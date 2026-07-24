## Function classes ------------------------------------------------------

#' dMod match function arguments
#' 
#' The function is exported for dependency reasons
#' 
#' @param arglist list
#' @param choices character
#' 
#' @export
match.fnargs <- function(arglist, choices) {

  # Catch the case of names == NULL
  if (is.null(names(arglist))) names(arglist) <- rep("", length(arglist))

  # exlude named arguments which are not in choices
  arglist <- arglist[names(arglist) %in% c(choices, "")]

  # determine available arguments
  available <- choices %in% names(arglist)

  if (!all(available)) names(arglist)[names(arglist) == ""] <- choices[!available]

  if (any(duplicated(names(arglist)))) stop("duplicate arguments in prdfn/obsfn/parfn function call")

  mapping <- match(choices, names(arglist))
  return(mapping)

}


## General concatenation of functions ------------------------------------------

#' Direct sum of objective functions
#'
#' @param x1 function of class `objfn`
#' @param x2 function of class `objfn`
#' @details The objective functions are evaluated and their results as added. Sometimes,
#' the evaluation of an objective function depends on results that have been computed
#' internally in a preceding objective function. Therefore, environments are forwarded
#' and all evaluations take place in the same environment. The first objective function
#' in a sum of functions generates a new environment.
#' @return Object of class `objfn`.
#' @seealso [normL2], [constraintL2], [priorL2], [datapointL2]
#' @aliases sumobjfn
#' @example inst/examples/objective.R
#' @export
"+.objfn" <- function(x1, x2) {

  if (is.null(x1)) return(x2)

  conditions.x1 <- attr(x1, "conditions")
  conditions.x2 <- attr(x2, "conditions")
  conditions12 <- union(conditions.x1, conditions.x2)

  parameters.x1 <- attr(x1, "parameters")
  parameters.x2 <- attr(x2, "parameters")
  parameters12 <- union(parameters.x1, parameters.x2)

  modelname.x1 <- attr(x1, "modelname")
  modelname.x2 <- attr(x2, "modelname")
  modelname12 <- union(modelname.x1, modelname.x2)


  # objfn + objfn
  if (inherits(x1, "objfn") & inherits(x2, "objfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = conditions12, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      # 1. If conditions.xi is null, always evaluate xi, but only once
      # 2. If not null, evaluate at intersection with conditions
      # 3. If not null & intersection is empty, don't evaluate xi at all
      v1 <- v2 <- NULL
      if (is.null(conditions.x1)) {
        v1 <- x1(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions.x1, env = env)
      } else if (any(conditions %in% conditions.x1)) {
        v1 <- x1(pars = pars, fixed = fixed, deriv = deriv, conditions = intersect(conditions, conditions.x1), env = env)
      }

      if (is.null(conditions.x2)) {
        v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions.x2, env = env)
      } else if (any(conditions %in% conditions.x2)) {
        v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, conditions = intersect(conditions, conditions.x2), env = attr(v1, "env"))
      }

      out <- v1 + v2
      attr(out, "env") <- attr(v1, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    # Propagate NLME reconstruction handles so a composed objective exposes its
    # model pieces (prdfn/data/errfn/omegaSpec) regardless of term order or
    # nesting. Coalesce from either operand. See .normalReconstruct() in nlmeNormal.R.
    for (.a in c("prdfn", "data", "errfn", "omegaSpec")) {
      .v <- attr(x1, .a, exact = TRUE)
      if (is.null(.v)) .v <- attr(x2, .a, exact = TRUE)
      if (!is.null(.v)) attr(outfn, .a) <- .v
    }
    # penaltySpec is MERGED (not first-wins): two constraintL1 terms combine
    # their penalty blocks under one shared lambda. See .mergePenaltySpec().
    .ps1 <- attr(x1, "penaltySpec", exact = TRUE)
    .ps2 <- attr(x2, "penaltySpec", exact = TRUE)
    .ps  <- .mergePenaltySpec(.ps1, .ps2)
    if (!is.null(.ps)) attr(outfn, "penaltySpec") <- .ps
    return(outfn)

  }


}


#' Multiplication of objective functions with scalars
#'
#' @description The `\%.*\%` operator allows to multiply objects of class objlist or objfn with
#' a scalar.
#'
#' @param x1 object of class objfn or objlist.
#' @param x2 numeric of length one.
#' @return An objective function or objlist object.
#'
#' @export
"%.*%" <- function(x1, x2) {

  if (inherits(x2, "objlist")) {

    out <- lapply(x2, function(x) {
      x1*x
    })
    # Multiply attributes
    out2.attributes <- attributes(x2)[sapply(attributes(x2), is.numeric)]
    attr.names <- names(out2.attributes)
    out.attributes <- lapply(attr.names, function(n) {
      x1*attr(x2, n)
    })
    attributes(out) <- attributes(x2)
    attributes(out)[attr.names] <- out.attributes

    return(out)


  } else if (inherits(x2, "objfn")) {

    conditions12 <- attr(x2, "conditions")
    parameters12 <- attr(x2, "parameters")
    modelname12 <- attr(x2, "modelname")
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = conditions12, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      v1 <- x1
      v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions, env = attr(v1, "env"))

      out <- v1 %.*% v2
      attr(out, "env") <- attr(v2, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    return(outfn)

  } else {

    x1*x2

  }

}


#' Direct sum of functions
#'
#' Used to add prediction function, parameter transformation functions or observation functions.
#'
#' @param x1 function of class `obsfn`, `prdfn` or `parfn`
#' @param x2 function of class `obsfn`, `prdfn` or `parfn`
#' @details Each prediction function is associated to a number of conditions. Adding functions
#' means merging or overwriting the set of conditions.
#' @return Object of the same class as `x1` and `x2` which returns results for the
#' union of conditions.
#' @aliases sumfn
#' @seealso [P], [Y], [Xs]
#' @example inst/examples/prediction.R
#' @export
"+.fn" <- function(x1, x2) {

  if (is.null(x1)) return(x2)

  mappings.x1 <- attr(x1, "mappings")
  mappings.x2 <- attr(x2, "mappings")

  conditions.x1 <- attr(x1, "conditions")
  conditions.x2 <- attr(x2, "conditions")
  overlap <- intersect(conditions.x1, conditions.x2)


  if (is.null(names(mappings.x1)) || is.null(names(mappings.x2))) stop("General transformations (NULL names) cannot be coerced.")

  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    mappings.x1 <- mappings.x1[!conditions.x1 %in% overlap]
    conditions.x1 <- conditions.x1[!conditions.x1 %in% overlap]
  }

  conditions.x12 <- c(conditions.x1, conditions.x2)
  mappings <- c(mappings.x1, mappings.x2)

  # prdfn + prdfn
  if (inherits(x1, "prdfn") & inherits(x2, "prdfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = names(mappings), env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]

      # Strip overlapping names so the kernels see disjoint (pars, fixed) sets.
      if (!is.null(fixed)) {
        pars  <- as.parvec(pars[setdiff(names(pars), names(fixed))])
        fixed <- as.parvec(fixed, deriv = FALSE, deriv2 = FALSE)
      }

      if (is.null(conditions)) {
        available <- names(mappings)
      } else {
        available <- intersect(names(mappings), conditions)
      }
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](times = times, pars = pars,
                                      fixed = fixed, deriv = deriv, deriv2 = deriv2)
      }

      out <- as.prdlist(outlist)
      return(out)

    }

    class(outfn) <- c("prdfn", "fn")

  }

  # obsfn + obsfn
  if (inherits(x1, "obsfn") & inherits(x2, "obsfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = names(mappings), env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]

      # Callers may pass overlapping `pars` / `fixed` (pinner contains
      # fixedinner); without this `c(pars, fixed)` inside the kernel would
      # duplicate names.
      if (!is.null(fixed)) {
        pars  <- as.parvec(pars[setdiff(names(pars), names(fixed))])
        fixed <- as.parvec(fixed, deriv = FALSE, deriv2 = FALSE)
      }

      if (is.null(conditions)) {
        available <- names(mappings)
      } else {
        available <- intersect(names(mappings), conditions)
      }
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](out = out, pars = pars,
                                      fixed = fixed, deriv = deriv, deriv2 = deriv2)
      }

      out <- as.prdlist(outlist)
      return(out)

    }

    class(outfn) <- c("obsfn", "fn")

  }


  # parfn + parfn
  if (inherits(x1, "parfn") & inherits(x2, "parfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = names(mappings), env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]


      if (is.null(conditions)) {
        available <- names(mappings)
      } else {
        available <- intersect(names(mappings), conditions)
      }
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2)
      }

      return(outlist)

    }

    class(outfn) <- c("parfn", "fn")

  }


  # detach outfn from the S3 dispatch frame; it only needs `mappings`
  environment(outfn) <- list2env(list(mappings = mappings),
                                 parent = asNamespace("dMod"))

  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- union(attr(x1, "parameters"), attr(x2, "parameters"))
  attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(x1, "compileInfo"), attr(x2, "compileInfo"))
  attr(outfn, "conditions") <- conditions.x12
  attr(outfn, "forcings") <- do.call(c, list(attr(x1, "forcings"), attr(x2, "forcings")))

  return(outfn)

}

#' Direct sum of datasets
#'
#' Used to merge datasets with overlapping conditions.
#'
#' @param data1 dataset of class `datalist`
#' @param data2 dataset of class `datalist`
#' @details Each data list contains data frames for a number of conditions.
#' The direct sum of datalist is meant as merging the two data lists and
#' returning the overarching datalist.
#' @return Object of class `datalist` for the
#' union of conditions.
#' @aliases sumdatalist
#' @example inst/examples/sumdatalist.R
#' @export
"+.datalist" <- function(data1, data2) {

  overlap <- names(data2)[names(data2) %in% names(data1)]
  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    data1 <- data1[!names(data1) %in% names(data2)]
  }

  conditions <- union(names(data1), names(data2))
  data <- lapply(conditions, function(C) rbind(data1[[C]], data2[[C]]))
  names(data) <- conditions

  grid1 <- attr(data1, "condition.grid")
  grid2 <- attr(data2, "condition.grid")

  grid <- combine(grid1, grid2)




  if (is.data.frame(grid)) grid <- grid[!duplicated(rownames(grid)), , drop = FALSE]

  out <- as.datalist(data)
  attr(out, "condition.grid") <- grid

  return(out)
}

out_conditions <- function(c1, c2) {

  if (!is.null(c1)) return(c1)
  if (!is.null(c2)) return(c2)
  return(NULL)

}

test_conditions <- function(c1, c2) {
  if (is.null(c1)) return(NULL)
  if (is.null(c2)) return(NULL)
  return(intersect(c1, c2))
}

#' Concatenation of functions
#'
#' Used to concatenate observation functions, prediction functions and parameter transformation functions.
#'
#' @param p1 function of class `obsfn`, `prdfn`, `parfn` or `idfn`
#' @param p2 function of class `obsfn`, `prdfn`, `parfn` or `idfn`
#' @return Object of the same class as `x1` and `x2`.
#' @aliases prodfn
#' @example inst/examples/prediction.R
#' @export
"*.fn" <- function(p1, p2) {
  
  # ============================================================
  # Global consistency check for condition handling
  #
  # Rules:
  # - A condition-unspecific function (conditions = NULL) may be
  #   combined with any other function.
  # - Two condition-specific functions must cover the same set
  #   of conditions.
  # - It is NOT allowed to combine a single-condition function
  #   with a multi-condition function.
  # ============================================================
  
  conditions.p1 <- attr(p1, "conditions")
  conditions.p2 <- attr(p2, "conditions")
  
  is_unspecific <- function(x) is.null(x)
  is_specific   <- function(x) !is.null(x) && length(x) == 1
  is_multiple   <- function(x) !is.null(x) && length(x) > 1
  
  if (!is_unspecific(conditions.p1) &&
      !is_unspecific(conditions.p2)) {
    
    # one specific, one multiple -> forbidden
    if ((is_specific(conditions.p1) && is_multiple(conditions.p2)) ||
        (is_specific(conditions.p2) && is_multiple(conditions.p1))) {
      
      stop(
        "Invalid composition of functions:\n",
        "Incompatible condition sets.\n\n",
        "Left-hand function conditions:  ",
        paste(conditions.p1, collapse = ", "), "\n",
        "Right-hand function conditions: ",
        paste(conditions.p2, collapse = ", "), "\n\n",
        "A function defined for a single condition cannot be\n",
        "combined with a function defined for multiple conditions.\n",
        "Either both functions must cover all conditions,\n",
        "or one function must be condition-unspecific."
      )
    }
  }

  # obsfn * obsfn -> obsfn
  if (inherits(p1, "obsfn") & inherits(p2, "obsfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]


      step1 <- p2(out = out, pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = attr(step1[[i]], "parameters"), fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = names(step1)[i])))


      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for observation function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
        outfn(out = out, pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])



      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings

    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"), attr(p2, "compileInfo"))
    class(outfn) <- c("obsfn", "fn", "composed")

    return(outfn)

  }


  # obsfn * parfn -> obsfn
  if (inherits(p1, "obsfn") & inherits(p2, "parfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = out, pars = step1[[i]], fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = names(step1)[i])))

      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for observation function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
        outfn(out = out, pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings

    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"), attr(p2, "compileInfo"))
    class(outfn) <- c("obsfn", "fn", "composed")

    return(outfn)

  }


  # obsfn * prdfn -> prdfn
  if (inherits(p1, "obsfn") & inherits(p2, "prdfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]

      step1 <- p2(times = times, pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) {
        pinner <- attr(step1[[i]], "parameters")
        fixedinner <- pinner[attr(pinner, "fixed")]
        p1(out = step1[[i]],
           pars = pinner,
           fixed = fixedinner,
           deriv = deriv,
           deriv2 = deriv2,
           conditions = names(step1)[i])
      }))


      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for prediction function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
        outfn(times = times, pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings

    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"), attr(p2, "compileInfo"))
    class(outfn) <- c("prdfn", "fn", "composed")

    return(outfn)

  }


  # prdfn * parfn -> prdfn
  if (inherits(p1, "prdfn") & inherits(p2, "parfn")) {


    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) {
        p1(times = times,
           pars = (step1[[i]])[setdiff(names(step1[[i]]), attr(step1[[i]], "fixed"))],
           fixed = (step1[[i]])[attr(step1[[i]], "fixed")],
           deriv = deriv,
           deriv2 = deriv2,
           conditions = names(step1)[i])
      }))

      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for prediction function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
        outfn(times = times, pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.out[i])[[1]]
      }
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings


    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"), attr(p2, "compileInfo"))
    class(outfn) <- c("prdfn", "fn", "composed")

    return(outfn)

  }

  # parfn * parfn -> parfn
  if (inherits(p1, "parfn") & inherits(p2, "parfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)


    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i)
        p1(pars = (step1[[i]])[setdiff(names(step1[[i]]), attr(step1[[i]], "fixed"))],
           fixed = (step1[[i]])[attr(step1[[i]], "fixed")],
           deriv = deriv, deriv2 = deriv2, conditions = names(step1)[i])))
      return(step2)

    }

    # Generate mappings for parameters function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
        outfn(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings


    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"), attr(p2, "compileInfo"))
    class(outfn) <- c("parfn", "fn", "composed")

    return(outfn)

  }


  # objfn * parfn -> objfn
  if (inherits(p1, "objfn") & inherits(p2, "parfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(...,  fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, "pars")]
      pars <- arglist[[1]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions)
      step2 <- Reduce("+", lapply(1:length(step1), function(i) p1(pars = step1[[i]], fixed = NULL, deriv = deriv, deriv2 = deriv2, conditions = names(step1)[i], env = env)))
      return(step2)


    }

    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"), attr(p2, "compileInfo"))
    class(outfn) <- c("objfn", "fn", "composed")

    return(outfn)


  }

  # idfn * fn -> fn
  if (inherits(p1, "idfn")) {
    return(p2)
  }

  # fn * idfn -> fn
  if (inherits(p2, "idfn")) {
    return(p1)
  }

}

