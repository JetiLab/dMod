
# Pharmacometric abbreviations used in the NLME plot helpers below:
#   DV    observed value
#   IPRED individual prediction (per subject, with eta_i*)
#   PRED  population prediction (eta = 0)
#   IRES  individual residual = DV - IPRED
#   IWRES individual weighted residual = (DV - IPRED) / sigma

# ggplot2 / dplyr NSE column references; declared so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c("IPRED", "PRED", "predicted", "observed",
                         "sd_est", "iter", "level",
                         # plot.sparsify ggplot2 aes() variables
                         "subject", "cluster", "centroid", "G", "score",
                         "selected"))


# Custom interface to ggplot2 ---

#' Open last plot in external pdf viewer
#' 
#' @description Convenience function to show last plot in an external viewer.
#' @param plot `ggplot2` plot object.
#' @param command character, indicatig which pdf viewer is started.
#' @param ... arguments going to `ggsave`.
#' @export
ggopen <- function(plot = last_plot(), command = "xdg-open", ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  system(command = paste(command, filename))
}




#' Standard plotting theme of dMod
#' 
#' @param base_size numeric, font-size
#' @param base_family character, font-name
#' @export
theme_dMod <- function(base_size = 12, base_family = "") {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]

  theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = "black"),
          axis.text = element_text(size = rel(1.0), colour = "black"),
          axis.text.x = element_text(margin = margin(t = 4, r = 4, b = 0, l = 4, unit = "mm")),
          axis.text.y = element_text(margin = margin(t = 4, r = 4, b = 4, l = 0, unit = "mm")),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = "black"),
          # panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(1.0)))
  
}

dMod_colors <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E",
                 "#4B1F6F", "#CC79A7", "#006400", "#F0E442", "#8B4513")

#' Generate `n` distinct colors from the dMod palette
#'
#' Returns the first `n` seed colors directly. For `n` larger than the seed
#' set, `Polychrome::createPalette()` extends the palette deterministically;
#' if Polychrome is not installed, the overflow falls back to gray.
#'
#' @param n integer, number of colors to produce.
#' @return Character vector of length `n` with hex color codes.
#' @export
dMod_palette <- function(n) {
  n <- as.integer(n)
  if (n <= 0L) return(character(0))
  if (n <= length(dMod_colors)) {
    return(unname(dMod_colors[seq_len(n)]))
  }
  if (requireNamespace("Polychrome", quietly = TRUE)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(123L)
    # createPalette() shifts the seeds slightly to maximize distinctness across
    # the whole set; keep the original seeds verbatim and only borrow the tail.
    extended <- unname(Polychrome::createPalette(n, seedcolors = dMod_colors))
    return(c(dMod_colors, extended[(length(dMod_colors) + 1L):n]))
  }
  c(dMod_colors, rep("gray", n - length(dMod_colors)))
}

#' Standard dMod color palette
#'
#' @param ... arguments forwarded to [ggplot2::discrete_scale()].
#' @export
#' @examples
#' library(ggplot2)
#' times <- seq(0, 2*pi, 0.1)
#' values <- sin(times)
#' data <- data.frame(
#'    time = times,
#'    value = c(values, 1.2*values, 1.4*values, 1.6*values),
#'    group = rep(c("C1", "C2", "C3", "C4"), each = length(times))
#' )
#' qplot(time, value, data = data, color = group, geom = "line") +
#'    theme_dMod() + scale_color_dMod()
#' @export
scale_color_dMod <- function(...) {
  ggplot2::discrete_scale(aesthetics = "colour", palette = dMod_palette, ...)
}


#' Standard dMod fill scheme
#'
#' @export
#' @param ... arguments forwarded to [ggplot2::discrete_scale()].
scale_fill_dMod <- function(...) {
  ggplot2::discrete_scale(aesthetics = "fill", palette = dMod_palette, ...)
}

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_dMod() + theme_dMod()




# Other ---------------------------------------------

#' Coordinate transformation for data frames
#' 
#' Applies a symbolically defined transformation to the `value`
#' column of a data frame. Additionally, if a `sigma` column is
#' present, those values are transformed according to Gaussian error
#' propagation.
#' @param data data frame with at least columns "name" (character) and
#' "value" (numeric). Can optionally contain a column "sigma" (numeric).
#' @param transformations character (the transformation) or named list of
#' characters. In this case, the list names must be a subset of those 
#' contained in the "name" column.
#' @return The data frame with the transformed values and sigma uncertainties.
#' @export
#' 
#' @examples
#' mydata1 <- data.frame(name = c("A", "B"), time = 0:5, value = 0:5, sigma = .1)
#' coordTransform(mydata1, "log(value)")
#' coordTransform(mydata1, list(A = "exp(value)", B = "sqrt(value)"))
coordTransform <- function(data, transformations) {
  
  mynames <- unique(as.character(data$name))
  
  # Replicate transformation if not a list
  if (!is.list(transformations))
    transformations <- as.list(structure(rep(transformations, length(mynames)), names = mynames))
  
  out <- do.call(rbind, lapply(mynames, function(n) {
    
    subdata <- subset(data, name == n)
    
    if (n %in% names(transformations)) {
      
      mysymbol <- getSymbols(transformations[[n]])[1]
      mytrafo <- replaceSymbols(mysymbol, "value", transformations[[n]])
      mytrafo <- parse(text = mytrafo)
      
      if ("sigma" %in% colnames(subdata))
        subdata$sigma <- abs(with(subdata, eval(D(mytrafo, "value")))) * subdata$sigma
      subdata$value <- with(subdata, eval(mytrafo))
      
    }
    
    return(subdata)
    
  }))
  
  
  return(out)
  
  
}


# Method dispatch for plotX functions -------------



#' Plot a list of model predictions
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param ... Further arguments going to `dplyr::filter`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet Either `"wrap"` or `"grid"`
#' @param transform list of transformation for the states, see [coordTransform].
#' @details The data.frame being plotted has columns `time`, `value`, `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @import ggplot2
#' @example inst/examples/plotting.R
#' @export
plotPrediction <- function(prediction,...) {
  UseMethod("plotPrediction", prediction)
}




#' Plot a list of model predictions and a list of data points in a combined plot
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param data Named list of data.frames as being used in [res], i.e. with columns `name`, `time`, 
#' `value` and `sigma`.
#' @param ... Further arguments going to `dplyr::filter`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet `"wrap"` or `"grid"`. Try `"wrap_plain"` for high amounts of conditions and low amounts of observables.
#' @param transform list of transformation for the states, see [coordTransform].
#' @param aesthetics Named list of aesthetic mappings, specified as character, e.g. `list(linetype = "name")`. 
#' Can refer to variables in the condition.grid
#' @details The data.frame being plotted has columns `time`, `value`, `sigma`,
#' `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @example inst/examples/plotting.R
#' @importFrom graphics par
#' @export
plotCombined <- function(prediction,...) {
  UseMethod("plotCombined", prediction)
}


#' Plot a list data points
#' 
#' @param data Named list of data.frames as being used in [res], i.e. with columns `name`, `time`, 
#' `value` and `sigma`.
#' @param ... Further arguments going to `subset`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet Either `"wrap"` or `"grid"`
#' @param transform list of transformation for the states, see [coordTransform].
#' @details The data.frame being plotted has columns `time`, `value`, `sigma`,
#' `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @example inst/examples/plotting.R
#' @export
plotData  <- function(data,...) {
  UseMethod("plotData", data)
}

#' @export
#' @rdname plotData
plotData.data.frame <- function(data, ...) {
  plotData.datalist(as.datalist(data), ...)
}

#' Profile likelihood plot
#' 
#' @param profs Lists of profiles as being returned by [profile].
#' @param ... logical going to subset before plotting.
#' @param maxvalue Numeric, the value where profiles are cut off.
#' @param parlist Matrix or data.frame with columns for the parameters to be added to the plot as points.
#' If a "value" column is contained, deltas are calculated with respect to lowest chisquare of profiles.
#' @param ncol Number of columns in the resulting plot grid.
#' @return A plot object of class `ggplot`.
#' @details See [profile] for examples.
#' @export
plotProfile <- function(profs,...) {
  UseMethod("plotProfile", profs)
}


#' Profile likelihood: plot of the parameter paths.
#' 
#' @param profs profile or list of profiles as being returned by [profile]
#' @param ... arguments going to subset
#' @param whichPar Character or index vector, indicating the parameters that are taken as possible reference (x-axis)
#' @param sort Logical. If paths from different parameter profiles are plotted together, possible
#' combinations are either sorted or all combinations are taken as they are.
#' @param relative logical indicating whether the origin should be shifted.
#' @param scales character, either `"free"` or `"fixed"`.
#' @return A plot object of class `ggplot`.
#' @details See [profile] for examples.
#' @export
plotPaths <- function(profs, ..., whichPar = NULL, sort = FALSE, relative = TRUE, scales = "fixed") {
  
  if ("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- as.data.frame(arglist[[i]])
    parameters <- attr(arglist[[i]], "parameters")
    
    if (is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    if (is.null(whichPar)) whichPar <- names(proflist)
    if (is.numeric(whichPar)) whichPar <- names(proflist)[whichPar]
    
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # matirx
      paths <- as.matrix(proflist[[n]][, parameters])
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      if (relative) 
        for(j in 1:ncol(paths)) paths[, j] <- as.numeric(paths[, j]) - as.numeric(paths[origin, j])
      
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if (sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- submatrix(combinations, cols = -which(combinations[1,] == combinations[2,]))
      combinations <- submatrix(combinations, cols = !duplicated(paste(combinations[1,], combinations[2,])))
      
      
      
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = profnames[i],
                   combination = paste(combinations[,j], collapse = " - \n"),
                   x = paths[, combinations[1,j]],
                   y = paths[, combinations[2,j]])
      }))
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  
  if (relative)
    axis.labels <- c(expression(paste(Delta, "parameter 1")), expression(paste(Delta, "parameter 2")))  
  else
    axis.labels <- c("parameter 1", "parameter 2")
  
  
  data <- droplevels(subset(data, ...))
  data$y <- as.numeric(data$y)
  data$x <- as.numeric(data$x)
  
  suppressMessages(
    p <- ggplot(data, aes(x = x, y = y, group = interaction(name, proflist), color = name, lty = proflist)) + 
      facet_wrap(~combination, scales = scales) + 
      geom_path() + #geom_point(aes=aes(size=1), alpha=1/3) +
      xlab(axis.labels[1]) + ylab(axis.labels[2]) +
      scale_linetype_discrete(name = "profile\nlist") +
      scale_color_dMod(name = "profiled\nparameter")
  )
  
  attr(p, "data") <- data
  return(p)
  
}



#' Plot Fluxes given a list of flux Equations
#'
#' @param pouter parameters
#' @param x The model prediction function `x(times, pouter, fixed, ...)`
#' @param fluxEquations list of chars containing expressions for the fluxes,
#' if names are given, they are shown in the legend. Easy to obtain via [subset.eqnlist], see Examples.
#' @param nameFlux character, name of the legend.
#' @param times Numeric vector of time points for the model prediction
#' @param ... Further arguments going to x, such as `fixed` or `conditions`
#'
#'
#' @return A plot object of class `ggplot`.
#' @examples
#' \dontrun{
#'
#' plotFluxes(bestfit, x, times, subset(f, "B"%in%Product)$rates, nameFlux = "B production")
#' }
#' @export
plotFluxes <- function(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", ...){

  if (is.null(names(fluxEquations))) names(fluxEquations) <- fluxEquations

  flux <- funCpp(fluxEquations, convenient = FALSE)$func
  prediction.all <- x(times, pouter, deriv = FALSE, ...)
  names.prediction.all <- names(prediction.all)
  if (is.null(names.prediction.all)) names.prediction.all <- paste0("C", 1:length(prediction.all))

  out <- lapply(1:length(prediction.all), function(cond) {
    prediction <- prediction.all[[cond]]
    pinner <- attr(prediction, "parameters")
    pinner.matrix <- matrix(pinner, nrow = length(pinner), ncol = nrow(prediction),
                            dimnames = list(names(pinner), NULL))
    fluxes <- cbind(time = prediction[, "time"], flux(cbind(prediction, t(pinner.matrix))))
    return(fluxes)
  }); names(out) <- names.prediction.all
  out <- wide2long(out)

  cbPalette <- c("#999999", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2",
                 "#D55E00", "#CC79A7","#CC6666", "#9999CC", "#66CC99","red", "blue", "green","black")

  P <- ggplot(out, aes(x = time, y = value, group = name, fill = name, log = "y")) +
    facet_wrap(~condition) + scale_fill_manual(values = cbPalette, name = nameFlux) +
    geom_density(stat = "identity", position = "stack", alpha = 0.3, color = "darkgrey", linewidth = 0.4) +
    xlab("time") + ylab("flux contribution")

  attr(P, "out") <- out

  return(P)

}


.stepDetect <- function(x, tol) {
  
  jumps <- 1
  while (TRUE) {
    i <- which(x - x[1] > tol)[1]
    if (is.na(i)) break
    jumps <- c(jumps, tail(jumps, 1) - 1 + i)
    x <- x[-seq(1, i - 1, 1)]
  }
  
  return(jumps)
  
  
}

#' Plotting objective values of a collection of fits
#' 
#' @param x data.frame with columns "value", "converged" and "iterations", e.g. 
#' a [parframe].
#' @param ... arguments for subsetting of x
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
#' @param showSteps logical, if `TRUE`, the detected steps are indicated by
#' dashed vertical lines and labelled by their index. Defaults to `FALSE`.
#' @export
plotValues <- function(x,...) {
  UseMethod("plotValues", x)
}



#' Plot parameter values for a fitlist
#' 
#' @param x parameter frame as obtained by as.parframe(mstrust)
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
#' @param ... arguments for subsetting of x
#' @export
plotPars <- function(x,...) {
  UseMethod("plotPars", x)
}



#' Plot residuals for a fitlist
#'
#' @description
#' Creates a plot of residuals from model fits, with flexible options for
#' grouping and faceting. Residuals can be summarized across different
#' dimensions (time, condition, observable, fit index).
#'
#' @param parframe Object of class \code{parframe}, e.g. returned by \link{mstrust}.
#' @param x Prediction function returning named list of data.frames with names 
#'   matching \code{data}.
#' @param data A \code{datalist} object, i.e. named list of data.frames with 
#'   columns \code{name}, \code{time}, \code{value}, and \code{sigma}.
#' @param split Character vector specifying how to summarize and display residuals.
#'   \itemize{
#'     \item \code{split[1]}: Variable for x-axis
#'     \item \code{split[2]}: Variable for grouping (color/line), defaults to \code{split[1]}
#'     \item \code{split[3+]}: Additional variables for \code{facet_wrap()}
#'   }
#' @param errmodel Optional error model function of type \code{prdfn}. If provided,
#'   residuals include the log-likelihood contribution from sigma.
#' @param ... Additional arguments passed to the prediction function \code{x}.
#'
#' @return A \code{ggplot} object with the summarized residual data frame
#'   attached as attribute \code{"out"}.
#'
#' @examples
#' \dontrun{
#' # Time on x-axis, faceted by condition and name
#' plotResiduals(myfitlist, g * x * p, data, 
#'               c("time", "index", "condition", "name"), 
#'               conditions = myconditions[1:4])
#'
#' # Condition on x-axis, residuals summed over time
#' plotResiduals(myfitlist, g * x * p, data, c("condition", "name", "index"))
#' }
#'
#' @export
#' @importFrom dplyr group_by summarise across
#' @importFrom rlang data_sym syms
plotResiduals <- function(parframe, x, data, split = "condition", errmodel = NULL, ...) {

  # Internal dispatch: NLME .fitNormal objects get their own residual diagnostic
  # (IWRES vs IPRED + vs TIME, see plotResiduals..fitNormal below). The classical
  # parframe/x/data path below is preserved unchanged for back-compat.
  if (inherits(parframe, "em")) return(plotResiduals.em(parframe))

  timesD <- sort(unique(c(0, unlist(lapply(data, function(d) d$time)))))
  
  if (!("index" %in% colnames(parframe))) {
    parframe$index <- seq_len(nrow(parframe))
  }
  
  # --- Compute residuals for all fits and conditions ---
  out <- do.call(rbind, lapply(seq_len(nrow(parframe)), function(j) {
    pred <- x(timesD, as.parvec(parframe, j), deriv = FALSE, ...)
    
    out_con <- do.call(rbind, lapply(names(pred), function(con) {
      err <- NULL
      if (!is.null(errmodel)) {
        err <- errmodel(out = pred[[con]], pars = getParameters(pred[[con]]), conditions = con)
      }
      out <- res(data[[con]], pred[[con]], err[[con]])
      cbind(out, condition = con)
    }))
    
    cbind(index = as.character(parframe[j, "index"]), out_con)
  }))
  
  # --- Summarize residuals ---
  out <- dplyr::group_by(out, across(all_of(split)))
  
  if (!is.null(errmodel)) {
    out <- dplyr::summarise(out, res = sum(weighted.residual^2 + log(sigma^2)), .groups = "drop")
  } else {
    out <- dplyr::summarise(out, res = sum(weighted.residual^2), .groups = "drop")
  }
  
  out <- as.data.frame(out)
  
  # --- Build aesthetics ---
  groupvar <- if (length(split) > 1) split[2] else split[1]
  
  p <- ggplot(out, aes(x = !!rlang::data_sym(split[1]), 
                       y = res, 
                       color = !!rlang::data_sym(groupvar), 
                       group = !!rlang::data_sym(groupvar))) + 
    theme_dMod() + 
    geom_point() + 
    geom_line()
  
  if (length(split) > 2) {
    facet_vars <- rlang::syms(split[3:length(split)])
    p <- p + facet_wrap(vars(!!!facet_vars))
  }
  
  attr(p, "out") <- out
  p
}


# NLME plot helpers ---------------------------------------------------------

#' Predictions from an EM object
#'
#' @description
#' Returns a long-format `data.frame` of observed values vs population (`PRED`,
#' eta = 0) and individual (`IPRED`, eta at posterior modes) predictions per
#' condition, plus residuals. Used directly by the diagnostic plot helpers.
#'
#' @param object An [EM] (from [EM]).
#' @param times Optional numeric vector of additional times for the smooth
#'   IPRED/PRED curves. Defaults to the union of observed times.
#' @param ... Ignored.
#'
#' @return A long-format `data.frame` with columns
#'   `condition, time, name, observed, sigma, IPRED, PRED, IRES, IWRES, PRES,
#'   PWRES, source` where `source = "obs"` for rows aligned with data
#'   observations and `source = "grid"` for the dense IPRED/PRED smooth grid
#'   (observed/sigma NA for grid rows).
#' @export
predict.em <- function(object, times = NULL, ...) {
  .emRequireOmega(object, "predict")
  if (is.null(object$prdfn) || is.null(object$data) || is.null(object$omega))
    stop("predict..fitNormal: fit is missing `prdfn`, `data`, or `omega` ",
         "(was the fit built by EM()?).")

  prdfn <- object$prdfn
  data  <- object$data
  om    <- object$omega
  etaModes <- object$etaModes
  pars      <- object$argument

  subjects <- rownames(om$subjectEtas)
  N <- length(subjects)
  obs_times <- sort(unique(unlist(lapply(data, `[[`, "time"))))
  grid_times <- if (is.null(times)) sort(unique(c(0, obs_times)))
                else sort(unique(c(0, obs_times, times)))

  eta_full <- as.vector(etaModes)
  names(eta_full) <- as.vector(om$subjectEtas)
  pars_ipred <- c(pars, eta_full)
  pars_pred  <- c(pars, setNames(rep(0, length(eta_full)), names(eta_full)))

  pred_ipred <- prdfn(times = grid_times, pars = pars_ipred, conditions = subjects)
  pred_pred  <- prdfn(times = grid_times, pars = pars_pred,  conditions = subjects)

  # Restrict to observables that actually appear in the data. Without this
  # filter, internal model states emitted by the prdfn (e.g. the depot Ag in a
  # PK model, present because the observation Y() is built with
  # attach.input = TRUE) would be plotted as ghost curves on the IPRED/PRED
  # axis next to the real observable, swamping the y-axis with the dose.
  data_names <- unique(unlist(lapply(data, `[[`, "name")))
  obs_names <- intersect(setdiff(colnames(pred_ipred[[1]]), "time"),
                         data_names)

  # If an errfn is attached and the data table has NA sigmas, evaluate
  # the errfn per-condition on the prediction so IWRES diagnostics get a
  # meaningful weight. Mirrors evalConditionResidual's contract.
  err <- object$errfn
  sigma_ipred <- NULL
  if (!is.null(err)) {
    sigma_ipred <- setNames(lapply(subjects, function(s) {
      pinner <- getParameters(pred_ipred[[s]])
      err(out = pred_ipred[[s]], pars = pinner, conditions = s)[[s]]
    }), subjects)
  }

  rows <- list()
  for (s in subjects) {
    d_s <- data[[s]]
    pi  <- pred_ipred[[s]]
    pp  <- pred_pred[[s]]
    si  <- if (!is.null(sigma_ipred)) sigma_ipred[[s]] else NULL
    # Coerce prdframe class away so plain matrix indexing works.
    if (!is.null(si)) si <- unclass(si)
    for (nm in obs_names) {
      d_nm <- d_s[d_s$name == nm, , drop = FALSE]
      get_sigma <- function(row_idx) {
        if (!is.null(si) && nm %in% colnames(si) && all(row_idx <= nrow(si)))
          as.numeric(si[row_idx, nm])
        else rep(NA_real_, length(row_idx))
      }
      if (nrow(d_nm)) {
        idx <- match(d_nm$time, pi[, "time"])
        sig_obs <- if (!is.null(si)) get_sigma(idx)
                   else if (all(is.na(d_nm$sigma))) rep(NA_real_, nrow(d_nm))
                   else d_nm$sigma
        rows[[length(rows) + 1L]] <- data.frame(
          condition = s, time = d_nm$time, name = nm,
          observed  = d_nm$value, sigma = sig_obs,
          IPRED     = pi[idx, nm], PRED  = pp[idx, nm],
          source    = "obs", stringsAsFactors = FALSE
        )
      }
      grid_idx <- setdiff(seq_len(nrow(pi)), match(d_nm$time, pi[, "time"]))
      if (length(grid_idx)) {
        rows[[length(rows) + 1L]] <- data.frame(
          condition = s, time = pi[grid_idx, "time"], name = nm,
          observed  = NA_real_, sigma = get_sigma(grid_idx),
          IPRED     = pi[grid_idx, nm], PRED = pp[grid_idx, nm],
          source    = "grid", stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$IRES  <- out$observed - out$IPRED
  out$PRES  <- out$observed - out$PRED
  out$IWRES <- out$IRES / out$sigma
  out$PWRES <- out$PRES / out$sigma
  rownames(out) <- NULL
  out
}



#' Per-subject individual fits (spaghetti plot)
#'
#' @description Faceted plot with one panel per subject: observed dots, IPRED
#'   curve, and (optionally) the population PRED curve overlaid dashed.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot.
#' @export
plotIndivs <- function(x, ...) UseMethod("plotIndivs", x)


#' Per-subject individual fits for an EM
#'
#' @description Per-subject IPRED curve with an IPRED +/- sigma ribbon derived
#'   from the fit's attached errfn, optional population PRED overlay dashed,
#'   and observed values as points. When the fit has a single observable the
#'   layout is `facet_wrap(~ condition)`; with multiple observables it switches
#'   to `facet_grid(name ~ condition)` (observables in rows, subjects in
#'   columns). Use `subjectsPerPage` to split very large cohorts into several
#'   plots.
#' @param x An [EM].
#' @param times Optional grid of additional times for the smooth IPRED/PRED.
#' @param ncol Facet column count for the single-observable
#'   `facet_wrap(~ condition)` layout (default 4). Ignored in
#'   `facet_grid(name ~ condition)` mode.
#' @param showPred Logical; overlay population PRED dashed (default TRUE).
#' @param showBand Logical; draw the IPRED +/- sigma ribbon if the fit carries
#'   an errfn (default TRUE).
#' @param subjectsPerPage Optional integer. If set, subjects are split into
#'   pages of at most `subjectsPerPage` each and the function returns a list
#'   of ggplots (one per page, page index appended to the title). `NULL`
#'   (default) keeps the single-plot behaviour.
#' @param ... Ignored.
#' @return A ggplot, or a list of ggplots when `subjectsPerPage` is set.
#' @export
plotIndivs.em <- function(x, times = NULL, ncol = 4L,
                              showPred = TRUE, showBand = TRUE,
                              subjectsPerPage = NULL, ...) {
  .emRequireOmega(x, "plotIndivs")
  fit <- x
  obs_times <- sort(unique(unlist(lapply(fit$data, `[[`, "time"))))
  if (is.null(times))
    times <- seq(min(obs_times), max(obs_times), length.out = 200L)

  pf <- predict(fit, times = times)
  # Trim grid rows to the observed time range. predict..fitNormal prepends t=0 to
  # the dense grid for the ODE solver, but if no subject has an observation at
  # that time the grid IPRED/PRED there can be far outside the data range
  # (e.g. log(Cc + eps) at Cc(0)=0 sits at log(eps)), which would dominate the
  # y-axis. Observation rows are kept verbatim.
  pf <- pf[pf$source == "obs" |
             (pf$time >= min(obs_times) & pf$time <= max(obs_times)),
           , drop = FALSE]

  n_obs_names <- length(unique(pf$name))
  cond_levels <- if (is.factor(pf$condition))
    levels(droplevels(pf$condition)) else unique(as.character(pf$condition))

  build_page <- function(page_pf, page_label = NULL) {
    obs <- page_pf[page_pf$source == "obs", , drop = FALSE]
    grd <- page_pf[order(page_pf$condition, page_pf$name, page_pf$time), ,
                   drop = FALSE]
    has_band <- showBand && any(is.finite(grd$sigma))

    p <- ggplot2::ggplot()
    if (has_band) {
      p <- p + ggplot2::geom_ribbon(
        data = grd,
        ggplot2::aes(x = time,
                     ymin = IPRED - sigma, ymax = IPRED + sigma),
        fill = dMod_colors[3], alpha = 0.2, linetype = 0)
    }
    p <- p +
      ggplot2::geom_line(data = grd,
                         ggplot2::aes(x = time, y = IPRED, color = "IPRED"),
                         linewidth = 0.7)
    if (showPred) {
      p <- p + ggplot2::geom_line(
        data = grd,
        ggplot2::aes(x = time, y = PRED, color = "PRED"),
        linewidth = 0.5, linetype = "dashed")
    }
    p <- p +
      ggplot2::geom_point(data = obs,
                          ggplot2::aes(x = time, y = observed),
                          size = 1.6, alpha = 0.85)

    if (n_obs_names > 1L) {
      p <- p + ggplot2::facet_grid(name ~ condition, scales = "free_y")
    } else {
      p <- p + ggplot2::facet_wrap(~ condition, ncol = ncol, scales = "free_y")
    }

    title <- sprintf("Individual fits (method = %s)", fit$method)
    if (!is.null(page_label))
      title <- paste0(title, " ", page_label)

    p +
      ggplot2::scale_color_manual(values = c(IPRED = dMod_colors[3],
                                             PRED  = dMod_colors[2])) +
      ggplot2::labs(x = "Time", y = "Value", color = NULL, title = title) +
      theme_dMod(base_size = 11)
  }

  if (is.null(subjectsPerPage))
    return(build_page(pf))

  subjectsPerPage <- as.integer(subjectsPerPage)
  if (length(subjectsPerPage) != 1L || is.na(subjectsPerPage) ||
      subjectsPerPage < 1L)
    stop("`subjectsPerPage` must be a positive integer or NULL.")

  pages <- split(cond_levels,
                 ceiling(seq_along(cond_levels) / subjectsPerPage))
  n_pages <- length(pages)
  lapply(seq_along(pages), function(i) {
    page_conds <- pages[[i]]
    sub <- pf[as.character(pf$condition) %in% page_conds, , drop = FALSE]
    if (is.factor(sub$condition))
      sub$condition <- factor(sub$condition, levels = page_conds)
    label <- sprintf("(page %d/%d)", i, n_pages)
    build_page(sub, page_label = label)
  })
}



#' Observed vs predicted scatter (DV vs IPRED and DV vs PRED)
#'
#' @description S3 plot method for [EM]. Two-panel scatter: DV vs IPRED
#'   on the left, DV vs PRED on the right, identity line shown dashed. With a
#'   multi-observable fit the observable becomes a faceting row
#'   (`name ~ panel`); otherwise one row with `facet_wrap(~ panel)`.
#' @param x An [EM].
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.em <- function(x, ...) {
  .emRequireOmega(x, "plot")
  fit <- x
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- rbind(
    data.frame(condition = pf$condition, name = pf$name,
               observed = pf$observed, predicted = pf$IPRED,
               panel = "IPRED", stringsAsFactors = FALSE),
    data.frame(condition = pf$condition, name = pf$name,
               observed = pf$observed, predicted = pf$PRED,
               panel = "PRED", stringsAsFactors = FALSE))
  multi_obs <- length(unique(long$name)) > 1L
  p <- ggplot2::ggplot(long, ggplot2::aes(x = predicted, y = observed,
                                          color = condition)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(alpha = 0.75)
  if (multi_obs) {
    # Per-observable scales differ by orders of magnitude (e.g. log-cp vs
    # linear PCA); free scales are required and coord_equal is incompatible
    # with that.
    p <- p + ggplot2::facet_grid(name ~ panel, scales = "free")
  } else {
    rng <- range(c(long$observed, long$predicted), na.rm = TRUE)
    p <- p +
      ggplot2::coord_equal(xlim = rng, ylim = rng) +
      ggplot2::facet_wrap(~ panel)
  }
  p +
    scale_color_dMod() +
    ggplot2::labs(x = "Prediction", y = "Observed",
                  title = sprintf("Observed vs predicted (method = %s)", fit$method)) +
    theme_dMod(base_size = 11) +
    ggplot2::theme(legend.position = "none")
}



#' Weighted-residual diagnostics for an EM
#'
#' @description Two-panel scatter: IWRES vs IPRED and IWRES vs TIME with a
#'   loess smoother. Sourced from `plotResiduals(fit, ...)` when the first
#'   argument inherits from `EM` (internal type dispatch, see
#'   [plotResiduals]).
#' @param fit An [EM].
#' @param ... Ignored.
#' @return A ggplot.
#' @keywords internal
plotResiduals.em <- function(fit, ...) {
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- rbind(
    data.frame(x = pf$IPRED, y = pf$IWRES, name = pf$name,
               panel = "IWRES vs IPRED", stringsAsFactors = FALSE),
    data.frame(x = pf$time,  y = pf$IWRES, name = pf$name,
               panel = "IWRES vs TIME",  stringsAsFactors = FALSE))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, color = name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                         color = dMod_colors[2], linewidth = 0.6)
  if (length(unique(long$name)) > 1L) {
    p <- p + ggplot2::facet_grid(name ~ panel, scales = "free_x")
  } else {
    p <- p + ggplot2::facet_wrap(~ panel, scales = "free_x")
  }
  p +
    scale_color_dMod() +
    ggplot2::labs(x = NULL, y = "IWRES",
                  title = "Weighted residual diagnostics") +
    theme_dMod(base_size = 11)
}



#' Random-effect distribution diagnostics
#'
#' @description Per-eta histogram against the estimated `N(0, Omega_kk)`
#'   density plus a QQ-plot against the estimated normal. Detects systematic
#'   shrinkage, bimodality, or distributional misfit. Generic to leave room
#'   for non-EM methods in the future.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot (or a list of ggplots if `cowplot` is unavailable).
#' @export
plotHistIndivs <- function(x, ...) UseMethod("plotHistIndivs", x)


#' @rdname plotHistIndivs
#' @param x An [EM].
#' @param ... Ignored.
#' @return A ggplot (or a list with `hist` and `qq` if cowplot is unavailable).
#' @export
plotHistIndivs.em <- function(x, ...) {
  .emRequireOmega(x, "plotHistIndivs")
  fit <- x
  etaModes <- fit$etaModes
  om <- fit$omega
  if (is.null(etaModes) || is.null(om))
    stop("plotHistIndivs..fitNormal: fit has no etaModes or omega.")
  Omega <- if (!is.null(fit$Omega)) fit$Omega else {
    L <- om$buildL(fit$argument[om$cholPars])
    tcrossprod(L)
  }
  K <- ncol(etaModes)
  eta_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k <- sqrt(Omega[k, k])
    data.frame(eta_name = om$eta[k],
               value    = etaModes[, k],
               sd_est   = sd_k,
               stringsAsFactors = FALSE)
  }))
  # Per-eta Gaussian curve on a wide grid so the full N(0, Omega_kk) shape is
  # visible; using stat_function with a single mean(sd_est) painted all panels
  # with the same density (wrong when omega varies across etas) and was
  # clipped to the data extent by facet_wrap(scales = "free"). geom_line on
  # the explicit grid both picks up the panel-specific SD and widens the
  # x-range to ~ +/- 3.5 sigma.
  dens_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k  <- sqrt(Omega[k, k])
    xlim  <- max(abs(etaModes[, k]), 3.5 * sd_k)
    xseq  <- seq(-xlim, xlim, length.out = 200L)
    data.frame(eta_name = om$eta[k],
               x        = xseq,
               density  = dnorm(xseq, 0, sd_k),
               stringsAsFactors = FALSE)
  }))
  p_hist <- ggplot2::ggplot(eta_long, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 12, fill = "grey80", color = "grey40") +
    ggplot2::geom_line(data = dens_long,
                       ggplot2::aes(x = x, y = density),
                       inherit.aes = FALSE,
                       color = dMod_colors[2], linewidth = 0.8) +
    ggplot2::facet_wrap(~ eta_name, scales = "free") +
    ggplot2::labs(x = "eta value", y = "density",
                  title = "Eta distribution vs N(0, Omega_kk)") +
    theme_dMod(base_size = 11)
  p_qq <- ggplot2::ggplot(eta_long,
                          ggplot2::aes(sample = value / sd_est)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = dMod_colors[2]) +
    ggplot2::facet_wrap(~ eta_name) +
    ggplot2::labs(x = "Theoretical quantile (N(0,1))",
                  y = "Standardised eta",
                  title = "Eta QQ vs N(0,1)") +
    theme_dMod(base_size = 11)
  if (requireNamespace("cowplot", quietly = TRUE))
    cowplot::plot_grid(p_hist, p_qq, ncol = 1)
  else list(hist = p_hist, qq = p_qq)
}



#' ECM convergence trace (OFV, |delta-psi|) per stage
#'
#' @description Four-panel trace of OFV, structural-parameter step
#'   `|delta psi|`, max softmax weight, and minimum effective node count
#'   across ECM iterations. Quadrature-method EM only.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot.
#' @export
plotTrace <- function(x, ...) UseMethod("plotTrace", x)


#' @rdname plotTrace
#' @param x An [EM] with `method = "quadrature"`. Errors otherwise.
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plotTrace.em <- function(x, ...) {
  .emRequireOmega(x, "plotTrace")
  fit <- x
  if (!fit$method %in% c("quadrature", "foceiQuadrature") ||
      is.null(fit$stageTrace))
    stop("plotTrace requires an .fitNormal fit with method 'quadrature' or 'foceiQuadrature'.")
  tr <- fit$stageTrace
  tr$iter <- seq_len(nrow(tr))
  long <- rbind(
    data.frame(iter = tr$iter, value = tr$OFV,        panel = "OFV",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$deltaPsi,  panel = "|delta psi|",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$maxSoftmax, panel = "max softmax",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$nEffMin,   panel = "n_eff min",
               level = factor(tr$level)))
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value, color = level)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    scale_color_dMod() +
    ggplot2::labs(x = "ECM iteration",
                  title = "ECM convergence trace") +
    theme_dMod(base_size = 11)
}


# Long-format samples for ggplot.
.bayes_samples_long <- function(samples, par_names = NULL) {
  if (!is.matrix(samples)) stop("samples must be a matrix.")
  if (is.null(colnames(samples))) {
    if (is.null(par_names))
      par_names <- paste0("p", seq_len(ncol(samples)))
    colnames(samples) <- par_names
  }
  N <- nrow(samples)
  data.frame(
    sample    = rep(seq_len(N), times = ncol(samples)),
    parameter = factor(rep(colnames(samples), each = N),
                       levels = colnames(samples)),
    value     = as.vector(samples),
    stringsAsFactors = FALSE)
}


#' Plot the marginal posterior distributions of an [mcmc()] result
#'
#' One panel per parameter: histogram of the post-warmup samples overlaid
#' with a kernel-density curve. The vertical dashed line marks the
#' posterior mean.
#'
#' @param x An object inheriting from `mcmcResult`.
#' @param parameters Optional character vector restricting which parameter
#'   panels are drawn.
#' @param bins Histogram bin count. Default 30.
#' @param ... Ignored.
#'
#' @return A ggplot.
#' @export
plot.mcmcresult <- function(x, parameters = NULL, bins = 30L, ...) {
  long <- .bayes_samples_long(x$samples)
  if (!is.null(parameters)) {
    bad <- setdiff(parameters, levels(long$parameter))
    if (length(bad))
      stop("Unknown parameter name(s): ", paste(bad, collapse = ", "))
    long <- long[long$parameter %in% parameters, , drop = FALSE]
    long$parameter <- droplevels(long$parameter)
  }
  means <- stats::aggregate(value ~ parameter, data = long, FUN = mean)
  ggplot2::ggplot(long, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = "grey80", color = "grey40") +
    ggplot2::geom_density(color = dMod_colors[3], linewidth = 0.7) +
    ggplot2::geom_vline(data = means,
                        ggplot2::aes(xintercept = value),
                        color = dMod_colors[2], linetype = "dashed",
                        linewidth = 0.6) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::labs(x = NULL, y = "density",
                  title = "Marginal posterior") +
    theme_dMod(base_size = 11)
}


#' Multi-chain density overlay with R-hat annotation
#'
#' @param x A `mcmcResultMulti` object.
#' @param parameters Optional character vector restricting parameters.
#' @param bins Histogram bin count (currently unused).
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.mcmcresultmulti <- function(x, parameters = NULL, bins = 30L, ...) {
  S <- x$samples
  ch <- factor(x$chainId)
  if (!is.null(parameters)) S <- S[, parameters, drop = FALSE]
  par_names <- colnames(S)
  long <- do.call(rbind, lapply(seq_along(par_names), function(j) {
    data.frame(parameter = par_names[j], chain = ch, value = S[, j],
               stringsAsFactors = FALSE)
  }))
  long$parameter <- factor(long$parameter, levels = par_names)

  finite_rhat <- if (!is.null(x$rHat)) x$rHat[is.finite(x$rHat)] else numeric()
  rhat_label <- if (length(finite_rhat))
    sprintf("Marginal posterior (%d chains, max R-hat = %.3f)",
            x$nChains, max(finite_rhat))
  else
    sprintf("Marginal posterior (%d chains)", x$nChains)

  ggplot2::ggplot(long, ggplot2::aes(x = value,
                                       color = .data$chain,
                                       fill  = .data$chain)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.6) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    scale_color_dMod() +
    ggplot2::scale_fill_manual(values = dMod_palette(nlevels(ch))) +
    ggplot2::labs(x = NULL, y = "density", title = rhat_label) +
    theme_dMod(base_size = 11)
}


#' @rdname plotTrace
#' @description Method for [mcmc()] sequential / SMC outputs: shows the
#'   adaptive tempering schedule (beta vs level), the post-resample ESS,
#'   and per-level acceptance.
#' @export
plotTrace.mcmcresultsequential <- function(x, ...) {
  diag_df <- rbind(
    data.frame(level = seq_along(x$betaPath) - 1L, value = x$betaPath,
               panel = "beta"),
    data.frame(level = seq_along(x$ESSPath) - 1L,  value = x$ESSPath,
               panel = "ESS"))
  if (length(x$acceptRates))
    diag_df <- rbind(diag_df,
                     data.frame(level = seq_along(x$acceptRates),
                                value = x$acceptRates,
                                panel = "accept rate"))
  ggplot2::ggplot(diag_df, ggplot2::aes(x = level, y = value)) +
    ggplot2::geom_line(color = dMod_colors[3], linewidth = 0.6) +
    ggplot2::geom_point(color = dMod_colors[1], size = 1.2) +
    ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 1L) +
    ggplot2::labs(x = "SMC level",
                  title = sprintf("SMC trace (log-evidence = %s)",
                                  format(x$logEvidence, digits = 4))) +
    theme_dMod(base_size = 11)
}


#' @rdname plotTrace
#' @description Method for [mcmc()] single-chain outputs: per-parameter
#'   trace of the chain across iterations.
#' @param parameters Optional character vector restricting which parameter
#'   panels are drawn.
#' @export
plotTrace.mcmcresultsingle <- function(x, parameters = NULL, ...) {
  S <- x$samples
  if (!is.null(parameters)) S <- S[, parameters, drop = FALSE]
  long <- .bayes_samples_long(S)
  long$iter <- long$sample
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value)) +
    ggplot2::geom_line(color = dMod_colors[3], linewidth = 0.4) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(x = "iteration", y = NULL, title = "Chain trace") +
    theme_dMod(base_size = 11)
}


#' @rdname plotTrace
#' @export
plotTrace.mcmcresultblocked <- function(x, parameters = NULL, ...) {
  plotTrace.mcmcresultsingle(x, parameters = parameters, ...)
}


#' Pair (corner-style) plot for [mcmc()] outputs
#'
#' Lower triangle: 2D scatter of post-warmup samples. Diagonal: 1D
#' marginal density.
#'
#' @param x An `mcmcResult` (or subclass) object.
#' @param parameters Optional character vector restricting parameters.
#' @param ... Method-specific arguments.
#' @return A ggplot.
#' @export
plotPairs <- function(x, ...) UseMethod("plotPairs", x)


#' @rdname plotPairs
#' @export
plotPairs.mcmcresult <- function(x, parameters = NULL, ...) {
  S <- x$samples
  if (!is.null(parameters)) S <- S[, parameters, drop = FALSE]
  par_names <- colnames(S)
  P <- ncol(S)
  if (P < 2L) stop("plotPairs needs >= 2 parameters.")

  combos <- expand.grid(xPar = par_names, yPar = par_names,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  combos$panel <- with(combos, ifelse(xPar == yPar, "diag",
                                       ifelse(match(yPar, par_names) >
                                              match(xPar, par_names),
                                              "lower", "upper")))
  combos <- combos[combos$panel %in% c("diag", "lower"), , drop = FALSE]
  long <- do.call(rbind, lapply(seq_len(nrow(combos)), function(idx) {
    xp <- combos$xPar[idx]; yp <- combos$yPar[idx]
    data.frame(xPar  = factor(xp, levels = par_names),
               yPar  = factor(yp, levels = par_names),
               x     = S[, xp],
               y     = if (combos$panel[idx] == "diag") NA_real_ else S[, yp],
               panel = combos$panel[idx], stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot()
  diag_df <- long[long$panel == "diag",  , drop = FALSE]
  off_df  <- long[long$panel == "lower", , drop = FALSE]
  if (nrow(off_df) > 0L)
    p <- p + ggplot2::geom_point(data = off_df,
                                  ggplot2::aes(x = x, y = y),
                                  alpha = 0.3, size = 0.5,
                                  color = dMod_colors[3])
  if (nrow(diag_df) > 0L)
    p <- p + ggplot2::geom_density(data = diag_df,
                                    ggplot2::aes(x = x),
                                    color = dMod_colors[2],
                                    linewidth = 0.7,
                                    inherit.aes = FALSE)
  p +
    ggplot2::facet_grid(yPar ~ xPar, scales = "free", switch = "both") +
    ggplot2::labs(x = NULL, y = NULL, title = "Posterior pair-plot") +
    theme_dMod(base_size = 10) +
    ggplot2::theme(strip.placement = "outside",
                   panel.spacing   = ggplot2::unit(0.1, "lines"))
}


## ---- profile / parameter-path plotting (moved from toolsSvenja.R) ---------
#' Plot an array of trajectories along the profile of a parameter
#' 
#' @param par Character of parameter name for which the array should be generated.
#' @param profs Lists of profiles as being returned by [profile]. 
#' @param prd Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param times Numeric vector of time points for the model prediction.
#' @param direction Character "up" or "down" indicating the direction the value should be traced along the profile starting at the bestfit value.
#' @param covtable Optional covariate table or condition.grid necessary if subsetting is required.
#' @param ... Further arguments for subsetting the plot.
#' @param nsimus Number of trajectories/ simulation to be calculated.
#' 
#' @return A plot object of class `ggplot`.
#' @author Svenja Kemmer, \email{svenja.kemmer@@fdm.uni-freiburg.de}
#' @examples
#' \dontrun{
#'  plotArray("myparameter", myprofiles, g*x*p, seq(0, 250, 1), 
#'     "up", condition.grid, name == "ProteinA" & condition == "c1") 
#' }
#' @export
#' @import data.table
plotArray <- function (par, profs, prd, times, direction = c("up", "down"), covtable, ..., nsimus = 4) {
  
  # select subframe from profiles
  mysub <- profs %>% as.data.table() %>% .[whichPar == par, ]
  mysub[, ID := 1:nrow(mysub)]
  
  # get ID of bestfit (constraint is 0 for bestfit)
  bestID <- mysub[constraint == 0.00]$ID
  if(direction == "up") mysubF <- mysub[ID >= bestID]  
  if(direction == "down") mysubF <- mysub[ID <= bestID]
  
  # select rows according to simulation number
  partable <- mysubF[seq(1, nrow(mysubF), (round(nrow(mysubF)/nsimus)))]
  
  # remove non_parameter names
  no_pars <- c("value", "constraint", "stepsize", "gamma", "whichPar", "data", "condition_obj", "AIC", "BIC", "prior", "ID", "chisquare")
  partable %>% .[, (no_pars) := NULL]
  
  # make predictions
  predictionDT <- .predictArray(prd, times, pars = partable, whichpar = par)
  out_plot <- copy(predictionDT)
  
  # use covtable for subsetting of the plot
  if(!is.null(covtable)) {
    if(!"condition" %in% names(covtable)){
      covtable <- as.data.table(covtable, keep.rownames = "condition")
    } else covtable <- as.data.table(covtable)
    out_plot <- merge(out_plot, covtable, by = "condition")
    out_plot <- out_plot[...]
  }
  
  # plot
  P <- ggplot(out_plot , aes(x = time, y = value, group = ParValue, color = ParValue)) +
    facet_grid(name~condition, scales = "free_y") +
    geom_line(size = 1) + 
    theme_dMod(base_size = 18) + scale_color_viridis_c() +
    theme(legend.position = "top", legend.key.size = unit(0.6,"cm")) + 
    theme(axis.line = element_line(colour = "black"), 
          panel.grid.major = element_line(colour = "grey97"), 
          panel.grid.minor = element_line(colour = "grey97"), 
          panel.background = element_blank()) +
    xlab("time") +
    ylab(paste0("value"))
  
  return(P)
}

.predictArray <- function (prd, times, pars = partable, whichpar = par, keep_names = NULL, FLAGverbose = FALSE, FLAGverbose2 = FALSE, FLAGbrowser = FALSE, ...) {
  .require_ns("purrr", ".predictArray()")
  if (FLAGverbose2) cat("Simulating", "\n")
  out <- lapply(1:nrow(pars), function(i) {
    if (FLAGverbose) cat("Parameter set", i, "\n")
    if (FLAGbrowser) browser()
    mypar <- pars[i,] %>% as.numeric()
    parval <- round(pars[i,][[whichpar]], digits = 2)
    names(mypar) <- names(pars)
    mypar <- as.parvec(mypar)
    prediction <- try(prd(times, mypar, deriv = FALSE, ...))
    if (inherits(prediction, "try-error")) {
      warning("parameter set ", i, " failed\n")
      return(NULL)
    }
    prediction <- purrr::imap(prediction, function(.x,.y){
      .x <- data.table(.x)
      if (!is.null(keep_names))
        .x[, (setdiff(names(.x), c(keep_names, "time"))) := NULL]
      .x[, `:=`(condition = .y, ParValue = parval)]
      .x
    })
    melt(rbindlist(prediction), variable.name = "name", value.name = "value", id.vars = c("time", "condition", "ParValue"))
  })
  if (FLAGverbose2) cat("postprocessing", "\n")
  out <- rbindlist(out[!is.null(out)])
  out
}


.findEmptyCorner <- function(x, y) {
  xmid <- (min(x, na.rm = TRUE) + max(x, na.rm = TRUE)) / 2
  ymid <- (min(y, na.rm = TRUE) + max(y, na.rm = TRUE)) / 2
  
  corners <- list(
    bottom_left  = c(0.05, 0.05),
    bottom_right = c(0.95, 0.05),
    top_left     = c(0.05, 0.95),
    top_right    = c(0.95, 0.95)
  )
  
  counts <- c(
    bottom_left  = sum(x <= xmid & y <= ymid, na.rm = TRUE),
    bottom_right = sum(x >  xmid & y <= ymid, na.rm = TRUE),
    top_left     = sum(x <= xmid & y >  ymid, na.rm = TRUE),
    top_right    = sum(x >  xmid & y >  ymid, na.rm = TRUE)
  )
  
  corners[[which.min(counts)]]
}

#' @keywords internal
#' @importFrom ggplot2 ggplot
PlotPaths <- function(profs=myprofiles, ..., whichPar, sort = FALSE, relative = TRUE, scales = "fixed", multi = TRUE, n_pars = 5, normalizePaths = FALSE) {
  
  if ("parframe" %in% class(profs)) {
    arglist <- list(profs)
  } else {
    arglist <- as.list(profs)
  }
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- as.data.frame(arglist[[i]])
    parameters <- attr(arglist[[i]], "parameters")
    
    if (is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    if (is.null(whichPar)) whichPar <- names(proflist)
    if (is.numeric(whichPar)) whichPar <- names(proflist)[whichPar]
    
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # matrix
      paths <- as.matrix(proflist[[n]][, parameters])
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      
      # Save absolute values of profiled parameter before relativizing
      abs_profiled <- as.numeric(paths[, n])
      
      if (relative) 
        for(j in 1:ncol(paths)) paths[, j] <- as.numeric(paths[, j]) - as.numeric(paths[origin, j])
      
      # Restore absolute values for the profiled parameter (x-axis always absolute)
      paths[, n] <- abs_profiled
      
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if (sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- submatrix(combinations, cols = -which(combinations[1,] == combinations[2,]))
      combinations <- submatrix(combinations, cols = !duplicated(paste(combinations[1,], combinations[2,])))
      
      
      
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = profnames[i],
                   combination = paste(combinations[,j], collapse = " - \n"),
                   x = paths[, combinations[1,j]],
                   y = paths[, combinations[2,j]])
      }))
      
      if(multi) path.data <- path.data %>% as.data.table %>% .[, partner := tstrsplit(as.character(combination), "\n", fixed=TRUE, keep = 2)]
      
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  if (relative){
    axis.labels <- c("parameter 1", expression(Delta ~ p[j]))
  } else {
    axis.labels <- c("parameter 1", "parameter 2")
  }
  
  data <- droplevels(subset(data, ...))
  removeBecauseNonsense <- c("value", "constraint", "stepsize", "chisquare", "data", "prior", "gamma", "whichPar")
  data <- data[!(partner %in% removeBecauseNonsense)]
  data$y <- as.numeric(data$y)
  data$x <- as.numeric(data$x)
  
  if (normalizePaths == TRUE) {
    data[, y := (ifelse(max(abs(y)) == 0, 0, y / abs(max(abs(y))))), by = combination] # if path is y, just return 0
    # data[, y := (2 * (y - min(y)) / (max(y) - min(y))) - 1, by = combination]
    removedCombinations <- unique(data[!is.finite(y), combination])
    data <- data[is.finite(y)]
    
    if(length(removedCombinations)>0) {warning(paste0("The following combinations have been removed due to failed paths:\n\t",paste(str_remove_all(removedCombinations, "\n"), collapse = "\n\t")))}
  }
  
  
  if(multi){
    
    # determine strength of change
    data[, max.dev := max(c(abs(max(as.numeric(y))), abs(min(as.numeric(y) )))), by = "partner"]
    setorder(data, name, -max.dev)
    # max.devis <- unique(data$max.dev)[1:n_pars]
    
    # create new column "label" only use to assign ploting colors
    data[,label := ifelse(max.dev %in% unique(max.dev)[1:n_pars], partner, "Others")]
    
    # Define the plotting colors
    species_colors <- c(
      setNames(dMod_palette(n_pars + 1L)[-1L], unique(data$partner)[1:n_pars]),
      "Others" = "gray"
    )
    
    # Automatically find the corner with the least data density
    legend_corner <- .findEmptyCorner(data$x, data$y)
    
    suppressMessages(
      p <- ggplot2::ggplot(data, aes(x = x, y = y, color = label, group = partner)) + 
        geom_line() +
        xlab(whichPar) + ylab(expression(Delta ~ p[j])) +
        scale_linetype_discrete(name = "profile\nlist") +
        scale_color_manual(values = species_colors) + theme_dMod() +
        theme(legend.position = "inside",
              legend.position.inside = legend_corner,
              legend.justification = legend_corner,
              legend.title = element_blank(),
              legend.background = element_rect(fill = alpha("white", 0.85), colour = "black", linewidth = 0.3),
              legend.key.size = unit(0.4, "cm"),
              legend.margin = margin(2, 4, 2, 4),
              legend.text = element_text(size = 7))
    )
  } else {
    suppressMessages(
      p <- ggplot2::ggplot(data, aes(x = x, y = y, group = interaction(name, proflist), color = name, lty = proflist)) + 
        facet_wrap(~combination, scales = scales) + 
        geom_path() + #geom_point(aes=aes(size=1), alpha=1/3) +
        xlab(axis.labels[1]) + ylab(axis.labels[2]) +
        scale_linetype_discrete(name = "profile\nlist") +
        scale_color_dMod(name = "profiled\nparameter")
    )
  }
  
  attr(p, "data") <- data
  return(p)
  
}

#' Profile likelihood: plot all parameter paths belonging to one profile in one plot
#' 
#' @param profs Lists of profiles as being returned by [profile]. 
#' @param whichpars Character vector of parameter names for which the profile paths should be generated.
#' @param npars Numeric vector of number of colored and named parameter paths.
#' @param normalizePaths Logical indicating whether the paths should be normalized to absolute values of 1. Default `FALSE`; `TRUE` only useful in corner cases when you know why to do so.
#' 
#' @return A plot object of class `ggplot` for length(whichpars) = 1 and otherwise an object of class `cowplot`.
#' @author Svenja Kemmer, \email{svenja.kemmer@@fdm.uni-freiburg.de}
#' @examples
#' \dontrun{
#'  plotPathsMulti(myprofiles, c("mypar1", "mypar2"), npars = 5) 
#' }
#' @export
#' @import data.table
plotPathsMulti <- function(profs, whichpars, npars = 5, normalizePaths = FALSE) {
  .require_ns("cowplot", "plotPathsMulti()")
  if(length(whichpars) == 1){
    p <- PlotPaths(profs=profs, whichPar = whichpars, n_pars = npars, normalizePaths = normalizePaths)
    return(p)
  } else {
    PlotList <- NULL
    for(i in 1:length(whichpars)){
      par <- whichpars[i]
      p <- PlotPaths(profs=profs, whichPar = par, n_pars = npars, normalizePaths = normalizePaths)
      PlotList[[i]] <- p
    }
    pl <- cowplot::plot_grid(plotlist = PlotList)
    return(pl)
  }
}



#' Profile likelihood: plot profiles along with their parameter paths
#' 
#' Generates combined plots of profile likelihoods and their parameter paths without a shared legend.
#' 
#' @param profs List of profiles as returned by [profile()].
#' @param whichpars Character vector of parameter names for which the profile paths should be generated.
#' @param npars Numeric indicating number of colored and named parameter paths.
#' @param ncols Number of columns in the resulting plot grid.
#' @param normalizePaths Logical indicating whether the paths should be normalized to absolute values of 1.
#'                       Default `FALSE`.
#' @param modes Character vector of profile modes to display in the profile plot.
#'              Default `c("data", "prior")`. Use e.g. `"data"` to show only the data contribution.
#' @param ... Additional arguments passed to `cowplot::plot_grid()`.
#' 
#' @return A combined `ggplot` object containing the profiles and paths (no shared legend).
#' 
#' @export
plotProfilesAndPaths <- function(profs, whichpars, npars = 5, ncols = 3, normalizePaths = FALSE, modes = c("data", "prior"), ...) {
  .require_ns("cowplot", "plotProfilesAndPaths()")

  # Save original obj.attributes before any subsetting drops them
  orig_oa <- attr(profs, "obj.attributes")
  filtered_oa <- if (!is.null(orig_oa)) intersect(orig_oa, modes) else NULL
  
  profs <- profs[profs$whichPar %in% whichpars]
  
  cleanProfilePlot <- function(prof_sub) {
    # Remove columns for modes we don't want, so plotProfile can't plot them
    cols_to_drop <- setdiff(orig_oa, modes)
    prof_sub <- prof_sub[, !(colnames(prof_sub) %in% cols_to_drop), drop = FALSE]
    attr(prof_sub, "obj.attributes") <- filtered_oa
    p <- plotProfile(prof_sub)
    
    # plotProfile always adds "total"; filter the underlying data to requested modes
    pdata <- attr(p, "data")
    pdata <- pdata[pdata$mode %in% modes, , drop = FALSE]
    
    threshold <- c(1, 2.7, 3.84)
    p_new <- ggplot(pdata, aes(x = par, y = delta, group = interaction(proflist, mode), 
                               color = proflist, linetype = mode)) +
      facet_wrap(~name, scales = "free_x") +
      geom_hline(yintercept = threshold, lty = 2, color = "gray") +
      geom_line() +
      geom_point(data = subset(pdata, is.zero)) +
      ylab(expression(paste("CL /", Delta * chi^2))) +
      scale_y_continuous(breaks = c(1, 2.7, 3.84), 
                         labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"),
                         limits = c(NA, 5)) +
      xlab("parameter value") +
      labs(title = NULL, x = NULL, linetype = "contrib") +
      theme(
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      guides(color = "none", fill = "none") +
      theme(legend.position = "none")
  }
  
  stacked_list <- vector("list", length(whichpars))
  
  for (z in seq_along(whichpars)) {
    prof_sub <- profs[profs$whichPar == whichpars[z]]
    
    p_prof_noleg <- cleanProfilePlot(prof_sub)
    
    p_paths <- plotPathsMulti(prof_sub, whichpars[z], npars, normalizePaths = normalizePaths) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    aligned_pair <- cowplot::align_plots(p_prof_noleg, p_paths, align = "v", axis = "tb")
    stacked_list[[z]] <- cowplot::plot_grid(aligned_pair[[1]], aligned_pair[[2]],
                                            ncol = 1, rel_heights = c(1, 0.7), align = "v", axis = "tb")
  }
  
  body <- cowplot::plot_grid(plotlist = stacked_list, ncol = ncols, ...)
  
  return(body)
}