
# Custom interface to ggplot2 ---

#' Open last plot in external pdf viewer
#' 
#' @description Convenience function to show last plot in an external viewer.
#' @param plot \code{ggplot2} plot object.
#' @param command character, indicatig which pdf viewer is started.
#' @param ... arguments going to \code{ggsave}.
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
  
  # theme_bw(base_size = base_size, base_family = base_family) + 
  #   theme(line = element_line(colour = black), 
  #         rect = element_rect(fill = "white", colour = NA), 
  #         text = element_text(colour = black), 
  #         axis.ticks = element_line(colour = black), 
  #         axis.text = element_text(color = black),
  #         legend.key = element_rect(colour = NA), 
  #         panel.border = element_rect(colour = black), 
  #         panel.grid = element_line(colour = "gray90", size = 0.2), 
  #         #panel.grid = element_blank(), 
  #         strip.background = element_rect(fill = "white", colour = NA))
  
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = "black"),
          axis.text = element_text(size = rel(1.0), colour = "black"),
          axis.text.x = element_text(margin=unit(c(4, 4, 0, 4), "mm")),
          axis.text.y = element_text(margin=unit(c(4, 4, 4, 0), "mm")),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = "black"),
          # panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(1.0)))
  
}

dMod_colors <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E", "#4B1F6F", "#CC79A7","#006400", "#F0E442", "#8B4513", rep("gray", 100))

#' Standard dMod color palette
#' 
#' @param ... arguments goint to code{scale_color_manual()}
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
  scale_color_manual(..., values = dMod_colors)
}


#' Standard dMod color scheme
#' 
#' @export
#' @param ... arguments goint to code{scale_color_manual()}
scale_fill_dMod <- function(...) {
  scale_fill_manual(..., values = dMod_colors)
}

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_dMod() + theme_dMod()




# Other ---------------------------------------------

#' Coordinate transformation for data frames
#' 
#' Applies a symbolically defined transformation to the \code{value}
#' column of a data frame. Additionally, if a \code{sigma} column is
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
#' as generated by \link{Xs}.
#' @param ... Further arguments going to \code{dplyr::filter}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @param transform list of transformation for the states, see \link{coordTransform}.
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @import ggplot2
#' @example inst/examples/plotting.R
#' @export
plotPrediction <- function(prediction,...) {
  UseMethod("plotPrediction", prediction)
}




#' Plot a list of model predictions and a list of data points in a combined plot
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by \link{Xs}.
#' @param data Named list of data.frames as being used in \link{res}, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param ... Further arguments going to \code{dplyr::filter}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet \code{"wrap"} or \code{"grid"}. Try \code{"wrap_plain"} for high amounts of conditions and low amounts of observables.
#' @param transform list of transformation for the states, see \link{coordTransform}.
#' @param aesthetics Named list of aesthetic mappings, specified as character, e.g. \code{list(linetype = "name")}. 
#' Can refer to variables in the condition.grid
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{sigma},
#' \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @example inst/examples/plotting.R
#' @importFrom graphics par
#' @export
plotCombined <- function(prediction,...) {
  UseMethod("plotCombined", prediction)
}


#' Plot a list data points
#' 
#' @param data Named list of data.frames as being used in \link{res}, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param ... Further arguments going to \code{subset}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @param transform list of transformation for the states, see \link{coordTransform}.
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{sigma},
#' \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
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
#' @param profs Lists of profiles as being returned by \link{profile}.
#' @param ... logical going to subset before plotting.
#' @param maxvalue Numeric, the value where profiles are cut off.
#' @param parlist Matrix or data.frame with columns for the parameters to be added to the plot as points.
#' If a "value" column is contained, deltas are calculated with respect to lowest chisquare of profiles.
#' @return A plot object of class \code{ggplot}.
#' @details See \link{profile} for examples.
#' @export
plotProfile <- function(profs,...) {
  UseMethod("plotProfile", profs)
}


#' Profile likelihood: plot of the parameter paths.
#' 
#' @param profs profile or list of profiles as being returned by \link{profile}
#' @param ... arguments going to subset
#' @param whichPar Character or index vector, indicating the parameters that are taken as possible reference (x-axis)
#' @param sort Logical. If paths from different parameter profiles are plotted together, possible
#' combinations are either sorted or all combinations are taken as they are.
#' @param relative logical indicating whether the origin should be shifted.
#' @param scales character, either \code{"free"} or \code{"fixed"}.
#' @return A plot object of class \code{ggplot}.
#' @details See \link{profile} for examples.
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
      scale_color_manual(name = "profiled\nparameter", values = dMod_colors)
  )
  
  attr(p, "data") <- data
  return(p)
  
}



#' Plot Fluxes given a list of flux Equations
#' 
#' @param pouter parameters
#' @param x The model prediction function \code{x(times, pouter, fixed, ...)} 
#' @param fluxEquations list of chars containing expressions for the fluxes, 
#' if names are given, they are shown in the legend. Easy to obtain via \link{subset.eqnlist}, see Examples.
#' @param nameFlux character, name of the legend.
#' @param times Numeric vector of time points for the model prediction
#' @param ... Further arguments going to x, such as \code{fixed} or \code{conditions}
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @examples 
#' \dontrun{
#' 
#' plotFluxes(bestfit, x, times, subset(f, "B"%in%Product)$rates, nameFlux = "B production")
#' }
#' @export
plotFluxes <- function(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", ...){
  
  if (is.null(names(fluxEquations))) names(fluxEquations) <- fluxEquations
  
  flux <- funC0(fluxEquations, convenient = FALSE)
  prediction.all <- x(times, pouter, deriv = FALSE, ...)
  names.prediction.all <- names(prediction.all)
  if (is.null(names.prediction.all)) names.prediction.all <- paste0("C", 1:length(prediction.all))
  
  out <- lapply(1:length(prediction.all), function(cond) {
    prediction <- prediction.all[[cond]]
    pinner <- attr(prediction, "parameters")
    pinner.matrix <- matrix(pinner, nrow = length(pinner), ncol = nrow(prediction), 
                            dimnames = list(names(pinner), NULL))
    #pinner.list <- as.list(pinner)
    #prediction.list <- as.list(as.data.frame(prediction))
    fluxes <- cbind(time = prediction[, "time"], flux(cbind(prediction, t(pinner.matrix))))
    return(fluxes)
  }); names(out) <- names.prediction.all
  out <- wide2long(out)
  
  cbPalette <- c("#999999", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2", 
                 "#D55E00", "#CC79A7","#CC6666", "#9999CC", "#66CC99","red", "blue", "green","black")
  
  P <- ggplot(out, aes(x = time, y = value, group = name, fill = name, log = "y")) + 
    facet_wrap(~condition) + scale_fill_manual(values = cbPalette, name = nameFlux) +
    geom_density(stat = "identity", position = "stack", alpha = 0.3, color = "darkgrey", size = 0.4) +
    xlab("time") + ylab("flux contribution")
   
  attr(P, "out") <- out
  
  return(P)
  
}


stepDetect <- function(x, tol) {
  
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
#' a \link{parframe}.
#' @param ... arguments for subsetting of x
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
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
#' @param parframe Object of class \code{parframe}, e.g. returned by \link{mstrust}
#' @param x Prediction function returning named list of data.frames with names as \code{data}.
#' @param data Named list of data.frames, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param split List of characters specifying how to summarise the residuals by \code{sqrt(res_i^2)}, 
#' \code{split[1]} used for x-axis, \code{split[2]} for grouping (color), and any additional for \code{facet_wrap()}
#' @param ... Additional arguments for x
#' @param errmodel object of type prdfn, the error model function.
#' 
#' @return A plot object of class \code{ggplot} with data.frame as attribute \code{attr(P,"out")}.
#' @examples
#' \dontrun{
#'  # time axis:
#'  plotResiduals(myfitlist, g*x*p, data, 
#'     c("time","index","condition","name"), 
#'     conditions = myconditions[1:4])
#'  # condition axis (residuals summed over time for each observable and condition):
#'  plotResiduals(myfitlist, g*x*p, data,  c("condition","name","index"))
#' }
#' @export
#' @importFrom plyr ddply summarise
plotResiduals <- function(parframe, x, data, split = "condition", errmodel = NULL, ...){
  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if(!("index" %in% colnames(parframe)))
    parframe$index <- 1:nrow(parframe)
  
  out <- do.call(rbind,lapply(1:nrow(parframe), function(j){
    pred <- x(timesD, as.parvec(parframe,j), deriv = FALSE, ...)
    
    out_con <- do.call(rbind,lapply(names(pred), function(con){
      err <- NULL
      if (!is.null(errmodel)) {
        err <- errmodel(out = pred[[con]], pars = getParameters(pred[[con]]), conditions = con)
      }
      out <- res(data[[con]], pred[[con]], err[[con]]) 
      return(cbind(out,condition = con))
    })
    )
    
    out_par <- cbind(index = as.character(parframe[j,"index"]), out_con)
    return(out_par)
  })
  )
  if (!is.null(errmodel)) {
    out <- plyr::ddply(out, split, summarise, res = sum(weighted.residual^2 + log(sigma^2))) 
  } else{
    out <- plyr::ddply(out, split, summarise, res = sum(weighted.residual^2)) 
  }
  groupvar <- split[1]
  if(length(split) > 1){
    groupvar <- split[2]
  }
  
  P <- ggplot2::ggplot(out, aes_string(x = split[1], y = "res", color = groupvar, group = groupvar)) + theme_dMod() + geom_point() + geom_line()
  if(length(split) > 2)
    P <- P + facet_wrap(split[3:length(split)]) 
  attr(P,"out") <- out
  return(P)
}
