#' plotFileScatters
#' 
#' Make a scatter plot per channel for all provided files
#'
#' @param  input      Either a flowSet, a flowFrame with a file ID column (e.g. 
#'                    output from the \code{\link{AggregateFlowFrames}} includes
#'                    a "File" column) or a vector of paths pointing to FCS files
#' @param  fileID     Name of the file ID column when the input is a flowFrame, 
#'                    default to "File" (File ID column in the 
#'                    \code{\link{AggregateFlowFrames}} flowFrame output).
#' @param  channels   Vector of channels or markers that need to be plotted, 
#'                    if NULL (default), all channels from the input will be 
#'                    plotted
#' @param  yLim       Optional vector of a lower and upper limit of the y-axis
#' @param  yLabel     Determines the label of the y-axis. Can be "marker" and\\or
#'                    "channel" or abbrevations. Default = "marker".
#' @param  quantiles  If provided (default NULL), a numeric vector with values
#'                    between 0 and 1. These quantiles are indicated on the plot
#' @param  names      Optional parameter to provide filenames. If \code{NULL} 
#'                    (default), the filenames will be numbers. Duplicated 
#'                    filenames will be made unique.
#' @param  groups     Optional parameter to specify groups of files, should have
#'                    the same length as the \code{input}. Id \code{NULL} 
#'                    (default), all files will be plotted in the same color
#' @param  color      Optional parameter to provide colors. Should have the same
#'                    lengths as the number of groups (or 1 if \code{groups} is 
#'                    \code{NULL})
#' @param  legend     Logical parameter to specify whether the group levels 
#'                    should be displayed. Default is \code{FALSE}
#' @param  maxPoints Total number of data points that will be plotted per 
#'                    channel, default is 50000
#' @param  silent     If FALSE, prints an update every time it starts processing 
#'                    a new file. Default = FALSE. 
#' @param  ncol       Number of columns in the final plot, optional
#' @param  nrow       Number of rows in the final plot, optional
#' @param  width      Width of png file. By default NULL the width parameter is 
#'                    estimated based on the input.
#' @param  height     Height of png file. By default NULL the width parameter is 
#'                    estimated based on the input.
#'  
#' @param  plotFile   Path to png file, default is "FileScatters.png". If 
#'                    \code{NULL}, the output will be a list of ggplots 
#' @param ...         Arguments for read.FCS (e.g. truncate_max_range)
#' 
#' @return List of ggplot objects if \code{plot} is \code{FALSE}, 
#'         otherwise \code{filePlot} with plot is created.
#'         
#' @examples 
#' # Preprocessing
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' 
#' flowCore::write.FCS(ff[1:1000, ], file = "ff_tmp1.fcs")
#' flowCore::write.FCS(ff[1001:2000, ], file = "ff_tmp2.fcs")
#' flowCore::write.FCS(ff[2001:3000, ], file = "ff_tmp3.fcs")
#'  
#' # Make plot
#' plotFileScatters(input = c("ff_tmp1.fcs", "ff_tmp2.fcs", "ff_tmp3.fcs"),
#'                  channels = c("Pacific Blue-A", 
#'                               "Alexa Fluor 700-A", 
#'                               "PE-Cy7-A"), 
#'                  maxPoints = 1000)
#' 
#' @import ggplot2
#' @importFrom methods is
#' @importFrom flowCore fsApply exprs
#' @importFrom dplyr tibble group_by summarise
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#' @importFrom stats quantile
#'  
#' @export 
plotFileScatters <- function(input, 
                             fileID = "File",
                             channels = NULL, 
                             yLim = NULL, 
                             yLabel = "marker",
                             quantiles = NULL,
                             names = NULL,
                             groups = NULL, 
                             color = NULL, 
                             legend = FALSE,
                             maxPoints = 50000, 
                             ncol = NULL, 
                             nrow = NULL,
                             width = NULL,
                             height = NULL,
                             silent = FALSE,
                             plotFile = "FileScatters.png",
                             ...){
  #----Warnings----
  if (!is.null(color) & !is.null(groups) & 
      length(unique(groups)) != length(color)){
    stop("Color vector length should be equal to the number of groups.")
  }
  
  if (!is.null(color) & is.null(groups) & length(color) != 1){
    stop("Color vector is too long for only 1 group.")
  }
  
  #---Read in data---
  if (is(input, "flowSet")) {
    data <- flowCore::fsApply(input, function(ff) {
      flowCore::exprs(ff)
    })
    cell_counts <- flowCore::fsApply(input, function(ff) {
      nrow(ff)
    })
    file_values <- unlist(sapply(seq_len(length(cell_counts)), 
                                 function(i) {
                                   rep(i, cell_counts[i])
                                 }))
    ff <- input[[1]]
  } else if (is(input, "flowFrame")) {
    ff <- input
    data <- flowCore::exprs(ff)
    data <- data[,c(channels, fileID)]
    file_values <- data[, fileID]
    input <- unique(file_values)
  } else {
    channels <- FlowSOM::GetChannels(read.FCS(input[1], ...), channels)
    ff <- FlowSOM::AggregateFlowFrames(input,
                                       cTotal = maxPoints, 
                                       channels = channels,
                                       silent = silent)
    data <- ff@exprs
    file_values <- data[, fileID]
  }
  
  subset <- sample(seq_len(nrow(data)), min(maxPoints, nrow(data)))
  if (is.null(channels)) {
    data <- data[subset, , drop = FALSE] 
  } else {
    data <- data[subset, channels, drop = FALSE]
  }
  file_values <- file_values[subset]
  channels <- colnames(data)
  if (fileID %in% channels){
    channels <- channels[-grep(fileID, channels)]
  }
  
  #----Additional warnings---
  if (!is.null(names) & length(unique(file_values)) != length(names)){
    stop("Names vector should have same length as number of files.")
  }
  if (!is.null(groups) & length(unique(file_values)) != length(groups)){
    stop("Groups vector should have same length as number of files.")
  }
  
  #----Organize file names and groups----
  if (is.null(names)) { # if no names are provided, the files will be numbered
    names <- as.character(seq_len(length(input)))
  }
  if (any(duplicated(names))){
    names <- make.unique(names)
  }
  if (is.null(groups)) { # if there are no groups, all files will be labeled "1"
    groups <- rep("1", length(unique(file_values)))
  }
  
  #----Generate plots----
  plots_list <- list()
  yLabel <- pmatch(yLabel, c("channel", "marker"))
  if (length(yLabel) > 2 | any(is.na(yLabel))){
    stop("\"yLabel\" should be \"marker\" and\\or \"channel\"")
  } else {
    yLabel <- c("channel", "marker")[yLabel]
  }
  for (channel in channels) {
    if ("marker" %in% yLabel && length(yLabel) == 1) {
      yLabs <- FlowSOM::GetMarkers(ff, channel)
    } else if ("channel" %in% yLabel && length(yLabel) == 1){
      yLabs <- channel
    } else if (all(c("channel", "marker") %in% yLabel) && length(yLabel) == 2){
      yLabs <- paste0(GetMarkers(ff, channel), " (", channel, ")")
    }
    
    
    df <- data.frame("intensity" = data[, channel],
                     "names" = factor(names[file_values], 
                                      levels = unique(names)),
                     "group" = factor(groups[file_values], 
                                      levels = unique(groups)))
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$names, .data$intensity)) +
      ggplot2::geom_jitter(position = position_jitter(width = 0.1), alpha = 0.5, 
                           ggplot2::aes(colour = .data$group), shape = ".") +
      ggplot2::ylab(yLabs) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         vjust = 0.5,
                                                         hjust = 0.95)) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 5, shape = 15, alpha = 1)))
    if (!is.null(color)) { # if manual colors are provided
      p <- p + ggplot2::scale_color_manual(values = color)
    }
    
    if (!is.null(yLim)) { # if y margins are provided
      p <- p + ggplot2::ylim(yLim)
    }
    
    if (!legend) { # if you don't want a legend on the plot
      p <- p + ggplot2::theme(legend.position = "none")
    }
    
    if(!is.null(quantiles)){
      my_quantile <- function(x, quantiles) {
        dplyr::tibble(intensity = stats::quantile(x, quantiles), 
                      quantile = quantiles)
      }
      
      quantile_intensities <- df %>%
        dplyr::group_by(names) %>% 
        dplyr::summarise(my_quantile(.data$intensity, quantiles))
      p <- p + ggplot2::geom_point(ggplot2::aes(x = .data$names, 
                                                y = .data$intensity), 
                                   col = "black", 
                                   shape = 3, #95,
                                   size = 3,
                                   data = quantile_intensities)
    }
    
    plots_list[[length(plots_list) + 1]] <- p
  }
  
  #----Return plots----
  if (!is.null(plotFile)) {
    if (is.null(nrow) & is.null(ncol)) {
      nrow <- floor(sqrt(length(channels)))
      ncol <- ceiling(length(channels) / nrow)
    } else if (!is.null(nrow) & !is.null(ncol)) {
      if(nrow * ncol < length(channels)) (stop("Too few rows/cols to make plot"))
    } else if (is.null(nrow)) {
      nrow <- ceiling(length(channels) / ncol)
    } else {
      ncol <- ceiling(length(channels) / nrow)
    }
    if (is.null(width)){
      width <-  ncol * (60 + 15 * length(unique(file_values)))
    }
    if (is.null(height)){
      height <-  250 * nrow
    }
    png(plotFile,
        width = width, 
        height = height)
    p <- ggpubr::annotate_figure(ggarrange(plotlist = plots_list,
                                           common.legend = legend, 
                                           ncol = ncol, nrow = nrow),
                                 bottom = ggpubr::text_grob("Files"))
    print(p)  
    dev.off()
  } else {
    return(plots_list)
  }
}