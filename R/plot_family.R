#' Create a Bubble Plot showing number of sequences and average/range length for each primer
#'
#' @param refdb (refdb) a refdb object.
#' @param level (character) a character vector indicating which taxonomic level
#' the x axis should refer to. It should be one of "phylum", "class", "order",
#' "family" and "genus".
#' @param select (character) this parameter allows to pre-filter the refdb
#' object based on a taxonomic name, which must be present in the records table
#' of the object.
#' @param size_range (integer vector) this allows to change the size of bubbles.
#' @param measure (character) either "range" or "average" (default). In the
#' first case, the range of length in the sequences of the records filtered for
#' each combination of primer and taxonomic name will be shown. The average
#' length will be shown otherwise.
#' @param tax.fct.levels (character) a vector of taxonomic names that allows to
#' re-order the x axis based on a custom order.
#' @param prim.fct.levels (character) a vector of primers' names that allows to
#' re-order the y axis based on a custom order.
#'
#' @return A Bubble plot, as a 'ggplot2' plot object, which can be saved and
#' modified.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples \dontrun{
#' (BubPlot <- plot_primers(total_ore))
#' }
plot_primers <- function(refdb, level = "phylum", select = NULL, size_range = NULL, measure = "average", tax.fct.levels = NULL, prim.fct.levels = NULL) {

  if (!(measure %in% c("average", "range")) | (length(measure) > 1)) {

    stop("The measure parameter must be one of 'average' or 'range'")

  }

  if (!(level %in% c("phylum", "class", "order", "family", "genus"))) {

    stop("The taxonomic level must be one of 'phylum', 'class', 'order', 'family' or 'genus'")

  }

  if (!is.null(select)) {

    selection_choices <- unique(c(refdb$phylum,
                                  refdb$class,
                                  refdb$order,
                                  refdb$family,
                                  refdb$genus,
                                  refdb$species))

    if (!(select %in% selection_choices) | (length(select) > 1)) {

      stop("'Select' parameter must be one of any taxonomic names included in the record table")

    } else {

      refdb <- refdb %>% dplyr::filter(.data$phylum == {{select}} |
                                       .data$class == {{select}} |
                                       .data$order == {{select}} |
                                       .data$family == {{select}} |
                                       .data$genus == {{select}} |
                                       .data$species == {{select}}
      )
    }

  }

  # filter required fields
  primers_tab_sel <- refdb %>% dplyr::filter(!is.na(.data$directionPrimers)) %>%
    dplyr::select({{level}},
                  directionPrimers,
                  PCR_primers,
                  lengthGene)

  # modify primer information. Only specific information on forward and reverse
  # primers can be easily computed, thus all combinations with multiple fwd and
  # rev primers will be excluded (e.g. "F|F|R|R").
  # Also, reverse the primers names if they are reported as R|F in "directionPrimers"

  primers_tab_filt <- primers_tab_sel %>%
    dplyr::filter(.data$directionPrimers == "F|R" |
                  .data$directionPrimers == "R|F" |
                  .data$directionPrimers == "f|r"|
                  .data$directionPrimers == "r|f",
                  .data$PCR_primers != "|",
                  !is.na(.data[[level]])) %>%
    dplyr::mutate(directionPrimers = stringr::str_replace(directionPrimers, "f", "F")) %>%
    dplyr::mutate(directionPrimers = stringr::str_replace(directionPrimers, "r", "R")) %>%
    dplyr::mutate(PCR_primers = ifelse(directionPrimers == "R|F",
                                       stringr::str_replace(PCR_primers, "(.*)\\|(.*)", "\\2|\\1"),
                                       PCR_primers)) %>%
    dplyr::mutate(directionPrimers = ifelse(directionPrimers == "R|F",
                                            stringr::str_replace(directionPrimers, "(.*)\\|(.*)", "\\2|\\1"),
                                            directionPrimers))

  # finally, prepare the table that will be processed by ggplot for bubble-plot
  # creation. This will include, for each level and primer combination, range of
  # length (difference between minimum and maximum length of sequences),
  # average length of sequences and number of records
  primers_tab_plot <- primers_tab_filt %>%
    dplyr::group_by(.data[[level]], .data$PCR_primers) %>%
    dplyr::summarise(range_length = (max(.data$lengthGene) - min(.data$lengthGene)),
                     avg_length = mean(.data$lengthGene),
                     n = dplyr::n())

  # taxonomic levels' order
  if (!is.null(tax.fct.levels)) {

    if (!all(tax.fct.levels %in% unique(primers_tab_plot[, level, drop = TRUE]))) {

      stop("Taxonomic levels indicated in 'tax.fct.levels' do not correspond to selected taxonomic names")

    } else {

      primers_tab_plot <- primers_tab_plot %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dplyr::across(.data[[level]],
                                    as.factor)) %>%
        dplyr::mutate({{level}} := forcats::fct_relevel(.[, level, drop = TRUE],
                                                        tax.fct.levels))

    }

  }

  # primers levels' order
  if (!is.null(prim.fct.levels)) {

    if (!all(prim.fct.levels %in% unique(primers_tab_plot[, "PCR_primers", drop = TRUE]))) {

      stop("Taxonomic levels indicated in 'prim.fct.levels' do not correspond to selected primers' names")

    } else {

      primers_tab_plot <- primers_tab_plot %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dplyr::across(.data$PCR_primers,
                                    as.factor)) %>%
        dplyr::mutate(PCR_primers = forcats::fct_relevel(PCR_primers,
                                                         prim.fct.levels))

    }

  }

  # set measure to plot
  if (measure == "average") {

    setMeasure <- "avg_length"
    setLegendTitle <- "Average\nsequence length (bp)"

  } else {

    setMeasure <- "range_length"
    setLegendTitle <- "Range of\nsequence length (bp)"

  }

  # Finally, print the bubble- plot

  bubblePlot <- ggplot2::ggplot(primers_tab_plot) +
    ggplot2::aes(.data[[level]], PCR_primers, color = .data[[setMeasure]], size = n) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(color = setLegendTitle,
                  size = "Number\nof sequences") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1,
                                                       vjust = 1,
                                                       size = 9))

  if (is.null(size_range)) {

    bubblePlot <- bubblePlot + ggplot2::scale_size_continuous(range = c(2, 30))

  } else {

    bubblePlot <- bubblePlot + ggplot2::scale_size_continuous(range = size_range)

  }

  bubblePlot

}

# set global variables from usage of dlpyr
utils::globalVariables(c("n", "directionPrimers", "PCR_primers"))

#' Create a barplot using the function geom_density_ridges from ggplot2, to compare
#' the final length of the sequences (field lengthGene) with the original length
#' of sequences (field lengthSource).
#'
#' @param refdb (refdb) a refdb object.
#' @param limit (numeric) the length in "base pairs" with which the sequences
#' longer than it will be grouped and reported as the last bar on the x axis.
#' @param breaks (numeric) the length range in "base pairs" of the x axis.
#' @param level (character) a character vector indicating which taxonomic level
#' the facet_wrap function from ggplot2 should refer to. It should be one of
#' "phylum", "class", "order", "family" and "genus".
#' @param select (character) this parameter allows to pre-filter the refdb
#' object based on a taxonomic name, which must be present in the records table
#' of the object.
#' @param scaling (logical) defaults to TRUE, it reports the height of the bars
#' according to the sequence count in the refdb object, thus, if used with the
#' level parameter, the heigth of bars in different wraps can be compared.
#'
#' @return a geom_density_ridges ggplot2 plot.
#'
#' @export
#'
#' @examples \dontrun{
#' (lengthPlot <- plot_length(refdb))
#' }
plot_length <- function(refdb, limit = 1000, breaks = 50, level = NULL, select = NULL, scaling = TRUE) {

  # set visible binding to variable
  Length <- NULL

  if (!is.null(select)) {

    selection_choices <- unique(c(refdb$phylum,
                                  refdb$class,
                                  refdb$order,
                                  refdb$family,
                                  refdb$genus,
                                  refdb$species))

    if (!(select %in% selection_choices) | (length(select) > 1)) {

      stop("'Select' parameter must be one of any taxonomic names included in the record table")

    } else {

      refdb <- refdb %>% dplyr::filter(.data$phylum == {{select}} |
                                       .data$class == {{select}} |
                                       .data$order == {{select}} |
                                       .data$family == {{select}} |
                                       .data$genus == {{select}} |
                                       .data$species == {{select}}
      )
    }

  }

  # set breaks bins' numbers (e.g. c("0", "0+breaks", "breaks+breaks" etc.)
  limits <- seq(from = 0, to = limit, by = breaks)

  # get all possible levels for future bin factor
  bin_levels <- sapply(seq(1:length(limits)), function(set) {

    if ((set) == max(seq(1:length(limits)))) {

      lev <- c("", paste0("over ", limit))
      return(lev)

    } else {

      lev <- paste0(limits[set], "-",
                    limits[set+1])
      return(lev)

    }})

  bin_levels <- do.call(c, bin_levels)

  # create table for plotting length distributions
  length_tab <- refdb %>% dplyr::select({{level}}, lengthSource, lengthGene) %>%
    tidyr::pivot_longer(., cols = c(lengthSource, lengthGene),
                        names_to = "LengthType",
                        values_to = "Length")

  # create bin categories
  bins <- lapply(seq(1:length(limits)), function(set) {

    if ((set) == max(seq(1:length(limits)))) {

      binned <- length_tab %>%
        dplyr::filter(.data$Length > limits[set] &
                        .data$Length <= max(c(refdb$lengthSource, refdb$lengthGene))) %>%
        dplyr::mutate(bin = paste0("over ", limit),
                      bin_num = length(bin_levels))

      return(binned)

    } else {

      binned <- length_tab %>%
        dplyr::filter(.data$Length > limits[set] &
                        .data$Length < limits[set+1]) %>%
        dplyr::mutate(bin = paste0(limits[set], "-", limits[set+1]),
                      bin_num = set)

      return(binned)

    }
  })

  # create final table for plotting
  length_tab_plot <- do.call(rbind, bins) %>%
    dplyr::select(-Length) %>%
    dplyr::mutate(dplyr::across(c({{level}}, LengthType), as.factor)) %>%
    dplyr::mutate(bin = forcats::fct(bin,
                                     bin_levels))

  if (scaling) {

    # set parameters to get heights of bins corresponding to actual count values
    height <- "count"
    scale <- FALSE

  } else {

    # set parameters to get heights of bins as default of ggridges
    height <- "density"
    scale <- TRUE

  }

  # plot histogram
  histPlot <- ggplot2::ggplot(length_tab_plot,
                              ggplot2::aes(x = bin_num,
                                           y = LengthType,
                                           fill = LengthType,
                                           height = ggplot2::after_stat(.data[[height]]))) +
    ggridges::geom_density_ridges(stat = "binline",
                                  panel_scaling = scale,
                                  binwidth = 1,
                                  alpha = .6) +
    ggplot2::scale_x_continuous(breaks = as.numeric(c("0", seq(1:(length(bin_levels)+2)))),
                                labels = c("", bin_levels, "", "")) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(.3, 1.9))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1,
                                                       vjust = 1),
                   panel.grid.minor = ggplot2::element_line(colour = "white")) +
    ggplot2::labs(fill = "Length type") +
    ggplot2::xlab("Length range (bp)") +
    ggplot2::ylab("") +
    ggplot2::scale_fill_brewer(palette="Set2",
                               labels = c("Length of chosen marker",
                                          "Length of whole sequence"),
                               guide = ggplot2::guide_legend(reverse=TRUE))

  if (!is.null(level)) {

    if (!all(level %in% c("phylum", "class", "order", "family", "genus", "species", "origin"))) {

      stop("The level parameter must be one of 'phylum', 'class', 'order', 'family', 'genus', 'species' or 'origin'")

    }

    histPlot <- histPlot + ggplot2::facet_wrap(level)

  }

  histPlot

}

# set global variables from usage of dlpyr
utils::globalVariables(c("bin", "bin_num", "lengthGene", "lengthSource", "LengthType"))
