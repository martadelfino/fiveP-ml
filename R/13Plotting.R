#' Function for plotting distribution of fiveP scores
#'
#' @param dis_plot A list of dfs of fiveP results to plot
#' @param bins Number of bins for histogram
#' @importFrom magrittr %>%
#' @return A ggplot2 plot of fiveP distribution scores
#' @export
fivep_distribution_plot <- function(dis_plots, bins = 20) {
  # Convert each to long format and tag with source name
  long_dfs <- lapply(names(dis_plots), function(name) {
    dis_plots[[name]] %>%
      pivot_longer(
        cols = -hgnc_id,
        names_to = "column",
        values_to = "value"
      ) %>%
      dplyr::mutate(source = name)
  }) %>% bind_rows()

  # Plot distributions with density curves
  ggplot(long_dfs, aes(x = value, fill = source)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = bins) +
    facet_wrap(~ column, scales = "free") +
    theme_minimal() +
    coord_cartesian(ylim = c(0, 11000)) +
    labs(
      x = "P threshold",
      y = "Count"
    )
}

#' Function for plotting distribution of fiveP scores with zoomed-in y-axis
#'
#' @param dis_plot A list of dfs of fiveP results to plot
#' @param bins Number of bins for histogram
#' @param ymax Maximum y-axis value for zooming in
#' @importFrom magrittr %>%
#' @return A ggplot2 plot of fiveP distribution scores with zoomed-in y-axis
#' @export
fivep_distribution_plot_zoomedin <- function(dis_plots, bins = 20, ymax = 2000) {
  long_dfs <- lapply(names(dis_plots), function(name) {
    dis_plots[[name]] %>%
      pivot_longer(
        cols = -hgnc_id,
        names_to = "column",
        values_to = "value"
      ) %>%
      dplyr::mutate(source = name)
  }) %>% bind_rows()
  # Precompute non-NA counts per column & dataset
  counts_df <- long_dfs %>%
    group_by(column, source) %>%
    summarise(non_na_count = sum(!is.na(value)), .groups = "drop") %>%
    mutate(
      label_y = 10500 - (row_number() - 1) * 500  # stagger labels
    )

  # Plot
  ggplot(long_dfs, aes(x = value, fill = source)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = bins) +
    facet_wrap(~ column, scales = "free_x") +  # only x is free
    geom_text(
      data = counts_df,
      aes(
        x = Inf,
        y = label_y,
        label = paste0(source, ": n=", non_na_count),
        fill = NULL
      ),
      inherit.aes = FALSE,
      hjust = 1.1,
      size = 4
    ) +
    theme(
      strip.text = element_text(size = 30),   # facet titles bigger
      plot.title = element_text(size = 30)    # main title bigger
    ) +
    coord_cartesian(ylim = c(0, ymax)) +  # fix y-axis for all facets
    theme_minimal() +
    labs(
      x = "P threshold",
      y = "Count"
    )
}

#' Function for plotting GO term similarity of input gene list vs random gene lists
#'
#' @param df A dataframe of GO similarity results to plot
#' @param sd_col The name of the column containing standard deviation values (as a string)
#' @param sd_series The name of the series to which the standard deviation applies (as a string)
#' @param title Optional title for the plot
#' @importFrom magrittr %>%
#' @return A ggplot2 plot of GO term similarity results
#' @export
plot_go_similarity <- function(df, sd_col = NULL, sd_series = NULL, title = NULL) {
  # If sd_col is provided, convert to symbol
  sd_col_sym <- if (!is.null(sd_col)) sym(sd_col) else NULL

  df_long <- df %>%
    dplyr::mutate(threshold = readr::parse_number(threshold_label)) %>%
    pivot_longer(
      cols = -c(threshold, threshold_label, if (!is.null(sd_col)) !!sd_col_sym),
      names_to = "series",
      values_to = "value"
    )

  p <- ggplot(df_long, aes(
    x = threshold,
    y = value,
    color = series,
    group = series
  )) +
    geom_line(size = 1) +
    geom_point()

  # Add error bars only if sd_col is provided and sd_series is a single value
  if (!is.null(sd_col) && !is.null(sd_series) && length(sd_series) == 1) {
    p <- p +
      geom_errorbar(
        data = df_long %>% filter(series == sd_series),
        aes(
          ymin = value - !!sd_col_sym,
          ymax = value + !!sd_col_sym
        ),
        width = 0.02,
        color = "black"
      )
  }

  p +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    coord_cartesian(ylim = c(0.1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = title,
      x = "Threshold",
      y = "Value"
    )
}
