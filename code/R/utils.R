
#' Given a data frame, provide the counts of groupings acording to
#' two columns. Also provide the count percentage relative to all entries.
#'
#' @param df Input data frame that contains to the two columns of interest
#' @param group_col1 string of first column name 
#' @param group_col2 string of second column name
#'
#' @return
#' @export
#'
#' @examples
two_group_count_percentage <- function(df, group_col1, group_col2) {
  n <- dim(df)[[1]]
  df %>%
    group_by(!!sym(group_col1), !!sym(group_col2)) %>%
    summarise(count = n()) %>%
    mutate(percentage = (count / n) * 100, 
           counts_and_percent = paste0(count, " (", round(percentage, 1), "%)")) %>%
    ungroup()
}
