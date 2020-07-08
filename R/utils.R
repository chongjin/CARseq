#' Permute case and control labels symmetrically
#' 
#' @param labels a character vector of labels with exactly two levels.
#' @export
permute_case_and_controls = function(labels) {
  # Use a more balanced way of resampling:
  # Same proportion of case and control samples are re-labelled as "case" sample after permutation.
  if (!is.factor(labels)) labels = as.factor(as.character(labels))
  label_levels = levels(labels)
  if (length(label_levels) <= 1) stop("There are not enough levels. Two levels are required in the vector of labels for the function to work.")
  if (length(label_levels) >= 3) stop("There are at least three levels. Two levels are required in the vector of labels for the function to work.")
  is_group1 = (labels == label_levels[1])
  proportion_of_group1 = sum(is_group1) / length(labels)
  labels[is_group1] = sample(c(rep(label_levels[1], round(sum(is_group1) * proportion_of_group1)),
                            rep(label_levels[2], sum(is_group1) - round(sum(is_group1) * proportion_of_group1))))
  labels[!is_group1] = sample(c(rep(label_levels[1], sum(is_group1) - round(sum(is_group1) * proportion_of_group1)),
                             rep(label_levels[2], sum(!is_group1) - sum(is_group1) + round(sum(is_group1) * proportion_of_group1))))
  labels
}

#' Calculate q-values (Storey 2003) from a p-value distribution with considerable amount of 1's
#' 
#' @param pvals a numeric vector of p-values.
#' 
#' @examples 
#' set.seed(1234)
#' pval = c(runif(100, 0, 0.1), runif(100,0, 1), rep(1, 100))
#' qval = get_qvalues_one_inflated(pval)
#' plot(sort(pval), sort(qval), xlim=c(0,1), ylim=c(0,1), type="l",
#'      xlab="p-values", ylab="q-values")
#' qval_range = seq(0.001, 0.25, by=0.001)
#' plot(qval_range,
#'      sapply(qval_range, function(x) sum(qval <= x, na.rm=TRUE)),
#'      type="l",
#'      xlab="q-value", ylab="number of significant genes")
#' @export
get_qvalues_one_inflated = function(pvals) {
  pvals_order = order(pvals)
  pvals_sorted = pvals[pvals_order]
  m = sum(!is.na(pvals_sorted))
  # pi0 = (proportion of p-value = 1) + 2*(proportion of p-value > 0.5 and < 1)
  pi0 = min(1, (sum(pvals_sorted == 1, na.rm=TRUE) + 2 * sum(pvals_sorted > 0.5 & pvals_sorted < 1, na.rm=TRUE)) / m)
  qvals = rep(NA, m)
  qvals[m] = pvals_sorted[m] * pi0
  for (i in rev(seq_len(m-1))) {
    qvals[i] = min(pi0*m*pvals_sorted[i]/i, qvals[i+1]) 
  }
  qvals[rank(pvals)]
}

#' Generate 3d barplot.
#' 
#' This function is used to illustrate the model-based mean expression in the CARseq paper.
#' 
#' @param x a numeric vector of proportions that will be shown on the x axis
#' @param y a numeric vector of heights that will be shown on the y axis
#' @param z a numeric vector of depths
#' @param col a numeric vector of colors
#' @param zvec a numeric vector of length 2 for the direction of depth. Both entries should be positive.
#' @param alpha a transparency value taking value between 0 and 1, where 1 is opaque and 0 is transparent
#' @param border the border color
#' @param xlim user-specified xlim
#' @param ylim user-specified ylim
#' 
#' @examples 
#' barplot_3d(x = 1:3, y = 1:3, z = 1)
#' 
#' @export
barplot_3d = function(x, y, z, col = NULL, zvec = c(0.1, 0.1), alpha = 0.8, border = "grey", xlim = NULL, ylim = NULL) {
  n = length(y)  # number of blocks
  if (length(x) == 1) {
    x = rep(1/n, n)
  } else if (length(x) != n) {
    stop("The length of x is not equal to the length of y!")
  }
  x = x / sum(x)
  if (length(z) == 1) {
    z = rep(z, n)
  } else if (length(z) != n) {
    stop("The length of z is not equal to the length of y!")
  }
  if (is.null(col)) {
    # by default, use rainbow palette
    col = hcl.colors(n, "Pastel 1")
  } else if (length(col) == 1) {
    col = rep(col, n)
  } else if (length(col) != n) {
    stop("The length of col is not equal to the length of y!")
  }
  col = rgb(t(col2rgb(col)/255), alpha = alpha)  # add alpha
  
  # empty plot
  if (is.null(xlim)) {
    xlim = c(0, 1 + max(z) * zvec[1])
  }
  if (is.null(ylim)) {
    ylim = c(0, max(y) + max(z) * zvec[2])
  }
  plot(0, type='n', axes=FALSE, ann=FALSE, xlim = xlim, ylim = ylim)
  x_start = c(0, cumsum(x)[-n])
  x_end = cumsum(x)
  for (i in seq_len(n)) {
    front_col = col[i]
    segments(x0 = x_start[i], y0 = 0,
             x1 = x_start[i] + z[i] * zvec[1], y1 = z[i] * zvec[2],
             col = border)
    segments(x0 = x_end[i] + z[i] * zvec[1], y0 = z[i] * zvec[2],
             x1 = x_start[i] + z[i] * zvec[1], y1 = z[i] * zvec[2],
             col = border)
    segments(x0 = x_start[i] + z[i] * zvec[1], y0 = y[i] + z[i] * zvec[2],
             x1 = x_start[i] + z[i] * zvec[1], y1 = z[i] * zvec[2],
             col = border)
    # front
    polygon(x = c(x_start[i], x_end[i], x_end[i], x_start[i]),
            y = c(0, 0, y[i], y[i]),
            border = border,
            col = front_col)
    # top
    polygon(x = c(x_start[i], x_end[i], x_end[i] + z[i] * zvec[1], x_start[i] + z[i] * zvec[1]),
            y = c(y[i], y[i], y[i] + z[i] * zvec[2], y[i] + z[i] * zvec[2]),
            border = border,
            col = colorRampPalette(c(front_col, "black"))(10)[2])  # slightly darker than col[i]
    # right side
    polygon(x = c(x_end[i], x_end[i] + z[i] * zvec[1], x_end[i] + z[i] * zvec[1], x_end[i]),
            y = c(0, z[i] * zvec[2], y[i] + z[i] * zvec[2], y[i]),
            border = border,
            col = colorRampPalette(c(front_col, "black"))(10)[4])  # much darker than col[i]
  }
}