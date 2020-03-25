#' @importFrom tidygraph as_tbl_graph
#' @importFrom dplyr group_by %>% filter ungroup
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @rdname plot.CGLassoFit
#' @export
as_tbl_graph.CGLassoFit <- function(x, ...){
  edge_arrays <- apply(x$edges, 3, function(X) which(lower.tri(X) & X, arr.ind = TRUE))

  edges <- do.call(rbind, lapply(seq_along(edge_arrays), function(i){
    ea <- edge_arrays[[i]]
    cbind(ea,
          lambda = rep(x$lambda[[i]], length.out = NROW(ea)),
          sparsity = rep(NROW(ea), length.out = NROW(ea)))
  }))

  colnames(edges) <- c("from", "to", "lambda", "sparsity")
  edges <- as_tibble(edges, rownames = NULL) %>% group_by(.data$sparsity) %>%
               filter(.data$lambda == min(.data$lambda)) %>%
               ungroup %>%
               ## Add back feature names (from original data matrix X)
               mutate(from = colnames(x$X)[.data$from],
                      to   = colnames(x$X)[.data$to])

  as_tbl_graph(edges, directed = FALSE)
}

#' Plot Complex Graphical Lasso Fits
#'
#' Given a complex graphical lasso path, as estimated by \code{\link{cglasso}},
#' plot several estimated graphs. Under the hood, this is a relatively thin wrapper
#' around the \code{\link[tidygraph]{as_tbl_graph}} and \code{\link[ggraph]{ggraph}}
#' functions, which can be used directly to produce customized visualizations.
#' @export
#' @importFrom ggraph ggraph geom_edge_fan geom_node_point geom_node_text facet_edges
#' @importFrom rlang .data
#' @importFrom dplyr %>% pull filter group_by ungroup select mutate summarize top_n
#' @importFrom ggplot2 theme aes element_rect coord_polar element_line rel
#' @importFrom tidygraph activate
#' @param x An object of class \code{CGLassoFit} as produced by \code{\link{cglasso}}
#' @param ... Additional arguments (currently ignored)
#' @param n_sparsity The number of graphs at different sparsities to display. For
#'                   some values of \code{x}, this might be an upper bound rather
#'                   than exact.
#' @examples
#' Theta <- matrix(c(5, 0, 0, 1, 0,
#'                   0, 5, 0, 0, 0,
#'                   0, 0, 5, 0, 0,
#'                   1, 0, 0, 5, 0,
#'                   0, 0, 0, 0, 5), nrow = 5)
#' Sigma <- solve(Theta) * 5
#' cov_spec <- cplx_cov_spec(Sigma = Sigma)
#'
#' X <- rmvcnorm(500, cov = cov_spec)
#' colnames(X) <- paste0("F", 1:5)
#'
#' cglasso_fit <- cglasso(X)
#' plot(cglasso_fit)
#' plot(cglasso_fit, n_sparsity = 9)
#'
plot.CGLassoFit <- function(x, ..., n_sparsity = 4){
    chkDots(...)

    x_gr <- as_tbl_graph(x) %>% activate(!!as.symbol("edges"))

    x_sparsities <- x_gr %>% pull(.data$sparsity) %>% unique

    x_sparsity_grid_fixed <- seq(min(x_sparsities), max(x_sparsities), length.out = n_sparsity + 1)

    ## Make groups based on "sparsity-buckets" and then get the densest element in each bucket
    x_gr <- x_gr %>% mutate(sparsity_factor = cut(.data$sparsity,
                                                  breaks = x_sparsity_grid_fixed,
                                                  include.lowest = TRUE)) %>%
                     group_by(.data$sparsity_factor) %>%
                     filter(.data$lambda == min(.data$lambda)) %>%
                     ungroup %>%
                     select(-.data$sparsity_factor) %>%
                     mutate(sparsity = naturalfactor(paste("Number of Estimated Edges:", .data$sparsity)))

    ggraph(x_gr, layout = "linear", circular = TRUE) +
      geom_edge_fan() + geom_node_point() +
      geom_node_text(aes(label = .data$name), repel = TRUE) +
      ggtitle("Complex Graphical Lasso Fit") +
      facet_edges(~sparsity) +
      # Adapted from the definition of theme_bw() in ggplot2
      # Normal + theme_bw() doesn't work as it adds back in x and y ticks
      theme(
        # white background and dark border
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border     = element_rect(fill = NA, colour = "grey20"),
        # make gridlines dark, same contrast with white as in theme_grey
        panel.grid = element_line(colour = "grey92"),
        panel.grid.minor = element_line(size = rel(0.5)),
        # contour strips to match panel contour
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        # match legend key to background
        legend.key       = element_rect(fill = "white", colour = NA),
      )
}
