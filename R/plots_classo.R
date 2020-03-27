#' Visualizing Complex Lasso Regularization Paths
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot lines text legend
#' @importFrom stats runif
#' @importFrom ggplot2 ggplot theme_bw scale_color_brewer geom_line
#' @importFrom ggplot2 xlab ylab theme_bw ggtitle aes coord_polar
#' @importFrom patchwork plot_layout
#' @importFrom dplyr filter mutate group_by ungroup %>% rename summarize top_n
#' @export
#' @param x An \code{CLassoFit} object produced by \code{\link{classo}}
#' @param xvar Value to use on the x-axis (ordinate). \itemize{
#'    \item \code{lambda}: The log of the regularization parameter \eqn{\lambda}{lambda}.
#'    \item \code{l1}: The \eqn{\ell_1}{l1}-norm of the solution
#'    }
#' @param yvar Value to use on the y-axis (abscissa).\itemize{
#'    \item \code{Mod}: the modulus (magnitude) of the estimated coefficients
#'    \item \code{Re}: the real part of the estimated regression coefficients
#'    \item \code{Im}: the imaginary part of the estimated regression coefficients
#'    \item \code{Arg}: the argument (phase) of the estimated coefficients
#' } Ignored if \code{gg = TRUE}.
#' @param label Should variables be labeled? Set \code{label=NULL} to label only
#'    those variables which are non-zero at the sparsest end of
#'    the regularization path. Only used if \code{gg = FALSE}
#' @param legend Should a legend of top variables (heuristically selcted) be displayed?
#'    Set \code{legend=NULL} to disable, else give the location of the legend.
#'    Only used if \code{gg = FALSE}
#' @param gg Should a \code{\link[ggplot2]{ggplot2}}-based visualization be returned?
#'           The default (\code{gg = FALSE}) plots are more standard and directly
#'           comparable to those produced by the \code{\link[glmnet]{glmnet}} package,
#'           while the \code{gg = TRUE} plots have are able to capture more of the
#'           complex structure and have more adaptive defaults.
#' @param ... Additional arguments passed to plotting functions
#' @param n_highlight Number of "top" variables to display when \code{gg = TRUE}.
#' @examples
#' n <- 200
#' p <- 500
#' groups <- rep(1:10, times=50)
#' beta <- numeric(p);
#' beta[1:10] <- 3 * exp(im * runif(10, 0, 2*pi))
#'
#' X <- matrix(rcnorm(n * p), ncol=p)
#' y <- X %*% beta + rcnorm(n)
#'
#' cfit <- classo(X, y)
#' plot(cfit)
#' plot(cfit, legend=NULL, xvar="lambda")
#' plot(cfit, gg = TRUE)
plot.CLassoFit <- function(x,
                           ...,
                           xvar = c("norm", "lambda"),
                           yvar = c("Mod", "Re", "Im", "Arg"),
                           label=FALSE,
                           legend = "topright",
                           gg = FALSE,
                           n_highlight = 10){

    chkDots(...)

    xvar <- match.arg(xvar)
    yvar <- match.arg(yvar, several.ok = gg)

    if(!gg){
        plot_data <- if(has_intercept(x)) coef(x)[-1,, drop=FALSE] else coef(x)

        xlab <- if(xvar == "norm") "L1 Norm" else "Log Lambda"
        ## Base Plotting
        xvar <- if(xvar == "norm") colSums(Mod(plot_data)) else x$lambda
        yvar <- match.fun(yvar)(plot_data)

        COLORS <- brewer.pal(9, "Set1")

        plot(xvar, yvar[1,], type="n", xlab=xlab,
             ylim=range(yvar),
             xlim=range(xvar),
             ylab=expression(hat(beta)))

        if(length(colnames(x$X))){
            labels <- colnames(x$X)
            do_legend <- !is.null(legend)

            if(is.null(label)){
                do_label <- TRUE
            } else{
                do_label <- label
            }
        } else {
            do_legend <- FALSE
            do_label <- FALSE
        }

        legend_labels <- character()
        legend_colors <- character()
        color_counter <- 0

        for(i in 1:NCOL(x$X)){
            y <- yvar[i,]

            if(max(abs(y)) == 0){
                next
            }

            color_counter <- color_counter + 1
            col <- COLORS[color_counter  %% 9 + 1]

            lines(xvar, y, col=col, ...)

            if(do_label || do_legend){
                if(do_legend){
                    legend_labels <- c(legend_labels, labels[i]);
                    legend_colors <- c(legend_colors, col);
                }

                if(do_label){
                    text_x <- range(xvar) * c(0.95, 1.05)
                    text_x <- min(text_x)

                    text(text_x, y[length(y)] * runif(1, 0.95, 1.05),
                         labels[i], col=col, offset=3)
                }
            }
        }

        if(do_legend){
            legend(legend, legend=legend_labels, col=legend_colors, lty=1, lwd=2, bg="white", ...)
        }
    } else {
        ## ggplot2 Plotting

        ## Get data ready to plot
        plot_data <- as.data.frame(x) %>% filter(.data$variable != "(Intercept)") %>%
            mutate(Mod = Mod(.data$coef),
                   Re  = Re(.data$coef),
                   Im  = Im(.data$coef),
                   Arg = Arg(.data$coef),
                   log_lambda = log(.data$lambda)) %>%
            ## Move Arg() to [0, 2 * pi)
            mutate(Arg = ifelse(.data$Arg < 0, .data$Arg + 2 * pi, .data$Arg)) %>%
            group_by(.data$variable) %>%
            filter(max(.data$Mod) != 0) %>%
            ungroup %>%
            group_by(.data$lambda) %>%
            mutate(norm = sum(.data$Mod)) %>%
            ungroup

        top_vars <- plot_data %>% group_by(.data$variable) %>%
                                  summarize(max_mod = max(.data$Mod)) %>%
                                  top_n(n_highlight, wt = .data$max_mod) %>%
                                  pull(.data$variable)

        plot_data <- plot_data %>% filter(.data$variable %in% top_vars) %>%
                                   mutate(variable = naturalfactor(.data$variable))

        ## Create four "basic" plots - we will build these all "manuall"
        ## and then use patchwork to combine them in a nice layout

        xlabel <- if(xvar == "norm") expression("\u2016"*beta*"\u2016"[1]) else expression(log(lambda))

        plot_base <- ggplot(plot_data, aes(x = !!as.symbol(xvar),
                                           group = .data$variable,
                                           color = .data$variable)) +
                                       theme_bw() +
                                       scale_color_brewer(name = paste("Top", n_highlight, "Variables"),
                                                          palette = "Paired")

        plot_mod <- plot_base + geom_line(aes(y = .data$Mod)) +
                                xlab(xlabel) + ylab(expression(abs(hat(beta[j])))) +
                                ggtitle("Modulus")
        plot_re  <- plot_base + geom_line(aes(y = .data$Re))  +
                                xlab(xlabel) + ylab(expression(Re(hat(beta[j])))) +
                                ggtitle("Real Part")
        plot_im  <- plot_base + geom_line(aes(y = .data$Im))  +
                                xlab(xlabel) + ylab(expression(Im(hat(beta[j])))) +
                                ggtitle("Imaginary Part")

        # plot_arg is a bit distinctive
        ## FIXME - This still looks a bit off - we can get weird "wraps" that shouldn't be there
        ##         compare to geom_point() on the same.
        ## geom_path gives better connections, but can still give weird wraps when
        ## Arg(beta_hat) fidgets around 0....
        plot_arg <- plot_base + geom_path(aes(x = .data$Arg, y = !!as.symbol(xvar))) +
                      coord_polar(start = - pi / 2, direction = -1) + ggtitle("Argument (Phase)") +
                      ylab(xlabel) + xlab(expression(Arg(hat(beta[j])))) +
                      scale_x_continuous(breaks = c(0, pi / 2, pi, 3 * pi / 2),
                                         labels = c(expression(0),
                                                    expression(pi / 2),
                                                    expression(pi),
                                                    expression(3 * pi / 2)),
                                         limits = c(0, 2 * pi))

        (plot_re + plot_im) / (plot_mod + plot_arg) + plot_layout(guides = "collect")
    }
}

#' @param x An \code{CLassoFit_cv} object as produced by
#'    \code{\link{cv.classo}}.
#' @param bar.width Width of error bars
#' @export
#' @rdname cv.classo
#' @importFrom graphics segments points abline
plot.CLassoFit_cv <- function(x, bar.width=0.01, ...){

    log_lambda <- log(x$lambda)
    bar.width <- bar.width * diff(range(log_lambda))

    plot(log_lambda, x$cvm,
         ylim=range(c(x$cvup, x$cvlo)),
         xlab=expression(log(lambda)),
         ylab=toupper(x$name),
         type="n", ...)

    points(log_lambda, x$cvm, pch=20, col="red")

    segments(log_lambda, x$cvlo,
             log_lambda, x$cvup,
             col="darkgrey")
    segments(log_lambda - bar.width, x$cvlo,
             log_lambda + bar.width, x$cvlo,
             col="darkgrey")
    segments(log_lambda - bar.width, x$cvup,
             log_lambda + bar.width, x$cvup,
             col="darkgrey")

    abline(v=log(x$lambda.min), lty = 3, lwd = 2)
    abline(v=log(x$lambda.1se), lty = 3, lwd = 2)
}
