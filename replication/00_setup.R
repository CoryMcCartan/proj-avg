suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scales)
    library(patchwork)
    library(wacolors)

    library(redist)
    library(alarmdata)
    library(ggredist)
    rlang::check_installed("rmapshaper", "for plotting")

    library(matrixStats)
    library(here)
})

theme_paper = function() {
    theme_void(base_size=10, base_family="Helvetica") +
        theme(legend.position="bottom",
              legend.key.width=unit(1.0, "cm"))
}

scale_contr = function(name, ...) {
    scale_fill_wa_c("vantage", name=name, midpoint=0, reverse=TRUE, ...)
}
scale_contr_b = function(name, ...) {
    scale_fill_wa_b("vantage", name=name, reverse=TRUE, ...,
                    rescaler=\(x, ...) rescale_mid(x, ..., mid=0))
}
scale_avg = function(name, ...) {
    scale_fill_wa_c("puget", name=name, ...)
}

scale_fill_party_c = function(
        name = "Vote share", midpoint = 0.5, limits = 0:1,
        labels = label_party_pct(), oob = scales::squish, reverse = FALSE, ...)  {
    pal = ggredist$partisan
    if (reverse)
        pal = rev(pal)
    if (is.null(midpoint)) {
        rescaler = scales::rescale
    } else {
        rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
            scales::rescale_mid(x, to, from, midpoint)
        }
    }
    ggplot2::scale_fill_gradientn(
        name = name, limits = limits, labels = labels,
        oob = oob, colours = pal, rescaler = rescaler, ...
    )
}


proj_distr <- function(plans, x, draws=NA) {
    plans_m <- get_plans_matrix(plans)
    n_ref <- 0

    if (!is.null(colnames(plans_m))) {
        refs <- which(nchar(colnames(plans_m)) > 0)
        n_ref <- length(unique(colnames(plans_m)[refs]))
    }
    if (is.null(draws)) {
        draw_idx <- seq_len(ncol(plans_m))
    } else if (length(draws) == 1 && is.na(draws)) {
        if (n_ref > 0) {
            draw_idx <- seq_len(ncol(plans_m))[-seq_len(n_ref)]
        } else {
            draw_idx <- seq_len(ncol(plans_m))
        }
    } else if (is.logical(draws)) {
        draw_idx <- which(draws)
    } else {
        draw_idx <- match(as.character(draws), levels(plans$draw))
    }

    plans <- arrange(plans, as.integer(.data$draw), .data$district)
    n_distr <- max(plans_m[, draw_idx[1]])
    m_val <- matrix(rlang::eval_tidy(rlang::enquo(x), plans),
                    nrow = n_distr)
    plans_m <- plans_m[, draw_idx, drop = FALSE]
    m_val <- m_val[, draw_idx, drop = FALSE]
    m_prec <- matrix(nrow = nrow(plans_m), ncol = ncol(plans_m))
    for (i in seq_len(ncol(plans_m))) {
        m_prec[, i] <- m_val[, i][plans_m[, i]]
    }
    m_prec
}

proj_contrast = function(plans, x, comp=1) {
    d_enac = rlang::eval_tidy(rlang::enquo(x), plans)[1:attr(plans, "ndists")]
    d_enac[as.matrix(plans)[, comp]] - rowMeans2(proj_distr(plans, {{ x }}))
}
proj_avg = function(plans, x) {
    rowMeans2(proj_distr(plans, {{ x }}))
}
proj_contr_norm = function(plans, x, comp=1) {
    d_enac = rlang::eval_tidy(rlang::enquo(x), plans)[1:attr(plans, "ndists")]
    pd = proj_distr(plans, {{ x }})
    (d_enac[as.matrix(plans)[, comp]] - rowMeans2(pd)) / rowSds(pd)
}


