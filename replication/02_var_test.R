V = nrow(map)
E = sum(lengths(map$adj))
N = 5000
n = attr(plans, "ndists")

# Sampled projective contrasts fig --------
set.seed(5118)
idx_shuf = sample(2:8)
idx = c(idx_shuf[1:3], 1, idx_shuf[4:7])
p = purrr::map(idx, function(i) {
    plot(map, proj_contrast(plans, ndshare, comp=i)) +
        scale_fill_party_c("Democratic\nvote share\nversus simulations\n",
                           labels = label_party_margin(midpoint=0.0),
                           midpoint = 0, limits=c(-0.15, 0.15)) +
        theme_paper()
}) |>
    wrap_plots(nrow=2) +
    plot_layout(guides="collect") +
    plot_annotation(theme = theme(plot.margin = margin())) &
    theme(legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"),
          legend.key.width=unit(0.6, "cm"))
ggsave(here("paper/figures/contr_samp.pdf"), plot=p, width=9.5, height=3.5)

# Multiple testing ------------
ee = rlang::expr(ndshare)
m = proj_distr(plans, !!ee) # V rows, 5000 cols
p_avg = rowMeans2(m)
m_pc = m - p_avg
pc_enac = rlang::eval_tidy(ee, plans)[as.matrix(plans)[, 1]] - p_avg

m_pc_pv = rowRanks(abs(m_pc), ties.method="max") / ncol(m)
pc_pv_enac = rowMeans2(abs(m_pc) >= abs(pc_enac))
hist(m_pc_pv)
hist(pc_pv_enac)

# check r-hats
rhats = map_dbl(1:V, function(i) {
    redist:::diag_rhat(pc_enac[i] - m[i, ], rep(1:2, each=N/2), split=FALSE)
})
hist(rhats)
if (any(rhats > 1.05)) stop("Projective contrast distribution did not converge.")


# Algorithm 1 of Storey & Tibshirani (2001)
est_pfdr <- function(p, ref=m_pc_pv, lambda=0.5, n_alpha=15) {
    alphas = quantile(p, probs=seq(0, 1, length.out=n_alpha), type=1, names=FALSE)
    m = length(p)
    pi0 = (m - sum(p <= lambda)) / (m - mean(colSums2(ref <= lambda)))
    pi0 = min(pi0, 1)
    fdrs = map_dbl(alphas, \(a) {
        pi0 *  mean(colSums2(ref <= a)) / sum(p <= a)
    })
    ok = !is.nan(fdrs)
    tibble::new_tibble(list(
        alpha = alphas,
        pfdr = coalesce(fdrs, 0.0)
    ))
}
# find highest (most powerful) threshold that controls FDR using linear interp.
find_thresh <- function(ests, fdr=0.05) {
    idx = nrow(ests) + 1 - which(rev(ests$pfdr) <= fdr)[1]
    if (is.na(idx)) return(0.0) # force no rejections
    if (idx == nrow(ests)) return(tail(ests$alpha, 1))

    adj = diff(ests$alpha[idx + 0:1]) / diff(ests$pfdr[idx + 0:1])
    ests$alpha[idx] + (fdr - ests$pfdr[idx]) * adj
}
# Produce q values with spline
qvalues <- function(ests, p) {
    splinefun(ests$alpha, rev(cummin(rev(ests$pfdr))), method="monoH.FC")(p)
}

# Make a (possibly normalized) projective contrast plot with FDR-adjusted
# discoveries marked
plot_contr_fdr <- function(map, plans, x, fdr=0.05, draw=1, norm=FALSE,
                           density=0.2, spacing=0.015) {
    x = rlang::enquo(x)
    m = proj_distr(plans, !!x) # V rows, 5000 cols
    p_avg = rowMeans2(m)
    m = m - p_avg
    n = attr(plans, "ndists")
    pc_draw = rlang::eval_tidy(x, plans)[(draw - 1)*n + 1:n][as.matrix(plans)[, draw]] - p_avg

    m_pv = rowRanks(abs(m), ties.method="max") / ncol(m)
    pv_draw = rowMeans2(abs(m) >= abs(pc_draw))

    q = est_pfdr(pv_draw, m_pv) |>
        qvalues(pv_draw)

    if (isTRUE(norm)) {
        pc_draw = pc_draw / sqrt(rowMeans2(m^2))
    }
    p = plot(map, pc_draw) + theme_paper()
    if (!any(q <= fdr)) return(p)

    geom_sel = map |>
        filter(q <= fdr) |>
        summarize(is_coverage=TRUE) |>
        suppressWarnings()

    p +
        geom_sf(data=geom_sel, inherit.aes=FALSE,
                fill=NA, linewidth=0.15, color="black") +
        ggpattern::geom_sf_pattern(data=geom_sel, inherit.aes=FALSE,
                                   fill=NA, linewidth=0.0, color=NA,
                                   pattern_fill="#00000070", pattern_color=NA,
                                   pattern_density=density, pattern_spacing=spacing)
}


p1 = plot_contr_fdr(map, plans, ndshare, fdr=0.01) +
    scale_fill_party_c("Democratic\nvote share\nversus simulations\n",
                       labels = label_party_margin(midpoint=0.0),
                       midpoint = 0, limits=c(-0.2, 0.2))
p2 = plot_contr_fdr(map, plans, ndshare, fdr=0.05) +
    scale_fill_party_c("Democratic\nvote share\nversus simulations\n",
                       labels = label_party_margin(midpoint=0.0),
                       midpoint = 0, limits=c(-0.2, 0.2))
p = (p1 + p2 &
         theme(plot.margin=margin(),
               legend.direction="vertical",
               legend.key.height=unit(1.0, "cm"),
               legend.key.width=unit(0.6, "cm"))
     ) + plot_layout(guides="collect")
ggsave(here("paper/figures/fdr_ctrl.pdf"), plot=p, width=11, height=3.75, device=cairo_pdf)


# FWER Appendix ----------
fwer_calc <- function(map, plans, x, n=10) {
    x = rlang::enquo(x)
    x_nm = rlang::as_label(x)
    m = proj_distr(plans, !!x) # V rows, 5000 cols
    m = m - rowMeans2(m)
    m_pv = rowRanks(abs(m), ties.method="max") / ncol(m)

    threshes = seq(0.1, 1, length.out=50)^2
    res = purrr::map(sample(N, n), function(i) {
        q = est_pfdr(m_pc_pv[, i], m_pc_pv[, -i]) |>
            qvalues(m_pc_pv[, i])
        map_int(threshes, \(t) sum(q <= t))
    }, .progress=TRUE) |>
        do.call(cbind, args=_)

    tibble::new_tibble(list(
        stat = rep_along(threshes, x_nm),
        thresh = threshes,
        fwer = rowMeans(res > 0)
    ))
}

n_sim = 250
d_thresh = bind_rows(
    fwer_calc(map, plans, ndshare, n=n_sim),
    fwer_calc(map, plans, comp_polsby, n=n_sim),
    fwer_calc(map, plans, vap_white/total_vap, n=n_sim),
    fwer_calc(map, plans, pop_aian/total_pop, n=n_sim),
)

d_thresh |>
    mutate(
        stat = c(
            ndshare = "Normal Democratic vote share",
            comp_polsby = "Polsby-Popper compactness",
            `pop_aian/total_pop` = "American Indian/Alaska Native population",
            `vap_white/total_vap` = "White voting-age population"
        )[stat]
    ) |>
ggplot(aes(thresh, fwer)) +
    facet_wrap(~ stat) +
    geom_abline(slope=1, color="#445566", linetype="dashed", linewidth=0.3) +
    geom_line() +
    labs(x="pFDR threshold", y="Family-wise Error Rate (FWER)") +
    coord_cartesian(expand=FALSE) +
    theme_bw(base_size=10, base_family="Helvetica") +
    theme(plot.margin=margin(t=2, b=2, l=2, r=10))
ggsave(here("paper/figures/fwer_null.pdf"), width=7, height=5)

# Other quantities appendix ----------
app_plot_fdr <- function(qty, name, ...) {
    plot_contr_fdr(map, plans, {{ qty }}, norm=T) +
        geom_district(aes(group=cd_2020), fill=NA, linewidth=0.3) +
        scale_contr(paste0(name, "\n"), ...) +
        labs(title=name)
}
p = list(
    app_plot_fdr(pop_aian/total_pop, "American Indian/Alaska Native population", labels=percent),
    app_plot_fdr(comp_polsby, "Polsby-Popper compactness"),
    app_plot_fdr(vap_white/total_vap, "White voting-age population", labels=percent),
    guide_area()
) |>
    wrap_plots(guides="collect")
ggsave(here("paper/figures/fdr_addl_ex.pdf"), width=10, height=8)

# Re-do Shuffled figure with FDR control ----------
p = purrr::map(idx, function(i) {
    plot_contr_fdr(map, plans, ndshare, draw=i, density=0.15, spacing=0.025, norm=TRUE) +
        scale_fill_party_c("Z-score of\nDemocratic\nvote share\nversus simulations\n",
                           labels = label_number(accuracy=0.1),
                           midpoint = 0, limits=c(-4, 4)) +
        theme_paper()
}) |>
    wrap_plots(nrow=2) +
    plot_layout(guides="collect") +
    plot_annotation(theme = theme(plot.margin = margin())) &
    theme(legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"),
          legend.key.width=unit(0.6, "cm"))
ggsave(here("paper/figures/contr_samp_norm.pdf"), plot=p, width=9.5, height=3.5)

