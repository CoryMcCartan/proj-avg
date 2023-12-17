V = nrow(map)
E = sum(lengths(map$adj))
N = 5000
n = attr(plans, "ndists")

# Sampled projective contrasts fig --------
idx_shuf = sample(2:8)
idx = c(idx_shuf[1:3], 1, idx_shuf[4:7])
p = lapply(idx, function(i) {
    plot(map, proj_contrast(plans, ndshare, comp=i)) +
        scale_fill_party_c("Democratic vote share\nversus simulations\n",
                           labels = label_party_margin(midpoint=0.0),
                           midpoint = 0, limits=c(-0.2, 0.2)) +
        theme_paper()
}) |>
    wrap_plots(nrow=2) +
    plot_layout(guides="collect") +
    plot_annotation(theme = theme(plot.margin = margin())) &
    theme(legend.direction="vertical")
ggsave(here("paper/figures/contr_samp.pdf"), plot=p, width=9.5, height=3.25)

# check r-hats
rhats = sapply(1:V, function(i) {
    redist:::diag_rhat(pc_enac[i] - m[i, ], rep(1:2, each=N/2), split=FALSE)
})
hist(rhats)
if (any(rhats > 1.05)) stop("Projective contrast distribution did not converge.")

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

# Algorithm 1 of Storey & Tibshirani (2001)
est_pfdr <- function(p, ref=m_pc_pv, lambda=0.5, n_alpha=15) {
    alphas = quantile(p, probs=seq(0, 1, length.out=n_alpha), type=1, names=FALSE)
    m = length(p)
    pi0 = (m - sum(p <= lambda)) / (m - mean(colSums2(ref <= lambda)))
    pi0 = min(pi0, 1)
    fdrs = vapply(alphas, \(a) {
        pi0 *  mean(colSums2(ref <= a)) / sum(p <= a)
    }, numeric(1))
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
qvalues <- function(ests, p) {
    splinefun(ests$alpha, rev(cummin(rev(ests$pfdr))), method="monoH.FC")(p)
}


plot_contr_fdr <- function(map, plans, x, fdr=0.05, draw=1, norm=FALSE) {
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
        geom_sf(data=sf::st_buffer(geom_sel, -1600), inherit.aes=FALSE,
                fill=NA, linewidth=0.7, color="#ffffff88") +
        geom_sf(data=geom_sel, inherit.aes=FALSE,
                fill=NA, linewidth=0.2, color="black")
}


p1 = plot_contr_fdr(map, plans, ndshare) +
    scale_fill_party_c("Democratic vote share\nversus simulations\n",
                       labels = label_party_margin(midpoint=0.0),
                       midpoint = 0, limits=c(-0.2, 0.2))
p2 = plot_contr_fdr(map, plans, ndshare, norm=TRUE) +
    scale_fill_party_c("Z-score of Democratic\nvote share versus simulations\n",
                       labels = label_number(accuracy=0.1),
                       midpoint = 0, limits=c(-4, 4))
p = p1 + p2 & theme(plot.margin=margin())
ggsave(here("paper/figures/fdr_norm.pdf"), plot=p, width=9, height=4)


ests = est_pfdr(pc_pv_enac)
plot(ests, pch=16, type='b', ylim=0:1); abline(a=0, b=1, col="red")
plot(pc_pv_enac, qvalues(ests, pc_pv_enac), cex=0.1, ylim=0:1); abline(a=0, b=1, col="red")

i = 425
ests = est_pfdr(m_pc_pv[, i], m_pc_pv[, -i])
plot(m_pc_pv[, i], qvalues(ests, m_pc_pv[, i]), cex=0.1, ylim=0:1); abline(a=0, b=1, col="red")

threshes = seq(0.1, 1, length.out=16)^2
res = purrr::map(sample(N, 100), function(i) {
    ests = est_pfdr(m_pc_pv[, i], m_pc_pv[, -i])
    sapply(threshes, \(t) sum(m_pc_pv[, i] <= find_thresh(ests, t)))
}, .progress=TRUE) |>
    do.call(cbind, args=_)
plot(threshes, rowMeans(res > 0)); abline(a=0, b=1, col="red")
