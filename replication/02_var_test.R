V = nrow(map)
E = sum(lengths(map$adj))
N = 5000

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

# Multiple testing experiments ------------
m = proj_distr(plans, ndshare) # V rows, 5000 cols
p_avg = rowMeans2(m)
m_pc = m - p_avg
pc_enac = plans$ndshare[1:attr(plans, "ndists")][as.matrix(plans)[, 1]] - p_avg

moran = geomander::global_morans(adj=map$adj, wts=p_avg)$moran / E
# moran = geomander::global_morans(adj=map$adj, wts=pc_enac)$moran / E
n_eff = V * (1 - moran) / (1 + moran)

n_effs = sample(5000, 100) |>
    sapply(function(i) {
        moran = geomander::global_morans(adj=map$adj, wts=m_pc[, i])$moran / E
        V * (1 - moran) / (1 + moran)
    })
hist(n_effs)

m_pc_pv = rowRanks(abs(m_pc), ties.method="max") / (ncol(m) + 1)
pc_pv_enac = rowMeans2(abs(m_pc) >= abs(pc_enac))
pc_z_enac = abs(pc_enac) / rowSds(m_pc)
m_pc_pv_cor = cor(t(m_pc_pv))
hist(m_pc_pv_cor, breaks=100)

fwer <- function(rule) {
    rule = rlang::as_function(rule)
    mean(apply(m_pc_pv, 2, \(x) any(rule(x))))
}

co = prec_cooccurrence(plans, sampled_only=TRUE, ncores=6)
mds = cmdscale(1 - co)
cl = hclust(as.dist(1 - co))

cl2 = cbind(sf::st_coordinates(sf::st_centroid(map$geometry)), pc_enac, p_avg) |>
    scale() |>
    dist() |>
    hclust()
plot(cl2, lwd=0.5, labels=FALSE)

plot(map, cutree(cl2, k=8))
plot(map, as.factor(cutree(cl2, k=30)))
plot(mds, pch=16, cex=1.0,
     col=scales::alpha(scales::colour_ramp(wacolors$sound_sunset)(-log10(pc_pv_enac)/3.7), alpha=0.2))

X = cbind(sf::st_coordinates(sf::st_centroid(map$geometry)), pc_enac, p_avg, rowMeans2(as.matrix(plans))) |>
    scale()
x = cluster::pam(X, k=25)$clustering / 20
# plot(mds, pch=16, cex=1.0,
#      col=scales::alpha(scales::colour_ramp(wacolors$sound_sunset)(x), alpha=0.2))
plot(map, as.factor(x))
    geom_district(aes(group=cd_2020), fill=NA, linewidth=0.4)

alpha = 0.05
alphas = c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15,
           0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

fwer_plot <- function(pv = pc_pv_enac, adj_proc = ~., alpha = alphas) {
    adj_proc = rlang::as_function(adj_proc)
    m_pv = apply(m_pc_pv, 2, adj_proc)
    pv = adj_proc(pv)
    fwers = sapply(alpha, \(a) mean(colSums(m_pv <= a) > 0))
    x = seq(0.05, 0.95, 0.05)
    y = sapply(x, \(t) (V - sum(pv <= t)) / (V - mean(colSums(m_pv <= t))))
    pi0 = median(y)
    print(pi0)
    print(y)
    fdrs = sapply(alphas, \(a) {
       pi0 *  mean(colSums(m_pv <= a)) / sum(pv <= a)
    })
    plot(alpha, fwers, ylim=0:1, type='b', pch=16, cex=0.8)
    points(alphas, fdrs, ylim=0:1, type='b', pch=17, cex=0.8, col="#305070")
    abline(a=0, b=1, col="red")
}

fwer_plot(pc_pv_enac)
fwer_plot(m_pc_pv[, 425])
fwer_plot(pc_pv_enac, ~ . * 4)
fwer_plot(pc_pv_enac, ~ . * n_eff)
fwer_plot(pc_pv_enac, ~ . * 100)
fwer_plot(pc_pv_enac, ~ . * mean(n_effs))
fwer_plot(pc_pv_enac, ~ p.adjust(., "BH")^1.0)
fwer_plot(pc_pv_enac, ~ p.adjust(., "BH")^1.5)

plot(sort(p.adjust(pc_pv_enac, "BH")^1.5), cex=0.2, ylim=0:1)
plot(sort(p.adjust(m_pc_pv[, 100], "BH")^1.5), cex=0.2, ylim=0:1)

plot(sort(m_pc_pv[, 100]), ylim=0:1, type='n')
for(i in sample(5000, 100)) lines(sort(m_pc_pv[, i]), ylim=0:1, lwd=0.1)
lines(sort(pc_pv_enac), ylim=0:1, col='red', lwd=1.5)
lines(rowMeans(apply(m_pc_pv, 2, sort)), ylim=0:1, col='blue', lwd=1.5)
# plot(sort(p.adjust2(pc_pv_enac, method="BY", n=round(mean(n_effs)))), ylim=0:1)

plot(map, p.adjust(pc_pv_enac, method="BH")^1.5 <= 0.1)
plot(map, pc_pv_enac*25 <= 0.1)
plot(map, p.adjust(m_pc_pv[, 25], method="BH")^2 <= 0.1)

apply(m_pc_pv, 2, \(x) any(x < alpha)) |> mean()
hist(m_pc_pv)

x = seq(0.05, 0.95, 0.05)
y = sapply(x, \(t) (V - sum(pc_pv_enac <= t)) / (V - mean(colSums(m_pc_pv <= t))))
pi0 = median(y)
plot(x, y, type='b')

n_rej = colSums(m_pc_pv <= 0.05)
n_rej_enac = sum(pc_pv_enac <= 0.05)
mean(n_rej / (n_rej + n_rej_enac))

pi0 * mean(colSums(m_pc_pv <= 0.05)) / sum(pc_pv_enac <= 0.05)


fwers = sapply(alphas, \(a) mean(colSums(m_pc_pv <= a) > 0))
plot(alphas, fwers, ylim=0:1, type='b', pch=16, cex=0.8)

abline(a=0, b=1, col="red")

geom_sel = map |>
    filter(p.adjust(pc_pv_enac, method="BH")^1.5 < 0.10) |>
    summarize() |>
    sf::st_buffer(-1200) |>
    suppressWarnings()

plot(map, pc_enac / rowSds(m_pc)) +
    geom_sf(data=sf::st_buffer(geom_sel, -1600), inherit.aes=FALSE,
            fill=NA, linewidth=0.8, color="#ffffff88") +
    geom_sf(data=geom_sel, inherit.aes=FALSE,
            # fill=NA, linewidth=0.30, color="yellow")
            fill=NA, linewidth=0.25, color="black") +
    # geom_district(aes(group=cd_2020), fill=NA, linewidth=0.4) +
    scale_fill_party_c("Democratic vote share\nversus simulations\n",
                       labels = label_party_margin(midpoint=0.0),
                       # midpoint = 0, limits=c(-0.2, 0.2)) +
                       midpoint = 0, limits=c(-4, 4)) +
    theme_paper()
