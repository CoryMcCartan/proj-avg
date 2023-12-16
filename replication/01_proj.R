map = alarm_50state_map("OR", year=2020)
map$geometry = rmapshaper::ms_simplify(map$geometry, keep=0.06, keep_shapes=TRUE)
plans = bind_cols(
    alarm_50state_plans("OR", stats=FALSE, year=2020),
    select(alarm_50state_stats("OR", year=2020), -draw:-pop_overlap)
)

# Workflow fig -------------

# Ensemble
idx_plot = which(!duplicated(as.matrix(plans), MARGIN=2))[1:6]
p1 = lapply(idx_plot, function(i) {
    out = ggplot(map, aes(group = as.matrix(plans)[, i])) +
        geom_district(linewidth=0.12, is_coverage=TRUE) +
        scale_fill_penn82() +
        theme_paper()
    if (i == 1) {
        out = out +
            labs(title="Enacted plan") +
            theme(
                plot.title = element_text(
                    hjust=0.5, size=8, face="bold", margin=margin(t=4, b=-2)
                ),
                plot.background = element_rect(
                    fill="#e7e7e7", color=NA, linewidth=0
                )
            )
    }
    out
}) |>
    wrap_plots(ncol=2, nrow=3) +
    plot_annotation(theme = theme(plot.margin = margin()))
ggsave(here("paper/figures/ensemble.pdf"), plot=p1, width=2.5, height=2.9)

# Boxplots
p2 = redist.plot.distr_qtys(
    plans, ndshare, geom=\(...) stat_boxplot(...),
    coef=100, outlier.shape=NA, linewidth=0.3
) +
    labs(x="Districts, ordered by Democratic vote") +
    scale_y_continuous("Normal Democratic vote share", labels=percent) +
    guides(color="none") +
    theme_bw(base_family="Helvetica")
ggsave(here("paper/figures/distr_dshare.pdf"), plot=p2, width=3.2, height=3.52)


# Projections
layout = lapply(1:4, function(i) {
    area(t=5-i, l=i, b=11-i, r=6+i)
}) |>
    do.call(c, args=_)
p3 = lapply(1:4, function(i) {
    cd_ex = as.matrix(plans)[, i + 1]
    plot(map, plans$ndshare[6*i + 1:6][cd_ex]) +
        geom_district(aes(group=cd_ex), fill=NA, linewidth=0.6) +
        scale_fill_party_c("Democratic vote share\n", limits=c(0.36, 0.64),
                           breaks=c(0.4, 0.5, 0.6)) +
        coord_sf(expand=FALSE) +
        theme_paper() +
        theme(plot.background = element_rect(fill="#ffffffcc", color="black", linewidth=0.8),
              plot.margin = margin(0, 8, 8, 8))
}) |>
    rev() |>
    wrap_plots(design=layout, heights=c(rep(2, 4), rep(1, 6)), widths=rep(1, 10))
ggsave(here("paper/figures/proj_stacked.pdf"), plot=p3, width=6.75, height=6.25)

p4 = plot(map, proj_avg(plans, ndshare)) +
    geom_district(aes(group=cd_2020), fill=NA, linewidth=0.4) +
    scale_fill_party_c("Normal Democratic\nvote share\n",
                       limits=c(0.4, 0.6), breaks=c(0.4, 0.5, 0.6)) +
    theme_paper()
ggsave(here("paper/figures/proj_avg.pdf"), plot=p4, width=4.25, height=4)

p5 = plot(map, proj_contrast(plans, ndshare)) +
    geom_district(aes(group=cd_2020), fill=NA, linewidth=0.4) +
    scale_fill_party_c("Democratic vote share\nversus simulations\n",
                       labels = label_party_margin(midpoint=0.0),
                       midpoint = 0, limits=c(-0.2, 0.2)) +
    theme_paper()
ggsave(here("paper/figures/proj_contr.pdf"), plot=p5, width=4.25, height=4)


