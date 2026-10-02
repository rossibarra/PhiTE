#!/usr/bin/env Rscript

# Deterministic versions of the two Phi-SFS explanatory figures (PDF, and PNG
# for the null figure).
# Run from the repository root with:
#   module load R/4.4.2
#   Rscript figures/make_phi_sfs_figures.R

orange <- "#E66101"
orange_fill <- "#F3C6AA"
teal <- "#257D7D"
teal_fill <- "#B9D7D7"
neutral <- "#9E9E9E"
neutral_fill <- "#D0D0D0"
ink <- "#202020"

normalize <- function(x) x / sum(x)

step_coordinates <- function(x, y) {
    list(
        x = c(x[1], rep(x[-1], each = 2)),
        y = c(y[1], as.vector(rbind(y[-length(y)], y[-1])))
    )
}

draw_sfs <- function(x, neutral_sfs, te_sfs, color, title, show_y = FALSE) {
    ymax <- max(c(neutral_sfs, te_sfs)) * 1.08
    plot(NA, xlim = c(0.02, 0.98), ylim = c(0, ymax), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")

    # The wider neutral bars remain a single uniform gray. Narrower opaque TE
    # bars are drawn in front, avoiding transparency-induced color changes.
    rect(x - 0.022, 0, x + 0.022, neutral_sfs,
         col = neutral_fill, border = neutral, lwd = 2.3)
    rect(x - 0.015, 0, x + 0.015, te_sfs,
         col = adjustcolor(color, alpha.f = 0.58), border = color, lwd = 3.1)

    axis(1, at = c(0.05, 0.95), labels = c("0.05", "0.95"),
         lwd = 3.0, lwd.ticks = 3.0, cex.axis = 1.45)
    if (show_y) {
        axis(2, at = c(0, ymax), labels = c("0", ""),
             lwd = 3.0, lwd.ticks = 3.0, cex.axis = 1.4)
    }
    box(bty = "l", lwd = 3.0)
    mtext("DAF", side = 1, line = 2.2, cex = 1.6, font = 2)
    title(main = title, col.main = color, cex.main = 1.55, font.main = 2,
          line = 0.45)
    legend("topright", legend = c("SNPs", "TEs"),
           fill = c(neutral_fill, adjustcolor(color, alpha.f = 0.58)),
           border = c(neutral, color), bty = "n", cex = 1.4,
           text.font = 2, x.intersp = 0.6, y.intersp = 0.9)
}

draw_sidebar <- function(label) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    rect(0.06, 0, 0.94, 1, col = "#F0F0F0", border = NA)
    text(0.5, 0.5, label, srt = 90, cex = 2.0, font = 2, col = ink)
}

draw_cdf <- function(x, neutral_sfs, te_sfs, color, fill,
                     show_y = FALSE, annotation_x = 0.48,
                     annotation_y = 0.34, arrow_x = 0.34) {
    neutral_cdf <- cumsum(neutral_sfs)
    te_cdf <- cumsum(te_sfs)
    neutral_step <- step_coordinates(x, neutral_cdf)
    te_step <- step_coordinates(x, te_cdf)

    plot(NA, xlim = c(0.02, 0.98), ylim = c(0, 1.02), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    polygon(c(neutral_step$x, rev(te_step$x)),
            c(neutral_step$y, rev(te_step$y)),
            col = fill, border = NA)
    lines(neutral_step$x, neutral_step$y, col = neutral, lwd = 4.1)
    lines(te_step$x, te_step$y, col = color, lwd = 4.4)

    axis(1, at = c(0.05, 0.95), labels = c("0.05", "0.95"),
         lwd = 3.0, lwd.ticks = 3.0, cex.axis = 1.45)
    if (show_y) {
        axis(2, at = c(0, 1), labels = c("0", "1"),
             lwd = 3.0, lwd.ticks = 3.0, cex.axis = 1.4)
    }
    box(bty = "l", lwd = 3.0)
    mtext("DAF", side = 1, line = 2.2, cex = 1.6, font = 2)

    target_bin <- max(which(x <= arrow_x))
    target_y <- (te_cdf[target_bin] + neutral_cdf[target_bin]) / 2
    text(annotation_x, annotation_y, expression(Phi[SFS]),
         cex = 1.75, font = 2, col = ink)
    arrows(annotation_x - 0.03, annotation_y + 0.08,
           arrow_x, target_y, length = 0.11, lwd = 3.2, col = ink)
}

draw_vertical_label <- function(label, cex = 1.5) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.55, 0.5, label, srt = 90, cex = cex, font = 2, xpd = NA)
}

make_definition_figure <- function(path) {
    bins <- 1:19
    x <- bins / 20
    neutral_sfs <- normalize(1 / bins)
    rare_sfs <- normalize(1 / bins^1.55)
    high_sfs <- normalize(0.55 / bins + 0.055 * exp(0.42 * (bins - 12)))

    pdf(path, width = 11.2, height = 8.6, family = "Helvetica",
        useDingbats = FALSE)
    layout(matrix(1:8, nrow = 2, byrow = TRUE),
           heights = c(1, 1), widths = c(0.10, 0.14, 1, 1))

    par(mar = rep(0, 4))
    draw_sidebar("SFS")
    par(mar = rep(0, 4))
    draw_vertical_label("Proportion")

    par(mar = c(3.8, 3.2, 2.5, 1.0), mgp = c(2.1, 0.65, 0),
        tcl = -0.35, las = 1)

    draw_sfs(x, neutral_sfs, rare_sfs, orange,
             "Excess rare variants", show_y = TRUE)
    draw_sfs(x, neutral_sfs, high_sfs, teal,
             "Excess high-frequency derived variants", show_y = FALSE)

    par(mar = rep(0, 4))
    draw_sidebar("CDF")
    par(mar = rep(0, 4))
    draw_vertical_label("Cumulative proportion", cex = 1.35)

    par(mar = c(3.8, 3.2, 0.8, 1.0), mgp = c(2.1, 0.65, 0),
        tcl = -0.35, las = 1)
    draw_cdf(x, neutral_sfs, rare_sfs, orange, orange_fill,
             show_y = TRUE, annotation_x = 0.49,
             annotation_y = 0.31, arrow_x = 0.30)
    draw_cdf(x, neutral_sfs, high_sfs, teal, teal_fill,
             show_y = FALSE, annotation_x = 0.69,
             annotation_y = 0.31, arrow_x = 0.57)
    dev.off()
}

# Null figure: Phi-SFS on its own scale, one category per column. The gray
# violin is the category's null distribution of Phi_i^0 (two neutral SNP sets of
# the same size M), so its height is the finite-sample floor, which shrinks as
# M grows. The point is Phi_obs, colored by its add-one P-value; the segment
# from the null mean mu_0 to Phi_obs is the effect size, Phi_obs - mu_0.
# Values are synthetic. Null means and SDs scale as 1/sqrt(M), anchored to the
# in-gene production null (mean 0.0050, SD 0.0018 at M = 4,067).
draw_violin <- function(values, center, width = 0.30) {
    d <- density(values, from = min(values), to = max(values), n = 512,
                 cut = 0, bw = "nrd0")
    half_width <- width * d$y / max(d$y)
    polygon(c(center - half_width, rev(center + half_width)),
            c(d$x, rev(d$x)), col = "#D0D0D0", border = "#B5B5B5",
            lwd = 1.8)
}

draw_null_figure <- function() {
    categories <- c("In gene", "0-2 kb", "2-5 kb", ">5 kb")
    site_count <- c(4000, 2500, 1200, 600)
    observed <- c(0.0300, 0.0185, 0.0150, 0.0160)
    null_count <- 500
    u <- ((1:null_count) - 0.5) / null_count
    shape <- (0.0050 / 0.0018)^2
    nulls <- lapply(site_count, function(m) {
        mu <- 0.0050 * sqrt(4067 / m)
        qgamma(u, shape = shape, rate = shape / mu)
    })
    p_values <- mapply(function(null, obs) (1 + sum(null >= obs)) / (length(null) + 1),
                       nulls, observed)
    p_cap <- log10(null_count + 1)
    minus_log10_p <- pmin(-log10(p_values), p_cap)
    blue <- colorRampPalette(c("#DEEBF7", "#6BAED6", "#08519C"))(301)
    point_colors <- blue[1 + round(minus_log10_p / p_cap * 300)]
    ylim <- c(0, 0.034)

    layout(matrix(c(1, 2, 3), nrow = 1), widths = c(0.38, 4.7, 1.35))
    par(oma = c(0, 0, 3.2, 0))
    par(mar = rep(0, 4))
    draw_vertical_label(expression(paste(Phi[SFS], " (DAF units)")), cex = 1.3)

    par(mar = c(6.2, 4.6, 1.0, 0.6), mgp = c(3.2, 0.75, 0),
        tcl = -0.35, las = 1)
    plot(NA, xlim = c(0.45, 4.55), ylim = ylim, axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    for (i in seq_along(nulls)) {
        draw_violin(nulls[[i]], i)
        mu <- mean(nulls[[i]])
        segments(i - 0.16, mu, i + 0.16, mu, col = "#6F6F6F", lwd = 2.6)
        segments(i, mu, i, observed[i], col = "#79A6D2", lwd = 3.1)
    }
    points(1:4, observed, pch = 21, bg = point_colors,
           col = "white", lwd = 2.0, cex = 2.4)
    axis(1, at = 1:4, labels = categories, lwd = 2.8,
         lwd.ticks = 2.8, cex.axis = 1.25)
    mtext(paste0("M = ", formatC(site_count, format = "d", big.mark = ",")), side = 1,
          at = 1:4, line = 2.4, cex = 0.95, col = "#5A5A5A")
    axis(2, at = seq(0, 0.03, by = 0.01), lwd = 2.8,
         lwd.ticks = 2.8, cex.axis = 1.25)
    box(bty = "l", lwd = 2.8)
    mtext("Distance to nearest gene", side = 1, line = 4.4,
          cex = 1.45, font = 2)

    par(mar = c(6.2, 0.4, 1.0, 0.5))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.5, 0.91, expression(-log[10](italic(P)*"-value")),
         cex = 1.3, font = 2)
    edges <- seq(0.10, 0.90, length.out = length(blue) + 1)
    rect(edges[-length(edges)], 0.81, edges[-1], 0.865,
         col = blue, border = NA)
    rect(0.10, 0.81, 0.90, 0.865, border = ink, lwd = 1.8)
    tick_value <- c(0, 1, 2, p_cap)
    tick_x <- 0.10 + tick_value / p_cap * 0.80
    segments(tick_x, 0.785, tick_x, 0.81, lwd = 1.7)
    text(tick_x, 0.75, labels = c("0", "1", "2", "2.7"), cex = 1.1)
    rect(0.12, 0.61, 0.24, 0.67, col = "#D0D0D0",
         border = "#B5B5B5", lwd = 1.8)
    text(0.31, 0.64, "Null distribution", adj = 0, cex = 1.2, font = 2)
    segments(0.12, 0.53, 0.24, 0.53, col = "#6F6F6F", lwd = 2.6)
    text(0.31, 0.53, expression(paste("Null mean ", mu[0])), adj = 0,
         cex = 1.2, font = 2)
    segments(0.18, 0.39, 0.18, 0.46, col = "#79A6D2", lwd = 3.1)
    text(0.31, 0.425, expression(Phi[obs] - mu[0]), adj = 0,
         cex = 1.2, font = 2)
    text(0.31, 0.375, "effect size", adj = 0, cex = 1.05)

    mtext("Illustrative example", side = 3, line = 1.0,
          outer = TRUE, cex = 1.85, font = 2)
    invisible(p_values)
}

make_null_figure <- function(stem) {
    pdf(paste0(stem, ".pdf"), width = 10.2, height = 7.2,
        family = "Helvetica", useDingbats = FALSE)
    p_values <- draw_null_figure()
    dev.off()
    png(paste0(stem, ".png"), width = 10.2, height = 7.2, units = "in",
        res = 180, family = "Helvetica")
    draw_null_figure()
    dev.off()
    p_values
}

make_definition_figure("figures/phi_sfs_definition_schematic.pdf")
print(make_null_figure("figures/phi_sfs_null_example"))
