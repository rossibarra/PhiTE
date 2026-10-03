#!/usr/bin/env Rscript

# Graphical abstract of the PhiTE pipeline (PDF and PNG). All values are
# synthetic and only illustrate the steps.
# Run from the repository root with:
#   module load R/4.4.2
#   Rscript figures/make_graphical_abstract.R

orange <- "#E66101"
orange_fill <- "#F3C6AA"
neutral <- "#9E9E9E"
neutral_fill <- "#D0D0D0"
blue <- "#79A6D2"
dark_blue <- "#08519C"
ink <- "#202020"
soft <- "#5A5A5A"

step_coordinates <- function(x, y) {
    list(
        x = c(x[1], rep(x[-1], each = 2)),
        y = c(y[1], as.vector(rbind(y[-length(y)], y[-1])))
    )
}

panel_title <- function(number, title) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    symbols(0.07, 0.45, circles = 0.05, inches = FALSE, add = TRUE,
            bg = ink, fg = ink)
    text(0.07, 0.45, number, col = "white", font = 2, cex = 1.5)
    text(0.15, 0.45, title, adj = 0, font = 2, cex = 1.6, col = ink)
}

panel_caption <- function(lines) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    y <- 0.92 - (seq_along(lines) - 1) * 0.26
    for (i in seq_along(lines)) {
        text(0.04, y[i], lines[[i]], adj = c(0, 1), cex = 1.18, col = soft)
    }
}

flow_arrow <- function() {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    polygon(c(0.1, 0.55, 0.55, 0.92, 0.55, 0.55, 0.1),
            c(0.46, 0.46, 0.40, 0.5, 0.60, 0.54, 0.54),
            col = "#BDBDBD", border = NA)
}

# 1. Posterior ARGs: three offset copies of one genealogy stand for the
# posterior draws; the mutation's branch gives the age interval.
draw_tree <- function(dx, dy, col, lwd, mutation = FALSE) {
    tips <- c(0.10, 0.24, 0.38, 0.52, 0.66, 0.80) + dx
    base <- 0.10 + dy
    n12 <- 0.34 + dy; n34 <- 0.28 + dy; n56 <- 0.40 + dy
    n1234 <- 0.55 + dy; root <- 0.80 + dy
    x12 <- mean(tips[1:2]); x34 <- mean(tips[3:4]); x56 <- mean(tips[5:6])
    x1234 <- mean(c(x12, x34)); xroot <- mean(c(x1234, x56))
    seg <- function(x0, y0, x1, y1) segments(x0, y0, x1, y1, col = col, lwd = lwd)
    for (i in 1:2) seg(tips[i], base, tips[i], n12)
    for (i in 3:4) seg(tips[i], base, tips[i], n34)
    for (i in 5:6) seg(tips[i], base, tips[i], n56)
    seg(tips[1], n12, tips[2], n12); seg(tips[3], n34, tips[4], n34)
    seg(tips[5], n56, tips[6], n56)
    seg(x12, n12, x12, n1234); seg(x34, n34, x34, n1234)
    seg(x12, n1234, x34, n1234)
    seg(x1234, n1234, x1234, root); seg(x56, n56, x56, root)
    seg(x1234, root, x56, root)
    if (mutation) {
        points(x12, (n12 + n1234) / 2, pch = 21, bg = orange, col = "white",
               cex = 2.6, lwd = 1.5)
        list(x = x12, lo = n12, hi = n1234)
    }
}

draw_args <- function() {
    par(mar = c(1.2, 1.2, 0.6, 0.6))
    plot(NA, xlim = c(0, 1.08), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    draw_tree(0.12, 0.12, "#DADADA", 2.4)
    draw_tree(0.06, 0.06, "#B8B8B8", 2.6)
    m <- draw_tree(0, 0, ink, 3.0, mutation = TRUE)
    # Age interval of the mutation's branch, and the time axis.
    segments(0.045, c(m$lo, m$hi), m$x, c(m$lo, m$hi), col = orange, lwd = 1.6, lty = 2)
    segments(0.03, m$lo, 0.03, m$hi, col = orange, lwd = 5)
    segments(0.015, c(m$lo, m$hi), 0.045, c(m$lo, m$hi), col = orange, lwd = 3)
    text(0.03, m$hi + 0.06, "age", col = orange, font = 2, cex = 1.15)
    arrows(-0.02, 0.10, -0.02, 0.98, length = 0.08, lwd = 2, col = soft, xpd = NA)
    text(-0.05, 0.55, "time", srt = 90, col = soft, cex = 1.1, xpd = NA)
    text(0.45, 0.03, "TEs and SNPs are sites in the same ARGs",
         cex = 1.05, col = soft)
}

# 2. Age-matched controls: the focal age distribution and matched SNP sets.
draw_matching <- function() {
    par(mar = c(3.6, 1.2, 0.6, 0.8), mgp = c(2.2, 0.6, 0), tcl = -0.3)
    grid <- seq(3, 6.5, length.out = 15)
    width <- diff(grid)[1]
    focal <- dnorm(grid, mean = log10(4.5e4), sd = 0.55)
    focal <- focal / sum(focal)
    plot(NA, xlim = range(grid) + c(-width, width) / 2, ylim = c(0, max(focal) * 1.35),
         axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    rect(grid - width / 2, 0, grid + width / 2, focal,
         col = orange_fill, border = orange, lwd = 2.2)
    set.seed(7)
    for (i in 1:6) {
        boot <- tabulate(sample(seq_along(grid), 900, replace = TRUE, prob = focal),
                         nbins = length(grid)) / 900
        s <- step_coordinates(c(grid - width / 2, max(grid) + width / 2),
                              c(boot, boot[length(boot)]))
        lines(s$x, s$y, col = adjustcolor(neutral, 0.9), lwd = 1.8)
    }
    axis(1, at = 3:6, labels = expression(10^3, 10^4, 10^5, 10^6),
         lwd = 2.4, lwd.ticks = 2.4, cex.axis = 1.15)
    box(bty = "l", lwd = 2.4)
    mtext("Allele age (generations)", side = 1, line = 2.4, cex = 1.05, font = 2)
    legend("topright", legend = c("A: focal TEs", expression(B[0] * ", ..., " * B[R] * ": SNPs")),
           fill = c(orange_fill, NA), border = c(orange, NA),
           lty = c(NA, 1), lwd = c(NA, 2), col = c(NA, neutral),
           bty = "n", cex = 1.05, x.intersp = 0.6, seg.len = 1.2)
}

# 3. Posterior-polarized SFS: CDFs of A and B0; Phi_SFS is the area between.
draw_phi <- function() {
    par(mar = c(3.6, 1.2, 0.6, 0.8), mgp = c(2.2, 0.6, 0), tcl = -0.3)
    bins <- 1:19
    x <- bins / 20
    b0 <- cumsum(1 / bins) / sum(1 / bins)
    a <- cumsum(1 / bins^1.45) / sum(1 / bins^1.45)
    sa <- step_coordinates(x, a)
    sb <- step_coordinates(x, b0)
    plot(NA, xlim = c(0.03, 0.97), ylim = c(0, 1.02), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    polygon(c(sb$x, rev(sa$x)), c(sb$y, rev(sa$y)), col = orange_fill, border = NA)
    lines(sb$x, sb$y, col = neutral, lwd = 3.4)
    lines(sa$x, sa$y, col = orange, lwd = 3.6)
    axis(1, at = c(0.05, 0.5, 0.95), labels = c("0.05", "0.5", "0.95"),
         lwd = 2.4, lwd.ticks = 2.4, cex.axis = 1.15)
    box(bty = "l", lwd = 2.4)
    mtext("Derived allele frequency", side = 1, line = 2.4, cex = 1.05, font = 2)
    text(0.55, 0.42, expression(Phi[obs] == W[1](A, B[0])), cex = 1.35, col = ink)
    arrows(0.42, 0.48, 0.27, 0.70, length = 0.09, lwd = 2.4, col = ink)
    legend("bottomright", legend = c("A", expression(B[0])),
           col = c(orange, neutral), lwd = 3.4, bty = "n", cex = 1.1,
           seg.len = 1.4, inset = c(0.02, 0.04))
}

# 4. Calibration: null distances between neutral sets, the floor mu_0, the
# observed distance and the floor-corrected effect size.
draw_null <- function() {
    par(mar = c(3.6, 1.2, 0.6, 0.8), mgp = c(2.2, 0.6, 0), tcl = -0.3)
    u <- (1:500 - 0.5) / 500
    shape <- 8
    mu <- 0.0114
    null <- qgamma(u, shape = shape, rate = shape / mu)
    d <- density(null, n = 512)
    obs <- 0.031
    phi_hat <- sqrt(obs^2 - mu^2)
    plot(NA, xlim = c(0, 0.036), ylim = c(0, max(d$y) * 1.45), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    polygon(d$x, d$y, col = neutral_fill, border = "#B5B5B5", lwd = 1.8)
    ytop <- max(d$y)
    segments(mu, 0, mu, ytop * 1.08, col = "#6F6F6F", lwd = 2.6)
    text(mu, ytop * 1.17, expression(mu[0] * " (floor)"), col = "#6F6F6F", cex = 1.1)
    segments(phi_hat, 0, phi_hat, ytop * 0.62, col = blue, lwd = 3.4)
    text(phi_hat - 0.0006, ytop * 0.72, expression(hat(Phi)[SFS]), col = dark_blue,
         cex = 1.25, adj = 1)
    points(obs, ytop * 0.30, pch = 21, bg = dark_blue, col = "white", cex = 2.6, lwd = 1.5)
    text(obs, ytop * 0.47, expression(Phi[obs]), col = dark_blue, cex = 1.25)
    text(0.0215, ytop * 1.36, expression(Phi[i]^0 == W[1](B[i], B[0])), col = soft,
         cex = 1.1)
    axis(1, at = c(0, 0.01, 0.02, 0.03), lwd = 2.4, lwd.ticks = 2.4, cex.axis = 1.15)
    box(bty = "l", lwd = 2.4)
    mtext(expression(bold(Phi[SFS] ~ "(DAF units)")), side = 1, line = 2.4, cex = 1.05)
}

draw_abstract <- function() {
    layout(matrix(c(1, 0, 2, 0, 3, 0, 4,
                    5, 6, 7, 8, 9, 10, 11,
                    12, 0, 13, 0, 14, 0, 15), nrow = 3, byrow = TRUE),
           widths = c(1, 0.16, 1, 0.16, 1, 0.16, 1),
           heights = c(0.16, 1, 0.32))
    par(oma = c(0.4, 1.2, 4.2, 0.6), mar = rep(0, 4), family = "Helvetica")

    panel_title("1", "Posterior ARGs")
    panel_title("2", "Age-matched controls")
    panel_title("3", "Compare spectra")
    panel_title("4", "Calibrate the floor")

    draw_args(); par(mar = rep(0, 4)); flow_arrow()
    draw_matching(); par(mar = rep(0, 4)); flow_arrow()
    draw_phi(); par(mar = rep(0, 4)); flow_arrow()
    draw_null()

    par(mar = rep(0, 4))
    panel_caption(list(
        "Every site gets an age interval",
        "and q = P(ALT is derived) from",
        "the posterior ARG draws."))
    panel_caption(list(
        "Each control set matches M SNPs",
        "to a bootstrap of A's ages:",
        "500 disjoint sets, frequency-blind."))
    panel_caption(list(
        "Project spectra to m = 20 and mix",
        "both orientations by q.",
        expression(Phi[SFS] * " is the area between the CDFs.")))
    panel_caption(list(
        "Null sets give the sampling floor.",
        "Report the add-one P-value and",
        expression(hat(Phi)[SFS] == sqrt(Phi[obs]^2 - mu[0]^2))))

    mtext("PhiTE: does a TE category's frequency spectrum differ from neutral SNPs of the same age?",
          side = 3, outer = TRUE, line = 1.4, cex = 1.55, font = 2, col = ink)
}

make_abstract <- function(stem) {
    pdf(paste0(stem, ".pdf"), width = 17, height = 5.3, family = "Helvetica",
        useDingbats = FALSE)
    draw_abstract()
    dev.off()
    png(paste0(stem, ".png"), width = 17, height = 5.3, units = "in", res = 160,
        family = "Helvetica")
    draw_abstract()
    dev.off()
}

make_abstract("figures/phite_graphical_abstract")
