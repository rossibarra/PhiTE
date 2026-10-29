"""Plot every tested matching set against the accepted sets and the observed TEs.

Reads one bootstrap_target_matcher.py output directory (the product of a data
store and one TE target set) and writes two figures:

  *_age_cdf.{png,pdf}  age CDF of every tested set (grey), the selected set of
                       each QC-passing replicate (blue) and the observed TEs
                       (vermillion), with a second panel of each CDF minus the
                       observed CDF.
  *_w1.{png,pdf}       the same three groups, as empirical CDFs of the W1
                       distance each set reached.

"Tested sets" are the best state of every restart of every replicate
(restart_best_cdfs.npy / restart_match_to_bootstrap_w1.npy). "Accepted" is the
selected restart of each replicate with qc_pass True. A replicate that fails QC
stays in the grey group only.

In the W1 figure the third group is bootstrap_to_observed_w1: how far each
bootstrap target lies from the observed TEs, i.e. the sampling scale a match
has to stay well below. Matches are made to the bootstrap target, so no match
can be compared with the observed W1 directly (see match_to_observed_w1 for
that).

An unmatched baseline can be added in green. --random-dir takes the output of
tools/random_control_sets.py (random SNP sets from the same store, target and
candidates) and adds it to both figures. --random-w1 takes just a .npy of W1
values and adds it to the W1 figure only; use --random-label to say where the
values came from.

Usage:
    python tools/plot_match_sets.py MATCH_DIR [--out-dir figures/match_sets]
        [--name NAME] [--linear-age] [--age-max GENERATIONS]
        [--random-dir DIR | --random-w1 FILE.npy] [--random-label TEXT]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

# Okabe-Ito blue and vermillion are CVD-safe against each other and the grey.
TESTED = "#9aa0a6"
ACCEPTED = "#0072B2"
OBSERVED = "#D55E00"
RANDOM = "#009E73"
INK = "#222222"


def load(match_dir: Path) -> dict[str, np.ndarray]:
    names = [
        "age_bins", "target_cdf", "cdfs", "restart_best_cdfs", "qc_pass",
        "restart_match_to_bootstrap_w1", "match_to_bootstrap_w1",
        "bootstrap_to_observed_w1",
    ]
    data = {n: np.load(match_dir / f"{n}.npy", mmap_mode="r") for n in names}
    reps = data["qc_pass"].shape[0]
    if data["cdfs"].shape[0] != reps or data["restart_best_cdfs"].shape[0] != reps:
        raise ValueError("match arrays are not aligned on the replicate axis")
    return data


def style(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#888888")
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(True, color="#e6e6e6", linewidth=0.6)
    ax.set_axisbelow(True)


def curves(ax, x, ys, color, width, alpha, z):
    ax.add_collection(LineCollection(
        [np.column_stack([x, y]) for y in ys],
        colors=color, linewidths=width, alpha=alpha, zorder=z,
    ))


def legend(ax, handles, loc):
    ax.legend(handles=handles, loc=loc, frameon=False, fontsize=9, labelcolor=INK)


def age_cdf_figure(d, title, linear_age, age_max, random_cdfs=None):
    ages = np.asarray(d["age_bins"], dtype=float)
    keep = ages > 0 if not linear_age else np.ones(ages.size, bool)
    if age_max is not None:
        keep &= ages <= age_max
    x = ages[keep]
    qc = np.asarray(d["qc_pass"], dtype=bool)
    tested = np.asarray(d["restart_best_cdfs"][:, :, keep]).reshape(-1, x.size)
    accepted = np.asarray(d["cdfs"][qc][:, keep])
    observed = np.asarray(d["target_cdf"][keep])
    random_cdfs = None if random_cdfs is None else random_cdfs[:, keep]

    # With a random baseline the matched differences are too small to see at
    # the random sets' scale, so the difference gets a zoomed and a full panel.
    rows = 3 if random_cdfs is not None else 2
    fig, axes = plt.subplots(
        rows, 1, figsize=(7, 6.4 if rows == 2 else 8.6), sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 2] if rows == 2 else [3, 2, 2]},
    )
    ax, dax = axes[0], axes[1]
    if random_cdfs is not None:
        curves(ax, x, random_cdfs, RANDOM, 0.8, 0.5, 0)
    curves(ax, x, tested, TESTED, 0.5, 0.35, 1)
    curves(ax, x, accepted, ACCEPTED, 0.8, 0.55, 2)
    ax.plot(x, observed, color=OBSERVED, linewidth=2.2, zorder=3)
    # The sets track the observed CDF closely, so the top panel cannot separate
    # them; the bottom panel shows each set's CDF minus the observed CDF.
    curves(dax, x, tested - observed, TESTED, 0.5, 0.35, 1)
    curves(dax, x, accepted - observed, ACCEPTED, 0.8, 0.55, 2)
    dax.axhline(0, color=OBSERVED, linewidth=2.0, zorder=3)
    if random_cdfs is not None:
        fax = axes[2]
        curves(fax, x, random_cdfs - observed, RANDOM, 0.8, 0.5, 0)
        curves(fax, x, tested - observed, TESTED, 0.5, 0.35, 1)
        curves(fax, x, accepted - observed, ACCEPTED, 0.8, 0.55, 2)
        fax.axhline(0, color=OBSERVED, linewidth=2.0, zorder=3)
        flim = float(np.abs(random_cdfs - observed).max()) * 1.05
        fax.set_ylim(-flim, flim)
        fax.set_ylabel("CDF minus observed\n(full scale)", color=INK)
        style(fax)
        dax.text(0.99, 0.04, "zoomed; random sets off scale",
                 transform=dax.transAxes, ha="right", va="bottom",
                 fontsize=8, color=INK)
    lim = float(np.abs(tested - observed).max()) * 1.05
    dax.set_ylim(-lim, lim)
    if not linear_age:
        ax.set_xscale("log")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1.01)
    axes[-1].set_xlabel("TE age (generations)", color=INK)
    ax.set_ylabel("cumulative fraction of set", color=INK)
    dax.set_ylabel("CDF minus observed", color=INK)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    style(ax)
    style(dax)
    handles = [
        Line2D([], [], color=TESTED, lw=1.5,
               label=f"tested, all restarts (n={tested.shape[0]})"),
        Line2D([], [], color=ACCEPTED, lw=1.5,
               label=f"accepted, final match (n={accepted.shape[0]})"),
        Line2D([], [], color=OBSERVED, lw=2.2, label="observed TEs"),
    ]
    if random_cdfs is not None:
        handles.append(Line2D([], [], color=RANDOM, lw=1.5,
                              label=f"random SNPs (n={random_cdfs.shape[0]})"))
    legend(ax, handles, "lower right")
    return fig


def ecdf(values):
    v = np.sort(np.asarray(values, dtype=float).ravel())
    return v, np.arange(1, v.size + 1) / v.size


def w1_figure(d, title, random_w1=None, random_label="random SNPs"):
    qc = np.asarray(d["qc_pass"], dtype=bool)
    tested = np.asarray(d["restart_match_to_bootstrap_w1"])
    accepted = np.asarray(d["match_to_bootstrap_w1"])[qc]
    observed = np.asarray(d["bootstrap_to_observed_w1"])
    groups = [
        (tested, TESTED, 1.6, f"tested, all restarts (n={tested.size})"),
        (accepted, ACCEPTED, 2.0, f"accepted, final match (n={accepted.size})"),
        (observed, OBSERVED, 2.0,
         f"bootstrap target vs observed (n={observed.size})"),
    ]
    if random_w1 is not None:
        random_w1 = np.asarray(random_w1, dtype=float).ravel()
        groups.append((random_w1, RANDOM, 2.0,
                       f"{random_label} (n={random_w1.size})"))
    fig, ax = plt.subplots(figsize=(7, 4.6), constrained_layout=True)
    for i, (values, color, width, _) in enumerate(groups):
        x, y = ecdf(values)
        ax.step(x, y, where="post", color=color, linewidth=width, zorder=i + 1)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("Wasserstein-1 distance (generations)", color=INK)
    ax.set_ylabel("cumulative fraction of sets", color=INK)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    style(ax)
    handles = [Line2D([], [], color=c, lw=w, label=lab)
               for _, c, w, lab in groups]
    if random_w1 is None:
        legend(ax, handles, "lower right")
    else:
        # Random sets sit orders of magnitude to the right; keep the legend
        # clear of their curve.
        ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.88, 0.03),
                  frameon=False, fontsize=9, labelcolor=INK)
    return fig


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("match_dir", type=Path)
    p.add_argument("--out-dir", type=Path, default=Path("figures/match_sets"))
    p.add_argument("--name", help="output prefix and title (default: match dir name)")
    p.add_argument("--linear-age", action="store_true",
                   help="linear age axis instead of log (drops no bins)")
    p.add_argument("--age-max", type=float,
                   help="drop age bins above this many generations")
    random = p.add_mutually_exclusive_group()
    random.add_argument("--random-dir", type=Path,
                        help="output of tools/random_control_sets.py")
    random.add_argument("--random-w1", type=Path,
                        help=".npy of random-set W1 values (W1 figure only)")
    p.add_argument("--random-label", help="legend text for the random baseline")
    args = p.parse_args()

    name = args.name or args.match_dir.name
    d = load(args.match_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    random_cdfs = random_w1 = None
    label = args.random_label or "random SNPs"
    if args.random_dir is not None:
        random_cdfs = np.load(args.random_dir / "random_cdfs.npy", mmap_mode="r")
        random_w1 = np.load(args.random_dir / "random_w1.npy")
        bins = np.load(args.random_dir / "age_bins.npy")
        if bins.shape != d["age_bins"].shape or not np.allclose(bins, d["age_bins"]):
            raise ValueError("random baseline age grid differs from the match's")
    elif args.random_w1 is not None:
        random_w1 = np.load(args.random_w1)
    figs = {
        "age_cdf": age_cdf_figure(d, f"{name}: age CDF of matched sets",
                                  args.linear_age, args.age_max,
                                  None if random_cdfs is None
                                  else np.asarray(random_cdfs)),
        "w1": w1_figure(d, f"{name}: W1 reached by matched sets",
                        random_w1, label),
    }
    for suffix, fig in figs.items():
        for ext in ("png", "pdf"):
            path = args.out_dir / f"{name}_{suffix}.{ext}"
            fig.savefig(path, dpi=200)
            print(path)
        plt.close(fig)


if __name__ == "__main__":
    main()
