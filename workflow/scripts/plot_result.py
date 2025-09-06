import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def manhattan_plot(
    data,
    chrom_col="CHR",
    pos_col="POS",
    p_col="P",
    title=None,
    threshold=None,
    sort=False,
    log=False,
    colors=["#B39CD0", "#FFC75F"],
    figsize=(15, 7),
    y_break=None,
    output_file="manhattan_plot.png",
):
    data = data.copy()

    data[chrom_col] = data[chrom_col].astype("category")
    if sort:
        data = data.sort_values(by=[chrom_col, pos_col])
    if log:
        data[p_col] = -np.log10(data[p_col])

    data_grouped = data.groupby(chrom_col, observed=True)

    gap = int(len(data) * 0.01)
    data["ind_gapped"] = 0
    x_labels = []
    x_labels_pos = []
    last_pos = 0
    last_pos += gap

    try:
        chrom_names_sorted = sorted(
            data_grouped.groups.keys(), key=lambda x: int(str(x).replace("chr", ""))
        )
    except ValueError:
        chrom_names_sorted = sorted(data_grouped.groups.keys())

    for name in chrom_names_sorted:
        group = data_grouped.get_group(name)
        data.loc[group.index, "ind_gapped"] = np.arange(last_pos, last_pos + len(group))
        x_labels.append(name)
        x_labels_pos.append(last_pos + len(group) / 2)
        last_pos += len(group) + gap

    axes = []
    if y_break and isinstance(y_break, (list, tuple)) and len(y_break) == 2:
        fig, (ax2, ax1) = plt.subplots(
            2, 1, sharex=True, figsize=figsize, gridspec_kw={"height_ratios": [1, 3]}
        )
        fig.subplots_adjust(hspace=0.1)
        axes = [ax1, ax2]
        ax1.set_ylim(0, y_break[0])
        ax2.set_ylim(
            y_break[1], data[p_col].replace([np.inf, -np.inf], np.nan).max() * 1.05
        )
        ax2.tick_params(axis="x", which="both", bottom=False)

        d = 0.002
        kwargs = dict(color="k", clip_on=False)
        ax2.plot((-d, +d), (-d, +d), transform=ax2.transAxes, **kwargs)
        ax1.plot((-d, +d), (1 - d, 1 + d), transform=ax1.transAxes, **kwargs)

        fig.text(
            0.08, 0.5, r"$-\log_{10}(P)$", va="center", rotation="vertical", fontsize=15
        )

        tick_label_size = 12
        ax1.tick_params(axis="y", labelsize=tick_label_size)
        ax1.tick_params(axis="x", labelsize=tick_label_size)
        ax2.tick_params(axis="y", labelsize=tick_label_size)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        axes = [ax]
        ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=15)

    for ax in axes:
        for num, name in enumerate(chrom_names_sorted):
            group = data_grouped.get_group(name)
            group.plot(
                kind="scatter",
                x="ind_gapped",
                y=p_col,
                color=colors[num % len(colors)],
                ax=ax,
                s=10,
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.xaxis.set_ticks_position("none")
        ax.set_xlim([0, last_pos])
        ax.set_ylabel("")
        ax.set_xlabel("")

        if threshold is not None:
            ax.axhline(y=threshold, color="grey", linestyle="--", linewidth=1)

    axes[0].set_xticks(x_labels_pos)
    axes[0].set_xticklabels(x_labels)

    if title is not None:
        axes[-1].set_title(title, fontsize=16)

    if output_file is not None:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")


data = pd.read_csv(snakemake.input.summary, sep=" ")

manhattan_plot(
    data,
    chrom_col="CHROM",
    pos_col="GENPOS",
    threshold=6,
    p_col="LOG10P",
    title=snakemake.params.title,
    output_file=snakemake.output.png,
    figsize=(18, 6),
)
