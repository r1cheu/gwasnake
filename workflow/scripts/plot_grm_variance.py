import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_variance_components(table: Path, label: str) -> dict[str, float | str]:
    variance_table = pd.read_csv(table, sep="\t")
    ratios = variance_table.set_index("component")["ratio"]
    ratio_ses = pd.to_numeric(
        variance_table.set_index("component")["ratio_se"], errors="coerce"
    ).fillna(0.0)
    return {
        "label": label,
        "qtl_add": float(ratios.get("qtl_add", 0.0)),
        "qtl_dom": float(ratios.get("qtl_dom", 0.0)),
        "background_add": float(ratios.get("background_add", 0.0)),
        "background_dom": float(ratios.get("background_dom", 0.0)),
        "qtl_add_se": float(ratio_ses.get("qtl_add", 0.0)),
        "qtl_dom_se": float(ratio_ses.get("qtl_dom", 0.0)),
        "background_add_se": float(ratio_ses.get("background_add", 0.0)),
        "background_dom_se": float(ratio_ses.get("background_dom", 0.0)),
    }


def plot_variance_components(data: pd.DataFrame, output: Path) -> None:
    if output.suffix == "":
        output = output.with_suffix(".svg")
    output.parent.mkdir(parents=True, exist_ok=True)
    labels = data["label"].tolist()
    y = np.arange(len(labels))
    height = 0.34

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18.0, max(7.0, len(labels) * 0.7)),
        sharex=True,
    )
    colors = {
        "qtl": "#0FA7A8",
        "background": "#E8B900",
        "total": "#E879F9",
        "se": "#FF6B4A",
    }

    qtl_add = data["qtl_add"].to_numpy()
    qtl_dom = data["qtl_dom"].to_numpy()
    background_add = data["background_add"].to_numpy()
    background_dom = data["background_dom"].to_numpy()
    qtl_add_se = data["qtl_add_se"].to_numpy()
    qtl_dom_se = data["qtl_dom_se"].to_numpy()
    background_add_se = data["background_add_se"].to_numpy()
    background_dom_se = data["background_dom_se"].to_numpy()
    panels = [
        ("Total", qtl_add + qtl_dom, background_add + background_dom, None, None),
        ("Additive", qtl_add, background_add, qtl_add_se, background_add_se),
        ("Dominance", qtl_dom, background_dom, qtl_dom_se, background_dom_se),
    ]

    all_totals = [
        qtl_values + background_values
        for _, qtl_values, background_values, _, _ in panels
    ]
    all_error_upper_limits = [
        qtl_values + qtl_se
        for _, qtl_values, _, qtl_se, _ in panels
        if qtl_se is not None
    ] + [
        background_values + background_se
        for _, _, background_values, _, background_se in panels
        if background_se is not None
    ]
    x_max = (
        max(
            [total.max() for total in all_totals]
            + [error_limit.max() for error_limit in all_error_upper_limits]
        )
        + 0.1
    )

    axes[0].set_ylabel("Trait")
    for panel_index, (
        ax,
        (title, qtl_values, background_values, qtl_se, background_se),
    ) in enumerate(zip(axes, panels)):
        total_values = qtl_values + background_values
        ax.barh(
            y - height / 2,
            qtl_values,
            height,
            color=colors["qtl"],
            label="QTL",
            xerr=qtl_se,
            error_kw={
                "ecolor": colors["se"],
                "elinewidth": 0.3,
                "capsize": 3,
                "capthick": 0.3,
                "clip_on": False,
            },
        )
        ax.barh(
            y + height / 2,
            background_values,
            height,
            color=colors["background"],
            label="Background",
            xerr=background_se,
            error_kw={
                "ecolor": colors["se"],
                "elinewidth": 0.3,
                "capsize": 3,
                "capthick": 0.3,
                "clip_on": False,
            },
        )
        for index, total_value in enumerate(total_values):
            ax.vlines(
                total_value,
                y[index] - height,
                y[index] + height,
                color=colors["total"],
                alpha=0.55,
                linestyle=(0, (2, 2)),
                linewidth=0.7,
            )
        ax.set_title(title)
        ax.set_yticks(y)
        ax.set_yticklabels(labels if panel_index == 0 else [])
        ax.set_ylim(-0.6, len(labels) - 0.4)
        ax.invert_yaxis()
        ax.set_xlim(0, x_max)
        ax.set_xlabel("Variance proportion")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.2)
        ax.spines["bottom"].set_linewidth(0.2)
        ax.tick_params(width=0.2, length=4)

    axes[0].plot(
        [],
        [],
        color=colors["total"],
        alpha=0.8,
        linestyle=(0, (2, 2)),
        linewidth=0.70,
        label="Total variance",
    )
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        frameon=False,
        ncol=3,
        loc="upper center",
    )

    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(output, dpi=300)


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot QTL and background variance proportions from variance.tsv tables."
    )
    parser.add_argument(
        "tables",
        nargs="+",
        type=Path,
        help="One or more results/{run_id}/{phenotype}/variance/variance.tsv tables.",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        type=Path,
        help="Output figure path.",
    )
    return parser.parse_args()


def main() -> None:
    args = get_args()
    # path is .../{phenotype}/variance/variance.tsv -> trait label is the phenotype dir
    rows = [
        read_variance_components(table, table.parent.parent.name)
        for table in args.tables
    ]
    plot_variance_components(pd.DataFrame(rows), args.out)


if __name__ == "__main__":
    main()
