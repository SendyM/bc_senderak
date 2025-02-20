#!/usr/bin/env python3
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Plot detailed histograms, boxplots, and a scatter plot for segment metrics (degree, centrality) with user-defined ranges."
    )
    parser.add_argument("table_file", help="Path to the input table file (CSV/TSV) with columns: segment, degree, centrality.")
    parser.add_argument("--delimiter", default="\t", help="Delimiter used in the table file (default: tab).")
    parser.add_argument("--bins", type=int, default=50, help="Number of bins for histograms (default: 50).")

    # Degree range arguments
    parser.add_argument("--degree-min", type=float, default=0,
                        help="Minimum degree value to display/filter (default: 0).")
    parser.add_argument("--degree-max", type=float, default=100,
                        help="Maximum degree value to display/filter (default: 100).")

    # Centrality range arguments
    parser.add_argument("--centrality-min", type=float, default=0.0,
                        help="Minimum centrality value to display/filter (default: 0.0).")
    parser.add_argument("--centrality-max", type=float, default=0.05,
                        help="Maximum centrality value to display/filter (default: 0.05).")

    # Output files
    parser.add_argument("--output-distributions", default="distributions_zoomed.png",
                        help="Output filename for the 2×2 subplot figure (default: distributions_zoomed.png).")
    parser.add_argument("--output-scatter", default="degree_vs_centrality_zoomed.png",
                        help="Output filename for the scatter plot (default: degree_vs_centrality_zoomed.png).")

    args = parser.parse_args()

    # Read the table file into a pandas DataFrame
    try:
        df = pd.read_csv(args.table_file, delimiter=args.delimiter)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Ensure required columns exist
    required_cols = {"segment", "degree", "centrality"}
    if not required_cols.issubset(set(df.columns)):
        print(f"Error: The input file must contain columns: {required_cols}")
        return

    # Convert columns to numeric
    df["degree"] = pd.to_numeric(df["degree"], errors="coerce")
    df["centrality"] = pd.to_numeric(df["centrality"], errors="coerce")

    # Drop rows with NaN
    df.dropna(subset=["degree", "centrality"], inplace=True)

    # Filter the DataFrame based on user-specified ranges
    df_filtered = df[
        (df["degree"] >= args.degree_min) &
        (df["degree"] <= args.degree_max) &
        (df["centrality"] >= args.centrality_min) &
        (df["centrality"] <= args.centrality_max)
    ].copy()

    # 1) Create a single figure with four subplots (2×2) for distributions
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Degree and Betweenness Centrality Distributions (Zoomed)", fontsize=16)

    # Top-left: Histogram + KDE for degree
    sns.histplot(data=df_filtered, x="degree", bins=args.bins, kde=True, edgecolor="black", ax=axes[0, 0])
    axes[0, 0].set_title("Degree (Histogram + KDE)")
    axes[0, 0].set_xlabel("Degree")
    axes[0, 0].set_ylabel("Frequency")
    axes[0, 0].grid(True, alpha=0.5)
    axes[0, 0].set_xlim(args.degree_min, args.degree_max)

    # Bottom-left: Boxplot for degree
    sns.boxplot(data=df_filtered, x="degree", ax=axes[1, 0])
    axes[1, 0].set_title("Boxplot for Degree")
    axes[1, 0].grid(True, alpha=0.5)
    axes[1, 0].set_xlim(args.degree_min, args.degree_max)

    # Top-right: Histogram + KDE for betweenness centrality
    sns.histplot(data=df_filtered, x="centrality", bins=args.bins, kde=True, edgecolor="black", ax=axes[0, 1])
    axes[0, 1].set_title("Betweenness Centrality (Histogram + KDE)")
    axes[0, 1].set_xlabel("Centrality")
    axes[0, 1].set_ylabel("Frequency")
    axes[0, 1].grid(True, alpha=0.5)
    axes[0, 1].set_xlim(args.centrality_min, args.centrality_max)

    # Bottom-right: Boxplot for betweenness centrality
    sns.boxplot(data=df_filtered, x="centrality", ax=axes[1, 1])
    axes[1, 1].set_title("Boxplot for Centrality")
    axes[1, 1].grid(True, alpha=0.5)
    axes[1, 1].set_xlim(args.centrality_min, args.centrality_max)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(args.output_distributions)
    plt.close(fig)
    print(f"Distribution figure saved as '{args.output_distributions}'.")

    # 2) Create a separate figure for scatter plot (degree vs. centrality)
    fig_scatter, ax_scatter = plt.subplots(figsize=(8, 6))
    ax_scatter.scatter(df_filtered["degree"], df_filtered["centrality"], alpha=0.7)
    ax_scatter.set_xlabel("Degree")
    ax_scatter.set_ylabel("Betweenness Centrality")
    ax_scatter.set_title("Scatter Plot: Degree vs. Centrality (Zoomed)")
    ax_scatter.grid(True, alpha=0.5)
    ax_scatter.set_xlim(args.degree_min, args.degree_max)
    ax_scatter.set_ylim(args.centrality_min, args.centrality_max)

    plt.tight_layout()
    plt.savefig(args.output_scatter)
    plt.close(fig_scatter)
    print(f"Scatter plot saved as '{args.output_scatter}'.")

if __name__ == "__main__":
    main()

