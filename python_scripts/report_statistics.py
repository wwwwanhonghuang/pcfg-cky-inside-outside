def analyze_dataframes(dataframes, key_pattern=None):
    """
    Combines and analyzes multiple DataFrames.

    Args:
        dataframes (list of pd.DataFrame): List of DataFrames to analyze.
        key_pattern (str, optional): A regex pattern to filter keys.

    Returns:
        pd.DataFrame: A combined DataFrame with descriptive statistics.
    """
    # Combine DataFrames
    combined_df = pd.concat(dataframes, axis=0)

    # Filter keys if a pattern is provided
    if key_pattern:
        combined_df = combined_df[combined_df['Key'].str.contains(key_pattern, regex=True)]

    # Perform statistical analysis
    stats = combined_df.groupby('Key')['Value'].describe()
    return stats


def plot_statistics(stats, title="Statistical Analysis", output_file=None):
    """
    Plots descriptive statistics.

    Args:
        stats (pd.DataFrame): Descriptive statistics DataFrame (result of analyze_dataframes).
        title (str): Title of the plot.
        output_file (str, optional): File path to save the plot. If None, displays the plot.
    """
    plt.figure(figsize=(10, 6))

    # Plot the mean with error bars (standard deviation)
    stats['mean'].plot(kind='bar', yerr=stats['std'], capsize=4, color='skyblue', alpha=0.8)

    plt.title(title)
    plt.ylabel("Mean Value with Std Dev")
    plt.xlabel("Key")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

def plot_comparison(stats_groups, group_labels, title="Statistical Comparison", output_file=None):
    """
    Plots statistical comparison of multiple groups.

    Args:
        stats_groups (list of pd.DataFrame): List of descriptive statistics DataFrames for each group.
        group_labels (list of str): Labels for each group.
        title (str): Title of the plot.
        output_file (str, optional): File path to save the plot. If None, displays the plot.
    """
    plt.figure(figsize=(12, 8))

    # Plot the mean values for each group
    for i, stats in enumerate(stats_groups):
        stats['mean'].plot(kind='line', label=group_labels[i], marker='o')

    plt.title(title)
    plt.ylabel("Mean Value")
    plt.xlabel("Key")
    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    else:
        plt.show()