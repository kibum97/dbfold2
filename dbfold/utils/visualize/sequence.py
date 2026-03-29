import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_sequence_subdomains(subdomains, protein_length, ax=None, height=0.2, colors=None):
    """
    Draws a schematic representation of protein subdomains.

    Parameters
    ----------
    subdomains : dict
        Dictionary where keys are subdomain labels and values are tuples of (start, end) indices.
    protein_length : int
        Total length of the protein sequence.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The Axes object with the drawn subdomains.
    """
    
    # Create a figure and axis for the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 2))
    else:
        fig = plt.gcf()
    #ax.add_patch(patches.Rectangle((0, -0.2), protein_length, 0.4, color='lightgray'))  # full protein
    ax.plot([0, protein_length], [0, 0], color='black', linewidth=1.5, zorder=0)

    for label, (start, end) in subdomains.items():
        width = end - start
        color = colors[label] if colors is not None else 'steelblue'
        rounded_box = patches.FancyBboxPatch(
            (start, -height / 2), width, height,
            boxstyle=f"round, pad=0.3, rounding_size={height / 2:.2f}",
            facecolor=color,
            edgecolor="black",
            linewidth=1
        )
        ax.add_patch(rounded_box)
        ax.text((start + end) / 2, 0.0, label, ha='center', va='center', fontsize=10)

    ax.set_xlim(-1, protein_length+1)
    ax.set_aspect('equal')
    ax.axis('off')

    return ax