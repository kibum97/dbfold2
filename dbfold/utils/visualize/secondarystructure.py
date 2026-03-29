from networkx import center
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

def _extract_secondary_segments(dssp_labels):
    """
    Extracts start and end indices of each secondary structure segment from the DSSP labels.
    
    Parameters
    ----------
    dssp_labels : list of str
        List of DSSP labels representing secondary structure elements.
    Returns
    -------
    segments : list of tuples
        List of tuples where each tuple contains the start and end indices of a secondary structure segment.
        Each tuple is of the form (ss_label, start_index, end_index).
    """
    segments = []
    current_segment = None

    for i, label in enumerate(dssp_labels):
        if label in ['H', 'E', 'C']:  # Helix, Strand, Coil
            if current_segment is None or current_segment[0] != label:
                if current_segment is not None:
                    segments.append(tuple(current_segment))
                current_segment = [label, i, i]
            elif current_segment[0] == label:
                # Extend the current segment
                current_segment[2] = i
        else:
            raise ValueError(f"Unexpected DSSP label: {label}. Please use the simplified DSSP codes.")

    if current_segment is not None:
        segments.append(tuple(current_segment))

    return segments

def draw_secondary_structure(dssp_labels, ax=None, skip_strand=3, y_pos=0., width=0.02):
    """
    Patched version of ss.draw_secondary_structure with an overlay option.
    
    Parameters
    ----------
    dssp_labels : list of str
        List of DSSP labels representing secondary structure elements.
    ax : matplotlib.axes.Axes, optional
        Matplotlib Axes object to draw on. If None, a new figure and axes will be created.
    skip_strand : int, optional
        Length of strand to skip when drawing. Default is 3.
    y_pos : float
        Vertical position to center the structure cartoon
    width : float
        Vertical thickness of the annotation
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        The Axes object with the drawn secondary structure.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 2))
    
    segments = _extract_secondary_segments(dssp_labels)
    
    for ss_label, start, end in segments:
        ax.plot([start, end + 1], [y_pos, y_pos], color='black', lw=1, zorder=10)
        if ss_label == 'H':
            x = np.linspace(start, end + 1, 100)
            y = np.sin(2 * np.pi * (x - start)) * width / 2 + y_pos
            ax.plot(x, y, color='black', lw=1.5, zorder=10)
        elif ss_label == 'E':
            if end - start + 1 >= skip_strand:
                ax.add_patch(
                    patches.FancyArrow(
                        start, y_pos, end - start + 1, 0, width=width,
                        length_includes_head=True, head_width=width*2,
                        head_length=0.6, color='black', zorder=10
                    )
                )

    if ax is None:
        ax.set_xlim(0, len(dssp_labels))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel('Residue Index')
        ax.set_title('Secondary Structure Representation')

    return ax