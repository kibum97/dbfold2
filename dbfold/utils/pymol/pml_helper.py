import numpy as np

def alter_b_factor(new_factor_dict):
    """
    Generate Script that alters the B-factor of atoms

    Parameters
    ----------
    new_factor_dict : dict
        Dictionary where keys are atom identifiers (e.g., 'resi {i} and name CA')
        and values are the new B-factors.
    Returns
    -------
    str
        A string containing the PyMOL script to alter the B-factors.
    """
    script_lines = []
    for atom, new_factor in new_factor_dict.items():
        script_lines.append(f"alter {atom}, b={new_factor}")
    return "\n".join(script_lines) + "\n"

def color_by_b_factor(b_factors):
    """
    Generate Script that colors atoms based on their B-factor

    Parameters
    ----------
    b_factors : dict
        Dictionary where keys are atom identifiers (e.g., 'resi {i} and name CA')
        and values are the B-factors.
    Returns
    -------
    str
        A string containing the PyMOL script to color atoms by B-factor.
    """
    return f"spectrum b, rainbow, minimum={min(b_factors.values())}, maximum={max(b_factors.values())}\n"

def sausage_by_b_factor(b_factors, radius=1.5):
    """
    Generate Script that creates putty representations based on B-factor

    Parameters
    ----------
    b_factors : dict
        Dictionary where keys are atom identifiers (e.g., 'resi {i} and name CA')
        and values are the B-factors.
    radius : float, optional
        Radius of the putty representation. Default is 0.5.
    Returns
    -------
    str
        A string containing the PyMOL script to create putty representations by B-factor.
    """
    script_lines = [
        'show cartoon',
        'cartoon putty',
        'set cartoon_putty_transform, 7',
        f'set cartoon_putty_radius, {radius}',
        'set cartoon_putty_scale_max, -1',
        'set cartoon_putty_quality, 20',
        'unset cartoon_smooth_loops',
        'unset cartoon_flat_sheets',
    ]   
    return "\n".join(script_lines) + "\n"