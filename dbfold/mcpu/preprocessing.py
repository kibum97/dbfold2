import mdtraj as md
import numpy as np
import pandas as pd
import warnings

def filter_atoms(traj, atom_indices, save_path=None):
    """
    Filter atoms from a trajectory based on the provided indices.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        The trajectory to filter.
    atom_indices : list
        List of atom indices to keep.

    Returns
    -------
    mdtraj.Trajectory
        The filtered trajectory.
    """
    processed_traj = traj.atom_slice(atom_indices)
    if save_path:
        processed_traj.save(save_path)
    return processed_traj

def load_atom_type(filename):
    """
    Load atom type from a file.

    Parameters
    ----------
    filename : str
        Path to the file containing atom type information.

    Returns
    -------
    dict
        Dictionary mapping (residue, atom) to atom type.
    """
    atom_type_df = pd.read_csv(filename)
    atom_type_dict = {
        (residue, atom): type
        for atom, residue, type in zip(atom_type_df["atom"], atom_type_df["residue"], atom_type_df["type"])
    }
    return atom_type_dict

def get_atom_indices(traj, atom_type_dict):
    """
    Get atom indices based on the atom type DataFrame.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        The trajectory to filter.
    atom_type_df : pandas.DataFrame
        DataFrame containing atom type information.

    Returns
    -------
    list
        List of atom indices to keep.
    """
    atom_indices = []
    
    table, bonds = traj.topology.to_dataframe()
    for i, row in table.iterrows():
        key = (row["resName"], row["name"])
        if key in atom_type_dict:
            atom_indices.append(i)
        elif ('XXX', row["name"]) in atom_type_dict:
            atom_indices.append(i)
        else:            
            warnings.warn(f"Atom {key} not found in atom type dictionary.", UserWarning)
    return atom_indices