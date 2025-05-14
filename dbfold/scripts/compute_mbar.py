import argparse
import os
import pickle
import tqdm
import pandas as pd
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns

import dbfold
import dbfold.utils
from dbfold.utils import substructures as subs
from dbfold import protein

def main(args):
    os.makedirs(args.savedir, exist_ok=True)

    for pdbroot in args.pdbroots:
        # Define the native file and data directory
        native_file = f'{args.basedir}/{pdbroot}/in_pdbs/{pdbroot}.pdb'
        datadir = f'{args.basedir}/{pdbroot}/MCPU_run/'

        # Initialize the Protein object
        prot = protein.Protein(native_file, pdbroot, args.homedir, args.savedir, datadir)
        if os.path.exists(f'{datadir}/logfiles_as_dataframe.pkl'):
            prot.log_df = pd.read_pickle(f'{datadir}/logfiles_as_dataframe.pkl')
            prot.eq_step = args.eq_step
        else:
            prot.create_log_dataframe(eq_step=args.eq_step)
        prot.log_df = prot.log_df[prot.log_df.index.get_level_values(2) >= args.eq_step]
        if args.max_step is not None:
            prot.log_df = prot.log_df[prot.log_df.index.get_level_values(2) < args.max_step]
        prot.get_xtc_list()
        prot.dist_cutoff = args.dist_cutoff

        # Create the MBAR and FES objects
        prot.create_fes(k_bias=0.02, recompute=True,solver_protocol="robust")
        prot.mbar = prot.fes.mbar
        if args.cache:
            os.makedirs(f'{args.savedir}/cache', exist_ok=True)
            with open(f'{args.savedir}/cache/{pdbroot}_fes.pkl', 'wb') as pickle_file:
                pickle.dump(prot.fes, pickle_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description="""Generate free energy landscapes using DBFOLD for given actin variants.

    Expected directory structure:
    <basedir>/
    ├── <pdbroot>/
    │   ├── in_pdbs/
    │   │   └── <pdbroot>.pdb (native structure)
    │   └── MCPU_run/ (trajectory files and log files)
    ├── <homedir>/
    │   └── <savedir>/ (output saved here)
    |       └── <cache>/ (results saved here)

    Example usage:
    python compute_mbar.py --pdbroots ACTA ACTB --eq_step 600000000 --cache
    """
    )
    parser.add_argument(
        '--pdbroots', nargs='+', required=True,
        help='List of actin variants to process (e.g., ACTB ACTA_cryo)'
    )
    parser.add_argument(
        '--eq_step', type=int, required=True,
        help='Minimum equilibration step for filtering log dataframe'
    )
    parser.add_argument(
        '--max_step', type=int, default=None,
        help='Maximum equilibration step for filtering log dataframe'
    )
    parser.add_argument(
        '--cache', action='store_true',
        help='Cache FES results to avoid recomputation'
    )
    cwd = os.getcwd()
    parser.add_argument(
        '--basedir', type=str, default=cwd,
        help='Root directory where per-pdbroot folders with PDB and MCPU_run are located (default: current working directory)'
    )
    parser.add_argument(
        '--homedir', type=str, default=os.path.join(cwd, 'Analysis'),
        help='Directory for analysis scripts, figures, and shared substructure definitions (default: ./Analysis)'
    )
    parser.add_argument(
        '--savedir', type=str, default=os.path.join(cwd, 'Results'),
        help='Path where output figures and results will be saved (default: ./Results)'
    )

    args = parser.parse_args()
    args.substructre_dict = os.path.abspath(args.substrcutre_dict)
    args.basedir = os.path.abspath(args.basedir)
    args.homedir = os.path.abspath(args.homedir)
    args.savedir = os.path.abspath(args.savedir)

    print("Using basedir:", args.basedir)
    print("Using homedir:", args.homedir)
    print("Using savedir:", args.savedir)

    main(args)