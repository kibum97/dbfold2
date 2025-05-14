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
        prot.get_xtc_list()
        prot.dist_cutoff = args.dist_cutoff

        # Create the MBAR and FES objects
        prot.create_fes(k_bias=0.02, recompute=True,solver_protocol="robust")
        prot.mbar = prot.fes.mbar
        if args.cache:
            os.makedirs(f'{args.savedir}/cache', exist_ok=True)
            with open(f'{args.savedir}/cache/{pdbroot}_fes.pkl', 'wb') as pickle_file:
                pickle.dump(prot.fes, pickle_file)

        # Load predefined substructure definitions
        substructure_dict = pickle.load(open(f'{args.substrcutre_dict}', 'rb'))
        native_distances, substructures = subs.load_substructures(native_file, substructure_dict, min_seq_separation=8)
        prot.min_seq_separation = 8
        prot.native_distances = native_distances
        prot.substructures = substructures
        subs.visualize_substructures(
            prot.native_distances, prot.substructures, prot.dist_cutoff, prot.min_seq_separation,
            onlylower=True, savepath=f'{args.savedir}/{pdbroot}_substructures.png'
        )

        # Compute the substructure labels and convert them to decimal
        prot.compute_score()
        prot.convert_score2subs(f=args.ratio)
        dec_subs = [subs.bin2dec(sub) for sub in prot.subs]

        for T in tqdm.tqdm(args.temperatures):
            # Compute the free energy landscape
            fes_result = prot.compute_fes(dec_subs, T, ftype='discrete')
            subs.plot_substructure_fes(fes_result, f"{args.savedir}/PMF_{pdbroot}_plot_{T}.png")
            for ymax in [5, 10, 15, 20, 25, 30]:
                subs.plot_substructure_fes(fes_result, f"{args.savedir}/PMF_{pdbroot}_plot_{T}_ymax{ymax}.png", ymax=ymax)
            with open(f'{args.savedir}/cache/PMF_{pdbroot}_plot_{T}.pkl', 'wb') as pickle_file:
                pickle.dump(fes_result, pickle_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description="""Generate free energy landscapes using DBFOLD for given actin variants.

    Expected directory structure:
    <basedir>/
    ├── <pdbroot>/
    │   ├── in_pdbs/
    │   │   └── <pdbroot>.pdb (native structure)
    │   └── MCPU_run/ (trajectory files and log files)
    ├── <substructure_definition>.pkl (substructure definitions)
    ├── <homedir>/
    │   └── <savedir>/ (output saved here)
    |       └── <cache>/ (results saved here; fes and pmf files)

    Example usage:
    python script.py --pdbroots ACTA ACTB --eq_step 600000000 --dist_cutoff 7.5
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
        '--substrcutre_dict', type=str, required=True,
        help='Path to the substructure definitions file (e.g., substructures_definition.pkl)'
    )
    parser.add_argument(
        '--dist_cutoff', type=float, default=7.5,
        help='Distance cutoff (in Å) used to define contacts for substructure analysis'
    )
    parser.add_argument(
        '--ratio', type=float, default=1.2,
        help='Ratio for native distances and sampled distances to define substructures'
    )
    parser.add_argument(
        '--temperatures', nargs='+', type=float, default=[0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000],
        help='List of simulation temperatures for free energy landscape calculations'
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