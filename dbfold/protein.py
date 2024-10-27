import os
import numpy as np
import pandas as pd
import mdtraj as md
import multiprocessing
import natsort
import argparse
import glob
import re
import collections
import joblib
import pickle
from functools import partial
from dbfold.utils import substructures as subs
from dbfold.utils import utils

class Protein:
    def __init__(self, native, pdbroot, homedir, savedir=None, datadir=None):
        self.native = native
        self.pdbroot = pdbroot
        self.homedir = homedir
        if savedir:
            self.savedir = savedir
        else:
            self.savedir = f'{homedir}/results/'
        if datadir:
            self.datadir = datadir
        else:
            self.datadir = f'{homedir}/{pdbroot}/MCPU_run/'

    def create_log_dataframe(self,eq_step):
        log_df = utils.logfiles_to_dataframe(self.datadir)
        self.log_df = log_df[log_df.index.get_level_values(2) >= eq_step]

    def create_mbar(self,k_bias):
        self.k_bias = k_bias
        if os.path.exists(f'{self.datadir}/mbar.pkl'):
            self.mbar = pickle.load(open(f'{self.datadir}/mbar.pkl','rb'))
        else:
            self.mbar = utils.initialize_mbar(self.log_df, self.k_bias, self.datadir)
    
    def create_fes(self,k_bias):
        self.k_bias = k_bias
        if os.path.exists(f'{self.datadir}/fes.pkl'):
            self.fes = pickle.load(open(f'{self.datadir}/fes.pkl','rb'))
        else:
            self.fes = utils.initialize_fes(self.log_df, self.k_bias, self.datadir)
        
    def get_xtc_list(self):
        xtc_list = glob.glob(f'{self.datadir}/*.xtc')
        xtc_list = natsort.natsorted(xtc_list)
        self.xtc_list = xtc_list
    
    def define_substructure(self, substructure_file=None):
        # Need to correct
        if substructure_file:
            substructures, native_distances, d_cutoff, min_seq_separation = joblib.load(substructure_file)
        else:
            min_seq_separation = 8
            contact_sep_thresh = 7
            min_clustersize = 5
            native_distances, substructures = subs.generate_substructures(self.native, self.dist_cutoff,
                                                                    min_seq_separation, contact_sep_thresh, min_clustersize)
        
        self.substructures = substructures
        self.native_distances = native_distances
    
    def compute_score(self):
        # substructure to pair
        pair_list = []
        for s in range(substructures.shape[2]):
            sub = substructures[:,:,s]
            temp_list = []
            for i in range(sub.shape[0]):
                for j in range(i,sub.shape[1]):
                    if sub[j,i] == 1.0:
                        temp_list.append([i,j])
            pair_list.append(temp_list)
        # compute score
        tot_score_list = []
        for xtc in prot.xtc_list:
            traj = md.load(xtc,top=native)
            traj = traj.atom_slice(traj.top.select('name CA'))
            for step in remove_list:
                print(f'Step {step} removed')
                traj = traj[traj.time != step]
            score_list = []
            for pair in pair_list:
                score_list.append(np.mean(md.compute_distances(traj,pair),axis=1))
            score_array = np.stack(score_list).T*10
            tot_score_list.append(score_array)
        score_arr = np.vstack(tot_score_list)
        self.score = score_arr
    
    def convert_score2subs(self,f=1.7):
        self.subs = subs.load_scores(self.score, self.native_distances, self.substructures, f)
        return self.subs

    def compute_features(self,f,save=None):
        """
        Compute features from trajectory
        Parameters
        ----------
        f : function
            Function to compute feature from trajectory; Need to be partial so that only input needed is MDtraj trajectory
        save : str, optional
        """
        # compute feature
        feature_list = []
        for xtc in self.xtc_list:
            traj = md.load(xtc,top=self.native)
            traj = traj.atom_slice(traj.top.select('name CA'))
            for step in self.remove_list:
                print(f'Step {step} removed')
                traj = traj[traj.time != step]
            feature_list = []
            feature_list.append(f(traj))
        feature_array = np.stack(feature_list)
        if save:
            np.save(feature_array,save)
        return feature_array

    def compute_fes(self,features,ftype='discrete'):
        self.log_df['feat'] = features
        self.log_df['mbar'] = mbar
        if ftype == 'discrete':
            print('d')
        elif ftype == 'continuous':
            print('c')
        return 0

    