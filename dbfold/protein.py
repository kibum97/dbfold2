import os
import numpy as np
from scipy.special import logsumexp
import pandas as pd
import mdtraj as md
import multiprocessing
import natsort
import argparse
import glob
import re
from collections import defaultdict
import joblib
import pickle
from functools import partial
from dbfold.utils import substructures as subs
from dbfold.utils import postprocessing
import seaborn as sns
import matplotlib.pyplot as plt

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
        log_df = postprocessing.logfiles_to_dataframe(self.datadir)
        self.eq_step = eq_step
        self.log_df = log_df[log_df.index.get_level_values(2) >= self.eq_step]

    def create_mbar(self,k_bias,recompute=False, **mbar_kwargs):
        self.k_bias = k_bias
        if os.path.exists(f'{self.datadir}/mbar.pkl') and not recompute:
            self.mbar = pickle.load(open(f'{self.datadir}/mbar.pkl','rb'))
        else:
            self.mbar = postprocessing.initialize_mbar(self.log_df, self.k_bias, self.datadir, **mbar_kwargs)
    
    def create_fes(self,k_bias,recompute=False, **fes_kwargs):
        self.k_bias = k_bias
        if os.path.exists(f'{self.datadir}/fes.pkl') and not recompute:
            self.fes = pickle.load(open(f'{self.datadir}/fes.pkl','rb'))
        else:
            self.fes = postprocessing.initialize_fes(self.log_df, self.k_bias, self.datadir, **fes_kwargs)
    """
    TODO: ADD MCPU feature to track replica and use auto detection
    def detect_equilibration(self):
        unique_trajs = self.log_df.index.droplevel("step").unique()
        for unique_traj in unique_trajs:
            traj = self.log_df.loc[unique_traj]
            g_k[k] = timeseries.statistical_inefficiency(u_kn[k, :], u_kn[k, 0 : N_k[k]])
            print(f"Correlation time for set {k:5d} is {g_k[k]:10.3f}")
            indices = timeseries.subsample_correlated_data(u_kn[k, 0 : N_k[k]])
            print(indices)
    """

    def compute_weights(self,temperature):
        u_n = self.log_df['energy'].values
        log_w_nb = self.mbar._computeUnnormalizedLogWeights(u_n/temperature)
        max_log_w_nb = np.max(log_w_nb)  # Use log to solve underflow
        w_nb = np.exp(log_w_nb - max_log_w_nb)
        w_nb = w_nb / np.sum(w_nb)
        self.weights = w_nb
        self.log_weights = log_w_nb

#    def compute_fes_kde(self,temperature):

    
    def compute_free_energy(self,indices):
        return -logsumexp(self.log_weights[indices])
    
    def compute_prob(self, indices):
        return np.sum(self.weights[indices])
    
    def get_xtc_list(self):
        xtc_list = glob.glob(f'{self.datadir}/*.xtc')
        xtc_list = natsort.natsorted(xtc_list)
        self.xtc_list = xtc_list
    
    #TODO: Move this to substructures.py and modify so that it can be used for compute_features
    def define_substructure(self, substructure_file=None, min_seq_separation=8, contact_sep_thresh=7, min_clustersize=5):
        self.min_seq_separation = min_seq_separation
        self.contact_sep_thresh = contact_sep_thresh
        self.min_clustersize = min_clustersize
        if substructure_file:
            substructures, native_distances, d_cutoff, min_seq_separation = joblib.load(substructure_file)
        else:
            native_distances, substructures = subs.generate_substructures(self.native, self.dist_cutoff,
                                                                    min_seq_separation, contact_sep_thresh, min_clustersize)
        self.substructures = substructures
        self.native_distances = native_distances
    
    def compute_score(self):
        # substructure to pair
        pair_list = []
        for s in range(self.substructures.shape[2]):
            sub = self.substructures[:,:,s]
            temp_list = []
            for i in range(sub.shape[0]):
                for j in range(i,sub.shape[1]):
                    if (sub[i,j] == 1.0) or (sub[j,i] == 1.0):
                        temp_list.append([i,j])
            pair_list.append(temp_list)
        # compute score
        tot_score_list = []
        for xtc in self.xtc_list:
            traj = md.load(xtc,top=self.native)
            traj = traj[traj.time >= self.eq_step]
            traj = traj[traj.time <= self.log_df.index.get_level_values(2).max()]
            traj = traj.atom_slice(traj.top.select('name CA'))
            #for step in remove_list:
            #    print(f'Step {step} removed')
            #    traj = traj[traj.time != step]
            score_list = []
            for pair in pair_list:
                score_list.append(np.mean(md.compute_distances(traj,pair),axis=1))
            score_array = np.stack(score_list).T*10
            tot_score_list.append(score_array)
        score_arr = np.vstack(tot_score_list)
        self.score = score_arr
    
    def convert_score2subs(self,f=1.7,mean_substructure_distances=None):
        if mean_substructure_distances:
            self.subs = subs.load_scores(self.score, self.native_distances, self.substructures, f, mean_substructure_distances=mean_substructure_distances)
        else:
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

    def compute_fes(self,features,temperature,ftype='discrete',nbins=30,n_bootstrap=None):
        self.log_df['feat'] = features
        #self.log_df['mbar'] = self.mbar
        if ftype == 'discrete':
            histogram_parameters = {}
            bin_edges = np.arange(self.log_df['feat'].values.min(),self.log_df['feat'].values.max()+2,1)
            bin_center_i = np.arange(self.log_df['feat'].values.min(),self.log_df['feat'].values.max()+1,1)
            hist_values, bin_edges = np.histogram(self.log_df['feat'].values, bins=bin_edges)
        elif ftype == 'continuous':
            histogram_parameters = {}
            nbins=nbins
            bin_center_i = np.zeros([nbins], np.float64)
            hist_values, bin_edges = np.histogram(self.log_df['feat'].values, bins=nbins)
            for i in range(nbins):
                bin_center_i[i] = 0.5 * (bin_edges[i] + bin_edges[i + 1])
        
        nconditions = self.log_df.index.droplevel("step").nunique()
        bin_edges[-1] += 1e-1 # To ensure np.digitize to work properly
        histogram_parameters["bin_edges"] = bin_edges
        if n_bootstrap is not None:
            self.fes.generate_fes(self.log_df['energy'].values.reshape(nconditions,-1)/temperature,self.log_df['feat'].values,fes_type="histogram",histogram_parameters=histogram_parameters,n_bootstraps=n_bootstrap)
        else:
            self.fes.generate_fes(self.log_df['energy'].values.reshape(nconditions,-1)/temperature,self.log_df['feat'].values,fes_type="histogram",histogram_parameters=histogram_parameters)
        uncertainty_method = "analytical"#"bootstrap"
        bin_center_i = [self.fes.histogram_data["bin_label"][key] for key in self.fes.histogram_data["nonzero_bins"]]
        results = self.fes.get_fes(
            bin_center_i, reference_point="from-lowest", uncertainty_method=uncertainty_method
        )
        #results = self.fes.get_fes(
        #    bin_center_i[hist_values != 0], reference_point="from-lowest", uncertainty_method=uncertainty_method
        #)
        return results, bin_center_i

    def merge_traj_from_indices(self,indices):
        """
        Create a single trajectory from snapshots of given indices
        Parameters
        ----------
        indices : list
            List of indices of snapshots
        Returns
        -------
        saving_traj : mdtraj.Trajectory
            Trajectory of snapshots
        """
        saving_snaps = self.log_df.loc[indices].index.values
        saving_dict = defaultdict(list)
        saving_traj = []
        for snap in saving_snaps:
            temperature, setpoint, step = snap
            saving_dict[f'{temperature:.3f}_{int(setpoint)}'].append(step)
        for key, value in saving_dict.items():
            traj = md.load(f'{self.datadir}/{self.pdbroot}_{key}.xtc',top=self.native)
            traj = traj[np.where(np.isin(traj.time, value))[0]]
            saving_traj.append(traj)
        saving_traj = md.join(saving_traj)
        return saving_traj
    
    def merge_traj_from_indices_nous(self,indices):
        """
        Create a single trajectory from snapshots of given indices
        Parameters
        ----------
        indices : list
            List of indices of snapshots
        Returns
        -------
        saving_traj : mdtraj.Trajectory
            Trajectory of snapshots
        """
        saving_snaps = self.log_df.loc[indices].index.values
        saving_dict = defaultdict(list)
        saving_traj = []
        for snap in saving_snaps:
            temperature, setpoint, step = snap
            saving_dict[f'{round(temperature,3):.3f}'].append(step)
        for key, value in saving_dict.items():
            traj = md.load(f'{self.datadir}/{self.pdbroot}_{key}.xtc',top=self.native)
            traj = traj[np.where(np.isin(traj.time, value))[0]]
            saving_traj.append(traj)
        saving_traj = md.join(saving_traj)
        return saving_traj
    
    def compute_contacts_frequency(self,indices=None):
        contacts_list = []
        for xtc in self.xtc_list:
            traj = md.load(xtc, top=self.native)
            traj = traj[traj.time >= self.eq_step]
            traj = traj[traj.time <= self.log_df.index.get_level_values(2).max()]
            traj = traj.atom_slice(traj.top.select('name CA'))
            results = postprocessing.compute_contacts(traj, mode='contacts',
                                              min_seq_separation=self.min_seq_separation,
                                              dist_cutoff=self.dist_cutoff,
                                              sqaureform=False)
            distances = results['contacts']
            pairs = results['pair']
            contacts_list.append(distances)
        contacts = np.vstack(contacts_list)
        weights = self.weights
        print(contacts.shape, pairs.shape)
        if indices is not None:
            contacts = contacts[indices]
            weights = weights[indices]
        self.contacts = contacts
        weights = weights / weights.sum()
        self.contacts_frequency = np.average(contacts, axis=0, weights=weights)
        self.contacts_pairs = pairs
        return self.contacts_pairs, self.contacts_frequency

    def compute_average_distance(self,indices=None):
        distances_list = []
        for xtc in self.xtc_list:
            traj = md.load(xtc, top=self.native)
            traj = traj[traj.time >= self.eq_step]
            traj = traj[traj.time <= self.log_df.index.get_level_values(2).max()]
            traj = traj.atom_slice(traj.top.select('name CA'))
            distances_list.append(postprocessing.compute_contacts(traj, mode='distances',
                                              min_seq_separation=self.min_seq_separation,
                                              dist_cutoff=self.dist_cutoff,
                                              sqaureform=True))
        distances = np.vstack(distances_list)
        weights = self.weights
        if indices is not None:
            contacts = contacts[indices]
            weights = weights[indices]
        weights = weights / weights.sum()
        self.average_distance = np.average(distances, axis=0, weights=weights)
        return self.average_distance

    def plot_contacts_frequency(self,savepath=None):
        fig, ax = plt.subplots()
        contacts_frequency_map = md.geometry.squareform(self.contacts_frequency.reshape(1,-1),self.contacts_pairs)[0]
        sns.heatmap(contacts_frequency_map,cmap='Reds', square=True, ax=ax, vmin=0, vmax=1,
                    linewidths=0, linecolor='lightgrey')
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
            spine.set_color('black')
        ticks = [i + 0.5 for i in range(0, len(contacts_frequency_map), 10)]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([str(int(i+0.5)) for i in ticks])
        ax.set_yticklabels([str(int(i+0.5)) for i in ticks])
        plt.tight_layout()
        if savepath:
            plt.savefig(savepath)
            plt.close()
        else:
            plt.show()
            plt.close()       
    
    def plot_average_distance(self,savepath=None):
        fig, ax = plt.subplots()
        sns.heatmap(self.average_distance,cmap='YlGnBu',square=True,ax=ax, linewidths=0.5, linecolor='lightgrey')
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
            spine.set_color('black')
        ticks = [i + 0.5 for i in range(0, len(self.average_distance), 10)]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([str(int(i+0.5)) for i in ticks])
        ax.set_yticklabels([str(int(i+0.5)) for i in ticks])
        plt.tight_layout()
        if savepath:
            plt.savefig(savepath)
            plt.close()
        else:
            plt.show()
            plt.close()