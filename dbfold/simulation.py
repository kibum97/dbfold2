import os
import sys
sys.path.append("/n/home01/kibumpark/pkg/dbfold2")
import time
import glob
import natsort
from dbfold.utils import *
import yaml

class Simulation:
    def __init__(self, yaml_file):
        self.config = self.load_yaml(yaml_file)
        
        # General simulation variables
        self.pdbroot = self.config['simulation'].get('pdbroot')
        self.simdir = self.config['simulation'].get('directory')
        self.mc_step = self.config['simulation'].get('mc_step')
        self.log_interval = self.config['simulation'].get('log_interval')
        self.use_cluster = self.config['simulation'].get('use_cluster')
        
        # Check if umbrella biasing is enabled and load related variables
        self.umbrella_biasing_enabled = self.config['simulation']['umbrella'].get('enabled', False)
        if self.umbrella_biasing_enabled:
            self.ca_cutoff = self.config['simulation']['umbrella'].get('ca_cutoff')
            self.num_replicas = self.config['simulation']['umbrella'].get('num_replicas')
        else:
            self.num_replicas = None
            self.temperature_min = None
            self.temperature_max = None
        
        # Check if temperature replica exchange is enabled and load related variables
        self.temperature_replica_enabled = self.config['simulation']['temperature_replica'].get('enabled', False)
        if self.temperature_replica_enabled:
            self.ca_cutoff = self.config['simulation']['temperature_replica'].get('ca_cutoff')
            self.num_replicas = self.config['simulation']['temperature_replica'].get('num_replicas')
        else:
            self.num_replicas = None
            self.temperature_min = None
            self.temperature_max = None
        self.temperature_min = self.config['simulation']['replica_exchange']['temperature_range'].get('min')
        self.temperature_max = self.config['simulation']['replica_exchange']['temperature_range'].get('max')

    def load_yaml(self, yaml_file):
        with open(yaml_file, 'r') as file:
            return yaml.safe_load(file)

    def run(self):
        print(f"Running simulation with var1: {self.var1}, var2: {self.var2}")
        
        if self.replica_exchange_enabled:
            print(f"Replica exchange is enabled with {self.num_replicas} replicas.")
            print(f"Temperature range: {self.temperature_min} to {self.temperature_max}")
        else:
            print("Replica exchange is not enabled.")
    
    def create_config_file(self,run_dir):

    def create_slurm_script(self,run_dir):

    def run(self,logfile=None):

# Usage
simulation = Simulation('simulation.yaml')
simulation.run()