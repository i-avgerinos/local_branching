from __future__ import division
import preprocessing # Loading instances from json files
import cp_model
import export # Exporting results to csv files
import numpy as np
import pandas as pd
from pyomo.environ import * # Import 'Pyomo'
import json
import logging
import time
import math

logging.getLogger('pyomo.core').setLevel(logging.ERROR)

if __name__ == "__main__":
	timelimit = 3600  # Timelimit for termination of the algorithm
	big_M = 100000  # A big numeric value

	file1 = open('Instances.json') # Import instances from 'file1' - appended in the GitHub directory
	datasets = json.load(file1)

	output_file = "CP_Results.csv"
	export.create_csv(output_file) # Create a csv file to export results

	for ID in datasets.keys():
		nJobs, nMachines, alpha, tau, rho, processing_times, setup_times, deadlines, setup_minus = preprocessing.input_data(datasets[ID], big_M) # Load a new instance 'ID'
		nSlots = nJobs # nSlots : Number of slots - equal to the number of jobs
		for R in range(2):
			resources = int(((R+2)/5)*nMachines) # Number of resources : 2/5 x |M| for R = 0, 3/5 x |M| for R = 1
			print(f"------------------------------------------------------------------")
			print(f"Instance {ID}")
			print(f"------------------------------------------------------------------")
			#algorithm:
			start_time = time.time()
			CP = cp_model.cp(nJobs, nSlots, nMachines, processing_times, big_M, setup_times, resources, deadlines)
			export.add_input(ID, nJobs, nMachines, resources, alpha, CP, int(time.time() - start_time), output_file)