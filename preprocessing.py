import math
import json
import numpy as np

# Load the new instance
def input_data(data, big_M):
    nJobs = int(data['Number_of_jobs'])  # Number of jobs |J|
    nMachines = int(data['Number_of_machines'])  # Number of machines |M|
    alpha = int(data['Alpha'])  # 'alpha' value
    tau = float(data['Tau'])  # 'tau' value
    rho = float(data['Rho'])  # 'rho' value
    processing_times = np.array(list(data['ProcessTimes'].values())).T  # Processing times of jobs to machines
    setup_times = np.array(list(data['SetupTimes'].values()))  # Sequence-dependent times of jobs to machines
    deadlines = np.array(data["Deadlines"])[0:]  # Due times of jobs
    setup_minus = np.where(setup_times >= 0, setup_times, np.inf).min(axis=1) # Lower bounds of setup times

    return nJobs, nMachines, alpha, tau, rho, processing_times, setup_times, deadlines, setup_minus