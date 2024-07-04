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
    if nJobs <= 100:
        setup_times = np.array(list(data['SetupTimes'].values()))  # Sequence-dependent times of jobs to machines
        processing_times = np.array(list(data['ProcessTimes'].values())).T  # Processing times of jobs to machines
    else:
        newDict = {}
        for i in data['SetupTimes'].keys():
            newDict[i] = []
            for j in data['SetupTimes'][i].keys():
                newDict[i].append([data['SetupTimes'][i][j][m] for m in data['SetupTimes'][i][j].keys()])
        setup_times = np.array(list(newDict.values()))  # Sequence-dependent times of jobs to machines
        newDict = {}
        for m in data['ProcessTimes'].keys():
            newDict[m] = []
            for j in data['ProcessTimes'][m].keys():
                newDict[m].append(data['ProcessTimes'][m][j])
        processing_times = np.array(list(newDict.values())).T  # Processing times of jobs to machines
    setup_minus = np.where(setup_times >= 0, setup_times, np.inf).min(axis=1) # Lower bounds of setup times

    return nJobs, nMachines, alpha, tau, rho, processing_times, setup_times, setup_minus