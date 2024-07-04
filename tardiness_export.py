from __future__ import division
import numpy as np
import pandas as pd
import json
import logging
import time
import math

# Create a new csv file and print the description of columns
def create_csv(filename):
    with open(filename, 'a') as outfile:
        outfile.write("Instance;Jobs;Machines;Resources;alpha;Algorithm;Number of solutions;LB;UB;M_time;S_time;Time;Utilization\n")
    outfile.close()
# Add new input in the csv file
def add_input(instance, nJobs, nMachines, resources, alpha, algorithm, iteration, lb, ub, elapsed_time, subproblem_time, util, filename):
    with open(filename, 'a') as outfile:
        outfile.write(f"{instance};{nJobs};{nMachines};{resources};{alpha};{algorithm+1};{iteration};{lb};{ub};{elapsed_time - subproblem_time};{subproblem_time};{elapsed_time};{util}\n")
    outfile.close()