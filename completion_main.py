from __future__ import division
import completion_preprocessing # Loading instances from json files
import completion_master # Constructing the master problem 'M'
import completion_subproblem # Solving subproblem, adding Benders cuts, constructing neighbourhoods
import completion_export # Exporting results to csv files
from gurobipy import GRB # Import Gurobi solver
import numpy as np
import pandas as pd
from pyomo.environ import * # Import 'Pyomo'
import json
import logging
import time
import math

logging.getLogger('pyomo.core').setLevel(logging.ERROR)

# A void model to neutralise output
def neutralize_log_output(solver):
    model = ConcreteModel()

    model.obj = Objective(rule = 0, sense = minimize)
    opt = SolverFactory(solver)
    opt.set_instance(model)
    opt.set_gurobi_param('OutputFlag', 0)

if __name__ == "__main__":
    timelimit = 7200  # Timelimit for termination of the algorithm
    big_M = 100000  # A big numeric value
    solver = 'gurobi_persistent'  # Selection of solver
    log_step = 3600 # Print output after a time interval of log_step seconds
    neutralize_log_output(solver)  # Neutralisation of log outputs

    file1 = open('Instances.json') # Import instances from 'file1' - appended in the GitHub directory
    datasets = json.load(file1)
    solutions = {}

    output_file = "Branch-and-Check_Results.csv"
    completion_export.create_csv(output_file) # Create a csv file to export results

    for ID in datasets.keys():
        nJobs, nMachines, alpha, tau, rho, processing_times, setup_times, setup_minus = completion_preprocessing.input_data(datasets[ID], big_M) # Load a new instance 'ID'
        '''
        nJobs : Number of jobs |J| (50, 100)
        nMachines : Number of machines |M| (5, 10)
        alpha : Group of setup times (0, 1, 2)
        tau : Subgroup of due-time strictness (0.5 - loose, 0.8 - strict)
        rho : Subgroup of due-time strictness (0.2)
        processing_times[j, m] :    Processing time of job 'j' on machine 'm' (p_{jm})
        setup_times[i, j, m] :  Setup time of job 'j' on machine 'm', if preceded by job 'i' (s_{ijm})
        setup_minus[j, m] : Lower bound of setup times (s^{-}_{jm})
        '''
        nSlots = nJobs # nSlots : Number of slots - equal to the number of jobs
        for R in range(2):
            resources = int(((R+2)/5)*nMachines) # Number of resources : 2/5 x |M| for R = 0, 3/5 x |M| for R = 1
            print(f"------------------------------------------------------------------")
            print(f"Instance {ID}")
            print(f"------------------------------------------------------------------")
            #algorithm:
            #0: Regular LBBD
            #1: Proposed Local Branching
            #2: Local Branching for all iterations
            #3: Local Branching w/o Domination rules
            for algorithm in range(4):
                best_lb = 0 # Initial value of incumbent lower bound
                best_ub = math.inf
                best_util = math.inf
                iteration = 1 # Initial counter of integer solutions
                next_log = 1 # After log_step x next_log seconds, the incumbent bounds are printed.
                master = completion_master.master_problem(nJobs, nSlots, nMachines, processing_times, big_M, setup_times, setup_minus) # Construction of the master problem
                start_time = time.time()
                subproblem_time = 0

                # Callback function
                def subproblem_Callback(cb_m, cb_opt, cb_where):
                    global best_ub, best_lb, iteration, next_log, subproblem_time, best_util

                    timer = time.time()
                    if cb_where == GRB.Callback.MIPNODE: # At each node, the callback updates the incumbent lower bound.
                        lb = cb_opt.cbGet(GRB.Callback.MIPNODE_OBJBND)
                        if lb > best_lb:
                            best_lb = lb

                    if cb_where == GRB.Callback.MIPSOL:  # If a new integer solution is found, the callback is triggered.
                        # var_list : All variables which will be used in the callback
                        var_list = [master.z]
                        for m in master.Machines:
                            for j in master.Jobs:
                                for i in master.Slots:
                                    var_list.append(master.x[i, j, m])

                        cb_opt.cbGetSolution(vars = var_list) # Get the values for all variables of interest

                        lb = cb_opt.cbGet(GRB.Callback.MIPSOL_OBJBND)  # Update the value of incumbent lower bound
                        if lb > best_lb:
                            best_lb = lb

                        if master.z.value < best_ub: # If the objective value of the master problem is greater/equal than the incumbent upper bound, no cuts are required.
                            solution = {} # A dictionary with the obtained solution
                            for m in master.Machines:
                                solution[m] = []
                                for i in master.Slots:
                                    for j in master.Jobs:
                                        if master.x[i, j, m].value > 0.9:
                                            solution[m].append(j)

                            x = [] # Group of variables 'x'
                            y = [] # Group of variables 'y'
                            for m in range(nMachines):
                                for i in range(len(solution[m])):
                                    x.append(master.x[nSlots - len(solution[m]) + i, solution[m][i], m])
                                    for k in range(nSlots):
                                        y.append(master.x[k, solution[m][i], m])

                            # Obtain the upper bound of the new integer solution
                            ub, util = completion_subproblem.solve_subproblem(nJobs, nMachines, processing_times, setup_times, resources, [solution], big_M, best_ub)
                            if algorithm == 0:
                                if ub < best_ub:
                                    best_ub = ub # The incumbent upper bound is updated
                                    best_util = util
                                gap = round((best_ub - best_lb) / best_ub, 4) * 100 # The incumbent value of optimality gap is updated 
                                cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, 0)) # Add Benders Cut of Algorithm 1
                                iteration += 1 # Increase the counter by 1
                            if algorithm == 1 or algorithm == 2:
                                neighbourhood = [solution] # ... construct a 4-OPT neighbourhood for Algorithm 2.

                                # Compute all solutions obtained from internal swaps
                                neighbourhood = completion_subproblem.internal_swaps(nJobs, nMachines, processing_times, setup_times, resources, solution, neighbourhood, big_M, best_ub, algorithm)
                                # Compute all solutions obtained from starting-jobs shifts
                                neighbourhood = completion_subproblem.starting_jobs(nJobs, nMachines, processing_times, setup_times, resources, solution, neighbourhood, big_M, best_ub, algorithm)
                                
                                thousands = math.ceil(len(neighbourhood)/1000) # To avoid unreasonably large subproblems, the neighbourhood is split in segments of 1000 solutions.
                                for t in range(thousands):
                                    if t < thousands-1:
                                        subneighbourhood = [neighbourhood[s] for s in range(t*1001, (t+1)*1000)]
                                    if t == thousands-1:
                                        subneighbourhood = [neighbourhood[s] for s in range(t*1001, len(neighbourhood)-1)]
                                     # For each segment of solution, the local optimum is obtained.
                                    ub, util = completion_subproblem.solve_subproblem(nJobs, nMachines, processing_times, setup_times, resources, subneighbourhood, big_M, best_ub)
                                    if ub < best_ub: # If the incumbent upper bound is improved, update its value.
                                        best_ub = ub
                                        best_util = util
                                gap = round((best_ub - best_lb) / best_ub, 4) * 100 # Update the value of incumbent optimality gap
                                cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, 1)) # Add the Benders cut for Algorithm 2
                                iteration += 1 # Increase the counter by 1
                            if algorithm == 3:
                                if ub < best_ub:
                                    best_ub = ub # The incumbent upper bound is updated
                                    gap = round((best_ub - best_lb) / best_ub, 4) * 100 # The incumbent value of optimality gap is updated 
                                    cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, 0)) # Add Benders Cut of Algorithm 1
                                    iteration += 1 # Increase the counter by 1
                                else:
                                    neighbourhood = [solution] # ... construct a 4-OPT neighbourhood for Algorithm 2.

                                    # Compute all solutions obtained from internal swaps
                                    neighbourhood = completion_subproblem.internal_swaps(nJobs, nMachines, processing_times, setup_times, resources, solution, neighbourhood, big_M, best_ub, 1)
                                    # Compute all solutions obtained from starting-jobs shifts
                                    neighbourhood = completion_subproblem.starting_jobs(nJobs, nMachines, processing_times, setup_times, resources, solution, neighbourhood, big_M, best_ub, 1)
                                    
                                    thousands = math.ceil(len(neighbourhood)/1000) # To avoid unreasonably large subproblems, the neighbourhood is split in segments of 1000 solutions.
                                    for t in range(thousands):
                                        if t < thousands-1:
                                            subneighbourhood = [neighbourhood[s] for s in range(t*1001, (t+1)*1000)]
                                        if t == thousands-1:
                                            subneighbourhood = [neighbourhood[s] for s in range(t*1001, len(neighbourhood)-1)]
                                         # For each segment of solution, the local optimum is obtained.
                                        ub, util = completion_subproblem.solve_subproblem(nJobs, nMachines, processing_times, setup_times, resources, subneighbourhood, big_M, best_ub)
                                        if ub < best_ub: # If the incumbent upper bound is improved, update its value.
                                            best_ub = ub
                                            best_util = util
                                    gap = round((best_ub - best_lb) / best_ub, 4) * 100 # Update the value of incumbent optimality gap
                                    cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, 1)) # Add the Benders cut for Algorithm 2
                                    iteration += 1 # Increase the counter by 1

                    subproblem_time += time.time() - timer
                    if time.time() - start_time > log_step*next_log: # If an interval of log_step seconds has passed, print output.
                        next_log += 1
                        completion_export.add_input(ID, nJobs, nMachines, resources, alpha, algorithm, iteration, math.ceil(best_lb), best_ub, int(time.time() - start_time), int(round(subproblem_time, 0)), best_util, output_file)                                   

                opt = SolverFactory(solver) # Add solver
                opt.set_instance(master) # Add the model
                opt.set_gurobi_param('OutputFlag', 0) # Avoid printing output
                opt.set_gurobi_param('TimeLimit', timelimit) # Set the timelimit
                opt.set_gurobi_param('PreCrush', 1) # Required for the Gurobi Callback
                opt.set_gurobi_param('LazyConstraints', 1) # Required for the Gurobi Callback
                opt.set_callback(subproblem_Callback) # Add the callback
                opt.solve(tee = False) # Solve the problem, considering a warm-start solution
                if best_lb > best_ub:
                    best_lb = best_ub

                gap = round((best_ub - best_lb) / best_ub, 4) * 100 # Final value of optimality gap
                # Print the results of the Algorithm after the timelimit or optimality is reached.
                print(f"Algorithm {algorithm+1}    |   {iteration} :       {round(best_lb, 0)}     {round(best_ub, 0)}     {round(gap, 2)}%        {int(time.time() - start_time)}")
                # Export the final results to csv
                completion_export.add_input(ID, nJobs, nMachines, resources, alpha, algorithm, iteration, math.ceil(best_lb), best_ub, int(time.time() - start_time), subproblem_time, best_util, output_file)