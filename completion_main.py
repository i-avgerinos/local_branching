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
    timelimit = 3600  # Timelimit for termination of the algorithm
    big_M = 100000  # A big numeric value
    solver = 'gurobi_persistent'  # Selection of solver
    log_step = 900 # Print output after a time interval of log_step seconds
    neutralize_log_output(solver)  # Neutralisation of log outputs

    file1 = open('Instances.json') # Import instances from 'file1' - appended in the GitHub directory
    datasets = json.load(file1)
    file2 = open("Warm_start.json") # Import warm start solutions from 'file2' - appended in the GitHub directory
    solutions = json.load(file2)

    output_file = "Branch-and-Check_Results.csv"
    completion_export.create_csv(output_file) # Create a csv file to export results

    for ID in datasets.keys():
        nJobs, nMachines, alpha, tau, rho, processing_times, setup_times, setup_minus, warm_start = completion_preprocessing.input_data(datasets[ID], solutions[ID], big_M) # Load a new instance 'ID'
        '''
        nJobs : Number of jobs |J| (50, 100)
        nMachines : Number of machines |M| (5, 10)
        alpha : Group of setup times (0, 1, 2)
        tau : Subgroup of due-time strictness (0.5 - loose, 0.8 - strict)
        rho : Subgroup of due-time strictness (0.2)
        processing_times[j, m] :    Processing time of job 'j' on machine 'm' (p_{jm})
        setup_times[i, j, m] :  Setup time of job 'j' on machine 'm', if preceded by job 'i' (s_{ijm})
        setup_minus[j, m] : Lower bound of setup times (s^{-}_{jm})
        warm_start : A solution provided as warm-start
        '''
        nSlots = nJobs # nSlots : Number of slots - equal to the number of jobs
        for R in range(2):
            resources = int(((R+2)/5)*nMachines) # Number of resources : 2/5 x |M| for R = 0, 3/5 x |M| for R = 1
            # first_bound : initial upper bound, obtained by the warm-start solution
            first_bound = completion_subproblem.solve_subproblem(nJobs, nMachines, processing_times, setup_times, resources, [warm_start], big_M, big_M)
            print(f"------------------------------------------------------------------")
            print(f"Instance {ID}")
            print(f"------------------------------------------------------------------")
            for algorithm in range(2):  # algorithm = 0: Algorithm 1 - Regular | algorithm = 1 : Algorithm 2 - Local Branching
                best_lb = 0 # Initial value of incumbent lower bound
                best_ub = first_bound # Initial value of incumbent upper bound
                iteration = 1 # Initial counter of integer solutions
                next_log = 1 # After log_step x next_log seconds, the incumbent bounds are printed.
                master = completion_master.master_problem(nJobs, nSlots, nMachines, processing_times, big_M, setup_times, setup_minus, warm_start) # Construction of the master problem
                start_time = time.time()

                # Callback function
                def subproblem_Callback(cb_m, cb_opt, cb_where):
                    global best_ub, best_lb, iteration, next_log

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
                            ub = completion_subproblem.solve_subproblem(nJobs, nMachines, processing_times, setup_times, resources, [solution], big_M, best_ub)
                            if ub < best_ub: # If the incumbent upper bound is improved, the regular cut is added.
                                best_ub = ub # The incumbent upper bound is updated
                                gap = round((best_ub - best_lb) / best_ub, 4) * 100 # The incumbent value of optimality gap is updated 
                                cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, 0)) # Add Benders Cut of Algorithm 1
                                iteration += 1 # Increase the counter by 1
                            else: # If the new upper bound is worse than the incumbent one...
                                if algorithm == 0:
                                    gap = round((best_ub - best_lb) / best_ub, 4) * 100
                                    cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, algorithm)) # ... add the regular Benders cut for Algorithm 1
                                    iteration += 1 # Increase the counter by 1
                                else:
                                    neighbourhood = [solution] # ... construct a 4-OPT neighbourhood for Algorithm 2.

                                    # Compute all solutions obtained from internal swaps
                                    neighbourhood = completion_subproblem.internal_swaps(nJobs, nMachines, processing_times, setup_times, resources, solution, neighbourhood, big_M, best_ub)
                                    # Compute all solutions obtained from starting-jobs shifts
                                    neighbourhood = completion_subproblem.starting_jobs(nJobs, nMachines, processing_times, setup_times, resources, solution, neighbourhood, big_M, best_ub)
                                    
                                    thousands = math.ceil(len(neighbourhood)/1000) # To avoid unreasonably large subproblems, the neighbourhood is split in segments of 1000 solutions.
                                    for t in range(thousands):
                                        if t < thousands-1:
                                            subneighbourhood = [neighbourhood[s] for s in range(t*1001, (t+1)*1000)]
                                        if t == thousands-1:
                                            subneighbourhood = [neighbourhood[s] for s in range(t*1001, len(neighbourhood)-1)]
                                         # For each segment of solution, the local optimum is obtained.
                                        ub = completion_subproblem.solve_subproblem(nJobs, nMachines, processing_times, setup_times, resources, subneighbourhood, big_M, best_ub)
                                        if ub < best_ub: # If the incumbent upper bound is improved, update its value.
                                            best_ub = ub
                                    gap = round((best_ub - best_lb) / best_ub, 4) * 100 # Update the value of incumbent optimality gap
                                    cb_opt.cbLazy(completion_subproblem.benders_cut(x, y, master, algorithm)) # Add the Benders cut for Algorithm 2
                                    iteration += 1 # Increase the counter by 1
                    if time.time() - start_time > log_step*next_log: # If an interval of log_step seconds has passed, print output.
                        next_log += 1
                        completion_export.add_input(ID, nJobs, nMachines, resources, alpha, algorithm, iteration, math.ceil(best_lb), best_ub, int(time.time() - start_time), output_file)                                   

                opt = SolverFactory(solver) # Add solver
                opt.set_instance(master) # Add the model
                opt.set_gurobi_param('OutputFlag', 0) # Avoid printing output
                opt.set_gurobi_param('TimeLimit', timelimit) # Set the timelimit
                opt.set_gurobi_param('PreCrush', 1) # Required for the Gurobi Callback
                opt.set_gurobi_param('LazyConstraints', 1) # Required for the Gurobi Callback
                opt.set_callback(subproblem_Callback) # Add the callback
                opt.solve(tee = False, warmstart = True) # Solve the problem, considering a warm-start solution

                gap = round((best_ub - best_lb) / best_ub, 4) * 100 # Final value of optimality gap
                # Print the results of the Algorithm after the timelimit or optimality is reached.
                print(f"{iteration} :       {round(best_lb, 0)}     {round(best_ub, 0)}     {round(gap, 2)}%        {int(time.time() - start_time)}")
                # Export the final results to csv
                completion_export.add_input(ID, nJobs, nMachines, resources, alpha, algorithm, iteration, math.ceil(best_lb), best_ub, int(time.time() - start_time), output_file)