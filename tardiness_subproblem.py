from __future__ import division
import sys
import time
import math
import itertools
import numpy as np
import pandas as pd
import copy
import cplex # Import CPLEX solver
from cplex.exceptions import CplexSolverError
from docplex.cp.model import CpoModel # Import DOCplex for the CP Optimizer module

# Solve a CP for resources-allocation for a given neighbourhood of solutions
def solve_subproblem(nJobs, nMachines, processing_times, setup_times, deadlines, resources, neighbourhood, big_M, best_ub):
    subproblem = CpoModel()

    # Variables
    selectedSolution = [subproblem.integer_var(0, 1, name = f"isSelected_{s}") for s in range(len(neighbourhood))] # Binary variable : 1 if solution 's' is the local optimum, 0 otherwise
    setupJob = {} # Interval variables 'setupJob[(i, m)]': Setup of the job of position 'i' on machine 'm'
    for s in range(len(neighbourhood)):
        for m in range(nMachines):
            for i in range(len(neighbourhood[s][m])):
                start = (0, big_M)
                end = (0, big_M)
                if i == 0:
                    size = 0
                else:
                    size = math.ceil(setup_times[neighbourhood[s][m][i-1], neighbourhood[s][m][i], m])
                setupJob[(s, i, m)] = subproblem.interval_var(start, end, size, optional = True, name = f"setup_{s},{neighbourhood[s][m][i]},{m}")
    processJob = {} # Interval variables 'processJob[(i, m)]': Processing of the job of position 'i' on machine 'm'
    for s in range(len(neighbourhood)):
        for m in range(nMachines):
            for i in range(len(neighbourhood[s][m])):
                start = (0, big_M)
                end = (0, big_M)
                size = math.ceil(processing_times[neighbourhood[s][m][i], m])
                processJob[(s, i, m)] = subproblem.interval_var(start, end, size, optional = True, name = f"process_{s},{neighbourhood[s][m][i]},{m}")
    sequence = {} # Sequence variables 'sequence[m]': Sequence of setup and processing interval variables on machine 'm'
    for s in range(len(neighbourhood)):
        for m in range(nMachines):
            intervals = []
            for i in range(len(neighbourhood[s][m])):
                intervals.append(setupJob[(s, i, m)])
                intervals.append(processJob[(s, i, m)])
            sequence[(s, m)] = subproblem.sequence_var(intervals, name = f"sequence_{s},{m}")
    tardiness = [[subproblem.integer_var(0, big_M, name = f"tardiness_{s},{j}") for j in range(nJobs)] for s in range(len(neighbourhood))] # Integer variables : tardiness

    # Constraints
    subproblem.add(sum(selectedSolution[s] for s in range(len(neighbourhood))) == 1.0) # Only one solution of the neighbourhood is the local optimal one.
    for s in range(len(neighbourhood)):
        for m in range(nMachines):
            subproblem.add(subproblem.no_overlap(sequence[(s, m)])) # Each machine can either setup or process one job at most at the same time.
            for i in range(len(neighbourhood[s][m])):
                subproblem.add(subproblem.presence_of(setupJob[(s, i, m)]) == selectedSolution[s]) # Setup and Process interval variables are present only for the local optimal solution.
                subproblem.add(subproblem.presence_of(processJob[(s, i, m)]) == selectedSolution[s])
                subproblem.add(subproblem.start_at_end(processJob[(s, i, m)], setupJob[(s, i, m)])) # The processing starts at the end of the setup.
                subproblem.add(tardiness[s][neighbourhood[s][m][i]] >= subproblem.end_of(processJob[(s, i, m)])*subproblem.presence_of(processJob[(s, i, m)]) - deadlines[neighbourhood[s][m][i]])
                if i > 0:
                    subproblem.add(subproblem.previous(sequence[(s, m)], processJob[(s, i-1, m)], setupJob[(s, i, m)])) # The sequence of jobs remains intact, as indicated by the provided solution.
        subproblem.add(sum(subproblem.pulse(setupJob[(s, i, m)], 1) for m in range(nMachines) for i in range(len(neighbourhood[s][m]))) <= resources) # No more than 'resources' setup tasks can be executed at the same time for all machines.
    
    # Objective function
    subproblem.add(subproblem.minimize(sum(tardiness[s][j] for s in range(len(neighbourhood)) for j in range(nJobs)))) # Objective function: minimisation of total tardiness
    
    sol = subproblem.solve(TimeLimit = 15, trace_log = False) # Solve the subproblem
    if sol: # If a feasible solution is found ...
        local_bound = 0
        for s in range(len(neighbourhood)):
            if sol[selectedSolution[s]] > 0.5:
                for j in range(nJobs):
                    local_bound += sol[tardiness[s][j]]
                break
        return local_bound # ...return the local bound...
    else:
        return big_M # ...otherwise, return a big numeric value (always greater than the incumbent upper bound)

# Add Benders cuts
'''
x_list : Group of variables x
y_list : Group of variables y
model : The model of the master problem
algorithm : 0 or 1 (for Algorithm 1 or 2)
'''
def benders_cut(x_list, y_list, model, algorithm):
    if algorithm == 0:
        return model.constraints.add(sum(1 - x for x in x_list) >= 1.0) # For Algorithm 1, the regular combinatorial cut is generated.
    else:
        return model.constraints.add(sum(1 - x for x in x_list) + sum(1 - y for y in y_list) >= 3.0) # For Algorithm 2, the Local Branching cut is generated.

# A relaxation for the upper bound of a solution: assuming infinite resources
def relaxation(solution, processing_times, setup_times, deadlines):
    rb = 0 # The bound of the relaxation
    for m in solution.keys():
        timeline = 0
        if len(solution[m]) > 0:
            timeline += processing_times[solution[m][0], m]
            rb += max([0, timeline - deadlines[solution[m][0]]])
        if len(solution[m]) > 1:
            old_job = solution[m][0]
            for j in solution[m][1:]:
                timeline += processing_times[j, m] + setup_times[old_job, j, m]
                rb += max([0, timeline - deadlines[j]])
                old_job = j
    return rb # Return the value of the relaxed bound

# Execution of all possible internal swaps over a given solution
def internal_swaps(nJobs, nMachines, processing_times, setup_times, deadlines, resources, solution, neighbourhood, big_M, best_ub):
    for m in solution.keys():
        for i in solution[m]:
            for j in solution[m]:
                if i != j:
                    new_solution = copy.deepcopy(solution)
                    pos_i, pos_j = solution[m].index(i), solution[m].index(j)
                    new_solution[m][pos_i], new_solution[m][pos_j] = j, i
                    new_rb = relaxation(new_solution, processing_times, setup_times, deadlines) # Obtain the relaxed bound for the generated solution
                    if new_solution not in neighbourhood and new_rb < best_ub: # If the relaxed bound is worse than the incumbent one, the solution can be ignored as suboptimal.
                        neighbourhood.append(new_solution)           
    return neighbourhood # Return the updated neighbourgood

# Execution of all possible starting-jobs shifts over a given solution
def starting_jobs(nJobs, nMachines, processing_times, setup_times, deadlines, resources, solution, neighbourhood, big_M, best_ub):
    for m in solution.keys():
        if len(solution[m]) > 0:
            job = solution[m][0]
            for k in solution.keys():
                if k != m:
                    new_solution = copy.deepcopy(solution)
                    new_solution[m].remove(job)
                    new_solution[k].insert(0, job)
                    new_rb = relaxation(new_solution, processing_times, setup_times, deadlines) # Obtain the relaxed bound for the generated solution
                    if new_solution not in neighbourhood and new_rb < best_ub: # If the relaxed bound is worse than the incumbent one, the solution can be ignored as suboptimal.
                        neighbourhood.append(new_solution)
    return neighbourhood # Return the updated neighbourgood