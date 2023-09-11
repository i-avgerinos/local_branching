from __future__ import division
import sys
import random
import time
import math
from pyomo.environ import *
from csv import reader
import logging

logging.getLogger('pyomo.core').setLevel(logging.ERROR)

# Set the values of variables for the warm-start solution
def warm_up(nJobs, nSlots, nMachines, processing_times, setup_times, solution):
	f_x = [[[0.0 for m in range(nMachines)] for j in range(nJobs)] for i in range(nSlots)] # f_x : Values of variables x
	f_P = [[0.0 for m in range(nMachines)] for j in range(nSlots)] # f_P : Values of variables P
	f_S = [[0.0 for m in range(nMachines)] for j in range(nSlots)] # f_S : Values of variables S
	f_C = [[0.0 for m in range(nMachines)] for j in range(nSlots)] # f_C : Values of variables C

	for m in range(nMachines):
		timeline = 0
		for i in range(len(solution[m])):
			f_x[nSlots - len(solution[m]) + i][solution[m][i]][m] = 1.0
			f_P[nSlots - len(solution[m]) + i][m] = processing_times[solution[m][i]][m]
			if i == 0:
				f_S[nSlots - len(solution[m]) + i][m] = 0
				f_C[nSlots - len(solution[m]) + i][m] = f_P[nSlots - len(solution[m]) + i][m] + f_S[nSlots - len(solution[m]) + i][m]
				timeline = f_C[nSlots - len(solution[m]) + i][m]
			else:
				f_S[nSlots - len(solution[m]) + i][m] = setup_times[solution[m][i-1], solution[m][i], m]
				f_C[nSlots - len(solution[m]) + i][m] = timeline + f_P[nSlots - len(solution[m]) + i][m] + f_S[nSlots - len(solution[m]) + i][m]
				timeline = f_C[nSlots - len(solution[m]) + i][m]

	return f_x, f_P, f_S, f_C

def master_problem(nJobs, nSlots, nMachines, processing_times, big_M, setup_times, setup_minus, warm_start):
	master = ConcreteModel(name = "Master Problem")

	# Sets
	master.Slots = range(nSlots) # Set of Slots: identical with set 'J'
	master.Jobs = range(nJobs) # Set of Jobs 'J'
	master.Machines = range(nMachines) # Set of Machines 'M'

	# Variables
	master.z = Var(within = NonNegativeReals) # Auxiliary variable 'z': objective function
	master.x = Var(master.Slots, master.Jobs, master.Machines, within = Binary) # Variables 'x': 1 if job 'j' is assigned to slot 'i' of machine 'm', 0 otherwise
	master.Completion = Var(master.Slots, master.Machines, within = NonNegativeReals) # Variables 'C': Completion time of slot 'i' in machine 'm'
	master.Processing = Var(master.Slots, master.Machines, within = NonNegativeReals) # Variables 'P': Processing time of slot 'i' in machine 'm'
	master.Setup = Var(master.Slots, master.Machines, within = NonNegativeReals) # Variables 'S': Setup time of slot 'i' in machine 'm'

	# Objective function
	master.obj = Objective(rule = master.z, sense = minimize)

	f_x, f_P, f_S, f_C = warm_up(nJobs, nSlots, nMachines, processing_times, setup_times, warm_start) # Import values from the warm-start solution.
	# Adjust values from warm-start solutions
	for m in master.Machines:
		for i in master.Slots:
			master.Processing[i, m] = f_P[i][m]
			master.Setup[i, m] = f_S[i][m]
			master.Completion[i, m] = f_C[i][m]
			for j in master.Jobs:
				master.x[i, j, m] = f_x[i][j][m]
	master.z = big_M

	# Constraints
	master.constraints = ConstraintList()

	def minObjective(model): # Original objective function: minimisation of sum of completion times
		return model.z >= sum(model.Completion[i, m] for i in model.Slots for m in model.Machines)
	master.c1 = Constraint(rule = minObjective)

	def scheduleJobs(model, j): # If job 'j' is assigned to machine 'm', then it will occupy exactly one slot of 'm'.
		return sum(model.x[i, j, m] for i in model.Slots for m in model.Machines) == 1.0
	master.c2 = Constraint(master.Jobs, rule = scheduleJobs)

	def slotCapacity(model, i, m): # Each slot can be occupied by one job at most.
		return sum(model.x[i, j, m] for j in model.Jobs) <= 1.0
	master.c3 = Constraint(master.Slots, master.Machines, rule = slotCapacity)

	def slotsContinuity(model, i, m): # A slot 'i' can be occupied only if slot 'i-1' is also occupied.
		if i > 0:
			return sum(model.x[i, j, m] for j in model.Jobs) - sum(model.x[i-1, j, m] for j in model.Jobs) >= 0.0
		else:
			return Constraint.Skip
	master.c4 = Constraint(master.Slots, master.Machines, rule = slotsContinuity)

	def processJobs(model, i, m): # Processing time of slot 'i' in machine 'm' is determined by the assigned job 'j'.
		return model.Processing[i, m] == sum(processing_times[j, m]*model.x[i, j, m] for j in model.Jobs)
	master.c5 = Constraint(master.Slots, master.Machines, rule = processJobs)
	
	def completeJobs(model, i, m): # The completion time of slot 'i' in machine 'm' is determined by the completion time of slot 'i-1', added by the processing and setup times of 'i'.
		if i > 0:
			return model.Completion[i, m] == model.Completion[i-1, m] + model.Processing[i, m] + model.Setup[i, m]
		else:
			return model.Completion[i, m] == model.Processing[i, m] + model.Setup[i, m]
	master.c6 = Constraint(master.Slots, master.Machines, rule = completeJobs)

	def setupJobs(model, i, j, m): # If job 'j' is assigned to slot 'i' of machine 'm', then the setup time of this slot is determined by the succession of jobs in slots 'i-1' --> 'i'.
		if i > 0:
			return sum(setup_times[k, j, m]*model.x[i-1, k, m] for k in model.Jobs) - model.Setup[i, m] <= big_M*(1 - model.x[i, j, m])
		else:
			return Constraint.Skip
	master.c7 = Constraint(master.Slots, master.Jobs, master.Machines, rule = setupJobs)

	def boundSetup(model, i, m): # Setup time lower bound
		if i > 0:
			return model.Setup[i, m] >= sum(setup_minus[j, m]*model.x[i-1, j, m] for j in model.Jobs)
		else:
			return Constraint.Skip
	master.c8 = Constraint(master.Slots, master.Machines, rule = boundSetup)

	return master