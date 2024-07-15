from __future__ import division
import sys
import random
import time
import math
from pyomo.environ import *
from csv import reader
import logging
import cplex # Import CPLEX solver
from cplex.exceptions import CplexSolverError
from docplex.cp.model import CpoModel # Import DOCplex for the CP Optimizer module

logging.getLogger('pyomo.core').setLevel(logging.ERROR)

def cp(nJobs, nSlots, nMachines, processing_times, big_M, setup_times, resources, deadlines):
	sTimes = [[[math.ceil(setup_times[i][j][m]) + int(bool(setup_times[i][j][m] < 0))*9999 for j in range(nJobs)] for i in range(nJobs)] for m in range(nMachines)]

	CP = CpoModel()

	processInt = {}
	for j in range(nJobs):
		for m in range(nMachines):
			start, end = (0, big_M), (0, big_M)
			size = math.ceil(processing_times[j][m])
			processInt[(j, m)] = CP.interval_var(start, end, size, optional = True, name = f"process_{j},{m}")
	setupInt = {}
	for j in range(nJobs):
		for m in range(nMachines):
			start, end = (0, big_M), (0, big_M)
			setupInt[(j, m)] = CP.interval_var(start, end, optional = True, name = f"setup_{j},{m}")
	sequence = {}
	for m in range(nMachines):
		sequence[m] = CP.sequence_var([processInt[(j, m)] for j in range(nJobs)], name = f"seq_{m}")

	tardiness = []
	for j in range(nJobs):
		tardiness.append(CP.integer_var(0, big_M, name = f"tardiness_{j}"))

	for j in range(nJobs):
		CP.add(sum(CP.presence_of(processInt[(j, m)]) for m in range(nMachines)) == 1)
		for m in range(nMachines):
			CP.add(tardiness[j] >= CP.end_of(processInt[(j, m)]) - deadlines[j])
			CP.add(CP.presence_of(processInt[(j, m)]) == CP.presence_of(setupInt[(j, m)]))
			CP.add(CP.start_at_end(processInt[(j, m)], setupInt[(j, m)]))
	for m in range(nMachines):
		CP.add(CP.no_overlap(sequence[m], sTimes[m]))
		CP.add(sum(CP.pulse(setupInt[(j, m)], 1) for j in range(nJobs)) <= resources)

	isNext = [[CP.integer_var(0, nJobs, name = f"next_{i}") for i in range(nJobs)] for j in range(nJobs)]

	for i in range(nJobs):
		CP.add(isNext[i] != i)
		for j in range(nJobs):
			for m in range(nMachines):
				CP.add(CP.if_then((CP.logical_and((CP.presence_of(processInt[(j, m)]) == 1), (CP.presence_of(processInt[(i, m)]) == 0))), (isNext[i] != j)))
				CP.add(CP.if_then((CP.start_of(processInt[(i, m)]) >= CP.start_of(processInt[(j, m)])), (isNext[i] != j)))
				CP.add(CP.if_then((isNext[i] == j), (CP.length_of(setupInt[(j, m)]) == setup_times[i][j][m])))

	CP.add(CP.minimize(sum(tardiness[j] for j in range(nJobs))))
	#CP.add(CP.minimize(sum(CP.end_of(processInt[(j, m)]) for j in range(nJobs) for m in range(nMachines))))
	
	sol = CP.solve(TimeLimit = 600, trace_log = False) # Solve the CP
	if sol:
		local_bound = 0
		for m in range(nMachines):
			for j in range(nJobs):
				if len(sol[processInt[(j, m)]]) > 0:
					local_bound += max([sol[processInt[(j, m)]][1] - deadlines[j], 0])
		return local_bound
	else:
		return big_M