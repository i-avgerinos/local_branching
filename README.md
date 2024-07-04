# main
'Instances.json' :    The datasets that have been used for the experiments.   

# tardiness branch
'tardiness_main.py' :          The file to be executed; the instances are imported one-by-one and solved by Algorithms 1 and 2  

'tardiness_preprocessing.py' : For a given instance, the values of |J|, |M|, p_{jm}, s_{ijm}, d_{j} are obtained.  

'tardiness_master.py' :        Construction of the master problem M and the respective warm-start solution  

'tardiness_subproblem.py' :    Including the CP subproblem S' and functions for the generation of Benders cuts and neighbourhoods of internal swaps - starting-jobs shifts  

'tardiness_export.py' :        Exporting results to csv files  


# completion branch
'completion_main.py' :          The file to be executed; the instances are imported one-by-one and solved by Algorithms 1 and 2  

'completion_preprocessing.py' : For a given instance, the values of |J|, |M|, p_{jm}, s_{ijm} are obtained.  

'completion_master.py' :        Construction of the master problem M and the respective warm-start solution  

'completion_subproblem.py' :    Including the CP subproblem S' and functions for the generation of Benders cuts and neighbourhoods of internal swaps - starting-jobs shifts  

'completion_export.py' :        Exporting results to csv files  

