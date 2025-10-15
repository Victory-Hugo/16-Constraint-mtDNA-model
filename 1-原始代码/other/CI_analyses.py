import datetime
import multiprocessing as mp
import os
import sys
from typing import Dict, Tuple
sys.path.append('../mitochondrial_constraint/')
import calculate_oe.oe_functions


def mp_function(
		obs: float, total: int, exp: float, vp: float,
		dict: Dict[Tuple[float, float, int, float], Tuple[float, float, float]]):
	"""Function that will be parallelized in iterate_oe().

	"""
	(lower_CI, upper_CI) = \
		calculate_oe.oe_functions.calculate_CI(obs_max_het=obs, total=total, exp_max_het=exp, max_parameter=vp)
	oe_ratio = obs / exp
	dict[(obs, exp, total, vp)] = (oe_ratio, lower_CI, upper_CI)
	print("Added to dictionary for obs, exp, total, vp: ", obs, exp, total, vp)
	return dict
	

def iterate_oe():
	"""Iterate through a range of variables for confidence interval calculation.

	"""
	# using parallel processing
	manager = mp.Manager()
	dict = manager.dict()  # new each iteration
	pool = mp.Pool(mp.cpu_count())
	
	for exp in (list(range(1, 1001))):
		for multiplier in [2.0, 1.5, 1, 0.5, 0.1]:  # to model fixed oe ratio 2.0, 1.5, 1.0, 0.5, 0.1
			obs = exp * multiplier
			for multiplier in [1, 2, 10]:  # to model range of prop of exp possible, i.e. exp/total = 1.0, 0.5, and 0.1
				total = round(exp * multiplier)
				if (total >= exp) and (total >= obs):
					for vp in [2.0, 2.5, 3.0]:  # to model different varying parameters
						pool.apply_async(mp_function, args=(obs, total, exp, vp, dict))
		print("Calculated all values when expected is: ", exp)
		
	pool.close()
	pool.join()  # postpones the execution of next line of code until all processes in the queue are done.
	
	# write to file
	file = open('output_files/other/CI_analyses.txt', "w")
	header = "obs	exp	total	ratio	lower_CI	upper_CI	vp"
	file.write(header + '\n')
	for key in dict:
		file.write(
			str(key[0]) + '\t' + str(key[1]) + '\t' + str(key[2]) + '\t'
			+ str(dict[key][0]) + '\t' + str(dict[key][1]) + '\t' + str(dict[key][2]) + '\t' + str(key[3]) + '\n')
	

if __name__ == "__main__":
	for path in ['output_files/other']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Complete analyses for confidence interval validation")
	
	# Iterating through variables for confidence interval construction
	iterate_oe()
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
