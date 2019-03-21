import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math

# SEQUENCE ANALYSIS
def calc_prob(seq_length, genome_length=12000000):
	"""
	Args: length of sequence, length of genome
	Returns: probability of sequence occuring arbitrarily in genome
	"""
	return 1 - (1 - (1/(4**seq_length)))**genome_length

def analyze_seq(seq):
	"""
	Args: string of sequence
	Returns: melting temp and prob of being found in genome
	"""
	# calc melting temp with exceptions
	numAT = seq.count("A") + seq.count("T")
	numGC = seq.count("G") + seq.count("C")
	if numAT + numGC != len(seq):
		raise ValueError("Invalid Sequence")
	tm = (numGC * 4) + (numAT * 2)
	# calc probabilitiy found in genome
	prob = calc_prob(len(seq))
	return tm, prob

# LOSS FUNCTIONS
def temp_loss(tm, temp_range=(65,75)):
	"""
	Args: current melting temp, optimal range tuple
	Returns: mean squared loss of current tm compared mean temp_range
	"""
	return (np.mean(temp_range) - tm) ** 2

def prob_loss(prob, scaling_constant=5000):
	"""
	Args: probability sequence of occurance
	Returns: scaled probability
	"""
	return prob * scaling_constant

# OPTIMIZATION FUNCTIONS
def optimize(pos_seq, prob_scale=10, vis=True):
	"""
	Args: string of max possible length for primer and scaling constant for probability weight
	Returns: primer that optimizes sum of temp_loss and probability_loss
	"""
	cur_seq = ""
	cur_loss = math.inf
	# stepwise descent and analysis
	for num_step, step_base in enumerate(pos_seq):
		test_seq = cur_seq + step_base
		step_tm, step_prob = analyze_seq(test_seq)
		step_loss = temp_loss(step_tm) + prob_loss(step_prob, prob_scale)
		# choose if cur_seq is optimized with step size 1
		if step_loss > cur_loss or num_step + 1 == len(pos_seq):
			return(cur_seq)
		else:
			cur_seq = test_seq
			cur_loss = step_loss

def visualize_descent(seq):
	"""
	Plots 3D gradient two loss funcitons
	"""
	temp_losses, prob_losses, total_losses, cur = [], [], [], ""
	for step in seq:
		# calc current losses
		tm, prob = analyze_seq(cur)
		cur_temp_loss, cur_prob_loss = temp_loss(tm), prob_loss(prob)
		# add to respective lists and progress a step
		temp_losses.append(cur_temp_loss)
		prob_losses.append(cur_prob_loss)
		total_losses.append(cur_temp_loss + cur_prob_loss)
		cur += step
	# plotting
	plt.plot(temp_losses)
	plt.plot(prob_losses)
	plt.plot(total_losses)
	plt.legend(["Temp Losses","Prob Losses","Total Losses"])
	plt.title("Optimization Metrics")
	plt.xlabel("Length Along Input")
	plt.ylabel("Loss")
	plt.show()

def find_primers(seq_dict, vis=True):
	"""
	Args: dict of possible primer sequence and reading direction of form {"ATGCA":"Forward"}
	Returns: optimized primers in list of form [forwards, backwards]
	"""
	primer_list = []
	for seq in seq_dict:
		if seq_dict[seq] == "Backward":
			checked_seq = seq[::-1]
		else:
			checked_seq = seq
		# optimize primers
		primer = optimize(checked_seq)
		primer_tm, primer_prob = analyze_seq(primer)
		# visualization
		if vis:
			print(f"{seq_dict[seq]} Primer:\n{primer}\nTm: {primer_tm} | Prob: {round(primer_prob, 10)} | Length: {len(primer)}\n{'-'*30}")
		primer_list.append(primer)
		return primer_list

visualize_descent("TGACGATGGTGATTATTTCGAACACGACGAATTGTAG")
