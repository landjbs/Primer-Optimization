import numpy as np
import matplotlib.pyplot as plt
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


def prob_loss(prob, scaling_constant):
  """
  Args: probability occurance
  Returns: scaled probability
  """
  return prob * scaling_constant
