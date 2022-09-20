import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from random import random

def f(x):
 return np.exp(-1.0/x**2 - x)			# unnormalized function used to sample
  
def get_norm():					# compute value of normalization constant
 (value,err) = quad(f, 1.0, np.inf)
 norm=1.0/value
 return norm

def p(x,norm):					# normalized distribution (not needed for sampling)
 return norm*f(x)
 
def sample_exp():				# sample exp(-x) over the domain (1,infinity). analytic.
 r = random()
 x = 1.0 + np.log(1.0/r)
 return x

def sample_dist():				# use rejection to sample the true distribution
 keep_going = True
 while keep_going:
  x = sample_exp()				# sample from exp(-x)
  r = random()
  if r<np.exp(-1.0/x**2):			# comparison function is exp(-1/x**2)
   keep_going = False
 return x

def make_histogram(norm):			# make a normalized histogram and compare to the normalized original distribution
 n=2**16
 d = []
 for i in range(n):
  d.append( sample_dist() )
 plt.figure()
 n, bins, patches = plt.hist(d,density=True,bins=100)
 plt.plot(bins,p(bins,norm))
 plt.show()
 plt.close()

def main():
 norm = get_norm()
 print("norm=",norm)
 make_histogram(norm) 
  

if __name__ == '__main__':
  main()
  
