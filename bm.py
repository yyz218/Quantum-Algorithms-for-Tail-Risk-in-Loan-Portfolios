import numpy as np
import csv
from scipy.special import erfcinv
from scipy.special import erfc
import itertools
import pandas as pd
import bisect
import time

# find the index of number in between of two values
def find_index(lst, value):
    index = bisect.bisect_right(lst, value)
    return index

# cdf function
def cdf(p):
    return 0.5 * erfc(-p / np.sqrt(2))

# inverse cdf function
def cdfinv(p):
    return -np.sqrt(2) * erfcinv(2*p)

# convert a single binary string b into real number
def br(b):
    s = 0
    for i, bit in enumerate(b):
        ai = int(bit)
        s += 2 ** (i) * ai
    return (2*s+1)/(2**(len(b)+1))

# get all possible binary string of length n
def all_b(n):
    bits = itertools.product([0,1], repeat=n)
    return list(bits)

# convert all binary strings in the bin into a list of real values
def binary_to_real(bin):
    reals=[]
    for i in bin:
        reals.append(br(i[::-1]))
    return reals

# compute the total number of counts of each possible loss value
def binary_counts(loss):
    n = len(loss)
    result = []
    for combo in itertools.product([0, 1], repeat=n):
        total = 1
        for i, (first, second) in enumerate(loss):
            if combo[i] == 0:
                total *= first
            else:
                total *= second
        result.append(total)
    return result

# compute all possible loss values 
def binary_sums(loss):
    n = len(loss)
    result = []
    for combo in itertools.product([0, 1], repeat=n):
        total = 0
        for i, (first, second) in enumerate(loss):
            if combo[i] == 0:
                total += first
            else:
                total += second
        result.append(total)
    return result

# compute all possible loss values 
def get_losses(banks):
    loss = []
    for bank in banks:
        L_i, v_i, P_i, r_i, gamma_i, λ_i, mu_i, sigma_i, t_i, Z_i = bank
        loss.append((P_i - gamma_i * P_i, - P_i * r_i))
    return loss

# simulation process
def simulation(real_z, reals, k, loans, market):
    Z, mu, sigma = market
    delta_Z = Z * np.exp(mu + sigma * cdfinv(real_z))
    out = []
    for b in range(len(loans)):
        loan = loans[b]
        L_i, v_i, P_i, r_i, gamma_i, λ_i, mu_i, sigma_i, t_i, Z_i = loan
        x = (t_i - λ_i * delta_Z) / (Z_i * (1 - λ_i))
        # if default
        if x>0:
            u = cdf((np.log(x)-mu_i)/sigma_i)
            f = find_index(reals,u)
            out.append((f,2**k-f))
        else:
            out.append((0,2**k))
    return out

start = time.time()

# load the loan portfolio
market = []
loans = []
with open('input_file.txt', 'r') as f:
    for line in f:
        loan = [float(x) for x in line.split()]
        if len(loan)==10:
            loans.append(loan)
        else:
            market = loan

# set the precision k
k = 15
bin = all_b(k)
l = [0] * (2**len(loans))

# compute all real values and start the simulation
reals = binary_to_real(bin)
for i in reals:
    out = simulation(i, reals, k, loans, market)
    p = binary_counts(out)
    for i in range(len(l)):
        l[i] += p[i]

s=sum(l)
for i in range(len(l)):
    l[i] = f"{l[i]/s:.20f}"

# get all possible loss values and compute the probability for each loss
binary_sum = binary_sums(get_losses(loans))
all_prob = dict(zip(binary_sum, l))
sorted_dict = dict(sorted(all_prob.items()))
losses = list(sorted_dict.keys())
prob = list(sorted_dict.values())
for i in range(len(prob)):
    prob[i] = float(prob[i])
all_m = [sum(prob[:i+1])/sum(prob) for i in range(len(prob))]

# write data into file
with open('prob.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['loss','prob'])
    for i in range(len(losses)):
        writer.writerow([losses[i],all_m[i]])

end = time.time()

print(f"Total time spent is {end - start:.4f} seconds for {len(loans)} loans")