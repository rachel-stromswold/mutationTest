import sys
import numpy as np
from matplotlib import pyplot as plt
import csv
import random

fname = sys.argv[1]
discard = 1
'''span = 1
if len(sys.argv) > 2:
    span = int(sys.argv[2])'''
if len(sys.argv) > 3:
    discard = int(sys.argv[3])

'''n_samples = 100
t_max = 10000
bitstrings = [[0 for i in range(t_max)] for i in range(n_samples)]
test_occ = [0 for i in range(t_max)]
for i in range(n_samples):
    bitstrings[i][0] = 1

for i in range(n_samples):
    for j in range(t_max):
        count = 0
        tmp_bitstring = bitstrings[i]
        for k in range(j+1):
            if (tmp_bitstring[k] == 1 and random.uniform(0,1) < 0.6447) or (k > 0 and (tmp_bitstring[k-1] == 1 and random.uniform(0,1) < 0.6447) ):
                bitstrings[i][k] = 1
                count += 1
            else:
                bitstrings[i][k] = 0
        test_occ[j] += float(count)/n_samples'''

names = ['average', 'bernoulli', 'poisson', 'binomial (old)', 'binomial (new)', 'finite digit (poisson)', 'finite digit (binomial)']
n_tests = len(names)
occ_arrs = [[] for i in range(n_tests)]
rho_arrs = [[] for i in range(n_tests)]
net_offset = 1
occ_offset = 0
rho_offset = 1
span = 3
with open(fname) as csvfile:
    read = csv.reader(csvfile)
    for row in read:
        if row[0] != "time":
            for i, x in enumerate(row):
                if (i-net_offset) % span == occ_offset and (i-net_offset)/span < n_tests:
                    occ_arrs[int((i-net_offset)/span)].append(float(x))
                if (i-occ_offset) % span == rho_offset and (i-net_offset)/span < n_tests:
                    rho_arrs[int((i-net_offset)/span)].append(float(x))

n = len(occ_arrs[0])

ln_time = np.log(np.arange(n-1))
ln_time[0] = 0

#ln_test_time = np.log(np.arange(t_max))
#ln_test_time = ln_time[:,np.newaxis]
#ln_test_occ = np.log(np.array[test_occ:])
#ln_time = np.log(np.array(time[discard:]))
ln_time = ln_time[:,np.newaxis]
ln_occ = [np.log(np.array(occ_arrs[i][discard:])) for i in range(n_tests)]
ln_rho = [np.log(np.array(rho_arrs[i][discard:])) for i in range(n_tests)]
#ln_rho = np.log(np.array(rho[discard:]))

print(len(ln_time),len(ln_occ[0]),len(ln_rho[0]))
theta, _, theta_sumsq, _ = np.linalg.lstsq(ln_time, ln_occ[0])
beta, _, beta_sumsq, _ = np.linalg.lstsq(ln_time, ln_rho[0])
#theta_test, _, theta_test_sumsq, _ = np.linalg.lstsq(ln_test_time, ln_test_occ)

theta_r2 = 1 - theta_sumsq/np.mean(ln_occ[0])
beta_r2 = 1 - beta_sumsq/np.mean(ln_rho[0])
#theta_test_r2 = 1 - theta_test_sumsq/np.mean(ln_test_occ)

print("data | power law\tr2     \tr2 adjusted")
print("theta: {0:.3f}\t\t{0:.3f}\t{0:.3f}".format(theta[0], theta_r2, 1-(1-theta_r2)*(n-1)/(n-2)))
print("beta:  {0:.3f}\t\t{0:.3f}\t{0:.3f}".format(beta[0], beta_r2, 1-(1-beta_r2)*(n-1)/(n-2)))
#print("thet2: {0:.3f}\t\t{0:.3f}\t{0:.3f}".format(theta_test[0], theta_test_r2, 1-(1-theta_test_r2)*(n-1)/(n-2)))

time = np.arange(n-1)
occ_fitted = np.power(time, theta)
rho_fitted = np.power(time, beta)
#occ_test_fitted = np.power(np.array(test_time), theta_test)

for i, name in enumerate(names):
    print(discard, len(rho_arrs[i]))
    plt.loglog(time, rho_arrs[i][discard:], label=name)
plt.loglog(time, rho_fitted, label='best fit')
plt.title("Average occupation density vs time")
plt.xlabel("time step")
plt.ylabel("rho(t)")
plt.legend(names)
plt.xlim([1.0, 110.0])
plt.ylim([0.1, 1.5])
plt.show()
plt.clf()

for i, name in enumerate(names):
    plt.loglog(time, occ_arrs[i][discard:], label=name)
plt.loglog(time, occ_fitted, label='best fit')
plt.title("Average occupation number vs time")
plt.xlabel("time step")
plt.ylabel("n(t)")
plt.legend(names)
plt.xlim([1.0, 110.0])
plt.ylim([0.1, 50])
plt.show()

'''plt.loglog(test_time, test_occ);
plt.loglog(test_time, occ_test_fitted);
plt.title("Average occupation number vs time");
plt.xlabel("time step");
plt.ylabel("n(t)");
plt.show("occ.png");'''
plt.clf();
