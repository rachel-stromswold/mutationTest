#! /bin/python

import os
import math
import sys
import csv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from scipy import stats

N = 100
N_GENERATORS = 5

#TRIAL_NUM = 100000
TRIAL_NUM = 100000
P = 0.6447
plot_traditional = False
plot_poisson = False
plot_binomial_old = False
plot_binomial_new = False
plot_hybrid_poisson = True
plot_hybrid_binomial = True

fname = "timeDat.csv"
if len(sys.argv) > 1:
    P = sys.argv[1]
if len(sys.argv) > 2:
    fname = sys.argv[2]

seedarr=[14983, 36709, 44207, 65604, 53050, 42595, 19764, 38077, 66747, 39341, 55352, 52370, 97828, 46377, 47277, 80069, 29518, 32172, 15958, 27783]
#lenarr=[8, 12, 16, 24, 32, 48, 64]
lenarr=[32, 64]

class Histogram:
    def __init__(self):
        self.data_points = []
        self.n_samples = [0 for i in range(N_GENERATORS)]
        self.labels = ["" for i in range(N_GENERATORS)]

    def add_data(self, fname):
        with open(fname) as csvfile:
            read = csv.reader(csvfile)
            this_file_samples = [0 for i in range(N_GENERATORS)]
            for row in read:
                if row[0] == "result":
                    for j in range(N_GENERATORS):
                        #print(row)
                        self.labels[j] = row[j+1]
                else:
                    ind = -1
                    #check to see if this data point is already in the array
                    for i, arr in enumerate(self.data_points):
                        if arr[0] == int(row[0]):
                            ind = i
                            break
                    if ind == -1:
                        self.data_points.append([0 for i in range(N_GENERATORS+2)])
                        self.data_points[ind][0] = int(row[0])
                        self.data_points[ind][1] = 0

                    for j in range(N_GENERATORS):
                        self.data_points[ind][j+2] += float(row[j+1])
                        self.data_points[ind][1] += self.data_points[ind][j+2]
                        this_file_samples[j] += self.data_points[ind][j+2]

            for j in range(N_GENERATORS):
                self.n_samples[j] += this_file_samples[j]

    def plot(self, n_bars=20):
        if len(self.data_points) < n_bars:
            n_bars = len(self.data_points)
        #sort by the total number of points
        self.data_points.sort(key=lambda x: x[1], reverse=True)
        axlabels = ["" for i in range(n_bars)]
        bernoulli_data = [0 for i in range(n_bars)]
        poisson_data = [0 for i in range(n_bars)]
        binomial_old_data = [0 for i in range(n_bars)]
        binomial_new_data = [0 for i in range(n_bars)]
        hybrid_poisson_data = [0 for i in range(n_bars)]

        for j in range(n_bars):
            axlabels[j] = "{:#06b}".format(self.data_points[j][0])
            bernoulli_data[j] = self.data_points[j][2]/self.n_samples[0]
            poisson_data[j] = self.data_points[j][3]/self.n_samples[1]
            binomial_old_data[j] = self.data_points[j][4]/self.n_samples[2]
            binomial_new_data[j] = self.data_points[j][5]/self.n_samples[3]
            hybrid_poisson_data[j] = self.data_points[j][6]/self.n_samples[4]
        print(axlabels)
        x = np.arange(len(axlabels))
        width=0.9

        fig, ax = plt.subplots()
        sw = width/N_GENERATORS
        rects1 = ax.bar(x-9*sw/4, bernoulli_data, sw, label="bernoulli", color="C0")
        rects2 = ax.bar(x-5*sw/4, poisson_data, sw, label="poisson", color="C1")
        rects2 = ax.bar(x-sw/4, binomial_old_data, sw, label="binomial (old)", color="C2")
        rects4 = ax.bar(x+3*sw/4, binomial_new_data, sw, label="binomial (new)", color="C3")
        rects5 = ax.bar(x+7*sw/4, hybrid_poisson_data, sw, label="hybrid (poisson)", color="C4")

        ax.set_ylabel('Frequency')
        ax.set_title('Histogram of bitstring occurrences for different generators')
        ax.set_xticks(x)
        ax.set_xticklabels(axlabels)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        ax.legend()
        plt.show()
        plt.clf()

berndat = [0 for i in range(N)]
poissdat = [0 for i in range(N)]
bindat = [0 for i in range(N)]
bindat_new = [0 for i in range(N)]
hyp_dat = [0 for i in range(N)]
hyb_dat = [0 for i in range(N)]

#0=mean, 1=var, 2=95% confidence range
berntimes = [[0 for i in range(len(lenarr))] for j in range(3)]
poisstimes = [[0 for i in range(len(lenarr))] for j in range(3)]
bintimes = [[0 for i in range(len(lenarr))] for j in range(3)]
bintimes_new = [[0 for i in range(len(lenarr))] for j in range(3)]
hyp_times = [[0 for i in range(len(lenarr))] for j in range(3)]
hyb_times = [[0 for i in range(len(lenarr))] for j in range(3)]

f = open(fname, 'w')
f.write("bitstream length,traditional time (ns),new time (ns),new time wo recursion (ns)\n")
for k, l in enumerate(lenarr):
    i = 0
    hist = Histogram()
    for s in seedarr:
        r_fd, w_fd = os.pipe()
        pidt = os.fork()
        if (pidt == 0):
            os.dup2(w_fd, 1)
            #P = 1.0/float(l)
            args = ["./main", "-l", str(l), "-s", str(s), "-n", str(TRIAL_NUM), "-p", str(P), "-q", "-d"]
            os.execv("./main", args)
        else:
            os.wait()
            os.close(w_fd)
            out = str(os.read(r_fd, 1000))
            datstr = out.rstrip()
            lst = datstr.split()
            if i < N:
                berndat[i] = float(lst[1].strip('\n \t'))
                poissdat[i] = float(lst[3].strip('\n \t'))
                bindat[i] = float(lst[5].strip('\n \t'))
                bindat_new[i] = float(lst[7].strip('\n n\\\'\t'))
                hyp_dat[i] = float(lst[9].strip('\n n\\\'\t'))
                hyb_dat[i] = float(lst[11].strip('\n n\\\'\t'))
                berntimes[0][k] += berndat[i]
                poisstimes[0][k] += poissdat[i]
                bintimes[0][k] += bindat[i]
                bintimes_new[0][k] += bindat_new[i]
                hyp_times[0][k] += hyp_dat[i]
                hyb_times[0][k] += hyb_dat[i]
                i+=1
            os.close(r_fd)
            if l == 4:
                hist.add_data("dat.csv")
    if l == 4:
        hist.plot()
    print('detected', i, 'samples')

    berntimes[0][k] /= i
    poisstimes[0][k] /= i
    bintimes[0][k] /= i
    bintimes_new[0][k] /= i
    hyp_times[0][k] /= i
    hyb_times[0][k] /= i
    for j in range(i):
        berntimes[1][k] += (berndat[j] - berntimes[0][k])*(berndat[j] - berntimes[0][k])
        poisstimes[1][k] += (poissdat[j] - poisstimes[0][k])*(poissdat[j] - berntimes[0][k])
        bintimes[1][k] += (bindat[j] - bintimes[0][k])*(bindat[j] - bintimes[0][k])
        bintimes_new[1][k] += (bindat_new[j] - bintimes_new[0][k])*(bindat_new[j] - bintimes_new[0][k])
        hyp_times[1][k] += (hyp_dat[j] - hyp_times[0][k])*(hyp_dat[j] - hyp_times[0][k])
        hyb_times[1][k] += (hyb_dat[j] - hyb_times[0][k])*(hyb_dat[j] - hyb_times[0][k])
    berntimes[1][k] /= (i-1)
    poisstimes[1][k] /= (i-1)
    bintimes[1][k] /= (i-1)
    bintimes_new[1][k] /= (i-1)
    hyp_times[1][k] /= (i-1)
    hyb_times[1][k] /= (i-1)

    s_p = math.sqrt((bintimes_new[1][k]+bintimes[1][k])/2)
    t = (bintimes_new[0][k]-bintimes[1][k])/(s_p*math.sqrt(2/i))

    #print out data about the averages variances etc.
    print('stream length:', l)
    print('  traditional:')
    print('\tavg:', berntimes[0][k])
    print('\tvar:', berntimes[1][k])
    print('\tdev:', math.sqrt(berntimes[1][k]))
    print('  poisson:')
    print('\tavg:', poisstimes[0][k])
    print('\tvar:', poisstimes[1][k])
    print('\tdev:', math.sqrt(poisstimes[1][k]))
    print('  binomial (old):')
    print('\tavg:', bintimes[0][k])
    print('\tvar:', bintimes[1][k])
    print('\tdev:', math.sqrt(bintimes[1][k]))
    print('  binomial (new):')
    print('\tavg:', bintimes_new[0][k])
    print('\tvar:', bintimes_new[1][k])
    print('\tdev:', math.sqrt(bintimes_new[1][k]))
    print('  hybrid (poisson):')
    print('\tavg:', hyp_times[0][k])
    print('\tvar:', hyp_times[1][k])
    print('\tdev:', math.sqrt(hyp_times[1][k]))
    print('  hybrid (binomial):')
    print('\tavg:', hyb_times[0][k])
    print('\tvar:', hyb_times[1][k])
    print('\tdev:', math.sqrt(hyb_times[1][k]))
    print('  s_p =', s_p, ', t =', t, ', d.f. =', 2*i-2)

    #update the error bars
    berntimes[2][k] = 1.95996*math.sqrt(berntimes[1][k])
    poisstimes[2][k] = 1.95996*math.sqrt(berntimes[1][k])
    bintimes[2][k] = 1.95996*math.sqrt(berntimes[1][k])
    bintimes_new[2][k] = 1.95996*math.sqrt(berntimes[1][k])
    hyp_times[2][k] = 1.95996*math.sqrt(hyp_times[1][k])
    hyb_times[2][k] = 1.95996*math.sqrt(hyb_times[1][k])

    mbps_list = []
    mbps_list.append(1000 * l / berntimes[0][k])   #bernoulli
    mbps_list.append(1000 * l / poisstimes[0][k])
    mbps_list.append(1000 * l / bintimes[0][k])    #binomial old
    mbps_list.append(1000 * l / bintimes_new[0][k])#binomial new
    mbps_list.append(1000 * l / hyp_times[0][k])
    mbps_list.append(1000 * l / hyb_times[0][k])
    axlabels = ["traditional", "poisson", "old binomial", "new binomial", "hybrid (poisson)", "hybrid (bin.)"]

    fig, ax = plt.subplots()
    bars = ax.bar(range(len(mbps_list)), mbps_list)
    bars[0].set_color("C0")
    bars[1].set_color("C1")
    bars[2].set_color("C2")
    bars[3].set_color("C3")
    bars[4].set_color("C4")
    bars[5].set_color("C5")

    ax.set_ylabel('Generation speed (MBPS)')
    ax.set_title('Speed for {} bits'.format(l))
    ax.set_xticks(range(len(mbps_list)))
    ax.set_xticklabels(axlabels)
    for tick in ax.get_xticklabels():
        tick.set_rotation(15)
    for i, bar in enumerate(bars):
        width = bar.get_width()
        height = bar.get_height()
        ax.text(bar.get_x() + width/2., height + 3, '%d' % mbps_list[i], ha='center', va='bottom')
    ax.legend()
    plt.show()
    plt.clf()

    wstr = str(l)+","+str(berntimes[0][k])+","+str(bintimes[0][k])+","+str(bintimes_new[0][k])+"\n"
    f.write(wstr)
f.close()

bern_slope, bern_int, bern_r, bern_p, bern_err = stats.linregress(lenarr, berntimes[0])
poiss_slope, poiss_int, poiss_r, poiss_p, poiss_err = stats.linregress(lenarr, poisstimes[0])
bin_slope, bin_int, bin_r, bin_p, bin_err = stats.linregress(lenarr, bintimes[0])
bin_slope_new, bin_int_new, bin_r_new, bin_p_new, bin_err_new = stats.linregress(lenarr, bintimes_new[0])
hyp_slope, hyp_int, hyp_r, hyp_p, hyp_err = stats.linregress(lenarr, hyp_times[0])
hyb_slope, hyb_int, hyb_r, hyb_p, hyb_err = stats.linregress(lenarr, hyb_times[0])

#print regression data
print("\n===== REGRESSION DATA =====")
print("traditional fit:")
print("\tslope:", str(bern_slope))
print("\ty int:", str(bern_int))
print("\tr:", str(bern_r))
print("\tp:", str(bern_p))
print("\terr:", str(bern_err))
print("poisson fit:")
print("\tslope:", str(poiss_slope))
print("\ty int:", str(poiss_int))
print("\tr:", str(poiss_r))
print("\tp:", str(poiss_p))
print("\terr:", str(poiss_err))
print("binomial (old) recursion fit:")
print("\tslope:", str(bin_slope))
print("\ty int:", str(bin_int))
print("\tr:", str(bin_r))
print("\tp:", str(bin_p))
print("\terr:", str(bin_err))
print("binomial (new) recursion fit:")
print("\tslope:", str(bin_slope_new))
print("\ty int:", str(bin_int_new))
print("\tr:", str(bin_r_new))
print("\tp:", str(bin_p_new))
print("\terr:", str(bin_err_new))
print("hybrid (poisson) recursion fit:")
print("\tslope:", str(hyp_slope))
print("\ty int:", str(hyp_int))
print("\tr:", str(hyp_r))
print("\tp:", str(hyp_p))
print("\terr:", str(hyp_err))
print("hybrid (binomial) recursion fit:")
print("\tslope:", str(hyb_slope))
print("\ty int:", str(hyb_int))
print("\tr:", str(hyb_r))
print("\tp:", str(hyb_p))
print("\terr:", str(hyb_err))

'''bernline = [0 for i in range(len(lenarr))]
binline = [0 for i in range(len(lenarr))]
poissonline = [0 for i in range(len(lenarr))]
binline_new = [0 for i in range(len(lenarr))]
for j in range(len(lenarr)):
    bernline[j] = bern_slope*lenarr[j] + bern_int
    binline[j] = bin_slope*lenarr[j] + bin_int
    poissonline[j] =
    binline_new[j] = bin2_slope*lenarr[j] + bin2_int'''
linexs = [lenarr[0], lenarr[-1]]
bernline = [bern_int + bern_slope*lenarr[0], bern_int + bern_slope*lenarr[-1]]
poissonline = [poiss_int + poiss_slope*lenarr[0], poiss_int + poiss_slope*lenarr[-1]]
binline = [bin_int + bin_slope*lenarr[0], bin_int + bin_slope*lenarr[-1]]
binline_new = [bin_int_new + bin_slope_new*lenarr[0], bin_int_new + bin_slope_new*lenarr[-1]]
hybline = [hyb_int + hyb_slope*lenarr[0], hyb_int + hyb_slope*lenarr[-1]]
hypline = [hyp_int + hyp_slope*lenarr[0], hyp_int + hyp_slope*lenarr[-1]]

if plot_traditional:
    plt.errorbar(lenarr, berntimes[0], berntimes[2], 0, 's', label='traditional method', color='C0')
    plt.plot(linexs, bernline, linestyle='-.', color='C0')

if plot_poisson:
    plt.errorbar(lenarr, poisstimes[0], poisstimes[2], 0, 's', label='poisson process', color='C2')
    plt.plot(linexs, poissonline, linestyle='-.', color='C1')

if plot_binomial_old:
    plt.errorbar(lenarr, bintimes[0], bintimes[2], 0, 's', label='old binomial shuffle', color='C1')
    plt.plot(linexs, binline, linestyle='-.', color='C2')

if plot_binomial_new:
    plt.errorbar(lenarr, bintimes_new[0], bintimes_new[2], 0, 's', label='new binomial shuffle', color='C3')
    plt.plot(linexs, binline_new, linestyle='-.', color='C3')

if plot_hybrid_binomial:
    plt.errorbar(lenarr, hyp_times[0], hyp_times[2], 0, 's', label='hybrid (poisson)', color='C4')
    plt.plot(linexs, hypline, linestyle='-.', color='C4')

if plot_hybrid_poisson:
    plt.errorbar(lenarr, hyb_times[0], hyb_times[2], 0, 's', label='hybrid (binomial)', color='C5')
    plt.plot(linexs, hybline, linestyle='-.', color='C5')

titlestr = 'Average time to generate '+str(TRIAL_NUM)+ ' bit strings'
plt.title(titlestr)
plt.xlabel('Length of Bitstream (bits)')
plt.xticks(range(0, 70, 4), ['%i'%(4*k) for k in range(18)])
plt.ylabel('Average Time (ns)')
plt.legend(loc='upper left')
plt.xlim([0, 70])
plt.ylim(bottom=0)
plt.savefig('time_vs_length.png')
plt.show()
