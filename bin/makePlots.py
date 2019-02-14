#! /bin/python

import os
import math
import sys
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from scipy import stats

N = 100

#TRIAL_NUM = 100000
TRIAL_NUM = 100000

fname = "timeDat.csv"
if len(sys.argv) > 1:
    fname = sys.argv[1]

seedarr=[14983, 36709, 44207, 65604, 53050, 42595, 19764, 38077, 66747, 39341, 55352, 52370, 97828, 46377, 47277, 80069, 29518, 32172, 15958, 27783]
lenarr=[4, 6, 8, 12, 16, 24, 32, 48, 64]
#lenarr=[4, 6, 8]

berndat = [0 for i in range(N)]
bindat = [0 for i in range(N)]
bin2dat = [0 for i in range(N)]

#0=mean, 1=var, 2=95% range
berntimes = [[0 for i in range(len(lenarr))] for j in range(3)]
bintimes = [[0 for i in range(len(lenarr))] for j in range(3)]
bin2times = [[0 for i in range(len(lenarr))] for j in range(3)]

k = 0
f = open(fname, 'w')
f.write("bitstream length,traditional time (ns),new time (ns),new time wo recursion (ns)\n")
for l in lenarr:
    i = 0
    for s in seedarr:
        r_fd, w_fd = os.pipe()
        pidt = os.fork()
        if (pidt == 0):
            os.dup2(w_fd, 1)
            #P = 1.0/float(l)
            P = 0.1
            args = ["./testPrg", "-l", str(l), "-s", str(s), "-n", str(TRIAL_NUM), "-p", str(P), "-q"]
            os.execv("./testPrg", args)
        else:
            os.wait()
            out = str(os.read(r_fd, 100))
            datstr = out.rstrip()
            lst = datstr.split()
            if i < N:
                berndat[i] = float(lst[2])
                bindat[i] = float(lst[4])
                bin2dat[i] = float(lst[6])
                berntimes[0][k] += berndat[i]
                bintimes[0][k] += bindat[i]
                bin2times[0][k] += bin2dat[i]
                i+=1
    print('detected', i, 'samples')

    berntimes[0][k] /= i
    bintimes[0][k] /= i
    bin2times[0][k] /= i
    for j in range(i):
        berntimes[1][k] += (berndat[j] - berntimes[0][k])*(berndat[j] - berntimes[0][k])
        bintimes[1][k] += (bindat[j] - bintimes[0][k])*(bindat[j] - bintimes[0][k])
        bin2times[1][k] += (bin2dat[j] - bin2times[0][k])*(bin2dat[j] - bin2times[0][k])
    berntimes[1][k] /= (i-1)
    bintimes[1][k] /= (i-1)
    bin2times[1][k] /= (i-1)

    s_p = math.sqrt((bin2times[1][k]+bintimes[1][k])/2)
    t = (bin2times[0][k]-bintimes[1][k])/(s_p*math.sqrt(2/i))

    #print out data about the averages variances etc.
    print('stream length:', l)
    print('  traditional:')
    print('\tavg:', berntimes[0][k])
    print('\tvar:', berntimes[1][k])
    print('\tdev:', math.sqrt(berntimes[1][k]))
    print('  optimized:')
    print('\tavg:', bintimes[0][k])
    print('\tvar:', bintimes[1][k])
    print('\tdev:', math.sqrt(bintimes[1][k]))
    print('  optimized2:')
    print('\tavg:', bin2times[0][k])
    print('\tvar:', bin2times[1][k])
    print('\tdev:', math.sqrt(bin2times[1][k]))
    print('  s_p =', s_p, ', t =', t, ', d.f. =', 2*i-2)

    #update the error bars
    berntimes[2][k] = 1.95996*math.sqrt(berntimes[1][k])
    bintimes[2][k] = 1.95996*math.sqrt(berntimes[1][k])
    bin2times[2][k] = 1.95996*math.sqrt(berntimes[1][k])

    wstr = str(l)+","+str(berntimes[0][k])+","+str(bintimes[0][k])+","+str(bin2times[0][k])+"\n"
    f.write(wstr)
    k += 1
f.close()

bern_slope, bern_int, bern_r, bern_p, bern_err = stats.linregress(lenarr, berntimes[0])
bin_slope, bin_int, bin_r, bin_p, bin_err = stats.linregress(lenarr, bintimes[0])
bin2_slope, bin2_int, bin2_r, bin2_p, bin2_err = stats.linregress(lenarr, bin2times[0])

#print regression data
print("\n===== REGRESSION DATA =====")
print("traditional fit:")
print("\tslope:", str(bern_slope))
print("\ty int:", str(bern_int))
print("\tr:", str(bern_r))
print("\tp:", str(bern_p))
print("\terr:", str(bern_err))
print("binomial w/ recursion fit:")
print("\tslope:", str(bin_slope))
print("\ty int:", str(bin_int))
print("\tr:", str(bin_r))
print("\tp:", str(bin_p))
print("\terr:", str(bin_err))
print("binomial wo recursion fit:")
print("\tslope:", str(bin2_slope))
print("\ty int:", str(bin2_int))
print("\tr:", str(bin2_r))
print("\tp:", str(bin2_p))
print("\terr:", str(bin2_err))

bernline = [0 for i in range(len(lenarr))]
binline = [0 for i in range(len(lenarr))]
bin2line = [0 for i in range(len(lenarr))]
for j in range(len(lenarr)):
    bernline[j] = bern_slope*lenarr[j] + bern_int
    binline[j] = bin_slope*lenarr[j] + bin_int
    bin2line[j] = bin2_slope*lenarr[j] + bin2_int

plt.errorbar(lenarr, berntimes[0], berntimes[2], 0, 's', label='traditional method', color='red')
plt.plot(lenarr, bernline, label="traditional fit", linestyle='-.', color='red')

plt.errorbar(lenarr, bintimes[0], bintimes[2], 0, 's', label='new method w/ recursion', color='blue')
plt.plot(lenarr, binline, label="w/ recursion fit", linestyle='-.', color='blue')

plt.errorbar(lenarr, bin2times[0], bin2times[2], 0, 's', label='new method wo recursion', color='green')
plt.plot(lenarr, bin2line, label="wo recursion fit", linestyle='-.', color='green')

titlestr = 'Average time to generate '+str(TRIAL_NUM)+ ' mutation masks'
plt.title(titlestr)
plt.xlabel('Length of Bitstream (bits)')
plt.xticks(range(0, 70, 4), ['%i'%(4*k) for k in range(18)])
plt.ylabel('Average Time (ns)')
plt.legend(loc='upper left')
plt.xlim([0, 70])
plt.ylim(bottom=0)
plt.savefig('time_vs_length.png')
plt.show()
