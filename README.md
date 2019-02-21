# mutationTest
This program tests a new method for generating single point mutations for genetic algorithms. The program runs both the traditional method which iterates over each bit and randomly decides whether it should be flipped, and the new method which first samples from a binomial and then uniform_int distribution to produce a random mask. In order to apply the mutations, the generated mask is simply xored with the existing genome.
## Usage
### -d
The behavior of the program differs based on whether the -d (distribution test) flag is passed. If set, the program will store the number of times each bitmask was created into a csv file. This file can then be used to produce histograms to confirm that both methods result in the same probability distribution
### -p [value between 0.0 and 1.0 inclusive]
Specify the probability of mutation (the probability that any given bit in a mask has value 1).
### -n [number of trials]
Specify the number of trials for each method to be performed.
### -l [length of bitmask]
Specify the length of the bitmask to be created. Note that this option cannot be used with the -d option, which sets the length to 8 bits automatically.
### -s [seed for RNG]
Specify the seed to be used for the random number generator.
### -q
Specify quiet execution. If passed, the only output data will be 6 fields with data as follows:

\[total time for trad. method\] \[avg time for trad.\] \[total time for new method w/ recursion] \[avg time for new method w/ recursion\] \[avg time for new method wo recursion\] \[avg time for new method wo recurion\]

## Example
```
[bin]$ ./testPrg -d
Running with:
        trials:            10000
        seed:              3141593
        length:            8
        prob:              0.01
        test distribution? yes
Bernoulli time              : 11566108  avg: 1156.61
binomial time (w/ recursion): 2661783   avg: 266.178
binomial time (wo recursion): 2421754   avg: 242.175
```

```
[bin]$ ./testPrg -s 2718282 -l 16 -n 100000 -p 0.05
Running with:
        trials:            100000
        seed:              2718282
        length:            16
        prob:              0.05
        test distribution? no
Bernoulli time              : 133176005 avg: 1331.76
binomial time (w/ recursion): 36425975  avg: 364.26
binomial time (wo recursion): 34986643  avg: 249.866
```
