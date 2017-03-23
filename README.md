# InvLaBP
Coalescent simulation of locally adapted inversions

This program simulates the coalescent process of genetic sequences that are in a polymorphic chromosomal inversion. The model behind it is the one described in Guerrero et al. (<a href="http://rstb.royalsocietypublishing.org/content/367/1587/430.long">2012</a>). Briefly, the inversion is assumed to be at migration-selection balance (i.e., the 'locally adapted breakpoints model in the paper) and its only effect is on recombination: heterokaryotypes (heterozygotes for the inversion) exchange genetic material ('gene flux') at a different rate than homozygotes. We're interested in inversions that have very low gene flux, but this rate is a parameter in the simulation, so I guess you can set it to whatever you want.

###Building
The simulation needs the Boost static libraries (mostly for shared pointers).
I've compiled it in multiple versions of MacOS (binary in the `tests` folder) and Linux, with the basic:
`g++ -O3 src/*.cpp -o labp_v17`

###Running
To run, the simulation needs a parameter file (example in `test/inLABP.pars`). If you run `./labp_v17` with no extra parameters, it will look fro the default `inLABP.pars`.The first argument passed will be taken as a run identifier and the program will look for `inLABP_MYRUNID.pars`, e.g. running `./labp_v17 foo` asks for `inLABP_foo.pars`. This run identifier will also appear in the output files (more below).

###Main parameters
- Number of populations
- Population sizes
- Frequency of inversion, per population
- Inversion size. Affects amount of recombination in homokaryotypes.
- Number of sampled chromosomes, per population (can do random draws weighted by inversion frequency)
- Gene flux
- Migration rate. Migration currently only works for two populations. Fixing that is not hard.
- Number of sites with ancestral material. Simulation of (recombining) sites can be very slow depending on sample size.
- Age of the inversion. A step function-style event. Inversion goes from arbitrary frequency to zero. If not set, the inversion is assumed to be infinetely old. 

###Other options
There are also some (primitive) options for more interesting evolutionary models:
- Simulate drift trajectory to loss of inversion
- Simulate speciation event (currently one event, will merge all existing populations).
- Simulate a change in population size (bottleneck or expansion, as a coefficient of population size).

###Output
The simulation will always output mean coalescent time (from n simulations) for each site that carries ancestral material. It can also output two files: `outLABP.sites` has the simulated data (in "ms" format), and `outLABP.stats` has some summary statistics (number of segregating sites, heterozygosity, heterozygosity per population, Fst, Dxy, and Tajima's D). As with the input file, if a run identifier is given the names of the output files will be `outLABP_foo.sites` and `outLABP_foo.stats`.
