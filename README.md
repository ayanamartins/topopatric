# topopatric

This is a reimplementation of the model described in
de Aguiar, M. A. M., M. Baranger, E. M. Baptestini, L. Kaufman, and Y. Bar-Yam. 2009
Global Patterns of Speciation and Diversity. Nature 460, no. 7253: 384–87. doi:10.1038/nature08168.

**General description:** individuals mate according to two threshols, genetic and spatial. The emergence
of reproductively isolated groups occurs under certain conditions

**General features:**
* Individuals are haploid and hermaphrodite
* Space is discrete (lattice)

##Notable differences in relation to the Aguiar's version:

1. Individuals disperse within the mating radius
2. Generations can be discrete or continuous
3. Genetic incompatibilites can be intralocus or interlocus
4. The genetic threshold is relative (0 to 1). If it is zero, inidividuals
will only mate to others with identical genotypes

##At least three files are required to run the program:

1. Compiled program (topopatric.f90)
2. input.in (input parameters)
3. seed.in (should cotain one random seed for each model run)

In addition, the file par.in is also required to run multiple parameter sets in a single run.
