# congruence-classes
An R package to explore congruence classes. Here we list what is still missing to do.

# Package

## Reading in empirical results

- **ACDC.read.RevBayes:** a single function to read in RevBayes HSMRF output and produce rate functions, to simplify user interaction.
- **ACDC.read.TreePar:** a single function to read in TreePar code and produce rate functions, to simplify user interaction. This function could possibly also perform the TreePar analysis, but I think that having at least one separate function that takes a TreePar result and then creates the functions would be ideal.
- **ACDC.read.TESS:**a single function to read in TESS output and produce rate functions, to simplify user interaction.

## Posterior samples instead of singe point estimate

- Update **ACDC.read.RevBayes** to construct a list of rate functions (one function per posterior sample) instead of using the mean/median rate.
- Update **ACDC.read.TESS** to construct a list of rate functions (one function per posterior sample) instead of using the mean/median rate.
- Update **congruence.class** to take a list of rate functions as an alternative. Then, the output should also be a list.
- Update **sample.congruence.class** to take a list of rate functions as an alternative. Then, we should both sample the reference rate as well as a congruent rate.
- When updating **sample.congruence.class**, **sample.rates** and **sample.basic.models** need adjusting to allow speciation rates to be set to more than one value. This probably means we have to make intermediate functions that take a vector and an index and plug them into the base functions. This could get ugly.

# Vignette

## Explanation of package structure and functions

This needs to be written. The basic idea is that
1. users know what the expected workflow is
2. users know what are the main functions and options within the function.
Overall, the idea is that we make clear what we had in mind when we designed the package, and why we designed it this way. So simply explain our strategy and later show this with the examples.


## Simpler functions to create rate samples

We should use and describe the rate functions implemented by Andy.


## Examples using posterior distributions instead of point estimates.

Once these features are available, we should show a simple example.
