# congruence-classes
An R package to explore congruence classes. Here we list what is still missing to do.

# Package

## Documentation
Not all arguments to all functions have explanations, and many functions have TODO as their example.
We may not need examples, but we do have primate data in the package now so it should be easy to generate these using code from the vignette.

## Reading in empirical results

- **read.TreePar:** a single function to read in TreePar code and produce rate functions, to simplify user interaction. This function could possibly also perform the TreePar analysis, but I think that having at least one separate function that takes a TreePar result and then creates the functions would be ideal.
- **read.TESS:** a single function to read in TESS output and produce rate functions, to simplify user interaction.

## Posterior samples instead of singe point estimate

- Update **read.TESS** to construct a list of rate functions (one function per posterior sample) instead of using the mean/median rate.

# Vignette

## Explanation of package structure and functions

This is in progress. The basic idea is that
1. users know what the expected workflow is
2. users know what are the main functions and options within the function.
Overall, the idea is that we make clear what we had in mind when we designed the package, and why we designed it this way. So simply explain our strategy and later show this with the examples.

AFM comment: If we swap chapters 3 (reading in data) and 4 (congruence class background) then chapters 4-6 show the workflow, which is as I see it,
1. Get in data
2. Explore the pulled rates and maybe a few alternatives
3. Automated exploration of many congruent scenarios.

AFM comment: I think the biggest TODO left in the vignette is documenting the basic functions.


## More motivation about generating rate functions

Why do you we have different methods to generate rates functions and not only draw rates uniformly? This needs to be explained more clearly. For example, this is not clear for section 6.2

## Examples using posterior distributions instead of point estimates.

Once these features are available, we should show a simple example.
