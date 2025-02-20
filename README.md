# Coarse-Graining Cascades Within Food Webs

Quantifying population dynamics is a fundamental challenge in ecology and evolutionary biology, particularly for species that are cryptic, microscopic, or extinct. Traditional approaches rely on continuous representations of population size, but in many cases, the precise number of individuals is unknowable. Here, we present a coarse-grained population model that simplifies population dynamics to binary states -- high or low -- determined by the balance of bottom-up resource availability and top-down predation pressure. This Boolean framework provides a minimal yet analytically tractable alternative to traditional Lotka-Volterra-based models, enabling direct insights into the role of food web structure in shaping community stability. Using this approach, we investigate how trophic interactions influence population persistence, cyclic dynamics, and extinction risk across model food webs. We find that top-down effects are a primary driver of cycling, aligning with theoretical expectations from traditional population models, and that trophic position strongly influences extinction risk, with higher-trophic species more prone to persistent low-population states. Additionally, we explore the role of trophic short-circuits -- direct interactions between apex predators and low-trophic prey -- and find that they can buffer cascades and alter extinction patterns in ways that are often overlooked in classical models. By simplifying population dynamics to a two-state system, this framework provides a powerful tool for disentangling the structural drivers of community stability. These results highlight the potential of coarse-grained approaches to complement existing models, offering new insights into trophic interactions, extinction risks, and the susceptibility of species to trophic cascades.

## Notes
The script script_figures reproduces the figures in the manuscript, but is not written incredibly efficiently, as these scripts were compiled from a larger assortment of scripts to investigate the model at different times... so there is a lot of duplication of computational effort.

Note also that Figures 5 and 6 are constructed using /math/analytical_solutions.nb as well as simulation export in the below script.

I have not included scripts for the supplementary figures.

I use the function src/smartpath() to direct figure and data files to a central directory. It's currently set for the path on my computer, but you will want to replace with the path to your file repository.
