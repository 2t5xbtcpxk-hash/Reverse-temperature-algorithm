# Reverse-temperature-algorithm

## The purpose of this algorithm

The purpose of this algorithm is to calculate hypotethical Fo91 crystallisation temperatures for olivines under Fo91 with the method detailed in (MATTHEWS, S., WONG, K., SHORTTLE, O., EDMONDS, M. & MACLENNAN, J. 2021. Do olivine crystallization temperatures faithfully record mantle temperature variability? Geochemistry, Geophysics, Geosystems, 22, e2020GC009157.). Slight variations were made compared to the calculations presented by Matthews and the principles in the calculations are given in detail at: (Link to Helda repository to be updated)


## How to cite

If you use this program or a modified version of it please cite the thesis (Link to Helda repository to be updated) according to the appropriate citation style of your work.

## How to use

There is a demo named `Fo91_extrapolation_demo.R` in the repository. It uses the `cat_olivine_demo.csv` file and calculates the Fo91 crystallisation temperature for a single sample. You can use this demo as a guideline as how this algorithm works. Any csv file that is given for the algorithm must contain the composition of olivines in cation proportions and have the same column names as in the demo file. If you miss data from e.g. some element just fill that column data with NA values. 
