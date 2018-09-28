# genesim
Like IKEA furniture, each person carries around their own building
instructions. Instead of a manual, we have DNA, a string of 3 billion
[nucleotides](https://www.google.com/search?hl=en&q=nucleotides). When
humans (or other
[diploid](https://en.wikipedia.org/wiki/Ploidy#Diploid) species) mate,
the parent strings get randomly interlaced and turned into building instructions
for the offspring.
A point at which we "snip out" a substring to interlace is called a *crossover*.
We can think about how chunks of genetic code, or
substrings, get passed around a family tree by this process, and from
looking at shared substrings, make guesses about how people are
related.

This program **genesim.py** simulates the random interlacing of
genetic code on a user-specified family tree structure.
We model random interlacing of parent strings as a spatial
[Poisson process](https://en.wikipedia.org/wiki/Poisson_point_process),
which can be simulated by:
1. choosing the total number of crossovers, which will be
[Poisson distributed](https://en.wikipedia.org/wiki/Poisson_distribution)
according to the
[genetic length](https://en.wikipedia.org/wiki/Centimorgan) of the
chromosome of interest;
2. distributing crossovers uniformly at random, once again with
respect to the genetic distance.

Note that genetic distance
