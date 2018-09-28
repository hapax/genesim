# genesim
Like IKEA furniture, each person carries around their own building
instructions. Instead of a manual, we have DNA, a string of 3 billion
[nucleotides](https://www.google.com/search?hl=en&q=nucleotides). When
humans (or other
[diploid](https://en.wikipedia.org/wiki/Ploidy#Diploid) species) mate,
the parent strings get randomly interlaced and turned into building instructions
for the offspring.
We can think about how chunks of genetic code, or
substrings, get passed around a family tree by this process, and from
looking at shared substrings, make guesses about how people are
related.

This program **genesim.py** simulates the random interlacing of
genetic code on a user-specified family tree structure.
A point at which we "snip out" a substring to interlace is called a *crossover*.
A good model for random interlacing is a spatial
[Poisson process](https://en.wikipedia.org/wiki/Poisson_point_process),
which can be simulated by:
1. choosing the total number of crossovers, which will be
[Poisson distributed](https://en.wikipedia.org/wiki/Poisson_distribution)
according to the
[genetic length](https://en.wikipedia.org/wiki/Centimorgan) of the
chromosome of interest;
2. distributing crossovers uniformly at random, once again with
respect to the genetic distance.

Note that genetic distance is different from physical distance in
nucleotides, since the degree of physical correlation changes as we
move along a chromosome.
Locally stretching and shrinking our "ruler" along the chromosome, to
ensure that crossovers are distributed uniformly, is exactly the
definition of genetic distance.

To complement the simulation, **kinship.py** provides some mathematical tools for analysing
relatedness between individuals on pedigrees.
These can be used to check the simulation, and vice versa!

### Pedigrees

A sample pedigree sample, **pedsample.txt**, is included.
The format is simple.
Each row is an individual, specified by four numbers:
1. A label. Usually, we number individuals consecutively, starting
   with one.
2. The label of the father of the individual.
3. The label of the mother.
4. The sex of the individual, with '1' corresponding to male and '2'
to female.

Note that *founders* of a pedigree are individuals whose parents are
not included.
We indicate this by placing a 0 in both parent columns.
Our sample pedigree reads:

| Label | Father | Mother | Sex |
| --- | --- | --- | --- |
|1	|0	|0	|1|
|2	|0	|0	|2|
|3	|1	|2	|1|
|4	|1	|2	|1|
|5	|1	|2	|2|
|6	|0	|0	|1|
|7	|0	|0	|2|
|8	|3	|7	|1|
|9	|6	|5	|2|

We can picture the pedigree (by convention, squares are males and
circle females:)

 <figure>
    <div style="text-align:center"><img src ="/pedigree.png" width="400" />
    <figcaption><i></i></figcaption>
	</div>
</figure>

### References

- [*Mathematical and Statistical Methods for Genetic Analysis*](https://www.springer.com/gp/book/9780387953892)
  (2002), Kenneth Lange. A clear and thorough book, with chapter 5
  relevant to the kinship measures here.
- ["Relatedness and the X chromosome"](https://hapax.github.io/assets/x-chromosome.pdf)
  (2012), David Wakeham. A lab presentation from my time in the
  [Bahlo lab](https://www.wehi.edu.au/people/melanie-bahlo) at [WEHI](https://www.wehi.edu.au/).
