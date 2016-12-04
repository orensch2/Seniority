# Local Seniority

The purpose of this script is to analyze the concentration effect of local seniority structure.

Local seniority structure is a vector of length k (# of neighbours)
Assumption: each neighbour is holding only 1 type of security
In case a structure has s < k seniority layers, zeros follow at the end of the vector

We also assume in this file:
All edges are of the same size y

We measure concentration with 2 measures: entropy based measure (E=1-entropy/log(k)) and standard deviation (StdEv)

We measure each seniority structure while varying y from x/k to x (most precise is to go until x*(k-1)

The result is a graph for each measure (subplot)
Each line in the subplot is pertaining to 1 seniority structure
Each point on the line is (y/x, measure)

Result:
Observe that the measures are non-decreasing on y
The measure a=(a_1,...,a_l, 0 ... 0) becomes flat for y/x > x/a_l
For 2 measures a,a' with same a_l if it is the case: measure(a) >= measure(a') at some point, then it is the case for all points
