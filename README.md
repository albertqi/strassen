# CS 124: Strassen's Algorithm 
## Albert Qi and Steve Dalla
### March 29, 2023

## Table of Contents
1. Introduction
2. Analytical Cross-Over Point
3. Experimental Cross-Over Point
4. Triangles in Random Graphs
5. Discussion of Experiments

## 1. Introduction



## 2. Analytical Cross-Over Point



## 3. Experimental Cross-Over Point



## 4. Triangles in Random Graphs

For each probability $p$, we create $10$ random graphs and calculate the average number of triangles. The following table compares these averages to the expected number of triangles, which is $\binom{1024}{3}p^3$.

| $p$    | Experimental Number of Triangles | Expected Number of Triangles |
| :-----:| :------------------------------: | :--------------------------: |
| $0.01$ | $175.8$                          | $178.4$                      |
| $0.02$ | $1419.8$                         | $1427.5$                     |
| $0.03$ | $4867.4$                         | $4817.7$                     |
| $0.04$ | $11395.3$                        | $11419.7$                    |
| $0.05$ | $22323.9$                        | $22304.1$                    |

We can then graph these results to better visualize how the experimental number of triangles compares to the expected number of triangles.

![Experimental and expected number of triangles vs. p](./assets/triangles.png)

As we can see by the graph, the experimental number of triangles aligns very closely with the expected number of triangles. Indeed, the experimental number only differs from the expected number of triangles by at most $1.46\%$ (with $p=0.01$).

## 5. Discussion of Experiments


