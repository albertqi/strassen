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

In this write-up, we determine the optimal cross-over point from Strassen multiplication to conventional multiplication. We will first find the analytical cross-over point and then compare that to our experimental cross-over point. Using this, we can also estimate the number of triangles in a random graph for various probabilities of including an edge. Finally, we will discuss our experiments in more depth, covering the optimizations of our algorithm and other intruiging details we discovered along the way.

## 2. Analytical Cross-Over Point

Before experimentally finding the cross-over point, it is helpful to first determine the cross-over point analytically. That is, we can mathematically calculate the value of $n_0$ that optimizes the runtime of Strassen's algorithm. We will begin by calculating the runtimes for both the conventional matrix multiplication algorithm and Strassen's multiplication algorithm.

Let $T_C(n)$ be the runtime of the conventional algorithm on an $n\times n$ matrix. When performing conventional multiplication on an $n\times n$ matrix, we need to update $n^2$ entries. Then, for each entry, there are $n$ multiplications and $n-1$ additions, meaning:

$$\begin{aligned}
T_C(n)&=n^2(n+(n-1))\\
&=n^2(2n-1)
\end{aligned}$$

Now, let $T_S(n)$ be the runtime of Strassen's algorithm on an $n\times n$ matrix. When performing Strassen multiplication on an $n\times n$ matrix, we need to solve $7$ subproblems, each of size $n/2$. Additionally, there are $18$ total additions and subtractions, each on $n/2\times n/2$ matrices. This means that:

$$\begin{aligned}
T_S(n)&=7\cdot T_S\left(\frac{n}{2}\right)+18\left(\frac{n}{2}\right)^2\\
&=7\cdot T_S\left(\frac{n}{2}\right)+\frac{9n^2}{2}
\end{aligned}$$

To then find the optimal cross-over point, we want to calculate the value of $n$ such that the time to perform conventional multiplication equals the time to perform Strassen multiplication with all subproblems using the conventional multiplication algorithm. For even values of $n$, we have that:

$$\begin{gathered}
T_C(n)=7\cdot T_C\left(\frac{n}{2}\right)+\frac{9n^2}{2}\\
n^2(2n-1)=7\left(\frac{n}{2}\right)^2\left(2\left(\frac{n}{2}\right)-1\right)+\frac{9n^2}{2}\\
n^2(2n-1)=\frac{7n^2(n-1)}{4}+\frac{9n^2}{2}\\
4n^2(2n-1)=7n^2(n-1)+18n^2\\
n^3-15n^2=0\\
n=0,15
\end{gathered}$$

We do not care about the degenerate case where $n=0$, so for even dimensions, we get that $n_0=15$.

For odd values of $n$, we first pad the matrix to be $(n+1)\times(n+1)$ in size. As such, we know that:

$$\begin{gathered}
T_C(n)=7\cdot T_C\left(\frac{n+1}{2}\right)+\frac{9(n+1)^2}{2}\\
n^2(2n-1)=7\left(\frac{n+1}{2}\right)^2\left(2\left(\frac{n+1}{2}\right)-1\right)+\frac{9(n+1)^2}{2}\\
n^2(2n-1)=\frac{7n(n+1)^2}{4}+\frac{9(n+1)^2}{2}\\
4n^2(2n-1)=7n(n+1)^2+18(n+1)^2\\
n^3-36n^2-43n-18=0\\
n\approx 37.170
\end{gathered}$$

Thus, for odd dimensions, we get that $n_0=37$.

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

As we can see by the graph, the experimental number of triangles aligns very closely with the expected number of triangles. Indeed, the experimental number only differs from the expected number of triangles by at most $1.46$ percent (with $p=0.01$).

## 5. Discussion of Experiments


