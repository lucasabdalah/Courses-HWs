<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 5 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Kronecker Product Singular Value Decomposition (KPSVD)](#kronecker-product-singular-value-decomposition-kpsvd)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Kronecker Product Singular Value Decomposition (KPSVD)

## Problem 1

---
Generate a block matrix according to the following structure

$$\begin{equation*}
\mathbf{X} = \begin{pmatrix}
\mathbf{X}_{1,1} & \dots & \mathbf{X}_{1,N} \\
\vdots & \ddots & \vdots \\
\mathbf{X}_{M,1} & \dots & \mathbf{X}_{M,N} \\
\end{pmatrix}, \mathbf{X}_{i,j} \in \mathbb{C}^{P \times Q}, 1 \leq i \leq M, \, 1 \leq j \leq N \, ,
\end{equation*}$$

Implement the KPSVD for the matrix $\mathbf{X}$ by computing $\sigma_{k}$, $\mathbf{U}_{k}$, and $\mathbf{V}_{k}$ such that

$$\begin{equation*}
\mathbf{X} = \sum_{k=1}^{r_{KP}} \sigma_{k} \mathbf{U}_{k} \otimes \mathbf{V}_{k}
\end{equation*}$$
	
---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that uses the SVD was applied to estimate the original data;
- $M = N = P = Q = 3$;
- Randomly generate $\mathbf{X}_{i,j} = \text{rand}(P,Q), 1 \leq i \leq M, \, 1 \leq j \leq N \,$
- Initializated from a Normal distribution $\mathcal{N}(0,\,1)\,$.

**Discussion**

We use the experiment with the real rank to validate the algorithm, by observing the NMSE between the given data and obtained as output to KPSVD.

NMSE with KPSVD

	Original Matrix vs KPSVD estimation (full rank): = -596.10 dB 

The output present a very low NMSE value what, what may be used as evidence to confirm the proper algorithm estimation.

[Problem 1 script][1].

</p>
</div>

<!-------------------------------------------------------------------------- -->

## Problem 2

In the above problem, set $M = N = P = Q = 3$ and randomly generate $\mathbf{X}_{i,j} = \text{rand}(P,Q), 1 \leq i \leq M, \, 1 \leq j \leq N \,$. Then compute the KPSVD and the Kronecker-rank $r_{KP}$ of $\mathbf{X}$ by using your KPSVD prototype function. Consider $r \leq r_{KP}$. Compute the nearest rank-$r$ for the matrix $\mathbf{X}$.
	
---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**


- 1000 Monte Carlo Runs;
- The algorithm that uses the SVD was applied to estimate the original data;
- $M = N = P = Q = 3$;
- Randomly generate $\mathbf{X}_{i,j} = \text{rand}(P,Q), 1 \leq i \leq M, \, 1 \leq j \leq N \,$
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the KPSVD to assess rank deficiency, for $R$ in range $\{1, 2, 3, 4, 5, 6, 7, 8, 9\}$, where the matrix presents its full rank for $R = 9$.

**Discussion**

The results are consistent with the proposed scenario, since that for randomly generated $\mathbf{X}$, the algorithm succeds to obtain the lower NMSE (dB) with the known full-rank. Futhermore, we can see that as the rank decreases, the NMSE increases abruptly.

- Original Matrix vs KPSVD estimation

	| rank | NMSE (dB)|
	| :--: | :------: |
	| 1	|	-4.72	|
	| 2	|	-9.37	|
	| 3	|	-14.33	|
	| 4	|	-19.99	|
	| 5	|	-26.64	|
	| 6	|	-35.08	|
	| 7	|	-46.19	|
	| 8	|	-66.89	|
	| 9	|	-597.43	|


As we can see results, the NMSE reduces as the rank increases, however it reaches the lowest point when the true rank is applied.

[Problem 2 script][2] and [Figures][3].

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw5/code/figures/hw5-problem2.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>


<!---------------------------------------------------------------------------->


[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw5/code/hw5.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw5/code/hw5_problem.m> (Problem 2 script)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw5/code/hw5.m> (Figures)