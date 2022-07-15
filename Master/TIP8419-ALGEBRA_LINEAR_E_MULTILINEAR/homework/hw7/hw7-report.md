<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 7 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [High Order Singular Value Decomposition (HOSVD)](#high-order-singular-value-decomposition-hosvd)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# High Order Singular Value Decomposition (HOSVD)

## Problem 1
For a third-order tensor $\mathbb{X} \in \mathbb{C}^{I\times J\times K}$ implement the truncated high-order singular value decomposition (HOSVD), using the following prototype function:

$$\begin{equation} 
\left[\mathcal{S}, \mathbf{U}^{(1)}, \mathbf{U}^{(2)}, \mathbf{U}^{(3)} \right] = \text{hosvd}(\mathcal{S}) \end{equation}$$

<u>Hint</u>: Use the file “hosvd_test.mat” to validate your results.

---

## Problem 2

Consider the two third-order tensors $\mathcal{X} \in \mathbb{C}^{8 \times 4 \times 10}$ and $\mathcal{Y} \in \mathbb{C}^{5 \times 5 \times 5}$ provided in the data file “hosvd_denoising.mat”. By using your HOSVD prototype function, find a low multilinear rank approximation for these tensors, defined as 

$\tilde{\mathcal{X}} \in \mathbb{C}^{R1 \times R2 \times R3}$ and $\tilde{\mathcal{Y}} \in \mathbb{C}^{P1 \times P2 \times P3}$. Then, calculate the normalized mean square error (NMSE) between the original tensor and its approximation, i.e.,:

$$\begin{equation*} 
\text{NMSE}(\tilde{\mathcal{X}}) =  \frac{|| \tilde{\mathcal{X}} - \mathcal{X} ||_{F}^{2}}{|| \mathcal{X} ||_{F}^{2}} \quad , \quad
\text{NMSE}(\tilde{\mathcal{Y}}) =  \frac{|| \tilde{\mathcal{Y}} - \mathcal{Y} ||_{F}^{2}}{|| \mathcal{Y} ||_{F}^{2}} 
\end{equation*}$$

<u>Hint</u>: The multilinear ranks of X and Y can be found by analysing the profile of the 1-mode, 2-mode and 3-mode singular values of these tensors.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that uses the SVD was applied to the initial factor matrices $\mathbf{A}_{0}$ and $\mathbf{B}_{0}$ initializated from a Normal distribution $\mathcal{N}(0,\,1)\,$;

**Discussion**

To compare the real data with the estimated factors, we may use two main results in the Experiment and Validation sections. The NMSE between the given data and obtained as output to LSKronF. As well the row/column factor scaling, i.e, apply the element-wise division between the given data and algorithm output for $\mathbf{A}$ vs $\mathbf{\hat{A}}$, and $\mathbf{B}$ vs $\mathbf{\hat{B}}$.

The NMSE value, with an emphasis to $\text{NMSE}(\mathbf{X_{0}, \hat{X}})$ value, with a very low SNR.

We can see that for all columns are composed by the same real value, for both $\mathbf{A}$ and $\mathbf{B}$. Hence, it presents the second evidence to confirm the proper algorithm estimation, since the columns differs only by a scale factor.

[Problem 1 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

---


### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 1000 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- SNR range $\{0, 5, 10, 15, 20, 25, 30\}$;
- Compute the LSKronF for each value, for $(I, J) = (2, 4)$, $(P, Q) = (3, 5)$ and $(I, J) = (4, 8)$, $(P, Q) = (3, 5)$.

**Discussion**

The results are consistent with the experiment perfomed, that for randomly generated $\mathbf{A}$ and $\mathbf{B}$, what confirmed as shown in the previous part, the columns from given to estimated data differs only by a scale factor.

From the figure results, we may assess the SNR gap between the NMSE curves.

For each value of SNR, respectively:

	Mean Diff d1 vs d2: 6.15 dB 
	Mean Diff d1 vs d2: 5.81 dB 
	Mean Diff d1 vs d2: 5.84 dB 
	Mean Diff d1 vs d2: 5.84 dB 
	Mean Diff d1 vs d2: 5.72 dB 
	Mean Diff d1 vs d2: 5.80 dB 
	Mean Diff d1 vs d2: 5.71 dB 

The mean value for the difference (gap) between the curves:
	
	Mean Diff: 5.84 dB 

[Problem 2 script][2] and [Figures][3].

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw4/code/figures/hw4-problem1.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>


<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw4/code/hw4.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw4/code/hw4_problem.m> (Problem 2 script)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw4/code/hw4.m> (Figures)