<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 10 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Multidimensional Least-Squares Kronecker Factorization (MLS-KronF)](#multidimensional-least-squares-kronecker-factorization-mls-kronf)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Multidimensional Least-Squares Kronecker Factorization (MLS-KronF)

## Problem 1

Let $\mathbf{X} \in \mathbb{C}^{I_1 I_2 \dots I_N \times J_1 J_2 \dots J_N}$ be a matrix that we wish to approximate as $\mathbf{X} \approx \mathbf{A}^{(1)} \otimes \mathbf{A}^{(2)} \otimes \dots \otimes \mathbf{A}^{(N)}$, that is, as Kronecker product of $N$ matrices $\mathbf{A}^{(n)} \in \mathbb{C}^{I_{n} \times J_{n}}$ with $n = 1, 2, \dots, N$. For $N = 3$ and arbitrary $I_{n}$ and $J_{n}$, implement the MLSKronF algorithm that estimates  $\mathbf{A}^{(1)}$, $\mathbf{A}^{(2)}$, $\mathbf{A}^{(3)}$ by solving the following problem:

$$\begin{equation*} 
(\mathbf{\hat{A}^{(1)}}, \mathbf{\hat{A}^{(2)}}, \mathbf{\hat{A}^{(3)}}) = \underset{\mathbf{\hat{A}^{(1)}}, \mathbf{\hat{A}^{(2)}}, \mathbf{\hat{A}^{(3)}}} \min \, || \mathbf{X} - \mathbf{A}^{(1)} \otimes \mathbf{A}^{(2)} \otimes \mathbf{A}^{(3)} ||_{F}^{2} \end{equation*}$$

using either the truncated HOSVD or the HOOI initialized with the HOSVD (you should implement both versions).

Test the algorithms on a matrix that exactly follows the model. Compare the estimated matrices $\mathbf{\hat{A}^{(1)}}$, $\mathbf{\hat{A}^{(2)}}$ and $\mathbf{\hat{A}^{(3)}}$ with the original ones. What can you conclude? Explain the results.

<u>Hint</u>: Use the file “Practice_10_kronf_matrix_3D.mat” to validate your result.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that uses the Kronecker Factorization was applied to the initial factor matrices, initializated from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- HOSVD initization
- HOOI initization

**Discussion**

To compare the real data with the estimated factors, we may use the experimental results for NMSE between the given data and obtained as output to MLSKRF with HOSVD and HOOI.

NMSE with MLSKronF (HOSVD)

	X and X_hat: -0.62 dB 
	A and A_hat: 4.07 dB 
	B and B_hat: 3.92 dB 
	C and C_hat: 1.91 dB 

NMSE with MLSKronF (HOOI)

	X and X_hat: -0.62 dB 
	A and A_hat: -9.49 dB 
	B and B_hat: -9.04 dB 
	C and C_hat: -3.97 dB 

The NMSE value, with an emphasis to $\text{NMSE}(\mathcal{X}, \mathcal{\hat{X}})$ value, suport the hypothesis of a propor implementation of the algorithm.

We can see that MLSKRF initialized with HOOI outperforms HOSVD for all values, presenting smaller NMSE values also for the factor matrices.

[Problem 1 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

## Problem 2 

Assuming 1000 Monte Carlo experiments, generate $\mathbf{X}_{0} = \mathbf{A} \otimes \mathbf{B} \otimes \mathbf{C} \in \mathbb{C}^{I_{1} I_{2} I_{3} \times J_{1} J_{2} J_{3}}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{I_{1} \times J_{1}}$, $\mathbf{B} \in \mathbb{C}^{I_{2} \times J_{2}}$ and $\mathbf{C} \in \mathbb{C}^{I_{3} \times J_{3}}$, whose elements are drawn from a normal distribution. Let $\mathbf{X} = \mathbf{X}_{0} + \alpha V$ be a noisy version of $\mathbf{X}_{0}$, where $V$ is the additive noise term, whose elements are drawn from a normal distribution. The parameter α controls the power (variance) of the noise term, and is defined as a function of the signal to noise ratio (SNR), in dB, as follows
$$\begin{equation} 
\text{SNR}_{\text{dB}} = 10\log_{10} \left( \frac{|| \mathbf{X}_{0} ||_{F}^{2}}{|| \alpha V ||_{F}^{2}} \right) 
\end{equation}$$

Assuming the SNR range $\{0, 5, 10, 15, 20, 25, 30\}$ dB, find the estimates $\mathbf{\hat{A}}$, $\mathbf{\hat{B}}$ and $\mathbf{\hat{C}}$ via the MLS-KronF algorithm, assuming:

1. $(I_1, I_2, I_3) = (J_1, J_2, J_3) = (2, 2, 2)$;
2. $(I_1, I_2, I_3) = (J_1, J_2, J_3) = (5, 5, 5)$;
3. $(I_1, I_2, I_3) = (J_1, J_2, J_3) = (2, 3, 4)$;
4. $(I_1, I_2, I_3) = (2, 3, 4)$ and $(J_1, J_2, J_3) = (5, 6, 7)$.

Let us define the normalized mean square error (NMSE) measure as follows

$$\begin{equation} 
\text{NMSE}(\mathbf{X}_{0} )= \frac{1}{1000} \sum_{i=1}^{1000}  \frac{|| \mathbf{\hat{X}}_{0}(i) - \mathbf{X}_{0}(i) ||_{F}^{2}}{|| \mathbf{X}_{0}(i) ||_{F}^{2}}
\end{equation}$$

where $\mathbf{X}_{0}(i)$ e $\mathbf{\hat{X}}_{0}(i)$ represent the original data matrix and the reconstructed one at the $i\text{th}$ experiment, respectively. For each SNR value and configuration, plot the NMSE vs. SNR curve. Discuss the obtained results. 

<u>Note</u>: For a given SNR (dB), the parameter α to be used in your experiment is determined from equation (1).

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 1000 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- SNR range $\{0, 5, 10, 15, 20, 25, 30\}$;
- Assuming four scenarios as proposed.

**Discussion**

The results are consistent with the experiment perfomed, that provide small SNR values for the randomly generated  factors $\mathbf{A}$, $\mathbf{B}$ and $\mathbf{C}$.

From the figure results, we may assess the SNR gap between the NMSE curves.

For each value of SNR, respectively:

| SNR 	| Scenario 1 NMSE | Scenario 2 NMSE | Scenario 3 NMSE | Scenario 4 NMSE |
| :---: | :---: | :---: | :---: | :---: |
| 0		| 8.64e-02 | 0.00e+00  | 5.22e-04  | 1.11e-16  |
| 5		| 1.23e-03 | -3.33e-16 | 4.07e-04  | -1.11e-16  |
| 10	|  0.00e+00  | -1.11e-16  | -1.23e-04  | 0.00e+00  |	
| 15	| 1.78e-15 |	0.00e+00  | -2.22e-16 | 1.11e-16  |
| 20	| -1.37e-05 |	0.00e+00  | 0.00e+00 |  1.11e-16 |
| 25	| 8.88e-16 |	6.66e-16   | 0.00e+00 |  0.00e+00  |
| 30	| 8.88e-16 | 0.00e+00 | 0.00e+00 | -4.44e-16  |
| Mean	| 1.25e-02| 3.17e-17 | 1.15e-04 | -3.17e-17  |
	 
We can see that the HOOI initialization outperfomrs HOSVD with a small advantage 
Each experiment is implemented in: [Problem 2 - Scenario 1][2], [Problem 2 - Scenario 2][3], [Problem 2 - Scenario 3][4], [Problem 2 - Scenario 4][5].

Code that generates the figures: [Problem 2 script][6], 

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/figures/hw10-problem1.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/figures/hw10-problem2.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/figures/hw10-problem3.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/figures/hw10-problem4.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/hw10.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/hw10_problem1.m> (Problem 2 - Scenario 1)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/hw10_problem2.m> (Problem 2 - Scenario 2)
[4]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/hw10_problem3.m> (Problem 2 - Scenario 3)
[5]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/hw10_problem4.m> (Problem 2 - Scenario 4)
[6]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/code/hw10.m> (Figures)