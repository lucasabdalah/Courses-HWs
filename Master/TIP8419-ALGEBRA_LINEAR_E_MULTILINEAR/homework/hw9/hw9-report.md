<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 9 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Multidimensional Least-Squares Khatri-Rao Factorization (MLS-KRF)](#multidimensional-least-squares-khatri-rao-factorization-mls-krf)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Multidimensional Least-Squares Khatri-Rao Factorization (MLS-KRF)

## Problem 1

Let $\mathbf{X} = \mathbf{A}^{(1)} \diamond \mathbf{A}^{(2)} \diamond \dots \diamond \mathbf{A}^{(N)} \in \mathbb{C}^{I_{1}I_{2}\dots I_{N} \times R}$ be a matrix generated from the Khatri-Rao product of $N$ matrices $\mathbf{A}^{(n)} \in \mathbb{C}^{I_{n} \times R}$, with $n = 1, 2, \dots, N$. Considering $N = 3$ and choosing your own values for $R$ and $I_{n}, n = 1, 2, 3$, implement the MLS-KRF algorithm
to find the estimates of $\mathbf{A}^{(1)}$, $\mathbf{A}^{(2)}$ and $\mathbf{A}^{(3)}$ by solving the following problem:

$$\begin{equation*} 
(\mathbf{\hat{A}^{(1)}}, \mathbf{\hat{A}^{(2)}}, \mathbf{\hat{A}^{(3)}}) = \underset{\mathbf{\hat{A}^{(1)}}, \mathbf{\hat{A}^{(2)}}, \mathbf{\hat{A}^{(3)}}} \min \, || \mathbf{X} - \mathbf{A}^{(1)} \diamond \mathbf{A}^{(2)} \diamond \mathbf{A}^{(3)} ||_{F}^{2} \end{equation*}$$

Compare the estimated matrices $\mathbf{\hat{A}^{(1)}}$, $\mathbf{\hat{A}^{(2)}}$ and $\mathbf{\hat{A}^{(3)}}$ with the original ones. What can you conclude? Explain the results.

<u>Hint</u>: Use the file “krf_matrix_3D.mat” to validate your result.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that uses the Khatri-Rao Factorization was applied to the initial factor matrices, initializated from a Normal distribution $\mathcal{N}(0,\,1)\,$;

**Discussion**

To compare the real data with the estimated factors, we may use two main results in the Experiment and Validation sections. The NMSE between the given data and obtained as output to MLSKRF. As well the row/column factor scaling, i.e, apply the element-wise division between the given data and algorithm output for $\mathcal{X}$ vs $\mathcal{\hat{X}}$, $\mathbf{A}$ vs $\mathbf{\hat{A}}$, $\mathbf{B}$ vs $\mathbf{\hat{B}}$ and $\mathbf{C}$ vs $\mathbf{\hat{C}}$.

NMSE with MLSKRF
	
	X and X_hat: -3.34 dB 
	A and A_hat: -0.92 dB 
	B and B_hat: 1.39 dB 
	C and C_hat: 2.22 dB 

Scale factor for X and X_hat with MLSKRF

	X_hat./X(1:160, 1): 0.23 
	X_hat./X(1:160, 2): 0.36 
	X_hat./X(1:160, 3): 0.2 
	X_hat./X(1:160, 4): 0.018 

Scale factor for A and A_hat with MLSKRF
	
	A_hat./A(1:5, 1): 0.28 
	A_hat./A(1:5, 2): 0.3 
	A_hat./A(1:5, 3): 0.27 
	A_hat./A(1:5, 4): -0.13

Scale factor for B and B_hat with MLSKRF

	B_hat./B(1:4, 1): 0.39 
	B_hat./B(1:4, 2): -0.47 
	B_hat./B(1:4, 3): 0.35 
	B_hat./B(1:4, 4): 0.15

Scale factor for B and B_hat with MLSKRF
	
	C_hat./C(1:8, 1): -0.27 
	C_hat./C(1:8, 2): 0.31 
	C_hat./C(1:8, 3): -0.26 
	C_hat./C(1:8, 4): 0.11 

The NMSE value, with an emphasis to $\text{NMSE}(\mathcal{X}, \mathcal{\hat{X}})$ value.

We can see that for all columns are composed by the same real value, for all matrices factors. Hence, it presents the second evidence to confirm the proper algorithm estimation, since the columns differs only by a scale factor.

[Problem 1 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

## Problem 2 

Assuming 1000 Monte Carlo experiments, generate $\mathbf{X}_{0} = \mathbf{A} \diamond \mathbf{B} \diamond \mathbf{C} \in \mathbb{C}^{I_{1} I_{2} I_{3} \times R}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{I_{1} \times R}$, $\mathbf{B} \in \mathbb{C}^{I_{2} \times R}$ and $\mathbf{C} \in \mathbb{C}^{I_{3} \times R}$, with $R = 4$, whose elements are drawn from a normal distribution. Let $\mathbf{X} = \mathbf{X}_{0} + \alpha V$ be a noisy version of $\mathbf{X}_{0}$, where $V$ is the additive noise term, whose elements are drawn from a normal distribution. The parameter α controls the power (variance) of the noise term, and is defined as a function of the signal to noise ratio (SNR), in dB, as follows
$$\begin{equation} 
\text{SNR}_{\text{dB}} = 10\log_{10} \left( \frac{|| \mathbf{X}_{0} ||_{F}^{2}}{|| \alpha V ||_{F}^{2}} \right) 
\end{equation}$$

Assuming the SNR range $\{0, 5, 10, 15, 20, 25, 30\}$ dB, find the estimates $\mathbf{\hat{A}}$, $\mathbf{\hat{B}}$ and $\mathbf{\hat{C}}$ via the MLS-KRF algorithm, assuming $I_{1} = 2$, $I_{2} = 3$ and $I_{3} = 4$.

Let us define the normalized mean square error (NMSE) measure as follows

$$\begin{equation} 
\text{NMSE}(\mathbf{X}_{0} )= \frac{1}{1000} \sum_{i=1}^{1000}  \frac{|| \mathbf{\hat{X}}_{0}(i) - \mathbf{X}_{0}(i) ||_{F}^{2}}{|| \mathbf{X}_{0}(i) ||_{F}^{2}}
\end{equation}$$

where $\mathbf{X}_{0}(i)$ e $\mathbf{\hat{X}}_{0}(i)$ represent the original data matrix and the reconstructed one at the $i\text{th}$ experiment, respectively. For each SNR value and configuration, plot the NMSE vs. SNR curve. Discuss the obtained results. 

<u>Note</u>: For a given SNR (dB), the parameter α to be used in your experiment is determined
from equation (1).

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 1000 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- SNR range $\{0, 5, 10, 15, 20, 25, 30\}$;
- For $I_{1} = 2$, $I_{2} = 3$ and $I_{3} = 4$;
- For $R = 4$.


**Discussion**

The results are consistent with the experiment perfomed, that for randomly generated $\mathbf{A}$, $\mathbf{B}$ and $\mathbf{C}$, what confirmed as shown in the previous part, the columns from given to estimated data differs only by a scale factor.

From the figure results, we may assess the SNR gap between the NMSE curves.

For each value of SNR, respectively:

| SNR 	| NMSE  |
| :---: | :---: |
| 0		| -3.34 |	
| 5		| -4.60 |	
| 10	| -4.80 |	
| 15	| -4.76 |	
| 20	| -4.83 |	
| 25	| -4.72 |	
| 30	| -4.84 |

[Problem 2 script][2] and [Figures][3].

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw9/code/figures/hw9-problem2.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>


<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw9/code/hw9.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw9/code/hw9_problem.m> (Problem 2 script)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw9/code/hw9.m> (Figures)