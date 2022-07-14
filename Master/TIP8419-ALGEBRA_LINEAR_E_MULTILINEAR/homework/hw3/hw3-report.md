<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 3 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Least-Squares Khatri-Rao Factorization (LSKRF)](#least-squares-khatri-rao-factorization-lskrf)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Least-Squares Khatri-Rao Factorization (LSKRF)

## Problem 1

---
Generate $\mathbf{X} = \mathbf{A} \diamond \mathbf{B} \in \mathbb{C}^{20 \times 4}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{5 \times 4}$ and $\mathbf{B} \in \mathbb{C}^{4 \times 4}$. Then, implement the Least-Squares Khatri-Rao Factorization (LSKRF) algorithm that estimate 
$\mathbf{A}$ and $\mathbf{B}$ by solving the following problem
$$\begin{equation*} 
(\mathbf{\hat{A}}, \mathbf{\hat{B}}) = \underset{\mathbf{A}, \mathbf{B}} \min || \mathbf{X} - \mathbf{A} \diamond \mathbf{B} ||_{F}^{2} \end{equation*}$$

Compare the estimated matrices $\mathbf{\hat{A}}$ and $\mathbf{\hat{B}}$ with the original ones. What can you conclude?
Explain the results.

<u>Hint</u>: Use the file “krf_matrix.mat” to validate your result.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that uses the SVD was applied to the initial factor matrices $\mathbf{A}_{0}$ and $\mathbf{B}_{0}$ initializated from a Normal distribution $\mathcal{N}(0,\,1)\,$;

**Discussion**

To compare the real data with the estimated factors, we may use two main results. The NMSE between the given data and obtained as output to LSKRF. As well the row/column factor scaling, i.e, apply the element-wise division between the given data and algorithm output for $\mathbf{A}$ vs $\mathbf{\hat{A}}$, and $\mathbf{B}$ vs $\mathbf{\hat{B}}$.

NMSE with LSKRF
	
    X and X_hat: -629.76 dB 
	A and A_hat: 10.05 dB 
	B and B_hat: 11.60 dB 

Scale factor for A and A_hat with LSKRF
	
    A_hat(:,1)./A(:,1): [-1.1; -1.1; -1.1; -1.1; -1.1]
    A_hat(:,2)./A(:,2): [0.66; 0.66; 0.66; 0.66; 0.66]
    A_hat(:,3)./A(:,3): [-1.1; -1.1; -1.1; -1.1; -1.1]
    A_hat(:,4)./A(:,4): [-0.69; -0.69; -0.69; -0.69; -0.69]

Scale factor for B and B_hat with LSKRF

	B_hat(:,1)./B(:,1): [-0.89; -0.89; -0.89; -0.89]
    B_hat(:,2)./B(:,2): [1.5; 1.5; 1.5; 1.5]
    B_hat(:,3)./B(:,3): [-0.88; -0.88; -0.88; -0.88]
    B_hat(:,4)./B(:,4): [-1.5; -1.5; -1.5; -1.5]

We can see that for all columns are composed by the same real value, for both $\mathbf{A}$ and $\mathbf{B}$. Hence, it presents the second evidence to confirm the proper algorithm estimation, since the columns differs only by a scale factor.

[Problem 1 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

## Problem 2 

Assuming 1000 Monte Carlo experiments, generate $\mathbf{X}_{0} = \mathbf{A} \diamond \mathbf{B} \in \mathbb{C}^{IJ \times R}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{I \times R}$ and $\mathbf{B} \in \mathbb{C}^{J \times R}$, with $R = 4$, whose elements are drawn from a normal distribution. 

Let $\mathbf{X} = \mathbf{X}_{0} + \alpha V$ be a noisy version of $\mathbf{X}_{0}$, where $V$ is the additive noise term, whose elements are drawn from a normal distribution. The parameter α controls the power (variance) of the noise term, and is defined as a function of the signal to noise ratio (SNR), in dB, as follows
$$\begin{equation} 
\text{SNR}_{\text{dB}} = 10\log_{10} \left( \frac{|| \mathbf{X}_{0} ||_{F}^{2}}{|| \alpha V ||_{F}^{2}} \right) 
\end{equation}$$

Assuming the SNR range $\{0, 5, 10, 15, 20, 25, 30\}$ dB, find the estimates 
$\hat{\mathbf{A}}$ and $\mathbf{\hat{B}}$ obtained with the LSKRF algorithm for the configurations $(I, J) = (10, 10)$ and $(I, J) = (30, 10)$.

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
- Compute the LSKRF for each value, for $I = \{10,30\}$.

**Discussion**

The results are consistent with the experiment perfomed, that for randomly generated $\mathbf{A}$ and $\mathbf{B}$, confirmed as shown in the previous part, the columns from given to estimated data differs only by a scale factor.

From the figure results, we may assess the SNR gap between the NMSE curves.

For each value of SNR, respectively:

	Mean Diff d1 vs d2: 3.86 dB 
	Mean Diff d1 vs d2: 3.38 dB 
	Mean Diff d1 vs d2: 3.25 dB 
	Mean Diff d1 vs d2: 3.27 dB 
	Mean Diff d1 vs d2: 3.25 dB 
	Mean Diff d1 vs d2: 3.27 dB 
	Mean Diff d1 vs d2: 3.31 dB 



The mean value for the difference (gap) between the curves:
	
	Mean Diff: 3.37 dB

[Problem 2 script][2] and [Figures][3].

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw3/code/figures/hw3-problem1.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>


<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw3/code/hw3.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw3/code/hw3_problem.m> (Problem 2 script)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw3/code/hw3.m> (Figures)