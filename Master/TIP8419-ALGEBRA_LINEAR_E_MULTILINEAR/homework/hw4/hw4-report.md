<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 4 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Least-Squares Kronnecker Product Factorization (LSKRONF)](#least-squares-kronnecker-product-factorization-lskronf)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Least-Squares Kronnecker Product Factorization (LSKronF)

## Problem 1

---
Generate $\mathbf{X} = \mathbf{A} \otimes \mathbf{B} \in \mathbb{C}^{24 \times 6}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{4 \times 2}$ and $\mathbf{B} \in \mathbb{C}^{6 \times 3}$. Then, implement the Least-Squares Kronnecker Product Factorization (LSKronF) algorithm that estimate 
$\mathbf{A}$ and $\mathbf{B}$ by solving the following problem
$$\begin{equation*} 
(\mathbf{\hat{A}}, \mathbf{\hat{B}}) = \underset{\mathbf{A}, \mathbf{B}} \min || \mathbf{X} - \mathbf{A} \otimes \mathbf{B} ||_{F}^{2} \end{equation*}$$

Compare the estimated matrices $\mathbf{\hat{A}}$ and $\mathbf{\hat{B}}$ with the original ones. What can you conclude?
Explain the results.

<u>Hint</u>: Use the file “kronf_matrix.mat” to validate your result.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that uses the SVD was applied to the initial factor matrices $\mathbf{A}_{0}$ and $\mathbf{B}_{0}$ initializated from a Normal distribution $\mathcal{N}(0,\,1)\,$;

**Discussion**

To compare the real data with the estimated factors, we may use two main results in the Experiment and Validation sections. The NMSE between the given data and obtained as output to LSKronF. As well the row/column factor scaling, i.e, apply the element-wise division between the given data and algorithm output for $\mathbf{A}$ vs $\mathbf{\hat{A}}$, and $\mathbf{B}$ vs $\mathbf{\hat{B}}$.

- Experiment

NMSE with LSKronF
	
    X and X_hat: -622.63 dB
	A and A_hat: 9.94 dB
	B and B_hat: 1.90 dB 

Scale factor for A and A_hat with LSKronF
	
    A_hat(:,1)./A(:,1): [0.19; 0.19; 0.19; 0.19]
    A_hat(:,2)./A(:,2): [0.19; 0.19; 0.19; 0.19]

Scale factor for B and B_hat with LSKronF

	B_hat(:,1)./B(:,1): [0.076; 0.076; 0.076; 0.076; 0.076; 0.076]
    B_hat(:,2)./B(:,2): [0.076; 0.076; 0.076; 0.076; 0.076; 0.076]
	B_hat(:,3)./B(:,3): [0.076; 0.076; 0.076; 0.076; 0.076; 0.076]

- Validation Data

NMSE with LSKronF
	
    X and X_hat: -608.15 dB	
	A and A_hat: 10.49 dB
	B and B_hat: 13.74 dB 

Scale factor for A and A_hat with LSKronF
	
    A_hat(:,1)./A(:,1): [-0.83; -0.83; -0.83; -0.83]
    A_hat(:,2)./A(:,2): [-0.83; -0.83; -0.83; -0.83]
	A_hat(:,2)./A(:,2): [-0.83; -0.83; -0.83; -0.83]

Scale factor for B and B_hat with LSKronF

	B_hat(:,1)./B(:,1): [-1.2; -1.2; -1.2; -1.2]
    B_hat(:,2)./B(:,2): [-1.2; -1.2; -1.2; -1.2]

The NMSE value, with an emphasis to $\text{NMSE}(\mathbf{X_{0}, \hat{X}})$ value, with a very low SNR.

We can see that for all columns are composed by the same real value, for both $\mathbf{A}$ and $\mathbf{B}$. Hence, it presents the second evidence to confirm the proper algorithm estimation, since the columns differs only by a scale factor.

[Problem 1 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

## Problem 2 

Assuming 1000 Monte Carlo experiments, generate $\mathbf{X}_{0} = \mathbf{A} \otimes \mathbf{B} \in \mathbb{C}^{IJ \times PQ}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{I \times P}$ and $\mathbf{B} \in \mathbb{C}^{J \times Q}$, whose elements are drawn from a normal distribution. 

Let $\mathbf{X} = \mathbf{X}_{0} + \alpha V$ be a noisy version of $\mathbf{X}_{0}$, where $V$ is the additive noise term, whose elements are drawn from a normal distribution. The parameter α controls the power (variance) of the noise term, and is defined as a function of the signal to noise ratio (SNR), in dB, as follows
$$\begin{equation} 
\text{SNR}_{\text{dB}} = 10\log_{10} \left( \frac{|| \mathbf{X}_{0} ||_{F}^{2}}{|| \alpha V ||_{F}^{2}} \right) 
\end{equation}$$

Assuming the SNR range $\{0, 5, 10, 15, 20, 25, 30\}$ dB, find the estimates 
$\hat{\mathbf{A}}$ and $\mathbf{\hat{B}}$ obtained with the LSKronF algorithm for the configurations $(I, J) = (2, 4)$, $(P, Q) = (3, 5)$ and $(I, J) = (4, 8)$, $(P, Q) = (3, 5)$

Let us define the normalized mean square error (NMSE) measure as follows

$$\begin{equation} 
\text{NMSE}(\mathbf{X}_{0} )= \frac{1}{1000} \sum_{i=1}^{1000}  \frac{|| \mathbf{\hat{X}}_{0}(i) - \mathbf{X}_{0}(i) ||_{F}^{2}}{|| \mathbf{X}_{0}(i) ||_{F}^{2}}
\end{equation}$$

where $\mathbf{X}_{0}(i)$ e $\mathbf{\hat{X}}_{0}(i)$ represent the original data matrix and the reconstructed one at the $i\text{th}$ experiment, respectively. For each SNR value and configuration, plot the NMSE vs. SNR curve. Discuss the obtained results. 

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