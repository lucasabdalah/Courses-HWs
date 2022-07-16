<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 11 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Alternating Least Squares (ALS) Algorithm](#alternating-least-squares-als-algorithm)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Alternating Least Squares (ALS) Algorithm

## Problem 1

For the third-order tensor $\mathcal{X} \in \mathbb{C}^{I \times J \times K}$ provided in the file “cpd_tensor.mat”, implement the plain-vanilla Alternating Least Squares (ALS) algorithm that estimates the factor matrices $\mathbf{A} \in \mathbb{C}^{I \times R}$, $\mathbf{B} \in \mathbb{C}^{J \times R}$ and $\mathbf{C} \in \mathbb{C}^{K \times R}$ by solving the following problem:

$$\begin{equation*} 
(\mathbf{\hat{A}}, \mathbf{\hat{B}}, \mathbf{\hat{C}}) = \underset{\mathbf{{A}}, \mathbf{{B}}, \mathbf{{C}}} \min \, || \mathcal{X} - \sum_{r=1}^{R} a_{r} \circ b_{r} \circ c_{r} ||_{F}^{2} \end{equation*}$$

where $\mathbf{A} = [a_1, \dots , a_R]$, $\mathbf{B} = [b_1, \dots , b_R]$, $\mathbf{C} = [c_1, \dots , c_R]$. 

Considering a successful run, compare the estimated matrices $\mathbf{\hat{A}}, \mathbf{\hat{B}}, \mathbf{\hat{C}}$ with the original ones (also provided in the
same Matlab file). Explain the results. 

<u>Hint</u>: An error measure at the ith iteration can be calculated from the following formula:

$$\begin{equation} 
e_{(i)} = \underset{\mathbf{{A}}, \mathbf{{B}}, \mathbf{{C}}} \min \, || [\mathcal{X}]_{(1)} - \mathbf{\hat{A}}_{(i)} \left(\mathbf{\hat{C}}) \diamond \mathbf{\hat{B}}_{(i)}\right)^{\top} ||_{F} \end{equation}$$

---

### Results

NMSE for ALS Validation
  
    X and X_hat: -123.33 dB 
    A and A_hat: 27.96 dB 
    B and B_hat: 0.64 dB 
    C and C_hat: 25.72 dB 
  


<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw11/code/figures/hw11-problem1.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<!--  -->

## Problem 2

In a Monte Carlo experiment with $M = 1000$ realizations, generate a tensor
$\mathcal{X}_{(0)} = \text{CPD}(\mathbf{{A}}, \mathbf{{B}}, \mathbf{{C}})$, where  $\mathbf{A} \in \mathbb{C}^{I \times R}$, $\mathbf{B} \in \mathbb{C}^{J \times R}$ and $\mathbf{C} \in \mathbb{C}^{K \times R}$ have unit norm columns
with elements randomly drawn from a normal distribution, with $R = 3$.

Let $\mathcal{X} = \mathcal{X}_{0} + \alpha \mathcal{V}$ be a noisy version of $\mathcal{X}_{0}$, where $\mathcal{V}$ is the additive noise term, whose elements are drawn from a normal distribution. The parameter $\alpha$ controls the power (variance) of the noise term, and is defined as a function of the signal to noise ratio (SNR), in dB, as follows
$$\begin{equation} 
\text{SNR}_{\text{dB}} = 10\log_{10} \left( \frac{|| \mathcal{X}_{0} ||_{F}^{2}}{|| \alpha \mathcal{V} ||_{F}^{2}} \right) 
\end{equation}$$

Assuming the SNR range $\{0, 5, 10, 15, 20, 25, 30\}$ dB, find the estimates $\mathbf{\hat{A}}$, $\mathbf{\hat{B}}$ and $\mathbf{\hat{C}}$ obtained with the ALS algorithm for $(I, J, K) = (10, 4, 2)$.

Let us define the normalized mean square error (NMSE) measure

$$\begin{equation} 
\text{NMSE}(\mathbf{Q})= \frac{1}{M} \sum_{m=1}^{M}  \frac{|| \mathbf{\hat{Q}}(m) - \mathbf{Q}(m) ||_{F}^{2}}{|| \mathbf{Q}(m) ||_{F}^{2}}
\end{equation}$$


where $\mathbf{Q}(m)$ and $\mathbf{\hat{Q}}$ denote the original data matrix and the reconstructed one at the $m$-th Monte Carlo experiment, respectively. For each SNR value, plot NMSE(A), NMSE(B) and NMSE(C) as a function of the SNR. Discuss the obtained results.

<u>Hint</u>Hint: Don’t forget to take into account the inherent ambiguities of the CP decomposition.

---

### Results

SNR (dB) Mean NMSE (dB) for 1000 MC rounds


    SNR 	|	NMSE 	
    ----------------------------------
      X 
    -----------------------------------
    0	|	-7.71	
    5	|	-30.53	
    10	|	-50.40	
    15	|	-69.78	
    20	|	-88.74	
    25	|	-108.26	
    30	|	-127.41	
    -----------------------------------
      A 
    -----------------------------------
    0	|	12.32	
    5	|	13.83	
    10	|	15.03	
    15	|	15.17	
    20	|	15.34	
    25	|	15.23	
    30	|	15.31	
    -----------------------------------
      B 
    -----------------------------------
    0	|	1.94	
    5	|	1.00	
    10	|	0.82	
    15	|	0.98	
    20	|	0.95	
    25	|	0.89	
    30	|	0.87	
    -----------------------------------
      C 
    -----------------------------------
    0	|	13.06	
    5	|	14.69	
    10	|	15.18	
    15	|	15.33	
    20	|	15.97	
    25	|	15.48	
    30	|	15.63	
    -----------------------------------

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw11/code/figures/hw11-problem2.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<!--  -->
