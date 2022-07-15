<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 8 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [High Order Orthogonal Iteration (HOOI)](#high-order-orthogonal-iteration-hooi)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# High Order Orthogonal Iteration (HOOI)

## Problem 1
For a third-order tensor $\mathcal{X} \in \mathbb{C}^{I\times J\times K}$ implement the High Order Orthogonal Iteration (HOOI) method, using the following prototype function:

$$\begin{equation} 
\left[\mathcal{S}, \mathbf{U}^{(1)}, \mathbf{U}^{(2)}, \mathbf{U}^{(3)} \right] = \text{hooi}(\mathcal{X}) \end{equation}$$

Compare the results with the HOSVD algorithm.

<u>Hint</u>: Use the file "hooi_test.mat” to validate your results.

--- 

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm that combines the SVD with a iterative estimation was applied to the given initial tensor $\mathcal{X} \in \mathbb{C}^{I\times J\times K}$;
- $I, J, K = 3, 4, 5$;
- Number of iterations: $100$.

**Discussion**

To compare the given data with the estimated factors, we may use two main  experiment results: the orthogonality between the subtensors (slices) of $\mathcal{S}$ and the NMSE between the given data and obtained as output to HOOI.

	-------------- Number of Iterations -------------- 
	it = 100

The orthogonality assessment consists in compute the function $f_{ort}(\mathcal{S})$, that acumulates the scalar product bewteen the slices, following the equation below.

$$\begin{equation*}
f_{ort}(\mathcal{S}) = \sum_{k_{M}=1}^{K} \sum_{k_{N}=1}^{K} \text{vec}(\mathcal{S}_{:, :, k_{M}})^{\top}\text{vec}(\mathcal{S}_{:, :, k_{N}}) \quad for \, \, k_{M} \neq k_{N}
\end{equation*}$$

We obtain $f_{ort}(\mathcal{S}) = 0$, as expected for successful HOOI.

	-------------- The tensor slices are orthogonal if = zero -------------- 
	f_ort = 0.00 

The values obtained for NMSE present low SNR, with an emphasis to $\text{NMSE}(\mathcal{X, \hat{X}})$ value, with a very low SNR. 

	--------------- NMSE between a given tensor X and estimation -------------- 
	NMSE: -629.27 dB
	--------------- NMSE between a given tensor core S and its estimation ----- 
	NMSE: 6.80 dB
	--------------- NMSE between the factor matrices U and their estimation ----- 
	NMSE between U1 and its estimation: 8.52 dB
	NMSE between U2 and its estimation: -17.23 dB
	NMSE between U3 and its estimation: 10.02 dB

We can see that both results, Orthogonality and NMSE, support the proper algorithm estimation hypothesis.

We assess also the HOOI vs HOSVD perfomance with the NMSE (dB). The HOOI algorithm outperforms HOSVD in four estimations, and has equal perfomance in one. 

**NMSE (dB)**
| $(\mathcal{X, \hat{X}})$ | $(\mathcal{S, \hat{S}})$ | $(\mathbf{U}^{(1)}, \hat{\mathbf{U}}^{(1)})$ | $(\mathbf{U}^{(2)}, \hat{\mathbf{U}}^{(2)})$ | $(\mathbf{U}^{(3)}, \hat{\mathbf{U}}^{(3)})$ |
| :---: | :---: | :---: | :---: | :---: |
| +10.65 | +0.29 | 0 | + 23.25 | +5.9 |

[Problem 1 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

---

## Problem 2

Consider the two third-order tensors $\mathcal{X} \in \mathbb{C}^{8 \times 4 \times 10}$ and $\mathcal{Y} \in \mathbb{C}^{5 \times 5 \times 5}$ provided in the data file “hosvd_denoising.mat”. By using your HOOI prototype function, find a low multilinear rank approximation for these tensors, defined as $\tilde{\mathcal{X}} \in \mathbb{C}^{R1 \times R2 \times R3}$ and $\tilde{\mathcal{Y}} \in \mathbb{C}^{P1 \times P2 \times P3}$. Then, calculate the normalized mean square error (NMSE) between the original tensor and its approximation, i.e.,:

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

- The algorithm that uses the SVD was applied to the given initial tensor $\mathcal{X} \in \mathbb{C}^{R1 \times R2 \times R3}$ and $\mathcal{Y} \in \mathbb{C}^{P1 \times P2 \times P3}$;
- $R1, R2, R3 = 8, 4, 10$;
- $P1, P2, P3 = 5, 5, 5$.

**Discussion**

To compare the both random tensor estimation with given multilinear ranks, we may use NMSE results between the given data and obtained as output to HOOI. We may assess also by comparing the multilinear rank obtained in the tensor core $\mathcal{S}_{\mathcal{X}}$ and $\mathcal{S}_{\mathcal{Y}}$ estimated with the given ones.

	--------------- NMSE between a given tensor X and its estimation ----- 
	NMSE: -603.32 dB
	--------------- NMSE between a given tensor Y and its estimation ----- 
	NMSE: -619.16 dB

The values obtained for NMSE present very low SNR, less than $-600 \, \text{dB}$.

As defined in the proposed problem, the given ranks of $\mathcal{X}$ $\mathcal{Y}$ are $R1, R2, R3 = 8, 4, 10$, and $P1, P2, P3 = 5, 5, 5$, respectively.

	Tensor X multilinear rank: [8 4 10] 
	Tensor Y multilinear rank: [5 5 5]

We can see that that the algorithm provide the expected result, with the given ranks equal to the estimated. 

We assess also the HOOI vs HOSVD perfomance with the NMSE (dB). The HOOI algorithm outperforms HOSVD in four estimations, and has equal perfomance in one. 

**NMSE (dB)**
| $(\mathcal{X, \hat{X}})$ | $(\mathcal{Y, \hat{Y}})$ |
| :---: | :---: |
| +2.83 | +8.52 |

In conclusion, both results, NMSE and ranks estimation using the tensor core, support the proper algorithm estimation hypothesis.

[Problem 2 script][1].

</p>
</div>

<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw8/code/hw8.m> (Problem 1 script)