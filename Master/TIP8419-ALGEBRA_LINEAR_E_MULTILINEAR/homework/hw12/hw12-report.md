<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 12 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Tensor Kronecker Product Singular Value Decomposition (TKPSVD)](#tensor-kronecker-product-singular-value-decomposition-tkpsvd)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)

# Tensor Kronecker Product Singular Value Decomposition (TKPSVD)

## Problem 1

On a previous homework we have implemented the KPSVD (Kronecker Product Singular Value Decomposition) algorithm. Now, we will implement the generalization of that to tensors, namely, the TKPSVD (Tensor Kronecker Product Singular Value Decomposition) algorithm. Consider the $N$-order tensor 
$\mathcal{X} \in \mathbb{R}^{I_1 \times I_2 \times \dots \times I_N}$. Then, a TKPSVD of $\mathcal{X}$ can be written as

$$\begin{equation*} 
\mathcal{X} = \sum_{j=1}^{R} \sigma_{j} \mathcal{A}^{(d)} \otimes \mathcal{A}^{(d-1)} \otimes \dots \otimes \mathcal{A}^{(1)} \end{equation*}$$

where the tensors $\mathcal{A}^{(i)}_{j} \in \mathbb{R}^{I_1^{(i)} \times I_2^{(i)} \times \dots \times I_N^{(i)}}$ satisfy

$$\begin{equation*} 
|| \mathcal{A}^{(i)}_{j} ||_{F} = 1 \, , \prod_{d}^{i=1} I_k^{(i)} = I_K \, , 1 \leq k \leq N
\end{equation*}$$

Set $R = 1$ and generate $\mathcal{X} = \sigma \mathcal{A}^{(3)} \otimes \mathcal{A}^{(2)} \otimes \mathcal{A}^{(1)}$, for randomly chosen $\sigma$, $\mathcal{A}^{(1)} \in \mathbb{R}^{10 \times 4 \times 2}$, $\mathcal{A}^{(2)} \in \mathbb{R}^{5 \times 2 \times 2}$, and $\mathcal{A}^{(3)} \in \mathbb{R}^{2 \times 2 \times 2}$.

To avoid scaling ambiguity issues, normalize the model so that $||\mathcal{A}^{(i)} ||_{F} = 1$ for $i = 1, 2, 3$. (In other words, the norm will be absorbed by $\sigma$). Then, implement the TKPSVD algorithm that estimate $\mathcal{A}^{(1)}$, $\mathcal{A}^{(2)}$, $\mathcal{A}^{(3)}$ by solving the following problem

$$\begin{equation*} 
(\mathbf{\hat{A}^{(1)}}, \mathbf{\hat{A}^{(2)}}, \mathbf{\hat{A}^{(3)}} ) = \underset{ \mathcal{A}^{(1)}, \mathcal{A}^{(2)}, \mathcal{A}^{(3)} } \min \, || \mathcal{X} - \sigma \mathcal{A}^{(3)} \otimes \mathcal{A}^{(2)} \otimes \mathcal{A}^{(1)} ||_{F}^{2} \end{equation*}$$

Compare the estimated tensors $\mathbf{\hat{A}^{(1)}}, \mathbf{\hat{A}^{(2)}}, \mathbf{\hat{A}^{(3)}}$ with the original (normalized) ones. What can you conclude? Explain the results.

---

### Results

<!--  -->

## Problem 2

In a Monte Carlo experiment with $M = 1000$ realizations, generate 

$$\begin{equation*}
\mathcal{X}_{(0)} = \sum_{r=1}^{R} \sigma_{r} \mathcal{A}^{(3)} \otimes \mathcal{A}^{(2)} \otimes \mathcal{A}^{(1)} 
\end{equation*}$$ 

for $R = 2$ and randomly chosen $\mathcal{A}^{(1)} \in \mathbb{R}^{10 \times 4 \times 2}$, $\mathcal{A}^{(2)} \in \mathbb{R}^{5 \times 2 \times 2}$, and $\mathcal{A}^{(3)} \in \mathbb{R}^{2 \times 2 \times 2}$ whose elements are drawn from a standard normal distribution. As in the previous case, the scalars $\sigma_{r}$ are meant to absorb the scaling of each term, while the random tensors will have unit Frobenius norm (you should first draw these tensors and then normalize them). 

Let $\mathcal{X} = \mathcal{X}_{0} + \alpha \mathcal{V}$ be a noisy version of $\mathcal{X}_{0}$, where $\mathcal{V}$ is the additive noise term, whose elements are drawn from a normal distribution. The parameter $\alpha$ controls the power (variance) of the noise term, and is defined as a function of the signal to noise ratio (SNR), in dB, as follows
$$\begin{equation} 
\text{SNR}_{\text{dB}} = 10\log_{10} \left( \frac{|| \mathcal{X}_{0} ||_{F}^{2}}{|| \alpha \mathcal{V} ||_{F}^{2}} \right) 
\end{equation}$$

Assuming the SNR range $\{0, 5, 10, 15, 20, 25, 30\}$ dB, find the estimates $\mathbf{\hat{A}^{(1)}}$, $\mathbf{\hat{A}^{(2)}}$ and $\mathbf{\hat{A}^{(3)}}$ obtained with the ALS algorithm for $(I, J, K) = (10, 4, 2)$.

Let us define the normalized mean square error (NMSE) measure

$$\begin{equation} 
\text{NMSE}(\mathbf{Q})= \frac{1}{M} \sum_{m=1}^{M}  \frac{|| \mathbf{\hat{Q}}(m) - \mathbf{Q}(m) ||_{F}^{2}}{|| \mathbf{Q}(m) ||_{F}^{2}}
\end{equation}$$

where $\mathbf{Q}(m)$ and $\mathbf{\hat{Q}}$ denote the original data matrix and the reconstructed one at the $m$-th Monte Carlo experiment, respectively. For each SNR value, plot NMSE $(\mathbf{A}^{(j)})$ as a function of the SNR for $j= 1, 2, 3$, where the matrix $\mathbf{A}^{(j)}$ is defined as

$$\begin{equation} 
\mathbf{A}^{(j)} = \left[\text{vec}(\mathcal{A}_{1}^{(j)}) \text{vec}(\mathcal{A}_{2}^{(j)}) \right]
\end{equation}$$

Discuss the obtained results.

<u>Hint</u>Hint: Don’t forget to take into account the inherent ambiguities of the decomposition when computing the NMSE.

---

### Results

<!--  -->
