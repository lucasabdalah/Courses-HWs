<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 0 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
1. [Problem 1](#problem-1)
2. [Problem 2](#problem-2)

# Problem 1

For randomly generated $\mathbf{A} \in \mathbb{C}^{N\times N}$ and $\mathbf{B} \in \mathbb{C}^{N\times N}$, evaluate the computational performance (run time) of the following matrix inversion formulas:

## (a)

**Method 1:**

$$\begin{equation*} 
(\mathbf{A}_{N \times N} \otimes \mathbf{B}_{N \times N} )^{-1}
\end{equation*}$$

**Method 2:**

$$\begin{equation*} 
(\mathbf{A}_{N \times N})^{-1} \otimes (\mathbf{B}_{N \times N})^{-1}
\end{equation*}$$

For $n \in \{2,4,6,8,16,32,64\}$.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 100 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value, for $N = 2,4,6,8,16,32,64$ ;
- Noiseless Scenario.

**Discussion**

We can see that for all values of $N$, the second method outperforms the first. For small values of $N$, the difference is more subtle, ten times faster. However as the $N$ increases, the performance gap increases up to a hundred times faster.

</p>
</div>


<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/homework_0/hw0-problem1-a.png" alt="Kronecker Product and Matrix Inversion Cost Figure" title="Kronecker Product and Matrix Inversion Cost" width="512" />
</p>

- - -

## (b)

**Method 1:**

$$\begin{equation*}  
\left(\mathbf{A}_{4 \times 4}^{(1)} \otimes \mathbf{A}_{4 \times 4}^{(2)} \otimes \mathbf{A}_{4 \times 4}^{(3)} \otimes \dots \otimes \mathbf{A}_{4 \times 4}^{(K)} \right)^{-1} = \left(\overset{K}{\underset{k=1}\otimes} \mathbf{A}_{4 \times 4}^{(k)} \right) ^{-1}
\end{equation*}$$

**Method 2:**

$$\begin{equation*}  
\left(\mathbf{A}_{4 \times 4}^{(1)}\right)^{-1}  \otimes \left(\mathbf{A}_{4 \times 4}^{(2)}\right)^{-1}  \otimes \left(\mathbf{A}_{4 \times 4}^{(3)}\right)^{-1}  \otimes \dots \otimes \left(\mathbf{A}_{4 \times 4}^{(K)}\right)^{-1} = \overset{K}{\underset{k=1}\otimes} \left(\mathbf{A}_{4 \times 4}^{(k)} \right) ^{-1}
\end{equation*}$$

For $k \in \{2,4,6,8,10\}$.

--- 

### Results


<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Discussion about the parameters:**

On the proposed scenario the 

The script to peform Problem B with matrices NxN and kron productory until K, with N = 4, and K = 8, results in an dimensional error, since it's required more the RAM memory available.

For $N = 4$ $k = 7$

Matrix Dimensions: 16384X16384 
N of elements: 268435456 
Memory use:  4 Gb

<div style="background-color:rgba(250, 0, 0, 0.25); text-align:left; padding:20px">
<p> 

Problem For $k = 8$ with $N = 4$

``` 
Requested 4x16384x4x16384 (64.0GB) array exceeds maximum array size preference (19.8GB). This might cause MATLAB to become unresponsive. 
``` 
</p>
</div>


**Simulation setup**

- 100 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value for $N = 2$ and $K = 2,4,6,8,10$ ;
- Noiseless Scenario.

**Discussion**

</p>
</div>





- - - 

- - -

# Problem 2
Let $eig(\mathbf{X})$ be the function that returns the matrix $\sum_{K \times K}$ of eigenvalues of $\mathbf{X}$. Show algebraically that $eig(\mathbf{A} \otimes \mathbf{B}) = eig(\mathbf{A}) \otimes eig(\mathbf{B})$.

<u>Hint</u>: Use the property $(\mathbf{A} \otimes \mathbf{B})(\mathbf{B} \otimes \mathbf{D})= \mathbf{A}\mathbf{C} \otimes \mathbf{B}\mathbf{D}$ 

[Hobbit lifestyles][1]


<!-- References -->

[1]: <https://en.wikipedia.org/wiki/Hobbit#Lifestyle> (Hobbit lifestyles)

<!-- See the section on [`Problem 1`](#Method1) -->


<!-- <span style="color:red">some **blue** text</span>. -->

<!-- <div style="background-color:rgba(250, 0, 0, 0.25); text-align:center; padding:20px">
<p> Alo teste </p>
</div>

<div style="background-color:rgba(0, 0, 250, 0.25); text-align:center; padding:20px">
<p> Alo teste </p>
</div>

<div style="background-color:rgba(0, 250, 0, 0.25); text-align:center; padding:20px">
<p> Alo teste </p>
</div> -->

<!-- <div style="background-color:rgba(0, 0, 0, 0.0470588); text-align:center; vertical-align: middle; padding:40px 0; margin-top:30px">
<a href="/blog">VIEW THE BLOG</a>
</div> -->