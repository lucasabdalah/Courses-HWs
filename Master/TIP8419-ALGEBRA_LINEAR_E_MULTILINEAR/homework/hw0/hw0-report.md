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
- Compute the mean for each value, for $N = 2,4,6,8,16,32,64$.

**Discussion**

We can see that for all values of $N$, the second method outperforms the first. For small values of $N$, the difference is more subtle, ten times faster. However as the $N$ increases, the performance gap increases up to a hundred times faster.

[Problem 1.a script][1]

</p>
</div>


<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw0/code/hw0-problem1-a.png" alt="Kronecker Product and Matrix Inversion Cost Figure" title="Kronecker Product and Matrix Inversion Cost" width="512" />
</p>

- - -

## (b)

**Method 1:**

$$\begin{equation*}  
\left(\mathbf{A}_{N \times N}^{(1)} \otimes \mathbf{A}_{N \times N}^{(2)} \otimes \mathbf{A}_{N \times N}^{(3)} \otimes \dots \otimes \mathbf{A}_{N \times N}^{(K)} \right)^{-1} = \left(\overset{K}{\underset{k=1}\otimes} \mathbf{A}_{N \times N}^{(k)} \right) ^{-1}
\end{equation*}$$

**Method 2:**

$$\begin{equation*}  
\left(\mathbf{A}_{N \times N}^{(1)}\right)^{-1}  \otimes \left(\mathbf{A}_{N \times N}^{(2)}\right)^{-1}  \otimes \left(\mathbf{A}_{N \times N}^{(3)}\right)^{-1}  \otimes \dots \otimes \left(\mathbf{A}_{N \times N}^{(K)}\right)^{-1} = \overset{K}{\underset{k=1}\otimes} \left(\mathbf{A}_{N \times N}^{(k)} \right) ^{-1}
\end{equation*}$$

For $k \in \{2,4,6,8,10\}$.

--- 

### Results


<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 200 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value for $N = 2$ and $K = 2,4,6,8,10$.

**Discussion**

On the scenario proposed, with $N = 4$, the amount of memory (ram) is up to greater than 64.0 Gb. Since a single complex element requires 16 bytes, the simulation using the homework setup fails from K = 8, since it's required more RAM memory than the available, 19.8 Gb. This value consider 100% of ram use, without taking into count the operational system (OS), backgroud scripts or matlab.

--- 

**Example**

To illustrate, the function [kron_dim][3] may be applied for the example with $N = 4$ $k = 7$:
```
Matrix Dimensions: 16384X16384 
N of elements: 268435456 
Memory use:  4 Gb
```
Since each matrix is $4 \times 4$, each Kronnecker product multiplies by 16 the amount of RAM required, hence the matrix product with $K = 8$ leads it to an error.

<div style="background-color:rgba(250, 0, 0, 0.25); text-align:left; padding:20px">
<p> 

``` 
Requested 4x16384x4x16384 (64.0GB) array exceeds maximum array size preference (19.8GB). This might cause MATLAB to become unresponsive. 
``` 
</p>
</div>

---

Finally, we set $N = 2$ for maximum usage when $K = 10$, since it leads to a $2^{10} \times 2^{10}$ matrix, with 1048576 elements and only 16 Mb of ram use. 
```
Matrix Dimensions: 1024X1024 
N of elements: 1048576 
Memory use: 16 Mb
```

We can see that for all values of $K$, the second method outperforms the first.
Both results support the hypothesis that the inversion of smaller matrices in Matlab is much more effective.

[Problem 1.b script][2]

</p>
</div>


<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw0/code/hw0-problem1-b.png" alt="Kronecker Product and Matrix Inversion Cost Figure" title="Kronecker Product and Matrix Inversion Cost" width="512" />
</p>

- - -

# Problem 2
Let $\text{eig}(\mathbf{X})$ be the function that returns the matrix $\sum_{K \times K}$ of eigenvalues of $\mathbf{X}$. Show algebraically that $\text{eig}(\mathbf{A} \otimes \mathbf{B}) = \text{eig}(\mathbf{A}) \otimes \text{eig}(\mathbf{B})$.

<u>Hint</u>: Use the property $(\mathbf{A} \otimes \mathbf{B})(\mathbf{C} \otimes \mathbf{D})= \mathbf{A}\mathbf{C} \otimes \mathbf{B}\mathbf{D}$ 

---

We write the SVD for each matrix, $\mathbf{A}$ and $\mathbf{B}$, as:

$$\begin{align*}
    \mathbf{A} &= \mathbf{U}_A \mathbf{\Sigma}_A \mathbf{V}_{A}^{H}\\
    \mathbf{B} &= \mathbf{U}_B \mathbf{\Sigma}_B \mathbf{V}_{B}^{H},
\end{align*}$$

We take advantage of the definitions to the equation 
$\text{eig}(\mathbf{A}\otimes\mathbf{B})$ and using two times the property suggested by the exercise, we have:

$$\begin{align*}
    \text{eig}\left(\mathbf{U}_A \mathbf{\Sigma}_A \mathbf{V}_{A}^{H}\otimes 
    \mathbf{U}_B \mathbf{\Sigma}_B \mathbf{V}_{B}^{H}\right) &=
    \text{eig}\left[(\mathbf{U}_A \otimes \mathbf{U}_B) ( \mathbf{\Sigma}_A \mathbf{V}_{A}^{H} \otimes \mathbf{\Sigma}_B \mathbf{V}_{B}^{H}) \right] \\
    
        &= \text{eig}\left[(\mathbf{U}_A \otimes \mathbf{U}_B \big)
        \big( \mathbf{\Sigma}_A \otimes \mathbf{\Sigma}_B \big)
        \big( \mathbf{V}_{A} \otimes \mathbf{V}_{B}\big)^{H}\right]\\

        &= \mathbf{\Sigma}_A \otimes \mathbf{\Sigma}_B = 
        \text{eig}\big( \mathbf{A}\big) \otimes \text{eig}\big( \mathbf{B} \big),
\end{align*}$$

 by applying the operator $\text{eig}(\cdot)$ that returns the 
eigenvalue matrix $\mathbf{\Sigma}_A \otimes \mathbf{\Sigma}_B$. 


--- 

<!-- References -->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/c185d153949c2784ac8e6e173d775dca0b3fef04/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw0/code/hw0_problem1_a.m#L4> (Problem 1.a script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/c185d153949c2784ac8e6e173d775dca0b3fef04/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw0/code/hw0_problem1_b.m#L4> (Problem 1.b script)
[3]: <https://github.com/lucasabdalah/basic-functions/blob/e8d92b3750a4041588375a4cc97c7636b9392165/functions/ram_use.m#L45> (kron_dim)

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