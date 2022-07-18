<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 2 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Problem 1](#problem-1)
- [Problem 2](#problem-2)

# Problem 1

Generate $\mathbf{X} = \mathbf{A} \diamond \mathbf{B} \in \mathbb{C}^{I \times R}$, for randomly chosen $\mathbf{A} \in \mathbb{C}^{I \times R}$ and $\mathbf{B} \in \mathbb{C}^{I \times R}$. Compute the left pseudo-inverse of $\mathbf{X}$ and obtain a graph that shows the run time vs. number of rows $(I)$ for the following methods.

**Method 1:**

Matlab/Octave function: $pinv(\mathbf{X}) = pinv(\mathbf{A} \diamond \mathbf{B})$

**Method 2:**

$\mathbf{X}^{\dagger} = (\mathbf{X}^{\top} \mathbf{X})^{-1} \mathbf{X}^{\top} = [(\mathbf{A} \diamond \mathbf{B})^{\top} (\mathbf{A} \diamond \mathbf{B})]^{-1} (\mathbf{A} \diamond \mathbf{B})^{\top}$

**Method 3:**

$\mathbf{X}^{\dagger} = [(\mathbf{A}^{\top} \mathbf{A}) \odot (\mathbf{B}^{\top} \mathbf{B})]^{-1} (\mathbf{A} \diamond \mathbf{B})^{\top}$

<u>Note</u>: Consider the range of values $I \in \{2, 4, 8, 16, 32, 64, 128, 256\}$ and plot the curves for $R = 2$ and $R = 4$.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 500 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value, for $N = \{2,4,6,8,16,32,64, 128, 256\}$.


**Discussion**

We can see that for all values of $I$, Matlab's method is outperformed by the methods 2 and 3. All methods present a subtle gap between their cost, approximately constant. Method 2 is two times faster then Matlab, while method 3 is ten times faster.

The experiment with $R = 4$ also supports the results presented for $R = 2$, with very similar plots.

[Problem 1 script][1] and [Figures][3].

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/code/figures/hw2-problem1a.png" alt="Left Pseudo Inverse Cost Figure (R=2)" title="Left Pseudo Inverse Cost Figure (R=2)" width="512" />
</p>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/code/figures/hw2-problem1b.png" alt="Left Pseudo Inverse Cost Figure (R=4)" title="Left Pseudo Inverse Cost Figure (R=4)" width="512" />
</p>

<!---------------------------------------------------------------------------->

# Problem 2 

Generate 
$\overset{N}{\underset{n=1}\diamond} \mathbf{A}_{(n)} = \mathbf{A}_{(1)} \diamond \dots \diamond \mathbf{A}_{(N)}$, where every $\mathbf{A}_{(n)}$ has dimensions $4 \times 2$, $n = 1, \dots , N$. Evaluate the run time associated with the computation of the Khatri-Rao product as a function of the number $N$ of matrices for the above methods.

<u>Note</u>: Consider the range of values $N \in \{2, 4, 6, 8, 10\}$. 

The symbols $\odot$ and $\diamond$  denotes the Hadamard and the Khatri-Rao Product, respectively.

---
### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 500 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Each matrix has $4 \times 2$ dimension;
- Compute the mean for each value, for $N = \{2,4,6,8,10\}$.

**Discussion**

The results are consistent with the experiment perfomed in [HW1](https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/hw1-report.pdf), that for randomly generated $\mathbf{A}$ and $\mathbf{B}$ $\in \mathbb{C}^{N\times N}$, an algorithm to compute the Khatri-Rao Product $\mathbf{A} \diamond \mathbf{B}$ was created according with the following prototype function: 
$$\begin{equation*} R = kr(A, B). \end{equation*}$$


[Problem 2 script][2]

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/code/figures/hw2-problem2.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>


<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/code/hw2_problem1.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/code/hw2_problem2.m> (Problem 2 script)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/code/hw2.m> (Figures)