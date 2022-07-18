<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 1 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Problem 1](#problem-1)
- [Problem 2](#problem-2)
- [Problem 3](#problem-3)

# Problem 1

For randomly generated $\mathbf{A}$ and $\mathbf{B}$ $\in \mathbb{C}^{N\times N}$, create an algorithm to compute the Hadamard Product $\mathbf{A} \odot \mathbf{B}$. Then, compare the run time of your algorithm with the operator A.*B of the software Octave/Matlab $^{\textregistered}$. Plot the run time curve as a function of the number of rows/columns $N \in \{2, 4, 8, 16, 32, 64, 128\}$.

---
### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 500 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value, for $N = \{2,4,6,8,16,32,64, 128\}$.

**Discussion**

We can see that for all values of $N$, Matlab's method outperforms the Author's. For small values of $N$, the gap between them, $6 \times 10^{-5}$s vs $6 \times 10^{-6}$s, approximately ten times faster. However as the $N$ increases, that performance gap becomes more subtle.

[Problem 1 script][1]


</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/code/figures/hw1-problem1.png" alt="Hadamard Product Cost Figure" title="Hadamard Product Cost Figure" width="512" />
</p>

<!---------------------------------------------------------------------------->

# Problem 2 
For randomly generated $\mathbf{A}$ and $\mathbf{B}$ $\in \mathbb{C}^{N\times N}$ , create an algorithm to compute the Kronecker Product $\mathbf{A} \otimes \mathbf{B}$. Then, compare the run time of your algorithm with the operator kron(A, B) of the software Octave/Matlab $^{\textregistered}$. Plot the run time curve as a function of the number of rows/columns $N \in \{2, 4, 8, 16, 32, 64, 128\}$.

---
### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 500 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value, for $N = \{2,4,6,8,16,32,64, 128\}$.

**Discussion**

We can see that for all values of $N$, Matlab's method outperforms the author's. There's a narrow performance gap between them, up to three times faster. The difference varies very little regardless the value of $N$ increase.

[Problem 2 script][2]

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/code/figures/hw1-problem2.png" alt="Hadamard Product Cost Figure" title="Hadamard Product Cost Figure" width="512" />
</p>

<!---------------------------------------------------------------------------->

# Problem 3
For randomly generated $\mathbf{A}$ and $\mathbf{B}$ $\in \mathbb{C}^{N\times N}$ , create an algorithm to compute the Khatri-Rao Product $\mathbf{A} \diamond \mathbf{B}$ according with the following prototype function: 
$$\begin{equation*} R = kr(A, B). \end{equation*}$$

---
### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- 500 Monte Carlo Runs;
- Each Monte Carlo iteration uses a new matrix initialization from a Normal distribution $\mathcal{N}(0,\,1)\,$;
- Compute the mean for each value, for $N = \{2,4,6,8,16,32,64, 128\}$.

**Discussion**

The method developed by the author present similar behavior to Kronnecker product and a predictable trend for all values of $N$. 

[Problem 3 script][3]

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/code/figures/hw1-problem3.png" alt="Hadamard Product Cost Figure" title="Hadamard Product Cost Figure" width="512" />
</p>


<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/code/hw1_problem1.m> (Problem 1 script)
[2]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/code/hw1_problem2.m> (Problem 2 script)
[3]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/code/hw1_problem3.m> (Problem 3 script)