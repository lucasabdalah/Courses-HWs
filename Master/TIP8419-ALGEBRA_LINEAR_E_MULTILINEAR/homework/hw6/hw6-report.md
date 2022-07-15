<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Homework 6 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Unfolding, folding, and $n$-mode product](#unfolding-folding-and-n-mode-product)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)
  - [Problem 3](#problem-3)

# Unfolding, folding, and $n$-mode product

<!-------------------------------------------------------------------------- -->

## Problem 1

For a third-order tensor $\mathbf{X} \in \mathbb{C}^{I \times J \times K}$, using the concept of $n$-mode fibers, implement the function unfold according to the following prototype 
$$\begin{equation*}
[\mathcal{X}]_{(n)} = \text{unfold}(\mathcal{X}, n)
\end{equation*}$$

<u>Hint</u>: Use the file “unfolding_folding.mat” to validate your function.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm was applied to reshape the original data into a $N$-mode tensor;
- $N$ in range $\{1, 2, 3\}$.

**Discussion**

- Experiment proposed in the example 2.6 of the book Multi-way Analysis With Applications in the Chemical Sciences (Smilde, 2004).

Tensor X

	X(:, :, 1)
	1  2  3;
	4  5  6;
	7  8  9;
	3  2  1;

	X(:, :, 2)
	5  6  7;
	8  9  4;
	5  3  2;
	4  5  6;

Tensor X (mode-1)

	X(4, 6)
	1  2  3  5  6  7;
	4  5  6  8  9  4;
	7  8  9  5  3  2;
	3  2  1  4  5  6;

Tensor X (mode-2)

	X(3, 8)
	1  4  7  3  5  8  5  4;
	2  5  8  2  6  9  3  5;
	3  6  9  1  7  4  2  6;

Tensor X (mode-3)

	X(2, 12)
	1  4  7  3  2  5  8  2  3  6  9  1;
	5  8  5  4  6  9  3  5  7  4  2  6;

- Validation

 Unfold difference

	sum(X1 - unfold(X, 1)) = 0.00 
	sum(X2 - unfold(X, 2)) = 0.00 
	sum(X3 - unfold(X, 3)) = 0.00 

We assess the difference between the given data and the algorithm output, and we can see that the residuals sum leads to zero.

[Problem script][1].

Fold experiment output log: [Unfold Txt File][2].

</p>
</div>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw6/ex1-img1.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw6/ex1-img2.png" alt="Khatri-Rao Product Cost Figure" title="Khatri-Rao Product Cost Figure" width="512" />
</p>

<!-------------------------------------------------------------------------- -->

## Problem 2

Implement the function fold that converts the unfolding $[\mathcal{X}]_{(n)}$ obtained with $\text{unfold}(\mathcal{X}, n)$ back to the tensor $\mathcal{X} \in \mathbb{C}^{I \times J \times K}$ (i.e., a 3-d array in Matlab/Octave), according to the following prototype:
$$\begin{equation*}
\mathcal{X} = \text{fold}([\mathcal{X}]_{(n)}, [I J K],  n)
\end{equation*}$$

<u>Hint</u>: Use the file “unfolding_folding.mat” to validate your function.

--- 

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm was applied to build a tensor from a $N$-mode tensor;
- $N$ in range $\{1, 2, 3\}$.

**Discussion**

- Experiment proposed in the example 2.6 of the book Multi-way Analysis With Applications in the Chemical Sciences (Smilde, 2004).

Tensor X (mode-1)

	X(4, 6)
	1  2  3  5  6  7;
	4  5  6  8  9  4;
	7  8  9  5  3  2;
	3  2  1  4  5  6;


Tensor X (mode-2)

	X(3, 8)
	1  4  7  3  5  8  5  4;
	2  5  8  2  6  9  3  5;
	3  6  9  1  7  4  2  6;


Tensor X (mode-3)

	X(2, 12)
	1  4  7  3  2  5  8  2  3  6  9  1;
	5  8  5  4  6  9  3  5  7  4  2  6;


Tensor X from (mode-1)

	X(:, :, 1)
	1  2  3;
	4  5  6;
	7  8  9;
	3  2  1;

	X(:, :, 2)
	5  6  7;
	8  9  4;
	5  3  2;
	4  5  6;


Tensor X from (mode-2)

	X(:, :, 1)
	1  2  3;
	4  5  6;
	7  8  9;
	3  2  1;

	X(:, :, 2)
	5  6  7;
	8  9  4;
	5  3  2;
	4  5  6;


Tensor X from (mode-3)

	X(:, :, 1)
	1  2  3;
	4  5  6;
	7  8  9;
	3  2  1;

	X(:, :, 2)
	5  6  7;
	8  9  4;
	5  3  2;
	4  5  6;

- Validation

 Fold difference

	sum(tenX - fold(X1)) = 0.00 
	sum(tenX - fold(X2)) = -0.00 
	sum(tenX - fold(X3)) = -0.00 


We assess the difference between the given data and the algorithm output, and we can see that the residuals sum leads to zero.

[Problem script][1].

Fold experiment output log: [Fold Txt File][3].

</p>
</div>

<!-------------------------------------------------------------------------- -->

## Problem 3

For given matrices $\mathbf{A} \in \mathbb{C}^{P \times I}$, $\mathbf{B} \in \mathbb{C}^{Q \times J}$, $\mathbf{C} \in \mathbb{C}^{R \times K}$
and tensor $\mathcal{X} \in \mathbb{C}^{I \times J \times K}$, calculate the tensor $\mathcal{Y} \in \mathbb{C}^{P \times Q \times R}$ via the following multilinear transformation:
$$\begin{equation*}
\mathcal{Y} = \mathcal{X} \times_{1} \mathbf{A} \times_{2} \mathbf{B} \times_{3} \mathbf{C}
\end{equation*}$$

<u>Hint</u>: Use the file “multilinear_product.mat” to validate your result.

---

### Results

<div style="background-color:rgba(0, 0, 200, 0.15); text-align:justify; padding:20px">
<p>

**Simulation setup**

- The algorithm was applied to compute the N-mode product between a given tensor and factor matrices.

**Discussion**

The results are consistent with the proposed scenario, since given data after the algorithm succeds to obtain a very low NMSE (dB) value.
	
	NMSE between a given tensor and its version afected by the N-mode product: -666.47 dB

[Problem script][1]

</p>
</div>

<!---------------------------------------------------------------------------->

[1]: <https://github.com/lucasabdalah/Courses-HWs/blob/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw6/code/hw6.m> (Problem script)
[2]: <https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw6/code/hw6_unfold.txt> (Unfold Txt File)
[3]: <https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw6/code/hw6_fold.txt> (Fold Txt File)
