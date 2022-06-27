<div style="background-color:rgb(100, 255, 100, 0.25); text-align:center; padding:20px">
<p> 
Lista 2 [TI8419 - Multilinear Algebra]

Lucas Abdalah

Professors: André Lima e Henrique Goulart

</p> 
</div>

- - - 

# Table of Contents
- [Problem 1](#problem-1)
- [Problem 2](#problem-2)
- [Problem 3](#problem-3)
- [Problem 4](#problem-4)

# Problem 1

<div style="background-color:rgb(100, 100, 100, 0.15); text-align:left; padding:20px">
<p>

By using the properties of the outer product, show that 

$$\begin{equation*}\mathcal{X} = \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{2} \circ \mathbf{c}_{2} \end{equation*}$$

is a rank-one tensor whenever $\mathbf{b}_{1} = \mathbf{b}_{2}$ and $\mathbf{c}_{1} = \mathbf{c}_{2}$. Is this also true in general
when $\mathbf{c}_{1} = \mathbf{c}_{2}$ but $\mathbf{b}_{1} \neq \mathbf{b}_{2}$ ?

</p> 
</div>

---

1.1) O tensor $\mathcal{X}$ é construido a partir do seguinte modelo:
$$\begin{equation*}\mathcal{X} = \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{2} \circ \mathbf{c}_{2} \end{equation*}$$

Assumindo $\mathbf{b}_{1} = \mathbf{b}_{2}$ e $\mathbf{c}_{1} = \mathbf{c}_{2}$

$$\begin{align*}
\mathcal{X} &= \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1}\\
    &= (\mathbf{a}_{1} + \mathbf{a}_{2}) \circ \mathbf{b}_{1} \circ \mathbf{c}_{1}
\end{align*}$$

1.2) É conveniente observar a soma como um novo vetor: $\mathbf{v}_{1} = \mathbf{a}_{1} + \mathbf{a}_{2}$

1.3) Finalmente, esta representação torna mais evidente que apenas 1 tensor de posto-1 é suficiente para representar $\mathcal{X}$, i.e, tem posto $R = 1$.

$$\begin{equation*} \boxed{\mathcal{X} = \mathbf{v}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1}} \end{equation*}$$

---

1.4) Assumindo $\mathbf{b}_{1} \neq \mathbf{b}_{2}$ e $\mathbf{c}_{1} = \mathbf{c}_{2}$

$$\begin{align*}
\mathcal{X} &= \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{2} \circ \mathbf{c}_{1}\\
    &= (\mathbf{a}_{1} \circ \mathbf{b}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{2}) \circ \mathbf{c}_{1}
\end{align*}$$

1.5) Com estas premissas, apenas a reorganização não permite afirmar que $\mathcal{X}$ é representado por um tensor de posto 1.

1.6) Entretanto, ao assumir que existe $\alpha$ tal que $\mathbf{b}_{2} = \alpha \mathbf{b}_{1}$, a equação fica mais simples.

$$\begin{align*}
\mathcal{X} &= \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{2} \circ \alpha \mathbf{b}_{1} \circ \mathbf{c}_{1}\\
    &= \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \alpha \mathbf{a}_{2} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1}\\
    &= (\mathbf{a}_{1} + \alpha \mathbf{a}_{2}) \circ \mathbf{b}_{1} \circ \mathbf{c}_{1}
\end{align*}$$

É conveniente observar a soma como um novo vetor: $\mathbf{v}_{2} = \mathbf{a}_{1} + \alpha \mathbf{a}_{2}$

1.7) Finalmente, esta representação torna mais evidente que apenas 1 tensor de posto-1 é suficiente para representar $\mathcal{X}$, i.e, tem posto $R = 1$.

$$\begin{equation*} \mathcal{X} = \mathbf{v}_{2} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} \end{equation*}$$

1.8) Entretanto, se não existe $\alpha$, ou seja $\mathbf{b}_{1}$ e $\mathbf{b}_{2}$ não são colineares, consequentemente $\mathbf{v}_{2}$ não existe e apenas 1 tensor de posto-1 é insuficiente para representar $\mathcal{X}$. Sendo utilizados 2 tensores de posto-1 para representação, implica em $posto(\mathcal{X}) = 2$.

$$\begin{equation*} \boxed{\mathcal{X} = (\mathbf{a}_{1} \circ \mathbf{b}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{2}) \circ \mathbf{c}_{1}} \end{equation*}$$

<!---------------------------------------------------------------------------->

# Problem 2

<div style="background-color:rgb(100, 100, 100, 0.15); text-align:left; padding:20px">
<p>

Show that the tensor rank is indeed a tensor property: in other words, it is invariant with respect to a multilinear transformation by nonsingular matrices, that is, if

$$\begin{equation*}\mathcal{X} = \mathcal{S} \times_{1} \mathbf{A}^{(1)} \dots \times_{N} \mathbf{A}^{(N)}  \end{equation*}$$

where $\mathbf{A}^{(n)} \in \mathbb{C}^{I_{n} \times I_{n}}$ is nonsingular for every $n$, then 

$$\begin{equation*} rank(\mathcal{X}) = rank(\mathcal{S}). \end{equation*}$$

(Hint: write $\mathcal{S}$ as a PD with a minimal number of terms, and then use the properties of the multilinear transformation to bound the rank of $\mathcal{X}$ ; similarly, use the invertibility of the multilinear transformation to bound the rank of $\mathcal{S}$. More generally, conclude that the same property holds for matrices $\mathbf{A}^{(n)} \in {\mathbb{C}}^{I_{n} \times R_{n}}$ having linearly independent columns (and thus $R_{n} \leq I_{n}$).

</p> 
</div>

---

2.1) O tensor *core* $\mathcal{S}$ reescrito em função de fatores da decomposição CP resulta em:

$$\begin{equation*} \mathcal{S} = \sum_{r=1}^{R} \mathbf{s}_{r}^{(1)} \circ \dots \circ \mathbf{s}_{r}^{(N)} \end{equation*}$$

2.2) Também é possível reorganizar a equação em função $\mathcal{S}$ de acordo com as equações (11) e (17) das notas de aula, utilizando as propriedades do operador transposto $\left(^{\top}\right)$, caso $\mathbf{A}^{(n)} \in \mathbb{R}$, e operador hermitiano/autoadjunto $\left(^{H}\right)$

---

$$\begin{equation*}\mathcal{S} = \mathcal{X} \times_{1} {\mathbf{A}^{(1)}}^{H} \dots \times_{N} {\mathbf{A}^{(N)}}^H \end{equation*}$$


$$\begin{equation*}\mathcal{X} = \mathcal{S} \times_{1} \mathbf{A}^{(1)} \dots \times_{N} \mathbf{A}^{(N)} \end{equation*}$$




<!---------------------------------------------------------------------------->

# Problem 3

Let X P CI1ˆI2ˆI3 be given by
X “ a1  ̋ b1  ̋ c1 ` a2  ̋ b2  ̋ c1 ` a1  ̋ b2  ̋ c2, (1)
where the vectors are assumed to satisfy the following:
• a1 is not collinear with a2;
• b1 is not collinear with b2;
• c1 is not collinear with c2.
The goal of this exercise is to show that any such tensor has rank three, that is,
it cannot be expressed as a sum of fewer terms. We will proceed by steps.
1
(i) First, show that
X “ S ˆ1 A ˆ2 B ˆ3 C,
where
A “ “a1 a2
‰ , B “ “b1 b2
‰ , C “ “c1 c2
‰ ,
and
S ̈ ̈1 “ I2 “
„1 0
0 1

, S ̈ ̈2 “
„0 1
0 0

.
Then, using the result of Exercise 2), conclude that X and S have the
same rank.
(ii) Hence, it suffices to show that rankpSq “ 3. Suppose, for a contradiction,
that rankpSq “ 2. Using the properties of the PARAFAC decomposition,
show that this imples the existence of matrices U, V, D1, D2 P C2ˆ2 such
that D1, D2 are diagonal and
S ̈ ̈1 “ UD1VT, S ̈ ̈2 “ UD2VT. (2)
(iii) Now, use the fact that X ̈ ̈1 “ I to show that (2) implies that S ̈ ̈2 can be
diagonalized by U, that is, there exists a diagonal matrix D P C2ˆ2 such
that
S “ UDU ́1.
(iv) Conclude that this leads to a contradiction, by taking into account the
Jordan form of S ̈ ̈2.

<!---------------------------------------------------------------------------->


# Problem 4

n this last exercise, we will show that, although the tensors of the form
considered in the last exercise have rank 3, they are limits of sequences of rank-
2 tensors. Thus, unlike happens for matrices, a sequence of rank-R tensors can
converge to a rank-S tensor with S ą R.
(i) First, show that the rank-1 tensor
Ym “ mpa1 ` m ́1b2q  ̋ pb2 ` m ́1b1q  ̋ pc1 ` m ́1c2q
is equal to X (as given by (1)) plus an Opmq term Zm and an Op1{mq
term.
(ii) Subtract the Opmq term to get:
Xm “ Ym  ́ Zm.
What is the rank of Xm?
(iii) Use the expression obtained for Xm to conclude that
lim
mÑ8 Xm “ X .