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

$$\begin{equation*}\mathcal{S} = \mathcal{X} \times_{1} {\mathbf{A}^{(1)}}^{H} \dots \times_{N} {\mathbf{A}^{(N)}}^H \end{equation*}$$

2.3) Desenvolvendo $\mathcal{X}$ em função de 2.1, obtém-se a representação em função do somatório ponderado pelas matrizes de mudança de base $\mathbf{A}^{(n)}$.

$$\begin{equation*}\mathcal{X} = \left(\sum_{r=1}^{R} \mathbf{s}_{r}^{(1)} \circ \dots \circ \mathbf{s}_{r}^{(N)} \right) \times_{1} \mathbf{A}^{(1)} \dots \times_{N} \mathbf{A}^{(N)} \end{equation*}$$

$$\begin{equation*}\mathcal{X} = \sum_{r=1}^{R} \mathbf{A}^{(1)} \mathbf{s}_{r}^{(1)} \circ \dots \circ \mathbf{A}^{(N)} \mathbf{s}_{r}^{(N)} \end{equation*}$$

2.4) Dado as considerações anteriores para $\mathcal{S}$, o tensor *core* $\mathcal{S}$ também pode ser reescrito como:

$$\begin{equation*}\mathcal{S} = \left(\sum_{r=1}^{R} \mathbf{A}^{(1)} \mathbf{s}_{r}^{(1)} \circ \dots \circ \mathbf{A}^{(N)} \mathbf{s}_{r}^{(N)} \right) \times_{1} {\mathbf{A}^{(1)}}^{H} \dots \times_{N} {\mathbf{A}^{(N)}}^H \end{equation*}$$

2.5) Dado as propriedas que relacionam os produtos de modo-$N$, o tensor $\mathcal{S}$, convenientemente é reorganizado de modo que as matrizes hermitianas multiplicam a matriz de transformação original de mesmo modo.

$$\begin{equation*}\mathcal{S} = \sum_{r=1}^{R} {\mathbf{A}^{(1)}}^{H} \mathbf{A}^{(1)} \mathbf{s}_{r}^{(1)} \circ \dots \circ {\mathbf{A}^{(N)}}^{H} \mathbf{A}^{(N)} \mathbf{s}_{r}^{(N)} \end{equation*}$$

2.6) Sendo ${\mathbf{A}^{(n)}}^{H} \mathbf{A}^{(n)} = \mathbf{I}$, obtemos:

$$\begin{equation*} \mathcal{S} = \sum_{r=1}^{R} \mathbf{s}_{r}^{(1)} \circ \dots \circ \mathbf{s}_{r}^{(N)} \end{equation*}$$

2.7) Consequentemente, $\mathcal{S}$ preserva o posto de $\mathcal{X}$, i.e, $\text{rank}(\mathcal{X}) = \text{rank}(\mathcal{S})$.

---

2.8) Para o primiro caso onde a matriz era quadrada, os operadores de inversa e hermitiano foram utilizdos. Já para o caso $\mathbf{A}^{(n)} \in {\mathbb{C}}^{I_{n} \times R_{n}}$, é necessário estender o conceito para pseudo-inversa de uma matriz. Assumindo a sua existência, obtém-se

$$\begin{equation*}{\mathbf{A}^{(n)}}^{\dagger} \mathbf{A}^{(n)} = \mathbf{I}\end{equation*}$$

E isso permite a extensão da demonstração para matrizes retangulares.

<!---------------------------------------------------------------------------->

# Problem 3

<div style="background-color:rgb(100, 100, 100, 0.15); text-align:left; padding:20px">
<p>

Let $\mathcal{X} \in {\mathbb{C}}^{I_{1} \times I_{2} \times I_{3}}$ be given by

$$\begin{equation*}\mathcal{X} = \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{2} \circ \mathbf{b}_{2} \circ \mathbf{c}_{1} + \mathbf{a}_{1} \circ \mathbf{b}_{2} \circ \mathbf{c}_{2} \end{equation*}$$

where the vectors are assumed to satisfy the following:
- a1 is not collinear with a2;
- b1 is not collinear with b2;
- c1 is not collinear with c2.

The goal of this exercise is to show that any such tensor has rank three, that is,
it cannot be expressed as a sum of fewer terms. We will proceed by steps.

(i) First, show that 

$$\begin{equation*}\mathcal{X} = \mathcal{S} \times_{1} \mathbf{A} \times_{2} \mathbf{B} \times_{3} \mathbf{C} \end{equation*}$$

where

$$\begin{equation*} \mathbf{A} = [\mathbf{a}_{1} \quad \mathbf{a}_{2}], \quad \mathbf{B} = [\mathbf{b}_{1} \quad \mathbf{b}_{2}], \quad \mathbf{C} = [\mathbf{c}_{1} \quad \mathbf{c}_{2}], \end{equation*}$$

and

$$\begin{equation*} \mathbf{S}_{..1} = \mathbf{I}_{2} = 
\begin{bmatrix}
1 & 0 \\
0 & 1 \\
\end{bmatrix}, \quad \mathbf{S}_{..2} = 
\begin{bmatrix}
0 & 1 \\
0 & 0 \\
\end{bmatrix}.
\end{equation*}$$

Then, using the result of Exercise 2), conclude that $\mathcal{X}$ and $\mathcal{S}$ have the
same rank.

(ii) Hence, it suffices to show that $\text{rank}(\mathcal{S}) = 3$. Suppose, for a contradiction, that $\text{rank}(\mathcal{S}) = 2$. Using the properties of the PARAFAC decomposition, show that this imples the existence of matrices $\mathbf{U}, \mathbf{V}, \mathbf{D}_{1}, \mathbf{D}_{2} \in \mathbb{C}^{2 \times 2}$ such
that $\mathbf{D}_{1}$, $\mathbf{D}_{2}$ are diagonal and

$$\begin{equation*} \mathbf{S}_{..1} = \mathbf{U} \mathbf{D}_{1} \mathbf{V}^{\top}, \quad \mathbf{S}_{..2} =  \mathbf{U} \mathbf{D}_{2} \mathbf{V}^{\top}
\end{equation*}$$

(iii) Now, use the fact that $\mathbf{X}_{..1} = \mathbf{I}$ to show that (2) implies that $\mathbf{S}_{..2}$ can be diagonalized by $\mathbf{U}$ , that is, there exists a diagonal matrix $\mathbf{D} such that

$$\begin{equation*} \mathbf{S} = \mathbf{U} \mathbf{D} \mathbf{U}^{-1} . \end{equation*}$$

(iv) Conclude that this leads to a contradiction, by taking into account the Jordan form of $\mathbf{S}_{..2}$.

</p> 
</div>

---

3.1) Para provar a que o tensor $\mathcal{X}$ tem posto 3 pode-se observar os *unfoldings* de modo 1, 2 e 3 do tensor *core* $\mathcal{S}$

$$\begin{equation*} \mathbf{[S]}_{(1)} = 
\begin{bmatrix}
1 & 0 & 0 & 1 \\
0 & 1 & 0 & 0\\
\end{bmatrix}
\end{equation*}$$

$$\begin{equation*} \mathbf{[S]}_{(2)} = 
\begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 1 & 0\\
\end{bmatrix}
\end{equation*}$$

$$\begin{equation*} \mathbf{[S]}_{(3)} = 
\begin{bmatrix}
1 & 0 & 0 & 1 \\
0 & 0 & 1 & 0\\
\end{bmatrix}
\end{equation*}$$

3.2) A partir da equação de representaçã o $\mathcal{X}$ em função da notação de matrizes *slices* (slide 162/244), obtém-se:

$$\begin{equation*} \mathbf{[X]}_{(1)} =  \mathbf{A} \mathbf{[S]}_{(1)} ( \mathbf{C} \otimes \mathbf{B})^{\top} \end{equation*}$$

3.3) Ao substituir os respectivos valores nas matrizes na expressão:

$$\begin{equation*} \mathbf{[X]}_{(1)} = [\mathbf{a}_{1} \quad \mathbf{a}_{2}] 
\begin{bmatrix}
1 & 0 & 0 & 1 \\
0 & 1 & 0 & 0\\
\end{bmatrix} 
\left[(\mathbf{c}_{1} \otimes \mathbf{b}_{1})^{\top} \quad (\mathbf{c}_{1} \otimes \mathbf{b}_{2})^{\top} \quad (\mathbf{c}_{2} \otimes \mathbf{b}_{1})^{\top} \quad (\mathbf{c}_{2} \otimes \mathbf{b}_{2})^{\top} \right] \end{equation*}$$

$$\begin{equation*} \mathbf{[X]}_{(1)} = [\mathbf{a}_{1} \quad \mathbf{a}_{2}] 
\begin{bmatrix}
1 (\mathbf{c}_{1} \otimes \mathbf{b}_{1})^{\top} + 0(\mathbf{c}_{1} \otimes \mathbf{b}_{2})^{\top} + 0(\mathbf{c}_{2} \otimes \mathbf{b}_{1})^{\top} + 1(\mathbf{c}_{2} \otimes \mathbf{b}_{2})^{\top} \\
0 (\mathbf{c}_{1} \otimes \mathbf{b}_{1})^{\top} + 1(\mathbf{c}_{1} \otimes \mathbf{b}_{2})^{\top} + 0(\mathbf{c}_{2} \otimes \mathbf{b}_{1})^{\top} + 0(\mathbf{c}_{2} \otimes \mathbf{b}_{2})^{\top} \\
\end{bmatrix} 
\end{equation*}$$

$$\begin{equation*} \mathbf{[X]}_{(1)} = [\mathbf{a}_{1} \quad \mathbf{a}_{2}] 
\begin{bmatrix}
(\mathbf{c}_{1} \otimes \mathbf{b}_{1})^{\top} + (\mathbf{c}_{2} \otimes \mathbf{b}_{2})^{\top} \\
(\mathbf{c}_{1} \otimes \mathbf{b}_{2})^{\top}\\
\end{bmatrix} 
\end{equation*}$$

$$\begin{equation*} \mathbf{[X]}_{(1)} = \mathbf{a}_{1} \left[(\mathbf{c}_{1} \otimes \mathbf{b}_{1})^{\top} + (\mathbf{c}_{2} \otimes \mathbf{b}_{2})^{\top} \right] + \mathbf{a}_{2}\left[(\mathbf{c}_{1} \otimes \mathbf{b}_{2})^{\top} \right]
\end{equation*}$$

$$\begin{equation*} \mathbf{[X]}_{(1)} = \mathbf{a}_{1} (\mathbf{c}_{1} \otimes \mathbf{b}_{1})^{\top} + \mathbf{a}_{1} (\mathbf{c}_{2} \otimes \mathbf{b}_{2})^{\top} + \mathbf{a}_{2} (\mathbf{c}_{1} \otimes \mathbf{b}_{2})^{\top}
\end{equation*}$$

3.4) A substituição avança, convenientemente para uma forma que pode ser justamente rescrita em função do produto externo dos vetores:

$$\begin{equation*}\mathcal{X} = \mathbf{a}_{1} \circ \mathbf{b}_{1} \circ \mathbf{c}_{1} + \mathbf{a}_{1} \circ \mathbf{b}_{2} \circ \mathbf{c}_{2}  + \mathbf{a}_{2} \circ \mathbf{b}_{2} \circ \mathbf{c}_{1}
\end{equation*}$$

3.5) A partir da demonstração e dos resultados do problema 2 (escrita do tensor *core* em função da decomposição CP), pode-se provar que $\mathcal{X}$ e $\mathcal{S}$ tem posto 3, como sugerido em (ii).

$$\begin{equation*} \mathcal{S} = \sum_{r=1}^{2} \mathbf{s}_{r}^{(1)} \circ \mathbf{s}_{r}^{(2)} \circ \mathbf{s}_{r}^{(3)} \end{equation*}$$

3.6) Com as ferramentas fornecidas, pode-se utilizar a equação (28) das notas de aula para aplicar a decomposição CP com a notação de *slices* frontais.

$$\begin{equation*} \mathcal{S}_{..1} = \sum_{r=1}^{2} \mathbf{s}_{1,r}^{(3)} \circ \mathbf{s}_{r}^{(1)} \circ {\mathbf{s}_{r}^{(2)}}^{\top} = \mathbf{S}^{(1)} \mathbf{D}_{1} \left( \mathbf{S}^{(3)} \right) {\mathbf{S}^{(2)}}^{\top} \end{equation*}$$

$$\begin{equation*} \mathcal{S}_{..2} = \sum_{r=1}^{2} \mathbf{s}_{2,r}^{(3)} \circ \mathbf{s}_{r}^{(1)} \circ {\mathbf{s}_{r}^{(2)}}^{\top} = \mathbf{S}^{(1)} \mathbf{D}_{2} \left( \mathbf{S}^{(3)} \right) {\mathbf{S}^{(2)}}^{\top} \end{equation*}$$

3.7) Levando em consideração o que é sugerido em (ii) e (iii), relacionando a expressão com a decomposição em valores singulares (SVD) 

$$\begin{equation*} \mathcal{S}_{..1} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{\top} = \mathbf{S}^{(1)} \mathbf{D}_{1} \left( \mathbf{S}^{(3)} \right) {\mathbf{S}^{(2)}}^{\top} = \mathbf{I} \end{equation*}$$

3.8) E como premissa: $\mathcal{S}_{..2} = \mathbf{I}$, isto implica que cada um dos termos da própria decomposição é uma matriz $2 \times 2$, tal que:

$$\begin{equation*} \mathbf{S}^{(1)} = \mathbf{D}_{1} \left( \mathbf{S}^{(3)} \right) = \mathbf{S}^{(2)} = \mathbf{U}  = \mathbf{\Sigma} = \mathbf{V}^{\top} = \mathbf{I} \end{equation*}$$

3.9) Como demonstrado acima, $\mathcal{S}_{..2}$ pode ser obtido também a partir dos resultados em função de $\mathcal{S}_{..1}$, com a multiplicação a esquerda por $\mathbf{S}^{(1)}$ e a direita ${\mathbf{S}^{(1)}}^{-1}$:

$$\begin{equation*} \mathcal{S}_{..2} =  \left[  \mathbf{S}^{(1)} \mathbf{D}_{2} \left( \mathbf{S}^{(3)} \right) {\mathbf{S}^{(2)}}^{\top} \right]  \end{equation*}$$

$$\begin{equation*} \mathcal{S}_{..2} = \mathbf{S}^{(1)} \left[  \mathbf{S}^{(1)} \mathbf{D}_{2} \left( \mathbf{S}^{(3)} \right) {\mathbf{S}^{(2)}}^{\top} \right] {\mathbf{S}^{(1)}}^{-1}  \end{equation*}$$

$$\begin{equation*} \mathcal{S}_{..2} =
\begin{bmatrix}
1 & 0\\
0 & 1\\ \end{bmatrix} 
\left\{ \begin{bmatrix}
1 & 0\\
0 & 1\\ \end{bmatrix} \mathbf{D}_{2} \left( \mathbf{S}^{(3)} \right) 
\begin{bmatrix}
1 & 0\\
0 & 1\\ \end{bmatrix} \right\}  \begin{bmatrix}
1 & 0\\
0 & 1\\ \end{bmatrix} \end{equation*}$$

3.10) Dado que a identidade é o elemento neutro na multiplicação de matrizes, podemos reescrever a equação omitindo-as:

$$\begin{equation*} \mathcal{S}_{..2} = \mathbf{D}_{2} \left( \mathbf{S}^{(3)} \right) \end{equation*}$$

3.11) Isso permite observar $\mathcal{S}_{..2}$ como diagonal, mas diferente da forma de Jordan, i.e, o posto do tensor é diferente de 2. 

$$\begin{equation*} \mathcal{S}_{..2} = \begin{bmatrix}
\mathbf{s}_{21} & 0\\
0 & \mathbf{s}_{22}\\ \end{bmatrix} \end{equation*}$$

3.12) Finalmente, o posto do tensor 
$$\begin{equation*}\text{posto}(\mathcal{X}) = \text{posto}(\mathcal{S}) = 3 \end{equation*}.$$

<!---------------------------------------------------------------------------->

# Problem 4

<div style="background-color:rgb(100, 100, 100, 0.15); text-align:left; padding:20px">
<p>

In this last exercise, we will show that, although the tensors of the form considered in the last exercise have rank $3$, they are limits of sequences of rank-$2$ tensors. Thus, unlike happens for matrices, a sequence of rank-$R$ tensors can converge to a rank-$S$ tensor with $S > R$.

(i) First, show that the rank-$1$ tensor

$$\begin{equation*}\mathcal{Y}_{m} = m(\mathbf{a}_{1} + m^{-1} \mathbf{b}_{2}) \circ (\mathbf{b}_{2} + m^{-1} \mathbf{b}_{1}) \circ (\mathbf{c}_{1} \circ m^{-1}\mathbf{c}_{2}) \end{equation*}$$

is equal to $\mathcal{X}$ (as given by (1)) plus an $O(m)$ term $\mathcal{Z}_m$ and an $O(1/m)$ term.

(ii) Subtract the $O(m)$ term to get:

$$\begin{equation*}\mathcal{X}_{m} = \mathcal{Y}_m - \mathcal{Z}_m . \end{equation*}$$

What is the rank of $\mathcal{X}_m$?

(iii) Use the expression obtained for $\mathcal{X}_m$ to conclude that
$\underset{m \rightarrow \infty}{\lim} \mathcal{X}_m = \mathcal{X}.$ 

</p> 
</div>

---

$$\begin{align*} \vdots
\end{align*}$$