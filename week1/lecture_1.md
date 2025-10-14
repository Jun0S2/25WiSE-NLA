# [0] 14.10.2025

소유자: Hyejun Park
생성 일시: 2025년 10월 14일 오전 10:19

### General Information

- Pure Math Lectures
- Exam Types
    - Portfolio Exam : Presentation, Oral
    - Num2Eng : Written
    - NLA : Oral Exam (1월에 알려준데 dates는)
    - Oral Exam은 멋이냐? mock exam이 잇을 것이다 +.+
    - 100 프로 Oral 시험으로만 점수가 정해진다.
- Components
    - programming lab
    - Exercise Session
    - Lecture
    - Ex Pres.
    - Exercise Section, EX sheet bi weekly
- No particular textbook
    
    ---
    

---

# Numerical Linear Algebra

### Contents

<aside>

We Study the three main problems of NLA

1. Solve Linear Systems `Ax = b`
2. Least Squares Problems `$min_{x}||b=Ax||$`
3. Eigenvalue Problems `$Ax = \lambda x$`

from a **computational** and **perturbation theory** perrspective.

- computational: specific algorithms
- perturbation : independent of an algorithms

We will see **direct** and **iterative** methods for **dense** and **sparse** matrices.

</aside>

# 1. A survey of matrix decompositions

### Idea :

Factorization of a matrix in *simpler* factors

### Major advantages…

- can be reused to solve several problems involving the same system matrix A.
- can be often updated with saving in computations
- are useful to solve very different problems (`useful` → highly effective software packages!)

We will discuss several important matrix decompositions regarding **existence** and **uniqueness**

### Later

- derivation of algorithms to compute matrix decompositions
- Analysis of their numerical stability

---

## Theorem 1.1 - LU Decomposition

Let $A \in C^{n*m}$ . Then the following are equivalent.

(1) There exists a unit lower triangular $L \in C^{n*m}$, a non singular diagonal $D \in C^{n*m}$, and a unit upper triangular matrix, $U \in C^{n*m}$

- `unit` : ones on the diagonal
- A = LDU

![image.png](%5B0%5D%2014%2010%202025%2028cf023a6c4580d2a46ddb08cb575b58/image.png)

(2) The matrix A (1:k, 1:k) is nonsingular for each k=1 ,…, n

Moreover (2) implies the uniqueness of the matrices in (1)

### Proof “(1) ⇒ (2)”

Let A = LDU. Let $k \in \set{1,...,n}$ . We consider

A = 

| **A11** | **A12** | k rows |
| --- | --- | --- |
| **A21** | **A22** | n-k rows |
| K columns | n-k Columns |  |

$$
= \begin{bmatrix}L_{11} & 0 \\L_{21} & L_{22}\end{bmatrix}\begin{bmatrix}D{n} & 0 \\0 & D_{22}\end{bmatrix}\begin{bmatrix}U_{11} & U_{12} \\0 & U_{22}\end{bmatrix}
$$

From the block entry (1,1), we get

- A(1:k, 1:k) → $A_{11}=L_{11}D_{11}U_{11}$

Since $L$ and $U$ are unit lower/upper triangular, also $L_{ii}$ and $U_{ii}$ are unit lower/upper triangular

$$
det(A) = det(LDU) = det(L)*det(D)*det(U)
$$

<aside>

- det(L) = 1
- det(u) = 1
- det(D) ≠ 0
</aside>

### “(2) ⇒ (1)”

We prove the assertion by induction over n.

Base case - N = 1. Let, 

<aside>

- L = [1] (lower triangular 1 by 1)
- D = [a11]
- U = [1] (upper tri. 1 by 1)
</aside>

Then A = [a11] = LDU

지금 우린 멀 증명하고 잇는건뎅? ㅇㅁㅇ

We assume that the assertion (1) is valid for (n-1) x (n-1) matrices.

Now let $A \in C^{n*m}.$ We consider the blocking

A =

| A11 | A12 | n-1 |
| --- | --- | --- |
| A21 | A22 | 1 |
| n-1 | 1 |  |

A22 : number

A12 : column vector

A21 : row vector

$$
\begin{bmatrix}I_{n-1} & 0 \\A_{21}*A^{-1}_{11} & 1\end{bmatrix}\begin{bmatrix}A{n} & 0 \\0 & s\end{bmatrix}\begin{bmatrix}I_{n-1} & A^{-1}_{11}*A_{12} \\0 & 1\end{bmatrix}
$$

where   $s = A_{22}-A_{21}A^{-1}_{11}A_{12} \in C$

is the Schur complement of An in A. 이거 백퍼 잘못 읽은거같은데 머라씨브리노

From

$$
0 != det(A)  = det(L_{n})det(D_{n}) det(U_{n}) = det (A_{11})*s
$$

- det(Ln)=1
- det(Un)=1
- det(A11) ≠ 0
- s ≠ 0

we deduce s ≠ 0.

deduce 뜻? 추론하다

By the induction hypothesis, the matrix A11 has a factorization $A_{11} = L_{n-1}D_{n-1}U_{n-1}$

Then, A =  $\begin{bmatrix}L_{n-1} & 0 \\0 & 1\end{bmatrix}\begin{bmatrix}D{n-1} & 0 \\0 & s\end{bmatrix}\begin{bmatrix}U_{n-1} & 0\\0 & 1\end{bmatrix}*U_{n}$

= $\begin{bmatrix}L_{n-1} & 0 \\0 & 1\end{bmatrix}\begin{bmatrix}D{n-1} & 0 \\0 & s\end{bmatrix}\begin{bmatrix}U_{n-1} & 0\\0 & 1\end{bmatrix}$

= $\begin{bmatrix}A_{11} & 0 \\0 & s\end{bmatrix}$

- 첫번째 행렬 : unit lower tirangular
- 세번째 행렬: unit upper triangular

→ Which yields the desired decomposition

$$
\begin{bmatrix}L_{n-1}*D_{n-1} & 0 \\0 & s\end{bmatrix}\begin{bmatrix}U_{n-1} & 0 \\0 & 1\end{bmatrix}\begin{bmatrix}L_{n-1}*D_{11}U_{n-1}& 0\\0 & s\end{bmatrix}
$$

### Uniqueness

Let $L1 = [l^{(1)}_{ij}], L2 = [l^{(2)}_{ij}], U1 = [u^{(1)}_{ij}],  U2 = [u^{(2)}_{ij}],$ 

$A = L_1D_1U_1 = L_{2}D_{2}U_2$

Then,   $L^{-1}_2L_1D_1 = D_{2}U_2u^{-1}_1$

<aside>

### Recall

The nonsingular lower/upper triangular matrices form a group w.r.l matrix multiplication

</aside>

- left : lower triangular matrix
- right upper triangular matrix

is both, lower and upper triangular (at the same time), and hence, diagonal matrix

Moreover, by `$l^{1}_{ii} = 1 = l^2_{ii}$`  and `$u^{1}_{ii} = 1 = u^2_{ii}$`  , we find $L^{-1}_2L_1 = I_n$ and $U^{-1}_2U_1 = I_1$

i.e.,  $L_1 = L_2$ and $U_1 = U_2$, and consequently,  $D_1 = D_2$

### Example

(1) For $A = \begin{bmatrix}0 & 1 \\0 & 0\end{bmatrix}$. `Theorem 1.1` does not hold since A(1:1, 1:1) = [0] is singular

However, there exist several LU- decompositions

$$
\begin{bmatrix}0 & 1 \\0 & 0\end{bmatrix} =\begin{bmatrix}1 & 0 \\ -\alpha\beta & 1\end{bmatrix} \begin{bmatrix}1 & 0 \\ 0 & \alpha \end{bmatrix} \begin{bmatrix}0 & 1 \\ 0 & \beta \end{bmatrix} , \alpha\beta \in C
$$

⚠️ beta : not unit upper triangular !!!

(2) For $A = \begin{bmatrix}0 & 1 \\1 & 1\end{bmatrix}$. `Theorem 1.1`  does not hold as well, but there exist **no** LU decomposition of the A=LDU

However, if we swap the rows (permutation), then $PA = \begin{bmatrix}1 & 1 \\0 & 1\end{bmatrix}= \begin{bmatrix}1 & 0 \\0 & 1\end{bmatrix}\begin{bmatrix}1 & 0 \\0 & 1\end{bmatrix} \begin{bmatrix}1 & 1 \\0 & 1\end{bmatrix}$ 

han an LU decomposition an in Thm 1.1 (? not sure what he wrote)

Allowing row permutations (”pivoting”) is also important for the numerical computation of LU-decomposition.

---

## Theorem 1.2 - Corollary ( $LDL^H$ decomposition)

Let $A = A^H \in C^{n*m}$, such that A(1:k, 1:k) is nonsingular for all k=1, n. 

Then there exist a unique unit lower triangular $L \in C^{n*m}$ and a unique nonsingular diagnoal $D \in C^{n*m}$ with $A = LDU^H$