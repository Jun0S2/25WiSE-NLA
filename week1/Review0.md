# Review Notes

소유자: Hyejun Park
생성 일시: 2025년 10월 14일 오후 1:32

## **Numerical Linear Algebra (NLA)**

수치선형대수는 **"컴퓨터로 행렬 계산을 어떻게 잘 할 것인가?"** 를 연구하는 학문

<aside>

1. **Solve Linear Systems**
    
    $$
    Ax = b
    $$
    

 ****

- 고등학교 연립방정식의 업그레이드
- **예시:** `2x + y = 5`, `x - y = 1` 이걸 행렬로 `A = [[2, 1], [1, -1]]`, `x = [x, y]^T`, `b = [5, 1]^T` 라고 쓰고  x 를 찾는 문제.
1. **Least Squares Problems**
    
    $$
    min_{x}||b-Ax||
    $$
    
    - `Ax = b`를 정확히 만족하는 해가 없을 때(예: 실험 데이터처럼 오차가 있는 경우), **"에러를 최소화하는"** 가장 나은 해 를 찾는 문제입니다.
    - **예시:** 흩어진 데이터 점들에 가장 잘 맞는 직선을 찾는 **회귀 분석**이 대표적입니다.
2. **Eigenvalue Problems**
    
    $$
    Ax = \lambda x
    $$
    
- 행렬 A 를 곱했을 때 **방향은 바꾸지 않고 크기만 lamda배**로 변하는 특별한 벡터x와 그 배수 lamda를 찾는 문제입니다.
- **예시:** 공학에서 구조물의 진동 모드, 데이터 사이언스에서 PCA(주성분 분석) 등에 폭넓게 쓰입니다.
</aside>

이 문제들을 바라보는 두 가지 관점:

- **Computational :** "어떤 **알고리즘**으로 문제를 풀 것인가?" (예: 가우스 소거법, 반복법)
- **Perturbation Theory :** "입력값에 작은 오차가 생기면 답은 얼마나 달라지는가?" 이는 **행렬 자체의 성질**을 분석하는 것이며, 사용하는 알고리즘과는 무관합니다.

---

풀어야 하는 행렬의 종류:

- Dense Matrices (밀집 행렬): 대부분의 원소가 0이 아님. (예: 일반적인 행렬)
- Sparse Matrices (희소 행렬): 대부분의 원소가 0임. (예: 대규모 편미분 방정식의 계수 행렬)
    - Direct Methods (직접법): 유한한 단계 내에 정확한 해를 찾는 방법 (예: LU 분해). Dense나 중소규모 Sparse에 좋음.
    - Iterative Methods (반복법): 점점 정답에 가까워지는 근사해를 찾는 방법. 대규모 Sparse에 필수.

---

## **Matrix Decompositions**

**아이디어:** "복잡한 행렬 A를, **다루기 쉬운 행렬들의 곱**으로 쪼개자!"

**왜 배우나요? (Major Advantages)**

1. **재사용성 (Reusability):** 한 번 행렬 A를 분해해두면, 다른  b에 대한 선형시스템을 풀 때 **다시 사용**할 수 있어 계산량이 줄어듭니다.
2. **효율성 (Efficiency):** A가 약간만 바뀌었을 때 분해 결과를 **쉽게 업데이트**할 수 있습니다.
3. **범용성 (Versatility):** 하나의 분해 방법으로 선형시스템, 최소제곱, 고유값 등 **아주 다른 문제들을 해결**할 수 있어서, MATLAB, NumPy 등 수치 계산 라이브러리의 기초가 됩니다.

**앞으로 우리가 할 일:**

- 여러 중요한 행렬 분해법의 **존재성과 유일성**을 논의할 것.
- 분해를 계산하는 **알고리즘**을 유도할 것.
- 알고리즘의 **수치적 안정성**을 분석할 것.

---

# **Theorem 1.1 - LU (LDU) Decomposition**

**목표:** 행렬 A를 A=LDU 형태로 분해하는 완벽한 조건을 찾기.

**용어 정의:**

- **Unit Lower Triangular Matrix (L):** 주대각선 아래쪽만 숫자가 있고, **주대각선은 모두 1**인 행렬.
    
    ```
    L = [ 1,  0,  0 ]
        [ ?,  1,  0 ]
        [ ?,  ?,  1 
    ```
    
- **Unit Upper Triangular Matrix (U):** 주대각선 위쪽만 숫자가 있고, **주대각선은 모두 1**인 행렬.
    
    ```
    U = [ 1,  ?,  ? ]
        [ 0,  1,  ? ]
        [ 0,  0,  1 ]
    ```
    
- **Nonsingular Diagonal Matrix (D):** 주대각선만 숫자가 있고, 그 숫자들이 **하나도 0이 아님**.
    
    ```
    D = [ d1, 0,  0 ]
        [ 0, d2,  0 ]
        [ 0,  0, d3 ] (d1, d2, d3 ≠ 0)
    ```
    

---

**Theorem 1.1:**

<aside>

### Lecture Note

Let $A \in C^{n*m}$ . Then the following are equivalent.

(1) There exists a unit lower triangular $L \in C^{n*m}$, a non singular diagonal $D \in C^{n*m}$, and a unit upper triangular matrix, $U \in C^{n*m}$

(2) The matrix A (1:k, 1:k) is `nonsingular` for each k=1 ,…, n

</aside>

**Theorem 1.1 조건:** $A \in C^{n*m}$ (n x n 정사각행렬)

다음 두 명제는 **Equivalent(동치)** 이다.

**(1)** 다음과 같은 분해가 존재한다: A=L⋅D⋅U

- L: Unit Lower Triangular (주대각선이 모두 1인 하삼각행렬)
- D: Non-singular Diagonal (대각성분이 0이 아닌 대각행렬, 즉 역행렬 존재)
- U: Unit Upper Triangular (주대각선이 모두 1인 상삼각행렬)

**(2)** 행렬 A의 모든 **Leading Principal Submatrix** A(1:k,1:k)*A*(1:*k*,1:*k*)가 
k=1,…,n*k*=1,…,*n*에 대해 **Non-singular**이다.

- Leading Principal Submatrix
    
    <aside>
    
    **Leading Principal Submatrix = "왼쪽 위 구석부터 차례로 잘라본 작은 행렬들"**
    
    마치 행렬의 '성장 앨범'이나 '단계별 스냅샷'
    
    ---
    
    ### **구체적인 예시 (Concrete Example)**
    
    다음 3x3 행렬 A가 있다고 가정해봅시다:
    
    ```
    A = [ 1,  2,  3 ]
        [ 4,  5,  6 ]
        [ 7,  8,  9 ]
    ```
    
    **이 행렬의 Leading Principal Submatrix는 3개입니다:**
    
    1. **k=1 (1x1 크기):**
        - **어떻게 자르나요?** 첫 1개 행 & 첫 1개 열을 취합니다.
        - `A(1:1, 1:1) = [ 1 ]`
        - **의미:** 그냥 맨 왼쪽 위 구석의 숫자 하나에요.
    2. **k=2 (2x2 크기):**
        - **어떻게 자르나요?** 첫 2개 행 & 첫 2개 열을 취합니다.
        - `A(1:2, 1:2) = [ 1, 2 ][ 4, 5 ]`
        - **의미:** "왼쪽 위" 2x2 블록입니다.
    3. **k=3 (3x3 크기):**
        - **어떻게 자르나요?** 첫 3개 행 & 첫 3개 열을 취합니다.
        - `A(1:3, 1:3) = [ 1, 2, 3 ][ 4, 5, 6 ][ 7, 8, 9 ]`
        - **의미:** 바로 원본 행렬 A 자신이에요.
    
    ---
    
    ### **왜 중요한가요? (Why It Matters)**
    
    **Theorem 1.1의 조건을 다시 보면:**
    
    > "행렬 A의 모든 Leading Principal Submatrix (k=1, 2, 3)가 Nonsingular (가역적, det ≠ 0)이어야 한다."
    > 
    
    우리 예시 행렬 A로 다시 돌아가보면:
    
    - `k=1`: `[1]` -> det = 1 (**≠ 0**, Good!)
    - `k=2`: `[1, 2; 4, 5]` -> det = (1*5 - 2*4) = -3 (**≠ 0**, Good!)
    - `k=3`: `[1, 2, 3; 4, 5, 6; 7, 8, 9]` -> det = 0 (**= 0**, Bad!❌)
    
    **결론:** 우리 예시 행렬 A는 `k=3`에서 조건을 만족하지 못하므로, Theorem 1.1이 보장하는 '완벽한' LDU 분해를 가질 수 *없습니다*.
    
    ---
    
    ### **다른 예시 (Another Example)**
    
    ```
    B = [ 2,  1,  0 ]
        [ 1,  2,  1 ]
        [ 0,  1,  2 ]
    
    ```
    
    - `k=1`: `[2]` -> det=2 ≠0
    - `k=2`: `[2, 1; 1, 2]` -> det= (4-1)=3 ≠0
    - `k=3`: 원본 B -> det=4 ≠0
    
    **이 행렬 B는 모든 Leading Principal Submatrix가 nonsingular이므로, LDU 분해가 유일하게 존재합니다.**
    
    ### **핵심 요약 (Key Takeaway)**
    
    - **Leading Principal Submatrix:** 그냥 "왼쪽 위 구석부터 `k x k` 크기로 자른 행렬"이에요.
    - **Theorem 1.1:** "작은 크기(1x1, 2x2, ...)부터 하나씩 다 뒤져봤을 때, 전부 역행렬이 존재해야(de
    </aside>
    
- Non-singular
    
    <aside>
    
    **Singular vs. Nonsingular Matrix:**
    
    - A **Nonsingular** matrix is a square matrix that *has an inverse* (is invertible). Another key point: Its determinant is **not zero** (`det(A) ≠ 0`). Its columns are linearly independent.
    - A **Singular** matrix is a square matrix that does *NOT* have an inverse. Its determinant **is zero** (`det(A) = 0`). Its columns are linearly dependent.
    </aside>
    

**추가 명제 :** 조건 (2)가 성립하면, (1)에서의 행렬 L,D,U는 **유일(Unique)** 하다.

---

### Flow

lecture note에 "Moreover (2) implies the uniqueness"라고 쓰여 있음 → 이 때문에 증명이 **세 부분**으로 나뉘어 설명해 준 것임.

---

# Before we start..

### **(1) -> (2)**

<aside>

**"(1)이 성립하면, (2)도 성립한다"** 는 뜻

수학/논리학에서 `->` 또는 `⇒` 기호는 **"If ... then ..."** (만약 ...라면, ...이다)를 의미

</aside>

### Block Partitioning (블록 분할)

<aside>

**큰 행렬을 작은 block으로 잘라서 생각하는 방법**

$$
A = \begin{bmatrix}A_{11} & A_{12} \\A_{21} & A_{22}\end{bmatrix} = \begin{bmatrix}L_{11} & 0 \\L_{21} & L_{22}\end{bmatrix}\begin{bmatrix}D_{11} & 0 \\0 & D_{22}\end{bmatrix}\begin{bmatrix}U_{11} & U_{12} \\0 & U_{22}\end{bmatrix}
$$

- A11 : A의 맨 위 k행 & 맨 왼쪽 k 열로 이루어진 부분 (Leading principal submatrix)
- A12 : 오른쪽 나머지 …
- L,D,U 도 마찬가지

**왜 분할을 하는지 ?** Proof 1에서는 A(1:k, 1:k) 가 non singular인지 확인 해야 함.

→ Simplify with block partitioning, you get $A_{11} = L_{11}D_{11}U_{11}$

→ $det(A_{11}) = det(L_{11})*det(D_{11})*det(U_{11}) = 1 * constant *1$

으로 determinant of A가 0 이 아님을 쉽게 증명할 수 있어서 사용

결론 : 행렬의 일부분만 따로 분석할때 유용한 도구!

**Complex Example** 

$$
A = \begin{bmatrix}1 & 2 & 3 & 4 \\ 5 &6&7&8 \\ 9& 10 &11& 12\\13&14&15&16 \\\end{bmatrix} = \begin{bmatrix} A_{11} & A_{12} \\ A_{21}& A_{22} \end{bmatrix} 
$$

- A11 = 1,2,5,6 (왼쪽 위), A12 = 3,4,7,8 (오른쪽위), A21 = 9,10,13,14 (왼쪽 아래), A22 = 11,12,15,16 (오른쪽 아래) 로 쪼개서 분석
</aside>

---

# 📝 Proves

## **(1) ⇒ (2)**

<aside>

- **가정:** A=LDU 분해가 존재한다.
- **결론:** 모든 Leading Principal Submatrix A(1:k,1:k)*A*(1:*k*,1:*k*)가 Non-singular이다.

즉, A=LDU 분해가 존재하면, 모든 Leading Principal Submatrix 가 Non-singular임을 보인다.

</aside>

### Step 0  Block Partitioning A11

A = LDU에서 k x k Leading Principal Submatrix를 블록으로 나눈다.

$$
A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21}& A_{22} \end{bmatrix}
$$

$A_{11}$이 바로 A(1:k, 1:k) 에 해당된다.

### Step 1  Block Partitioning LDUs

$$
L = \begin{bmatrix} L_{11} & 0 \\ L_{21}& L_{22} \end {bmatrix}, D= \begin{bmatrix} D_{11} & 0 \\ 0& D_{22} \end{bmatrix} , U = \begin{bmatrix} U_{11} & U_{12} \\ 0& U_{22} \end{bmatrix}
$$

### Step 2  Multiplication

$$
A_{11} = L_{11}*D_{11}*U_{11}
$$

### Step 3  Calculation

$$
det(A_{11}) = det(L_{11}D_{11}U_{11}) = det(L_{11})*det(D_{11})*det(U_{11})
$$

- L, U 는 unit triangular → 행렬식 = 1
- D11 는 대각 성분이 0 이 아님 → 행렬식 ≠ 0

### Conclusion  Thus, $A_{11}$ ≠ 0 → $A_{11}$ is Nonsingular

---

## **(2) ⇒ (1)**

<aside>

- **가정:** 모든 Leading Principal Submatrix A(1:k,1:k)*A*(1:*k*,1:*k*)가 Non-singular이다.
- **결론** : A = LDU 분해가 존재한다

즉,  모든 Leading Principal Submatrix 가 Non-singular이면, A = LDU 분해가 존재한다 (by Induction)

</aside>

### Step 0  Base Case

- Let $A = [a_{11}]$
- A(1:1, 1:1) = $[a_{11}]$ is Nonsingular → then, $a_{11}$ ≠ 0
- Thus, $A = [1] * [a_{11}]*1$ 로 분해 가능
- → 각각, $L = [1] , D = [a_{11}], U = [1]$ (왜냐면, A= LDU니까)

### Step 1  Inductive Step

n-1 에서 성립 → n 에서 성립

**Block Partition matrix A (n x n)**

$$
A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21}& A_{22} \end{bmatrix}
$$

 → A₁₁ = (n-1)(n-1), A₂₂ = 1×1

<aside>

**A₁₁ = (n-1)(n-1), A₂₂ = 1×1 은 어디서 왔는가..?**

- 대상 : n x n 행렬 A
- 방법 : n x n 행렬에서, **맨 마지막 행과 맨 마지막 열**을 떼어내서 블록으로 만듦

**Example (n=4):**

$$
A = \begin{bmatrix} a_{11} & a_{12} & a_{13} & a_{14} \\ a_{21}&a_{22}&a_{23}&a_{24}\\a_{31}&a_{32}&a_{33}&a_{34} \\ a_{41} & a_{42} & a_{43} & a_{44} \end{bmatrix}
$$

**Block Partition:**

- $A_{11}$ : 맨 아래 행과 오른쪽 열 하나를 제외한 부분

$$
A_{11} = \begin{bmatrix} a_{11}& a_{12}&a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31}&a_{32}&a_{33} \end{bmatrix}
$$

3 x 3 → 즉, (n-1) x (n-1) array!

- $A_{22}$ : 맨 오른쪽 구석

$$
A_{22} = [a_{44}]
$$

1 x 1 → 즉 그냥 상수

- $A_{12}$  : 오른쪽 위 3 x 1 (column vector - 세로)
- $A_{21}$ : 왼쪽 아래 1 x 3 (row vector - 가로)

**[정리] (n-1) x (n-1) array를 가장 단순한 2x2 행렬로 간소화 하려면 (nx1) x 1, 1x (nx1) 행렬과 상수 하나가 남게 됨. 그 녀석이 바로 a22가 된다.**

</aside>

### Above we got,

$$
A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21}& A_{22} \end{bmatrix}
$$

### Now, Let’s convert L,D,U as 2 x 2 Array

**L (Lower Triangular)** 

**위쪽 오른쪽 블록은 0 이어야함**

$$
L = \begin{bmatrix} L_{11} & 0 \\ L_{21}& L_{22} \end{bmatrix}
$$

**D(Diagonal)** 

**대각선 외에는 0 이어야함**

$$
D = \begin{bmatrix} D_{11} &0\\ 0& D_{22} \end{bmatrix}
$$

**U(Upper Triangular)**

**아래 왼쪽 블록은 0 이어야함**

$$
U = \begin{bmatrix} U_{11} & U_{12} \\ 0& U_{22} \end{bmatrix}
$$

### Now, Compute A = LDU

<aside>

$$
A = LDU =\begin{bmatrix} L_{11} & 0 \\ L_{21}& L_{22} \end{bmatrix} \begin{bmatrix} {D_{11}} & 0 \\ 0 & D_{22} \end{bmatrix}\begin{bmatrix} {U_{11}} & U_{12} \\ 0 & U_{22} \end{bmatrix}
$$

### **1.  L x D**

$$
LD = \begin{bmatrix} L_{11}D_{11} & 0 \\ L_{21}D_{21}& L_{22}D_{22} \end{bmatrix}
$$

### **2.  LD x U**

$$
LDU = LD = \begin{bmatrix} L_{11}D_{11}U_{11} & L_{11}D_{11}U_{12} \\ L_{21}D_{21}U_{11}& L_{21}D_{11}U_{12}+L_{22}D_{22}U_{22} \end{bmatrix}
$$

### 3. A = LDU

$$
A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21}& A_{22} \end{bmatrix}= \begin{bmatrix} L_{11}D_{11}U_{11} & L_{11}D_{11}U_{12} \\ L_{21}D_{21}U_{11}& L_{21}D_{11}U_{12}+L_{22}D_{22}U_{22} \end{bmatrix}
$$

</aside>

We now get 3 equations

- $A_{11} = L_{11}D_{11}U_{11}$
- $A_{12} = L_{11}D_{11}U_{12}$
- $A_{21} = L_{21}D_{21}U_{11}$
- $A_{22} = L_{21}D_{11}U_{12}+L_{22}D_{22}U_{22}$

### Step 2  Schur Complement

<aside>

**Schulr Complement 란 ?**

$$
s = A_{22} - A_{21}*A^{-1}_{11}*A_{12}
$$

- $A_{22}$:  그냥 오른쪽 아래 구석값
- $A_{21}A^{-1}_{11}A_{12} : A_{11}$을 통해 $A_{22}$에 간접적으로 영향을 미치는 값
- s: $A_{22}$ 값에서 $A_{11}$의 간접 영향을 빼면 **진짜 고유한 값**이 뭔가?

**중요한 점:** s는 **그냥 숫자 하나**(1×1 행렬 = 스칼라)

- If $A_{11}$ is Nonsingular, there exist $A^{-1}_{11}$
- $det(A) = det(A_{11}) * s \not = 0$  이고, $det(A_{11}) \not=0$  → $s\not=0$

걍 D_22 를 이 값으로 설정하면 모든 방정식이 깔끔하게 풀린다.

큰 행렬 A 에서 A11의 영향을 제외한 후 A22에 남는 고유한 값을 obtain 하는 공식

$$
A_{11} = L_{n-1}D_{n-1}U_{n-1}
$$

$$
L_{11} = L_{n-1}, D_{11} = D_{n-1}, U_{11} = U_{n-1}
$$

</aside>

### Step 3  Decomposition

By Induction, one can decompose $A_{11} = L_{n-1}*D_{n-1}*U_{n-1}$

[참고] 여기서 n-1 하는 이유는,,, A11로 쪼개면서 맨 마지막 행 열 요소가 날라감. 즉 다른 애들도 n-1 행렬에 맞게 조각내야햄.

$$
L = \begin{bmatrix}L_{n-1} & 0 \\A_{21}*U^{-1}_{n-1}*D^{-1}_{n-1} & 1\end{bmatrix}, D = \begin{bmatrix}D{n-1} & 0 \\0 & s\end{bmatrix}, U = \begin{bmatrix} U_{n-1} & L^{-1}_{n-1}*A_{12} \\ 0 & 1 \end{bmatrix}
$$

- 이 L,D,U를 곱하면 원래 A 가 나옴을 확인할 수 있음

### 방정식 풀기

$A_{12} = L_{11}D_{11}U_{12}$ 풀이

<aside>

$$
A_{12} = L_{n-1}D_{n-1}U_{12}
$$

$$
(L_{n-1}D_{n-1})^{-1}*A_{12} = U_{12}
$$

$$
U_{12}=D^{-1}_{n-1}L^{-1}_{n-1}A_{12}
$$

</aside>

$A_{21} = L_{21}D_{21}U_{11}$  풀이

<aside>

$$
L_{21} = A_{21}U^{-1}_{n-1}D^{-1}_{n-1}
$$

</aside>

$A_{22} = L_{21}D_{11}U_{12}+L_{22}D_{22}U_{22}$ 풀이

<aside>

$$
L = \begin{bmatrix} L_{11}&0\\L_{21}&L_{22} \end{bmatrix}
$$

</aside>

- $L_{22}$ 는 맨 오른쪽 아래에 있는 block, L 은 2x2 행렬이고 L11 은 (n-1) x (n-1) 이므로 $L_{22}$는 1x1 행렬
- U 도 마찬가지
- Unit Triangular 가 되고 싶어서 의도적으로 1을 L과, U 에 집어 넣는다.

<aside>

$$
A_{22} = L_{21}D_{n-1}U_{12} + 1 * D_{22}*1
$$

$$
D_{22} = A_{22} -L_{21}D_{n-1}U_{12}   ... (1)
$$

</aside>

**Left Side** $L_{21}D_{n-1}U_{12}$를 계산해보면,

<aside>

$$
L_{21}D_{n-1}U_{12} = (A_{21}U^{-1}_{n-1}D^{-1}_{n-1})*D_{n-1}*(D^{-1}_{n-1}L^{-1}_{n-1}A_{12})
$$

$$
L_{21}D_{n-1}U_{12} = (A_{21}U^{-1}_{n-1}*I)*(D^{-1}_{n-1}L^{-1}_{n-1}A_{12})
$$

$$
L_{21}D_{n-1}U_{12} = A_{21}U^{-1}_{n-1}*D^{-1}_{n-1}L^{-1}_{n-1}A_{12}
$$

$$
L_{21}D_{n-1}U_{12} = A_{21}(U^{-1}_{n-1}*D^{-1}_{n-1}L^{-1}_{n-1})A_{12}
$$

</aside>

역행렬에 의해,

Tipp ) $A_{11} = L_{n-1}D_{n-1}U_{n-1}$ . 자꾸 11 이 랜덤 숫자로 착각하게 되네.. 기억하자..

$$
L_{21}D_{n-1}U_{12} = A_{21}(A^{-1}_{11})A_{12}
$$

다시 (1) 식으로 되돌아 가서 위 식을 대입하면,

$$
D_{22} = A_{22} - A_{21}A^{-1}_{11}A_{12}
$$

근데? 이게 바로 Schur Complement 임 (D22가 s) 임을 발견 → 이걸로 이제, 행렬의 마지막 성분을 s 로 설정해도 무관하다!

<aside>

**블록 행렬 공식 (외워야함)**

$det(A) = det(A_{11})*s$

→ det(determinant) 는 역행렬이 있는지 알려주는 숫자

- det(A) ≠ 0 → non singular
- det(A) = 0 → singular

우리가 지금 모든 Leading Principal Submatrix가 nonsingular임을 증명하고 있어서 튀어나옴.

</aside>

$$
det(\begin{bmatrix} A_{11}&A_{12} \\ A_{21}& A_{22} \end{bmatrix}) = det(A_{11})det(A_{22} - A_{21}A^{-1}_{11}A_{12})
$$

- $A_{22}-A_{21}A^{-1}_{11}A_{12} = s$  → Schur complement
- $A_{11}$은 nonsingular(가정) → 역행렬 존재
- $det(s) = s$ (s 는 1 x1 행렬이므로)
- Thus, $det(A) = det(A_{11})*s$ ≠ 0 이고, $det(A_{11})$≠ 0 → s≠ 0

### Conclusion  A = LDU

LDU를 곱하면 A 가 나옴을 확인 가능.

<aside>

$$
A_{11} = L_{n-1}D_{n-1}U_{n-1}
$$

$$
L = \begin{bmatrix}L_{n-1}& 0\\ A_{21}U^{-1}_{n-1}D^{-1}_{n-1} & 1 \end{bmatrix} , D = \begin{bmatrix}D_{n-1}& 0\\ 0 & s\end{bmatrix}, U = \begin{bmatrix}U_{n-1}& U_{12}\\ 0 & 1\end{bmatrix}
$$

$$
LD = \begin{bmatrix} L_{n-1}D_{n-1}  & 0 \\ A_{21}U^{01}_{n-1} & {s} \end{bmatrix}
$$

$$
LDU = \begin{bmatrix} L_{n-1}D_{n-1}U_{n-1} & L_{n-1}D_{n-1}U_{12} \\ A_{21}U^{-1}_{n-1}U_{n-1} & A_{21}U^{-1}_{n-1}U_{12}+s \end{bmatrix}
$$

</aside>

### Equations

- $L_{n-1}D_{n-1}U_{n-1} = A_{11}$  (assumption from the induction)
- $L_{n-1}D_{n-1}U_{12} = A_{12}$ (Definition of equation 2)
- $A_{21}U^{-1}_{n-1}U_{n-1} = A_{21}$
- $A_{21}U^{-1}_{n-1}U_{12}+s$
    - We proved, $L_{21}D_{n-1}U_{12} = A_{21}(A^{-1}_{11})A_{12}$
    - 방정식 풀이 2  $U_{12} = D^{-1}_{n-1}L^{-1}_{n-1}A_{12}$ 대입
    - $A_{21}U^{-1}_{n-1}U_12 = A_{21}U^{-1}_{n-1}(D^{-1}_{n-1}L^{-1}_{n-1}A_{12})$
    - $A_{21}(U^{-1}_{n-1}D^{-1}_{n-1}L^{-1}_{n-1})A_{12}$
    - $A_{21}A^{-1}_{11}A_{12}$
    - 다시 대입 및 Schur Complement 실행하면, 
    $A_{21}A^{-1}_{11}A_{12} + (A_{22} - A_{21}A^{-1}_{11}A_{12}) = A_{22}$

→ 모든 블록이 일치함 ! 따라서, A = LDU

---

## Proof 3 - Uniqueness

<aside>

- 가정 : $A = L_1D_1U_1 = L_2D_2U_2$ (두개의 분해가 있다고 가정)
- **결론:** 그 분해를 구성하는 L,D,U는 **오직 하나뿐**이다.
</aside>

### Step 0 식 정리

$$
L^{-1}_2L_1D_1 = D_2U_2U^{-1}_1
$$

- `left side`: $L^{-1}_2L_1$ (unit lower triangular) x $D_1$ (diagonal) ⇒ Lower Triangular
- `right side` : $D_2$ (diagonal) x $U_2U^{-1}_1$ (unit upper triangle) ⇒ Upper Triangular

### Step 1 | Observation

Lower triangular = Upper triangular → 대각 행렬이 되어야함

### Step 2 | Apply Unit req.

- $L^{-1}_2L_1$ 는 unit lower tri. 면서 대각 행렬 → 단위 행렬
- $U_2U^{-1}_1$ 는 unit upper tri. 면서 대각 행렬 → 단위 행렬

$$
L^{-1}_2L_1 = I \implies L_1 = L_2
$$

$$
U_2U^{-1}_1 = I \implies U_1 = U_2
$$

### Conclusion | Unique

L1 = L2, U1 = U2. Thus, D1 = D2

→ Decomposition is unique.

<aside>

What is “Unique”

**"행렬 A를 LDU로 분해하는 방법은 오직 한 가지뿐이다"** 

L,U를 단위 행렬로 강제로 만듦 (LDU 분해의 핵심)

따라서, D도 고정됨 ⇒ 유니크하다.

</aside>

---

# Pivoting

(1) For $A = \begin{bmatrix}0 & 1 \\0 & 0\end{bmatrix}$. `Theorem 1.1` does not hold since A(1:1, 1:1) = [0] is singular

However, there exist several LU- decompositions

$$
\begin{bmatrix}0 & 1 \\0 & 0\end{bmatrix} =\begin{bmatrix}1 & 0 \\ -\alpha\beta & 1\end{bmatrix} \begin{bmatrix}1 & 0 \\ 0 & \alpha \end{bmatrix} \begin{bmatrix}0 & 1 \\ 0 & \beta \end{bmatrix} , \alpha\beta \in C
$$

⚠️ $\beta$ : not unit upper triangular !!!

(2) For $A = \begin{bmatrix}0 & 1 \\1 & 1\end{bmatrix}$. `Theorem 1.1`  does not hold as well, but there exist **no** LU decomposition of the A=LDU

However, if we swap the rows (permutation), then $PA = \begin{bmatrix}1 & 1 \\0 & 1\end{bmatrix}= \begin{bmatrix}1 & 0 \\0 & 1\end{bmatrix}\begin{bmatrix}1 & 0 \\0 & 1\end{bmatrix} \begin{bmatrix}1 & 1 \\0 & 1\end{bmatrix}$ 

Allowing row permutations (”pivoting”) is also important for the numerical computation of LU-decomposition.

---

# Theorem 1.2 - Corollary ( $LDL^H$ decomposition)

Let $A = A^H \in C^{n*m}$, such that A(1:k, 1:k) is nonsingular for all k=1, n. 

Then there exist a unique unit lower triangular $L \in C^{n*m}$ and a unique nonsingular diagnoal $D \in C^{n*m}$ with $A = LDU^H$

---

<aside>

따라서, D도 자동으로 ㄱ

</aside>