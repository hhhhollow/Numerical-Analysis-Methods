# Numerical Analysis Lab Experiments

## Overview
This repository contains a collection of Python scripts implementing fundamental numerical analysis algorithms. The focus of these implementations is on **first principles**, computational stability, and rigorous step-by-step verification (验算). 

These algorithms are built from scratch using only the standard Python `math` library, intentionally avoiding high-level computational libraries like `SciPy` or `NumPy` to fully demonstrate the underlying mathematical mechanics and error handling.

## Project Structure

### Problem 1: Computational Stability & Overflow Avoidance
**Objective:** Calculate the function values of $g(x) = \frac{e^{-x}}{1 + e^{-x}}$ at extreme values ($x = \pm 750$) and analyze computational errors.
* **Challenge:** Standard evaluation at $x = -750$ causes an `OverflowError` due to exceeding the double-precision floating-point limit ($~1.8 \times 10^{308}$).
* **Solution:** Implemented a numerically robust piecewise algorithm using algebraic identity transformation ($g(x) = \frac{1}{e^x + 1}$ for $x < 0$) to completely eliminate overflow while handling underflow safely.

### Problem 2: Root Finding via Newton-Raphson & Secant Methods
**Objective:** Find the positive root of the non-linear equation $x e^x - x^2 = 1$.
* **Features:**
  * Implements both **Newton's Method** (using analytical derivatives) and the **Secant Method** (using finite difference approximations).
  * Includes strict mathematical verification at each iteration: checking for zero-derivatives (to prevent division by zero) and tracking the residual error $|f(x_k)|$ to observe convergence rates.

### Problem 3: Strict Root Finding via Bisection Method
**Objective:** Find the root of the polynomial $x^5 - 5x^2 + 1 = 0$ in the interval $[-1, 0]$ with a strict absolute error tolerance of $|x_k - x^*| \le 10^{-5}$.
* **Features:**
  * Implements the textbook Bisection Algorithm (Algorithm 2.2.1).
  * Continuously evaluates the interval length $|b-a|$ and the absolute function value $|f(x)|$ as dual stopping criteria.
  * Outputs detailed step-by-step verification, proving that the theoretical minimum iterations align with the practical computational output.

### Problem 4: Large-Scale Sparse Matrix Computation & Cholesky Decomposition
**Objective:** Solve a symmetric positive definite linear system $Ax = B$ utilizing experimental data (`data2.mat`), transitioning from a $100 \times 100$ principal submatrix to a full $6461 \times 6461$ sparse matrix.
* **Features:**
  * Implements standard Cholesky decomposition ($A = LL^T$) and Modified Cholesky decomposition ($A = LDL^T$) purely from scratch, adhering to first principles without calling pre-built matrix factorization functions.
  * Strictly verifies numerical accuracy by calculating the infinity norm of the residual vector $||Ax - B||_{\infty}$ after deriving the solution.
  * Evaluates algorithm scalability and memory efficiency, contrasting the custom loop-based implementations with optimized sparse direct solvers when handling the full-scale sparse dataset.

### Problem 5: Iterative Methods for Linear Systems
**Objective:** Solve a $4 \times 4$ linear system $Ax = b$ using classical stationary iterative methods and compare their convergence behavior under the relative residual stopping criterion $\frac{||b - Ax_k||}{||b||} < 10^{-6}$.
* **Features:**
  * Implements the **Jacobi Method**, **Gauss-Seidel Method**, and **SOR Method** based on the matrix splitting $A = D - L - U$.
  * Uses the zero vector as the initial approximation and tracks the residual norm after each iteration to determine whether the numerical solution has reached the required accuracy.
  * Includes divergence protection for unstable iterations and selects an under-relaxation factor $\omega = 0.1$ for SOR to improve convergence reliability on the given coefficient matrix.
  
## Prerequisites
* Python 3.6 or higher.
* MATLAB for running `problem4.m` and `problem5.m`.
* No external Python libraries required for the Python-based problems.

## How to Run
Clone the repository and run the scripts directly from your terminal or IDE:
```bash
git clone <your-repository-url>
cd <your-repository-directory>
python problem1_overflow.py
python problem2_newton_secant.py
python problem3_bisection.py
matlab -batch "problem4"
matlab -batch "problem5"
```
