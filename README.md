# recursive-nystrom -- Recursive important sampling for Nyström approximation
MATLAB code implementing the recursive ridge leverage score sampling algorithm developed in: [Recursive Sampling for the Nyström Method](https://arxiv.org/abs/1605.07583)

## Installataion

Download `recursiveNystrom.m`, [add to MATLAB path](https://www.mathworks.com/help/matlab/ref/addpath.html), or include directly in project directory. For an example of usage see `exampleApplication.m`.

### Usage
**Input:**
`recursiveNystrom(A,s,kernelFunc,accelerated_flag)`

- `X` : matrix with n rows (data points) and d columns (features)
- `s` : the number of samples used to construct the Nystrom approximation. default = $\sqrt{n}$.
- `kernelFunc` : A function that can compute arbitrary submatrices of A's kernel matrix for some positive semidefinite kernel. For implementation specifics, see the provided example `gaussianKernel.m`. Default = gaussian kernel with $\gamma = 1$.
- `accelerated_flag`: either 0 or 1, default = 0. If the flag is set to 1, the code uses an accelerated version of the algorithm as described in Section 5.2.1 of [NIPS paper.](https://arxiv.org/abs/1605.07583) This version will output a lower quality Nystrom approximation, but will run more quickly. We recommend setting accelerated_flag = 0 (the default) unless the standard version of the algorithm runs too slowly for your purposes.

**Output:**
Rank s Nyström approximation, in factored form.

- `C` : A subset of s columns from A's n x n kernel matrix
- `W` : An s x s positive semidefinite matrix

`C*W*C` approximates `X`'s full kernel matrix, `K`.


In learning applications, it is natural to compute `F = C*chol(W)'`. `F` has n rows and each row can be supplied as a data point to a linear algorithm (regression, SVM, etc.) to approximate the kernel version of the algorithm. Caveat: the accelerated version of our algorithm runs in O(ns) time. Computing `F = C*chol(W)'` takes O(ns^2) time, so it may be more prudent to access the matrix implicitly.

### Example

**Compute a Nyström approximation for a Gaussian kernel matrix**
-- Compute approximation for Gaussian kernel matrix with variance parameter $gamma = 40$, i.e., kernel function for points $x,y$ equal to $e^-(40*\|x - y\|^2)$.

```
gamma = 40;
kFunc = @(X,rowInd,colInd) gaussianKernel(X,rowInd,colInd,gamma);
[C,W] = recursiveNystrom(X,s,kFunc);
```