# recursive-nystrom: Recursive Importance Sampling for the Nyström Method
MATLAB code implementing the recursive ridge leverage score sampling algorithm developed in: [Recursive Sampling for the Nyström Method](https://arxiv.org/abs/1605.07583) (NIPS 2017).

## Installation

Download `recursiveNystrom.m`, [add to MATLAB path](https://www.mathworks.com/help/matlab/ref/addpath.html), or include directly in project directory. For an example of usage see `exampleApplication.m`.

## Usage
**Input:**
`recursiveNystrom(X,s,kernelFunc,accelerated_flag)`

- `X` : matrix with n rows (data points) and d columns (features)
- `s` : number of samples used in Nyström approximation. default = sqrt(n). Generally should set s < n.
- `kernelFunc` : A function that can compute arbitrary submatrices of `X`'s kernel matrix for some positive semidefinite kernel. For implementation specifics, see the provided example `gaussianKernel.m`. Default = gaussian kernel, i.e. e<sup>-&gamma;||x - y||<sup>2</sup></sup>, with width parameter &gamma; = 1.
- `accelerated_flag`: 0 or 1, default = 0. If the flag is set to 1, the code uses an accelerated version of the algorithm as described in Section 5.2.1 of the [NIPS paper.](https://arxiv.org/abs/1605.07583) This version will output a lower quality Nyström approximation, but run more quickly. We recommend setting accelerated_flag = 0 (the default) unless the standard version of the algorithm runs too slowly for your purposes.

**Output:**
Rank s Nyström approximation, in factored form.

- `C` : A subset of s columns from `X`'s n x n kernel matrix
- `W` : An s x s positive semidefinite matrix

`C*W*C'` approximates `X`'s full kernel matrix, `K`.

In learning applications, it is natural to compute `F = C*chol(W)'`. `F` has n rows and s columns. Each row can be supplied as a data point to a linear algorithm (regression, SVM, etc.) to approximate the kernel version of the algorithm. Caveat: the accelerated version of our algorithm runs in O(ns) time. Computing `F = C*chol(W)'` takes O(ns<sup>2</sup>) time, so it may be more prudent to access the matrix implicitly.

### Example

**Compute a Nyström approximation for a Gaussian kernel matrix with variance parameter &gamma; = 20**

```
% generate random test matrix
X = randn(2000,100); X = normc(X);

% define function for computing kernel dot product
gamma = 20;
kFunc = @(X,rowInd,colInd) gaussianKernel(X,rowInd,colInd,gamma);

% compute factors of Nyström approximation
[C,W] = recursiveNystrom(X,500,kFunc);
```
