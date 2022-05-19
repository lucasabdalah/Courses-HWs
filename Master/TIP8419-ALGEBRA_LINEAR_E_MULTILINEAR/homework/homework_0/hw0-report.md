# TI8419 - Multilinear Algebra
# Lucas Abdalah
## Homework 0

See the section on [`Problem 1`](#Method1)

## Table of Contents
1. [Problem 1](#problem-1)
2. [Problem 2](#problem-2)


## Problem 1
### Item (a) 

- Method 1
  
$$ (\mathbf{A}_{N \times N} \otimes \mathbf{B}_{N \times N} )^{-1}$$

- Method 2

$$ (\mathbf{A}_{N \times N})^{-1} \otimes (\mathbf{B}_{N \times N})^{-1}$$

for $n \in \{2,4,6,8,16,32,64\}$

<p align="center">
<img src="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/master/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/homework_0/hw0-problem1-a.png" alt="alt text" title="image Title" width="512" />
</p>

## Problem 2


### (a) Method 1
$$ (\mathbf{A}_{N \times N} \otimes \mathbf{B}_{N \times N} )^{-1}$$

### (b) Method 2

The script to peform Problem B with matrices NxN and kron productory until K, with N = 4, and K = 8, results in an dimensional error, since it's required more the RAM memory available.

```
Requested 4x16384x4x16384 (64.0GB) array exceeds maximum array size preference (19.8GB). This might cause MATLAB to become unresponsive.

Error in kron (line 35)
   K = reshape(A.*B, ma*mb, na*nb);

Error in problem1_itemb>b_method1 (line 43)
            C = kron(C, A{ii});

Error in problem1_itemb (line 16)
[time, ~] = b_method1(8)
```
<div style="background-color:rgba(250, 0, 0, 0.25); text-align:left; padding:20px">
<p> 

``` 
Requested 4x16384x4x16384 (64.0GB) array exceeds maximum array size preference (19.8GB). This might cause MATLAB to become unresponsive. 

Error in kron (line 35)
   K = reshape(A.*B, ma*mb, na*nb);

Error in problem1_itemb>b_method1 (line 43)
            C = kron(C, A{ii});

Error in problem1_itemb (line 16)
[time, ~] = b_method1(8)
``` 
</p>
</div>



with N = 4, and K = 8


```
problem1_itemb
       16384       16384

time =

  486.4500
```

- - - 

[Hobbit lifestyles][1]


<!-- References -->

[1]: <https://en.wikipedia.org/wiki/Hobbit#Lifestyle> (Hobbit lifestyles)


<!-- <span style="color:red">some **blue** text</span>. -->

<!-- <div style="background-color:rgba(250, 0, 0, 0.25); text-align:center; padding:20px">
<p> Alo teste </p>
</div>

<div style="background-color:rgba(0, 0, 250, 0.25); text-align:center; padding:20px">
<p> Alo teste </p>
</div>

<div style="background-color:rgba(0, 250, 0, 0.25); text-align:center; padding:20px">
<p> Alo teste </p>
</div> -->

<!-- <div style="background-color:rgba(0, 0, 0, 0.0470588); text-align:center; vertical-align: middle; padding:40px 0; margin-top:30px">
<a href="/blog">VIEW THE BLOG</a>
</div> -->