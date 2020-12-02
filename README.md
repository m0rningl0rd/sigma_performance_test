# Singular value decomposition performance test
A simple FORTRAN test to assess performance of several methods to compute singular values of a 3x3 matrix.

Sample output with ```ifort``` 18.0.5 using ```-O2 -mkl=sequential -fp-model precise -static-intel```
```
 *** Computing singular values of matrix A_ij: ***
 ---------------------------------------------------
| 0.39208680E-06 | 0.66691446E+00 | 0.33535504E+00 |
| 0.25480442E-01 | 0.96305549E+00 | 0.91532719E+00 |
| 0.35251614E+00 | 0.83828819E+00 | 0.79586363E+00 |
 ---------------------------------------------------
     1000000 times execution for each method
 ---------------------------------------------------
*** METHOD 1: Directly with LAPACK sgesvd
[ 0.19166092E+01, 0.28137150E+00, 0.17612383E+00]
Time elapsed:  1102.5 ms.

*** METHOD 2: From sqrt of eigenvalues via LAPACK ssyev
[ 0.19166092E+01, 0.28137144E+00, 0.17612405E+00]
Time elapsed:   560.2 ms.

*** METHOD 3: Directly with LAPACK sgesdd
[ 0.19166092E+01, 0.28137150E+00, 0.17612383E+00]
Time elapsed:   941.5 ms.

*** METHOD 4: Self contained algebraic method
[ 0.19166092E+01, 0.28137100E+00, 0.17612576E+00]
Time elapsed:  63.284 ms.
```

Sample output with ```ifort``` 18.0.5 using ```-O3 -mkl=sequential -fp-model fast=2 -no-prec-div -static-intel```
```
  *** Computing singular values of matrix A_ij: ***
 ---------------------------------------------------
| 0.39208680E-06 | 0.66691446E+00 | 0.33535504E+00 |
| 0.25480442E-01 | 0.96305549E+00 | 0.91532719E+00 |
| 0.35251614E+00 | 0.83828819E+00 | 0.79586363E+00 |
 ---------------------------------------------------
     1000000 times execution for each method
 ---------------------------------------------------
*** METHOD 1: Directly with LAPACK sgesvd
[ 0.19166092E+01, 0.28137150E+00, 0.17612383E+00]
Time elapsed:  1100.3 ms.

*** METHOD 2: From sqrt of eigenvalues via LAPACK ssyev
[ 0.19166092E+01, 0.28137147E+00, 0.17612405E+00]
Time elapsed:   549.4 ms.

*** METHOD 3: Directly with LAPACK sgesdd
[ 0.19166092E+01, 0.28137150E+00, 0.17612383E+00]
Time elapsed:  1026.3 ms.

*** METHOD 4: Self contained algebraic method
[ 0.19166092E+01, 0.28136674E+00, 0.17613320E+00]
Time elapsed:   0.010 ms.
```
