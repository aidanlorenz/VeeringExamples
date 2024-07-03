# Small dilatation pseudo-Anosovs via veering triangulations

The code in this repository was written as part of the author's PhD studies.

Dependencies include [Sage 10.1](https://www.sagemath.org/index.html), [Mathematica](https://www.wolfram.com/mathematica/?src=google&416&gad_source=1&gclid=CjwKCAiAopuvBhBCEiwAm8jaMVFCN6LFjhjLvtFWyESeESXPOsPQ2CPjl-eNtUbDTQMpBW2gK9Nl5RoCPJ4QAvD_BwE) (installed with a command line interface that runs when given the `math` command), [Regina](https://regina-normal.github.io/), and [Parlak, Segerman, and Schleimer's Veering repository](https://github.com/henryseg/Veering).

The most important folder is `General` which contains folders with 2-dimensional manifold examples, 3-dimensional manifold examples, and the most important file `VeeringGeneralFunctions.ipynb` where all the methods are written.

`DifferentFace.ipynb` has some examples of things that might be malfunctioning.

## Example usage
In a Jupyer notebook running a Sage kernel is a user-friendly environment to experiment with the code here.
Consider the manifold with the following signature from the veering census:

    sig = 'ivLLQQccdgffhghhvvaaaaavv_01111220'

Call the manifold $M$. The veering structure specifies a fibered face $F$ of $M$ for which we can obtain a Hilbert basis:

```
get_fibered_cone(sig).Hilbert_basis()
```
```
N(0, 0,  1),
N(1, 0, -1),
N(1, 1,  2),
N(2, 1,  0)
in 3-d lattice N
```

To get an expression for the Thurston norm in this fibered cone:
```
get_Thurston_norm(sig)
```
```
2*a0 - 2*a1 + a2
```
where `(a0, a1, a2)` is an arbitrary point in the cone over $F$.

We can easily find an approximation of the minimal normalized dilatation in $F$:
```
get_min_norm_dila_log_mathematica_approx(sig)
```
```
3.6538111277006253
```
and we can also find the vector pointing in the direction of the minimal normalized dilatation:
```
get_minimal_direction(sig)
```
```
[(2.2319678240391077, 0.81012999779366, 0.81012999779366)]
```
The manifold $M$ has $3$ boundary components $b_0, b_1, b_2$ and we can find the maps $\partial_i \colon H_2(M, \partial M) \to H_1(b_i)$ for each $i$ coming from the long exact sequence for the pair $(M, b_i)$:
```
get_LES_boundary_matrices(sig)
```
```
[
[-3  6 -2]    [-1  3  0]    [ 0  0 -1]
[ 0  1  0],   [ 1  0  0],   [ 0 -2  0]
]
```
indicating that $\partial_0$ is represented by the matrix $\begin{pmatrix} -3 & 6 & -2 \\ 0 & 1 & 0 \end{pmatrix}$ and so on.

We can also figure out a lot about any given class.
Pick, for example, the class $(2,1,3)$.
We know this class is in the interior of the fibered cone over $F$ because it can be written as a positive linear combination of three of the hilbert basis vectors:
$$(2,1,3) = 1*(1,1,2) + 1*(1,0,-1) + 2*(0,0,1).$$
Thus we can use the formula obtained above for the Thurston norm:
$$ \|(2,1,3)\| = 2*2 - 2*1 + 3 = 5$$
and we can determine the dilatation of this class:
```
get_dila_mathematica(sig, (2,1,3))
```
```
2.296630262886538
```
Therefore we can compute the normalized dilatation of the class:
$$\overline{\lambda}(2,1,3) = \|(2,1,3)\|\log \lambda(2,1,3) = 5\log(2.296630262886538) \approx 4.157214727646552$$
which, unsurprisingly, is larger than the minimum we found previously.
We can even determine the number of singularities on the minimal surface representative for the class in the orbits corresponding to each boundary component of $M$.
This is just the $\gcd$ of the image under the boundary maps:
$$
\gcd\left(\partial_0\left(2,1,3\right)\right) = \gcd\left(\begin{pmatrix} -3 & 6 & -2 \\ 0 & 1 & 0 \end{pmatrix} \cdot \begin{pmatrix} 2 \\ 1 \\ 3 \end{pmatrix}\right) = \gcd\left(-6, 1\right) = 1
\\
\gcd\left(\partial_1\left(2,1,3\right)\right) = \gcd\left(\begin{pmatrix} -1 & 3 & 0 \\ 1 & 0 & 0 \end{pmatrix} \cdot \begin{pmatrix} 2 \\ 1 \\ 3 \end{pmatrix}\right) = \gcd\left(1, 2\right) = 1
\\
\gcd\left(\partial_2\left(2,1,3\right)\right) = \gcd\left(\begin{pmatrix} 0 & 0 & -1 \\ 0 & -2 & 0 \end{pmatrix} \cdot \begin{pmatrix} 2 \\ 1 \\ 3 \end{pmatrix}\right) = \gcd\left(-3, -2\right) = 1
$$
and hence the minimal surface representative for $(2,1,3)$ has $3$ singularities; $1$ corresponding to each of the $3$ boundary components of $M$.

One can even determine the number of prongs of the invariant foliations for the monodromy associated to a fibered class at each singularity; one might investigate the example files in this repository to see such examples.

It is also worth noting that most functions are designed to work with variables so that one can obtain *formulas* for the quantities discussed above (and others) and hence study whole classes of examples at once instead of looking at individual classes and obtaining individual examples.

## Associated works
This codebase was developed as part of the author's PhD thesis from Vanderbilt University titled "Small dilatation pseudo-Anosovs coming from Dehn fillings of hyperbolic fibered $3$-manifolds" which should be available to the public by August of 2025. Some of the algorithms and examples in this codebase are discussed in detail in the thesis.

A paper which distills the results in said thesis is in preparation and will be linked to here when complete.