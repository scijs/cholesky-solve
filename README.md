# cholesky-solve[WIP]

This module solves sparse symmetric positive definite linear systems,
by finding the Cholesky decomposition, and then doing forward
substitution and backward substitution. It is basically a Javascript
port of the paper "Algorithm 8xx: a concise sparse Cholesky
factorization package". This kind of solver has many applications in
digital geometry processing.

## Install

    npm install cholesky-solve

## Example

```javascript
var choleskySolve = require('cholesky-solve')

// matrix dimension.
const n = 10

// sparse matrix on left-hand side
var M = [
  [2, 2, 1.5],
  [1, 1, 1.0],
  [1, 4, 0.02],
  [5, 5, 1.2],
  [7, 7, 1.6],
  [4, 4, 2.6],
  [3, 3, 1.1],
  [4, 7, 0.09],
  [4, 6, 0.16],
  [0, 0, 1.7],
  [4, 8, 0.52],
  [0, 8, 0.13],
  [6, 6, 1.3],
  [7, 8, 0.11],
  [4, 9, 0.53],
  [8, 8, 1.4],
  [9, 9, 3.1],
  [1, 9, 0.01],
  [6, 9, 0.56]
]

// right-hand side
var b = [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]

var P = require('cuthill-mckee')(M, n)
// solve the equation
// Mx = b
// and print x
console.log(choleskySolve(M, b, n, P))
```

## API

### `require("cholesky-solve")(M, b, n, [P])`
Solves the equation `Mx = b` by using the Cholesky decomposition.

* `M` a list of the matrix coefficients of the sparse matrix `M`.
* `b` the vector on the right-hand side. A regular array of length `n`
* `n` the dimension of the matrix `M`
* `P` encodes a permutation matrix that preconditions `M` before the Cholesky decomposition is solved for. A possible algorithm for finding a good permutation is
[Cuthill–McKee](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm). See
the module [cuthill-mckee](https://github.com/mikolalysenko/cuthill-mckee) for a
Javascript implementation.

**Returns** An array encoding the solution to the equation `Mx = b`.
