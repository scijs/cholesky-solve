var choleskySolve = require('./')
var test        = require('tape')
var almostEqual = require("almost-equal")

// deterministic RNG for generating test data.
var rng = new require('xorshift').constructor([1, 0, 2, 0]);

var eps = 1e-5

function solveAndAssert(t, n, M, b, P, expectedSolution) {
  var foundSolution = choleskySolve(M, b, n, P)

  for(var i=0; i< n; ++i) {
    t.assert(almostEqual(expectedSolution[i], foundSolution[i], eps, eps), "solution element " + i + ": "+ expectedSolution[i] + " = " + foundSolution[i])
  }
}

test('solve10x10matrix', function(t) {
  const n = 10
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

  var b = [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]

  var expectedSolution =[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

  var P = null // identity permutation

  solveAndAssert(t, n, M, b, P, expectedSolution)

  t.end();
})

test('solve10x10matrix_cuthill_mckee', function(t) {
  const n = 10
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

  var b = [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]

  var expectedSolution =[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

  var P = require('cuthill-mckee')(M, n)

  solveAndAssert(t, n, M, b, P, expectedSolution)

  t.end();
})

test('solve1000x1000matrix', function(t) {
  var L = []
  var n = 1000

  // first generate a lower-triangular matrix L
  for(var i = 0; i < n; ++i) {
    L[i] = []

    for(var j = 0; j <n; ++j) {
      L[i][j] = 0
    }

    for(var j = 0; j <= i; ++j) {
      if(rng.random() > 0.9 || i === j) {
        L[i][j] = Math.floor(rng.random() * 10)+1
      } else {
        L[i][j] = 0
      }
    }
  }

  // next, we simply multiply L with its transpose, and put the result in A.
  // the resulting matrix is symmetric, and positive definite,
  // thus it must have a cholesky decomposition.
  var A = []
  for(var i = 0; i < n; ++i) {
    A[i] = []
  }
  for(var i = 0; i < n; ++i) {
    for(var j = 0; j < n; ++j) {
      var s = 0.0
      for(var k = 0; k < n; ++k) {
        s += L[i][k] * L[j][k]
      }
      A[i][j] = s
    }
  }

  // now store A as a sparse matrix M.
  var M = []
  for(var row = 0; row < n; ++row) {
    for(var col = row; col < n; ++col) {
      //      console.log(row, col)
      if(A[row][col] > 0.0001) {
        M.push([row, col, A[row][col]])
      }
    }
  }

  // In our test, we shall solve the equation
  // Mx = b
  // so randomly generate x.
  var x = []
  var b = []
  for(var i = 0; i < n; ++i) {
    x[i] = Math.floor(rng.random() * 9)
  }

  // Now compute b as b = Mx
  for(var i = 0; i < n; ++i) {
    var s = 0.0
    for(var k = 0; k < n; ++k) {
      s += A[i][k] * x[k]
    }
    b[i] = s
  }


  var P = require('cuthill-mckee')(M, n)
  // solve.
  var foundSolution = choleskySolve(M, b, n, P)

  // check that the residual vector is 0.
  for(var i = 0; i < n; ++i) {
    var s = 0.0
    for(var k = 0; k < n; ++k) {
      s += A[i][k] * foundSolution[k]
    }
    var res = b[i] - s
    t.assert(almostEqual(0.0, res, eps, eps), "residual #" + i + ":" +  "0.0 = " + res)
  }

  t.end();
})
