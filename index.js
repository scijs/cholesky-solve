function ldl_symbolic
(
  n, /* A and L are n-by-n, where n >= 0 */
  Ap, /* input of size n + 1, not modified */
  Ai, /* input of size nz=Ap[n], not modified */
  Lp, /* output of size n + 1, not defined on input */
  Parent, /* output of size n, not defined on input */
  Lnz, /* output of size n, not defined on input */
  Flag /* workspace of size n, not defn. on input or output */
) {

  var i, k, p, kk, p2

  for (k = 0; k < n; k++) {
    /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
    Parent[k] = -1 /* parent of k is not yet known */
    Flag[k] = k /* mark node k as visited */
    Lnz[k] = 0 /* count of nonzeros in column k of L */
    kk = (k) /* kth original, or permuted, column */
    p2 = Ap[kk + 1]
    for (p = Ap[kk]; p < p2; p++) {
      /* A (i,k) is nonzero (original or permuted A) */
      i = (Ai[p])

      if (i < k) {
        /* follow path from i to root of etree, stop at flagged node */
        for (; Flag[i] !== k; i = Parent[i]) {
          /* find parent of i if not yet determined */
          if (Parent[i] === -1) Parent[i] = k
          Lnz[i]++ /* L (k,i) is nonzero */
          Flag[i] = k /* mark i as visited */
        }
      }
    }
  }
  /* construct Lp index array from Lnz column counts */
  Lp[0] = 0
  for (k = 0; k < n; k++) {
    Lp[k + 1] = Lp[k] + Lnz[k]
  }
}

function ldl_numeric /* returns n if successful, k if D (k,k) is zero */
(
  n, /* A and L are n-by-n, where n >= 0 */
  Ap, /* input of size n+1, not modified */
  Ai, /* input of size nz=Ap[n], not modified */
  Ax, /* input of size nz=Ap[n], not modified */
  Lp, /* input of size n+1, not modified */
  Parent, /* input of size n, not modified */
  Lnz, /* output of size n, not defn. on input */
  Li, /* output of size lnz=Lp[n], not defined on input */
  Lx, /* output of size lnz=Lp[n], not defined on input */
  D, /* output of size n, not defined on input */
  Y, /* workspace of size n, not defn. on input or output */
  Pattern, /* workspace of size n, not defn. on input or output */
  Flag /* workspace of size n, not defn. on input or output */
) {

  var yi, l_ki
  var i, k, p, kk, p2, len, top
  for (k = 0; k < n; k++) {
    /* compute nonzero Pattern of kth row of L, in topological order */
    Y[k] = 0.0 /* Y(0:k) is now all zero */
    top = n /* stack for pattern is empty */
    Flag[k] = k /* mark node k as visited */
    Lnz[k] = 0 /* count of nonzeros in column k of L */
    kk = (k) /* kth original, or permuted, column */
    p2 = Ap[kk + 1]
    for (p = Ap[kk]; p < p2; p++) {
      i = (Ai[p]) /* get A(i,k) */
      if (i <= k) {
        Y[i] += Ax[p] /* scatter A(i,k) into Y (sum duplicates) */
        for (len = 0; Flag[i] !== k; i = Parent[i]) {
          Pattern[len++] = i /* L(k,i) is nonzero */
          Flag[i] = k /* mark i as visited */
        }
        while (len > 0) Pattern[--top] = Pattern[--len]
      }
    }
    /* compute numerical values kth row of L (a sparse triangular solve) */
    D[k] = Y[k] /* get D(k,k) and clear Y(k) */
    Y[k] = 0.0
    for (; top < n; top++) {
      i = Pattern[top] /* Pattern[top:n-1] is pattern of L(:,k) */
      yi = Y[i] /* get and clear Y(i) */
      Y[i] = 0.0
      p2 = Lp[i] + Lnz[i]
      for (p = Lp[i]; p < p2; p++) {
        Y[Li[p]] -= Lx[p] * yi
      }
      l_ki = yi / D[i] /* the nonzero entry L(k,i) */
      D[k] -= l_ki * yi
      Li[p] = k /* store L(k,i) in column form of L */
      Lx[p] = l_ki
      Lnz[i]++ /* increment count of nonzeros in col i */
    }

    if (D[k] === 0.0) return (k) /* failure, D(k,k) is zero */
  }

  return (n) /* success, diagonal of D is all nonzero */
}

function ldl_lsolve
(
  n, /* L is n-by-n, where n >= 0 */
  X, /* size n. right-hand-side on input, soln. on output */
  Lp, /* input of size n+1, not modified */
  Li, /* input of size lnz=Lp[n], not modified */
  Lx  /* input of size lnz=Lp[n], not modified */
) {
  var j, p, p2
  for (j = 0; j < n; j++) {
    p2 = Lp[j + 1]
    for (p = Lp[j]; p < p2; p++) {
      X[Li[p]] -= Lx[p] * X[j]
    }
  }
}
function ldl_dsolve
(
  n, /* D is n-by-n, where n >= 0 */
  X, /* size n. right-hand-side on input, soln. on output */
  D /* input of size n, not modified */
) {
  var j
  for (j = 0; j < n; j++) {
    X[j] /= D[j]
  }
}
function ldl_ltsolve
(
  n, /* L is n-by-n, where n >= 0 */
  X, /* size n. right-hand-side on input, soln. on output */
  Lp, /* input of size n+1, not modified */
  Li, /* input of size lnz=Lp[n], not modified */
  Lx  /* input of size lnz=Lp[n], not modified */
) {
  var j, p, p2
  for (j = n - 1; j >= 0; j--) {
    p2 = Lp[j + 1]
    for (p = Lp[j]; p < p2; p++) {
      X[j] -= Lx[p] * X[Li[p]]
    }
  }
}

function ldl_perm
(
  n,		/* size of X, B, and P */
  X,	/* output of size n. */
  B,	/* input of size n. */
  P	/* input permutation array of size n. */
) {
  var j ;
  for (j = 0; j < n; j++)
  {
    X[j] = B [P[j]]
  }
}

function ldl_permt(
  n,		/* size of X, B, and P */
  X,	/* output of size n. */
  B,	/* input of size n. */
  P	/* input permutation array of size n. */
) {
  var j
  for (j = 0; j < n; j++)
  {
    X [P[j]] = B[j]
  }
}

function prepare (M, n, P) {
  const ANZ = M.length

  // if a permutation was specified, apply it.
  if (P) {
    var Pinv = new Array(n)

    for (k = 0; k < n; k++) {
      Pinv[P[k]] = k
    }

    var Mt = [] // scratch memory
    // Apply permutation. We make M into P*M*P^T
    for(var a = 0; a < M.length; ++a) {
      var ar = Pinv[M[a][0]]
      var ac = Pinv[M[a][1]]

      // we only store the upper-diagonal elements(since we assume matrix is symmetric, we only need to store these)
      // if permuted element is below diagonal, we simply transpose it.
      if(ac < ar) {
        var t = ac
        ac = ar
        ar = t
      }

      Mt[a] = []
      Mt[a][0] = ar
      Mt[a][1] = ac
      Mt[a][2] = M[a][2]
    }

    M = Mt // copy scratch memory.
  } else {
    // if P argument is null, we just use an identity permutation.
    var P = []
    for(var i = 0; i < n; ++i) {
      P[i] = i
    }
  }

  // The sparse matrix we are decomposing is A.
  // Now we shall create A from M.
  var Ap = new Array(n + 1)
  var Ai = new Array(M.length)
  var Ax = new Array(M.length)

  // count number of non-zero elements in columns.
  var LNZ = []
  for(var i = 0; i < n; ++i) {
    LNZ[i] = 0
  }
  for(var a = 0; a < M.length; ++a) {
    LNZ[M[a][1]]++
  }

  Ap[0] = 0
  for(var i = 0; i < n; ++i) {
    Ap[i+1] = Ap[i] + LNZ[i]
  }

  var coloffset = []
  for(var a = 0; a < n; ++a) {
    coloffset[a] = 0
  }

  // go through all elements in M, and add them to sparse matrix A.
  for(var i = 0; i < M.length; ++i) {
    var e = M[i]
    var col = e[1]

    var adr = Ap[col] + coloffset[col]
    Ai[adr] = e[0]
    Ax[adr] = e[2]

    coloffset[col]++
  }

  var D = new Array(n)
  var Y = new Array(n)
  var Lp = new Array(n + 1)
  var Parent = new Array(n)
  var Lnz = new Array(n)
  var Flag = new Array(n)
  var Pattern = new Array(n)
  var bp1 = new Array(n)
  var x = new Array(n)
  var d

  ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag)

  var Lx = new Array(Lp[n])
  var Li = new Array(Lp[n])

  d = ldl_numeric(n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern, Flag)

  if (d === n) {
    return function(b) {
      ldl_perm(n, bp1, b, P);
      ldl_lsolve(n, bp1, Lp, Li, Lx)
      ldl_dsolve(n, bp1, D)
      ldl_ltsolve(n, bp1, Lp, Li, Lx)
      ldl_permt(n, x, bp1, P);

      return x
    }

  } else {
    return null
  }
}

module.exports.prepare = prepare
