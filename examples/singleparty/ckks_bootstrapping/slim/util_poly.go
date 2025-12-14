package main

import (
	"math"
	"math/cmplx"
)

// Input : Ciphertext after EvalMod
// Output : Ciphertext exp(x) to x with cleaning
// Note : deg = 2*n -1
func HermiteInterpolation_exp2id(n int, deg int) []complex128 {
	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(n), 0))

	A := make([][]complex128, 2*n)
	b := make([]complex128, 2*n)

	for i := 0; i < n; i++ {
		zi := cmplx.Pow(z, complex(float64(i), 0))

		A[i] = make([]complex128, 2*n)
		for j := 0; j < 2*n; j++ {
			A[i][j] = cmplx.Pow(zi, complex(float64(j), 0))
		}
		b[i] = complex(float64(i), 0)

		A[n+i] = make([]complex128, 2*n)
		for j := 1; j < 2*n; j++ {
			A[n+i][j] = complex(float64(j), 0) * cmplx.Pow(zi, complex(float64(j-1), 0))
		}
		b[n+i] = 0
	}

	list := solveLinearSystem(A, b)
	coeffs := make([]complex128, deg+1)
	for i := 0; i < 2*n; i++ {
		coeffs[i] = list[i]
	}
	for i := 2 * n; i <= deg; i++ {
		coeffs[i] = complex(0.0, 0.0)
	}

	return coeffs
}

func Coefficient_Decompose_to_four_Polynomial(coeff []complex128) (complex128, [][]complex128) {
	rtn := make([][]complex128, 4)
	t := len(coeff)
	t_h := int(math.Floor(float64(t) / 2.0))
	P_x := make([]complex128, t)
	P_x_bar := make([]complex128, t)
	Q_x := make([]complex128, t_h)
	Q_x_bar := make([]complex128, t_h)

	coeff_rev := make([]complex128, len(coeff))
	for k := 1; k < len(coeff); k++ {
		coeff_rev[k] = coeff[t-k]
	}

	for k := 1; k < t; k++ {
		if k <= t/2 {
			P_x[k] = coeff[k] / complex(float64(t-k)/float64(t)*float64(k+1), 0.0)
		} else {
			P_x[k] = coeff[k] / complex(float64(t-k)/float64(t), 0.0)
		}
	}

	for k := 1; k < t; k++ {
		if k <= t/2 {
			P_x_bar[k] = coeff_rev[k] / complex(float64(t-k)/float64(t)*float64(k+1), 0.0)
		} else {
			P_x_bar[k] = coeff_rev[k] / complex(float64(t-k)/float64(t), 0.0)
		}
	}

	// scaling i
	for k := 1; k < t; k++ {
		if k%4 == 0 {
			P_x_bar[k] *= 1.0 * P_x_bar[k]
		} else if k%4 == 1 {
			P_x_bar[k] *= 1i * P_x_bar[k]
		} else if k%4 == 2 {
			P_x_bar[k] *= -1 * P_x_bar[k]
		} else if k%4 == 3 {
			P_x_bar[k] *= -1i * P_x_bar[k]
		}
	}

	for k := 1; k < t_h; k++ {
		Q_x[k] = coeff[k] / complex(float64(t-k)/float64(-t)*float64(k), 0.0)
	}

	for k := 1; k < t_h; k++ {
		Q_x_bar[k] = coeff_rev[k] / complex(float64(t-k)/float64(-t)*float64(k), 0.0)
	}

	// scaling i
	for k := 1; k < t_h; k++ {
		if k%4 == 0 {
			Q_x_bar[k] *= 1.0 * Q_x_bar[k]
		} else if k%4 == 1 {
			Q_x_bar[k] *= 1i * Q_x_bar[k]
		} else if k%4 == 2 {
			Q_x_bar[k] *= -1 * Q_x_bar[k]
		} else if k%4 == 3 {
			Q_x_bar[k] *= -1i * Q_x_bar[k]
		}
	}

	rtn[0] = P_x
	rtn[1] = P_x_bar
	rtn[2] = Q_x
	rtn[3] = Q_x_bar

	return coeff[0], rtn
}

func HermiteInterpolation_exp2symbol(n int, deg int, base int) []complex128 {
	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(n), 0))

	A := make([][]complex128, 2*n)
	b := make([]complex128, 2*n)

	for i := 0; i < n; i++ {
		zi := cmplx.Pow(z, complex(float64(i), 0))

		A[i] = make([]complex128, 2*n)
		for j := 0; j < 2*n; j++ {
			A[i][j] = cmplx.Pow(zi, complex(float64(j), 0))
		}

		if i >= base {
			b[i] = complex(0, float64(1))
		} else if i == base-1 {
			b[i] = complex(float64(0.5), 0)
		} else {
			b[i] = complex(float64(0), 0)
		}

		A[n+i] = make([]complex128, 2*n)
		for j := 1; j < 2*n; j++ {
			A[n+i][j] = complex(float64(j), 0) * cmplx.Pow(zi, complex(float64(j-1), 0))
		}
		b[n+i] = 0
	}

	list := solveLinearSystem(A, b)
	coeffs := make([]complex128, deg+1)
	for i := 0; i < 2*n; i++ {
		coeffs[i] = list[i]
	}
	for i := 2 * n; i <= deg; i++ {
		coeffs[i] = complex(0.0, 0.0)
	}
	return coeffs
}

func HermiteInterpolation_exp2negatesymbol(n int, deg int, base int) []complex128 {
	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(n), 0))

	A := make([][]complex128, 2*n)
	b := make([]complex128, 2*n)

	for i := 0; i < n; i++ {
		zi := cmplx.Pow(z, complex(float64(i), 0))

		A[i] = make([]complex128, 2*n)
		for j := 0; j < 2*n; j++ {
			A[i][j] = cmplx.Pow(zi, complex(float64(j), 0))
		}

		if i >= base {
			b[i] = complex(0, float64(1))
		} else if i == 0 {
			b[i] = complex(float64(0.5), 0)
		} else {
			b[i] = complex(float64(0), 0)
		}

		A[n+i] = make([]complex128, 2*n)
		for j := 1; j < 2*n; j++ {
			A[n+i][j] = complex(float64(j), 0) * cmplx.Pow(zi, complex(float64(j-1), 0))
		}
		b[n+i] = 0
	}

	list := solveLinearSystem(A, b)
	coeffs := make([]complex128, deg+1)
	for i := 0; i < 2*n; i++ {
		coeffs[i] = list[i]
	}
	for i := 2 * n; i <= deg; i++ {
		coeffs[i] = complex(0.0, 0.0)
	}
	return coeffs
}

func HermiteInterpolation_exp2bin(n int, deg int, base int) []complex128 {
	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(n), 0))

	A := make([][]complex128, 2*n)
	b := make([]complex128, 2*n)

	for i := 0; i < n; i++ {
		zi := cmplx.Pow(z, complex(float64(i), 0))

		A[i] = make([]complex128, 2*n)
		for j := 0; j < 2*n; j++ {
			A[i][j] = cmplx.Pow(zi, complex(float64(j), 0))
		}

		if i == 1 {
			b[i] = complex(float64(0), 0)
		} else {
			b[i] = complex(float64(1), 0)
		}

		A[n+i] = make([]complex128, 2*n)
		for j := 1; j < 2*n; j++ {
			A[n+i][j] = complex(float64(j), 0) * cmplx.Pow(zi, complex(float64(j-1), 0))
		}
		b[n+i] = 0
	}

	list := solveLinearSystem(A, b)
	coeffs := make([]complex128, deg+1)
	for i := 0; i < 2*n; i++ {
		coeffs[i] = list[i]
	}
	for i := 2 * n; i <= deg; i++ {
		coeffs[i] = complex(0.0, 0.0)
	}
	return coeffs
}

// n = 3 : 0 0.5 i
func HermiteInterpolation_symbol2symbol() []complex128 {

	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(32), 0))

	A := make([][]complex128, 4)
	b := make([]complex128, 4)

	for i := 0; i < 2; i++ {
		zi := cmplx.Pow(z, complex(float64(i), 0))

		A[i] = make([]complex128, 4)
		for j := 0; j < 4; j++ {
			A[i][j] = cmplx.Pow(zi, complex(float64(j), 0))
		}

		if i == 0 {
			b[i] = complex(float64(0), 0)
		} else {
			b[i] = complex(float64(0.5), 0)
		}

		A[2+i] = make([]complex128, 4)
		for j := 1; j < 4; j++ {
			A[2+i][j] = complex(float64(j), 0) * cmplx.Pow(zi, complex(float64(j-1), 0))
		}
		b[2+i] = 0
	}

	list := solveLinearSystem(A, b)
	//fmt.Println(len(list))
	coeffs := make([]complex128, 4)
	for i := 0; i < 4; i++ {
		coeffs[i] = list[i]
	}

	return coeffs
}

func HermiteInterpolation_symbol2symbol_imag() []complex128 {

	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(16), 0))

	A := make([][]complex128, 4)
	b := make([]complex128, 4)

	for i := 0; i < 2; i++ {
		zi := cmplx.Pow(z, complex(float64(i), 0))

		A[i] = make([]complex128, 4)
		for j := 0; j < 4; j++ {
			A[i][j] = cmplx.Pow(zi, complex(float64(j), 0))
		}

		if i == 0 {
			b[i] = complex(float64(0), 0)
		} else {
			b[i] = complex(float64(0), 1.0)
		}

		A[2+i] = make([]complex128, 4)
		for j := 1; j < 4; j++ {
			A[2+i][j] = complex(float64(j), 0) * cmplx.Pow(zi, complex(float64(j-1), 0))
		}
		b[2+i] = 0
	}

	list := solveLinearSystem(A, b)
	//fmt.Println(len(list))
	coeffs := make([]complex128, 4)
	for i := 0; i < 4; i++ {
		coeffs[i] = list[i]
	}

	return coeffs
}

func solveLinearSystem(A [][]complex128, b []complex128) []complex128 {
	n := len(b)
	for i := 0; i < n; i++ {
		maxRow := i
		for k := i + 1; k < n; k++ {
			if cmplx.Abs(A[k][i]) > cmplx.Abs(A[maxRow][i]) {
				maxRow = k
			}
		}
		A[i], A[maxRow] = A[maxRow], A[i]
		b[i], b[maxRow] = b[maxRow], b[i]

		pivot := A[i][i]
		for j := i; j < n; j++ {
			A[i][j] /= pivot
		}
		b[i] /= pivot

		for k := i + 1; k < n; k++ {
			factor := A[k][i]
			for j := i; j < n; j++ {
				A[k][j] -= factor * A[i][j]
			}
			b[k] -= factor * b[i]
		}
	}

	x := make([]complex128, n)
	for i := n - 1; i >= 0; i-- {
		x[i] = b[i]
		for j := i + 1; j < n; j++ {
			x[i] -= A[i][j] * x[j]
		}
	}

	return x
}
