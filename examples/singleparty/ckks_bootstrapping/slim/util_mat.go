package main

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"math/big"
	"math/cmplx"
	"os"
	"sort"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func Test_Gen_matrix() {

	N := 64
	f_coeffs := []complex128{1, 1, 0, 0, 0, 0, 0, 0}
	g_coeffs := []complex128{1, 1, 0, 0, 0, 0, 0, 0}

	logN := int(math.Log2(float64(N)))
	mid := logN / 2

	P := GenerateBitReverseMatrix(N)

	// Head = B_{mid-1} * ... * B_0
	dftHead := Identity(N)
	for i := 0; i < mid; i++ {
		B := GenerateButterflyStage(N, i)
		dftHead = MatMul(B, dftHead)
	}

	// Tail = B_{logN-1} * ... * B_{mid}
	dftTail := Identity(N)
	for i := mid; i < logN; i++ {
		B := GenerateButterflyStage(N, i)
		dftTail = MatMul(B, dftTail)
	}

	scale := 1.0 / math.Sqrt(float64(N))
	idftTail := ScaleMatrix(ConjugateTranspose(dftTail), complex(scale, 0))
	idftHead := ScaleMatrix(ConjugateTranspose(dftHead), complex(scale, 0))

	Row_To_Diagonal(MatMul(dftHead, P))
	Row_To_Diagonal(dftTail)

	// DFT(f) = dftTail * dftHead * f
	f_coeffs = MatvecMul(P, f_coeffs)
	g_coeffs = MatvecMul(P, g_coeffs)

	f_dft := MatvecMul(dftTail, MatvecMul(dftHead, f_coeffs))
	// DFT(g) = dftTail * dftHead * g
	g_dft := MatvecMul(dftTail, MatvecMul(dftHead, g_coeffs))

	h_dft := make([]complex128, N)
	for i := 0; i < N; i++ {
		h_dft[i] = f_dft[i] * g_dft[i]
	}

	h_coeffs := MatvecMul(idftHead, MatvecMul(idftTail, h_dft))

	Row_To_Diagonal(idftTail)
	Row_To_Diagonal(MatMul(P, idftHead))

	fmt.Println("f(x) = 1 + x")
	fmt.Println("g(x) = 1 + x")
	fmt.Println("---------------------------------")
	fmt.Println("Result h(x) = f(x) * g(x) Coefficients:")
	fmt.Println("Expected: 1.00 + 2.00x + 1.00x^2")
	fmt.Print("Actual:   ")
	for _, v := range h_coeffs {
		fmt.Printf("%.2f ", real(v))
	}
	fmt.Println()

	Row_To_Diagonal(dftHead)
	Row_To_Diagonal(dftTail)
	Row_To_Diagonal(idftHead)
	Row_To_Diagonal(idftTail)

}

func GenerateBitReverseMatrix(N int) [][]complex128 {
	P := Zeros(N, N)
	logN := uint(math.Log2(float64(N)))

	for i := 0; i < N; i++ {
		rev_i := 0
		for j := uint(0); j < logN; j++ {
			if (i>>j)&1 == 1 {
				rev_i |= 1 << (logN - 1 - j)
			}
		}
		P[rev_i][i] = 1
	}
	return P
}

func GenerateButterflyStage(N, stage int) [][]complex128 {
	B := Zeros(N, N)
	groupSize := 1 << (stage + 1)
	numGroups := N / groupSize
	twiddleStep := 1 << (uint(math.Log2(float64(N))) - 1 - uint(stage))

	for g := 0; g < numGroups; g++ {
		for k := 0; k < groupSize/2; k++ {
			base := g*groupSize + k
			pair := base + groupSize/2

			// Twiddle factor: w = e^(-2*pi*i*k*step / N)
			angle := -2.0 * math.Pi * float64(k*twiddleStep) / float64(N)
			twiddle := cmplx.Exp(complex(0, angle))

			B[base][base] = 1
			B[base][pair] = twiddle
			B[pair][base] = 1
			B[pair][pair] = -twiddle
		}
	}
	return B
}

func MatMul(A, B [][]complex128) [][]complex128 {
	n, m, p := len(A), len(B[0]), len(B)
	C := Zeros(n, m)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			for k := 0; k < p; k++ {
				C[i][j] += A[i][k] * B[k][j]
			}
		}
	}
	return C
}

func MatvecMul(A [][]complex128, v []complex128) []complex128 {
	n, m := len(A), len(v)
	res := make([]complex128, n)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			res[i] += A[i][j] * v[j]
		}
	}
	return res
}

func ConjugateTranspose(A [][]complex128) [][]complex128 {
	n, m := len(A), len(A[0])
	B := Zeros(m, n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			B[i][j] = cmplx.Conj(A[j][i])
		}
	}
	return B
}

func ScaleMatrix(A [][]complex128, alpha complex128) [][]complex128 {
	n, m := len(A), len(A[0])
	B := Zeros(n, m)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			B[i][j] = alpha * A[i][j]
		}
	}
	return B
}

func Zeros(r, c int) [][]complex128 {
	m := make([][]complex128, r)
	for i := range m {
		m[i] = make([]complex128, c)
	}
	return m
}

func InverseHeadTail(head, tail [][]complex128) ([][]complex128, [][]complex128) {
	n := len(head)
	scale := 1.0 / math.Sqrt(float64(n))

	headInv := ScaleMatrix(ConjugateTranspose(head), complex(scale, 0))
	tailInv := ScaleMatrix(ConjugateTranspose(tail), complex(scale, 0))

	return headInv, tailInv
}

func MatVecScale(x []complex128, alpha complex128) []complex128 {
	n := len(x)
	y := make([]complex128, n)
	for i := 0; i < n; i++ {
		y[i] = x[i] * alpha
	}
	return y
}

func Generate_baseREDCmat(base int, slice_length int, mod *big.Int) [][]complex128 {

	rtn := make([][]complex128, slice_length)
	for i := 0; i < slice_length; i++ {
		rtn[i] = make([]complex128, slice_length)
	}

	radix := make([][]int64, slice_length/2)
	for i := 0; i < slice_length/2; i++ {
		radix[i] = make([]int64, slice_length/2)
	}

	bigbase := big.NewInt(int64(base))

	for i := 0; i < slice_length/2; i++ {
		powerbase := big.NewInt(1.0)
		powerbase.Exp(bigbase, big.NewInt(int64(slice_length/2+i)), mod)
		radix[i] = DecomposeInt64Base(powerbase, int64(base), slice_length/2)
	}

	for i := 0; i < slice_length/2; i++ {
		for j := 0; j < slice_length/2; j++ {
			rtn[i+slice_length/2][j+slice_length/2] = complex(float64(radix[i][j]), 0)
		}
	}

	rtn = Transpose(rtn)

	return rtn
}

func checkRectComplex(A [][]complex128) (m, n int, err error) {
	m = len(A)
	if m == 0 {
		return 0, 0, errors.New("empty matrix")
	}
	n = len(A[0])
	if n == 0 {
		return 0, 0, errors.New("empty row")
	}
	for i := 1; i < m; i++ {
		if len(A[i]) != n {
			return 0, 0, errors.New("ragged rows: matrix must be rectangular")
		}
	}
	return m, n, nil
}

func ComputeMaxNorm(A [][]complex128, base int, mod *big.Int) *big.Float {
	size := len(A[0])

	max_int := big.NewInt(0.0)
	basepower := big.NewInt(1)

	decomposed_mod_square_m1 := DecomposeInt64Base(big.NewInt(0.0).Mul(big.NewInt(0.0).Sub(mod, big.NewInt(1.0)), big.NewInt(0.0).Sub(mod, big.NewInt(1.0))), int64(base), len(A[0]))
	idx := len(A[0]) - 1
	for i := 0; i < len(A[0]); i++ {

		if decomposed_mod_square_m1[len(A[0])-1-i] != 0 {
			break
		}
		idx -= 1
	}
	idx += 1
	fmt.Println(idx)

	max_coeff := big.NewInt(0)
	for i := size / 2; i < size; i++ {
		basetemp := big.NewInt(int64(base - 1))
		coeff := big.NewInt(0.0)
		for j := size / 2; j < idx; j++ {
			if j == idx-1 {
				entry := big.NewInt(int64(math.Round(real(A[i][j]))))
				entry.Mul(entry, big.NewInt(int64(decomposed_mod_square_m1[j])))
				coeff.Add(coeff, entry)
			} else {
				entry := big.NewInt(int64(math.Round(real(A[i][j]))))
				entry.Mul(entry, basetemp)
				coeff.Add(coeff, entry)
			}
		}
		coeff.Add(coeff, basetemp)
		if max_coeff.Cmp(coeff) == -1 {
			max_coeff = coeff
		}
		temp := big.NewInt(0.0)
		temp.Mul(coeff, basepower)
		max_int.Add(max_int, temp)
		basepower.Mul(big.NewInt(int64(base)), basepower)
	}

	max_float := IntToFloatExact(max_int)
	fmt.Println("max_int : ", max_int)
	fmt.Println("mod : ", mod)
	fmt.Println("max_digit : ", max_coeff)
	//fmt.Println("mod_square : ", big.NewInt(1).Mul(mod, mod))

	mod_float := IntToFloatExact(mod)
	rtn := max_float.Quo(max_float, mod_float)

	fmt.Println("float : ", rtn)

	decomposed_max_int := DecomposeInt64Base(max_int, int64(base), len(A[0]))
	fmt.Println(decomposed_max_int)

	idx = len(A[0]) - 1
	for i := 0; i < len(A[0]); i++ {

		if decomposed_max_int[len(A[0])-1-i] != 0 {
			break
		}
		idx -= 1
	}
	idx += 1
	fmt.Println(idx)

	max_int = big.NewInt(0.0)
	basepower = big.NewInt(1)
	max_coeff = big.NewInt(0)
	for i := size / 2; i < size; i++ {
		basetemp := big.NewInt(int64(base - 1))
		coeff := big.NewInt(0.0)
		for j := size / 2; j < idx; j++ {
			if j == idx-1 {
				entry := big.NewInt(int64(math.Round(real(A[i][j]))))
				entry.Mul(entry, big.NewInt(int64(decomposed_max_int[j])))
				coeff.Add(coeff, entry)
			} else {
				entry := big.NewInt(int64(math.Round(real(A[i][j]))))
				entry.Mul(entry, basetemp)
				coeff.Add(coeff, entry)
			}

		}
		coeff.Add(coeff, basetemp)
		if max_coeff.Cmp(coeff) == -1 {
			max_coeff = coeff
		}
		temp := big.NewInt(0.0)
		temp.Mul(coeff, basepower)
		max_int.Add(max_int, temp)
		basepower.Mul(big.NewInt(int64(base)), basepower)
	}

	max_float = IntToFloatExact(max_int)
	fmt.Println("max_int : ", max_int)
	fmt.Println("mod : ", mod)
	fmt.Println("max_digit : ", max_coeff)
	//fmt.Println("mod_square : ", big.NewInt(1).Mul(mod, mod))

	mod_float = IntToFloatExact(mod)
	rtn = max_float.Quo(max_float, mod_float)

	fmt.Println("float : ", rtn)

	decomposed_max_int = DecomposeInt64Base(max_int, int64(base), len(A[0]))
	fmt.Println(decomposed_max_int)

	idx = len(A[0]) - 1
	for i := 0; i < len(A[0]); i++ {

		if decomposed_max_int[len(A[0])-1-i] != 0 {
			break
		}
		idx -= 1
	}
	idx += 1
	fmt.Println(idx)

	max_int = big.NewInt(0.0)
	basepower = big.NewInt(1)
	max_coeff = big.NewInt(0)
	for i := size / 2; i < size; i++ {
		basetemp := big.NewInt(int64(base - 1))
		coeff := big.NewInt(0.0)
		for j := size / 2; j < idx; j++ {
			if j == idx-1 {
				entry := big.NewInt(int64(math.Round(real(A[i][j]))))
				entry.Mul(entry, big.NewInt(int64(decomposed_max_int[j])))
				coeff.Add(coeff, entry)
			} else {
				entry := big.NewInt(int64(math.Round(real(A[i][j]))))
				entry.Mul(entry, basetemp)
				coeff.Add(coeff, entry)
			}

		}
		coeff.Add(coeff, basetemp)
		if max_coeff.Cmp(coeff) == -1 {
			max_coeff = coeff
		}
		temp := big.NewInt(0.0)
		temp.Mul(coeff, basepower)
		max_int.Add(max_int, temp)
		basepower.Mul(big.NewInt(int64(base)), basepower)
	}

	max_float = IntToFloatExact(max_int)
	fmt.Println("max_int : ", max_int)
	fmt.Println("mod : ", mod)
	fmt.Println("max_digit : ", max_coeff)
	//fmt.Println("mod_square : ", big.NewInt(1).Mul(mod, mod))

	mod_float = IntToFloatExact(mod)
	rtn = max_float.Quo(max_float, mod_float)

	fmt.Println("float : ", rtn)

	return rtn
}

var magic = [6]byte{'P', 'T', 'M', 'A', 'P', '1'}

type maxLevelQer interface{ MaxLevelQ() int }
type maxLeveler interface{ MaxLevel() int }

func SavePlaintextMap(path string, m map[int]*rlwe.Plaintext) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	if _, err := w.Write(magic[:]); err != nil {
		return err
	}

	if err := binary.Write(w, binary.LittleEndian, uint64(len(m))); err != nil {
		return err
	}

	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	for _, k := range keys {
		pt := m[k]
		if pt == nil {
			return fmt.Errorf("key %d: plaintext is nil", k)
		}
		blob, err := pt.MarshalBinary()
		if err != nil {
			return fmt.Errorf("key %d: marshal failed: %w", k, err)
		}
		// key
		if err := binary.Write(w, binary.LittleEndian, int64(k)); err != nil {
			return err
		}
		// length
		if len(blob) > int(^uint32(0)) {
			return fmt.Errorf("key %d: serialized size too large", k)
		}
		if err := binary.Write(w, binary.LittleEndian, uint32(len(blob))); err != nil {
			return err
		}
		// payload
		if _, err := w.Write(blob); err != nil {
			return err
		}
	}
	return nil
}

func LoadPlaintextMap(path string, params ckks.Parameters) (map[int]*rlwe.Plaintext, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := bufio.NewReader(f)

	var got [6]byte
	if _, err := io.ReadFull(r, got[:]); err != nil {
		return nil, err
	}
	if got != magic {
		return nil, errors.New("invalid file format: bad magic")
	}

	var count uint64
	if err := binary.Read(r, binary.LittleEndian, &count); err != nil {
		return nil, err
	}

	levelMax := params.MaxLevel()

	out := make(map[int]*rlwe.Plaintext, int(count))
	for i := uint64(0); i < count; i++ {
		// key (int64)
		var key int64
		if err := binary.Read(r, binary.LittleEndian, &key); err != nil {
			return nil, fmt.Errorf("read key[%d]: %w", i, err)
		}

		// payload length (uint32)
		var n uint32
		if err := binary.Read(r, binary.LittleEndian, &n); err != nil {
			return nil, fmt.Errorf("read length for key %d: %w", key, err)
		}
		if n == 0 {
			return nil, fmt.Errorf("key %d: zero-length payload", key)
		}

		blob := make([]byte, n)
		if _, err := io.ReadFull(r, blob); err != nil {
			return nil, fmt.Errorf("read payload for key %d: %w", key, err)
		}

		var pt *rlwe.Plaintext
		var umErr error
		ok := false
		for L := levelMax; L >= 0; L-- {
			pt = ckks.NewPlaintext(params, L)
			if umErr = pt.UnmarshalBinary(blob); umErr == nil {
				ok = true
				break
			}
		}
		if !ok {
			return nil, fmt.Errorf("unmarshal key %d failed: %w", key, umErr)
		}

		out[int(key)] = pt
	}
	return out, nil
}

func Identity(n int) [][]complex128 {
	mat := make([][]complex128, n)
	for i := range mat {
		mat[i] = make([]complex128, n)
		mat[i][i] = 1
	}
	return mat
}

func isPow2(n int) bool { return n > 0 && (n&(n-1)) == 0 }

func dftMatrixScaled(n int, scale float64, sign int) [][]complex128 {
	out := make([][]complex128, n)
	w := float64(sign) * 2 * math.Pi / float64(n)
	for r := 0; r < n; r++ {
		row := make([]complex128, n)
		base := cmplx.Rect(1, w*float64(r)) // e^{iwr}
		val := complex(scale, 0)
		for c := 0; c < n; c++ {
			row[c] = val
			val *= base
		}
		out[r] = row
	}
	return out
}

func kronIK_FM(K, M int, FM [][]complex128) [][]complex128 {
	N := K * M
	out := make([][]complex128, N)
	for blk := 0; blk < K; blk++ {
		for i := 0; i < M; i++ {
			row := make([]complex128, N)
			copy(row[blk*M:blk*M+M], FM[i])
			out[blk*M+i] = row
		}
	}
	return out
}

func kronFK_IM(FK [][]complex128, K, M int) [][]complex128 {
	N := K * M
	out := make([][]complex128, N)
	for a := 0; a < K; a++ {
		for i := 0; i < M; i++ {
			row := make([]complex128, N)
			for b := 0; b < K; b++ {
				row[b*M+i] = FK[a][b]
			}
			out[a*M+i] = row
		}
	}
	return out
}

func twiddleDiagKM(K, M int) []complex128 {
	N := K * M
	out := make([]complex128, N)
	w := -2 * math.Pi / float64(N)
	for a := 0; a < K; a++ {
		for b := 0; b < M; b++ {
			out[a*M+b] = cmplx.Rect(1, w*float64(a*b))
		}
	}
	return out
}

func scaleRows(A [][]complex128, diag []complex128) {
	for r := range A {
		s := diag[r]
		row := A[r]
		for c := range row {
			row[c] *= s
		}
	}
}

func scaleCols(A [][]complex128, diag []complex128) {
	R, C := len(A), len(A[0])
	for c := 0; c < C; c++ {
		s := diag[c]
		for r := 0; r < R; r++ {
			A[r][c] *= s
		}
	}
}

func applyColPerm(A [][]complex128, new2old []int) [][]complex128 {
	R, C := len(A), len(A[0])
	B := make([][]complex128, R)
	for r := 0; r < R; r++ {
		dst := make([]complex128, C)
		src := A[r]
		for c := 0; c < C; c++ {
			dst[c] = src[new2old[c]]
		}
		B[r] = dst
	}
	return B
}

func applyRowPermLeft(A [][]complex128, new2old []int) [][]complex128 {
	R := len(A)
	B := make([][]complex128, R)
	for r := 0; r < R; r++ {
		B[r] = append([]complex128(nil), A[new2old[r]]...)
	}
	return B
}

func invertPerm(p []int) []int {
	n := len(p)
	inv := make([]int, n)
	for old, neu := range p { // p: old->new
		inv[neu] = old // inv: new->old
	}
	return inv
}

// stride Π_{K,M}: old->new  (c=aM+b → bK+a)
func strideOld2New(K, M int) []int {
	N := K * M
	p := make([]int, N)
	for a := 0; a < K; a++ {
		for b := 0; b < M; b++ {
			old := a*M + b
			new := b*K + a
			p[old] = new
		}
	}
	return p
}

func blockShiftColsOld2New(N, K, M, s int) []int {
	p := make([]int, N)
	for old := 0; old < N; old++ {
		a, i := old/M, old%M
		new := ((a+s)%K)*M + i
		p[old] = new
	}
	return p
}

func withinShiftRowsOld2New(N, K, M, t int) []int {
	p := make([]int, N)
	for old := 0; old < N; old++ {
		a, i := old/M, old%M
		new := a*M + ((i + t) % M)
		p[old] = new
	}
	return p
}

func GenerateDFT_SparseDiag(N int) (head, tail [][]complex128) {
	if !isPow2(N) {
		panic("N must be power of 2")
	}
	logN := int(math.Log2(float64(N)))
	mid := logN / 2
	M := 1 << mid
	K := N / M
	Ntot := N

	FM := dftMatrixScaled(M, 1.0, -1)
	FK := dftMatrixScaled(K, 1.0, -1)
	H0 := kronIK_FM(K, M, FM)

	A := kronFK_IM(FK, K, M)
	pi_old2new := strideOld2New(K, M)
	pi_new2old := invertPerm(pi_old2new) // for right-multiply
	T0 := applyColPerm(A, pi_new2old)    // A * Π
	D := twiddleDiagKM(K, M)
	scaleRows(T0, D) // D · (A*Π)

	Q_old2new := blockShiftColsOld2New(Ntot, K, M, 1)
	P_old2new := withinShiftRowsOld2New(Ntot, K, M, 1)
	Q_new2old := invertPerm(Q_old2new)
	P_new2old := invertPerm(P_old2new)

	H1 := applyColPerm(H0, Q_new2old)
	head = applyRowPermLeft(H1, P_new2old) // **diag(H') ≡ 0**

	// T' = Q^{-1} * T0 * P^{-1}
	T1 := applyRowPermLeft(T0, Q_old2new)
	tail = applyColPerm(T1, P_old2new) // **diag(T') has exactly K nonzeros**
	return
}

func GenerateIDFT_SparseDiag(N int) (headInv, tailInv [][]complex128) {
	if !isPow2(N) {
		panic("N must be power of 2")
	}
	logN := int(math.Log2(float64(N)))
	mid := logN / 2
	M := 1 << mid
	K := N / M
	Ntot := N

	Q_old2new := blockShiftColsOld2New(Ntot, K, M, 1)
	P_old2new := withinShiftRowsOld2New(Ntot, K, M, 1)
	Q_new2old := invertPerm(Q_old2new)
	P_new2old := invertPerm(P_old2new)

	// H^{-1} = I_K ⊗ (1/M) F_M^H
	FM_H := dftMatrixScaled(M, 1.0/float64(M), +1)
	Hinv0 := kronIK_FM(K, M, FM_H)
	// (H')^{-1} = Q^{-1} * H^{-1} * P^{-1}
	Hb := applyRowPermLeft(Hinv0, Q_old2new)
	tailInv = applyColPerm(Hb, P_old2new)

	// T^{-1} = Π^{-1} · ((1/K)F_K^H⊗I_M) · D^{-1}
	FK_H := dftMatrixScaled(K, 1.0/float64(K), +1)
	A := kronFK_IM(FK_H, K, M)
	D := twiddleDiagKM(K, M)
	for i := range D {
		D[i] = cmplx.Conj(D[i])
	} // D^{-1}=D^H
	scaleCols(A, D)
	pi_old2new := strideOld2New(K, M)
	Tinv0 := applyRowPermLeft(A, pi_old2new)

	// (T')^{-1} = P * T^{-1} * Q
	Ta := applyRowPermLeft(Tinv0, P_new2old)
	headInv = applyColPerm(Ta, Q_new2old)
	return
}

func GenerateDFT_SparseDiagUnitary(N int) (headU, tailU [][]complex128) {
	if !isPow2(N) {
		panic("N must be power of 2")
	}
	logN := int(math.Log2(float64(N)))
	mid := logN / 2
	M := 1 << mid
	K := N / M
	Ntot := N

	FM_u := dftMatrixScaled(M, 1.0/math.Sqrt(float64(M)), -1)
	FK_u := dftMatrixScaled(K, 1.0/math.Sqrt(float64(K)), -1)
	H0 := kronIK_FM(K, M, FM_u)

	A := kronFK_IM(FK_u, K, M)
	pi_old2new := strideOld2New(K, M)
	pi_new2old := invertPerm(pi_old2new)
	T0 := applyColPerm(A, pi_new2old)
	D := twiddleDiagKM(K, M)
	scaleRows(T0, D)

	Q_old2new := blockShiftColsOld2New(Ntot, K, M, 1)
	P_old2new := withinShiftRowsOld2New(Ntot, K, M, 1)
	Q_new2old := invertPerm(Q_old2new)
	P_new2old := invertPerm(P_old2new)

	H1 := applyColPerm(H0, Q_new2old)
	headU = applyRowPermLeft(H1, P_new2old)

	T1 := applyRowPermLeft(T0, Q_old2new)
	tailU = applyColPerm(T1, P_old2new)
	return
}

func GenerateIDFT_SparseDiagUnitary(N int) (headInvU, tailInvU [][]complex128) {
	if !isPow2(N) {
		panic("N must be power of 2")
	}
	logN := int(math.Log2(float64(N)))
	mid := logN / 2
	M := 1 << mid
	K := N / M
	Ntot := N

	Q_old2new := blockShiftColsOld2New(Ntot, K, M, 1)
	P_old2new := withinShiftRowsOld2New(Ntot, K, M, 1)
	Q_new2old := invertPerm(Q_old2new)
	P_new2old := invertPerm(P_old2new)

	// (Ĥ)^{-1} = I_K ⊗ (F̂_M)^H
	FM_uH := dftMatrixScaled(M, 1.0/math.Sqrt(float64(M)), +1)
	HinvU0 := kronIK_FM(K, M, FM_uH)
	Hb := applyRowPermLeft(HinvU0, Q_old2new)
	tailInvU = applyColPerm(Hb, P_old2new)

	// (T̂)^{-1} = Π^{-1} · ((F̂_K)^H⊗I_M) · D^H
	FK_uH := dftMatrixScaled(K, 1.0/math.Sqrt(float64(K)), +1)
	A := kronFK_IM(FK_uH, K, M)
	D := twiddleDiagKM(K, M)
	for i := range D {
		D[i] = cmplx.Conj(D[i])
	}
	scaleCols(A, D)
	pi_old2new := strideOld2New(K, M)
	TinvU0 := applyRowPermLeft(A, pi_old2new)

	Ta := applyRowPermLeft(TinvU0, P_new2old)
	headInvU = applyColPerm(Ta, Q_new2old)
	return
}

func GenerateTwoStageDFT(N int) ([][]complex128, [][]complex128) {
	logN := int(math.Log2(float64(N)))
	if 1<<logN != N {
		panic("N must be power of 2")
	}

	mid := logN / 2
	//fmt.Println(mid)
	//mid = logN

	// Head = B_0 to B_{mid-1}
	head := Identity(N)
	for i := 0; i < mid; i++ {
		B := GenerateButterflyStage(N, i)
		head = MatMul(B, head)
	}

	// Tail = B_{mid} to B_{logN-1}
	tail := Identity(N)
	for i := mid; i < logN; i++ {
		B := GenerateButterflyStage(N, i)
		tail = MatMul(B, tail)
	}

	return head, tail
}

func GenerateSpecialNormalizedDFT(N int) [][]complex128 {
	matrix := make([][]complex128, N)
	scale := 1.0 / math.Sqrt(float64(N))
	twoPi := 2 * math.Pi
	omega := complex(math.Cos(-twoPi/float64(N)), math.Sin(-twoPi/float64(N))) // e^{-2πi/N}

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			exp := float64(j * k)
			matrix[j][k] = cmplx.Pow(omega, complex(exp, 0)) * complex(scale, 0)
		}
	}
	return matrix
}

func GenerateNormalizedIDFT(N int) [][]complex128 {
	matrix := make([][]complex128, N)
	scale := 1.0 / math.Sqrt(float64(N))
	//twoPi := 2 * math.Pi
	omegaInv := complex(math.Cos(2*math.Pi/float64(N)), math.Sin(2*math.Pi/float64(N))) // e^{2πi/N}

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			exp := float64(j * k)
			matrix[j][k] = cmplx.Pow(omegaInv, complex(exp, 0)) * complex(scale, 0)
		}
	}
	return matrix
}

func GenerateDFT(N int) [][]complex128 {
	matrix := make([][]complex128, N)
	scale := 1.0 // math.Sqrt(float64(N))
	twoPi := 2 * math.Pi
	omega := complex(math.Cos(-twoPi/float64(N)), math.Sin(-twoPi/float64(N))) // e^{-2πi/N}

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			exp := float64(j * k)
			matrix[j][k] = cmplx.Pow(omega, complex(exp, 0)) * complex(scale, 0)
		}
	}
	return matrix
}

func GenerateNormalizedInvDFT(N int) [][]complex128 {
	matrix := make([][]complex128, N)
	scale := 1.0 // math.Sqrt(float64(N))
	//twoPi := 2 * math.Pi
	omegaInv := complex(math.Cos(2*math.Pi/float64(N)), math.Sin(2*math.Pi/float64(N))) // e^{2πi/N}

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			exp := float64(j * k)
			matrix[j][k] = cmplx.Pow(omegaInv, complex(exp, 0)) * complex(scale, 0)
		}
	}

	return matrix
}

func GenerateNormalizedInvDFT_with_masking(N int) [][]complex128 {
	matrix := make([][]complex128, N)
	scale := 1.0 // math.Sqrt(float64(N))
	//twoPi := 2 * math.Pi
	omegaInv := complex(math.Cos(2*math.Pi/float64(N)), math.Sin(2*math.Pi/float64(N))) // e^{2πi/N}

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			exp := float64(j * k)
			matrix[j][k] = cmplx.Pow(omegaInv, complex(exp, 0)) * complex(scale, 0)
		}
	}

	for j := N / 2; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			matrix[j][k] = complex(0.0, 0.0)
		}
	}

	return matrix
}

func GenerateInvDFT(N int) [][]complex128 {
	matrix := make([][]complex128, N)
	scale := 1.0 // math.Sqrt(float64(N))
	//twoPi := 2 * math.Pi
	omegaInv := complex(math.Cos(2*math.Pi/float64(N)), math.Sin(2*math.Pi/float64(N))) // e^{2πi/N}

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
		for k := 0; k < N; k++ {
			exp := float64(j * k)
			matrix[j][k] = cmplx.Pow(omegaInv, complex(exp, 0)) * complex(scale, 0)
		}
	}
	return matrix
}

func MatrixMult(A [][]complex128, B [][]complex128) [][]complex128 {
	N := len(B[0])
	matrix := make([][]complex128, N)

	for j := 0; j < N; j++ {
		matrix[j] = make([]complex128, N)
	}

	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			for k := 0; k < N; k++ {
				matrix[i][j] += A[i][k] * B[k][j]
			}
		}
	}

	return matrix
}

func LinTrans(A [][]complex128, x []complex128) []complex128 {
	N := len(A[0])
	b := make([]complex128, N)

	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			b[i] += A[i][j] * x[j]
		}
	}

	return b
}

func HadmardMult(x []complex128, y []complex128) []complex128 {
	N := len(x)
	b := make([]complex128, N)

	for i := 0; i < N; i++ {
		b[i] = x[i] * y[i]
	}
	return b
}

func PrintMatrix(mat [][]complex128, name string) {
	fmt.Printf("=== %s ===\n", name)

	//for i := 0; i < 10; i++ {
	//	fmt.Print(mat[0][i])
	//}

	for _, row := range mat {
		for _, v := range row {
			fmt.Printf("%6.3f + %6.3fi\t", real(v), imag(v))
		}
		fmt.Println()
	}
	fmt.Println()
}

func PrintVec(vec []complex128, name string) {
	fmt.Printf("=== %s ===\n", name)
	for i := 0; i < len(vec); i++ {
		fmt.Printf("%6.3f + %6.3fi\t", real(vec[i]), imag(vec[i]))
	}
	fmt.Println()
}

func ConcatMatrix(A [][]complex128, n int, N int) [][]complex128 {
	ConcatA := make([][]complex128, N)

	for i := 0; i < N; i++ {
		ConcatA[i] = make([]complex128, N)
	}

	for k := 0; k < N/n; k++ {
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				ConcatA[i+n*k][j+n*k] = A[i][j]
			}
		}
	}

	return ConcatA
}

func TwistedMatrix(A [][]complex128, n int, N int) [][]complex128 {
	batch := N / n
	twistedA := make([][]complex128, N)

	for i := 0; i < N; i++ {
		twistedA[i] = make([]complex128, N)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < batch; k++ {
				twistedA[i*batch+k][j*batch+k] = A[i][j]
			}
		}
	}

	return twistedA
}

func MatrixPadding(A [][]complex128, n int, N int) [][]complex128 {
	doubleA := make([][]complex128, n)
	for i := 0; i < n; i++ {
		doubleA[i] = make([]complex128, n)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i < len(A[0]) && j < len(A[0]) {
				doubleA[i][j] = A[i][j]
			} else {
				doubleA[i][j] = complex(0, 0)
			}
		}
	}

	return doubleA
}

func TwistedVec(vec []complex128, n int, N int) []complex128 {
	batch := N / n
	twistedVec := make([]complex128, N)

	for i := 0; i < n; i++ {
		for j := 0; j < batch; j++ {
			twistedVec[i*batch+j] = vec[j*n+i]
		}
	}

	return twistedVec
}

func InvTwistedVec(vec []complex128, n int, N int) []complex128 {
	batch := N / n
	InvtwistedVec := make([]complex128, N)

	for i := 0; i < n; i++ {
		for j := 0; j < batch; j++ {
			InvtwistedVec[j*n+i] = vec[i*batch+j]
		}
	}

	return InvtwistedVec
}

func Rotation(z []complex128, r int) []complex128 {
	m := len(z)
	rotz := make([]complex128, m)
	r += m // Handling negative r
	for i := 0; i < m; i++ {
		rotz[i] = z[(i+r)%m]
	}
	return rotz
}

func Row_To_Diagonal(A [][]complex128) [][]complex128 {
	N := len(A[0])
	diag := make([][]complex128, N)

	for i := 0; i < N; i++ {
		diag[i] = make([]complex128, N)
	}

	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			diag[i][j] = A[j][(i+j)%N]
		}
	}

	cnt := 0
	valid_idx := []int{}
	for i := 0; i < N; i++ {
		sum := complex(0.0, 0.0)
		for j := 0; j < N; j++ {
			sum += diag[i][j]
		}

		if sum != complex(0.0, 0.0) {
			cnt += 1
			valid_idx = append(valid_idx, i)
		}
	}

	//fmt.Println(cnt)
	//fmt.Println(valid_idx)

	return diag
}

// ilog2Pow2 returns log2(N) when N is a power of two; panics otherwise.
func ilog2Pow2(N int) int {
	if N <= 0 || (N&(N-1)) != 0 {
		panic("N must be a positive power of two")
	}
	// compute k such that N == 1<<k
	k := 0
	for (1 << k) != N {
		k++
	}
	return k
}

// BitReverseMatrixDense returns an N x N 0/1 matrix P such that y = P*x
// applies bit-reversal permutation: y[brev(i)] = x[i].
func BitReverseMatrixDense(N int) [][]complex128 {
	k := ilog2Pow2(N)
	P := make([][]complex128, N)
	for r := 0; r < N; r++ {
		P[r] = make([]complex128, N)
	}
	for i := 0; i < N; i++ {
		j := int(reverseBits(uint(i), k))
		P[j][i] = 1
	}
	return P
}

func Transpose(A [][]complex128) [][]complex128 {
	N := len(A[0])
	AT := make([][]complex128, N)
	for i := 0; i < N; i++ {
		AT[i] = make([]complex128, N)
	}
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			AT[i][j] = A[j][i]
		}
	}

	return AT

}

// COO is a simple sparse coordinate list: ones at (Rows[t], Cols[t]).
type COO struct {
	N    int
	Rows []int
	Cols []int
}

// BitReverseMatrixCOO builds the same permutation in sparse-COO form.
// Exactly N ones: for each column i, a single row j=brev(i).
func BitReverseMatrixCOO(N int) *COO {
	k := ilog2Pow2(N)
	rows := make([]int, 0, N)
	cols := make([]int, 0, N)
	for i := 0; i < N; i++ {
		j := int(reverseBits(uint(i), k))
		rows = append(rows, j)
		cols = append(cols, i)
	}
	return &COO{N: N, Rows: rows, Cols: cols}
}

// ApplyDense multiplies P (0/1) by x.
func ApplyDense(P [][]int, x []float64) []float64 {
	N := len(P)
	if len(x) != N {
		panic("dimension mismatch")
	}
	y := make([]float64, N)
	for r := 0; r < N; r++ {
		row := P[r]
		var acc float64
		for c := 0; c < N; c++ {
			if row[c] != 0 {
				acc += x[c] // row has a single 1, but this is general
			}
		}
		y[r] = acc
	}
	return y
}

// ApplyCOO multiplies the COO permutation by x in O(N).
func ApplyCOO(P *COO, x []float64) []float64 {
	if len(x) != P.N {
		panic("dimension mismatch")
	}
	y := make([]float64, P.N)
	for t := 0; t < P.N; t++ {
		r := P.Rows[t]
		c := P.Cols[t]
		y[r] = y[r] + x[c] // exactly one 1 per column → simple move
	}
	return y
}

func MaxAbsDiff(A, B [][]complex128) float64 {
	n, m := len(A), len(A[0])
	if n != len(B) || m != len(B[0]) {
		panic("dimension mismatch")
	}
	maxd := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			d := cmplx.Abs(A[i][j] - B[i][j])
			if d > maxd {
				maxd = d
			}
		}
	}
	return maxd
}

func reverseBits(x uint, k int) uint {
	var r uint
	for i := 0; i < k; i++ {
		r = (r << 1) | (x & 1)
		x >>= 1
	}
	return r
}

func BitReverseMatrix(N int) [][]complex128 {
	k := ilog2Pow2(N)
	P := make([][]complex128, N)
	for r := 0; r < N; r++ {
		P[r] = make([]complex128, N)
	}
	for i := 0; i < N; i++ {
		j := int(reverseBits(uint(i), k))
		P[j][i] = 1
	}
	return P
}

func GenerateStageMatrix(N, q int) [][]complex128 {
	if !isPow2(N) {
		panic("N must be power of two")
	}
	m := ilog2Pow2(N)
	if q < 1 || q > m {
		panic("q must be in 1..log2(N)")
	}

	L := 1 << q
	r := N / L
	A := make([][]complex128, N)
	for i := 0; i < N; i++ {
		A[i] = make([]complex128, N)
	}

	// twiddle base
	w := complex(math.Cos(-2*math.Pi/float64(L)), math.Sin(-2*math.Pi/float64(L)))

	for b := 0; b < r; b++ {
		start := b * L
		for t := 0; t < L/2; t++ {
			d := cmplx.Pow(w, complex(float64(t), 0)) // w^t

			r0 := start + t
			r1 := start + t + L/2
			c0 := start + t
			c1 := start + t + L/2

			// [ [I,  D],
			//   [I, -D] ]
			A[r0][c0] += 1
			A[r0][c1] += d
			A[r1][c0] += 1
			A[r1][c1] += -d
		}
	}
	return A
}

type FFTFactorization struct {
	P      [][]complex128   // bit-reversal
	Stages [][][]complex128 // A_1, A_2, ..., A_m
}

func FactorizeDFT(N int) FFTFactorization {
	m := ilog2Pow2(N)
	P := BitReverseMatrix(N)
	stages := make([][][]complex128, 0, m)
	for q := 1; q <= m; q++ {
		Aq := GenerateStageMatrix(N, q)
		stages = append(stages, Aq)
	}
	return FFTFactorization{P: P, Stages: stages}
}

func ExpandFactors(f FFTFactorization) [][]complex128 {
	N := len(f.P)
	M := Identity(N)
	M = MatMul(f.P, M)
	for _, Aq := range f.Stages {
		M = MatMul(Aq, M)
	}
	return M
}

func HeadTail(f FFTFactorization, s int) (head, tail [][]complex128) {
	if s < 0 || s > len(f.Stages) {
		panic("0 <= s <= m")
	}
	N := len(f.P)
	head = Identity(N)
	for i := 0; i < s; i++ {
		head = MatMul(f.Stages[i], head) // A_{i+1} ... A_1
	}
	tail = Identity(N)
	for i := len(f.Stages) - 1; i >= s; i-- {
		tail = MatMul(f.Stages[i], tail) // A_m ... A_{s+1}
	}
	return
}
func Inverse(A [][]complex128) [][]complex128 {
	n := len(A)
	if n == 0 || len(A[0]) != n {
		panic("Inverse: square matrix required")
	}
	M := make([][]complex128, n)
	for i := 0; i < n; i++ {
		M[i] = make([]complex128, 2*n)
		copy(M[i][:n], A[i])
		M[i][n+i] = 1
	}
	const eps = 1e-15
	for col := 0; col < n; col++ {
		pivot := col
		maxAbs := 0.0
		for r := col; r < n; r++ {
			v := cmplx.Abs(M[r][col])
			if v > maxAbs {
				maxAbs = v
				pivot = r
			}
		}
		if maxAbs < eps {
			panic("Inverse: singular or ill-conditioned matrix (pivot ~ 0)")
		}
		if pivot != col {
			M[col], M[pivot] = M[pivot], M[col]
		}
		p := M[col][col]
		for j := 0; j < 2*n; j++ {
			M[col][j] /= p
		}
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := M[r][col]
			if f == 0 {
				continue
			}
			for j := 0; j < 2*n; j++ {
				M[r][j] -= f * M[col][j]
			}
		}
	}
	Inv := make([][]complex128, n)
	for i := 0; i < n; i++ {
		Inv[i] = make([]complex128, n)
		copy(Inv[i], M[i][n:])
	}
	return Inv
}
