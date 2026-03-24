// Package main implements an example showcasing slim for bootstrapping for fixed-point approximate
// arithmetic over the reals/complexes numbers.
// This re-ordering of the bootstrapping steps was first proposed for the BFV/BGV schemes by Chen and Han
// in Homomorphic Lower Digits Removal and Improved FHE Bootstrapping (https://eprint.iacr.org/2018/067).
// It was also used by Kim and Guyot in Optimized Privacy-Preserving CNN Inference With Fully Homomorphic
// Encryption (https://ieeexplore.ieee.org/document/10089847) to efficiently perform the convolution in
// the coefficient domain.
//
// This example assumes that the user is already familiar with the bootstrapping and its different steps.
// See the basic example `lattigo/single_party/applications/reals_bootstrapping/basics` for an introduction into the
// bootstrapping.
//
// The usual order of the bootstrapping operations is:
//
// 0) User defined circuit in the slots domain
// 1) ScaleDown: Scale the ciphertext to q0/|m|
// 2) ModUp: Raise modulus from q0 to qL
// 3) CoeffsToSlots: Homomorphic encoding
// 4) EvalMod: Homomorphic modular reduction
// 5) SlotsToCoeffs (and go back to 0): Homomorphic Decoding
//
// This example instantiates a custom order of the circuit evaluating:
//
// 0) User defined circuit in the slots domain
// 1) SlotsToCoeffs: Homomorphic Decoding
// 2) User defined circuit in the coeffs domain
// 3) ScaleDown: Scale the ciphertext to q0/|m|
// 4) ModUp: Raise modulus from q0 to qL
// 5) CoeffsToSlots: Homomorphic encoding
// 6) EvalMod (and to back to 0): Homomorphic modular reduction
//
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

	// Set setting = false to pre-generate the required plaintexts.
	// After that, set setting = true and run the benchmark/experiments.
	setting := false

	// PolyMult for Bk arithmetic!
	Test_Bk_bit_length := []int{16, 32, 64, 128, 256, 512, 1024, 2048, 4096}
	for i := 0; i < 8; i++ {
		// Bbox Init
		var bit_length, logbase int
		bit_length = Test_Bk_bit_length[i]
		logbase = 4

		var Bbox RParams
		Bbox.Init(false, bit_length, logbase)

		fmt.Println("=======================================================")
		fmt.Printf("============ Bk arithmetic Test : %d bit ============\n", bit_length)
		fmt.Println("=======================================================")
		Test_Bk_arthmetic_B16(Bbox, setting)
	}

	Test_Bk_bit_length = []int{256, 512, 2048}
	for i := 0; i < 3; i++ {
		// Bbox Init
		var bit_length, logbase int
		bit_length = Test_Bk_bit_length[i]
		logbase = 4

		var radix_params RParams
		radix_params.Init(true, bit_length, logbase)

		// Test Start
		if i == 0 {
			fmt.Println("=======================================================")
			fmt.Printf("=============== Curve25519 arithmetic Test ==============\n")
			fmt.Println("=======================================================")
			Test_curve25519_arthmetic(radix_params, setting)
		} else if i == 1 {
			fmt.Println("=======================================================")
			fmt.Printf("================ P-384 arithmetic Test== ==============\n")
			fmt.Println("=======================================================")
			Test_P384_baseREDC(radix_params, setting)
		} else {
			fmt.Println("=======================================================")
			fmt.Printf("============== Montogmery arithmetic Test ==============\n")
			fmt.Println("=======================================================")
			Test_montREDC(radix_params, setting)
		}
	}

	return

}

type RParams struct {
	// type
	ismod bool

	// radix parameter
	bit_length   int
	B            int
	logB         int
	K            int
	slice_length int
	logK         int

	lazy_iter int
}

func (radix_params *RParams) Init(ismod bool, bit_length int, logB int) {

	radix_params.ismod = ismod

	if radix_params.ismod == true {
		radix_params.bit_length = 2 * bit_length
	} else {
		radix_params.bit_length = bit_length
	}

	radix_params.logB = logB
	radix_params.B = 1 << logB

	radix_params.K = radix_params.bit_length / radix_params.logB
	radix_params.slice_length = 2 * radix_params.K
	radix_params.logK = int(math.Log2(float64(radix_params.K)))

	compute_lazycarry_iter := func(K int, B int, ismod bool) int {

		if ismod == true {
			K = K / 2
		}

		max_digit := B - 1
		output_digit := K * max_digit * max_digit
		iter := 0
		for output_digit >= 2*B-1 {
			output_digit = output_digit/B + max_digit
			iter += 1
		}

		return iter
	}

	radix_params.lazy_iter = compute_lazycarry_iter(radix_params.K, radix_params.B, ismod)
}

type Context struct {
	params    *ckks.Parameters
	encoder   *ckks.Encoder
	encryptor *rlwe.Encryptor
	decryptor *rlwe.Decryptor
	eval      *ckks.Evaluator
	btp       *bootstrapping.Evaluator
}
