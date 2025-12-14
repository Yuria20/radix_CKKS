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

	// PolyMult for Bk arithmetic!
	Test_Bk_bit_length := []int{16, 32, 64, 128, 256, 512, 1024, 2048, 4096}
	for i := 0; i < 8; i++ {
		// Bbox Init
		var bit_length, base int
		bit_length = Test_Bk_bit_length[i]
		if bit_length == 16 {
			base = 16
		} else {
			base = 256
		}

		var Bbox Bk_Arithmetic_toolbox
		Bbox.Init(bit_length, base)
		for j := 0; j < 1; j++ {
			fmt.Println("=======================================================")
			fmt.Printf("============ Bk arithmetic Test : %d bit ============\n", bit_length)
			fmt.Println("=======================================================")
			Test_Bk_arthmetic(Bbox)
		}
	}

	return

	Test_Bk_bit_length = []int{512, 1024, 4096}
	for i := 0; i < 3; i++ {
		// Bbox Init
		var bit_length, base int
		bit_length = Test_Bk_bit_length[i]
		if bit_length == 16 {
			base = 16
		} else {
			base = 256
		}

		var Bbox Bk_Arithmetic_toolbox
		//Bbox.Init(bit_length, base)
		Bbox.InitModular(bit_length, base)
		for j := 0; j < 1; j++ {
			// Test Start
			if i == 0 {
				Test_curve25519_arthmetic(Bbox)
			} else if i == 1 {
				Test_P384_baseREDC(Bbox)
			} else {
				Test_montREDC(Bbox)
			}
		}
	}

	return

	// Modular Aritmetic!

	//fmt.Println("=======================================================")
	//fmt.Printf("=========== Zp arithmetic Test : Curve25519 ===========\n")
	//fmt.Println("=======================================================")
	var Bbox Bk_Arithmetic_toolbox
	Bbox.InitModular(512, 16)
	Test_curve25519_arthmetic(Bbox)

	fmt.Println("=======================================================")
	fmt.Printf("=========== Zp arithmetic Test : RSA2048 ===========\n")
	fmt.Println("=======================================================")
	//var Bbox Bk_Arithmetic_toolbox
	//Bbox.InitRSA(2048, 16)
	Test_montREDC(Bbox)

	return

}

type Bk_Arithmetic_toolbox struct {
	bit_length     int
	base           int
	base_bit       int
	slice_length   int
	log_slice_half int

	lazy_iter  map[int]int
	carry_iter map[int]int
}

func (Bbox *Bk_Arithmetic_toolbox) Init(bit_length int, base int) {

	Bbox.bit_length = bit_length
	if bit_length == 16 {
		Bbox.base_bit = 4
	} else {
		Bbox.base_bit = 8
	}

	Bbox.base_bit = 4

	Bbox.base = 1 << Bbox.base_bit

	Bbox.lazy_iter = make(map[int]int, 9)
	Bbox.carry_iter = make(map[int]int, 9)

	if Bbox.base_bit == 4 {
		Bbox.lazy_iter[16] = 2
		Bbox.lazy_iter[32] = 2
		Bbox.lazy_iter[64] = 2
		Bbox.lazy_iter[128] = 3
		Bbox.lazy_iter[256] = 3
		Bbox.lazy_iter[512] = 3
		Bbox.lazy_iter[1024] = 3
		Bbox.lazy_iter[2048] = 4
		Bbox.lazy_iter[4096] = 4

		Bbox.carry_iter[16] = 3
		Bbox.carry_iter[32] = 4
		Bbox.carry_iter[64] = 5
		Bbox.carry_iter[128] = 6
		Bbox.carry_iter[256] = 7
		Bbox.carry_iter[512] = 8
		Bbox.carry_iter[1024] = 9
		Bbox.carry_iter[2048] = 10
		Bbox.carry_iter[4096] = 11
	}

	if Bbox.base_bit == 8 {
		Bbox.lazy_iter[16] = 2
		Bbox.lazy_iter[32] = 2
		Bbox.lazy_iter[64] = 2
		Bbox.lazy_iter[128] = 3
		Bbox.lazy_iter[256] = 3
		Bbox.lazy_iter[512] = 3
		Bbox.lazy_iter[1024] = 3
		Bbox.lazy_iter[2048] = 4
		Bbox.lazy_iter[4096] = 4

		Bbox.carry_iter[16] = 3
		Bbox.carry_iter[32] = 4
		Bbox.carry_iter[64] = 5
		Bbox.carry_iter[128] = 6
		Bbox.carry_iter[256] = 7
		Bbox.carry_iter[512] = 8
		Bbox.carry_iter[1024] = 9
		Bbox.carry_iter[2048] = 10
		Bbox.carry_iter[4096] = 11
	}

	Bbox.slice_length = 2 * (Bbox.bit_length / Bbox.base_bit)
	Bbox.log_slice_half = int(math.Log2(float64((Bbox.bit_length / Bbox.base_bit))))
}

func (Bbox *Bk_Arithmetic_toolbox) InitModular(bit_length int, base int) {

	Bbox.bit_length = bit_length
	if bit_length == 16 {
		Bbox.base_bit = 4
	} else {
		Bbox.base_bit = 8
	}

	Bbox.base_bit = 4

	Bbox.base = 1 << Bbox.base_bit

	Bbox.lazy_iter = make(map[int]int, 9)
	Bbox.carry_iter = make(map[int]int, 9)

	Bbox.lazy_iter[16] = 2
	Bbox.lazy_iter[32] = 2
	Bbox.lazy_iter[64] = 2
	Bbox.lazy_iter[128] = 3
	Bbox.lazy_iter[256] = 3
	Bbox.lazy_iter[512] = 3
	Bbox.lazy_iter[1024] = 3
	Bbox.lazy_iter[2048] = 4
	Bbox.lazy_iter[4096] = 4

	Bbox.carry_iter[16] = 3
	Bbox.carry_iter[32] = 4
	Bbox.carry_iter[64] = 5
	Bbox.carry_iter[128] = 6
	Bbox.carry_iter[256] = 7
	Bbox.carry_iter[512] = 8
	Bbox.carry_iter[1024] = 9
	Bbox.carry_iter[2048] = 10
	Bbox.carry_iter[4096] = 11

	Bbox.slice_length = 2 * (Bbox.bit_length / Bbox.base_bit)
	Bbox.log_slice_half = int(math.Log2(float64((Bbox.bit_length / Bbox.base_bit))))
}

func (Bbox *Bk_Arithmetic_toolbox) InitCurve25519(base int) {

	Bbox.bit_length = 1024
	Bbox.base_bit = 5
	Bbox.base = 1 << Bbox.base_bit

	Bbox.lazy_iter = make(map[int]int, 6)
	Bbox.carry_iter = make(map[int]int, 6)

	Bbox.lazy_iter[1024] = 3
	Bbox.carry_iter[1024] = 8

	Bbox.slice_length = 256
	Bbox.log_slice_half = 7
}

type Gmn_Arithmetic_toolbox struct {
	bit_length     int
	base           int
	base_bit       int
	slice_length   int
	log_slice_half int

	lazy_iter  map[int]int
	carry_iter map[int]int
}

func (Gbox *Gmn_Arithmetic_toolbox) Init(bit_length int, base int) {

	Gbox.bit_length = bit_length
	Gbox.base_bit = 8
	Gbox.base = 1 << Gbox.base_bit

	Gbox.lazy_iter = make(map[int]int, 6)
	Gbox.carry_iter = make(map[int]int, 6)

	Gbox.lazy_iter[64] = 2
	Gbox.lazy_iter[128] = 2
	Gbox.lazy_iter[256] = 2
	Gbox.lazy_iter[512] = 2
	Gbox.lazy_iter[1024] = 2
	Gbox.lazy_iter[2048] = 3

	Gbox.carry_iter[64] = 4
	Gbox.carry_iter[128] = 5
	Gbox.carry_iter[256] = 6
	Gbox.carry_iter[512] = 7
	Gbox.carry_iter[1024] = 8
	Gbox.carry_iter[2048] = 9

	Gbox.slice_length = 2 * (Gbox.bit_length / Gbox.base_bit)
	Gbox.log_slice_half = int(math.Log2(float64((Gbox.bit_length / Gbox.base_bit))))
}

type Context struct {
	params    *ckks.Parameters
	encoder   *ckks.Encoder
	encryptor *rlwe.Encryptor
	decryptor *rlwe.Decryptor
	eval      *ckks.Evaluator
	btp       *bootstrapping.Evaluator
}
