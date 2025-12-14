package main

import (
	"fmt"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func LazyMult(ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext, Plain_DFT map[int]*rlwe.Plaintext, Plain_InvDFT map[int]*rlwe.Plaintext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {

	bit_length := Bbox.bit_length
	base := Bbox.base

	slice_length := Bbox.slice_length

	values := make([]complex128, cc.params.MaxSlots())
	var ciphertext *rlwe.Ciphertext

	// LazyCarry

	//fmt.Print("PolyMult... ")
	//PrintDebug(slice_length, *cc.params, ciphertext1, values, cc.decryptor, cc.encoder)
	ciphertext = PolyMult(ciphertext1, ciphertext2, Plain_DFT, Plain_InvDFT, Bbox, cc)
	PrintDebug(slice_length, *cc.params, ciphertext, values, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *cc.params, ciphertext, values, cc.decryptor, cc.encoder)

	for iter := 0; iter < Bbox.lazy_iter[bit_length]; iter++ {
		fmt.Print("LazyCarry... ")
		start := time.Now()
		ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
		PrintDebug(slice_length, *cc.params, ciphertext, values, cc.decryptor, cc.encoder)
		elapsed_time := time.Since(start)
		fmt.Println(elapsed_time)
	}
	return ciphertext
}

func LazyCMult(ciphertext1 *rlwe.Ciphertext, DFT_ptxt *rlwe.Plaintext, Plain_DFT map[int]*rlwe.Plaintext, Plain_InvDFT map[int]*rlwe.Plaintext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {

	bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	//log_slice_half := Bbox.log_slice_half

	var ciphertext *rlwe.Ciphertext

	// LazyCarry
	ciphertext = PolyCMult(ciphertext1, DFT_ptxt, Plain_DFT, Plain_InvDFT, Bbox, cc)

	for iter := 0; iter < Bbox.lazy_iter[bit_length]; iter++ {
		ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
	}
	return ciphertext
}

func LazyCarry2Carry(ciphertext *rlwe.Ciphertext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {
	//bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	log_slice_half := Bbox.log_slice_half

	var ct_symbol, ct_carry *rlwe.Ciphertext
	ct_symbol = ciphertext.CopyNew()
	ct_carry = ciphertext.CopyNew()

	// LazyCarry To Carry
	fmt.Print("SI Boot... ")
	start := time.Now()
	ct_symbol = IntToSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
	elapsed_time := time.Since(start)
	fmt.Println(elapsed_time)

	//values := make([]complex128, cc.params.MaxSlots())
	//PrintDebug(slice_length, *cc.params, ciphertext, values, cc.decryptor, cc.encoder)
	start = time.Now()

	ct_carry = LazyCarryToCarry(base, log_slice_half, slice_length, ct_symbol, ct_carry, cc)

	ciphertext = ct_carry.CopyNew()

	return ciphertext
}

func ExactCMult(ciphertext1 *rlwe.Ciphertext, DFT_ptxt *rlwe.Plaintext, Plain_DFT map[int]*rlwe.Plaintext, Plain_InvDFT map[int]*rlwe.Plaintext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {

	bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	log_slice_half := Bbox.log_slice_half

	var ciphertext *rlwe.Ciphertext

	// LazyCarry
	ciphertext = PolyCMult(ciphertext1, DFT_ptxt, Plain_DFT, Plain_InvDFT, Bbox, cc)

	for iter := 0; iter < Bbox.lazy_iter[bit_length]; iter++ {
		fmt.Print("LazyCarry... ")
		start := time.Now()
		ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
		elapsed_time := time.Since(start)
		fmt.Println(elapsed_time)
	}

	var ct_symbol, ct_carry *rlwe.Ciphertext
	ct_symbol = ciphertext.CopyNew()
	ct_carry = ciphertext.CopyNew()

	// LazyCarry To Carry
	fmt.Print("SI Boot... ")
	start := time.Now()
	ct_symbol = IntToSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
	elapsed_time := time.Since(start)
	fmt.Println(elapsed_time)

	ct_carry = LazyCarryToCarry(base, log_slice_half, slice_length, ct_symbol, ct_carry, cc)
	ciphertext = ct_carry.CopyNew()

	return ciphertext
}

func ExactMult(ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext, Plain_DFT map[int]*rlwe.Plaintext, Plain_InvDFT map[int]*rlwe.Plaintext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {

	bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	log_slice_half := Bbox.log_slice_half

	var ciphertext *rlwe.Ciphertext

	// LazyCarry

	ciphertext = PolyMult(ciphertext1, ciphertext2, Plain_DFT, Plain_InvDFT, Bbox, cc)

	for iter := 0; iter < Bbox.lazy_iter[bit_length]; iter++ {
		fmt.Print("LazyCarry... ")
		start := time.Now()
		ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
		elapsed_time := time.Since(start)
		fmt.Println(elapsed_time)
	}

	var ct_symbol, ct_carry *rlwe.Ciphertext
	ct_symbol = ciphertext.CopyNew()
	ct_carry = ciphertext.CopyNew()

	// LazyCarry To Carry
	fmt.Print("SI Boot... ")
	start := time.Now()
	ct_symbol = IntToSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
	elapsed_time := time.Since(start)
	fmt.Println(elapsed_time)

	//fmt.Print("Carry... ")
	//start = time.Now()
	ct_carry = LazyCarryToCarry(base, log_slice_half, slice_length, ct_symbol, ct_carry, cc)
	//elapsed_time = time.Since(start)
	//fmt.Println(elapsed_time)
	ciphertext = ct_carry.CopyNew()

	return ciphertext
}

func Comparison(ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {

	//bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	log_slice_half := Bbox.log_slice_half
	MAX_SLOT := cc.params.MaxSlots()
	BATCH := MAX_SLOT / slice_length
	var ciphertext *rlwe.Ciphertext

	ciphertext, _ = cc.eval.SubNew(ciphertext1, ciphertext2)

	masking := make([]complex128, MAX_SLOT)
	for i := 0; i < BATCH; i++ {
		for j := 0; j < slice_length; j++ {
			if j >= slice_length/2 {
				masking[i*slice_length+j] = complex(1.0, 0.0)
			} else {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			}
		}
	}
	twisted_masking := TwistedVec(masking, slice_length, MAX_SLOT)

	cc.eval.Add(ciphertext, twisted_masking, ciphertext)

	var ct_symbol, ct_carry *rlwe.Ciphertext
	ct_symbol = ciphertext.CopyNew()
	ct_carry = ciphertext.CopyNew()

	ct_symbol = IntToSymbol_negate_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)

	//fmt.Print("Carry Start...")
	ct_carry = LazyCarryToCarry_negate(base, log_slice_half, slice_length, ct_symbol, ct_carry, cc)

	masking = make([]complex128, MAX_SLOT)
	for i := 0; i < BATCH; i++ {
		for j := 0; j < slice_length; j++ {
			if j == slice_length/2 {
				masking[i*slice_length+j] = complex(1.0, 0.0)
			} else {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			}
		}
	}
	twisted_masking = TwistedVec(masking, slice_length, MAX_SLOT)
	cc.eval.MulRelin(ct_carry, twisted_masking, ct_carry)
	cc.eval.Rescale(ct_carry, ct_carry)
	cc.eval.Rotate(ct_carry, MAX_SLOT/2, ct_carry)

	return ct_carry
}

func ConditonalSub(ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext, Bbox Bk_Arithmetic_toolbox, cc Context, tag string) *rlwe.Ciphertext {

	//bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	log_slice_half := Bbox.log_slice_half
	MAX_SLOT := cc.params.MaxSlots()
	BATCH := MAX_SLOT / slice_length
	var ciphertext *rlwe.Ciphertext

	ciphertext, _ = cc.eval.SubNew(ciphertext1, ciphertext2)

	masking := make([]complex128, MAX_SLOT)
	for i := 0; i < BATCH; i++ {
		for j := 0; j < slice_length; j++ {
			if j == slice_length/2 {
				masking[i*slice_length+j] = complex(1.0, 0.0)
			} else {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			}
		}
	}
	twisted_masking := TwistedVec(masking, slice_length, MAX_SLOT)

	cc.eval.Add(ciphertext, twisted_masking, ciphertext)

	var ct_symbol, ct_carry *rlwe.Ciphertext
	ct_symbol = ciphertext.CopyNew()
	ct_carry = ciphertext.CopyNew()

	fmt.Print("SIBoot... ")
	start := time.Now()
	//ct_symbol = IntToSymbol_negate_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
	ct_symbol = IntToSymbol_negate(base, slice_length, ct_symbol, cc)
	//PrintDebug(slice_length, *cc.params, ct_symbol, twisted_masking, cc.decryptor, cc.encoder)
	elapsed_time := time.Since(start)
	fmt.Println(elapsed_time)

	ct_carry = LazyCarryToCarry_negate(base, log_slice_half, slice_length, ct_symbol, ct_carry, cc)
	//PrintDebug(slice_length, *cc.params, ct_carry, twisted_masking, cc.decryptor, cc.encoder)
	fmt.Print("Cleaning... ")
	start = time.Now()
	//ct_carry = Cleaning(ct_carry, Bbox, cc)
	ciphertext = ct_carry.CopyNew()

	// using modular arithmetic
	//debug := ct_carry.CopyNew()

	elapsed_time = time.Since(start)
	fmt.Println(elapsed_time)

	masking = make([]complex128, MAX_SLOT)
	for i := 0; i < BATCH; i++ {
		for j := 0; j < slice_length; j++ {
			if j == slice_length/2 {
				masking[i*slice_length+j] = complex(1.0, 0.0)
			} else {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			}
		}
	}

	twisted_masking = TwistedVec(masking, slice_length, MAX_SLOT)
	cc.eval.MulRelin(ct_carry, twisted_masking, ct_carry)
	cc.eval.Rescale(ct_carry, ct_carry)
	cc.eval.Rotate(ct_carry, MAX_SLOT/2, ct_carry)

	var one *rlwe.Ciphertext
	var err error
	one = ct_carry.CopyNew()
	for i := 0; i < log_slice_half; i++ {
		var temp *rlwe.Ciphertext
		temp, err = cc.eval.RotateNew(one, -1*(1<<i)*BATCH)
		if err != nil {
			fmt.Println(err)
		}

		cc.eval.Add(one, temp, one)
	}

	//PrintDebug(slice_length, *cc.params, one, masking, cc.decryptor, cc.encoder)

	var rtn *rlwe.Ciphertext
	cc.eval.MulRelin(ciphertext, one, ciphertext)
	cc.eval.Rescale(ciphertext, ciphertext)

	cc.eval.Sub(one, 1, one)
	cc.eval.Mul(one, -1, one)

	var temp *rlwe.Ciphertext
	temp = ciphertext1.CopyNew()
	cc.eval.MulRelin(temp, one, temp)
	cc.eval.Rescale(temp, temp)

	rtn, _ = cc.eval.AddNew(ciphertext, temp)

	if tag == "curve25519" {
		/*masking = make([]complex128, MAX_SLOT)
		for i := 0; i < BATCH; i++ {
			for j := 0; j < slice_length; j++ {
				if j == slice_length/2-1 {
					masking[i*slice_length+j] = complex(1.0, 0.0)
				} else {
					masking[i*slice_length+j] = complex(0.0, 0.0)
				}
			}
		}

		twisted_masking = TwistedVec(masking, slice_length, MAX_SLOT)
		cc.eval.MulRelin(debug, twisted_masking, debug)
		cc.eval.Rescale(debug, debug)
		cc.eval.Add(rtn, debug, rtn)*/
	}

	return rtn
}
