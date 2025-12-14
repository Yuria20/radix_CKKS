package main

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func PolyMult(ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext, Diag_DFT map[int]*rlwe.Plaintext, Diag_InvDFT map[int]*rlwe.Plaintext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {
	slice_length := Bbox.slice_length
	params := cc.params
	eval := cc.eval
	batch := params.MaxSlots() / slice_length
	var ctxt, DFT_ctxt, DFT_ctxt1, DFT_ctxt2 *rlwe.Ciphertext

	// 1. DFT transform
	fmt.Print("DFT transform... ")
	start := time.Now()
	DFT_ctxt1 = BSGS_FFT(batch, slice_length, ctxt1, Diag_DFT, cc)
	elapsed_time := time.Since(start)
	fmt.Println(elapsed_time)

	fmt.Print("DFT transform... ")
	start = time.Now()
	DFT_ctxt2 = BSGS_FFT(batch, slice_length, ctxt2, Diag_DFT, cc)
	elapsed_time = time.Since(start)
	fmt.Println(elapsed_time)

	// 2. CKKS Mult
	DFT_ctxt, _ = eval.MulRelinNew(DFT_ctxt1, DFT_ctxt2)
	eval.Rescale(DFT_ctxt, DFT_ctxt)

	// 3. Inverse DFT transform
	fmt.Print("DFT transform... ")
	start = time.Now()
	ctxt = BSGS_FFT(batch, slice_length, DFT_ctxt, Diag_InvDFT, cc)
	elapsed_time = time.Since(start)
	fmt.Println(elapsed_time)

	return ctxt
}

func PolyCMult(ctxt1 *rlwe.Ciphertext, DFT_ptxt *rlwe.Plaintext, Diag_DFT map[int]*rlwe.Plaintext, Diag_InvDFT map[int]*rlwe.Plaintext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {
	slice_length := Bbox.slice_length
	params := cc.params
	eval := cc.eval
	batch := params.MaxSlots() / slice_length
	var ctxt, DFT_ctxt, DFT_ctxt1 *rlwe.Ciphertext

	// 1. DFT transform
	fmt.Print("DFT transform... ")
	start := time.Now()
	DFT_ctxt1 = BSGS_FFT(batch, slice_length, ctxt1, Diag_DFT, cc)
	elapsed_time := time.Since(start)
	fmt.Println(elapsed_time)
	//DFT_ctxt2 = BSGS_FFT(batch, slice_length, ctxt2, Diag_DFT, cc)

	// 2. CKKS Mult
	DFT_ctxt, _ = eval.MulRelinNew(DFT_ctxt1, DFT_ptxt)
	eval.Rescale(DFT_ctxt, DFT_ctxt)

	// 3. Inverse DFT transform
	fmt.Print("DFT transform... ")
	start = time.Now()
	ctxt = BSGS_FFT(batch, slice_length, DFT_ctxt, Diag_InvDFT, cc)
	elapsed_time = time.Since(start)
	fmt.Println(elapsed_time)

	return ctxt
}

func Compute_Shift_Index(params ckks.Parameters, slice_length int, batch int) []int {

	shift := []int{}
	// compute ceil(root(slice))
	root_slice := int(math.Ceil(math.Sqrt(float64(slice_length))))

	shiftMap := make(map[int]struct{})
	for i := 0; i < root_slice; i++ {
		index := i * root_slice
		shiftMap[batch*index] = struct{}{}

		shift = append(shift, int(batch*index))
		//fmt.Println(index * batch)
	}

	shiftIndex := make([]int, 0, len(shiftMap))
	for key := range shiftMap {
		shiftIndex = append(shiftIndex, key)
	}

	for j := 0; j < root_slice; j++ {

		for i := 0; i < root_slice; i++ {
			if i*root_slice+j >= slice_length {
				continue
			}
		}

		// j rotation for babystep
		shift = append(shift, int(batch*j))
		//fmt.Println(batch * j)

	}

	log_slice_half := int(math.Log2(float64(slice_length) / 2))

	for i := 0; i < log_slice_half; i++ {

		shift = append(shift, int(-1*(1<<i)*(params.MaxSlots())/slice_length))
		//shift = append(shift, int((1<<i)*(params.MaxSlots())/slice_length))
	}

	return shift
}

func Compute_Shift_Index_special(params ckks.Parameters, slice_length int, batch int) []int {

	shift := []int{}
	// compute ceil(root(slice))
	root_slice := int(math.Ceil(math.Sqrt(float64(slice_length))))

	shiftMap := make(map[int]struct{})
	for i := 0; i < root_slice; i++ {
		index := i * root_slice
		shiftMap[batch*index] = struct{}{}

		shift = append(shift, int(batch*index))
		//fmt.Println(index * batch)
	}

	shiftIndex := make([]int, 0, len(shiftMap))
	for key := range shiftMap {
		shiftIndex = append(shiftIndex, key)
	}

	for j := 0; j < root_slice; j++ {

		for i := 0; i < root_slice; i++ {
			if i*root_slice+j >= slice_length {
				continue
			}
		}

		// j rotation for babystep
		shift = append(shift, int(batch*j))
		//fmt.Println(batch * j)

	}

	log_slice_half := int(math.Log2(float64(2*slice_length) / 2))

	for i := 0; i < log_slice_half; i++ {

		shift = append(shift, int(-(1<<i)*(params.MaxSlots())/(2*slice_length)))
	}

	return shift
}

func CalculateL2NormStats(a, b []complex128) (avgLog2, maxLog2 float64, err error) {
	if len(a) != len(b) {
		return 0, 0, errors.New("input arrays must have the same length")
	}

	if len(a) == 0 {
		return 0, 0, nil
	}

	var sumLog2 float64
	maxLog2 = math.Inf(-1)

	for i := 0; i < len(a); i++ {
		diff := a[i] - b[i]
		norm := cmplx.Abs(diff)

		log2Norm := math.Log2(norm)

		sumLog2 += log2Norm
		if log2Norm > maxLog2 {
			maxLog2 = log2Norm
		}
	}

	avgLog2 = sumLog2 / float64(len(a))

	return avgLog2, maxLog2, nil
}

func BSGS_plain_Gen(params ckks.Parameters, slice_length int, batch int, lv int, Diag_mat [][]complex128, cc Context) map[int]*rlwe.Plaintext {

	pre_pt := make(map[int]*rlwe.Plaintext, params.MaxSlots())

	// compute ceil(root(slice))
	root_slice := int(math.Ceil(math.Sqrt(float64(slice_length))))
	valuesDebug := make([]complex128, params.MaxSlots())

	for j := 0; j < root_slice; j++ {

		for i := 0; i < root_slice; i++ {

			if i*root_slice+j >= slice_length {
				continue
			}
			rot_vec := Rotation(Diag_mat[batch*(i*root_slice+j)], -batch*j)
			//Giantstep[i*batch*root_slice].Scale = rlwe.NewScale(cc.params.Q()[Giantstep[i*batch*root_slice].Level()])

			pt := ckks.NewPlaintext(params, lv)
			pt.Scale = rlwe.NewScale(cc.params.Q()[lv])
			cc.encoder.Encode(rot_vec, pt)

			if i == 0 && j == 0 {
				cc.encoder.Decode(pt, valuesDebug)

			}

			pre_pt[batch*(i*root_slice+j)] = pt

			//baby_step, err := cc.eval.MulRelinNew(Giantstep[i*batch*root_slice], pt)

		}
	}

	return pre_pt
}

func BSGS_FFT(batch int, slice_length int, ctxt *rlwe.Ciphertext, Diag_mat map[int]*rlwe.Plaintext, cc Context) *rlwe.Ciphertext {
	var rtn *rlwe.Ciphertext
	//params := cc.params

	// compute ceil(root(slice))
	root_slice := int(math.Ceil(math.Sqrt(float64(slice_length))))

	shiftMap := make(map[int]struct{})
	for i := 0; i < root_slice; i++ {
		index := i * root_slice
		shiftMap[batch*index] = struct{}{}

	}

	shiftIndex := make([]int, 0, len(shiftMap))
	for key := range shiftMap {
		shiftIndex = append(shiftIndex, key)
	}

	Giantstep, err := cc.eval.RotateHoistedNew(ctxt, shiftIndex)
	if err != nil {
		fmt.Println(err)
	}

	for j := 0; j < root_slice; j++ {
		var temp *rlwe.Ciphertext
		start := true

		for i := 0; i < root_slice; i++ {

			if i*root_slice+j >= slice_length {
				continue
			}
			//rot_vec := Rotation(Diag_mat[batch*(i*root_slice+j)], -batch*j)
			//Giantstep[i*batch*root_slice].Scale = rlwe.NewScale(cc.params.Q()[Giantstep[i*batch*root_slice].Level()])

			//pt := ckks.NewPlaintext(*params, params.MaxLevel())
			//pt.Scale = rlwe.NewScale(cc.params.Q()[Giantstep[i*batch*root_slice].Level()])
			//cc.encoder.Encode(rot_vec, pt)
			baby_step, err := cc.eval.MulRelinNew(Giantstep[i*batch*root_slice], Diag_mat[batch*(i*root_slice+j)])

			if err != nil {
				fmt.Println(err)
			}

			if start {
				temp = baby_step.CopyNew()
				start = false
			} else {
				err = cc.eval.Add(temp, baby_step, temp)
				if err != nil {
					fmt.Println(err)
				}
			}
			baby_step = nil
		}

		// j rotation for babystep
		err = cc.eval.Rotate(temp, batch*j, temp)
		//fmt.Println(batch * j)
		if err != nil {
			fmt.Println(err)
		}

		// add rtn ciphertext
		if j == 0 {
			rtn = temp.CopyNew()
		} else {
			err = cc.eval.Add(rtn, temp, rtn)
			if err != nil {
				fmt.Println(err)
			}
		}
		temp = nil
	}

	if err = cc.eval.Rescale(rtn, rtn); err != nil {
		panic(err)
	}

	Giantstep = nil
	runtime.GC()

	return rtn
}

func BSGS_baseREDC(batch int, slice_length int, ctxt *rlwe.Ciphertext, Diag_mat map[int]*rlwe.Plaintext, cc Context) *rlwe.Ciphertext {
	var rtn *rlwe.Ciphertext
	//params := cc.params

	// compute ceil(root(slice))
	root_slice := int(math.Ceil(math.Sqrt(float64(slice_length))))

	shiftMap := make(map[int]struct{})
	for i := 0; i < root_slice; i++ {
		index := i * root_slice
		shiftMap[batch*index] = struct{}{}

	}

	shiftIndex := make([]int, 0, len(shiftMap))
	for key := range shiftMap {
		shiftIndex = append(shiftIndex, key)
	}

	Giantstep, err := cc.eval.RotateHoistedNew(ctxt, shiftIndex)
	if err != nil {
		fmt.Println(err)
	}

	for j := 0; j < root_slice; j++ {
		var temp *rlwe.Ciphertext
		start := true

		for i := 0; i < root_slice; i++ {

			if i*root_slice+j >= slice_length {
				continue
			}
			//rot_vec := Rotation(Diag_mat[batch*(i*root_slice+j)], -batch*j)
			//Giantstep[i*batch*root_slice].Scale = rlwe.NewScale(cc.params.Q()[Giantstep[i*batch*root_slice].Level()])

			//pt := ckks.NewPlaintext(*params, params.MaxLevel())
			//pt.Scale = rlwe.NewScale(cc.params.Q()[Giantstep[i*batch*root_slice].Level()])
			//cc.encoder.Encode(rot_vec, pt)
			baby_step, err := cc.eval.MulRelinNew(Giantstep[i*batch*root_slice], Diag_mat[batch*(i*root_slice+j)])

			if err != nil {
				fmt.Println(err)
			}

			if start {
				temp = baby_step.CopyNew()
				start = false
			} else {
				err = cc.eval.Add(temp, baby_step, temp)
				if err != nil {
					fmt.Println(err)
				}
			}
			//rot_vec = nil
			baby_step = nil
		}

		// j rotation for babystep
		err = cc.eval.Rotate(temp, batch*j, temp)
		//fmt.Println(batch * j)
		if err != nil {
			fmt.Println(err)
		}

		// add rtn ciphertext
		if j == 0 {
			rtn = temp.CopyNew()
		} else {
			err = cc.eval.Add(rtn, temp, rtn)
			if err != nil {
				fmt.Println(err)
			}
		}
		temp = nil
	}

	if err = cc.eval.Rescale(rtn, rtn); err != nil {
		panic(err)
	}

	Giantstep = nil
	runtime.GC()

	return rtn
}

func BSGS_FFT_Optim(
	batch int,
	sliceLength int,
	ctxt *rlwe.Ciphertext,
	DiagMat map[int]*rlwe.Plaintext,
	cc Context,
) *rlwe.Ciphertext {

	root := int(math.Sqrt(float64(sliceLength)))
	if root*root < sliceLength {
		root++
	}

	gsCount := (sliceLength + root - 1) / root

	strideRot := batch * root
	shiftIndex := make([]int, gsCount)
	for i := 0; i < gsCount; i++ {
		shiftIndex[i] = i * strideRot // = batch*(i*root)
	}

	//    GiantMap: map[rotation]*Ciphertext
	GiantMap, err := cc.eval.RotateHoistedNew(ctxt, shiftIndex)
	if err != nil {
		panic(err)
	}

	Giant := make([]*rlwe.Ciphertext, gsCount)
	for i := 0; i < gsCount; i++ {
		Giant[i] = GiantMap[shiftIndex[i]]
	}
	GiantMap = nil

	var rtn *rlwe.Ciphertext

	for j := 0; j < root; j++ {
		if j >= sliceLength {
			break
		}

		var temp *rlwe.Ciphertext

		baseDiag := batch * j
		stepDiag := batch * root

		for i := 0; i < gsCount; i++ {
			idx := i*root + j
			if idx >= sliceLength {
				break
			}

			// Giant[i] * DiagMat[batch*(i*root + j)]
			diagIdx := baseDiag + i*stepDiag
			prod, err := cc.eval.MulRelinNew(Giant[i], DiagMat[diagIdx])
			if err != nil {
				panic(err)
			}

			if temp == nil {
				temp = prod
			} else {
				if err = cc.eval.Add(temp, prod, temp); err != nil {
					panic(err)
				}
			}
		}

		if temp == nil {
			continue
		}

		if err = cc.eval.Rotate(temp, batch*j, temp); err != nil {
			panic(err)
		}

		if rtn == nil {
			rtn = temp
		} else {
			if err = cc.eval.Add(rtn, temp, rtn); err != nil {
				panic(err)
			}
		}
	}

	if err = cc.eval.Rescale(rtn, rtn); err != nil {
		panic(err)
	}

	Giant = nil

	return rtn
}

func Base_Reduction() {

}

func Make_BaseReduction_mat() {

}
