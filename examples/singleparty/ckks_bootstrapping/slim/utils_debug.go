package main

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func ComputePrec(base int, slice_length int, params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128, min_err float64, avg_err float64) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	valuesTest = InvTwistedVec(valuesTest, slice_length, slots)

	//fmt.Println()
	//fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	min_L2_err := 100.0
	avg_L2_err := 0.0
	//max_L2_err := 0.0
	is_carry := true

	for i := 0; i < params.MaxSlots(); i++ {

		// compute real err
		real_round := math.Round(real(valuesTest[i]))
		real := real(valuesTest[i])
		real_err := math.Abs(real_round - real)

		// compute imag err
		imag_round := math.Round(imag(valuesTest[i]))
		imag := imag(valuesTest[i])
		imag_err := math.Abs(imag_round - imag)

		L2_err := math.Sqrt(real_err*real_err + imag_err*imag_err)
		log_L2_err := -1 * math.Log2(L2_err)

		avg_L2_err += log_L2_err
		min_L2_err = min(min_L2_err, log_L2_err)
	}

	for i := 0; i < params.MaxSlots()/slice_length; i++ {
		for j := 0; j < slice_length/2; j++ {
			real_round := int(math.Round(real(valuesTest[i*slice_length+j])))
			if real_round < 0 || real_round >= base {
				is_carry = false
				//fmt.Println(j, real_round)
			}
		}
	}

	//fmt.Println("Max Precision(real) :", max_real_err)
	//fmt.Println("Min Precision(real) :", -1*min_L2_err)
	//fmt.Println("Avg Noise(real) :", -2*avg_L2_err/float64(len(valuesTest)))

	//fmt.Println(valuesTest[:8])
	//fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	//fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	//fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	//precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)
	fmt.Println("is carry? : ", is_carry)
	//fmt.Println(precStats.String())
	fmt.Println()

	return valuesTest, -1 * min_L2_err, -2 * avg_L2_err / float64(len(valuesTest))
}

func ComputePrec_vertical(slice_length int, params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128, min_err float64, avg_err float64) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	//valuesTest = InvTwistedVec(valuesTest, slice_length, slots)

	//fmt.Println()
	//fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	min_L2_err := 100.0
	avg_L2_err := 0.0
	//max_L2_err := 0.0

	for i := 0; i < params.MaxSlots(); i++ {

		// compute real err
		real_round := math.Round(2 * real(valuesTest[i]))
		real := real(2 * valuesTest[i])
		real_err := math.Abs(real_round - real)

		// compute imag err
		imag_round := math.Round(imag(valuesTest[i]))
		imag := imag(valuesTest[i])
		imag_err := math.Abs(imag_round - imag)

		L2_err := math.Sqrt(real_err*real_err + imag_err*imag_err)
		log_L2_err := -1 * math.Log2(L2_err)

		avg_L2_err += log_L2_err
		min_L2_err = min(min_L2_err, log_L2_err)

		if i < 2*slice_length {
			//fmt.Println(valuesTest[i])
		}
	}

	//fmt.Println("Max Precision(real) :", max_real_err)
	//fmt.Println("Min Precision(real) :", -1*min_L2_err)
	//fmt.Println("Avg Noise(real) :", -2*avg_L2_err/float64(len(valuesTest)))

	//fmt.Println(valuesTest[:8])
	//fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	//fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	//fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	//precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	//fmt.Println(precStats.String())
	fmt.Println()

	return valuesTest, -1 * min_L2_err, -2 * avg_L2_err / float64(len(valuesTest))
}

func ComputePrecSymbol(slice_length int, params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	min_real_err := 100.0
	avg_real_err := 0.0
	max_real_err := 0.0

	min_imag_err := 100.0
	avg_imag_err := 0.0
	max_imag_err := 0.0

	for i := 0; i < len(valuesTest)/slice_length; i++ {
		for j := 0; j < slice_length/2; j++ {
			//out_round := math.Round(real(valuesTest[i]))
			real_err := -1 * math.Log2(math.Abs(math.Round(2*real(valuesTest[i*slice_length+j]))-2*real(valuesTest[i*slice_length+j])))
			avg_real_err += real_err
			min_real_err = min(min_real_err, real_err)
			max_real_err = max(max_real_err, real_err)

			imag_err := -1 * math.Log2(math.Abs(math.Round(imag(valuesTest[i*slice_length+j]))-imag(valuesTest[i*slice_length+j])))
			avg_imag_err += imag_err
			min_imag_err = min(min_imag_err, imag_err)
			max_imag_err = max(max_imag_err, imag_err)
		}

	}
	fmt.Println()

	//fmt.Println("Max Precision(real) :", max_real_err)
	fmt.Println("Min Precision(real) :", min_real_err)
	//fmt.Println("Avg Precision(real) :", 2*avg_real_err/float64(len(valuesTest)))

	//fmt.Println("Max Precision(imag) :", max_imag_err)
	fmt.Println("Min Precision(imag) :", min_imag_err)
	//fmt.Println("Avg Precision(imag) :", 2*avg_imag_err/float64(len(valuesTest)))

	return
}

func PrintDebug(slice_length int, params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	valuesTest = InvTwistedVec(valuesTest, slice_length, slots)

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("LogScale: %d \n", ciphertext.LogScale())

	avg_L2_err := 0.0
	min_L2_err := 1000.0

	//min_idx := 0
	is_carry := true

	for i := 0; i < params.MaxSlots(); i++ {

		// compute real err
		real_round := math.Round(2 * real(valuesTest[i]))
		realr := real(2 * valuesTest[i])
		real_err := math.Abs(real_round - realr)

		// compute imag err
		imag_round := math.Round(imag(valuesTest[i]))
		imag := imag(valuesTest[i])
		imag_err := math.Abs(imag_round - imag)

		L2_err := math.Sqrt(real_err*real_err + imag_err*imag_err)
		log_L2_err := -1 * math.Log2(L2_err)

		avg_L2_err += log_L2_err
		min_L2_err = min(min_L2_err, log_L2_err)
		if log_L2_err == min_L2_err {
			//min_idx = i
		}
	}

	for i := 0; i < params.MaxSlots()/slice_length; i++ {
		for j := 0; j < slice_length; j++ {
			real_round := int(math.Round(real(valuesTest[i*slice_length+j])))
			if real_round < 0 || real_round >= 16 {
				is_carry = false
				//fmt.Println(j, real_round)
			}
		}
	}

	//fmt.Println(min_idx, valuesTest[min_idx])
	fmt.Println("is carry? : ", is_carry)
	//fmt.Println("B_in :", math.Log2(B_in/float64(len(valuesTest))))
	fmt.Println("Avg Precision :", avg_L2_err/float64(len(valuesTest)))
	fmt.Println("Max Precision :", min_L2_err)

	//fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
