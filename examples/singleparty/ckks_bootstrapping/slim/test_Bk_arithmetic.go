package main

import (
	"fmt"
	"math"
	"math/big"
	"runtime"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
	"gonum.org/v1/gonum/dsp/fourier"
)

// Made by ChatGPT
func warmupCPU(duration time.Duration) {
	numCPU := runtime.NumCPU()
	fmt.Printf("Warming up %d CPU cores for %v...\n", numCPU, duration)

	var wg sync.WaitGroup
	wg.Add(numCPU)

	for i := 0; i < numCPU; i++ {
		go func(id int) {
			defer wg.Done()
			end := time.Now().Add(duration)
			x := 0.0
			for time.Now().Before(end) {
				// floating point
				x += math.Sin(float64(id) + x)
				if x > 1e6 {
					x = 0 // overflow
				}
			}
			//fmt.Printf("Core %d done\n", id)
		}(i)
	}

	wg.Wait()
	fmt.Println("Warm-up complete.")
}

func Test_Bk_arthmetic(Bbox Bk_Arithmetic_toolbox) {

	//==============================================
	//=== 0) Set Large CKKS paramter  ==============
	//==============================================

	// We use Single Thread !!
	runtime.GOMAXPROCS(1)

	fmt.Println()
	fmt.Println("=== Scheme Parameter ===")

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 13

	bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	//log_slice_half := Bbox.log_slice_half

	var LogDefaultScale int

	var Lv int
	var LogQ []int
	LogDefaultScale = 48

	// SAFE-Guard for overflow
	q0 := []int{52, 48}
	qiSlotsToCoeffs := []int{48, 48, 48}
	qiCircuits := []int{48, 48, 55, 48, 55}

	var qiLookUpTable []int
	if Bbox.base_bit == 4 {
		qiLookUpTable = []int{48, 48, 48, 48, 48, 48}
	} else {
		qiLookUpTable = []int{48, 48, 48, 48, 48, 48, 48, 48, 48, 48}
	}

	qiEvalMod := []int{52, 52, 52, 52, 52, 52, 52, 52}
	qiCoeffsToSlots := []int{52, 52, 52}

	// Additonal modulus for bit length
	if 128 <= Bbox.bit_length {
		qiCircuits = append([]int{48}, qiCircuits...)
	}
	if 2048 <= Bbox.bit_length {
		qiCircuits = append([]int{48}, qiCircuits...)
	}

	LogQ = append(q0, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiCircuits...)

	// start level
	Lv = len(LogQ)
	LogQ = append(LogQ, qiLookUpTable...)
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	var err error
	var params ckks.Parameters

	params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,                      // Log2 of the ring degree
		LogQ:            LogQ,                      // Log2 of the ciphertext modulus
		LogP:            []int{52, 52, 52, 52, 52}, // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: LogDefaultScale,           // Log2 of the scale
		Xs:              ring.Ternary{H: 32768},
	})

	if err != nil {
		panic(err)
	}

	//====================================
	//=== BOOTSTRAPPING PARAMETERS ===
	//====================================

	// CoeffsToSlots parameters (homomorphic encoding)
	CtS_scale := *big.NewFloat(float64(0.5)) // We use real bootstrapping
	CoeffsToSlotsParameters := dft.MatrixLiteral{
		Type:         dft.HomomorphicEncode,
		Format:       dft.RepackImagAsReal, // Returns the real and imaginary part into separate ciphertexts
		LogSlots:     params.LogMaxSlots(),
		LevelQ:       params.MaxLevelQ(),
		LevelP:       params.MaxLevelP(),
		LogBSGSRatio: 1,
		Scaling:      &CtS_scale,
		Levels:       []int{1, 1, 1}, //qiCoeffsToSlots
	}

	// Parameters of the homomorphic modular reduction x mod 1
	var Mod1ParametersLiteral mod1.ParametersLiteral
	Mod1ParametersLiteral = mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - CoeffsToSlotsParameters.Depth(true),
		LogScale:        52,               // Matches qiEvalMod
		Mod1Type:        mod1.CosDiscrete, // Multi-interval Chebyshev interpolation
		Mod1Degree:      31,               // Depth iter_num
		DoubleAngle:     3,                // Depth 3
		K:               16,               // With EphemeralSecretWeight = 32 and 2^{155 slots, ensures < 2^{-138.7} failure probability
		LogMessageRatio: 4,                // q/|m| = 2^10
		Mod1InvDegree:   0,                // Depth 0
	}

	// SlotsToCoeffs parameters (homomorphic decoding)
	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type:         dft.HomomorphicDecode,
		LogSlots:     params.LogMaxSlots(),
		LogBSGSRatio: 1,
		LevelP:       params.MaxLevelP(),
		Levels:       []int{1, 1, 1}, // qiSlotsToCoeffs
	}

	SlotsToCoeffsParameters.LevelQ = 4

	// Custom bootstrapping.Parameters.
	// All fields are public and can be manually instantiated.
	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		Mod1ParametersLiteral:   Mod1ParametersLiteral,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		EphemeralSecretWeight:   32, // > 128bit secure for LogN=16 and LogQP = 1155
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}

	// We pring some information about the bootstrapping parameters (which are identical to the residual parameters in this example).
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1555, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	//fmt.Print("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}

	var eval *ckks.Evaluator
	eval = ckks.NewEvaluator(params, evk)
	rot := -1 * params.MaxSlots() / slice_length

	//fmt.Println(Compute_Shift_Index(params, slice_length, params.MaxSlots()/slice_length))
	galEls := []uint64{
		// The galois element for the cyclic rotations by 55positions to the left.
		params.GaloisElement(rot),
		// The galois element for the complex conjugatation.
		params.GaloisElementForComplexConjugation(),
	}

	//params.GaloisElementOrderTwoOrthogonalSubgroup()

	shift_idx := Compute_Shift_Index(params, slice_length, params.MaxSlots()/slice_length)

	for i := 0; i < len(shift_idx); i++ {
		galEls = append(galEls, []uint64{params.GaloisElement(shift_idx[i])}...)
	}

	//fmt.Println("Required Galois Key")
	//fmt.Println(galEls)

	// We then generate the `rlwe.GaloisKey`s element that corresponds to these galois elements.
	// And we update the evaluator's `rlwe.EvaluationKeySet` with the new keys.
	eval = eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...))

	// Instantiates the bootstrapper
	var btp *bootstrapping.Evaluator

	if btp, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	var cc Context
	cc.params = &params
	cc.encoder = encoder
	cc.encryptor = encryptor
	cc.decryptor = decryptor
	cc.eval = eval
	cc.btp = btp

	//==============================================
	//=== 1) Encrypt Large Integer  ================
	//==============================================

	batch := int(float64(params.MaxSlots()) / float64(slice_length))
	fmt.Printf("Large Integer parameter : bit_length=%d, batch=%d, lazycarry_iter=%d, carry_iter=%d", bit_length, batch, Bbox.lazy_iter[bit_length], Bbox.carry_iter[bit_length])
	fmt.Println()
	fmt.Println()

	fmt.Printf("make DFT matrix...")
	if file, _ := FileExists("DFT_" + strconv.Itoa(slice_length)); file == false {
		// Generate FFT matrix
		Normalized_DFT := GenerateSpecialNormalizedDFT(slice_length)
		//Normalized_DFT = MatrixPadding(Normalized_DFT, slice_length, params.MaxSlots())
		//Normalized_DFT2 = MatrixPadding(Normalized_DFT2, slice_length, params.MaxSlots())
		Twisted_DFT := TwistedMatrix(Normalized_DFT, slice_length, params.MaxSlots())
		//Twisted_DFT2 := TwistedMatrix(Normalized_DFT2, slice_length, params.MaxSlots())
		Diag_DFT := Row_To_Diagonal(Twisted_DFT)
		//Diag_DFT2 := Row_To_Diagonal(Twisted_DFT2)
		//fmt.Println(Diag_DFT2[0][0])
		Normalized_DFT = nil
		Twisted_DFT = nil
		runtime.GC()

		Plain_DFT := BSGS_plain_Gen(params, slice_length, batch, Lv-1, Diag_DFT, cc)

		if err := SavePlaintextMap("DFT_"+strconv.Itoa(slice_length), Plain_DFT); err != nil {
			panic(err)
		}

		if err != nil {
			panic(err)
		}
		Diag_DFT = nil
		runtime.GC()
	}

	Plain_DFT, err := LoadPlaintextMap("DFT_"+strconv.Itoa(slice_length), params)

	if file, _ := FileExists("InvDFT_" + strconv.Itoa(slice_length)); file == false {
		InvDFT := GenerateNormalizedInvDFT_with_masking(slice_length)
		//InvDFT = MatrixPadding(InvDFT, slice_length, params.MaxSlots())
		Twisted_InvDFT := TwistedMatrix(InvDFT, slice_length, params.MaxSlots())
		Diag_InvDFT := Row_To_Diagonal(Twisted_InvDFT)

		InvDFT = nil
		Twisted_InvDFT = nil
		runtime.GC()

		Plain_InvDFT := BSGS_plain_Gen(params, slice_length, batch, Lv-3, Diag_InvDFT, cc)

		if err := SavePlaintextMap("InvDFT_"+strconv.Itoa(slice_length), Plain_InvDFT); err != nil {
			panic(err)
		}

		if err != nil {
			panic(err)
		}
		Diag_InvDFT = nil
		runtime.GC()
	}
	fmt.Println("Done!")

	Plain_InvDFT, err := LoadPlaintextMap("InvDFT_"+strconv.Itoa(slice_length), params)

	//==============================================
	//=== 2) Mult Large Integer : MultPoly =========
	//==============================================

	iter_num := 1
	lazy_time := make([]float64, iter_num)
	exact_time := make([]float64, iter_num)
	min_err := make([]float64, iter_num)
	avg_err := make([]float64, iter_num)

	warmupCPU(30 * time.Second)
	fmt.Println("=== Multiplication Start ===")
	for i := 0; i < iter_num; i++ {

		// == 1-1 : packing
		value1 := make([]complex128, params.MaxSlots())
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			for j := 0; j < slice_length/2; j++ {
				nBig := sampling.RandInt(big.NewInt(int64(base)))
				nInt := nBig.Int64()
				value1[i*slice_length+j] = complex(float64(nInt), 0)
			}
			for j := 0; j < slice_length/2; j++ {
				value1[i*slice_length+j+slice_length/2] = complex(0, 0)
			}
		}

		value2 := make([]complex128, params.MaxSlots())
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			for j := 0; j < slice_length/2; j++ {
				nBig := sampling.RandInt(big.NewInt(int64(base)))
				nInt := nBig.Int64()
				value2[i*slice_length+j] = complex(float64(nInt), 0)
			}
			for j := 0; j < slice_length/2; j++ {
				value2[i*slice_length+j+slice_length/2] = complex(0, 0)
			}
		}

		// Encrypt value 1

		twisted_values1 := TwistedVec(value1, slice_length, params.MaxSlots())
		twisted_values2 := TwistedVec(value2, slice_length, params.MaxSlots())

		plaintext1 := ckks.NewPlaintext(params, Lv-1)
		if err := encoder.Encode(twisted_values1, plaintext1); err != nil {
			panic(err)
		}

		ciphertext1, err := encryptor.EncryptNew(plaintext1)
		if err != nil {
			panic(err)
		}

		// Encrypt value 2
		plaintext2 := ckks.NewPlaintext(params, Lv-1)
		if err := encoder.Encode(twisted_values2, plaintext2); err != nil {
			panic(err)
		}

		ciphertext2, err := encryptor.EncryptNew(plaintext2)
		if err != nil {
			panic(err)
		}

		// == 1-2 : Compute Want
		fft := fourier.NewCmplxFFT(slice_length)
		valuesWant_cpx := make([]complex128, params.MaxSlots())
		valuesWant_int := make([]*big.Int, params.MaxSlots())
		valuesWant_carry := make([]*big.Int, params.MaxSlots())

		valuesWant := make([]complex128, params.MaxSlots())
		valuesTest := make([]complex128, params.MaxSlots())
		valuesWant_big := make([]*big.Int, batch)
		valuesTest_big := make([]*big.Int, batch)

		// FFT transform
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			start := i * slice_length
			end := start + slice_length
			seg1 := value1[start:end]
			seg2 := value2[start:end]

			fft.Coefficients(value1[start:end], seg1)
			fft.Coefficients(value2[start:end], seg2)
		}

		// Mult over FFT form
		for i := 0; i < params.MaxSlots(); i++ {
			valuesWant_cpx[i] = value1[i] * value2[i]
			valuesWant_cpx[i] /= complex(float64(slice_length), 0)
		}

		// Inverse FFT transform
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			start := i * slice_length
			end := start + slice_length
			seg := valuesWant_cpx[start:end]

			fft.Sequence(valuesWant_cpx[start:end], seg)
		}

		// Convert Compolex to Int
		for i := 0; i < params.MaxSlots(); i++ {
			valuesWant_int[i] = big.NewInt(int64(math.Round(real(valuesWant_cpx[i]))))
		}

		// Carry
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			for j := 0; j < slice_length; j++ {
				if j != slice_length-1 {
					var temp *big.Int
					temp = new(big.Int).Set(valuesWant_int[i*slice_length+j])
					temp.Quo(temp, big.NewInt(int64(base)))
					valuesWant_carry[i*slice_length+j] = valuesWant_int[i*slice_length+j].Mod(valuesWant_int[i*slice_length+j], big.NewInt(int64(base)))
					valuesWant_int[i*slice_length+j+1] = valuesWant_int[i*slice_length+j+1].Add(valuesWant_int[i*slice_length+j+1], temp)
				} else {
					valuesWant_carry[i*slice_length+j] = valuesWant_int[i*slice_length+j].Mod(valuesWant_int[i*slice_length+j], big.NewInt(int64(base)))
				}
			}
		}

		// Carry
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			for j := 0; j < slice_length/2; j++ {
				valuesWant_carry[i*slice_length+j+slice_length/2] = big.NewInt(0)
			}
		}

		// Convert Int to Complex
		for i := 0; i < params.MaxSlots(); i++ {
			valuesWant[i] = complex(float64(valuesWant_carry[i].Uint64()), 0)
		}

		mult_time := time.Now()
		var ciphertext *rlwe.Ciphertext

		// lazy mult
		ciphertext = LazyMult(ciphertext1, ciphertext2, Plain_DFT, Plain_InvDFT, Bbox, cc)
		lazy_elapsed := time.Since(mult_time)
		fmt.Println("lazy time : ", lazy_elapsed)
		//PrintDebug(slice_length, params, ciphertext, valuesWant, decryptor, encoder)

		// excat mult
		ciphertext = LazyCarry2Carry(ciphertext, Bbox, cc)

		total_elapsed := time.Since(mult_time)
		fmt.Println("exact time : ", total_elapsed)
		fmt.Println()

		//==============================================
		//=== 53 Decrypt and Validation ================
		//==============================================

		plaintext := decryptor.DecryptNew(ciphertext)
		encoder.Decode(plaintext, valuesTest)

		valuesTest = InvTwistedVec(valuesTest, slice_length, params.MaxSlots())

		// Time
		fmt.Println("=== Experiment Result ===")

		// Accuracy

		// Want : Convert Radix int to Int
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			base_temp := big.NewInt(1)
			valuesWant_big[i] = big.NewInt(0)
			for j := 0; j < slice_length/2; j++ {
				valuesWant_big[i].Add(valuesWant_big[i], big.NewInt(int64(math.Round(real(valuesWant[i*slice_length+j])))).Mul(big.NewInt(int64(math.Round(real(valuesWant[i*slice_length+j])))), base_temp))
				base_temp.Mul(base_temp, big.NewInt(int64(base)))
			}

			Mod := big.NewInt(int64(base))
			Mod = Mod.Exp(Mod, big.NewInt(int64(slice_length/2)), nil)

			valuesWant_big[i] = valuesWant_big[i].Mod(valuesWant_big[i], Mod)
		}

		// Test : Convert Radix int to Int
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			base_temp := big.NewInt(1)
			valuesTest_big[i] = big.NewInt(0)
			for j := 0; j < slice_length/2; j++ {
				valuesTest_big[i].Add(valuesTest_big[i], big.NewInt(int64(math.Round(real(valuesTest[i*slice_length+j])))).Mul(big.NewInt(int64(math.Round(real(valuesTest[i*slice_length+j])))), base_temp))
				base_temp.Mul(base_temp, big.NewInt(int64(base)))
			}

			Mod := big.NewInt(int64(base))
			Mod = Mod.Exp(Mod, big.NewInt(int64(slice_length/2)), nil)

			valuesTest_big[i] = valuesTest_big[i].Mod(valuesTest_big[i], Mod)
		}

		acc := 0.0
		for i := 0; i < len(valuesTest_big); i++ {
			if valuesTest_big[i].Cmp(valuesWant_big[i]) == 0 {
				acc += 1.0
			}
		}
		fmt.Printf("Large Integer accuracy : %0.4f\n", acc/float64(batch))

		//fmt.Print("Test(Big) : ")
		//fmt.Printf("[%d, %d, %d, ... , %d]\n", valuesTest_big[0], valuesTest_big[1], valuesTest_big[2], valuesTest_big[len(valuesTest_big)-1])
		//fmt.Print("Want(Big) : ")
		//fmt.Printf("[%d, %d, %d, ... , %d]\n", valuesWant_big[0], valuesWant_big[1], valuesWant_big[2], valuesWant_big[len(valuesWant_big)-1])

		// precision

		lazy_time[i] = float64(lazy_elapsed.Seconds())
		exact_time[i] = float64(total_elapsed.Seconds())

		_, min, avg := ComputePrec(base, slice_length, params, ciphertext, valuesWant, decryptor, encoder)

		min_err[i] = min
		avg_err[i] = avg

		fmt.Println()
	}

	var lazy_mean, exact_mean, min_mean, avg_mean float64
	for i := 0; i < iter_num; i++ {
		lazy_mean += lazy_time[i]
		exact_mean += exact_time[i]
		min_mean += min_err[i]
		avg_mean += avg_err[i]
	}
	lazy_mean /= float64(iter_num)
	exact_mean /= float64(iter_num)
	min_mean /= float64(iter_num)
	avg_mean /= float64(iter_num)

	fmt.Println("lazy lat. : ", lazy_mean, " s")
	fmt.Println("lazy amot.. : ", 1000*lazy_mean/float64(batch), "ms")
	fmt.Println("exact lat. : ", exact_mean, " s")
	fmt.Println("exact amot. : ", 1000*exact_mean/float64(batch), " ms")
	fmt.Println("min err : ", min_mean)
	fmt.Println("avg err : ", avg_mean)

}
