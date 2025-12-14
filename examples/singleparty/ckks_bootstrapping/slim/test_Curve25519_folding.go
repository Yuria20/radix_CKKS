package main

import (
	"crypto/rand"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"strconv"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/dsp/fourier"
)

// RandomBigInt returns a uniform random integer with the
// top bit set, so the result has EXACTLY 'bits' bits.
func RandomBigInt_mod(mod *big.Int) (*big.Int, error) {
	// max = 1 << bits
	max := mod

	// 0 ≤ n < max
	n, err := rand.Int(rand.Reader, max) // crypto/rand.Int :contentReference[oaicite:1]{index=1}
	if err != nil {
		return nil, err
	}

	return n, nil
}

func DecomposeInt64Base(x *big.Int, base int64, k int) []int64 {
	if base < 2 {
		panic("base must be ≥ 2")
	}
	if k <= 0 {
		panic("k must be positive")
	}

	B := big.NewInt(base)
	z := new(big.Int).Set(x) // copy; keep caller’s Int intact
	zero := big.NewInt(0)
	digits := make([]int64, 0, k) // pre-allocate

	quot, rem := new(big.Int), new(big.Int)

	for z.Cmp(zero) > 0 {
		quot, rem = quot.QuoRem(z, B, rem) // z/B, z%B  :contentReference[oaicite:0]{index=0}
		digits = append(digits, rem.Int64())
		z.Set(quot)

		if len(digits) > k {
			panic("number requires more than k digits; increase k")
		}
	}

	if len(digits) == 0 {
		digits = append(digits, 0)
	}

	// Zero-pad up to k
	for len(digits) < k {
		digits = append(digits, 0)
	}

	return digits
}

func Test_curve25519_arthmetic(Bbox Bk_Arithmetic_toolbox) {

	//==============================================
	//=== 0) Set Large CKKS paramter  ==============
	//==============================================

	// We use Single Thread !!
	runtime.GOMAXPROCS(1)

	fmt.Println()
	fmt.Println("=== Scheme Parameter ===")

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 16

	bit_length := Bbox.bit_length
	base := Bbox.base
	//base_bit := Bbox.base_bit
	slice_length := Bbox.slice_length
	//log_slice_half := Bbox.log_slice_half

	LogDefaultScale := 48

	q0 := []int{52, 48}                  // 3) ScaleDown & 4) ModUp :
	qiSlotsToCoeffs := []int{48, 48, 48} // 1) SlotsToCoeffs
	qiCircuits := []int{48, 48, 48, 55, 48, 55}
	qiLookUpTable := []int{48, 48, 48, 48, 48, 48}     // 0) Circuit in the slot domain
	qiEvalMod := []int{52, 52, 52, 52, 52, 52, 52, 52} // 6) EvalMod : 8
	qiCoeffsToSlots := []int{52, 52, 52}               // 52 CoeffsToSlots

	if Bbox.bit_length == 2048 {
		qiCircuits = append([]int{48}, qiCircuits...)
	}

	fmt.Println("base : ", q0)
	fmt.Println("StC : ", qiSlotsToCoeffs)
	fmt.Println("Circuit", qiCircuits)
	fmt.Println("LUT", qiLookUpTable)
	fmt.Println("EvalMod : ", qiEvalMod)
	fmt.Println("CtS : ", qiCoeffsToSlots)

	LogQ := append(q0, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiCircuits...)

	// start level
	Lv := len(LogQ)
	LogQ = append(LogQ, qiLookUpTable...)
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
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
	CtS_scale := *big.NewFloat(float64(0.5))
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
	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - CoeffsToSlotsParameters.Depth(true),
		LogScale:        52,               // Matches qiEvalMod
		Mod1Type:        mod1.CosDiscrete, // Multi-interval Chebyshev interpolation
		Mod1Degree:      31,               // Depth 5
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
		params.GaloisElement((slice_length / 4) * (params.MaxSlots()) / slice_length),
		params.GaloisElement(-255 * (params.MaxSlots()) / slice_length),
	}

	params.GaloisElementOrderTwoOrthogonalSubgroup()

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

	//btp.Bootstrap()

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
	big_two := big.NewInt(int64(2))
	mod := big_two.Exp(big_two, big.NewInt(255), nil)
	mod = mod.Sub(mod, big.NewInt(19))
	fmt.Printf("Large Integer parameter : mod=Curve25519, slice_length=%d, bit_length=%d, batch=%d, lazycarry_iter=%d, carry_iter=%d", slice_length, bit_length, batch, Bbox.lazy_iter[bit_length], Bbox.carry_iter[bit_length])
	fmt.Println()
	fmt.Println()

	if file, _ := FileExists("DFT_mod_" + strconv.Itoa(slice_length)); file == false {
		// Generate FFT matrix
		Normalized_DFT := GenerateSpecialNormalizedDFT(slice_length)
		Normalized_DFT = MatrixPadding(Normalized_DFT, slice_length, params.MaxSlots())
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

		if err := SavePlaintextMap("DFT_mod_"+strconv.Itoa(slice_length), Plain_DFT); err != nil {
			panic(err)
		}

		if err != nil {
			panic(err)
		}
		Diag_DFT = nil
		runtime.GC()
	}

	Plain_DFT, err := LoadPlaintextMap("DFT_mod_"+strconv.Itoa(slice_length), params)

	if file, _ := FileExists("InvDFT_mod_" + strconv.Itoa(slice_length)); file == false {
		InvDFT := GenerateNormalizedInvDFT(slice_length)
		InvDFT = MatrixPadding(InvDFT, slice_length, params.MaxSlots())
		Twisted_InvDFT := TwistedMatrix(InvDFT, slice_length, params.MaxSlots())
		Diag_InvDFT := Row_To_Diagonal(Twisted_InvDFT)

		InvDFT = nil
		Twisted_InvDFT = nil
		runtime.GC()

		Plain_InvDFT := BSGS_plain_Gen(params, slice_length, batch, Lv-3, Diag_InvDFT, cc)

		if err := SavePlaintextMap("InvDFT_mod_"+strconv.Itoa(slice_length), Plain_InvDFT); err != nil {
			panic(err)
		}

		if err != nil {
			panic(err)
		}
		Diag_InvDFT = nil
		runtime.GC()
	}

	Plain_InvDFT, err := LoadPlaintextMap("InvDFT_mod_"+strconv.Itoa(slice_length), params)

	//fmt.Println(Lv-1, Lv-2)
	runtime.GC()

	//fmt.Println("hit!")
	lazy_time := make([]float64, 5)
	exact_time := make([]float64, 6)
	min_err := make([]float64, 6)
	avg_err := make([]float64, 6)

	for k := 0; k < 5; k++ {

		values1_big := make([]*big.Int, batch)
		values2_big := make([]*big.Int, batch)
		values1_int := make([]int64, params.MaxSlots())
		values2_int := make([]int64, params.MaxSlots())

		for i := 0; i < batch; i++ {
			values1_big[i], _ = RandomBigInt_mod(mod)
			values2_big[i], _ = RandomBigInt_mod(mod)
			digits1 := DecomposeInt64Base(values1_big[i], int64(base), slice_length)
			digits2 := DecomposeInt64Base(values2_big[i], int64(base), slice_length)
			//fmt.Println(digits1)
			//fmt.Println(digits2)
			for j := 0; j < slice_length; j++ {
				values1_int[i*slice_length+j] = digits1[j]
				values2_int[i*slice_length+j] = digits2[j]
			}
		}

		// == 1-1 : packing
		value1 := make([]complex128, params.MaxSlots())
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			for j := 0; j < slice_length; j++ {
				value1[i*slice_length+j] = complex(float64(values1_int[i*slice_length+j]), 0)
			}
		}

		value2 := make([]complex128, params.MaxSlots())
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			for j := 0; j < slice_length/2; j++ {
				value2[i*slice_length+j] = complex(float64(values2_int[i*slice_length+j]), 0)
			}
		}

		//fmt.Println(value1[:12])
		//fmt.Println(value2[:12])

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

		valueN_int := make([]int64, params.MaxSlots())
		valueN_cpx := make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			valueN_int = DecomposeInt64Base(mod, int64(base), slice_length)
			//value2N_int = DecomposeInt64Base(big.NewInt(0.0).Mul(mod, big.NewInt(2.0)), int64(base), slice_length)

			for j := 0; j < slice_length; j++ {
				valueN_cpx[i*slice_length+j] = complex(float64(valueN_int[j]), 0)
				//value2N_cpx[i*slice_length+j] = complex(float64(value2N_int[j]), 0)
			}
		}
		//==============================================
		//=== 2) Mult Large Integer : MultPoly =========
		//==============================================

		fmt.Println("=== Multiplication Start ===")

		mult_time := time.Now()
		var ciphertext *rlwe.Ciphertext

		ciphertext = ExactMult(ciphertext1, ciphertext2, Plain_DFT, Plain_InvDFT, Bbox, cc)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)
		fmt.Print("Cleaning Start...")
		start := time.Now()
		ciphertext = Cleaning_with_vectorized_evaluation(ciphertext, Bbox, cc)
		elapsed := time.Since(start)
		fmt.Println(elapsed)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		//fmt.Println(ciphertext.Scale.Value.Float64())
		masking := make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			for j := 0; j < slice_length; j++ {
				if slice_length/4 <= j && j < slice_length/2 {
					masking[i*slice_length+j] = complex(38.0, 0.0)
				} else {
					masking[i*slice_length+j] = complex(0.0, 0.0)
				}
			}
		}

		twisted_masking := TwistedVec(masking, slice_length, params.MaxSlots())

		pt := ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
		cc.encoder.Encode(twisted_masking, pt)

		var temp *rlwe.Ciphertext
		temp, _ = cc.eval.MulRelinNew(ciphertext, pt)
		cc.eval.Rescale(temp, temp)
		//fmt.Println(ciphertext.Scale.Value.Float64())
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		err = cc.eval.Rotate(temp, (slice_length/4)*(params.MaxSlots())/slice_length, temp)
		if err != nil {
			fmt.Println(err)
		}

		masking = make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			for j := 0; j < slice_length; j++ {
				if j < slice_length/4 {
					masking[i*slice_length+j] = complex(1.0, 0.0)
				} else {
					masking[i*slice_length+j] = complex(0.0, 0.0)
				}
			}
		}

		twisted_masking = TwistedVec(masking, slice_length, params.MaxSlots())

		pt = ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
		cc.encoder.Encode(twisted_masking, pt)

		ciphertext, _ = cc.eval.MulRelinNew(ciphertext, pt)
		cc.eval.Rescale(ciphertext, ciphertext)
		cc.eval.Add(temp, ciphertext, ciphertext)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		for iter := 0; iter < 2; iter++ {
			fmt.Printf("LazyCarry Start...")
			start = time.Now()
			ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
			elapsed = time.Since(start)
			fmt.Printf(" %s\n", elapsed)
			//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)
		}

		ciphertext = LazyCarry2Carry(ciphertext, Bbox, cc)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		masking = make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			for j := 0; j < slice_length; j++ {
				if slice_length/4 <= j && j < slice_length/2 {
					masking[i*slice_length+j] = complex(38.0, 0.0)
				} else {
					masking[i*slice_length+j] = complex(0.0, 0.0)
				}
			}
		}

		twisted_masking = TwistedVec(masking, slice_length, params.MaxSlots())

		pt = ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
		cc.encoder.Encode(twisted_masking, pt)

		temp, _ = cc.eval.MulRelinNew(ciphertext, pt)
		cc.eval.Rescale(temp, temp)
		//fmt.Println(ciphertext.Scale.Value.Float64())
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		err = cc.eval.Rotate(temp, (slice_length/4)*(params.MaxSlots())/slice_length, temp)
		if err != nil {
			fmt.Println(err)
		}

		masking = make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			for j := 0; j < slice_length; j++ {
				if j < slice_length/4 {
					masking[i*slice_length+j] = complex(1.0, 0.0)
				} else {
					masking[i*slice_length+j] = complex(0.0, 0.0)
				}
			}
		}

		twisted_masking = TwistedVec(masking, slice_length, params.MaxSlots())

		pt = ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
		cc.encoder.Encode(twisted_masking, pt)

		ciphertext, _ = cc.eval.MulRelinNew(ciphertext, pt)
		cc.eval.Rescale(ciphertext, ciphertext)
		cc.eval.Add(temp, ciphertext, ciphertext)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		for iter := 0; iter < 2; iter++ {
			fmt.Printf("LazyCarry Start...")
			start = time.Now()
			ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
			elapsed = time.Since(start)
			fmt.Printf(" %s\n", elapsed)
			//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)
		}

		lazy_elapsed := time.Since(mult_time)

		ciphertext = LazyCarry2Carry(ciphertext, Bbox, cc)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)
		ciphertext = Cleaning_with_vectorized_evaluation(ciphertext, Bbox, cc)
		//PrintDebug(slice_length, params, ciphertext, twisted_values1, cc.decryptor, cc.encoder)

		var cipherN *rlwe.Ciphertext
		valueN_cpx = TwistedVec(valueN_cpx, slice_length, params.MaxSlots())
		plain_N := ckks.NewPlaintext(params, ciphertext.Level())
		//plain_N_rad.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
		encoder.Encode(valueN_cpx, plain_N)
		cipherN, err = cc.encryptor.EncryptNew(plain_N)
		if err != nil {
			fmt.Println(err)
		}

		ciphertext = ConditonalSub(ciphertext, cipherN, Bbox, cc, "curve25519")

		//fmt.Println(ciphertext.Level())

		total_elapsed := time.Since(mult_time)

		plaintext := decryptor.DecryptNew(ciphertext)
		encoder.Decode(plaintext, valuesTest)

		valuesTest = InvTwistedVec(valuesTest, slice_length, params.MaxSlots())

		//fmt.Println("mod : ", mod)
		//fmt.Println("mod : ", big.NewInt(0.0).Mul(mod, big.NewInt(2.0)))

		// Time
		fmt.Println("=== Experiment Result ===")
		fmt.Printf("Latency(Lazy): %0.3f s\n", lazy_elapsed.Seconds())
		fmt.Printf("Amortized(Lazy): %0.3f ms\n", float64(lazy_elapsed.Milliseconds())*float64(slice_length)/(float64(params.MaxSlots())))
		fmt.Printf("Latency: %0.3f s\n", total_elapsed.Seconds())
		fmt.Printf("Amortized : %0.3f ms\n", float64(total_elapsed.Milliseconds())*float64(slice_length)/(float64(params.MaxSlots())))

		// Accuracy

		// Want : Convert Radix int to Int
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			base_temp := big.NewInt(1)
			valuesWant_big[i] = big.NewInt(0)
			for j := 0; j < slice_length/2; j++ {
				valuesWant_big[i].Add(valuesWant_big[i], big.NewInt(int64(math.Round(real(valuesWant[i*slice_length+j])))).Mul(big.NewInt(int64(math.Round(real(valuesWant[i*slice_length+j])))), base_temp))
				base_temp.Mul(base_temp, big.NewInt(int64(base)))
			}

			Mod := mod
			valuesWant_big[i] = valuesWant_big[i].Mod(valuesWant_big[i], Mod)

			if i == 0 {
				//fmt.Println(valuesWant_big[0])
			}

		}

		//fmt.Println(valuesWant_big)
		// Test : Convert Radix int to Int
		max_rat := big.NewFloat(0.0)

		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			base_temp := big.NewInt(1)
			valuesTest_big[i] = big.NewInt(0)
			for j := 0; j < slice_length/2; j++ {
				valuesTest_big[i].Add(valuesTest_big[i], big.NewInt(int64(math.Round(real(valuesTest[i*slice_length+j])))).Mul(big.NewInt(int64(math.Round(real(valuesTest[i*slice_length+j])))), base_temp))
				base_temp.Mul(base_temp, big.NewInt(int64(base)))
			}

			big_Test := IntToFloatExact(valuesTest_big[i])
			mod_float := IntToFloatExact(mod)
			mod_float.Quo(big_Test, mod_float)

			if max_rat.Cmp(mod_float) == -1 {
				max_rat = mod_float
			}
			if i == 0 {
				//fmt.Println(valuesTest_big[0])
			}
			//Mod := mod
			//valuesTest_big[i] = valuesTest_big[i].Mod(valuesTest_big[i], Mod)
		}

		//fmt.Println(max_rat)
		//fmt.Println(valuesTest_big)

		acc := 0.0
		for i := 0; i < len(valuesTest_big); i++ {
			if valuesTest_big[i].Cmp(valuesWant_big[i]) == 0 {
				acc += 1.0
			}
		}
		fmt.Printf("Large Integer accuracy : %0.4f", acc/float64(batch))

		//fmt.Print("Test(Big) : ")
		//fmt.Printf("[%d, %d, %d, ... , %d]\n", valuesTest_big[0], valuesTest_big[1], valuesTest_big[2], valuesTest_big[len(valuesTest_big)-1])
		//fmt.Print("Want(Big) : ")
		//fmt.Printf("[%d, %d, %d, ... , %d]\n", valuesWant_big[0], valuesWant_big[1], valuesWant_big[2], valuesWant_big[len(valuesWant_big)-1])
		lazy_time[k] = float64(lazy_elapsed.Seconds())
		exact_time[k] = float64(total_elapsed.Seconds())
		// precision
		_, min, avg := ComputePrec(base, slice_length, params, ciphertext, valuesWant, decryptor, encoder)
		min_err[k] = min
		avg_err[k] = avg
		//ComputePrec(params, ct, valuesWant, decryptor, encoder)
		fmt.Println()
	}

	var exact_mean, lazy_mean, min_mean, avg_mean float64
	for i := 0; i < 5; i++ {
		lazy_mean += lazy_time[i]
		exact_mean += exact_time[i]
		min_mean += min_err[i]
		avg_mean += avg_err[i]
	}
	lazy_mean /= 5.0
	exact_mean /= 5.0
	min_mean /= 5.0
	avg_mean /= 5.0

	fmt.Println("lazy lat. : ", lazy_mean, " s")
	fmt.Println("lazy amot.. : ", 1000*lazy_mean/float64(batch), "ms")
	fmt.Println("exact lat. : ", exact_mean, " s")
	fmt.Println("exact amot. : ", 1000*exact_mean/float64(batch), " ms")
	fmt.Println("min err : ", min_mean)
	fmt.Println("avg err : ", avg_mean)

}
