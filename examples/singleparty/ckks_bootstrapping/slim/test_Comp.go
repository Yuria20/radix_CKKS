package main

import (
	"fmt"
	"math"
	"math/big"
	"runtime"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func Test_Comp(Bbox Bk_Arithmetic_toolbox) {

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

	var LogDefaultScale int

	var Lv int
	var LogQ []int

	if Bbox.base == 16 {
		LogDefaultScale = 48

		q0 := []int{52, 48}                  // 3) ScaleDown & 4) ModUp :
		qiSlotsToCoeffs := []int{48, 48, 48} // 1) SlotsToCoeffs
		qiCircuits := []int{48, 48, 55, 48, 55}

		qiLookUpTable := []int{48, 48, 48, 48, 48, 48} // 0) Circuit in the slot domain

		qiEvalMod := []int{52, 52, 52, 52, 52, 52, 52, 52} // 6) EvalMod : 8
		qiCoeffsToSlots := []int{52, 52, 52}               // 53 CoeffsToSlots

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

	} else {
		LogDefaultScale = 45

		q0 := []int{53, 45}                  // 3) ScaleDown & 4) ModUp :
		qiSlotsToCoeffs := []int{45, 45, 45} // 1) SlotsToCoeffs
		qiCircuits := []int{45, 45, 53, 45, 53}

		qiLookUpTable := []int{45, 45, 45, 45, 45, 45, 45, 45, 45} // 0) Circuit in the slot domain

		qiEvalMod := []int{53, 53, 53, 53, 53, 53, 53, 53} // 6) EvalMod : 8
		qiCoeffsToSlots := []int{53, 53, 53}               // 53 CoeffsToSlots

		if Bbox.bit_length == 2048 {
			qiCircuits = append([]int{45}, qiCircuits...)
		}

		//fmt.Println("base : ", q0)
		//fmt.Println("StC : ", qiSlotsToCoeffs)
		//fmt.Println("Circuit", qiCircuits)
		//fmt.Println("LUT", qiLookUpTable)
		//fmt.Println("EvalMod : ", qiEvalMod)
		//fmt.Println("CtS : ", qiCoeffsToSlots)

		LogQ = append(q0, qiSlotsToCoeffs...)
		LogQ = append(LogQ, qiCircuits...)

		// start level
		Lv = len(LogQ)
		LogQ = append(LogQ, qiLookUpTable...)
		LogQ = append(LogQ, qiEvalMod...)
		LogQ = append(LogQ, qiCoeffsToSlots...)
	}

	var err error
	var params ckks.Parameters
	if Bbox.base == 16 {
		params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            LogN,                      // Log2 of the ring degree
			LogQ:            LogQ,                      // Log2 of the ciphertext modulus
			LogP:            []int{53, 53, 53, 53, 53}, // Log2 of the key-switch auxiliary prime moduli
			LogDefaultScale: LogDefaultScale,           // Log2 of the scale
			Xs:              ring.Ternary{H: 192},
		})
	} else {
		params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            LogN,                      // Log2 of the ring degree
			LogQ:            LogQ,                      // Log2 of the ciphertext modulus
			LogP:            []int{53, 53, 53, 53, 53}, // Log2 of the key-switch auxiliary prime moduli
			LogDefaultScale: LogDefaultScale,           // Log2 of the scale
			Xs:              ring.Ternary{H: 32768},
		})
	}

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
	var Mod1ParametersLiteral mod1.ParametersLiteral
	if Bbox.base == 16 {
		Mod1ParametersLiteral = mod1.ParametersLiteral{
			LevelQ:          params.MaxLevel() - CoeffsToSlotsParameters.Depth(true),
			LogScale:        52,               // Matches qiEvalMod
			Mod1Type:        mod1.CosDiscrete, // Multi-interval Chebyshev interpolation
			Mod1Degree:      31,               // Depth 5
			DoubleAngle:     3,                // Depth 3
			K:               16,               // With EphemeralSecretWeight = 32 and 2^{155 slots, ensures < 2^{-138.7} failure probability
			LogMessageRatio: 4,                // q/|m| = 2^10
			Mod1InvDegree:   0,                // Depth 0
		}

	} else {
		Mod1ParametersLiteral = mod1.ParametersLiteral{
			LevelQ:          params.MaxLevel() - CoeffsToSlotsParameters.Depth(true),
			LogScale:        53,               // Matches qiEvalMod
			Mod1Type:        mod1.CosDiscrete, // Multi-interval Chebyshev interpolation
			Mod1Degree:      31,               // Depth 5
			DoubleAngle:     3,                // Depth 3
			K:               16,               // With EphemeralSecretWeight = 32 and 2^{155 slots, ensures < 2^{-138.7} failure probability
			LogMessageRatio: 8,                // q/|m| = 2^10
			Mod1InvDegree:   0,                // Depth 0
		}
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
		params.GaloisElement(params.MaxSlots() / 2),
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
	fmt.Printf("Large Integer parameter : bit_length=%d, batch=%d, lazycarry_iter=%d, carry_iter=%d", bit_length, batch, Bbox.lazy_iter[bit_length], Bbox.carry_iter[bit_length])
	fmt.Println()
	fmt.Println()

	//fmt.Println(Lv-1, Lv-2)
	//Diag_DFT = nil
	//Diag_InvDFT = nil
	//runtime.GC()

	mod := big.NewInt(int64(base))
	mod = mod.Exp(mod, big.NewInt(int64(slice_length/2-1)), nil)
	mod = mod.Add(mod, big.NewInt(0.0).Exp(big.NewInt(int64(base)), big.NewInt(int64(slice_length/2-2)), nil))
	primeStr := "20001030500114365578276172632409225656496377638593185792480199378537058893741203912381988890014108377553928167405028407253810959620542493980998363450945858435827230425360559865510737689270866756515652866551976510105545537970228013696617490642319004379799051093313365944568373216127132036998920433135243233125809277778869287991992773691437819178746812300604717045712177799445705659577853791653666000889181040230645767816513604639136008269498511864003611922378089262260835535463621902062288442797058009859720284517772755190922161553252699460863880033665341687790834810005136108173760509308004664781377968493112574251293"
	// p-384
	//primeStr := "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"
	mod = new(big.Int)
	mod.SetString(primeStr, 10)

	// == 1-1 : packing
	values1_big := make([]*big.Int, batch)
	values2_big := make([]*big.Int, batch)
	values1_int := make([]int64, params.MaxSlots())
	values2_int := make([]int64, params.MaxSlots())

	for i := 0; i < batch; i++ {
		b := "20339409470481469145876538288174964040612543182783178761894482178651659055521856271732705234358884807705206971189306301040769101050802363680531997488074738297008992287308782219277760591532235283585148351951621996424642128638902869665385225110477050590350041199575415996541349275116302294455582188283042297515958470768906139287768859394591477063772212926573015227566133099886845530519494981343285430642332770328324265590466374227021623410161616931685560000157684592045842840853370619318484472261837445210284357694302274233657545943095958609506835822446051340525392221645072426073978591781607379922290434699102781889934"
		temp := new(big.Int)
		temp.SetString(b, 10)
		values1_big[i] = temp
		//intstr := ""
		//values1_big[i], _ = RandomBigInt_mod(big.NewInt(0.0).Mul(mod, big.NewInt(2.0)))
		//values1_big[i], _ = RandomBigInt_mod(mod)
		//values2_big[i], _ = RandomBigInt_mod(mod)
		values2_big[i] = mod
		//values1_big[i].Add(values1_big[i], mod)

		//fmt.Println(values1_big[i])
		//fmt.Println(values2_big[i])

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

	// Encrypt value 1
	fmt.Println(Lv)
	twisted_values1 := TwistedVec(value1, slice_length, params.MaxSlots())
	twisted_values2 := TwistedVec(value2, slice_length, params.MaxSlots())

	plaintext1 := ckks.NewPlaintext(params, 11)
	if err := encoder.Encode(twisted_values1, plaintext1); err != nil {
		panic(err)
	}

	ciphertext1, err := encryptor.EncryptNew(plaintext1)
	if err != nil {
		panic(err)
	}

	// Encrypt value 2
	plaintext2 := ckks.NewPlaintext(params, 11)
	if err := encoder.Encode(twisted_values2, plaintext2); err != nil {
		panic(err)
	}

	ciphertext2, err := encryptor.EncryptNew(plaintext2)
	if err != nil {
		panic(err)
	}

	// == 1-2 : Compute Want
	//fft := fourier.NewCmplxFFT(slice_length)
	//valuesWant_cpx := make([]complex128, params.MaxSlots())
	//valuesWant_int := make([]*big.Int, params.MaxSlots())
	//valuesWant_carry := make([]*big.Int, params.MaxSlots())

	valuesWant := make([]complex128, params.MaxSlots())
	valuesTest := make([]complex128, params.MaxSlots())
	//values1_big := make([]*big.Int, batch)
	//values2_big := make([]*big.Int, batch)

	valuesWant_big := make([]*big.Int, batch)
	valuesTest_big := make([]*big.Int, batch)

	for i := 0; i < batch; i++ {
		if values1_big[i].Cmp(values2_big[i]) == -1 {
			valuesWant_big[i] = values1_big[i]
		} else {
			rtn := big.NewInt(1).Sub(values1_big[i], values2_big[i])
			valuesWant_big[i] = rtn
		}
	}

	//==============================================
	//=== 2) Comp Large Integer : MultPoly =========
	//==============================================

	fmt.Println("=== Comparison Start ===")

	//lazy_time := make([]float64, 5)
	exact_time := make([]float64, 5)
	min_err := make([]float64, 5)
	avg_err := make([]float64, 5)

	for i := 0; i < 1; i++ {
		mult_time := time.Now()
		var ciphertext *rlwe.Ciphertext
		ciphertext = ciphertext1.CopyNew()
		//PrintDebug(slice_length, params, ciphertext1, value1, decryptor, encoder)
		//PrintDebug(slice_length, params, ciphertext2, value1, decryptor, encoder)
		//ciphertext = Comparison(ciphertext1, ciphertext2, Bbox, cc)
		ciphertext = ConditonalSub(ciphertext, ciphertext2, Bbox, cc, "std")
		//ciphertext = ConditonalSub_with_MSB(ciphertext1, ciphertext2, Bbox, cc)
		total_elapsed := time.Since(mult_time)
		fmt.Println()
		//PrintDebug(slice_length, params, ciphertext, value1, decryptor, encoder)
		//==============================================
		//=== 53 Decrypt and Validation ================
		//==============================================

		// Decrypt
		//if ciphertext.Level() != 0 {
		//	eval.DropLevel(ciphertext, ciphertext.Level())
		//}
		plaintext := decryptor.DecryptNew(ciphertext)
		encoder.Decode(plaintext, valuesTest)

		valuesTest = InvTwistedVec(valuesTest, slice_length, params.MaxSlots())

		// Time
		fmt.Println("=== Experiment Result ===")
		//fmt.Printf("Latency(Lazy): %0.3f s\n", lazy_elapsed.Seconds())
		//fmt.Printf("Amortized(Lazy): %0.3f ms\n", float64(lazy_elapsed.Milliseconds())*float64(slice_length)/(float64(params.MaxSlots())))
		//fmt.Printf("Latency: %0.3f s\n", total_elapsed.Seconds())
		//fmt.Printf("Amortized : %0.3f ms\n", float64(total_elapsed.Milliseconds())*float64(slice_length)/(float64(params.MaxSlots())))

		// Test : Convert Radix int to Int
		for i := 0; i < params.MaxSlots()/slice_length; i++ {
			base_temp := big.NewInt(1)
			valuesTest_big[i] = big.NewInt(0)
			for j := 0; j < slice_length/2; j++ {
				valuesTest_big[i].Add(valuesTest_big[i], big.NewInt(int64(math.Round(real(valuesTest[i*slice_length+j])))).Mul(big.NewInt(int64(math.Round(real(valuesTest[i*slice_length+j])))), base_temp))
				base_temp.Mul(base_temp, big.NewInt(int64(base)))
			}

			//Mod := big.NewInt(int64(base))
			//Mod = Mod.Exp(Mod, big.NewInt(int64(slice_length/2)), nil)

			//valuesTest_big[i] = valuesTest_big[i].Mod(valuesTest_big[i], Mod)
		}
		fmt.Println(valuesTest_big[0])
		acc := 0.0
		for i := 0; i < len(valuesTest_big); i++ {
			if valuesTest_big[i].Cmp(valuesWant_big[i]) == 0 {
				acc += 1.0
				//fmt.Print("index :", i, valuesTest_big[i], valuesWant_big[i])
				//PrintDebug(slice_length, params, ct_carry, twisted_values1, cc.decryptor, cc.encoder)
			} else {
				//fmt.Print("index :", i, valuesTest[i*slice_length], valuesWant_big[i])
			}
		}
		fmt.Printf("Large Integer accuracy : %0.4f", acc/float64(batch))

		//fmt.Print("Test(Big) : ")
		//fmt.Printf("[%d, %d, %d, ... , %d]\n", valuesTest_big[0], valuesTest_big[1], valuesTest_big[2], valuesTest_big[len(valuesTest_big)-1])
		//fmt.Print("Want(Big) : ")
		//fmt.Printf("[%d, %d, %d, ... , %d]\n", valuesWant_big[0], valuesWant_big[1], valuesWant_big[2], valuesWant_big[len(valuesWant_big)-1])

		//lazy_time[i] = float64(lazy_elapsed.Seconds())
		exact_time[i] = float64(total_elapsed.Seconds())

		_, min, avg := ComputePrec(base, slice_length, params, ciphertext, valuesWant, decryptor, encoder)

		min_err[i] = min
		avg_err[i] = avg
		//ComputePrec(params

	}

	var exact_mean, min_mean, avg_mean float64
	for i := 0; i < 5; i++ {
		//lazy_mean += lazy_time[i]
		exact_mean += exact_time[i]
		min_mean += min_err[i]
		avg_mean += avg_err[i]
	}
	//lazy_mean /= 5.0
	exact_mean /= 5.0
	min_mean /= 5.0
	avg_mean /= 5.0

	//fmt.Println("lazy lat. : ", lazy_mean, " s")
	//fmt.Println("lazy amot.. : ", 1000*lazy_mean/float64(batch), "ms")
	fmt.Println("exact lat. : ", exact_mean, " s")
	fmt.Println("exact amot. : ", 1000*exact_mean/float64(batch), " ms")
	fmt.Println("min err : ", min_mean)
	fmt.Println("avg err : ", avg_mean)
}
