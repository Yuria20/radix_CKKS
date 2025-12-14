package main

import (
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

func IntToFloatExact(x *big.Int) *big.Float {
	prec := uint(x.BitLen())
	if prec == 0 {
		prec = 1
	} // x==0 보호
	return new(big.Float).SetPrec(prec).SetMode(big.ToNearestEven).SetInt(x)
}

func Test_P384_baseREDC(Bbox Bk_Arithmetic_toolbox) {

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

	}

	var err error
	var params ckks.Parameters
	if Bbox.base == 16 {
		params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            LogN,                      // Log2 of the ring degree
			LogQ:            LogQ,                      // Log2 of the ciphertext modulus
			LogP:            []int{52, 52, 52, 52, 52}, // Log2 of the key-switch auxiliary prime moduli
			LogDefaultScale: LogDefaultScale,           // Log2 of the scale
			Xs:              ring.Ternary{H: 32768},
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
	batch := int(float64(params.MaxSlots()) / float64(slice_length))

	//fmt.Println(Compute_Shift_Index(params, slice_length, params.MaxSlots()/slice_length))
	galEls := []uint64{
		// The galois element for the cyclic rotations by 55positions to the left.
		params.GaloisElement(rot),
		params.GaloisElement(slice_length / 4 * batch),
		params.GaloisElement(96 * batch),
		// The galois element for the complex conjugatation.
		params.GaloisElementForComplexConjugation(),
		params.GaloisElement(params.MaxSlots() / 2),
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

	//batch := int(float64(params.MaxSlots()) / float64(slice_length))
	//mod1, _ := GenPrime(1024)
	//mod2, _ := GenPrime(1024)
	//mod := big.NewInt(0.0).Mul(mod1, mod2)

	//fmt.Println(mod)

	//primeStr := "28057595600765895243607060750207179866082002587411688425150475127313997461021553774559567496072065457119562022632013018804912241256785149185912900393066214448895544554146993039771469568911245723703124380899723110776536815504192700147941524695921628054591696543610699669469187758736947007037477107117417644889380003346927707969870087324621347765390422879505467815523091080136196634558818534660716563492607325360006900142335935990378449285144564176765791638185088071563116084995998490999991893631560515168011944277666697044347488763409505184936079541099564139815133305098656604674174396469051359407274511309384955839319"
	// p-384
	primeStr := "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"
	mod := new(big.Int)
	mod.SetString(primeStr, 10)

	//fmt.Println("hit!")

	fmt.Printf("Large Integer parameter : slice_length=%d, bit_length=%d, batch=%d, lazycarry_iter=%d, carry_iter=%d", slice_length, bit_length, batch, Bbox.lazy_iter[bit_length], Bbox.carry_iter[bit_length])
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

	// 로드(동일한 params 필요)
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

	if file, _ := FileExists("P384_baseREDC"); file == false {
		baseREDCmat := Generate_baseREDCmat(base, 96*2, mod)
		inf_norm := ComputeMaxNorm(baseREDCmat, base, mod)

		baseREDCmat = MatrixPadding(baseREDCmat, slice_length, params.MaxSlots())
		//fmt.Println(baseREDCmat)
		//Normalized_DFT = MatrixPadding(Normalized_DFT, slice_length, params.MaxSlots())
		Twisted_baseREDCmat := TwistedMatrix(baseREDCmat, slice_length, params.MaxSlots())
		Diag_baseREDCmat := Row_To_Diagonal(Twisted_baseREDCmat)
		Twisted_baseREDCmat = nil
		runtime.GC()

		Plain_baseREDCmat := BSGS_plain_Gen(params, slice_length, batch, Lv-1, Diag_baseREDCmat, cc)
		//fmt.Println(Plain_baseREDCmat[0].IsBatched)

		log_inf_norm, _ := inf_norm.Float64()
		log_inf_norm = math.Log2(log_inf_norm)
		cond_iter := int(math.Ceil(log_inf_norm))
		fmt.Println(inf_norm, cond_iter)

		//fmt.Println(baseREDCmat)
		fmt.Println(inf_norm)
		baseREDCmat = nil

		if err := SavePlaintextMap("P384_baseREDC", Plain_baseREDCmat); err != nil {
			panic(err)
		}

		if err != nil {
			panic(err)
		}
		Diag_baseREDCmat = nil
		runtime.GC()
	}

	Plain_baseREDCmat, _ := LoadPlaintextMap("P384_baseREDC", params)

	fmt.Println("hit!")
	lazy_time := make([]float64, 5)
	exact_time := make([]float64, 6)
	min_err := make([]float64, 6)
	avg_err := make([]float64, 6)

	for k := 0; k < 5; k++ {

		values1_big := make([]*big.Int, batch)
		values2_big := make([]*big.Int, batch)
		values1_int := make([]int64, params.MaxSlots())
		values2_int := make([]int64, params.MaxSlots())
		valueN_int := make([]int64, params.MaxSlots())
		value2N_int := make([]int64, params.MaxSlots())
		valueN_cpx := make([]complex128, params.MaxSlots())
		value2N_cpx := make([]complex128, params.MaxSlots())

		//fmt.Println(mod)

		//fmt.Println("mod : ", mod)
		//fmt.Println("mod : ", big.NewInt(0.0).Mul(mod, big.NewInt(2.0)))

		for i := 0; i < batch; i++ {
			valueN_int = DecomposeInt64Base(mod, int64(base), slice_length)
			value2N_int = DecomposeInt64Base(big.NewInt(0.0).Mul(mod, big.NewInt(2.0)), int64(base), slice_length)

			for j := 0; j < slice_length; j++ {
				valueN_cpx[i*slice_length+j] = complex(float64(valueN_int[j]), 0)
				value2N_cpx[i*slice_length+j] = complex(float64(value2N_int[j]), 0)
			}
		}

		//fmt.Println("mod : ", mod)
		//fmt.Println("mod : ", big.NewInt(0.0).Mul(mod, big.NewInt(2.0)))

		for i := 0; i < batch; i++ {
			values1_big[i], _ = RandomBigInt_mod(mod)
			values2_big[i], _ = RandomBigInt_mod(mod)

			//a := "4782023829109761722472798430788939709906878331370176896010083845243703667195716390422561705445533716818391830077633304753132856536606244374064944512442211381330539364108232275298850362727549976192629747912046784031567511263006900319250524374688664238773001156851250239633673553001508095167591904915817489764589351635139479019675126045901894885451493124291867257091031896522945041347311189338409488073779555089803866665026242543778177792656358643057635883902239807286315547153333440779555696639536679624979650640159884619749574126162336140961414485017905961651285080634588109808594547887738051443179450842851107500716"
			//b := "19678016127222884972896136671114436891386376113641182645059574433738579923312170080874082864904953960012218474080417368012685499252287949068873353300714173855363692522094013903887774739716768484345999909118549570553848868492994897969204774609536119948282572239490594821884855210590256668294703266748158241793858295590117880259433947370570057179737892899365292999291517095321844583672221482409397774183311110719404687022511982996346752941192811208660024038141714391450163534946616705033203844337993762068661721048373156446562875197068880610230787544614695226791113113772406507923272847292085602070442977711080899909660"
			//values1_big[i].SetString(a, 10)
			//values2_big[i].SetString(b, 10)

			//values1_big[i] = big.NewInt(1)
			//values2_big[i] = big.NewInt(1)

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
		// TODO : 교차 패킹으로 수정 필요.
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
		//fmt.Println(value1[:8])
		//fmt.Println(value2[:8])

		// Encrypt value 1

		twisted_values1 := TwistedVec(value1, slice_length, params.MaxSlots())
		twisted_values2 := TwistedVec(value2, slice_length, params.MaxSlots())

		plaintext1 := ckks.NewPlaintext(params, Lv-1)
		if err := encoder.Encode(twisted_values1, plaintext1); err != nil {
			panic(err)
		}

		//fmt.Println(Lv)

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

		// Convert Complex to Int
		for i := 0; i < params.MaxSlots(); i++ {
			valuesWant_int[i] = big.NewInt(int64(math.Round(real(valuesWant_cpx[i]))))
		}
		//fmt.Println(valuesWant_int[:8])

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

		//==============================================
		//=== 2) Mult Large Integer : MultPoly =========
		//==============================================
		if k == 0 {
			fmt.Println("warm-up round")
		} else {
			fmt.Printf("%d th iteration\n", k)
		}
		fmt.Println("Multiplication Start!")

		mult_time := time.Now()
		var ciphertext *rlwe.Ciphertext
		//var debug *rlwe.Ciphertext

		// T = A*B
		var T *rlwe.Ciphertext

		// large mult
		T = ExactMult(ciphertext1, ciphertext2, Plain_DFT, Plain_InvDFT, Bbox, cc)

		//fmt.Println("T lv : ", T.Level())

		fmt.Print("Cleaning... ")
		start := time.Now()
		T = Cleaning_with_vectorized_evaluation(T, Bbox, cc)
		elapsed_time := time.Since(start)
		fmt.Println(elapsed_time)

		//ciphertext = T.CopyNew()
		//Plain_DFT = nil
		//Plain_InvDFT = nil
		runtime.GC()

		fmt.Print("baseREDC... ")
		start = time.Now()
		temp := BSGS_FFT(batch, slice_length, T, Plain_baseREDCmat, cc)
		elapsed_time = time.Since(start)
		fmt.Println(elapsed_time)

		err = eval.Rotate(temp, 96*batch, temp)
		if err != nil {
			fmt.Println(err)
		}
		masking := make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			for j := 0; j < slice_length; j++ {
				if j < 96 {
					masking[i*slice_length+j] = complex(1.0, 0.0)
				}
			}
		}
		masking = TwistedVec(masking, slice_length, params.MaxSlots())

		pt := ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(params.Q()[temp.Level()])
		encoder.Encode(masking, pt)

		eval.MulRelin(T, masking, T)
		err = eval.Rescale(T, T)
		if err != nil {
			fmt.Println(err)
		}

		eval.Add(temp, T, T)
		ciphertext = T.CopyNew()
		runtime.GC()

		for iter := 0; iter < 3; iter++ {
			fmt.Printf("LazyCarry... ")
			lazy_time := time.Now()
			ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
			elapsed_time := time.Since(lazy_time)
			fmt.Println(elapsed_time)
		}

		ciphertext = LazyCarry2Carry(ciphertext, Bbox, cc)
		ciphertext = Cleaning_with_vectorized_evaluation(ciphertext, Bbox, cc)
		//fmt.Println(ciphertext.Level())

		fmt.Print("baseREDC... ")
		start = time.Now()
		temp = BSGS_FFT(batch, slice_length, ciphertext, Plain_baseREDCmat, cc)
		elapsed_time = time.Since(start)
		fmt.Println(elapsed_time)

		err = eval.Rotate(temp, 96*batch, temp)
		if err != nil {
			fmt.Println(err)
		}
		masking = make([]complex128, params.MaxSlots())
		for i := 0; i < batch; i++ {
			for j := 0; j < slice_length; j++ {
				if j < 96 {
					masking[i*slice_length+j] = complex(1.0, 0.0)
				}
			}
		}
		masking = TwistedVec(masking, slice_length, params.MaxSlots())

		pt = ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(params.Q()[temp.Level()])
		encoder.Encode(masking, pt)

		eval.MulRelin(ciphertext, masking, ciphertext)
		err = eval.Rescale(ciphertext, ciphertext)
		if err != nil {
			fmt.Println(err)
		}

		eval.Add(temp, ciphertext, ciphertext)

		runtime.GC()

		for iter := 0; iter < 2; iter++ {
			fmt.Printf("LazyCarry... ")
			lazy_time := time.Now()
			ciphertext = LazyCarry_with_vectorized_evaluation(base, slice_length, ciphertext, cc)
			elapsed_time := time.Since(lazy_time)
			fmt.Println(elapsed_time)
		}

		lazy_elapsed := time.Since(mult_time)

		ciphertext = LazyCarry2Carry(ciphertext, Bbox, cc)

		var cipherN *rlwe.Ciphertext
		valueN_cpx = TwistedVec(valueN_cpx, slice_length, params.MaxSlots())
		plain_N := ckks.NewPlaintext(params, ciphertext.Level())
		//plain_N_rad.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
		encoder.Encode(valueN_cpx, plain_N)
		cipherN, err = cc.encryptor.EncryptNew(plain_N)
		if err != nil {
			fmt.Println(err)
		}

		ciphertext = ConditonalSub(ciphertext, cipherN, Bbox, cc, "p384")
		//fmt.Println(ciphertext.Level())

		total_elapsed := time.Since(mult_time)
		fmt.Println()
		//==============================================
		//=== 53 Decrypt and Validation ================
		//==============================================

		// Decrypt
		//if ciphertext.Level() != 0 {
		//	eval.DropLevel(ciphertext, ciphertext.Level())
		//}
		//ciphertext = debug.CopyNew()

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
		fmt.Printf("Large Integer accuracy : %0.4f\n", acc/float64(batch))

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
