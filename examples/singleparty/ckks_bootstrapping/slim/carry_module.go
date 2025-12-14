package main

import (
	"fmt"
	"math"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func LazyCarry(base int, slice_length int, ciphertext *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	var err error
	var ct_mod, ct_Q *rlwe.Ciphertext
	var temp *rlwe.Ciphertext
	temp = ciphertext.CopyNew()
	//value := make([]complex128, params.MaxSlots())

	//fmt.Print("Before Lazycarry")
	//PrintDebug(*params, ciphertext, value, cc.decryptor, cc.encoder)

	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}

	//fmt.Println("(In)LogScale : ", ciphertext.LogScale())

	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	//fmt.Println("(StC)LogScale : ", ciphertext.LogScale())
	//fmt.Println("(ScaleDown)LogScale : ", ciphertext.Level())

	// Step 2: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(base))))

	// Step 3 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext

	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}
	//fmt.Println("(CtS)LogScale : ", ciphertext.LogScale())
	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalExp (Homomorphic modular reduction)
	start := time.Now()
	if imag, err = btp.EvalModAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}
	elapsed := time.Since(start)
	fmt.Println("EvalExp(cos) :", elapsed)

	if real, err = btp.EvalElseAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = btp.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if err = btp.Evaluator.Add(real, imag, ciphertext); err != nil {
		panic(err)
	}
	//fmt.Println("(EvalExp)LogScale : ", ciphertext.LogScale())

	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.Polynomial
	var eval_poly bignum.Polynomial

	//start := time.Now()
	eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_exp2id(base, 2*base-1), nil)
	//t := time.Since(start)
	//fmt.Println(t)

	if polys = polynomial.NewPolynomial(eval_poly); err != nil {
		panic(err)
	}

	// [z]_b is computed
	start = time.Now()
	ct_mod, _ = polyEval.Evaluate(ciphertext, polys, temp.Scale)
	elapsed = time.Since(start)
	fmt.Println("LUT :", elapsed)

	//return ct_mod
	//fmt.Println("(LUT)LogScale : ", ct_mod.LogScale())

	// Q_b(z) is computed
	//fmt.Println(temp.Scale.Float64())
	//fmt.Println(ct_mod.Scale.Float64())
	ct_Q, _ = eval.SubNew(temp, ct_mod)

	eval.Rotate(ct_Q, -1*params.MaxSlots()/slice_length, ct_Q)
	//eval.Rotate(ct_Q, -1, ct_Q)

	// make masking vector
	masking := make([]complex128, params.MaxSlots())
	for i := 0; i < params.MaxSlots()/slice_length; i++ {
		for j := 0; j < slice_length; j++ {
			if j == 0 {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			} else if j < slice_length/2 {
				masking[i*slice_length+j] = complex(1.0/float64(base)*float64(params.DefaultScale().Float64())/ct_Q.Scale.Float64(), 0.0)
			} else {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			}
		}
	}

	masking = TwistedVec(masking, slice_length, params.MaxSlots())
	pt := ckks.NewPlaintext(*params, params.MaxLevel())
	pt.Scale = rlwe.NewScale(cc.params.Q()[ct_Q.Level()])
	cc.encoder.Encode(masking, pt)

	eval.MulRelin(ct_Q, pt, ct_Q)
	eval.Rescale(ct_Q, ct_Q)
	ct_Q.Scale.Value.Set(big.NewFloat(float64(params.DefaultScale().Float64())))

	eval.SetScale(ct_mod, params.DefaultScale())
	eval.Add(ct_mod, ct_Q, ciphertext)

	//eval.MulRelin(ciphertext, masking, ciphertext)
	//eval.Rescale(ciphertext, ciphertext)
	//fmt.Println("(mod + Q )LogScale : ", ct_mod.LogScale())

	//fmt.Print("after LazyCarry")
	//PrintDebug(*params, ciphertext, value, cc.decryptor, cc.encoder)

	//fmt.Println("Level after bootstrapping : ", ciphertext.Level())
	//fmt.Println("Done")

	return ciphertext
}

func LazyCarry_with_vectorized_evaluation(base int, slice_length int, ciphertext *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	var err error
	var ct_mod, ct_Q *rlwe.Ciphertext
	var temp *rlwe.Ciphertext
	temp = ciphertext.CopyNew()
	//value := make([]complex128, params.MaxSlots())

	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}

	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	fmt.Println("base : ", base)
	btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(base))))

	// Step 3 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext

	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalExp (Homomorphic modular reduction)
	if ciphertext, err = btp.EvalExpAndScale_with_vectroized_evaluation(real, 2*math.Pi/complex(float64(base), 0.0), params); err != nil {
		panic(err)
	}

	polyEval := polynomial.NewEvaluator(*params, eval)
	var polys polynomial.PolynomialVector
	eval_poly_vec := make([]bignum.Polynomial, 2)

	eval_poly_vec[0] = bignum.NewPolynomial(0, HermiteInterpolation_exp2id(base, 2*base-1), nil)
	eval_poly_vec[1] = bignum.NewPolynomial(0, make([]complex128, 2*base), nil)

	mapping := make(map[int][]int, 2)
	for i := 0; i < 2; i++ {
		mapping[i] = make([]int, params.MaxSlots()/2)
		for j := 0; j < params.MaxSlots()/2; j++ {
			mapping[i][j] = j + i*params.MaxSlots()/2
		}
	}

	if polys, err = polynomial.NewPolynomialVector(eval_poly_vec, mapping); err != nil {
		panic(err)
	}

	// [z]_b is computed
	ct_mod, err = polyEval.Evaluate(ciphertext, polys, temp.Scale)
	if err != nil {
		fmt.Println(err)
	}

	// Q_b(z) is computed
	ct_Q, _ = eval.SubNew(temp, ct_mod)
	eval.Rotate(ct_Q, -1*params.MaxSlots()/slice_length, ct_Q)

	// make masking vector
	masking := make([]complex128, params.MaxSlots())
	for i := 0; i < params.MaxSlots()/slice_length; i++ {
		for j := 0; j < slice_length; j++ {
			if j == 0 {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			} else if j < slice_length/2 {
				masking[i*slice_length+j] = complex(1.0/float64(base)*float64(params.DefaultScale().Float64())/ct_Q.Scale.Float64(), 0.0)
			} else {
				masking[i*slice_length+j] = complex(0.0, 0.0)
			}
		}
	}

	masking = TwistedVec(masking, slice_length, params.MaxSlots())
	pt := ckks.NewPlaintext(*params, params.MaxLevel())
	pt.Scale = rlwe.NewScale(cc.params.Q()[ct_Q.Level()])
	cc.encoder.Encode(masking, pt)

	eval.MulRelin(ct_Q, pt, ct_Q)
	eval.Rescale(ct_Q, ct_Q)
	ct_Q.Scale.Value.Set(big.NewFloat(float64(params.DefaultScale().Float64())))
	//PrintDebug(slice_length, *params, ct_Q, value, cc.decryptor, cc.encoder)

	eval.SetScale(ct_mod, params.DefaultScale())
	eval.Add(ct_mod, ct_Q, ciphertext)

	return ciphertext
}

func IntToSymbol(base int, slice_length int, ciphertext *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	var err error

	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}
	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	//fmt.Println(ciphertext.Level())

	value := make([]complex128, params.MaxSlots())

	for i := 0; i < params.MaxSlots(); i++ {
		value[i] = complex(float64(base)/float64(2*base-1)*params.DefaultScale().Float64()/ciphertext.Scale.Float64(), 0.0)
	}

	pt := ckks.NewPlaintext(*params, params.MaxLevel())
	pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
	cc.encoder.Encode(value, pt)

	eval.MulRelin(ciphertext, pt, ciphertext)
	eval.Rescale(ciphertext, ciphertext)
	//btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(2*base-1))))
	//fmt.Println(ciphertext.Scale)
	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalMod (Homomorphic modular reduction)
	if imag, err = btp.EvalModAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	if real, err = btp.EvalElseAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = btp.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if err = btp.Evaluator.Add(real, imag, ciphertext); err != nil {
		panic(err)
	}

	//fmt.Println("Eval")
	//value := make([]complex128, params.MaxSlots())
	//PrintDebug(slice_length, *params, ciphertext, value, cc.decryptor, cc.encoder)

	//opeval := ckks.NewEvaluator(*params, evk)
	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.Polynomial
	var eval_poly bignum.Polynomial

	//start := time.Now()
	//eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_exp2id(base, 2*base-1), nil)
	eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_exp2symbol(2*base-1, 4*base-3, base), nil)
	//t := time.Since(start)
	//fmt.Println(t)

	if polys = polynomial.NewPolynomial(eval_poly); err != nil {
		panic(err)
	}

	ciphertext, _ = polyEval.Evaluate(ciphertext, polys, params.DefaultScale())

	//fmt.Println(ciphertext.Scale.Value.Float64())

	//fmt.Println(ciphertext.Level())

	return ciphertext
}

func IntToSymbol_with_vectorized_evaluation(base int, slice_length int, ciphertext *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	var err error

	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}
	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	//fmt.Println(ciphertext.Level())

	value := make([]complex128, params.MaxSlots())

	for i := 0; i < params.MaxSlots(); i++ {
		value[i] = complex(float64(base)/float64(2*base-1)*params.DefaultScale().Float64()/ciphertext.Scale.Float64(), 0.0)
	}

	pt := ckks.NewPlaintext(*params, params.MaxLevel())
	pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
	cc.encoder.Encode(value, pt)

	eval.MulRelin(ciphertext, pt, ciphertext)
	eval.Rescale(ciphertext, ciphertext)
	//btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(2*base-1))))
	//fmt.Println(ciphertext.Scale)
	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalExp (Homomorphic modular reduction)
	//start := time.Now()
	if ciphertext, err = btp.EvalExpAndScale_with_vectroized_evaluation(real, 2*math.Pi/complex(float64(base), 0.0), params); err != nil {
		panic(err)
	}
	//elapsed := time.Since(start)
	//fmt.Println("EvalExp(cos) :", elapsed)

	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.PolynomialVector
	eval_poly_vec := make([]bignum.Polynomial, 2)

	//start := time.Now()
	eval_poly_vec[0] = bignum.NewPolynomial(0, HermiteInterpolation_exp2symbol(2*base-1, 4*base-3, base), nil)
	eval_poly_vec[1] = bignum.NewPolynomial(0, make([]complex128, 4*base-2), nil)

	mapping := make(map[int][]int, 2)
	for i := 0; i < 2; i++ {
		mapping[i] = make([]int, params.MaxSlots()/2)
		for j := 0; j < params.MaxSlots()/2; j++ {
			mapping[i][j] = j + i*params.MaxSlots()/2
		}
	}
	//t := time.Since(start)
	//fmt.Println(t)

	if polys, err = polynomial.NewPolynomialVector(eval_poly_vec, mapping); err != nil {
		panic(err)
	}

	ciphertext, _ = polyEval.Evaluate(ciphertext, polys, params.DefaultScale())

	//fmt.Println(ciphertext.Scale.Value.Float64())

	//fmt.Println(ciphertext.Level())

	return ciphertext
}

func IntToSymbol_negate(base int, slice_length int, ciphertext *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	var err error

	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}
	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	//fmt.Println(ciphertext.Level())

	value := make([]complex128, params.MaxSlots())

	for i := 0; i < params.MaxSlots(); i++ {
		value[i] = complex(float64(base)/float64(2*base-1)*params.DefaultScale().Float64()/ciphertext.Scale.Float64(), 0.0)
	}

	pt := ckks.NewPlaintext(*params, params.MaxLevel())
	pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
	cc.encoder.Encode(value, pt)

	eval.MulRelin(ciphertext, pt, ciphertext)
	eval.Rescale(ciphertext, ciphertext)
	//btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(2*base-1))))
	//fmt.Println(ciphertext.Scale)
	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalMod (Homomorphic modular reduction)
	if imag, err = btp.EvalModAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	if real, err = btp.EvalElseAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = btp.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if err = btp.Evaluator.Add(real, imag, ciphertext); err != nil {
		panic(err)
	}

	//fmt.Println("Eval")
	//value := make([]complex128, params.MaxSlots())
	//PrintDebug(slice_length, *params, ciphertext, value, cc.decryptor, cc.encoder)

	//opeval := ckks.NewEvaluator(*params, evk)
	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.Polynomial
	var eval_poly bignum.Polynomial

	//start := time.Now()
	//eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_exp2id(base, 2*base-1), nil)
	eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_exp2negatesymbol(2*base-1, 4*base-3, base), nil)
	//t := time.Since(start)
	//fmt.Println(t)

	if polys = polynomial.NewPolynomial(eval_poly); err != nil {
		panic(err)
	}

	ciphertext, _ = polyEval.Evaluate(ciphertext, polys, params.DefaultScale())

	//fmt.Println(ciphertext.Scale.Value.Float64())

	//fmt.Println(ciphertext.Level())

	return ciphertext
}

func IntToSymbol_negate_with_vectorized_evaluation(base int, slice_length int, ciphertext *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	var err error

	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}
	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	//fmt.Println(ciphertext.Level())

	value := make([]complex128, params.MaxSlots())

	for i := 0; i < params.MaxSlots(); i++ {
		value[i] = complex(float64(base)/float64(2*base-1)*params.DefaultScale().Float64()/ciphertext.Scale.Float64(), 0.0)
	}

	pt := ckks.NewPlaintext(*params, params.MaxLevel())
	pt.Scale = rlwe.NewScale(cc.params.Q()[ciphertext.Level()])
	cc.encoder.Encode(value, pt)

	eval.MulRelin(ciphertext, pt, ciphertext)
	eval.Rescale(ciphertext, ciphertext)
	//btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(2*base-1))))
	//fmt.Println(ciphertext.Scale)
	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalExp (Homomorphic modular reduction)
	//start := time.Now()
	/*
	 */
	if ciphertext, err = btp.EvalExpAndScale_with_vectroized_evaluation(real, 2*math.Pi/complex(float64(base), 0.0), params); err != nil {
		panic(err)
	}
	//elapsed := time.Since(start)
	//fmt.Println("EvalExp(cos) :", elapsed)

	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.PolynomialVector
	eval_poly_vec := make([]bignum.Polynomial, 2)

	//start := time.Now()
	eval_poly_vec[0] = bignum.NewPolynomial(0, HermiteInterpolation_exp2negatesymbol(2*base-1, 4*base-3, base), nil)
	eval_poly_vec[1] = bignum.NewPolynomial(0, make([]complex128, 4*base-2), nil)

	mapping := make(map[int][]int, 2)
	mapping[0] = make([]int, params.MaxSlots()/2)
	for j := 0; j < params.MaxSlots()/2+params.MaxSlots()/slice_length; j++ {
		mapping[0][j] = j
	}

	mapping[1] = make([]int, params.MaxSlots()/2-params.MaxSlots()/(slice_length))
	for j := 0; j < params.MaxSlots()/2-params.MaxSlots()/slice_length; j++ {
		mapping[1][j] = j + params.MaxSlots()/2 + params.MaxSlots()/slice_length
	}

	//t := time.Since(start)
	//fmt.Println(t)

	if polys, err = polynomial.NewPolynomialVector(eval_poly_vec, mapping); err != nil {
		panic(err)
	}

	ciphertext, _ = polyEval.Evaluate(ciphertext, polys, params.DefaultScale())

	//fmt.Println(ciphertext.Scale.Value.Float64())

	//fmt.Println(ciphertext.Level())

	return ciphertext
}

func LazyCarryToCarry(base int, log_slice_half int, slice_length int, ct_symbol *rlwe.Ciphertext, ct_carry *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	//btp := cc.btp
	//value := make([]complex128, params.MaxSlots())
	encoder := cc.encoder

	//ct_symbol = CleaningSymbol(base, slice_length, ct_symbol, cc)
	prec := 30

	for i := 0; i < log_slice_half; i++ {
		var ct_rot *rlwe.Ciphertext
		var err error
		//ct_rot, _ = eval.RotateNew(ct_symbol, -1*(1<<i))
		ct_rot, err = eval.RotateNew(ct_symbol, -1*(1<<i)*(params.MaxSlots())/slice_length)
		//ct_rot, err = eval.RotateNew(ct_symbol, (1<<i)*(params.MaxSlots())/slice_length)
		if err != nil {
			fmt.Println(err)
		}

		prec -= 2

		carry_part1 := func() {
			var temp, temp2 *rlwe.Ciphertext
			temp = ct_symbol.CopyNew()

			// compute (x+conj(x))
			eval.Conjugate(temp, temp)
			eval.Add(temp, ct_symbol, temp)

			// compute y-x
			temp2, _ = eval.SubNew(ct_rot, ct_symbol)

			// compute (x+conj(x))*(y-x)
			eval.MulRelin(temp, temp2, temp)
			eval.Rescale(temp, temp)

			eval.Add(temp, ct_symbol, ct_symbol)

			if ct_symbol.Level() <= 4 {
				if ct_symbol.Level() <= 4 {
					fmt.Print("Symbol Cleaning... ")
					start := time.Now()
					ct_symbol = CleaningSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
					elapsed_time := time.Since(start)
					fmt.Println(elapsed_time)

					prec = 30
				}
			}
		}
		carry_part1()
	}

	if ct_symbol.Level() < 5 {
		fmt.Print("Symbol Cleaning... ")
		start := time.Now()
		ct_symbol = CleaningSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
		elapsed_time := time.Since(start)
		fmt.Println(elapsed_time)

		prec = 30
	}

	carry_part2 := func() {
		var temp *rlwe.Ciphertext
		temp = ct_symbol.CopyNew()

		eval.Conjugate(temp, temp)

		eval.Sub(temp, ct_symbol, temp)

		values := make([]complex128, params.MaxSlots())
		for i := 0; i < params.MaxSlots(); i++ {
			values[i] = complex(0, float64(params.DefaultScale().Float64())*0.5/ct_symbol.Scale.Float64())
		}

		pt := ckks.NewPlaintext(*params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(params.Q()[temp.Level()])
		encoder.Encode(values, pt)

		eval.MulRelin(temp, pt, ct_symbol)
		eval.Rescale(ct_symbol, ct_symbol)
		ct_symbol.Scale.Value.Set(big.NewFloat(float64(params.DefaultScale().Float64())))

		prec -= 1
	}

	carry_part2()

	var temp *rlwe.Ciphertext

	temp = ct_symbol.CopyNew()
	eval.MulRelin(ct_symbol, base, temp)

	eval.Sub(ct_carry, temp, ct_carry)

	err := eval.Rotate(ct_symbol, -1*params.MaxSlots()/slice_length, ct_symbol)
	if err != nil {
		fmt.Println(err)
	}
	eval.Add(ct_carry, ct_symbol, ct_carry)

	return ct_carry
}

func LazyCarryToCarry_negate(base int, log_slice_half int, slice_length int, ct_symbol *rlwe.Ciphertext, ct_carry *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	//btp := cc.btp
	//value := make([]complex128, params.MaxSlots())
	encoder := cc.encoder

	//ct_symbol = CleaningSymbol(base, slice_length, ct_symbol, cc)
	prec := 30

	for i := 0; i < log_slice_half; i++ {
		var ct_rot *rlwe.Ciphertext
		var err error
		//ct_rot, _ = eval.RotateNew(ct_symbol, -1*(1<<i))
		ct_rot, err = eval.RotateNew(ct_symbol, -1*(1<<i)*(params.MaxSlots())/slice_length)
		//ct_rot, err = eval.RotateNew(ct_symbol, (1<<i)*(params.MaxSlots())/slice_length)
		if err != nil {
			fmt.Println(err)
		}
		//fmt.Println("ct_symbol ")
		//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

		//fmt.Println("rot : ", -1*(1<<i))
		//PrintDebug(slice_length, *params, ct_rot, value, cc.decryptor, cc.encoder)

		prec -= 2

		carry_part1 := func() {
			var temp, temp2 *rlwe.Ciphertext
			temp = ct_symbol.CopyNew()

			// compute (x+conj(x))
			eval.Conjugate(temp, temp)
			eval.Add(temp, ct_symbol, temp)

			// compute y-x
			temp2, _ = eval.SubNew(ct_rot, ct_symbol)

			// compute (x+conj(x))*(y-x)
			eval.MulRelin(temp, temp2, temp)
			eval.Rescale(temp, temp)

			eval.Add(temp, ct_symbol, ct_symbol)

			if ct_symbol.Level() <= 4 {
				if ct_symbol.Level() <= 4 {
					//ct_symbol = CleaningSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
					ct_symbol = CleaningSymbol(base, slice_length, ct_symbol, cc)

					prec = 30
				}
			}
		}
		carry_part1()
	}

	if ct_symbol.Level() < 5 {
		//ct_symbol = CleaningSymbol_with_vectorized_evaluation(base, slice_length, ct_symbol, cc)
		ct_symbol = CleaningSymbol(base, slice_length, ct_symbol, cc)
		prec = 30
	}

	//fmt.Println("symbol ciphertext")
	//ComputePrecSymbol(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

	carry_part2 := func() {
		var temp *rlwe.Ciphertext
		temp = ct_symbol.CopyNew()

		// compute 1/2*i*(conj(x)-x)
		eval.Conjugate(temp, temp)
		//fmt.Println("conj(x)")
		//PrintDebug(slice_length, *params, temp, value, cc.decryptor, cc.encoder)
		eval.Sub(temp, ct_symbol, temp)
		//fmt.Println("conj(x)-x")
		//PrintDebug(slice_length, *params, temp, value, cc.decryptor, cc.encoder)

		values := make([]complex128, params.MaxSlots())
		for i := 0; i < params.MaxSlots(); i++ {
			values[i] = complex(0, float64(params.DefaultScale().Float64())*0.5/ct_symbol.Scale.Float64())
		}

		pt := ckks.NewPlaintext(*params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(params.Q()[temp.Level()])
		encoder.Encode(values, pt)

		//fmt.Println(float64())
		eval.MulRelin(temp, pt, ct_symbol)
		eval.Rescale(ct_symbol, ct_symbol)
		ct_symbol.Scale.Value.Set(big.NewFloat(float64(params.DefaultScale().Float64())))

		//fmt.Println("i/2*(conj(x)-x)")
		//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

		prec -= 1
	}

	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

	carry_part2()

	//fmt.Println("carry2")
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	//fmt.Println("symbol ciphertext")
	//ComputePrecSymbol(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	//PrintDebug(*params, ct_symbol, value, cc.decryptor, cc.encoder)

	var temp *rlwe.Ciphertext

	temp = ct_symbol.CopyNew()
	eval.MulRelin(ct_symbol, base, temp)

	//fmt.Println("B * s")
	//PrintDebug(slice_length, *params, temp, value, cc.decryptor, cc.encoder)

	eval.Add(ct_carry, temp, ct_carry)

	//fmt.Println("x - B * s")
	//PrintDebug(slice_length, *params, ct_carry, value, cc.decryptor, cc.encoder)

	err := eval.Rotate(ct_symbol, -1*params.MaxSlots()/slice_length, ct_symbol)
	if err != nil {
		fmt.Println(err)
	}
	eval.Sub(ct_carry, ct_symbol, ct_carry)

	//fmt.Println("x - B * s + rho(s)")
	//PrintDebug(slice_length, *params, ct_carry, value, cc.decryptor, cc.encoder)

	return ct_carry
}

func CleaningSymbol(base int, slice_length int, ct_symbol *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	// Step 0: Some circuit in the slots domain
	var err error

	//value := make([]complex128, params.MaxSlots())
	//eval.Add(ct_symbol, 0.5, ct_symbol)
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

	if ct_symbol.Level() > 4 {
		btp.DropLevel(ct_symbol, ct_symbol.Level()-4)
	}

	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ct_symbol, err = btp.SlotsToCoeffs(ct_symbol, nil); err != nil {
		panic(err)
	}

	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ct_symbol, _, err = btp.ScaleDown(ct_symbol); err != nil {
	//	panic(err)
	//}
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	btp.SetScale(ct_symbol, rlwe.NewScale(float64(params.Q()[0])/(float64(base))))
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

	// Step 4 : Extend the basis from q to Q
	if ct_symbol, err = btp.ModUp(ct_symbol); err != nil {
		panic(err)
	}
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	var real2, imag2 *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ct_symbol); err != nil {
		panic(err)
	}

	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)
	eval.Mul(real, 2.0, real)
	eval.Mul(imag, 2.0, imag)
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)
	//fmt.Println("After CoeffsToSlot Scale : ", real.Scale.Uint64())

	// Step 6 : EvalMod (Homomorphic modular reduction)

	if imag2, err = btp.EvalModAndScale(imag, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	if imag, err = btp.EvalElseAndScale(imag, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	if real2, err = btp.EvalModAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	if real, err = btp.EvalElseAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = btp.Evaluator.Mul(imag2, 1i, imag2); err != nil {
		panic(err)
	}

	if err = btp.Evaluator.Add(imag, imag2, imag); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = btp.Evaluator.Mul(real2, 1i, real2); err != nil {
		panic(err)
	}

	if err = btp.Evaluator.Add(real, real2, real); err != nil {
		panic(err)
	}

	//fmt.Println("dddd")
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.Polynomial
	var eval_poly bignum.Polynomial
	eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_symbol2symbol(), nil)

	if polys = polynomial.NewPolynomial(eval_poly); err != nil {
		panic(err)
	}

	// [z]_b is computed
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	real, _ = polyEval.Evaluate(real, polys, params.DefaultScale())
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)

	eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_symbol2symbol_imag(), nil)

	if polys = polynomial.NewPolynomial(eval_poly); err != nil {
		panic(err)
	}

	// [z]_b is computed
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)
	imag, _ = polyEval.Evaluate(imag, polys, params.DefaultScale())
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	//value := make([]complex128, params.MaxSlots())
	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)

	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(*params, imag, value, cc.decryptor, cc.encoder)

	//fmt.Println("Scale : ", real.Scale.Uint64())
	//fmt.Println("Scale : ", imag.Scale.Uint64())

	// Recombines the real and imaginary part
	//if err = btp.Evaluator.Mul(imag, 1i, imag); err != nil {
	//	panic(err)
	//}

	if err = btp.Evaluator.Add(real, imag, ct_symbol); err != nil {
		panic(err)
	}

	//fmt.Println("Scale : ", ciphertext.Scale.Uint64())

	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(*params, ct_symbol, value, cc.decryptor, cc.encoder)

	return ct_symbol
}

func CleaningSymbol_with_vectorized_evaluation(base int, slice_length int, ct_symbol *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	params := cc.params
	eval := cc.eval
	btp := cc.btp

	// Step 0: Some circuit in the slots domain
	var err error

	//value := make([]complex128, params.MaxSlots())
	//eval.Add(ct_symbol, 0.5, ct_symbol)
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

	if ct_symbol.Level() > 4 {
		btp.DropLevel(ct_symbol, ct_symbol.Level()-4)
	}

	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ct_symbol, err = btp.SlotsToCoeffs(ct_symbol, nil); err != nil {
		panic(err)
	}

	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ct_symbol, _, err = btp.ScaleDown(ct_symbol); err != nil {
	//	panic(err)
	//}
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	btp.SetScale(ct_symbol, rlwe.NewScale(float64(params.Q()[0])/(float64(base))))
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)

	// Step 4 : Extend the basis from q to Q
	if ct_symbol, err = btp.ModUp(ct_symbol); err != nil {
		panic(err)
	}
	//PrintDebug(slice_length, *params, ct_symbol, value, cc.decryptor, cc.encoder)
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	//var real2, imag2 *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ct_symbol); err != nil {
		panic(err)
	}

	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	eval.Mul(real, 2.0, real)
	eval.Mul(imag, 2.0, imag)

	//eval.Add(real, imag, real)

	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)
	//fmt.Println("After CoeffsToSlot Scale : ", real.Scale.Uint64())

	// Step 6 : EvalMod (Homomorphic modular reduction)

	if imag, err = btp.EvalExpAndScale_with_vectroized_evaluation(imag, 2*math.Pi/complex(float64(base), 0.0), params); err != nil {
		panic(err)
	}

	if real, err = btp.EvalExpAndScale_with_vectroized_evaluation(real, 2*math.Pi/complex(float64(base), 0.0), params); err != nil {
		panic(err)
	}

	//fmt.Println("dddd")
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.PolynomialVector
	eval_poly_vec := make([]bignum.Polynomial, 2)

	//start := time.Now()
	eval_poly_vec[0] = bignum.NewPolynomial(0, HermiteInterpolation_symbol2symbol(), nil)
	eval_poly_vec[1] = bignum.NewPolynomial(0, make([]complex128, 4), nil)

	mapping := make(map[int][]int, 2)
	for i := 0; i < 2; i++ {
		mapping[i] = make([]int, params.MaxSlots()/2)
		for j := 0; j < params.MaxSlots()/2; j++ {
			mapping[i][j] = j + i*params.MaxSlots()/2
		}
	}

	if polys, err = polynomial.NewPolynomialVector(eval_poly_vec, mapping); err != nil {
		panic(err)
	}

	real, _ = polyEval.Evaluate(real, polys, params.DefaultScale())
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)

	eval_poly_vec[0] = bignum.NewPolynomial(0, HermiteInterpolation_symbol2symbol_imag(), nil)
	eval_poly_vec[1] = bignum.NewPolynomial(0, make([]complex128, 4), nil)

	if polys, err = polynomial.NewPolynomialVector(eval_poly_vec, mapping); err != nil {
		panic(err)
	}
	// [z]_b is computed
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)
	imag, _ = polyEval.Evaluate(imag, polys, params.DefaultScale())
	//PrintDebug(slice_length, *params, imag, value, cc.decryptor, cc.encoder)

	//value := make([]complex128, params.MaxSlots())
	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)

	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(*params, imag, value, cc.decryptor, cc.encoder)

	//fmt.Println("Scale : ", real.Scale.Uint64())
	//fmt.Println("Scale : ", imag.Scale.Uint64())

	// Recombines the real and imaginary part
	//if err = btp.Evaluator.Mul(imag, 1i, imag); err != nil {
	//	panic(err)
	//}

	if err = btp.Evaluator.Add(real, imag, ct_symbol); err != nil {
		panic(err)
	}

	//fmt.Println("Scale : ", ciphertext.Scale.Uint64())

	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(slice_length, *params, real, value, cc.decryptor, cc.encoder)
	//fmt.Print("Before Symbol bootstrapping")
	//PrintDebug(*params, ct_symbol, value, cc.decryptor, cc.encoder)

	return ct_symbol
}

func Cleaning(ciphertext *rlwe.Ciphertext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {
	var err error
	btp := cc.btp
	params := cc.params
	eval := cc.eval
	base := Bbox.base

	// Step 0: Some circuit in the slots domain
	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}
	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(base))))
	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}
	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	//fmt.Println("After CoeffsToSlot Scale : ", real.Scale.Uint64())

	// Step 6 : EvalMod (Homomorphic modular reduction)
	if imag, err = btp.EvalModAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	if real, err = btp.EvalElseAndScale(real, 2*math.Pi/complex(float64(base), 0.0)); err != nil {
		panic(err)
	}

	//fmt.Println("Scale : ", real.Scale.Uint64())
	//fmt.Println("Scale : ", imag.Scale.Uint64())

	// Recombines the real and imaginary part
	if err = btp.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if err = btp.Evaluator.Add(real, imag, ciphertext); err != nil {
		panic(err)
	}

	//fmt.Println("Scale : ", ciphertext.Scale.Uint64())

	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.Polynomial
	var eval_poly bignum.Polynomial

	//start := time.Now()
	eval_poly = bignum.NewPolynomial(0, HermiteInterpolation_exp2id(base, 2*base-1), nil)
	//t := time.Since(start)
	//fmt.Println(t)

	if polys = polynomial.NewPolynomial(eval_poly); err != nil {
		panic(err)
	}

	// [z]_b is computed
	ciphertext, _ = polyEval.Evaluate(ciphertext, polys, params.DefaultScale())

	return ciphertext
}

func Cleaning_with_vectorized_evaluation(ciphertext *rlwe.Ciphertext, Bbox Bk_Arithmetic_toolbox, cc Context) *rlwe.Ciphertext {
	var err error
	btp := cc.btp
	params := cc.params
	eval := cc.eval
	base := Bbox.base

	// Step 0: Some circuit in the slots domain
	if ciphertext.Level() > 4 {
		btp.DropLevel(ciphertext, ciphertext.Level()-4)
	}
	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = btp.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturally given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = btp.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	btp.SetScale(ciphertext, rlwe.NewScale(float64(params.Q()[0])/(float64(base))))
	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = btp.ModUp(ciphertext); err != nil {
		panic(err)
	}
	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = btp.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}
	eval.Conjugate(real, imag)
	eval.Add(real, imag, real)

	// Step 6 : EvalExp (Homomorphic modular reduction)
	//start := time.Now()
	if ciphertext, err = btp.EvalExpAndScale_with_vectroized_evaluation(real, 2*math.Pi/complex(float64(base), 0.0), params); err != nil {
		panic(err)
	}
	//elapsed := time.Since(start)
	//fmt.Println("EvalExp(cos) :", elapsed)

	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(*params, eval)

	var polys polynomial.PolynomialVector
	eval_poly_vec := make([]bignum.Polynomial, 2)

	//start := time.Now()
	eval_poly_vec[0] = bignum.NewPolynomial(0, HermiteInterpolation_exp2id(base, 2*base-1), nil)
	eval_poly_vec[1] = bignum.NewPolynomial(0, make([]complex128, 2*base), nil)

	mapping := make(map[int][]int, 2)
	for i := 0; i < 2; i++ {
		mapping[i] = make([]int, params.MaxSlots()/2)
		for j := 0; j < params.MaxSlots()/2; j++ {
			mapping[i][j] = j + i*params.MaxSlots()/2
		}
	}
	//t := time.Since(start)
	//fmt.Println(t)

	if polys, err = polynomial.NewPolynomialVector(eval_poly_vec, mapping); err != nil {
		panic(err)
	}

	ciphertext, _ = polyEval.Evaluate(ciphertext, polys, params.DefaultScale())

	return ciphertext
}
