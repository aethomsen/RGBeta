Begin["Betas`"]

(*#####################################*)
(*----------Utility functions----------*)
(*#####################################*)
(*Symmetrizing in indices*)
Sym[i1_, i2_][expr_] := expr /2 + ReplaceAll[expr, {i1 -> i2, i2 -> i1}] /2 ;
AntiSym[i1_, i2_][expr_] := expr /2 - ReplaceAll[expr, {i1 -> i2, i2 -> i1}] /2 ;
Sym[a_, b_, c_, d_][expr_] := 
	Block[{perm},
		perm = {a -> #[[1]], b -> #[[2]], c -> #[[3]], d -> #[[4]]} & /@ Permutations[{a, b, c, d}];
		Mean[expr/.perm]	
	];
Sym4[a_, b_, c_, d_][expr_] := 
	Block[{perm},
		perm = {a -> #, # -> a} & /@ {a, b, c, d};
		Mean[expr/.perm]	
	];	
Sym[a_, b_, c_][expr_] := 
	Block[{perm},
		perm = {a -> #[[1]], b -> #[[2]], c -> #[[3]]} & /@ Permutations[{a, b, c}];
		Mean[expr/.perm]	
	];
	

(*Functions that speed up evaluation by applying expand succesively to each couple of terms in the evaluation.*)
Ttimes[a_, b_, c___] := Ttimes[Expand[a b], c];
Ttimes[a_] = a;
Tdot[a_, b_, c___] := Tdot[Expand[a.b], c];
Tdot[a_] = a;


(*##################################*)
(*----------Beta functions----------*)
(*##################################*)

(*Function returns the beta function of the given coupling and loop order*)
BetaTerm::gaugeLoops = "The gauge beta function is only implemented to 3 loops."
BetaTerm::yukawaLoops = "The Yukawa beta function is only implemented to 2 loops."
BetaTerm::quarticLoops = "The quartic beta function is only implemented to 2 loops."
BetaTerm::loopNumber = "The `1` beta function is only implemented up to `2` loops."
BetaTerm::unkown = "The coupling `1` has not been defined."
BetaTerm[coupling_, loop_Integer] :=
	Module[{beta, group, tensor, C1, C2},
		(*Determines the correct tensor structure based on coupling type*)
		Switch[$couplings @ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			If[loop > 3, 
				Message[BetaTerm::loopNumber, "gauge", 3];
				Return[Null];
			];
			group = $couplings @ coupling;
			beta = Ttimes[G2Matrix[C1, $A] ,GaugeTensors[loop], G2Matrix[$B, C2], $gaugeGroups[group, Projector][C1, C2] ];
			
			(*Puts U1-mixing matrix on matrix form*)
			If[MatchQ[$gaugeGroups[group, LieGroup], U1[n_] /; n >1],
				beta = beta /. Matrix[l1_List][group[adj] @ v1] Matrix[l2_List][group[adj] @ v2] :> Outer[Times, l1, l2];
			];
		,Yukawa,
			If[loop > 2, 
				Message[BetaTerm::loopNumber, "Yukawa", 2];
				Return[Null];
			];
			
			Switch[$yukawas[coupling, Chirality]
			,Left,
				tensor = YukawaTensors[loop][[1, 1]];
			,Right,
				tensor = YukawaTensors[loop][[2, 2]];
			];
			
			beta = tensor $yukawas[coupling, Projector][$a, $i, $j] // Expand // Expand;
		,Quartic,
			If[loop > 2, 
				Message[BetaTerm::loopNumber, "quartic", 2];
				Return[Null];
			];
			
			beta = QuarticTensors[loop] $quartics[coupling, Projector][$a, $b, $c, $d] // Expand // Expand;
		,FermionMass,
			If[loop > 2, 
				Message[BetaTerm::loopNumber, "fermion mass", 2];
				Return[Null];
			];
			
			Switch[$fermionMasses[coupling, Chirality]
			,Left,
				tensor = YukawaTensors[loop][[1, 1]];
			,Right,
				tensor = YukawaTensors[loop][[2, 2]];
			];
			
			beta = tensor $fermionMasses[coupling, Projector][$a, $i, $j] // Expand // Expand;
		,Trilinear,
			If[loop > 2, 
				Message[BetaTerm::loopNumber, "trilinear scalar", 2];
				Return[Null];
			];
			
			beta = QuarticTensors[loop] $trilinears[coupling, Projector][$a, $b, $c, $d] // Expand // Expand;
		,ScalarMass,
			If[loop > 2, 
				Message[BetaTerm::loopNumber, "scalar mass", 2];
				Return[Null];
			];
			
			beta = QuarticTensors[loop] $scalarMasses[coupling, Projector][$a, $b, $c, $d] // Expand // Expand;
		,_Missing,
			Message[BetaTerm::unkown, coupling];
			Return[Null];
		];
		
		Return @ beta;
	];

(*Function that produces the beta function for the requested coupling*)
BetaFunction::unkown = "The coupling `1` has not been defined."
BetaFunction[coupling_Symbol, loop_Integer, OptionsPattern[{RescaledCouplings -> False, FourDimensions -> True}] ] :=
	Block[{coef = 4 Pi, firstTerm = 0, l},
		If[Head @ $couplings @ coupling === Missing, 
			Message[BetaFunction::unkown, coupling];
			Return @ Null;
		];
		
		If[OptionValue @ RescaledCouplings, coef = 1; ];
		If[OptionValue @ FourDimensions, firstTerm = 1; ];
		
		Sum[ Power[coef, -2 l] BetaTerm[coupling, l], {l, firstTerm, loop}]
	];

(*Fuction for diagonalizing the quartic beta functions. It inherits the options from BetaFunction*)
QuarticBetaFunctions[loop_Integer, opt:OptionsPattern[]] :=
	Block[{betaFunctions, couplings, qProjections, invMatrix},
		couplings = Keys @ $quartics;
		Print["The quartic couplings are ", couplings];
		
		(*Finds inversion matrix for the quartic projectors*)
		qProjections = CheckProjection /@ couplings;
		invMatrix = Inverse @ Transpose @ Table[Simplify @ D[qProjections, c], {c, couplings}];
		
		(*Extracts beta functions*)
		betaFunctions = Monitor[
							Table[BetaFunction[c, loop, opt], {c, couplings}]
						,StringForm["Evaluating the `` \[Beta]-function", c] ];
		
		Return[invMatrix . betaFunctions];
	];


(*Function for finalizing a betafunction, bringing it from a nice compact output to a form more suitable for 
further Mathematica manipulations. Can also be used to specify particular cases for coupling matrices.*)
Finalize[expr_, OptionsPattern[{Parametrizations -> {}}] ] :=
	Internal`InheritedBlock[{out, Bar, Trans, Matrix},
		out = expr /. OptionValue @ Parametrizations;
			Bar[a_List] := Bar /@ a;
			Bar[Times[a__]] := Bar /@ Times[a];
			Bar[Plus[a__]] := Bar /@ Plus[a];
			Bar[a_Symbol] := Conjugate @ a;
			Bar[a_] /; NumberQ[a] := Conjugate @ a;
			Trans[a_List] /; MatrixQ[a] === True := Transpose @ a;
		Matrix[y__][a_[f1_], b_[f2_]] /; !OrderedQ[{f1, f2}] :=
			Matrix[Sequence @@ Reverse[Trans /@ List@ y]][b[f2], a[f1]];
		Matrix[y__][a_[f1], b_[f2]] := Dot[y];
		(*Matrix /: Matrix[u_List][_[v1]] Matrix[w_List][_[v2]] := List /@ u . List @ w;*)
		out
	];

(*Function to check the projected value of a coupling*)
CheckProjection::unkown = "The coupling `1` has not been defined."
CheckProjection[coupling_Symbol] :=
	Module[{cop, tensor, A, B, a, i, j, b, c, d},
		Switch[$couplings @ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			cop = Ttimes[G2Matrix[A, B], $gaugeGroups[$couplings @ coupling, Projector][B, A] ] // Expand;
		,Yukawa,			
			Switch[$yukawas[coupling, Chirality]
			,Left,
				tensor = Yuk[a, i, j][[1, 1]];
			,Right,
				tensor = Yuk[a, i, j][[2, 2]];
			];
			cop = tensor $yukawas[coupling, Projector][a, i, j] // Expand // Expand;
		,Quartic,
			cop = Lam[a, b, c, d] $quartics[coupling, Projector][a, b, c, d] // Expand // Expand;
		,FermionMass,
			Switch[$fermionMasses[coupling, Chirality]
			,Left,
				tensor = Yuk[a, i, j][[1, 1]];
			,Right,
				tensor = Yuk[a, i, j][[2, 2]];
			];
			cop = tensor $fermionMasses[coupling, Projector][a, i, j] // Expand // Expand;
		,Trilinear,
			cop = Lam[a, b, c, d] $trilinears[coupling, Projector][a, b, c, d] // Expand // Expand;
		,ScalarMass,
			cop = Lam[a, b, c, d] $scalarMasses[coupling, Projector][a, b, c, d] // Expand // Expand;
		,_Missing,
			Message[CheckProjection::unkown, coupling];
			Return @ Null;
		];
		Return @ cop;
	];



(*#####################################*)
(*----------Beta Coefficients----------*)
(*#####################################*)
$quarticCoefficients = {
	(* 1-loop *)
	Bcoef[3, 1, 1] -> 36,
	Bcoef[3, 1, 2] -> -12,
	Bcoef[3, 1, 3] -> 3,
	Bcoef[3, 1, 4] -> 2,
	Bcoef[3, 1, 5] -> -12,
	(*2-loop*)
	Bcoef[3, 2, 1] -> 324,
	Bcoef[3, 2, 2] -> -684,
	Bcoef[3, 2, 3] -> 646,
	Bcoef[3, 2, 4] -> -28,
	Bcoef[3, 2, 5] -> -32,
	Bcoef[3, 2, 6] -> 12,
	Bcoef[3, 2, 7] -> 60,
	Bcoef[3, 2, 8] -> 0,
	Bcoef[3, 2, 9] -> 6,
	Bcoef[3, 2, 10] -> -143/3,
	Bcoef[3, 2, 11] -> 11/3,
	Bcoef[3, 2, 12] -> 10/3,
	Bcoef[3, 2, 13] -> -18,
	Bcoef[3, 2, 14] -> 24,
	Bcoef[3, 2, 15] -> -18,
	Bcoef[3, 2, 16] -> 1/3,
	Bcoef[3, 2, 17] -> -6,
	Bcoef[3, 2, 18] -> 0,
	Bcoef[3, 2, 19] -> -144,
	Bcoef[3, 2, 20] -> 60,
	Bcoef[3, 2, 21] -> 10,
	Bcoef[3, 2, 22] -> 0,
	Bcoef[3, 2, 23] -> -3,
	Bcoef[3, 2, 24] -> 0,
	Bcoef[3, 2, 25] -> 24,
	Bcoef[3, 2, 26] -> -48,
	Bcoef[3, 2, 27] -> 12,
	Bcoef[3, 2, 28] -> 0,
	Bcoef[3, 2, 29] -> -2,
	Bcoef[3, 2, 30] -> -3,
	Bcoef[3, 2, 31] -> 48,
	Bcoef[3, 2, 32] -> 24,
	Bcoef[3, 2, 33] -> 24
}; 

$yukawaCoefficients = {
	(* 1-loop *)
	Bcoef[2, 1, 1] -> 0,
	Bcoef[2, 1, 2] -> -6,
	Bcoef[2, 1, 3] -> 2,
	Bcoef[2, 1, 4] -> 1,
	Bcoef[2, 1, 5] -> 1/2,
	(* 2-loop *)
	Bcoef[2, 2, 1] -> -21/2,
	Bcoef[2, 2, 2] -> 12, 
	Bcoef[2, 2, 3] -> 0,
	Bcoef[2, 2, 4] -> -3,
	Bcoef[2, 2, 5] -> 49/4,
	Bcoef[2, 2, 6] -> -1/4,
	Bcoef[2, 2, 7] -> -1/2,
	Bcoef[2, 2, 8] -> -97/3,
	Bcoef[2, 2, 9] -> 11/6,
	Bcoef[2, 2, 10] -> 5/3,
	Bcoef[2, 2, 11] -> 1/12,
	Bcoef[2, 2, 12] -> 12,
	Bcoef[2, 2, 13] -> 0,
	Bcoef[2, 2, 14] -> 6,
	Bcoef[2, 2, 15] -> -12,
	Bcoef[2, 2, 16] -> 10,
	Bcoef[2, 2, 17] -> 6,
	Bcoef[2, 2, 18] -> 5/2,
	Bcoef[2, 2, 19] -> 9,
	Bcoef[2, 2, 20] -> -1/2,
	Bcoef[2, 2, 21] -> -7/2,
	Bcoef[2, 2, 22] -> 0,
	Bcoef[2, 2, 23] -> -2,
	Bcoef[2, 2, 24] -> 2,
	Bcoef[2, 2, 25] -> 0,
	Bcoef[2, 2, 26] -> -2,
	Bcoef[2, 2, 27] -> 0,
	Bcoef[2, 2, 28] -> -1/2,
	Bcoef[2, 2, 29] -> -2,
	Bcoef[2, 2, 30] -> -1/4,
	Bcoef[2, 2, 31] -> -3/4,
	Bcoef[2, 2, 32] -> -1,
	Bcoef[2, 2, 33] -> -3/4
};

$gaugeCoefficients = {
	(* 1-loop *)
	Bcoef[1, 1, 1] -> -22/3,
	Bcoef[1, 1, 2] -> 2/3,
	Bcoef[1, 1, 3] -> 1/3,
	(* 2-loop *)
	Bcoef[1, 2, 1] -> 2,
	Bcoef[1, 2, 2] -> 4,
	Bcoef[1, 2, 3] -> -68/3,
	Bcoef[1, 2, 4] -> 10/3,
	Bcoef[1, 2, 5] -> 2/3,
	Bcoef[1, 2, 6] -> -1,
	Bcoef[1, 2, 7] -> 0,
	(* 3-loop *)
	Bcoef[1, 3, 1] -> -1,
	Bcoef[1, 3, 2] -> 29/2,
	Bcoef[1, 3, 3] -> 133/18,
	Bcoef[1, 3, 4] -> 679/36,
	Bcoef[1, 3, 5] -> -11/18,
	Bcoef[1, 3, 6] -> -25/18,
	Bcoef[1, 3, 7] -> -23/36,
	Bcoef[1, 3, 8] -> -49/36,
	Bcoef[1, 3, 9] -> 4,
	Bcoef[1, 3, 10] -> 25/2,
	Bcoef[1, 3, 11] -> -2857/27,
	Bcoef[1, 3, 12] -> -79/108,
	Bcoef[1, 3, 13] -> 1/54,
	Bcoef[1, 3, 14] -> 1415/54,
	Bcoef[1, 3, 15] -> 545/108,
	Bcoef[1, 3, 16] -> -29/54,
	Bcoef[1, 3, 17] -> 1,
	Bcoef[1, 3, 18] -> -1/12,
	Bcoef[1, 3, 19] -> -5/4,
	Bcoef[1, 3, 20] -> -1/4,
	Bcoef[1, 3, 21] -> -1,
	Bcoef[1, 3, 22] -> -7,
	Bcoef[1, 3, 23] -> -7/2,
	Bcoef[1, 3, 24] -> -6,
	Bcoef[1, 3, 25] -> 9/4,
	Bcoef[1, 3, 26] -> 1,
	Bcoef[1, 3, 27] -> -1,
	Bcoef[1, 3, 28] -> 3/2,
	Bcoef[1, 3, 29] -> 7/8,
	Bcoef[1, 3, 30] -> 1/2,
	Bcoef[1, 3, 31] -> 1/8,
	Bcoef[1, 3, 32] -> 3/8,
	Bcoef[1, 3, 33] -> -1/8
};

Protect[$gaugeCoefficients, $quarticCoefficients, $yukawaCoefficients];

End[]