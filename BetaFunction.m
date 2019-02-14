<< GroupsAndIndices`
<< FieldsAndCouplings`
<< Tensors`

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
(*Global couplings variable*)
$couplings = <||>;

(*Function returns the beta function of the given coupling and loop order*)
BetaFunction::gaugeLoops = "The gauge beta function is only implemented to 3 loops."
BetaFunction::yukawaLoops = "The Yukawa beta function is only implemented to 2 loops."
BetaFunction::quarticLoops = "The quartic beta function is only implemented to 1 loops."
BetaFunction::unkown = "The coupling `1` has not been defined."
BetaFunction[coupling_Symbol, loop_Integer] := 
	Block[{beta},
		Switch[$couplings @ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			If[loop > 3, 
				Message[BetaFunction::gaugeLoops];
				Return[Null];
			];
			
			beta = Ttimes[GaugeTensors[loop], G2[B, C, 3/2], $gaugeGroups[$couplings @ coupling, Projector][C, A] ] /. gaugeCoefficients;
		,Yukawa,
			If[loop > 2, 
				Message[BetaFunction::yukawaLoops];
				Return[Null];
			];
			
			Switch[$yukawas[coupling, Chirality]
			,Left,
				tensor = YukawaTensors[loop][[1, 1]];
			,Right,
				tensor = YukawaTensors[loop][[2, 2]];
			];
			
			beta = tensor $yukawas[coupling, Projector][a, i, j] /. yukawaCoefficients // Expand;
		,Quartic,
			If[loop > 1, 
				Message[BetaFunction::quarticLoops];
				Return[Null];
			];
			
			beta = QuarticTensors[loop] $quartics[coupling, Projector][a, b, c, d] /. quarticCoefficients // Expand;
		,_,
			Message[BetaFunction::unkown, coupling];
			Return[Null];
		];
		
		Return @ beta;
	];

(*Function for finalizing a betafunction, bringing it from a nice compact output to a form more suitable for 
further Mathematica manipulations. Can also be used to specify particular cases for coupling matrices.*)
Finalize[expr_, OptionsPattern[{Parametrizations -> {} }] ] :=
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
		out
	];


(*#####################################*)
(*----------Beta Coefficients----------*)
(*#####################################*)

quarticCoefficients = {
	(* 1-loop *)
	B[3, 1, 1] ->  36,
	B[3, 1, 2] -> -12,
	B[3, 1, 3] -> 3,
	B[3, 1, 4] -> 2,
	B[3, 1, 5] -> -12
}; 

yukawaCoefficients = {
	(* 1-loop *)
	B[2, 1, 1] -> 0,
	B[2, 1, 2] -> -6,
	B[2, 1, 3] -> 2,
	B[2, 1, 4] -> 1,
	B[2, 1, 5] -> 1/2,
	(* 2-loop *)
	B[2, 2, 1] -> -21/2,
	B[2, 2, 2] -> 12, 
	B[2, 2, 3] -> 0,
	B[2, 2, 4] -> -3,
	B[2, 2, 5] -> 49/4,
	B[2, 2, 6] -> -1/4,
	B[2, 2, 7] -> -1/2,
	B[2, 2, 8] -> -97/3,
	B[2, 2, 9] -> 11/6,
	B[2, 2, 10] -> 5/3,
	B[2, 2, 11] -> 1/12,
	B[2, 2, 12] -> 12,
	B[2, 2, 13] -> 0,
	B[2, 2, 14] -> 6,
	B[2, 2, 15] -> -12,
	B[2, 2, 16] -> 10,
	B[2, 2, 17] -> 6,
	B[2, 2, 18] -> 5/2,
	B[2, 2, 19] -> 9,
	B[2, 2, 20] -> -1/2,
	B[2, 2, 21] -> -7/2,
	B[2, 2, 22] -> 0,
	B[2, 2, 23] -> -2,
	B[2, 2, 24] -> 2,
	B[2, 2, 25] -> 0,
	B[2, 2, 26] -> -2,
	B[2, 2, 27] -> 0,
	B[2, 2, 28] -> -1/2,
	B[2, 2, 29] -> -2,
	B[2, 2, 30] -> -1/4,
	B[2, 2, 31] -> -3/4,
	B[2, 2, 32] -> -1,
	B[2, 2, 33] -> -3/4
};

gaugeCoefficients = {
	(* 1-loop *)
	B[1, 1, 1] -> -11/3,
	B[1, 1, 2] -> 1/3,
	B[1, 1, 3] -> 1/6,
	(* 2-loop *)
	B[1, 2, 1] -> 1,
	B[1, 2, 2] -> 2,
	B[1, 2, 3] -> -34/3,
	B[1, 2, 4] -> 5/3,
	B[1, 2, 5] -> 1/3,
	B[1, 2, 6] -> -1/2,
	B[1, 2, 7] -> 0,
	(* 3-loop *)
	B[1, 3, 1] -> -1/2,
	B[1, 3, 2] -> 29/4,
	B[1, 3, 3] -> 133/36,
	B[1, 3, 4] -> 679/72,
	B[1, 3, 5] -> -11/36,
	B[1, 3, 6] -> -25/36,
	B[1, 3, 7] -> -23/72,
	B[1, 3, 8] -> -49/72,
	B[1, 3, 9] -> 2,
	B[1, 3, 10] -> 25/4,
	B[1, 3, 11] -> -2857/54,
	B[1, 3, 12] -> -79/216,
	B[1, 3, 13] -> 1/108,
	B[1, 3, 14] -> 1415/108,
	B[1, 3, 15] -> 545/216,
	B[1, 3, 16] -> -29/108,
	B[1, 3, 17] -> 1/2,
	B[1, 3, 18] -> -1/24,
	B[1, 3, 19] -> -5/8,
	B[1, 3, 20] -> -1/8,
	B[1, 3, 21] -> -1/2,
	B[1, 3, 22] -> -7/2,
	B[1, 3, 23] -> -7/4,
	B[1, 3, 24] -> -3,
	B[1, 3, 25] -> 9/8,
	B[1, 3, 26] -> 1/2,
	B[1, 3, 27] -> -1/2,
	B[1, 3, 28] -> 3/4,
	B[1, 3, 29] -> 7/16,
	B[1, 3, 30] -> 1/4,
	B[1, 3, 31] -> 1/16,
	B[1, 3, 32] -> 3/16,
	B[1, 3, 33] -> -1/16
};