(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageExport["AnomalousDimension"]
PackageExport["AnomalousDimTerm"]
PackageExport["BetaFunction"]
PackageExport["BetaTerm"]
PackageExport["CheckProjection"]
PackageExport["Finalize"]
PackageExport["QuarticBetaFunctions"]

PackageScope["AntiSym"]
PackageScope["ResetBetas"]
PackageScope["Sym"]
PackageScope["Sym4"]
PackageScope["Tdot"]
PackageScope["Ttimes"]
PackageScope["TsSym"]
PackageScope["TsSym4"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

AnomalousDimension::usage =
	"AnomalousDimTerm[field, loop] gives the anomalous dimension of the given field up to the l-loop order."
AnomalousDimTerm::usage =
	"AnomalousDimTerm[field, loop] gives the l-loop contribution to the anomalous dimension of the given field."
BetaFunction::usage =
	"BetaFunction[coupling, loop] computes the entire beta function of the coupling to the given loop order."
BetaTerm::usage =
	"BetaTerm[coupling, loop] computes the l-loop contribution to the beta function of the coupling."
CheckProjection::usage =
	"CheckProjection[coupling] returns the result of the automatic projection operator of the coupling on the corresponding generalized coupling."
Finalize::usage =
	"Finalize[expr] performs additional refinement of beta function expressions. Only works for expressions with up to 2 nontrivial flavor indices."
QuarticBetaFunctions::usage =
	"QuarticBetaFunctions[loop] returns all quartic beta functions to the given loop order using diagonalized projectors."

AntiSym::usage =
	"AntiSym[a, b][expr] is an internal function for antisymmetrizing expr in a and b."
ResetBetas::usage =
	"ResetBetas[] is a function used to dump all internally stored beta computations from the kernel."
Sym::usage =
	"Sym[a1, a2, a3, a4][expr] is an internal function for symmetrizing expr over up to four dummy indices."
Sym4::usage =
	"Sym4[a1, a2, a3, a4][expr] averages the expression over the 4 ways of switching a1 with one of the indices."
Tdot::usage =
	"Tdot[a, b,...] is an internal function the sequantially expands the matrix product of all the arguments."
Ttimes::usage =
	"Ttimes[a, b,...] is an internal function the sequantially expands the product of all the arguments."
TsSym::usage = TsSym4::usage =
	"TsSym and TsSym4 are functions used to symmetrize TStructures in their indices."

(*#####################################*)
(*----------Utility functions----------*)
(*#####################################*)
(*Symmetrizing tensor structures*)
SymmetrizeTS[TStructure[tsInds__][sa_], indices_, permutations_] :=
	Module[{elem, indPositions, ind, n, numPerm, posPermutations, rules, subs},
		numPerm = Length@ permutations;

		indPositions = Flatten@ Table[FirstPosition[{tsInds}, ind][[1]], {ind, indices}];
		posPermutations = Table[n, {n, Length@ {tsInds} }] /. MapThread[Rule, {indPositions, indPositions[[#]] }] & /@ permutations;
		subs = MapThread[Rule, {indices, indices[[#]] }] & /@ permutations;

		rules = Delete[ArrayRules@ sa, -1];
		rules = Table[#1[[Ordering@ posPermutations[[n]] ]] -> (#2 /. subs[[n]]), {n, numPerm}] & @@@ rules;
		TStructure[tsInds]@ SparseArray[GatherToList@ rules, Dimensions@ sa] /numPerm
	];

TsSym[i1_, i2_][expr_] := expr /. TStructure[inds__][sa_] :> SymmetrizeTS[TStructure[inds]@ sa, {i1, i2}, {{1, 2}, {2, 1}}];

TsSym[i1_, i2_, i3_][expr_] := expr /. TStructure[inds__][sa_] :> SymmetrizeTS[TStructure[inds]@ sa,
 	{i1, i2, i3}, Permutations@ {1, 2, 3}];

TsSym[i1_, i2_, i3_, i4_][expr_] := expr /. TStructure[inds__][sa_] :> SymmetrizeTS[TStructure[inds]@ sa,
 	{i1, i2, i3, i4}, Permutations@ {1, 2, 3, 4}];

TsSym4[i1_, i2_, i3_, i4_][expr_] := expr /. TStructure[inds__][sa_] :> SymmetrizeTS[TStructure[inds]@ sa,
	{i1, i2, i3, i4}, {{1, 2, 3, 4}, {2, 1, 3, 4}, {3, 2, 1, 4}, {4, 2, 3, 1}}];

(* Symmetrizing dummy index expressions: *)
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
AntiSym[indices__][expr_] :=
	Block[{indList, subs, s},
		indList = List[indices];
		subs = MapThread[Rule, {indList, #}] & /@ Permutations @ indList;
		Signature[indList] Sum[Signature[s[[;;, 2]] ] expr /.s, {s, subs}] / Factorial @Length @ indList
	];

CanonizeMatrices[expr_] := ReplaceAll[Matrix[y__][h1_[i1_], h2_[i2_]] /; !OrderedQ[{i1, i2}] :>
	Matrix[Sequence @@ Reverse[Trans /@ List@ y]][h1 @ i2, h2 @ i1] ] @ expr;

(*Functions that speed up evaluation by applying expand succesively to each couple of terms in the evaluation.*)
Ttimes[a_, b_, c___] := Ttimes[Expand[a b], c];
Ttimes[a_] = a;
Tdot[a_, b_, c___] := Tdot[Expand[a.b], c];
Tdot[a_] = a;


(*##################################*)
(*----------Beta functions----------*)
(*##################################*)

(*Function returns the l-loop contribution to the beta function of a given coupling*)
BetaTerm::loopNumber = "The `1` beta function has only been implemented up to `2` loops."
BetaTerm::unkown = "The coupling `1` has not been defined."
BetaTerm[coupling_, loop_Integer] :=
	Module[{beta, group, tensor, C1, C2},
		(*Determines the correct tensor structure based on coupling type*)
		Switch[$couplings @ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			If[loop > 3 || loop < 0,
				Message[BetaTerm::loopNumber, "gauge", 3];
				Abort[];
			];
			group = $couplings @ coupling;
			beta = GaugeTensors[coupling, loop] /. $gaugeCoefficients // Expand;

			(*Puts U1-mixing matrix on matrix form*)
			If[MatchQ[$gaugeGroups[group, LieGroup], U1[n_] /; n >1],
				beta = beta /. {Matrix[l1_List][group[adj] @ v1] Matrix[l2_List][group[adj] @ v2] :> Outer[Times, l1, l2],
					Matrix[mat_][group[adj] @ v1, group[adj] @ v2] /; MatrixQ @ mat :> mat};
			];
		,Yukawa,
			If[loop > 2 || loop < 0,
				Message[BetaTerm::loopNumber, "Yukawa", 2];
				Abort[];
			];
			beta = YukawaTensors[coupling, loop] /. $yukawaCoefficients // Expand;

		,Quartic,
			If[loop > 2 || loop < 0,
				Message[BetaTerm::loopNumber, "quartic", 2];
				Abort[];
			];
			beta = QuarticTensors[coupling, loop] /. $quarticCoefficients // Expand;

		,FermionMass,
			If[loop > 2 || loop < 0,
				Message[BetaTerm::loopNumber, "fermion mass", 2];
				Abort[];
			];
			beta = FermionMassTensors[coupling, loop] /. $yukawaCoefficients // Expand;

		,Trilinear,
			If[loop > 2 || loop < 0,
				Message[BetaTerm::loopNumber, "trilinear scalar", 2];
				Abort[];
			];
			beta = ScalarMassiveTensors[coupling, loop] /. $quarticCoefficients // Expand;

		,ScalarMass,
			If[loop > 2 || loop < 0,
				Message[BetaTerm::loopNumber, "scalar mass", 2];
				Abort[];
			];
			If[loop === 0,
				Return @ 0;
			];
			beta = ScalarMassiveTensors[coupling, loop] /. $quarticCoefficients // Expand;

		,_Missing,
			Message[BetaTerm::unkown, coupling];
			Abort[];
		];

		(*Canonically order coupling indices if any*)
		CanonizeMatrices @ RefineGroupStructures @ beta
	];

(*Function that produces the beta function for the requested coupling*)
Options[BetaFunction] = {RescaledCouplings -> False, FourDimensions -> True};
BetaFunction::unkown = "The coupling `1` has not been defined."
BetaFunction::loopNumber = "The `1` beta function has only been implemented up to `2` loops."
BetaFunction[coupling_Symbol, loop_Integer, OptionsPattern[] ] ? OptionsCheck :=
	Module[{coef = 4 Pi, firstTerm = 0, l},
		Switch[$couplings @ coupling
		,gr_ /; MemberQ[Keys @ $gaugeGroups, gr],
			If[loop > 3 || loop < 0,
				Message[BetaFunction::loopNumber, "gauge", 3];
				Abort[];
			];
		,Yukawa|Quartic|ScalarMass|Trilinear|FermionMass,
			If[loop > 2 || loop < 0,
				Message[BetaFunction::loopNumber, $couplings @ coupling, 2];
				Abort[];
			];
		,_,
			Message[BetaFunction::unkown, coupling];
			Abort[];
		];

		If[OptionValue @ RescaledCouplings, coef = 1; ];
		If[OptionValue @ FourDimensions, firstTerm = 1; ];

		Sum[ Power[coef, -2 l] BetaTerm[coupling, l], {l, firstTerm, loop}]
	];

(*Fuction for diagonalizing the quartic beta functions. It inherits the options from BetaFunction*)
QuarticBetaFunctions::singular = "The projection matrix is singular. Some of the couplings may be redundant."
QuarticBetaFunctions::loopNumber = "The quartic beta function has only been implemented up to 2 loops."
QuarticBetaFunctions[loop_Integer, opt:OptionsPattern[]] ? OptionsCheck :=
	Module[{betaFunctions, couplings, qProjections, invMatrix, c},
		If[loop > 2,
			Message[QuarticBetaFunctions::loopNumber];
			Abort[];
		];

		couplings = Keys @ $quartics;
		Print["The quartic couplings are ", couplings];

		(*Finds inversion matrix for the quartic projectors*)
		qProjections = CheckProjection /@ couplings;
		invMatrix = Transpose @ Table[Simplify @ D[qProjections, c], {c, couplings}];
		If[Det @ invMatrix === 0,
			Message[QuarticBetaFunctions::singular];
			Abort[];
		];
		invMatrix = Inverse @ invMatrix;

		(*Extracts beta functions*)
		betaFunctions = Monitor[
							Table[BetaFunction[c, loop, opt], {c, couplings}]
						,StringForm["Evaluating the `` \[Beta]-function", c] ];

		invMatrix . betaFunctions // Expand
	];

(*Function returns the l-loop contribution to the beta function of a given coupling*)
AnomalousDimTerm::loopNumber = "The `1` scalar anomalous dimension has only been implemented up to `2` loops."
AnomalousDimTerm::unkown = "The field `1` has not been defined."
AnomalousDimTerm[field_, loop_Integer] :=
	Module[{anomDim},
		(*Determines the correct tensor structure based on the field type*)
		Switch[field
		,x_ /; MemberQ[Keys @ $fermions, x],
			If[loop > 2 || loop < 1,
				Message[AnomalousDimTerm::loopNumber, "fermion", 2];
				Abort[];
			];
			anomDim = FermionAnomalousTensors[field, loop] /. $fermionAnomalousCoefficients;

		,x_ /; MemberQ[Keys @ $scalars, x],
			If[loop > 2 || loop < 1,
				Message[AnomalousDimTerm::loopNumber, "scalar", 2];
				Abort[];
			];
			anomDim = ScalarAnomalousTensors[field, loop] /. $scalarAnomalousCoefficients;

		,_,
			Message[AnomalousDimTerm::unkown, field];
			Abort[];
		];

		(*Canonically order coupling indices if any*)
		CanonizeMatrices @ RefineGroupStructures @ anomDim
	];

(*Function that produces the anomalous dimension for a given field*)
Options[AnomalousDimension] = {RescaledCouplings -> False};
AnomalousDimension::unkown = "The coupling `1` has not been defined."
AnomalousDimension::loopNumber = "The anomalous dimension has only been implemented up to `1` loops."
AnomalousDimension[field_, loop_Integer, OptionsPattern[] ] ? OptionsCheck :=
	Module[{coef = 4 Pi, l},
		(* Checks if the input request has been implemented.*)
		Switch[field
		,x_ /; MemberQ[Keys @ $fermions, x] || MemberQ[Keys @ $scalars, x],
			If[loop > 2 || loop < 1,
				Message[AnomalousDimension::loopNumber, 2];
				Abort[];
			];
		,_,
			Message[AnomalousDimension::unkown, field];
			Abort[];
		];

		If[OptionValue @ RescaledCouplings, coef = 1;];

		Sum[ Power[coef, -2 l] AnomalousDimTerm[field, l], {l, loop}]
	];

(*###################################*)
(*----------Extra functions----------*)
(*###################################*)

(*Function for finalizing a betafunction, bringing it from a nice compact output to a form more suitable for
further Mathematica manipulations. Can also be used to specify particular cases for coupling matrices.*)
Options[Finalize] = {Parameterizations -> {}, BarToConjugate -> False};
Finalize[expr_, opt:OptionsPattern[]] ? OptionsCheck :=
	Internal`InheritedBlock[{out, Bar, Trans, Matrix},
		out = CanonizeMatrices @ expr /. OptionValue @ Parameterizations;
			Bar[a_List] := Bar /@ a;
			Bar[Times[a_, b__]] := Bar /@ Times[a, b];
			Bar[Plus[a_, b__]] := Bar /@ Plus[a, b];
			If[OptionValue @ BarToConjugate,
				Bar[a_Symbol] := Conjugate @ a;
			];
			Bar[a_] /; NumberQ[a] := Conjugate @ a;
			Trans[a_List] /; MatrixQ[a] := Transpose @ a;
			Trans[a_List] /; VectorQ[a] := a;
		Matrix[y__][__] := Dot[y];
		out/. {del[_, $i, $j]-> 1}
	];

(*Function to check the projected value of a coupling*)
CheckProjection::unkown = "The coupling `1` has not been defined."
CheckProjection[coupling_Symbol] :=
	Module[{tensor, A, B, a, i, j, b, c, d},
		Switch[$couplings @ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			G2Matrix[$A, $B] $gaugeGroups[$couplings @ coupling, Projector][$B, $A]
		,Yukawa,
			Tr[Yuk[$a, $i, $j]. $yukawas[coupling, Projector][$a, $i, $j] ]
		,Quartic,
			Lam[$a, $b, $c, $d] $quartics[coupling, Projector][$a, $b, $c, $d]
		,FermionMass,
			Tr[Yuk[$a, $i, $j, True]. $fermionMasses[coupling, Projector][$a, $i, $j] ]
		,Trilinear,
			Lam[$a, $b, $c, $d, True] $trilinears[coupling, Projector][$a, $b, $c, $d]
		,ScalarMass,
			Lam[$a, $b, $c, $d, True] $scalarMasses[coupling, Projector][$a, $b, $c, $d]
		,_Missing,
			Message[CheckProjection::unkown, coupling];
			$Failed
		]// RefineGroupStructures// CanonizeMatrices// Simplify
	];
