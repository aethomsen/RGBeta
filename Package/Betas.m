(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageExport["Gauge"]
PackageExport["Yukawa"]
PackageExport["Quartic"]
PackageExport["Trilinear"]
PackageExport["ScalarMass"]
PackageExport["FermionMass"]

PackageExport["AnomalousDimension"]
PackageExport["AnomalousDimTerm"]
PackageExport["BetaFunction"]
PackageExport["BetaTerm"]
PackageExport["CheckProjection"]
PackageExport["Finalize"]
PackageExport["QuarticBetaFunctions"]
PackageExport["CheckQuarticMixing"]
PackageExport["UpsilonFunction"]
PackageExport["UpsilonTerm"]

PackageScope["AntiSym"]
PackageScope["Sym"]
PackageScope["Sym4"]
PackageScope["SymmetrizeTS"]
PackageScope["Tdot"]
PackageScope["Ttimes"]
PackageScope["TsSym"]
PackageScope["TsSym4"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

Gauge::usage= Yukawa::usage= Quartic::usage= Trilinear::usage= ScalarMass::usage= FermionMass::usage=
	"Is a symbol used to denote a class of couplings in BetaTerm or BetaFunction."

AnomalousDimension::usage =
	"AnomalousDimension[field, loop] gives the anomalous dimension of the given field(s) up to the l-loop order. Use {f1, f2} for field for mixing term between multiple fields of similar quantum numbers."
AnomalousDimTerm::usage =
	"AnomalousDimTerm[field, loop] gives the l-loop contribution to the anomalous dimension of the given field. Use {f1, f2} for field for mixing term between multiple fields of similar quantum numbers."
BetaFunction::usage =
	"BetaFunction[coupling, loop] computes the entire beta function of the coupling to the given loop order.
	Option FlavorImproved-> True (default) uses the flavor-improved \[Beta]-function (relevant from 3-loop order).
	Option FourDimensions-> True (default) removes the 0-order \[Epsilon]-dependent part.
	Option RescaledCouplings-> True (not default) will absorb the loop factors of 1/16\[Pi]^2 into the couplings.
	Rather than a single coupling, Gauge, Yukawa, Quartic, Trilinear, ScalarMass, and FermionMass can be used for a class of couplings."
BetaTerm::usage =
	"BetaTerm[coupling, loop] computes the l-loop contribution to the beta function of the coupling.
	Option FlavorImproved-> True (default) uses the flavor-improved \[Beta]-function (relevant from 3-loop order).
	Rather than a single coupling, Gauge, Yukawa, Quartic, Trilinear, ScalarMass, and FermionMass can be used for a class of couplings."
CheckProjection::usage =
	"CheckProjection[coupling] returns the result of the automatic projection operator of the coupling on the corresponding generalized coupling."
Finalize::usage =
	"Finalize[expr] performs additional refinement of beta function expressions. Only works for expressions with up to 2 nontrivial flavor indices."
QuarticBetaFunctions::usage =
	"QuarticBetaFunctions[loop] returns all quartic beta functions to the given loop order using diagonalized projectors."
UpsilonFunction::usage =
	"UpsilonFunction[field, loop] gives the anomalous dimension of the given field up to the l-loop order. Use {f1, f2} for field for mixing term between multiple fields of similar quantum numbers."
UpsilonTerm::usage =
	"UpsilonTerm[field, loop] gives the l-loop contribution to the anomalous dimension of the given field. Use {f1, f2} for field for mixing term between multiple fields of similar quantum numbers."

AntiSym::usage =
	"AntiSym[a, b][expr] is an internal function for antisymmetrizing expr in a and b."
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

$betaTypes= {Gauge, Yukawa, Quartic, Trilinear, ScalarMass, FermionMass};
$betaOrders= <|Gauge-> 4, Yukawa-> 3, Quartic-> 3 , Trilinear-> 2, ScalarMass-> 2, FermionMass-> 2|>;

CouplingType::unkown = "The coupling `1` has not been defined."
CouplingType[coupling_] := Switch[$couplings@ coupling,
	_Missing,
		Message[CouplingType::unkown, coupling];
		Abort[],
	x_ /; MemberQ[Keys @ $gaugeGroups, x], Gauge,
	_, $couplings@ coupling
]

CheckKnownBeta::loopNumber = "The `1` beta function has been implemented only up to `2` loops."
CheckKnownBeta::loopint= "Loop number `1` should be a non-negative integer."
CheckKnownBeta[type_, loop_] :=
	Module[{},
		If[!IntegerQ@ loop || Negative@ loop,
			Message[CheckKnownBeta::loopint, loop];
			Abort[];
		];
		betaType= If[MemberQ[$betaTypes, type],
				type,
				CouplingType@ type
			];
		If[loop > $betaOrders@ betaType,
			Message[CheckKnownBeta::loopNumber, betaType, $betaOrders@ betaType];
			Abort[];
		]
	]

(*Function returns the l-loop contribution to the beta function of a given coupling*)
Options[BetaTerm] = {
		FlavorImproved-> True (*Determines if returning the flavor-improved beta function*)
	};
BetaTerm[coupling_, loop_, opt:OptionsPattern[] ]? OptionsCheck :=
	Module[{beta, group, tensor, C1, C2},
		CheckKnownBeta[coupling, loop];

		If[MemberQ[$betaTypes, coupling], Return@ ProjectionToUnmixedBetas[coupling, loop, BetaTerm, opt];];

		CheckPotentialCouplingMixing[coupling, BetaTerm];
		(*Determines the correct tensor structure based on coupling type*)
		Switch[CouplingType @ coupling
		,Gauge,
			group = $couplings @ coupling;
			beta = GaugeTensors[coupling, loop] /. $gaugeCoefficients // Expand;

			(*Puts U1-mixing matrix on matrix form*)
			If[MatchQ[$gaugeGroups[group, LieGroup], U1[n_] /; n >1],
				beta = beta /. {Matrix[l1_List][group[adj] @ v1] Matrix[l2_List][group[adj] @ v2] :> Outer[Times, l1, l2],
					Matrix[mat_][group[adj] @ v1, group[adj] @ v2] /; MatrixQ @ mat :> mat};
			];
		,Yukawa,
			beta = YukawaTensors[coupling, loop] + If[OptionValue@ FlavorImproved && loop > 0,
					UpsilonYukawaTensors[coupling, loop]/. $fermionUpsilonCoefficients/. $scalarUpsilonCoefficients, 0
				]/. $yukawaCoefficients // Expand;
		,Quartic,
			CheckQuarticMixing[coupling, BetaTerm];
			beta = QuarticTensors[coupling, loop] /. $quarticCoefficients // Expand;

		,FermionMass,
			beta = FermionMassTensors[coupling, loop] /. $yukawaCoefficients // Expand;

		,Trilinear,
			beta = ScalarMassiveTensors[coupling, loop] /. $quarticCoefficients // Expand;

		,ScalarMass,
			If[loop === 0, Return @ 0; ];
			beta = ScalarMassiveTensors[coupling, loop] /. $quarticCoefficients // Expand;
		];

		(*Canonically order coupling indices if any*)
		CanonizeMatrices @ beta
		(* CanonizeMatrices @ RefineGroupStructures @ beta *)
	];

(*Function that produces the beta function for the requested coupling*)
Options[BetaFunction] = {
		RescaledCouplings-> False,
		FourDimensions-> True,
		FlavorImproved-> True
	};
BetaFunction[coupling_, loop_, opt:OptionsPattern[] ] ? OptionsCheck :=
	Module[{coef = 4 Pi, firstTerm = 0, l},
		CheckKnownBeta[coupling, loop];
		If[MemberQ[$betaTypes, coupling], Return@ ProjectionToUnmixedBetas[coupling, loop, BetaFunction, opt];];

		CheckPotentialCouplingMixing[coupling, BetaFunction];

		If[OptionValue @ RescaledCouplings, coef = 1; ];
		If[OptionValue @ FourDimensions, firstTerm = 1; ];

		Sum[ Power[coef, -2 l] Quiet@ BetaTerm[coupling, l, FilterRules[{opt}, Options@ BetaTerm]],
			{l, firstTerm, loop}]
	];

(* Function for checking possible coupling mixings *)
CheckPotentialCouplingMixing::mixes= "The quartic couplings `1` involve the same fields. This function returns an unpredictable linear combination of all the associated \[Beta]-functions. Please use the type call `2`[`3`, <loop>] instead."
CheckPotentialCouplingMixing[coupling_, func_]:=
	Module[{fields, couplings, couplingsAssoc, type},
		type= CouplingType@ coupling;
		If[MatchQ[type, Gauge| FermionMass| ScalarMass], Return[]; ];

		couplingsAssoc= CouplingInfo@ type;

		fields= Sort@ couplingsAssoc[coupling, Fields];
		couplings= Keys@ Select[couplingsAssoc, (Sort@ #[Fields] === fields &)];

		If[Length@ couplings> 1,
			Message[CheckPotentialCouplingMixing::mixes, couplings, func, type];
		];
	];

CouplingInfo@ couplingType_:= Switch[couplingType,
		Gauge, $gaugeGroups,
		Yukawa, $yukawas,
		Quartic, $quartics,
		Trilinear, $trilinears,
		ScalarMass, $scalarMasses,
		FermionMass, $fermionMasses
	];

(*Fuction returning the diagonalized list of all beta functions (terms) of a given type. It inherits the options from BetaFunction or BetaTerm*)
ProjectionToUnmixedBetas::singular = "The projection matrix of `1` is singular. Some of the couplings are redundant."
ProjectionToUnmixedBetas[type_, loop_, func_, opt:OptionsPattern[]]:=
	Module[{beta, couplings, projections, c},
		couplings= If[type === Gauge,
			List @@ Query[All, Key@ Coupling]@ $gaugeGroups,
			Keys@ CouplingInfo@ type
		];
		Print["The ", type, " couplings are ", couplings, " and the associated \[Beta]-functions(terms) are"];

		If[Length@ couplings === 0, Return@ {}; ];

		(*Finds inversion matrix for the projectors*)
		If[MatchQ[type, Yukawa| Quartic| Trilinear],
			projections= CheckProjection /@ couplings/. {Matrix[c_][__]-> c, Tensor[c_][__]-> c};
			invMatrix= (Simplify@ D[projections, #]&)/@ couplings// Transpose;
			If[Simplify@ Det@ invMatrix === 0,
				Message[ProjectionToUnmixedBetas::singular, type];
				Abort[];
			];
			invMatrix= Simplify@ Inverse @ invMatrix;
		,
			invMatrix= IdentityMatrix@ Length@ couplings;
		];

		beta = Monitor[
				Table[Quiet@ func[c, loop, opt], {c, couplings}]
			, StringForm["Evaluating the `` \[Beta]-function", c] ];

		invMatrix . beta // Expand

	];


(*Fuction for diagonalizing the quartic beta functions. It inherits the options from BetaFunction*)
QuarticBetaFunctions::singular = "The projection matrix is singular. Some of the couplings may be redundant."
QuarticBetaFunctions::loopNumber = "The quartic beta function has only been implemented up to 3 loops."
QuarticBetaFunctions[loop_Integer, opt:OptionsPattern[]] ? OptionsCheck :=
	Module[{betaFunctions, couplings, qProjections, invMatrix, c},
		If[loop > 3,
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
							Table[Quiet@ BetaFunction[c, loop, opt], {c, couplings}]
						,StringForm["Evaluating the `` \[Beta]-function", c] ];

		invMatrix . betaFunctions // Expand
	];

(*Function returns the l-loop contribution to the anomalous dimension of a given matter field*)
AnomalousDimTerm::loopNumber = "The `1` scalar anomalous dimension has only been implemented up to `2` loops."
AnomalousDimTerm::incompatible = "The two input fields do not share quantum numbers."
AnomalousDimTerm::unkown = "The field `1` has not been defined."
AnomalousDimTerm[inp_/; MemberQ[{0, 2}, Length@ inp], loop_Integer] :=
	Module[{anomDim, fields},
		fields= If[Head@ inp === List, inp, {inp, inp}];

		(*Determines the correct tensor structure based on the field type*)
		Switch[First@ fields
		,x_ /; MemberQ[Keys @ $fermions, x],
			If[loop > 2 || loop < 1,
				Message[AnomalousDimTerm::loopNumber, "fermion", 2];
				Abort[];
			];
			(* Checks that the qunatum numbers agree *)
			If[Sort@ $fermions[fields[[1]], GaugeRep] =!= Sort@ $fermions[fields[[2]], GaugeRep],
				Message[AnomalousDimTerm::incompatible];
				Abort[];
			];
			anomDim = FermionAnomalousTensors[fields, loop] /. $fermionAnomalousCoefficients;

		,x_ /; MemberQ[Keys @ $scalars, x],
			If[loop > 2 || loop < 1,
				Message[AnomalousDimTerm::loopNumber, "scalar", 2];
				Abort[];
			];
			(* Checks that the qunatum numbers agree *)
			If[Sort@ $scalars[fields[[1]], GaugeRep] =!= Sort@ $scalars[fields[[2]], GaugeRep],
				Message[AnomalousDimTerm::incompatible];
				Abort[];
			];
			anomDim = ScalarAnomalousTensors[fields, loop] /. $scalarAnomalousCoefficients;

		,_,
			Message[AnomalousDimTerm::unkown, First@ fields];
			Abort[];
		];

		(*Canonically order coupling indices if any*)
		CanonizeMatrices @ RefineGroupStructures @ anomDim
	];

(*Function that produces the anomalous dimension for a given field*)
Options[AnomalousDimension] = {RescaledCouplings -> False};
AnomalousDimension::unkown = "The coupling `1` has not been defined."
AnomalousDimension::loopNumber = "The anomalous dimension has only been implemented up to `1` loops."
AnomalousDimension[inp_/; MemberQ[{0, 2}, Length@ inp], loop_Integer, OptionsPattern[] ] ? OptionsCheck :=
	Module[{coef = 4 Pi, l},
		(* Checks if the input request has been implemented.*)
		Switch[If[Head@ inp === List, First@ inp, inp]
		,x_ /; MemberQ[Keys @ $fermions, x] || MemberQ[Keys @ $scalars, x],
			If[loop > 2 || loop < 1,
				Message[AnomalousDimension::loopNumber, 2];
				Abort[];
			];
		,_,
			Message[AnomalousDimension::unkown, inp];
			Abort[];
		];

		If[OptionValue @ RescaledCouplings, coef = 1;];

		Sum[ Power[coef, -2 l] AnomalousDimTerm[inp, l], {l, loop}]
	];

(*Function returns the l-loop contribution to the upsilon function of a given matter field*)
UpsilonTerm::loopNumber= "The `1` scalar anomalous dimension has only been implemented up to `2` loops."
UpsilonTerm::incompatible = "The two input fields do not share quantum numbers."
UpsilonTerm::unkown= "The field `1` has not been defined."
UpsilonTerm[inp_/; MemberQ[{0, 2}, Length@ inp], loop_Integer]:=
	Module[{ups, fields},
		fields= If[Head@ inp === List, inp, {inp, inp}];

		(*Determines the correct tensor structure based on the field type*)
		Switch[First@ fields
		,x_ /; MemberQ[Keys@ $fermions, x],
			If[loop > 3 || loop < 1,
				Message[UpsilonTerm::loopNumber, "fermion", 3];
				Abort[];
			];
			(* Checks that the qunatum numbers agree *)
			If[Sort@ $fermions[fields[[1]], GaugeRep] =!= Sort@ $fermions[fields[[2]], GaugeRep],
				Message[UpsilonTerm::incompatible];
				Abort[];
			];
			ups= FermionUpsilonTensors[fields, loop]/. $fermionUpsilonCoefficients;

		,x_ /; MemberQ[Keys@ $scalars, x],
			If[loop > 3 || loop < 1,
				Message[UpsilonTerm::loopNumber, "scalar", 3];
				Abort[];
			];
			(* Checks that the qunatum numbers agree *)
			If[Sort@ $scalars[fields[[1]], GaugeRep] =!= Sort@ $scalars[fields[[2]], GaugeRep],
				Message[UpsilonTerm::incompatible];
				Abort[];
			];
			ups= ScalarUpsilonTensors[fields, loop]/. $scalarUpsilonCoefficients;

		,_,
			Message[UpsilonTerm::unkown, First@ fields];
			Abort[];
		];

		(*Canonically order coupling indices if any*)
		CanonizeMatrices @ RefineGroupStructures @ ups
	];

(*Function that produces the upsilon function for a given field*)
Options@ UpsilonFunction= {RescaledCouplings-> False};
UpsilonFunction::unkown = "The coupling `1` has not been defined."
UpsilonFunction::loopNumber = "The anomalous dimension has only been implemented up to `1` loops."
UpsilonFunction[inp_/; MemberQ[{0, 2}, Length@ inp], loop_Integer, OptionsPattern[] ] ? OptionsCheck :=
	Module[{coef = 4 Pi, l, type},
		(* Checks if the input request has been implemented.*)
		Switch[If[Head@ inp === List, First@ inp, inp]
		,x_ /; MemberQ[Keys @ $fermions, x] || MemberQ[Keys @ $scalars, x],
			If[loop > 3 || loop < 1,
				Message[UpsilonFunction::loopNumber, 3];
				Abort[];
			];
		,_,
			Message[UpsilonFunction::unkown, inp];
			Abort[];
		];

		If[OptionValue @ RescaledCouplings, coef = 1;];

		Sum[Power[coef, -2 l] UpsilonTerm[inp, l], {l, loop}]
	];

(*###################################*)
(*----------Extra functions----------*)
(*###################################*)

(*Function for finalizing a betafunction, bringing it from a nice compact output to a form more suitable for
further Mathematica manipulations. Can also be used to specify particular cases for coupling matrices.*)
Options[Finalize] = {Parameterizations -> {}, BarToConjugate -> False};
Finalize[expr_, opt:OptionsPattern[]] ? OptionsCheck :=
	Module[{out},
		out = CanonizeMatrices @ expr /. OptionValue @ Parameterizations;
		If[OptionValue @ BarToConjugate,
			out = out/. Bar[a_Symbol] :> Conjugate@ a;
		];
		out/. Matrix[y__][__] :> Dot[y] /. {del[_, $i, $j]-> 1}
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
