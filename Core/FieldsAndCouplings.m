(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Begin["FieldsAndCouplings`"]

(*Structure deltas*)
sDelS /: sDelS[field1_, ind_, f_] sDelS[field2_, ind_, g_] := 0;
sDelF /: sDelF[field1_, ind_, f_] sDelF[field2_, ind_, g_] := 0;
sDelV /: sDelV[field1_, ind_, f_] sDelV[field2_, ind_, g_] := 0;
(*For masses*)
sDelS /: sDelS[$vev, ind_, f_] sDelS[$vevSelect, ind_, g_] := 1;
Bar@ $vev = $vev;


(*Initiates a scalar field*)
Options[AddScalar] = {SelfConjugate -> False, GaugeRep -> {}, FlavorIndices -> {}, Mass -> None};
AddScalar::failure = "Failed to add scalar field.";
AddScalar[field_, OptionsPattern[] ] ? OptionsCheck :=
	Block[{rep, massTerm = 1, projector},
		(*Checks gauge representations*)
		If[! And @@ Table[RepresentationCheck @ rep, {rep, OptionValue[GaugeRep]}],
			Message[AddScalar::failure];
			Return @ $Failed;
		];

		(* To project out the anomalous dimension of the field from the general structure *)
		projector = Evaluate[sDelS[Bar @ field, #1, s1] sDelS[field, #2, s2] Product[del[rep, s1, s2] / Dim @ rep, {rep, OptionValue[GaugeRep]}] ]&;

		(*Adds field options to the list of scalar fields*)
		AppendTo[$scalars, field -> <|
			GaugeRep -> OptionValue[GaugeRep],
			FlavorIndices -> OptionValue[FlavorIndices],
			Projector -> projector,
			SelfConjugate -> OptionValue[SelfConjugate],
			Mass -> OptionValue @ Mass|>];

		If[OptionValue @ Mass =!= None,
			massTerm = UnitStep[Global`t - OptionValue @ Mass];
		];

		(*Initializes the structure deltas for the field*)
		If[OptionValue[SelfConjugate],
			sDelS/: sDelS[field, ind_, s1_] sDelS[field, ind_, s2_] = 1 * massTerm
				* Product[del[rep, s1, s2], {rep, OptionValue[GaugeRep]}]
				* Product[del[rep, s1, s2], {rep, OptionValue[FlavorIndices]}];
			Bar[field] = field;
		,
			sDelS/: sDelS[field, ind_, s1_] sDelS[Bar[field], ind_, s2_] = 2 * massTerm
				* Product[del[rep, s1, s2], {rep, OptionValue[GaugeRep]}]
				* Product[del[rep, s1, s2], {rep, OptionValue[FlavorIndices]}];
		];

		ResetBetas[];
	];

(*Initiates a fermion field*)
Options[AddFermion] = {GaugeRep -> {}, FlavorIndices -> {}, Mass -> None};
AddFermion::failure = "Failed to add fermion field.";
AddFermion[field_, OptionsPattern[] ] ? OptionsCheck :=
	Block[{massTerm = 1, rep},
		(*Checks gauge representations*)
		If[! And @@ Table[RepresentationCheck @ rep, {rep, OptionValue[GaugeRep]}],
			Message[AddFermion::failure];
			Return @ $Failed;
		];

		(* To project out the anomalous dimension of the field from the general structure *)
		projector = Evaluate[sDelF[Bar @ field, #1, f1] sDelF[field, #2, f2] Product[del[rep, f1, f2] / Dim @ rep, {rep, OptionValue[GaugeRep]}] ]&;

		(*Adds field options to the list of fermion fields*)
		AppendTo[$fermions, field -> <|
			GaugeRep -> OptionValue[GaugeRep],
			FlavorIndices -> OptionValue[FlavorIndices],
			Projector -> projector,
			Mass -> OptionValue @ Mass|>];

		If[OptionValue @ Mass =!= None,
			massTerm = UnitStep[Global`t - OptionValue @ Mass];
		];

		(*Initializes the structure deltas for the field*)
		sDelF/: sDelF[field, ind_, f1_] sDelF[Bar[field], ind_, f2_] = massTerm
			* Product[del[rep, f1, f2], {rep, OptionValue[GaugeRep]}]
			* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];

		ResetBetas[];
	];

(*Initiates a vector field*)
AddVector[name_, group_] := (sDelV/: sDelV[name, ind_, v1_] sDelV[name, ind_, v2_] = del[group[adj], v1, v2]);


(*###################################*)
(*----------Other functions----------*)
(*###################################*)
(* Determines all unique permutations of fields with a given invariant *)
CouplingPermutations[fields_List, invariant_Function] :=
	Block[{arg, fs, inv, number, perms, uniformSymbols, uniqueArrangements},
		uniformSymbols = s_Symbol /; StringMatchQ[SymbolName@s, "*$*"] :> Symbol @ StringSplit[SymbolName@s, "$"][[1]];
		number = Length @ fields;
		perms = Permutations @ Table[n, {n, number}];
		uniqueArrangements = {};
		Reap[
			Do[
				fs = fields[[Ordering @ p]];
				inv = (invariant @@ arg /@ p) /.uniformSymbols /.
					{del[rep_, a__] :> del[rep, Sequence @@ Sort @ List @ a],
					eps[rep_, a__] :> eps[rep, Sequence @@ Sort @ List @ a],
					lcSymb[rep_, a__] :> lcSymb[rep, Sequence @@ Sort @ List @ a],
					delA2[rep_, i_, a__] :> delA2[rep, i, Sequence @@ Sort @ List @ a],
					delS2[rep_, i_, a__] :> delS2[rep, i, Sequence @@ Sort @ List @ a]};
				If[MemberQ[uniqueArrangements, {fs, inv}],
					Continue[];
				];
				AppendTo[uniqueArrangements, {fs, inv}];
				Sow @ p;
			,{p, perms}];
		][[2,1]]
	];

AveragePermutations[indices_List, permutations_List][expr_] :=
	Block[{subs},
		subs = MapThread[Rule, {indices, indices[[#]]}] & /@ permutations;
		Mean[expr/.subs]
	];

(* Puts a tensor on canonical form---Brute force *)
CanonizeTensor[symbForm_List][expr_] :=
	Block[{symb, newSymbs, indexAssignments, temp},
		symb = DeleteDuplicates @ Cases[expr, _Symbol?(StringMatchQ[SymbolName@ #, Alternatives @@ symbForm] &), Infinity];
		newSymbs = Table[Symbol["internal$" <> ToString@ n], {n, Length @ symb}];
		indexAssignments = (Thread@Rule[#, newSymbs] &) /@ Permutations @ symb;
		Sort[Flatten @ Table[
			expr /. assign /. {eps[rep_, a_, c_] tGen[rep_, A_, c_, b_] /; !OrderedQ[{a, b}] -> eps[rep, b, c] tGen[rep, A, c, a],
			eps[rep_, c_, a_] tGen[rep_, A_, c_, b_] /; !OrderedQ@{a, b} -> eps[rep, c, b] tGen[rep, A, c, a]} /.
			{eps[rep_, a_, b_] /; ! OrderedQ@{a, b} -> -eps[rep, b, a],
			del[rep_, a_, b_] /; ! OrderedQ@{a, b} -> del[rep, b, a] }
		, {assign, indexAssignments}]][[1]]
	];
CanonizeTensor[symbForm_List][Plus[x_, y__]] := CanonizeTensor[symbForm]@ x + CanonizeTensor[symbForm] @ Plus@ y;

(*###################################*)
(*----------Gauge couplings----------*)
(*###################################*)
(*Function for adding gauge groups to the model*)
AddGaugeGroup::unkown = "`1` is not a reckognized Lie group.";
AddGaugeGroup::nonstring = "Use a string for the field or leave it as the default value.";
AddGaugeGroup::dimensions = "The dimension of the coupling matrix does not match the dimension of
	the kineitc mixing term";
AddGaugeGroup::automatic = "Automatic naming of the coupling matrix only suitable for U(1)^p for p <= 9";
Options[AddGaugeGroup] = {CouplingMatrix -> Automatic, Field -> Automatic};
AddGaugeGroup[coupling_Symbol, groupName_Symbol, U1, opts:OptionsPattern[]] :=
	AddGaugeGroup[coupling, groupName, U1[1], opts];
AddGaugeGroup[coupling_Symbol, groupName_Symbol, lieGroup_Symbol[n_Integer|n_Symbol], OptionsPattern[]] ? OptionsCheck :=
	Block[{cMatrix, invalid, projector, fieldName},
		(*Checks for mismatch between U1 power and coupling matrix*)
		If[MatchQ[lieGroup[n], U1[x_] /; x > 1] && OptionValue @ CouplingMatrix =!= Automatic,
			If[!SymmetricMatrixQ @ OptionValue @ CouplingMatrix || Length @ OptionValue @ CouplingMatrix =!= n,
				Message[AddGaugeGroup::dimensions];
				Return @ $Failed;
			];
		];

		invalid = DefineLieGroup[groupName, lieGroup[n]];
		If[invalid === $Failed,
			Return @ $Failed;
		];

		(*Decides on the field name*)
		Switch[OptionValue @ Field
		,Automatic,
			fieldName = "A_" <> ToString @ groupName;
		,_,
			fieldName = OptionValue @ Field;
		];

		(*Sets up the gauge fields and 2-point projection*)
		AddVector[fieldName, groupName];
		projector = Evaluate[If[lieGroup =!= U1, del[groupName[adj], v1, v2] / Dim @ groupName[adj], 1] *
			sDelV[fieldName, #1, v1] sDelV[fieldName, #2, v2] ] &;

		(*Adds the group information and the coupling to the repsective lists*)
		AppendTo[$gaugeGroups, groupName ->
			<|Coupling -> coupling,
			Field -> fieldName,
			LieGroup -> lieGroup[n],
			Projector -> projector|>];
		AppendTo[$couplings, coupling -> groupName];

		(*Setting up the coupling matrix for Abelian groups with Kinetic mixing*)
		If[lieGroup === U1 && n > 1,
			If[OptionValue @ CouplingMatrix === Automatic,
				If[n > 9,
					Message[AddGaugeGroup::automatic];
				];
				cMatrix = Table[Symbol[ToString @ coupling <>
						Which[i < j, ToString @ i <> ToString @ j,
							i===j, ToString @ i,
							j < i, ToString @ j <> ToString @ i] ]
					,{i, n}, {j, n} ];
			,
				cMatrix = OptionValue @ CouplingMatrix;
			];
			AppendTo[$gaugeGroups @ groupName, CouplingMatrix -> cMatrix];
			coupling /: Matrix[a_List, coupling, b___] = Matrix[Dot[a, cMatrix], b];
			coupling /: Matrix[a___, coupling, b_List] = Matrix[a, Dot[cMatrix, b]];
			Trans[coupling] = coupling;
		];

		ResetBetas[];
	];

(*Checks the validity of a given representation. Returns True iff rep is valid given gauge groups
 of the present model*)
RepresentationCheck::invalid = "`1` is not a reckognized format for a representation.";
RepresentationCheck::representation = "`1` is not a reckognized representation of a `2` group.";
RepresentationCheck::noGroup = "The `1` gauge group has not been defined.";
RepresentationCheck[rep_] :=
	Block[{group, gName},
		If[! MatchQ[rep, _@_], (*Checks form*)
			Message[RepresentationCheck::invalid, rep];
			Return @ False;
		];

		gName = Head@ If[Head@ rep === Bar, rep[[1]], rep];
		group = $gaugeGroups[gName, LieGroup];
		Switch[group
		,SU[_],
			If[MemberQ[Join[gName/@{fund, adj, S2, A2}, Bar/@ gName/@ {fund, S2, A2}], rep],
				Return @ True;
			];
			Message[RepresentationCheck::representation, rep, group];
		,SO[_],
			If[MemberQ[gName/@{fund, adj, S2}, rep],
				Return @ True;
			];
			Message[RepresentationCheck::representation, rep, group];
		,Sp[_],
			If[MemberQ[Join[gName/@{fund, adj, A2}, Bar/@ gName/@ {fund, A2}], rep],
				Return @ True;
			];
			Message[RepresentationCheck::representation, rep, group];
		,U1[1],
			If[Head @ rep =!= Bar && Head @ rep[[1]] =!= List &&
					!MemberQ[{adj, fund, A2, S2}, rep[[1]] ],
				Return @ True;
			];
			Message[RepresentationCheck::representation, rep, group];
		,U1[_],
			If[Head @ rep =!= Bar && VectorQ @ rep[[1]] && Length @ rep[[1]] === group[[1]],
				Return @ True;
			];
			Message[RepresentationCheck::representation, rep, group];
		,_,
			Message[RepresentationCheck::noGroup, Head @rep];
			Return @ False;
		];
		Return @ False;
	];

(*The general gauge coupling matrix G^2_{AB} used in the computation of the beta function tensors*)
G2Matrix[A_, B_] :=
	Module[{gauge, v1, v2},
		Sum[
			sDelV[$gaugeGroups[gauge, Field], A, v1] sDelV[$gaugeGroups[gauge, Field], B, v2]
				* If[MatchQ[$gaugeGroups[gauge, LieGroup], U1[n_] /; n > 1],
					Matrix[$gaugeGroups[gauge, Coupling]][gauge[adj] @v1, gauge[adj] @v2]
				,
					$gaugeGroups[gauge, Coupling]^2 del[gauge[adj], v1, v2]
				]
		,{gauge, Keys @ $gaugeGroups}]
	];

(*Defining the gauge generators for the left-handed spinors, \bar{\psi} T^A \psi.*)
TfLeft[A_, i_, j_] :=
	Module[{ferm, group, rep, gRep1, gRep2, f1, f2, v1},
		Sum[
			sDelF[Bar @ ferm, i, f1] sDelF[ferm, j, f2] * Product[del[rep, f1, f2], {rep, $fermions[ferm, FlavorIndices]}]
			* Sum[group = Head@ If[Head@ gRep1 === Bar, gRep1[[1]], gRep1];
					sDelV[$gaugeGroups[group, Field], A, v1] tGen[gRep1, v1, f1, f2]
					* Product[del[gRep2, f1, f2], {gRep2, DeleteCases[$fermions[ferm, GaugeRep], gRep1]}],
				{gRep1, $fermions[ferm, GaugeRep]}]
		,{ferm, Keys @ $fermions}]
	];

(*The general fermion gauge generators used in the computation of the beta function tensors*)
Tferm[A_, i_, j_] := {{TfLeft[A, i, j], 0}, {0, -TfLeft[A, j, i]}};
TfermTil[A_, i_, j_] := {{-TfLeft[A, j, i], 0}, {0, TfLeft[A, i, j]}};

(*The general scalar gauge generators used in the computation of the beta function tensors*)
Tscal[A_, a_, b_] :=
	Module[{scal, group, rep, gRep1, gRep2, s1, s2, v1},
		Sum[
			AntiSym[a, b][sDelS[Bar @ scal, a, s1] sDelS[scal, b, s2]]
			* Product[del[rep, s1, s2], {rep, $scalars[scal, FlavorIndices]}]
			* Sum[group = Head@ If[Head@ gRep1 === Bar, gRep1[[1]], gRep1];
					sDelV[$gaugeGroups[group, Field], A, v1] tGen[gRep1, v1, s1, s2]
					* Product[del[gRep2, s1, s2], {gRep2, DeleteCases[$scalars[scal, GaugeRep], gRep1]}],
				{gRep1, $scalars[scal, GaugeRep]}]
		,{scal, Keys @ $scalars}]
	];

(*The general gauge structure constants G^{-2}_{AD} f^{DBC} used in the computation of the beta function tensors*)
FGauge[A_, B_, C_] :=
	Module[{gauge, v1, v2, v3},
		Sum[
			sDelV[$gaugeGroups[gauge, Field], A, v1] sDelV[$gaugeGroups[gauge, Field], B, v2] sDelV[$gaugeGroups[gauge, Field], C, v3]
				* Power[$gaugeGroups[gauge, Coupling], -2] fStruct[gauge, v1, v2, v3]
		,{gauge, Keys @ $gaugeGroups}]
	];

(*####################################*)
(*----------Yukawa couplings----------*)
(*####################################*)
(*Function for defining the Yukawa couplings of the theory*)
Options[AddYukawa] = {CouplingIndices -> (Null &),
	GroupInvariant -> (1 &),
	Chirality -> Left,
	CheckInvariance -> False};
AddYukawa::unkown = "`1` does not match any of the `2`s.";
AddYukawa::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors."
AddYukawa::gaugeInv = "Could not verify invariance under the `1` group. Plase check the GroupInvariant for errors.";
AddYukawa[coupling_, {phi_, psi1_, psi2_}, OptionsPattern[]] ? OptionsCheck:=
	Block[{g, group, normalization, projection, symmetryFactor, temp, test, yuk, yukawa, yukbar, y},
		(*Tests if the fields have been defined*)
		If[!MemberQ[Keys @ $scalars, phi /. Bar[x_] -> x],
			Message[AddYukawa::unkown, phi /. Bar[x_] -> x, "scalar"];
			Return @ $Failed;
		];
		If[!MemberQ[Keys@ $fermions, psi1],
			Message[AddYukawa::unkown, psi1, "fermion"];
			Return @ $Failed;
		];
		If[!MemberQ[Keys@ $fermions, psi2],
			Message[AddYukawa::unkown, psi2, "fermion"];
			Return @ $Failed;
		];

		(*If the chirality is right handed, the coupling is written with the barred fields*)
		Switch[OptionValue @ Chirality
		,Left,
			g = coupling;
		,Right,
			g = Bar @ coupling;
		];

		(*Constructs the coupling structure*)
		normalization = 2; (*To account for the structure deltas working with real/complex scalars in Yuk[a, i, j]*)
		If[!$scalars[phi, SelfConjugate]|| Head @ phi === Bar, normalization *= 1/Sqrt[2];];

		If[Length @ OptionValue[CouplingIndices][s1, f1, f2] <= 2,
			yuk = Evaluate[normalization Matrix[g] @ ##] &;
			yukbar = Evaluate[normalization Matrix[Bar @ g] @ ##] &
		,
			yuk = Evaluate[normalization Tensor[g] @ ##] &;
			yukbar = Evaluate[normalization Tensor[Bar @ g] @ ##] &
		];

		(*Adds the Yukawa coupling to the association*)
		AppendTo[$yukawas, coupling ->
			<|Chirality -> OptionValue @ Chirality,
			Coupling -> yuk,
			CouplingBar -> yukbar,
			Fields -> {phi, psi1, psi2},
			Indices -> OptionValue @ CouplingIndices,
			Invariant -> OptionValue @ GroupInvariant|>];

		(*Constructs the projection operator*)
		projection = Evaluate[ ReplaceAll[tGen[rep_, A_, a_, b_] -> tGen[Bar@rep, A, a, b]] @ OptionValue[GroupInvariant][s1, f1, f2] *
			Switch[OptionValue @ Chirality
			,Left,
				sDelS[Bar@ phi, #1, s1] sDelF[Bar@ psi1, #2, f1] sDelF[Bar@ psi2, #3, f2]
			,Right,
				sDelS[phi, #1, s1] sDelF[psi1, #2, f1] sDelF[psi2, #3, f2]
			] ] &;
		symmetryFactor = ReplaceAll[Rule[#, 0] & /@ DeleteCases[Keys @ $yukawas, coupling]] @ RefineGroupStructures @
			Expand[projection[$a, $i, $j] Switch[OptionValue @ Chirality, Left, YukawaLeft[$a, $i, $j], Right, YukawaRight[$a, $i, $j]	] ] /.
				(Matrix|Tensor)[x_][__] -> x /. coupling -> 1 // Simplify;
		If[symmetryFactor === 0,
			Message[AddYukawa::projection0];
			KeyDropFrom[$yukawas, coupling];
			Return @ $Failed;
		];
		projection = Evaluate[projection[#1, #2, #3] / symmetryFactor] &;
		AppendTo[$yukawas @ coupling, Projector -> projection];

		(*Tests whether the coupling satisfy gauge invariance*)
		If[OptionValue @ CheckInvariance,
			Block[{A,i1,i2,i3,a1,a2},
				yukawa[$da_, $di_, $dj_] = sDelS[phi, $da, s1$1] sDelF[psi1, $di, f1$1] sDelF[psi2, $dj, f2$1] OptionValue[GroupInvariant][s1$1, f1$1, f2$1];
				test = yukawa[a2, i1, i2] Tscal[A, a2, a1] + yukawa[a1, i3, i2] TfLeft[A, i3, i1] + yukawa[a1, i1, i3] TfLeft[A, i3, i2] //Expand;
				(* test = TfLeft[A, k, i] YukawaLeft[a, k ,j] + YukawaLeft[a, i, k] TfLeft[A, k, j] + YukawaLeft[b, i, j] Tscal[A, b, a] //Expand; *)
				test *= sDelS[Bar@phi, a1, s1] sDelF[Bar@psi1, i1, f1] sDelF[Bar@psi2, i2, f2] //Expand;
				Do[
					temp = CanonizeTensor[{"s1$*","s2$*","s3$*","s4$*","f1$*","f2$*"}][ test sDelV[$gaugeGroups[group, Field], A, v1] //Expand ] /. Matrix[x_][_[v1]] -> x;
					If [temp =!= 0 && !MatchQ[temp, {0 ..}],
						Message[AddYukawa::gaugeInv, group];
						Print["0 != ", temp];
					];
				,{group, Keys@ $gaugeGroups}];
			];
		];

		AppendTo[$couplings, coupling -> Yukawa];
		ResetBetas[];
	];


(*Function for defining the fermion masses of the theory*)
Options[AddFermionMass] = {MassIndices -> (Null &),
	GroupInvariant -> (1 &),
	Chirality -> Left};
AddFermionMass::unkown = "`1` does not match any of the `2`s.";
AddFermionMass::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors."
AddFermionMass[mass_, {psi1_, psi2_}, OptionsPattern[]] ? OptionsCheck:=
	Block[{g, group, projection, symmetryFactor, temp, test, yuk, yukbar, y},
		(*Tests if the fields have been defined*)
		If[!MemberQ[Keys@ $fermions, psi1],
			Message[AddFermionMass::unkown, psi1, "fermion"];
			Return @ Null;
		];
		If[!MemberQ[Keys@ $fermions, psi2],
			Message[AddFermionMass::unkown, psi2, "fermion"];
			Return @ Null;
		];

		(*If the chirality is right handed, the coupling is written with the barred fields*)
		Switch[OptionValue @ Chirality
		,Left,
			g = mass;
		,Right,
			g = Bar @ mass;
		];

		(*Constructs the coupling structure*)
		(*The normalization is 2 to account for the symmetrization in Yukawa[a, i, j]*)
		If[Length @ OptionValue[MassIndices][f1, f2] <= 2,
			yuk = Evaluate[2 Matrix[g] @ ##] &;
			yukbar = Evaluate[2 Matrix[Bar @ g] @ ##] &
		,
			yuk = Evaluate[2 Tensor[g] @ ##] &;
			yukbar = Evaluate[2 Tensor[Bar @ g] @ ##] &
		];

		(*Adds the Yukawa coupling to the association*)
		AppendTo[$fermionMasses, mass ->
			<|Chirality -> OptionValue @ Chirality,
			Coupling -> yuk,
			CouplingBar -> yukbar,
			Fields -> {psi1, psi2},
			Indices -> OptionValue @ MassIndices,
			Invariant -> OptionValue @ GroupInvariant|>];

		(*Constructs the projection operator*)
		projection = Evaluate[ ReplaceAll[tGen[rep_, A_, a_, b_] -> tGen[Bar@rep, A, a, b]] @ OptionValue[GroupInvariant][f1, f2] *
			Switch[OptionValue @ Chirality
			,Left,
				sDelS[$vevSelect, #1, s1] sDelF[Bar@ psi1, #2, f1] sDelF[Bar@ psi2, #3, f2]
			,Right,
				sDelS[$vevSelect, #1, s1] sDelF[psi1, #2, f1] sDelF[psi2, #3, f2]
			] ] &;
		symmetryFactor = ReplaceAll[Rule[#, 0] & /@ DeleteCases[Keys @ $fermionMasses, mass]] @ RefineGroupStructures @
			Expand[projection[$a, $i, $j] Switch[OptionValue @ Chirality, Left, YukawaLeft[$a, $i, $j, True],
				Right, YukawaRight[$a, $i, $j, True] ] ] /. (Matrix|Tensor)[x_][__] -> x /. mass -> 1 // Simplify;
		If[symmetryFactor === 0,
			Message[AddFermionMass::projection0];
			KeyDropFrom[$fermionMasses, mass];
			Return @ $Failed;
		];
		projection = Evaluate[projection[#1, #2, #3] / symmetryFactor] &;
		AppendTo[$fermionMasses @ mass, Projector -> projection];
		AppendTo[$couplings, mass -> FermionMass];

		ResetBetas[];
	];


(*Chiral Yukawa couplings*)
YukawaLeft[$da_, $di_, $dj_, massive_:False] :=
	Module[{f1, f2, yu, s1},
		Sym[$di, $dj][
			Sum[sDelS[yu[Fields][[1]], $da, s1] sDelF[yu[Fields][[2]], $di, f1] sDelF[yu[Fields][[3]], $dj, f2]
					* yu[Invariant][s1, f1, f2] yu[Coupling][Sequence @@ yu[Indices][s1, f1, f2]]
				,{yu, $yukawas}] +
			If[massive,
				Sum[sDelS[$vev, $da, s1] sDelF[yu[Fields][[1]], $di, f1] sDelF[yu[Fields][[2]], $dj, f2]
					* yu[Invariant][f1, f2] yu[Coupling][Sequence @@ yu[Indices][f1, f2]]
				,{yu, $fermionMasses}]
			,0]
		]
	];

YukawaRight[$da_, $di_, $dj_, massive_:False] :=
	Module[{f1, f2, yu, s1},
		Sym[$di, $dj][
			Sum[sDelS[Bar @ yu[Fields][[1]], $da, s1] sDelF[Bar @ yu[Fields][[2]], $di, f1] sDelF[Bar @ yu[Fields][[3]], $dj, f2]
					* yu[Invariant][s1, f1, f2] yu[CouplingBar][Sequence @@ yu[Indices][s1, f1, f2]]
				,{yu, $yukawas}] +
			If[massive,
				Sum[sDelS[$vev, $da, s1] sDelF[Bar @ yu[Fields][[1]], $di, f1] sDelF[Bar @ yu[Fields][[2]], $dj, f2]
					* yu[Invariant][f1, f2] yu[CouplingBar][Sequence @@ yu[Indices][f1, f2]]
				,{yu, $fermionMasses}]
			,0]
		]
	];

(*Genral Yukawa couplings used in the computation of the beta function tensors.*)
Yuk[a_, i_, j_, massive_:False] := {{YukawaLeft[a, i, j, massive], 0}, {0, YukawaRight[a, i, j, massive]}};
YukTil[a_, i_, j_, massive_:False] := {{YukawaRight[a, i, j, massive], 0}, {0, YukawaLeft[a, i, j, massive]}};


(*#####################################*)
(*----------Quartic couplings----------*)
(*#####################################*)

(*Function for defining the quartic couplings of the theory*)
AddQuartic::unkown = "`1` does not match any of the scalars.";
AddQuartic::gaugeInv = "Could not verify invariance under the `1` group. Plase check the GroupInvariant for errors.";
AddQuartic::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors.";
Options[AddQuartic] = {CouplingIndices -> (Null &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True,
	CheckInvariance -> False};
AddQuartic [coupling_, {phi1_, phi2_, phi3_, phi4_}, OptionsPattern[]] ? OptionsCheck :=
	Block[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
		(*Tests if the fields have been defined*)
		Do[
			If[!MemberQ[Keys@ $scalars, temp /. Bar[x_] -> x],
				Message[AddQuartic::unkown, temp /. Bar[x_] -> x];
				Return[$Failed];
			];
		,{temp, {phi1, phi2, phi3, phi4}}];

		If[OptionValue @ SelfConjugate,
			SetReal @ coupling;
		];

		normalization = 24; (*To account for the structure deltas working with real/complex scalars in Lam[a, b, c, d]*)
		Do[
			If[!$scalars[phi, SelfConjugate]|| Head @ phi === Bar, normalization *= 1/Sqrt[2];];
		,{phi, {phi1, phi2, phi3, phi4}}];

		If[Length @ OptionValue[CouplingIndices][s1, s2, s3, s4] <= 2,
			lam = Evaluate[normalization Matrix[coupling] @ ##] &;
			lambar = Evaluate[normalization Matrix[Bar @ coupling] @ ## ] &
		,
			lam = Evaluate[normalization Tensor[coupling] @ ##] &;
			lambar = Evaluate[normalization Tensor[Bar @ coupling] @ ## ] &
		];

		(*Adds the quartic coupling to the association*)
		AppendTo[$quartics, coupling ->
			<|Coupling -> lam,
			CouplingBar -> lambar,
			Fields -> {phi1, phi2, phi3, phi4},
			Indices -> OptionValue @ CouplingIndices,
			Invariant -> OptionValue @ GroupInvariant,
			SelfConjugate -> OptionValue @ SelfConjugate,
			UniqueArrangements -> CouplingPermutations[{phi1, phi2, phi3, phi4}, OptionValue @ GroupInvariant]|>];

		(*Constructs the projection operator*)
		projection = Evaluate[ ReplaceAll[tGen[rep_, A_, a_, b_] -> tGen[Bar@rep, A, a, b]] @ OptionValue[GroupInvariant][s1, s2, s3, s4] *
			sDelS[Bar@phi1, #1, s1] sDelS[Bar@phi2, #2, s2] sDelS[Bar@phi3, #3, s3] sDelS[Bar@phi4, #4, s4]] &;
		symmetryFactor = ReplaceAll[Rule[#, 0] & /@ DeleteCases[Keys @ $quartics, coupling]] @ RefineGroupStructures @
			Expand[projection[$a, $b, $c, $d] Lam[$a, $b, $c, $d] ] /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 // Simplify;
		If[symmetryFactor === 0,
			Message[AddQuartic::projection0];
			KeyDropFrom[$quartics, coupling];
			Return @ $Failed;
		];
		projection = Evaluate[projection[#1, #2, #3, #4] / symmetryFactor] &;
		AppendTo[$quartics @ coupling, Projector -> projection];

		(*Tests whether the coupling satisfy gauge invariance*)
		If[OptionValue @ CheckInvariance,
			Block[{A,a1,a2,a3,a4,a5},
				lambda[$da_, $db_, $dc_, $dd_] = sDelS[phi1, $da, s1$1] sDelS[phi2, $db, s2$1] sDelS[phi3, $dc, s3$1] sDelS[phi4, $dd, s4$1] OptionValue[GroupInvariant][s1$1, s2$1, s3$1, s4$1];
				test = lambda[a5,a2,a3,a4] Tscal[A, a5, a1] + lambda[a1,a5,a3,a4] Tscal[A, a5, a2] + lambda[a1,a2,a5,a4] Tscal[A, a5, a3] + lambda[a1,a2,a3,a5] Tscal[A, a5, a4] //Expand;
				(* test = Tscal[A, a1, a5] Lam[a5,a2,a3,a4] //Expand //Sym4[a1,a2,a3,a4]; *)
				test *= sDelS[Bar@phi1, a1, s1] sDelS[Bar@phi2, a2, s2] sDelS[Bar@phi3, a3, s3] sDelS[Bar@phi4, a4, s4] //Expand;
				Do[
					temp = CanonizeTensor[{"s1$*","s2$*","s3$*","s4$*","f1$*","f2$*"}][ test sDelV[$gaugeGroups[group, Field], A, v1] //Expand ] /. Matrix[x_][_[v1]] -> x;
					If [temp =!= 0 && !MatchQ[temp, {0 ..}],
						Message[AddQuartic::gaugeInv, group];
						Print["0 != ", temp];
					];
				,{group, Keys@ $gaugeGroups}];
			];
		];

		AppendTo[$couplings, coupling -> Quartic];
		ResetBetas[];
	];


(*Function for defining the trilinear scalar couplings of the theory*)
AddTrilinear::unkown = "`1` does not match any of the scalars.";
AddTrilinear::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors.";
Options[AddTrilinear] = {CouplingIndices -> (Null &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True};
AddTrilinear [coupling_, {phi1_, phi2_, phi3_}, OptionsPattern[]] ? OptionsCheck :=
	Block[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
		(*Tests if the fields have been defined*)
		Do[
			If[!MemberQ[Keys@ $scalars, temp /. Bar[x_] -> x],
				Message[AddQuartic::unkown, temp /. Bar[x_] -> x];
				Return[Null];
			];
		,{temp, {phi1, phi2, phi3}}];

		If[OptionValue @ SelfConjugate,
			SetReal @ coupling;
		];

		normalization = 24; (*To account for the structure deltas working with real/complex scalars in Lam[a, b, c, d]*)
		Do[
			If[!$scalars[phi, SelfConjugate] || Head @ phi === Bar, normalization *= 1/Sqrt[2];];
		,{phi, {phi1, phi2, phi3}}];

		If[Length @ OptionValue[CouplingIndices][s1, s2, s3] <= 2,
			lam = Evaluate[normalization Matrix[coupling] @ ##] &;
			lambar = Evaluate[normalization Matrix[Bar @ coupling] @ ## ] &
		,
			lam = Evaluate[normalization Tensor[coupling] @ ##] &;
			lambar = Evaluate[normalization Tensor[Bar @ coupling] @ ## ] &
		];

		(*Adds the quartic coupling to the association*)
		AppendTo[$trilinears, coupling ->
			<|Coupling -> lam,
			CouplingBar -> lambar,
			Fields -> {phi1, phi2, phi3},
			Indices -> OptionValue @ CouplingIndices,
			Invariant -> OptionValue @ GroupInvariant,
			SelfConjugate -> OptionValue @ SelfConjugate,
			UniqueArrangements -> CouplingPermutations[{phi1, phi2, phi3, $vev}, OptionValue @ GroupInvariant]|>];

		(*Constructs the projection operator*)
		projection = Evaluate[ ReplaceAll[tGen[rep_, A_, a_, b_] -> tGen[Bar@rep, A, a, b]] @ OptionValue[GroupInvariant][s1, s2, s3] *
			sDelS[Bar@phi1, #1, s1] sDelS[Bar@phi2, #2, s2] sDelS[Bar@phi3, #3, s3] sDelS[$vevSelect, #4, s4]] &;
		symmetryFactor = ReplaceAll[Rule[#, 0] & /@ DeleteCases[Keys @ $trilinears, coupling]] @ RefineGroupStructures @
			Expand[projection[$a, $b, $c, $d] Lam[$a, $b, $c, $d, True] ] /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 // Simplify;
		If[symmetryFactor === 0,
			Message[AddTrilinear::projection0];
			KeyDropFrom[$trilinears, coupling];
			Return @ $Failed;
		];
		projection = Evaluate[projection[#1, #2, #3, #4] / symmetryFactor] &;
		AppendTo[$trilinears @ coupling, Projector -> projection];
		AppendTo[$couplings, coupling -> Trilinear];

		ResetBetas[];
	];

(*Function for defining the Scalar mass terms of the theory*)
AddScalarMass::unkown = "`1` does not match any of the scalars.";
AddScalarMass::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors.";
Options[AddScalarMass] = {MassIndices -> (Null &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True};
AddScalarMass [coupling_, {phi1_, phi2_}, OptionsPattern[]] ? OptionsCheck:=
	Block[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
		(*Tests if the fields have been defined*)
		Do[
			If[!MemberQ[Keys@ $scalars, temp /. Bar[x_] -> x],
				Message[AddQuartic::unkown, temp /. Bar[x_] -> x];
				Return[Null];
			];
		,{temp, {phi1, phi2}}];

		If[OptionValue @ SelfConjugate,
			SetReal @ coupling;
		];

		normalization = 24; (*To account for the structure deltas working with real/complex scalars in Lam[a, b, c, d]*)
		Do[
			If[!$scalars[phi, SelfConjugate] || Head @ phi === Bar, normalization *= 1/Sqrt[2];];
		,{phi, {phi1, phi2}}];

		If[Length @ OptionValue[MassIndices][s1, s2] <= 2,
			lam = Evaluate[normalization Matrix[coupling] @ ##] &;
			lambar = Evaluate[normalization Matrix[Bar @ coupling] @ ## ] &
		,
			lam = Evaluate[normalization Tensor[coupling] @ ##] &;
			lambar = Evaluate[normalization Tensor[Bar @ coupling] @ ## ] &
		];

		(*Adds the quartic coupling to the association*)
		AppendTo[$scalarMasses, coupling ->
			<|Coupling -> lam,
			CouplingBar -> lambar,
			Fields -> {phi1, phi2},
			Indices -> OptionValue @ MassIndices,
			Invariant -> OptionValue @ GroupInvariant,
			SelfConjugate -> OptionValue @ SelfConjugate,
			UniqueArrangements -> CouplingPermutations[{phi1, phi2, $vev, $vev}, OptionValue @ GroupInvariant]|>];

		(*Constructs the projection operator*)
		projection = Evaluate[ ReplaceAll[tGen[rep_, A_, a_, b_] -> tGen[Bar@rep, A, a, b]] @ OptionValue[GroupInvariant][s1, s2] *
			sDelS[Bar@phi1, #1, s1] sDelS[Bar@phi2, #2, s2] sDelS[$vevSelect, #3, s3] sDelS[$vevSelect, #4, s4]]& ;
		symmetryFactor = ReplaceAll[Rule[#, 0] & /@ DeleteCases[Keys @ $scalarMasses, coupling]] @ RefineGroupStructures @
			Expand[projection[$a, $b, $c, $d] Lam[$a, $b, $c, $d, True] ] /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 // Simplify;
		If[symmetryFactor === 0,
			Message[AddScalarMass::projection0];
			KeyDropFrom[$scalarMasses, coupling];
			Return @ $Failed;
		];
		projection = Evaluate[projection[#1, #2, #3, #4] / symmetryFactor] &;
		AppendTo[$scalarMasses @ coupling, Projector -> projection];
		AppendTo[$couplings, coupling -> ScalarMass];

		ResetBetas[];
	];

(*Genral quartic coupling used in the computation of the beta function tensors.*)
Lam[$da_, $db_, $dc_, $dd_, massive_:False] :=
	Module[{l, s1, s2, s3, s4},
		(*Quartic couplings*)
		Sum[
			AveragePermutations[{$da, $db, $dc, $dd}, l[UniqueArrangements] ][
				sDelS[l[Fields][[1]], $da, s1] sDelS[l[Fields][[2]], $db, s2] sDelS[l[Fields][[3]], $dc, s3] sDelS[l[Fields][[4]], $dd, s4]
				* l[Invariant][s1, s2, s3, s4] l[Coupling][Sequence @@ l[Indices][s1, s2, s3, s4]]
			+ If[! l@SelfConjugate,
				sDelS[Bar@ l[Fields][[1]], $da, s1] sDelS[Bar@ l[Fields][[2]], $db, s2] sDelS[Bar@ l[Fields][[3]], $dc, s3]
				* sDelS[Bar@ l[Fields][[4]], $dd, s4] l[Invariant][s1, s2, s3, s4] l[CouplingBar][Sequence @@ l[Indices][s1, s2, s3, s4]]
			,0]
			]
		,{l, $quartics}] +
		If[massive,
			(*Trilinear couplings*)
			Sum[AveragePermutations[{$da, $db, $dc, $dd}, l[UniqueArrangements] ][
					sDelS[l[Fields][[1]], $da, s1] sDelS[l[Fields][[2]], $db, s2] sDelS[l[Fields][[3]], $dc, s3] sDelS[$vev, $dd, s4]
					* l[Invariant][s1, s2, s3] l[Coupling][Sequence @@ l[Indices][s1, s2, s3]]
				+ If[! l@SelfConjugate,
					sDelS[Bar@ l[Fields][[1]], $da, s1] sDelS[Bar@ l[Fields][[2]], $db, s2] sDelS[Bar@ l[Fields][[3]], $dc, s3]
					* sDelS[$vev, $dd, s4] l[Invariant][s1, s2, s3] l[CouplingBar][Sequence @@ l[Indices][s1, s2, s3]]
				,0]
				]
			,{l, $trilinears}] +
			(*Scalar masses*)
			Sum[AveragePermutations[{$da, $db, $dc, $dd}, l[UniqueArrangements] ][
					sDelS[l[Fields][[1]], $da, s1] sDelS[l[Fields][[2]], $db, s2] sDelS[$vev, $dc, s3] sDelS[$vev, $dd, s4]
					* l[Invariant][s1, s2] l[Coupling][Sequence @@ l[Indices][s1, s2]]
				+ If[! l@SelfConjugate,
					sDelS[Bar@ l[Fields][[1]], $da, s1] sDelS[Bar@ l[Fields][[2]], $db, s2] sDelS[$vev, $dc, s3] sDelS[$vev, $dd, s4]
					* l[Invariant][s1, s2] l[CouplingBar][Sequence @@ l[Indices][s1, s2]]
				,0]
				]
			,{l, $scalarMasses}]
		,0]
	];



(*######################################*)
(*---------------Clean up---------------*)
(*######################################*)
(*Removes all model information.*)
ResetModel[] :=
	Block[{},
		(*Resets tensor dummy notation*)
		ReInitializeSymbols[];

		(*Resets all information of the current model.*)
		(*Global couplings variable*)
		$couplings = <||>;
		(*Association with all information on the gauge groups: fields, couplings etc.*)
		$gaugeGroups = <||>;
		(*Associationwith all information on the fermion fields: representations etc.*)
		$fermions = <||>;
		(*Associationwith all information on the scalar fields: representations, etc.*)
		$scalars = <||>;
		(*Associationwith all information on the quartic couplings.*)
		$quartics = <||>;
		(*Associationwith all information on the Yukawa couplings.*)
		$yukawas = <||>;
		(*Associationwith all information on the fermion masses.*)
		$fermionMasses = <||>;
		(*Associationwith all information on the trilinear scalar couplings.*)
		$trilinears = <||>;
		(*Associationwith all information on the scalar masses.*)
		$scalarMasses = <||>;

		(*Removes stored computations of the beta tensors.*)
		ResetBetas[];
	];

(*Function for removing a scalar or fermion field (and associated Yukawa and quartic couplings) from the model.*)
RemoveField::unkown = "The field `1` has not been defined."
RemoveField[field_] :=
	Module[{coupling},
		Switch[field
		,f_ /; MemberQ[Keys @ $fermions, f],
			$fermions = Delete[$fermions, Key @ field];
			(*Remove yukawa interactions the fermion is involved in*)
			Do[
				If[MemberQ[$yukawas[coupling, Fields][[2;;3]], field],
					RemoveInteraction[coupling];
				];
			,{coupling, Keys @ $yukawas}];
		,f_ /; MemberQ[Keys @ $scalars, f],
			If[$scalars[field, SelfConjugate],
				Bar[field] =.;
			];

			$scalars = Delete[$scalars, Key @ field];
			(*Remove yukawa interactions the scalar is involved in*)
			Do[
				If[MemberQ[{#, Bar @ #}& @ $yukawas[coupling, Fields][[1]], field],
					RemoveInteraction[coupling];
				];
			,{coupling, Keys @ $yukawas}];
			(*Remove quartic interactions the scalar is involved in*)
			Do[
				If[MemberQ[$quartics[coupling, Fields], field] || MemberQ[Bar /@ $quartics[coupling, Fields], field],
					RemoveInteraction[coupling];
				];
			,{coupling, Keys @ $quartics}];
		,_,
			Message[RemoveField::unkown, field];
			Return @ $Failed;
		];
	];
(*For simmultaneous removal of multiple fields*)
RemoveField[field_, fields__] :=
	Block[{},
		RemoveField[field];
		RemoveField[fields];
	];

(*Function for removing an interaction from the model.*)
RemoveInteraction::unkown = "The coupling `1` has not been defined."
RemoveInteraction[coupling_] :=
	Module[{group, lieG},
		Switch[$couplings @ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			(*For gauge groups the corresponding gauge field is removed.*)
			sDelV /: sDelV[$gaugeGroups[$couplings @ coupling, Field], ind_, v1_] *
				sDelV[$gaugeGroups[$couplings @ coupling, Field], ind_, v2_] =. ;
			$gaugeGroups = Delete[$gaugeGroups, Key @ $couplings @ coupling];
			(*The tensor symbols are reset, and reloaded for all other gauge groups.*)
			ReInitializeSymbols[];
			Do[
				lieG = $gaugeGroups[group, LieGroup];
				Switch[Head @ lieG
				,SO,
					DefineSOGroup[group, lieG[[1]] ];
				,Sp,
					DefineSpGroup[group, lieG[[1]]];
				,SU,
					DefineSUGroup[group, lieG[[1]]];
				,U,
					DefineU1Group[group];
				];
			,{group, Keys @ $gaugeGroups}];
		,Yukawa,
			$yukawas = Delete[$yukawas, Key @ coupling];
		,Quartic,
			$quartics = Delete[$quartics, Key @ coupling];
		,FermionMass,
			$fermionMasses = Delete[$fermionMasses, Key @ coupling];
		,Trilinear,
			$trilinears = Delete[$trilinears, Key @ coupling];
		,ScalarMass,
			$scalarMasses = Delete[$scalarMasses, Key @ coupling];
		,_Missing,
			Message[RemoveInteraction::unkown, coupling];
			Return @ $Failed;
		];
		$couplings = Delete[$couplings, Key @ coupling];
		ResetBetas[];
	];
(*For simmultaneous removal of multiple interactions*)
RemoveInteraction[coupling_, couplings__] :=
	Block[{},
		RemoveInteraction[coupling];
		RemoveInteraction[couplings];
	];



End[]
