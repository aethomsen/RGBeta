(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageExport["AddFermion"]
PackageExport["AddFermionMass"]
PackageExport["AddGaugeGroup"]
PackageExport["AddQuartic"]
PackageExport["AddScalar"]
PackageExport["AddScalarMass"]
PackageExport["AddTrilinear"]
PackageExport["AddVector"]
PackageExport["AddYukawa"]
PackageExport["SetSymmetries"]

PackageExport["ResetModel"]

PackageExport["$couplings"]
PackageExport["$fermionMasses"]
PackageExport["$fermions"]
PackageExport["$gaugeGroups"]
PackageExport["$quartics"]
PackageExport["$scalarMasses"]
PackageExport["$scalars"]
PackageExport["$trilinears"]
PackageExport["$yukawas"]

PackageScope["FGauge"]
PackageScope["G2Matrix"]
PackageScope["Lam"]
PackageScope["Tferm"]
PackageScope["TfermTil"]
PackageScope["Tscal"]
PackageScope["Yuk"]
PackageScope["YukTil"]

PackageScope["$fieldIndexMap"]
PackageScope["$flavorReps"]

PackageScope["AveragePermutations"]
PackageScope["CouplingPermutations"]
PackageScope["GatherToList"]
PackageScope["RepresentationCheck"]

PackageScope["UpdateProjectors"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

AddFermion::usage =
	"AddFermion[field] is a function used to define a fermion field in the model."
AddFermionMass::ussage =
	"AddFermionMass[mass, {ferm1, ferm2}] defines a mass term between the two fermion fields."
AddGaugeGroup::usage =
	"AddGaugeGroup[coupling, groupName, lieGroup[n]] is a function used to define a gauge group in the model."
AddQuartic::usage =
	"AddQuartic[coupling, {scal1, scal2, scal3, scal4}] defines a quartic interaction between 4 scalar fields."
AddScalar::usage =
	"AddScalar[field] is a function used to define a scalar field in the model."
AddScalarMass::ussage =
	"AddScalarMass[mass, {scal1, scal2}] defines a mass term between the two scalar fields."
AddTrilinear::usage =
	"AddTrilinear[coupling, {scal1, scal2, scal3}] defines a Yukawa interaction between a scalar and two fermion fields."
AddVector::usage =
	"AddVector[field, group] is an internal function used to define a vector field of a gauge group."
AddYukawa::usage =
	"AddYukawa[coupling, {scal, ferm1, ferm2}] defines a Yukawa interaction between a scalar and two fermion fields."
SetSymmetries::usage =
	"SetSymmetries[coupling, symmetries, antisymmetries] defines (anti)symmetric indices of a coupling. Format is e.g. {{1,4}, {2, ..}}, where the numbers refer to the coupling indices."

RemoveField::usage =
	"RemoveField[field, ...] is a function that removes one or more scalar and fermion fields from the model."
RemoveInteraction::usage =
	"RemoveCoupling[coupling, ...] is a function that removes one or more couplings and corresponding interactions/gauge groups from the model."
ResetModel::usage =
	"ResetModel[] resets all tensors and removes all fields and coupling definitions made in the current instance of RGBeta."

$couplings::usage =
	"$couplings is an association containing the types of all the couplings that have been defined in the model."
$fermionMasses::usage =
	"$fermionMasses is an association with all the fermion masses that have been declared in the model."
$fermions::usage =
	"$fermions is an association containing all the internal information on the fermions that have been declared in the model."
$gaugeGroups::usage =
	"$gaugeGroups is an association containing all the internal information on the gauge groups that have been declared in the model."
$quartics::usage =
	"$quartics is an association containing all internal information on the quartic couplings."
$scalarMasses::usage =
	"$scalarMasses is an association with all the scalar masses that have been declared in the model."
$scalars::usage =
	"$scalars is an association containing all the internal information on the scalars that have been declared in the model."
$trilinears::usage =
	"$trilinears is an association contraining all the internal informaiton on the trilinear scalar couplings."
$yukawas::usage =
	"$yukawas is an association containing all internal information on the yukawa couplings."


FGauge::usage =
	"FGauge[A, B, C] is a function that generates the general gauge structure constants."
G2Matrix::usage =
	"G2Matrix[a, b] is a function that generates the general gauge coupling matrix."
Lam::usage =
	"Lam[a, b, c, d] is a function that generates the general quartic coupling."
Tferm::usage =
	"Tferm[A, i, j] is a function that generates the general gauge generator matrix for the fermions."
TfermTil::usage =
	"TfermTil[A, i, j] is a function that generates the general gauge generator matrix for the fermions with an applied tilde."
Tscal::usage =
	"Tscal[A, i, j] is a function that generates the general gauge generator matrix for the scalars."
Yuk::usage =
	"Yuk[a, i, j] is a function that generates the general yukawa coupling."
YukTil::usage =
	"YukTil[a, i, j] is a function that generates the general yukawa coupling with an applied tilde."

AveragePermutations::usage =
	"AveragePermutations[indices, permutations][expr] gives the average of expr with indices (List) exchanged in all way perscribed by permutations."
CouplingPermutations::usage =
	"CouplingPermutations[fields, invariant] provides at list of all unique permutations of a coupling of the fields (List) and with the given group invariance function."
RepresentationCheck::usage =
	"RepresentationCheck[rep] returns True if rep is a valid representation given the gauge groups of the present model and False otherwise."

(* Function that updates the mapping of fields to indices in the $fieldIndexMap *)
UpdateFieldIndexMap[] :=
	Module[{count, temp, field},
		$fieldIndexMap = <||>;
		count = 1;
		temp = <||>;
		Do[
			AppendTo[temp, field -> count++];
			If[ ! $scalars[field, SelfConjugate], AppendTo[temp, Bar@ field -> count++]; ];
		, {field, Keys@ $scalars}];
		AppendTo[temp, $vev -> count++];	AppendTo[temp, $vevSelect -> count++];
		AppendTo[$fieldIndexMap, "Scalars" -> Association@ temp];

		$scalarContraction = SparseArray[ Table[ Switch[field,
				$vev, temp /@ {$vev, $vevSelect} -> 1,
				$vevSelect, temp /@ {$vevSelect, $vev} -> 1,
				_, temp /@ {field, Bar@ field} -> If[field === Bar@ field, 1, 2] ],
			{field, Keys@ temp}], Length@ temp {1, 1}];

		count = 1;
		AppendTo[$fieldIndexMap,
			"Fermions" -> Association@ Table[field -> count++, {field, Keys@ $fermions} ]];

		count = 1;
		AppendTo[$fieldIndexMap, "GaugeBosons" -> Association@ Table[group -> count++, {group, Keys@ $gaugeGroups} ]];

		(* The association with all fields *)
		AppendTo[$fieldIndexMap, "All" -> Join @@ $fieldIndexMap /@ {"Scalars", "Fermions", "GaugeBosons"}];
	];

(* Constructs projection operators for the various types of couplings*)
UpdateProjectors::projection0 = "The projcetion operator does not pick out the `1` coupling. Please check the GroupInvariant of the coupling for errors and/or set the symmetries of the couplings with the SymmetricIndices/AntisymmetricIndices options."
UpdateProjectors[couplings_List] := Scan[UpdateProjectors, couplings];
UpdateProjectors[coupling_] :=
	Module[{couplingInfo, dim, group, projector, symmetryFactor},
		Switch[$couplings@ coupling
		,x_ /; MemberQ[Keys @ $gaugeGroups, x],
			dim = Length@ $fieldIndexMap["GaugeBosons"] {1, 1};
			group = $couplings@ coupling;
			(* As of Mathematica v13 SparseArray behaves atomically w.r.t. the Slots and has to be made Inactive *)
			AppendTo[$gaugeGroups@ group, Projector ->
				Activate@* (Evaluate[ TStructure[$gauge@ #1, $gauge@ #2]@
				Inactive[SparseArray][$fieldIndexMap["GaugeBosons"] /@ (group {1, 1}) ->
					If[Head@ $gaugeGroups[group, LieGroup] =!= U1,
					del[group@ adj, #1, #2]/Dim@ group@ adj,
					del[group@ adj, #1, v1] del[group@ adj, #2, v2] ],
				dim] ] &)
			];

		,Yukawa,
			dim = {Length@ $fieldIndexMap["Scalars"], Length@ $fieldIndexMap["Fermions"], Length@ $fieldIndexMap["Fermions"]};
			couplingInfo = $yukawas@ coupling;
			(* Constructs the proto-projector *)
			projector = Evaluate[ DiagonalMatrix[ TStructure[$scalar@ #1, $fermion@ #2, $fermion@ #3] /@ Switch[couplingInfo@ Chirality
					,Left,
						{Inactive[SparseArray][$fieldIndexMap["All"] /@
							MapAt[Bar, couplingInfo@ Fields, 1] -> GroupInvBar@ couplingInfo[Invariant][#1, #2, #3]
						, dim], 0}
					,Right,
						{0, Inactive[SparseArray][$fieldIndexMap["All"] /@
							couplingInfo@ Fields -> couplingInfo[Invariant][#1, #2, #3]
						, dim]}
					] ] ] &;
			(* Determines the correct symmetry factor *)
			symmetryFactor = ReplaceAll[Rule[#, 0] & /@ Keys @ $yukawas]@ Simplify@ RefineGroupStructures[
				Tr[Activate@ projector[$a, $i, $j] . Yuk[$a, $i, $j] ]  /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1] // Simplify;
			If[symmetryFactor === 0,
				Message[UpdateProjectors::projection0, coupling];
				KeyDropFrom[$yukawas, coupling];
				Return@ $Failed
			];
			projector = Activate@* (Evaluate[projector[#1, #2, #3] / symmetryFactor] &);
			AppendTo[$yukawas@ coupling, Projector -> projector];

		,Quartic,
			dim = Length@ $fieldIndexMap["Scalars"] {1, 1, 1, 1};
			couplingInfo = $quartics@ coupling;
			(* Constructs the proto-projector *)
			projector = Function[{$da, $db, $dc, $dd}, Evaluate[
				TStructure[$scalar@ $da, $scalar@ $db, $scalar@ $dc, $scalar@ $dd] @ Inactive[SparseArray][
				$fieldIndexMap["Scalars"] /@ Bar /@ couplingInfo@ Fields ->
				GroupInvBar@ couplingInfo[Invariant][$da, $db, $dc, $dd], dim] ] ];
			(* Determines the correct symmetry factor *)
			symmetryFactor = ReplaceAll[Rule[#, 0] & /@ Keys @ $quartics] @ RefineGroupStructures[
				Activate@ projector[$a, $b, $c, $d] Lam[$a, $b, $c, $d]  /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 ] // Simplify;
			If[symmetryFactor === 0,
				Message[UpdateProjectors::projection0, coupling];
				KeyDropFrom[$quartics, coupling];
				Return@ $Failed
			];
			projector = Activate@* (Evaluate[projector[#1, #2, #3, #4] / symmetryFactor] &);
			AppendTo[$quartics@ coupling, Projector -> projector];

		,FermionMass,
			dim = {Length@ $fieldIndexMap["Scalars"], Length@ $fieldIndexMap["Fermions"], Length@ $fieldIndexMap["Fermions"]};
			couplingInfo = $fermionMasses@ coupling;
			(* Constructs the proto-projector *)
			projector = Function[{$da, $di, $dj}, Evaluate[ DiagonalMatrix[
				TStructure[$scalar@ $da, $fermion@ $di, $fermion@ $dj] /@ Switch[couplingInfo@ Chirality
					,Left,
						{Inactive[SparseArray][$fieldIndexMap["All"] /@ Join[{$vevSelect}, couplingInfo@ Fields] ->
							GroupInvBar@ couplingInfo[Invariant][$di, $dj], dim], 0}
					,Right,
						{0, Inactive[SparseArray][$fieldIndexMap["All"] /@ Join[{$vevSelect}, couplingInfo@ Fields] ->
							couplingInfo[Invariant][$di, $dj], dim]}
				] ] ] ];
			(* Determines the correct symmetry factor *)
			symmetryFactor = ReplaceAll[Rule[#, 0] & /@ Keys @ $fermionMasses] @ RefineGroupStructures @
				Tr[Activate@ projector[$a, $i, $j] . Yuk[$a, $i, $j, True]  /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 ]// Simplify;
			If[symmetryFactor === 0,
				Message[UpdateProjectors::projection0, coupling];
				KeyDropFrom[$fermionMasses, coupling];
				Return@ $Failed
			];
			projector = Activate@* (Evaluate[projector[#1, #2, #3] / symmetryFactor] &);
			AppendTo[$fermionMasses@ coupling, Projector -> projector];

		,Trilinear,
			dim = Length@ $fieldIndexMap["Scalars"] {1, 1, 1, 1};
			couplingInfo = $trilinears@ coupling;
			(* Constructs the proto-projector *)
			projector = Function[{$da, $db, $dc, $dd}, Evaluate[
				TStructure[$scalar@ $da, $scalar@ $db, $scalar@ $dc, $scalar@ $dd] @ Inactive[SparseArray][
				$fieldIndexMap["Scalars"] /@ Join[Bar /@ couplingInfo@ Fields, {$vevSelect}] ->
				GroupInvBar@ couplingInfo[Invariant][$da, $db, $dc], dim] ] ];
			(* Determines the correct symmetry factor *)
			symmetryFactor = ReplaceAll[Rule[#, 0] & /@ Keys @ $trilinears] @ RefineGroupStructures[
				Activate@ projector[$a, $b, $c, $d] Lam[$a, $b, $c, $d, True]  /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 ] // Simplify;
			If[symmetryFactor === 0,
				Message[UpdateProjectors::projection0, coupling];
				KeyDropFrom[$trilinears, coupling];
				Return@ $Failed
			];
			projector = Activate@* (Evaluate[projector[#1, #2, #3, #4] / symmetryFactor] &);
			AppendTo[$trilinears@ coupling, Projector -> projector];

		,ScalarMass,
			dim = Length@ $fieldIndexMap["Scalars"] {1, 1, 1, 1};
			couplingInfo = $scalarMasses@ coupling;
			(* Constructs the proto-projector *)
			projector = Function[{$da, $db, $dc, $dd}, Evaluate[
				TStructure[$scalar@ $da, $scalar@ $db, $scalar@ $dc, $scalar@ $dd] @ Inactive[SparseArray][
				$fieldIndexMap["Scalars"] /@ Join[Bar /@ couplingInfo@ Fields, {$vevSelect, $vevSelect}] ->
				GroupInvBar@ couplingInfo[Invariant][$da, $db], dim] ] ];
			(* Determines the correct symmetry factor *)
			symmetryFactor = ReplaceAll[Rule[#, 0] & /@ Keys @ $scalarMasses] @ RefineGroupStructures[
				Activate@ projector[$a, $b, $c, $d] Lam[$a, $b, $c, $d, True]  /. (Matrix|Tensor)[x_][__] -> x /. coupling -> 1 ] // Simplify;
			If[symmetryFactor === 0,
				Message[UpdateProjectors::projection0, coupling];
				KeyDropFrom[$scalarMasses, coupling];
				Return@ $Failed
			];
			projector = Activate@* (Evaluate[projector[#1, #2, #3, #4] / symmetryFactor] &);
			AppendTo[$scalarMasses@ coupling, Projector -> projector];
		];
	];

(*Initiates a scalar field*)
Options[AddScalar] = {SelfConjugate -> False, GaugeRep -> {}, FlavorIndices -> {}, Mass -> None};
AddScalar::failure = "Failed to add scalar field.";
AddScalar[field_, OptionsPattern[] ] ? OptionsCheck :=
	Module[{rep, projector, scal, dim},
		(*Checks gauge representations*)
		If[! And @@ Table[RepresentationCheck @ rep, {rep, OptionValue[GaugeRep]}],
			Message[AddScalar::failure];
			Return @ $Failed;
		];

		If[OptionValue@ SelfConjugate,
			SetReal @ field;
		];

		(*Adds field options to the list of scalar fields*)
		AppendTo[$scalars, field -> <|
			GaugeRep -> OptionValue[GaugeRep],
			FlavorIndices -> OptionValue[FlavorIndices],
			SelfConjugate -> OptionValue[SelfConjugate]|>];
		UpdateFieldIndexMap[];

		(* Update flavor reps list *)
		$flavorReps = Union[$flavorReps, OptionValue[FlavorIndices]];

		(* The projectors must be updated to accomodate the changes in the $fieldIndexMap *)
		UpdateProjectors[Join[Keys@ $yukawas, Keys@ $quartics, Keys@ $fermionMasses, Keys@ $trilinears, Keys@ $scalarMasses] ];
		ResetBetas[];
	];

(*Initiates a fermion field*)
Options[AddFermion] = {GaugeRep -> {}, FlavorIndices -> {}};
AddFermion::failure = "Failed to add fermion field.";
AddFermion[field_, OptionsPattern[] ] ? OptionsCheck :=
	Module[{rep, projector, ferm, dim},
		(*Checks gauge representations*)
		If[! And @@ Table[RepresentationCheck @ rep, {rep, OptionValue[GaugeRep]}],
			Message[AddFermion::failure];
			Return @ $Failed;
		];

		(*Adds field options to the list of fermion fields*)
		AppendTo[$fermions, field -> <|
			GaugeRep -> OptionValue[GaugeRep],
			FlavorIndices -> OptionValue[FlavorIndices]
			|>];
		UpdateFieldIndexMap[];

		(* Update flavor reps list *)
		$flavorReps = Union[$flavorReps, OptionValue[FlavorIndices]];

		(* The projectors must be updated to accomodate the changes in the $fieldIndexMap *)
		UpdateProjectors[Join[Keys@ $yukawas, Keys@ $fermionMasses] ];
		ResetBetas[];
	];

(*###################################*)
(*----------Other functions----------*)
(*###################################*)
(* Determines all unique permutations of fields with a given invariant *)
CouplingPermutations[fields_List, invariant_Function] :=
	Module[{arg, fs, inv, number, perms, uniformSymbols, uniqueArrangements},
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
					delS2[rep_, i_, a__] :> delS2[rep, i, Sequence @@ Sort @ List @ a],
					dSym[gr_, inds__]:> dSym[gr, Sequence @@ Sort @ List @ inds]};
				If[MemberQ[uniqueArrangements, {fs, inv}],
					Continue[];
				];
				AppendTo[uniqueArrangements, {fs, inv}];
				Sow @ p;
			,{p, perms}];
		][[2,1]]
	];

AveragePermutations[indices_List, permutations_List][expr_] :=
	Module[{subs},
		subs = MapThread[Rule, {indices, indices[[#]]}] & /@ permutations;
		Mean[expr/.subs]
	];

(* Puts a tensor on canonical form---Brute force *)
CanonizeTensor[symbForm_List][expr_] :=
	Module[{symb, newSymbs, indexAssignments, temp},
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

(* Collects a list with overlapping entires to make it suitbale for SparseArray *)
GatherToList[expr_] := Plus @@@ GroupBy[Flatten@ expr, (#[[1]] &) -> (#[[2]] &)] // Normal;

(* Function for complex conjugating an expression made up of the group invariants defined in GroupsAndIndices *)
GroupInvBar[expr_] := ReplaceAll[{tGen[rep_, A_, a_, b_] -> tGen[rep, A, b, a], x_Complex-> Conjugate@ x}]@ expr;

(* Function for setting symmetries of couplings and behavior in Matrix or Tensor functions and Trans. *)
SetSymmetries::invalid= "`1` is not a valid format for symmetrized couplings."
SetSymmetries::unkown= "`1` is not a Yukawa or Quartic coupling."
SetSymmetries::inds= "`1` has but `2` indices and no symmetry can be assigned."
SetSymmetries::dupl= "Index appearing multiple times accros the (anti)symmetry lists."
SetSymmetries::outbound= "(Anti)symmetry lists refers to indices of higher rank than the `1` tensor."
SetSymmetries[coupling_Symbol, {}, {}]:= Null;
SetSymmetries[coupling_Symbol, symInds_List, asymInds_List]:= Module[{couplingInds, allInds, set},
	{sym, asym}= PrepareSymList/@ {symInds, asymInds};

	Switch[$couplings@ coupling
		,Yukawa,
			couplingInds= Length@ $yukawas[coupling, Indices][$a, $i, $j];
		,Quartic,
			couplingInds= Length@ $quartics[coupling, Indices][$a, $b, $c, $d];
		,_,
			Message[SetSymmetries::unkown, coupling];
			Abort[];
	];
	If[couplingInds <= 1,
		Message[SetSymmetries::inds, coupling, couplingInds];
		Abort[];
	];

	allInds= sym~Join~asym;
	If[allInds=== {}, Return@ Null;];
	If[!DuplicateFreeQ@ allInds,
		Message[SetSymmetries::dupl];
		Abort[];
	];
	If[Max@ allInds > couplingInds,
		Message[SetSymmetries::outbound, coupling];
		Abort[];
	];

	(* If coupling is a matrix *)
	If[couplingInds === 2,
		If[sym=== {{1,2}},
			Matrix[coupling][i_, j_] /; ! OrderedQ@ {i, j} := Matrix[coupling][j, i];
			Trans@ coupling := coupling;
			Trans@ Bar@ coupling := Bar@ coupling;
		];
		If[asym=== {{1,2}},
			Matrix[coupling][i_, j_] /; ! OrderedQ@ {i, j} := -Matrix[coupling][j, i];
			Trans@ coupling := -coupling;
			Trans@ Bar@ coupling := -Bar@ coupling;
		];
	];

	(* If coupling is a higher rank tensor *)
	If[couplingInds > 2,
		Do[
			With[{setHere= set},
				Tensor[coupling][inds__] /; !OrderedQ@ {inds}[[setHere]]:= Tensor[coupling]@@ SubSort[{inds}, setHere];
			]
		, {set, sym}];
		Do[
			With[{setHere= set},
				Tensor[coupling][inds__] /; !OrderedQ@ {inds}[[setHere]]:= Signature@ {inds}[[setHere]] Tensor[coupling]@@ SubSort[{inds}, setHere];
			];
		, {set, asym}];
	];

];
PrepareSymList[inds_List]:= Module[{},
	If[MatchQ[inds, {_Integer..}], Return@ {inds}; ];
	If[inds=== {} || inds=== {{}}, Return@{}; ];
	If[!MatchQ[inds, {{_Integer..}..}],
		Message[SetSymmetries::invalid, inds];
		Abort[];
	];
	Sort/@ DeleteCases[inds, (Length@ # <2 &)]
];
SubSort[list_List, set_List]:= Block[{out=list},
	out[[set]]= Sort@ out[[set]];
	out
];

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
	Module[{cMatrix, invalid, projector, fieldName},
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
		(* AddVector[fieldName, groupName];
		projector = Evaluate[If[lieGroup =!= U1, del[groupName[adj], v1, v2] / Dim @ groupName[adj], 1] *
			sDelV[fieldName, #1, v1] sDelV[fieldName, #2, v2] ] &; *)

		(*Adds the group information and the coupling to the repsective lists*)
		AppendTo[$gaugeGroups, groupName ->
			<|Coupling -> coupling,
			Field -> fieldName,
			LieGroup -> lieGroup[n]|>];
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

		(* The construction of the projector of the coupling has been relegated to UpdateFieldIndexMap[] *)
		UpdateFieldIndexMap[];
		UpdateProjectors[$gaugeGroups[#, Coupling] & /@ Keys@ $gaugeGroups];
		ResetBetas[];
	];

(*Checks the validity of a given representation. Returns True iff rep is valid given gauge groups
 of the present model*)
RepresentationCheck::invalid = "`1` is not a reckognized format for a representation.";
RepresentationCheck::representation = "`1` is not a reckognized representation of a `2` group.";
RepresentationCheck::noGroup = "The `1` gauge group has not been defined.";
RepresentationCheck[rep_] :=
	Module[{group, gName},
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


(*####################################*)
(*----------Yukawa couplings----------*)
(*####################################*)
(*Function for defining the Yukawa couplings of the theory*)
Options[AddYukawa] = {
	CouplingIndices -> ({} &),
	GroupInvariant -> (1 &),
	Chirality -> Left,
	CheckInvariance -> False,
	AntisymmetricIndices -> {},
	SymmetricIndices -> {}
};
AddYukawa::unkown = "`1` does not match any of the `2`s.";
AddYukawa::gaugeInv = "Could not verify invariance under the `1` group. Plase check the GroupInvariant for errors.";
AddYukawa[coupling_, {phi_, psi1_, psi2_}, OptionsPattern[]] ? OptionsCheck:=
	Module[{g, group, normalization, projection, symmetryFactor, temp, test, yuk, yukawa, yukbar, y},
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
			,Left, g = coupling;
			,Right,	g = Bar @ coupling;
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
			Invariant -> OptionValue @ GroupInvariant,
			UniqueArrangements -> Cases[CouplingPermutations[{phi, psi1, psi2}, OptionValue @ GroupInvariant], {1, __}]
			|>];

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
		SetSymmetries[coupling, OptionValue@ SymmetricIndices, OptionValue@ AntisymmetricIndices];
		UpdateProjectors[coupling];
		ResetBetas[];
	];


(*Function for defining the fermion masses of the theory*)
Options[AddFermionMass] = {
	MassIndices -> ({} &),
	GroupInvariant -> (1 &),
	Chirality -> Left,
	AntisymmetricIndices -> {},
	SymmetricIndices -> {}
};
AddFermionMass::unkown = "`1` does not match any of the `2`s.";
AddFermionMass[coupling_, {psi1_, psi2_}, OptionsPattern[]] ? OptionsCheck:=
	Module[{g, group, projection, symmetryFactor, temp, test, yuk, yukbar, y},
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
			,Left, g = coupling;
			,Right, g = Bar @ coupling;
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
		AppendTo[$fermionMasses, coupling ->
			<|Chirality -> OptionValue @ Chirality,
			Coupling -> yuk,
			CouplingBar -> yukbar,
			Fields -> {psi1, psi2},
			Indices -> OptionValue @ MassIndices,
			Invariant -> OptionValue @ GroupInvariant,
			UniqueArrangements -> Prepend[1]/@ (1+ CouplingPermutations[{psi1, psi2}, OptionValue @ GroupInvariant])
			|>];

		AppendTo[$couplings, coupling -> FermionMass];
		SetSymmetries[coupling, OptionValue@ SymmetricIndices, OptionValue@ AntisymmetricIndices];
		UpdateProjectors[coupling];
		ResetBetas[];
	];


(*#####################################*)
(*----------Quartic couplings----------*)
(*#####################################*)

(*Function for defining the quartic couplings of the theory*)
AddQuartic::unkown = "`1` does not match any of the scalars.";
AddQuartic::gaugeInv = "Could not verify invariance under the `1` group. Plase check the GroupInvariant for errors.";
Options[AddQuartic] = {
	CouplingIndices -> ({} &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True,
	CheckInvariance -> False,
	AntisymmetricIndices -> {},
	SymmetricIndices -> {}
};
AddQuartic [coupling_, {phi1_, phi2_, phi3_, phi4_}, OptionsPattern[]] ? OptionsCheck :=
	Module[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
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
		SetSymmetries[coupling, OptionValue@ SymmetricIndices, OptionValue@ AntisymmetricIndices];
		UpdateProjectors[coupling];
		ResetBetas[];
	];


(*Function for defining the trilinear scalar couplings of the theory*)
AddTrilinear::unkown = "`1` does not match any of the scalars.";
AddTrilinear::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors.";
Options[AddTrilinear] = {
	CouplingIndices -> ({} &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True,
	AntisymmetricIndices -> {},
	SymmetricIndices -> {}
};
AddTrilinear [coupling_, {phi1_, phi2_, phi3_}, OptionsPattern[]] ? OptionsCheck :=
	Module[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
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

		AppendTo[$couplings, coupling -> Trilinear];
		SetSymmetries[coupling, OptionValue@ SymmetricIndices, OptionValue@ AntisymmetricIndices];
		UpdateProjectors[coupling];
		ResetBetas[];
	];

(*Function for defining the Scalar mass terms of the theory*)
AddScalarMass::unkown = "`1` does not match any of the scalars.";
AddScalarMass::projection0 = "The projcetion operator does not pick out the coupling. Please check the GroupInvariant for errors.";
Options[AddScalarMass] = {
	MassIndices -> ({} &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True,
	AntisymmetricIndices -> {},
	SymmetricIndices -> {}
};
AddScalarMass [coupling_, {phi1_, phi2_}, OptionsPattern[]] ? OptionsCheck:=
	Module[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
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

		AppendTo[$couplings, coupling -> ScalarMass];
		SetSymmetries[coupling, OptionValue@ SymmetricIndices, OptionValue@ AntisymmetricIndices];
		UpdateProjectors[coupling];
		ResetBetas[];
	];


(*#####################################*)
(*----------Couplings tensors----------*)
(*#####################################*)

(* Function for averaging a pos -> coupling over permutations of the indices *)
AvgPermutations[couplingPos_ -> coupling_, indices_, permutations_] :=
	Module[{subs, len = Length@ permutations, i},
		subs = MapThread[Rule, {indices, indices[[#]] }] & /@ permutations;
		Table[couplingPos[[Ordering@ permutations[[i]] ]] -> coupling/len /. subs[[i]], {i, len}]
	];

(* Generic structure for the coupling tensor creation *)
CouplingTensorMap[info_, fields_, coupInds_, allInds_, permutations_, conj_:False] :=
	AvgPermutations[
		(* If the conjugate coupling is used *)
			If[conj,
				$fieldIndexMap["All"] /@ fields ->
				GroupInvBar[info[Invariant] @@ coupInds] info[CouplingBar] @@ (info[Indices] @@ coupInds)
			,
				$fieldIndexMap["All"] /@ fields ->
				info[Invariant] @@ coupInds info[Coupling] @@ (info[Indices] @@ coupInds)
			]
		, allInds, permutations];

(*The general gauge coupling matrix G^2_{AB} used in the computation of the beta function tensors*)
G2Matrix[$dA_, $dB_] :=
	Module[{temp, group, groupInfo, dim = Length@ $gaugeGroups {1, 1} },
		If[MemberQ[dim, 0], Return@ 0;];
		temp = Table[groupInfo = $gaugeGroups@ group;
				$fieldIndexMap["GaugeBosons"] /@ (group {1, 1}) ->
				If[MatchQ[groupInfo@ LieGroup, U1[n_] /; n > 1],
					Matrix[groupInfo@ Coupling][group[adj]@ $dA, group[adj]@ $dB]
				,
					groupInfo[Coupling]^2 del[group@ adj, $dA, $dB]
				]
			, {group, Keys@$gaugeGroups}];
		TStructure[$gauge@ $dA, $gauge@ $dB]@ SparseArray[GatherToList @ temp, dim]
	];

(*The general fermion gauge generators used in the computation of the beta function tensors*)
Tferm[$dA_, $di_, $dj_] :=
	Module[{tLeft, ferm, gRep, group, otherInds,
		dim = {Length@ $fieldIndexMap["GaugeBosons"], Length@ $fieldIndexMap["Fermions"], Length@ $fieldIndexMap["Fermions"]}},
		If[MemberQ[dim, 0], Return@ {{0, 0}, {0, 0}};];
		tLeft = Evaluate@ Flatten@ Table[
			group = Head@ If[Head@ gRep === Bar, gRep[[1]], gRep];
			otherInds = DeleteCases[$fermions[ferm, GaugeRep], gRep] ~ Join ~ $fermions[ferm, FlavorIndices];
			$fieldIndexMap["All"] /@ {group, ferm, ferm} -> tGen[gRep, $dA, #1, #2] * Product[del[rep, #1, #2], {rep, otherInds}]
		,{ferm, Keys@ $fermions}, {gRep, $fermions[ferm, GaugeRep]}] &;

		TStructure[$gauge@ $dA, $fermion@ $di, $fermion@ $dj] /@ {SparseArray[GatherToList@ tLeft[$di, $dj], dim],
			- SparseArray[GatherToList@ tLeft[$dj, $di], dim]} //DiagonalMatrix
	];

TfermTil[A_, i_, j_] := {{0, 1}, {1, 0}}.Tferm[A, i, j].{{0, 1}, {1, 0}};

(*The general scalar gauge generators used in the computation of the beta function tensors*)
Tscal[$dA_, $da_, $db_] :=
	Module[{temp, scal, gRep, group, entries,
		dim = {Length@ $fieldIndexMap["GaugeBosons"], Length@ $fieldIndexMap["Scalars"], Length@ $fieldIndexMap["Scalars"]}},
		If[MemberQ[dim, 0], Return@ 0;];
		(* A factor 1/2 included for the anti symmetrization of the external scalar indices *)
		entries = Flatten @ Table[
			group = Head@ If[Head@ gRep === Bar, gRep[[1]], gRep];
			temp = Product[del[rep, $da, $db], {rep,  DeleteCases[$scalars[scal, GaugeRep], gRep] ~ Join ~ $scalars[scal, FlavorIndices]}];
			{$fieldIndexMap["All"] /@ {group, Bar@ scal, scal} -> temp tGen[gRep, $dA, $da, $db] /2,
				$fieldIndexMap["All"] /@ {group, scal, Bar@ scal} -> -temp tGen[gRep, $dA, $db, $da] /2 }
		,{scal, Keys@ $scalars}, {gRep, $scalars[scal, GaugeRep]}];

		TStructure[$gauge@ $dA, $scalar@ $da, $scalar@ $db] @ SparseArray[GatherToList@ entries, dim]
	];

(*The general gauge structure constants G^{-2}_{AD} f^{DBC} used in the computation of the beta function tensors*)
FGauge[$dA_, $dB_, $dC_] :=
	Module[{entries, group, groupInfo, dim = Length@ $fieldIndexMap["GaugeBosons"] {1, 1, 1}},
		If[MemberQ[dim, 0], Return@ 0;];
		entries = Table[groupInfo = $gaugeGroups@ group;
				If[MatchQ[groupInfo@ LieGroup, U1[n_]],
					Nothing
				,
					$fieldIndexMap["GaugeBosons"] /@ (group {1, 1, 1}) -> Power[groupInfo[Coupling], -2] fStruct[group, $dA, $dB, $dC]
				]
			, {group, Keys@$gaugeGroups}];
		TStructure[$gauge@ $dA, $gauge@ $dB, $gauge@ $dC]@ SparseArray[GatherToList@ entries, dim]
	];

(*Genral Yukawa couplings used in the computation of the beta function tensors.*)
Yuk[$da_, $di_, $dj_, massive_:False] :=
	Module[{couplingInfo, yL, yR,
		dim = {Length@ $fieldIndexMap["Scalars"], Length@ $fieldIndexMap["Fermions"], Length@ $fieldIndexMap["Fermions"]} },
		If[MemberQ[dim, 0], Return@ {{0, 0}, {0, 0}};];
		yL = Table[
				CouplingTensorMap[couplingInfo, couplingInfo@ Fields, {$da, $di, $dj},
					{$da, $di, $dj}, couplingInfo@ UniqueArrangements]
			, {couplingInfo, $yukawas}];
		yR =  Table[
				CouplingTensorMap[couplingInfo, MapAt[Bar, couplingInfo@ Fields, 1], {$da, $di, $dj},
					{$da, $di, $dj}, couplingInfo@ UniqueArrangements, True]
			, {couplingInfo, $yukawas}];
		If[massive,
			yL = Join[yL, Table[
					CouplingTensorMap[couplingInfo, Join[{$vev}, couplingInfo@ Fields], {$di, $dj},
						{$da, $di, $dj}, couplingInfo@ UniqueArrangements]
				, {couplingInfo, $fermionMasses}] ];
			yR = Join[yR, Table[
					CouplingTensorMap[couplingInfo, Join[{$vev}, couplingInfo@ Fields], {$di, $dj},
						{$da, $di, $dj}, couplingInfo@ UniqueArrangements, True]
				, {couplingInfo, $fermionMasses}] ];
		];
		(* The chiral couplings are organized in full Yukawa matrix *)
		TStructure[$scalar@ $da, $fermion@ $di, $fermion@ $dj] /@
			{SparseArray[GatherToList@ yL, dim], SparseArray[GatherToList@ yR, dim]} //DiagonalMatrix
	];

YukTil[$da_, $di_, $dj_, massive_:False] := {{0, 1}, {1, 0}}.Yuk[$da, $di, $dj, massive].{{0, 1}, {1, 0}};

(*Genral quartic coupling used in the computation of the beta function tensors.*)
Lam[$da_, $db_, $dc_, $dd_, massive_:False] :=
	Module[{couplingInfo, temp, dim = Length@ $fieldIndexMap["Scalars"] {1, 1, 1, 1}},
		If[MemberQ[dim, 0], Return@ 0;];
		temp = Table[
				Join[CouplingTensorMap[couplingInfo, couplingInfo@ Fields, {$da, $db, $dc, $dd},
						{$da, $db, $dc, $dd}, couplingInfo@ UniqueArrangements],
					If[couplingInfo@ SelfConjugate, {}, (*Adds conjugate coupling if relevant *)
						CouplingTensorMap[couplingInfo, Bar/@ couplingInfo@ Fields, {$da, $db, $dc, $dd},
							{$da, $db, $dc, $dd}, couplingInfo@ UniqueArrangements, True]
					] ]
			, {couplingInfo, $quartics}];
		If[massive,
			temp = Join[temp,
				Table[
					Join[CouplingTensorMap[couplingInfo, Join[couplingInfo@ Fields, {$vev}], {$da, $db, $dc},
							{$da, $db, $dc, $dd}, couplingInfo@ UniqueArrangements],
						If[couplingInfo@ SelfConjugate, {}, (*Adds conjugate coupling if relevant *)
							CouplingTensorMap[couplingInfo, Join[Bar/@ couplingInfo@ Fields, {$vev}], {$da, $db, $dc},
								{$da, $db, $dc, $dd}, couplingInfo@ UniqueArrangements, True]
						] ], {couplingInfo, $trilinears}],
				Table[
					Join[CouplingTensorMap[couplingInfo, Join[couplingInfo@ Fields, {$vev, $vev}], {$da, $db},
							{$da, $db, $dc, $dd}, couplingInfo@ UniqueArrangements],
						If[couplingInfo@ SelfConjugate, {}, (*Adds conjugate coupling if relevant *)
							CouplingTensorMap[couplingInfo, Join[Bar/@ couplingInfo@ Fields, {$vev, $vev}], {$da, $db},
								{$da, $db, $dc, $dd}, couplingInfo@ UniqueArrangements, True]
						] ], {couplingInfo, $scalarMasses}]
			];
		];
		TStructure[$scalar@$da, $scalar@$db, $scalar@$dc, $scalar@$dd] @ SparseArray[GatherToList@ temp, dim]
	];


(*######################################*)
(*---------------Clean up---------------*)
(*######################################*)
(*Removes all model information.*)
ResetModel[] :=
	Module[{},
		(*Resets tensor dummy notation*)
		ReInitializeSymbols[];
		(* Resets coupling symmetries and reality *)
		ReInitializeCouplingBehavior[];

		(*Resets all information of the current model.*)
		(*Global couplings variable*)
		$couplings = <||>;
		(*List of all flavor representations*)
		$flavorReps = {};
		(*Association with all information on the gauge groups: fields, couplings etc.*)
		$gaugeGroups = <||>;
		(*Association with all information on the fermion fields: representations etc.*)
		$fermions = <||>;
		(*Association with all information on the scalar fields: representations, etc.*)
		$scalars = <||>;
		(*Association with all information on the quartic couplings.*)
		$quartics = <||>;
		(*Association with all information on the Yukawa couplings.*)
		$yukawas = <||>;
		(*Association with all information on the fermion masses.*)
		$fermionMasses = <||>;
		(*Association with all information on the trilinear scalar couplings.*)
		$trilinears = <||>;
		(*Association with all information on the scalar masses.*)
		$scalarMasses = <||>;
		(*Association with all information on the fermion fields: representations etc.*)
		$LieGroups = <||>;
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
	Module[{},
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
	Module[{},
		RemoveInteraction[coupling];
		RemoveInteraction[couplings];
	];
