(*Structure deltas*)
Clear[SdelS, SdelF, SdelV]
SdelS /: SdelS[field1_, ind_, f_] SdelS[field2_, ind_, g_] := 0; 
SdelF /: SdelF[field1_, ind_, f_] SdelF[field2_, ind_, g_] := 0;
SdelV /: SdelV[field1_, ind_, f_] SdelV[field2_, ind_, g_] := 0;

(*Associationwith all information on the scalar fields: representations, etc.*)
$scalars = <||>;
(*Initiates a scalar field*)
Options[AddScalar] = {SelfConjugate -> False, GaugeRep -> {}, FlavorIndices -> {}, Mass -> 0};
AddScalar[field_String, OptionsPattern[] ] :=
	Block[{rep, massTerm = 1},
		AppendTo[$scalars, field -> <|
			GaugeRep -> OptionValue[GaugeRep], 
			FlavorIndices -> OptionValue[FlavorIndices],
			SelfConjugate -> OptionValue[SelfConjugate], 
			Mass -> OptionValue @ Mass|>];
		If[OptionValue @ Mass =!= 0,
			massTerm = UnitStep[t - OptionValue @ Mass];
		];
		
		If[OptionValue[SelfConjugate],
			SdelS/: SdelS[field, ind_, s1_] SdelS[field, ind_, s2_] = 1 * massTerm 
				* Product[del[rep, s1, s2], {rep, OptionValue[GaugeRep]}]
				* Product[del[rep, s1, s2], {rep, OptionValue[FlavorIndices]}];
			Bar[field] = field;
		,
			SdelS/: SdelS[field, ind_, s1_] SdelS[Bar[field], ind_, s2_] = 2 * massTerm
				* Product[del[rep, s1, s2], {rep, OptionValue[GaugeRep]}]
				* Product[del[rep, s1, s2], {rep, OptionValue[FlavorIndices]}];
		];
	];

(*Associationwith all information on the fermion fields: representations etc.*)
$fermions = <||>;
(*Initiates a fermion field*)	
Options[AddFermion] = {GaugeRep -> {}, FlavorIndices -> {}, Mass -> 0};
AddFermion[field_String, OptionsPattern[] ] :=
	Block[{massTerm = 1, rep},
		If[OptionValue @ Mass =!= 0,
			massTerm = UnitStep[t - OptionValue @ Mass];
		];
		AppendTo[$fermions, field -> <|
			GaugeRep -> OptionValue[GaugeRep], 
			FlavorIndices -> OptionValue[FlavorIndices], 
			Mass -> OptionValue @ Mass|>]; 
		SdelF/: SdelF[field, ind_, f1_] SdelF[Bar[field], ind_, f2_] = massTerm 
			* Product[del[rep, f1, f2], {rep, OptionValue[GaugeRep]}]
			* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];
	];

(*Initiates a vector field*)	
AddVector[name_String, group_] :=
	Block[{},
		SdelV/: SdelV[name, ind_, v1_] SdelV[name, ind_, v2_] = del[group[adj], v1, v2];
	];

(*###################################*)
(*----------Gauge couplings----------*)
(*###################################*)
(*Association with all information on the gauge groups: fields, couplings etc.*)
$gaugeGroups = <||>;

(*Function for adding gauge groups to the model*)
AddGaugeGroup::unkown = "`1` is not a reckognized Lie group."
AddGaugeGroup::nonstring = "Use a string for the field or leave it as the default value."
AddGaugeGroup[coupling_Symbol, groupName_Symbol, lieGroup_[n_], OptionsPattern[{Field -> Automatic}]] :=
	Block[{projector, fieldName},
		(*Sets the field name*)
		Switch[OptionValue @ Field
		,Automatic,
			fieldName = "A_" <> ToString @ groupName;
		,_String,
			fieldName = OptionValue @ Field;
		,_,
			Message[AddGaugeGroup::nonstring];
			Return @ Null;
		];
		
		(*Sets up the group symbols*)
		Switch[lieGroup
		,SU,
			DefineSUGroup[groupName, n];
		,SO,
			DefineSOGroup[groupName, n];
		,U,
			If[n =!= 1,
				Message[AddGaugeGroup::unkown, lieGroup[n] ];
				Return @ Null;
			];
			DefineU1Group[groupName];
		,_,
			Message[AddGaugeGroup::unkown, lieGroup[n] ];
			Return @ Null;
		];
		
		(*Sets up the gauge fields and 2-point projection*)
		AddVector[fieldName, groupName];
		projector = With[{V = fieldName, gStruct = del[groupName[adj], v1, v2] / Dim @ groupName[adj]}, 
			SdelV[V, #1, v1] SdelV[V, #2, v2] gStruct &]; 
		
		(*Adds the group information and the coupling to the repsective lists*)
		AppendTo[$gaugeGroups, groupName -> 
			<|Field -> fieldName, 
			Coupling -> coupling, 
			Projector -> projector|>];
		AppendTo[$couplings, coupling -> groupName];
	];



(*The gauge coupling matrix G^2_{AB}*)
G2[A_, B_, power_: 1] := 
	Module[{gauge, v1, v2},
		Sum[
			SdelV[$gaugeGroups[gauge, Field], A, v1] SdelV[$gaugeGroups[gauge, Field], B, v2] 
				*Power[$gaugeGroups[gauge, Coupling], 2 power] del[gauge[adj], v1, v2]  
		,{gauge, Keys @ $gaugeGroups}]
	];
	
(*Defining the gauge generators for the left-handed spinors, \bar{\psi} T^A \psi.*)
TfLeft[A_, i_, j_] := 
	Module[{ferm, group, rep, gRep1, gRep2, f1, f2, v}, 
		Sum[
			SdelF[Bar @ ferm, i, f1] SdelF[ferm, j, f2] * Product[del[rep, f1, f2], {rep, $fermions[ferm, FlavorIndices]}]
			* Sum[group = Head@ If[Head@ gRep1 === Bar, gRep1[[1]], gRep1]; 
					SdelV[$gaugeGroups[group, Field], A, v] tGen[gRep1, v, f1, f2] 
					* Product[del[gRep2, f1, f2], {gRep2, DeleteCases[$fermions[ferm, GaugeRep], gRep1]}], 
				{gRep1, $fermions[ferm, GaugeRep]}] 
		,{ferm, Keys @ $fermions}]
	];
(*And for the Majorana-spinor*)
Tferm[A_, i_, j_] := {{TfLeft[A, i, j], 0}, {0, -TfLeft[A, j, i]}};
TfermTil[A_, i_, j_] := {{-TfLeft[A, j, i], 0}, {0, TfLeft[A, i, j]}};  

(*Defining the anti-symmetric gauge generators for the scalars*)
Tscal[A_, a_, b_] := 
	Module[{scal, group, rep, gRep1, gRep2, s1, s2, v}, 
		Sum[
			AntiSym[a, b][SdelS[Bar @ scal, a, s1] SdelS[scal, b, s2]]  
			* Product[del[rep, s1, s2], {rep, $scalars[scal, FlavorIndices]}]
			* Sum[group = Head@ If[Head@ gRep1 === Bar, gRep1[[1]], gRep1];
					SdelV[$gaugeGroups[group, Field], A, v] tGen[gRep1, v, s1, s2] 
					* Product[del[gRep2, s1, s2], {gRep2, DeleteCases[$scalars[scal, GaugeRep], gRep1]}],
				{gRep1, $scalars[scal, GaugeRep]}] 
		,{scal, Keys @ $scalars}]
	];

(*Defining the gauge structure constant multiplied with G^-2*)
FGauge[A_, B_, C_] := 
	Module[{gauge, v1, v2, v3},
		Sum[
			SdelV[$gaugeGroups[gauge, Field], A, v1] SdelV[$gaugeGroups[gauge, Field], B, v2] SdelV[$gaugeGroups[gauge, Field], C, v3]
				* Power[$gaugeGroups[gauge, Coupling], -2] fStruct[gauge, v1, v2, v3]  
		,{gauge, Keys @ $gaugeGroups}]
	];

(*####################################*)
(*----------Yukawa couplings----------*)
(*####################################*)
(*Associationwith all information on the Yukawa couplings.*)
$yukawas = <||>;
(*Function for defining the Yukawa couplings of the theory*)
Options[AddYukawa] = {CouplingIndices -> (Null &),
	GroupInvariant -> (1 &),
	Chirality -> Left, 
	CheckInvariance -> False};
AddYukawa::chirality = "The chirality `1` is invalid: Left or Right expected."; 
AddYukawa::unkown = "`1` does not match any of the `2`s.";
AddYukawa::nonfunction = "The value given in `1` is not a function."; 
AddYukawa[coupling_, {phi_, psi1_, psi2_}, OptionsPattern[]] :=
	Block[{g, group, normalization, projection, symmetryFactor, temp, test, yuk, yukbar, y}, 
		(*Tests if the fields have been defined*)
		If[!MemberQ[Keys @ $scalars, phi /. Bar[x_] -> x],
			Message[AddYukawa::unkown, phi /. Bar[x_] -> x, "scalar"];
			Return @ Null;
		];
		If[!MemberQ[Keys@ $fermions, psi1],
			Message[AddYukawa::unkown, psi1, "fermion"];
			Return @ Null;
		];
		If[!MemberQ[Keys@ $fermions, psi2],
			Message[AddYukawa::unkown, psi2, "fermion"];
			Return @ Null;
		];

		(*Tests if the fields have been defined*)
		If[! Head @ OptionValue @ CouplingIndices === Function,
			Message[AddYukawa::nonfunction, CouplingIndices];
			Return @ Null;
		]; 
		If[! Head @ OptionValue @ GroupInvariant === Function,
			Message[AddYukawa::nonfunction, GroupInvariant];
			Return @ Null;
		]; 
		
		(*If the chirality is right handed, the coupling is written with the barred fields*)
		Switch[OptionValue @ Chirality
		,Left,
			g = coupling;
		,Right, 	
			g = Bar @ coupling;
		,_,
			Message[AddYukawa::chirality, OptionValue @ Chirality];
			Return @ Null;
		];
		
		(*Constructs the coupling structure*)
		normalization = 2; (*To account for the symmetrization in Yukawa[a, i, j]*)
		If[!$scalars[phi, SelfConjugate]|| Head @ phi === Bar, normalization *= 1/Sqrt[2];];
		With[{n = normalization, g0 = g},
			If[Length @ coupling <= 2,
				yuk = n Matrix[g0]@ ## &;
				yukbar = n Matrix[Bar @g0]@ ## &
			,
				yuk = n * g0 @ ## &;
				yukbar = n Bar[g0]@ ## &
			];
		];
		
		(*Tests whether the Yukawa coupling satisfy gauge invariance*)
		If[OptionValue @ CheckInvariance,
			y = With[{y1 = yuk, ind = Sequence@@ OptionValue[CouplingIndices][s, f1, f2], 
				gi = OptionValue[GroupInvariant][s, f1, f2]},
				Sym[#2, #3][y1[ind] SdelS[phi, #1, s] SdelF[psi1, #2, f1] SdelF[psi2, #3, f2] gi ] &];
			test = TfLeft[A, k, i] y[a, k ,j] + y[a, i, k] TfLeft[A, k, j] + y[b, i, j] Tscal[A, b, a] //Expand;
			test *= SdelS[Bar@phi, a, scal] SdelF[Bar@psi1, i, ferm1] SdelF[Bar@psi2, j, ferm2] //Expand;
			Do[
				temp = test SdelV[group @ Field, A, vec1] //Expand;
				If [temp =!= 0, 
					Print[coupling,"---Gauge invarinace check inconclusive for the ", group @ Field, " field:"];
					Print["0 = ", temp];
				];
			,{group, $gaugeGroups}];
		];
		
		(*Defines the projection operator for extracting out the particular Yukawa coupling.*)
		symmetryFactor = If[psi1 === psi2, 2, 1];
		projection = With[{c = normalization/2 /symmetryFactor/ Expand[OptionValue[GroupInvariant][a,b,c] 
				*(OptionValue[GroupInvariant][a,b,c] /. tGen[rep_, A_, a_, b_] -> tGen[Bar @ rep, A, a, b])], 
				gInv = OptionValue[GroupInvariant][s, f1, f2] /. tGen[rep_, A_, a_, b_] -> tGen[Bar @ rep, A, a, b]}, 
			Switch[OptionValue @ Chirality
			,Left,
				c SdelS[Bar@ phi, #1, s] SdelF[Bar@ psi1, #2, f1] SdelF[Bar@ psi2, #3, f2] gInv &	
			,Right,
				c SdelS[phi, #1, s] SdelF[psi1, #2, f1] SdelF[psi2, #3, f2] gInv &
			]
		];
		
		(*Adds the Yukawa coupling to the association*)
		AppendTo[$yukawas, coupling -> 
			<|Chirality -> OptionValue @ Chirality,
			Coupling -> yuk,
			CouplingBar -> yukbar, 
			Fields -> {phi, psi1, psi2},
			Indices -> OptionValue @ CouplingIndices,
			Invariant -> OptionValue @ GroupInvariant,
			Projector -> projection|>];
		AppendTo[$couplings, coupling -> Yukawa];
	];

YukawaLeft[a_, i_, j_] :=
	Module[{f1, f2, yu, s1},
		Sum[SdelS[yu[Fields][[1]], a, s1] SdelF[yu[Fields][[2]], i, f1] SdelF[yu[Fields][[3]], j, f2]
				* yu[Invariant][s1, f1, f2] yu[Coupling][Sequence @@ yu[Indices][s1, f1, f2]]  
			,{yu, $yukawas}] //Sym[i, j]
	];

YukawaRight[a_, i_, j_] :=
	Module[{f1, f2, yu, s1},
		Sum[SdelS[Bar @ yu[Fields][[1]], a, s1] SdelF[Bar @ yu[Fields][[2]], i, f1] SdelF[Bar @ yu[Fields][[3]], j, f2]
				* yu[Invariant][s1, f1, f2] yu[CouplingBar][Sequence @@ yu[Indices][s1, f1, f2]]  
			,{yu, $yukawas}] //Sym[i, j]
	];

Yuk[a_, i_, j_] := {{YukawaLeft[a, i, j], 0}, {0, YukawaRight[a, i, j]}};
YukTil[a_, i_, j_] := {{YukawaRight[a, i, j], 0}, {0, YukawaLeft[a, i, j]}};


(*#####################################*)
(*----------Quartic couplings----------*)
(*#####################################*)
(*Associationwith all information on the quartic couplings.*)
$quartics = <||>;
(*Function for defining the Yukawa couplings of the theory*)
AddQuartic::unkown = "`1` does not match any of the scalars.";
AddQuartic::nonfunction = "The value given in `1` is not a function.";
Options[AddQuartic] = {CouplingIndices -> (Null &),
	GroupInvariant -> (1 &),
	SelfConjugate -> True, 
	InvarianceCheck -> False}; 
AddQuartic [coupling_, {phi1_, phi2_, phi3_, phi4_}, OptionsPattern[]] :=
	Block[{group, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
		(*Tests if the fields have been defined*)
		Do[
			If[!MemberQ[Keys@ $scalars, temp /. Bar[x_] -> x],
				Message[AddQuartic::unkown, temp /. Bar[x_] -> x];
				Return[Null];
			];
		,{temp, {phi1, phi2, phi3, phi4}}];
		
		(*Tests if the fields have been defined*)
		If[! Head @ OptionValue @ CouplingIndices === Function,
			Message[AddQuartic::nonfunction, CouplingIndices];
			Return @ Null;
		]; 
		If[! Head @ OptionValue @ GroupInvariant === Function,
			Message[AddQuartic::nonfunction, GroupInvariant];
			Return @ Null;
		]; 
		
		normalization = 24; (*To account for the symmetrization in Lam[a, b, c, d]*)
		Do[
			If[!$scalars[phi, SelfConjugate]|| Head @ phi === Bar, normalization *= 1/Sqrt[2];];
		,{phi, {phi1, phi2, phi3, phi4}}];
		With[{n = normalization},
			If[Length @ coupling <= 2,
				lam = n Matrix[coupling]@ ## &;
				lambar = n Matrix[Bar @ coupling]@ ## &
			,
				lam = n coupling @ ## &;
				lambar = n Bar[coupling] @ ## &
			];
		];
		(*Tests whether the Yukawa coupling satisfy gauge invariance*)
		If[OptionValue @ InvarianceCheck,
			lambda = With[{l1 = lam, ind = Sequence@@ OptionValue @ CouplingIndices[s1, s2, s3, s4], 
				gi = OptionValue[GroupInvariant][s1, s2, s3, s4]},
				Sym[#1, #2, #3, #4][l1[ind] SdelS[phi1, #1, s1] SdelS[phi2, #2, s2] SdelS[phi3, #3, s3] SdelS[phi4, #4, s4] gi] &];
			test = Tscal[A, a, e] lambda[e, b, c, d] + Tscal[A, b, e] lambda[a, e, c, d] 
				+ Tscal[A, c, e] lambda[a, b, e, d] + Tscal[A, d, e] lambda[a, b, c, e]//Expand;
			test = test SdelS[Bar@phi1, a, scal1] SdelS[Bar@phi2, b, scal2] SdelS[Bar@phi3, c, scal3] SdelS[Bar@phi4, d, scal4] //Expand;
			Do[
				temp = test SdelV[group @ Field, A, vec1] //Expand;
				If [temp =!= 0, 
					Print[coupling,"---Gauge invarinace check inconclusive for the ", group @ Field, " field:"];
					Print["0 = ", temp//S];
				];
			,{group, $gaugeGroups}];
		];
		
		(*Defines the projection operator for extracting out the particular quartic coupling.*)
		symmetryFactor = 24 / Length @ DeleteDuplicates @ Permutations @ {phi1, phi2, phi3, phi4};	
		projection = With[{c = normalization /24 /symmetryFactor / Expand[OptionValue[GroupInvariant][a, b, c, d] 
				* (OptionValue[GroupInvariant][a, b, c, d] /. tGen[rep_, A_, a_, b_] -> tGen[Bar @ rep, A, a, b])],
				gInv = OptionValue[GroupInvariant][s1, s2, s3, s4] /. tGen[rep_, A_, a_, b_] -> tGen[Bar @ rep, A, a, b] },
			c SdelS[Bar@phi1, #1, s1] SdelS[Bar@phi2, #2, s2] SdelS[Bar@phi3, #3, s3] * SdelS[Bar@phi4, #4, s4] gInv &
			];
		
		(*Adds the quartic coupling to the association*)
		AppendTo[$quartics, coupling -> 
			<|Coupling -> lam,
			CouplingBar -> lambar,
			Fields -> {phi1, phi2, phi3, phi4},
			Indices -> OptionValue @ CouplingIndices,
			Invariant -> OptionValue @ GroupInvariant,
			Projector -> projection,
			SelfConjugate -> OptionValue @ SelfConjugate|>];
		AppendTo[$couplings, coupling -> Quartic];
	];

Lam[a_, b_, c_, d_] :=
	Module[{l, s1, s2, s3, s4},
		Sum[SdelS[l[Fields][[1]], a, s1] SdelS[l[Fields][[2]], b, s2] SdelS[l[Fields][[3]], c, s3] SdelS[l[Fields][[4]], d, s4]
				* l[Invariant][s1, s2, s3, s4] l[Coupling][Sequence @@ l[Indices][s1, s2, s3, s4]]
			+If[! l@SelfConjugate,
				SdelS[Bar@ l[Fields][[1]], a, s1] SdelS[Bar@ l[Fields][[2]], b, s2] SdelS[Bar@ l[Fields][[3]], c, s3] 
				* SdelS[Bar@ l[Fields][[4]], d, s4] l[Invariant][s1, s2, s3, s4] l[CouplingBar][Sequence @@ l[Indices][s1, s2, s3, s4]]
				,0] 
			,{l, $quartics}] //Sym[a, b, c, d]
	];




