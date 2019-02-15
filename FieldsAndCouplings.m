(*Structure deltas*)
Clear[SdelS, SdelF, SdelV]
SdelS /: SdelS[field1_, ind_, f_] SdelS[field2_, ind_, g_] := 0; 
SdelF /: SdelF[field1_, ind_, f_] SdelF[field2_, ind_, g_] := 0;
SdelV /: SdelV[field1_, ind_, f_] SdelV[field2_, ind_, g_] := 0;

(*Associationwith all information on the scalar fields: representations, etc.*)
$scalars = <||>;
(*Initiates a scalar field*)
Options[CreateScalar] = {SelfConjugate -> False, GaugeRep -> {}, FlavorIndices -> {}};
CreateScalar[field_, OptionsPattern[] ] :=
	Block[{rep},
		AppendTo[$scalars, field -> <|GaugeRep -> OptionValue[GaugeRep], FlavorIndices -> OptionValue[FlavorIndices],
			SelfConjugate -> OptionValue[SelfConjugate]|>];
		If[OptionValue[SelfConjugate],
			SdelS/: SdelS[field, ind_, s1_] SdelS[field, ind_, s2_] = 1 
				* Product[del[rep, s1, s2], {rep, OptionValue[GaugeRep]}]
				* Product[del[rep, s1, s2], {rep, OptionValue[FlavorIndices]}];
			Bar[field] = field;
		,
			SdelS/: SdelS[field, ind_, s1_] SdelS[Bar[field], ind_, s2_] = 2 
				* Product[del[rep, s1, s2], {rep, OptionValue[GaugeRep]}]
				* Product[del[rep, s1, s2], {rep, OptionValue[FlavorIndices]}];
		];
	];

(*Associationwith all information on the fermion fields: representations etc.*)
$fermions = <||>;
(*Initiates a fermion field*)	
Options[CreateFermion] = {GaugeRep -> {}, FlavorIndices -> {}};
CreateFermion[field_, OptionsPattern[] ] :=
	Block[{rep},
		AppendTo[$fermions, field -> <|GaugeRep -> OptionValue[GaugeRep], FlavorIndices -> OptionValue[FlavorIndices]|>]; 
		SdelF/: SdelF[field, ind_, f1_] SdelF[Bar[field], ind_, f2_] = 
			Product[del[rep, f1, f2], {rep, OptionValue[GaugeRep]}]
			* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];
	];

(*Initiates a vector field*)	
CreateVector[name_, group_] :=
	Block[{},
		SdelV/: SdelV[name, ind_, v1_] SdelV[name, ind_, v2_] = del[group[adj], v1, v2];
	];

(*###################################*)
(*----------Gauge couplings----------*)
(*###################################*)
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
	Module[{ferm, rep, gRep1, gRep2, f1, f2, v}, 
		Sum[
			SdelF[Bar @ ferm, i, f1] SdelF[ferm, j, f2] * Product[del[rep, f1, f2], {rep, $fermions[ferm, FlavorIndices]}]
			* Sum[SdelV[$gaugeGroups[Head @ gRep1, Field], A, v] TGen[gRep1, v, f1, f2] * Product[del[gRep2, f1, f2],
				{gRep2, DeleteCases[$fermions[ferm, GaugeRep], gRep1]}],{gRep1, $fermions[ferm, GaugeRep]}] 
		,{ferm, Keys @ $fermions}]
	];
(*And for the Majorana-spinor*)
Tferm[A_, i_, j_] := {{TfLeft[A, i, j], 0}, {0, -TfLeft[A, j, i]}};
TfermTil[A_, i_, j_] := {{-TfLeft[A, j, i], 0}, {0, TfLeft[A, i, j]}};  

(*Defining the anti-symmetric gauge generators for the scalars*)
Ts[A_, a_, b_] := 
	Module[{scal, rep, gRep1, gRep2, s1, s2, v}, 
		Sum[
			AntiSym[a, b][SdelS[Bar @ scal, a, s1] SdelS[scal, b, s2]]  
			* Product[del[rep, s1, s2], {rep, $scalars[scal, FlavorIndices]}]
			* Sum[SdelV[$gaugeGroups[Head @ gRep1, Field], A, v] TGen[gRep1, v, s1, s2] * Product[del[gRep2, s1, s2],
				{gRep2, DeleteCases[$scalars[scal, GaugeRep], gRep1]}],{gRep1, $scalars[scal, GaugeRep]}] 
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
Options[AddYukawa] = {OverallFactor -> 1, Chirality -> Left, InvarianceCheck -> False};
AddYukawa::chirality = "The chirality `1` is invalid: Left or Right expected. Returning Null"; 
AddYukawa[coupling_, {phi_, psi1_, psi2_}, indices_Function, groupInvariant_Function, OptionsPattern[]] :=
	Block[{g, group, invariance, normalization, projection, symmetryFactor, temp, test, yuk, yukbar, y}, 
		normalization = 2 OptionValue @ OverallFactor;
		
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
		If[OptionValue @ InvarianceCheck,
			y = With[{y1 = yuk, ind = Sequence@@ indices[s, f1, f2], gi = groupInvariant[s, f1, f2]},
				Sym[#2, #3][y1[ind] SdelS[phi, #1, s] SdelF[psi1, #2, f1] SdelF[psi2, #3, f2] gi ] &];
			test = TfLeft[A, k, i] y[a, k ,j] + y[a, i, k] TfLeft[A, k, j] + y[b, i, j] Ts[A, b, a]//Expand;
			test = test SdelS[Bar@phi, a, scal] SdelF[Bar@psi1, i, ferm1] SdelF[Bar@psi2, j, ferm2] //Expand;
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
		Switch[OptionValue @ Chirality
		,Left,
			projection = With[{c = normalization/2 /symmetryFactor / Expand @ Power[OptionValue[OverallFactor] groupInvariant[a,b,c], 2] },
			c SdelS[Bar@phi, #1, s] SdelF[Bar@psi1, #2, f1] SdelF[Bar@psi2, #3, f2] groupInvariant[s, f1, f2] &];	
		,Right,
			projection = With[{c = normalization/2 /symmetryFactor / Expand @ Power[OptionValue[OverallFactor] groupInvariant[a,b,c], 2] },
			c SdelS[phi, #1, s] SdelF[psi1, #2, f1] SdelF[psi2, #3, f2] groupInvariant[s, f1, f2] &];
		];
		
		(*Adds the Yukawa coupling to the association*)
		AppendTo[$yukawas, coupling -> 
			<|Chirality -> OptionValue @ Chirality,
			Coupling -> yuk,
			CouplingBar -> yukbar, 
			Fields -> {phi, psi1, psi2},
			Indices -> indices,
			Invariant -> groupInvariant,
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
Options[AddQuartic] = {OverallFactor -> 1, SelfConjugate -> True, InvarianceCheck -> False}; 
AddQuartic [coupling_, {phi1_, phi2_, phi3_, phi4_}, indices_Function, groupInvariant_Function, OptionsPattern[]] :=
	Block[{group, invariance, lam, lambar, lambda,  normalization, phi, projection, symmetryFactor, temp},
		normalization = 24 OptionValue[OverallFactor];
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
			lambda = With[{l1 = lam, ind = Sequence@@ indices[s1, s2, s3, s4], gi = groupInvariant[s1, s2, s3, s4]},
				Sym[#1, #2, #3, #4][l1[ind] SdelS[phi1, #1, s1] SdelS[phi2, #2, s2] SdelS[phi3, #3, s3] SdelS[phi4, #4, s4] gi] &];
			test = Ts[A, a, e] lambda[e, b, c, d] + Ts[A, b, e] lambda[a, e, c, d] 
				+ Ts[A, c, e] lambda[a, b, e, d] + Ts[A, d, e] lambda[a, b, c, e]//Expand;
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
		projection = With[{c = normalization /24 /symmetryFactor 
			/ Expand @ Power[OptionValue[OverallFactor] groupInvariant[a, b, c, d], 2] },
			c SdelS[Bar@phi1, #1, s1] SdelS[Bar@phi2, #2, s2] SdelS[Bar@phi3, #3, s3] 
			* SdelS[Bar@phi4, #4, s4] groupInvariant[s1, s2, s3, s4] &];
		
		(*Adds the quartic coupling to the association*)
		AppendTo[$quartics, coupling -> 
			<|Coupling -> lam,
			CouplingBar -> lambar,
			Fields -> {phi1, phi2, phi3, phi4},
			Indices -> indices,
			Invariant -> groupInvariant,
			Projector -> projection,
			SelfConjugate -> OptionValue[SelfConjugate]|>];
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




