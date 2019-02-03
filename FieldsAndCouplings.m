Clear[Bar];
Bar[Bar[field_]] := field; 

(*Structure deltas*)
Clear[SdelS, SdelF, SdelV]
SdelS /: SdelS[field1_, ind_, f_] SdelS[field2_, ind_, g_] := 0; 
SdelF /: SdelF[field1_, ind_, f_] SdelF[field2_, ind_, g_] := 0;
SdelV /: SdelV[field1_, ind_, f_] SdelV[field2_, ind_, g_] := 0;

(*Associationwith all information on the scalar fields: representations, etc.*)
scalars = <||>;
(*Initiates a scalar field*)
Options[CreateScalar] = {SelfConjugate -> False, GaugeRep -> {}, FlavorIndices -> {}};
CreateScalar[field_, OptionsPattern[] ] :=
	Block[{rep},
		AppendTo[scalars, field -> <|GaugeRep -> OptionValue[GaugeRep], FlavorIndices -> OptionValue[FlavorIndices],
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
fermions = <||>;
(*Initiates a fermion field*)	
Options[CreateFermion] = {GaugeRep -> {}, FlavorIndices -> {}};
CreateFermion[field_, OptionsPattern[] ] :=
	Block[{rep},
		AppendTo[fermions, field -> <|GaugeRep -> OptionValue[GaugeRep], FlavorIndices -> OptionValue[FlavorIndices]|>]; 
		SdelF/: SdelF[field, ind_, f1_] SdelF[Bar[field], ind_, f2_] = 
			Product[del[rep, f1, f2], {rep, OptionValue[GaugeRep]}]
			* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];
	];

(*Initiates a vector field*)	
CreateVector[name_, group_] :=
	Block[{},
		SdelV/: SdelV[name, ind_, v1_] SdelV[name, ind_, v2_] = del[group[adj], v1, v2];
	];


(*----------Gauge couplings----------*)
(*The gauge coupling matrix G^2_{AB}*)
G2[A_, B_, power_: 1] := 
	Module[{gauge, v1, v2},
		Sum[
			SdelV[gaugeGroups[gauge, Field], A, v1] SdelV[gaugeGroups[gauge, Field], B, v2] 
				*Power[gaugeGroups[gauge, Coupling], 2 power] del[gauge[adj], v1, v2]  
		,{gauge, Keys @ gaugeGroups}]
	];
	
(*Defining the gauge generators for the left-handed spinors, \bar{\psi} T^A \psi.*)
TfLeft[A_, i_, j_] := 
	Module[{ferm, rep, gRep1, gRep2, f1, f2, v}, 
		Sum[
			SdelF[Bar @ ferm, i, f1] SdelF[ferm, j, f2] * Product[del[rep, f1, f2], {rep, fermions[ferm, FlavorIndices]}]
			* Sum[SdelV[gaugeGroups[Head @ gRep1, Field], A, v] TGen[gRep1, v, f1, f2] * Product[del[gRep2, f1, f2],
				{gRep2, DeleteCases[fermions[ferm, GaugeRep], gRep1]}],{gRep1, fermions[ferm, GaugeRep]}] 
		,{ferm, Keys @ fermions}]
	];
(*And for the Majorana-spinor*)
Tf[A_, i_, j_] := {{TfLeft[A, i, j], 0}, {0, -TfLeft[A, j, i]}};
Tft[A_, i_, j_] := {{-TfLeft[A, j, i], 0}, {0, TfLeft[A, i, j]}};  

(*Defining the anti-symmetric gauge generators for the scalars*)
Ts[A_, a_, b_] := 
	Module[{scal, rep, gRep1, gRep2, s1, s2, v}, 
		Sum[
			AntiSym[a, b][SdelS[Bar @ scal, a, s1] SdelS[scal, b, s2]]  
			* Product[del[rep, s1, s2], {rep, scalars[scal, FlavorIndices]}]
			* Sum[SdelV[gaugeGroups[Head @ gRep1, Field], A, v] TGen[gRep1, v, s1, s2] * Product[del[gRep2, s1, s2],
				{gRep2, DeleteCases[scalars[scal, GaugeRep], gRep1]}],{gRep1, scalars[scal, GaugeRep]}] 
		,{scal, Keys @ scalars}]
	];

(*Defining the gauge structure constant multiplied with G^-2*)
FGauge[A_, B_, C_] := 
	Module[{gauge, v1, v2, v3},
		Sum[
			SdelV[gaugeGroups[gauge, Field], A, v1] SdelV[gaugeGroups[gauge, Field], B, v2] SdelV[gaugeGroups[gauge, Field], C, v3]
				* Power[gaugeGroups[gauge, Coupling], -2] fStruct[gauge, v1, v2, v3]  
		,{gauge, Keys @ gaugeGroups}]
	];


(*----------Yukawa couplings----------*)

(*Associationwith all information on the Yukawa couplings.*)
yukawas = <||>;
AddYukawa[coupling_, {phi_, psi1_, psi2_}, indices_Function, groupInvariant_Function] :=
	Block[{group, invariance, normalization = 1, projection, temp, test, yuk, y}, 
		If[!scalars[phi, SelfConjugate]|| Head @ phi === Bar, normalization = 1/Sqrt[2];];
		If[Length @ coupling <= 2,
			yuk = Matrix @ coupling;
		,
			yuk = coupling;
		];
		
		(*Tests whether the Yukawa coupling satisfy gauge invariance*)
		y = With[{c = normalization, y1 = yuk, ind = Sequence@@ indices[s, f1, f2], gi = groupInvariant[s, f1, f2]},
			Sym[#2, #3][c * y1[ind] SdelS[phi, #1, s] SdelF[psi1, #2, f1] SdelF[psi2, #3, f2] gi ] &];
		test = TfLeft[A, k, i] y[a, k ,j] + y[a, i, k] TfLeft[A, k, j] + y[b, i, j] Ts[A, b, a]//Expand;
		test = test SdelS[Bar@phi, a, scal] SdelF[Bar@psi1, i, ferm1] SdelF[Bar@psi2, j, ferm2] //Expand;
		Do[
			temp = test SdelV[group @ Field, A, vec1] //Expand;
			If [temp =!= 0, 
				Print[coupling,"---Gauge invarinace check inconclusive for the ", group @ Field, " field:"];
				Print["0 = ", temp];
			];
		,{group, gaugeGroups}];
		
		(*Defines the projection operator for extracting out the particular Yukawa coupling.*)	
		projection = With[{c = normalization / Expand @ Power[groupInvariant[a,b,c], 2] },
			c SdelS[Bar@phi, #1, s] SdelF[Bar@psi1, #2, f1] SdelF[Bar@psi2, #3, f2] groupInvariant[s, f1, f2] &];	
			
		AppendTo[yukawas, coupling -> <|Coupling -> yuk normalization, Fields -> {phi, psi1, psi2}, Indices -> indices,
			Invaraint -> groupInvariant, Projector -> projection|>];
	];





