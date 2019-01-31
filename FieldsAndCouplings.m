(*Structure deltas*)
Clear[SdelS, SdelF, SdelV]
SdelS /: SdelS[field1_, ind_, f_] SdelS[field2_, ind_, g_] := 0; 
SdelF /: SdelF[field1_, ind_, f_] SdelF[field2_, ind_, g_] := 0;
SdelV /: SdelV[field1_, ind_, f_] SdelV[field2_, ind_, g_] := 0;

(*Initiates a scalar field*)
Options[CreateScalar] = {SelfConjugate -> False, GaugeRepresentation -> {}, FlavorIndices -> {}};
CreateScalar[name_, OptionsPattern[] ] :=
	Block[{},
		If[OptionValue[SelfConjugate],
			SdelS/: SdelS[name, ind_, f1_] SdelS[name, ind_, f2_] = 1 
				* Product[del[rep, f1, f2], {rep, OptionValue[GaugeRepresentation]}]
				* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];
		,
			SdelS/: SdelS[name, ind_, f1_] SdelS[Bar[name], ind_, f2_] = 2 
				* Product[del[rep, f1, f2], {rep, OptionValue[GaugeRepresentation]}]
				* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];
		];
	];

(*Associationwith all information on the fermion fields: representations etc.*)
fermions = <||>;
(*Initiates a fermion field*)	
Options[CreateFermion] = {GaugeRep -> {}, FlavorIndices -> {}};
CreateFermion[field_, OptionsPattern[] ] :=
	Block[{},
		AppendTo[fermions, field -> <|GaugeRep -> OptionValue[GaugeRep], FlavorIndices -> OptionValue[FlavorIndices]|>]; 
		SdelF/: SdelF[field, ind_, f1_] SdelF[Bar[field], ind_, f2_] = 
			Product[del[rep, f1, f2], {rep, OptionValue[GaugeRep]}]
			* Product[del[rep, f1, f2], {rep, OptionValue[FlavorIndices]}];
	];

(*Initiates a vector field*)	
CreateVector[name_, group_] :=
	Block[{},
		SdelV/: SdelV[name, ind_, f1_] SdelV[name, ind_, f2_] = del[group[adj], f1, f2];
	];