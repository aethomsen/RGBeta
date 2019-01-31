(*Structure deltas*)
Clear[SdelS, SdelF, SdelV]
SdelS /: SdelS[field1_, ind_, f_] SdelS[field2_, ind_, g_] := 0; 
SdelF /: SdelF[field1_, ind_, f_] SdelF[field2_, ind_, g_] := 0;
SdelV /: SdelV[field1_, ind_, f_] SdelV[field2_, ind_, g_] := 0;

Options[CreateField] = {Type -> Scalar, SelfConjugate -> False, GaugeRepresentation -> <||>, FlavorIndices -> <||>};
CreateField[name_, OptionsPattern[] ] :=
	Module[{},0
		
	];