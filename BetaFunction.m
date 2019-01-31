<< GroupsAndIndices`
<< FieldsAndCouplings`

Tf[A_, i_, j_] := Module[{ferm, rep, gRep, f1, f2, v}, 
	Sum[
		SdelF[ferm, i, f1] SdelF[Bar[ferm], j, f2] 
		* Sum[SdelV[gaugeGroups[Head @ gRep][Field], A, v] TGen[gRep, v, f1, f2] ,{gRep, fermions[ferm][GaugeRep]}] 
		* Product[del[rep, f1, f2], {rep, fermions[ferm][FlavorIndices]}]
	,{ferm, Keys @ fermions}]
];
