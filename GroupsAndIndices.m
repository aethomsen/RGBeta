(*----------Generic tensor properties----------*)
(*del[rep, a, b] is the symbol for Kronecker delta \delta_{a,b} belong to the indices specified by the representation.*)
Clear[del];
	del /: del[rep_, a___, x_, b___] del[rep_, c___, x_, d___] := del[rep, c, a, b, d];
	del /: Power[del[rep_, a_, b_], 2] := Dim[rep];
	del[rep_, a_, a_] := Dim[rep];

(*Default properties for anti-symmetric inavariants.*)
Clear[eps];
	eps /: del[rep_, fund_, a___, x_, b___] eps[rep_, c___, x_, d___] := eps[rep, c, a, b, d];
	eps /: eps[rep_ ,a_, b_] eps[rep_ ,b_, c_] := -Dim[rep];
	eps /: eps[rep_ ,a_, b_] eps[rep_ ,c_, b_] := Dim[rep];
	eps /: eps[rep_ ,b_, a_] eps[rep_ ,b_, c_] := Dim[rep];
	eps /: Power[eps[rep_, a_, b_], 2] := Dim[rep];
	eps[rep_, a_, a_] := 0;

(*Default group generator properties.*)
Clear[TGen];
	TGen /: del[rep_, a___, x_, b___] TGen[rep_, A_, c___, x_, d___] := TGen[rep, A, c, a, b, d];
	TGen /: del[group_[adj], A___, X_, B___] TGen[group_[rep_], X_, c__] := TGen[group[rep], A, B, c];
	TGen /: TGen[rep_, A_, a_, b_] TGen[rep_, A_, b_, c_] := Casimir2[rep] del[rep, a, c];
	TGen /: TGen[rep_, A_, a_, b_] TGen[rep_, B_, b_, a_] := TraceNormalization[rep] del[Head[rep][adj], A, B];
	TGen[rep_, A_, a_, a_] := 0;
	
(*Default structure constant properties.*)
CasimirSig[group_, a_, b_] := 
	Block[{common, ind1, ind2},
		common = Intersection[a, b];
		ind1 = Complement[a, common][[1]];
		ind2 = Complement[b, common][[1]];
		Signature[DeleteCases[a, ind1]] Signature[DeleteCases[b, ind2]] del[group[adj], ind1, ind2]
	];
Clear[fStruct];
	fStruct /: del[group_[adj], A___, X_, B___] fStruct[group_, C___, X_, D___] := fStruct[group, C, A, B, D];
	fStruct /: fStruct[group_, A___] fStruct[group_, B___] /; Length @ Intersection[List @A, List @B] === 2 := 
		Casimir2@ group[adj] CasimirSig[group, List@A, List@B];
	fStruct[group_, a___, x_, b___, x_, c___] := 0;

(*Matrix head for symbolic matrix manipulations*)
Clear[Trans];
	Trans[Trans[x_]] := x;
Clear[Matrix];
	Matrix /: del[ind_, a___, x_, b___] Matrix[m__][c___, ind_[x_], d___] := Matrix[m][c, ind[a, b], d];
	Matrix /: Matrix[m1__][ind_[b_]] Matrix[m2__][ind_[b_], c_] := Matrix[Sequence @@ Reverse[Trans /@ List@m1], m2][c];
	Matrix /: Matrix[m1__][a___, ind_[b_]] Matrix[m2__][ind_[b_], c___] := Matrix[m1, m2][a, c];
	Matrix /: Matrix[m1__][a_, ind_[b_]] Matrix[m2__][c_, ind_[b_]] := Matrix[m1, Sequence @@ Reverse[Trans /@ List@m2]][a, c];
	Matrix /: Matrix[m1__][ind_[b_], a_] Matrix[m2__][ind_[b_], c_] := Matrix[Sequence @@ Reverse[Trans /@ List@m1], m2][a, c];
	Matrix[m__][ind_[a_], ind_[a_]] := Tr @ Dot[m];
	Matrix[m__][] := m;


(*----------Gauge group definitions----------*)
(*Association with all information on the gauge groups: fields, couplings etc.*)
gaugeGroups = <||>;

(*Initialization for an SU(n) gauge group.*)
SUGroup[group_Symbol, n_Integer, OptionsPattern[{Field -> A[group], Coupling -> g[group]}] ] := 
	Block[{},		
		AppendTo[gaugeGroups, group -> <|Field -> OptionValue[Field], Coupling -> OptionValue[Coupling] |>];
		
		(*Fundamental*)
		Dim[group[fund]] = n;
		TraceNormalization[group[fund]] = 1/2;
		Casimir2[group[fund]] = (n^2 - 1) / (2 n);
		(*Fierz identitiy*)
		TGen /: TGen[group[fund], A_, a_, b_] TGen[group[fund], A_, c_, d_] = TraceNormalization[group[fund]] *
			(del[group[fund], a, d] del[group[fund], c, b] - del[group[fund], a, b] del[group[fund], c, d] / n);
		
		(*Anti-Fundamental*)
		del[group[afund], a_, b_] = del[group[fund], a, b];
		TGen[group[afund], A_, a_, b_] = - TGen[group[fund], A, b, a];
		
		(*Adjoint*)
		Dim[group[adj]] = n^2 - 1;
		TraceNormalization[group[adj]] = n;
		Casimir2[group[adj]] = n;
		
		(*Sets up the gauge fields and 2-point projection*)
		CreateVector[OptionValue[Field], group];
		GaugeProjection[group, A_, B_] := SdelV[OptionValue[Field], A, v1] SdelV[OptionValue[Field], B, v2] *
			del[group[adj], v1, v2] / Dim @ group[adj]; 
	];

(*Initialization for a U(1) gauge group.*)	
U1Group[group_Symbol, OptionsPattern[{Field -> A[group], Coupling -> g[group]}] ] := 
	Block[{},
		AppendTo[gaugeGroups, group -> <|Field -> OptionValue[Field], Coupling -> OptionValue[Coupling] |>];
		
		del[group[_],___] = 1;
		Dim[group[x_]] = x;
		TGen[group[x_],___] = x; 
		Dim[group[adj]] = 1; 
		fStruct[group, __] = 0; 
		
		(*Sets up the gauge fields and 2-point projection*)
		CreateVector[OptionValue[Field], group];
		GaugeProjection[group, A_, B_] := SdelV[OptionValue[Field], A, v1] SdelV[OptionValue[Field], B, v2];
	];	
	
	