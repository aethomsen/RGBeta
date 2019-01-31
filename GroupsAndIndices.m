(*del[group, rep, a, b] is the symbol for Kronecker delta \delta_{a,b} belong to the indices specified by the group and representation.*)
Clear[del];
	del /: del[rep_, a___, x_, b___] del[rep_, c___, x_, d___] := del[rep, c, a, b, d];
	del /: Power[del[rep_, a_, b_], 2] := Dim[rep];
	del[rep_, a_, a_] := Dim[rep];

Clear[eps];
	eps /: del[rep_, fund_, a___, x_, b___] eps[rep_, c___, x_, d___] := eps[rep, c, a, b, d];
	eps /: eps[rep_ ,a_, b_] eps[rep_ ,b_, c_] := -Dim[rep, fund];
	eps /: eps[rep_ ,a_, b_] eps[rep_ ,c_, b_] := Dim[rep, fund];
	eps /: eps[rep_ ,b_, a_] eps[rep_ ,b_, c_] := Dim[rep, fund];
	eps /: Power[eps[rep_, a_, b_], 2] := Dim[rep];
	eps[rep_, a_, a_] := 0;

Clear[TGen];
	TGen /: del[rep_, a___, x_, b___] TGen[rep_, A_, c___, x_, d___] := TGen[rep, A, c, a, b, d];
	TGen /: del[group_[adj], A___, X_, B___] TGen[group_, rep_, X_, c__] := TGen[rep, A, B, c];
	TGen /: TGen[rep_, A_, a_, b_] TGen[rep_, A_, b_, c_] := Casimir2[rep] del[rep, a, c];
	TGen /: TGen[rep_, A_, a_, b_] TGen[rep_, B_, b_, a_] := TraceNormalization[rep] del[Head[rep][adj], A, B];
	TGen[rep_, A_, a_, a_] := 0;

(*Initiation for an SU(n) gauge group*)
SUGroup[group_, n_Integer] := 
	Module[{},
		(*Fundamental*)
		Dim[group[fund]] = n;
		TraceNormalization[group[fund]] = 1/2;
		Casimir2[group[fund]] = (n^2 - 1) / (2 n);
		(*Fierz identitiy*)
		TGen /: TGen[group[fund], A_, a_, b_] TGen[group[fund], A_, c_, d_] = TraceNormalization[group[fund]] *
			(del[group[fund], a, d] del[group[fund], c, b] - del[group[fund], a, b] del[group[fund], c, d] / n);
		
		(*Adjoint*)
		Dim[group[adj]] = n^2 - 1;
		TraceNormalization[group[adj]] = n;
		Casimir2[group[adj]] = n;
	];
