(*del[index,a,b] is the symbol for Kronecker delta \delta_{a,b} belong to the indices specified by index.*)
Clear[del];
del /: del[group_, rep_, a___, x_, b___] del[group_, rep_, c___, x_, d___] := del[group, rep, c, a, b, d];
del /: Power[del[group_, rep_, a_, b_], 2] := Dim[group, rep];
del /: Power[del[group_, rep_, a_, b_], 2] := Dim[group, rep];
del[group_, rep_, a_, a_] := Dim[group, rep];

Clear[TGen];
TGen /: del[group_, rep_, a___, x_, b___] TGen[group_, rep_, A_, c___, x_, d___] := TGen[group, rep, c, a, b, d];
TGen /: del[group_, adj_, a___, x_, b___] TGen[group_, rep_, x_, c__] := TGen[group, rep, a, b, c];
TGen /: TGen[group_, rep_, A_, a_, b_] TGen[group_, rep_, A_, b_, c_] := Casimir2[group, rep] del[group, rep, a, c];
TGen /: TGen[group_, rep_, A_, a_, b_] TGen[group_, rep_, B_, b_, a_] := TraceNormalization[group, rep] del[group, adj, A, B];
TGen[group_, rep_, A_, a_, a_] := 0;

SUGroup[group_, n_Integer] := 
	Block[{},
		Dim[group, rep] = n;
		Dim[group, adj] = n^2 - 1;
		TraceNormalization[group, fund] = 1/2;
		TraceNormalization[group, adj] = n;
		Casimir2[group, fund] = (n^2 - 1) / (2 n);
		Casimir2[group, adj] = n;
		(*Fierz identitiy*)
		TGen /: TGen[group, fund, A_, a_, b_] TGen[group, fund, A_, c_, d_] = TraceNormalization[group, fund] *
			(del[group, fund, a, d] del[group, fund, c, b] - del[group, fund, a, b] del[group, fund, c, d] / n); 
	];
