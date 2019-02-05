sig1 = {{0, 1}, {1, 0}};

(*Gauge coupling multiplied on generators*)
TfG2[A_, i_, j_] := Module[{B}, Tf[B, i, j] G2[A, B] // Expand]
TfG2t[A_, i_, j_] := Module[{B}, Tft[B, i, j] G2[A, B] // Expand]
TsG2[A_, a_, b_] := Module[{B}, Ts[B, a, b] G2[A, B] // Expand]
(*--------------------------------------*)
(*----------2-point structures----------*)
(*--------------------------------------*)

(*1-loop*)
C2G[A_, B_] := Module[{C1, C2, D1, D2},
		Ttimes[FGauge[A, C1, D1], G2[C1, C2], G2[D1, D2], FGauge[B, C2, D2]]
	]; 
S2F[A_, B_] := Module[{i, j},
		Tr[Tf[A, i, j].Tf[B, j, i]] // Expand
	];
S2S[A_, B_] := Module[{i, j},
		Ts[A, i, j] Ts[B, j, i] // Expand
	];
C2F[i_, j_] := Module[{A, k},
		TfG2[A, i, k].Tf[A, k, j] // Expand
	];
C2Ft[i_, j_] :=	Module[{A, k},
		TfG2t[A, i, k].Tft[A, k, j] // Expand
	];
C2S[a_, b_] := Module[{A, c},
	TsG2[A, a, c] Ts[A, c, b] // Expand
	];
Y2S[a_, b_] := Module[{i, j},
		Tr[y[a, i, j].yt[b, j, i]] // Expand
	];
Y2F[i_, j_] := Module[{a, k},
		y[a, i, k].yt[a, k, j] // Expand
	];
Y2Ft[i_, j_] := Module[{a, k},
		yt[a, i, k].y[a, k, j] // Expand
	];

(*----------2-loop----------*)
(*g^4*)
S2SC2S[A_, B_] := Module[{a, b, c},
		Ttimes[Ts[B, c, a], Ts[A, a, b], C2S[b, c]]
	];
S2FC2F[A_, B_] := Module[{i, j, k},
		Tr[Tdot[Tf[B, k, i], Tf[A, i, j], C2F[j, k]]]
	];
C2FC2G[i_, j_] := Module[{A, B, k},
		Tdot[TfG2[A, i, k], TfG2[B, k, j]] C2G[A, B] // Expand
	];
C2FS2F[i_, j_] := Module[{A, B, k},
		Tdot[TfG2[A, i, k], TfG2[B, k, j]] S2F[A, B] // Expand
	];
C2FS2S[i_, j_] := Module[{A, B,k},
		Tdot[TfG2[A, i, k], TfG2[B, k, j]] S2S[A, B] // Expand
	];
C2FC2Gt[i_, j_] := sig1.C2FC2G[i, j].sig1;
C2FS2Ft[i_, j_] := sig1.C2FS2F[i, j].sig1;
C2FS2St[i_, j_] := sig1.C2FS2S[i, j].sig1;
C2SC2G[a_, b_] := Module[{A, B, c},
		Ttimes[TsG2[A, a, c], TsG2[B, c, b], C2G[A, B]]
	];
C2SS2F[a_, b_] := Module[{A, B, c},
		Ttimes[TsG2[A, a, c], TsG2[B, c, b], S2F[A, B]]
	];
C2SS2S[a_, b_] := Module[{A, B, c},
		Ttimes[TsG2[A, a, c], TsG2[B, c, b], S2S[A, B]]
	];

(*g^2 y^2*)
Y2SC2F[a_, b_] := Module[{i, j, k},
		Tr @ Tdot[y[a, i, j], C2F[j, k], yt[b, k, i]]
	];
Y2FC2F[i_, j_] := Module[{a, k, l},
		Tdot[y[a, i, k], C2F[k, l], yt[a, l, j]]
	];
Y2FC2S[i_, j_] := Module[{a, b, k},
		Tdot[y[a, i, k], yt[b, k, j]] C2S[a, b] // Expand
	];
Y2FC2Ft[i_, j_] := sig1.Y2FC2F[i, j].sig1;
Y2FC2St[i_, j_] := sig1.Y2FC2S[i, j].sig1;
S2SY2S[A_, B_] := Module[{a, b, c},
		Ttimes[Ts[B, c, a], Ts[A, a, b], Y2S[b, c]]
	];
S2FY2F[A_, B_] := Module[{i, j, k},
		Tr @ Tdot[Tft[B, k, i], Tft[A, i, j], Y2F[j, k]]
	];

(*y^4*)
Y2SY2F[a_, b_] := Module[{k1, k2, k3},
		Tr @ Tdot[y[a, k1, k2], yt[b, k2, k3], Y2F[k3, k1]]
	];
Y4cS[a_, b_] := Module[{c, k1, k2, k3, k4},
		Tr @ Tdot[y[a, k1, k2], yt[c, k2, k3], y[b, k3, k4], yt[c, k4, k1]]
	];
Y4cF[i_, j_] := Module[{a, b, k1, k2, k3},
		Ttimes[y[a, i, k1], yt[b, k1, k2], y[a, k2, k3], yt[b, k3, j]]
	];
Y2FY2F[i_, j_] := Module[{a, k1, k2},
		Tdot[y[a, i, k1], Y2Ft[k1, k2], yt[a, k2, j]]
	];
Y2FY2S[i_, j_] := Module[{a, b, k1},
		Tdot[y[a, i, k1], yt[b, k1, j]] Y2S[a, b] // Expand
	];
Y4cFt[i_, j_] := sig1.Y4cF[i, j].sig1;
Y2FY2Ft[i_, j_] := sig1.Y2FY2F[i, j].sig1;
Y2FY2St[i_, j_] := sig1.Y2FY2S[i, j].sig1;

(*lambda^2*)
Lam2[a_, b_] := Module[{c, d, e},
		Lam[a, c, d, e] Lam[c, d, e, b] // Expand
	];

(*---------------------------------*)
(*----------Gauge tensors----------*)
(*---------------------------------*)
(*Gauge tensors at 1-loop order*)
GaugeTensors[1] := 
	Block[{bg, n},
		bg[1, 1] = C2G[A, B];
		bg[1, 2] = S2F[A, B];
		bg[1, 3] = S2S[A, B];
		Sum[bg[1, n], {n, 3}]
	];
	
	
	
	
	
	
	