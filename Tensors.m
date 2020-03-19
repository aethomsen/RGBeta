(*
	Author: Anders Eller Thomsen 
	Released under the MIT license (see 'MIT_license.txt').
*)
Begin["Tensors`"]

sig1 = {{0, 1}, {1, 0}};

(*Gauge coupling multiplied on generators*)
TfG2[A_, i_, j_] := Module[{B}, Tferm[B, i, j] G2Matrix[A, B] // Expand]
TfG2t[A_, i_, j_] := Module[{B}, TfermTil[B, i, j] G2Matrix[A, B] // Expand]
TsG2[A_, a_, b_] := Module[{B}, Tscal[B, a, b] G2Matrix[A, B] // Expand]

(*######################################*)
(*----------2-point structures----------*)
(*######################################*)
(*1-loop*)
C2G[A_, B_] := Module[{C1, C2, D1, D2},
		Ttimes[FGauge[A, C1, D1], G2Matrix[C1, C2], G2Matrix[D1, D2], FGauge[B, C2, D2]]
	]; 
S2F[A_, B_] := Module[{i, j},
		Tr[Tferm[A, i, j].Tferm[B, j, i]] // Expand
	];
S2S[A_, B_] := Module[{i, j},
		Tscal[A, i, j] Tscal[B, j, i] // Expand
	];
C2F[i_, j_] := Module[{A, k},
		TfG2[A, i, k].Tferm[A, k, j] // Expand
	];
C2Ft[i_, j_] :=	Module[{A, k},
		TfG2t[A, i, k].TfermTil[A, k, j] // Expand
	];
C2S[a_, b_] := Module[{A, c},
	TsG2[A, a, c] Tscal[A, c, b] // Expand
	];
Y2S[a_, b_] := Module[{i, j},
		Tr[Yuk[a, i, j].YukTil[b, j, i]] // Expand
	];
Y2F[i_, j_] := Module[{a, k},
		Yuk[a, i, k].YukTil[a, k, j] // Expand
	];
Y2Ft[i_, j_] := Module[{a, k},
		YukTil[a, i, k].Yuk[a, k, j] // Expand
	];

(*----------2-loop----------*)
(*g^4*)
S2SC2S[A_, B_] := Module[{a, b, c},
		Ttimes[Tscal[B, c, a], Tscal[A, a, b], C2S[b, c]]
	];
S2FC2F[A_, B_] := Module[{i, j, k},
		Tr[Tdot[Tferm[B, k, i], Tferm[A, i, j], C2F[j, k]]]
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
		Tr @ Tdot[Yuk[a, i, j], C2F[j, k], YukTil[b, k, i]]
	];
Y2FC2F[i_, j_] := Module[{a, k, l},
		Tdot[Yuk[a, i, k], C2F[k, l], YukTil[a, l, j]]
	];
Y2FC2S[i_, j_] := Module[{a, b, k},
		Tdot[Yuk[a, i, k], YukTil[b, k, j]] C2S[a, b] // Expand
	];
Y2FC2Ft[i_, j_] := sig1.Y2FC2F[i, j].sig1;
Y2FC2St[i_, j_] := sig1.Y2FC2S[i, j].sig1;
S2SY2S[A_, B_] := Module[{a, b, c},
		Ttimes[Tscal[B, c, a], Tscal[A, a, b], Y2S[b, c]]
	];
S2FY2F[A_, B_] := Module[{i, j, k},
		Tr @ Tdot[TfermTil[B, k, i], TfermTil[A, i, j], Y2F[j, k]]
	];

(*y^4*)
Y2SY2F[a_, b_] := Module[{k1, k2, k3},
		Tr @ Tdot[Yuk[a, k1, k2], YukTil[b, k2, k3], Y2F[k3, k1]]
	];
Y4cS[a_, b_] := Module[{c, k1, k2, k3, k4},
		Tr @ Tdot[Yuk[a, k1, k2], YukTil[c, k2, k3], Yuk[b, k3, k4], YukTil[c, k4, k1]]
	];
Y4cF[i_, j_] := Module[{a, b, k1, k2, k3},
		Ttimes[Yuk[a, i, k1], YukTil[b, k1, k2], Yuk[a, k2, k3], YukTil[b, k3, j]]
	];
Y2FY2F[i_, j_] := Module[{a, k1, k2},
		Tdot[Yuk[a, i, k1], Y2Ft[k1, k2], YukTil[a, k2, j]]
	];
Y2FY2S[i_, j_] := Module[{a, b, k1},
		Tdot[Yuk[a, i, k1], YukTil[b, k1, j]] Y2S[a, b] // Expand
	];
Y4cFt[i_, j_] := sig1.Y4cF[i, j].sig1;
Y2FY2Ft[i_, j_] := sig1.Y2FY2F[i, j].sig1;
Y2FY2St[i_, j_] := sig1.Y2FY2S[i, j].sig1;

(*lambda^2*)
Lam2[a_, b_] := Module[{c, d, e},
		Lam[a, c, d, e] Lam[c, d, e, b] // Expand
	];

(*1-loop vertices*)
FourGLam[a_, b_, c_, d_] := Module[{A1, A2, b1, b2},
	Ttimes[TsG2[A1, a, b1], TsG2[A2, b1, b], Tscal[A1, c, b2], Tscal[A2, b2, d]]
]

(*#################################*)
(*----------Gauge tensors----------*)
(*#################################*)

(*Evaluates the tensors involved in the gauge beta function*)	
GaugeTensors[loop_Integer] := GaugeTensors[loop] =
	Module[{bg, n, i, j, k, l, a, b, c, d, C1, C2, C3, C4, e, f},
		Switch[loop
		,0,
			Return[- \[Epsilon] / 2 G2Matrix[$A, $B, -1]];
		,1,
			bg[1, 1] := C2G[$A, $B];
			bg[1, 2] := S2F[$A, $B];
			bg[1, 3] := S2S[$A, $B];
		,2,
			bg[2, 1] := S2FC2F[$A, $B];
			bg[2, 2] := S2SC2S[$A, $B];
			bg[2, 3] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], C2G[C2, $B]];
			bg[2, 4] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], S2F[C2, $B]];
			bg[2, 5] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], S2S[C2, $B]];
			bg[2, 6] := S2FY2F[$A, $B];
			bg[2, 7] := 0 (*S2SY2S[$A, $B]*);
		,3,
		(*y^0 terms*)
			bg[3, 1] := Tr@Tdot[Tferm[$B, l, i], Tferm[$A, i, j], C2F[j, k], C2F[k, l]];
			bg[3, 2] := Ttimes[Tscal[$B, d, a], Tscal[$A, a, b], C2S[b, c], C2S[c, d]];
			bg[3, 3] := Tr@Tdot[Tferm[$B, k, i], Tferm[$A, i, j], C2FC2G[j, k]];
			bg[3, 4] := Ttimes[Tscal[$B, c, a], Tscal[$A, a, b], C2SC2G[b, c]];
			bg[3, 5] := Tr@Tdot[Tferm[$B, k, i], Tferm[$A, i, j], C2FS2F[j, k]];
			bg[3, 6] := Ttimes[Tscal[$B, c, a], Tscal[$A, a, b], C2SS2F[b, c]];
			bg[3, 7] := Tr@Tdot[Tferm[$B, k, i], Tferm[$A, i, j], C2FS2S[j, k]];
			bg[3, 8] := Ttimes[Tscal[$B, c, a], Tscal[$A, a, b], C2SS2S[b, c]];		
			bg[3, 9] := Ttimes[S2FC2F[$A, C1], G2Matrix[C1, C2], C2G[C2, $B]];
			bg[3, 10] := Ttimes[S2SC2S[$A, C1], G2Matrix[C1, C2], C2G[C2, $B]];
			bg[3, 11] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], C2G[C2, C3], G2Matrix[C3, C4], C2G[C4, $B]];
			bg[3, 12] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], S2F[C2, C3], G2Matrix[C3, C4], S2F[C4, $B]];
			bg[3, 13] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], S2S[C2, C3], G2Matrix[C3, C4], S2S[C4, $B]];
			bg[3, 14] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], C2G[C2, C3], G2Matrix[C3, C4], S2F[C4, $B]];
			bg[3, 15] := Ttimes[C2G[$A, C1], G2Matrix[C1, C2], C2G[C2, C3], G2Matrix[C3, C4], S2S[C4, $B]];
			bg[3, 16] := Ttimes[S2F[$A, C1], G2Matrix[C1, C2], C2G[C2, C3], G2Matrix[C3, C4], S2S[C4, $B]];
			bg[3, 17] := Ttimes[Tscal[$A, a, b], TsG2[C1, b, c], Lam[a, c, d, f], Tscal[$B, d, e], Tscal[C1, e, f]];
			bg[3, 18] := Ttimes[Tscal[$A, a, b], Lam2[b, c], Tscal[$B, c, a]];
		(*y^2 terms*)
			bg[3, 19] := Tr @ Tdot[Tferm[$B, l, i], Tferm[$A, i, j], Y2Ft[j, k], C2F[k, l]];
			bg[3, 20] := Tr @ Tdot[Tferm[$B, k, i], Tferm[$A, i, j], Y2FC2Ft[j, k]];
			bg[3, 21] := Ttimes[Tscal[$B, c, a], Tscal[$A, a, b], Y2SC2F[b, c]];
			bg[3, 22] := Tr @ Tdot[Tferm[$B, k, i], Tferm[$A, i, j], Y2FC2St[j, k]];
			bg[3, 23] := Ttimes[Tscal[$B, d, a], Tscal[$A, a, b], Y2S[b, c], C2S[c, d]];
			bg[3, 24] := Ttimes[S2FY2F[$A, C1], G2Matrix[C1, C2], C2G[C2, $B]];
			bg[3, 25] := Ttimes[S2SY2S[$A, C1], G2Matrix[C1, C2], C2G[C2, $B]];
		(*y^4 terms*)
			bg[3, 26] := Tr @ Tdot[Yuk[a, f, i], Tferm[$A, i, j], YukTil[a, j, k] Yuk[b, k, l], Tferm[$B, l, e], YukTil[b, e, f]];
			bg[3, 27] := Ttimes[Tscal[$B, c, a], Tscal[$A, a, b], Y4cS[b, c]];
			bg[3, 28] := Tr @ Tdot[Tferm[$B, k, i], Tferm[$A, i, j], Y4cFt[j, k]];
			bg[3, 29] := Tr @ Tdot[Tferm[$B, k, i], Tferm[$A, i, j], Y2FY2St[j, k]];
			bg[3, 30] := Ttimes[Tscal[$B, c, a], Tscal[$A, a, b], Y2SY2F[b, c]];
			bg[3, 31] := Tr @ Tdot[Tferm[$B, l, i], Tferm[$A, i, j], Y2Ft[j, k], Y2Ft[k, l]];
			bg[3, 32] := Tr @ Tdot[Tferm[$B, k, i], Tferm[$A, i, j], Y2FY2Ft[j, k]];
			bg[3, 33] := Ttimes[Tscal[$B, d, a], Tscal[$A, a, b], Y2S[b, c], Y2S[c, d]];
		];
		
		With[{diagrams = {3, 7, 33}[[loop]]},
			Monitor[
				Return[Sum[Bcoef[1, loop, n] bg[loop, n], {n, diagrams}] /. $gaugeCoefficients]
			,StringForm["Evaluating term `` / ``", n, diagrams]]
		]
	];


(*##################################*)
(*----------Yukawa tensors----------*)
(*##################################*)

(*Evaluates the tensors involved in the yukawa beta function*)
YukawaTensors[loop_Integer] := YukawaTensors[loop] = 
	Module[{bYuk, n, b, c, d, k, k1, k2, k3, k4, A},
		Switch[loop
		,0,
			Return[- \[Epsilon] / 2 Yuk[$a, $i, $j] ];
		,1,
			bYuk[1, 1] := 0 (*Yuk[b, $i, $j] C2S[$a, b] // Expand*);
			bYuk[1, 2] := Tdot[C2Ft[$i, k1], Yuk[$a, k1, $j]] // Sym[$i, $j];
			bYuk[1, 3] := Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2], Yuk[b, k2, $j]];
			bYuk[1, 4] := Tdot[Y2F[$i, k1], Yuk[$a, k1, $j]] // Sym[$i, $j];
			bYuk[1, 5] := Yuk[b, $i, $j] Y2S[b, $a] // Expand;
		,2,
		(* y^1 terms*)
			bYuk[2, 1] := Ttimes[Yuk[c, $i, $j], C2S[c, b], C2S[b, $a]];
			bYuk[2, 2] := Tdot[C2Ft[$i, k], Yuk[b, k, $j]] C2S[b, $a] // Expand // Sym[$i, $j];
			bYuk[2, 3] := 0 (*Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2F[k2, $j]]*);
			bYuk[2, 4] := Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Yuk[$a, k2, $j]] // Sym[$i, $j];
			bYuk[2, 5] := Yuk[b, $i, $j] C2SC2G[b, $a] // Expand;
			bYuk[2, 6] := Yuk[b, $i, $j] C2SS2S[b, $a] // Expand;
			bYuk[2, 7] := Yuk[b, $i, $j] C2SS2F[b, $a] // Expand;
			bYuk[2, 8] := C2FC2Gt[$i, k]. Yuk[$a, k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 9] := C2FS2St[$i, k]. Yuk[$a, k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 10] := C2FS2Ft[$i, k]. Yuk[$a, k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 11] := Yuk[b, $i, $j] Lam2[b, $a] // Expand;
		(*y^3 terms*)
			bYuk[2, 12] := Tdot[TfermTil[A, $i, k1], Yuk[$a, k1, k2], YukTil[b, k2, k3], TfG2t[A, k3, k4], Yuk[b, k4, $j]] // Sym[$i, $j];
			bYuk[2, 13] := 0 (*Tdot[Y2F[$i, k1], TfermTil[A, k1, k2], Yuk[$a, k2, k3], TfG2[A, k3, $j]] // Sym[$i, $j]*);
			bYuk[2, 14] := Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2], Yuk[c, k2, $j]] C2S[b, c] // Expand;
			bYuk[2, 15] := Tdot[Yuk[b, $i, k1], YukTil[c, k1, k2], Yuk[b, k2, $j]] C2S[c, $a] // Expand;
			bYuk[2, 16] := Tdot[Yuk[b, $i, k1], C2F[k1, k2], YukTil[$a, k2, k3], Yuk[b, k3, $j]] // Sym[$i, $j];
			bYuk[2, 17] := Tdot[C2Ft[$i, k1], Yuk[b, k1, k2], YukTil[$a, k2, k3], Yuk[b, k3, $j]] // Sym[$i, $j];
			bYuk[2, 18] := Yuk[b, $i, $j] Y2SC2F[b, $a] // Expand;
			bYuk[2, 19] := Yuk[$a, $i, k].Y2FC2St[k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 20] := Yuk[$a, $i, k].Y2FC2Ft[k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 21] := Tdot[Yuk[$a, $i, k1], Y2Ft[k1, k2], C2F[k2, $j]] // Sym[$i, $j];
			bYuk[2, 22] := 0 (*Ttimes[Yuk[c, $i, $j], Y2S[c, b], C2S[b, $a]]*);
			bYuk[2, 23] := Tdot[Yuk[b, $i, k1], YukTil[c, k1, k2], Yuk[d, k2, $j]] Lam[$a, b, c, d] // Expand;
		(*y^5 terms*)
			bYuk[2, 24] := Tdot[Yuk[b, $i, k1], YukTil[c, k1, k2], Yuk[$a, k2, k3], YukTil[b, k3, k4], Yuk[c, k4, $j]];
			bYuk[2, 25] := 0 (*Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2], Yuk[c, k2, k3], YukTil[b, k3, k4], Yuk[c, k4, $j]] // Sym[$i, $j]*);
			bYuk[2, 26] := Tdot[Yuk[c, $i, k1], YukTil[b, k1, k2], Yuk[$a, k2, k3], YukTil[b, k3, k4], Yuk[c, k4, $j]];
			bYuk[2, 27] := 0 (*Tdot[Yuk[$a, $i, k], Y4cFt[k, $j]] // Sym[$i, $j]*);
			bYuk[2, 28] := Yuk[b, $i, $j] Y4cS[b, $a] // Expand;
			bYuk[2, 29] := Tdot[Yuk[b, $i, k1], Y2Ft[k1, k2], YukTil[$a, k2, k3], Yuk[b, k3, $j]] // Sym[$i, $j];
			bYuk[2, 30] := Tdot[Yuk[$a, $i, k], Y2FY2Ft[k, $j]] // Sym[$i, $j];
			bYuk[2, 31] := Yuk[b, $i, $j] Y2SY2F[b, $a] // Expand;
			bYuk[2, 32] := Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2], Yuk[c, k2, $j]] Y2S[b, c] // Expand;
			bYuk[2, 33] := Tdot[Yuk[$a, $i, k], Y2FY2St[k, $j]] // Sym[$i, $j];
		];
		
		With[{diagrams = {5, 33}[[loop]]},
			Monitor[
				Return[Sum[Bcoef[2, loop, n] bYuk[loop, n], {n, diagrams}] /. $yukawaCoefficients]
			,StringForm["Evaluating term `` / ``", n, diagrams]]
		]
	];
	
(*Evaluates the tensors involved in the Fermion mass beta function*)
FermionMassTensors[loop_Integer] := FermionMassTensors[loop] = 
	Module[{bYuk, n, b, c, d, k, k1, k2, k3, k4, A},
		Switch[loop
		,0,
			Return[0];
		,1,
			bYuk[1, 1] := 0;
			bYuk[1, 2] := Tdot[C2Ft[$i, k1], Yuk[$a, k1, $j, True]] // Sym[$i, $j];
			bYuk[1, 3] := Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2, True], Yuk[b, k2, $j]];
			bYuk[1, 4] := Tdot[Y2F[$i, k1], Yuk[$a, k1, $j, True]] // Sym[$i, $j];
			bYuk[1, 5] := 0;
		,2,
		(* y^1 terms*)
			bYuk[2, 1] := 0;
			bYuk[2, 2] := 0;
			bYuk[2, 3] := 0 (*Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2F[k2, $j]]*);
			bYuk[2, 4] := Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Yuk[$a, k2, $j, True]] // Sym[$i, $j];
			bYuk[2, 5] := 0;
			bYuk[2, 6] := 0;
			bYuk[2, 7] := 0;
			bYuk[2, 8] := C2FC2Gt[$i, k]. Yuk[$a, k, $j, True] // Expand // Sym[$i, $j];
			bYuk[2, 9] := C2FS2St[$i, k]. Yuk[$a, k, $j, True] // Expand // Sym[$i, $j];
			bYuk[2, 10] := C2FS2Ft[$i, k]. Yuk[$a, k, $j, True] // Expand // Sym[$i, $j];
			bYuk[2, 11] := 0;
		(*y^3 terms*)
			bYuk[2, 12] := Tdot[TfermTil[A, $i, k1], Yuk[$a, k1, k2, True], YukTil[b, k2, k3], TfG2t[A, k3, k4], Yuk[b, k4, $j]] // Sym[$i, $j];
			bYuk[2, 13] := 0 (*Tdot[Y2F[$i, k1], TfermTil[A, k1, k2], Yuk[$a, k2, k3], TfG2[A, k3, $j]] // Sym[$i, $j]*);
			bYuk[2, 14] := Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2, True], Yuk[c, k2, $j]] C2S[b, c] // Expand;
			bYuk[2, 15] := 0;
			bYuk[2, 16] := Tdot[Yuk[b, $i, k1], C2F[k1, k2], YukTil[$a, k2, k3, True], Yuk[b, k3, $j]] // Sym[$i, $j];
			bYuk[2, 17] := Tdot[C2Ft[$i, k1], Yuk[b, k1, k2], YukTil[$a, k2, k3, True], Yuk[b, k3, $j]] // Sym[$i, $j];
			bYuk[2, 18] := 0;
			bYuk[2, 19] := Yuk[$a, $i, k, True].Y2FC2St[k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 20] := Yuk[$a, $i, k, True].Y2FC2Ft[k, $j] // Expand // Sym[$i, $j];
			bYuk[2, 21] := Tdot[Yuk[$a, $i, k1, True], Y2Ft[k1, k2], C2F[k2, $j]] // Sym[$i, $j];
			bYuk[2, 22] := 0;
			bYuk[2, 23] := Tdot[Yuk[b, $i, k1], YukTil[c, k1, k2], Yuk[d, k2, $j]] Lam[$a, b, c, d, True] // Expand;
		(*y^5 terms*)
			bYuk[2, 24] := Tdot[Yuk[b, $i, k1], YukTil[c, k1, k2], Yuk[$a, k2, k3, True], YukTil[b, k3, k4], Yuk[c, k4, $j]];
			bYuk[2, 25] := 0 (*Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2], Yuk[c, k2, k3], YukTil[b, k3, k4], Yuk[c, k4, $j]] // Sym[$i, $j]*);
			bYuk[2, 26] := Tdot[Yuk[c, $i, k1], YukTil[b, k1, k2], Yuk[$a, k2, k3, True], YukTil[b, k3, k4], Yuk[c, k4, $j]];
			bYuk[2, 27] := 0 (*Tdot[Yuk[$a, $i, k], Y4cFt[k, $j]] // Sym[$i, $j]*);
			bYuk[2, 28] := 0;
			bYuk[2, 29] := Tdot[Yuk[b, $i, k1], Y2Ft[k1, k2], YukTil[$a, k2, k3, True], Yuk[b, k3, $j]] // Sym[$i, $j];
			bYuk[2, 30] := Tdot[Yuk[$a, $i, k, True], Y2FY2Ft[k, $j]] // Sym[$i, $j];
			bYuk[2, 31] := 0;
			bYuk[2, 32] := Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2, True], Yuk[c, k2, $j]] Y2S[b, c] // Expand;
			bYuk[2, 33] := Tdot[Yuk[$a, $i, k, True], Y2FY2St[k, $j]] // Sym[$i, $j];
		];
		
		With[{diagrams = {5, 33}[[loop]]},
			Monitor[
				Return[Sum[Bcoef[2, loop, n] bYuk[loop, n], {n, diagrams}] /. $yukawaCoefficients]
			,StringForm["Evaluating term `` / ``", n, diagrams]]
		]
	];


(*##################################*)
(*----------Quartic tensors----------*)
(*##################################*)

(*Evaluates the tensors involved in the quartic beta function*)
QuarticTensors[loop_Integer] := QuarticTensors[loop] = 
	Module[{bl, n, b1, b2, b3, b4, A1, A2, A3, k1, k2, k3, k4, k5, k6},
		Switch[loop
		,0,
			Return[- \[Epsilon] Lam[$a, $b, $c, $d] ];
		,1,
			bl[1, 1] := FourGLam[$a, $b, $c, $d] // Sym[$b, $c, $d];
			bl[1, 2] := C2S[$a, b1] Lam[b1, $b, $c, $d] // Expand // Sym4[$a, $b, $c, $d];
			bl[1, 3] := Lam[$a, $b, b1, b2] Lam[b1, b2, $c, $d] // Expand // Sym[$b, $c, $d];
			bl[1, 4] := Y2S[$a, b1] Lam[b1, $b, $c, $d] // Expand // Sym4[$a, $b, $c, $d];
			bl[1, 5] := Tr @ Tdot[Yuk[$a, k4, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4]] // Sym[$b, $c, $d];
		,2,		
		(* y^0, lam^0 terms*)
			bl[2, 1] := Ttimes[Tscal[A1, $a, b1], FourGLam[b1, b2, $b, $c], TsG2[A1, b2, $d]] // Sym[$a, $b, $c, $d]; 	
			bl[2, 2] := Ttimes[C2S[$a, b1], FourGLam[b1, $b, $c, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 3] := Ttimes[TsG2[A1, $a, b1], TsG2[A2, b1, $b], C2G[A2, A3], Tscal[A1, $c, b2], TsG2[A3, b2, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 4] := Ttimes[TsG2[A1, $a, b1], TsG2[A2, b1, $b], S2S[A2, A3], Tscal[A1, $c, b2], TsG2[A3, b2, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 5] := Ttimes[TsG2[A1, $a, b1], TsG2[A2, b1, $b], S2F[A2, A3], Tscal[A1, $c, b2], TsG2[A3, b2, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 6] := Ttimes[FourGLam[$a, b1, $d, b2], Lam[b1, b2, $b, $c]] // Sym[$a, $b, $c, $d];
			bl[2, 7] := Ttimes[FourGLam[$a, $b, b1, b2], Lam[b1, b2, $c, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 8] := 0 (*Ttimes[Lam[$a, $b, b1, b2], C2S[b1, $c], C2S[b2, $d]] // Sym[$a, $b, $c, $d]*);
			bl[2, 9] := Ttimes[C2S[$a, b1], C2S[b1, b2], Lam[b2, $b, $c, $d] ]  // Sym4[$a, $b, $c, $d];
			bl[2, 10] := Ttimes[C2SC2G[$a, b1], Lam[b1, $b, $c, $d] ]  // Sym4[$a, $b, $c, $d];
			bl[2, 11] := Ttimes[C2SS2S[$a, b1], Lam[b1, $b, $c, $d] ]  // Sym4[$a, $b, $c, $d];
			bl[2, 12] := Ttimes[C2SS2F[$a, b1], Lam[b1, $b, $c, $d] ]  // Sym4[$a, $b, $c, $d];
		(* y^0, lam^{n>0} terms*)
			bl[2, 13] := Ttimes[Tscal[A1, $a, b1], Lam[b1, b2, b3, b4], TsG2[A1, $b, b2], Lam[b3, b4, $c, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 14] := Ttimes[Lam[$a, $b, b1, b2], C2S[b2, b3], Lam[b3, b1, $c, $d]] // Sym[$b, $c, $d];
			bl[2, 15] := Ttimes[C2S[$a, b1], Lam[b1, $b, b2, b3], Lam[b2, b3, $c, $d]]  // Sym[$a, $b, $c, $d];
			bl[2, 16] := Ttimes[Lam2[$a, b1], Lam[b1,$b, $c, $d]]  // Sym4[$a, $b, $c, $d];
			bl[2, 17] := Ttimes[Lam[$a, b1, b2, b3], Lam[$b, b4, b2, b3], Lam[b1, b4, $c, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 18] := 0 (*Ttimes[Lam[$a, $b, b1, b2], Lam[b1, b2, b3, b4], Lam[b3, b4, $c, $d]] // Sym[$b, $c, $d]*);
		(*y^2 terms*)
			bl[2, 19] := Ttimes[Tr @ Tdot[Yuk[$a, k4, k1], Tferm[A1, k1, k2], Tferm[A2, k2, k3], YukTil[$b, k3, k4]], 
				TsG2[A1, $c, b1], TsG2[A2, b1, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 20] := Ttimes[Y2S[$a, b1], FourGLam[b1, $b, $c, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 21] := Ttimes[Y2SC2F[$a, b1], Lam[b1, $b, $c, $d]] // Sym4[$a, $b, $c, $d];
			bl[2, 22] := 0 (*Ttimes[C2S[$a, b1], Y2S[b1, b2], Lam[b2, $b, $c, $d] ]  // Sym4[$a, $b, $c, $d]*);
			bl[2, 23] := Ttimes[Lam[$a, $b, b1, b2], Y2S[b2, b3], Lam[b3, b1, $c, $d]] // Sym[$b, $c, $d];
		(*y^4 terms*)
			bl[2, 24] := 0 (*Tr @ Tdot[Yuk[$a, k6, k1], Tferm[A1, k1, k2], YukTil[$b, k2, k3], Yuk[$c, k3, k4], 
				TfG2[A1, k4, k5], YukTil[$d, k5, k6]] // Sym[$b, $c, $d]*);
			bl[2, 25] := Tr @ Tdot[Yuk[b1, k4, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4]] 
				* C2S[$a, b1] // Expand // Sym[$a, $b, $c, $d];
			bl[2, 26] := Tr @ Tdot[Yuk[$a, k5, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4], 
				C2Ft[k4, k5]] // Sym[$a, $b, $c, $d];
			bl[2, 27] := Tr @ Tdot[Yuk[$a, k4, k1], YukTil[b1, k1, k2], Yuk[$b, k2, k3], YukTil[b2, k3, k4]] 
				*Lam[b1, b2, $c, $d] // Sym[$a, $b, $c, $d];
			bl[2, 28] := 0 (*Tr @ Tdot[Yuk[$a, k4, k1], YukTil[$b, k1, k2], Yuk[b1, k2, k3], YukTil[b2, k3, k4]] 
				*Lam[b1, b2, $c, $d] // Sym[$a, $b, $c, $d]*);
			bl[2, 29] := Y4cS[$a, b1] Lam[b1, $b, $c, $d] // Sym4[$a, $b, $c, $d];
			bl[2, 30] := Y2SY2F[$a, b1] Lam[b1, $b, $c, $d] // Sym4[$a, $b, $c, $d];
		(*y^6 terms*)
			bl[2, 31] := Tr @ Tdot[Yuk[b1, k6, k1], YukTil[$a, k1, k2], Yuk[b1, k2, k3], YukTil[$b, k3, k4], 
				Yuk[$c, k4, k5], YukTil[$d, k5, k6]] // Sym[$a, $b, $c, $d];
			bl[2, 32] := Tr @ Tdot[Yuk[$a, k6, k1], YukTil[b1, k1, k2], Yuk[$b, k2, k3], YukTil[$c, k3, k4], 
				Yuk[b1, k4, k5], YukTil[$d, k5, k6]] // Sym[$b, $c, $d];
			bl[2, 33] := Tr @ Tdot[Yuk[$a, k5, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4], 
				Y2F[k4, k5]] // Sym[$a, $b, $c, $d];
		];
		
		With[{diagrams = {5, 33}[[loop]]},
			Monitor[
				Return[Sum[Bcoef[3, loop, n] bl[loop, n], {n, diagrams}] /. $quarticCoefficients]
			,StringForm["Evaluating term `` / ``", n, diagrams]]
		]
	];

(*Evaluates the tensors involved in the scalar mass and trillinear beta function*)
ScalarMassiveTensors[loop_Integer] := ScalarMassiveTensors[loop] = 
	Module[{bl, n, b1, b2, b3, b4, A1, A2, k1, k2, k3, k4, k5, k6},
		Switch[loop
		,0,
			Return[- \[Epsilon]/2 Lam[$a, $b, $c, $d, True] ];
		,1,
			bl[1, 1] := 0;
			bl[1, 2] := C2S[$a, b1] Lam[b1, $b, $c, $d, True] // Expand // Sym4[$a, $b, $c, $d];
			bl[1, 3] := Lam[$a, $b, b1, b2, True] Lam[b1, b2, $c, $d, True] // Expand // Sym[$b, $c, $d];
			bl[1, 4] := Y2S[$a, b1] Lam[b1, $b, $c, $d, True] // Expand // Sym4[$a, $b, $c, $d];
			bl[1, 5] := Tr @ Tdot[Yuk[$a, k4, k1, True], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], 
				YukTil[$d, k3, k4, True]] // Sym[$b, $c, $d];
		,2,		
		(* y^0, lam^0 terms*)
			bl[2, 1] := 0; 	
			bl[2, 2] := 0;
			bl[2, 3] := 0;
			bl[2, 4] := 0;
			bl[2, 5] := 0;
			bl[2, 6] := Ttimes[FourGLam[$a, b1, $d, b2], Lam[b1, b2, $b, $c, True]] // Sym[$a, $b, $c, $d];
			bl[2, 7] := Ttimes[FourGLam[$a, $b, b1, b2], Lam[b1, b2, $c, $d, True]] // Sym[$a, $b, $c, $d];
			bl[2, 8] := 0 (*Ttimes[Lam[$a, $b, b1, b2], C2S[b1, $c], C2S[b2, $d]] // Sym[$a, $b, $c, $d]*);
			bl[2, 9] := Ttimes[C2S[$a, b1], C2S[b1, b2], Lam[b2, $b, $c, $d, True] ]  // Sym4[$a, $b, $c, $d];
			bl[2, 10] := Ttimes[C2SC2G[$a, b1], Lam[b1, $b, $c, $d, True] ]  // Sym4[$a, $b, $c, $d];
			bl[2, 11] := Ttimes[C2SS2S[$a, b1], Lam[b1, $b, $c, $d, True] ]  // Sym4[$a, $b, $c, $d];
			bl[2, 12] := Ttimes[C2SS2F[$a, b1], Lam[b1, $b, $c, $d, True] ]  // Sym4[$a, $b, $c, $d];
		(* y^0, lam^{n>0} terms*)
			bl[2, 13] := Ttimes[Tscal[A1, $a, b1], Lam[b1, b2, b3, b4], TsG2[A1, $b, b2], Lam[b3, b4, $c, $d, True]] // Sym[$a, $b, $c, $d];
			bl[2, 14] := Ttimes[Lam[$a, $b, b1, b2, True], C2S[b2, b3], Lam[b3, b1, $c, $d, True]] // Sym[$b, $c, $d];
			bl[2, 15] := Ttimes[C2S[$a, b1], Lam[b1, $b, b2, b3, True], Lam[b2, b3, $c, $d, True]]  // Sym[$a, $b, $c, $d];
			bl[2, 16] := Ttimes[Lam2[$a, b1], Lam[b1,$b, $c, $d, True]]  // Sym4[$a, $b, $c, $d];
			bl[2, 17] := Ttimes[Lam[$a, b1, b2, b3, True], Lam[$b, b4, b2, b3, True], Lam[b1, b4, $c, $d, True]] // Sym[$a, $b, $c, $d];
			bl[2, 18] := 0 (*Ttimes[Lam[$a, $b, b1, b2], Lam[b1, b2, b3, b4], Lam[b3, b4, $c, $d]] // Sym[$b, $c, $d]*);
		(*y^2 terms*)
			bl[2, 19] := Ttimes[Tr @ Tdot[Yuk[$a, k4, k1, True], Tferm[A1, k1, k2], Tferm[A2, k2, k3], YukTil[$b, k3, k4, True]], 
				TsG2[A1, $c, b1], TsG2[A2, b1, $d]] // Sym[$a, $b, $c, $d];
			bl[2, 20] := 0;
			bl[2, 21] := Ttimes[Y2SC2F[$a, b1], Lam[b1, $b, $c, $d, True]] // Sym4[$a, $b, $c, $d];
			bl[2, 22] := 0 (*Ttimes[C2S[$a, b1], Y2S[b1, b2], Lam[b2, $b, $c, $d] ]  // Sym4[$a, $b, $c, $d]*);
			bl[2, 23] := Ttimes[Lam[$a, $b, b1, b2, True], Y2S[b2, b3], Lam[b3, b1, $c, $d, True]] // Sym[$b, $c, $d];
		(*y^4 terms*)
			bl[2, 24] := 0 (*Tr @ Tdot[Yuk[$a, k6, k1], Tferm[A1, k1, k2], YukTil[$b, k2, k3], Yuk[$c, k3, k4], 
				TfG2[A1, k4, k5], YukTil[$d, k5, k6]] // Sym[$b, $c, $d]*);
			bl[2, 25] := Tr @ Tdot[Yuk[b1, k4, k1], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True]] 
				* C2S[$a, b1] // Expand // Sym[$a, $b, $c, $d];
			bl[2, 26] := Tr @ Tdot[Yuk[$a, k5, k1, True], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True], 
				C2Ft[k4, k5]] // Sym[$a, $b, $c, $d];
			bl[2, 27] := Tr @ Tdot[Yuk[$a, k4, k1, True], YukTil[b1, k1, k2], Yuk[$b, k2, k3, True], YukTil[b2, k3, k4]] 
				*Lam[b1, b2, $c, $d, True] // Sym[$a, $b, $c, $d];
			bl[2, 28] := 0 (*Tr @ Tdot[Yuk[$a, k4, k1], YukTil[$b, k1, k2], Yuk[b1, k2, k3], YukTil[b2, k3, k4]] 
				*Lam[b1, b2, $c, $d] // Sym[$a, $b, $c, $d]*);
			bl[2, 29] := Y4cS[$a, b1] Lam[b1, $b, $c, $d, True] // Sym4[$a, $b, $c, $d];
			bl[2, 30] := Y2SY2F[$a, b1] Lam[b1, $b, $c, $d, True] // Sym4[$a, $b, $c, $d];
		(*y^6 terms*)
			bl[2, 31] := Tr @ Tdot[Yuk[b1, k6, k1], YukTil[$a, k1, k2, True], Yuk[b1, k2, k3], YukTil[$b, k3, k4, True], 
				Yuk[$c, k4, k5, True], YukTil[$d, k5, k6, True]] // Sym[$a, $b, $c, $d];
			bl[2, 32] := Tr @ Tdot[Yuk[$a, k6, k1, True], YukTil[b1, k1, k2], Yuk[$b, k2, k3, True], YukTil[$c, k3, k4, True], 
				Yuk[b1, k4, k5], YukTil[$d, k5, k6, True]] // Sym[$b, $c, $d];
			bl[2, 33] := Tr @ Tdot[Yuk[$a, k5, k1, True], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True], 
				Y2F[k4, k5]] // Sym[$a, $b, $c, $d];
		];
		
		With[{diagrams = {5, 33}[[loop]]},
			Monitor[
				Return[Sum[Bcoef[3, loop, n] bl[loop, n], {n, diagrams}] /. $quarticCoefficients]
			,StringForm["Evaluating term `` / ``", n, diagrams]]
		]
	];


(*#####################################*)
(*----------Beta Coefficients----------*)
(*#####################################*)
$quarticCoefficients = {
	(* 1-loop *)
	Bcoef[3, 1, 1] -> 36,
	Bcoef[3, 1, 2] -> -12,
	Bcoef[3, 1, 3] -> 3,
	Bcoef[3, 1, 4] -> 2,
	Bcoef[3, 1, 5] -> -12,
	(*2-loop*)
	Bcoef[3, 2, 1] -> 324,
	Bcoef[3, 2, 2] -> -684,
	Bcoef[3, 2, 3] -> 646,
	Bcoef[3, 2, 4] -> -28,
	Bcoef[3, 2, 5] -> -32,
	Bcoef[3, 2, 6] -> 12,
	Bcoef[3, 2, 7] -> 60,
	Bcoef[3, 2, 8] -> 0,
	Bcoef[3, 2, 9] -> 6,
	Bcoef[3, 2, 10] -> -143/3,
	Bcoef[3, 2, 11] -> 11/3,
	Bcoef[3, 2, 12] -> 10/3,
	Bcoef[3, 2, 13] -> -18,
	Bcoef[3, 2, 14] -> 24,
	Bcoef[3, 2, 15] -> -18,
	Bcoef[3, 2, 16] -> 1/3,
	Bcoef[3, 2, 17] -> -6,
	Bcoef[3, 2, 18] -> 0,
	Bcoef[3, 2, 19] -> -144,
	Bcoef[3, 2, 20] -> 60,
	Bcoef[3, 2, 21] -> 10,
	Bcoef[3, 2, 22] -> 0,
	Bcoef[3, 2, 23] -> -3,
	Bcoef[3, 2, 24] -> 0,
	Bcoef[3, 2, 25] -> 24,
	Bcoef[3, 2, 26] -> -48,
	Bcoef[3, 2, 27] -> 12,
	Bcoef[3, 2, 28] -> 0,
	Bcoef[3, 2, 29] -> -2,
	Bcoef[3, 2, 30] -> -3,
	Bcoef[3, 2, 31] -> 48,
	Bcoef[3, 2, 32] -> 24,
	Bcoef[3, 2, 33] -> 24
}; 

$yukawaCoefficients = {
	(* 1-loop *)
	Bcoef[2, 1, 1] -> 0,
	Bcoef[2, 1, 2] -> -6,
	Bcoef[2, 1, 3] -> 2,
	Bcoef[2, 1, 4] -> 1,
	Bcoef[2, 1, 5] -> 1/2,
	(* 2-loop *)
	Bcoef[2, 2, 1] -> -21/2,
	Bcoef[2, 2, 2] -> 12, 
	Bcoef[2, 2, 3] -> 0,
	Bcoef[2, 2, 4] -> -3,
	Bcoef[2, 2, 5] -> 49/4,
	Bcoef[2, 2, 6] -> -1/4,
	Bcoef[2, 2, 7] -> -1/2,
	Bcoef[2, 2, 8] -> -97/3,
	Bcoef[2, 2, 9] -> 11/6,
	Bcoef[2, 2, 10] -> 5/3,
	Bcoef[2, 2, 11] -> 1/12,
	Bcoef[2, 2, 12] -> 12,
	Bcoef[2, 2, 13] -> 0,
	Bcoef[2, 2, 14] -> 6,
	Bcoef[2, 2, 15] -> -12,
	Bcoef[2, 2, 16] -> 10,
	Bcoef[2, 2, 17] -> 6,
	Bcoef[2, 2, 18] -> 5/2,
	Bcoef[2, 2, 19] -> 9,
	Bcoef[2, 2, 20] -> -1/2,
	Bcoef[2, 2, 21] -> -7/2,
	Bcoef[2, 2, 22] -> 0,
	Bcoef[2, 2, 23] -> -2,
	Bcoef[2, 2, 24] -> 2,
	Bcoef[2, 2, 25] -> 0,
	Bcoef[2, 2, 26] -> -2,
	Bcoef[2, 2, 27] -> 0,
	Bcoef[2, 2, 28] -> -1/2,
	Bcoef[2, 2, 29] -> -2,
	Bcoef[2, 2, 30] -> -1/4,
	Bcoef[2, 2, 31] -> -3/4,
	Bcoef[2, 2, 32] -> -1,
	Bcoef[2, 2, 33] -> -3/4
};

$gaugeCoefficients = {
	(* 1-loop *)
	Bcoef[1, 1, 1] -> -22/3,
	Bcoef[1, 1, 2] -> 2/3,
	Bcoef[1, 1, 3] -> 1/3,
	(* 2-loop *)
	Bcoef[1, 2, 1] -> 2,
	Bcoef[1, 2, 2] -> 4,
	Bcoef[1, 2, 3] -> -68/3,
	Bcoef[1, 2, 4] -> 10/3,
	Bcoef[1, 2, 5] -> 2/3,
	Bcoef[1, 2, 6] -> -1,
	Bcoef[1, 2, 7] -> 0,
	(* 3-loop *)
	Bcoef[1, 3, 1] -> -1,
	Bcoef[1, 3, 2] -> 29/2,
	Bcoef[1, 3, 3] -> 133/18,
	Bcoef[1, 3, 4] -> 679/36,
	Bcoef[1, 3, 5] -> -11/18,
	Bcoef[1, 3, 6] -> -25/18,
	Bcoef[1, 3, 7] -> -23/36,
	Bcoef[1, 3, 8] -> -49/36,
	Bcoef[1, 3, 9] -> 4,
	Bcoef[1, 3, 10] -> 25/2,
	Bcoef[1, 3, 11] -> -2857/27,
	Bcoef[1, 3, 12] -> -79/108,
	Bcoef[1, 3, 13] -> 1/54,
	Bcoef[1, 3, 14] -> 1415/54,
	Bcoef[1, 3, 15] -> 545/108,
	Bcoef[1, 3, 16] -> -29/54,
	Bcoef[1, 3, 17] -> 1,
	Bcoef[1, 3, 18] -> -1/12,
	Bcoef[1, 3, 19] -> -5/4,
	Bcoef[1, 3, 20] -> -1/4,
	Bcoef[1, 3, 21] -> -1,
	Bcoef[1, 3, 22] -> -7,
	Bcoef[1, 3, 23] -> -7/2,
	Bcoef[1, 3, 24] -> -6,
	Bcoef[1, 3, 25] -> 9/4,
	Bcoef[1, 3, 26] -> 1,
	Bcoef[1, 3, 27] -> -1,
	Bcoef[1, 3, 28] -> 3/2,
	Bcoef[1, 3, 29] -> 7/8,
	Bcoef[1, 3, 30] -> 1/2,
	Bcoef[1, 3, 31] -> 1/8,
	Bcoef[1, 3, 32] -> 3/8,
	Bcoef[1, 3, 33] -> -1/8
};

Protect[$gaugeCoefficients, $quarticCoefficients, $yukawaCoefficients];



End[]





	