(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Begin["TensorCalculations`"]
<< BetaTensors`; (*Loads all the tensors from ancillary file*)

sig1 = {{0, 1}, {1, 0}};
sig3 = {{1, 0}, {0, -1}};
sig3Til = {{-1, 0}, {0, 1}};

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
];

(* Extra identities required to handle the more complicated combinations of struture constants appearing in the 4-loop gauge beta functions. *)
NewIdentities[expr_] := Block[{fStruct, tGen, n},
	fStruct /: fStruct[gr_, x : OrderlessPatternSequence[A_, B_, C_]] tGen[gr_[rep_], B_, a_, b_] tGen[gr_[rep_], C_, b_, c_] :=
    	I/2 Signature[List@x] Signature[{A, B, C}] Casimir2[gr@adj] tGen[gr@rep, A, a, c];
	fStruct /: fStruct[gr_, OrderlessPatternSequence[A_, C_, E_]] fStruct[gr_, OrderlessPatternSequence[B_, D_, E_]] *
		tGen[gr_@rep_, A_, a_,b_] tGen[gr_@rep_, B_, b_, c_] tGen[gr_@rep_, C_, c_, d_] tGen[gr_@rep_, D_, d_, a_] := 0;
	fStruct /: fStruct[gr_, x1 : OrderlessPatternSequence[E_, A_, F_]] fStruct[gr_, x2 : OrderlessPatternSequence[F_, B_, G_]] *
		fStruct[gr_, x3 : OrderlessPatternSequence[G_, C_, H_]] fStruct[gr_, x4 : OrderlessPatternSequence[H_, D_, E_]] fStruct[gr_, y1 : OrderlessPatternSequence[I_, A_, J_]] *
		fStruct[gr_, y2 : OrderlessPatternSequence[J_, B_, K_]] fStruct[gr_, y3 : OrderlessPatternSequence[K_, C_, L_]] fStruct[gr_, y4 : OrderlessPatternSequence[L_, D_, I_]] /;
		Head@$gaugeGroups[gr, LieGroup] == SU :=
		(n = $gaugeGroups[gr, LieGroup][[1]]; Signature[List@x1] Signature[List@x2] Signature[List@x3] Signature[List@x4] Signature[List@y1]  *
			Signature[List@y2] Signature[List@y3] Signature[List@y4] Signature[{E, A, F}] Signature[{F, B, G}] Signature[{G, C, H}] Signature[{H, D, E}] *
			Signature[{I, A, J}] Signature[{J, B, K}] Signature[{K, C, L}] Signature[{L, D, I}] n^2 (n^2 - 1) (n^2 + 12)/8);
	expr
];

(*##################################*)
(*------------RG tensors------------*)
(*##################################*)
(* These functions project out specific coupling beta functions from the general tensors.*)

GaugeTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {3, 7, 33, 202}[[loop]], n, C1, C2, proj, bTerm},
		(* The 0-loop contribution in d = 4 - \epsilon dimensions *)
		If[loop === 0,
			bTerm = Expand[ - Global`\[Epsilon] $gaugeGroups[$couplings @ coupling, Projector][C1, C2] G2Matrix[C1, C2] ];
			If[ MatchQ[ $gaugeGroups[$couplings @ coupling, LieGroup], U1[m_] /; m > 1],
				bTerm = bTerm /. coupling -> $gaugeGroups[$couplings @ coupling, CouplingMatrix];
			];
			Return @ bTerm;
		];

		proj = Ttimes[G2Matrix[$A, C1], $gaugeGroups[$couplings @ coupling, Projector][C1, C2], G2Matrix[C2, $B]];
		Monitor[
			bTerm = Sum[
				Bcoef[1, loop, n] Ttimes[proj, BetaTensor[1, loop, n] ],
			{n, diagrams}];
		,StringForm["Evaluating term `` / ``", n, diagrams]];

		(* At 4-loops additional identities are required.*)
		If[loop > 3,
			NewIdentities @ bTerm
		,
			bTerm
		]
	];

YukawaTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33, 308}[[loop + 1]], n, proj},
		proj = $yukawas[coupling, Projector][$a, $i, $j];
		(* The BetaTensors hav been explicitly symmetrized up to 2-loop order to save time. At 3-loop symmetrization is done here. *)
		If[loop > 2,
			proj = TsSym[$i, $j] @ proj;
		];
		Monitor[
			(* The projection point depends on the "chirlality ascribed to the given Yukawa coupling/" *)
			Switch[$yukawas[coupling, Chirality]
			,Left,
				Sum[
					Bcoef[2, loop, n] Ttimes[proj, BetaTensor[2, loop, n][[1, 1]] ],
				{n, diagrams}]
			,Right,
				Sum[
					Bcoef[2, loop, n] Ttimes[proj, BetaTensor[2, loop, n][[2, 2]] ],
				{n, diagrams}]
			]
		,StringForm["Evaluating term `` / ``", n, diagrams] ]
	];

FermionMassTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33}[[loop + 1]], n},
		Monitor[
			(* The projection point depends on the "chirlality ascribed to the given Yukawa coupling/" *)
			Switch[$fermionMasses[coupling, Chirality]
			,Left,
				Sum[
					Bcoef[2, loop, n] Ttimes[$fermionMasses[coupling, Projector][$a, $i, $j],
					BetaTensor[4, loop, n][[1, 1]] ],
				{n, diagrams}]
			,Right,
				Sum[
					Bcoef[2, loop, n] Ttimes[$fermionMasses[coupling, Projector][$a, $i, $j],
					BetaTensor[4, loop, n][[2, 2]] ],
				{n, diagrams}]
			]
		,StringForm["Evaluating term `` / ``", n, diagrams] ]
	];

QuarticTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33}[[loop + 1]], n},
		Monitor[
			Sum[
				Bcoef[3, loop, n] Ttimes[$quartics[coupling, Projector][$a, $b, $c, $d],
				BetaTensor[3, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

ScalarMassiveTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33}[[loop + 1]], n},
		Monitor[
			Switch[$couplings @ coupling
			,Trilinear,
				Sum[
					Bcoef[3, loop, n] Ttimes[$trilinears[coupling, Projector][$a, $b, $c, $d],
					BetaTensor[5, loop, n] ],
				{n, diagrams}]
			,ScalarMass,
				Sum[
					Bcoef[3, loop, n] Ttimes[$scalarMasses[coupling, Projector][$a, $b, $c, $d],
					BetaTensor[5, loop, n] ],
				{n, diagrams}]
			]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

FermionAnomalousTensors[field_, loop_Integer] :=
	Module[{diagrams = {2, 9}[[loop]], n},
		Monitor[
			Sum[
				Acoef[1, loop, n] Ttimes[$fermions[field, Projector][$i, $j], (AnomalousTensor[1, loop, n] = AnomalousTensor[1, loop, n])[[1,1]] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

ScalarAnomalousTensors[field_, loop_Integer] :=
	Module[{diagrams = {2, 8}[[loop]], n},
		Monitor[
			Sum[
				Acoef[2, loop, n] Ttimes[$scalars[field, Projector][$a, $b], AnomalousTensor[2, loop, n] = AnomalousTensor[2, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];



(*#####################################*)
(*----------Beta Coefficients----------*)
(*#####################################*)
$quarticCoefficients = {
	Bcoef[3, 0, 1] -> 1,
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
	Bcoef[2, 0, 1] -> 1,
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
	Bcoef[1, 0, 1] -> 1,
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

$fermionAnomalousCoefficients = {
	(* 1-loop *)
	Acoef[1, 1, 1] -> Global`\[Xi],
	Acoef[1, 1, 2] -> 1/2,
	(* 2-loop *)
	Acoef[1, 2, 1] -> -3/2,
	Acoef[1, 2, 2] -> 25/4 + 2 Global`\[Xi] + Global`\[Xi]^2/4,
	Acoef[1, 2, 3] -> -1/4,
	Acoef[1, 2, 4] -> -1/2,
	Acoef[1, 2, 5] -> 9/2,
	Acoef[1, 2, 6] -> -1/4,
	Acoef[1, 2, 7] -> -7/4,
	Acoef[1, 2, 8] -> -1/8,
	Acoef[1, 2, 9] -> -3/8
};

$scalarAnomalousCoefficients = {
	(* 1-loop *)
	Acoef[2, 1, 1] -> Global`\[Xi] - 3,
	Acoef[2, 1, 2] -> 1/2,
	(* 2-loop *)
	Acoef[2, 2, 1] -> 3/2,
	Acoef[2, 2, 2] -> Global`\[Xi]^2/4 + 2 Global`\[Xi] - 35/3,
	Acoef[2, 2, 3] -> 11/12,
	Acoef[2, 2, 4] -> 5/6,
	Acoef[2, 2, 5] -> 1/12,
	Acoef[2, 2, 6] -> 5/2,
	Acoef[2, 2, 7] -> -1/2,
	Acoef[2, 2, 8] -> -3/4
};


Protect[$gaugeCoefficients, $quarticCoefficients, $yukawaCoefficients, $fermionAnomalousCoefficients, $scalarAnomalousCoefficients];

End[]
