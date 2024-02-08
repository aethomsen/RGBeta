(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageExport["Acoef"]
PackageExport["Bcoef"]
PackageExport["Ucoef"]

PackageScope["FermionAnomalousTensors"]
PackageScope["FermionMassTensors"]
PackageScope["GaugeTensors"]
PackageScope["QuarticTensors"]
PackageScope["ResetBetas"]
PackageScope["ScalarAnomalousTensors"]
PackageScope["ScalarMassiveTensors"]
PackageScope["YukawaTensors"]

PackageScope["UpsilonYukawaTensors"]
PackageScope["FermionUpsilonTensors"]
PackageScope["ScalarUpsilonTensors"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

Acoef::usage =
	"Acoef[i, j, k] is used internally to denote the coefficients of the differnet tensor contraction appering in the anomalous dimensions."
Bcoef::usage =
	"Bcoef[i, j, k] is used internally to denote the coefficients of the differnet tensor contraction appering in the beta functions."
Ucoef::usage =
	"Ucoef[i, j, k] is used internally to denote the coefficients of the differnet tensor contraction appering in the upsilon functions."

FermionAnomalousTensors::uasge =
	"FermionAnomalousTensors[field, loop] evaluates all tensor contractions of the general anomalous dimension tensor and the field projector."
FermionMassTensors::usage =
	"FermionMassTensors[coupling, loop] is a function that computes all the tensor contractions used the general fermion mass beta function at the given loop order."
GaugeTensors::usage =
	"GaugeTensors[coupling, loop] is a function that computes all the tensor contractions used the general gauge beta function at the given loop order."
QuarticTensors::usage =
	"QuarticTensors[coupling, loop] is a function that computes all the tensor contractions used the general quartic beta function at the given loop order."
ResetBetas::usage =
	"ResetBetas[] is a function used to dump all internally stored beta computations from the kernel."
ScalarAnomalousTensors::uasge =
	"ScalarAnomalousTensors[field, loop] evaluates all tensor contractions of the general anomalous dimension tensor and the field projector."
ScalarMassiveTensors::usage =
	"ScalarMassiveTensors[coupling, loop] is a function that computes all the tensor contractions used in the the trillinear and scalar mass beta functions at the given loop order."
YukawaTensors::usage =
	"YukawaTensors[coupling, loop] is a function that computes all the tensor contractions used the general yukawa beta function at the given loop order."

(*########################################*)
(*----------All kinds of tensors----------*)
(*########################################*)

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
	Module[{diagrams = {3, 7, 33, 202}[[loop]], n= 0, C1, C2, proj, bTerm},
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
				Bcoef[1, loop, n] RefineGroupStructures @ Ttimes[proj, BetaTensor[1, loop, n] ],
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
	Module[{diagrams = {1, 5, 33, 308}[[loop + 1]], n= 0},
		Monitor[
			Sum[
				Bcoef[2, loop, n] RefineGroupStructures @ Tr@ Tdot[$yukawas[coupling, Projector][$a, $i, $j], BetaTensor[2, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams] ]
	];

FermionMassTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33}[[loop + 1]], n= 0},
		Monitor[
			Sum[
				Bcoef[2, loop, n] RefineGroupStructures @ Tr@ Tdot[$fermionMasses[coupling, Projector][$a, $i, $j], BetaTensor[4, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams] ]
	];

QuarticTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33, 315}[[loop + 1]], n= 0},
		Monitor[
			Sum[
				Bcoef[3, loop, n] RefineGroupStructures @ Ttimes[$quartics[coupling, Projector][$a, $b, $c, $d],
				BetaTensor[3, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

ScalarMassiveTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {1, 5, 33}[[loop + 1]], n= 0},
		Monitor[
			Switch[$couplings @ coupling
			,Trilinear,
				Sum[
					Bcoef[3, loop, n] RefineGroupStructures @ Ttimes[$trilinears[coupling, Projector][$a, $b, $c, $d],
					BetaTensor[5, loop, n] ],
				{n, diagrams}]
			,ScalarMass,
				Sum[
					Bcoef[3, loop, n] RefineGroupStructures @ Ttimes[$scalarMasses[coupling, Projector][$a, $b, $c, $d],
					BetaTensor[5, loop, n] ],
				{n, diagrams}]
			]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

(* For the projection of scalar and fermion fields *)
ScalarFieldProjector[{phi1_, phi2_}, $da_, $db_] :=
	Module[{dim},
		dim = Length@ $fieldIndexMap["Scalars"];
		TStructure[$scalar@ $da, $scalar@ $db] @ SparseArray[
			$fieldIndexMap["Scalars"] /@ {Bar@ phi1, phi2} ->
				Product[del[rep, $da, $db]/ Dim@ rep, {rep, $scalars[phi1, GaugeRep]} ] /
					Product[If[$scalars[scal, SelfConjugate], 1, Sqrt@ 2], {scal, {phi1, phi2}}], (*Normalization for complex scalars*)
			dim]
	]

FermionFieldProjector[{psi1_, psi2_}, $di_, $dj_] :=
	Module[{dim},
		dim = Length@ $fieldIndexMap["Fermions"];
		DiagonalMatrix@ {
			TStructure[$fermion@ $di, $fermion@ $dj] @ SparseArray[
				$fieldIndexMap["Fermions"] /@ {psi1, psi2} -> Product[del[rep, $di, $dj]/ Dim@ rep, {rep, $fermions[psi1, GaugeRep]} ],
				dim], 0}
	]

(* Anomalous dimensions *)
FermionAnomalousTensors[fields_, loop_Integer] :=
	Module[{diagrams = {2, 9}[[loop]], n= 0},
		Monitor[
			Sum[
				Acoef[1, loop, n] RefineGroupStructures @ Tr@ Tdot[FermionFieldProjector[fields, $i, $j], AnomalousTensor[1, loop, n] = AnomalousTensor[1, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

ScalarAnomalousTensors[fields_, loop_Integer] :=
	Module[{diagrams = {2, 8}[[loop]], n= 0},
		Monitor[
			Sum[
				Acoef[2, loop, n] RefineGroupStructures @ Ttimes[ScalarFieldProjector[fields, $a, $b], AnomalousTensor[2, loop, n] = AnomalousTensor[2, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

(* For the evlauation of upsilon tensors*)
MinusHc[i1_, i2_][expr_]:= expr - (sig1. expr. sig1/. TStructure[inds__][sa_] :>
	SymmetrizeTS[TStructure[inds]@ sa, {i1, i2}, {{2, 1}}]);
MinusTrans[i1_, i2_][expr_]:= expr - (expr/. TStructure[inds__][sa_] :>
	SymmetrizeTS[TStructure[inds]@ sa, {i1, i2}, {{2, 1}}]);

CompleteReplace[expr_, repRules_]:= expr/. repRules/. x_SparseArray:>
	SparseArray[ArrayRules@ x/. repRules, Dimensions@ x];

(* Returns (\upsilon g)_aij for use in B = \beta - (\upsilon g)*)
UpsilonYukawaTensors[coupling_Symbol, loop_Integer] :=
	Module[{diagrams = {{0, 0}, {0, 0}, {6, 3}}[[loop]], n= 0, d= 0, kind},
		Monitor[
			Sum[
				d++;
				Ucoef[1, loop, n]( Tdot[$yukawas[coupling, Projector][$a, $i, $j], CompleteReplace[UpsilonTensor[1, loop, n], $j-> kind], Yuk[$a, kind, $j]]-
					Tdot[$yukawas[coupling, Projector][$a, $i, $j], Yuk[$a, $i, kind], sig1. CompleteReplace[UpsilonTensor[1, loop, n], $i-> kind]. sig1])// Tr// RefineGroupStructures
			, {n, diagrams[[1]]}] +
			Sum[
				d++;
				Ucoef[2, loop, n] RefineGroupStructures @ Tr@ Tdot[Ttimes[$yukawas[coupling, Projector][$a, $i, $j], UpsilonTensor[2, loop, n]], Yuk[$b, $i, $j]]
			, {n, diagrams[[2]]}]
		,StringForm["Evaluating term `` / ``", d, Plus@@ diagrams]]
	];

FermionUpsilonTensors[fields_, loop_Integer] :=
	Module[{diagrams = {0, 0, 6}[[loop]], n= 0},
		Monitor[
			Sum[
				Ucoef[1, loop, n] Tr@ Tdot[FermionFieldProjector[fields, $i, $j], UpsilonTensor[1, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];

ScalarUpsilonTensors[fields_, loop_Integer] :=
	Module[{diagrams = {0, 0, 3}[[loop]], n= 0},
		Monitor[
			Sum[
				Ucoef[2, loop, n] Ttimes[ScalarFieldProjector[fields, $a, $b], UpsilonTensor[2, loop, n] ],
			{n, diagrams}]
		,StringForm["Evaluating term `` / ``", n, diagrams]]
	];




(*###########################################*)
(*-----------All tensor structures-----------*)
(*###########################################*)



ResetBetas[] := Module[{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,
    k1,k2,k3,k4,k5,k6,k7,k8,b1,b2,b3,b4,b5,b6,b7,b8},
(*###################################*)
(*-----------Gauge Tensors-----------*)
(*###################################*)
    BetaTensor[1, 1, 1] := BetaTensor[1, 1, 1] = C2G[$A, $B];
    BetaTensor[1, 1, 2] := BetaTensor[1, 1, 2] = S2F[$A, $B];
    BetaTensor[1, 1, 3] := BetaTensor[1, 1, 3] = S2S[$A, $B];

    BetaTensor[1, 2, 1] := BetaTensor[1, 2, 1] = S2FC2F[$A, $B];
    BetaTensor[1, 2, 2] := BetaTensor[1, 2, 2] = S2SC2S[$A, $B];
    BetaTensor[1, 2, 3] := BetaTensor[1, 2, 3] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], C2G[A2, $B]];
    BetaTensor[1, 2, 4] := BetaTensor[1, 2, 4] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], S2F[A2, $B]];
    BetaTensor[1, 2, 5] := BetaTensor[1, 2, 5] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], S2S[A2, $B]];
    BetaTensor[1, 2, 6] := BetaTensor[1, 2, 6] = S2FY2F[$A, $B];
    BetaTensor[1, 2, 7] := BetaTensor[1, 2, 7] = 0 (*S2SY2S[$A, $B]*);
(*y^0 terms*)
    BetaTensor[1, 3, 1] := BetaTensor[1, 3, 1] = Tr@Tdot[Tferm[$B, k4, k1], Tferm[$A, k1, k2], C2F[k2, k3], C2F[k3, k4]];
    BetaTensor[1, 3, 2] := BetaTensor[1, 3, 2] = Ttimes[Tscal[$B, b4, b1], Tscal[$A, b1, b2], C2S[b2, b3], C2S[b3, b4]];
    BetaTensor[1, 3, 3] := BetaTensor[1, 3, 3] = Tr@Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], C2FC2G[k2, k3]];
    BetaTensor[1, 3, 4] := BetaTensor[1, 3, 4] = Ttimes[Tscal[$B, b3, b1], Tscal[$A, b1, b2], C2SC2G[b2, b3]];
    BetaTensor[1, 3, 5] := BetaTensor[1, 3, 5] = Tr@Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], C2FS2F[k2, k3]];
    BetaTensor[1, 3, 6] := BetaTensor[1, 3, 6] = Ttimes[Tscal[$B, b3, b1], Tscal[$A, b1, b2], C2SS2F[b2, b3]];
    BetaTensor[1, 3, 7] := BetaTensor[1, 3, 7] = Tr@Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], C2FS2S[k2, k3]];
    BetaTensor[1, 3, 8] := BetaTensor[1, 3, 8] = Ttimes[Tscal[$B, b3, b1], Tscal[$A, b1, b2], C2SS2S[b2, b3]];
    BetaTensor[1, 3, 9] := BetaTensor[1, 3, 9] = Ttimes[S2FC2F[$A, A1], G2Matrix[A1, A2], C2G[A2, $B]];
    BetaTensor[1, 3, 10] := BetaTensor[1, 3, 10] = Ttimes[S2SC2S[$A, A1], G2Matrix[A1, A2], C2G[A2, $B]];
    BetaTensor[1, 3, 11] := BetaTensor[1, 3, 11] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], C2G[A2, A3], G2Matrix[A3, A4], C2G[A4, $B]];
    BetaTensor[1, 3, 12] := BetaTensor[1, 3, 12] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], S2F[A2, A3], G2Matrix[A3, A4], S2F[A4, $B]];
    BetaTensor[1, 3, 13] := BetaTensor[1, 3, 13] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], S2S[A2, A3], G2Matrix[A3, A4], S2S[A4, $B]];
    BetaTensor[1, 3, 14] := BetaTensor[1, 3, 14] = Ttimes[C2G[$A, A1], G2Matrix[A1, A2], C2G[A2, A3], G2Matrix[A3, A4], S2F[A4, $B]];
    BetaTensor[1, 3, 15] := BetaTensor[1, 3, 15] =  Ttimes[C2G[$A, A1], G2Matrix[A1, A2], C2G[A2, A3], G2Matrix[A3, A4], S2S[A4, $B]];
    BetaTensor[1, 3, 16] := BetaTensor[1, 3, 16] = Ttimes[S2F[$A, A1], G2Matrix[A1, A2], C2G[A2, A3], G2Matrix[A3, A4], S2S[A4, $B]];
    BetaTensor[1, 3, 17] := BetaTensor[1, 3, 17] = Ttimes[Tscal[$A, b1, b2], TsG2[A1, b2, b3], Lam[b1, b3, b4, b6], Tscal[$B, b4, b5], Tscal[A1, b5, b6]];
    BetaTensor[1, 3, 18] := BetaTensor[1, 3, 18] = Ttimes[Tscal[$A, b1, b2], Lam2[b2, b3], Tscal[$B, b3, b1]];
(*y^2 terms*)
    BetaTensor[1, 3, 19] := BetaTensor[1, 3, 19] = Tr @ Tdot[Tferm[$B, k4, k1], Tferm[$A, k1, k2], Y2Ft[k2, k3], C2F[k3, k4]];
    BetaTensor[1, 3, 20] := BetaTensor[1, 3, 20] = Tr @ Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], Y2FC2Ft[k2, k3]];
    BetaTensor[1, 3, 21] := BetaTensor[1, 3, 21] = Ttimes[Tscal[$B, b3, b1], Tscal[$A, b1, b2], Y2SC2F[b2, b3]];
    BetaTensor[1, 3, 22] := BetaTensor[1, 3, 22] = Tr @ Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], Y2FC2St[k2, k3]];
    BetaTensor[1, 3, 23] := BetaTensor[1, 3, 23] = Ttimes[Tscal[$B, b4, b1], Tscal[$A, b1, b2], Y2S[b2, b3], C2S[b3, b4]];
    BetaTensor[1, 3, 24] := BetaTensor[1, 3, 24] = Ttimes[S2FY2F[$A, A1], G2Matrix[A1, A2], C2G[A2, $B]];
    BetaTensor[1, 3, 25] := BetaTensor[1, 3, 25] = Ttimes[S2SY2S[$A, A1], G2Matrix[A1, A2], C2G[A2, $B]];
(*y^4 terms*)
    BetaTensor[1, 3, 26] := BetaTensor[1, 3, 26] = Tr @ Tdot[Yuk[b1, b6, k1], Tferm[$A, k1, k2], YukTil[b1, k2, k3] Yuk[b2, k3, k4], Tferm[$B, k4, b5], YukTil[b2, b5, b6]];
    BetaTensor[1, 3, 27] := BetaTensor[1, 3, 27] = Ttimes[Tscal[$B, b3, b1], Tscal[$A, b1, b2], Y4cS[b2, b3]];
    BetaTensor[1, 3, 28] := BetaTensor[1, 3, 28] = Tr @ Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], Y4cFt[k2, k3]];
    BetaTensor[1, 3, 29] := BetaTensor[1, 3, 29] = Tr @ Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], Y2FY2St[k2, k3]];
    BetaTensor[1, 3, 30] := BetaTensor[1, 3, 30] = Ttimes[Tscal[$B, b3, b1], Tscal[$A, b1, b2], Y2SY2F[b2, b3]];
    BetaTensor[1, 3, 31] := BetaTensor[1, 3, 31] = Tr @ Tdot[Tferm[$B, k4, k1], Tferm[$A, k1, k2], Y2Ft[k2, k3], Y2Ft[k3, k4]];
    BetaTensor[1, 3, 32] := BetaTensor[1, 3, 32] = Tr @ Tdot[Tferm[$B, k3, k1], Tferm[$A, k1, k2], Y2FY2Ft[k2, k3]];
    BetaTensor[1, 3, 33] := BetaTensor[1, 3, 33] = Ttimes[Tscal[$B, b4, b1], Tscal[$A, b1, b2], Y2S[b2, b3], Y2S[b3, b4]];

 (* 4-loop *)
    BetaTensor[1, 4, 1] := BetaTensor[1, 4, 1] = Ttimes[FGauge[A1, A2, A3], FGauge[A4, A5, A6], FGauge[A7, A8, A9], FGauge[A10, A11, A12], FGauge[A13, A14, A15], FGauge[A16, A17, $B], FGauge[A18, A19, A20], FGauge[$A, A21, A22], G2Matrix[A1, A4], G2Matrix[A2, A7], G2Matrix[A3, A10], G2Matrix[A5, A19], G2Matrix[A6, A21], G2Matrix[A8, A14], G2Matrix[A9, A22], G2Matrix[A11, A15], G2Matrix[A12, A20], G2Matrix[A13, A16], G2Matrix[A17, A18]];
    BetaTensor[1, 4, 2] := BetaTensor[1, 4, 2] = Ttimes[Tscal[A1, b1, b2], Tscal[A2, b3, b4], Tscal[A3, b2, b3], Tscal[$A, b4, b1], Tscal[$B, b5, b6], TsG2[A1, b7, b5], TsG2[A2, b6, b8], TsG2[A3, b8, b7]];
    BetaTensor[1, 4, 3] := BetaTensor[1, 4, 3] = Ttimes[FGauge[A1, A2, A3], FGauge[$A, A4, A5], G2Matrix[A1, A4], Tscal[A6, b1, b2], Tscal[$B, b3, b1], TsG2[A2, b4, b3], TsG2[A3, b5, b6], TsG2[A5, b2, b5], TsG2[A6, b6, b4]];
    BetaTensor[1, 4, 4] := BetaTensor[1, 4, 4] = Ttimes[Tr[Tdot[TfermTil[A1, k1, k2], TfermTil[A2, k2, k3], TfermTil[A3, k3, k4], TfermTil[$A, k4, k1]]], Tscal[$B, b1, b2], TsG2[A1, b2, b3], TsG2[A2, b3, b4], TsG2[A3, b4, b1]];
    BetaTensor[1, 4, 5] := BetaTensor[1, 4, 5] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], Tr[Tdot[TfermTil[A1, k1, k2], TfermTil[A3, k2, k3], TfermTil[A5, k3, k4], TfermTil[$B, k4, k1]]], Tr[Tdot[TfermTil[A4, k5, k6], TfermTil[A2, k6, k7], TfermTil[$A, k7, k8], TfermTil[A6, k8, k5]]]];
    BetaTensor[1, 4, 6] := BetaTensor[1, 4, 6] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], Tr[Tdot[TfermTil[A1, k1, k2], TfermTil[A3, k2, k3], TfermTil[$A, k3, k4], TfermTil[A2, k4, k5], TfermTil[A5, k5, k6], TfermTil[$B, k6, k7], TfermTil[A4, k7, k8], TfermTil[A6, k8, k1]]]];
    BetaTensor[1, 4, 7] := BetaTensor[1, 4, 7] = Ttimes[S2SC2S[A1, A2], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A1, b3, b4], TsG2[A2, b4, b1]];
    BetaTensor[1, 4, 8] := BetaTensor[1, 4, 8] = Ttimes[S2FC2F[A1, A2], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A1, b4, b1], TsG2[A2, b3, b4]];
    BetaTensor[1, 4, 9] := BetaTensor[1, 4, 9] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], S2SC2S[A1, A3], Tr[Tdot[TfermTil[A4, k1, k2], TfermTil[A2, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 10] := BetaTensor[1, 4, 10] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], S2FC2F[A1, A3], Tr[Tdot[TfermTil[$B, k1, k2], TfermTil[A4, k2, k3], TfermTil[A2, k3, k4], TfermTil[$A, k4, k1]]]];
    BetaTensor[1, 4, 11] := BetaTensor[1, 4, 11] = Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], C2Ft[k3, k4], C2Ft[k4, k5], TfermTil[$A, k5, k1]]];
    BetaTensor[1, 4, 12] := BetaTensor[1, 4, 12] = Ttimes[C2S[b1, b2], C2S[b2, b3], C2S[b3, b4], Tscal[$A, b4, b5], Tscal[$B, b5, b1]];
    BetaTensor[1, 4, 13] := BetaTensor[1, 4, 13] = Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], C2FC2Gt[k3, k4], TfermTil[$A, k4, k1]]];
    BetaTensor[1, 4, 14] := BetaTensor[1, 4, 14] = Ttimes[C2S[b1, b2], C2SC2G[b2, b3], Tscal[$A, b3, b4], Tscal[$B, b4, b1]];
    BetaTensor[1, 4, 15] := BetaTensor[1, 4, 15] = Tr[Tdot[C2FS2Ft[k1, k2], C2Ft[k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]];
    BetaTensor[1, 4, 16] := BetaTensor[1, 4, 16] = Ttimes[C2S[b1, b2], C2SS2F[b2, b3], Tscal[$A, b3, b4], Tscal[$B, b4, b1]];
    BetaTensor[1, 4, 17] := BetaTensor[1, 4, 17] = Tr[Tdot[C2FS2St[k1, k2], C2Ft[k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]];
    BetaTensor[1, 4, 18] := BetaTensor[1, 4, 18] = Ttimes[C2S[b1, b2], C2SS2S[b3, b1], Tscal[$A, b4, b3], Tscal[$B, b2, b4]];
    BetaTensor[1, 4, 19] := BetaTensor[1, 4, 19] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[TfermTil[$A, k1, k2], C2Ft[k2, k3], C2Ft[k3, k4], TfermTil[A2, k4, k1]]]];
    BetaTensor[1, 4, 20] := BetaTensor[1, 4, 20] = Ttimes[C2G[A1, $A], C2S[b1, b2], C2S[b2, b3], Tscal[$B, b4, b1], TsG2[A1, b3, b4]];
    BetaTensor[1, 4, 21] := BetaTensor[1, 4, 21] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], S2F[A1, A3], S2F[A5, A4], Tr[Tdot[TfermTil[A2, k1, k2], TfermTil[A6, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 22] := BetaTensor[1, 4, 22] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], S2F[A1, A3], S2S[A5, A4], Tr[Tdot[TfermTil[A6, k1, k2], TfermTil[A2, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 23] := BetaTensor[1, 4, 23] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], S2S[A1, A3], S2S[A5, A4], Tr[Tdot[TfermTil[A2, k1, k2], TfermTil[A6, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 24] := BetaTensor[1, 4, 24] = Ttimes[G2Matrix[A1, A2], S2S[A3, A1], S2S[A4, A2], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A3, b4, b1], TsG2[A4, b3, b4]];
    BetaTensor[1, 4, 25] := BetaTensor[1, 4, 25] = Ttimes[C2G[A1, A2], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2F[A3, A5], Tr[Tdot[TfermTil[A6, k1, k2], TfermTil[A4, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 26] := BetaTensor[1, 4, 26] = Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2S[A3, A4], Tscal[$A, b1, b2], Tscal[$B, b3, b1], TsG2[A2, b4, b3], TsG2[A4, b2, b4]];
    BetaTensor[1, 4, 27] := BetaTensor[1, 4, 27] = Ttimes[G2Matrix[A1, A2], S2F[A3, A1], S2S[A2, A4], Tscal[$A, b1, b2], Tscal[$B, b3, b1], TsG2[A3, b4, b3], TsG2[A4, b2, b4]];
    BetaTensor[1, 4, 28] := BetaTensor[1, 4, 28] = Ttimes[C2G[A1, A2], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2S[A3, A5], Tr[Tdot[TfermTil[A6, k1, k2], TfermTil[A4, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 29] := BetaTensor[1, 4, 29] = Ttimes[C2G[A1, A2], C2G[A3, A4], G2Matrix[A1, A3], G2Matrix[A2, A5], G2Matrix[A4, A6], Tr[Tdot[TfermTil[$B, k1, k2], TfermTil[A6, k2, k3], TfermTil[A5, k3, k4], TfermTil[$A, k4, k1]]]];
    BetaTensor[1, 4, 30] := BetaTensor[1, 4, 30] = Ttimes[C2G[A1, A2], C2G[A3, A4], G2Matrix[A1, A3], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A2, b4, b1], TsG2[A4, b3, b4]];
    BetaTensor[1, 4, 31] := BetaTensor[1, 4, 31] = Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2F[A3, A4], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A2, b3, b4], TsG2[A4, b4, b1]];
    BetaTensor[1, 4, 32] := BetaTensor[1, 4, 32] = Ttimes[G2Matrix[A1, A2], S2F[A1, A3], S2F[A2, A4], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A3, b4, b1], TsG2[A4, b3, b4]];
    BetaTensor[1, 4, 33] := BetaTensor[1, 4, 33] = Ttimes[C2G[A1, A2], C2G[A3, $B], G2Matrix[A1, A3], G2Matrix[A2, A4], S2FC2F[A4, $A]];
    BetaTensor[1, 4, 34] := BetaTensor[1, 4, 34] = Ttimes[C2G[A1, A2], C2G[A3, $A], G2Matrix[A1, A3], G2Matrix[A2, A4], S2SC2S[$B, A4]];
    BetaTensor[1, 4, 35] := BetaTensor[1, 4, 35] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[TfermTil[$A, k1, k2], TfermTil[A2, k2, k3], C2FC2Gt[k3, k1]]]];
    BetaTensor[1, 4, 36] := BetaTensor[1, 4, 36] = Ttimes[C2G[A1, $A], C2SC2G[b1, b2], Tscal[$B, b2, b3], TsG2[A1, b3, b1]];
    BetaTensor[1, 4, 37] := BetaTensor[1, 4, 37] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[C2FS2Ft[k1, k2], TfermTil[A2, k2, k3], TfermTil[$A, k3, k1]]]];
    BetaTensor[1, 4, 38] := BetaTensor[1, 4, 38] = Ttimes[C2G[A1, $A], C2SS2F[b1, b2], Tscal[$B, b2, b3], TsG2[A1, b3, b1]];
    BetaTensor[1, 4, 39] := BetaTensor[1, 4, 39] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[C2FS2St[k1, k2], TfermTil[A2, k2, k3], TfermTil[$A, k3, k1]]]];
    BetaTensor[1, 4, 40] := BetaTensor[1, 4, 40] = Ttimes[C2G[A1, $A], C2SS2S[b1, b2], Tscal[$B, b3, b1], TsG2[A1, b2, b3]];
    BetaTensor[1, 4, 41] := BetaTensor[1, 4, 41] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], G2Matrix[A3, A4], S2F[A3, $A], S2FC2F[A2, A4]];
    BetaTensor[1, 4, 42] := BetaTensor[1, 4, 42] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], G2Matrix[A3, A4], S2F[A3, $A], S2SC2S[A4, A2]];
    BetaTensor[1, 4, 43] := BetaTensor[1, 4, 43] = Ttimes[C2G[A1, $A], G2Matrix[A1, A2], G2Matrix[A3, A4], S2FC2F[A2, A3], S2S[A4, $B]];
    BetaTensor[1, 4, 44] := BetaTensor[1, 4, 44] = Ttimes[C2G[A1, $A], G2Matrix[A1, A2], G2Matrix[A3, A4], S2S[A3, $B], S2SC2S[A4, A2]];
    BetaTensor[1, 4, 45] := BetaTensor[1, 4, 45] = Ttimes[C2G[A1, A2], C2G[A3, $B], C2G[A4, A5], C2G[A6, $A], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6]];
    BetaTensor[1, 4, 46] := BetaTensor[1, 4, 46] = Ttimes[C2G[A1, A2], C2G[A3, A4], C2G[A5, $B], G2Matrix[A1, A3], G2Matrix[A2, A5], G2Matrix[A4, A6], S2F[A6, $A]];
    BetaTensor[1, 4, 47] := BetaTensor[1, 4, 47] = Ttimes[C2G[A1, A2], C2G[A3, A4], C2G[A5, $A], G2Matrix[A1, A3], G2Matrix[A2, A5], G2Matrix[A4, A6], S2S[A6, $B]];
    BetaTensor[1, 4, 48] := BetaTensor[1, 4, 48] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], S2F[A2, A3], S2F[A4, A5], S2F[A6, $A]];
    BetaTensor[1, 4, 49] := BetaTensor[1, 4, 49] = Ttimes[C2G[A1, A2], C2G[A3, $B], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2F[A4, A5], S2F[A6, $A]];
    BetaTensor[1, 4, 50] := BetaTensor[1, 4, 50] = Ttimes[C2G[A1, $A], G2Matrix[A1, A2], G2Matrix[A3, A4], G2Matrix[A5, A6], S2S[A2, A3], S2S[A4, A5], S2S[A6, $B]];
    BetaTensor[1, 4, 51] := BetaTensor[1, 4, 51] = Ttimes[C2G[A1, A2], C2G[A3, $A], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2S[A4, A5], S2S[A6, $B]];
    BetaTensor[1, 4, 52] := BetaTensor[1, 4, 52] = Ttimes[C2G[A1, A2], C2G[A3, $A], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2F[A4, A5], S2S[A6, $B]];
    BetaTensor[1, 4, 53] := BetaTensor[1, 4, 53] = Ttimes[C2G[A1, A2], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2F[A3, A5], S2F[A4, $A], S2S[A6, $B]];
    BetaTensor[1, 4, 54] := BetaTensor[1, 4, 54] = Ttimes[C2G[A1, A2], G2Matrix[A1, A3], G2Matrix[A2, A4], G2Matrix[A5, A6], S2F[A3, $A], S2S[A4, A5], S2S[A6, $B]];
    BetaTensor[1, 4, 55] := BetaTensor[1, 4, 55] = Ttimes[FGauge[A1, A2, A3], FGauge[$A, A4, A5], G2Matrix[A1, A4], Lam[b1, b2, b3, b4], Tscal[$B, b5, b4], TsG2[A2, b1, b5], TsG2[A3, b2, b6], TsG2[A5, b6, b3]];
    BetaTensor[1, 4, 56] := BetaTensor[1, 4, 56] = Ttimes[C2S[b1, b2], Lam[b2, b3, b4, b5], Tscal[A1, b5, b6], Tscal[$A, b6, b1], Tscal[$B, b7, b4], TsG2[A1, b3, b7]];
    BetaTensor[1, 4, 57] := BetaTensor[1, 4, 57] = Ttimes[Lam[b1, b2, b3, b4], Tscal[A1, b5, b3], Tscal[A2, b4, b6], Tscal[$A, b7, b8], Tscal[$B, b8, b2], TsG2[A1, b6, b7], TsG2[A2, b1, b5]];
    BetaTensor[1, 4, 58] := BetaTensor[1, 4, 58] = Ttimes[C2G[A1, A2], Lam[b1, b2, b3, b4], Tscal[$A, b5, b3], Tscal[$B, b6, b4], TsG2[A1, b2, b5], TsG2[A2, b1, b6]];
    BetaTensor[1, 4, 59] := BetaTensor[1, 4, 59] = Ttimes[Lam[b1, b2, b3, b4], S2S[A1, A2], Tscal[$A, b1, b5], Tscal[$B, b2, b6], TsG2[A1, b5, b4], TsG2[A2, b6, b3]];
    BetaTensor[1, 4, 60] := BetaTensor[1, 4, 60] = Ttimes[C2G[A1, $A], Lam[b1, b2, b3, b4], Tscal[A2, b2, b5], Tscal[$B, b5, b3], TsG2[A1, b1, b6], TsG2[A2, b6, b4]];
    BetaTensor[1, 4, 61] := BetaTensor[1, 4, 61] = Ttimes[Lam[b1, b2, b3, b4], S2F[A1, A2], Tscal[$A, b5, b3], Tscal[$B, b6, b4], TsG2[A1, b2, b5], TsG2[A2, b1, b6]];
    BetaTensor[1, 4, 62] := BetaTensor[1, 4, 62] = Ttimes[Lam[b1, b2, b3, b4], Lam[b3, b5, b6, b4], Tscal[A1, b7, b2], Tscal[$A, b5, b7], Tscal[$B, b8, b6], TsG2[A1, b1, b8]];
    BetaTensor[1, 4, 63] := BetaTensor[1, 4, 63] = Ttimes[Lam[b1, b2, b3, b4], Lam[b3, b5, b6, b4], Tscal[A1, b7, b6], Tscal[$A, b5, b7], Tscal[$B, b8, b2], TsG2[A1, b1, b8]];
    BetaTensor[1, 4, 64] := BetaTensor[1, 4, 64] = Ttimes[C2S[b1, b2], Lam2[b3, b2], Tscal[$A, b4, b1], Tscal[$B, b3, b4]];
    BetaTensor[1, 4, 65] := BetaTensor[1, 4, 65] = Ttimes[C2S[b1, b2], Lam[b1, b3, b4, b5], Lam[b2, b4, b6, b5], Tscal[$A, b3, b7], Tscal[$B, b7, b6]];
    BetaTensor[1, 4, 66] := BetaTensor[1, 4, 66] = Ttimes[C2G[A1, $A], Lam2[b1, b2], Tscal[$B, b1, b3], TsG2[A1, b3, b2]];
    BetaTensor[1, 4, 67] := BetaTensor[1, 4, 67] = Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b5, b6, b3], Lam[b2, b7, b8, b4], Tscal[$A, b7, b5], Tscal[$B, b8, b6]];
    BetaTensor[1, 4, 68] := BetaTensor[1, 4, 68] = Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b5, b6, b3], Lam[b2, b6, b7, b4], Tscal[$A, b7, b8], Tscal[$B, b8, b5]];
    BetaTensor[1, 4, 69] := BetaTensor[1, 4, 69] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], Tr[Tdot[TfermTil[A1, k1, k2], TfermTil[A3, k2, k3], TfermTil[$A, k3, k4], TfermTil[A2, k4, k5], Yuk[b1, k6, k5], Tferm[$B, k6, k7], Tferm[A4, k7, k8], YukTil[b1, k8, k1]]]];
    BetaTensor[1, 4, 70] := BetaTensor[1, 4, 70] = Ttimes[S2FY2F[A1, A2], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A1, b4, b1], TsG2[A2, b3, b4]];
    BetaTensor[1, 4, 71] := BetaTensor[1, 4, 71] = Ttimes[S2SY2S[A1, A2], Tscal[$A, b1, b2], Tscal[$B, b2, b3], TsG2[A1, b4, b1], TsG2[A2, b3, b4]];
    BetaTensor[1, 4, 72] := BetaTensor[1, 4, 72] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], S2FY2F[A1, A3], Tr[Tdot[TfermTil[$B, k1, k2], TfermTil[A4, k2, k3], TfermTil[A2, k3, k4], TfermTil[$A, k4, k1]]]];
    BetaTensor[1, 4, 73] := BetaTensor[1, 4, 73] = Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], S2SY2S[A1, A3], Tr[Tdot[TfermTil[A2, k1, k2], TfermTil[A4, k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]]];
    BetaTensor[1, 4, 74] := BetaTensor[1, 4, 74] = Tr[Tdot[C2Ft[k1, k2], C2Ft[k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k5], Y2F[k5, k1]]];
    BetaTensor[1, 4, 75] := BetaTensor[1, 4, 75] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k4, k3], C2F[k4, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 76] := BetaTensor[1, 4, 76] = Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], Yuk[b1, k4, k3], C2F[k4, k5], YukTil[b1, k5, k6], TfermTil[$A, k6, k1]]];
    BetaTensor[1, 4, 77] := BetaTensor[1, 4, 77] = Ttimes[C2S[b1, b2], C2S[b2, b3], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]]];
    BetaTensor[1, 4, 78] := BetaTensor[1, 4, 78] = Tr[Tdot[C2Ft[k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k6], YukTil[b1, k1, k6]]];
    BetaTensor[1, 4, 79] := BetaTensor[1, 4, 79] = Ttimes[Tr[Tdot[C2Ft[k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k1, k4]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 80] := BetaTensor[1, 4, 80] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k2]]]];
    BetaTensor[1, 4, 81] := BetaTensor[1, 4, 81] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b3, k3, k2]]], Tscal[$A, b4, b3], Tscal[$B, b2, b4]];
    BetaTensor[1, 4, 82] := BetaTensor[1, 4, 82] = Ttimes[C2S[b1, b2], Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k1]]]];
    BetaTensor[1, 4, 83] := BetaTensor[1, 4, 83] = Ttimes[C2S[b1, b2], C2S[b2, b3], Tscal[$A, b4, b5], Tscal[$B, b5, b1], Y2S[b4, b3]];
    BetaTensor[1, 4, 84] := BetaTensor[1, 4, 84] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], C2FC2Gt[k3, k4], TfermTil[$A, k4, k1]]];
    BetaTensor[1, 4, 85] := BetaTensor[1, 4, 85] = Tr[Tdot[C2FS2Ft[k1, k2], Y2F[k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]];
    BetaTensor[1, 4, 86] := BetaTensor[1, 4, 86] = Tr[Tdot[C2FS2St[k1, k2], Y2F[k2, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k1]]];
    BetaTensor[1, 4, 87] := BetaTensor[1, 4, 87] = Tr[Tdot[YukTil[b1, k1, k2], C2FC2Gt[k1, k3], Yuk[b1, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k2]]];
    BetaTensor[1, 4, 88] := BetaTensor[1, 4, 88] = Ttimes[C2SC2G[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]]];
    BetaTensor[1, 4, 89] := BetaTensor[1, 4, 89] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2FC2Gt[k1, k3], Yuk[b2, k3, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 90] := BetaTensor[1, 4, 90] = Tr[Tdot[YukTil[b1, k1, k2], C2FS2Ft[k2, k3], Yuk[b1, k4, k3], Tferm[$A, k4, k5], Tferm[$B, k5, k1]]];
    BetaTensor[1, 4, 91] := BetaTensor[1, 4, 91] = Ttimes[C2SS2F[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]]];
    BetaTensor[1, 4, 92] := BetaTensor[1, 4, 92] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2FS2Ft[k1, k3], Yuk[b2, k3, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 93] := BetaTensor[1, 4, 93] = Tr[Tdot[YukTil[b1, k1, k2], C2FS2St[k2, k3], Yuk[b1, k4, k3], Tferm[$A, k4, k5], Tferm[$B, k5, k1]]];
    BetaTensor[1, 4, 94] := BetaTensor[1, 4, 94] = Ttimes[C2SS2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]]];
    BetaTensor[1, 4, 95] := BetaTensor[1, 4, 95] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2FS2St[k1, k3], Yuk[b2, k3, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 96] := BetaTensor[1, 4, 96] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[TfermTil[$A, k1, k2], C2Ft[k2, k3], Y2F[k3, k4], TfermTil[A2, k4, k1]]]];
    BetaTensor[1, 4, 97] := BetaTensor[1, 4, 97] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], Tferm[$A, k4, k5], Tferm[A2, k5, k2]]]];
    BetaTensor[1, 4, 98] := BetaTensor[1, 4, 98] = Ttimes[C2G[A1, $A], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k2]]], Tscal[$B, b3, b2], TsG2[A1, b1, b3]];
    BetaTensor[1, 4, 99] := BetaTensor[1, 4, 99] = Ttimes[C2G[A1, $B], C2S[b1, b2], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A2, k1, k3], TfermTil[$A, k3, k4], Yuk[b2, k4, k2]]]];
    BetaTensor[1, 4, 100] := BetaTensor[1, 4, 100] = Ttimes[C2SC2G[b1, b2], Tscal[$A, b3, b4], Tscal[$B, b4, b1], Y2S[b3, b2]];
    BetaTensor[1, 4, 101] := BetaTensor[1, 4, 101] = Ttimes[C2SS2F[b1, b2], Tscal[$A, b3, b4], Tscal[$B, b4, b1], Y2S[b3, b2]];
    BetaTensor[1, 4, 102] := BetaTensor[1, 4, 102] = Ttimes[C2SS2S[b1, b2], Tscal[$A, b3, b4], Tscal[$B, b2, b3], Y2S[b1, b4]];
    BetaTensor[1, 4, 103] := BetaTensor[1, 4, 103] = Ttimes[C2G[A1, $A], C2S[b1, b2], Tscal[$B, b3, b4], TsG2[A1, b4, b1], Y2S[b3, b2]];
    BetaTensor[1, 4, 104] := BetaTensor[1, 4, 104] = Ttimes[C2G[A1, A2], C2G[A3, $B], G2Matrix[A1, A3], G2Matrix[A2, A4], S2FY2F[A4, $A]];
    BetaTensor[1, 4, 105] := BetaTensor[1, 4, 105] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], G2Matrix[A3, A4], S2F[A3, $A], S2FY2F[A2, A4]];
    BetaTensor[1, 4, 106] := BetaTensor[1, 4, 106] = Ttimes[C2G[A1, $A], G2Matrix[A1, A2], G2Matrix[A3, A4], S2FY2F[A2, A3], S2S[A4, $B]];
    BetaTensor[1, 4, 107] := BetaTensor[1, 4, 107] = Ttimes[C2G[A1, A2], C2G[A3, $A], G2Matrix[A1, A3], G2Matrix[A2, A4], S2SY2S[A4, $B]];
    BetaTensor[1, 4, 108] := BetaTensor[1, 4, 108] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], G2Matrix[A3, A4], S2F[A3, $A], S2SY2S[A2, A4]];
    BetaTensor[1, 4, 109] := BetaTensor[1, 4, 109] = Ttimes[C2G[A1, $A], G2Matrix[A1, A2], G2Matrix[A3, A4], S2S[A3, $B], S2SY2S[A2, A4]];
    BetaTensor[1, 4, 110] := BetaTensor[1, 4, 110] = Ttimes[Lam[b1, b2, b3, b4], Tr[Tdot[YukTil[b3, k1, k2], TfermTil[A1, k1, k3], TfermTil[$A, k3, k4], Yuk[b4, k4, k2]]], Tscal[$B, b5, b2], TsG2[A1, b1, b5]];
    BetaTensor[1, 4, 111] := BetaTensor[1, 4, 111] = Ttimes[Lam[b1, b2, b3, b4], Tscal[A1, b4, b5], Tscal[$A, b5, b6], Tscal[$B, b7, b3], TsG2[A1, b2, b7], Y2S[b6, b1]];
    BetaTensor[1, 4, 112] := BetaTensor[1, 4, 112] = Ttimes[Lam2[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]]];
    BetaTensor[1, 4, 113] := BetaTensor[1, 4, 113] = Ttimes[Lam2[b1, b2], Tscal[$A, b3, b4], Tscal[$B, b2, b3], Y2S[b4, b1]];
    BetaTensor[1, 4, 114] := BetaTensor[1, 4, 114] = Ttimes[Lam[b1, b2, b3, b4], Lam[b3, b5, b6, b4], Tscal[$A, b2, b7], Tscal[$B, b7, b6], Y2S[b1, b5]];
    BetaTensor[1, 4, 115] := BetaTensor[1, 4, 115] = Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[$A, k3, k4], Yuk[b1, k4, k5], YukTil[b2, k6, k5], TfermTil[$B, k6, k7], TfermTil[A2, k7, k8], Yuk[b2, k8, k2]]]];
    BetaTensor[1, 4, 116] := BetaTensor[1, 4, 116] = Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[$A, k3, k4], Yuk[b2, k4, k2]]], Tr[Tdot[YukTil[b1, k5, k6], TfermTil[A2, k5, k7], TfermTil[$B, k7, k8], Yuk[b2, k8, k6]]]];
    BetaTensor[1, 4, 117] := BetaTensor[1, 4, 117] = Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[$A, k3, k4], Yuk[b2, k4, k5], Tferm[A2, k5, k6], YukTil[b1, k6, k7], Yuk[b2, k8, k7], Tferm[$B, k8, k2]]]];
    BetaTensor[1, 4, 118] := BetaTensor[1, 4, 118] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k2], Tferm[$A, k3, k4], YukTil[b3, k4, k5], Yuk[b2, k6, k5], Tferm[$B, k6, k1]]]];
    BetaTensor[1, 4, 119] := BetaTensor[1, 4, 119] = Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[b2, k3, k4], YukTil[b3, k5, k4], TfermTil[A2, k5, k6], Yuk[b3, k6, k2]]], Tscal[$A, b1, b4], Tscal[$B, b4, b2]];
    BetaTensor[1, 4, 120] := BetaTensor[1, 4, 120] = Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], YukTil[b2, k5, k6], Yuk[b1, k7, k6], Tferm[$B, k7, k2]]];
    BetaTensor[1, 4, 121] := BetaTensor[1, 4, 121] = Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k6], Yuk[b2, k6, k7], YukTil[b1, k1, k7]]];
    BetaTensor[1, 4, 122] := BetaTensor[1, 4, 122] = Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[b1, k4, k3], Tferm[$A, k4, k5], Tferm[$B, k5, k6], YukTil[b2, k6, k7], TfermTil[A2, k7, k8], Yuk[b2, k8, k2]]]];
    BetaTensor[1, 4, 123] := BetaTensor[1, 4, 123] = Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[$A, k3, k4], TfermTil[$B, k4, k5], Yuk[b1, k5, k6], YukTil[b2, k7, k6], TfermTil[A2, k7, k8], Yuk[b2, k8, k2]]]];
    BetaTensor[1, 4, 124] := BetaTensor[1, 4, 124] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k1], YukTil[b2, k3, k4], Yuk[b4, k4, k2]]], Tscal[$A, b4, b5], Tscal[$B, b5, b3]];
    BetaTensor[1, 4, 125] := BetaTensor[1, 4, 125] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k1], YukTil[b4, k3, k4], Yuk[b3, k4, k2]]], Tscal[$A, b5, b4], Tscal[$B, b2, b5]];
    BetaTensor[1, 4, 126] := BetaTensor[1, 4, 126] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], YukTil[b3, k4, k5], Yuk[b2, k2, k5]]], Tscal[$A, b3, b4], Tscal[$B, b4, b1]];
    BetaTensor[1, 4, 127] := BetaTensor[1, 4, 127] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b2, k5, k6], Yuk[b3, k6, k2]]]];
    BetaTensor[1, 4, 128] := BetaTensor[1, 4, 128] = Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], YukTil[b1, k5, k4], TfermTil[$A, k5, k6], TfermTil[$B, k6, k7], Yuk[b2, k7, k2]]];
    BetaTensor[1, 4, 129] := BetaTensor[1, 4, 129] = Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k6], YukTil[b1, k6, k7], Yuk[b2, k7, k2]]];
    BetaTensor[1, 4, 130] := BetaTensor[1, 4, 130] = Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k6, k5], YukTil[b2, k7, k6], TfermTil[$A, k7, k1]]];
    BetaTensor[1, 4, 131] := BetaTensor[1, 4, 131] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k2]]], Tr[Tdot[YukTil[b2, k4, k5], Yuk[b1, k6, k5], Tferm[$A, k6, k7], Tferm[$B, k7, k4]]]];
    BetaTensor[1, 4, 132] := BetaTensor[1, 4, 132] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]], Y2S[b3, b2]];
    BetaTensor[1, 4, 133] := BetaTensor[1, 4, 133] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b4], Y2S[b2, b4]];
    BetaTensor[1, 4, 134] := BetaTensor[1, 4, 134] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k2]]], Y2S[b1, b2]];
    BetaTensor[1, 4, 135] := BetaTensor[1, 4, 135] = Ttimes[Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k1]]], Y2S[b1, b2]];
    BetaTensor[1, 4, 136] := BetaTensor[1, 4, 136] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b3, k3, k2]]], Tscal[$A, b4, b3], Tscal[$B, b2, b4]];
    BetaTensor[1, 4, 137] := BetaTensor[1, 4, 137] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k1], YukTil[b4, k3, k4], Yuk[b2, k4, k2]]], Tscal[$A, b4, b5], Tscal[$B, b5, b3]];
    BetaTensor[1, 4, 138] := BetaTensor[1, 4, 138] = Ttimes[Tr[Tdot[Y2Ft[k1, k2], C2F[k2, k3], YukTil[b1, k3, k4], Yuk[b2, k1, k4]]], Tscal[$A, b2, b3], Tscal[$B, b3, b1]];
    BetaTensor[1, 4, 139] := BetaTensor[1, 4, 139] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], Y2Ft[k4, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 140] := BetaTensor[1, 4, 140] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k2]]]];
    BetaTensor[1, 4, 141] := BetaTensor[1, 4, 141] = Ttimes[C2S[b1, b2], Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k1]]]];
    BetaTensor[1, 4, 142] := BetaTensor[1, 4, 142] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[b3, k2, k5]]], Tscal[$A, b2, b4], Tscal[$B, b4, b3]];
    BetaTensor[1, 4, 143] := BetaTensor[1, 4, 143] = Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b3, k5, k6], Yuk[b2, k6, k2]]]];
    BetaTensor[1, 4, 144] := BetaTensor[1, 4, 144] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], C2Ft[k3, k4], Y2F[k4, k5], TfermTil[$A, k5, k1]]];
    BetaTensor[1, 4, 145] := BetaTensor[1, 4, 145] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Yuk[b1, k4, k3], C2F[k4, k5], YukTil[b1, k5, k6], TfermTil[$A, k6, k1]]];
    BetaTensor[1, 4, 146] := BetaTensor[1, 4, 146] = Tr[Tdot[Y2Ft[k1, k2], C2F[k2, k3], YukTil[b1, k3, k4], TfermTil[$A, k4, k5], TfermTil[$B, k5, k6], Yuk[b1, k1, k6]]];
    BetaTensor[1, 4, 147] := BetaTensor[1, 4, 147] = Tr[Tdot[TfermTil[$B, k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], Y2Ft[k4, k5], YukTil[b1, k5, k6], TfermTil[$A, k6, k1]]];
    BetaTensor[1, 4, 148] := BetaTensor[1, 4, 148] = Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k6], TfermTil[$B, k6, k7], Yuk[b2, k7, k2]]];
    BetaTensor[1, 4, 149] := BetaTensor[1, 4, 149] = Ttimes[C2S[b1, b2], Tscal[$A, b3, b4], Tscal[$B, b5, b3], Y2S[b4, b2], Y2S[b5, b1]];
    BetaTensor[1, 4, 150] := BetaTensor[1, 4, 150] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A2, k1, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k6], Yuk[b2, k6, k2]]]];
    BetaTensor[1, 4, 151] := BetaTensor[1, 4, 151] = Ttimes[C2G[A1, $A], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k4, k2]]], Tscal[$B, b4, b3], TsG2[A1, b1, b4]];
    BetaTensor[1, 4, 152] := BetaTensor[1, 4, 152] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A2, k1, k3], TfermTil[$A, k3, k4], Yuk[b2, k4, k5], YukTil[b1, k5, k6], Yuk[b2, k2, k6]]]];
    BetaTensor[1, 4, 153] := BetaTensor[1, 4, 153] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A2, k1, k3], TfermTil[$A, k3, k4], Yuk[b2, k4, k2]]], Y2S[b1, b2]];
    BetaTensor[1, 4, 154] := BetaTensor[1, 4, 154] = Ttimes[C2G[A1, $A], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k2]]], Tscal[$B, b3, b2], TsG2[A1, b1, b3]];
    BetaTensor[1, 4, 155] := BetaTensor[1, 4, 155] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[TfermTil[$A, k1, k2], Y2F[k2, k3], Y2F[k3, k4], TfermTil[A2, k4, k1]]]];
    BetaTensor[1, 4, 156] := BetaTensor[1, 4, 156] = Ttimes[C2G[A1, $B], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k4], Tferm[$A, k4, k5], Tferm[A2, k5, k2]]]];
    BetaTensor[1, 4, 157] := BetaTensor[1, 4, 157] = Ttimes[C2G[A1, $A], Tscal[$B, b1, b2], TsG2[A1, b2, b3], Y2S[b4, b1], Y2S[b4, b3]];
    BetaTensor[1, 4, 158] := BetaTensor[1, 4, 158] = Ttimes[Lam[b1, b2, b3, b4], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k2], Tferm[$A, k3, k4], YukTil[b2, k4, k5], Yuk[b3, k6, k5], Tferm[$B, k6, k1]]]];
    BetaTensor[1, 4, 159] := BetaTensor[1, 4, 159] = Ttimes[Lam[b1, b2, b3, b4], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b5, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k4, k2]]], Tscal[$A, b5, b6], Tscal[$B, b6, b1]];
    BetaTensor[1, 4, 160] := BetaTensor[1, 4, 160] = Ttimes[Lam[b1, b2, b3, b4], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b2, k5, k6], Yuk[b3, k6, k2]]]];
    BetaTensor[1, 4, 161] := BetaTensor[1, 4, 161] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], YukTil[b3, k4, k5], Yuk[b1, k5, k6], YukTil[b2, k7, k6], TfermTil[$B, k7, k8], Yuk[b3, k8, k2]]];
    BetaTensor[1, 4, 162] := BetaTensor[1, 4, 162] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b4, k4, k5], YukTil[b2, k6, k5], Yuk[b3, k6, k2]]], Tscal[$A, b1, b5], Tscal[$B, b5, b4]];
    BetaTensor[1, 4, 163] := BetaTensor[1, 4, 163] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b3, k5, k6], Yuk[b1, k6, k7], YukTil[b2, k7, k8], Yuk[b3, k2, k8]]];
    BetaTensor[1, 4, 164] := BetaTensor[1, 4, 164] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], YukTil[b1, k4, k5], Yuk[b3, k6, k5], Tferm[$B, k6, k7], YukTil[b3, k7, k8], Yuk[b2, k8, k2]]];
    BetaTensor[1, 4, 165] := BetaTensor[1, 4, 165] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], YukTil[b3, k4, k5], Yuk[b2, k5, k6], YukTil[b1, k7, k6], Yuk[b3, k8, k7], Tferm[$B, k8, k1]]];
    BetaTensor[1, 4, 166] := BetaTensor[1, 4, 166] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k5, k6], YukTil[b3, k7, k6], TfermTil[$B, k7, k8], Yuk[b3, k8, k2]]];
    BetaTensor[1, 4, 167] := BetaTensor[1, 4, 167] = Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], YukTil[b2, k5, k6], Yuk[b1, k7, k6], Tferm[$B, k7, k2]]];
    BetaTensor[1, 4, 168] := BetaTensor[1, 4, 168] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k6], Yuk[b2, k6, k7], YukTil[b1, k1, k7]]];
    BetaTensor[1, 4, 169] := BetaTensor[1, 4, 169] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], YukTil[b2, k4, k5], Yuk[b3, k6, k5], Tferm[$B, k6, k1]]], Y2S[b1, b3]];
    BetaTensor[1, 4, 170] := BetaTensor[1, 4, 170] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b4, k5, k4], YukTil[b2, k6, k5], Yuk[b4, k6, k2]]], Tscal[$A, b1, b5], Tscal[$B, b5, b3]];
    BetaTensor[1, 4, 171] := BetaTensor[1, 4, 171] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b4, k5, k4], YukTil[b3, k6, k5], Yuk[b2, k6, k2]]], Tscal[$A, b1, b5], Tscal[$B, b5, b4]];
    BetaTensor[1, 4, 172] := BetaTensor[1, 4, 172] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b1, k5, k6], Yuk[b3, k6, k7], YukTil[b2, k8, k7], Yuk[b3, k2, k8]]];
    BetaTensor[1, 4, 173] := BetaTensor[1, 4, 173] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b3, k5, k6], Yuk[b1, k6, k7], YukTil[b3, k7, k8], Yuk[b2, k8, k2]]];
    BetaTensor[1, 4, 174] := BetaTensor[1, 4, 174] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k4, k2]]], Tr[Tdot[YukTil[b3, k5, k6], Yuk[b1, k7, k6], Tferm[$A, k7, k8], Tferm[$B, k8, k5]]]];
    BetaTensor[1, 4, 175] := BetaTensor[1, 4, 175] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], YukTil[b3, k5, k4], Yuk[b2, k5, k2]]], Tscal[$A, b1, b4], Tscal[$B, b4, b3]];
    BetaTensor[1, 4, 176] := BetaTensor[1, 4, 176] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[b3, k5, k6], Yuk[b4, k6, k2]]], Tscal[$A, b1, b5], Tscal[$B, b5, b4]];
    BetaTensor[1, 4, 177] := BetaTensor[1, 4, 177] = Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], YukTil[b1, k5, k4], TfermTil[$A, k5, k6], TfermTil[$B, k6, k7], Yuk[b2, k7, k2]]];
    BetaTensor[1, 4, 178] := BetaTensor[1, 4, 178] = Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k6], YukTil[b1, k6, k7], Yuk[b2, k7, k2]]];
    BetaTensor[1, 4, 179] := BetaTensor[1, 4, 179] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], Yuk[b1, k6, k5], YukTil[b2, k7, k6], TfermTil[$A, k7, k1]]];
    BetaTensor[1, 4, 180] := BetaTensor[1, 4, 180] = Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b2, k5, k6], Yuk[b3, k6, k7], YukTil[b1, k8, k7], Yuk[b3, k8, k2]]];
    BetaTensor[1, 4, 181] := BetaTensor[1, 4, 181] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k4, k2]]], Tscal[$A, b1, b4], Tscal[$B, b4, b5], Y2S[b3, b5]];
    BetaTensor[1, 4, 182] := BetaTensor[1, 4, 182] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b4, k4, k2]]], Tscal[$A, b1, b5], Tscal[$B, b5, b3], Y2S[b4, b2]];
    BetaTensor[1, 4, 183] := BetaTensor[1, 4, 183] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b3, k5, k6], Yuk[b2, k6, k2]]], Y2S[b1, b3]];
    BetaTensor[1, 4, 184] := BetaTensor[1, 4, 184] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k2]]], Tr[Tdot[YukTil[b2, k4, k5], Yuk[b1, k6, k5], Tferm[$A, k6, k7], Tferm[$B, k7, k4]]]];
    BetaTensor[1, 4, 185] := BetaTensor[1, 4, 185] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k4, k3], Y2Ft[k4, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 186] := BetaTensor[1, 4, 186] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[b3, k2, k5]]], Tscal[$A, b2, b4], Tscal[$B, b4, b3]];
    BetaTensor[1, 4, 187] := BetaTensor[1, 4, 187] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Yuk[b1, k4, k3], Y2Ft[k4, k5], YukTil[b1, k5, k6], TfermTil[$A, k6, k1]]];
    BetaTensor[1, 4, 188] := BetaTensor[1, 4, 188] = Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k6], TfermTil[$B, k6, k7], Yuk[b2, k7, k2]]];
    BetaTensor[1, 4, 189] := BetaTensor[1, 4, 189] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b4], Y2S[b2, b4]];
    BetaTensor[1, 4, 190] := BetaTensor[1, 4, 190] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b4, k4, k2]]], Tscal[$A, b1, b5], Tscal[$B, b5, b4], Y2S[b3, b2]];
    BetaTensor[1, 4, 191] := BetaTensor[1, 4, 191] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], Tferm[$A, k4, k5], Tferm[$B, k5, k2]]], Y2S[b1, b2]];
    BetaTensor[1, 4, 192] := BetaTensor[1, 4, 192] = Ttimes[Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfermTil[$A, k5, k1]]], Y2S[b1, b2]];
    BetaTensor[1, 4, 193] := BetaTensor[1, 4, 193] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], Tferm[$A, k3, k4], Tferm[$B, k4, k5], YukTil[b2, k5, k6], Yuk[b3, k6, k2]]], Y2S[b1, b3]];
    BetaTensor[1, 4, 194] := BetaTensor[1, 4, 194] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[$A, k3, k4], Tferm[$B, k4, k1]]], Y2S[b1, b3], Y2S[b2, b3]];
    BetaTensor[1, 4, 195] := BetaTensor[1, 4, 195] = Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Y2F[k3, k4], Yuk[b2, k4, k2]]], Tscal[$A, b1, b3], Tscal[$B, b3, b2]];
    BetaTensor[1, 4, 196] := BetaTensor[1, 4, 196] = Tr[Tdot[TfermTil[$B, k1, k2], Y2F[k2, k3], Y2F[k3, k4], Y2F[k4, k5], TfermTil[$A, k5, k1]]];
    BetaTensor[1, 4, 197] := BetaTensor[1, 4, 197] = Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Y2F[k3, k4], Yuk[b1, k4, k5], Tferm[$A, k5, k6], Tferm[$B, k6, k2]]];
    BetaTensor[1, 4, 198] := BetaTensor[1, 4, 198] = Ttimes[Tscal[$A, b1, b2], Tscal[$B, b3, b1], Y2S[b4, b2], Y2S[b5, b3], Y2S[b5, b4]];
(* gamma_5 terms*)
    BetaTensor[1, 4, 199] := BetaTensor[1, 4, 199] = Tr@Tdot[Tferm[$A, k1, k2], TfG2[A1, k2, k3], sig3Til, YukTil[b1, k3, k4], Yuk[b2, k4, k1]] Tr@Tdot[Tferm[$B, k5, k6], Tferm[A1, k6, k7], sig3Til, YukTil[b1, k7, k8], Yuk[b2, k8, k5]] // Expand;
    BetaTensor[1, 4, 200] := BetaTensor[1, 4, 200] = Tr@Tdot[Tferm[$A, k1, k2], sig3Til, YukTil[b1, k2, k3], TfG2t[A1, k3, k4], Yuk[b2, k4, k1]] Tr@Tdot[Tferm[$B, k5, k6], sig3Til, YukTil[b1, k6, k7], TfermTil[A1, k7, k8], Yuk[b2, k8, k5]] // Expand;
    BetaTensor[1, 4, 201] := BetaTensor[1, 4, 201] = Tr@Tdot[Tferm[$A, k1, k2], TfG2[A1, k2, k3], sig3Til, YukTil[b1, k3, k4], Yuk[b2, k4, k1]] Tr@Tdot[Tferm[$B, k5, k6], Tferm[A1, k6, k7], YukTil[b2, k7, k8], sig3, Yuk[b1, k8, k5]] // Expand;
    BetaTensor[1, 4, 202] := BetaTensor[1, 4, 202] = Tr@Tdot[Tferm[$A, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k3, k4], sig3, Yuk[b2, k4, k1]] Tr@Tdot[Tferm[$B, k5, k6], YukTil[b1, k6, k7], TfermTil[A1, k7, k8], sig3, Yuk[b2, k8, k5]] // Expand;

(*####################################*)
(*-----------Yukawa Tensors-----------*)
(*####################################*)
    BetaTensor[2, 0, 1] := BetaTensor[2, 0, 1] = - Global`\[Epsilon] / 2 Yuk[$a, $i, $j];

    BetaTensor[2, 1, 1] := BetaTensor[2, 1, 1] = {{0, 0}, {0, 0}} (*Yuk[b2, $i, $j] C2S[$a, b2] // Expand*);
    BetaTensor[2, 1, 2] := BetaTensor[2, 1, 2] = Tdot[C2Ft[$i, k1], Yuk[$a, k1, $j]] // TsSym[$i, $j];
    BetaTensor[2, 1, 3] := BetaTensor[2, 1, 3] = Tdot[Yuk[b1, $i, k1], YukTil[$a, k1, k2], Yuk[b1, k2, $j]];
    BetaTensor[2, 1, 4] := BetaTensor[2, 1, 4] = Tdot[Y2F[$i, k1], Yuk[$a, k1, $j]] // TsSym[$i, $j];
    BetaTensor[2, 1, 5] := BetaTensor[2, 1, 5] = Yuk[b1, $i, $j] Y2S[b1, $a] // Expand;
(* y^1 terms*)
    BetaTensor[2, 2, 1] := BetaTensor[2, 2, 1] = Ttimes[Yuk[b3, $i, $j], C2S[b3, b2], C2S[b2, $a]];
    BetaTensor[2, 2, 2] := BetaTensor[2, 2, 2] = Tdot[C2Ft[$i, k1], Yuk[b2, k1, $j]] C2S[b2, $a] // Expand // TsSym[$i, $j];
    BetaTensor[2, 2, 3] := BetaTensor[2, 2, 3] = {{0, 0}, {0, 0}} (*Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2F[k2, $j]]*);
    BetaTensor[2, 2, 4] := BetaTensor[2, 2, 4] = Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Yuk[$a, k2, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 5] := BetaTensor[2, 2, 5] = Yuk[b2, $i, $j] C2SC2G[b2, $a] // Expand;
    BetaTensor[2, 2, 6] := BetaTensor[2, 2, 6] = Yuk[b2, $i, $j] C2SS2S[b2, $a] // Expand;
    BetaTensor[2, 2, 7] := BetaTensor[2, 2, 7] = Yuk[b2, $i, $j] C2SS2F[b2, $a] // Expand;
    BetaTensor[2, 2, 8] := BetaTensor[2, 2, 8] = C2FC2Gt[$i, k1]. Yuk[$a, k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[2, 2, 9] := BetaTensor[2, 2, 9] = C2FS2St[$i, k1]. Yuk[$a, k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[2, 2, 10] := BetaTensor[2, 2, 10] = C2FS2Ft[$i, k1]. Yuk[$a, k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[2, 2, 11] := BetaTensor[2, 2, 11] = Yuk[b2, $i, $j] Lam2[b2, $a] // Expand;
(*y^3 terms*)
    BetaTensor[2, 2, 12] := BetaTensor[2, 2, 12] = Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b2, k2, k3], TfG2t[A1, k3, k4], Yuk[b2, k4, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 13] := BetaTensor[2, 2, 13] = {{0, 0}, {0, 0}} (*Tdot[Y2F[$i, k1], TfermTil[A1, k1, k2], Yuk[$a, k2, k3], TfG2[A1, k3, $j]] // TsSym[$i, $j]*);
    BetaTensor[2, 2, 14] := BetaTensor[2, 2, 14] = Tdot[Yuk[b2, $i, k1], YukTil[$a, k1, k2], Yuk[b3, k2, $j]] C2S[b2, b3] // Expand;
    BetaTensor[2, 2, 15] := BetaTensor[2, 2, 15] = Tdot[Yuk[b2, $i, k1], YukTil[b3, k1, k2], Yuk[b2, k2, $j]] C2S[b3, $a] // Expand;
    BetaTensor[2, 2, 16] := BetaTensor[2, 2, 16] = Tdot[Yuk[b2, $i, k1], C2F[k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 17] := BetaTensor[2, 2, 17] = Tdot[C2Ft[$i, k1], Yuk[b2, k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 18] := BetaTensor[2, 2, 18] = Yuk[b2, $i, $j] Y2SC2F[b2, $a] // Expand;
    BetaTensor[2, 2, 19] := BetaTensor[2, 2, 19] = Yuk[$a, $i, k1].Y2FC2St[k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[2, 2, 20] := BetaTensor[2, 2, 20] = Yuk[$a, $i, k1].Y2FC2Ft[k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[2, 2, 21] := BetaTensor[2, 2, 21] = Tdot[Yuk[$a, $i, k1], Y2Ft[k1, k2], C2F[k2, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 22] := BetaTensor[2, 2, 22] = {{0, 0}, {0, 0}} (*Ttimes[Yuk[b3, $i, $j], Y2S[b3, b2], C2S[b2, $a]]*);
    BetaTensor[2, 2, 23] := BetaTensor[2, 2, 23] = Tdot[Yuk[b2, $i, k1], YukTil[b3, k1, k2], Yuk[b4, k2, $j]] Lam[$a, b2, b3, b4] // Expand;
(*y^5 terms*)
    BetaTensor[2, 2, 24] := BetaTensor[2, 2, 24] = Tdot[Yuk[b2, $i, k1], YukTil[b3, k1, k2], Yuk[$a, k2, k3], YukTil[b2, k3, k4], Yuk[b3, k4, $j]];
    BetaTensor[2, 2, 25] := BetaTensor[2, 2, 25] = {{0, 0}, {0, 0}} (*Tdot[Yuk[b2, $i, k1], YukTil[$a, k1, k2], Yuk[b3, k2, k3], YukTil[b2, k3, k4], Yuk[b3, k4, $j]] // TsSym[$i, $j]*);
    BetaTensor[2, 2, 26] := BetaTensor[2, 2, 26] = Tdot[Yuk[b3, $i, k1], YukTil[b2, k1, k2], Yuk[$a, k2, k3], YukTil[b2, k3, k4], Yuk[b3, k4, $j]];
    BetaTensor[2, 2, 27] := BetaTensor[2, 2, 27] = {{0, 0}, {0, 0}} (*Tdot[Yuk[$a, $i, k1], Y4cFt[k1, $j]] // TsSym[$i, $j]*);
    BetaTensor[2, 2, 28] := BetaTensor[2, 2, 28] = Yuk[b2, $i, $j] Y4cS[b2, $a] // Expand;
    BetaTensor[2, 2, 29] := BetaTensor[2, 2, 29] = Tdot[Yuk[b2, $i, k1], Y2Ft[k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 30] := BetaTensor[2, 2, 30] = Tdot[Yuk[$a, $i, k1], Y2FY2Ft[k1, $j]] // TsSym[$i, $j];
    BetaTensor[2, 2, 31] := BetaTensor[2, 2, 31] = Yuk[b2, $i, $j] Y2SY2F[b2, $a] // Expand;
    BetaTensor[2, 2, 32] := BetaTensor[2, 2, 32] = Tdot[Yuk[b2, $i, k1], YukTil[$a, k1, k2], Yuk[b3, k2, $j]] Y2S[b2, b3] // Expand;
    BetaTensor[2, 2, 33] := BetaTensor[2, 2, 33] = Tdot[Yuk[$a, $i, k1], Y2FY2St[k1, $j]] // TsSym[$i, $j];

 (* 3-loop *)
    BetaTensor[2, 3, 1] := BetaTensor[2, 3, 1] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], TfermTil[A2, k1, k2], Yuk[$a, k2, k3], Tferm[A3, k3, k4], TfG2[A2, k4, k5], TfG2[A1, k5, k6], TfG2[A3, k6, $j]];
    BetaTensor[2, 3, 2] := BetaTensor[2, 3, 2] = TsSym[$i, $j]@ Ttimes[S2SC2S[A1, A2], TsG2[A1, b1, b2], TsG2[A2, $a, b1], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 3] := BetaTensor[2, 3, 3] = TsSym[$i, $j]@ Ttimes[S2SC2S[A1, A2], Tdot[Yuk[$a, $i, k1], TfG2[A2, k1, k2], TfG2[A1, k2, $j]]];
    BetaTensor[2, 3, 4] := BetaTensor[2, 3, 4] = TsSym[$i, $j]@ Ttimes[S2FC2F[A1, A2], TsG2[A1, $a, b1], TsG2[A2, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 5] := BetaTensor[2, 3, 5] = TsSym[$i, $j]@ Ttimes[S2FC2F[A1, A2], Tdot[Yuk[$a, $i, k1], TfG2[A1, k1, k2], TfG2[A2, k2, $j]]];
    BetaTensor[2, 3, 6] := BetaTensor[2, 3, 6] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], C2S[b3, b1], C2S[$a, b3], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 7] := BetaTensor[2, 3, 7] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], C2S[b2, $a], Tdot[C2Ft[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 8] := BetaTensor[2, 3, 8] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2Ft[$i, k1], Yuk[b1, k2, k1], C2F[k2, $j]]];
    BetaTensor[2, 3, 9] := BetaTensor[2, 3, 9] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 10] := BetaTensor[2, 3, 10] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2F[k2, k3], C2F[k3, $j]];
    BetaTensor[2, 3, 11] := BetaTensor[2, 3, 11] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], C2Ft[k1, k2], C2Ft[k2, k3], Yuk[$a, $j, k3]];
    BetaTensor[2, 3, 12] := BetaTensor[2, 3, 12] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], C2SC2G[b2, b1], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 13] := BetaTensor[2, 3, 13] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], C2SS2S[b2, b1], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 14] := BetaTensor[2, 3, 14] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], C2SS2F[b2, b1], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 15] := BetaTensor[2, 3, 15] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2FC2Gt[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 16] := BetaTensor[2, 3, 16] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2FS2St[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 17] := BetaTensor[2, 3, 17] = TsSym[$i, $j]@ Ttimes[C2SC2G[b1, $a], Tdot[C2Ft[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 18] := BetaTensor[2, 3, 18] = TsSym[$i, $j]@ Ttimes[C2SS2S[b1, $a], Tdot[C2Ft[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 19] := BetaTensor[2, 3, 19] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2FS2Ft[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 20] := BetaTensor[2, 3, 20] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2FC2G[k2, $j]];
    BetaTensor[2, 3, 21] := BetaTensor[2, 3, 21] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2FS2S[k2, $j]];
    BetaTensor[2, 3, 22] := BetaTensor[2, 3, 22] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], C2FC2Gt[k1, k2], Yuk[$a, $j, k2]];
    BetaTensor[2, 3, 23] := BetaTensor[2, 3, 23] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], C2FS2St[k1, k2], Yuk[$a, $j, k2]];
    BetaTensor[2, 3, 24] := BetaTensor[2, 3, 24] = TsSym[$i, $j]@ Ttimes[C2SS2F[b1, $a], Tdot[C2Ft[$i, k1], Yuk[b1, k1, $j]]];
    BetaTensor[2, 3, 25] := BetaTensor[2, 3, 25] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], C2FS2F[k2, $j]];
    BetaTensor[2, 3, 26] := BetaTensor[2, 3, 26] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], C2FS2Ft[k1, k2], Yuk[$a, $j, k2]];
    BetaTensor[2, 3, 27] := BetaTensor[2, 3, 27] = TsSym[$i, $j]@ Ttimes[G2Matrix[A1, A2], S2S[A3, A1], S2S[A4, A2], TsG2[A3, $a, b1], TsG2[A4, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 28] := BetaTensor[2, 3, 28] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2S[A3, A4], TsG2[A2, $a, b1], TsG2[A4, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 29] := BetaTensor[2, 3, 29] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], C2G[A3, A4], G2Matrix[A1, A3], TsG2[A2, $a, b1], TsG2[A4, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 30] := BetaTensor[2, 3, 30] = TsSym[$i, $j]@ Ttimes[G2Matrix[A1, A2], S2F[A3, A1], S2S[A2, A4], TsG2[A3, $a, b1], TsG2[A4, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 31] := BetaTensor[2, 3, 31] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2F[A3, A4], TsG2[A2, $a, b1], TsG2[A4, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 32] := BetaTensor[2, 3, 32] = TsSym[$i, $j]@ Ttimes[G2Matrix[A1, A2], S2S[A3, A1], S2S[A4, A2], Tdot[Yuk[$a, $i, k1], TfG2[A3, k1, k2], TfG2[A4, k2, $j]]];
    BetaTensor[2, 3, 33] := BetaTensor[2, 3, 33] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2S[A3, A4], Tdot[Yuk[$a, $i, k1], TfG2[A4, k1, k2], TfG2[A2, k2, $j]]];
    BetaTensor[2, 3, 34] := BetaTensor[2, 3, 34] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], C2G[A3, A4], G2Matrix[A1, A3], Tdot[Yuk[$a, $i, k1], TfG2[A2, k1, k2], TfG2[A4, k2, $j]]];
    BetaTensor[2, 3, 35] := BetaTensor[2, 3, 35] = TsSym[$i, $j]@ Ttimes[G2Matrix[A1, A2], S2F[A1, A3], S2F[A2, A4], TsG2[A3, $a, b1], TsG2[A4, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 36] := BetaTensor[2, 3, 36] = TsSym[$i, $j]@ Ttimes[G2Matrix[A1, A2], S2F[A3, A1], S2S[A4, A2], Tdot[Yuk[$a, $i, k1], TfG2[A4, k1, k2], TfG2[A3, k2, $j]]];
    BetaTensor[2, 3, 37] := BetaTensor[2, 3, 37] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2F[A3, A4], Tdot[Yuk[$a, $i, k1], TfG2[A4, k1, k2], TfG2[A2, k2, $j]]];
    BetaTensor[2, 3, 38] := BetaTensor[2, 3, 38] = TsSym[$i, $j]@ Ttimes[G2Matrix[A1, A2], S2F[A3, A1], S2F[A4, A2], Tdot[Yuk[$a, $i, k1], TfG2[A3, k1, k2], TfG2[A4, k2, $j]]];
    BetaTensor[2, 3, 39] := BetaTensor[2, 3, 39] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[TfermTil[A1, $i, k1], TfermTil[A2, k1, k2], Yuk[b3, k2, $j]], TsG2[A1, b4, b2], TsG2[A2, b1, b4]];
    BetaTensor[2, 3, 40] := BetaTensor[2, 3, 40] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tscal[A1, b4, b2], Tscal[A2, b3, b5], TsG2[A1, b5, b6], TsG2[A2, b1, b4], Yuk[b6, $i, $j]];
    BetaTensor[2, 3, 41] := BetaTensor[2, 3, 41] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Lam2[b1, $a], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 42] := BetaTensor[2, 3, 42] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Lam[b1, b3, b4, $a], Lam[b2, b5, b3, b4], Yuk[b5, $i, $j]];
    BetaTensor[2, 3, 43] := BetaTensor[2, 3, 43] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b2, b5, b4, b6], Lam[b3, b5, b6, $a], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 44] := BetaTensor[2, 3, 44] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k2, k1], TfG2[A2, k2, k3], TfG2[A1, k3, k4], YukTil[$a, k4, k5], TfermTil[A2, k5, k6], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 45] := BetaTensor[2, 3, 45] = TsSym[$i, $j]@ Ttimes[Tdot[TfG2t[A1, $i, k1], TfG2t[A2, k1, k2], Yuk[b1, k2, $j]], Tr[Tdot[YukTil[b1, k3, k4], TfermTil[A2, k3, k5], TfermTil[A1, k5, k6], Yuk[$a, k6, k4]]]];
    BetaTensor[2, 3, 46] := BetaTensor[2, 3, 46] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], TfermTil[A2, k1, k2], Yuk[$a, k2, k3], YukTil[b1, k4, k3], TfG2t[A2, k4, k5], TfG2t[A1, k5, k6], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 47] := BetaTensor[2, 3, 47] = TsSym[$i, $j]@ Ttimes[S2FY2F[A1, A2], TsG2[A1, $a, b1], TsG2[A2, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 48] := BetaTensor[2, 3, 48] = TsSym[$i, $j]@ Ttimes[S2FY2F[A1, A2], Tdot[Yuk[$a, $i, k1], TfG2[A1, k1, k2], TfG2[A2, k2, $j]]];
    BetaTensor[2, 3, 49] := BetaTensor[2, 3, 49] = TsSym[$i, $j]@ Ttimes[S2SY2S[A1, A2], TsG2[A1, $a, b1], TsG2[A2, b1, b2], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 50] := BetaTensor[2, 3, 50] = TsSym[$i, $j]@ Ttimes[S2SY2S[A1, A2], Tdot[Yuk[$a, $i, k1], TfG2[A1, k1, k2], TfG2[A2, k2, $j]]];
    BetaTensor[2, 3, 51] := BetaTensor[2, 3, 51] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 52] := BetaTensor[2, 3, 52] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], YukTil[b2, k3, k2], TfG2t[A1, k3, k4], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 53] := BetaTensor[2, 3, 53] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k2, k1], C2F[k2, k3], YukTil[b1, k3, k4], TfG2t[A1, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 54] := BetaTensor[2, 3, 54] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], C2Ft[k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 55] := BetaTensor[2, 3, 55] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], Tferm[A1, k2, k3], YukTil[b1, k3, k4], Yuk[$a, k5, k4], TfG2[A1, k5, $j]];
    BetaTensor[2, 3, 56] := BetaTensor[2, 3, 56] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[$a, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k4, k3], TfermTil[A1, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 57] := BetaTensor[2, 3, 57] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], TfermTil[A1, k1, k2], Yuk[$a, k2, k3], YukTil[b1, k4, k3], TfG2t[A1, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 58] := BetaTensor[2, 3, 58] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Y2F[$i, k1], TfermTil[A1, k1, k2], Yuk[b1, k3, k2], TfG2[A1, k3, $j]]];
    BetaTensor[2, 3, 59] := BetaTensor[2, 3, 59] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], TfG2[A1, k2, k3], YukTil[b2, k3, k4], Yuk[b1, k4, $j]]];
    BetaTensor[2, 3, 60] := BetaTensor[2, 3, 60] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k4, k3], C2Ft[k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 61] := BetaTensor[2, 3, 61] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], TfermTil[A1, k1, k2], Yuk[$a, k2, k3], TfG2[A1, k3, k4], Y2Ft[k4, $j]];
    BetaTensor[2, 3, 62] := BetaTensor[2, 3, 62] = TsSym[$i, $j]@ Ttimes[-1, Tdot[C2Ft[$i, k1], TfermTil[A1, k1, k2], Yuk[b1, k2, $j]], TsG2[A1, b1, b2], Y2S[b2, $a]];
    BetaTensor[2, 3, 63] := BetaTensor[2, 3, 63] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], C2S[b2, b3], Tdot[Yuk[b3, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 64] := BetaTensor[2, 3, 64] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], C2S[b2, $a], Tdot[Yuk[b3, $i, k1], YukTil[b1, k1, k2], Yuk[b3, $j, k2]]];
    BetaTensor[2, 3, 65] := BetaTensor[2, 3, 65] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], C2S[b2, b3], Tdot[Yuk[b3, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 66] := BetaTensor[2, 3, 66] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b2, k1, $i], C2F[k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]]];
    BetaTensor[2, 3, 67] := BetaTensor[2, 3, 67] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Yuk[b2, k1, $i], C2F[k1, k2], YukTil[b1, k2, k3], Yuk[b2, k3, $j]]];
    BetaTensor[2, 3, 68] := BetaTensor[2, 3, 68] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[C2Ft[$i, k1], Yuk[b2, k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]]];
    BetaTensor[2, 3, 69] := BetaTensor[2, 3, 69] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2Ft[$i, k1], Yuk[b2, k1, k2], YukTil[b1, k2, k3], Yuk[b2, $j, k3]]];
    BetaTensor[2, 3, 70] := BetaTensor[2, 3, 70] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[$a, k2, k3], C2Ft[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 71] := BetaTensor[2, 3, 71] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], C2F[k2, k3], YukTil[$a, k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 72] := BetaTensor[2, 3, 72] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b1, k4, k3], C2F[k4, $j]];
    BetaTensor[2, 3, 73] := BetaTensor[2, 3, 73] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k2, k1], C2F[k2, k3], YukTil[$a, k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 74] := BetaTensor[2, 3, 74] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Yuk[b1, k2, k3], YukTil[$a, k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 75] := BetaTensor[2, 3, 75] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k3, k2], C2Ft[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 76] := BetaTensor[2, 3, 76] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k2]]], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 77] := BetaTensor[2, 3, 77] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k4, k3], C2F[k4, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 78] := BetaTensor[2, 3, 78] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[C2Ft[k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[$a, k1, k4]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 79] := BetaTensor[2, 3, 79] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], C2S[b2, b3], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b3, $j, k2]]];
    BetaTensor[2, 3, 80] := BetaTensor[2, 3, 80] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], C2Ft[k2, k3], Yuk[b2, k3, $j]]];
    BetaTensor[2, 3, 81] := BetaTensor[2, 3, 81] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[C2Ft[$i, k1], Y2F[k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 82] := BetaTensor[2, 3, 82] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[C2Ft[$i, k1], Yuk[b2, k1, k2], YukTil[b1, k3, k2], Yuk[$a, $j, k3]]];
    BetaTensor[2, 3, 83] := BetaTensor[2, 3, 83] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], C2Ft[k2, k3], C2Ft[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 84] := BetaTensor[2, 3, 84] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Y2F[k2, k3], Yuk[$a, $j, k3]];
    BetaTensor[2, 3, 85] := BetaTensor[2, 3, 85] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k2, k1], C2F[k2, k3], YukTil[b1, k3, k4], Yuk[$a, $j, k4]];
    BetaTensor[2, 3, 86] := BetaTensor[2, 3, 86] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[$a, k2, k1], C2F[k2, k3], Y2Ft[k3, $j]];
    BetaTensor[2, 3, 87] := BetaTensor[2, 3, 87] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], C2S[$a, b3], Y2S[b2, b3], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 88] := BetaTensor[2, 3, 88] = TsSym[$i, $j]@ Ttimes[C2S[$a, b1], Tdot[C2Ft[$i, k1], Yuk[b2, k1, $j]], Y2S[b2, b1]];
    BetaTensor[2, 3, 89] := BetaTensor[2, 3, 89] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], Tdot[TfG2t[A2, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b1, k4, $j]]];
    BetaTensor[2, 3, 90] := BetaTensor[2, 3, 90] = TsSym[$i, $j]@ Ttimes[S2S[A1, A2], Tdot[TfG2t[A2, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b1, k4, $j]]];
    BetaTensor[2, 3, 91] := BetaTensor[2, 3, 91] = TsSym[$i, $j]@ Ttimes[S2F[A1, A2], Tdot[TfG2t[A2, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b1, k4, $j]]];
    BetaTensor[2, 3, 92] := BetaTensor[2, 3, 92] = TsSym[$i, $j]@ Ttimes[C2G[A1, A2], Tdot[Y2F[$i, k1], TfG2t[A1, k1, k2], Yuk[$a, k3, k2], TfG2[A2, k3, $j]]];
    BetaTensor[2, 3, 93] := BetaTensor[2, 3, 93] = TsSym[$i, $j]@ Ttimes[S2S[A1, A2], Tdot[Y2F[$i, k1], TfG2t[A1, k1, k2], Yuk[$a, k3, k2], TfG2[A2, k3, $j]]];
    BetaTensor[2, 3, 94] := BetaTensor[2, 3, 94] = TsSym[$i, $j]@ Ttimes[S2F[A1, A2], Tdot[Y2F[$i, k1], TfG2t[A1, k1, k2], Yuk[$a, k3, k2], TfG2[A2, k3, $j]]];
    BetaTensor[2, 3, 95] := BetaTensor[2, 3, 95] = TsSym[$i, $j]@ Ttimes[C2SC2G[b1, b2], Tdot[Yuk[b2, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 96] := BetaTensor[2, 3, 96] = TsSym[$i, $j]@ Ttimes[C2SS2S[b1, b2], Tdot[Yuk[b2, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 97] := BetaTensor[2, 3, 97] = TsSym[$i, $j]@ Ttimes[C2SC2G[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 98] := BetaTensor[2, 3, 98] = TsSym[$i, $j]@ Ttimes[C2SS2S[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 99] := BetaTensor[2, 3, 99] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2FC2G[k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 100] := BetaTensor[2, 3, 100] = TsSym[$i, $j]@ Ttimes[C2SS2F[b1, b2], Tdot[Yuk[b2, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 101] := BetaTensor[2, 3, 101] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2FS2S[k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 102] := BetaTensor[2, 3, 102] = TsSym[$i, $j]@ Ttimes[C2SS2F[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 103] := BetaTensor[2, 3, 103] = TsSym[$i, $j]@ Tdot[C2FC2Gt[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 104] := BetaTensor[2, 3, 104] = TsSym[$i, $j]@ Tdot[C2FS2St[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 105] := BetaTensor[2, 3, 105] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2FS2F[k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 106] := BetaTensor[2, 3, 106] = TsSym[$i, $j]@ Tdot[C2FS2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 107] := BetaTensor[2, 3, 107] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2FC2Gt[k1, k3], Yuk[$a, k3, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 108] := BetaTensor[2, 3, 108] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2FS2St[k1, k3], Yuk[$a, k3, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 109] := BetaTensor[2, 3, 109] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2FS2Ft[k1, k3], Yuk[$a, k3, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 110] := BetaTensor[2, 3, 110] = TsSym[$i, $j]@ Ttimes[C2SC2G[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 111] := BetaTensor[2, 3, 111] = TsSym[$i, $j]@ Ttimes[C2SS2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 112] := BetaTensor[2, 3, 112] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], C2FC2Gt[k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 113] := BetaTensor[2, 3, 113] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], C2FS2St[k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 114] := BetaTensor[2, 3, 114] = TsSym[$i, $j]@ Ttimes[C2SS2F[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 115] := BetaTensor[2, 3, 115] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], Y2Ft[k1, k2], C2FC2G[k2, $j]];
    BetaTensor[2, 3, 116] := BetaTensor[2, 3, 116] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], Y2Ft[k1, k2], C2FS2S[k2, $j]];
    BetaTensor[2, 3, 117] := BetaTensor[2, 3, 117] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], C2FS2Ft[k2, k3], Yuk[b1, k3, $j]];
    BetaTensor[2, 3, 118] := BetaTensor[2, 3, 118] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], Y2Ft[k1, k2], C2FS2F[k2, $j]];
    BetaTensor[2, 3, 119] := BetaTensor[2, 3, 119] = TsSym[$i, $j]@ Ttimes[C2SC2G[b1, b2], Y2S[b2, $a], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 120] := BetaTensor[2, 3, 120] = TsSym[$i, $j]@ Ttimes[C2SS2S[b1, b2], Y2S[b2, $a], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 121] := BetaTensor[2, 3, 121] = TsSym[$i, $j]@ Ttimes[C2SS2F[b1, b2], Y2S[b2, $a], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 122] := BetaTensor[2, 3, 122] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[TfermTil[A1, $i, k1], Yuk[b3, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 123] := BetaTensor[2, 3, 123] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Lam[b2, b3, b4, $a], Tdot[Yuk[b3, k1, $i], YukTil[b1, k1, k2], Yuk[b4, k2, $j]]];
    BetaTensor[2, 3, 124] := BetaTensor[2, 3, 124] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Lam[b2, b3, b4, $a], Tdot[Yuk[b1, $i, k1], YukTil[b4, k1, k2], Yuk[b3, $j, k2]]];
    BetaTensor[2, 3, 125] := BetaTensor[2, 3, 125] = TsSym[$i, $j]@ Ttimes[C2S[$a, b1], Lam[b1, b2, b3, b4], Tdot[Yuk[b3, k1, $i], YukTil[b2, k1, k2], Yuk[b4, k2, $j]]];
    BetaTensor[2, 3, 126] := BetaTensor[2, 3, 126] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b2, k1, $i], C2F[k1, k2], YukTil[b1, k2, k3], Yuk[b3, k3, $j]]];
    BetaTensor[2, 3, 127] := BetaTensor[2, 3, 127] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[b3, k2, k3], Yuk[b2, $j, k3]]];
    BetaTensor[2, 3, 128] := BetaTensor[2, 3, 128] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Lam[b2, b4, b5, b3], Tdot[Yuk[b4, k1, $i], YukTil[b1, k1, k2], Yuk[b5, k2, $j]]];
    BetaTensor[2, 3, 129] := BetaTensor[2, 3, 129] = TsSym[$i, $j]@ Ttimes[Lam2[b1, b2], Tdot[Yuk[b2, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k2, $j]]];
    BetaTensor[2, 3, 130] := BetaTensor[2, 3, 130] = TsSym[$i, $j]@ Ttimes[Lam2[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]]];
    BetaTensor[2, 3, 131] := BetaTensor[2, 3, 131] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b3, b5, b4, $a], Tdot[Yuk[b5, k1, $i], YukTil[b1, k1, k2], Yuk[b2, k2, $j]]];
    BetaTensor[2, 3, 132] := BetaTensor[2, 3, 132] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Lam[b2, b4, b5, b3], Y2S[b1, b5], Yuk[b4, $i, $j]];
    BetaTensor[2, 3, 133] := BetaTensor[2, 3, 133] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], YukTil[b2, k3, k2], TfG2t[A1, k3, k4], Yuk[$a, k4, k5], YukTil[b1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 134] := BetaTensor[2, 3, 134] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b3, $i, k1], YukTil[b2, k2, k1], Yuk[$a, k3, k2], YukTil[b3, k3, k4], Yuk[b1, k4, $j]]];
    BetaTensor[2, 3, 135] := BetaTensor[2, 3, 135] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Yuk[b2, k1, $i], YukTil[b3, k1, k2], Yuk[b1, k2, k3], YukTil[b2, k4, k3], Yuk[b3, k4, $j]]];
    BetaTensor[2, 3, 136] := BetaTensor[2, 3, 136] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[b2, k2, k3], Yuk[$a, k3, k4], YukTil[b1, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 137] := BetaTensor[2, 3, 137] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], C2Ft[k2, k3], Yuk[$a, k3, k4], YukTil[b1, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 138] := BetaTensor[2, 3, 138] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[b2, k2, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 139] := BetaTensor[2, 3, 139] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], TfG2[A1, k1, k2], YukTil[b1, k2, k3], Yuk[$a, k3, k4], YukTil[b2, k5, k4], TfermTil[A1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 140] := BetaTensor[2, 3, 140] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], Tferm[A1, k1, k2], YukTil[$a, k2, k3], Yuk[b2, k4, k3], TfG2[A1, k4, k5], YukTil[b2, k5, k6], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 141] := BetaTensor[2, 3, 141] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], Tferm[A1, k1, k2], YukTil[b1, k2, k3], Yuk[b2, k4, k3], TfG2[A1, k4, k5], YukTil[$a, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 142] := BetaTensor[2, 3, 142] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], TfermTil[A1, k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], TfG2t[A1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 143] := BetaTensor[2, 3, 143] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], YukTil[b2, k3, k2], TfG2t[A1, k3, k4], Yuk[b2, k4, k5], YukTil[$a, k5, k6], Yuk[b1, $j, k6]];
    BetaTensor[2, 3, 144] := BetaTensor[2, 3, 144] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b2, k4, k3], YukTil[b1, k5, k4], TfG2t[A1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 145] := BetaTensor[2, 3, 145] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b1, k4, k3], YukTil[b2, k5, k4], TfG2t[A1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 146] := BetaTensor[2, 3, 146] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k3, k4], Yuk[b2, k4, k5], YukTil[$a, k6, k5], Yuk[b2, $j, k6]];
    BetaTensor[2, 3, 147] := BetaTensor[2, 3, 147] = TsSym[$i, $j]@ Ttimes[Tdot[TfG2t[A1, $i, k1], Yuk[b1, k1, $j]], Tr[Tdot[YukTil[b2, k2, k3], TfermTil[A1, k2, k4], Yuk[b2, k4, k5], YukTil[b1, k5, k6], Yuk[$a, k3, k6]]]];
    BetaTensor[2, 3, 148] := BetaTensor[2, 3, 148] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[b1, k1, k2], YukTil[b2, k3, k2], TfG2t[A1, k3, k4], Yuk[$a, k4, k5], YukTil[b2, k5, k6], Yuk[b1, $j, k6]];
    BetaTensor[2, 3, 149] := BetaTensor[2, 3, 149] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], Yuk[b2, k4, k3], TfG2[A1, k4, k5], YukTil[b2, k5, k6], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 150] := BetaTensor[2, 3, 150] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b2, k4, k5], YukTil[b1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 151] := BetaTensor[2, 3, 151] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k3, k4], Yuk[b2, k4, k5], YukTil[b1, k6, k5], Yuk[b2, $j, k6]];
    BetaTensor[2, 3, 152] := BetaTensor[2, 3, 152] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k2, k3], Yuk[b2, k4, k3], YukTil[b1, k5, k4], TfG2t[A1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 153] := BetaTensor[2, 3, 153] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], Tferm[A1, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k6, k5], TfG2[A1, k6, $j]];
    BetaTensor[2, 3, 154] := BetaTensor[2, 3, 154] = TsSym[$i, $j]@ Ttimes[-1, Tdot[TfermTil[A1, $i, k1], Y2F[k1, k2], Yuk[$a, k2, k3], YukTil[b1, k4, k3], TfG2t[A1, k4, k5], Yuk[b1, k5, $j]]];
    BetaTensor[2, 3, 155] := BetaTensor[2, 3, 155] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k2, k1], Y2Ft[k2, k3], YukTil[b1, k3, k4], TfG2t[A1, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 156] := BetaTensor[2, 3, 156] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k4, k3], Y2F[k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 157] := BetaTensor[2, 3, 157] = TsSym[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Y2F[k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 158] := BetaTensor[2, 3, 158] = TsSym[$i, $j]@ Ttimes[-1, Tdot[TfermTil[A1, $i, k1], Y2F[k1, k2], Yuk[b1, k2, k3], TfG2[A1, k3, k4], YukTil[$a, k4, k5], Yuk[b1, k5, $j]]];
    BetaTensor[2, 3, 159] := BetaTensor[2, 3, 159] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], Y2Ft[k1, k2], YukTil[b1, k2, k3], TfermTil[A1, k3, k4], Yuk[b1, k5, k4], TfG2[A1, k5, $j]];
    BetaTensor[2, 3, 160] := BetaTensor[2, 3, 160] = TsSym[$i, $j]@ Ttimes[-1, Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], TfermTil[A1, k2, k3], Yuk[b1, k4, k3], Y2Ft[k4, k5], TfG2[A1, k5, $j]]];
    BetaTensor[2, 3, 161] := BetaTensor[2, 3, 161] = TsSym[$i, $j]@ Tdot[Y2F[$i, k1], TfermTil[A1, k1, k2], Yuk[b1, k2, k3], YukTil[$a, k3, k4], Yuk[b1, k5, k4], TfG2[A1, k5, $j]];
    BetaTensor[2, 3, 162] := BetaTensor[2, 3, 162] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], Tferm[A1, k1, k2], YukTil[b1, k2, k3], Yuk[b2, k3, $j]], TsG2[A1, b2, b3], Y2S[b3, $a]];
    BetaTensor[2, 3, 163] := BetaTensor[2, 3, 163] = TsSym[$i, $j]@ Ttimes[Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], TfG2[A1, k2, k3], YukTil[b1, k3, k4], Yuk[b2, k4, $j]], Y2S[b2, b1]];
    BetaTensor[2, 3, 164] := BetaTensor[2, 3, 164] = TsSym[$i, $j]@ Ttimes[Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2], YukTil[b1, k3, k2], TfG2t[A1, k3, k4], Yuk[b2, k4, $j]], Y2S[b1, b2]];
    BetaTensor[2, 3, 165] := BetaTensor[2, 3, 165] = TsSym[$i, $j]@ Tdot[Y2F[$i, k1], TfermTil[A1, k1, k2], Y2F[k2, k3], Yuk[$a, k3, k4], TfG2[A1, k4, $j]];
    BetaTensor[2, 3, 166] := BetaTensor[2, 3, 166] = TsSym[$i, $j]@ Tdot[Y2F[$i, k1], TfermTil[A1, k1, k2], Yuk[$a, k2, k3], TfG2[A1, k3, k4], Y2Ft[k4, $j]];
    BetaTensor[2, 3, 167] := BetaTensor[2, 3, 167] = TsSym[$i, $j]@ Ttimes[-1, Tdot[TfermTil[A1, $i, k1], Y2F[k1, k2], Yuk[b1, k2, $j]], TsG2[A1, b1, b2], Y2S[b2, $a]];
    BetaTensor[2, 3, 168] := BetaTensor[2, 3, 168] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b3, $i, k1], YukTil[b1, k1, k2], Yuk[$a, k3, k2], YukTil[b2, k3, k4], Yuk[b3, $j, k4]]];
    BetaTensor[2, 3, 169] := BetaTensor[2, 3, 169] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b3, $i, k1], YukTil[b1, k2, k1], Yuk[b3, k3, k2], YukTil[$a, k3, k4], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 170] := BetaTensor[2, 3, 170] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b3, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k3, k2], YukTil[b3, k4, k3], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 171] := BetaTensor[2, 3, 171] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b2, k1, $i], YukTil[b3, k2, k1], Yuk[$a, k3, k2], YukTil[b3, k3, k4], Yuk[b1, k4, $j]]];
    BetaTensor[2, 3, 172] := BetaTensor[2, 3, 172] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b3, k3, k2], YukTil[b2, k3, k4], Yuk[b3, $j, k4]]];
    BetaTensor[2, 3, 173] := BetaTensor[2, 3, 173] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b3, k1, k2], Yuk[b1, k3, k2], YukTil[b3, k4, k3], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 174] := BetaTensor[2, 3, 174] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b1, k2, k1], Yuk[b3, k3, k2], YukTil[b2, k3, k4], Yuk[b3, k4, $j]]];
    BetaTensor[2, 3, 175] := BetaTensor[2, 3, 175] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b3, k1, k2], Yuk[b1, k3, k2], YukTil[b3, k4, k3], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 176] := BetaTensor[2, 3, 176] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b3, k3, k2], YukTil[b2, k4, k3], Yuk[$a, k4, k1]]], Yuk[b3, $i, $j]];
    BetaTensor[2, 3, 177] := BetaTensor[2, 3, 177] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k4, k2]]], Yuk[b3, $i, $j]];
    BetaTensor[2, 3, 178] := BetaTensor[2, 3, 178] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 179] := BetaTensor[2, 3, 179] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[b2, k2, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 180] := BetaTensor[2, 3, 180] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[b2, k2, k3], Yuk[$a, k3, k4], YukTil[b2, k5, k4], Yuk[b1, $j, k5]];
    BetaTensor[2, 3, 181] := BetaTensor[2, 3, 181] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[$a, k2, k1], C2Ft[k2, k3], Yuk[b2, k3, k4], YukTil[b1, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 182] := BetaTensor[2, 3, 182] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[$a, k3, k2], C2F[k3, k4], YukTil[b2, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 183] := BetaTensor[2, 3, 183] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, k4], YukTil[b1, k5, k4], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 184] := BetaTensor[2, 3, 184] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b1, k4, k3], C2Ft[k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 185] := BetaTensor[2, 3, 185] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], C2F[k3, k4], YukTil[b1, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 186] := BetaTensor[2, 3, 186] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], C2Ft[k2, k3], Yuk[b2, k3, k4], YukTil[b1, k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 187] := BetaTensor[2, 3, 187] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], YukTil[$a, k4, k5], Yuk[b2, k2, k5]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 188] := BetaTensor[2, 3, 188] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[$a, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k2, k5]]], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 189] := BetaTensor[2, 3, 189] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, k4], YukTil[b1, k4, k5], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 190] := BetaTensor[2, 3, 190] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[b2, k2, k3], Yuk[b1, k4, k3], YukTil[$a, k5, k4], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 191] := BetaTensor[2, 3, 191] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[b2, k2, k3], Yuk[$a, k4, k3], YukTil[b2, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 192] := BetaTensor[2, 3, 192] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[b2, k2, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], Yuk[$a, $j, k5]];
    BetaTensor[2, 3, 193] := BetaTensor[2, 3, 193] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[$a, k1, k2], Yuk[b2, k2, $j]], Tr[Tdot[YukTil[b1, k3, k4], C2Ft[k3, k5], Yuk[b2, k5, k4]]]];
    BetaTensor[2, 3, 194] := BetaTensor[2, 3, 194] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]], Tr[Tdot[YukTil[b2, k3, k4], C2Ft[k3, k5], Yuk[b1, k5, k4]]]];
    BetaTensor[2, 3, 195] := BetaTensor[2, 3, 195] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b2, k1, $i], YukTil[$a, k2, k1], Y2F[k2, k3], Yuk[b1, k3, $j]]];
    BetaTensor[2, 3, 196] := BetaTensor[2, 3, 196] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], Y2F[k2, k3], Yuk[b2, k3, $j]]];
    BetaTensor[2, 3, 197] := BetaTensor[2, 3, 197] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b1, k2, k1], Y2F[k2, k3], Yuk[b2, k3, $j]]];
    BetaTensor[2, 3, 198] := BetaTensor[2, 3, 198] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k2]]], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 199] := BetaTensor[2, 3, 199] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b3, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[b3, k4, $j]]];
    BetaTensor[2, 3, 200] := BetaTensor[2, 3, 200] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b3, k1, k2], Yuk[b2, k2, k3], YukTil[b1, k4, k3], Yuk[b3, $j, k4]]];
    BetaTensor[2, 3, 201] := BetaTensor[2, 3, 201] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k4, k2]]], Yuk[b3, $i, $j]];
    BetaTensor[2, 3, 202] := BetaTensor[2, 3, 202] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], Y2Ft[k2, k3], YukTil[$a, k3, k4], Yuk[b1, $j, k4]];
    BetaTensor[2, 3, 203] := BetaTensor[2, 3, 203] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[$a, k2, k3], Y2F[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 204] := BetaTensor[2, 3, 204] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], C2Ft[k2, k3], Y2F[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 205] := BetaTensor[2, 3, 205] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[Y2Ft[k1, k2], C2F[k2, k3], YukTil[$a, k3, k4], Yuk[b1, k1, k4]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 206] := BetaTensor[2, 3, 206] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], Y2Ft[k4, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 207] := BetaTensor[2, 3, 207] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Yuk[b2, k3, k2], C2F[k3, k4], YukTil[b2, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 208] := BetaTensor[2, 3, 208] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], C2F[k3, k4], YukTil[b2, k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 209] := BetaTensor[2, 3, 209] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[$a, k2, k5]]], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 210] := BetaTensor[2, 3, 210] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], Y2Ft[k2, k3], YukTil[$a, k3, k4], Yuk[b1, $j, k4]];
    BetaTensor[2, 3, 211] := BetaTensor[2, 3, 211] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k3, k2], Y2F[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 212] := BetaTensor[2, 3, 212] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Y2F[k1, k2], Yuk[b1, k2, k3], YukTil[$a, k4, k3], Yuk[b1, $j, k4]];
    BetaTensor[2, 3, 213] := BetaTensor[2, 3, 213] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], Y2Ft[k2, k3], YukTil[b1, k3, k4], Yuk[$a, $j, k4]];
    BetaTensor[2, 3, 214] := BetaTensor[2, 3, 214] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[b3, k1, $i], YukTil[$a, k1, k2], Yuk[b1, k2, $j]], Y2S[b3, b2]];
    BetaTensor[2, 3, 215] := BetaTensor[2, 3, 215] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b3, $j, k2]], Y2S[b3, b2]];
    BetaTensor[2, 3, 216] := BetaTensor[2, 3, 216] = TsSym[$i, $j]@ Ttimes[C2S[b1, $a], Tdot[Yuk[b2, $i, k1], YukTil[b1, k1, k2], Yuk[b3, $j, k2]], Y2S[b3, b2]];
    BetaTensor[2, 3, 217] := BetaTensor[2, 3, 217] = TsSym[$i, $j]@ Ttimes[C2S[$a, b1], Tdot[Yuk[b2, $i, k1], YukTil[b3, k1, k2], Yuk[b2, $j, k2]], Y2S[b3, b1]];
    BetaTensor[2, 3, 218] := BetaTensor[2, 3, 218] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], C2F[k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, $j]], Y2S[b2, b1]];
    BetaTensor[2, 3, 219] := BetaTensor[2, 3, 219] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], C2Ft[k2, k3], Yuk[b2, k3, $j]], Y2S[b1, b2]];
    BetaTensor[2, 3, 220] := BetaTensor[2, 3, 220] = TsSym[$i, $j]@ Ttimes[Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3], Yuk[b2, k3, $j]], Y2S[b2, b1]];
    BetaTensor[2, 3, 221] := BetaTensor[2, 3, 221] = TsSym[$i, $j]@ Ttimes[Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[b2, k3, k2], Yuk[$a, $j, k3]], Y2S[b2, b1]];
    BetaTensor[2, 3, 222] := BetaTensor[2, 3, 222] = TsSym[$i, $j]@ Tdot[C2Ft[$i, k1], Y2F[k1, k2], Y2F[k2, k3], Yuk[$a, $j, k3]];
    BetaTensor[2, 3, 223] := BetaTensor[2, 3, 223] = TsSym[$i, $j]@ Ttimes[C2S[b1, b2], Y2S[b3, b1], Y2S[b3, $a], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 224] := BetaTensor[2, 3, 224] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b1, k1, $i], YukTil[b4, k1, k2], Yuk[b3, k3, k2], YukTil[b4, k4, k3], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 225] := BetaTensor[2, 3, 225] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b4, $i, k1], YukTil[b3, k2, k1], Yuk[b4, k3, k2], YukTil[b2, k3, k4], Yuk[b1, $j, k4]]];
    BetaTensor[2, 3, 226] := BetaTensor[2, 3, 226] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b4, $i, k1], YukTil[b1, k1, k2], Yuk[b3, k3, k2], YukTil[b4, k4, k3], Yuk[b2, k4, $j]]];
    BetaTensor[2, 3, 227] := BetaTensor[2, 3, 227] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b4, $i, k1], YukTil[b3, k1, k2], Yuk[b1, k3, k2], YukTil[b2, k4, k3], Yuk[b4, $j, k4]]];
    BetaTensor[2, 3, 228] := BetaTensor[2, 3, 228] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, b4], Tdot[Yuk[b1, $i, k1], YukTil[b4, k2, k1], Yuk[$a, k3, k2], YukTil[b3, k3, k4], Yuk[b2, $j, k4]]];
    BetaTensor[2, 3, 229] := BetaTensor[2, 3, 229] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, b4], Tdot[Yuk[b2, k1, $i], YukTil[b1, k1, k2], Yuk[b4, k3, k2], YukTil[$a, k4, k3], Yuk[b3, k4, $j]]];
    BetaTensor[2, 3, 230] := BetaTensor[2, 3, 230] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, b4], Tdot[Yuk[$a, $i, k1], YukTil[b4, k1, k2], Yuk[b1, k3, k2], YukTil[b3, k3, k4], Yuk[b2, $j, k4]]];
    BetaTensor[2, 3, 231] := BetaTensor[2, 3, 231] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b4, k3, k1], YukTil[b2, k3, k4], Yuk[b1, k4, k2]]], Yuk[b4, $i, $j]];
    BetaTensor[2, 3, 232] := BetaTensor[2, 3, 232] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, b4], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[$a, k4, k1]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 233] := BetaTensor[2, 3, 233] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b1, $i, k1], YukTil[b3, k2, k1], Y2F[k2, k3], Yuk[b2, k3, $j]]];
    BetaTensor[2, 3, 234] := BetaTensor[2, 3, 234] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b2, k1, $i], YukTil[b4, k1, k2], Yuk[b3, k2, $j]], Y2S[b4, b1]];
    BetaTensor[2, 3, 235] := BetaTensor[2, 3, 235] = TsSym[$i, $j]@ Ttimes[Lam[b1, b2, b3, $a], Tdot[Yuk[b4, $i, k1], YukTil[b3, k1, k2], Yuk[b2, $j, k2]], Y2S[b4, b1]];
    BetaTensor[2, 3, 236] := BetaTensor[2, 3, 236] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[b3, k2, k3], YukTil[$a, k4, k3], Yuk[b3, k5, k4], YukTil[b1, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 237] := BetaTensor[2, 3, 237] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[b3, k2, k3], YukTil[$a, k4, k3], Yuk[b1, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 238] := BetaTensor[2, 3, 238] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[$a, k2, k3], YukTil[b3, k4, k3], Yuk[b1, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 239] := BetaTensor[2, 3, 239] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[b3, k2, $j]], Tr[Tdot[YukTil[b2, k3, k4], Yuk[b3, k5, k3], YukTil[b1, k5, k6], Yuk[$a, k6, k4]]]];
    BetaTensor[2, 3, 240] := BetaTensor[2, 3, 240] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[b1, k3, k2], YukTil[b3, k3, k4], Yuk[$a, k5, k4], YukTil[b2, k5, k6], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 241] := BetaTensor[2, 3, 241] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[b3, k3, k2], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[b2, k5, k6], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 242] := BetaTensor[2, 3, 242] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[b3, k3, k2], YukTil[b1, k3, k4], Yuk[$a, k5, k4], YukTil[b2, k5, k6], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 243] := BetaTensor[2, 3, 243] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[b1, k5, k4], YukTil[b2, k5, k6], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 244] := BetaTensor[2, 3, 244] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[b3, k3, k2], YukTil[$a, k4, k3], Yuk[b2, k4, k5], YukTil[b3, k5, k6], Yuk[b1, $j, k6]];
    BetaTensor[2, 3, 245] := BetaTensor[2, 3, 245] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[b1, k4, k5], YukTil[b2, k5, k6], Yuk[b3, $j, k6]];
    BetaTensor[2, 3, 246] := BetaTensor[2, 3, 246] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[$a, k4, k5], YukTil[b2, k6, k5], Yuk[b3, k6, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 247] := BetaTensor[2, 3, 247] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[$a, k3, k2], Y2Ft[k3, k4], YukTil[b1, k4, k5], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 248] := BetaTensor[2, 3, 248] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[$a, k2, k3], YukTil[b1, k4, k3], Y2F[k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 249] := BetaTensor[2, 3, 249] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[$a, k3, k2], YukTil[b1, k3, k4], Yuk[b3, k4, $j]], Y2S[b3, b2]];
    BetaTensor[2, 3, 250] := BetaTensor[2, 3, 250] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[b3, k2, k3], YukTil[b1, k4, k3], Yuk[$a, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 251] := BetaTensor[2, 3, 251] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[b1, k2, k3], YukTil[b3, k3, k4], Yuk[$a, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 252] := BetaTensor[2, 3, 252] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[b3, k2, $j]], Tr[Tdot[YukTil[b2, k3, k4], Yuk[b3, k5, k3], YukTil[$a, k5, k6], Yuk[b1, k6, k4]]]];
    BetaTensor[2, 3, 253] := BetaTensor[2, 3, 253] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[b1, k2, k3], YukTil[$a, k4, k3], Yuk[b3, k4, k5], YukTil[b2, k5, k6], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 254] := BetaTensor[2, 3, 254] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[b3, k2, k3], YukTil[b1, k4, k3], Yuk[b3, k4, k5], YukTil[$a, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 255] := BetaTensor[2, 3, 255] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Yuk[b2, k2, k3], YukTil[b3, k3, k4], Yuk[b1, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, $j]];
    BetaTensor[2, 3, 256] := BetaTensor[2, 3, 256] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Yuk[b2, k2, k3], YukTil[b1, k4, k3], Yuk[b3, k4, k5], YukTil[b2, k5, k6], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 257] := BetaTensor[2, 3, 257] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[$a, k1, k2], Yuk[b2, k2, $j]], Tr[Tdot[YukTil[b3, k3, k4], Yuk[b2, k5, k3], YukTil[b3, k5, k6], Yuk[b1, k6, k4]]]];
    BetaTensor[2, 3, 258] := BetaTensor[2, 3, 258] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[b3, k3, k2], YukTil[b2, k4, k3], Yuk[$a, k4, k5], YukTil[b3, k6, k5], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 259] := BetaTensor[2, 3, 259] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[b3, k3, k2], YukTil[$a, k3, k4], Yuk[b3, k5, k4], YukTil[b2, k6, k5], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 260] := BetaTensor[2, 3, 260] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Yuk[b2, k2, k3], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[b3, k6, k5], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 261] := BetaTensor[2, 3, 261] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[b3, k6, k5], Yuk[b1, k6, $j]];
    BetaTensor[2, 3, 262] := BetaTensor[2, 3, 262] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]], Tr[Tdot[YukTil[b3, k3, k4], Yuk[b1, k5, k3], YukTil[b3, k5, k6], Yuk[b2, k6, k4]]]];
    BetaTensor[2, 3, 263] := BetaTensor[2, 3, 263] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b1, k4, k3], Yuk[b3, k4, k5], YukTil[b2, k6, k5], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 264] := BetaTensor[2, 3, 264] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[b1, k5, k4], YukTil[b3, k5, k6], Yuk[b2, $j, k6]];
    BetaTensor[2, 3, 265] := BetaTensor[2, 3, 265] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[b1, k6, k5], Yuk[b3, k6, $j]];
    BetaTensor[2, 3, 266] := BetaTensor[2, 3, 266] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[b3, k6, k5], Yuk[$a, k6, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 267] := BetaTensor[2, 3, 267] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[$a, k6, k5], Yuk[b3, k6, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 268] := BetaTensor[2, 3, 268] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[$a, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 269] := BetaTensor[2, 3, 269] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Y2F[k2, k3], Yuk[b1, k3, k4], YukTil[$a, k5, k4], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 270] := BetaTensor[2, 3, 270] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, $i, k1], YukTil[$a, k1, k2], Yuk[b2, $j, k2]], Tr[Tdot[YukTil[b1, k3, k4], Y2F[k3, k5], Yuk[b2, k5, k4]]]];
    BetaTensor[2, 3, 271] := BetaTensor[2, 3, 271] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Y2F[k2, k3], Yuk[b2, k3, k4], YukTil[$a, k5, k4], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 272] := BetaTensor[2, 3, 272] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[b1, k3, k2], Y2Ft[k3, k4], YukTil[$a, k4, k5], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 273] := BetaTensor[2, 3, 273] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Y2F[k2, k3], Yuk[$a, k3, k4], YukTil[b2, k5, k4], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 274] := BetaTensor[2, 3, 274] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Yuk[b2, k2, k3], YukTil[b1, k4, k3], Y2F[k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 275] := BetaTensor[2, 3, 275] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k1, k2], Yuk[$a, k2, k3], YukTil[b2, k4, k3], Y2F[k4, k5], Yuk[b1, k5, $j]];
    BetaTensor[2, 3, 276] := BetaTensor[2, 3, 276] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[b1, k2, k3], YukTil[$a, k4, k3], Y2F[k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 277] := BetaTensor[2, 3, 277] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], Y2F[k2, k3], Yuk[b2, k3, k4], YukTil[b1, k5, k4], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 278] := BetaTensor[2, 3, 278] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], Y2Ft[k3, k4], YukTil[b1, k4, k5], Yuk[b2, $j, k5]];
    BetaTensor[2, 3, 279] := BetaTensor[2, 3, 279] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]], Tr[Tdot[YukTil[b2, k3, k4], Y2F[k3, k5], Yuk[b1, k5, k4]]]];
    BetaTensor[2, 3, 280] := BetaTensor[2, 3, 280] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], Y2Ft[k3, k4], YukTil[b2, k4, k5], Yuk[b1, $j, k5]];
    BetaTensor[2, 3, 281] := BetaTensor[2, 3, 281] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b1, k4, k3], Y2F[k4, k5], Yuk[b2, k5, $j]];
    BetaTensor[2, 3, 282] := BetaTensor[2, 3, 282] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], YukTil[$a, k5, k4], Yuk[b2, k5, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 283] := BetaTensor[2, 3, 283] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k4], YukTil[b2, k5, k4], Yuk[$a, k5, k2]]], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 284] := BetaTensor[2, 3, 284] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], Yuk[b2, k5, k2]]], Yuk[b2, $i, $j]];
    BetaTensor[2, 3, 285] := BetaTensor[2, 3, 285] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[$a, k3, k2], YukTil[b3, k3, k4], Yuk[b1, $j, k4]], Y2S[b2, b3]];
    BetaTensor[2, 3, 286] := BetaTensor[2, 3, 286] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[b3, k3, k2], YukTil[$a, k4, k3], Yuk[b1, k4, $j]], Y2S[b2, b3]];
    BetaTensor[2, 3, 287] := BetaTensor[2, 3, 287] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, $i, k1], YukTil[b2, k2, k1], Yuk[b1, k3, k2], YukTil[$a, k3, k4], Yuk[b3, k4, $j]], Y2S[b2, b3]];
    BetaTensor[2, 3, 288] := BetaTensor[2, 3, 288] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[$a, k1, k2], Yuk[b2, k3, k2], YukTil[b1, k4, k3], Yuk[b3, k4, $j]], Y2S[b2, b3]];
    BetaTensor[2, 3, 289] := BetaTensor[2, 3, 289] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[b2, k2, k1], Yuk[$a, k3, k2], YukTil[b2, k3, k4], Yuk[b3, k4, $j]], Y2S[b3, b1]];
    BetaTensor[2, 3, 290] := BetaTensor[2, 3, 290] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k3, k4], Yuk[b2, $j, k4]], Y2S[b1, b3]];
    BetaTensor[2, 3, 291] := BetaTensor[2, 3, 291] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[b1, $j, k4]], Y2S[b3, b2]];
    BetaTensor[2, 3, 292] := BetaTensor[2, 3, 292] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b1, k3, k4], Yuk[b3, k4, $j]], Y2S[b3, b2]];
    BetaTensor[2, 3, 293] := BetaTensor[2, 3, 293] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k2]]], Y2S[b3, b2], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 294] := BetaTensor[2, 3, 294] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[$a, k4, k1]]], Y2S[b1, b2], Yuk[b3, $i, $j]];
    BetaTensor[2, 3, 295] := BetaTensor[2, 3, 295] = TsSym[$i, $j]@ Tdot[Yuk[b1, k1, $i], Y2Ft[k1, k2], YukTil[$a, k2, k3], Y2F[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 296] := BetaTensor[2, 3, 296] = TsSym[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Y2F[k2, k3], Y2F[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 297] := BetaTensor[2, 3, 297] = TsSym[$i, $j]@ Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], Y2F[k2, k3], Y2F[k3, k4], Yuk[b1, k4, $j]];
    BetaTensor[2, 3, 298] := BetaTensor[2, 3, 298] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k4, k3], Y2Ft[k4, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 299] := BetaTensor[2, 3, 299] = TsSym[$i, $j]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Y2F[k3, k4], Yuk[$a, k4, k2]]], Yuk[b1, $i, $j]];
    BetaTensor[2, 3, 300] := BetaTensor[2, 3, 300] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, $i, k1], YukTil[$a, k2, k1], Y2F[k2, k3], Yuk[b2, k3, $j]], Y2S[b2, b1]];
    BetaTensor[2, 3, 301] := BetaTensor[2, 3, 301] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k2, k1], Y2F[k2, k3], Yuk[b2, k3, $j]], Y2S[b1, b2]];
    BetaTensor[2, 3, 302] := BetaTensor[2, 3, 302] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[b1, k1, $i], YukTil[$a, k1, k2], Yuk[b2, k2, $j]], Y2S[b1, b3], Y2S[b2, b3]];
    BetaTensor[2, 3, 303] := BetaTensor[2, 3, 303] = TsSym[$i, $j]@ Ttimes[Tdot[Yuk[$a, $i, k1], YukTil[b1, k1, k2], Yuk[b2, $j, k2]], Y2S[b1, b3], Y2S[b2, b3]];
(* gamma_5 terms *)
    BetaTensor[2, 3, 304] := BetaTensor[2, 3, 304] = TsSym[$i, $j]@ Ttimes[Tdot[sig3, Yuk[b1, $i, k1], TfG2[A1, k1, k2], TfG2[A2, k2, $j]], Tr@Tdot[sig3, Yuk[b1, k3, k4], Tferm[A1, k4, k5], Tferm[A2, k5, k6], YukTil[$a, k6, k3]] ];
    BetaTensor[2, 3, 305] := BetaTensor[2, 3, 305] = TsSym[$i, $j]@ Ttimes[Tdot[TfG2t[A1, $i, k1], sig3, Yuk[b1, k1, k2], TfG2[A2, k2, $j]], Tr@Tdot[sig3, Yuk[b1, k3, k4], Tferm[A2, k4, k5], YukTil[$a, k5, k6], TfermTil[A1, k6, k3]] ];
    BetaTensor[2, 3, 306] := BetaTensor[2, 3, 306] = TsSym[$i, $j]@ Ttimes[Tdot[sig3, Yuk[b1, $i, k1], TfG2[A2, k1, k2], TfG2[A1, k2, $j]], Tr@Tdot[sig3, Yuk[b1, k3, k4], Tferm[A1, k4, k5], Tferm[A2, k5, k6], YukTil[$a, k6, k3]] ];
    BetaTensor[2, 3, 307] := BetaTensor[2, 3, 307] = TsSym[$i, $j]@ Ttimes[Tdot[sig3, Yuk[b1, $i, k1], TfG2[A1, k1, k2], TfG2[A2, k2, $j]], Tr@Tdot[sig3, Yuk[b1, k3, k4], Tferm[A2, k4, k5], YukTil[$a, k5, k6], TfermTil[A1, k6, k3]] ];
    BetaTensor[2, 3, 308] := BetaTensor[2, 3, 308] = TsSym[$i, $j]@ Ttimes[Tdot[TfG2t[A1, $i, k1], sig3, Yuk[b1, k1, k2], TfG2[A2, k2, $j]], Tr@Tdot[sig3, Yuk[b1, k3, k4], Tferm[A1, k4, k5], Tferm[A2, k5, k6], YukTil[$a, k6, k3]] ];

(*#####################################*)
(*-----------Quartic Tensors-----------*)
(*#####################################*)
    BetaTensor[3, 0, 1] := BetaTensor[3, 0, 1] = - Global`\[Epsilon] Lam[$a, $b, $c, $d];

    BetaTensor[3, 1, 1] := BetaTensor[3, 1, 1] = FourGLam[$a, $b, $c, $d] // TsSym[$b, $c, $d];
    BetaTensor[3, 1, 2] := BetaTensor[3, 1, 2] = C2S[$a, b1] Lam[b1, $b, $c, $d] // Expand // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 1, 3] := BetaTensor[3, 1, 3] = Lam[$a, $b, b1, b2] Lam[b1, b2, $c, $d] // Expand // TsSym[$b, $c, $d];
    BetaTensor[3, 1, 4] := BetaTensor[3, 1, 4] = Y2S[$a, b1] Lam[b1, $b, $c, $d] // Expand // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 1, 5] := BetaTensor[3, 1, 5] = Tr @ Tdot[Yuk[$a, k4, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4]] // TsSym[$b, $c, $d];

(* y^0, lam^0 terms*)
    BetaTensor[3, 2, 1] := BetaTensor[3, 2, 1] = Ttimes[Tscal[A1, $a, b1], FourGLam[b1, b2, $b, $c], TsG2[A1, b2, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 2] := BetaTensor[3, 2, 2] = Ttimes[C2S[$a, b1], FourGLam[b1, $b, $c, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 3] := BetaTensor[3, 2, 3] = Ttimes[TsG2[A1, $a, b1], TsG2[A2, b1, $b], C2G[A2, A3], Tscal[A1, $c, b2], TsG2[A3, b2, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 4] := BetaTensor[3, 2, 4] = Ttimes[TsG2[A1, $a, b1], TsG2[A2, b1, $b], S2S[A2, A3], Tscal[A1, $c, b2], TsG2[A3, b2, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 5] := BetaTensor[3, 2, 5] = Ttimes[TsG2[A1, $a, b1], TsG2[A2, b1, $b], S2F[A2, A3], Tscal[A1, $c, b2], TsG2[A3, b2, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 6] := BetaTensor[3, 2, 6] = Ttimes[FourGLam[$a, b1, $d, b2], Lam[b1, b2, $b, $c]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 7] := BetaTensor[3, 2, 7] = Ttimes[FourGLam[$a, $b, b1, b2], Lam[b1, b2, $c, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 8] := BetaTensor[3, 2, 8] = 0 (*Ttimes[Lam[$a, $b, b1, b2], C2S[b1, $c], C2S[b2, $d]] // TsSym[$a, $b, $c, $d]*);
    BetaTensor[3, 2, 9] := BetaTensor[3, 2, 9] = Ttimes[C2S[$a, b1], C2S[b1, b2], Lam[b2, $b, $c, $d] ]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 2, 10] := BetaTensor[3, 2, 10] = Ttimes[C2SC2G[$a, b1], Lam[b1, $b, $c, $d] ]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 2, 11] := BetaTensor[3, 2, 11] = Ttimes[C2SS2S[$a, b1], Lam[b1, $b, $c, $d] ]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 2, 12] := BetaTensor[3, 2, 12] = Ttimes[C2SS2F[$a, b1], Lam[b1, $b, $c, $d] ]  // TsSym4[$a, $b, $c, $d];
(* y^0, lam^{n>0} terms*)
    BetaTensor[3, 2, 13] := BetaTensor[3, 2, 13] = Ttimes[Tscal[A1, $a, b1], Lam[b1, b2, b3, b4], TsG2[A1, $b, b2], Lam[b3, b4, $c, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 14] := BetaTensor[3, 2, 14] = Ttimes[Lam[$a, $b, b1, b2], C2S[b2, b3], Lam[b3, b1, $c, $d]] // TsSym[$b, $c, $d];
    BetaTensor[3, 2, 15] := BetaTensor[3, 2, 15] = Ttimes[C2S[$a, b1], Lam[b1, $b, b2, b3], Lam[b2, b3, $c, $d]]  // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 16] := BetaTensor[3, 2, 16] = Ttimes[Lam2[$a, b1], Lam[b1,$b, $c, $d]]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 2, 17] := BetaTensor[3, 2, 17] = Ttimes[Lam[$a, b1, b2, b3], Lam[$b, b4, b2, b3], Lam[b1, b4, $c, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 18] := BetaTensor[3, 2, 18] = 0 (*Ttimes[Lam[$a, $b, b1, b2], Lam[b1, b2, b3, b4], Lam[b3, b4, $c, $d]] // TsSym[$b, $c, $d]*);
(*y^2 terms*)
    BetaTensor[3, 2, 19] := BetaTensor[3, 2, 19] = Ttimes[Tr @ Tdot[Yuk[$a, k4, k1], Tferm[A1, k1, k2], Tferm[A2, k2, k3], YukTil[$b, k3, k4]], TsG2[A1, $c, b1], TsG2[A2, b1, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 20] := BetaTensor[3, 2, 20] = Ttimes[Y2S[$a, b1], FourGLam[b1, $b, $c, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 21] := BetaTensor[3, 2, 21] = Ttimes[Y2SC2F[$a, b1], Lam[b1, $b, $c, $d]] // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 2, 22] := BetaTensor[3, 2, 22] = 0 (*Ttimes[C2S[$a, b1], Y2S[b1, b2], Lam[b2, $b, $c, $d] ]  // TsSym4[$a, $b, $c, $d]*);
    BetaTensor[3, 2, 23] := BetaTensor[3, 2, 23] = Ttimes[Lam[$a, $b, b1, b2], Y2S[b2, b3], Lam[b3, b1, $c, $d]] // TsSym[$b, $c, $d];
(*y^4 terms*)
    BetaTensor[3, 2, 24] := BetaTensor[3, 2, 24] = 0 (*Tr @ Tdot[Yuk[$a, k6, k1], Tferm[A1, k1, k2], YukTil[$b, k2, k3], Yuk[$c, k3, k4], TfG2[A1, k4, k5], YukTil[$d, k5, k6]] // TsSym[$b, $c, $d]*);
    BetaTensor[3, 2, 25] := BetaTensor[3, 2, 25] = Tr @ Tdot[Yuk[b1, k4, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4]] C2S[$a, b1] // Expand // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 26] := BetaTensor[3, 2, 26] = Tr @ Tdot[Yuk[$a, k5, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4], C2Ft[k4, k5]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 27] := BetaTensor[3, 2, 27] = Tr @ Tdot[Yuk[$a, k4, k1], YukTil[b1, k1, k2], Yuk[$b, k2, k3], YukTil[b2, k3, k4]] Lam[b1, b2, $c, $d] //Expand // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 28] := BetaTensor[3, 2, 28] = 0 (*Tr @ Tdot[Yuk[$a, k4, k1], YukTil[$b, k1, k2], Yuk[b1, k2, k3], YukTil[b2, k3, k4]] Lam[b1, b2, $c, $d] // TsSym[$a, $b, $c, $d]*);
    BetaTensor[3, 2, 29] := BetaTensor[3, 2, 29] = Y4cS[$a, b1] Lam[b1, $b, $c, $d] // TsSym4[$a, $b, $c, $d];
    BetaTensor[3, 2, 30] := BetaTensor[3, 2, 30] = Y2SY2F[$a, b1] Lam[b1, $b, $c, $d] // TsSym4[$a, $b, $c, $d];
(*y^6 terms*)
    BetaTensor[3, 2, 31] := BetaTensor[3, 2, 31] = Tr @ Tdot[Yuk[b1, k6, k1], YukTil[$a, k1, k2], Yuk[b1, k2, k3], YukTil[$b, k3, k4], Yuk[$c, k4, k5], YukTil[$d, k5, k6]] // TsSym[$a, $b, $c, $d];
    BetaTensor[3, 2, 32] := BetaTensor[3, 2, 32] = Tr @ Tdot[Yuk[$a, k6, k1], YukTil[b1, k1, k2], Yuk[$b, k2, k3], YukTil[$c, k3, k4], Yuk[b1, k4, k5], YukTil[$d, k5, k6]] // TsSym[$b, $c, $d];
    BetaTensor[3, 2, 33] := BetaTensor[3, 2, 33] = Tr @ Tdot[Yuk[$a, k5, k1], YukTil[$b, k1, k2], Yuk[$c, k2, k3], YukTil[$d, k3, k4], Y2F[k4, k5]] // TsSym[$a, $b, $c, $d];

(* 3-loop *)
	BetaTensor[3, 3, 1] := BetaTensor[3, 3, 1] = TsSym[$a, $b, $c, $d]@ Ttimes[Tscal[A1, b1, b2], Tscal[A2, b3, b4], Tscal[A3, $a, b1], Tscal[A4, $b, b3], TsG2[A1, b4, b5], TsG2[A2, b2, b6], TsG2[A3, b5, $d], TsG2[A4, b6, $c]];
	BetaTensor[3, 3, 2] := BetaTensor[3, 3, 2] = TsSym[$a, $b, $c, $d]@ Ttimes[Tscal[A1, b1, b2], Tscal[A2, b3, b4], Tscal[A3, $a, b3], Tscal[A4, $b, b1], TsG2[A1, b4, b5], TsG2[A2, b2, b6], TsG2[A3, b5, $d], TsG2[A4, b6, $c]];
	BetaTensor[3, 3, 3] := BetaTensor[3, 3, 3] = TsSym[$a, $b, $c, $d]@ Ttimes[Tscal[A1, b1, $d], Tscal[A2, b2, $c], Tscal[A3, $a, b1], Tscal[A4, $b, b2], TsG2[A1, b3, b4], TsG2[A2, b5, b6], TsG2[A3, b6, b3], TsG2[A4, b4, b5]];
	BetaTensor[3, 3, 4] := BetaTensor[3, 3, 4] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[TfermTil[A1, k1, k2], TfermTil[A2, k2, k3], TfermTil[A3, k3, k4], TfermTil[A4, k4, k1]]], TsG2[A1, $a, b1], TsG2[A2, b1, $d], TsG2[A3, $b, b2], TsG2[A4, b2, $c]];
	BetaTensor[3, 3, 5] := BetaTensor[3, 3, 5] = TsSym[$a, $b, $c, $d]@ Ttimes[S2SC2S[A1, A2], Tscal[A3, b1, $c], TsG2[A1, $a, b2], TsG2[A2, $b, b1], TsG2[A3, b2, $d]];
	BetaTensor[3, 3, 6] := BetaTensor[3, 3, 6] = TsSym[$a, $b, $c, $d]@ Ttimes[S2FC2F[A1, A2], Tscal[A3, b1, $c], TsG2[A1, $b, b1], TsG2[A2, $a, b2], TsG2[A3, b2, $d]];
	BetaTensor[3, 3, 7] := BetaTensor[3, 3, 7] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $a], Tscal[A1, b2, $d], Tscal[A2, $b, b3], Tscal[A3, $c, b2], TsG2[A1, b4, b5], TsG2[A2, b5, b1], TsG2[A3, b3, b4]];
	BetaTensor[3, 3, 8] := BetaTensor[3, 3, 8] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $b], C2S[b2, $a], Tscal[A1, $d, b3], Tscal[A2, b3, b2], TsG2[A1, $c, b4], TsG2[A2, b4, b1]];
	BetaTensor[3, 3, 9] := BetaTensor[3, 3, 9] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $b], C2S[$a, b2], Tscal[A1, $c, b3], Tscal[A2, b3, $d], TsG2[A1, b4, b1], TsG2[A2, b2, b4]];
	BetaTensor[3, 3, 10] := BetaTensor[3, 3, 10] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], Tscal[A3, $c, b1], Tscal[A4, b2, b3], TsG2[A1, b3, $b], TsG2[A2, $a, b4], TsG2[A3, b4, b2], TsG2[A4, b1, $d]];
	BetaTensor[3, 3, 11] := BetaTensor[3, 3, 11] = TsSym[$a, $b, $c, $d]@ Ttimes[S2S[A1, A2], Tscal[A3, b1, $b], Tscal[A4, $a, b1], TsG2[A1, $d, b2], TsG2[A2, b3, $c], TsG2[A3, b2, b4], TsG2[A4, b4, b3]];
	BetaTensor[3, 3, 12] := BetaTensor[3, 3, 12] = TsSym[$a, $b, $c, $d]@ Ttimes[S2F[A1, A2], Tscal[A3, $c, b1], Tscal[A4, b2, b3], TsG2[A1, b3, $b], TsG2[A2, $a, b4], TsG2[A3, b4, b2], TsG2[A4, b1, $d]];
	BetaTensor[3, 3, 13] := BetaTensor[3, 3, 13] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], C2S[b1, $a], Tscal[A3, b2, $d], TsG2[A1, $c, b2], TsG2[A2, $b, b3], TsG2[A3, b3, b1]];
	BetaTensor[3, 3, 14] := BetaTensor[3, 3, 14] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SC2G[b1, $a], Tscal[A1, $b, b2], Tscal[A2, $c, b3], TsG2[A1, b3, $d], TsG2[A2, b2, b1]];
	BetaTensor[3, 3, 15] := BetaTensor[3, 3, 15] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $a], S2S[A1, A2], Tscal[A3, $b, b2], TsG2[A1, b2, $d], TsG2[A2, $c, b3], TsG2[A3, b3, b1]];
	BetaTensor[3, 3, 16] := BetaTensor[3, 3, 16] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2S[$d, b1], Tscal[A1, b2, $b], Tscal[A2, $a, b2], TsG2[A1, b1, b3], TsG2[A2, b3, $c]];
	BetaTensor[3, 3, 17] := BetaTensor[3, 3, 17] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $a], S2F[A1, A2], Tscal[A3, b2, $d], TsG2[A1, $c, b2], TsG2[A2, $b, b3], TsG2[A3, b3, b1]];
	BetaTensor[3, 3, 18] := BetaTensor[3, 3, 18] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2F[b1, $a], Tscal[A1, $b, b2], Tscal[A2, $c, b3], TsG2[A1, b3, $d], TsG2[A2, b2, b1]];
	BetaTensor[3, 3, 19] := BetaTensor[3, 3, 19] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], C2G[A3, A4], TsG2[A1, b1, $d], TsG2[A2, $a, b2], TsG2[A3, b2, $c], TsG2[A4, $b, b1]];
	BetaTensor[3, 3, 20] := BetaTensor[3, 3, 20] = TsSym[$a, $b, $c, $d]@ Ttimes[S2S[A1, A2], S2S[A3, A4], TsG2[A1, b1, $d], TsG2[A2, $a, b2], TsG2[A3, b2, $c], TsG2[A4, $b, b1]];
	BetaTensor[3, 3, 21] := BetaTensor[3, 3, 21] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], S2S[A3, A1], S2S[A4, A2], Tscal[A5, $b, b1], TsG2[A3, b2, $d], TsG2[A4, b1, $c], TsG2[A5, $a, b2]];
	BetaTensor[3, 3, 22] := BetaTensor[3, 3, 22] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2S[A3, A4], Tscal[A5, $b, b1], TsG2[A2, $a, b2], TsG2[A4, b1, $d], TsG2[A5, b2, $c]];
	BetaTensor[3, 3, 23] := BetaTensor[3, 3, 23] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], C2G[A3, A4], G2Matrix[A1, A3], Tscal[A5, b1, $c], TsG2[A2, $b, b1], TsG2[A4, $a, b2], TsG2[A5, b2, $d]];
	BetaTensor[3, 3, 24] := BetaTensor[3, 3, 24] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], S2S[A3, A4], TsG2[A1, $b, b1], TsG2[A2, $a, b2], TsG2[A3, b2, $d], TsG2[A4, b1, $c]];
	BetaTensor[3, 3, 25] := BetaTensor[3, 3, 25] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], S2F[A3, A1], S2S[A2, A4], Tscal[A5, $b, b1], TsG2[A3, $a, b2], TsG2[A4, b1, $d], TsG2[A5, b2, $c]];
	BetaTensor[3, 3, 26] := BetaTensor[3, 3, 26] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], S2F[A3, A4], Tscal[A5, b1, $c], TsG2[A2, $a, b2], TsG2[A4, $b, b1], TsG2[A5, b2, $d]];
	BetaTensor[3, 3, 27] := BetaTensor[3, 3, 27] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], S2F[A3, A4], TsG2[A1, $b, b1], TsG2[A2, $a, b2], TsG2[A3, b2, $d], TsG2[A4, b1, $c]];
	BetaTensor[3, 3, 28] := BetaTensor[3, 3, 28] = TsSym[$a, $b, $c, $d]@ Ttimes[S2F[A1, A2], S2S[A3, A4], TsG2[A1, $b, b1], TsG2[A2, $a, b2], TsG2[A3, b2, $d], TsG2[A4, b1, $c]];
	BetaTensor[3, 3, 29] := BetaTensor[3, 3, 29] = TsSym[$a, $b, $c, $d]@ Ttimes[S2F[A1, A2], S2F[A3, A4], TsG2[A1, b1, $d], TsG2[A2, $a, b2], TsG2[A3, b2, $c], TsG2[A4, $b, b1]];
	BetaTensor[3, 3, 30] := BetaTensor[3, 3, 30] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], S2F[A1, A3], S2F[A2, A4], Tscal[A5, b1, $c], TsG2[A3, $b, b1], TsG2[A4, $a, b2], TsG2[A5, b2, $d]];
	BetaTensor[3, 3, 31] := BetaTensor[3, 3, 31] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Tscal[A1, $a, b4], Tscal[A2, $b, b5], Tscal[A3, b1, b6], TsG2[A1, b5, $c], TsG2[A2, b6, b2], TsG2[A3, b4, b3]];
	BetaTensor[3, 3, 32] := BetaTensor[3, 3, 32] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Tscal[A1, b3, b1], Tscal[A2, $c, b3], Tscal[A3, b2, b4], TsG2[A1, b4, b5], TsG2[A2, b5, b6], TsG2[A3, b6, $d]];
	BetaTensor[3, 3, 33] := BetaTensor[3, 3, 33] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Tscal[A1, b3, b4], Tscal[A2, $c, b3], Tscal[A3, b2, b5], TsG2[A1, b5, b6], TsG2[A2, b6, $d], TsG2[A3, b4, b1]];
	BetaTensor[3, 3, 34] := BetaTensor[3, 3, 34] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Tscal[A1, b3, b2], Tscal[A2, $c, b4], Tscal[A3, b1, b3], TsG2[A1, b5, b6], TsG2[A2, b6, $d], TsG2[A3, b4, b5]];
	BetaTensor[3, 3, 35] := BetaTensor[3, 3, 35] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $c], Lam[b2, b3, $a, $b], Tscal[A1, $d, b4], Tscal[A2, b2, b5], TsG2[A1, b5, b1], TsG2[A2, b4, b3]];
	BetaTensor[3, 3, 36] := BetaTensor[3, 3, 36] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b2, b3, $a, $b], Tscal[A1, $c, b4], Tscal[A2, b5, b1], TsG2[A1, b3, b5], TsG2[A2, b4, $d]];
	BetaTensor[3, 3, 37] := BetaTensor[3, 3, 37] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$a, b1], Lam[b1, b2, b3, $d], Tscal[A1, $b, b4], Tscal[A2, b2, b5], TsG2[A1, b5, b3], TsG2[A2, b4, $c]];
	BetaTensor[3, 3, 38] := BetaTensor[3, 3, 38] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $a, $b, $c], S2SC2S[A1, A2], TsG2[A1, $d, b2], TsG2[A2, b2, b1]];
	BetaTensor[3, 3, 39] := BetaTensor[3, 3, 39] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $c], Lam[b2, b3, $a, $b], Tscal[A1, $d, b4], Tscal[A2, b2, b5], TsG2[A1, b5, b3], TsG2[A2, b4, b1]];
	BetaTensor[3, 3, 40] := BetaTensor[3, 3, 40] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $a, $b, $c], S2FC2F[A1, A2], TsG2[A1, b2, b1], TsG2[A2, $d, b2]];
	BetaTensor[3, 3, 41] := BetaTensor[3, 3, 41] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], C2S[b2, $c], Lam[b3, b4, $a, $b], Tscal[A1, b4, b1], TsG2[A1, b3, b2]];
	BetaTensor[3, 3, 42] := BetaTensor[3, 3, 42] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$a, b1], C2S[$b, b2], C2S[$c, b3], Lam[b1, b2, b3, $d]];
	BetaTensor[3, 3, 43] := BetaTensor[3, 3, 43] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], C2S[b2, b1], C2S[$c, b3], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 44] := BetaTensor[3, 3, 44] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[b2, $a], C2S[b3, b1], Lam[b3, $b, $c, $d]];
	BetaTensor[3, 3, 45] := BetaTensor[3, 3, 45] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], Lam[b1, b2, $a, $b], Tscal[A3, b3, b1], TsG2[A1, $d, b3], TsG2[A2, $c, b4], TsG2[A3, b4, b2]];
	BetaTensor[3, 3, 46] := BetaTensor[3, 3, 46] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], S2S[A1, A2], Tscal[A3, b3, b1], TsG2[A1, $d, b3], TsG2[A2, $c, b4], TsG2[A3, b4, b2]];
	BetaTensor[3, 3, 47] := BetaTensor[3, 3, 47] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], S2F[A1, A2], Tscal[A3, b3, b1], TsG2[A1, $d, b3], TsG2[A2, $c, b4], TsG2[A3, b4, b2]];
	BetaTensor[3, 3, 48] := BetaTensor[3, 3, 48] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], Lam[b1, b2, $a, $b], Tscal[A3, b3, b2], TsG2[A1, b1, b3], TsG2[A2, $c, b4], TsG2[A3, b4, $d]];
	BetaTensor[3, 3, 49] := BetaTensor[3, 3, 49] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], S2S[A1, A2], Tscal[A3, $c, b3], TsG2[A1, b4, b2], TsG2[A2, b3, $d], TsG2[A3, b1, b4]];
	BetaTensor[3, 3, 50] := BetaTensor[3, 3, 50] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], S2F[A1, A2], Tscal[A3, b3, b2], TsG2[A1, b1, b3], TsG2[A2, $c, b4], TsG2[A3, b4, $d]];
	BetaTensor[3, 3, 51] := BetaTensor[3, 3, 51] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], C2SC2G[$d, b2], Lam[b1, b2, $a, $b]];
	BetaTensor[3, 3, 52] := BetaTensor[3, 3, 52] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], C2SS2S[$d, b2], Lam[b1, b2, $a, $b]];
	BetaTensor[3, 3, 53] := BetaTensor[3, 3, 53] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $a], C2SC2G[b2, b1], Lam[b2, $b, $c, $d]];
	BetaTensor[3, 3, 54] := BetaTensor[3, 3, 54] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $a], C2SS2S[b2, b1], Lam[b2, $b, $c, $d]];
	BetaTensor[3, 3, 55] := BetaTensor[3, 3, 55] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], C2SS2F[$d, b2], Lam[b1, b2, $a, $b]];
	BetaTensor[3, 3, 56] := BetaTensor[3, 3, 56] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $a], C2SS2F[b2, b1], Lam[b2, $b, $c, $d]];
	BetaTensor[3, 3, 57] := BetaTensor[3, 3, 57] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, $a, $b, $c], S2S[A3, A1], S2S[A4, A2], TsG2[A3, b2, b1], TsG2[A4, $d, b2]];
	BetaTensor[3, 3, 58] := BetaTensor[3, 3, 58] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], Lam[b1, $a, $b, $c], S2S[A3, A4], TsG2[A2, $d, b2], TsG2[A4, b2, b1]];
	BetaTensor[3, 3, 59] := BetaTensor[3, 3, 59] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], C2G[A3, A4], G2Matrix[A1, A3], Lam[b1, $a, $b, $c], TsG2[A2, b2, b1], TsG2[A4, $d, b2]];
	BetaTensor[3, 3, 60] := BetaTensor[3, 3, 60] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, $a, $b, $c], S2F[A3, A1], S2S[A2, A4], TsG2[A3, $d, b2], TsG2[A4, b2, b1]];
	BetaTensor[3, 3, 61] := BetaTensor[3, 3, 61] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], Lam[b1, $a, $b, $c], S2F[A3, A4], TsG2[A2, $d, b2], TsG2[A4, b2, b1]];
	BetaTensor[3, 3, 62] := BetaTensor[3, 3, 62] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, $a, $b, $c], S2F[A1, A3], S2F[A2, A4], TsG2[A3, b2, b1], TsG2[A4, $d, b2]];
	BetaTensor[3, 3, 63] := BetaTensor[3, 3, 63] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b3, b4, $a, $b], Tscal[A1, $c, b5], Tscal[A2, b1, b6], TsG2[A1, b6, b2], TsG2[A2, b5, b4]];
	BetaTensor[3, 3, 64] := BetaTensor[3, 3, 64] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b3, b4, $a, $b], Tscal[A1, $c, b5], Tscal[A2, b1, b6], TsG2[A1, b6, b4], TsG2[A2, b5, b2]];
	BetaTensor[3, 3, 65] := BetaTensor[3, 3, 65] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b2, b4, b3, $c], Tscal[A1, $b, b5], Tscal[A2, b5, b1], TsG2[A1, $a, b6], TsG2[A2, b6, b4]];
	BetaTensor[3, 3, 66] := BetaTensor[3, 3, 66] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $a, $b], Tscal[A1, $d, b5], Tscal[A2, b5, b2], TsG2[A1, $c, b6], TsG2[A2, b6, b3]];
	BetaTensor[3, 3, 67] := BetaTensor[3, 3, 67] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $a, $b, $c], Lam[b2, b3, b4, $d], Tscal[A1, b5, b4], Tscal[A2, b3, b5], TsG2[A1, b2, b6], TsG2[A2, b6, b1]];
	BetaTensor[3, 3, 68] := BetaTensor[3, 3, 68] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam2[b1, $a], Tscal[A1, $b, b2], Tscal[A2, $c, b3], TsG2[A1, b3, $d], TsG2[A2, b2, b1]];
	BetaTensor[3, 3, 69] := BetaTensor[3, 3, 69] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b2, b4, b3, $c], Tscal[A1, $a, b5], Tscal[A2, b1, b6], TsG2[A1, b6, b4], TsG2[A2, b5, $b]];
	BetaTensor[3, 3, 70] := BetaTensor[3, 3, 70] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $a, $b], Tscal[A1, $c, b5], Tscal[A2, b2, b6], TsG2[A1, b6, b3], TsG2[A2, b5, $d]];
	BetaTensor[3, 3, 71] := BetaTensor[3, 3, 71] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b3, b4, $a, $b], Tscal[A1, b5, b4], Tscal[A2, b3, b5], TsG2[A1, b1, b6], TsG2[A2, b6, b2]];
	BetaTensor[3, 3, 72] := BetaTensor[3, 3, 72] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b2, b3, $a, $b], Lam[b3, b4, b5, $d], Tscal[A1, $c, b4], TsG2[A1, b5, b1]];
	BetaTensor[3, 3, 73] := BetaTensor[3, 3, 73] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $c], Lam[b2, b3, b4, b5], Lam[b2, b5, $a, $b], Tscal[A1, $d, b3], TsG2[A1, b4, b1]];
	BetaTensor[3, 3, 74] := BetaTensor[3, 3, 74] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[b3, b4], Lam[b1, b3, $c, $d], Lam[b2, b4, $a, $b]];
	BetaTensor[3, 3, 75] := BetaTensor[3, 3, 75] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[b3, b1], Lam[b2, b4, $c, $d], Lam[b3, b4, $a, $b]];
	BetaTensor[3, 3, 76] := BetaTensor[3, 3, 76] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[$c, b3], Lam[b1, b3, b4, $d], Lam[b2, b4, $a, $b]];
	BetaTensor[3, 3, 77] := BetaTensor[3, 3, 77] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$a, b1], C2S[$b, b2], Lam[b1, b3, b4, $c], Lam[b2, b3, b4, $d]];
	BetaTensor[3, 3, 78] := BetaTensor[3, 3, 78] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], C2S[$d, b2], Lam[b1, b2, b3, b4], Lam[b3, b4, $a, $b]];
	BetaTensor[3, 3, 79] := BetaTensor[3, 3, 79] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $c], C2S[b2, b1], Lam[b2, b3, b4, $d], Lam[b3, b4, $a, $b]];
	BetaTensor[3, 3, 80] := BetaTensor[3, 3, 80] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], Lam[b1, b2, b3, b4], Lam[b1, b4, $a, $b], TsG2[A1, $d, b3], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 81] := BetaTensor[3, 3, 81] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $a, $b], S2S[A1, A2], TsG2[A1, $d, b3], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 82] := BetaTensor[3, 3, 82] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $a, $b], S2F[A1, A2], TsG2[A1, $d, b3], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 83] := BetaTensor[3, 3, 83] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SC2G[b1, b2], Lam[b1, b3, $c, $d], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 84] := BetaTensor[3, 3, 84] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2S[b1, b2], Lam[b1, b3, $c, $d], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 85] := BetaTensor[3, 3, 85] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2F[b1, b2], Lam[b1, b3, $c, $d], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 86] := BetaTensor[3, 3, 86] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SC2G[$c, b1], Lam[b1, b2, b3, $d], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 87] := BetaTensor[3, 3, 87] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2S[$c, b1], Lam[b1, b2, b3, $d], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 88] := BetaTensor[3, 3, 88] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2F[$c, b1], Lam[b1, b2, b3, $d], Lam[b2, b3, $a, $b]];
	BetaTensor[3, 3, 89] := BetaTensor[3, 3, 89] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b3, b5, b4, $d], Lam[b6, b5, $a, $b], Tscal[A1, $c, b1], TsG2[A1, b6, b2]];
	BetaTensor[3, 3, 90] := BetaTensor[3, 3, 90] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b3, b4, b5, $c], Lam[b6, b4, b5, $b], Tscal[A1, $a, b1], TsG2[A1, b2, b6]];
	BetaTensor[3, 3, 91] := BetaTensor[3, 3, 91] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b3, $a, $b], Lam[b2, b5, b6, b4], Tscal[A1, $d, b6], TsG2[A1, $c, b5]];
	BetaTensor[3, 3, 92] := BetaTensor[3, 3, 92] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, b3, b4, $d], Lam[b2, b4, b5, $c], Lam[b3, b5, $a, $b]];
	BetaTensor[3, 3, 93] := BetaTensor[3, 3, 93] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$a, b1], Lam[b2, $b, $c, $d], Lam2[b1, b2]];
	BetaTensor[3, 3, 94] := BetaTensor[3, 3, 94] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, b3, b4, b5], Lam[b2, b4, b5, $d], Lam[b3, $a, $b, $c]];
	BetaTensor[3, 3, 95] := BetaTensor[3, 3, 95] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, b3, $a, $b], Lam[b2, b4, b5, $c], Lam[b3, b4, b5, $d]];
	BetaTensor[3, 3, 96] := BetaTensor[3, 3, 96] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, b3, b4, b5], Lam[b2, b4, $c, $d], Lam[b3, b5, $a, $b]];
	BetaTensor[3, 3, 97] := BetaTensor[3, 3, 97] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], Lam[b1, b2, b3, b4], Lam[b2, b5, $a, $b], Lam[b3, b5, b4, $d]];
	BetaTensor[3, 3, 98] := BetaTensor[3, 3, 98] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$a, b1], Lam[b1, b2, b3, $d], Lam[b2, b4, b5, $b], Lam[b3, b4, b5, $c]];
	BetaTensor[3, 3, 99] := BetaTensor[3, 3, 99] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], Lam[b1, b2, b3, $d], Lam[b2, b4, b3, b5], Lam[b4, b5, $a, $b]];
	BetaTensor[3, 3, 100] := BetaTensor[3, 3, 100] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $a], Lam[b1, b4, b5, $b], Lam[b2, b4, b6, $c], Lam[b3, b5, b6, $d]];
	BetaTensor[3, 3, 101] := BetaTensor[3, 3, 101] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b2, b3, $a, $b], Lam2[b1, b3]];
	BetaTensor[3, 3, 102] := BetaTensor[3, 3, 102] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b5, $a, $b], Lam[b2, b6, $c, $d], Lam[b3, b5, b6, b4]];
	BetaTensor[3, 3, 103] := BetaTensor[3, 3, 103] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b5, $a, $b], Lam[b2, b6, b4, $c], Lam[b3, b5, b6, $d]];
	BetaTensor[3, 3, 104] := BetaTensor[3, 3, 104] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, $a, $b, $c], Lam[b2, b5, b4, b6], Lam[b3, b5, b6, $d]];
	BetaTensor[3, 3, 105] := BetaTensor[3, 3, 105] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b5, b3, $c], Lam[b2, b6, b4, $d], Lam[b5, b6, $a, $b]];
	BetaTensor[3, 3, 106] := BetaTensor[3, 3, 106] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $a], Lam[b1, b4, b5, $c], Lam[b2, b6, b3, $d], Lam[b4, b6, b5, $b]];
	BetaTensor[3, 3, 107] := BetaTensor[3, 3, 107] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $a, $b], Lam[b2, b5, b6, $c], Lam[b3, b5, b6, $d]];
	BetaTensor[3, 3, 108] := BetaTensor[3, 3, 108] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b3, $c, $d], Lam[b2, b5, b4, b6], Lam[b5, b6, $a, $b]];
	BetaTensor[3, 3, 109] := BetaTensor[3, 3, 109] = TsSym[$a, $b, $c, $d]@ Ttimes[S2FY2F[A1, A2], Tscal[A3, b1, $c], TsG2[A1, $b, b1], TsG2[A2, $a, b2], TsG2[A3, b2, $d]];
	BetaTensor[3, 3, 110] := BetaTensor[3, 3, 110] = TsSym[$a, $b, $c, $d]@ Ttimes[S2SY2S[A1, A2], Tscal[A3, $b, b1], TsG2[A1, b2, $d], TsG2[A2, b1, $c], TsG2[A3, $a, b2]];
	BetaTensor[3, 3, 111] := BetaTensor[3, 3, 111] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[$a, k4, k2]]], Tscal[A3, $d, b1], TsG2[A1, b2, b3], TsG2[A2, b1, b2], TsG2[A3, b3, $c]];
	BetaTensor[3, 3, 112] := BetaTensor[3, 3, 112] = TsSym[$a, $b, $c, $d]@ Ttimes[Tscal[A1, b1, b2], Tscal[A2, $d, b3], Tscal[A3, $b, b1], TsG2[A1, b3, b4], TsG2[A2, b2, b5], TsG2[A3, b5, $c], Y2S[b4, $a]];
	BetaTensor[3, 3, 113] := BetaTensor[3, 3, 113] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $c], Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[$a, k4, k2]]], TsG2[A1, b2, b1], TsG2[A2, $d, b2]];
	BetaTensor[3, 3, 114] := BetaTensor[3, 3, 114] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $b], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[$a, k4, k2]]], TsG2[A1, b2, $d], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 115] := BetaTensor[3, 3, 115] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$b, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], Tferm[A1, k4, k5], Tferm[A2, k5, k2]]], TsG2[A1, $c, b1], TsG2[A2, b1, $d]];
	BetaTensor[3, 3, 116] := BetaTensor[3, 3, 116] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[TfermTil[A1, k1, k2], C2Ft[k2, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k1]]], TsG2[A1, b1, $d], TsG2[A2, $c, b1]];
	BetaTensor[3, 3, 117] := BetaTensor[3, 3, 117] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k2]]], Tscal[A1, b2, $c], Tscal[A2, $b, b2], TsG2[A1, b1, b3], TsG2[A2, b3, $d]];
	BetaTensor[3, 3, 118] := BetaTensor[3, 3, 118] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $b], Tscal[A1, b2, b1], Tscal[A2, $d, b3], TsG2[A1, b3, b4], TsG2[A2, $c, b2], Y2S[b4, $a]];
	BetaTensor[3, 3, 119] := BetaTensor[3, 3, 119] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $b], Tscal[A1, b2, $d], Tscal[A2, $c, b2], TsG2[A1, b3, b4], TsG2[A2, b4, b1], Y2S[b3, $a]];
	BetaTensor[3, 3, 120] := BetaTensor[3, 3, 120] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], S2S[A1, A3], Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A4, k1, k3], TfermTil[A2, k3, k4], Yuk[$a, k4, k2]]], TsG2[A3, b1, $d], TsG2[A4, $c, b1]];
	BetaTensor[3, 3, 121] := BetaTensor[3, 3, 121] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A3, k1, k3], TfermTil[A4, k3, k4], Yuk[$a, k4, k2]]], TsG2[A2, $c, b1], TsG2[A4, b1, $d]];
	BetaTensor[3, 3, 122] := BetaTensor[3, 3, 122] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], S2F[A1, A3], Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A2, k1, k3], TfermTil[A4, k3, k4], Yuk[$a, k4, k2]]], TsG2[A3, b1, $d], TsG2[A4, $c, b1]];
	BetaTensor[3, 3, 123] := BetaTensor[3, 3, 123] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], Tscal[A3, b1, b2], TsG2[A1, $c, b1], TsG2[A2, $b, b3], TsG2[A3, b3, $d], Y2S[b2, $a]];
	BetaTensor[3, 3, 124] := BetaTensor[3, 3, 124] = TsSym[$a, $b, $c, $d]@ Ttimes[S2S[A1, A2], Tscal[A3, $b, b1], TsG2[A1, $d, b2], TsG2[A2, b1, $c], TsG2[A3, b2, b3], Y2S[b3, $a]];
	BetaTensor[3, 3, 125] := BetaTensor[3, 3, 125] = TsSym[$a, $b, $c, $d]@ Ttimes[S2F[A1, A2], Tscal[A3, b1, b2], TsG2[A1, $c, b1], TsG2[A2, $b, b3], TsG2[A3, b3, $d], Y2S[b2, $a]];
	BetaTensor[3, 3, 126] := BetaTensor[3, 3, 126] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $b, $c], Tr[Tdot[YukTil[b2, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[$a, k4, k2]]], TsG2[A1, b3, b1], TsG2[A2, $d, b3]];
	BetaTensor[3, 3, 127] := BetaTensor[3, 3, 127] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Tr[Tdot[YukTil[b2, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[b1, k4, k2]]], TsG2[A1, b3, $d], TsG2[A2, $c, b3]];
	BetaTensor[3, 3, 128] := BetaTensor[3, 3, 128] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Tscal[A1, b3, $d], Tscal[A2, $c, b3], TsG2[A1, b4, b5], TsG2[A2, b2, b4], Y2S[b5, b1]];
	BetaTensor[3, 3, 129] := BetaTensor[3, 3, 129] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $b, $c], Tscal[A1, $d, b3], Tscal[A2, b2, b4], TsG2[A1, b4, b5], TsG2[A2, b3, b1], Y2S[b5, $a]];
	BetaTensor[3, 3, 130] := BetaTensor[3, 3, 130] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[$a, k4, k2]]], TsG2[A1, b3, b2], TsG2[A2, b1, b3]];
	BetaTensor[3, 3, 131] := BetaTensor[3, 3, 131] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $a, $b, $c], S2FY2F[A1, A2], TsG2[A1, b2, b1], TsG2[A2, $d, b2]];
	BetaTensor[3, 3, 132] := BetaTensor[3, 3, 132] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $a, $b, $c], S2SY2S[A1, A2], TsG2[A1, b2, b1], TsG2[A2, $d, b2]];
	BetaTensor[3, 3, 133] := BetaTensor[3, 3, 133] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $b, $c], Tscal[A1, b3, b2], Tscal[A2, $d, b4], TsG2[A1, b4, b5], TsG2[A2, b1, b3], Y2S[b5, $a]];
	BetaTensor[3, 3, 134] := BetaTensor[3, 3, 134] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Lam[b2, b3, $b, $c], Tscal[A1, b3, b1], TsG2[A1, b2, b4], Y2S[b4, $a]];
	BetaTensor[3, 3, 135] := BetaTensor[3, 3, 135] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Lam[b2, $a, $b, $c], Tr[Tdot[YukTil[b2, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k2]]]];
	BetaTensor[3, 3, 136] := BetaTensor[3, 3, 136] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k4, k3], C2F[k4, k2]]]];
	BetaTensor[3, 3, 137] := BetaTensor[3, 3, 137] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[C2Ft[k1, k2], C2Ft[k2, k3], Yuk[b1, k3, k4], YukTil[$a, k1, k4]]]];
	BetaTensor[3, 3, 138] := BetaTensor[3, 3, 138] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[$a, b3], Lam[b1, $b, $c, $d], Y2S[b3, b2]];
	BetaTensor[3, 3, 139] := BetaTensor[3, 3, 139] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[$d, b3], Lam[b1, b3, $b, $c], Y2S[b2, $a]];
	BetaTensor[3, 3, 140] := BetaTensor[3, 3, 140] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2FC2Gt[k1, k3], Yuk[$a, k3, k2]]]];
	BetaTensor[3, 3, 141] := BetaTensor[3, 3, 141] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2FS2St[k1, k3], Yuk[$a, k3, k2]]]];
	BetaTensor[3, 3, 142] := BetaTensor[3, 3, 142] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2FS2Ft[k1, k3], Yuk[$a, k3, k2]]]];
	BetaTensor[3, 3, 143] := BetaTensor[3, 3, 143] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SC2G[b1, b2], Lam[b1, $b, $c, $d], Y2S[b2, $a]];
	BetaTensor[3, 3, 144] := BetaTensor[3, 3, 144] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2S[b1, b2], Lam[b1, $b, $c, $d], Y2S[b2, $a]];
	BetaTensor[3, 3, 145] := BetaTensor[3, 3, 145] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2F[b1, b2], Lam[b1, $b, $c, $d], Y2S[b2, $a]];
	BetaTensor[3, 3, 146] := BetaTensor[3, 3, 146] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Lam[b2, b3, b4, $d], Tscal[A1, $c, b4], TsG2[A1, b3, b5], Y2S[b5, b1]];
	BetaTensor[3, 3, 147] := BetaTensor[3, 3, 147] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $b, $c], Tscal[A1, $d, b3], TsG2[A1, b2, b5], Y2S[b5, $a]];
	BetaTensor[3, 3, 148] := BetaTensor[3, 3, 148] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b2, b3, $a, $b], Tr[Tdot[YukTil[b3, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k2]]]];
	BetaTensor[3, 3, 149] := BetaTensor[3, 3, 149] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, b3, $c, $d], Lam[b4, b3, $a, $b], Y2S[b2, b4]];
	BetaTensor[3, 3, 150] := BetaTensor[3, 3, 150] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, b3, $c, $d], Lam[b2, b4, $a, $b], Y2S[b3, b4]];
	BetaTensor[3, 3, 151] := BetaTensor[3, 3, 151] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], Lam[b1, b2, b3, $d], Lam[b3, b4, $a, $b], Y2S[b2, b4]];
	BetaTensor[3, 3, 152] := BetaTensor[3, 3, 152] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], Lam[b2, b3, b4, $d], Lam[b3, b4, $a, $b], Y2S[b1, b2]];
	BetaTensor[3, 3, 153] := BetaTensor[3, 3, 153] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b2, b4, $a, $b], Lam[b3, b5, b4, $c], Y2S[b1, b5]];
	BetaTensor[3, 3, 154] := BetaTensor[3, 3, 154] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b2, $a, $b, $c], Lam[b3, b5, b4, $d], Y2S[b1, b5]];
	BetaTensor[3, 3, 155] := BetaTensor[3, 3, 155] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $a, $b], Lam[b2, b3, b4, $d], Lam[b5, b3, b4, $c], Y2S[b1, b5]];
	BetaTensor[3, 3, 156] := BetaTensor[3, 3, 156] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b2, b4, $a, $b], Lam[b3, b5, $c, $d], Y2S[b1, b5]];
	BetaTensor[3, 3, 157] := BetaTensor[3, 3, 157] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[$a, k2, k5]]], Tscal[A2, $c, b2], TsG2[A1, $b, b3], TsG2[A2, b3, $d]];
	BetaTensor[3, 3, 158] := BetaTensor[3, 3, 158] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A1, k1, k3], Yuk[b1, k3, k4], Tferm[A2, k4, k5], YukTil[$a, k5, k6], Yuk[b1, k6, k2]]], TsG2[A1, b2, $d], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 159] := BetaTensor[3, 3, 159] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[b2, k3, k4], YukTil[$a, k5, k4], TfermTil[A2, k5, k6], Yuk[$b, k6, k2]]], Tscal[A3, $d, b2], TsG2[A3, $c, b1]];
	BetaTensor[3, 3, 160] := BetaTensor[3, 3, 160] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[A3, k3, k4], Yuk[$a, k4, k5], YukTil[$b, k6, k5], TfermTil[A2, k6, k7], Yuk[$c, k7, k2]]], TsG2[A3, $d, b1]];
	BetaTensor[3, 3, 161] := BetaTensor[3, 3, 161] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A1, k1, k3], TfermTil[A3, k3, k4], Yuk[$a, k4, k5], YukTil[$b, k6, k5], TfermTil[A4, k6, k7], TfermTil[A2, k7, k8], Yuk[$c, k8, k2]]]];
	BetaTensor[3, 3, 162] := BetaTensor[3, 3, 162] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]], Tscal[A1, b3, $d], Tscal[A2, $c, b3], TsG2[A1, b1, b4], TsG2[A2, b4, b2]];
	BetaTensor[3, 3, 163] := BetaTensor[3, 3, 163] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b2, k4, k2]]], Tscal[A1, b3, $c], Tscal[A2, $b, b3], TsG2[A1, b1, b4], TsG2[A2, b4, $d]];
	BetaTensor[3, 3, 164] := BetaTensor[3, 3, 164] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[b1, k4, k5], YukTil[$a, k5, k6], Yuk[$b, k2, k6]]], TsG2[A1, b2, $d], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 165] := BetaTensor[3, 3, 165] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[b1, k4, k5], YukTil[$a, k5, k6], Yuk[b1, k2, k6]]], TsG2[A1, b2, $d], TsG2[A2, $c, b2]];
	BetaTensor[3, 3, 166] := BetaTensor[3, 3, 166] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k2]]], Tscal[A1, b2, $c], Tscal[A2, $b, b2], TsG2[A1, b1, b3], TsG2[A2, b3, $d]];
	BetaTensor[3, 3, 167] := BetaTensor[3, 3, 167] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$b, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], Tferm[A1, k4, k5], Tferm[A2, k5, k2]]], TsG2[A1, $c, b1], TsG2[A2, b1, $d]];
	BetaTensor[3, 3, 168] := BetaTensor[3, 3, 168] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[TfermTil[A1, k1, k2], Y2F[k2, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k1]]], TsG2[A1, $c, b1], TsG2[A2, b1, $d]];
	BetaTensor[3, 3, 169] := BetaTensor[3, 3, 169] = TsSym[$a, $b, $c, $d]@ Ttimes[Tscal[A1, b1, $c], Tscal[A2, $d, b2], TsG2[A1, b2, b3], TsG2[A2, $b, b1], Y2S[b4, b3], Y2S[b4, $a]];
	BetaTensor[3, 3, 170] := BetaTensor[3, 3, 170] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], Tr[Tdot[YukTil[$b, k1, k2], TfermTil[A1, k1, k3], TfermTil[A3, k3, k4], Yuk[$a, k4, k2]]], Tr[Tdot[YukTil[$d, k5, k6], TfermTil[A4, k5, k7], TfermTil[A2, k7, k8], Yuk[$c, k8, k6]]]];
	BetaTensor[3, 3, 171] := BetaTensor[3, 3, 171] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$c, k1, k2], TfermTil[A1, k1, k3], TfermTil[A2, k3, k4], Yuk[$b, k4, k2]]], TsG2[A1, $d, b1], TsG2[A2, b1, b2], Y2S[b2, $a]];
	BetaTensor[3, 3, 172] := BetaTensor[3, 3, 172] = TsSym[$a, $b, $c, $d]@ Ttimes[Tscal[A1, b1, b2], Tscal[A2, $d, b3], TsG2[A1, b3, b4], TsG2[A2, $c, b1], Y2S[b2, $a], Y2S[b4, $b]];
	BetaTensor[3, 3, 173] := BetaTensor[3, 3, 173] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[$c, k6, k2]]]];
	BetaTensor[3, 3, 174] := BetaTensor[3, 3, 174] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], Tferm[A1, k4, k5], YukTil[$b, k5, k6], Yuk[$c, k7, k6], Tferm[A2, k7, k2]]]];
	BetaTensor[3, 3, 175] := BetaTensor[3, 3, 175] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[TfermTil[A1, k1, k2], C2Ft[k2, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[$c, k6, k7], YukTil[$d, k1, k7]]]];
	BetaTensor[3, 3, 176] := BetaTensor[3, 3, 176] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], C2S[b2, $d], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 177] := BetaTensor[3, 3, 177] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], C2S[b2, $c], Tr[Tdot[YukTil[$b, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k4, k2]]]];
	BetaTensor[3, 3, 178] := BetaTensor[3, 3, 178] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], C2S[b2, $c], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 179] := BetaTensor[3, 3, 179] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k2, k5]]]];
	BetaTensor[3, 3, 180] := BetaTensor[3, 3, 180] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[$c, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k4, k5], Yuk[$b, k2, k5]]]];
	BetaTensor[3, 3, 181] := BetaTensor[3, 3, 181] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], C2Ft[k5, k6], Yuk[$c, k6, k2]]];
	BetaTensor[3, 3, 182] := BetaTensor[3, 3, 182] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k4, k3], C2F[k4, k5], YukTil[$b, k5, k6], Yuk[$c, k6, k2]]];
	BetaTensor[3, 3, 183] := BetaTensor[3, 3, 183] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[C2Ft[k1, k2], C2Ft[k2, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k6, k5], YukTil[$d, k1, k6]]];
	BetaTensor[3, 3, 184] := BetaTensor[3, 3, 184] = TsSym[$a, $b, $c, $d]@ Ttimes[C2G[A1, A2], G2Matrix[A1, A3], G2Matrix[A2, A4], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A4, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A3, k5, k6], Yuk[$c, k6, k2]]]];
	BetaTensor[3, 3, 185] := BetaTensor[3, 3, 185] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], S2S[A1, A3], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A4, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[$c, k6, k2]]]];
	BetaTensor[3, 3, 186] := BetaTensor[3, 3, 186] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], G2Matrix[A3, A4], S2F[A1, A3], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A4, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[$c, k6, k2]]]];
	BetaTensor[3, 3, 187] := BetaTensor[3, 3, 187] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SC2G[b1, $d], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 188] := BetaTensor[3, 3, 188] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2S[b1, $d], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 189] := BetaTensor[3, 3, 189] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2FC2Gt[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k2, k5]]];
	BetaTensor[3, 3, 190] := BetaTensor[3, 3, 190] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2FS2St[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k2, k5]]];
	BetaTensor[3, 3, 191] := BetaTensor[3, 3, 191] = TsSym[$a, $b, $c, $d]@ Ttimes[C2SS2F[b1, $d], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 192] := BetaTensor[3, 3, 192] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2FS2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k2, k5]]];
	BetaTensor[3, 3, 193] := BetaTensor[3, 3, 193] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], TfermTil[A2, k5, k6], Yuk[$b, k6, k2]]]];
	BetaTensor[3, 3, 194] := BetaTensor[3, 3, 194] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b2, b3, $c, $d], Tr[Tdot[YukTil[$b, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 195] := BetaTensor[3, 3, 195] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], Lam[b1, b2, b3, $d], Tr[Tdot[YukTil[$b, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 196] := BetaTensor[3, 3, 196] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Lam[b2, b3, $b, $c], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b1, k3, k2], YukTil[b2, k4, k3], Yuk[$a, k4, k1]]]];
	BetaTensor[3, 3, 197] := BetaTensor[3, 3, 197] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[b2, k4, k5], Yuk[$b, k2, k5]]]];
	BetaTensor[3, 3, 198] := BetaTensor[3, 3, 198] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], TfermTil[A1, k1, k3], Yuk[b1, k3, k4], YukTil[$a, k5, k4], TfermTil[A2, k5, k6], Yuk[$b, k6, k2]]]];
	BetaTensor[3, 3, 199] := BetaTensor[3, 3, 199] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[b1, k6, k2]]]];
	BetaTensor[3, 3, 200] := BetaTensor[3, 3, 200] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[b2, k5, k4], TfermTil[A2, k5, k6], Yuk[b2, k6, k2]]]];
	BetaTensor[3, 3, 201] := BetaTensor[3, 3, 201] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $b, $c], Tscal[A1, b2, b3], TsG2[A1, b4, $d], Y2S[b3, $a], Y2S[b4, b1]];
	BetaTensor[3, 3, 202] := BetaTensor[3, 3, 202] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b2, b3, $c, $d], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 203] := BetaTensor[3, 3, 203] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$c, b1], Lam[b1, b2, b3, $d], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 204] := BetaTensor[3, 3, 204] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b3, $b, $c, $d], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k4, k2]]]];
	BetaTensor[3, 3, 205] := BetaTensor[3, 3, 205] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Lam[b2, b3, $b, $c], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 206] := BetaTensor[3, 3, 206] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Lam[b2, $a, $b, $c], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[b1, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 207] := BetaTensor[3, 3, 207] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[$b, k2, k5]]]];
	BetaTensor[3, 3, 208] := BetaTensor[3, 3, 208] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], YukTil[$a, k4, k5], Yuk[b2, k2, k5]]]];
	BetaTensor[3, 3, 209] := BetaTensor[3, 3, 209] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[b2, k2, k5]]]];
	BetaTensor[3, 3, 210] := BetaTensor[3, 3, 210] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[$b, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k2, k5]]]];
	BetaTensor[3, 3, 211] := BetaTensor[3, 3, 211] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[$a, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], YukTil[b1, k4, k5], Yuk[b2, k2, k5]]]];
	BetaTensor[3, 3, 212] := BetaTensor[3, 3, 212] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b3, $b, $c, $d], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b2, k3, k2], YukTil[b1, k4, k3], Yuk[$a, k4, k1]]]];
	BetaTensor[3, 3, 213] := BetaTensor[3, 3, 213] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Lam[b2, $a, $b, $c], Tr[Tdot[YukTil[b2, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k2]]]];
	BetaTensor[3, 3, 214] := BetaTensor[3, 3, 214] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[Y2Ft[k1, k2], C2F[k2, k3], YukTil[$a, k3, k4], Yuk[b1, k1, k4]]]];
	BetaTensor[3, 3, 215] := BetaTensor[3, 3, 215] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], Y2Ft[k4, k2]]]];
	BetaTensor[3, 3, 216] := BetaTensor[3, 3, 216] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], C2Ft[k1, k3], Yuk[b2, k3, k4], YukTil[b1, k4, k5], Yuk[$a, k2, k5]]]];
	BetaTensor[3, 3, 217] := BetaTensor[3, 3, 217] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Lam[b1, $b, $c, $d], Y2S[b3, b2], Y2S[b3, $a]];
	BetaTensor[3, 3, 218] := BetaTensor[3, 3, 218] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b3, b4, $a, $b], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k4, k2]]]];
	BetaTensor[3, 3, 219] := BetaTensor[3, 3, 219] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $b, $c], Lam[b2, b3, b4, $d], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k2], YukTil[b3, k4, k3], Yuk[$a, k4, k1]]]];
	BetaTensor[3, 3, 220] := BetaTensor[3, 3, 220] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b2, b4, b3, $c], Tr[Tdot[YukTil[$b, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b4, k4, k2]]]];
	BetaTensor[3, 3, 221] := BetaTensor[3, 3, 221] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $c, $d], Tr[Tdot[YukTil[$b, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 222] := BetaTensor[3, 3, 222] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b3, b4, $a, $b], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k1], YukTil[b2, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 223] := BetaTensor[3, 3, 223] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b2, b3, $a, $b], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b4, k3, k1], YukTil[b1, k3, k4], Yuk[b4, k4, k2]]]];
	BetaTensor[3, 3, 224] := BetaTensor[3, 3, 224] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $a, $b, $c], Lam[b2, b3, b4, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b4, k3, k4], Yuk[b3, k4, k2]]]];
	BetaTensor[3, 3, 225] := BetaTensor[3, 3, 225] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b3, b4, $b, $c], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k2], YukTil[b2, k4, k3], Yuk[$a, k4, k1]]]];
	BetaTensor[3, 3, 226] := BetaTensor[3, 3, 226] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[$a, k4, k1]]]];
	BetaTensor[3, 3, 227] := BetaTensor[3, 3, 227] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b2, b3, $a, $b], Tr[Tdot[YukTil[b3, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k2]]]];
	BetaTensor[3, 3, 228] := BetaTensor[3, 3, 228] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Lam[b2, b4, b3, $c], Tr[Tdot[YukTil[b4, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 229] := BetaTensor[3, 3, 229] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, b4], Lam[b1, b4, $c, $d], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]]];
	BetaTensor[3, 3, 230] := BetaTensor[3, 3, 230] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b3, b4, $a, $b], Y2S[b1, b4], Y2S[b2, b3]];
	BetaTensor[3, 3, 231] := BetaTensor[3, 3, 231] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Lam[b2, b3, $a, $b], Y2S[b4, b1], Y2S[b4, b3]];
	BetaTensor[3, 3, 232] := BetaTensor[3, 3, 232] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[$a, k5, k4], YukTil[b2, k6, k5], Yuk[$b, k6, k2]]], Tscal[A1, $d, b3], TsG2[A1, $c, b1]];
	BetaTensor[3, 3, 233] := BetaTensor[3, 3, 233] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$c, k1, k2], TfermTil[A1, k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[$a, k6, k5], YukTil[$b, k6, k7], Yuk[b1, k2, k7]]], TsG2[A1, $d, b2]];
	BetaTensor[3, 3, 234] := BetaTensor[3, 3, 234] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$c, k1, k2], TfermTil[A1, k1, k3], Yuk[b1, k3, k4], YukTil[b2, k4, k5], Yuk[$a, k6, k5], YukTil[b1, k6, k7], Yuk[$b, k2, k7]]], TsG2[A1, $d, b2]];
	BetaTensor[3, 3, 235] := BetaTensor[3, 3, 235] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Yuk[b1, k6, k5], Tferm[A2, k6, k7], YukTil[b1, k7, k8], Yuk[$c, k8, k2]]]];
	BetaTensor[3, 3, 236] := BetaTensor[3, 3, 236] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k6, k5], YukTil[b1, k7, k6], TfermTil[A2, k7, k8], Yuk[b1, k8, k2]]]];
	BetaTensor[3, 3, 237] := BetaTensor[3, 3, 237] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[b1, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[b1, k6, k7], YukTil[$c, k7, k8], Yuk[$d, k2, k8]]]];
	BetaTensor[3, 3, 238] := BetaTensor[3, 3, 238] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[$d, k1, k2], TfermTil[A1, k1, k3], Yuk[$a, k3, k4], YukTil[b1, k4, k5], Yuk[$b, k6, k5], YukTil[b1, k7, k6], TfermTil[A2, k7, k8], Yuk[$c, k8, k2]]]];
	BetaTensor[3, 3, 239] := BetaTensor[3, 3, 239] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[YukTil[$d, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], Tferm[A1, k4, k5], YukTil[$b, k5, k6], Yuk[$c, k7, k6], Tferm[A2, k7, k2]]]];
	BetaTensor[3, 3, 240] := BetaTensor[3, 3, 240] = TsSym[$a, $b, $c, $d]@ Ttimes[G2Matrix[A1, A2], Tr[Tdot[TfermTil[A1, k1, k2], Y2F[k2, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], TfermTil[A2, k5, k6], Yuk[$c, k6, k7], YukTil[$d, k1, k7]]]];
	BetaTensor[3, 3, 241] := BetaTensor[3, 3, 241] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$b, k4, k3], Yuk[$c, k4, k1]]], Tscal[A1, $d, b2], TsG2[A1, b1, b3], Y2S[b3, $a]];
	BetaTensor[3, 3, 242] := BetaTensor[3, 3, 242] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[$d, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k3, k4], Yuk[b1, k4, k5], YukTil[$b, k6, k5], Yuk[$c, k6, k1]]]];
	BetaTensor[3, 3, 243] := BetaTensor[3, 3, 243] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[$d, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k5], YukTil[b1, k6, k5], Yuk[$c, k1, k6]]]];
	BetaTensor[3, 3, 244] := BetaTensor[3, 3, 244] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[$b, k6, k5], Yuk[b2, k6, k2]]]];
	BetaTensor[3, 3, 245] := BetaTensor[3, 3, 245] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k4, k5], YukTil[b2, k6, k5], Yuk[$c, k6, k2]]]];
	BetaTensor[3, 3, 246] := BetaTensor[3, 3, 246] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b2, k5, k4], YukTil[$b, k6, k5], Yuk[$c, k6, k2]]]];
	BetaTensor[3, 3, 247] := BetaTensor[3, 3, 247] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k5, k4], YukTil[$c, k5, k6], Yuk[b2, k6, k2]]]];
	BetaTensor[3, 3, 248] := BetaTensor[3, 3, 248] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[$b, k6, k5], YukTil[$c, k6, k7], Yuk[b1, k2, k7]]];
	BetaTensor[3, 3, 249] := BetaTensor[3, 3, 249] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[b1, k6, k5], YukTil[$b, k7, k6], Yuk[$c, k2, k7]]];
	BetaTensor[3, 3, 250] := BetaTensor[3, 3, 250] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[$b, k5, k6], YukTil[b1, k6, k7], Yuk[$c, k2, k7]]];
	BetaTensor[3, 3, 251] := BetaTensor[3, 3, 251] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k4, k5], Yuk[$b, k6, k5], YukTil[$c, k6, k7], Yuk[b1, k2, k7]]];
	BetaTensor[3, 3, 252] := BetaTensor[3, 3, 252] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k4, k5], Yuk[$b, k6, k5], YukTil[b1, k6, k7], Yuk[$c, k2, k7]]];
	BetaTensor[3, 3, 253] := BetaTensor[3, 3, 253] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Yuk[$c, k5, k2]]]];
	BetaTensor[3, 3, 254] := BetaTensor[3, 3, 254] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, $d], Tr[Tdot[YukTil[$c, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], Yuk[$b, k2, k5]]]];
	BetaTensor[3, 3, 255] := BetaTensor[3, 3, 255] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[b1, b2], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k5], YukTil[$c, k6, k5], Yuk[$d, k6, k1]]]];
	BetaTensor[3, 3, 256] := BetaTensor[3, 3, 256] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[Y2Ft[k1, k2], C2F[k2, k3], YukTil[$a, k3, k4], Yuk[$b, k4, k5], YukTil[$c, k5, k6], Yuk[$d, k1, k6]]];
	BetaTensor[3, 3, 257] := BetaTensor[3, 3, 257] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Y2F[k5, k6], Yuk[$c, k6, k2]]];
	BetaTensor[3, 3, 258] := BetaTensor[3, 3, 258] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], C2Ft[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k6, k5], Y2Ft[k6, k2]]];
	BetaTensor[3, 3, 259] := BetaTensor[3, 3, 259] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], C2Ft[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k4, k5], Yuk[$b, k6, k5], YukTil[$c, k6, k7], Yuk[$d, k2, k7]]];
	BetaTensor[3, 3, 260] := BetaTensor[3, 3, 260] = TsSym[$a, $b, $c, $d]@ Ttimes[C2S[$d, b1], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]], Y2S[b2, b1]];
	BetaTensor[3, 3, 261] := BetaTensor[3, 3, 261] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Tr[Tdot[YukTil[$c, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k5, k4], YukTil[$b, k6, k5], Yuk[b2, k6, k2]]]];
	BetaTensor[3, 3, 262] := BetaTensor[3, 3, 262] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k5], YukTil[b3, k6, k5], Yuk[$c, k1, k6]]]];
	BetaTensor[3, 3, 263] := BetaTensor[3, 3, 263] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k5], YukTil[b1, k6, k5], Yuk[$b, k6, k2]]]];
	BetaTensor[3, 3, 264] := BetaTensor[3, 3, 264] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[$b, k6, k5], Yuk[b3, k6, k2]]]];
	BetaTensor[3, 3, 265] := BetaTensor[3, 3, 265] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k4, k5], YukTil[b3, k6, k5], Yuk[$b, k6, k2]]]];
	BetaTensor[3, 3, 266] := BetaTensor[3, 3, 266] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[$a, k4, k5], YukTil[b2, k6, k5], Yuk[b3, k6, k2]]]];
	BetaTensor[3, 3, 267] := BetaTensor[3, 3, 267] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], Yuk[$b, k5, k2]]]];
	BetaTensor[3, 3, 268] := BetaTensor[3, 3, 268] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[$b, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[b2, k4, k2]]], Y2S[b3, b1]];
	BetaTensor[3, 3, 269] := BetaTensor[3, 3, 269] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, b3, $d], Tr[Tdot[YukTil[b3, k1, k2], Yuk[b1, k3, k2], YukTil[b2, k4, k3], Yuk[$a, k5, k4], YukTil[$b, k6, k5], Yuk[$c, k6, k1]]]];
	BetaTensor[3, 3, 270] := BetaTensor[3, 3, 270] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[b1, k3, k4], Yuk[$a, k5, k4], YukTil[b3, k6, k5], Yuk[$b, k6, k2]]]];
	BetaTensor[3, 3, 271] := BetaTensor[3, 3, 271] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k5, k4], YukTil[$b, k5, k6], Yuk[b1, k6, k2]]]];
	BetaTensor[3, 3, 272] := BetaTensor[3, 3, 272] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[b1, k3, k4], Yuk[b3, k5, k4], YukTil[$a, k6, k5], Yuk[$b, k6, k2]]]];
	BetaTensor[3, 3, 273] := BetaTensor[3, 3, 273] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k4, k5], YukTil[b3, k5, k6], Yuk[b1, k6, k2]]]];
	BetaTensor[3, 3, 274] := BetaTensor[3, 3, 274] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[b3, k6, k5], Yuk[$a, k6, k2]]]];
	BetaTensor[3, 3, 275] := BetaTensor[3, 3, 275] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[b2, k5, k4], YukTil[$a, k6, k5], Yuk[b3, k6, k2]]]];
	BetaTensor[3, 3, 276] := BetaTensor[3, 3, 276] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[b3, k3, k4], Yuk[$a, k5, k4], YukTil[b3, k5, k6], Yuk[b2, k6, k2]]]];
	BetaTensor[3, 3, 277] := BetaTensor[3, 3, 277] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k5, k4], Yuk[$b, k5, k2]]]];
	BetaTensor[3, 3, 278] := BetaTensor[3, 3, 278] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[b1, k5, k2]]]];
	BetaTensor[3, 3, 279] := BetaTensor[3, 3, 279] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[$b, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k2, k5]]]];
	BetaTensor[3, 3, 280] := BetaTensor[3, 3, 280] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], YukTil[$a, k5, k4], Yuk[b2, k5, k2]]]];
	BetaTensor[3, 3, 281] := BetaTensor[3, 3, 281] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Y2F[k1, k3], Yuk[b2, k3, k4], YukTil[b1, k5, k4], Yuk[$a, k2, k5]]]];
	BetaTensor[3, 3, 282] := BetaTensor[3, 3, 282] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b2, k4, k5], Yuk[b1, k5, k2]]]];
	BetaTensor[3, 3, 283] := BetaTensor[3, 3, 283] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, b2, $c, $d], Tr[Tdot[YukTil[b2, k1, k2], Yuk[b3, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k1]]], Y2S[b3, b1]];
	BetaTensor[3, 3, 284] := BetaTensor[3, 3, 284] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[b3, k4, k3], Yuk[$a, k4, k1]]], Y2S[b3, b2]];
	BetaTensor[3, 3, 285] := BetaTensor[3, 3, 285] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b3, k4, k2]]], Y2S[b3, b2]];
	BetaTensor[3, 3, 286] := BetaTensor[3, 3, 286] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k4, k3], Y2Ft[k4, k2]]]];
	BetaTensor[3, 3, 287] := BetaTensor[3, 3, 287] = TsSym[$a, $b, $c, $d]@ Ttimes[Lam[b1, $b, $c, $d], Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Y2F[k3, k4], Yuk[$a, k4, k2]]]];
	BetaTensor[3, 3, 288] := BetaTensor[3, 3, 288] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$b, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b2, k4, k2]]], Tr[Tdot[YukTil[$d, k5, k6], Yuk[b2, k7, k5], YukTil[$c, k7, k8], Yuk[b1, k8, k6]]]];
	BetaTensor[3, 3, 289] := BetaTensor[3, 3, 289] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$c, k4, k3], Yuk[$d, k4, k1]]], Tr[Tdot[YukTil[$b, k5, k6], Yuk[b1, k7, k5], YukTil[$a, k7, k8], Yuk[b2, k8, k6]]]];
	BetaTensor[3, 3, 290] := BetaTensor[3, 3, 290] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b2, k5, k4], YukTil[$b, k5, k6], Yuk[b1, k6, k7], YukTil[$c, k8, k7], Yuk[b2, k8, k2]]];
	BetaTensor[3, 3, 291] := BetaTensor[3, 3, 291] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[$b, k5, k6], Yuk[b2, k7, k6], YukTil[$c, k8, k7], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 292] := BetaTensor[3, 3, 292] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[$b, k5, k6], Yuk[$c, k6, k7], YukTil[b2, k8, k7], Yuk[$d, k2, k8]]];
	BetaTensor[3, 3, 293] := BetaTensor[3, 3, 293] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k4, k5], YukTil[b2, k5, k6], Yuk[$b, k7, k6], YukTil[$c, k7, k8], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 294] := BetaTensor[3, 3, 294] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k5, k4], YukTil[b1, k5, k6], Yuk[b2, k6, k7], YukTil[$c, k8, k7], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 295] := BetaTensor[3, 3, 295] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Yuk[b1, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[$b, k6, k5], Yuk[b2, k6, k7], YukTil[$c, k8, k7], Yuk[b2, k8, k2]]];
	BetaTensor[3, 3, 296] := BetaTensor[3, 3, 296] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k4, k5], YukTil[b2, k5, k6], Yuk[$c, k7, k6], YukTil[b1, k7, k8], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 297] := BetaTensor[3, 3, 297] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b2, k5, k4], YukTil[$b, k6, k5], Yuk[$c, k7, k6], YukTil[b1, k8, k7], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 298] := BetaTensor[3, 3, 298] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[b1, k5, k4], YukTil[$b, k5, k6], Yuk[$c, k7, k6], YukTil[$d, k7, k8], Yuk[b2, k8, k2]]];
	BetaTensor[3, 3, 299] := BetaTensor[3, 3, 299] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k5, k4], YukTil[$c, k5, k6], Yuk[b2, k7, k6], YukTil[b1, k8, k7], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 300] := BetaTensor[3, 3, 300] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k5, k4], YukTil[b1, k5, k6], Yuk[$c, k6, k7], YukTil[$d, k7, k8], Yuk[b2, k8, k2]]];
	BetaTensor[3, 3, 301] := BetaTensor[3, 3, 301] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k4, k5], YukTil[b2, k6, k5], Yuk[b1, k7, k6], YukTil[$c, k7, k8], Yuk[$d, k8, k2]]];
	BetaTensor[3, 3, 302] := BetaTensor[3, 3, 302] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$c, k3, k4], Yuk[$d, k4, k2]]], Tr[Tdot[YukTil[b2, k5, k6], Yuk[b1, k7, k5], YukTil[$a, k7, k8], Yuk[$b, k8, k6]]]];
	BetaTensor[3, 3, 303] := BetaTensor[3, 3, 303] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k1], YukTil[$a, k3, k4], Yuk[$b, k5, k4], YukTil[$c, k5, k6], Yuk[$d, k7, k6], YukTil[b1, k7, k8], Yuk[b2, k8, k2]]];
	BetaTensor[3, 3, 304] := BetaTensor[3, 3, 304] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[b1, k3, k4], YukTil[$a, k5, k4], Yuk[$b, k5, k6], YukTil[$c, k7, k6], Yuk[$d, k7, k2]]];
	BetaTensor[3, 3, 305] := BetaTensor[3, 3, 305] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Yuk[$c, k6, k5], YukTil[b1, k7, k6], Yuk[$d, k7, k2]]];
	BetaTensor[3, 3, 306] := BetaTensor[3, 3, 306] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Yuk[b1, k5, k6], YukTil[$c, k6, k7], Yuk[$d, k7, k2]]];
	BetaTensor[3, 3, 307] := BetaTensor[3, 3, 307] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[b1, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], Yuk[$b, k5, k6], YukTil[$c, k7, k6], Yuk[$d, k7, k2]]];
	BetaTensor[3, 3, 308] := BetaTensor[3, 3, 308] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Yuk[b1, k5, k6], YukTil[$c, k7, k6], Yuk[b1, k7, k2]]];
	BetaTensor[3, 3, 309] := BetaTensor[3, 3, 309] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[b1, k5, k4], Yuk[$b, k6, k5], YukTil[$c, k7, k6], Yuk[b1, k7, k2]]];
	BetaTensor[3, 3, 310] := BetaTensor[3, 3, 310] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$d, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k3, k4], Yuk[b2, k4, k5], YukTil[$b, k6, k5], Yuk[$c, k6, k1]]], Y2S[b2, b1]];
	BetaTensor[3, 3, 311] := BetaTensor[3, 3, 311] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[$d, k1, k2], Yuk[b1, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k5], YukTil[b2, k6, k5], Yuk[$c, k1, k6]]], Y2S[b2, b1]];
	BetaTensor[3, 3, 312] := BetaTensor[3, 3, 312] = TsSym[$a, $b, $c, $d]@ Ttimes[Tr[Tdot[YukTil[b1, k1, k2], Yuk[b2, k3, k2], YukTil[$a, k4, k3], Yuk[$b, k4, k5], YukTil[$c, k6, k5], Yuk[$d, k6, k1]]], Y2S[b1, b2]];
	BetaTensor[3, 3, 313] := BetaTensor[3, 3, 313] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k4, k5], Yuk[$c, k6, k5], Y2Ft[k6, k2]]];
	BetaTensor[3, 3, 314] := BetaTensor[3, 3, 314] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Y2F[k1, k3], Yuk[$a, k3, k4], YukTil[$b, k5, k4], Y2F[k5, k6], Yuk[$c, k6, k2]]];
	BetaTensor[3, 3, 315] := BetaTensor[3, 3, 315] = TsSym[$a, $b, $c, $d]@ Tr[Tdot[YukTil[$d, k1, k2], Y2F[k1, k3], Y2F[k3, k4], Yuk[$a, k4, k5], YukTil[$b, k6, k5], Yuk[$c, k6, k2]]];


(*########################################*)
(*---------Massive Yukawa Tensors---------*)
(*########################################*)
    BetaTensor[4, 0, 1] := BetaTensor[4, 0, 1] = {{0, 0}, {0, 0}};

    BetaTensor[4, 1, 1] := BetaTensor[4, 1, 1] = {{0, 0}, {0, 0}};
    BetaTensor[4, 1, 2] := BetaTensor[4, 1, 2] = Tdot[C2Ft[$i, k1], Yuk[$a, k1, $j, True]] // TsSym[$i, $j];
    BetaTensor[4, 1, 3] := BetaTensor[4, 1, 3] = Tdot[Yuk[b, $i, k1], YukTil[$a, k1, k2, True], Yuk[b, k2, $j]];
    BetaTensor[4, 1, 4] := BetaTensor[4, 1, 4] = Tdot[Y2F[$i, k1], Yuk[$a, k1, $j, True]] // TsSym[$i, $j];
    BetaTensor[4, 1, 5] := BetaTensor[4, 1, 5] = {{0, 0}, {0, 0}};
(* y^1 terms*)
    BetaTensor[4, 2, 1] := BetaTensor[4, 2, 1] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 2] := BetaTensor[4, 2, 2] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 3] := BetaTensor[4, 2, 3] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 4] := BetaTensor[4, 2, 4] = Tdot[C2Ft[$i, k1], C2Ft[k1, k2], Yuk[$a, k2, $j, True]] // TsSym[$i, $j];
    BetaTensor[4, 2, 5] := BetaTensor[4, 2, 5] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 6] := BetaTensor[4, 2, 6] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 7] := BetaTensor[4, 2, 7] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 8] := BetaTensor[4, 2, 8] = C2FC2Gt[$i, k1]. Yuk[$a, k1, $j, True] // Expand // TsSym[$i, $j];
    BetaTensor[4, 2, 9] := BetaTensor[4, 2, 9] = C2FS2St[$i, k1]. Yuk[$a, k1, $j, True] // Expand // TsSym[$i, $j];
    BetaTensor[4, 2, 10] := BetaTensor[4, 2, 10] = C2FS2Ft[$i, k1]. Yuk[$a, k1, $j, True] // Expand // TsSym[$i, $j];
    BetaTensor[4, 2, 11] := BetaTensor[4, 2, 11] = {{0, 0}, {0, 0}};
(*y^3 terms*)
    BetaTensor[4, 2, 12] := BetaTensor[4, 2, 12] = Tdot[TfermTil[A1, $i, k1], Yuk[$a, k1, k2, True], YukTil[b1, k2, k3], TfG2t[A1, k3, k4], Yuk[b1, k4, $j]] // TsSym[$i, $j];
    BetaTensor[4, 2, 13] := BetaTensor[4, 2, 13] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 14] := BetaTensor[4, 2, 14] = Tdot[Yuk[b1, $i, k1], YukTil[$a, k1, k2, True], Yuk[b2, k2, $j]] C2S[b1, b2] // Expand;
    BetaTensor[4, 2, 15] := BetaTensor[4, 2, 15] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 16] := BetaTensor[4, 2, 16] = Tdot[Yuk[b1, $i, k1], C2F[k1, k2], YukTil[$a, k2, k3, True], Yuk[b1, k3, $j]] // TsSym[$i, $j];
    BetaTensor[4, 2, 17] := BetaTensor[4, 2, 17] = Tdot[C2Ft[$i, k1], Yuk[b1, k1, k2], YukTil[$a, k2, k3, True], Yuk[b1, k3, $j]] // TsSym[$i, $j];
    BetaTensor[4, 2, 18] := BetaTensor[4, 2, 18] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 19] := BetaTensor[4, 2, 19] = Yuk[$a, $i, k1, True].Y2FC2St[k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[4, 2, 20] := BetaTensor[4, 2, 20] = Yuk[$a, $i, k1, True].Y2FC2Ft[k1, $j] // Expand // TsSym[$i, $j];
    BetaTensor[4, 2, 21] := BetaTensor[4, 2, 21] = Tdot[Yuk[$a, $i, k1, True], Y2Ft[k1, k2], C2F[k2, $j]] // TsSym[$i, $j];
    BetaTensor[4, 2, 22] := BetaTensor[4, 2, 22] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 23] := BetaTensor[4, 2, 23] = Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[b3, k2, $j]] Lam[$a, b1, b2, b3, True] // Expand;
(*y^5 terms*)
    BetaTensor[4, 2, 24] := BetaTensor[4, 2, 24] = Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[$a, k2, k3, True], YukTil[b1, k3, k4], Yuk[b2, k4, $j]];
    BetaTensor[4, 2, 25] := BetaTensor[4, 2, 25] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 26] := BetaTensor[4, 2, 26] = Tdot[Yuk[b2, $i, k1], YukTil[b1, k1, k2], Yuk[$a, k2, k3, True], YukTil[b1, k3, k4], Yuk[b2, k4, $j]];
    BetaTensor[4, 2, 27] := BetaTensor[4, 2, 27] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 28] := BetaTensor[4, 2, 28] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 29] := BetaTensor[4, 2, 29] = Tdot[Yuk[b1, $i, k1], Y2Ft[k1, k2], YukTil[$a, k2, k3, True], Yuk[b1, k3, $j]] // TsSym[$i, $j];
    BetaTensor[4, 2, 30] := BetaTensor[4, 2, 30] = Tdot[Yuk[$a, $i, k1, True], Y2FY2Ft[k1, $j]] // TsSym[$i, $j];
    BetaTensor[4, 2, 31] := BetaTensor[4, 2, 31] = {{0, 0}, {0, 0}};
    BetaTensor[4, 2, 32] := BetaTensor[4, 2, 32] = Tdot[Yuk[b1, $i, k1], YukTil[$a, k1, k2, True], Yuk[b2, k2, $j]] Y2S[b1, b2] // Expand;
    BetaTensor[4, 2, 33] := BetaTensor[4, 2, 33] = Tdot[Yuk[$a, $i, k1, True], Y2FY2St[k1, $j]] // TsSym[$i, $j];

(*#######################################*)
(*--------Massive Quartic Tensors--------*)
(*#######################################*)
    BetaTensor[5, 0, 1] := BetaTensor[5, 0, 1] = - Global`\[Epsilon]/2 Lam[$a, $b, $c, $d, True];
    BetaTensor[5, 1, 1] := BetaTensor[5, 1, 1] = 0;
    BetaTensor[5, 1, 2] := BetaTensor[5, 1, 2] = C2S[$a, b1] Lam[b1, $b, $c, $d, True] // Expand // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 1, 3] := BetaTensor[5, 1, 3] = Lam[$a, $b, b1, b2, True] Lam[b1, b2, $c, $d, True] // Expand // TsSym[$b, $c, $d];
    BetaTensor[5, 1, 4] := BetaTensor[5, 1, 4] = Y2S[$a, b1] Lam[b1, $b, $c, $d, True] // Expand // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 1, 5] := BetaTensor[5, 1, 5] = Tr @ Tdot[Yuk[$a, k4, k1, True], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True]] // TsSym[$b, $c, $d];
(* y^0, lam^0 terms*)
    BetaTensor[5, 2, 1] := BetaTensor[5, 2, 1] = 0;
    BetaTensor[5, 2, 2] := BetaTensor[5, 2, 2] = 0;
    BetaTensor[5, 2, 3] := BetaTensor[5, 2, 3] = 0;
    BetaTensor[5, 2, 4] := BetaTensor[5, 2, 4] = 0;
    BetaTensor[5, 2, 5] := BetaTensor[5, 2, 5] = 0;
    BetaTensor[5, 2, 6] := BetaTensor[5, 2, 6] = Ttimes[FourGLam[$a, b1, $d, b2], Lam[b1, b2, $b, $c, True]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 7] := BetaTensor[5, 2, 7] = Ttimes[FourGLam[$a, $b, b1, b2], Lam[b1, b2, $c, $d, True]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 8] := BetaTensor[5, 2, 8] = 0;
    BetaTensor[5, 2, 9] := BetaTensor[5, 2, 9] = Ttimes[C2S[$a, b1], C2S[b1, b2], Lam[b2, $b, $c, $d, True] ]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 2, 10] := BetaTensor[5, 2, 10] = Ttimes[C2SC2G[$a, b1], Lam[b1, $b, $c, $d, True] ]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 2, 11] := BetaTensor[5, 2, 11] = Ttimes[C2SS2S[$a, b1], Lam[b1, $b, $c, $d, True] ]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 2, 12] := BetaTensor[5, 2, 12] = Ttimes[C2SS2F[$a, b1], Lam[b1, $b, $c, $d, True] ]  // TsSym4[$a, $b, $c, $d];
(* y^0, lam^{n>0} terms*)
    BetaTensor[5, 2, 13] := BetaTensor[5, 2, 13] = Ttimes[Tscal[A1, $a, b1], Lam[b1, b2, b3, b4], TsG2[A1, $b, b2], Lam[b3, b4, $c, $d, True]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 14] := BetaTensor[5, 2, 14] = Ttimes[Lam[$a, $b, b1, b2, True], C2S[b2, b3], Lam[b3, b1, $c, $d, True]] // TsSym[$b, $c, $d];
    BetaTensor[5, 2, 15] := BetaTensor[5, 2, 15] = Ttimes[C2S[$a, b1], Lam[b1, $b, b2, b3, True], Lam[b2, b3, $c, $d, True]]  // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 16] := BetaTensor[5, 2, 16] = Ttimes[Lam2[$a, b1], Lam[b1,$b, $c, $d, True]]  // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 2, 17] := BetaTensor[5, 2, 17] = Ttimes[Lam[$a, b1, b2, b3, True], Lam[$b, b4, b2, b3, True], Lam[b1, b4, $c, $d, True]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 18] := BetaTensor[5, 2, 18] = 0;
(*y^2 terms*)
    BetaTensor[5, 2, 19] := BetaTensor[5, 2, 19] = Ttimes[Tr @ Tdot[Yuk[$a, k4, k1, True], Tferm[A1, k1, k2], Tferm[A2, k2, k3], YukTil[$b, k3, k4, True]], TsG2[A1, $c, b1], TsG2[A2, b1, $d]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 20] := BetaTensor[5, 2, 20] = 0;
    BetaTensor[5, 2, 21] := BetaTensor[5, 2, 21] = Ttimes[Y2SC2F[$a, b1], Lam[b1, $b, $c, $d, True]] // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 2, 22] := BetaTensor[5, 2, 22] = 0;
    BetaTensor[5, 2, 23] := BetaTensor[5, 2, 23] = Ttimes[Lam[$a, $b, b1, b2, True], Y2S[b2, b3], Lam[b3, b1, $c, $d, True]] // TsSym[$b, $c, $d];
(*y^4 terms*)
    BetaTensor[5, 2, 24] := BetaTensor[5, 2, 24] = 0;
    BetaTensor[5, 2, 25] := BetaTensor[5, 2, 25] = Tr @ Tdot[Yuk[b1, k4, k1], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True]] C2S[$a, b1] // Expand // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 26] := BetaTensor[5, 2, 26] = Tr @ Tdot[Yuk[$a, k5, k1, True], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True], C2Ft[k4, k5]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 27] := BetaTensor[5, 2, 27] = Tr @ Tdot[Yuk[$a, k4, k1, True], YukTil[b1, k1, k2], Yuk[$b, k2, k3, True], YukTil[b2, k3, k4]] Lam[b1, b2, $c, $d, True] // Expand // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 28] := BetaTensor[5, 2, 28] = 0;
    BetaTensor[5, 2, 29] := BetaTensor[5, 2, 29] = Y4cS[$a, b1] Lam[b1, $b, $c, $d, True] // Expand // TsSym4[$a, $b, $c, $d];
    BetaTensor[5, 2, 30] := BetaTensor[5, 2, 30] = Y2SY2F[$a, b1] Lam[b1, $b, $c, $d, True] // Expand // TsSym4[$a, $b, $c, $d];
(*y^6 terms*)
    BetaTensor[5, 2, 31] := BetaTensor[5, 2, 31] = Tr @ Tdot[Yuk[b1, k6, k1], YukTil[$a, k1, k2, True], Yuk[b1, k2, k3], YukTil[$b, k3, k4, True], Yuk[$c, k4, k5, True], YukTil[$d, k5, k6, True]] // TsSym[$a, $b, $c, $d];
    BetaTensor[5, 2, 32] := BetaTensor[5, 2, 32] = Tr @ Tdot[Yuk[$a, k6, k1, True], YukTil[b1, k1, k2], Yuk[$b, k2, k3, True], YukTil[$c, k3, k4, True], Yuk[b1, k4, k5], YukTil[$d, k5, k6, True]] // TsSym[$b, $c, $d];
    BetaTensor[5, 2, 33] := BetaTensor[5, 2, 33] = Tr @ Tdot[Yuk[$a, k5, k1, True], YukTil[$b, k1, k2, True], Yuk[$c, k2, k3, True], YukTil[$d, k3, k4, True], Y2F[k4, k5]] // TsSym[$a, $b, $c, $d];


(*#################################################*)
(*-----------Anomalous dimension tensors-----------*)
(*#################################################*)
(* Fermion field anomalous dimension *)
    AnomalousTensor[1, 1, 1] := C2Ft[$i, $j];
    AnomalousTensor[1, 1, 2] := Y2F[$i, $j];

    AnomalousTensor[1, 2, 1] := Tdot[C2Ft[$i, k1], C2Ft[k1, $j] ];
    AnomalousTensor[1, 2, 2] := C2FC2Gt[$i, $j];
    AnomalousTensor[1, 2, 3] := C2FS2St[$i, $j];
    AnomalousTensor[1, 2, 4] := C2FS2Ft[$i, $j];
    AnomalousTensor[1, 2, 5] := Y2FC2S[$i, $j];
    AnomalousTensor[1, 2, 6] := Y2FC2F[$i, $j];
    AnomalousTensor[1, 2, 7] := Tdot[C2Ft[$i, k1], Y2F[k1, $j]];
    AnomalousTensor[1, 2, 8] := Y2FY2F[$i, $j];
    AnomalousTensor[1, 2, 9] := Y2FY2S[$i, $j];

(* Fermion field anomalous dimension *)
    AnomalousTensor[2, 1, 1] := C2S[$a, $b];
    AnomalousTensor[2, 1, 2] := Y2S[$a, $b];

    AnomalousTensor[2, 2, 1] := Ttimes[C2S[$a, b1], C2S[b1, $b]];
    AnomalousTensor[2, 2, 2] := C2SC2G[$a, $b];
    AnomalousTensor[2, 2, 3] := C2SS2S[$a, $b];
    AnomalousTensor[2, 2, 4] := C2SS2F[$a, $b];
    AnomalousTensor[2, 2, 5] := Lam2[$a, $b];
    AnomalousTensor[2, 2, 6] := Y2SC2F[$a, $b];
    AnomalousTensor[2, 2, 7] := Y4cS[$a, $b];
    AnomalousTensor[2, 2, 8] := Y2SY2F[$a, $b];


(*#####################################*)
(*-----------Upsilon tensors-----------*)
(*#####################################*)
(* Fermion upsilon tensors *)
    UpsilonTensor[1, 3, 1] := UpsilonTensor[1, 3, 1] = MinusHc[$i, $j]@ Tdot[TfermTil[A1, $i, k1], Y2F[k1, k2], Yuk[b1, k2, k3], TfG2[A1, k3, k4], YukTil[b1, k4, $j]];
    UpsilonTensor[1, 3, 2] := UpsilonTensor[1, 3, 2] = MinusHc[$i, $j]@ Tdot[Ttimes[Yuk[b1, $i, k1], C2S[b1, b2]], YukTil[b3, k1, k2], Yuk[b2, k2, k3], YukTil[b3, k3, $j]];
    UpsilonTensor[1, 3, 3] := UpsilonTensor[1, 3, 3] = MinusHc[$i, $j]@ Tdot[Yuk[b1, $i, k1], C2F[k1, k2], YukTil[b2, k2, k3], Yuk[b1, k3, k4], YukTil[b2, k4, $j]];
    UpsilonTensor[1, 3, 4] := UpsilonTensor[1, 3, 4] = MinusHc[$i, $j]@ Tdot[Yuk[b1, $i, k1], YukTil[b2, k1, k2], Yuk[b3, k2, k3], YukTil[b2, k3, k4], Yuk[b1, k4, k5], YukTil[b3, k5, $j]];
    UpsilonTensor[1, 3, 5] := UpsilonTensor[1, 3, 5] = MinusHc[$i, $j]@ Tdot[Yuk[b1, $i, k1], Y2Ft[k1, k2], YukTil[b2, k2, k3], Yuk[b1, k3, k4], YukTil[b2, k4, $j]];
    UpsilonTensor[1, 3, 6] := UpsilonTensor[1, 3, 6] = MinusHc[$i, $j]@ Tdot[Ttimes[Yuk[b1, $i, k1], Y2S[b1, b2]], YukTil[b3, k1, k2], Yuk[b2, k2, k3], YukTil[b3, k3, $j]];

    UpsilonTensor[2, 3, 1] := UpsilonTensor[2, 3, 1] = MinusTrans[$a, $b]@ Tr@ Tdot[Yuk[$a, k1, k2], C2F[k2, k3], YukTil[b1, k3, k4], Yuk[$b, k4, k5], YukTil[b1, k5, k1]];
    UpsilonTensor[2, 3, 2] := UpsilonTensor[2, 3, 2] = MinusTrans[$a, $b]@ Ttimes[Tr@ Tdot[Yuk[$a, k1, k2], YukTil[b1, k2, k3], Yuk[b2, k3, k4], YukTil[b3, k4, k1]], Lam[b1, b2, b3, $b]];
    UpsilonTensor[2, 3, 3] := UpsilonTensor[2, 3, 3] = MinusTrans[$a, $b]@ Tr@ Tdot[Yuk[$a, k1, k2], Y2Ft[k2, k3], YukTil[b1, k3, k4], Yuk[$b, k4, k5], YukTil[b1, k5, k1]];



];
