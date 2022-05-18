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
ScalarAnomalousTensors::uasge =
	"ScalarAnomalousTensors[field, loop] evaluates all tensor contractions of the general anomalous dimension tensor and the field projector."
ScalarMassiveTensors::usage =
	"ScalarMassiveTensors[coupling, loop] is a function that computes all the tensor contractions used in the the trillinear and scalar mass beta functions at the given loop order."
YukawaTensors::usage =
	"YukawaTensors[coupling, loop] is a function that computes all the tensor contractions used the general yukawa beta function at the given loop order."

(*########################################*)
(*----------All kinds of tensors----------*)
(*########################################*)

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
	Module[{diagrams = {1, 5, 33}[[loop + 1]], n= 0},
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
