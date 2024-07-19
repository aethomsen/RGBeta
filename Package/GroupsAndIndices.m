(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageExport["del"]
PackageExport["delA2"]
PackageExport["delIndex"]
PackageExport["delS2"]
PackageExport["dSym"]
PackageExport["eps"]
PackageExport["fStruct"]
PackageExport["lcSymb"]
PackageExport["tGen"]
PackageExport["Dim"]
PackageExport["Casimir2"]
PackageExport["TraceNormalization"]

PackageExport["SO"]
PackageExport["Sp"]
PackageExport["SU"]
PackageExport["U1"]

PackageExport["adj"]
PackageExport["fund"]
PackageExport["A2"]
PackageExport["S2"]

PackageExport["Bar"]
PackageExport["DefineLieGroup"]
PackageExport["Matrix"]
PackageExport["SetReal"]
PackageExport["Tensor"]
PackageExport["Trans"]

PackageExport["$LieGroups"]

PackageScope["AdjTrFactor"]
PackageScope["PerformColorAlg"]
PackageScope["RefineGroupStructures"]
PackageScope["ReInitializeSymbols"]
PackageScope["ReInitializeCouplingBehavior"]
PackageScope["TStructure"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

del::usage =
	"del[rep, a, b] represents a Kronecker delta in \"rep space\" with indices a and b."
delA2::usage =
	"delA2[group, a, i, j] represents the singlet contraction between the A2 index a and the two fundamental indices i and j."
delIndex::usage =
	"del[rep, a, b] represents a Kronecker delta function in \"rep space\" with either a or b being an Integer fixing the value of the other summation index."
delS2::usage =
	"delS2[group, a, i, j] represents the singlet contraction between the S2 index a and the two fundamental indices i and j."
dSym::usage =
	"dSym[group, A, B, C] represents the symmetric 3-index tensor of the adjoint representation of SU(N) with indices A, B, and C."
eps::usage =
	"eps[rep, a, b] represents the 2-index antisymmetric tensor in \"rep space\" with indices a and b."
fStruct::usage =
	"fStruct[group, A, B, C] represents the structure constant of the group with indices A, B, and C."
lcSymb::usage =
	"Levi-Civita symbol"
tGen::usage =
	"tGen[rep, A, a, b] represents a group generator of the representation \"rep\" with adjoint index A. a and b are the two indices of rep. "

adj::usage = fund::usage = S2::usage = A2::usage =
	"Representation name defined defined by the Define(SO/Sp/SU)Group functions."
SO::usage = Sp::usage = SU::usage = U1::usage =
	"SO, Sp, SU, and U1 are used to specify different Lie groups."

Bar::usage =
	"Bar[x] rerpresents the conjugate of x. Used both for fields, couplings, and representations."
DefineLieGroup::usage =
	"DefineLieGroup[groupName, lieGroup[n]] sets groupName to be a lie group of type SU(n), SO(n), Sp(n), or U(1)^n and defines invariants for several common representations of said group."
Matrix::usage =
	"Matrix[x,...][i, j] represents the matrix product of couplings x,... with open indices i and j."
SetReal::usage =
	"SetReal[x1,...] makes x1,... behave as real parameters under complex conjugation (the Bar function)."
Tensor::usage =
	"Tensor[x][i,...] represents the tensor x with open indices i,..."
Trans::usage =
	"Trans[coupling] represents the transposed quantity of a coupling with two indices."

Casimir2::usage =
	"Casimir2[rep] sets the quadratic casimir of a given representation."
Dim::usage =
	"Dim[rep] sets the dimension of a given representation."
TraceNormalization::usage =
	"TraceNormalization[rep] sets the trace normalization of a given representation."


RefineGroupStructures::usage =
	"RefineGroupStructures[expr] decomposes generators of non-fundamental representations of the groups to the fundamental ones, whereby identities can be applied."
ReInitializeSymbols::usage =
	"ReInitializeSymbols[] is a function which when called flushes all previous definitions for symbol behaviour under implicit summation."
TStructure::usage =
	"TStructure[ind1, ind2,...][sparseArray] is used as a wrapper for a sparse array containing a particular tensor structure, while providing the indices of said struture."

(*#############################################*)
(*----------Generic tensor properties----------*)
(*#############################################*)

(*An overall function which sets up the properties of symbols wrt. implicit summation. It flushes all previous definitions.*)
ReInitializeSymbols[] :=
	Module[{},
		(*Dim[rep] is the dimension of a given representation*)
		Clear@ Dim;
		Dim[Bar[rep_]] := Dim[rep];

		(*Group constants for certain representations.*)
		Clear@ TraceNormalization;
		Clear@ Casimir2;

		(*del[rep, a, b] is the symbol for Kronecker delta \delta_{a,b} belong to the indices specified by the representation.*)
		Clear @ del;
		del[rep_, a_, b_] /; !OrderedQ@ {a, b}:= del[rep, b, a];
		(* del /: del[rep_, a___, x_, b___] del[rep_, c___, x_, d___] := del[rep, c, a, b, d]; *)
		del /: del[rep_, OrderlessPatternSequence[x_, a_]] del[rep_, OrderlessPatternSequence[x_, b_]] := del[rep, a, b];
		del /: Power[del[rep_, a_, b_], 2] := Dim[rep];
		del[rep_, a_, a_] = Dim[rep];
		del[Bar[rep_], a_, b_] = del[rep, a, b];

		(*del for taking specific index*)
		Clear @ delIndex;
		delIndex /: del[rep_, a___, x:Except[_Integer], b___] delIndex[rep_, c___, x:Except[_Integer], d___] = delIndex[rep, c, a, b, d];
		delIndex /: delIndex[rep_, a___, x:Except[_Integer], b___] delIndex[rep_, c___, x:Except[_Integer], d___] = delIndex[rep, c, a, b, d];
		delIndex /: Power[delIndex[rep_, OrderlessPatternSequence[Except[_Integer], _]], 2] = 1;
		delIndex[rep_, a_Integer, b_Integer] := KroneckerDelta[a, b];

		(*Default properties for 2-index anti-symmetric inavariants.*)
		Clear @ eps;
		(* eps /: del[rep_, a___, x_, b___] eps[rep_, c___, x_, d___] := eps[rep, c, a, b, d]; *)
		eps /: del[rep_, OrderlessPatternSequence[x_, a_]] eps[rep_, b___, x_, c___] := eps[rep, b, a, c];
		eps /: eps[rep_, a_, b_] eps[rep_, b_, c_] := -del[rep, a, c];
		eps /: eps[rep_, a_, b_] eps[rep_, c_, b_] := del[rep, a, c];
		eps /: eps[rep_, b_, a_] eps[rep_, b_, c_] := del[rep, a, c];
		eps /: Power[eps[rep_, a_, b_], 2] := Dim[rep];
		eps[rep_, a_, a_] := 0;

		(*Levi-Civita symbols*)
		Clear @ lcSymb;
		(* lcSymb /: del[rep_, a___, x_, b___] lcSymb[rep_, c___, x_, d___] := lcSymb[rep, c, a, b, d]; *)
		lcSymb /: del[rep_, OrderlessPatternSequence[x_, a_]] lcSymb[rep_, b___, x_, c___] := lcSymb[rep, b, a, c];
		lcSymb[rep_, a___, x_, b___, x_, c___] := 0;
		lcSymb /: lcSymb[rep_, x:OrderlessPatternSequence[k__, i___]] lcSymb[rep_, y:OrderlessPatternSequence[k__, j___]] :=
			Factorial @ Length @ List[k] * Signature@ List[j] *
			Signature[List@ x] Signature[List@ y] Signature[{k, i}] Signature[{k, j}] *
			Sum[Signature @ perm * Times @@ MapThread[ del[rep, #1, #2] &, {List[i], perm}], {perm, Permutations @ List[j] }];
		(*lcSymb /: lcSymb[rep_, i__] lcSymb[rep_, j__] := 	Signature@ List[j] *
			Sum[Signature @ perm Times @@ Thread @ del[rep, List[i], perm], {perm, Permutations @ List[j] }];*)
		lcSymb /: Power[lcSymb[rep_, i__], 2] = Factorial @ Dim @ rep;

		(*Invariant of S2 and fundamentals*)
		Clear @ delS2;
		delS2 /: del[group_[S2], a___, x_ , b___] delS2[group_, x_, i_, j_]=  delS2[group, a, b, i, j];
		delS2 /: del[group_[fund], k___, x_, l___] delS2[group_, a_, i___, x_, j___]=  delS2[group, a, i, k, l, j];
		delS2 /: delS2[group_, a_, i_ , j_] delS2[group_, b_, i_, j_] = del[group@ S2, a, b];
		delS2 /: delS2[group_, a_, i_ , j_] delS2[group_, b_, j_, i_] = del[group@ S2, a, b];
		delS2 /: Power[delS2[group_, a_, i_ , j_], 2] = Dim[group@ S2];

		(*Invariant of A2 and fundamentals*)
		Clear @ delA2;
		delA2 /: del[group_[A2], a___, x_ , b___] delA2[group_, x_, i_, j_]=  delA2[group, a, b, i, j];
		delA2 /: del[group_[fund], k___, x_, l___] delA2[group_, a_, i___, x_, j___]=  delA2[group, a, i, k, l, j];
		delA2 /: delA2[group_, a_, i_ , j_] delA2[group_, b_, i_, j_] = del[group@ A2, a, b];
		delA2 /: delA2[group_, a_, i_ , j_] delA2[group_, b_, j_, i_] = - del[group@ A2, a, b];
		delA2 /: Power[delA2[group_, a_, i_ , j_], 2] = Dim[group@ A2];

		(*Default group generator properties.*)
		Clear @ tGen;
		tGen /: del[rep_, OrderlessPatternSequence[x_, a_]] tGen[rep_, A_, b___, x_, c___] := tGen[rep, A, b, a, c];
		tGen /: del[group_[adj], OrderlessPatternSequence[X_, A_]] tGen[group_[rep_], X_, c__] := tGen[group[rep], A, c];
		tGen /: tGen[rep_, A_, a_, b_] tGen[rep_, A_, b_, c_] := Casimir2[rep] del[rep, a, c];
		tGen /: tGen[rep_, A_, a_, b_] tGen[rep_, B_, b_, a_] := TraceNormalization[rep] del[Head[rep][adj], A, B];
		(*For the real representation we need*)
		tGen /: Power[tGen[rep_, A_, a_, b_], 2] := -Casimir2[rep] * Dim[rep];
		tGen /: tGen[rep_, A_, a_, b_] tGen[rep_, A_, c_, b_] := - Casimir2[rep] del[rep, a, c];
		tGen /: tGen[rep_, A_, b_, a_] tGen[rep_, A_, b_, c_] := - Casimir2[rep] del[rep, a, c];
		tGen[rep_, A_, a_, a_] = 0;
		tGen[Bar[rep_], A_, a_, b_ ] = - tGen[rep, A, b, a];

		(*Default structure constant properties.*)
		Clear @ fStruct;
		fStruct /: del[group_[adj], OrderlessPatternSequence[X_, A_]] fStruct[group_, B___, X_, C___] := fStruct[group, B, A, C];
		fStruct /: fStruct[group_, x:OrderlessPatternSequence[A_, B_, C_] ] fStruct[group_, y:OrderlessPatternSequence[D_, B_, C_] ] :=
			Signature[List@ x] Signature[List@ y] Signature[{A, B, C}] Signature[{D, B, C}] Casimir2[group@ adj] del[group@ adj, A, D];
		fStruct /: Power[fStruct[group_, __], 2] := Dim[group@ adj] Casimir2[group@ adj];
		fStruct /: fStruct[group_, a_, b_, c_] fStruct[group_, p1:OrderlessPatternSequence[a_, b_, c_]] := 
			Signature @ {a, b, c} * Signature @ {p1} * Dim[group@ adj] Casimir2[group@ adj];
		fStruct[group_, a___, x_, b___, x_, c___] := 0;

		(*Symmetric tensor*)
		Clear @ dSym;
		dSym[group_, OrderlessPatternSequence[A_, A_, _]]:= 0;
		dSym/: del[group_@ adj, OrderlessPatternSequence[A_, X_]] dSym[group_, OrderlessPatternSequence[X_, B_, C_]] := dSym[group, A, B, C];
		dSym/: fStruct[group_, OrderlessPatternSequence[A_, X_, Y_]] dSym[group_, OrderlessPatternSequence[X_, Y_, B_]] := 0; 

		(*Two-index delta*)
		Clear @ twoIndexRepDelta;
		twoIndexRepDelta[rep_, a_, b_] = 0;
	];

(*Function applying the group *)
	replace = {
		(* tGen[group_[adj], A_, B_, C_] :> - Module[{a, b, c}, 2 AntiSym[A, B][tGen[group@ fund, A, a, b] tGen[group@ fund, B, b, c]] tGen[group@ fund, C , c, a]
			] / TraceNormalization @ group @ fund, *)
		tGen[group_[adj], A_, B_, C_]:> - I fStruct[group, A, B, C],
		del[group_[S2], a_ ,b_] :> twoIndexRepDelta[group @ S2, a ,b],
		del[group_[A2], a_ ,b_] :> twoIndexRepDelta[group @ A2, a ,b],
		tGen[group_[S2], A_, a_, b_] -> 2 Sym[a@1, a@2] @ Sym[b@1, b@2][tGen[group@ fund, A, a@1, b@1] del[group@ fund, a@2, b@2]],
		delS2[group_, a_, i_, j_] -> Sym[a@1, a@2][ del[group@ fund, a@1, i] del[group@ fund, a@2, j]],
		tGen[group_[A2], A_, a_, b_] -> 2 AntiSym[a@1, a@2] @ AntiSym[b@1, b@2][tGen[group@ fund, A, a@1, b@1] del[group@ fund, a@2, b@2]],
		delA2[group_, a_, i_, j_] -> AntiSym[a@1, a@2][ del[group@ fund, a@1, i] del[group@ fund, a@2, j]]
		(*del[group_[A2], a_ ,b_] :> AntiSym[a @ 1, a @ 2 ][del[group@fund, a@1, b@1] del[group@fund, a@2, b@2] ]*)
	};
RefineGroupStructures[expr_] := expr /. replace // PerformColorAlg // EvalTrReductionFactors // Expand;

(*Adds case to the system built in function Tr and Dot, to deal with substituting couplings for 0.*)
	Unprotect[Tr, Dot];
		Dot[___, 0, ___] = 0;
		Dot[a___, x_? NumericQ b_, c___] := x Dot[a, b, c];
		Tr[0] = 0;
		Tr[x_? NumericQ a_ ] = x Tr@ a;
	Protect[Tr, Dot];

(*An overall function which sets up the properties of couplings in tensor and matrix contractions etc. It flushes all previous definitions.*)
ReInitializeCouplingBehavior[]:= Module[{},
	Clear[Bar, Trans, Matrix, Tensor];

	(*Formating*)
		Format[Trans[x_], StandardForm] := HoldForm[x^Global`T];
		Format[Bar[x_], StandardForm] := OverBar @ x;
		Format[Trans[Bar[x_]], StandardForm] := x^Style[Global`\[Dagger], Bold, 12];
		Format[Matrix[x__][h1_[i1_] ], StandardForm] := Subscript[Dot[x], i1];
		Format[Matrix[x__][h1_[i1_], h2_[i2_]], StandardForm] := Subscript[Dot[x], i1, i2];
		Format[Tensor[x_][inds__], StandardForm] := Subscript[x, Sequence@@ {inds}[[;;, 1]]];

	(*Complex conjugation and transposition of matrices*)
		Bar[Bar[x_]] = x;
		Bar[0] = 0;
		Bar[a:(_List| _Times| _Plus)] := Bar /@ a;
		Bar[a_] /; NumericQ[a] := Conjugate @ a;

		Trans[Trans[x_]] := x;
		Trans[0] = 0;
		Trans[x_? NumericQ a_]:= x Trans@ a;
		Trans[a_Plus]:= Trans /@ a;
		Trans[a_List] /; VectorQ[a] := a;
		Trans[a_List] /; MatrixQ[a] := Transpose @ a;

	(*Matrix head for symbolic matrix manipulations*)
		Matrix /: del[ind_, a___, x_, b___] Matrix[m__][c___, ind_[x_], d___] := Matrix[m][c, ind[a, b], d];
		Matrix /: Matrix[m1__][a_] Matrix[m2__][a_] := Matrix[Sequence @@ Reverse[Trans /@ List@m1], m2][];
		Matrix /: Matrix[m1__][b_] Matrix[m2__][b_, c_] := Matrix[Sequence @@ Reverse[Trans /@ List@m2], m1][c];
		Matrix /: Matrix[m1__][a_, b_] Matrix[m2__][b_] := Matrix[m1, m2][a];
		Matrix /: Matrix[m1__][a_, b_] Matrix[m2__][b_, c_] := Matrix[m1, m2][a, c];
		Matrix /: Matrix[m1__][a_, b_] Matrix[m2__][c_, b_] := Matrix[m1, Sequence @@ Reverse[Trans /@ List@m2]][a, c];
		Matrix /: Matrix[m1__][b_, a_] Matrix[m2__][b_, c_] := Matrix[Sequence @@ Reverse[Trans /@ List@m1], m2][a, c];
		Matrix /: Power[Matrix[m1__][a_, b_], 2] := Matrix[m1, Sequence @@ Reverse[Trans /@ List@m1]][a, a];
		Matrix /: Power[Matrix[m1__][a_], 2] := Matrix[m1, Sequence @@ Reverse[Trans /@ List@m1]][a, a];
		Matrix[a_List, b_List] = Matrix[Dot[a, b]];
		Matrix[m_][] = m;
		Matrix[m_][Null] = m;
		Matrix[___, 0, ___][___] = 0;
		Matrix[a___, x_?NumericQ b_, c___][ind__]:= x Matrix[a, b, c]@ ind;
		Matrix[m__][a_, a_] :=
			Block[{matrices, permutations},
				matrices = List @ m;
				permutations = NestList[RotateLeft, matrices, Length@ matrices];
				permutations = Join[permutations, Map[Trans, Reverse[permutations, 2], {2}] ];
				matrices = permutations[[Ordering[permutations, 1][[1]] ]];
				Tr @ Dot[Sequence @@ matrices]
			];
		Matrix[m__][] :=
			Block[{forms},
				forms = {Dot[m], Trans /@ Reverse @ Dot[m]};
				Sort[forms][[1]]
			];

	(*Tensor head for tensor coupling contractions*)
		Tensor /: del[ind_, a___, x_, b___] Tensor[t_][c___, ind_[x_], d___] := Tensor[t][c, ind[a, b], d];
];

SetReal[x_] := (Bar @ x = x;);
SetReal[x_, y__] := (SetReal @ x; SetReal @ y );

(* Tensor structure head used for all coupling contraction*)
TStructure /: TStructure[ind1__][ar1_] TStructure[ind2__][ar2_] = TStructure[ind1, ind2][TensorProduct[ar1, ar2]];
TStructure /: TStructure[ind__][ar_]^2 := TStructure[ind, ind]@ TensorProduct[ar, ar]
TStructure /: TStructure[ind__][ar1_] + TStructure[ind__][ar2_] = TStructure[ind][TensorProduct[ar1 + ar2]];
TStructure[][expr_] = expr;
TStructure[___][0] = 0;
(* Duplicate indices are contracted as per dummy index convention *)
TStructure[ind__][ar_] /; ! DuplicateFreeQ@ List@ ind :=
	Module[{cont, indices, temp},
		indices = List@ ind;
		temp = ar;
		While[!DuplicateFreeQ@ indices,
			(* For scalars the contraction is done via the $scalarContraction *)
			cont = First@ Cases[PositionIndex@ indices, _?(Length@ # === 2 &)];
			Switch[Head@ indices[[cont[[1]] ]]
			, $scalar,
				temp = TensorContract[TensorProduct[temp, $scalarContraction],
				Transpose@ {cont, Length@ indices + {1, 2}}];
			, $fermion | $gauge,
				temp = TensorContract[temp, {cont}];
			];
			indices = Delete[indices, List /@ cont];
		];
		(* The entries of the resulting array are expanded to allow the dummyindices to contract. *)
		If[Head@temp === SparseArray,
			temp = SparseArray[ArrayRules@ temp // Expand, Dimensions@ temp];
		,
			temp = Expand@ temp;
		];
		TStructure[Sequence @@ indices][temp]
	];

(*###########################################*)
(*----------Gauge group definitions----------*)
(*###########################################*)
DefineLieGroup::unkown = "`1` is not a Lie group or has not been implemented.";
DefineLieGroup::integerU1 = "U(1) groups should have a positive integer power";
DefineLieGroup[groupName_Symbol, U1] = DefineLieGroup[groupName, U1 @ 1 ];
DefineLieGroup[groupName_Symbol, lieGroup_Symbol[n_Integer|n_Symbol] ] :=
	Module[{},
		Switch[lieGroup
		,SO,
			DefineSOGroup[groupName, n];
		,Sp,
			If[IntegerQ @n && !EvenQ @ n,
				Message[DefineLieGroup::unkown, lieGroup[n] ];
				Return @ $Failed;
			];
			DefineSpGroup[groupName, n];
		,SU,
			DefineSUGroup[groupName, n];
		,U1,
			If[!IntegerQ @ n ||  n <= 0,
				Message[DefineLieGroup::integerU1];
				Return @ $Failed;
			];
			DefineU1Group[groupName, n];
		,_,
			Message[DefineLieGroup::unkown, lieGroup[n] ];
			Return @ $Failed;
		];
		$LieGroups@ groupName= lieGroup @ n;
	];

(*Initialization for an SO(n) gauge group.*)
DefineSOGroup[group_Symbol, n_Integer|n_Symbol] :=
	Module[{},
		(*Fundamental*)
		Dim[group[fund]] = n;
		TraceNormalization[group[fund]] = 1/2;
		Casimir2[group[fund]] = (n - 1) / 4;
		(*Fierz identitiy*)
		tGen /: tGen[group[fund], A_, a_, b_] tGen[group[fund], A_, c_, d_] = TraceNormalization[group[fund]] / 2 *
			(del[group[fund], a, d] del[group[fund], c, b] - del[group[fund], a, c] del[group[fund], b, d]);

		(*Symmetric traceless*)
		Dim[group[S2]] = (n - 1) (n + 2) / 2 ;
		TraceNormalization[group[S2]] = (n + 2) / 2;
		Casimir2[group[S2]] =  n / 2;
		delS2 /: delS2[group, a_, i_ , j_] delS2[group, a_, k_, l_] = Sym[k, l][del[group@ fund, i, k] del[group@ fund, j, l]]
			- del[group@ fund, i, j] del[group@ fund, k, l] / n;
		twoIndexRepDelta[group @ S2, a_, b_] = Sym[a @ 1, a @ 2 ][del[group@ fund, a@ 1, b@ 1] del[group@ fund, a@ 2, b@ 2] ] -
			del[group@ fund, a@ 1, a@ 2] del[group@ fund, b@ 1, b@ 2] / n;

		(*Adjoint*)
		Dim[group[adj]] = n (n - 1) / 2;
		TraceNormalization[group[adj]] = (n - 2) / 2;
		Casimir2[group[adj]] = (n - 2) / 2;
	];

(*Initialization for an Sp(n) gauge group.*)
DefineSpGroup[group_Symbol, n_Integer|n_Symbol] :=
	Module[{},
		(*Fundamental*)
		Dim[group[fund]] = n;
		TraceNormalization[group[fund]] = 1/2;
		Casimir2[group[fund]] = (n + 1) / 4;
		(*Fierz identitiy*)
		tGen /: tGen[group[fund], A_, a_, b_] tGen[group[fund], A_, c_, d_] = TraceNormalization[group[fund]] / 2 *
			(del[group[fund], a, d] del[group[fund], c, b] - eps[group[fund], a, c] eps[group[fund], b, d]);

		(*Antisymmetric "traceless"*)
		Dim[group[A2]] = (n - 2) (n + 1) / 2 ;
		TraceNormalization[group[A2]] = (n - 2) / 2;
		Casimir2[group[A2]] =  n / 2;
		delA2 /: delA2[group, a_, i_ , j_] delA2[group, a_, k_, l_] = AntiSym[k, l][del[group@ fund, i, k] del[group@ fund, j, l]] -
			eps[group@ fund, i, j] eps[group@ fund, k, l] / n;
		twoIndexRepDelta[group @ A2, a_, b_] = AntiSym[a @ 1, a @ 2 ][del[group@ fund, a@ 1, b@ 1] del[group@ fund, a@ 2, b@ 2] ] -
			eps[group@ fund, a@ 1, a@ 2] eps[group@ fund, b@ 1, b@ 2] / n;

		(*Adjoint*)
		Dim[group[adj]] = n (n + 1) / 2;
		TraceNormalization[group[adj]] = (n + 2) / 2;
		Casimir2[group[adj]] = (n + 2) / 2;
	];

(*Initialization for an SU(n) Lie group.*)
DefineSUGroup[group_Symbol, n_Integer|n_Symbol] :=
	Module[{},
		(*Fundamental*)
		Dim[group[fund]] = n;
		TraceNormalization[group[fund]] = 1/2;
		Casimir2[group[fund]] = (n^2 - 1) / (2 n);
		(*Fierz identitiy*)
		tGen /: tGen[group[fund], A_, a_, b_] tGen[group[fund], A_, c_, d_] = TraceNormalization[group[fund]] *
			(del[group[fund], a, d] del[group[fund], c, b] - del[group[fund], a, b] del[group[fund], c, d] / n);

		(*Symmetric*)
		Dim[group[S2]] = n (n +1) /2;
		TraceNormalization[group[S2]] = (n + 2) / 2;
		Casimir2[group[S2]] =  (n - 1) (n + 2) / n;
		delS2 /: delS2[group, a_, i_ , j_] delS2[group, a_, k_, l_] = Sym[k, l][del[group@ fund, i, k] del[group@ fund, j, l]];
		twoIndexRepDelta[group @ S2, a_, b_] = Sym[a @ 1, a @ 2 ][del[group@fund, a@1, b@1] del[group@fund, a@2, b@2] ];

		(*Anti-Symmetric*)
		Dim[group[A2]] = n (n -1) /2;
		TraceNormalization[group[A2]] = (n - 2) / 2;
		Casimir2[group[A2]] =  (n - 2) (n + 1) / n;
		delA2 /: delA2[group, a_, i_ , j_] delA2[group, a_, k_, l_] = AntiSym[k, l][del[group@ fund, i, k] del[group@ fund, j, l]];
		twoIndexRepDelta[group @ A2, a_, b_] = AntiSym[a @ 1, a @ 2 ][del[group@fund, a@1, b@1] del[group@fund, a@2, b@2] ];
		(*del/: del[group @ A2, a_, b_] del[group @ fund, a_[[1]], b_]*)

		(*Symmetric tensor*)
		dSym/: dSym[group, OrderlessPatternSequence[A_, X_, Y_]] dSym[group, OrderlessPatternSequence[B_, X_, Y_]] = AdjTrFactor[n, 4]* del[group@ adj, A, B];
		dSym/: Power[dSym[group, __], 2] = (n^2 - 1)* AdjTrFactor[n, 4];

		(*Adjoint*)
		Dim[group[adj]] = n^2 - 1;
		TraceNormalization[group[adj]] = n;
		Casimir2[group[adj]] = n;
	];

(*Initialization for a U(1) gauge group.*)
DefineU1Group[group_Symbol, power_Integer:1] :=
	Module[{},
		Switch[power
		,1,
			del[group[_], ___] = 1;
			Dim[group[_]] = 1;
			tGen[group[x_], ___] = x;
		,n_ /; n > 1,
			tGen[group[x_List], A_, a_, b_] = Matrix[x][group[adj][A]];
			Dim[group[adj]] = power;
			del[group[_List], ___] = 1;
			Dim[group[_List]] = 1;
		];
		fStruct[group, __] = 0;
	];


(*##########################################*)
(*----------Advanced color algebra----------*)
(*##########################################*)
(* SU(N) algebra for the adjoint representation is evaluated using the *)


PerformColorAlg[expr_] := 
	Module[{out = Expand @ expr, cgs, cgGr},
		If[Head @ out === Plus, PerformColorAlg /@ out // Return; ];
		out = If[Head @ out === Times, List @@ out, {out} ];

		{cgs, out} = SelectAndDelteCases[out, _fStruct | _dSym];
		cgs = GatherBy[cgs, First];
		Times @@ out* Product[
				If[Length@ cgGr > 2, AdjContraction[Times @@ cgGr], Times @@ cgGr]
			, {cgGr, cgs}]
	];

AdjContraction @ prod:Times[prefact___, cgs:Longest[(_fStruct | _dSym) ..]] /; Length@{cgs} > 2 := 
	Module[{out},
		out = EliminateShortestCycle[Times@ cgs];
		If[FreeQ[out, AdjTrace], Return@ prod;];
		AdjContraction[Times@prefact* (out /. x_AdjTrace :> EvaluateAdjTrace@x) // Expand]	
	];
AdjContraction @ sum_Plus := AdjContraction /@ sum;
AdjContraction @ x_ := x;

EliminateShortestCycle[prod_Times] := Replace[prod, {
    	(*3-cycles*)
		Times[rest___, cyc:OrderlessPatternSequence[_[_, OrderlessPatternSequence[_, x1_, x2_]], 
	   		_[_, OrderlessPatternSequence[_, x2_, x3_]], _[_, OrderlessPatternSequence[_, x3_, x1_]] ]] :>
     	Times @ rest * CycleToTrace3 @ cyc
    ,
		(*4-cycles*)
    	Times[rest___, cyc:OrderlessPatternSequence[_[_, OrderlessPatternSequence[_, x1_, x2_]], 
			_[_, OrderlessPatternSequence[_, x2_, x3_]], _[_, OrderlessPatternSequence[_, x3_, x4_]], 
			_[_, OrderlessPatternSequence[_, x4_, x1_]] ]] :>
     	Times @ rest * CycleToTrace4 @ cyc
	,
		(*5-cycles*)
    	Times[rest___, cyc:OrderlessPatternSequence[_[_, OrderlessPatternSequence[_, x1_, x2_]], 
			_[_, OrderlessPatternSequence[_, x2_, x3_]], _[_, OrderlessPatternSequence[_, x3_, x4_]], 
			_[_, OrderlessPatternSequence[_, x4_, x5_]], _[_, OrderlessPatternSequence[_, x5_, x1_]] ]] :>
     	Times @ rest * CommuteTrace5 @ Times@ cyc * AdjTrace[]
    } ];

CycleToTrace3 @ OrderlessPatternSequence[h1_[gr_, p1:OrderlessPatternSequence[a1_, x1_, x2_]],
		h2_[gr_, p2 : OrderlessPatternSequence[a2_, x2_, x3_]], 
		h3_[gr_, p3 : OrderlessPatternSequence[a3_, x3_, x1_]] ] := 
 	Module[{pref},
		pref = If[h1 === fStruct, Signature @ {p1} * Signature @ {a1, x1, x2}, 1]* 
			If[h2 === fStruct, Signature @ {p2} * Signature @ {a2, x2, x3}, 1]*
			If[h3 === fStruct, Signature @ {p3} * Signature @ {a3, x3, x1}, 1]*
			Power[I, Count[{h1, h2, h3}, fStruct]];
		pref* AdjTrace[gr, {h1, h2, h3} /. {dSym -> DMat, fStruct -> FMat}, {a1, a2, a3} ]
	];

CycleToTrace4 @ OrderlessPatternSequence[h1_[gr_, p1 : OrderlessPatternSequence[a1_, x1_, x2_]], 
		h2_[gr_, p2 : OrderlessPatternSequence[a2_, x2_, x3_]], 
		h3_[gr_, p3 : OrderlessPatternSequence[a3_, x3_, x4_]], 
		h4_[gr_, p4 : OrderlessPatternSequence[a4_, x4_, x1_]] ] := 
	Module[{pref},
		pref = 
		If[h1 === fStruct, Signature @ {p1} * Signature @ {a1, x1, x2}, 1] * 
			If[h2 === fStruct, Signature@ {p2}* Signature@ {a2, x2, x3}, 1] *
				If[h3 === fStruct, Signature@ {p3}* Signature@ {a3, x3, x4}, 1] *
				If[h4 === fStruct, Signature@ {p4}* Signature@ {a4, x4, x1}, 1] *
				Power[I, Count[{h1, h2, h3, h4}, fStruct]];
		pref* AdjTrace[gr, {h1, h2, h3, h4} /. {dSym -> DMat, fStruct -> FMat}, {a1, a2, a3, a4} ]
	];

(* There is only one possible full contraction of 10 tensors with no 3- or 4-cycle. Any commutation will reduce it. *)
 CommuteTrace5 @ Times[fStruct[gr_, p1:OrderlessPatternSequence[a_, b_, e_]], 
 		fStruct[gr_, p2:OrderlessPatternSequence[c_, d_, e_]], rest__]:= 
	Signature @ {a, b, e} * Signature @ {c, d, e} * Signature @ {p1} * Signature @ {p2} * Times @ rest *
		(-fStruct[gr, c, b, e] fStruct[gr, a, e, d] - fStruct[gr, d, b, e] fStruct[gr, a, c, e]);

CommuteTrace5 @ Times[fStruct[gr_, p1:OrderlessPatternSequence[a_, b_, e_]], 
 		dSym[gr_, OrderlessPatternSequence[c_, d_, e_]], rest__]:= 
	Signature @ {a, b, e} * Signature @ {p1} * Times @ rest *
		(-fStruct[gr, a, c, e] dSym[gr, b, d, e] - fStruct[gr, a, d, e] fStruct[gr, b, c, e]);

CommuteTrace5 @ Times[dSym[gr_, OrderlessPatternSequence[b_, c_, e_]], dSym[gr_, OrderlessPatternSequence[a_, d_, e_]], rest__]:= 
	Module[{n = First@ $LieGroups@ gr},
		Times @ rest * (
			dSym[gr, a, c, e] dSym[gr, b, d, e] - fStruct[gr, a, b, e] fStruct[gr, c, d, e] + 
			2/n (del[gr@ adj, a, c] del[gr@ adj, b, d] - del[gr@ adj, a, d] del[gr@ adj, b, c]) 
		)
	];

CanonizeAdjTrace @ AdjTrace[gr_, types_, inds_] := 
	Block[{can, canTransp, ind, indTransp, ord, rots = Length @ types - 1},
		can = NestList[RotateLeft, types, rots];
		ind = NestList[RotateLeft, inds, rots];
		canTransp = Reverse /@ can; 
		indTransp = Reverse /@ ind;
		ord = First @ Ordering @ can;
		can = can[[ord]]; ind = ind[[ord]];
		ord = First @ Ordering @ canTransp;
		canTransp = canTransp[[ord]]; indTransp = indTransp[[ord]];
		If[OrderedQ @ {can, canTransp},
			AdjTrace[gr, can, ind]
		,
			Power[-1, Count[types, FMat]] * AdjTrace[gr, canTransp, indTransp]
		]
	];

EvaluateAdjTrace::unkwn = "The cadjoint trace has not been implemented: `1`";
EvaluateAdjTrace @ adjTr:AdjTrace[gr_, types_, _] := 
	Module[{out, n = First @ $LieGroups @ gr, e},
		out = CanonizeAdjTrace @ adjTr;
		out /. Switch[Length @ types
			, 3, {	
				AdjTrace[gr, {FMat, FMat, FMat}, {a_, b_, c_}] :> 
					I* n/ 2* fStruct[gr, a, b, c],
				AdjTrace[gr, {DMat, FMat, FMat}, {a_, b_, c_}] :> 
					n/ 2* dSym[gr, a, b, c],
				AdjTrace[gr, {DMat, DMat, FMat}, {a_, b_, c_}] :> 
					I*  AdjTrFactor[n, 4]/2* fStruct[gr, a, b, c],
				AdjTrace[gr, {DMat, DMat, DMat}, {a_, b_, c_}] :> 
					AdjTrFactor[n, 12]/2* dSym[gr, a, b, c]
				}
			, 4, {
				AdjTrace[gr, {FMat, FMat, FMat, FMat}, {a_, b_, c_, d_}] :> (
					del[gr@ adj, a, d]  del[gr@ adj, b, c] +
					(del[gr@ adj, a, b]  del[gr@ adj, c, d] + del[gr@ adj, a, c]  del[gr@ adj, b, d])/2 +
					n/4* ( fStruct[gr, a, d, e]  fStruct[gr, b, c, e] + dSym[gr, a, d, e]  dSym[gr, b, c, e])
					),
				AdjTrace[gr, {DMat, FMat, FMat, FMat}, {d_, a_, b_, c_}] :> 
					1/4  I  n (dSym[gr, a, d, e]  fStruct[gr, b, c, e] - fStruct[gr, a, d, e]  dSym[gr, b, c, e] ),
				AdjTrace[gr, {DMat, DMat, FMat, FMat}, {c_, d_, a_, b_}] :> (
					(del[gr@ adj, a, b]  del[gr@ adj, c, d] - del[gr@ adj, a, c]  del[gr@ adj, b, d])/2 +
					AdjTrFactor[n, 8]/4  fStruct[gr, a, d, e]  fStruct[gr, b, c, e] +
					n/4   dSym[gr, a, d, e]  dSym[gr, b, c, e]),
				AdjTrace[gr, {DMat, FMat, DMat, FMat}, {d_, a_, b_, c_}] :> (
					-(del[gr@ adj, a, b]  del[gr@ adj, c, d] - del[gr@ adj, a, c]  del[gr@ adj, b, d])/2 +
					n/4  fStruct[gr, a, d, e]  fStruct[gr, b, c, e] +
					n/4   dSym[gr, a, d, e]  dSym[gr, b, c, e]),
				AdjTrace[gr, {DMat, DMat, DMat, FMat}, {b_, c_, d_, a_}] :> (
					2 I/n * fStruct[gr, a, d, e]  dSym[gr, b, c, e] +
					I * AdjTrFactor[n, 8]/4  fStruct[gr, a, b, e]  dSym[gr, c, d, e] +
					I/4 * n * dSym[gr, a, b, e]  fStruct[gr, c, d, e]),
				AdjTrace[gr, {DMat, DMat, DMat, DMat}, {a_, b_, c_, d_}] :> (
					AdjTrFactor[n, 4]/n (del[gr@ adj, a, b]  del[gr@ adj, c, d] + del[gr@ adj, a, d]  del[gr@ adj, b, c]) +
					AdjTrFactor[n, 16]/4 (dSym[gr, a, b, e]  dSym[gr, c, d, e] + dSym[gr, a, d, e]  dSym[gr, b, c, e]) -
					n/4  dSym[gr, a, c, e]  dSym[gr, b, d, e])
				}
			, _,
				Message[EvaluateAdjTrace::unkwn, out]; Abort[];
		]
	];	
EvaluateAdjTrace@ AdjTrace[]= 1;

(*Introduced to minimize the no. of terms at intermediate steps when N is kept symbolic*)
EvalTrReductionFactors@ expr_ := expr /. AdjTrFactor[n_, a_] :> (n^2 - a)/n; 


