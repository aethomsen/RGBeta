Begin["GroupsAndIndices`"]

(*#############################################*)
(*----------Generic tensor properties----------*)
(*#############################################*)
CasimirSig[group_, a_, b_] := 
	Block[{common, ind1, ind2},
		common = Intersection[a, b];
		ind1 = Complement[a, common][[1]];
		ind2 = Complement[b, common][[1]];
		Signature[DeleteCases[a, ind1]] Signature[DeleteCases[b, ind2]] del[group[adj], ind1, ind2]
	];

(*An overall function which sets up the properties of symbols wrt. implicit summation. It flushes all previous definitions.*)
ReInitializeSymbols[] :=
	Block[{},
		(*Dim[rep] is the dimension of a given representation*)
		Clear @ Dim;
		Dim[Bar[rep_]] := Dim[rep];
		
		(*Group constants for certain representations.*)
		Clear @ TraceNormalization;
		Clear @ Casimir2;
		
		(*del[rep, a, b] is the symbol for Kronecker delta \delta_{a,b} belong to the indices specified by the representation.*)
		Clear @ del;
		del /: del[rep_, a___, x_, b___] del[rep_, c___, x_, d___] := del[rep, c, a, b, d];
		del /: Power[del[rep_, a_, b_], 2] := Dim[rep];
		del[rep_, a_, a_] = Dim[rep];
		del[Bar[rep_], a_, b_] = del[rep, a, b];
		
		(*del for taking specific index*)
		Clear @ delIndex;
		delIndex /: del[rep_, a___, x_Symbol, b___] delIndex[rep_, c___, x_Symbol, d___] = delIndex[rep, c, a, b, d];
		delIndex /: delIndex[rep_, a___, x_Symbol, b___] delIndex[rep_, c___, x_Symbol, d___] = delIndex[rep, c, a, b, d];
		delIndex /: Power[delIndex[rep_, a_, b_], 2] = 1;
		delIndex[rep_, a_Integer, b_Integer] := KroneckerDelta[a, b];
		
		(*Default properties for anti-symmetric inavariants.*)
		Clear @ eps;
		eps /: del[rep_, a___, x_, b___] eps[rep_, c___, x_, d___] := eps[rep, c, a, b, d];
		eps /: eps[rep_ ,a_, b_] eps[rep_ ,b_, c_] := -del[rep, a, c];
		eps /: eps[rep_ ,a_, b_] eps[rep_ ,c_, b_] := del[rep, a, c];
		eps /: eps[rep_ ,b_, a_] eps[rep_ ,b_, c_] := del[rep, a, c];
		eps /: Power[eps[rep_, a_, b_], 2] := Dim[rep];
		eps[rep_, a_, a_] := 0;
		
		(*Invariant of S2 and fundamentals*)
		Clear @ delS2;
		delS2 /: del[group_[S2], a___, x_ , b___] delS2[group_, x_, i_, j_]=  delS2[group, a, b, i, j];
		delS2 /: del[group_[fund], k___, x_ , l___] delS2[group_, a_, i___ x_, j___]=  delS2[group, a, i, k, l, j];
		delS2 /: delS2[group_, a_, i_ , j_] delS2[group_, a_, k_, l_] = Sym[k, l][del[group@ fund, i, k] del[group@ fund, j, l]];
		delS2 /: delS2[group_, a_, i_ , j_] delS2[group_, b_, i_, j_] = del[group@ S2, a, b];
		delS2 /: delS2[group_, a_, i_ , j_] delS2[group_, b_, j_, i_] = del[group@ S2, a, b];
		
		(*Invariant of A2 and fundamentals*)
		Clear @ delA2;
		delA2 /: del[group_[A2], a___, x_ , b___] delA2[group_, x_, i_, j_]=  delA2[group, a, b, i, j];
		delA2 /: del[group_[fund], k___, x_ , l___] delA2[group_, a_, i___ x_, j___]=  delA2[group, a, i, k, l, j];
		delA2 /: delA2[group_, a_, i_ , j_] delA2[group_, a_, k_, l_] = AntiSym[k, l][del[group@ fund, i, k] del[group@ fund, j, l]];
		delA2 /: delA2[group_, a_, i_ , j_] delA2[group_, b_, i_, j_] = del[group@ A2, a, b];
		delA2 /: delA2[group_, a_, i_ , j_] delA2[group_, b_, j_, i_] = - del[group@ A2, a, b];
		
		(*Default group generator properties.*)
		Clear @ tGen;
		tGen /: del[rep_, a___, x_, b___] tGen[rep_, A_, c___, x_, d___] := tGen[rep, A, c, a, b, d];
		tGen /: del[group_[adj], A___, X_, B___] tGen[group_[rep_], X_, c__] := tGen[group[rep], A, B, c];
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
		fStruct /: del[group_[adj], A___, X_, B___] fStruct[group_, C___, X_, D___] := fStruct[group, C, A, B, D];
		(*fStruct /: fStruct[group_, A___] fStruct[group_, B___] /; Length @ Intersection[List @A, List @B] === 2 := 
			Casimir2 @ group[adj] CasimirSig[group, List@A, List@B];*)
		fStruct /: fStruct[group_, x:OrderlessPatternSequence[A_, B_, C_] ] fStruct[group_, y:OrderlessPatternSequence[D_, B_, C_] ] := 
			Signature[List@ x] Signature[List@ y] Signature[{A, B, C}] Signature[{D, B, C}] Casimir2[group@ adj] del[group@ adj, A, D];
		fStruct[group_, a___, x_, b___, x_, c___] := 0;
	];

(*Function applying the group *)
RefineGroupStructures[expr_] := Block[{replace},
	replace = {
			tGen[group_[adj], A_, B_, C_] :> - Module[{a, b, c}, 2 AntiSym[A, B][tGen[group@ fund, A, a, b] tGen[group@ fund, B, b, c]] tGen[group@ fund, C , c, a] 
				] / TraceNormalization @ group @ fund,
			tGen[group_[S2], A_, a_, b_] :> 2 Sym[a@1, a@2] @ Sym[b@1, b@2][tGen[group@ fund, A, a@1, b@1] del[group@ fund, a@2, b@2]],
			delS2[group_, a_, i_, j_] :> Sym[a@1, a@2][ del[group@ fund, a@1, i] del[group@ fund, a@2, j]],
			tGen[group_[A2], A_, a_, b_] :> 2 AntiSym[a@1, a@2] @ AntiSym[b@1, b@2][tGen[group@ fund, A, a@1, b@1] del[group@ fund, a@2, b@2]],
			delA2[group_, a_, i_, j_] :> AntiSym[a@1, a@2][ del[group@ fund, a@1, i] del[group@ fund, a@2, j]]  
		};
	Return[expr /. replace];	
];

(*Adds case to the system built in function Tr and Dot, to deal with substituting couplings for 0.*)
	Unprotect[Tr, Dot];
		Dot[___, 0, ___] = 0;
		Tr[0] = 0;
	Protect[Tr, Dot];

(*Complex conjugation and transposition of matrices*)
	Bar[Bar[x_]] = x;
	Bar[0] = 0;	
	Trans[Trans[x_]] := x;
	Trans[0] = 0;
	Trans[a_List] /; VectorQ[a] := a;
(*Matrix head for symbolic matrix manipulations*)
	Matrix /: del[ind_, a___, x_, b___] Matrix[m__][c___, ind_[x_], d___] := Matrix[m][c, ind[a, b], d];
	Matrix /: Matrix[m1__][a_] Matrix[m2__][a_] := Matrix[Sequence @@ Reverse[Trans /@ List@m1], m2][];
	Matrix /: Matrix[m1__][b_] Matrix[m2__][b_, c_] := Matrix[Sequence @@ Reverse[Trans /@ List@m2], m1][c];
	Matrix /: Matrix[m1__][a_, b_] Matrix[m2__][b_] := Matrix[m1, m2][a];
	Matrix /: Matrix[m1__][a_, b_] Matrix[m2__][b_, c_] := Matrix[m1, m2][a, c];
	Matrix /: Matrix[m1__][a_, b_] Matrix[m2__][c_, b_] := Matrix[m1, Sequence @@ Reverse[Trans /@ List@m2]][a, c];
	Matrix /: Matrix[m1__][b_, a_] Matrix[m2__][b_, c_] := Matrix[Sequence @@ Reverse[Trans /@ List@m1], m2][a, c];
	Matrix[a_List, b_List] = Matrix[Dot[a, b]];
	Matrix[m_][] := m;
	Matrix[m_][Null] := m;
	Matrix[___, 0, ___][___] = 0;
	(*Matrix[m__][a_, a_] := Tr @ Dot[m];*)
	Matrix[m__][a_, a_] :=
		Block[{matrices, permutations},
			matrices = List @ m;
			permutations = NestList[RotateLeft, matrices, Length@ matrices];
			permutations = Join[permutations, Map[Trans, Reverse[permutations, 2], {2}] ]; 
			matrices = permutations[[Ordering[permutations, 1][[1]] ]];
			Tr @ Dot[Sequence @@ matrices]
		];
		
(*Formating*)
	Format[Trans[x_]] := HoldForm[x^Global`T];
	Format[Bar[x_]] := OverBar @ x;
	Format[Trans[Bar[x_]] ] := x^Style[Global`\[Dagger], Bold, 12];
	Format[Matrix[x_, y__][] ] := HoldForm @Dot[x, y];
	Format[Matrix[x__][h1_[i1_] ] ] := Subscript[Dot[x], i1];
	Format[Matrix[x__][h1_[i1_], h2_[i2_]] ] := Subscript[Dot[x], i1, i2];
	

(*###########################################*)
(*----------Gauge group definitions----------*)
(*###########################################*)
(*Initialization for an SO(n) gauge group.*)
DefineSOGroup[group_Symbol, n_Integer|n_Symbol] := 
	Block[{projection},		
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
		
		(*Adjoint*)
		Dim[group[adj]] = n (n - 1) / 2;
		TraceNormalization[group[adj]] = (n - 2) / 2;
		Casimir2[group[adj]] = (n - 2) / 2;
	];

(*Initialization for an Sp(n) gauge group.*)
DefineSpGroup[group_Symbol, n_Integer|n_Symbol] := 
	Block[{projection},		
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
		
		(*Adjoint*)
		Dim[group[adj]] = n (n + 1) / 2;
		TraceNormalization[group[adj]] = (n + 2) / 2;
		Casimir2[group[adj]] = (n + 2) / 2;
	];

(*Initialization for an SU(n) Lie group.*)
DefineSUGroup[group_Symbol, n_Integer|n_Symbol] :=
	Block[{},
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
		
		(*Anti-Symmetric*)
		Dim[group[A2]] = n (n -1) /2;
		TraceNormalization[group[A2]] = (n - 2) / 2;
		Casimir2[group[A2]] =  (n - 2) (n + 1) / n;  
		
		(*Adjoint*)
		Dim[group[adj]] = n^2 - 1;
		TraceNormalization[group[adj]] = n;
		Casimir2[group[adj]] = n;
	];

(*Initialization for a U(1) gauge group.*)	
DefineU1Group[group_Symbol, power_Integer:1] := 
	Block[{},
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


End[]



	