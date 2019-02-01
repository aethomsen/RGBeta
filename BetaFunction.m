<< GroupsAndIndices`
<< FieldsAndCouplings`

(*Symmetrizing in indices*)
Sym[i1_, i2_][expr_] := expr /2 + ReplaceAll[expr, {i -> j, j -> i}] /2 ;
AntiSym[i1_, i2_][expr_] := expr /2 + ReplaceAll[expr, {i -> j, j -> i}] /2 ;

(*Functions that speed up evaluation by applying expand succesively to each couple of terms in the evaluation.*)
Ttimes[a_, b_, c___] := Ttimes[Expand[a b], c];
Ttimes[a_] = a;
Tdot[a_, b_, c___] := Tdot[Expand[a.b], c];
Tdot[a_] = a;

