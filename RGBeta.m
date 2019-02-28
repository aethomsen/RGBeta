(*################################################################*)
(*	
This is RGBeta, Version 0.1. 
Last modefied xx-xx-2019

Copyright...

Anders Eller Thomsen 
*)
(*################################################################*)


BeginPackage["RGBeta`"]

(*################################################################*)
(*----------Defines symbols and functions used by RGBeta----------*)
(*################################################################*)  
$couplings::usage = 
	"$coupling is an association containing the types of all the couplings that have been defined in the model."

$fermions::usage =
	"$fermions is an association containing all the internal information on the fermions that have been declared in the model."

$gaugeGroups::usage =
	"$gaugeGroups is an association containing all the internal information on the gauge groups that have been declared in the model."

$scalars::usage =
	"$scalars is an association containing all the internal information on the scalars that have been declared in the model."

$RGBetaVersion::usage =
	"CurrentVersion of the RGBeta package."

del::usage = 
	"del[rep, a, b] represents a Kronecker delta in \"rep space\" with indices a and b."

delIndex::usage = 
	"del[rep, a, b] represents a Kronecker delta function in \"rep space\" with either a or b being an Integer fixing the value of the other summation index."

eps::usage = 
	"eps[rep, a, b] represents the 2-index antisymmetric tensor in \"rep space\" with indices a and b."

fStruct::usage = 
	"fStruct[group, A, B, C] represents the structure constant of the group with indices A, B, and C."

tGen::usage =
	"tGen[rep, A, a, b] represents a group generator of the representation \"rep\" with adjoint index A. a and b are the two indices of rep. "

Bar::usage =
	"Bar[x] rerpresents the conjugate of x. Used both for fields, couplings, and representations."

Casimir2::usage = 
	"Casimir2[rep] sets the quadratic casimir of a given representation."

CasimirSig::usage =
	"CasimirSig is internal function used to determined the relative sign of the contraction of two structure constants."

DefineSOGroup::usage =
	"DefineSOGroup[group, n] defines group to be an SO(n) group and sets defines invariants for several common representations accordingly."

DefineSUGroup::usage =
	"DefineSUGroup[group, n] defines group to be an SU(n) group and sets defines invariants for several common representations accordingly."

DefineU1Group::usage =
	"DefineU1Group[group, n] defines group to be a U(1) group and defines representations accordingly."

Dim::usage =
	"Dim[rep] sets the dimension of a given representation."

Matrix::usage =
	"Matrix[x,...][i, j] represents the matrix product of couplings x,... with open indices i and j."

TraceNormalization::usage = 
	"TraceNormalization[rep] sets the trace normalization of a given representation."

Trans::usage = 
	"Trans[coupling] represents the transposed quantity of a coupling with two indices."


Begin["`Private`"] (* Begin Private Context *) 
	$RGBetaVersion = "RGBeta 0.1";
	Print @ $RGBetaVersion;
	Print @ "By Anders Eller Thomsen"
	
	dir = DirectoryName[$InputFileName];
	Block[{$Path = {dir}},
		<< GroupsAndIndices`;
	];
	
End[] (* End Private Context *)

EndPackage[]