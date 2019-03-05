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
$a::usage = $b::usage = $c::usage = $d::usage = $i::usage = $j::usage = $A::usage = $B::usage = 
	"Used as a global dummy index."

$couplings::usage = 
	"$coupling is an association containing the types of all the couplings that have been defined in the model."

$fermions::usage =
	"$fermions is an association containing all the internal information on the fermions that have been declared in the model."

$gaugeCoefficients::usage = 
	"$gaugeCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the quartic beta function."

$gaugeGroups::usage =
	"$gaugeGroups is an association containing all the internal information on the gauge groups that have been declared in the model."

$quarticCoefficients::usage = 
	"$quarticCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the quartic beta function."

$scalars::usage =
	"$scalars is an association containing all the internal information on the scalars that have been declared in the model."

$yukawaCoefficients::usage = 
	"$yukawaCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the yukawa beta function."

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

sDelF::usage = 
	"sDelF[field, generalIndex, specificIndex] represents the structure delta for a fermionic field."

sDelS::usage = 
	"sDelS[field, generalIndex, specificIndex] represents the structure delta for a scalar field."

sDelV::usage = 
	"sDelV[field, generalIndex, specificIndex] represents the structure delta for a vector field."

tGen::usage =
	"tGen[rep, A, a, b] represents a group generator of the representation \"rep\" with adjoint index A. a and b are the two indices of rep. "

Coupling::usage = CouplingBar::usage = Field::usage = Fields::usage = FlavorIndices::usage = 
	GaugeRep::usage = Indices::usage = Invariant::usage = Projector::usage = SelfConjugate::usage =
	"Key used in global association lists of fields and/or couplings."

AddFermion::usage =
	"AddFermion[field] is a function used to define a fermion field in the model."

AddGaugeGroup::usage =
	"AddGaugeGroup[coupling, groupName, lieGroup[n]] is a function used to define a gauge group in the model."
	
AddScalar::usage =
	"AddScalar[field] is a function used to define a scalar field in the model."

AddVector::usage =
	"AddVector[field, group] is an internal function used to define a vector field of a gauge group."

AntiSym::usage = 
	"AntiSym[a, b][expr] is an internal function for antisymmetrizing expr in a and b."

Bar::usage =
	"Bar[x] rerpresents the conjugate of x. Used both for fields, couplings, and representations."

BetaTerm::usage =
	"BetaTerm[coupling, loop] computes the contribution to the beta function of the coupling at the given loop order."

BetaFunction::usage =
	"BetaFunction[coupling, loop] computes the entire beta function of the coupling to the given loop order."

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

Finalize::usage = 
	"Finalize[expr] performs additional refinement of beta function expressions."

Matrix::usage =
	"Matrix[x,...][i, j] represents the matrix product of couplings x,... with open indices i and j."

ProjectionCheck::usage =
	"ProjectionCheck[coupling] returns the result of the automatic projection operator of the coupling on the corresponding generalized coupling."

Sym::usage =
	"Sym[a1, a2, a3, a4][expr] is an internal function for symmetrizing expr over up to four dummy indices."

TraceNormalization::usage = 
	"TraceNormalization[rep] sets the trace normalization of a given representation."

Trans::usage = 
	"Trans[coupling] represents the transposed quantity of a coupling with two indices."

Ttimes::usage = 
	"Tdot[a, b,...] is an internal function the sequantially expands the matrix product of all the arguments."

Ttimes::usage = 
	"Ttimes[a, b,...] is an internal function the sequantially expands the product of all the arguments."

(*##############################################*)
(*---------------Loads components---------------*)
(*##############################################*) 
Begin["`Private`"] (* Begin Private Context *) 
	$RGBetaVersion = "RGBeta 0.1";
	Print[$RGBetaVersion, " by Anders Eller Thomsen"];
	
	dir = DirectoryName @ $InputFileName;
	Block[{$Path = {dir}},
		<< GroupsAndIndices`;
		<< FieldsAndCouplings`;
	];
	
End[] (* End Private Context *)

EndPackage[]