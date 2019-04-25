(*
################################################################
This is RGBeta, Version 0.2. 
Last modefied xx-03-2019

Copyright...

Anders Eller Thomsen 
################################################################
*)

BeginPackage["RGBeta`"]

(*################################################################*)
(*----------Defines symbols and functions used by RGBeta----------*)
(*################################################################*)  
$a::usage = $b::usage = $c::usage = $d::usage = $i::usage = $j::usage = $A::usage = $B::usage = $vev::usage =  $vevSelect::usage =  
	"Used as a global dummy index."

$couplings::usage = 
	"$coupling is an association containing the types of all the couplings that have been defined in the model."

$fermionMasses::usage =
	"$fermionMasses is an association with all the fermion masses that have been declared in the model."

$fermions::usage =
	"$fermions is an association containing all the internal information on the fermions that have been declared in the model."

$gaugeCoefficients::usage = 
	"$gaugeCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the quartic beta function."

$gaugeGroups::usage =
	"$gaugeGroups is an association containing all the internal information on the gauge groups that have been declared in the model."

$quarticCoefficients::usage = 
	"$quarticCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the quartic beta function."

$quartics::usage = 
	"$quartics is an association containing all internal information on the quartic couplings."

$scalars::usage =
	"$scalars is an association containing all the internal information on the scalars that have been declared in the model."

$trilinears::usage =
	"$trilinears is an association contraining all the internal informaiton on the trilinear scalar couplings." 

$yukawaCoefficients::usage = 
	"$yukawaCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the yukawa beta function."

$yukawas::usage = 
	"$yukawas is an association containing all internal information on the yukawa couplings."

$RGBetaVersion::usage =
	"CurrentVersion of the RGBeta package."

adj::usage = fund::usage = S2::usage = A2::usage =
	"Representation name defined defined by the Define(SO/Sp/SU)Group functions."

f1::usage = f2::usage =  s1::usage =  s2::usage =  s3::usage =  s4::usage =
	"Index name used for flavor indices on yukawa or quartic couplings."  

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

Chirality::usage = Coupling::usage = CouplingBar::usage = FermionMass::usage = Field::usage = Fields::usage = FlavorIndices::usage = 
	GaugeRep::usage = Indices::usage = Invariant::usage = LieGroup::usage = Mass::usage= Projector::usage = Quartic::usage =
	SelfConjugate::usage = Trilinear::usage = Yukawa::usage = 
	"Key used in global association lists of fields and/or couplings."

SO::usage = Sp::usage = SU::usage = U::usage =
	"SO, Sp, SU, and U are used to specify different Lie groups."

AddFermion::usage =
	"AddFermion[field] is a function used to define a fermion field in the model."

AddFermionMass::ussage =
	"AddFermionMass[mass, {ferm1, ferm2}] defines a mass term between the two fermion fields."

AddGaugeGroup::usage =
	"AddGaugeGroup[coupling, groupName, lieGroup[n]] is a function used to define a gauge group in the model."

AddQuartic::usage = 
	"AddQuartic[coupling, {scal1, scal2, scal3, scal4}] defines a quartic interaction between 4 scalar fields."
	
AddScalar::usage =
	"AddScalar[field] is a function used to define a scalar field in the model."

AddTrilinear::usage = 
	"AddTrilinear[coupling, {scal1, scal2, scal3}] defines a Yukawa interaction between a scalar and two fermion fields."

AddYukawa::usage = 
	"AddYukawa[coupling, {scal, ferm1, ferm2}] defines a Yukawa interaction between a scalar and two fermion fields."

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

Bcoef::usage = 
	"Bcoef[i, j, k] is used internally to denote the coefficients of the differnet tensor contraction appering in the beta functions."

Casimir2::usage = 
	"Casimir2[rep] sets the quadratic casimir of a given representation."

CasimirSig::usage =
	"CasimirSig is internal function used to determined the relative sign of the contraction of two structure constants."

CheckProjection::usage =
	"CheckProjection[coupling] returns the result of the automatic projection operator of the coupling on the corresponding generalized coupling."

DefineU1Group::usage =
	"DefineU1Group[group, n] defines group to be a U(1) group and defines representations accordingly."

DefineSpGroup::usage =
	"DefineSpGroup[group, n] defines group to be an Sp(n) group and sets defines invariants for several common representations accordingly."

DefineSOGroup::usage =
	"DefineSOGroup[group, n] defines group to be an SO(n) group and sets defines invariants for several common representations accordingly."

DefineSUGroup::usage =
	"DefineSUGroup[group, n] defines group to be an SU(n) group and sets defines invariants for several common representations accordingly."

Dim::usage =
	"Dim[rep] sets the dimension of a given representation."

Finalize::usage = 
	"Finalize[expr] performs additional refinement of beta function expressions."

FGauge::usage = 
	"FGauge[A, B, C] is a function that generates the general gauge structure constants."

G2Matrix::usage = 
	"G2Matrix[a, b] is a function that generates the general gauge coupling matrix."

GaugeTensors::usage = 
	"GaugeTensors[loop] is a function that computes all the tensor contractions used the general gauge beta function at the given loop order."

Lam::usage = 
	"Lam[a, b, c, d] is a function that generates the general quartic coupling."

Matrix::usage =
	"Matrix[x,...][i, j] represents the matrix product of couplings x,... with open indices i and j."

QuarticTensors::usage = 
	"QuarticTensors[loop] is a function that computes all the tensor contractions used the general quartic beta function at the given loop order."

RemoveField::usage =
	"RemoveField[field, ...] is a function that removes one or more scalar and fermion fields from the model."

RemoveInteraction::usage =
	"RemoveCoupling[coupling, ...] is a function that removes one or more couplings and corresponding interactions/gauge groups from the model."

ResetModel::usage = 
	"ResetModel[] resets all tensors and removes all fields and coupling definitions made in the current instance of RGBeta."

ReInitializeSymbols::usage = 
	"ReInitializeSymbols[] is a function which when called flushes all previous definitions for symbol behaviour under implicit summation."

Sym::usage =
	"Sym[a1, a2, a3, a4][expr] is an internal function for symmetrizing expr over up to four dummy indices."

Sym4::usage =
	"Sym4[a1, a2, a3, a4][expr] averages the expression over the 4 ways of switching a1 with one of the indices."

Tferm::usage = 
	"Tferm[A, i, j] is a function that generates the general gauge generator matrix for the fermions."

TfermTil::usage = 
	"TfermTil[A, i, j] is a function that generates the general gauge generator matrix for the fermions with an applied tilde."

TraceNormalization::usage = 
	"TraceNormalization[rep] sets the trace normalization of a given representation."

Trans::usage = 
	"Trans[coupling] represents the transposed quantity of a coupling with two indices."

Tscal::usage = 
	"Tscal[A, i, j] is a function that generates the general gauge generator matrix for the scalars."

Tdot::usage = 
	"Tdot[a, b,...] is an internal function the sequantially expands the matrix product of all the arguments."

Ttimes::usage = 
	"Ttimes[a, b,...] is an internal function the sequantially expands the product of all the arguments."

Yuk::usage = 
	"Yuk[a, i, j] is a function that generates the general yukawa coupling."

YukawaTensors::usage = 
	"YukawaTensors[loop] is a function that computes all the tensor contractions used the general yukawa beta function at the given loop order."

YukTil::usage = 
	"YukTil[a, i, j] is a function that generates the general yukawa coupling with an applied tilde."

(*##############################################*)
(*---------------Loads components---------------*)
(*##############################################*) 
Begin["`Private`"] (* Begin Private Context *) 
	$RGBetaVersion = "RGBeta 0.2";
	Print[$RGBetaVersion, " by Anders Eller Thomsen"];
	
	(*Loads package files*)
	dir = DirectoryName @ $InputFileName;
	Block[{$Path = {dir}},
		<< Betas`;
		<< FieldsAndCouplings`;
		<< GroupsAndIndices`;
		<< Tensors`;
	];
	
	(*Initializes dummy index notation and global model information*)
	ResetModel[];

	(*Protects global symbols*)	
	Protect[$a, $b, $c, $d, $i, $j, $A, $B];
	Protect[SO, Sp, SU, U];
	Protect[f1, f2, s1, s2, s3, s4];
	
End[] (* End Private Context *)

EndPackage[]