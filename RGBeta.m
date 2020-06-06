(*
################################################################
This is RGBeta, Version 0.4.
Last modefied 23-04-2020

Released under the MIT license (see 'MIT_license.txt').

Author: Anders Eller Thomsen
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

$scalarMasses::usage =
	"$scalarMasses is an association with all the scalar masses that have been declared in the model."

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

f1::usage = f2::usage = s1::usage = s2::usage = s3::usage = s4::usage = v1::usage = v2::usage =
	"Index name used for flavor indices on kinetic-mixing, Yukawa, or quartic couplings."

del::usage =
	"del[rep, a, b] represents a Kronecker delta in \"rep space\" with indices a and b."

delA2::usage =
	"delA2[group, a, i, j] represents the singlet contraction between the A2 index a and the two fundamental indices i and j."

delIndex::usage =
	"del[rep, a, b] represents a Kronecker delta function in \"rep space\" with either a or b being an Integer fixing the value of the other summation index."

delS2::usage =
	"delS2[group, a, i, j] represents the singlet contraction between the S2 index a and the two fundamental indices i and j."

eps::usage =
	"eps[rep, a, b] represents the 2-index antisymmetric tensor in \"rep space\" with indices a and b."

fStruct::usage =
	"fStruct[group, A, B, C] represents the structure constant of the group with indices A, B, and C."

lcSymb::usage =
	"Levi-Civita symbol"

sDelF::usage =
	"sDelF[field, generalIndex, specificIndex] represents the structure delta for a fermionic field."

sDelS::usage =
	"sDelS[field, generalIndex, specificIndex] represents the structure delta for a scalar field."

sDelV::usage =
	"sDelV[field, generalIndex, specificIndex] represents the structure delta for a vector field."

tGen::usage =
	"tGen[rep, A, a, b] represents a group generator of the representation \"rep\" with adjoint index A. a and b are the two indices of rep. "

BarToConjugate::usage = Chirality::usage = Coupling::usage = CouplingBar::usage = CouplingIndices::usage = CouplingMatrix::usage = FermionMass::usage = Field::usage =
Fields::usage = FlavorIndices::usage = FourDimensions::usage = GaugeRep::usage = GroupInvariant::usage = Indices::usage = Invariant::usage = LieGroup::usage = Mass::usage=
MassIndices::usage = Parameterizations::usage = Projector::usage = Quartic::usage = RescaledCouplings::usage = ScalarMass::usage = SelfConjugate::usage = Trilinear::usage = UniqueArrangements::ussage = Yukawa::usage =
	"Function option and/or key used in global association lists of fields and/or couplings."

SO::usage = Sp::usage = SU::usage = U1::usage =
	"SO, Sp, SU, and U1 are used to specify different Lie groups."

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

AddScalarMass::ussage =
	"AddScalarMass[mass, {scal1, scal2}] defines a mass term between the two scalar fields."

AddTrilinear::usage =
	"AddTrilinear[coupling, {scal1, scal2, scal3}] defines a Yukawa interaction between a scalar and two fermion fields."

AddYukawa::usage =
	"AddYukawa[coupling, {scal, ferm1, ferm2}] defines a Yukawa interaction between a scalar and two fermion fields."

AddVector::usage =
	"AddVector[field, group] is an internal function used to define a vector field of a gauge group."

AntiSym::usage =
	"AntiSym[a, b][expr] is an internal function for antisymmetrizing expr in a and b."

AveragePermutations::usage =
	"AveragePermutations[indices, permutations][expr] gives the average of expr with indices (List) exchanged in all way perscribed by permutations."

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

CheckProjection::usage =
	"CheckProjection[coupling] returns the result of the automatic projection operator of the coupling on the corresponding generalized coupling."

CouplingPermutations::usage =
	"CouplingPermutations[fields, invariant] provides at list of all unique permutations of a coupling of the fields (List) and with the given group invariance function."

DefineLieGroup::usage =
	"DefineLieGroup[groupName, lieGroup[n]] sets groupName to be a lie group of type SU(n), SO(n), Sp(n), or U(1)^n and defines invariants for several common representations of said group."

Dim::usage =
	"Dim[rep] sets the dimension of a given representation."

FermionMassTensors::usage =
	"FermionMassTensors[loop] is a function that computes all the tensor contractions used the general fermion mass beta function at the given loop order."

Finalize::usage =
	"Finalize[expr] performs additional refinement of beta function expressions. Only works for expressions with up to 2 nontrivial flavor indices."

FGauge::usage =
	"FGauge[A, B, C] is a function that generates the general gauge structure constants."

ResetBetas::usage =
	"FlsuhBetas[] is a function used to dump all internally stored beta computations from the kernel."

G2Matrix::usage =
	"G2Matrix[a, b] is a function that generates the general gauge coupling matrix."

GaugeTensors::usage =
	"GaugeTensors[loop] is a function that computes all the tensor contractions used the general gauge beta function at the given loop order."

Lam::usage =
	"Lam[a, b, c, d] is a function that generates the general quartic coupling."

Matrix::usage =
	"Matrix[x,...][i, j] represents the matrix product of couplings x,... with open indices i and j."

OptionsCheck::usage =
	"Function internally used to check validity of options passed as arguments."

QuarticBetaFunctions::usage =
	"QuarticBetaFunctions[loop] returns all quartic beta functions to the given loop order using diagonalized projectors."

QuarticTensors::usage =
	"QuarticTensors[loop] is a function that computes all the tensor contractions used the general quartic beta function at the given loop order."

RemoveField::usage =
	"RemoveField[field, ...] is a function that removes one or more scalar and fermion fields from the model."

RemoveInteraction::usage =
	"RemoveCoupling[coupling, ...] is a function that removes one or more couplings and corresponding interactions/gauge groups from the model."

ResetModel::usage =
	"ResetModel[] resets all tensors and removes all fields and coupling definitions made in the current instance of RGBeta."

RefineGroupStructures::usage =
	"RefineGroupStructures[expr] decomposes generators of non-fundamental representations of the groups to the fundamental ones, whereby identities can be applied."

ReInitializeSymbols::usage =
	"ReInitializeSymbols[] is a function which when called flushes all previous definitions for symbol behaviour under implicit summation."

RepresentationCheck::usage =
	"RepresentationCheck[rep] returns True if rep is a valid representation given the gauge groups of the present model and False otherwise."

ScalarMassiveTensors::usage =
	"ScalarMassiveTensors[loop] is a function that computes all the tensor contractions used in the the trillinear and scalar mass beta functions at the given loop order."

SetReal::usage =
	"SetReal[x1,...] makes x1,... behave as real parameters under complex conjugation (the Bar function)."

Sym::usage =
	"Sym[a1, a2, a3, a4][expr] is an internal function for symmetrizing expr over up to four dummy indices."

Sym4::usage =
	"Sym4[a1, a2, a3, a4][expr] averages the expression over the 4 ways of switching a1 with one of the indices."

Tensor::usage =
	"Tensor[x][i,...] represents the tensor x with open indices i,..."

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
	$RGBetaVersion = "RGBeta v0.4";
	Print[$RGBetaVersion, " by Anders Eller Thomsen"];

	(*Loads package files*)
	dir = DirectoryName @ $InputFileName;
	Block[{$Path = {FileNameJoin@ {dir,"Core"} }},
		<< Betas`;
		<< FieldsAndCouplings`;
		<< GroupsAndIndices`;
		<< TensorCalculations`;
	];

	(*Initializes dummy index notation and global model information*)
	ResetModel[];

	(*Protects global symbols*)
	Protect[$a, $b, $c, $d, $i, $j, $A, $B];
	Protect[SO, Sp, SU, U1];
	Protect[f1, f2, s1, s2, s3, s4];

End[] (* End Private Context *)


(*###########################################*)
(*---------------Options check---------------*)
(*###########################################*)
General::invalidopt = "Option `1` for function `2` received invalid value `3`.";
General::optexpectsval = "Option `1` for function `2` received invalid value `3`. A `4` is expected.";
OptionMessage[opt_, func_, val_] := Message[General::invalidopt, opt, func, val];
OptionMessage[Chirality, func_, val_] := Message[General::optexpectsval, Chirality, func, val, "Left or a Right"];
OptionMessage[CouplingIndices, func_, val_] := Message[General::optexpectsval, CouplingIndices, func, val, Function];
OptionMessage[FlavorIndices, func_, val_] := Message[General::optexpectsval, FlavorIndices, func, val, List];
OptionMessage[GaugeRep, func_, val_] := Message[General::optexpectsval, GaugeRep, func, val, List];
OptionMessage[GroupInvariant, func_, val_] := Message[General::optexpectsval, GroupInvariant, func, val, Function];
OptionMessage[MassIndices, func_, val_] := Message[General::optexpectsval, MassIndices, func, val, Function];
OptionMessage[Parameterizations, func_, val_] := Message[General::optexpectsval, Parameterizations, func, val, "List of replacement rules"];

(*Tests for specific options. Form expected is OptionTest[function, options] *)
OptionTest[_, BarToConjugate] = BooleanQ;
OptionTest[_, Chirality] = MatchQ[Left|Right];
OptionTest[_, CouplingIndices] = MatchQ[_Function];
OptionTest[_, CouplingMatrix] = MatchQ[Automatic|_List];
OptionTest[_, FlavorIndices] = MatchQ[_List];
OptionTest[_, FourDimensions] = BooleanQ;
OptionTest[_, GaugeRep] = MatchQ[_List];
OptionTest[_, GroupInvariant] = MatchQ[_Function];
OptionTest[_, MassIndices] = MatchQ[_Function];
OptionTest[_, Parameterizations] = MatchQ[_List];
OptionTest[_, RescaledCouplings] = BooleanQ;
OptionTest[_, SelfConjugate] = BooleanQ;
OptionTest[_, CheckInvariance] = BooleanQ;

Attributes @ OptionsCheck = {HoldFirst};
OptionsCheck @ func_[___, opts : OptionsPattern[]] :=
	And @@ (OptionTest[func, #1][#2] || OptionMessage[#1, func, #2] &) @@@ FilterRules[List[opts], Options @ func];




EndPackage[]
