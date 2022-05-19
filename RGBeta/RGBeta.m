(*
################################################################
This is RGBeta, Version 1.1.2.
Last modefied 19-05-2022

Released under the MIT license (see 'LICENSE').

Author: Anders Eller Thomsen
################################################################
*)

Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageExport["$a"]
PackageExport["$b"]
PackageExport["$c"]
PackageExport["$d"]
PackageExport["$i"]
PackageExport["$j"]
PackageExport["$A"]
PackageExport["$B"]
PackageExport["$vev"]
PackageExport["$vevSelect"]

PackageExport["f1"]
PackageExport["f2"]
PackageExport["s1"]
PackageExport["s2"]
PackageExport["s3"]
PackageExport["s4"]
PackageExport["v1"]
PackageExport["v2"]

PackageExport["$RGBetaVersion"]

PackageExport["AntisymmetricIndices"]
PackageExport["BarToConjugate"]
PackageExport["Chirality"]
PackageExport["CouplingIndices"]
PackageExport["CouplingMatrix"]
PackageExport["FlavorIndices"]
PackageExport["FlavorImproved"]
PackageExport["FourDimensions"]
PackageExport["GaugeRep"]
PackageExport["GroupInvariant"]
PackageExport["MassIndices"]
PackageExport["Parameterizations"]
PackageExport["RescaledCouplings"]
PackageExport["SelfConjugate"]
PackageExport["SymmetricIndices"]

PackageScope["$fermion"]
PackageScope["$gauge"]
PackageScope["$scalar"]
PackageScope["$scalarContraction"]

PackageScope["OptionsCheck"]

PackageScope["Coupling"]
PackageScope["CouplingBar"]
PackageScope["FermionMass"]
PackageScope["Field"]
PackageScope["Fields"]
PackageScope["Indices"]
PackageScope["Invariant"]
PackageScope["LieGroup"]
PackageScope["Mass"]
PackageScope["Projector"]
PackageScope["Quartic"]
PackageScope["ScalarMass"]
PackageScope["Trilinear"]
PackageScope["UniqueArrangements"]
PackageScope["Yukawa"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

$a::usage = $b::usage = $c::usage = $d::usage = $i::usage = $j::usage = $A::usage = $B::usage = $vev::usage =  $vevSelect::usage =
	"Used as a global dummy index."

$RGBetaVersion::usage =
	"CurrentVersion of the RGBeta package."

OptionsCheck::usage =
	"Function internally used to check validity of options passed as arguments."

(*################################################################*)
(*----------Defines symbols and functions used by RGBeta----------*)
(*################################################################*)

$fermion::usage = $gauge::usage = $scalar::usage =
	"Denotes index types for TStructure."

$scalarContraction::usage =
	"Contraction array for the scalar fields."

f1::usage = f2::usage = s1::usage = s2::usage = s3::usage = s4::usage = v1::usage = v2::usage =
	"Index name used for flavor indices on kinetic-mixing, Yukawa, or quartic couplings."

(* sDelF::usage =
	"sDelF[field, generalIndex, specificIndex] represents the structure delta for a fermionic field."

sDelS::usage =
	"sDelS[field, generalIndex, specificIndex] represents the structure delta for a scalar field."

sDelV::usage =
	"sDelV[field, generalIndex, specificIndex] represents the structure delta for a vector field." *)

BarToConjugate::usage = Chirality::usage = CouplingIndices::usage = CouplingMatrix::usage =  FlavorIndices::usage = FlavorImproved::usage = FourDimensions::usage =
GaugeRep::usage = GroupInvariant::usage = MassIndices::usage = Parameterizations::usage = RescaledCouplings::usage = SelfConjugate::usage =
	"Function option and/or key used in global association lists of fields and/or couplings."

Coupling::usage= CouplingBar::usage = FermionMass::usage = Field::usage = Fields::usage = Indices::usage = Invariant::usage =
LieGroup::usage = Mass::usage = Projector::usage = Quartic::usage = ScalarMass::usage = Trilinear::usage = UniqueArrangements::ussage =
Yukawa::usage =
	"Function option and/or key used in global association lists of fields and/or couplings."

AntisymmetricIndices::usage= SymmetricIndices::usage= "Coupling option specifying sets of completely (anti)symmetric indices. Option is given on the {{2, 3}, {1, ...}, ...}, where each set is a set of (anti)symmetrized couplings."

(*Protects global symbols*)
Protect[$a, $b, $c, $d, $i, $j, $A, $B];
Protect[SO, Sp, SU, U1];
Protect[f1, f2, s1, s2, s3, s4];

$RGBetaVersion= "v1.1.2"


(*###########################################*)
(*---------------Options check---------------*)
(*###########################################*)
General::invalidopt = "Option `1` for function `2` received invalid value `3`.";
General::optexpectsval = "Option `1` for function `2` received invalid value `3`. A `4` is expected.";
OptionMessage[opt_, func_, val_] := Message[General::invalidopt, opt, func, val];
OptionMessage[AntisymmetricIndices, func_, val_] := Message[General::optexpectsval, SymmetricIndices, func, val, "list of integers or a list of lists of integers"];
OptionMessage[Chirality, func_, val_] := Message[General::optexpectsval, Chirality, func, val, "Left or a Right"];
OptionMessage[CouplingIndices, func_, val_] := Message[General::optexpectsval, CouplingIndices, func, val, Function];
OptionMessage[FlavorIndices, func_, val_] := Message[General::optexpectsval, FlavorIndices, func, val, List];
OptionMessage[GaugeRep, func_, val_] := Message[General::optexpectsval, GaugeRep, func, val, List];
OptionMessage[GroupInvariant, func_, val_] := Message[General::optexpectsval, GroupInvariant, func, val, Function];
OptionMessage[MassIndices, func_, val_] := Message[General::optexpectsval, MassIndices, func, val, Function];
OptionMessage[Parameterizations, func_, val_] := Message[General::optexpectsval, Parameterizations, func, val, "list of replacement rules"];
OptionMessage[SymmetricIndices, func_, val_] := Message[General::optexpectsval, SymmetricIndices, func, val, "list of integers or a list of lists of integers"];

(*Tests for specific options. Form expected is OptionTest[function, options] *)
OptionTest[_, AntisymmetricIndices] = MatchQ[{} | {_Integer ..} | {{_Integer ..} ..}];
OptionTest[_, BarToConjugate] = BooleanQ;
OptionTest[_, Chirality] = MatchQ[Left|Right];
OptionTest[_, CouplingIndices] = MatchQ[_Function];
OptionTest[_, CouplingMatrix] = MatchQ[Automatic|_List];
OptionTest[_, FlavorIndices] = MatchQ[_List];
OptionTest[_, FlavorImproved] = BooleanQ;
OptionTest[_, FourDimensions] = BooleanQ;
OptionTest[_, GaugeRep] = MatchQ[_List];
OptionTest[_, GroupInvariant] = MatchQ[_Function];
OptionTest[_, MassIndices] = MatchQ[_Function];
OptionTest[_, Parameterizations] = MatchQ[_List];
OptionTest[_, RescaledCouplings] = BooleanQ;
OptionTest[_, SelfConjugate] = BooleanQ;
OptionTest[_, SymmetricIndices] = MatchQ[{} | {_Integer ..} | {{_Integer ..} ..}];
OptionTest[_, CheckInvariance] = BooleanQ;

Attributes @ OptionsCheck = {HoldFirst};
OptionsCheck @ func_[___, opts : OptionsPattern[]] :=
	And @@ (OptionTest[func, #1][#2] || OptionMessage[#1, func, #2] &) @@@ FilterRules[List[opts], Options @ func];
