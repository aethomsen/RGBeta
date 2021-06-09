(*
	Author: Anders Eller Thomsen
	Released under the MIT license (see 'LICENSE').
*)
Package["RGBeta`"]

(*##################################*)
(*----------Package Export----------*)
(*##################################*)

PackageScope["$fermionAnomalousCoefficients"]
PackageScope["$gaugeCoefficients"]
PackageScope["$quarticCoefficients"]
PackageScope["$scalarAnomalousCoefficients"]
PackageScope["$yukawaCoefficients"]

(*#####################################*)
(*----------Usage Definitions----------*)
(*#####################################*)

$fermionAnomalousCoefficients::usage =
	"$fermionAnomalousCoefficients is an internal replacement lsit containing the coefficients of all the tensors structures in the fermion anomalous dimension."
$gaugeCoefficients::usage =
	"$gaugeCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the quartic beta function."
$quarticCoefficients::usage =
	"$quarticCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the quartic beta function."
$scalarAnomalousCoefficients::usage =
	"$scalarAnomalousCoefficients is an internal replacement lsit containing the coefficients of all the tensors structures in the fermion anomalous dimension."
$yukawaCoefficients::usage =
	"$yukawaCoefficients is an internal replacement list containing the coefficient of all tensor constractions used in the yukawa beta function."

(*#####################################*)
(*----------Beta Coefficients----------*)
(*#####################################*)
$quarticCoefficients = {
	Bcoef[3, 0, 1] -> 1,
	(* 1-loop *)
	Bcoef[3, 1, 1] -> 36,
	Bcoef[3, 1, 2] -> -12,
	Bcoef[3, 1, 3] -> 3,
	Bcoef[3, 1, 4] -> 2,
	Bcoef[3, 1, 5] -> -12,
	(*2-loop*)
	Bcoef[3, 2, 1] -> 324,
	Bcoef[3, 2, 2] -> -684,
	Bcoef[3, 2, 3] -> 646,
	Bcoef[3, 2, 4] -> -28,
	Bcoef[3, 2, 5] -> -32,
	Bcoef[3, 2, 6] -> 12,
	Bcoef[3, 2, 7] -> 60,
	Bcoef[3, 2, 8] -> 0,
	Bcoef[3, 2, 9] -> 6,
	Bcoef[3, 2, 10] -> -143/3,
	Bcoef[3, 2, 11] -> 11/3,
	Bcoef[3, 2, 12] -> 10/3,
	Bcoef[3, 2, 13] -> -18,
	Bcoef[3, 2, 14] -> 24,
	Bcoef[3, 2, 15] -> -18,
	Bcoef[3, 2, 16] -> 1/3,
	Bcoef[3, 2, 17] -> -6,
	Bcoef[3, 2, 18] -> 0,
	Bcoef[3, 2, 19] -> -144,
	Bcoef[3, 2, 20] -> 60,
	Bcoef[3, 2, 21] -> 10,
	Bcoef[3, 2, 22] -> 0,
	Bcoef[3, 2, 23] -> -3,
	Bcoef[3, 2, 24] -> 0,
	Bcoef[3, 2, 25] -> 24,
	Bcoef[3, 2, 26] -> -48,
	Bcoef[3, 2, 27] -> 12,
	Bcoef[3, 2, 28] -> 0,
	Bcoef[3, 2, 29] -> -2,
	Bcoef[3, 2, 30] -> -3,
	Bcoef[3, 2, 31] -> 48,
	Bcoef[3, 2, 32] -> 24,
	Bcoef[3, 2, 33] -> 24
};

$yukawaCoefficients = {
	Bcoef[2, 0, 1] -> 1,
	(* 1-loop *)
	Bcoef[2, 1, 1] -> 0,
	Bcoef[2, 1, 2] -> -6,
	Bcoef[2, 1, 3] -> 2,
	Bcoef[2, 1, 4] -> 1,
	Bcoef[2, 1, 5] -> 1/2,
	(* 2-loop *)
	Bcoef[2, 2, 1] -> -21/2,
	Bcoef[2, 2, 2] -> 12,
	Bcoef[2, 2, 3] -> 0,
	Bcoef[2, 2, 4] -> -3,
	Bcoef[2, 2, 5] -> 49/4,
	Bcoef[2, 2, 6] -> -1/4,
	Bcoef[2, 2, 7] -> -1/2,
	Bcoef[2, 2, 8] -> -97/3,
	Bcoef[2, 2, 9] -> 11/6,
	Bcoef[2, 2, 10] -> 5/3,
	Bcoef[2, 2, 11] -> 1/12,
	Bcoef[2, 2, 12] -> 12,
	Bcoef[2, 2, 13] -> 0,
	Bcoef[2, 2, 14] -> 6,
	Bcoef[2, 2, 15] -> -12,
	Bcoef[2, 2, 16] -> 10,
	Bcoef[2, 2, 17] -> 6,
	Bcoef[2, 2, 18] -> 5/2,
	Bcoef[2, 2, 19] -> 9,
	Bcoef[2, 2, 20] -> -1/2,
	Bcoef[2, 2, 21] -> -7/2,
	Bcoef[2, 2, 22] -> 0,
	Bcoef[2, 2, 23] -> -2,
	Bcoef[2, 2, 24] -> 2,
	Bcoef[2, 2, 25] -> 0,
	Bcoef[2, 2, 26] -> -2,
	Bcoef[2, 2, 27] -> 0,
	Bcoef[2, 2, 28] -> -1/2,
	Bcoef[2, 2, 29] -> -2,
	Bcoef[2, 2, 30] -> -1/4,
	Bcoef[2, 2, 31] -> -3/4,
	Bcoef[2, 2, 32] -> -1,
	Bcoef[2, 2, 33] -> -3/4
};

$gaugeCoefficients = {
	Bcoef[1, 0, 1] -> 1,
	(* 1-loop *)
	Bcoef[1, 1, 1] -> -22/3,
	Bcoef[1, 1, 2] -> 2/3,
	Bcoef[1, 1, 3] -> 1/3,
	(* 2-loop *)
	Bcoef[1, 2, 1] -> 2,
	Bcoef[1, 2, 2] -> 4,
	Bcoef[1, 2, 3] -> -68/3,
	Bcoef[1, 2, 4] -> 10/3,
	Bcoef[1, 2, 5] -> 2/3,
	Bcoef[1, 2, 6] -> -1,
	Bcoef[1, 2, 7] -> 0,
	(* 3-loop *)
	Bcoef[1, 3, 1] -> -1,
	Bcoef[1, 3, 2] -> 29/2,
	Bcoef[1, 3, 3] -> 133/18,
	Bcoef[1, 3, 4] -> 679/36,
	Bcoef[1, 3, 5] -> -11/18,
	Bcoef[1, 3, 6] -> -25/18,
	Bcoef[1, 3, 7] -> -23/36,
	Bcoef[1, 3, 8] -> -49/36,
	Bcoef[1, 3, 9] -> 4,
	Bcoef[1, 3, 10] -> 25/2,
	Bcoef[1, 3, 11] -> -2857/27,
	Bcoef[1, 3, 12] -> -79/108,
	Bcoef[1, 3, 13] -> 1/54,
	Bcoef[1, 3, 14] -> 1415/54,
	Bcoef[1, 3, 15] -> 545/108,
	Bcoef[1, 3, 16] -> -29/54,
	Bcoef[1, 3, 17] -> 1,
	Bcoef[1, 3, 18] -> -1/12,
	Bcoef[1, 3, 19] -> -5/4,
	Bcoef[1, 3, 20] -> -1/4,
	Bcoef[1, 3, 21] -> -1,
	Bcoef[1, 3, 22] -> -7,
	Bcoef[1, 3, 23] -> -7/2,
	Bcoef[1, 3, 24] -> -6,
	Bcoef[1, 3, 25] -> 9/4,
	Bcoef[1, 3, 26] -> 1,
	Bcoef[1, 3, 27] -> -1,
	Bcoef[1, 3, 28] -> 3/2,
	Bcoef[1, 3, 29] -> 7/8,
	Bcoef[1, 3, 30] -> 1/2,
	Bcoef[1, 3, 31] -> 1/8,
	Bcoef[1, 3, 32] -> 3/8,
	Bcoef[1, 3, 33] -> -1/8
};

$fermionAnomalousCoefficients = {
	(* 1-loop *)
	Acoef[1, 1, 1] -> Global`\[Xi],
	Acoef[1, 1, 2] -> 1/2,
	(* 2-loop *)
	Acoef[1, 2, 1] -> -3/2,
	Acoef[1, 2, 2] -> 25/4 + 2 Global`\[Xi] + Global`\[Xi]^2/4,
	Acoef[1, 2, 3] -> -1/4,
	Acoef[1, 2, 4] -> -1/2,
	Acoef[1, 2, 5] -> 9/2,
	Acoef[1, 2, 6] -> -1/4,
	Acoef[1, 2, 7] -> -7/4,
	Acoef[1, 2, 8] -> -1/8,
	Acoef[1, 2, 9] -> -3/8
};

$scalarAnomalousCoefficients = {
	(* 1-loop *)
	Acoef[2, 1, 1] -> Global`\[Xi] - 3,
	Acoef[2, 1, 2] -> 1/2,
	(* 2-loop *)
	Acoef[2, 2, 1] -> 3/2,
	Acoef[2, 2, 2] -> Global`\[Xi]^2/4 + 2 Global`\[Xi] - 35/3,
	Acoef[2, 2, 3] -> 11/12,
	Acoef[2, 2, 4] -> 5/6,
	Acoef[2, 2, 5] -> 1/12,
	Acoef[2, 2, 6] -> 5/2,
	Acoef[2, 2, 7] -> -1/2,
	Acoef[2, 2, 8] -> -3/4
};

Protect[$gaugeCoefficients, $quarticCoefficients, $yukawaCoefficients, $fermionAnomalousCoefficients, $scalarAnomalousCoefficients];
