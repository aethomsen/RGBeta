(* ::Package:: *)

(* ::Title:: *)
(*Initialization [RGBeta`]*)


(* Generate warning for Mathematica versions earlier than 11.0.0 *)
If[$VersionNumber < 11.0, 
CellPrint[{
TextCell["RGBeta` was developed for Mathematica 11.0.0 and later. 
Your current Mathematica version [" <> ToString@$Version <> "] might not be compatible. 
In case you are experiencing problems with SuperTracer`, please update your Mathematica version.",
"Text",Background->LightRed]
}]
]


(* Loading the package *)
If[MemberQ[$Packages,"RGBeta`"],
	(* Avoid double loading the package *)
	Print[Style["The package RGBeta` is already loaded. Please restart the kernel to reload it.",RGBColor[.8,.4706,0.2573]]],
	
	(* Set directory of AutoEFT package *)
	$DirectoryRGBeta= DirectoryName[$InputFileName, 2];

	(* Loading *)
	Check[
		Get[FileNameJoin[{$DirectoryRGBeta, "RGBeta.m"}]],
		Print[Style["Loading failed!",RGBColor[.6,.0706,0.1373]]];
		Abort[]
	];
	
	(*Print the packagename*)
	CellPrint[TextCell["RGBeta","Text",Blue,Background->LightGreen,FontSize->14]];

	Print[
	"by Anders Eller Thomsen \n",
	"Reference: ",Hyperlink["arXiv:2021.xxxx","https://arxiv.org/abs/2021.xxxxx"],"\n",
	"Website: ",Hyperlink["https://github.com/aethomsen/RGBeta","https://github.com/aethomsen/RGBeta"]
	];
];
