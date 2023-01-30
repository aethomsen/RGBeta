(* Generate warning for Mathematica versions earlier than 11.0.0 *)
If[$VersionNumber < 11.0,
	CellPrint[{
		TextCell["RGBeta` was developed for Mathematica 11.0.0 and later.
		Your current Mathematica version [" <> ToString@$Version <> "] might not be compatible.
		In case you are experiencing problems with RGBeta`, please update your Mathematica version.",
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

	(*Initializes dummy index notation and global model information*)
	RGBeta`ResetModel[];

	(*Print the packagename*)
	CellPrint@ Cell[StandardForm@ RowBox@ {
		AdjustmentBox[
			StyleBox["\t    R", FontColor-> RGBColor[.6, .0706, 0.1373]],
		BoxMargins-> {{0, -.2}, {0, 0}}],
		AdjustmentBox[
			StyleBox["G", FontColor-> RGBColor[0.3, 0.55, 0.2]],
		BoxMargins -> {{0, -.2}, {0, 0}}],
		StyleBox["Beta", FontColor-> RGBColor[0, 0.396, 0.741]]
		}, "text", FontSize-> 20, FontWeight-> Bold,
	Background-> White];

	Print[
	"by Anders Eller Thomsen \n",
	"Reference: ",Hyperlink["arXiv:2101.08265","https://arxiv.org/abs/2101.08265"],"\n",
	"Website: ",Hyperlink["https://github.com/aethomsen/RGBeta","https://github.com/aethomsen/RGBeta"]
	];
];
