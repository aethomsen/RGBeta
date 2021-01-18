InstallRGBeta::version="Warning: The package structure of `1` is only supported by Mathematica versions \[GreaterEqual] `2`. You are using Mathematica `3`.";


InstallRGBeta[]:= Module[{
        pkgDir= FileNameJoin[{$UserBaseDirectory, "Applications", "RGBeta"}],
        pkgLink= "https://github.com/aethomsen/RGBeta/archive/master.zip",
        pkgName= "RGBeta",
        minVersion= 11.0,
        questionOverwrite, tmpFile, unzipDir, zipDir},


	(* Messages *)
	questionOverwrite= pkgName<> " is already installed. Do you want to replace the content of "<> pkgDir<> " with a newly downloaded version?";

	(* Check Mathematica version *)
	If[$VersionNumber< minVersion,
		Message[InstallRGBeta::version, pkgName, ToString@ minVersion, $VersionNumber];
	];

	(* Check if SuperTracer has already been installed *)
	If[
		DirectoryQ[pkgDir],
		If[
			ChoiceDialog[questionOverwrite, {"Yes"-> True, "No"-> False},
                WindowFloating-> True,
                WindowTitle-> pkgName<> " installation detected"],
			Quiet@ DeleteDirectory[pkgDir, DeleteContents-> True],
			Abort[]
		];
	];

	(* Download SuperTracer *)
	Print["Downloading "<> pkgName<> " from ", pkgLink<> "."];

	tmpFile= Quiet@ URLSave[pkgLink];

	If[tmpFile=== $Failed,
		Print["Failed to download "<> pkgName<> ".\nInstallation aborted!"];
		Abort[]
	];

	(* Unzip SuperTracer file *)
	Print["Extracting "<> pkgName<> " zip file."];

	unzipDir= tmpFile<>".dir";
	ExtractArchive[tmpFile, unzipDir];

	(* Move files to the Mathematica packages folder *)
	Print["Copying "<> pkgName<> " to "<> pkgDir<> "."];

	zipDir= FileNames["RGBeta.m", unzipDir, Infinity];
	CopyDirectory[DirectoryName[zipDir[[1]], 2], pkgDir];

	(* Delete the extracted archive *)
	Quiet@ DeleteDirectory[unzipDir, DeleteContents-> True];
	Print["Installation complete!"];
];

InstallRGBeta[];
