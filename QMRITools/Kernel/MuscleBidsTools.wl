(* ::Package:: *)

(* ::Title:: *)
(*QMRITools MuscleBidsTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`MuscleBidsTools`", Join[{"Developer`"}, 
	Complement[QMRITools`$Contexts, {"QMRITools`MuscleBidsTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


ImportJSON::usage = 
"ImportJSON[file] imports a JSON file as rawJSON."

GetJSONPosition::usage = 
"GetJSONPosition[{json..}, {{key, value}..}] gets the position from a list of JSON association lists where keys have the given value.
GetJSONPosition[{json..}, {{key, value}..}, sortKey] same but finally sorts the positions for the value of the sortKey."

MergeJSON::usage = 
"MergeJSON[{json..}] merges a list of JSON association lists where duplicate keys with the same values are removed and duplicate keys with different values are merged."


SelectSubjects::usage = 
"SelectSubjects[dir] selects the subjects in the given data directory which has a config file.
SelectSubjects[dir, Dynamic[sel]] does the same but stores the selected subjects in the variable sel."

ViewConfig::usage =
"ViewConfig[config] shows a config file for Muscle Bids processing.
ViewConfig[folder] same but imports the config file from the given folder.
ViewConfig[\"default\"] shows the default config."

ViewProtocolNames::usage =
"ViewProtocolNames[folder] shows the protocol names in the radData folder specified by the config file with
the series number and protocol name from the json files in the raw folders for each subject and session.
Duplicate protocol names are highlighted and can be shown only using the checkbox. The config file only needs to define the folders."

GetProtocolNames::usage =
"GetProtocolNames[config, subs] gets the list of protocol names in the specified subject folder."

CheckConfigLabels::usage =
"CheckConfigLabels[dir] checks for each subject and session if the protocol names that are defined as Label in the config
file in dir are present. Missing labels are highlighted and can be shown only using the checkbox.
Only subjects in the rawData folder are shown, so run it after BidsDcmToNii.
CheckConfigLabels[dir, subs] does the same for the subjects subs."

CheckCustomConfigs::usage =
"CheckCustomConfigs[dir] lists for each subject that has a custom config, the settings that differ from the config file in dir, one row per setting.
It also checks that the custom config was carried along by the steps: the copy in a later stage that is missing, unreadable or different from the first copy is highlighted.
The stages are checked from the source folder via rawData and derivedData to mergeData. Subjects without a custom config are not shown.
The option ProcessSubjects limits the subjects."

CheckBidsOutput::usage =
"CheckBidsOutput[dir, step] checks for each subject if the output expected from the config file in dir is present after the given step.
Currently step can be \"Convert\", which checks for each Label if the json and nii files are present in the rawData folder.
Missing files are highlighted and can be shown only using the checkbox."

MuscleBidsConfig::usage =
"MuscleBidsConfig[dir] shows the config of the config file in dir and its checks as tabs: ViewConfig, CheckCustomConfigs, CheckConfigLabels and ViewProtocolNames.
The data is scanned when the function is called, the Reload buttons scan again. The option ProcessSubjects is passed on."

MuscleBidsCheck::usage =
"MuscleBidsCheck[dir] shows CheckBidsOutput and ManageCheckFiles for the config file in dir as tabs, the checks of the processing steps.
The Report tab exports the problems of a step, or of all steps, as md or pdf report.
The data is scanned when the function is called, the Reload buttons scan again. The options ProcessSubjects and BidsOutputImages are passed on."

MuscleBidsReport::usage =
"MuscleBidsReport[dir, step, file] writes the problem report of the config file in dir to file, for the step or for all steps when step is \"All\".
The steps are \"Protocol names\", \"Config labels\", \"Convert\", \"Custom configs\", \"Process\", \"Merge\", \"Segment\", \"Tractography\" and \"Analysis\".
Problems are grouped by class: error, unreadable, missing, duplicate and outdated. The extension of file selects the format, \".md\" or \".pdf\".
The options are ProcessSubjects to limit the subjects and BidsOutputImages to define the expected Analysis images. Returns the full path of file."

ManageCheckFiles::usage =
"ManageCheckFiles[dir] shows per processing step the check files expected from the config file in dir and their status.
Selected check files can be deleted such that the next run of that step processes the data again. Reload scans again.
The options are ProcessSubjects to limit the subjects and BidsOutputImages to define the expected Analysis images."

MuscleBidsGUI::usage =
"MuscleBidsGUI[dir] shows a GUI to run the Muscle Bids processing steps for the config file in dir, with the options of each step and the subjects to process.
Each step runs in a separate kernel, nothing is printed in the notebook and the progress is shown from the log files.
The tabs Config and Check show MuscleBidsConfig and MuscleBidsCheck."


GetConfig::usage =
"GetConfig[folder] imports a Muscle Bids config file from the given folder.
GetConfig[folder, name] same but for a config file with the given name."

MergeConfig::usage =
"MergeConfig[assoc, replace] merges the replace association with the assoc association."


PartitionBidsName::usage = 
"PartitionBidsName[name] converts a Bids name to the Bids labels as an association, i.e. 
{\"sub\",\"ses\",\"vol\",\"stk\",\"rep\",\"chunk\",\"acq\",\"part\",\"type\",\"suf\"}."

PartitionBidsFolderName::usage = 
"PartitionBidsFolderName[fol] partitions the Bids folder and file name. It returns the bids root folder and the 
label parts using PartitionBidsName."

GenerateBidsName::usage = 
"GenerateBidsName[parts] generates a Bids file name from the Bids labels association which can be generated by
PartitionBidsName. The labels are joined in the order {\"sub\",\"ses\",\"vol\",\"stk\",\"rep\",\"chunk\",\"acq\",\"part\",\"type\",\"suf\"}."

GenerateBidsFolderName::usage =
"GenerateBidsFolderName[parts] generates a Bids folder name from the Bids labels association which can be
generated by PartitionBidsFolderName.
GenerateBidsFolderName[fol, parts] the same but with a custom root folder."

GenerateBidsFileName::usage = 
"GenerateBidsFileName[parts] generates a Bids file name from the Bids labels association which can be generated 
by PartitionBidsName.
GenerateBidsFileName[fol, parts] the same but with a custom root folder."

SelectBidsFolders::usage =
"SelectBidsFolders[fol, tag] selects all folders in the fol with the name tag."

SelectBidsSubjects::usage =
"SelectBidsSubjects[fol] selects all subjects in the bids folder."

SelectBidsSessions::usage =
"SelectBidsSessions[fol] selects all sessions in the bids subject folder."


BidsDcmToNii::usage =
"BidsDcmToNii[dir] converts the bids sourceFolder with dicom files to raw nii files based on the config file."

MuscleBidsConvert::usage =
"MuscleBidsConvert[dir] converts all raw nii data in the to Muscle-Bids named nii based on the config 
file in the bids sourceFolder dir."

MuscleBidsProcess::usage = 
"MuscleBidsProcess[dir] processes all Muscle-Bids named nii based on the config file in the bids sourceFolder dir."

MuscleBidsMerge::usage = 
"MuscleBidsMerge[dir] merges multiple stack data for all Muscle-Bids named nii based on the config file in the 
bids sourceFolder dir."

MuscleBidsSegment::usage = 
"MuscleBidsSegment[dir] segments the data of Muscle-Bids named nii based on the config file in the bids 
sourceFolder dir. The segmentations are generated by the function SegmentData."

MuscleBidsTractography::usage =
"MuscleBidsTractography[dir] performs tractography on the Muscle-Bids named nii based on the config 
file in the bids sourceFolder dir. If a segmentation is present it is used as a mask for per muscle segmentation."

MuscleBidsAnalysis::usage = 
"MuscleBidsAnalysis[dir] performs analysis on the Muscle-Bids named nii based on the config file in the bids 
sourceFolder dir. If a segmentation is present it is used to calculate the mean per segmentation."

CheckDataDescription::usage =
"CheckDataDescription[description, method] checks the data description config file used in BidsDcmToNii,
MuscleBidsConvert, MuscleBidsProcess, and MuscleBidsMerge."



(* ::Subsection::Closed:: *)
(*Options*)


BidsIncludeSession::usage = 
"BidsIncludeSession is an option for BidsDcmToNii. If True session folders will be used in output even if 
they are not specified."

DeleteAfterConversion::usage = 
"DeleteAfterConversion is an option for MuscleBidsConvert. If set True all files that have been converted 
will be deleted."

BidsTractographyMethod::usage =
"BidsTractographyMethod is an option for MuscleBidsTractography and can be \"Full\", \"Tractography\" 
or \"Segmentation\". With \"Tractography\" the tractography is performed without segmentation.
With \"Segmentation\" only the segmentation is performed without tractography. With Full both are performed."

BidsOutputImages::usage = 
"BidsOutputImages is an option for MuscleBidsAnalysis and ManageCheckFiles. It sets which images are saved in the output folder,
one of \"All\", \"Quantitative\", \"Segmentation\", \"Tractography\" or \"None\"."

ProcessSubjects::usage = 
"ProcessSubjects is an option for Bids functions. Can be a list of bids subject names else it is All."

VersionCheck::usage = 
"VersionCheck is an option for all Bids functions. If set True data processed with an old version is reprocessed."


(* ::Subsection::Closed:: *)
(*Error Messages*)


CheckDataDescription::key = "Datasets have duplicate names which is not allowed."

CheckDataDescription::type = "Unknown Muscle-BIDS type: `1`, using folder \"miss\"."

CheckDataDescription::class = "Unknown Muscle-BIDS Class: `1`. Must be \"Volume\", \"Stacks\", 
\"Repetitions\", \"Chunks\", \"Acquisitions\"."

CheckDataDescription::lab = "Invalid combination of Class and Label: `1` with `2` is not allowed."

CheckDataDescription::man = "Mandatory values `1` are not in the data description."

CheckDataDescription::stk = "Class \"stacks\" or \"Chunk\" is used but overlap is not defined, assuming overlap 0."

GetConfig::conf = "Could not find config file in given folder."

MuscleBidsReport::dir = "The study folder `1` does not exist."

MuscleBidsReport::step = "Unknown step: `1`. Use a Muscle-BIDS step name or \"All\"."

MuscleBidsReport::format = "Unsupported file format, use .md or .pdf."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*General Definitions*)


(* ::Subsubsection::Closed:: *)
(*Debug and logging*)


debugBids[x___] := DebugPrint["Bids", x];


dataToLog[data_] := If[KeyExistsQ[data, $Failed], 
    "Wrong data description: " <> data[$Failed], 
    StringRiffle[KeyValueMap[ToString[#1] <> ": " <> ToString[#2] &, 
        KeyDrop[data, {"Process", "Merging", "Segment", "Tractography"}]], "; "]
]


compress = ($OperatingSystem === "Windows");


(* ::Subsubsection::Closed:: *)
(*Bids Definitions*)


bidsTypes = <|
	(*anatomy types*)
	"T1w"->"anat", "T1w-FS"->"anat", "T2w"->"anat", "T2w-FS"->"anat",
	(*dixon*)
	"megre"->"dix", "tse"->"dix",
	(*quant types*)
	"mese"->"quant",
	"T1"->"quant", "T2"->"quant", "wT2"->"quant",
	(*diff types*)
	"dwi"->"dwi",
	(*seg types*)
	"seg"->"seg"
|>;


bidsClass = {"Volume", "Volumes", "Stacks", "Repetitions", "Chunks", "Acquisitions", "Mixed"};


(* ::Subsubsection::Closed:: *)
(*Helper Functions*)


WipStrip = StringDelete[#, {"WIP ", "WIP_", "wip ", "wip_"}, IgnoreCase -> True] &;


StringStrip = StringDelete[#, {"-", "_", ".", " "}] &;


(* ::Subsubsection::Closed:: *)
(*BidsType*)


BidsType[type_?StringQ] := Lookup[bidsTypes, type, "miss"]

BidsType[parts_?AssociationQ] := Lookup[bidsTypes, parts["type"], "miss"]


(* ::Subsubsection::Closed:: *)
(*BidsValue*)


BidsValue[parts_, val_?ListQ] := Flatten[BidsValue[parts, #] &/@ val]

BidsValue[parts_, val_?StringQ] := Lookup[parts, val, Nothing]


(* ::Subsubsection::Closed:: *)
(*BidsString*)


BidsString[parts_, val_?ListQ] := BidsString[parts, #] &/@ val

BidsString[parts_, val_?StringQ] := Block[{str},
	str = BidsValue[parts, val];
	If[StringQ[str], val<>"-"<>str, str]
]


(* ::Subsubsection::Closed:: *)
(*GetClassName*)


GetClassName[class_, nameIn_] := Switch[class,
	"Volume"|"Volumes"|"Stacks"|"Repetitions"|"Acquisitions"|"Chunks",
	Switch[class, 
		"Volume"|"Volumes", "vol",
		"Stacks", "stk",
		"Chunks", "chunk", 
		"Repetitions", "rep", 
		"Acquisitions", "acq"
	] -> StringStrip[nameIn]
]


(* ::Subsubsection::Closed:: *)
(*SubNameToBids*)


Options[SubNameToBids] = {BidsIncludeSession -> True};

SubNameToBids[nameIn_?ListQ, opts : OptionsPattern[]] := SubNameToBids[#, "", opts] & /@ nameIn

SubNameToBids[nameIn_?ListQ, met_, opts : OptionsPattern[]] := SubNameToBids[#, met, opts] & /@ nameIn

SubNameToBids[nameIn_?StringQ, opts : OptionsPattern[]] := SubNameToBids[nameIn, "", opts]

SubNameToBids[nameIn_?StringQ, met_, OptionsPattern[]] := Block[{ass, keys, name, ses},
	(*get the names*)
	ass = Switch[met, 
		"Sub", PartitionBidsName, 
		"BidsDcmToNii", PartitionBidsName[FileNameTake[#, {2, -1}]]&, 
		_, PartitionBidsFolderName[#][[-1]]&
	]@nameIn;
	keys = Keys[ass];

	(*if bids take sub key else assume first suf is name*)
	name = "sub" -> If[MemberQ[keys, "sub"], ass["sub"], First[ass["suf"]]];

	(*if bids take ses key else assume last suf is session*)
	ses = "ses" -> If[MemberQ[keys, "ses"],
		(*session is present take session*)
		ass["ses"],
		(*more than one suf last is session*)
		If[Length[ass["suf"]] > 1, Last[ass["suf"]],
			(*no session,see if need to be forced*)
			If[OptionValue[BidsIncludeSession], "001", ""]]];
	Association[{name, ses, "suf" -> {}}]
]


(* ::Subsection:: *)
(*BIDS name and select*)


(* ::Subsubsection::Closed:: *)
(*PartitionBidsName*)


SyntaxInformation[PartitionBidsName] = {"ArgumentsPattern" -> {_}};

PartitionBidsName[list:{_?StringQ ..}] := PartitionBidsName/@list

PartitionBidsName[string_?StringQ] := Block[{parts, entity, suffix, suf},
	(*split entities and suffixes*)
	parts = StringSplit[string, "_"];
	entity = StringCases[#, k__ ~~ "-" ~~ v__ :> k -> v] & /@ parts;
	suf = Select[parts, ! StringContainsQ[#, Flatten@Keys@entity] &];

	(*see if type is part of suffixes*)
	suf = Which[
		suf==={}, {"parts"->{}},
		MemberQ[Keys[bidsTypes], First@suf], {"type"->First@suf,"suf"->Rest@suf},
		True, {"suf"->suf}
	];

	(*give as association*)
	Association[Join[entity, suf]]
]


(* ::Subsubsection::Closed:: *)
(*PartitionBidsFolderName*)


SyntaxInformation[PartitionBidsFolderName] = {"ArgumentsPattern" -> {_}};

PartitionBidsFolderName[fol:{_?StringQ ..}] := PartitionBidsFolderName/@fol

PartitionBidsFolderName[fol_?StringQ] := {
	First@StringSplit[fol, "sub-"], 
	PartitionBidsName@StringJoin@Riffle[Select[FileNameSplit[fol],StringContainsQ[#,"-"]&],"_"]
}


(* ::Subsubsection::Closed:: *)
(*GenerateBidsName*)


SyntaxInformation[GenerateBidsName] = {"ArgumentsPattern" -> {_}};

GenerateBidsName[list:{_?AssociationQ ..}] := GenerateBidsName/@list

GenerateBidsName[parts_?AssociationQ] := StringRiffle[Join[
	BidsString[parts, {"sub", "ses", "vol", "stk", "rep", "chunk", "acq", "part"}], 
	DeleteCases[BidsValue[parts, {"type", "suf"}],""]
], "_"]


(* ::Subsubsection::Closed:: *)
(*GenerateBidsFolderName*)


SyntaxInformation[GenerateBidsFolderName] = {"ArgumentsPattern" -> {_, _.}};

GenerateBidsFolderName[parts_?AssociationQ]:=GenerateBidsFolderName["", parts]

GenerateBidsFolderName[fol_?StringQ, list:{_?AssociationQ ..}] := GenerateBidsFolderName[fol,#]& /@ list

GenerateBidsFolderName[{fol_?StringQ, parts_?AssociationQ}] := GenerateBidsFolderName[fol, parts]

GenerateBidsFolderName[fol_?StringQ, parts_?AssociationQ] := FileNameJoin[{
	fol, BidsString[parts, "sub"], BidsString[parts, "ses"]
}]


(* ::Subsubsection::Closed:: *)
(*GenerateBidsFileName*)


SyntaxInformation[GenerateBidsFileName] = {"ArgumentsPattern" -> {_, _.}};

GenerateBidsFileName[list:{_?AssociationQ ..}] := GenerateBidsFileName["", #]& /@ list

GenerateBidsFileName[fol_?StringQ, list:{_?AssociationQ ..}] := GenerateBidsFileName[fol, #]& /@ list

GenerateBidsFileName[parts_?AssociationQ] := GenerateBidsFileName["", parts]

GenerateBidsFileName[{fol_?StringQ, parts_?AssociationQ}] := GenerateBidsFileName[fol, parts]

GenerateBidsFileName[fol_?StringQ, parts_?AssociationQ] := FileNameJoin[{
	GenerateBidsFolderName[fol, parts], BidsType[parts], GenerateBidsName[parts]
}]


(* ::Subsection:: *)
(*SelectBids*)


(* ::Subsubsection::Closed:: *)
(*SelectBids*)


SelectBids[folder_?ListQ, entity_?StringQ] := Catenate[SelectBids[#, entity]& /@ folder]

SelectBids[folder_?StringQ, entity_?StringQ] := Block[{
		baseName, start, end, what
	},
	(*see if the folder is a sub- or ses- folder*)
	baseName = StringStartsQ[FileNameTake[folder], #]&;

	(*levels below the folder that are still to be searched: sub-, ses-, entity*)
	start = Which[baseName["ses-"], 3, baseName["sub-"], 2, True, 1];
	end = Switch[entity, "sub", 1, "ses", 2, _, 3];
	what = {"sub-", "ses-", entity}[[start;;end]];

	(*step down level by level and keep the folders*)
	Select[Catenate[{Fold[FileNames[#2 <> "*", #1, {1}] &, folder, what]}], DirectoryQ]
]


(* ::Subsubsection::Closed:: *)
(*FilterBidsSubjects*)


FilterBidsSubjects[fols_, subs_] := FilterBidsSubjects[fols, subs, ""]

FilterBidsSubjects[fols_, subs_, met_] := If[subs === All || subs === "All", fols,
	Select[fols, MemberQ[SubNameToBids[Flatten[{subs}], "Sub"], SubNameToBids[#, met]]&]
]


(* ::Subsubsection::Closed:: *)
(*SelectBidsFolders*)


SyntaxInformation[SelectBidsFolders] = {"ArgumentsPattern" -> {_, _}};

SelectBidsFolders[fol_, tag_] := SelectBids[fol, tag]


(* ::Subsubsection::Closed:: *)
(*SelectBidsSubjects*)


SyntaxInformation[SelectBidsSubjects] = {"ArgumentsPattern" -> {_}};

SelectBidsSubjects[fol_] := SelectBids[fol, "sub"]


(* ::Subsubsection::Closed:: *)
(*SelectBidsSessions*)


SyntaxInformation[SelectBidsSessions] = {"ArgumentsPattern" -> {_}};

SelectBidsSessions[fol_] := SelectBids[fol, "ses"]


(* ::Subsection:: *)
(*Config*)


(* ::Subsubsection::Closed:: *)
(*Default config values*)


defaultConfig = <|
	"folders" -> <|
		"dicomData" -> "01_sourcedata", 
		"rawData" -> "02_rawdata", 
		"derivedData" -> "03_derivatives",
		"mergeData" -> "04_merged", 
		"analysis" -> "05_analysis"
	|>,
	"conversion"-> <|
		"Version"-> If[$OperatingSystem === "MacOSX", 1, "Own"]
	|>,
	"Process" -> <|
		"Masking" -> 5,
		"FlipPermute" -> {{1, 1, 1}, {"x", "y", "z"}},
		"SplitRegistration" -> False,
		"RegistrationDimension" -> "3D",
		"RegistrationMethod" -> "bspline",
		"RicianCorrection" -> True,
		"IVIMCorrection" -> True,
		"GradientCorrection" -> False,
		"FasciculationDetection" -> False
	|>,
	"Merging" -> <|
		"Overlap" -> 0,
		"Padding" -> 0,
		"Motion" -> True,
		"Reverse" -> False,
		"SplitRegistration" -> True
	|>,
	"Segment" -> <|
		"Device" -> "CPU",
		"VoxSize" -> Automatic,
		"Dimensions" -> "3D",
		"Location" -> "Body",
		"Method" -> Automatic
	|>,
	"Tractography" -><|
		"FlipPermute" -> {{1, 1, 1}, {"x", "y", "z"}},
		"TractLength"-> {15, 500},
		"TractAngle" -> 25,
		"TractSeed" -> 0.66,
		"TractStep" -> "Automatic",
		"SegmentLength" -> {15, 500},
		"BoneLabel" -> 100,
		"HarmonicDenoise" -> False
	|>,
	"Analysis" -> <|
		"MaskErosion" -> True,
		"TractWeighting" -> False,
		"UseFilter" -> "None"
	|>
|>;


(* ::Subsubsection::Closed:: *)
(*ConfigLookup*)


ConfigLookup[config_?StringQ, key_?StringQ] := ConfigLookup[GetConfig[config], key]

ConfigLookup[config_?AssociationQ, key_?StringQ] := Lookup[config, key, defaultConfig[key]]

ConfigLookup[config_?StringQ, key_?StringQ, def_?StringQ] := ConfigLookup[GetConfig[config], key, def]

ConfigLookup[config_?AssociationQ, key_?StringQ, def_?StringQ] := Lookup[ConfigLookup[config, key], def, defaultConfig[key, def]];


(* ::Subsubsection::Closed:: *)
(*CheckConfig*)


CheckConfig[fol_?StringQ] := CheckConfig[fol, fol, False]

CheckConfig[inFol_?StringQ, outFol_?StringQ] := CheckConfig[inFol, outFol, True]

CheckConfig[inFol_?StringQ, outFol_?StringQ, export_?BooleanQ] := Block[{conf, nam, path},
	(*subject-session name and the config file name made from it*)
	nam = GenerateBidsName[Last@PartitionBidsFolderName[outFol]];
	path = FileNameJoin[{outFol, nam<>"_config.json"}];
	(*find the config for this subject, a general or a subject specific one*)
	conf = Quiet@GetConfig[inFol, nam];

	debugBids[path];

	(*optionally write the config to the output folder, return if it was found*)
	If[conf =!= $Failed,
		If[export, Export[path, conf]]; {True, conf},
		{False, <||>}
	]
]


(* ::Subsubsection::Closed:: *)
(*BidsConfigFile*)


BidsConfigFile[folder_, name_] := SelectFirst[{
	folder,
	FileNameJoin[{folder, "config.json"}],
	FileNameJoin[{folder, name <> "_config.json"}]
}, FileExistsQ[#] && ! DirectoryQ[#] &, Missing[]]


(* ::Subsubsection::Closed:: *)
(*GetConfig*)


GetConfig[folder_?StringQ] := GetConfig[folder, ""]

GetConfig[folder_?StringQ, name_?StringQ] := Block[{file},
	(*find the config file in the folder*)
	file = BidsConfigFile[folder, name];
	(*import it or give a message*)
	If[!MissingQ[file],
		Import[file, "RawJSON"],
		Message[GetConfig::conf]; $Failed
	]
]


(* ::Subsubsection::Closed:: *)
(*BidsFolder*)


BidsFolder[dir_?StringQ, key_?StringQ] := BidsFolder[Replace[Quiet@GetConfig[dir], Except[_Association] -> <||>], dir, key]

BidsFolder[config_?AssociationQ, dir_?StringQ, key_?StringQ] := FileNameJoin[{dir, ConfigLookup[config, "folders", key]}]

BidsConfigQ[dir_] := AssociationQ[dir] || AssociationQ[Quiet@GetConfig[dir]]


(* ::Subsubsection::Closed:: *)
(*MergeConfig*)


MergeConfig[assoc_?ListQ, replace_?AssociationQ] := Normal@MergeConfig[Association@assoc, replace]

MergeConfig[assoc_?AssociationQ, replace_?ListQ] := MergeConfig[assoc, Association@replace]

MergeConfig[assoc_?AssociationQ, replace_?AssociationQ] := Merge[{assoc, replace}, 
	Replace[{
		(*Both are Association, go deeper*)
		{b_Association, p_Association} :> MergeConfig[b, p],
		(*Overwrite with patch if not both Associations*)
		{_, p_} :> p,
		(*Key exists in only one Association*)
		{v_} :> v 
	}]
]


(* ::Subsubsection::Closed:: *)
(*BuildBidsNameFromConfig*)


BuildBidsNameFromConfig[{fol_, parts_}, datType_] := BuildBidsNameFromConfig[{fol, parts}, datType, ""]

BuildBidsNameFromConfig[{fol_, parts_}, datType_, con_] := Block[{key, movStacks},
	(*Filename for the dataType it self*)
	
	If[con ==="",
		(*If no contrast generate generic output file name*)
		key = If[!datType["HasDuplicate"], Nothing,
			If[datType["Class"] === "Stacks", "stk", "chunk"] -> StringStrip@datType["Key"]];
		GenerateBidsFileName[fol, <|
			parts, key, "type" -> datType["Type"], "suf" -> datType["Suffix"]
		|>],

		(*if contrast specified generate the input file stack*)
		movStacks = StringStrip /@ Flatten[{datType["Label"]}];
		If[datType["Class"] === "Volumes", movStacks = movStacks[[{1}]]];
		(*output the stack file names*)
		GenerateBidsFileName[fol, <|parts, 
			GetClassName[datType["Class"], StringStrip[#]], 
			"type" -> datType["Type"], "suf" -> {datType["Suffix"], con}|>
		] <> ".nii" & /@ movStacks
	]
]


BuildBidsNameFromConfig[{folOut_, parts_}, {datType_, all_}, target_?ListQ] := BuildBidsNameFromConfig[{folOut, "", parts}, {datType, all}, target]

BuildBidsNameFromConfig[{folOut_, folIn_, parts_}, {datType_, all_}, tar_?ListQ] := Block[{
		target, isKey, key, type, suf, datTar, needKey, keyOut, tarType
	},
	(*Filename for remote target specified as list*)
	target = DeleteDuplicates[tar];

	(*see if target description has key and split target*)
	isKey = MemberQ[all[[All, "Key"]], First[target]];
	{key, type} = If[isKey, {First[target], Rest[target]}, {datType["Key"], target}];
	
	{type, suf} = {First[type], Rest[type]};

	(*make correct output key, quits if key is needed but not present*)
	keyOut = If[!datType["HasDuplicate"], Nothing,
		Switch[datType["Class"],
			"Stacks", "stk", 
			(*"Volume", "vol",*)
			_, "chunk"
		] -> StringStrip@key];

	(*If folder in is not specified only the output file is needed*)
	If[folIn === "", Return[GenerateBidsFileName[folOut, 
		<|parts, keyOut, "type" -> type, "suf" -> suf|>]<>".nii"]];

	(*if folIn is defined look for input stack file names*)
	tarType = First@If[isKey,
		Select[all, #["Key"] === key &],
		Select[all, #["InFolder"] === First[suf] &]
	];

	(*output the stack file names*)
	GenerateBidsFileName[folIn, <|parts, 
		GetClassName[tarType["Class"], StringStrip[#]], "type" -> type, "suf" -> suf|>
	] <> ".nii" & /@ Flatten[{tarType["Label"]}]
]


(* ::Subsection:: *)
(*JSON*)


(* ::Subsubsection::Closed:: *)
(*ImportJSON*)


SyntaxInformation[ImportJSON] = {"ArgumentsPattern" -> {_}};

ImportJSON[file_] := Replace[Quiet[Developer`ReadRawJSONFile[file]], Except[_Association | _List] :> Import[file, "RawJSON"]]


(* ::Subsubsection::Closed:: *)
(*GetJSONPosition*)


GetJSONPosition[json_, selection_] := GetJSONPosition[json, selection, ""]

GetJSONPosition[json_, selection_, sort_] := Block[{selIndex, selFunc, list, key, val, pos},
	(*selection functions*)
	selIndex = WipStrip[ToLowerCase[Last[Flatten[{#1 /. #3}]]]] === ToLowerCase[#2] &;
	selFunc = (list=#1; key=#2[[1]]; val=#2[[2]]; Select[list, selIndex[key, val, json[[#]]]&])&;

	(*get the file positions*)
	pos = Fold[selFunc, Range[Length[json]], selection];
	(*if sort is empty return the positions or if pos is empty return {}*)
	pos = If[sort===""||pos==={}, pos, pos[[Ordering[sort /. json[[pos]]]]]]
]


(* ::Subsubsection::Closed:: *)
(*MergeJSON*)


MergeJSON[json:{_?AssociationQ..}] := Block[{keys},
	keys = DeleteDuplicates[Flatten[Keys /@ json]];
	Association[If[#[[2]]==={}, Nothing, #]& /@ Thread[
		keys->(If[Length[#]===1,First@#,#]& /@ ((DeleteDuplicates /@ Transpose[(# /@ keys)& /@ json]) /. Missing[___]->Nothing))
	]]
]


(* ::Subsubsection::Closed:: *)
(*BidsDcmConsistentQ*)


(*dcm2niix builds that write one file per echo in the same layout for every Philips release: the own build (dcm2niix-own.exe)
and the official release that contains it (dcm2niix-20261029.exe) and later ones, not the official builds in between*)
bidsDcmOwn = 20260915;

bidsDcmConsistent = 20261029;

BidsDcmDate[json_] := Block[{date},
	(*build date in the dcm2niix version of the json, 0 when there is none*)
	date = StringCases[ToString[Lookup[json, "ConversionSoftwareVersion", ""]], "v1.0." ~~ d:Repeated[DigitCharacter, {8}] :> ToExpression[d]];
	If[date === {}, 0, First[date]]
]

BidsDcmConsistentQ[json_List] := AllTrue[json, (BidsDcmDate[#] === bidsDcmOwn || BidsDcmDate[#] >= bidsDcmConsistent)&]


(* ::Subsection:: *)
(*Bids Support Functions*)


(* ::Subsubsection::Closed:: *)
(*ViewConfig*)


SyntaxInformation[ViewConfig] = {"ArgumentsPattern" -> {_}};

ViewConfig["default"] := ViewConfig[defaultConfig]

ViewConfig[folder_?StringQ] := ViewConfig[Join[GetConfig[folder], <|"default"->defaultConfig|>]]

ViewConfig[config_?AssociationQ] := If[Length[config]>10, MenuView, TabView][# -> Which[
	# === "datasets"||# === "default", ViewConfig[config[#]],
	# === "analysis", If[KeyExistsQ[config["analysis"], "Analysis"], MakeTable[config[#]], ViewConfig[config[#]]], 
	True, MakeTable[config[#]]
]& /@ Keys[config], ControlPlacement -> Left]


MakeTable[association_] := Block[{key, value, bCol},
	bCol = LightDarkV[White, Darker@Gray];
	Grid[KeyValueMap[{
		Style[key = #1,Bold], 
		value = #2;
		Which[
			AssociationQ[value], MakeTable[value],
			ListQ[value] && Length[value] > 0, Which[
				AssociationQ[First[value]], Column[MakeTable /@ value],
				MatchQ[value, {__String}] && Length[value] > 6, Column[Row[#, ", "]&/@Partition[value, 6, 6, 1,{}]],
				VectorQ[value], Row[value, ", "], 
				key==="QuantImages", TabView[#[[1]] -> MakeTable[Association[#[[2]]]] & /@ ProcessImage[value], ControlPlacement -> Left],
				True, Column@value
			], 
			True, value
		]
	} &, association], 
	Frame -> All, Alignment -> Left, Background -> {{Gray, {bCol}}, bCol}, Spacings -> {1, 0.5}]
]


ProcessImage[an_] := Normal[GroupBy[an, Most[#[[1]]] &, <|#[[1, -1]] -> <|
	"Color map [Range]" -> StringRiffle@Flatten@{If[Length[#] >= 2, #[[2]], "Automatic"], "[", If[Length[#] >= 3, #[[3]], "Automatic"], "]"},
	"Label" -> StringRiffle@Flatten@If[Length[#] >= 4, {#[[4, 1]], "  ",#[[4, -1]], "  ", #[[4, 2]]}, {"Automatic"}]|>
|> & /@ # &]]


(* ::Subsubsection::Closed:: *)
(*ViewProtocolNames*)


Options[ViewProtocolNames] = {
	ProcessSubjects->All
}

SyntaxInformation[ViewProtocolNames] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ViewProtocolNames[dir:(_?StringQ | _?AssociationQ), OptionsPattern[]] :=
	BidsCheckView[dir, OptionValue[ProcessSubjects], "Protocol", {"Series", "Protocol name", "Duplicate"}]


(* ::Subsubsection::Closed:: *)
(*GetProtocolNames*)


GetProtocolNames[fold_?ListQ] := GetProtocolNames /@ fold

GetProtocolNames[fold_?StringQ] := Block[{list, json, name},
	(*series number and cleaned protocol name of every json in the folder*)
	list = Sort@DeleteDuplicates[(
		Quiet[json = ImportJSON[#]];
		{
			json["SeriesNumber"], 
			name = json["ProtocolName"]; 
			If[StringQ[name], StringTrim[WipStrip@name], name]
		}
	) & /@ FileNames["*.json", fold, 2]]; 

	(*keep valid entries and skip the names with a date and time stamp*)
	Select[list,
		Head[#[[1]]] =!= Missing && Head[#[[1]]] =!= $Failed && StringQ[#[[2]]] &&
		!StringContainsQ[#[[2]], RegularExpression["_\\d{6}\\.\\d{3}"]] &
	]
]


GetProtocolNames[dir_, sub_] := Block[{config},
	config = If[AssociationQ[dir], dir, GetConfig[dir]];
	GetProtocolNames[FilterBidsSubjects[SelectBids[ConfigLookup[config, "folders", "rawData"], "ses"], sub]][[If[StringQ[sub], 1, All]]]
]


(* ::Subsubsection::Closed:: *)
(*BidsScan*)


bidsTableSteps = {"Protocol names", "Config labels", "Convert", "Custom configs"};

bidsCheckSteps = {"Process", "Merge", "Segment", "Tractography", "Analysis"};

bidsReportSteps = Join[bidsTableSteps, bidsCheckSteps];

bidsReportClasses = {"error", "unreadable", "missing", "duplicate", "outdated"};


$BidsScan = <||>;

(*a lazy view has no scan result until Reload, set by the GUI for its tabs*)
$BidsLazy = False;

BidsJsonFiles[fol_, labels_] := Flatten[FileNames["*"<>StringStrip[#]<>"*.json", fol]& /@ Flatten[{labels}]]

BidsProtocolNames[fol_] := Block[{names},
	(*scan the folder only if it is not in the cache, is cleared by Reload and by a new view*)
	names = Lookup[$BidsScan, fol, Missing[]];
	If[MissingQ[names],
		names = GetProtocolNames[fol];
		AssociateTo[$BidsScan, fol -> names]
	];
	names
]


(* ::Subsubsection::Closed:: *)
(*SelectSubjects*)


SyntaxInformation[SelectSubjects] = {"ArgumentsPattern" -> {_, _.}};

SelectSubjects[dir_?StringQ] := DynamicModule[{sel = {}}, SelectSubjects[dir, Dynamic[sel]]]

SelectSubjects[dir_?StringQ, Dynamic[sel_]] := DynamicModule[{list},
	(*subject folders of the dicom folder*)
	list = Sort[FileBaseName/@ Select[FileNames["*", BidsFolder[dir, "dicomData"]], DirectoryQ]];
	(*checkbox list with buttons to select all, none or copy*)
	Column[{
		CheckboxBar[Dynamic[sel], list, Appearance -> "Vertical" -> {Min[{15, Length@list}], Automatic}, Method -> "Active"],
		Row[{
			Button["Select All", sel = list],
			Button["Deselect All", sel = {}],
			Button["Copy selected list to clipboard", CopyToClipboard[sel]]
		}]
	}]
]


(* ::Subsubsection::Closed:: *)
(*CheckConfigLabels*)


SyntaxInformation[CheckConfigLabels] = {"ArgumentsPattern" -> {_, _.}};

CheckConfigLabels[dir_?StringQ] := CheckConfigLabels[dir, All]

CheckConfigLabels[dir_?StringQ, subs_] := BidsCheckView[dir, subs, "Labels", {"Dataset", "Label", "Present"}]


(* ::Subsubsection::Closed:: *)
(*BidsSubjectList*)


BidsSubjectList[dir_, subs_] := Block[{config, root, dis},
	(*config from association or folder, default if invalid*)
	{config, root} = If[AssociationQ[dir], {dir, Directory[]}, {Quiet@GetConfig[dir], dir}];
	config = If[AssociationQ[config], config, <||>];
	
	(*per subject folder the merged datasets*)
	{config, {#, Last@PartitionBidsFolderName[#],
		dis = Quiet@CheckDataDescription[MergeConfig[Lookup[config, "datasets", <||>], Last@CheckConfig[#]], "MuscleBidsProcess"];
		If[ListQ[dis], Select[dis, !KeyExistsQ[#, $Failed]&], {}]
	}& /@ FilterBidsSubjects[SelectBids[BidsFolder[config, root, "rawData"], "ses"], subs]}
]


(* ::Subsubsection::Closed:: *)
(*BidsTable*)


BidsTable[rows_, head_] := With[{hrows = If[MatrixQ[head], head, {head}]}, Pane[
	Grid[Join[Map[Style[#, Bold]&, hrows, {2}], First /@ rows],
		Alignment -> Left, Spacings -> {1.5, 0.6}, ItemSize -> {Automatic, 1.6},
		(*bad rows get a light red row, the header a grey one*)
		Background -> {None, Join[ConstantArray[LightDarkV[GrayLevel[.9], GrayLevel[.26]], Length[hrows]], MapIndexed[Which[
					Last[#1], LightDarkV[RGBColor[1, .9, .9], RGBColor[.38, .15, .15]],
					EvenQ[First[#2]], LightDarkV[GrayLevel[.96], GrayLevel[.19]],
					True, LightDarkV[White, GrayLevel[.13]]
				]&, rows]]
			},
		Dividers -> {None, All}, FrameStyle -> LightDarkV[GrayLevel[.8], GrayLevel[.32]]
	],
	
	(*scroll when more than 15 rows*)
	If[Length[rows] > 15, {Automatic, 430}, Automatic], 
	Scrollbars -> Automatic, ImageSizeAction -> "Scrollable", 
	Background -> LightDarkV[White, GrayLevel[.13]]
]]


(* ::Subsubsection::Closed:: *)
(*BidsCheckView*)


BidsCheckView[dir_?StringQ, subs_, what_, head_] /; ExpandFileName[dir] =!= dir := BidsCheckView[ExpandFileName[dir], subs, what, head]

BidsCheckView[dir_, subs_, what_, head_] := With[{
		data = If[$BidsLazy, None, BidsCheckRows[dir, subs, what]],
		valid0 = BidsConfigQ[dir],
		t0 = If[$BidsLazy, "not scanned yet", DateString["Time"]]
	},
	DynamicModule[{rows, valid, time, only, show, cur, bad},
		(*scan result and time, a lazy view (GUI tab) has no scan result until Reload*)
		rows = data;
		valid = valid0;
		time = t0;
		(*the filter (on for all but the protocol names) and the selected subject*)
		only = !MemberQ[{"Protocol", "Custom"}, what];
		cur = None;

		Column[{
			Row[{
				(*rescan with an empty protocol names cache*)
				Button["Reload",
					$BidsScan = <||>;
					rows = BidsCheckRows[dir, subs, what];
					valid = BidsConfigQ[dir];
					time = DateString["Time"]
				, Method -> "Queued"], "  ",
				Checkbox[Dynamic[only]], " only problems  ", 
				Dynamic["scanned at "<>time]
			}],
			Dynamic[
				(*only the bad rows and their subjects*)
				show = If[rows === None, None, If[only, Select[(First[#] -> Select[Last[#], Last])& /@ rows, Last[#] =!= {}&], rows]];
				If[show === {} || show === None, 
					Which[
						show === None, "Not scanned yet, press Reload",
						only, "No problems found",
						!valid, "No subjects found, the config file is missing or invalid",
						what === "Custom", "No custom configs found",
						True, "No subjects found"
					],
					MenuView[
						Function[{sub, list},
							(*menu entry with number of problems*)
							bad = Count[list[[All, 2]], True];
							If[bad > 0, Style[sub<>" ("<>ToString[bad]<>")", Red], sub] ->
							BidsTable[list, head]
						] @@@ show,
						(*keep selected subject*)
						Dynamic[
							Replace[FirstPosition[show[[All, 1]], cur], {{i_} :> i, _ -> 1}], 
							(cur = show[[#, 1]])&
						],
						ControlPlacement -> Top
					]
				]
			]
		}]
	]
]


BidsCheckRows[dir_, subs_, what_] := Block[{config, sets, rows},
	(*a new protocol names scan starts with an empty cache*)
	If[what === "Protocol", $BidsScan = <||>];
	(*config and the subjects with their datasets*)
	{config, sets} = BidsSubjectList[dir, subs];
	(*rows per subject, the custom configs need the whole config*)
	rows = (
		GenerateBidsName[#[[2]]] -> If[what === "Custom",
			BidsCustomRows[config, dir, #],
			BidsCheckRows[#, what]
		]
	)& /@ sets;
	(*custom configs only keep the subjects with a row*)
	If[what === "Custom", DeleteCases[rows, _ -> {}], rows]
]

BidsCheckRows[{fol_, _, datDis_}, what_] := Block[{list, dup, names, labs, files},
	(*Labels of each dataset*)
	labs = Flatten[{Lookup[#, "Label", {}]}]& /@ datDis;

	Switch[what,
		(*-------------- protocol names --------------*)
		"Protocol",
		list = BidsProtocolNames[fol];
		(*duplicate protocol names are bad*)
		dup = Keys[Select[Counts[list[[All, 2]]], # > 1 &]];
		{
			Append[#, If[MemberQ[dup, #[[2]]], "yes", ""]], 
			MemberQ[dup, #[[2]]]
		} & /@ list
		,
		(*-------------- config labels --------------*)
		"Labels",
		names = BidsProtocolNames[fol][[All, 2]];
		(*per Label if it is in the protocol names*)
		list = MapThread[Function[{dat, labels},
			{dat["Key"], #, MemberQ[names, #]}& /@ labels
		], {datDis, labs}];
		BidsPresentRows@DeleteDuplicates@Flatten[list, 1]
		,
		(*-------------- converted output --------------*)
		"Convert",
		list = MapThread[Function[{dat, labels}, Table[
			(*per Label the converted json and if its nii exists*)
			files = BidsJsonFiles[FileNameJoin[{fol, dat["OutFolder"]}], lab];
			(*no json, an empty json (failed conversion) or a missing nii is bad*)
			If[files === {},
				{{dat["Key"], lab, "no json found", False}},
				{
					dat["Key"], lab, FileBaseName[#],
					Which[
						!NiiFileExistQ[#], False,
						FileByteCount[#] === 0, "json empty",
						True, True
					]
				}& /@ files
			]
			(*Volumes: first Label only*)
			, {lab, If[dat["Class"] === "Volumes", Take[labels, UpTo[1]], labels]}
		]], {datDis, labs}];
		BidsPresentRows@Flatten[list, 2]
	]
]


BidsPresentRows[rows_] := {
	Append[Most[#], Switch[Last[#], 
		True, "yes", 
		False, Style["missing", Red, Bold], 
		_, Style[Last[#], Red, Bold]
	]],	Last[#] =!= True
}& /@ rows


(* ::Subsubsection::Closed:: *)
(*BidsCustomConfig*)


BidsCustomStages[config_, dir_, parts_] := Block[{dicDir, srcNames, stages},
	(*dicom data folder, the plain subject folder is also tried for session 001*)
	dicDir = BidsFolder[config, dir, "dicomData"];
	srcNames = Join[{parts["sub"]<>"_"<>parts["ses"]}, If[parts["ses"] === "001", {parts["sub"]}, {}]];
	(*source stage first, then the raw, derived and merged stages*)
	stages = Join[
		{{ConfigLookup[config, "folders", "dicomData"], SelectFirst[FileNameJoin[{dicDir, #}]& /@ srcNames, DirectoryQ, Missing[]]}},
		({ConfigLookup[config, "folders", #], GenerateBidsFolderName[BidsFolder[config, dir, #], parts]})& /@ {"rawData", "derivedData", "mergeData"}
	];
	(*keep the stages whose subject folder exists*)
	Select[stages, DirectoryQ[#[[2]]]&]
]

BidsCustomRead[folder_, name_] := Block[{file, conf},
	(*config file of the stage, Missing when there is none*)
	file = BidsConfigFile[folder, name];
	(*parse quietly, $Failed unless it is an association*)
	conf = If[MissingQ[file], Missing[], Quiet@Import[file, "RawJSON"]];
	If[MissingQ[conf] || AssociationQ[conf], conf, $Failed]
]

BidsCustomLeaves[assoc_?AssociationQ, path_] := Join @@ KeyValueMap[BidsCustomLeaves[#2, Append[path, #1]]&, assoc]

BidsCustomLeaves[val_, path_] := {{path, val}}

BidsCustomText[val_] := Block[{str},
	str = TextString[val];
	If[StringLength[str] > 80, StringTake[str, 77]<>"...", str]
]

BidsCustomRow[datasets_, path_, val_] := Block[{dat, study},
	(*dataset name and the value at the same path in the study config*)
	dat = First[path];
	study = Fold[If[AssociationQ[#1] && KeyExistsQ[#1, #2], #1[#2], Missing[]]&, datasets, path];
	{
		If[KeyExistsQ[datasets, dat], dat, dat<>" (not in config)"],
		StringRiffle[Rest[path], "."],
		If[MissingQ[study], "", BidsCustomText[study]],
		BidsCustomText[val]
	}
]

BidsCustomSettings[custom_?AssociationQ, datasets_?AssociationQ] :=
	(*one row per leaf of the custom patch (the first key of its path is the dataset), sorted by dataset and setting*)
	Sort[BidsCustomRow[datasets, #1, #2]& @@@ BidsCustomLeaves[custom, {}]]

BidsCustomProblem[stage_, text_] := {{stage, "custom config", "", Style[text, Red, Bold]}, True}

BidsCustomNorm[assoc_?AssociationQ] := KeySort[BidsCustomNorm /@ assoc]

BidsCustomNorm[val_] := val

BidsCustomRows[config_, dir_, {_, parts_, _}] := Block[{name, stages, reads, k, ref, refName},
	(*config file name of the subject and the config read of each stage*)
	name = GenerateBidsName[parts];
	stages = BidsCustomStages[config, dir, parts];
	reads = BidsCustomRead[#[[2]], name]& /@ stages;
	(*the first stage with a config is the reference*)
	k = Length[TakeWhile[reads, MissingQ]] + 1;
	If[k > Length[stages], Return[{}, Block]];
	ref = reads[[k]];
	refName = stages[[k, 1]];
	If[ref === $Failed, Return[{BidsCustomProblem[refName, "unreadable"]}, Block]];
	(*settings of the reference, then the later stages without the config, unreadable or with other content*)
	Join[
		{#, False}& /@ BidsCustomSettings[ref, Replace[Lookup[config, "datasets", <||>], Except[_?AssociationQ] -> <||>]],
		BidsCustomProblem[#1, Which[
			MissingQ[#2], "not copied",
			#2 === $Failed, "unreadable",
			True, "differs from "<>refName
		]]& @@@ Select[Drop[Thread[{stages[[All, 1]], reads}], k], BidsCustomNorm[#[[2]]] =!= BidsCustomNorm[ref]&]
	]
]


(* ::Subsubsection::Closed:: *)
(*CheckBidsOutput*)


Options[CheckBidsOutput] = {ProcessSubjects -> All}

SyntaxInformation[CheckBidsOutput] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

CheckBidsOutput[dir_?StringQ, step:"Convert", OptionsPattern[]] :=
	BidsCheckView[dir, OptionValue[ProcessSubjects], step, {"Dataset", "Label", "File", "Present"}]


(* ::Subsubsection::Closed:: *)
(*CheckCustomConfigs*)


Options[CheckCustomConfigs] = {ProcessSubjects -> All}

SyntaxInformation[CheckCustomConfigs] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

CheckCustomConfigs[dir_?StringQ, OptionsPattern[]] :=
	BidsCheckView[dir, OptionValue[ProcessSubjects], "Custom", 
		{{"", "", "Normal config", "Custom config"}, {"Dataset", "Setting", "Value", "Value"}}]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsConfig*)


Options[MuscleBidsConfig] = {ProcessSubjects -> All}

SyntaxInformation[MuscleBidsConfig] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

MuscleBidsConfig[dir_?StringQ, OptionsPattern[]] := Block[{subs},
	subs = OptionValue[ProcessSubjects];
	(*one tab per config tool, each scans its data now unless the views are lazy*)
	TabView[{
		"Full config" -> If[BidsConfigQ[dir], ViewConfig[dir], "No valid config file found in the folder"],
		"Custom config" -> CheckCustomConfigs[dir, ProcessSubjects -> subs],
		"Config labels" -> CheckConfigLabels[dir, subs],
		"Protocol names" -> ViewProtocolNames[dir, ProcessSubjects -> subs]
	}, ControlPlacement -> Top]
]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsCheck*)


Options[MuscleBidsCheck] = {ProcessSubjects -> All, BidsOutputImages -> "All"}

SyntaxInformation[MuscleBidsCheck] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

MuscleBidsCheck[dir_?StringQ, OptionsPattern[]] := Block[{subs, imOut},
	{subs, imOut} = OptionValue[{ProcessSubjects, BidsOutputImages}];
	(*one tab per check tool, each scans its data now unless the views are lazy, the report scans on export*)
	TabView[{
		"Convert output" -> CheckBidsOutput[dir, "Convert", ProcessSubjects -> subs],
		"Check files" -> ManageCheckFiles[dir, ProcessSubjects -> subs, BidsOutputImages -> imOut],
		"Report" -> BidsReportView[dir, subs, imOut]
	}, ControlPlacement -> Top]
]


(* ::Subsubsection::Closed:: *)
(*ManageCheckFiles*)


Options[ManageCheckFiles] = {ProcessSubjects -> All, BidsOutputImages -> "All"}

SyntaxInformation[ManageCheckFiles] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ManageCheckFiles[dir_?StringQ, OptionsPattern[]] := BidsCheckFilesView[dir, OptionValue[ProcessSubjects], OptionValue[BidsOutputImages]]


BidsCheckFilesView[dir_?StringQ, subs_, imOut_] /; ExpandFileName[dir] =!= dir := BidsCheckFilesView[ExpandFileName[dir], subs, imOut]

BidsCheckFilesView[dir_, subs_, imOut_] := With[{
		scan0 = If[$BidsLazy, {None, <||>}, BidsScanChecks[dir, subs, imOut]],
		t0 = If[$BidsLazy, "not scanned yet", DateString["Time"]]
	},
	DynamicModule[{checks, stat, time, sel, upd, ind, sub, scan, exist, row, delFiles, buttons, tab},
		(*scan result and time, a lazy view (GUI tab) has no scan result until Reload*)
		{checks, stat} = scan0;
		time = t0;
		(*the state of the tabs and the selection*)
		sel = {};
		upd = 0;
		ind = 1;
		sub = <||>;

		(*scan the expected check files and read the status of each*)
		scan = (
			{checks, stat} = BidsScanChecks[dir, subs, imOut];
			time = DateString["Time"];
			upd++
		)&;

		(*check file exists, its status was read during the scan*)
		exist = (Lookup[stat, #, None] =!= None)&;

		(*cells of one check file and if it is bad: missing, unreadable or error*)
		row = Function[{nam, base}, Block[{st, cls, vals},
			st = Lookup[stat, base, None];
			cls = BidsCheckClass[st];
			vals = If[AssociationQ[st], Lookup[st, {"Check", "Version", "Date"}, ""], {"", "", ""}];
			(*bad status in red*)
			If[cls =!= "", vals[[1]] = Style[cls, Red, Bold]];
			{Join[{
				If[exist[base],
					Checkbox[Dynamic[
						MemberQ[sel, base],
						(sel = If[#, Union[sel, {base}], DeleteCases[sel, base]])&
					]],
					""
				],
				nam, BidsCheckName[nam, base]
			}, vals], cls =!= ""}
		]];

		(*confirm with the file list and delete the selected check files of one tab*)
		delFiles = Function[{step, list}, Block[{del, msg},
			del = Intersection[sel, list[[All, 2]]];
			msg = Column[{
				"Delete "<>ToString[Length[del]]<>" check files? The next run of "<>step<>" processes them again.",
				Pane[Column[FileNameTake /@ del], ImageSize -> {UpTo[600], UpTo[300]}, Scrollbars -> Automatic]
			}];
			If[del =!= {} && ChoiceDialog[msg],
				DeleteFile[#<>"_check.json"]& /@ del;
				stat = Join[stat, AssociationThread[del -> None]];
				sel = Complement[sel, del];
				upd++
			]
		]];

		(*buttons act on the check files of one tab*)
		buttons = Function[{step, list}, Row[{
			Button["Select all", sel = Union[sel, Select[list[[All, 2]], exist]]],
			Button["Deselect all", sel = Complement[sel, list[[All, 2]]]],
			Button["Delete selected", delFiles[step, list], Method -> "Queued"],
			(*Reload scans the expected check files again, Refresh status only reads the known ones*)
			Button["Reload", scan[], Method -> "Queued"],
			Button["Refresh status",
				stat = BidsCheckStatus[checks];
				time = DateString["Time"];
				upd++
			, Method -> "Queued"],
			Dynamic["  "<>ToString[Length[Intersection[sel, list[[All, 2]]]]]<>" selected, scanned at "<>time]
		}]];

		(*buttons above the table of one tab*)
		tab = Column[{
			buttons[#1, #2],
			Dynamic[upd; BidsTable[row @@@ #2, {"", "Subject", "Check file", "Status", "Version", "Date"}]]
		}]&;

		(*tabs per step and data type*)
		Dynamic[
			If[checks === None,
				Row[{Button["Reload", scan[], Method -> "Queued"], "  Not scanned yet, press Reload"}],
				TabView[KeyValueMap[Function[{step, list}, step -> If[list === {}, "No check files expected",
					TabView[KeyValueMap[#1 -> tab[step, #2]&, GroupBy[list, First -> Rest]],
						Dynamic[Lookup[sub, step, 1], (sub = Append[sub, step -> #])&], ControlPlacement -> Top]
				]], checks], Dynamic[ind], ControlPlacement -> Top]
			]
		]
	]
]


BidsScanChecks[dir_, subs_, imOut_] := Block[{checks},
	(*expected check files and the status of each*)
	checks = BidsCheckFiles[dir, subs, imOut];
	{checks, BidsCheckStatus[checks]}
]


BidsCheckStatus[checks_] := Block[{bases},
	(*file base of every expected check file*)
	bases = DeleteDuplicates[Last /@ Catenate[Values[checks]]];
	(*status of each*)
	AssociationThread[bases, BidsReadCheck /@ bases]
]

BidsReadCheck[base_] := Block[{file, json},
	(*None is missing, $Failed is unreadable, else the content*)
	file = base<>"_check.json";
	If[!FileExistsQ[file], Return[None, Block]];
	json = Quiet@ImportJSON[file];
	If[AssociationQ[json], json, $Failed]
]


BidsCheckClass[st_] := Which[
	(*problem class of a check file status: missing, unreadable or error, "" when ok*)
	st === None, "missing",
	!AssociationQ[st], "unreadable",
	Lookup[st, "Check", ""] === "error", "error",
	True, ""
]

BidsCheckName[subj_, base_] := StringReplace[
	(*check file name without the subject prefix*)
	FileNameTake[base], StartOfString ~~ subj ~~ "_" -> ""
]


BidsCheckFiles[dir_, subs_, imOut_] := Block[{config, sets, der, mer, ana, anDis},
	{config, sets} = BidsSubjectList[dir, subs];
	(*output folders and analysis datasets*)
	{der, mer, ana} = BidsFolder[config, dir, #]& /@ {
			"derivedData", "mergeData", "analysis"
		};
	anDis = If[KeyExistsQ[config, "analysis"], 
		Quiet@CheckDataDescription[config["analysis"], "MuscleBidsAnalysis"], 
		{}
	];
	anDis = If[ListQ[anDis], 
		Select[anDis, !KeyExistsQ[#, $Failed]&], 
		{}
	];

	(*expected check files of all subjects grouped per step*)
	AssociationThread[bidsCheckSteps, Lookup[GroupBy[
		Flatten[BidsCheckFiles[{der, mer, ana}, #, anDis, imOut]& /@ sets, 1], First -> Rest
	], bidsCheckSteps, {}]]
]


BidsCheckFiles[{der_, mer_, ana_}, {fol_, parts_, datDis_}, anDis_, imOut_] := Block[{
		nam, typ, files, sets, base, partsO, dixMet
	},
	nam = GenerateBidsName[parts];
	dixMet = {"Dixon-S", "Dixon-P", "Dixon", "Dixon-B", "Dixon-A"};
	Join[
		Flatten[Table[
			typ = dat["Type"];
			(*raw input files of the dataset as sets*)
			files = BidsJsonFiles[FileNameJoin[{fol, dat["InFolder"]}], Lookup[dat, "Label", {}]];
			sets = If[files === {}, {}, PartitionBidsName[FileBaseName /@ files]];
			sets = DeleteDuplicates@If[MemberQ[{"megre", "tse"}, typ],
				Append[#, "suf" -> {First[#["suf"]]}]& /@ sets,
				sets
			];
			base = GenerateBidsFileName[der, #]& /@ sets;
			Join[
				(*process, dwi prep and fit if the method is processed*)
				If[typ === "dwi",
					{"Process", typ, nam, #<>"_prep"}& /@ base,
					{}
				],
				If[MemberQ[Lookup[<|"megre" -> dixMet, "tse" -> dixMet, "dwi" -> {"DTI"}, "mese" -> {"EPGT2", "Exp"}|>, 
						typ, {}], dat["Process", "Method"]],
					{"Process", typ, nam, #}& /@ base,
					{}
				],
				(*merge, segment and tractography output*)
				If[KeyExistsQ[dat, "Merging"],
					{{"Merge", typ, nam, BuildBidsNameFromConfig[{mer, parts}, dat]}},
					{}
				],
				If[KeyExistsQ[dat, "Segment"],
					{{"Segment", typ, nam,
						ConvertExtension@BuildBidsNameFromConfig[{mer, parts}, {dat, datDis}, 
							{dat["Key"], "seg", "auto", typ}]
					}},
					{}
				],
				If[KeyExistsQ[dat, "Tractography"],
					{{"Tractography", typ, nam,
						ConvertExtension@BuildBidsNameFromConfig[{mer, parts}, {dat, datDis}, 
							Join[dat["Tractography", "Target"][[;;-2]], {"trk"}]]
					}},
					{}
				]
			]
		, {dat, datDis}], 1](*close dataset loop*),
		(*analysis xls and images per analysis dataset*)
		Flatten[(
			typ = #["Key"];
			partsO = If[typ === "Default",
				parts,
				Join[<|If[#["Analysis", "Class"] === "Stacks", "stk", "chunk"] -> StringStrip@typ|>, parts]
			];
			{
				"Analysis", If[typ === "Default", #, typ<>" "<>#], nam,
				FileNameJoin[{GenerateBidsFolderName[ana, partsO], GenerateBidsName[<|partsO, "suf" -> {#}|>]}]
			}& /@ If[imOut === "None", {"xls"}, {"xls", "img"}]
		)& /@ anDis, 1]
	]
]


(* ::Subsubsection::Closed:: *)
(*BidsReportData*)


BidsReportScope[step_] := If[step === "All", bidsReportSteps, {step}]


BidsReportData[dir_, subs_, imOut_, step_] := Block[{
		config, sets, configOK, datasets, names, scope, checked,
		tables, checkSteps, custom, tabRows, checks, stat, checkRows, rows
	},
	(*config and subject folders of the study*)
	{config, sets} = BidsSubjectList[dir, subs];
	(*config validity and the largest number of datasets of a subject*)
	configOK = BidsConfigQ[dir];
	datasets = Max[Prepend[Length[Last[#]] & /@ sets, 0]];
	(*dataset names of all subjects*)
	names = DeleteDuplicates[Flatten[Lookup[Last[#], "Key", {}]& /@ sets]];

	(*steps to check, without a valid config only the protocol names*)
	scope = BidsReportScope[step];
	checked = If[configOK,
		scope,
		Select[scope, # === "Protocol names" &]
	];
	tables = Select[checked, MemberQ[bidsTableSteps, #] &];
	checkSteps = Select[checked, MemberQ[bidsCheckSteps, #] &];
	(*reset the protocol names cache*)
	$BidsScan = <||>;
	(*custom config scan of the subjects, once for the step*)
	custom = If[MemberQ[tables, "Custom configs"],
		(GenerateBidsName[#[[2]]] -> BidsCustomRows[config, dir, #])& /@ sets,
		{}
	];
	(*problem rows of the table steps*)
	tabRows = Association[(# -> BidsReportRows[#, sets, dir, custom])& /@ tables];
	(*expected check files of the check steps and their status*)
	checks = If[checkSteps === {},
		<||>,
		KeyTake[BidsCheckFiles[dir, subs, imOut], checkSteps]
	];
	stat = BidsCheckStatus[checks];
	(*problem rows of the check steps*)
	checkRows = Association[
		Function[st,
			st -> Catenate[BidsReportCheckRow[st, #, stat, dir]& /@ checks[st]]
		] /@ checkSteps
	];

	(*problem rows per step, grouped by class in the fixed class order*)
	rows = Map[KeyTake[GroupBy[#, First -> Last], bidsReportClasses]&, Join[tabRows, checkRows]];
	<|
		"Subjects" -> Length[sets],
		"Datasets" -> datasets,
		"DatasetNames" -> names,
		"ConfigOK" -> configOK,
		"Checked" -> checked,
		"Info" -> <|"Custom configs" -> BidsCustomInfo[custom]|>,
		"CustomSubjects" -> Length[Select[custom, Length[Last[#]] > 0 &]],
		"Rows" -> rows
	|>
]


BidsReportRows[step_, sets_, dir_, custom_] := Block[{what},
	(*custom configs come from the scan of the subjects*)
	If[step === "Custom configs",
		Return[BidsCustomItems[custom], Block]
	];
	(*code of the table step in BidsCheckRows*)
	what = Switch[step, "Protocol names", "Protocol", "Config labels", "Labels", "Convert", "Convert"];

	(*problem rows of each subject, with the session folder of its dataset entry*)
	Catenate@Table[
		BidsReportRow[what, GenerateBidsName[entry[[2]]], #[[1]], entry[[1]], dir]& /@
			Select[BidsCheckRows[entry, what], Last]
	, {entry, sets}]
]


BidsCustomItems[custom_] := Catenate@Table[
	BidsCustomItem[pair[[1]], #]& /@ Select[pair[[2]], Last]
, {pair, custom}]


BidsCustomInfo[custom_] := Catenate@Table[
	BidsCustomInfoRow[pair[[1]], #[[1]]]& /@ Select[pair[[2]], !Last[#] &]
, {pair, custom}]


BidsCustomItem[subj_, {cells_, True}] := Block[{plain, text, cls},
	(*plain text cells, without Style*)
	plain = cells /. Style[x_, ___] :> x;
	text = plain[[4]];

	(*unreadable stage or a stage that differs from the reference*)
	cls = If[text === "unreadable", "unreadable", "outdated"];
	BidsReportItem[cls, subj, plain[[1]], "custom config", "", <|"Detail" -> text|>]
]


BidsCustomInfoRow[subj_, cells_] := <|
	"Subject" -> subj,
	"Dataset" -> cells[[1]],
	"Setting" -> cells[[2]],
	"Study value" -> cells[[3]],
	"Custom value" -> cells[[4]]
|>


BidsReportItem[cls_, subj_, type_, name_, log_, extra_] := {
	cls, 
	Join[<|"Subject" -> subj, "Type" -> type, "Name" -> name, "Log" -> log|>, extra]
}


BidsReportRow[what_, subj_, rowCells_, folder_, dir_] := Block[{cells},
	(*plain text cells, without Style*)
	cells = rowCells /. Style[x_, ___] :> x;

	(*class, data type, name and the extra values of the row*)
	Switch[what,
		"Protocol", BidsReportItem["duplicate", subj, "", cells[[2]], "", <|"Series" -> ToString[cells[[1]]]|>],
		"Labels", BidsReportItem["missing", subj, cells[[1]], cells[[2]], "", <||>],
		"Convert", BidsReportItem["missing", subj, cells[[1]], cells[[2]], BidsLogFile[dir, folder, subj, "Convert"],
			<|"Detail" -> Which[
				cells[[3]] === "no json found", "no json found",
				cells[[4]] === "json empty", "json empty: "<>cells[[3]],
				True, "nii missing: "<>cells[[3]]
			]|>
		]
	]
]


BidsReportCheckRow[step_, {type_, subj_, base_}, stat_, dir_] := Block[{st, cls, date},
	(*status of the check file and its problem class*)
	st = Lookup[stat, base, None];
	cls = BidsCheckClass[st];
	(*date of an error check file, empty for the other classes*)
	date = If[cls === "error", Lookup[st, "Date", ""], ""];

	(*only problem files get a row*)
	If[cls === "",
		{},
		{BidsReportItem[
			cls, subj, type, BidsCheckName[subj, base],
			BidsLogFile[dir, BidsSessionFolder[base, subj], subj, step],
			<|"Date" -> date|>
		]}
	]
]


BidsSessionFolder[base_, subj_] := Block[{parts, pos},
	(*folder of the subject session above the file, "" when there is none*)
	parts = FileNameSplit[base];

	(*last path part named as the subject session*)
	pos = Flatten[Position[parts, Last[StringSplit[subj, "_"]]]];
	If[pos === {}, "", FileNameJoin[Take[parts, Last[pos]]]]
]


BidsLogFile[dir_, folder_, subj_, step_] := Block[{file},
	(*log file of the step for the subject, "" when the folder is unknown*)
	file = If[StringQ[folder] && folder =!= "", FileNameJoin[{folder, subj<>"_BIDS"<>step<>".log"}], ""];

	(*path relative to the study folder when the file exists, else ""*)
	If[FileExistsQ[file], FileNameDrop[ExpandFileName[file], FileNameDepth[ExpandFileName[dir]]], ""]
]


BidsReportSection[level_, title_, step_, cls_, rows_] := <|
	"Level" -> level, 
	"Title" -> title, 
	"Step" -> step, 
	"Class" -> cls,
	"Columns" -> If[rows === None, None, BidsReportColumns[step, cls]],
	"Rows" -> rows
|>


BidsReportAllSections[byStep_, cls_] := Block[{steps, n},
	(*steps with rows of the class and the total of the class*)
	steps = Select[Keys[byStep], KeyExistsQ[byStep[#], cls] &];
	n = Total[Length[byStep[#][cls]] & /@ steps];

	(*class heading, then one table per step*)
	Prepend[
		BidsReportSection[3, #<>" ("<>ToString[Length[byStep[#][cls]]]<>")", #, cls, byStep[#][cls]] & /@ steps,
		BidsReportSection[2, cls<>" ("<>ToString[n]<>")", "", cls, None]
	]
]


BidsReportStepSections[step_, rowsByClass_, classes_] :=
	(*step heading, then one table per class*)
	Prepend[
		BidsReportSection[3, #<>" ("<>ToString[Length[rowsByClass[#]]]<>")", step, #, rowsByClass[#]] & /@ classes,
		BidsReportSection[2, step, step, "", None]
	]


BidsReportLayout[data_, step_] := Block[{
		head, scope, checked, unchecked, byStep, classes, summary, sections, infoRows, infoSections, fields, notes
	},
	head = data["Head"];
	(*steps in the report, the one step or all of them*)
	scope = BidsReportScope[step];
	(*checked steps are summarised, the other steps only get a note*)
	checked = Select[scope, MemberQ[data["Checked"], #] &];
	unchecked = Select[scope, !MemberQ[data["Checked"], #] &];
	(*problem rows of each checked step*)
	byStep = Association[(# -> Lookup[data["Rows"], #, <||>])& /@ checked];
	(*classes that occur, in fixed order*)
	classes = Select[bidsReportClasses, MemberQ[Catenate[Keys /@ Values[byStep]], #] &];
	
	(*summary: column names and one row per checked step with its count per class and in total*)
	summary = If[classes === {},
		<|"Title" -> "Summary", "Head" -> {}, "Rows" -> {}|>,
		<|
			"Title" -> "Summary",
			"Head" -> Join[{"Step"}, classes, {"Total"}],
			"Rows" -> KeyValueMap[
				Join[
					{#1},
					Length /@ Lookup[#2, classes, {}],
					{Total[Length /@ Values[#2]]}
				]&,
				byStep
			]
		|>
	];

	(*sections: class then step for all steps, step then class for one step*)
	sections = If[classes === {},
		{},
		If[step === "All",
			Catenate[BidsReportAllSections[byStep, #] & /@ classes],
			BidsReportStepSections[step, byStep[step], classes]
		]
	];

	(*info section of the custom config settings, not a problem*)
	infoRows = Lookup[data["Info"], "Custom configs", {}];
	infoSections = If[infoRows === {},
		{},
		{BidsReportSection[2, "Custom configs ("<>ToString[Length[infoRows]]<>" settings)", "Custom configs", "info", infoRows]}
	];

	(*header fields, the custom configs only when they are checked*)
	fields = Join[
		{
			"Study" -> head["Study"],
			"Date" -> head["Date"],
			"QMRITools version" -> head["Version"],
			"Subjects" -> head["Subjects"],
			"Datasets" -> If[head["DatasetNames"] === {}, "none", StringRiffle[head["DatasetNames"], ", "]]
		},
		If[MemberQ[checked, "Custom configs"],
			{"Custom configs" -> ToString[data["CustomSubjects"]]<>" of "<>ToString[head["Subjects"]]<>" subjects"},
			{}
		]
	];

	(*notes: config problem, unchecked steps, and the reason when nothing is listed*)
	notes = Join[
		If[head["ConfigOK"] === False,
			{"The config file is missing or invalid, the default folders were used."},
			{}
		],
		If[unchecked === {},
			{},
			{"Not checked: the config file is missing or invalid. Steps: "<>StringRiffle[unchecked, ", "]<>"."}
		],
		Which[
			classes =!= {}, {},
			head["Subjects"] === 0, {"Not checked: no subjects found."},
			head["Datasets"] === 0, {"Not checked: no datasets in the config."},
			True, {"No problems found."}
		]
	];

	<|
		"Title" -> "Muscle-BIDS problem report",
		"Fields" -> fields,
		"Notes" -> notes,
		"Summary" -> summary,
		"Sections" -> Join[sections, infoSections]
	|>
]


BidsReportMarkdown[layout_] := Block[{lines},
	(*title, fields and notes, then the summary and the sections*)
	lines = Join[
		BidsReportHeader[layout],
		BidsReportSummary[layout["Summary"]],
		Catenate[BidsReportTable /@ layout["Sections"]]
	];
	StringRiffle[lines, "\n"]
]


BidsReportHeader[layout_] :=
	(*title, one line per field, then the notes as quotes*)
	Join[
		{"# "<>layout["Title"], ""},
		("- "<>ToString[First[#]]<>": "<>BidsReportCell[ToString[Last[#]], True]) & /@ layout["Fields"],
		{""},
		Catenate[({"> "<>#, ""}) & /@ layout["Notes"]]
	]


BidsReportSummary[summary_] :=
	(*no lines when the summary has no rows*)
	If[summary["Rows"] === {}, {},
		{
			"## "<>summary["Title"], 
			"", 
			BidsReportGrid[summary["Head"], 
			Map[BidsReportCell[#, False]&, summary["Rows"], {2}]], 
			""
		}
	]


BidsReportTable[section_] := Block[{heading, cols, cellRow, grid},
	(*heading of the section*)
	heading = StringRepeat["#", section["Level"]]<>" "<>section["Title"];

	(*heading-only section has no table*)
	If[section["Rows"] === None,
		{heading, ""},
		cols = section["Columns"];
		cellRow = Function[row,
			MapThread[BidsReportCell[Lookup[row, #1, ""], #2]&, {cols["Keys"], cols["Code"]}]
		];
		grid = BidsReportGrid[cols["Head"], cellRow /@ section["Rows"]];
		{heading, "", grid, ""}
	]
]


BidsReportColumns[step_, cls_] := Block[{spec, head, keys, codes},
	(*column names, keys of the row and the keys shown as code*)
	spec = Switch[step,
		(*-------------- protocol names --------------*)
		"Protocol names", {
			{"Subject", "Series", "Protocol name"}, 
			{"Subject", "Series", "Name"}, 
			{}
		},
		(*-------------- config labels --------------*)
		"Config labels", {
			{"Subject", "Dataset", "Label"},
			{"Subject", "Type", "Name"},
			{}
		},
		(*-------------- converted output --------------*)
		"Convert", {
			{"Subject", "Dataset", "Label", "Detail", "Log"},
			{"Subject", "Type", "Name", "Detail", "Log"},
			{"Log"}
		},
		(*-------------- custom configs, problems and settings --------------*)
		"Custom configs", If[cls === "info",
			{
				{{"", "", "", "Normal config", "Custom config"}, {"Subject", "Dataset", "Setting", "Value", "Value"}},
				{"Subject", "Dataset", "Setting", "Study value", "Custom value"},
				{"Setting"}
			},
			{
				{"Subject", "Stage", "Problem"},
				{"Subject", "Type", "Detail"},
				{}
			}
		],
		(*-------------- check files, error with its date --------------*)
		_, If[cls === "error",
			{
				{"Subject", "Type", "Check file", "Date", "Log"},
				{"Subject", "Type", "Name", "Date", "Log"},
				{"Name", "Log"}
			},
			{
				{"Subject", "Type", "Check file", "Log"},
				{"Subject", "Type", "Name", "Log"},
				{"Name", "Log"}
			}
		]
	];
	{head, keys, codes} = spec;
	<|
		"Head" -> head,
		"Keys" -> keys,
		"Code" -> (MemberQ[codes, #] & /@ keys)
	|>
]


BidsReportCell[val_, isCode_] := Block[{txt},
	(*text of the cell on one line*)
	txt = StringReplace[ToString[val], "\n" -> " "];
	
	(*code cells: backtick becomes a quote and pipes are escaped, other cells escape the markdown characters*)
	If[isCode,
		txt = StringReplace[txt, {"`" -> "'", "|" -> "\\|"}],
		txt = StringReplace[txt, {
			"\\" -> "\\\\", "`" -> "\\`", "*" -> "\\*", "_" -> "\\_", "<" -> "\\<",
			">" -> "\\>", "[" -> "\\[", "]" -> "\\]", "|" -> "\\|"
		}]
	];
	
	(*empty stays empty, code goes in backticks*)
	Which[txt === "", "", isCode, "`"<>txt<>"`", True, txt]
]


BidsReportGrid[head_?MatrixQ, rows_] := BidsReportGrid[
	(*markdown has one header row, the group label goes in front of the column name*)
	MapThread[If[#1 === "", #2, #1<>" "<>ToLowerCase[#2]]&, head], rows
]

BidsReportGrid[head_, rows_] := StringRiffle[Join[
	{"| "<>StringRiffle[head, " | "]<>" |", "|"<>StringRepeat["---|", Length[head]]},
	("| "<>StringRiffle[#, " | "]<>" |")& /@ rows
], "\n"]


(* ::Subsubsection::Closed:: *)
(*BidsReportPdf*)


BidsReportPdf[layout_, file_] := Block[{cells, nb},
	(*title, fields, notes, summary and the sections as cells*)
	cells = Join[
		{Cell[layout["Title"], "Title", FontFamily -> "Helvetica", FontSize -> 16, FontColor -> Black]},
		BidsReportGridCells[{}, List @@@ layout["Fields"], {False, False}],
		(Cell[#, "Text", FontSize -> 9, FontColor -> Black]& /@ layout["Notes"]),
		BidsReportSummaryCells[layout["Summary"]],
		Catenate[BidsReportSectionCells /@ layout["Sections"]]
	];
	(*margins of the cells*)
	cells = Append[Append[#, CellMargins -> {{0, 0}, {2, 2}}], LineSpacing -> {1, 0}]& /@ cells;

	(*A4 portrait with 50 pt margins and a white background*)
	nb = Notebook[cells,
		Background -> White,
		PageWidth -> 495,
		PrintingOptions -> {"PaperSize" -> {595, 842}, "PrintingMargins" -> {{50, 50}, {50, 50}}}
	];

	(*export the notebook as PDF*)
	Export[file, nb, "PDF"]
]


BidsReportSummaryCells[summary_] := Block[{head},
	(*no cells when the summary has no rows*)
	head = summary["Head"];
	If[summary["Rows"] === {},
		{},
		Prepend[
			BidsReportGridCells[head, summary["Rows"], ConstantArray[False, Length[head]]],
			Cell[summary["Title"], "Section", FontFamily -> "Helvetica", FontSize -> 12, FontColor -> Black]
		]
	]
]


BidsReportSectionCells[section_] := Block[{cols, body},
	(*columns of the table, none for a heading-only section*)
	cols = section["Columns"];
	(*cells of the rows, keyed by the column keys*)
	body = If[cols === None,
		{},
		BidsReportGridCells[
			cols["Head"],
			Map[Lookup[#, cols["Keys"], ""]&, section["Rows"]],
			cols["Code"]
		]
	];
	(*heading in the style of its level, then the table*)
	Prepend[body, Cell[
		section["Title"],
		Switch[section["Level"], 2, "Section", 3, "Subsection"],
		FontFamily -> "Helvetica",
		FontSize -> Switch[section["Level"], 2, 12, 3, 10],
		FontColor -> Black
	]]
]


BidsReportWrap[txt_String, n_] := Block[{tokens, lines},
	(*words split after a separator, words longer than a line are cut*)
	tokens = Catenate[
		If[StringLength[#] > n, StringJoin /@ Partition[Characters[#], UpTo[n]], {#}]& /@
			StringCases[txt, RegularExpression["[^\\\\/_ -]*[\\\\/_ -]|[^\\\\/_ -]+"]]
	];
	(*greedy fill of the lines*)
	lines = If[tokens === {},
		{""},
		Fold[
			If[StringLength[Last[#1]] + StringLength[#2] <= n,
				Append[Most[#1], Last[#1] <> #2],
				Append[#1, #2]
			]&,
			{First[tokens]},
			Rest[tokens]
		]
	];
	lines
]


BidsReportGridCells[head_, rows_, codes_] := Block[{
		hrows, all, body, adv, lens, fracs, chars, cellOf, rowOf, headBoxes, bgs, groups
	},
	(*header rows, none, one or a group row above the column names*)
	hrows = Which[head === {}, {}, MatrixQ[head], head, True, {head}];
	(*cell texts of the head and the rows on one line*)
	all = Map[StringReplace[ToString[#], "\n" -> " "]&, Join[hrows, rows], {2}];
	(*body rows without the head*)
	body = Drop[all, Length[hrows]];
	(*width of a character per column, wider for the monospace code columns*)
	adv = If[#, 5.6, 5.0]& /@ codes;
	(*longest text of each column in characters, at least 10 and at most 20 for code and 24 for text columns*)
	lens = MapThread[Max[10, Min[If[#2, 20, 24], Max[StringLength /@ #1]]]&, {Transpose[all], codes}];
	(*column fractions of the text width in proportion to the width of the longest text*)
	fracs = lens adv/Total[lens adv];
	(*characters per line of each column*)
	chars = Clip[Floor[(495 fracs - 6)/adv], {1, Infinity}];
	
	(*cell box of the wrapped lines as strings, monospace for code cells*)
	cellOf = Function[{txt, n, code, weight},
		(*one line is a plain string, more lines a grid without row spacing*)
		With[{lines = ToBoxes /@ BidsReportWrap[txt, n]},
			StyleBox[
				If[Length[lines] === 1,
					First[lines],
					GridBox[List /@ lines,
						GridBoxAlignment -> {"Columns" -> {{Left}}},
						GridBoxSpacings -> {"Rows" -> {{0}}}
					]
				],
				FontFamily -> If[code, "Courier", "Helvetica"],
				FontSize -> 8,
				FontColor -> Black,
				FontWeight -> weight,
				LineSpacing -> {1, 0},
				ShowStringCharacters -> False
			]
		]
	];

	(*one row of cell boxes with its font weight*)
	rowOf = Function[{cells, weight}, MapThread[cellOf[#1, #2, #3, weight]&, {cells, chars, codes}]];

	(*header row in bold*)
	headBoxes = rowOf[#, "Bold"]& /@ Take[all, Length[hrows]];
	(*background of the whole rows: grey header, then white and light grey alternating*)
	bgs = Table[If[EvenQ[i], GrayLevel[0.95], White], {i, Length[body]}];

	(*groups of 30 rows, one notebook cell each with the header, so that long tables paginate*)
	groups = Partition[Thread[{body, bgs}], UpTo[30]];
	Function[group,
		Cell[BoxData[GridBox[
			Join[headBoxes, (rowOf[#[[1]], "Plain"]&) /@ group],
			GridBoxAlignment -> {"Columns" -> {{Left}}, "Rows" -> {{Top}}},
			GridBoxSpacings -> {"Columns" -> {{1}}, "Rows" -> {{0.2}}},
			GridBoxBackground -> {"Rows" -> {Join[ConstantArray[GrayLevel[0.85], Length[hrows]], group[[All, 2]]]}},
			GridBoxItemSize -> {"Columns" -> {Scaled /@ fracs}}
		]], "Text"]
	] /@ groups
]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsReport*)


Options[MuscleBidsReport] = {ProcessSubjects -> All, BidsOutputImages -> "All"}

SyntaxInformation[MuscleBidsReport] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

MuscleBidsReport[dir_?StringQ, step_?StringQ, file_?StringQ, opts:OptionsPattern[]] /; ExpandFileName[dir] =!= dir := MuscleBidsReport[ExpandFileName[dir], step, file, opts]

MuscleBidsReport[dir_?StringQ, step_?StringQ, file_?StringQ, OptionsPattern[]] := Block[{
		subs, imOut, data, info, layout, path
	},
	(*study folder must exist*)
	If[!DirectoryQ[dir],
		Message[MuscleBidsReport::dir, dir];
		Return[$Failed, Block]
	];
	(*step must be a step name or All*)
	If[!MemberQ[Append[bidsReportSteps, "All"], step],
		Message[MuscleBidsReport::step, step];
		Return[$Failed, Block]
	];
	(*file must be .md or .pdf*)
	If[!MemberQ[{"md", "pdf"}, ToLowerCase[FileExtension[file]]],
		Message[MuscleBidsReport::format];
		Return[$Failed, Block]
	];
	(*subjects and output images from the options*)
	{subs, imOut} = OptionValue[{ProcessSubjects, BidsOutputImages}];
	(*problem rows of the step and the subject count*)
	data = BidsReportData[dir, subs, imOut, step];
	(*study details for the header*)
	info = <|
		"Study" -> dir,
		"Date" -> DateString["ISODate"],
		"Version" -> QMRITools`$InstalledVersion,
		"Subjects" -> data["Subjects"],
		"Datasets" -> data["Datasets"],
		"DatasetNames" -> data["DatasetNames"],
		"ConfigOK" -> data["ConfigOK"]
	|>;
	(*neutral layout of the report*)
	layout = BidsReportLayout[Append[data, "Head" -> info], step];
	(*full path of the file, folder created when missing*)
	path = ExpandFileName[file];
	If[!DirectoryQ[DirectoryName[path]], CreateDirectory[DirectoryName[path]]];
	(*write md as UTF-8 text, pdf as notebook*)
	Switch[ToLowerCase[FileExtension[file]],
		"md", Export[path, BidsReportMarkdown[layout], "Text", CharacterEncoding -> "UTF-8"],
		"pdf", BidsReportPdf[layout, path]
	]
]


(* ::Subsubsection::Closed:: *)
(*BidsReportView*)


BidsReportFile[dir_, folder_, step_, format_] := Block[{name, stepName},
	(*study folder name without drive and path characters*)
	name = StringDelete[FileNameTake[ExpandFileName[dir]], {":", "\\", "/"}];
	(*study when nothing is left*)
	If[name === "", name = "study"];
	(*step name, All gives all*)
	stepName = If[step === "All", "all", StringReplace[step, " " -> "_"]];
	(*report file in the folder with the date of today*)
	FileNameJoin[{folder, name<>"_report_"<>stepName<>"_"<>DateString["ISODate"]<>"."<>format}]
]


BidsReportExport[dir_, subs_, imOut_, folder_, step_, format_] := Block[{file},
	(*file name of the report*)
	file = BidsReportFile[dir, If[folder === "", dir, folder], step, format];
	(*write the report, returns the file or $Failed*)
	MuscleBidsReport[dir, step, file, ProcessSubjects -> subs, BidsOutputImages -> imOut]
]


BidsReportView[dir_?StringQ, subs_, imOut_] /; ExpandFileName[dir] =!= dir := BidsReportView[ExpandFileName[dir], subs, imOut]

BidsReportView[dir_, subs_, imOut_] := DynamicModule[{
		format, folder, status, last, file, pick, export, busy
	},
	(*output format, output folder, status line and the last exported file*)
	format = "md";
	folder = dir;
	status = "";
	last = None;
	busy = "Working ...";

	(*export button of one step, disabled while an export runs*)
	export = Function[{label, step}, Button[label,
		status = busy;
		file = CheckAbort[BidsReportExport[dir, subs, imOut, folder, step, format], $Failed];
		(*status and last file of the export*)
		If[StringQ[file],
			status = "Exported: "<>file;
			last = file
			,
			status = "Failed, see the messages"
		]
		, Method -> "Queued", ImageSize -> 200, Enabled -> Dynamic[status =!= busy]]];

	Column[{
		(*format and folder rows*)
		Row[{Style["Format", Bold], "  ", RadioButtonBar[Dynamic[format], {"md", "pdf"}]}],
		Row[{Style["Folder", Bold], "  ", InputField[Dynamic[folder], String, FieldSize -> 50], "  ",
			Button["Browse",
				(*pick the output folder*)
				pick = SystemDialogInput["Directory", folder];
				If[StringQ[pick], folder = pick]
				, Method -> "Queued"]}],
		(*export of all steps and of the steps in two groups*)
		export["Export all steps", "All"],
		Row[{
			Column[Prepend[
				export[#, #]& /@ {"Protocol names", "Config labels", "Convert", "Custom configs"},
				Style["Set-up checks", Bold]
			], Spacings -> 1],
			"    ",
			Column[Prepend[
				export[#, #]& /@ {"Process", "Merge", "Segment", "Tractography", "Analysis"},
				Style["Pipeline steps", Bold]
			], Spacings -> 1]
		}, Alignment -> Top],
		(*status line with the open button of the last file*)
		Row[{Dynamic[status], "  ",
			Button["Open", SystemOpen[last], Enabled -> Dynamic[last =!= None], Method -> "Queued"]
		}]
	}, Spacings -> 1]
]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsGUI*)


SyntaxInformation[MuscleBidsGUI] = {"ArgumentsPattern" -> {_}};

MuscleBidsGUI[dir_?StringQ] /; ExpandFileName[dir] =!= dir := MuscleBidsGUI[ExpandFileName[dir]]

MuscleBidsGUI[dir_?StringQ] := With[{
		config = Block[{$BidsLazy = True}, MuscleBidsConfig[dir]],
		check = Block[{$BidsLazy = True}, MuscleBidsCheck[dir]],
		(*step, folder of its logs and the options in the GUI*)
		steps = {
			{BidsDcmToNii, "rawData", {}},
			{MuscleBidsConvert, "rawData", {DeleteAfterConversion}},
			{MuscleBidsProcess, "derivedData", {VersionCheck}},
			{MuscleBidsMerge, "mergeData", {VersionCheck}},
			{MuscleBidsSegment, "mergeData", {VersionCheck}},
			{MuscleBidsTractography, "mergeData", {VersionCheck, BidsTractographyMethod}},
			{MuscleBidsAnalysis, "analysis", {VersionCheck, BidsOutputImages}}
		},
		optVals = <|
			BidsTractographyMethod -> {"Full", "Tractography", "Segmentation"},
			BidsOutputImages -> {"All", "Quantitative", "Segmentation", "Tractography", "None"}
		|>,
		msgFile = FileNameJoin[{$TemporaryDirectory, "MuscleBidsGUI_"<>IntegerString[Hash[dir], 16]<>"_"<>#<>".txt"}]&
	},
	DynamicModule[{
		link, t0, tick, state, step, root, now, total, touched, tail, msgs,
		all, sel, opt, ctrl, run, poll
	},

	(*no step running, empty status and all subjects selected*)
	link = None;
	t0 = 0;
	tick = 0;
	state = "Idle";
	step = "";
	root = "";
	now = "";
	total = 0;
	touched = 0;
	tail = {};
	msgs = {};
	all = True;
	sel = {};

	(*option values per step start at the defaults of the step*)
	opt = Association[Flatten[Function[row, ({SymbolName[row[[1]]], #} -> (# /. Options[row[[1]]]))& /@ row[[3]]] /@ steps, 1]];

	(*checkbox for True/False options, popup for the others*)
	ctrl = Row[{ToString[#2], ": ",
		If[KeyExistsQ[optVals, #2],
			PopupMenu[Dynamic[opt[{#1, #2}]], optVals[#2]],
			Checkbox[Dynamic[opt[{#1, #2}]]]
		], "   "}]&;

	(*start a step in a separate kernel*)
	run = Function[row, Which[
		link =!= None, Null,
		!all && sel === {}, now = "Select subjects or tick All subjects",
		True,
		step = SymbolName[row[[1]]];
		t0 = AbsoluteTime[];
		tail = {};
		now = "";
		touched = 0;
		msgs = {};
		root = BidsFolder[dir, row[[2]]];
		total = If[row[[1]] === BidsDcmToNii, 0,
			Length@FilterBidsSubjects[SelectBids[BidsFolder[dir, "rawData"], "ses"], If[all, All, sel]]];
		link = BidsStartStep[dir, row[[1]], Join[{ProcessSubjects -> If[all, All, sel]}, (# -> opt[{step, #}])& /@ row[[3]]], msgFile[step]];
		(*no link means the kernel did not start*)
		If[Head[link] === LinkObject,
			state = "Running",
			link = None;
			state = "Crashed";
			now = "Could not start the kernel"
		]
	]];

	(*read the state of the running step*)
	poll = (If[link =!= None,
		{state, now, touched, tail, msgs} = Lookup[BidsStepStatus[link, root, step, t0, msgFile[step]], 
			{"State", "Now", "Touched", "Tail", "Messages"}];
		tick++;
		If[state =!= "Running",
			Quiet@LinkClose[link];
			link = None;
			Quiet@DeleteFile[msgFile[step]]
		]
	])&;

	Column[{
		TabView[{
			"Run" -> Column[{
				Style["Subjects", Bold],
				Row[{Checkbox[Dynamic[all]], " All subjects"}],
				Dynamic[If[all, "", SelectSubjects[dir, Dynamic[sel]]]],
				Style["Steps", Bold],
				Grid[Function[row, {
					Button[SymbolName[row[[1]]], run[row], ImageSize -> 200, Enabled -> Dynamic[link === None], Method -> "Queued"],
					Row[ctrl[SymbolName[row[[1]]], #]& /@ row[[3]]]
				}] /@ steps, Alignment -> Left],
				Button["Stop",
					If[link =!= None && ChoiceDialog["Stop "<>step<>"? The subjects that are not finished are processed again the next time."],
						(*the step can finish while the dialog is open, then there is nothing to stop*)
						If[link =!= None,
							LinkClose[link];
							link = None;
							state = "Stopped"
						]
					],
					ImageSize -> 200, Enabled -> Dynamic[link =!= None], Method -> "Queued"],
				Style["Status", Bold],
				Row[{
					(*moving bar while a step runs*)
					Dynamic[If[state === "Running", ProgressIndicator[Appearance -> "Percolate", ImageSize -> {120, 12}], ""], TrackedSymbols :> {state}],
					Dynamic[tick; Style[If[step === "", "", step<>": "]<>state<>If[state === "Running", " ("<>ToString[Round[AbsoluteTime[] - t0]]<>" s)", ""], Bold,
						Sequence @@ Which[state === "Crashed", {Red}, state === "Done with messages", {Orange}, True, {}]]]
				}, "  "],
				Dynamic[tick; Column[{
					If[total > 0, "subjects started "<>ToString[touched]<>" of "<>ToString[total], Nothing],
					Style[now, Italic]
				}]],
				OpenerView[{"Log", Dynamic[Pane[Style[StringRiffle[tail, "\n"], FontFamily -> "Courier"], ImageSize -> {600, 200}, Scrollbars -> Automatic]]}],
				OpenerView[{"Messages", Dynamic[Pane[Style[StringRiffle[msgs, "\n"], FontFamily -> "Courier"], ImageSize -> {600, 100}, Scrollbars -> Automatic]]}]
			}, Spacings -> 1],
			"Config" -> config,
			"Check" -> check
		}, ControlPlacement -> Top],
		(*poll the running step, only while a link is open*)
		Dynamic[If[link =!= None, Refresh[poll[]; "", UpdateInterval -> 2, TrackedSymbols :> {}], ""], TrackedSymbols :> {link}, SynchronousUpdating -> False]
	}],

	Initialization :> (link = None; If[state === "Running", state = "Interrupted"]),
	UnsavedVariables :> {link},
	Deinitialization :> If[link =!= None, Quiet@LinkClose[link]]
]]


BidsStartStep[dir_, st_, opts_, msg_] := Block[{link, load},
	If[FileExistsQ[msg], DeleteFile[msg]];
	link = LinkLaunch[First[$CommandLine] <> " -wstp"];
	(*dev or release paclet as in LaunchTrainingKernels*)
	load = If[StringContainsQ[First[PacletFind["QMRITools"]]["Location"], "workspace"], 
		"QMRIToolsDev`", "QMRITools`"];
	(*run the step with the messages in a file*)
	With[{
			loadContext = load,
			stepName = Context[st] <> SymbolName[st],
			folder = ExpandFileName[dir],
			stepOptions = Replace[opts, (k_ -> v_) :> {Context[k] <> SymbolName[k], v}, {1}],
			messageFile = msg
		},
		LinkWrite[link, Unevaluated[Get[loadContext]; Block[{stream, $Output, $Messages},
			stream = OpenWrite[messageFile];
			$Output = {};
			$Messages = {stream};
			Quiet[ToExpression[stepName][folder, 
				Sequence @@ (ToExpression[First[#]] -> Last[#]& /@ stepOptions)], FrontEndObject::notavail];
			Close[stream]]]]];
	link
]


BidsStepStatus[link_, root_, st_, t0_, msg_] := Block[{pkt, state, pat, logs, lines, msgs},
	(*drain the link, a return packet means the step is finished*)
	state = "Running";
	While[state === "Running" && Quiet@LinkReadyQ[link],
		pkt = Quiet@LinkRead[link];
		state = Which[pkt === $Failed, "Crashed", MatchQ[pkt, _ReturnPacket | _ReturnExpressionPacket], "Done", True, "Running"]
	];
	If[state === "Running" && !Quiet@LinkConnectedQ[link], state = "Crashed"];

	(*the logs of this run, newest last*)
	pat = If[st === "BidsDcmToNii", "DcmToNii_*.log", "*_BIDS"<>StringReplace[st, "MuscleBids" -> ""]<>".log"];
	logs = {#, AbsoluteTime[FileDate[#]]}& /@ Quiet@FileNames[pat, root, 3];
	logs = First /@ SortBy[Select[logs, Last[#] > t0 - 2&], Last];
	lines = If[logs === {}, {}, Quiet@Import[Last[logs], "Lines"]];
	lines = If[ListQ[lines], Select[lines, StringTrim[#] =!= ""&], {}];

	(*messages of the step*)
	msgs = If[FileExistsQ[msg], Quiet@Import[msg, "Lines"], {}];
	msgs = If[ListQ[msgs], msgs, {}];

	<|
		"State" -> If[state === "Done" && msgs =!= {}, "Done with messages", state],
		(*last line without the time stamp, the slash and the dash*)
		"Now" -> If[lines === {}, "", StringTrim@StringReplace[Last[lines], RegularExpression["^\\s*(\\d\\d-\\d\\d-\\d\\d \\d\\d:\\d\\d:\\d\\d)?\\s*/\\s*(- )?"] -> ""]],
		"Touched" -> Length[logs],
		"Tail" -> Take[lines, -Min[30, Length[lines]]],
		"Messages" -> msgs
	|>
]


(* ::Subsubsection::Closed:: *)
(*CheckDataDescription*)


SyntaxInformation[CheckDataDescription] = {"ArgumentsPattern" -> {_, _}};

CheckDataDescription[disIn_Association, met_] := Block[{
		dis, val, vals, keys, duplicate, disOut
	},

	(*all data names should be unique*)
	If[!DuplicateFreeQ[Keys[disIn]],
		(*datasets cannot have duplicate names*)
		Return[Message[CheckDataDescription::key]; $Failed]
		,
		dis = If[met === "MuscleBidsAnalysis",
        	(* If nested like, grab all sub-keys; if flat, use the block directly *)
			If[KeyExistsQ[disIn, "Analysis"], {"Default"->disIn}, disIn], 
			disIn
		];
		keys = Keys[dis];
		vals = Values[dis];
		(*check if there are duplicated datasets, i.e. same type and suffix*)
		duplicate = !DuplicateFreeQ[{#["Type"], #["Suffix"]} & /@ vals];
		debugBids["data Descriptions duplicates: ", duplicate];

		(*Add the data key and if it has duplicates*)
		disOut = MapThread[Join[#1, <|"Key" -> StringStrip[#2], "HasDuplicate" -> duplicate|>] &
			, {vals, keys}];
		If[met === "MuscleBidsAnalysis",
			disOut,
			Flatten[CheckDataDescription[Normal[#], met]& /@ disOut]
		]
	]
]

CheckDataDescription[dis:{_Rule..}, met_] := Block[{
		manLab, ass, key, man, cls, typ, fail
	},

	(*Get the data Description keys*)
	ass = Association[dis];
	key = Keys[ass];
	(*fail output*)
	fail = Association[$Failed->ToString[Normal[ass]]];

	manLab = Switch[met,
		"MuscleBidsConvert", {"Label", "Type"},
		"MuscleBidsProcess", {"Type"},
		"MuscleBidsMerge", {"Type", "Merging"},
		"MuscleBidsSegment", {"Type"},
		"MuscleBidsTractography", {"Type"}
	];

	(*Check if mandatory keys are present*)
	man = ContainsAll[key, manLab];

	If[!man,
		Return[Message[CheckDataDescription::man, manLab]; fail],

		(*Check if type is valid*)
		If[!MemberQ[Keys[bidsTypes], ass["Type"]], Message[CheckDataDescription::type, ass["Type"]]];

		(*Check if class is present*)
		If[KeyExistsQ[ass, "Class"],
			(*if present add check class is valid*)
			If[!MemberQ[bidsClass, ass["Class"]], Return[Message[CheckDataDescription::class,ass["Class"]]; fail]],
			(*add class if not present*)
			ass = Association[ass, "Class"->"Volume"]
		];

		(*check if labels match class*)
		cls = Switch[ass["Class"],
			"Volume", StringQ[ass[["Label"]]],
			"Stacks"|"Volumes"|"Repetitions"|"Acquisitions"|"Mixed", ListQ[ass["Label"] && Length[ass]>1]
		];
		If[!cls, Return[Message[CheckDataDescription::lab, ass["Class"], ass["Label"]]; fail]];

		(*check suffix, in and out folder*)
		If[!KeyExistsQ[ass, "Suffix"], ass = Association[ass, "Suffix"->""]];
		Switch[met,
			"MuscleBidsConvert",
			If[!KeyExistsQ[ass, "InFolder"], ass = Association[ass, "InFolder"->"raw"]];
			ass = Association[ass, "OutFolder" -> BidsType[ass["Type"]]];
			,
			"MuscleBidsProcess"|"MuscleBidsMerge"|"MuscleBidsSegment"|"MuscleBidsTractography",
			If[!KeyExistsQ[ass, "InFolder"], ass = Association[ass, "InFolder"->BidsType[ass["Type"]]]];
			ass = Association[ass, "OutFolder" -> BidsType[ass["Type"]]];			
		];

		(*output the completed data Description*)
		{KeySort@ass}
	]
]


(* ::Subsection::Closed:: *)
(*0. BidsFolderLoop*)


Options[BidsFolderLoop] = {
	(*loop method*)
	Method->"MuscleBidsConvert", 
	(*general options*)
	ProcessSubjects->All, 
	VersionCheck->False, 
	(*method specific options*)	
	DeleteAfterConversion->True, 
	BidsTractographyMethod->"Full",
	BidsOutputImages->"All"
};

SyntaxInformation[BidsFolderLoop] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

BidsFolderLoop[inFol_?StringQ, outFol_?StringQ, opts:OptionsPattern[]] := BidsFolderLoop[inFol, outFol, {""}, opts]

BidsFolderLoop[inFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := BidsFolderLoop[inFol, inFol, datDis, opts]

BidsFolderLoop[inFol_?StringQ, outFol_?StringQ, datDisIn_?AssociationQ, OptionsPattern[]] := Block[{
		met, datType, folders, subs, logFile, ass, nam, filesSl, versCheck, delete, tractMet,
		files, nTyp, pat, inputFolders, customConfig, out, datDis, imOut
	},

	{met, subs, versCheck, delete, tractMet, imOut} = OptionValue[{
		Method, ProcessSubjects, VersionCheck, DeleteAfterConversion, BidsTractographyMethod, BidsOutputImages}];

	debugBids["Enter BidsFolderLoop for method: "<>met];
	debugBids[inFol];

	(*select the subjects and folders to be processed*)
	folders =Switch[met,
		"BidsDcmToNii", ResetLog[]; Select[FileNames[All, inFol], DirectoryQ],
		_, SelectBids[inFol, "ses"]
	];

	subs = FilterBidsSubjects[folders, subs, met];

	debugBids[Column@subs];

	(*loop over the subjects*)
	Table[
		debugBids["input folder:", fol];

		(* -------------- Config and naming --------------*)

		(*get the BidsName*)
		ass = SubNameToBids[fol, met];
		nam = GenerateBidsName[ass];
		out = GenerateBidsFolderName[outFol, ass];

		debugBids["output folder:", out];

		(*check for custom config - merge if config exists in input folder and copy it to output folder*)
		(* Load config: use local/subject config for all but Analysis *)
		{customConfig, datDis} = If[met === "MuscleBidsAnalysis",
			{False, datDisIn}, 
			CheckConfig[fol, out]
		];
		debugBids[{customConfig, datDis}];

		(* Merge global/local and inject Keys/Duplicate logic*)
		datDis = CheckDataDescription[If[met === "MuscleBidsAnalysis", datDis, MergeConfig[datDisIn, datDis]], met];

		(*-------------- Logging --------------*)
		debugBids[{"starting logging: ", met}];
		SetLogFile@Switch[met,
			"BidsDcmToNii", FileNameJoin[{outFol, "DcmToNii_"<>DateName[]<>".log"}],
			"MuscleBidsConvert", FileNameJoin[{fol, nam<>"_BIDSConvert.log"}],
			"MuscleBidsProcess", FileNameJoin[{out, nam<>"_BIDSProcess.log"}],
			"MuscleBidsMerge", FileNameJoin[{out, nam<>"_BIDSMerge.log"}],
			"MuscleBidsSegment", FileNameJoin[{out, nam<>"_BIDSSegment.log"}],
			"MuscleBidsTractography", FileNameJoin[{out, nam<>"_BIDSTractography.log"}],
			"MuscleBidsAnalysis", FileNameJoin[{out, nam<>"_BIDSAnalysis.log"}]
		];
		If[met=!= "BidsDcmToNii", ImportLog[]]; 
		ShowLog[];
		
		(*-----*)AddToLog[{"Starting "<>met<>" for directory: ", fol}, True, 0];
		(*-----*)If[customConfig, AddToLog["********** ----- Using custom config ----- **********", 0]];
		If[met === "BidsDcmToNii",
			(*-----*)AddToLog["Using Chris Rorden's dcm2niix.exe (https://github.com/rordenlab/dcm2niix)", 1];
		];

		(* -------------- The actual processing loops -------------- *)
		debugBids[{"starting method: ", met}];
		Switch[met,
			"BidsDcmToNii", BidsDcmToNiiI[fol, out, datDisIn],
			"MuscleBidsAnalysis", MuscleBidsAnalysisI[{fol, outFol}, #, versCheck, imOut]& /@ datDis
			,
			_,
			(*loop over the data Descriptions*)
			Table[
				(*check if datDis is valid*)
				If[KeyExistsQ[type, $Failed],
					(*-----*)AddToLog[dataToLog@type, 2, True];
					(*-----*)AddToLog["Skipping", 3],
					(*if valid perform conversion*)
					(*-----*)AddToLog[dataToLog@type, 2, True];
					inputFolders = SelectBids[fol, type["InFolder"]];
					(*method specific scripts: loop over all folders in subject/session folder*)
					Table[
						(*loop needs feedback to show where bugs are.*)
						MonitorFunction[DateString[], folIn<>" - "<>type["Key"]<>":"];
						Switch[met,
							"MuscleBidsConvert", MuscleBidsConvertI[folIn, type, delete],
							"MuscleBidsProcess", MuscleBidsProcessI[{folIn, outFol}, type, versCheck],
							"MuscleBidsMerge", MuscleBidsMergeI[{folIn, outFol}, {type, datDis}, versCheck],
							"MuscleBidsSegment", MuscleBidsSegmentI[{folIn, outFol}, {type, datDis}, versCheck],
							"MuscleBidsTractography", MuscleBidsTractographyI[{folIn, outFol}, {type, datDis}, versCheck, tractMet]
						];(*close method switch*)
					, {folIn, inputFolders}];(*Close sub folders loop*)		
				];(*close type check*)
			, {type, datDis}];(*close datatype loop*)
		]; (*close method switch*)

		(*clear log filename*)
		SetLogFile[];

	, {fol, subs}];(*close subject/ses folder loop*)
]


(* ::Subsection:: *)
(*1. BidsDcmToNii*)


(* ::Subsubsection::Closed:: *)
(*BidsDcmToNii*)


Options[BidsDcmToNii]={BidsIncludeSession->True, ProcessSubjects->All}

SyntaxInformation[BidsDcmToNii] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

BidsDcmToNii[folder_?StringQ, opts:OptionsPattern[]] := BidsDcmToNii[folder, GetConfig[folder], opts];


BidsDcmToNii[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir},
	debugBids["starting BidsDcmToNii"];
	dir = Directory[]; SetDirectory[folder];
	BidsDcmToNii[
		ConfigLookup[config, "folders", "dicomData"],(*the input folder of the dcm data*)
		ConfigLookup[config, "folders", "rawData"],(*the output folder for conversion*)
		ConfigLookup[config, "conversion"],(*the conversion settings*)
		opts];
	SetDirectory[dir];
]


BidsDcmToNii[inFol_?StringQ, outFol_?StringQ, settings_, opts:OptionsPattern[]] := BidsFolderLoop[inFol, outFol, settings, Method->"BidsDcmToNii", opts]


(* ::Subsubsection::Closed:: *)
(*BidsDcmToNiiI*)


BidsDcmToNiiI[fol_, outI_, settings_] := Block[{out},
	(*define the out folder*)
	out = FileNameJoin[{outI, "raw"}];
	(*-----*)AddToLog[{"Output folder: ", out}, 1];
	Quiet[CreateDirectory[out]];

	(*perform the conversions only when output folder is empty*)
	If[EmptyDirectoryQ[out],
		(*perform conversion*)			
		(*-----*)AddToLog["Starting the conversion", 1, True];
		DcmToNii[FileNameJoin[{Directory[],#}]&/@{fol,out}, 
			MonitorCalc->False, UseVersion->settings["Version"]];
		(*-----*)AddToLog["Folder was converted", 1],
		(*-----*)AddToLog["Folder was skipped since output folder already exists", 1];
	];
]


(* ::Subsection:: *)
(*2. MuscleBidsConvert*)


(* ::Subsubsection::Closed:: *)
(*MuscleBidsConvert*)


Options[MuscleBidsConvert] = {DeleteAfterConversion->True, ProcessSubjects->All};

SyntaxInformation[MuscleBidsConvert] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

MuscleBidsConvert[folder_?StringQ, opts:OptionsPattern[]] := MuscleBidsConvert[folder, GetConfig[folder], opts];


MuscleBidsConvert[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir},
	debugBids["starting MuscleBidsConvert"];
	dir = Directory[]; SetDirectory[folder];
	MuscleBidsConvert[
		ConfigLookup[config, "folders", "rawData"],(*the input and output folder for the data*)
		ConfigLookup[config, "folders", "rawData"],
		config["datasets"],(*what data for conversion*)
		opts];
	SetDirectory[dir];
]


MuscleBidsConvert[niiFol_?StringQ, outFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := BidsFolderLoop[niiFol, outFol, datDis, Method->"MuscleBidsConvert", opts]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsConvertI*)


MuscleBidsConvertI[folIn_, datType_, del_] := Block[{
		type, fol, parts, files, json, infoExtra, pos, posIn, info, data, vox,
		grad, val, suffix, outFile, echo, nEch, fit, labels, class, types,
		vx, vy, vz, dx, dy, dz, sx, sy, sz, off, diffFile, hasBval, hdr, nSl, len,
		labs, noFiles, sel, imgTags, recTags, consistent
	},

	debugBids["Starting MuscleBidsConvertI"];
	debugBids[folIn];

	(*see if one label or session|repetition*)
	{fol, parts} = PartitionBidsFolderName[folIn];

	(*user defined data description*)
	labels = Flatten[{datType["Label"]}];
	type = datType["Type"];
	class = datType["Class"];
	
	(*for joining dti data*)
	labels = If[class === "Volumes", {labels}, labels];

	(*-----*)AddToLog[{"Converting", ToString[Length[labels]], class}, 2, True];
	(*-----*)AddToLog[StringJoin@@Riffle[labels,", "], 3];

	(*loop over stack names*)
	Table[
		(*import the json belonging to name*)
		(*-----*)AddToLog[{"Converting", nameIn, "as", type,":"}, True, 3];
		files = Flatten@If[class === "Volumes",
			FileNames["*"<>StringReplace[#," "->"_"]<>"*.json", folIn]& /@ nameIn,
			FileNames["*"<>StringReplace[nameIn," "->"_"]<>"*.json", folIn]
		];

		(*add settings to output json*)
		infoExtra = <|
			"ConversionSoftware"->"QMRITools.com", 
			"ConversionSoftwareVersion"->QMRITools`$InstalledVersion,
			If[class==="Volume", "Volume"->nameIn, Nothing],
			If[class==="Volumes", "Volumes"->nameIn, Nothing],
			If[class==="Repetitions", "Repetition"->nameIn, Nothing],
			If[class ==="Acquisitions", "Acquisition"->nameIn, Nothing],
			If[class==="Stacks" || class==="Chunks" || class==="Mixed", "Stack"->nameIn, Nothing],
			If[class==="Stacks" || class==="Chunks" || class==="Mixed", "OverLap"->datType["Overlap"], Nothing]
		|>;

		(*check if there are files to do something with*)
		If[Length@Flatten[files]===0,
			(*no json files found*)
			debugBids["!!!!!!!!! No json files found !!!!!!!!!!"];
			(*-----*)AddToLog[{"No json files found with label ", nameIn , " skipping conversion"}, 4],

			(*if json files found import them*)
			debugBids[files];
			json = ImportJSON/@files;

			(*see which data type is expected*)
			debugBids["Converting ",datType["Type"], " as ", datType["Process", "Method"]];
			Switch[datType["Type"],

				(*-------------------------------------------------*)
				(*-------- DIXON conversion script megre ----------*)
				(*-------------------------------------------------*)
				"megre"|"tse",

				Switch[datType["Process", "Method"],
					"Dixon-P",
					(*philips data with online recon*)
					suffix = datType["Process", "Types"];
					debugBids[suffix];

					If[ListQ[suffix]&&!MatrixQ[suffix],
						(*data is a 4D dataset and the order is in types*)

						pos = posIn = GetJSONPosition[json, {{"ProtocolName", nameIn}}];
						debugBids["Converting Dix data, json position: ", pos];

						If[pos=!={},
							pos = CheckPos[pos];
							
							(*get the json and data*)
							(*-----*)AddToLog[{"Importing dataset with properties: ", nameIn}, 4];
							info = json[[pos]];
							{data, vox} = ImportNii[ConvertExtension[files[[pos]],".nii"], NiiScaling->False];
							(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];

							(*loop over types and export them*)
							Table[(*export to the correct folder*)
								outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
									"suf"->Flatten@{datType["Suffix"], suffix[[i]]}|>];
								(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
								ExportNii[data[[All, i]], vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
								Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];
							,{i, 1, Length[suffix]}];

							(*export used files*)
							DelFiles[files[[posIn]], del];
						];
					];

					If[ListQ[suffix]&&MatrixQ[suffix],
						(*data is multiple files and the order and extension is in types*)
						{suffix, types} = Transpose[suffix];
						
						pos = posIn = GetJSONPosition[json, {{"ProtocolName", nameIn}}];
						debugBids["Converting Dix data, json position: ", pos];
						(*get the json and data*)
						(*-----*)AddToLog[{"Importing ", Length[pos], "dataset with properties: ", nameIn}, 4];

						labs = Last[StringSplit[FileBaseName[#], "_"]] & /@files;
						debugBids["File labels: ", labs];

						(*loop over types and export them*)
						Table[(*export to the correct folder*)
							pos = Flatten[Position[labs, types[[i]]]];

							If[pos=!={},
								pos = CheckPos[pos];
								outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
									"suf"->Flatten@{datType["Suffix"], suffix[[i]]}|>];
								(*-----*)AddToLog[{"Importing:", files[[pos]]}, 4];
								info = json[[pos]];
								{data, vox} = ImportNii[ConvertExtension[files[[pos]],".nii"], NiiScaling->False];
								(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];
								ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
								Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];
								(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
							]
						, {i, 1, Length[suffix]}];

						(*export used files*)
						DelFiles[files[[posIn]], del];
					];

					,
					"Dixon-S",
					(*siemens data with online recon*)

					{suffix, types} = Transpose@datType["Process", "Types"];

					debugBids[types];

					Table[
						pos = posIn = GetJSONPosition[json, {{"SeriesDescription", nameIn<>"_"<>types[[i]]}}];
						If[pos=!={},
							debugBids[{"Dixon-S", nameIn<>"_"<>types[[i]], pos}];
							(*-----*)AddToLog[{"Importing dataset with properties: ", {nameIn, i}}, 4];
							(*get the json and data*)
							pos = CheckPos[pos];
							info = json[[pos]];
							{data, vox} = ImportNii[ConvertExtension[files[[pos]],".nii"], NiiScaling->False];

							(*export to the correct folder*)
							outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
								"suf"->Flatten@{datType["Suffix"], suffix[[i]]}|>];
							debugBids[{outFile, GetClassName[class, nameIn]}];
							(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
							ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
							Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];
							
							(*export used files*)
							DelFiles[files[[posIn]], del];
						]
					, {i, Length@types}];

					,
					"Dixon-B",

					(*non default with data with 4D nii, where data correction is needed*)
					pos = posIn = GetJSONPosition[json, {{"ProtocolName", nameIn}}];
					debugBids["Converting Dix data, json position: ", pos];
					(*get the json and data*)
					(*-----*)AddToLog[{"Importing dataset with properties: ", nameIn}, 4];

					(*get the json and data*)
					pos = CheckPos[pos];
					info = json[[pos]];
					{data, vox} = ImportNii[ConvertExtension[files[[pos]],".nii"], NiiScaling->False];
					(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];

					(*assuming data has source and echos figure out the echos and slices*)
					{nSl, nEch} = Dimensions[data][[1 ;; 2]];
					
					debugBids[Last[info["ImageType"]]];
					(*extract and convert the relevant data*)
					(*check for data with only B0 (SMA AMC)*)
					If[Last[info["ImageType"]] === "FFE",
						(*-----*)AddToLog[{"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"}, 4];
						nEch = (nEch - 1)/4;
						data = Partition[#, nEch] & /@ Partition[Flatten[data, 1][[;; 4 (nSl*nEch)]], nSl nEch];
						vox = vox / {4, 1, 1};
						data[[4]] = data[[4]]/1000.;
						info = KeyDrop[info, "EchoTime"];
						,
						nEch = (nEch - 5)/4;
						data = Partition[#, nEch] & /@ Partition[Flatten[Transpose[data], 1][[;; 4 (nSl*nEch)]], nSl nEch];
						(*mag, real, imag, phase*)
						data = {1000. data[[1]]/2047., 1000. (data[[2]] - 2047.)/2047., 1000. (data[[3]] - 2047.)/2047., Pi (data[[4]] - 2047.)/2047.};
					];
					(*-----*)AddToLog[{"Slices:", nSl, "; Echos:", nEch}, 4];

					echo = datType["Process", "EchoTime"];
					infoExtra = Join[infoExtra, <|
						"EchoNumber" -> Range[nEch], 
						"EchoTime" -> (echo[[1]] + Range[0, nEch - 1] echo[[2]])/1000.,
						"ForthDimension"->"EchoTime",
						"DataClass"->class
					|>];

					Table[
						{i, suffix} = Switch[dixType, "Mixed", {1, ""}, "Phase", {4, "ph"}, "Real", {2, "real"}, "Imaginary", {3, "imag"}];
						(*export to the correct folder*)
						outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
							"suf"->Flatten@{datType["Suffix"], suffix}|>];
						(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
						ExportNii[data[[i]], vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
						Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];

					, {dixType, {"Mixed", "Phase", "Real", "Imaginary"}}];

					(*Delete used files*)
					DelFiles[files[[posIn]], del];
					,
					"Dixon-A",
					(*philips data, auto detect online recon/tse vs raw source data from the ImageType tag*)
					imgTags = ToLowerCase[Last[#["ImageType"]]]& /@ json;
					recTags = Select[{"WATER", "FAT", "IN_PHASE", "OUT_OF_PHASE"}, MemberQ[imgTags, ToLowerCase[#]]&];
					debugBids[recTags];

					If[recTags=!={},
						(*online recon/tse data: one file per dixon type*)
						types = {"WATER", "FAT", "IN_PHASE", "OUT_OF_PHASE"};
						suffix = {"wat", "fat", "inph", "outph"};
						posIn = {};

						(*loop over types and export them*)
						Table[(*get the position of the files needed*)
							pos = GetJSONPosition[json, {{"ProtocolName", nameIn}, {"ImageType", types[[i]]}}];

							If[pos==={},
								(*-----*)AddToLog[{"No json files found with label ", nameIn , " and type ", suffix[[i]], " skipping conversion"}, 4],

								pos = CheckPos[pos];
								posIn = Join[posIn, Flatten[{pos}]];

								(*get the json and data*)
								(*-----*)AddToLog[{"Importing dataset with properties: ", nameIn, suffix[[i]]}, 4];
								info = json[[pos]];
								{data, vox} = ImportNii[ConvertExtension[files[[pos]], ".nii"], NiiScaling -> False];
								(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];

								(*export to the correct folder*)
								outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn],
									"suf"->Flatten@{datType["Suffix"], suffix[[i]]}|>];
								(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
								ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
								Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];
							]
						, {i, 1, Length[suffix]}];

						(*export used files*)
						If[posIn=!={}, DelFiles[files[[posIn]], del]];
						,
						(*raw source data: magnitude, phase, real and imaginary per echo*)
						types = {"MIXED", "PHASE", "REAL", "IMAGINARY"};
						suffix = {"", "ph", "real", "imag"};

						(*loop over types and export them*)
						Table[(*get the position of the files needed*)
							pos = GetJSONPosition[json, {{"ProtocolName", nameIn}, {"ImageType", types[[i]]}}, "EchoNumber"];

							If[pos==={},
								(*-----*)AddToLog[{"No json files found with label ", nameIn , " and type ", suffix[[i]], " skipping conversion"}, 4],

								(*get the json and data*)
								(*-----*)AddToLog[{"Importing", Length[pos], "datasets with properties: ", {nameIn, suffix[[i]]}}, 4];
								info = MergeJSON[json[[pos]]];
								{data, vox} = Transpose[ImportNii[#, NiiScaling -> False]& /@ ConvertExtension[files[[pos]], ".nii"]];
								vox = First@vox;
								data = Transpose[data];

								(*correct data for different types*)
								data = Switch[suffix[[i]],
									"", 1000. data / 2047.,
									"ph", Pi (data - 2047.) / 2047.,
									"real" | "imag", 1000. (data - 2047.) / 2047.
								];
								(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];

								(*make the additional mandatory bids json values*)
								infoExtra = Join[infoExtra, <|"ForthDimension"->"EchoTime", "DataClass"->class|>];

								(*export to the correct folder*)
								outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn],
									"suf"->Flatten@{datType["Suffix"], suffix[[i]]}|>];
								(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
								ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
								Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];

								(*Delete used files*)
								DelFiles[files[[pos]], del];
							]
						, {i, 1, Length[suffix]}]
					];
					,
					"Dixon",
					(*default script with bids standard of each echo in one file*)
					(*get the position of the files needed loop over dixon data types*)
					Table[
						(*get the position of the files needed*)
						pos = GetJSONPosition[json, {{"ProtocolName", nameIn}, {"ImageType", dixType}}, "EchoNumber"];

						(*-----*)AddToLog[{"Importing", Length[pos], "datasets with properties: ", {nameIn, dixType}}, 4];
						debugBids["Converting Dix data of type ", dixType, ", json position: ", pos];

						If[pos==={},
							(*-----*)AddToLog[{"No json files found with label ", nameIn , " and type ", dixType, " skipping conversion"}, 4],
							(*if json files found import them*)
							info = MergeJSON[json[[pos]]];
							{data, vox} = Transpose[ImportNii[#]&/@ConvertExtension[files[[pos]], ".nii"]];

							If[ToExpression[First@StringSplit[info["SoftwareVersions"], "\\" | "."]] >= 12,
								(*philips release 12*)
								Which[
									(*data export of GE echo data (MOTION)*)
									First[info["ImageType"]] === "ORIGINAL",
									data = Partition[Flatten[data, 1], Length[data]];
									,
									(*data export of DIXON echo data (MOTOR/TWITCH)*)
									First[info["ImageType"]] === "DERIVED",
									{e, s} = info /@ {"EchoTrainLength", "PhaseEncodingStepsOutOfPlane"};

									debugBids["!!!!!! 12.3 Derived: data dimensions:", Dimensions@data];
									debugBids["Echo / Slices:", {e, s, e/2, s/2}];
									
									Switch[dixType, 
										"Mixed",
										data = Partition[Flatten[data, 1], e];
										Table[data[[i]] = RotateRight[data[[i]]], {i, 1, s, e/2}];
										info["EchoTime"] = Sort[info["EchoTime"]];
										,
										"Phase",
										phase1 = Partition[Flatten[Flatten[data, {1, e/2}][[Range[1, e, 2]]], 1], e];
										Table[phase1[[i]] = RotateRight[phase1[[i]]], {i, 2, s/2, e/2}];
										phase2 = Partition[Flatten[Flatten[data, {1, e/2}][[Range[2, e, 2]]], 1], e];
										Table[phase2[[i]] = RotateRight[phase2[[i]]], {i, 2, s/2, e/2}];
										data = Join[phase1, phase2];
									]
								];
								, 
								(*philips older release*)
								data = Transpose[data]
							];

							vox = First@vox;
							(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];

							(*correct data for different types*)
							{data, suffix} = Switch[dixType,
								"Mixed", {1000. data / 2047.,""},
								"Phase", {Pi (data - 2047.) / 2047.,"ph"},
								"Real", {1000.(data - 2047.) / 2047.,"real"},
								"Imaginary", {1000. (data - 2047.) / 2047.,"imag"}
							];

							(*make the additional mandatory bids json values*)
							infoExtra = Join[infoExtra, <|
								"ForthDimension"->"EchoTime",
								"DataClass"->class
							|>];

							(*export to the correct folder*)
							outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
								"suf"->Flatten@{datType["Suffix"], suffix}|>];
							(*-----*)AddToLog[{"Exporting to file:", outFile}, 4];
							ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
							Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];

							(*Delete used files*)
							DelFiles[files[[pos]], del];
						]
					(*Close loop over dixon data types*)
					, {dixType, {"Mixed", "Phase", "Real", "Imaginary"}}]
				],

				(*-------------------------------------------*)
				(*---------- DWI conversion script ----------*)
				(*-------------------------------------------*)
				"dwi",

				noFiles = False; 

				Switch[class,
					"Volumes",
					(*if volumes get the position of the files needed*)
					pos = posIn = GetJSONPosition[json, {{"ProtocolName", #}}] & /@ nameIn;
					debugBids["Converting DWI data ", nameIn, ", json position: ", pos];
					
					If[Flatten[pos]==={}, 
						noFiles = True;
						(*-----*)AddToLog[{"No json files found with label ", nameIn , " skipping conversion"}, 4],

						(*get the json and data*)
						(*-----*)AddToLog[{"Importing datasets with properties: ", nameIn}, 4];
						posIn = Flatten[posIn];
						pos = CheckPos/@pos;
						info = json[[First@pos]];
						diffFile = ConvertExtension[files[[#]],".nii"]&/@pos;
						hdr = Last@ImportNii[First@diffFile, NiiMethod -> "header"];

						hasBval = (FileExistsQ[ConvertExtension[#,".bval"]]&&FileExistsQ[ConvertExtension[#,".bvec"]])&/@diffFile;
						hasBval = AllTrue[hasBval, # === True &];

						If[hasBval, 
							{data, grad, val, vox} = Transpose[ImportNiiDiff[#, FlipBvec->False]&/@ diffFile];
							{data, grad, val, vox} = ConcatenateDiffusionData[data, grad, val, vox];
							,
							(*-----*)AddToLog[{"!!!!!!!!!!!!!!! WARNING NO BVAL OR BVEC FILE !!!!!!!!!!!!!!!!!!!"}, 4];
						]
					];

					,_ ,
					(*get the position of the files needed*)
					sel = If[First[json]["Manufacturer"]==="Siemens",
						{{"SeriesDescription", nameIn}},
						{{"ProtocolName", nameIn}}
					];
					pos = posIn = GetJSONPosition[json, sel];
					debugBids["Converting DWI data ", nameIn, ", json position: ", pos];

					If[pos==={}, 
						noFiles = True;
						(*-----*)AddToLog[{"No json files found with label ", nameIn , " skipping conversion"}, 4],

						(*get the json and data*)
						(*-----*)AddToLog[{"Importing dataset with properties: ", nameIn}, 4];
						pos = CheckPos[pos];
						info = json[[pos]];
						diffFile = ConvertExtension[files[[pos]],".nii"];
						hdr = Last@ImportNii[diffFile, NiiMethod -> "header"];
						
						hasBval = FileExistsQ[ConvertExtension[diffFile,".bval"]] && FileExistsQ[ConvertExtension[diffFile,".bvec"]];

						If[hasBval, 
							{data, grad, val, vox} = ImportNiiDiff[diffFile, FlipBvec->False],
							(*-----*)AddToLog[{"!!!!!!!!!!!!!!! WARNING NO BVAL OR BVEC FILE !!!!!!!!!!!!!!!!!!!"}, 4];
						];
						(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];
					]
				];

				If[!noFiles,
					(*get the offset to add to header*)
					{vz, vy, vx} = vox;
					{sz, sy, sx} = vox Dimensions@data[[All, 1]];
					{dz, dy, dx} = {"qOffsetZ", "qOffsetY", "qOffsetX"} /. hdr;
					off = -{(sz - vz)/2, (dy + sy) - vy, dx};

					(*make the additional mandatory bids json values*)
					infoExtra = Join[infoExtra, <|
						"Offset"->off,
						"ForthDimension"->"Diffusion",
						"DataClass"->class
					|>];

					If[class==="Volumes",
						files = Flatten[files];
						pos = Flatten[pos];
						nameIn = First@nameIn;
					];

					(*export to the correct folder*)
					debugBids[{parts, type, GetClassName[class, nameIn], datType["Suffix"]}];
					outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
						"suf"->{datType["Suffix"]}|>];
					(*-----*)AddToLog[{"Exporting to file:", outFile}, 5];
					If[hasBval,
						ExportBval[val, ConvertExtension[outFile, ".bval"]];
						ExportBvec[grad, ConvertExtension[outFile, ".bvec"]];
					];
					ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
					Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];

					(*Delete used files*)
					DelFiles[files[[posIn]], del];
				],

				(*-------------------------------------------*)
				(*----------- T2 conversion script ----------*)
				(*-------------------------------------------*)
				"mese",

				(*skip when no json has exactly this protocol name*)
				If[GetJSONPosition[json, {{"ProtocolName", nameIn}}] === {},
				(*-----*)AddToLog[{"No json with protocol name ", nameIn, " found, skipping conversion"}, 4],

				(*newer dcm2niix writes one file per echo in the same layout for every release, older ones need the old handling*)
				noFiles = False;
				consistent = BidsDcmConsistentQ[json[[GetJSONPosition[json, {{"ProtocolName", nameIn}}]]]];
				debugBids["dcm2niix with consistent output: ", consistent];

				Which[
(*TODO
I have to properly implement the logic that figures out if its siemens of Philips and if its individual files of 4D and if the correct 
echo time is in the json. and then spit out warnings and skip when things are going wrong
*)
					(*if echo time exists assume 4D nii without correct echo times, only older dcm2niix writes 4D nii*)
					KeyExistsQ[Lookup[datType, "Process", <||>], "EchoTime"] && !consistent,

					pos = posIn = GetJSONPosition[json, {{"ProtocolName", nameIn}}];
					debugBids["Converting MESE data, json position: ", pos];

					(*get the json and data*)
					(*-----*)AddToLog[{"Importing dataset with properties: ", nameIn}, 4];
					pos = CheckPos[pos];
					(*echo times - get the echo time from the data json if present*)
					info = json[[pos]];

					If[KeyExistsQ[info, "EchoTime"] && KeyExistsQ[info, "EchoTrainLength"] || info["Manufacturer"]==="Siemens",
						{data, vox} = ImportNii[ConvertExtension[files[[pos]],".nii"]];
						,
						{data, fit, vox} = ImportNiiT2[ConvertExtension[files[[pos]],".nii"]];
					];
					(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];
					(*echo = Lookup[info, "EchoTime", datType["Process", "EchoTime"]/1000.];*)
					echo = datType["Process", "EchoTime"]/1000.;
					nEch = Lookup[info, "EchoTrainLength", Length[data[[1]]]];
					info = KeyDrop[info, "EchoTime"];
					echo = <|
						"EchoNumber" -> Range[nEch], 
						"EchoTime" -> (echo + Range[0, nEch - 1] echo)
					|>;
					debugBids[echo];

					,
					(*Else its a list of single echo nii files*)
					True,

					If[consistent && KeyExistsQ[Lookup[datType, "Process", <||>], "EchoTime"],
						(*-----*)AddToLog[{"The EchoTime of the config is ignored, the echo times are in the json files"}, 4]
					];

					(*get the position of the files needed*)
					pos = posIn = GetJSONPosition[json, {{"ProtocolName", nameIn}}, "EchoNumber"];
					debugBids["Converting MESE data ", nameIn, ", json position: ", pos];
					info = MergeJSON[json[[pos]]];

					(*select only echos*)
					len = info["AcquisitionNumber"];
					len = If[ListQ[len], Max[len], Lookup[info, "EchoTrainLength", Length[pos]]];
					pos = pos[[;; Min[len, Length[pos]]]];
					info = MergeJSON[json[[pos]]];

					(*-----*)AddToLog[{"Importing ", Length[pos], "dataset with properties: ", nameIn}, 4];

					(*get the json and data*)
					{data, vox} = Transpose[ImportNii /@ ConvertExtension[files[[pos]],".nii"]];

					(*files of the newer dcm2niix have to be complete and in order: all echos, equal sizes and increasing echo times*)
					noFiles = consistent && !And[
						Length[pos] === Lookup[info, "EchoTrainLength", Length[pos]],
						Length[DeleteDuplicates[Dimensions /@ data]] === 1,
						OrderedQ[Flatten[{Lookup[info, "EchoTime", {}]}], Less]
					];

					If[noFiles,
						(*-----*)AddToLog[{"!!!!!!!!!!!!!!! SKIPPING - INCOMPLETE OR UNORDERED ECHO FILES, convert the dicom data again !!!!!!!!!!!!!!!!!!!"}, 4];
						debugBids["!!!!!!!!!!!!! Echo files: ", Length[pos], Dimensions /@ data, info["EchoTime"]];
						,
						(*older dcm2niix writes release 12 data in another order, newer always one file per echo*)
						data = If[!consistent && ToExpression[First@StringSplit[info["SoftwareVersions"], "\\" | "."]] >= 12,
							AssociateTo[info, "EchoTime" -> Sort[info["EchoTime"]]];
							Partition[Flatten[data, 1], Length@data],(*philips release 12*)
							Transpose[data](*philips older release*)];
						vox = First@vox;
						(*-----*)AddToLog[{"Dimensions:", Dimensions@data, "; Voxel size:", vox}, 4];
					];

					echo = <||>;
				];

				If[!noFiles,

				(*make the additional mandatory bids json values*)
				infoExtra = Join[infoExtra, echo, <|
					"ForthDimension"->"EchoTime",
					"DataClass"->class
				|>];

				(*export to the correct folder*)
				outFile = GenerateBidsFileName[fol, <|parts, "type"->type, GetClassName[class, nameIn], 
					"suf"->{datType["Suffix"]}|>];
				(*-----*)AddToLog[{"Exporting to file:", outFile}, 5];
				ExportNii[data, vox, ConvertExtension[outFile, ".nii"], CompressNii -> compress];
				Export[ConvertExtension[outFile, ".json"], MergeJSON[{info, infoExtra}]];

				(*Delete used files*)
				DelFiles[files[[posIn]], del];
				];
				]

				,

				(*-------------------------------------------*)
				(*-------- Other processing script ----------*)
				(*-------------------------------------------*)
				_,
				(*-----*)AddToLog["Unknown datatype for conversion", 4];

			(*Close Type switch*)
			];
		(*close file check*)	
		];
	(*close loop over stack names*)
	,{nameIn, labels}];

	(*compress the nii files if compression during ExportNii -> False*)
	If[!compress, CompressNiiFiles[fol]];
]


CheckPos[pos_] := Block[{posOut},
	posOut = Last[pos];
	If[Length[pos] > 1,
		(*-----*)AddToLog[{"!!!!!!!!!!!!!!! More than one file found using last !!!!!!!!!!!!!!!!!!!"}, 4];
		debugBids["!!!!!!!!!!!!!!! More than one file found using last: ", posOut];
	];
	posOut
]


DelFiles[files_, del_]:=Quiet@If[del, 
	DeleteFile[ConvertExtension[files,#]]&/@{".nii", ".nii.gz", ".json", ".bval", ".bvec"}
];


(* ::Subsection:: *)
(*3. MuscleBidsProcess*)


(* ::Subsubsection::Closed:: *)
(*MuscleBidsProcess*)


Options[MuscleBidsProcess] = {ProcessSubjects->All, VersionCheck->False};

SyntaxInformation[MuscleBidsProcess] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

MuscleBidsProcess[folder_?StringQ, opts:OptionsPattern[]] := MuscleBidsProcess[folder, GetConfig[folder], opts];


MuscleBidsProcess[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir},
	debugBids["starting MuscleBidsProcess"];
	dir = Directory[]; SetDirectory[folder];
	MuscleBidsProcess[
		ConfigLookup[config, "folders", "rawData"],(*the input folder for the data*)
		ConfigLookup[config, "folders", "derivedData"],(*the output folder for processing*)
		config["datasets"],(*what data for processing*)
		opts];
	SetDirectory[dir];	
]


MuscleBidsProcess[niiFol_?StringQ, outFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := BidsFolderLoop[niiFol, outFol, datDis, Method->"MuscleBidsProcess", opts]


(* ::Subsubsection::Closed:: *)
(*OptionsToJSON*)


(*option rules used in a processing call as json entries, so the exported json shows what was run*)
OptionsToJSON[opts_] := JSONSafe@KeyMap[SymbolName, Association[opts]]


(*values json cannot hold (Automatic, Missing, other symbols) become strings*)
JSONSafe[a_Association] := JSONSafe /@ a

JSONSafe[l_List] := JSONSafe /@ l

JSONSafe[x : (_String | _Integer | _Real | True | False | Null)] := x

JSONSafe[x_] := ToString[x]


(*the json that belongs to a data file, empty if there is none*)
SourceJSON[file_] := If[FileExistsQ[#], ImportJSON[#], <||>]&@ConvertExtension[file, ".json"]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsProcessI*)


MuscleBidsProcessI[{folIn_, folOut_}, datType_, verCheck_] := Block[{
		con, fol, parts, type, files, sets, diffFile, nfile, keys, dixFiles, jsonFile, nFiles, phbpt, dbond,
		outFile, json, echos, mag, ph, real, imag, dixVox, magM, B0mask, ph0i, pos, e1, e2, hz, b0i,
		t2stari, watfr, fatfr, wat, fat, inph, outph, b0, t2star, r2star, phi, itt, res, outTypes, preProc, 
		niiFileP, jsonFileP, resi, data, grad, val, diffVox, mask, den, sig, snr, snr0, reg, valU, mean, 
		fiti, s0i, fri, adci, pD, tens, s0, out, l1, l2, l3, md, fa, rd, t2vox, t2w, t2f, b1, n, 
		angle, ex, ref, thk, phii, phbpi, phbp, ta, filt, field, settingPre, settingPro, regF, coil, off, 
		int, dint, suffix, types, flip, per, bmat, magph, split, ivim, shift, t2, gradField,
		meanV, coor, valV, fasc, fascm, sel, norm, rdim, raw, rmet, rician, sigI, bv, init, exist,
		denOpts, dixOpts, epgOpts, phaseEchoes, procInfo, precRatio, precResult
	},

	debugBids["Starting MuscleBidsProcessI"];
	debugBids[{folIn, folOut}];

	(*get the context for exporting*)
	con = Context[con];
	outTypes = {};

	(*get the information needed for processing, e.g. session|repetition*)
	{fol, parts} = PartitionBidsFolderName[folIn];
	type = datType["Type"];
	keys = {"EchoTime", "ForthDimension", "DataClass", "Stack",
		"OverLap", "SliceThickness", "SpacingBetweenSlices"};

	(*config settings and version that are added to every output json*)
	procInfo = JSONSafe@<|"ProcessSettings" -> datType["Process"], "ProcessVersion" -> QMRITools`$InstalledVersion|>;

	(*see what needs to be processed*)
	files = BidsJsonFiles[folIn, datType["Label"]];
	If[files === {},
		(*-----*)AddToLog[{"!!!!!!!!!!!!!!! SKIPPING - NO FILES FOUND !!!!!!!!!!!!!!!!!!!"}, 2];
		debugBids["No files found"];
		Return[];	
	];

	sets = If[type==="megre"||type==="tse",
		DeleteDuplicates[(ta = #;AssociateTo[ta, "suf"->{First[ta["suf"]]}])&/@PartitionBidsName[FileBaseName/@files]],
		DeleteDuplicates[PartitionBidsName[FileBaseName/@files]]];
	(*-----*)AddToLog[{"Found", ToString[Length[sets]], datType["Class"], "that will be processed:"}, 2];

	debugBids[Column@files];
	debugBids[Column@sets];

	(*loop over sets*)
	Table[
		(*-----*)AddToLog[dataToLog@set, 2];
		debugBids[{set, type, datType["Process", "Method"]}];
		(*see which data type is expected*)
		Switch[type,

			"megre" | "tse",

			(*-------------------------------------------*)
			(*-------- Dixon processing scripts ---------*)
			(*-------------------------------------------*)

			Switch[datType["Process", "Method"],
				"Dixon-S" | "Dixon-P" | "Dixon" | "Dixon-B" | "Dixon-A",

				outFile = GenerateBidsFileName[folOut, set];
				debugBids[outFile];

				(*check if files are already done*)
				If[CheckFile[outFile, "done", verCheck],
					(*if check file has label done and version is recent skip*)
					(*-----*)AddToLog["Processing already done for: ", True, 3];
					(*-----*)AddToLog[outFile, 4],

					(*-----*)AddToLog["Starting processing for data:", 3, True];
					outTypes = {};
					echos = {};
					settingPro = <||>;

					(*----------------------------------------------*)
					(*-------- TSE/GE online recon scripts ---------*)
					(*----------------------------------------------*)

					(*branch 1: data with online reconstruction, default to the standard dixon types if not configured*)
					(*results get overwritten by branch 2 if dixon processing is feasible.*)
					If[MemberQ[{"Dixon-S", "Dixon-P", "Dixon-A"}, datType["Process", "Method"]],
						suffix = Lookup[datType["Process"], "Types", {"wat", "fat", "inph", "outph"}];
						If[MatrixQ[suffix], {suffix, types} = Transpose@suffix];

						dixFiles = GenerateBidsFileName[fol, <|set, "suf"->{datType["Suffix"], #}|>]&/@suffix;
						nFiles = ConvertExtension[dixFiles, ".nii"];
						(*-----*)AddToLog[First@dixFiles, 4];

						exist = NiiFileExistQ /@ nFiles;
						If[!Or@@exist,
							(*-----*)AddToLog[{"Could not find any of the ", First@dixFiles}, 4],
							If[!AllTrue[exist, TrueQ], (*-----*)AddToLog[{"Missing types for ", set, ": ", Pick[suffix, exist, False]}, 4]];
							{suffix, nFiles} = {Pick[suffix, exist], Pick[nFiles, exist]};

							(*-----*)AddToLog["Importing the data", 4];
							{data, dixVox} = Transpose[ImportNii/@nFiles];
							dixVox = First@dixVox;
							json = ImportJSON[ConvertExtension[First[nFiles], ".json"]];

							pos = Flatten[Position[suffix, #] & /@ {"inph", "outph", "fat", "wat"}];
							debugBids[pos];
							mask = If[pos==={}, 1,
								Mask[NormalizeMeanData[Transpose[NormalizeData /@ data[[pos]]]], 15,
									MaskSmoothing -> True, MaskComponents -> 2, MaskClosing -> 2, MaskDilation -> 1]
							];
							debugBids[{Dimensions@mask, Dimensions@data}];

							If[MemberQ[suffix, "wat"], wat = mask data[[Position[suffix, "wat"][[1,1]]]]];
							If[MemberQ[suffix, "fat"], fat = mask data[[Position[suffix, "fat"][[1,1]]]]];
							If[MemberQ[suffix, "inph"], inph = mask data[[Position[suffix, "inph"][[1,1]]]]];
							If[MemberQ[suffix, "outph"], outph = mask data[[Position[suffix, "outph"][[1,1]]]]];

							If[MemberQ[suffix, "t2star"],
								t2star = mask data[[Position[suffix, "t2star"][[1,1]]]] / 10000.;
								r2star = DivideNoZero[1, t2star];
								AppendTo[suffix, "r2star"]];
							If[MemberQ[suffix, "fatfr"],
								fatfr = mask data[[Position[suffix, "fatfr"][[1,1]]]] / 1000.;
								If[!MemberQ[suffix, "watfr"], watfr = mask - fatfr; AppendTo[suffix, "watfr"]];
							];

							If[!MemberQ[suffix, "fatfr"] && MemberQ[suffix, "wat"] && MemberQ[suffix, "fat"],
								{watfr, fatfr} = MaskData[DixonToPercent[wat, fat], mask];
								AppendTo[suffix, "fatfr"];
								If[!MemberQ[suffix, "watfr"], AppendTo[suffix, "watfr"]];
							];

							outTypes = suffix;
							settingPro = <|"DixonSource" -> "ScannerRecon"|>;
						];
					(*close branch 1*)
					];

					(*-------------------------------------------*)
					(*-------- GE offline recon scripts ---------*)
					(*-------------------------------------------*)

					(*branch 2: philips data with offline reconstruction not tested for Siemens*)
					If[MemberQ[{"Dixon", "Dixon-B", "Dixon-A"}, datType["Process", "Method"]],
						(*input file names*)
						dixFiles = GenerateBidsFileName[fol, <|set, "suf"->{datType["Suffix"], #}|>]& /@ {"real", "imag"};
						jsonFile = ConvertExtension[First@dixFiles, ".json"];
						nFiles = ConvertExtension[dixFiles, ".nii"];
						magph = False;

						(*if not real imag check for mag phase*)
						If[!FileExistsQ[jsonFile],
							dixFiles = GenerateBidsFileName[fol, <|set, "suf"->{datType["Suffix"], #}|>]& /@ {"", "ph"};
							jsonFile = ConvertExtension[First@dixFiles, ".json"];
							nFiles = ConvertExtension[dixFiles, ".nii"];
							magph = True;
						];

						(*-----*)AddToLog[First@dixFiles, 4];

						(*Check if needed json Exist*)
						If[!FileExistsQ[jsonFile],
							(*-----*)AddToLog["Could not find the needed JSON file", 4];
							,
							(*Check if needed nii Exist*)
							If[!AllTrue[nFiles, NiiFileExistQ],
								(*-----*)AddToLog[{"Could not find all the ", First@dixFiles}, 4];
								,
								(*-----*)AddToLog["Importing the data", 4];

								(*import the data*)
								json = ImportJSON[jsonFile];
								{echos, field} = json /@ {"EchoTime", "MagneticFieldStrength"};

								If[magph,
									{{mag, ph}, dixVox} = Transpose[ImportNii/@nFiles];
									real = mag Cos[ph];
									imag = mag Sin[ph];
									,
									{{real, imag}, dixVox} = Transpose[ImportNii/@nFiles];
									{mag, ph} = Through[{Abs, Arg}[real + I imag]];
								];
								dixVox = First@dixVox;
								debugBids[{Length[echos], Dimensions[real]}];

								(*check if number of echos matches number of volumes*)
								If[Length[echos] =!= Dimensions[real][[2]],
									(*-----*)AddToLog[{"!!!!! Dimensions echos not match !!!!! echos:", Length[echos], "volumes:", Dimensions[real][[2]], "- skipping Dixon fit"}, 4],

									(*Apply background mask*)
									magM = NormalizeMeanData@mag;
									B0mask = Mask[magM, 15, MaskSmoothing->True, MaskComponents->2, MaskClosing->2, MaskDilation->1];
									{real, imag} = MaskData[#, B0mask] &/@ {real, imag};

									(*-----*)AddToLog["Starting denoising and SNR calculation", 4];
									denOpts = {PCAKernel -> 5, Method -> "Patch", PCAComplex -> True};
									{{real, imag}, sig} = PCADeNoise[{real, imag}, Sequence @@ denOpts];
									{mag, ph} = Through[{Abs, Arg}[real + I imag]];
									snr = SNRCalc[Mean@Transpose@mag, sig];

									Switch[datType["Process", "Method"],
										(*Dixon processing scrip for multi echo gradient echo complex data as used in motion study*)
										"Dixon" | "Dixon-A",

										(*see if there are dixon flips*)
										(*{{mag, ph, real, imag}, pos} = FixDixonFlips[{mag, ph, real, imag}];
										(*-----*)If[pos=!={}, AddToLog[{"Found complex flips in volumes: ", pos}, 4]];*)

										(*calculated field maps*)
										(*-----*)AddToLog[{"Starting field map calculation"}, 4];
										{{b0i, t2stari, phii, phbpi}, {e1, e2, n}} = DixonPhase[{real, imag}, echos];
										phaseEchoes = {e1, e2};
										(*-----*)AddToLog[{"used echo ", ToString[e1], "(", 1000 echos[[e1]],"ms ) and", ToString[e2], "(", 1000 echos[[e2]], "ms )"}, 5];

										(*perform the IDEAL dixon fit, with double bond fat model if enough echos*)
										(*-----*)AddToLog["Starting Dixon reconstruction", 4];
										dixOpts = {
											DixonPrecessions -> 1, DixonPhases -> {True, True, True, True, True},
											DixonFixT2 -> False, DixonFieldStrength -> field,
											DixonAmplitudes -> If[Length[echos] > 6, "CallDB", "Fixed"], DixonTolerance -> 1
										};

										If[Length[echos] > 6,
											{{watfr, fatfr}, {wat, fat, dbond}, {inph, outph},
												{{b0, phbp, phi, phbpt}, {t2star, r2star}}, itt, res} = DixonReconstruct[
													{real, imag}, echos, {b0i, t2stari, phii, phbpi}, Sequence @@ dixOpts];
											outTypes = {"dbond", "phbp", "phi", "phbpt", "phii", "phbpi"};
											,
											{{watfr, fatfr}, {wat, fat}, {inph, outph},
												{{b0, phbp, phi, phbpt}, {t2star, r2star}}, itt, res} = DixonReconstruct[
													{real, imag}, echos, {b0i, t2stari, phii, phbpi}, Sequence @@ dixOpts];
											outTypes = {"phbp", "phi", "phbpt", "phii", "phbpi"};
										];

										,
										(*Dixon processing scrip for multi echo gradient echo complex data as used in Bochum cohort*)
										"Dixon-B",

										(*uwrap and convert B0 to hz*)
										(*-----*)AddToLog[{"Starting field map calculation"}, 4];
										b0i = UnwrapSplit[ph[[All, -1]] - ph[[All, 1]], mag, UnwrapDimension -> "3D", MonitorUnwrap -> False];
										b0i = b0i/(2 Pi (Last[echos] - First[echos]));
										phaseEchoes = {1, Length[echos]};
										(*calculate the t2 star from the two in phase images*)
										t2stari = T2Fit[mag, echos][[2]];

										debugBids[Dimensions/@{real,imag, b0i, t2stari}];
										debugBids[echos];

										(*perform the IDEAL dixon fit*)
										(*-----*)AddToLog["Starting Dixon reconstruction", 4];
										dixOpts = {DixonClipFraction -> True, DixonFieldStrength -> field};
										{{watfr, fatfr}, {wat, fat}, {inph, outph},
											{{b0}, {t2star, r2star}}, itt, res} = DixonReconstruct[
												{real, imag}, echos, {b0i, t2stari}, Sequence @@ dixOpts];
									];

									(*check the precession sense of the fat model on the fitted data*)
									{precRatio, precResult} = DixonPrecessionCheck[{real, imag}, echos, fatfr,
										If[datType["Process", "Method"] === "Dixon-B", 0., phbp],
										DixonFieldStrength -> field, DixonPrecessions -> Lookup[dixOpts, DixonPrecessions, 1]];
									(*-----*)AddToLog[{"Precession check:", precResult, "- residual ratio", Round[precRatio, 0.01]}, 5];

									(*settings used for the fit, exported with the data*)
									settingPro = <|
										"DixonSource" -> "QMRITools", "Denoise" -> OptionsToJSON[denOpts],
										"DixonPhaseEchos" -> phaseEchoes, "DixonSettings" -> OptionsToJSON[dixOpts],
										"DixonPrecessionCheck" -> <|"Ratio" -> Round[precRatio, 0.001], "Result" -> precResult|>
									|>;

									(*correct water and fat*)
									{wat, fat} = Abs[{wat, fat}];

									(*add outtypes present in both branches of dixon fitting*)
									outTypes = Join[{"real", "imag", "mag", "ph", "b0i", "t2stari", 
										"b0", "t2star", "r2star", "inph", "outph", "wat", "fat", 
										"watfr", "fatfr", "itt", "res", "snr", "sig"}, outTypes];
								] (*close echo check*)
							] (*close nii check*)
						] (*close json check*)
					];(*close branch 2*)

					(*export the calculated data, check file and compress once, for whichever branch(es) produced data*)
					If[outTypes == {},
						(*no data to process found*)
						(*-----*)AddToLog[{"No data found for either branch, skipping: ", set}, 4]
						,
						(*data was found and processing done export what is needed *)
						(*-----*)AddToLog["Exporting the calculated data to:", 4];
						(*-----*)AddToLog[outFile,5];
						(
							ExportNii[ToExpression[con<>#], dixVox, outFile<>"_"<>#<>".nii", CompressNii -> compress];
							Export[ConvertExtension[outFile <> "_"<>#, ".json"], MergeJSON[{json, procInfo, settingPro}]]
						) & /@ outTypes;

						MakeCheckFile[outFile, Sort@Join[
							{"Check"->"done", "Outputs" -> outTypes, "SetProperties"->set},
							If[echos=!={}, {"EchoTimes"->echos}, {}]
						]];

						If[!compress, CompressNiiFiles[DirectoryName[outFile]]];
					];

					(*-----*)AddToLog["Finished processing", 3, True];
				];(*close dixon processing*)
			,

				_,
				(*-------------------------------------------*)
				(*-------------- Unknown megre --------------*)
				(*-------------------------------------------*)
				(*-----*)AddToLog[{"Unknown processing ", datType["Process"], "for datatype", type}, True, 3];
			],

			"dwi",
			(*-------------------------------------------*)
			(*--------- dwi processing script -----------*)
			(*-------------------------------------------*)

			(*input file names*)
			diffFile = GenerateBidsFileName[fol, set];
			jsonFile = ConvertExtension[diffFile, ".json"];
			nfile = ConvertExtension[diffFile, ".nii"];

			(*output file names*)
			outFile = GenerateBidsFileName[folOut, set];
			debugBids[outFile];

			(*default settings*)
			settingPre = settingPro = <||>;

			(*-------------------------------------------*)
			(*------- dwi pre -processing script --------*)
			(*-------------------------------------------*)

			(*check if pre-processing is already done*)
			preProc = False;
			If[CheckFile[outFile<>"_prep", "done", verCheck],
				(*if check file has label done and version is recent skip*)
				(*-----*)AddToLog["Pre-processing already done for: ", True, 3];
				(*-----*)AddToLog[outFile, 4],
				(*-----*)AddToLog["Starting pre-processing for data:", 3, True];
				(*-----*)AddToLog[diffFile, 4];

				If[!FileExistsQ[jsonFile],
					(*-----*)AddToLog["Could not find the needed JSON file",4],
					(*Check if needed nii Exist*)
					If[!(NiiFileExistQ[nfile] && FileExistsQ[ConvertExtension[nfile,".bval"]] && 
						FileExistsQ[ConvertExtension[nfile,".bvec"]]),
						(*-----*)AddToLog[{"Skipping, could not find .nii, .bval and .bvec"}, 4],
						(*-----*)AddToLog["Importing the data", 4];

						(*import the data*)
						json = ImportJSON[jsonFile];
						{raw, grad, val, diffVox} = ImportNiiDiff[nfile, FlipBvec->False];
						{raw, grad, val} = SortDiffusionData[NormalizeData[raw], grad, val];
						debugBids[{"Raw data dimensions", Dimensions[raw]}];

						(*gradient flip correction*)
						{flip, per} = ConfigLookup[datType, "Process", "FlipPermute"];
						settingPre = MergeConfig[settingPre, <| "FlipPermute" -> {flip, per} |>];
						(*-----*)AddToLog[{"Gradient flips used: ", flip, per}, 4];
						grad = FlipGradientOrientation[grad, flip, per];

						(*Denoise and SNR*)
						(*-----*)AddToLog["Starting dwi denoising", 4];
						mask = Mask[NormalizeMeanData[raw], ConfigLookup[datType, "Process", "Masking"], 
							MaskSmoothing->True, MaskComponents->2, MaskDilation->1];
						rician = ConfigLookup[datType, "Process", "RicianCorrection"];
						settingPre = MergeConfig[settingPre, <| "RicianCorrection" -> rician |>];
						{den, sig} = PCADeNoise[raw, mask, PCAOutput->False, PCATolerance->0, PCAKernel->5,
							PCARicianCorrect -> rician];
						snr = SNRCalc[den, sig];
						snr0 = Mean@Transpose@First@SelectBvalueData[{snr, val}, {0, Max[{2, Min[val]}]}];
						debugBids[{"de-noised data dimensions", Dimensions[den]}];

						(*register data - 2D or 3D - each side separate*)
						(*-----*)AddToLog["Starting dwi motion and eddy correction", 4];
						rdim = ConfigLookup[datType, "Process", "RegistrationDimension"];
						rmet = ConfigLookup[datType, "Process", "RegistrationMethod"];
						split = ConfigLookup[datType, "Process", "SplitRegistration"];

						settingPre = MergeConfig[settingPre, <|
							"RegistrationDimension" -> rdim,
							"RegistrationMethod" -> rmet,
							"SplitRegistration" -> split
						|>];
						split = ConfigLookup[datType, "Process", "SplitRegistration"];
						regF = If[rdim==="2D", split = False; RegisterCardiacData, 
							If[split, RegisterDiffusionDataSplit, RegisterDiffusionData]];
						

						debugBids[{"Registration function", regF}];
						reg = regF[{den, mask, diffVox}, Iterations->300, NumberSamples->5000, 
							PrintTempDirectory->False, MethodReg -> rmet, InterpolationOrderReg ->1, 
							BsplineSpacing -> {30, 30, 30}, BsplineDirections -> {0, 1, 1}
						];

						(*anisotropic filtering*)
						(*-----*)AddToLog["Starting anisotropic data smoothing", 4];
						filt = Unitize[reg] P2SDenoise[reg, mask];

						(*export all the calculated data*)
						(*-----*)AddToLog["Exporting the calculated data to:", 4];
						(*-----*)AddToLog[outFile, 5];
						outTypes = {"raw", "den", "reg", "filt", "sig", "snr0", "snr"};

						(
							ExportNii[ToExpression[con<>#], diffVox, outFile<>"_"<>#<>".nii", CompressNii -> compress];
							If[MemberQ[{"reg", "filt"}, #], ExportBvalvec[{val, grad}, outFile <> "_"<>#]];
							Export[ConvertExtension[outFile<>"_"<>#, ".json"], MergeJSON[{json, procInfo, settingPre}]];
						) & /@ outTypes;

						(*export the check file*)
						MakeCheckFile[outFile<>"_prep", Sort@Join[
							{"Check"->"done", "Bvalue" -> val, "Gradient" -> grad, 
							"Outputs" -> outTypes, "SetProperties"->set},
							Normal@KeyTake[json, keys]
						]];

						(*compress the nii files if compression during ExportNii -> False*)
						If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

						(*-----*)AddToLog["Finished pre-processing", 3, True];

						(*Set pre processing true, overrules check file for processing*)
						preProc = True;
					]
				]
			(*close preprocessing*)
			];

			Switch[datType["Process", "Method"],

				"DTI",
				(*-------------------------------------------*)
				(*---------- dwi processing script ----------*)
				(*-------------------------------------------*)

				(*input file for processing*)
				niiFileP = ConvertExtension[GenerateBidsFileName[folOut, <|set, "suf"->{datType["Suffix"], "filt"}|>],".nii"];
				jsonFileP = ConvertExtension[niiFileP, ".json"];

				(*check if processing is already done, redo is prep is done*)					
				If[If[!preProc, CheckFile[outFile, "done", verCheck], False],
					(*if check file has label done and version is recent skip*)
					(*-----*)AddToLog["Processing already done for: ", True, 3];
					(*-----*)AddToLog[outFile, 4],
					(*-----*)AddToLog["Starting processing for data:", 3, True];
					(*-----*)AddToLog[niiFileP, 4];				

					If[!FileExistsQ[jsonFileP],
						(*-----*)AddToLog["Could not find the needed JSON file", 4];,

						(*Check if the needed nii files exist*)
						If[!(NiiFileExistQ[niiFileP]&&FileExistsQ[ConvertExtension[niiFileP,".bval"]]&&FileExistsQ[ConvertExtension[niiFileP,".bvec"]]),
							(*-----*)AddToLog[{"Skipping, could not find .nii, .bval and .bvec"}, 4],
							(*-----*)AddToLog["Importing the data", 4];

							(*import the data*)
							json = ImportJSON[jsonFileP];
							{data, grad, val, diffVox} = ImportNiiDiff[niiFileP, FlipBvec->False];
							mask = Mask[NormalizeMeanData[data], ConfigLookup[datType, "Process", "Masking"], 
								MaskSmoothing->True, MaskComponents->2, MaskClosing->2];
							data = MaskData[data, mask];

							(*perform gradient nonlinearity correction*)
							coil = ConfigLookup[datType, "Process", "GradientCorrection"];
							settingPro = Join[settingPro, <|"GradientCorrection"->coil|>];
							off = Lookup[json, "Offset", False];
							If[!StringQ[coil], coil = False,
								(*-----*)AddToLog[{"Using Gradient correction: ", coil}, 4];
								If[!VectorQ[off], coil = False;
									(*-----*)AddToLog["No offset defined in json, skipping grad correction" , 5],
									{int, dint} = MakeGradientDerivatives[diffVox, coil];
									coil = {diffVox, off, dint};
									outTypes = Join[{"field"}, outTypes];
								];
							];

							(*perform ivim correction*)
							ivim = ConfigLookup[datType, "Process", "IVIMCorrection"];
							settingPro = Join[settingPro, <|"IVIMCorrection"->ivim|>];
							If[ivim, (*-----*)AddToLog["Starting ivim calculation", 4];
								(*make initial guess for ivim*)
								{mean, bv} = MeanBvalueSignal[data, val];
								sigI = MeanSignal[mean];
								init = IVIMCalc[sigI, bv, {1, .05, .0015, .015}, IVIMFixed -> False, IVIMConstrained -> False];
								debugBids[{"IVIM Init", init}];
								(*actual ivim fitting*)
								{s0i, fri, adci, pD} = If[coil === False,
									(*use normal b-value*)
									IVIMCalc[data, val, init, IVIMConstrained -> True, Parallelize -> True, IVIMFixed -> True],
									(*use gradient field corrected b-value*)
									IVIMCalc[data, {val, grad}, init, coil, IVIMConstrained -> True, Parallelize -> True, IVIMFixed -> True]
								];
								debugBids[{"IVIM fit results", MinMax/@{s0i, fri, adci, pD}}];
								data = First@IVIMCorrectData[data, {s0i, fri, pD}, val, FilterMaps->False];
								adci = 1000 adci;
								outTypes = Join[{"mean", "adci", "fri", "s0i"}, outTypes];
							];

							(*perform the actual tensor calculation*)
							(*-----*)AddToLog["Starting tensor calculation", 4];
							{tens, s0, out} = Quiet@TensorCalc[data, grad, val, coil, FullOutput->True, 
								Method->"iWLLS", RobustFit->True, Parallelize->True, MonitorCalc->False];
							If[coil === False, {out, res} = out, {out, res, field} = out];
							out = Total@Transpose@out;
							(*calculate tensor parameters*)
							{l1, l2, l3, md, fa} = ParameterCalc[tens];
							rd = Mean[{l2, l3}];

							(*perform fasciculation analysis*)
							fasc = ConfigLookup[datType, "Process", "FasciculationDetection"];
							settingPro = Join[settingPro, <|"FasciculationDetection"->fasc|>];
							If[fasc,
								(*-----*)AddToLog["Starting fasciculation analysis", 4];
								(*normalize the diffusion data based on the tensor and select b>200*)
								reg = First@ImportNii[StringReplace[niiFileP, "filt" -> "reg"]];
								sel = First@SelectBvalueData[reg, val, 1];
								norm = NormalizeFasciculationData[sel, mask, {tens, grad, val}];
								(*perform activation analysis*)
								fasc = First@FindActivations[norm, IgnoreSlices -> {0, 0}, 
									ActivationThreshold -> {3., .65}, ActivationOutput -> All, MaskDilation -> 2, 
									ActivationIterations -> 10, ActivationBackground -> 10];
								{fasc, fascm} = SelectActivations[fasc];
								outTypes = Join[{"fasc", "fascm", "norm"}, outTypes];
							];

							(*export all the calculated data*)
							(*-----*)AddToLog["Exporting the calculated data to:", 4];
							(*-----*)AddToLog[outFile, 5];
							debugBids[Column[{json, settingPro}]];
							
							tens = Transpose[tens];
							outTypes = Join[{"data", "tens", "res", "out", "s0", 
								"l1", "l2", "l3", "md", "fa", "rd"}, outTypes];

							(
								ExportNii[ToExpression[con<>#], diffVox, outFile<>"_"<>#<>".nii", CompressNii -> compress];
								If[MemberQ[{"data"}, #], ExportBvalvec[{val, grad}, outFile <> "_"<>#]];
								Export[ConvertExtension[outFile<>"_"<>#, ".json"], MergeJSON[{json, procInfo, settingPro}]];
							) & /@ outTypes;

							(*export the check file*)
							MakeCheckFile[outFile, Sort@Join[
								{"Check"->"done", "Bvalue" -> val, "Gradient" -> grad, 
									"Outputs" -> outTypes, "SetProperties"->set},
								Normal@KeyTake[json, keys]
							]];

							(*compress the nii files if compression during ExportNii -> False*)
							If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

							(*-----*)AddToLog["Finished processing", 3, True];				
						]
					]
				(*close dti processing*)
				],

				_,
				(*-------------------------------------------*)
				(*--------------- Unknown dti ---------------*)
				(*-------------------------------------------*)
				(*-----*)AddToLog[{"Unknown processing ", datType["Process"], "for datatype", type}, True, 3];
			],

			"mese",
			(*-------------------------------------------*)
			(*-------- mese processing scripts ----------*)
			(*-------------------------------------------*)

			Switch[datType["Process", "Method"],				

				"EPGT2"|"Exp",
				(*-------------------------------------------*)
				(*---------- EPG processing script ----------*)
				(*-------------------------------------------*)

				(*input file names*)
				diffFile = GenerateBidsFileName[fol, set];
				jsonFile = ConvertExtension[diffFile, ".json"];
				nfile = ConvertExtension[diffFile, ".nii"];

				(*output file names*)
				outFile = GenerateBidsFileName[folOut, set];

				(*check if files are already done*)
				If[CheckFile[outFile, "done", verCheck],
					(*if check file has label done and version is recent skip*)
					(*-----*)AddToLog["Processing already done for: ", True, 3];
					(*-----*)AddToLog[outFile, 4],
					(*-----*)AddToLog["Starting processing for data:", 3, True];
					(*-----*)AddToLog[diffFile, 4];

					(*Check if needed json Exist*)
					If[!FileExistsQ[jsonFile],
						(*-----*)AddToLog[{"Could not find the needed JSON file of", jsonFile}, 4];,
						(*Check if needed nii Exist*)
						If[!NiiFileExistQ[nfile],
							(*-----*)AddToLog[{"Could not find the data of", diffFile}, 4],
							(*-----*)AddToLog["Importing the data", 4];

							(*import the data*)
							json = ImportJSON[jsonFile];
							echos = json["EchoTime"];
							debugBids["importing: ", nfile];
							{data, t2vox} = ImportNii[nfile];

							If[ArrayDepth[data]=!=4,
								(*-----*)AddToLog["!!!!!!!!!!!!!!! SKIPPING - WRONG DATA DIMENSIONS !!!!!!!!!!!!!!!!!!!", 4];	
								debugBids["!!!!!!!!!!!!! Wrong data dimensions: ", Dimensions@data];
								,

								(*mask the background*)
								mask = Mask[NormalizeMeanData[data], ConfigLookup[datType, "Process", "Masking"], MaskSmoothing -> True, MaskComponents -> 2, MaskClosing -> 2];
								data = NormalizeData@MaskData[data, mask];
								settingPro = <||>;

								Switch[datType["Process", "Method"],

									"EPGT2",

									(*determine the pulse profiles*)
									(*-----*)AddToLog["Calculating the slice profiles", 4];	
									{ex, ref} = ConfigLookup[datType, "Process", "Settings"];
									debugBids[{ex, ref, echos}];
									(*check if the profiles are given as numbers or as files*)
									angle = If[NumberQ[ex] && NumberQ[ref],
										shift = 0.;
										{ex, ref},
										thk = json["SliceThickness"];
										angle = GetPulseProfile[ex, ref, SliceRange -> 3 thk, SliceRangeSamples -> 6 thk][[1;;2]];
										shift = ConfigLookup[datType, "Process", "Shift"];
										shift = If[shift=!=True, 0.,
											shift = (1/ref[[3, 1]] - 1/ex[[3, 1]]) (3(*filed times ppm*) 3.4);
											(*-----*)AddToLog[{"Shifting fat profile with: ", shift}, 4];
											Round[shift/((3 thk/2)/(6 thk))]
										];
										angle
									];

									debugBids[Dimensions/@{data, t2vox}];

									(*calculate the water t2 map, the fit settings are exported with the data*)
									(*-----*)AddToLog["Starting EPG T2 calculation", 4];
									epgOpts = {
										DictT2IncludeWater -> True, EPGFitFat -> False, EPGCalibrate -> True, EPGFatShift -> shift,
										DictT2fValue -> 150, DictT2fRange -> {150, 250, 5},
										DictB1Range -> {0.5, 1.4, 0.02}, DictT2Range -> {15, 45, 0.2}
									};
									{{t2w, t2f, b1}, {wat, fat, fatfr}, res} = EPGT2Fit[data, 1000 echos, angle, MonitorCalc -> False, Sequence @@ epgOpts];
									settingPro = <|"EPGSettings" -> OptionsToJSON[epgOpts]|>;

									outTypes = {"data", "t2w", "t2f", "b1", "wat", "fat", "fatfr", "res"};

									,
									"Exp",

									(*-----*)AddToLog["Starting exponential T2 fitting", 4];
									{s0, t2} = T2Fit[data, 1000 echos];

									outTypes = {"data", "t2", "s0"};
								];

								(*export all the calculated data*)
								(*-----*)AddToLog["Exporting the calculated data to:", 4];
								(*-----*)AddToLog[outFile, 5];		

								(
									ExportNii[ToExpression[con<>#], t2vox, outFile<>"_"<>#<>".nii", CompressNii -> compress];
									Export[ConvertExtension[outFile<>"_"<>#, ".json"], MergeJSON[{json, procInfo, settingPro}]];
								) & /@ outTypes;

								(*export the check file*)
								MakeCheckFile[outFile, Sort@Join[
									{"Check"->"done", "EchoTimes" -> echos, "Outputs" -> outTypes, "SetProperties"->set},
									Normal@KeyTake[json, keys]
								]];

								(*compress the nii files if compression during ExportNii -> False*)
								If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

								(*-----*)AddToLog["Finished processing", 3, True];
							]
						]
					]
				(*close t2 processing*)
				],

				_,
				(*-------------------------------------------*)
				(*--------------- Unknown mese ---------------*)
				(*-------------------------------------------*)
				(*-----*)AddToLog[{"Unknown processing ", datType["Process"], "for datatype", type}, True, 3];
			],
			(*-------------------------------------------*)
			(*------------------ Other ------------------*)
			(*-------------------------------------------*)
			_,
			(*-----*)AddToLog["Unknown datatype for processing", 4];

		(*Close Type switch*)
		];

	(*close loop over sets*)
	, {set, sets}]
]


(* ::Subsection:: *)
(*4. MuscleBidsMerge*)


(* ::Subsubsection::Closed:: *)
(*MuscleBidsMerge*)


Options[MuscleBidsMerge] = {ProcessSubjects->All, VersionCheck->False};

SyntaxInformation[MuscleBidsMerge] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

MuscleBidsMerge[folder_?StringQ, opts:OptionsPattern[]] := MuscleBidsMerge[folder, GetConfig[folder], opts];


MuscleBidsMerge[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir}, 
	debugBids["starting MuscleBidsMerge"];
	dir = Directory[]; SetDirectory[folder];
	MuscleBidsMerge[
		ConfigLookup[config, "folders", "derivedData"],(*the input folder for the data*)
		ConfigLookup[config, "folders", "mergeData"],(*the output folder for merging*)
		config["datasets"],(*what data merging*)
		opts];
	SetDirectory[dir];	
]


MuscleBidsMerge[datFol_?StringQ, merFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := BidsFolderLoop[datFol, merFol, datDis, Method->"MuscleBidsMerge", opts]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsMergeI*)


MuscleBidsMergeI[{folIn_, folOut_}, {datType_, allType_}, verCheck_] := Block[{
		nonQuant, multDim, native, fol, parts, merge, fileCheck, 
		outFile, tarFile, tarStack, movStack, nStack, files, sameType,
		process, processMD, processAll, nSet, processStacks,
		overT, overM, pad, motion, reverse, split, metReg, settings,
		target, vox, voxTar, moving, movingMD, length, lengthMD, voxMov, names,
		posAll, posNat, posScale, voxTarR, voxF, json, 
		im, targetR, mskm, mskt, msktR, func, reg, movp
	},

	(*get the input settings*)
	debugBids["Starting MuscleBidsMergeI"];
	debugBids[{folIn, folOut}];
	debugBids[datType];

	(*muscle bids that are non quantitative of multi dimensional and or have to stay native*)
	nonQuant = {"inph", "outph", "wat", "fat", "s0", "mean"};
	multDim = {"tens", "data", "filt", "fasc", "data", "reg"};
	native = {"tens", "fasc"};

	(*get the outFile names*)
	{fol, parts} = PartitionBidsFolderName[folIn];
	merge = datType["Merging"];

	(*------------------ figure out the file names and check them ------------------*)

	outFile = BuildBidsNameFromConfig[{folOut, parts}, datType];
	debugBids[{parts, outFile}];

	(*target out file and input stack*)
	tarFile = BuildBidsNameFromConfig[{folOut, parts}, {datType, allType}, merge["Target"]];
	tarStack = BuildBidsNameFromConfig[{folOut, fol, parts}, {datType, allType}, merge["Target"]];
	fileCheck = AllTrue[NiiFileExistQ /@ tarStack, # &];
	debugBids[{datType["Class"], tarFile}];
	debugBids["target stack: ", Column@tarStack];
	debugBids[fileCheck];

	If[!fileCheck,
		(*-----*)AddToLog[{"Skipping merging since the target files could not be found"}, 3];
		Return[]];

	(*moving input stack*)
	movStack = BuildBidsNameFromConfig[{fol, parts}, datType, merge["Moving"]];
	fileCheck = AllTrue[NiiFileExistQ /@ tarStack, # &];
	nStack = Length @ movStack;
	debugBids["moving stack: ", Column@movStack];

	If[!fileCheck,
		(*-----*)AddToLog[{"Skipping merging since the moving files could not be found"}, 3];
		Return[]];

	(*process stack*)
	process = Select[merge["Process"], !MemberQ[multDim, #] &];
	processMD = Select[merge["Process"], MemberQ[multDim, #] &];
	length = Length @ process;
	processAll = Join[process, processMD];
	nSet = Length[processAll];
	processStacks = BuildBidsNameFromConfig[{fol, parts}, datType, #]& /@ processAll;
	fileCheck = (files = #; AllTrue[NiiFileExistQ /@ files, # &])& /@ processStacks;
	debugBids["process stack: ", Column[Column/@processStacks]];
	debugBids[fileCheck];

	(*check if number of stacks are consistent*)
	If[!AllTrue[fileCheck, #&],
		(*-----*)AddToLog[{"Not all process files exist:", StringRiffle[ToString/@fileCheck, ", "]}, 3];
		Return[]
	];

	(*start the merging, if check file has label done and version is recent skip*)
	(*-----*)AddToLog[{"Start joining ", nStack, "stacs for", nSet, "dataypes"}, 3];
	(*-----*)AddToLog[{"The types that will be merged are: "}, 3];
	(*-----*)AddToLog[{StringJoin@Riffle[processAll,", "]}, 4];
	If[CheckFile[outFile, "done", verCheck],
		(*-----*)AddToLog[{"Processing already done."}, 3];
		Return[]];

	(*------------------ figure out settings ------------------*)

	(*get the settings for the merging*) 
	overT = ConfigLookup[datType, "Merging", "Overlap"];
	pad = ConfigLookup[datType, "Merging", "Padding"];
	motion = ConfigLookup[datType, "Merging", "Motion"];
	reverse = ConfigLookup[datType, "Merging", "Reverse"];
	split = ConfigLookup[datType, "Merging", "SplitRegistration"];

	If[pad > 0, 
		(*-----*)AddToLog[{"Padding overlap with: ", pad}, 5]];

	{overT, overM} = If[IntegerQ[overT], {overT, overT}, overT];
	If[overT === 0, pad = 0; motion = False];

	(*registration settings*)
	metReg = Switch[datType["InFolder"], 
		"dix", "rigid", 
		"quant", {"rigid", "affine"}, 
		_, {"rigid","affine","bspline"}
	];

	settings = <|
		"Merge Overlap" -> {overT, overM},
		"Merge Padding" -> pad,
		"Merge Motion" -> motion,
		"Merge Reverse" -> reverse,
		"Merge Split" -> split,
		"Merge Method" -> metReg
	|>;

	debugBids["overlap and split: ", {overT, pad, split}];

	(*------------------ Importing off all data after checks ------------------*)

	(*import the Target data, import merged target if it exists else merge it*)
	(*-----*)AddToLog[{"Importing and processing the target data"}, 4];

	(*if target is same as moving never load target from disk*)
	sameType = tarStack === movStack;
	If[NiiFileExistQ[tarFile] && !sameType,
		{target, voxTar} = ImportNii[tarFile];
		(*-----*)AddToLog[{"Splitting the primary datatype that already existed"}, 4];
		If[nStack=!=1, 
			target = SplitSets[target, nStack, overT, ReverseSets->reverse, PadOverlap->pad], 
			target = {target}
		];
		,
		(*remake target*)
		{target, vox} = Transpose[ImportNii /@ tarStack];
		voxTar = First@vox;
		debugBids["voxel size per stack: ", {vox, voxTar}];
		(* make data real valued*)
		target = If[RealValuedNumberQ[target[[1, 1, 1, 1]]], target, Abs[target]];

		(*-----*)AddToLog[{"Joining the primary datatype", If[motion,"with","without"], "motion correction"}, 4];
		If[nStack=!=1,
			target = JoinSets[target, overT, voxTar, ReverseSets->reverse, MotionCorrectSets->motion, 
				NormalizeSets->True, NormalizeOverlap->True, MonitorCalc->False];
			debugBids["target dimensions: ", Dimensions@target];
			target = SplitSets[target, nStack, overT, ReverseSets->reverse, PadOverlap->pad];
		];
	];
	debugBids["target dimensions: ", Dimensions/@target];

	(*import the Moving data, import the multi dimm separately*)
	(*-----*)AddToLog[{"Importing and processing the moving data"}, 4];

	(*import the 3D moving data *)
	{moving, vox} = Transpose[(files=#; Transpose[ImportNii /@ files]) &/@ processStacks];

	voxMov = First@vox; (*voxel size per stack*)
	(*check voxel sizes of moving*)
	If[!Equal@@voxMov,
		(*-----*)AddToLog[{"********** The voxel size is not the same for all stacks **********"}, 0];
		(*-----*)AddToLog[{"Voxel size per stack: ", voxMov}, 4];
		(*-----*)AddToLog[{"", (Dimensions /@ First@moving) voxMov}, 4];
	];

	(*multi dim data has to be flattened*)
	movingMD = If[processMD==={}, lengthMD=0; {},	
		movingMD = moving[[length+1;;]];
		lengthMD = Length /@ movingMD[[All, 1, 1]]; (*4th dimensions of the multi dims*)
		Flatten[movingMD, {1, 4}]
	];
	names = Flatten[{process, ConstantArray[#[[1]], #[[2]]] & /@ Thread[{processMD, lengthMD}]}];

	(*redefine moving - make data real valued to prevent errors*)
	moving = Join[moving[[;;length]], movingMD];
	moving = If[RealValuedNumberQ[N@#[[1,1,1,1]]], #, Abs[#]] &/@ moving;

	debugBids["moving dimensions: ", Dimensions/@{moving, moving[[;;length]], movingMD}];
	debugBids["number of sets: ", {length, lengthMD}];
	debugBids["names: ", names];

	(*figure out the scaling of each contrast*)
	posAll = Range[Length[names]];
	posNat = Flatten[Position[names, _?(MemberQ[native, #] &), {1}]];
	posScale = Complement[posAll, posNat];
	voxTarR = First@voxMov;
	voxF = If[MemberQ[native, names[[#]]], voxTarR, voxTar]&;

	debugBids["moving before registr: ", Column[Dimensions /@ # & /@ moving]];
	debugBids["positions: ", {posScale, posNat}];

	(*import the json information*)
	json = (files=#; MergeJSON[ImportJSON[ConvertExtension[#, ".json"]]&/@files]) &/@ processStacks;
	debugBids["json dimensions: ", Dimensions@json];

	(*------------------ performing registration to target for all contrasts ------------------*)

	(*perform motion correction after target merging*)
	(*If motion correction for joning is False and target is of same type no need for motion correction*)
	debugBids["go into registration: ", {motion, sameType, !(!motion && sameType)}];

	If[!(!motion && sameType), 
		(*-----*)AddToLog[{"Performing the registration for the all the datasets"}, 4];
		im = First@FirstPosition[names, merge["Moving"]];
		debugBids["location moving contrast: ", im];

		moving = Table[
			(*-----*)AddToLog[{"Stack: ", i}, 5];
			debugBids["start registration: ", i];

			(*make target in native space*)
			targetR = RescaleData[target[[i]], {voxTar, voxTarR}];
			debugBids[Dimensions/@{moving[[im, i]], target[[i]], targetR}];

			(*make masks*)
			mskm = DilateMask[Mask[NormalizeData[moving[[im, i]]], 5], 5];
			mskt = DilateMask[Mask[NormalizeData[target[[i]]], 5], 5];
			msktR = DilateMask[Mask[NormalizeData[targetR], 5], 5];

			(*only register if not the same contrast and not the first stack*)
			debugBids["check for registration: ", {i===If[reverse, nStack, 1], sameType, i===If[reverse, nStack, 1]&& sameType}];
			
			If[i===If[reverse, nStack, 1] && sameType,
				(*just select data and do nothing*)
				debugBids["no registration: ", i];
				moving[[All, i]]
				,
				(*move the target from anatomical to native space*)
				func = If[i===If[reverse, nStack, 1] || split===False, 
					RegisterData, RegisterDataSplit];
				reg = ToPackedArray@N@Chop@func[
					{moving[[im, i]], mskm, voxMov[[i]]}, {target[[i]], voxTar}, 
					Iterations->300, BsplineSpacing->20 voxTar, InterpolationOrderReg->1, 
					NumberSamples -> 20000, PrintTempDirectory->False, 
					MethodReg->metReg, HistogramBins -> 128];

				(*if padding enlarge the moving files*)
				If[pad > 0,					
					reg = ArrayPad[reg, Ceiling@{{pad, pad}/2, 0, 0}];
					movp = Transpose[ArrayPad[#, Ceiling@{{pad, pad}/2, 0, 0}]&/@moving[[All, i]]],
					movp = Transpose[moving[[All, i]]]
				];
				debugBids["Dimension before registration: ", Dimensions@movp, Dimensions/@Transpose[movp]];

				(*register back the target from native space to anatomy and transform the rest*)
				func = If[i===If[reverse, nStack, 1] || split===False, 
					RegisterDataTransform, RegisterDataTransformSplit];
				movp = Flatten[(
					debugBids["registering: ", #[[1]]];
					If[#[[1]]==={}, {},
						Transpose@ToPackedArray@N@Chop@Last@func[
							#[[2]], {reg, voxMov[[i]]}, {movp[[All, #[[1]]]], voxMov[[i]]},
							Iterations->300, BsplineSpacing->30 voxTar, InterpolationOrderReg->1, 
							NumberSamples -> 20000, PrintTempDirectory->False, DeleteTempDirectory->True, 
							MethodReg->metReg, BsplineDirections -> {1, 1, 0}, HistogramBins -> 128
						]
					])& /@ {{posScale, {target[[i]], mskt, voxTar}}, {posNat, {targetR, msktR, voxTarR}}}
				, 1][[Ordering[Join[posScale, posNat]]]];
				debugBids["Dimensions after registration: ", Dimensions@movp, Dimensions/@movp];
				movp
			]
		, {i, 1, nStack}];

		(*extract all parameters after registration*)
		debugBids["after dimensions before transpose: ", Dimensions/@moving];
		moving = Transpose[moving];
	];(*close motion moving*)

	(*------------------ actually join the data and reorder for multi dim ------------------*)

	debugBids["after dimensions: ", Dimensions@moving];
	debugBids["after registration: ", Column[Dimensions /@ # & /@ moving]];

	(*join the moving types*)
	(*-----*)AddToLog[{"Joining the data"}, 4];
	moving = If[nStack===1,
		moving[[All, 1]],
		debugBids["joining: ", posAll];
		JoinSets[
			If[MemberQ[nonQuant, names[[#]]], Ramp[moving[[#]]], moving[[#]]], 
			If[MemberQ[posNat, #], overM, overT], voxF[#], 
			MonitorCalc->False, MotionCorrectSets->False, 
			PadOverlap->pad, ReverseSets->reverse, 
			NormalizeSets->MemberQ[nonQuant, names[[#]]], 
			NormalizeOverlap->MemberQ[nonQuant, names[[#]]]
		]&/@posAll
	];

	debugBids["after merging: ", Column[Dimensions /@ moving]];		

	(*collect the multi dim data and unflatten them*)
	If[movingMD =!= {},
		movingMD = moving[[length+1;;]];
		movingMD = Transpose[movingMD[[#[[1]];;#[[2]]]]] & /@ (
			{1, 0} + # & /@ Partition[Prepend[Accumulate[lengthMD], 0], 2, 1]);
		moving = Join[moving[[1;;length]], movingMD];
	];

	(*------------------ exporting the data ------------------*)

	(*export the joined data with the merged json*)
	(*-----*)AddToLog["Exporting the calculated data to:", 4];
	(*-----*)AddToLog[outFile, 5];
	(
		debugBids["Exporting: ", {processAll[[#]], voxF[#]}];
		ExportNii[moving[[#]], voxF[#], outFile<>"_"<>processAll[[#]]<>".nii", CompressNii -> compress];
		Export[outFile<>"_"<>processAll[[#]]<>".json", MergeJSON[{json[[#]], settings}]];
	)&/@ Range[nSet];

	(*make the check file*)
	MakeCheckFile[outFile, Sort@Join[{"Check"->"done"}, Normal@datType]];

	(*compress the nii files if compression during ExportNii -> False*)
	If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

	(*-----*)AddToLog["Finished merging", 3, True];
]


(* ::Subsection:: *)
(*5. MuscleBidsSegment*)


(* ::Subsubsection::Closed:: *)
(*MuscleBidsSegment*)


Options[MuscleBidsSegment] = {ProcessSubjects->All, VersionCheck->False};

SyntaxInformation[MuscleBidsSegment] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

MuscleBidsSegment[folder_?StringQ, opts:OptionsPattern[]] := MuscleBidsSegment[folder, GetConfig[folder], opts];


MuscleBidsSegment[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir},
	debugBids["starting MuscleBidsSegment"]; 
	dir = Directory[]; SetDirectory[folder];
	MuscleBidsSegment[
		ConfigLookup[config, "folders", "mergeData"],(*the input folder for the data, output is in same folder*)
		ConfigLookup[config, "folders", "mergeData"],
		config["datasets"],(*what data for segmentation*)
		opts];
	SetDirectory[dir];	
]


MuscleBidsSegment[datFol_?StringQ, outFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := BidsFolderLoop[datFol, outFol, datDis, Method->"MuscleBidsSegment", opts]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsSegmentI*)


MuscleBidsSegmentI[{folIn_, folOut_}, {datType_, allType_}, verCheck_] := Block[{
		segType, segTypeLab, checkFile, fol, device, location, segment,
		parts, outFile, segFile, out, vox, seg, status, segDim,
		voxS, tarI, movI, segI, tar, mov, dim, voxTar, voxMov, voxSeg, mask,
		procInfo, segOpts, regOpts, srcFiles
	},

	status = "done";

	(*config settings and version that are added to every output json*)
	procInfo = JSONSafe@<|"SegmentSettings" -> datType["Segment"], "ProcessVersion" -> QMRITools`$InstalledVersion|>;

	debugBids["Starting MuscleBidsSegmentI"];
	debugBids[datType];
	{fol, parts} = PartitionBidsFolderName[folIn];
	debugBids[folIn, folOut];
	debugBids[{fol, parts}];

	(*Check if segmentation is defined*)
	If[datType["Segment"] === Missing["KeyAbsent", "Segment"],	
		(*-----*)AddToLog[{"No Segmentations defined for this data"}, 3];
		debugBids["No segmentation defined"];
		Return[]
	];

	(*Check if segmentation needs to be performed*)
	checkFile = ConvertExtension@BuildBidsNameFromConfig[{folOut, parts}, {datType, allType}, 
		{datType["Key"], "seg", "auto", datType["Type"]}];
	debugBids[checkFile];
	If[CheckFile[checkFile, "done", verCheck],
		(*-----*)AddToLog[{"Segmentation already done for:", datType["key"]}, 3];
		Return[]
	];

	(*Get the segmentation targets and its names*)
	segType = ConfigLookup[datType, "Segment", "Target"];
	debugBids[segType];

	If[ArrayDepth[segType]===1, segType = {segType}];
	Switch[ConfigLookup[datType, "Segment", "Method"],

		(*Automatic segmentations using NN*)
		Automatic,

		(*Loop over the segmentation types if more are specified*)
		Table[
			(*-----*)AddToLog[{"Performing segmentation for ", StringRiffle[segI, "_"]}, 3];

			(*get the correct filenames*)
			outFile = BuildBidsNameFromConfig[{folOut, parts}, {datType, allType},
				Join[{datType["Key"], "seg", "auto"}, segI]];
			segFile = BuildBidsNameFromConfig[{fol, parts}, {datType, allType}, segI];

			(*check if target file exists if so perform the segmentation*)
			If[!NiiFileExistQ[segFile],
				(*-----*)AddToLog[{"The segmentation file does not exist: ", segFile}, 4];
				status = "error"
				,
				{out, vox} = ImportNii[segFile];

				segDim = ConfigLookup[datType, "Segment", "Dimensions"];
				If[segDim=!="2D" && segDim=!="3D", segDim="3D";];
				location = ConfigLookup[datType, "Segment", "Location"];
				(*-----*)AddToLog[{"Segmenting location using dimensions: ", location, " - ", segDim}, 4];
				debugBids[{segDim, location, vox}];

				segOpts = {SegmentationDimension -> segDim, TargetDevice -> ConfigLookup[datType, "Segment", "Device"]};
				seg = SegmentData[{out, vox}, location, Monitor -> False, Sequence @@ segOpts];

				(*export the segmentation with a json of its source, settings and networks*)
				ExportNii[seg, vox, outFile, CompressNii -> compress];
				Export[ConvertExtension[outFile, ".json"], MergeJSON[{SourceJSON[segFile], procInfo, <|
						"SegmentationMethod" -> "Automatic", "SegmentationSource" -> FileNameTake[ConvertExtension[segFile, ".nii"]],
						"SegmentData" -> OptionsToJSON[segOpts],
						"SegmentationNetworks" -> JSONSafe@QMRITools`SegmentationTools`Private`SegmentationInfo[location, segDim]
				|>}]];
			];
			, {segI, segType}			
		],

		(*Register existing segmentation*)
		"Registration",

		(*figure out if duplicate handling is needed.*)
		segment = datType["Segment"];
		outFile = BuildBidsNameFromConfig[{folOut, parts}, {datType, allType},
				Join[{datType["Key"], "seg", "reg"}, segment["Target"]]];
		tar = BuildBidsNameFromConfig[{fol, parts}, {datType, allType}, segment["Target"]];
		mov = BuildBidsNameFromConfig[{fol, parts}, {datType, allType}, segment["Moving"]];
		seg = BuildBidsNameFromConfig[{fol, parts}, {datType, allType}, segment["Segmentation"]];
		srcFiles = FileNameTake[ConvertExtension[#, ".nii"]]& /@ <|"Target" -> tar, "Moving" -> mov, "Segmentation" -> seg|>;
		debugBids[Column[{outFile, seg, tar, mov}]];

		(*Import and prepare the data*)
		{tar, voxTar} = ImportNii[tar];
		tar = NormalizeData[Clip[tar, {0, Infinity}]];
		{mov, voxMov} = ImportNii[mov];
		mov = NormalizeData[Clip[mov, {0, Infinity}]];
		{seg, voxSeg} = ImportNii[seg];
		mask = Mask[NormalizeData[tar], 5, MaskSmoothing -> True, MaskDilation -> 5];
		voxS = ConfigLookup[datType, "Segment", "VoxSize"];
		If[voxS =!= Automatic, mov = RescaleData[mov, {voxMov, voxS}]; voxMov = voxS;];
		debugBids[Column[Dimensions/@{seg, tar, mov}]];

		(*perform the segmentation registration*)
		regOpts = {MethodReg -> {"rigid", "affine"}, NumberSamples -> 50000, HistogramBins -> 128, Resolutions -> 1, Iterations -> 200};
		seg = Table[RegisterDataTransform[
			{tar[[{i}]], mask[[{i}]], voxTar}, {mov, voxMov}, {seg, voxSeg}, Sequence @@ regOpts,
			PrintTempDirectory -> False, DeleteTempDirectory -> True, TransformMethod -> "Segmentation"
		][[2]], {i, 1, 3}];

		(*export the segmentation with a json of its sources and settings*)
		ExportNii[seg, voxTar, outFile, CompressNii -> compress];
		Export[ConvertExtension[outFile, ".json"], MergeJSON[{procInfo, <|
			"SegmentationMethod" -> "Registration", "SegmentationSource" -> srcFiles,
			"Registration" -> OptionsToJSON[regOpts], "VoxSize" -> voxMov
		|>}]];
	];

	(*make the check file*)
	MakeCheckFile[checkFile, Sort@Join[{"Check" -> status}, Normal@datType]];

	(*compress the nii files if compression during ExportNii -> False*)
	If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

	(*-----*)AddToLog["Finished the segmentation", 3, True];
]


(* ::Subsection:: *)
(*6. MuscleBidsTractography*)


(* ::Subsubsection::Closed:: *)
(*MuscleBidsTractography*)


Options[MuscleBidsTractography] = {ProcessSubjects->All, VersionCheck->False, BidsTractographyMethod->"Full"};

SyntaxInformation[MuscleBidsTractography] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

MuscleBidsTractography[folder_?StringQ, opts:OptionsPattern[]] := MuscleBidsTractography[folder, GetConfig[folder], opts];


MuscleBidsTractography[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir},
	debugBids["starting MuscleBidsTractography"];
	dir = Directory[]; SetDirectory[folder];
	MuscleBidsTractography[
		ConfigLookup[config, "folders", "mergeData"],(*the input folder for the data, output is in same folder*)
		ConfigLookup[config, "folders", "mergeData"],
		config["datasets"],(*what data for tractography*)
		opts];
	SetDirectory[dir];	
]


MuscleBidsTractography[datFol_?StringQ, outFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := BidsFolderLoop[datFol, outFol, datDis, Method->"MuscleBidsTractography", opts]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsTractographyI*)


MuscleBidsTractographyI[{folIn_, folOut_}, {datType_, allType_}, verCheck_, met_] := Block[{
		tracto, tractType, tractSeg, tractStopLab, tractStopVal, outFile,
		tractTypeLab, fol, parts, checkFile, seed, lenS, segBone, tractSegLab,
		datfile, stopFile, tens, vox, dim, stop, ang, step, tracts, seeds, len, seg, curv,
		segFile, muscles, muscLabs, mus, bones, context, leng, dens, flip, per, 
		voxSeg, dims, harm, tensH, con, amp, trkFileF,
		procInfo, trackOpts, harmOpts, tsegOpts, trackInfo, harmInfo, tsegInfo
	},

	debugBids["Starting MuscleBidsTractographyI"];
	debugBids[datType];
	{fol, parts} = PartitionBidsFolderName[folIn];
	debugBids[{folIn, folOut}];
	debugBids[{fol, parts}];

	(* Extract tractography settings*)
	tracto = datType["Tractography"]; 
	If[tracto===Missing["KeyAbsent", "Tractography"],
		(*-----*)AddToLog[{"No tractography defined for this data"}, 3];
		Return[];
	];

	(*config settings and version that are added to every output json*)
	procInfo = JSONSafe@<|"TractographySettings" -> tracto, "ProcessVersion" -> QMRITools`$InstalledVersion|>;

	(* Extract tractography settings*)
	{flip, per} = ConfigLookup[datType, "Tractography", "FlipPermute"];
	len = ConfigLookup[datType, "Tractography", "TractLength"];
	ang = ConfigLookup[datType, "Tractography", "TractAngle"];
	step = ConfigLookup[datType, "Tractography", "TractStep"];
	seed = ConfigLookup[datType, "Tractography", "TractSeed"];
	lenS = ConfigLookup[datType, "Tractography", "SegmentLength"];
	segBone = ConfigLookup[datType, "Tractography", "BoneLabel"];
	harm = ConfigLookup[datType, "Tractography", "HarmonicDenoise"];

	(*tract files and nii*)
	tractType = tracto["Target"];
	tractTypeLab = StringRiffle[tractType, "_"];

	checkFile = outFile = ConvertExtension@BuildBidsNameFromConfig[{folOut, parts}, 
		{datType, allType},	Join[tractType[[;;-2]], {"trk"}]];
	datfile = BuildBidsNameFromConfig[{folOut, parts}, {datType, allType}, tractType];

	(* stopping, ensure tractStopLab and tractStopVal are lists *)
	{tractStopLab, tractStopVal} = Transpose@tracto["Stopping"];
	If[ArrayDepth[tractStopLab]===1, tractStopLab = {tractStopLab}];
	If[ArrayDepth[tractStopVal]===1, tractStopVal = {tractStopVal}];
	stopFile = BuildBidsNameFromConfig[{folOut, parts}, {datType, allType}, #]& /@ tractStopLab;

	(*optional segmentation file*)
	tractSeg = tracto["Segmentation"];
	segFile = If[!ListQ[tractSeg], "No Seg File",
		tractSegLab = StringRiffle[tractSeg, "_"];
		BuildBidsNameFromConfig[{folOut, parts}, {datType, allType}, tractSeg]];

	(* If tractography and segmentation is already done, log the event *)
	debugBids[Column@{outFile, datfile, Column@stopFile, segFile}];
	If[CheckFile[checkFile, "done", verCheck],
		(*-----*)AddToLog[{"Tractography and segmentation already done for:", tractTypeLab}, 3];
		Return[]];

	(* Check if tractography needs to be performed *)
	Which[
		CheckFile[checkFile, "track", verCheck],
		(* If tractography is already done, log the event *)
		(*-----*)AddToLog[{"Tractography already done for:", tractTypeLab}, 3];
		,
		!(met === "Full" || met==="Tractography"),
		(*-----*)AddToLog[{"Skipping tractography because of method:", met}, 3];
		,
		True,
		(* If tractography is not done, log the event and proceed with the processing*)
		(*-----*)AddToLog[{"Starting the whole volume tractography"}, 3, True];
		(*-----*)AddToLog[{"The type that will be tracted is: ", tractTypeLab}, 4];

		(* Check if the tensor file and stop files exist *)
		Which[
			!NiiFileExistQ[datfile],
			(*-----*)AddToLog[{"The tensor file does not exist", datfile}, 4];
			,
			!And@@(NiiFileExistQ/@stopFile),
			(*-----*)AddToLog[{"Not all stop files exist not exist", stopFile}, 4];
			,
			True,
			(* If all files exist, proceed with the tractography processing*)
			(*-----*)AddToLog[{"Importing the needed data"}, 4];
			{tens, vox} = ImportNii[datfile];
			tens = Transpose@ToPackedArray@N@tens;
			dim = Rest@Dimensions@tens;

			(* Import stop files *)
			stop = (
				{stop, voxSeg} = ImportNii[#];
				debugBids["Voxel sizes: ", {voxSeg, vox , Dimensions@stop}];
				If[voxSeg===vox, stop, PadToDimensions[RescaleData[stop, {voxSeg, vox}], dim]]
			)& /@ stopFile;
			debugBids["Data dimensions: ", Dimensions/@Join[stop, {tens}]];
			stop = Transpose[{stop, tractStopVal}];
			

			(* Perform tractography, the settings are exported with the tracts *)
			(*-----*)AddToLog[{"Starting the whole volume tractography"}, 4];
			trackOpts = {
				InterpolationOrder -> 0, StepSize -> step, TractMethod -> "RK4",
				MaxSeedPoints -> If[seed<1, Scaled[seed], seed], FiberLengthRange -> len, FiberAngle -> ang
			};
			{tracts, seeds} = FiberTractography[tens, vox, stop, Sequence @@ trackOpts,
				TensorFlips -> flip, TensorPermutations -> per, TractMonitor -> False, Parallelization -> True
			];
			trackInfo = <|
				"TractSource" -> FileNameTake[ConvertExtension[datfile, ".nii"]],
				"StopFiles" -> (FileNameTake[ConvertExtension[#, ".nii"]]& /@ stopFile), "StopValues" -> tractStopVal,
				"Tractography" -> OptionsToJSON[Join[trackOpts, {TensorFlips -> flip, TensorPermutations -> per}]]
			|>;

			(* Export the tractography results *)
			(*-----*)AddToLog[{"Exporting the whole volume tractography"}, 4];
			ExportTracts[outFile<>".trk", tracts, vox, dim, seeds];
			Export[outFile<>".json", MergeJSON[{SourceJSON[datfile], procInfo, trackInfo}]];

			(*perform harmonic denoising if needed*)
			Which[
				!NiiFileExistQ[segFile] && harm,
				(*-----*)AddToLog[{"The segmentation file does not exist which is needed for harmonic denoise: ", segFile}, 4];
				,
				harm,
				(*-----*)AddToLog[{"Performing harmonic tensor denoising"}, 4];
				(* get the segmentation data *)			
				{seg, voxSeg} = ImportNii[segFile];
				seg = RescaleSegmentation[seg, {voxSeg, vox}];

				{muscles, muscLabs} = SelectSegmentations[seg, Range[segBone], False];
				harmOpts = {MaxIterations -> 250, GradientStepSize -> {2, 1}, RadialBasisKernel -> 12, TensorFlips -> flip, TensorPermutations -> per};
				{tensH, con, amp} = HarmonicDenoiseTensor[tens, seg, vox, muscLabs, Sequence @@ harmOpts, Monitor -> False, Parallelize -> True];

				(*-----*)AddToLog[{"Starting the whole volume tractography on harmonic de-noised data"}, 4];
				stop = {{Dilation[Normal@Total@Transpose@muscles, 1], {0.9, 1.1}}};
				{tracts, seeds} = FiberTractography[tensH, vox, stop, Sequence @@ trackOpts, TractMonitor -> False, Parallelization -> True];
				harmInfo = <|
					"TractSource" -> FileNameTake[ConvertExtension[datfile, ".nii"]], "SegmentationSource" -> FileNameTake[ConvertExtension[segFile, ".nii"]],
					"HarmonicDenoise" -> OptionsToJSON[harmOpts], "StopMask" -> "muscle segmentation dilated by 1 voxel, range {0.9, 1.1}",
					"Tractography" -> OptionsToJSON[trackOpts]
				|>;

				(* Export the tractography results with their settings *)
				(*-----*)AddToLog[{"Exporting the whole volume tractography"}, 4];
				context = Context[context];
				tensH = Transpose@tensH;
				(
					ExportNii[ToExpression[context<>#], vox, outFile<>"_"<>#<>".nii.gz", CompressNii -> compress];
					Export[outFile<>"_"<>#<>".json", MergeJSON[{SourceJSON[datfile], procInfo, harmInfo}]]
				)& /@ {"con", "amp", "tensH"};
				ExportTracts[outFile<>"_har.trk", tracts, vox, dim, seeds];
				Export[outFile<>"_har.json", MergeJSON[{SourceJSON[datfile], procInfo, harmInfo}]];
			];

			MakeCheckFile[checkFile, Sort@Join[{"Check" -> "track"}, Normal@datType]];

			(*compress the nii files if compression during ExportNii -> False*)
			If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

			(*-----*)AddToLog["Finished the tractography", 3, True];
		];
	];

	(* Check if segmentation needs to be performed *)
	trkFileF = If[harm, outFile<>"_har.trk", outFile<>".trk"];
	Which[
		CheckFile[checkFile, "seg", verCheck],
		(* If segmentation or tractography is already done, log the event *)
		(*-----*)AddToLog[{"Segmentation of tractography already done for:", tractSegLab}, 3];
		,
		!(met === "Full" || met==="Segmentation"),
		(*-----*)AddToLog[{"Skipping tractography segmentation because of method:", met}, 3];
		,
		True,
		(* If segmentation or tractography is not done, log the event and proceed with the processing *)
		(*-----*)AddToLog[{"Starting the tractography segmentation"}, 3, True];
		(*-----*)AddToLog[{"The tractography will be segmented using: ", tractSegLab}, 4];

		Which[
			!NiiFileExistQ[segFile],
			(*-----*)AddToLog[{"The segmentation file does not exist: ", segFile}, 4];
			,
			!FileExistsQ[trkFileF],
			(*-----*)AddToLog[{"The tracts file does not exist: ", trkFileF}, 4];
			,
			True,

			(*-----*)AddToLog[{"Importing the needed data"}, 4];
			(*import trk file if needed, if processing was done in same run this is skipped*)
			If[Dimensions[tracts] === {}, {tracts, vox, dim, seeds} = ImportTracts[trkFileF]];
			(* get the segmentation data *)			
			{seg, voxSeg} = ImportNii[segFile];
			{muscles, muscLabs} = SelectSegmentations[seg, Range[segBone], False];
			bones = Unitize[SelectSegmentations[seg, Range[segBone + 1, segBone + 30]]];
			dims = Dimensions@seg;

			debugBids[{{dims,voxSeg, dims voxSeg}, {dim, vox, dim vox}}];

			(*-----*)AddToLog[{"Segmenting the tracts"}, 4];
			(*perform fitting and segmentations of the tracts, the settings are exported with the maps*)
			tsegOpts = {FiberLengthRange -> lenS, FitTractSegments -> True};
			tracts = SegmentTracts[tracts, muscles, voxSeg, dims, Sequence @@ tsegOpts];
			tsegInfo = <|
				"TractSource" -> FileNameTake[trkFileF], "SegmentationSource" -> FileNameTake[ConvertExtension[segFile, ".nii"]],
				"BoneLabel" -> segBone, "SegmentTracts" -> OptionsToJSON[tsegOpts]
			|>;

			(*-----*)AddToLog[{"Analyzing the tracts"}, 4];
			(*Calculate tract parameters*)
			seed = SeedDensityMap[seeds, voxSeg, dims];
			dens = TractDensityMap[tracts, voxSeg, dims];
			leng = TractLengthMap[tracts, voxSeg, dims];
			ang = TractAngleMap[tracts, voxSeg, dims];
			curv = TractCurvatureMap[tracts, voxSeg, dims];

			(*-----*)AddToLog[{"Exporting the results and maps"}, 4];
			(*export stuff*)
			context = Context[context];
			(
				ExportNii[ToExpression[context<>#], voxSeg, outFile<>"_"<>#<>".nii.gz", CompressNii -> compress];
				Export[outFile<>"_"<>#<>".json", MergeJSON[{SourceJSON[segFile], procInfo, tsegInfo}]]
			)& /@ {"dens", "leng", "ang", "seed","curv"};
			ExportTracts[outFile<>"_seg.trk", tracts, voxSeg, dims, seeds];
			Export[outFile<>"_seg.json", MergeJSON[{SourceJSON[segFile], procInfo, tsegInfo}]];

			(*export plot scene*)
			(*-----*)AddToLog[{"Exporting the scene"}, 4];
			Export[outFile<>"_plot.wxf",
				PlotSegmentedTracts[tracts, muscles, bones, dims, voxSeg, 
					OutputForm -> "All", Method -> "tube", MaxTracts -> 10000]
			];

			MakeCheckFile[checkFile, Sort@Join[{"Check"->"seg"}, Normal@datType]];

			(*compress the nii files if compression during ExportNii -> False*)
			If[!compress, CompressNiiFiles[DirectoryName[outFile]]];

			(*-----*)AddToLog["Finished the tractography segmentation", 3, True];
		];
	];

	If[CheckFile[checkFile, "seg", verCheck],
		MakeCheckFile[checkFile, Sort@Join[{"Check"->"done"}, Normal@datType]];

		(*compress the nii files if compression during ExportNii -> False*)
		If[!compress, CompressNiiFiles[DirectoryName[outFile]]];
	];
]


(* ::Subsection:: *)
(*7. MuscleBidsAnalysis*)


(* ::Subsubsection::Closed:: *)
(*MuscleBidsAnalysis*)


Options[MuscleBidsAnalysis] = {ProcessSubjects->All, VersionCheck->False, BidsOutputImages->"All"};

SyntaxInformation[MuscleBidsAnalysis] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

MuscleBidsAnalysis[folder_?StringQ, opts:OptionsPattern[]] := MuscleBidsAnalysis[folder, GetConfig[folder], opts];


MuscleBidsAnalysis[folder_?StringQ, config_?AssociationQ, opts:OptionsPattern[]] := Block[{dir}, 
	debugBids["starting MuscleBidsAnalysis"];
	dir = Directory[]; SetDirectory[folder];
	MuscleBidsAnalysis[
		ConfigLookup[config, "folders", "mergeData"],(*the input folder for the data*)
		ConfigLookup[config, "folders", "analysis"],(*the output folder for analysis*)
		config["analysis"],(*what data for analysis*)
		opts];
	SetDirectory[dir];	
]


MuscleBidsAnalysis[datFol_?StringQ, anFol_?StringQ, datDis_?AssociationQ, opts:OptionsPattern[]] := Block[{data, name},
	(*loop over all folders*)
	BidsFolderLoop[datFol, anFol, datDis, Method->"MuscleBidsAnalysis", opts];

	(*processing for joining all generated datafiles*)
	data = Join @@ (Import /@ Select[FileNames["*.wxf", anFol, Infinity], 
		First[StringSplit[FileBaseName[#], "_"]] =!= "All" &]);
	name = FileNameJoin[{anFol, "All_"<>DateName[]}];

	(*export to summary data file*)
	Export[name <> ".xlsx", data];
	Export[name <> ".wxf", data];
]


(* ::Subsubsection::Closed:: *)
(*MuscleBidsAnalysisI*)


bidsLabelRange = <|"Body" -> 799, "Legs" -> 90, "Shoulder" -> 30, "Arm" -> 12|>;

BidsAnalysisLabels[lab_] := Block[{items, nr, wh},
	(*the labels setting is a number, a word or both, in any order*)
	items = Flatten[{lab}];
	nr = Cases[items, _Integer];
	wh = Cases[items, _String];
	If[items === {} || Length[nr] > 1 || Length[wh] > 1 || Length[nr] + Length[wh] =!= Length[items],
		Return[$Failed, Block]
	];
	(*the word defaults to Body, the number to the range of the word or to all labels*)
	wh = If[wh === {}, "Body", First[wh]];
	{If[nr === {}, Lookup[bidsLabelRange, wh, Automatic], First[nr]], wh}
]


BidsMaskSeg[seg_, erosion_, trMask_] := Block[{segE},
	(*erode the masks*)
	segE = If[erosion, DilateMask[seg, -1], seg];
	(*masks and masks restricted to the tract density*)
	{segE, If[trMask =!= 1, MaskSegmentation[segE, trMask], segE]}
]


BidsImportSeg[file_, musNr_, dim_] := Block[{seg, vox, masks, labs, zero, mat, pos},
	(*import the segmentation*)
	{seg, vox} = ImportNii[file];
	seg = If[ArrayDepth[seg] === 4, seg[[All, 1]], seg];
	(*grid has to be the same as the default segmentation*)
	If[Dimensions[seg] =!= dim, Return[$Failed, Block]];
	(*select the labels of the default segmentation*)
	{masks, labs} = SelectSegmentations[seg, musNr, False];
	mat = If[labs === {}, {}, Transpose[masks]];
	zero = SparseArray[{}, dim];
	(*the mask of each default label, empty for the labels this segmentation does not have*)
	pos = FirstPosition[labs, #, {0}][[1]] & /@ musNr;
	{Transpose[If[# === 0, zero, mat[[#]]] & /@ pos], Pick[musNr, pos, 0]}
]


MuscleBidsAnalysisI[{folIn_, folOut_}, datDis_, verCheck_, imOut_] := Block[{
		maskErosion, tractWeighting, anaSeg, fol, parts, segFile, fileName, partsO, fileNameO, checkFileX, checkFileI,
		n, what, seg, vox, vol, musNr, musName, sideName, sideNr, dataLabs, anaType, densLab, str,
		densFile, trType, trMask, segT,	datfile, data, scale, tract, outFile, meanType, hasKey,
		quantIm, segIm, tractIm, imRef, ref, crp, refC, size, pos, sliceData, make3DImage, make2DImage,
		cols, cFun, ran, clip, type, imFile, imDat, voxi, voxSeg, segPl, imTrk, trkfile, reffile,
		addLabel, img, lab, label, opts, filt, filtMask, cross, volF, crossF, typeSuf,
		whatName, whatNr, segF, labAll, labCfg, nMax, typeSeg, segFileO, segMiss, segO, segDim, segKey, sg, sgT, sgF
	},

	debugBids["Starting MuscleBidsAnalysisI"];
	debugBids[{folIn, folOut}];

	(*Options*)
	maskErosion = ConfigLookup[datDis, "Analysis", "MaskErosion"];
	tractWeighting = False;

	(*----------- make the xls files -------------*)

	(*get the segmentation settings*)

	debugBids[datDis];

	hasKey = !(datDis["Key"]==="Default");
	{fol, parts} = PartitionBidsFolderName[folIn];
	typeSuf = If[datDis["Analysis", "Class"] === "Stacks", "stk", "chunk"];
	
	partsO = If[hasKey, Join[<|typeSuf -> StringStrip@datDis["Key"]|>, parts], parts];
	debugBids[{parts, partsO, hasKey, typeSuf}];

	If[hasKey,
		(*-----*)AddToLog[{"Multiple analysis - starting:", datDis["Key"]}, 2, True];
	];

	(*file name functions*)
	fileName = If[hasKey,
		GenerateBidsFileName[fol, <|parts, typeSuf -> StringStrip@#[[1]], "type" -> StringStrip@#[[2]], "suf" -> #[[3;;]]|>]&,
		GenerateBidsFileName[fol, <|parts, "type" -> First[#], "suf" -> Rest[#]|>]&
	];
	fileNameO = FileNameJoin[{GenerateBidsFolderName[folOut, #], GenerateBidsName[#]}]&;

	(*checkfiles image and xls*)
	partsO["suf"] = {"xls"}; 
	checkFileX = fileNameO[partsO];
	partsO["suf"] = {"img"}; 
	checkFileI = fileNameO[partsO];

	anaSeg = datDis["Segmentation", "Type"];
	(*the labels setting is a number, a word or both, anything else is a config error*)
	labCfg = BidsAnalysisLabels[datDis["Segmentation", "Labels"]];
	If[labCfg === $Failed,
		(*-----*)AddToLog[{"Segmentation Labels in the analysis config should be a number, a word or both, for example 799, \"Body\" or [799, \"Body\"], it is:", 
			ToString[datDis["Segmentation", "Labels"], InputForm], ". Skipping the analysis."}, 3, True];
		MakeCheckFile[checkFileX, Sort@Join[{"Check" -> "error"}, Normal@datDis]];
		Return[Null, Block]
	];
	{n, what} = labCfg;
	debugBids[{anaSeg, n, what}];
	(*-----*)AddToLog[{"Segmentation file used for analysis is:", StringRiffle[StringStrip/@anaSeg, "_"]}, 3];

	(*Perform the segmentation analysis, what are the label names and volumes*)
	debugBids[{"Segmentation to xls analysis", parts}];
	segFile = fileName[anaSeg]<>".nii";

	Which[
		(*segmentation is already done*)
		CheckFile[checkFileX, "done", verCheck],
		(*-----*)AddToLog[{"Skipping: the segmentation to xls analysis is already done "}, 4],

		(*no segmentation file exists*)
		!NiiFileExistQ[segFile],
		(*-----*)AddToLog[{"The segmentation file does not exist: ", segFile}, 4],

		(*segmentation file existes perform the analysis*)
		True,
		(*-----*)AddToLog[{"Importing and processing the needed segmentation"}, 4];
		{seg, vox} = ImportNii[segFile];
		seg = If[ArrayDepth[seg]===4, seg[[All, 1]], seg];
		segDim = Dimensions[seg];
		labAll = GetSegmentationLabels[seg];
		(*without a number all labels of the segmentation are used*)
		nMax = If[n === Automatic, Max[Append[labAll, 0]], n];
		{seg, musNr} = SelectSegmentations[seg, Range[nMax], False];

		(*no label in the range of the config, for example Legs labels 1-90 for a Body segmentation (101-910)*)
		If[musNr === {},
			(*-----*)AddToLog[{"No labels from 1 to", nMax, "found in the segmentation:",
				If[labAll === {}, "it is empty",
					"it has "<>ToString[Length[labAll]]<>" labels from "<>ToString[Min[labAll]]<>" to "<>ToString[Max[labAll]]
				]<>". Check Segmentation Labels in the analysis config and the Segment Location. Skipping the analysis."}, 3, True];
			MakeCheckFile[checkFileX, Sort@Join[{"Check" -> "error"}, Normal@datDis]];
			Return[Null, Block]
		];

		debugBids["Type segmentations"];
		(*segmentations that replace the default for some data types*)
		typeSeg = Lookup[datDis["Segmentation"], "TypeSegmentation", <||>];
		segFileO = fileName[#] <> ".nii" & /@ typeSeg;
		segMiss = Select[segFileO, !NiiFileExistQ[#] &];
		If[segMiss =!= <||>,
			(*-----*)AddToLog[{"The TypeSegmentation files do not exist:", StringRiffle[Values[segMiss], ", "], ". Skipping the analysis."}, 3, True];
			MakeCheckFile[checkFileX, Sort@Join[{"Check" -> "error"}, Normal@datDis]];
			Return[Null, Block]
		];
		(*-----*)AddToLog[{"Importing the TypeSegmentation files"}, 4];
		segO = BidsImportSeg[#, musNr, segDim] & /@ segFileO;
		If[MemberQ[segO, $Failed],
			(*-----*)AddToLog[{"The TypeSegmentation files for", StringRiffle[Keys@Select[segO, # === $Failed &], ", "], "do not have the grid of the default segmentation", ToString[segDim], ". Skipping the analysis."}, 3, True];
			MakeCheckFile[checkFileX, Sort@Join[{"Check" -> "error"}, Normal@datDis]];
			Return[Null, Block]
		];
		KeyValueMap[
			If[#2[[2]] =!= {},
				(*-----*)AddToLog[{"The segmentation for", #1, "does not have the labels", StringRiffle[ToString /@ #2[[2]], ", "], ", the values will be zero."}, 3, True]
			]&
		, segO];

		debugBids["Use Filter"];
		filt = ConfigLookup[datDis, "Analysis", "UseFilter"];
		If[filt =!= "None",
			(*-----*)AddToLog[{"Using filtering for parameter calculation"}, 4];
			{ran, filt} = filt;
			{filt, vox} = ImportNii[fileName[filt]<>".nii"];
			filtMask = Mask[filt, ran];
			segF = MaskSegmentation[seg, filtMask];
			filt = True;
			,
			filtMask = 1;
			filt = False;
		];

		debugBids["Volumes and cross"];
		(*-----*)AddToLog[{"Calculating the volume of the segmentation"}, 4];
		vol = SegmentationVolume[seg, vox];
		cross = SegmentationCrossSection[seg, vox];

		If[filt, 
			volF = SegmentationVolume[segF, vox];
			crossF = SegmentationCrossSection[segF, vox];
		];

		debugBids["names and labels"];

		(*switch to the correct segmentation label*)
		Which[
			what === "Legs", 
				(*-----*)AddToLog[{"Using the Legs for muscle labeling"}, 4];
				whatName = GetAssetLocation["MuscleLegLabels"];
				whatNr = GetAssetLocation["MusclesLegAllLabels"],
			what === "Arm", 
				(*-----*)AddToLog[{"Using the Arms for muscle labeling"}, 4];
				whatName = GetAssetLocation["MuscleArmLabels"];
				whatNr = GetAssetLocation["MusclesAllLabels"],
			what === "Body", 
				(*-----*)AddToLog[{"Using the Body for muscle labeling"}, 4];
				whatName = GetAssetLocation["MuscleLabels"];
				whatNr = GetAssetLocation["MusclesAllLabels"],
			what === "Shoulder", 
				(*-----*)AddToLog[{"Using the Shoulder for muscle labeling"}, 4];
				whatName = GetAssetLocation["MuscleShoulderLabels"];
				whatNr = GetAssetLocation["MusclesAllLabels"],
			FileExistsQ[what], 
				(*-----*)AddToLog[{"Using custom file for muscle labeling"}, 4];
				whatName = what;
				whatNr = None,
			True, 
				(*-----*)AddToLog[{"Unknown Label type: ", what}, 4];
				whatName = whatNr = None
		];

		(*convert segmentation nr to names, sides and label numbers*)
		If[whatName =!= None,
			musName = MuscleLabelToName[musNr, whatName];
			{musName, sideName} = Transpose[(str = StringSplit[#, "_"];
				If[Last[str] == "Left" || Last[str] == "Right",
					{StringRiffle[Most@str, "_"], Last@str},
					{StringRiffle[str, "_"], "Both"}
				]) & /@ musName];
			sideNr = sideName /. Thread[{"Left", "Right", "Both"} -> {1, 2, 3}];
			,
			musName = ConstantArray["Unknown", Length[musNr]];
			sideName = ConstantArray["Unknown", Length[musNr]];
			sideNr = ConstantArray[0, Length[musNr]];
		];
		If[whatNr =!= None,
			musNr = MuscleNameToLabel[musName, whatNr];
		];

		(*summarize the data labels for export later*)
		dataLabs = Join[
			Thread[{"subject", "subjectID", "session", "sessionID"} -> Transpose@ConstantArray[{
				parts["sub"], ToExpression[Last@StringCases[parts["sub"], NumberString]],
				parts["ses"], ToExpression[Last@StringCases[parts["ses"], NumberString]]}, Length[musNr]
			]],
			If[hasKey, {"type" -> ConstantArray[datDis["Key"], Length[musNr]]}, {}],
			{
				"muscle" -> musName, "muscleID" -> musNr, "side" -> sideName, "sideID" -> sideNr, 
				"volume" -> vol, "cross" -> cross
			},
			If[filt, {"volume_F" -> volF, "cross_F" -> crossF}, {}]
		];

		(*get the labels for analysis, see if and which need to be done using tract based analysis*)
		debugBids["data analysis"];
		(*labels to be analysed*)
		anaType = datDis["Analysis", "Types"];
		anaType = Flatten[Thread /@ anaType, 1];
		(*-----*)AddToLog[{"Analsysis will be performed for:", StringRiffle[StringStrip/@#, "_"]&/@anaType}, 3];

		(*Figure out the tract analysis*)
		debugBids["tract mask"];
		If[KeyExistsQ[datDis["Analysis"], "TractBased"],
			densLab = {If[hasKey, StringStrip@datDis["Analysis", "TractBased"][[1,1]], Nothing], "dwi", "dti", "trk", "dens"};
			densFile = fileName[densLab]<>".nii";
			trType = Flatten[Thread /@ datDis["Analysis", "TractBased"], 1];
			trType = Select[trType, MemberQ[anaType, #] &];
			(*-----*)AddToLog[{"The types with tract weighting will be:", StringRiffle[StringStrip/@#, "_"]&/@trType}, 3];
			(*-----*)AddToLog[{"Import tract mask:", StringRiffle[densLab, "_"]}, 4];
			trMask = ImportNii[densFile][[1]];
			trMask = If[tractWeighting, trMask, Unitize@trMask];
			,
			(*-----*)AddToLog[{"No tract bases analysis given"}, 4];
			trMask = 1;
		];

		(*perform the actual data analysis *)
		(*-----*)AddToLog[{"Starting the data analysis:"}, 3, True];
		(*make the correcet masks*)
		{seg, segT} = BidsMaskSeg[seg, maskErosion, trMask];
		(*the same for the segmentations of the data types, with the filter on the uneroded masks*)
		segO = Association@KeyValueMap[
			#1 -> Join[BidsMaskSeg[#2[[1]], maskErosion, trMask], {If[filt, MaskSegmentation[#2[[1]], filtMask], 0]}] &
		, segO];

		(*loop over all datatypes and perform the mask analysis*)
		data = Flatten[Table[
			datfile = fileName[datType]<>".nii";
			debugBids["data file: ",datfile];
			Which[
				(*data does not exist so skip*)
				!NiiFileExistQ[datfile],
				(*-----*)AddToLog[{"The data does not exist: ", StringRiffle[datType, "_"]}, 4];
				Nothing,

				(*data exists so perform the analysis*)
				True,
				debugBids[datType];
				{data, vox} = ImportNii[datfile];
				(*figure out how to handle this type *)
				tract = MemberQ[trType, datType];
				scale = Switch[datType[[-2;;]], {"dix", "fatfr"} | {"t2", "fatfr"}, 100, {"dix", "t2star"}, 1000, _, 1];
				meanType = (datType[[-2;;]] === {"trk", "seed"} || datType[[-2;;]] === {"trk", "dens"});
				(*-----*)AddToLog[{"Processing file "<>If[tract, "with", "without"]<>" tract weighting:", StringRiffle[StringStrip@datType, "_"]}, 4];

				(*the segmentation of this data type or the default*)
				segKey = SelectFirst[Keys[segO], MemberQ[datType, #] &, None];
				{sg, sgT, sgF} = If[segKey === None, {seg, segT, segF}, segO[segKey]];
				(*-----*)AddToLog[{"Segmentation used for", StringRiffle[StringStrip@datType, "_"], "is:",
					StringRiffle[StringStrip /@ If[segKey === None, anaSeg, typeSeg[segKey]], "_"]}, 4];

				(*mask based analysis*)
				label = StringRiffle[datType[[-2;;]], "_"];
				Join[
					Thread[{label, label<>"_IQR"} -> Transpose[scale GetMaskData[data, If[tract, sgT, sg],
						GetMaskOutput -> If[meanType, "MeanSTD", "MedianIQR"],
						GetMaskOnly -> If[meanType, True, False]
					]]]
					,
					If[filtMask===1, {},
						Thread[{label<>"_F", label<>"_F_IQR"} -> Transpose[scale GetMaskData[data, If[tract, sgT, sgF],
							GetMaskOutput -> If[meanType, "MeanSTD", "MedianIQR"],
							GetMaskOnly -> If[meanType, True, False]
						]]]
					]
				]

			]
		, {datType, anaType}], 1];

		(*merge the data and export*)
		partsO["suf"] = {};
		outFile = fileNameO[partsO];
		debugBids[{"exporting", outFile}];
		(*-----*)AddToLog[{"Data will be exported to:", DirectoryName@outFile}, 3, True];
		data = Dataset[Association /@ Transpose[Thread[#] & /@ Join[dataLabs, data]]];
		(*the check file is only made when both files are exported*)
		If[AllTrue[{Export[outFile<>".xlsx", data], Export[outFile<>".wxf", data]}, StringQ],
			MakeCheckFile[checkFileX, Sort@Join[{"Check"->"done"}, Normal@datDis]],
			(*-----*)AddToLog[{"Export failed, no check file made for:", outFile}, 3]
		];

		(*compress the nii files if compression during ExportNii -> False*)
		If[!compress, CompressNiiFiles[DirectoryName[outFile]]];
	];

	(*----------- make the images -------------*)

	(*figure out which images to make based on setting*)
	{quantIm, segIm, tractIm} = Switch[imOut, 
		"All", {True, True, True},
		"Quantitative", {True, False, False},
		"Segmentation", {False, True, False},
		"Tractography", {False, False, True},
		_, {False, False, False}
	];
	debugBids[{"Image Analysis", {quantIm, segIm, tractIm}}];

	(*make images if needed*)
	Which[
		(*segmentation is already done*)
		CheckFile[checkFileI, "done", verCheck],
		(*-----*)AddToLog[{"Skipping: the iamges are already done"}, 4],

		(*No images need to be made*)
		!AnyTrue[{quantIm, segIm, tractIm}, # &],
		(*-----*)AddToLog[{"No images will be made since option is set to None:"}, 3],

		(*making the images*)
		True,
		(*-----*)AddToLog[{"Starting making the images:"}, 3];

		(*get the reference file and figure out image slice posisions*)
		imRef = datDis["Images", "Reference"];
		(*-----*)AddToLog[{"Checking the reference file for 2D images: ", StringRiffle[imRef, "_"]}, 4];
		reffile = fileName[imRef]<>".nii";
		debugBids[{reffile, NiiFileExistQ[reffile]}];

		(*quantIm and SegIm need refffile, check if there to get needed information.*)
		If[!NiiFileExistQ[reffile],
			(*-----*)AddToLog[{"No reference file skipping 2D quant and seg images."}, 4]; 
			quantIm = segIm = False;
			,
			(*get the ref data for the slice posisions and background of seg images*)
			{ref, vox} = ImportNii[reffile];
			crp = FindCrop[ref, CropPadding -> 10];
			refC = ApplyCrop[ref, crp];
			size = Dimensions[refC] vox;

			pos = GetSlicePositions[GaussianFilter[refC, 15], vox, MakeCheckPlot -> False, 
				DropSlices -> {1, 1, 1}, PeakNumber -> {0, 1, 2}];
			debugBids[pos];
			pos[[1]] = Reverse[Range[0., 1., 1/(Ceiling[Divide @@ size[[;; 2]]] + 1)][[2 ;; -2]] size[[1]]];
			If[Length[pos[[3]]]===1, pos[[3]] = {1./3, 2/3} size[[3]]];
			debugBids[pos];

			(*Function to extract slice data for 2D images*)
			sliceData = Block[{slDat},
				slDat = GetSliceData[ApplyCrop[#, crp, {vox, voxi}], pos, voxi];
				{slDat[[1]], {slDat[[3, 1]], slDat[[2, 1]], Reverse[slDat[[3, 2]], 2]}, {}}
			] &;
		];

		(*3D image function needed for segmentation and tractography*)
		make3DImage = With[{gc = #, sc = 0.75}, 
			ImageResize[ImagePad[ImageCrop@Image[Graphics3D[Table[
				Translate[Rotate[(First@gc), -i 90 Degree, {0, 0, 1}, 0.5 Options[gc, PlotRange][[1, 2, All, 2]]], 
					{sc i Options[gc, PlotRange][[1, 2, 1, 2]], 0, 0}], 
				{i, 0, 3}], Background -> Lighter@Gray, ViewPoint -> {0, -2, .1}, ##] & @@ Join[{
					BoxRatios -> {sc 4, 1, 1} Options[gc, BoxRatios][[1, 2]], 
					PlotRange -> {sc 4, 1, 1} Options[gc, PlotRange][[1, 2]], 
					Options@gc
					}], ImageSize -> {Automatic, 2400}, ImageResolution -> 300], 
			{{60, 60}, {60, 60}}, Lighter@Gray], {Automatic, 2000}]
		]&;

		(*2D image function*)
		make2DImage = With[{di = Max[ImageDimensions[#][[2]]&/@#[[2]]]/4},
		ImagePad[ImageAssemble[ImageResize[#, {Automatic, di}] & /@ Join[
			{ImageAssemble[Transpose@{ImagePad[#, -5] & /@ #[[1]]}, Spacings -> 20, Background -> White]},
			ImagePad[#, -5] & /@ #[[2]]], Spacings -> 20, Background -> White, 
			ImageResolution -> 300], 20, White
		]]&;

		addLabel = ImageAssemble[{
			{ImageCrop[#2, ImageDimensions[#2] - {0, 20}, {Left, Bottom}]}, 
			{LegendImage[#1, First[ImageDimensions@#2], #3]}}
		]&;

		(*------------ quantiatavite images ------------*)

		If[!quantIm,
			(*-----*)AddToLog[{"Not making Quant images since setting is False."}, 4],
			(*-----*)AddToLog[{"Start making Quant images:"}, 4]; 

			debugBids["Making quantitative map images:"];
			Table[
				(*get the color and styling function for the image*)
				cols = Rest@im;
				{cFun, ran, lab} = Switch[Length@cols,
					0, {"BlackToWhite", Automatic, None},
					1, {cols[[1]], Automatic, None},
					2, {cols[[1]], cols[[2]], None},
					3, cols];
				clip = If[cFun === "BlackToWhite", Automatic, Black];

				(*get the filenames for import and export*)
				type = First@im;				
				imFile = fileName[type]<>".nii";
				partsO["suf"] = If[hasKey, type[[2;;]], type];

				(*check if the data file exist*)				
				If[!NiiFileExistQ[imFile],
					(*-----*)AddToLog[{"Cant make image for because data does not exist: ", StringRiffle[StringStrip/@type, "_"]}, 5],
					(*-----*)AddToLog[{"Making image for: ", StringRiffle[type, "_"]}, 5];

					(*import data and make image and export*)
					{imDat, voxi} = ImportNii[imFile];
					img = make2DImage@MakeSliceImages[sliceData@imDat, voxi, 
							ColorFunction -> cFun, PlotRange -> ran, ClippingStyle -> clip, ImageSize -> 2400];
					img = If[lab =!= None, addLabel[cFun, img, lab], img];
					Export[fileNameO[partsO]<>".jpg", img, ImageResolution -> 300];
				]
			, {im, datDis["Images", "QuantImages"]}];
		];

		(*------------ segmentation images ------------*)

		(*check if segmentation can be done*)
		debugBids[{segFile, NiiFileExistQ[segFile]}];
		If[!NiiFileExistQ[segFile], 
			(*-----*)AddToLog[{"Segmentation file does not exist."}, 4];
			segIm = False
		];

		If[!segIm,
			(*-----*)AddToLog[{"Not making Segment images since setting is False."}, 4],
			(*-----*)AddToLog[{"Start making segmentation images:"}, 4]; 
			debugBids["Making segment images:"];

			(*import the segmation*)
			{seg, voxi} = ImportNii[segFile];

			(*make the 2D segmentation image*)
			(*-----*)AddToLog[{"Making 2D Segment image"}, 5]; 
			partsO["suf"] = If[hasKey, anaSeg[[2;;4]], anaSeg[[;;3]] ];
			Export[fileNameO[partsO]<>".jpg", 
				make2DImage@MakeSliceImages[sliceData@ref, {sliceData@seg, GetSegmentationLabels[seg]}, vox,
					ColorFunction -> {"BlackToWhite","RomaO"}, PlotRange -> Automatic, ClippingStyle -> Automatic, ImageSize -> 2400]
			, ImageResolution -> 300];

			(*make the 3D segmentation image*)
			(*-----*)AddToLog[{"Making 3D Segment image"}, 5];
			partsO["suf"] = Join[partsO["suf"], {"vol"}];
			segPl = PlotSegmentations[SelectSegmentations[seg, Range[n]], SelectSegmentations[seg, Range[n+1, n+30]], 
				voxi, ContourResolution -> 2 voxi];
			Export[fileNameO[partsO]<>".jpg", make3DImage@segPl, ImageResolution -> 300];

			(*make the grid segmentation image*)
			(*-----*)AddToLog[{"Making 2D Segment grid image"}, 5];
			partsO["suf"] = If[hasKey, anaSeg[[2;;4]], anaSeg[[;;3]] ];
			partsO["suf"] = Join[partsO["suf"], {"grid"}];
			{ref, crp} = AutoCropData[ref];
			{segPl, lab} = SplitSegmentations[ApplyCrop[seg, crp]];
			SeedRandom[12345];
			segPl = MergeSegmentations[segPl, Join[RandomSample[Select[lab, # <= n &]], Select[lab, # > n &]]];
			Export[fileNameO[partsO]<>".jpg", 
				MakeChannelClassGrid[{ref}, segPl,
					Which[Length[ref] > 32, {4, 8}, Length[ref] > 18, {3, 6}, Length[ref] > 8, {2, 4}, True, {1, 3}]
				]
			, ImageResolution -> 300, ImageSize->{Automatic, 2000}];
		];

		(*------------ tractography images ------------*)

		(*check if tract image can be done*)
		imTrk = datDis["Images", "TractImages"];
		trkfile = fileName[imTrk]<>".wxf";
		debugBids[{trkfile, FileExistsQ[trkfile]}];
		If[!FileExistsQ[trkfile], 
			(*-----*)AddToLog[{"Tract file does not exist."}, 4];
			tractIm = False;
		];

		If[!tractIm,
			(*-----*)AddToLog[{"Not making Tract images since setting is False."}, 4],
			(*-----*)AddToLog[{"Start making tractography images:"}, 4]; 

			(*import the tractography, make the image and export*)
			(*-----*)AddToLog[{"Making 3D tract image"}, 5];
			debugBids["Making tract images:"];
			partsO["suf"] = Join[If[hasKey, imTrk[[2;;4]], imTrk[[;;3]]], {"vol"}];
			Export[fileNameO[partsO]<>".jpg", make3DImage@Import@trkfile, ImageResolution -> 300];
		];

		(*finalize image making*)
		(*-----*)AddToLog[{"Finished making the images"}, 3, True];
		MakeCheckFile[checkFileI, Sort@Join[{"Check"->"done"}, Normal@datDis]];

		(*compress the nii files if compression during ExportNii -> False*)
		If[!compress, CompressNiiFiles[DirectoryName[outFile]]];
	]
];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
