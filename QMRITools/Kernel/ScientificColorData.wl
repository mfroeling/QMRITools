(* ::Package:: *)

(* ::Title:: *)
(*QMRITools ScientificColorData*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`ScientificColorData`"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


AddScientificColours::usage = 
"AddScientificColours[dir] adds the scientific colour data (https://zenodo.org/records/8409685) from the specified folder."

ExtractColorData::usage = 
"ExtractColorData[] Extracts the scientific colordata archive (https://zenodo.org/records/8409685)."


(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*AddScientificColours*)


AddScientificColours[dir_]:=Block[{
		grads,gradDiv,gradMulti,cyclic,all,swatches,colName,groupName,getCol,colRange,allDef,
		nameDef,groupDef,newGroupsPattern,newShemeNames,newShemes,scientificColorMapsGroups,gr
	},
	
	(*activated ColorDataDump*)
	ColorData[];
	
	(*skip if either the files are not there or the colors have been difined already*)
	If[Quiet[Head[ColorData["ScientificColorMapsGroups"]]===ColorData]&&DirectoryQ[dir],
		
		(*color functions defined in scientific colour maps v8*)
		grads={"acton","bamako","batlow","batlowK","batlowW","bilbao","buda","davos","devon","glasgow","grayC","hawaii","imola","lajolla","lapaz","lipari","navia","nuuk","oslo","tokyo","turku"};
		gradDiv={"bam","berlin","broc","cork","lisbon","managua","roma","tofino","vanimo","vik"};
		gradMulti={"fes","bukavu","oleron"};
		cyclic={"bamO","brocO","corkO","romaO","vikO"};
		
		all=Join[grads,gradDiv,gradMulti,cyclic];
		
		(*number of swatches for each*)
		swatches={"10","25","50","100"};
		
		(*functions that generate the information needed for the named colorfunctions*)
		colName=Switch[#2,
			1, {Capitalize[#1],#1,{}},
			2, {Capitalize[#1]<>"Discrete"<>#3,#1<>" discrete "<>#3,{}},
			3, {Capitalize[#1]<>"Categorical",#1<>" categorical",{}}
		]&;
		
		groupName={
			If[#2===1,"Gradients","Indexed"],
			"ScientificColorMaps",
			Which[
				MemberQ[grads,#1],"Sequential",
				MemberQ[gradDiv,#1],"Diverging",
				MemberQ[gradMulti,#1],"MultiSequential",
				MemberQ[cyclic,#1],"Cyclic"
			]<>"Gradients"<>Switch[#2,1,"",2,"Discrete",3,"Categorical"]
		}&;
		
		getCol=Switch[#2,
			1, If[!MemberQ[gradMulti,#1],
				RGBColor/@Import[FileNameJoin[{dir,#1,#1}]<>".txt","Data"],
				Transpose@{Join[Subdivide[0.,0.5,127],Subdivide[0.5,1.,127]], RGBColor/@Import[FileNameJoin[{dir,#1,#1}]<>".txt","Data"]}],
			2, RGBColor[ToExpression@StringSplit[#][[1;;3]]/255]&/@Import[FileNameJoin[{dir,#1,"DiscretePalettes",#1}]<>#3<>".txt","Lines"][[3;;]],
			3, RGBColor/@Import[FileNameJoin[{dir,#1,"CategoricalPalettes",#1}]<>"S.txt","Data"]
		]&;
		
		colRange=Switch[#1,
			1, {0, 1},
			2, {1, ToExpression[#2], 1},
			3, {1, 100, 1} 
		]&;
		
		(*generate all color functions*)
		allDef=Flatten[Table[
			{colName[name,i,j], groupName[name,i,j], 1, colRange[i,j], getCol[name,i,j], ""}
		,{name,all}, {i, If[MemberQ[grads,name], {1,2,3}, {1,2}]}, {j, If[i==2, swatches, {""}]}], 2];
		
		(*modify ColorDataDump*)
		nameDef=allDef[[All,1,1]];
		groupDef=DeleteCases[DeleteCases[DeleteDuplicates[Flatten[allDef[[All,2]]]],"Gradients"],"Indexed"];
		newGroupsPattern=Alternatives@@Join[DataPaclets`ColorDataDump`colorSchemeGroupsPattern/.Alternatives->List,groupDef];
		newShemeNames=Join[DataPaclets`ColorDataDump`colorSchemeNames,nameDef];
		newShemes=Join[DataPaclets`ColorDataDump`colorSchemes,allDef];
		
		DataPaclets`ColorDataDump`colorSchemeGroupsPattern=newGroupsPattern;
		DataPaclets`ColorDataDump`colorSchemeNames=newShemeNames;
		DataPaclets`ColorDataDump`colorSchemes=newShemes;
		
		(*Modify ColorData such it can display the Scientific color map groups*)
		scientificColorMapsGroups={"SequentialGradients","SequentialGradientsDiscrete","SequentialGradientsCategorical","DivergingGradients",
			"DivergingGradientsDiscrete","MultiSequentialGradients","MultiSequentialGradientsDiscrete",
			"CyclicGradients","CyclicGradientsDiscrete"};
			
		Unprotect[ColorData];
		ColorData["ScientificColorMapsGroups"]=scientificColorMapsGroups;
		(gr=#;ColorData[gr]=Sort[Pick[DataPaclets`ColorDataDump`colorSchemeNames,MemberQ[#,gr]&/@DataPaclets`ColorDataDump`colorSchemes[[All,2]],True]])&/@Prepend[scientificColorMapsGroups,"ScientificColorMaps"];
		Protect[ColorData];
	]
]


ExtractColorData[dir_] := Block[{file},
	file = dir<>".zip";
	If[! DirectoryQ[FileNameJoin[{DirectoryName[file], "ColorData"}]],
		If[FileExistsQ[file], Quiet@ExtractArchive[file, DirectoryName[file]],
		Print["DemoData archive does not exist"]]
	];
]

(* ::Section:: *)
(*End Package*)


End[](* End Private Context *)

EndPackage[]
