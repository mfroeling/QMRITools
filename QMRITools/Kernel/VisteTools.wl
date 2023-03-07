(* ::Package:: *)

(* ::Title:: *)
(*QMRITools VisteTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`VisteTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`VisteTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DatWrite::usage = 
"DatWrite[file, data] exports data to *.dat file as binary data using Real32 format."

DTItoolExpInd::usage = 
"DTItoolExpInd[data, file] exports a 3D array data to the file filename DTItool format (*.dat) using DatWrite.
DTItoolExpInd[data, file ,folder] exports data to given file and folder.
DTItoolExpInd[data, file ,folder, add] exports data to given file and folder and adds -add to the filename.";

DTItoolExpTens::usage =
"DTItoolExpTens[tensor] exports a diffustion tensor array to the DTItool format (*.dat).
DTItoolExpTens[tensor, add] exports tensor and adds - add to the filenames.
DTItoolExpTens[tensor, add, folder] exports tensor to the given folder and adds - add to the filenames."

DTItoolExpFile::usage = 
"DTItoolExpFile[file, background, add, voxsize] exports a *.dti text file."

DTItoolExp::usage =
"DTItoolExp[tensor, voxsize] exports tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} and uses XX.dat as background and generates corresponding *dti files.
DTItoolExp[tensor, voxsize, folder] exports tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and uses XX.dat as background and generates corresponding *dti files. 
DTItoolExp[tensor, voxsize, folder, add] exports tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and uses XX.dat as background and generates corresponding *dti files adds - add to the filenames.
DTItoolExp[back, tensor, voxsize] exports background to back.dat and tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} and generates corresponding *dti files. 
DTItoolExp[back, tensor, voxsize, folder] exports background to back.dat and tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and generates corresponding *dti files.
DTItoolExp[back, tensor, voxsize, folder, add] exports background to back.dat and tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and generates corresponding *dti files and adds - add to the filenames."

ExportVol::usage =
"ExportVol[filename, data, voxsize] exports a .vol and .raw file which can be loaded in DTItool 3.0."


ImportDTI::usage = 
"ImportDTI[folder] imports xx.dat, yy.dat, zz.dat, xy.dat, xz.dat and yz.dat from the given folder.
ImportDTI[folder, add] imports xx-add.dat, yy-add.dat, zz-add.dat, xy-add.dat, xz-add.dat and yz-add.dat from the given folder.
ImportDTI[{file1, file2, ..}] imports the given  *.dat files."

DatRead::usage =
"DatRead[file] imports data from file (dtitool *.dat format) as binary data using Real32 format."


ImportVol::usage = 
"ImportVol[] promts for a vol file to open.
ImportVol[\"file\"] inpormts the file.
the function returns data and voxsize."


LoadFiberTracts::usage = 
"LoadFiberTracts[] promts for a .fbs to open.
LoadFiberTracts[\"file\"] imports the file."


(* ::Subsection::Closed:: *)
(*Options*)


BinaryType::usage = 
"BinaryType is an option for ExportVol and must be \"Integer16\" for an integer array and \"Real32\" for a Double array."


(* ::Subsection::Closed:: *)
(*Error Messages*)


DTItoolExpInd::dim = "data is a `1` dimensional array and must be 3 dimensional {slices, x, y}."

DTItoolExpTens::dim = "data is a `1` dimensional array and must be 4 dimensional {6 tensor elements, slices, x, y}."

DatWrite::dim = "Dat write is used for exporting data. The data is a `1` dimensional array and must be 3 dimensional {slices, x, y}."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*DTItool Export*)


(* ::Subsubsection::Closed:: *)
(*DatWrite*)


SyntaxInformation[DatWrite] = {"ArgumentsPattern" -> {_, _}};

DatWrite[file_?StringQ,data_?ArrayQ]:=
Module[{strm},
	If[!ArrayQ[data,3],
		Message[DatWrite::dim,ArrayDepth[data]],
		strm=OpenWrite[file,BinaryFormat->True];
		BinaryWrite[file,Reverse[Dimensions[data]],"UnsignedInteger16",ByteOrdering->-1];
		BinaryWrite[file,N[Flatten[data]],"Real32",ByteOrdering->-1];
		Close[strm];]
		]


(* ::Subsubsection::Closed:: *)
(*DTItoolExpInd*)


SyntaxInformation[DTItoolExpInd] = {"ArgumentsPattern" -> {_, _, _., _.}};

DTItoolExpInd[data_?ArrayQ,filename_String,folder_String:"",add_String:""]:=
Module[{folderi,addi,file},
	If[!ArrayQ[data,3],
		Message[DTItoolExpInd::dim,ArrayDepth[data]],
		addi=If[add!="","-"<>add,""];
		folderi=If[folder!="",folder<>$PathnameSeparator,""];
		file=folder<>filename<>add<>".dat";
		DatWrite[file,N[data]];
		Print["Exported: ",file]]
		];


(* ::Subsubsection::Closed:: *)
(*DTItoolExpTens*)


SyntaxInformation[DTItoolExpTens] = {"ArgumentsPattern" -> {_, _., _.}};

DTItoolExpTens[data_?ArrayQ,add_String:"",folder_String:""]:=
Module[{dir,addi,folderi,file,files},
	If[!ArrayQ[data,4] || Length[data]!=6,
		Message[DTItoolExpTens::dim,ArrayDepth[data]],
		dir={"XX","YY","ZZ","XY","XZ","YZ"};
		addi=If[add!="","-"<>add,""];
		folderi=If[folder!="",folder<>$PathnameSeparator,""];
		files=MapThread[(file=folder<>#1<>add<>".dat";DatWrite[file,N[#2]];file)&,{dir,data}];
		Print["Exported: ",files];
		]];


(* ::Subsubsection::Closed:: *)
(*DTItoolExp*)


SyntaxInformation[DTItoolExp] = {"ArgumentsPattern" -> {_, _, _., _., _.}};

DTItoolExp[tens_?ArrayQ,vox_?ListQ,folder_String:"",add_String:""]:=DTItoolExpi[{},tens,vox,folder,add];

DTItoolExp[back:{_?ArrayQ,_?StringQ},tens_?ArrayQ,vox_?ListQ,folder_String:"",add_String:""]:=DTItoolExpi[{back},tens,vox,folder,add];

DTItoolExp[back:{{_?ArrayQ,_?StringQ}..},tens_?ArrayQ,vox_?ListQ,folder_String:"",add_String:""]:=DTItoolExpi[back,tens,vox,folder,add];


(* ::Subsubsection::Closed:: *)
(*DTItoolExpi*)


DTItoolExpi[back_,tens_,vox_,folder_,add_]:=
Module[{addi,folderi,file,files,fileb,filesb,print},
	addi=If[add!="","-"<>add,""];
	folderi=If[folder!="",folder<>$PathnameSeparator,""];
	
	print=Reap[
		If[back=={},
			fileb=folderi<>"tens"<>addi;
			DTItoolExpFile[fileb<>".dti","XX",addi,vox];
			Sow["No background given, using XX as background.\nExported: "<>fileb<>".dti\n"];
			,
			filesb=Map[(
				fileb=folderi<>#[[2]]<>addi;
				DatWrite[fileb<>".dat",N[#[[1]]]];
				DTItoolExpFile[fileb<>".dti",#[[2]],addi,vox];
				fileb
				)&,back];
			Sow[StringJoin[
				{"Exported: ",Riffle[(#<>".dat"&/@filesb),", "],"as background.\nExported: ",
				Riffle[(#<>".dti"&/@filesb),", "],"\n"
				}]];
			];
		files=MapThread[(
			file=folderi<>#1<>addi<>".dat";
			DatWrite[file,N[#2]];file)&,{{"XX","YY","ZZ","XY","XZ","YZ"},tens}];
			Sow[StringJoin[{"Exported: ",Riffle[files,", "]}]];
	][[2,1]];
	StringJoin[print]
]


(* ::Subsubsection::Closed:: *)
(*DTItoolExpFile*)


SyntaxInformation[DTItoolExpFile] = {"ArgumentsPattern" -> {_, _, _, _}};

DTItoolExpFile[file_String,back_String,add_String,vox_List]:= Export[file,
"/* DTI BMT format */
T
float
XX XX"<>add<>".dat
YY YY"<>add<>".dat
ZZ ZZ"<>add<>".dat
XY XY"<>add<>".dat
XZ XZ"<>add<>".dat
YZ YZ"<>add<>".dat
I "<>back<>add<>".dat
"<>ToString[vox[[2]]]<>" "<>ToString[vox[[3]]]<>" "<>ToString[vox[[1]]]
		,"Text"]


(* ::Subsection:: *)
(*DTItool Import*)


(* ::Subsubsection::Closed:: *)
(*ImpotDTI*)


SyntaxInformation[ImportDTI] = {"ArgumentsPattern" -> {_,_.}};

ImportDTI[folder_String,add_String:""]:=DatRead[folder<>$PathnameSeparator<>#<>If[add=="","","-"<>add]<>".dat"]&/@{"xx","yy","zz","xy","xz","yz"};

ImportDTI[files:{_?StringQ..}]:=Module[{filesc},
filesc=If[StringLength[#]<4,#<>".dat",If[StringTake[#,-4]==".dat",#,#<>".dat"]]&/@files;
DatRead[#]&/@filesc];


(* ::Subsubsection::Closed:: *)
(*DatRead*)


SyntaxInformation[DatRead] = {"ArgumentsPattern" -> {_}};

DatRead[file_?StringQ]:=Module[{dims,res,strm},
	strm=OpenRead[file,BinaryFormat->True];
	dims=BinaryReadList[strm,"UnsignedInteger16",3,ByteOrdering->-1];
	res=BinaryReadList[strm,"Real32",ByteOrdering->-1];
	Close[strm];
	Return[Reverse[First[Fold[Partition,res,dims]],1]]
	];


(* ::Subsection:: *)
(*Export Volume files *)


(* ::Subsubsection::Closed:: *)
(*ExportVol*)


Options[ExportVol]={BinaryType->"Integer16"};

SyntaxInformation[ExportVol] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ExportVol[file_String,dat_?ArrayQ,voxsize_List,OptionsPattern[]]:=
Module[{data,type},
	{data,type}=If[Element[Flatten[dat],Integers],
		{dat,
			If[Max[dat] > 32768 || Min[dat] < -32768,
				"Integer32",
				"Integer16"
				]},
		{N[dat],"Real32"}
		];
	ExportVolRaw[data,file,type];
	ExportVolFile[file,Dimensions[data],voxsize,type];
	]


(* ::Subsubsection::Closed:: *)
(*ExportVolRaw*)


ExportVolRaw[data_,file_,bit_]:=
Module[{strm},
	strm=OpenWrite[file<>".raw",BinaryFormat->True];
	BinaryWrite[strm,Flatten[data],bit];
	Close[strm];
	]


(* ::Subsubsection::Closed:: *)
(*ExportVolFile*)


ExportVolFile[file_,dim_,vox_,bit_]:=Export[file<>".vol",
"Data.FileName     = "<>Last[StringSplit[file,$PathnameSeparator]]<>".raw
Data.Type         = raw
Data.Dimensions   = "<>StringReplace[StringTrim[ToString[Reverse[dim]],("{"|"}")...],","->""]<>"
Data.PixelSpacing = "<>StringReplace[StringTrim[ToString[Reverse[vox]],("{"|"}")...],","->""]<>"
Data.NrBits       = "<>StringCases[bit,NumberString]<>"
Data.NrComponents = 1
"
,"Text"]


(* ::Subsection::Closed:: *)
(*Import Volume files*)


Options[ImportVol] = {};

SyntaxInformation[ImportVol] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

ImportVol[fil_String: "", OptionsPattern[]] := Module[{file,hdr,datafile,dim,vox,bits,rbit,strm,dat,res},
  file = If[fil == "", FileSelect["FileOpen", {"*.vol"}, WindowTitle -> "Select the *.vol file to import"], fil];
  If[file == Null || file === $Canceled, Return[]];
  
  hdr = #[[1]] -> #[[2]] & /@ (StringTrim[StringSplit[#, "="]] & /@ 
      Import[file, "Lines"]);
  
  datafile = DirectoryName[file] <> ("Data.FileName" /. hdr);
  dim = ToExpression[
    "{" <> (StringReplace["Data.Dimensions" /. hdr, " " -> ","]) <> 
     "}"];
  vox = ToExpression[
    "{" <> (StringReplace["Data.PixelSpacing" /. hdr, " " -> ","]) <> 
     "}"];
  bits = "Data.NrBits" /. hdr;
  rbit = Switch[bits, "32", "Real32", "16", "Integer16"];
  
  strm = OpenRead[datafile, BinaryFormat -> True];
  res = BinaryReadList[strm, rbit, ByteOrdering -> -1];
  Close[strm];
  
  dat = Reverse[First[Fold[Partition, res, dim]], 1];
  {dat, Reverse[vox]}
  ]


(* ::Subsection::Closed:: *)
(*LoadFiberTracts*)


SyntaxInformation[LoadFiberTracts] = {"ArgumentsPattern" -> {_.}};


LoadFiberTracts[]:=Module[{file},
	file = FileSelect["FileOpen", {"*.fbs"}, WindowTitle -> "Select the *.fbs file to import"];
	If[file == Null || file === $Canceled, 
		Return[file],
		LoadFiberTracts[file]
	]
]

LoadFiberTracts[filename_] :=  Import[filename, {{"VertexData", "LineData"}}]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
