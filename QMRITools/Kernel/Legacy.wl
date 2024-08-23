(* ::Package:: *)

(* ::Title:: *)
(*QMRITools VisteTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`Legacy`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`Legacy`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
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

FiberDensityMap::usage =
"FiberDensityMap[fiberPoins, dim, vox] generates a fiber density map for the fiberPoins which are imported by LoadFiberTracts. 
The dimensions dim should be the dimensions of the tracked datasets van vox its volxel size."

FiberLengths::usage =
"FiberLengths[fpoints, flines] calculates the fiber lenght using the output from LoadFiberTacts.
FiberLengths[{fpoints, flines}] calculates the fiber lenght using the output from LoadFiberTacts."


ReadVoxSize::usage = 
"ReadVoxSize[filename] imports the voxelsize from a .dcm file. filename must be a string. 
Imports the pixel and slice spacing from the dicom header. Output is a list containg the voxels size {slice thickness, x, y}."

ReadDicom::usage = 
"ReadDicom[folder] imports all dicom files from the given folder.
ReadDicom[{file1, file2,...}] imports all the given filenames.
ReadDicom[folder, {file1, file2,...}] imports all the given filenames from the given folder.
ReadDicom[folder, partsize] imports all dicom files from the given folder and partions them in given partsize.
ReadDicom[{file1, file2, ...}, partsize] imports all the given filenames and partions them in given partsize.
ReadDicom[folder, {file1, file2, ...}, partsize] imports all the given filenames from the given folder and partions them in given partsize."

ReadDicomDiff::usage = 
"ReadDicomDiff[folder, part] imports all dicom files from the given folder and the corresponding diffusion parameters. 

part is the number of diffusion images per slice including the unweighted images."

ReadDicomDir::usage = 
"ReadDicomDir[file] reads the image data from a dicom directory.";

ReadDicomDirDiff::usage = 
"ReadDicomDirDiff[file] reads the image data and relevant diffuison parameters from a dicom directory."

ReadGradients::usage = 
"ReadGradients[folder, nr] imports the diffusion gradient directions from the dicom header of the first nr of files in de given folder.

folder must be a string, nr must be a int. Uses GradRead."

GradRead::usage = 
"GradRead[filename] imports the diffusion gradient direction from a .dcm file.
filename must be a string.";

ReadBvalue::usage = 
"ReadBvalue[folder,nr] imports the gradient directions from the dicom header of the first nr of files in de given folder.
folder must be a string, nr must be a int. Uses BvalRead."

BvalRead::usage = 
"BvalRead[file] imports the bvalue from a .dcm file. file must be a string."

ShiftPar::usage = 
"ShiftPar[B0file.dcm,DTIfile.dcm] imports the parameters from the dicom headeand and calculates the needed values to preform B0 field map correction.
Needs a B0 dicom file and a diffusion dicom file."

ReadBrukerDiff::usage = 
"ReadBrukerDiff[\"\"] imports the bruker diffusion data selected by the input dialog.
ReadBrukerDiff[\"file\"] imports the bruker diffusion data from \"file\", file must be location of 2dseq."


AlignRespLog::usage =
"AlignRespLog[physLog, respirect, scanTime] aligns respirect and physlog data. physLog is output from ImportPhyslog.
resirect is the first output from ImportRespirect.";

ImportPhyslog::usage =
"ImportPhyslog[] imports all physlog files from the folder selcted.
ImportPhyslog[\"forder\"] imports all physlog files from \"folder\" selcted."

ImportRespirect::usage = 
"ImportRespirect[] impors all the respirect log files from the folder selcted.
ImportRespirect[\"folder\"] impors all the respirect log files from the \"folder\" selcted."

PlotPhyslog::usage =
"PlotPhyslog[{time, resp}, {start, stop}] plots the physlog from ImportPhyslog.
PlotPhyslog[{time, resp}, {start, stop}, scanTime] plots the physlog from ImportPhyslog."

PlotRespiract::usage =
"PlotRespiract[data, dataP, scantimes] plots the respirect data to correct peaks. data and dataP are the first outputs of ImportResirect. scantimes is the output from AlignRespLog. 
PlotRespiract[data, dataP, scantimes, steps]."


ROIMask::usage = 
"ROIMask[maskdim, {name->{{{x,y},slice}..}..}] crates mask from coordinates x and y at slice. 
maskdim is the dimensions of the output {zout,xout,yout}."


SetupDataStructure::usage = 
"SetupDataStructure[dcmFolder] makes nii folders and generates nii files for a directory of dmc data where the data is structured per subject."


(* ::Subsection::Closed:: *)
(*Options*)


SeedDensity::usage = 
"SeedDensity is an option for FiberDensityMap. The seedpoint spacing in mm."

BinaryType::usage = 
"BinaryType is an option for ExportVol and must be \"Integer16\" for an integer array and \"Real32\" for a Double array."


RotateGradient::usage = 
"RotateGradient is an option for ReadDicomDirDiff. If False it will also output the gradient direction as stored in the dicom header."

ScaleCorrect::usage = 
"ScaleCorrect is an option for ReadDicom, ReadDicomDiff, ReadDicomDir and ReadDicomDirDiff. \
The dicom image values are corrected for rescale slope, scale slope and rescale intercept."

BmatrixOut::usage = 
"BmatrixOut is a option for ImportBrukerData if True the bmatrix is given, if false the gradients and bvec are given."

ConvertDcm::usage = 
"ConvertDcm is an option for GradRead."


OutputMethod::usage = "OutputMethod can be \"val\" or \"plot\"."

SampleStep::usage= "SampleStep is an option for AlignRespiract."


(* ::Subsection::Closed:: *)
(*Error Messages*)


DTItoolExpInd::dim = "data is a `1` dimensional array and must be 3 dimensional {slices, x, y}."

DTItoolExpTens::dim = "data is a `1` dimensional array and must be 4 dimensional {6 tensor elements, slices, x, y}."

DatWrite::dim = "Dat write is used for exporting data. The data is a `1` dimensional array and must be 3 dimensional {slices, x, y}."


ReadDicom::unknown = "Unknown filename: `1`."

ReadDicom::imp = "Warning: Some files were not imported."

ReadDicom::part = "Warning: Total number of unpartitioned slices is not a multiple of the partition size."

ReadDicom::fls = "No files to import."

ShiftPar::file = "`1` does not exist."

ReadVoxSize::file = "`1` does not exist."

ReadBrukerDiff::seq = "File 2dseq not found at: `1`."

ReadBrukerDiff::proc = "File d3proc not found at: `1`."

ReadBrukerDiff::meth = "File meth not found at: `1`."

ROIMask::war = "there are more slices in the roi set than in the given dimensions."


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


(* ::Subsection::Closed:: *)
(*FiberDensityMap*)


Options[FiberDensityMap] = {SeedDensity -> Automatic};

SyntaxInformation[FiberDensityMap] = {"ArgumentsPattern" -> {_, _, _,OptionsPattern[]}};

FiberDensityMap[fibers_, dim_, vox_, OptionsPattern[]] := 
 Module[{pixindex, density, dens,densi},
  pixindex = GetFiberCoor[fibers, vox];
  pixindex = Transpose[MapThread[Clip[#1, {1, #2}] &, {Transpose[pixindex], dim}]];
  density = CountVoxels[ConstantArray[0, dim], pixindex];
  densi = OptionValue[SeedDensity];
  (*Print[{(Times @@ vox)/0.75,Median[DeleteCases[Flatten[density], 0]]}];*)
  dens = If[NumberQ[densi],
    Times @@ (vox/densi),
    Median[Cases[Flatten[density], Except[0]]]
    ];
  Clip[NormalizeDens[density, dens], {0., 10.}]
  ]

GetFiberCoor = Compile[{{fibcor, _Real, 1}, {vox, _Real, 1}},
   Round[Reverse[fibcor + vox]/vox],
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", 
   Parallelization -> True];

CountVoxels = Compile[{{const, _Integer, 3}, {pix, _Integer, 2}}, Block[{out = const},
    (out[[#[[1]], #[[2]], #[[3]]]] += 1) & /@ pix;
    out
    ]];

NormalizeDens = Compile[{{dens, _Integer, 3}, {n, _Real, 0}}, dens/n];


(* ::Subsection::Closed:: *)
(*FiberLengths*)


SyntaxInformation[FiberLengths] = {"ArgumentsPattern" -> {_, _.}};

FiberLengths[fpoints_, flines_] := FiberLengths[{fpoints, flines}]
FiberLengths[{fpoints_, flines_}] := Module[{len, mpos},
   len = (Length /@ flines) - 1;
   mpos = First@First@Position[len, Max[len]];
   len Mean[EuclideanDistance @@@ Partition[fpoints[[flines[[mpos]]]], 2, 1]]
];


(* ::Subsection::Closed:: *)
(*ReadVoxSize*)


SyntaxInformation[ReadVoxSize] = {"ArgumentsPattern" -> {_}};

ReadVoxSize[fil_String] := Module[{meta, slice, pixel,file},
file = If[fil == "", FileSelect["FileSave",WindowTitle->"Select the dicom file to read the vox size from"], fil];
If[FileExistsQ[file], meta = Import[file, "MetaInformation"];
slice = "SlicesSpacing" /. meta;
pixel = "PixelSpacing" /. meta;
If[StringQ[pixel], 
pixel = "PixelSpacing" /. ("(0028,9110)" /. First["(5200,9230)" /. meta])];
Flatten[{slice, pixel}] 
,
Message[ReadVoxSize::file, file]]]


(* ::Subsection::Closed:: *)
(*ReadGradients*)


SyntaxInformation[ReadGradients] = {"ArgumentsPattern" -> {_,_}};

ReadGradients[fol_String,part_Integer:1]:=
Module[{files,folder},
folder = If[fol == "", FileSelect["Directory",WindowTitle->"Select the directory containing the *.dcm files"], fol];
If[folder == Null, Return[]];
If[!DirectoryQ[folder],Return[]];
files=FileNames["*.dcm",folder];
-Round[If[part==1,
GradRead[files[[1]]],
GradRead/@files[[;;part]]
] . Inverse[(GradRotMat[files[[1]]])],.00001]
]


(* ::Subsection::Closed:: *)
(*GradRead*)


Options[GradRead] = {ConvertDcm->True}

SyntaxInformation[GradRead] = {"ArgumentsPattern" -> {_,OptionsPattern[]}};

GradRead[file_String,OptionsPattern[]] := 
Module[{impfile, grad}, 
impfile = If[StringTake[file, -4] == ".dcm"||!OptionValue[ConvertDcm], file, file <> ".dcm"];
grad = ({"(2005,10B0)", "(2005,10B1)", "(2005,10B2)"} /. 
Import[impfile, "MetaInformation"]);
If[StringQ[grad[[1]]],
Chop[(ImportString[#, "Real32"] // First) & /@ ({"(2005,10B0)", 
"(2005,10B1)", "(2005,10B2)"} /. 
Import[impfile, "MetaInformation"]), 10^-9],
grad]]


(* ::Subsection:: *)
(*ReadBvalue*)


(* ::Subsubsection::Closed:: *)
(*ReadBvalue*)


SyntaxInformation[ReadBvalue] = {"ArgumentsPattern" -> {_,_}};

ReadBvalue[fol_String,part_Integer:1]:=
Module[{files,folder},
folder = If[fol == "", FileSelect["Directory",WindowTitle->"Select the directory containing the *.dcm files"], fol];
If[folder == Null, Return[]];
If[!DirectoryQ[folder],Return[]];
files=FileNames["*.dcm",folder];
If[part==1,
BvalRead[files[[1]]],
BvalRead/@files[[;;part]]
]
]


(* ::Subsubsection::Closed:: *)
(*BvalRead*)


SyntaxInformation[BvalRead] = {"ArgumentsPattern" -> {_}};

BvalRead[file_String]:=
Module[{impfile,out,meta},
impfile=If[StringTake[file,-4]==".dcm",file,file<>".dcm"];
meta=Import[file,"MetaInformation"];
out="(0018,9087)"/.meta;
If[StringQ[out],
If[out=="(0018,9087)",out="(0018,9087)"/.("(2005,140F)"/.meta)];
Return[ImportString[out,"Real64"]//First],
out]
]


(* ::Subsection::Closed:: *)
(*ReadDocomDiff*)


Options[ReadDicomDiff]={ScaleCorrect->False};

SyntaxInformation[ReadDicomDiff] = {"ArgumentsPattern" -> {_,_,OptionsPattern[]}};

ReadDicomDiff[fol_String,part_Integer,OptionsPattern[]]:=
Module[{files,dti,grad,vox,bvec,folder},
folder = If[fol == "", FileSelect["Directory",WindowTitle->"Select the directory containing the *.dcm files"], fol];
If[folder == Null, Return[]];
If[!DirectoryQ[folder],Return[]];
files=FileNames["*.dcm",folder];
dti=iReadDicom[files,part,OptionValue[ScaleCorrect]];
Print["Importing gradients"];
grad=ReadGradients[folder,part];
bvec=ReadBvalue[folder,part];
vox=ReadVoxSize[files[[1]]];
{dti,grad,bvec,vox}
]


(* ::Subsection:: *)
(*ReadDicom*)


(* ::Subsubsection::Closed:: *)
(*ReadDicom*)


Options[ReadDicom]={ScaleCorrect->False};

SyntaxInformation[ReadDicom] = {"ArgumentsPattern" -> {_,_.,_.,OptionsPattern[]}};

(* Input is "folder", patitioning default 0 *)
ReadDicom[fol_String,part_Integer:0,OptionsPattern[]]:=Module[{folder},
folder = If[fol == "", FileSelect["Directory",WindowTitle->"Select the directory containing the *.dcm files"], fol];
If[folder == Null, Return[]];
iReadDicom[FileNames["*.dcm",folder],part,OptionValue[ScaleCorrect]]
]

(* Input is "folder" and {"file1","file2",...}, patitioning default 0 *)
ReadDicom[folder_String,files:{_?StringQ..},part_Integer:0,OptionsPattern[]]:=
iReadDicom[If[StringTake[folder,-2]== $PathnameSeparator,folder,folder<>$PathnameSeparator]<>#&/@files,part,OptionValue[ScaleCorrect]]

(* Input is {"file1","file2",...}, patitioning default 0 *)
ReadDicom[files:{_?StringQ..},part_Integer:0,OptionsPattern[]]:=
iReadDicom[files,part,OptionValue[ScaleCorrect]]


(* ::Subsubsection::Closed:: *)
(*iReadDicom*)


(* Input is list of filnames {"file1","file2",...}, patitioning *)
iReadDicom[files_List, part_Integer, scor_] := 
Module[{imported, output, err = False, filesc, ss, rs, ri, file, nr},
(*check if Filenames have.dcm extention if nog add the extention*)

filesc = If[StringMatchQ[#, ___ ~~ ".dcm"], #, # <> ".dcm"] & /@ files;

If[Length[filesc]==0,Return[Message[ReadDicom::fls]]];

(*Import files,check if file exists,
if not reprot and set error to true*)
imported = DeleteCases[
Monitor[
MapIndexed[(
file = #1;
nr = #2[[1]];
If[FileExistsQ[#],
file = #; Import[#, "Data"][[1]],
err = True;
Message[ReadDicom::unknown, #]]
) &, filesc]
, Column[{Row[{"Importing file: ", file}], 
ProgressIndicator[nr, {0, Length[filesc]}]}]
]
, Null];
If[err, Message[ReadDicom::imp]];

(*see if importe file need partitioning,
if so check is of the number of files i a multiple of the \
partition number,if not report error*)

output = If[part <= 1, imported,
If[Mod[Length[imported], part] != 0,
Message[ReadDicom::part]];
Partition[imported, part]
];

Print["Done importing " <> ToString[Length[filesc]] <> " files!"];

If[scor,
{ss, rs, ri} = 
ReleaseHold[{Hold[
If[StringQ["(2005,100E)"], Unit82Num["(2005,100E)"], 
"(2005,100E)"]], "RescaleSlope", "RescaleIntercept"} /. 
Import[filesc[[1]], "MetaInformation"]];

(*Return[ToPackedArray[(output*rs + ri)/(rs*ss)]];,*)
Return[(output*rs + ri)/(rs*ss)];
,
(*Return[ToPackedArray[Round[output]]];*)
Return[Round[output]];
];
];

Unit82Num[unit8_]:=Module[{bin,sign,exp,sum,indices,out},
bin=Flatten[Reverse[IntegerDigits[ToCharacterCode[unit8],2,8]]];
sign=If[bin[[1]]==1,-1,1];
exp=(FromDigits[bin[[2;;9]],2]-127);
indices=Flatten[Position[bin[[10;;]],1]];
sum=If[exp==-127,Total[2^(-(indices-1))],1+Total[2^-indices]]//N;
out=sign*(2^exp)*sum
]


(* ::Subsection::Closed:: *)
(*ReadDicomDir*)


SyntaxInformation[ReadDicomDir] = {"ArgumentsPattern" -> {_,_.}};

ReadDicomDir[fil_,part_:0]:=
Module[{meta,slice,directions,groups,sliceSpacing,pixelSpacing,data,output,file},
file = If[fil == "", FileSelect["FileOpen",WindowTitle->"Select the enhanced DICOM file"], fil];
If[file == Null, Return[]];
If[!FileExistsQ[file],Return[]];
meta = Import[file,"MetaInformation"];
slice = "(2001,1018)"/.meta;
directions = "(2005,1415)"/.meta;
groups = ("(5200,9230)"/.meta);
sliceSpacing = {"SlicesSpacing"}/.meta;
pixelSpacing = "PixelSpacing"/.("(0028,9110)"/.groups[[1]]);
data = Import[file,"Data"];
output=If[part<=1,data,
If[Mod[Length[data],part]!=0,Message[ReadDicom::part]];
Partition[data,part]];

{output,Flatten[{sliceSpacing,pixelSpacing}]}
]


(* ::Subsection::Closed:: *)
(*ReadDicomDiff*)


Options[ReadDicomDirDiff]={RotateGradient->True}

SyntaxInformation[ReadDicomDirDiff] = {"ArgumentsPattern" -> {_,OptionsPattern[]}};

ReadDicomDirDiff[fil_,OptionsPattern[]]:=
Module[{meta, slice, directions, groups, sliceSpacing, pixelSpacing, 
orientation, gradRotmat, grads, gradsRot, data, bvec,file}, 
file = If[fil == "", FileSelect["FileOpen",WindowTitle->"Select the enhanced DICOM file"], fil];
If[file == Null, Return[]];
If[!FileExistsQ[file],Return[]];
meta = Import[file, "MetaInformation"];
slice = "(2001,1018)" /. meta;
directions = "(2005,1415)" /. meta;

(*groups = ("(5200,9230)" /. meta)[[1 ;; directions]];*)
groups = If[("FrameCount"/slice /. meta) == directions,
("(5200,9230)" /. meta)[[1 ;; directions]],
("(5200,9230)" /. meta)[[1 ;;("FrameCount"/"(2001,1018)") /. meta]]
];

sliceSpacing = {"SlicesSpacing"} /. meta;
pixelSpacing = "PixelSpacing" /. ("(0028,9110)" /. groups[[1]]);
orientation = "ImageOrientation" /. ("(0020,9116)" /. groups[[1]]);
gradRotmat = Transpose[{orientation[[1 ;; 3]], orientation[[4 ;; 6]], Cross[orientation[[1 ;; 3]], orientation[[4 ;; 6]]]}];
grads = If[("(0018,9075)" /. #) == "NONE" || ("(0018,9075)" /. #) == 
"ISOTROPIC", {0, 0, 
0}, ("(0018,9089)" /. ("(0018,9076)" /. #))] & /@ (("(0018,9117)" /. #) & /@ groups);
bvec = "(0018,9087)" /. ("(0018,9117)" /. #[[2]]) & /@ groups;
gradsRot = Round[grads . gradRotmat, .0001];
data = If[("FrameCount"/slice /. meta) == directions,
Partition[Import[file, "Data"], directions],
Partition[Import[file, "Data"],("FrameCount"/"(2001,1018)") /. meta]
(*Drop[#, -1] & /@ Partition[Import[file, "Data"], directions + 1]*)
];
If[OptionValue[RotateGradient],
{data, gradsRot, bvec, Flatten[{sliceSpacing, pixelSpacing}]},
{data, {gradsRot,grads}, bvec, Flatten[{sliceSpacing, pixelSpacing}]}
]
]


(* ::Subsection::Closed:: *)
(*ShiftPar*)


SyntaxInformation[ShiftPar] = {"ArgumentsPattern" -> {_,_}};

ShiftPar[pFile_,dFile_]:=Block[{metaDTI, metaB0, dDir, phaseM, freqN, phaseSteps, echoLength, bw, te, TEB0}, 
If[! FileExistsQ[pFile], Return[Message[ShiftPar::file, pFile]]];
If[! FileExistsQ[dFile], Return[Message[ShiftPar::file, dFile]]];
metaDTI = Import[dFile, "MetaInformation"];
metaB0 = Import[pFile, "MetaInformation"];
If[StringTake[dFile, -4] == StringTake[pFile, -4] == ".dcm",
dDir = "PhaseEncodingDirection" /. metaDTI;
{phaseM, freqN, phaseSteps, echoLength, bw} = 
Join[If[dDir == "COL", {"Columns", "Rows"}, {"Rows", 
"Columns"}], {"PhaseEncodingSteps", "EchoTrainLength", 
"PixelBandwidth"}] /. metaDTI;
,
{phaseM, freqN, phaseSteps, echoLength, bw} = Join[
If[(dDir = ("PhaseEncodingDirection" /. ("(0018,9125)" /. \
("(5200,9229)" /. metaDTI)))) == "COLUMN", {"Columns", 
"Rows"}, {"Rows", "Columns"}] /. 
metaDTI, {"PhaseEncodingSteps"} /. ("(2005,140E)" /. \
("(5200,9229)" /. metaDTI)),
{"EchoTrainLength"} /. ("(0018,9112)" /. ("(5200,9229)" /. 
metaDTI)),
{"PixelBandwidth"} /. ("(0018,9006)" /. ("(5200,9229)" /. metaDTI))
]
];
te = "(2001,1025)" /. metaB0;
TEB0 = ToExpression[StringTake[te, -3] <> "-" <> StringTake[te, 3]] 10^-3;
{(1/(2 Pi TEB0))*(freqN*echoLength/(bw*phaseSteps)), dDir} /. "COLUMN" -> "COL"]


(* ::Subsection::Closed:: *)
(*ReadBrukkerDiff*)


Clear[ReadBrukerDiff]
Options[ReadBrukerDiff] = {BmatrixOut -> True};

SyntaxInformation[
ReadBrukerDiff] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ReadBrukerDiff[file_, OptionsPattern[]] := 
Module[{type, name, title, file1, file2, file3, imstr, destr, mestr, mestrr, binType, flipDim, 
xdim, ydim, zdim, diffdir, sdim, dat, bmat, grad, bvec, vox},

type = {"2dseq"};
name = "Image files";
title = "Select the 2dseq file";

file1 = If[file == "",
SystemDialogInput["FileOpen", {Directory[], {name -> type}}, 
WindowTitle -> title],
file];
If[file1 === $Canceled, Return["Canceled!"]];
file2 = DirectoryName[file1] <> "d3proc";
file3 = DirectoryName[file1, 3] <> "method";

If[! FileExistsQ[file1], 
Return[Message[ReadBrukerDiff::seq, file1]]];
If[! FileExistsQ[file2], 
Return[Message[ReadBrukerDiff::proc, file2]]];
If[! FileExistsQ[file3], 
Return[Message[ReadBrukerDiff::meth, file3]]];

destr = Import[file2, "Text"];
mestr = Import[file3, "Text"];
mestrr = StringReplace[mestr, "\n" -> " "];

xdim = ToExpression[
StringCases[destr, 
"##$IM_SIX=" ~~ (x : DigitCharacter ..) -> x][[1]]];
ydim = ToExpression[
StringCases[destr, 
"##$IM_SIY=" ~~ (x : DigitCharacter ..) -> x][[1]]];
zdim = ToExpression[
StringCases[destr, 
"##$IM_SIZ=" ~~ (x : DigitCharacter ..) -> x][[1]]];

diffdir = 
ToExpression[
StringCases[mestr, 
"##$PVM_DwNDiffExp=" ~~ (x : DigitCharacter ..) -> x][[1]]];
sdim = zdim/diffdir;

binType = 
StringCases[mestrr, 
"##$RECO_wordtype=_" ~~ (x : NumberString) ~~ "BIT_SGN_INT" -> 
x][[1]];
imstr = Import[file1, "Integer" <> binType, ByteOrdering -> -1];

flipDim = Switch[
StringCases[mestrr, 
"##$PVM_SPackArrReadOrient=( 1 ) " ~~ x__ ~~ 
" ##$PVM_SPackArrReadOffset=(" -> x][[1]],
"L_R", False,
"H_F", True,
_, False
];

dat = If[flipDim,
Transpose[Fold[Partition, imstr, {ydim, xdim, sdim}]],
Transpose[Fold[Partition, imstr, {xdim, ydim, sdim}]]
];

bmat = Append[-{1, 1, 1, 2, -2, -2} TensVec[#], 1] & /@ 
Fold[Partition, 
ToExpression[
"{" <> StringReplace[
StringCases[mestrr, 
"##$PVM_DwBMat=( " <> ToString[diffdir] <> ", 3, 3 ) " ~~ 
x__ ~~ " ##$PVM_DwEffBval=( " <> ToString[diffdir] <> 
" )" -> x], {"  " -> ",", " " -> ","}] <> "}"], {3, 3}];

vox = Flatten[{
ToExpression[
StringCases[mestrr, 
"##$PVM_SliceThick=" ~~ x : NumberString -> x]][[1]],
Reverse[
ToExpression[
StringCases[mestrr, 
"##$PVM_SpatResol=( 2 ) " ~~ x : NumberString ~~ " " ~~ 
y : NumberString -> {x, y}]][[1]]]
}];

If[OptionValue[BmatrixOut],
{dat, bmat, vox},
{grad, bvec} = BmatrixInv[bmat[[All, ;; 6]]];
{dat, grad, bvec, vox}
]
]


(* ::Subsection::Closed:: *)
(*AlignRespLog*)


Options[AlignRespLog] = Options[AlignRespLogi] = {OutputMethod -> "val", SampleStep -> 0.005};

SyntaxInformation[AlignRespLog] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

AlignRespLog[physLog_, respirect_, scanTime_, opts:OptionsPattern[]]:=AlignRespLog[physLog, respirect, scanTime, Range[Length[respirect]], opts]

AlignRespLog[physLog_, respirect_, scanTime_, order_, opts:OptionsPattern[]]:= Block[{n},
Switch[OptionValue[OutputMethod],
"val",
Transpose[AlignRespLogi[physLog, respirect[[#]], scanTime[[#]], opts] & /@ order],
"plot",
Manipulate[AlignRespLogi[physLog, respirect[[order[[n]]]], scanTime[[order[[n]]]], opts] ,{{n,1,"dataset"},1,Length[order],1}]
]
];

AlignRespLogi[physLog_, respirect_, scanTime_, OptionsPattern[]] := 
Block[{ptime, pdat, start, stop, rtime, rdat, rstart, rend, pstart, co2dataSel,
pend, pran, samp, rtimei, rint, ptimei, pint, corr, poff, len,corrAll,n,nsel,
startscan, stopscan},

co2dataSel=respirect[[1]];

{rtime, rdat} = co2dataSel[[All, {1, 4}]] // Transpose;
samp = OptionValue[SampleStep];
len = Length[physLog];

corrAll = {corr, rstart, samp, pint, ptime, pdat, stop} = Transpose[
(
n = #;
{ptime, pdat} = physLog[[n, 1]];
{start, stop} = physLog[[n, 2]];

{rstart, rend} = StartEnd[rtime];
{pstart, pend} = StartEnd[ptime];

pran = pend - pstart;

rtimei = Range[rstart, rend, samp];
rint = Interpolation[{rtime, rdat} // Transpose][rtimei];
ptimei = Range[pstart, pend, samp] + rstart;
pint = 
Interpolation[{ptime + rstart, pdat} // Transpose][ptimei];
corr = Abs[ListCorrelate[pint, rint, {-1, 1}, 0]];

{corr, rstart, samp, pint, ptime, pdat, stop}
) & /@ Range[len]];

nsel = First[Flatten[Position[corrAll[[1]], Max[corrAll[[1]]]]]];
{corr, rstart, samp, pint, ptime, pdat, stop} = corrAll[[All, nsel]];
poff = rstart + samp (First@Flatten[Position[corr, Max[corr]]] - Length[pint]);
{startscan, stopscan} = {stop - scanTime, stop} + poff;

Switch[OptionValue[OutputMethod],
"val", {{startscan, stopscan},physLog[[nsel]],respirect},
"plot", ListLinePlot[
{{rtime, 1 - (rdat - Min[rdat])/(Max[rdat] - Min[rdat])} // Transpose,
{ptime + poff, (pdat - Min[pdat])/(Max[pdat] - Min[pdat])} // Transpose
},
AspectRatio -> 0.2, PlotRange -> {Full, Full}, ImageSize -> 1000, 
PlotStyle -> {Red, Black}, 
GridLines -> {{{startscan, 
Directive[{Thickness[.005], Green}]}, {stopscan, 
Directive[{Thickness[.005], Red}]}}, None}]]
]

StartEnd = {First[#], Last[#]} &;


(* ::Subsection::Closed:: *)
(*AlignRespLog*)


SyntaxInformation[ImportPhyslog] = {"ArgumentsPattern" -> {_.}};

ImportPhyslog[] := ImportPhyslog[FileSelect["Directory" ,WindowTitle->"Select directory that contains the physlogs"]]
ImportPhyslog[folder_] := Block[
{files, file, resp, mark, time, start, stop, sel},
files = Sort[FileNames["*.log", {folder}]];
(
{resp, mark} = 
Drop[ToExpression[
"{" <> StringReplace[#, " " .. -> ","] <> "}"] & /@ 
DeleteCases[Import[files[[#]], "Lines"][[6 ;;]], 
"#"], -1][[All, {6, 10}]] // Transpose;
time = Range[0, Length[resp] - 1] (1/500/60.);

start = time[[First[Flatten[Position[mark, 10]]]]];

stop = Flatten[Position[mark, 20]];
stop = time[[If[stop == {}, Length[mark], Last[stop]]]];

sel = Range[1, Length[time], 50];

{{time, resp}[[All, sel]], {start, stop}}

) & /@ Range[Length[files]]
]


(* ::Subsection::Closed:: *)
(*PlotPhyslog*)


SyntaxInformation[PlotPhyslog] = {"ArgumentsPattern" -> {_, _, _.}};

PlotPhyslog[{time_, resp_}, {start_, stop_}, scanTime__: 0] := 
Block[{start2},
start2 = If[scanTime == 0, start, stop - scanTime];
ListLinePlot[Transpose[{time, resp}], GridLines -> {
{{start, Directive[{Dashed, Thickness[0.005], Green}]}, {start2, 
Directive[{Thickness[0.005], Green}]}, {stop, 
Directive[{Thickness[0.005], Red}]}}
, None}, PlotStyle -> Black, AspectRatio -> 0.1,
Axes->False,Frame->{{True,False},{True,False}}, 
ImageSize -> 1000, PlotLabel -> {stop - start, stop - start2}]
]


(* ::Subsection::Closed:: *)
(*ImportRespirect*)


SyntaxInformation[ImportRespirect] = {"ArgumentsPattern" -> {_.}};

ImportRespirect[] := ImportRespirect[FileSelect["Directory",WindowTitle->"Select directory that contains the respiract data"]]
ImportRespirect[folder_] := 
Module[{co2data, co2dataP, events, events2, checks, xsel, timesel, sel, start, end, a, b},

If[folder === Null, 

Return[],

co2dataP = (Import[FileNames[{"bbb_*.txt"}, {folder}][[1]], "Data"])[[2 ;;, 1 ;; 3]];
co2data = (DeleteCases[Import[FileNames[{"raw_*.txt"}, {folder}][[1]],"Data"], {""}])[[3 ;;, {1, 4, 5, 2}]];
events = (Import[FileNames[{"events_*.txt"}, {folder}][[1]], "Data"])[[2 ;;]];

events2 = Delete[events, Position[events, "END sequence"][[All, {1}]]];

DynamicModule[{
x = ConstantArray[False, Length[events2]],
time = 6},

checks = {Checkbox[Dynamic[x[[#]]]]} & /@ Range[Length[events2]];

{xsel, timesel} = DialogInput[
{
TextCell["Choose events: "],
Row[{Thread[{checks, events2[[All, 2]]}] // TableForm}],
TextCell["Enter the respirect experiment duration: "],
InputField[Dynamic[time], Number],
DefaultButton[DialogReturn[{x, time}]]
}, Modal -> True];
];

If[Total[Boole[xsel]] == 0,

Return[Print["Error, more than one event need to be selected"]],

sel = Flatten[Position[xsel, True]];

((
start = #;
end = start + timesel;
a = Select[co2data, end > #[[1]] > start &];
b = Select[co2dataP, end > #[[1]] > start &];
{a, b}
) & /@ events2[[sel, 1]]) 
]
]
]


(* ::Subsection:: *)
(*ImportRespirect*)


(* ::Subsubsection::Closed:: *)
(*ImportRespirect*)


SyntaxInformation[PlotRespiract] = {"ArgumentsPattern" -> {_, _, _.}};

PlotRespiract[dataAll_, scantimes_, steps__: 10] := PlotRespiracti[dataAll, scantimes, steps][[2]]

PlotRespiracti[dataAll_, scantimes_, steps_] := DynamicModule[
{datas, samplesO2, allx, samplesCO2, pos, CO2, O2, pt, pd, peak,data,dataP,
allCO2, allO2,  max, min, scant, r, but, temp,range, allPmouth, samplesx,
rr,  CO2tot, O2tot, dataout, len, shift, n, pointsall, x,cent, span, pmin, pmax, pran},

{data, dataP}=Transpose[dataAll];

range = Range[1, Length[data[[1, All, 1]]], steps];
len = Length[data];
shift = {0, 1.2, 0, 0};

CO2tot = O2tot = ConstantArray[0, len];
dataout = Transpose /@ dataP;

datas = Transpose /@ data[[All, range]];

{min, max} = 
Transpose[Transpose[{Min[#], Max[#]} & /@ #] & /@ datas];
rr = max - min; rr[[All, 1]] = 1;

{allx, allCO2, allO2, allPmouth} = Transpose[# + shift & /@ ((datas - min)/rr)];

{samplesx, samplesCO2, samplesO2} = 
Transpose[# + shift[[1 ;; 3]] & /@ ((dataout - min[[All, 1 ;; 3]])/
rr[[All, 1 ;; 3]])];

scant = scantimes - min[[All, 1]];
pointsall = MapThread[Transpose[{#1, #2}] &, {samplesx, samplesCO2}];
temp = pointsall = MapThread[(r = #2; Select[#1, r[[1]] < #[[1]] < r[[2]] &]) &, {pointsall, scant}];
n = 1; but = False;

{pmin,pmax}={Min[allx],Max[allx]};
pran=pmax-pmin;

DialogInput[{
Manipulate[
cent = Clip[cent,{pmin+span,pmax-span}];
LocatorPane[
Dynamic[pointsall[[n]]],
Dynamic[
n = Clip[n, {1, len}, {1, len}];

pos = 
Position[allx[[n]], First[Nearest[allx[[n]], #, 1]]][[1, 
1]] & /@ pointsall[[n, All, 1]];

x = allx[[n, pos]];
CO2 = pointsall[[n, All, 2]] = allCO2[[n, pos]];
O2 = allO2[[n, pos]];

dataout[[n]] = ({x, CO2 - shift[[2]], O2} rr[[n, 1 ;; 3]]) + min[[n, 1 ;; 3]];

CO2tot[[n]] = Round[Mean[dataout[[n, 2]]], .1];
O2tot[[n]] = Round[Mean[dataout[[n, 3]]], .1];


Graphics[
{
Thickness[.0015],
Darker[Darker[Red]], Line[Transpose[{allx[[n]], allCO2[[n]]}]],
Darker[Darker[Blue]], Line[Transpose[{allx[[n]], allO2[[n]]}]]
,
PointSize[Large],
Red, Point[Transpose[{x, CO2}]],
Blue, Point[Transpose[{x, O2}]]
,
Thickness[.001],
Red, Line[SortBy[Transpose[{x, CO2}], First]],
Blue, Line[SortBy[Transpose[{x, O2}], First]],

Opacity[.5], Thickness[.0025], Green,
Line[{{scant[[n, 1]], -0.1}, {scant[[n, 1]], 2.3}}],
Red,
Line[{{scant[[n, 2]], -0.1}, {scant[[n, 2]], 2.3}}]
},

Background -> White,
PlotRange->{{cent-span,cent+span},{-0.1,2.4}},
ImageSize -> 1200, AspectRatio -> 0.4, Frame -> True, 
FrameTicks -> {True, False}, LabelStyle -> Large,
GridLines->{x,None},
LabelStyle->Black,
PlotLabel -> Column[{
Style["Dataset "<>ToString[n],Bold,Large],
Style[LabFun[x, CO2tot[[n]], O2tot[[n]]], Bold, Medium], 
TableForm[{CO2tot, O2tot}, TableHeadings -> {{"CO2", "O2"}, Range[len]}]}, Alignment -> Center]
]

], LocatorAutoCreate -> True, Appearance -> None]

,{{cent,pran/2,"Center"},Dynamic[pmin+span],Dynamic[pmax-span],.1}
,{{span,pran/2,"Span"},1,pran/2,.1}
]
,
Dynamic[Grid[{
{
Button["   Back   ", If[n > 1, n--]],
Button["   Next   ", If[n < len, n++]],
CancelButton[]

},
{
Button["Peak Detect", (
temp[[n]] = pointsall[[n]];
{pt, pd} = 
Transpose[
Select[Transpose[{allx[[n]], allCO2[[n]]}], 
scant[[n, 1]] < #1[[1]] < scant[[n, 2]] &]];

peak = (Partition[
PeakDetect[pd/Mean[pd], 20/steps, 0, 0.5, Padding -> 10], 
4, 1, 4] /. {{0, 0, 0, 1} -> 1, {_, _, _, _} -> 0});

pos = Position[allx[[n]], #][[1, 1]] & /@ (pt = 
DeleteCases[(peak pt), 0.]);
pointsall[[n]] = Transpose[{pt, allCO2[[n, pos]]}];
)]
,
Button["    Undo    ", pointsall[[n]] = temp[[n]]],
If[n == len || but, but = True; 
DefaultButton[
DialogReturn[{CO2tot, O2tot, Transpose[{data, SortBy[Transpose[#],1]& /@ dataout}]}]]]
}
}]]
}, Modal -> True, WindowFloating -> True]
]


(* ::Subsubsection::Closed:: *)
(*LabFun*)


LabFun[x_, co2_, o2_] := Block[{st, en},
StringJoin@(ToString /@ {
"CO2:  ", co2, "   ",
"O2:  ", o2, "\n",
"start:  ", st = Round[Min[x], .1], "  ",
"end:  ", en = Round[Max[x], .1], "  ",
"span:  ", en - st, "  "
})]


(* ::Subsection::Closed:: *)
(*PlotData3D Old*)

(*
SyntaxInformation[PlotData3D] = {"ArgumentsPattern" -> {_, _.}};

PlotData3D[data_, vox:{_, _, _}:{1, 1, 1}] := Module[{
   tran, depth, mind, maxd, dimd, ratio, dim, setp, slicep, columnp, rowp, line,
  scal, dats, pdat, im3D,colf, linax, lincor, linsag, plax, plcor, plsag, imall, sel, slices3D, plsl3D, plim3D, dat,
  imax, imcor, imsag, size1, size2, size3, size1b, mins,
     maxs
  },
 
 NotebookClose[plotwindow3D];
 ClearTemporaryVariables[];
 
 tran = False;
 depth = ArrayDepth[data];
 mind = Min[data];
 maxd = Max[data];
 dat = If[depth == 4, Reverse@ToByte[data, {mind, maxd}], 
    Reverse@ToByte[data, {mind, maxd}]];
 dimd = Dimensions[dat];
 
   size1b = 1;
   
   pan = Manipulate[
     
     t0 = (
         (*determine scaling and dimensions*)
         dim = {size1, size2, size3} = If[depth == 4, 
            If[trans, 
            	size1b = dimd[[1]]; dimd[[{2, 3, 4}]],
            	size1b = dimd[[2]]; dimd[[{1, 3, 4}]]
            	], 
            dimd];
         ratio = Reverse[vox*dim]/Max[(vox*dim)];
         
         (*correctly clip the slice numbers and mirror slices if needed*)
         If[depth == 4 && trans != tran, {tran, slice, set} = {trans, set, slice}];
         
         set = If[NumericQ[set], Clip[set, {1, size1b}], 1];
         slice = If[NumericQ[slice], Clip[slice, {1, size1}], 1];
         column = If[NumericQ[column], Clip[column, {1, size2}], 1];
         row = If[NumericQ[row], Clip[row, {1, size3}], 1];
         
         {setp, slicep, columnp, rowp} = If[trans && depth == 4, 
         	{size1b - set + 1, slice, column, size3 - row + 1}, 
         	{set, size1 - slice + 1, column, size3 - row + 1}];
         
         (*determine to draw lines and adjust pot scaling for all pannels*)
         {line, scal} = 
          If[show != 4, {False, scale}, {lines, 0.5 scale}];
         ) // AbsoluteTiming // First;
     
     (*create im3D*)
     t1 = (
         (*rescale 3D values for 3D image and select correct dataset*)

         
         dats = 
          If[depth == 4, 
           If[trans, {dat[[setp]]}, {dat[[All, setp]]}], {dat}];
         mins = Min[dats];
         maxs = Max[dats];
         min3D = Clip[min3D, {mins, 1}];
         max3D = Clip[max3D, {0, 2 maxs}];
         
         (*create im3D*)
         pdat = 
          If[show == 6, ToByte[##, {min3D, max3D}], ##] & @@
            dats;
         im3D = 
          Image3D[##, ColorFunction -> col3D] & @@ 
           If[reverse, {Reverse[pdat]}, {pdat}];
         ) // AbsoluteTiming // First;
     
     (*Create dynamic color function, 1000 values, 
     with clipping for min and max vals*)
     t2 = (
         If[show != 6,
         	(*colf = LookUpTable3[cfs, {lstyle, color}, {minclip, maxclip},If[cfs,{mind,maxd,mind,maxd},{mind,maxd,min,max}]]*)
         	colf = LookUpTable[{lstyle, color}, {minclip, maxclip}]         	
         	]
         ) // AbsoluteTiming // First;
     
     (*create the lices for the pannel all view*)
     t3 = (If[show == 4 && lines,
          
          linax = Graphics[{Red, Thickness[.01], 
             Line[{{column, 0}, {column, size2}}], Green, 
             Thickness[.01], Line[{{0, row}, {size3, row}}]}];
          
          lincor = 
           Graphics[{Red, Thickness[.01], 
             Line[{{column, 0}, {column, size3}}], Blue, 
             Thickness[.01], Line[{{0, slice}, {size3, slice}}]}];
          
          linsag = 
           Graphics[{Green, Thickness[.01], 
             Line[{{row, 0}, {row, size3}}], Blue, Thickness[.01], 
             Line[{{0, slice}, {size3, slice}}]}];
          ]
         ) // AbsoluteTiming // First;
     
     (*create the axial coronal and sagital images*)
     t4 = (
         If[MemberQ[If[planez, {1, 4, 5}, {1, 4}], show],
          
          imax = Colorize[Image3DSlices[im3D, {slicep}, 1][[1]], 
            ColorFunction -> colf, ColorFunctionScaling -> cfs];
          
          plax = Show[##, ImageSize -> scal {ratio[[1]], ratio[[2]]}, 
              AspectRatio -> Full] & @@ 
            If[line, {imax, linax}, {imax}];
          ];
         If[MemberQ[If[planey, {2, 4, 5}, {2, 4}], show],
          
          imcor = Colorize[Image3DSlices[im3D, {rowp}, 2][[1]], 
            ColorFunction -> colf, ColorFunctionScaling -> cfs];
          
          plcor = Show[##, ImageSize -> scal {ratio[[1]], ratio[[3]]},
               AspectRatio -> Full] & @@ 
            If[line, {imcor, lincor}, {imcor}];
          ];
         If[MemberQ[If[planex, {3, 4, 5}, {3, 4}], show],
          
          imsag = ImageReflect[
            Colorize[Image3DSlices[im3D, {columnp}, 3][[1]], 
             ColorFunction -> colf, ColorFunctionScaling -> cfs], 
            Left -> Right];
          
          plsag = Show[##, ImageSize -> scal {ratio[[2]], ratio[[3]]},
               AspectRatio -> Full] & @@ 
            If[line, {imsag, linsag}, {imsag}];
          ];
         ) // AbsoluteTiming // First;
     
     (*create the pannel all view with event handles to click and \
select*)
     t5 = (
         If[show == 4,
           imall = Grid[{{
               
               EventHandler[
                plax, {"MouseDown" :> ({column, row} = 
                    Abs[Round[
                    MousePosition["Graphics"]] - {-1, -1}])}]}, {
               
               EventHandler[
                plcor, {"MouseDown" :> ({column, slice} = 
                    Abs[Round[MousePosition["Graphics"]] - {-1, 1}])}],
               
               EventHandler[
                
                plsag, {"MouseDown" :> ({row, 
                    slice} = (Abs[
                    Round[MousePosition["Graphics"]] - {-1, 1}]))}]
               }}, Background -> White, Spacings -> {0, 0}, 
             Frame -> All, 
             FrameStyle -> Directive[{Thickness[6], White}]]];
         ) // AbsoluteTiming // First;
     
     If[MemberQ[views[[All, 1]], vp], vv = {0, 0, 1}];
     
     t6 = (If[show == 5,
           
           sel = DeleteCases[{If[planez, 1], If[planey, 2], 
              If[planex, 3]}, Null];
           
           slices3D = 
            If[sel == {}, {}, {Opacity[{opz, opy, opx}[[#]]], 
                Dynamic[Texture[{imax, imcor, imsag}[[#]]]], Polygon[{
                   {{1, 1, slice}, {size3, 1, slice}, {size3, size2, 
                    slice}, {1, size2, slice}},
                   {{1, row, 1}, {size2, row, 1}, {size2, row, 
                    size1}, {1, row, size1}},
                   {{column, 1, 1}, {column, size3, 1}, {column, 
                    size3, size1}, {column, 1, size1}}
                   }[[#]], 
                 VertexTextureCoordinates -> {{0, 0}, {1, 0}, {1, 
                    1}, {0, 1}}]} & /@ sel
             ];
           plsl3D = Show[Graphics3D[slices3D,
              BoxRatios -> ratio, ImageSize -> scale, 
              SphericalRegion -> True, Background -> back, 
              Lighting -> "Neutral",
              ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
              ViewAngle -> Dynamic[va],
              Boxed -> box, Axes -> axes, 
              AxesStyle -> Thread[List[{Red, Green, Blue}, Thick]]
              ], 
             PlotRange -> {{-1, size3 + 1}, {-1, size2 + 1}, {-1, 
                size1 + 1}}]
           ];
         ) // AbsoluteTiming // First;
     
     If[show == 6, plim3D = Show[
         im3D,
         BoxRatios -> ratio, ImageSize -> scale, 
         SphericalRegion -> True, Background -> back, 
         Lighting -> "Neutral",
         ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
         ViewAngle -> Dynamic[va],
         Boxed -> box, Axes -> axes, 
         AxesStyle -> Thread[List[{Red, Green, Blue}, Thick]],
         PlotRange -> {{-1, size3 + 1}, {-1, size2 + 1}, {-1, 
            size1 + 1}}
         ];
      ];
     
     Switch[show,
      1, plax,
      2, plcor,
      3, plsag,
      4, imall,
      5, plsl3D,
      6, plim3D
      ]
     
     ,
     
     (*show what*)
     {{show, 4, "Plot Mode"}, {1 -> "Axial", 2 -> "Coronal", 
       3 -> "Sagital", 4 -> "All Planes", 5 -> "Planes 3D", 
       6 -> "Volume 3D"}},
     {{scale, 500, "Plot Size"}, psizes},
     {{back, Gray, "BackGround"}, 
      ColorSlider[#, ImageSize -> {Automatic, 15}] &},
     
     (*general 4D*)
     {{set, 1, "Set (4D)"}, 1, Dynamic[size1b], 1},
     (*general 3D 4D*)
     {{trans, False, "Transpose 4D"}, {True, False}},
     {{reverse, False, "Reverse slices"}, {True, False}},
     {{slice, Round[size1/2], "Axial"}, 1, Dynamic[size1], 1},
     {{row, Round[size2/2], "Coronal"}, 1, Dynamic[size2], 1},
     {{column, Round[size3/2], "Sagital"}, 1, Dynamic[size3], 1},
     
     (*all planes 4*)
     {{lines, True, "Show lines"}, {True -> "On", False -> "Off"}},
     
     (*planes color 1-5*)
     {{color, "BlackToWhite", "ColorFunction"}, colors},
     {{lstyle, 1, "lstyle"}, colfuncs},
     {{cfs, False, "Auto Scaling"}, {True -> "On", False -> "Off"}},
     (*on or off by auto scale*)
     {{min, mind, "Min"}, mind, 0.9 max},
     {{minclip, RGBColor[{0, 0, 0}], "MinClip"}, 
      ColorSlider[#, ImageSize -> {Automatic, 15}] &},
     {{max, maxd, "Max"}, 1.1 min, maxd},
     {{maxclip, RGBColor[{255, 255, 255}], "MaxClip"}, 
      ColorSlider[#, ImageSize -> {Automatic, 15}] &},
     
     (*3D general 5-6*)
     {{box, True, "Show box"}, {True, False}},
     {{axes, True, "Show axis"}, {True, False}},
     {{vp, 3.5 {0.384, 0.709, 0.591}, "Viewpoint"}, views, 
      ControlType -> SetterBar},
     
     (*3D planes 5*)
     {{planex, True, "Show plane x"}, {True, False}},
     {{opx, 1, "Opacity plane x"}, 0, 1, 0.1},
     {{planey, True, "Show plane y"}, {True, False}},
     {{opy, 1, "Opacity plane y"}, 0, 1, 0.1},
     {{planez, True, "Show plane z"}, {True, False}},
     {{opz, 1, "Opacity plane z"}, 0, 1, 0.1},
     
     (*3Dvol 6*)
     {{col3D, Automatic, "Colorfunction 3D"}, colors3D},
     {{min3D, mins, "min 3D"}, Dynamic[mins], max3D},
     {{max3D, maxs, "max 3D"}, min3D, Dynamic[2 maxs]},
     
     {{vp, 3.5 {0.384, 0.709, 0.591}, "ViewPoint"}, Dynamic[vp] &, 
      ControlType -> None},
     {{vv, {0, 0, 1}, "ViewVertical"}, Dynamic[vv] &, 
      ControlType -> None},
     {{va, 25 Degree, "ViewAngle"}, Dynamic[va] &, 
      ControlType -> None},
     
     ControlPlacement -> Right,
     SynchronousInitialization -> False
     ];

   plotwindow3D = 
    CreateWindow[
     DialogNotebook[{CancelButton["Close",Clear[data]; DialogReturn[]], pan}, 
      WindowSize -> All, WindowTitle -> "Plot data window"]];
 ]
*)

(*
PointsFunc = 
  Compile[{{qual, _Real, 0}, {dx, _Integer, 0}, {dy, _Integer, 
     0}, {dz, _Integer, 0}, {size, _Integer, 0}, {alpha, _Real, 
     0}, {beta, _Real, 0}, {or, _Real, 1}},
   Block[{pts, pt = {1, 1, 1}, ptls, blank = {{1, 1, 1}}, test = 0, 
     step = Round[size/(qual*size)]},
    pts = Table[
      Round[(({{1, 0, 0}, {0, Cos[alpha], -Sin[alpha]}, {0, 
              Sin[alpha], Cos[alpha]}}.{{Cos[beta], 0, Sin[beta]}, {0,
               1, 0}, {-Sin[beta], 0, Cos[beta]}}).({x, y, or[[3]]} - 
            or)) + or],
      {y, -Round[(size - dy)/2], Round[(size - dy)/2] + dy, step},
      {x, -Round[(size - dx)/2], Round[(size - dx)/2] + dx, step}];
    Do[pts = DeleteCases[(
          ptls = #;
          blank = ConstantArray[0, Dimensions[ptls]];
          
          test = Total[(pt = #; 
               If[1 - step <= pt[[1]] <= dx + step && 
                 1 - step <= pt[[2]] <= dy + step && 
                 1 - step <= pt[[3]] <= dz + step, 1, 0]) & /@ ptls];
          If[test > 0, ptls, blank]
          ) & /@ pts, blank];
     pts = Transpose[DeleteCases[(
           ptls = #;
           blank = ConstantArray[0, Dimensions[ptls]];
           
           test = Total[(pt = #; 
                If[1 - step <= pt[[1]] <= dx + step && 
                  1 - step <= pt[[2]] <= dy + step && 
                  1 - step <= pt[[3]] <= dz + step, 1, 0]) & /@ ptls];
           If[test > 0, ptls, blank]
           ) & /@ Transpose[pts], blank]];
     , {2}];
    pts
    ]];

AngCor = Function[{ang, cor}, 
   If[-45 <= ang <= 45, ang*(cor/45), 
     If[Positive[ang], 
      cor + (ang - 45)*(90 - cor)/45, -cor + (ang + 45)*(90 - cor)/
         45]] Degree];


SyntaxInformation[PlotData3D] = {"ArgumentsPattern" -> {_, _.}};

PlotData3D[data_, vox : {_, _, _} : {1, 1, 1}] := 
 Module[{tab1, tab2, tab3, tab4, control, dim, dz, dx, dy, planes, 
   qual, xx, yy, zz, planex, opx, planey, opy, planez, opz, ratio, or,
    clip, slicex, slicey, slicez, 
   slicea, sliceao, surf, vp, va, vv, plot, exp, diag, dorig, dqual, 
   opd, diagx, diagy, diagz, alpha, beta, box, axes, label, ps, color,
    lstyle, bcol, legend, min, max, minclip, maxclip, 
   pxmin, pxmax, pymin, pymax, pzmin, pzmax, iso, isoval, surfqual, 
   surfop, surfcol, fileType, size, pannel, dimq, pts, pol, gr, vec, 
   a, b, c, ang, angx, angy}, 
  If[(! ArrayQ[data, _, NumericQ]) || (! ArrayDepth[data] == 3), Return[Message[PlotData3D::data]]];
  
  dim = {dz, dy, dx} = Dimensions[data];
  size = Round[Sqrt[Total[#^2 & /@ Drop[Sort[dim], 1]]] // N];
  
  tab1 = Column[{ManPannel[
      "Planes", {
      	{"Show Planes",Control@{{planes, True, ""}, {True, False}}}, 
      	{"Plane Quality (%)", Control@{{qual, .5, ""}, .2, 1, .2}}, 
      	{Style["Plane Navigation", Bold], ""}, 
        {"Corronal Slice (x)",Control@{{xx, Round[dx/2], ""}, 1, dx, 1}}, 
        {"Saggital Slice (y)",Control@{{yy, Round[dy/2], ""}, 1, dy, 1}}, 
        {"Axial Slice (z)",Control@{{zz, Round[dz/2], ""}, 1, dz, 1}}, 
        {Style["Plane Settings", Bold], ""}, 
        {"Show Corronal (x)",Control@{{planex, True, ""}, {True, False}}}, 
        {"Plane Opacity (x)",Control@{{opx, 1, ""}, 0, 1, 0.1}}, 
        {"Show Saggital (y)",Control@{{planey, True, ""}, {True, False}}}, 
        {"Plane Opacity (y)",Control@{{opy, 1, ""}, 0, 1, 0.1}}, 
        {"Show Axial (z)",Control@{{planez, True, ""}, {True, False}}}, 
        {"Plane Opacity (z)",Control@{{opz, 1, ""}, 0, 1, 0.1}}}], 
     ManPannel[
      "Diagonal plane", {
      	{"Show Diagonal",  Control@{{diag, False, ""}, {True,False}}}, 
        {"Show Plane Origin",Control@{{dorig, False, ""}, {True,False}}}, 
        {"Plane Quality (%)", Control@{{dqual, .5, ""},.2,1,.2}}, 
        {"Plane Opacity", Control@{{opd, 1, ""}, 0, 1, 0.1}}, 
        {Style["Position", Bold], ""}, 
        {"Diagonal x Position", Control@{{diagx, Round[dx/2], ""}, 1, dx, 1}}, 
        {"Diagonal y Position", Control@{{diagy, Round[dy/2], ""}, 1, dy,1}}, 
        {"Diagonal z Position",  Control@{{diagz, Round[dz/2], ""}, 1, dz, 1}}, 
        {Style["Rotation", Bold], ""}, 
        {"Rotation x-axis (\[Degree])", Control@{{alpha, 15, ""}, -90, 90,  1}}, 
        {"Rotation y-axis (\[Degree])",  Control@{{beta, 15, ""}, -90, 90, 1}}}, False]}];
  tab2 = Column[{ManPannel[
      "Plot Style", {
      	{"Show Box", Control@{{box, False, ""}, {True, False}}}, 
        {"Show Axes", Control@{{axes, False, ""}, {True, False}}}, 
        {"Plot Title", Control@{{label, "", ""}, InputField[#, String] &}}, 
        {"Plot Size", Control@{{ps, 400, ""}, sizes, ControlType -> PopupMenu}},
        {"Color function", Control@{{color, "GrayTones", ""}, colors, ControlType -> PopupMenu}}, 
        {"Color style", Control@{{lstyle, 1, ""}, colfuncs}}, 
        {"Background Color", Control@{{bcol, Gray, ""}, ColorSlider[#,  ImageSize -> {Automatic, 15}] &}}, 
        {"Legend on/off", Control@{{legend, False, ""}, {True, False}}}}], 
     ManPannel[
      "Plot Range", {
      	{"Min value", Control@{{min, Min[data], ""}, Min[data], max, (max - Min[data])/100}}, 
        {"Max value",  Control@{{max, Max[data], ""}, min, Max[data], (Max[data] - min)/100}}, 
        (*{"Transparent Clipping", Control@{{transclip, False, ""}, {True, False}}},*) 
        {"Min Clipping",  Control@{{minclip, Black, ""},  ColorSlider[#,  ImageSize -> {Automatic, 15}] &}}, 
        {"Max Clipping", Control@{{maxclip, White, ""},   ColorSlider[#, ImageSize -> {Automatic, 15}] &}}}], 
     ManPannel[
      "Slice Range", {
      	{"Minimal x value", Control@{{pxmin, 1, ""}, 1, pxmax, 1}}, 
        {"Maximal x value", Control@{{pxmax, dx, ""}, pxmin + 1, dx,1}}, 
        {"Minimal y value", Control@{{pymin, 1, ""}, 1, pymax, 1}}, 
        {"Maximal y value", Control@{{pymax, dy, ""}, pymin + 1, dy, 1}}, 
        {"Minimal z value", Control@{{pzmin, 1, ""}, 1, pzmax, 1}}, 
        {"Maximal z value", Control@{{pzmax, dz, ""}, pzmin + 1, dz, 1}}}]}];
  tab3 = Column[{ManPannel[
      "Iso Surface", {
      	{"Show IsoSurface", Control@{{iso, False, ""}, {True, False}}}, 
        {"Iso Value",  Control@{{isoval, {Round[Max[data]/2]}, ""}, InputField[#] &}}, 
        {"Surface Quality (%)", Control@{{surfqual, 50, ""}, 20, 125, 1}}, 
        {"Surface Opacity", Control@{{surfop, 1, ""}, 0, 1, 0.1}}, 
        {"Surface Color", Control@{{surfcol, Darker[Red], ""}, ColorSlider[#, ImageSize -> {Automatic, 15}] &}}}]}];
  tab4 = Column[{ManPannel[
      "Export plot", {
      	{"File Type", Control@{{fileType, ".jpg", ""}, files}}, 
        {"Export Size", Control@{{size, 400, ""}, sizes, ControlType -> PopupMenu}}, 
        {"Export", Button["Save Plot", FileSave[exp, fileType, size],  Method -> "Queued", ImageSize -> 150]}}]}
        ];
  control = {{{pannel, 1, ""}, {1 -> "Planes", 
      2 -> "Plotting Options", 3 -> "IsoSurface", 4 -> "Export"}}, 
    Delimiter, 
    PaneSelector[{1 -> tab1, 2 -> tab2, 3 -> tab3, 4 -> tab4}, 
     pannel]};
  
  (*mind=If[Min[data//N]==0.,0.0001 Max[data],0.];*)
  PrintTemporary["Initializing plot window, please wait"];
  
  Manipulate[If[! ListQ[data], Return[]];
    If[! ArrayQ[data], Return[]];
    (*define box ratio*)
    
    ratio = {(pxmax - pxmin), (pymax - pymin), (pzmax - 
         pzmin)} Reverse[vox];
    ang = {angx, angy} = 
      N[{ArcTan[(ratio[[3]]/vox[[1]])/(ratio[[3]]/vox[[2]])], 
         ArcTan[(ratio[[3]]/vox[[1]])/(ratio[[3]]/vox[[3]])]}/Degree];
    
    (*diagonal slice parameters*)
    
    or = {diagx, diagy, diagz};
    vec = {a, b, c} = 
      Normalize[
       Reverse[vox] {Sin[
          AngCor[beta, angy]], -Cos[AngCor[beta, angy]] Sin[
           AngCor[alpha, angx]], 
         Cos[AngCor[alpha, angx]] Cos[AngCor[beta, angy]]}];
    
    clip = {minclip,maxclip};
     (*If[transclip, {Transparent, Transparent}, {minclip, maxclip}];*)
    With[{
      
      (*slice plot function*)
      
      SlicePlot = 
       Function[{n, op}, 
        dimq = Round[qual*dim[[Drop[{1, 2, 3}, {n}]]]];
        {Opacity[op], 
          Texture[Graphics[
            Raster[Clip[
              Rescale[
               RescaleImg[{data[[zz]], data[[All, yy, All]], 
                  data[[All, All, xx]]}[[n]], dimq], {min, max}], {0, 
               1}, {-1, -2}], 
             ColorFunction -> (ColSelC[#, clip, {lstyle, color}] &)], 
            PlotRange -> {{0, dimq[[2]]}, {0, dimq[[1]]}}]], 
          Polygon[{{{1, 1, zz}, {dx, 1, zz}, {dx, dy, zz}, {1, dy, 
               zz}}, {{1, yy, 1}, {dx, yy, 1}, {dx, yy, dz}, {1, yy, 
               dz}}, {{xx, 1, 1}, {xx, dy, 1}, {xx, dy, dz}, {xx, 1, 
               dz}}}[[n]], 
           VertexTextureCoordinates -> {{0, 0}, {1, 0}, {1, 1}, {0, 
              1}}]}]
      
     (*, SlicePlotAng = 
       Function[{col, op}, 
        pts = PointsFunc[dqual, dx, dy, dz, size, AngCor[alpha, angx],
           AngCor[beta, angy], {diagx, diagy, diagz}];
        pol = {First[First[pts]], First[Last[pts]], Last[Last[pts]], 
          Last[First[pts]]};
          
        gr = 
         Graphics[
          Raster[Clip[
            Rescale[
             Map[data[[Clip[#[[3]], {1, dz}], Clip[#[[2]], {1, dy}], 
                Clip[#[[1]], {1, dx}]]] &, pts, {2}], {min, max}], {0,
              1}, {-1, -2}], 
           ColorFunction -> (ColSelC[#, clip, {lstyle, col}] &)], 
          PlotRange -> {{0, Length[pts[[1]]]}, {0, Length[pts]}}];
        {Opacity[op], Texture[gr], 
          Polygon[pol, 
           VertexTextureCoordinates -> {{0, 0}, {0, 1}, {1, 1}, {1, 
              0}}]}]*)},
     
     (*Draw planes*)
     
     slicex = If[#1 && #2, SlicePlot[2, opx], {}] &;
     slicey = If[#1 && #2, SlicePlot[3, opy], {}] &;
     slicez = If[#1 && #2, SlicePlot[1, opz], {}] &;
     slicea = If[#, 
     	pts = PointsFunc[dqual, dx, dy, dz, size, AngCor[alpha, angx],
           AngCor[beta, angy], {diagx, diagy, diagz}];
        pol = {First[First[pts]], First[Last[pts]], Last[Last[pts]], 
          Last[First[pts]]};
        gr = Graphics[Raster[Clip[Rescale[
             Map[data[[Clip[#[[3]], {1, dz}], Clip[#[[2]], {1, dy}], 
                Clip[#[[1]], {1, dx}]]] &, pts, {2}], {min, max}], {0,
              1}, {-1, -2}], 
           ColorFunction -> (ColSelC[#, clip, {lstyle, color}] &)], 
          PlotRange -> {{0, Length[pts[[1]]]}, {0, Length[pts]}}];
        {Opacity[opd], Texture[gr], 
          Polygon[pol, 
           VertexTextureCoordinates -> {{0, 0}, {0, 1}, {1, 1}, {1, 
              0}}]}
  , {}] &;];
    
    (*Draw diagonal slice marker*)
    
    sliceao = 
     If[#, Dynamic[{Red, 
         Scale[Sphere[or, (1/20 Min[dim*vox])], Min[vox]/Reverse[vox],
           or], Green, Thick, Arrowheads[0.05], 
         Scale[Arrow[
           Tube[{or, 
             or + (.5 Min[dim*vox]) (Normalize[vec/Reverse[vox]])}, 
            0.8]], Min[vox]/Reverse[vox], or]}], {}] &;
    
    (*Draw iso surface*)
    
    surf = 
     If[#, Dynamic[
        ListContourPlot3D[data, Contours -> Cases[isoval, _?NumberQ], 
          Mesh -> False, Axes -> False, 
          ContourStyle -> Directive[Opacity[surfop], surfcol], 
          MaxPlotPoints -> Round[0.75 surfqual], 
          BoundaryStyle -> None][[1]]], {}] &;
    
    (*Generate Plot*)
    
    plot = 
     Dynamic[Graphics3D[{sliceao[dorig], surf[iso], slicez[planes, planez], 
       slicey[planes, planey], slicex[planes, planex], slicea[diag]}, 
      Lighting -> "Neutral", BoxRatios -> ratio, 
      ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
      ViewAngle -> Dynamic[va], ImageSize -> ps, Background -> bcol, 
      SphericalRegion -> True, Boxed -> box, Axes -> axes, 
      AxesStyle -> Thread[List[{Red, Green, Blue}, Thick]], 
      BaseStyle -> {FontWeight -> Bold, FontFamily -> "sans-serif", 
        28}, LabelStyle -> 14, AxesLabel -> {"X", "Y", "Z"}, 
      PlotRange -> {{pxmin - 1, pxmax + 1}, {pymin - 1, 
         pymax + 1}, {pzmin - 1, pzmax + 1}}, 
      ContentSelectable -> True, PlotLabel -> label, 
      ImagePadding -> {{5, 5}, {5, 5}}]];
    
    (*Display Plot*)
    
    exp = If[legend, Legendi[plot, {lstyle, color, bcol}, min, max, ps], plot]
    
    (*Insert control pannels*)
    , ## ,
    (*Hidden manipulation parameters*)
    {{vp, {1.3, -2.4, 2}, "ViewPoint"}, Dynamic[vp] &, ControlType -> None}, 
    {{vv, {0, 0, 1}, "ViewVertical"}, Dynamic[vv] &, ControlType -> None}, 
    {{va, 25 Degree, "ViewAngle"}, Dynamic[va] &, ControlType -> None}, 
    Deployed -> False, 
    SynchronousInitialization -> False, 
    SynchronousUpdating -> False, 
    ContinuousAction -> False,
    ControlPlacement -> Right
    ]&@@ control]
*)


(*
PlotData3D[data_, vox:{_,_,_}:{1,1,1}] := 
 Module[{tab1, tab2, tab3, tab4, control, dim, dz, dx, dy,
   planes, qual, xx, yy, zz, planex, opx, planey, opy, planez, opz,
   ratio,or,clip,CR=Clip[Round[#1], {1, #2}]&,
   slicex,slicey,slicez,slicea,sliceao,surf,vp,va,vv,plot,exp,
   diag, dorig, dqual, opd, diagx, diagy, diagz, alpha, beta,
   box, axes, label, ps, color, lstyle, bcol, legend,
   min, max, transclip, minclip, maxclip,
   pxmin, pxmax, pymin, pymax, pzmin, pzmax,
   iso, isoval, surfqual, surfop, surfcol,
   fileType, size, pannel,dimq,
   vec, a, b, c, ang, angx, angy, afunc, arang, asel, v1, v2
   },
  
  If[(! ArrayQ[data, _, NumericQ]) || (! ArrayDepth[data] == 3),Return[Message[PlotData3D::data]]];
  
  dim = {dz, dy, dx} = Dimensions[data];
  
  tab1 = Column[{
     ManPannel["Planes", {
       {"Show Planes", Control@{{planes, True, ""}, {True, False}}},
       {"Plane Quality (%)", Control@{{qual, .5, ""}, .1, 1, .1}},
       {Style["Plane Navigation", Bold], ""},
       {"Corronal Slice (x)", 
        Control@{{xx, Round[dx/2], ""}, 1, dx, 1}},
       {"Saggital Slice (y)", 
        Control@{{yy, Round[dy/2], ""}, 1, dy, 1}},
       {"Axial Slice (z)", Control@{{zz, Round[dz/2], ""}, 1, dz, 1}},
       {Style["Plane Settings", Bold], ""},
       {"Show Corronal (x)", 
        Control@{{planex, True, ""}, {True, False}}},
       {"Plane Opacity (x)", Control@{{opx, 1, ""}, 0, 1, 0.1}},
       {"Show Saggital (y)", 
        Control@{{planey, True, ""}, {True, False}}},
       {"Plane Opacity (y)", Control@{{opy, 1, ""}, 0, 1, 0.1}},
       {"Show Axial (z)", Control@{{planez, True, ""}, {True, False}}},
       {"Plane Opacity (z)", Control@{{opz, 1, ""}, 0, 1, 0.1}}
       }]
     ,
     ManPannel["Diagonal plane", {
       {"Show Diagonal", Control@{{diag, False, ""}, {True, False}}},
       {"Show Plane Origin", 
        Control@{{dorig, False, ""}, {True, False}}},
       {"Plane Quality (%)", Control@{{dqual, .5, ""}, .1, 1, .1}},
       {"Plane Opacity", Control@{{opd, 1, ""}, 0, 1, 0.1}},
       {Style["Position", Bold], ""},
       {"Diagonal x Position", 
        Control@{{diagx, Round[dx/2], ""}, 1, dx, 1}},
       {"Diagonal y Position", 
        Control@{{diagy, Round[dy/2], ""}, 1, dy, 1}},
       {"Diagonal z Position", 
        Control@{{diagz, Round[dz/2], ""}, 1, dz, 1}},
       {Style["Rotation", Bold], ""},
       {"Rotation x-axis (\[Degree])", 
        Control@{{alpha, 15, ""}, -90, 90, 1}},
       {"Rotation y-axis (\[Degree])", 
        Control@{{beta, 15, ""}, -90, 90, 1}}
       }, False]
     }];
  
  tab2 = Column[{
     ManPannel["Plot Style", {
       {"Show Box", Control@{{box, False, ""}, {True, False}}},
       {"Show Axes", Control@{{axes, False, ""}, {True, False}}},
       {"Plot Title", 
        Control@{{label, "", ""}, InputField[#, String] &}},
       {"Plot Size", 
        Control@{{ps, 400, ""}, sizes, ControlType -> PopupMenu}},
       {"Color function", 
        Control@{{color, "GrayTones", ""}, colors, 
          ControlType -> PopupMenu}},
       {"Color style", Control@{{lstyle, 1, ""}, colfuncs}},
       {"Background Color", 
        Control@{{bcol, Gray, ""}, 
          ColorSlider[#, ImageSize -> {Automatic, 15}] &}},
       {"Legend on/off", Control@{{legend, False, ""}, {True, False}}}
       }]
     ,
     ManPannel["Plot Range", {
       {"Min value", 
        Control@{{min, Min[data], ""}, Min[data], 
          max, (max - Min[data])/100}},
       {"Max value", 
        Control@{{max, Max[data], ""}, min, 
          Max[data], (Max[data] - min)/100}},
       {"Transparent Clipping", 
        Control@{{transclip, False, ""}, {True, False}}},
       {"Min Clipping", 
        Control@{{minclip, Black, ""}, 
          ColorSlider[#, ImageSize -> {Automatic, 15}] &}},
       {"Max Clipping", 
        Control@{{maxclip, White, ""}, 
          ColorSlider[#, ImageSize -> {Automatic, 15}] &}}
       }]
     ,
     ManPannel["Slice Range", {
       {"Minimal x value", Control@{{pxmin, 1, ""}, 1, pxmax, 1}},
       {"Maximal x value", 
        Control@{{pxmax, dx, ""}, pxmin + 1, dx, 1}},
       {"Minimal y value", Control@{{pymin, 1, ""}, 1, pymax, 1}},
       {"Maximal y value", 
        Control@{{pymax, dy, ""}, pymin + 1, dy, 1}},
       {"Minimal z value", Control@{{pzmin, 1, ""}, 1, pzmax, 1}},
       {"Maximal z value", Control@{{pzmax, dz, ""}, pzmin + 1, dz, 1}}
       }]
     }];
  
  tab3 = Column[{
     ManPannel["Iso Surface", {
       {"Show IsoSurface", Control@{{iso, False, ""}, {True, False}}},
       {"Iso Value", 
        Control@{{isoval, {Round[Max[data]/2]}, ""}, InputField[#] &}},
       {"Surface Quality (%)", 
        Control@{{surfqual, 50, ""}, 20, 125, 1}},
       {"Surface Opacity", Control@{{surfop, 1, ""}, 0, 1, 0.1}},
       {"Surface Color", 
        Control@{{surfcol, Darker[Red], ""}, 
          ColorSlider[#, ImageSize -> {Automatic, 15}] &}}
       }]
     }];
  
  tab4 = Column[{
     ManPannel["Export plot", {
       {"File Type", Control@{{fileType, ".jpg", ""}, files}},
       {"Export Size", 
        Control@{{size, 400, ""}, sizes, ControlType -> PopupMenu}},
       {"Export", 
        Button["Save Plot", FileSave[exp, fileType, size], 
         Method -> "Queued", ImageSize -> 150]}
       }]
     }];
  
  control = {
    {{pannel, 1, ""}, {1 -> "Planes", 2 -> "Plotting Options", 
      3 -> "IsoSurface", 4 -> "Export"}},
    Delimiter,
    PaneSelector[{1 -> tab1, 2 -> tab2, 3 -> tab3, 4 -> tab4}, pannel]
    };
  
  (*mind = If[Min[data // N] == 0., 0.0001 Max[data], 0.];*)
  
  PrintTemporary["Initializing plot window, please wait"];
  
  Manipulate[
  	
  	If[!ArrayQ[data],Return[]];
     
     (*define box ratio*)
     ratio = {(pymax - pymin), (pxmax - pxmin), (pzmax - pzmin)} Reverse[vox];
     
     (*diagonal slice parameters*)
     or = {diagx, diagy, diagz};
     vec = {a, b, c} = Normalize[Reverse[vox] {Sin[alpha Degree] Cos[beta Degree], Sin[alpha Degree] Sin[beta Degree], Cos[alpha Degree]}];
     If[diag, 
     	ang = {angx, angy} = N[{ArcTan[dz/dx], ArcTan[dx/dy]}/Degree];
      	asel = If[Abs[alpha] > angx && Abs[beta] > angy, 1, If[Abs[alpha] > angx, 2, 3]];
      	afunc = {
      		{v1,(-a (v1 - diagx) - b (-diagy) - c (v2 - diagz))/If[b == 0 || b == 0., 1, b], v2}, 
          	{(-a (-diagx) - b (v1 - diagy) - c (v2 - diagz))/If[a == 0 || a == 0., 1, a], v1, v2}, 
          	{v1, v2, (-a (v1 - diagx) - b (v2 - diagy) - c (-diagz))/If[c == 0 || c == 0., 1, c]}
          	}[[asel]]// N;
      	arang = {{dx, dz}, {dy, dz}, {dx, dy}}[[asel]];
      	;
      ];
     
     clip=If[transclip, {Transparent,Transparent}, {minclip,maxclip}];
     
     With[{
        (*slice plot function*)
        SlicePlot = Function[{vecf, rang, op, qualf},
        	ParametricPlot3D[vecf, {v1, 1, rang[[1]]}, {v2, 1, rang[[2]]},
        		PlotStyle -> Opacity[op], PlotPoints -> Round[qualf rang], Mesh -> False, ColorFunctionScaling -> False,
        		ColorFunction -> (ColSelC[Clip[Rescale[data[[CR[#3,dz],CR[#1,dx],CR[#2,dy]]],{min,max}],{0,1},{-1,2}],clip,{lstyle,color}]&)
        	]]
        ,
    	SlicePlot2 = 
  Function[{n, op}, dimq = Round[qual*dim[[Drop[{1, 2, 3}, {n}]]]];
   Graphics3D[{Opacity[op], 
     Texture[Graphics[
       Raster[Clip[
         Rescale[RescaleImg[{data[[zz]], data[[All, All, yy]], 
             data[[All, xx, All]]}[[n]], dimq], {min, max}], {0, 
          1}, {-1, -2}], 
        ColorFunction -> (ColSelC[#, clip, {lstyle, color}] &)], 
       PlotRange -> {{0, dimq[[2]]}, {0, dimq[[1]]}}, 
       AspectRatio -> (Divide @@ (dimq*Drop[vox, {n}]))]],
     Polygon[{
        {{1, 1, zz}, {dim[[2]], 1, zz}, {dim[[2]], dim[[3]], zz}, {1, 
          dim[[3]], zz}},
        {{1, xx, 1}, {dim[[2]], xx, 1}, {dim[[2]], xx, dim[[1]]}, {1, 
          xx, dim[[1]]}},
        {{yy, 1, 1}, {yy, dim[[3]], 1}, {yy, dim[[3]], dim[[1]]}, {yy,
           1, dim[[1]]}}
        }[[n]], 
      VertexTextureCoordinates -> {{{0, 0}, {0, 1}, {1, 1}, {1, 
           0}}, {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, {{0, 0}, {1, 0}, {1,
            1}, {0, 
           1}}}[[n]]]}]]},
           (*Draw planes*)
		(*slicex=If[#1&&#2,SlicePlot[{v1,xx,v2},{dy,dz},opx,qual],Graphics3D[]]&;
		slicey=If[#1&&#2,SlicePlot[{yy,v1,v2},{dx,dz},opy,qual],Graphics3D[]]&;
		slicez=If[#1&&#2,SlicePlot[{v1,v2,zz},{dx,dy},opz,qual],Graphics3D[]]&;*)
		slicex = If[#1 && #2, SlicePlot2[2, opx], Graphics3D[]] &;
		slicey = If[#1 && #2, SlicePlot2[3, opy], Graphics3D[]] &;
		slicez = If[#1 && #2, SlicePlot2[1, opz], Graphics3D[]] &;
      slicea = If[#, SlicePlot[afunc, arang, opd, dqual], Graphics3D[]]&;
      ];
     
     (*Draw diagonal slice marker*)
     sliceao = If[#, Graphics3D[{Red,Scale[Sphere[or, 3], Min[vox]/Reverse[vox], or], Green, Thick, Arrowheads[0.05], 
         Scale[Arrow[Tube[{or, or + (1/5 Min[dim*vox]) (Normalize[vec/Reverse[vox]])}, 0.8]], Min[vox]/Reverse[vox], or]}],
       Graphics3D[]]&;
     
     (*Draw iso surface*)
     surf = If[#, ListContourPlot3D[Transpose[data, {1, 3, 2}], Contours -> Cases[isoval, _?NumberQ],Mesh -> False, Axes -> False,
        ContourStyle -> Directive[Opacity[surfop], surfcol], MaxPlotPoints -> Round[0.75 surfqual], BoundaryStyle -> None],
       Graphics3D[]]&;
     
     (*Generate Plot*)
     plot = Show[
       sliceao[dorig], surf[iso], slicez[planes,planez], slicey[planes,planey], slicex[planes,planex], slicea[diag],
       Lighting -> "Neutral", BoxRatios -> ratio, 
       ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
       ViewAngle -> Dynamic[va], ImageSize -> ps, Background -> bcol, 
       SphericalRegion -> True, Boxed -> box, Axes -> axes, 
       AxesStyle -> Thread[List[{ Green, Red, Blue}, Thick]], 
       BaseStyle -> {FontWeight -> Bold, FontFamily -> "sans-serif", 
         28}, LabelStyle -> 14, AxesLabel -> { "Y", "X", "Z"}, 
       PlotRange -> {{pymin - 1, pymax + 1}, {pxmin - 1, 
          pxmax + 1}, {pzmin - 1, pzmax + 1}}, 
       ContentSelectable -> True, PlotLabel -> label, 
       ImagePadding -> {{5, 5}, {5, 5}}
       ];
     
     (*Display Plot*)  
     exp=If[legend,Dynamic[Legendi[plot,{lstyle,color,bcol},min,max,ps]],Dynamic[plot]]
       
     (*Insert control pannels*)  
     , ##,
     (*Hidden manipulation parameters*)
     {{vp, {1.3, -2.4, 2}, "ViewPoint"}, Dynamic[vp] &, ControlType -> None},
     {{vv, {0, 0, 1}, "ViewVertical"}, Dynamic[vv] &, ControlType -> None},
     {{va, 25 Degree, "ViewAngle"}, Dynamic[va] &, ControlType -> None},
     
     Deployed->True,
     SynchronousInitialization -> False,
     ControlPlacement -> Right,
     SynchronousUpdating -> False,
     ContinuousAction -> False] & @@ control
  ]
*)


(* ::Subsection:: *)
(*MakeUnet - old*)

(*

Options[MakeUnet] = {
	BlockType -> "ResNet", 
	DropoutRate -> 0.2, 
	NetworkDepth -> 5, 
	DownsampleSchedule -> Automatic, 
	FeatureSchedule -> Automatic,
	InputFilters -> 32, 
	ActivationType -> "GELU"
}

SyntaxInformation[MakeUnet] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeUnet[nClass_, dimIn_, opts:OptionsPattern[]] :=MakeUnet[1, nClass, dimIn, opts]

MakeUnet[nChan_, nClass_, dimIn_, OptionsPattern[]] := Block[{
		dep, dep1, drop, type, dim, nDim, filt, feat, enc, dec, stride, filtIn, actType
	},

	(*Get the options*)
	{dep, drop, type, stride, feat, filtIn, actType} = OptionValue[
		{NetworkDepth, DropoutRate, BlockType, DownsampleSchedule, FeatureSchedule, InputFilters, ActivationType}
	];

	(*Define UNET properties*)
	enc ="enc_" <> ToString[#]&;
	dec ="dec_" <> ToString[#]&;
	dep1 = dep-1;
	nDim = Length@dimIn;
	dim = Switch[nDim, 2, "2D", 3, "3D"];
	feat = If[feat===Automatic,
		Switch[type, "DenseNet" | "UDenseNet", {1, 2, 4, 6, 8}, _, {1, 2, 4, 8, 16}],
		feat];
	feat = PadRight[feat, dep, Last@feat];
	filt = Switch[type, "DenseNet" | "UDenseNet", Table[{filtIn, 1 + i}, {i, feat}], _, filtIn feat];
	stride = Prepend[If[stride===Automatic, ConstantArray[2, {dep-1, nDim}], stride], {1, 1, 1}[[;;nDim]]];

	(*make the UNET*)
	NetGraph[
		Association@Join[
			Table[
				enc[i] -> ConvNode[filt[[i]], "Dropout" -> drop, "Dimensions" -> dim, "Stride" -> stride[[i]],
					"ConvType" -> type, "NodeType" -> "Encode", "ActivationType" -> actType]
			, {i, 1, dep}],
			Table[
				dec[i] -> ConvNode[filt[[i]], "Dropout" -> drop, "Dimensions" -> dim, "Stride" -> stride[[i+1]],
					"ConvType" -> type, "NodeType" -> "Decode", "ActivationType" -> actType]
			, {i, 1, dep1}],
			{"start" -> UNetStart[filt[[1]], nChan, dimIn, actType]},
			{"map" -> UNetMap[dimIn, nClass]}
		],

		Join[
			{NetPort["Input"] -> "start" -> enc[1], {enc[dep], enc[dep1]} -> dec[dep1], dec[1] -> "map"},
			Table[enc[i - 1] -> enc[i], {i, 2, dep}],
			Table[{dec[i + 1], enc[i]} -> dec[i], {i, 1, dep-2}]
		]
	]
]


UNetMap[dim_, nClass_] :=  Flatten[{
	ConvolutionLayer[nClass, 1], If[nClass > 1,	
		{TransposeLayer[Switch[Length@dim, 2, {3, 1, 2}, 3, {4, 1, 2, 3}]], SoftmaxLayer[]},
		{LogisticSigmoid, FlattenLayer[1]}
	]
}]


UNetStart[filt_, nChan_, dimIn_, actType_] := {ConvolutionLayer[If[IntegerQ[filt],filt,First@filt], 1, "Input" -> Prepend[dimIn, nChan]], BatchNormalizationLayer[], ActivationLayer[actType]}


Options[ConvNode] = {
	"Dimensions" -> "3D",
	"ActivationType" -> "GELU",
	"Dropout" -> 0.2,
	"ConvType" -> "ResNet",
	"NodeType" -> "Encode",(*encode, decode, start*)
	"Stride" -> Automatic
};

ConvNode[chan_, OptionsPattern[]] := Block[{
		convType, nodeType, actType, mode, node, drop, dim, stride
	},

	(*get the options*)
	{convType, nodeType, actType, drop, dim, stride} = OptionValue[
		{"ConvType", "NodeType", "ActivationType", "Dropout", "Dimensions", "Stride"}
	];

	(*mode is encoding or decoding, decoding is solved later and treated as normal here*)
	mode = If[nodeType === "Encode", "down", "normal"];

	(*make convblocks for various convolution types*)
	node = Switch[convType,	
		"UResNet", 
		Flatten[{
			ConvBlock[chan/2, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], 
			ConvBlock[chan, "ActivationType" -> actType, "Stride"->stride]
		}],

		"ResNet",
		{<|
			"con" -> Join[
				ConvBlock[chan/2, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], 
				ConvBlock[chan, "ActivationType" -> "None"]], 
			"skip" -> ConvBlock[chan, "ConvMode" -> mode<>"S", "ActivationType" -> "None", "Stride"->stride],
			"tot" -> {TotalLayer[], ActivationLayer[actType]}
		|>, {
			{"con", "skip"} -> "tot"
		}},

		"DenseNet",
		With[{n = chan[[1]], dep = chan[[2]], layName = "lay_" <> ToString[#] &},{
			Join[
				<|If[mode === "down", "down" -> ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], Nothing]|>,
				Association@Table[If[rep==dep, "lay_end", layName[rep]] -> ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> "catenate"], {rep, 1, dep}]
			],
			Table[Table[If[rr == 0, If[mode==="down", "down", NetPort["Input"]], layName[rr]], {rr, 0, rep - 1}] -> If[rep==dep, "lay_end", layName[rep]], {rep, 1, dep}]
		}],

		"UDenseNet", 
		Flatten[{If[mode === "down", ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], Nothing], 
			ConstantArray[ConvBlock[chan[[1]], "ActivationType" -> actType], chan[[2]]]}],

		_,
		Flatten[{ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], ConvBlock[chan, "ActivationType" -> actType]}]
	];


	(*Add dropout and upconv for deconding block*)
	NetFlatten@If[nodeType === "Decode",

		(*convert to decoding block and add dropout*)
		NetGraph[<|
			"upconv" -> ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> "up", Dimensions -> dim, "Stride"->stride],
			"conv" -> If[convType==="ResNet"||convType==="DenseNet",
				NetGraph[
					Join[node[[1]], <|"cat"->CatenateLayer[],"drop"->DropoutLayer[drop]|>],
					Switch[convType,
						"ResNet", Join[node[[2]], {"cat"->{"con","skip"}, "tot"->"drop"}],
						"DenseNet", Join[node[[2]] /. NetPort["Input"]->"cat", {"lay_end"->"drop"}]
					]
				],
				Flatten[{CatenateLayer[], node, DropoutLayer[drop]}]
			]
		|>, {
			{NetPort["Input2"] -> "upconv", NetPort["Input1"]} -> "conv"
		}]
		,

		(*add dropout to encoding block*)
		If[convType==="ResNet"||convType==="DenseNet",
				NetGraph[
					Join[node[[1]], <|"drop"->DropoutLayer[drop]|>],
					Join[node[[2]], {Switch[convType,"ResNet", "tot", "DenseNet", "lay_end"]->"drop"}]
				],
				NetChain[Flatten@{node, DropoutLayer[drop]}]
			]
		
	]
]


Options[ConvBlock] = {
	"Dimensions" -> "3D",
	"ActivationType" -> "GELU",
	"ConvMode" -> "normal"(*normal, up, down, catenate*),
	"Stride" -> 2
};

ConvBlock[channels_, OptionsPattern[]] := Block[{
		chan, kern,  actType, pad, actLayer, convMode, dim, str
	},

	{actType, convMode, dim, str} = OptionValue[{"ActivationType", "ConvMode", "Dimensions", "Stride"}];
	chan = Round@First@Flatten@{channels};
	
	Switch[convMode,
		"up", 
		{ResizeLayer[Scaled/@str, Resampling -> "Nearest"], ConvolutionLayer[chan, 2, "PaddingSize" -> ConstantArray[{0,1},Length[str]], "Stride" -> 1]},
		"down"|"downS", 
		{ConvolutionLayer[chan, str, "PaddingSize" -> 0, "Stride" -> str], BatchNormalizationLayer[], ActivationLayer[actType]},
		"normal", 
		{ConvolutionLayer[chan, 3, "PaddingSize" -> 1, "Stride" -> 1], BatchNormalizationLayer[], ActivationLayer[actType]},
		"normalS", 
		{ConvolutionLayer[chan, 1, "PaddingSize" -> 0, "Stride" -> 1], BatchNormalizationLayer[], ActivationLayer[actType]},
		"catenate", 
		{CatenateLayer[], ConvolutionLayer[chan, 3, "PaddingSize" -> 1, "Stride" -> 1], BatchNormalizationLayer[], ActivationLayer[actType]}
	]
]


ActivationLayer[actType_] := If[StringQ[actType],
	Switch[actType, "LeakyRELU", ParametricRampLayer[], "None", Nothing, _, ElementwiseLayer[actType]],
	actType
]


*)



(* ::Subsection::Closed:: *)
(*ROIMask*)


SyntaxInformation[ROIMask] = {"ArgumentsPattern" -> {_, _, _.}};

ROIMask[roiDim_, maskdim_,ROI:{(_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..})..}]:=
Module[{output},
	output=Map[#[[1]]->ROIMask[roiDim,maskdim,#[[2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[roiDim_,maskdim_,ROI:{{_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}}..}]:=
Module[{output},
	output=Map[#[[1,1]]->ROIMask[roiDim,maskdim,#[[1,2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[roiDim_,maskdim_,ROI:{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}]:=
Module[{output,roiCor,roiSlice,msk},
	output=ConstantArray[0,Join[{roiDim[[1]]},maskdim]];
	If[ROI[[All,1]]!={{{0,0}}},
		roiCor=Round[ROI[[All,1]]];
		roiSlice=Clip[ROI[[All,2]],{1,roiDim[[1]]}];
		msk=1-ImageData[Image[Graphics[Polygon[#],PlotRange->{{0,roiDim[[3]]},{0,roiDim[[2]]}}],"Bit",ColorSpace->"Grayscale",ImageSize->maskdim]]&/@roiCor;
		MapIndexed[output[[#1]]=msk[[First[#2]]];&,roiSlice];
		];
	Return[output];
	]

ROIMask[maskdim_,ROI:{(_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..})..}]:=
Module[{output},
	output=Map[#[[1]]->ROIMask[maskdim,#[[2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[maskdim_,ROI:{{_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}}..}]:=
Module[{output},
	output=Map[#[[1,1]]->ROIMask[maskdim,#[[1,2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[maskdim_,ROI:{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}]:=
Module[{output, roiCor, roiSlice, msk},
 output = ConstantArray[0, maskdim];
 If[ROI[[All, 1]] != {{{0, 0}}},
  roiCor = Round[Map[Reverse[maskdim[[2 ;; 3]]]*# &, ROI[[All, 1]], {2}]];
  If[Max[ROI[[All, 2]]] > maskdim[[1]], Message[ROIMask::war]];
  roiSlice = Clip[ROI[[All, 2]], {1, maskdim[[1]]}];
  msk = 1 - 
      ImageData[
       Image[Graphics[Polygon[#], 
         PlotRange -> {{0, maskdim[[3]]}, {0, maskdim[[2]]}}], "Bit", 
        ColorSpace -> "Grayscale", 
        ImageSize -> maskdim[[2 ;; 3]]]] & /@ roiCor;
  MapIndexed[output[[#1]] = msk[[First[#2]]]; &, roiSlice];
  ];
 Return[output];]


(* ::Subsection::Closed:: *)
(*SetupDataStructure*)


SetupDataStructure[dcmFolder_] := 
 Module[{folderdcm, foldernii, folderout, folders,fol, niiFolder, outFolder},
  folderdcm = Directory[] <> $PathnameSeparator <> # & /@ Select[FileNames["*", "dcm"], DirectoryQ];

  foldernii = StringReplace[#, "dcm" -> "nii"] & /@ folderdcm;
  folderout = StringReplace[#, "dcm" -> "out"] & /@ folderdcm;
  folders = Transpose[{folderdcm, foldernii, folderout}];
  
  fol = Last@FileNameSplit[dcmFolder];
  niiFolder = StringReplace[dcmFolder, fol -> "nii"];
  outFolder = StringReplace[dcmFolder, fol -> "out"];
  If[! DirectoryQ[niiFolder], CreateDirectory[niiFolder]];
  If[! DirectoryQ[outFolder], CreateDirectory[outFolder]];
  
  (*create nii files*)
  If[! DirectoryQ[#[[2]]], CreateDirectory[#[[2]]]; DcmToNii[#[[1 ;; 2]]]] & /@ folders;
  
  folders
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
