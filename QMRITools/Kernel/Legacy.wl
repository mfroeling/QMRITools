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



(* ::Subsection::Closed:: *)
(*Options*)


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
].Inverse[(GradRotMat[files[[1]]])],.00001]
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
gradsRot = Round[grads.gradRotmat, .0001];
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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
