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
"DTItoolExpTens[tensor] exports a diffusion tensor array to the DTItool format (*.dat).
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
"ReadDicomDirDiff[file] reads the image data and relevant diffusion parameters from a dicom directory."

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



BayesianIVIMFit2::usage = 
"BayesianIVIMFit2[data, bval, init, mask] performs bayesian IVIM fit of data.

data is the data which should be {slice, Ndiff, x, y}.
bval is the bvector would be length Ndiff.
init is the initialization of the bayesian fit which comes from IVIMCalc, (without s0 using 2 compartments).
mask is the region in which the bayesian fit is performed.

output is {f1, dc, pdc1}. The fraction is defined between 0 and 1, the dc, pdc1 is in mm^2/s."

BayesianIVIMFit3::usage = 
"BayesianIVIMFit3[data, bval, init, mask] performs bayesian IVIM fit of data.

data is the data which should be {slice, Ndiff, x, y}.
bval is the bvector would be length Ndiff.
init is the initialization of the bayesian fit which comes from IVIMCalC, (without s0 using 3 compartments).
mask is the region in which the bayesian fit is performed.

output is {f1, f2, dc, pdc1, pdc2}. The fractions f1 and f2 are defined between 0 and 1, the dc, pdc1 and pdc1 is in mm^2/s."

FracCorrect::usage = 
"FracCorrect[fraction, time] corrects the signal fraction calculated with the IVIM model for tissue relaxation and acquisition parameters.
After correction the signal fraction can be regarded as volume fraction.
FracCorrect[{fraction1, fraction2}, time] corrects the signal fraction1 and fraction2 from a 3 compartment IVIM model. 

time is {{te, tr}, {t2t, t21}, {t1t, t11}} or {{te, tr}, {t2t, t21, t22}, {t1t, t11, t12}}.
where t2t and t1t are \"tissue\" relaxation times and t11 t12, t21 and t22 the \"fluid\" relaxation times.

The te and tr as well as the relaxation times T2 and T1 can be defines in any time unit as long as they are consistant for all, e.g. all in ms.

output is the corrected fraction maps.";

ThetaConv::usage = 
"ThetaConv[{f1, Fc, pDc}] converts the parameters from Log space to normal space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConv[{f1, f2, dc, pDc1}] converts the parameters from Log space to normal space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConv[{f1, f2, dc, pDc1, pDc2}] converts the parameters from Log space to normal space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3."

ThetaConvi::usage = 
"ThetaConvi[{f, dc, pdc}] converts the parameters from Normal space to Log space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConvi[{f1, f2, dc, pdc1}] converts the parameters from Normal space to Log space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConvi[{f1, f2, dc, pdc1, pdc2}] converts the parameters from Normal space to Log space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3."

FConvert::usage = 
"FConvert[f] convers the fraction f from log space."

FConverti::usage = 
"FConverti[f] converts the fraction f to log space."

CorrectParMap::usage = 
"CorrectParMap[par, constraints, mask] removes the IVIM parameters outside the constraints within the mask.

par is {f1, dc, pdc1} or {f1, f2, dc, pdc1, pdc2}.
constraints are the lower and upper constraints for each parameters {{min, max},...}.
mask has the same dimensions as the parameter maps. 

output are the corrected paremeter maps."

HistogramPar::usage = 
"HistogramPar[data, {constraints, Nbins}, style, color, range] plots histograms of IVIM solution.
HistogramPar[data, {constraints, Nbins, mu, conv}, components, color, range] plots histograms of IVIM solution.

data is {f1, dc, pdc1} or {f1, f2, dc, pdc1, pdc2}.
constraints are the ranges of the x-axes for the plots.
Nbins are the number of histogram bins.
style is the plot type, can be 1, 2, or 3.
color is the color of the histogram.
range are the ranges of the y-axes.  

output is a row of histograms."


FitSpectra::usage = 
"FitSpectra[specBasis, spec, {st,end}, dt, {lwVals,lwAmps}] Fits the basis spectra from GetSpectraBasisFunctions to the spec overt the ppm range {st, end} and dt the dwell time."

FindSpectraPpmShift::usage = 
"FindSpectraPpmShift[spectra, {dw, gyro}, peaks] finds the ppm value that aligns the spectra with the given peak positions peaks which is a list of ppm values. 
FindSpectraPpmShift[spectra, {dw, gyro}, {peaks, amps}] finds the ppm value that aligns the spectra with the given peak positions peaks which is a list of ppm values and amps are their relative amplitudes. 
FindSpectraPpmShift[spectra, {dw, gyro}, specTar] finds the ppm value that aligns the spectra with the given target spectra specTar." 

FitSpectraResultTable::usage = 
"FitSpectraResultTable[parFit, parsF, names, ref, out] function not done."

CompareSpectraFitPlot::usage =
"CompareSpectraFitPlot[ppmPl, specPlot, fitPlot] function not done."

CompareFidFitPlot::usage =
"CompareFidFitPlot[time, fidPlot, fitPlot] function not done."

MakeSpectraResultPlot::usage = 
"MakeSpectraResultPlot[ppmF, specF, {fit, basisFit}, names, sc, met] function not done."

SpectraFitResult::usage = 
"SpectraFitResult[spec, {fit, basisFit}, te, {dw, gyro}, {pars, names, metRef, log}, plots, OptionsPattern[]] function not done."

CSIInterface::usage = 
"CSIInterface[] opens the CSI interface. Function not done.
CSIInterface[te, bw] opens the CSI interface with known te and bw.
CSIInterface[file] opens the CSI interface with the data from file loaded.
CSIInterface[file, {tei, bwi}] opens the CSI interface with the data from file loaded with known te and bw."


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



ChainSteps::usage = 
"ChainSteps is an option for BayesianIVIMFit2 and BayesianIVIMFit3. It determines how long the algorithm runs.
three values must be given {iterations, burn steps, sample density}."

UpdateStep::usage = 
"UpdateStep is an option for BayesianIVIMFit2 and BayesianIVIMFit3. It determines how often the parameters are updated. Is optimized during the first 500 burn steps.";

FixPseudoDiff::usage = 
"FixPseudoDiff is an option for BayesianIVIMFit2 and BayesianIVIMFit3. If the pDc1 and pD2 were fixed in IVIMCalc this value should be True."

FixPseudoDiffSD::usage = 
"FixPseudoDiffSD is an option for BayesianIVIMFit2 and BayesianIVIMFit3. Gives the standard deviation of pDc1 and pD2 if FixPseudoDiff is True."

CorrectPar::usage = 
"CorrectPar is an option for BayesianIVIMFit2 and BayesianIVIMFit3. If True it removes the values outside the constraints using CorrectParMap."

FitConstrains::usage = 
"FitConstrains is an option for BayesianIVIMFit2 and BayesianIVIMFit3. Gives the constraints of the parameters. 
The values are used for displaying the histograms and for the initialization if CorrectPar is True."

OutputSamples::usage = 
"OutputSamples is an option for BayesianIVIMFit2 and BayesianIVIMFit3. If set True the full marcov chain is given as an additionaln output."



SplineSpacingFactor::usage = 
"SplineSpacingFactor is an option for FitSpectra and defines the distance between the b-spline points relative the the mean linewidth of the peaks."

FineTuneFit::usage = 
"FineTuneFit is an option for FitSpectra and when True it performs a second fitting run where for each peak is an individual linewidth, line shape and shift are fitted."

InitializeFit::usage = 
"InitializeFit is an option for FitSpectra and is used to set initial values for the global fit {gam, eps, {phi0, phi1}, line shape}."

FitLineShape::usage = 
"FitLineShape is an option for FitSpectra and when True allows to fit the line shape. If False a voigt line shape is used."

SpectraOutputPlots::usage = 
"SpectraOutputPlots is an option for FitSpectra. If True the automatic calibration plot for the initial fit are generated."

ReadoutType::usage = 
"ReadoutType is an option for FitSpectra and padding and apodization functions. Value can be \"Fid\" or \"Echo\"."


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
   RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}, 
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
       3 -> "Sagittal", 4 -> "All Planes", 5 -> "Planes 3D", 
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
     {{column, Round[size3/2], "Sagittal"}, 1, Dynamic[size3], 1},
     
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
  size = Round[Norm[#^2 & /@ Drop[Sort[dim], 1]] // N];
  
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


(* ::Subsection::Closed:: *)
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


(* ::Subsection:: *)
(*OLD IVIM functions*)


(* ::Subsection:: *)
(*Bayesian Functions*)


(* ::Subsubsection::Closed:: *)
(*BayesianIVIMFit2*)


Options[BayesianIVIMFit2] = {ChainSteps -> {20000, 1000, 10}, UpdateStep -> {0.5, 0.2, 0.5}, 
	FixPseudoDiff -> False, CorrectPar->True, FixPseudoDiffSD -> 0.5, OutputSamples->False,
	FitConstrains -> ThetaConv[{{-7.6, 7.6}, {-10.0, -5.7}, {-7.0, 0.}}]
	};

SyntaxInformation[BayesianIVIMFit2] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

BayesianIVIMFit2[data_, bval_, fitpari_, maski_, opts : OptionsPattern[]] := Module[{
	useDat, thetai, fix, ynf, fixSD, out1, out2, h1, solution,mask,
	fitpar, deviation, con2, con2e, mui, covi, mmu, mcov,post, dep
	},

	con2 = OptionValue[FitConstrains];
	con2e = ThetaConvi[con2];
	fix = OptionValue[FixPseudoDiff];
	fixSD = OptionValue[FixPseudoDiffSD];

	dep=ArrayDepth[data];

	mask = Mask[Switch[dep,3,Mean[data],4,Mean@Transpose@data], 0.000001]maski;

	fitpar=ThetaConvi[MapThread[N[mask Clip[#1, #2]] &, {fitpari, con2}]];
	fitpar=If[OptionValue[CorrectPar], CorrectParMap[fitpar, con2e, mask], fitpar];

	useDat = MaskData[data, mask];

	{thetai,post} = DataToVector[fitpar,mask];
	thetai=Transpose@thetai;

	ynf = First@DataToVector[Switch[dep,3,data,4,Transpose@data],mask];

	If[fix, 
		thetai[[3]] = RandomVariate[NormalDistribution[Mean[thetai[[3]]], fixSD], Length[thetai[[3]]]]
		];

	{mui, covi} = MeanCov[thetai];

	(*show the pre fit distribution*)
	Print[Dimensions[ynf],Dimensions[thetai]];
	Print[h1 = HistogramPar[thetai, {con2e, 75, mui, covi}, 3, Gray, {0.1, 0.1, 0.1}]];    

	out2 = BayesianIVIMFitI2[thetai, bval, ynf, FilterRules[{opts}, Options[BayesianIVIMFitI2]]];
	solution = out2[[2]];
	deviation = out2[[3]];
	out1 = Chop[VectorToData[Transpose@ThetaConv[solution], post]];

	{mmu, mcov} = MeanCov[solution];
	Print[Column[{
		(*h1,*)
		HistogramPar[solution, {con2e, 75, mmu, mcov}, 3, Blue, {0.1, 0.2, 0.2}],
		UncertainPlot[solution, deviation, con2e, 5 Median[#] & /@ deviation]
		}]
		];

	If[OptionValue[OutputSamples],{out1, out2},out1]
	]

Options[BayesianIVIMFitI2] = {ChainSteps -> {20000, 1000, 10}, UpdateStep -> {0.5, 0.1, 0.5}};

BayesianIVIMFitI2[thetai_, bval_, yn_, OptionsPattern[]] := Block[{
		j, w, w1, w2, w3, wup1, wup2, wup3, yty, t1, t2, mu, cov, theta, nvox, nbval,
		t2s, ttot, t2m, fj, fjt, dj, djt, pdj, pdjt, muj, covj, icovj, gj, gjt, 
		bool1, bool2, bool3, boolf, rU, steps, wstart, nit, burn, sow
	},

	{nit, burn, sow} = OptionValue[ChainSteps];
	steps = nit + burn;
	wstart = OptionValue[UpdateStep];
	nvox = Length[thetai[[1]]];
	nbval = Length[bval];
	ttot={};   

	t1 = First[Timing[
		j = 0;
		(*number of voxels and bvals*)
		Print[nbval, " bvalues x ", nvox, " voxels"];
		(*define yn*)
		yty = Dotc1[yn];
		(*rU := RandomReal[1, nvox];*)

		(*step 2 - initialize mu(j) and cov(j) for j=1 - thetaj={fj,dj,pdj}*)
		{fj, dj, pdj} = thetai;
		{muj, covj} = MeanCov[thetai];
		(*initialize loop pars*)
		gj = Transpose@FunceC2l[fj, dj, pdj, bval];
		(* define Nfr(i), Ndc(i), 
		Npdc(i) needed to update w in first 500 burn steps*)
		{w1, w2, w3} = Transpose[ConstantArray[wstart, {nvox}]];
		wup1 = wup2 = wup3 = ConstantArray[0, {nvox}];

		(*step 3 - further steps of the MCMC j= 2, 3, ... *)
		Monitor[
		{mu, cov, theta, t2s} = Last[Reap[
			Do[t2 = First[Timing[
					j++;

					(*step 3a - Sampel mu(j) [A2] for j=2*)
					(*step 3b - Sampel covu(j) [A3] for j=2*)
					{muj, covj, icovj} = RandomGibsSample[{fj, dj, pdj}, covj, nvox];
										
					(*steps 3c "loop" over the voxels i, perform as vector for each of the parameters *)

					(*step 3c-i - define theta(j) - {fj,dj,pdj}=thetaj;*)

					(*step 3c-ii - ramom sample frtmp*)
					(*comp 1*)
					fjt = RandomNormalCf[fj, w1];
					gjt = Transpose@FunceC2l[fjt, dj, pdj, bval];
					bool1 = Quiet@AlphaC[{fj, dj, pdj}, {fjt, dj, pdj}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool1, gj, gjt];
					fj = BoolAdd[bool1, fj, fjt];

					(*step 3c-iii - ramom sample dctmp*)
					djt = RandomNormalCd[dj, w2];
					gjt = Transpose@FunceC2l[fj, djt, pdj, bval];
					bool2 = Quiet@AlphaC[{fj, dj, pdj}, {fj, djt, pdj}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool2, gj, gjt];
					dj = BoolAdd[bool2, dj, djt];

					(*step 3c-iv - ramom sample pdctmp *)
					pdjt = RandomNormalCd[pdj, w3];
					gjt = Transpose@FunceC2l[fj, dj, pdjt, bval];
					bool3 = Quiet@AlphaC[{fj, dj, pdj}, {fj, dj, pdjt}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool3, gj, gjt];
					pdj = BoolAdd[bool3, pdj, pdjt];

					(*flip and clip if dc>pdc1 *)
					boolf = BooleC[dj, pdj];
					{fj, dj, pdj} = {(-2 boolf + 1) fj, BoolAdd[boolf, dj, pdj], BoolAdd[boolf, pdj, dj]};

					(*check update wup*)
					If[j < 500, 
						wup1 += bool1; wup2 += bool2; wup3 += bool3;
					(*each 100th step reset wup*)
					If[MemberQ[{100, 200, 300, 400, 500}, j],
						w = {w1, w2, w3} = {w1, w2, w3} 101/(2(101 - {wup1, wup2, wup3}));
						wup1 = wup2 = wup3 = ConstantArray[0, {nvox}];
					]];
					]];(*close timing2*)

				(*sow solution every x steps*)
				If[Mod[j - 1, sow] == 0 && j > burn,
				Sow[muj, 1]; 
				Sow[covj, 2]; 
				Sow[{fj, dj, pdj}, 3]; 
				Sow[t2, 4];
				];

				t2m = If[j < 20, Mean[AppendTo[ttot, t2]], Mean[Drop[AppendTo[ttot, t2], 1]]];

				, {steps}];(*close Do loop*)
			]];(*close Reap*)

		(*monitor stuff*)
		, Row[{
			{(steps) - j, NumberForm[t2, {4, 3}], 
			NumberForm[Round[(((steps) - j) t2m)/60, .1], {4, 1}]},
			"     mu:  ", NumberForm[#, {5, 2}] & /@ muj // MatrixForm, 
			",", 
			NumberForm[#, {5, 2}] & /@ ({100, 1000, 1000} ThetaConv[muj]) // MatrixForm,
			"     cov:  ", NumberForm[#, {5, 2}] & /@ # & /@ covj // MatrixForm
			}]
		](*close monitor*)
		]];(*close Timing1*)

	Print[PlotPerformance[{steps, nvox}, {t1, t2s}, {w, wstart}]];

	{
	theta,
	Table[Mean[theta[[All, i]]], {i, Length[theta[[1]]]}],
	Table[StandardDeviation[theta[[All, i]]], {i, Length[theta[[1]]]}],
	cov, mu, Length[theta]}
];


(* ::Subsubsection::Closed:: *)
(*BayesianIVIMFit3*)


Options[BayesianIVIMFit3] = {ChainSteps -> {20000, 1000, 10}, 
	UpdateStep -> {0.5, 0.5, 0.1, 0.5, 0.5}, FixPseudoDiff -> False, 
	CorrectPar->True, OutputSamples->False, FixPseudoDiffSD -> 0.5, 
	FitConstrains ->ThetaConv[{{-7.6, 7.6}, {-7.6, 7.6}, {-10., -5.5}, {-6.5, -2.3}, {-5.2, 0.}}]};

SyntaxInformation[BayesianIVIMFit3] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

BayesianIVIMFit3[data_, bval_, fitpari_, maski_, opts : OptionsPattern[]] := Module[{fitpar,con3,mask,
		useDat, thetai, ynf, fix, fixSD, out1, out2, h1, dep,
		solution, deviation, con3e, mui, covi, mmu, mcov,post
	},

	con3 = OptionValue[FitConstrains];
	con3e = ThetaConvi[con3];
	fix = OptionValue[FixPseudoDiff];
	fixSD = OptionValue[FixPseudoDiffSD];

	dep=ArrayDepth[data];

	mask = Mask[Switch[dep,3,Mean[data],4,Mean@Transpose@data], 0.000001]maski;

	fitpar=ThetaConvi[MapThread[N[mask Clip[#1, #2]] &, {fitpari, con3}]];
	fitpar=If[OptionValue[CorrectPar], CorrectParMap[fitpar, con3e, mask], fitpar];

	useDat = MaskData[data, mask];

	{thetai,post} = DataToVector[fitpar,mask];
	thetai=Transpose@thetai;

	ynf = First@DataToVector[Switch[dep,3,data,4,Transpose@data],mask];

	If[fix,
		thetai[[4]] = RandomVariate[NormalDistribution[Mean[thetai[[4]]], fixSD], Length[thetai[[4]]]];
		thetai[[5]] = RandomVariate[NormalDistribution[Mean[thetai[[5]]], fixSD], Length[thetai[[5]]]];
	];

	{mui, covi} = MeanCov[thetai];

	(*show the pre fit distribution*)
	Print[Dimensions[ynf],Dimensions[thetai]];
	Print[h1 = HistogramPar[thetai, {con3e, 75, mui, covi}, 3, Gray, {0.1, 0.1, 0.1, 0.1, 0.1}]];

	out2 = BayesianIVIMFitI3[thetai, bval, ynf, FilterRules[{opts}, Options[BayesianIVIMFitI3]]];
	solution = out2[[2]];
	deviation = out2[[3]];
	out1 = Chop[VectorToData[Transpose@ThetaConv[solution], post]];

	{mmu, mcov} = MeanCov[solution];
	Print[Column[{
		(*h1,*)
		HistogramPar[solution, {con3e, 75, mmu, mcov}, 3, Blue, {0.1, 0.1, 0.1, 0.1, 0.1}],
		UncertainPlot[solution, deviation, con3e, 5 Median[#] & /@ deviation]
		}]
	];

	If[OptionValue[OutputSamples],{out1, out2},out1]
	]

Options[BayesianIVIMFitI3] = {ChainSteps -> {20000, 1000, 10}, UpdateStep -> {0.5, 0.5, 0.2, 0.5, 0.5}};

BayesianIVIMFitI3[thetai_, bval_, yn_, OptionsPattern[]] := Block[{
		j, w, w1, w2, w3, w4, w5, wup1, wup2, wup3, wup4, wup5, yty, 
		t1, t2, ttot, t2m, mu, cov, muj, covj, icovj, theta, rU, t2s, nvox, nbval,
		f1jc, f2jc, f1j, f2j, f1jt, f2jt, dj, djt, pd1j, pd2j, pd1jt, pd2jt, gj, 
		gjt, bool1, bool2, bool3, bool4, bool5, boolf, steps, wstart, nit, burn, sow
	},

	{nit, burn, sow} = OptionValue[ChainSteps];
	steps = nit + burn;
	wstart = OptionValue[UpdateStep];
	nvox = Length[thetai[[1]]];
	nbval = Length[bval];
	ttot={};

	t1 = First[Timing[
		j = 0;
		(*number of voxels and bvals*)
		Print[nbval, " bvalues x ", nvox, " voxels"];
		(*define yn*)
		yty = Dotc1[yn];
		(*rU := RandomReal[1, nvox];*)

		(*step 2 - initialize mu(j) and cov(j) for j=1 - thetaj={fj,dj,pdj}*)
		{f1j, f2j, dj, pd1j, pd2j} = thetai;
		{muj, covj} = MeanCov[thetai];
		(*initialize loop pars*)
		gj = Transpose@FunceC3l[f1j, f2j, dj, pd1j, pd2j, bval];
		(* define Nfr(i), Ndc(i), 
		Npdc(i) needed to update w in first 500 burn steps*)
		{w1, w2, w3, w4, w5} = Transpose[ConstantArray[wstart, {nvox}]];
		wup1 = wup2 = wup3 = wup4 = wup5 = ConstantArray[0, {nvox}];

		(*step 3 - further steps of the MCMC j= 2, 3, ... *)
		Monitor[
			{mu, cov, theta, t2s} = Last[Reap[
				Do[t2 = First[Timing[
					j++;

					(*step 3a - Sampel mu(j) [A2] for j=2*)
					(*step 3b - Sampel covu(j) [A3] for j=2*)
					{muj, covj, icovj} = RandomGibsSample[{f1j, f2j, dj, pd1j, pd2j}, covj, nvox];

					(*steps 3c "loop" over the voxels i, perform as vector for each of the parameters *)

					(*step 3c-i - define theta(j) - {fj,dj,pdj}=thetaj;*)

					(*step 3c-ii - ramom sample frtmp*)
					(*comp 1*)
					f1jt = RandomNormalCf[f1j, w1];
					gjt = Transpose@FunceC3l[f1jt, f2j, dj, pd1j, pd2j, bval];
					bool1 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1jt, f2j, dj,  pd1j, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool1, gj, gjt];
					f1j = BoolAdd[bool1, f1j, f1jt];

					(*comp 2*)
					f2jt = RandomNormalCf[f2j, w2];
					gjt = Transpose@FunceC3l[f1j, f2jt, dj, pd1j, pd2j, bval];
					bool2 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2jt, dj, pd1j, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool2, gj, gjt];
					f2j = BoolAdd[bool2, f2j, f2jt];

					(*step 3c-iii - ramom sample dctmp*)
					djt = RandomNormalCd[dj, w3];
					gjt = Transpose@FunceC3l[f1j, f2j, djt, pd1j, pd2j, bval];
					bool3 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2j, djt, pd1j, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool3, gj, gjt];
					dj = BoolAdd[bool3, dj, djt];

					(*step 3c-iv - ramom sample pdctmp *)
					(*comp 1*)
					pd1jt = RandomNormalCd[pd1j, w4];
					gjt = Transpose@FunceC3l[f1j, f2j, dj, pd1jt, pd2j, bval];
					bool4 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2j, dj, pd1jt, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool4, gj, gjt];
					pd1j = BoolAdd[bool4, pd1j, pd1jt];

					(*comp 2*)
					pd2jt = RandomNormalCd[pd2j, w5];
					gjt = Transpose@FunceC3l[f1j, f2j, dj, pd1j, pd2jt, bval];
					bool5 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2j, dj, pd1j, pd2jt}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
					gj = BoolAdd[bool5, gj, gjt];
					pd2j = BoolAdd[bool5, pd2j, pd2jt];

					(*flip if dc>pdc1 and if pdc1>pdc2*)
					boolf = BooleC[dj, pd1j];
					f1jc = FConvf[f1j];
					f2jc = FConvf[f2j];
					{f1j, dj, pd1j} = {FConvif[Clip[boolf (1 - 2 f1jc - f2jc) + f1jc, {0.0005, 0.9995}]], BoolAdd[boolf, dj, pd1j], BoolAdd[boolf, pd1j, dj]};
					boolf = BooleC[pd1j, pd2j];
					{f1j, f2j, pd1j, pd2j} = {BoolAdd[boolf, f1j, f2j], BoolAdd[boolf, f2j, f1j], BoolAdd[boolf, pd1j, pd2j], BoolAdd[boolf, pd2j, pd1j]};

					(*check update wup*)
					If[j < 500,
						wup1 += bool1; wup2 += bool2; wup3 += bool3; wup4 += bool4; 
						wup5 += bool5;
						(*each 100th step reset wup*)
						If[MemberQ[{100, 200, 300, 400, 500}, j],
							w = {w1, w2, w3, w4, w5} = {w1, w2, w3, w4, w5} 101/(2 (101 - {wup1, wup2, wup3, wup4, wup5}));
							wup1 = wup2 = wup3 = wup4 = wup5 = ConstantArray[0, {nvox}];
						]
					];
				]];(*close Timing2*)

				(*sow solution every x steps*)
				If[Mod[j - 1, sow] == 0 && j > burn, 
					Sow[muj, 1]; 
					Sow[covj, 2];
					Sow[{f1j, f2j, dj, pd1j, pd2j}, 3]; 
					Sow[t2, 4];
				];

				t2m = If[j < 20, Mean[AppendTo[ttot, t2]], Mean[Drop[AppendTo[ttot, t2], 1]]];

				, {steps}];(*close Do loop*)
			]];(*close Reap*)

			(*monitor stuff*)
			, Row[{
				{steps - j, NumberForm[t2, {4, 3}], 
				NumberForm[Round[(((steps) - j) t2m)/60, .1], {4, 1}]},
				" mu:  ", NumberForm[#, {5, 2}] & /@ muj // MatrixForm, 
				",", 
				NumberForm[Round[#, .01], {5, 2}] & /@ ({100, 100, 1000, 1000, 1000} ThetaConv[muj]) // MatrixForm,
				" cov:  ", NumberForm[#, {5, 2}] & /@ # & /@ covj // MatrixForm
			}]
		](*close monitor*)
	]];(*close Timing*)

	Print[PlotPerformance[{steps, nvox}, {t1, t2s}, {w, wstart}]];

	{
		theta,
		Table[Mean[theta[[All, i]]], {i, Length[theta[[1]]]}],
		Table[StandardDeviation[theta[[All, i]]], {i, Length[theta[[1]]]}],
		cov, mu, Length[theta]
	}
];


(* ::Subsubsection::Closed:: *)
(*Bayesian Core Functions*)


BooleC = Compile[{{val, _Real, 1}, {ru, _Real, 1}}, 
	UnitStep[val - ru]
, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

BooleC1 = Compile[{{val, _Real, 1}, {ru, _Real, 0}}, 
	UnitStep[val - ru]
, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

BooleC2 = Compile[{{val, _Real, 1}, {min, _Real, 0}, {max, _Real, 0}}, 
	UnitStep[val - min] (1. - UnitStep[val - max])
, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

BoolAdd = N[#2 - #1 #2 + #1 #3] &;


ClipC = Compile[{{theta, _Real, 2}, {trans, _Integer, 0}}, Block[{chk1, out},
	If[Length[theta] == 3,
		chk1 = BooleC2[theta[[1]], -7.0, 7.0];
		chk1 = chk1*BooleC2[theta[[2]], -9.5, -5.0];
		(*chk1=chk1*(1-BooleC1[theta[[3]],-0.001]);*)
		chk1 = chk1*BooleC2[theta[[3]], -5.25, -0.001];
		out = DeleteCases[chk1*Transpose[theta], {0., 0., 0.}];
		,
		chk1 = BooleC2[theta[[1]], -7.0, 7.0];
		chk1 = chk1*BooleC2[theta[[2]], -7.0, 7.0];
		chk1 = chk1*BooleC2[theta[[3]], -9.5, -5.5];
		chk1 = chk1*BooleC2[theta[[5]], -7.5, -0.001];
		(*chk1=chk1*(1-BooleC1[theta[[5]],-0.0001]);*)
		out = DeleteCases[chk1*Transpose[theta], {0., 0., 0., 0., 0.}];
	];

	If[trans == 1, Transpose[out], out]
], Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

MeanCov = Block[{inp = ClipC[#, 0]}, {Mean[inp], Covariance[inp]}] &;

(*random Gibs Samplers*)

PosSym[mati_] := Block[{mat = Round[mati, 10.^-20]},
	mat = If[PositiveDefiniteMatrixQ[mat], mat, PosDef[mat]];
	(mat + Transpose[mat])/2
];

PosDef[mat_, tol_: 10.^-5] := Block[{eigsys},
	(*make matrix posdef*)
	NestWhile[(
		eigsys = Eigensystem[#];
		(Eigensystem[#][[2]] . DiagonalMatrix[
			Max[#, tol] & /@ (eigsys[[1]])] . Transpose[eigsys[[2]]])
		) &, N[mat], (! PositiveDefiniteMatrixQ[#] &)]
];

RandomGibsSample[theta_, cov_, m_] := Block[{munew, tm, icov, mat,tmt,mi},
	(*PosSym[cov/m]*)
	munew = N[RandomVariate[MultinormalDistribution[Mean /@ N[theta],N[cov/m]]]];
	(*munew = N[(1 + munew) - 1];*)
	tm = ClipC[theta, 1] - munew;
	tmt=Chop[N[tm . Transpose[tm]], 10^-5];
	(*mi=m-3;*)
	icov = N[RandomVariate[InverseWishartMatrixDistribution[m-3, tmt]]];
	{munew, icov, N@PseudoInverse[icov]}
];

(*random Normal Samplers*)
RandomNormalC = Compile[{{m, _Real, 1}, {s, _Real, 1}},
	Chop[MapThread[RandomVariate[NormalDistribution[#1, #2^2]] &, {m, s}]]
, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];
RandomNormalCf = Compile[{{m, _Real, 1}, {s, _Real, 1}},
	Chop[MapThread[RandomVariate[NormalDistribution[#1, #2^2]] &, {m, s}]]
, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];
RandomNormalCd = Compile[{{m, _Real, 1}, {s, _Real, 1}},
	Chop[Clip[MapThread[RandomVariate[NormalDistribution[#1, #2^2]] &, {m, s}],{-15.,0.4}]]
, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

(*calulated fitted points g(fr, dc, pdc)*)
FunceC2 = Compile[{{fr, _Real, 1}, {dc, _Real, 1}, {pdc, _Real, 1}, {bm, _Real, 1}},Block[{fre=Exp[fr]},
	Chop[Transpose[Map[((Exp[Exp[dc] #] + fre Exp[Exp[pdc] #])/(1 + fre)) &, -bm]]]
], Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

FunceC2l = Compile[{{fr, _Real, 1}, {dc, _Real, 1}, {pdc, _Real, 1}, {bm, _Real, 0}}, 
	Chop[((Exp[-bm Exp[dc]] + Exp[fr] Exp[-bm Exp[pdc]])/(1 + Exp[fr]))]
, Parallelization -> True, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

(*calulated fitted points g(fr1, fr2, dc, pdc1, pdc2)*)
FunceC3 = Compile[{{fr1, _Real, 1}, {fr2, _Real, 1}, {dc, _Real, 1}, {pdc1, _Real, 1}, {pdc2, _Real, 1}, {bm, _Real, 1}},
	Block[{fr1e = Exp[fr1], fr2e = Exp[fr2]},
		Chop[Transpose[Map[((
			(Exp[Exp[pdc1] #] fr1e)/(1 + fr1e) +
			(Exp[Exp[pdc2] #] fr2e)/(1 + fr2e) -
			(Exp[Exp[dc] #] (-1 + fr1e fr2e))/((1 + fr1e) (1 + fr2e))
			)) &, -bm]]
		]
], Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

FunceC3l = Compile[{{fr1, _Real, 1}, {fr2, _Real, 1}, {dc, _Real, 1}, {pdc1, _Real, 1}, {pdc2, _Real, 1}, {bm, _Real, 0}},
	Block[{fr1e = Exp[fr1], fr2e = Exp[fr2]},
		Chop[(
			(Exp[-bm Exp[pdc1]] fr1e)/(1 + fr1e) +
			(Exp[-bm Exp[pdc2]] fr2e)/(1 + fr2e) -
			(Exp[-bm Exp[dc]] (-1 + fr1e fr2e))/((1 + fr1e) (1 + fr2e))
		)]
], Parallelization -> True, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

(*calculate probability*)
DotC = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}}, ((vec1 . vec2)^2)/(vec2 . vec2),
	RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];
Dotc1 = Compile[{{vec, _Real, 1}}, vec . vec,
	RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];
MatDot2 = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}, {mat, _Real, 2}}, (vec1 . mat . vec1) - (vec2 . mat . vec2),
	RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

AlphaC = Compile[{
		{theta, _Real, 2}, {thetat, _Real, 2}, {mu, _Real, 1},
		{icov, _Real, 2}, {y, _Real, 2}, {yty, _Real, 1}, {g, _Real, 2}, 
		{gt, _Real, 2}, {nb, _Real, 0}, {nvox, _Real, 0}
	}, Block[{gttgt, gtg, ytg, ytgt, pt, pd, alpha,pdpt,rand,bool,top,bot},
		(*probability 1*)
		pt = Exp[0.5 (MatDot2[Transpose[theta - mu], Transpose[thetat - mu], icov])];
		(*probability 2*)
		pd = Chop[((yty - DotC[y, gt])/(yty - DotC[y, g]))^(-nb/2)];
		(*bool=alpha-RU*)
		rand = RandomReal[1, nvox];
		UnitStep[(pd*pt) - rand]
],	Parallelization -> True, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


(* ::Subsection:: *)
(*Bayesian Support Functions*)


(* ::Subsubsection::Closed:: *)
(*FracCorrect*)


SyntaxInformation[FracCorrect] = {"ArgumentsPattern" -> {_, _, _.}};

FracCorrect[f1_, time_] := Block[{te, tr, t2t, t21, t1t, t11, st, s1},
	{{te, tr}, {t2t, t21}, {t1t, t11}} = time;
	st = Sigval[{1, t1t, t2t}, tr, te] // N;
	s1 = Sigval[{1, t11, t21}, tr, te] // N;
	((f1*st)/(s1 - f1*s1 + f1*st))
]

FracCorrect[{f1_, f2_?VectorQ}, time_] := Block[{te, tr, t2t, t21, t22, t1t, t11, t12, st, s1, s2},
	{{te, tr}, {t2t, t21, t22}, {t1t, t11, t12}} = time;
	st = Sigval[{1, t1t, t2t}, tr, te] // N;
	s1 = Sigval[{1, t11, t21}, tr, te] // N;
	s2 = Sigval[{1, t12, t22}, tr, te] // N;
	{(f1*s2*st)/(s1*s2 - f1*s1*s2 - f2*s1*s2 + f2*s1*st + f1*s2*st),
	(f2*s1*st)/(s1*s2 - f1*s1*s2 - f2*s1*s2 + f2*s1*st + f1*s2*st)}
]

(*correct fraction for T2 relaxation*)
Sigval[par_, tr_, te_] := par[[1]] (1 - Exp[-tr/par[[2]]]) Exp[-te/par[[3]]]=






(* ::Subsubsection::Closed:: *)
(*ThetaConv*)


SyntaxInformation[ThetaConv] = {"ArgumentsPattern" -> {_}};

ThetaConv[{f1_, dc_, pdc_}] := {Exp[f1]/(1 + Exp[f1]), Exp[dc], Exp[pdc]};
ThetaConv[{f1_, f2_, dc_, pdc1_}] := {Exp[f1]/(1 + Exp[f1]), Exp[f2]/(1 + Exp[f2]), Exp[dc], Exp[pdc1]};
ThetaConv[{f1_, f2_, dc_, pdc1_, pdc2_}] := {Exp[f1]/(1 + Exp[f1]), Exp[f2]/(1 + Exp[f2]), Exp[dc], Exp[pdc1], Exp[pdc2]};


(* ::Subsubsection::Closed:: *)
(*ThetaConvi*)


SyntaxInformation[ThetaConvi] = {"ArgumentsPattern" -> {_}};

ThetaConvi[{f_, dc_, pDc_}] := N[{Log[f] - Log[1 - f], Log[dc], Log[pDc]}] /. {-Infinity -> 0., Indeterminate -> 0.};
ThetaConvi[{f1_, f2_, dc_, pDc1_}] := N[{Log[f1] - Log[1 - f1], Log[f2] - Log[1 - f2], Log[dc], Log[pDc1]}] /. {-Infinity -> 0., Indeterminate -> 0.};
ThetaConvi[{f1_, f2_, dc_, pDc1_, pDc2_}] := N[{Log[f1] - Log[1 - f1], Log[f2] - Log[1 - f2], Log[dc], Log[pDc1], Log[pDc2]}] /. {-Infinity -> 0., Indeterminate -> 0.};


(* ::Subsubsection::Closed:: *)
(*FConvert*)


SyntaxInformation[FConvert] = {"ArgumentsPattern" -> {_}};

FConvert[f_] := If[VectorQ[f], FConvf[f], FConv[f]]
FConv = Compile[{{f1, _Real, 3}}, Exp[f1]/(1 + Exp[f1]), Parallelization -> True];
FConvf = Compile[{{f1, _Real, 1}}, Exp[f1]/(1 + Exp[f1]), Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*FConverti*)


SyntaxInformation[FConverti] = {"ArgumentsPattern" -> {_}};

FConverti[f_] := If[VectorQ[f], FConvif[f], FConvi[f]];
FConvi = Compile[{{f, _Real, 3}}, Log[f] - Log[1 - f], Parallelization -> True];
FConvif = Compile[{{f, _Real, 1}}, Log[f] - Log[1 - f], Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*CorrectParMap*)


SyntaxInformation[CorrectParMap] = {"ArgumentsPattern" -> {_, _, _}};

CorrectParMap[par_, con_, mask_] := Module[{dim, mean, cov, sig, clipmap, rand, clippar, parnew},

	{mean, cov} = MeanCov[Transpose@DataToVector[par,mask][[1]]];
	sig = Diagonal[cov];
	dim = Dimensions[First@par];

	MapThread[(
		clipmap = mask ((1 - Mask[#1, #2[[1]] + .001]) + Mask[#1, #2[[2]] - .001]);
		rand = clipmap RandomVariate[NormalDistribution[#3, #4], dim];
		clippar = (1 - clipmap) #1;
		parnew = rand + clippar
	) &, {par, con, mean, sig}]
]


(* ::Subsubsection::Closed:: *)
(*HistogramPar*)


SyntaxInformation[HistogramPar] = {"ArgumentsPattern" -> {_, _, _, _, _.}};

HistogramPar[dat_, {con_, bin_}, sel_, col_, ran_: .5] :=Block[{mu,cov},
	{mu, cov} = MeanCov[dat]; 
	HistogramPar[dat, {con, bin, mu, cov}, sel, col, ran]
	];

HistogramPar[dat_, {con_, bin_, mu_, cov_}, sel_, col_, ran_: .5] := Module[{label, data, ticks, tickst, tickste, len, ss, hist, pdf, binsize,x},
	If[Length[dat] != Length[con], Return[]];

	binsize = (-Subtract @@ #/bin) & /@ con;

	tickst = {{0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}, {0.1, 0.2, 0.5,
			1, 2, 3}, {0.5, 1, 2, 5, 10, 20, 50, 100, 500}};
	tickste = ThetaConvi[{1., .001, .001} tickst];

	data = DeleteCases[Flatten[#], 0.] & /@ dat;
	(*If[sel == 1, 
		DeleteCases[Flatten[#], 0.] & /@ dat, 
		DeleteCases[Flatten[#], 0.(*Indeterminate*)] & /@ dat
		];*)

	label = {{"f", "D", "pD"}, {"f", "d", "pd"}, {"f  (no units)", 
			"D  (\!\(\*SuperscriptBox[\(10\), \(-3\)]\) \
	\!\(\*SuperscriptBox[\(mm\), \(2\)]\)/s)", 
			"pD  (\!\(\*SuperscriptBox[\(10\), \(-3\)]\) \
	\!\(\*SuperscriptBox[\(mm\), \(2\)]\)/s)"}}[[sel]];
	len = Length[data];

	hist = (
			ss = If[len == 5, {1, 1, 2, 3, 3}[[#]], #];

			ticks = (If[sel == 3, Thread[{tickste[[ss]], tickst[[ss]]}], 
				Automatic]);

			Histogram[
			data[[#]],
			HistRange[con[[#]], bin], "Probability",
			AxesOrigin -> {con[[#, 1]], 0},
			Frame -> {{True, False}, {True, False}}, 
			FrameLabel -> {label[[ss]], "Probability"},
			FrameTicks -> {{Automatic, Automatic}, {ticks, Automatic}}, 
			PlotRange -> {0, ran[[#]]},
			ChartStyle -> Directive[{col, EdgeForm[None], Opacity[0.9]}], 
			LabelStyle -> Directive[{Bold, 12, FontFamily -> "Helvetica"}],
			PerformanceGoal -> "Speed", ImageSize -> 300
			]
			) & /@ Range[len];

	pdf = If[mu === 0 && cov === 0, 
		ConstantArray[Graphics[{}], {len}],
		Plot[PDF[NormalDistribution[mu[[#]], Sqrt[cov[[#, #]]]], x]* binsize[[#]], {x, con[[#, 1]], con[[#, 2]]}, PlotStyle -> {Thick, Red}, PlotRange -> Full] & /@ Range[len]
		];

	GraphicsRow[MapThread[Show[#1, #2] &, {hist, pdf}], 
		ImageSize -> len*300, Spacings -> 0]
	]

HistRange[rri_, n_] := Module[{rr}, rr = If[rri[[2]] > 0, {1, 1} rri, {1, 1} rri]; {rr[[1]], rr[[2]], (rr[[2]] - rr[[1]])/n}];


(* ::Subsubsection::Closed:: *)
(*PlotPerformance*)


SyntaxInformation[PlotPerformance] = {"ArgumentsPattern" -> {_, _, _}};

PlotPerformance[{nit_, nvox_}, {t1_, t2s_}, {w_, wstart_}] := Column[{
	Row[{nit, " steps: ", Round[t1/60, .1], 
		" min - each step takes: ", Round[t1/nit, .001], 
		" s - full chain (21000) takes: ", 
		Round[(21000/nit) (t1/60), .1], " min"}]
	,
	(*GraphicsRow[MapThread[Show[
			ListPlot[#1, AspectRatio -> .1 Length[wstart], 
			PlotRange -> All, PlotStyle -> {Black}],
			Plot[#2, {x, 0, nvox}, PlotStyle -> {Red, Thick}]
			] &, {w, wstart}], ImageSize -> 1000]
	,*)
		{w, wstart,nvox};
	Show[
		ListPlot[t2s, AspectRatio -> 0.075, ImageSize -> 1000, 
		PlotStyle -> {Black, PointSize[Medium]}],
		ListLinePlot[{
		GaussianFilter[t2s, 20, Padding -> "Reversed"],
		GaussianFilter[t2s, Length[t2s], Padding -> "Reversed"]
		}, AspectRatio -> 0.075, ImageSize -> 1000, 
		PlotStyle -> {Directive[{Red, Thick}], 
			Directive[{Red, Dashed, Thick}]}, PlotRange -> Full
		]
	]
}]


(* ::Subsubsection::Closed:: *)
(*UncertainPlot*)


SyntaxInformation[UncertainPlot] = {"ArgumentsPattern" -> {_, _, _, _.}};

UncertainPlot[mn_, sig_, con_, ran_:.1] := Module[{tickst, tickste, label, ticks, len, ss},
	tickst = {{0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}, {0.1, 0.2, 0.5, 1, 2, 3}, {0.5, 1, 2, 5, 10, 20, 50, 100, 500}};
	tickste = ThetaConvi[{1., .001, .001} tickst];

	label = {{"f  (no units)", 
		"\!\(\*SubscriptBox[\(\[Sigma]\), \(f\)]\)"}, {"D  \
	(\!\(\*SuperscriptBox[\(10\), \(-3\)]\) \!\(\*SuperscriptBox[\(mm\), \
	\(2\)]\)/s)", 
		"\!\(\*SubscriptBox[\(\[Sigma]\), \(d\)]\)"}, {"pD  \
	(\!\(\*SuperscriptBox[\(10\), \(-3\)]\) \!\(\*SuperscriptBox[\(mm\), \
	\(2\)]\)/s)", "\!\(\*SubscriptBox[\(\[Sigma]\), \(pd\)]\)"}};
	len = Length[mn];

	GraphicsRow[(
		ss = If[len == 5, {1, 1, 2, 3, 3}[[#]], #];

		ticks = Thread[{tickste[[ss]], tickst[[ss]]}];

		ListPlot[{mn[[#]], sig[[#]]} // Transpose,
			Frame -> {{True, False}, {True, False}}, 
			FrameLabel -> label[[ss]], PlotStyle -> Red,
			FrameTicks -> {{Automatic, Automatic}, {ticks, Automatic}}, 
			Axes -> False,
			PlotRange -> {con[[#]], {0, ran[[#]]}}, ImageSize -> 300,
			LabelStyle -> Directive[{Bold, 12, FontFamily -> "Helvetica"}]
		]
		) & /@ Range[len]
	, ImageSize -> len*300, Spacings -> 0]
]


(* ::Subsection:: *)
(*Spectra Fitting Functions*)


(* ::Subsubsection::Closed:: *)
(*FitSpectra*)


Options[FitSpectra]={
	SpectraNucleus->"1H",
	SpectraPpmShift->4.65,
	SpectraFieldStrength->3,
	PaddingFactor->2,
	SplineSpacingFactor->1.5,
	FineTuneFit->True,
	InitializeFit->Automatic,
	FitLineShape->False,
	SpectraOutputPlots->False,
	ReadoutType->"Fid"
};

SyntaxInformation[FitSpectra] = {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}}

FitSpectra[specBasisIn_, specIn_, {st_,end_}, dTime_, lwVals_?VectorQ, opts : OptionsPattern[]]:=FitSpectra[specBasisIn,specIn,{st,end},dTime,{lwVals, 0 lwVals + 1.},opts]

FitSpectra[specBasisIn_, specIn_, {st_,end_}, dTime_, {lwVals_?VectorQ, lwAmps_?VectorQ}, OptionsPattern[]]:=Block[{
	tTotal,log,pad,spFac,field,nuc,shift,plots,init,scale,nBasis,len,
	timeBasis,specFull,timeFull,ppmFull,nSamp,gyro,indSt,indEnd,
	gamI,epsI,phi0i,phi1i,lineI,phiI,plLine,plShift, readout, varF,
	splineSpace,cpn,var,phi0f,phi1f,gamF,epsF,phiF,lineF,gam,eps,line,phi,sigI,
	tFit1,fit1,sol,output,tFit2,fit2,fit,timeBasisIn,time,ppm,spline,basis,error,errors,specFit
	},

	(*time the fitting*)
	tTotal=AbsoluteTiming[
		(*turn of error messages of FinMinimum*)
		Off[FindMinimum::cvmit];Off[FindMinimum::lstol];

		(*logging*)
		log={};(*Print[Dynamic[Column[log]]];*)

		(*get options*)
		pad = OptionValue[PaddingFactor];
		spFac = OptionValue[SplineSpacingFactor];
		field = OptionValue[SpectraFieldStrength];
		nuc = OptionValue[SpectraNucleus];
		shift = OptionValue[SpectraPpmShift];
		plots = OptionValue[SpectraOutputPlots];
		init = OptionValue[InitializeFit];
		readout = OptionValue[ReadoutType];

		(*set general parameters*)
		scale = 1000/Max[Abs[specIn]];
		nBasis = Length[specBasisIn];

		
		(*-------------------------------------------------------------------*)
		(*create the basis functions fid and pad*)
		timeBasis = CashBasisTime[specBasisIn, pad, readout];

		(*pad and normalize the spectra spectra*)
		specFull = scale ApodizePadSpectra[specIn, PaddingFactor->pad, ReadoutType->readout];

		(*get the time and ppm axes*)
		{timeFull,ppmFull}=GetTimePpmRange[specFull,dTime,field,nuc];
		gyro = GetGyro[nuc,field];
		nSamp = Length[specFull];

		(*find the positions of the fit range*)
		{indSt,indEnd}=Flatten[Position[ppmFull,Nearest[ppmFull,#][[1]]]&/@{st,end}];

		(*logging of input parameters*)
		AppendTo[log,Style["Spectral Properties",Bold]];
		AppendTo[log,"    - Number of samples:          "<>ToString[nSamp]];
		AppendTo[log,"    - Gyro magnetic ratio:        "<>ToString[gyro]<>" Hz"];
		AppendTo[log,"    - Number of basis functions:  "<>ToString[nBasis]];
		AppendTo[log,"    - Number of fit samples:      "<>ToString[indEnd-indSt+1]];
		AppendTo[log,""];

		
		(*-------------------------------------------------------------------*)
		If[init=!=Automatic,
			(*in not automatic it is an initial fit result*)
			{gamI,epsI,{phi0i,phi1i},lineI} = init;
			phiI={phi0i,phi1i};
			,
			(*find initial linewidth and spec shift*)
			{gamI,epsI,plLine}=EstimateLineWidth[{ppmFull,specFull},{lwVals,lwAmps},gyro,{st,end},plots];
			(*find initial phase estimate*)
			{{phi0i,phi1i},plShift}=EstimatePhaseShift[{ppmFull,specFull},{timeFull,timeBasis},{gamI,epsI},gyro,{indSt,indEnd},readout,plots];
			(*define the initial line shape*)
			lineI=0.5;

			(*logging of parameters*)
			AppendTo[log,Style["Estimating linewidth",Bold]];
			AppendTo[log,"    - spectral linewidth:         "<>ToString[Round[gamI / gyro,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[gamI,.0001]]<>" Hz"];
			AppendTo[log,"    - base spectra shift:         "<>ToString[Round[epsI,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[epsI gyro,.0001]]<>" Hz"];
			AppendTo[log,""];
			AppendTo[log,Style["Estimating phase correction",Bold]];
			AppendTo[log,"    - zeroth order phase:         "<>ToString[Round[phi0i,.0001]]<>" rad"];
			AppendTo[log,"                                  "<>ToString[Round[phi0i/Degree,.0001]]<>" deg"];
			AppendTo[log,"    - first order phase:          "<>ToString[Round[2Pi phi1i,.0001]]<>" rad/kHz"];
			AppendTo[log,"                                  "<>ToString[Round[phi1i,.0001]]<>" ms"];
			AppendTo[log,""];
		];

		
		(*-------------------------------------------------------------------*)
		(*get the spline basis function parameters*)
		splineSpace=spFac Mean[Flatten[{gamI}]]/gyro;
		cpn =Clip[ Round[Subtract@@Reverse[Sort[{st,end}]]/splineSpace],{4,Round[nSamp/10]}];

		(*logging of parameters*)
		AppendTo[log,Style["Estimating spline smoothness",Bold]];
		AppendTo[log,"    - spline spacing:             "<>ToString[splineSpace]<>" ppm"];
		AppendTo[log,"    - spline control points:      "<>ToString[cpn]];
		AppendTo[log,""];

		
		(*-------------------------------------------------------------------*)
		(*perform the first run minimization*)
		tFit1=0;
		If[NumberQ[gamI]&&NumberQ[epsI]&&NumberQ[lineI],
			(*define fit parameters and initialization*)
			Clear[phi0f,phi1f,gamF,epsF,lineF,var,varF];

			(*define the fit variables*)
			var={{gamF,gamI},{epsF,epsI},{lineF,lineI},{phi0f,phi0i},{phi1f,phi1i}};
			varF={gamF,epsF,{phi0f,phi1f}, lineF};
			(*get the std of initial fit*)
			init={gamI, epsI, phi1i};

			(*perform the fit*)
			{tFit1,fit1}=AbsoluteTiming[FindMinimum[FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},varF, init, Output->"Error", ReadoutType ->readout], var][[2]]];
							
			(*Get the fit results and output, wrap phi between -pi and pi*)
			sol={gamI,epsI,phiI,lineI}={Clip[gamF,{1,500}],epsF,{2ArcTan[Tan[phi0f/2]],phi1f},lineF}/.fit1;

			(*logging of parameters*)
			AppendTo[log,Style["Performing spectra First run",Bold]];
			AppendTo[log,"    - line shape:                 "<>ToString[Round[lineI,.0001]]];
			AppendTo[log,"    - spectral linewidth:         "<>ToString[Round[gamI/gyro,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[gamI,.0001]]<>" Hz"];
			AppendTo[log,"    - base spectra shift:         "<>ToString[Round[epsI,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[epsI gyro,.0001]]<>" Hz"];
			AppendTo[log,"    - zeroth order phase:         "<>ToString[Round[phiI[[1]],.0001]]<>" rad"];
			AppendTo[log,"                                  "<>ToString[Round[phiI[[1]]/Degree,.0001]]<>" deg"];
			AppendTo[log,"    - first order phase:          "<>ToString[Round[2Pi phiI[[2]],.0001]]<>" rad/kHz"];
			AppendTo[log,"                                  "<>ToString[Round[ phiI[[2]],.0001]]<>" ms"];
			AppendTo[log,""];
		];

		(*redefine the spline spacing *)
		splineSpace = spFac Mean[Flatten[{gamI}]]/gyro;
		cpn = Clip[Round[Subtract@@Reverse[MinMax[ppmFull]]/splineSpace],{4,Round[nSamp/10]}];

		(*make the output*)
		{fit,sigI}=FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},sol, Output->"Fit", ReadoutType ->readout];

		
		(*-------------------------------------------------------------------*)
		(*perform the second run minimization*)
		tFit2=0;
		If[OptionValue[FineTuneFit],
			(*prepare the fit parameters*)
			Clear[phi0f,phi1f,phi,gam,eps,line,var,varF,init];

			gamF=Table[Unique[gam],{i,1,nBasis}];
			epsF=Table[Unique[eps],{i,1,nBasis}];
			lineF=If[OptionValue[FitLineShape],Table[Unique[line],{i,1,nBasis}],Unique[line]];

			(*define the fit variables*)
			var=Join[MakeVars[gamF,gamI,1],MakeVars[epsF,epsI,1],MakeVars[{phi0f,phi1f},phiI,0],MakeVars[lineF,lineI,1]];
			varF={gamF,epsF,{phi0f,phi1f}, lineF};
			(*get the std of initial fit*)
			init={gamI,epsI,phiI[[2]]};

			(*perform the minimization*)
			{tFit2,fit2}=AbsoluteTiming[FindMinimum[FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},varF, init, Output->"Error", ReadoutType ->readout], var][[2]]];

			(*get the solution and output, wrap phi between -pi and pi*)
			sol={gamI, epsI, phiI, lineI}={Clip[gamF,{1,500}],epsF,{2ArcTan[Tan[phi0f/2]],phi1f},lineF}/.fit2;

			(*recalculate the spline spacings*)
			splineSpace = spFac Mean[Flatten[{gamI}]]/gyro;
			cpn = Clip[Round[Subtract@@Reverse[MinMax[ppmFull]]/splineSpace],{4,Round[nSamp/10]}];

			(*generate the output*)
			{fit,sigI}=FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},sol, Output->"Fit", ReadoutType ->readout];

			(*logging of parameters*)
			AppendTo[log,Style["Performing spectra: Second run",Bold]];
			AppendTo[log,"    - line shape:                 "<>ToString[Round[Mean[Flatten[{lineI}]],.0001]]];
			AppendTo[log,"                                  "<>ToString[OptionValue[FitLineShape]]];
			AppendTo[log,"    - mean spectral linewidth:    "<>ToString[Round[Mean[gamI]/gyro,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[Mean[gamI],.0001]]<>" Hz"];
			AppendTo[log,"    - mean base spectra shift:    "<>ToString[Round[Mean[epsI],.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[Mean[epsI] gyro,.0001]]<>" Hz"];
			AppendTo[log,"    - zeroth order phase:         "<>ToString[Round[phiI[[1]],.0001]]<>" rad"];
			AppendTo[log,"                                  "<>ToString[Round[phiI[[1]]/Degree,.0001]]<>" deg"];
			AppendTo[log,"    - first order phase:          "<>ToString[Round[2Pi phiI[[2]],.0001]]<>" rad/kHz"];
			AppendTo[log,"                                  "<>ToString[Round[phiI[[2]],.0001]]<>" ms"];
			AppendTo[log,""];
		]
	(*close timing*)
	][[1]];

	
	(*-------------------------------------------------------------------*)
	(*apply the results on the raw input*)
	{time,ppm}=GetTimePpmRange[specIn,dTime,field,nuc];
	gyro=GetGyro[nuc,field];
	nSamp=Length[specFull];

	timeBasisIn = ShiftedInverseFourier[#, readout]&/@specBasisIn;
	basis = BasisSpectraApply[{ppm, time, timeBasisIn}, sol, gyro, readout];
	fit = fit/scale;
	specFit = fit . basis;

	(*fit a spline through the residuals*)
	spline=BSplineCurveFit[specIn-specFit, SplineKnotsNumber-> cpn, SplineRegularization->0, SplineDegree-> 2];

	(*calculate the error*)
	error=specIn-specFit-spline;
	errors=100error/Max[Abs[specIn]];
	errors=ToString[Round[Abs[Mean[errors]],.01]]<>" \[PlusMinus] "<>ToString[Round[StandardDeviation[errors],.01]];

	(*logging*)
	AppendTo[log,Style["Fit performance",Bold]];
	AppendTo[log,"    - Total computation time:     "<>ToString[tTotal]<>" s"];
	AppendTo[log,"         - fit1 time:             "<>ToString[tFit1]<>" s"];
	AppendTo[log,"         - fit2 time:             "<>ToString[tFit2]<>" s"];
	AppendTo[log,"    - Residual error (mn \[PlusMinus] std):  "<>errors<>" %"];

	(*turn on error messages of FinMinimum*)
	On[FindMinimum::cvmit];On[FindMinimum::lstol];

	
	(*-------------------------------------------------------------------*)
	(*give the output,scale data to original values*)
	{Prepend[fit,1],Prepend[basis,spline],error,sol,log,{plLine,plShift}}
]


MakeVars[par_,val_,0]:=If[Length[par]===0,
	{{par,val}},
	If[Length[par]==Length[val],Transpose[{par,val}],Transpose[{par,ConstantArray[val,Length[par]]}]]
]

MakeVars[par_,val_,1]:=If[Length[par]===0,
	{{par,val}},
	If[Length[par]==Length[val],
		Transpose[{par,val}],
		Transpose[{par,RandomReal[{0.97,1.03},Length[par]]ConstantArray[val,Length[par]]}]
	]
]

CashBasisTime[specBasisIn_, pad_, readout_] := CashBasisTime[specBasisIn, pad, readout] = Block[{func},
	func = Switch[readout,"Fid",ApodizePadFid,"Echo",ApodizePadEcho];
	func[ShiftedInverseFourier[#, readout], PaddingFactor->pad]&/@specBasisIn
]


(* ::Subsubsection::Closed:: *)
(*FitSpectraError*)


Options[FitSpectraError] = {Output -> "Error", ReadoutType -> "Fid"};

FitSpectraError[{ppmFull_, spec_}, {timeFull_, timeBasis_}, {indSt_, indEnd_}, {cpn_, gyro_},
	{gam_ /; AllTrue[gam, NumericQ], eps_ /; AllTrue[eps, NumericQ], phi_ /; AllTrue[phi, NumericQ], f_ /; AllTrue[f, NumericQ]}, 
	init___ : 0, OptionsPattern[]] := Block[{
		specF, fidF, specBasisF, fit, fErr, gErr, pErr, errorF, errorS, err, gamI, epsI,
		sigI, phiI, specBasis, fidBasisF, specFit, spline, readout, eErr
	},

	(*get the readout type*)
	readout = OptionValue[ReadoutType];

	(*give output either error or fit results*)
	Switch[OptionValue[Output],
		(* ------------ Calculate Error for minimization --------------*)
		"Error",

		(*----------- apply parameters -----------*)
		(*apply phase to target instead of basis functions (faster)*)
		specF = PhaseShiftSpectra[spec, ppmFull, gyro, -phi];
		(*also convert to fid for error calculations*)
		fidF= ShiftedInverseFourier[specF, readout];
		(*select spectra range*)
		specF = specF[[indSt ;; indEnd]];

		(*generate basis spectra from time domain by applying gam, eps and line shape*)
		{fidBasisF,specBasisF} = BasisSpectraApply[{timeFull, timeBasis}, {gam, eps, f}, gyro, {indSt, indEnd}, readout];

		(*----------- Perform Fit and calculate error -------------*)
		(*perform Fit of basis spectra*)
		fit = Quiet@NNLeastSquares[Join[Transpose[Re[specBasisF]],Transpose[Im[specBasisF]]], Join[Re[specF],Im[specF]]];			

		(*define errors fid and spectra*)
		errorS = specF - fit . specBasisF;
		errorF = fidF - fit . fidBasisF;

		(*Re and Im error normalized for number of points*)
		err = Mean[Join[Re[errorS]^2, Re[errorF]^2, Im[errorS]^2, Im[errorF]^2]];

		If[init === 0,
			(*constrain f between 0 and 1*)
			fErr = Total[ConFuncC[f, 0, 1, 4]];
			(*constrain phase to be small*)
			pErr = ConFuncC[phi[[2]], -0.5, 0.5, 5];
			(*constrain gam to be positive*)
			gErr = Total[ConFuncC[gam, 1, 500, 3]];
			(*no initial values, minimize RMSE with f, gam and phase constraint*)
			err + fErr + gErr + pErr
			,

			(*get initial values to constrain fine tune fit, init = {gamI, epsI, phiI}*)
			(*constrain f between 0 and 1*)
			fErr = Total[ConFuncC[f, 0, 1, 4]];
			(*constrain lw gam to initial value*)
			gErr = Total[ConFuncSC[gam,init[[1]], 10, 2]];
			(*constrain shift eps to initial value*)
			eErr = Total[ConFuncSC[eps,init[[2]], 1, 2]];
			(*constrain phase to initial value*)
			pErr = ConFuncC[phi[[2]]-init[[3]], -0.1, 0.1, 5];

			(*calculate error, minimize RMSE with f, gam, eps and phase constraint*)
			err + fErr + gErr + eErr + pErr
		],

		(* ------------ Calculate fitted spectra for output --------------*)
		"Fit",

		(*generate basis spectra from time domain by applying gam, eps and line shape*)
		specBasis = BasisSpectraApply[{ppmFull, timeFull, timeBasis}, {gam, eps, phi, f}, gyro, readout];

		(*perform Fit of basis spectra*)
		fit = Quiet@Clip[NNLeastSquares[Join[Transpose[Re[specBasis]],Transpose[Im[specBasis]]], Join[Re[spec],Im[spec]]],{0,Infinity}];

		specFit = fit . specBasis;

		(*fit a spline through the residuals*)
		spline = BSplineCurveFit[spec - specFit, SplineKnotsNumber -> cpn, SplineRegularization -> 0, SplineDegree -> 2];

		(*recalculate the error*)
		error = spec - specFit - spline;
		{fit, StandardDeviation[error]^2}
	]
]

ConFuncSC = Compile[{{par, _Real, 0}, {ref, _Real, 0}, {nor, _Real, 0}, {sc, _Real, 0}}, 
	10^sc ((par-ref)/nor)^2
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]

ConFuncUC = Compile[{{par, _Real, 0}, {min, _Real, 0}, {sc, _Real, 0}}, Block[{off = -par + min},
	10^sc (UnitStep[off] (off))^2
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]

ConFuncLC = Compile[{{par, _Real, 0}, {max, _Real, 0}, {sc, _Real, 0}}, Block[{off = par - max},
	10^sc (UnitStep[off] (off))^2
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]

ConFuncC = Compile[{{par, _Real, 0}, {min, _Real, 0}, {max, _Real, 0}, {sc, _Real, 0}}, Block[{off1 = -par + min, off2 = par - max},
	10^sc (UnitStep[off1] (off1)^2 + UnitStep[off2] (off2)^2)
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]


(* ::Subsubsection::Closed:: *)
(*FindSpectraPpmShift*)


FindSpectraPpmShift[spec_, {dw_, gyro_}, peaks_] := FindSpectraPpmShift[spec, {dw, gyro}, {peaks, 0}]

FindSpectraPpmShift[spec_, {dw_, gyro_}, {peaks_, amp_}] := Block[{ppm, dPpm, tar, amps, weight, s, corrF, maxPos, sft},
	ppm = GetPpmRange[spec, dw, gyro];
	dPpm = (ppm[[1]] - ppm[[2]]);

	tar = If[Length[spec] === Length[peaks],
		Abs[peaks],
		amps = If[amp === 0, 0. peaks + 1, amp];
		tar = 0 ppm;
		tar[[Flatten[Position[ppm, First@Nearest[ppm, #]] & /@ peaks]]] = amps/Max[amps];
		tar
	];

	(*perform correlation of spectra with delta function*)
	s = Max[ppm]/4;
	weight = 1/(Exp[(ppm^2/(2*s^2))]*(Sqrt[2*Pi]*s));
	weight = 1;
	corrF = ListCorrelate[tar, Abs[spec], Round[Length[tar]/2], 0];
	corrF = weight DivideNoZero[corrF, Max[corrF]];
	maxPos = Position[corrF, Max[corrF]][[1, 1]];

	(*constrain shift to 5ppm*)
	sft = -dPpm (maxPos - (Length[ppm]/2.));
	If[-5 < sft < 5, sft, 0]
]


(* ::Subsubsection::Closed:: *)
(*BasisSpectraApply*)


BasisSpectraApply[{timeFull_,timeBasis_},{gam_,eps_,f_},gyro_,readout_]:=BasisSpectraApply[{0,timeFull,timeBasis},{gam,eps,0,f},gyro, {1,-1},readout]

BasisSpectraApply[{timeFull_,timeBasis_},{gam_,eps_,f_},gyro_,{st_,end_},readout_]:=BasisSpectraApply[{0,timeFull,timeBasis},{gam,eps,0,f},gyro,{st,end},readout]

BasisSpectraApply[{ppmFull_,timeFull_,timeBasis_},{gam_,eps_,phi_,f_},gyro_,readout_]:=BasisSpectraApply[{ppmFull,timeFull,timeBasis},{gam,eps,phi,f},gyro,{1,-1},readout]

BasisSpectraApply[{ppmFull_,timeFull_,timeBasis_},{gam_,eps_,phi_,f_},gyro_,{st_,end_},readout_]:=Block[{fidBasis, specBasis, vec, gamV, epsV, fv,func},
	(*make all parameters vectors*)
	{gamV, epsV, fv} = ConstantArray[1, {3,Length[timeBasis]}] {gam, eps, f};

	(*get the readout function*)
	func = Switch[readout,"Fid",TimeShiftFid,"Echo",TimeShiftEcho];
	(*generate basis spectra from time domain by applying gam,eps and line shape*)
	fidBasis = MapThread[func[#1, timeFull, gyro, {#2, #3, #4}] &, {timeBasis, gamV, epsV, fv}];
	specBasis = ShiftedFourier[#, readout]&/@fidBasis;

	(*apply phase to the basis spectra, phi is only 0 if used in fitting and there also the fid is needed in all other cases only spectra is needed*)
	If[phi === 0, {fidBasis, specBasis[[All, st ;; end]]}, PhaseShiftSpectra[#, ppmFull, gyro, phi] & /@ specBasis]
]


(* ::Subsubsection::Closed:: *)
(*Estimate Line width*)


(*Function to estimate linewidth*)
EstimateLineWidth[{ppm_,spec_},{peaks_,amps_},gyro_,ran_,plot_:True]:=Block[{
	dPpm,deltaF,corrF,weight, s, max,pts,pos,ppmC,maxF,block,line,sol,x,ppmCf,pl,lw,sft},

	(*define delta ppm and the delta function for correlation*)
	dPpm=(ppm[[1]]-ppm[[2]]);
	deltaF=0 ppm;
	deltaF[[Flatten[Position[ppm,First@Nearest[ppm,#]]&/@peaks]]]=amps/Max[amps];

	(*perform correlation of spectra with delta function*)
	corrF=ListCorrelate[deltaF,Abs[spec],Round[Length[deltaF]/2],0];
	corrF=DivideNoZero[Length[corrF]corrF,Max[corrF]];

	(*Find max correlation and position*)
	s = Max[ppm]/4;
	weight = 1/(Exp[(ppm^2/(2*s^2))]*(Sqrt[2*Pi]*s));
	max=Max[weight corrF];
	pos=Position[weight corrF,max][[1,1]];
	max=corrF[[pos]];

	(*constrain shift to 5ppm*)
	sft=-dPpm(pos-(Length[ppm]/2.));
	sft=If[-5<sft<5,sft,0];
	pos=If[-5<sft<5,pos,0];

	(*find two points closest to FWHM*)
	maxF=0.5max;
	block=If[#<maxF,1,0]&/@corrF;
	pts={-FirstPosition[Reverse@block[[;;pos]],1][[1]]+pos+2,FirstPosition[block[[pos;;]],1][[1]]+pos-1};

	(*prevent not found errors in pts*)
	If[AllTrue[pts,NumericQ],
		line=Transpose@{{-1,0,1}+#,corrF[[{-1,0,1}+#]]}&/@pts;
		sol=Flatten[x/.Solve[Fit[#,{1,x},x]==maxF]&/@line];

		(*constrain the lw*)
		lw=gyro dPpm(Subtract@@(Reverse@Sort@sol));
		lw=If[1<lw<200,lw,50];
		sol=If[1<lw<200,{1,2},sol];
		,
		{lw,sol}={50,{1,2}};
	];
	(*debugging plots*)
	pl=If[plot,
		ppmC=dPpm(Range[Length[corrF]]-Length[corrF]/2);
		ppmCf=ListInterpolation[ppmC];
		FlipView[{Show[
				PlotSpectra[ppm,Max[Abs[spec]]deltaF,Method->"Abs",GridLineSpacing->5,PlotRange->{ran,Full},PlotColor->Red, PlotLabel->"Raw signal and calibration metabolites"],
				PlotSpectra[ppm,spec,Method->"Abs",GridLineSpacing->10,PlotRange->{ran,{0,Max[Abs[spec]]}}]
			],Show[
				PlotSpectra[ppmC,corrF,GridLineSpacing->5,PlotRange->{ran,Full}, PlotLabel->"Convolution signal"],
				ListLinePlot[{Transpose[{ppmCf[sol],{maxF,maxF}}],{{ppmC[[pos]],0},{ppmC[[pos]],max}}},PlotStyle->Directive[{Thick,Red}],ScalingFunctions->{"Reverse",Automatic}]
		]}], Null];

	(*calculate the estimated lw and shift*)
	{lw,sft,pl}
]


(* ::Subsubsection::Closed:: *)
(*EstimatePhaseShift*)


(*function to estimate phase form abs fitting of basis spectra*)
EstimatePhaseShift[{ppm_,spec_},{time_,fids_},{gam_,eps_},gyro_,{st_,en_},readout_,plot_:True]:=Block[{
	phi1,sol1,phi2,sol2,specsC,fit,phi0f,phi1f,phi,specC,ran,pl,lim,specF,ppmF,func
	},

	specF=spec[[st;;en]];
	ppmF=ppm[[st;;en]];
	lim=.1;

	(*convert basis fids in spectra find function based on fid or echo*)
	func=Switch[readout,"Fid",TimeShiftFid,"Echo",TimeShiftEcho];
	specsC=Transpose[ShiftedFourier[func[#,time,gyro,{gam,eps,.5}], readout][[st;;en]]&/@fids];

	(*Fit absolute basis spectra to absolute spectrum*)
	fit=specsC . (NNLeastSquares[Abs[specsC],Abs[specF]]);
	(*minimize error with the target spectra*)
	sol1=Quiet@NMinimize[{PhaseError[ppmF,fit,specF,{phi0f,0 },gyro],-Pi<phi0f<Pi},{phi0f},MaxIterations->25][[2]];
	phi1={phi0f, 0}/.sol1;

	(*apply the zeroth order phase to the basis spectra*)
	specsC=Transpose[PhaseShiftSpectra[#,ppmF,gyro,phi1]&/@Transpose[specsC]];
	(*calculate the fit based on the imaginary part of the spectra*)
	fit=specsC . (NNLeastSquares[Re@specsC,Re@specF]);

	(*minimize error with the target spectra*)
	sol2=Quiet@NMinimize[{PhaseError[ppmF,fit,specF,{phi0f,phi1f },gyro],-Pi<phi0f<Pi,-lim<phi1f<lim},{phi0f,phi1f},MaxIterations->25][[2]];
	phi2={phi0f, phi1f}/.sol2;
	phi=phi1+phi2;

	(*debugging plots*)
	pl = If[plot,
		specC=PhaseShiftSpectra[fit,ppmF,gyro,phi2];
		ran={-1,1}Max[Abs[{specF,specC}]];
		FlipView[{
			Show[PlotSpectra[ppmF,specC,GridLineSpacing->5,PlotRange->{Full,ran}], PlotLabel->"Calibrated fit signal"],
			Show[PlotSpectra[ppmF,specF,GridLineSpacing->5,PlotRange->{Full,ran}], PlotLabel->"Raw signal"]
		}]
		,
		Null
	];

	{phi,pl}
];


PhaseError[ppm_,specI_,specT_,{phi0_?NumericQ,phi1_?NumericQ},gyro_]:=PhaseErrorC[ppm,specI,specT,phi0,phi1,gyro]

PhaseErrorC=Compile[{{ppm,_Real,1},{specI,_Complex,1},{specT,_Complex,1},{phi0,_Real,0},{phi1,_Real,0},{gyro,_Real,0}},Block[{spec},
	spec=Exp[-I(phi0+2Pi (phi1/1000) gyro ppm)]specI;
	Total[(Re[specT]-Re[spec])^2] + Total[(Im[specT]-Im[spec])^2]
],RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


(* ::Subsection::Closed:: *)
(*FitSpectraResultTable*)


SyntaxInformation[FitSpectraResultTable] = {"ArgumentsPattern" -> {_, _, _, _, _.}}

FitSpectraResultTable[parFit_, parsF_, names_, ref_, out_:"tab"] := Block[{
	par, phi, amp, lw, ls, shift, sc, rowName, colName,tabDat, tab, dat
	},

	par = parFit[[2 ;;]];
	phi = {
		{"", "", "", ""},
		Flatten@Thread[{Style[#, Bold,Black] & /@ {"\!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\) [deg]", "\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\) [ms]"},
			Round[{parsF[[3, 1]]/Degree, parsF[[3, 2]]}, .001]}
		]
	};
	sc = If[ref =!= "" && MemberQ[names, ref], par[[Position[names, ref][[1, 1]]]] /. {0. -> Max[par]}, Max[par]];
	amp = If[sc === 1, par, 100 par/sc];

	{lw, ls, shift}=ConstantArray[1,{3,Length[names]}]parsF[[{1, 4, 2}]];

	rowName = Join[names, {"", "phase"}];
	colName = {"Amp.", "LW [Hz]", "shift [ppm]", "LS [L<>G]"};
	tabDat = Join[Transpose[Round[{amp, lw, shift, ls}, .001]], phi];

	tab = TableForm[tabDat, TableHeadings -> {Style[#, Bold] & /@ rowName, Style[#, Bold] & /@ colName}, TableSpacing -> {2, 2}, TableAlignments -> Center];
	dat = Transpose[Prepend[Transpose[Prepend[tabDat, colName]], Prepend[rowName, ""]]] /. {
		Style["\!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\) [deg]", Bold] -> "\[Theta]0 [deg]", Style["\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\) [ms]", Bold] -> "\[Theta]1 [ms]"};

	Switch[out,
		"tab", tab,
		"dat", dat,
		"both", {tab, dat}
	]
]


(* ::Subsection::Closed:: *)
(*CompareSpectraFitPlot*)


SyntaxInformation[CompareSpectraFitPlot] = {"ArgumentsPattern" -> {_, _, _, _.}}

CompareSpectraFitPlot[ppmPl_, specPlot_, fitPlot_, ranPpm_:Full] := Block[{ran, sp, error, errorF},
	ran = {-1, 1} Max[Abs[specPlot], Abs[fitPlot]];
	sp = 2;
	error = specPlot - fitPlot;
	errorF = GaussianFilter[error,Round[Length[error]/10]];

	Column[{
		FlipView[{
			Column[{
				PlotSpectra[ppmPl, specPlot - errorF, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "Abs", PlotLabel->"Raw signal"],
				PlotSpectra[ppmPl, specPlot - errorF, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "ReIm"]
			}],
			Column[{
				PlotSpectra[ppmPl, fitPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "Abs", PlotLabel->"Fitted signal"],
				PlotSpectra[ppmPl, fitPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "ReIm"]
			}]
		}],
		FlipView[{
			PlotSpectra[ppmPl, error - errorF, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "ReIm", PlotLabel->"Fit error full"],
			PlotSpectra[ppmPl, error - errorF, GridLineSpacing -> sp, PlotRange -> {ranPpm, Full}, Method -> "ReIm", PlotLabel->"Fit error scaled"]
		}]
	}]
]


(* ::Subsection::Closed:: *)
(*CompareFidFitPlot*)


SyntaxInformation[CompareFidFitPlot] = {"ArgumentsPattern" -> {_,_,_}}

CompareFidFitPlot[time_, fidPlot_, fitPlot_] := Block[{ran, sp, error, errorF},
	sp = 2;
	ran = {-1, 1} Max[Abs[fidPlot], Abs[fitPlot]];

	error = fidPlot - fitPlot;
	errorF= GaussianFilter[error,Round[Length[error]/10]];

	Column[{
		FlipView[{
			Column[{
				PlotFid[time, fidPlot - errorF, GridLineSpacing -> sp, PlotRange -> {Full,ran}, Method -> "Abs", PlotLabel->"Raw signal"],
				PlotFid[time, fidPlot - errorF, GridLineSpacing -> sp, PlotRange -> {Full,ran}, Method -> "ReIm"]
			}],
			Column[{
				PlotFid[time, fitPlot, GridLineSpacing -> sp, PlotRange -> {Full,ran}, Method -> "Abs", PlotLabel->"Fitted signal"],
				PlotFid[time, fitPlot, GridLineSpacing -> sp, PlotRange -> {Full,ran}, Method -> "ReIm"]
			}]
		}]
		,
		FlipView[{
			PlotFid[time, error - errorF, GridLineSpacing -> sp, PlotRange -> {Full,ran}, Method -> "ReIm", PlotLabel->"Fit error full"],
			PlotFid[time, error - errorF, GridLineSpacing -> sp, PlotRange -> {Full,Full}, Method -> "ReIm", PlotLabel->"Fit error scaled"]
		}]
	}]
]


(* ::Subsection::Closed:: *)
(*MakeSpectraResultPlot*)


SyntaxInformation[MakeSpectraResultPlot] = {"ArgumentsPattern" -> {_, _, _, _, _.}}

MakeSpectraResultPlot[ppmF_, specF_, {fit_, basisFit_}, names_, ppmRan_] := Block[{
	sp, specFit, resTotPl, errPl, fitPl, resPl, outPl, pRan, pMax,
	lab1, lab2, resFitRI, resFit, resBasPl, met},

	sp = 2;
	met = "ReIm";
	specFit = fit . basisFit;
	pMax = Max[Abs[specFit], Abs[specF]];
	pRan = {-pMax, pMax};

	resTotPl = Column[{
		FlipView[errPl = {
			PlotSpectra[ppmF, specF - specFit, Method -> met, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp, PlotLabel->"Fit error full"],
			PlotSpectra[ppmF, specF - specFit, Method -> met, PlotRange -> {ppmRan, Full}, GridLineSpacing -> sp, PlotLabel->"Fit error scaled"]
		}],
		FlipView[fitPl = {
			Show[
				PlotSpectra[ppmF, specF, Method -> met /. "ReIm" -> "Re", PlotColor -> Red, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp, PlotLabel->"Real signal"],
				PlotSpectra[ppmF, basisFit[[1]], Method -> met /. "ReIm" -> "Re", PlotColor -> Green, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, specFit, Method -> met /. "ReIm" -> "Re", PlotColor -> Black, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp]
			],
			Show[
				PlotSpectra[ppmF, specF, Method -> met /. "ReIm" -> "Im", PlotColor -> Red, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp, PlotLabel->"Imaginary signal"],
				PlotSpectra[ppmF, basisFit[[1]], Method -> met /. "ReIm" -> "Im", PlotColor -> Green, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, specFit, Method -> met /. "ReIm" -> "Im", PlotColor -> Black, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp]
			]
		}]
		,
		FlipView[resPl = {
			PlotSpectra[ppmF, specF, Method -> met, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp, PlotLabel->"Raw signal"],
			PlotSpectra[ppmF, specFit, Method -> met, PlotRange -> {ppmRan, pRan}, GridLineSpacing -> sp, PlotLabel->"Fitted signal"]
		}]
	}, Alignment -> Center];

	lab1 = Style[#, Bold, Black, Large] & /@ {"Error scaled to signal", "Error scaled to Max"};
	lab2 = Style[#, Bold, Black, Large] & /@ {"Fit and signal Re", "Fit and signal Im"};
	outPl = Column[Flatten[{Thread[{lab1, errPl}], Thread[{lab2, fitPl}]}], Alignment -> Center];

	resBasPl = FlipView[{
		resFitRI = PlotSpectra[ppmF, fit basisFit, Method -> "ReIm", PlotColor -> Red, SpectraSpacing -> 0.2, GridLines -> {}, GridLineSpacing -> sp, 
			PlotLabels -> Prepend[names, "spline"],PlotRange->{ppmRan,Full}, PlotLabel->"Fitted basis spectra real and imaginary"],
		resFit = PlotSpectra[ppmF, fit basisFit, Method -> "Abs", PlotColor -> Red, SpectraSpacing -> 0.2, GridLines -> {}, GridLineSpacing -> sp, 
			PlotLabels -> Prepend[names, "spline"],PlotRange->{ppmRan,Full}, PlotLabel->"Fitted basis spectra absolute"]
	}];

	{resTotPl, resBasPl, {errPl, fitPl, resPl, outPl, resFitRI, resFit}}
]


(* ::Subsection::Closed:: *)
(*SpectraFitResult*)


Options[SpectraFitResult] = {PlotRange -> Full};

SyntaxInformation[SpectraFitResult] = {"ArgumentsPattern" -> {_, _, _, _, _, _, _., OptionsPattern[]}};

SpectraFitResult[specF_, {fit_, basisFit_}, te_, {dw_, gyro_}, {pars_, names_, metRef_, log_}, ops : OptionsPattern[]] := SpectraFitResult[specF, {fit, basisFit}, te, {dw, gyro}, {pars, names, metRef, log}, {}, ops]

SpectraFitResult[specF_, {fit_, basisFit_}, te_, {dw_, gyro_}, {pars_, names_, metRef_, log_}, plots_, OptionsPattern[]] := Block[{resTab, ppm, phase, resTotPl, resBasPl, specFit, fidF, fidFit, time, comFidPlot, specPlot, fitPlot, comSpecPlot, logging, ran},
	NotebookClose[resultWindow];
	resTab = FitSpectraResultTable[fit, pars, names, metRef];
	ran = OptionValue[PlotRange];

	(*make result plot*)
	ppm = GetPpmRange[specF, dw, gyro];
	phase = -pars[[3]];
	{resTotPl, resBasPl} = MakeSpectraResultPlot[ppm, specF, {fit, basisFit}, names, ran][[1 ;; 2]];

	(*make fitted spectra and fids*)
	specFit = fit . basisFit;
	fidF = ShiftedInverseFourier[specF];
	fidFit = ShiftedInverseFourier[specFit];

	(*make fid plot*)
	time = GetTimeRange[fidF, dw];
	comFidPlot = CompareFidFitPlot[time, fidF, fidFit];

	(*correct fitted spectra*)
	specPlot = CorrectTESpec[ApodizePadSpectra[PhaseShiftSpectra[specF, ppm, gyro, phase]], dw, te, gyro, ran];
	fitPlot = CorrectTESpec[ApodizePadSpectra[PhaseShiftSpectra[specFit, ppm, gyro, phase]], dw, te, gyro, ran];

	(*make spectra plot*)
	ppm = GetPpmRange[specPlot, dw, gyro];
	comSpecPlot = CompareSpectraFitPlot[ppm, specPlot, fitPlot, ran];

	(*make logging*)
	logging = MakeLogging[log];
	(*make result window*)
	resultWindow = CreateDialog[{
		CancelButton["Close", DialogReturn[]],
		TabView[{
			"Fit result spectra" -> comSpecPlot, 
			"Fit result fid" -> comFidPlot, 
			"Result table" -> resTab, 
			"Fit log" -> logging,
			If[plots === {}, Nothing, "Fit Calibration" -> Column[plots]], 
			"Raw fit result" -> resTotPl, 
			"Raw fitted basis" -> resBasPl
		}, Alignment -> Center]
	}, Background -> White];
]

MakeLogging[log_] := Block[{tmp},
	tmp = DeleteCases[Flatten[{If[StringQ[#], StringSplit[#, "  "], #]}], ""] & /@ log;
	tmp = If[Length[#] == 1 && Head[#[[1]]] === String, Prepend[#, ""], #] & /@ tmp;
	Grid[{{Grid[tmp[[1 ;; 22]], Alignment -> Left], Grid[tmp[[23 ;;]], Alignment -> Left]}}, Alignment -> Top, Spacings -> 4]
]





(* ::Subsection::Closed:: *)
(*CSIInterface31P*)


Options[CSIInterface] = {SpectraFieldStrength -> 7,SpectraNucleus->"31P"};

SyntaxInformation[CSIInterface] = {"ArgumentsPattern" -> {_., _., OptionsPattern[]}};

CSIInterface[opts : OptionsPattern[]] := CSIInterface["", {0, 0}, opts]

CSIInterface[tei_?NumberQ, bwi_?NumberQ, opts : OptionsPattern[]] := CSIInterface["", {tei, bwi}, opts]

CSIInterface[file_?StringQ, opts : OptionsPattern[]] := CSIInterface[file, {0, 0}, opts]

CSIInterface[file_?StringQ, {tei_?NumberQ, bwi_?NumberQ}, OptionsPattern[]] := Module[{
	f, te, bw, nuc, field, gyro, names, met, metSel, metRef, method, plot, kLoad, rec, spectraC, dec, line, fine, den, z, y, x, 
	sPhase, status, statusP, kspace, noise,	header, type, ham, spectra, spectraR, spec, proc, shift, fids, specs, 
	table, fit, basisFit, errorFit, pars, log, plots, specF, fitted, xm, ym, zm, dn, dc, mr, dw, nSamp, filt, teu,
	fileSave, spectraPlot, lab, fovZ, fovY, fovX, coils, nCoils, teE,
	statPart, loadPart, reconPart, plotPart, fitPart, closePart
	},

	NotebookClose[interfaceWindow];
	NotebookClose[resultWindow];
	NotebookClose[plotWindow];

	(*convert input to editable parameters*)
	f = file; te = tei; bw = bwi;

	(*fixed parameters and default values*)
	nuc = OptionValue[SpectraNucleus]; 
	field = OptionValue[SpectraFieldStrength];
	gyro = GetGyro[nuc, field];

	(*nucleus specific settings*)
	names = met = metSel = {"PE", "PC", "Piex", "Piin", "GPE", "GPC", "PCr", "ATP", "NAD", "UDPG"};
	metRef = "PCr"; 

	(*initial manipulator values*)
	method = "WSVD"; plot = "Cor";
	dec = False; 
	line = fine = den = filt = True;

	(*initial button options*)
	kLoad = rec = sPhase = fitted = proc = dn = dc = False;
	mr = ""; 
	spectraC = shift = coils= 0;

	(*initialize fit coors*)
	z = y = x = fovZ = fovY = fovX = zm = ym = xm = 0;

	(*monitoring*)
	status = ""; statusP = False;

	(*------------ STATUS -------------*)
	statPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell["Status "], Dynamic[Row[{TextCell[status, Bold, Black], If[statusP, ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate], TextCell[""]]}, " "]]}
	};

	(*------------ LOADING -------------*)
	loadPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Import Data", Bold, 14]},

		(*select file button*)
		{TextCell["Select CSI list data ", FontFamily -> "Helvetica"], FileNameSetter[Dynamic[f], "Open", {"List Data" -> {"*.list", "*.data"}}, ImageSize -> 175]},

		(*load button*)
		{Button["Load data",
			NotebookClose[resultWindow];
			NotebookClose[plotWindow];
			status = "Loading data"; statusP = True;
			rec = proc = dn = dc = sPhase = False;
			{{kspace, noise}, {header, type}} = ReadListData[f, False];
			If[filt,
				ham = Total@Unitize[Abs[OrderKspace[kspace, type, {"N_aver", "N_kz", "N_ky", "N_kx"}][[1]]]];
				kspace = OrderKspace[TotalType[kspace, type, {"N_aver"}], {"N_chan", "N_samp", "N_kz", "N_ky", "N_kx"}][[1]];
				,
				kspace = OrderKspace[MeanType[kspace, type, {"N_aver"}], {"N_chan", "N_samp", "N_kz", "N_ky", "N_kx"}][[1]];
				ham = MakeHammingFilter[Dimensions[kspace][[-3;;]]];
				kspace = Map[ham #&,kspace,{2}];
			];
			coils = nCoils = Range[Length[kspace]];
			kLoad = True;
			status = "Done loading data!"; statusP = False;
			,
			Background -> Dynamic[If[kLoad, RGBColor[0.86, 0.97, 0.77], Automatic]], Enabled -> Dynamic[FileExistsQ[f]], Method -> "Queued", ImageSize -> 175
		],Dynamic[f]},
		{TextCell[" Filter data  "], SetterBar[Dynamic[filt], {True -> " k-space weighting ", False -> " Hamming "}]}
	};

	(*------------ ACQUISITION / RECONSTRUCTION -------------*)
	reconPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["acquisition parameters", Bold, 14]},

		(*acquisition parameters*)
		{TextCell["Echo time [ms] "], InputField[Dynamic[te], Number]},
		{TextCell["Bandwidth [Hz] "], InputField[Dynamic[bw], Number]},

		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Data reconstruction", Bold, 14]},

		(*reconstruct button*)
		{Button["Reconstruct data",
			NotebookClose[resultWindow];
			NotebookClose[plotWindow];
			status = "Reconstructing using " <> method <> ""; statusP = True;
			rec = proc = dn = dc = sPhase = False; mr = "";
			mr = method;
			dw = 1./bw;
			teu = te/1000;
			coils = If[coils==={},All,coils];
			spectra = spectraR = CoilWeightedReconCSI[kspace[[coils]], noise[[coils]], header, Method -> method];
			rec = True; 
			If[den, spectra = DenoiseCSIdata[spectra]; proc = dn = True; ];
			If[dec, spectra = DeconvolveCSIdata[spectra, ham]; proc = dc = True;];
			{zm, ym, xm, nSamp} = Dimensions[spectra];

			(*find shift*)
			spec = Mean@Flatten[spectra, 2];
			shift = 0(*new find shift function to be added*); 
			status = "Done reconstructing!"; statusP = False;
			,
			Background -> Dynamic[If[rec, RGBColor[0.86, 0.97, 0.77], Automatic]], Enabled -> Dynamic[kLoad && bw =!= 0], Method -> "Queued", ImageSize -> 175
		],
		Row[{TextCell[" Denoise data  "], Checkbox[Dynamic[den]], TextCell["   Deconvolve data  "], Checkbox[Dynamic[dec]]}, ""]},
		{TextCell[" Select coils  "], Dynamic[If[coils===0,TextCell[""], TogglerBar[Dynamic[coils], nCoils]]]},

		(*select method*)
		{TextCell[" Reconstruction method  "], SetterBar[Dynamic[method], {"Roemer" -> " Roemer ", "WSVD" -> " WSVD "}]},

		(*report*)
		{Dynamic[TextCell[If[rec, "PCr shift is " <> ToString[Round[shift, .01]] <> " ppm ", ""]]], Dynamic[TextCell[If[rec, "Reconstructed using: " <> mr <> If[dn, ", de-noised", ""] <> If[dc, ", de-convolved", ""], ""]]]}
		};

	(*------------ PLOTTING -------------*)
	plotPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Data Plotting", Bold, 14]},

		(*plotting data button*)
		{Button["Plot CSI data", 
			NotebookClose[resultWindow];
			NotebookClose[plotWindow];
			spectraPlot = Switch[plot,
				"Raw", spectraR,
				"Proc", spectra,
				"Cor",
				If[!sPhase,
					status = "Correcting phase of spectra"; statusP = True;
					spectraC = Map[PhaseCorrectSpectra[ApodizePadSpectra[ShiftSpectra[#, {dw, gyro}, -shift]], dw, teu, gyro, {10, -20}]&, spectra, {-2}];
					status = "Done phase correcting spectra!"; statusP = False; sPhase = True
				];
				spectraC
			];
			status = "Plotting the CSI data"; 
			PlotCSIData[spectraPlot, {dw, gyro}, PlotRange -> {10, -20}]
			,
			Enabled -> Dynamic[rec], Method -> "Queued", ImageSize -> 175
		],

		(*select plot method*)
		Dynamic[SetterBar[Dynamic[plot], {"Raw" -> " Raw ", If[proc, "Proc" -> " Processed ", Nothing], "Cor" -> " Phase corrected "}]]},
		{Button["Save spar/sdat",
			(*saving make better!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
			spectraPlot = Switch[plot,
				"Raw", teE=te; spectraR,				
				"Proc",teE=te; spectra,
				"Cor",
				If[!sPhase,
					status = "Correcting phase of spectra"; statusP = True;
					spectraC = Map[PhaseCorrectSpectra[ApodizePadSpectra[ShiftSpectra[#, {dw, gyro}, -shift]], dw, teu, gyro, {10, -20}, True]&, spectra, {-2}];
					status = "Done phase correcting spectra!"; statusP = False;	sPhase = True
				];
				teE=0.;
				spectraC
			];
			status = "Saving the CSI data";
			fileSave = SystemDialogInput["FileSave"]; 
			lab = mr <> If[plot==="Raw","",If[dn, ", de-noised", ""] <> If[dc, ", deconvolve", ""]]<>If[plot ==="Cor"," ,phase corrected",""];
			If[fileSave =!= $Canceled, ExportSparSdat[fileSave, spectraPlot, {bw, teE}, {gyro, nuc}, {fovZ, fovY, fovX}]];
			(*saving make better!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
			,
			Enabled -> Dynamic[rec&&fovY>0&&fovZ>0&&fovX>0], Method -> "Queued", ImageSize -> 175],
		Row[{
			{fovZ, fovY, fovX} = Round[{fovZ, fovY, fovX}];
			TextCell["  vox z [mm] "], InputField[Dynamic[fovZ, (fovZ = Round[#]) &], Number, FieldSize -> 3, ContinuousAction -> True, Background -> Dynamic[If[fovZ>0, White, RGBColor[1, 0.9, 0.9]]]],
			TextCell["   vox y [mm] "], InputField[Dynamic[fovY, (fovY = Round[#]) &], Number, FieldSize -> 3, ContinuousAction -> True, Background -> Dynamic[If[fovY>0, White, RGBColor[1, 0.9, 0.9]]]],
			TextCell["   vox x [mm]  "], InputField[Dynamic[fovX, (fovX = Round[#]) &], Number, FieldSize -> 3, ContinuousAction -> True, Background -> Dynamic[If[fovX>0, White, RGBColor[1, 0.9, 0.9]]]]
		}, " "]}
	};

	(*------------ FITTING -------------*)
	fitPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Spectra Fitting", Bold, 14]},
		(*fitting button*)
		{Button["Fit selected spectra",
			NotebookClose[resultWindow];
			NotebookClose[plotWindow];
			(*basis spectra*)
			status = "Generating basis spectra"; statusP = True;
			{names, fids, specs, table} = GetSpectraBasisFunctions[metSel, {"ATP"}, BasisSequence -> {"PulseAcquire", teu}, 
				SpectraSamples -> nSamp, SpectraBandwidth -> bw, SpectraPpmShift -> 0, SpectraFieldStrength -> field, SpectraNucleus -> nuc];
			status = "Done generating basis spectra!"; statusP = False;

			(*fitting*)
			status = "Start fitting spectra"; statusP = True;
			specF = ShiftSpectra[spectra[[z, y, x]], {dw, gyro}, -shift];
			{fit, basisFit, errorFit, pars, log, plots} = FitSpectra[specs, specF, {10, -20}, dw, { {0, -2.52, -7.56, -16.15}, {2, 1, 1, 1}},
				PaddingFactor -> 2, SpectraPpmShift -> 0, SplineSpacingFactor -> 2.5, SpectraOutputPlots -> True, 
				SpectraNucleus -> nuc, SpectraFieldStrength -> field, FitLineShape -> line, FineTuneFit -> fine];
			status = "Done fitting spectra!"; statusP = False;
			fitted = True;
			,
			Background -> Dynamic[If[fitted, RGBColor[0.86, 0.97, 0.77], Automatic]], Enabled -> Dynamic[rec && 1 <= x <= xm && 1 <= y <= ym && 1 <= z <= zm], Method -> "Queued", ImageSize -> 175],

			(*size reporting*)
			Dynamic[TextCell[If[rec, "The CSI grid size is " <> ToString[Dimensions[spectra][[;; -2]]], ""]]]
		},

		(*input fit parameters*)
		{TextCell["Input coordinate to fit "], Row[{
			{z, y, x} = Round[{z, y, x}];
			TextCell["  z  "], InputField[Dynamic[z, (z = Round[#]) &], Number, FieldSize -> 2, ContinuousAction -> True, Background -> Dynamic[If[1 <= z <= zm, White, RGBColor[1, 0.9, 0.9]]]],
			TextCell["   y  "], InputField[Dynamic[y, (y = Round[#]) &], Number, FieldSize -> 2, ContinuousAction -> True, Background -> Dynamic[If[1 <= y <= ym, White, RGBColor[1, 0.9, 0.9]]]],
			TextCell["   x  "], InputField[Dynamic[x, (x = Round[#]) &], Number, FieldSize -> 2, ContinuousAction -> True, Background -> Dynamic[If[1 <= x <= xm, White, RGBColor[1, 0.9, 0.9]]]]
		}, " "]},

		(*fitting options*)
		{TextCell["Fitting options "], Row[{TextCell[" Fine tune fit  "], Checkbox[Dynamic[fine]], TextCell["   Fit line shape  "], Checkbox[Dynamic[line]]}, ""]},
		(*select metabolites*)
		{TextCell["Select metabolites to fit "], TogglerBar[Dynamic[metSel], met]},

		(*------------ RESULTS -------------*)
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Fitting Results", Bold, 14]},

		(*show fit results button*)
		{TextCell["Reference metabolite "], Dynamic[SetterBar[Dynamic[metRef], names]]},
		{Button["Show results",
			NotebookClose[resultWindow];
			NotebookClose[plotWindow];
			status = "Generating fit results"; statusP = True;
			SpectraFitResult[specF, {fit, basisFit}, teu, {dw, gyro}, {pars, names, metRef, log}, plots, PlotRange -> {10, -20}];
			status = "Done with generating fitting results!"; statusP = False;
			,
			Enabled -> Dynamic[fitted], Method -> "Queued", ImageSize -> 175]
		}
	};

	(*------------ CLOSING -------------*)
	closePart ={
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], CancelButton["Close GUI", DialogReturn[
			NotebookClose[resultWindow];
			NotebookClose[plotWindow];
			Clear[kspace, spectra]
			], ImageSize -> 175]},{TextCell[""]
		},
		{TextCell[""]}
	};

	grid = Join[statPart, loadPart, reconPart, plotPart, If[nuc==="31P",fitPart,{}], closePart];

	dialogGrid = Manipulate[(*Column[{{f,te,bw},{den,dec,method,plot},{fine,line,{x,y,z}},{metSel},{metRef},{Dimensions[kspace]};}],*)"",
		Grid[grid, Alignment -> {{Right, Left}}, ItemSize -> {{18, 37}, Automatic}],
		AppearanceElements -> None,
		SynchronousUpdating -> False,
		Paneled -> False,
		ControlPlacement -> Right
	];

	(*------------ DIALOG -------------*)
	interfaceWindow = CreateDialog[dialogGrid, Background -> White, WindowSize -> All, WindowTitle -> "CSI Processing GUI"];
];


(* ::Subsubsection::Closed:: *)
(*DKI*)


(*TensMinDKI[s_,ls_,bmat_,bmatI_] := bmatI.ls*)
TensMinDKI = Compile[{{s, _Real, 1}, {bmatI, _Real, 2}},
	If[Total[s]==0.,
		{0.,0.,0.,0.,0.,0.,0.},
		bmatI . s
	]
,RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


(* ::Subsubsection::Closed:: *)
(*NLS*)


TensMinNLS[s_,ls_,bmat_,bmatI_] := 
Module[{v,xx,yy,zz,xy,xz,yz,init,tens,sol},
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},
		tens,
		v={xx,yy,zz,xy,xz,yz,tens[[7]]};
		init=Thread[{v[[1;;6]],tens[[1;;6]]}];
		sol=FindMinimum[.5 Total[(s-Exp[bmat . v])^2],init][[2]];
		v/.sol
	]
]


(* ::Subsubsection::Closed:: *)
(*NLS*)


TensMinGMM[s_,ls_,bmat_,bmatI_] := 
Module[{v,xx,yy,zz,xy,xz,yz,init,tens,res,w},
	s;
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={xx,yy,zz,xy,xz,yz,tens[[7]]};
		init=Thread[{v[[1;;6]],tens[[1;;6]]}];
		v/.FindMinimum[(
			res=ls-bmat . v;
			w=1/(res^2+Mean[res]^2);
			.5 Total[(w/Mean[w])*(res)^2]
		),init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CLLS*)


TensMinCLLS[s_,ls_,bmat_,bmatI_] := 
Module[{v,r0,r1,r2,r3,r4,r5,init,tens},
	s;
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={r0^2,r1^2+r3^2,r2^2+r4^2+r5^2,r0 r3,r0 r4,r3 r4+r1 r5,tens[[7]]};
		init=Thread[{{r0,r1,r2,r3,r4,r5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(ls-bmat . v)^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CWLLS*)


TensMinCWLLS[s_,ls_,bmat_,bmatI_] := 
Module[{v,r0,r1,r2,r3,r4,r5,init,tens,std=1,wMat},
	bmatI;
	wMat=Transpose[bmat] . DiagonalMatrix[s^2/std^2];
	tens = LinearSolve[wMat . bmat, wMat . ls];
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={r0^2,r1^2+r3^2,r2^2+r4^2+r5^2,r0 r3,r0 r4,r3 r4+r1 r5,tens[[7]]};
		init=Thread[{{r0,r1,r2,r3,r4,r5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(s^2/std^2)*(ls-bmat . v)^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CNLS*)


TensMinCNLS[s_,ls_,bmat_,bmatI_] := 
Module[{v,r0,r1,r2,r3,r4,r5,init,tens},
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={r0^2,r1^2+r3^2,r2^2+r4^2+r5^2,r0 r3,r0 r4,r3 r4+r1 r5,tens[[7]]};
		init=Thread[{{r0,r1,r2,r3,r4,r5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(s-Exp[bmat . v])^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*ExtendedCholeskyDecomposition*)


ExtendedCholeskyDecomposition[tm_] := Block[{n,beta,theta,cm,lm,dm,em,j},
	n=Length[tm];
	beta=Max[{Max[Diagonal[tm]],Max[UpperTriangularize[tm,1]]/Sqrt[n^2-1],10^-15}];
	cm=DiagonalMatrix[Diagonal[tm]];
	lm=dm=em=ConstantArray[0,{n,n}];
	Table[
		If[j==1,
			(*j=1 make first column cm equal to tm*)
			cm[[j+1;;,j]]=tm[[j+1;;,j]];
			,
			(*j>1 fill lm matrix*)
			lm[[j,;;j-1]]=cm[[j,;;j-1]]/(Diagonal[dm][[;;j-1]]/.(0.->Infinity));
			If[j<n,
				cm[[j+1;;,j]]=tm[[j+1;;,j]]-lm[[j,j-1;;]] . Transpose[cm[[j+1;;,j-1;;]]]
				];
			];
		theta=If[j==n,0,Max[Abs[cm[[j+1;;,j]]]]];
		dm[[j,j]]=Max[{Abs[cm[[j,j]]],theta^2/beta}];
		em[[j,j]]=dm[[j,j]]-cm[[j,j]];
		cm=cm-DiagonalMatrix[PadLeft[(1/(dm[[j,j]]/.(0.->Infinity)))*cm[[j+1;;,j]]^2,n]];
	,{j,1,3}];
	lm=lm+IdentityMatrix[n];
	Transpose[lm . MatrixPower[dm,.5]]
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
