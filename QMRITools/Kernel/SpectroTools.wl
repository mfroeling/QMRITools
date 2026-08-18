(* ::Package:: *)

(* ::Title:: *)
(*QMRITools SpectroTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`SpectroTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`SpectroTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


ReadjMRUI::usage = 
"ReadjMRUI[file] read a jMRUI spectrum file. 
Output is the {time, spec, {beginTime, samplingInterval}}."

PadFid::usage = 
"PadFid[fid] pads the fid with zeros to increase its length."

PadEcho::usage =
"PadEcho[echo] pads the echo with zeros to increase its length."

ApodizeFid::usage = 
"ApodizeFid[fid] performs apodization on the fid. The apodization function is set with the option ApodizationFunction."

ApodizeEcho::usage =
"ApodizeEcho[echo] performs apodization on the echo. The apodization function is set with the option ApodizationFunction."

ApodizePadFid::usage = 
"ApodizePadFid[fid] performs apodization on the fid and pads the fid with zeros to increase its length."

ApodizePadEcho::usage =
"ApodizePadEcho[echo] performs apodization on the echo and pads the echo with zeros to increase its length."

PadSpectra::usage = 
"PadSpectra[spec] doubles the number of spectral points while maintaining the dwell time."

ApodizeSpectra::usage = 
"ApodizeSpectra[spec] performs apodization of the spectra. The apodization function is set with the option ApodizationFunction."

ApodizePadSpectra::usage = 
"ApodizePadSpectra[spec] and doubles the number of spectral points while maintaining the dwell time."


GetTimePpmRange::usage = 
"GetTimePpmRange[spec, {dt, field, nuc}] get the timing of the fid and the ppm values of the spec where dt is the well time in ms, field the field strength in Tesla and nuc the nucleus available in GyromagneticRatio. 
GetTimePpmRange[spec, dt, field, nuc] get the timing of the fid and the ppm values of the spec. 
GetTimePpmRange[spec, dt, gyro] get the timing of the fid and the ppm values of the spec."

GetTimeRange::usage = 
"GetTimeRange[fid, dt] get the timing of the fid where dt is the well time in ms."

GetPpmRange::usage = 
"GetPpmRange[spec, {dt, field, nuc}] get the ppm values of the spec where dt is the well time in ms, field the field strength in Tesla and nuc the nucleus available in GyromagneticRatio. 
GetPpmRange[spec, dt, field, nuc] get the ppm values of the spec. 
GetPpmRange[spec, dt, gyro] get the ppm values of the spec."

GetGyro::usage = 
"GetGyro[nuc, field] geth the gyromagnetic ratio with field the field strength in Tesla and nuc the nucleus available in GyromagneticRatio."


PhaseShiftSpectra::usage = 
"PhaseShiftSpectra[spectra, phi0] applies the 0th order phase phi0 to the spectra. 
PhaseShiftSpectra[spectra, ppm, gyro, phi1] applies the 1st order phase phi1 to the spectra. The ppm can be obtained using GetPpmRange and gyro with GetGyro. 
PhaseShiftSpectra[spec, ppm, gyro, {phi0, phi1}] applies the 0th and 1st order phases {phi0, phi1} to the spectra. The ppm can be obtained using GetPpmRange and gyro with GetGyro.

The 0th order phase phi0 is in radians and the 1st order phase phi1 is in ms."

PhaseShiftFid::usage = 
"PhaseShiftFid[spectra, phi0] applies the 0th order phase phi0 to the FID."

TimeShiftFid::usage = 
"TimeShiftFid[fid, time, gam] applies a line broadening with linewidth gam and a Voigt line shape to the fid. The time can be obtained using GetTimeRange.
TimeShiftFid[fid, time, {gamL, gamG}] applies a line broadening with linewidth gamG \"Gaussian\" and gamL \"Lorentzian\".
TimeShiftFid[fid, time, gyro, {gam, eps}] applies a line broadening with linewidth gam to the fid and a phase eps that results in eps ppm shift of the spectra. The gyro can be obtained with GetGyro.
TimeShiftFid[fid, time, gyro, {{gamL, gamG}, eps}] applies a line broadening with linewidth linewidth gamG \"Gaussian\" and gamL \"Lorentzian\" to the fid and a phase eps that results in eps ppm shift of the spectra.

The linewidth gam is given in ms and the spectra shift eps is given in ppm."

TimeShiftEcho::usage =
"TimeShiftEcho[fid, time, gam] applies a line broadening with linewidth gam and a Voigt line shape to the fid. The time can be obtained using GetTimeRange.
TTimeShiftEcho[fid, time, {gam, f}] applies a line broadening with linewidth gam and a custom line shape f to the fid (f=0, \"Gaussian\", f=1 \"Lorentzian\").
TTimeShiftEcho[fid, time, gyro, {gam, eps}] applies a line broadening with linewidth gam to the fid and a phase eps that results in eps ppm shift of the spectra. The gyro can be obtained with GetGyro.
TTimeShiftEcho[fid, time, gyro, {gam, eps, f}] applies a line broadening with linewidth gam using a custom line shape f to the fid and a phase eps that results in eps ppm shift of the spectra.

The linewidth gam is given in ms and the spectra shift eps is given in ppm."

ShiftSpectra::usage = 
"ShiftSpectra[spectra, {dw, gyro}, shift] shifts the spectra by shift. The shift is in ppm." 


ChangeDwellTimeFid::usage = 
"ChangeDwellTimeFid[fid, dt, dtNew] changes the sampling time of an fid from dwell time dt to dwell time dtNew."

PhaseCorrectSpectra::usage =
"PhaseCorrectSpectra[spec] performs 0th order phase correction of the spectra by minimizing the difference between the real and absolute spectra value.
PhaseCorrectSpectra[spec, dw] performs 0th order phase correction of the spectra using Hankel matrix SVD fitting.
PhaseCorrectSpectra[spec, dw, te] := performs 0th and 1st order phase correction of the spectra using Hankel matrix SVD fitting. The first order phase is corrected by padding the fid with the missing values in the time before the TE.
PhaseCorrectSpectra[spec, dw, te, gyro, ppmRan] performs 0th and 1st order phase correction of the spectra using Hankel matrix SVD fitting. Only the part of the spectra in the ppmRan is used for optimization."

CorrectTESpec::usage = 
"CorrectTESpec[spectra, dw, te] corrects the spectra for 1st order phase by extrapolating the missing FID samples in the TE using Hankel matrix SVD analysis.
CorrectTESpec[spectra, dw, te, gyro, ppmRan] corrects the spectra for 1st order phase by extrapolating the missing FID samples in the TE using Hankel matrix SVD analysis. Only the part of the spectra in the ppmRan is used for optimization."

CorrectTEFid::usage = 
"CorrectTEFid[fid, dw, te] corrects the fid for 1st order phase by extrapolating the missing FID samples in the TE using Hankel matrix SVD analysis.
CorrectTEFid[fid, dw, te, gyro, ppmRan] corrects the fid for 1st order phase by extrapolating the missing FID samples in the TE using Hankel matrix SVD analysis. Only the part of the spectra in the ppmRan is used for optimization."


FitSpectra::usage = 
"FitSpectra[specBasis, spec, {st,end}, dt, {lwVals,lwAmps}] Fits the basis spectra from GetSpectraBasisFunctions to the spec overt the ppm range {st, end} and dt the dwell time."

GetSpectraBasisFunctions::usage = 
"GetSpectraBasisFunctions[{met1, ..., metN}] generates a list of spectra basis functions with names met1 to metN. The names are strings and are the metabolites available in GetSpinSystem.
GetSpectraBasisFunctions[{{prop1}, ..., {propN}}] generates a list of spectra basis functions with properties prop1 to propN. The properties are those specified in MakeSpinSystem.
GetSpectraBasisFunctions[inp, split] generates a list of spectra basis functions. Each metabolite name present in the list split wil be split in individual spectra per peak."


PlotCSIData::usage =
"PlotCSIData[spectra, {dwell, gyro}] plots the CSI spectra which has dimensions {z,y,x,nSamp}. The ppm axes is determined by dwell and gyro. Gyro can be obtained with GetGyro.
PlotCSIData[spectra, {dwell, field, nuc}] plots the CSI spectra which has dimensions {z,y,x,nSamp}. The ppm axes is determined by dwell and field and nuc."

PlotFid::usage = 
"PlotFid[fid, dwell] plots the fid assuming dwell as the sampling time.
PlotFid[time, fid] plot the fid where time is the timing of the fid which can be obtained with GetTimeRange."

PlotSpectra::usage = 
"PlotSpectra[spectra, {dwell, gyro}] plots the spectra, the ppm axes is determined by dwell and gyro. Gyro can be obtained with GetGyro.
PlotSpectra[spectra {dwell, field, nuc}] plots the spectra, the ppm axes is determined by dwell field and nuc.
PlotSpectra[ppm, spectra] plots the spectra where ppm is the pmm range of the spectra which can be obtained with GetPpmRange." 


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


SpectraSimulator::usage = 
"SpectraSimulator[] opens a interactive spectra simulator."


CSIInterface::usage = 
"CSIInterface[] opens the CSI interface. Function not done.
CSIInterface[te, bw] opens the CSI interface with known te and bw.
CSIInterface[file] opens the CSI interface with the data from file loaded.
CSIInterface[file, {tei, bwi}] opens the CSI interface with the data from file loaded with known te and bw."

SpectraFitResult::usage = 
"SpectraFitResult[spec, {fit, basisFit}, te, {dw, gyro}, {pars, names, metRef, log}, plots, OptionsPattern[]] function not done."


ImportSparSdat::usage = 
"ImportSparSdat[fileSpar, fileSdat] imports fileSpar and fileSdat files. Function not done."

ExportSparSdat::usage =
"ExportSparSdat[file, specs, {bw ,te}, {gyro ,nuc}] exports specs to file. Function not done."


(* ::Subsection::Closed:: *)
(*Options*)


ApodizationFunction::usage = 
"ApodizationFunction is an options for ApodizeFid, ApodizeSpectra, ApodizePadFid, and ApodizePadSpectra. Values can be \"Hanning\", \"Hamming\", \"Gaussian\", \"Lorentzian\", and \"Voigt\"."

PaddingFactor::usage = 
"PaddingFactor is an option for PadFid, PadSpectra, ApodizePadFid, ApodizePadSpectra and FitSpectra. It Specifies with which factor to lengthen the fid."

BasisSequence::usage = 
"BasisSequence is an option for GetSpectraBasisFunctions and specifies which sequence to use."

SpectraSamples::usage =
"SpectraSamples is an option for GetSpectraBasisFunctions and sets the number of samples in the spectra."

SpectraBandwidth::usage =
"SpectraBandwidth is an option for GetSpectraBasisFunctions and sets the bandwidth of the spectra." 

SpectraNucleus::usage =
"SpectraNucleus is an option for GetSpectraBasisFunctions and FitSpectra and specifies which nucleus to Simulate or fit, see GyromagneticRatio."

SpectraPpmShift::usage = 
"SpectraPpmShift is an option for GetSpectraBasisFunctions and FitSpectra and defines how much the center frequency is shifted, default is water at 4.65 ppm."

SpectraFieldStrength::usage = 
"SpectraFieldStrength is an option for GetSpectraBasisFunctions and FitSpectra and sets the field strength at which the simulations and fitting is performed."

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

SpectraSpacing::usage = 
"SpectraSpacing is an option for PlotSpectra and defines the amount of spacing between spectra when multiple spectra are plotted."

SparName::usage = 
"SparName is an option for ExportSparSdat."

SparOrientation::usage = 
"SparOrientation is an option for ExportSparSdat."

SparID::usage = 
"SparID is an option for ExportSparSdat."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*Import data*)


(* ::Subsubsection::Closed:: *)
(*ReadJMRUI file*)


SyntaxInformation[ReadjMRUI]={"ArgumentsPattern"->{_}}

ReadjMRUI[file_]:=Block[{imp,data,head,series,pts,spec,time},
	imp=Import[file,"Data"];
	data=Select[imp,AllTrue[#,NumericQ]&&#=!={}&];
	head=Flatten[Select[imp,!AllTrue[#,NumericQ]&&#=!={}&]];
	head=(StringTrim/@StringSplit[#<>" ",":"])&/@Select[head,StringContainsQ[#,":"]&];
	head[[2;;,2]]=(ToExpression[StringReplace[#,"E"->" 10^" ]]&/@head[[2;;,2]])/.Null->0;
	head=Thread[head[[All,1]]->head[[All,2]]];
	series="DatasetsInFile"/.head;
	pts="PointsInDataset"/.head;
	data=Partition[data,pts];
	spec=Reverse/@(data[[All,All,3]]+data[[All,All,4]]I);
	time=(data[[All,All,1]]+data[[All,All,2]]I);
	{time,spec,N@{"BeginTime","SamplingInterval"}/1000/.head}
]


(* ::Subsection::Closed:: *)
(*PhaseCorrection*)


(* ::Subsubsection::Closed:: *)
(*PhaseCorrectSpectra*)


SyntaxInformation[PhaseCorrectSpectra]={"ArgumentsPattern"->{_,_.,_.,_.,_.,_.}}

PhaseCorrectSpectra[spec_] := Exp[-I Quiet[Last[FindMinimum[PhaseCorrectError[spec, phi0], {phi0}]]][[1,2]]] spec


PhaseCorrectError[specI_, phi0_?NumericQ] := PhaseCorrectErrorC[specI, phi0]

PhaseCorrectErrorC = Compile[{{specI, _Complex, 1}, {phi0, _Real, 0}},Block[{specR},
	specR = Re[Exp[-I phi0] specI];
	Total[(Abs[specI] - specR)^2] - Total[specR]
],RuntimeOptions -> "Speed", Parallelization -> True];


PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ] := PhaseCorrectSpectra[spec, dw, 0., 0., Full]

PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, te_?NumberQ] := PhaseCorrectSpectra[spec, dw, te, 0., Full]

PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, gyro_?NumberQ, ppmRan_?ListQ] := PhaseCorrectSpectra[spec, dw, 0., gyro, ppmRan]

PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, te_?NumberQ, gyro_?NumberQ, ppmRan_] := Block[{
		fid, specOut, missing, full, hankelSpec, phi, phi0
	},
	(*create the fid*)
	fid = ShiftedInverseFourier[spec];

	(*perform the HankelFit*)
	{missing, full} = HankelFit[fid ,dw, te, gyro, ppmRan];

	(*create the full fid*)
	specOut = ShiftedFourier[Join[missing, fid][[;;Length[fid]]]];

	(*create the Hankel spectra*)
	hankelSpec = ShiftedFourier[full];
	If[ppmRan =!= Full, 
		hankelSpec = Pick[hankelSpec, Unitize[Clip[GetPpmRange[hankelSpec, dw, gyro], Sort[ppmRan], {0, 0}]], 1]
	];

	(*find the first order phase*)
	Clear[phi0];
	phi = Quiet[Last[FindMaximum[PhaseErrorH[hankelSpec, phi0], phi0]][[1,2]]];
	phi = If[NumberQ[phi], phi, phi[[2]]];
	specOut Exp[-I phi]
]


PhaseErrorH[specI_, phi0_?NumericQ] := PhaseErrorHC[specI, phi0]

PhaseErrorHC = Compile[{{specI, _Complex, 1}, {phi0, _Real, 0}}, 
	Total[Re[Exp[-I phi0] specI]], 
	RuntimeOptions -> "Speed", Parallelization -> True
];


(* ::Subsubsection::Closed:: *)
(*CorrectTESpec*)


SyntaxInformation[CorrectTESpec]={"ArgumentsPattern"->{_,_,_,_.,_.}}

CorrectTESpec[spec_, dw_, te_] := CorrectTESpec[spec, dw, te, 0, Full]

CorrectTESpec[spec_, dw_, te_, gyro_, ppmRan_] := ShiftedFourier[CorrectTEFid[ShiftedInverseFourier[spec], dw, te, gyro, ppmRan]]


(* ::Subsubsection::Closed:: *)
(*CorrectTEFid*)


SyntaxInformation[CorrectTEFid]={"ArgumentsPattern"->{_,_,_,_.,_.}}

CorrectTEFid[fid_, dw_, te_] := CorrectTEFid[fid, dw, te, 0, Full]

CorrectTEFid[fid_, dw_, te_, gyro_, ppmRan_] := Block[{missing},
	(*get the missing time points for the fid*)
	missing = HankelFit[fid ,dw, te, gyro, ppmRan][[1]];
	(*output full fid*)
	Join[missing, fid][[;;Length[fid]]]
]


(* ::Subsubsection::Closed:: *)
(*HankelFit*)


HankelFit[fid_ ,dw_, te_, gyro_, ppmRan_]:=Block[{timeOr, timeMis, henkel, fit},
	(*
	10.1016/0022-2364(87)90023-0
	10.1006/jmra.1995.1116
	10.1016/0022-2364(92)90188-D
	10.1006/jmra.1993.1249
	*)
	(*get the correct timing of the fid and missing values*)
	timeOr = dw Range[1, Length[fid]];
	timeMis = Reverse@Range[te - dw, If[Mod[te, dw] > 0.5 dw, -0.5 dw, 0], -dw];

	(*get the hankel values*)
	henkel = HankelSVDFid[fid, dw, gyro, ppmRan];
	fit = (PseudoInverse[HankelSVDBasisC[timeOr, henkel]] . fid);

	(*missing and full hankel fid*)
	{If[timeMis =!= {}, HankelSVDBasisC[timeMis, henkel].fit, {}], HankelSVDBasisC[timeOr, henkel] . fit}
]


(* ::Subsubsection::Closed:: *)
(*HankelSVDFid*)


HankelSVDFid[fid_, dw_] := HankelSVDFid[fid, dw, 0, Full]

HankelSVDFid[fid_, dw_, gyro_] := HankelSVDFid[fid, dw, gyro, Full]

HankelSVDFid[fid_, dw_, gyro_, ppmRan_] := Block[{
		lMax, mMax, H, U, Utr, Uup, q, freq, decay, ppmVal, sel, ppmMax, ppmMin
	},

	(*get the matrix size*)
	lMax = Round[0.4 Length[fid]];
	mMax = Length[fid] + 1 - lMax;

	(*create the hankel matrix and the singular values*)
	H = fid[[Range[mMax] + #]] & /@ Range[0, lMax - 1];
	U = First@SingularValueDecomposition[H, 32];
	q = Log[Eigenvalues[PseudoInverse[U[[;; -2]]] . U[[2 ;;]]]];

	(*get the frequencies and delay times*)
	decay = Re[q]/dw;
	freq = Im[q]/(2 Pi dw);

	(*select the the correct frequency range and only long decay times*)
	sel = If[ppmRan === Full, 
		UnitStep[Abs[1000/decay] - Pi], 
		ppmVal = -(freq/gyro);
		{ppmMin, ppmMax} = Sort[ppmRan];
		(1 - UnitStep[ppmVal + ppmMin]) UnitStep[ppmVal + ppmMax] UnitStep[Abs[1000/decay] - Pi]
	];
	sel = If[Total[sel] === 0, sel + 1, sel];

	(*give the hankel output*)
	{Pick[decay, sel, 1], Pick[freq, sel, 1]}
]


HankelSVDBasisC = Compile[{{time, _Real, 1}, {hankel, _Real, 2}}, 
	Transpose[Exp[# time] & /@ (hankel[[1]] + 2 Pi I hankel[[2]])], 
	RuntimeOptions -> "Speed"
];


(* ::Subsection:: *)
(*Padding and apodization*)


(* ::Subsubsection::Closed:: *)
(*PadFid*)


Options[PadFid] = {PaddingFactor -> 2}

SyntaxInformation[PadFid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

PadFid[fid_, OptionsPattern[]] := PadRight[fid, Round[OptionValue[PaddingFactor] Length[fid]]]


(* ::Subsubsection::Closed:: *)
(*PadEcho*)


Options[PadEcho] = {PaddingFactor -> 2}

SyntaxInformation[PadEcho] = {"ArgumentsPattern" -> {_,OptionsPattern[]}}

PadEcho[echo_, OptionsPattern[]] := ArrayPad[echo, {Floor[(OptionValue[PaddingFactor]-1) (Length[echo]/2)], Ceiling[(OptionValue[PaddingFactor] -1)(Length[echo]/2)]}]


(* ::Subsubsection::Closed:: *)
(*ApodizeFid*)


Options[ApodizeFid] = {ApodizationFunction -> "Hanning"}

SyntaxInformation[ApodizeFid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeFid[fid_, OptionsPattern[]] := ApodizeFun[Length[fid], OptionValue[ApodizationFunction]] fid


(* ::Subsubsection::Closed:: *)
(*ApodizeEcho*)


Options[ApodizeEcho] = {ApodizationFunction -> "Hanning"}

SyntaxInformation[ApodizeEcho] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeEcho[echo_, OptionsPattern[]] := ApodizeFun[Length[echo], OptionValue[ApodizationFunction], "Echo"] echo


(* ::Subsubsection::Closed:: *)
(*ApodizeFun*)


ApodizeFun[length_, apM_:"Hanning", type_ : "Fid"] := ApodizeFun[length, apM, type] = Block[{app, xdat, xmax},
	xdat = Switch[type,
		"Fid", Range[0, length - 1],
		"Echo", Abs[Round[Range[-length/2, length/2 - 1]]]
	];
	xmax = Max[Abs[xdat]];
	app = Switch[apM,
		"Hanning", 0.5 + 0.5 Cos[xdat Pi/xmax],
		"Hamming", 0.54 + 0.46 Cos[xdat Pi/xmax],
		"Gaussian", Exp[-(3./xmax) xdat],
		"Lorentzian", Exp[-(2./xmax)^2 xdat^2],
		"Voigt", 0.5 Exp[-(3./xmax) xdat] + 0.5 Exp[-(2./xmax)^2 xdat^2]
	];
	app = app/Max[app];
	app
]


(* ::Subsubsection::Closed:: *)
(*ApodizePadFid*)


Options[ApodizePadFid] = {ApodizationFunction -> "Hanning", PaddingFactor -> 2}

SyntaxInformation[ApodizePadFid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadFid[fid_, OptionsPattern[]] := PadFid[ApodizeFid[fid, ApodizationFunction->OptionValue[ApodizationFunction]], PaddingFactor->OptionValue[PaddingFactor]]


(* ::Subsubsection::Closed:: *)
(*ApodizePadEcho*)


Options[ApodizePadEcho] = {ApodizationFunction -> "Hanning", PaddingFactor -> 2}

SyntaxInformation[ApodizePadEcho] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadEcho[echo_, OptionsPattern[]] := PadEcho[ApodizeEcho[echo, ApodizationFunction->OptionValue[ApodizationFunction]],PaddingFactor->OptionValue[PaddingFactor]]


(* ::Subsubsection::Closed:: *)
(*PadSpectra*)


Options[PadSpectra] = {PaddingFactor -> 2, ReadoutType -> "Fid"}

SyntaxInformation[PadSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

PadSpectra[spec_, OptionsPattern[]] := Block[{func,readout},
	readout = OptionValue[ReadoutType];
	func = Switch[readout, "Fid", PadFid, "Echo", PadEcho];
	ShiftedFourier[func[ShiftedInverseFourier[spec, readout],PaddingFactor-> OptionValue[PaddingFactor]], readout]
]


(* ::Subsubsection::Closed:: *)
(*ApodizeSpectra*)


Options[ApodizeSpectra] = {ApodizationFunction -> "Hanning", ReadoutType -> "Fid"}

SyntaxInformation[ApodizeSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeSpectra[spec_, OptionsPattern[]] :=  Block[{func,readout},
	readout = OptionValue[ReadoutType];
	func = Switch[readout, "Fid", ApodizeFid, "Echo", ApodizeEcho];
	ShiftedFourier[func[ShiftedInverseFourier[spec, readout], ApodizationFunction-> OptionValue[ApodizationFunction]], readout]
]


(* ::Subsubsection::Closed:: *)
(*ApodizePadSpectra*)


Options[ApodizePadSpectra] = {ApodizationFunction -> "Hanning", PaddingFactor -> 2, ReadoutType -> "Fid"}

SyntaxInformation[ApodizePadSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadSpectra[spec_, OptionsPattern[]] := Block[{func,readout},
	readout = OptionValue[ReadoutType];
	func = Switch[readout, "Fid", ApodizePadFid, "Echo", ApodizePadEcho];
	ShiftedFourier[func[ShiftedInverseFourier[spec, readout], ApodizationFunction->OptionValue[ApodizationFunction], PaddingFactor->OptionValue[PaddingFactor]],readout]
]


(* ::Subsection:: *)
(*Time and ppm*)


(* ::Subsubsection::Closed:: *)
(*GetTimePpmRange*)


SyntaxInformation[GetTimePpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetTimePpmRange[spec_, {dt_, field_, nuc_}] := GetTimePpmRange[spec, dt, GetGyro[nuc, field]]

GetTimePpmRange[spec_, dt_, field_, nuc_] := GetTimePpmRange[spec, dt, GetGyro[nuc, field]]

GetTimePpmRange[spec_, dt_, gyro_] := GetTimePpmRange[spec, {dt, 0}, gyro]

GetTimePpmRange[spec_, {dt_, te_}, gyro_] := {GetTimeRange[spec, {dt, te}], GetPpmRange[spec, dt, gyro]}


(* ::Subsubsection::Closed:: *)
(*GetPpmRange*)


SyntaxInformation[GetPpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetPpmRange[spec_?VectorQ, {dt_, field_, nuc_}] := GetPpmRange[Length[spec], dt, GetGyro[nuc, field]]

GetPpmRange[spec_?VectorQ, dt_, field_, nuc_] := GetPpmRange[Length[spec], dt, GetGyro[nuc, field]]

GetPpmRange[spec_?VectorQ, dt_, gyro_] := GetPpmRange[Length[spec], dt, gyro]


GetPpmRange[len_?IntegerQ, {dt_, field_, nuc_}] := GetPpmRange[len, dt, GetGyro[nuc, field]]

GetPpmRange[len_?IntegerQ, dt_, field_, nuc_] := GetPpmRange[len, dt, GetGyro[nuc, field]]

GetPpmRange[len_?IntegerQ, dt_, gyro_] := Block[{ppmBw = 1./(dt gyro)}, Reverse@Range[-ppmBw/2, ppmBw/2, ppmBw/(len - 1)]]


(* ::Subsubsection::Closed:: *)
(*GetTimeRange*)


SyntaxInformation[GetTimeRange] = {"ArgumentsPattern" -> {_, _}};

GetTimeRange[fid_?VectorQ, dt_] :=GetTimeRange[fid, {dt, 0}]

GetTimeRange[fid_?VectorQ, {dt_, te_}] := GetTimeRange[Length[fid], {dt, te}]

GetTimeRange[len_?IntegerQ, dt_]:=GetTimeRange[len, {dt, 0.}]

GetTimeRange[len_?IntegerQ, {dt_, te_}] := dt (Range[len]-1) + te


(* ::Subsubsection::Closed:: *)
(*GetGyro*)


SyntaxInformation[GetGyro] = {"ArgumentsPattern" -> {_, _}};

GetGyro[nuc_?StringQ, field_?NumberQ] := GyromagneticRatio[nuc] field

GetGyro[field_?NumberQ,nuc_?StringQ] := GyromagneticRatio[nuc] field


(* ::Subsubsection::Closed:: *)
(*PhaseShiftSpectra*)


SyntaxInformation[PhaseShiftSpectra] = {"ArgumentsPattern" -> {_, _, _., _.}}

PhaseShiftSpectra[spec_, phi0_] := PhaseShiftSpectraC0[spec, phi0];

PhaseShiftSpectra[spec_, ppm_, gyro_, phi1_] := PhaseShiftSpectraC[spec, ppm, gyro, 0., phi1];

PhaseShiftSpectra[spec_, ppm_, gyro_, {phi0_, phi1_}] := PhaseShiftSpectraC[spec, ppm, gyro, phi0, phi1];


PhaseShiftSpectraC0 = Compile[{{spec, _Complex, 1}, {phi0, _Real, 0}},
	Exp[-I phi0] spec,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]

PhaseShiftSpectraC = Compile[{{spec, _Complex, 1}, {ppm, _Real, 1}, {gyro, _Real, 0}, {phi0, _Real, 0}, {phi1, _Real, 0}},
	Exp[-I (phi0 + 2 Pi (phi1 / 1000) gyro ppm)] spec,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


(* ::Subsubsection::Closed:: *)
(*PhaseShiftFid*)


SyntaxInformation[PhaseShiftFid] = {"ArgumentsPattern" -> {_, _, _., _.}};

PhaseShiftFid[fid_, phi0_] := PhaseShiftFidC0[fid, phi0];

PhaseShiftFid[fid_, time_, phi1_] := PhaseShiftFidC[fid, time, 0., phi1];

PhaseShiftFid[fid_, time_, {phi0_, phi1_}] := PhaseShiftFidC[fid, time, phi0, phi1];

PhaseShiftFidC0 = Compile[{{fid, _Complex, 1}, {phi0, _Real, 0}}, 
	Exp[-I phi0] fid, 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

PhaseShiftFidC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {phi0, _Real, 0}, {phi1, _Real, 0}}, 
	Exp[-I (phi0 + 2 Pi phi1 time)] fid, 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*TimeShiftFid*)


SyntaxInformation[TimeShiftFid] = {"ArgumentsPattern" -> {_, _, _., _.}}

TimeShiftFid[fid_, time_, gam_] := TimeShiftFidC[fid, time, 0, gam, gam, 0.];

TimeShiftFid[fid_, time_, {gamL_, gamG_}] := TimeShiftFidC[fid, time, 0., gamL, gamG, 0.];

TimeShiftFid[fid_, time_, gyro_, {gam_, eps_}] := TimeShiftFidC[fid, time, gyro, gam, gam, eps];

TimeShiftFid[fid_, time_, gyro_, {{gamL_, gamG_}, eps_}] := TimeShiftFidC[fid, time, gyro, gamL, gamG, eps];

TimeShiftFidC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {gamL, _Real, 0}, {gamG, _Real, 0}, {eps, _Real, 0}},
	Exp[- Pi (gamL time + 2 Pi gamG^2 time^2) + 2 Pi eps gyro I time] fid, 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


(* ::Subsubsection::Closed:: *)
(*TimeShiftEcho*)


SyntaxInformation[TimeShiftEcho] = {"ArgumentsPattern" -> {_, _, _., _.}}

TimeShiftEcho[fid_, time_, gam_] := TimeShiftEchoC[fid, time, 0., gam, 0., .5];

TimeShiftEcho[fid_, time_, {gam_, f_}] := TimeShiftEchoC[fid, time, 0., gam, 0., f];

TimeShiftEcho[fid_, time_, gyro_, {gam_, eps_}] := TimeShiftEchoC[fid, time, gyro, gam, eps, .5];

TimeShiftEcho[fid_, time_, gyro_, {gam_, eps_, f_}] := TimeShiftEchoC[fid, time, gyro, gam, eps, f];

TimeShiftEchoC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {gam, _Real, 0}, {eps, _Real, 0}, {f, _Real, 0}},Block[{timeNew},
		timeNew = time - (time[[-1]]/2);
		Exp[- Pi (gamL time + 2 Pi gamG^2 time^2) + 2 Pi eps gyro I time] fid
	],RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


(* ::Subsubsection::Closed:: *)
(*ShiftSpectra*)


Options[ShiftSpectra] = {ReadoutType -> "Fid"}

SyntaxInformation[ShiftSpectra] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}}

ShiftSpectra[spec_, {dw_, gyro_}, shift_, OptionsPattern[]] := Block[{readout, func},
	readout = OptionValue[ReadoutType]; 
	func = If[readout==="Fid",ShiftFidC,ShiftEchoC];
	ShiftedFourier[func[ShiftedInverseFourier[spec,readout], GetTimeRange[spec, dw], gyro, shift],readout]
]

ShiftFidC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {eps, _Real, 0}}, 
	Exp[2 Pi eps gyro I time] fid,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]

ShiftEchoC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {eps, _Real, 0}}, Block[{timeNew},
	timeNew = time - (time[[-1]]/2);
	Exp[2 Pi eps gyro I timeNew] fid
	],RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


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
(*ChangeDwellTimeFid*)


SyntaxInformation[ChangeDwellTimeFid] = {"ArgumentsPattern" -> {_, _, _}}

ChangeDwellTimeFid[time_, dwOrig_, dwTar_] := Block[{timeOrig, timeTar, sc},
	sc = dwOrig / dwTar;
	(*get time of original signal*)
	timeOrig = dwOrig(Range[Round[Length@time]]-1);
	timeTar = dwTar(Range[Round[sc Length@time]]-1);
	(*Interpolate the time to the new timescale*)
	Interpolation[Transpose[{timeOrig, time}], InterpolationOrder -> 1,"ExtrapolationHandler"->{(0.0 &),"WarningMessage"->False}][timeTar]
]


(* ::Subsection::Closed:: *)
(*GetBasisFunctions*)


Options[GetSpectraBasisFunctions] = {
	BasisSequence -> {"PulseAcquire", 0},
	SpectraSamples -> 2046,
	SpectraBandwidth -> 2000,
	SpectraNucleus -> "1H",
	SpectraPpmShift -> 4.65,
	SpectraFieldStrength -> 3
};

SyntaxInformation[GetSpectraBasisFunctions] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

GetSpectraBasisFunctions[inp_, opts:OptionsPattern[]] := GetSpectraBasisFunctions[inp, {""}, opts]

GetSpectraBasisFunctions[inp_, split_, OptionsPattern[]] := Block[{
		cf, bw, nSamp, lw, lws, read, sys, din, struct, dr, te, tm, sVals, seq, readout, 
		t1, t2, nEcho, sysAll, spinTab, field, nuc, nam, labs, names, 
		times, fids, specs
	},

	(*get the option values*)
	{seq, sVals} = OptionValue[BasisSequence];
	readout = Switch[seq, "PulseAcquire"|"STEAM","Fid", "SpaceEcho","Echo"];

	nSamp = OptionValue[SpectraSamples];
	bw = OptionValue[SpectraBandwidth];
	cf = OptionValue[SpectraPpmShift];
	field = OptionValue[SpectraFieldStrength];
	nuc = OptionValue[SpectraNucleus];

	(*get or make the spin systems*)
	sysAll = If[StringQ[#], GetSpinSystem, MakeSpinSystem][#, CenterFrequency -> cf] & /@ inp;
	sysAll = Reverse@SortBy[sysAll,First[#[[2]]] &];
	spinTab = SysTable[sysAll];

	(*loop over systems to perform simulation of basis functions*)
	{names, fids, specs} = Transpose[Flatten[
		Map[(
			(*get the system*)
			sys = #;
			nam = sys[[7]];
			labs = sys[[5]];

			(*see if molecule has to be split into individual peaks*)
			read = If[MemberQ[split, sys[[-1]]], "each", "all"];

			(*generate hamiltonian*)
			{din, struct} = SimHamiltonian[sys, FieldStrength -> field, SimNucleus -> nuc];

			(*perform sequence*)
			Switch[seq,
				"PulseAcquire",
				te = sVals;
				dr = SequencePulseAcquire[din, struct, te],
				"STEAM",
				{te, tm} = sVals;
				dr = SequenceSteam[din, struct, {te, tm}],
				"SpaceEcho",
				{t1, t2, nEcho} = sVals;
				dr = SequenceSpaceEcho[din, struct, t1, t2, nEcho, 1]
			];

			(*simulate readout*)
			fids = 4 First@SimReadout[dr, struct, 
				ReadoutSamples -> nSamp, ReadoutBandwidth -> bw, CenterFrequency -> cf, 
				Linewidth -> 0, ReadoutOutput -> read
			];

			(*check if multiple fids for metabolite and act accordingly*)
			If[VectorQ[fids],
				specs = ShiftedFourier[fids, readout];
				{{nam, fids, specs}}
				,
				(*multiple fids*)
				MapIndexed[(
					fids = (Re[#1] + Im[#1] I);
					specs = ShiftedFourier[fids, readout];
					{nam <> "-" <> labs[[#2]], fids, specs}
				) &, fids]
			]

		(*close mapping function*)
		) &, sysAll]
	, 1]];

	{names, fids, specs, spinTab}
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
, RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}]

ConFuncUC = Compile[{{par, _Real, 0}, {min, _Real, 0}, {sc, _Real, 0}}, Block[{off = -par + min},
	10^sc (UnitStep[off] (off))^2
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}]

ConFuncLC = Compile[{{par, _Real, 0}, {max, _Real, 0}, {sc, _Real, 0}}, Block[{off = par - max},
	10^sc (UnitStep[off] (off))^2
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}]

ConFuncC = Compile[{{par, _Real, 0}, {min, _Real, 0}, {max, _Real, 0}, {sc, _Real, 0}}, Block[{off1 = -par + min, off2 = par - max},
	10^sc (UnitStep[off1] (off1)^2 + UnitStep[off2] (off2)^2)
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}]


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
],RuntimeOptions->"Speed"];


(* ::Subsection::Closed:: *)
(*PlotSpectra*)


Options[PlotSpectra] = {
	PlotRange -> Full, 
	Method -> "All", 
	GridLines -> {},
	PlotColor -> Automatic,
	GridLineSpacing -> 1, 
	SpectraSpacing -> 0.2, 
	PlotLabels -> None,
	AspectRatio -> .2, 
	ImageSize -> 750, 
	PlotLabel -> None,
	CenterFrequency -> 0.,
	Filling->False
};

SyntaxInformation[PlotSpectra] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

PlotSpectra[spec_, {dwell_?NumberQ, gyro_?NumberQ}, opts : OptionsPattern[]] := PlotSpectra[GetPpmRange[If[MatrixQ[spec],spec[[1]],spec], dwell, gyro], spec, opts]

PlotSpectra[spec_, {dwell_?NumberQ, field_?NumberQ, nuc_?StringQ}, opts : OptionsPattern[]] := PlotSpectra[GetPpmRange[If[MatrixQ[spec],spec[[1]],spec], dwell, field, nuc], spec, opts]

PlotSpectra[ppm_?VectorQ, spec_, OptionsPattern[]] := Block[{
		fun, plot, plot2, grid, gridS, or, rr, col, space, cols, cols2, pl1, pl2, labels, min, max, shift
	},

	(*get the plot range*)
	rr = OptionValue[PlotRange] /. Automatic -> Full;
	Switch[rr,
		{_, Full}, rr[[2]] = {Min[{Re[spec], Im[spec]}], Max[Abs[spec]]},
		Full, rr = {Full, {Min[{Re[spec], Im[spec]}], Max[Abs[spec]]}}
	];

	(*get grid line options*)
	gridS = OptionValue[GridLineSpacing];
	{min,max} = If[rr[[1]]===Full,Round[MinMax[ppm]],MinMax[rr[[1]]]];

	grid = Sort@DeleteDuplicates@Join[
		If[gridS === 0, {}, Range[0, Round[max], gridS]],
		If[gridS === 0, {}, -Range[0, -Round[min], gridS]],
		OptionValue[GridLines]
	];
	shift = OptionValue[CenterFrequency];

	(*if only one spectra plot normal else plot expanded list of spectra*)
	If[VectorQ[spec],
		(*plot single spectra*)
		(*get the plot functions*)
		fun = Switch[OptionValue[Method], "Abs", {Abs}, "Re", {Re}, "Im", {Im}, "ReIm", {Im, Re}, "All", {Im, Re, Abs}];
		(*plot single spectra*)
		plot = Transpose[{ppm + shift, #}] & /@ (#@spec & /@ fun);
		(*get the plot color*)
		col = If[OptionValue[PlotColor] === Automatic, 
			({{LightDarkV[Gray, White], AbsoluteThickness[1]}, {Darker[Red], AbsoluteThickness[1]}, {LightDarkV[]}}[[-Length[fun] ;;]]), 
			OptionValue[PlotColor]
		];

		(*Make the plot*)
		ListLinePlot[plot, PlotStyle -> col, PlotRange -> rr, GridLines -> {grid, {0}}, 
			AspectRatio -> OptionValue[AspectRatio], ImageSize -> OptionValue[ImageSize], 
			PlotLabel -> OptionValue[PlotLabel], Frame -> {{False, False}, {True, False}}, 
			FrameStyle -> Directive[{Thick, LightDarkV[]}], FrameLabel -> {"PPM", None},
			LabelStyle -> {Bold, 14, LightDarkV[]}, PlotHighlighting -> False,
			PerformanceGoal->"Speed", MaxPlotPoints->Infinity, Background->None, 
			Filling->OptionValue[Filling],
			ScalingFunctions -> {"Reverse", Automatic}
		]

		,
		(*plot List of spectra*)
		(*get the plot functions*)
		fun = Switch[OptionValue[Method], "All", Abs, "Abs", Abs, "Re", Re, "ReIm", Re, "Im", Im, _, Return[]];

		(*space the spectra over the y axes*)
		space = Reverse@Range[0, Length[spec]] Max[Abs[spec]] OptionValue[SpectraSpacing];

		(*plot the spectra*)
		plot = Transpose[{ppm + shift, #}] & /@ (fun[Append[spec, Total[spec]]] + space);

		(*correct the plot range*)
		If[rr[[2]] =!= Full, rr[[2, 2]] = 1.1 Max[plot[[All, All, 2]]]];
		If[rr[[2]] =!= Full, rr[[2, 1]] = Min[plot[[All, All, 2]]] - .1 Max[plot[[All, All, 2]]]];

		If[rr[[1]]=!=Full, plot=Select[#,(Min[rr[[1]]]<=#[[1]]<=Max[rr[[1]]])&]&/@plot];

		(*get the plot colors*)
		cols = Thread[{Append[ConstantArray[LightDarkV[], Length[plot] - 1], Red], Thick}];
		labels = OptionValue[PlotLabels];

		(*make the plot*)
		pl1 = ListLinePlot[plot, Frame -> {{False, False}, {True, False}}, FrameStyle -> Thickness[.003], FrameTicksStyle -> Thickness[.003],
			PlotRange -> rr, PlotRangeClipping -> True, PlotStyle -> cols, ScalingFunctions -> {"Reverse", Automatic}, PlotLabel -> OptionValue[PlotLabel], 
			PlotLabels ->If[(OptionValue[Method] === "ReIm") || (labels === None), None, (Style[#, LightDarkV[], Bold, 14] & /@ Append[labels, "All"])],
			GridLines -> {grid, None}, AspectRatio -> .5, ImageSize -> 1000, FrameLabel -> {"PPM", None}, LabelStyle -> {Bold, 14, LightDarkV[]}, 
			PerformanceGoal->"Speed", MaxPlotPoints->Infinity, PlotHighlighting -> False
		];

		(*make Im plot if needed*)
		If[OptionValue[Method] === "ReIm",
			plot2 = Transpose[{ppm + shift, #}] & /@ (Im[Append[spec, Total[spec]]] + space);
			cols2 = {Append[ConstantArray[Gray, Length[plot] - 1], Gray]};

			pl2 = ListLinePlot[plot2, Frame -> {{False, False}, {True, False}}, FrameStyle -> Thickness[.003], FrameTicksStyle -> Thickness[.003], 
				PlotRange -> rr, PlotRangeClipping -> True, PlotStyle -> cols2, ScalingFunctions -> {"Reverse", Automatic},
				PlotLabels -> (Style[#, LightDarkV[], Bold, 14] & /@ Append[OptionValue[PlotLabels], "All"]), PlotLabel -> OptionValue[PlotLabel],
				GridLines -> {grid, None}, PlotRange -> rr, AspectRatio -> .5, 
				ImageSize -> 1000, FrameLabel -> {"PPM", None}, MaxPlotPoints->Infinity,
				LabelStyle -> {Bold, 14, LightDarkV[]}, PerformanceGoal->"Speed", PlotHighlighting -> False
			];

			Show[pl2, pl1]
			,
			pl1
		]
	]
]


(* ::Subsection::Closed:: *)
(*PlotFid*)


Options[PlotFid] = {
	PlotRange -> Full, 
	Method -> "All", 
	GridLines -> {}, 
	PlotColor -> Automatic, 
	GridLineSpacing -> 1, 
	AspectRatio -> 0.2, 
	ImageSize -> 750, 
	PlotLabel -> None
};

SyntaxInformation[PlotFid] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

PlotFid[fid_?VectorQ, dwell_?NumberQ, opts : OptionsPattern[]] := PlotFid[
	GetTimeRange[fid, dwell], fid, opts]

PlotFid[time_?VectorQ, fid_?VectorQ, OptionsPattern[]] := Block[{fun, plot, grid, gridS, rr, col},
	gridS = OptionValue[GridLineSpacing];
	grid = DeleteDuplicates@Join[If[gridS === 0, {},Range[0, Max[time], gridS]],OptionValue[GridLines]];

	fun = Switch[OptionValue[Method], "Abs", {Abs}, "Re", {Re}, "Im", {Im}, "ReIm", {Im, Re}, "All", {Im, Re, Abs}];

	plot = Transpose[{time, #}] & /@ (#@fid & /@ fun);

	rr = OptionValue[PlotRange];
	Switch[rr,
		{_, Full}, rr[[2]] = {-Max[Abs[fid]], Max[Abs[fid]]},
		Full, rr = {Full, {-Max[Abs[fid]], Max[Abs[fid]]}}
	];
	col = If[OptionValue[PlotColor] === Automatic, 
		{
			{LightDarkV[Gray, White], AbsoluteThickness[1]}, 
			{Darker[Red], AbsoluteThickness[1]}, 
			{LightDarkV[]}
		}[[-Length[fun] ;;]], 
		OptionValue[PlotColor]];

	ListLinePlot[plot, PlotStyle -> col, PlotRange -> rr, GridLines -> {grid, {0.}}, 
		AspectRatio -> OptionValue[AspectRatio], ImageSize -> OptionValue[ImageSize], 
		PlotLabel -> OptionValue[PlotLabel], Frame -> {{False, False}, {True, False}},
		FrameStyle -> Directive[{Thick, LightDarkV[]}], FrameLabel -> {"time [s]", None},
		LabelStyle -> {Bold, 14, LightDarkV[]}, PlotHighlighting -> False,
		PerformanceGoal->"Speed", MaxPlotPoints->Infinity, Background->None
	]
]


(* ::Subsection::Closed:: *)
(*PlotCSIData*)


Options[PlotCSIData] = {PlotRange -> Full};

SyntaxInformation[PlotCSIData] = {"ArgumentsPattern" -> {_, _, _., _., OptionsPattern[]}};

PlotCSIData[dataInp_, dw_?NumberQ, field_?NumberQ, nuc_?StringQ, opts : OptionsPattern[]] := PlotCSIData[dataInp, {dw, GetGyro[nuc, field]}, opts];

PlotCSIData[dataInp_, {dw_?NumberQ, field_?NumberQ, nuc_?StringQ}, opts : OptionsPattern[]] := PlotCSIData[dataInp, {dw, GetGyro[nuc, field]}, opts];

PlotCSIData[dataInp_, dw_?NumberQ, gyro_?NumberQ, opts : OptionsPattern[]] := PlotCSIData[dataInp, {dw, gyro}, opts];

PlotCSIData[dataInp_, {dw_?NumberQ, gyro_?NumberQ}, OptionsPattern[]] := DynamicModule[{
		data, dim, nMax, xdat, xmin, xmax, xdatR, datN, sp, coors, maxCol, totCol, 
		lab, specP, size, ccs, dataPlot, backColM, backColT
	},

	data = N@dataInp;

	If[! ListQ[data], Return[];];

	dim = Dimensions[data];
	nMax = dim[[1]];

	xdat = GetPpmRange[data[[1, 1, 1]], dw, gyro];
	{xmin, xmax} = MinMax[xdat];
	xdatR = Reverse[Rescale[xdat[[;; ;; 4]], {xmin, xmax}, {0.1, 0.9}]];

	datN = Map[(sp = Abs[#1[[;; ;; 4]]]; {xdatR, Rescale[sp, MinMax[sp], {0.1, .9}] 1}) &, data, {-2}];

	coors = {
		Table[{j, i}, {i, dim[[2]]}, {j, dim[[3]]}],
		Table[{j, i}, {i, dim[[1]]}, {j, dim[[3]]}],
		Table[{j, i}, {i, dim[[1]]}, {j, dim[[2]]}]
	};

	maxCol = Map[GrayLevel, Rescale[Map[Max, Abs[data], {-2}]], {3}];
	totCol = Map[GrayLevel, Rescale[Map[Total, Abs[data], {-2}]], {3}];

	pan = Manipulate[
		nMax = dim[[or]];
		sz = {dx, dy} = Most[Drop[dim, {or}]];

		size = Reverse[sz 25];

		Column[{
			Dynamic[
				lab = Switch[or, 1, {n, -coor[[1]] + dx + 1, coor[[2]]}, 2, {coor[[1]], n, coor[[2]]}, 3, {coor[[1]], coor[[2]], n}];
				specP = PhaseShiftSpectra[
						Switch[or, 1, Reverse[data][[n, coor[[1]], coor[[2]]]], 2, data[[coor[[1]], n, coor[[2]]]], 3, data[[coor[[1]], coor[[2]], n]]]
					, xdat, gyro, {ph, te}];
				
				FlipView[{PlotSpectra[If[app, ApodizePadSpectra[specP], specP], {gyro, dw},
						PlotRange -> {{pMin, pMax}, Full}, Method -> functions, 
						PlotLabel -> lab,ImageSize -> .5 First@size, AspectRatio -> 0.5
					],
					PlotFid[ShiftedInverseFourier[specP], dw, Method -> functions, 
						PlotLabel -> lab, ImageSize -> .5 First@size, AspectRatio -> 0.5]
				}]
			, TrackedSymbols :> {n, or, ph, te, app, pMin, pMax, functions, coor}]
			,
			Dynamic[
				ccs = Flatten[coors[[or]], 1];
				dataPlot = 
				Flatten[Switch[or, 1, Reverse@datN[[n]], 2, datN[[All, n]], 3, 
					datN[[All, All, n]]], 1];
				backColM = 
				Flatten[Switch[or, 1, Reverse@maxCol[[n]], 2, maxCol[[All, n]], 
					3, maxCol[[All, All, n]]], 1];
				backColT = 
				Flatten[Switch[or, 1, Reverse@totCol[[n]], 2, totCol[[All, n]], 
					3, totCol[[All, All, n]]], 1];
				
				EventHandler[
					Deploy@Graphics[{
							Flatten[
								With[{idx = Reverse[#[[1]]]}, {
									EdgeForm[If[coor == idx, Directive[Thick, Red], None]],
									#[[2]], Rectangle[#[[1]] + 0.05, #[[1]] + 0.95]
									}
								] & /@ Thread[{ccs, Switch[backScale, "Max", backColM, "Total", backColT]}]
							],
							{StandardBlue, Line[Transpose /@ (dataPlot + ccs)]}
						},
						ImageSize -> size, Axes -> True, AxesOrigin -> {0, 0}
					],
					{"MouseDown" :> With[{p = MousePosition["Graphics"]}, If[p =!= None, coor = Reverse[Floor[p]]]]}
				]
			, TrackedSymbols :> {n, or, backScale, coor}]
			}, Alignment -> Center]

		, {{or, 1, "Orientation"},
			{1 -> "Transversal", 2 -> "Coronal", 3 -> "Sagittal"},
			ControlType -> SetterBar,
			TrackingFunction -> (With[{o = #},
				or = o;
				n = Min[n, dim[[o]]];
				
				coor = MapThread[
				Min, {coor, {dim[[{2, 3}]], dim[[{1, 3}]], 
					dim[[{1, 2}]]}[[o]]}];
				] &)}
		, {{n, Ceiling[Length[data]/2], "Slice"}, 1, Dynamic[nMax], 1}
		, {{backScale, "Max", "Background"}, {"Max", "Total"}}
		, {{app, True, "Apodize and Pad"}, {True, False}}
		, {{ph, 0, "Phase spectra"}, -2 Pi, 2 Pi}
		, {{te, 0, "EchoTime"}, 0, 1, .1}
		, Delimiter
		, {{functions, "ReIm", "Spectra function"}, {"Re", "Im", "ReIm", "Abs", "All"}}
		, {{pMin, xmin, "Min pmm"}, xmin, Dynamic[pMax - 1]}
		, {{pMax, xmax, "Min pmm"}, Dynamic[pMin + 1], xmax}

		, {{coor, {1, 1}}, ControlType -> None}
		,
		ControlPlacement->Right,
		SynchronousInitialization -> False
	];

	NotebookClose[plotWindow];
	plotWindow = CreateWindow[DialogNotebook[{CancelButton["Close", Clear[data]; DialogReturn[]], pan}, 
		WindowSize -> All, WindowTitle -> "Plot data window"]];
];


(* ::Subsection:: *)
(*SpectraSimulator*)


SpectraSimulator[]:=DynamicModule[{
		use, vals, names, namD, nSamp, bwS, fieldS, nuc, nucI, dw, gyro, field, base, grids, ppmI
	},

	vals = {
		{100, 50, 20, 10}, 
		{20, 30, 40, 100, 0, 75, 500, 250, 40, 60}, 
		{100}
	};
	use = (0 vals + 1) True;

	names = {{"HDO","GlucoseD","GlutamateD","LactateD"}, {"PE","PC","Piex","Piin","GPE","GPC","PCr","ATP","NAD","UDPG"}, {"Na"}};
	namD = {{"Wat","Gluc","Glx","Lact"}, {"PE","PC","Piex","Piin","GPE","GPC","PCr","ATP","NAD","UDPG"}, {"Na"}};
	

	grids={
		{
			{Row[{namD[[1,1]]<>": ",InputField[Dynamic[vals[[1,1]]],FieldSize->3],Checkbox[Dynamic[use[[1,1]]]]}],Row[{namD[[1,2]]<>": ",InputField[Dynamic[vals[[1,2]]],FieldSize->3],Checkbox[Dynamic[use[[1,2]]]]}],
			Row[{namD[[1,3]]<>": ",InputField[Dynamic[vals[[1,3]]],FieldSize->3],Checkbox[Dynamic[use[[1,3]]]]}],Row[{namD[[1,4]]<>": ",InputField[Dynamic[vals[[1,4]]],FieldSize->3],Checkbox[Dynamic[use[[1,4]]]]}]},
			{}
		},{
			{Row[{names[[2,1]]<>": ",InputField[Dynamic[vals[[2,1]]],FieldSize->3],Checkbox[Dynamic[use[[2,1]]]]}],Row[{names[[2,2]]<>": ",InputField[Dynamic[vals[[2,2]]],FieldSize->3],Checkbox[Dynamic[use[[2,2]]]]}],
			Row[{names[[2,3]]<>": ",InputField[Dynamic[vals[[2,3]]],FieldSize->3],Checkbox[Dynamic[use[[2,3]]]]}],Row[{names[[2,4]]<>": ",InputField[Dynamic[vals[[2,4]]],FieldSize->3],Checkbox[Dynamic[use[[2,4]]]]}],
			Row[{names[[2,5]]<>": ",InputField[Dynamic[vals[[2,5]]],FieldSize->3],Checkbox[Dynamic[use[[2,5]]]]}]},
			{Row[{names[[2,6]]<>": ",InputField[Dynamic[vals[[2,6]]],FieldSize->3],Checkbox[Dynamic[use[[2,6]]]]}],
			Row[{names[[2,7]]<>": ",InputField[Dynamic[vals[[2,7]]],FieldSize->3],Checkbox[Dynamic[use[[2,7]]]]}],Row[{names[[2,8]]<>": ",InputField[Dynamic[vals[[2,8]]],FieldSize->3],Checkbox[Dynamic[use[[2,8]]]]}],
			Row[{names[[2,9]]<>": ",InputField[Dynamic[vals[[2,9]]],FieldSize->3],Checkbox[Dynamic[use[[2,9]]]]}],Row[{names[[2,10]]<>": ",InputField[Dynamic[vals[[2,10]]],FieldSize->3],Checkbox[Dynamic[use[[2,10]]]]}]}
		},{
			{Row[{names[[3,1]]<>": ",InputField[Dynamic[vals[[3,1]]],FieldSize->3],Checkbox[Dynamic[use[[3,1]]]]}]},
			{""}
		}
	};

	{nSamp, bwS, fieldS, nucI} = {8096, 10000, 7, 2};
	
	nuc = {"2H", "31P", "23Na"}[[nucI]];
	{dw, gyro, field, base} = simulatedFid[{nSamp, bwS, fieldS, nuc}, names[[nucI]]];
	ppmI = {{7.5, -2.5}, {8, -18}, {10, -10}}[[nucI]];

	man = Manipulate[
		simTrigger;

		(*see which metabolites to use and with which amplitude*)
		useI = If[AllTrue[use[[nucI]], #===False&], (0 vals[[nucI]] + 1) True, use[[nucI]]];
		valI = Replace[vals[[nucI]], x_ /; ! NumberQ[x] :> 0, {1}];
		peakSel = Pick[Thread[{names[[nucI]] , valI}], useI];
		fidsT = peakSel[[All, 2]] . (peakSel[[All, 1]] /. base);

		(*resample fid and generate noise*)
		dw2 = 1./bw;
		ni = Round[(te / 1000) / dw];
		ns = Round[Min[{ns, nSamp (dw / dw2) - ni}]];
		fid1 = ChangeDwellTimeFid[fidsT[[ni+1;;]], dw, dw2][[;;ns]];

		(* Define reference values *)
		bwRef = 5000.;
		nRef = 512.;
		(* Calculate scaled sigma *)
		(* We use the Max of the original fidsT to keep the SNR relative to the source signal *)
		sig = (First[Abs[fidsT]] / snr) * Sqrt[bw / bwRef] * Sqrt[nRef / ns] / Sqrt[2];

		(*recalculate noise only if needed*)
		If[lFid =!= Length[fid1] || sigI =!= sig,
			lFid = Length[fid1];
			sigI = sig;
			noise1 = Complex@@@RandomReal[NormalDistribution[0., sig], {lFid, 2}]
		];
		
		(*get the time of the acquired signal for manipulations*)
		{timeI, ppmI} = GetTimePpmRange[fid1, {dw2, field, nuc}];

		(*add shift and linewidth fo the FID*)
		fid2 = If[timeShift, 
			PhaseShiftFid[TimeShiftFid[fid1, timeI + te / 1000, gyro, {{gamma, sigma}, eps}], ph0], 
			fid1
		];
		lw = effectiveLinewidth[gamma, sigma, gyro];

		(*add phase to spectra as post processing and visualization*)
		fid3 = ShiftedInverseFourier[PhaseShiftSpectra[ShiftedFourier[fid2], ppmI, gyro, {ph0s, ph1s}]];

		(*apply apodization and padding and figure out time of spectra to plot for visualization*)
		{fid4, noise4} = If[ap==="None", {fid3, noise1}, 
			ApodizeFid[#, ApodizationFunction->ap]& /@ {fid3, noise1}];
		{fid5, noise5} = If[pad===1, {fid4, noise4}, PadFid[#, 
			PaddingFactor->pad]& /@ {fid4, noise4}];

		(*prepare for plotting*)
		fidTot = fid5 + noise5;
		noiseF = ShiftedFourier[noise5];
		specF = ShiftedFourier[fid5];
		specTot = noiseF + specF;
		fidP = Switch[show, "both", fidTot, "signal", fid5, "noise", noise5];
		specP = Switch[show, "both", specTot, "signal", specF, "noise", noiseF];

		Grid[{
			{Style["Effective linewidth [Hz]: "<> ToString[Round[lw[[1]]]], LightDarkV[]]},
			{
				timeP = GetTimeRange[fidP, dw2] + (te - ph1s) / 1000;
				maxFid = Max[Abs[fidTot]];
				Show[
					PlotFid[timeP, TimeShiftFid[fidP, timeP, gyro, {0, off}], 
						Method -> met, PlotRange -> {{-0.001, Max[timeP] + 0.001}, {-1.2, 1.2} maxFid}, 
						GridLineSpacing -> 0.05, ImageSize -> pSize
					],
					If[!timeShift, Graphics[], ListLinePlot[Transpose@{timeP, maxFid shiftFun[timeP, eps, gyro]}, PlotStyle->{Darker@Blue, Thin}]],
					If[!timeShift, Graphics[], ListLinePlot[Transpose@{timeP, maxFid relaxFun[timeP, {gamma, sigma}]}, PlotStyle->{Darker@Orange, Dashed}, PlotRange->Full]],
					If[ap==="None", Graphics[], ListLinePlot[Transpose@{timeI, maxFid apodizeFun[timeI, ap]}, PlotStyle->{Darker@Green, Dashed}]]
				]
			},{
				maxSpec = Max[Abs[specTot]];
				PlotSpectra[ShiftSpectra[specP, {dw2, gyro}, off], {dw2, gyro}, 
					PlotRange -> {If[pRan==="Automatic", Full, {pMax, pMin}], {-0.3, 1.2} maxSpec},
					GridLineSpacing -> 5, AspectRatio -> .5, Method -> met, ImageSize -> pSize]
			}
		}]

		,
		Style["Spectra simulation", 14, Bold],
		{{nucleus, 2, "nucleus"}, {1->"2H", 2->"31P", 3->"23Na"}},
		{{fieldS, 7, "field"}, {1, 1.5, 3, 7}},
		Button["Simulate fid for nucleus and field",
			simTrigger = False;
			nucI = nucleus;
			nuc = {"2H", "31P", "23Na"}[[nucI]];
			{dw, gyro, field, base} = simulatedFid[{nSamp, bwS, fieldS, nuc}, names[[nucI]]];
			{pMin, pMax} = {{7.5, -2.5}, {8, -18}, {10, -10}}[[nucI]];
			simTrigger = True;
		, Method->"Queued"],
		
		Delimiter,
		Dynamic[Grid[grids[[nucI]]]],

		Delimiter,
		TabView[{
			Style["acquisition and signal", 12, Bold] -> Column[{
				"",
				Style["acquisition settings", 10, Bold],
				Control[{{off, 0., "F0 [ppm]"}, -10, 10}],
				Control[{{bw, 5000., "bandwidth [Hz]"}, 500, 10000}],
				Control[{{ns, 512, "n samples"}, 32, 1024, 16}],
				Control[{{te, 0., "echo time [ms]"}, 0,5}],
				"",
				Style["FID properties", 10, Bold],
				Control[{{snr, 40., "snr (t = te)"}, 2, 100}],
				Control[{{ph0, 0, "phase [rad]"}, -2 Pi, 2 Pi}],
				Control[{{eps, 0, "shift [ppm]"}, -5, 5}],
				Control[{{gamma , 10, "Lorentzian linewidth [Hz]"}, 0, 100}],
				Control[{{sigma, 0, "Gaussian linewidth [Hz]"}, 0, 100}]
			}]
			,
			Style["plotting",12,Bold] -> Column[{
				"",
				Style["spectra phasing", 10, Bold],
				Control[{{ph0s, 0., "0th order phase [rad]"}, -2 Pi, 2 Pi}],
				Control[{{ph1s, 0., "1st order phase [ms]"}, -5., 5.}],
				"",
				Style["padding and apodization", 10, Bold],
				Control[{{pad,1,"padding factor"}, 1,4, 0.5}],
				Control[{{ap, "None", "apodization kernel"}, 
					{"None", "Hanning", "Hamming", "Gaussian", "Lorentzian", "Voigt"},
					ControlType->SetterBar}],
				"",
				Style["what to plot", 10, Bold],
				Control[{{timeShift, True, "apply linewidth and shift"}, {False, True}}],
				Control[{{show,"both","signal and noise"}, {"both","signal","noise"}}],
				Control[{{met,"ReIm","plot method"}, {"All","Abs","Re","Im","ReIm"}}],
				"",
				Style["options", 10, Bold],
				Control[{{pRan, "Automatic", "ppm plot range"}, {"Automatic","Manual"}}],
				Control[{{pMin, 10, "ppm min"}, Dynamic[pMax], Dynamic[Max[ppmI]]}],
				Control[{{pMax, -20, "ppm max"}, Dynamic[Min[ppmI]], Dynamic[pMin]}],
				Button["reset ppm range", {pMin, pMax} = {{7.5, 2.5}, {10, -20}, {-10,10}}[[nucleus]]],
				Control[{{pSize, 600, "plot size"}, {200, 300, 400, 500, 600, 700}}]
			}]
		}],

		{{bw, 5000.}, ControlType -> None},
		{{ns, 512}, ControlType -> None},
		{{te, 0.}, ControlType -> None},
		{{snr, 40.}, ControlType -> None},

		{{timeShift, True}, ControlType -> None},
		{{ph0, 0}, ControlType -> None},
		{{eps, 0}, ControlType -> None},
		{{gamma, 10}, ControlType -> None},
		{{sigma, 0.}, ControlType -> None},
		{{ph0s, 0.}, ControlType -> None},
		{{ph1s, 0.}, ControlType -> None},

		{{pad, 1}, ControlType -> None},
		{{ap, "None"}, ControlType -> None},
		{{show, "both"}, ControlType -> None},
		{{pRan, "Manual"}, ControlType -> None},
		{{pMin, 8}, ControlType -> None},
		{{pMax, -18}, ControlType -> None},
		{{met, "ReIm"}, ControlType -> None},
		{{pSize, 600}, ControlType -> None},

		(*hidden dynamic parameters*)
		{dw2, ControlType -> None},
		{ni, ControlType -> None},
		{fid1, ControlType -> None},
		{fid2, ControlType -> None},
		{fid3, ControlType -> None},
		{fid4, ControlType -> None},
		{fid5, ControlType -> None},
		{noise1, ControlType -> None},
		{noise4, ControlType -> None},
		{noise5, ControlType -> None},
		{noiseF, ControlType -> None},
		{fidTot, ControlType -> None},
		{fidP, ControlType -> None},
		{specF, ControlType -> None},
		{specTot, ControlType -> None},
		{specP, ControlType -> None},
		{sig, ControlType -> None},
		{n1, ControlType -> None},
		{timeI, ControlType -> None},

		{timeP, ControlType -> None},
		{maxFid, ControlType -> None},
		{maxSpec, ControlType -> None},
		{lFid, ControlType -> None},
		{sigI, ControlType -> None},

		{fidsT, ControlType -> None},
		{valI, ControlType -> None},
		{useI, ControlType -> None},

		{simTrigger, ControlType -> None},

		ControlPlacement->Left,
		SaveDefinitions->True,
		TrackedSymbols:>{bw, ns, te, timeShift, eps, gamma, sigma, ph0, 
			ph0s, ph1s, off, pSize, snr, pad, ap, show, pRan, 
			pMin, pMax, met, simTrigger, vals, use
		}
	];

	NotebookClose[spectraSimulator];
	spectraSimulator = CreateWindow[DialogNotebook[{CancelButton["Close", DialogReturn[]], man},
		WindowSize -> All, WindowTitle -> "Spectra simulator"]];
];


simulatedFid[{nSamp_, bw_, field_, nuc_}, namesI_]:=Block[{dw, gyro, names, fids, specs, table},
	dw = 1. / bw;
	gyro = GetGyro[field, nuc];
	{names, fids, specs, table} = GetSpectraBasisFunctions[namesI,
		BasisSequence->{"PulseAcquire", 0}, SpectraSamples->nSamp, SpectraBandwidth->bw,
		SpectraPpmShift->0, SpectraFieldStrength->field, SpectraNucleus->nuc];
	{dw, gyro, field, Thread[names -> fids]}
]


apodizeFun[time_,apM_:"Hanning"]:=Block[{length,app,xdat,xmax},
	xdat = time;
	length = Length@time;
	xmax = Max[Abs[xdat]];
	app = Switch[apM,
		"Hanning", 0.5 + 0.5 Cos[xdat Pi/xmax],
		"Hamming", 0.54 + 0.46 Cos[xdat Pi/xmax],
		"Lorentzian", Exp[-(3./xmax) xdat],
		"Gaussian", Exp[-(2./xmax)^2 xdat^2],
		"Voigt", 0.5 Exp[-(3./xmax) xdat] + 0.5 Exp[-(2./xmax)^2 xdat^2]
	];
	app = app/Max[app];
	app
]


relaxFun[time_, {gamma_, sigma_}] := Exp[- Pi (gamma time + 2 Pi sigma^2 time^2)]
shiftFun[time_, eps_, gyro_] := Re@Exp[2 Pi eps gyro I time] 


effectiveLinewidth[gamma_, sigma_, gyro_] := Block[{fwhmHz},
	(* Olivero Voigt Approximation *)
	fwhmHz = 0.5346 * gamma + Sqrt[0.2166 gamma^2 + 2.3548 sigma^2];
	{fwhmHz, fwhmHz / gyro}
]


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
			shift = FindSpectraPpmShift[spec,{1./bw,gyro}, {{0, -2.52, -7.56, -16.15}, {2, 1, 1, 1}}]; 
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


(* ::Subsection:: *)
(*Spar and Sdat*)


(* ::Subsubsection::Closed:: *)
(*ImportSparSdat*)


SyntaxInformation[ImportSparSdat] = {"ArgumentsPattern" -> {_, _}};

ImportSparSdat[fSpar_,fSdat_]:=Block[{numType,head,depth,row,nSamp,x,y,z,bw,te,nuc,gyro,field,fids,specs},
	(*read raw data and header*)
	numType=FromVaxD[BinaryReadList[fSdat, "UnsignedInteger32"]];
	head=ParseSpar[Import[fSpar, "Lines"]];

	(*get header information*)
	{bw,te,nuc,gyro}={"sample_frequency","spectrum_echo_time","nucleus","synthesizer_frequency"}/.head;
	gyro = gyro 10^-6;
	field=Round[gyro /GyromagneticRatio[nuc],.1];

	{depth,row,nSamp,x,y,z}={"num_dimensions","spec_num_row","dim1_pnts","dim2_pnts","dim3_pnts","dim4_pnts"}/.head;
	(*partition the data to the correct size*)
	fids=Switch[depth,
		2,
		If[row==1,
			(*single voxel*)
			RotateDimensionsRight[Fold[Partition,numType,{2}]],
			(*dynamic series*)
			RotateDimensionsRight[Fold[Partition,numType,{2,nSamp}]]
		],
		3,
		RotateDimensionsRight[Fold[Partition,numType,{2,nSamp,x}]],
		4,(*3D CSI*)
		RotateDimensionsRight[Fold[Partition,numType,{2,nSamp,x,y}]]
	];
	fids=fids[[1]]+fids[[2]]I;

	(*convert fid to spec*)
	specs=If[ArrayDepth[fids]>1, Map[ShiftedFourier,fids,{-2}],	ShiftedFourier[fids]];

	(*output data*)
	{specs,{bw,te},{gyro,nuc},head}
]


(* ::Subsubsection::Closed:: *)
(*ParseSpar*)


ParseSpar[sparI_]:=Block[{textIn,values,lab,val,spar},
	spar=Select[sparI,#=!=""&];
	textIn=Select[spar,StringTake[#,1]==="!"&];
	values=StringSplit[#," :"]&/@Select[spar,StringTake[#,1]=!="!"&];
	{lab,val}=Transpose[If[Length[#]===1,Join[#,{""}],StringTrim/@#]&/@values];
	val=If[StringMatchQ[#,NumberString],ToExpression[#],#]&/@val;
	Thread[lab->val]
]


(* ::Subsubsection::Closed:: *)
(*FromVaxD*)


FromVaxD=Compile[{{int,_Integer,0}},Block[{bin,sign ,fraction, exponent},
	(*pad to 32 binary numbers*)
	bin=IntegerDigits[int,2,32];
	(*sign is bit 17*)
	sign=(-1)^bin[[17]];
	(*fraction is stored at the end the beginning of the binary number*)
	fraction=(FromDigits[Join[bin[[26;;-1]],bin[[1;;16]]],2]/16777216.)+0.5;
	(*the exponent is the center*)
	exponent=2.^(FromDigits[bin[[18;;25]],2]-128.);
	(*create the number*)
	sign fraction exponent
],RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


(* ::Subsubsection::Closed:: *)
(*ExportSparSdat*)


Options[ExportSparSdat]={SparName->"QMRITools", SparOrientation->{0,0},SparID->""}

SyntaxInformation[ExportSparSdat] = {"ArgumentsPattern" -> {_, _, {_,_}, {_,_}, _., OptionsPattern[]}};

ExportSparSdat[file_, specs_, {bw_, te_}, {gyro_, nuc_}, opts:OptionsPattern[]] := ExportSparSdat[file, specs, {bw, te}, {gyro, nuc}, {1,1,1}, opts]

ExportSparSdat[file_, specs_, {bw_, te_}, {gyro_, nuc_}, vox_, opts:OptionsPattern[]]:=Block[{fidsOut, datOut, fileOut, headOut},
	(*export data*)
	fidsOut=Map[ShiftedInverseFourier,specs,{-2}];
	datOut=ToVaxD[Flatten[RotateDimensionsLeft[{Re@fidsOut,Im@fidsOut}]]];
	fileOut=file<>".SDAT";
	If[FileExistsQ[fileOut],DeleteFile[fileOut]];
	BinaryWrite[fileOut,datOut,"UnsignedInteger32"];
	Close[fileOut];

	(*export header*)
	headOut=MakeSpar[specs, {bw, te}, {gyro, nuc}, vox, opts];
	Export[file<>".SPAR", headOut, "text"];
]


(* ::Subsubsection::Closed:: *)
(*ToVaxD*)


ToVaxD=Compile[{{num,_Real,0}},Block[{signBin,numA,exp,expBin,frac,fracBin},
	If[num===0.,
		FromDigits[ConstantArray[0,32],2],
		(*get the sign bit ans remove sign*)
		signBin={UnitStep[Sign[-num]]};
		numA=Abs[num];
		(*get the exponent*)
		exp=Ceiling[1./Log[numA,2.]];
		expBin=PadLeft[IntegerDigits[exp+128,2],8];
		(*get the fraction*)
		fracBin=PadLeft[IntegerDigits[Round[16777216.(numA/(2.^exp)-0.5)],2],23];
		(*output the integer*)
		FromDigits[Join[fracBin[[8;;-1]],signBin,expBin,fracBin[[1;;7]]],2]
	]
],RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


(* ::Subsubsection::Closed:: *)
(*MakeSpar*)


Options[MakeSpar] = Options[ExportSparSdat];

MakeSpar[specs_, {bw_, te_}, {gyro_, nuc_}, vox_, OptionsPattern[]]:=Block[{
		dimZO,dimYO,dimXO, nSampO, gyroO, nucO, bwO, teO, nameO, hf, ps, id,
		text,lab,filHeader,fixedHeader, vals, head,row,depth
	},
	(*switch between data dimensions*)
	depth=ArrayDepth[specs];
	Switch[depth,
		1,(*single voxel*)
		nSampO=Length[specs];
		row=dimZO=dimYO=dimXO=1;
		,
		2,(*dynamic data*)
		{row,nSampO}=Dimensions[specs];
		dimZO=dimYO=dimXO=1;
		,
		3,(*2D CSI*)
		{dimYO,dimXO,nSampO}=Dimensions[specs];
		row=dimXO dimYO;
		dimZO=1;
		,
		4,(*3D CSI*)
		{dimZO,dimYO,dimXO,nSampO}=Dimensions[specs];
		row=dimXO dimYO;
	];

	(*mandatory input parameters*)
	{gyroO,nucO} = {10^6 gyro,nuc};
	{bwO, teO} = {bw, te};

	(*optional input parameters*)
	nameO = OptionValue[SparName];
	{hf,ps} = OptionValue[SparOrientation];(*head or feat first / prone or supine*)
	id = OptionValue[SparID];(*prefer the processing steps*)

	(*all needed fixed header information orders of text and lab are important*)
	text={
		"!--------------------------------------------------------------------","!","!","!      CAUTION - Investigational device.","!      Limited by Federal Law to investigational use.","!","!",
		"!      GYROSCAN spectro parameter file. ","!      Last revised 05-July-2007.","!--------------------------------------------------------------------","!   This file contains time domain data in the spectral dimension.",
		"!   S15/ACS: set of *.SPAR and *.SDAT files is created, (data format: VAX CPX floats)","!--------------------------------------------------------------------","!---------------------------------------------",
		"!--------------------","! Column parameters ","!--------------------","!-----------------","! Row parameters ","!-----------------","!-------------------------------------------------",
		"! Extra parameters in order to make data transfer ","! possible between S15/ACS and SUNSPEC: ","!-------------------------------------------------","!-------------------------------------------------",
		"!-----------------------------------------------------","!-----------------------------------------------------","!-----------------------------------------------------",
		"!---------------------------------------------------","!   Additional parameters","!---------------------------------------------------","!-----------------------------------------------------"
	};

	lab={
		"examination_name","scan_id","scan_date","patient_name","patient_birth_date","patient_position","patient_orientation","samples","rows","synthesizer_frequency","offset_frequency",
		"sample_frequency","echo_nr","mix_number","nucleus","t0_mu1_direction","echo_time","repetition_time","averages","volume_selection_enable","t1_measurement_enable","t2_measurement_enable",
		"time_series_enable","phase_encoding_enable","nr_phase_encoding_profiles","si_ap_off_center","si_lr_off_center","si_cc_off_center","si_ap_off_angulation","si_lr_off_angulation",
		"si_cc_off_angulation","t0_kx_direction","t0_ky_direction","nr_of_phase_encoding_profiles_ky","phase_encoding_direction","phase_encoding_fov","slice_thickness",
		"image_plane_slice_thickness","slice_distance","nr_of_slices_for_multislice","Spec.image in plane transf","spec_data_type","spec_sample_extension","spec_num_col",
		"spec_col_lower_val","spec_col_upper_val","spec_col_extension","spec_num_row","spec_row_lower_val","spec_row_upper_val","spec_row_extension","num_dimensions","dim1_ext",
		"dim1_pnts","dim1_low_val","dim1_step","dim1_direction","dim1_t0_point","dim2_ext","dim2_pnts","dim2_low_val","dim2_step","dim2_direction","dim2_t0_point","dim3_ext",
		"dim3_pnts","dim3_low_val","dim3_step","dim3_direction","dim3_t0_point","dim4_ext","dim4_pnts","dim4_low_val","dim4_step","dim4_direction","dim4_t0_point","echo_acquisition",
		"TSI_factor","spectrum_echo_time","spectrum_inversion_time","image_chemical_shift","resp_motion_comp_technique","de_coupling","equipment_sw_versions","placeholder1","placeholder2"
	};

	(*from input*)
	filHeader={
		(*general acquisition names*)
		"patient_position"->Switch[hf, 0, "\"head_first\"", 1, "\"feat_first\""],
		"patient_orientation" ->Switch[ps, 0, "\"supine\"", 1, "\"prone\""],
		(*get from input window*)
		"examination_name"->"Generated by QMRITools",
		"patient_name"->nameO,
		"phase_encoding_fov"->vox[[2]] dimYO,(*mm fov in freq*)
		"slice_thickness"->vox[[1]] dimZO,(*mm fov in phase*)
		"slice_distance"->If[depth>2, vox[[1]], 0],(*mm slice thickness*)
		"repetition_time"->0,(*ms could be an input parameter but not relevant for now*)
		(*save date*)
		"scan_date"->StringRiffle[ToString/@DateList[Today][[1;;3]],"."]<>" "<>StringRiffle[ToString/@Round[DateList[Now][[-3;;]]],":"],
		(*get from gui with processing settings*)
		"scan_id"->id (*a string describing the data*),

		(*get from processing tool*)
		"synthesizer_frequency"->gyroO,(*MHz*)
		"sample_frequency"->bwO,(*Hz*)
		"nucleus"->nucO,(*string*)
		"echo_time"->teO,(*ms*)
		"spectrum_echo_time"->teO,

		(*get from data dimensions*)
		"num_dimensions"->Clip[depth,{2,4}],
		(*dynamics if needed*)
		"rows"->row,
		"spec_row_upper_val"->row,
		"spec_num_row"->row,
		(*fid parameters and time*)
		"samples" ->nSampO,
		"spec_num_col"->nSampO,
		"dim1_pnts"->nSampO,
		"dim1_step"-> 1./bwO,(*s*)
		"spec_col_upper_val"->(nSampO-1)(1./bw),(*s*)
		(*data dimensions*)
		"dim2_pnts"->dimXO,
		"dim3_pnts"->dimYO,
		"nr_of_slices_for_multislice"->dimZO,
		"nr_phase_encoding_profiles"->dimXO,
		"nr_of_phase_encoding_profiles_ky"->dimYO,
		"phase_encoding_enable"-> If[depth<=2,"\"no\"","\"yes\""],
		"dim4_pnts" ->dimZO
	};

	(*fixed parameters that are default*)
	fixedHeader={
		(*possible non fixed parameters*)
		"offset_frequency" ->0,"averages"->1,"echo_nr"-> 1,"mix_number"-> 1,
		"spectrum_inversion_time"->0,"image_chemical_shift"->0,
		(*start fixed parameters*)
		"patient_birth_date"->"1980.11.25","image_plane_slice_thickness"->0,
		(*general fixed parameters*)
		"t0_mu1_direction"->0,
		"volume_selection_enable"->"\"no\"","t1_measurement_enable"->"\"no\"","t2_measurement_enable"->"\"no\"","time_series_enable"->"\"no\"",

		
		"si_ap_off_center"->0,"si_lr_off_center"->0,"si_cc_off_center"->0,"si_ap_off_angulation"->0,
		"si_lr_off_angulation"->0,"si_cc_off_angulation"->0,
		"t0_kx_direction"->50,"t0_ky_direction"->50,
		"phase_encoding_direction"->"\"trans\"",

		
		"Spec.image in plane transf"->"\"minA-minB\"",
		(*column and row fixed parameters*)
		"spec_data_type"->"cf","spec_sample_extension"->"[V]",
		"spec_col_lower_val"->0,"spec_col_extension"->"[sec]",
		"spec_row_lower_val"->1,"spec_row_extension"->"[index]",
		(*per dim fixed parameters*)
		"dim1_ext"->"[sec]","dim1_low_val"->0,"dim1_direction"->"mu1","dim1_t0_point"->0,
		"dim2_ext"->"[num]","dim2_low_val"->1.0,"dim2_step"->1.0,"dim2_direction"->"x","dim2_t0_point"->50,
		"dim3_ext"->"[num]","dim3_low_val"->1.0,"dim3_step"->1.0,"dim3_direction"->"y","dim3_t0_point"->50,
		"dim4_ext"->"[index]","dim4_low_val"->1.0,"dim4_step"->1.0,"dim4_direction"->"slice","dim4_t0_point"->"-",
		(*closing values*)
		"echo_acquisition"->"NO","TSI_factor"->0,"resp_motion_comp_technique"->"NONE","de_coupling"->"NO",
		"equipment_sw_versions"->"QMRITools","placeholder1"->"","placeholder2"->""
	};

	(*generate the header values*)
	vals=lab/.Join[filHeader,fixedHeader];
	head=Thread[lab->vals];
	head=StringReplace[#[[1]]<>" : "<>If[StringQ[#[[2]]],#[[2]],If[IntegerQ[#[[2]]],ToString[#[[2]]],ToString[DecimalForm[#[[2]],{20,15}]]]],": ["->":["]&/@head;

	StringRiffle[Join[text[[;;13]],head[[;;41]],text[[{14}]],head[[42;;43]],
		text[[15;;17]],head[[44;;47]],
		text[[18;;20]],head[[48;;51]],
		text[[21;;24]],head[[{52}]],
		text[[{25}]],head[[53;;58]],
		text[[{26}]],head[[59;;64]],
		text[[{27}]],head[[65;;70]],
		text[[{28}]],head[[71;;76]],
		text[[29;;31]],head[[77;;86]],
		text[[{32}]],{""}
	],"\n\n"]
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
