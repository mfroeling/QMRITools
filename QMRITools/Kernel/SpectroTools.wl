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
Output is the {time, spec, {begintime, samplingInterval}}."

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
"PadSpectra[spec] doubles the number of spectral points while maintainig the dwell time."

ApodizeSpectra::usage = 
"ApodizeSpectra[spec] performs apodization of the spectra. The apodization function is set with the option ApodizationFunction."

ApodizePadSpectra::usage = 
"ApodizePadSpectra[spec] and doubles the number of spectral points while maintainig the dwell time."


GetTimePpmRange::usage = 
"GetTimePpmRange[spec, {dt, field, nuc}] get the timing of the fid and the ppm values of the spec where dt is the well time in ms, field the field strength in Tesla and nuc the nucleus availible in GyromagneticRatio. 
GetTimePpmRange[spec, dt, field, nuc] get the timing of the fid and the ppm values of the spec. 
GetTimePpmRange[spec, dt, gyro] get the timing of the fid and the ppm values of the spec."

GetTimeRange::usage = 
"GetTimeRange[fid, dt] get the timing of the fid where dt is the well time in ms."

GetPpmRange::usage = 
"GetPpmRange[spec, {dt, field, nuc}] get the ppm values of the spec where dt is the well time in ms, field the field strength in Tesla and nuc the nucleus availible in GyromagneticRatio. 
GetPpmRange[spec, dt, field, nuc] get the ppm values of the spec. 
GetPpmRange[spec, dt, gyro] get the ppm values of the spec."

GetGyro::usage = 
"GetGyro[nuc, field] geth the gyromagnetic ratio with field the field strength in Tesla and nuc the nucleus availible in GyromagneticRatio."


PhaseShiftSpectra::usage = 
"PhaseShiftSpectra[spectra, phi0] aplies the 0th order phase phi0 to the spectra. 
PhaseShiftSpectra[spectra, ppm, gyro, phi1] aplies the 1st order phase phi1 to the spectra. The ppm can be obtained using GetPpmRange and gyro with GetGyro. 
PhaseShiftSpectra[spec, ppm, gyro, {phi0, phi1}] aplies the 0th and 1st order phases {phi0, phi1} to the spectra. The ppm can be obtained using GetPpmRange and gyro with GetGyro.

The 0th order phase phi0 is in radians and the 1st order phase phi1 is in ms."

PhaseShiftFid::usage = 
"PhaseShiftFid[spectra, phi0] aplies the 0th order phase phi0 to the FID."

TimeShiftFidV::usage = 
"TimeShiftFidV[fid, time, gam] aplies a line broadening with linewidth gam and a Voigt line shape to the fid. The time can be obtained using GetTimeRange.
TimeShiftFidV[fid, time, {gamL, gamG}] aplies a line broadening with linewidth gamG \"Gaussian\" and gamL \"Lorentzian\".
TimeShiftFidV[fid, time, gyro, {gam, eps}] aplies a line broadening with linewidth gam to the fid and a phase eps that results in eps ppm shift of the spectra. The gyro can be obtained with GetGyro.
TimeShiftFidV[fid, time, gyro, {{gamL, gamG}, eps}] aplies a line broadening with linewidth linewidth gamG \"Gaussian\" and gamL \"Lorentzian\" to the fid and a phase eps that results in eps ppm shift of the spectra.

The linewidth gam is given in ms and the spectra shift eps is given in ppm."

TimeShiftFid::usage = 
"TimeShiftFid[fid, time, gam] aplies a line broadening with linewidth gam and a Voigt line shape to the fid. The time can be obtained using GetTimeRange.
TimeShiftFid[fid, time, {gam, f}] aplies a line broadening with linewidth gam and a custom line shape f to the fid (f=0, \"Gaussian\", f=1 \"Lorentzian\").
TimeShiftFid[fid, time, gyro, {gam, eps}] aplies a line broadening with linewidth gam to the fid and a phase eps that results in eps ppm shift of the spectra. The gyro can be obtained with GetGyro.
TimeShiftFid[fid, time, gyro, {gam, eps, f}] aplies a line broadening with linewidth gam using a custom line shape f to the fid and a phase eps that results in eps ppm shift of the spectra.

The linewidth gam is given in ms and the spectra shift eps is given in ppm."

TimeShiftEcho::usage =
"TimeShiftEcho[fid, time, gam] aplies a line broadening with linewidth gam and a Voigt line shape to the fid. The time can be obtained using GetTimeRange.
TTimeShiftEcho[fid, time, {gam, f}] aplies a line broadening with linewidth gam and a custom line shape f to the fid (f=0, \"Gaussian\", f=1 \"Lorentzian\").
TTimeShiftEcho[fid, time, gyro, {gam, eps}] aplies a line broadening with linewidth gam to the fid and a phase eps that results in eps ppm shift of the spectra. The gyro can be obtained with GetGyro.
TTimeShiftEcho[fid, time, gyro, {gam, eps, f}] aplies a line broadening with linewidth gam using a custom line shape f to the fid and a phase eps that results in eps ppm shift of the spectra.

The linewidth gam is given in ms and the spectra shift eps is given in ppm."

ShiftSpectra::usage = 
"ShiftSpectra[spectra, {dw, gyro}, shift] shifts the spectra by shift. The shift is in ppm." 


ChangeDwellTimeFid::usage = 
"ChangeDwellTimeFid[fid, dt, dtnew] changes the sampling time of an fid from dwelltime dt to dwelltime dtnew."

PhaseCorrectSpectra::usage =
"PhaseCorrectSpectra[spec] performs 0th order phase correction of the spectra by minimizing the difference between the real and absolute spectra velaue.
PhaseCorrectSpectra[spec, dw] performs 0th order phase correction of the spectra using Henkel matrix SVD fitting.
PhaseCorrectSpectra[spec, dw, te] := performs 0th and 1st order phase correction of the spectra using Henkel matrix SVD fitting. The first order phase is corrected by padding the fid with the missing values in the time befroe the TE.
PhaseCorrectSpectra[spec, dw, te, gyro, ppmRan] performs 0th and 1st order phase correction of the spectra using Henkel matrix SVD fitting. Only the part of the spectra in the ppmRan is used for optimization."

CorrectTESpec::usage = 
"CorrectTESpec[spectra, dw, te] corrects the spectra for 1st order phase by extrapolating the missing FID samples in the TE using Henkel matrix SVD analsis.
CorrectTESpec[spectra, dw, te, gyro, ppmRan] corrects the spectra for 1st order phase by extrapolating the missing FID samples in the TE using Henkel matrix SVD analsis. Only the part of the spectra in the ppmRan is used for optimization."

CorrectTEFid::usage = 
"CorrectTEFid[fid, dw, te] corrects the fid for 1st order phase by extrapolating the missing FID samples in the TE using Henkel matrix SVD analsis.
CorrectTEFid[fid, dw, te, gyro, ppmRan] corrects the fid for 1st order phase by extrapolating the missing FID samples in the TE using Henkel matrix SVD analsis. Only the part of the spectra in the ppmRan is used for optimization."


FitSpectra::usage = 
"FitSpectra[specBasis, spec, {st,end}, dt, {lwvals,lwamsp}] Fits the basis spectra from GetSpectraBasisFunctions to the spec overt the ppm range {st, end} and dt the dweltime."

GetSpectraBasisFunctions::usage = 
"GetSpectraBasisFunctions[{met1, ..., metn}] generates a list of spectra baisis functions with names met1 to metn. The names are strings and are the metabolites availible in GetSpinSystem.
GetSpectraBasisFunctions[{{props1}, ..., {propsn}}] generates a list of spectra baisis functions with properties prop1 to propn. The properties are those specified in MakeSpinSystem.
GetSpectraBasisFunctions[inp, split] generates a list of spectra basisfunctions. Each metabolite name present in the list split wil be split in individual spectra per peak."


PlotCSIData::usage =
"PlotCSIData[spectra, {dwell, gyro}] plots the CSI spectra which has dimensions {z,y,x,nsamp}. The ppm axes is determined by dwell and gyro. Gyro can be obtained with GetGyro.
PlotCSIData[spectra, {dwell, field, nuc}] plots the CSI spectra which has dimensions {z,y,x,nsamp}. The ppm axes is determined by dwell and field and nuc."

PlotFid::usage = 
"PlotFid[fid, dwell] plots the fid assuming dwell as the sampling time.
PlotFid[time, fid] plot the fid where time is the timing of the fid which can be obtained with GetTimeRange."

PlotSpectra::usage = 
"PlotSpectra[spectra, {dwell, gyro}] plots the spectra, the ppm axes is determined by dwell and gyro. Gyro can be obtained with GetGyro.
PlotSpectra[spectra {dwell, field, nuc}] plots the spectra, the ppm axes is determined by dwell field and nuc.
PlotSpectra[ppm, spectra] plots the spectra where ppm is the pmm range of the spectra which can be obtained with GetPpmRange." 


FindSpectraPpmShift::usage = 
"FindSpectraPpmShift[spectra, {dw, gyro}, peaks] finds the ppm value that aligns the spectra with the given peak positions peaks wich is a list of ppm values. 
FindSpectraPpmShift[spectra, {dw, gyro}, {peaks, amps}] finds the ppm value that aligns the spectra with the given peak positions peaks wich is a list of ppm values and amps are ther relative amplitudes. 
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
"ImportSparSdat[fspar, fsdat] imports sfpar and fsdata file. Function not done."

ExportSparSdat::usage =
"ExportSparSdat[file, specs, {bw ,te}, {gyro ,nuc}] exports specs to file. Function not done."


(* ::Subsection::Closed:: *)
(*Options*)


ApodizationFunction::usage = 
"ApodizationFunction is an options for ApodizeFid, ApodizeSpectra, ApodizePadFid, and ApodizePadSpectra. Values can be \"Hanning\", \"Hamming\", \"Gaussian\", \"Lorentzian\", and \"Voigt\"."

PaddingFactor::usage = 
"PaddingFactor is an option for PadFid, PadSpectra, ApodizePadFid, ApodizePadSpectra and FitSpectra. It Specifies with which factro to lengthen the fid."

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
"SpectraFieldStrength is an option for GetSpectraBasisFunctions and FitSpectra and sets the field strenght at which the simulations and fitting is perforemd."

SplineSpacingFactor::usage = 
"SplineSpacingFactor is an option for FitSpectra and defines the distance between the bsplien points relative the the mean linewithd of the peaks."

FineTuneFit::usage = 
"FineTuneFit is an option for FitSpectra and when True it performs a second fitting run where for each peak is an individual linewidth, line shape and shift are fitted."

InitializeFit::usage = 
"InitializeFit is an option for FitSpectra and is used to set initila values for the global fit {gami,epsi,{phi0i,phi1i},line shape}."

FitLineShape::usage = 
"FitLineShape is an option for FitSpectra and when True allows to fit the lineshap. If False a voigt line shape is used."

SpectraOutputPlots::usage = 
"SpectraOutputPlots is an option for FitSpectra. If True the automatica calibration plot for the initial fit are generated."

ReadoutType::usage = 
"ReadoutType is an option for FitSpectra and padding and apodization fucntions. Value can be \"Fid\" or \"Echo\"."

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


PhaseCorrectError[speci_, phi0_?NumericQ] := PhaseCorrectErrorC[speci, phi0]

PhaseCorrectErrorC = Compile[{{speci, _Complex, 1}, {phi0, _Real, 0}},Block[{specR},
	specR= Re[Exp[-I phi0] speci];
	Total[(Abs[speci] - specR)^2]-Total[specR]
],RuntimeOptions -> "Speed", Parallelization -> True];


PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, out_:True] := PhaseCorrectSpectra[spec, dw, 0, 0, Full, out]

PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, te_?NumberQ, out_:True] := PhaseCorrectSpectra[spec, dw, te, 0, Full, out]

PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, gyro_?NumberQ, ppmRan_, out_:True] := PhaseCorrectSpectra[spec, dw, 0, gyro, ppmRan, out]

PhaseCorrectSpectra[spec_?ListQ, dw_?NumberQ, te_?NumberQ, gyro_?NumberQ, ppmRan_, out_:True] := Module[{fid, specOut, missing, full, henkelSpec, phi, phi0},
	(*create the fid*)
	fid = ShiftedInverseFourier[spec];

	(*perform the HenkelFit*)
	{missing, full} = HenkelFit[fid ,dw, te, gyro, ppmRan];

	(*create the full fid*)
	specOut = ShiftedFourier[Join[missing, fid][[;;Length[fid]]]];

	(*create the Henkle spectra*)
	henkelSpec = ShiftedFourier[full];
	If[ppmRan =!= Full, henkelSpec = Pick[henkelSpec, Unitize[Clip[GetPpmRange[henkelSpec, dw, gyro], Sort[ppmRan], {0, 0}]], 1]];

	(*find the first order phase*)
	Clear[phi0];
	phi = Quiet[Last[FindMaximum[PhaseErrorH[henkelSpec, phi0], phi0]][[1,2]]];
	phi = If[NumberQ[phi], phi, phi[[2]]];
	If[out,
		specOut Exp[-I phi],
		phi
	]
]


PhaseErrorH[speci_, phi0_?NumericQ] := PhaseErrorHC[speci, phi0]

PhaseErrorHC = Compile[{{speci, _Complex, 1}, {phi0, _Real, 0}}, 
	Total[Re[Exp[-I phi0] speci]], 
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
	missing = HenkelFit[fid ,dw, te, gyro, ppmRan][[1]];
	(*output full fid*)
	Join[missing, fid][[;;Length[fid]]]
]


(* ::Subsubsection::Closed:: *)
(*HenkelFit*)


HenkelFit[fid_ ,dw_, te_, gyro_, ppmRan_]:=Block[{timeOr, timeMis, henk, fit},
	(*get the correct timing of the fid and missing values*)
	timeOr = dw Range[1, Length[fid]];
	timeMis = Reverse@Range[te - dw, If[Mod[te, dw] > 0.5 dw, -0.5 dw, 0], -dw];

	(*get the henkel values*)
	henk = HenkelSVDFid[fid, dw, gyro, ppmRan];
	fit = (PseudoInverse[HenkelSVDBasisC[timeOr, henk]] . fid);

	(*missing and full henkle fid*)
	{If[timeMis =!= {}, HenkelSVDBasisC[timeMis,henk] . fit,{}],	HenkelSVDBasisC[timeOr, henk] . fit}
]


(* ::Subsubsection::Closed:: *)
(*HenkelSVDFid*)


HenkelSVDFid[fid_, dw_] := HenkelSVDFid[fid, dw, 0, Full]

HenkelSVDFid[fid_, dw_, gyro_] := HenkelSVDFid[fid, dw, gyro, Full]

HenkelSVDFid[fid_, dw_, gyro_, ppmRan_] := Block[{
	lmax, mmax, H, U, Utr, Uup, Udown, q, freq, decay, ppmval, sel, ppmMax, ppmMin
	},

	(*get the matrix size*)
	lmax = Round[0.4 Length[fid]];
	mmax = Length[fid] + 1 - lmax;

	(*create the henkel matrix and the singula values*)
	H = fid[[Range[mmax] + #]] & /@ Range[0, lmax - 1];
	U = First@SingularValueDecomposition[H, 32];
	q = Log[Eigenvalues[PseudoInverse[U[[;; -2]]] . U[[2 ;;]]]];

	(*get the frequencies and delay times*)
	decay = Re[q]/dw;
	freq = Im[q]/(2 Pi dw);

	(*select the the correct frequency range and only long decay times*)
	sel = If[ppmRan === Full, 
		UnitStep[Abs[1000/decay] - Pi], 
		ppmval = -(freq/gyro);
		{ppmMin, ppmMax} = Sort[ppmRan];
		(1 - UnitStep[ppmval + ppmMin]) UnitStep[ppmval + ppmMax] UnitStep[Abs[1000/decay] - Pi]
	];
	sel = If[Total[sel] === 0, sel + 1, sel];

	(*give the henkel output*)
	{Pick[decay, sel, 1], Pick[freq, sel, 1]}
]


HenkelSVDBasisC = Compile[{{time, _Real, 1}, {henkel, _Real, 2}}, 
	Transpose[Exp[# time] & /@ (henkel[[1]] + 2 Pi I henkel[[2]])], 
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
	func = Switch[readout,"Fid",PadFid,"Echo",PadEcho];
	ShiftedFourier[func[ShiftedInverseFourier[spec, readout],PaddingFactor-> OptionValue[PaddingFactor]], readout]
]


(* ::Subsubsection::Closed:: *)
(*ApodizeSpectra*)


Options[ApodizeSpectra] = {ApodizationFunction -> "Hanning", ReadoutType -> "Fid"}

SyntaxInformation[ApodizeSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeSpectra[spec_, OptionsPattern[]] :=  Block[{func,readout},
	readout = OptionValue[ReadoutType];
	func = Switch[readout,"Fid",ApodizeFid,"Echo",ApodizeEcho];
	ShiftedFourier[func[ShiftedInverseFourier[spec, readout], ApodizationFunction-> OptionValue[ApodizationFunction]], readout]
]


(* ::Subsubsection::Closed:: *)
(*ApodizePadSpectra*)


Options[ApodizePadSpectra] = {ApodizationFunction -> "Hanning", PaddingFactor -> 2, ReadoutType -> "Fid"}

SyntaxInformation[ApodizePadSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadSpectra[spec_, OptionsPattern[]] := Block[{func,readout},
	readout = OptionValue[ReadoutType];
	func = Switch[readout,"Fid",ApodizePadFid,"Echo",ApodizePadEcho];
	ShiftedFourier[func[ShiftedInverseFourier[spec, readout], ApodizationFunction->OptionValue[ApodizationFunction], PaddingFactor->OptionValue[PaddingFactor]],readout]
]


(* ::Subsection:: *)
(*Time and ppm*)


(* ::Subsubsection::Closed:: *)
(*GetTimePpmRange*)


SyntaxInformation[GetTimePpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetTimePpmRange[spec_, {dt_, field_, nuc_}] := GetTimePpmRange[spec, dt, GetGyro[nuc, field]]

GetTimePpmRange[spec_, dt_, field_, nuc_] := GetTimePpmRange[spec, dt, GetGyro[nuc, field]]

GetTimePpmRange[spec_, dt_, gyro_] := {GetTimeRange[spec, dt], GetPpmRange[spec, dt, gyro]}


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

GetTimeRange[fid_?VectorQ, dt_] := GetTimeRange[Length[fid], dt]

GetTimeRange[len_?IntegerQ, dt_] := dt (Range[len]-1)


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
	Exp[-I (phi0 + 2 Pi (phi1/1000) gyro ppm)] spec,
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

TimeShiftFid[fid_, time_, gam_] := TimeShiftFidC[fid, time, 0., gam, 0., .5];

TimeShiftFid[fid_, time_, {gam_, f_}] := TimeShiftFidC[fid, time, 0., gam, 0., f];

TimeShiftFid[fid_, time_, gyro_, {gam_, eps_}] := TimeShiftFidC[fid, time, gyro, gam, eps, .5];

TimeShiftFid[fid_, time_, gyro_, {gam_, eps_, f_}] := TimeShiftFidC[fid, time, gyro, gam, eps, f];

TimeShiftFidC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {gam, _Real, 0}, {eps, _Real, 0}, {f, _Real, 0}},
	Exp[-(f gam time + (1 - f) gam^2 time^2) + 2 Pi eps gyro I time] fid, 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


SyntaxInformation[TimeShiftFidV] = {"ArgumentsPattern" -> {_, _, _., _.}}

TimeShiftFidV[fid_, time_, gam_] := TimeShiftFidC2[fid, time, 0, gam, gam, 0.];

TimeShiftFidV[fid_, time_, {gamL_, gamG_}] := TimeShiftFidC2[fid, time, 0., gamL, gamG, 0.];

TimeShiftFidV[fid_, time_, gyro_, {gam_, eps_}] := TimeShiftFidC2[fid, time, gyro, gam, gam, eps];

TimeShiftFidV[fid_, time_, gyro_, {{gamL_, gamG_}, eps_}] := TimeShiftFidC2[fid, time, gyro, gamL, gamG, eps];

TimeShiftFidC2 = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {gamL, _Real, 0}, {gamG, _Real, 0}, {eps, _Real, 0}},
	Exp[-(gamL time + (gamG time)^2)] Exp[2 Pi eps gyro I time] fid, 
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
	(f Exp[-gam Abs[timeNew]] + (1 - f) Exp[- gam^2 Abs[timeNew]^2]) Exp[2 Pi eps gyro I timeNew] fid
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

FindSpectraPpmShift[spec_, {dw_, gyro_}, {peaks_, amp_}] := Block[{ppm, dppm, tar, amps, wgth, s, corrf, maxPos, sft},
	ppm = GetPpmRange[spec, dw, gyro];
	dppm = (ppm[[1]] - ppm[[2]]);

	tar = If[Length[spec] === Length[peaks],
		Abs[peaks],
		amps = If[amp === 0, 0. peaks + 1, amp];
		tar = 0 ppm;
		tar[[Flatten[Position[ppm, First@Nearest[ppm, #]] & /@ peaks]]] = amps/Max[amps];
		tar
	];

	(*perfomr correlation of spectra with delta function*)
	s = Max[ppm]/4;
	wgth = 1/(Exp[(ppm^2/(2*s^2))]*(Sqrt[2*Pi]*s));
	wgth = 1;
	corrf = ListCorrelate[tar, Abs[spec], Round[Length[tar]/2], 0];
	corrf = wgth DivideNoZero[corrf, Max[corrf]];
	maxPos = Position[corrf, Max[corrf]][[1, 1]];

	(*contrain shift to 5ppm*)
	sft = -dppm (maxPos - (Length[ppm]/2.));
	If[-5 < sft < 5, sft, 0]
]


(* ::Subsubsection::Closed:: *)
(*ChangeDwellTimeFid*)


SyntaxInformation[ChangeDwellTimeFid] = {"ArgumentsPattern" -> {_, _, _}}

ChangeDwellTimeFid[time_, dwOrig_, dwTar_] := Block[{NsampOrig, timeOrig, NsampTar, timeTar, sc},
	sc = dwOrig/dwTar;
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
	cf, bw, nsamp, lw, lws, read, sys, din, struct, dr, te, svals, seq, readout, t1, t2, necho,
	sysAll, spinTab, field, nuc, nam, labs, names, times, fids, ppms, specs, douts},

	(*get the option values*)
	{seq, svals} = OptionValue[BasisSequence];
	readout = Switch[seq, "PulseAcquire"|"STEAM","Fid", "SpaceEcho","Echo"];

	nsamp = OptionValue[SpectraSamples];
	bw = OptionValue[SpectraBandwidth];
	cf = OptionValue[SpectraPpmShift];
	field = OptionValue[SpectraFieldStrength];
	nuc = OptionValue[SpectraNucleus];

	(*get or make the spin systems*)
	sysAll = If[StringQ[#], GetSpinSystem, MakeSpinSystem][#, CenterFrequency -> cf] & /@ inp;
	sysAll = Reverse@SortBy[sysAll,First[#[[2]]] &];
	spinTab = SysTable[sysAll];

	(*loop over systems to perform sysmulation of basis functions*)
	{names, fids, specs} = Transpose[Flatten[
		Map[(
			(*get the system*)
			sys = #;
			nam = sys[[7]];
			labs = sys[[5]];

			(*see if molicule has to be splitted into individual peaks*)
			read = If[MemberQ[split, sys[[-1]]], "each", "all"];

			(*generate hamiltonian*)
			{din, struct} = SimHamiltonian[sys, FieldStrength -> field, SimNucleus -> nuc];

			(*perform sequence*)
			Switch[seq,
				"PulseAcquire",
				te = svals;
				dr = SequencePulseAcquire[din, struct, te],
				"STEAM",
				{te, tm} = svals;
				dr = SequenceSteam[din, struct, {te, tm}],
				"SpaceEcho",
				{t1, t2, necho} = svals;
				dr = SequenceSpaceEcho[din, struct, t1, t2, necho, 1]
			];

			(*simulate readout*)
			fids = 4 First@SimReadout[dr, struct, ReadoutSamples -> nsamp, ReadoutBandwidth -> bw, CenterFrequency -> cf, Linewidth -> 0, ReadoutOutput -> read];

			(*check if multiple fids for metabolite and act acordingly*)
			If[VectorQ[fids],
				specs = ShiftedFourier[fids, readout];
				{{nam, fids, specs}}
				,
				(*mulitple fids*)
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

FitSpectra[specBasisIn_, specIn_, {st_,end_}, dtime_, lwvals_?VectorQ, opts : OptionsPattern[]]:=FitSpectra[specBasisIn,specIn,{st,end},dtime,{lwvals, 0lwvals + 1.},opts]

FitSpectra[specBasisIn_, specIn_, {st_,end_}, dtime_, {lwvals_?VectorQ, lwamsp_?VectorQ}, OptionsPattern[]]:=Block[{
	ttotal,log,pad,spfac,field,nuc,shift,plots,init,scale,nbas,len,
	timeBasis,specFull,timeFull,ppmFull,nsamp,gyro,indSt,indEnd,
	gami,epsi,phi0i,phi1i,linei,phii,plLine,plShift, readout, varf,
	splineSpace,cpn,var,phi0f,phi1f,gamf,epsf,phif,linef,gam,eps,line,phi,sigi,
	tfit1,fit1,sol,output,tfit2,fit2,fit,timeBasisIn,time,ppm,spline,basis,error,errors,specFit
	},

	(*time the fitting*)
	ttotal=AbsoluteTiming[
		(*turn of error messages of FinMinimum*)
		Off[FindMinimum::cvmit];Off[FindMinimum::lstol];

		(*logging*)
		log={};(*Print[Dynamic[Column[log]]];*)

		(*get options*)
		pad = OptionValue[PaddingFactor];
		spfac = OptionValue[SplineSpacingFactor];
		field = OptionValue[SpectraFieldStrength];
		nuc = OptionValue[SpectraNucleus];
		shift = OptionValue[SpectraPpmShift];
		plots = OptionValue[SpectraOutputPlots];
		init = OptionValue[InitializeFit];
		readout = OptionValue[ReadoutType];

		(*set general parameters*)
		scale = 1000/Max[Abs[specIn]];
		nbas = Length[specBasisIn];

		
		(*-------------------------------------------------------------------*)
		(*create the basis functions fid and pad*)
		timeBasis = CashBasisTime[specBasisIn, pad, readout];

		(*pad and normalize the spectra spectra*)
		specFull = scale ApodizePadSpectra[specIn, PaddingFactor->pad, ReadoutType->readout];

		(*get the time and ppm axes*)
		{timeFull,ppmFull}=GetTimePpmRange[specFull,dtime,field,nuc];
		gyro = GetGyro[nuc,field];
		nsamp = Length[specFull];

		(*find the positions of the fit range*)
		{indSt,indEnd}=Flatten[Position[ppmFull,Nearest[ppmFull,#][[1]]]&/@{st,end}];

		(*logging of input parameters*)
		AppendTo[log,Style["Spectral Properties",Bold]];
		AppendTo[log,"    - Number of samples:          "<>ToString[nsamp]];
		AppendTo[log,"    - Gyro magnetic ratio:        "<>ToString[gyro]<>" Hz"];
		AppendTo[log,"    - Number of basis functions:  "<>ToString[nbas]];
		AppendTo[log,"    - Number of fit samples:      "<>ToString[indEnd-indSt+1]];
		AppendTo[log,""];

		
		(*-------------------------------------------------------------------*)
		If[init=!=Automatic,
			(*in not automatic it is an initial fit result*)
			{gami,epsi,{phi0i,phi1i},linei} = init;
			phii={phi0i,phi1i};
			,
			(*find initial linewidth and spec shift*)
			{gami,epsi,plLine}=EstimateLineWidth[{ppmFull,specFull},{lwvals,lwamsp},gyro,{st,end},plots];
			(*find initial phase estimate*)
			{{phi0i,phi1i},plShift}=EstimatePhaseShift[{ppmFull,specFull},{timeFull,timeBasis},{gami,epsi},gyro,{indSt,indEnd},readout,plots];
			(*define the initial line shape*)
			linei=0.5;

			(*logging of parameters*)
			AppendTo[log,Style["Estimating linewidth",Bold]];
			AppendTo[log,"    - spectral linewidth intit:   "<>ToString[Round[gami / gyro,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[gami,.0001]]<>" Hz"];
			AppendTo[log,"    - base spectra shift:         "<>ToString[Round[epsi,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[epsi gyro,.0001]]<>" Hz"];
			AppendTo[log,""];
			AppendTo[log,Style["Estimating phase correction",Bold]];
			AppendTo[log,"    - zeroth order phase:         "<>ToString[Round[phi0i,.0001]]<>" rad"];
			AppendTo[log,"                                  "<>ToString[Round[phi0i/Degree,.0001]]<>" deg"];
			AppendTo[log,"    - first order phase:          "<>ToString[Round[2Pi phi1i,.0001]]<>" rad/kHz"];
			AppendTo[log,"                                  "<>ToString[Round[phi1i,.0001]]<>" ms"];
			AppendTo[log,""];
		];

		
		(*-------------------------------------------------------------------*)
		(*get the spline basis fucntion paramters*)
		splineSpace=spfac Mean[Flatten[{gami}]]/gyro;
		cpn =Clip[ Round[Subtract@@Reverse[Sort[{st,end}]]/splineSpace],{4,Round[nsamp/10]}];

		(*logging of parameters*)
		AppendTo[log,Style["Estimating spline smoothness",Bold]];
		AppendTo[log,"    - spline spacing:             "<>ToString[splineSpace]<>" ppm"];
		AppendTo[log,"    - spline control points:      "<>ToString[cpn]];
		AppendTo[log,""];

		
		(*-------------------------------------------------------------------*)
		(*perform the first run minimization*)
		tfit1=0;
		If[NumberQ[gami]&&NumberQ[epsi]&&NumberQ[linei],
			(*define fit parameters and initialization*)
			Clear[phi0f,phi1f,gamf,epsf,linef,var,varf];

			(*define the fit varables*)
			var={{gamf,gami},{epsf,epsi},{linef,linei},{phi0f,phi0i},{phi1f,phi1i}};
			varf={gamf,epsf,{phi0f,phi1f}, linef};
			(*get the std of initial fit*)
			init={gami, epsi, phi1i};

			(*perform the fit*)
			{tfit1,fit1}=AbsoluteTiming[FindMinimum[FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},varf, init, Output->"Error", ReadoutType ->readout], var][[2]]];
							
			(*Get the fit results and output, wrap phi between -pi and pi*)
			sol={gami,epsi,phii,linei}={Clip[gamf,{1,500}],epsf,{2ArcTan[Tan[phi0f/2]],phi1f},linef}/.fit1;

			(*logging of parameters*)
			AppendTo[log,Style["Performing spectra First run",Bold]];
			AppendTo[log,"    - line shape:                 "<>ToString[Round[linei,.0001]]];
			AppendTo[log,"    - spectral linewidth:         "<>ToString[Round[gami/gyro,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[gami,.0001]]<>" Hz"];
			AppendTo[log,"    - base spectra shift:         "<>ToString[Round[epsi,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[epsi gyro,.0001]]<>" Hz"];
			AppendTo[log,"    - zeroth order phase:         "<>ToString[Round[phii[[1]],.0001]]<>" rad"];
			AppendTo[log,"                                  "<>ToString[Round[phii[[1]]/Degree,.0001]]<>" deg"];
			AppendTo[log,"    - first order phase:          "<>ToString[Round[2Pi phii[[2]],.0001]]<>" rad/kHz"];
			AppendTo[log,"                                  "<>ToString[Round[ phii[[2]],.0001]]<>" ms"];
			AppendTo[log,""];
		];

		(*redifine the spline spacing *)
		splineSpace = spfac Mean[Flatten[{gami}]]/gyro;
		cpn = Clip[Round[Subtract@@Reverse[MinMax[ppmFull]]/splineSpace],{4,Round[nsamp/10]}];

		(*make the output*)
		{fit,sigi}=FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},sol, Output->"Fit", ReadoutType ->readout];

		
		(*-------------------------------------------------------------------*)
		(*perform the second run minimization*)
		tfit2=0;
		If[OptionValue[FineTuneFit],
			(*prepare the fit parameters*)
			Clear[phi0f,phi1f,phi,gam,eps,line,var,varf,init];

			gamf=Table[Unique[gam],{i,1,nbas}];
			epsf=Table[Unique[eps],{i,1,nbas}];
			linef=If[OptionValue[FitLineShape],Table[Unique[line],{i,1,nbas}],Unique[line]];

			(*define the fit varables*)
			var=Join[MakeVars[gamf,gami,1],MakeVars[epsf,epsi,1],MakeVars[{phi0f,phi1f},phii,0],MakeVars[linef,linei,1]];
			varf={gamf,epsf,{phi0f,phi1f}, linef};
			(*get the std of initial fit*)
			init={gami,epsi,phii[[2]]};

			(*perform the minimization*)
			{tfit2,fit2}=AbsoluteTiming[FindMinimum[FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},varf, init, Output->"Error", ReadoutType ->readout], var][[2]]];

			(*get the solution and output, wrap phi between -pi and pi*)
			sol={gami, epsi, phii, linei}={Clip[gamf,{1,500}],epsf,{2ArcTan[Tan[phi0f/2]],phi1f},linef}/.fit2;

			(*recalculate the spline spacings*)
			splineSpace = spfac Mean[Flatten[{gami}]]/gyro;
			cpn = Clip[Round[Subtract@@Reverse[MinMax[ppmFull]]/splineSpace],{4,Round[nsamp/10]}];

			(*generate the output*)
			{fit,sigi}=FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},sol, Output->"Fit", ReadoutType ->readout];

			(*logging of parameters*)
			AppendTo[log,Style["Performing spectra: Second run",Bold]];
			AppendTo[log,"    - line shape:                 "<>ToString[Round[Mean[Flatten[{linei}]],.0001]]];
			AppendTo[log,"                                  "<>ToString[OptionValue[FitLineShape]]];
			AppendTo[log,"    - mean spectral linewidth:    "<>ToString[Round[Mean[gami]/gyro,.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[Mean[gami],.0001]]<>" Hz"];
			AppendTo[log,"    - mean base spectra shift:    "<>ToString[Round[Mean[epsi],.0001]]<>" ppm"];
			AppendTo[log,"                                  "<>ToString[Round[Mean[epsi] gyro,.0001]]<>" Hz"];
			AppendTo[log,"    - zeroth order phase:         "<>ToString[Round[phii[[1]],.0001]]<>" rad"];
			AppendTo[log,"                                  "<>ToString[Round[phii[[1]]/Degree,.0001]]<>" deg"];
			AppendTo[log,"    - first order phase:          "<>ToString[Round[2Pi phii[[2]],.0001]]<>" rad/kHz"];
			AppendTo[log,"                                  "<>ToString[Round[phii[[2]],.0001]]<>" ms"];
			AppendTo[log,""];
		]
	(*close timeing*)
	][[1]];

	
	(*-------------------------------------------------------------------*)
	(*apply the resuluts on the raw input*)
	{time,ppm}=GetTimePpmRange[specIn,dtime,field,nuc];
	gyro=GetGyro[nuc,field];
	nsamp=Length[specFull];

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
	AppendTo[log,"    - Total computation time:     "<>ToString[ttotal]<>" s"];
	AppendTo[log,"         - fit1 time:             "<>ToString[tfit1]<>" s"];
	AppendTo[log,"         - fit2 time:             "<>ToString[tfit2]<>" s"];
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
		specF, fidF, specBasisF, fit, ferr, gerr, perr, errorF, errorS, err, gami, epsi,
		sigi, phii, specBasis, fidBasisF, specFit, spline, readout, eerr
	},

	(*get the readout type*)
	readout = OptionValue[ReadoutType];

	(*give ouput either error or fit results*)
	Switch[OptionValue[Output],
		(* ------------ Calculate Error for minimization --------------*)
		"Error",

		(*----------- Aply parameters -----------*)
		(*apply phase to target instead of basis functions (faster)*)
		specF = PhaseShiftSpectra[spec, ppmFull, gyro, -phi];
		(*also convert to fid for error calcualtion*)
		fidF= ShiftedInverseFourier[specF, readout];
		(*select spectra range*)
		specF = specF[[indSt ;; indEnd]];

		(*generate basis spectra from time domain by applying gam, eps and line shape*)
		{fidBasisF,specBasisF} = BasisSpectraApply[{timeFull, timeBasis}, {gam, eps, f}, gyro, {indSt, indEnd}, readout];

		(*----------- Perform Fit and calculate errro -------------*)
		(*perform Fit of basis spectra*)
		fit = Quiet@NNLeastSquares[Join[Transpose[Re[specBasisF]],Transpose[Im[specBasisF]]], Join[Re[specF],Im[specF]]];			

		(*define errors fid and spectra*)
		errorS = specF - fit . specBasisF;
		errorF = fidF - fit . fidBasisF;

		(*Re and Im error normalized for number of points*)
		err = Mean[Join[Re[errorS]^2, Re[errorF]^2, Im[errorS]^2, Im[errorF]^2]];

		If[init === 0,
			(*constrain f between 0 and 1*)
			ferr = Total[ConFuncC[f, 0, 1, 4]];
			(*constrain phase to be small*)
			perr = ConFuncC[phi[[2]], -0.5, 0.5, 5];
			(*constrain gam to be positive*)
			gerr = Total[ConFuncC[gam, 1, 500, 3]];
			(*no initial values, minimize RMSE with f, gam and phase contraint*)
			err + ferr + gerr + perr
			,

			(*get initial values to constrain finetune fit, init ={gami, epsi, phii}*)
			(*constrain f between 0 and 1*)
			ferr = Total[ConFuncC[f, 0, 1, 4]];
			(*constrain lw gam to initial value*)
			gerr = Total[ConFuncSC[gam,init[[1]], 10, 2]];
			(*constrain shift eps to initial value*)
			eerr = Total[ConFuncSC[eps,init[[2]], 1, 2]];
			(*constrain phase to initial value*)
			perr = ConFuncC[phi[[2]]-init[[3]], -0.1, 0.1, 5];

			(*calculate error, minimize RMSE with f, gam, eps and phase contraint*)
			err + ferr + gerr + eerr + perr
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

BasisSpectraApply[{ppmFull_,timeFull_,timeBasis_},{gam_,eps_,phi_,f_},gyro_,{st_,end_},readout_]:=Block[{fidBasis, specBasis, vec, gamv, epsv, fv,func},
	(*make all parameters vectors*)
	{gamv, epsv, fv} = ConstantArray[1, {3,Length[timeBasis]}] {gam, eps, f};

	(*get the reaout function*)
	func = Switch[readout,"Fid",TimeShiftFid,"Echo",TimeShiftEcho];
	(*generate basis spectra from time domain by applying gam,eps and line shape*)
	fidBasis = MapThread[func[#1, timeFull, gyro, {#2, #3, #4}] &, {timeBasis, gamv, epsv, fv}];
	specBasis = ShiftedFourier[#, readout]&/@fidBasis;

	(*apply phase to the basis spectra, phi is only 0 if used in fitting and there also the fid is needed in all other cases only spectra is needed*)
	If[phi === 0, {fidBasis, specBasis[[All, st ;; end]]}, PhaseShiftSpectra[#, ppmFull, gyro, phi] & /@ specBasis]
]


(* ::Subsubsection::Closed:: *)
(*Estimate Line width*)


(*Function to estimate linewidth*)
EstimateLineWidth[{ppm_,spec_},{peaks_,amps_},gyro_,ran_,plot_:True]:=Block[{
	dppm,deltaf,corrf,wgth, s, max,corrf2,pts,pos,ppmC,maxf,block,line,sol,x,ppmCf,pl,lw,sft},

	(*define delta ppm and the delta function for correlation*)
	dppm=(ppm[[1]]-ppm[[2]]);
	deltaf=0 ppm;
	deltaf[[Flatten[Position[ppm,First@Nearest[ppm,#]]&/@peaks]]]=amps/Max[amps];

	(*perfomr correlation of spectra with delta function*)
	corrf=ListCorrelate[deltaf,Abs[spec],Round[Length[deltaf]/2],0];
	corrf=DivideNoZero[Length[corrf]corrf,Max[corrf]];

	(*Find max correlation and position*)
	s = Max[ppm]/4;
	wgth = 1/(Exp[(ppm^2/(2*s^2))]*(Sqrt[2*Pi]*s));
	max=Max[wgth corrf];
	pos=Position[wgth corrf,max][[1,1]];
	max=corrf[[pos]];

	(*contrain shift to 5ppm*)
	sft=-dppm(pos-(Length[ppm]/2.));
	sft=If[-5<sft<5,sft,0];
	pos=If[-5<sft<5,pos,0];

	(*find two points closest to FWHM*)
	maxf=0.5max;
	block=If[#<maxf,1,0]&/@corrf;
	pts={-FirstPosition[Reverse@block[[;;pos]],1][[1]]+pos+2,FirstPosition[block[[pos;;]],1][[1]]+pos-1};

	(*prevent not found errors in pts*)
	If[AllTrue[pts,NumericQ],
		line=Transpose@{{-1,0,1}+#,corrf[[{-1,0,1}+#]]}&/@pts;
		sol=Flatten[x/.Solve[Fit[#,{1,x},x]==maxf]&/@line];

		(*constrain the lw*)
		lw=gyro dppm(Subtract@@(Reverse@Sort@sol));
		lw=If[1<lw<200,lw,50];
		sol=If[1<lw<200,{1,2},sol];
		,
		{lw,sol}={50,{1,2}};
	];
	(*debugging plots*)
	pl=If[plot,
		ppmC=dppm(Range[Length[corrf]]-Length[corrf]/2);
		ppmCf=ListInterpolation[ppmC];
		FlipView[{Show[
				PlotSpectra[ppm,Max[Abs[spec]]deltaf,Method->"Abs",GridLineSpacing->5,PlotRange->{ran,Full},PlotColor->Red, PlotLabel->"Raw singal and calibration metabolites"],
				PlotSpectra[ppm,spec,Method->"Abs",GridLineSpacing->10,PlotRange->{ran,{0,Max[Abs[spec]]}}]
			],Show[
				PlotSpectra[ppmC,corrf,GridLineSpacing->5,PlotRange->{ran,Full}, PlotLabel->"Convolution signal"],
				ListLinePlot[{Transpose[{ppmCf[sol],{maxf,maxf}}],{{ppmC[[pos]],0},{ppmC[[pos]],max}}},PlotStyle->Directive[{Thick,Red}],ScalingFunctions->{"Reverse",Automatic}]
		]}], Null];

	(*calculate the estimated lw and shift*)
	{lw,sft,pl}
]


(* ::Subsubsection::Closed:: *)
(*EstimatePhaseShift*)


(*function to estimate phase form abs fitting of baiss spectra*)
EstimatePhaseShift[{ppm_,spec_},{time_,fids_},{gam_,eps_},gyro_,{st_,en_},readout_,plot_:True]:=Block[{
	phi1,sol1,phi2,sol2,specsC,fit,phi0f,phi1f,phi,specc,ran,pl,lim,specf,ppmf,func
	},

	specf=spec[[st;;en]];
	ppmf=ppm[[st;;en]];
	lim=.1;

	(*convert basis fids in spectra find function based on fid or echo*)
	func=Switch[readout,"Fid",TimeShiftFid,"Echo",TimeShiftEcho];
	specsC=Transpose[ShiftedFourier[func[#,time,gyro,{gam,eps,.5}], readout][[st;;en]]&/@fids];

	(*Fit absolute basis spectra to absolute spectrum*)
	fit=specsC . (NNLeastSquares[Abs[specsC],Abs[specf]]);
	(*minimize error with the target spectra*)
	sol1=Quiet@NMinimize[{PhaseError[ppmf,fit,specf,{phi0f,0 },gyro],-Pi<phi0f<Pi},{phi0f},MaxIterations->25][[2]];
	phi1={phi0f, 0}/.sol1;

	(*apply the zeroth order phase to the basis spectra*)
	specsC=Transpose[PhaseShiftSpectra[#,ppmf,gyro,phi1]&/@Transpose[specsC]];
	(*calculate the fit based on the imaginary part of the spectra*)
	fit=specsC . (NNLeastSquares[Re@specsC,Re@specf]);

	(*minimize error with the target spectra*)
	sol2=Quiet@NMinimize[{PhaseError[ppmf,fit,specf,{phi0f,phi1f },gyro],-Pi<phi0f<Pi,-lim<phi1f<lim},{phi0f,phi1f},MaxIterations->25][[2]];
	phi2={phi0f, phi1f}/.sol2;
	phi=phi1+phi2;

	(*debugging plots*)
	pl = If[plot,
		specc=PhaseShiftSpectra[fit,ppmf,gyro,phi2];
		ran={-1,1}Max[Abs[{specf,specc}]];
		FlipView[{
			Show[PlotSpectra[ppmf,specc,GridLineSpacing->5,PlotRange->{Full,ran}], PlotLabel->"Calibrated fit signal"],
			Show[PlotSpectra[ppmf,specf,GridLineSpacing->5,PlotRange->{Full,ran}], PlotLabel->"Raw signal"]
		}]
		,
		Null
	];

	{phi,pl}
];


PhaseError[ppm_,speci_,spect_,{phi0_?NumericQ,phi1_?NumericQ},gyro_]:=PhaseErrorC[ppm,speci,spect,phi0,phi1,gyro]

PhaseErrorC=Compile[{{ppm,_Real,1},{speci,_Complex,1},{spect,_Complex,1},{phi0,_Real,0},{phi1,_Real,0},{gyro,_Real,0}},Block[{spec},
	spec=Exp[-I(phi0+2Pi (phi1/1000) gyro ppm)]speci;
	Total[(Re[spect]-Re[spec])^2] + Total[(Im[spect]-Im[spec])^2]
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
		fun, plot, plot2, grid, gridS, or, rr, col, space, cols, cols2, pl1, pl2, lables, min, max, shift
	},

	(*get the plot range*)
	rr = OptionValue[PlotRange] /. Automatic -> Full;
	Switch[rr,
		{_, Full}, rr[[2]] = {Min[{Re[spec], Im[spec]}], Max[Abs[spec]]},
		Full, rr = {Full, {Min[{Re[spec], Im[spec]}], Max[Abs[spec]]}}
	];

	(*get gridline options*)
	gridS = OptionValue[GridLineSpacing];
	{min, max} = If[rr[[1]]===Full, MinMax[ppm], MinMax[rr[[1]]]];

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
		plot = Transpose[{Reverse[ppm + shift], #}] & /@ (#@spec & /@ fun);
		(*get the plot color*)
		col = If[OptionValue[PlotColor] === Automatic, 
			({{Gray, Thin}, {Red, Thin}, {Black}}[[-Length[fun] ;;]]), 
			OptionValue[PlotColor]
		];

		(*Make the plot*)
		{min, max} = {Ceiling@min, Floor@max};
		tck = Reverse@Range[min, max, Round[(max - min)/5]];
		ListLinePlot[plot, PlotStyle -> col, PlotRange -> rr, GridLines -> {grid, {0.}}, AspectRatio -> OptionValue[AspectRatio],
			ImageSize -> OptionValue[ImageSize], PlotLabel -> OptionValue[PlotLabel], (*ScalingFunctions -> {"Reverse", Automatic},*)
			Frame -> {{False, False}, {True, False}}, FrameStyle -> Directive[{Thick, Black}], FrameLabel -> {"PPM", None},
			LabelStyle -> {Bold, 14, Black}, PerformanceGoal->"Speed", MaxPlotPoints->Infinity, Filling->OptionValue[Filling],
			FrameTicks -> {Thread[{-tck, tck}], None}
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
		cols = Thread[{Append[ConstantArray[Black, Length[plot] - 1], Red], Thick}];
		lables = OptionValue[PlotLabels];

		(*make the plot*)
		pl1 = ListLinePlot[plot, Frame -> {{False, False}, {True, False}}, FrameStyle -> Thickness[.003], FrameTicksStyle -> Thickness[.003],
			PlotRange -> rr, PlotRangeClipping -> True, PlotStyle -> cols, ScalingFunctions -> {"Reverse", Automatic}, PlotLabel -> OptionValue[PlotLabel], 
			PlotLabels ->If[(OptionValue[Method] === "ReIm") || (lables === None), None, (Style[#, Black, Bold, 14] & /@ Append[lables, "All"])],
			GridLines -> {grid, None}, AspectRatio -> .5, ImageSize -> 1000, FrameLabel -> {"PPM", None}, LabelStyle -> {Bold, 14, Black}, 
			PerformanceGoal->"Speed", MaxPlotPoints->Infinity
		];

		(*make Im plot if needed*)
		If[OptionValue[Method] === "ReIm",
			plot2 = Transpose[{ppm + shift, #}] & /@ (Im[Append[spec, Total[spec]]] + space);
			cols2 = {Append[ConstantArray[Gray, Length[plot] - 1], Gray]};

			pl2 = ListLinePlot[plot2, Frame -> {{False, False}, {True, False}}, FrameStyle -> Thickness[.003], FrameTicksStyle -> Thickness[.003], 
				PlotRange -> rr, PlotRangeClipping -> True, PlotStyle -> cols2, ScalingFunctions -> {"Reverse", Automatic},
				PlotLabels -> (Style[#, Black, Bold, 14] & /@ Append[OptionValue[PlotLabels], "All"]), PlotLabel -> OptionValue[PlotLabel],
				GridLines -> {grid, None}, PlotRange -> rr, AspectRatio -> .5, 
				ImageSize -> 1000, FrameLabel -> {"PPM", None}, MaxPlotPoints->Infinity,
				LabelStyle -> {Bold, 14, Black}, PerformanceGoal->"Speed"
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

PlotFid[fid_?VectorQ, dwell_?NumberQ, opts : OptionsPattern[]] := PlotFid[GetTimeRange[fid, dwell], fid, opts]

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
	col = If[OptionValue[PlotColor] === Automatic, ({{Gray, Thin}, {Red, Thin}, {Black}}[[-Length[fun] ;;]]), OptionValue[PlotColor]];

	ListLinePlot[plot, PlotStyle -> col, PlotRange -> rr, GridLines -> {grid, {0.}}, AspectRatio -> OptionValue[AspectRatio],
		ImageSize -> OptionValue[ImageSize], PlotLabel -> OptionValue[PlotLabel], Frame -> {{False, False}, {True, False}},
		FrameStyle -> Directive[{Thick, Black}], FrameLabel -> {"time [s]", None}, LabelStyle -> {Bold, 14, Black}
	]
]


(* ::Subsection::Closed:: *)
(*PlotCSIData*)


Options[PlotCSIData]={PlotRange->Full}

SyntaxInformation[PlotCSIData] = {"ArgumentsPattern" -> {_, _, _., _.,OptionsPattern[]}}

PlotCSIData[datainp_, dw_?NumberQ, field_?NumberQ, nuc_?StringQ, opts:OptionsPattern[]] := PlotCSIData[datainp, {dw, GetGyro[nuc, field]}, opts]

PlotCSIData[datainp_, {dw_?NumberQ, field_?NumberQ, nuc_?StringQ}, opts:OptionsPattern[]] := PlotCSIData[datainp, {dw, GetGyro[nuc, field]}, opts]

PlotCSIData[datainp_, dw_?NumberQ, gyro_?NumberQ, opts:OptionsPattern[]] := PlotCSIData[datainp, {dw, gyro}, opts]

PlotCSIData[datainp_, {dw_?NumberQ, gyro_?NumberQ}, OptionsPattern[]] := Module[{
	data,datai,fun,nmax,dim,or,n,yran,dataPlot,maxPlot,maxAll,totAll,yrans,totPlot,colp,back,
	col,spec,coor,xdat,pmin,pmax,tdat,size,scale,leg,xmin,xmax,ymax,backScale, funs, c1, c2, lab},

	NotebookClose[plotwindow];

	data = N@datainp;

	If[! NumberQ[dw] && ! NumberQ[gyro],
		(*error dw and gyro are not OK*)
		$Failed
		,

		pan = Manipulate[
			(*prepare the data*)
			datai = fun@data;

			(*clip the range*)
			nmax = dim[[or]];
			n = Round[Clip[n, {1, nmax}]];
			{c1,c2} = Drop[dim[[;;-2]],{or}];

			(*get the correct data and ranges*)
			yran = {Min[{-0.5 Max[datai], Min[datai]}], 1.5 Max[datai]};
			dataPlot = Switch[or, 1, datai[[n]], 2, datai[[All, n]], 3, datai[[All, All, n]]];
			maxPlot = Switch[or, 1, maxAll[[n]], 2, maxAll[[All, n]], 3, maxAll[[All, All, n]]];
			totPlot = Switch[or, 1, totAll[[n]], 2, totAll[[All, n]], 3, totAll[[All, All, n]]];

			yrans = {Min[{-0.5 Max[dataPlot], Min[dataPlot]}], 1.5 Max[dataPlot]};

			Column[{
				Dynamic[
					spec = If[coor === {0, 0}, 
						0. data[[1, 1, 1]], 
						coor = {Clip[coor[[1]],{1,c1}], Clip[coor[[2]],{1,c2}]};
						Switch[or, 1, data[[n, coor[[1]], coor[[2]]]], 2, data[[coor[[1]], n, coor[[2]]]], 3, data[[coor[[1]], coor[[2]], n]]]
					];
					lab = Switch[or, 1, {n, coor[[1]], coor[[2]]}, 2, {coor[[1]], n, coor[[2]]}, 3, {coor[[1]], coor[[2]], n}];

					FlipView[{
						PlotSpectra[If[app,ApodizePadSpectra[spec],spec], {dw,gyro}, PlotRange -> {{pmin, pmax}, Full}, Method -> funs, PlotLabel->lab, ImageSize->Length[First[dataPlot]] size],
						PlotFid[tdat, ShiftedInverseFourier[spec], Method -> funs, PlotLabel->lab, ImageSize->Length[First[dataPlot]] size]
					}]
				]
				,
				"", 
				Style["Slice "<>ToString[n], Bold, Black, 24], 
				gridPlot = Grid[
					MapIndexed[(
						Item[
							EventHandler[
								(*Tooltip[*)
							(*the images*)
							Graphics[{Directive[{Thick, ColorData[{"DarkRainbow", "Reverse"}][#[[2]]]}], Line[Thread[{xdat(*[[;; ;; 2]]*), #1[[1]](*[[1, ;; ;; 2]]*)}]]},
								AspectRatio -> 1, ImageSize -> size, 
								Background -> If[back, Switch[backScale,"Max",GrayLevel[#[[2]]],"Total",GrayLevel[#[[3]]]], White],
								PlotRange -> {-{pmin, pmax}, Switch[scale, "Max", {Min[{-0.5 Max[#1[[1]]], 1.5 Min[#1[[1]]]}], 1.5 Max[#1[[1]]]}, "Full", yran, "Slice", yrans]}
							],
							(*(*the coordinate tooltip*)
							Switch[or, 1, {n, #2[[1]], #2[[2]]}, 2, {#2[[1]], n, #2[[2]]}, 3, {#2[[1]], #2[[2]], n}], 
							TooltipStyle -> {Directive[Black, Bold, Medium], Background -> White, CellFrameColor -> None, CellFrame -> None}],*)
							(*create the popup window*)
							"MouseClicked" :> (coor = #2)],
							(*highlight selected plot*)
							Background -> Dynamic[If[#2 == coor, Red, col]], Frame -> True, FrameStyle -> Dynamic[If[#2 == coor, Red, col]]
						]

					(*loop over all voxesl*)
					) &, RotateDimensionsLeft[{dataPlot, maxPlot, totPlot}], {2}]
					, Spacings -> {0.3,0.35}, Alignment -> Center, Background -> If[back, col, White], ItemSize -> Full, Frame -> All, FrameStyle -> If[back, col, White]
				],
				leg
			}, Alignment -> Center]

			(*active manipulate parametes*)
			, {{or, 1, "Orientation"}, {1 -> "Transversal", 2 -> "Coronal", 3 -> "Sagittal"}}
			, {{n, Ceiling[Length[data]/2], "Slice"}, 1, Dynamic[nmax], 1}
			, Delimiter
			, {{funs, "ReIm", "Function vox"}, {"Re", "Im", "ReIm", "Abs", "All"}}
			, {{app, False, "Apodize and Pad"}, {True,False}}
			, {{fun, Abs, "Function CSI"}, {Abs -> "Absolute", Re -> "Real", Im -> "Imaginary"}}
			, {{size, 40, "Plot size"}, {20 -> "Small", 40 -> "Medium", 60 -> "Large", 80 -> "Extra large"}}
			, Delimiter
			, {{pmin, xmin, "Min pmm"}, xmin, Dynamic[pmax - 1]}
			, {{pmax, xmax, "Min pmm"}, Dynamic[pmin + 1], xmax}
			, {{scale, "Max", "Plot scale"}, {"Max", "Full","Slice"}}
			, Delimiter
			, {{backScale, "Max", "Background"}, {"Max", "Total"}}
			, {{back, True, "Magnitude background"}, {True,False}}
			, {{col, Black, "Grid Color"}, (*ColorSlider*)
Button[
Dynamic[Graphics[{col, Rectangle[]}, ImageSize -> {20, 20}]],
new = SystemDialogInput["Color", #];
col = If[new === $Canceled, col, new]
, Background -> White, Method -> "Queued", FrameMargins -> 0, Appearance -> "Frameless"]&
}

			(* hidden manipulate paramterrs *)
			, {{coor, {0, 0}}, ControlType -> None}
							
			, {datai, ControlType -> None}, {dim, ControlType -> None}, {nmax, ControlType -> None}, {dataPlot, ControlType -> None} , {maxPlot, ControlType -> None}, {totPlot, ControlType -> None}
			, {yran, ControlType -> None}, {totAll, ControlType -> None}, {maxAll, ControlType -> None}, {ymax, ControlType -> None}, {xdat, ControlType -> None}, {tdat, ControlType -> None}
			, {xmin, ControlType -> None}, {xmax, ControlType -> None}, {spec, ControlType -> None}, {colp, ControlType -> None}

			, TrackedSymbols :> {or, n, fun, size, pmin, pmax, scale, col, back}
			, Initialization :> (
				dim = Dimensions[data];
				or = 1;
				nmax = dim[[1]];

				maxAll = Map[Max, Abs[data], {-2}];
				ymax = Max[maxAll];
				maxAll = maxAll/Max[maxAll];

				totAll = Map[Total, Abs[data], {-2}];
				totAll = totAll/Max[totAll];

				xdat = -GetPpmRange[data[[1, 1, 1]], dw, gyro];
				tdat = GetTimeRange[data[[1, 1, 1]], dw];

				{pmin,pmax} = {xmin, xmax} = If[OptionValue[PlotRange]===Full, MinMax[xdat], OptionValue[PlotRange]];

				leg = BarLegend[{{"DarkRainbow", "Reverse"}, {0, 100}}, LegendLayout -> "Row", LegendMarkerSize -> 400, LabelStyle -> Directive[Large, Black, Bold]];
			)

			, ControlPlacement -> Right
		];

		plotwindow = CreateWindow[DialogNotebook[{CancelButton["Close", Clear[data]; DialogReturn[]], pan}, WindowSize -> All, WindowTitle -> "Plot data window"]];
	]
]


(* ::Subsection:: *)
(*SpectraSimulator*)


 SpectraSimulator[]:=DynamicModule[{
use,vals,names, nsamp,bwS,fieldS,nuc,peakSel,dw,gyro,field,fidsT,time,ppm,grids
},

vals={{100,50,20,10},{20,30,40,100,0,75,500,250,40,60},{100}};
names={{"HDO","GlucoseD","GlutamateD","LactateD"},{"PE","PC","Piex","Piin","GPE","GPC","PCr","ATP","NAD","UDPG"},{"Na"}};
use=(0{{100,50,20,10},{20,30,40,100,0,75,500,250,40,60},{100}}+1)True;

grids={
{
{Row[{names[[1,1]]<>": ",InputField[Dynamic[vals[[1,1]]],FieldSize->3],Checkbox[Dynamic[use[[1,1]]]]}],Row[{names[[2,2]]<>": ",InputField[Dynamic[vals[[2,2]]],FieldSize->3],Checkbox[Dynamic[use[[1,2]]]]}],
Row[{names[[1,3]]<>": ",InputField[Dynamic[vals[[1,3]]],FieldSize->3],Checkbox[Dynamic[use[[1,3]]]]}],Row[{names[[2,4]]<>": ",InputField[Dynamic[vals[[2,4]]],FieldSize->3],Checkbox[Dynamic[use[[1,4]]]]}]
},{}},
{
{Row[{names[[2,1]]<>": ",InputField[Dynamic[vals[[2,1]]],FieldSize->3],Checkbox[Dynamic[use[[2,1]]]]}],Row[{names[[2,2]]<>": ",InputField[Dynamic[vals[[2,2]]],FieldSize->3],Checkbox[Dynamic[use[[2,2]]]]}],
Row[{names[[2,3]]<>": ",InputField[Dynamic[vals[[2,3]]],FieldSize->3],Checkbox[Dynamic[use[[2,3]]]]}],Row[{names[[2,4]]<>": ",InputField[Dynamic[vals[[2,4]]],FieldSize->3],Checkbox[Dynamic[use[[2,4]]]]}],
Row[{names[[2,5]]<>": ",InputField[Dynamic[vals[[2,5]]],FieldSize->3],Checkbox[Dynamic[use[[2,5]]]]}]},
{Row[{names[[2,6]]<>": ",InputField[Dynamic[vals[[2,6]]],FieldSize->3],Checkbox[Dynamic[use[[2,6]]]]}],
Row[{names[[2,7]]<>": ",InputField[Dynamic[vals[[2,7]]],FieldSize->3],Checkbox[Dynamic[use[[2,7]]]]}],Row[{names[[2,8]]<>": ",InputField[Dynamic[vals[[2,8]]],FieldSize->3],Checkbox[Dynamic[use[[2,8]]]]}],
Row[{names[[2,9]]<>": ",InputField[Dynamic[vals[[2,9]]],FieldSize->3],Checkbox[Dynamic[use[[2,9]]]]}],Row[{names[[2,10]]<>": ",InputField[Dynamic[vals[[2,10]]],FieldSize->3],Checkbox[Dynamic[use[[2,10]]]]}]}
},
{{Row[{names[[3,1]]<>": ",InputField[Dynamic[vals[[3,1]]],FieldSize->3],Checkbox[Dynamic[use[[3,1]]]]}]},{""}}
};

{nsamp,bwS,fieldS,nuc}={4096,10000,7,"31P"};

peakSel=Pick[Thread[names[[2]]->vals[[2]]],use[[2]]];
{dw,gyro,field,fidsT}=simulatedFid[{nsamp,bwS,fieldS,nuc},peakSel];

man = Manipulate[
simTrigger;

(*resample fid and generate noise*)
dw2=1./bw;
ni=Round[(te/1000)/dw];
ns=Round[Min[{ns,nsamp dw/dw2-ni}]];

fid1=ChangeDwellTimeFid[fidsT[[ni+1;;]],dw,dw2][[;;ns]];
sig=(Max[Abs[fid1]]/snr)/Sqrt[2];

If[lfid=!=Length[fid1]||sigi=!=sig,
lfid=Length[fid1];
sigi=sig;
noise1=Complex@@@RandomReal[NormalDistribution[0.,sig],{lfid,2}]
];

{timei,ppmi}=GetTimePpmRange[fid1,{dw2,field,nuc}];
tpi=timei+te/1000;

(*add shift and linewidth*)
fid2=If[timeshift,PhaseShiftFid[TimeShiftFid[fid1,tpi,gyro,{gam,eps,f}],ph0],fid1];

(*add phase to spectra*)
fid3=ShiftedInverseFourier[PhaseShiftSpectra[ShiftedFourier[fid2],ppmi,gyro,{ph0s,ph1s}]];

(*apply apodization and padding and figure out time of spectra to plot*)
{fid4,noise4}=If[ap==="None",{fid3,noise1}, ApodizeFid[#,ApodizationFunction->ap]&/@{fid3,noise1}];
{fid5,noise5}=If[pad===1,{fid4,noise4},PadFid[#,PaddingFactor->pad]&/@{fid4,noise4}];

(*prepare for plotting*)
fidTot=fid5+noise5;
noiseF=ShiftedFourier[noise5];
specF=ShiftedFourier[fid5];
specTot=noiseF+specF;

maxFid=Max[Abs[fidTot]];
maxSpec=Max[Abs[specTot]];

fidP=Switch[show,"both",fidTot,"signal",fid5,"noise",noise5];
specP=Switch[show,"both",specTot,"signal",specF,"noise",noiseF];

{time,ppm}=GetTimePpmRange[fidP,{dw2,field,nuc}];
timep=time+(te-ph1s)/1000;
{pmin,pmax}={Min[{pmin,Max[ppm]}],Max[{pmax,Min[ppm]}]};

Column[{
(*{Total@Re@fidP,Re@First@fidP,Total@Abs@ShiftedFourier@fidP,sig,simTrigger},
peakSel,*)
Show[
PlotFid[timep,fidP,Method->met,PlotRange->{{Min[{0,Min[timep]}],Max[timep]},{-1.5,1.5}maxFid},
	GridLineSpacing->0.05,ImageSize->psize],
If[!timeshift,Graphics[],ListLinePlot[Transpose@{timep,maxFid shiftFun[timep,eps,gyro]},PlotStyle->{Darker@Blue,Thin}]],
If[!timeshift,Graphics[],ListLinePlot[Transpose@{timep,maxFid relaxFun[timep,{gam,gam},f]},PlotStyle->{Darker@Orange,Dashed},PlotRange->Full]],
If[ap==="None",Graphics[],ListLinePlot[Transpose@{timei,maxFid apodizeFun[timei,ap]},PlotStyle->{Darker@Green,Dashed}]]
],
PlotSpectra[ppm,specP,PlotRange->{If[pran==="Automatic",Full,{pmin,pmax}],{-0.2,1.2}maxSpec},GridLineSpacing->5,AspectRatio->.5,Method->met, ImageSize->psize]
}]

,Style["Spectra simulation",14,Bold]
,{{nucs,2,"Nucleus"},{1->"2H",2->"31P",3->"23Na"}}
,{{fieldS,7,"Field"},{1,1.5,3,7}}
,Dynamic[Grid[grids[[nucs]]]]
,Button["Simulate fid",
simTrigger=False;
peakSel=Pick[Thread[names[[nucs]]->vals[[nucs]]],use[[nucs]]];
If[
AllTrue[use,#===False&],
Print["Select at least one peak"],
nuc={"2H","31P","23Na"}[[nucs]];
{dw,gyro,field,fidsT}=simulatedFid[{nsamp,bwS,fieldS,nuc},peakSel];
{time,ppm}=GetTimePpmRange[fidsT,{dw,field,nuc}];
simTrigger=True;
];,Method->"Queued"],

Delimiter,
OpenerView[{Style["acquisition parameters",12,Bold],Column[{
Control[{{bw,5000.,"bandwidth"},500,10000}],
Control[{{ns,512,"n samples"},64,1024,16}],
Control[{{te,0.,"echo time"},0,5}],
Control[{{snr,40.,"snr"},2,100}]
}]}],
{{bw,5000.},ControlType->None},
{{ns,512},ControlType->None},
{{te,0.},ControlType->None},
{{snr,40.},ControlType->None},

Delimiter,
OpenerView[{Style["linewidth, shift and phase",12,Bold],Column[{
Control[{{timeshift,True,"apply linewidth and shift"},{False,True}}],
Control[{{ph0,0,"phase [rad]"},-Pi,Pi}],
Control[{{eps,0,"shift [ppm]"},-5,5}],
Control[{{gam,40,"linewidth [Hz]"},0,200}],
Control[{{f,0.5,"lineshape"},0,1}],
Style["spectra phase",10,Bold],
Control[{{ph0s,0.,"0th order phase [rad]"},-Pi,Pi}],
Control[{{ph1s,0.,"1st order phase [ms]"},-10.,10.}]
}]}],
{{timeshift,True},ControlType->None},
{{ph0,0},ControlType->None},
{{eps,0},ControlType->None},
{{gam,40},ControlType->None},
{{f,0.5},ControlType->None},
{{ph0s,0.},ControlType->None},
{{ph1s,0.},ControlType->None},

Delimiter,
OpenerView[{Style["plotting options",12,Bold],Column[{
Control[{{pad,1,"padding"},1,4,0.5}]"",
Control[{{ap,"None","apodization"},{"None","Hanning","Hamming","Gaussian","Lorentzian","Voigt"},ControlType->SetterBar}],
Style["what to plot",10,Bold],
Control[{{show,"both","signal and noise"},{"both","signal","noise"}}],
Control[{{met,"ReIm","plot method"},{"All","Abs","Re","Im","ReIm"}}],
Style["options",10,Bold],
Control[{{pran,"Manual","ppm plot range"},{"Automatic","Manual"}}],
Control[{{pmin,8,"ppm min"},Dynamic[pmax],Dynamic[Max[ppm]]}],
Control[{{pmax,-18, "ppm max"},Dynamic[Min[ppm]],Dynamic[pmin]}],
Button["reset ppm range",{pmin,pmax}={8,-18}],
Control[{{psize,400,"plot size"},{200,300,400,500,600,700}}]
}]}],
{{pad,1},ControlType->None},
{{ap,"None"},ControlType->None},
{{show,"both"},ControlType->None},
{{pran,"Manual"},ControlType->None},
{{pmin,8},ControlType->None},
{{pmax,-18},ControlType->None},
{{met,"ReIm"},ControlType->None},
{{psize,400},ControlType->None},

(*hidden dynamic parameters*)
{dw2,ControlType->None},
{ni,ControlType->None},
{fid1,ControlType->None},
{fid2,ControlType->None},
{fid3,ControlType->None},
{fid4,ControlType->None},
{fid5,ControlType->None},
{noise1,ControlType->None},
{noise4,ControlType->None},
{noise5,ControlType->None},
{noiseF,ControlType->None},
{fidTot,ControlType->None},
{fidP,ControlType->None},
{specF,ControlType->None},
{specTot,ControlType->None},
{specP,ControlType->None},
{sig,ControlType->None},
{n1,ControlType->None},
{timei,ControlType->None},
{ppmi,ControlType->None},
{tpi,ControlType->None},
{timep,ControlType->None},
{maxFid,ControlType->None},
{maxSpec,ControlType->None},
{lfid,ControlType->None},
{sigi,ControlType->None},

{simTrigger,ControlType->None},

ControlPlacement->Left,
SaveDefinitions->True,
TrackedSymbols:>{bw,ns,te,timeshift,eps,gam,f,ph0,ph0s,ph1s,snr,pad,ap,show,pran,pmin,pmax,met,simTrigger}
];

NotebookClose[specsim];
specsim=CreateWindow[DialogNotebook[{CancelButton["Close",DialogReturn[]],man},
WindowSize->All,WindowTitle->"Spectra simulator"]];

];

simulatedFid[{nsamp_,bw_,field_,nuc_},peaks_]:=Block[{dw,gyro,names,fids,specs,table},
dw=1./bw;
gyro=GetGyro[field,nuc];
{names,fids,specs,table}=GetSpectraBasisFunctions[peaks[[All,1]],
BasisSequence->{"PulseAcquire",0},SpectraSamples->nsamp,SpectraBandwidth->bw,
SpectraPpmShift->0,SpectraFieldStrength->field,SpectraNucleus->nuc];
{dw,gyro,field,(names/.peaks) . fids}
]

apodizeFun[time_,apM_:"Hanning"]:=Block[{length,app,xdat,xmax},
xdat=time;
length=Length@time;
xmax=Max[Abs[xdat]];
app=Switch[apM,
"Hanning",0.5+0.5 Cos[xdat Pi/xmax],
"Hamming",0.54+0.46 Cos[xdat Pi/xmax],
"Gaussian",Exp[-(3./xmax) xdat],
"Lorentzian",Exp[-(2./xmax)^2 xdat^2],
"Voigt",0.5 Exp[-(3./xmax) xdat]+0.5 Exp[-(2./xmax)^2 xdat^2]
];
app=app/Max[app];
app
]

relaxFun[time_,{gam1_,gam2_},f_]:=Exp[-(f gam1 time + (1-f)gam2^2 time^2)]
shiftFun[time_,eps_,gyro_]:=Re@Exp[2 Pi eps gyro I time] 


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

MakeSpectraResultPlot[ppmF_, specF_, {fit_, basisFit_}, names_, ppmran_] := Block[{
	sp, specFit, resTotPl, errPl, fitPl, resPl, outPl, pran, pmax,
	lab1, lab2, resfitRI, resfit, resBasPl, met},

	sp = 2;
	met = "ReIm";
	specFit = fit . basisFit;
	pmax = Max[Abs[specFit], Abs[specF]];
	pran = {-pmax, pmax};

	resTotPl = Column[{
		FlipView[errPl = {
			PlotSpectra[ppmF, specF - specFit, Method -> met, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp, PlotLabel->"Fit error full"],
			PlotSpectra[ppmF, specF - specFit, Method -> met, PlotRange -> {ppmran, Full}, GridLineSpacing -> sp, PlotLabel->"Fit error scaled"]
		}],
		FlipView[fitPl = {
			Show[
				PlotSpectra[ppmF, specF, Method -> met /. "ReIm" -> "Re", PlotColor -> Red, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp, PlotLabel->"Real signal"],
				PlotSpectra[ppmF, basisFit[[1]], Method -> met /. "ReIm" -> "Re", PlotColor -> Green, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, specFit, Method -> met /. "ReIm" -> "Re", PlotColor -> Black, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp]
			],
			Show[
				PlotSpectra[ppmF, specF, Method -> met /. "ReIm" -> "Im", PlotColor -> Red, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp, PlotLabel->"Imaginary signal"],
				PlotSpectra[ppmF, basisFit[[1]], Method -> met /. "ReIm" -> "Im", PlotColor -> Green, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, specFit, Method -> met /. "ReIm" -> "Im", PlotColor -> Black, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp]
			]
		}]
		,
		FlipView[resPl = {
			PlotSpectra[ppmF, specF, Method -> met, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp, PlotLabel->"Raw signal"],
			PlotSpectra[ppmF, specFit, Method -> met, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp, PlotLabel->"Fitted signal"]
		}]
	}, Alignment -> Center];

	lab1 = Style[#, Bold, Black, Large] & /@ {"Error scaled to signal", "Error scaled to Max"};
	lab2 = Style[#, Bold, Black, Large] & /@ {"Fit and signal Re", "Fit and signal Im"};
	outPl = Column[Flatten[{Thread[{lab1, errPl}], Thread[{lab2, fitPl}]}], Alignment -> Center];

	resBasPl = FlipView[{
		resfitRI = PlotSpectra[ppmF, fit basisFit, Method -> "ReIm", PlotColor -> Red, SpectraSpacing -> 0.2, GridLines -> {}, GridLineSpacing -> sp, 
			PlotLabels -> Prepend[names, "spline"],PlotRange->{ppmran,Full}, PlotLabel->"Fitted basis spectra real and imaginary"],
		resfit = PlotSpectra[ppmF, fit basisFit, Method -> "Abs", PlotColor -> Red, SpectraSpacing -> 0.2, GridLines -> {}, GridLineSpacing -> sp, 
			PlotLabels -> Prepend[names, "spline"],PlotRange->{ppmran,Full}, PlotLabel->"Fitted basis spectra absolute"]
	}];

	{resTotPl, resBasPl, {errPl, fitPl, resPl, outPl, resfitRI, resfit}}
]


(* ::Subsection::Closed:: *)
(*CSIInterface31P*)


Options[CSIInterface] = {SpectraFieldStrength -> 7,SpectraNucleus->"31P"};

SyntaxInformation[CSIInterface] = {"ArgumentsPattern" -> {_., _., OptionsPattern[]}};

CSIInterface[opts : OptionsPattern[]] := CSIInterface["", {0, 0}, opts]

CSIInterface[tei_?NumberQ, bwi_?NumberQ, opts : OptionsPattern[]] := CSIInterface["", {tei, bwi}, opts]

CSIInterface[file_?StringQ, opts : OptionsPattern[]] := CSIInterface[file, {0, 0}, opts]

CSIInterface[file_?StringQ, {tei_?NumberQ, bwi_?NumberQ}, OptionsPattern[]] := Module[{
	f, te, bw, nuc, field, gyro, names, met, metSel, metRef, method, plot, kload, rec, spectraC, dec, line, fine, den, z, y, x, 
	sphase, status, statusP, kspace, noise,	header, type, ham, spectra, spectraR, spec, proc, shift, fids, specs, 
	table, fit, basisFit, errorFit, pars, log, plots, specf, fitted, xm, ym, zm, dn, dc, mr, dw, nsamp, filt, teu,
	fileSave, spectraPlot, lab, fovz, fovy, fovx, coils, ncoils, teE,
	statPart, loadPart, reconPart, plotpart, fitPart, closePart
	},

	NotebookClose[interfaceWindow];
	NotebookClose[reswindow];
	NotebookClose[plotwindow];

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
	kload = rec = sphase = fitted = proc = dn = dc = False;
	mr = ""; 
	spectraC = shift = coils= 0;

	(*initialize fit coors*)
	z = y = x = fovz = fovy = fovx = zm = ym = xm = 0;

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
			NotebookClose[reswindow];
			NotebookClose[plotwindow];
			status = "Loading data"; statusP = True;
			rec = proc = dn = dc = sphase = False;
			{{kspace, noise}, {header, type}} = ReadListData[f, False];
			If[filt,
				ham = Total@Unitize[Abs[OrderKspace[kspace, type, {"N_aver", "N_kz", "N_ky", "N_kx"}][[1]]]];
				kspace = OrderKspace[TotalType[kspace, type, {"N_aver"}], {"N_chan", "N_samp", "N_kz", "N_ky", "N_kx"}][[1]];
				,
				kspace = OrderKspace[MeanType[kspace, type, {"N_aver"}], {"N_chan", "N_samp", "N_kz", "N_ky", "N_kx"}][[1]];
				ham = MakeHammingFilter[Dimensions[kspace][[-3;;]]];
				kspace = Map[ham #&,kspace,{2}];
			];
			coils = ncoils = Range[Length[kspace]];
			kload = True;
			status = "Done loading data!"; statusP = False;
			,
			Background -> Dynamic[If[kload, RGBColor[0.86, 0.97, 0.77], Automatic]], Enabled -> Dynamic[FileExistsQ[f]], Method -> "Queued", ImageSize -> 175
		],Dynamic[f]},
		{TextCell[" Filter data  "], SetterBar[Dynamic[filt], {True -> " k-space weighting ", False -> " Hamming "}]}
	};

	(*------------ ACQUISITION / RECONSTRUCTION -------------*)
	reconPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Acquistion parameters", Bold, 14]},

		(*acquistion paramters*)
		{TextCell["Echo time [ms] "], InputField[Dynamic[te], Number]},
		{TextCell["Bandwidth [Hz] "], InputField[Dynamic[bw], Number]},

		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Data reconstruction", Bold, 14]},

		(*reconstruct button*)
		{Button["Reconstruct data",
			NotebookClose[reswindow];
			NotebookClose[plotwindow];
			status = "Reconstructing using " <> method <> ""; statusP = True;
			rec = proc = dn = dc = sphase = False; mr = "";
			mr = method;
			dw = 1./bw;
			teu = te/1000;
			coils = If[coils==={},All,coils];
			spectra = spectraR = CoilWeightedReconCSI[kspace[[coils]], noise[[coils]], header, Method -> method];
			rec = True; 
			If[den, spectra = DenoiseCSIdata[spectra]; proc = dn = True; ];
			If[dec, spectra = DeconvolveCSIdata[spectra, ham]; proc = dc = True;];
			{zm, ym, xm, nsamp} = Dimensions[spectra];

			(*find shift*)
			spec = Mean@Flatten[spectra, 2];
			shift = FindSpectraPpmShift[spec,{1./bw,gyro}, {{0, -2.52, -7.56, -16.15}, {2, 1, 1, 1}}]; 
			status = "Done reconstructing!"; statusP = False;
			,
			Background -> Dynamic[If[rec, RGBColor[0.86, 0.97, 0.77], Automatic]], Enabled -> Dynamic[kload && bw =!= 0], Method -> "Queued", ImageSize -> 175
		],
		Row[{TextCell[" Denoise data  "], Checkbox[Dynamic[den]], TextCell["   Deconvolve data  "], Checkbox[Dynamic[dec]]}, ""]},
		{TextCell[" Select coils  "], Dynamic[If[coils===0,TextCell[""], TogglerBar[Dynamic[coils], ncoils]]]},

		(*select method*)
		{TextCell[" Reconstruction method  "], SetterBar[Dynamic[method], {"Roemer" -> " Roemer ", "WSVD" -> " WSVD "}]},

		(*report*)
		{Dynamic[TextCell[If[rec, "PCr shift is " <> ToString[Round[shift, .01]] <> " ppm ", ""]]], Dynamic[TextCell[If[rec, "Reconstructed using: " <> mr <> If[dn, ", denoised", ""] <> If[dc, ", deconvolved", ""], ""]]]}
		};

	(*------------ PLOTTING -------------*)
	plotpart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Data Plotting", Bold, 14]},

		(*plotting data button*)
		{Button["Plot CSI data", 
			NotebookClose[reswindow];
			NotebookClose[plotwindow];
			spectraPlot = Switch[plot,
				"Raw", spectraR,
				"Proc", spectra,
				"Cor",
				If[!sphase,
					status = "Correcting phase of spectra"; statusP = True;
					spectraC = Map[PhaseCorrectSpectra[ApodizePadSpectra[ShiftSpectra[#, {dw, gyro}, -shift]], dw, teu, gyro, {10, -20}]&, spectra, {-2}];
					status = "Done phase correcting spectra!"; statusP = False; sphase = True
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
				If[!sphase,
					status = "Correcting phase of spectra"; statusP = True;
					spectraC = Map[PhaseCorrectSpectra[ApodizePadSpectra[ShiftSpectra[#, {dw, gyro}, -shift]], dw, teu, gyro, {10, -20}, True]&, spectra, {-2}];
					status = "Done phase correcting spectra!"; statusP = False;	sphase = True
				];
				teE=0.;
				spectraC
			];
			status = "Saving the CSI data";
			fileSave = SystemDialogInput["FileSave"]; 
			lab = mr <> If[plot==="Raw","",If[dn, ", denoised", ""] <> If[dc, ", deconvolved", ""]]<>If[plot ==="Cor"," ,phase corrected",""];
			If[fileSave =!= $Canceled, ExportSparSdat[fileSave, spectraPlot, {bw, teE}, {gyro, nuc}, {fovz, fovy, fovx}]];
			(*saving make better!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
			,
			Enabled -> Dynamic[rec&&fovy>0&&fovz>0&&fovx>0], Method -> "Queued", ImageSize -> 175],
		Row[{
			{fovz, fovy, fovx} = Round[{fovz, fovy, fovx}];
			TextCell["  vox z [mm] "], InputField[Dynamic[fovz, (fovz = Round[#]) &], Number, FieldSize -> 3, ContinuousAction -> True, Background -> Dynamic[If[fovz>0, White, RGBColor[1, 0.9, 0.9]]]],
			TextCell["   vox y [mm] "], InputField[Dynamic[fovy, (fovy = Round[#]) &], Number, FieldSize -> 3, ContinuousAction -> True, Background -> Dynamic[If[fovy>0, White, RGBColor[1, 0.9, 0.9]]]],
			TextCell["   vox x [mm]  "], InputField[Dynamic[fovx, (fovx = Round[#]) &], Number, FieldSize -> 3, ContinuousAction -> True, Background -> Dynamic[If[fovx>0, White, RGBColor[1, 0.9, 0.9]]]]
		}, " "]}
	};

	(*------------ FITTING -------------*)
	fitPart = {
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Spectra Fitting", Bold, 14]},
		(*fitting button*)
		{Button["Fit selected spectra",
			NotebookClose[reswindow];
			NotebookClose[plotwindow];
			(*basis spectra*)
			status = "Generating basis spectra"; statusP = True;
			{names, fids, specs, table} = GetSpectraBasisFunctions[metSel, {"ATP"}, BasisSequence -> {"PulseAcquire", teu}, 
				SpectraSamples -> nsamp, SpectraBandwidth -> bw, SpectraPpmShift -> 0, SpectraFieldStrength -> field, SpectraNucleus -> nuc];
			status = "Done generating basis spectra!"; statusP = False;

			(*fitting*)
			status = "Start fitting spectra"; statusP = True;
			specf = ShiftSpectra[spectra[[z, y, x]], {dw, gyro}, -shift];
			{fit, basisFit, errorFit, pars, log, plots} = FitSpectra[specs, specf, {10, -20}, dw, { {0, -2.52, -7.56, -16.15}, {2, 1, 1, 1}},
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
		{TextCell["Fitting options "], Row[{TextCell[" Finetune fit  "], Checkbox[Dynamic[fine]], TextCell["   Fit line shape  "], Checkbox[Dynamic[line]]}, ""]},
		(*select metabolites*)
		{TextCell["Select metabolites to fit "], TogglerBar[Dynamic[metSel], met]},

		(*------------ RESULTS -------------*)
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], TextCell["Fitting Results", Bold, 14]},

		(*show fit results button*)
		{TextCell["Reference metabolite "], Dynamic[SetterBar[Dynamic[metRef], names]]},
		{Button["Show results",
			NotebookClose[reswindow];
			NotebookClose[plotwindow];
			status = "Generating fit results"; statusP = True;
			SpectraFitResult[specf, {fit, basisFit}, teu, {dw, gyro}, {pars, names, metRef, log}, plots, PlotRange -> {10, -20}];
			status = "Done with generating fitting results!"; statusP = False;
			,
			Enabled -> Dynamic[fitted], Method -> "Queued", ImageSize -> 175]
		}
	};

	(*------------ CLOSING -------------*)
	closePart ={
		{Item[TextCell[""], Background -> Automatic], SpanFromLeft},
		{TextCell[""], CancelButton["Close GUI", DialogReturn[
			NotebookClose[reswindow];
			NotebookClose[plotwindow];
			Clear[kspace, spectra]
			], ImageSize -> 175]},{TextCell[""]
		},
		{TextCell[""]}
	};

	grid = Join[statPart, loadPart, reconPart, plotpart, If[nuc==="31P",fitPart,{}], closePart];

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

SpectraFitResult[specf_, {fit_, basisFit_}, te_, {dw_, gyro_}, {pars_, names_, metRef_, log_}, ops : OptionsPattern[]] := SpectraFitResult[specf, {fit, basisFit}, te, {dw, gyro}, {pars, names, metRef, log}, {}, ops]

SpectraFitResult[specf_, {fit_, basisFit_}, te_, {dw_, gyro_}, {pars_, names_, metRef_, log_}, plots_, OptionsPattern[]] := Block[{resTab, ppm, phase, resTotPl, resBasPl, specFit, fidf, fidFit, time, comFidPlot, specPlot, fitPlot, comSpecPlot, logging, ran},
	NotebookClose[reswindow];
	resTab = FitSpectraResultTable[fit, pars, names, metRef];
	ran = OptionValue[PlotRange];

	(*make result plot*)
	ppm = GetPpmRange[specf, dw, gyro];
	phase = -pars[[3]];
	{resTotPl, resBasPl} = MakeSpectraResultPlot[ppm, specf, {fit, basisFit}, names, ran][[1 ;; 2]];

	(*make fitted spectra and fids*)
	specFit = fit . basisFit;
	fidf = ShiftedInverseFourier[specf];
	fidFit = ShiftedInverseFourier[specFit];

	(*make fid plot*)
	time = GetTimeRange[fidf, dw];
	comFidPlot = CompareFidFitPlot[time, fidf, fidFit];

	(*correct fitted spectra*)
	specPlot = CorrectTESpec[ApodizePadSpectra[PhaseShiftSpectra[specf, ppm, gyro, phase]], dw, te, gyro, ran];
	fitPlot = CorrectTESpec[ApodizePadSpectra[PhaseShiftSpectra[specFit, ppm, gyro, phase]], dw, te, gyro, ran];

	(*make spectra plot*)
	ppm = GetPpmRange[specPlot, dw, gyro];
	comSpecPlot = CompareSpectraFitPlot[ppm, specPlot, fitPlot, ran];

	(*make logging*)
	logging = MakeLogging[log];
	(*make result window*)
	reswindow = CreateDialog[{
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

ImportSparSdat[fspar_,fsdat_]:=Block[{nums,head,depth,row,nsamp,x,y,z,bw,te,nuc,gyro,field,fids,specs},
	(*read raw data and header*)
	nums=FromVaxD[BinaryReadList[fsdat, "UnsignedInteger32"]];
	head=ParseSpar[Import[fspar, "Lines"]];

	(*get header infor*)
	{bw,te,nuc,gyro}={"sample_frequency","spectrum_echo_time","nucleus","synthesizer_frequency"}/.head;
	gyro = gyro 10^-6;
	field=Round[gyro /GyromagneticRatio[nuc],.1];

	{depth,row,nsamp,x,y,z}={"num_dimensions","spec_num_row","dim1_pnts","dim2_pnts","dim3_pnts","dim4_pnts"}/.head;
	(*partition the data to the correct size*)
	fids=Switch[depth,
		2,
		If[row==1,
			(*single voxel*)
			RotateDimensionsRight[Fold[Partition,nums,{2}]],
			(*dynamic series*)
			RotateDimensionsRight[Fold[Partition,nums,{2,nsamp}]]
		],
		3,
		RotateDimensionsRight[Fold[Partition,nums,{2,nsamp,x}]],
		4,(*3DCSI*)
		RotateDimensionsRight[Fold[Partition,nums,{2,nsamp,x,y}]]
	];
	fids=fids[[1]]+fids[[2]]I;

	(*convert fid to spec*)
	specs=If[ArrayDepth[fids]>1, Map[ShiftedFourier,fids,{-2}],	ShiftedFourier[fids]];

	(*ouput data*)
	{specs,{bw,te},{gyro,nuc},head}
]


(* ::Subsubsection::Closed:: *)
(*ParseSpar*)


ParseSpar[spari_]:=Block[{textIn,values,lab,val,spar},
	spar=Select[spari,#=!=""&];
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

ExportSparSdat[file_, specs_, {bw_, te_}, {gyro_, nuc_}, vox_, opts:OptionsPattern[]]:=Block[{fidsOut,numsOut,fileOut,datOut,headOut},
	(*export data*)
	fidsOut=Map[ShiftedInverseFourier,specs,{-2}];
	numsOut=ToVaxD[Flatten[RotateDimensionsLeft[{Re@fidsOut,Im@fidsOut}]]];
	fileOut=file<>".SDAT";
	If[FileExistsQ[fileOut],DeleteFile[fileOut]];
	BinaryWrite[fileOut,numsOut,"UnsignedInteger32"];
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
		(*get the fracqtion*)
		fracBin=PadLeft[IntegerDigits[Round[16777216.(numA/(2.^exp)-0.5)],2],23];
		(*output the integer*)
		FromDigits[Join[fracBin[[8;;-1]],signBin,expBin,fracBin[[1;;7]]],2]
	]
],RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


(* ::Subsubsection::Closed:: *)
(*MakeSpar*)


Options[MakeSpar] = Options[ExportSparSdat];

MakeSpar[specs_, {bw_, te_}, {gyro_, nuc_}, vox_, OptionsPattern[]]:=Block[{
		dimzO,dimyO,dimxO,nsampO,gyroO,nucO,bwO,teO,nameO, hf, ps, id,
		text,lab,filHeader,fixedHeader, vals, head,row,depth
	},
	(*swith between data dimensions*)
	depth=ArrayDepth[specs];
	Switch[depth,
		1,(*single voxel*)
		nsampO=Length[specs];
		row=dimzO=dimyO=dimxO=1;
		,
		2,(*dynamic data*)
		{row,nsampO}=Dimensions[specs];
		dimzO=dimyO=dimxO=1;
		,
		3,(*2D CSI*)
		{dimyO,dimxO,nsampO}=Dimensions[specs];
		row=dimxO dimyO;
		dimzO=1;
		,
		4,(*3D CSI*)
		{dimzO,dimyO,dimxO,nsampO}=Dimensions[specs];
		row=dimxO dimyO;
	];

	(*manditory input paramters*)
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
		"!   S15/ACS: set of *.SPAR and *.SDAT files is created, (dataformat: VAX CPX floats)","!--------------------------------------------------------------------","!---------------------------------------------",
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
		"TSI_factor","spectrum_echo_time","spectrum_inversion_time","image_chemical_shift","resp_motion_comp_technique","de_coupling","equipment_sw_verions","placeholder1","placeholder2"
	};

	(*from input*)
	filHeader={
		(*general acquistion names*)
		"patient_position"->Switch[hf, 0, "\"head_first\"", 1, "\"feat_first\""],
		"patient_orientation" ->Switch[ps, 0, "\"supine\"", 1, "\"prone\""],
		(*get from input window*)
		"examination_name"->"Generated by QMRITools",
		"patient_name"->nameO,
		"phase_encoding_fov"->vox[[2]] dimyO,(*mm fov in freq*)
		"slice_thickness"->vox[[1]] dimzO,(*mm fov in phase*)
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
		"samples" ->nsampO,
		"spec_num_col"->nsampO,
		"dim1_pnts"->nsampO,
		"dim1_step"-> 1./bwO,(*s*)
		"spec_col_upper_val"->(nsampO-1)(1./bw),(*s*)
		(*data dimensions*)
		"dim2_pnts"->dimxO,
		"dim3_pnts"->dimyO,
		"nr_of_slices_for_multislice"->dimzO,
		"nr_phase_encoding_profiles"->dimxO,
		"nr_of_phase_encoding_profiles_ky"->dimyO,
		"phase_encoding_enable"-> If[depth<=2,"\"no\"","\"yes\""],
		"dim4_pnts" ->dimzO
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
		"equipment_sw_verions"->"QMRITools","placeholder1"->"","placeholder2"->""
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
