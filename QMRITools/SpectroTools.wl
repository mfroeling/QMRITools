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
"ReadJMRUI[file] read a jMRUI spectrum file. 
Output is the {time, spec, {begintime, samplingInterval}}."


PadFid::usage = 
"PadFid[fid] pads the fid with zeros to double its length."

ApodizeFid::usage = 
"ApodizeFid[fid] performs apodization on the fid. The apodization function is set with the option ApodizationFunction."

ApodizePadFid::usage = 
"ApodizePadFid[fid] performs apodization on the fid and pads the fid with zeros to double its length"

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


PhaseCorrectSpectra::usage =
"PhaseCorrectSpectra[spect] performs 0th order phase correction of the spectra.
PhaseCorrectSpectra[spect, ppm, gyro] performs 0th and 1st order phase correction of the spectra.
PhaseCorrectSpectra[spect, ppm, gyro, phi1] performs 0th order phase correction of the spectra with a known 1st order phase.

0th order phase is in radians and the 1st order phase is in ms."


PhaseShiftSpectra::usage = 
"PhaseShiftSpectra[spectra, phi0] aplies the 0th order phase phi0 to the spectra. 
PhaseShiftSpectra[spectra, ppm, gyro, phi1] aplies the 1st order phase phi1 to the spectra. The ppm can be obtained using GetPpmRange and gyro with GetGyro. 
PhaseShiftSpectra[spec, ppm, gyro, {phi0, phi1}] aplies the 0th and 1st order phases {phi0, phi1} to the spectra. The ppm can be obtained using GetPpmRange and gyro with GetGyro.

The 0th order phase phi0 is in radians and the 1st order phase phi1 is in ms."

TimeShiftFid::usage = 
"TimeShiftFid[fid, time, gam] aplies a linebroadening with linewidth gam and a Voight lineshape to the fid. The time can be obtained using GetTimeRange.
TimeShiftFid[fid, time, {gam, f}] aplies a linebroadening with linewidth gam and a custom lineshape f to the fid (f=0, \"Gaussinan\", f=1 \"Laurentian\").
TimeShiftFid[fid, time, gyro, {gam, eps}] aplies a linebroadening with linewidth gam to the fid and a phase eps that results in eps ppm shift of the spectra. The gyro can be obtained with GetGyro
TimeShiftFid[fid, time, gyro, {gam, eps, f}] aplies a linebroadening with linewidth gam using a custom lineshape f to the fid and a phase eps that results in eps ppm shift of the spectra.

The linewidth gam is given in ms and the spectra shift eps is given in ppm."

ChangeDwellTimeFid::usage = 
"ChangeDwellTimeFid[fid, dt, dtnew] changes the sampleling time of an fid from dwelltime dt to dwelltime dtnew."


(* ::Subsection:: *)
(*Options*)


ApodizationFunction::usage = 
"ApodizationFunction is an options for ApodizeFid, ApodizeSpectra, ApodizePadFid, and ApodizePadSpectra. Values can be \"Han\", \"Ham\", \"Gaus\", \"Laur\", and \"GausLaur\"."


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


(* ::Subsection:: *)
(*PhaseCorrectSpectra*)


SyntaxInformation[PhaseCorrectSpectra]={"ArgumentsPattern"->{_, _., _., _.}}

PhaseCorrectSpectra[spect_] := Block[{phi0, sol},
	sol = Quiet[Last[FindMinimum[PhaseCorrectError[spect, phi0], {phi0, 0}]]];
	Exp[-I phi0] spect /. sol
]

PhaseCorrectSpectra[spect_, ppm_, gyro_] := Block[{phi0, phi1, sol},
	sol = Quiet[Last[FindMinimum[PhaseCorrectError[ppm, spect, phi0, phi1, gyro], {{phi0, 0}, {phi1, 0}}]]];
	Exp[-I (phi0 + 2 Pi (phi1/1000)  gyro ppm)] spect /. sol
]

(*correct zeroth phase with known first order phase*)
PhaseCorrectSpectra[spect_, ppm_, gyro_, phi1_] := Block[{phi0, sol},
	sol = Quiet[Last[FindMinimum[PhaseCorrectError[ppm, spect, phi0, phi1, gyro], {phi0, 0}]]];
	Exp[-I (phi0 + 2 Pi (phi1/1000)  gyro ppm)] spect /. sol
]


PhaseCorrectError[speci_, phi0_?NumericQ] := PhaseCorrectErrorC1[speci, phi0]

PhaseCorrectError[ppm_, speci_, phi0_?NumericQ, phi1_?NumericQ, gyro_] := PhaseCorrectErrorC2[ppm, speci, phi0, phi1, gyro]

PhaseCorrectErrorC1 = Compile[{{speci, _Complex, 1}, {phi0, _Real, 0}}, Block[{spec},
	spec = Exp[-I phi0] speci;
	Total[(Abs[speci] - Re[spec])^2] (*+ Total[Im@spec]^2*)
], RuntimeOptions -> "Speed", Parallelization -> True];

PhaseCorrectErrorC2 = Compile[{{ppm, _Real, 1}, {speci, _Complex, 1}, {phi0, _Real, 0}, {phi1, _Real, 0}, {gyro, _Real, 0}}, Block[{spec},
	spec = Exp[-I (phi0 + 2 Pi (phi1/1000)  gyro ppm)] speci;
	Total[(Abs[speci] - Re[spec])^2] (*+ Total[Im@spec]^2*)
], RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsection:: *)
(*Padding and apodization*)


(* ::Subsubsection:: *)
(*PadFid*)


SyntaxInformation[PhaseCorrectSpectra] = {"ArgumentsPattern" -> {_}}

PadFid[fid_] := PadRight[fid, 2 Length[fid]]


(* ::Subsubsection:: *)
(*ApodizeFid*)


Options[ApodizeFid] = {ApodizationFunction -> "Ham"}

SyntaxInformation[PhaseCorrectSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeFid[fid_, OptionsPattern[]] := ApodizeFun[Length[fid], OptionValue[ApodizationFunction]] fid


(* ::Subsubsection:: *)
(*ApodizePadFid*)


Options[ApodizePadFid] = {ApodizationFunction -> "Ham"}

SyntaxInformation[ApodizePadFid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadFid[fid_, opts : OptionsPattern[]] := PadFid[ApodizeFid[fid, opts]]


(* ::Subsubsection:: *)
(*PadFid*)


SyntaxInformation[PadSpectra] = {"ArgumentsPattern" -> {_}}

PadSpectra[spec_] := ShiftedFourier[PadFid[ShiftedInverseFourier[spec]]]


(* ::Subsubsection:: *)
(*ApodizeFid*)


Options[ApodizeSpectra] = {ApodizationFunction -> "Ham"}

SyntaxInformation[ApodizeSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeSpectra[spec_, opts : OptionsPattern[]] := ShiftedFourier[ApodizeFid[ShiftedInverseFourier[spec], opts]]


(* ::Subsubsection:: *)
(*ApodizePadFid*)


Options[ApodizeSpectra] = {ApodizationFunction -> "Ham"}

SyntaxInformation[ApodizeSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadSpectra[spec_, opts : OptionsPattern[]] := ShiftedFourier[ApodizePadFid[ ShiftedInverseFourier[spec], opts]]


(* ::Subsubsection:: *)
(*ApodizeFun*)


ApodizeFun[length_, apM_ : "GausLaur"] := ApodizeFun[length, apM] = Block[{app},
	app = Switch[apM,
		"Han", 0.5 + 0.5 Table[Cos[x Pi/length], {x, 1, length}],
		"Ham", 0.54 + 0.46 Table[Cos[x Pi/length], {x, 1, length}],
		"Gaus", Table[N@PDF[HalfNormalDistribution[3/length], x], {x, 1, length}],
		"Laur", Table[Exp[-3 x/(length)], {x, 1, length}],
		"GausLaur", Table[Exp[2 x/(length)], {x, 1, length}] Table[N@PDF[HalfNormalDistribution[4/length], x], {x, 1, length}]
	];
	app = app/Max[app]
]


(* ::Subsection:: *)
(*Time and ppm*)


(* ::Subsubsection:: *)
(*ApodizeFun*)


SyntaxInformation[GetTimePpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetTimePpmRange[spec_, {dt_, field_, nuc_}] := GetTimePpmRange[spec, dt, field, nuc]

GetTimePpmRange[spec_, dt_, field_, nuc_] := GetTimePpmRange[spec, dt, GyromagneticRatio[nuc] field]

GetTimePpmRange[spec_, dt_, gyro_] := {GetTimeRange[spec, dt], GetPpmRange[spec, dt, gyro]}


(* ::Subsubsection:: *)
(*ApodizeFun*)


SyntaxInformation[GetPpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetPpmRange[spec_, {dt_, field_, nuc_}] := GetPpmRange[spec, dt, field, nuc]

GetPpmRange[spec_, dt_, field_, nuc_] := GetPpmRange[spec, dt, GetGyro[nuc, field]]

GetPpmRange[spec_, dt_, gyro_] := Block[{ppmBw},
	ppmBw = 1./(dt gyro);
	Reverse@Range[-ppmBw/2, ppmBw/2, ppmBw/(Length[spec] - 1)]
]


(* ::Subsubsection:: *)
(*GetTimeRange*)


SyntaxInformation[GetTimeRange] = {"ArgumentsPattern" -> {_, _}};

GetTimeRange[fid_, dt_] := N@Range[dt, Length[fid] dt, dt ]


(* ::Subsubsection:: *)
(*GetGyro*)


SyntaxInformation[GetGyro] = {"ArgumentsPattern" -> {_, _}};

GetGyro[nuc_, field_] := GyromagneticRatio[nuc] field


(* ::Subsubsection:: *)
(*PhaseShiftSpectra*)


SyntaxInformation[PhaseShiftSpectra] = {"ArgumentsPattern" -> {_, _, _., _.}}

PhaseShiftSpectra[spec_, phi0_] := PhaseShiftSpectraC[spec, 0., 0., phi0, 0.];

PhaseShiftSpectra[spec_, ppm_, gyro_, phi1_] := PhaseShiftSpectraC[spec, ppm, gyro, 0., phi1];

PhaseShiftSpectra[spec_, ppm_, gyro_, {phi0_, phi1_}] := PhaseShiftSpectraC[spec, ppm, gyro, phi0, phi1];

PhaseShiftSpectraC = Compile[{{spec, _Complex, 1}, {ppm, _Real, 1}, {gyro, _Real, 0}, {phi0, _Real, 0}, {phi1, _Real, 0}},
	Exp[-I (phi0 + 2 Pi (phi1/1000) gyro ppm)] spec,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


(* ::Subsubsection:: *)
(*TimeShiftFid*)


SyntaxInformation[TimeShiftFid] = {"ArgumentsPattern" -> {_, _, _., _.}}

TimeShiftFid[fid_, time_, gam_] := TimeShiftFidC[fid, time, 0., gam, 0., .5];

TimeShiftFid[fid_, time_, {gam_, f_}] := TimeShiftFidC[fid, time, 0., gam, 0., f];

TimeShiftFid[fid_, time_, gyro_, {gam_, eps_}] := TimeShiftFidC[fid, time, gyro, gam, eps, .5];

TimeShiftFid[fid_, time_, gyro_, {gam_, eps_, f_}] := TimeShiftFidC[fid, time, gyro, gam, eps, f];

TimeShiftFidC = Compile[{{fid, _Complex, 1}, {time, _Real, 1}, {gyro, _Real, 0}, {gam, _Real, 0}, {eps, _Real, 0}, {f, _Real, 0}},
	(f Exp[-gam time] + (1 - f) Exp[- gam^2 time^2]) Exp[2 Pi eps gyro I time] fid, 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]

(* ::Subsubsection:: *)
(*ChangeDwellTimeFid*)


SyntaxInformation[ChangeDwellTimeFid] = {"ArgumentsPattern" -> {_, _, _}}

ChangeDwellTimeFid[time_, dwOrig_, dwTar_] := Block[{NsampOrig, timeOrig, NsampTar, timeTar},
	(*get time of original signal*)
	NsampOrig = Length@time;
	timeOrig = dwOrig (Range[NsampOrig] - 1);
	(*define time of new signal*)
	NsampTar = Round[Max[timeOrig]/dwTar];
	timeTar = dwTar (Range[NsampTar] - 1);
	(*Interpolate the time to the new timescale*)
	Interpolation[Transpose[{timeOrig, time}], InterpolationOrder -> 1][timeTar]
]

(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
