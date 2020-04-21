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


(* ::Subsection::Closed:: *)
(*Functions*)


ReadjMRUI::usage = 
"ReadjMRUI[file] read a jMRUI spectrum file. 
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


GetSpectraBasisFunctions::usage = 
"GetSpectraBasisFunctions[{met1, ..., metn}] generates a list of spectra baisis functions with names met1 to metn. The names are strings and are the metabolites availible in GetSpinSystem.
GetSpectraBasisFunctions[{{props1}, ..., {propsn}}] generates a list of spectra baisis functions with properties prop1 to propn. The properties are those specified in MakeSpinSystem.
GetSpectraBasisFunctions[inp, split] generates a list of spectra basisfunctions. Each metabolite name present in the list split wil be split in individual spectra per peak."

FitSpectra::usage = 
"FitSpectra[specBasis, spec, {st,end}, dt, {lwvals,lwamsp}] Fits the basis spectra from GetSpectraBasisFunctions to the spec overt the ppm range {st, end} and dt the dweltime."


PlotFid::usage = 
"PlotFid[fid, dwell] 
PlotFid[time, fid]"

PlotSpectra::usage = 
"PlotSpectra[spec, {dwell, field, nuc}]
PlotSpectra[ppm, spec]" 

MakeSpectraGrid::usage = 
"MakeSpectraGrid[spectra]
MakeSpectraGrid[spectra, ppm]" 

FitSpectraResultTable::usage
"FitSpectraResultTable[parFit, parsF, names, ref, out] function not done"

CompareSpectraFitPlot::usage
"CompareSpectraFitPlot[ppmPl, specPlot, fitPlot] function not done"

MakeSpectraResultPlot::usage
"MakeSpectraResultPlot[ppmF_, specF_, {fit_, basisFit_}, names_, sc__ : 1, met__ : "ReIm"] function not done"


(* ::Subsection::Closed:: *)
(*Options*)


ApodizationFunction::usage = 
"ApodizationFunction is an options for ApodizeFid, ApodizeSpectra, ApodizePadFid, and ApodizePadSpectra. Values can be \"Han\", \"Ham\", \"Gaus\", \"Laur\", and \"GausLaur\"."

PaddingFactor::usage = 
"PaddingFactor is an option for PadFid, PadSpectra, ApodizePadFid, ApodizePadSpectra and FitSpectra. It Specifies with which factro to lengthen the fid."

BasisSequence::usage = 
"BasisSequence is an option for GetBasisFunction and specifies which sequence to use."

SpectraSamples::usage="SpectraSamples is an option for GetSpectraBasisFunctions."
SpectraBandwith::usage="SpectraBandwith is an option for GetSpectraBasisFunctions." 
SpectraLinewidth::usage="SpectraLinewidth is an option for GetSpectraBasisFunctions."
SpectraLinewidthShape::usage="SpectraLinewidthShape is an option for GetSpectraBasisFunctions."

SpectraNucleus::usage="SpectraNucleus is an option for GetSpectraBasisFunctions and FitSpectra."
SpectraPpmShift::usage="SpectraPpmShift is an option for GetSpectraBasisFunctions and FitSpectra."
SpectraFieldStrength::usage="SpectraFieldStrength is an option for GetSpectraBasisFunctions and FitSpectra."


SplineSpacingFactor::usage = 
"SplineSpacingFactor is an option for FitSpectra."
FineTuneFit::usage = 
"FineTuneFit is an option for FitSpectra."
InitializeFit::usage = 
"InitializeFit is an option for FitSpectra."
FitLineShape::usage = 
"FitLineShape is an option for FitSpectra."
SpectraOutputPlots::usage = 
"SpectraOutputPlots is an option for FitSpectra."

SpectraSpacing::usage = "SpectraSpacing is an option for PlotSpectra."

PlotScaling::usage = 
"PlotScaling is an option for MakeSpectraGrid"

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


(* ::Subsubsection::Closed:: *)
(*PadFid*)


Options[PadFid] = {PaddingFactor -> 2}

SyntaxInformation[PadFid] = {"ArgumentsPattern" -> {_}}

PadFid[fid_, OptionsPattern[]] := PadRight[fid, Round[OptionValue[PaddingFactor] Length[fid]]]


(* ::Subsubsection::Closed:: *)
(*ApodizeFid*)


Options[ApodizeFid] = {ApodizationFunction -> "Ham"}

SyntaxInformation[ApodizeFid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeFid[fid_, OptionsPattern[]] := ApodizeFun[Length[fid], OptionValue[ApodizationFunction]] fid


(* ::Subsubsection::Closed:: *)
(*ApodizePadFid*)


Options[ApodizePadFid] = {ApodizationFunction -> "Ham", PaddingFactor -> 2}

SyntaxInformation[ApodizePadFid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadFid[fid_, OptionsPattern[]] := PadFid[ApodizeFid[fid, ApodizationFunction->OptionValue[ApodizationFunction]],PaddingFactor->OptionValue[PaddingFactor]]


(* ::Subsubsection::Closed:: *)
(*PadFid*)


Options[PadSpectra] = {PaddingFactor -> 2}

SyntaxInformation[PadSpectra] = {"ArgumentsPattern" -> {_}}

PadSpectra[spec_, opts : OptionsPattern[]] := ShiftedFourier[PadFid[ShiftedInverseFourier[spec], opts]]


(* ::Subsubsection::Closed:: *)
(*ApodizeFid*)


Options[ApodizeSpectra] = {ApodizationFunction -> "Ham"}

SyntaxInformation[ApodizeSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizeSpectra[spec_, opts : OptionsPattern[]] := ShiftedFourier[ApodizeFid[ShiftedInverseFourier[spec], opts]]


(* ::Subsubsection::Closed:: *)
(*ApodizePadFid*)


Options[ApodizePadSpectra] = {ApodizationFunction -> "Ham", PaddingFactor -> 2}

SyntaxInformation[ApodizePadSpectra] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ApodizePadSpectra[spec_, opts : OptionsPattern[]] := ShiftedFourier[ApodizePadFid[ShiftedInverseFourier[spec], opts]]


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*GetTimePpmRange*)


SyntaxInformation[GetTimePpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetTimePpmRange[spec_, {dt_, field_, nuc_}] := GetTimePpmRange[spec, dt, field, nuc]

GetTimePpmRange[spec_, dt_, field_, nuc_] := GetTimePpmRange[spec, dt, GetGyro[nuc, field]]

GetTimePpmRange[spec_, dt_, gyro_] := {GetTimeRange[spec, dt], GetPpmRange[spec, dt, gyro]}


(* ::Subsubsection::Closed:: *)
(*GetPpmRange*)


SyntaxInformation[GetPpmRange] = {"ArgumentsPattern" -> {_, _, _., _.}};

GetPpmRange[spec_, {dt_, field_, nuc_}] := GetPpmRange[spec, dt, field, nuc]

GetPpmRange[spec_, dt_, field_, nuc_] := GetPpmRange[spec, dt, GetGyro[nuc, field]]

GetPpmRange[spec_, dt_, gyro_] := Block[{ppmBw},
	ppmBw = 1./(dt gyro);
	Reverse@Range[-ppmBw/2, ppmBw/2, ppmBw/(Length[spec] - 1)]
]


(* ::Subsubsection::Closed:: *)
(*GetTimeRange*)


SyntaxInformation[GetTimeRange] = {"ArgumentsPattern" -> {_, _}};

GetTimeRange[fid_, dt_] := N@Range[0, (Length[fid]-1) dt, dt ]


(* ::Subsubsection::Closed:: *)
(*GetGyro*)


SyntaxInformation[GetGyro] = {"ArgumentsPattern" -> {_, _}};

GetGyro[nuc_, field_] := GyromagneticRatio[nuc] field


(* ::Subsubsection::Closed:: *)
(*PhaseShiftSpectra*)


SyntaxInformation[PhaseShiftSpectra] = {"ArgumentsPattern" -> {_, _, _., _.}}

PhaseShiftSpectra[spec_, phi0_] := PhaseShiftSpectraC[spec, 0., 0., phi0, 0.];

PhaseShiftSpectra[spec_, ppm_, gyro_, phi1_] := PhaseShiftSpectraC[spec, ppm, gyro, 0., phi1];

PhaseShiftSpectra[spec_, ppm_, gyro_, {phi0_, phi1_}] := PhaseShiftSpectraC[spec, ppm, gyro, phi0, phi1];

PhaseShiftSpectraC = Compile[{{spec, _Complex, 1}, {ppm, _Real, 1}, {gyro, _Real, 0}, {phi0, _Real, 0}, {phi1, _Real, 0}},
	Exp[-I (phi0 + 2 Pi (phi1/1000) gyro ppm)] spec,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
]


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*GetBasisFunctions*)


Options[GetSpectraBasisFunctions] = {
	BasisSequence -> {"PulseAquire", 0},

	SpectraSamples -> 2046,
	SpectraBandwith -> 2000,
	
	SpectraLinewidth -> 5,
	SpectraLinewidthShape -> "LG",
	
	SpectraNucleus -> "1H",
	SpectraPpmShift -> 4.65,
	SpectraFieldStrength -> 3
};

SyntaxInformation[GetSpectraBasisFunctions] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

GetSpectraBasisFunctions[inp_, opts:OptionsPattern[]] := GetSpectraBasisFunctions[inp, {""}, opts]

GetSpectraBasisFunctions[inp_, split_, OptionsPattern[]] := Block[{
	cf, bw, nsamp, lw, lws, read, sys, din, struct, dr, te, svals, seq,
	sysAll, spinTab, field, nuc, nam, labs, names, times, fids, ppms, specs, douts},
	
	(*get the option values*)
	{seq, svals} = OptionValue[BasisSequence];
	
	nsamp = OptionValue[SpectraSamples];
	bw = OptionValue[SpectraBandwith];
	
	cf = OptionValue[SpectraPpmShift];
	lw = OptionValue[SpectraLinewidth];
	lws = OptionValue[SpectraLinewidthShape] /. {"Lorentzian"->"L", "Gaussian"->"G", "Voigt"->"LG", _->"LG"};
	
	field = OptionValue[SpectraFieldStrength];
	nuc = OptionValue[SpectraNucleus];
	
	(*get or make the spin systems*)
	sysAll = If[StringQ[#], GetSpinSystem, MakeSpinSystem][#, CenterFrequency -> cf] & /@ inp;
	spinTab = SysTable[sysAll];
	
	(*loop over systems to perform sysmulation of basis functions*)
	{names, times, fids, ppms, specs} = Transpose[Flatten[
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
				"PulseAquire",
				te = svals;
				dr = SequencePulseAcquire[din, struct, te];
			];
			
			(*simulate readout*)
			{times, fids, ppms, specs, douts} = SimReadout[dr, struct, ReadoutSamples -> nsamp, ReadoutBandwith -> bw, CenterFrequency -> cf, 
				Linewidth -> lw, LinewidthShape -> lws, ReadoutOutput -> read];
				
			(*check if multiple fids for metabolite and act acordingly*)
			If[VectorQ[fids],
				(*one fid*)
				fids = (Re[fids] + Im[fids] I);
				specs = ShiftedFourier[Re[fids] + Im[fids] I] Exp[-Pi I];
				{{nam, times, fids, Reverse[ppms], specs}}
				,
				(*mulitple fids*)
				MapIndexed[(
					fids = (Re[#1] + Im[#1] I);
					specs = ShiftedFourier[Re[fids] + Im[fids] I] Exp[-Pi I];
					{nam <> "-" <> labs[[#2]], times, fids, Reverse[ppms], specs}
				) &, fids]
			]
		
		(*close mapping function*)
		) &, sysAll]
	, 1]];
	
	{names, times, fids, ppms, specs, spinTab}
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
	SpectraOutputPlots->False
};

SyntaxInformation[FitSpectra] = {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}}

FitSpectra[specBasisIn_, specIn_, {st_,end_}, dtime_, lwvals_?VectorQ, opts : OptionsPattern[]]:=FitSpectra[specBasisIn,specIn,{st,end},dtime,{lwvals, 0lwvals + 1.},opts]

FitSpectra[specBasisIn_, specIn_, {st_,end_}, dtime_, {lwvals_?VectorQ, lwamsp_?VectorQ}, OptionsPattern[]]:=Block[{
	ttotal,log,pad,spfac,field,nuc,shift,plots,init,scale,nbas,len,
	timeBasis,specFull,timeFull,ppmFull,nsamp,gyro,indSt,indEnd,
	gami,epsi,phi0i,phi1i,linei,phii,plLine,plShift,
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
		pad=OptionValue[PaddingFactor];
		spfac=OptionValue[SplineSpacingFactor];
		field= OptionValue[SpectraFieldStrength];
		nuc =OptionValue[SpectraNucleus];
		shift = OptionValue[SpectraPpmShift];
		plots=OptionValue[SpectraOutputPlots];
		init=OptionValue[InitializeFit];
		
		(*set general parameters*)
		scale=1000/Max[Abs[specIn]];
		nbas=Length[specBasisIn];
		
		
		(*-------------------------------------------------------------------*)
		(*create the basis functions fid and pad*)
		timeBasis=CashBasisTime[specBasisIn,pad];
		
		(*pad and normalize the spectra spectra*)
		specFull=scale ShiftedFourier[ApodizePadFid[ShiftedInverseFourier[specIn],PaddingFactor->pad]];
		
		(*get the time and ppm axes*)
		{timeFull,ppmFull}=GetTimePpmRange[specFull,dtime,field,nuc];
		gyro=GetGyro[nuc,field];
		nsamp=Length[specFull];
		
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
			{gami,epsi,{phi0i,phi1i},linei}=init;
			phii={phi0i,phi1i};
			,
			(*find initial linewidth and spec shift*)
			{gami,epsi,plLine}=EstimateLineWidth[{ppmFull,specFull},{lwvals,lwamsp},gyro,{st,end},plots];
			(*find initial phase estimate*)
			{{phi0i,phi1i},plShift}=EstimatePhaseShift[{ppmFull,specFull},{timeFull,timeBasis},{gami,epsi},gyro,{indSt,indEnd},plots];
			(*define the initial line shape*)
			linei=0.5;
			
			(*logging of parameters*)
			AppendTo[log,Style["Estimating linewidth",Bold]];
			AppendTo[log,"    - spectral linewidth intit:   "<>ToString[Round[gami /gyro,.0001]]<>" ppm"];
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
			Clear[phi0f,phi1f,gamf,epsf,linef];
			var={{gamf,gami},{epsf,epsi},{linef,linei},{phi0f,phi0i},{phi1f,phi1i}};
			init={gami,epsi,1,phi1i};
			
			(*perform the fit*)
			{tfit1,fit1}=AbsoluteTiming[FindMinimum[
				FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},{gamf,epsf,{phi0f,phi1f}, linef},init,Output->"Error"], 
				var,MaxIterations->50,Method->"QuasiNewton"][[2]]];
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
		splineSpace=spfac Mean[Flatten[{gami}]]/gyro;
		cpn = Clip[Round[Subtract@@Reverse[MinMax[ppmFull]]/splineSpace],{4,Round[nsamp/10]}];
		
		(*make the output*)
		{fit,sigi}=FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},sol,Output->"Fit"];
		
		
		(*-------------------------------------------------------------------*)
		(*perform the second run minimization*)
		tfit2=0;
		If[OptionValue[FineTuneFit],
			(*prepare the fit parameters*)
			Clear[phi0f,phi1f,phi,gam,eps,line,var,init];
			gamf=Table[Unique[gam],{i,1,nbas}];
			epsf=Table[Unique[eps],{i,1,nbas}];
			linef=If[OptionValue[FitLineShape],Table[Unique[line],{i,1,nbas}],Unique[line]];
			
			(*define the fit varables*)
			var=Join[MakeVars[gamf,gami,0],MakeVars[epsf,epsi,0],MakeVars[{phi0f,phi1f},phii,0],MakeVars[linef,linei,0]];
			
			(*get the std of initial fit*)
			init={gami,epsi,sigi,phii[[2]]};
			
			(*perform the minimization*)
			{tfit2,fit2}=AbsoluteTiming[FindMinimum[
				FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},{gamf,epsf,{phi0f,phi1f},linef},init, Output->"Error"],
				var, MaxIterations->500, Method->"QuasiNewton"][[2]]];
			
			(*get the solution and output, wrap phi between -pi and pi*)
			sol={gami, epsi, phii, linei}={Clip[gamf,{1,500}],epsf,{2ArcTan[Tan[phi0f/2]],phi1f},linef}/.fit2;
			
			(*recalculate the spline spacings*)
			splineSpace=spfac Mean[Flatten[{gami}]]/gyro;
			cpn = Clip[Round[Subtract@@Reverse[MinMax[ppmFull]]/splineSpace],{4,Round[nsamp/10]}];
			
			(*generate the output*)
			{fit,sigi}=FitSpectraError[{ppmFull,specFull},{timeFull,timeBasis},{indSt,indEnd},{cpn,gyro},sol,Output->"Fit"];
			
			(*logging of parameters*)
			AppendTo[log,Style["Performing spectra: Second run",Bold]];
			AppendTo[log,"    - line shape:                 "<>ToString[Round[Mean[Flatten[{linei}]],.0001]]];
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
	
	timeBasisIn=(ShiftedInverseFourier[#]&/@specBasisIn);
	basis=BasisSpectraApply[{ppm,time,timeBasisIn},sol,gyro];
	fit=fit/scale;
	specFit=fit.basis;
	
	(*fit a spline through the residuals*)
	spline=BSplineCurveFit[specIn-specFit,SplineKnotsNumber-> cpn,SplineRegularization->0,SplineDegree-> 2];
	
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
	(*give the output,  scale data to original values*)
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
		Transpose[{par,RandomReal[{0.99,1.01},Length[par]]ConstantArray[val,Length[par]]}]
	]
]

CashBasisTime[specBasisIn_,pad_]:=CashBasisTime[specBasisIn,pad]=ApodizePadFid[ShiftedInverseFourier[#],PaddingFactor->pad]&/@specBasisIn


(* ::Subsubsection::Closed:: *)
(*Fit Basis spectra*)


Options[FitSpectraError] = {Output -> "Error"};

FitSpectraError[{ppmFull_, spec_}, {timeFull_, timeBasis_}, {indSt_, indEnd_}, {cpn_, gyro_},
	{gam_ /; AllTrue[gam, NumericQ], eps_ /; AllTrue[eps, NumericQ], phi_ /; AllTrue[phi, NumericQ], f_ /; AllTrue[f, NumericQ]},
	init___ : 0, OptionsPattern[]] := Block[{
	specF, specBasisF, fit, ferr, gerr, error, reerr, imerr, gami, epsi, sigi, phii, specBasis, specFit, spline
	},
	(*give ouput either error or fit results*)
	Switch[OptionValue[Output],
		(*output error of fit calculated on selected ppm range*)
		"Error",
		(*apply phase to target instead of basis functions (faster)*)
		specF = PhaseShiftSpectra[spec[[indSt ;; indEnd]], ppmFull[[indSt ;; indEnd]], gyro, -phi];
		(*generate basis spectra from time domain by applying gam, eps and lineshape*)
		specBasisF = BasisSpectraApply[{timeFull, timeBasis}, {gam, eps, f}, gyro, {indSt, indEnd}];
		
		(*perform Fit of basis spectra*)
		fit = Quiet@NNLeastSquares[Transpose[Re[specBasisF]], Re[specF]];
		(*constrain f between 0 and 1 using power function*)
		ferr = If[NumberQ[f],10( ((f - 0.5)/.6)^100), Total[10(((# - 0.6)/.6)^100) & /@ f]];
		(*ferr = If[NumberQ[f], (2 (Ramp[f - 0.4] + Ramp[-(f - 0.6)]))^8, Total[(2 (Ramp[# - 0.4] + Ramp[-(# - 0.6)]))^8 &/@f]];*)
		(*ferr=0;*)
		(*constrain gam to be positive*)
		gerr = If[NumberQ[gam], (UnitStep[-(gam - 2)] (gam - 2))^4, Total[(UnitStep[-(# - 2)] (# - 2))^4 & /@ gam]];
		(*define errors*)
		error = specF - fit.specBasisF;
		
		(*rel and Im error normalized for number of points*)
		reerr = Total[Re[error]^2]/Length[error];
		imerr = Total[Im[error]^2]/Length[error];
		
		If[init === 0,
			(*no initial values only minimize RMSE*)
			reerr + imerr + ferr + gerr,
			(*get initial values to constrain finetune fit*)
			{gami, epsi, sigi, phii} = init;
			(*calculate error, normalize for sigma contrain gam en eps*)
			reerr/sigi + imerr/sigi + Total[((gam - gami)/40)^4] + Total[((eps - epsi)/1)^4] + ferr + gerr
		],
		
		(*ouput the fit results calculated on full ppm range*)
		"Fit",
		(*generate basis spectra from time domain by applying gam, eps and lineshape*)
		specBasis = BasisSpectraApply[{ppmFull, timeFull, timeBasis}, {gam, eps, phi, f}, gyro];
		(*perform Fit of basis spectra*)
		fit = Quiet@NNLeastSquares[Transpose[Re[specBasis]], Re[spec]];
		specFit = fit.specBasis;
		(*fit a spline through the residuals*)
		spline = BSplineCurveFit[spec - specFit, SplineKnotsNumber -> cpn, SplineRegularization -> 0, SplineDegree -> 2];
		(*recalculate the error*)
		error = spec - specFit - spline;
		{fit, StandardDeviation[error]^2}
	]
]


(* ::Subsubsection::Closed:: *)
(*BasisSpectraApply*)


BasisSpectraApply[{timeFull_,timeBasis_},{gam_,eps_,f_},gyro_]:=BasisSpectraApply[{0,timeFull,timeBasis},{gam,eps,0,f},gyro,{1,Length[timeFull]}]

BasisSpectraApply[{timeFull_,timeBasis_},{gam_,eps_,f_},gyro_,{st_,end_}]:=BasisSpectraApply[{0,timeFull,timeBasis},{gam,eps,0,f},gyro,{st,end}]

BasisSpectraApply[{ppmFull_,timeFull_,timeBasis_},{gam_,eps_,phi_,f_},gyro_]:=BasisSpectraApply[{ppmFull,timeFull,timeBasis},{gam,eps,phi,f},gyro,{1,Length[timeFull]}]

BasisSpectraApply[{ppmFull_,timeFull_,timeBasis_},{gam_,eps_,phi_,f_},gyro_,{st_,end_}]:=Block[{specBasis},
	(*generate basis spectra from time domain by applying gam, eps and lineshape*)
	specBasis=If[NumberQ[gam]&&NumberQ[eps]&&NumberQ[f],
		(*global gamma and epsilon*)
		Map[ShiftedFourier[TimeShiftFid[#1,timeFull,gyro,{gam,eps,f}]][[st;;end]]&,timeBasis],
		If[VectorQ[gam]&&VectorQ[eps]&&NumberQ[f],
			(*basis function specific gamma and epsilon*)
			MapThread[ShiftedFourier[TimeShiftFid[#1,timeFull,gyro,{#2,#3,f}]][[st;;end]]&,{timeBasis,gam,eps}],
			(*basis function specific gamma, epsilon and shape*)
			MapThread[ShiftedFourier[TimeShiftFid[#1,timeFull,gyro,{#2,#3,#4}]][[st;;end]]&,{timeBasis,gam,eps,f}]
		]
	];
	
	(*apply phase to the basis spectra*)
	If[phi===0,specBasis,PhaseShiftSpectra[#,ppmFull,gyro,phi]&/@specBasis]
]


(* ::Subsubsection::Closed:: *)
(*Estimate Line width*)


(*Function to estimate linewidth*)
EstimateLineWidth[{ppm_,spec_},{peaks_,amps_},gyro_,ran_,plot_:True]:=Block[{
	dppm,deltaf,corrf,wgth, max,corrf2,pts,pos,ppmC,maxf,block,line,sol,x,ppmCf,pl,lw,sft},
	
	(*define delta ppm and the delta function for correlation*)
	dppm=(ppm[[1]]-ppm[[2]]);
	deltaf=0 ppm;
	deltaf[[Flatten[Position[ppm,First@Nearest[ppm,#]]&/@peaks]]]=amps/Max[amps];
	
	(*perfomr correlation of spectra with delta function*)
	corrf=ListCorrelate[deltaf,Abs[spec],Round[Length[deltaf]/2],0];
	corrf=DevideNoZero[Length[corrf]corrf,Max[corrf]];
	
	(*Find max correlation and position*)
	wgth=Table[N@PDF[NormalDistribution[0,Max[ppm]/4],x],{x,ppm}];
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
		FlipView[{
			Show[
				PlotSpectra[ppm,Max[Abs[spec]]deltaf,Method->"Abs",GridLineSpacing->5,PlotRange->{ran,Full},PlotColor->Red],
				PlotSpectra[ppm,spec,Method->"Abs",GridLineSpacing->10,PlotRange->{ran,{0,Max[Abs[spec]]}}]
			,ImageSize->1000],
			Show[
				PlotSpectra[ppmC,corrf,GridLineSpacing->5,PlotRange->{ran,Full}],
				ListLinePlot[{Transpose[{ppmCf[sol],{maxf,maxf}}],{{ppmC[[pos]],0},{ppmC[[pos]],max}}},PlotStyle->Directive[{Thick,Red}],ScalingFunctions->{"Reverse",Automatic}]
			,ImageSize->1000]
		}]
		,
		Null
	];
	
	(*calculate the estimated lw and shift*)
	{lw,sft,pl}
]


(* ::Subsubsection::Closed:: *)
(*EstimatePhaseShift*)


(*function to estimate phase form abs fitting of baiss spectra*)
EstimatePhaseShift[{ppm_,spec_},{time_,fids_},{gam_,eps_},gyro_,{st_,en_},plot_:True]:=Block[{
	phi1,sol1,phi2,sol2,specsC,fit,phi0f,phi1f,phi,specs,ran,pl,lim,specf,ppmf
	},
	
	specf=spec[[st;;en]];
	ppmf=ppm[[st;;en]];
	lim=.1;
	
	(*convert basis fids in spectra*)
	specsC=Transpose[ShiftedFourier[TimeShiftFid[#,time,gyro,{gam,eps,.5}]][[st;;en]]&/@fids];
	(*Fit absolute basis spectra to absolute spectrum*)
	fit=specsC.(NNLeastSquares[Abs[specsC],Abs[specf]]);
	
	(*minimize error with the target spectra*)
	sol1=Quiet@NMinimize[{PhaseError[ppmf,fit,specf,{phi0f,0 },gyro],-Pi<phi0f<Pi},{phi0f},MaxIterations->25][[2]];
	phi1={phi0f, 0}/.sol1;
	
	(*apply the zeroth order phase to the basis spectra*)
	specsC=Transpose[PhaseShiftSpectra[#,ppmf,gyro,phi1]&/@Transpose[specsC]];
	(*calculate the fit based on the imaginary part of the spectra*)
	fit=specsC.(NNLeastSquares[Re@specsC,Re@specf]);
	
	(*minimize error with the target spectra*)
	sol2=Quiet@NMinimize[{PhaseError[ppmf,fit,specf,{phi0f,phi1f },gyro],-Pi<phi0f<Pi,-lim<phi1f<lim},{phi0f,phi1f},MaxIterations->25][[2]];
	phi2={phi0f, phi1f}/.sol2;
	phi=phi1+phi2;
	
	(*debugging plots*)
	pl = If[plot,
		specs=PhaseShiftSpectra[fit,ppmf,gyro,phi2];
		ran={-1,1}Max[Abs[{specf,specs}]];
		FlipView[{
			Show[PlotSpectra[ppmf,specs,GridLineSpacing->5,PlotRange->{Full,ran}],ImageSize->1000],
			Show[PlotSpectra[ppmf,specf,GridLineSpacing->5,PlotRange->{Full,ran}],ImageSize->1000]
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
	PlotLabel -> None
};

SyntaxInformation[PlotSpectra] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

PlotSpectra[spec_, {dwell_?NumberQ, field_?NumberQ, nuc_?StringQ}, opts : OptionsPattern[]] := PlotSpectra[GetPpmRange[spec, dwell, field, nuc], spec, opts]

PlotSpectra[ppm_?VectorQ, spec_, OptionsPattern[]] := Block[{
	fun, plot, plot2, grid, gridS, or, rr, col, space, cols, cols2, pl1, pl2, lables
	},
	
	(*get gridline options*)
	gridS = OptionValue[GridLineSpacing];
	grid = Sort@DeleteDuplicates@Join[
		If[gridS === 0, {}, Range[0, Round[Max[ppm]], gridS]],
		If[gridS === 0, {}, -Range[0, -Round[Min[ppm]], gridS]],
		OptionValue[GridLines]
	];
	
	(*get the plot range*)
	rr = OptionValue[PlotRange] /. Automatic -> Full;
	Switch[rr,
		{_, Full}, rr[[2]] = {Min[{Re[spec], Im[spec]}], Max[Abs[spec]]},
		Full, rr = {Full, {Min[{Re[spec], Im[spec]}], Max[Abs[spec]]}}
	];
	
	(*if only one spectra plot normal else plot expanded list of spectra*)
	If[VectorQ[spec],
		(*plot single spectra*)
		(*get the plot functions*)
		fun = Switch[OptionValue[Method], "Abs", {Abs}, "Re", {Re}, "Im", {Im}, "ReIm", {Im, Re}, "All", {Im, Re, Abs}];
		(*plot single spectra*)
		plot = Transpose[{ppm, #}] & /@ (#@spec & /@ fun);
		(*get the plot color*)
		col = If[OptionValue[PlotColor] === Automatic, ({Gray, {Red, Thin}, {Thick, Black}}[[-Length[fun] ;;]]), OptionValue[PlotColor]];
		
		(*Make the plot*)
		ListLinePlot[plot, PlotStyle -> col, PlotRange -> rr, GridLines -> {grid, None}, AspectRatio -> OptionValue[AspectRatio],
			ImageSize -> OptionValue[ImageSize], PlotLabel -> OptionValue[PlotLabel], ScalingFunctions -> {"Reverse", Automatic},
			Frame -> {{False, False}, {True, False}}, FrameStyle -> Directive[{Thick, Black}], FrameLabel -> {"PPM", None},
			LabelStyle -> {Bold, 14, Black}
		]
		
		,
		(*plot List of spectra*)
		(*get the plot functions*)
		fun = Switch[OptionValue[Method], "Abs", Abs, "Re", Re, "ReIm", Re, "Im", Im, _, Return[]];
		(*space the spectra over the y axes*)
		space = Reverse@Range[0, Length[spec]] Max[Abs[spec]] OptionValue[SpectraSpacing];
		
		(*plot the spectra*)
		plot = Transpose[{ppm, #}] & /@ (fun[Append[spec, Total[spec]]] + space);
		(*correct the plot range*)
		If[rr[[2]] =!= Full, rr[[2, 2]] = 1.1 Max[plot[[All, All, 2]]]];
		If[rr[[2]] =!= Full, rr[[2, 1]] = Min[plot[[All, All, 2]]] - .1 Max[plot[[All, All, 2]]]];
		(*get the plot colors*)
		cols = Thread[{Append[ConstantArray[Black, Length[plot] - 1], Red], Thick}];
		lables = OptionValue[PlotLabels];
		
		(*make the plot*)
		pl1 = ListLinePlot[plot, Frame -> {{False, False}, {True, False}}, FrameStyle -> Thickness[.003], FrameTicksStyle -> Thickness[.003],
			PlotRange -> rr, PlotStyle -> cols, ScalingFunctions -> {"Reverse", Automatic}, 
			PlotLabels ->If[(OptionValue[Method] === "ReIm") || (lables === None), None, (Style[#, Black, Bold, 14] & /@ Append[lables, "All"])],
			GridLines -> {grid, None}, PlotRange -> rr, AspectRatio -> .5, ImageSize -> 1000, FrameLabel -> {"PPM", None}, LabelStyle -> {Bold, 14, Black}
		];
		
		(*make Im plot if needed*)
		If[OptionValue[Method] === "ReIm",
			plot2 = Transpose[{ppm, #}] & /@ (Im[Append[spec, Total[spec]]] + space);
			cols2 = {Append[ConstantArray[Gray, Length[plot] - 1], Gray]};
			
			pl2 = ListLinePlot[plot2, Frame -> {{False, False}, {True, False}}, FrameStyle -> Thickness[.003], FrameTicksStyle -> Thickness[.003], 
				PlotRange -> rr, PlotStyle -> cols2, ScalingFunctions -> {"Reverse", Automatic},
				PlotLabels -> (Style[#, Black, Bold, 14] & /@ Append[OptionValue[PlotLabels], "All"]),
				GridLines -> {grid, None}, PlotRange -> rr, AspectRatio -> .5,
				ImageSize -> 1000, FrameLabel -> {"PPM", None},
				LabelStyle -> {Bold, 14, Black}
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
	
	fun = Switch[OptionValue[Method], "Abs", {Abs}, "ReIm", {Im, Re}, "All", {Im, Re, Abs}];
	
	plot = Transpose[{time, #}] & /@ (#@fid & /@ fun);
	
	rr = OptionValue[PlotRange];
	Switch[rr,
		{_, Full}, rr[[2]] = {-Max[Abs[fid]], Max[Abs[fid]]},
		Full, rr = {Full, {-Max[Abs[fid]], Max[Abs[fid]]}}
	];
	col = If[OptionValue[PlotColor] === Automatic, ({Gray, {Red, Thin}, {Thick, Black}}[[-Length[fun] ;;]]), OptionValue[PlotColor]];
	
	ListLinePlot[plot, PlotStyle -> col, PlotRange -> rr, GridLines -> {grid, {0.}}, AspectRatio -> OptionValue[AspectRatio],
		ImageSize -> OptionValue[ImageSize], PlotLabel -> OptionValue[PlotLabel], Frame -> {{False, False}, {True, False}},
		FrameStyle -> Directive[{Thick, Black}], FrameLabel -> {"time [s]", None}, LabelStyle -> {Bold, 14, Black}
	]
]



(* ::Subsection::Closed:: *)
(*MakeSpectraGrid*)


Options[MakeSpectraGrid] = {
	Method -> "Abs", 
	PlotScaling -> "Max",
	PlotRange -> Full,
	ImageSize->50
};

SyntaxInformation[MakeSpectraGrid] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}}

MakeSpectraGrid[spectra_, opts : OptionsPattern[]] := MakeSpectraGrid[spectra, 0, opts]

MakeSpectraGrid[spectra_, ppm_, OptionsPattern[]] := Block[{
	met, sc, max, spectraP, maxs, xdat, mid, leg, ran, grid, xran, speci, maxi, plots
	},
	
	(*get options*)
	met = OptionValue[Method];
	sc = OptionValue[PlotScaling];
	xran = OptionValue[PlotRange];
	
	spectraP = Switch[met, "Abs", Abs, "Im", Im, "Re", Re][spectra];
	
	(*prepare values*)
	max = Max[Abs[spectraP]];
	spectraP = spectraP/max;(*scale all between 0 and 1*)
	maxs = Map[Max, spectraP, {-2}];
	
	xdat = If[ppm === 0,
		mid = Dimensions[spectra][[-1]]/2;
		Range[Dimensions[spectra][[-1]]]-mid
		,
		-ppm
	];
	
	(*make legend bar*)
	leg = BarLegend[{{"DarkRainbow", "Reverse"}, {0, 100}}, LegendLayout -> "Row", LegendMarkerSize -> 400, LabelStyle -> Directive[Large, Black, Bold]];
	
	ran = {If[met === "Abs", -0.1, -1], 1} If[NumberQ[sc], sc, Switch[sc, "Max", 1, "Full", 1, "Half", 0.5, "Third", 0.3]];
	If[sc === "Max", spectraP = spectraP/maxs];
	
	grid = If[ppm === 0,
		grid = {{{0, Black}}, {{0, Black}}},
		xran = -If[xran === Full, MinMax[xdat], xran];
		grid = {#, Lighter@Gray} & /@ Range[Round[xran[[1]], 10], Round[xran[[2]], 10], 10];
		{Join[grid, {{0, Black}}], {{0, Black}}}
	];
	
	plots = Map[(
		{speci, maxi} = #;
		speci = Transpose[{xdat, speci}];
		Graphics[{Directive[{(*Thick,*)ColorData[{"DarkRainbow", "Reverse"}][1.25 Sqrt[maxi]]}], Line[speci]}, AspectRatio -> 0.7, ImageSize -> OptionValue[ImageSize],
			PlotRange -> {xran, ran}, GridLines -> grid, PlotRegion -> {{0.05, 0.95}, {.1, 0.9}}
			(*,Frame\[Rule]{{False, False},{True,False}},FrameStyle\[Rule]{Thick,Black}*)
		]
	) &, TransData[{spectraP, maxs}, "l"], {3}];
	
	Column[{Grid[#, Spacings -> 0], leg}, Alignment -> Center] & /@ plots
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
		Flatten@Thread[{Style[#, Bold] & /@ {"\!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\) [deg]", "\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\) [ms]"},
			Round[{parsF[[3, 1]]/Degree, parsF[[3, 2]]}, .001]}
		]
	};
	sc = If[ref === "", 1, Clip[par[[Position[names, ref][[1, 1]]]], {Max[DeleteCases[par, 0.]],Infinity}]];
	amp = If[sc === 1, par, 100 par/sc];
		
	{lw, ls, shift} = parsF[[{1, 4, 2}]];
	If[NumberQ[ls], ls = ConstantArray[ls, Length[lw]]];
	If[Length[lw] == 0, {lw, ls, shift} = ConstantArray[#, Length[names]] & /@ {lw, ls, shift}];
	
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
(*FitSpectraResultTable*)


SyntaxInformation[CompareSpectraFitPlot] = {"ArgumentsPattern" -> {_, _, _, _.}}

CompareSpectraFitPlot[ppmPl_, specPlot_, fitPlot_, ranPpm_:Full] := Block[{ran, sp},
  ran = {-1, 1} Max[Abs[specPlot], Abs[fitPlot]];
  sp = 2;
  Column[{FlipView[{
      Column[{
        PlotSpectra[ppmPl, specPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "Abs"],
        PlotSpectra[ppmPl, specPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "ReIm"]
        }],
      Column[{
        PlotSpectra[ppmPl, fitPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "Abs"],
        PlotSpectra[ppmPl, fitPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "ReIm"]
        }]
      }]
    ,
    FlipView[{
    	PlotSpectra[ppmPl, fitPlot - specPlot, GridLineSpacing -> sp, PlotRange -> {ranPpm, ran}, Method -> "ReIm"],
      PlotSpectra[ppmPl, fitPlot - specPlot, GridLineSpacing -> 2, PlotRange -> {ranPpm, Full}, Method -> "ReIm"]
      }]
    }]
  ]


(* ::Subsection::Closed:: *)
(*FitSpectraResultTable*)


SyntaxInformation[MakeSpectraResultPlot] = {"ArgumentsPattern" -> {_, _, _, _, _.}}

MakeSpectraResultPlot[ppmF_, specF_, {fit_, basisFit_}, names_, ppmran_] := Block[{
	sp, specFit, resTotPl, errPl, fitPl, resPl, outPl, pran, pmax,
	lab1, lab2, resfitRI, resfit, resBasPl, met},
	
	sp = 2;
	met = "ReIm";
	specFit = fit.basisFit;
	pmax = Max[Abs[specFit], Abs[specF]];
	pran = {-pmax, pmax};
	
	resTotPl = Column[{
		FlipView[errPl = {
			PlotSpectra[ppmF, specF - specFit, Method -> met, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
			PlotSpectra[ppmF, specF - specFit, Method -> met, PlotRange -> {ppmran, Full}, GridLineSpacing -> sp]
		}],
		FlipView[fitPl = {
			Show[
				PlotSpectra[ppmF, specF, Method -> met /. "ReIm" -> "Re", PlotColor -> Red, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, basisFit[[1]], Method -> met /. "ReIm" -> "Re", PlotColor -> Green, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, specFit, Method -> met /. "ReIm" -> "Re", PlotColor -> Black, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp]
			],
			Show[
				PlotSpectra[ppmF, specF, Method -> met /. "ReIm" -> "Im", PlotColor -> Red, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, basisFit[[1]], Method -> met /. "ReIm" -> "Im", PlotColor -> Green, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
				PlotSpectra[ppmF, specFit, Method -> met /. "ReIm" -> "Im", PlotColor -> Black, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp]
			]
		}]
		,
		FlipView[resPl = {
			PlotSpectra[ppmF, specF, Method -> met, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp],
			PlotSpectra[ppmF, specFit, Method -> met, PlotRange -> {ppmran, pran}, GridLineSpacing -> sp]
		}]
	}, Alignment -> Center];
	
	lab1 = Style[#, Bold, Black, Large] & /@ {"Error scaled to signal", "Error scaled to Max"};
	lab2 = Style[#, Bold, Black, Large] & /@ {"Fit and signal Re", "Fit and signal Im"};
	outPl = Column[Flatten[{Thread[{lab1, errPl}], Thread[{lab2, fitPl}]}], Alignment -> Center];
	
	resBasPl = FlipView[{
		resfitRI = PlotSpectra[ppmF, fit basisFit, Method -> "ReIm", PlotColor -> Red, SpectraSpacing -> 0.2, GridLines -> {}, GridLineSpacing -> sp, PlotLabels -> Prepend[names, "spline"],PlotRange->{ppmran,Full}],
		resfit = PlotSpectra[ppmF, fit basisFit, Method -> "Abs", PlotColor -> Red, SpectraSpacing -> 0.2, GridLines -> {}, GridLineSpacing -> sp, PlotLabels -> Prepend[names, "spline"],PlotRange->{ppmran,Full}]
	}];
	
	{resTotPl, resBasPl, {errPl, fitPl, resPl, outPl, resfitRI, resfit}}
]

(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
