(* ::Package:: *)

(* ::Title:: *)
(*QMRITools DixonTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`DixonTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`DixonTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DixonToPercent::usage = 
"DixonToPercent[water, fat] converts the dixon water and fat data to percent maps.

Output is {waterFraction, fatFraction}.
The values of water and fat are arbitraty units and the ouput fractions are between 0 and 1."

DixonReconstruct::usage = 
"DixonReconstruct[real, imag, echo] reconstruxt Dixon data with initital guess b0 = 0 and T2star = 0.
DixonReconstruct[real, imag, echo, b0] reconstructs Dixon data with intitial guess T2star = 0.
DixonReconstruct[real, imag, echo, b0, t2] reconstructs Dixon data.

real is the real data in radials.
imag is the imaginary data in radians.
B0 can be estimated from two phase images using Unwrap.
T2 can be estimated from multiple echos using T2fit.

Output is {{watF,fatF},{watSig,fatSig},{inphase,outphase},{B0,T2star},itterations}.

The fractions are between 0 and 1, the B0 field map is in Hz and the T2start map is in ms.

DixonReconstruct[] is based on DOI: 10.1002/mrm.20624 and 10.1002/mrm.21737."

SimulateDixonSignal::usage = 
"SimulateDixonSignal[echo, fr, B0, T2] simulates an Dixon gradient echo sequence with echotimes.
Echotimes echo in ms, fat fraction fr between 0 and 1, field of resonance B0 in Hz and relaxation T2 in ms."


Unwrap::usage = 
"Unwrap[data] unwraps the given dataset.

The data should be between -Pi and Pi. 

Unwrap[] is based on DOI: 10.1364/AO.46.006623 and 10.1364/AO.41.007437."

UnwrapSplit::usage = 
"UnwrapSplit[phase, data] unwarps the give phase dataset but splits the data into left and right using SplitData based in the data and performs the unwrapping seperately.

The data should be between -Pi and Pi.

UnwrapSplit[] is based on DOI: 10.1364/AO.46.006623 and 10.1364/AO.41.007437."


(* ::Subsection::Closed:: *)
(*Options*)


DixonPrecessions::usage = 
"DixonPrecessions is an options for DixonReconstruct. Defines the rotation of the signal {-1,1} default is -1."

DixonFieldStrength::usage = 
"DixonFieldStrength is an options for DixonReconstruct. Defines the fieldstrengths in Tesla on which the data was acquired."

DixonFrequencies::usage = 
"DixonFrequencies is an options for DixonReconstruct. Defines the frequencies in ppm of the fat peaks being used."

DixonAmplitudes::usage = 
"DixonAmplitudes is an options for DixonReconstruct. Defines the relative amplitudes of the fat peaks being used."

DixonNucleus::usage = 
"DixonNucleus is an option for DixonReconstruct. Defines the nucleus for which the reconstruction is performed."

DixonTollerance::usage = 
"DixonTollerance is an options for DixonReconstruct. Defines at which change per itteration of b0 and R2star the ittarative methods stops. Default value is 0.1."

DixonMaskThreshhold::usage = 
"DixonMaskThreshhold is an options for DixonReconstruct. Defines at which threshhold the dixon reconstruction considers a voxel to be background noise. Defualt values is 0.05."

DixonFilterInput::usage = 
"DixonFilterInput is an options for DixonReconstruct. If True the input b0 and T2star values are smoothed using a gaussian kernel."

DixonFilterOutput::usage = 
"DixonFilterOutput is an options for DixonReconstruct. If True the out b0 and T2star values are smoothed Median filter and lowpassfiltering after which the water and fat maps are recomputed."

DixonFilterSize::usage = 
"DixonFilterSize is an options for DixonReconstruct. Defines the number of voxel with which the input b0 and T2star values are smoothed."

DixonIterations::usage = 
"DixonIterations is an options for DixonReconstruct. Defines the maximum itterations the fit can use."


MonitorUnwrap::usage = 
"MonitorUnwrap is an option for Unwrap. Monitor the unwrapping progress."

UnwrapDimension::usage = 
"UnwrapDimension is an option for Unwrap. Can be \"2D\" or \"3D\". 2D is for unwarpping 2D images or unwrapping the individual images from a 3D dataset \
(does not unwrap in the slice direction). 3D unwraps a 3D dataset in all dimensions."

UnwrapThresh::usage = 
"UnwrapThresh is an option for Unwrap. Is a value between 0.6 and 0.9, and defines when to unwrap, the higher the value the less unwrapping will be done."


(* ::Subsection::Closed:: *)
(*Error Messages*)


Unwrap::data2D = "Unwrapping Dimensions is 2D and this can only be preformed on 2D or 3D data. This data is `1`D."

Unwrap::data3D = "Unwrapping Dimensions is 3D and this can only be preformed on 3D data. This data is `1`D."

Unwrap::dim = "Unwrapping Dimensions can be \"2D\" or \"3D\", current value is `1`."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*DixonToPercent*)


SyntaxInformation[DixonToPercent] = {"ArgumentsPattern" -> {_, _}};

DixonToPercent[water_, fat_] := Block[{atot, fatMap, waterMap, fmask, wmask, back, afat, awater},
	 afat = Abs[fat];
	 awater = Abs[water];
	 atot = Abs[fat + water];
	 (*define water and fat maps*)
	 fatMap = Clip[Chop[DevideNoZero[afat, atot]], {0., 1.}];
	 waterMap = Clip[Chop[DevideNoZero[awater, atot]], {0., 1.}];
	 (*see where fat>50%*)
	 fmask = Mask[fatMap, .5];
	 wmask = 1 - fmask;
	 (*define background*)
	 back = Mask[fatMap + waterMap, .01];
	 fatMap = (wmask fatMap + fmask (1 - waterMap));
	 waterMap = (fmask waterMap + wmask (1 - fatMap));
	 (*output*)
	 N[{back waterMap, back fatMap}]
 ]

(* ::Subsection:: *)
(*DixonReconstruct*)


(* ::Subsubsection::Closed:: *)
(*DixonReconstruct*)


Options[DixonReconstruct] = {
	DixonPrecessions -> -1, 
	DixonFieldStrength -> 3, 
	DixonNucleus -> "1H",
	DixonFrequencies -> {{0}, {3.8, 3.4, 3.13, 2.67, 2.46, 1.92, 0.57, -0.60}},
	DixonAmplitudes -> {{1}, {0.089, 0.598, 0.047, 0.077, 0.052, 0.011, 0.035, 0.066}},
	DixonIterations -> 15, 
	DixonTollerance -> 0.001, 
	DixonMaskThreshhold -> 0.01,
	DixonFilterInput -> False, 
	DixonFilterOutput -> True, 
	DixonFilterSize -> 1};

SyntaxInformation[DixonReconstruct] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

DixonReconstruct[real_, imag_, echo_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, 0, 0, opts]

DixonReconstruct[real_, imag_, echo_, b0i_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, b0i, 0, opts]

DixonReconstruct[real_, imag_, echoi_, b0i_, t2_, OptionsPattern[]] := Block[{
	freqs, amps, gyro, precession, field, sigFW, sigPhi, eta, maxItt, r2Star,
	thresh, complex, ydat, result, input, b0f, b0, iopPhase, Amat, Amat2,
	cWater, cFat, b0fit, t2Star, fraction, signal, fit, itt, dim, mask,
	msk, t2i, t2f, echo, iop, ioAmat, phiEst, phiIn, phiInit, res, r2star, fsize, 
	r2, r2f ,dep, range, settings},
	
	(*algorithems are base on: *)	
	(*Triplett WT et.al. MRM 2014;72:8-19 doi 10.1002/mrm.23917 - fat peaks*)
	(*Huanzhou et.al. DOI 10.1002/jmri.21090 - T2* correction*)
	(*Reeder et.al. doi.org/10.1002/mrm.20624 - iDEAL*)
	(*Bydder et.al. doi:10.1016/j.mri.2010.08.011 - initial phase*)
	
	(*fixed setting*)
	echo = echoi;
	precession = OptionValue[DixonPrecessions](*-1,1*);
	field = OptionValue[DixonFieldStrength];
	freqs = precession field GyromagneticRatio[OptionValue[DixonNucleus]] OptionValue[DixonFrequencies];
	amps = #/Total[#] &/@ OptionValue[DixonAmplitudes];
	eta = OptionValue[DixonTollerance];
	maxItt = OptionValue[DixonIterations];
	thresh = OptionValue[DixonMaskThreshhold];
	fsize = OptionValue[DixonFilterSize];
	
	(*define in out phase*)
	Amat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ echo;
	iop = {0, 0.5}/Abs[Total[amps[[2]] freqs[[2]]]];
	ioAmat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ iop;
		
	(*create complex data for fit*)
	complex = N[real + imag I];
	If[ArrayDepth[real] === 4, complex = Transpose[complex]];
	range = {0., 2 Max[Abs[complex]]};
	dim = Dimensions[complex][[2;;]];
	dep = {ArrayDepth[complex]-1};
	
	(*find background voxels*)
	(*mask = Closing[Times @@ UnitStep[Abs[complex] - thresh], 1];*)
	mask = Closing[UnitStep[Abs[First@complex] - thresh], 1];
	
	(*prepare b0map and r2 map*)
	b0 = If[b0i === 0, ConstantArray[0., dim], b0i];
	r2 = If[t2 === 0,  ConstantArray[0., dim], DevideNoZero[1., t2]];
	
	(*smooth maps if needed*)
	{b0f, r2f} = If[OptionValue[DixonFilterInput],
		PrintTemporary["Filtering input B0 and T2* maps "];
		{If[b0i=!=0, LapFilter[b0], b0], If[t2=!=0, LapFilter[r2], r2]},
		{b0 ,r2}
	];
	
	(*perform the dixon reconstruction*)
	PrintTemporary["performing dixon iDEAL reconstruction"];
	
	(*define complex field map*)
	phiInit = mask(2 Pi I b0f - r2f);
	complex=TransData[complex,"l"];
	input = TransData[{complex, phiInit, mask}, "l"];
	
	(*Perform the dixon fit*)
	Quiet@Monitor[ii=0;result =Map[(ii++;DixonFiti[#, echo, Amat, {eta, maxItt}])&, input, dep];,ProgressIndicator[ii, {0, Times @@ dim}]];
 	{cWater, cFat, phiEst ,phiIn, res, itt} = TransData[Chop[result],"r"];

	(*filter the output*) 
	 If[OptionValue[DixonFilterOutput],
	 	PrintTemporary["Filtering field estimation and recalculating signal fractions"];
	 	(*smooth b0 field and R2star maps*)
	 	phiEst = mask(LapFilter[Im[phiEst]] I - LapFilter[Abs[-Re[phiEst]]]);
	 	phiIn = MedianFilter[#, 1] & /@ phiIn;
	 	
	 	(*recalculate the water fat signals*)
	 	input = TransData[{complex, phiEst, phiIn, mask}, "l"];
	 	Monitor[jj=0;result = Map[(jj++;DixonFiti[#, echo, Amat])&, input, dep];,ProgressIndicator[jj, {0, Times @@ dim}]];
		
		{cWater, cFat, res} = TransData[Chop[result],"r"];
	 ]; 	 
	 
	 
	 
	 (*create the output*)
	 PrintTemporary["performing water fat calculation"];
	 fraction = DixonToPercent[cWater, cFat];

	 (*signal and in/out phase data *)
	 signal = 1000 (Clip[Abs[{cWater, cFat}],range]/range[[2]]);
	 iopPhase = 1000 (Clip[TransData[InOutPhase[phiEst, iop, ioAmat, cWater, cFat], "r"],range]/range[[2]]);
	 res = 1000 Abs[res]/range[[2]];
	 
	 (*estimate b0 and t2star*)
	 b0fit = Im[phiEst]/(2 Pi);
	 r2Star = -Re[phiEst];
	 t2Star = DevideNoZero[1,r2Star];
	 fit = {
	 	Clip[b0fit, {-400., 400.}, {-400., 400.}], 
	 	Clip[t2Star, {0., 0.25}, {0., 0.25}], 
	 	Clip[r2Star, {0., 1000.}, {0., 1000.}],
	 	phiIn
	 	};

	 (*give the output*)
	 {fraction, signal, iopPhase, fit, itt, res}

 ]


(* ::Subsubsection::Closed:: *)
(*DixonFit*)


DixonFiti[{ydat_, phiInit_, mask_}, echo_, Amat_, {eta_, maxItt_}] := Block[
	{continue, phiEst, phiMat, pAmat, phivec, phi0, yfit, cFrac, Bmat, deltaPhi, i, iophiMat, iopImag, sol, res},
	If[mask>0,
		(*initialize fit*)
		deltaPhi = 0.;
		phiEst = phiInit;
		(*perform itterative fit*)
		i = 0;
		continue = True;
		While[continue,
			(*update the field map*)
			phiEst = phiEst + deltaPhi;
			(*find solution for complex fractions*)
			pAmat = Chop[Exp[phiEst echo] Amat];
			(*estimate the initial phase*)
			phi0 = DixPhaseEstimate[pAmat, ydat];
			(*find solution for complex fractions*)
			yfit = ydat Exp[-I phi0];
			cFrac = DixLeastSquaresP[pAmat, yfit];
			(*calculate solution and residuals*)
			sol = Chop[pAmat.cFrac];
			res = yfit - sol;
			(*calculate field map error*)
			Bmat = Join[Transpose[{echo sol}], pAmat, 2];
			deltaPhi = First@DixLeastSquaresC[Bmat, res];
			(*chech for continue*)
			i++;
			i++; continue = ! (Abs[deltaPhi] < eta || i >= maxItt);
		];
		
		(*give output*)
		{cFrac[[1]], cFrac[[2]], phiEst, phi0, RootMeanSquare[res], i}
		,
		{0.,0.,0.,0.,0.,0.}
	]
  ]


DixonFiti[{ydat_, phiInit_, phi0_, mask_}, echo_, Amat_] := Block[{pAmat, cFrac, yfit, res},
	If[mask > 0,
		(*find solution for complex fractions with smooth phase map*)
		pAmat = Chop[Exp[phiInit echo] Amat];
		yfit = ydat Exp[-I phi0];
		cFrac = DixLeastSquaresC[pAmat, yfit];
		(*calculate the residuals*)
		res = RootMeanSquare[yfit - Chop[pAmat.cFrac]];
		(*ouput*)
		{cFrac[[1]], cFrac[[2]], res},
		{0., 0., 0.}
	]
]


InOutPhase = Compile[{{phi, _Complex, 0}, {iop, _Real, 1}, {ioAmat, _Complex, 2}, {cWat, _Complex, 0}, {cFat, _Complex, 0}},
   Abs[(Exp[phi iop] ioAmat).{cWat, cFat}], 
   RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}, Parallelization -> True];


DixLeastSquaresP = Compile[{{A, _Complex, 2}, {y, _Complex, 1}}, Block[{AT, ATA, IATA},
	    AT = ConjugateTranspose[A];
	    ATA = Chop[AT.A];
	    IATA = If[Total[Flatten[Abs[Chop[ATA, 10^-6]]]] > 0., Chop[Inverse[ATA]], 0. ATA];
	    Re[IATA].Re[AT.y]
    ],
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


DixLeastSquaresC = Compile[{{A, _Complex, 2}, {y, _Complex, 1}}, Block[{AT, ATA, IATA},
	    AT = ConjugateTranspose[A];
	    ATA = Chop[AT.A];
	    IATA = If[Total[Flatten[Abs[Chop[ATA, 10^-6]]]] > 0., Chop[Inverse[ATA]], 0. ATA];
	    IATA.(AT.y)
    ], 
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


DixPhaseEstimate = Compile[{{A, _Complex, 2}, {y, _Complex, 1}}, Block[{AT, ATA, IATA},
		AT = ConjugateTranspose[A];
	    ATA = Chop[AT.A];
	    IATA = If[Total[Flatten[Abs[Chop[ATA, 10^-6]]]] > 0., Chop[Inverse[ATA]], 0. ATA];
	    .5 Arg[(AT.y).Re[IATA].(AT.y)]
    ], 
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*SimulateDixonSignal*)


Clear[SimulateDixonSignal]

Options[SimulateDixonSignal] = {
	DixonNucleus -> "1H", 
	DixonPrecessions -> -1, 
	DixonFieldStrength -> 3, 
	DixonFrequencies -> {{0}, {3.8, 3.4, 3.13, 2.67, 2.46, 1.92, 0.57, -0.60}}, 
	DixonAmplitudes -> {{1}, {0.089, 0.598, 0.047, 0.077, 0.052, 0.011, 0.035, 0.066}}
}

SyntaxInformation[SimulateDixonSignal] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

SimulateDixonSignal[echo_, fr_, B0_, T2_, OptionsPattern[]] := 
 Block[{precession,field,freqs,amps, Amat, phi, sig},
  precession = OptionValue[DixonPrecessions](*-1,1*);
  field = OptionValue[DixonFieldStrength];
  freqs = precession field GyromagneticRatio[OptionValue[DixonNucleus]] OptionValue[DixonFrequencies];
  amps = #/Total[#] & /@ OptionValue[DixonAmplitudes];
  
  Amat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ echo;
  phi = N@2 Pi B0 I - 1./T2;
  
  sig = Exp[phi echo] Amat;
  sig = sig.{fr, 1 - fr};
  
  {Re[sig], Im[sig]}
  ]


(* ::Subsection:: *)
(*Phase unwrap*)


(* ::Subsubsection::Closed:: *)
(*Unwrap*)


Options[Unwrap]={MonitorUnwrap->True,UnwrapDimension->"2D", UnwrapThresh->0.75};

SyntaxInformation[Unwrap] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

(* Phase unwrapping algorithem based on M.A. Herraez et al 2002 and Abdul-Rahman 2007.
unwraping algorithem. Unwraps one image, needs functions: Diff, SecDiff, EdgeReliability and PairCoor*)

Unwrap[dat_,OptionsPattern[]]:= Block[{data, ind, undim, out, mon, thresh, len},
	undim = OptionValue[UnwrapDimension];
	mon = OptionValue[MonitorUnwrap];
	thresh = Clip[OptionValue[UnwrapThresh],{0.6, 0.9}];
	
	Switch[undim,
		
		"2D",
		data = N[dat];
		len = Length[data];
		Which[
			MatrixQ[data],
			If[mon,PrintTemporary["Unwrapping one image using 2D algorithm."]];
			out = Unwrapi[data, thresh];
			,
			ArrayQ[data,3],
			If[mon,PrintTemporary["Unwrapping ",Length[data]," images using 2D algorithm"]];
			Monitor[
				out = MapIndexed[( ind = First[#2]; Unwrapi[#1, thresh] )&, data, {ArrayDepth[data]-2} ];
				out = UnwrapZi[out, thresh];
			 ,If[mon,ProgressIndicator[ind, {0, len}],""]
			]
		];,
		
		"3D",
		data = N[ArrayPad[dat,1]];
		If[ArrayQ[data,3],
			If[mon,PrintTemporary["Unwrapping 3D data using 3D algorithm"]];
			out = ArrayPad[Unwrapi[data, thresh, mon],-1];
			,
			Message[Unwrap::data3D,ArrayDepth[data]]
			];,
		_,
		Message[Unwrap::dim,undim]
		];
		
	(*center around 0*)
	out
]


(* ::Subsubsection::Closed:: *)
(*UnwrapSplit*)


Options[UnwrapSplit] = Options[UnwrapSplit]

SyntaxInformation[UnwrapSplit] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

UnwrapSplit[phase_, mag_,opts:OptionsPattern[]] := Block[{cutVal, phaseSplit, B0split},
	cutVal = CutData[mag][[3]];
	phaseSplit = CutData[phase, cutVal][[1 ;; 2]];
	B0split = Unwrap[#, opts] & /@ phaseSplit;
	StichData @@ B0split
]


(* ::Subsubsection::Closed:: *)
(*UnwrapZi*)


UnwrapZi[data_, thresh_]:= Block[{mask,slice,diff,meandiff,steps,off,unwrap,dat,num2,Roundi},
	
	Roundi = If[Negative[#], 
		num2=#-Ceiling[#]; If[-1<num2<-thresh,Floor[#],Ceiling[#]],
		num2=#-Floor[#]; If[1>num2>thresh,Ceiling[#],Floor[#]]]&;
	
	mask=Unitize[data];
	slice=Round[0.5Length[data]];
	diff=(#[[1]]-#[[2]]&/@ Partition[data/(2Pi),2,1]);
	
	meandiff=Median[#]&/@ Map[DeleteCases[Flatten[N[#]],0.]&,diff];
	steps=FoldList[Plus,0,Map[Roundi[#]&,meandiff]];
	off=Round[Median[DeleteCases[Flatten[N[data[[slice]]/(2Pi)]],0.]]];
	
	unwrap=steps-(steps[[slice]]+off);
	dat=(2Pi unwrap+data)mask//N;
	(dat - mask Round[MeanNoZero[Flatten[dat]],2Pi])
]


(* ::Subsubsection::Closed:: *)
(*Unwrapi*)


Unwrapi[dat_, thresh_] := Unwrapi[dat, thresh, False]

Unwrapi[dat_, thresh_, mon_] := Block[{data,datai,mask, crp, dimi, sorted,groups,groupsize,groupnr,task},
	(*monitor*)
	task="Preclustering data.";
	If[mon,PrintTemporary[Dynamic[task]]];
	
	(*rescale the dat to integers allows for faster matrix replacements*)
	datai = Round[10000. * dat / (2. Pi)];
	dimi = Dimensions[datai];
	
	(*remove zeros*)
	data = If[ArrayDepth[datai] == 3, crp = FindCrop[datai]; ApplyCrop[datai, crp], datai];
	
	(*make mask to pervent unwrapping in background*)
	mask = ArrayPad[Closing[ArrayPad[Mask[Ceiling[Abs@data], 1], 5], 1], -5];
	
	(*Get the edges sotrted for reliability and precluster groups*)
	sorted = GetEdgeList[data, mask];
	{groups, groupsize, groupnr} = MakeGroups[data, mask];

	(*make 2D data 3D and define shifts in add*)
	If[ArrayDepth[data] == 2,
		groups = {groups}; data = {data};
		sorted = {#[[1]] + 1, 1, #[[2]], #[[3]]} & /@ sorted;
		];
	
	(*Unwrap the data*)
	task="Unwrapping edges.";
	data = UnWrapC[sorted, data, groups, groupsize, groupnr, Round[thresh 10000]];
	
	(*make output in rad*)	
	If[ArrayDepth[dat] == 2, 
		(*output the 2D in rad*)
		2 Pi data[[1]]/10000.,
		(*align to zero and ouput 3D in rad*)
		data = 2 Pi (data - mask Round[MeanNoZero[Flatten[data]],10000])/10000.;
		ReverseCrop[data, dimi, crp]
		]
]


UnWrapC = Compile[{{sorted, _Integer, 2}, {datai, _Integer, 3}, {groupsi, _Integer, 3}, {groupsizei, _Integer, 1}, {groupnri, _Integer, 0}, {thresh,_Integer,0}},
	Block[{data, const, dir, dim, groups, group1, group2, groupsize, groupnr, z1, z2, x1, x2, y1, y2, wrap, wrapT, pos, g1, g2, out, adds,add, diff},
   	
    groups = groupsi;
    groupsize = groupsizei;
    data = datai;
    groupnr = groupnri;

    dim = Dimensions[data];
    adds = {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};
    z1=0;x1=0;y1=0;z2=0;x2=0;y2=0;
    group1=0;group2=0;g1=0;g2=0;

    out = Map[(

        (*get the voxel corrdinates and the neighbour, contrain to dimensions*)
        add = adds[[#[[1]]]];
        z1=#[[2]]; z2=z1+add[[1]]; If[z2>dim[[1]],z2=dim[[1]]];
        x1=#[[3]]; x2=x1+add[[2]]; If[x2>dim[[2]],x2=dim[[2]]];
        y1=#[[4]]; y2=y1+add[[3]]; If[y2>dim[[3]],y2=dim[[3]]];
                
        (*Get the group numbers*)
        group1 = groups[[z1, x1, y1]];
        group2 = groups[[z2, x2, y2]];
        
        (*unwrapping logic*)
        (*0. check if both are not in background if one is background skip*)
        If[! (group1 == 0 || group2 == 0),
         (*1. both already in same group but not zero (most cases 60%) do nothing*)
         If[(! (group1 > 1 && group1 == group2)),
          (*If not in the same group determine the wrap of the edge and determine how to unwrap*)
          diff = data[[z1, x1, y1]] - data[[z2, x2, y2]];
          (*wrap = Round[Round[diff,5000]-Sign[diff], 10000];*)
          wrap = Sign[diff] Ceiling[Abs[diff] - thresh, 10000];
          wrapT = (wrap != 0);
          
          (*get group sizes*)
          g1 = groupsize[[group1]];
          g2 = groupsize[[group2]];
          
          Which[
           (*2. one of two pixels already in group othter not in group (second most cases 30%)*)
           (*2A. group 1 existst add group 2 to group 1*)
           group2 == 1 && group1 > 1, 
           (*add the non group voxel to the existing group*)
           groups[[z2, x2, y2]] = group1;
           groupsize[[group1]] += 1;
           If[wrapT, data[[z2, x2, y2]] += wrap];
           ,
           (*2B. group 2 existst add group 1 to group 2*)
           group1 == 1 && group2 > 1, 
           (*add the non group voxel to the existing group*)
           groups[[z1, x1, y1]] = group2;
           groupsize[[group2]] += 1;
           If[wrapT, data[[z1, x1, y1]] -= wrap];
           ,
           (*3. both belong to no group (third most cases 6%)*)
           group1 == group2 == 1, 
           (*unwrap right or botom pixel and assign both to group*)
           groupnr++;
           groups[[z1, x1, y1]] = groups[[z2, x2, y2]] = groupnr;
           groupsize[[groupnr]] = 2;
           If[wrapT, data[[z2, x2, y2]] += wrap];
           ,
           (*4. both already in a group (least cases 4%)*)
           (*4A. group 1 is greather than or eaqual to group 2,
           add group 2 to group 1*)
           g1 >= g2, 
           (*unwrap group 2 with respect to group 1,
           only do if wrap\[NotEqual]0*)
           pos = 1 - Unitize[groups - group2];
           If[wrapT, data += wrap pos];
           groups += ((group1 - group2) pos);
           groupsize[[group1]] += g2;
           groupsize[[group2]] = 0;
           ,
           (*4B. group 2 is greather than or eaqual to group 1,
           add group 1 to group 2*)
           g2 > g1, 
           (*unwrap group 1 with respect to group 2,
           only do if wrap\[NotEqual]0*)
           pos = 1 - Unitize[groups - group1];
           If[wrapT, data -= wrap pos];
           groups += ((group2 - group1) pos);
           groupsize[[group2]] += g1;
           groupsize[[group1]] = 0;
           ]]];
        1
        ) &, sorted];
    data],
    {{dim, _Integer, 1},{add, _Integer, 1}, {group1, _Integer, 0}, {group2, _Integer, 0}, {g1, _Integer, 0}, {g2, _Integer, 0}, 
    	{z1, _Integer, 0}, {x1, _Integer, 0}, {y1, _Integer, 0}, {z2, _Integer, 0}, {x2, _Integer, 0}, {y2, _Integer, 0}},
    RuntimeOptions -> "Speed", Parallelization -> True];



(* ::Subsubsection::Closed:: *)
(*GetKernels*)


GetKernels[dep_] := Block[{ker, i, j, k, keri, kers},
   Switch[dep,
    2,
    ker = ConstantArray[0, {3, 3}];
    ker[[2, 2]] = -1;
    kers = ({i, j} = #; keri = ker; keri[[i, j]] = 1; keri) & /@ {
       {2, 1}, {2, 3}, {1, 2}, {3, 2}, {1, 1}, {3, 3}, {1, 3}, {3, 1}
       }(*H,V,D1,D2*),
    3,
    ker = ConstantArray[0, {3, 3, 3}];
    ker[[2, 2, 2]] = -1;
    kers = ({i, j, k} = #; keri = ker; keri[[i, j, k]] = 1; 
        keri) & /@ {
       {2, 1, 2}, {2, 3, 2}, {1, 2, 2}, {3, 2, 2}, {2, 2, 1}, {2, 2, 3},(*H,V,N*)
       {1, 1, 2}, {3, 3, 2}, {1, 3, 2}, {3, 1, 2},(*plane 1*)
       {1, 2, 1}, {3, 2, 3}, {1, 2, 3}, {3, 2, 1},(*plane 2*)
       {2, 1, 1}, {2, 3, 3}, {2, 1, 3}, {2, 3, 1},(*plane 3*)
       {1, 1, 1}, {3, 3, 3}, {1, 1, 3}, {3, 3, 1}, {1, 3, 1}, {3, 1, 3}, {3, 1, 1}, {1, 3, 3}(*diagonals*)
       }];
   Total /@ Partition[kers, 2]
   ];


(* ::Subsubsection::Closed:: *)
(*GetEdgeList*)


GetEdgeList[data_] := GetEdgeList[data, 1]

GetEdgeList[data_, maski_] := Block[{dep, diff, mask, edge, coor, fedge, ord, pos},
	dep = ArrayDepth[data];
	(*maske a mask if needed*)
	mask = If[maski === 1, Closing[Mask[Ceiling[Abs@data], 1], 1], maski];
	(*calculate the second order diff*)
	diff = ListConvolve[#, data, ConstantArray[2, dep], 0] & /@ GetKernels[dep];
	diff = If[dep == 2, DiffC2[diff], DiffC3[diff]];
	
	(*get the edge reliability*)
	edge = Switch[dep,
		2(*2D data*), 
		N@{(RotateLeft[#] & /@ diff) + diff, RotateLeft[diff] + diff},
		3(*3D data*), 
		N@{RotateLeft[diff] + diff, ((RotateLeft[#] & /@ #) & /@ diff) + diff, (RotateLeft[#] & /@ diff) + diff}
	];
	edge = mask # & /@ edge;
	 
	(*sort the edges for reliability*)
	coor = MapIndexed[#2 &, edge, {dep + 1}];
	fedge = Flatten[edge, dep];
	ord = Ordering[fedge];
	pos = Position[Unitize[fedge[[ord]]], 1, 1, 1][[1, 1]];
	Flatten[coor, dep][[ord]][[pos ;;]]
  ]

DiffC2 = Compile[{{diff, _Real, 3}}, Total[(diff - Round[diff, 10000])^2]];
DiffC3 = Compile[{{diff, _Real, 4}}, Total[(diff - Round[diff, 10000])^2]];


(* ::Subsubsection::Closed:: *)
(*MakeGroups*)


MakeGroups[data_] := MakeGroups[data, 1]
MakeGroups[data_, maski_]:=Block[{dep,dim,fun,min,max,part,dat,masks,mclus,clus,groupsize,groups,groupnr,mask},
	(*get data properties*)
	dep=ArrayDepth[data];
	dim=Dimensions[data];
	fun=If[dep==2,Image,Image3D];
	
	(*maske a mask if needed*)
	mask=If[maski===1,Closing[Mask[Ceiling[Abs@data],1],1],maski];
	
	(*find mask ranges*)
	{min,max}=MinMax[data];
	part=Partition[Range[min,max,(max-min)/12]//N,2,1];
	
	(*remove background form masks, and create masks*)
	dat=data/. 0->(-2min);
	masks=Mask[dat,#]&/@part;
	
	(*make groups from masks*)
	mclus=0;
	groups=Total[(
		(*make groups from masks and only keeps groups*)
		clus=MorphologicalComponents[DeleteSmallComponents[fun[#]]];
		clus=Clip[clus,{0,1}](mclus+clus);
		(*find max group number and export*)
		mclus=Max[clus];clus
	)&/@masks];
	groups=groups+mask;
	
	(*create outputs, the size vector and group nrs*)
	groupsize=ConstantArray[0,Count[Flatten[groups],1]];
	(groupsize[[#]]=Count[groups,#,2])&/@Range[1,Max[groups]];
	groupnr=Max[groups];
	{groups,groupsize,groupnr}
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
