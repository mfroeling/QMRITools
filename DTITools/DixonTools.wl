(* ::Package:: *)

(* ::Title:: *)
(*DTITools DixonTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`DixonTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DixonToPercent::usage = 
"DixonToPercent[water, fat] converts the dixon water and fat data to percent maps.

Output is {waterFraction, fatFraction}."

DixonReconstruct::usage = 
"DixonReconstruct[real, imag, echo] reconstruxt Dixon data with initital guess b0 = 0 and T2star = 0.
DixonReconstruct[real, imag, echo, b0] reconstructs Dixon data with intitial guess T2star = 0.
DixonReconstruct[real, imag, echo, b0, t2] reconstructs Dixon data.

real is the real data in radials.
imag is the imaginary data in radians.
b0 can be estimated from two phase images using Unwrap.
t2 can be estimated from multiple echos using T2fit.

Output is {{watF,fatF},{watSig,fatSig},{inphase,outphase},{b0,t2star},itterations}."


(* ::Subsection::Closed:: *)
(*Options*)


DixonPrecessions::usage = 
"DixonPrecessions is an options for DixonReconstruct. Defines the rotation of the signal {-1,1} default is -1."

DixonFieldStrength::usage = 
"DixonFieldStrength is an options for DixonReconstruct. Defines the fieldstrengths on which the data was acquired."

DixonFrequencies::usage = 
"DixonFrequencies is an options for DixonReconstruct. Defines the frequencies of the fat peaks being used."

DixonAmplitudes::usage = 
"DixonAmplitudes is an options for DixonReconstruct. Defines the amplitudes of the fat peaks being used."

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


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*DixonToPercent*)


SyntaxInformation[DixonToPercent] = {"ArgumentsPattern" -> {_, _}};

DixonToPercent[water_, fat_] := Block[{atot, fatMap, waterMap, fmask, wmask, back, afat, awater},
 	afat = Abs[fat];
 	awater = Abs[water];
 	atot = afat + awater;
 	
 	(*define water and fat maps*)
 	fatMap =  DevideNoZero[afat, atot];
 	waterMap =  DevideNoZero[awater, atot];
 	
 	(*see where fat > 50%*)
 	fmask = Mask[fatMap, .5];
 	wmask = 1 - fmask;

 	(*define background*)
 	back = Mask[fatMap + waterMap, .1];
 	(*define water and fat fractions*)
 	fatMap = (fmask fatMap + wmask (1-waterMap));
 	Clip[N[Chop[{back (1 - fatMap), back fatMap},10^-6]],{0.,1.},{0.,1.}]
  ]


(* ::Subsection:: *)
(*DixonReconstruct*)


(* ::Subsubsection::Closed:: *)
(*DixonReconstruct*)


Options[DixonReconstruct] = {DixonPrecessions -> -1, DixonFieldStrength -> 3, 
  DixonFrequencies -> {{0}, {3.8, 3.4, 3.13, 2.67, 2.46, 1.92, 0.57, -0.60}}, 
  DixonAmplitudes -> {{1}, {0.089, 0.598, 0.047, 0.077, 0.052, 0.011, 0.035, 0.066}}, 
  DixonIterations -> 5, DixonTollerance -> 0.01, DixonMaskThreshhold -> 0.05, 
  DixonFilterInput -> True, DixonFilterOutput -> True, DixonFilterSize -> 1};

SyntaxInformation[DixonReconstruct] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

DixonReconstruct[real_, imag_, echo_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, 0, 0, opts]

DixonReconstruct[real_, imag_, echo_, b0i_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, b0i, 0, opts]

DixonReconstruct[real_, imag_, echoi_, b0i_, t2_, OptionsPattern[]] := Block[{
	freqs, amps, gyro, precession, field, sigFW, sigPhi, eta, maxItt, r2Star,
	thresh, complex, ydat, result, input, b0f, b0, iopPhase, Amat, Amat2,
	cWater, cFat, b0fit, t2Star, fraction, signal, fit, itt, dim, mask,
	msk, t2i, t2f, echo, iop, ioAmat, phiEst, phiInit, res, r2star, fsize, 
	r2, r2f ,dep, range, settings},
	
	(*{3.80,3.40,2.60,1.94,0.39,-0.60} and {0.087,0.693,0.128,0.004,0.039,0.048}*)
	(*{3.8,3.4,3.11,2.67,2.45,-0.61} and {0.088,0.635,0.071,0.096,0.068,0.042};*)
	(*{3.8,3.4,3.13,2.67,2.46,1.92,0.57,-0.60} and {0.089,0.598,0.048,0.077,0.052,0.011,0.035,0.066};*)
	(*Triplett WT et.al. MRM 2014;72:8-19 doi 10.1002/mrm.23917*)
	
	(*fixed setting*)
	echo = echoi;
	precession = OptionValue[DixonPrecessions](*-1,1*);
	field = OptionValue[DixonFieldStrength];
	freqs = precession field 42.58 OptionValue[DixonFrequencies];
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
	complex = N[real + imag I];If[ArrayDepth[real] === 4, complex = Transpose[complex]];
	range = {0., 1.5Max[Abs[complex]]};
	dim = Dimensions[complex][[2;;]];
	dep = {ArrayDepth[complex]-1};
	
	(*find background voxels*)
	mask = Closing[Times @@ UnitStep[Abs[complex] - thresh], 1];
	
	(*prepare b0map and r2 map*)
	b0 = If[b0i === 0, ConstantArray[0., dim], b0i];
	r2 = If[t2 === 0,  ConstantArray[0., dim], DevideNoZero[1., t2]];
	(*smooth maps if needed*)
	{b0f, r2f} = If[OptionValue[DixonFilterInput],
		PrintTemporary["Filtering input B0 and T2* maps "];
		{If[b0i=!=0, LapFilt[b0], b0], If[t2=!=0, LapFilt[r2], r2]},
		{b0 ,r2}
	];
	
	(*perform the dixon reconstruction*)
	PrintTemporary["performing dixon iDEAL reconstruction"];
	
	(*define complex field map*)
	phiInit = mask(2 Pi I b0f - r2f);
	complex=TransData[complex,"l"];
	input = TransData[{complex, phiInit, mask}, "l"];
	(*Perform the dixon fit*)
	Monitor[ii=0;result =Map[(ii++;DixonFiti[#, echo, Amat, {eta, maxItt}])&, input, dep];,ProgressIndicator[ii, {0, Times @@ dim}]];
 	{cWater, cFat, phiEst ,res, itt} = TransData[Chop[result],"r"];

	(*filter the output*) 
	 If[OptionValue[DixonFilterOutput],
	 	PrintTemporary["Filtering field estimation and recalculating signal fractions"];
	 	(*smooth b0 field and R2star maps*)
	 	phiEst = mask(LapFilt[Im[phiEst]] I - LapFilt[Clip[-Re[phiEst],{0,500}]]);
	 	(*recalculate the water fat signals*)
	 	input = TransData[{complex, phiEst, mask}, "l"];
	 	Monitor[jj=0;result = Map[(jj++;DixonFiti[#, echo, Amat])&, input, dep];,ProgressIndicator[jj, {0, Times @@ dim}]];
		{cWater, cFat, res} = TransData[Chop[result],"r"];
	 ]; 	 
	 
	 (*create the output*)
	 PrintTemporary["performing water fat calculation"];
	 fraction = DixonToPercent[cWater, cFat];

	 (*signal and in/out phase data *)
	 signal = 1000 Clip[Abs[{cWater, cFat}],range];
	 iopPhase = 1000 Clip[TransData[InOutPhase[phiEst, iop, ioAmat, cWater, cFat], "r"],range];
	 res = 1000 Abs[res];
	 
	 (*estimate b0 and t2star*)
	 b0fit = Im[phiEst]/(2 Pi);
	 r2Star = -Re[phiEst];
	 t2Star = DevideNoZero[1,r2Star];
	 fit = {Clip[b0fit, {-400., 400.}, {-400., 400.}], Clip[t2Star, {0., 0.25}, {0., 0.25}], Clip[r2Star, {0., 1000.}, {0., 1000.}]};

	 (*give the output*)
	 {fraction, signal, iopPhase, fit, itt, res}

 ]

LapFilt[data_, fil_:0.8] := Clip[Chop[ImageData[TotalVariationFilter[Image3D[N@data, "Real"], fil, Method -> "Laplacian", MaxIterations -> 15]]], MinMax[data]]


(* ::Subsubsection::Closed:: *)
(*DixonFit*)


DixonFiti[{ydat_, phiInit_, mask_}, echo_, Amat_, {eta_, maxItt_}] := Block[
	{continue, phiEst, phiMat, pAmat, phivec, cFrac, Bmat, deltaPhi, i, iophiMat, iopImag, sol, res},
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
			cFrac = LeastSquares[pAmat, ydat, Tolerance -> 10^-6];
			(*calculate solution and residuals*)
			sol = Chop[pAmat.cFrac];
			res = ydat - sol;
			(*calculate field map error*)
			Bmat = Join[Transpose[{echo sol}], pAmat, 2];
			deltaPhi = First@LeastSquares[Bmat, res];
			(*chech for continue*)
			i++;
			i++; continue = ! (Abs[deltaPhi] < eta || i >= maxItt);
		];
		
		(*give output*)
		{cFrac[[1]], cFrac[[2]], phiEst, RootMeanSquare[res], i}
		,
		{0.,0.,0.,0.,0.}
	]
  ]


DixonFiti[{ydat_, phiInit_, mask_}, echo_, Amat_] := Block[{pAmat, cFrac},
	If[mask > 0,
		(*find solution for complex fractions with smooth phase map*)
		pAmat = Chop[Exp[phiInit echo] Amat];
		cFrac = LeastSquares[pAmat, ydat];
		(*calculate the residuals*)
		{cFrac[[1]], cFrac[[2]], RootMeanSquare[ydat - (pAmat.cFrac)]},
		{0., 0., 0.}
	]
]


InOutPhase = Compile[{{phi, _Complex, 0}, {iop, _Real, 1}, {ioAmat, _Complex, 2}, {cWat, _Complex, 0}, {cFat, _Complex, 0}},
   Abs[(Exp[phi iop] ioAmat).{cWat, cFat}], 
   RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}, Parallelization -> True];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
