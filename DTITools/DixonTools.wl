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
  DixonIterations -> 20, DixonTollerance -> 0.1, 
  DixonMaskThreshhold -> 0.05, DixonFilterInput -> False, DixonFilterOutput -> True, 
  DixonFilterSize -> 1};

SyntaxInformation[DixonReconstruct] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

DixonReconstruct[real_, imag_, echo_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, 0, 0, opts]

DixonReconstruct[real_, imag_, echo_, b0i_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, b0i, 0, opts]

DixonReconstruct[real_, imag_, echoi_, b0i_, t2_, OptionsPattern[]] := Block[{
	freqs, amps, gyro, precession, field, sigFW, sigPhi, eta, maxItt, r2Star,
	thresh, complex, ydat, result, input, b0f, b0, inphase, outphase, Amat, Amat2,
	cWater, cFat, b0fit, t2Star, fraction, signal, fit, itt, dim, mask,
	msk, t2i, t2f, echo, iop, ioAmat, phiEst, phiInit, res, r2star,
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
	iop = {1, 1.5}/Abs[freqs[[2, First@First@Position[amps[[2]], Max@amps[[2]]]]]];
	ioAmat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ iop;
		
	(*create complex data for fit*)
	complex = If[ArrayDepth[real] === 4, TransData[Transpose[N[real + imag I]], "l"], TransData[N[real + imag I], "l"]];
	range = {0., 1.5Max[Abs[complex]]};
	dim = Dimensions[complex][[;; -2]];
	dep = {ArrayDepth[complex] - 1};
	
	(*find background voxels*)
	mask = Times @@ TransData[UnitStep[Abs[complex] - thresh], "r"];
	
	(*prepare b0map and r2 map*)
	b0 = If[b0i === 0, ConstantArray[0., dim], b0i];
	r2 = If[t2 === 0,  ConstantArray[0., dim], DevideNoZero[1., t2]];
	{b0f, r2f} = If[OptionValue[DixonFilterInput],
		PrintTemporary["Filtering input B0 and T2* maps "];
		{MedianFilter[b0,1],MedianFilter[r2,1]},
		{b0 ,r2}];

	(*perform the dixon reconstruction*)
	PrintTemporary["performing dixon iDEAL reconstruction"];
	(*define complex field map*)
	phiInit = 2 Pi I b0f - r2f;
	input = TransData[{complex, phiInit, mask}, "l"];
	settings= {{echo, iop}, {Amat, ioAmat}, {eta, maxItt}};
	
	System`SetSystemOptions["CheckMachineUnderflow" -> False];
	(*monitor Calculations*)
	ii=j=0;
	PrintTemporary[ProgressIndicator[Dynamic[ii], {0, Times @@ dim}]];
	result =Map[(j++;If[j>500,ii+=j;j=1;];DixonFiti[#, settings])&, input, dep];
	System`SetSystemOptions["CheckMachineUnderflow" -> True];
 	
 	(*get the result*)
 	{cWater, cFat, phiEst, inphase, outphase, itt ,res} = TransData[Chop[result],"r"];

	(*filter the output*) 
	 If[OptionValue[DixonFilterOutput],
	 	PrintTemporary["Filtering field estimation and recalculating fractions"];
	 	(*smooth b0 field and R2star maps*)
	 	phiEst = MedianFilter[Im[phiEst],1] I - MedianFilter[Clip[-Re[phiEst],{0,500}],1];
	 	(*recalculate phi and dixon estimation*)
	 	input = TransData[{complex, phiEst, mask}, "l"];
	 	(*recalculate the water fat signals*)
	 	result = Map[Dixoni[#, {{echo, iop}, {Amat, ioAmat}}]&, input, dep];
		{cWater, cFat, inphase, outphase ,res} = TransData[Chop[result],"r"];
	 ];
	 
	 (*create the output*)
	 PrintTemporary["performing water fat calculation"];
	 fraction = DixonToPercent[cWater, cFat];
	 
	 (*estimate b0 and t2star*)
	 b0fit = Im[phiEst]/(2 Pi);
	 r2Star = -Re[phiEst];
	 t2Star = DevideNoZero[1,r2Star];
	 fit = {Clip[b0fit, {-400., 400.}, {-400., 400.}], Clip[t2Star, {0., 0.25}, {0., 0.25}], Clip[r2Star, {0., 1000.}, {0., 1000.}]};
	 
	 (*signal and in/out phase data *)
	 signal = 1000 Clip[Abs[{cWater, cFat}],range];
	 {inphase, outphase} = 1000 Clip[{inphase, outphase},range];
	 res = 1000 Abs[res];

	 (*give the output*)
	 {fraction, signal,  {inphase, outphase}, fit, itt, res}

 ]


(* ::Subsubsection::Closed:: *)
(*DixonFit*)


DixonFiti[{ydat_, phiInit_, mask_}, {{echo_, iop_}, {Amat_, ioAmat_}, {eta_, maxItt_}}] := Block[
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
			pAmat = Exp[phiEst echo] Amat;
			cFrac = LeastSquares[pAmat, ydat, Tolerance -> 10^-6];
			sol = pAmat.cFrac;
			(*calculate field map error*)
			Bmat = Join[Transpose[{echo sol}], pAmat, 2];
			deltaPhi = First@LeastSquares[Bmat, (ydat - sol), Tolerance -> 10^-6];
			(*chech for continue*)
			i++;
			continue = ! ((Abs[Re[deltaPhi]] < eta && (Abs[Im[deltaPhi]]) < eta) || i >= maxItt);
		];
		
		(*calculate the final values*)
		(*calculate in out phase images*)
		iopImag = Abs[(Exp[phiEst iop] ioAmat).cFrac];
		res = Sqrt[Mean[(ydat-sol)^2]];
		
		(*give the output,comp_Water,comp_Fat,b0fit,t2star*)
		Flatten[{cFrac, phiEst, iopImag, i, res}]
		,
		{0.,0.,0.,0.,0.,0.,0.}
	]
  ]

Dixoni[{ydat_, phiInit_, mask_}, {{echo_, iop_}, {Amat_, ioAmat_}}] := Block[
	{pAmat, pAmat2, cFrac, iophiMat, iopImag, sol, res},
	If[mask>0,
		(*find solution for complex fractions with smooth phase map*)
		pAmat = Exp[phiInit echo] Amat;
		cFrac = LeastSquares[pAmat, ydat, Tolerance -> 10^-6];	
		(*calculate the residuals*)
		res = Sqrt[Mean[(ydat-(pAmat.cFrac))^2]];
		(*calculate in out phase images*)
		iopImag = Abs[(Exp[phiInit iop] ioAmat).cFrac];
		(*give the output,comp_Water,comp_Fat,b0fit,t2star*)
		Flatten[{cFrac, iopImag, res}]
		,
		{0.,0.,0.,0.,0.}
	]
  ]

(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
