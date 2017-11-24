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

DixonFilterInputSize::usage = 
"DixonFilterInputSize is an options for DixonReconstruct. Defines the number of voxel with which the input b0 and T2star values are smoothed."

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

DixonToPercent[water_, fat_] := 
 Block[{atot, fatMap, waterMap, fmask, wmask, back, afat, awater},
 	
 	afat = Abs[fat];
 	awater = Abs[water];
 	atot = Abs[fat + water];
 	
 	(*calcualte fat fractions*)
 	fatMap =  DevideNoZero[afat, afat + awater];
 	
 	(*see where fat > 50%*)
 	fmask = Mask[fatMap, .5];
 	wmask = 1 - fmask;
 	
 	(*define water and fat maps*)
 	fatMap =  DevideNoZero[afat, atot];
 	waterMap =  DevideNoZero[awater, atot];
 	(*define background*)
 	back = Mask[fatMap + waterMap, .1];
 	(*define water and fat fractions*)
 	fatMap = (fmask fatMap + wmask (1-waterMap));
 	Clip[N[{back (1 - fatMap), back fatMap}],{0.,1.}]
  ]


(* ::Subsection:: *)
(*DixonReconstruct*)


(* ::Subsubsection::Closed:: *)
(*DixonReconstruct*)


Options[DixonReconstruct] = {DixonPrecessions -> -1, DixonFieldStrength -> 3, 
  DixonFrequencies -> {{0}, {3.8, 3.4, 3.13, 2.67, 2.46, 1.92, 0.57, -0.60}}, 
  DixonAmplitudes -> {{1}, {0.089, 0.598, 0.048, 0.077, 0.052, 0.011, 0.035, 0.066}}, 
  DixonIterations -> 50, DixonTollerance -> 0.1, 
  DixonMaskThreshhold -> 0.05, DixonFilterInput -> False, DixonFilterOutput -> True, 
  DixonFilterInputSize -> 2};

SyntaxInformation[DixonReconstruct] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

DixonReconstruct[real_, imag_, echo_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, 0, 0, opts]

DixonReconstruct[real_, imag_, echo_, b0_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, b0, 0, opts]

DixonReconstruct[real_, imag_, echoi_, b0_, t2_, OptionsPattern[]] := Block[{
	freqs, amps, gyro, precession, field, sigFW, sigPhi, eta, maxItt, r2Star,
	thresh, complex, ydat, result, input, b0f, b0i, inphase, outphase, Amat,
	cWater, cFat, b0fit, t2Star, fraction, signal, fit, itt, dim, mask,
	msk, t2i, t2f, echo, iop, ioAmat, phiEst, phiInit, res, r2star},
	(*{3.80,3.40,2.60,1.94,0.39,-0.60} and {0.087,0.693,0.128,0.004,0.039,0.048}*)
	(*{3.8,3.4,3.11,2.67,2.45,-0.61} and {0.088,0.635,0.071,0.096,0.068,0.042};*)
	(*{3.8,3.4,3.13,2.67,2.46,,1.92,0.57,-0.60} and {0.089,0.598,0.048,0.077,0.052,0.011,0.035,0.066};*)
	(*Triplett WT et.al. MRM 2014;72:8-19 doi 10.1002/mrm.23917*)
	
	(*fixed setting*)
	echo = echoi;
	precession = OptionValue[DixonPrecessions](*-1,1*);
	gyro = 42.58(*Hz*);
	field = OptionValue[DixonFieldStrength];
	{freqs, amps} = {precession field gyro OptionValue[DixonFrequencies], OptionValue[DixonAmplitudes]};
	Amat = Transpose[Map[Total, Transpose[(amps Exp[freqs (2 Pi I) #] & /@ echo)], {2}]];
	
	(*define in out phase*)
	iop = {1, 1.5}/Abs[freqs[[2, First@First@Position[amps[[2]], Max@amps[[2]]]]]];
	ioAmat = Transpose[Map[Total, Transpose[(amps Exp[freqs (2 Pi I) #] & /@ iop)], {2}]];
	
	(*define stop criterea*)
	eta = OptionValue[DixonTollerance];
	maxItt = OptionValue[DixonIterations];
	thresh = OptionValue[DixonMaskThreshhold];
	
	(*create complex data for fit*)
	complex = If[ArrayDepth[real] === 4,
		TransData[Transpose[N[real + imag I]], "l"],
		TransData[N[real + imag I], "l"]
		];
	dim = Dimensions[complex][[;; -2]];
	mask = Times @@ TransData[UnitStep[Abs[complex] - thresh], "r"];
	
	(*prepare b0map*)
	b0f = If[b0 === 0,
		ConstantArray[0., dim],
		If[OptionValue[DixonFilterInput], LowpassFilter[#,.75]&/@b0, b0]
		];
	(*prepare t2Star map*)
	t2f = If[t2 === 0,
		ConstantArray[0., dim],
		If[OptionValue[DixonFilterInput], LowpassFilter[#,.75]&/@t2, t2]
		];
		
	(*define complex field map*)
	phiInit = b0f + I DevideNoZero[1, t2f]/(2 Pi);
	
	(*monitor Calculations*)
	PrintTemporary["performing dixon iDEAL reconstruction"];
	i=j=0;
	SetSharedVariable[i];ParallelEvaluate[j=0];
	PrintTemporary[ProgressIndicator[Dynamic[i], {0, Times @@ dim}]];
	
	(*make parallel*)
	DistributeDefinitions[echo, iop, Amat, ioAmat, eta, maxItt, DixonFiti, Dixoni];
	(*perform the dixon reconstruction*)
	
	input = TransData[{complex, phiInit, mask}, "l"];
	result = ParallelMap[(
		j++;If[j>1000,i+=j;j=1;];
		DixonFiti[#, {echo, iop}, {Amat, ioAmat}, {eta, maxItt}]
		) &, input, {ArrayDepth[complex] - 1}];
 	
 	{cWater, cFat, phiEst, inphase, outphase, itt ,res} = TransData[result,"r"];

	(*filter the output*) 
	 If[OptionValue[DixonFilterOutput],
	 	PrintTemporary["Filtering field estimation and recalculating fractions"];
	 	(*smooth b0 field map*)
	 	b0fit = MedianFilter[Re[phiEst],1];
	 	b0fit = LowpassFilter[#, 0.75] & /@ b0fit;
	 	(*smooth R2star map*)
	 	r2star = Clip[-Re[2 Pi I result[[All, All, All, 3]]], {0., 1000.}, {0., 0.}];
	 	r2star = MedianFilter[r2star,1];
	 	r2star = LowpassFilter[#, 0.75] & /@ r2star;
	 	(*recalculate phi and dixon estimation*)
	 	phiEst = b0fit + I r2star/(2 Pi); 
	 	input = TransData[{complex, phiEst, mask}, "l"];
	 	result = ParallelMap[(Dixoni[#, {echo, iop}, {Amat, ioAmat}]) &, input, {ArrayDepth[complex] - 1}];
	 	(*output fit*)
		{cWater, cFat, inphase, outphase ,res} = TransData[result,"r"];
	 ];
	 
	 (*create the output*)
	 PrintTemporary["performing water fat calculation"];
	 fraction = DixonToPercent[cWater, cFat];
	 (*estimate b0 and t2star*)
	 b0fit = Re[phiEst];
	 r2Star = (-Re[2 Pi I phiEst]);
	 t2Star = DevideNoZero[1,r2Star];
	 fit = {Clip[b0fit, {-400., 400.}, {0., 0.}], Clip[t2Star, {0., 0.25}, {0., 0.}], Clip[r2Star, {0., 1000.}, {0., 0.}]};
	  (*signal and in/out phase data *)
	 signal = Clip[Abs[{cWater, cFat}],{1,1.5}MinMax[Abs[complex]]];
	 {inphase, outphase} = Clip[{inphase, outphase},{1,1.5}MinMax[Abs[complex]],{0.,0.}];
	 
	 (*give the output*)
	 {fraction, 1000 signal, 1000 {inphase, outphase}, fit, itt, 1000 Abs[res]}

 ]


(* ::Subsubsection::Closed:: *)
(*DixonFit*)


DixonFiti[{ydat_, phiInit_, mask_}, {echo_, iop_}, {Amat_, ioAmat_}, {eta_, maxItt_}] := Block[
	{continue, phiEst, phiMat, pAmat, phivec, cFrac, Bmat, deltaPhi, i, iophiMat, iopImag, sol, res},
	If[mask>0,
		(*initialize fit*)
		phiEst = phiInit;
		deltaPhi = 0.;
		
		(*perform itterative fit*)
		i = 0;
		continue = True;
		While[continue,
			(*update the field map*)
			phiEst = phiEst + deltaPhi;
			(*find solution for complex fractions*)
			pAmat = DiagonalMatrix[Exp[2 Pi I phiEst echo]].Amat;
			cFrac = LeastSquares[pAmat, ydat, Tolerance -> 10^-8];
			sol = pAmat.cFrac;
			(*calculate field map error*)
			phivec = (2 Pi I echo (sol));
			Bmat = Transpose[{phivec, pAmat[[All, 1]], pAmat[[All, 2]]}];
			deltaPhi = First@LeastSquares[Bmat, (ydat - sol), Tolerance -> 10^-8];
			(*chech for continue*)
			i++;
			continue = ! ((Abs[Re[deltaPhi]] < eta && (Abs[Im[deltaPhi]]) < eta) || i >= maxItt);
		];
		
		(*calculate the final values*)
		(*calculate in out phase images*)
		iophiMat = DiagonalMatrix[Exp[2 Pi I phiEst iop]];
		iopImag = Abs[iophiMat.ioAmat.cFrac];
		
		res = Sqrt[Mean[(ydat-sol)^2]];
		
		(*give the output,comp_Water,comp_Fat,b0fit,t2star*)
		Flatten[{cFrac, phiEst, iopImag, i ,res}]
		,
		{0.,0.,0.,0.,0.,0.,0.}
	]
  ]

Dixoni[{ydat_, phiInit_, mask_}, {echo_, iop_}, {Amat_, ioAmat_}] := Block[
	{pAmat, cFrac, iophiMat, iopImag, sol, res},
	If[mask>0,
		(*find solution for complex fractions*)
		pAmat = DiagonalMatrix[Exp[2 Pi I phiInit echo]].Amat;
		cFrac = LeastSquares[pAmat, ydat, Tolerance -> 10^-8];
		sol = pAmat.cFrac;		
		
		(*calculate in out phase images*)
		iophiMat = DiagonalMatrix[Exp[2 Pi I phiInit iop]];
		iopImag = Abs[iophiMat.ioAmat.cFrac];
		
		(*calculate the residuals*)
		res = Sqrt[Mean[(ydat-sol)^2]];
		
		(*give the output,comp_Water,comp_Fat,b0fit,t2star*)
		Flatten[{cFrac, iopImag ,res}]
		,
		{0.,0.,0.,0.,0.}
	]
  ]

(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
