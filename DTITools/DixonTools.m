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

Unprotect @@ Names["DTITools`DixonTools`*"];
ClearAll @@ Names["DTITools`DixonTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


DixonToPercent::usage = 
"DixonToPercent[water, fat] converts the dixon water and fat data to percent maps."

DixonReconstruct::usage = 
"DixonReconstruct[real, imag, echo, b0]."


(* ::Subsection:: *)
(*General Options*)


DixonPrecessions::usage = 
"DixonPrecessions"

DixonFieldStrength::usage = 
"DixonFieldStrength"

DixonFrequencies::usage = 
"DixonFrequencies"

DixonAmplitudes::usage = 
"DixonAmplitudes"

DixonTollerance::usage = 
"DixonTollerance"

DixonMaskThreshhold::usage = 
"DixonMaskThreshhold"

DixonFilterB0::usage = 
"DixonFilterB0"

DixonFilterB0Size::usage = 
"DixonFilterB0Size"



(* ::Subsection:: *)
(*Error Messages*)





(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*DixonToPercent*)


SyntaxInformation[DixonToPercent] = {"ArgumentsPattern" -> {_, _}};

DixonToPercent[water_, fat_] := 
 Block[{tot, fatMap, waterMap, fmask, wmask, back, afat, awater},
 	(*define water and fat signals*)
 	afat=Abs[fat];
 	awater=Abs[water];
 	tot = Abs[fat + water];
 	
 	(*calcualte water and fat fractions*)
 	fatMap =  DevideNoZero[afat, tot];
 	waterMap = (1 - DevideNoZero[awater, tot]);
 	
 	(*see where fat > 50%*)
 	fmask = Mask[afat, .5];
 	wmask = 1 - fmask;
 	back = 1 - Mask[waterMap + wmask, 1.9];
 	
 	(*define water and fat fractions*)
 	fatMap = Clip[(fmask fatMap + wmask waterMap),{0,1}];
 	N[{back (1 - fatMap), back fatMap}]
  ]


(* ::Subsection::Closed:: *)
(*DixonReconstruct*)


(* ::Subsubsection::Closed:: *)
(*DixonReconstruct*)


Options[DixonReconstruct] = {
	DixonPrecessions -> -1, DixonFieldStrength -> 3,
	DixonFrequencies -> {{0}, {3.8, 3.4, 3.11, 2.67, 2.45, -0.61}},
	DixonAmplitudes -> {{1}, {0.088, 0.635, 0.071, 0.096, 0.068, 0.042}},
	MaxIterations -> 300, DixonTollerance -> 0.01,
	DixonMaskThreshhold -> 0.05, DixonFilterB0 -> True,
	DixonFilterB0Size -> 2
	};

SyntaxInformation[DixonReconstruct] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

DixonReconstruct[real_, imag_, echo_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, 0, 0, opts]

DixonReconstruct[real_, imag_, echo_, b0_, opts : OptionsPattern[]] := DixonReconstruct[real, imag, echo, b0, 0, opts]

DixonReconstruct[real_, imag_, echo_, b0_, t2_, OptionsPattern[]] := Block[{freqs, amps, gyro, precession, field,
  sigFW, sigPhi, eta, maxItt, thresh, complex, ydat, result, b0f, b0i, inphase, outphase, freqsAmps, Amat,
  cWater, cFat, b0fit, t2Star, fraction, signal, fit, itt, dim, mask, msk,t2i,t2f},
 
 (*{3.80,3.40,2.60,1.94,0.39,-0.60} or {3.8,3.4,3.11,2.67,2.45,-0.61}*)
 (*{0.087,0.693,0.128,0.004,0.039,0.048} or {0.088,0.635,0.071,0.096,0.068,0.042};*)
 
 (*fixed setting*)
 precession = OptionValue[DixonPrecessions](*-1,1*);
 gyro = 42.58(*Hz*);
 field = OptionValue[DixonFieldStrength];
 freqsAmps = {freqs, amps} = {precession field gyro OptionValue[DixonFrequencies], OptionValue[DixonAmplitudes]};
 Amat = Transpose[Map[Total, Transpose[(amps Exp[freqs (2 Pi I) # ] & /@ echo)], {2}]];
 
 (*define stop criterea*)
 eta = OptionValue[DixonTollerance];
 maxItt = OptionValue[MaxIterations];
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
 	If[OptionValue[DixonFilterB0], GaussianFilter[b0, OptionValue[DixonFilterB0Size]], b0]
 	];
 (*prepare t2Star map*)
 t2f = If[t2 === 0, 
 	ConstantArray[0., dim],
 	If[OptionValue[DixonFilterB0], GaussianFilter[t2, OptionValue[DixonFilterB0Size]], t2]
 	];
 
 (*monitor Calculations*)
 j = 0;
 PrintTemporary["performing dixon iDEAL reconstruction"];
 PrintTemporary[ProgressIndicator[Dynamic[j], {0, Times @@ dim}]];
 
 (*make parallel*)
 ParallelEvaluate[j=0];
 SetSharedVariable[j];
 DistributeDefinitions[echo, Amat, freqsAmps, eta, maxItt, DixonFiti];
 
 (*perform the dixon reconstruction*)
 {cWater, cFat, b0fit, t2Star, inphase, outphase, itt} = TransData[Map[(
      j++; 
      If[#[[4]] > 0, 
      	DixonFiti[#[[1]], #[[2]], #[[3]], echo, Amat, freqsAmps, {eta, maxItt}], 
      	{0., 0., 0., 0., 0., 0., 0.}]
      ) &, TransData[{complex, b0f, t2f, mask}, "l"], {ArrayDepth[complex] - 1}],"r"];
 
 (*create the output*)
 fraction = DixonToPercent[cWater, cFat];
 signal = Abs[{cWater, cFat}];
 fit = {Clip[b0fit, {-400., 400.}, {0., 0.}], 
   Clip[1000 t2Star, {0., 200.}, {0., 0.}]};
 
 (*five the output*)
 {fraction, 1000 signal, 1000 {inphase, outphase}, fit, itt}
 ]


(* ::Subsubsection::Closed:: *)
(*DixonFit*)


DixonFiti[ydat_, b0_, t2_, echo_, Amat_, {freqs_, amps_}, {eta_, maxItt_}] := 
 Block[{continue, phiEst, phiMat, phivec, cFrac, Bmat, deltaPhi, i,iop, iophiMat, ioAmat, iopImag},
  
  (*initialize fit*)
  continue = True;
  phiEst = If[t2>0, b0 + I (1/t2)/(2 Pi), b0];
  i = 0;
  
  (*perform itterative fit*)
  While[continue,
   (*find solution for complex fractions*)
   phiMat = Quiet@Inverse[DiagonalMatrix[Exp[2 Pi I phiEst echo]]];
   cFrac = LeastSquares[Amat, phiMat.ydat];
   (*update the field map*)
   phivec = (2 Pi I echo (cFrac[[1]] + cFrac[[2]] Amat[[All, 2]]));
   Bmat = Transpose[{phivec, Amat[[All, 1]], Amat[[All, 2]]}];
   deltaPhi = LeastSquares[Bmat, (phiMat.ydat - Amat.cFrac)][[1]];
   phiEst = phiEst + deltaPhi;
   (*chech for continue*)
   i++;
   continue = ! ((Abs[Re[deltaPhi]] < eta  && Abs[Im[deltaPhi]] < eta) || i >= maxItt);
   ];
  (*calculate the final values*)
  phiMat = Inverse[DiagonalMatrix[Exp[2 Pi I phiEst echo]]];
  cFrac = LeastSquares[Amat, phiMat.ydat];
  
  (*calculate in out phase images*)
  iop = {1, 1.5}/Abs[freqs[[2, First@First@Position[amps[[2]], Max@amps[[2]]]]]];
  iophiMat = DiagonalMatrix[Exp[2 Pi I phiEst iop]];
  ioAmat = Transpose[Map[Total,Transpose[(amps Exp[freqs (2 Pi I) # ] & /@ iop)], {2}]];
  iopImag = Abs[cFrac.iophiMat.ioAmat];
  
  (*give the output, comp_Water,comp_Fat,b0fit,t2star*)
  Flatten[{cFrac, Re[phiEst], 1/(Abs[2 Pi Im[phiEst]]), iopImag, i}]
  ]



(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`DixonTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`DixonTools`*"];

EndPackage[]
