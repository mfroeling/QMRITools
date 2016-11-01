(* ::Package:: *)

(* ::Title:: *)
(*DTITools Tools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`IVIMTools`"];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`IVIMTools`*"];
ClearAll @@ Names["DTITools`IVIMTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


IVIMCalc::usage=
"IVIMCalc[data, binp, init] calculates the IVIM fit.
 
data should be 1D ,2D, 3D or 4D. 
binp should be full bmatrix which can be calculated from the bvecs en bvals using Bmatrix. 
init should are the initialization parameters for 2 components this is {S0, f, D, Dp} for 3 componentes this is {S0, f1, f2, D, Dp1, Dp2}.

output is {S0, f1, D, pD1} or {S0, f1, f2, D, pD1, pD2}."

IVIMFunction::usage =
"IVIMFunction[] gives the IVIM function with 2 comps.
IVIMFunction[components] gives the IVIM function.
IVIMFunction[components, type] gives the IVIM function. 

type can be \"Normal\" or \"Exp\".
componenets can be 2 or 3.

output is the function with b, S0, f1, f2, D, pD1, pD2 as parameters"

FracCorrect::usage = 
"FracCorrect[fraction, time] corrects the signal fraction calculated with the IVIM model for tissue relaxation and acquisition parameters.
After correction the signal fraction can be regarded as volume fraction.
FracCorrect[{fraction1, fraction2}, time] corrects the signal fraction1 and fraction2 from a 3 compartement IVIM model. 

time is {{te, tr}, {t2t, t21}, {t1t, t11}} or {{te, tr}, {t2t, t21, t22}, {t1t, t11, t12}}
where t2t and t1t are \"tissue\" relaxation times and t11 t12, t21 and t22 the \"fluid\" relaxation times 

output is the corrected fraction maps";

MeanNoZero::usage = 
"MeanNoZero[data] calculates the mean of the data ignoring the zeros."

ThetaConv::usage = 
"ThetaConv[{F1, Fc, pDc}] converts the parameters from Log space to normal space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConv[{F1, F2, Dc, pDc1}] converts the parameters from Log space to normal space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConv[{F1, F2, Dc, pDc1, pDc2}] converts the parameters from Log space to normal space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3."

ThetaConvi::usage = 
"ThetaConvi[{f, dc, pdc}] converts the parameters from Normal space to Log space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConvi[{f1, f2, dc, pdc1}] converts the parameters from Normal space to Log space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3.
ThetaConvi[{f1, f2, dc, pdc1, pdc2}] converts the parameters from Normal space to Log space. Is used in BayesianIVIMFit2 and BayesianIVIMFit3."

FConvert::usage = 
"FConvert[F] convers the fraction F from log space."

FConverti::usage = 
"FConverti[f] converts the fraction f to log space."

BayesianIVIMFit2::usage = 
"BayesianIVIMFit2[data, bval, init, mask] performs bayesian IVIM fit of data.

data is the data which should be {slice, Ndiff, x, y}.
bval is the bvector whould be length Ndiff.
init is the initalization of the bayesian fit which comes from IVIMCals, (without S0 using 2 compartments).
mask is the region in which the bayesian fit is performed.

output is {f1, dc, pdc1}."

BayesianIVIMFit3::usage = 
"BayesianIVIMFit3[data, bval, init, mask] performs bayesian IVIM fit of data.

data is the data which should be {slice, Ndiff, x, y}.
bval is the bvector whould be length Ndiff.
init is the initalization of the bayesian fit which comes from IVIMCals, (without S0 using 3 compartments).
mask is the region in which the bayesian fit is performed.

output is {f1, f2, dc, pdc1, pdc2}."

IVIMResiduals::usage = 
"IVIMResiduals[data, binp, pars] calculates the root mean square residuals of an IVIM fit ussing IVIMCalc, BayesianIVIMFit2 or BayesianIVIMFit3."

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

CorrectParMap::usage = 
"CorrectParMap[par, constraints, mask] removes the IVIM parameters outside the constraints within the mask.

par is {f1, dc, pdc1} or {f1, f2, dc, pdc1, pdc2}.
constraints are the lower and upper constraints for each parameters {{min, max},...}
mask has the same dimensions as the parameter maps. 

output are the corrected paremeter maps."

NNLeastSquares::usage = 
"NNLeastSquares[A, y] performs a Non Negative Linear Least Squares fit.
finds an x that solves the linear least-squares problem for the matrix equation A.x==y.

output is the solution x."

IVIMCorrectData::usage = 
"IVIMCorrectData[data, {S0, f, pdc}, bval] removes the ivim signal from the data.

data is the original data.
{S0, f, pdc} are the solution to a 2 compartment IVIM fit using IVIMCalc or BayesianIVIMFit2.
bval are the bvalues.

output is the corrected data."


(* ::Subsection:: *)
(*General Options*)


IVIMConstrained::usage = 
"IVIMConstrained is an option for IVIMCalc. When set True the fit wil be constrained to the values given in IVIMConstrains."

IVIMTensFit::usage =
"IVIMTensFit is an option for IVIMCalc. When set True the tissue diffusion component wil be calculated as a tensor."

IVIMComponents::usage =
"IVIMComponents is an option for IVIMCalc. Default value is 2, the tissue and the blood component. can also be set to 3."

IVIMConstrains::usage = 
"IVIMConstrains is an option for IVIMCalc. \
Default values are: {{0.8, 1.2}, {0, 1}, {0.0005, 0.0035}, {0.005, 0.5}, {0.002, 0.015}}. \
Where {{S0 in percentage},{fractions},{tissue diffusion},{blood compartment Dp},{third compartment}}."

IVIMFixed::usage =
"IVIMFixed is an option for IVIMCalc and the default value is False. \
When set True the pseudo diffusion wil be fixed to the parameter given as init.
When set to \"One\" only the fast component of a 3 compartment fit is fixed."

MonitorIVIMCalc::usage = 
"MonitorIVIMCalc is an option for IVIMCalc. When true the proceses of the calculation is shown."

ChainSteps::usage = 
"ChainSteps is an option for BayesianIVIMFit2 and BayesianIVIMFit3. It determines how long the algorithm runs.
three values must be given {itterations, burn steps, sample density}."

UpdateStep::usage = 
"UpdateStep is an option for BayesianIVIMFit2 and BayesianIVIMFit3. It determines how often the parameters are updated. Is optimized during the first 500 burn steps.";

FixPseudoDiff::usage = 
"FixPseudoDiff is an option for BayesianIVIMFit2 and BayesianIVIMFit3. If the pDc1 and pD2 were fixed in IVIMCalc this value should be True."

FixPseudoDiffSD::usage = 
"FixPseudoDiffSD is an option for BayesianIVIMFit2 and BayesianIVIMFit3. Gives the standard deviation of pDc1 and pD2 if FixPseudoDiff is True"

CorrectPar::usage = 
"CorrectPar is an option for BayesianIVIMFit2 and BayesianIVIMFit3. If True it removes the values outside the contraints using CorrectParMap"

FitConstrains::usage = 
"FitConstrains is an option for BayesianIVIMFit2 and BayesianIVIMFit3. Gives the contraints of the parameters. 
The values are used for displaying the histograms and for the initialization if CorrectPar is True"

OutputSamples::usage = 
"OutputSamples is an option for BayesianIVIMFit2 and BayesianIVIMFit3. If set True the full marcov chain is given as an additionaln output."

FilterMaps::usage = 
"FilterMaps is an option for IVIMCorrectData. If True the IVIM parameter maps are filtered before signal correction"

FilterType::usage = 
"FilterType is an option for IVIMCorrectData. If FilterMaps is True it tells which filter to use. can be \"Median\" of \"Gausian\""

FilterSize::usage = 
"FilterSize is an option for IVIMCorrectData. If FilterMaps is True it gives the kernel size."


(* ::Subsection:: *)
(*Error Messages*)


IVIMCalc::init = "Number of initialization values is `1` and should be `2`"

IVIMCalc::time = "Number of components is `1` but the number of relaxations times is `2`"

IVIMCalc::comp = "Number of components should be 2 or 3 not `1`"

IVIMCalc::bvec = "The length of the data (`1`) should be the same as the length of the bmatrix (`2`)"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*IVIMCalc*)


Options[IVIMCalc] = {Method -> Automatic, Parallelize->False, MonitorIVIMCalc -> True, 
	IVIMFixed -> False, IVIMConstrained -> True, IVIMTensFit -> False, IVIMComponents -> 2,
	IVIMConstrains -> {{0.8, 1.2}, {0, 1}, {0.0005, 0.0035}, {0.001, 0.5}, {0.001, 0.5}}
   };

SyntaxInformation[IVIMCalc] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

IVIMCalc[data_, binp_, init_, OptionsPattern[]] := 
 Module[{tensFit, components, fixed, constrained, method, bmdc, bin,fcon,
   s0min, s0max, fmin, fmax, dcmin, dcmax, pdc1min, pdc1max, pdc2min, pdc2max,
   depthD,dirD,dirB,func,dat,dat0,datn,rl,rr,ivim, mdat, sol, fitd, start, funcf,
   S0s, frin1, frin2, dcin, pdcin1, pdcin2, funcin, fixrule,cons,dccon,pdccon, fpars, mapfun, out
   },
  
  Clear[bm, bm1, bm2, bm3, bm4, bm5, bm6, S0, f1, f2, dc, pdc1, pdc2];
  
  (*data checks*)
  depthD = ArrayDepth[data];
  dirD = If[depthD == 4, Length[data[[1]]], Length[data]];
  dirB = Length[binp];
  
  (*check if data is 4D,3D,2D or 1D*)
  If[depthD > 4, Return[Message[IVIMCalc::data, ArrayDepth[data]]]];
  (*check if bmat is the same lengt as data*)
  If[dirB != dirD, Return[Message[IVIMCalc::bvec, dirD, dirB]]];
  
  (*perform tensor fit*)
  tensFit = OptionValue[IVIMTensFit];
  components = Clip[OptionValue[IVIMComponents], {1, 3}];
  fixed = OptionValue[IVIMFixed];
  constrained = OptionValue[IVIMConstrained];
  method = OptionValue[Method];
  
  {{s0min, s0max}, {fmin, fmax}, {dcmin, dcmax}, {pdc1min, pdc1max}, {pdc2min, pdc2max}} = OptionValue[IVIMConstrains];
  
  (*running parameters*)
  bmdc = If[tensFit, 
  	bm = bm1 + bm2 + bm3; {bm1, bm2, bm3, bm4, bm5, bm6}.{xx, yy, zz, xy, xz, yz}, 
  	Clear[bm]; bm* dc
  	];
  
  dat = N[If[depthD == 4, Transpose[data], data]];
  dat0 = DeleteCases[Flatten[dat[[{1}]]], 0.];
  mdat=Mean[dat0];
  datn = If[dat0 === {}, dat, dat/mdat];
  
  rl = RotateRight[Range[depthD]];
  rr = RotateLeft[Range[depthD]];
  
  (*contruct bvals for fit*)
  bin = If[!tensFit, If[VectorQ[binp], binp, Abs[Total[#[[1 ;; 3]]]] & /@ binp], binp];
  
 (*initial values for fit*)
  Switch[components,
   1, f1 = f2 = 0; {S0s,dcin} = init;,
   2, f2 = 0; {S0s, frin1, dcin, pdcin1} = init;,
   3, Clear[f2]; {S0s, frin1, frin2, dcin, pdcin1, pdcin2} = init;
  ];
  
  (*initial fit values*)
  funcin = Join[
  	If[components==1,{},{{f1, frin1}, {f2, frin2}}[[1 ;; components-1]]],
    If[tensFit, Thread[{{xx, yy, zz, xy, xz, yz}, dcin}], {{dc, dcin}}],
     (*if not fixed give fit start parameters for pdc values*)
     If[components==1,
     	{},
     	Switch[fixed,
     	False, {{pdc1, pdcin1}, {pdc2, pdcin2}}[[1 ;; components - 1]],
     	"One", {{pdc1, pdcin1}},
     	_, {}
     	]] 
  ];
 
  (*fix fixed parameters*)
  fixrule = If[components==1,
  	{},
  	Switch[fixed, 
  		True,{pdc1 -> pdcin1, pdc2 -> pdcin2}[[1 ;; components - 1]],
  		"One", {pdc2 -> pdcin2},
  		_, {}]];
  
  (*generate fix fuctions*)
  func =Chop[Simplify[(S0*((((1 - f1 - f2)*Exp[-bmdc]) + (f1* Exp[-bm pdc1]) + (f2* Exp[-bm pdc2]))))]] /. fixrule;
  
  (*constrains for fit*)
  cons = If[constrained,
	  (*constrains dc and tens*)
	  dccon = If[tensFit, {dcmin < xx < dcmax, dcmin < yy < dcmax, dcmin < zz < dcmax, dcmin < (xx + yy + zz)/3 < dcmax}, {dcmin < dc < dcmax}];
	  (*if not fixed dc and pdc need to be constrained*)
	  pdccon = If[components==1,
	  	{},
	  	Switch[fixed,
	  	False,{dc < pdc1, pdc1min < pdc1 < pdc1max, pdc2min < pdc2 < pdc2max}[[1 ;; components]],
	  	"One",{dc < pdc1, pdc1min < pdc1 < pdc1max},
	  	_,{}
	  ]];
	  (*if 3 components pdc1 and/or pdc2 also need to be constrained dc < pdc1 < pdc(in)2*)
	  If[components == 3 && (fixed === False), AppendTo[pdccon, pdc1 < pdc2]];
	  If[components == 3 && (fixed === "One"), AppendTo[pdccon, pdc1 < pdcin2]];
	  
	  fcon=If[components==1,{},{fmin < f1 < fmax, fmin < f2 < fmax}[[1 ;; components - 1]]];
	  (*all constrains together*)
	  Join[{(f1 + f2) < 1}, fcon, dccon, pdccon],
	  (*no constrains*)
	  {}
	];
	
   (*running parameters*)
   fpars = If[tensFit, {bm1, bm2, bm3, bm4, bm5, bm6}, {bm}];
   
   (*fit function with of without constrains*)
   funcf = If[constrained,{func,cons},func];
   
   (*define output*)
   out=Join[
   	{S0, f1, f2}[[1 ;; components]],
   	If[tensFit, {{xx, yy, zz, xy, xz, yz}}, {dc}], 
   	If[components==1,{},{pdc1, pdc2}[[1 ;; components - 1]]]
   	] /. fixrule;
  
  (*perform fit*)
  mapfun=If[OptionValue[Parallelize],
  	DistributeDefinitions[bin,funcin,funcf,fpars,method,out];
  	ParallelMap,Map];
  	 ivim = Transpose[mapfun[
  		(
  		S0s = #[[1]];
  		If[N[#] == #*0. || S0s == 0.,
		   (*masked voxel*)
		   0. out
		   ,
		   (*data voxel*)
		   fitd = Flatten /@ ({bin, #} // Transpose);
		   start=Prepend[funcin,{S0,S0s}];
		   sol = Quiet[FindFit[fitd, funcf, start , fpars, Method -> method, MaxIterations -> 150]];
		   out /. sol		    
  		]
  		
  		)&, 
	  	Transpose[datn, rl], {depthD - 1}], rr];  	
	  	
  ivim[[1]]=ivim[[1]]*mdat;
  ivim
  ]


(* ::Subsection::Closed:: *)
(*IVIMFunction*)


SyntaxInformation[IVIMFunction] = {"ArgumentsPattern" -> {_, _}};

IVIMFunction[]:=IVIMFunction[2, "Normal"]

IVIMFunction[pars_]:=IVIMFunction[pars, "Normal"]

IVIMFunction[pars_, fun_] := 
 Block[{func, Global`S0, Global`f1, Global`bm, Global`dc, Global`pdc1, Global`f2, Global`pdc2, ff1, ff2, fdc, fpdc1, 
   fpdc2},
  Switch[fun,
   "Normal", ff1 = Global`f1; ff2 = Global`f2; fdc = Global`dc; fpdc1 = Global`pdc1; fpdc2 = Global`pdc2;,
   "Exp", ff1 = Exp[Global`f1]/(1 + Exp[Global`f1]); ff2 = Exp[Global`f2]/(1 + Exp[Global`f2]); 
   fdc = Exp[Global`dc]; fpdc1 = Exp[Global`pdc1]; fpdc2 = Exp[Global`pdc2];
   ];
  
  func = Switch[pars,
    2, Simplify[Global`S0*(((1 - ff1)*Exp[-Global`bm fdc]) + (ff1*Exp[-Global`bm fpdc1]))],
    3, Simplify[
     Global`S0*( (1 - ff1 - ff2)*Exp[-Global`bm*fdc] + ff1*Exp[-Global`bm*fpdc1] + 
        ff2*Exp[-Global`bm*fpdc2] )]
    ]]


(* ::Subsection::Closed:: *)
(*FracCorrect*)


SyntaxInformation[FracCorrect] = {"ArgumentsPattern" -> {_, _, _.}};

FracCorrect[f1_, time_] := Module[{te, tr, t2t, t21, t1t, t11, st, s1},
  {{te, tr}, {t2t, t21}, {t1t, t11}} = time;
  st = Sigval[{1, t1t, t2t}, tr, te] // N;
  s1 = Sigval[{1, t11, t21}, tr, te] // N;
  ((f1*st)/(s1 - f1*s1 + f1*st))
  ]
  
FracCorrect[{f1_, f2_?VectorQ}, time_] := 
 Module[{te, tr, t2t, t21, t22, t1t, t11, t12, st, s1, s2},
  {{te, tr}, {t2t, t21, t22}, {t1t, t11, t12}} = time;
  st = Sigval[{1, t1t, t2t}, tr, te] // N;
  s1 = Sigval[{1, t11, t21}, tr, te] // N;
  s2 = Sigval[{1, t12, t22}, tr, te] // N;
  {(f1*s2*st)/(s1*s2 - f1*s1*s2 - f2*s1*s2 + f2*s1*st + f1*s2*st),
   (f2*s1*st)/(s1*s2 - f1*s1*s2 - f2*s1*s2 + f2*s1*st + f1*s2*st)}
  ]

(*correct fraction for T2 relaxation*)
Sigval[par_, TR_, TE_] := par[[1]] (1 - Exp[-TR/par[[2]]]) Exp[-TE/par[[3]]]


(* ::Subsection::Closed:: *)
(*MeanNoZero*)


SyntaxInformation[MeanNoZero] = {"ArgumentsPattern" -> {_, _.}};

Default[MeanNoZero] = 0;
MeanNoZero[data_, cor_.] := 
 Mean[DeleteCases[Flatten[N[data], ArrayDepth[data] - (cor + 1)], 0.]]


(* ::Subsection::Closed:: *)
(*ThetaConv*)


SyntaxInformation[ThetaConv] = {"ArgumentsPattern" -> {_}};

ThetaConv[{f1_, dc_, pdc_}] := {Exp[f1]/(1 + Exp[f1]), Exp[dc], Exp[pdc]};
ThetaConv[{f1_, f2_, dc_, pdc1_}] := {Exp[f1]/(1 + Exp[f1]), Exp[f2]/(1 + Exp[f2]), Exp[dc], Exp[pdc1]};
ThetaConv[{f1_, f2_, dc_, pdc1_, pdc2_}] := {Exp[f1]/(1 + Exp[f1]), Exp[f2]/(1 + Exp[f2]), Exp[dc], Exp[pdc1], Exp[pdc2]};


(* ::Subsection::Closed:: *)
(*ThetaConvi*)


SyntaxInformation[ThetaConvi] = {"ArgumentsPattern" -> {_}};

ThetaConvi[{F_, Dc_, pDc_}] := N[{Log[F] - Log[1 - F], Log[Dc], Log[pDc]}] /. {-Infinity -> 0., Indeterminate -> 0.};
ThetaConvi[{F1_, F2_, Dc_, pDc1_}] := N[{Log[F1] - Log[1 - F1], Log[F2] - Log[1 - F2], Log[Dc], Log[pDc1]}] /. {-Infinity -> 0., Indeterminate -> 0.};
ThetaConvi[{F1_, F2_, Dc_, pDc1_, pDc2_}] := N[{Log[F1] - Log[1 - F1], Log[F2] - Log[1 - F2], Log[Dc], Log[pDc1], Log[pDc2]}] /. {-Infinity -> 0., Indeterminate -> 0.};


(* ::Subsection::Closed:: *)
(*FConvert*)


SyntaxInformation[FConvert] = {"ArgumentsPattern" -> {_}};

FConvert[f_] := If[VectorQ[f], FConvf[f], FConv[f]]
FConv = Compile[{{f1, _Real, 3}}, Exp[f1]/(1 + Exp[f1]), Parallelization -> True];
FConvf = Compile[{{f1, _Real, 1}}, Exp[f1]/(1 + Exp[f1]), Parallelization -> True];


(* ::Subsection::Closed:: *)
(*FConverti*)


SyntaxInformation[FConverti] = {"ArgumentsPattern" -> {_}};

FConverti[F_] := If[VectorQ[F], FConvif[F], FConvi[F]];
FConvi = Compile[{{F, _Real, 3}}, Log[F] - Log[1 - F], Parallelization -> True];
FConvif = Compile[{{F, _Real, 1}}, Log[F] - Log[1 - F], Parallelization -> True];


(* ::Subsection::Closed:: *)
(*CorrectParMap*)


SyntaxInformation[CorrectParMap] = {"ArgumentsPattern" -> {_, _, _}};

CorrectParMap[par_, con_, mask_] := 
 Module[{dim, mean, cov, sig, clipmap, rand, clippar, parnew},

  {mean, cov} = MeanCov[Transpose@Data3DToVector[par,mask][[1]]];
  sig = Diagonal[cov];
  dim = Dimensions[First@par];
  
  MapThread[(
     clipmap = mask ((1 - Mask[#1, #2[[1]] + .001]) + Mask[#1, #2[[2]] - .001]);
     rand = clipmap RandomVariate[NormalDistribution[#3, #4], dim];
     clippar = (1 - clipmap) #1;
     parnew = rand + clippar
     ) &, {par, con, mean, sig}]
  ]


(* ::Subsection::Closed:: *)
(*BayesianIVIMFit2*)


Options[BayesianIVIMFit2] = {ChainSteps -> {20000, 1000, 10}, UpdateStep -> {0.5, 0.2, 0.5}, 
	FixPseudoDiff -> False, CorrectPar->True, 
   FixPseudoDiffSD -> 0.5, OutputSamples->False,
   FitConstrains -> ThetaConv[{{-7.6, 7.6}, {-10.0, -5.7}, {-7.0, 0.}}]
   };

SyntaxInformation[BayesianIVIMFit2] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

BayesianIVIMFit2[data_, bval_, fitpari_, maski_, opts : OptionsPattern[]] := Module[{
	useDat, thetai, fix, ynf, fixSD, out1, out2, h1, solution,mask,
	fitpar, deviation, con2, con2e, mui, covi, mmu, mcov,post
   },
  
  con2 = OptionValue[FitConstrains];
  con2e = ThetaConvi[con2];
  fix = OptionValue[FixPseudoDiff];
  fixSD = OptionValue[FixPseudoDiffSD];
  
  mask = Mask[data[[All, 1]], 0.000001]maski;

  fitpar=ThetaConvi[MapThread[N[mask Clip[#1, #2]] &, {fitpari, con2}]];
  fitpar=If[OptionValue[CorrectPar], CorrectParMap[fitpar, con2e, mask], fitpar];
    	
  useDat = MaskDTIdata[data, mask];
  {thetai,post} = Data3DToVector[fitpar,mask];
  thetai=Transpose@thetai;
  ynf = First@Data3DToVector[Transpose[data],mask];
  
  If[fix, 
  	thetai[[3]] = RandomVariate[NormalDistribution[Mean[thetai[[3]]], fixSD], Length[thetai[[3]]]]
  	];
  
  {mui, covi} = MeanCov[thetai];
  
  (*show the pre fit distribution*)
  Print[Dimensions[ynf],Dimensions[thetai]];
  PrintTemporary[h1 = HistogramPar[thetai, {con2e, 75, mui, covi}, 3, Gray, {0.1, 0.1, 0.1}]];    

  out2 = BayesianIVIMFitI2[thetai, bval, ynf, FilterRules[{opts}, Options[BayesianIVIMFitI2]]];
  solution = out2[[2]];
  deviation = out2[[3]];
  out1 = Chop[VectorToData[Transpose@ThetaConv[solution], post]];
  
  {mmu, mcov} = MeanCov[solution];
  Print[Column[{
     h1,
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
      gj = FunceC2[fj, dj, pdj, bval];
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
               (*{muj, covj} = MeanCov[{fj, dj, pdj}];
               icovj=N@PseudoInverse[covj];*)
               
               (*steps 3c "loop" over the voxels i, perform as vector for each of the parameters *)
               (*step 3c-i - define theta(j) - {fj,dj,pdj}=thetaj;*)
               
               (*step 3c-ii - ramom sample frtmp*)
               fjt = RandomNormalCf[fj, w1];
               gjt = FunceC2[fjt, dj, pdj, bval];
               bool1 = Quiet@AlphaC[{fj, dj, pdj}, {fjt, dj, pdj}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
               gj = BoolAdd[bool1, gj, gjt];
               fj = BoolAdd[bool1, fj, fjt];

               (*step 3c-iii - ramom sample dctmp*)
               djt = RandomNormalCd[dj, w2];
               gjt = FunceC2[fj, djt, pdj, bval];
               bool2 = Quiet@AlphaC[{fj, dj, pdj}, {fj, djt, pdj}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
               gj = BoolAdd[bool2, gj, gjt];
               dj = BoolAdd[bool2, dj, djt];
               
               (*step 3c-iv - ramom sample pdctmp *)
               pdjt = RandomNormalCd[pdj, w3];
               gjt = FunceC2[fj, dj, pdjt, bval];
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



(* ::Subsection::Closed:: *)
(*BayesianIVIMFit3*)


Options[BayesianIVIMFit3] = {ChainSteps -> {20000, 1000, 10}, 
	UpdateStep -> {0.5, 0.5, 0.1, 0.5, 0.5}, FixPseudoDiff -> False, 
	CorrectPar->True, OutputSamples->False, FixPseudoDiffSD -> 0.5, 
	FitConstrains ->ThetaConv[{{-7.6, 7.6}, {-7.6, 7.6}, {-10., -5.5}, {-6.5, -2.3}, {-5.2, 0.}}]};

SyntaxInformation[BayesianIVIMFit3] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

BayesianIVIMFit3[data_, bval_, fitpari_, maski_, opts : OptionsPattern[]] := 
 Module[{fitpar,con3,mask,
   useDat, thetai, ynf, fix, fixSD, out1, out2, h1,
   solution, deviation, con3e, mui, covi, mmu, mcov,post},
  
  con3 = OptionValue[FitConstrains];
  con3e = ThetaConvi[con3];
  fix = OptionValue[FixPseudoDiff];
  fixSD = OptionValue[FixPseudoDiffSD];
  
  mask = Mask[data[[All, 1]], 0.000001]maski;

  fitpar=ThetaConvi[MapThread[N[mask Clip[#1, #2]] &, {fitpari, con3}]];
  fitpar=If[OptionValue[CorrectPar], CorrectParMap[fitpar, con3e, mask], fitpar];
  
  useDat = MaskDTIdata[data, mask];
  {thetai,post} = Data3DToVector[fitpar,mask];
  thetai=Transpose@thetai;
  ynf = Data3DToVector[Transpose[data],mask][[1]];
  
  If[fix,
   thetai[[4]] = RandomVariate[NormalDistribution[Mean[thetai[[4]]], fixSD], Length[thetai[[4]]]];
   thetai[[5]] = RandomVariate[NormalDistribution[Mean[thetai[[5]]], fixSD], Length[thetai[[5]]]];
   ];
  
  {mui, covi} = MeanCov[thetai];
  
  (*show the pre fit distribution*)
  Print[Dimensions[ynf],Dimensions[thetai]];
  PrintTemporary[h1 = HistogramPar[thetai, {con3e, 75, mui, covi}, 3, Gray, {0.1, 0.1, 0.1, 0.1, 0.1}]];
    
  out2 = BayesianIVIMFitI3[thetai, bval, ynf, FilterRules[{opts}, Options[BayesianIVIMFitI3]]];
  solution = out2[[2]];
  deviation = out2[[3]];
  out1 = Chop[VectorToData[Transpose@ThetaConv[solution], post]];

  {mmu, mcov} = MeanCov[solution];
  Print[Column[{
     h1,
     HistogramPar[solution, {con3e, 75, mmu, mcov}, 3, 
      Blue, {0.1, 0.1, 0.1, 0.1, 0.1}],
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
      gj = FunceC3[f1j, f2j, dj, pd1j, pd2j, bval];
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
               {muj, covj, icovj} = 
                RandomGibsSample[{f1j, f2j, dj, pd1j, pd2j}, covj, nvox];
               
               (*steps 3c "loop" over the voxels i, perform as vector for each of the parameters *)
               
               (*step 3c-i - define theta(j) - {fj,dj,pdj}=thetaj;*)

			   (*step 3c-ii - ramom sample frtmp*)
               (*comp 1*)
               f1jt = RandomNormalCf[f1j, w1];
               gjt = FunceC3[f1jt, f2j, dj, pd1j, pd2j, bval];
               bool1 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1jt, f2j, dj,  pd1j, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
               gj = BoolAdd[bool1, gj, gjt];
               f1j = BoolAdd[bool1, f1j, f1jt];
               
               (*comp 2*)
               f2jt = RandomNormalCf[f2j, w2];
               gjt = FunceC3[f1j, f2jt, dj, pd1j, pd2j, bval];
               bool2 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2jt, dj, pd1j, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
               gj = BoolAdd[bool2, gj, gjt];
               f2j = BoolAdd[bool2, f2j, f2jt];
               
               (*step 3c-iii - ramom sample dctmp*)
               djt = RandomNormalCd[dj, w3];
               gjt = FunceC3[f1j, f2j, djt, pd1j, pd2j, bval];
               bool3 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2j, djt, pd1j, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
               gj = BoolAdd[bool3, gj, gjt];
               dj = BoolAdd[bool3, dj, djt];
               
               (*step 3c-iv - ramom sample pdctmp *)
               (*comp 1*)
               pd1jt = RandomNormalCd[pd1j, w4];
               gjt = FunceC3[f1j, f2j, dj, pd1jt, pd2j, bval];
               bool4 = AlphaC[{f1j, f2j, dj, pd1j, pd2j}, {f1j, f2j, dj, pd1jt, pd2j}, muj, icovj, yn, yty, gj, gjt, nbval, nvox];
               gj = BoolAdd[bool4, gj, gjt];
               pd1j = BoolAdd[bool4, pd1j, pd1jt];
               
               (*comp 2*)
               pd2jt = RandomNormalCd[pd2j, w5];
               gjt = FunceC3[f1j, f2j, dj, pd1j, pd2jt, bval];
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
                 ]];
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
         "     mu:  ", NumberForm[#, {5, 2}] & /@ muj // MatrixForm, 
         ",", 
         NumberForm[Round[#, .01], {5, 2}] & /@ ({100, 100, 1000, 1000, 1000} ThetaConv[muj]) // MatrixForm,
         "     cov:  ", NumberForm[#, {5, 2}] & /@ # & /@ covj // MatrixForm
         }]
       ](*close monitor*)
      ]];(*close Timing*)
      
   Print[PlotPerformance[{steps, nvox}, {t1, t2s}, {w, wstart}]];
   
   {
    theta,
    Table[Mean[theta[[All, i]]], {i, Length[theta[[1]]]}],
    Table[StandardDeviation[theta[[All, i]]], {i, Length[theta[[1]]]}],
    cov, mu, Length[theta]}
   
   ];



(* ::Subsection::Closed:: *)
(*IVIM Core Functions*)


BooleC = Compile[{{val, _Real, 1}, {ru, _Real, 1}}, UnitStep[val - ru], Parallelization -> True, RuntimeOptions -> "Speed"];
BooleC1 = Compile[{{val, _Real, 1}, {ru, _Real, 0}}, UnitStep[val - ru], Parallelization -> True, RuntimeOptions -> "Speed"];
BooleC2 = Compile[{{val, _Real, 1}, {min, _Real, 0}, {max, _Real, 0}}, UnitStep[val - min] (1. - UnitStep[val - max]), Parallelization -> True, RuntimeOptions -> "Speed"];
BoolAdd = N[#2 - #1 #2 + #1 #3] &;

ClipC = Compile[{{theta, _Real, 2}, {trans, _Integer, 0}}, Block[{chk1, out},
    If[Length[theta] == 3,
     chk1 = BooleC2[theta[[1]], -7.0, 7.0];
     chk1 = chk1*BooleC2[theta[[2]], -9.5, -5.5];
     (*chk1=chk1*(1-BooleC1[theta[[3]],-0.001]);*)
     chk1 = chk1*BooleC2[theta[[3]], -7.5, -0.001];
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
    ], Parallelization -> True, RuntimeOptions -> "Speed"];

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
      (Eigensystem[#][[2]].DiagonalMatrix[
         Max[#, tol] & /@ (eigsys[[1]])].Transpose[eigsys[[2]]])
      ) &, N[mat], (! PositiveDefiniteMatrixQ[#] &)]
   ];

(*
RandomGibsSample[theta_, cov_, m_] := Block[{munew, tm, icov,mat},
   
   	munew = N[RandomVariate[System`MultinormalDistribution[Mean /@ N[theta], N[PosSym[cov/m]]]]]; 
    munew = N[(1 + munew) - 1];
    tm = N[(ClipC[theta, 1] - munew)];
    mat=N[PosSym[PseudoInverse[tm.Transpose[tm]]]];
    icov = N[RandomReal[WishartDistribution[mat, m - 3]]];
    {munew, N[PseudoInverse[icov]], icov}
   ];
*)

RandomGibsSample[theta_, cov_, m_] := Block[{munew, tm, icov, mat,tmt,mi},
   munew = N[RandomVariate[MultinormalDistribution[Mean /@ N[theta],N[(*PosSym[cov/m]*)cov/m]]]];
   (*munew = N[(1 + munew) - 1];*)
   tm = Chop[N[(ClipC[theta, 1] - munew)]];
   tmt=Chop[N[tm.Transpose[tm]], 10^-5];
   (*mi=m-3;*)
   icov = N[RandomVariate[InverseWishartMatrixDistribution[m-3, tmt]]];
   {munew, icov, N@PseudoInverse[icov]}
];

(*random Normal Samplers*)
RandomNormalC = Compile[{{m, _Real, 1}, {s, _Real, 1}},
   Chop[MapThread[RandomVariate[NormalDistribution[#1, #2^2]] &, {m, s}]],
   Parallelization -> True, RuntimeOptions -> "Speed"];
RandomNormalCf = Compile[{{m, _Real, 1}, {s, _Real, 1}},
   Chop[Clip[MapThread[RandomVariate[NormalDistribution[#1, #2^2]] &, {m, s}],{-9.,9.}]],
   Parallelization -> True, RuntimeOptions -> "Speed"];
RandomNormalCd = Compile[{{m, _Real, 1}, {s, _Real, 1}},
   Chop[Clip[MapThread[RandomVariate[NormalDistribution[#1, #2^2]] &, {m, s}],{-15.,0.}]],
   Parallelization -> True, RuntimeOptions -> "Speed"];

(*calulated fitted points g(fr, dc, pdc)*)
FunceC2 = Compile[{{fr, _Real, 1}, {dc, _Real, 1}, {pdc, _Real, 1}, {bm, _Real, 1}},
   Chop[Transpose[Map[((Exp[Exp[dc] #] + Exp[Exp[pdc] # + fr])/(1 + Exp[fr])) &, -bm]]]
   , Parallelization -> True, RuntimeOptions -> "Speed"];

(*calulated fitted points g(fr1, fr2, dc, pdc1, pdc2)*)
FunceC3 = Compile[{{fr1, _Real, 1}, {fr2, _Real, 1}, {dc, _Real, 1}, {pdc1, _Real, 1}, {pdc2, _Real, 1}, {bm, _Real, 1}},
   Block[{fr1e = Exp[fr1], fr2e = Exp[fr2]},
    Chop[Transpose[Map[((
          (Exp[Exp[pdc1] #] fr1e)/(1 + fr1e) +
           (Exp[Exp[pdc2] #] fr2e)/(1 + fr2e) -
           (Exp[Exp[dc] #] (-1 + fr1e fr2e))/((1 + fr1e) (1 + fr2e))
          )) &, -bm]]]
    ], Parallelization -> True, RuntimeOptions -> "Speed"];

(*calculate probability*)
DotC = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}}, 
	((vec1.vec2)^2)/(vec2.vec2),
   RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> "Speed"];
Dotc1 = Compile[{{vec, _Real, 1}}, 
	vec.vec,
   RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> "Speed"];
MatDot2 = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}, {mat, _Real, 2}}, 
	(vec1.mat.vec1) - (vec2.mat.vec2),
   RuntimeAttributes -> {Listable}, Parallelization -> True, RuntimeOptions -> "Speed"];
AlphaC = Compile[{
	{theta, _Real, 2}, {thetat, _Real, 2}, {mu, _Real, 1},
	{icov, _Real, 2}, {y, _Real, 2}, {yty, _Real, 1}, {g, _Real, 2}, 
	{gt, _Real, 2}, {nb, _Real, 0}, {nvox, _Real, 0}}, 
   Block[{gttgt, gtg, ytg, ytgt, pt, pd, alpha,pdpt,rand,bool,top,bot},
    (*probability 1*)
    pt = Exp[0.5 (MatDot2[Transpose[theta - mu], Transpose[thetat - mu], icov])];
    (*probability 2*)
    pd = Chop[((yty - DotC[y, gt])/(yty - DotC[y, g]))^(-nb/2)];
    (*bool=alpha-RU*)
    rand = RandomReal[1, nvox];
    pdpt=(pd*pt);
    bool = pdpt - rand;
    UnitStep[bool]
    ],
   Parallelization -> True, RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*IVIM Plot Functions*)


(* ::Subsubsection::Closed:: *)
(*HistogramPar*)


SyntaxInformation[HistogramPar] = {"ArgumentsPattern" -> {_, _, _, _, _.}};

HistogramPar[dat_, {con_, bin_}, sel_, col_, ran_: .5] :=Block[{mu,cov},
	{mu, cov} = MeanCov[dat]; 
	HistogramPar[dat, {con, bin, mu, cov}, sel, col, ran]
	];

HistogramPar[dat_, {con_, bin_, mu_, cov_}, sel_, col_, ran_: .5] := 
 Module[{label, data, ticks, tickst, tickste, len, ss, hist, pdf, binsize,x},
  If[Length[dat] != Length[con], Return[]];
  
  binsize = (-Subtract @@ #/bin) & /@ con;
  
  tickst = {{0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}, {0.1, 0.2, 0.5,
      1, 2, 3}, {0.5, 1, 2, 5, 10, 20, 50, 100, 500}};
  tickste = ThetaConvi[{1., .001, .001} tickst];
  
  data = If[sel == 1, 
  	DeleteCases[Flatten[#], 0.] & /@ dat, 
    DeleteCases[Flatten[#], 0.(*Indeterminate*)] & /@ dat
    ];
  
  label = {{"F", "D", "pD"}, {"f", "d", "pd"}, {"F  (no units)", 
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
     " min   -   each step takes: ", Round[t1/nit, .001], 
     " s   -   full chain (21000) takes: ", 
     Round[(21000/nit) (t1/60), .1], " min"}]
   ,
   GraphicsRow[MapThread[Show[
       ListPlot[#1, AspectRatio -> .1 Length[wstart], 
        PlotRange -> All, PlotStyle -> {Black}],
       Plot[#2, {x, 0, nvox}, PlotStyle -> {Red, Thick}]
       ] &, {w, wstart}], ImageSize -> 1000]
   ,
   Show[
    ListPlot[t2s, AspectRatio -> 0.075, ImageSize -> 1000, 
     PlotStyle -> {Black, PointSize[Medium]}],
    ListLinePlot[{
      GaussianFilter[t2s, 20, Padding -> "Reversed"],
      GaussianFilter[t2s, Length[t2s], Padding -> "Reversed"]
      }, AspectRatio -> 0.075, ImageSize -> 1000, 
     PlotStyle -> {Directive[{Red, Thick}], 
       Directive[{Red, Dashed, Thick}]}, PlotRange -> Full
     ]]
   }]


(* ::Subsubsection::Closed:: *)
(*UncertainPlot*)


SyntaxInformation[UncertainPlot] = {"ArgumentsPattern" -> {_, _, _, _.}};

UncertainPlot[mn_, sig_, con_, ran_:.1] := 
 Module[{tickst, tickste, label, ticks, len, ss},
  tickst = {{0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}, {0.1, 0.2, 0.5,
      1, 2, 3}, {0.5, 1, 2, 5, 10, 20, 50, 100, 500}};
  tickste = ThetaConvi[{1., .001, .001} tickst];
  
  label = {{"F  (no units)", 
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



(* ::Subsubsection::Closed:: *)
(*NNLeastSquares*)


SyntaxInformation[NNLeastSquares] = {"ArgumentsPattern" -> {_, _}};

(*Main function*)
NNLeastSquares[A_, f_] := With[{
    toIndicesZ = Compile[{{v, _Real, 1}}, Flatten[Position[v, 1.]]],
    toIndicesP = 
     Compile[{{v, _Real, 1}}, Flatten[Position[1 - v, 1.]]],
    wCalc = 
     Compile[{{AA, _Real, 2}, {ff, _Real, 1}, {x, _Real, 1}}, 
      Transpose[AA].(ff - AA.x)],
    zCalc = 
     Compile[{{R, _Real, 2}, {Q, _Real, 2}, {ff, _Real, 1}}, 
      Inverse[R].Q.ff],
    xCalc = 
     Compile[{{x, _Real, 1}, {z, _Real, 1}, {posSet, _Integer, 1}},
      Block[{negz, alpha},
       negz = Select[posSet, z[[#]] < 0. &];
       alpha = Min[x[[negz]]/(x[[negz]] - z[[negz]])];
       N@Chop[x + alpha (z - x), 10.^-10]]]
    },
   Quiet[
    Block[{zerSet, posSet, maxPos, setz, z, x, zeroed, w, zeros, i, j,
       l},
     i = j = 1;
     zeros = 0 A[[1]];
     l = Length[zeros];
     (*fuctions to get posSet and zerSet and to calc z*)
     zerSet := toIndicesZ[zeroed];
     posSet := toIndicesP[zeroed];
     maxPos := Last[Ordering[zeroed w]];
     setz := Block[{Q, R},
       {Q, R} = QRDecomposition[A[[All, posSet]]];
       z = zeros; z[[posSet]] = zCalc[R, Q, f];
       ];
     (*initialize x, w, posSet and zerSet*)
     x = zeros;
     w = wCalc[A, f, x];
     zeroed = zeros + 1.;
     (*fist While loop*)
     While[j < l && zerSet != {} && Max[w[[zerSet]]] > 10.^-10,
      j++;
      (*index of max(w) and add j to posSet*)
      zeroed[[maxPos]] = 0.;
      (*calc z values for posSet*)
      setz;
      (*second While loop*)
      While[Min[z] < 0.,
       i++;
       (*calc alpha and recalculate x*)
       x = xCalc[x, z, posSet];
       (*move to zerSet all posSet vals were x=0*)
       zeroed[[Select[posSet, x[[#]] == 0. &]]] = 1.;
       (*calc z values for posSet*)
       setz;
       ];
      (*recalculate x and w*)
      w = wCalc[A, f, x = z];
      ];
     (*output solution*)
     Return[x];
     ]]];

(*support functions*)
toIndicesZ = Compile[{{v, _Real, 1}}, Flatten[Position[v, 1.]]];
toIndicesP = Compile[{{v, _Real, 1}}, Flatten[Position[1 - v, 1.]]];
wCalc = Compile[{{A, _Real, 2}, {f, _Real, 1}, {x, _Real, 1}}, 
   Transpose[A].(f - A.x)];
zCalc = Compile[{{R, _Real, 2}, {Q, _Real, 2}, {f, _Real, 1}}, 
   Inverse[R].Q.f];
xCalc = Compile[{{x, _Real, 1}, {z, _Real, 1}, {posSet, _Integer, 1}},
   Block[{negz, alpha},
    negz = Select[posSet, z[[#]] < 0 &];
    alpha = Min[x[[negz]]/(x[[negz]] - z[[negz]])];
    N@Chop[x + alpha (z - x), 10.^-10]]];


(* ::Subsubsection::Closed:: *)
(*IVIMCorrectData*)


Options[IVIMCorrectData] = {FilterMaps -> True, FilterType -> "Median", FilterSize -> 1};

SyntaxInformation[IVIMCorrectData] = {"ArgumentsPattern" -> {_, {_, _,_} , _ , OptionsPattern[]}};

IVIMCorrectData[data_, {S0_, f_, pdc_}, bval_, OptionsPattern[]] := Module[{ff, pdcf, filt, dataSyn, dataCor},
    
  {ff, pdcf} = If[OptionValue[FilterMaps],
    filt = Switch[OptionValue[FilterType],
    	"Median", MedianFilter,
    	_, GaussianFilter
    	];
    {filt[f, OptionValue[FilterSize]],filt[pdc, OptionValue[FilterSize]]}
    ,
    {f, pdc}
    ];
  
  dataSyn = Round@Clip[SynDatai[S0, ff, pdcf, bval], {0, Infinity}];
  dataCor = Clip[data - dataSyn, {0, 1.1 Max[data]}];
  
  {dataCor, dataSyn}
  ]

SynDatai = Compile[{{s0, _Real, 3}, {f, _Real, 3}, {pdc, _Real, 3}, {bval, _Real, 1}}, Transpose[Map[(f s0 Exp[-# pdc]) &, bval]]];


(* ::Subsubsection::Closed:: *)
(*IVIMResiduals*)


IVIMResiduals[data_, binp_, pars_] := 
 Module[{depthD,depthP,dat,par,res},
  
  (*data checks*)
  depthD = ArrayDepth[data];
  depthP = ArrayDepth[pars];
  
  dat = N[If[depthD == 4, Transpose[data], data]];
  dat = If[depthD > 1, TransData[dat, "l"], dat];
  par = If[depthP > 1, TransData[pars, "l"], pars];
  
  res = IVIMResCalcC[dat, binp, par];
  
  res = If[depthD > 1, TransData[res, "r"], res];
  
  Sqrt[Mean[Drop[res, 1]^2]] // N
  
  ];
  
  IVIMResCalcC = 
  Block[{S0, f1, f2, dc, pdc1, pdc2, out}, 
   Compile[{{dat, _Real, 1}, {binp, _Real, 1}, {pars, _Real, 1}},
    out = Switch[Length[pars],
      2,
      {S0, dc} = pars;
      (S0*(((Exp[-binp dc])))),
      4,
      {S0, f1, dc, pdc1} = pars;
      (S0*((((1 - f1)*Exp[-binp dc]) + (f1*Exp[-binp pdc1])))),
      6,
      {S0, f1, f2, dc, pdc1, pdc2} = pars;
      (S0*((((1 - f1 - f2)*Exp[-binp dc]) + (f1*
             Exp[-binp pdc1]) + (f2*Exp[-binp pdc2]))))
      ];
    dat - out
    , {{out, _Real, 1}}, RuntimeAttributes -> {Listable}, 
    RuntimeOptions -> "Speed"
    ]];



(* ::Section:: *)
(*End Package*)


End[]

SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`IVIMTools`*"];

EndPackage[]
