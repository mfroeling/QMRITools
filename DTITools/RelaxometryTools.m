(* ::Package:: *)

(* ::Title:: *)
(*DTITools RelaxometryTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`RelaxometryTools`", {"Developer`"}];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`RelaxometryTools`*"];
ClearAll @@ Names["DTITools`RelaxometryTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


T1rhoFit::usage = "
T1rhoFit[data, EchoTimes]"

T2Fit::usage = "
T2Fit[data, EchoTimes]"

TriExponentialT2Fit::usage = "
TriExponentialT2Fit[data, EchoTimes] fits the T2 based on Azzabou N et.al. Validation of a generic approach to muscle water 
T2 determination at 3T in fat-infiltrated skeletal muscle. J. Magn. Reson. 2015."

EPGT2Fit::usage = "
EPGT2Fit[data, EchoTimes, angle, relax] fits the T2 based on Marty B et.al. Simultaneous muscle water T2 and fat fraction mapping using transverse relaxometry with stimulated echo compensation.
angle is the refocussing angle in degree."

EPGSignal::usage = "
EPGSignal[Necho, echoSpace, T1, T2, angle, B1] generates a EPG T2 curve with stimulated echos. T1, T2 and echoSpace are in ms, angel is in degree, B1 is between 0 and 1."

CreateT2Dictionary::usage = "
CreateT2Dictionary[{T1m, T1f, T2f}, {Necho, echoSpace, angle}] Creates a EPG signal dictionary used for EPGT2fit."

DictionaryMinSearch::usage = "
DictionaryMinSearch[dictionary, y] performs dictionary minimization of data y. dictionary is generated with CreateT2Dictionary."


(* ::Subsection:: *)
(*General Options*)


DictT2Range::usage = "DictT2Range is an option for CreateT2Dictionary and EPGT2Fit. is specifies the range and step of the T2 values in the dictionary {min, max, step} in ms."

DictB1Range::usage = "DictB1Range is an option for CreateT2Dictionary and EPGT2Fit. It specifies the range and step of the B1 values in the dictionary {min, max, step}."

EPGRelaxPars::usage = "EPGRelaxPars is and option for EPGT2Fit. Needs to be {T1muscl, T1Fat, T2Fat} in ms, defaul is {1400,365,137}."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*T1rhoFit*)


Options[T1rhoFit] = {Method -> "Linear"};

SyntaxInformation[T1rhoFit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

T1rhoFit[datan_, times_, OptionsPattern[]] := 
 Switch[OptionValue[Method],
  "Linear", LinFit[datan, times],
  _, LogFit[datan, times]
  ]


(* ::Subsection::Closed:: *)
(*T2Fit*)


Options[T2Fit] = {Method -> "Linear"};

SyntaxInformation[T2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

T2Fit[datan_, times_, OptionsPattern[]] := Switch[OptionValue[Method],
  "Linear", LinFit[N[datan], times],
  _, LogFit[N[datan], times]
  ]


(* ::Subsection:: *)
(*linFit and LogFit*)


(* ::Subsubsection::Closed:: *)
(*LinFit*)


LinFit[datan_, times_] := 
 Module[{datal, result, fdat, offset, T1r, t, ad,off,t1rho},
  ad = ArrayDepth[datan];
  datal = Log[datan] /. {Indeterminate -> 0, -Infinity -> 0} // N;
  datal = Switch[ad,
    3, Transpose[datal, {3, 1, 2}],
    4, Transpose[datal, {1, 4, 2, 3}]
    ];
  
  result = ParallelMap[If[Total[#]==0.,
  		{0.,0.},
  		fdat = Transpose[{times, #}];
  		(*Quiet[LinearModelFit[fdat, t, t]["BestFitParameters"]]*)
  		{off, t1rho} /. Quiet[FindFit[fdat, off + t1rho t, {off, t1rho}, t]]
  	]&, datal, {ad - 1}];
  
  {offset, T1r} = Switch[ad,
    3, Transpose[result, {2, 3, 1}],
    4, Transpose[result, {2, 3, 4, 1}]
    ];
  
  {offset, T1r} = {Exp[offset], -1/(T1r /. 0. -> Infinity)};
  T1r = Clip[T1r, {0, 500}, {0, 500}];
  {offset, T1r}
  ]


(* ::Subsubsection::Closed:: *)
(*LogFit*)


LogFit[datan_, times_] := 
 Module[{result, fdat, offset, T1r, off, t1rho, t, ad, datal},
  ad = ArrayDepth[datan];
  datal = Switch[ad,
    3, Transpose[datan, {3, 1, 2}],
    4, Transpose[datan, {1, 4, 2, 3}]
    ];
  
  result = ParallelMap[
    (
      fdat = Transpose[{times, #}];
      {off, t1rho} /. Quiet[NonlinearModelFit[
          fdat,
          {off Exp[-t/t1rho], 10 < t1rho < 301},
          {{off, Last@Last@fdat}, {t1rho, 100}},
          t]["BestFitParameters"]]
      ) &, datal, {ad - 1}];
  
  {offset, T1r} = Switch[ad,
    3, Transpose[result, {2, 3, 1}],
    4, Transpose[result, {2, 3, 4, 1}]
    ];
  
  T1r = Clip[T1r, {0, 500}, {0, 500}];
  {offset, T1r}
  ]


(* ::Subsection::Closed:: *)
(*EPGT2Fit*)


Options[EPGT2Fit]= {DictT2Range -> {20., 80., 0.3}, DictB1Range -> {0.4, 1., 0.02}, EPGRelaxPars->{1400., 365., 137.}}

SyntaxInformation[EPGT2Fit]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}}

EPGT2Fit[datan_, times_, angle_, opts:OptionsPattern[]]:=Block[{Necho,echoSpace,T1m, T1f, T2f, ad, datal, sol,
	wat, fat, fatMap, T2map, B1Map,
	dictf, valsf,ydat,fwf,residualError
	},
  
  {T1m, T1f, T2f} = OptionValue[EPGRelaxPars];
  
  Necho = Length[times];
  echoSpace = First[times];
 
  {dictf, valsf} = CreateT2Dictionary[{T1m, T1f, T2f}, {Necho, echoSpace, angle}, DictB1Range->OptionValue[DictB1Range], DictT2Range->OptionValue[DictT2Range]];
  
  ad = ArrayDepth[datan];
  datal = N[Switch[ad,
    3, Transpose[datan, {3, 1, 2}],
    4, Transpose[datan, {1, 4, 2, 3}]
    ]];
  
  PrintTemporary["starting dictionary min search"];  
  DistributeDefinitions[dictf, valsf,LeastSquaresC,ErrorC];
  sol = ParallelMap[(
  	ydat = N@#;
  	If[Total[ydat] == 0.,
  	(*skip if background*)
  	{{0., 0.}, {0., 0.}, 0.},
  	(*calcualte dictionary error*)
  	fwf = LeastSquaresC[dictf, ydat];
  	residualError = ErrorC[ydat, fwf, dictf];
  	(*find Min value*)
  	{valsf, fwf, residualError}[[All, First@Ordering[residualError, 1]]]
   ]
   ) &, datal, {ad - 1}];
  
  Switch[ad,
  	3,
  		{wat, fat} = Clip[{sol[[All, All, 2, 1]], sol[[All, All, 2, 2]]}, {0, Max[sol[[All, All, 2]]]}];
  		fatMap = Clip[sol[[All, All, 2, 2]]/((sol[[All, All, 2, 2]] + sol[[All, All, 2, 1]]) /. 0. -> Infinity), {0, 1}];
  		T2map = sol[[All, All, 1, 1]];
  		B1Map = MedianFilter[sol[[All, All, 1, 2]], 1];
  	,4,
		{wat, fat} = Clip[{sol[[All, All, All, 2, 1]], sol[[All, All, All, 2, 2]]}, {0, Max[sol[[All, All, All, 2]]]}];
  		fatMap = Clip[sol[[All, All, All, 2, 2]]/((sol[[All, All, All, 2, 2]] + sol[[All, All, All, 2, 1]]) /. 0. -> Infinity), {0, 1}];
  		T2map = sol[[All, All, All, 1, 1]];
  		B1Map = MedianFilter[sol[[All, All, All, 1, 2]], 1];
  ];
  
  {{T2map,B1Map},{wat, fat, fatMap}}
]

LeastSquaresC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]}, 
	(Inverse[T.A].T).y], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget->System`$DTIToolsCompiler];

ErrorC = Compile[{{y, _Real, 1}, {f, _Real, 1}, {A, _Real, 2}}, 
	Total[(y - A.f)^2], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget->System`$DTIToolsCompiler];


(* ::Subsection::Closed:: *)
(*DictionaryMinSearch*)


SyntaxInformation[DictionaryMinSearch]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

DictionaryMinSearch[{dictf_, valsf_}, yi_] := Block[{fwf, residualError, ydat},
  ydat = N@yi;
  If[Total[ydat] == 0.,
  	(*skip if background*)
  	{{0., 0.}, {0., 0.}, 0.},
  	(*calcualte dictionary error*)
  	fwf = LeastSquaresC[dictf, ydat];
  	residualError = ErrorC[ydat, fwf, dictf];
  	(*find Min value*)
  	{valsf, fwf, residualError}[[All, First@Ordering[residualError, 1]]]
   ]]



(* ::Subsection:: *)
(*EPGSignal*)


(* ::Subsubsection::Closed:: *)
(*EPGSignal*)


SyntaxInformation[EPGSignal]= {"ArgumentsPattern" -> {_, _, _, _, _, _, OptionsPattern[]}}

EPGSignal[Nechoi_, echoSpace_, T1_, T2_, angle_, B1_] := Block[
  {tau, T0, R0, alpha, Smat, Tmat, Rmat, Pmat, Emat, xvec, 
   out, t2r, t1r,Necho},
  (*define internal paramters*)
  Necho=Round[Nechoi];
  alpha = N[B1 angle Degree];
  tau = echoSpace/2.;
  t2r = Exp[-tau/T2];
  t1r = Exp[-tau/T1];
  
  (*Selection matrix to move all traverse states up one coherence level*)
  Smat = MixMatrix[Necho];
  (*RF mixing matrix*)
  T0 = RotMatrixT[alpha];
  Tmat = MakeDiagMat[T0, Necho, 1];
  (* Relaxation matrix*)
  R0 = DiagonalMatrix[N[{t2r, t2r, t1r}]];
  Rmat = MakeDiagMat[R0, Necho, t2r];
  (*Precession and relaxation matrix*)
  Pmat = Rmat.Smat;
  (*Matrix representing the inter-echo duration*)
  Emat = Pmat.Tmat.Pmat;
  (*Recursively apply matrix to get echo amplitudes*)
  xvec = ConstantArray[0, 3*Necho + 1]; xvec[[1]] = 1;
  (*create outpu*)
  out = Table[xvec = Emat.xvec; First[xvec], {i, Necho}];
  out/out[[1]]
  ]


(* ::Subsubsection::Closed:: *)
(*MakeDiagMat*)


MakeDiagMat[mat_, Necho_, first_] := Block[{out},
   out = ArrayPad[Flatten[Flatten[(DiagonalMatrix[ConstantArray[1, Necho]]*ConstantArray[mat, {Necho, Necho}]), {2, 4}], {2, 3}], {{1, 0}, {1, 0}}];
   out[[1, 1]] = first;
   out
   ];


(* ::Subsubsection::Closed:: *)
(*MixMatrix*)


(*if run once with Necho definition is stored*)
MixMatrix[Necho_] := MixMatrix[Necho] = Block[{len, Smat, off1, off2},
   len = 3*Necho + 1;
   (*mixing matirx*)
   Smat = ConstantArray[0, {len, len}];
   Smat[[1, 3]] = Smat[[2, 1]] = Smat[[3, 6]] = Smat[[4, 4]] = 1;
   Table[
    off1 = ((o - 1) - 1)*3 + 2;
    If[off1 <= len, Smat[[3*o - 1, off1]] = 1];
    off2 = ((o + 1) - 1)*3 + 3;
    If[off2 <= len, Smat[[3*o, off2]] = 1];
    Smat[[3*o + 1, 3*o + 1]] = 1;
    , {o, 2, Necho}];
   (*output mixing matrix*)
   Smat
   ]


(* ::Subsubsection::Closed:: *)
(*RotMatrixT*)


(*if run once with alpha definition is stored*)
RotMatrixT[alpha_] := RotMatrixT[alpha] = {
    {Cos[alpha/2]^2, Sin[alpha/2]^2, Sin[alpha]},
    {Sin[alpha/2]^2, Cos[alpha/2]^2, -Sin[alpha]},
    {-0.5 Sin[alpha], 0.5 Sin[alpha], Cos[alpha]}
    };


(* ::Subsection::Closed:: *)
(*CreateT2Dictionary*)


Options[CreateT2Dictionary] = {DictT2Range -> {20., 80., 0.3}, DictB1Range -> {0.4, 1., 0.02}};

SyntaxInformation[CreateT2Dictionary]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

CreateT2Dictionary[relax_, ang_, opts : OptionsPattern[]] := CreateT2Dictionaryi[N[relax], N[ang], N[OptionValue[DictT2Range]],N[OptionValue[DictB1Range]]]

(*save each unique dictionary*)
CreateT2Dictionaryi[relax_, ang_, T2range_, B1range_] := CreateT2Dictionaryi[relax, ang, T2range, B1range] = Block[
   {dict, dictf, valsf, T1m, T1f, T2f, Necho, echoSpace, angle, t2s, t2e, t2i, b1s, b1e, b1i},
   (*set parameters*)
   {T1m, T1f, T2f} = relax;
   {Necho, echoSpace, angle} = ang;
   (*get dictionary values*)
   {t2s, t2e, t2i} = T2range;
   {b1s, b1e, b1i} = B1range;
   Print["Creating dictionary for new values : ", {T1m, T1f, T2f}, "  ", {Necho, echoSpace, angle}, "  ", {t2s, t2e, t2i}, "  ", {b1s, b1e, b1i}];
   (*create dictionary*)
   dict = Table[
     {{EPGSignal[Necho, echoSpace, T1m, T2m, angle, B1], EPGSignal[Necho, echoSpace, T1f, T2f, angle, B1]}, {T2m, B1}}, 
     	{T2m, t2s, t2e, t2i}, {B1, b1s, b1e, b1i}];
   (*flatten dictionary*)  	
   valsf = Flatten[dict[[All, All, 2]], 1];
   dictf = Transpose /@ Flatten[dict[[All, All, 1]], 1];
   Print["The dictionary contains "<>ToString[Length[dictf]]<>" values."];
   {dictf, valsf}
   ]


(* ::Subsection::Closed:: *)
(*TriExponentialT2Fit*)


SyntaxInformation[TriExponentialT2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

TriExponentialT2Fit[datan_, times_] := 
 Block[{result, fdat, offset, T1r, off, t1rho, t, ad, datal,
   model, cs, cf, T2s, T2f, Af, Am, T2m, x, Aff, Amm, T2mf, S0, ffr, 
   mfr, T2, model2,
   maskT2, dataT2, fmask, fitData,
   Afi, Ami, csi, T2mi, T2fi, T2si, sci
   },
  ad = ArrayDepth[datan];
  
  model = Af *(cs*Exp[-x/T2s] + (1 - cs)*Exp[-x/T2f]) + Am * Exp[-x/T2m];
  
  Switch[ad,
   3,(*single slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[{Mask[Mean[datan], {2}]}, MaskComponents -> 2, MaskClosing -> 1][[1]];
   dataT2 = maskT2 # & /@ datan;
   dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[-1]], {0.5}];
   fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ dataT2}];
   
   datal = Transpose[dataT2, {3, 1, 2}],
   4,(*mulit slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[Mask[Mean[Transpose[datan]], {2}], MaskComponents -> 2, MaskClosing -> 1]; 
   dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[All, -1]], {0.5}];
   fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
   
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ Transpose[dataT2]}];
   
   datal = Transpose[dataT2, {1, 4, 2, 3}]
   ];
  
  (*perform callibration fit*)
  {Afi, Ami, csi, T2mi, T2fi, T2si} = {Af, Am, cs, T2m, T2f, T2s} /. 
    FindFit[fitData,
     {model, {0 <= cs <= 1, T2m < T2f, T2f < T2s, Am >= 0, Af > Am}},
     {{Af, 1.25 0.8 fitData[[1, 2]]}, {Am, 1.25 0.2 fitData[[1, 2]]}, {cs, 0.33}, {T2m, 35}, {T2f, 81}, {T2s, 250}}, 
     x, Method -> "NMinimize"];
  sci = Afi + Ami;
  {Afi, Ami} = {Afi, Ami}/sci;
  
  (*Print the callibration results*)
  Print[Row[{Column[Round[{sci, 100 Afi, {100 Ami, T2mi}, {100 csi, T2si, T2fi}}, .1], Alignment -> Center],
     Show[
     	ListLinePlot[fitData, PlotStyle -> Directive[Thick, Red], PlotRange -> {{0, 150}, {0, sci}}, ImageSize -> 150],
     	Plot[sci model /. Thread[{Af, Am, cs, T2m, T2f, T2s} -> {Afi, Ami, csi, T2mi, T2fi, T2si}], {x, 0, 200}, PlotStyle -> Directive[{Black, Dashed}]]
      ]}]];
  
  (*define the fat values in the model*)
  model2 = Af *(csi*Exp[-x/T2si] + (1 - csi)*Exp[-x/T2fi]) + Am * Exp[-x/T2m];
  DistributeDefinitions[model2];
  
  result = ParallelMap[(
      If[N@Total[#] === 0.,
       {0., 0., 0., 0.},
       
       (*fit the muscle T2*)
       fdat = Transpose[{times, #}];
       {Aff, Amm, T2mf} = {Af, Am, T2m} /. Quiet[FindFit[fdat, model2, {{Af, 1.25 0.1 fdat[[1, 2]]}, {Am, 1.25 0.9 fdat[[1, 2]]}, {T2m, 35}}, x]];
       
       S0 = Aff + Amm;
       {Aff, Amm} = {Aff, Amm}/S0;
       {S0, Aff, Amm, T2mf}
       
       ]) &, datal, {ad - 1}];
  
  {S0, ffr, mfr, T2} = Switch[ad, 
  	3, Transpose[result, {2, 3, 1}], 
  	4, Transpose[result, {2, 3, 4, 1}]];
  
  T2 = Clip[T2, {0, 500}, {0, 500}];
  S0 = Clip[S0, {Min[datan], 1.5 Max[datan]}, {0, 0}];
  mfr = Clip[mfr, {0, 1}, {0, 1}];
  ffr = Clip[ffr, {0, 1}, {0, 1}];
  
  N@{S0, ffr, mfr, T2}
  ]


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`RelaxometryTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`RelaxometryTools`*"];

EndPackage[]
