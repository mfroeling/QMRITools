(* ::Package:: *)

(* ::Title:: *)
(*DTITools Tools*)


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
T2 determination at 3T in fat-infiltrated skeletal muscle. J. Magn. Reson. 2015"

(* ::Subsection:: *)
(*General Options*)




(* ::Subsection:: *)
(*Error Messages*)




(* ::Section:: *)
(*Functions*)


Begin["`Private`"]

(* ::Subsection::Closed:: *)
(*General local definitions*)


LinFit[datan_, times_] := 
 Module[{datal, result, fdat, offset, T1r, t, ad},
  ad = ArrayDepth[datan];
  datal = Log[datan] /. {Indeterminate -> 0, -Infinity -> 0} // N;
  datal = Switch[ad,
    3, Transpose[datal, {3, 1, 2}],
    4, Transpose[datal, {1, 4, 2, 3}]
    ];
  
  result = ParallelMap[If[Total[#]==0.,
  	{0.,0.},
  	fdat = Transpose[{times, #}];
  	Quiet[LinearModelFit[fdat, t, t]["BestFitParameters"]]
  	]
      &, datal, {ad - 1}];
  
  {offset, T1r} = Switch[ad,
    3, Transpose[result, {2, 3, 1}],
    4, Transpose[result, {2, 3, 4, 1}]
    ];
  
  {offset, T1r} = {Exp[offset], -1/(T1r /. 0. -> Infinity)};
  T1r = Clip[T1r, {0, 500}, {0, 500}];
  {offset, T1r}
  ]

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


(* ::Subsection::Closed:: *)
(*TriExponentialT2Fit*)

SyntaxInformation[TriExponentialT2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

TriExponentialT2Fit[datan_, times_] := 
 Block[{result, fdat, offset, T1r, off, t1rho, t, ad, datal,
   model, cs, cf, T2s, T2f, Af, Am, T2m, x, Aff, Amm, T2mf, S0, ffr, 
   mfr, T2,
   maskT2, dataT2, fmask, fitData,
   Afi, Ami, csi, T2mi, T2fi, T2si, sci
   },
  ad = ArrayDepth[datan];
  
  model = 
   Af *(cs*Exp[-x/T2s] + (1 - cs)*Exp[-x/T2f]) + Am * Exp[-x/T2m];
  
  Switch[ad,
   3,(*single slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = 
    SmoothMask[{Mask[Mean[datan], {2}]}, MaskComponents -> 2, 
      MaskClosing -> 1][[1]];
   dataT2 = maskT2 # & /@ datan;
   dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[-1]], {0.5}];
   fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = 
    Transpose[{times, 
      Mean[Flatten[GetMaskData[#, fmask]]] & /@ dataT2}];
   
   datal = Transpose[dataT2, {3, 1, 2}],
   4,(*mulit slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = 
    SmoothMask[Mask[Mean[Transpose[datan]], {2}], MaskComponents -> 2,
      MaskClosing -> 1]; 
   dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[All, -1]], {0.5}];
   fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
   
   (*data for calibration fit*)
   fitData = 
    Transpose[{times, 
      Mean[Flatten[GetMaskData[#, fmask]]] & /@ Transpose[dataT2]}];
   
   datal = Transpose[dataT2, {1, 4, 2, 3}]
   ];
  
  (*perform callibration fit*)
  {Afi, Ami, csi, T2mi, T2fi, T2si} = {Af, Am, cs, T2m, T2f, T2s} /. 
    FindFit[fitData,
     {model, {0 <= cs <= 1, T2m < T2f, T2f < T2s, Am >= 0, Af > Am}},
     {{Af, 1.25 0.8 fitData[[1, 2]]}, {Am, 
       1.25 0.2 fitData[[1, 2]]}, {cs, 0.33}, {T2m, 35}, {T2f, 
       81}, {T2s, 250}}, x, Method -> "NMinimize"];
  sci = Afi + Ami;
  {Afi, Ami} = {Afi, Ami}/sci;
  
  (*Print the callibration results*)
  Print[Row[{
     Column[
      Round[{sci, 
        100 Afi, {100 Ami, T2mi}, {100 csi, T2si, T2fi}}, .1], 
      Alignment -> Center],
     Show[
      ListLinePlot[fitData, PlotStyle -> Directive[Thick, Red], 
       PlotRange -> {{0, 150}, {0, sci}}, ImageSize -> 150],
      Plot[
       sci model /. 
        Thread[{Af, Am, cs, T2m, T2f, T2s} -> {Afi, Ami, csi, T2mi, 
           T2fi, T2si}], {x, 0, 200}, 
       PlotStyle -> Directive[{Black, Dashed}]]
      ]}]];
  
  (*define the fat values in the model*)
  cs = csi; T2s = T2si; T2f = T2fi;
  
  result = ParallelMap[(
      If[N@Total[#] === 0.,
       {0., 0., 0., 0.},
       
       fdat = Transpose[{times, #}];
       {Aff, Amm, T2mf} = Quiet[{Af, Am, T2m} /. FindFit[fdat, model,
           {{Af, sci Afi}, {Am, sci Ami}, {T2m, T2mi}}, x]];
       
       S0 = Aff + Amm;
       {Aff, Amm} = {Aff, Amm}/S0;
       {S0, Aff, Amm, T2mf}
       
       ]) &, datal, {ad - 1}];
  
  Print[Dimensions[result]];
  
  {S0, ffr, mfr, T2} = 
   Switch[ad, 3, Transpose[result, {2, 3, 1}], 4, 
    Transpose[result, {2, 3, 4, 1}]];
  
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
