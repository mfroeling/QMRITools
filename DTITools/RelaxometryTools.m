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

T1rhoFit::usage="T1rhoFit[data, EchoTimes]"

T2Fit::usage="T2Fit[data, EchoTimes]"

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


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`RelaxometryTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`RelaxometryTools`*"];

EndPackage[]
