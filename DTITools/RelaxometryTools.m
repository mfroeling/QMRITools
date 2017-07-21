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


T1rhoFit::usage = 
"T1rhoFit[data, EchoTimes] fits the T1rho value to the data using linear or nonlinear methdos.

Output is {S(0), T1rhomap}."

T2Fit::usage = 
"T2Fit[data, EchoTimes] fits the T1rho value to the data using linear or nonlinear methdos.

Output is {S(0), T2}."

TriExponentialT2Fit::usage = 
"TriExponentialT2Fit[data, EchoTimes] fits the T2 based on Azzabou N et.al. Validation of a generic approach to muscle water T2 determination at 3T in fat-infiltrated skeletal muscle. J. Magn. Reson. 2015.

The fat T2 parameters are automatically estimated from the high signal voxels from the last echo.

Output is {{S(0), fatFraction, muscleFraction, T2map},callibration} or {S(0), fatFraction, muscleFranction, T2map}."

EPGSignal::usage = 
"EPGSignal[Necho, echoSpace, T1, T2, angle, B1] generates a EPG T2 curve with stimulated echos. 
T1, T2 and echoSpace are in ms, angel is in degree, B1 is between 0 and 1.

Output is the EPG Signal vector."

EPGT2Fit::usage = 
"EPGT2Fit[data, {Necho, detlaTE}, {exitation, refoucs}] fits the T2 based on Marty B et.al. Simultaneous muscle water T2 and fat fraction mapping using transverse relaxometry with stimulated echo compensation.

exitation and refocus are the RF pulse angles e.g. 90,180. They can also be a range of angeles over the slice profile as defined by GetSliceProfile.

Output is {{{T2map,B1Map},{wat, fat, fatMap}},callibration} or {{T2map,B1Map},{wat, fat, fatMap}}"

CreateT2Dictionary::usage = 
"CreateT2Dictionary[{T1m, T1f, T2f}, {Necho, echoSpace, angle}] Creates a EPG signal dictionary used for EPGT2fit.
Every dictionary that is defined is cached.

Output is {dictionary, vals}"

DictionaryMinSearch::usage = 
"DictionaryMinSearch[dictionary, y] performs dictionary minimization of data y. dictionary is generated with CreateT2Dictionary.

Output is {{T2, B1}, fwfraction, residualError}."

NonLinearEPGFit::usage = 
"NonLinearEPGFit[{vals, T2cons}, y] performs dictionary minimization of data y. 

vals = {{T1muscle, T1fat, T2fat}, {Necho, echoSpace, angle}}.

Output is {{T2, B1}, fwfraction, residualError}."

CalibrateEPGT2Fit::usage = 
"CalibrateEPGT2Fit[datan, times, angle] calculates the Fat T2 ralaxation that will be used in the EPGT2fit.

Outputs the fat T2 value."


(* ::Subsection:: *)
(*General Options*)


DictT2Range::usage = 
"DictT2Range is an option for CreateT2Dictionary and EPGT2Fit. is specifies the range and step of the T2 values in the dictionary {min, max, step} in ms."

DictB1Range::usage = 
"DictB1Range is an option for CreateT2Dictionary and EPGT2Fit. It specifies the range and step of the B1 values in the dictionary {min, max, step}."

EPGRelaxPars::usage = 
"EPGRelaxPars is and option for EPGT2Fit. Needs to be {T1muscl, T1Fat, T2Fat} in ms, defaul is {1400,365,137}."

MonitorEPGFit::usage = 
"MonitorEPGFit show waitbar during EPGT2Fit."

EPGFitPoints::usage = 
"EPGFitPoints is a option for CalibrateEPGT2Fit and EPGT2Fit. Number of points is 200 by default."

EPGCalibrate::usage = 
"EPGCalibrate is an option for EPGT2Fit. If set to True it does autmatic callibration of the T2 fat relaxation time."

OutputCalibration::usage = 
"OutputCalibration is an option for EPGT2Fit and TriExponentialT2Fit. If true it outputs the calibartion values."


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
 Block[{datal, result, fdat, offset, T1r, t, ad,off,t1rho},
  ad = ArrayDepth[datan];
  datal = LogNoZero[datan];
  datal = Switch[ad,
    3, Transpose[datal, {3, 1, 2}],
    4, Transpose[datal, {1, 4, 2, 3}]
    ];
    
  PrintTemporary["performing linear T2 fit"];
  
  result = ParallelMap[(
  	If[Total[#]==0.,
  		{0.,0.},
  		{off, t1rho} /. Quiet[FindFit[Transpose[{times, #}], off + t1rho t, {off, t1rho}, t]]]
  		)&, datal, {ad - 1}];
  
  {offset, T1r} = TransData[result,"r"];
  
  {offset, T1r} = {Exp[offset], DevideNoZero[-1,T1r]};
  T1r = Clip[T1r, {0, 500}, {0, 500}];
  {offset, T1r}
  ]


(* ::Subsubsection::Closed:: *)
(*LogFit*)


LogFit[datan_, times_] := 
 Block[{result, fdat, offset, T1r, off, t1rho, t, ad, datal},
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


(* ::Subsection:: *)
(*EPGSignal*)


(* ::Subsubsection::Closed:: *)
(*EPGSignal*)


SyntaxInformation[EPGSignal] = {"ArgumentsPattern" -> {_, _, _, _}};

EPGSignal[{Nechoi_, echoSpace_}, {T1_, T2_}, {ex_, angle_}, B1_] := EPGSignali[{Nechoi, echoSpace}, {T1, T2}, {ex, angle}, B1]

EPGSignali[{Nechoi_, echoSpace_}, {T1_, T2_}, {ex_?ListQ, angle_?ListQ}, B1_] := Block[{sig},
  sig = Map[EPGSignali[{Nechoi, echoSpace}, {T1, T2}, #, B1] &, Transpose[{ex, angle}]];
  sig=Total@Join[sig,sig[[2;;]]];
  sig/sig[[1]]
  ]
  
EPGSignali[{Nechoi_, echoSpace_}, {T1_, T2_}, {ex_, angle_}, B1_] := Block[{
	tau, T0, R0, alpha, Smat, Tmat, Rmat, Pmat, Emat, xvec, out, t2r, t1r, Necho, Mzex
  },
 (*define internal paramters*)
 Necho = Round[Nechoi];
 Mzex = Sin[B1 ex (Pi/180.)];
 alpha = B1 angle (Pi/180.);
 tau = echoSpace/2;
 
 (*define relaxation*)
 t2r = N@Exp[-tau/T2];
 t1r = N@Exp[-tau/T1];
   
 (*Selection matrix to move all traverse states up one coherence level*)
 {Smat, xvec} = MixMatrix[Necho];
 (*RF mixing matrix*)
 Tmat = MakeDiagMat[RotMatrixT[alpha], Necho, 1.];
   (* Relaxation matrix*)
 Rmat = MakeDiagMat[DiagonalMatrix[{t2r, t2r, t1r}], Necho, t2r];
   
 (*Precession and relaxation matrix*)
 Pmat = Rmat.Smat;
 (*Matrix representing the inter-echo duration*)
 Emat = Pmat.Tmat.Pmat;
 
 (*create outpu*)
 out = NestList[Dot[Emat, #] &, xvec, Necho][[2 ;;, 1]];
 Mzex out
 ]


(* ::Subsubsection::Closed:: *)
(*MakeDiagMat*)


MakeDiagMat[mat_, Necho_, first_] := Block[{out},
  out = ArrayPad[ArrayFlatten[IdentityMatrix[Necho] ConstantArray[mat, {Necho, Necho}]], {{1, 0}, {1, 0}}];
  out[[1, 1]] = first;
  out
  ];


(* ::Subsubsection::Closed:: *)
(*MixMatrix*)


(*if run once with Necho definition is stored*)
MixMatrix[Necho_] := MixMatrix[Necho] = Block[{len, Smat, off1, off2},
      len = 3*Necho + 1;
      (*mixing matirx*)
      Smat = ConstantArray[0., {len, len}];
      Smat[[1, 3]] = Smat[[2, 1]] = Smat[[3, 6]] = Smat[[4, 4]] = 1.;
      Table[
        off1 = ((o - 1) - 1)*3 + 2.;
        If[off1 <= len, Smat[[3*o - 1, off1]] = 1.];
        off2 = ((o + 1) - 1)*3 + 3.;
        If[off2 <= len, Smat[[3*o, off2]] = 1.];
        Smat[[3*o + 1, 3*o + 1]] = 1.;
        , {o, 2, Necho}];
      (*output mixing matrix*)
      N@{Smat, Prepend[ ConstantArray[0, 3*Necho ], 1.]}
      ]


(* ::Subsubsection::Closed:: *)
(*RotMatrixT*)


RotMatrixT = Compile[{{alpha, _Real, 0}}, {
        {Cos[alpha/2]^2, Sin[alpha/2]^2, Sin[alpha]},
        {Sin[alpha/2]^2, Cos[alpha/2]^2, -Sin[alpha]},
        {-0.5 Sin[alpha], 0.5 Sin[alpha], Cos[alpha]}
        }];


(* ::Subsection:: *)
(*EPGT2Fit*)


(* ::Subsubsection::Closed:: *)
(*EPGT2Fit*)


Options[EPGT2Fit]= {
	DictB1Range -> {0.5, 1.5, 0.025}, DictT2Range -> {10, 70, 0.3}, 
	EPGRelaxPars -> {1400., 365., 150}, Method -> "dictionaryM", MonitorEPGFit -> True,
	EPGCalibrate -> True, EPGFitPoints -> 50, OutputCalibration -> False}

SyntaxInformation[EPGT2Fit]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}}

EPGT2Fit[datan_, echoi_, angle_, OptionsPattern[]]:=Block[{
	echo, T1m, T1f, T2f, ad, datal, sol, wat, fat, fatMap, T2map, B1Map, clip, fatc,
	B1i, T2i, T2s, B1s, soli, error, dictf, valsf, ydat, fwf, residualError, T2mc, B1c, cal,
	ran,dict,vals,cons,start, b1ran, b1step, t2ran, t2step, points
	},
	
	{T1m, T1f, T2f} = OptionValue[EPGRelaxPars];
	clip = OptionValue[DictT2Range][[1 ;; 2]];
	echo = If[Length[echoi]===2, echoi, {Length[echoi],First[echoi]}];
	
	(*clibrate the fat signal from the real data*)
	If[OptionValue[EPGCalibrate],
	  Print["Callibrating EPG fat T2."];
	  cal = {T2mc, T2f, B1c, fatc} = CalibrateEPGT2Fit[datan, echo, angle, EPGRelaxPars -> {clip, {50, 300}, {T1m, T1f}}, EPGFitPoints -> OptionValue[EPGFitPoints]];
	  T2f = N@Round[T2f,5];(*make whole ms such that it is more likely for same values*)
	  Print["EPG fat callibration:  ", T2f, " ms"];
	  ];
  
  	(*get data *)
  	ad = ArrayDepth[datan];
	datal = N@Switch[ad, 3, Transpose[datan, {3, 1, 2}], 4, Transpose[datan, {1, 4, 2, 3}]];

	(*monitor calculation*)
	i = j = 0;
	SetSharedVariable[i]; ParallelEvaluate[j = 0];
	If[OptionValue[MonitorEPGFit], PrintTemporary[ProgressIndicator[Dynamic[i], {0, Times @@ Dimensions[datan[[All, 1]]]}]]];
	
	(*switch between fitting algorithms*)
	sol = Switch[OptionValue[Method],
		"NLLS",
		(*define fit values*)
		valsf = {echo, {T1m, T1f, T2f}, angle};
		(*monitor calculation*)
		PrintTemporary["starting NLLS fitting: ", DateString[]];
		(*perform the fit using parallel kernels*)
		DistributeDefinitions[ErrorFunc, LeastSquaresC, EPGSignali, NonLinearEPGFiti, valsf, clip];
		ParallelMap[((*monitor calculation*)
		   j++; If[j > 100, i += j; j = 1;];
		   NonLinearEPGFiti[{valsf, clip}, #]
		   ) &, datal, {ad - 1}]
		   
		,"dictionaryM",
		b1ran=N[OptionValue[DictB1Range]];
		t2ran=N[OptionValue[DictT2Range]];
	  	b1step = (b1ran[[2]] - b1ran[[1]])/6;
	  	t2step = (t2ran[[2]] - t2ran[[1]])/6;
	  	
	  	(*create the dictionary*)
		{dictf, valsf} = CreateT2Dictionary[{T1m, T1f, T2f}, echo, angle, DictB1Range -> b1ran, DictT2Range -> t2ran];
		(*get the range and make the dictionary*)
		ran = Length[Range @@ b1ran];
		dict = Partition[dictf, ran];
		vals = Partition[valsf, ran];
		
		(*define the constraints and intial search points*)
		cons = Dimensions[vals][[;; -2]];
		points = Flatten[Table[{i, j},
    		{i, t2ran[[1]] + 0.5 t2step, t2ran[[2]]-t2step, t2step},
    		{j, b1ran[[1]] + 0.5 b1step, b1ran[[2]]-b1step, b1step}], 1];
		start = First@Position[vals, First@Nearest[Flatten[vals, 1], #, 1]] & /@ points;
		
		(*monitor calculation*)
		PrintTemporary["starting dictionary min search (Nmin): ", DateString[]];
		(*perform the fit using parallel kernels*)
		DistributeDefinitions[dict, vals, cons, start, LeastSquaresC, DictionaryMinSearchi];
		ParallelMap[((*monitor calculation*)
		   j++; If[j > 100, i += j; j = 1;];
		   DictionaryMinSearchi[{dict, vals}, #, {cons,start}]
		   ) &, datal, {ad - 1}]
		
		,_,
	  	(*create the dictionary*)
		b1ran=N[OptionValue[DictB1Range]];
		t2ran=N[OptionValue[DictT2Range]];
		{dictf, valsf} = CreateT2Dictionary[{T1m, T1f, T2f}, echo, angle, DictB1Range -> b1ran, DictT2Range -> t2ran];
		(*monitor calculation*)
		PrintTemporary["starting dictionary min search (Brute): ", DateString[]];
		(*perform the fit using parallel kernels*)
		DistributeDefinitions[dictf, valsf, LeastSquaresC, ErrorC, DictionaryMinSearchi];
		ParallelMap[((*monitor calculation*)
		   j++; If[j > 100, i += j; j = 1;];
		   DictionaryMinSearchi[{dictf, valsf}, #]
		   ) &, datal, {ad - 1}]
	];

	(*restructure fit solution*)
	sol = TransData[sol, "r"];
	{wat, fat} = Clip[{sol[[3]], sol[[4]]}, {0, Max[sol[[3 ;; 4]]]}];
	fatMap = Clip[fat/((wat + fat) /. 0. -> Infinity), {0, 1}];
	T2map = sol[[1]];
	(*B1Map = MedianFilter[sol[[2]], 1];*)
	B1Map = sol[[2]];
	error = sol[[5]];
	
	(*if neede also output callibaration*)
	If[OptionValue[OutputCalibration],
		{{{T2map, B1Map}, {wat, fat, fatMap},Sqrt[error]/echo[[1]]}, cal},
		{{T2map, B1Map}, {wat, fat, fatMap},Sqrt[error]/echo[[1]]}
	]
]


(* ::Subsubsection::Closed:: *)
(*NonLinearEPGFit*)


SyntaxInformation[NonLinearEPGFit]= {"ArgumentsPattern" -> {_, _}}

NonLinearEPGFit[{valsf_,cons_}, yi_] := NonLinearEPGFiti[{valsf, cons}, N[yi]]

NonLinearEPGFiti[{_,_}, {0. ..}] = {0., 0., 0., 0., 0.};

NonLinearEPGFiti[{valsf_, cons_}, ydat_] := Block[{
   fwf, residualError, soli, T1m, T1f, T2f, echo, angle, T2s, B1s
   },
  
	(*perform the fit*)
	{residualError, soli} = Quiet@FindMinimum[{
    	 ErrorFunc[ydat, T2i, B1i, valsf],
     	{0.4 <= B1i, B1i <= 1.6, cons[[1]] <= T2i, T2i <= cons[[2]]}
     	}, {{T2i, 35}, {B1i, 0.9}},
     	AccuracyGoal -> 5, PrecisionGoal -> 5, MaxIterations -> 25];
     {T2s, B1s} = soli[[All, 2]];

	(*get corresponding fat fractions*)
	{echo, {T1m, T1f, T2f}, angle} = valsf;
	fwf = LeastSquaresC[Transpose[{
    	EPGSignali[echo, {T1m, T2s}, angle, B1s],
    	EPGSignali[echo, {T1f, T2f}, angle, B1s]
    }], ydat];
	
	(*export paramters*)
	Flatten[{{T2s, B1s}, fwf, residualError/Total[fwf]^2}]
   ]


(* ::Subsubsection::Closed:: *)
(*DictionaryMinSearch*)


SyntaxInformation[DictionaryMinSearch]= {"ArgumentsPattern" -> {_, _, _.}}

(*brute force dictionary min search*)
DictionaryMinSearch[{dictf_, valsf_}, ydat_] := DictionaryMinSearchi[{dictf, valsf}, ydat]

DictionaryMinSearchi[_, {0. ..}] = {0., 0., 0., 0., 0.}

DictionaryMinSearchi[{dictf_, valsf_}, ydat_] := Block[{
	fwf, residualError
	},
	(*calcualte dictionary error*)
 	fwf = LeastSquaresC[dictf, ydat];
 	residualError = ErrorC[ydat, fwf, dictf];
 	(*find Min value*)
 	Flatten[{valsf, fwf, residualError}[[All, First@Ordering[residualError, 1]]]]
  ]

(*minimization dictionary min search*)  
DictionaryMinSearch[{dict_, vals_}, ydat_, {cons_,start_}] := DictionaryMinSearchi[{dict, vals}, ydat, {cons,start}]

DictionaryMinSearchi[_, {0. ..},_] = {0., 0., 0., 0., 0.}

DictionaryMinSearchi[{dict_, vals_}, ydat_, {{maxx_,maxy_},input_}] := Block[{
	err, coor, ErrorFuncDic
	},
	ErrorFuncDic[sig_, {x_Integer, y_Integer}] := LeastSquaresErrorC[dict[[x, y]], sig];
	(*minimize the dictionary value*)
	{err, coor} = Quiet@NMinimize[{ErrorFuncDic[ydat, {x, y}], 1 <= x <= maxx && 1 <= y <= maxy}, {x, y}, Integers, 
		Method -> {"NelderMead", "InitialPoints" -> input}
		,PrecisionGoal -> 6, AccuracyGoal -> 6, MaxIterations -> 50];
	(*find the min vals*)
	Flatten[{vals[[##]], LeastSquaresC[dict[[##]], ydat], err}] & @@ coor[[All, 2]]
	]


(* ::Subsubsection::Closed:: *)
(*ErrorFunc*)


ErrorFunc[y_, T2m_Real, B1_Real, vals_] := Quiet@Block[{
    sig, T1m, T1f, T2f, echo, angle
    },
   {echo, {T1m, T1f, T2f}, angle} = vals;
   sig = Transpose[{
      EPGSignali[echo, {T1m, T2m}, angle, B1],
      EPGSignali[echo, {T1f, T2f}, angle, B1]
      }];
   LeastSquaresErrorC[sig, y]]

ErrorFunc[y_, T2m_Real, T2f_Real, B1_Real, vals_] := Quiet@Block[{
    sig, T1m, T1f, echo, angle
    },
   {echo, {T1m, T1f}, angle} = vals;
   sig = Transpose[{
      EPGSignali[echo, {T1m, T2m}, angle, B1],
      EPGSignali[echo, {T1f, T2f}, angle, B1]
      }];
   LeastSquaresErrorC[sig, y]
   ]


(* ::Subsubsection::Closed:: *)
(*LeastSquaresC*)


LeastSquaresC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]},
    (Inverse[T.A].T).y
    ], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget -> System`$DTIToolsCompiler];


(* ::Subsubsection::Closed:: *)
(*LeastSquaresErrorC*)


LeastSquaresErrorC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]},
    Total[(y - A.(Inverse[T.A].T).y)^2]
    ], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget -> System`$DTIToolsCompiler];


(* ::Subsubsection::Closed:: *)
(*ErrorC*)


ErrorC = Compile[{{y, _Real, 1}, {f, _Real, 1}, {A, _Real, 2}}, 
   Total[((y - A.f))^2]
   , RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget -> System`$DTIToolsCompiler];


(* ::Subsubsection::Closed:: *)
(*CalibrateEPGT2Fit*)


Options[CalibrateEPGT2Fit] = {EPGRelaxPars -> {{0, 100}, {20, 300}, {1400., 365.}}, EPGFitPoints -> 50};

SyntaxInformation[CalibrateEPGT2Fit]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

CalibrateEPGT2Fit[datan_, echo_, angle_, OptionsPattern[]] := Block[{
	ad, Necho, echoSpace, maskT2, dataT2, fmask, fitData, step, T2mmin, T2mmax, T2fmin, T2fmax, 
	T1m, T1f, fatFr, watFr, cons, valsf, wcons, soli, fits, residualError
 	},
  
	ad = ArrayDepth[datan];
  
	(*Swtich between 3D and 4D data*)
	(*define thet fat mask and get the fat only signals*)
	Switch[ad, 3,
	  (*single slice*)(*make mask an normalize data to first echo*)
	  maskT2 = 
	   SmoothMask[{Mask[Mean[datan], {2}]}, MaskComponents -> 2, MaskClosing -> 1][[1]];
	  dataT2 = maskT2 # & /@ datan;
	  dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
	  (*create mask selecting fat*)
	  fmask = Mask[dataT2[[-1]], {0.5}];
	  fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
	  (*data for calibration fit*)
	  fitData = Transpose[Flatten[GetMaskData[#, fmask]] & /@ dataT2],
	  4,
	  (*mulit slice*)
	  (*make mask an normalize data to first echo*)
	  maskT2 = Mask[Mean[Transpose[datan]], {2}];
	  dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
	  
	  (*create mask selecting fat*)
	  fmask = Mask[dataT2[[All, -1]], {0.5}];
	  fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
	  
	  (*data for calibration fit*)
	  fitData = 
	   Transpose[Flatten[GetMaskData[#, fmask]] & /@ Transpose@dataT2]
	  ];
  
	(*select random fit points to calibrate fat signal and get the boundries*)
	step = Ceiling[Length[fitData]/OptionValue[EPGFitPoints]];
	{{T2mmin, T2mmax}, {T2fmin, T2fmax}, {T1m, T1f}} = OptionValue[EPGRelaxPars];

	(*define fit values*)
	valsf = {echo, {T1m, T1f}, angle};
	  
	(*perform the fit using parallel kernels*)
	DistributeDefinitions[ErrorFunc, LeastSquaresC, LeastSquaresErrorC, EPGSignali, valsf, echo, T2mmin, T2mmax, T2fmin, T2fmax];
	fits = ParallelMap[(
		(*calcualte NLLS error*)
		{residualError, soli} = Quiet@FindMinimum[{
	         ErrorFunc[#, T2mi, T2fi, B1i, valsf], 
	         {0.5 <= B1i, B1i <= 1.5, T2mmin <= T2mi, T2mi <= T2mmax, T2fmin <= T2fi, T2fi <= T2fmax}
	         }, {{T2mi, 10.}, {T2fi, 135.}, {B1i, 1}}, AccuracyGoal -> 5, PrecisionGoal -> 5, MaxIterations -> 25];
	     {T2mif, T2fif, B1if} = soli[[All, 2]];
	     (*Find the fat fraction*)
	     {watFr, fatFr} = LeastSquaresC[Transpose[{
	         EPGSignali[echo, {T1m, T2mif}, angle, B1if],
	         EPGSignali[echo, {T1f, T2fif}, angle, B1if]
	         }], #];
	     (*give the output*)
	     {T2mif, T2fif, B1if, fatFr/(watFr + fatFr)}
	     ) &, fitData[[1 ;; ;; step]]];
  
  Mean[fits]
  ]


(* ::Subsubsection::Closed:: *)
(*CreateT2Dictionary*)


Options[CreateT2Dictionary] = {DictB1Range -> {0.5, 1.5, 0.02}, DictT2Range -> {10, 70, 0.25}};

SyntaxInformation[CreateT2Dictionary]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

CreateT2Dictionary[relax_, echo_, ang_, OptionsPattern[]] := CreateT2Dictionaryi[relax, echo, ang, OptionValue[DictT2Range], OptionValue[DictB1Range]]

(*save each unique dictionary*)
CreateT2Dictionaryi[relax_, echo_, ang_, T2range_, B1range_] := CreateT2Dictionaryi[relax, echo, ang, T2range, B1range] = Block[
   {dict, dictf, valsf, T1m, T1f, T2f, Necho, echoSpace, angle, t2s, t2e, t2i, b1s, b1e, b1i},
	(*set parameters*)
	{T1m, T1f, T2f} = N@relax;
	(*get dictionary values*)
	{t2s, t2e, t2i} = N@T2range;
	{b1s, b1e, b1i} = N@B1range;

Print["Creating new dictionary"];
(*
Print["Creating new dictionary for new values :",
  "\n	{T1muscle, T1fat, T2fat} = ", {T1m, T1f, T2f},
  "\n	{Necho, echoSpace}       = ", echo,
  "\n	{T2min, T2max, T2step}   = ", {t2s, t2e, t2i},
  "\n	{B1min, B1max, B1step}   = ", {b1s, b1e, b1i}];
*)
	DistributeDefinitions[EPGSignali, echo, T1m, ang, T1f,T2f, t2s, t2e, t2i, b1s, b1e, b1i];
	dict = ParallelTable[{{
	      EPGSignali[echo, {T1m, T2m}, ang, B1],
	      EPGSignali[echo, {T1f, T2f}, ang, B1]
	      }, 
	      {T2m, B1}
	      }, {T2m, t2s, t2e, t2i}, {B1, b1s, b1e, b1i}];
	 
	 (*flatten dictionary*)
	 valsf = Flatten[dict[[All, All, 2]], 1];
	 dictf = Transpose /@ Flatten[dict[[All, All, 1]], 1];
	 
	 Print["The dictionary contains " <> ToString[Length[dictf]] <> 
	    " values."];
	 {dictf, valsf}
	]


(* ::Subsection::Closed:: *)
(*TriExponentialT2Fit*)


Options[TriExponentialT2Fit]={OutputCalibration->False}

SyntaxInformation[TriExponentialT2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

TriExponentialT2Fit[datan_, times_,OptionsPattern[]] := 
 Block[{result, fdat, offset, T1r, off, t1rho, t, ad, datal,
   model, cs, cf, T2s, T2f, Af, Am, T2m, x, Aff, Amm, T2mf, S0, ffr, 
   mfr, T2, model2,cal,
   maskT2, dataT2, fmask, fitData,
   Afi, Ami, csi, T2mi, T2fi, T2si, sci
   },
  
  ad = ArrayDepth[datan];
  
  Switch[ad,
   
   3,(*single slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[{Mask[Mean[datan], {2}]}, MaskComponents -> 2, MaskClosing -> 1][[1]];
   dataT2 = maskT2 # & /@ datan;
   dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[-1]], {0.4}];
   fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ dataT2}];
   datal = Transpose[dataT2, {3, 1, 2}],
   
   4,(*mulit slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[Mask[Mean[Transpose[datan]], {2}], MaskComponents -> 2, MaskClosing -> 1]; 
   dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[All, -1]], {0.4}];
   fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ Transpose[dataT2]}];
   datal = Transpose[dataT2, {1, 4, 2, 3}]
   ];
 
  model = Af *(cs*Exp[-x/T2s] + (1 - cs)*Exp[-x/T2f]) + Am * Exp[-x/T2m];
  
  (*perform callibration fit*)
 {Afi, Ami, csi, T2mi, T2fi, T2si} = {Af, Am, cs, T2m, T2f, T2s} /. 
   FindFit[fitData, {model, 
   	{0.1 <= cs <= 0.9, 30 < T2m < 50, T2m < T2f, T2f < T2s, 0 <= Am, Am < Af}},
   	{{Af, 1 fitData[[1, 2]]}, {Am, 0.25 fitData[[1, 2]]},{cs, 0.33}, {T2m, 35}, {T2f, 81}, {T2s, 250}}, 
   	x, Method -> "NMinimize"];
  (*normalize signal fractions*)  	
  sci = Afi + Ami;
  {Afi, Ami} = {Afi, Ami}/sci;
  cal = Round[{sci, {100 Ami, T2mi}, {100 csi, T2si, T2fi}}, .1];
  
  (*Print the callibration results*)
  Print[Row[{
    Column[Row /@ {
       {"signal: ", Round[sci, .1]},
       {"f mus: ", Round[100 Ami, .1] , " - T2 mus: ", 
        Round[T2mi, .1] },
       {""},
       {"f fat-slow: ", Round[100 csi, .1]}, 
       {" T2 slow: ", Round[T2si, .1], " - T2 fast: ", Round[T2fi, .1]}
       }, Alignment -> Center], 
    Show[
    	ListLinePlot[fitData, PlotStyle -> Directive[Thick, Red], PlotRange -> {{0, 150}, {0, sci}}, ImageSize -> 150], 
    	Plot[sci model /. Thread[{Af, Am, cs, T2m, T2f, T2s} -> {Afi, Ami, csi, T2mi, T2fi, T2si}], {x, 0, 200}, PlotStyle -> Directive[{Black, Dashed}]]
    	]
    	}, "   "]];
  
  (*define the fat values in the model*)
  model2 = Af *(csi*Exp[-x/T2si] + (1 - csi)*Exp[-x/T2fi]) + Am * Exp[-x/T2m];
  DistributeDefinitions[model2];
  
  (*perform the per voxel fit*)
  result = ParallelMap[(
      If[N@Total[#] === 0.,
       {0., 0., 0., 0.},
       
       (*fit the muscle T2*)
       fdat = Transpose[{times, #}];
       {Aff, Amm, T2mf} = {Af, Am, T2m} /. Quiet[FindFit[fdat, model2, {{Af, 0.125 fdat[[1, 2]]}, {Am, 1.125 fdat[[1, 2]]}, {T2m, 35}}, x]];
       
       S0 = Aff + Amm;
       {Aff, Amm} = {Aff, Amm}/S0;
       {S0, Aff, Amm, T2mf}
       
       ]) &, datal, {ad - 1}];
  
  (*generate output*)
  {S0, ffr, mfr, T2} = Switch[ad, 
  	3, Transpose[result, {2, 3, 1}], 
  	4, Transpose[result, {2, 3, 4, 1}]];
  	
  T2 = Clip[T2, {0, 500}, {0, 500}];
  S0 = Clip[S0, {Min[datan], 1.5 Max[datan]}, {0, 0}];
  mfr = Clip[mfr, {0, 1}, {0, 1}];
  ffr = Clip[ffr, {0, 1}, {0, 1}];
  
  (*in needed also output the calibration values*)
  If[OptionValue[OutputCalibration],{N@{S0, ffr, mfr, T2},cal},N@{S0, ffr, mfr, T2}]
  ]


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`RelaxometryTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`RelaxometryTools`*"];

EndPackage[]
