(* ::Package:: *)

(* ::Title:: *)
(*QMRITools RelaxometryTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`RelaxometryTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$QMRIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


T1rhoFit::usage = 
"T1rhoFit[data, EchoTimes] fits the T1rho value to the data using linear or nonlinear methdos.

Output is {S(0), T1rhomap}."

T2Fit::usage = 
"T2Fit[data, EchoTimes] fits the T1rho value to the data using linear or nonlinear methods.

Output is {S(0), T2}."

TriExponentialT2Fit::usage = 
"TriExponentialT2Fit[data, EchoTimes] fits the T2 based on Azzabou N et.al. Validation of a generic approach to muscle water T2 determination at 3T in fat-infiltrated skeletal muscle. J. Magn. Reson. 2015.
The fat T2 parameters are automatically estimated from the high signal voxels from the last echo.

Output is {{S(0), fatFraction, muscleFraction, T2map},callibration} or {S(0), fatFraction, muscleFranction, T2map}."


EPGSignal::usage = 
"EPGSignal[{Necho, echoSpace}, {T1, T2}, {ex_angle,ref_angle}, B1] generates a EPG T2 curve with stimulated echos. 
T1, T2 and echoSpace are in ms, angel is in degree, B1 is between 0 and 1.

Output is the EPG Signal vector."

EPGT2Fit::usage = 
"EPGT2Fit[data, {Necho, detlaTE}, {exitation, refoucs}] fits the T2 based on Marty B et.al. Simultaneous muscle water T2 and fat fraction mapping using transverse relaxometry with stimulated echo compensation.
Exitation and refocus are the RF pulse angles e.g. 90,180. They can also be a range of angeles over the slice profile as defined by GetSliceProfile.

Output is {{{T2map,B1Map},{wat, fat, fatMap}},callibration} or {{T2map,B1Map},{wat, fat, fatMap}}"

CalibrateEPGT2Fit::usage = 
"CalibrateEPGT2Fit[datan, times, angle] calculates the Fat T2 ralaxation that will be used in the EPGT2fit.

Outputs the fat T2 value."

CreateT2Dictionary::usage = 
"CreateT2Dictionary[{T1m, T1f, T2f}, {Necho, echoSpace, angle}] Creates a EPG signal dictionary used for EPGT2fit.
Every dictionary that is defined is cached.

Output is {dictionary, vals}"

NonLinearEPGFit::usage = 
"NonLinearEPGFit[{vals, T2cons}, y] performs dictionary minimization of data y. vals = {{T1muscle, T1fat, T2fat}, {Necho, echoSpace, angle}}.

Output is {{T2, B1}, fwfraction, residualError}."

DictionaryMinSearch::usage = 
"DictionaryMinSearch[dictionary, y] performs dictionary minimization of data y. dictionary is generated with CreateT2Dictionary.

Output is {{T2, B1}, fwfraction, residualError}."


(* ::Subsection::Closed:: *)
(*General Options*)


EPGMethod::usage =
"EPGMethod is an optionf for EPGT2Fit. Values can be \"NLLS\", \"dictionary\" or \"dictionaryM\"."

MonitorEPGFit::usage = 
"MonitorEPGFit show waitbar during EPGT2Fit."


EPGRelaxPars::usage = 
"EPGRelaxPars is and option for EPGT2Fit. Needs to be {T1muscl, T1Fat, T2Fat} in ms, defaul is {1400,365,137}."

EPGCalibrate::usage = 
"EPGCalibrate is an option for EPGT2Fit. If set to True it does autmatic callibration of the T2 fat relaxation time."

EPGSmoothB1::usage = 
"EPGSmoothB1 is an options for EPGT2Fit. If set to True the B1 map of the fit will be smoothed after which the minimization if perfomed again but with a fixed B1."

OutputCalibration::usage = 
"OutputCalibration is an option for EPGT2Fit and TriExponentialT2Fit. If true it outputs the calibartion values."

DictT2Range::usage = 
"DictT2Range is an option for CreateT2Dictionary and EPGT2Fit. is specifies the range and step of the T2 values in the dictionary {min, max, step} in ms."

DictB1Range::usage = 
"DictB1Range is an option for CreateT2Dictionary and EPGT2Fit. It specifies the range and step of the B1 values in the dictionary {min, max, step}."

DictT2fRange::usage = 
"DictT2fRange is an option for CreateT2Dictionary and EPGT2Fit. is specifies the range and step of the T2 fat values in the dictionary {min, max, step} in ms. 
If a single value is given this fixed value is used a long as EPGCalibrate is False."

EPGFitPoints::usage = 
"EPGFitPoints is a option for CalibrateEPGT2Fit and EPGT2Fit. Number of points is 200 by default."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Exponential Fitting*)


(* ::Subsubsection::Closed:: *)
(*T1rhoFit*)


Options[T1rhoFit] = {Method -> "Linear"};

SyntaxInformation[T1rhoFit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

T1rhoFit[datan_, times_, OptionsPattern[]] := 
 Switch[OptionValue[Method],
  "Linear", LinFit[datan, times],
  _, LogFit[datan, times]
  ]


(* ::Subsubsection::Closed:: *)
(*T2Fit*)


Options[T2Fit] = {Method -> "Linear"};

SyntaxInformation[T2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

T2Fit[datan_, times_, OptionsPattern[]] := Switch[OptionValue[Method],
  "Linear", LinFit[N[datan], times],
  _, LogFit[N[datan], times]
  ]


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
   maskT2 = Mask[Mean[datan], {2}, MaskSmoothing ->True, MaskComponents -> 2, MaskClosing -> 1];
   dataT2 = MaskData[datan, maskT2];
   dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[-1]], {0.4}];
   fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ dataT2}];
   datal = Transpose[dataT2, {3, 1, 2}],
   
   4,(*mulit slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = Mask[Mean[Transpose[datan]], {2}, MaskSmoothing->True, MaskComponents -> 2, MaskClosing -> 1]; 
   dataT2 = NormalizeData[MaskData[datan, maskT2]];
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


(* ::Subsection:: *)
(*EPGSignal*)


(* ::Subsubsection::Closed:: *)
(*EPGSignal*)


SyntaxInformation[EPGSignal] = {"ArgumentsPattern" -> {_, _, _, _, _.}};

EPGSignal[{Nechoi_, echoSpace_}, {T1_, T2_}, {ex_, ref_}, B1_, f_: 0] := EPGSignali[{Nechoi, echoSpace}, {T1, T2}, {ex, ref}, B1, f]

EPGSignali[{Nechoi_, echoSpace_}, {T1_, T2_}, {ex_?ListQ, ref_?ListQ}, B1_, f_:0.] := Block[{sig},
	sig = Map[EPGSignali[{Nechoi, echoSpace}, {T1, T2}, #, B1, f] &, Transpose[{ex, ref}]];
	sig = Mean@Join[sig, sig[[2 ;;]]]
  ]

EPGSignali[{Necho_, echoSpace_}, {T1_, T2_}, {exi_, refi_}, B1_, f_:0.] := Block[
	{tau, T0, R0, ex, ref, Smat, Tmat, Rmat, Rvec, svec, t2r, t1r, states, w, funRot,funMove},
	(*define internal paramters*)
	states = Round[If[Necho >= 10, Max[{Necho/2, 10}], Necho]];
	(*convert to Rad*)
	ex = N[B1 exi Degree];
	ref = N[B1 refi Degree];
	tau = echoSpace/2.;
	
	(*if use off ressonance then use complex matrix*)
	w = -tau 2. Pi f/1000.;
	{funRot,funMove}=If[w==0.,{RotMatrixT, MoveStates}, {RotMatrixTI, MoveStatesI}];
	
	(*Selection matrix to move all traverse states up one coherence Level*)
	Smat = MixMatrix[states];
	svec = Rvec = ConstantArray[0., Length[Smat]];
	(*define relaxation*)
	t2r = Chop[Exp[-tau/T2 + w I]];
	t1r = Chop[Exp[-tau/T1 - w I]];
	(*Relaxation matrix*)
	Rmat = MakeDiagMat[DiagonalMatrix[{t2r, t2r, t1r}], states];
	Rvec[[3]] = (1. - t1r);
	(*RF mixing matrix*)
	Tmat = MakeDiagMat[funRot[ref, 0], states];
	(*Create Initial state*)
	svec[[1 ;; 3]] = funRot[ex, 90].{0., 0., 1.};
	(*combined relax and gradient and create output*)
	Abs[funMove[Rmat, Rvec, Smat, Tmat, svec, Round@Necho][[2 ;;, 1]]]
  ]


(* ::Subsubsection::Closed:: *)
(*MakeDiagMat*)


MakeDiagMat[mat_, Necho_] := ArrayFlatten[IdentityMatrix[Necho] ConstantArray[mat, {Necho, Necho}]]


(* ::Subsubsection::Closed:: *)
(*MixMatrix*)


(*if run once with Necho definition is stored*)
MixMatrix[Necho_] := MixMatrix[Necho] = Block[{len, Smat, vec, off1, off2},
   (*mixing matirx*)
   len = 3*(Necho);
   Smat = ConstantArray[0, {len, len}];
   (*define state transitions*)
   Smat[[1, 5]] = 1;(*yi-1\[Rule]xi*)
   Smat[[len, len]] = 1;(*zn*)
   Table[
    Smat[[o - 1, o + 2]] = 1;(*yi\[Rule]yi-1*)
    Smat[[o + 1, o - 2]] = 1;(*xi\[Rule]xi+1*)
    Smat[[o, o]] = 1;(*zi\[Rule]zi*)
    , {o, 3, len - 3, 3}];
   Smat
   ]
   



(* ::Subsubsection::Closed:: *)
(*RotMatrixT*)


RotMatrixT[alpha_, ___] := RotMatrixTC[alpha];

(*using CPMG condition*)
RotMatrixTC = Compile[{{alpha, _Real, 0}}, Chop[{
     {Cos[alpha/2]^2, Sin[alpha/2]^2, Sin[alpha]},
     {Sin[alpha/2]^2, Cos[alpha/2]^2, -Sin[alpha]},
     {-0.5 Sin[alpha], 0.5 Sin[alpha], Cos[alpha]}
     }], RuntimeOptions -> "Speed"];

RotMatrixTI[alpha_, phi_: 90] := RotMatrixTCI[alpha, phi];

(*Specify angle and phase*)
RotMatrixTCI = Compile[{{alpha, _Real, 0}, {phi, _Real, 0}}, Chop[{
     {Cos[alpha/2]^2, Exp [2 phi I] Sin[alpha/2]^2, -I Exp [phi I] Sin[alpha]},
     {Exp [-2 phi I] Sin[alpha/2]^2, Cos[alpha/2]^2, I Exp [-phi I] Sin[alpha]},
     {-0.5 I Exp [-phi I] Sin[alpha], 0.5 I Exp [phi I] Sin[alpha], Cos[alpha]}
     }], RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MoveStates*)


MoveStates = Compile[{{Rmat, _Real, 2}, {Rvec, _Real, 1}, {Smat, _Real, 2}, {Tmat, _Real, 2}, {svec, _Real, 1}, {Necho, _Integer, 0}}, 
	(*Rmat = relaxation; Rvec = Mz recovery; Tmat = Rf pulse;*)
    (*1. Relaxation - 2. Mz-rec - 3. Change states - 4. RF pulse - 5. Relaxation - 6. Mz-rec - 7. Change states*)
    NestList[Chop[Smat.(Rmat.(Tmat.(Smat.((Rmat.#) + Rvec))) + Rvec)] &, svec, Necho]
    , RuntimeOptions -> "Speed"];

MoveStatesI = Compile[{{Rmat, _Complex, 2}, {Rvec, _Complex, 1}, {Smat, _Integer, 2}, {Tmat, _Complex, 2}, {svec, _Complex, 1}, {Necho, _Integer, 0}},
   (*Rmat = relaxation; Rvec = Mz recovery; Tmat = Rf pulse;*)
   (*1. Relaxation - 2. Mz-rec - 3. Change states - 4. RF pulse - 5. Relaxation - 6. Mz-rec - 7. Change states*)
   NestList[Chop[Smat.(Rmat.(Tmat.(Smat.((Rmat.#) + Rvec))) + Rvec)] &, svec, Necho]
   , RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*EPGT2Fit*)


(* ::Subsubsection::Closed:: *)
(*EPGT2Fit*)


Options[EPGT2Fit]= {
	DictB1Range -> {0.5, 1.4, 0.01}, DictT2Range -> {10., 60., 0.2}, DictT2fRange -> {120., 170., 1.}, 
	EPGRelaxPars -> {1400., 365.}, 
	EPGCalibrate -> False, EPGFitPoints -> 50, EPGMethod -> "dictionaryM", 
	MonitorEPGFit -> True, OutputCalibration -> False, EPGSmoothB1 -> True}

SyntaxInformation[EPGT2Fit]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}}

EPGT2Fit[datan_, echoi_, angle_, OptionsPattern[]]:=Block[{
	echo, T1m, T1f, T2f, ad, datal, sol, wat, fat, fatMap, T2map, B1Map, clip, fatc,
	B1i, T2i, T2s, B1s, soli, error, dictf, valsf, ydat, fwf, residualError, T2mc, B1c, cal,
	ran,dict,vals,cons,start, b1ran, b1step, t2ran, t2step, points, dim, size,
	b1rule, B1Int, dataf, b1, b1vals, b1len, t2rule, t2Int, t2, t2vals, t2len,
	sig, dictMat, dictfMat, t2fran, S0c, met, val, out, T2fmap, wMat
	},
	
	(*Get Input*)
	echo = If[Length[echoi]===2, echoi, {Length[echoi],First[echoi]}];
  	ad = ArrayDepth[datan];
	datal = N@Switch[ad, 1, datan, 3, Transpose[datan, {3, 1, 2}], 4, Transpose[datan, {1, 4, 2, 3}]];
	
	(*Get Options*)
	{T1m, T1f} = N@OptionValue[EPGRelaxPars];
	t2ran=N[OptionValue[DictT2Range]];
	b1ran=N[OptionValue[DictB1Range]];
	t2fran=N[OptionValue[DictT2fRange]];
	clip = t2ran[[1 ;; 2]];
	
	(*clibrate the fat signal from the real data*)
	If[OptionValue[EPGCalibrate]&&!VectorQ[datan],
	  Print["Callibrating EPG fat T2."];
	  cal = CalibrateEPGT2Fit[datan, echo, angle, EPGRelaxPars -> {clip, {50, 300}, {T1m, T1f}}, EPGFitPoints -> OptionValue[EPGFitPoints]];
	  (*{T2mc, T2f, B1c, fatc} = cal[[1]];*)
	  {T2f,B1c,S0c} = cal[[1]];
	  T2f = N@Round[T2f,5]+10; (* better to have to high T2fat than to Low *)
	  (*make whole ms such that it is more likely for same values*)
	  Print["EPG fat callibration:  ", T2f, " ms"];
	  ,
	  T2f = If[NumberQ[t2fran], Print["DictT2frange is not a number using 150ms"]; 150, t2fran];
	  ];
  
	(*monitor calculation*)
	i=0; SetSharedVariable[i]; ParallelEvaluate[j = 0];
	If[OptionValue[MonitorEPGFit]&&!VectorQ[datan], 
		dim = Times@@Dimensions[datal][[;;-2]];
		size = Round[dim/20];
		PrintTemporary[ProgressIndicator[Dynamic[i], {0, dim}]]];
	
	(*find the correct method*)
	met = If[OptionValue[EPGCalibrate] || NumberQ[t2fran], 
		PrintTemporary["Fitting using single Fat value: ",T2f]; val = 2; OptionValue[EPGMethod], 
		PrintTemporary["Fitting using dictionary of Fat values: ",T2f]; val = 3; "dictionaryM"
		];
	
	(*switch between fitting algorithms*)
	sol = Switch[met,
		
		"NLLS", (*non linear least sqyares*)
		(*define fit values*)
		valsf = {echo, {T1m, T1f, T2f}, angle};
		
		If[VectorQ[datal],
			(*single voxel*)
			NonLinearEPGFiti[{valsf, clip}, datal],
			(*monitor calculation*)
			PrintTemporary["Starting NLLS fitting: ", DateString[]];
			(*perform the fit using parallel kernels*)
			DistributeDefinitions[ErrorFunc, LeastSquaresC, NonLinearEPGFiti, EPGSignali, MixMatrix, MakeDiagMat, RotMatrixT, RotMatrixTI, MoveStates, MoveStatesI valsf, clip];
			ParallelMap[(
				j++; If[j > size, i += j; j = 1;]; 
				NonLinearEPGFiti[{valsf, clip}, #])&, datal, {ad - 1}]
		]
		   
		,"dictionaryM", (*NealderMead Nmnimize*)
	  	(*create the dictionary*)
		{dict, vals} = CreateT2Dictionaryi[{T1m, T1f}, echo, angle, t2ran, b1ran, T2f];
		(*extra weight on first two echos to get correct b1*)
		wMat = ConstantArray[1., echo[[1]]]; wMat[[1 ;; 2]] = 2 echo[[1]];wMat = DiagonalMatrix[wMat];
		dictMat = PseudoInverseWC[dict, wMat];
		cons = Dimensions[vals][[;; -2]];
		
		If[VectorQ[datal],
			(*single voxel*)
			DictionaryMinSearchi[{dict, dictMat, vals}, datal, cons],
			(*monitor calculation*)
			PrintTemporary["Starting dictionary min search (Nmin): ", DateString[]];
			(*perform the fit using parallel kernels*)
			DistributeDefinitions[dict, dictMat, vals, cons, size, LeastSquaresC, LeastSquaresError2C, DictionaryMinSearchi];
			ParallelMap[(
				j++; If[j > size, i += j; j = 1;]; 
				DictionaryMinSearchi[{dict, dictMat, vals}, #, cons])&, datal, {ad - 1}]
		]
		
		,_, (*Brute force min search*)
	  	(*create the dictionary*)
		{dict, vals} = CreateT2Dictionaryi[{T1m, T1f}, echo, angle, t2ran, b1ran, T2f];
		cons = Dimensions[vals][[;; -2]];
		dictf = Flatten[dict,1];
		dictfMat = PseudoInverseC[dictf];
		valsf = Flatten[vals,1];
		
		If[VectorQ[datal],
			(*single voxel*)
			DictionaryMinSearchi[{dict, dictfMat, vals}, datal],
			(*monitor calculation*)
			PrintTemporary["Starting dictionary min search (Brute): ", DateString[]];
			(*perform the fit using parallel kernels*)
			DistributeDefinitions[dictf, dictfMat, valsf, size, LeastSquaresC, LeastSquaresErrorC, ErrorC, DictionaryMinSearchi];
			ParallelMap[(
				j++; If[j > size, i += j; j = 1;]; 
				DictionaryMinSearchi[{dictf, dictfMat, valsf}, #])&, datal, {ad - 1}]
		]
	];

	(*restructure fit solution*)
	sol = If[VectorQ[datal], sol, TransData[sol, "r"]];
	
	(*perform B1/T2fat smoothing if needed*)
	B1Map = sol[[val]];
	If[val==3, T2fmap = sol[[2]]];
	
	(* check if B1map needs to be smoothed, only works for dictionary methods *)
	If[(OptionValue[EPGSmoothB1] && (OptionValue[EPGMethod]=!="NLLS") && !VectorQ[datal]),

		PrintTemporary["Starting B1 smoothing and refit with smooth B1: ", DateString[]];
					
		(*smooth the B1 map*)
		B1Map = Clip[N[Round[LapFilt[B1Map] - b1ran[[1]], b1ran[[3]]] + b1ran[[1]]], b1ran[[1;;2]], b1ran[[1;;2]]];
		(*convert B1 values to dictionaly integers*)
		b1vals = N@Range @@ b1ran; b1len = Length[b1vals];
		b1rule = Thread[b1vals -> Range[b1len]];
		B1Int = Clip[Round[B1Map /. b1rule ],{1,cons[[val]]}];
		
		(*monitor calculation*)
		i = 0; SetSharedVariable[i]; ParallelEvaluate[j = 0];
		DistributeDefinitions[dict, dictMat, vals, size, LeastSquaresC, ErrorC, DictionaryMinSearchi];
		
		(*perfomr the brute force fit of water only *)
		sol =If[val==2,
			(*definet fit data*)
			dataf = TransData[{datal, B1Int}, "l"];
			(*only fit T2 water*)
			 ParallelMap[(
				j++; If[j > size, i += j; j = 1;];
				sig=#[[1]]; b1=#[[2]];
				DictionaryMinSearchi[{dict[[All,b1]], dictMat[[All,b1]], vals[[All,b1]]}, sig]) &, dataf, {ad - 1}]
			,
			(*smooth the T2fmat*)
			T2fmap = Clip[N@Round[LapFilt[T2fmap] - t2fran[[1]], t2fran[[3]]] + t2fran[[1]], t2fran[[1 ;; 2]], t2fran[[1 ;; 2]]];
			t2vals = N@Range @@ t2fran; t2len = Length[t2vals];
			t2rule = Thread[t2vals -> Range[t2len]];
			t2Int = Clip[Round[T2fmap /. t2rule], {1, cons[[2]]}];
					
			(*definet fit data*)
			dataf = TransData[{datal, B1Int, t2Int}, "l"];
			(*only fit T2 water*)
			ParallelMap[(
				j++; If[j > size, i += j; j = 1;];
				sig=#[[1]];b1=#[[2]];t2=#[[3]];
				DictionaryMinSearchi[{dict[[All, t2, b1]], dictMat[[All, t2 ,b1]], vals[[All, t2, b1]]}, sig]) &, dataf, {ad - 1}]
		];	
		(*update the solution*)
		sol = TransData[sol,"r"];
	];
	
	(*get the outputs*)
	T2map = sol[[1]];
	{wat, fat} = N@Clip[{sol[[val+1]], sol[[val+2]]}, {0., Infinity}];
	fatMap = DevideNoZero[fat, (wat + fat)];
	error = Sqrt[sol[[val+3]]];
	
	(*if needed also output callibaration*)
	out = If[val==2,
		{{T2map, B1Map}, {wat, fat, fatMap}, error},
		{{T2map, T2fmap, B1Map}, {wat, fat, fatMap}, error}
	];
	
	If[OptionValue[OutputCalibration]&&OptionValue[EPGCalibrate], {out, cal[[1]]}, out]
]

LapFilt[data_, fil_:0.8] := Clip[Chop[ImageData[TotalVariationFilter[Image3D[N@data, "Real"], fil, Method -> "Laplacian", MaxIterations -> 15]]], MinMax[data]]


(* ::Subsubsection::Closed:: *)
(*CalibrateEPGT2Fit*)


Options[CalibrateEPGT2Fit] = {EPGRelaxPars -> {{0, 100}, {20, 300}, {1400., 365.}}, EPGFitPoints -> 50};

SyntaxInformation[CalibrateEPGT2Fit]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

CalibrateEPGT2Fit[datan_, echoi_, angle_, OptionsPattern[]] := Block[{
	ad, Necho, echoSpace, maskT2, dataT2, fmask, fitData, step, T2mmin, T2mmax, T2fmin, T2fmax, 
	T1m, T1f, fat, wat, cons, valsf, wcons, soli, fits, residualError, echo},
	
	echo = If[Length[echoi]===2, echoi, {Length[echoi],First[echoi]}];
	ad = ArrayDepth[datan];
  
	(*Swtich between 3D and 4D data*)
	(*define thet fat mask and get the fat only signals*)
	Switch[ad, 
		
	  3,
	  (*single slice*)
	  (*make mask an normalize data to first echo*)
	  maskT2 = Mask[Mean[datan]];
	  dataT2 = NormalizeData[maskT2 # & /@ datan];
	  (*create mask selecting fat*)
	  fmask = Mask[dataT2[[-1]], {50}];
	  fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
	  (*data for calibration fit*)
	  fitData = Transpose[Flatten[GetMaskData[#, fmask]] & /@ (dataT2+10.^-10)]-10.^-10;
	  
	  ,4,
	  (*mulit slice*)
	  (*make mask an normalize data to first echo*)
	  maskT2 = Mask[Mean[Transpose[datan]]];
	  dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
	   (*create mask selecting fat*)
	  fmask = Mask[dataT2[[All, -1]], {50}];
	  fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
	  (*data for calibration fit*)
	  fitData = Transpose[Flatten[GetMaskData[#, fmask]] & /@ Transpose[dataT2+10.^-10]]-10.^-10;
	  ];
  
	(*select random fit points to calibrate fat signal and get the boundries*)
	step = Ceiling[Length[fitData]/OptionValue[EPGFitPoints]];
	{{T2mmin, T2mmax}, {T2fmin, T2fmax}, {T1m, T1f}} = OptionValue[EPGRelaxPars];
	
	valsf = {echo, T1f, angle};
	
	DistributeDefinitions[ErrorFunc1, EPGSignali, MixMatrix, MakeDiagMat, RotMatrixT, RotMatrixTI, MoveStates, MoveStatesI, valsf];
	fits = ParallelMap[(
		S0 = #[[1]];
		{residualError, soli} = Quiet@FindMinimum[{
		    ErrorFunc1[#, T2fi, B1i, S0i, valsf],
		    {0.5 <= B1i <= 1.5, 20. <= T2fi <= 300., 0 < S0i}
		    }, {{T2fi, 50.}, {B1i, 1}, {S0i, 5 S0}}, MaxIterations -> 25];
		out = {T2fif, B1if, S0if} = soli[[All, 2]];
		out
		) &, fitData[[1 ;; ;; step]]];     
  
  	{Median[fits],StandardDeviation[fits]}
  ]


(* ::Subsubsection::Closed:: *)
(*CreateT2Dictionary*)


Options[CreateT2Dictionary] = {DictB1Range -> {0.5, 1.4, 0.01}, DictT2Range -> {10., 70., 0.2}, DictT2fRange -> {100., 200., 2.}};

SyntaxInformation[CreateT2Dictionary]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

CreateT2Dictionary[relax_, echo_, ang_, OptionsPattern[]] := CreateT2Dictionaryi[relax, echo, ang, OptionValue[DictT2Range], OptionValue[DictB1Range], OptionValue[DictT2fRange]]

(*save each unique dictionary*)
CreateT2Dictionaryi[relax_, echo_, ang_, t2range_, b1range_, t2frange_] := CreateT2Dictionaryi[relax, echo, ang, t2range, b1range, t2frange] = Block[{
	T1m, T1f, t2Mvals, t2Mlen, b1vals, b1len, t2val, t2Fvals, t2Flen, time, fatSig, watSig, dict, vals
	},
	
	(*set parameters*)
	{T1m, T1f} = N@relax;
	(*get dictionary values*)
	t2Mvals = N@Range @@ t2range;
	t2Mlen = Length[t2Mvals];
	b1vals = N@Range @@ b1range;
	b1len = Length[b1vals];
	(*check in number or range*)
	If[NumberQ[t2frange],
		t2val = t2frange;
		,
		t2Fvals = N@Range @@ t2frange;
		t2Flen = Length[t2Fvals];
	];
	
	DistributeDefinitions[EPGSignali, MixMatrix, MakeDiagMat, RotMatrixT, RotMatrixTI, MoveStates, MoveStatesI, echo, T1m, ang, T1f, b1vals, t2Mvals, t2Fvals, t2val];
	If[NumberQ[t2frange],
		PrintTemporary["Creating new dictionary with fixed T2 fat value"];
		(*fixed T2 value*)
		time = AbsoluteTiming[
			fatSig = ParallelTable[EPGSignali[echo, {T1f, t2val}, ang, B1], {B1, b1vals}];
			watSig = ParallelTable[EPGSignali[echo, {T1m, T2m}, ang, B1], {B1, b1vals}, {T2m, t2Mvals}];
			dict = Table[Transpose@{watSig[[b1i, t2mi]], fatSig[[b1i]]}, {t2mi, 1, t2Mlen}, {b1i, 1, b1len}];
			vals = Table[{t2m, b1}, {t2m, t2Mvals}, {b1, b1vals}];
		][[1]];
		
		Print["The dictionary contains "<>ToString[Times @@ Dimensions[dict][[;; -3]]]<>" values, with T2fat = "<>ToString[Round[t2val]]<>" ms, and took "<>ToString[Round[time, .1]]<>" seconds to generate."];
		,
		PrintTemporary["Creating new dictionary with range of T2 fat values"];
		(*range of T2 values*)
		time = AbsoluteTiming[
			fatSig = ParallelTable[EPGSignali[echo, {T1f, T2f}, ang, B1], {B1, b1vals}, {T2f, t2Fvals}];
			watSig = ParallelTable[EPGSignali[echo, {T1m, T2m}, ang, B1], {B1, b1vals}, {T2m, t2Mvals}];
			dict = Table[Transpose@{watSig[[b1i, t2mi]], fatSig[[b1i, t2fi]]}, {t2mi, 1, t2Mlen}, {t2fi, 1, t2Flen}, {b1i, 1, b1len}];
			vals = Table[{t2m, t2f, b1}, {t2m, t2Mvals}, {t2f, t2Fvals}, {b1, b1vals}];
		][[1]];
		
		Print["The dictionary contains "<>ToString[Times @@ Dimensions[dict][[;; -3]]]<>" values, and took "<> ToString[Round[time, .1]]<>" seconds to generate."];
	];
	(*output*)
	{dict,vals}
]


(* ::Subsubsection::Closed:: *)
(*NonLinearEPGFit*)


SyntaxInformation[NonLinearEPGFit]= {"ArgumentsPattern" -> {_, _}}

NonLinearEPGFit[{valsf_,cons_}, yi_] := NonLinearEPGFiti[{valsf, cons}, N[yi]]

NonLinearEPGFiti[{_,_}, {0. ..}] = {0., 0., 0., 0., 0.};

NonLinearEPGFiti[{valsf_, cons_}, ydat_] := Block[{fwf, residualError, soli, T1m, T1f, T2f, echo, angle, T2s, B1s},
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
	Flatten[{{T2s, B1s}, fwf, residualError}]
   ]


(* ::Subsubsection::Closed:: *)
(*DictionaryMinSearch*)


SyntaxInformation[DictionaryMinSearch]= {"ArgumentsPattern" -> {_, _, _.}}


(*brute force dictionary min search*)
(*usging a normal dictionary*)
DictionaryMinSearch[{dictf_, valsf_}, ydat_] := DictionaryMinSearchi[{dictf, valsf}, ydat]

DictionaryMinSearchi[{_, valsf_}, {0. ..}] := Join[0. valsf[[1]], {0., 0., 0.}]

DictionaryMinSearchi[{dictf_, valsf_}, ydat_] := Block[{fwf, residualError, sol},
	(*calcualte dictionary error*)
 	fwf = LeastSquaresC[dictf, ydat];
 	residualError = ErrorC[ydat, fwf, dictf];
 	sol = First@Ordering[residualError, 1];
 	(*find Min value*)
 	Flatten[{valsf, fwf, residualError}[[All, sol]]]
  ]


(*brute force dictionary min search*)
(*usging a the pseudo inverse dictionary*)
DictionaryMinSearch[{dictf_, dictfMat_, valsf_}, ydat_] := DictionaryMinSearchi[{dictf, dictfMat, valsf}, ydat]

DictionaryMinSearchi[{_, _, valsf_}, {0. ..}] := Join[0. valsf[[1]], {0., 0., 0.}]

DictionaryMinSearchi[{dictf_, dictfMat_, valsf_}, ydat_] := Block[{fwf, residualError, sol},
	(*calcualte dictionary error*)
 	fwf = LeastSquares2C[dictfMat, ydat];
 	residualError = ErrorC[ydat, fwf, dictf];
 	sol = First@Ordering[residualError, 1];
 	(*find Min value*)
 	Flatten[{valsf, fwf, residualError}[[All, sol]]]
  ]
  
  
(*minimization dictionary min search*)
(*usging a normal dictionary*)  

dicmet = {"NelderMead", "PostProcess" -> False, "ExpandRatio" -> 2, "ShrinkRatio" -> .75, "ReflectRatio" -> .85, "ContractRatio" -> .75};

DictionaryMinSearch[{dict_, vals_}, ydat_, cons_] := DictionaryMinSearchi[{dict, vals}, ydat, cons]

(*if background return zeors*)
DictionaryMinSearchi[{_, _}, {0. ..}, {_, _, _}] = {0., 0., 0., 0., 0., 0.}

DictionaryMinSearchi[{_, _}, {0. ..}, {_, _}] = {0., 0., 0., 0., 0.}

(*Fit without the predefined pseudo inverse*)
DictionaryMinSearchi[{dict_, vals_}, ydat_, {maxx_,maxy_}] := Block[{err, coor, ErrorFuncDic},
	
	(*define the cost function*)
	ErrorFuncDic[sig_, {x_Integer, y_Integer}] := LeastSquaresErrorC[dict[[x, y]], sig];
	(*minimize the dictionary value*)
	{err, coor} = Quiet@NMinimize[
		{ErrorFuncDic[ydat, {x, y}], 1 <= x <= maxx && 1 <= y <= maxy}, {x, y}, Integers, 
		Method -> Append[dicmet, "InitialPoints" -> Round[{maxx, maxy}/2]]];
	(*find the min vals*)
	Flatten[{vals[[##]], LeastSquaresC[dict[[##]], ydat], err}] & @@ coor[[All, 2]]
	]

(*Fit with the predifined pseudo inverse*)
DictionaryMinSearchi[{dict_, vals_}, ydat_, {maxx_, maxy_, maxz_}] := Block[{err, coor, ErrorFuncDic},
	(*define the cost function*)
	ErrorFuncDic[sig_, {x_Integer, y_Integer, z_Integer}] := LeastSquaresErrorC[dict[[x, y, z]], sig];
	(*minimize the dictionary value*)
	{err, coor} = Quiet@NMinimize[
		{ErrorFuncDic[ydat, {x, y, z}], 1 <= x <= maxx && 1 <= y <= maxy && 1 <= z <= maxz}, {x, y, z}, Integers,
		Method -> Append[dicmet, "InitialPoints" -> Round[{maxx, maxy, maxz}/2]]];
	(*find the min vals*)
	Flatten[{vals[[##]], LeastSquaresC[dict[[##]], ydat], err}] & @@ coor[[All, 2]]
  ]


(*minimization dictionary min search*)
(*usging a the pseudo inverse dictionary*)  
DictionaryMinSearch[{dict_, dictMat_, vals_}, ydat_, cons_] := DictionaryMinSearchi[{dict, dictMat, vals}, ydat, cons]

DictionaryMinSearchi[{_, _, _}, {0. ..}, {_, _, _}] = {0., 0., 0., 0., 0., 0.}

DictionaryMinSearchi[{_, _, _}, {0. ..}, {_, _}] = {0., 0., 0., 0., 0.}

DictionaryMinSearchi[{dict_, dictMat_, vals_}, ydat_, {maxx_,maxy_}] := Block[{err, coor, ErrorFuncDic},
	(*define the cost function*)
	ErrorFuncDic[sig_, {x_Integer, y_Integer}] := LeastSquaresError2C[dict[[x, y]], dictMat[[x, y]], sig];
	(*minimize the dictionary value*)
	{err, coor} = Quiet@NMinimize[
		{ErrorFuncDic[ydat, {x, y}], 1 <= x <= maxx && 1 <= y <= maxy}, {x, y}, Integers, 
		Method -> Append[dicmet, "InitialPoints" -> Round[{maxx, maxy}/2]]];
	(*find the min vals*)
	Flatten[{vals[[##]], LeastSquaresC[dict[[##]], ydat], err}] & @@ coor[[All, 2]]
	]

DictionaryMinSearchi[{dict_, dictMat_, vals_}, ydat_, {maxx_, maxy_, maxz_}] := Block[{err, coor, ErrorFuncDic},
	(*define the cost function*)
	ErrorFuncDic[sig_, {x_Integer, y_Integer, z_Integer}] := LeastSquaresError2C[dict[[x, y, z]], dictMat[[x, y, z]], sig];
	(*minimize the dictionary value*)
		{err, coor} = Quiet@NMinimize[
			{ErrorFuncDic[ydat, {x, y, z}], 1 <= x <= maxx && 1 <= y <= maxy && 1 <= z <= maxz}, {x, y, z}, Integers,
			Method -> Append[dicmet, "InitialPoints" -> Round[{maxx, maxy, maxz}/2]]];
	(*find the min vals*)
	Flatten[{vals[[##]], LeastSquaresC[dict[[##]], ydat], err}] & @@ coor[[All, 2]]
  ]


(* ::Subsubsection::Closed:: *)
(*ErrorFunc*)


ErrorFunc1[y_, T2f_Real, B1_Real, S0_, vals_] :=  Quiet@Block[{sig, T1m, T1f, echo, angle},   
	{echo, T1f, angle} = vals;
	Total[(y - S0 EPGSignali[echo, {T1f, T2f}, angle, B1])^2]]


ErrorFunc[y_, T2m_Real, B1_Real, vals_] := Quiet@Block[{sig, T1m, T1f, T2f, echo, angle},
   {echo, {T1m, T1f, T2f}, angle} = vals;
   sig = Transpose[{
      Abs[EPGSignali[echo, {T1m, T2m}, angle, B1]],
      Abs[EPGSignali[echo, {T1f, T2f}, angle, B1]]
      }];
   LeastSquaresErrorC[sig, y]]


ErrorFunc[y_, T2m_Real, T2f_Real, B1_Real, vals_] := Quiet@Block[{sig, T1m, T1f, echo, angle},
   {echo, {T1m, T1f}, angle} = vals;
   sig = Transpose[{
      EPGSignali[echo, {T1m, T2m}, angle, B1],
      EPGSignali[echo, {T1f, T2f}, angle, B1]
      }];
   LeastSquaresErrorC[sig, y]
   ]


(* ::Subsubsection::Closed:: *)
(*LeastSquaresC*)


(*calculate the pseudoinverse Ai*)
PseudoInverseC = Compile[{{A, _Real, 2}}, Block[{T = Transpose[A]}, (Inverse[T.A].T)], 
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

PseudoInverseWC = Compile[{{A, _Real, 2}, {W, _Real, 2}}, Block[{T = Transpose[A]}, (Inverse[T.W.A].T.W)], 
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];
   
LeastSquaresC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]}, (Inverse[T.A].T).y], 
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

(*Ai is Inverse[T.A].T*)   
LeastSquares2C = Compile[{{Ai, _Real, 2}, {y, _Real, 1}},  Ai.y, 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

(*f is Ai.y*)
ErrorC = Compile[{{y, _Real, 1}, {f, _Real, 1}, {A, _Real, 2}}, Total[((y - A.f))^2], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

LeastSquaresErrorC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]}, Total[(y - A.(Inverse[T.A].T).y)^2]], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

(*Ai is Inverse[T.A].T*)	
LeastSquaresError2C = Compile[{{A, _Real, 2}, {Ai, _Real, 2}, {y, _Real, 1}}, Total[(y - A.Ai.y)^2], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
