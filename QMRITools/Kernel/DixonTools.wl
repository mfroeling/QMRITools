(* ::Package:: *)

(* ::Title:: *)
(*QMRITools DixonTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`DixonTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`DixonTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


DixonToPercent::usage = 
"DixonToPercent[water, fat] converts the dixon water and fat data to percent maps.

Output is {waterFraction, fatFraction}.
The values of water and fat are arbitraty units and the ouput fractions are between 0 and 1."

DixonReconstruct::usage = 
"DixonReconstruct[{real, imag}, echo] reconstruxt Dixon data with initital guess b0 = 0 and T2star = 0.
DixonReconstruct[{real, imag}, echo, {b0}] reconstructs Dixon data with intitial guess T2star = 0.
DixonReconstruct[{real, imag}, echo, {b0, t2}] reconstructs Dixon data with t2star and B0.
DixonReconstruct[{real, imag}, echo, {b0, t2, ph0}] reconstructs Dixon data with initial phase.
DixonReconstruct[{real, imag}, echo, {b0, t2, ph0, phb}] reconstructs Dixon data with bipolar phase .

Output is {{watF,fatF},{watSig,fatSig},{inphase,outphase},{{b0, ph0, phb},{r2, t2}},itterations}.

The fractions are between 0 and 1, the B0 field map is in Hz and the T2start map is in ms.

DixonReconstruct[] is based on DOI: 10.1002/mrm.20624 and 10.1002/mrm.21737 (10.1002/nbm.3766)."


FixDixonFlips::usage = 
"FixDixonFlips[{mag, phase, real, imag}] checks if any volumes are 180 degrees out of phase and corrects them."


SimulateDixonSignal::usage = 
"SimulateDixonSignal[echo, fr, B0, T2] simulates an Dixon gradient echo sequence with echotimes.
Echotimes echo in ms, fat fraction fr between 0 and 1, field of resonance B0 in Hz and relaxation T2 in ms."

OptimizeDixonEcho::usage = 
"OptimizeDixonEcho[] shows a manipulate pannel which allos to optimize the dixon echos.
OptimizeDixonEcho[echos] shows a manipulate pannel which allos to optimize the predifined dixon echos."


FindInPhaseEchos::usage = 
"FindInPhaseEchos[echos, iop] finds the two nearest echos to inphase which are best used for unwrapping using the iop time."

Wrap::usage = 
"Wrap[data] wraps phase values between -Pi and Pi."

Unwrap::usage = 
"Unwrap[data] unwraps the given dataset. The data should be between -Pi and Pi. 
Unwrap[] is based on DOI: 10.1364/AO.46.006623 and 10.1364/AO.41.007437."

UnwrapSplit::usage = 
"UnwrapSplit[phase, data] unwarps the give phase dataset but splits the data into left and right using SplitData based in the data and performs the unwrapping seperately. The data should be between -Pi and Pi.
UnwrapSplit[] is based on DOI: 10.1364/AO.46.006623 and 10.1364/AO.41.007437."

UnwrapList::usage = 
"UnwrapList[list] unwraps a 1D list of values between -Pi and Pi."

UnwrapDCT::usage = 
"UnwrapDCT[data] unwraps the given dataset using DCT transform . The data should be between -Pi and Pi. 
UnwrapDCT[] is based on DOI: 10.1364/JOSAA.11.000107."

DixonPhase::usage = 
"DixonPhase[real, imag, echos] calculates the b0 and ph0 maps."


(* ::Subsection::Closed:: *)
(*Options*)


DixonPrecessions::usage = 
"DixonPrecessions is an options for DixonReconstruct. Defines the rotation of the signal {-1,1} default is -1."

DixonFieldStrength::usage = 
"DixonFieldStrength is an options for DixonReconstruct. Defines the fieldstrengths in Tesla on which the data was acquired."

DixonFrequencies::usage = 
"DixonFrequencies is an options for DixonReconstruct. Defines the frequencies in ppm of the fat peaks being used."

DixonAmplitudes::usage = 
"DixonAmplitudes is an options for DixonReconstruct. Defines the relative amplitudes of the fat peaks being used."

DixonNucleus::usage = 
"DixonNucleus is an option for DixonReconstruct. Defines the nucleus for which the reconstruction is performed."

DixonTollerance::usage = 
"DixonTollerance is an options for DixonReconstruct. Defines at which change per itteration of b0 and R2star the ittarative methods stops. Default value is 0.1."

DixonMaskThreshhold::usage = 
"DixonMaskThreshhold is an options for DixonReconstruct. Defines at which threshhold the dixon reconstruction considers a voxel to be background noise. Defualt values is 0.05."

DixonFilterInput::usage = 
"DixonFilterInput is an options for DixonReconstruct. If True the input b0 and T2star values are smoothed using a gaussian kernel."

DixonFilterOutput::usage = 
"DixonFilterOutput is an options for DixonReconstruct. If True the out b0 and T2star values are smoothed Median filter and lowpassfiltering after which the water and fat maps are recomputed."

DixonComplexOutput::usage = 
"DixonComplexOutput is an options for DixonReconstruct. If set True DixonReconstruct will return the complex data instead of the magnitude data."


DixonFilterSize::usage = 
"DixonFilterSize is an options for DixonReconstruct. Defines the number of voxel with which the input b0 and T2star values are smoothed."

DixonIterations::usage = 
"DixonIterations is an options for DixonReconstruct. Defines the maximum itterations the fit can use."

DixonBipolar::usage = 
"DixonBipolar is an option for DixonReconstruct. If set to true it assumes alternating readout directions."

DixonClipFraction::usage =
"DixonClipFraction is an option for DixonReconstruct. If set true the fat fraction is clipped between 0 and 1."

DixonFilterType::usage = 
"DixonFilterType is an option for DixonReconstruct. FilterType can me \"Median\" or \"Laplacian\"."

DixonCorrectT1::usage = 
"DixonCorrectT1 is an option for DixonReconstruct. To perform T1 correction provide the TR and FA as a list, {TR, FA}. TR is in ms and FA in degrees."

MonitorUnwrap::usage = 
"MonitorUnwrap is an option for Unwrap. Monitor the unwrapping progress."

UnwrapDimension::usage = 
"UnwrapDimension is an option for Unwrap. Can be \"2D\" or \"3D\". 2D is for unwarpping 2D images or unwrapping the individual images from a 3D dataset \
(does not unwrap in the slice direction). 3D unwraps a 3D dataset in all dimensions."

UnwrapThresh::usage = 
"UnwrapThresh is an option for Unwrap. Is a value between 0.6 and 0.9, and defines when to unwrap, the higher the value the less unwrapping will be done."


(* ::Subsection::Closed:: *)
(*Error Messages*)


Unwrap::data2D = "Unwrapping Dimensions is 2D and this can only be preformed on 2D or 3D data. This data is `1`D."

Unwrap::data3D = "Unwrapping Dimensions is 3D and this can only be preformed on 3D data. This data is `1`D."

Unwrap::dim = "Unwrapping Dimensions can be \"2D\" or \"3D\", current value is `1`."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*DixonToPercent*)


SyntaxInformation[DixonToPercent] = {"ArgumentsPattern" -> {_, _}};

DixonToPercent[water_, fat_]:=DixonToPercent[water, fat, True]

DixonToPercent[water_, fat_, clip_?BooleanQ] := Block[{atot, fatMap, waterMap, fMask, wMask, afat, awater},
	afat = Abs[fat];
	awater = Abs[water];
	atot = Abs[fat + water];
	
	(*define water and fat maps*)
	waterMap = Chop[DevideNoZero[awater, atot, "Comp"]];
	fatMap = Chop[DevideNoZero[afat, atot, "Comp"]];
	
	(*noise bias correction*)
	wMask = Mask[waterMap, .5, MaskSmoothing->False];
	fMask = (1 - wMask) Mask[fatMap, .5, MaskSmoothing->False];
	
	If[clip,
		fatMap = fMask (1 - waterMap) + wMask fatMap;
		waterMap = (wMask + fMask) - fatMap;
		Clip[{waterMap, fatMap}, {0, 1}],
		waterMap = wMask waterMap + (fMask - fMask fatMap);
		fatMap = (wMask + fMask) - waterMap;
		waterMap = (wMask + fMask) - fatMap;
		Clip[{waterMap, fatMap}, {-0.1, 1.1}, {-0.1, 1.1}]
 	]
]


(* ::Subsection:: *)
(*DixonReconstruct*)


(* ::Subsubsection::Closed:: *)
(*FixDixonFlips*)


SyntaxInformation[FixDixonFlips] = {"ArgumentsPattern" -> {{_, _, _, _}}};

FixDixonFlips[{mag_,phase_,real_,imag_}]:=Block[{p,r,i,c},
	p = FindComplexFlips[real,imag];
	{If[p==={},
		{mag,phase,real,imag}, 
		r = real;r[[All,p]]=-r[[All,p]];
		i = imag;i[[All,p]]=-i[[All,p]];
		Through[{Abs,Arg,Re,Im}[r+I i]]
	],p}
]


FindComplexFlips[real_,imag_]:=Block[{comp,diff,diffM,means},
	comp=real+I imag;
	diff=Differences[Transpose[comp]];
	diffM=Median@Abs@diff;
	means=MeanNoZero[Flatten[#]]&/@Abs[#-diffM&/@Abs[diff]];
	Flatten[Position[Partition[Append[Boole[# > 2 Median[means] & /@ means],1], 2, 1], {1, 1}]] + 1
]


(* ::Subsubsection:: *)
(*DixonReconstruct*)


Options[DixonReconstruct] = {
	DixonPrecessions -> -1, 
	DixonFieldStrength -> 3, 
	DixonNucleus -> "1H",
	(*DixonFrequencies -> {{0}, {3.8, 3.4, 3.13, 2.67, 2.46, 1.92, 0.57, -0.60}},
	DixonAmplitudes -> {{1}, {0.089, 0.598, 0.047, 0.077, 0.052, 0.011, 0.035, 0.066}},*)
	DixonFrequencies -> {{0}, {3.8, 3.4, 3.1, 2.7, 2.5, 1.95, 0.5, -0.5, -0.6}},
	DixonAmplitudes -> {{1}, {0.088, 0.628, 0.059, 0.064, 0.059, 0.01, 0.039, 0.01, 0.042}},
	DixonIterations -> 20, 
	DixonTollerance -> .1, 
	DixonMaskThreshhold -> 0.1,
	DixonFilterInput -> True, 
	DixonFilterOutput -> True,
	DixonFilterSize -> 1,
	DixonFilterType -> "Median",
	DixonBipolar -> False,
	DixonClipFraction -> False,
	DixonCorrectT1->False,
	DixonComplexOutput->False,
	MonitorCalc -> False
	};


SyntaxInformation[DixonReconstruct] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

DixonReconstruct[{real_, imag_}, echo_, opts : OptionsPattern[]] := DixonReconstruct[{real, imag}, echo, {0, 0, 0, 0}, opts]

DixonReconstruct[{real_, imag_}, echo_, {b0i_}, opts : OptionsPattern[]] := DixonReconstruct[{real, imag}, echo, {b0i, 0, 0, 0}, opts]

DixonReconstruct[{real_, imag_}, echo_, {b0i_, t2i_}, opts : OptionsPattern[]] := DixonReconstruct[{real, imag}, echo, {b0i, t2i, 0, 0}, opts]

DixonReconstruct[{real_, imag_}, echo_, {b0i_, t2i_, ph0i_}, opts : OptionsPattern[]] := DixonReconstruct[{real, imag}, echo, {b0i, t2i, ph0i, 0}, opts]


DixonReconstruct[{real_, imag_}, echo_, {b0i_, t2i_, ph0i_, phbi_}, OptionsPattern[]] := Block[{
		freqs, amps, eta, maxItt, thresh, filti, filto, filtFunc, mat, iop, Amat, Amati, ioAmat,
		complex, mask, range, zero, b0, r2, phi, ph0, phb, result, cWat, cFat, phEst ,ph0Est, phbEst, res, itt,
		fraction, signal, ioPhase, fit, t2, func, mon, n, t1c, tr, fa, t1f, t1m, sf, sw
	},
	
	(*algorithem is base on: *)	
	(*Triplett WT et.al. 10.1002/mrm.23917 - fat peaks*)
	
	(*Reeder et.al. 10.1002/mrm.20624 - iDEAL*)
	(*Huanzhou Yu et.al. 10.1002/jmri.21090 - iDEAL algorithm and T2* correction*)
	(*Bydder et.al. 10.1016/j.mri.2010.08.011 - initial phase*)
	(*Peterson et.al.10.1002/mrm.24657 - bipolar*)
	
	(*Chebrolu et.al. 10.1002/mrm.22300 - two t2* species*)
	(*Byder et.al. 10.1016/j.mri.2011.07.004 - a matrix with bonds*)
	
	mon = OptionValue[MonitorCalc];
	
	(*metabolite frequencies and amplitudes*)
	freqs = OptionValue[DixonPrecessions] OptionValue[DixonFieldStrength] OptionValue[DixonFrequencies] GyromagneticRatio[OptionValue[DixonNucleus]];
	amps = #/Total[#] &/@ OptionValue[DixonAmplitudes];
	iop = {0, 0.5}/Abs[Total[(amps[[2]]^2) freqs[[2]]]/Total[amps[[2]]^2]];
	
	(*optimization settings*)
	eta = OptionValue[DixonTollerance];
	maxItt = OptionValue[DixonIterations];
	thresh = OptionValue[DixonMaskThreshhold];
	
	(*define filter*)
	filti = OptionValue[DixonFilterInput];
	filto = OptionValue[DixonFilterOutput];
	func = Switch[OptionValue[DixonFilterType],"Median", MedFilter, "Laplacian",LapFilter];
	filtFunc = func[#, OptionValue[DixonFilterSize]]&;
	
	(*Get the T1 correction Factor*)
	t1c = If[OptionValue[DixonCorrectT1]===False, 1,
		{tr, fa} = OptionValue[DixonCorrectT1];
		{t1f, t1m} = {400, 1400};
		sf = ((1. - Exp[-tr/t1f]) Sin[fa Degree])/(1 - Exp[-tr/t1f] Cos[fa Degree]);
		sw = ((1. - Exp[-tr/t1m]) Sin[fa Degree])/(1 - Exp[-tr/t1m] Cos[fa Degree]);
		sf/sw
	];

	(*get alternating readout signs for bipolar acquisition*)
	(*mat = (2 Pi I) Transpose[{echo, (If[OptionValue[DixonBipolar], -1, 1])^Range[Length[echo]]}];*)
	n = If[OptionValue[DixonBipolar], 3, 2];
	mat = (2. Pi I) Transpose[{echo, ConstantArray[1, Length[echo]], (-1)^Range[Length[echo]]}][[All,;;n]];

	(*define the water fat matrixes*)
	Amat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ echo;
	Amati = Inverse[ConjugateTranspose[Amat] . Amat] . ConjugateTranspose[Amat];
	ioAmat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ iop;
			
	(*create complex data for fit*)
	If[mon, PrintTemporary["Prepairing data and field maps"]];
	complex = N[real + imag I];

	If[ArrayDepth[complex] === 4, complex = Transpose[complex]];
	mask = Closing[UnitStep[Abs[First@complex] - thresh], 1] Unitize[Total[Abs[complex]]];
	range = {0., 2 Max[Abs[complex]]};
	zero = 0. Re[First@complex];
	
	(*define complex field map and background phase*)
	{b0, t2, ph0, phb} = mask If[# === 0, zero, If[filti, filtFunc[#], #]]& /@ {b0i, t2i, ph0i, phbi};
	
	(*perform the dixon reconstruction*)
	If[mon, PrintTemporary["Performing dixon iDEAL reconstruction"]];
	complex = RotateDimensionsLeft@MaskData[complex, mask];
	phi = RotateDimensionsLeft[{(b0 + I DevideNoZero[1., t2i] / (2 Pi)), ph0, phb}[[;;n]]];
	result = DixonFitiC[complex, phi, mask, mat, Amat, Amati, eta, maxItt, n];
 	If[n==2,
 		{cWat, cFat, phEst ,ph0Est, res, itt} = RotateDimensionsRight[Chop[result]];
 		phbEst = zero;
 		,
 		{cWat, cFat, phEst ,ph0Est, phbEst, res, itt} = RotateDimensionsRight[Chop[result]];
 	];

	(*filter the output*)
	If[filto,
		If[mon, PrintTemporary["Filtering field estimation and recalculating signal fractions"]];
		(*smooth b0 field and R2star maps*)
		(*phiEst = mask (filtFunc[Re[phiEst]] + DevideNoZero[1. ,filtFunc[Ramp[DevideNoZero[1., Im[phiEst]]]]] I);*)
		phEst = mask (filtFunc[Re[phEst]] + I DevideNoZero[1., filtFunc[DevideNoZero[1.,filtFunc[Ramp[Im[phEst]]]]]]);
		ph0Est = mask filtFunc[Re[ph0Est]](*(filtFunc[Re[ph0Est]] + I filtFunc[Im[ph0Est]])*);
		phbEst = mask filtFunc[Re[phbEst]](*(filtFunc[Re[phbEst]] + I filtFunc[Im[phbEst]])*);
		phi = RotateDimensionsLeft[{phEst ,ph0Est, phbEst}[[;;n]]];
		
		(*recalculate the water fat signals*)
		result = DixonFitiC2[complex, phi, mask, mat, Amat, Amati];
		{cWat, cFat, res} = RotateDimensionsRight[Chop[result]];
	];
	
	(*create the output*)
	If[mon, PrintTemporary["Generating all the needed outputs"]];
	itt = Round@Re@itt;
	mask = Times @@ (Mask[Abs[#], 2 range, MaskSmoothing -> False] & /@ {cWat, cFat});
	{cWat, cFat} = {t1c mask cWat, mask cFat};
	fraction = DixonToPercent[cWat, cFat, OptionValue[DixonClipFraction]];

	(*signal and in/out phase data *)
	signal = If[OptionValue[DixonComplexOutput],
		1000 {cWat, cFat} / range[[2]],
		1000 Clip[Abs[{cWat, cFat}], range] / range[[2]]
	];	
	res = 1000 Clip[Abs[res], range] / range[[2]];
	ioPhase = 1000 Clip[RotateDimensionsRight[InOutPhase[cWat, cFat, ioAmat]], range]  / range[[2]];
	
	(*estimate b0 and t2star*)
	b0 = Re[phEst];
	r2 = 2 Pi Im[phEst];
	t2 = DevideNoZero[1, r2];
	phi = {Clip[b0, {-400., 400.}, {-400., 400.}], ph0Est, phbEst}; 
	t2 = {Clip[t2, {0., 0.25}, {0., 0.25}],	Clip[r2, {0., 1000.}, {0., 1000.}]};
	
	(*give the output*)
	{fraction, signal, ioPhase, {phi, t2}, itt, res}
]


InOutPhase = Compile[{{cWat, _Complex, 0}, {cFat, _Complex, 0}, {ioAmat, _Complex, 2}},
	Abs[ioAmat . {cWat, cFat}]
, RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}];


(* ::Subsubsection:: *)
(*DixonFit*)


DixonFitiC = Compile[{
		{ydat, _Complex, 1}, {phi, _Complex, 1}, {mask, _Real, 0},
		{mat,_Complex,2}, {Amat, _Complex, 2}, {Amati, _Complex, 2}, 
		{eta, _Real, 0}, {maxItt, _Integer, 0}, {n, _Integer, 0}
	}, Block[{
		i, continue, phiEst, phi0Est, dPhi, dPhi0, rms, res, rho, sigd, sol, Bi, sc
	},
	
	(*initialize variables*)
	res = ydat 0.;
	rho = {0. + 0. I, 0. + 0. I};
	phiEst = phi;
	dPhi = 0 phi;
	rms = 0. I;

	(*loop parameters*)
	i = 0;
	continue = True;
	
	(*perform itterative fitting*)
	If[mask > 0, 
		While[continue,	i++;
			(*update the field map*)
			phiEst += dPhi;
			
			(*define complex field map P(-phi) or (E D)^-1 and demodulate signal phase*)
			sigd = Exp[-mat . phiEst] ydat;
			
			(*perform A.2 form 10.1002/jmri.21090, the water fat fraction*)
			(*A.rho is needed for Matrix B eq A.4 and eq A.5, calculat the fitted demodulated signal*)
			rho = (Amati . sigd);
			sol = Amat . rho;
			res = sigd - sol;
			
			(*Define the matrix B including bipolar 10.1002/mrm.24657 and initial phase 10.1016/J.MRI.2010.08.011*)
			(*Obtain the error terms eq A.5*)
			Bi = PseudoInverse[Join[mat sol, Amat, 2]][[1 ;; n]];
			dPhi = Bi . res;
			(*only look at phase for initial and bipolar*)
			dPhi[[2;;n]] = Re@dPhi[[2;;n]];

			(*chech for continue*)
			continue = ! ((Abs[Re@dPhi[[1]]] < eta && Abs[(Im@dPhi[[1]])] < eta && Abs[10 Re@dPhi[[2]]] < eta && Abs[100 Re@dPhi[[n]]] < eta) || i >= maxItt);
		];
		
		rms = Sqrt[Mean[res^2]];
	];

	(*output*)
	Join[rho, phiEst, {rms, i}]
	
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


DixonFitiC2 = Compile[{
	{ydat, _Complex, 1}, {phi, _Complex, 1}, {mask, _Real, 0},
	{mat, _Complex, 2}, {Amat, _Complex, 2}, {Amati, _Complex, 2}}, 
	Block[{sigd, rho, res, rms, r1, r2},
		If[mask > 0,
			(*initialize values*)
			rho = {0. + 0. I, 0. + 0. I};
			res = ydat 0.;
			
			(*define complex field map P(-phi) or (E D)^-1, demodulate signal and find complex fraction*)
			sigd = Exp[-mat . phi] ydat;
			rho = Amati . sigd;
			
			(*calculate the residuals*)
			res = sigd - (Amat . rho);
						
			(*output*)
   			Join[rho, {Sqrt[Mean[res^2]]}],
   			ConstantArray[0., Length[Amat[[1]]] + 1]
		]
	]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection::Closed:: *)
(*SimulateDixonSignal*)


Options[SimulateDixonSignal] = {
	DixonNucleus -> "1H", 
	DixonPrecessions -> -1, 
	DixonFieldStrength -> 3, 
	DixonFrequencies -> {{0}, {3.8, 3.4, 3.13, 2.67, 2.46, 1.92, 0.57, -0.60}}, 
	DixonAmplitudes -> {{1}, {0.089, 0.598, 0.047, 0.077, 0.052, 0.011, 0.035, 0.066}}
}

SyntaxInformation[SimulateDixonSignal] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

SimulateDixonSignal[echo_, fr_, b0_, t2_, OptionsPattern[]] := Block[{precession,field,freqs,amps, Amat, phi, sig},
	precession = OptionValue[DixonPrecessions](*-1,1*);
	field = OptionValue[DixonFieldStrength];
	freqs = precession field GyromagneticRatio[OptionValue[DixonNucleus]] OptionValue[DixonFrequencies];
	amps = #/Total[#] & /@ OptionValue[DixonAmplitudes];
	
	Amat = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ echo;
	phi = N@2 Pi b0 I - 1./t2;
	sig = Exp[If[Length[t2] === 2, Transpose[echo # & /@ phi], phi echo]] Amat;
	
	sig = sig . {fr, 1 - fr};
	{Re[sig], Im[sig]}
]


(* ::Subsection::Closed:: *)
(*OptimizeDixonEcho*)


SyntaxInformation[UnwrapList] = {"ArgumentsPattern" -> {}};

Options[OptimizeDixonEcho] = {
	DixonPrecessions -> -1,
	DixonFieldStrength -> 3,
	DixonNucleus -> "1H",
	DixonFrequencies -> {{0}, {3.8, 3.4, 3.1, 2.7, 2.5, 1.95, 0.5, -0.5, -0.6}},
	DixonAmplitudes -> {{1}, {0.088, 0.628, 0.059, 0.064, 0.059, 0.01, 0.039, 0.01, 0.042}}
}
	
OptimizeDixonEcho[ops:OptionsPattern[]]:=OptimizeDixonEcho[1,1,ops]

OptimizeDixonEcho[fi_,di_,OptionsPattern[]]:=Manipulate[
	fr = 297.466/(field GyromagneticRatio["1H"]);
	echos = first + Range[0,necho-1]delta;
	pts = Transpose[Through[{Im,Re}[Exp[-2Pi echos/fr I]]]];
	e = FindInPhaseEchos[echos, fr,DixonBipolar -> bip];
	pre = Table[1/n -> Ceiling[n/2] fr/n, {n, 2, 7}];
	
	(*define the water and fat frequencies and amplitudes to calcluate the condition number*)
	freqs = OptionValue[DixonPrecessions] OptionValue[DixonFieldStrength] OptionValue[DixonFrequencies] GyromagneticRatio[OptionValue[DixonNucleus]];
	amps = #/Total[#] & /@ OptionValue[DixonAmplitudes];
	A = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ echos;
	condNum = Divide @@ SingularValueList[A];
	
	Column[{
		Switch[vis,
			"path",
			Show[
				Graphics[{Black,PointSize[.02],Point[pts]},PlotRange->{{-1.2,1.2},{-1.2,1.2}},AspectRatio->1,Axes->True,Ticks->None,ImageSize->300,
					PlotLabel->Column[{"Inphase echos for unwrapping: "<>ToString[e],"Condition Number: "<>ToString[Round[condNum, 0.001]]}]],
				ListLinePlot[
					Table[Callout[pts[[i]], i, LabelStyle -> {FontSize -> 14, Bold}, CalloutStyle -> None], {i, Length[pts]}],
					PlotRange->{{-1.2,1.2},{-1.2,1.2}}, PlotStyle->Black, Mesh->All],
				Graphics[{{Green,PointSize[.04],Point@First@pts},{Red,PointSize[.04],Point@Last@pts}}],
					If[io,Graphics[{{PointSize[.02],Blue,Point[pts[[e]]]}}],Graphics[]],
				Graphics[{Lighter@Gray, Dashed, Opacity[.5], Circle[], Line[{{0, 0}, RotationMatrix[# Degree] . {0, 1}}] & /@ Range[0, 360, 60]}]
			],
			"vectors",
			Show[
				Graphics[{Black,PointSize[.02],Point[pts]},PlotRange->{{-1.2,1.2},{-1.2,1.2}},AspectRatio->1,Axes->True,Ticks->None,ImageSize->300,PlotLabel->"Inphase echos for unwrapping: "<>ToString[e]],
				ListLinePlot[Table[Callout[pts[[i]], i, LabelStyle -> {FontSize -> 14, Bold}, CalloutStyle -> None], {i, Length[pts]}], PlotStyle -> Transparent,PlotRange->{{-1.2,1.2},{-1.2,1.2}}],
				Graphics[{{PointSize[.04],Green,Point@First@pts}}],
				Graphics[{{PointSize[.04],Red,Point@Last@pts}}],
				Graphics[{Thick,Arrow[{{0,0},#}&/@pts[[2;;-2]]]}],
				Graphics[{Thick,Green,Arrow[{{0,0},First@pts}]}],
				Graphics[{Thick,Red,Arrow[{{0,0},Last@pts}]}],
				If[io,Graphics[{Blue,Arrow[{{0,0},#}&/@pts[[e]]]}],Graphics[]],
				If[io,Graphics[{{PointSize[.02],Blue,Point[pts[[e]]]}}],Graphics[]],
				Graphics[{Lighter@Gray, Dashed, Opacity[.5], Circle[], Line[{{0, 0}, RotationMatrix[# Degree] . {0, 1}}] & /@ Range[0, 360, 60]}]
			]
		],
		ListLinePlot[Transpose@pts, Mesh -> All, PlotStyle -> {Red, Black}, ImageSize -> 300, Ticks -> None]
	}]
	,
	{{first, fi(*0.95*)(*fr*),"first echo"},0,2fr,0.005},
	{{delta, di(*1.3*),"echo spacing"},0.0,2fr,.005},
	{{delta, di, "echo spacing"}, Table[Ceiling[n/2] fr/n -> 1/n, {n, 2, 7}], ControlType -> SetterBar},
	{{necho,10,"number of echos"},2,25,1},
	Delimiter,
	{{field,3,"field strength"},{1,1.5,3,7,9.4}},
	Delimiter,
	{{vis,"path","visualization"},{"path","vectors"}},
	Delimiter,
	{{io,False,"show inphase echos"},{True,False}},
	{{bip,False,"bipolar"},{True,False}},
	
	{fr, ControlType->None},
	{echos, ControlType->None},
	{pts, ControlType->None},
	{freqs, ControlType->None},
	{amps, ControlType->None},
	{A, ControlType->None},
	{condNum, ControlType->None},
	
	
	ControlPlacement->Left,
	Initialization:>{
		field=3;
		fr=300/(field GyromagneticRatio["1H"])},
	SaveDefinitions->True
]


(* ::Subsection:: *)
(*Phase unwrap*)


(* ::Subsubsection::Closed:: *)
(*FindInPhaseEchos*)


Options[FindInPhaseEchos] = {DixonBipolar -> False}

SyntaxInformation[FindInPhaseEchos] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

FindInPhaseEchos[echos_, iop_, OptionsPattern[]]:=Block[{phase, bip, sel, list},
	phase = Arg@Exp[-I 2 Pi echos/iop];
	phase = Flatten[Table[{
		Abs@Nearest[phase[[{i, j}]], 0, 2],
		Abs[phase[[i]] - phase[[j]]], 
		Abs[i - j],
		Min[echos[[{i, j}]]]
		,{i, j}
	}, {i, 1, Length@echos}, {j, i + 1, Length@echos}], 1];
	
	bip = OptionValue[DixonBipolar];
	(*step 1 select all echos pairs where first is more inphase then out phase*)
	sel = Select[phase, (#[[1, 1]] < 0.5 Pi && #[[1, 2]] < 0.75 Pi) &];
	(*step 2 select all echo pairs that are less than 0.5 out phase*)
	sel = Select[sel, (#[[2]] < 0.5 Pi) &];
	(*step 3 If bipolar only select pairs with even difference*)
	If[bip, sel = Select[sel, EvenQ[#[[3]]] &]];
	
	If[sel=!={},
		(*get ranking for each of the values, smallest phase, phase diff, echo diff, first echo*)
		list = Range[Length@sel];
		sel = sel[[Ordering[Total /@ Transpose[(list /. Thread[Ordering[#] -> list]) & /@ Transpose[sel[[All, 1 ;; 4]]]]]]];
		sel[[1, -1]],
		{1,Length@echos}
	]
	
	(*phase=If[#>0.5 ,#-1,#]&/@FractionalPart[echos/iop];
	ord=Flatten[Position[phase,#]&/@Nearest[phase,0,Length[phase]]];
	Sort[If[OptionValue[DixonBipolar],
		Select[ord,If[OddQ[First@ord],OddQ,EvenQ]][[1;;2]],
	ord[[;;2]]]]*)
]


(* ::Subsubsection::Closed:: *)
(*UnwrapList*)


SyntaxInformation[UnwrapList] = {"ArgumentsPattern" -> {_}};

UnwrapList[list_]:=Block[{jumps,lst,diff,out},
	lst=If[Head[First@list]===Complex,Arg@list,list]/Pi;
	diff=Differences[lst];
	jumps=2 Prepend[Accumulate[(-Sign[diff]) Round[Chop[Abs[diff], 1.25]/2]],0];
	out=jumps+lst;
	Pi(Round[Subtract[Mean[list],Mean[out]],2]+out)
]


(* ::Subsubsection::Closed:: *)
(*Unwrap*)


Options[Unwrap]={MonitorUnwrap->True, UnwrapDimension->"2D", UnwrapThresh->0.5};

SyntaxInformation[Unwrap] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

Unwrap[dat_,OptionsPattern[]]:= Block[{data, ind, undim, out, mon, thresh, len},
	(* Phase unwrapping algorithem based on *)
	(*M.A. Herraez et al 2002 - 2D phase unwrapping using noncontinuous path -  DOI: 10.1364/ao.41.007437*)
	(*Abdul-Rahman 2007. - 3D phase unwrapping usiong noncontinuous path- DOI: 10.1364/AO.46.006623*)
	
	(*get options*)
	undim = OptionValue[UnwrapDimension];
	mon = OptionValue[MonitorUnwrap];
	thresh = Clip[OptionValue[UnwrapThresh],{0.3, 0.9}];
	
	(*swithc between 2D and 3D methods*)
	Switch[undim,
		
		(*2D using noncontinuous path*)
		"2D",
		data = ToPackedArray@Wrap[N[dat]];
		Which[
			MatrixQ[data],
			If[mon,PrintTemporary["Unwrapping one image using 2D algorithm."]];
			out = Unwrapi[data, thresh];
			,
			ArrayQ[data,3],
			len = Length[data];
			If[mon, PrintTemporary["Unwrapping ", len," images using 2D algorithm"]];
			out = Map[Unwrapi[#1, thresh]&, data, {ArrayDepth[data]-2}];
			out = UnwrapZi[out, thresh];
		],
		
		(*3D using noncontinuous path*)
		"3D",
		data = ToPackedArray@ArrayPad[Wrap[N[dat]], 1, 0.];
		If[ArrayQ[data,3],
			If[mon,PrintTemporary["Unwrapping 3D data using 3D algorithm"]];
			out = ArrayPad[Unwrapi[data, thresh],-1];
			,
			Message[Unwrap::data3D,ArrayDepth[data]]
		],
		
		(*Unknown option*)
		_,
		Message[Unwrap::dim,undim]
	];
		
	(*center around 0*)
	ToPackedArray@N@out
]


(* ::Subsubsection::Closed:: *)
(*UnwrapSplit*)


Options[UnwrapSplit] = Options[Unwrap]

SyntaxInformation[UnwrapSplit] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

UnwrapSplit[phase_, mag_,opts:OptionsPattern[]] := Block[{cutVal, phaseSplit, B0split},
	cutVal = CutData[mag][[3]];
	phaseSplit = CutData[phase, cutVal][[1 ;; 2]];
	B0split = Unwrap[#, opts] & /@ phaseSplit;
	StichData @@ B0split
]


(* ::Subsubsection::Closed:: *)
(*UnwrapZi*)


UnwrapZi[data_, thresh_]:= Block[{mask,slice,diff,meandiff,steps,off,unwrap,dat,num2,Roundi},
	Roundi = If[Negative[#], 
		num2=#-Ceiling[#]; If[-1<num2<-thresh,Floor[#],Ceiling[#]],
		num2=#-Floor[#]; If[1>num2>thresh,Ceiling[#],Floor[#]]
	]&;
	
	mask = Unitize[data];
	slice = Round[0.5Length[data]];
	diff = (#[[1]]-#[[2]]&/@ Partition[data/(2Pi),2,1]);
	
	meandiff = Median[#]&/@ Map[DeleteCases[Flatten[N[#]],0.]&,diff];
	steps = FoldList[Plus,0,Map[Roundi[#]&,meandiff]];
	off = Round[Median[DeleteCases[Flatten[N[data[[slice]]/(2Pi)]],0.]]];
	
	unwrap = steps-(steps[[slice]]+off);
	dat = (2Pi unwrap+data)mask//N;
	(dat - mask Round[MeanNoZero[Flatten[dat]],2Pi])
]


(* ::Subsubsection::Closed:: *)
(*Unwrapi*)


Unwrapi[dat_, thresh_] := Block[{data, mask, crp, dimi, sorted, groups, groupsize, groupnr, task},
	If[MinMax[dat]==={0.,0.},
		dat,
		
		(*rescale the data to 2Pi = 1, makes it easyer to process, removes needs for 2 PI checks*)
		data = dat / (2. Pi);
		dimi = Dimensions[data];
		
		(*remove zeros in back ground to reduce datasize in 3D*)
		data = If[ArrayDepth[data] == 3, crp = FindCrop[data]; ApplyCrop[data, crp], data];
		
		(*make mask to pervent unwrapping in 0 values*)
		mask = Mask[Ceiling[Abs@data], 1, MaskSmoothing ->False];
		
		(*Get the edges sotrted for reliability and precluster groups*)
		sorted = GetEdgeList[data, mask];
		{groups, groupsize} = MakeGroups[data, mask];
	
		(*make 2D data 3D and define shifts in add*)
		If[ArrayDepth[data] == 2,
			groups = {groups}; data = {data};
			sorted = Transpose[{#[[1]] + 1, 0 #[[1]] + 1, #[[2]], #[[3]]} &[Transpose[sorted]]];
		];
		
		(*Unwrap the data*)
		data = UnWrapC[sorted, data, groups, groupsize, thresh];
		
		(*make output in rad*)	
		If[ArrayDepth[dat] == 2, 
			(*output the 2D in rad*)
			2 Pi data[[1]],
			(*align to zero and ouput 3D in rad*)
			data = 2 Pi (data - mask Round[MeanNoZero[Flatten[data]]]);
			ReverseCrop[data, dimi, crp]
		]
	]
]


UnWrapC = Compile[{{sorted, _Integer, 2}, {datai, _Real, 3}, {groupsi, _Integer, 3}, {groupsizei, _Integer, 1}, {thresh,_Real,0}},
	Block[{data, const, dir, dim, groups, group1, group2, groupsize, groupnr, z1, z2, x1, x2, y1, y2, wrap, wrapT, pos, g1, g2, out, adds,add, diff},
		
		data = datai;
		groups = groupsi;
		groupsize = groupsizei;
		
		(*initialize parameters*)
		dim = Dimensions[data];
		groupnr = Length[groupsize];
		adds = {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};
		z1=x1=y1=z2=x2=y2=group1=group2=g1=g2=0;
		
		(*loop over all edges*)
		out = Map[(
			
			(*Get the voxel corrdinates and the neighbour, contrain to dimensions*)
			add = adds[[#[[1]]]];
			z1=#[[2]]; z2 = If[z1==dim[[1]], dim[[1]], z1+add[[1]]];
			x1=#[[3]]; x2 = If[x1==dim[[2]], dim[[2]], x1+add[[2]]];
			y1=#[[4]]; y2 = If[y1==dim[[3]], dim[[3]], y1+add[[3]]];
			        
			(*Get the group numbers*)
			group1 = groups[[z1, x1, y1]];
			group2 = groups[[z2, x2, y2]];
			
			(*Unwrapping logic*)
			(*0. Only process if both are not in background*)
			If[group1 > 0 && group2 > 0,
				
				(*1. Only process if both are not in same group or no group*)
				If[group1 != group2 || group1 == group2 == 1,
					
				(*Determine the wrap of the edge*)
				diff = data[[z1, x1, y1]] - data[[z2, x2, y2]];
				wrap = Sign[diff] Ceiling[Abs[diff] - thresh];
				wrapT = (wrap != 0);
				
				(*2. both already in a group*)
				If[group1 > 1 && group2 > 1,
					g1 = groupsize[[group1]];
					g2 = groupsize[[group2]];
					If[g1 >= g2,
						(*2A. group 1 is larger, add group 2 to group 1*)
						pos = BitXor[1, Unitize[groups - group2]];
						groups += ((group1 - group2) pos);
						groupsize[[group1]] = g1 + g2;
						groupsize[[group2]] = 0;
						If[wrapT, data += wrap pos];
						,
						(*2B. group 2 larger, add group 1 to group 2*)
						pos = BitXor[1, Unitize[groups - group1]];
						groups += ((group2 - group1) pos);
						groupsize[[group1]] = 0;
						groupsize[[group2]] = g1 + g2;
						If[wrapT, data -= wrap pos];
					],
					
					(*3. one of two pixels not in group*)
					Which[
						(*3A. only group 1 existst, add group 2 to group 1*)
						group1 > 1,
						groups[[z2, x2, y2]] = group1;
						groupsize[[group1]] += 1;
						If[wrapT, data[[z2, x2, y2]] += wrap];
						,
						(*3B. only group 2 existst, add group 1 to group 2*)
						group2 > 1,
						groups[[z1, x1, y1]] = group2;
						groupsize[[group2]] += 1;
						If[wrapT, data[[z1, x1, y1]] -= wrap];
						,
						(*3C. both belong to no group, make new group*)
						True,
						groupnr++;
						groups[[z1, x1, y1]] = groups[[z2, x2, y2]] = groupnr;
						AppendTo[groupsize,2];
						If[wrapT, data[[z2, x2, y2]] += wrap];
					]
				]
			]
		];
		(*close the map fucntion*)
		1) &, sorted];
		
	(*output the unwraped data*)
	data],
RuntimeOptions -> "Speed", Parallelization -> True];



(* ::Subsubsection::Closed:: *)
(*GetEdgeList*)


GetEdgeList[data_] := GetEdgeList[data, 1]

GetEdgeList[data_, maski_] := Block[{dep, diff, diff1, diff2, mask, edge, coor, fedge, ord, pos},
	dep = ArrayDepth[data];
	(*maske a mask if needed*)
	mask = If[maski === 1, Mask[Ceiling[Abs@data], 1, MaskSmoothing -> False], maski];
	
	(*calculate the second order diff*)
	{diff1, diff2} = Transpose[Partition[DiffU[ListConvolve[#, data, {2, -2}]] & /@ GetKernels[dep], 2]];
	diff = Total[(diff1 + diff2)^2];
	
	(*get the edge reliability*)
	edge = (RotateLeft[diff, #] + diff) & /@ Switch[dep, 2, {{0, 1}, {1, 0}}, 3, {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}];
	edge = MaskData[edge, mask];
	 
	(*sort the edges for reliability*)
	coor = MapIndexed[#2 &, edge, {dep + 1}];
	fedge = Flatten[edge, dep];
	ord = Ordering[fedge];
	pos = Position[Unitize[fedge[[ord]]], 1, 1, 1];
	pos=If[pos==={},1,pos[[1,1]]];
	
	Flatten[coor, dep][[ord]][[pos ;;]]
]


GetKernels[dep_] := Block[{ker, i, j, k, keri},
	Switch[dep,
		2,
		ker = ConstantArray[0, {3, 3}];
		ker[[2, 2]] = -1;
		({i, j} = #; keri = ker; keri[[i, j]] = 1; keri) & /@ {
			{2, 1}, {2, 3}, {1, 2}, {3, 2}, {1, 1}, {3, 3}, {1, 3}, {3, 1}
		}(*H,V,D1,D2*),
		3,
		ker = ConstantArray[0, {3, 3, 3}];
		ker[[2, 2, 2]] = -1;
		({i, j, k} = #; keri = ker; keri[[i, j, k]] = 1; keri) & /@ {
			{2, 1, 2}, {2, 3, 2}, {1, 2, 2}, {3, 2, 2}, {2, 2, 1}, {2, 2, 3},(*H,V,N*)
			{1, 1, 2}, {3, 3, 2}, {1, 3, 2}, {3, 1, 2},(*plane 1*)
			{1, 2, 1}, {3, 2, 3}, {1, 2, 3}, {3, 2, 1},(*plane 2*)
			{2, 1, 1}, {2, 3, 3}, {2, 1, 3}, {2, 3, 1},(*plane 3*)
			{1, 1, 1}, {3, 3, 3}, {1, 1, 3}, {3, 3, 1}, {1, 3, 1}, {3, 1, 3}, {3, 1, 1}, {1, 3, 3}(*diagonals*)
		}
	]
];


DiffU = Compile[{{diff, _Real, 0}}, 
	diff - Sign[diff] Ceiling[ Abs[diff] - 0.5]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MakeGroups*)


MakeGroups[data_] := MakeGroups[data, 1]

MakeGroups[data_, maski_]:=Block[{dep,dim,fun,min,max,part,dat,masks,small,nclus,clus,groupsize,groups,groupnr,mask},

	(*maske a mask if needed*)
	mask = If[maski === 1, Mask[Ceiling[Abs@data], 1, MaskSmoothing -> False], maski];
	dat = (mask data) - 2 (1 - mask);
		
	(*find mask ranges*)
	{min,max} = MinMax[data];
	part = {#[[1]] + 0.001, #[[2]] - 0.001} & /@ Partition[Range[-1, 1, 0.2] // N, 2, 1];
		
	(*make groups from masks*)
	clus = DeleteSmallComponents[MorphologicalComponents[
		Mask[dat, #, MaskSmoothing -> False], CornerNeighbors -> False], 
		If[ArrayDepth[data]===3, 15, 3], CornerNeighbors -> False] & /@ part;
		
	nclus = Prepend[Drop[Accumulate[Max /@ clus], -1], 0];
	groups = Total[MapThread[#2 Unitize[#1] + #1 &, {clus, nclus}]] + mask;
	
	(*create outputs, the size vector and group nrs*)
	groupsize = ConstantArray[0, Max[groups]];
	(groupsize[[#[[1]]]] = #[[2]]) & /@ Sort[Tally[Flatten[groups]]][[2 ;;]];
	{groups, groupsize}
]


(* ::Subsection::Closed:: *)
(*Wrap*)


SyntaxInformation[Wrap] = {"ArgumentsPattern" -> {_}};

Wrap[dat_]:= Mod[dat + Pi, 2 Pi] - Pi


(* ::Subsection:: *)
(*UnwrapDCT*)


(* ::Subsubsection::Closed:: *)
(*UnwrapDCT*)


SyntaxInformation[UnwrapDCT] = {"ArgumentsPattern" -> {_, _.}};

UnwrapDCT[psi_]:=UnwrapDCT[psi, None]

UnwrapDCT[psii_, wi_]:=Block[{
		psi, a, d, itt, w, alpha, rhoi,  norm, normi, phi, phii, 
		Qphii, maxi , i, soli, num, dena, denb
	},
	
	(*Phase unwrapping algorithem based on Ghiglia,Dennis C.,and Louis A.Romero. 10.1364/JOSAA.11.000107.*)
	
	(*prepare data*)
	psi = ToPackedArray@N@psii;
	a = ArrayDepth[psi];
	d = Dimensions[psi];
	
	(*make weights, w is min of weights in each direction eq 36 paper*)
	itt = If[wi===None, True, False];
	w = MakeWeights[psi, wi, d, a];
	
	(*initialize values*)
	rhoi = GetDifference[psi, w, True];(*should not be w*)
	
	(*no weigths do instant solve has weigts defined to itterative solver*)
	If[itt,
		(*instan solve*)
		SolvePoisson[rhoi],
		
		(*step 1: initialize parameters for loop*)
		i = 0;
		phi = 0.psi;
		norm = 10^-6 Norm@Flatten@rhoi;
		maxi = 100 (*Round[0.1 Times@@d]*);
		
		(*run loop*)
		(*If[a===3, PrintTemporary[Dynamic[i]," / ", maxi, "   ", norm, " < ", Dynamic[normi]]];*)
		
		While[True,(*should check for rhoi is all zero*)
			
			(*step 2: find solution calculate the phi update*)
			soli = SolvePoisson[rhoi];
			
			(*step 3: update k*)
			i += 1;
			
			(*step 4 or 5: define initial phi or update phi*)
			num = Total[rhoi soli, -1];
			phii = If[i===1, soli, soli + (num/denb) phii];
			
			(*store current value as i-1 value*)
			denb = num; If[denb===0., Break[]];
			
			(*step 6: perform one scalar and two vectors update*)
			Qphii = GetDifference[phii, w, False];
			dena = Total[phii Qphii, -1]; If[dena===0., Break[]];
			alpha = num/dena;
			rhoi -= alpha Qphii;
			phi += alpha phii;
						
			(*step 7: check for continue*)
			(*calculate norm*)
			normi = Norm@Flatten@rhoi;
			
			If[i > maxi || normi < norm, Break[]]
		];
		
		(*Print[i," / ",maxi,"   ", norm," < ",normi,"   "];*)
		
		phi
	]
]


(* ::Subsubsection::Closed:: *)
(*MakeWeights*)


MakeWeights[psi_, wi_, d_, a_]:=Block[{w},
	(*weighting for unwrapping which can be None, Automatic or predifined*)
	w = ToPackedArray@N@Switch[wi,
		None, ConstantArray[1., d],
		Automatic, 1. - Rescale[MedianFilter[N[StandardDeviationFilter[Sin[psi], 1] + StandardDeviationFilter[Cos[psi], 1]], 1]],
		_, If[d===Dimensions[wi], wi, Return[Message[UnwrapDCT::dim, d, Dimensions[w]]]]
	];
	
	(*calculates the min of w of paired voxels in all dimensions*)
	MinAt[w(*^2*), #]& /@ Range[a]
]


MinAt[arr_, lev_]:=Block[{k},
	k = UnitStep[DifferenceAt[arr, lev]];
	ToPackedArray@N@Switch[lev,
		1, k arr[[2;;]] + (1-k) arr[[;;-2]],
		2, k arr[[All, 2;;]] + (1-k) arr[[All, ;;-2]],
		3, k arr[[All, All, 2;;]] + (1-k) arr[[All, All, ;;-2]]
	]
]


(* ::Subsubsection::Closed:: *)
(*GetDifference*)


GetDifference[psi_, w_, True]:=GetDifference[psi, w, Wrap[#]&]

GetDifference[psi_, w_, False]:=GetDifference[psi, w, #&]

GetDifference[psi_, w_, wrap_]:=Block[{wrapF, pad, a, out},
	(*calculate the difference in all dimensions apply wrapping and weigtheing and calculate rho*)
	a = ArrayDepth[psi];
	pad = {{{1,1},0,0}, {0,{1,1},0}, {0,0,{1,1}}};
	out = w(wrap[DifferenceAt[psi, #]&/@Range[a]]);
	Total[DifferenceAt[ArrayPad[out[[#]], pad[[#,;;a]], 0.], #]&/@Range[a]]
]


DifferenceAt[arr_, lv_]:= -RotateDimensionsRight[Differences[RotateDimensionsLeft[arr, lv-1]], lv-1]


(* ::Subsubsection::Closed:: *)
(*SolvePoisson*)


SolvePoisson[rho_]:=SolvePoisson[rho, ArrayDepth[rho], Dimensions[rho]]

SolvePoisson[rho_, a_, d_]:=Block[{dctRho, dev, dctPhi},
	(* solve the poisson equation using DCT cash the divisor and handle the /0 for first index*)
	dctRho = FourierDCT[rho];
	dev = GetDev[a, d];
	Switch[a, 1, dev[[1]] = 1., 2, dev[[1,1]] = 1., 3, dev[[1,1,1]] = 1.];
	dctPhi = dctRho / dev;
	Switch[a, 1, dctPhi[[1]] = 0., 2, dctPhi[[1,1]] = 0., 3, dctPhi[[1,1,1]] = 0.];
	FourierDCT[dctPhi, 3]
]


GetDev[a_, d_] := GetDev[a, d] = a (Total@Cos[Pi RotateDimensionsRight[N[Array[{##}&, d]] - 1.] / N[d]] - a);


(* ::Subsection:: *)
(*DixonPhase*)


Options[DixonPhase] = {
   DixonBipolar -> True,
   DixonPrecessions -> -1,
   DixonFieldStrength -> 3,
   DixonNucleus -> "1H",
   DixonFrequencies -> {{0}, {3.8, 3.4, 3.1, 2.7, 2.5, 1.95, 0.5, -0.5, -0.6}},
   DixonAmplitudes -> {{1}, {0.088, 0.628, 0.059, 0.064, 0.059, 0.01, 0.039, 0.01, 0.042}},
   MonitorCalc->False
   };

SyntaxInformation[DixonPhase] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

DixonPhase[{real_, imag_}, echos_, OptionsPattern[]] := Block[{
	freqs, amps, iop, A, Ah, Ai, e1, e2, de, dt, hz, bip, mat, msk, comp, compi, ph, ph1, ph0, phi, 
	ph1i, ph0i, norm, i, itt, l, hl, sw, sw1, sw0, f0, f1, normN, UF, ni, n, t2s, Ac},
		
	(*initial phase 10.1016/J.MRI.2010.08.011*)
	(*bipolar phase 10.1002/mrm.24657*)
	
	(*define the water and fat frequencies and amplitudes*)
	freqs = OptionValue[DixonPrecessions] OptionValue[DixonFieldStrength] OptionValue[DixonFrequencies] GyromagneticRatio[OptionValue[DixonNucleus]];
	amps = #/Total[#] & /@ OptionValue[DixonAmplitudes];
	iop = 1/Abs[Total[(amps[[2]]^2) freqs[[2]]]/Total[amps[[2]]^2]];
	
	(*define the solution matrix*)
	A = (Total /@ (amps Exp[freqs (2 Pi I) #])) & /@ echos;
	Ah = ConjugateTranspose[A];
	Ai = Inverse[Re[Ah . A]];
	
	(*get the two first inphase locations and calculate timing*)
	{e1, e2} = FindInPhaseEchos[echos, iop, DixonBipolar -> False];
	de = e2 - e1;
	dt = (echos[[e2]] - echos[[e1]]);
	hz = {Abs[1/dt], 0.5, 0.25};
	
	(*make the time matrix*)
	bip = (-1)^Range[Length[echos]];
	mat = -I Transpose[{echos, Abs[bip], bip}];
	
	(*prepare signals for *)
	comp = compi = Transpose[real + imag I];
	msk = Unitize@Total@Abs@compi;
	ph = ph1 = ph0 = phi = ph1i = ph0i = 0 Re@First@compi;
	
	(*prepare itterative optimization*)
	norm = {{1., 1., 1.}};
	i = 0;
	itt = 50;
	l = Length[echos];
	hl = Round[.5 l];
	sw = sw1 = sw0 = f0 = f1 = False;
	
	t1=t2=t3=0.;
	If[OptionValue[MonitorCalc], PrintTemporary[Dynamic[{i, Last@Transpose[100 Transpose[norm]/norm[[1]]],{t1,t2,t3}}]]];
		
	(*start optimization*)
	Do[
		If[i =!= 0, AppendTo[norm, norm[[-1]]]];
		i++;
		
		t1=First@AbsoluteTiming[
		(*b0 phase*)
		phi = Mean[(UnwrapDCT /@ (Arg[compi[[# + de]]] - Arg[compi[[#]]])) & /@ Range[e1, Min[{e1 + hl, l - de}], Min[{hl,3}]]];
		ph += msk phi;
		compi = ApplyPhase[comp, hz {ph, ph1, ph0}, mat];
		
		(*clac norm*)
		norm[[-1, 1]] = Norm@Flatten@Pick[phi, msk, 1];
		normN = 100 (#/norm[[1]]&/@norm);
		];
		
		t3=First@AbsoluteTiming[
		(*biplolar phase*)
		If[normN[[-1, 3]] >= .1,
			UF = If[normN[[-1, 3]] < 25 || sw0, sw0 = True; # &, UnwrapDCT];
			ph0i = Mean[(bip[[#]] UF[Arg@DevideNoZero[compi[[#]]^2, compi[[# - 1]] compi[[# + 1]], "Comp"]] & /@ Range[2, l-1, 3])];
						
			ph0 += msk ph0i;
			compi = ApplyPhase[comp, hz {ph, ph1, ph0}, mat];
			
			(*clac norm*)
			norm[[-1, 3]] = Norm@Flatten@Pick[ph0i, msk, 1],
			norm[[-1, 3]] = norm[[-2, 3]]
		];
		(*clac norm*)
		normN = 100 (#/norm[[1]]&/@norm);
		];
		
		
		t2=First@AbsoluteTiming[
		(*initial phase*)
		Ac = Ah.compi;
		ph1i = UnwrapDCT /@ Arg[DotAc[RotateDimensionsLeft[Ac], RotateDimensionsLeft[Ai . Ac]]];
		ph1 += msk ph1i;
		compi = ApplyPhase[comp, hz {ph, ph1, ph0}, mat];
		
		(*clac norm*)
		norm[[-1, 2]] = Norm@Flatten@Pick[ph1i, msk, 1];
		normN = 100 (#/norm[[1]]&/@norm);
		];
		
		If[AllTrue[Last[normN], # < 4 &], Break[]]
	, {itt}];
		
	(*fix the initial phase*)
	phi = Mean[(UnwrapDCT[(Arg[compi[[# + de]]] - Arg[compi[[#]]])]) & /@ Range[e1, Min[{e1 + hl, l - de}], Min[{hl,3}]]];
	ph += msk phi;
	compi = ApplyPhase[comp, hz {ph, ph1, ph0}, mat];
		
	ph1i = UnwrapDCT[Arg@compi[[e2]]];
	ph1 += msk ph1i;
	compi = ApplyPhase[comp, hz {ph, ph1, ph0}, mat];

	(*get R2 star*)
	n = First@ FirstPosition[echos, First[Select[echos, # > 1.5 iop &]]];
	t2s = Last@T2Fit[Abs[Transpose[comp[[n ;;]]]], echos[[n ;;]]];
	
	(*give output*)
	{ph, ph1, ph0} = hz {ph, ph1, ph0} / (2 Pi);
	{{ph, t2s, ph1, ph0}, {e1, e2, n}}
]


ApplyPhase = Compile[{{comp, _Complex, 4}, {phase, _Complex, 4}, {mat, _Complex, 2}},
	comp Exp[mat . phase],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


DotAc = Compile[{{v1, _Complex, 1}, {v2, _Complex, 1}}, 
	v1 . v2,
RuntimeOptions -> "Speed", Parallelization -> True, RuntimeAttributes -> {Listable}];

FitBipolar[ph0_, msk_] := Block[{m, ydat, xdat, dat, fit, vals},
	m = UnitStep[Rescale[Total@Total@msk] - 0.5];
	ydat = MeanNoZero@MeanNoZero[ph0];
	xdat = Range[Length[ydat]];
	dat = Pick[Transpose[{xdat, ydat}], m, 1];
	fit = Fit[dat, {1, x}, x];
	vals = fit /. x -> xdat;
	msk ConstantArray[vals, Dimensions[msk][[1 ;; 2]]]
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
