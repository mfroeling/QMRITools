(* ::Package:: *)

(* ::Title:: *)
(*QMRITools DenoiseTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`DenoiseTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`DenoiseTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DeNoise::usage = 
"DeNoise[data,sigma,filtersize] removes Rician noise with standard deviation \"sigma\" from the given dataset using a kernel with size \"filtersize\" a gaussian kernel.
DeNoise[data,sigma,filtersize, Kernel->\"kerneltype\"] removes Rician noise with standard deviation \"sigma\" from the given dataset using a kernel with size \"filtersize\" and type \"kerneltype\".

Output is data denoised.

DeNoise[] is based on DOI: 10.1109/TMI.2008.920609"


PCADeNoise::usage = 
"PCADeNoise[data] removes rician noise from the data with PCA.
PCADeNoise[data, mask] removes rician noise from the data with PCA only withing the mask.
PCADeNoise[data, mask, sig] removes rician noise from the data with PCA only withing the mask using sig as prior knowledge or fixed value.

Output is de {data denoise, sigma map} by default if PCAOutput is Full then fitted {data dnoise , {sigma fit, average sigma}, {number components, number of fitted voxesl, number of max fits}, total fit -time per 500 ittt}.

PCADeNoise[] is based on DOI: 10.1016/j.neuroimage.2016.08.016 and 10.1002/mrm.26059"

PCAFitHist::usage =  
"PCAFitHist[data] fits the marchencopasteur distribution to the PCA of the data using hist fit.
PCAFitHist[data, sig] fits the marchencopasteur distribution to the PCA of the data using sig as start value or fixed value using hist fit.

Output is {simga, number of noise comp, and denoised matrix, itterations}."

PCAFitEq::usage = 
"PCAFitEq[data] fits the marchencopasteur distribution to the PCA of the data using grid search.
PCAFitEq[data, sig] fits the marchencopasteur distribution to the PCA of the data using sig as start value or fixed value using grid search.

Output is {simga, number of noise comp, and denoised matrix}."


AnisoFilterTensor::usage = 
"AnisoFilterTensor[tens, diffdata] Filter the tensor tens using an anisotropic diffusion filter (Perona-Malik). 
It uses the diffusion weighted data diffdata to find edges that are not visible in the tensor.
Edge weights based on the diffusion data are averaged over all normalized diffusion direction.

Output is the smoothed tensor.

AnisoFilterTensor[] is based on DOI: 10.1109/ISBI.2006.1624856"

WeightMapCalc::usage =  
"WeightMapCalc[diffdata] calculates a weight map which is used in AnisoFilterTensor.

Output is a weight map of the diffdata which is high in isotropic regions and low at edges."


(* ::Subsection::Closed:: *)
(*Options*)


DeNoiseKernel::usage = 
"DeNoiseKernel is and option for DeNoise. Values can be \"Disk\", \"Box\" or \"Gaussian\"."

DeNoiseMonitor::usage = 
"DeNoiseMonitor is and option for DeNoise. Monitor the denoising progres."

DeNoiseIterations::usage = 
"DeNoiseIterations is and option for DeNoise. Specifies the number of the denoising iterations."


PlotSolution::usage = 
"PlotSolution is an option for PCAFitHist, if set true it dispays the fitting itterations."

FitSigma::usage = 
"FitSigma is an option of PCAFitHist, PCAFitEq and PCADeNoise, if set True sig is fitted if set False sigma is fixed to input value."

PCAFitParameters::usage = 
"PCAFitParameters is an option of PCADeNoise and PCAFitHist. {nb, pi, maxit} = bins, initial signal components, maximum number of itterations."

PCAKernel::usage = 
"PCAKernel is an option of PCADeNoise. It sets the kernel size."

PCAOutput::usage = 
"PCAOutput is an option of PCADeNoise. If output is full the output is {datao, {output[[1]], sigmat}, {output[[2]], output[[3]], j}, timetot}.
Else the output is {datao, sigmat}."

PCATollerance::usage = 
"PCATollerance is an option of PCADeNoise and shuld be an integer > 0. Default value is 0. When increased the denoise method removes less noise."

PCAWeighting::usage = 
"PCAWeighting is an option of PCADeNoise and can be True of False. Default value is False. When True the weights of the per voxel result are calculated based on the number of non noise components."

PCAClipping::usage = 
"PCAClipping is an option of PCADeNoise and can be True of False. If True the output is clipped between 0 and the max absolute value of the input data."

AnisoStepTime::usage =
"AnisoStepTime is an option for AnisoFilterTensor and defines the diffusion time, when small more step are needed."

AnisoFilterSteps::usage =
"AnisoFilterSteps is an option for AnisoFilterTensor and defines the amoutn of diffusin steps taken. Higher is more smoothing"

AnisoWeightType::usage =
"AnisoWeightType is an option for AnisoFilterTensor and WeightMapCalc and defines the weighting, eigher 1, the exponent of (-g/kappa) or 2, 1/(1+g/kappa)."

AnisoKappa::usage =
"AnisoKappa is an option for AnisoFilterTensor and WeightMapCalc and defines the weighting strenght, all data is normalize to 100 before filetering."


(* ::Subsection::Closed:: *)
(*Error Messages*)


DeNoise::data =
"Error: not able to proces the combination of this data set (Dimensions `2` ) with the given size kernel `1`, posibilities:
- 2D data, 2D kernel, sigma = single value
- 3D data, 2D kernel, sigma = single value or list with value specified for each slice
- 3D data, 3D kernel, sigma = single value
- 4D data, 2D kernel, sigma = single value or 2D array with value specified for each slice and diffusion direction
- 4D data, 3D kernel, sigma = single value or list with value specified for each diffusion direction"

DeNoise::filt =
"Error: The dimension of the kernel (`1`D) is of higher order than the dimension of the dataset (`2`D)."

DeNoise::dim =
"Error: The dimension of the sigmap (`1`) is not the same as the dimension of the dataset (`2`)."

DeNoise::kern = 
"Error: Unknown kernel type:`1`, use \"Gaussian\", \"Box\" or \"Disk\"."

DeNoise::sig =
"Data and simga are of unequal dimensions. Data: `1`, Sigma: `2`."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Denoise*)


(* ::Subsubsection::Closed:: *)
(*Denoise*)


Options[DeNoise] = {DeNoiseKernel -> "Gaussian", DeNoiseMonitor -> False, DeNoiseIterations -> 1};

SyntaxInformation[DeNoise] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

DeNoise[dat_, sigi_, filt_, OptionsPattern[]] := 
 Module[{kern, out, type, dimsig, dimdat, data, sig},
  sig = N[sigi];
  data = ToPackedArray@N@dat;
  (*Check dimensions,must be of lower order than data*)
  If[
   ArrayQ[sig] && (ArrayDepth[data] == 3 || ArrayDepth[data] == 4),
   dimsig = Dimensions[sig];
   dimdat = Dimensions[data];
   If[
    (ArrayDepth[data] == 3 && dimsig != dimdat) || 
     ArrayDepth[data] == 4 && dimsig != dimdat[[{1, 3, 4}]], 
    Return[Message[DeNoise::dim, dimsig, dimdat]]
    ];
   ];
  
  (*Check dimensions,filter must be of lower order than data*)
  If[Length[filt] > ArrayDepth[data], 
   Return[Message[DeNoise::filt, Length[filt], ArrayDepth[data]]]];
  
  (*Create filter*)
  kern = Switch[type = OptionValue[DeNoiseKernel],
    "Box", kern = BoxMatrix[filt]/Total[Flatten[BoxMatrix[filt]]],
    "Disk", kern = DiskMatrix[filt]/Total[Flatten[DiskMatrix[filt]]],
    "Gaussian", kern = GaussianMatrix[{filt}],
    _, Return[Message[DeNoise::kern, type]];];
  
  If[OptionValue[DeNoiseMonitor], 
   PrintTemporary["Using " <> type <> " kernel."]];
  
  out = ToPackedArray@N@Nest[DeNoisei[#, sig, filt, kern] &, data, OptionValue[DeNoiseIterations]];
  
  If[
   ArrayQ[out], Return[Clip[out, {0.9 Min[data], 1.1 Max[data]}]], 
   Return[Message[DeNoise::data, Dimensions[sig], Dimensions[data]]]
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*DeNoisei*)


DeNoisei[data_?MatrixQ,sig_,filt_,kern_?MatrixQ]:=
Module[{},
	(*PrintTemporary["Appling 2D kernel on 2D dataset."];*)
	DeNoiseApp[data,sig,filt,kern]]

DeNoisei[data:{_?MatrixQ..},sig_,filt_,kern:{_?MatrixQ..}]:=
Module[{},
	(*PrintTemporary["Appling 3D kernel on 3D dataset."];*)
	DeNoiseApp[data,sig,filt,kern]]

DeNoisei[data:{_?MatrixQ..},sig_,filt_,kern_?MatrixQ]:=
Module[{},
	(*PrintTemporary["Appling 2D kernel on 3D dataset."];*)
	(*Using one sigma value for all slices."];*)
	If[NumberQ[sig],
		Map[DeNoiseApp[#,sig,filt,kern]&,data],
		MapThread[DeNoiseApp[#1,#2,filt,kern]&,{data,sig}]
		]
	]

DeNoisei[data:{{_?MatrixQ..}..},sig_,filt_,kern_?MatrixQ]:=
Module[{},
	(*PrintTemporary["Appling 2D kernel on 4D dataset."];*)
	(*Using one sigma value for all slices and directions."];*)
	If[NumberQ[sig],
		Map[DeNoiseApp[#,sig,filt,kern]&,data,{2}],
		Transpose[Map[MapThread[DeNoiseApp[#1,#2,filt,kern]&,{#,sig}]&,Transpose[data]]]
		]
	]

DeNoisei[data:{{_?MatrixQ..}..},sig_,filt_,kern:{_?MatrixQ..}]:=
Module[{},
	(*PrintTemporary["Appling 3D kernel on 4D dataset."];*)
	(*Using one sigma value for all diffusion directions."];*)
	Transpose[Map[DeNoiseApp[#1,sig,filt,kern]&,Transpose[data]]]
	]


(* ::Subsubsection::Closed:: *)
(*DeNoiseApp*)


DeNoiseApp[data_, sig_, filt_, kern_] := 
 Module[{secmod, quadmod},
  secmod = ListConvolve[kern, data^2, Transpose[{filt + 1, -(filt + 1)}], 0.];
  quadmod = ListConvolve[kern, data^4, Transpose[{filt + 1, -(filt + 1)}], 0.];
  If[NumberQ[sig],
  	NoiseAppCN[secmod, quadmod, data, sig],
  	NoiseAppC[secmod, quadmod, data, sig]
  ]
  ]

NoiseAppCN = Compile[{{secmod, _Real, 2}, {quadmod, _Real, 2}, {data, _Real, 2}, {sig, _Real, 0}},
   Block[{top, div, K, deb},
    top = (4 sig^2 (secmod - sig^2));
    div = (quadmod - secmod^2) + 10^-10;
    K = (1 - top/div);
    deb = Sqrt[Clip[(secmod - 2 sig^2) + (K (data^2 - secmod)), {0., Infinity}]]
    ]];

NoiseAppC = Compile[{{secmod, _Real, 2}, {quadmod, _Real, 2}, {data, _Real, 2}, {sig, _Real, 2}},
   Block[{top, div, K, deb},
    top = (4 sig^2 (secmod - sig^2));
    div = (quadmod - secmod^2) + 10^-10;
    K = (1 - top/div);
    deb = Sqrt[Clip[(secmod - 2 sig^2) + (K (data^2 - secmod)), {0., Infinity}]]
    ]];

NoiseAppCN = Compile[{{secmod, _Real, 3}, {quadmod, _Real, 3}, {data, _Real, 3}, {sig, _Real, 0}},
   Block[{top, div, K, deb},
    top = (4 sig^2 (secmod - sig^2));
    div = (quadmod - secmod^2) + 10^-10;
    K = (1 - top/div);
    deb = Sqrt[Clip[(secmod - 2 sig^2) + (K (data^2 - secmod)), {0., Infinity}]]
    ]];
    
NoiseAppC = Compile[{{secmod, _Real, 3}, {quadmod, _Real, 3}, {data, _Real, 3}, {sig, _Real, 3}},
   Block[{top, div, K, deb},
    top = (4 sig^2 (secmod - sig^2));
    div = (quadmod - secmod^2) + 10^-10;
    K = (1 - top/div);
    deb = Sqrt[Clip[(secmod - 2 sig^2) + (K (data^2 - secmod)), {0., Infinity}]]
    ]];



(* ::Subsection:: *)
(*PCADenoise*)


(* ::Subsubsection::Closed:: *)
(*PCADeNoise*)


Options[PCADeNoise] = {PCAKernel -> 5, PCAFitParameters -> {10, 6, 10}, FitSigma -> False, PCAOutput -> Full, Method->"Equation", PCATollerance->0, PCAWeighting->True, PCAClipping->True};

SyntaxInformation[PCADeNoise] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

PCADeNoise[data_, opts : OptionsPattern[]] := PCADeNoise[data, 1, 0., opts];

PCADeNoise[data_, mask_, opts : OptionsPattern[]] := PCADeNoise[data, mask, 0., opts];

PCADeNoise[datai_, maski_, sigmai_, OptionsPattern[]] := Block[{
	data, mask, maskd, sigm, ker, off, datao, weights, sigmat, dim, zdim,
	ydim, xdim, ddim, nb, pi, maxit, output, time1, time, weigth,
	timetot, sigi, sigf, zm, zp, xm, xp, ym, yp, fitdata, filt, sigo,
	Nes, datn, it, m, n, fitsig, step, totalItt, j, max,
	start, maxIttN, tol
	},
	
	(*make everything numerical to speed up*)
	data = Transpose[ToPackedArray[N@datai]];
	maskd = Unitize@Total@data;
	data = TransData[data,"l"];
	mask = ToPackedArray[N@maski];
	sigm = ToPackedArray[N@sigmai];
	
	(*get data dimensions*)
	dim = {zdim, ydim, xdim, ddim} = Dimensions[data];
	
	(*get options for fit*)
	(*initial compoenents and itterations*)
	{nb,pi,maxit} = OptionValue[PCAFitParameters];
	nb = If[NumberQ[nb], nb, 5];
	pi = If[NumberQ[pi], pi, 10];
	maxit = If[NumberQ[maxit], maxit, 10];
	(*kernel size*)
	ker = OptionValue[PCAKernel];
	ker = If[EvenQ[ker], ker - 1, ker];
	(*fit sigma*)
	fitsig = OptionValue[FitSigma];
	(*tollerane if >0 more noise components are kept*)
	tol=OptionValue[PCATollerance];
	
	(*define runtime parameters*)
	{m, n} = MinMax[{ddim, ker^3}];
	off = Round[(ker - 1)/2];
	datao = ConstantArray[0., dim];
	weights = sigmat = datao[[All, All, All, 1]];
	sigm = If[NumberQ[sigm], ConstantArray[sigm, {zdim, ydim, xdim}], sigm];
	
	(*if mask is a number make it 1 for all voxels*)
	mask = If[NumberQ[mask], maskd weights + 1, maskd mask];
	
	(*parameters for monitor*)
	start = off + 1;
	totalItt = Total[Flatten[mask[[start;;zdim-off, start;;ydim-off, start;;xdim-off]]]];
	j=0;
	
	Monitor[output = Table[
		(*Check if masked voxel*)
		If[mask[[z, y, x]] == 0.||AllTrue[data[[z, y, x]], # === 0. &],
			(*background do nothing*)
			{0., 0., 0.}
			,
			j++;
			(*define initial sigma*)
			sigi = sigm[[z, y, x]];
			(*get pixel range and data*)
			{{zm, ym, xm}, {zp, yp, xp}} = {{z, y, x} - off, {z, y, x} + off};
			(*kernel box to vector*)
			fitdata = Flatten[data[[zm ;; zp, ym ;; yp, xm ;; xp]], 2];
			
			(*perform the fit and reconstruct the noise free data*)
			Switch[OptionValue[Method],
				"Equation",
				sigi = If[fitsig, 0., sigi];
				{sigo, Nes, datn} = PCAFitEqi[fitdata, {m, n}, sigi , tol];
				it = 1;
				, _,
				{sigo, Nes, datn, it} = PCAFitHisti[fitdata, {m, n}, sigi, tol, FitSigma -> fitsig, PCAFitParameters -> {nb, pi, maxit}];
			];
			
			(*reshape the vector into kernel box*)
			datn = Fold[Partition, datn, {ker, ker}];
			
			(*collect the noise free data and weighting matrix*)
			(*weight the signal for number of components*)
			weigth = If[OptionValue[PCAWeighting], 1. / (m-Nes+1), 1.];
			
			(*sum data and sigma and weight for numer of components*)
			datao[[zm ;; zp, ym ;; yp, xm ;; xp, All]] += (weigth datn);
			sigmat[[zm ;; zp, ym ;; yp, xm ;; xp]] += weigth sigo;
			
			(*count the total weighting*)
			weights[[zm ;; zp, ym ;; yp, xm ;; xp]] += weigth;
			
			(*output sig,Nest and itterations*)
			{sigo, Nes, it}
		], {z, start, zdim - off}, {y, start, ydim - off}, {x, start, xdim - off}];
		,
		(*monitor*)
		ProgressIndicator[j, {0, totalItt}]
	];
	
	max = 1.1 Max[Abs[data]];
	
	(*correct output data for weightings*)
	datao = Transpose@TransData[DevideNoZero[datao, weights],"r"];
	If[OptionValue[PCAClipping], datao = Clip[datao,{0, max}]];
	
	sigmat = DevideNoZero[sigmat, weights];
	output = ArrayPad[#, off] & /@ TransData[output, "r"];
	
	(*define output*)
	If[OptionValue[PCAOutput] === Full,
		(*fitted dta,average sigma,{sigma fit,number components, number of fitted voxesl,number of max fits}*)
		{ToPackedArray@N@datao, sigmat, output},
		{ToPackedArray@N@datao, sigmat}
	]
]


(* ::Subsubsection::Closed:: *)
(*MarchenkoPasturC*)


(*compiled marchenco pastur distribution function*)
MarchenkoPasturC[lab_,Q_?NumericQ,sig_?NumericQ]:=MarchenkoPasturCi[lab,Q,sig];
MarchenkoPasturCi=Compile[{{lab,_Real,0},{Q,_Real,0},{sig,_Real,0}},
	Block[{labm,sig2=sig^2,labp,Qs=Sqrt[Q]},
		(*define parameters for function*)
		labm=sig2 (1-Qs)^2;
		labp=sig2 (1+Qs)^2;
		(*define piecewise function*)
		Piecewise[{{Sqrt[(labp-lab) (lab-labm)]/(2 Pi sig2 Q lab),labm<lab<labp}},0]
],RuntimeAttributes->{Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SVD*)


(*singular ValueDecomposition of matrix and eigenval normalisation*)
(*determine dimension,make sure that m<n*)
SVD[mat_] := SVD[mat, MinMax[Dimensions[mat]]];
(*if m and n are known*)
SVD[mat_, {m_, n_}] := Block[{u, w, v, eig},
	(*no need for transpose,eig of mat and mat` are equal*)
	(*perform singular value decomposition*)
	{u, w, v} = SingularValueDecomposition[mat];
	(*normalize eigenvalues from SVD*)
	eig = Diagonal[w]^2/n;
	{u, w, Transpose[v], eig, m, n}
]


(* ::Subsubsection::Closed:: *)
(*CalcSigFunc*)


(*Fit sig on MarchenkoPatur distribution*)
CalcSigFunc[dat_, Q_, sigi_] := Block[{sig, lab},
   (Quiet@FindFit[dat, MarchenkoPasturC[lab, Q, sig], {{sig, sigi}}, lab,
       AccuracyGoal -> 5, PrecisionGoal -> 5,Method -> "LevenbergMarquardt", MaxIterations -> 50]
       )[[1, 2]]];


(* ::Subsubsection::Closed:: *)
(*HistListC*)


(*Make Histogram bins, is fast then using normal mathematica function*)
HistListC = Compile[{{dat, _Real, 1}, {nbins, _Integer, 0}}, Block[
    {min, max, maxmin, binw, bins, cpos, binst, comp, ydat, xdat, tall, miss},
    (*get histogram range*)
    min = Min[dat]; max = Max[dat];
    {min, max} = Chop[{min, max}];
    (*check if range > 0*)
    If[min === max,
     (*no range make delta function*)
     binw = 1/nbins;
     cpos = Ceiling[(nbins + 0.001)/2];
     bins = ConstantArray[0, nbins];
     bins[[cpos]] = nbins;
     Transpose[{Range[-0.5 + binw/2, 0.5, binw], bins}]
     ,
     (*range is greate then zero, thus calculate the bin width *)
     maxmin = max - min;
     If[maxmin<=0.,Print[{"error",{min,max}}];Print[dat]];
     binw = maxmin/nbins;
     (*count in number of bins*)
     tall = Tally[Floor[nbins (dat - min)/(1.001 maxmin)]];
     (*check if bins are empty*)
     comp = Complement[Range[0, nbins], tall[[All, 1]]];
     If[comp[[1]] != nbins,
      (*fil in missing bins*)
      miss = {#, 0} & /@ comp[[;; -2]];
      binst = Sort[Join[tall, miss]],
      (*no bins missing*)
      binst = Sort[tall]
      ];
     (*calcualte the y and x data*)
     bins = binst[[All, 2]];
     ydat = N@bins/Total[bins]/binw;
     xdat = Range[min + binw/2, max, binw];
     (*output the data*)
     Transpose[{xdat, ydat}]
     ]], 
	RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*ErrorFunc*)


(*Errorfunction only used for plotting the minimization*)
ErrorFunc[data_, Q_, sig_] := Block[{xdata, ydata, vals, tvals},
	xdata = data[[1]];
	ydata = data[[2]];
	vals = MarchenkoPasturC[xdata, Q, sig];
	Total[(ydata - vals)^2]
];


(* ::Subsubsection::Closed:: *)
(*PCAFitHist*)


Options[PCAFitHist] = {PlotSolution -> False, FitSigma -> True, PCAFitParameters -> {10, 6, 10}};

SyntaxInformation[PCAFitHist] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

(*no initial sigma given*)
PCAFitHist[data_, {m_, n_}, opts : OptionsPattern[]] := PCAFitHisti[data, {m, n}, 0., 0, opts]
(*no initial sigma given*)
PCAFitHist[data_, {m_, n_}, sigi_, opts : OptionsPattern[]] := PCAFitHisti[data, {m, n}, sigi, 0, opts]
(*initial sigma is given*)
PCAFitHist[data_, {m_, n_}, sigi_, toli_, opts :OptionsPattern[]] := PCAFitHisti[data, {m, n}, sigi, toli, opts]

(*internal function*)
Options[PCAFitHisti] = Options[PCAFitHist];

PCAFitHisti[data_, {mi_, ni_}, sigii_, toli_, OptionsPattern[]] := Block[
  {nb, pi, maxit, u, w, v, eig, m, n, i, pi1, pi0, Nes, Q, Qs, sigi, 
   sig, hlist, eigp, tol},
  (*get options,number of bins,initial p and max itterations*)
  {nb, pi, maxit} = OptionValue[PCAFitParameters];
  (*perform svd*)
  {u, w, v, eig, m, n} = SVD[data, {mi, ni}];
  (*perform heuristic ittarative fitting*)
  i = pi1 = pi0 = 0;
  Do[(*count fit and how often max fit*)
   i++;
   (*number of nois comp,Q and Qs*)
   Nes = (m - pi);
   Q = N[Nes/n];
   Qs = Sqrt[Q];
   (*custom histogram list function for speed*)
   hlist = HistListC[eig[[pi + 1 ;;]], nb];
   
   (*calcualte initial sig from data or use input for i=1*)
   sigi = 
    If[i == 1, 
     If[sigii === 0., Sqrt[Last[eig]]/Sqrt[(1 - Qs)^2], sigii], sig];
   (*perform the fit,data from histogramlist*)
   (*fit MP function to data,returns sig if fitsimgam is true, if sigma is fixed no fit*)
   sig = If[OptionValue[FitSigma], CalcSigFunc[hlist, Q, sigi], sigi];
   
   If[sig < 0,
    (*if sig becomes negative increase the number of bins reset and try again*)
    nb += 2; i = pi1 = pi0 = 0;
    ,
    (*determine number of noise components with given sig*)
    eigp = sig^2 (1 + Qs)^2;
    pi1 = Clip[Length[Select[eig, # > eigp &]], {0, m - 1}];
    (*this ends if the same solution or the same solution as the previous itteration is found*)
    If[pi == pi1 || pi1 == pi0, Break[]];
    (*updata pi values*)
    {pi0, pi} = {pi, pi1};
    ];
   , {maxit}];
   
   (*constartin pi plus tol*)
  tol=Round[Clip[pi+toli,{0,n}]];
  (*set the noise components to zero*)
  w[[tol ;;, tol ;;]] = 0.;
  (*give output, number of noise comp and sigma and number of itterations*)
  {sig, Nes, u.w.v, i}
  ]


(* ::Subsubsection::Closed:: *)
(*PCAFitEq*)


(*PCAfit using set of equations*)
SyntaxInformation[PCAFitEq]={"ArgumentsPattern"->{_,_,_.,_.}};

(*no initial sigma given*)
PCAFitEq[data_, {m_, n_}] := PCAFitEqi[data, {m, n}, 0., 0]
(*no initial normal tolarance*)
PCAFitEq[data_, {m_, n_}, sigi_] := PCAFitEqi[data, {m, n}, sigi, 0]
(*initial sigma is given*)
PCAFitEq[data_, {m_, n_}, sigi_, toli_] :=  PCAFitEqi[data, {m, n}, sigi, toli]

(*internal function*)
PCAFitEqi[data_, {mi_, ni_}, sigi_, toli_] := Block[{u, w, v, m, n, eig, pi, sig,tol},
	(*perform svd*)
	{u, w, v, eig, m, n} = SVD[data, {mi, ni}];
	(*if sigma is given perform with fixed sigma,else fit both*)
	{pi, sig} = If[N[sigi] != 0.,
		GridSearchSig[eig, m, n, sigi],
		GridSearch[eig, m, n]
	];
	
	(*constartin pi plus tol*)
	tol=Round[Clip[pi+toli,{0,n}]];
	(*set the noise components to zero*)
	w[[tol ;;, tol ;;]] = 0.;
	(*give output,simga,number of noise comp,and denoised matrix*)
	{sig, m - (tol), u.w.v}
]


(* ::Subsubsection::Closed:: *)
(*GridSearch*)


(*gird search to find p at which sig is almost equal*)
GridSearch = Compile[{{eig, _Real, 1}, {m, _Integer, 0}, {n, _Integer, 0}}, Block[{Nes, llab, eq1, eq2, pi, sig, eigl},
	(*initialize values*)
	eigl = Last[eig];
	pi = 1; eq1 = 0.; eq2 = 10.;

	(*find p for which eq1 and eq2 is equal to given sig*)
	While[eq2 > eq1 && pi < m,
		(*/Max[sig,1.0],sig is<1 per definition*)
		eq1 = (Mean[eig[[pi ;; m]]]);
		eq2 = ((eig[[pi]] - eigl)/(4 Sqrt[((m - pi)/n)]));
		pi++;
	];
	(*give output,number of noise comp and sigma*)
	{pi, Sqrt[(eq1 + eq2)/2]}],
	RuntimeOptions -> "Speed",  Parallelization -> True
];


(* ::Subsubsection::Closed:: *)
(*GridSearchSig*)


(*gird search to find p with a given sig,get mean p of both equations*)
GridSearchSig = Compile[{{eig, _Real, 1}, {m, _Integer, 0}, {n, _Integer, 0}, {sig, _Real, 0}},	Block[{pi, eq1, eq2, sigm, sig2, eigl},
	(*initialize values*)
	eigl = Last[eig];
	sig2 = sig^2;
	pi = 1; eq1 = 2 sig2; eq2 = 2 sig2;
	(*find p for which eq1 and eq2 is equal to given sig*)
	While[(eq1 - sig2 > 0 || eq2 - sig2 > 0) && pi < m,
		(*/Max[sig,1.0],sig is<1 per definition*)
		eq1 = (Mean[eig[[pi ;; m]]]);
		eq2 = ((eig[[pi]] - eigl)/(4 Sqrt[((m - pi)/n)]));
		pi++;
	];
	(*give output,number of noise comp and sigma*)
	{pi, sig}],
	RuntimeOptions -> "Speed",  Parallelization -> True
];


(* ::Subsection:: *)
(*AnisotropicFilterTensor*)


(* ::Subsubsection::Closed:: *)
(*AnisotropicFilterTensor*)


Options[AnisoFilterTensor] = {AnisoWeightType->2, AnisoKappa->5., AnisoStepTime->1, AnisoFilterSteps->5};

SyntaxInformation[AnisoFilterTensor] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

AnisoFilterTensor[tensi_,dat_,OptionsPattern[]]:=Block[{
weights,kernels,mn,j,datf,kers,wts,lambda,finDiff,wtsI,
itt,time,kappa,type, data, tens
},
(*get the options*)
itt=OptionValue[AnisoFilterSteps];
time=OptionValue[AnisoStepTime];
kappa=N@OptionValue[AnisoKappa];
type=Clip[Round@OptionValue[AnisoWeightType],{1,2}];

(*calculate the edges based on the diffusion images*)
PrintTemporary["Determaning the weights based on the data."];
data=ToPackedArray@N@dat;
tens=ToPackedArray@N@tensi;
weights=WeightMapCalc[data, AnisoKappa->kappa, AnisoWeightType->type];
(*get the fixed parameters*)
mn=Mean[tens[[1;;3]]];
{kers,wts}=KernelWeights[];
lambda=1/Length[kers];

(*filter the tensor*)
PrintTemporary["Anisotropic filtering of the tensor."];
j=0;PrintTemporary[ProgressIndicator[Dynamic[j],{0,itt 6}]];
Table[
(*Normalize the diffusion tensor*)
datf = 100 DevideNoZero[tens[[tt]],mn];
(*perform the diffusion smoothing itterations*)
Do[
j++;
finDiff=FinDiffCalc[datf,kers];
wtsI=weights WeightCalc[finDiff,wts,kappa,type];
datf=datf+time lambda Total@(wtsI finDiff);
,itt];
(*revert tensor normalization*)
datf=mn datf/100
,{tt,1,6}](*loop over tensor*)
]


(* ::Subsubsection::Closed:: *)
(*WeightMapCalc*)


Options[WeightMapCalc]={AnisoWeightType->2, AnisoKappa->10.};

SyntaxInformation[WeightMapCalc] = {"ArgumentsPattern" -> {_,  OptionsPattern[]}};

WeightMapCalc[data_,OptionsPattern[]]:=Block[{kers,wts,weights,finDiff,dat,dim,len, kappa, type },
(*get the options*)
kappa=N@OptionValue[AnisoKappa];
type=Clip[Round@OptionValue[AnisoWeightType],{1,2}];
(*get the kernerl and weights*)
{kers,wts}=KernelWeights[];
(*prepare output *)
dim=Dimensions[data];
len=dim[[2]];dim=Drop[dim,{2}];
weights=ConstantArray[0,Prepend[dim,Length[wts]]];
(*get the weighting for all diffusion images*)
i=0;PrintTemporary[ProgressIndicator[Dynamic[i],{0,len}]];
(
i++;
(*normalize the data*)
dat=100#/Max[Abs[#]];
(*add to the weights*)
weights+=WeightCalc[FinDiffCalc[dat,kers],wts,kappa,type];
)&/@Transpose[data];

(*normalize the weights between 0 and 1*)
(*weights=Mean[weights];
weights=weights/Max[weights];*)
(#/Max[#])&/@weights
]


(* ::Subsubsection::Closed:: *)
(*KernelWeights*)


KernelWeights[]:=Block[{cent,ker,keri,wtsi},
ker=ConstantArray[0,{3,3,3}];
ker[[2,2,2]]=-1;
cent={2,2,2};
Transpose[Flatten[Table[
If[{i,j,k}==cent,
Nothing,
keri=ker;keri[[i,j,k]]=1;
wtsi=N@Norm[cent-{i,j,k}];
{keri,1/wtsi^2}
],{i,1,3},{j,1,3},{k,1,3}],2]]
];


(* ::Subsubsection::Closed:: *)
(*WeightCalc*)


WeightCalc[finDiff_,wts_,kappa_,type_]:=wts Switch[type,1,Exp[-((finDiff/kappa)^2)],2,1./(1.+(finDiff/kappa)^2)];


(* ::Subsubsection::Closed:: *)
(*FinDiffCalc*)


FinDiffCalc[dat_,kers_]:=ParallelMap[ListConvolve[#,dat,{2,2,2},0]&,kers]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
