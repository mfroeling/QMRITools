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

Output is data denoised.

DeNoise[] is based on DOI: 10.1109/TMI.2008.920609."

PCADeNoise::usage = 
"PCADeNoise[data] removes rician noise from the data with PCA.
PCADeNoise[data, mask] removes rician noise from the data with PCA only withing the mask.
PCADeNoise[data, mask, sig] removes rician noise from the data with PCA only withing the mask using sig as prior knowledge or fixed value.

Output is de {data denoise, sigma map} by default if PCAOutput is Full then fitted {data dnoise , {sigma fit, average sigma}, {number components, number of fitted voxesl, number of max fits}, total fit -time per 500 ittt}.

PCADeNoise[] is based on DOI: 10.1016/j.neuroimage.2016.08.016 and 10.1002/mrm.26059."

NNDeNoise::usage = 
"NNDeNoise[data] removes rician noise from the data using self supersized neural net.
NNDeNoise[data, mask] removes rician noise from the data with PCA using self supersized neural net withing the mask.

PCADeNoise[] is based on DOI:10.48550/arXiv.2011.01355."

DenoiseCSIdata::usage = 
"DenoiseCSIdata[spectra] perfroms PCA denoising of the complex values spectra, data has to be 3D and the spectral dimensions is last, {x,y,z,spectra}."

DenoiseDynamicSpectraData::usage = 
"DenoiseDynamicSpectraData[spectra] perfroms PCA denoising of the complex values spectra, The data is given as a list of dynamicly acquired spectra {dynamic ,spectra}."


AnisoFilterTensor::usage = 
"AnisoFilterTensor[tens, diffdata] Filter the tensor tens using an anisotropic diffusion filter (Perona-Malik). 
It uses the diffusion weighted data diffdata to find edges that are not visible in the tensor.
Edge weights based on the diffusion data are averaged over all normalized diffusion direction.

AnisoFilterTensor[tens] Same but does not use the data for edge identification.

Output is the smoothed tensor.

AnisoFilterTensor[] is based on DOI: 10.1109/ISBI.2006.1624856."

WeightMapCalc::usage = 
"WeightMapCalc[diffdata] calculates a weight map which is used in AnisoFilterTensor.

Output is a weight map of the diffdata which is high in isotropic regions and low at edges."


AnisoFilterData::usage = 
"AnisoFilterData[data] Filter the diffusion tensor data using an anisotropic filter based on the strucure tensor of the data. 

Output is the smoothed data.

AnisoFilterData[] is based on DOI: 10.1016/j.jbiomech.2021.110540 and 10.1016/j.mri.2009.10.001 and 10.1371/journal.pone.0126953."


HarmonicDenoiseTensor::usage =
"HarmonicDenoiseTensor[tens, mask, vox] uses the harmonic denoising method to denoise the tensor within the mask.
Values for which the tensor is not defined but are within the mask will be filled in.
HarmonicDenoiseTensor[tens, seg, vox] will do the same for each segmentation label in seg.
HarmonicDenoiseTensor[tens, seg, vox, labs] will do the same for each segmentation label in seg that is specified in labs.

HarmonicDenoiseTensor[] is based on 10.1016/j.media.2011.01.005."


(* ::Subsection::Closed:: *)
(*Options*)


DeNoiseKernel::usage = 
"DeNoiseKernel is and option for DeNoise. Values can be \"Disk\", \"Box\" or \"Gaussian\"."

DeNoiseMonitor::usage = 
"DeNoiseMonitor is and option for DeNoise. Monitor the denoising progres."

DeNoiseIterations::usage = 
"DeNoiseIterations is and option for DeNoise. Specifies the number of the denoising iterations."


PCAKernel::usage = 
"PCAKernel is an option of PCADeNoise. It sets the kernel size."

PCAOutput::usage = 
"PCAOutput is an option of PCADeNoise. If output is full the output is {datao, {output[[1]], sigmat}, {output[[2]], output[[3]], j}, timetot}.
Else the output is {datao, sigmat}."

PCATolerance::usage = 
"PCATolerance is an option of PCADeNoise and shuld be an integer > 0. Default value is 0. When increased the denoise method removes less noise."

PCAWeighting::usage = 
"PCAWeighting is an option of PCADeNoise and can be True of False. Default value is False. When True the weights of the per voxel result are calculated based on the number of non noise components."

PCAClipping::usage = 
"PCAClipping is an option of PCADeNoise and can be True of False. If True the output is clipped between 0 and the max absolute value of the input data."

PCANoiseSigma::usage = 
"PCANoiseSigma is an option of DenoiseCSIdata and can be \"Corners\" or \"Automatic\"."

PCAComplex::usage = 
"PCAComplex is an option of PCADeNoise and can be True of False. If set true the input data is expexted to be {real, imag}."


NNThreshold::usage = 
"NNThreshold is an options for NNDeNoise and specifies the automated back ground masking value."


AnisoStepTime::usage =
"AnisoStepTime is an option for AnisoFilterTensor and defines the diffusion time, when small more step are needed."

AnisoFilterSteps::usage =
"AnisoFilterSteps is an option for AnisoFilterTensor and defines the amoutn of diffusin steps taken. Higher is more smoothing."

AnisoWeightType::usage =
"AnisoWeightType is an option for AnisoFilterTensor and WeightMapCalc and defines the weighting, eigher 1, the exponent of (-g/kappa) or 2, 1/(1+g/kappa)."

AnisoKappa::usage =
"AnisoKappa is an option for AnisoFilterTensor and WeightMapCalc and defines the weighting strenght, all data is normalize to 100 before filetering."

AnisoIterations::usage = 
"AnisoIterations is an options for AnisoFilterData. It specifies the amount of denoising iterations."

AnisoKernel::usage = 
"AnisoKernel is an options for AnisoFilterData. It defines the kernel size."


RadialBasisKernel::usage =
"RadialBasisKernel is an option for HarmonicDenoiseTensor. It defines the kernel size of the radial basis functions in mm."

GradientStepSize::usage =
"GradientStepSize is an option for HarmonicDenoiseTensor. It defines the step size of the gradient descent for the harmonic and radial parts."

RangeFA::usage =
"RangeFA is an option for HarmonicDenoiseTensor. It defines the range of the FA values of voxels to include in the minization."

RangeMD::usage =
"RangeMD is an option for HarmonicDenoiseTensor. It defines the range of the MD values of voxels to include in the minization."



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


debugDenoise[x___] := If[$debugDenoise, MonitorFunction[x]];


(* ::Subsection:: *)
(*Denoise*)


(* ::Subsubsection::Closed:: *)
(*Denoise*)


Options[DeNoise] = {DeNoiseKernel -> "Gaussian", DeNoiseMonitor -> False, DeNoiseIterations -> 1};

SyntaxInformation[DeNoise] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

DeNoise[dat_, sigi_, filt_, OptionsPattern[]] := Module[{kern, out, type, dimsig, dimdat, data, sig},
	sig = N@sigi;
	data = ToPackedArray@N@dat;
	(*Check dimensions,must be of lower order than data*)
	If[ArrayQ[sig] && (ArrayDepth[data] == 3 || ArrayDepth[data] == 4),
		dimsig = Dimensions[sig];
		dimdat = Dimensions[data];
		If[(ArrayDepth[data] == 3 && dimsig != dimdat) || ArrayDepth[data] == 4 && dimsig != dimdat[[{1, 3, 4}]],
			Return[Message[DeNoise::dim, dimsig, dimdat]]
		];
	];

	(*Check dimensions,filter must be of lower order than data*)
	If[Length[filt] > ArrayDepth[data], Return[Message[DeNoise::filt, Length[filt], ArrayDepth[data]]]];

	(*Create filter*)
	kern = Switch[type = OptionValue[DeNoiseKernel],
		"Box", kern = BoxMatrix[filt]/Total[Flatten[BoxMatrix[filt]]],
		"Disk", kern = DiskMatrix[filt]/Total[Flatten[DiskMatrix[filt]]],
		"Gaussian", kern = GaussianMatrix[{filt}],
		_, Return[Message[DeNoise::kern, type]];
	];
	If[OptionValue[DeNoiseMonitor], PrintTemporary["Using " <> type <> " kernel."]];

	out = ToPackedArray@N@Nest[DeNoisei[#, sig, filt, kern] &, data, OptionValue[DeNoiseIterations]];
	If[
		ArrayQ[out], Return[Clip[out, {0.9 Min[data], 1.1 Max[data]}]],
		Return[Message[DeNoise::data, Dimensions[sig], Dimensions[data]]]
	]
]


(* ::Subsubsection::Closed:: *)
(*DeNoisei*)


(*Appling 2D kernel on 2D dataset*)
DeNoisei[data_?MatrixQ,sig_,filt_,kern_?MatrixQ]:=DeNoiseApp[data,sig,filt,kern]

(*Appling 3D kernel on 3D dataset*)
DeNoisei[data:{_?MatrixQ..},sig_,filt_,kern:{_?MatrixQ..}]:=DeNoiseApp[data,sig,filt,kern]

(*Appling 2D kernel on 3D dataset, Using one sigma value for all slices*)
DeNoisei[data:{_?MatrixQ..},sig_,filt_,kern_?MatrixQ]:=If[NumberQ[sig],
	Map[DeNoiseApp[#,sig,filt,kern]&,data],
	MapThread[DeNoiseApp[#1,#2,filt,kern]&,{data,sig}]
]

(*Appling 2D kernel on 4D dataset, Using one sigma value for all slices and directions.*)
DeNoisei[data:{{_?MatrixQ..}..},sig_,filt_,kern_?MatrixQ]:=If[NumberQ[sig],
	Map[DeNoiseApp[#,sig,filt,kern]&,data,{2}],
	Transpose[Map[MapThread[DeNoiseApp[#1,#2,filt,kern]&,{#,sig}]&,Transpose[data]]]
]

(*Appling 3D kernel on 4D dataset, Using one sigma value for all diffusion directions*)
DeNoisei[data:{{_?MatrixQ..}..},sig_,filt_,kern:{_?MatrixQ..}]:=Transpose[Map[DeNoiseApp[#1,sig,filt,kern]&,Transpose[data]]]


(* ::Subsubsection::Closed:: *)
(*DeNoiseApp*)


DeNoiseApp[data_, sig_, filt_, kern_] := Module[{secmod, quadmod},
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


(* ::Subsubsection:: *)
(*PCADeNoise*)


Options[PCADeNoise] = {
	PCAKernel -> 5, 
	PCAOutput -> False,
	PCATolerance -> 0, 
	PCAWeighting -> False, 
	PCAClipping -> True,
	PCAComplex->False,
	Method -> "Similarity",
	MonitorCalc -> False,
	Parallelize -> True
};

SyntaxInformation[PCADeNoise] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

PCADeNoise[data_, opts : OptionsPattern[]] := PCADeNoise[data, 1, 0., opts];

PCADeNoise[data_, mask_, opts : OptionsPattern[]] := PCADeNoise[data, mask, 0., opts];

PCADeNoise[datai_, maski_, sigmai_, OptionsPattern[]] := Block[{
		weight, ker, tol, mon, data, min, max, maskd, mask, sigm, dim, zdim, ydim, xdim, ddim, 
		k, m, n, off, datao, weights, sigmat, start, sigmati, nmati, clip, comp, len,
		output, j, sigi, zm, ym, xm, zp, yp, xp, fitdata, sigo, Nes, datn, 
		weightFun, posV, leng, nearPos, p, pi, pos, np, met, trans, comps,
		mid, pad, fun, par, dimv, coors
	},

	(*tollerane if>0 more noise components are kept*)
	{met, mon, weight, tol, ker, clip, comp, par} = OptionValue[{Method, MonitorCalc, PCAWeighting, 
		PCATolerance, PCAKernel, PCAClipping, PCAComplex, Parallelize}];
	If[par, weight = False];

	(*concatinate complex data and make everything numerical to speed up*)
	comps = False;
	data = ToPackedArray@N@If[comp, 
		clip = False;
		data = If[Length@datai===2, datai, comps = True; Through[{Re, Im}[datai]]];
		Join[data[[1]], data[[2]], 2],
		datai
	];

	(*prep data and masks*)
	{min, max} = 1.1 MinMax[Abs[data]];
	maskd = Unitize@Total@Transpose[data];
	mask = ToPackedArray[N@(maski maskd)];
	sigm = ToPackedArray[N@sigmai];

	(*get data dimensions*)
	dim = {zdim, ddim, ydim, xdim} = Dimensions[data];
	ker = If[EvenQ[ker] && met =!= "Similarity", ker - 1, ker];
	k = Round[ker^3];
	{trans, m, n} = If[ddim < k, {False, ddim, k}, {True, k, ddim}];
	mid = Ceiling[ker^3/2];

	(*define sigma*)
	sigm = If[NumberQ[sigm], ConstantArray[sigm, {zdim, ydim, xdim}], sigm];

	(*define weighting function*)
	weightFun = If[weight, 1./(m - #)&, 1.&];
	
	(*padd the data with kernel offset*)
	off = Round[(ker - 1)/2];
	pad = (off + 1){{1, 1}, {1, 1}, {1, 1}};
	data = ArrayPad[data, Insert[pad, {0, 0}, 2], 0.];
	mask = ArrayPad[mask, pad, 0.];
	sigm = ArrayPad[sigm, pad, 0.];
	
	Switch[met,

		(*--------------use similar signals--------------*)
		"Similarity",
		(*vectorize data*)
		{data, posV} = DataToVector[data, mask];
		sigm = First@DataToVector[sigm, mask];

		(*parameters for monitor*)
		leng = Length[data];
		pos = Range@leng;

		(*get positions of similar signals but in random batches such that nearest has speed*)
		SeedRandom[1234];
		np = 50000;
		pos = If[leng <= np, {pos}, 
			np = Ceiling[leng/Ceiling[leng/np]];
			Partition[RandomSample[pos], np, np, {1, 1}, {}]
		];
		
		(*for each random batch find the nearest voxels and make pos array*)
		If[mon, PrintTemporary["Preparing data similarity"]];
		nearPos = Flatten[Nearest[data[[#]] -> #, data[[#]], 2 k, 
			DistanceFunction -> EuclideanDistance, Method -> "Scan", 
			WorkingPrecision -> MachinePrecision] & /@ pos, 1];
		nearPos[[All, 1]] = Flatten[pos];
		nearPos = Prepend[RandomSample[Rest@#, k - 1], First@#] & /@ Sort[nearPos];

		(*perform denoising*)
		j = 0;
		If[mon, PrintTemporary["Performing the denoising"]];
		If[mon&&!par, PrintTemporary[ProgressIndicator[Dynamic[j], {0, leng}]]];
		
		If[weight,
			(*--------------use weighting---------------*)
			
			(*ouput data*)
			datao = 0. data;
			weights = sigmat = sigmati = nmati = datao[[All, 1]];
		
			Map[(
				If[!par,j++];
				p = #;
				pi = First@p;
	
				(*perform the fit and reconstruct the noise free data*)
				{sigo, Nes, datn} = PCADeNoiseFit[data[[p]], sigm[[pi]], {trans, m, n},  tol];
				weight = weightFun[Nes];
	
				(*sum data and sigma and weight for numer of components*)
				datao[[p]] += weight datn;
				sigmat[[p]] += weight sigo;
				weights[[p]] += weight;
	
				(*output sig, Nest and iterations*)
				sigmati[[pi]] = sigo;
				nmati[[pi]] = Nes;
			) &, nearPos];
			
			(*make everything in arrays*)
			datao = VectorToData[DivideNoZero[datao, weights], posV];
			sigmat = VectorToData[DivideNoZero[sigmat, weights], posV];
			output = {VectorToData[sigmati, posV], VectorToData[nmati, posV]};
			
			,
			
			(*--------------no weighting allows parallel ---------------*)

			fun = If[par, DistributeDefinitions[data, sigm, trans, m, n, tol, par, PCADeNoiseFit]; ParallelMap, Map];
			{datao, sigmat, sigmati, nmati} = Transpose@fun[(If[!par, j++];
				p = #;
				pi = First@p;
				{sigo, Nes, datn} = PCADeNoiseFit[data[[p]], sigm[[pi]], {trans, m, n}, tol];
				{First@datn, sigo, sigo, Nes}
			) &, nearPos];

			(*make everything in arrays*)
			datao = VectorToData[datao, posV];
			sigmat = VectorToData[sigmat, posV];
			output = {VectorToData[sigmati, posV], VectorToData[nmati, posV]};
		];

		,
		(*--------------use patch---------------*)

		_,
		(*vectorize problem*)
		{dimv, coors} = DataToVector[data, mask][[2]];
		
		(*define runtime parameters and data*)
		pad = off + 1;
		data[[;;pad]] = Reverse[data[[-2 pad;;-pad-1]]];
		data[[-pad;;]] = Reverse[data[[pad+1;;2 pad]]];
		data = RotateDimensionsLeft[Transpose[data]];
		
		(*parameters for monitor*)
		j = 0;
		leng = Length@coors;
		If[mon&&!par, PrintTemporary[ProgressIndicator[Dynamic[j], {0, leng}]]];
		
		If[weight,
		
			(*ouput data*)
			datao = 0. data;
			weights = sigmat = datao[[All, All, All, 1]];
	
			(*perform denoising*)
			output = Transpose@Map[(j++;
					{z, y, x} = #;
					
					(*define initial sigma and get pixel range and data*)
					{{zm, ym, xm}, {zp, yp, xp}} = {{z, y, x} - off, {z, y, x} + off};
					fitdata = ArrayReshape[data[[zm ;; zp, ym ;; yp, xm ;; xp]], {k, ddim}];
					sigi = sigm[[z, y, x]];
	
					(*perform the fit and reconstruct the noise free data*)
					{sigo, Nes, datn} = PCADeNoiseFit[fitdata, sigi, {trans, m, n}, tol];
	
					(*reshape the vector into kernel box and get the weightes*)
					datn = Fold[Partition, datn, {ker, ker}];
					weight = weightFun[Nes];
	
					(*sum data and sigma and weight for numer of components*)
					datao[[zm ;; zp, ym ;; yp, xm ;; xp, All]] += (weight datn);
					sigmat[[zm ;; zp, ym ;; yp, xm ;; xp]] += weight sigo;
					weights[[zm ;; zp, ym ;; yp, xm ;; xp]] += weight;
	
					(*output sig, Nest and iterations*)
					{sigo, Nes}
				)&, coors];
	
			(*correct output data for weightings*)
			datao = Transpose@RotateDimensionsRight[Re@DivideNoZero[datao, weights]];
			sigmat = DivideNoZero[sigmat, weights];
			output = VectorToData[#, {dimv, coors}] & /@ output;
			
			,
			
			fun = If[par, DistributeDefinitions[data, sigm, trans, m, n, tol, par, k, ddim, mid, off, pad, PCADeNoiseFit];ParallelMap, Map];
			{datao, sigmat, sigmati, nmati} = Transpose@fun[(If[!par, j++];
				{z, y, x} = #;
				{{zm, ym, xm}, {zp, yp, xp}} = {{z, y, x} - off, {z, y, x} + off};
				fitdata = ArrayReshape[data[[zm ;; zp, ym ;; yp, xm ;; xp]], {k, ddim}];
				{sigo, Nes, datn} = PCADeNoiseFit[fitdata, sigm[[z, y, x]], {trans, m, n}, tol];
				{datn[[mid]], sigo, sigo, Nes}
			)&, coors];
			
			datao = VectorToData[datao, {dimv, coors}];
			sigmat = VectorToData[sigmat, {dimv, coors}];
			output = {VectorToData[sigmati, {dimv, coors}], VectorToData[nmati, {dimv, coors}]};
		];
	];
	
	pad = off + 2;
	datao = datao[[pad;;-pad, All, pad;;-pad, pad;;-pad]];
	sigmat = sigmat[[pad;;-pad, pad;;-pad, pad;;-pad]];
	output = output[[All, pad;;-pad, pad;;-pad, pad;;-pad]];

	(*define output, split it if data is complex*)
	datao = Which[
		comp, len = Round[Length[datao[[1]]]/2]; {datao[[All,1;;len]], datao[[All,len+1;;-1]]},
		clip, Clip[datao, {min, max}], 
		True, datao
	];

	If[comps, datao = datao[[1]]+ I datao[[2]]];

	If[OptionValue[PCAOutput],
		(*fitted dta,average sigma,{sigma fit,number components, number of fitted voxesl,number of max fits}*)
		{datao, sigmat, output},
		{datao, sigmat}
	]
]


(* ::Subsubsection:: *)
(*PCADeNoiseFit*)


(*internal function*)
PCADeNoiseFit[data_, sigi_?NumberQ, {trans_, m_, n_}, toli_] := Block[{
		xmat, xmatT, val, mat, pi, sig, xmatN, tol, out
	},
	
	(*perform decomp*)
	{xmat, xmatT} = If[trans, {data, Transpose@data}, {Transpose@data, data}];
	{val, mat} = Reverse /@ Eigensystem[xmat . xmatT];

	(*if sigma is given perform with fixed sigma,else fit both*)
	{pi, sig} = GridSearch[Re[val], m, n, sigi];

	(*constartin pi plus tol*)
	tol = Round[Clip[pi - toli, {1, m}]];

	(*give output,simga,number of noise comp,and denoised matrix*)
	out = Transpose[mat[[tol ;;]]] . (mat[[tol ;;]] . xmat);
	out = If[trans, out, Transpose@out];
	
	{sig, tol, out}
]


GridSearch = Compile[{{val, _Real, 1}, {m, _Integer, 0}, {n, _Integer, 0}, {sig, _Real, 0}}, Block[
	{valn, sigq1, sigq2, p, pi, si},

	(*calculate all possible values for eq1 and eq2*)
	valn = val[[;; -2]] / n;
	p = Range[m - 1];
	sigq1 = Accumulate[valn] / p;
	sigq2 = (valn - First[valn]) / (4 Sqrt[N[p/(n - (m - (p + 1)))]]);

	(*find at which value eq1>eq2*)
	pi = Max[1, Total[1 - UnitStep[If[sig === 0.,
		sigq2 - sigq1,
		Mean[{sigq1, sigq2}] - sig^2
	]]]];
	si = Sqrt[Max[0, (sigq1[[pi]] + sigq2[[pi]])/2]];

	(*give output*)
	{pi, si}
], RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*NNDeNoise*)


Options[NNDeNoise] = {NNThreshold -> 2};

SyntaxInformation[NNDeNoise] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

NNDeNoise[data_, opts : OptionsPattern[]] := NNDeNoise[data, 1, opts];

NNDeNoise[data_, mask_, opts : OptionsPattern[]] := Block[{
		back, dat, coor, n, ran, dati, train, i
	},
	(*make selection mask and vectorize data*)
	back = Round[mask Mask[NormalizeMeanData[data], OptionValue[NNThreshold]]];
	{dat, coor} = DataToVector[data, back];
	dat = ToPackedArray@N@dat;

	(*Get dimensions, training sample and define network*)
	n = Range@Length@First@dat;
	ran = 1.1 MinMax[dat];

	(*trian network per volume and generate denoised data*)
	(*DistributeDefinitions[dat, n];*)
	dat = Monitor[Table[
		train = Thread[(dati = dat[[All, Complement[n, {i}]]]) -> dat[[All, i]]];
		trained = Predict[train, ValidationSet -> RandomSample[train, Round[.1 Length[dat]]],
			Method -> {"LinearRegression", "OptimizationMethod" -> "StochasticGradientDescent","L2Regularization"->0.01},
			PerformanceGoal -> {"DirectTraining"},
			AnomalyDetector -> None, TrainingProgressReporting -> None, MissingValueSynthesis -> None
		];
		trained[dati]
	, {i, n}], ProgressIndicator[Dynamic[i], {0, Max[n]}]];

	ToPackedArray@N@Clip[Transpose[VectorToData[#, coor] & /@ dat], ran]
]


(* ::Subsubsection::Closed:: *)
(*DenoiseCSIdata*)


Options[DenoiseCSIdata] = {PCAKernel -> 5, PCANoiseSigma->"Corners", PCAWeighting ->True, Method->"Patch"}

SyntaxInformation[DenoiseCSIdata]={"ArgumentsPattern"->{_, OptionsPattern[]}}

DenoiseCSIdata[spectra_, OptionsPattern[]] := Block[{sig, out, hist, len, spectraDen, nn ,sel},
	(* assusmes data is (x,y,z,spectra)*)
	len = Dimensions[spectra][[-1]];

	(*get the corner voxels to calcluate the noise standard deviation or automatic estimation*)

	sig = Switch[OptionValue[PCANoiseSigma],
		"Corners2",
		0.9 StandardDeviation[Flatten[spectra[[{1, -1}, {1, -1}, {1, -1}]]]],
		"Corners", 
		nn = Flatten[spectra[[{1, -1}, {1, -1}, {1, -1}]]];
		sel = FindOutliers[Re@nn, OutlierRange -> 5] FindOutliers[Im@nn, OutlierRange -> 5];
		nn = Pick[nn, sel, 1];
		StandardDeviation[nn]/Sqrt[2]
		,
		"Automatic", 0
	];

	(*Denoise the spectra data*)
	{spectraDen, sig} = PCADeNoise[Transpose[Join[Re@#, Im@#]]&[RotateDimensionsRight[spectra]], 1, sig, 
		PCAClipping -> False, MonitorCalc->False, Parallelize -> True, 
		PCAKernel -> OptionValue[PCAKernel], PCAWeighting ->OptionValue[PCAWeighting], Method->OptionValue[Method]
	];

	ToPackedArray@N@RotateDimensionsLeft[Transpose[spectraDen][[1 ;; len]] + Transpose[spectraDen][[len + 1 ;;]] I]
]


(* ::Subsubsection::Closed:: *)
(*DenoiseDynamicSpectraData*)


SyntaxInformation[DenoiseDynamicSpectraData]={"ArgumentsPattern"->{_}}

DenoiseDynamicSpectraData[spectra_] := Block[{len, data, sig, comp, trans, m, n},
	(*merge Re and Im data*)
	len = Dimensions[spectra][[-1]];
	data = Join[Re@#, Im@#] &[Transpose[spectra]];

	{m, n} = Dimensions[data];
	{trans, m, n} = If[m < n, {True, m, n}, {False, n, m}];

	(*perform denoising*)	
	{sig, comp, data} = PCADeNoiseFit[data, 0., {trans, m, n}, 0];

	(*reconstruct complex spectra*)
	data = Transpose[data[[;;len]] + I data[[len+1;;]]];

	(*output data and sigma*)
	{ToPackedArray@N@data, sig}
]


(* ::Subsection:: *)
(*AnisotropicFilterTensor*)


(* ::Subsubsection::Closed:: *)
(*AnisotropicFilterTensor*)


Options[AnisoFilterTensor] = {AnisoWeightType->2, AnisoKappa->5., AnisoStepTime->1, AnisoFilterSteps->5};

SyntaxInformation[AnisoFilterTensor] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

AnisoFilterTensor[tensi_,opts:OptionsPattern[]]:=AnisoFilterTensor[tensi, 1, opts]

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
	weights = If[dat===1,
		1,
		PrintTemporary["Determaning the weights based on the data."];
		WeightMapCalc[ToPackedArray@N@dat, AnisoKappa->kappa, AnisoWeightType->type]
	];

	(*get the fixed parameters*)
	tens = ToPackedArray@N@tensi;
	mn = Mean[tens[[1;;3]]];
	{kers, wts} = KernelWeights[];
	lambda = 1/Length[kers];

	(*filter the tensor*)
	PrintTemporary["Anisotropic filtering of the tensor."];
	j=0;PrintTemporary[ProgressIndicator[Dynamic[j],{0,itt 6}]];
	Table[
		(*Normalize the diffusion tensor*)
		datf = 100 DivideNoZero[tens[[tt]],mn];
		(*perform the diffusion smoothing iterations*)
		Do[
			j++;
			finDiff = FinDiffCalc[datf,kers];
			wtsI = weights WeightCalc[finDiff,wts,kappa,type];
			datf = datf+time lambda Total@(wtsI finDiff);
			,itt
		];
		(*revert tensor normalization*)
		datf=mn datf/100
		(*loop over tensor*)
	,{tt,1,6}]
]


(* ::Subsubsection::Closed:: *)
(*WeightMapCalc*)


Options[WeightMapCalc]={AnisoWeightType->2, AnisoKappa->10.};

SyntaxInformation[WeightMapCalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

WeightMapCalc[data_,OptionsPattern[]]:=Block[{
	kers,wts,weights,finDiff,dat,dim,len, kappa, type
	},
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
		dat=100 #/Max[Abs[#]];
		(*add to the weights*)
		weights+=WeightCalc[FinDiffCalc[dat,kers],wts,kappa,type];
	)&/@Transpose[ToPackedArray@N@data];

	(*normalize the weights between 0 and 1*)
	(*weights=Mean[weights];
	weights=weights/Max[weights];*)
	(#/Max[#])&/@weights
]


(* ::Subsubsection::Closed:: *)
(*KernelWeights*)


KernelWeights[]:=Block[{cent, ker, keri, wtsi},
	ker = ConstantArray[0, {3, 3, 3}];
	ker[[2, 2, 2]] = -1;
	cent = {2, 2, 2};
	Transpose[Flatten[Table[
		If[{i, j, k}==cent,
			Nothing,
			keri = ker; keri[[i, j, k]] = 1;
			wtsi = N@Norm[cent - {i, j, k}];
			{keri, 1/wtsi^2}
		]
	, {i, 1, 3}, {j, 1, 3}, {k, 1, 3}], 2]]
];


(* ::Subsubsection::Closed:: *)
(*WeightCalc*)


WeightCalc[finDiff_,wts_,kappa_,type_] := wts Switch[type,1,Exp[-((finDiff/kappa)^2)],2,1./(1.+(finDiff/kappa)^2)];


(* ::Subsubsection::Closed:: *)
(*FinDiffCalc*)


FinDiffCalc[dat_,kers_] := ParallelMap[ListConvolve[#,dat,{2,2,2},0]&,kers]


(* ::Subsection::Closed:: *)
(*AnisoFilterData*)


Options[AnisoFilterData] = {
	AnisoStepTime -> 1, 
	AnisoIterations -> 1, 
	AnisoKernel -> {0.25, 0.5}};

SyntaxInformation[AnisoFilterData] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

AnisoFilterData[data_, opts:OptionsPattern[]] := AnisoFilterData[data, {1, 1, 1}, opts]

AnisoFilterData[data_, vox_, opts:OptionsPattern[]] := Block[{
	dd, grads, k, jacTot, tMat, eval, evec, div, dati,
	sig, rho , step, itt, sc, max, tr, cr, ind
	},
	(*for implementation see 10.1002/mrm.20339*)

	(*crop background for speed*)
	{dati, cr} = AutoCropData[data];

	(*get the 4th dimensions on first place *)
	dati = ToPackedArray@N@Transpose[dati];
	max = Max[dati];
	tr = max/10000;
	ind = IdentityMatrix[3];
	sc = Max[vox]/vox;

	(*aniso filter kernel size *)
	{sig, rho} = OptionValue[AnisoKernel];
	step = OptionValue[AnisoStepTime];
	itt = OptionValue[AnisoIterations];

	Do[(*loop over itterrations*)
		(*get the data gradients*)
		grads = ToPackedArray@N[(
			dd = GaussianFilter[#, {1, sc sig}];
			sc (GaussianFilter[dd, 1, #1] & /@ ind)
		) & /@ dati];

		(*calculate the jacobian (aka structure tensor) and smooth it*)
		jacTot = ToPackedArray@N@TensMat[
			GaussianFilter[#, {1, sc rho}] & /@ Total[
				{#[[1]]^2, #[[2]]^2, #[[3]]^2, #[[1]] #[[2]], #[[1]] #[[3]], #[[2]] #[[3]]} &[Transpose[grads]]
			, {2}]
		];

		(*get the step matrix*)
		tMat = ToPackedArray@N@Map[If[Max[#] < tr, 
			{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
			{eval, evec} = Eigensystem[#];
			eval = eval = 1/(eval + 10.^-16);
			Transpose[evec] . DiagonalMatrix[3 eval/Total[eval]] . evec
		] &, jacTot, {3}];

		(*get the time step and calculate divergence of vector field*)
		div = Total[MapThread[GaussianFilter[#1, 1, #2] &, {#, ind}] & /@ 
			RotateDimensionsRight[DivDot[tMat, RotateDimensionsLeft[grads, 2]], 2], {2}];

		(*perform the smoothing step*)
		dati = ToPackedArray@N@Clip[dati + step div, {0, 2} max];
	, {itt}];

	(*output the data*)
	ReverseCrop[Transpose@dati, Dimensions@data, cr]
]


DivDot = Compile[{{t, _Real, 2}, {gr, _Real, 2}}, 
	gr . t, 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


StrucTensCalc = Compile[{{eval, _Real, 1}, {evec, _Real, 2}, {tr, _Real, 0}},
	If[Max[eval] < tr, {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
		Transpose[evec] . DiagonalMatrix[3 eval/Total[eval]] . evec
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection:: *)
(*HarmonicDenoiseTensor*)


(* ::Subsubsection:: *)
(*HarmonicDenoiseTensor*)


Options[HarmonicDenoiseTensor] = {
	RadialBasisKernel -> 12,
	TensorFlips -> {1, 1, 1},
	TensorPermutations -> {"x", "y", "z"},
	MaxIterations -> 150,
	GradientStepSize -> {0.5,0.5},
	Tolerance -> 10.^-5,
	RangeFA -> {0.05, 0.4},
	RangeMD -> {1., 2.5},
	Monitor -> False
};


SyntaxInformation[HarmonicDenoiseTensor] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

HarmonicDenoiseTensor[tensI__?ArrayQ, seg_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, opts:OptionsPattern[]]:=
	HarmonicDenoiseTensor[tensI, seg, vox, 0, opts]

HarmonicDenoiseTensor[tensI__?ArrayQ, segI_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, labs_, OptionsPattern[]]:=Block[{
		sigma, flip, per, itt, step, tol, rFA, rMD, seg, lab, mon, tensL, crops, denoise,
		tensO, dimT, conO, dimC, ampO, dimA, mus, crp, tens, con, amp
	},

	(*get options*)
	{sigma, flip, per, itt, step, tol, rFA, rMD, mon}=OptionValue[{RadialBasisKernel, TensorFlips, TensorPermutations,
		MaxIterations, GradientStepSize, Tolerance, RangeFA, RangeMD, Monitor}];

	(*figure out how to loop over the masks, seg is 4D with first index the segmentations*)
	seg = Which[
		ArrayDepth[segI]==3 && labs===0,
		{segI},
		ArrayDepth[segI]==3 && ListQ[labs],
		Transpose@First@SelectSegmentations[segI, labs, False],
		ArrayDepth[segI]==4 && labs===0,
		Transpose[segI],
		ArrayDepth[segI]==4 && ListQ[labs],
		Transpose[segI][[labs]]
	];

	(*make cropped data for each muscle*)
	If[mon, MonitorFunction["Making cropped muscles"]];
	tensL = Transpose[FlipTensorOrientation[tensI, per, flip]];
	crops = Table[crp = FindCrop[vol, CropPadding -> sigma];
		{Normal@Transpose[ApplyCrop[tensL, crp]], Normal@ApplyCrop[vol, crp], crp}
	, {vol, Transpose[seg[[All, pos]]]}];

	(*paralell map over all muscles and denoise*)
	If[mon, MonitorFunction["Starting parallel denoising"]];
	DistributeDefinitions[HarmonicDenoiseTensorI, sigma, itt, step, tol, rFA, rMD, vox];
	denoise = ParallelMap[HarmonicDenoiseTensorI[#[[1]], #[[2]], vox,
		RadialBasisKernel->sigma, MaxIterations->itt, GradientStepSize->step, 
		Tolerance->tol,	RangeFA->rFA, RangeMD->rMD
	]&, crops, ProgressReporting -> mon];

	(*prepare the output*)
	If[mon, MonitorFunction["Reasembling data"]];
	tensO = Transpose[SparseArray[0. tensI]];
	dimT = Dimensions@tensO;
	conO = Transpose[SparseArray[0. tensI[[1;;3]]]];
	dimC = Dimensions@conO;
	ampO = SparseArray[0. First@tensI];
	dimA = Dimensions@ampO;

	(*add result to output*)
	Map[(
		{{tens, con, amp}, crp} = #;
		conO += ReverseCrop[con, dimC, crp];
		ampO += ReverseCrop[amp, dimA, crp];
		tensO += ReverseCrop[Transpose[tens], dimT, crp];
    ) &, Thread[{denoise, crops[[All, 3]]}]];
	tensO = Transpose[tensO];

	(*give the output*)
	{tensO, conO, ampO}
]


Options[HarmonicDenoiseTensorI]=Options[HarmonicDenoiseTensor]

HarmonicDenoiseTensorI[tens_, mask_, vox_, OptionsPattern[]]:=Block[{
		mon, sigma, md, fa, msel, rFA, rMD, G, L, R,coor, sel, map, vecN, 
		vecH, sol, tensN, maps, itt, step, tol, tgl, tr, tt
	},

	{sigma, itt, step, tol, rFA, rMD} = OptionValue[{RadialBasisKernel, MaxIterations,
		GradientStepSize, Tolerance, RangeFA, RangeMD}];
	debugDenoise[{sigma, itt, step, tol, rFA, rMD}, "Settings"];

	(*make gradien,laplace and RBF matrix functions*)
	tt = First@AbsoluteTiming[
		tgl = First@AbsoluteTiming[
			{{G, L}, coor, sel, map} = MakeGradientLaplacian[mask, vox, sigma];
		];
		debugDenoise[tgl, "Making G, L matrix"];
		
		tr = First@AbsoluteTiming[
			R = MakeRBF[coor, sigma, vox];
		];
		debugDenoise[tr, "Making R matrix"];

		(*becouse of the implementation coordiante system the tensor needs to be reversed*)
		tensN = FlipTensorOrientation[MaskData[tens, mask], {"z", "y", "x"}];
		(*remove unreliable voxel*)
		{md, fa} = ParameterCalc[tensN][[4;;5]];
		msel = Mask[{{fa, rFA}, {md, rMD}}];

		(*perform the recon*)
		debugDenoise["starting fitting"];
		{vecN, vecH, sol} = FitHarmonicBasis[MaskData[tensN, msel], sel, {G, L, R},
			MaxIterations->itt, GradientStepSize->step, Tolerance->tol];
	];

	debugDenoise[tt, "Total denoise time"];

	(*generate output*)
	maps = MakeSolutionMaps[sol, map];
	tensN = FlipTensorOrientation[ReconstrucTensor[vecN, tensN, sel, coor], {"z", "y", "x"}];
	{tensN, Transpose[maps[[1;;3]]], maps[[4]]}
]


(* ::Subsubsection::Closed:: *)
(*MakeGradientLaplacian*)


MakeGradientLaplacian[maski_]:=MakeGradientLaplacian[maski,{1,1,1},0]

MakeGradientLaplacian[maski_,vox_?VectorQ]:=MakeGradientLaplacian[maski,vox,0]

MakeGradientLaplacian[maski_,sig_?NumberQ]:=MakeGradientLaplacian[maski,{1,1,1},sig]

MakeGradientLaplacian[mask_,vox_?VectorQ,sig_?NumberQ]:=Block[{aDepth,dim,const,matI,maskDilated,maskPadded,di,
		z, y, x, coors, cx, cy, cz, coorsPad, selMask, selMaskC, selCent, selGrad, selLap, coorMask, coorCent, coorGrad, coorLap,
		ix0, ix1, iy0, iy1, iz0, iz1, gx, gy, gz, lx1, lx0, ly1, ly0, lz1, lz0, xx, xy, yx, yy, zx, zy, g, l, lx, ly, lz, ix, iy, iz, G, L
	},

	(*get mask properties*)
	aDepth = ArrayDepth@mask;
	dim = Dimensions@mask;
	const = ConstantArray[0, aDepth];
	matI = IdentityMatrix[aDepth];

	(*extend matrix in all directions and dilate mask in each direction if needed*)
	maskDilated = If[sig === 0, mask, Dilation[mask, BoxMatrix[Round[sig/vox]]]];
	maskPadded = (di = ArrayPad[ArrayPad[maskDilated, Thread[{const, #}]], 1];
	Unitize[di + RotateRight[di,#]])& /@ Reverse[matI];

	(*define which points can be selected for which application*)
	selMask = Flatten[mask];(*where the data is*)
	selCent = Flatten[maskDilated];(*where the padded centers are*)
	selMaskC = Pick[selMask, selCent,1];
	selGrad = Flatten[ArrayPad[#, -1]]& /@ maskPadded; (*where the concentrations for the gradients are*)
	selLap = Table[Flatten[ArrayPad[Times@@(RotateRight[di, #]&/@Join[matI, -matI]),-1]], {di, maskPadded}];

	(*make the coordinates for the centers, the x shifted and y shifted cells*)
	Switch[aDepth,
		2,
		{y, x} = dim;
		coors = Flatten[Table[{j, i }, {j, 1, y}, {i, 1, x}], 1];
		cx = Flatten[Table[{j, i }, {j, 1, y}, {i, 1/2, x+1/2}], 1];
		cy = Flatten[Table[{j, i }, {j, 1/2, y+1/2}, {i, 1, x}], 1];
		coorsPad = {cx, cy};,
		3,
		{z, y, x} = dim;
		coors = Flatten[Table[{k, j, i}, {k, 1, z}, {j, 1, y}, {i, 1, x}], 2];
		cx = Flatten[Table[{k, j, i}, {k, 1, z},{j, 1, y},{i, 1/2, x+1/2}], 2];
		cy = Flatten[Table[{k, j, i}, {k, 1, z},{j, 1/2, y+1/2},{i, 1, x}], 2];
		cz = Flatten[Table[{k, j, i}, {k, 1/2, z+1/2},{j, 1, y},{i, 1, x}], 2];
		coorsPad = {cx, cy, cz};
	];

	(*selecet the correct coordinates for visualization*)
	coorMask = Pick[coors, selMask];
	coorCent = Pick[coors, selCent, 1];
	coorGrad = Pick[#[[1]], #[[2]], 1]&/@Transpose[{coorsPad, selGrad}];
	coorLap = Pick[#[[1]], #[[2]], 1]&/@Transpose[{coorsPad, selLap}];

	(*make the Gradient and Laplacien matrix for 2D and 3D*)
	Switch[aDepth

		(*----2D case----*)
		, 2,
		(*band and identiy matrixes for gradient and laplacian matrixes*)
		{ix0, iy0} = SparseArray[IdentityMatrix[#]]&/@{x, y};
		{ix1, iy1} = SparseArray[IdentityMatrix[#+1]]&/@{x, y};
		{gx, gy} = SparseArray[{Band[{1, 1}]->-1., Band[{1, 2}]->1.}, {#, #+1}]&/@{x, y};
		{lx1, ly1} = SparseArray[{Band[{1, 1}]->-2., Band[{1, 2}]->1., Band[{2, 1}]->1.}, {#+1, #+1}]&/@{x, y};
		{lx0, ly0} = SparseArray[{Band[{1, 1}]->-2., Band[{1, 2}]->1., Band[{2, 1}]->1.}, {#, #}]&/@{x, y};

		(*make the matrix for calculating gradient of concentrations for x, y, and z shifted cells*)
		g = Table[
			g = KroneckerProduct[If[i == 2, gy, iy0], If[i == 1, gx, ix0]];
			Pick[Transpose[Pick[Transpose[g], selGrad[[i]], 1]], selMask, 1]
		, {i, 1, 2}];

		(*combine into one matrix, rows are xyz... for easyer vector calculation later*)
		{{xx, xy}, {yx, yy}} = Dimensions/@g;
		G = Flatten[Transpose[SparseArray[ArrayPad[#[[1]], {{0, 0}, #[[2]]}]&/@Transpose[{g, {{0, yy}, {xy, 0}}}]]], 1];

		(*make the laplacian matrixes for all directions of the shifted cells*)
		{lx, ly} = Table[
			{ly, iy, lx, ix} = {{ly0, iy0, lx1, ix1}, {ly1, iy1, lx0, ix0}}[[j]];
			l = Total[Table[KroneckerProduct[If[i == 2, ly, iy], If[i == 1, lx, ix]], {i, 1, 2}]];
			Pick[Transpose[Pick[Transpose[l], selGrad[[j]], 1]], selLap[[j]], 1]
		, {j, 1, 2}];

		(*join all the Laplacians *)
		{{xx, xy}, {yx, yy}} = Dimensions/@{lx, ly};
		L = ArrayPad[lx, {{0, yx}, {0, yy}}]+ArrayPad[ly, {{xx, 0}, {xy, 0}}];

		(*----3D case----*)
		, 3,
		(*band and identiy matrixes for gradient and laplacian matrixes*)
		{ix0, iy0, iz0} = SparseArray[IdentityMatrix[#]]&/@{x, y, z};
		{ix1, iy1, iz1} = SparseArray[IdentityMatrix[#+1]]&/@{x, y, z};
		{gx, gy, gz} = SparseArray[{Band[{1, 1}]->-1., Band[{1, 2}]->1.}, {#, #+1}]&/@{x, y, z};
		{lx1, ly1, lz1} = SparseArray[{Band[{1, 1}]->-2., Band[{1, 2}]->1., Band[{2, 1}]->1.}, {#+1, #+1}]&/@{x, y, z};
		{lx0, ly0, lz0} = SparseArray[{Band[{1, 1}]->-2., Band[{1, 2}]->1., Band[{2, 1}]->1.}, {#, #}]&/@{x, y, z};

		(*make the matrix for calculating gradient of concentrations for x, y, and z shifted cells*)
		g = Table[
			g = KroneckerProduct[If[i == 3, gz, iz0], KroneckerProduct[If[i == 2, gy, iy0], If[i == 1, gx, ix0]]];
			Pick[Transpose[Pick[Transpose[g], selGrad[[i]], 1]], selMask, 1]
		, 	{i, 1, 3}];

		(*combine into one matrix, rows are xyz... for easyer vector calculation later*)
		{{xx, xy}, {yx, yy}, {zx, zy}} = Dimensions/@g;
		G = Flatten[Transpose[SparseArray[ArrayPad[#[[1]], {{0, 0}, #[[2]]}]&/@Transpose[{g, {{0, yy+zy}, {xy, zy}, {xy+yy, 0}}}]]], 1];

		(*make the laplacian matrixes for all directions of the shifted cells*)
		{lx, ly, lz} = Table[
			{lz, iz, ly, iy, lx, ix} = {{lz0, iz0, ly0, iy0, lx1, ix1}, {lz0, iz0, ly1, iy1, lx0, ix0}, {lz1, iz1, ly0, iy0, lx0, ix0}}[[j]];
			l = Total[Table[KroneckerProduct[If[i == 3, lz, iz], KroneckerProduct[If[i == 2, ly, iy], If[i == 1, lx, ix]]], {i, 1, 3}]];
			Pick[Transpose[Pick[Transpose[l], selGrad[[j]], 1]], selLap[[j]], 1]
		, {j, 1, 3}];

		(*join all the Laplacians *)
		{{xx, xy}, {yx, yy}, {zx, zy}} = Dimensions/@{lx, ly, lz};
		L = ArrayPad[lx, {{0, yx+zx}, {0, yy+zy}}]+ArrayPad[ly, {{xx, zx}, {xy, zy}}]+ArrayPad[lz, {{xx+yx, 0}, {xy+yy, 0}}];
	];

	{{G, L}, {coorCent, selMaskC} ,selMask, {coorCent, coorGrad, dim}}
]


(* ::Subsubsection::Closed:: *)
(*MakeRBF*)


MakeRBF[{coor_, sel_}, rad_]:=MakeRBF[{coor, sel}, rad, 1]

MakeRBF[{coor_, sel_}, rad_, vox_]:=Block[{seed, target, n, m, nr, index, indexFinal, rbf, rows, cols, vals},
	seed = Transpose[Reverse[vox Transpose[ToPackedArray@N@coor]]];
	target = Pick[seed, sel, 1];

	n = Length[First@target];
	m = Length@seed;
	nr = n Total@sel;

	index = Nearest[target -> "Index", seed, {All, 2.5 rad}, 
	DistanceFunction -> EuclideanDistance];
	indexFinal = IndexSwitch[index, n];
	rbf = gaussianRBFGradientC[seed, target, index, 1./rad^2];

	rows = ToPackedArray@Join @@ MapThread[ConstantArray, {Range[Length@indexFinal], Length /@ indexFinal}];
	cols = ToPackedArray@Join @@ indexFinal;
	vals = ToPackedArray@Join @@ rbf;

	SparseArray[Transpose[{rows, cols}] -> vals, {m, nr}, 0.]
]


gaussianRBFGradientC = Compile[{{center, _Real, 1}, {points, _Real, 2}, {ind, _Integer, 1}, {r2, _Real, 0}}, Block[{c},
    c = Transpose[points[[ind]]] - center;
    Flatten[2 r2 Exp[-Total[c^2] r2] Transpose[c], 1]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


IndexSwitch = Compile[{{i, _Integer, 1}, {dm, _Integer, 0}}, 
	Flatten[If[dm === 2, Transpose[{2 i - 1, 2 i}], Transpose[{3 i - 2, 3 i - 1, 3 i}]], 1]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*FitHarmonicBasis*)


Options[FitHarmonicBasis] = {
	MaxIterations -> 250,
	GradientStepSize -> {0.5,0.5},
	Tolerance -> 10.^-4
};

FitHarmonicBasis[tv_, sel_, {matG_ ,matL_, matR_}, opts:OptionsPattern[]]:=Block[{
		itt, step, matGt, h0, vec0, v0, solver, dvh, dvg, i, j, thp, th, tgp, tg, vsh, vsg,
		h1, vec1, v1, vh, a0, a1, vha, d, dim, tens, vec, val, sth, stg, tol, diff, mon
	},

	(*get the options*)
	{itt, step, tol} = OptionValue[{MaxIterations, GradientStepSize, Tolerance}];

	(*vectorize the data*)
	Switch[Length[tv],
		2,
		{tens, vec} = tv;
		dim = Drop[Dimensions@vec, -1];
		d = Length@dim;
		tens = SelectVector[tens, sel, d];
		vec = SelectVector[vec, sel, d];,

		_,
		dim = Dimensions@tv[[1]];
		d = Length@dim;
		tens = SelectVector[#, sel]& /@ tv;
		{val, vec} = EigensysCalc[tens][[All,All,1]];
		vec = Unitize[val] vec;
		tens = TensMat[tens];
	];

	diff = Max@Abs[tens];
	tol = diff tol;

	thp = First@AbsoluteTiming[
	(*initialize the h phase*)
		matGt = Transpose[matG];
		solver = MakeSolver[matL];
		h0 = LinearSolve[matGt . matG, matGt . Flatten[vec]];
		h0 = NullSpaceProjection[matL, h0, solver];
		vec0 = Partition[matG . h0, d];
		v0 = GetDiffusionValues[tens, vec0];
		sth = Norm[GetObjectiveGradient[tens, vec0, matGt]]
	];
	debugDenoise[thp, "Prep time h-phase: "];

	(*h phase loop*)
	dvh = v0; i = 0;
	debugDenoise[Dynamic[{i, dvh/tol, sth}], "Start h-phase: "];

	th = First@AbsoluteTiming[
		vsh = Reap[While[dvh>tol && i<itt,
			i++; Sow[{v0, dvh}];
			h1 = h0 + (100 step[[1]] / sth) GetObjectiveGradient[tens, vec0, matGt];
			h1 = NullSpaceProjection[matL, h1, solver];
			vec1 = Partition[matG . h1, d];
			v1 = GetDiffusionValues[tens, vec1];
			dvh = Abs[v0 - v1];
			h0 = h1; v0 = v1;
			vec0 = vec1;
		]]
	];
	debugDenoise[th, "Fit time h-phase: "];

	(*normalize h0 to max h0 for g phase*)
	vh = vec1 / Median[(Norm /@ vec1)];

	tgp = First@AbsoluteTiming[
		(*initialize the g phase*)
		a0 = ConstantArray[0.,Length@matR];
		vec0 = vh + Partition[a0 . matR, d];
		v0 = GetDiffusionValues[tens, vh];
		stg = Norm[GetObjectiveGradient[tens, vec0, matR]];
	];
	debugDenoise[tgp, "Prep time g-phase: "];

	(*g phase loop*)
	dvg = v0; j = 0;
	debugDenoise[Dynamic[{j,dvg/tol,stg}], "Start g-phase: "];
	tg = First@AbsoluteTiming[
		vsg = Reap[While[dvg>tol && j<Round[itt],
			j++; Sow[{v0, dvg}];
			a1 = a0 + (0.05 step[[2]]/stg) GetObjectiveGradient[tens, vec0, matR];
			vec1 = vh + Partition[a1 . matR, d];
			v1 = GetDiffusionValues[tens, vec1];
			dvg = Abs[v1 - v0];
			a0 = a1; v0 = v1;
			vec0 = vec1;
		]]
	];
	debugDenoise[tg, "Fit time g-phase: "];

	(*normalize vectors after g phase*)
	vha = Normalize/@vec1;

	debugDenoise[thp+th+tgp+tg, "Total fit time: "];

	debugDenoise[Grid[{{
		ListLinePlot[{vsh[[2,1,All,1]], vsg[[2,1,All,1]]} / diff, ImageSize->200, PlotRange->Full,PlotLabel->"Difference"],
		ListLinePlot[{vsh[[2,1,All,2]], vsg[[2,1,All,2]]} / tol, ImageSize->200,PlotLabel->"Tollerance"]
	}}], "Gradient functions: "];

	{vha, Normalize/@vh, {h1,a1}}
]


SelectVector[data_, sel_]:=SelectVector[data, sel, 3]

SelectVector[data_, sel_, d_]:=Pick[Flatten[data, d-1], sel, 1]


(* ::Subsubsection::Closed:: *)
(*GetDiffusionValues*)


GetDiffusionValues[tens_,vec_]:=Mean[GetDiffC[tens,vec]]


GetDiffC=Compile[{{t,_Real,2},{v,_Real,1}},Block[{vv},
	vv = Dot[v, v];
	If[vv == 0., 0., v . t . v/vv]
],RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


(* ::Subsubsection::Closed:: *)
(*GetObjectiveGradient*)


GetObjectiveGradient[tens_, vec_, obj_] := -obj . Flatten[GetGradC[tens, vec]]


GetGradC=Compile[{{t, _Real, 2}, {v, _Real, 1}}, Block[{vv, tv},
	vv = Dot[v, v];
	tv = Dot[t, v];
	If[vv == 0., v, ((v . tv)v - (vv tv))/vv^2]
],RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


(* ::Subsubsection::Closed:: *)
(*NullSpaceProjection*)


NullSpaceProjection[matL_, x_] := x - Transpose[matL] . LinearSolve[matL . Transpose[matL], matL . x, Method->"Pardiso"];
NullSpaceProjection[matL_, x_, sol_] := x - Transpose[matL] . sol[matL . x];


MakeSolver[matL_] := LinearSolve[matL . Transpose[matL], Method->"Pardiso"];


(* ::Subsubsection::Closed:: *)
(*MakeSolutionMaps*)


MakeSolutionMaps[{con_, amp_}, {coorCent_, coorGrad_, dim_}]:=Block[{
		cc, ca, l1, l2, l3, z, y, x, c1, c2, c3, amap
	},

	cc = Quantile[Abs[con], .75];
	ca = Quantile[Abs[amp], .75];

	Switch[Length@coorGrad,
		(*reconstruct for 2D maps*)
		2,
		{l1, l2} = Length/@coorGrad;
		{y, x} = dim;
		c1 = PadRight[SparseArray[Thread[Ceiling[coorGrad[[1]]]->con[[;;l1]]/cc], {y, x+1}], {y, x}];
		c2 = PadRight[SparseArray[Thread[Ceiling[coorGrad[[2]]]->con[[l1+1;;]]/cc], {y+1, x}], {y, x}];
		amap = SparseArray[Thread[coorCent->amp/ca], {y, x}];
		{c1,c2,amap}
		(*reconstruct for 3D maps *)
		,3,
		{l1, l2, l3}=Length/@coorGrad;
		{z, y, x}=dim;
		c1 = PadRight[SparseArray[Thread[Ceiling[coorGrad[[1]]]->con[[;;l1]]/cc], {z, y, x+1}], {z, y, x}];
		c2 = PadRight[SparseArray[Thread[Ceiling[coorGrad[[2]]]->con[[l1+1;;l1+l2]]/cc], {z, y+1, x}], {z, y, x}];
		c3 = PadRight[SparseArray[Thread[Ceiling[coorGrad[[3]]]->con[[l1+l2+1;;l1+l2+l3]]/cc], {z+1, y, x}], {z, y, x}];
		amap = SparseArray[Thread[coorCent->amp/ca], {z, y, x}];

		{c1, c2, c3, amap}
	]
]


(* ::Subsubsection::Closed:: *)
(*ReconstrucTensor*)


ReconstrucTensor[vecN1_, tens_, sel_, coor_]:=Block[{cors, dim, vecO, vecN2, vecN3, vecN, tensN},
	cors = Pick[coor[[1]], coor[[2]],1];
	dim = Dimensions@First@tens;

	vecO = Transpose@EigenvecCalc[RotateDimensionsRight[SelectVector[RotateDimensionsLeft[tens], sel]]];

	vecN2 = MakePerpendicular[vecN1, vecO[[2]]];
	vecN3 = CrossC[vecN1, vecN2];
	tensN = TensVec@MakeTens[Transpose[{vecN1, vecN2, vecN3}],DiagonalMatrix[{2.5, 1.75, 1}/1000]];
	Normal[SparseArray[Thread[cors->#],dim,0.]& /@ tensN]
]


MakePerpendicular = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}},
	Normalize[vec2-(vec2 . vec1) vec1]
,RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"];

CrossC = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}},
	Cross[vec1, vec2]
,RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"];

MakeTens = Compile[{{vec, _Real, 2},{val, _Real, 2}},
	Transpose[vec] . val . vec
, RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
