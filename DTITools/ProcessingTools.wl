(* ::Package:: *)

(* ::Title:: *)
(*DTITools ProcessingTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`ProcessingTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


TensorCalc::usage = 
"TensorCalc[data, gradients, bvalue] calculates the diffusion tensor for the given dataset. Allows for one unweighted image and one b value. 
Gradient directions must be in the form {{x1,y1,z1}, ..., {xn,yn,zn}} without the unweighted gradient direction. 
bvalue is a singe number indicating the b-value used.
TensorCalc[data, gradients, bvec] calculates the diffusion tensor for the given dataset. allows for multiple unweighted images and multiple bvalues. 
allows for differnt tensor fitting methods. gradient directions must be in the form {{x1,y1,z1}, ..., {xn,yn,zn}} with the unweighted direction as {0,0,0}. \
bvec the bvector, with a bvalue defined for each gradient direction. b value for unweighted images is 0.
TensorCalc[data, bmatix] calculates the diffusion tensor for the given dataset. allows for multiple unweighted images and multiple bvalues. 
bmat is the bmatrix which can be generated usiong Bmatrix."

EigenvalCalc::usage = 
"EigenvalCalc[tensor] caculates the eigenvalues for the given tensor."

EigenvecCalc::usage =
"EigenvecCalc[tensor] caculates the eigenvectors for the given tensor."

EigensysCalc::usage = 
"EigensysCalc[tensor] caculates the eigensystem for the given tensor."

ADCCalc::usage =
"ADCCalc[eigenvalues] caculates the ADC from the given eigenvalues."

FACalc::usage =
"FACalc[eigenvalues] caculates the FA from the given eigenvalues."

ECalc::usage =
"ECalc[eigenvalues] caculates the e from the given eigenvalues."

ParameterCalc::usage = "ParameterCalc[tensor] caculates the eigenvalues and MD and FA from the given tensor."

SigmaCalc::usage = 
"SigmaCalc[DTI,grad,bvec] calculates the noise sigma based on the tensor residual, using a blur factor of 10.
SigmaCalc[DTI,tens,grad,bvec] calculates the noise sigma based on the tensor residual, using a blur factor of 10.
SigmaCalc[DTI,grad,bvec,blur] calculates the noise sigma based on the tensor residual, If blur is 1 ther is no blurring.
SigmaCalc[DTI,tens,grad,bvec,blur] calculates the noise sigma based on the tensor residual. If blur is 1 ther is no blurring."

SNRCalc::usage = 
"SNRCalc[data,masksig,masknoise] calculates the Signal to noise ratio of the signal selected by masksig and the noise selected by masknoise."

SNRMapCalc::usage =
"SNRMapCalc[data1,noisemap] calcualtes the signal to noise ratio of the data using MN[data]/(1/sqrt[pi/2] sigma), \
where sigma is the local mean of the noise map assuming it is a rician distribution.

SNRMapCalc[{data1,data2}] calcualtes the signal to noise ratio from two identical images using \
MN[data1,data2] / (.5 SQRT[2] STDV[data2-data1]).

SNRMapCalc[{data1, .. dataN}] calcualtes the signal to noise ratio of the data using MN/sigma where the mean signal MN is the average voxe \
value over all dynamics N and the sigma is the standard deviation over all dynamics N."

CoilSNRCalc::usage = 
"CoilSNRCalc[coils, noise] calculates the sensitivity weighted snr of multiple coil elements using magnitude signal and noise."

PhaseCalc::usage = 
"PhaseCalc[B0data] unwraps the two B0 phase maps and calculates the phase difference between the two sets. Output is in radials."

AngleCalc::usage = 
"AngleCalc[data, vector] calculates the angel between the vector and the data. Data shoud be an array of dimensions {xxx,3}."

AngleMap::usage = 
"AngleMap[data] calculates the zennith and azimuth angles of a 3D dataset (z,x,y,3) containing vectors relative to the slice direction."

ResidualCalc::usage =
"ResidualCalc[DTI,{tensor,S0},gradients,bvector] calculates the tensor residuals for the given dataset.
ResidualCalc[DTI,{tensor,S0},outlier,gradients,bvector] calculates the tensor residuals for the given dataset taking in account the outliers.
ResidualCalc[DTI,{tensor,S0},bmat] calculates the tensor residuals for the given dataset.
ResidualCalc[DTI,{tensor,S0},outlier,bmat] calculates the tensor residuals for the given dataset taking in account the outliers.
ResidualCalc[DTI,tensor,gradients,bvector] calculates the tensor residuals for the given dataset. Tensor must contain Log[S0].
ResidualCalc[DTI,tensor,outlier,gradients,bvector] calculates the tensor residuals for the given dataset taking in account the outliers. Tensor must contain Log[S0].
ResidualCalc[DTI,tensor,bmat] calculates the tensor residuals for the given dataset. Tensor must contain Log[S0].
ResidualCalc[DTI,tensor,outlier,bmat] calculates the tensor residuals for the given dataset taking in account the outliers. Tensor must contain Log[S0]."

FitData::usage = 
"FitData[data,range] converts the data into 100 bins within the +/- range around the mean. Function is used in ParameterFit."

ParameterFit::usage = 
"ParameterFit[data] fits a (skew)Normal probability density function to the data.
ParameterFit[{data1, data2,...}] fits a (skew)Normal probability density function to each of the datasets. Is used in Hist."

ParameterFit2::usage = 
"ParameterFit2[data] fits two skewNormal probaility density fucntions to the data. Assuming two compartments, \
one for fat and one for muscle. Is used in SmartMask2 and Hist2."

DatTot::usage = 
"DatTot[{data1, data2, ..}, name, vox] calculates the parameter table conating the volume, mean, std and 95 CI for each of the diffusion parameters."

DatTotXLS::usage = 
"DatTotXLS[{data1, data2, ..}, name, vox] is the same as DatTot, but gives the parameters as strings for easy export to excel."

SliceData::usage = 
"SliceData[data] calculates the mean and std of the diffuison parameters per slice of data."

MeanSignal::usage = 
"MeanSignal[diffdata] calculates the mean signal per volume of the diff data.
MeanSignal[diffdata, pos] calculates the mean signal per volume of the diff data cor the given positions."

GetMaskMeans::usage = 
"GetMaskMeans[dat, mask, name] calculates the mean, std, 5,50 and 95% CI form the given data for each of the given masks. 
Mask can be genereated by SplitSegmentations. name is a string that is added to the header."

FiberDensityMap::usage =
"FiberDensityMap[fiberPoins, dim, vox] generates a fiber density map for the fiberPoins which are imported by LoadFiberTracts. \
The dimensions dim should be the dimensions of the tracked datasets van vox its volxel size."

FiberLengths::usage =
"FiberLengths[fpoints,flines] calculates the fiber lenght using the output from LoadFiberTacts.
FiberLengths[{fpoints,flines}] calculates the fiber lenght using the output from LoadFiberTacts."

(*FindOutliers::usage = "";*)

(* ::Subsection::Closed:: *)
(*Options*)


SmoothPhase::usage = 
"SmoothPhase is an option for PhaseCalc. Defines how the fasemap is smoothed. Default setting is \"Smooth\". Only works when a mask is also given as input. 
Possible values are \"None\", \"Mask\", \"Median\", \"Smooth\", \"Grow\""

Mara::usage = 
"Mara is an option for PhaseCalc. When True it uses a different phase unwrapping and phasemap calculation approach to cope with two legs. Default value is False."

PhaseCorrect::usage = 
"PhaseCorrect is an option for PhaseCalc. Sometimes the enitre dataset is unwraped to the wrong baseline. 
Shifts the entire phasemap with the given value. Default value is 0." 

MonitorCalc::usage = 
"MonitorCalc is an option for all Calc fucntions. When true the proceses of the calculation is shown."

FullOutput::usage = 
"FullOutput is an option for TensorCalc when using bvector. When True also the S0 is given as output."

RobustFit::usage = 
"RobustFit is an option for TensorCalc. If true outliers will be rejected in the fit, only works with WLLS.
If FullOutput is given the outlier map is given.";

RobustFitParameters::usage =
"RobustFitParameters is an option for TensorCalc. gives the threshold for stopping the itterations and the kappa for the outlier marging, {tr,kappa}."

RejectMap::usage = 
"RejectMap is an option for EigenvalCalc. If Reject is True and RejectMap is True both the eigenvalues aswel as a map showing je rejected values is returned."

Reject::usage = 
"Reject is an option for EigenvalCalc. It True then voxels with negative eigenvalues are rejected and set to 0."

Distribution::usage = 
"Distribution is an option for AngleCalc. values can be \"0-180\", \"0-90\" and \"-90-90\"."

MeanRes::usage = 
"MeanRes is an option for ResidualCalc. When True the root mean square of the residual is calculated."

NormResidual::usage = 
"NormResidual is an option for ResidualCalc. When True the residuals are normalize to the S0 image."

FilterShape::usage = 
"FilterShape is an option for SigmaCalc. Can be \"Gaussian\" of \"Median\"."

FitFunction::usage = 
"FitFunction is an option for ParameterFit. Options are \"Normal\" or \"SkewNormal\". Indicates which function wil be fitted."

FitOutput::usage = 
"FitOutput is an option for ParameterFit and ParameterFit2. Option can be \"Parameters\", \"Function\" or \"BestFitParameters\"."

UseMask::usage = 
"UseMask is a function for MeanSignal and DriftCorrect"

OutputSNR::usage = 
"OutputSNR is an option for SNRMapCalc."

SmoothSNR::usgae = 
"SmoothSNR is an option for SNRMapCalc"

SeedDensity::usage = 
"SeedDensity is an option for FiberDensityMap. The seedpoint spacing in mm."

MeanMethod::usage = 
"MeanMethod is an option for GetMaskMeans. The option can be  \"NormalDist\", \"SkewNormalDist\", or \"Mean\"."

(* ::Subsection::Closed:: *)
(*Error Messages*)


TensorCalc::grad = "The `2` gradient directions defined do not match the `1` gradients directions in the data set."

TensorCalc::data =
"Data set dimensions (`1`D) unknown, posibilities:
- Multiple slices (4D)-> {slices, directions, x,y}
- Single slice (3D)-> {directions, x, y}
- Multiple voxels (2D)-> {directions, voxels}
- Single voxel (1D)-> {directions}"

TensorCalc::bvec = "the `1` gradient directions do not match the `2` b values in the b vector."

TensorCalc::met = "The method specified (`1`) is not a valid method, please use: \"LLS\",\"WLLS\",\"NLS\"."

PhaseCalc::inp =  "The given data is nog a 4D array, or the mask is not a 3D array or the dimensions of the data and the mask are not the same \
Dimensions data: `1`; Dimensions mask: `2`."

ParameterFit::func = "Unknow fit function: `1`. options are SkewNormal or Normal."

ParameterFit::outp = "Unknow output format: `1`. options are Parameters or Function."

ResidualCalc::datdim = "DTIdata (`1`) and tensor data (`2`) are not the same dimensions."

AngleCalc::dist = "Unknown option (`1`), options can be. \"0-180\", \"0-90\" or \"-90-90\"."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*TensorCalc*)


(* ::Subsubsection::Closed:: *)
(*TensorCalc*)


Options[TensorCalc]= {MonitorCalc->True, Method->"iWLLS", FullOutput->False, RobustFit->True,Parallelize->True , RobustFitParameters->{10^-4,6}};

SyntaxInformation[TensorCalc] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

(*bvalue, only one number, gr does not have b=0*)
TensorCalc[data_,gr_,bvalue:_?NumberQ,opts:OptionsPattern[]]:=
Block[{depthD,dirD,dirG,grad,bvec},
	
	depthD=ArrayDepth[data];
	dirD=If[depthD==4,Length[data[[1]]],Length[data]];
	dirG=Length[gr];
	
	(*check if data is 4D, 3D, 2D or 1D*)
	If[depthD>4,Return[Message[TensorCalc::data,ArrayDepth[data]]]];
	(*check if gradient dimensions are the same in the data and grad vector*)
	If[(dirD-1)!=dirG,Return[Message[TensorCalc::grad,dirD,dirG]]];
	
	bvec=Prepend[ConstantArray[bvalue,{dirG}],0];
	grad=N[Prepend[gr,{0,0,0}]];
	
	If[OptionValue[Method]!="DKI",
		TensorCalc[data,Bmatrix[bvec,grad],opts],
		TensorCalc[data,Bmatrix[bvec,grad,Method->"DKI"],opts]
	]
]

(*bvector*)
TensorCalc[data_,grad_,bvec:{_?NumberQ ..},opts:OptionsPattern[]]:=
Block[{depthD,dirD,dirG,dirB},
	
	depthD=ArrayDepth[data];
	dirD=If[depthD==4,Length[data[[1]]],Length[data]];
	dirG=Length[grad];
	dirB=Length[bvec];
	
	(*check if data is 4D, 3D, 2D or 1D*)
	If[depthD>4,Return[Message[TensorCalc::data,ArrayDepth[data]]]];
	(*check if gradient dimensions are the same in the data and grad vector*)
	If[dirD!=dirG,Return[Message[TensorCalc::grad,dirD,dirG]]];
	(*check if bvec is the same lengt as gradient vec*)
	If[dirB!=dirG,Return[Message[TensorCalc::bvec,dirG,dirB]]];
	
	If[OptionValue[Method]!="DKI",
		TensorCalc[data,Bmatrix[bvec,grad],opts],
		TensorCalc[data,Bmatrix[bvec,grad,Method->"DKI"],opts]
	]
]


(*bmatrix*)
TensorCalc[dat_,bmat:{_?ListQ ..},OptionsPattern[]]:=
Block[{dirD,dirB,tensor,rl,rr,TensMin,out,tenscalc,x,data,depthD,xx,bmatI,fout,method,output,robust,dataL,func,con,kappa},
	
	(*get output form*)
	output=OptionValue[FullOutput];
	robust=OptionValue[RobustFit];
	{con,kappa}=OptionValue[RobustFitParameters];
	
	(*chekc method*)
	method=OptionValue[Method];
	If[!MemberQ[{"LLS","WLLS","iWLLS"(*,"NLS","GMM","CLLS","CWLLS","CNLS","DKI"*)},method],
		Return[Message[TensorCalc::met,method];$Failed]
		];
	
	data=N[Clip[Chop[dat],{0.,Infinity}]];
	(*get the data dimensions*)	
	depthD=ArrayDepth[data];
	dirD=If[depthD==4,Length[data[[1]]],Length[data]];
	dirB=Length[bmat];
	
	(*check if data is 4D, 3D, 2D or 1D*)
	If[depthD>4,Return[Message[TensorCalc::data,ArrayDepth[data]];$Failed]];
	(*check if bmat is the same lengt as data*)
	If[dirB!=dirD,Return[Message[TensorCalc::bvec,dirD,dirB];$Failed]];
	
	(*define data*)
	dataL=LogNoZero[Chop[data]];	
	(*calculate the inverse bmat*)
	bmatI=PseudoInverse[bmat];
	
	(*if data is 4D handle as multiple 3D sets (saves memory and calculation time)*)
	If[depthD==4,
		
		xx=0;
		
		func=If[OptionValue[Parallelize],
			SetSharedVariable[xx,data,dataL];
			DistributeDefinitions[data,dataL,bmat,bmatI,method,output,TensorCalci];
			ParallelTable,
			Table
		];	
		
		If[OptionValue[MonitorCalc],PrintTemporary[ProgressIndicator[Dynamic[xx], {0, Length[data]}]]];
		tensor = func[
			xx++;
			TensorCalci[data[[x]],dataL[[x]],bmat,bmatI,Method->method,FullOutput->output,RobustFit->robust,RobustFitParameters->{con,kappa}]
			,{x,1,Length[data],1}];
		
		(*full output returns {tens,S0,(outliers),residuals}*)
		If[output,
			tensor = Transpose[tensor];
			tensor[[1]] = Transpose[tensor[[1]]]
			,
			tensor = Transpose[tensor]
		];	
		
		,(*1D,2D,3D*)
		tensor=TensorCalci[data,dataL,bmat,bmatI,Method->method,FullOutput->output,RobustFit->robust,RobustFitParameters->{con,kappa}];
	];
	
	tensor
]


(* ::Subsubsection::Closed:: *)
(*TensorCalci*)


Options[TensorCalci] = Options[TensorCalc];

TensorCalci[data_, dataL_, bmat_, bmatI_,OptionsPattern[]]:=Block[
	{l,r,depthD,TensMin,tensor,w, method,outliers,S0,fitresult,robust,residual,con,kappa},
	
	(*transpose the data*)
	depthD = ArrayDepth[data];
	l = RotateRight[Range[depthD]];
	r = RotateLeft[Range[depthD]];
	
	method = OptionValue[Method];
	robust = (OptionValue[RobustFit] && method =!= "LLS");
	{con,kappa}=OptionValue[RobustFitParameters];
	
	outliers = If[robust,Transpose[FindOutliers[Transpose[dataL,l], bmat, con, kappa], r], ConstantArray[0.,Dimensions[data]]];
	
	fitresult = Switch[method,
		"LLS", Transpose[TensMinLLS[Transpose[dataL,l], bmatI], r],
		"WLLS", Transpose[TensMinWLLS[Transpose[(1-outliers) data,l],Transpose[(1-outliers) dataL,l], bmat], r],
		"iWLLS", Transpose[TensMiniWLLS[Transpose[(1-outliers) data,l],Transpose[(1-outliers) dataL,l], bmat], r]
		];
	
	If[OptionValue[FullOutput],residual = ResidualCalc[data,fitresult,outliers,bmat,MeanRes->"MAD"]];
		
	S0 = ExpNoZero[Last[fitresult]];
	tensor = Clip[Drop[fitresult,-1],{-0.1,0.1}];

	If[OptionValue[FullOutput],If[robust,{tensor,S0,outliers,residual},{tensor,S0,residual}],tensor]
]


(* ::Subsubsection:: *)
(*FindOutliers*)


FindOutliers = Block[{ittA,itt,contA,cont,sol,solA,soli,res,weigths, wmat,fitE,LS2,bmat2,mad,out}, 
  	Compile[{{LS, _Real, 1}, {bmat, _Real, 2}, {con, _Real, 0}, {kappa, _Real, 0}},
  		If[AllTrue[LS, 0. === # &]||Total[Unitize[LS]]<7,
  			(*skip if background*)
  			out = 0. LS;
  			,
  			(*Find the outliers*)
  			(*initialize*)
  			ittA = 0; contA = 1;
  			(*Step1: initial LLS fit*)
  			sol = LeastSquares[bmat,LS];
  			
  			(*check if LLS fit is plausable*)
  			If[Negative[Last[sol]],
  				out = 0. LS;
  				,
	  			While[contA == 1,(*init itteration values*)
	  				ittA++;
	  				solA = sol;
	  				
	  				itt = 0; cont = 1;
	  				(*Step2: Compute a robust estimate for homoscedastic regression using IRLS.*)
	  				While[cont == 1,
	  					itt++;
	  					soli = sol;
	  					(*a. Calculate the residuals e* in the linear domain*)
	  					res = LS - bmat.sol;
	  					(*b. Obtain an estimate of the dispersion of the residuals by calculating the median absolute deviation (MAD).*)
	  					mad = 1.4826 MedianDeviation[res];
	  					(*prevent calculation with 0*)
	  					If[AllTrue[res, (0. === #) &] || mad === 0.,
	  						cont = 0,
	  						(*c. Recompute the weights according to Eq. [13].*)
	  						weigths = 1/(1 + (res/mad)^2)^2;
	  						(*d. Perform WLLS fit with new weights*)
	  						wmat = DiagonalMatrix[weigths];
	  						sol = LeastSquares[Transpose[bmat].wmat.bmat, Transpose[bmat].wmat.LS];
	  						(*e. Check convergence*)
	  						If[!AnyTrue[Abs[sol - soli] - con (Max /@ Transpose[{Abs[sol], Abs[soli]}]),Positive] || itt === 5, cont = 0];
	  					];
	  				];(*end first while*)
	   
					itt = 0; cont = 1;
	   
					(*Step 3: Transform variables for heteroscedasticity*)
					fitE = Exp[bmat.sol]+10^-6;
					LS2 = LS / fitE;
					bmat2 = bmat / fitE;
									
					(*Step 4: Initial LLS fit in * domain*)
					sol = LeastSquares[bmat2,LS2];
	   
					(*Step 5: Compute a robust estimate for homoscedastic regression using IRLS.*)
					While[cont == 1,
						itt++;
						soli = sol;
						(*a. Calculate the residuals e* in the linear domain*)
						res = LS2 - bmat2.sol;
						(*b. Obtain an estimate of the dispersion of the residuals by calculating the median absolute deviation (MAD).*)
						mad = 1.4826 MedianDeviation[res];
						(*prevent calculation with 0*)
						If[AllTrue[res, (0. === #) &] || mad === 0.,
							cont = 0,
							(*c. Recompute the weights according to Eq. [13].*)
							weigths = 1/(1 + (res/mad)^2)^2;
							(*d. Perform WLLS fit with new weights*)
							wmat = DiagonalMatrix[weigths];
							sol = LeastSquares[Transpose[bmat2].wmat.bmat2, Transpose[bmat2].wmat.LS2];
							(*e. Check convergence*)
							If[!AnyTrue[Abs[sol - soli] - con (Max /@ Transpose[{Abs[sol], Abs[soli]}]),Positive] || itt === 5, cont = 0];
						];
					];(*end second while*)
					
					(*Step 6: Check convergence overall loop*)
					If[! AnyTrue[Abs[sol - solA] - con (Max /@ Transpose[{Abs[sol], Abs[solA]}]), Positive] || ittA === 10, contA = 0];
				];(*end main while*)
  			
				
				(*Step 7: Identify and exclude outliers*)
				res = LS2 - bmat2.sol;
				out = UnitStep[Abs[res] - (kappa 1.4826 MedianDeviation[res])];
			
			];(*close if negative S0*)
		];(*close if background*)
		
		out
		
		,{{wmat, _Real, 2}, {bmat2, _Real, 2}, {out, _Real, 1}},
		RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
	](*close compile*)
];


(* ::Subsubsection::Closed:: *)
(*LLS*)


TensMinLLS = Compile[{{LS, _Real, 1}, {bmatI, _Real, 2}},
	bmatI.LS,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*WLLS*)


TensMinWLLS = Block[{wmat,mvec,sol},
	Compile[{{S, _Real, 1},{LS, _Real, 1},{bmat, _Real, 2}}, 
	    If[AllTrue[LS, 0. === # &]||Total[Unitize[LS]]<7,
	    	sol = ConstantArray[0., Length@First@bmat]
	    	,
	    	mvec = UnitStep[LS] Unitize[LS]; (*if 0 then it is not used because w=0*)
	    	wmat = DiagonalMatrix[mvec S^2];
	    	sol = LeastSquares[Transpose[bmat].wmat.bmat,Transpose[bmat].wmat.LS];
	    ];
	    sol
    ,{{wmat,_Real,2}, {sol, _Real, 1}}, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]];


(* ::Subsubsection::Closed:: *)
(*iWLLS*)


TensMiniWLLS = Block[{wmat, mat, cont, itt, mvec, soli, max, sol, w},
   Compile[{{S, _Real, 1}, {LS, _Real, 1}, {bmat, _Real, 2}},
    mvec = UnitStep[S] Unitize[S];
    max = Max[mvec S];
    If[AllTrue[LS, 0. === # &] || Total[mvec] <= 7,
     (*skip background or not enough data for fit*)
     sol = 0. First@bmat;
     ,
     (*initialize*)itt = 0;
     cont = 1;
     (*initialize using LLS*)
     sol = LeastSquares[bmat, LS];
     (*check for implausabole solution (negative S0)*)
     If[Last[sol] >= 3*max || Last[sol] <= 0,
      sol = 0. First@bmat;
      ,
      (*itterative reweighting*)
      While[cont == 1,
       (*init itteration values*)
       itt++;
       soli = sol;
       (*perform WLLS*)
       w = (mvec Exp[bmat.sol])^2;
       wmat = DiagonalMatrix[w];
       sol = LeastSquares[Transpose[bmat].wmat.bmat, Transpose[bmat].wmat.LS];
       (*update weight*)
       (*see if to quit loop*)
       If[ Last[sol] >= 3*max || Last[sol] <= 0, cont=0;  sol = 0. First@bmat];
       If[! AnyTrue[Abs[sol - soli] - 0.0001 (Max /@ Transpose[{Abs[sol], Abs[soli]}]), Positive] || itt === 10 , cont = 0];
       ](*close while*)
      ];(*close if S0*)
     ];(*close if back*)
    sol, 
    {{mat, _Real, 2}, {wmat, _Real, 2}, {sol, _Real, 1}},
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]];
        

(* ::Subsubsection::Closed:: *)
(*DKI*)


(*TensMinDKI[S_,LS_,bmat_,bmatI_]:=bmatI.LS*)
TensMinDKI = Compile[{{S, _Real, 1}, {bmatI, _Real, 2}},
	If[Total[S]==0.,
    	{0.,0.,0.,0.,0.,0.,0.},
    	bmatI.S
	],RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"(*, Parallelization -> True*)];


(* ::Subsubsection::Closed:: *)
(*NLS*)


TensMinNLS[S_,LS_,bmat_,bmatI_]:=
Module[{v,xx,yy,zz,xy,xz,yz,init,tens,sol},
	tens=bmatI.LS;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},
		tens,
		v={xx,yy,zz,xy,xz,yz,tens[[7]]};
		init=Thread[{v[[1;;6]],tens[[1;;6]]}];
		sol=FindMinimum[.5 Total[(S-Exp[bmat.v])^2],init][[2]];
		v/.sol
	]
]


(* ::Subsubsection::Closed:: *)
(*NLS*)


TensMinGMM[S_,LS_,bmat_,bmatI_]:=
Module[{v,xx,yy,zz,xy,xz,yz,init,tens,res,w},
	S;
	tens=bmatI.LS;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={xx,yy,zz,xy,xz,yz,tens[[7]]};
		init=Thread[{v[[1;;6]],tens[[1;;6]]}];
		v/.FindMinimum[(
			res=LS-bmat.v;
			w=1/(res^2+Mean[res]^2);
			.5 Total[(w/Mean[w])*(res)^2]
		(*w=1/(res^2+(1.4826*Median[Abs[res-Median[res]]])^2);*)
		),init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CLLS*)


TensMinCLLS[S_,LS_,bmat_,bmatI_]:=
Module[{v,R0,R1,R2,R3,R4,R5,init,tens},
	S;
	tens=bmatI.LS;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={R0^2,R1^2+R3^2,R2^2+R4^2+R5^2,R0 R3,R0 R4,R3 R4+R1 R5,tens[[7]]};
		init=Thread[{{R0,R1,R2,R3,R4,R5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(LS-bmat.v)^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CWLLS*)


TensMinCWLLS[S_,LS_,bmat_,bmatI_]:=
Module[{v,R0,R1,R2,R3,R4,R5,init,tens,std=1,wmat},
	bmatI;
	wmat=Transpose[bmat].DiagonalMatrix[S^2/std^2];
	tens=PseudoInverse[wmat.bmat].wmat.LS;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={R0^2,R1^2+R3^2,R2^2+R4^2+R5^2,R0 R3,R0 R4,R3 R4+R1 R5,tens[[7]]};
		init=Thread[{{R0,R1,R2,R3,R4,R5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(S^2/std^2)*(LS-bmat.v)^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CNLS*)


TensMinCNLS[S_,LS_,bmat_,bmatI_]:=
Module[{v,R0,R1,R2,R3,R4,R5,init,tens},
	tens=bmatI.LS;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={R0^2,R1^2+R3^2,R2^2+R4^2+R5^2,R0 R3,R0 R4,R3 R4+R1 R5,tens[[7]]};
		init=Thread[{{R0,R1,R2,R3,R4,R5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(S-Exp[bmat.v])^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*ExtendeCholeskyDecomposition*)


ExtendedCholeskyDecomposition[Tm_]:=
Module[{n,beta,theta,Cm,Lm,Dm,Em,j},
	n=Length[Tm];
	beta=Max[{Max[Diagonal[Tm]],Max[UpperTriangularize[Tm,1]]/Sqrt[n^2-1],10^-15}];
	Cm=DiagonalMatrix[Diagonal[Tm]];
	Lm=Dm=Em=ConstantArray[0,{n,n}];
	For[j=1,j<=3,j++,
		If[j==1,
			(*j=1 maak eerste colom Cm gelijk aan Tm*)
			Cm[[j+1;;,j]]=Tm[[j+1;;,j]];
			,
			(*j>1 vul Lm matrix*)
			Lm[[j,;;j-1]]=Cm[[j,;;j-1]]/(Diagonal[Dm][[;;j-1]]/.(0.->Infinity));
			If[j<n,
				Cm[[j+1;;,j]]=Tm[[j+1;;,j]]-Lm[[j,j-1;;]].Transpose[Cm[[j+1;;,j-1;;]]]
				];
			];
		theta=If[j==n,0,Max[Abs[Cm[[j+1;;,j]]]]];
		Dm[[j,j]]=Max[{Abs[Cm[[j,j]]],theta^2/beta}];
		Em[[j,j]]=Dm[[j,j]]-Cm[[j,j]];
		Cm=Cm-DiagonalMatrix[PadLeft[(1/(Dm[[j,j]]/.(0.->Infinity)))*Cm[[j+1;;,j]]^2,n]];
		];
	Lm=Lm+IdentityMatrix[n];
	Transpose[Lm.MatrixPower[Dm,.5]]
	]


(* ::Subsection::Closed:: *)
(*EigenvalCalc*)


(* ::Subsubsection::Closed:: *)
(*EigenvalCalc*)


Options[EigenvalCalc]= {MonitorCalc->True,RejectMap->False,Reject->True};

SyntaxInformation[EigenvalCalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

EigenvalCalc[tensor_?ArrayQ,OptionsPattern[]]:=
Module[{output},
	Switch[ArrayDepth[tensor],
		4,
		output=EigenvalCalci[Transpose[tensor,{4,1,2,3}],RejectMap->OptionValue[RejectMap],Reject->OptionValue[Reject]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for ",Length[tensor[[1]]]," slices!"]];,
		3,
		output=EigenvalCalci[Transpose[tensor,{3,1,2}],RejectMap->OptionValue[RejectMap],Reject->OptionValue[Reject]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 slice!"]];,
		1,
		output=EigenvalCalci[tensor,RejectMap->OptionValue[RejectMap],Reject->OptionValue[Reject]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 voxel!"]];	
	];
		
	(*If[ArrayQ[tensor,4],
		output=EigenvalCalci[Transpose[tensor,{4,1,2,3}],RejectMap->OptionValue[RejectMap],Reject->OptionValue[Reject]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for ",Length[tensor[[1]]]," slices!"]];,
		If[ArrayQ[tensor,3],
			output=EigenvalCalci[Transpose[tensor,{3,1,2}],RejectMap->OptionValue[RejectMap],Reject->OptionValue[Reject]];
			If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 slice!"]];,
			If[VectorQ[tensor],
				output=EigenvalCalci[tensor,RejectMap->OptionValue[RejectMap],Reject->OptionValue[Reject]];
				If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 voxel!"]];
				]
			]		
		];*)
	Return[N[output]];
	]


(* ::Subsubsection::Closed:: *)
(*EigenvalCalci*)


Options[EigenvalCalci]={RejectMap->False,Reject->True}

EigenvalCalci[tensor_, OptionsPattern[]] := Module[{rejectmap, values},
  (*calculate eigenvalues*)
  values = Map[Eigenvalues2, tensor, {ArrayDepth[tensor] - 1}];
  (*reject negative eigenvalues*)
  values = If[OptionValue[Reject], RejectEig[values], values];
  (*create rejectmap*)
  If[OptionValue[RejectMap],
  rejectmap = Map[RejectMapi, values, {ArrayDepth[tensor] - 1}];
   {values, rejectmap}
   ,
   values
   ]
  ]

Eigenvalues2[{0., 0., 0., 0., 0., 0.}] := {0., 0., 0.};
Eigenvalues2[{0, 0, 0, 0, 0, 0}] := {0., 0., 0.};
Eigenvalues2[tens_] := Eigenvalues[{{tens[[1]],tens[[4]],tens[[5]]},{tens[[4]],tens[[2]],tens[[6]]},{tens[[5]],tens[[6]],tens[[3]]}}];
RejectEig = Compile[{{eig, _Real, 1}}, 
	If[Total[1 - UnitStep[eig]] > 0, {0., 0., 0.}, eig], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];
RejectMapi[{0., 0., 0.}] := 0;
RejectMapi[{_, _, _}] := 1;


(* ::Subsection::Closed:: *)
(*EigenvecCalc*)


(* ::Subsubsection::Closed:: *)
(*EigenvecCalc*)


Options[EigenvecCalc]= {MonitorCalc->False};

SyntaxInformation[EigenvecCalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

EigenvecCalc[tensor_?ArrayQ,OptionsPattern[]]:=
Module[{output},
	
	Switch[ArrayDepth[tensor],
		4,
		output=If[OptionValue[MonitorCalc],
			Monitor[EigenvecCalci[Transpose[tensor,{4,1,2,3}]],ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate]],
			EigenvecCalci[Transpose[tensor,{4,1,2,3}]]
			];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for "<>ToString[Length[tensor[[1]]]]<>" slices!"]];,
		3,
		output=EigenvecCalci[Transpose[tensor,{3,1,2}]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 slice!"]];,
		1,
		output=EigenvecCalci[tensor];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 voxel!"]];	
	];
	(*
	If[ArrayQ[tensor,4],
		Monitor[
			output=EigenvecCalci[Transpose[tensor,{4,1,2,3}]];,If[OptionValue[MonitorCalc],ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate],""]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for "<>ToString[Length[tensor[[1]]]]<>" slices!"]];,
		If[ArrayQ[tensor,3],
			output=EigenvecCalci[Transpose[tensor,{3,1,2}]];
			If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 slice!"]];,
			If[VectorQ[tensor],
				output=EigenvecCalci[tensor];
				If[OptionValue[MonitorCalc],Print["Done calculating eigenvalues for 1 voxel!"]];
				]
			]
		];*)
	Return[output]
	]


(* ::Subsubsection::Closed:: *)
(*EigenvecCalci*)


EigenvecCalci[tensor_]:=
Map[
	If[#[[1]]==#[[2]]==#[[3]]==#[[4]]==#[[5]]==#[[6]],
		{{0,0,1},{0,1,0},{1,0,0}},
		Eigenvectors[{{#[[1]],#[[4]],#[[5]]},{#[[4]],#[[2]],#[[6]]},{#[[5]],#[[6]],#[[3]]}}]
		]&,tensor,{ArrayDepth[tensor]-1}
	]


(* ::Subsection::Closed:: *)
(*EigensysCalc*)


(* ::Subsubsection::Closed:: *)
(*EigensysCalc*)


Options[EigensysCalc]= {MonitorCalc->False};

SyntaxInformation[EigensysCalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

EigensysCalc[tensor_?ArrayQ,OptionsPattern[]]:=
Module[{output},
	
	Switch[ArrayDepth[tensor],
		4,
		output=If[OptionValue[MonitorCalc],
			Monitor[EigensysCalci[Transpose[tensor,{4,1,2,3}]],
				ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate]],
			EigensysCalci[Transpose[tensor,{4,1,2,3}]]
			];
		If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for "<>ToString[Length[tensor[[1]]]]<>" slices!"]];,
		3,
		output=EigensysCalci[Transpose[tensor,{3,1,2}]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for 1 slice!"]];,
		2,
		output=EigensysCalci[Transpose[tensor]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for Row of voxels!"]];,
		1,
		output=EigensysCalci[tensor];
		If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for 1 voxel!"]];	
	];
	(*
	If[ArrayQ[tensor,4],
		output=EigensysCalci[Transpose[tensor,{4,1,2,3}]];
		If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for "<>ToString[Length[tensor[[1]]]]<>" slices!"]];,
		If[ArrayQ[tensor,3],
			output=EigensysCalci[Transpose[tensor,{3,1,2}]];
			If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for 1 slice!"]];,
			If[VectorQ[tensor],
				output=EigensysCalci[tensor];
				If[OptionValue[MonitorCalc],Print["Done calculating eigensystem for 1 voxel!"]];
				]
			]
		];*)
	Return[output]
	]


(* ::Subsubsection::Closed:: *)
(*EigensysCalc*)


EigensysCalci[tensor_]:=
Module[{eig},
	Map[
		If[#[[1]]==#[[2]]==#[[3]]==#[[4]]==#[[5]]==#[[6]],
			{{0,0,0},{{0,0,1},{0,1,0},{1,0,0}}},
			eig=Eigensystem[TensMat[#]];
			If[Positive[eig[[1,1]]]&&Positive[eig[[1,2]]]&&Positive[eig[[1,3]]],eig,{{0,0,0},{{0,0,1},{0,1,0},{1,0,0}}}]
			]&,tensor,{ArrayDepth[tensor]-1}
		]
	]


(* ::Subsection::Closed:: *)
(*ADCCalc*)


SyntaxInformation[ADCCalc] = {"ArgumentsPattern" -> {_}};

ADCCalc[eig_] := ADCCalci[eig]

ADCCalci = Compile[{{eig, _Real, 1}}, Mean[eig], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsection::Closed:: *)
(*FACalc*)


SyntaxInformation[FACalc] = {"ArgumentsPattern" -> {_}};

FACalc[eig_] := FACalci[eig]

FACalci = Compile[{{eig, _Real, 1}}, Block[{l1, l2, l3, teig},
   l1 = eig[[1]]; l2 = eig[[2]]; l3 = eig[[3]];
   teig = Sqrt[2.*Total[eig^2]];
   If[teig == 0., 0. ,Sqrt[(l1 - l2)^2 + (l2 - l3)^2 + (l1 - l3)^2]/teig]
   ], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection::Closed:: *)
(*ECalc*)


Options[ECalc]= {MonitorCalc->True};

SyntaxInformation[ECalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ECalc[eigen_,OptionsPattern[]]:=
Module[{output,slices,x},
	slices=Length[eigen];
	If[ArrayQ[eigen,4],
		Monitor[output=Table[ECalci[eigen[[x]]],{x,1,slices,1}];,If[OptionValue[MonitorCalc],Column[{"Calculation E for Multiple slices",ProgressIndicator[x,{0,slices}]}],""]];
		If[OptionValue[MonitorCalc],Print["Done calculating e for "<>ToString[slices]<>" slices!"]];
		,
		output=ECalci[eigen];
		If[ArrayQ[eigen,3],If[OptionValue[MonitorCalc],Print["Done calculating e for 1 slice!"]],
			If[VectorQ[eigen],If[OptionValue[MonitorCalc],Print["Done calculating e for 1 voxel!"]]]
			]
		];
	Return[output];
	]


ECalci[eigen_]:=
Module[{EC},
	EC=Compile[{l1,l3},Sqrt[1-(l3/l1)]];
	Map[If[#[[3]]!=0,
		EC[#[[1]],#[[3]]],
		0]&,eigen,{ArrayDepth[eigen]-1}]
	]


(* ::Subsection::Closed:: *)
(*ParameterCalc*)


Options[ParameterCalc]={Reject->True,MonitorCalc->False}

SyntaxInformation[ParameterCalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ParameterCalc[tensor_,OptionsPattern[]]:=
Module[{eig,adc,fa},
	eig=1000.*EigenvalCalc[tensor,Reject->OptionValue[Reject],MonitorCalc->OptionValue[MonitorCalc]];
	adc=ADCCalc[eig];
	fa=FACalc[eig];
	Join[TransData[eig,"r"],{adc,fa}]
	]


(* ::Subsection::Closed:: *)
(*PhaseCalc*)


Options[PhaseCalc]={SmoothPhase->"Smooth",BackgroundFilter->6,MonitorUnwrap->True,UnwrapDimension->"2D"};

SyntaxInformation[PhaseCalc] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

PhaseCalc[dat_?ArrayQ, opts : OptionsPattern[]] :=PhaseCalc[dat, 1, opts]

PhaseCalc[dat_?ArrayQ, mask_,OptionsPattern[]] :=
Module[{B0data, B0unw, phase, data},
  If[mask != 1 && Dimensions[dat[[All, 1]]] != Dimensions[mask], 
   Return["error"]];
   data=If[Min[dat]>=-Pi&&Max[dat]<=Pi,
   	Transpose[dat],
   	N[(((1.53455433455433 Transpose[dat]) - 3142)/1000)]
   ];
  B0data = mask*(data[[2]] - data[[1]]);
  B0unw = Unwrap[B0data, BackgroundFilter->OptionValue[BackgroundFilter],MonitorUnwrap->OptionValue[MonitorUnwrap],UnwrapDimension->OptionValue[UnwrapDimension]];
  phase = 
   Switch[OptionValue[SmoothPhase], 
   	"None", Chop[B0unw], 
   	"Median",Chop[MedianFilter[B0unw, 1]], 
    "Smooth", Chop[GaussianFilter[MedianFilter[B0unw, 1], 3]]];
  Return[phase]]


(* ::Subsection::Closed:: *)
(*AngleCalc*)


Options[AngleCalc]={Distribution->"0-180"};

SyntaxInformation[AngleCalc] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};


AngleCalc[data_?ArrayQ,vec_?VectorQ,OptionsPattern[]]:=
Module[{angles},
	angles=Map[If[Re[#]==#,ArcCos[#.vec],"no"]&,data,{Depth[data]-2}];
	
	Switch[
		OptionValue[Distribution],
		"0-180",
		(180/Pi)angles,
		"0-90",
		Map[If[#!="no",If[#>=1/2Pi,(180/Pi)(Pi-#),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
		"-90-90",
		Map[If[#!="no",If[#>=1/2Pi,(180/Pi)(#-Pi),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
		_,
		Message[AngleCalc::dist,OptionValue[Distribution]]
		]
	]


AngleCalc[data_?ArrayQ,vec_?ArrayQ,OptionsPattern[]]:=
Module[{angles},
	If[Dimensions[data]!=Dimensions[vec],
		Print["Error"],
		
		angles=MapThread[If[Re[#1]==#1,ArcCos[#1.#2],"no"]&,{data,vec},ArrayDepth[vec]-1];
		
		Switch[
			OptionValue[Distribution],
			"0-180",
			(180/Pi)angles,
			"0-90",
			Map[If[#!="no",If[#>=1/2Pi,(180/Pi)(-(#-Pi)),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
			"-90-90",
			Map[If[#!="no",If[#>=1/2Pi,(180/Pi)(#-Pi),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
			_,
			Message[AngleCalc::dist,OptionValue[Distribution]]
			]
		]
	]


(* ::Subsection::Closed:: *)
(*AngleMap*)


SyntaxInformation[AngleMap] = {"ArgumentsPattern" -> {_}};

AngleMap[vec_]:=
Module[{az,zen},
	Transpose[Map[If[#=={0,0,1},
		{0,0},
		If[Negative[#[[3]]],v=-#,v=#];
		az=ArcCos[v[[3]]]/Degree;
		zen=ArcTan[v[[2]],v[[1]]]/Degree;
		(*zen=If[zen<-90,zen+180,If[zen>90,zen-180,zen]];*)
		zen=If[Negative[zen],zen+180,zen];
		{az,zen}
		]&,vec,{3}],{2,3,4,1}]
	]


(* ::Subsection::Closed:: *)
(*ResidualCalc*)


Options[ResidualCalc] = {MeanRes -> "All"};

SyntaxInformation[ResidualCalc] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

(*b0, no outliers, bval, bvec*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, S0_?ArrayQ}, grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}, bval_, opts : OptionsPattern[]] := 
 ResidualCalc[data, Join[tensor, {LogNoZero[S0]}], ConstantArray[0., Dimensions[data]], Bmatrix[bval, grad], opts]

(*b0, no outliers bmat*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, S0_?ArrayQ}, bmat_?ArrayQ, opts : OptionsPattern[]] :=
 ResidualCalc[data, Join[tensor, {LogNoZero[S0]}], ConstantArray[0., Dimensions[data]], bmat, opts]

(*b0, outliers, bval, bvec*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, S0_?ArrayQ}, outlier_?ArrayQ, grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}, bval_, opts : OptionsPattern[]] :=
 ResidualCalc[data, Join[tensor, {LogNoZero[S0]}], outlier, Bmatrix[bval, grad], opts]

(*b0, outliers, bmat*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, S0_?ArrayQ}, outlier_?ArrayQ, bmat_?ArrayQ, opts : OptionsPattern[]] :=
 ResidualCalc[data, Join[tensor, {LogNoZero[S0]}], outlier, bmat, opts]

(*no outliers, bval, bvec*)
ResidualCalc[data_?ArrayQ, tensor_?ArrayQ, grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}, bval_, opts : OptionsPattern[]] := 
 ResidualCalc[data, tensor, ConstantArray[0., Dimensions[data]], Bmatrix[bval, grad], opts]

(*no outliers bmat*)
ResidualCalc[data_?ArrayQ, tensor_?ArrayQ, bmat_?ArrayQ, opts : OptionsPattern[]] :=
 ResidualCalc[data, tensor, ConstantArray[0., Dimensions[data]], bmat, opts]

(*outliers, bval, bvec*)
ResidualCalc[data_?ArrayQ, tensor_?ArrayQ, outlier_?ArrayQ, grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}, bval_, opts : OptionsPattern[]] :=
 ResidualCalc[data, tensor, outlier, Bmatrix[bval, grad], opts]

ResidualCalc[data_?ArrayQ, tensor_?ArrayQ, outlier_?ArrayQ, bmat_?ArrayQ, OptionsPattern[]] := Block[{
	fit, dat, err, dimD,dimT
	},
  dat = N[data];
  (*check data and tensor dimensions*)
  dimD=Dimensions[If[ArrayDepth[dat] == 4,dat[[All,1]],dat[[1]]]];
  dimT=Dimensions[tensor[[1]]];
  
  If[dimD != dimT || Length[tensor]!=7, Return[Message[ResidualCalc::datdim, Dimensions[data], Dimensions[tensor]]]];

  (*remove ouliers*)
  
  fit = Clip[ExpNoZero[bmat.tensor], {-1.5, 1.5} Max[dat], {0., 0.}];
 
  err = If[ArrayDepth[dat] == 4,
    Transpose[(1 - outlier)] (Transpose[dat] - fit),
    (1 - outlier) (dat - fit)
    ];
  
  Switch[OptionValue[MeanRes],
   "RMSE", RMSNoZero[err],
   "MAD", MADNoZero[err],
   _, If[ArrayDepth[dat] == 4, Transpose[err], err] // N
   ]
  ]


(* ::Subsection:: *)
(*ParameterFit*)


(* ::Subsubsection::Closed:: *)
(*ParameterFit*)


Options[ParameterFit] = {FitFunction -> "SkewNormal", FitOutput -> "Parameters", Method -> Automatic}

SyntaxInformation[ParameterFit] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ParameterFit[dat : {_?ListQ ..}, opts : OptionsPattern[]] := ParameterFit[Flatten[#], opts] & /@ dat

ParameterFit[dat_List, OptionsPattern[]] := Module[{mod, out, met, data, mdat, sdat, fdat},
  
  (*get option values*)
  mod = OptionValue[FitFunction];
  out = OptionValue[FitOutput];
  met = OptionValue[Method];
  
  (*prepare data*)
  data = DeleteCases[Flatten[dat] // N, 0.];
  (*initialization for mean and std*)
  mdat = Mean[data];
  sdat = StandardDeviation[data];
  (*fit data*)
  fdat = FitData[data];
  
  Off[NonlinearModelFit::"cvmit"]; Off[NonlinearModelFit::"sszero"];
  (*perform the fit*)
  If[Length[data] <= 10,
   Print["Not Enough data in the ROI"];
   ,
   Switch[mod,
    (*SkewNormal dist parameter fit*)
    "SkewNormal",
    sol = NonlinearModelFit[fdat,  PDF[SkewNormalDistribution[Mu, Sigma, Alpha], x], {{Mu, mdat}, {Sigma, sdat}, {Alpha, 0}}, x, Method -> met];
    par = {Mu, Sigma, Alpha} /. sol["BestFitParameters"];
    fun = SkewNormalDistribution[Mu, Sigma, Alpha] /. sol["BestFitParameters"];
    ,
    (*Normal dist parameter fit*)
    "Normal",
    sol = NonlinearModelFit[fdat, PDF[NormalDistribution[Mu, Sigma], x], {{Mu, mdat}, {Sigma, mdat/2}}, x];
    par = {Mu, Sigma} /. sol["BestFitParameters"];
    fun = NormalDistribution[Mu, Sigma] /. sol["BestFitParameters"];
    ,
    _,
    Message[ParameterFit::func, mod]]
   ];
  On[NonlinearModelFit::"cvmit"]; On[NonlinearModelFit::"sszero"];
  
  (*generate Output*)
  Switch[out,
   "Parameters",
   par,
   "ParametersExtra",
   Flatten[{Mean[fun], StandardDeviation[fun], Quantile[fun, {.5, .05, .95}]}],
   "Function",
   sol,
   "BestFitParameters",
   sol["BestFitParameters"],
   _, Message[ParameterFit::outp, out]
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*ParameterFit2*)


Options[ParameterFit2]={FitOutput->"BestFitParameters"}

SyntaxInformation[ParameterFit2] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ParameterFit2[dat_List,OptionsPattern[]]:=
Module[{i,datf,init,out,sol,par,
	Omega1i,Omega2i,Alpha1i,Alpha2i,Xi1i,Xi2i,
	Omega1,Xi1,Alpha1,Omega2,Xi2,Alpha2},
	Off[NonlinearModelFit::cvmit];Off[NonlinearModelFit::eit];Off[NonlinearModelFit::sszero];
	
	init={
		{.2,.2,0,0,.6,2.1},
		{.2,.2,0,0,.6,1.6},
		{.2,.2,0,0,.6,1.3},
		{.2,.2,0,0,.6,1.7},
		{.2,.2,0,0,.3,.2}
		};
	datf=DeleteCases[DeleteCases[Flatten[#]//N,0.],1.]&/@dat;
	out=OptionValue[FitOutput];
	i=0;
	sol=MapThread[
		(
		i++;
		{Omega1i,Omega2i,Alpha1i,Alpha2i,Xi1i,Xi2i}=#2;
		NonlinearModelFit[
			FitData[#1,3.5],
			{f SkewNorm[x,Omega1,Xi1,Alpha1]+(1-f)SkewNorm[x,Omega2,Xi2,Alpha2],
			0<=f<=1,
			0<Omega1,
			0<Omega2,
			0<Mn[Omega1,Xi1,Alpha1]<1,
			If[i==5,Mn[Omega1,Xi1,Alpha1]>1.2Mn[Omega2,Xi2,Alpha2],Mn[Omega1,Xi1,Alpha1]<Mn[Omega2,Xi2,Alpha2]],
			If[i==5,-2<Alpha1<0,-1.5<Alpha1<1.5],
			If[i==5,0<Alpha2<2,-1.5<Alpha2<1.5]},
			{{f,0.5},{Omega1,Omega1i},{Omega2,Omega2i},{Alpha1,Alpha1i},{Alpha2,Alpha2i},{Xi1,Xi1i},{Xi2,Xi2i}},
			x,MaxIterations->1000]
		)&,{datf,init}];
	
	par={Mn[Omega1,Xi1,Alpha1],Sqrt[Var[Omega1,Alpha1]]}/.#["BestFitParameters"]&/@sol;
	
	On[NonlinearModelFit::cvmit];On[NonlinearModelFit::eit];On[NonlinearModelFit::sszero];
	
	Switch[
		out,
		"Parameters",
		{Mn[Omega1,Xi1,Alpha1],Sqrt[Var[Omega1,Alpha1]],Mn[Omega2,Xi2,Alpha2],Sqrt[Var[Omega2,Alpha2]]}/.#["BestFitParameters"]&/@sol,
		"Function",
		sol,
		"BestFitParameters",
		{f,Omega1,Omega2,Xi1,Xi2,Alpha1,Alpha2}/.#["BestFitParameters"]&/@sol,
		_,
		Message[ParameterFit::outp,out]
		]
	]


(* ::Subsubsection::Closed:: *)
(*FitData*)


SyntaxInformation[FitData] = {"ArgumentsPattern" -> {_, _.}};

FitData[dat_,sdr_:2]:=
Module[{m, s, min, max, range, step, xdat, data, out}, 
  If[dat == {} || Length[dat] == 1, {}, 
  	m = Mean[dat];
  	s = StandardDeviation[dat];
  	min = (m - sdr s); max = (m + sdr s);
  	range = max - min;
  	step = range/100;
  	data = BinCounts[dat, {min, max, step}];
  	xdat = Range[min + 0.5 step, max - 0.5 step, step];
  	out = Transpose[{xdat, data/Length[dat]/step}];
  	DeleteCases[out, {_, 0.}]
   ]
  ];


(* ::Subsubsection::Closed:: *)
(*GetMaskMeans*)


Options[GetMaskMeans] = {MeanMethod -> "SkewNormalDist"}

SyntaxInformation[GetMaskMeans] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

GetMaskMeans[dat_, mask_, opts:OptionsPattern[]] := GetMaskMeans[dat, mask, "", opts]

GetMaskMeans[dat_, mask_, name_, OptionsPattern[]] := 
 Block[{labels, out, fl},
  labels = If[name==="",
  	{"mean", "std", "Median", "5%", "95%"},
  	name <> " " <> # & /@ {"mean", "std", "Median", "5%", "95%"}
  ];
  out = (
      fl = Flatten[GetMaskData[dat, #]];
      Switch[OptionValue[MeanMethod],
       "NormalDist",
       ParameterFit[fl, FitOutput -> "ParametersExtra", FitFunction -> "Normal"],
       "SkewNormalDist",
       ParameterFit[fl, FitOutput -> "ParametersExtra", FitFunction -> "SkewNormal"],
       _,
       Flatten[{Mean[fl], StandardDeviation[fl], Quantile[fl, {.5, .05, .95}]}]
       ]
      ) & /@ Transpose[mask];
  Prepend[out, labels]
  ]


(* ::Subsubsection::Closed:: *)
(*RegNorm*)


RegNorm[x_,Mu_,Sigma_]:=1/(E^((x - Mu)^2/(2*Sigma^2))*(Sqrt[2*Pi]*Sigma));


(* ::Subsubsection::Closed:: *)
(*SkewNorm*)


Phi[x_]:=1/(E^(x^2/2)*Sqrt[2*Pi]);
CapitalPhi[x_]:=.5(1+Erf[(x)/Sqrt[2]]);
SkewNorm[x_,Omega_,Xi_,Alpha_]:=(2/Omega)Phi[(x-Xi)/Omega]CapitalPhi[Alpha (x-Xi)/Omega];
Delta[a_]:=a/Sqrt[1+a^2];
Mn[w_,e_,a_]:=e+w Delta[a] Sqrt[2/Pi];
Var[w_,a_]:=w^2(1-(2Delta[a]^2/Pi));

SkewNormC=Compile[{{x, _Real},{Omega, _Real},{Xi, _Real},{Alpha, _Real}},
Chop[(2/Omega)(1/(E^(((x-Xi)/Omega)^2/2)*Sqrt[2*Pi]))(.5(1+Erf[((Alpha (x-Xi)/Omega))/Sqrt[2]]))]
];


(* ::Subsection::Closed:: *)
(*SliceData*)


SyntaxInformation[SliceData] = {"ArgumentsPattern" -> {_}};

SliceData[data_]:=
Module[{fitdat},
	fitdat=ParameterFit[DeleteCases[#,Null,{2}]]&/@data;
	Prepend[MapIndexed[Flatten[{#2,#1}]&,Transpose[fitdat]],
	{"Slice Number","first-Mean","first-SD","second-Mean","second-SD","third-Mean","third-SD","ADC-Mean","ADC-SD","FA-Mean","FA-SD"}]
	]


(* ::Subsection::Closed:: *)
(*SigmaCalc*)


Options[SigmaCalc] = {FilterShape -> "Median"};

SyntaxInformation[SigmaCalc] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

SigmaCalc[DTI_?ArrayQ, grad : {{_, _, _} ..}, bvalue_, blur_: 2, OptionsPattern[]] := Module[
	{tens, res,len,sig}, 
  tens = TensorCalc[DTI, grad, bvalue, MonitorCalc -> False];
  res = ResidualCalc[DTI, tens, grad, bvalue, MeanRes -> "MAD"];
  len = Length[grad];
  sig = Sqrt[len/(len - 7)]*res;
   PrintTemporary["Filtering noisemap"];
  Switch[OptionValue[FilterShape],
   "Gaussian",
   GaussianFilter[sig, blur],
   "Median",
   MedianFilter[sig, blur]
   ]
  ]


SigmaCalc[DTI_?ArrayQ, tens_?ArrayQ, grad : {{_, _, _} ..}, bvalue_, blur_: 2, OptionsPattern[]] := Module[
	{res,len,sig},
  res = ResidualCalc[DTI, tens, grad, bvalue, MeanRes -> "MAD"];
  len = Length[grad];
  sig = Sqrt[len/(len - 7)]*res;
   PrintTemporary["Filtering noisemap"];
  Switch[OptionValue[FilterShape],
   "Gaussian",
   GaussianFilter[sig, blur],
   "Median",
   MedianFilter[sig, blur]
   ]
  ]


(* ::Subsection::Closed:: *)
(*SNRCalc*)


SyntaxInformation[SNRCalc] = {"ArgumentsPattern" -> {_, _, _}};

SNRCalc[data_?ArrayQ,mask1_?ArrayQ,mask2_?ArrayQ]:=
Module[{noise,signal,Msignal,Mnoise},
	signal=GetMaskData[data,mask1];
	noise=GetMaskData[data,mask2];
	Mnoise=Mean[Flatten[noise]];
	Msignal=N[Map[Mean[#]&,signal]]/.Mean[{}]->0;
	Msignal/(0.8Mnoise)/.(1/Mean[{}])->0
	]



(* ::Subsection::Closed:: *)
(*SNRMapCalc*)


Options[SNRMapCalc] = {OutputSNR -> "SNR", SmoothSNR->2};

SyntaxInformation[SNRMapCalc] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};

SNRMapCalc[data_?ArrayQ, noise_?ArrayQ, opts:OptionsPattern[]] := SNRMapCalc[data, noise, 1, opts]
SNRMapCalc[data_?ArrayQ, noise_?ArrayQ, k_?NumberQ, OptionsPattern[]] := Module[{sigma, sigmac, snr, depthD, depthN},
	
 	sigma = N[GaussianFilter[noise, 4]];
 	sigmac = (sigma/Sqrt[Pi/2.]) /. 0. -> Infinity;
 	
 	depthD=ArrayDepth[data];
 	depthN=ArrayDepth[noise];
 	snr = If[k>=1,
 		If[depthD==depthN,
 		GaussianFilter[data/(sigmac), k],
 		If[depthD==depthN+1&&k>=1,
 			If[depthD==4,
 				Transpose[GaussianFilter[#/sigmac, k]&/@Transpose[data]],
 				GaussianFilter[#/sigmac, k]&/@data
 				] 				
 			]
 		]
 		,
 		If[depthD==depthN,
 		data/(sigmac),
 		If[depthD==depthN+1&&k>=1,
 			If[depthD==4,
 				Transpose[(#/sigmac)&/@Transpose[data]],
 				#/sigmac&/@data
 				] 				
 			]
 		]
 	];
  
  Switch[OptionValue[OutputSNR],
	 "Sigma", sigma,
	 "Both", {snr, sigma},
	 _, snr
	 ]
  ]

SNRMapCalc[{data1_?ArrayQ, data2_?ArrayQ}, opts:OptionsPattern[]] := SNRMapCalc[{data1, data2}, 2, opts]
SNRMapCalc[{data1_?ArrayQ, data2_?ArrayQ}, k_?NumberQ, OptionsPattern[]] := 
 Module[{noise, signal, sigma, snr},
  noise = (data1 - data2);
  signal = Mean[{data1, data2}];
  sigma = ConstantArray[StandardDeviation[DeleteCases[Flatten[noise] // N, 0.]],Dimensions[signal]];
  snr = GaussianFilter[signal/(.5 Sqrt[2] sigma), k];
  Switch[OptionValue[OutputSNR],
	 "Sigma", sigma,
	 "Both", {snr, sigma},
	 _, snr
	 ]
 ]

SNRMapCalc[data : {_?ArrayQ ...}, opts:OptionsPattern[]] := SNRMapCalc[data, 2, opts]
SNRMapCalc[data : {_?ArrayQ ...}, k_?NumberQ, OptionsPattern[]] := 
 Module[{signal, sigma, snr,div},
  signal = Mean[data];
  sigma = Chop[StandardDeviation[data]]-10^-15;
  div=N@Clip[signal / sigma, {0, Infinity}];
  div=Clip[div, {0., 100 Median[DeleteCases[Flatten[div], 0.]]}];
  snr = GaussianFilter[div, k];
  
  Switch[OptionValue[OutputSNR],
	 "Sigma", sigma,
	 "Both", {snr, sigma},
	 _, snr
	 ]
 ]


(* ::Subsection::Closed:: *)
(*MeanSignal*)


SyntaxInformation[CoilSNRCalc] = {"ArgumentsPattern" -> {_, _}};

(*calculate the combineds snr form multiple coils images*)
CoilSNRCalc[coils_, noise_] := 
 Block[{mn, sigmap, coilsN, noiseN, sumSquares, weights, snr},
  (*get mean noise*)
  mn = MeanNoZero@Flatten@N@noise;
  (*normalize all coils to constant noise level*)
  coilsN = 10. coils/mn;
  noiseN = 10. noise/mn;
  (*calcualte the sum of squares signal*)
  {sumSquares, weights} = SumOfSquares[coilsN];
  (*calculated the weitghted noise addition*)
  {snr, sigmap} = WeigthedSNR[coilsN, noiseN, weights];
  
  {coilsN, noiseN, weights, sumSquares, sigmap, snr}
  ]


WeigthedSNR[signal_, noise_, weights_] := Block[{sigmap, sigtot, snr},
  sigtot = Total[signal weights];
  sigmap = 
   Sqrt[Total[weights^2 (Sqrt[2./Pi] GaussianFilter[noise, 3])^2]];
  snr = DevideNoZero[sigtot, sigmap];
  {snr, sigmap}
  ]


(* ::Subsection::Closed:: *)
(*MeanSignal*)


Options[MeanSignal]={UseMask->True}

SyntaxInformation[MeanSignal] = {"ArgumentsPattern" -> {_, _.,OptionsPattern[]}};

MeanSignal[data_, opts : OptionsPattern[]] := MeanSignal[data, 1, opts];

MeanSignal[data_, pos_ ,OptionsPattern[]] := Block[{datat, mean, mask},
  datat = Transpose[data];
  
  If[ListQ[pos],
   mask = If[OptionValue[UseMask],
   	Round@GaussianFilter[Mask[datat[[First@pos]]/MeanNoZero[Flatten[datat[[First@pos]]]], .5], 5],
   	ConstantArray[1,Dimensions[data[[All,1]]]]];   
   mean = MeanNoZero[Flatten[#]] & /@ Transpose[MaskDTIdata[data[[All, pos]], mask]];
   ,
   mask = If[OptionValue[UseMask],
   	Round@GaussianFilter[Mask[datat[[pos]]/MeanNoZero[Flatten[datat[[pos]]]], .5], 5],
   	ConstantArray[1,Dimensions[data[[All,1]]]]];
   mean = MeanNoZero[Flatten[#]] & /@ Transpose[MaskDTIdata[data, mask]];
   ];
  N@mean
  ]


(* ::Subsection::Closed:: *)
(*FiberDensityMap*)


Options[FiberDensityMap] = {SeedDensity -> Automatic};

SyntaxInformation[FiberDensityMap] = {"ArgumentsPattern" -> {_, _, _,OptionsPattern[]}};

FiberDensityMap[fibers_, dim_, vox_, OptionsPattern[]] := 
 Module[{pixindex, density, dens,densi},
  pixindex = GetFiberCoor[fibers, vox];
  pixindex = Transpose[MapThread[Clip[#1, {1, #2}] &, {Transpose[pixindex], dim}]];
  density = CountVoxels[ConstantArray[0, dim], pixindex];
  densi = OptionValue[SeedDensity];
  (*Print[{(Times @@ vox)/0.75,Median[DeleteCases[Flatten[density], 0]]}];*)
  dens = If[NumberQ[densi],
    Times @@ (vox/densi),
    Median[DeleteCases[Flatten[density], 0]]
    ];
  Clip[NormalizeDens[density, dens], {0., 10.}]
  ]

GetFiberCoor = Compile[{{fibcor, _Real, 1}, {vox, _Real, 1}},
   Round[Reverse[fibcor + vox]/vox],
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", 
   Parallelization -> True];

CountVoxels = Compile[{{const, _Integer, 3}, {pix, _Integer, 2}}, Block[{out = const},
    (out[[#[[1]], #[[2]], #[[3]]]] += 1) & /@ pix;
    out
    ]];

NormalizeDens = Compile[{{dens, _Integer, 3}, {n, _Real, 0}}, dens/n];


(* ::Subsection::Closed:: *)
(*FiberLengths*)


SyntaxInformation[FiberLengths] = {"ArgumentsPattern" -> {_, _.}};

FiberLengths[fpoints_, flines_] := FiberLengths[{fpoints, flines}]
FiberLengths[{fpoints_, flines_}] := Module[{len, mpos},
   len = (Length /@ flines) - 1;
   mpos = First@First@Position[len, Max[len]];
   len Mean[EuclideanDistance @@@ Partition[fpoints[[flines[[mpos]]]], 2, 1]]
];


(* ::Subsection::Closed:: *)
(*DataTot and DataTotXLS*)


SyntaxInformation[DatTot] = {"ArgumentsPattern" -> {_, _, _}};

DatTot[data_,name_,vox_]:=
Module[{fitdat},
	fitdat=ParameterFit[DeleteCases[Flatten[#],Null]&/@data];
	With[{Quant=Function[dat,{dat[[1]],dat[[2]],100dat[[2]]/dat[[1]]}]},
		Flatten[{name,vox[[1]],vox[[2]],Quant[fitdat[[1]]],Quant[fitdat[[2]]],Quant[fitdat[[3]]],Quant[fitdat[[4]]],Quant[fitdat[[5]]]}]
		]
	]

SyntaxInformation[DatTotXLS] = {"ArgumentsPattern" -> {_, _, _}};

DatTotXLS[data_,name_,vox_]:=
Module[{fitdat},
	fitdat=ParameterFit[DeleteCases[Flatten[#],Null]&/@data];
	With[{Quant=Function[dat,ToString[Round[dat[[1]],.01]]<>" \[PlusMinus] "<>ToString[Round[dat[[2]],.01]]]},
		Flatten[{name,vox[[1]],vox[[2]],Quant[fitdat[[1]]],Quant[fitdat[[2]]],Quant[fitdat[[3]]],Quant[fitdat[[4]]],Quant[fitdat[[5]]]}]
		]
	]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
