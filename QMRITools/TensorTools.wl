(* ::Package:: *)

(* ::Title:: *)
(*QMRITools TensorTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`TensorTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`TensorTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DriftCorrect::usage = 
"DriftCorrect[data, bval] dirft corrects the data using the signals of the lowest bvalue that has 6 or more unique volumes.
For the function to work optimal it is best to have these volumes evenly spread througout thet data \
and for the first and last volume to have this low bvalue.

DriftCorrect[] is based on DOI: 10.1002/mrm.26124." 

ConcatenateDiffusionData::usage=
"ConcatenateDiffusionData[{{data1, .., dataN}, {grad1, .., gradN}, {bval, .., bvalN}, {vox, .., voxN}}] concatenates the diffusion data sets.
ConcatenateDiffusionData[{data1, .., dataN}, {grad1, .., gradN}, {bval, .., bvalN}, {vox, .., voxN}] concatenates the diffusion data sets."

SortDiffusionData::usage = 
"SortDiffusionData[data, grad, bval] sorts the diffusion datasets grad and bval for magnitude of bvalue."

RemoveIsoImages::usage = 
"RemoveIsoImages[data, grad, bval] Romoves the ISO images from the philips scanner from the data. ISO images have g={0,0,0} and b>0."


TensorCalc::usage = 
"TensorCalc[data, gradients, bvalue] calculates the diffusion tensor for the given dataset. Allows for one unweighted image and one b value. 
Gradient directions must be in the form {{x1,y1,z1}, ..., {xn,yn,zn}} without the unweighted gradient direction. 
bvalue is a singe number indicating the b-value used.
TensorCalc[data, gradients, bvec] calculates the diffusion tensor for the given dataset. allows for multiple unweighted images and multiple bvalues. 
allows for differnt tensor fitting methods. gradient directions must be in the form {{x1,y1,z1}, ..., {xn,yn,zn}} with the unweighted direction as {0,0,0}. \
bvec the bvector, with a bvalue defined for each gradient direction. b value for unweighted images is 0.
TensorCalc[data, bmatix] calculates the diffusion tensor for the given dataset. allows for multiple unweighted images and multiple bvalues. 
bmat is the bmatrix which can be generated usiong Bmatrix.

The bvalue assumed to be is in s/mm^2 and therfore the output is in mm^2/2

TensorCalc[] is based on DOI: 10.1016/j.neuroimage.2013.05.028 and 10.1002/mrm.25165."


ResidualCalc::usage =
"ResidualCalc[DTI,{tensor,S0},gradients,bvector] calculates the tensor residuals for the given dataset.
ResidualCalc[DTI,{tensor,S0},outlier,gradients,bvector] calculates the tensor residuals for the given dataset taking in account the outliers.
ResidualCalc[DTI,{tensor,S0},bmat] calculates the tensor residuals for the given dataset.
ResidualCalc[DTI,{tensor,S0},outlier,bmat] calculates the tensor residuals for the given dataset taking in account the outliers.
ResidualCalc[DTI,tensor,gradients,bvector] calculates the tensor residuals for the given dataset. Tensor must contain Log[S0].
ResidualCalc[DTI,tensor,outlier,gradients,bvector] calculates the tensor residuals for the given dataset taking in account the outliers. Tensor must contain Log[S0].
ResidualCalc[DTI,tensor,bmat] calculates the tensor residuals for the given dataset. Tensor must contain Log[S0].
ResidualCalc[DTI,tensor,outlier,bmat] calculates the tensor residuals for the given dataset taking in account the outliers. Tensor must contain Log[S0]."

SigmaCalc::usage = 
"SigmaCalc[DTI,grad,bvec] calculates the noise sigma based on the tensor residual, using a blur factor of 10.
SigmaCalc[DTI,tens,grad,bvec] calculates the noise sigma based on the tensor residual, using a blur factor of 10.
SigmaCalc[DTI,grad,bvec,blur] calculates the noise sigma based on the tensor residual, If blur is 1 ther is no blurring.
SigmaCalc[DTI,tens,grad,bvec,blur] calculates the noise sigma based on the tensor residual. If blur is 1 ther is no blurring."


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

ParameterCalc::usage = "ParameterCalc[tensor] caculates the eigenvalues and MD and FA from the given tensor. The parameters are l1, l2, l3, MD and FA. l1, l2, l3, MD are in (10^-3 mm^2/s)"


AngleCalc::usage = 
"AngleCalc[data, vector] calculates the angel between the vector and the data. Data shoud be an array of dimensions {xxx,3}."

AngleMap::usage = 
"AngleMap[data] calculates the zennith and azimuth angles of a 3D dataset (z,x,y,3) containing vectors relative to the slice direction."


Correct::usage =
"Correct[data, phase, shiftpar] corrects the dataset data using the phasemap and the shiftpar and interpolation order 1.
Correct[data, phase, shiftpar, int] corrects the dataset data using the phasemap and the shiftpar and interpolation order int."

TensorCorrect::usage=
"TensorCorrect[tensor, phase, shift, vox] corrects the tensor based on B0 field map. Can perform both translation and rotation of tensor."

Deriv::usage = 
"Deriv[disp, vox] calculates the derivative of the displacement along the three main axes. disp is the displacement field, vox is the voxel size.
Deriv[disp, vox, mask] calculates the derivative of the displacement along the three main axes. Sharp edges between the background en disp are solved by the mask. mask is a mask delining the edge of the displacement field.";


ColorFAPlot::usage = 
"ColorFAPlot[tenor] create a color coded FA map from the tensor for l1, l2 and l3."


(* ::Subsection::Closed:: *)
(*Options*)


NormalizeSignal::usage = 
"NormalizeSignal is an option for DriftCorrect."


MonitorCalc::usage = 
"MonitorCalc is an option for all Calc fucntions. When true the proceses of the calculation is shown."

FullOutput::usage = 
"FullOutput is an option for TensorCalc when using bvector. When True also the S0 is given as output."

RobustFit::usage = 
"RobustFit is an option for TensorCalc. If true outliers will be rejected in the fit, only works with WLLS.
If FullOutput is given the outlier map is given.";

RobustFitParameters::usage =
"RobustFitParameters is an option for TensorCalc. gives the threshold for stopping the itterations and the kappa for the outlier marging, {tr,kappa}."


FilterShape::usage = 
"FilterShape is an option for SigmaCalc. Can be \"Gaussian\" of \"Median\"."

RejectMap::usage = 
"RejectMap is an option for EigenvalCalc. If Reject is True and RejectMap is True both the eigenvalues aswel as a map showing je rejected values is returned."

Reject::usage = 
"Reject is an option for EigenvalCalc. It True then voxels with negative eigenvalues are rejected and set to 0."

Distribution::usage = 
"Distribution is an option for AngleCalc. values can be \"0-180\", \"0-90\" and \"-90-90\"."

MeanRes::usage = 
"MeanRes is an option for ResidualCalc. When True the root mean square of the residual is calculated."


RotationCorrect::usage =
"RotationCorrect is an option for TensorCorrect. Default is False. Is a tensor is deformed setting to True also the shear is accounted for by local rotation of the tensor"


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

ResidualCalc::datdim = "DTIdata (`1`) and tensor data (`2`) are not the same dimensions."

AngleCalc::dist = "Unknown option (`1`), options can be. \"0-180\", \"0-90\" or \"-90-90\"."

ConcatenateDiffusionData::dim= "data, grad and bval should be the same length:  data `1` / grad `2` / bval `2`"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*DriftCorrect*)


Options[DriftCorrect]={NormalizeSignal->True, UseMask->True}

SyntaxInformation[DriftCorrect] = {"ArgumentsPattern" -> {_, _,_., OptionsPattern[]}}

DriftCorrect[data_, bi_, opts:OptionsPattern[]] := 
 Block[{bval,pos},

  bval = If[ArrayDepth[bi] == 2, BmatrixInv[bi][[1]], bi];
  pos = First@UniqueBvalPosition[bval, 5][[2]];

	DriftCorrect[data, bval, pos, opts]

  ];
  
DriftCorrect[data_, bi_, pos_, OptionsPattern[]] := 
 Block[{sig, cor, bval, sol1, sol2, sol3, a, b, c, x,outp},
  bval = If[ArrayDepth[bi] == 2, BmatrixInv[bi][[1]], bi];
  sig = MeanSignal[data, pos, UseMask->OptionValue[UseMask]];
  
  dat = Transpose[{pos, sig}];
  
  {sol1, sol2, sol3} = {a, b, c} /. FindFit[dat, {c + b x + a x^2}, {a, b, c}, x];
  cor = sol3/Table[sol3 + sol2 x + sol1 x^2, {x, 1, Length[bi]}];
  
  outp = ConstantArray[cor, Length[data]] data;
  
  If[OptionValue[NormalizeSignal], 100 outp / (sig[[1]] cor[[1]]) , outp ]
  ];


(* ::Subsection::Closed:: *)
(*ConcatenateDiffusionData*)


SyntaxInformation[ConcatenateDiffusionData] = {"ArgumentsPattern" -> {_, _., _., _.}};

ConcatenateDiffusionData[data_?ListQ] :=If[Length[data] == 4,ConcatenateDiffusionData[data[[1]], data[[2]], data[[3]], data[[4]]]]

ConcatenateDiffusionData[data_, grad_, val_, vox_] := 
  Module[{dataout, gradout, valout, voxout},
   If[Length[data] == Length[grad] == Length[val],
    dataout = Transpose@Flatten[Transpose[NormalizeData[#]] & /@ data, 1];
    gradout = Flatten[grad, 1];
    valout = Flatten[val];
    
    {dataout, gradout, valout} = RemoveIsoImages[dataout, gradout, valout];
    {dataout, gradout, valout} = SortDiffusionData[dataout, gradout, valout];
    ,
    Return[Message[ConcatenateDiffusionData::dim, Length[data],Length[grad], Length[val]]];
    ];
   
   voxout = If[ListQ[vox] && ! ListQ[vox[[1]]], vox, vox[[1]]];
   
   {dataout, gradout, valout, voxout}
   ];


(* ::Subsection::Closed:: *)
(*SortDiffusionData*)


SyntaxInformation[SortDiffusionData] = {"ArgumentsPattern" -> {_, _, _}};

SortDiffusionData[data_, grad_, val_] := Module[{pos, valu, sel},
  {valu, pos} = UniqueBvalPosition[val];
  sel = Flatten[pos];
  {data[[All, sel]], grad[[sel]], val[[sel]]}
  ]


(* ::Subsection::Closed:: *)
(*RemoveIsoImages*)


SyntaxInformation[RemoveIsoImages] = {"ArgumentsPattern" -> {_, _, _}};

RemoveIsoImages[data_, grad_, val_] := Module[{sel},
  sel = Complement[Range[Length[val]], 
    Complement[Flatten[Position[grad, {0., 0., 0.}]], 
     Flatten[Position[val, 0.]]]];
  {data[[All, sel]], grad[[sel]], val[[sel]]}
  ]


(* ::Subsection:: *)
(*TensorCalc*)


(* ::Subsubsection::Closed:: *)
(*TensorCalc*)


Options[TensorCalc]= {MonitorCalc->True, Method->"iWLLS", FullOutput->True, RobustFit->True, Parallelize->True , RobustFitParameters->{10.^-4,6}};

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
	dataL=Chop[LogNoZero[Chop[data]]];	
	(*calculate the inverse bmat*)
	bmatI=PseudoInverse[bmat];
	
	(*if data is 4D handle as multiple 3D sets (saves memory and calculation time)*)
	If[depthD==4,
		
		xx=0;
		
		func=If[OptionValue[Parallelize],
			SetSharedVariable[xx];
			DistributeDefinitions[bmat, bmatI, method, output ,robust, con ,kappa, 
				TensorCalci, FindTensOutliers, TensMinLLS, TensMinWLLS, TensMiniWLLS, ResidualCalc,
				ExpNoZero, RMSNoZero, MADNoZero, Bmatrix];
			ParallelMap,
			Map
		];	
		
		If[OptionValue[MonitorCalc],PrintTemporary[ProgressIndicator[Dynamic[xx], {0, Length[data]}]]];
		result = func[(
			xx++;
			TensorCalci[#[[1]], #[[2]], bmat, bmatI, Method->method, FullOutput->output, RobustFit->robust, RobustFitParameters->{con,kappa}]
			)&,Transpose[{data, dataL}]];
			
		result = Transpose[result];
		
		(*full output returns {tens,S0,(outliers),residuals}*)
		If[output, result[[1]] = Transpose[result[[1]]]];	
		
		,(*1D,2D,3D*)
		result = TensorCalci[data,dataL,bmat,bmatI, Method->method, FullOutput->output, RobustFit->robust, RobustFitParameters->{con,kappa}];
	];

	result
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
	
	
	outliers = If[robust,
		Transpose[FindTensOutliers[Transpose[dataL,l], bmat, con, kappa], r]
		, 
		ConstantArray[0.,Dimensions[data]]
		];
	
	fitresult = Switch[method,
		"LLS", Transpose[TensMinLLS[Transpose[dataL,l], bmatI], r],
		"WLLS", Transpose[TensMinWLLS[Transpose[(1-outliers) data,l],Transpose[(1-outliers) dataL,l], bmat], r],
		"iWLLS", Transpose[TensMiniWLLS[Transpose[(1-outliers) data,l],Transpose[(1-outliers) dataL,l], bmat], r]
		];
	
	If[OptionValue[FullOutput], residual = ResidualCalc[data,fitresult,outliers,bmat, MeanRes->"MAD"]];
		
	S0 = N@Clip[ExpNoZero[N@Chop[Last[fitresult]]],{0, 1.5 Max[data]}];
	tensor = N@Clip[Drop[fitresult,-1],{-0.1,0.1}];

	If[OptionValue[FullOutput],If[robust,{tensor,S0,outliers,residual},{tensor,S0,residual}],tensor]
]


(* ::Subsubsection::Closed:: *)
(*FindOutliers*)


FindTensOutliers = Quiet@Compile[{{LS, _Real, 1}, {bmat, _Real, 2}, {con, _Real, 0}, {kappa, _Real, 0}},	Block[
		{ittA, contA, solA, itt, cont, soli, res, mad, wts, wmat, fitE, LS2, bmat2, out},
		(*initialize some values*)
		out = 0. LS; LS2 = LS; bmat2 = bmat;
		
		(*skip if background*)
	  	If[Total[Unitize[LS]] >= 7 ,
	  		(*If not background find the outliers*)
	  		(*moniotr the overall while loop*)
	  		ittA = 0; contA = 1;
	  			
	  		(*Step1: initial LLS fit*)
	  		sol = LeastSquares[bmat,LS];
	  			
	  		(*check if LLS fit is plausable, i.e. S0 > 0*)
	  		If[Last[sol]>0,
	  			While[contA == 1,ittA++;
	  			(*init the solution*)
	  			solA = sol;
	  				
	  			(*Step2: Compute a robust estimate for homoscedastic regression using IRLS.*)
	  			itt = 0; cont = 1;
	  			While[cont == 1, itt++;
	  				soli = sol;
	  				(*a. Calculate the residuals e* in the linear domain*)
	  				res = LS - bmat.sol;
	  				(*b. Obtain an estimate of the dispersion of the residuals by calculating the median absolute deviation (MAD).*)
	  				mad = N@Chop[1.4826 MedianDeviation[res]];
	  				(*prevent calculation with 0*)
	  				If[AllTrue[res, (0. === #) &] || mad === 0. ,
	  					cont = 0,
	  					(*c. Recompute the weights according to Eq. [13].*)
	  					wts = 1 / (1 + (res/mad)^2)^2;
	  					(*d. Perform WLLS fit with new weights*)
	  					wmat = Transpose[bmat].DiagonalMatrix[wts];
	  					sol = Chop[LeastSquares[wmat.bmat, wmat.LS]];
	  					(*e. Check convergence*)
	  					If[!AnyTrue[Abs[sol - soli] - con (Max /@ Transpose[{Abs[sol], Abs[soli]}]),Positive] || itt === 5, cont = 0];
	  				];
	  			];(*end first while*)
	
				(*Step 3: Transform variables for heteroscedasticity*)
				fitE = Exp[bmat.sol] + 10^-10;
				LS2 = LS / fitE;
				bmat2 = bmat / fitE;
				
				(*Step 4: Initial LLS fit in * domain*)
				sol = LeastSquares[bmat2,LS2];
	   
				(*Step 5: Compute a robust estimate for homoscedastic regression using IRLS.*)
				itt = 0; cont = 1;
				While[cont == 1, itt++;
					soli = sol;
					(*a. Calculate the residuals e* in the linear domain*)
					res = LS2 - bmat2.sol;
					(*b. Obtain an estimate of the dispersion of the residuals by calculating the median absolute deviation (MAD).*)
					mad = N@Chop[1.4826 MedianDeviation[res]];
					(*prevent calculation with 0*)
					If[AllTrue[res, (0. === #) &] || mad === 0.,
						cont = 0,
						(*c. Recompute the weights according to Eq. [13].*)
						wts = 1 / (1 + (res/mad)^2)^2;
						(*d. Perform WLLS fit with new weights*)
						wmat = Transpose[bmat2].DiagonalMatrix[wts];
						sol = Chop[LeastSquares[wmat.bmat2, wmat.LS2]];
						(*e. Check convergence*)
						If[!AnyTrue[Abs[sol - soli] - con (Max /@ Transpose[{Abs[sol], Abs[soli]}]),Positive] || itt === 5, cont = 0];
					];
				];(*end second while*)
			
				(*Step 6: Check convergence overall loop*)
				If[!AnyTrue[Abs[sol - solA] - con (Max /@ Transpose[{Abs[sol], Abs[solA]}]), Positive] || ittA === 10, contA = 0];
			];(*end main while*)
	  			
			(*Step 7: Identify and exclude outliers*)
			res = LS2 - bmat2.sol;
			out = UnitStep[Abs[res] - (kappa 1.4826 MedianDeviation[res])];
			
			];(*close if negative S0*)
		];(*close if background*)
		
		out
	],
	
	{{sol, _Real, 1}},
	RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages"->False}, 
	CompilationOptions -> {"ExpressionOptimization" -> False}
]



(* ::Subsubsection::Closed:: *)
(*LLS*)


TensMinLLS = Compile[{{LS, _Real, 1}, {bmatI, _Real, 2}},
	bmatI.LS,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*WLLS*)


TensMinWLLS = Compile[{{S, _Real, 1},{LS, _Real, 1},{bmat, _Real, 2}}, 
	Block[{wmat,mvec,sol},
		sol = 0. First[bmat];
		If[!(AllTrue[LS, 0. === # &] || Total[Unitize[LS]] < 7),
	    	mvec = UnitStep[LS] Unitize[LS]; 
	    	(*if 0 then it is not used because w=0*)
	    	wmat = DiagonalMatrix[mvec S^2];
	    	sol = LeastSquares[Transpose[bmat].wmat.bmat,Transpose[bmat].wmat.LS]
	    ];
	    sol]
    ,{{wmat,_Real,2}, {sol, _Real, 1}}, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*iWLLS*)


TensMiniWLLS = Quiet@Compile[{{S, _Real, 1}, {LS, _Real, 1}, {bmat, _Real, 2}},
	Block[{wmat, mat, cont, itt, mvec, soli, sol0, max, sol, w},
    mvec = UnitStep[S] Unitize[S];
    max = Max[mvec S];
    sol0 = 0. First[bmat];
    sol = sol0;
    (*skip background or not enough data for fit*)
    If[!(AllTrue[LS, 0. === # &] || Total[mvec] <= 7),
     (*initialize*)
     itt = 0;
     cont = 1;
     (*initialize using LLS*)
     sol = N@Chop@LeastSquares[bmat, LS];
     (*check for implausabole solution (negative S0 or high S0)*)
     If[Last[sol] >= 3*max || Last[sol] <= 0.,
      sol = sol0;
      ,
      (*itterative reweighting*)
      While[cont == 1,
       (*init itteration values*)
       itt++;
       soli = sol;
       (*perform WLLS*)
       w = (mvec Exp[bmat.sol])^2;
       wmat =Transpose[bmat].DiagonalMatrix[w];
       sol = LeastSquares[wmat.bmat, wmat.LS];
       (*update weight*)
       (*see if to quit loop*)
       If[(Last[sol] >= 3*max || Last[sol] <= 0), 
       	cont = 0.; sol = sol0
       	];
       If[! AnyTrue[Abs[sol - soli] - 0.0001 (Max /@ Transpose[{Abs[sol], Abs[soli]}]), Positive] || itt === 10 , cont = 0.];
       ]]
       ];
    sol], 
    {{mat, _Real, 2}, {wmat, _Real, 2}, {sol, _Real, 1}},
    RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages"->False}];
        

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


(* ::Subsection:: *)
(*Tensor Parameters*)


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


(* ::Subsubsection::Closed:: *)
(*ADCCalc*)


SyntaxInformation[ADCCalc] = {"ArgumentsPattern" -> {_}};

ADCCalc[eig_] := ADCCalci[eig]

ADCCalci = Compile[{{eig, _Real, 1}}, Mean[eig], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*FACalc*)


SyntaxInformation[FACalc] = {"ArgumentsPattern" -> {_}};

FACalc[eig_] := FACalci[eig]

FACalci = Compile[{{eig, _Real, 1}}, Block[{l1, l2, l3, teig},
   l1 = eig[[1]]; l2 = eig[[2]]; l3 = eig[[3]];
   teig = Sqrt[2.*Total[eig^2]];
   If[teig == 0., 0. ,Sqrt[(l1 - l2)^2 + (l2 - l3)^2 + (l1 - l3)^2]/teig]
   ], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


(* ::Subsection:: *)
(*Correct*)


(* ::Subsubsection::Closed:: *)
(*Correct*)


SyntaxInformation[Correct] = {"ArgumentsPattern" -> {_, _, _, _.}};

Correct[data_?MatrixQ,phase_?MatrixQ,shift_]:=
Correcti[data,phase,shift,1]

Correct[data_?MatrixQ,phase_?MatrixQ,shift_,int_]:=
Correcti[data,phase,shift,int]

Correct[data:{_?MatrixQ..},phase:{_?MatrixQ..},shift_]:=
MapThread[Correcti[#1,#2,shift,1]&,{data,phase}]

Correct[data:{_?MatrixQ..},phase:{_?MatrixQ..},shift_,int_]:=
MapThread[Correcti[#1,#2,shift,int]&,{data,phase}]

Correct[data:{{_?MatrixQ..}..},phase:{_?MatrixQ..},shift_]:=
Transpose[Map[MapThread[Correcti[#1,#2,shift,1]&,{#,phase}]&,Transpose[data,{2,1}]],{2,1}]

Correct[data:{{_?MatrixQ..}..},phase:{_?MatrixQ..},shift_,int_]:=
Transpose[Map[MapThread[Correcti[#1,#2,shift,int]&,{#,phase}]&,Transpose[data,{2,1}]],{2,1}]


(* ::Subsubsection::Closed:: *)
(*Correcti*)


Correcti[dat_,ph_,shift_,int_]:=
Module[{pos,acpos,shiftpx,data,phase,output},
	If[shift[[2]]=="COL",
		data=Transpose[dat];phase=Transpose[ph];,
		data=dat;phase=ph;
		];
	shiftpx=phase*shift[[1]];
	output=Round[MapThread[(
		pos=Range[Length[#1]];
		acpos=pos-#1;
		ListInterpolation[#2,InterpolationOrder->int][acpos]
		)&,{shiftpx,data}]];
	If[shift[[2]]=="COL",Return[Transpose[output]],Return[output]]
	]


(* ::Subsection:: *)
(*TensorCorrect*)


(* ::Subsubsection::Closed:: *)
(*TensorCorrect*)


Options[TensorCorrect]={RotationCorrect->False};

SyntaxInformation[TensorCorrect] = {"ArgumentsPattern" -> {_, _, _, _, _., OptionsPattern[]}};

(* zonder masker, dus met sprongen in de afgeleide by grens tussen deformatie veld en achtergrond *)
TensorCorrect[tens_,phase_,shift_,vox_,OptionsPattern[]]:=
TensorCorrect[tens,phase,0,shift,vox];

(* met masker, dus zonder sprongen in de afgeleide by grens tussen deformatie veld en achtergrond *)
TensorCorrect[tens_,phase_,mask_,shift_,vox_,OptionsPattern[]]:=
	Module[{dim,pxshift,der,F,tensM,tensC,tensCV,tensT},
	
	dim=Dimensions[phase];
	(*deformation expessed in pixels*)
	pxshift=phase*shift[[1]];
	
	If[OptionValue[RotationCorrect]==True,
		PrintTemporary["Cacluating Derivative"];
		(*local derivative of the displacement in the slice direction*)
		der=If[!ArrayQ[mask],
		Deriv[pxshift,vox],
		Deriv[pxshift,vox,mask]
		];
		F=Fmat[der,shift[[2]]];
		
		PrintTemporary["Rotation Correction"];
		(*rotation correction of matrix*)
		(*tensor to matrixform*)
		tensM=TensMat[tens];
		(*rotation correct tensor matrix*)
		tensC=MapThread[DRot[#1,#2]&,{tensM,F},3];
		(*corrected tensor back to vector form*)
		tensCV=TensVec[tensC];
		,
		tensCV=tens;
		];
	
	PrintTemporary["Translation Correction"];
	(*Translation correction of the rotation corrected Tensor*)
	tensT=Map[(
		MapThread[
			TransCorrect[#1,#2,shift[[2]],1]
			&,{#,pxshift}
			]
		)&,tensCV]
	];


(* ::Subsubsection::Closed:: *)
(*TransCorrect*)


(* Translation correct one slice*)
TransCorrect[dat_,sh_,dir_,int_]:=
Module[{data,shift,pos,acpos,out},
	(*Transpose the data zo the deformation is always in the "ROW" direction*)
	If[dir=="COL",
		data=Transpose[dat];shift=Transpose[sh];,
		data=dat;shift=sh;
		];
	(*{dims,dimx,dimy}=Dimensions[data];*)
	(*deformation Correction*)
	out=MapThread[(
		pos=Range[Length[#1]];
		acpos=pos-#1;
		ListInterpolation[#2,InterpolationOrder->int][Clip[acpos,{1,Length[acpos]}]]
		(*[{Clip[acpos[[1]],{1,dimx}],Clip[acpos[[2]],{1,dimy}]}]*)
		)&,{shift,data}];
	
	(*If deformation was in the "COL" direction rotate back*)
	If[dir=="COL",Return[Transpose[Chop[out]]],Return[Chop[out]]]
	];


(* ::Subsubsection::Closed:: *)
(*FMat*)


Fmat[der_,shift_]:=
Module[{Dx,Dy,Dz,dim,zero,ones,F},
	{Dx,Dy,Dz}=der;
	dim=Dimensions[Dx];
	zero=ConstantArray[0,dim];
	ones=ConstantArray[1,dim];
	If[shift=="COL",
		F=Transpose[{{ones,zero,zero},{Dx,Dy+1,Dz},{zero,zero,ones}},{4,5,1,2,3}],
		If[shift=="ROW",
			F=Transpose[{{Dx+1,Dy,Dz},{zero,ones,zero},{zero,zero,ones}},{4,5,1,2,3}]
			],
		Print["error, unknown direction"]
		]
	];


(* ::Subsubsection::Closed:: *)
(*Drot*)


DRot[D_,F_]:=Module[{val,e1,e2,e3,n1,n2,n3,NN},
{val,{e1,e2,e3}}=Eigensystem[D];
n1=Normalize[F.e1];
n2=Normalize[F.e2-(n1.(F.e2))*n1]//N;
n3=Normalize[Cross[n1,n2]]//N;
NN=Transpose[{n1,n2,n3}];
Chop[NN.(IdentityMatrix[3]val).Transpose[NN]]
];


(* ::Subsection:: *)
(*Deriv*)


(* ::Subsubsection::Closed:: *)
(*Deriv*)


SyntaxInformation[Deriv] = {"ArgumentsPattern" -> {_, _, _}};

Deriv[disp_,vox_]:=
Module[{dim,Dx,Dy,Dz},
	dim=Dimensions[disp];
	Dx=Transpose[DerivFunc[Transpose[disp,{1,3,2}],dim[[2]],vox[[2]]],{1,3,2}];
	Dy=DerivFunc[disp,dim[[3]],vox[[3]]];
	Dz=Transpose[DerivFunc[Transpose[disp,{3,2,1}],dim[[1]],vox[[1]]],{3,2,1}];
	{Dx,Dy,Dz}
	];

Deriv[disp_,vox_,mask_]:=
Module[{dim,Dx,Dy,Dz},
	dim=Dimensions[disp];
	Dx=Transpose[DerivFunc[Transpose[disp,{1,3,2}],Transpose[mask,{1,3,2}],dim[[2]],vox[[2]]],{1,3,2}];
	Dy=DerivFunc[disp,mask,dim[[3]],vox[[3]]];
	Dz=Transpose[DerivFunc[Transpose[disp,{3,2,1}],Transpose[mask,{3,2,1}],dim[[1]],vox[[1]]],{3,2,1}];
	{Dx,Dy,Dz}
	];


(* ::Subsubsection::Closed:: *)
(*DerivFunc*)


DerivFunc[disp_,length_,step_]:=
Module[{coor,f},
	coor=Range[length]*step;
	Map[(
		f=Interpolation[Transpose[{coor,#}],InterpolationOrder->1];
		Table[f'[st],{st,coor}]
		)&,disp,{2}]
	];

DerivFunc[disp_,mask_,length_,step_]:=
Module[{coor,f,fr,df,dfr},
	coor=Range[length]*step;
	MapThread[(
		f=Interpolation[Transpose[{coor,#1}],InterpolationOrder->1];
		fr=Interpolation[Transpose[{coor,Reverse[#1]}],InterpolationOrder->1];
		df=RotateRight[#2,1]*#2*Table[f'[st],{st,coor}];
		dfr=RotateLeft[#2,1]*#2*-Reverse[Table[fr'[st],{st,coor}]];
		MapThread[If[#1==0,#2,#1]&,{df,dfr}]
		)&,{disp,mask},2]
	];


(* ::Subsection:: *)
(*Angle maps*)


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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
(*ColorFAPlot*)


SyntaxInformation[ColorFAPlot] = {"ArgumentsPattern" -> {_}};  
  
ColorFAPlot[tens_] := Block[{FA, eig, eigv, mid, eigFA, mask},
  {eig, eigv} = Transpose[EigensysCalc[tens], {2, 3, 4, 1, 5}];
  mask = Mask[tens[[1]], 10^-6];
  
  eigv = mask Abs[eigv];
  FA = FACalc[eig];
  eigFA = mask FA eigv;
  
  DynamicModule[{colEigFA, colEig, im},
   colEigFA = Table[Image[eigFA[[j, All, All, i]], ColorSpace -> "RGB"], {j, 1, Length[eigv]}, {i, 1, 3}];
   colEig = Table[Image[eigv[[j, All, All, i]], ColorSpace -> "RGB"], {j, 1, Length[eigv]}, {i, 1, 3}];
   
   Manipulate[
    im = GraphicsRow[{colEig, colEigFA}[[i, j, sel]], ImageSize -> Length[sel]*size],
    {{j, Round[Length[colEig]/2], "slice"}, 1, Length[colEig], 1},
    {{i, 2, "method"}, {1 -> "raw", 2 -> "FA"}},
    {{sel, {1, 2, 3}, "eigenvectors"}, {1 -> "first", 2 -> "second", 3 -> "third"}, ControlType -> TogglerBar},
    Button["save image", SaveImage[im], Method -> "Queued"],
    {{size, 300, "image size"}, {200, 300, 400, 600, 1000}}]
   ]
  ]  


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
