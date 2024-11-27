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


(* ::Subsection:: *)
(*Functions*)


TensorCalc::usage = 
"TensorCalc[data, gradients, bvalue] calculates the diffusion tensor for the given dataset. Allows for one unweighted image and one b value. 
Gradient directions must be in the form {{x1,y1,z1}, ..., {xn,yn,zn}} without the unweighted gradient direction.
bvalue is a singe number indicating the b-value used.
TensorCalc[data, gradients, bvec] calculates the diffusion tensor for the given dataset. allows for multiple unweighted images and multiple bvalues.
allows for differnt tensor fitting methods. gradient directions must be in the form {{x1,y1,z1}, ..., {xn,yn,zn}} with the unweighted direction as {0,0,0}.
bvec the bvector, with a bvalue defined for each gradient direction. b value for unweighted images is 0.
TensorCalc[data, bmatix] calculates the diffusion tensor for the given dataset. allows for multiple unweighted images and multiple bvalues.
bmat is the bmatrix which can be generated usiong Bmatrix.

The bvalue assumed to be is in s/mm^2 and therfore the output is in mm^2/2.

TensorCalc[] is based on DOI: 10.1016/j.neuroimage.2013.05.028 and 10.1002/mrm.25165."


FlipTensorOrientation::usage = 
"FlipTensorOrientation[tens, perm] permutes the internal orientation of the tensor, perm can be any permutation of {\"x\",\"y\",\"z\"}.
FlipTensorOrientation[tens, flip] flips the internal orientation of the tensor, flip can be {1,1,1}, {-1,1,1}, {1,-1,1} or {1,1,-1}.
FlipTensorOrientation[tens, flip, perm] flips and permuter the internal orientation of the tensor.
FlipTensorOrientation[tens, perm, flip]flips and permuter the internal orientation of the tensor."

FlipGradientOrientation::usage = 
"FlipGradientOrientation[grad, perm] permutes the internal orientation of the gradients, perm can be any permutation of {\"x\",\"y\",\"z\"}.
FlipGradientOrientation[grad, flip] flips the internal orientation of the gradients, flip can be {1,1,1}, {-1,1,1}, {1,-1,1} or {1,1,-1}.
FlipGradientOrientation[grad, flip, perm] flips and permuter the internal orientation of the gradients.
FlipGradientOrientation[grad, perm, flip]flips and permuter the internal orientation of the gradients."


EigenvalCalc::usage = 
"EigenvalCalc[tensor] calculates the eigenvalues for the given tensor."

EigenvecCalc::usage =
"EigenvecCalc[tensor] calculates the eigenvectors for the given tensor."

EigensysCalc::usage = 
"EigensysCalc[tensor] calculates the eigensystem for the given tensor."

ADCCalc::usage =
"ADCCalc[eigenvalues] calculates the ADC from the given eigenvalues."

FACalc::usage =
"FACalc[eigenvalues] calculates the FA from the given eigenvalues."

ECalc::usage =
"ECalc[eigenvalues] calculates the E from the given eigenvalues."

WestinMeasures::usage = 
"WestinMeasures[eigenvalues] calculates the westin measures."

ParameterCalc::usage = 
"ParameterCalc[tensor] calculates the eigenvalues and MD and FA from the given tensor. The parameters are l1, l2, l3, MD and FA. l1, l2, l3, MD are in (10^-3 mm^2/s)."


LogTensor::usage = 
"LogTensor[tensor] transforms the tensor to LogEuclidian space.

LogTensor[] is based on DOI: 10.1109/42.963816."

ExpTensor::usage = 
"ExpTensor[tensor] transforms the tensor from LogEuclidian space.

ExpTensor[] is based on DOI: 10.1109/42.963816."


AngleCalc::usage = 
"AngleCalc[data, vector] calculates the angel between the vector and the data. Data shoud be an array of dimensions {xxx,3}."

AngleMap::usage = 
"AngleMap[data] calculates the zennith and azimuth angles of a 3D dataset (z,x,y,3) containing vectors relative to the slice direction."


DriftCorrect::usage = 
"DriftCorrect[data, bval] dirft corrects the data using the signals of the lowest bvalue that has 6 or more unique volumes.
For the function to work optimal it is best to have these volumes evenly spread throughout thet data and for the first and last volume to have this low bvalue.

DriftCorrect[] is based on DOI: 10.1002/mrm.26124."

ConcatenateDiffusionData::usage =
"ConcatenateDiffusionData[{{data1, .., dataN}, {grad1, .., gradN}, {bval, .., bvalN}, {vox, .., voxN}}] concatenates the diffusion data sets.
ConcatenateDiffusionData[{data1, .., dataN}, {grad1, .., gradN}, {bval, .., bvalN}, {vox, .., voxN}] concatenates the diffusion data sets."

SortDiffusionData::usage = 
"SortDiffusionData[data, grad, bval] sorts the diffusion datasets grad and bval for magnitude of bvalue."

RemoveIsoImages::usage = 
"RemoveIsoImages[data, grad, bval] Romoves the ISO images from the philips scanner from the data. ISO images have g={0,0,0} and b>0."


ResidualCalc::usage =
"ResidualCalc[dti,{tensor,s0},gradients,bvector] calculates the tensor residuals for the given dataset.
ResidualCalc[dti,{tensor,s0},outlier,gradients,bvector] calculates the tensor residuals for the given dataset taking in account the outliers.
ResidualCalc[dti,{tensor,s0},bmat] calculates the tensor residuals for the given dataset.
ResidualCalc[dti,{tensor,s0},outlier,bmat] calculates the tensor residuals for the given dataset taking in account the outliers.
ResidualCalc[dti,tensor,gradients,bvector] calculates the tensor residuals for the given dataset. Tensor must contain Log[s0].
ResidualCalc[dti,tensor,outlier,gradients,bvector] calculates the tensor residuals for the given dataset taking in account the outliers. Tensor must contain Log[s0].
ResidualCalc[dti,tensor,bmat] calculates the tensor residuals for the given dataset. Tensor must contain Log[s0].
ResidualCalc[dti,tensor,outlier,bmat] calculates the tensor residuals for the given dataset taking in account the outliers. Tensor must contain Log[s0]."

SigmaCalc::usage = 
"SigmaCalc[dti,grad,bvec] calculates the noise sigma based on the tensor residual, using a blur factor of 10.
SigmaCalc[dti,tens,grad,bvec] calculates the noise sigma based on the tensor residual, using a blur factor of 10.
SigmaCalc[dti,grad,bvec,blur] calculates the noise sigma based on the tensor residual, If blur is 1 ther is no blurring.
SigmaCalc[dti,tens,grad,bvec,blur] calculates the noise sigma based on the tensor residual. If blur is 1 ther is no blurring."


TransformTensor::usage = 
"TransformTensor[tensor, disp, vox] corrects the tensor with voxel size vox based on the displacement field disp. The displacement field is te displacement in mm
for each voxel location in x, y and z.

TransformTensor[] is based on DOI: 10.1109/42.963816."


Correct::usage =
"Correct[data, phase, shiftpar] corrects the dataset data using the phasemap and the shiftpar and interpolation order 1.
Correct[data, phase, shiftpar, int] corrects the dataset data using the phasemap and the shiftpar and interpolation order int."

TensorCorrect::usage=
"TensorCorrect[tensor, phase, shift, vox] corrects the tensor based on B0 field map. Can perform both translation and rotation of tensor."

Deriv::usage = 
"Deriv[disp, vox] calculates the derivative of the displacement along the three main axes. disp is the displacement field, vox is the voxel size.
Deriv[disp, vox, mask] calculates the derivative of the displacement along the three main axes. Sharp edges between the background en disp are solved by the mask. mask is a mask delining the edge of the displacement field.";


(* ::Subsection::Closed:: *)
(*Options*)


NormalizeSignal::usage = 
"NormalizeSignal is an option for DriftCorrect."


FullOutput::usage = 
"FullOutput is an option for TensorCalc when using bvector. When True also the s0 is given as output."

RobustFit::usage = 
"RobustFit is an option for TensorCalc. If true outliers will be rejected in the fit, only works with WLLS.
If FullOutput is given the outlier map is given.";

RobustFitParameters::usage =
"RobustFitParameters is an option for TensorCalc. gives the threshold for stopping the iterations and the kappa for the outlier marging, {tr,kappa}."


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
"RotationCorrect is an option for TensorCorrect. Default is False. Is a tensor is deformed setting to True also the shear is accounted for by local rotation of the tensor."


(* ::Subsection::Closed:: *)
(*Error Messages*)


TensorCalc::grad = 
"The `2` gradient directions defined do not match the `1` gradients directions in the data set."

TensorCalc::data =
"Data set dimensions (`1`D) unknown, posibilities:
- Multiple slices (4D)-> {slices, directions, x,y}
- Single slice (3D)-> {directions, x, y}
- Multiple voxels (2D)-> {directions, voxels}
- Single voxel (1D)-> {directions}"

TensorCalc::bvec = 
"the `1` gradient directions do not match the `2` b values in the b vector."

TensorCalc::met = 
"The method specified (`1`) is not a valid method, please use: \"LLS\",\"WLLS\",\"NLS\"."

ResidualCalc::datdim = 
"DTIdata (`1`) and tensor data (`2`) are not the same dimensions."

AngleCalc::dist = 
"Unknown option (`1`), options can be. \"0-180\", \"0-90\" or \"-90-90\"."

ConcatenateDiffusionData::dim = 
"data, grad and bval should be the same length:  data `1` / grad `2` / bval `2`."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*TensorCalc*)


(* ::Subsubsection::Closed:: *)
(*TensorCalc*)


Options[TensorCalc]= {
	MonitorCalc->True, 
	Method->"iWLLS", 
	FullOutput->True, 
	RobustFit->True, 
	Parallelize->True , 
	RobustFitParameters->{10.^-1, 5}
};

SyntaxInformation[TensorCalc] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

(*bvalue, only one number, gr does not have b=0*)
TensorCalc[data_, gr_, bvalue:_?NumberQ, opts:OptionsPattern[]]:=
Block[{depthD,dirD,dirG,grad,bvec},

	depthD = ArrayDepth[data];
	dirD = Length[If[depthD==4, data[[1]], data]];
	dirG = Length[gr];

	(*check if data is 4D, 3D, 2D or 1D*)
	If[depthD > 4, Return[Message[TensorCalc::data, ArrayDepth[data]]]];
	(*check if gradient dimensions are the same in the data and grad vector*)
	If[(dirD-1) != dirG,Return[Message[TensorCalc::grad,dirD,dirG]]];

	bvec = Prepend[ConstantArray[bvalue,{dirG}],0];
	grad = N[Prepend[gr,{0,0,0}]];

	If[OptionValue[Method]!="DKI",
		TensorCalc[data, Bmatrix[bvec, grad], opts],
		TensorCalc[data, Bmatrix[bvec, grad, Method->"DKI"], opts]
	]
]


(*bvector*)
TensorCalc[data_, grad_?MatrixQ, bvec:{_?NumberQ ..}, opts:OptionsPattern[]]:=
Block[{depthD,dirD,dirG,dirB},

	depthD=ArrayDepth[data];
	dirD = Length[If[depthD==4, data[[1]], data]];
	dirG = Length[grad];
	dirB = Length[bvec];

	(*check if data is 4D, 3D, 2D or 1D*)
	If[depthD>4,Return[Message[TensorCalc::data,ArrayDepth[data]]]];
	(*check if gradient dimensions are the same in the data and grad vector*)
	If[dirD!=dirG,Return[Message[TensorCalc::grad,dirD,dirG]]];
	(*check if bvec is the same length as gradient vec*)
	If[dirB!=dirG,Return[Message[TensorCalc::bvec,dirG,dirB]]];

	If[OptionValue[Method]!="DKI",
		TensorCalc[data,Bmatrix[bvec, grad], opts],
		TensorCalc[data,Bmatrix[bvec, grad, Method->"DKI"],opts]
	]
]


(*bmatrix*)
TensorCalc[dat_, bmati_?MatrixQ, OptionsPattern[]]:=Block[{
		dirD, dirB, tensor, rl, rr, TensMin, out, tenscalc, x, data, depthD, bmat, 
		fout, method, output, robust, func, con, kappa, result, dataL, outliers, parallel, 
		mon, fitFun, outFit, fitresult, residual, dataFit, s0, outFun
	},

	(*get output form*)
	{output, robust, {con,kappa}, parallel, mon, method} = OptionValue[{
		FullOutput, RobustFit, RobustFitParameters, Parallelize, MonitorCalc, Method}];

	(*chekc method*)
	If[!MemberQ[{"LLS", "WLLS", "iWLLS"(*,"NLS","GMM","CLLS","CWLLS","CNLS","DKI"*)}, method],
		Return[Message[TensorCalc::met, method];$Failed]
	];

	bmat = bmati;
	fitFun = Switch[method, "LLS", TensMinLLS, "WLLS", TensMinWLLS, "iWLLS", TensMiniWLLS];
	outFun = FindTensOutliers[#, bmat, con, kappa]&;

	(*get the data dimensions*)	
	depthD = ArrayDepth[dat];
	dirD = Length[If[depthD==4, dat[[1]], dat]];
	dirB = Length[bmat];

	(*make diff direction last dimension*)
	data = ToPackedArray@Ramp@N@Round[dat, .000001];
	data = RotateDimensionsLeft@If[depthD==4, Transpose@data, data];
	dataL = ToPackedArray@N@LogNoZero[data];

	(*check if data is 4D, 3D, 2D or 1D*)
	If[depthD>4, Return[Message[TensorCalc::data, depthD];$Failed]];
	(*check if bmat is the same length as data*)
	If[dirB!=dirD, Return[Message[TensorCalc::bvec, dirD, dirB];$Failed]];

	(*prepare for parallel computing if needed*)
	func = If[parallel && method =!= "LLS",
		DistributeDefinitions[bmat, con, kappa, fitFun, outFun,
			FindTensOutliers, TensMinLLS, TensMinWLLS, TensMiniWLLS];
		ParallelMap, 
		Map];	

	(*define outliers if needed*)	
	outliers = If[robust && method =!= "LLS",
		If[mon, PrintTemporary["Finding tensor outliers"]]; 
		If[depthD == 1,	outFun[dataL], func[outFun, dataL]],
		SparseArray[{}, Dimensions@data, 0.]
	];
	outFit = ToPackedArray@N@(1. - outliers);

	If[mon, PrintTemporary["Fitting tensor"]];
	fitresult = If[method === "LLS", 
		RotateDimensionsLeft[PseudoInverse[bmat].RotateDimensionsRight[dataL]],
		If[depthD == 1,
			(*single voxel fit*)
			fitFun[outFit data, outFit dataL, bmat],
			(*2-4D data fit*)
			dataFit = RotateDimensionsLeft[{RotateDimensionsRight[outFit data], RotateDimensionsRight[outFit dataL]}, 2];
			func[fitFun[#, bmat]&, dataFit]
		]
	];

	If[mon, PrintTemporary["Finalizing output tensor"]]; 

	fitresult = RotateDimensionsRight[fitresult];
	outliers = RotateDimensionsRight[outliers];
	If[depthD == 4,	outliers = Transpose[outliers]];

	If[OptionValue[FullOutput], 
		residual = ResidualCalc[dat, fitresult, outliers, bmat, MeanRes->"RMSE"]
	];

	s0 = N@Clip[ExpNoZero[N@Chop[Last[fitresult]]],{0., 1.5 Max[data]}];
	tensor = N@Clip[Drop[fitresult, -1],{-0.1,0.1}];

	If[OptionValue[FullOutput],If[robust,{tensor, s0, outliers, residual}, {tensor, s0, residual}], {tensor, s0}]
]


(* ::Subsubsection::Closed:: *)
(*FindOutliers*)


FindTensOutliers = Quiet@Compile[{{ls, _Real, 1}, {bmat, _Real, 2}, {con, _Real, 0}, {kappa, _Real, 0}}, Block[{
	sol, ittA, contA, solA, itt, cont, soli, res, mad, wts, wmat, fitE, LS2, bmat2, out},

	(*based on DOI: 10.1002/mrm.25165*)
	(*initialize some values*)
	out = (0. ls); 
	LS2 = ls; 
	bmat2 = bmat;

	(*If not background find the outliers*)
	If[Total[ls] >1 ,
		(*Step1: initial LLS fit*)
		sol = PseudoInverse[bmat] . ls;
		(*check if LLS fit is plausable, i.e. s0 > 0*)
		If[Last[sol] > 0,

			(*Start outer loop*)
			ittA = 0; contA = 1;
			While[contA == 1,
				(*init the solution*)
				solA = sol;
				ittA++;

				(*Step2: Compute a robust estimate for homoscedastic regression using IRLS.*)
				itt = 0; cont = 1;
				While[cont == 1, itt++;
					(*init the solution*)
					soli = sol;
					(*a. Calculate the residuals e* in the linear domain*)
					res = ls - bmat . sol;
					(*b. Obtain an estimate of the dispersion of the residuals by calculating the median absolute deviation (MAD).*)
					mad = 1.4826 MedianDeviation[res];
					(*prevent calculation with 0*)
					If[mad === 0., cont = 0,
						(*c. Recompute the weights according to Eq. [13].*)
						wts = Normalize[1 / (1 + (res/mad)^2)^2];
						(*d. Perform WLLS fit with new weights*)
						wmat = Transpose[bmat] . DiagonalMatrix[wts];
						sol = PseudoInverse[wmat . bmat] . wmat . ls;
						(*e. Check convergence*)
						If[Total[UnitStep[Abs[sol - soli] - con Abs[soli]]] === 0 || itt === 3, cont = 0];
					];
				];(*end first while*)

				(*Step 3: Transform variables for heteroscedasticity*)
				fitE = Exp[-bmat . Chop[sol]] + 10^-10;
				LS2 = ls / fitE;
				bmat2 = bmat / fitE;

				(*Step 4: Initial LLS fit in * domain*)
				sol = PseudoInverse[bmat2] . LS2;

				(*Step 5: Compute a robust estimate for homoscedastic regression using IRLS.*)
				itt = 0; cont = 1;
				While[cont == 1, 
					(*init the solution*)
					soli = sol;
					itt++;
					(*a. Calculate the residuals e* in the linear domain*)
					res = LS2 - bmat2 . sol;
					(*b. Obtain an estimate of the dispersion of the residuals by calculating the median absolute deviation (MAD).*)
					mad = 1.4826 MedianDeviation[res];
					(*prevent calculation with 0*)
					If[mad === 0., 
						cont = 0,
						(*c. Recompute the weights according to Eq. [13].*)
						wts = Normalize[1 / (1 + (res/mad)^2)^2];
						(*d. Perform WLLS fit with new weights*)
						wmat = Transpose[bmat2] . DiagonalMatrix[wts];
						sol = PseudoInverse[wmat . bmat2] . wmat . LS2;
						(*e. Check convergence*)
						If[Total[UnitStep[Abs[sol - soli] - con Abs[soli]]] === 0 || itt === 3, cont = 0];
					];
				];(*end second while*)

				(*Step 6: Check convergence overall loop*)
				If[Total[UnitStep[Abs[sol - solA] - con Abs[solA]]] === 0 || ittA === 5, contA = 0];
			];(*end main while*)

			(*Step 7: Identify and exclude outliers*)
			res = LS2 - bmat2 . sol;
			out = UnitStep[Abs[res] - (kappa 1.4826 MedianDeviation[res])];

			];(*close if negative s0*)
		];(*close if background*)

		out
	],
	RuntimeAttributes -> {Listable}, 
	RuntimeOptions -> {"Speed", "WarningMessages" -> False}
]


(* ::Subsubsection::Closed:: *)
(*LLS*)


TensMinLLS = Compile[{{dat, _Real, 2}, {bmatI, _Real, 2}},
	bmatI . dat[[2]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]


(* ::Subsubsection::Closed:: *)
(*WLLS*)


TensMinWLLS = Compile[{{dat, _Real, 2}, {bmat, _Real, 2}}, 
	Block[{wmat, mvec, sol, s, ls},
		{s, ls} = dat;
		mvec = UnitStep[s] Unitize[s];
		sol = First[bmat];
		wmat = 0. bmat;
		If[! (Total[ls]=== 0. || Total[mvec] < 7),
			wmat = Transpose[bmat] . DiagonalMatrix[mvec s^2];
			sol = PseudoInverse[wmat . bmat] . wmat . ls;
		];
	sol]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]


(* ::Subsubsection::Closed:: *)
(*iWLLS*)


TensMiniWLLS = Compile[{{dat, _Real, 2}, {bmat, _Real, 2}}, 
	Block[{s, ls, wmat, mat, cont, itt, mvec, soli, sol0, max, sol, w},
		{s, ls} = dat;
		mvec = UnitStep[s] Unitize[s];
		max = 3. Max[mvec s];
		sol0 = 0. First[bmat];
		sol = sol0;
		wmat = 0. bmat;

		(*skip background or not enough data for fit*)
		If[!Total[ls]=== 0.||Total[mvec] > 7,
			(*initialize*)
			itt = 0;
			cont = 1.;
			(*initialize using LLS*)
			sol = PseudoInverse[bmat] . ls;
			(*check for implausabole solution (negative s0 or high s0)*)
			If[Last[sol] >= max || Last[sol] <= 0.,
				sol = sol0;
				,
				(*itterative reweighting*)
				While[cont == 1,
					(*init iteration values*)
					itt++;
					soli = sol;
					(*perform WLLS*)
					w = mvec Exp[2 bmat . sol];
					wmat =Transpose[bmat] . DiagonalMatrix[w];
					sol = PseudoInverse[wmat . bmat] . wmat . ls;
					(*update weight*)
					(*see if to quit loop*)
					If[(Last[sol] >= max || Last[sol] <= 0), cont = 0.;	sol = sol0];
					If[! AnyTrue[Abs[sol - soli] - 0.001 Abs[soli], Positive] || itt === 10 , cont = 0.];
				]
			]
		];
		sol
	]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages"->False}]


(* ::Subsubsection::Closed:: *)
(*DKI*)


(*TensMinDKI[s_,ls_,bmat_,bmatI_]:=bmatI.ls*)
TensMinDKI = Compile[{{s, _Real, 1}, {bmatI, _Real, 2}},
	If[Total[s]==0.,
		{0.,0.,0.,0.,0.,0.,0.},
		bmatI . s
	],RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


(* ::Subsubsection::Closed:: *)
(*NLS*)


TensMinNLS[s_,ls_,bmat_,bmatI_]:=
Module[{v,xx,yy,zz,xy,xz,yz,init,tens,sol},
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},
		tens,
		v={xx,yy,zz,xy,xz,yz,tens[[7]]};
		init=Thread[{v[[1;;6]],tens[[1;;6]]}];
		sol=FindMinimum[.5 Total[(s-Exp[bmat . v])^2],init][[2]];
		v/.sol
	]
]


(* ::Subsubsection::Closed:: *)
(*NLS*)


TensMinGMM[s_,ls_,bmat_,bmatI_]:=
Module[{v,xx,yy,zz,xy,xz,yz,init,tens,res,w},
	s;
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={xx,yy,zz,xy,xz,yz,tens[[7]]};
		init=Thread[{v[[1;;6]],tens[[1;;6]]}];
		v/.FindMinimum[(
			res=ls-bmat . v;
			w=1/(res^2+Mean[res]^2);
			.5 Total[(w/Mean[w])*(res)^2]
		(*w=1/(res^2+(1.4826*Median[Abs[res-Median[res]]])^2);*)
		),init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CLLS*)


TensMinCLLS[s_,ls_,bmat_,bmatI_]:=
Module[{v,r0,r1,r2,r3,r4,r5,init,tens},
	s;
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={r0^2,r1^2+r3^2,r2^2+r4^2+r5^2,r0 r3,r0 r4,r3 r4+r1 r5,tens[[7]]};
		init=Thread[{{r0,r1,r2,r3,r4,r5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(ls-bmat . v)^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CWLLS*)


TensMinCWLLS[s_,ls_,bmat_,bmatI_]:=
Module[{v,r0,r1,r2,r3,r4,r5,init,tens,std=1,wmat},
	bmatI;
	wmat=Transpose[bmat] . DiagonalMatrix[s^2/std^2];
	tens=PseudoInverse[wmat . bmat] . wmat . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={r0^2,r1^2+r3^2,r2^2+r4^2+r5^2,r0 r3,r0 r4,r3 r4+r1 r5,tens[[7]]};
		init=Thread[{{r0,r1,r2,r3,r4,r5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(s^2/std^2)*(ls-bmat . v)^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*CNLS*)


TensMinCNLS[s_,ls_,bmat_,bmatI_]:=
Module[{v,r0,r1,r2,r3,r4,r5,init,tens},
	tens=bmatI . ls;
	If[tens=={0.,0.,0.,0.,0.,0.,0.},tens,
		v={r0^2,r1^2+r3^2,r2^2+r4^2+r5^2,r0 r3,r0 r4,r3 r4+r1 r5,tens[[7]]};
		init=Thread[{{r0,r1,r2,r3,r4,r5},TensVec[ExtendedCholeskyDecomposition[TensMat[tens]]]}];
		v/.FindMinimum[.5Total[(s-Exp[bmat . v])^2],init][[2]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*ExtendeCholeskyDecomposition*)


ExtendedCholeskyDecomposition[tm_]:= Block[{n,beta,theta,cm,lm,dm,em,j},
	n=Length[tm];
	beta=Max[{Max[Diagonal[tm]],Max[UpperTriangularize[tm,1]]/Sqrt[n^2-1],10^-15}];
	cm=DiagonalMatrix[Diagonal[tm]];
	lm=dm=em=ConstantArray[0,{n,n}];
	Tabel[
		If[j==1,
			(*j=1 maak eerste colom cm gelijk aan tm*)
			cm[[j+1;;,j]]=tm[[j+1;;,j]];
			,
			(*j>1 vul lm matrix*)
			lm[[j,;;j-1]]=cm[[j,;;j-1]]/(Diagonal[dm][[;;j-1]]/.(0.->Infinity));
			If[j<n,
				cm[[j+1;;,j]]=tm[[j+1;;,j]]-lm[[j,j-1;;]] . Transpose[cm[[j+1;;,j-1;;]]]
				];
			];
		theta=If[j==n,0,Max[Abs[cm[[j+1;;,j]]]]];
		dm[[j,j]]=Max[{Abs[cm[[j,j]]],theta^2/beta}];
		em[[j,j]]=dm[[j,j]]-cm[[j,j]];
		cm=cm-DiagonalMatrix[PadLeft[(1/(dm[[j,j]]/.(0.->Infinity)))*cm[[j+1;;,j]]^2,n]];
	,{j,1,3}];
	lm=lm+IdentityMatrix[n];
	Transpose[lm . MatrixPower[dm,.5]]
]


(* ::Subsection:: *)
(*Reorient Tensor*)


(* ::Subsubsection::Closed:: *)
(*FlipTensorOrientation*)


SyntaxInformation[FlipTensorOrientation] = {"ArgumentsPattern" -> {_, _, _.}};

FlipTensorOrientation[tensor_, p_] /; AllTrue[p, NumberQ] := FlipTensorOrientation[tensor, {"x", "y", "z"}, p]

FlipTensorOrientation[tensor_, v_] /; AllTrue[v, StringQ] := FlipTensorOrientation[tensor, v, {1, 1, 1}]

FlipTensorOrientation[tensor_, p_, v_]/; (AllTrue[v, StringQ] && AllTrue[p, NumberQ]):=FlipTensorOrientation[tensor, v, p]

FlipTensorOrientation[tensor_, v_, p_]/; (AllTrue[v, StringQ] && AllTrue[p, NumberQ]):= Block[{times, transp},
	If[DeleteDuplicates[Abs[p]] === {1} && Sort[v] === {"x", "y", "z"},
		times = Join[{1, 1, 1}, Flatten[Table[p[[i]] p[[j]], {i, 1, 3}, {j, i + 1, 3}]]];
		transp = (StringJoin[Sort[Characters[#]]] & /@ (
			StringReplace[{"xx", "yy", "zz", "xy","xz", "yz"}, Thread[{"x", "y", "z"} -> v]])
		) /. Thread[{"xx", "yy", "zz", "xy", "xz", "yz"} -> {1, 2, 3, 4, 5, 6}];
		(times tensor)[[transp]]
		,
		$Failed
	]
]


(* ::Subsubsection::Closed:: *)
(*FlipTensorOrientation*)


FlipGradientOrientation[grad_, p_] /; AllTrue[p, NumberQ] := FlipGradientOrientation[grad, {"x", "y", "z"}, p]

FlipGradientOrientation[grad_, v_] /; AllTrue[v, StringQ] := FlipGradientOrientation[grad, v, {1, 1, 1}]

FlipGradientOrientation[grad_, p_, v_] /; (AllTrue[v, StringQ] && AllTrue[p, NumberQ]) := FlipGradientOrientation[grad, v, p]

FlipGradientOrientation[grad_, v_, p_] /; (AllTrue[v, StringQ] && AllTrue[p, NumberQ]) := Block[{transp, times},
	If[DeleteDuplicates[Abs[p]] === {1} && Sort[v] === {"x", "y", "z"},
		times = ConstantArray[p, Length[grad]];
		transp = v /. Thread[{"x", "y", "z"} -> {1, 2, 3}];
		(times grad)[[All, transp]]
		,
		$Failed
	]
]


(* ::Subsection:: *)
(*Tensor Parameters*)


(* ::Subsubsection::Closed:: *)
(*EigensysCalc*)


Options[EigensysCalc]={RejectMap->False,Reject->True, PerformanceGoal->"Quality"};

SyntaxInformation[EigensysCalc]={"ArgumentsPattern"->{_,OptionsPattern[]}};

EigensysCalc[tens_,opts:OptionsPattern[]]:=EigenSys[tens,"all",opts]


(* ::Subsubsection::Closed:: *)
(*EigenvalCac*)


Options[EigenvalCalc]=Options[EigensysCalc];

SyntaxInformation[EigenvalCalc]={"ArgumentsPattern"->{_,OptionsPattern[]}};

EigenvalCalc[tens_,opts:OptionsPattern[]]:=EigenSys[tens,"val",opts]


(* ::Subsubsection::Closed:: *)
(*EigenvecCalc*)


Options[EigenvecCalc]=Options[EigensysCalc];

SyntaxInformation[EigenvecCalc]={"ArgumentsPattern"->{_,OptionsPattern[]}};

EigenvecCalc[tens_,opts:OptionsPattern[]]:=EigenSys[tens,"vec",opts]


(* ::Subsubsection:: *)
(*EigenSys*)


Options[EigenSys]=Options[EigensysCalc];

EigenSys[tens_,out_,OptionsPattern[]]:=Block[{t, met, val, vec,reject, sel},
	met = OptionValue[PerformanceGoal];

	t=Which[
		VectorQ[tens], tens,
		MatrixQ[tens]&&Dimensions[tens]==={3,3}, TensVec[tens],
		(*ArrayQ[tens]*)True, RotateDimensionsLeft[tens]
	];

	{val,vec} = If[met==="Speed", EigenSysC[t,out], EigenSysQ[t,out]];

	If[OptionValue[Reject],
		reject = SelectEig[val];
		sel = 1-reject;
		val = sel val;
		If[vec=!=0,vec=sel vec+reject ConstantArray[{{0.,0.,1.},{0.,1.,0.},{1.,0.,0.}},Dimensions[reject]]];
	];

	If[OptionValue[RejectMap]&&OptionValue[Reject],
		Switch[out,
			"val",{val,reject},
			"vec",{vec,reject},
			"all",{val,vec,reject}
		],
		Switch[out,
			"val",val,
			"vec",vec,
			"all",{val,vec}
		]
	]
]


SelectEig=Compile[{{eig,_Real,1}},1-UnitStep[Last[eig]],RuntimeAttributes->{Listable},RuntimeOptions->"Speed",Parallelization->True];


(* ::Subsubsection::Closed:: *)
(*EigenSysQ*)


(*slow precice method*)
EigenSysQ[tens_,out_]:=Block[{val,vec},
	If[out=!="val",
		If[VectorQ[tens],
			(*tensor is just one value*)
			EigenSysi[tens],
			(*calculate the eigensystem*)
			val=Map[EigenSysi,tens,{-2}];
			vec=Switch[ArrayDepth[tens]-1,1,val[[All,2]],2,val[[All,All,2]],3,val[[All,All,All,2]]];
			val=Switch[ArrayDepth[tens]-1,1,val[[All,1]],2,val[[All,All,1]],3,val[[All,All,All,1]]];
			{val,vec}
		],
		If[VectorQ[tens],
			(*tensor is just one value*)
			{EigenVali[tens], 0},
			(*calculate the eigensystem*)
			{Map[EigenVali, tens, {-2}], 0}
		]
	]
];


EigenSysi[{0.,0.,0.,0.,0.,0.}]:={{0.,0.,0.},{{0.,0.,1.},{0.,1.,0.},{1.,0.,0.}}}
EigenSysi[tensor_?VectorQ]:=Eigensystem[TensMat[tensor]]


EigenVali[{0.,0.,0.,0.,0.,0.}]:={0.,0.,0.}
EigenVali[tensor_?VectorQ]:=Eigenvalues[TensMat[tensor]]


(* ::Subsubsection::Closed:: *)
(*EigenSysC*)


(*fast direct method*)
EigenSysC[tens_, out_]:=Block[{val},
	val=EigenValC[tens];
	If[out=!="val",{val,EigenVecC[tens,val]},{val,0}]
]


(* ::Subsubsection::Closed:: *)
(*EigenValC*)


EigenValC = Compile[{{tens,_Real,1}},Block[{
		dxx,dyy,dzz,dxy,dxz,dyz,dxy2,dxz2,dyz2,i1,i2,i3,i,v,s,p3,v2,phi
	},
	(*method https://doi.org/10.1016/j.mri.2009.10.001*)
	If[Total[Abs[tens]]<10.^-15,
		{0.,0.,0.},

		{dxx,dyy,dzz,dxy,dxz,dyz}=tens;
		{dxy2,dxz2,dyz2}={dxy,dxz,dyz}^2;

		i1=dxx+dyy+dzz;
		i2=dxx dyy+dxx dzz+dyy dzz-dxy2-dxz2-dyz2;
		i3=dxx dyy dzz+2 dxy dxz dyz-dzz dxy2-dyy dxz2-dxx dyz2;

		i=i1/3;
		v=i^2-i2/3;
		s=i^3-(i1 i2)/6+i3/2;
		p3=Pi/3;
		v2=Sqrt[v];

		phi=Re[ArcCos[If[(v v2)===0,0,s/(v v2)]]/3];
		{i+2 v2 Cos[phi],i-2 v2 Cos[p3+phi],i-2 v2 Cos[p3-phi]}
	]
], RuntimeAttributes->{Listable}, RuntimeOptions -> {"Speed", "WarningMessages"->False}];


(* ::Subsubsection::Closed:: *)
(*EigenVecC*)


EigenVecC=Compile[{{tens,_Real,1},{eig,_Real,1}},Block[{dxx,dyy,dzz,dxy,dxz,dyz,a,b,c,norm},
	If[Total[Abs[tens]]<10.^-15,
		{{0,0,1},{0,1,0},{1,0,0}},

		{dxx,dyy,dzz,dxy,dxz,dyz}=tens;
		(
			{a,b,c}={dxz dxy,dxy dyz,dxz dyz}-{dyz,dxz,dxy} ({dxx,dyy,dzz}-#);
			{a,b,c}={b c,a c,a b};
			norm = Sqrt[Abs[a]^2 + Abs[b]^2 + Abs[c]^2];
			If[norm===0,{a,b,c}/norm,{0.,0.,0.}]
		)&/@eig
	]
], RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"];


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
		If[OptionValue[MonitorCalc], Print["Done calculating e for "<>ToString[slices]<>" slices!"]];
		,
		output=ECalci[eigen];
		If[ArrayQ[eigen,3],If[OptionValue[MonitorCalc], Print["Done calculating e for 1 slice!"]],
			If[VectorQ[eigen],If[OptionValue[MonitorCalc], Print["Done calculating e for 1 voxel!"]]]
			]
		];
	Return[output];
	]


ECalci[eigen_]:= Block[{ec},
	ec=Compile[{l1,l3},Sqrt[1-(l3/l1)]];
	Map[If[#[[3]]!=0, ec[#[[1]],#[[3]]], 0]&,eigen,{ArrayDepth[eigen]-1}]
]


(* ::Subsubsection::Closed:: *)
(*WestinMeasures*)


WestinMeasures[eig_]:=Block[{l1, l2, l3},
	{l1,l2,l3} = RotateDimensionsRight[eig];
	{DivideNoZero[l1-l2,l1], DivideNoZero[l2-l3,l1], DivideNoZero[l3,l1]}
]


(* ::Subsubsection::Closed:: *)
(*ParameterCalc*)


Options[ParameterCalc] = {Reject->False, PerformanceGoal -> "Quality"}

SyntaxInformation[ParameterCalc] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ParameterCalc[tensor_,OptionsPattern[]]:= Block[{eig,adc,fa},
	eig = 1000 EigenvalCalc[tensor, Reject->OptionValue[Reject], RejectMap->False, PerformanceGoal->OptionValue[PerformanceGoal]];
	adc = ADCCalc[eig];
	fa = FACalc[eig];
	Join[RotateDimensionsRight[eig], {adc,fa}]
]


(* ::Subsection:: *)
(*LogEuclidian *)


(* ::Subsubsection::Closed:: *)
(*LogTensor*)


SyntaxInformation[LogTensor]={"ArgumentsPattern"->{_}};

LogTensor[tens_]:=Block[{t,v,e},
	t=TensMat[tens];
	t=Map[(
		If[Total[Flatten[#]]===0.,#,{v,e}=Eigensystem[#];
		If[AnyTrue[v,#<=0.&],0.#,
		Transpose[e] . DiagonalMatrix[Log[v]] . e]]
	)&,t,{-3}];
	{1.,1.,1.,Sqrt[2.],Sqrt[2.],Sqrt[2.]}TensVec[t]
]


(* ::Subsubsection::Closed:: *)
(*ExpTensor*)


SyntaxInformation[ExpTensor]={"ArgumentsPattern"->{_}};

ExpTensor[tens_]:=Block[{t,e,v},
	t=TensMat[{1.,1.,1.,1./Sqrt[2.],1./Sqrt[2.],1./Sqrt[2.]}tens];
	t=Map[(
		If[Total[Flatten[#]]===0.,#,{v,e}=Eigensystem[#];
		Transpose[e] . DiagonalMatrix[Exp[v]] . e]
	)&,t,{-3}];
	TensVec[t]
]


(* ::Subsection:: *)
(*Angle maps*)


(* ::Subsubsection::Closed:: *)
(*AngleCalc*)


Options[AngleCalc]={Distribution->"0-180"};

SyntaxInformation[AngleCalc] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};


AngleCalc[data_?ArrayQ,vec_?VectorQ,OptionsPattern[]]:=
Module[{angles},
	angles=Map[If[Re[#]==#,ArcCos[# . vec],"no"]&,data,{Depth[data]-2}];

	Switch[
		OptionValue[Distribution],
		"0-180",
		(180/Pi)angles,
		"0-90",
		Map[If[#=!="no",If[#>=1/2Pi,(180/Pi)(Pi-#),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
		"-90-90",
		Map[If[#=!="no",If[#>=1/2Pi,(180/Pi)(#-Pi),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
		_,
		Message[AngleCalc::dist,OptionValue[Distribution]]
		]
	]


AngleCalc[data_?ArrayQ,vec_?ArrayQ,OptionsPattern[]]:=
Module[{angles},
	If[Dimensions[data]!=Dimensions[vec],
		Print["Error"],

		angles=MapThread[If[Re[#1]==#1,ArcCos[#1 . #2],"no"]&,{data,vec},ArrayDepth[vec]-1];

		Switch[
			OptionValue[Distribution],
			"0-180",
			(180/Pi)angles,
			"0-90",
			Map[If[#=!="no",If[#>=1/2Pi,(180/Pi)(-(#-Pi)),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
			"-90-90",
			Map[If[#=!="no",If[#>=1/2Pi,(180/Pi)(#-Pi),(180/Pi)(#)]]&,angles,{ArrayDepth[angles]}],
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


(* ::Subsection:: *)
(*Diffusion data functions*)


(* ::Subsubsection::Closed:: *)
(*DriftCorrect*)


Options[DriftCorrect]={NormalizeSignal->True, UseMask->True}

SyntaxInformation[DriftCorrect] = {"ArgumentsPattern" -> {_, _,_., OptionsPattern[]}}

DriftCorrect[data_, bi_, opts:OptionsPattern[]] := Block[{bval,pos},
	bval = If[ArrayDepth[bi] == 2, BmatrixInv[bi][[1]], bi];
	pos = First@UniqueBvalPosition[bval, 5][[2]];
	DriftCorrect[data, bval, pos, opts]
];

DriftCorrect[data_, bi_, pos_, OptionsPattern[]] := Block[{
	sig, cor, bval, sol1, sol2, sol3, a, b, c, x, outp, dat
	},
	bval = If[ArrayDepth[bi] == 2, BmatrixInv[bi][[1]], bi];
	sig = MeanSignal[data, pos, UseMask->OptionValue[UseMask]];
	dat = Transpose[{pos, sig}];

	{sol1, sol2, sol3} = {a, b, c} /. FindFit[dat, {c + b x + a x^2}, {a, b, c}, x];
	cor = sol3/Table[sol3 + sol2 x + sol1 x^2, {x, 1, Length[bi]}];

	outp = ConstantArray[cor, Length[data]] data;

	If[OptionValue[NormalizeSignal], 100 outp / (sig[[1]] cor[[1]]) , outp]
];


(* ::Subsubsection::Closed:: *)
(*ConcatenateDiffusionData*)


SyntaxInformation[ConcatenateDiffusionData] = {"ArgumentsPattern" -> {_, _., _., _.}};

ConcatenateDiffusionData[data_?ListQ] :=If[Length[data] == 4,ConcatenateDiffusionData[data[[1]], data[[2]], data[[3]], data[[4]]]]

ConcatenateDiffusionData[data_, grad_, val_, vox_] := Module[{dataout, gradout, valout, voxout},
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

	{ToPackedArray@N@dataout, gradout, valout, voxout}
];


(* ::Subsubsection::Closed:: *)
(*SortDiffusionData*)


SyntaxInformation[SortDiffusionData] = {"ArgumentsPattern" -> {_, _, _}};

SortDiffusionData[data_, grad_, val_] := Module[{pos, valu, sel},
	{valu, pos} = UniqueBvalPosition[val];
	sel = Flatten[pos];
	{ToPackedArray@N@data[[All, sel]], grad[[sel]], val[[sel]]}
]


(* ::Subsubsection::Closed:: *)
(*RemoveIsoImages*)


SyntaxInformation[RemoveIsoImages] = {"ArgumentsPattern" -> {_, _, _}};

RemoveIsoImages[data_, grad_, val_] := Module[{sel},
	sel = Complement[
		Range[Length[val]],Complement[Flatten[Position[grad, {0., 0., 0.}]], Flatten[Position[val, 0.]]]
	];
	{ToPackedArray@N@data[[All, sel]], grad[[sel]], val[[sel]]}
]


(* ::Subsection:: *)
(*Tensor Residuals*)


(* ::Subsubsection::Closed:: *)
(*ResidualCalc*)


Options[ResidualCalc] = {MeanRes -> "All"};

SyntaxInformation[ResidualCalc] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

(*b0, no outliers, bval, bvec*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, s0_?ArrayQ}, grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}, bval_, opts : OptionsPattern[]] := 
	ResidualCalc[data, Join[tensor, {LogNoZero[s0]}], ConstantArray[0., Dimensions[data]], Bmatrix[bval, grad], opts]

(*b0, no outliers bmat*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, s0_?ArrayQ}, bmat_?ArrayQ, opts : OptionsPattern[]] :=
	ResidualCalc[data, Join[tensor, {LogNoZero[s0]}], ConstantArray[0., Dimensions[data]], bmat, opts]

(*b0, outliers, bval, bvec*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, s0_?ArrayQ}, outlier_?ArrayQ, grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}, bval_, opts : OptionsPattern[]] :=
	ResidualCalc[data, Join[tensor, {LogNoZero[s0]}], outlier, Bmatrix[bval, grad], opts]

(*b0, outliers, bmat*)
ResidualCalc[data_?ArrayQ, {tensor_?ArrayQ, s0_?ArrayQ}, outlier_?ArrayQ, bmat_?ArrayQ, opts : OptionsPattern[]] :=
	ResidualCalc[data, Join[tensor, {LogNoZero[s0]}], outlier, bmat, opts]

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
	dimD = Dimensions[If[ArrayDepth[dat] == 4,dat[[All,1]],dat[[1]]]];
	dimT = Dimensions[tensor[[1]]];

	If[dimD != dimT || Length[tensor]!=7, Return[Message[ResidualCalc::datdim, Dimensions[data], Dimensions[tensor]]]];

	(*remove ouliers*)
	fit = Clip[ExpNoZero[bmat . tensor], {-1.5, 1.5} Max[dat], {0., 0.}];

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


(* ::Subsubsection::Closed:: *)
(*SigmaCalc*)


Options[SigmaCalc] = {FilterShape -> "Median"};

SyntaxInformation[SigmaCalc] = {"ArgumentsPattern" -> {_, _, _, _., _., OptionsPattern[]}};

SigmaCalc[dti_?ArrayQ, grad : {{_, _, _} ..}, bvalue_, blur_: 2, OptionsPattern[]] := Module[
	{tens, res,len,sig}, 
	tens = TensorCalc[dti, grad, bvalue, MonitorCalc -> False];
	res = ResidualCalc[dti, tens, grad, bvalue, MeanRes -> "MAD"];
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


SigmaCalc[dti_?ArrayQ, tens_?ArrayQ, grad : {{_, _, _} ..}, bvalue_, blur_: 2, OptionsPattern[]] := Module[
	{res,len,sig},
	res = ResidualCalc[dti, tens, grad, bvalue, MeanRes -> "MAD"];
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
(*TensorCorrect*)


(* ::Subsubsection::Closed:: *)
(*TransformTensor*)


SyntaxInformation[TransformTensor] = {"ArgumentsPattern" -> {_, _, _}};

TransformTensor[tens_, disp_, vox_]:=Block[{imat, jac},
	imat=IdentityMatrix[3];
	jac=Chop[imat+Table[GaussianFilter[disp[[i]],1,imat[[j]]]/vox[[i]],{i,1,3},{j,1,3}]];

	TensVec[Apply[TensorRotate,RotateDimensionsLeft[{RotateDimensionsRight[TensMat[tens],2],jac},3],{-4}]]
]


(* ::Subsubsection::Closed:: *)
(*TensorRotate*)


TensorRotate[tens_,f_]:=Block[{val,e1,e2,e3,n1,n2,n3,nMat,fMat},
	If[tens[[1,1]]==0.,
		tens,
		{val, {e1,e2,e3}}=Eigensystem[tens];
		fMat=PseudoInverse[f];
		n1=Normalize[fMat . e1];
		n2=Normalize[fMat . e2-(n1 . (fMat . e2))*n1]//N;
		n3=Normalize[Cross[n1,n2]]//N;
		nMat={n1,n2,n3};
		Chop[Transpose[nMat] . DiagonalMatrix[val] . nMat]
	]
];


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


Correcti[dat_,ph_,shift_,int_]:= Module[{pos,acpos,shiftpx,data,phase,output},
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
	Module[{dim,pxshift,der,f,tensM,tensC,tensCV,tensT},

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
		f=Fmat[der,shift[[2]]];

		PrintTemporary["Rotation Correction"];
		(*rotation correction of matrix*)
		(*tensor to matrixform*)
		tensM=TensMat[tens];
		(*rotation correct tensor matrix*)
		tensC=MapThread[DRot[#1,#2]&,{tensM,f},3];
		(*corrected tensor back to vector form*)
		tensCV=TensVec[tensC];
		,
		tensCV=tens;
		];

	PrintTemporary["Translation Correction"];
	(*Translation correction of the rotation corrected Tensor*)
	tensT=Map[(MapThread[TransCorrect[#1,#2,shift[[2]],1]&, {#,pxshift}])&,tensCV]
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
Module[{dx,dy,dz,dim,zero,ones,f},
	{dx,dy,dz}=der;
	dim=Dimensions[dx];
	zero=ConstantArray[0,dim];
	ones=ConstantArray[1,dim];
	If[shift=="COL",
		f=Transpose[{{ones,zero,zero},{dx,dy+1,dz},{zero,zero,ones}},{4,5,1,2,3}],
		If[shift=="ROW",
			f=Transpose[{{dx+1,dy,dz},{zero,ones,zero},{zero,zero,ones}},{4,5,1,2,3}]
			],
		Print["error, unknown direction"]
		]
	];


(* ::Subsubsection::Closed:: *)
(*Drot*)


DRot[tens_,f_]:=Module[{val,e1,e2,e3,n1,n2,n3,nn},
	{val,{e1,e2,e3}}=Eigensystem[tens];
	n1=Normalize[f . e1];
	n2=Normalize[f . e2-(n1 . (f . e2))*n1]//N;
	n3=Normalize[Cross[n1,n2]]//N;
	nn=Transpose[{n1,n2,n3}];
Chop[nn . (IdentityMatrix[3]val) . Transpose[nn]]
];


(* ::Subsection:: *)
(*Deriv*)


(* ::Subsubsection::Closed:: *)
(*Deriv*)


SyntaxInformation[Deriv] = {"ArgumentsPattern" -> {_, _, _}};

Deriv[disp_,vox_]:=
Module[{dim,dx,dy,dz},
	dim=Dimensions[disp];
	dx=Transpose[DerivFunc[Transpose[disp,{1,3,2}],dim[[2]],vox[[2]]],{1,3,2}];
	dy=DerivFunc[disp,dim[[3]],vox[[3]]];
	dz=Transpose[DerivFunc[Transpose[disp,{3,2,1}],dim[[1]],vox[[1]]],{3,2,1}];
	{dx,dy,dz}
	];

Deriv[disp_,vox_,mask_]:=
Module[{dim,dx,dy,dz},
	dim=Dimensions[disp];
	dx=Transpose[DerivFunc[Transpose[disp,{1,3,2}],Transpose[mask,{1,3,2}],dim[[2]],vox[[2]]],{1,3,2}];
	dy=DerivFunc[disp,mask,dim[[3]],vox[[3]]];
	dz=Transpose[DerivFunc[Transpose[disp,{3,2,1}],Transpose[mask,{3,2,1}],dim[[1]],vox[[1]]],{3,2,1}];
	{dx,dy,dz}
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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
