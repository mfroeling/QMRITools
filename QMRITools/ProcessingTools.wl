(* ::Package:: *)

(* ::Title:: *)
(*QMRITools ProcessingTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`ProcessingTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`ProcessingTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


SetupDataStructure::usage = 
"SetupDataStructure[dcmFolder] makes nii folders and generates nii files for a directory of dmc data where the data is structured per subject."

SNRCalc::usage = 
"SNRCalc[data,masksig,masknoise] calculates the Signal to noise ratio of the signal selected by masksig and the noise selected by masknoise."

SNRMapCalc::usage =
"SNRMapCalc[data1,noisemap] calcualtes the signal to noise ratio of the data using MN[data]/(1/sqrt[pi/2] sigma), \
where sigma is the local mean of the noise map assuming it is a rician distribution.

SNRMapCalc[{data1,data2}] calcualtes the signal to noise ratio from two identical images using \
MN[data1,data2] / (.5 SQRT[2] STDV[data2-data1]).

SNRMapCalc[{data1, .. dataN}] calcualtes the signal to noise ratio of the data using MN/sigma where the mean signal MN is the average voxe \
value over all dynamics N and the sigma is the standard deviation over all dynamics N."

FitData::usage = 
"FitData[data,range] converts the data into 100 bins within the +/- range around the mean. Function is used in ParameterFit."

ParameterFit::usage = 
"ParameterFit[data] fits a (skew)Normal probability density function to the data.
ParameterFit[{data1, data2,...}] fits a (skew)Normal probability density function to each of the datasets. Is used in Hist."

ParameterFit2::usage = 
"ParameterFit2[data] fits two skewNormal probaility density fucntions to the data. Assuming two compartments, \
one for fat and one for muscle."

DatTot::usage = 
"DatTot[{data1, data2, ..}, name, vox] calculates the parameter table conating the volume, mean, std and 95 CI for each of the diffusion parameters."

DatTotXLS::usage = 
"DatTotXLS[{data1, data2, ..}, name, vox] is the same as DatTot, but gives the parameters as strings for easy export to excel."

GetMaskMeans::usage = 
"GetMaskMeans[dat, mask, name] calculates the mean, std, 5,50 and 95% CI form the given data for each of the given masks. 
Mask can be genereated by SplitSegmentations. name is a string that is added to the header."

FiberDensityMap::usage =
"FiberDensityMap[fiberPoins, dim, vox] generates a fiber density map for the fiberPoins which are imported by LoadFiberTracts. \
The dimensions dim should be the dimensions of the tracked datasets van vox its volxel size."

FiberLengths::usage =
"FiberLengths[fpoints,flines] calculates the fiber lenght using the output from LoadFiberTacts.
FiberLengths[{fpoints,flines}] calculates the fiber lenght using the output from LoadFiberTacts."


JoinSets::usage =
"JoinSets[{dat1,dat2,...}, over] joins dat1, dat2, ... with over slices overlap.
JoinSets[{dat1,dat2,dat3...},{over1,over2,...}] joins dat1 and dat2 with over1 slices overlap, Joins dat2 and dat3 with over2 slices overlap and so on.
JoinSets[{dat1,dat2,...},{{over,drop1,drop2},...}] joins dat1, dat2 with over slices overlap and drops drop1 slices for dat1 and drop2 from drop 2.

DOI: 10.1148/radiol.14140702."

SplitSets::usage = 
"SplitSets[data, Nsets, Nover] splits the data in Nsets with Nover slices overlap."

CorrectJoinSetMotion::usage =
"CorrectJoinSetMotion[[{dat1,dat2,...}, vox, over] motion correts multiple sets with overlap. Over is the number of slices overlap between stes. A Translation registration is performed."

DataTransformation::usage = 
"DataTransformation[data,vox,w] transforms a 3D dataset accordint to the affine transformation vector w"

InvertDataset::usage = 
"InvertDataset[data] inverts the data along the x y and z axes. In other words it is rotated aroud the origin such that (x,y,z)=(-x,-y,-z) and (0,0,0)=(0,0,0)"


Hist::usage = 
"Hist[data, range] plots a probability density histogram of the data from xmin to xmax with a fitted (skew)normal distribution. Uses ParameterFit.
Hist[data, range, label] plots a probability density histogram of the data from xmin to xmax with a fitted (skew)normal distribution and label as x-axis label.
Hist[{data1..,data2,..}, {range1,range2,..}] plots a probability density histogram of the data from xmin to xmax with a fitted (skew)normal distribution. Uses ParameterFit.
Hist[{data1,data2,..}, {range1,range2,..}, {label1,label2,..}] plots a probability density histogram of the data from xmin to xmax with a fitted (skew)normal distribution and label as x-axis label."

Hist2::usage = 
"Hist2[pars, range] plots a probability density histogram of the data over range with two fitted (skew)normal distribution. Uses ParameterFit2.
Hist2[pars, range, label] plots a probability density histogram of the data over range with two fitted (skew)normal distribution. Uses ParameterFit2."

ErrorPlot::usage = 
"ErrorPlot[data, xdata] plots a errorplot of the data where the first dim of the data is the xrange which matches the xdata list. 
ErrorPlot[data, xdata, range] similar with a given y range."


NumberTableForm::usage = 
"NumberTableForm[data] makes a right aligned table of the numbers with 3 decimal percision.
NumberTableForm[data, n] makes a right aligned table of the numbers with n decimal percision.";

MeanStd::usage = 
"MeanStd[data] calculates the mean and standard deviation and reports it as a string."

MeanRange::usage = 
"MeanRange[Range] calculates the medain (50%) and standard deviation (14% and 86%) range and reports it as a string."


FindOutliers::usage =
"FindOutliers[data] finds the outliers of a list of data."

MedCouple::usage = 
"MedCouple[data] calculates the medcouple of a list of data."


SmartMask::usage = 
"SmartMask[input] crates a smart mask of input, which is either the tensor or the tensor parameters calculated using ParameterCalc.
SmartMask[input, mask] crates a smart mask of input and used the mask as a prior selection of the input."


(* ::Subsection::Closed:: *)
(*Options*)


FitFunction::usage = 
"FitFunction is an option for ParameterFit. Options are \"Normal\" or \"SkewNormal\". Indicates which function wil be fitted."

FitOutput::usage = 
"FitOutput is an option for ParameterFit and ParameterFit2. Option can be \"Parameters\", \"Function\" or \"BestFitParameters\"."

OutputSNR::usage = 
"OutputSNR is an option for SNRMapCalc."

SmoothSNR::usage = 
"SmoothSNR is an option for SNRMapCalc."

SeedDensity::usage = 
"SeedDensity is an option for FiberDensityMap. The seedpoint spacing in mm."

MeanMethod::usage = 
"MeanMethod is an option for GetMaskMeans. The option can be  \"NormalDist\", \"SkewNormalDist\", or \"Mean\"."


ReverseData::usage =
"ReverseData is an option for JoinSets. Reverses each individual datset given as input for the JoinSets function. True by default."

ReverseSets::usage =
"ReverseSets is an option for JoinSets. Reverses the order of the datsets, False by default."

NormalizeSets::usage = 
"NormalizeSets is an option for JoinSets. True normalizes the individual stacs before joining."

NormalizeOverlap::usage = 
"NormalizeOverlap is an option for JoinSets. True removes strong signal dropoff at the end of a stack."

MotionCorrectSets::usage = 
"MotionCorrectSets is an option for JoinSets. True motion corrects the individual stacs before joining using CorrectJoinSetMotion."

JoinSetSplit::usage = 
"JoinSetSplit is an option ofr CorrectJoinSetMotion. If True RegisterDataTransformSplit is used else RegisterDataTransform is used."

PaddOverlap::usage = 
"PaddOverlap is an option of CorrectJoinSetMotion and JoinSets. it allows for extra motion in the z direction."


OutlierMethod::usage = 
"OutlierMethod is an option for FindOutliers. values can be \"IQR\", \"SIQR\" or \"aIQR\". \"IRQ\" is used for normly distributed data, \"SIQR\" or \"aIQR\" are better for skewed distributions."

OutlierOutput::usage = 
"OutlierOutput is an option for FindOutliers. If value is \"Mask\" it gives a list of 1 for data and 0 for outliers. Else the output is {data, outliers}."

OutlierIterations::usage = 
"OutlierIterations is an option for FindOutliers. Specifies how many iterations are used to find the outliers. 
Each itteration the outliers are reevaluated on the data with the previously found outliers alread rejected."

OutlierRange::usage = 
"OutlierRange is an option for FindOutliers. Specifies how many times the IQR is considred an oulier."

OutlierIncludeZero::usage = 
"OutlierIncludeZero is an option for FindOutliers. If set to True all values that are zero are ignored and considered outliers."


ColorValue::usage = 
"ColorValue is an option for Hist and ErrorPlot. Default {Black, Red}."

Scaling::usage = 
"Scaling is an option for Hist2. Scales the individual fits of the fat and muscle compartment."


Strictness::usage = 
"Strictness is an option for SmartMask value between 0 and 1. Higer values removes more data."

MaskCompartment::usage = 
"MaskCompartment is an option for SmartMask. Can be \"Muscle\" or \"Fat\"."

SmartMethod::usage = 
"SmartMethod is an option for SmartMask. This specifies how the mask is generated. Can be \"Continuous\" or \"Catagorical\""

SmartMaskOutput::usage = 
"SmartMaskOutput is an option for Smartmask. Can be set to \"mask\" to output only the mask or \"full\" to also output the probability mask."


TableMethod::usage = 
"TableMethod is an option for NumberTableForm. It specifies which number form to uses. Values can be NumberForm, ScientificForm or EngineeringForm"


(* ::Subsection::Closed:: *)
(*Error Messages*)


ParameterFit::func = "Unknow fit function: `1`. options are SkewNormal or Normal."

ParameterFit::outp = "Unknow output format: `1`. options are Parameters or Function."

JoinSets::over = "Error: The overlap must be a number or a list which gives the overlap and how many slice must be droped. Not: `1`."

Hist::size = "Length of data (`1`)must be the same as the length of the range (`2`) and labels (`3`)."

ErrorPlot::size = "Length of data (`1`)must be the same as the length of the range (`2`) and labels (`3`)."

Hist2::size = "Length of data (`1`), labels (`2`) and range (`3`) must be 5."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*SetupDataStructure*)


SetupDataStructure[dcmFolder_] := 
 Module[{folderdcm, foldernii, folderout, folders,fol, niiFolder, outFolder},
  folderdcm = Directory[] <> $PathnameSeparator <> # & /@ Select[FileNames["*", "dcm"], DirectoryQ];
  foldernii = StringReplace[#, "dcm" -> "nii"] & /@ folderdcm;
  folderout = StringReplace[#, "dcm" -> "out"] & /@ folderdcm;
  folders = Transpose[{folderdcm, foldernii, folderout}];
  
  fol = Last@FileNameSplit[dcmFolder];
  niiFolder = StringReplace[dcmFolder, fol -> "nii"];
  outFolder = StringReplace[dcmFolder, fol -> "out"];
  If[! DirectoryQ[niiFolder], CreateDirectory[niiFolder]];
  If[! DirectoryQ[outFolder], CreateDirectory[outFolder]];
  
  (*create nii files*)
  If[! DirectoryQ[#[[2]]], CreateDirectory[#[[2]]]; DcmToNii[#[[1 ;; 2]]]] & /@ folders;
  
  folders
]


(* ::Subsection:: *)
(*ParameterFit*)


(* ::Subsubsection::Closed:: *)
(*ParameterFit*)


Options[ParameterFit] = {FitFunction -> "SkewNormal", FitOutput -> "Parameters", Method -> Automatic}

SyntaxInformation[ParameterFit] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ParameterFit[dat : {_?ListQ ..}, opts : OptionsPattern[]] := ParameterFit[Flatten[#], opts] & /@ dat

ParameterFit[dat_List, OptionsPattern[]] := Module[{mod, out, met, data, mdat, sdat, fdat, sol ,par, fun},
  
  (*get option values*)
  mod = OptionValue[FitFunction];
  out = OptionValue[FitOutput];
  met = OptionValue[Method];
  
  (*prepare data*)
  data = dat;
  data=Pick[data, Unitize[data], 1];
  
  (*initialization for mean and std*)
  mdat = Mean[data];
  sdat = StandardDeviation[data];

  Off[NonlinearModelFit::"cvmit"]; Off[NonlinearModelFit::"sszero"];
  
  nodat = Length[data] <= 10;
  (*perform the fit for one compartment*)
  If[nodat,
   Print["Not Enough data in the ROI"];
   ,
   (*fit data*)
   fdat = FitData[data];
   
   Switch[mod,
    (*SkewNormal dist parameter fit*)
    "SkewNormal",
    sol = NonlinearModelFit[fdat,  PDF[SkewNormalDistribution[Mu, Sigma, Alpha], x], {{Mu, mdat}, {Sigma, sdat}, {Alpha, 0}}, x, Method -> met];
    par = sol["BestFitParameters"];
    fun = SkewNormalDistribution[Mu, Sigma, Alpha] /. par;
    ,
    (*Normal dist parameter fit*)
    "Normal",
    sol = NonlinearModelFit[fdat, {PDF[NormalDistribution[Mu, Sigma], x],Sigma>0}, {{Mu, mdat}, {Sigma, sdat}}, x];
    par = sol["BestFitParameters"];
    fun = NormalDistribution[Mu, Sigma] /. par;
    ,
    _,
    Message[ParameterFit::func, mod]]
   ];
  On[NonlinearModelFit::"cvmit"]; On[NonlinearModelFit::"sszero"];
  
  (*generate Output*)
  Switch[out,
   "Parameters",
   If[nodat,
   	{0,0},
   	{Mean[fun], StandardDeviation[fun]}
   ],
   "ParametersExtra",
   If[nodat,
   	{0.,0.,0.,0.,0.},
   	Flatten[{Mean[fun], StandardDeviation[fun], Quantile[fun, {.5, .05, .95}]}]
   ],
   "Function",
   If[nodat,0.,sol],
   "BestFitParameters",
   If[nodat,0.,par[[All,2]]],
   _, Message[ParameterFit::outp, out]
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*ParameterFit2*)


Options[ParameterFit2]={FitOutput->"BestFitParameters"}

SyntaxInformation[ParameterFit2] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ParameterFit2[dat_List, OptionsPattern[]]:=
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
	datf=Cases[Cases[Flatten[#]//N,Except[0.]],Except[1.]]&/@dat;
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
  labels = If[name==="", {"mean", "std", "Median", "5%", "95%"}, name <> " " <> # & /@ {"mean", "std", "Median", "5%", "95%"}
  ];
  out = If[Total[Flatten[#]]<=10,
  	Print["Less than 10 voxels, output will be 0."];{0.,0.,0.,0.,0.}
  	,
      fl = GetMaskData[dat, #, GetMaskOutput -> All];
      Switch[OptionValue[MeanMethod],
       "NormalDist",
       ParameterFit[fl, FitOutput -> "ParametersExtra", FitFunction -> "Normal"],
       "SkewNormalDist",
       ParameterFit[fl, FitOutput -> "ParametersExtra", FitFunction -> "SkewNormal"],
       _,
       Flatten[{Mean[fl], StandardDeviation[fl], Quantile[fl, {.5, .05, .95}]}]
       ]
  ]&/@ Transpose[mask];
      
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


(* ::Subsection:: *)
(*FindOutliers*)


(* ::Subsubsection::Closed:: *)
(*FindOutliers*)


Options[FindOutliers] = {OutlierMethod -> "IQR", OutlierOutput -> "Mask", OutlierIterations -> 1, OutlierRange -> 1.5, OutlierIncludeZero -> True}

SyntaxInformation[FindOutliers] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

FindOutliers[datai_?VectorQ, opts:OptionsPattern[]]:=FindOutliers[datai,1,opts]

FindOutliers[datai_?VectorQ, ignore_, OptionsPattern[]] :=  Block[{
	data, maxIt, diff, it, out, outI, outNew, q1, q2, q3, sc, iqr, dataQ, up, low, mc, met, incZero, output
	},
  (*make numeric*)
  data = N@datai;
  
  (*get options*)
  met = OptionValue[OutlierMethod];
  output = OptionValue[OutlierOutput];
  maxIt = OptionValue[OutlierIterations];
  sc = OptionValue[OutlierRange];
  incZero = OptionValue[OutlierIncludeZero];
  
  (*initialize*)
  diff = it = 1;
  outI = out = N@If[incZero, 0 data + 1, Unitize[data]];
  
  (*perform itterative outlier detection*)
  While[(diff != 0.) && it <= maxIt,
   
   (*get the data quantiles and iqr*)
   dataQ = Pick[data, ignore out, 1.];
   {q1, q2, q3} = Quantile[dataQ, {.25, .50, .75}];
   iqr = (q3 - q1);
   
   (*switch methods*)
   (*IQR-inter quantile range, SIQR-skewed iql, aIQR-
   adjusted iqr using medcouple for skewness*)
   {low, up} = Switch[OptionValue[OutlierMethod],
     "IQR", {q1 - sc iqr, q3 + sc iqr},
     "SIQR", {q1 - sc 2 (q2 - q1), q3 + sc 2 (q3 - q2)},
     "aIQR",
     mc = MedCouple[dataQ, q2];
     If[mc >= 0,
      {q1 - sc iqr Exp[-4 mc], q3 + sc iqr Exp[3 mc]},
      {q1 - sc iqr Exp[-3 mc], q3 + sc iqr Exp[4 mc]}
      ]
     ];
   (*make the oulier mask*)
   outNew = N[outI (If[(# < low || # > up), 0, 1] & /@ N[data])];
   (*update ouliers and itteration*)
   diff = Total[out - outNew];
   out = outNew;
   it++
   ];
  
  (*make the output*)
  If[output === "Mask",
   Round[out],
   {Pick[datai, out, 1.], Pick[datai, out, 0.]}
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*MedCouple*)


MedCouple[data_] := MedCouple[data,Median[data]];

MedCouple[data_, q2_] := Block[{xi, li, xj, lj, pi, hxixj},
  xi = Select[data, # >= q2 &];
  li = Range[pi = Length[xi]];
  xj = Select[data, # <= q2 &];
  lj = Range[Length[xj]];
  hxixj = Flatten@Table[
     If[xi[[i]] > xj[[j]],
      ((xi[[i]] - q2) - (q2 - xj[[j]]))/(xi[[i]] - xj[[j]]),
      pi - 1 - i - j
      ], {i, li}, {j, lj}];
  Median[hxixj]
  ]


(* ::Subsection:: *)
(*Number Functions*)


(* ::Subsubsection::Closed:: *)
(*NumberTableForm*)


Options[NumberTableForm] = Join[{TableMethod -> NumberForm}, Options[TableForm]];

SyntaxInformation[NumberTableForm] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

NumberTableForm[dat_, opts : OptionsPattern[]] := NumberTableForm[dat, 3, opts];

NumberTableForm[dat_, depth_, opts : OptionsPattern[]] := 
  Block[{opt, met},
   met = OptionValue[TableMethod];
   met = If[MemberQ[{NumberForm, ScientificForm, EngineeringForm}, met], met, NumberForm];
   TableForm[
   	Map[met[#, {20*depth, depth}] &, dat, {ArrayDepth[dat]}],
     FilterRules[{opts}, Options[TableForm]], TableAlignments -> Right]
   ];


(* ::Subsubsection::Closed:: *)
(*MeanStd*)

MeanStd[inp_]:=MeanStd[inp, 1]

MeanStd[inp_,n_] := Block[{dat}, 
	dat = inp /. {Mean[{}] -> Nothing, 0. -> Nothing};
	Quiet@Row[{NumberForm[Round[Mean[dat], .001], {3, n}], NumberForm[Round[StandardDeviation[dat], .001], {7, n}]}, "\[PlusMinus]"]
  ]


(* ::Subsubsection::Closed:: *)
(*MeanRange*)


MeanRange[inp_] := Block[{q1, q2, q3},
  {q1, q2, q3} = Quantile[inp /. {Mean[{}] -> Nothing, 0. -> Nothing}, {.14, .5, .86}];
  Quiet@Row[{NumberForm[Round[q2, .0001], {7, 2}], "  (", NumberForm[Round[q1, .0001], {7, 2}], " - ", NumberForm[Round[q3, .001], {7, 2}], ")"}]
  ]
  
MeanRange[inp_,quant_] := Block[{q1, q2, q3},
  {q1, q2, q3} = Quantile[inp /. {Mean[{}] -> Nothing, 0. -> Nothing}, {quant[[1]],.5,quant[[2]]}];
  Quiet@Row[{NumberForm[Round[q2, .0001], {7, 2}], "  (", NumberForm[Round[q1, .0001], {7, 2}], " - ", NumberForm[Round[q3, .001], {7, 2}], ")"}]
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

SNRMapCalc[data_?ArrayQ, noise_?ArrayQ, opts:OptionsPattern[]] := SNRMapCalc[data, noise, OptionValue[SmoothSNR], opts]
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
  sigma = ConstantArray[StandardDeviation[Cases[Flatten[noise] // N, Except[0.]]],Dimensions[signal]];
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
  div=Clip[div, {0., 100 Median[Cases[Flatten[div], Except[0.]]]}];
  snr = GaussianFilter[div, k];
  
  Switch[OptionValue[OutputSNR],
	 "Sigma", sigma,
	 "Both", {snr, sigma},
	 _, snr
	 ]
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
    Median[Cases[Flatten[density], Except[0]]]
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


(* ::Subsection:: *)
(*TransformData*)


(* ::Subsubsection::Closed:: *)
(*TransformData*)


Options[DataTransformation]={InterpolationOrder->1}

SyntaxInformation[DataTransformation]={"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};

DataTransformation[data_, vox_, wi_, OptionsPattern[]] := 
 Block[{coor, rot, coorR, interFunc, interFuncC, w},
  w = If[Length[wi]==3,Join[wi,{0,0,0,1,1,1,0,0,0}],wi];
  
  coor = GetCoordinates[data, vox];
  rot = ParametersToTransformFull[w, "Inverse"];
  coorR = ApplyRotC[coor, rot];
  interFunc = Interpolation[Transpose[{Flatten[coor, ArrayDepth[coor] - 2], Flatten[data]}], InterpolationOrder -> OptionValue[InterpolationOrder], "ExtrapolationHandler" -> {0. &, "WarningMessage" -> False}];
  interFuncC = Compile[{{coor, _Real, 1}}, interFunc[coor[[1]], coor[[2]], coor[[3]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];
  interFuncC[coorR]
  ]


(* ::Subsubsection::Closed:: *)
(*GetCoordinates*)


GetCoordinates[data_, vox_] := Block[{dim, off, coor},
   off = Dimensions[data]/2;
   coor = MapIndexed[#2 &, data, {ArrayDepth[data]}] - 0.5;
   CoordC[coor, off, vox]
   ];


(* ::Subsubsection::Closed:: *)
(*CoordC*)


   
CoordC = Compile[{{coor, _Real, 1}, {off, _Real, 1}, {vox, _Real, 1}},
    vox (coor - off),
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*ApplyRotC*)


ApplyRotC = Compile[{{coor, _Real, 1}, {rot, _Real, 2}}, 
  	(rot.Append[coor, 1])[[1 ;; 3]],
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*ParametersToTransformFull*)


ParametersToTransformFull[w_] := ParametersToTransformFull[w, "Normal"]

ParametersToTransformFull[w_, opt_] := Block[{
	tx, ty, tz, rx, ry, rz, sx, sy, sz, gx, gy, gz, 
	T, R, G, S, Rx, Ry, Rz, Gx, Gy, Gz, 
	mat, rMat, tMat}, 
	
	{rz, rx, ry, tz, tx, ty, sz, sx, sy, gz, gx, gy} = w;
	
	(*translation*)
	T = N@{{1, 0, 0, tz}, {0, 1, 0, tx}, {0, 0, 1, ty}, {0, 0, 0, 1}};
	
	(*rotation*)
	rx = rx Degree; ry = ry Degree; rz = rz Degree;
	Rz = {{1, 0, 0, 0}, {0, Cos[rz], -Sin[rz], 0}, {0, Sin[rz], Cos[rz], 0}, {0, 0, 0, 1}};
	Rx = {{Cos[rx], 0, Sin[rx], 0}, {0, 1, 0, 0}, {-Sin[rx], 0, Cos[rx], 0}, {0, 0, 0, 1}};
	Ry = {{Cos[ry], -Sin[ry], 0, 0}, {Sin[ry], Cos[ry], 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	R = N@Rx.Ry.Rz;
	
	(*scaling*)
	If[N@{sx,sy,sz}=={0.,0.,0.},{sx,sy,sz}={1.,1.,1.}];
	S = N@{{sz, 0, 0, 0}, {0, sx, 0, 0}, {0, 0, sy, 0}, {0, 0, 0, 1}};
	
	(*skew*)
	Gz = {{1, gz, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	Gx = {{1, 0, 0, 0}, {0, 1, gx, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	Gy = {{1, 0, gy, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	G = N@Gx.Gy.Gz;

	mat = T.R.G.S;
	
	Switch[opt,
		"Normal",
		mat,
		"Inverse",
		rMat = Inverse[mat[[1 ;; 3, 1 ;; 3]]];
		tMat = -rMat.mat[[1 ;; 3, 4]];
		Append[Flatten /@ Thread[{rMat, tMat}], {0, 0, 0, 1}]
		]
	]


(* ::Subsubsection::Closed:: *)
(*InvertDataset*)


InvertDataset[data_] := Module[{dep},
  dep = ArrayDepth[data];
  Switch[dep,
   4, Transpose[Inverse3Di /@ Transpose[data]],
   _, Inverse3Di[data]
   ]
  ]

Inverse3Di[data_] := Block[{out},
  out = data;
  (out = Reverse[out, #]) & /@ {1, 2, 3};
  out]


(* ::Subsection:: *)
(*Join sets*)


(* ::Subsubsection::Closed:: *)
(*JoinSets*)


Options[JoinSets]={ReverseSets->True,ReverseData->True, NormalizeOverlap->False, NormalizeSets -> True, MotionCorrectSets -> False, PaddOverlap -> 2, JoinSetSplit -> True};

SyntaxInformation[JoinSets] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

JoinSets[data_?ArrayQ,over_,opts:OptionsPattern[]]:=JoinSets[data,over,{1,1,1},opts]

JoinSets[data_?ArrayQ,over_,vox_,OptionsPattern[]]:=Block[
	{dat, overlap, motion, pad, normalize, depth, meth, target},
	
	(*get the options*)
	motion = OptionValue[MotionCorrectSets];
	pad = OptionValue[PaddOverlap];
	normalize=OptionValue[NormalizeSets];
	normover=OptionValue[NormalizeOverlap];
	depth=ArrayDepth[data];
	overlap = If[ListQ[over],First@over,over];
	
	(*normalize the data*)
	dat=If[normalize, PrintTemporary["normalizing data"]; NormalizeData/@data, data];

	(*reverse the order of the sets if needed*)
	dat=If[OptionValue[ReverseSets],Reverse[dat],dat];
	
	If[motion,
		Switch[depth,
			5,
			motion=False;
			(*define the moving data*)
			Print["motion correct is only for 3D volues"]
			,
			4,
			PrintTemporary["motion correcting data"];
			dat = CorrectJoinSetMotion[dat, vox, over, PaddOverlap->pad, JoinSetSplit->OptionValue[JoinSetSplit]];
			overlap = overlap + 2*pad;
		]
	];
	
	(*reverse the order of the slices if needed*)
	dat=N@If[OptionValue[ReverseData],Reverse[dat,2],dat];
	
	PrintTemporary["Joining data"];	
	dat = Switch[depth,
		5,Transpose[(JoinSetsi[dat[[All, All, #]],overlap,normover]) & /@ Range[Length[dat[[1, 1]]]]],
		4,JoinSetsi[dat,overlap,normover],
		_,$Failed
	];
	
	(*give output*)	
	dat = If[motion, ArrayPad[dat, Prepend[ConstantArray[{0, 0}, ArrayDepth[dat] - 1], {-pad, -pad}]],dat];
	
	Return[If[OptionValue[ReverseData],Reverse[dat],dat]]
]


JoinSetsi[data_?ArrayQ,overlap_?IntegerQ,norm_:False]:=
Module[{sets,set1,set2,step,set1over,set2over,joined,mn1,mn2},
	
	sets=Length[data];
	step=1/(overlap+1);
	
	(*perform the join*)
	For[i=1,i<sets,i++,
		If[i==1,
			set1=Drop[data[[i]],{-overlap,-1}];
			set1over=Take[data[[i]],{-overlap,-1}];
			,
			set1=Drop[joined,{-overlap,-1}];
			set1over=Take[joined,{-overlap,-1}];
			];
		set2=Drop[data[[i+1]],{1,overlap}];
		set2over=Take[data[[i+1]],{1,overlap}];
		
		If[norm,
			mn1 = MeanNoZero[Flatten[#]] & /@ set1over;
			mn2 = MeanNoZero[Flatten[#]] & /@ set2over;
			mn1 = mn1[[1]]/mn1;
			mn2 = mn2[[-1]]/mn2;
			set1over = mn1 set1over;
			set2over = mn2 set2over;
		];

		joined = Joini[{set1, set2}, {set1over, set2over}, overlap];
		];
	
	joined	

	];


JoinSetsi[data_?ArrayQ,overlap_?ListQ,OptionsPattern[]]:=
Module[{sets,set1,set2,i,step,set1over,set2over,joined,overSet,data1,data2,drop1,drop2,overl},
	
	sets=Length[data];
	
	(*perform the join*)
	For[i=1,i<sets,i++,
		overSet=overlap[[i]];
		If[i==1,
			data1=data[[i]];,
			data1=joined;
			];
		If[Length[overSet]!=3&&!IntegerQ[overSet],
			Return[Message[JoinSets::over,overSet]];
			,
			If[IntegerQ[overSet],
				step=1/(overSet+1);
				data2=data[[i+1]];
				overl=overSet;
				,
				If[Length[overSet]==3,
					overl=overSet[[1]];
					step=1/(overl+1);
					drop1=overSet[[2]];
					drop2=overSet[[3]];
					If[drop1!=0,data1=Drop[data1,{-drop1,-1}]];
					If[drop2!=0,data2=Drop[data[[i+1]],{1,drop2}];,data2=data[[i+1]];];
					]
				]
			];
		set1=Drop[data1,{-overl,-1}];
		set1over=Take[data1,{-overl,-1}];
		set2=Drop[data2,{1,overl}];
		set2over=Take[data2,{1,overl}];
		
		joined=Joini[{set1,set2},{set1over,set2over},overl];
		];
		
	joined

	];
	



(* ::Subsubsection::Closed:: *)
(*Joini*)


Joini[sets_, setover_, step_] := Module[{over,dato,unit,noZero,tot},
  (*define the overlapping voxels*)
  unit = Unitize[setover];
  noZero = Times @@ unit;
  tot = Total[noZero];
  (*prepare the data for listable compliled function*)
  noZero = TransData[noZero, "l"];
  dato = TransData[TransData[setover, "l"], "l"];
  (*merge the overlapping data*)
  over = TransData[JoinFuncC[dato, noZero, tot, step], "r"];
  (*merge the non ovelap with the overlap*)
  Chop[Join[sets[[1]], over, sets[[2]]]]
  ]


JoinFuncC = Compile[{{dat, _Real, 2}, {noZero, _Integer, 1}, {tot, _Integer, 0}, {steps, _Integer, 0}},
	Block[{ran, unit, tot1, out},
    If[tot === 0,
     (*all zeros, no overlap of signals so just the sum of signals*)
     out = Total[(0 dat + 1) dat];
     ,
     (*overlap of signals*)
     (*define the range needed*)
     tot1 = 1./(tot + 1.);
     ran = Reverse@Range[tot1, 1. - tot1, tot1];
     
     (*replace with gradient*)
     If[tot === steps,
      (*full overlap*)
      out = Total[{ran, 1. - ran} dat];
      ,
      (*partial overlap*)
      (*summ all signals*)
      unit = (0 dat + 1);
      (*replace the overlapping signals with a gradient*)
      unit[[All, Flatten[Position[noZero, 1]]]] = {ran, 1 - ran};
      (*sum the signals*)
      out = Total[unit dat]
      ]];
    (*give the output*)
    out],
    {{out, _Real, 1}}, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*CorrectJoinSetMotion*)


Options[CorrectJoinSetMotion] = {JoinSetSplit -> True, PaddOverlap -> 2}

SyntaxInformation[CorrectJoinSetMotion] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

CorrectJoinSetMotion[input_, vox_, over_, OptionsPattern[]] := Module[
	{sets, nmax, dim, d1, d2, maskd1, maskd2, samp, overp, pad, regFunc, depth},
 	
 	(*get the input*)
 	pad = OptionValue[PaddOverlap];
 	depth = ArrayDepth[input];
 	
	(*data which will be joined, make all data sets 4D*)
	sets = Switch[depth,
		5,input,
		4,Transpose[{#}]&/@input
	];
	
	(*add z padding to allow more overlap*)
	sets = ArrayPad[#, Prepend[ConstantArray[{0, 0}, ArrayDepth[#] - 1], {pad, pad}]] & /@ sets;
	
	(*set needed values*)
	nmax = Length[sets];
	overp = over + 2 pad;
	dim = Dimensions[sets[[1,All,1]]];
	
	(*define the registration function*)
	regFunc = If[OptionValue[JoinSetSplit],RegisterDataTransformSplit,RegisterDataTransform];
	
	i=0;
	PrintTemporary[Dynamic[i]];
	
	(*perform the motion correction*)
	Table[
		i=n;
		(*get the seconds overlap stac*)
		d1 = sets[[n, ;; overp,1]];
		maskd1 = Dilation[#, 10] & /@ Mask[d1];
		(*pad to allow motion*)
		d1 = PadLeft[d1, dim];
		maskd1 = PadLeft[maskd1, dim];
		
		(*get the seconds overlap stac*)
		d2 = sets[[n + 1, -overp ;;,1]];
		maskd2 = Dilation[#, 10] & /@ Mask[d2];
		(*pad to allow motion*)
		d2 = PadLeft[d2, dim];
		maskd2 = PadLeft[maskd2, dim];
		
		maskd1 = maskd2 = Dilation[maskd1 maskd2, 1];
		
		(*get the number of samples for the registration*)
		samp = Round[((Total@Flatten@maskd1)+(Total@Flatten@maskd2))/20];
		
		(*perform the registration*)
		sets[[n + 1]] = Last@regFunc[{d1, maskd1, vox}, {d2, maskd2, vox}, {sets[[n + 1]], vox},
				MethodReg -> "translation", Iterations -> 100, NumberSamples -> samp, PrintTempDirectory -> False, InterpolationOrderReg -> 0];
		
		, {n, 1, nmax - 1}
	];
	
	(*output the data, make the 3D data 3D again*)
	Switch[depth,
		5,sets,
		4,sets[[All,All,1]]
		]
  
  ]


(* ::Subsection::Closed:: *)
(*SplitSets*)


Options[SplitSets] = {ReverseSets -> False, ReverseData -> True, PaddOverlap->0};

SyntaxInformation[SplitSets] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

SplitSets[data_, sets_, overlap_, OptionsPattern[]] := Module[{lengthSet, sels, start, end, dat, over, pad},
  
  dat = If[OptionValue[ReverseData], Reverse[data], data];
  
  pad = OptionValue[PaddOverlap];  
  dat = ArrayPad[dat,{{pad,pad},{0,0},{0,0}}];
  
  over=overlap+2pad;
  
  lengthSet = Round[(Length[dat] + (sets - 1)*over)/sets];
  sels = Table[
    start = (i lengthSet + 1) - i over;
    end = start + lengthSet - 1;
    Range[start, end]
    , {i, 0, sets - 1}];
  
  dat = (dat[[#]] & /@ sels);
  
  dat = If[OptionValue[ReverseData], Reverse[dat, 2], dat];
  dat = If[OptionValue[ReverseSets], Reverse[dat], dat];
  
  dat
  ]


(* ::Subsection::Closed:: *)
(*Hist*)


labStyle=Directive[Bold,FontFamily->"Helvetica",14,Black];


Options[Hist] = {ColorValue -> {{Black,White}, Red, Green, Blue}, Method -> "SkewNormal", PlotLabel -> "", AxesLabel -> "", ImageSize -> 300}

SyntaxInformation[Hist] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

Hist[dat_, ops : OptionsPattern[]] := Hist[dat, 0, ops]

Hist[dat_, range_, OptionsPattern[]] := Module[{sol, line, hist, x, colbar, coledge, color2, color3, color4, data, r1, r2, sdr, m, s, title, label,mn,std ,fdat},
  
  {{colbar,coledge}, color2, color3, color4} = OptionValue[ColorValue];
  
  title = OptionValue[PlotLabel];
  title = If[title === "", None, title];
  label = OptionValue[AxesLabel];
  label = If[label === "", None, label];
  
  fdat = N@Flatten[dat];
  data = Pick[fdat, Unitize[fdat], 1];
  
  {r1, r2} = If[range===0,
  	Quantile[data, {0.01, .99}],
  	If[IntegerQ[range],
  		sdr = range;
  		{m, s} = {Mean[data], StandardDeviation[data]};
  		{(m - sdr s), (m + sdr s)},
  		range
    ]];

  Switch[OptionValue[Method],
   "None",
   line = Graphics[{}],
   "Normal",
   sol = ParameterFit[data, FitOutput -> "Function", FitFunction -> "Normal"];
   line = Plot[sol[x], {x, r1, r2}, PlotStyle -> {Thick, color2}, PlotRange -> Full],
   "SkewNormal",
   sol = ParameterFit[data, FitOutput -> "Function", FitFunction -> "SkewNormal"];
   line = Plot[sol[x], {x, r1, r2}, PlotStyle -> {Thick, color2}, PlotRange -> Full],
   "Both",
   sol = {
     ParameterFit[data, FitOutput -> "Function", FitFunction -> "SkewNormal"], 
     ParameterFit[data, FitOutput -> "Function", FitFunction -> "Normal"]
     };
   line = Plot[{sol[[1]][x], sol[[2]][x]}, {x, r1, r2}, PlotStyle -> {Directive[Thick, color2], Directive[Thick, color3]}, PlotRange -> Full],
   "All",
   sol = {
     ParameterFit[data, FitOutput -> "Function", FitFunction -> "SkewNormal"], 
     ParameterFit[data, FitOutput -> "Function", FitFunction -> "Normal"]
     };
     mn=Mean[data];
     std=StandardDeviation[data];
     line = Plot[{sol[[1]][x], sol[[2]][x],PDF[NormalDistribution[mn, std], x]}, {x, r1, r2}, 
     	PlotStyle -> {Directive[Thick, color2], Directive[Thick, color3], Directive[Thick, color4]}, PlotRange -> Full]
   ];
  
  hist = Histogram[
    Select[data, (r1 < # < r2) &], {r1, r2, (r2 - r1)/30}, 
    "ProbabilityDensity",
    PerformanceGoal -> "Speed", PlotRange -> {{r1, r2}, All},
    PlotLabel -> title, LabelStyle -> labStyle, Axes -> False, 
    FrameStyle -> Thick,
    FrameLabel -> {label, "Probability Density"}, 
    Frame -> {True, True, False, False}, 
    ChartBaseStyle -> EdgeForm[coledge], ChartStyle -> colbar];
  
  Show[hist, line, ImageSize -> OptionValue[ImageSize]]
  ]


(* ::Subsection::Closed:: *)
(*Hist2*)


Options[Hist2]={Scaling->False}

SyntaxInformation[Hist2] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

Hist2[dat_?ArrayQ,range:{{_,_}..},label:{_String..},OptionsPattern[]]:=
Module[{data,line,line1,line2,hist,x,f,Omega1,Omega2,Xi1,Xi2,Alpha1,Alpha2,r1,r2,sol},
	
	If[!(Length[range]==Length[dat]==Length[label]==5),Return[Message[Hist2::size,Length[dat],Length[range],Length[label]]]];
	data=DeleteCases[Flatten[#]//N,0.]&/@dat;
	sol=ParameterFit2[data];
	Map[
		(
		{r1,r2}=range[[#]];
		
		{f,Omega1,Omega2,Xi1,Xi2,Alpha1,Alpha2}=sol[[#]];
			
			line1=Plot[If[OptionValue[Scaling],(1-f),1]*SkewNorm[x,Omega2,Xi2,Alpha2],
				{x,r1,r2},PlotStyle->{Thick,Red},PlotRange->Full];
			
			line2=Plot[If[OptionValue[Scaling],(f),1]*SkewNorm[x,Omega1,Xi1,Alpha1],
				{x,r1,r2},PlotStyle->{Thick,Blue},PlotRange->Full];
			
			line=Plot[f SkewNorm[x,Omega1,Xi1,Alpha1]+(1-f)SkewNorm[x,Omega2,Xi2,Alpha2],
				{x,r1,r2},PlotStyle->{Thick,Green},PlotRange->Full];
			
			hist=Histogram[
				Select[data[[#]],(r1<#<r2)&],{Range[r1,r2,(r2-r1)/50]},"ProbabilityDensity",
				PerformanceGoal->"Speed",
				PlotRange->{{r1,r2},All},LabelStyle->labStyle,Axes->False,FrameStyle->Thick,
				FrameLabel->{label[[#]],"Probability Density"},Frame->{True,True,False,False},
				ChartBaseStyle->EdgeForm[White],ChartStyle->Black
				];
				
			Show[hist,line1,line2,line])&,Range[Length[range]]
		]
	]


Phi[x_]:=1/(E^(x^2/2)*Sqrt[2*Pi]);
CapitalPhi[x_]:=.5(1+Erf[(x)/Sqrt[2]]);
Delta[a_]:=a/Sqrt[1+a^2];
Mn[w_,e_,a_]:=e+w Delta[a] Sqrt[2/Pi];
Var[w_,a_]:=w^2(1-(2Delta[a]^2/Pi));
SkewNorm[x_,Omega_,Xi_,Alpha_]:=(2/Omega)Phi[(x-Xi)/Omega]CapitalPhi[Alpha (x-Xi)/Omega];


(* ::Subsection::Closed:: *)
(*ErrorPlot*)


Options[ErrorPlot] = {ColorValue -> {Black, Red}, PlotLabel -> "", AxesLabel -> "", ImageSize -> 300, Method->"median"}

SyntaxInformation[ErrorPlot] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

ErrorPlot[dat_, xdat_, ops : OptionsPattern[]] := ErrorPlot[dat, xdat, 0, ops]

ErrorPlot[dat_, xdat_, range_, OptionsPattern[]] := 
 Block[{color1, color2, title, label, fdat, mn, sd, er1, er2, sdr, m, s, plr},
  {color1, color2} = OptionValue[ColorValue];
  
  title = OptionValue[PlotLabel];
  title = If[title === "", None, title];
  label = OptionValue[AxesLabel];
  label = If[label === "", None, label];
  
  fdat = DeleteCases[N@Flatten[#], 0.]&/@dat;
  
  Switch[OptionValue[Method],
  	"mean",
  	{mn,sd} = Transpose[{Mean[#],StandardDeviation[#]}&/@fdat];
  	{er1, er2} = {mn-sd,mn+sd};
  	,"median",
  	{er1, mn, er2} = Transpose[Quantile[#,{.25,.5,.75}]&/@fdat];
  	,"fit",
  	{mn, sd} = Transpose[ParameterFit[fdat]];
  	{er1, er2} = {mn - sd, mn + sd} /. 0 -> Null;
  	,_,
  	{er1, mn, er2} = Transpose[Quantile[#,{.25,.5,.75}]&/@fdat];
  ];
    
  plr = If[range === 0,
    Quantile[Flatten[fdat], {0.01, .99}],
    If[IntegerQ[range],
     sdr = range;
     {m, s} = {Mean[Flatten[fdat]], StandardDeviation[Flatten[fdat]]};
     {(m - sdr s), (m + sdr s)},
     range]
    ];
  
  ListPlot[Transpose[{xdat, #}] & /@ {mn, er1, er2},
   PlotRange -> {MinMax[xdat], plr},
   PlotLabel -> title, FrameLabel -> label,
   Axes -> False, Frame -> {True, True, False, False},
   FrameStyle -> Thick, 
   PlotStyle -> {{color1, Thick}, {Dashed, Thick, color2}, {Dashed, 
      Thick, color2}},
   Joined -> {True, True, True}, Filling -> {2 -> {3}}, 
   FillingStyle -> Directive[Opacity[0.2], color2], LabelStyle -> labStyle, ImageSize->OptionValue[ImageSize]
   ]
  ]


(* ::Subsection::Closed:: *)
(*SmartMask*)


Options[SmartMask]={Strictness->0.50, MaskCompartment->"Muscle", SmartMethod->"Continuous", SmartMaskOutput->"mask"};

SyntaxInformation[SmartMask] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SmartMask[input_,ops:OptionsPattern[]]:=SmartMask[input, 0, ops]

SmartMask[input_,maski_,OptionsPattern[]]:=Module[{
	sol,func,range,map,mask,pmask,pars
	},
	
	(*get the parameter from the tensor else use input parameters*)
	pars = If[Length[input]==6,
		PrintTemporary["Caculating Parameters"];
		ParameterCalc[input],
		input];
	
	pmask = Mask[pars[[4]] , {0.1, 4}];
	
	(*find the histogram solution*)
	sol=If[maski===0,
		Switch[
			OptionValue[MaskCompartment],
			"Muscle",
			ParameterFit2[pars][[All,{3,5,7}]],
			"Fat",
			ParameterFit2[pars][[All,{2,4,6}]]
			]
			,
			ParameterFit[GetMaskData[#,maski pmask]&/@pars,FitOutput->"BestFitParameters"]
		];
	
	
	Switch[OptionValue[SmartMethod],
		"Catagorical",
		range = (func = SkewNormalDistribution[#2[[1]], #2[[2]], #2[[3]]]; Quantile[func, {.02, .98}]) & /@ sol;
		range = Clip[range, {0, Infinity}, {10^-3, 0.}];
		map = Total[MapThread[Mask[#1,#2]&,{pars,range}]]/5;
		mask = pmask * Mask[TotalVariationFilter[map,.15],{OptionValue[Strictness]}];
		,
		"Continuous",
		map = MapThread[PDF[SkewNormalDistribution[#2[[1]], #2[[2]], #2[[3]]], #1] &, {pars, sol}];
		map = Total[{1, 1, 1, 1, 2}*(#/Max[#] & /@ map)]/6;
		mask = pmask * Mask[TotalVariationFilter[map, .25], {OptionValue[Strictness]}];
		];
		
		If[OptionValue[SmartMaskOutput]==="mask",mask,{mask,map}]
	]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
