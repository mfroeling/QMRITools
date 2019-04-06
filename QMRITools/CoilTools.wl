(* ::Package:: *)

(* ::Title:: *)
(*QMRITools CoilTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`CoilTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`CoilTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


LoadCoilSetup::usage = 
"LoadCoilSetup[file] load a very specific type of coil experiment, a dynmic scan with a setup of which the second dynamic is a noise measurement.
The input file is the Nii file that conatins the individualy reconstruted coil images and the noise data.
Internaly it uses CoilSNRCalc and SumOfSquares. 

Output is the coil data with coil noise data and snrmap based on the SumOfSquares addition, the SOS reconstruction and the SOS weights.
{dataC, noiseC, sosC, snrC, sigmapC, weights, vox}."

LoadCoilTarget::usage =
"LoadCoilTarget[file] loads a very specific typ of experiment, a dynamic scan with with the second dynmaic is a noise measuremnt.
The input file is the Nii file that conatins the scanner reconstruction and the noise data.
Internaly it uses SNRMapCalc, 

Output is the reconstructed data with noise data and snrMap {dataC, noiseC, sosC, snrC, sigmapC, weights, vox}."

CoilSNRCalc::usage = 
"CoilSNRCalc[coils, noise] calculates the sensitivity weighted snr of multiple coil elements using magnitude signal and noise.

Output is {data, noise, sos, snr, sigmap, weights}."

MakeWeightMask::usage =
"MakeWeightMask[weights] creates a mask of homogeneous regions of weightmaps removing the noise."

FindCoilPosition::usage = 
"FindCoilPosition[weights] finds the coil posision by locating the highest intensity location in the coil weight map, which can be obtianed by LoadCoilSetup or SumOfSquares.
Internally it uses MakeWeightMask to remove the noise of the weightmasks.
FindCoilPosition[weights, mask] limits the search region to the provided mask."

NoiseCorrelation::usage = 
"NoiseCorrelation[noise] calculates the noise correlation matrix, noise is {nrCoils, noise Samples}."

NoiseCovariance::usage = 
"NoiseCovariance[noise] calculates the noise covariance matrix, noise is {nrCoils, noise Samples}."

MakeNoisePlots::usage =
"MakeNoisePlots[noise] returns a grid of plots of the noise per channel
MakeNoisePlots[noise, {met, prt}] met can be \"Grid\" with prt a number or Automatic. Else all plots will be returend as a list of plots.
MakeNoisePlots[noise, {met, prt}, sub] sub defines how much the noise is subsampled, default is 40 (every 40th sample is used in plot)."

MakeCoilLayout::usage = 
"MakeCoilLayout[{name, size, number}] makes a coil grid with label name, partioned in size rows and with label number.
MakeCoilLayout[{name, size, number}, val] makes a coil grid with label name, partioned in size rows and with label the val at location number.
MakeCoilLayout[{coils..}] same but for multile coils grids. Each coil grid is defined as {name, size, number}.
MakeCoilLayout[{coils..}, val] savem but for multiple coil grids."


(* ::Subsection:: *)
(*Options*)


OutputCoilSurface::usage = 
"OutputCoilSurface is an option for FindCoilPosition. If set true it will also output a SurfacePlot of the coil location volume."

CoilSurfaceVoxelSize::usage =
"CoilSurfaceVoxelSize is an option for FindCoilPosition. Specifies the voxel size used for OutputCoilSurface."

CoilArrayPlot::usage = 
"CoilArrayPlot is an option for MakeCoilLayout. If True and values are provided it makes an arrayplot of the coil layouts"

(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*LoadCoilSetup*)


SyntaxInformation[LoadCoilSetup] = {"ArgumentsPattern" -> {_}};

 (*load a sereis of linked coils measurements.*)
LoadCoilSetup[fileC_?StringQ]:=Block[{dataC, vox, len, noiseC, sosC, snrC, sigmapC, weights, mn},
	(*import the data*)
	{dataC, vox} = ImportNii[fileC];
	(*split the noise dan data*)
	len = Round[Length[dataC[[1]]]/2];
	noiseC = Transpose@dataC[[All, (len + 1) ;;]];
	dataC = Transpose@dataC[[All, 1 ;; len]];
	
	mn =  MeanNoZero@Flatten@N@noiseC;
	dataC = 10. dataC/mn;
	noiseC = 10. noiseC/mn;
	
	(*calculate snr*)
	{dataC, noiseC, sosC, snrC, sigmapC, weights} = CoilSNRCalc[dataC, noiseC];
	(*output*)
	{dataC, noiseC, sosC, snrC, sigmapC, weights, vox}
]


(* ::Subsection:: *)
(*LoadCoilTarget*)


SyntaxInformation[LoadCoilTarget] = {"ArgumentsPattern" -> {_}};

(*load a series of body coil images*)
LoadCoilTarget[file_?StringQ] := Block[{
	dataR, vox, noiseR, mn, mask, snr, sigma
	},
	(*import target1*)
	{dataR, vox} = ImportNii[file];
	(*split the noise dan data*)
	noiseR = dataR[[All, 2]];
	dataR = dataR[[All, 1]];
	(*normalize the data*)
	mn =  MeanNoZero@Flatten@N@noiseR;
	dataR = 10. dataR/mn;
	noiseR = 10. noiseR/mn;
	(*calculate the SNR*)
	{snr, sigma} = SNRMapCalc[dataR, noiseR, OutputSNR -> "Both"];
	(*output*)
	{dataR, noiseR, snr, sigma, vox}
]


(* ::Subsection::Closed:: *)
(*CoilSNRCalc*)


SyntaxInformation[CoilSNRCalc] = {"ArgumentsPattern" -> {_, _}};

(*calculate the combineds snr form multiple coils images*)
CoilSNRCalc[coils_, noise_] := Block[{
	mn, sigmap, coilsN, noiseN, sumSquares, weights, snr
	},
	(*get mean noise*)
	mn = MeanNoZero@Flatten@N@noise;
	(*normalize all coils to constant noise level*)
	coilsN = 10. coils/mn;
	noiseN = 10. noise/mn;
	(*calcualte the sum of squares signal*)
	{sumSquares, weights} = SumOfSquares[coilsN];
	(*calculated the weitghted noise addition*)
	{snr, sigmap} = WeigthedSNR[coilsN, noiseN, weights];
	(*output*)
	{coilsN, noiseN, sumSquares, snr, sigmap, weights}
]


WeigthedSNR[signal_, noise_, weights_] := Block[{sigmap, sigtot, snr},
	sigtot = Total[signal weights];
	sigmap = Sqrt[Total[weights^2 (Sqrt[2./Pi] GaussianFilter[noise, 3])^2]];
	snr = DevideNoZero[sigtot, sigmap];
	{snr, sigmap}
  ]


(* ::Subsection::Closed:: *)
(*FindCoilPosision*)


SyntaxInformation[MakeWeightMask] = {"ArgumentsPattern" -> {_}};

MakeWeightMask[weights_?ArrayQ] := Block[{out, back},
	out = Total[StdFilter[#, 2] & /@ weights];
	out = out/Max[out];
	back = 1 - Mask[Total@weights, {0, 0.0001}];
	out = back Mask[out, {0, 0.4}];
	out=ImageData[SelectComponents[Image3D[out], "Count", -1]];
	{out, ArrayPad[Closing[ArrayPad[out,30],15],-30]}
]


(* ::Subsection::Closed:: *)
(*FindCoilPosision*)


Options[FindCoilPosition] = {OutputCoilSurface->False, CoilSurfaceVoxelSize->{1,1,1}}

SyntaxInformation[FindCoilPosition]={"ArgumentsPattern" -> {_, _., OptionsPattern[]}}

FindCoilPosition[weights_, opts:OptionsPattern[]] := FindCoilPosition[weights, 1, opts]

FindCoilPosition[weights_, mask_, OptionsPattern[]] := Block[{
	dim, ss, max, mm, pdat, comp, pl, pts, maski, weightsM, maska, ncoils
	},
	
	(*mask the noise in the weights*)
	PrintTemporary["Masking weights"];
	{maski,maska} = MakeWeightMask[weights];
	
	(*dimensions*)
	dim = Reverse[Dimensions[weights[[1]]]];
	
	(*smooth data*)
	weightsM = N[(maski - Erosion[maski, 4]) GaussianFilter[#, 4]] &/@ weights;
	ncoils=Length[weightsM];
	
	(*find the coil posisions*)
	
	{pts, pl} = Monitor[Transpose@Table[
		(*get coil signal mask*)
		pdat = mask weightsM[[i]];
		max = Quantile[DeleteCases[Flatten[pdat], 0.], .995];
		comp = SelectComponents[Image3D[ArrayPad[Mask[pdat, max], 1]], "Count", -1];
		pdat = ArrayPad[GaussianFilter[ImageData[comp], 2], -1];
	
		(*get coil coordinate*)
		comp = ComponentMeasurements[comp, "Centroid"];
		pts = Round[If[comp === {}, {-10, -10, -10}, comp[[1, 2]]]];
		pts = Abs[pts - {0, 1, 1} (dim + 1)];
		
		If[OptionValue[OutputCoilSurface],
			pl = PlotContour[pdat, OptionValue[CoilSurfaceVoxelSize], ContourStyle -> {Red, 1}],
			pl = 1;
		];
		
		{pts, pl}
		,{i,1,ncoils}], ProgressIndicator[i,{1,ncoils}]];
	
	If[OptionValue[OutputCoilSurface],{pts, pl},pts]
]


(* ::Subsection::Closed:: *)
(*NoiseCorrelation*)


SyntaxInformation[NoiseCorrelation]={"ArgumentsPattern" -> {_}}

NoiseCorrelation[noise_] := Block[{nrCoils, corr},
  nrCoils = Length[noise];
  corr = Table[Which[i == j, .5, i < j, 0., True, Correlation[noise[[i]], noise[[j]]]], {i, 1, nrCoils}, {j, 1, nrCoils}];
  corr + Transpose[corr]
  ]


(* ::Subsection::Closed:: *)
(*NoiseCovariance*)


SyntaxInformation[NoiseCovariance]={"ArgumentsPattern" -> {_}}

NoiseCovariance[noise_] := Block[{nrCoils, cova},
  nrCoils = Length[noise];
  cova = Table[Which[i == j, .5 Covariance[noise[[i]], noise[[j]]], i < j, 0., True, Covariance[noise[[i]], noise[[j]]]], {i, 1, nrCoils}, {j, 1, nrCoils}];
  cova + Transpose[cova]
  ]


(* ::Subsection::Closed:: *)
(*MakeNoisePlots*)


SyntaxInformation[MakeNoisePlots]={"ArgumentsPattern" -> {_, _. ,_.}}

MakeNoisePlots[noise_] := MakeNoisePlots[noise, {"Grid", Automatic}, 40];

MakeNoisePlots[noise_, met_] := MakeNoisePlots[noise, met, 40];

MakeNoisePlots[noise_, {output_, met_}, sub_] := Block[{nrCoils, mean, max, noiseC, col, part, plots, ran},
   nrCoils = Length[noise];
   mean = Mean@Abs@Flatten@noise;
   
   plots = Table[
     noiseC = noise[[i, 1 ;; ;; sub]];
     max = Max@Abs@Flatten@noiseC;
     
     (*define plot range and color*)
     Which[
      max < mean, col = Red; ran = 5 mean;,
      max > 10 mean, col = Green; ran = 50 mean;,
      True, col = Blue; ran = 5 mean;
      ];
     (*make plot*)
     ListLinePlot[Re@noiseC, ImageSize -> 200, PlotStyle -> col, 
      PlotRange -> {-ran, ran}, PerformanceGoal -> "Speed",
      PlotLabel -> Style["Coil Nr: " <> ToString[i], Black, 12], 
      AspectRatio -> 1, Ticks -> None]
     , {i, 1, nrCoils}];
   
   If[output == "Grid",
    part = 
     If[met =!= Automatic && NumberQ[met], Round[met], 
      Ceiling@Sqrt[nrCoils], met];
    GraphicsGrid[Partition[plots, part, part, 1, Graphics[]], 
     ImageSize -> 1200]
    ,
    plots
    ]
   ];


(* ::Subsection::Closed:: *)
(*MakeCoilLayout*)


(* ::Subsubsection::Closed:: *)
(*MakeCoilLayout*)


Options[MakeCoilLayout] = {PlotRange -> Automatic, ColorFunction -> "SunsetColors", ImageSize -> 100, CoilArrayPlot -> False};

SyntaxInformation[MakeNoisePlots]={"ArgumentsPattern" -> {_, _. , OptionsPattern[]}}

MakeCoilLayout[{name_?StringQ, size_?NumberQ, number_?ListQ}, opts : OptionsPattern[]] := MakeCoilLayouti[{name, size, number}, 0, opts]

MakeCoilLayout[{name_?StringQ, size_?NumberQ, number_?ListQ}, val_, opts : OptionsPattern[]] := MakeCoilLayouti[{name, size, number}, val, opts]

MakeCoilLayout[coils_?ListQ, val_?ListQ, opts : OptionsPattern[]] := MakeCoilLayouti[#, val, opts] & /@ coils

MakeCoilLayout[coils_, opts : OptionsPattern[]] := MakeCoilLayouti[#, 0, opts] & /@ coils


(* ::Subsubsection::Closed:: *)
(*MakeCoilLayouti*)


Options[MakeCoilLayouti] = Options[MakeCoilLayout];

MakeCoilLayouti[{name_?StringQ, size_?NumberQ, number_?ListQ}, val_, OptionsPattern[]] := Block[{grid, out},
  If[val === 0,
   grid = Transpose[Partition[number, size]];
   out = Column[{name, Grid[grid, Frame -> All]}, Alignment -> Center];
   ,
   grid = Transpose[Partition[val[[number]], size]] /. _?Negative -> Green;
   If[OptionValue[CoilArrayPlot],
    out = ArrayPlot[grid, ImageSize -> OptionValue[ImageSize],
       PlotRange -> OptionValue[PlotRange], 
       PlotLabel -> Style[name, Black, Bold, 12],
       ColorFunction -> OptionValue[ColorFunction], Frame -> False];
    ,
    out = Column[{name, Grid[grid, Frame -> All]}, Alignment -> Center];
    ]
   ];
  out
  ]


(* ::Section:: *)
(*End Package*)


End[](* End Private Context *)

EndPackage[]
