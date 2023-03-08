(* ::Package:: *)

(* ::Title:: *)
(*QMRITools CardiacTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`CardiacTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`CardiacTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


HelixAngleCalc::usage = 
"HelixAngleCalc[eigenvectors, mask, vox] calculates the helix angle matrix of cardiac data using only a left ventricle mask.
HelixAngleCalc[eigenvectors, mask, maskp, vox] calculates the helix angle matrix of cardiac data using only a left ventricle mask, and a maskp for visualization.
HelixAngleCalc[eigenvectors, mask, centerpoint, vec, inout, vox]  calculates the helix angle matrix of cardiac data using only a left ventricle mask.
HelixAngleCalc[eigenvectors, mask, maskp, centerpoint, vec, inout, vox] calculates the helix angle matrix of cardiac data using a left vantricle mask and a maskp for visualization.

eigenvectors are the tensor eigenvectors calculated with EigenvecCalc.
mask is a mask of the left ventricle.
maskp is a mask used for visualization.
vox is the voxels size, {slice, x, y}.

The folowing values are calculated automaticlay Using CentralAxes but can also be provided as an input.
centerpoint is the center of each slice calculated with CentralAxes.
inout is the inner and outer radius calcualted with CentralAxes.
vec is the vector describin the central axes of the heart, calculated with CentralAxes.

Output is het fiber angle matrix FAM = {9, slice, x, y} or {FAM, plot}.
The angles are in degrees. 

HelixAngleCalc[] is based on DOI: 10.1186/1532-429X-17-S1-P15."

CardiacCoordinateSystem::usage = 
"CardiacCoordinateSystem[mask, vox] creates the cardiac coordinate system within the mask and is used in HelixAngleCalc. 
Output is a set of vectors {radvecn, norvecc, cirvec}, being the radial, normal and circular axes of each voxel respectivley.
If the option showPlot is true the output is {{radvecn, norvecc, cirvec}, plots}."

CentralAxes::usage = 
"CentralAxes[mask, vox] calculates the center of the lumen from a mask of the left ventricle. vox is the voxels size, {slice, x, y}.
CentralAxes[mask, maskp, vox] allows for fancy visualization of the other structures using maskp.

Output is {centerpoints, normalvecs, inout} or {centerpoints, normalvecs, inout, fit}."

CalculateWallMap::usage = 
"CalculateWallMap[mask,vox] calculates the wall distance map and the wall derivative.

Output is {wallmap, wallDerivative}."


GetMaskSegmentPoints::usage = 
"GetMaskSegmentPoints[mask] get the attacthment points from a cardiac segmentation where the heart has label 1, and the attachment points have label 2 and 3.

Output is {maks, pts} where now in mask the points are removed."

GetSegmentSlices::usage = 
"GetSegmentSlices[mask] based on the mask it gives back the slice numbers of the apex, apical, mid-ventircal, and basal slices.
GetSegmentSlices[points] does the same but then based on the points obtained form GetMaskSegmentPoints."

SegmentsPerSlice::usage = 
"SegmentsPerSlice[points] gives the number of segments per slice where the slice distribution is determined by GetSegmentSlices.
SegmentsPerSlice[slices, points] does the same but the slices are given manually."

MaskToLines::usage = 
"MaskToLines[mask, vox] calculates lines perpendicular to the heart wall per slice within the mask. Internally it uses CalculateWallMap and CentralAxes to obtain the cardiac geometry from mask. 
MaskToLines[mask, wall, cent] where mask is the first output of CalculateWallMap and cent is the first output of CentralAxes."

LinesToSegmentIndex::usage = 
"LinesToSegmentIndex[lines, points, segments] finds the lines indeces correspoinding to the points and the segments borders. Additionally it finds all the lines indeces for all lines within each segment.
The lines are comupted by MaskToLines, the points are cumputed by GetMaskSegmentPoints, and the segments is the output of SegmentsPerSlices.

Output {pointIndex, segmentIndex, lineIndex}."

GetSegmentLines::usage = 
"GetSegmentLines[lines, lineIndex, segments] groups the transmural lines per segment."

SegmentLinesToMask::usage = 
"SegmentLinesToMask[mask, segLines] cuts the mask based one the tranmural lines per segments which can be obtained by GetGesmentLines."

MakeLineImage::usage = 
"MakeLineImage[back, segLines, pts] makes an image of the cardiac segmentation lines."

MakeMaskImage::usage = 
"MakeMaskImage[back, mask] maskes an image of the cardiac segmentation mask."

CardiacSegment::usage = 
"CardiacSegment[mask, vox, pts] segments the mask in the AHA17 segmenation using pts to indicate the attachemnts.
CardiacSegment[mask, back, vox, pts] the same where back is used for image generation.
CardiacSegment[mask, vox, pts, seg] does the same but seg can be an alternate segmentation to the AHA17.
CardiacSegment[mask, back, vox, pts, seg] does the same but seg can be an alternate segmentation to the AHA17 where back is used for image generation."

CardiacSegmentGUI::usage = 
"CardiacSegmentGUI[data, mask, vox] allows to segment the heart in 1, 4, 6 or AHA-17 segements for each slice 360 radial samples are generated.

data is a background image on which all overlays are projected. 
mask is the mask of the left ventricle (same as used for CentralAxes) and defines the area in which the data is sampled.
off is the centerpoints generated by CentralAxes.

Output is {points, slices , {rev, flip}}."


MaskHelix::usage = 
"MaskHelix[helix, mask] masks helix angle data, sets the background to -100 and allows for Median filter of the helix mask.
helix can be a singel map or the FAM.

Output is the masked helix angle data."


PlotSegments::usage = 
"PlotSegments[mask, data, segang] shows how the heart wil be sampled by RadialSample. 

mask is a mask the left ventricle that was used in the CardiacSegment.
function and the segang is the output of the cardaic SegmentFunction.

Output is a plot window."

PlotSegmentMask::usage = 
"PlotSegmentMask[mask, segmask, vox] plots the mask segements created by CardiacSegment.

mask is a mask the left ventricle that was used in the CardiacSegment.
segmask is the output of CardiacSegemnt.
vox is the voxels size, {slice, x, y}.

Output is a plot window."


RadialSample::usage = 
"RadialSample[mask, data, segang] radialy samples the provided parametermap data. 

The mask should be a mask of the left ventricle that was used in the CardiacSegment.
segang is the output of the cardaic SegmentFunction.

Output is {points, vals} which are orderd as indicated by the user."

TransmuralPlot::usage = 
"TransmuralPlot[data] plots transmural profiles of the data which are created by RadialSample.

data can be a single profile or a list of profiles. In the second case the mean and standardeviations are plotted.

Output is a plot of the transmural profile."


BullseyePlot::usage = 
"BullseyePlot[data, segmask] generates a AHA-17 segement bullseye plot. 
BullseyePlot[list] generates a AHA-17 segement bullseye plot of the lists (which needs to have 17 values) provide.

data is a 3D volume used for the plot. 
segmask is the AHA-17 segmentation resulting form the CardiacSegment function when AHA17 is selected.

Output is a bullseye plot or a plotwindow, depending on the Method which can be \"Dynamic\" else it will be static.

BullseyePlot[] is based on DOI: 10.1161/hc0402.102975."


ExcludeSlices::usage = 
"ExcludeSlices[data] excludes slices that do not look like the others based on various distance measures.

Output is an array with 1 or 0 with the dimensiosn {slices, diff dirs}."


MakeECVBloodMask::usage = 
"MakeECVBloodMask[T1pre, T1post] makes a bloodpool mask based on the T1pre and T1post images. It assumes that the hart is cropped with the blood in the center.

The T1pre and T1post maps are assuemed to be in ms."


ECVCalc::usage = 
"ECVCalc[T1pre, T1post, hema] calculates the ECVmap using MakeECVBloodMask.
ECVCalc[T1pre, T1post, bloodMask, hema] calculates the ECVmap using bloodMask.

The T1pre and T1post maps are assuemed to be in ms."


CreateHeart::usage = 
"CreateHeart[] creates a simulated left ventricle shape.
CreateHeart[pars] creates a simulated left ventricle shape with predifined parameters pars.

Output is the heart shape, the voxel size and the parameters needed to generate the heart, {mask, vox, pars}."


(* ::Subsection:: *)
(*Options*)


LCMMethod::usage = 
"LCMMethod is an option for HelixAngleCalc and LMCSytemCalc. Can be \"CentralAxes\" or \"WallMap\". 
\"CentralAxes\" uses wall distance calculation using projection of the centarl axes and circular approximation of the ventricle. This method is fairly fast and uses CentralAxes internaly.
\"WallMap\" uses wall distance interpolation and subsequential gradient calculation. Can take long for high res datasets but is most accurate. Uses CalculateWallMap internaly."

AxesMethod::usage = 
"AxesMethod is an option for HelixAngleCalc and CentralAxes and CardiacCoordinateSystem. Can be \"Linear\", \"Quadratic\", \"Cubic\"."

RowSize::usage =
"RowSize is an option for CentralAxes. defines the number or images per showing the segmentation.
Can be \"Automatic\" of an integer." 

ShowPlot::usage = 
"ShowPlot is an option for CentralAxes, HelixAngleCalc and CardiacCoordinateSystem. True shows the fit of the central axes and outpu the plot as extra output."

MaskWallMap::usage = 
"MaskWallMap is an option for CalculateWallMap. if True or False."


GroupPerSegment::usage =
"GroupPerSegment is an option for SegmentsPerSlice. If set False segements are grouped per slice and not per segment."

SegmentationMethod::usage =
"SegmentationMethod is an option for SegmentsPerSlice. Values can be \"AHA\", \"AHA+\", 1, 2, 3, 6 or 8."


ReversePoints::usage = 
"ReversePoints is an option for LinesToSegmentIndex, CardiacSegment. Defines at which point to start, can be True or False."

ReverseDirection::usage = 
"ReverseDirection is an option for LinesToSegmentIndex, CardiacSegment. Defines the direction of rotiation, clockwise or anti-clockwise, can be True of False."

MakeSegmentPlots::usage = 
"MakeSegmentPlots is an option for CardiacSegment. If True plots of the sementation are made."

StartPoints::usage = 
"StartPoints is an option for CardiacSegmentGUI. Value is \"Default\" or the point list given by CardiacSegment."

StartSlices::usage = 
"StartSlices is an option for CardiacSegmentGUI. Value is \"Default\" or the list given by CardiacSegment."


RadialSamples::usage =
"RadialSamples is an option for RadialSample and PlotSegments. Defines how manny transmural samples are taken."

DropSamples::usage = 
"DropSamples is an option for RadialSample and PlotSegments. Defines how manny samples are droped form star and end. Can be an number or set (strat, end) of numbers."


GridLineSpacing::usage = 
"GridLineSpacing is an option of TransmuralPlot. It defines the spacing of the gridlines."

ShowOutliers::usage = "ShowOutliers is an option for ExcludeSlices."

SmoothHelix::usage = 
"SmoothHelix is an option for MaskHelix, sets the kernelsize for the MedianFilter." 

BackgroundValue::usage = 
"BackgroundValue is an option for MaskHelix. Sets the backgroud value (default is -100)."

TextOffset::usage = 
"TextOffset is an option for BullseyePlot. Determines where the text is placed, can be 0 to 1."

TextSize::usage = 
"TextSize is an option for BullseyePlot. Determines the text size."

TextNumberForm::usage =
"TextNumberForm is an option for BullseyePlot. Specifies how many number and decimals to use like in NumberForm."

BullPlotMethod::usage = 
"BullPlotMethod is an option for BullseyePlot. Can be \"Dynamic\" of \"Normal\". 
\"Dynamic\" allows to change plotting parameters in Manipulation window."

CutOffMethod::usage =
"CutOffMethod is an option for ExcludeSlices. Default value is \"Auto\" or it can be a fixed percentage (value between 0 and .5)."

DistanceMeasure::usage = 
"DistanceMeasure is an option for ExcludeSlices. Defaul value is 5. (1 ManhattanDistance, 2 SquaredEuclideanDistance, 3 EuclideanDistance, 4 Correlation, 5 SpearmanRho."

BloodMaskRange::usage = 
"BloodMaskRange is an option for MakeECVBloodMask." 

OutputCheckImage::usage =
"OutputCheckImage is an option for MakeECVBloodMask." 	


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*General functions*)



(* ::Subsubsection::Closed:: *)
(*DotC*)


DotC = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}}, 
	vec1 . vec2, 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MakePerpendicular*)


MakePerpendicular = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}}, 
	Normalize[vec1 - (vec1 . vec2) vec2], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*NormalizeC*)


NormalizeC = Compile[{{vec, _Real, 1}}, 
	Normalize[vec], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*CrossC*)


CrossC = Compile[{{vec1, _Real, 1}, {vec2, _Real, 1}}, 
	Cross[vec1, vec2], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*PlotMaskVolume*)


Options[PlotMaskVolume] = {Filter->True}

PlotMaskVolume[mask_,vox_,color_:Darker[Gray],OptionsPattern[]] := Block[{pmask,dim ,xi,yi,zi,pmask2},
	pmask = ArrayPad[Reverse[If[OptionValue[Filter], GaussianFilter[Clip[mask], 1],	Clip[mask]], 2], 1];
	dim=Dimensions[pmask];
		
	ListContourPlot3D[
		(*fix for 13.1*)
		RescaleData[pmask,Reverse@dim], 
		MaxPlotPoints->Infinity, Contours -> {.4}, ContourStyle -> {color, Opacity[.5]}, Mesh -> False,
		Lighting -> "Neutral", BoundaryStyle -> None, PlotRange -> (Thread[{{0, 0, 0}, Reverse@dim - 1}]), DataRange->(Thread[{{0, 0, 0}, dim - 1}]),
		BoxRatios -> Reverse[(vox (dim + 2))], Axes -> True, ImageSize -> 400, SphericalRegion -> True
	]
]


(* ::Subsection::Closed:: *)
(*HelixAngleCalc*)


Options[HelixAngleCalc]={ShowPlot->True, LCMMethod->"WallMap", AxesMethod->"Quadratic"};

SyntaxInformation[HelixAngleCalc]={"ArgumentsPattern"->{_,_,_,_.,OptionsPattern[]}};

HelixAngleCalc[data_?ArrayQ, mask_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, opts:OptionsPattern[]]:=HelixAngleCalc[data,mask,0,vox,opts]

HelixAngleCalc[data_?ArrayQ, mask_?ArrayQ, maskp_, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, opts:OptionsPattern[]]:=Block[
	{out, evec, projection, inp, helix, sign, norvec,radvec, cirvec, coors, plots, i, j},

	coors = CardiacCoordinateSystem[mask, maskp, vox, 
		ShowPlot->OptionValue[ShowPlot], 
		LCMMethod->OptionValue[LCMMethod], 
		AxesMethod->OptionValue[AxesMethod]
	];
	
	{radvec, norvec, cirvec} = If[OptionValue[ShowPlot], plots = coors[[2]]; coors[[1]], coors];

	(*create helix angle maps*)
	out = Flatten[Table[
		evec = Map[Reverse, data[[All,All,All,All,i]], {3}];
		(*align vector with projection vector*)
		evec = SignNoZero[DotC[{norvec,radvec,cirvec}[[j]],evec]] evec;
		(*Helix,Transverse,Sheet}*)
		projection = MakePerpendicular[evec, {radvec,cirvec,norvec}[[j]]];
		inp = DotC[projection,{cirvec,norvec,radvec}[[j]]];
		helix = ArcCos[Abs[inp]]/Degree;
		
		(*keep sign only for helix[[1,1]]*)
		If[(j==1 && i==1), SignNoZero[inp], 1] helix;
	, {i, 3}, {j, 3}],1];

	If[OptionValue[ShowPlot], {Re@out,plots}, Re@out]
]


(* ::Subsection:: *)
(*CardiacCoordinateSystem*)


(* ::Subsubsection::Closed:: *)
(*CardiacCoordinateSystem*)


Options[CardiacCoordinateSystem] = {ShowPlot -> False, LCMMethod->"WallMap", AxesMethod->"Quadratic"}

SyntaxInformation[CardiacCoordinateSystem] = {"ArgumentsPattern" -> {_,_, OptionsPattern[]}};

CardiacCoordinateSystem[mask_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, opts:OptionsPattern[]]:= CardiacCoordinateSystem[mask, 0, vox, opts]

CardiacCoordinateSystem[mask_?ArrayQ, maskp_, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, OptionsPattern[]] := Block[{
		dim, axesout, off, vec, inout, pla, met, radvec, norvec, cirvec,norvecc, radvecn, wall, der,
		sp, spz, spxy, maskCont, vectorField,n, z, y, x, coo, rav, nov, rov, vec1 ,vec2, vec3, plot, plw 
	},
	
	dim = Dimensions[mask];
	
	(*get the cardiac center line*)
	axesout = CentralAxes[mask, maskp, vox, AxesMethod -> OptionValue[AxesMethod], ShowPlot -> OptionValue[ShowPlot]];
	{off, vec, inout} = If[OptionValue[ShowPlot], pla = axesout[[4]]; axesout[[1 ;; 3]], axesout];
	
	PrintTemporary["LMCS caclulation start"];
	Switch[OptionValue[LCMMethod],
		"CentralAxes",
		plw = Nothing;
		(*calculate the wall angle map*)
		wall = N[WallAngleMap[mask, vox, inout] Degree];
		(*define te rad vector using center points*)
		radvec = RadVecC[dim, off];
		(*define norvecs using centerline vectors*)
		norvec = ConstantArray[vec[[#]],dim[[2;;]]]&/@Range[dim[[1]]];
		(*define rotation vectors perpendicualr to radvecn and norvec*)
		cirvec = NormalizeC[CrossC[radvec,norvec]];
		(*correct the norvec for the wall curvature by rotation around rotvec*)
		norvecc = NorVecR[wall,cirvec,norvec];
		(*make radvec purpendicular to corrected norvec*)
		radvecn = MakePerpendicular[radvec,norvecc];
		,
		"WallMap",
		(*Calculate the wall distance map, and the wall direction*)
		wall = CalculateWallMap[mask, vox, ShowPlot -> OptionValue[ShowPlot]];
		{wall, der} = If[OptionValue[ShowPlot], plw = wall[[3]]; wall[[1 ;; 2]], wall];
		(*make all vectors perpendicual*)
		radvecn = NormalizeC[Transpose[der/vox, {4, 1, 2, 3}]];
		norvec = ConstantArray[#, dim[[2 ;;]]] & /@ vec;
		norvecc = MakePerpendicular[norvec, radvecn];
		cirvec = NormalizeC[CrossC[radvecn,norvecc]];
	];
	
	PrintTemporary["LMCS is done"];
	If[OptionValue[ShowPlot],
		(*plot hear geo*)
		sp = Ceiling[Dimensions[mask]/{12, 24, 24}];
		{spz, spxy} = {sp[[1]], Min[sp[[2 ;; 3]]]};
		maskCont = PlotContour[Reverse[Clip[mask + maskp,{0,1},{0,1}],2], vox, ContourColor->Gray];
		n = (spz 0.6 vox[[1]]) {1, -1, 1}/vox;
		vectorField = Table[
			If[mask[[z, y, x]] == 0,
				{None, None, None},
				coo = {x, -y + dim[[2]] + 1, z};
				rav = Reverse[n radvecn[[z, y, x]]];
				nov = Reverse[n norvecc[[z, y, x]]];
				rov = Reverse[n cirvec[[z, y, x]]];
				{
					{Darker[Green], Thick, Line[{coo(*-rav*), coo + rav}]},
					{Darker[Blue], Thick, Line[{coo(*-nov*), coo + nov}]}, 
					{Darker[Red], Thick, Line[{coo(*-rov*), coo + rov}]}
				}
			], 
		{z, 1, dim[[1]], spz}, {y, 1, dim[[2]], spxy}, {x, 1, dim[[3]], spxy}];
		
		{vec1, vec2, vec3} =   DeleteCases[Flatten[#, 2], None] & /@ Transpose[vectorField, {2, 3, 4, 1}];
	    plot = Show[maskCont, Graphics3D[vec1], Graphics3D[vec2], Graphics3D[vec3]];
	    
	    PrintTemporary[plot];
	    Pause[1];
	];
	
	If[OptionValue[ShowPlot],
		{{radvecn, norvecc, cirvec}, {pla, plw, plot}},
		{radvecn, norvecc, cirvec}
	]
]


(* ::Subsubsection::Closed:: *)
(*RadVecC*)


RadVecC = Compile[{{dim, _Real, 1}, {off, _Real, 2}}, 
	Table[Normalize[{i, j, k}-off[[i]]], {i, dim[[1]]}, {j, dim[[2]]}, {k, dim[[3]]}], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*NorVecR*)


NorVecR = Compile[{{angmap, _Real, 0}, {rotvec, _Real, 1}, {norvec, _Real, 1}}, Block[{v1, v2, v3, W, iden},
	{v1, v2, v3} = rotvec;
	iden = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	W = -{{0, -v3, v2}, {v3, 0, -v1}, {-v2, v1, 0}};
	(iden + Sin[angmap] W + (2 Sin[angmap/2]^2 MatrixPower[W, 2])) . norvec], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*WallAngleMap*)


WallAngleMap[mask_, vox_, inout_] := Block[{dim, cent, len, edge1, edge2, fout, fin, walldir, wallang, wallvec, sign, wallangfunc, dist, in, out},
	dim = Dimensions[mask];
	len = dim[[1]];
	
	(*fit the wall profile with Quadratic function*)
	{in, out} = (FitWall[#, "Quartic"]) & /@ {inout[[1, 2]], inout[[2, 2]]};
	fout = First@First@Position[out[[1]], _Real?(# != 0. &)];
	fin = First@First@Position[in[[1]], _Real?(# != 0. &)];
	
	walldir = {in, out};
	wallang = Transpose[Map[(
		wallvec = Normalize[vox {1, #[[2]], 0}];
		sign = SignNoZero[#[[2]]];
		{#[[1]], sign VectorAngle[{1, 0, 0}, wallvec]/Degree}
	) &, Transpose[walldir, {2, 3, 1}], {2}], {3, 1, 2}];
	
	cent = {ConstantArray[0, len], wallang[[1, 2]]};

	If[fout > 1, wallang[[2, 2, Range[fout - 1]]] = 90];
	If[fin > 1, wallang[[1, 2, Range[fin - 1]]] = 90];
	cent = {ConstantArray[0, len], wallang[[1, 2]]};

	edge1 = {ConstantArray[ Max[walldir[[All, 1]]], len], wallang[[2, 2]]};
	edge2 = {ConstantArray[ Max[vox dim], len], wallang[[2, 2]]};
	wallang = Transpose[{Transpose[vox[[1 ;; 2]] {Range[len], #[[1]]}], #[[2]]}] & /@ Join[{cent}, wallang, {edge1, edge2}];
	
	wallang = Sort[DeleteDuplicates[Flatten[wallang, 1]]];
	wallangfunc = Interpolation[wallang, InterpolationOrder -> 1];
	
	Quiet@Table[
		dist = EuclideanDistance[{j, k}, inout[[1, 1, z, 2 ;;]]];
		wallangfunc[vox[[1]] z, vox[[2]] dist]
	, {z, dim[[1]]}, {j, dim[[2]]}, {k, dim[[3]]}
   ]
]


(* ::Subsubsection::Closed:: *)
(*FitWall*)


(*fit the waal profile for the normalized central axes*)
FitWall[data_, met_] := Block[{fun, pf, fdata, xdat, points, pos},
	(*define function and data*)
	fun = Switch[met,
		"Quadratic", {1, t, t^2},
		"Cubic", {1, t, t^2, t^3},
		"Quartic", {1, t, t^2, t^3, t^4},
		_, {1, t}];
	xdat = Range[Length[data]];
	fdata = DeleteCases[Transpose[{xdat, data}], {_, {}}];
	
	(*perform the fit and gererate the fitted poins *)
	pf = Simplify@Chop@Fit[fdata, fun, t];
	points = Chop[pf /. t -> # & /@ xdat];
	
	(*set the radius to zero for slices without a mask*)
	pos = Flatten[Position[data, {}]];
	If[! (pos === {}), points[[pos]] = 0];
	{points, (Unitize[points[[#]]] (D[pf, t] /. t -> #) /.  0. -> 10) & /@ xdat}
]


(* ::Subsection:: *)
(*CentralAxes*)


(* ::Subsubsection::Closed:: *)
(*CentralAxes*)


Options[CentralAxes]={ShowPlot->False, RowSize->"Automatic", AxesMethod->"Cubic"(*,Output->True*)};

SyntaxInformation[CentralAxes]={"ArgumentsPattern"->{_,_,_.,OptionsPattern[]}};

CentralAxes[mask_,vox_,opts:OptionsPattern[]]:=CentralAxes[mask,0,vox,opts]

CentralAxes[mask_,maskp_,vox_,OptionsPattern[]]:=Module[{
	rad,met,row,dim,half,minmaxr,inner,outer,inout,offi,offo,vecsi,vecso,
	offouti,offouto,pl1,off,vecs,pl2,offout,vecsout,fit,last},
	
	(*get option values*)
	rad={0.01,1};
	met=OptionValue[AxesMethod];
	row=If[OptionValue[RowSize]==="Automatic"||!IntegerQ[OptionValue[RowSize]],
		Round[Sqrt[Length[mask]]],
		OptionValue[RowSize]
	];
	
	(*get data dimensions*)
	dim=Dimensions[mask];
	(*half=CenterPoint[mask];*)
	half=Drop[dim,1]/2.;
	minmaxr= rad Max[(Drop[dim,1]/1)];
	
	(*get inner and outer radius*)
	{inner,outer} = GetRadius[mask];

	(*finde the upper most closed outer*)
	last = First@Last@Position[Unitize[inner[[3]] /. {} -> 0] + Unitize[outer[[3]] /. {} -> 0], 2];
	outer[[All, last + 1 ;;]] = Transpose@ConstantArray[{{}, {}, {}}, Length[outer[[3]]] - last];
	
	(*fit off centers*)
	{off, vecs} = FitCenterLine[inner[[1]], outer[[1]], vox, met];
	{offi, vecsi} = FitCenterLine[inner[[1]], vox, met];
	{offo, vecso} = FitCenterLine[outer[[1]], vox, met];
	
	{off, vecs} = BoundCorrect[Min /@ Transpose[{inner[[3]], outer[[3]]} /. {{} -> 0}], off, vecs];
	{offi, vecsi} = BoundCorrect[inner[[3]]/. {{} -> 0}, offi, vecsi];
	{offo, vecso} = BoundCorrect[outer[[3]]/. {{} -> 0}, offo, vecso];
	  
	(*generate plots*)
	If[OptionValue[ShowPlot],
		pl1 = PlotRadius[Clip[2 mask + maskp,{0,2}], inner, outer];
		pl2 = PlotSegmentation[mask + maskp, inner, outer, {off, offi, offo}, vox];
		fit = Row[{GraphicsGrid[Partition[pl1, row, row, 1, {}], ImageSize -> row*100], pl2}];
	];
	
	(*create output*)
	offout = Reverse[{1, -1, 1} (# + {0, -(dim[[2]]), 0})] & /@ off;
	offouti = Reverse[{1, -1, 1} (# + {0, -(dim[[2]]), 0})] & /@ offi;
	offouto = Reverse[{1, -1, 1} (# + {0, -(dim[[2]]), 0})] & /@ offo;
	vecsout=Reverse[{1,-1,1} #]&/@vecs;
	inout = {{offouti, inner[[3]]}, {offouto, outer[[3]]}};
	
	(*If[OptionValue[ShowPlot]==True,Print[fit]];*)
	If[OptionValue[ShowPlot],{offout,vecsout,inout,fit},{offout,vecsout,inout}]
]


(* ::Subsubsection::Closed:: *)
(*GetRadius*)


GetRadius[mask_] := Block[{comps, in, out, fout, fin, seg, min ,mout},
	(*create the inner and outer volume*)
	seg = MorphologicalComponents[N[#]] & /@ (1. - mask);
	seg = If[#[[1, 1]] == 2, # /. {2 -> 1, 1 -> 2}, #] & /@ seg;
	min = Unitize[Clip[seg - 1, {0, 1}, {0, 0}]];
	mout = Unitize[1 - Clip[seg, {0, 1}, {0, 0}]];
	comps =Reverse /@ ( (mout - (Erosion[#, 1] & /@ mout)) + 2 ((Dilation[#, 1] & /@ min) - min));
			
	{in, out} = Transpose[Switch[Max[#],
		1, {{}, Reverse /@ Position[#, 1]},
		2, SortBy[{Reverse /@ Position[#, 1], Reverse /@ Position[#, 2]}, (Length[#] &)],
		_, {{}, {}}
	] & /@ comps];
	
	fout = Transpose[MapIndexed[If[#1 === {}, {{}, {}, {}}, FitEllipse[#1, #2]] &, out]];
	fin = Transpose[MapIndexed[If[#1 === {}, {{}, {}, {}}, FitEllipse[#1, #2]] &, in]];
	
	{fin, fout}
]


(* ::Subsubsection::Closed:: *)
(*FitEllipse*)


FitEllipse[pts_, i_] := Block[{lsMat, a, b, c, d, e, f, s, val ,vec,  den, num ,fac, r0, r1,r2,ang},
	(*transform coordinates*)
	lsMat = N[Function[{x, y}, {x^2, 2 x y, y^2, 2 x, 2 y, 1}] @@@ pts];
	(*solve eclipse equation*)
	(*{a, b, c, d, e} = Chop[LeastSquares[lsMat, ConstantArray[1., Length[pts]]]];*)
	d = Transpose[lsMat] . lsMat;
	c = {{0, 0, -2, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {-2, 0, 0, 0, 0, 0}, 
		{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
	{val, vec} = Eigensystem[{d, c}];
	{a, b, c, d, e, f} = vec[[Position[Boole[Negative /@ val], 1][[1, 1]]]];
		
	den = b^2 - a c;
	num = 2*(a e^2 + c d^2 + f b^2 - 2 b d e - a c f);
	fac = Sqrt[((a - c)^2 + 4 b^2)];
	
	If[(b^2 - a c) < 0, 
		(*solve for centers*)
		r0 = {(c*d - b*e)/den, (a*e - b*d)/den};
		(*get the radisu*)
		r1 = Sqrt[num/den/(fac - a - c)];
		r2 = Sqrt[num/den/(-fac - a - c)];
		(*get the angle*)
		ang = ArcTan[(2 b)/(a - c)]/2;
		If[ a > c, ang += Pi/2];
		{Join[r0, i], {r1, r2, ang}, Mean[{r1, r2}]},
		{{}, {}, {}}
	]
]


(* ::Subsubsection::Closed:: *)
(*FitCenterLine*)


(*fit the center line for the average inner and outer center points*)
FitCenterLine[inner_, outer_, vox_, met_] := Block[{data},
	data = Mean[DeleteCases[#, {}]] & /@ Transpose[{inner, outer}];
	FitCenterLine[data, vox, met]
]

FitCenterLine[datai_, vox_, met_] := Block[{fun, pf, xdat, data, fdata},
	data = datai /. {Mean[{}] -> {{}, {}, {}}, {} -> {{}, {}, {}}};
	
	fun = Switch[met,
		"Quadratic", {1, t, t^2},
		"Cubic", {1, t, t^2, t^3},
		"Quartic", {1, t, t^2, t^3, t^4},
	_, {1, t}];
	
	xdat = Range[Length[data]];
	fdata = DeleteCases[Transpose[{Range[Length[#]], #}], {_, {}}] & /@ Transpose[data];
	pf = Simplify@Chop@Fit[#1, fun, t] & /@ fdata;
	
	{Transpose[If[NumberQ[#], ConstantArray[#, Length[xdat]], # /. t -> xdat] & /@ pf],
	(Normalize[D[Reverse[vox] pf, t] /. t -> #] & /@ xdat)}
]


(* ::Subsubsection::Closed:: *)
(*BoundCorrect*)


BoundCorrect[points_,off_,vec_]:=Block[{vv,first,last,offo,veco},
	vv = Unitize[points];
		
	first = (First@Position[vv, 1])[[1]];
	last = (Last@Position[vv, 1])[[1]];
	offo = Join[
		{off[[first, 1]] + (first - #) (off[[first, 1]] - off[[first + 1, 1]]), off[[first, 2]] + (first - #) (off[[first, 2]] - off[[first + 1, 2]]), #} & /@ Range[1, first - 1],
		off[[Range[first, last]]],
	  	{off[[last, 1]] + (# - last) (off[[last, 1]] - off[[last - 1, 1]]), off[[last, 2]] + (# - last) (off[[last, 2]] - off[[last - 1, 2]]), #} & /@ Range[last + 1, Length[vv]]
	  ];
	veco=Join[ConstantArray[vec[[first]],first-1],vec[[first;;last]],ConstantArray[vec[[last]],Length[vv]-last]];
	{offo,veco}
]


(* ::Subsubsection::Closed:: *)
(*PlotRadius*)


PlotRadius[mask_, inner_, outer_] := MapThread[(
	Show[
		Image[#1/Max[mask], ImageSize -> 100],
		If[#2[[1]] != {}, ParametricPlot[Ellipse2D[#2[[{1, 2}]], t] - 0.5, {t, 0, 2 Pi}, PlotStyle -> Directive[Thick, Red]], Graphics[]],
		If[#2[[1]] != {}, Graphics[{Red, PointSize[Medium], Point[#2[[1, 1 ;; 2]] - 0.5]}], Graphics[]],
		If[#3[[1]] != {}, ParametricPlot[Ellipse2D[#3[[{1, 2}]], t] - 0.5, {t, 0, 2 Pi}, PlotStyle -> Directive[Thick, Blue]], Graphics[]],
		If[#3[[1]] != {}, Graphics[{Blue, PointSize[Medium], Point[#3[[1, 1 ;; 2]]]}], Graphics[]]
	]
) &, {mask, Transpose@inner, Transpose@outer}]


Ellipse2D[{{x0_, y0_, z0_}, {a_, b_, alpha_}}, theta_] := {
	x0 + a*Cos[theta]*Cos[alpha] - b*Sin[theta]*Sin[alpha],
	y0 + a*Cos[theta]*Sin[alpha] + b*Sin[theta]*Cos[alpha]}


Ellipse3D[{{x0_, y0_, z0_}, {a_, b_, alpha_}}, theta_] := {
	x0 + a*Cos[theta]*Cos[alpha] - b*Sin[theta]*Sin[alpha],
	y0 + a*Cos[theta]*Sin[alpha] + b*Sin[theta]*Cos[alpha], 
	z0}


(* ::Subsubsection::Closed:: *)
(*PlotSegmentation*)


PlotSegmentation[mask_, inner_, outer_, {off_, offi_, offo_}, vox_] := Block[{voxl, offp, offip, offop},
	voxl = Reverse[vox];
	
	offp = {.5, .5, 0} + {1, 1, 1} # & /@ DeleteCases[off, {}];
	offip = {.5, .5, 0} + {1, 1, 1} # & /@ DeleteCases[offi, {}];
	offop = {.5, .5, 0} + {1, 1, 1} # & /@ DeleteCases[offo, {}];
	
	Show[
		PlotContour[Reverse[mask,2], vox],
		(*center points*)
		ListPointPlot3D[offp, PlotStyle -> Directive[{Thick, Black, PointSize[Large]}]],
		Graphics3D[{Thick, Black, Line[offp]}],
		
		ListPointPlot3D[Delete[outer[[1]], Position[outer[[2]], {}]], PlotStyle -> Directive[{Thick, Blue, PointSize[Large]}]],
		Graphics3D[{Thick, Blue, Line[Delete[offop, Position[outer[[2]], {}]]]}],
		
		ListPointPlot3D[Delete[offip, Position[inner[[2]], {}]], PlotStyle -> Directive[{Thick, Red, PointSize[Large]}]],
		Graphics3D[{Thick, Red, Line[Delete[offip, Position[inner[[2]], {}]]]}],
		
		(*Plot the segmented outlines*)
		If[outer[[2, #]] === {}, 
			Graphics3D[],
			ParametricPlot3D[Ellipse3D[outer[[{1, 2}, #]], u] - {0.5, .5, 0}, {u, 0, 2 Pi}, PlotStyle -> Directive[{Thick, Blue, Opacity[.5]}]] 
			(*ParametricPlot3D[{outer[[3, #]] Sin[u], outer[[3, #]] voxl[[1]]/voxl[[2]] Cos[u],0} + offop[[#]], {u, 0, 2 Pi}, PlotStyle -> Directive[{Thick, Blue, Opacity[.5]}]]*)
		] & /@ Range[Length[mask]],
		
		If[inner[[2, #]] === {}, 
			Graphics3D[], 
			ParametricPlot3D[Ellipse3D[inner[[{1, 2}, #]], u] - {0.5, .5, 0}, {u, 0, 2 Pi}, PlotStyle -> Directive[{Thick, Red, Opacity[.5]}]]
			(*ParametricPlot3D[{inner[[3, #]] Sin[u], inner[[3, #]] voxl[[1]]/voxl[[2]] Cos[u], 0} + offip[[#]], {u, 0, 2 Pi}, PlotStyle -> Directive[{Thick, Red, Opacity[.5]}]]*)
		] & /@ Range[Length[mask]],
		
	ImageSize->400]
]


(* ::Subsection:: *)
(*CalculateWallMap*)


(* ::Subsubsection::Closed:: *)
(*CalculateWallMap*)


Options[CalculateWallMap] = {ShowPlot -> True, MaskWallMap->False};

SyntaxInformation[CalculateWallMap] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

CalculateWallMap[maski_, vox_, OptionsPattern[]] := Module[{
		seg, min, mout, mtot, pts, ptsi, pos, ptso, x, y, z, plane, mask, dim,
		surfpl, pointspl, planepl, planefit, fit, planem, d1, d2, d3, zc,
		min2, mout2, surfin, surfout, ptsin, ptsout, ptspl, inpnt, outpnt,
		wall, pt, pt1, pt2, dist1, ptsm, ptsmv, dist, der, mn ,mx, ind, max
	},
	
	(*prepare the mask*)
	mask = SelectComponents[#, "Count", -1] & /@ Round[maski];
	dim = Dimensions[mask];
	
	(*create the inner and outer volume*)
	seg = MorphologicalComponents[N[#]] & /@ (1. - mask);
	seg = If[#[[1, 1]] == 2, # /. {2 -> 1, 1 -> 2}, #] & /@ seg;
	min = Unitize[Clip[seg - 1, {0, 1}, {0, 0}]];
	mout = Unitize[1 - Clip[seg, {0, 1}, {0, 0}]];
	mtot = (mout + 1) - 2 min;
	
	(*get the top points to fit top plane*)
	ind = Range[Length[mask]] - 1;
	pts = ptsi = Flatten[MapIndexed[(If[(pos = Pick[ind, #1, 1]) != {}, Reverse@Flatten[{Last@pos, #2}], Nothing]) &, RotateDimensionsLeft[mask], {2}], 1];
	max = Max[pts[[All, 3]]] - Ceiling[Length[mask]/10];
	ptso = Select[pts, #[[3]] > max &];
	
	(*fit the top plane*)
	Clear[x, y, z];
	plane = Fit[ptso, {1, x, y}, {x, y}];
	While[ptsi != ptso,
		ptsi = ptso;
		ptso = Select[pts, (((plane - z) /. Thread[{x, y, z} -> #]) < 1) &];
		plane = Fit[ptso, {1, x, y}, {x, y}];
	];
	
	(*plane fit visualisation*)
	If[OptionValue[ShowPlot],
		(*surfpl = ListContourPlot3D[
			(*fix for 13.1*)
			RescaleData[GaussianFilter[mask, 1], Reverse@dim], 
			
			Contours -> {0.6}, Mesh -> False, PlotRange -> Transpose[{{0, 0, 0}, Reverse@(dim)}],
			BoxRatios -> Reverse@(vox dim), ContourStyle -> Directive[Gray, Opacity[0.5]], Lighting -> "Neutral", Axes -> False];*)
		surfpl = PlotContour[mask, vox, ContourColor->Gray];
		pointspl = ListPointPlot3D[
			ptsi, PlotRange -> Transpose[{{0, 0, 0}, Reverse@dim}], BoxRatios -> 1, PlotStyle -> Red];
		planepl = Plot3D[plane,  {x, 0, dim[[3]]}, {y, 0, dim[[2]]}, Mesh -> False, PlotStyle -> Directive[Red, Opacity[.2]], BoundaryStyle -> Darker[Red]];
		planefit = Show[surfpl, planepl, pointspl, PerformanceGoal -> "Speed", ImageSize->350];
	];
	
	(*make mask from plan*)
	planem = 0 mout;
	{d1, d2, d3} = Dimensions[planem];
	Table[
		zc = (Round[plane] + o) /. {x -> xc, y -> yc};
		If[zc <= d1, planem[[zc, yc, xc]] = 1]
		, {xc, 1, d3}, {yc, 1, d2}, {o, 1, 2}];
	planem = (1 - mout) planem;
	
	(*close inner mask to plane*)
	min2 = UnitStep[RotateDimensionsRight[FillFun[RotateDimensionsLeft[planem + mout + min]]] - 2];
	mout2 = Clip[mout + min2, {0, 1}];
	
	ptsin = Position[min2 - (Erosion[#, 1] & /@ min2), 1];
	ptsout = Position[(Dilation[#, 1] & /@ mout2) - mout2, 1];
	ptsm = Position[mask, 1];
	
	If[OptionValue[ShowPlot],
		(*Create Inner and outer surfaces*)
		{surfout, surfin} = PlotContour[{mout2, min2}[[#]], vox, ContourColor->{Red, Blue}[[#]]] & /@ {1, 2};
			(*{surfout, surfin} = ListContourPlot3D[
			(*fix for 13.1*)
			RescaleData[GaussianFilter[min2 + mout2, 1], Reverse@dim], 
			
			Contours -> {0.3, 1.5}[[#]], Mesh -> False, 
			PlotRange -> Transpose[{{0, 0, 0}, Reverse@(dim)}], Lighting -> "Neutral", 
			BoxRatios -> Reverse@(vox dim), ContourStyle -> Directive[{Red, Blue}[[#]], Opacity[0.4]], 
			MaxPlotPoints -> {50, 100}] & /@ {1, 2};*)
		
		ptspl = Show[surfin, surfout, ListPointPlot3D[{Reverse[#]-1&/@ptsin[[;; ;;2]], Reverse[#]-1&/@ptsout[[;; ;;2]]}, PlotStyle -> {Blue, Red}], PerformanceGoal -> "Speed",ImageSize->350]
	];
	
	inpnt = Nearest[N[vox # & /@ ptsin]];
	outpnt = Nearest[N[vox # & /@ ptsout]];
	ptsmv = N[vox # & /@ ptsm];
	
	(*generate the wall distance function*)
	min = 0.1 Norm[Drop[vox,1]];
	dist = Map[(
		pt = #;
		pt1 = Mean@inpnt[pt, 1];
		pt2 = Mean@outpnt[pt, 1];
		dist1 = Norm[pt1 - pt2];
		If[dist1 < min, 0., Mean[{Norm[pt1 - pt]/dist1, 1 - (Norm[pt2 - pt]/dist1)}]]
	) &, ptsmv];
	
	(*create the wall distance map*)
	wall = GaussianFilter[Normal[SparseArray[Thread[ptsm -> dist], dim]] + (1 - mout2), 1];
	der = GaussianFilter[wall, {{3,3,3}}(*{2 (1/vox)/(1/vox[[2]])}*), #] & /@ (IdentityMatrix[3]);
    
    If[OptionValue[MaskWallMap],
    	wall = mask wall;
    	der = mask # & /@ der;
    	{mn, mx} = Quantile[Flatten[GetMaskData[wall, mask]],{0.01,0.99}];
    	wall = Clip[mask ((wall - mn)/(mx-mn)),{0,1}];
    ];
    
    If[OptionValue[ShowPlot], PrintTemporary[Row[fit = {planefit, ptspl}]]; {wall, der, fit}, {wall, der}]
]


(* ::Subsubsection::Closed:: *)
(*FillFun*)


FillFun = Compile[{{v, _Integer, 1}}, Block[{out, in, i, m},
	out = v;
	in = 0 v;
	m = 1 - UnitStep[v - 2];
	i = 0;
	If[Last[out] =!= 2,
		While[in =!= out && i < 10,
			i++;
			in = out;
			out = out + m 2 Prepend[UnitStep[Abs[Differences[out]] - 2], 0]
		];
	out, out]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*GetMaskSegmentPoints*)


(* ::Subsubsection::Closed:: *)
(*GetMaskSegmentPoints*)


SyntaxInformation[GetMaskSegmentPoints] = {"ArgumentsPattern" -> {_}};

GetMaskSegmentPoints[maski_]:=Block[{mask, lab, m1, m2, points},
	(*split in main mask and points*)
	{mask, lab} = SplitSegmentations[maski];
	If[Round[lab]=!={1,2,3},
		
		(*not valid mask*)
		Print["not valid input mask"];
		{maski,$Failed},
		
		(*valid mask*)
		{mask, m1, m2} = Transpose[Normal[mask]];
		
		(*get point coordinates*)
		points = Transpose[{PointCenter/@m1, PointCenter/@m2}];
		
		(*remove points from mask*)
		mask = RemovePoint[m2, RemovePoint[m1, mask]];
		
		(*output*)
		{mask, points}
	]
]


(* ::Subsubsection::Closed:: *)
(*PointCenter*)


PointCenter = (1 /. ComponentMeasurements[SelectComponents[Image[Reverse[#]], "Count", -1], "Centroid"]) /. {1 -> None, {x_, y_} -> {y, x}}&


(* ::Subsubsection::Closed:: *)
(*RemovePoint*)


RemovePoint[m_, msk_] := Block[{v}, 
	MapThread[(
		v = GetMaskData[#2, Dilation[#1, 1] - #1, GetMaskOnly -> True];
		v = If[v === {}, 0, Round[Mean[v]]];
		If[v === 1, #2 + #1, #2]
	) &, {m, msk}]
]


(* ::Subsection::Closed:: *)
(*GetSegmentSlices*)


SyntaxInformation[GetSegmentSlices] = {"ArgumentsPattern" -> {_}};

GetSegmentSlices[mask_?ArrayQ] := GetSegmentSlicesi[Position[Max[#] & /@ Clip[Round[mask], {0, 1}], 1]]

GetSegmentSlices[pts_] := GetSegmentSlicesi[Position[If[#[[1]] =!= None && #[[2]] =!= None, 1, 0] & /@ pts, 1]]

GetSegmentSlicesi[pos_] := Block[{st, en, app, seg},
	st = First@First[pos - 1];
	en = First@Last[pos];
	app = Round[(en - st)/7];
	seg = Ceiling@Range[st + app, en - (en - (app + st))/3, (en - (app + st))/3];
	Flatten[{st, seg, en}]
]


(* ::Subsection::Closed:: *)
(*SegmentsPerSlice*)


Options[SegmentsPerSlice]={GroupPerSegment ->True, SegmentationMethod -> "AHA"}

SyntaxInformation[SegmentsPerSlice] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SegmentsPerSlice[slices_?VectorQ, ops:OptionsPattern[]]:= SegmentsPerSlice[slices, None, ops]

SegmentsPerSlice[points_, ops:OptionsPattern[]]:= SegmentsPerSlice[GetSegmentSlices[points], points, ops]

SegmentsPerSlice[slices_?VectorQ, points_, OptionsPattern[]]:=Block[
	{segmi, grp, loc, segs},
	
	segmi = OptionValue[SegmentationMethod];
	grp = OptionValue[GroupPerSegment]||StringQ[segmi];
	
	(*segment ranges*)
	loc = Range@@(#+{1,0})&/@Partition[slices, 2, 1];
	(*select only sices that have 2 point markers*)
	If[points=!=None, loc = Select[#, (points[[#,1]]=!=None && points[[#,2]]=!=None)&]&/@loc];
	loc = Append[loc, _];
	
	(*get segments*)
	segs = Which[
		segmi==="AHA",{1,4,6,6,0},
		segmi==="AHA+",{2,6,8,8,0},
		MemberQ[{1,2,4,6,8},segmi],{1,1,1,1,0}segmi,
		True, $Failed
	];
	
	(*make the segmentations*)
	segs = Thread[loc->segs];
	If[!StringQ[segmi],segs = DeleteCases[segs,{}->_]];
	If[grp, segs, Flatten[Thread/@segs]]
]


(* ::Subsection:: *)
(*MaskToLines*)


(* ::Subsubsection::Closed:: *)
(*MaskToLines*)


SyntaxInformation[MaskToLines] = {"ArgumentsPattern" -> {_, _, _.}};

MaskToLines[mask_, vox_] := Block[{wall, off}, 
	off = CentralAxes[mask, vox, ShowPlot -> False][[1]]; 
	wall = CalculateWallMap[mask, vox, MaskWallMap -> False, ShowPlot -> False][[1]]; 
	MaskToLines[mask, wall, off]
]

MaskToLines[mask_?ArrayQ, wall_?ArrayQ, cent_?ListQ] := Block[{
	mid, offi, maski, midi, walli, dim, dimx, dimy, ptsmi, ptsm, angles,
	lines, tresh, seg, min, mout, mmid, int, int1, int2, intDer, step, centVal,
	ord, vec, ptsol, ptso, ptsil, ptsi, ptsml, x, y, val, vali, k, dis, 
	ii, pp},
	
	mid = Mask[wall, {0, .5}];
	(*Print[Dimensions[mask]];*)
	Table[
		(*Print[i];*)(*values that need indexing*)
		offi = cent[[i, {2, 3}]];
			
		maski = Round @ mask[[i]];
		midi = mid[[i]];
		walli = wall[[i]];
		(*static values*)
		dim = {dimx, dimy} = Dimensions[maski];
		tresh = 0.3;
		(*points if no mask*)
		ptsmi = Table[RotationMatrix[-i Degree] . {0., 1.}, {i, 0, 359, 1.}];
		If[Total[Flatten[maski]] < 5,
			
			(*if no maks do nothing*)
			(*Print["no Mask"];*)
			ptsm = ptsmi; lines = ConstantArray[{}, 360];
			,
			(*get the inner and/or outer and/or mid mask*)
			seg = If[#[[1, 1]] == 2, # /. {2 -> 1, 1 -> 2}, #]&[MorphologicalComponents[N[1 - maski]]];
			min = Unitize[Clip[seg, {1.5, 2}, {0, 0}]];
			mout = Unitize[1 - Clip[seg, {0, 1}, {0, 0}]];
			seg = If[#[[1, 1]] == 2, # /. {2 -> 1, 1 -> 2}, #]&[MorphologicalComponents[N[1 - midi]]];
			mmid = Unitize[1 - Clip[seg, {0, 1}, {0, 0}]];
			
			(*check which method to use*)
			(*1. both inner and outer mask*)
			(*2. apex,only other mask with cent in the mask*)
			(*3. base,no closed maks with cent outside the mask*)
			If[Total[Flatten[min]] =!= 0,
				
				(*Print["inner and outer mask"];*)
				(*there is an inner and outer mask-here lines are formed perpendicular to the mid wall derivative*)
				
				(*get outer and inner points-720 points equally distributed radially*)
				{ptsol, ptso} = PerimiterPoints[(*Dilation[mout,1]*)mout, dim, .5];
				{ptsil, ptsi} = PerimiterPoints[Dilation[min, 1], dim, .5];
				
				(*get mid points-360 points equally distributed radially*)
				{ptsml, ptsm} = PerimiterPoints[mmid, dim, 1];
				dis = Flatten[(pp = #; EuclideanDistance[pp, #]& /@ Nearest[ptsol, pp, 2])& /@ ptsil];
				k = Floor[.8 Mean[dis]];
				dis = .4 Min[dis];
				
				(*get the wall derivatives*)
				(*make the derivative funciton-vector pointing perpendicularto the wall at wall mid point*)
				{int1, int2} = ListInterpolation[GaussianFilter[walli, {k, k}, #], InterpolationOrder -> 1]& /@ IdentityMatrix[2];
				intDer = {int1[#[[1]], #[[2]]], int2[#[[1]], #[[2]]]}&;
				
				(*projects line in and outward to closest wall point*)
				lines = (step = dis Normalize[intDer[#]]; {FindPoint[#, -step, Nearest[ptsi]], FindPoint[#, step, Nearest[ptso]]})& /@ ptsm;
				,
				
				(*there is only an outer mask*)
				centVal = maski[[Round[offi[[1]]], Round[offi[[2]]]]];
				
				(*see if center is in mask*)
				If[centVal === 1,
					(*Print["center in the mask"];*)
					(*the center falls within the mask-apex region*)
					{ptsml, ptsm} = PerimiterPoints[mout, dim, 1];
					
					(*get the line points and angles*)
					lines = Thread[{ConstantArray[offi, Length[ptsm]], ptsm}];
					,
					
					(*Print["center outside the mask"];*)
					(*the center falls outside the mask*)
					step = 0.75;
					ptsm = ptsmi;
					int = ListInterpolation[maski, InterpolationOrder -> 1];
					{x, y} = offi;
					vali = int[x, y];
					
					(*get the line points and angles*)
					lines = Map[(
						vec = #;
						ii = 0;
						val = vali;
						{x, y} = offi;
						
						(*find the first position for which point is inside the mask*)	
						While[(1 + step <= y <= dimy - step - 1 && 1 + step <= x <= dimx - step - 1 && val <= 0.), 
							{x, y} = offi + (ii step) vec; 
							val = int[x, y]; ii++;
						];
						ptsi = offi + ((ii - 1) step) vec;
						
						(*Then find the first position for which point is outside the																																	 																mask again*)
						While[(1 + step <= y <= dimy - step - 1 && 1 + step <= x <= dimx - step - 1 && val > 0.), 
							{x, y} = offi + (ii step) vec; 
							val = int[x, y]; ii++;
						];
						ptso = offi + ((ii - 2) step) vec;
						
						(*if none found return {} else return cent and firt point*)
						If[val > tresh, {ptsi, ptso}, {}]
					)&, ptsm];
				]
			]
		];
		(*put everythin in order*)
		lines
	,{i, 1, Length[mask]}]
]


(* ::Subsubsection::Closed:: *)
(*PerimiterPoints*)


PerimiterPoints[im_, {x_, _}, deg_] := Block[
	{ptsp, pts, len, int}
	,
	(*get perimiter points and tranfer to array index coordinates *)
	ptsp = Reverse /@ ComponentMeasurements[Image[im], "PerimeterPositions"][[1, 2, 1]];
	ptsp = Transpose[{Floor @ ptsp[[All, 1]], Ceiling @ ptsp[[All, 2]]}];
	ptsp = N[{x, 0} + {-1, 1} #& /@ ptsp];
	(*interpolate perimiter points 360 equally spaced points*)
	pts = Append[ptsp, First[ptsp]];
	len = Prepend[Accumulate[ArcLength[Line[#]]& /@ Partition[pts, 2, 1]], 0];
	int = Interpolation[Thread[{len, pts}], InterpolationOrder -> 1];
	pts = int /@ Range[0, Max[len], Max[len] / ((360 / deg) - 1)];
	(*output*)
	{ptsp, pts}
]


(* ::Subsubsection::Closed:: *)
(*FindPoint*)


FindPoint[pmid_, step_, fun_] := Block[
	{ptar, pnear, ni, nd, n, fnear, v, u, p, j}
	,
	(*define line perpendicular to wall and find nearest wall point*)
	ptar = pmid + step;
	pnear = First @ fun[ptar, 1];
	ni = nd = Sqrt[Total[(pnear - ptar) ^ 2]];
	j = 0;
	(*project point on line and find the new point closest to new line*)
		
	While[nd > 10 ^ -3 && j < 10, 
		j++; 
		v = ptar - pmid; 
		u = pnear - pmid; 
		p = (v . u / v . v) v; 
		ptar = pmid + Sign[p . step] p; 
		pnear = First @ fun[ptar, 1]; n = Sqrt[Total[(pnear - ptar) ^ 2]]; 
		nd = Abs[n - ni]; 
		ni = n;
	];
	pnear
]


(* ::Subsection:: *)
(*LinesToSegmentIndex*)


(* ::Subsubsection::Closed:: *)
(*LinesToSegmentIndex*)


Options[LinesToSegmentIndex]= {ReversePoints->True, ReverseDirection->False};

SyntaxInformation[LinesToSegmentIndex] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

LinesToSegmentIndex[lines_, pts_, seg_, OptionsPattern[]]:=Block[{outLine,ptIndex,segIndex,lineIndex,flip,rev,segI},
	flip = OptionValue[ReversePoints];
	rev = If[OptionValue[ReverseDirection], -1, 1];
	segI = Flatten[Thread/@seg];
	
	{ptIndex, segIndex, lineIndex} = Transpose@Table[
	outLine = DeleteCases[lines[[i]], {}];
	If[Length[outLine]=!=360||pts[[i,1]]===None||pts[[i,2]]===None,
		{None, None, None}
		,
		outLine = outLine[[All,2]];
		ptIndex = FindLine[outLine,lines[[i]],pts[[i]]];
		segIndex = SegmentIndex[i/.segI,If[flip,Reverse@ptIndex,ptIndex],rev];
		lineIndex = LinesIndex[segIndex,rev];
		{ptIndex, segIndex, lineIndex}
	],{i,1,Length@lines}]
]


(* ::Subsubsection::Closed:: *)
(*FindLine*)


FindLine[line_,lines_,pts:{{_,_},{_,_}}]:=FindLine[line,lines,#]&/@pts

FindLine[line_,lines_,pt:{_,_}]:=Block[{ind,dist},
	ind=Sort@DeleteDuplicates@Flatten[Position[line,#]&/@Nearest[line,pt,10]];
	dist=DistanceToLine[lines[[#]],pt]&/@ind;
	ind[[Position[dist,Min@dist][[1,1]]]]
]


(* ::Subsubsection::Closed:: *)
(*DistanceToLine*)


DistanceToLine[{a1_,a2_},p_]:=Block[{v,u,pp},
	If[a1===a2,
		Nothing,
		(*vectors to project point on line*)
		v=a2-a1;
		u=p-a1;
		(*distance of point to projection*)
		EuclideanDistance[p,a1+(v . u/v . v)v]
	]
]


(* ::Subsubsection::Closed:: *)
(*SegmentIndex*)


SegmentIndex[num_,{p1_,p2_},rev_]:=Block[{dwall,dsept,inds},
	(*index differences between points over the free wall*)
	dwall=WrapIndex[If[rev==1,p1-p2,p2-p1]];
	(*index differences between points over the septum*)
	dsept=360-dwall;
	
	(*get all the possible step sizes for the septum (2,3) or free wall (3,4,5)*)
	inds=Transpose[Round@{p1,p2,1./2  dsept,1./4  dwall,
	1./3  dwall,1./3  dsept,1./5  dwall}/.Round[None]->None];
	
	(*based on number of segments get the segment start indexes*)
	WrapIndex/@Switch[num,
		0,{None},
		1,{inds[[1]]},
		2,{inds[[1]],inds[[2]]},
		4,{inds[[1]]-(rev inds[[5]]),inds[[1]],inds[[2]],inds[[2]]+(rev inds[[5]])},
		6,{inds[[1]]-(rev inds[[4]]),inds[[1]],inds[[1]]+(rev inds[[3]]),inds[[2]],inds[[2]]+(rev inds[[4]]),inds[[2]]+(rev 2 inds[[4]])},8,{inds[[1]]-(rev 2 inds[[7]]),inds[[1]]-(rev inds[[7]]),inds[[1]],inds[[1]]+(rev inds[[6]]),inds[[1]]+(rev 2 inds[[6]]),inds[[2]],inds[[2]]+(rev inds[[7]]),inds[[2]]+(rev 2 inds[[7]])},_,{inds[[1]]}
	]
]


(* ::Subsubsection::Closed:: *)
(*WrapIndex*)


WrapIndex[ind_?ListQ]:=WrapIndex/@ind
WrapIndex[None]:=None
WrapIndex[ind_]:=If[ind<1,ind+360,If[ind>360,ind-360,ind]]


(* ::Subsubsection::Closed:: *)
(*LinesdIndex*)


LinesIndex[inds_,rev_]:=Block[{i1,i2},
	(*convert indexes to lists*)
	If[Length[inds]===1,
		(*if only 1 index None or full circle*)
		If[inds[[1]]===None,{},{Range[360]}]
		,(
		{i1,i2}=#;
		(*depending of direction range between indexes*)
		If[Sign[i2-i1]===rev,
			Range[i1,i2,rev],
			If[rev===1,
			Join[Range[i1,360,rev],Range[1,i2,rev]],
			Join[Range[i1,1,rev],Range[360,i2,rev]]
			]
		])&/@Partition[Join[inds,{inds[[1]]}],2,1]
	]
]


(* ::Subsection:: *)
(*GetSegmentLines*)


(* ::Subsubsection::Closed:: *)
(*GetSegmentLines*)


SyntaxInformation[GetSegmentLines] = {"ArgumentsPattern" -> {_, _, _}};

GetSegmentLines[lines_, lineIndex_, seg_]:=Block[{sl,ss},
	Table[(
		{sl,ss}=#;
		If[Length[lines[[sl]]]===360&&lineIndex[[sl]]=!=None,
		{sl,lines[[sl,lineIndex[[sl,ss]]]]},
		{sl,None}
		]
	)&/@loc,{loc, SegmentLocations[seg]}]
];


(* ::Subsubsection::Closed:: *)
(*SegmentLocations*)


SegmentLocations[segments_]:=Block[{location,trans,sl,ind},
	(*get the locations*)
	location = If[segments=!=$Failed,
		location = If[Head[#]===List,
			trans = True;
			(
				sl=#[[1]];
				ind=#[[2]];
				{sl,#}&/@Range[ind]
			)&/@#,
			trans=False;
			sl=#[[1]];
			ind=#[[2]];
			{sl,#}&/@Range[ind]
		]&/@Select[(Thread/@segments),Head[First[#]]=!=Blank&];
		If[trans,Flatten[Transpose/@location,1],Transpose[location]]
	];
	Reverse@location
]


(* ::Subsection::Closed:: *)
(*SegmentLinesToMask*)


SyntaxInformation[SegmentLinesToMask] = {"ArgumentsPattern" -> {_, _}};

SegmentLinesToMask[smsk_, segLines_]:=Block[{out,ran,segMask,sl,ln,pol,msk},
	out=0 smsk;
	ran=Reverse@Dimensions[smsk][[2;;]];
	
	segMask=Table[
		out=0 smsk;
		(
			{sl,ln}=#;
			If[ln=!=None,
				pol=Polygon@Join[DeleteDuplicates[RevRound/@ln[[All,2]]],
				Reverse@DeleteDuplicates[RevRound/@ln[[All,1]]]];
				msk=Reverse@Closing[ImageData@Binarize[1-Rasterize[Graphics[pol,PlotRange->Thread@{1,ran}],RasterSize->ran]],1];
				out[[sl]]=Round[smsk[[sl]] msk];
			]
		)&/@segment;
		out
	,{segment,segLines}];
	MergeSegmentations[RemoveMaskOverlaps[Transpose@segMask],Range[Length@segMask]]
]

RevRound={Floor[#[[2]]],Ceiling[#[[1]]]}&;


(* ::Subsection:: *)
(*MakeImages*)


(* ::Subsubsection::Closed:: *)
(*MakeLineImage*)


SyntaxInformation[MakeLineImage]={"ArgumentsPattern"->{_,_,_}};

MakeLineImage[back_,segLines_,pts_]:=Block[{colRule,pos,slines,p1,p2},
	colRule=Thread[Range[17]->(SeedRandom[113];RandomSample[Ncol[17]])];
	ImageCollage[Table[
		pos=Position[segLines[[All,All,1]],i];
		If[pos=!={},
			{p1,p2}=Switch[Length@pos,4,{2,3},6,{2,4}+1,_,{None,None}];
			slines=segLines[[First@#,Last@#,2]]&/@pos;
			Show[
				ArrayPlot[N[Reverse@back[[i]]],ColorRules->{0.->White,1.->Gray,2.->Lighter@Lighter@Red,3.->Lighter@Lighter@Green},ImageSize->300,PlotLabel->Style[i,Large],Frame->False],
				Graphics[{PointSize[Large],Red,Point@Reverse@pts[[i,1]],Green,Point@Reverse@pts[[i,2]]}],
				If[AnyTrue[slines,#===None&],
					Graphics[],
					Show[
						If[Length@slines=!=1,Graphics[{Thickness[.025],Darker@Gray,Line[RevShift/@#&/@slines[[All,1]]]}],Graphics[]],
						Graphics[Thread[{pos[[All,1]]/.colRule,Map[Line,Map[RevShift,slines[[All,;;;;3]],{3}],{2}]}]],
						If[p1=!=None,Graphics[{{Red,Line[RevShift/@slines[[p1,1]]]},{Green,Line[RevShift/@slines[[p2,1]]]}}],Graphics[]]
					]
				]
			],
			Nothing
		]
	,{i,1,Length[back]}],Background->White]
]


RevShift=Reverse[#1-{0.5,0.5}]&;


(* ::Subsubsection::Closed:: *)
(*MakeMaskImage*)


SyntaxInformation[MakeMaskImage]={"ArgumentsPattern"->{_,_}};

MakeMaskImage[back_,mask_]:=ImageCollage[Table[
	If[DeleteDuplicates[Flatten[mask[[i]]]]==={0},
		Nothing,
		Show[
			ArrayPlot[N[Reverse@back[[i]]],ColorRules->{0.->White,1.->Gray,2.->Lighter@Lighter@Red,3.->Lighter@Lighter@Green},ImageSize->300,PlotLabel->Style[i,Black,Large],Frame->False],
			ArrayPlot[Reverse@mask[[i]],ColorRules->Join[{0->Transparent},Thread[Range[17]->(SeedRandom[113];RandomSample[Ncol[17]])]]]
		]
	]
,{i,1,Length[mask]}],Background->White]


(* ::Subsubsection::Closed:: *)
(*Ncol*)


Ncol=If[#===1,{ColorData["DarkRainbow"][0.]},ColorData["DarkRainbow"]/@(Range[0,#-1]/(#-1))]&;


(* ::Subsection::Closed:: *)
(*CardiacSegment*)


Options[CardiacSegment] = Options[CardiacSegmenti] = {ReversePoints->True, ReverseDirection->False, MakeSegmentPlots->True}

SyntaxInformation[CardiacSegment]={"ArgumentsPattern"->{_,_,_,_.,_.,OptionsPattern[]}};

CardiacSegment[mask_?ArrayQ,vox:{_?NumberQ,_?NumberQ,_?NumberQ}, pts_, ops:OptionsPattern[]]:=CardiacSegmenti[mask, mask, vox, pts, SegmentsPerSlice[pts], ops]

CardiacSegment[mask_?ArrayQ,vox:{_?NumberQ,_?NumberQ,_?NumberQ}, pts_, seg_?VectorQ, ops:OptionsPattern[]]:=CardiacSegmenti[mask, mask, vox, pts, SegmentsPerSlice[seg],  ops]

CardiacSegment[mask_?ArrayQ, back_?ArrayQ, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, pts_, ops:OptionsPattern[]]:=CardiacSegmenti[mask, back, vox, pts, SegmentsPerSlice[pts], ops]

CardiacSegment[mask_?ArrayQ, back_?ArrayQ, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, pts_, seg_?VectorQ, ops:OptionsPattern[]]:=CardiacSegmenti[mask, back, vox, pts, SegmentsPerSlice[seg], ops]

CardiacSegmenti[mask_?ArrayQ, back_?ArrayQ, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, pts_, seg_, OptionsPattern[]]:=Block[{
		lines, ptIndex, segIndex, lineIndex, segLines, segMask, plotLines, plotMask
	},
	
	Print[seg];
	
	lines = MaskToLines[mask, vox];
	{ptIndex, segIndex, lineIndex} = LinesToSegmentIndex[lines, pts, seg, 
		ReversePoints->OptionValue[ReversePoints], ReverseDirection->OptionValue[ReverseDirection]];
	segLines = GetSegmentLines[lines, lineIndex, seg];
	segMask = SegmentLinesToMask[mask, segLines];
	
	If[OptionValue[MakeSegmentPlots],
		plotMask = MakeMaskImage[back, segMask];
		plotLines = MakeLineImage[back, segLines, pts];
		{segLines, segMask, {plotLines, plotMask}},
		{segLines, segMask}
	]
]


(* ::Subsection:: *)
(*CardiacSegmentGUI*)


Options[CardiacSegmentGUI]={StartPoints->"Default", StartSlices->"Default", ReversePoints->True, ReverseDirection->False};

SyntaxInformation[CardiacSegmentGUI]={"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};

CardiacSegmentGUI[data_, maski_, vox_, OptionsPattern[]]:=DialogInput[{
	DynamicModule[{
			off, mask, radiusStart, centers, startPoint, defPoints, revi, flipi, pointsIn, centerpl,
			slices, sti, api, midi, basi, endi, carPl, colsOver, car,lm, allpli, radi
		},
		
		mask=Round[maski];
		radiusStart=ConstantArray[radi=Max[Dimensions[data]]/8,Length[data]];
		
		off = CentralAxes[mask, vox, ShowPlot -> False][[1]];
		centers=Reverse/@off[[All,2;;3]];
		
		startPoint=OptionValue[StartPoints];
		
		defPoints=Transpose[{centers+Transpose@{radiusStart,radiusStart},centers-Transpose@{radiusStart,radiusStart}},{2,1,3}];
		pointsIn=If[startPoint==="Default",defPoints,Rev@startPoint];
		revi=OptionValue[ReverseDirection];
		flipi=OptionValue[ReversePoints];
		
		centerpl=Graphics[{Red,Disk[#,1]}]&/@centers;
		slices=Length[data];
		
		{sti,api,midi,basi,endi}=If[OptionValue[StartSlices]==="Default"||OptionValue[StartSlices]==={},GetSegmentSlices[mask],OptionValue[StartSlices]];
		
		carPl=ArrayPlot[car=Reverse@Total[mask,{2}],ColorFunction->"GrayTones",Frame->False,ImageSize->{200,200},AspectRatio->1];
		colsOver=0 car;
		lm=Dimensions[car];
		allpli=Graphics[MapThread[{#2,Thick,Line[{{0,#1},{lm[[2]],#1}}]}&,{{sti,api,midi,basi,endi},{Gray,Orange,Blue,Red,Purple}}]];
		
		Manipulate[
			n = Clip[n,{1,slices}];
			rev = If[revb,-1,1];
			
			(*slice sements*)
			end = Clip[end,{bas,slices}];
			bas = Clip[bas,{mid,end}];
			ap = Clip[ap,{start,mid}];
			mid = Clip[mid,{ap,bas}];
			start = Clip[start,{0,ap}];
			
			allpl = Graphics[MapThread[{#2,Thick,Line[{{0,#1},{lm[[2]],#1}}]}&,{{start,ap,mid,bas,end},{Gray,Orange,Blue,Red,Purple}}]];
			
			(*segments*)
			segments = SegmentsPerSlice[{start,ap,mid,bas,end}, SegmentationMethod->segmi, GroupPerSegment->slcgrp];
			segm = n/.Flatten[Thread/@segments];
			
			(*Colors*)
			backcols=Flatten[If[Head[#[[1]]]===List,Thread[#],#]&/@Thread[segments[[All,1]]->{Orange,Blue,Red,Purple,Gray}]];
			colls=backcols[[All,2]];
			backcol=n/.backcols;
			
			(*plots*)
			cent=centers[[n]];
			ptSel=(points[[n,1]]=!=None&&points[[n,2]]=!=None);
			If[!ptSel&&segm>0,points[[n]]=defPoints[[n]];ptsSel=True];
			
			If[ptSel,
				rad=Clip[RadCalc[points[[n]],cent],{5,Min[Dimensions[mask[[1]]]/2]}];
				angs=VecAngleC[If[flip,Reverse@points[[n]],points[[n]]],cent,segm,rev];
				
				dsks=RegionDisk[cent,angs,rad,rev];
				polplot=Graphics[Table[{Opacity[0.5],ColFun[i,segm,backcol],dsks[[i]]},{i,1,segm}]];
				
				arplot=Switch[segm,
					0,Graphics[],
					_,
					segpt=(RotateRad[rad,cent,#]&/@angs);
					Graphics[{Thick,backcol,Arrow[{cent,#}]&/@segpt}]
				];
				,
				rad=radi;
				polplot=arplot=Graphics[];
			];
			
			anatplot=ArrayPlot[Reverse[data[[n]]],ColorFunction->"GrayTones",FrameTicks->Automatic,ImageSize->400, 
				LabelStyle->Directive[{Bold,Black,FontFamily->"Helvetica"}],FrameTicksStyle->Directive[{Black,Thick}],
				FrameStyle->Directive[{Thick,backcol}]
			];
			maskplot=ArrayPlot[Reverse[1-mask[[n]]],ColorRules->{1->Transparent,0->GrayLevel[1,mop]}];
			circplot=Graphics[{backcol,Thick,Circle[cent,rad]}];
			
			(*locator pane*)
			Column[{
				(*sclice buttons*)
				Row[{
					Button["<<<",n=1],
					Button["<<",n=start+1],
					Button["<",n=n-1],
					Button[">",n=n+1],
					Button[">>",n=end],
					Button[">>>",n=slices]
				}],
				(*locator pane*)
				LocatorPane[
					(*dynamic points*)
					If[segm===0,None,Dynamic[points[[n]]]],
					(*dynamic plot*)
					Dynamic[Show[anatplot, maskplot, circplot, polplot, arplot, centerpl[[n]], ImageSize->600]]
					,
					(*locator appearance*)
					Appearance->If[segm===0,None,{Graphics[{Green,Disk[]},ImageSize->15],Graphics[{Blue,Disk[]},ImageSize->15]}]
				]
			},Alignment->Center]
			
			,
			(*manipulate controls*)
			{{n,Round[slices/2],"Slice"},1,slices,1},
			{{mop,.2,"Mask opacity"},0,1},
			
			Delimiter,
			{{segmi,"AHA","Number of segments"},{1->"1 per sliec",2->"2 per sliec",4->"4 per slice",6->"6 per slice",8->"8 per slice","AHA"->"AHA-17","AHA+"->"AHA-24"}},
			{{revb,revi,"Reverse direction"},{True,False}},
			{{flip,flipi,"Flip points"},{True,False}},
			
			Delimiter,
			(*slice segmentation*)
			Row[{
				VerticalSlider[Dynamic[start],{0,slices,1},Background->Lighter@Gray,Appearance->{Vertical,Tiny}],"  ",
				VerticalSlider[Dynamic[ap],{0,slices,1},Background->Lighter@Orange,Appearance->{Vertical,Tiny}],"  ",
				VerticalSlider[Dynamic[mid],{0,slices,1},Background->Lighter@Blue,Appearance->{Vertical,Tiny}],"  ",
				VerticalSlider[Dynamic[bas],{0,slices,1},Background->Lighter@Red,Appearance->{Vertical,Tiny}],"  ",
				VerticalSlider[Dynamic[end],{0,slices,1},Background->Lighter@Purple,Appearance->{Vertical,Tiny}],"      ",
				Dynamic[Show[carPl,allpl]]
			}],
			
			(*close buttons*)
			Delimiter,
			Row[{
				DefaultButton["Done",DialogReturn[
					out=ConstantArray[{None,None},Length[points]];
					sl=Range[start+1,end];
					out[[sl]]=points[[sl]];
					{(Rev@out), {start,ap,mid,bas,end}, {revb, flip}}
				]],
				CancelButton["Cancel",DialogReturn[$Canceled]]
			}],
			
			(*hidden controls*)
			{{ap,api},ControlType->None},
			{{mid,midi},ControlType->None},
			{{bas,basi},ControlType->None},
			{{end,endi},ControlType->None},
			{{start,sti},ControlType->None},
			{{allpl,allpli},ControlType->None},
			
			{{points,pointsIn},ControlType->None},
			{{rev,revi},ControlType->None},
			
			{ptsSel,ControlType->None},
			{rad,ControlType->None},
			{cent,ControlType->None},
			{angs,ControlType->None},
			{segpt,ControlType->None},
			{dsks,ControlType->None},
			
			{segments,ControlType->None},
			{segm,ControlType->None},
			{{slcgrp,True},ControlType->None},
			{backcols,ControlType->None},
			{backcol,ControlType->None},
			{colls,ControlType->None},
			
			{anatplot,ControlType->None},
			{maskplot,ControlType->None},
			{circplot,ControlType->None},
			{polplot,ControlType->None},
			{arplot,ControlType->None},
			
			{out,ControlType->None},
			{sl,ControlType->None},
			{ptSel,ControlType->None},
			
			
			(*initialization*)
			ControlPlacement->Left,
			SynchronousUpdating->True,
			SynchronousInitialization->False
		](*close Manipulate*)
	](*close dynamic module*)
(*close dialog input*)
},WindowTitle->"Segement the heart"(*,WindowFloating->True,Modal->False*)];


Rev[pts_]:=Map[If[#=!=None,Reverse@#,#]&,pts,{2}]


RadCalc[pnt_, cent_]:=Mean[Norm[#-cent]&/@pnt];


VecAngleC=Compile[{{pnt,_Real,2},{cent,_Real,1},{num,_Integer,0},{rev,_Real,0}},Block[{
		x1,x2,y1,y2,ang1,ang2,angs,angs2,angsout,diff,dsep,dwall
	},
	{{x1,y1},{x2,y2}}=(#-cent&/@pnt);
	ang1=ArcTan[y1,x1];
	ang2=ArcTan[y2,x2];
	
	diff=If[rev==1,ang2-ang1,ang1-ang2];
	diff=If[Negative[diff],2Pi+diff,diff];
	{dsep,dwall}={diff,2Pi-diff};
	
	angs = {ang1,ang2,1./2 dsep,1./4 dwall,1./3 dwall,1./3 dsep,1./5 dwall};
	
	(*define angels based on number of segments*)
	angs2=Switch[num,
		0,{0.123},
		1,{angs[[1]]},
		2,{angs[[1]],angs[[2]]},
		4,{angs[[1]]-(rev angs[[5]]),angs[[1]],angs[[2]],angs[[2]]+(rev angs[[5]])},
		6,{angs[[1]]-(rev angs[[4]]),angs[[1]],angs[[1]]+(rev angs[[3]]),angs[[2]],angs[[2]]+(rev angs[[4]]),angs[[2]]+(rev 2 angs[[4]])},
		8,{angs[[1]]-(rev 2 angs[[7]]),angs[[1]]-(rev angs[[7]]),angs[[1]],angs[[1]]+(rev angs[[6]]),angs[[1]]+(rev 2 angs[[6]]),angs[[2]],angs[[2]]+(rev angs[[7]]),angs[[2]]+(rev 2 angs[[7]])},
		_,{angs[[1]]}
	];

	(*export angles that are always between 0 and 2Pi*)
	If[Negative[#],2Pi+#,#]&/@angs2
]];


RotateRad[_,_,0.123]={};

RotateRad[_,_,{}]={};

RotateRad[rad_,cent_,ang_]:=(cent+{{Cos[ang],Sin[ang]},{-Sin[ang],Cos[ang]}} . {0,rad});


RegionDisk[cent_,ang_,ran_,rev_]:=Block[{ans},
	ans=Partition[If[Negative[#],#+2Pi,#]&/@(-ang+0.5Pi),2,1,1];
	If[rev==1,
		Disk[cent,ran,If[#[[1]]>#[[2]],#,#+{0,-2Pi}]]&/@ans,
		Disk[cent,ran,If[#[[2]]>#[[1]],Reverse@#,Reverse[#+{0, 2Pi}]]]&/@ans
	]
]


ColFun[i_,seg_,back_]:=If[back===Gray,back,If[seg==1,ColorData["Rainbow"][0.],ColorData["Rainbow"][((i-1)/(seg-1.))]]]

ColFun[i_,seg_]:=If[seg==1,ColorData["Rainbow"][0.],ColorData["Rainbow"][((i-1)/(seg-1.))]]


(* ::Subsection:: *)
(*PlotSegments*)


(* ::Subsubsection::Closed:: *)
(*PlotSegments*)


Options[PlotSegments] = {RadialSamples -> 10};

SyntaxInformation[PlotSegments] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

PlotSegments[data_, mask_, angs_, OptionsPattern[]] := Block[{pan}, 
 	pan = Manipulate[ 		
 		slices=angs[[n,All,1]];
 		pts=DeleteCases[#,{}]&/@angs[[n,All,2]];
 		nmr=Range[Length[slices]];
 		
 		datpl=ArrayPlot[data[[#]], DataReversed -> True, ColorFunction -> "GrayTones"]&/@slices;
 		maskpl=mask[[slices]];
 		
 		size = Length[angs[[n]]] 250;
 		(*switch between views*)
 		GraphicsRow[Show[
 		Switch[m,
 			(*mask with inner and outer pionts only*) 
 			1, {datpl[[#]],
 				If[NumberQ[Mean[Flatten[pts[[#]]]]],
 					ListPlot[{pts[[#,All,1]], pts[[#,All,2]]}, AspectRatio -> 1, PlotStyle -> {Red, Orange}],
 					Graphics[]
 				]},
 			(*mask with inner and outer points plus radial samples*)
 			2, {datpl[[#]],
 				If[NumberQ[Mean[Flatten[pts[[#]]]]],
 					ListPlot[Transpose@PointRange[pts[[#]], OptionValue[RadialSamples]], AspectRatio -> 1],
 					Graphics[]
 				]},
 			3, {
 				If[NumberQ[Mean[Flatten[pts[[#]]]]],
 					ArrayPlot[Fun[maskpl[[#]], Round[pts[[#,All,1]]], Round[pts[[#,All,2]]]], DataReversed -> True, ColorFunction -> "GrayTones"],
 					ArrayPlot[maskpl[[#]], DataReversed -> True, ColorFunction -> "GrayTones"]
 				]}
     		]
 		] &/@ nmr, ImageSize->size], 
	    {{n, 1, "segment"}, 1, Length[angs], 1}, 
	    {{m, 1, "plot type"}, {1 -> "start stop point", 2 -> "radial samples", 3 -> "mask"}}, 
	    {size, ControlType -> None},
	    {pts, ControlType -> None},
	    {nmr, ControlType -> None},
	    {slices, ControlType -> None},
	    {datpl, ControlType -> None},
	    {maskpl, ControlType -> None},
	    
	    SaveDefinitions->True
     ];
  NotebookClose[plotwindow];
  plotwindow = CreateWindow[DialogNotebook[{CancelButton["Close", DialogReturn[]], pan}, WindowSize -> All, WindowTitle -> "Plot data window"]];
  ]


(* ::Subsubsection::Closed:: *)
(*Fun*)


Fun[mask_,pts1_,pts2_]:=Block[{tmp=mask},
	(tmp[[#[[2]],#[[1]]]]=10)&/@pts1;
	(tmp[[#[[2]],#[[1]]]]=-10)&/@pts2;
	tmp]


(* ::Subsubsection::Closed:: *)
(*PointRange*)


PointRange[pts_,steps_]:=Block[{step=steps-1,pt1,pt2},
	(
	pt1=#[[2,{1,2}]];
	pt2=#[[1,{1,2}]];
	Table[pt2+(i (pt1-pt2)/step),{i,0,step,1}]
	)&/@pts
]

PointRange[pts_,steps_,drop_]:=Block[{step=steps-1,pt1,pt2},
	(
	pt1=#[[2,{1,2}]];
	pt2=#[[1,{1,2}]];
	Take[Table[pt2+(i (pt1-pt2)/step),{i,0,step,1}],{1+drop,steps-drop}]
	)&/@pts
]


(* ::Subsection::Closed:: *)
(*PlotSegmentsMask*)


SyntaxInformation[PlotSegmentMask] = {"ArgumentsPattern" -> {_, _, _}};

PlotSegmentMask[maski_, segmaski_, vox_] := Block[{heart, seg,pan},
	heart = PlotMaskVolume[maski, vox];
	seg = PlotMaskVolume[segmaski[[#]], vox, Red,Filter->False] & /@ Range[Length[segmaski]];
	pan=Manipulate[
		GraphicsGrid[{{
			Show[heart, seg[[n]]],
			Show[heart, seg[[n]], ViewPoint -> Front, Method -> {"RotationControl" -> None}]
			},{
			Show[heart, seg[[n]], ViewPoint -> Top, Method -> {"RotationControl" -> None}],
			Show[heart, seg[[n]], ViewPoint -> Left, Method -> {"RotationControl" -> None}]
		}}, ImageSize -> 600]
   , {{n,1,"Segment"}, 1, Length[segmaski], 1},SaveDefinitions->True];
   
   NotebookClose[plotwindow];
   plotwindow = CreateWindow[DialogNotebook[{CancelButton["Close", DialogReturn[]], pan}, WindowSize -> All, WindowTitle -> "Plot data window"]];
]


(* ::Subsection::Closed:: *)
(*RadialSample*)


Options[RadialSample]={RadialSamples->10, DropSamples->0};

SyntaxInformation[RadialSample] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

RadialSample[mask_,data_,segang_,OptionsPattern[]]:=Block[{
	pos,dat,val,intdat,exdat,intfunc, output,slice,pts,int,ptsr,vals
	},

	(*slice by slice interpolation function*)
	intfunc = MapThread[(
		pos = Position[Round[#1], 1];
		dat = #2;
		val = dat[[#[[1]], #[[2]]]] & /@ pos;
		If[pos =!= {},
			intdat = N@Thread[{pos, val}];
			exdat = N@Thread[pos -> val];
			With[{nearF = Nearest[exdat]}, Interpolation[intdat, InterpolationOrder -> 1,
				"ExtrapolationHandler" -> {Mean[nearF[{#1, #2}, 3]] &, "WarningMessage" -> False}]
			],
			0 &
		]
	) &, {mask, data}, 1];
	
	output = Map[(
		slice = #[[1]];
		pts = #[[2]];
		int = intfunc[[slice]];
		
		If[pts=!=None,
			(*get radial pts steps*)
			ptsr = PointRangeC[pts, OptionValue[RadialSamples],OptionValue[DropSamples]];
			vals = int @@ Transpose[ptsr, {2, 3, 1}];
			{ptsr, vals},
			{{}, {}}
		]
	) &, segang, {2}];
	
	{output[[All,All,1]],output[[All,All,2]]}
]

PointRangeC = Compile[{{pts, _Real, 2}, {steps, _Integer, 0}, {drop, _Integer, 0}}, Block[{n, step},
	n = steps - 1;
	step = (pts[[2]] - pts[[1]])/n;
	pts[[1]] + # step & /@ Range[0 + drop, n - drop]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection::Closed:: *)
(*TransmuralPlot*)


Options[TransmuralPlot] = {GridLineSpacing -> 10, PlotStyle -> Red, PlotRange -> Automatic, ImageSize->300, Method->"Median", PlotLabel->None}

SyntaxInformation[TransmuralPlot] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

TransmuralPlot[data_, OptionsPattern[]] := Block[{mn, std, steps, stdM, stdP, xdata, min, max, col,pdat,dat,fil,style},
	dat=If[ArrayDepth[data]==1,
		{data},
		Switch[OptionValue[Method],
			"Median", Transpose[Quantile[data, {.5, .35, .65}]],
			"MedianQ", Transpose[Quantile[data, {.5, .25, .75}]],
			"MedianSD", Transpose[Quantile[data, {.5, .16, .84}]],
			"Median95", Transpose[Quantile[data, { .5, .05, .95}]],
			"Median0", Transpose[Quantile[data, {.5}]],
			_, 
			mn = Mean[data];
			std = StandardDeviation[data];
			{stdM, stdP} = {mn - std, mn + std};
			{mn,stdM,stdP}
		]
	];
	
	steps = Length[dat[[1]]];
	xdata = (Range[0, steps - 1])/(steps - 1);
	pdat=(Thread[{xdata, #}]&/@ dat);
	
	fil=If[OptionValue[Method]=="Median0",None,{2 -> {3}}];
	{min, max} = If[OptionValue[PlotRange] === Automatic, MinMax[dat], OptionValue[PlotRange]];
	
	style=If[ArrayDepth[data]==1,
		OptionValue[PlotStyle],
		col = OptionValue[PlotStyle];
		(Directive[#, col,Thick] & /@ {Dashing[None], Dashed, Dashed})
	];
	
	col = OptionValue[PlotStyle];
	ListLinePlot[pdat, PlotStyle ->style, Filling -> fil, FillingStyle -> Directive[Opacity[0.2], col],
		Axes -> False, Frame -> {{True, False}, {True, False}}, FrameStyle -> Directive[{Thick, Black}],
		LabelStyle -> Directive[{Bold, Black, 14, FontFamily -> "Helvetica"}], PlotRange -> {{0, 1}, {min, max}},
		ImageSize->OptionValue[ImageSize],PlotLabel->OptionValue[PlotLabel], 
		GridLines -> {{{.5, Directive[Thick, Black, Dashed]}}, Join[System`FindDivisions[{min, max},Round[(max-min)/ OptionValue[GridLineSpacing]]], {{0, Directive[{Thick, Black}]}}]},
		FrameTicks -> {{{0, "Endo"}, {.5, "Mid"}, {1, "Epi"}}, Automatic}
	]
]


(* ::Subsection::Closed:: *)
(*MaskHelix*)


Options[MaskHelix]={BackgroundValue-> -100, SmoothHelix->False};

SyntaxInformation[MaskHelix] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

MaskHelix[helix_, mask_,OptionsPattern[]] := Block[{fun},
	fun=Switch[
		OptionValue[SmoothHelix],
		1, N[mask MedianFilter[#, 1]] /. 0. -> OptionValue[BackgroundValue],
		_, N[mask #] /. 0. -> OptionValue[BackgroundValue]
		]&;
	If[ArrayDepth[helix]==4,fun/@helix,fun[helix]]	
]


(* ::Subsection::Closed:: *)
(*BullseyePlot*)


Options[BullseyePlot]={
	TextOffset->.5,
	TextSize->12,
	PlotRange->Automatic, 
	ColorFunction->"TemperatureMap", 
	BullPlotMethod->"Dynamic", 
	TextNumberForm->{5,2},
	ImageSize->200
};

SyntaxInformation[BullseyePlot] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

BullseyePlot[data_?ArrayQ,segMask_?ArrayQ,opts:OptionsPattern[]]:=Block[{fdata}, 
	fdata=(Flatten@GetMaskData[data,#])&/@segMask;
	BullseyePlot[fdata,opts]
]

BullseyePlot[dati_?ListQ,OptionsPattern[]]:=Block[{number, radius, datat, min, max, cols, sdata, pts, disks, textv, textn, plfun, col},
	number = {6, 6, 4, 1};
	radius = {4.9, 3.6, 2.3, 1};
	datat = If[# === {} || NumberQ[#], #, Round[Mean[#], .01]] & /@ dati;
	{min, max} = If[OptionValue[PlotRange] === Automatic, MinMax[DeleteCases[datat, 0.]], OptionValue[PlotRange]];
	
	cols = (# -> Show[ColorData[#, "Image"], ImageSize -> 100]) & /@ {"GrayTones", "Rainbow", "ThermometerColors", "SunsetColors", "TemperatureMap", "GrayYellowTones", "BlueGreenYellow", "AvocadoColors", "SouthwestColors"};
	
	max = If[max <= min, min + 0.1, max];
	sdata = (datat - min)/(max - min);
	
	pts = Flatten[Table[If[i == 4, {{0, 0}}, RotationMatrix[# Degree] . {0,radius[[i]] - 1.3 OptionValue[TextOffset]} & /@ Range[0, 359, 360/number[[i]]]], {i, 4}], 1];
	
	disks = Flatten[{
		Table[{col, Disk[{0, 0}, radius[[1]], Pi/3 {i, i + 1}]}, {i, 1, 6}],
		Table[{col, Disk[{0, 0}, radius[[2]], Pi/3 {i, i + 1}]}, {i, 1, 6}],
		Table[{col, Disk[{0, 0}, radius[[3]], (i - 1) Pi/2 + {Pi/4, 3 Pi/4}]}, {i, 1, 4}],
		{{col, Disk[{0, 0}, radius[[4]], {0, 2 Pi}]}}
	}, 1];
	
	textv = Table[Text[Style[If[datat[[i]] === {} || sdata[[i]] < 0, "", NumberForm[datat[[i]],OptionValue[TextNumberForm]]], Bold, FontFamily -> "Helvetica", Black, FontSize -> OptionValue[TextSize]], pts[[i]]],{i, 17}]; 
	textn = Table[Text[Style[If[datat[[i]] === {} || sdata[[i]] < 0, "", i], Bold, FontFamily -> "Helvetica", Black, FontSize -> OptionValue[TextSize]], pts[[i]]], {i, 17}];
	
	plfun[{colf_, cstyle_}, {pText_, textVal_}] := Block[{blcol, colfunc}, 
		Legended[
			blcol = If[colf == "GrayTones", Darker[Red], Gray];
			colfunc = ColorData[colf][If[cstyle, 1 - #, #]] &;
			Graphics[{EdgeForm[{Thick, Black}], MapThread[#1 /. If[#2 === {} || #2 < 0, col -> blcol, col -> colfunc[#2]] &, {disks, sdata}], If[pText, Switch[textVal, 1, textn, 2, textv]]},
				ImageSize -> OptionValue[ImageSize]], 
			BarLegend[{colf, {min, max}}, LabelStyle -> Directive[{Bold, Black, FontFamily -> "Helvetica", FontSize -> OptionValue[TextSize]}], LegendMarkerSize -> 0.8 OptionValue[ImageSize]]
		]
	];
	
	pan = Manipulate[
		plfun[{colf, cstyle}, {pText, textVal}]
		,
		{{pText, True, "Show labels"}, {True, False}},
		{{textVal, 2, "Label"}, {1 -> "Segment", 2 -> "Value"}},
		{{colf, OptionValue[ColorFunction], "Color function"}, cols},
		{{cstyle, False, "Reverse color"}, {True, False}},
	SaveDefinitions->True, Deployed->False, SynchronousInitialization -> False];
	
	If[OptionValue[BullPlotMethod] === "Dynamic",
		NotebookClose[plotwindow];
		plotwindow = CreateWindow[DialogNotebook[{CancelButton["Close", DialogReturn[]], pan}, WindowSize -> All, WindowTitle -> "Plot data window"]];
		,
		plfun[{OptionValue[ColorFunction], False}, {True, 2}]
	]
]


(* ::Subsection:: *)
(*ExcludeSlices*)


(* ::Subsubsection::Closed:: *)
(*ExcludeSlices*)


Options[ExcludeSlices] = {CutOffMethod -> "Auto",DistanceMeasure->5, ShowOutliers->False};

SyntaxInformation[ExcludeSlices] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

ExcludeSlices[data_, OptionsPattern[]] := Block[{measure, selmask, cutoff, std, mn, bin, type,q1,q2,q3},
	type = OptionValue[DistanceMeasure];
	cutoff = OptionValue[CutOffMethod];
	(*get similarity measure*)
	measure = CalculateMeasure[data, type];
	(*calculate cutoff and selection mask*)
	cutoff = If[NumberQ[cutoff] && cutoff < .5, Quantile[Flatten@measure, cutoff],
		{q1, q2, q3} = Quantile[Flatten[measure], {0.25, 0.5, .75}];
		q2 - 2 (q3 - q1)
	];
	selmask = Mask[measure, cutoff];
	(*report the outliers*)
	If[OptionValue[ShowOutliers],ShowOutlierDistribution[measure, selmask, cutoff]];
	(*output mask*)
	selmask
]


(* ::Subsubsection::Closed:: *)
(*ShowOutlierDistribution*)


ShowOutlierDistribution[measure_, selmask_, cutoff_] := Block[{mn, fmeas, minmax},
	
	fmeas = Flatten@measure;
	mn = Mean@fmeas;
	minmax = MinMax[fmeas];
	
	(*Plot Outlier distributions*)
	Print[GraphicsRow[{
		ListPlot[fmeas, GridLines -> {None, Flatten[{mn, cutoff}]}, GridLinesStyle -> Directive@{Thick, Red}, PlotStyle -> Directive@{PointSize[0.03], Black}, PlotRange -> {All, minmax}],
		Histogram[fmeas, {(minmax[[2]] - minmax[[1]])/20}, ChartStyle -> Black, PerformanceGoal -> "Speed", GridLines -> {{mn, cutoff}, None}, GridLinesStyle -> Directive@{Thick, Red}, PlotRange -> {minmax, Full}]
	}, PlotLabel -> {mn, StandardDeviation[fmeas], cutoff}]];
	
	(*report outliers*)
	Print@Row[Flatten[{"% slices excluded: ", Round[100 Total[Flatten[(1 - selmask)]]/Length[Flatten@selmask]], "  /  ", "% directions per slice exlcuded: ", Round[100 (Total[(1 - #)]/Length[#] & /@ selmask)]}], "  "];
]


(* ::Subsubsection::Closed:: *)
(*CalculateMeasure*)


CalculateMeasure[data_, type_] := Block[{
	target, datan, fun, measure, slice, dirs, mask,mm,q1,q2,q3,iqr
	},
	
	target = Median /@ data;
	target = Flatten /@ (target/Median[Flatten[target]]);
	datan = Map[Flatten[#/Median[Flatten[#]]] &, data, {2}];
	{slice, dirs} = Dimensions[datan][[1 ;; 2]];
	
	(*select distance measure*)
	fun = Switch[type,
		1, ManhattanDistance,
		2, SquaredEuclideanDistance,
		3, EuclideanDistance,
		4, Correlation,
		5, SpearmanRho,
		_, SpearmanRho
	];
	(*calculate measure*)
	measure = Table[
		mask=Unitize[target[[i]] datan[[i, j]]];
		fun[Pick[target[[i]], mask, 1], Pick[datan[[i, j]], mask, 1]],
   	 {i, 1, slice, 1}, {j, 1, dirs, 1}];
   	 
   	 (*normalize measure*)
   	 measure = (
   	 	mm = #;
   	 	{q1, q2, q3} = Quantile[#, {0.25, 0.5, .75}];
   	 	iqr = q3 - q1;
   	 	mm = Select[mm, ((q1 - 1 iqr) < # < (q3 + 1 iqr)) &];
   	 	{q1, q2, q3} = Quantile[mm, {0.25, 0.5, .75}];
   	 	iqr = q3 - q1;
   	 	mm = Select[mm, ((q1 - 1 iqr) < # < (q3 + 1 iqr)) &];
   	 	{q1, q2, q3} = Quantile[mm, {0.25, 0.5, .75}];
   	 	1 + (# - q2)/(10 (q3 - q1))
   	 ) & /@ measure;
   	 
   	 (*make low value bad*)
   	 If[type <= 3, 2 - measure, measure]
]


(* ::Subsection::Closed:: *)
(*MakeECVBloodMask*)


Options[MakeECVBloodMask] = {BloodMaskRange -> {1400, {0, 700}}, OutputCheckImage -> True}

SyntaxInformation[MakeECVBloodMask] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

MakeECVBloodMask[pre_, post_, OptionsPattern[]] := Block[{
	mask1, mask2, bloodMask, meas, cent, preM, postM, slice, ims
	},
	
	{preM, postM} = OptionValue[BloodMaskRange];
	
	mask1 = Mask[pre, preM];
	mask2 = Mask[post, postM];
	
	bloodMask = Image3D[ImageData[
		SelectComponents[Erosion[Image[#], 2], "Count", -2]] & /@ (mask1 mask2)];
	meas = ComponentMeasurements[bloodMask, "IntensityCentroid"];
	
	cent = First@Nearest[meas[[All, 2]], Reverse@Dimensions[mask1]/2.];
	bloodMask = Erosion[#, 1] & /@ ImageData[SelectComponents[bloodMask, #IntensityCentroid == cent &]];
	
	If[OptionValue[OutputCheckImage],
		slice = Ceiling[Length[pre]/2];
		ims = (Image[((pre[[slice]]/Max[pre[[slice]]]) + bloodMask[[slice]])/2]);
		{bloodMask, ims},
		bloodMask
	]
]


(* ::Subsection::Closed:: *)
(*MakeECVBloodMask*)


ECVCalc[mappre_, mappost_, hema_?RealQ] := Block[{z, x, y, mask},
	mask = MakeECVBloodMask[mappre, mappost, OutputCheckImage -> False];
	ECVCalc[mappre, mappost, mask, hema]
]

ECVCalc[mappre_, mappost_, bloodMask_, hema_] := Block[{deltaR1, deltaR1b, rPre, rPost},
	deltaR1 = Clip[DevideNoZero[1, mappost] - DevideNoZero[1, mappre], {0, Infinity}];
	(*deltaR1b = Median@Flatten[GetMaskData[deltaR1, bloodMask]];*)
	rPre = 1./ Median@Select[GetMaskData[mappre, bloodMask], 2000 > # > 1500 &];
	rPost = 1./ Median@Select[GetMaskData[mappost, bloodMask], 500 > # > 200 &];
	deltaR1b = rPost - rPre;
	Clip[100 (deltaR1/deltaR1b) (1 - hema), {0, 100}]
]  


(* ::Subsection::Closed:: *)
(*CreateHeart*)


SyntaxInformation[CreateHeart] = {"ArgumentsPattern" -> {_.}};

CreateHeart[] := CreateHeart[0.]

CreateHeart[setin_] := Block[{
	set, col, contin, contout, shape, seto, shapeplot, topline, shapeout, con, out, shapeoutC, shout, shin
	},
	
	set = out = If[setin === 0. || ! ListQ[setin], {{73, 2.85, 2.8, 0.018, 0}, {67, 3.42, 5.76, 0.078}, 110}, setin];
	col = Gray;
	
	NotebookClose[cardiacWindow];
	cardiacWindow = DialogInput[{
		CancelButton["Generate", DialogReturn[out = seto]],
		Manipulate[
			contin = ContourPlot[With[
				{zi = 0.06 (zp - higi), xi = 0.06 (xp - 59), yi = 0},
				((xi - shifti)^2/widthi*(1 - cupi zi) + (yi)^2/widthi*(1 - cupi zi) + zi^2/lengthi^2)
			], {xp, 0, 120}, {zp, 0, 120}, Contours -> {1}, ContourStyle -> {Thick, Red}, ContourShading -> None];
			
			contout = ContourPlot[With[
				{zo = 0.06 (zp - higo), xo = 0.06 (xp - 59), yo = 0},
				(xo^2/widtho*(1 - cupo zo) + yo^2/widtho*(1 - cupo zo) + zo^2/lengtho^2)
			], {xp, 0, 120}, {zp, 0, 120}, Contours -> {1}, ContourStyle -> {Thick, Blue}, ContourShading -> None];
			
			shape = Table[With[{
				zi = 0.06 (zp - higi), xi = 0.06 (xp - 59.5), yi = 0,
				zo = 0.06 (zp - higo), xo = 0.06 (xp - 59.5), yo = 0},
				shout = If[((xo^2/widtho*(1 - cupo zo) + yo^2/widtho*(1 - cupo zo) + zo^2/lengtho^2) > 1), 0, 1];
				shin = If[(((xi - shifti)^2/widthi*(1 - cupi zi) + yi^2/widthi*(1 - cupi zi) + zi^2/lengthi^2) > 1), 0, 1];
				shout - shin
			], {zp, 120, 1, -1}, {xp, 1, 120}];
			
			seto = {{higi, lengthi, widthi, cupi, shifti}, {higo, lengtho, widtho, cupo}, top};
			shapeplot = ArrayPlot[shape];
			topline = Graphics[{
				{Green, Thick, Line[{{0, top}, {120, top}}]},
				{White, Polygon[{{0, top}, {120, top}, {120, Length[shape] + 1}, {0, Length[shape] + 1}}]}
			}];
			
			Show[shapeplot, topline, contin, contout]
			
			, Delimiter
			, {{higi, set[[1, 1]], "inner hight"}, 60, 90}
			, {{lengthi, set[[1, 2]], "inner length"}, 2, 4}
			, {{widthi, set[[1, 3]], "inner width"}, 1, 10}
			, {{cupi, set[[1, 4]], "inner cup"}, 0, 0.25}
			, {{shifti, set[[1, 5]], "inner shift"}, -1, 1}
			, Delimiter
			, {{higo, set[[2, 1]], "outer hight"}, 60, 90}
			, {{lengtho, set[[2, 2]], "outer length"}, 2, 4}
			, {{widtho, set[[2, 3]], "outer width"}, 1, 10}
			, {{cupo, set[[2, 4]], "outer cup"}, 0, 0.25}
			, Delimiter
			, {{top, set[[3]],"top loation"}, 90, 120, 1}
			, Button["set 1", {{higi, lengthi, widthi, cupi}, {higo, lengtho, widtho, cupo}, top} = {{73, 2.85, 2.8, 0.018}, {67, 3.42, 5.76, 0.078}, 110}]
			, Button["set 2", {{higi, lengthi, widthi, cupi}, {higo, lengtho, widtho, cupo}, top} = {{69, 3.1, 2.8, 0.16}, {67.5, 3.6, 7.3, 0.12}, 105}],
			
			SynchronousUpdating -> True, Method -> "Queued"
		]
	}, WindowSize -> All, WindowTitle -> "Plot data window", WindowFloating -> True, Modal -> True];
	
	shapeoutC = Compile[{{seti, _Real, 1}, {seto, _Real, 1}}, 
		Table[With[{
			zi = 0.06 (zp - seti[[1]]),
			xi = 0.06 (xp - 59.5),
			yi = 0.06 (yp - 59.5),
			zo = 0.06 (zp - seto[[1]]),
			xo = 0.06 (xp - 59.5),
			yo = 0.06 (yp - 59.5)},
			If[((xo^2/seto[[3]]*(1 - seto[[4]] zo) + yo^2/seto[[3]]*(1 - seto[[4]] zo) + zo^2/seto[[2]]^2) > 1), 0, 1] - If[(((xi - seti[[5]])^2/seti[[3]]*(1 - seti[[4]] zi) + yi^2/seti[[3]]*(1 - seti[[4]] zi) + zi^2/seti[[2]]^2) > 1), 0, 1]
		], {zp, 1, 120, 1}, {xp, 1, 120, 1}, {yp, 1, 120, 1}]
	];
	
	shapeout = shapeoutC[out[[1]], out[[2]]];
	con = ConstantArray[0, Dimensions[shapeout]];
	shapeout[[out[[3]] ;;]] = con[[out[[3]] ;;]];
	Return[{ArrayPad[shapeout, 10], {0.7, 0.7, 0.7}, seto}];
]


(* ::Section:: *)
(*End Package*)


End[](* End Private Context *)

EndPackage[]
