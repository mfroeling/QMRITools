(* ::Package:: *)

(* ::Title:: *)
(*QMRITools ShapeTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`ShapeTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`ShapeTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


MakeRegionMesh::usage = 
"MakeRegionMesh[data]
MakeRegionMesh[data, n]
MakeRegionMesh[data, vox]
MakeRegionMesh[data, vox, n]
MakeRegionMesh[data, vox, {n, l}]..."

SplitRegionMesh::usage =
"SplitRegionMesh..."

MaskToDistanceMap::usage = 
"MaskToDistanceMap[mask]
MaskToDistanceMap[mask, method]..."

MaskFromDistanceMap::usage =
"MaskFromDistanceMap[dist]
MaskFromDistanceMap[dist, method]..."

PlotMesh::usage =
"PlotMesh[mesh]
PlotMesh[points, cells]..."

ScaleToVolume::usage =
"ScaleToVolume[data, vox]
ScaleToVolume[data, vox, vol]..."

MakeMuscleTemplate::usage =
"MakeMuscleTemplate[masks, vox]..."

TemplateToVolume::usage =
"TemplateToVolume[{masks, voxM}, {tempI, voxT}]..."

MakeTemplatePlot::usage = 
"MakeTemplatePlot[meshes, points]..."

MakeShapeModel::usage = 
"MakeShapeModel[points]..."

FitShapeModel::usage = 
"FitShapeModel[{mean_, mat_}, points_, n_] fit the eigensystem to points using n eigenvectors and make mesh points..."

EvaluateModel::usage = 
"EvaluateModel[points_, cells_]..."

MeshGridPlot::usage = 
"MeshGridPlot[]..."

(* ::Subsection::Closed:: *)
(*Options*)


MeshOutput::usage = 
"MeshOutput is an option for MakeRegionMesh. Values can be \"Mesh\", \"Points\", \"Cells\", \"PointsCells\", \"All\"."

MeshOpacity::usage =
"MeshOpacity is an option for PlotMesh..."

MeshColor::usage =
"MeshColor is an option for PlotMesh..."

MeshEdgeColor::usage = 
"MeshEdgeColor is an option for PlotMesh..."

MeshPointColor::usage =
"MeshPointColor is an option for PlotMesh..."

MeshPointSize::usage =
"MeshPointSize is an option for PlotMesh..."


MeshPoints::usage =
"MeshPoints is an option for TemplateToVolume..."

TemplatePadding ::usage =
"TemplatePadding is an option for TemplateToVolume..."

TemplateDilation ::usage =
"TemplateDilation is an option for TemplateToVolume..."

VolumeRescale ::usage =
"VolumeRescale is an option for TemplateToVolume..."

SplineSpacing ::usage =
"SplineSpacing is an option for TemplateToVolume..."

TemplateOutput::usage =
"TemplateOutput is an option for TemplateToVolume..."

MeshesPerRow::usage = 
"MeshesPerRow is and option for MeshGridPlot..."

(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*MaskToDistanceMap*)


SyntaxInformation[MaskToDistanceMap] = {"ArgumentsPattern" -> {_, _.}};

MaskToDistanceMap[mask_] := MaskToDistanceMap[mask, "both"]

MaskToDistanceMap[mask_, met_] := Block[{im},
	If[ArrayDepth[mask] === 3,
		MaskToDistanceMapI[mask, met],
		DistributeDefinitions[MaskToDistanceMapI, met];
		Transpose[ParallelMap[MaskToDistanceMapI[#, met]&, Normal@Transpose[mask]]]
	]
];


MaskToDistanceMapI[mask_, met_]:= Block[{im},
	im = Image3D[mask, "Bit"];
	Switch[met,
		"both", 
		ImageData[DistanceTransform[im]] - 
		ImageData[DistanceTransform[1 - im]],
		"out", -ImageData[DistanceTransform[1 - im]],
		_, ImageData[DistanceTransform[im]]
	]
]

(* ::Subsection::Closed:: *)
(*MaskFromDistanceMap*)


SyntaxInformation[MaskFromDistanceMap] = {"ArgumentsPattern" -> {_, _.}};

MaskFromDistanceMap[dist_] := MaskFromDistanceMap[dist, "both"]

MaskFromDistanceMap[distI_, met_] := Block[{dist},
	dist = N[distI] /. 0. -> -0.1;
	dist = Switch[met,
		"both", UnitStep[-Ramp[-dist]],
		"out", UnitStep[-Ramp[-dist]],
		_, 1 - UnitStep[-Ramp[dist]]
	];
	If[ArrayDepth[dist] === 3,
		MaskFromDistanceMapI[dist],
		DistributeDefinitions[MaskFromDistanceMapI];
		Transpose[ParallelMap[MaskFromDistanceMapI, Transpose[dist]]]
	]
]


MaskFromDistanceMapI[dist_] := SparseArray@ImageData@SelectComponents[Image3D[dist, "Bit"], "Count", -1]


(* ::Subsection::Closed:: *)
(*MakeRegionMesh*)


Options[MakeRegionMesh] ={
	MeshOutput->"Mesh"
}


SyntaxInformation[MakeRegionMesh] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeRegionMesh[data_, opts : OptionsPattern[]]:=MakeRegionMesh[data, {1, 1, 1}, 1000, opts]

MakeRegionMesh[data_, n_?IntegerQ, opts : OptionsPattern[]]:=MakeRegionMesh[data, {1, 1, 1}, n, opts]

MakeRegionMesh[data_, vox_?VectorQ, opts : OptionsPattern[]]:=MakeRegionMesh[data, vox, 1000, opts]

MakeRegionMesh[data_, vox_?VectorQ, n_?IntegerQ, opts : OptionsPattern[]] := MakeRegionMesh[data, vox, {n, 0}, opts]

MakeRegionMesh[data_, vox_?VectorQ, {n_?IntegerQ, l_?IntegerQ}, opts : OptionsPattern[]] := Block[{mesh, meshS, length, lengthMax, meshOut},
	(*create initial mesh and smooth out the pixilation*)
	mesh = DiscretizeGraphics[PlotContour[data, vox, ContourSmoothRadius -> 2, ContourResolution -> vox]];
	(*remesch to aproximately 2x needed cells *)
	mesh = Remesh[mesh, Method -> {"Adaptive", "MinEdgeLength" -> EstimateEdgeLength[mesh, 2 n]}];

	(*initialize exact mesh count loop*)
	length = EstimateEdgeLength[mesh, 0.8 n];
	lengthMax = 2 length;
	(*increase min edgelength till exact count is reached*)
	While[length < lengthMax, length++;
		meshOut = SimplifyMesh[mesh, {{"TriangleQuality", 4}, {"MinEdgeLength", length}, {"MaxVertexCount", n}, {"MinTriangleArea", 0.5 (length^2)}}];
		If[Length[MeshCoordinates[meshOut]] === n, Break[]]
	];

	mesh = MeshRegion[meshOut, SphericalRegion -> True, PlotTheme -> "Default"];

	Switch[OptionValue[MeshOutput], 
		"Mesh", mesh,
		"Points", First@SplitRegionMesh[mesh],
		"Cells", Last@SplitRegionMesh[mesh],
		"PointsCells", SplitRegionMesh[mesh],
		"All", Join[{mesh}, SplitRegionMesh[mesh]],
		_, mesh
	]
]


EstimateEdgeLength[mesh_, n_] := Block[{length, coors},
	length = Mean[MeshEdgeLength[mesh]];
	coors = Length[MeshCoordinates[mesh]];
	length N[Sqrt[coors/n]]
]


MeshEdgeLength[mesh_MeshRegion] := Block[{edges, pts},
	edges = MeshCells[mesh, 1][[All, 1]];
	pts = MeshCoordinates[mesh];
	If[edges === {}, {}, Norm /@ Subtract @@@ (pts[[#1]] &) /@ edges]
]


(* ::Subsection::Closed:: *)
(*SplitRegionMesh*)


SyntaxInformation[SplitRegionMesh] = {"ArgumentsPattern" -> {_}};

SplitRegionMesh[mesh_]:={MeshCoordinates@mesh, MeshCells[mesh, 2]}


(* ::Subsection::Closed:: *)
(*PlotMesh*)


(*efficient way of plotting mesh with or without points*)
Options[PlotMesh] = {
	MeshOpacity -> 1,
	MeshColor -> StandardRed,
	MeshEdgeColor -> None,
	MeshPointColor -> None,
	MeshPointSize -> 1
}


SyntaxInformation[PlotMesh] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

PlotMesh[region_MeshRegion, opts : OptionsPattern[]] := PlotMesh[SplitRegionMesh[region], opts]

PlotMesh[pts_?MatrixQ, cell_, opts : OptionsPattern[]]:=PlotMesh[{pts, cell}, opts]

PlotMesh[{pts_?MatrixQ, cell_}, opts : OptionsPattern[]] := Block[{
		meshCol, pointCol, edgeCol, ptSize, npts, opacity
	},
	npts = Length[pts];

	{meshCol, pointCol, edgeCol, ptSize, opacity} = OptionValue[{
			MeshColor, MeshPointColor, MeshEdgeColor, MeshPointSize, MeshOpacity
		}];

	(*figure out the mesh color*)
	meshCol = Which[
		meshCol === None, None,
		ColorQ[meshCol], ConstantArray[meshCol, npts],
		VectorQ[meshCol], meshCol,
		True, ConstantArray[meshCol, npts]
	];

	(*figure out the point color*)
	pointCol = Which[
		pointCol === None, False,
		pointCol === RandomColor, SeedRandom[1234]; RandomColor[npts],
		ColorQ[pointCol], Darker@pointCol,
		VectorQ[pointCol], pointCol,
		True, Darker@StandardRed
	];

	(*figure out edge color*)
	edgeCol = Which[
		ColorQ[edgeCol], edgeCol,
		True, None
	];

	(*make the plot*)
	Show[
		If[meshCol =!= None || edgeCol =!= None,
			Graphics3D[{Opacity[opacity], EdgeForm[edgeCol],
				If[meshCol === None,
					{edgeCol, GraphicsComplex[pts, cell /. Polygon -> Line]},
					{Opacity[opacity], EdgeForm[edgeCol], 
						GraphicsComplex[pts, cell, VertexNormals -> MakeVertexNormals[pts, cell], VertexColors -> meshCol]}
				]
			}], 
			Graphics3D[]
		],
		If[pointCol =!= False, 
			ListSpherePlot[pts, SphereSize -> ptSize, SphereColor -> pointCol], 
			Graphics3D[]
		], 
		Lighting -> "Neutral", SphericalRegion -> True, ImageSize -> 400, Boxed -> False
	]
]


MakeVertexNormals[p_List, cells_] := Block[{t0, t1, t2, P, ux, uy, uz, vx, vy, vz},
	{t0, t1, t2} = Transpose@cells[[All, 1]];
	P = SparseArray[Join @@ (Transpose@{#, Range@Length@#} & /@ {t0, t1, t2}) -> ConstantArray[1, 3 Length@t0]];
	{ux, uy, uz} = Transpose[p[[t1]] - p[[t0]]];
	{vx, vy, vz} = Transpose[p[[t2]] - p[[t0]]];
	Normalize /@ (P . Transpose@{uy vz - uz vy, uz vx - ux vz, ux vy - uy vx})
]


(* ::Subsection::Closed:: *)
(*ScaleToVolume*)


SyntaxInformation[ScaleToVolume] = {"ArgumentsPattern" -> {_, _, _.}};

ScaleToVolume[dat_, vox_]:=ScaleToVolume[dat, {vox, vox}, 0.]

ScaleToVolume[dat_, vox_, vol_?NumberQ]:=ScaleToVolume[dat, {vox, vox}, vol]

ScaleToVolume[dat_, {vox_, voxT_}]:=ScaleToVolume[dat, {vox, voxT}, 0.]

ScaleToVolume[dat_, {vox_, voxT_}, volI_?NumberQ] := Block[{data, vol},
	Switch[ArrayDepth[dat],
		3, 
		If[volI===0., $Failed, ScaleToVolumeI[dat, {vox, voxT}, volI]],
		4,
		data = Transpose[dat];
		vol = If[volI===0., MaskVolume[data, vox], volI];
		data = Transpose[PadToDimensions[ScaleToVolumeI[#, {vox, voxT}, vol]& /@ data]];
		First[AutoCropData[data]],
		_,$Failed
	]
]


ScaleToVolumeI[dat_, {vox_, voxT_}, vol_] := Block[{scale},
	scale = N[MaskVolume[dat, vox] / vol]^(1/3);
	Round[RescaleData[GaussianFilter[dat, 1], {vox, voxT scale}, InterpolationOrder -> 1]]
]


(* ::Subsection::Closed:: *)
(*MakeMuscleTemplate*)


SyntaxInformation[MakeMuscleTemplate] = {"ArgumentsPattern" -> {_, _}};

MakeMuscleTemplate[masks_, vox_] := MakeMuscleTemplate[masks, {vox, vox}]

MakeMuscleTemplate[masksI_, {vox_, voxT_}] := Block[{masks, sel, vol, mean, reg},
	(*rescale all muscles to same volume*)
	masks = ScaleToVolume[masksI, {vox, voxT}];
	vol = MaskVolume[Transpose@masks, voxT];
	
	(*make the distance maps for each mask*)
	sel = MaskToDistanceMap[masks];
	
	(*Make the muscle template*)
	(*register all individual volumes to mean template *)
	mean = MaskToDistanceMap[ScaleToVolume[Round@Mean@Transpose@masks, voxT, vol]];
	reg = RegisterData[{mean, voxT}, {sel, voxT}, MethodReg -> {"rigid"}, 
		NumberSamples -> 10000, HistogramBins -> 128, InterpolationOrderReg -> 1, 
		Iterations -> 250, DeleteTempDirectory -> True, PrintTempDirectory -> False];

	(*refine registration with affine transform*)
	mean = MaskToDistanceMap[ScaleToVolume[MaskFromDistanceMap@Mean@Transpose@reg, voxT, vol]];
	reg = RegisterData[{mean, voxT}, {sel, voxT}, MethodReg -> {"affineMask"}, 
		NumberSamples -> 10000, HistogramBins -> 128, InterpolationOrderReg -> 1, 
		Iterations -> 250, DeleteTempDirectory -> True, PrintTempDirectory -> False];

	(*make and rescale mean template to original mean volume*)
	mean = MaskFromDistanceMap@Mean@Transpose@reg;
	First@AutoCropData[ScaleToVolume[mean, voxT, vol]]
]


(* ::Subsection::Closed:: *)
(*TemplateToVolume*)


Options[TemplateToVolume] = {
	MeshPoints -> {500, 4000},
	TemplatePadding -> 20,
	TemplateDilation -> 2,
	VolumeRescale -> True,
	SplineSpacing -> 5,
	TemplateOutput->"Points",
	Monitor->False
}

SyntaxInformation[TemplateToVolume] = {"ArgumentsPattern" -> {{_, _}, {_, _}, OptionsPattern[]}};

TemplateToVolume[{masks_, voxM_}, {tempI_, voxT_}, OptionsPattern[]] :=Block[{
		nLow, nHigh, padData, padReg, volRescale, space, output, volume,
		dim, template, vol, templateMesh, templatePoints, templateCells,
		moving, templateDist, movingDist, movingRigid, templateReg, deformation,
		movingMesh, pointsTemp, pointsLocal, meshes, points, 
		meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid
	},

	(*get options*)
	{{nLow, nHigh}, padData, padReg, volRescale, space, output, mon} = OptionValue[{MeshPoints, TemplatePadding, 
		TemplateDilation, VolumeRescale, SplineSpacing, TemplateOutput, Monitor}];

	(*single or multiple volumes*)
	volume = ArrayDepth[masks] === 3;
	If[! 3 <= volume <= 4, Return[$Failed]];

	If[mon, Echo[DateString[], "Making template: "]];
	(*Prepare template*)
	dim = Dimensions[tempI] + 2 padData;
	template = PadToDimensions[tempI, dim];
	vol = MaskVolume[template, voxT];
	{templateMesh, templatePoints, templateCells} = MakeRegionMesh[template, voxT, nLow, MeshOutput -> "All"];

	If[mon, Echo[DateString[], "Rescalling moving data: "]];
	(*choose to rescle or volume scale to target resolution and pad to template dimensions*)
	moving = If[volRescale, 
		ScaleToVolume[masks, {voxM, voxT}, MaskVolume[template, voxT]],
		RescaleData[masks, {voxM, voxT}, InterpolationOrder -> 0]
	];
	moving = If[volume,
		PadToDimensions[moving, dim],
		Transpose[PadToDimensions[Transpose[moving], dim]]
	];

	If[mon, Echo[DateString[], "Making distance maps: "]];
	(*make distance map*)
	templateDist = MaskToDistanceMap@template;
	movingDist = MaskToDistanceMap@moving;

	(*rigid alig moving Volume to template space*)
	If[mon, Echo[DateString[], "Rigid alignment: "]];
	movingRigid = RegisterData[{templateDist, voxT}, {movingDist, voxT}, 
		Resolutions -> 1, MethodReg -> {"rigid"}, Iterations -> 250, NumberSamples -> 10000, 
		HistogramBins -> 128, InterpolationOrderReg -> 1, DeleteTempDirectory -> True, 
		PrintTempDirectory -> False];
	movingRigid = MaskFromDistanceMap@movingRigid;
	
	(*warp template to local aligned muscle and get displacement*)
	If[mon, Echo[DateString[], "Template warping: "]];
	movingDist = MaskToDistanceMap[movingRigid] + padReg;
	{templateReg, deformation} = Last@RegisterDataTransform[
		{movingDist, voxT}, {templateDist, voxT}, {template, voxT}, 
		Resolutions -> 1, MethodReg -> {"affine", "bsplineMask"}, Iterations -> 250, NumberSamples -> 10000, 
		HistogramBins -> 128, InterpolationOrderReg -> 1, PrintTempDirectory -> False, 
		DeleteTempDirectory -> True, BsplineSpacing -> space voxT, ImportDeformation -> True];
	templateReg = Round@templateReg;

	(*define the mesh with template points in local space*)
	If[mon, Echo[DateString[], "Moving mesh points: "]];
	{meshMovingRigid, pointsTemp, pointsLocal} = If[volume,
		MovePoints[{movingRigid, voxT, nHigh}, Transpose[deformation, {1, 4, 2, 3}], templatePoints],
		Transpose[Map[
			MovePoints[{#[[1]], voxT, nHigh}, #[[2]], templatePoints] &, 
			Transpose[{Transpose[movingRigid], Transpose[deformation, {2, 1, 5, 3, 4}]}]
		]]
	];

	Switch[output,
		"Meshes",
		
		If[mon, Echo[DateString[], "Making meshes: "]];
		meshTemplate = MakeRegionMesh[template, voxT, nHigh];
		If[volume,
			meshTemplateReg = MakeRegionMesh[templateReg, voxT, nHigh];
			meshMoving = MakeRegionMesh[moving, voxT, nHigh];
			,
			{meshTemplateReg, meshMoving} = Monitor[Transpose@Table[{
				MakeRegionMesh[templateReg[[All, ni]], voxT, nHigh],
				MakeRegionMesh[moving[[All, ni]], voxT, nHigh]
			}, {ni, 1, Length[pointsTemp]}], ProgressIndicator[ni, {1, Length[pointsTemp]}]];
		];

		{
			{meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid},
			{pointsTemp, pointsLocal, templatePoints, templateCells}
		},

		"PointsCells", 
		{pointsLocal, templateCells},

		_, pointsLocal]
]


MovePoints[{mask_, vox_, nHigh_}, delta_, points_] := Block[{meshH, temp, local},
	meshH = MakeRegionMesh[mask, vox, nHigh];
	temp = points - (delta[[#[[3]], #[[2]], #[[1]]]] & /@ Round[#/vox & /@ points]);
	local = RegionNearest[meshH, temp];
	{meshH, temp, local}
]


(* ::Subsection::Closed:: *)
(*MakeTemplatePlot*)


SyntaxInformation[MakeTemplatePlot] = {"ArgumentsPattern" -> {_, _}};

MakeTemplatePlot[meshes_, points_] := Block[{
		meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid,
		pointsTemp, pointsLocal, templatePoints, templateCells
	},
	{meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid} = meshes;
	{pointsTemp, pointsLocal, templatePoints, templateCells} = points;
	If[Length[meshTemplateReg] > 0,
		Table[MakeTemplatePlot[{meshTemplate, meshTemplateReg[[ni]], meshMoving[[ni]], meshMovingRigid[[ni]], 
			pointsTemp[[ni]], pointsLocal[[ni]], templatePoints}], {ni, 1, Length@meshTemplateReg}],
		MakeTemplatePlot[{meshTemplate, meshTemplateReg, meshMoving, 
			meshMovingRigid, pointsTemp, pointsLocal, templatePoints}]
	]
]


MakeTemplatePlot[{meshTemplate_, meshTemplateReg_, meshMoving_, meshMovingRigid_, 
		pointsTemp_, pointsLocal_, templatePoints_}] := Block[{
		range, templateRegPlot, movingRigPlot, templatePlot, movingPlot, col
	},

	(*data range in mm after padding*)
	range = MinMax[#] + {-40, 40} & /@ Transpose[Flatten[{pointsTemp, pointsLocal, templatePoints}, 1]];

	(*make the mesh plots*)
	templateRegPlot = PlotMesh[meshTemplateReg, MeshOpacity -> 0.5, MeshColor -> StandardGreen];
	movingRigPlot = PlotMesh[meshMovingRigid, MeshOpacity -> 0.5, MeshColor -> Gray];
	templatePlot = PlotMesh[meshTemplate, MeshOpacity -> 0.5, MeshColor -> StandardGreen];
	movingPlot = PlotMesh[meshMoving, MeshOpacity -> 0.5, MeshColor -> Gray];

	(*generate stadard random colors*)
	SeedRandom[1234];
	col = RandomColor[Length@templatePoints];

	(*Make plot grid*)
	Grid[{{
		Link3DGraphic@Show[templatePlot, movingPlot, ImageSize -> 300, 
			PlotLabel -> "Before template alignment", PlotRange -> range],
		Link3DGraphic@Show[templatePlot, movingRigPlot, ImageSize -> 300, 
			PlotLabel -> "After rigid template alignment", PlotRange -> range], 
		Link3DGraphic@Show[templateRegPlot, movingRigPlot, ImageSize -> 300, 
			PlotLabel -> "After bspline template alignment", PlotRange -> range]
	}, {
		Link3DGraphic@Show[movingRigPlot, ListSpherePlot[templatePoints, SphereColor -> col], 
			ImageSize -> 300, PlotLabel -> "Template points", PlotRange -> range],
		Link3DGraphic@Show[movingRigPlot, ListSpherePlot[pointsTemp, SphereColor -> col], 
			ImageSize -> 300, PlotLabel -> "Registered template points", PlotRange -> range],
		Link3DGraphic@Show[movingRigPlot, ListSpherePlot[pointsLocal, SphereColor -> col], 
			ImageSize -> 300, PlotLabel -> "Registered template points on mesh", PlotRange -> range]
	}}]
]


(* ::Subsection:: *)
(*ShapeModel*)


(* ::Subsubsection::Closed:: *)
(*MakeShapeModel*)


MakeShapeModel[points_] := MakeShapeModel[points, False]

MakeShapeModel[points_, mon_] := Block[{
		meanPoints, variation, val, vec, nvec, vecMat, std, model, plot
	},

	(*calculate the mean muscle shape*)
	meanPoints = Mean[points];
	(*subtract the mean coordinates to get the variation per point and v	ectorize*)
	variation = Transpose@Flatten[# - meanPoints & /@ points, {2, 3}];
	(*get shape covariance eigensystem*)
	{val, vec} = Eigensystem[Covariance[variation]];
	(*find the maximum number of relevant vectors*)
	nvec = First@Last@Position[UnitStep[(Accumulate[val]/Total[val]) - 0.99], 0];

	(*make the model*)
	vecMat = vec[[;; nvec]];
	std = StandardDeviation /@ Transpose[(FitVariationVec[variation, vecMat, nvec])];
	model = std vecMat;

	plot = ListLinePlot[
		Transpose[{Range[0, nvec], (Prepend[Accumulate[val[[;; nvec]]], 0]/Total[val])}], Mesh -> Full, 
		PlotRange -> {0, 1.1}, PlotLabel -> "Number of PCs for 99%", ImageSize->300, 
		GridLines -> {None, Range[0, 1.1, 0.1]}, $plotOptions];

	{meanPoints, model, If[mon, plot, Nothing]}
]


(* ::Subsubsection::Closed:: *)
(*FitShapeModel*)


FitShapeModel[{mean_, mat_}, points_] := FitShapeModel[{mean, mat}, points, All]

FitShapeModel[{mean_, mat_}, points_, ni_] := Block[{n, var, fit, fitPoints},
	n = If[ni===All, Length@mat, Min[{ni, Length@mat}]];
	If[MatrixQ[points],
		var = points - mean;
		fit = FitVariationVec[var, mat, n];
		fitPoints = mean + Partition[fit . mat[[;; n]], 3],
		var = Transpose@Flatten[# - mean & /@ points, {2, 3}];
		fit = FitVariationVec[var, mat, n];
		fitPoints = (mean + Partition[#, 3]) & /@ (fit . mat[[;; n]])
	];
	{fit, fitPoints}
]


(* ::Subsubsection::Closed:: *)
(*FitVariationVec*)


(*fit the eigensystem to points using n eigenvectors*)
FitVariationVec[var_, mat_] := FitVariationVec[var, mat, All]

FitVariationVec[var_, mat_, n_] := var . PseudoInverse[mat[[;; n]]]


(* ::Subsection:: *)
(*EvaluateModel*)


(* ::Subsubsection::Closed:: *)
(*EvaluateModel*)


EvaluateModel[points_, cells_] := EvaluateModel[MakeShapeModel[points, True], points, cells]

EvaluateModel[{mean_, mat_}, points_, cells_] := EvaluateModel[{mean, mat, None}, points, cells]

EvaluateModel[{mean_, mat_, pc_}, points_, cells_] := Manipulate[
	range = MinMax[#] + {-40, 40} & /@ Transpose[Flatten[points, 1]];
	{fit, pointsFit} = FitShapeModel[{mean, mat}, points, nvecs];
	plot = FindClusters[DimensionReduce[fit, 2], 3, Method -> "KMeans"];
	row = Grid[{{
		If[pc =!= None, pc, Nothing],
		ListPlot[plot, ImageSize -> 210, PlotRange -> {{-5, 5}, {-5, 5}}, Axes -> True, ImageSize -> 300,
			AspectRatio -> 1, PlotStyle -> PointSize[Large], PlotLabel -> "2D projection of " <> ToString[nvecs] <> " PCs", 
			$plotOptions],
		SmoothHistogram[Transpose[fit], PlotRange -> {{-5, 5}, {0, .7}}, ImageSize -> 300, AspectRatio -> 0.7, 
			PlotLegends -> ("PC: " <> ToString[#] & /@ Range[nvecs]), PlotLabel -> "Distribution of PCs", $plotOptions]
	}}, Spacings -> {3, 3}];

	MakeEvalPlot[nvecs, row, {mean, std, mat, cells}, range]
	,
	{{nvecs, 5, "PC"}, 1, 10, 1, ControlType -> SetterBar},
	{{std, 0, "range"}, -5, 5},
	Button["Make animation",
		fl = FileSelect["FileSave", {"*.gif"}];
		If[fl =!= $Canceled,
			anim = Table[MakeEvalPlot[nvecs, row, {mean, stdi, mat, cells}, range], {stdi, -5, 5, 1}];
			anim = Map[Rasterize[#, ImageSize -> 1000, ImageResolution -> 150] &, Join[anim, Reverse@anim[[2 ;; -2]]]];
			Export[ConvertExtension[fl, "gif"], anim, "DisplayDuration" -> 0.1, AnimationRepetitions -> Infinity]
		], Method -> "Queued"],
	{fit, ControlType -> None},
	{pointsFit, ControlType -> None},
	{plot, ControlType -> None},
	{row, ControlType -> None},
	{anim, ControlType -> None},
	{fl, ControlType -> None},
	{range, ControlType -> None}
]


(* ::Subsubsection::Closed:: *)
(*MakeEvalPlot*)


MakeEvalPlot[nvevs_, row_, {mean_, std_, mat_, cells_}, range_] := Grid[{
	{
		row
	}, {
		Grid[
			Partition[Table[
				Link3DGraphic@Show[
						PlotMesh[mean + Partition[std  mat[[j]], 3], cells, MeshPointColor -> Darker@Red], 
						If[std===0,Graphics3D[],PlotMesh[mean, cells, MeshOpacity -> 0.2, MeshColor -> Gray]], 
					ImageSize -> 200, PlotLabel -> Style["PC: " <> ToString[j], Black, Bold, 20],
					PlotRange->range]
			, {j, 1, nvevs, 1}], 5, 5, 1, {}]
		]
	}
}, Spacings -> {3, 3}, Alignment -> Center]


(* ::Subsubsection::Closed:: *)
(*MakeEvalPlot*)



Options[MeshGridPlot] = {
	MeshesPerRow -> 5
}

MeshGridPlot[meshesP_, pointsP_, sel_, opts : OptionsPattern[]] := MeshGridPlot[meshesP, pointsP, 15, opts]

MeshGridPlot[meshesP_, pointsP_, sel_, OptionsPattern[]] := Block[{col, range, part},
	SeedRandom[1234];
	col = RandomColor[Length@pointsP[[2, 1]]];
	range = MinMax[#] + {-40, 40} & /@ Transpose[Flatten[pointsP[[1 ;; 2]], 2]];

	sel = Which[
		ListQ[sel], sel,
		IntegerQ[sel], RandomSample[Range@Length@pointsP[[2]], sel],
		True, RandomSample[Range@Length@pointsP[[2]], 15]
	];

	part = OptionValue[MeshesPerRow];

	Grid[Partition[Table[Link3DGraphic@Show[
		PlotMesh[meshesP[[4, i]], MeshOpacity -> 0.5, MeshColor -> Gray],
		ListSpherePlot[pointsP[[2, i]], SphereColor -> col], ImageSize -> 300, PlotRange -> range
	], {i, sel}], part, part, 1, {}], Spacings -> {0, 0}]
]




(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
