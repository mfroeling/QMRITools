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
"MakeRegionMesh[data] makes a MeshRegion from a binary mask assuming an 1x1x1mm isotropic voxel size and 1000 mesh coordinates.
MakeRegionMesh[data, n] makes a MeshRegion with exactly n mesh coordinates. 
MakeRegionMesh[data, vox] makes a MeshRegion taking in account the given voxel size vox with 1000 mesh coordinates.
MakeRegionMesh[data, vox, n] makes a MeshRegion taking in account the given voxel size vox with exactly n mesh coordinates."

SplitRegionMesh::usage =
"SplitRegionMesh[mesh] splits the mesh in the MeshCoordinates and the mesh cells."

MaskToDistanceMap::usage = 
"MaskToDistanceMap[mask] converts a binary mask into a mask distance map.
MaskToDistanceMap[mask, method] converts a binary mask into a distance mask. The method can be \"both\", \"out\" or \"in\"
which calculates the mask distance inside and outside, outside, or inside the mask."

MaskFromDistanceMap::usage =
"MaskFromDistanceMap[dist] converts a mask distance map created with MaskToDistanceMap back to a binary mask.
MaskFromDistanceMap[dist, method] converts a mask distance map created with MaskToDistanceMap back to a binary mask. The method can be \"both\", \"out\" or \"in\"
which calculates the mask distance inside and outside, outside, or inside the mask."

PlotMesh::usage =
"PlotMesh[mesh] Plots a MeshRegion as a GraphicsComplex. 
PlotMesh[points, cells] plots a MeshRegion which has been split into the coordinates and cells using SplitRegionMesh."


ScaleToVolume::usage =
"ScaleToVolume[mask, vox] rescales a series of mask all to the same volume, which will be the median volume of all the masks.
ScaleToVolume[data, vox, vol] does the same but then rescales to the target volume vol."

MakeMuscleTemplate::usage =
"MakeMuscleTemplate[masks, vox] makes a muscle template form a series of masks. The output is the mean template muscle mask as a binary volume."

TemplateToVolume::usage =
"TemplateToVolume[{masks, voxM}, {template, voxT}] registers a mask or a series of masks to the given muscle template.
The default output are the mesh points per mask volume and the cell needed to generate the mesh as {points, mesh}."


MakeShapeModel::usage = 
"MakeShapeModel[points] makes the statistical shape model for the points given. The points per muscle must have the same number of coordinates.
The points can be obtained form TemplateToVolume. The output is {meanPoints, model} which are needed in FitShapeModel."

FitShapeModel::usage = 
"FitShapeModel[{meanPoints, model}, points] fits the shape model to the points, the shape model is defined by meanPoints and model which are obtained form MakeShapeModel.
FitShapeModel[{meanPoints, model}, points, n] fits the model using the first n shape vectors from the model."

ApplyShapeModel::usage = 
"ApplyShapeModel[model, std] applies the std to the shape model to create the points. The model is made by MakeShapeModel, the std is 
obtained by FitShapeModel. The std can be a number, vector or matrix. 
ApplyShapeModel[model, std, n] does the same but std is a number and n the integer indicating which vector to be used."


MakeTemplatePlot::usage = 
"MakeTemplatePlot[meshes, points] visualizes the template alignment. The input meshes and points are made by TemplateToVolume with
the TemplateOutput set to \"Meshes\"."

EvaluateModel::usage = 
"EvaluateModel[points, mesh] return a dynamic plot showing an evaluation of the shape model. 
The points and mesh are obtained from TemplateToVolume.
EvaluateModel[{mean, mat}, points, cells] does the same but for a explicit shape model given by mean and mat. 
the shape model can be obtained from MakeShapeModel."

MeshGridPlot::usage = 
"MeshGridPlot[meshes, points] gives a grid of shape meshes defined by meshes and points. These can be obtained from
with the TemplateOutput set to \"Meshes\"."

PlotShapeVariation::usage =
"PlotShapeVariation[points, mesh] plots the template mesh which shows where the model has variation for each principal
component. The needed points and mesh are obtained from TemplateToVolume.
PlotShapeVariation[points, mesh, n] does the same but only shows the first n principal components."


(* ::Subsection::Closed:: *)
(*Options*)


MeshOutput::usage = 
"MeshOutput is an option for MakeRegionMesh. Values can be \"Mesh\", \"Points\", \"Cells\", \"PointsCells\", \"All\"."

MeshOpacity::usage =
"MeshOpacity is an option for PlotMesh. It specifies how much opacity the mesh has."

MeshColor::usage =
"MeshColor is an option for PlotMesh. It specifies the color of the mesh. The color can be a single color or a list of colors
with the same length as the mesh coordinates."

MeshEdgeColor::usage = 
"MeshEdgeColor is an option for PlotMesh. It specifies the color for the border of the mesh vertices."

MeshPointColor::usage =
"MeshPointColor is an option for PlotMesh. It specifies the color of the mesh coordinates. The points need to be one and the points are 
shown as spheres."

MeshPointSize::usage =
"MeshPointSize is an option for PlotMesh. It specifies how large the points will be, its defined in mm."


MeshPoints::usage =
"MeshPoints is an option for TemplateToVolume. It specifies how many mesh coordinates will be used in the template."

TemplatePadding ::usage =
"TemplatePadding is an option for TemplateToVolume. It specifies how much the template will be padded to prevent cut off of 
larger muscles."

TemplateDilation ::usage =
"TemplateDilation is an option for TemplateToVolume. It specifies how much the template is dilated during registration.
The value is specified in mm and it helps to properly align the points."

VolumeRescale ::usage =
"VolumeRescale is an option for TemplateToVolume. If this option is set to True the volume of each muscle will be normalized 
to the mean volume of the template volumes or the template volume."

SplineSpacing ::usage =
"SplineSpacing is an option for TemplateToVolume. It defines how big the b-splines are during the registration process. Its value
is multiplied with the voxel size which then becomes the spline spacing."

TemplateOutput::usage =
"TemplateOutput is an option for TemplateToVolume. The value can be \"Meshes\" \"PointsCells\" or \"Points\". The \"Meshes\" \"PointsCells\" outputs
are needed for template and model visualization. The \"Points\" is the only thing needed for model fitting."

MeshesPerRow::usage = 
"MeshesPerRow is and option for MeshGridPlot and PlotShapeVariation. It specifies how many volumes are shown per row."


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


Options[MakeRegionMesh] = Options[MakeRegionMeshI] = {
	MeshOutput->"Mesh"
}


SyntaxInformation[MakeRegionMesh] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeRegionMesh[data_, opts : OptionsPattern[]]:=MakeRegionMesh[data, {1, 1, 1}, 1000, opts]

MakeRegionMesh[data_, n_?IntegerQ, opts : OptionsPattern[]]:=MakeRegionMesh[data, {1, 1, 1}, n, opts]

MakeRegionMesh[data_, vox_?VectorQ, opts : OptionsPattern[]]:=MakeRegionMesh[data, vox, 1000, opts]


MakeRegionMeshI[data_, vox_?VectorQ, n_?IntegerQ, opts : OptionsPattern[]]:=Block[{
		dim, size, mesh, meshS, length, lengthMax, meshOut
	},
	(*create initial mesh and smooth out the pixilation*)
	dim = Dimensions@data;
	size = Reverse[vox dim];
	mesh = DiscretizeGraphics[ListContourPlot3D[GaussianFilter[data, 2], Contours -> {0.5}, MaxPlotPoints -> Reverse[dim],
		BoxRatios -> size, DataRange -> Thread[{0, size}], PlotRange -> Thread[{0, size}], Mesh -> False]];

	(*re-mesh to approximately 2x needed cells *)
	mesh = meshOut = Remesh[mesh, Method -> {"Adaptive", "MinEdgeLength" -> EstimateEdgeLength[mesh, 2 n]}];
	(*initialize exact mesh count loop*)
	length = EstimateEdgeLength[mesh, 0.8 n];
	lengthMax = 2 length;
	(*increase min edge length till exact count is reached*)
	While[length < lengthMax, length++;
		meshOut = SimplifyMesh[mesh, {{"TriangleQuality", 4}, {"MinEdgeLength", length}, 
		{"MaxVertexCount", n}, {"MinTriangleArea", 0.5 (length^2)}}];
		If[Length[MeshCoordinates[meshOut]] === n, Break[]]
	];

	(*make the correct output*)
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


EstimateEdgeLength[mesh_, n_] := Block[{length, edges, pts},
	edges = MeshCells[mesh, 1][[All, 1]];
	pts = MeshCoordinates[mesh];
	length = If[edges === {}, {}, Norm /@ Subtract @@@ (pts[[#1]] &) /@ edges];
	Mean[length] N[Sqrt[Length[pts]/n]]
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
		meshCol, pointCol, edgeCol, ptSize, nPoints, opacity
	},
	nPoints = Length[pts];

	{meshCol, pointCol, edgeCol, ptSize, opacity} = OptionValue[{
			MeshColor, MeshPointColor, MeshEdgeColor, MeshPointSize, MeshOpacity
		}];

	(*figure out the mesh color*)
	meshCol = Which[
		meshCol === None, None,
		ColorQ[meshCol], ConstantArray[meshCol, nPoints],
		VectorQ[meshCol], meshCol,
		True, ConstantArray[meshCol, nPoints]
	];

	(*figure out the point color*)
	pointCol = Which[
		pointCol === None, False,
		pointCol === RandomColor, SeedRandom[1234]; RandomColor[nPoints],
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
						GraphicsComplex[pts, cell, VertexNormals -> MakeVertexNormals[pts, cell], 
							VertexColors -> meshCol]}
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
	MeshPoints -> {500, 2000},
	TemplatePadding -> 20,
	TemplateDilation -> 2,
	VolumeRescale -> True,
	SplineSpacing -> 5,
	TemplateOutput->"Points",
	Monitor->False
}


SyntaxInformation[TemplateToVolume] = {"ArgumentsPattern" -> {{_, _}, {_, _}, OptionsPattern[]}};

TemplateToVolume[{masks_, voxM_}, {tempI_, voxT_}, OptionsPattern[]] :=Block[{
		nLow, nHigh, padData, padReg, volRescale, space, output, volume, mon, 
		dim, template, vol, templateMesh, templatePoints, templateCells,
		moving, templateDist, movingDist, movingRigid, templateReg, deformation,
		movingMesh, pointsTemp, pointsLocal, meshes, points, 
		meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid
	},

	(*get options*)
	{{nLow, nHigh}, padData, padReg, volRescale, space, output, mon} = OptionValue[{MeshPoints, TemplatePadding, 
		TemplateDilation, VolumeRescale, SplineSpacing, TemplateOutput, Monitor}];
	mon = If[mon, MonitorFunction, List];

	(*single or multiple volumes*)
	volume = ArrayDepth[masks] === 3;
	If[! 3 <= volume <= 4, Return[$Failed]];
	
	(*Prepare template*)
	mon[DateString[], "Making template: "];
	dim = Dimensions[tempI] + 2 padData;
	template = PadToDimensions[tempI, dim];
	{templateMesh, templatePoints, templateCells} = MakeRegionMesh[template, voxT, nLow, MeshOutput -> "All"];

	(*choose to rescale or volume scale to target resolution and pad to template dimensions*)
	mon[DateString[], "Rescaling moving data: "];
	moving = If[volRescale, 
		ScaleToVolume[masks, {voxM, voxT}],
		RescaleData[masks, {voxM, voxT}, InterpolationOrder -> 0]
	];
	moving = If[volume,
		PadToDimensions[moving, dim],
		Transpose[PadToDimensions[Transpose[moving], dim]]
	];

	(*make distance map*)
	mon[DateString[], "Making distance maps: "];
	templateDist = MaskToDistanceMap@template;
	movingDist = MaskToDistanceMap@moving;

	(*rigid align moving Volume to template space*)
	mon[DateString[], "Rigid alignment: "];
	movingRigid = RegisterData[{templateDist, voxT}, {movingDist, voxT}, 
		Resolutions -> 1, MethodReg -> {"rigid"}, Iterations -> 250, NumberSamples -> 10000, 
		HistogramBins -> 128, InterpolationOrderReg -> 1, DeleteTempDirectory -> True, 
		PrintTempDirectory -> False];
	movingRigid = MaskFromDistanceMap@movingRigid;
	
	(*warp template to local aligned muscle and get displacement*)
	mon[DateString[], "Template warping: "];
	movingDist = MaskToDistanceMap[movingRigid] + padReg;
	{templateReg, deformation} = Last@RegisterDataTransform[
		{movingDist, voxT}, {templateDist, voxT}, {template, voxT}, 
		Resolutions -> 1, MethodReg -> {"affine", "bsplineMask"}, Iterations -> 250, NumberSamples -> 10000, 
		HistogramBins -> 128, InterpolationOrderReg -> 1, PrintTempDirectory -> False, 
		DeleteTempDirectory -> True, BsplineSpacing -> space voxT, ImportDeformation -> True];
	templateReg = Round@templateReg;

	(*define the mesh with template points in local space*)
	mon[DateString[], "Moving mesh points: "];
	{meshMovingRigid, pointsTemp, pointsLocal} = If[volume,
		MovePoints[{movingRigid, voxT, nHigh}, Transpose[deformation, {1, 4, 2, 3}], templatePoints],
		DistributeDefinitions[MovePoints, MakeRegionMeshI, SplitRegionMesh, EstimateEdgeLength, MeshEdgeLength,
			voxT, nHigh, templatePoints];
		Transpose[ParallelMap[(
			MovePoints[{#[[1]], voxT, nHigh}, #[[2]], templatePoints]
			)&, Transpose[{Transpose[movingRigid], Transpose[deformation, {2, 1, 5, 3, 4}]}]
		]]
	];

	Switch[output,
		"Meshes",
		mon[DateString[], "Making meshes: "];
		meshTemplate = MakeRegionMesh[template, voxT, nHigh];
		{meshTemplateReg, meshMoving} = If[volume,
			{
				MakeRegionMesh[templateReg, voxT, nHigh],
				MakeRegionMesh[moving, voxT, nHigh]
			},
			DistributeDefinitions[MakeRegionMeshI, SplitRegionMesh, EstimateEdgeLength, MeshEdgeLength, voxT, nHigh];
			{
				Transpose[ParallelMap[MakeRegionMeshI[#, voxT, nHigh]&,Transpose[templateReg]]],
				Transpose[ParallelMap[MakeRegionMeshI[#, voxT, nHigh]&,Transpose[moving]]]
			}
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
	meshH = MakeRegionMeshI[mask, vox, nHigh];
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

	(*generate standard random colors*)
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


SyntaxInformation[MakeShapeModel] = {"ArgumentsPattern" -> {_, _.}};

MakeShapeModel[points_] := MakeShapeModel[points, False]

MakeShapeModel[points_, mon_] := Block[{
		meanPoints, variation, val, vec, nVec, vecMat, std, model, plot
	},

	(*calculate the mean muscle shape*)
	meanPoints = Mean[points];
	(*subtract the mean coordinates to get the variation per point and vectorize*)
	variation = Transpose@Flatten[# - meanPoints & /@ points, {2, 3}];
	(*get shape covariance Eigensystem*)
	{val, vec} = Eigensystem[Covariance[variation]];
	(*find the maximum number of relevant vectors*)
	nVec = First@Last@Position[UnitStep[(Accumulate[val]/Total[val]) - 0.99], 0];

	(*make the model*)
	vecMat = vec[[;; nVec]];
	std = StandardDeviation /@ LeastSquares[Transpose[vecMat], Transpose[variation]];
	model = Partition[#, 3] & /@ (std vecMat);

	plot = ListLinePlot[
		Transpose[{Range[0, nVec], (Prepend[Accumulate[val[[;; nVec]]], 0]/Total[val])}], Mesh -> Full, 
		PlotRange -> {0, 1.1}, PlotLabel -> "Number of PCs for 99%", ImageSize->300, 
		GridLines -> {None, Range[0, 1.1, 0.1]}, $plotOptions];

	{meanPoints, model, If[mon, plot, Nothing]}
]


(* ::Subsubsection::Closed:: *)
(*FitShapeModel*)


SyntaxInformation[FitShapeModel] = {"ArgumentsPattern" -> {_, _, _.}};

FitShapeModel[{mean_, mat_}, points_] := FitShapeModel[{mean, mat}, points, All]

FitShapeModel[{mean_, mat_}, points_, ni_] := Block[{n, var, fit, fitPoints},
	n = If[ni===All, Length@mat, Min[{ni, Length@mat}]];

	If[MatrixQ[points],
		fit = LeastSquares[Flatten[mat[[;; n]], {2,3}], Flatten[points - mean]];
		fitPoints = mean + fit . mat[[;; n]]
		,
		fit = Transpose[LeastSquares[Flatten[mat[[;; n]], {2,3}], Transpose[Flatten[# - mean] & /@ points]]];
		fitPoints = (mean + #) & /@ (fit . mat[[;; n]])
	];
	{fit, fitPoints}
]


(* ::Subsubsection::Closed:: *)
(*ApplyShapeModel*)


SyntaxInformation[ApplyShapeModel] = {"ArgumentsPattern" -> {_, _, _.}};

ApplyShapeModel[{mean_, mat_}, std_] := ApplyShapeModel[{mean, mat}, std, 0]

ApplyShapeModel[{mean_, mat_}, std_, n_] := Which[
	NumberQ[std], mean + std If[IntegerQ[n] && n>0, mat[[n]], First[mat]],
	VectorQ[std], mean + std . mat[[;;Length[std]]],
	MatrixQ[std], (mean + # . mat[[;;Length[#]]])& /@ std
]


(* ::Subsection:: *)
(*EvaluateModel*)


(* ::Subsubsection::Closed:: *)
(*EvaluateModel*)


EvaluateModel[points_, cells_] := EvaluateModel[MakeShapeModel[points, True], points, cells]

EvaluateModel[{mean_, mat_}, points_, cells_] := EvaluateModel[{mean, mat, None}, points, cells]

EvaluateModel[{mean_, mat_, pc_}, points_, cells_] := Manipulate[
	range = MinMax[#] + {-40, 40} & /@ Transpose[Flatten[points, 1]];
	{fit, pointsFit} = FitShapeModel[{mean, mat}, points, nVec];
	plot = FindClusters[DimensionReduce[fit, 2], 3, Method -> "KMeans"];
	row = If[srow, 
		Grid[{{
			If[pc =!= None, pc, Nothing],
			ListPlot[plot, ImageSize -> 210, PlotRange -> {{-5, 5}, {-5, 5}}, Axes -> True, ImageSize -> 300,
				AspectRatio -> 1, PlotStyle -> PointSize[Large], PlotLabel -> "2D projection of " <> ToString[nVec] <> " PCs", 
				$plotOptions],
			SmoothHistogram[Transpose[fit], PlotRange -> {{-5, 5}, {0, .7}}, ImageSize -> 300, AspectRatio -> 0.7, 
				PlotLegends -> ("PC: " <> ToString[#] & /@ Range[nVec]), PlotLabel -> "Distribution of PCs", $plotOptions]
		}}, Spacings -> {3, 3}]
		, Nothing
	];

	(*make the plot*)
	MakeEvalPlot[nVec, row, {{mean, mat}, std, cells}, range, part, {col, temp}]
	,
	(*configure the plot controls*)
	{{nVec, 3, "PC"}, 1, 12, 1, ControlType -> SetterBar},
	{{part, 3, "plots per row"}, 2, 6, 1, ControlType -> SetterBar},
	{{std, 0, "range"}, -5, 5},
	{{col, "No", "distance color"},{"No","Constant","Scaled"}},
	{{temp, False, "show template"}, {True, False}},
	{{srow, False, "show pcs"}, {True, False}},
	{{dark, "Light", "export mode"}, {"Light", "Dark"}},

	(*export button*)
	Button["Make animation",
		fl = FileSelect["FileSave", {"*.gif"}];
		If[fl =!= $Canceled,
			anim = Table[MakeEvalPlot[nVec, row,  {{mean, mat}, std, cells}, range, part, {col, temp}], {std, -5, 5, 1}];
			anim = Map[Rasterize[#, ImageSize -> 1000, ImageResolution -> 150, LightDark -> dark] &, Join[anim, Reverse@anim[[2 ;; -2]]]];
			Export[ConvertExtension[fl, "gif"], anim, "DisplayDuration" -> 0.1, AnimationRepetitions -> Infinity]
		], Method -> "Queued"],

	(*hidden controls*)
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


MakeEvalPlot[nVec_, row_, {{mean_,mat_}, std_, cells_}, range_, part_, {col_, temp_}] := Block[{
		pts, pcol, step
	},
	step = Rescale[Map[Norm, mat[[;;nVec]], {-2}]];
	Grid[{
		{
			row
		}, {
			Grid[Partition[Table[
				pcol = Switch[col, 
					"No",StandardRed, 
					"Constant",ColorData["Lipari"]/@step[[j]],
					_,ColorData["Lipari"]/@Rescale[step[[j]]]
				];
				Link3DGraphic[Show[
					PlotMesh[ApplyShapeModel[{mean,mat}, std, j], cells, MeshColor->pcol, MeshPointColor -> Darker@pcol], 
					If[std===0||temp,Graphics3D[], PlotMesh[mean, cells, MeshOpacity -> 0.2, MeshColor -> Gray]], 
				PlotLabel -> Style["PC: " <> ToString[j], LightDarkSwitched[Black, White], Bold, 20],
					ImageSize -> Round[1000/part], PlotRange->range]]
			, {j, 1, nVec, 1}]	, part, part, 1, {}]]
		}
	}, Spacings -> {3, 3}, Alignment -> Center]
]


(* ::Subsubsection::Closed:: *)
(*MakeEvalPlot*)


Options[MeshGridPlot] = {
	MeshesPerRow -> 5
}

MeshGridPlot[meshesP_, pointsP_, opts : OptionsPattern[]] := MeshGridPlot[meshesP, pointsP, 15, opts]

MeshGridPlot[meshesP_, pointsP_, selI_, OptionsPattern[]] := Block[{col, range, sel, part},
	SeedRandom[1234];
	col = RandomColor[Length@pointsP[[2, 1]]];
	range = MinMax[#] + {-40, 40} & /@ Transpose[Flatten[pointsP[[1 ;; 2]], 2]];

	sel = selI;
	sel = Which[
		ListQ[sel], sel,
		IntegerQ[sel], RandomSample[Range@Length@pointsP[[2]], sel],
		True, RandomSample[Range@Length@pointsP[[2]], 15]
	];

	part = OptionValue[MeshesPerRow];

	Grid[Partition[Table[Link3DGraphic@Show[
		PlotMesh[meshesP[[4, i]], MeshOpacity -> 0.5, MeshColor -> Gray],
		ListSpherePlot[pointsP[[2, i]], SphereColor -> col], ImageSize -> Round[(1000/part)], PlotRange -> range
	], {i, sel}], part, part, 1, {}], Spacings -> {0, 0}]
]


(* ::Subsubsection::Closed:: *)
(*PlotShapeVariation*)


Options[PlotShapeVariation] = {
	MeshesPerRow -> 5
};

PlotShapeVariation[points_, cells_, opts : OptionsPattern[]] := PlotShapeVariation[points, cells, 10, opts]

PlotShapeVariation[points_, cells_, nVec_, OptionsPattern[]] := Block[{
		mean, mat, plot, part, std1, std2, pt1, pt2, diff, col
	},

	{mean, mat, plot} = MakeShapeModel[points, True];
	part = OptionValue[MeshesPerRow];
	
	Grid[Partition[Table[
		{std1, std2} = {-1, 1} 2;
		{pt1, pt2} = {mean + std1 mat[[j]], mean + std2 mat[[j]]};
		diff = Norm /@ (pt1 - pt2);
		col = ColorData["Lipari"] /@ (diff/Ceiling[Max[diff], 5]);
		Link3DGraphic@Show[PlotMesh[mean, cells, MeshColor -> col], ImageSize -> 200, 
			PlotLabel -> Style["PC: " <> ToString[j], LightDarkSwitched[Black, White], Bold, 20]]
	, {j, 1, nVec, 1}], part, part, 1, {}], Spacings -> {0, 0}, Alignment -> Center]
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
