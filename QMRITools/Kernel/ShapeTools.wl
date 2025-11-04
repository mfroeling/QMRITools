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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*MaskToDistanceMap*)


SyntaxInformation[MaskToDistanceMap] = {"ArgumentsPattern" -> {_, _.}};

MaskToDistanceMap[mask_] := MaskToDistanceMap[mask, "both"]

MaskToDistanceMap[mask_, met_] := Block[{im},
	If[ArrayDepth[mask] === 3,
		im = Image3D[mask, "Bit"];
		Switch[met,
			"both", 
			ImageData[DistanceTransform[im]] - 
			ImageData[DistanceTransform[1 - im]],
			"out", -ImageData[DistanceTransform[1 - im]],
			_, ImageData[DistanceTransform[im]]
		],
	Transpose[ParallelMap[MaskToDistanceMap, Normal@Transpose[mask]]]
	]
];


(* ::Subsection:: *)
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
		SparseArray@ImageData@SelectComponents[Image3D[dist, "Bit"], "Count", -1],
		Transpose[ParallelMap[SparseArray@ImageData[SelectComponents[Image3D[#, "Bit"], "Count", -1]] &, Transpose[dist]]]
	]
];


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*SplitRegionMesh*)


SyntaxInformation[SplitRegionMesh] = {"ArgumentsPattern" -> {_}};

SplitRegionMesh[mesh_]:={MeshCoordinates@mesh, MeshCells[mesh, 2]}


(* ::Subsection:: *)
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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*ScaleToVolume*)


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


(* ::Subsection:: *)
(*TemplateToVolume*)


Options[TemplateToVolume] = {
	MeshPoints -> {500, 4000},
	TemplatePadding -> 20,
	TemplateDilation -> 2,
	VolumeRescale -> True,
	SplineSpacing -> 5
}


SyntaxInformation[TemplateToVolume] = {"ArgumentsPattern" -> {{_, _}, {_, _}, OptionsPattern[]}};

TemplateToVolume[{masks_, voxM_}, {tempI_, voxT_}, OptionsPattern[]] :=Block[{
		nLow, nHigh, padData, padReg, volRescale, space, volume,
		dim, template, vol, templateMesh, templatePoints, templateCells,
		moving, templateDist, movingDist, movingRigid, templateReg, 
		deformation,
		movingMesh, pointsTemp, pointsLocal, meshes, points,
		meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid
	},

	(*get options*)
	{{nLow, nHigh}, padData, padReg, volRescale, space} = OptionValue[{MeshPoints, TemplatePadding, 
		TemplateDilation, VolumeRescale, SplineSpacing}];

	(*single or multiple volumes*)
	volume = ArrayDepth[masks] === 3;
	If[! 3 <= volume <= 4, Return[$Failed]];

	(*Prepare template*)
	dim = Dimensions[tempI] + 2 padData;
	template = PadToDimensions[tempI, dim];
	vol = MaskVolume[template, voxT];
	{templateMesh, templatePoints, templateCells} = MakeRegionMesh[template, voxT, nLow, MeshOutput -> "All"];

	(*choose to rescle or volume scale to target resolution and pad to t\
	emplate dimensions*)
	moving = If[volRescale, 
		ScaleToVolume[masks, {voxM, voxT}, MaskVolume[template, voxT]],
		RescaleData[masks, {voxM, voxT}, InterpolationOrder -> 0]
	];
	moving = If[volume,
		PadToDimensions[moving, dim],
		Transpose[PadToDimensions[Transpose[moving], dim]]
	];

	(*make distance map*)
	templateDist = MaskToDistanceMap@template;
	movingDist = MaskToDistanceMap@moving;

	(*rigid alig moving Volume to template space*)
	movingRigid = MaskFromDistanceMap@RegisterData[{templateDist, voxT}, {movingDist, voxT}, 
		Resolutions -> 1, MethodReg -> {"rigid"}, Iterations -> 250, NumberSamples -> 10000, 
		HistogramBins -> 128, InterpolationOrderReg -> 1, DeleteTempDirectory -> True, 
		PrintTempDirectory -> False];

	(*warp template to local aligned muscle and get displacement*)
	{templateReg, deformation} = Last@RegisterDataTransform[
		{MaskToDistanceMap[movingRigid] + padReg, voxT}, {templateDist, voxT}, {template, voxT}, 
		Resolutions -> 1, MethodReg -> {"affine", "bsplineMask"}, Iterations -> 250, NumberSamples -> 10000, 
		HistogramBins -> 128, InterpolationOrderReg -> 1, PrintTempDirectory -> False, 
		DeleteTempDirectory -> True, BsplineSpacing -> space voxT, ImportDeformation -> True];
	templateReg = Round@templateReg;

	(*define the mesh with template points in local space*)
	{meshMovingRigid, pointsTemp, pointsLocal} = If[volume,
		MovePoints[{movingRigid, voxT, nHigh}, Transpose[deformation, {1, 4, 2, 3}], templatePoints],
		Transpose[Map[
			MovePoints[{#[[1]], voxT, nHigh}, #[[2]], templatePoints] &, 
			Transpose[{Transpose[movingRigid], Transpose[deformation, {2, 1, 5, 3, 4}]}]
		]]
	];

	meshTemplate = MakeRegionMesh[template, voxT, nHigh];
	If[volume,
		meshTemplateReg = MakeRegionMesh[templateReg, voxT, nHigh];
		meshMoving = MakeRegionMesh[moving, voxT, nHigh];
		,
		{meshTemplateReg, meshMoving} = Transpose@Table[{
			MakeRegionMesh[templateReg[[All, ni]], voxT, nHigh],
			MakeRegionMesh[moving[[All, ni]], voxT, nHigh]
		}, {ni, 1, Length[pointsTemp], 3}];
	];

	meshes = {meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid};
	points = {pointsTemp, pointsLocal, templatePoints, templateCells};

	{meshes, points}
]

MovePoints[{mask_, vox_, nHigh_}, delta_, points_] := Block[{meshH, temp, local},
	meshH = MakeRegionMesh[mask, vox, nHigh];
	temp = points - (delta[[#[[3]], #[[2]], #[[1]]]] & /@ Round[#/vox & /@ points]);
	local = RegionNearest[meshH, temp];
	{meshH, temp, local}
]


(* ::Subsection:: *)
(*MakeTemplatePlot*)


SyntaxInformation[MakeTemplatePlot] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

MakeTemplatePlot[meshes_, points_] := Block[{
		meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid,
		pointsTemp, pointsLocal, templatePoints, templateCells
	},
	{meshTemplate, meshTemplateReg, meshMoving, meshMovingRigid} = 
	meshes;
	{pointsTemp, pointsLocal, templatePoints, templateCells} = points;
	If[Length[meshTemplateReg] > 0,
		Table[MakeTemplatePlot[{meshTemplate, meshTemplateReg[[ni]], meshMoving[[ni]], meshMovingRigid[[ni]], 
			pointsTemp[[ni]], pointsLocal[[ni]], templatePoints}], {ni, 1, Length@meshTemplateReg}],
		MakeTemplatePlot[{meshTemplate, meshTemplateReg, meshMoving, 
			meshMovingRigid, pointsTemp, pointsLocal, templatePoints}]
	]
]

MakeTemplatePlot[{meshTemplate_, meshTemplateReg_, meshMoving_, meshMovingRigid_, pointsTemp_, pointsLocal_, templatePoints_}] := Block[{
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



(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
