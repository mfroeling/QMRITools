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
"MakeRegionMesh..."

MaskToDistanceMap::usage = 
"MaskToDistanceMap..."

MaskFromDistanceMap::usage =
"MaskFromDistanceMap..."

MaskFromDistanceMap::usage =
"MaskFromDistanceMap..."

PlotMesh::usage =
"PlotMesh..."


(* ::Subsection:: *)
(*Options*)


MeshOpacity ::usage =
"MeshOpacity is an option for PlotMes..."

MeshColor ::usage =
"MeshColor is an option for PlotMes..."

MeshEdgeColor ::usage = 
"MeshEdgeColor is an option for PlotMes..."

MeshPointColor ::usage =
"MeshPointColor is an option for PlotMes..."

MeshPointSize ::usage =
"MeshPointSize is an option for PlotMes..."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*MaskToDistanceMap*)


SyntaxInformation[MaskToDistanceMap] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

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
	Transpose[ParallelMap[MaskToDistanceMap, Transpose[mask]]]
	]
];


(* ::Subsection:: *)
(*MaskFromDistanceMap*)


SyntaxInformation[MaskFromDistanceMap] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

MaskFromDistanceMap[dist_] := MaskFromDistanceMap[dist, "both"]

MaskFromDistanceMap[distI_, met_] := Block[{dist},
	dist = N[distI] /. 0. -> -0.1;
	dist = Switch[met,
		"both", UnitStep[-Ramp[-dist]],
		"out", UnitStep[-Ramp[-dist]],
		_, 1 - UnitStep[-Ramp[dist]]
	];
	SparseArray@If[ArrayDepth[dist] === 3,
		SparseArray@ImageData@SelectComponents[Image3D[dist, "Bit"], "Count", -1],
		Transpose[ParallelMap[SparseArray@ImageData[SelectComponents[Image3D[#, "Bit"], "Count", -1]] &, Transpose[dist]]]
	]
];


(* ::Subsection:: *)
(*MakeRegionMesh*)


SyntaxInformation[MaskFromDistanceMap] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeRegionMesh[data_]:=MakeRegionMesh[data, {1, 1, 1}, 1000]

MakeRegionMesh[data_, n_?IntegerQ]:=MakeRegionMesh[data, {1, 1, 1}, n]

MakeRegionMesh[data_, vox_?VectorQ]:=MakeRegionMesh[data, vox, 1000]

MakeRegionMesh[data_, vox_?VectorQ, n_?IntegerQ] := MakeRegionMesh[data, vox, {n, 0}]

MakeRegionMesh[data_, vox_?VectorQ, {n_?IntegerQ, l_?IntegerQ}] := Block[{mesh, meshS, length, lengthMax, meshOut},
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

	MeshRegion[meshOut, SphericalRegion -> True, PlotTheme -> "Default"]
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
(*PlotMesh*)


(*efficient way of plotting mesh with or without points*)
Options[PlotMesh] = {
	MeshOpacity -> 1,
	MeshColor -> StandardRed,
	MeshEdgeColor -> None,
	MeshPointColor -> None,
	MeshPointSize -> 1
}


PlotMesh[region_MeshRegion, opts : OptionsPattern[]] := PlotMesh[MeshCoordinates@region, MeshCells[region, 2], opts]

PlotMesh[pts_, cell_, opts : OptionsPattern[]] := Block[{
		meshCol, pointCol, edgeCol, ptSize, npts, opacity
	},
	npts = Length[pts];

	{meshCol, pointCol, edgeCol, ptSize, opacity} = OptionValue[{MeshColor, MeshPointColor, 
		MeshEdgeColor, MeshPointSize, MeshOpacity}];

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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
