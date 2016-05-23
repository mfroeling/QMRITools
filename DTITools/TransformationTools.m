(* ::Package:: *)

(* ::Title:: *)
(*DTITools TransformationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`TransformationTools`"];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`TransformationTools`*"];
ClearAll @@ Names["DTITools`TransformationTools"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Fuctions*)


TransformData::usage = "TransformData[data,vox,w] transforms a 3D dataset accordint to the affine transformation vector w";


(* ::Subsection::Closed:: *)
(*Options*)



(* ::Subsection::Closed:: *)
(*Error Messages*)



(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*TransformData*)

(* ::Subsubsection:: *)
(*TransformData*)


Options[TransformData]={InterpolationOrder->1}

SyntaxInformation[TransformData]={"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};

TransformData[data_, vox_, w_,OptionsPattern[]] := 
 Block[{coor, rot, coorR, interFunc, interFuncC},
  coor = GetCoordinates[data, vox];
  rot = ParametersToTransformFull[w, "Inverse"];
  coorR = ApplyRotC[coor, rot];
  interFunc = 
   Interpolation[
    Transpose[{Flatten[coor, ArrayDepth[coor] - 2], Flatten[data]}], 
    InterpolationOrder -> OptionValue[InterpolationOrder], 
    "ExtrapolationHandler" -> {0. &, "WarningMessage" -> False}];
  interFuncC = 
   Compile[{{coor, _Real, 1}}, 
    interFunc[coor[[1]], coor[[2]], coor[[3]]], 
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];
  interFuncC[coorR]
  ]


(* ::Subsubsection:: *)
(*GetCoordinates*)


GetCoordinates[data_, vox_] := Block[{dim, off, coor},
   off = Dimensions[data]/2;
   coor = MapIndexed[#2 &, data, {ArrayDepth[data]}] - 0.5;
   CoordC[coor, off, vox]
   ];


(* ::Subsubsection:: *)
(*CoordC*)

   
CoordC = Compile[{{coor, _Real, 1}, {off, _Real, 1}, {vox, _Real, 1}},
    vox (coor - off),
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection:: *)
(*ApplyRotC*)


ApplyRotC = Compile[{{coor, _Real, 1}, {rot, _Real, 2}}, 
  	(rot.Append[coor, 1])[[1 ;; 3]],
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection:: *)
(*ParametersToTransformFull*)


ParametersToTransformFull[w_] := ParametersToTransformFull[w, "Normal"]

ParametersToTransformFull[w_, opt_] := Block[{
	tx, ty, tz, rx, ry, rz, sx, sy, sz, gx, gy, gz, 
	T, R, G, S, Rx, Ry, Rz, Gx, Gy, Gz, 
	mat, rMat, tMat}, 
	
	{rx, ry, rz, tx, ty, tz, sx, sy, sz, gx, gy, gz} = w;
	rx = -rx Degree; ry = -ry Degree; rz = -rz Degree;
	T = {{1, 0, 0, tx}, {0, 1, 0, ty}, {0, 0, 1, tz}, {0, 0, 0, 1}};
	
	Rx = {{1, 0, 0, 0}, {0, Cos[rx], Sin[rx], 0}, {0, -Sin[rx], Cos[rx], 0}, {0, 0, 0, 1}};
	Ry = {{Cos[ry], 0, -Sin[ry], 0}, {0, 1, 0, 0}, {Sin[ry], 0, Cos[ry], 0}, {0, 0, 0, 1}};
	Rz = {{Cos[rz], Sin[rz], 0, 0}, {-Sin[rz], Cos[rz], 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	R = Rx.Ry.Rz;
	
	Gx = {{1, 0, gx, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	Gy = {{1, 0, 0, 0}, {gy, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	Gz = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, gz, 1, 0}, {0, 0, 0, 1}};
	G = Gx.Gy.Gz;
	
	S = {{sx, 0, 0, 0}, {0, sy, 0, 0}, {0, 0, sz, 0}, {0, 0, 0, 1}};
	
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


(* ::Section:: *)
(*End Package*)


End[]

SetAttributes[#,{Protected, ReadProtected}]&/@ Names["DTITools`TransformationTools`*"];

EndPackage[]
