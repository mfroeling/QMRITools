(* ::Package:: *)

(* ::Title:: *)
(*DTITools RegistrationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`RegistrationTools`"];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`RegistrationTools`*"];
ClearAll @@ Names["DTITools`RegistrationTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Fuctions*)


GridData::usage = 
"GridData[data] returns a array of coordinates, the coordinates are between 1 and de dimensions of the dataset. the spacing between the coordinates is 1.
GridData[data,sp] returns a array of coordinates, the coordinates are between 1 and de dimensions of the dataset. the spacing between the coordinates is sp.
GridData[data,{s1,s2}] returns a array of coordinates, the coordinates are between 1 and de dimensions of the dataset. \
For a 2D dataset the spacing between the coordinates is s1 and s2 for each dimensions respectively. \
For a 3D dataset the spacing between the coordinates is s1 for the slice direction and s2 in plane.
GridData[data,{s1,s2,s3}] returns a array of coordinates, the coordinates are between 1 and de dimensions of the dataset. \
For a 3D dataset the spacing between the coordinates is \"s1\", \"s2\" and \"s3\" for each dimensions respectively."

TransformGrid2D::usage = 
"TransformGrid2D[gridlist, vox, w, dim] transforms the 2D coordinates in gridlist acoording to w = (\[Theta], dx, dy, sx, sy, q). The transformation is arround the coordinate 0.5 dim = 0.5 (dimx, dimy). \
Returns the transformed list of 2D coordinates."

TransformGrid3D::usage = 
"TransformGrid3D[gridlist, vox, w, dim] transforms the 3D coordinates in gridlist acoording to w = (\[Theta]x, \[Theta]y, \[Theta]z, dx, dy, dz, sx, sy, sz, q1, q2, q3). The transformation is arround the coordinate 0.5 dim = 0.5 (dimx, dimy, dimz) \
Returns the transformed list of 3D coordinates."

GridInterpolation2D::usage = 
"GridInterpolation2D[image, grid] interpolates a image for the given grid of 2D coordinates (GridData and TransformGrid2D). \
The images is a 2D matrix of values, the grid is a 2D matrix with 2D coordinates. \
The dimensions of the grid can be different than the dimensions of the image (both 2D), \
but the 2D coordinates given in the grid can not exceed the 2D dimensions of the image."

GridInterpolation3D::usage = 
"GridInterpolation3D[data, grid] interpolates a 3D dataset data for the given grid of 3D coordinates (GridData and TransformGrid3D). \
The images is a 3D matrix of values, the grid is a 3D matrix with 3D coordinates. \
The dimensions of the grid can be different than the dimensions of the image (both 3D), \
but the 3D coordinates given in the grid can not exceed the 3D dimensions of the image."

DefImage::usage = 
"DefImage[image, vox, w] deforms a image according to the 6 degrees of freedom w = (\[Theta], dx, dy, sx, sy, q). \
With \[Theta] = rotation, dx and dy = translation, sx and sy = scale and q = skew. \
The transformation is arround the coordinate 0.5 dim = 0.5 (dimx, dimy)."

DefData::usage = 
"DefData[data, vox, w] deforms a 3D data set according to the 12 degrees of freedom w = (\[Theta]x, \[Theta]y, \[Theta]z, dx, dy, dz, sx, sy, sz, q1, q2, q3). \
With \[Theta]x, \[Theta]y and \[Theta]z = rotation, dx, dy and dz = translation, sx, sy and sz = scale and q1, q2 and q3 = skew. \
The transformation is arround the coordinate 0.5 dim = 0.5 (dimx, dimy, dimz), with dim the dimensions of the dataset."

RemovePeaks::usage = 
"RemovePeaks[data] removes the bright spots in the data and smooths out the image if there is a gradient.
RemovePeaks[data,perc] removes the bright spots in the data. Perc is the percentage of the data that is considered a peak, default is 0.96.
RemovePeaks[data,perc,kern] removes the bright spots in the data. Perc is the percentage of the data that is considered a peak, default is 0.96. \
kern defines how sharp the peaks are removed. Large kernel very smooth removal, small kernel verry sharp removal."

JointHistogram::usage=
"JointHistogram[data, {range,step}] creates a joint histgoram of diffusion data.";


(* ::Subsection::Closed:: *)
(*Options*)



(* ::Subsection::Closed:: *)
(*Error Messages*)


GridData::dim = "Dataset should be 2D of 3D, current data is `1`D";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*General*)


(* ::Subsubsection::Closed:: *)
(*JointHistogram*)


SyntaxInformation[JointHistogram] = {"ArgumentsPattern" -> {_, {_,_}}};

JointHistogram[data_, {r1_, s1_}] := 
 Module[{dat = Transpose[data] // N, fdat, hist},
  fdat = DeleteCases[
      DeleteCases[
       Transpose[{Flatten[dat[[1]]], Flatten[#]}], {_, 
        0.}], {0., _}] & /@ dat[[2 ;;]];
  hist = Reverse[BinCounts[#, {0, r1, s1}, {0, r1, s1}]] & /@ fdat;
  Log[hist /. {0 -> 1, 1 -> 1.1}] // N
  ]


(* ::Subsubsection::Closed:: *)
(*RemovePeaks*)


SyntaxInformation[RemovePeaks] = {"ArgumentsPattern" -> {_, _., _.}};

RemovePeaks[data_]:=RemovePeaks[data,0.96,2];

RemovePeaks[data_,perc_]:=RemovePeaks[data,perc,2];

RemovePeaks[data_,perc_,kern_]:=
Module[{quan,tops},
	quan=Round[Quantile[Flatten[GaussianFilter[data,kern],{2,3}],perc]];
	tops=GaussianFilter[Clip[data-quan,{0,Infinity}],kern];
	Clip[data-tops,{0,Infinity}]
	];


(* ::Subsubsection::Closed:: *)
(*Grid Data 2 D and 3 D*)


SyntaxInformation[GridData] = {"ArgumentsPattern" -> {_, _.}};

GridData[data_,step_:1]:=
Module[{s1,s2,s3,idim,jdim,kdim},
	If[!(ArrayDepth[data]==2||ArrayDepth[data]==3),
		Return[Message[GridData::dim,ArrayDepth[data]]];
		,
		(* define grid spacing in each of the directions, one number - same for all directions 2D and 3D, two numbers - 2D, three numbers - 3D*)
		Switch[
			Length[step],
			0,s1=s2=s3=step;,
			1,s1=s2=s3=step[[1]];,
			2,s1=step[[1]];s2=s3=step[[2]];,
			3,s1=step[[1]];s2=step[[2]];s3=step[[3]];
			];
		(* grid creations, 2D or 3D*)
		Switch[
			ArrayDepth[data],
			2,(*2D dataset*)
			{idim,jdim}=Dimensions[data];
			Table[{i,j},{i,1,idim,s1},{j,1,jdim,s2}],
			3,(*3D dataset*)
			{idim,jdim,kdim}=Dimensions[data];
			Table[{i,j,k},{i,1,idim,s1},{j,1,jdim,s2},{k,1,kdim,s3}]
			]
		]
	];


(* ::Subsection:: *)
(*2D Functions*)


(* ::Subsubsection::Closed:: *)
(*Grid Transformation 2D*)


SyntaxInformation[TransformGrid2D] = {"ArgumentsPattern" -> {_, _, _, _}};

TransformGrid2D[grid_List, vox_List, w_List, dim_List] := 
 Module[{off, part},
  off = Round[.5 dim];
  part = Dimensions[grid][[2]];
  Partition[Transpose[(Transpose[GridTransform2D[Transpose[(Transpose[Flatten[grid, 1]] - off)*vox], w]]/vox) + off], part]
  ]

GridTransform2D = (*input parameters {w, dx, dy, s1, s2, a}*)
Compile[
	{{pts,_Real,2},
	{p,_Real,1}},
	With[
		{
		(*p1, p2: scaling factors; p2, p3: skewing factors; p5, p6: translation factors
		rotation is a combination of skew and scale *)
		dx = p[[2]],dy = p[[3]],
		S={{1,0},{0,1}}+(0.1{{p[[4]],p[[6]]},{0,p[[5]]}}),
		R=N[{{Cos[#],-Sin[#]},{Sin[#],Cos[#]}}]&@@{p[[1]]Degree}
		},
		Map[N[R.S.#+{dx,dy}]&,pts]
		]
	];


(* ::Subsubsection::Closed:: *)
(*Grid Interpolation 2 D*)


(* ::Text:: *)
(*These functions make fast interpolation possible of an image on a grid.*)


Options[GridInterpolation2D] = {InterpolationOrder -> 1};

SyntaxInformation[GridInterpolation2D] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

GridInterpolation2D[dat_List, grd_List, opts: OptionsPattern[]] := 
 Module[{interFunction, data, grid, ddim, gdim, cgrid},
  grid = Flatten[grd, 1];
  gdim = Dimensions[grd];
  
  data = ArrayPad[dat, 1] // N;
  ddim = Dimensions[data];
  
  With[{
    (*Calls the interpolation function for grid coordinates*)
    InterpolateOnGrid = Compile[{{pts, _Real, 2}}, Map[interFunction[#[[1]], #[[2]]] &, pts]]
    },
   (*Define interpolation function*)
   interFunction = ListInterpolation[data,opts];
   (*actual grid interpolation*)
   cgrid = Transpose[Clip[grid[[All, #]] + 1, {1, ddim[[#]]}] & /@ {1, 2}];
   Partition[InterpolateOnGrid[cgrid], gdim[[2]]]
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*Deform image*)


SyntaxInformation[DefImage] = {"ArgumentsPattern" -> {_, _, _}};

DefImage[image_, vox_, w_] := GridInterpolation2D[image,TransformGrid2D[GridData[image, 1], vox, w, Dimensions[image]]]


(* ::Subsection:: *)
(*3D Functions*)


(* ::Subsubsection::Closed:: *)
(*Grid Transformation 3D*)


(* ::Text:: *)
(*The function MyRigid2D transforms a grid with parameters w where w=(\[Theta], \[Gamma], \[Psi],Subscript[t, x],Subscript[t, y], Subscript[t, z], Subscript[s, x], Subscript[s, y],Subscript[s, z], Subscript[s, xy],Subscript[s, xz],Subscript[s, yz]), *)


TransformGrid3D[grid_List,vox_List,w_List,dim_]:=
Module[{off,part1,part2},
	off=Round[.5dim];
	part1=Dimensions[grid][[3]];
	part2=Dimensions[grid][[2]];
	Partition[Partition[Transpose[(Transpose[GridTransform3D[Transpose[(Transpose[Flatten[grid,2]]-off)*vox],w]]/vox)+off],part1],part2]
]

GridTransform3D = (*input parameters {w1, w2, w3, dx, dy, dz, sx, sy, sz, sxy, sxz, syz}*)
Compile[
	{{pts,_Real,2},
	{p,_Real,1}},
	Block[
		{
			(*p1, p2, p3 - rotation on 3 axis; p4, p5 ,p6 - translation on 3 axis; p7, p8, p9 - scaling on 3 axis; p10, p11, p12 - skew*)
			dx = p[[5]],
			dy = p[[6]],
			dz = p[[4]],
			S={{1,0, 0},{0,1, 0}, {0, 0, 1}}+(0.1{{p[[7]],0,0},{0,p[[8]],0},{0,0,p[[9]]}}),
			G={{1,0.1 p[[10]] 0.1 p[[12]], 0.1 p[[10]]},{ 0.1 p[[11]],1,0},{0, 0.1 p[[12]],1}},
			p1=p[[3]]Degree,
			p2=p[[1]]Degree,
			p3=p[[2]]Degree,
			R
		},
		R={
			{Cos[p2] Cos[p3],-Cos[p2] Sin[p3],Sin[p2]},
			{Cos[p3] Sin[p1] Sin[p2]+Cos[p1] Sin[p3],Cos[p1] Cos[p3]-Sin[p1] Sin[p2] Sin[p3],-Cos[p2] Sin[p1]},
			{-Cos[p1] Cos[p3] Sin[p2]+Sin[p1] Sin[p3],Cos[p3] Sin[p1]+Cos[p1] Sin[p2] Sin[p3],Cos[p1] Cos[p2]}
		};
		Map[N[R.G.S.#+{dz,dx,dy}]&,pts]
		]
	];


(* ::Subsubsection::Closed:: *)
(*Grid Interpolation 3D*)


(* ::Text:: *)
(*These functions make fast interpolation possible of an image on a grid.*)


Options[GridInterpolation3D]={InterpolationOrder->1};

SyntaxInformation[GridInterpolation3D] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

GridInterpolation3D[dat_List, grd_List, opts:OptionsPattern[]]:=
Module[{interFunction,data,grid,cgrid,ddim,gdim},
	data=ArrayPad[dat,1]//N;
	grid=Flatten[grd,2];
	ddim=Dimensions[data];
	gdim=Dimensions[grd];
	With[
		{
			(*Calls the interpolation function for grid coordinates*)
			InterpolateOnGrid=Compile[{{pts,_Real,2}},Map[interFunction[#[[1]],#[[2]],#[[3]]]&,pts]]
		},
		(*Define interpolation function*)
		interFunction=ListInterpolation[data,opts];
		
		(*actual grid interpolation*)
		cgrid=Transpose[Clip[grid[[All,#]]+1,{1,ddim[[#]]}]&/@{1,2,3}];
		Partition[Partition[InterpolateOnGrid[cgrid],gdim[[3]]], gdim[[2]]]
		]
	]


(* ::Subsubsection::Closed:: *)
(*Deform Data*)


SyntaxInformation[DefData] = {"ArgumentsPattern" -> {_, _, _}};

DefData[data_,vox_,w_]:=GridInterpolation3D[data,TransformGrid3D[GridData[data,1],vox,w,Dimensions[data]]]


(* ::Section:: *)
(*End Package*)


End[]

SetAttributes[#,{Protected, ReadProtected}]&/@ Names["DTITools`RegistrationTools`*"];

EndPackage[]
