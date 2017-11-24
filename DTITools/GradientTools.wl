(* ::Package:: *)

(* ::Title:: *)
(*DTITools GradientTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`GradientTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


GenerateGradients::usage = 
"GenerateGradients[numb] optimizes a set with numb gradients, numb mus be an integer.
GenerateGradients[{numb, fixed}] optimizes a set with numb gradients, numb must ba an integer and fixed a list of 3D coordiantes e.g. {{0,0,1},{0,1,0}}.\
The fixed gradients will not be moved.
GenerateGradients[{numb1, numb2 ...}, alpha] optimizes a multi shel gradient set with numb gradients per shel. If alpha is set to 0.5 equal importance is given to\
the optimal distribution of each shell en the enitre set. if alpha is 0 only the sub shels will be optimized, if alpha is set to 1 only the global set wil be optimized."

ConditionNumberCalc::usage = 
"ConditionNumberCalc[grads] calcualtes the condition number of the gradient set."

EnergyCalc::usage = 
"EnergyCalc[grads] calcualtes the total Energy of the gradient set."

OverPlusCalc::usage =
"OverPlusCalc[grads] determines the minimal overplus factor of of the gradient set."

GenerateGradientsGUI::usage =
"GenerateGradientsGUI[] runs the GenerateGradients function in GUI with output for the philips system."

FinalGrads::usage = 
"FinalGrads[grtxt,{int,intn},{rand,order}] finalizes the gradient txt file. 
grtxt is the output from the function ConvertGrads, which convert the grad to txt format.
int is True or False, if set to True it interleaves b=0 gradients every intn directions.
rand indicates if the gradients need to be randomized, for this it uses the order which is the output of FindOrder."

ConvertGrads::usage = 
"ConvertGrads[grad, bv] converts the gradients to txt format, which is needed for FinalGrads."

FindOrder::usage = 
"FindOrder[grad,bv] finds the optimal order of the gradient directions which minimizes the duty cycle.
The output is needed for FinalGrads.
grad is a list of gradient sets and bv is a list of b-values with the same number as the list of gradient sets."

GradRotMat::usage =
"GradRotMat[file] loads the rotation matrix from a .dcm file based on the imageorientation dicom field."

GradRot::usage =
"GradRot[\"DTIfile.dcm\",gradients] rotates the gradients with the rotationmatrix obtained from the header of the given dicomfile. uses GradRotMat."

GradDirs::usage = 
"GradDirs[grad] imports the gradient directions implemented on the philips defauls mri console and the patch by martijn froeling."

Bmatrix::usage =
"Bmatrix[bvec,grad] creates bmatrix form grad and bvec in form {-bxx, -byy, -bzz, -bxy, -bxz, -byz ,1}.
Bmatrix[{bvec,grad}] creates bmatrix form grad and bvec in form {bxx, byy, bzz, bxy, bxz, byz}."

BmatrixInv::usage = 
"BmatrixInv[bm] generates a bvecotr and gradiens directions form a given bmatrx.
BmatrixInv[bm, bvi] generates a bvecotr and gradiens directions form a given bmatrx using the given bvalues bvi."

BmatrixConv::usage = 
"BmatrixConv[bm] converts the bmatrix form 7 to 6 or from 6 to 7."

BmatrixRot::usage = 
"BmatrixRot[bmat, rotmat] Rotates the B-matrix."

BmatrixCalc::usage = 
"BmatrixCalc[\"folder\", grads] calculates the true bmatrix from the exported sequence parameters from the philips scanner that are stored in \"folder\" for each of the gradient directions grads."

BmatrixToggle::usage = 
"BmatrixToggle[bmat, axes, flip], axes can be any order of {\"x\",\"y\",\"z\"}. flip should be {1,1,1},{1,1,-1},{1,-1,1} or {-1,1,1}."

GradBmatrix::usage = 
"GradBmatrix[Gt, hw, te, t] Calculates the true bmatrix from the sequence created by GradSeq."

GradSeq::usage = 
"GradSeq[pars, t, grad] Creates a sequence from the gradient pars imported by ImportGradObj."

ImportGradObj::usage = 
"ImportGradObj[folder] Imports the gradient par files exported from the philips scanner."

GetSliceNormal::usage = 
"GetSliceNormal[file] imports the slice normal from a dicom image."

GetSliceNormalDir::usage = 
"GetSliceNormalDir[file] imports the slice normal from a enhanced dicom image."

CalculateMoments::usage = 
"CalculateMoments[{Gt, hw, te}, t] calculates the 0th to 3th order moments of the sequence created by GradSeq. Output is {{Gt, M0, M1, M2, M3}, vals}."

UniqueBvalPosition::usage = 
"UniqueBvalPosition[bval] generates a list of all the unique bvalues and their positions.
UniqueBvalPosition[bval, num] generates a list of all the unique bvalues and their positions that are present in the \
dataset equal or more than num times"

GetGradientScanOrder::usage = 
"GetGradientScanOrder[grad, bval] determines the scanorder based on the txt file provided to the scanner as input. 
GetGradientScanOrder[file, grad, bval] determines the scanorder based on the txt file provided to the scanner as input."

BvecCreate::usage = 
"BvecCreate[bval, grad] creates a bvector."


(* ::Subsection::Closed:: *)
(*Options*)


Steps::usage = 
"Steps is the number of step that is used in Generate Grads."

Runs::usage = 
"Runs is an option for GenerateGradients. Set how often the minimalization function is run. The best solution of all runs is the output. Default value is 1."

VisualOpt::usage = 
"VisualOpt is an option for GenerateGradients. Show the minimalization proces of eacht calculation step. Default is False."

GradType::usage = 
"GradType is what type of gradient set wil be produced in GenerateGradients \"Normal\" or \"OverPlus\"."

ConditionCalc::usage = 
"ConditionCalc is an option for GenerateGradients if set to true GenerateGradients will also give the condition number evolution of the system."

PhaseEncoding::usage = 
"PhaseEncoding is an options of GradSeq. Values can be \"A\", \"P\", \"R\" and \"L\"."

FlipAxes::usage =
"FlipAxes is an option for GradSeq. Defaul value is {{1,1,1},{1,1,1}}. First three values are for diffusion gradients last three are for the acquisition gradients."

SwitchAxes::usage =
"SwitchAxes is an option for GradSeq. Defaul value is {{1,2,3},{1,2,3}}. First three values are for diffusion gradients last three are for the acquisition gradients."

FullGrad::usage = 
"FullGrad is an option for Grad. Default is True. When true the gradient directions wil be loaded with the first gradient {0,0,0}."

UseGrad::usage = 
"UseGrad is an option for GradSeq. The default value is {0, 1, {1, 0}, 1} where {grex, gr180, {grepi1, grepi2}, grdiff, grflow}."

FlipGrad::usage = 
"FlipGrad is an option for GradSeq. When FlipGrad is true the gr180 is fliped."

UnitMulti::usage= 
"UnitMulti is an option for GradSeq. Defaul value is 10^-3. Defines the scaling of the gradient strength."

OutputType::usage = 
"OutputType is an option for BmatrixCalc. Values can be \"Matrix\" of \"Gradients\"."

OutputPlot::usage =
"OutputPlot is an option for GradBmatrix. It specifies if the plots of the gradients should also be exported."

StepSizeI::usage = 
"StepSizeI is an option for GradBmatrix. Specifies the integration stepsize is Method -> \"Numerical\" is used."

FullSphere::usage = 
"FullSphere is an option for GenerateGradients. If set True the gradients will be optimized on a full sphere rather than half a sphere."

OrderSpan::usage = 
"OrderSpan is an options for FindOrder."


(* ::Subsection::Closed:: *)
(*Error Messages*)


GradDirs::grad = 
"`1` is not a valid gradient specification. Posible values are:
	Philips normal: \"gr6\", \"gr15\" or \"gr32\".
	Philips over plus mode: \"gr6op\", \"gr15op\" or \"gr32op\".
	Patch normal mode: \"opt6\",\"opt10\",\"opt15\",\"opt24\",\"opt32\",\"opt48\" or \"opt60\".
	Patch over plus mode: \"opt6op\",\"opt10op\",\"opt15op\",\"opt24op\",\"opt32op\",\"opt48op\" or \"opt60op\"."

BmatrixToggle::axes = "input: `1` should contain permuations of {\"x\",\"y\",\"z\"}.";

BmatrixToggle::flip = "input:`1` should be be {1,1,1},{1,1,-1},{1,-1,1} or {-1,1,1}.";



(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*GenerateGradients*)


(* ::Subsubsection::Closed:: *)
(*GenerateGradients*)


Options[GenerateGradients] = {Steps -> 1000, Runs -> 1, 
   VisualOpt -> False, GradType -> "Normal", ConditionCalc -> False, 
   FullSphere -> False};
(*type\[Rule]Normal of OverPlus*)

SyntaxInformation[GenerateGradients] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

(*default*)
GenerateGradients[numb_Integer, opts : OptionsPattern[]] := 
 If[OptionValue[GradType] == "OverPlus", 
  GenerateGradientsi[{numb}, {}, 0, 0, opts, Method -> "OverPlus"], 
  If[OptionValue[GradType] == "Normal2", 
   GenerateGradientsi[{numb}, {}, 0, 0, opts, Method -> "Default2"], 
   GenerateGradientsi[{numb}, {}, 0, 0, opts, Method -> "Default"]]]

(*use fixed*)
GenerateGradients[inp : {_Integer, _List}, opts : OptionsPattern[]] :=
  If[OptionValue[GradType] == "OverPlus", 
  GenerateGradientsi[{inp[[1]]}, inp[[2]], 0, 0, opts, 
   Method -> "OverPlus"], 
  GenerateGradientsi[{inp[[1]]}, inp[[2]], 0, 0, opts, 
   Method -> "Fixed"]]

(*use shels*)
GenerateGradients[numbs : {_Integer ..}, opts : OptionsPattern[]] := 
 GenerateGradientsi[numbs, 0, 0.5, 0, opts, Method -> "Shels"]
GenerateGradients[numbs : {_Integer ..}, alpha_, 
  opts : OptionsPattern[]] := 
 GenerateGradientsi[numbs, 0, Clip[alpha, {0, 1}], 0, opts, 
  Method -> "Shels"]

(*Non rand*)
GenerateGradients[grds : {{_, _, _} ..}, opts : OptionsPattern[]] := 
 If[OptionValue[GradType] == "OverPlus", 
  GenerateGradientsi[{Length[grds]}, 0, 0, grds, opts, 
   Method -> "OverPlusNon"], 
  GenerateGradientsi[{Length[grds]}, 0, 0, grds, opts, 
   Method -> "NonRandom"]]


(* ::Subsubsection::Closed:: *)
(*GenerateGradientsi*)


Options[GenerateGradientsi] = {Steps -> 1000, Runs -> 1, 
   VisualOpt -> True, Method -> "Default", GradType -> "Normal", 
   ConditionCalc -> False, FullSphere -> False};

GenerateGradientsi[numbs_, fixed_, alph_, initp_, OptionsPattern[]] :=

 
 Block[{steps, runs, method, vis, cond, condtot, nf, ni, condnr, 
   velocity0, velocity1, parti, part, ns, plot, tempplot, vel, charge,
    tempc, sph, cols, output, half},
  DynamicModule[{vp, vv, va, tempp, points},
   (*Initialisation *)
   cols = {Red, Green, Blue, Yellow, Pink, Darker[Red], 
     Darker[Green], Darker[Blue], Darker[Yellow], Darker[Pink]};
   
   (*get options*)
   steps = OptionValue[Steps];
   runs = OptionValue[Runs];
   method = OptionValue[Method];
   vis = OptionValue[VisualOpt];
   cond = OptionValue[ConditionCalc];
   half = 1 - Boole[OptionValue[FullSphere]];
   
   (*determine length*)
   nf = Length[fixed];
   ni = Total[numbs] - nf;
   
   
   If[method == "OverPlus" || method == "OverPlusNon", half = 1];
   
   condnr = Infinity;
   If[method == "NonRandom" || method == "OverPlusNon",
    points = tempp = initp;
    runs = 1;
    ,
    points = tempp = ConstantArray[{0, 0, 0}, ni];
    ];
     
   (*get velocity matrix for multishell*)
   If[method == "Shels",
    ns = Length[numbs];
    {vel, part} = cc = Prepare[numbs, half, {}, alph][[2 ;;]]];
   
   (*visualisation*)
   If[OptionValue[VisualOpt],
    
    (*Initialize graphics*)
    vp = {1.3, -2.4, 2}; vv = {0, 0, 1}; va = 30 Degree;
    sph = Graphics3D[{White, Sphere[{0, 0, 0}, 0.95]}, 
      Lighting -> "Neutral", ImageSize -> 400, 
      PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}, {-1.1, 1.1}},
      ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
      ViewAngle -> Dynamic[va]];
    
    If[! (method == "Shels"),
     (*single shell*)
     plot = Row[{
         Show[sph, ListSpherePloti[Dynamic[points], Black, 0.05], 
          If[half == 1, ListSpherePloti[Dynamic[-points], Red, 0.05], 
           Graphics3D[]]],
         Show[sph, ListSpherePloti[Dynamic[tempp], Black, 0.05], 
          If[half == 1, ListSpherePloti[Dynamic[-tempp], Red, 0.05], 
           Graphics3D[]]]
         }];
     ,
     (*multi shell*)
     plot = Row[{
         Show[sph,
          
          MapThread[
           ListSpherePloti[#1, #2, 0.05] &, {Dynamic[points[[#]]] & /@ 
             part, cols[[1 ;; ns]]}],
          
          If[half == 1, ListSpherePloti[Dynamic[-points], Gray, 0.05], 
           Graphics3D[]]],
         Show[sph,
          
          MapThread[
           ListSpherePloti[#1, #2, 0.05] &, {Dynamic[tempp[[#]]] & /@ 
             part, cols[[1 ;; ns]]}],
          
          If[half == 1, ListSpherePloti[Dynamic[-tempp], Gray, 0.05], 
           Graphics3D[]]]
         }];
     ];
    tempplot = PrintTemporary[plot];
    ];
   
   (*Do number of runs*)
   condtot = Reap[Do[
       (*initialize*)
       If[! (method == "NonRandom" || method == "OverPlusNon"),
        tempp = RandInit[ni, half];
        If[! (method == "Fixed" || method == "OverPlus" || 
            method == "Default2"),
         tempp[[1]] = {0, 0, 1};
         ]
        ];
       
       tempp = Switch[
         OptionValue[Method],
         
         "Default2",(* default *)
         Do[tempp = GradOptimize1C[tempp, half];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         tempp
         ,
         
         "Default",(* default *)
         Do[tempp = GradOptimize2C[tempp, 1, half];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         tempp
          ,
         
         "NonRandom",(* non rand *)
         Do[tempp = GradOptimize1C[tempp, half];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         tempp
         ,
         
         "Fixed",(* fixed *)
         tempp = Join[fixed, tempp];
         Do[tempp = GradOptimize2C[tempp, nf, half];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         tempp
         ,
         
         "Shels",(* shells *)
         Do[tempp = GradOptimize4C[tempp, vel, half];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         tempp
         ,
         
          "OverPlus",(* overplus default *)
         tempp = Join[{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}, fixed, tempp];
         charge = 
          Join[ConstantArray[(.5 ni)^(1.2), 3], 
           ConstantArray[1, ni + nf]];
         
         Do[
          tempp = GradOptimize2C[tempp, nf + 3, 0], {Round[steps/10]}];
         Do[tempp = GradOptimize3C[tempp, charge, nf + 3];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         tempp = Drop[tempp, 3];
         Do[
          tempp = Normalize /@ 
            Clip[tempp, {-1/Sqrt[2], 1/Sqrt[2]}, {-1/Sqrt[2], 
              1/Sqrt[2]}], {25}];
         tempp
         ,
         
         "OverPlusNon",(* overplus default *)
         tempp = Join[{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}, tempp];
         charge = 
          Join[ConstantArray[(.5 ni)^(1.2), 3], 
           ConstantArray[1, ni + nf]];
         
         Do[tempp = GradOptimize3C[tempp, charge, nf + 3];
          If[cond, Sow[ConditionNumberCalc[tempp]]];
          , {steps}];
         
         tempp = Drop[tempp, 3];
         Do[
          tempp = Normalize /@ 
            Clip[tempp, {-1/Sqrt[2], 1/Sqrt[2]}, {-1/Sqrt[2], 
              1/Sqrt[2]}], {25}];
         tempp
         
         (*end switch*)];
       
       tempc = ConditionNumberCalc[tempp];
       
       If[tempc < condnr, points = tempp; condnr = tempc;];
       
       (*end runs do loop*)
       Pause[0.5];
       , {runs}]][[2]];
   NotebookDelete[tempplot];
   
   output = Chop[If[method == "Shels", points[[#]] & /@ part, points]];
   ];
  If[cond, {output, condtot[[1]]}, output]
  ]


(* ::Subsubsection::Closed:: *)
(*RandInit*)


(*random seed points on sphere*)
RandInit[ni_, half_] := If[half == 1,
   Sign[#[[3]] + 10.^-16] Normalize[#] & /@ 
    RandomReal[NormalDistribution[], {ni, 3}],
   Normalize[#] & /@ RandomReal[NormalDistribution[], {ni, 3}]
   ];


(* ::Subsubsection::Closed:: *)
(*ListSpherePloti*)


ListSpherePloti[pts_, col_, size_] := Graphics3D[{col, Sphere[#, size] & /@ pts}]


(* ::Subsubsection::Closed:: *)
(*ListStickPlot*)


ListStickPlot[pts_, size_] := Graphics3D[{Gray, Tube[0.98 {#, -#}, size] & /@ pts}]


(* ::Subsubsection::Closed:: *)
(*Prepare*)


(*initialization of points*)
Prepare[numbs_, half_] := Prepare[numbs, half, {}, 0];
Prepare[numbs_, half_, fixed_] := Prepare[numbs, half, fixed, 0];
Prepare[numbs_, half_, fixed_, alph_] := 
 Block[{nf, ni, ni2, parti, part, velocity0, velocity1, vel},
  nf = Length[fixed];
  ni = Total[numbs] - nf;
  
  (*only for multi shell numbs vector is longer than 1*)
  If[Length[numbs] > 1,
   (*case 1 multi shell*)
   parti = {1, 0} + # & /@ Drop[Partition[Prepend[Accumulate[numbs], 0], 2, 1, 1], -1];
   part = Range[#[[1]], #[[2]]] & /@ parti;
   
   ni2 = If[half == 0, ni, 2 ni];
   
   velocity0 = ConstantArray[0, {ni2, ni2, 3}];
   velocity1 = ConstantArray[1, {ni2, ni2, 3}];
   If[half == 1,
    Table[
        velocity0[[i, j]] = {1, 1, 1};
        velocity0[[i, j + ni]] = {1, 1, 1};
        velocity0[[i + ni, j]] = {1, 1, 1};
        velocity0[[i + ni, j + ni]] = {1, 1, 1};
        , {i, #1[[1]], #1[[2]]}, {j, #1[[1]], #1[[2]]}] & /@ parti;,
    Table[
        velocity0[[i, j]] = {1, 1, 1};
        , {i, #1[[1]], #1[[2]]}, {j, #1[[1]], #1[[2]]}] & /@ parti;
    ];
   vel = (1 - alph) velocity0 + alph velocity1;
   {RandInit[ni, half], vel, part}
   ,
   (*case 2 single shell*)
   If[fixed === {},
    (*no fixed*)
    RandInit[ni, half],
    (*fixed*)
    Join[fixed, RandInit[ni, half]]
    ]
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*GradOptimize functions*)


(*Generate cartesian grid*)
GradGrid[n_, full_] := Block[{points},
  points = If[EvenQ[n],
    If[full,
     (*even grid full*)
     Flatten[
       Table[{i, j, k}, {i, -1, 1, 2/(n - 1)}, {j, -1, 1, 
         2/(n - 1)}, {k, 1/(n - 1), 1, 2/(n - 1)}], 2] // N
     ,
     (*even spaced in between odd grid*)
     Flatten[
       Table[{i, j, k}, {i, -1 + 1/n, 1, 2/n}, {j, -1 + 1/n, 1, 
         2/n}, {k, 1/n, 1, 2/n}], 2] // N
     ]
    ,
    (*odd spaced*)
    points = 
     Flatten[Table[{i, j, k}, {i, -1, 1, 2/(n - 1)}, {j, -1, 1, 
         2/(n - 1)}, {k, 0, 1, 2/(n - 1)}], 2] // N;
    DeleteCases[
     If[(#[[1]] < 0. && ##[[3]] == 0.) || (#[[1]] == 0. && #[[2]] < 
            0 && ##[[3]] == 0.), Null, #] & /@ points, Null]
    ];
  points = Sort[points, Norm[#1] < Norm[#2] &] // N
  ]

(*optimize singel shell no fixed gradients*)
GradOptimize1C = Compile[{{points, _Real, 2}, {half, _Integer, 0}},
   Block[{n, n2, pointsi, pointsmat, distmatxyz, diag, distmat, 
     velocity, pointsnew, sign},
    n = Length[points];
    If[half == 1,
     pointsi = Join[points, -points]; n2 = 2 n;,
     pointsi = points; n2 = n
     ];
    pointsmat = ConstantArray[pointsi, n2];
    distmatxyz = pointsmat - Transpose[pointsmat];
    diag = DiagonalMatrix[ConstantArray[10.^32, n2]];
    distmat = Total[Transpose[distmatxyz, {2, 3, 1}]^2] + diag;
    velocity = Total[(distmatxyz/distmat)];
    pointsnew = (Normalize[#] & /@ (pointsi + velocity));
    If[half == 1, 
     pointsnew = Sign[pointsnew[[All, 3]] + 10.^-16] pointsnew;];
    pointsnew[[;; n]]]
   ];

(*optimize singel shell with fixed gradients*)
GradOptimize2C = 
  Compile[{{points, _Real, 2}, {nf, _Real, 0}, {half, _Integer, 0}},
   Block[{n, n2, pointsi, pointsmat, distmatxyz, diag, distmat, 
     velocity, pointsnew, rang},
    n = Length[points];
    pointsi = Join[points, -points];
    If[half == 1,
     pointsi = Join[points, -points]; n2 = 2 n;,
     pointsi = points; n2 = n
     ];
    pointsmat = ConstantArray[pointsi, n2];
    distmatxyz = pointsmat - Transpose[pointsmat];
    diag = DiagonalMatrix[ConstantArray[10.^16, n2]];
    distmat = Total[Transpose[distmatxyz, {2, 3, 1}]^2] + diag;
    velocity = Total[(distmatxyz/distmat)];
    If[half == 1,
     rang = Round[Join[Range[1, nf], Range[1 + n, n + nf]]];,
     rang = Round[Range[1, nf]];
     ];
    (velocity[[#]] = {0., 0., 0.}) & /@ rang;
    pointsnew = (Normalize[#] & /@ (pointsi + velocity));
    If[half == 1, 
     pointsnew = Sign[pointsnew[[All, 3]] + 10.^-16] pointsnew;];
    pointsnew[[;; n]]],
   {{rang, _Real, 1}}  
   ];

(*optimize singe shell Overplus*)
GradOptimize3C = 
  Compile[{{points, _Real, 2}, {char, _Real, 1}, {nf, _Real, 0}},
   Block[{n, pointsi, pointsmat, distmatxyz, diag, distmat, velocity, 
     pointsnew, chari, charmat, chars, rang},
    n = Length[points];
    pointsi = Join[points, -points];
    chari = Join[char, char];
    pointsmat = ConstantArray[pointsi, 2 n];
    charmat = ConstantArray[chari, 2 n];
    distmatxyz = pointsmat - Transpose[pointsmat];
    chars = charmat*Transpose[charmat];
    diag = DiagonalMatrix[ConstantArray[10.^32, 2 n]];
    distmat = Total[Transpose[distmatxyz, {2, 3, 1}]^2] + diag;
    velocity = Total[chars (distmatxyz/distmat)];
    rang = Round[Join[Range[1, nf], Range[1 + n, n + nf]]];
    velocity[[rang]] = {0., 0., 0.};
    pointsnew = (Sign[#[[3]] + 10.^-16]*Normalize[#] & /@ (pointsi + 
         velocity));
    pointsnew[[;; n]]],
   {{rang, _Real, 1}}
   ];

(*optimize multi shell*)
GradOptimize4C = 
  Compile[{{points, _Real, 2}, {vel, _Real, 3}, {half, _Integer, 0}},
   Block[{n, n2, pointsi, pointsmat, distmatxyz, diag, distmat, 
     velocity, pointsnew},
    n = Length[points];
    If[half == 1,
     pointsi = Join[points, -points]; n2 = 2 n;,
     pointsi = points; n2 = n;
     ];
    pointsmat = ConstantArray[pointsi, n2];
    distmatxyz = pointsmat - Transpose[pointsmat];
    diag = DiagonalMatrix[ConstantArray[10.^16, n2]];
    distmat = Total[Transpose[distmatxyz, {2, 3, 1}]^2] + diag;
    velocity = Total[vel (distmatxyz/distmat)];
    pointsnew = (Normalize[#] & /@ (pointsi + velocity));
    If[half == 1, 
     pointsnew = Sign[pointsnew[[All, 3]] + 10.^-16] pointsnew;];
    pointsnew[[;; n]]]
   ];


(* ::Subsection::Closed:: *)
(*GenerateGradientsGUI*)


GenerateGradientsGUI[] := Block[{pan},
  pan = Manipulate[
    (*covert bval and names*)
    bvall = If[NumberQ[bvall], {bvall}, bvall];
    names = StringReplace[#, " " -> "_"] & /@ names;
    
    (*constrain one shel to 124 directions*)
    dirs1 = Clip[dirs1, {3, 128}];
    
    (*Constrain, multi shel to 60 per shell*)
    {dirs21, dirs22, dirs23, dirs24, dirs25, dirs26} = 
     Clip[{dirs21, dirs22, dirs23, dirs24, dirs25, dirs26}, {3, 128}];
    dirs2 = {dirs21, dirs22, dirs23, dirs24, dirs25, dirs26}[[
      1 ;; nshels]];
    
    (*convert gradients to readable format*)
    out = names[[mult]] <> Switch[mult,
       (*DWI*)4, 
       If[gradd === "", gradd, 
        FinalGrads[outd, {inter, int}, {random, orderd}]],
       (*cartesian grid*)3, 
       If[gradc === "", gradc, 
        FinalGrads[outc, {inter, int}, {random, orderc}]],
       (*multi shel*)2, 
       If[gradm === "", gradm, 
        FinalGrads[outm, {inter, int}, {random, orderm}]],
       (*single shel*)1, 
       If[grads === "", grads, 
        FinalGrads[outs, {inter, int}, {random, orders}]]];
    
    (*The display, show the plots or show the Gradient directions txt*)

    
    Dynamic[
     Switch[disp,
      (*display as 3D plot*)
      1,
      Switch[mult,
       1,(*show single shell gradient directions*)
       If[points === {},
        (*no gradients present only show sphere*)
        Column[{"", 
          Show[SpherePlot[1, opacity], ImageSize -> size, 
           PlotLabel -> ""]}],
        pointspl = 
         If[proj, 
          Map[Sign[#[[3]] + 10.^-16]*Normalize[#] &, points, {1}], 
          points];
        Column[{"",
          Show[
           SpherePlot[1, opacity],
           (*Shell positive z*)
           ListSpherePloti[pointspl, Red, 0.05],
           (*Shells Negative z*)
           
           If[mirror, ListSpherePloti[-pointspl, Gray, 0.05], 
            Graphics3D[{}]],
           (*sticks*)
           If[sticks, ListStickPlot[pointspl, 0.01], Graphics3D[{}]],
           ImageSize -> size, PlotLabel -> "", Background -> app
           ]
          }, Alignment -> Center]
        ],
       2,(*show multi shell gradient directions*)
       If[ppoints === {},
        (*no gradients present only show sphere*)
        Column[{"", 
          Show[SpherePlot[1, opacity], ImageSize -> size, 
           PlotLabel -> ""]}]
        ,
        (*gradients are defiended *)
        ppointspl = 
         If[proj, 
          Map[Sign[#[[3]] + 10.^-16]*Normalize[#] &, ppoints, {2}], 
          ppoints];
        len = Length[ppoints];
        rlen = Range[len];
        If[running, show = rlen];
        Column[{
          
          SetterBar[Dynamic[show], 
           Join[{0 -> "multi", rlen -> "all"}, rlen]],
          If[show === 0,
           (*show all gradients multi shell*)
           Show[
            (*Multi shells positive z*)
            MapThread[{
                SpherePlot[#3, opacity],
                ListSpherePloti[#1 #3, #2, 0.05]} &, {Reverse@
                ppointspl, 
               Reverse@{Red, Green, Blue, Yellow, Pink, Purple}[[
                 1 ;; len]], Range[1, .5, -.5/(len - 1)]
               }] // Flatten,
            (*multi shells negative z*)
            
            If[mirror, 
             MapThread[
              ListSpherePloti[-#1 #2, Gray, 0.05] &, {ppointspl, 
               Range[1, .5, -.5/(len - 1)]}], Graphics3D[{}]],
            ImageSize -> size, Background -> app
            ]
           ,
           (*show all gradients single shell*)
           Show[
            SpherePlot[1, opacity],
            (*Multi shells positive z*)
            
            MapThread[
              ListSpherePloti[#1, #2, 
                0.05] &, {ppointspl, {Red, Green, Blue, Yellow, Pink, 
                 Purple}[[1 ;; len]]}][[Clip[show, {1, len}]]],
            (*multi shells negative z*)
            
            If[mirror, (ListSpherePloti[-#, Gray, 0.05] & /@ 
                ppointspl)[[Clip[show, {1, len}]]], Graphics3D[{}]],
            (*sticks*)
            
            If[sticks, (ListStickPlot[#, 0.01] & /@ ppointspl)[[
              Clip[show, {1, len}]]], Graphics3D[{}]],
            ImageSize -> size, Background -> app
            ]
           ]
          }, Alignment -> Center]
        ],
       3,(*show cartesian grid directions*)
       Column[{"",
         Show[
          SpherePlot[0, opacity],
          (*Shell positive z*)
          ListSpherePloti[pointsc, Red, 0.05],
          (*Shells Negative z*)
          
          If[mirror, ListSpherePloti[-pointsc, Gray, 0.05], 
           Graphics3D[{}]]
          ,
          ImageSize -> size, PlotLabel -> "", Background -> app
          ]
         }],
       4,(*DWI*)
       Show[
        ListLinePlot[{{0, #}, {1, #}} & /@ Prepend[bvald, 0],
         PlotStyle -> Directive[{Gray, Thick}], 
         PlotRange -> {{0, 2}, Full},
         AxesStyle -> Thick, Axes -> {False, True}, 
         AspectRatio -> 1.4, 
         AxesLabel -> {None, Style["b-value", Bold, Black]}],
        ListPlot[{1, #} & /@ Prepend[bvald, 0], PlotStyle -> Red]
        ]
       ],
      
      (*display as polar plot*)
      2,
      Column[{
        (*controls*)
        Row[{
          
          Slider2D[Dynamic[viewvec], {{-1, -1}, {1, 1}}, 
           ContinuousAction -> True],
          
          SetterBar[
           Dynamic[ctype], {1 -> "polar", 2 -> "even grid", 
            3 -> "scaled grid"}]
          }, ImageSize -> 400, Alignment -> Center],
        If[mult == 2, 
         SetterBar[Dynamic[showc], Join[{All -> "all"}, rlenc]], ""],
        (*the plot*)
        Dynamic@Show[
          charts[[ctype]],
          (*change between single and multi slice*)
          Switch[mult,
           1,
           If[points === {},
            Graphics[],
            
            pointspl = 
             If[proj, 
              Map[Sign[#[[3]] + 10.^-16]*Normalize[#] &, points, {1}],
               points];
            PlotChartPoints[pointspl, {mirror, Red}]
            ],
           2,
           If[ppoints === {},
            Graphics[],
            
            ppointspl = 
             If[proj, 
              Map[Sign[#[[3]] + 10.^-16]*Normalize[#] &, 
               ppoints, {2}], ppoints];
            len = Length[ppoints];
            rlenc = Range[len];
            
            PlotChartPoints[
             ppointspl[[
              showc]], {mirror, {Red, Green, Blue, Yellow, Pink, 
                Purple}[[showc]]}]
            ],
           _,
           Graphics[]
           ]
          ]
        }, Alignment -> Center],
      
      (*display as text*)
      3,
      out,
      
      (*display G load*)
      4,
      Switch[mult,
       1, 
       If[orders === "", Graphics[], 
        PlotDuty[{ConstantArray[grads, Length[bvall]], bvall, orders},
          random]],
       2, 
       If[orderm === "", Graphics[], 
        PlotDuty[{gradm, bvals[[;; nshels]], orderm}, random]],
       3, 
       If[orderc === "", Graphics[], 
        PlotDuty[{gradc, bvalc, orderc}, random]],
       4, 
       If[orderd === "", Graphics[], 
        PlotDuty[{ConstantArray[gradd, Length[bvald]], bvald, orderd},
          random]]
       ]
      ]
     ]
    ,
    
    (*Controls*)
    (*set Name*)
    Row[{"      Set Name     ", 
      InputField[Dynamic[names[[mult]]], String]}],
    
    Delimiter,
    
    (*display controls*)
    {{disp, 1, "display gradients"}, {1 -> "graphics", 2 -> "chart", 
      3 -> "text", 4 -> "G load"}},
    {{opacity, 0.5, "sphere opacity"}, 0, 1, .1, 
     ControlType -> Slider},
    Row[{
      " sticks: ", Checkbox[Dynamic[sticks]],
      "   mirror grad.: ", Checkbox[Dynamic[mirror]],
      "   project grad. on half: ", Checkbox[Dynamic[proj]]
      }],
    Grid[{
      {
       Button["top", vp = {0, 0, 3.38}, ImageSize -> {50, 20}, 
        FrameMargins -> 0],
       Button["right", vp = {3.38, 0, 0}, ImageSize -> {50, 20}, 
        FrameMargins -> 0],
       Button["front", vp = {0, 3.38, 0}, ImageSize -> {50, 20}, 
        FrameMargins -> 0],
       Button[
        "reset view point", {vp, vv, va} = {{1.3, -2.4, 2}, {0, 0, 1},
           30. Degree}, ImageSize -> {100, 20}, FrameMargins -> 0]}
      (*,
      {
      Button["bottom",vp={0,0,-3.38},ImageSize\[Rule]{50,20},
      FrameMargins\[Rule]0],
      Button["left",vp={-3.38,0,0},ImageSize\[Rule]{50,20},
      FrameMargins\[Rule]0],
      Button["back",vp={0,-3.38,0},ImageSize\[Rule]{50,20},
      FrameMargins\[Rule]0]
      }*)
      }],
    
    Delimiter,
    
    (*multi or single shell*)
    {{half, 1, "Full or half sphere"}, {1 -> "half sphere", 
      0 -> "full sphere"}},
    {{mult, 1, "shells"},
     {1 -> "single shell", 2 -> "multi shell", 3 -> "cartesian", 
      4 -> "DWI"},
     ControlType -> PopupMenu, FieldSize -> {13, 0.7}},
    
    (*single multi and cartesian controls*)
    PaneSelector[{
      (*Singel shell pannel*)
      1 -> Column[{
         Control[{{type, "normal", "type"}, {"normal", 
            "normal fixed z", "normal fixed x, y and z",
            "over-plus", "over-plus fixed z", 
            "over-plus fixed x, y and z"}, ControlType -> PopupMenu, 
           FieldSize -> {13, 0.7}}],
         Control[{{dirs1, 30, "number of gradients"}, 6, 128, 1, 
           Appearance -> "Labeled", ImageSize -> Tiny, 
           AppearanceElements -> {"InputField"}}]
         }]
      ,
      (*Multi shell pannel*)
      2 -> Column[{
         Row[{Control[{{nshels, 2, "number of shells"}, {2, 3, 4, 5, 
              6}}],
           
           Control[{{shel, 1,     "     shell"}, 
             Dynamic[Range[nshels]], ControlType -> SetterBar}]
           }],
         PaneSelector[{
           
           Control[{{dirs21, 15, "number of gradients shell 1"}, 3, 
             128, 1, Appearance -> "Labeled", ImageSize -> Tiny, 
             AppearanceElements -> {"InputField"}}],
           
           Control[{{dirs22, 15, "number of gradients shell 2"}, 3, 
             128, 1, Appearance -> "Labeled", ImageSize -> Tiny, 
             AppearanceElements -> {"InputField"}}],
           
           Control[{{dirs23, 15, "number of gradients shell 3"}, 3, 
             128, 1, Appearance -> "Labeled", ImageSize -> Tiny, 
             AppearanceElements -> {"InputField"}}],
           
           Control[{{dirs24, 15, "number of gradients shell 4"}, 3, 
             128, 1, Appearance -> "Labeled", ImageSize -> Tiny, 
             AppearanceElements -> {"InputField"}}],
           
           Control[{{dirs25, 15, "number of gradients shell 5"}, 3, 
             128, 1, Appearance -> "Labeled", ImageSize -> Tiny, 
             AppearanceElements -> {"InputField"}}],
           
           Control[{{dirs26, 15, "number of gradients shell 6"}, 3, 
             128, 1, Appearance -> "Labeled", ImageSize -> Tiny, 
             AppearanceElements -> {"InputField"}}]
           }, Dynamic[shel]],
         Control[{{weight, .5, "shell weighting"}, 0, 1, .05, 
           ControlType -> Slider, Appearance -> "Labeled", 
           ImageSize -> Tiny, AppearanceElements -> {"InputField"}}]
         }]
      ,
      (*cartesian grid*)
      3 -> Column[{
         Control[{{grid, 9, "cartesian grid size"}, 5, 15, 1, 
           Appearance -> "Labeled", ImageSize -> Tiny}],
         Control[{{gridf, True, 
            "full even grid"}, {False -> "no (in between odd grid)", 
            True -> "yes"}}]
         }],
      (*DWI pannel*)
      4 -> 
       Control[{{typed, "normal", "               type"}, {"normal", 
          "over-plus"}, ControlType -> SetterBar}]
      }, mult],
    
    Delimiter,
    
    (*input bvals*)
    PaneSelector[
     {
      1 -> 
       Row[{"      b-value:   ", InputField[Dynamic[bvall], Expression,
          
          Background -> 
           Dynamic[
            If[((AllTrue[bvall, NumberQ] && ListQ[bvall]) || 
               NumberQ[bvall]), None, Lighter[Lighter[Red]]]]
          ]}],
      2 -> 
       Dynamic[Grid[
         Partition[
          PadRight[
           Row[{"b-val" <> ToString[#] <> ":", 
               InputField[Dynamic[bvals[[#]]], Number, 
                FieldSize -> 5]}] & /@ Range[1, nshels], 6, ""], 3]]],
      3 -> 
       Row[{"max b (corner):   ", InputField[Dynamic[bvalc], Number]}],
      4 -> 
       Row[{"      b-value:   ", InputField[Dynamic[bvald], Expression,
          
          Background -> 
           Dynamic[
            If[((AllTrue[bvald, NumberQ] && ListQ[bvald]) || 
               NumberQ[bvald]), None, Lighter[Lighter[Red]]]]
          ]}]
      }, mult],
    
    Delimiter,
    
    (*interleave b=0 and optimize gradient load*)
    Row[{"      interleave b=0: " Checkbox[Dynamic[inter]], 
      "      Optimize G load: ", Checkbox[Dynamic[random]]}],
    {{int, 10, "interleave b=0 every: "}, Range[3, 20]},
    
    Delimiter,
    
    (*quality and multi or single shell*)
    {{steps, 1000, "quality (iterations)"}, {500 -> "poor (500)", 
      1000 -> "normal (1000)", 2500 -> "excellent (2500)", 
      5000 -> "perfect (5000)", 10000 -> "extreme (10000)"}
     , ControlType -> PopupMenu, FieldSize -> {9, 0.7}},
    (*generate gradietns button*)
    (*mult:1-single shell; 2-multi shell; 3-cartesian;4-DWI*)
    Row[{Button["generate gradients",
       (*initiate gradient generation*)
       app = Lighter[Lighter[LightGray]];
       disp = 1; running = True;
       mirror = If[half == 0, False, True];
       proj = False;
       Pause[.2];
       
       (*switch to correct method*)
       Switch[mult
        , 4,(*DWI*)
        gradd = If[typed == "normal",
           {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
           {{-0.707107, -0.5, 0.5}, {0.707107, -0.5, 0.5}, {0., 
             0.707107, 0.707107}}
           ];
        , 3, (*cartesian*)
        gradc = pointsc = GradGrid[grid, gridf];
        
        , 2,(*multi shell*)
        {mpoints, vel, part} = Prepare[dirs2, half, {}, weight];
        Pause[.5];
        Do[
         mpoints = GradOptimize4C[mpoints, vel, half];
         ppoints = Chop[mpoints[[#]]] & /@ part;
         , {steps}];
        gradm = Chop[mpoints[[#]]] & /@ part;
        
        , 1,(*single shell, normal or overplus*)
        grads = Switch[type
          , "normal",
          points = Prepare[{dirs1}, half];
          Pause[.5];
          Do[points = GradOptimize1C[points, half], {steps}];
          Chop[points]
          
          , "normal fixed z",
          points = Prepare[{dirs1}, half, {{0, 0, 1}}];
          Pause[.5];
          Do[points = GradOptimize2C[points, 1, half], {steps}];
          Chop[points]
          
          , "normal fixed x, y and z",
          
          points = 
           Prepare[{dirs1}, half, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}];
          Pause[.5];
          Do[points = GradOptimize2C[points, 3, half], {steps}];
          Chop[points]
          
          (*all over plus options only work on half shell*)
          , "over-plus",
          half = 1;
          points = Prepare[{dirs1}, half];
          points = Join[{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}, points];
          Pause[.5];
          (*quick distribution for init*)
          
          Do[points = 
            GradOptimize2C[points, 3, half], {Round[steps/10]}];
          (*optimize overplys*)
          
          charge = 
           Join[ConstantArray[(.5 dirs1)^(1.2), 3], 
            ConstantArray[1, dirs1]];
          Do[points = GradOptimize3C[points, charge, 3], {steps}];
          points = Chop[Drop[points, 3]]
          
          , "over-plus fixed z",
          half = 1;
          points = Prepare[{dirs1 - 1}, half];
          points = Join[{
             {0, 0, 1}, {0, 1, 0}, {1, 0, 0},
             {0., 0.707107, 0.707107}
             }, points];
          Pause[.5];
          (*quick distribution for init*)
          
          Do[points = 
            GradOptimize2C[points, 4, half], {Round[steps/10]}];
          (*optimize overplys*)
          
          charge = 
           Join[ConstantArray[(.5 dirs1)^(1.2), 3], 
            ConstantArray[1, dirs1]];
          Do[points = GradOptimize3C[points, charge, 4], {steps}];
          points = Chop[Drop[points, 3]]
          
          , "over-plus fixed x, y and z",
          half = 1;
          points = Prepare[{dirs1 - 3}, half];
          points = Join[{
             {0, 0, 1}, {0, 1, 0}, {1, 0, 0},
             {-0.707107, -0.5, 0.5}, {0.707107, -0.5, 0.5}, {0., 
              0.707107, 0.707107}
             }, points];
          Pause[.5];
          (*quick distribution for init*)
          
          Do[points = 
            GradOptimize2C[points, 6, half], {Round[steps/10]}];
          (*optimize overplys*)
          
          charge = 
           Join[ConstantArray[(.5 dirs1)^(1.2), 3], 
            ConstantArray[1, dirs1]];
          Do[points = GradOptimize3C[points, charge, 6], {steps}];
          points = Chop[Drop[points, 3]]
          
          ];(*end single shel options*)
        
        ];(*end generate grads*)
       
       (*covert gradients to txt and find optimal order dutycycle*)
       (*mult:1-single shell; 2-multi shell; 3-cartesian;4-DWI*)
       Switch[mult,
        1,
        outs = 
         ConvertGrads[ConstantArray[grads, Length[bvall]], bvall];
        orders = FindOrder[ConstantArray[grads, Length[bvall]], bvall];
        , 2,
        outm = ConvertGrads[gradm, bvals[[;; nshels]]];
        orderm = FindOrder[gradm, bvals[[;; nshels]]];
        , 3,
        outc = ConvertGrads[gradc, {bvalc}];
        orderc = FindOrder[gradc, bvalc];
        , 4,
        outd = 
         ConvertGrads[ConstantArray[gradd, Length[bvald]], bvald];
        outd[[2]] = Join[
          If[typed == "normal",
           {" 1.00000    0.00000    0.00000       0.1", 
            " 0.00000    1.00000    0.00000       0.1", 
            " 0.00000    0.00000    1.00000       0.1"},
           {"-0.70711   -0.50000    0.50000       0.1", 
            " 0.70711   -0.50000    0.50000       0.1", 
            " 0.00000    0.70711    0.70711       0.1"}],
          {
           " 0.02704    0.79706    0.60330       0.1", 
           "-0.09999   -0.59783    0.79536       0.1", 
           " 0.23191   -0.77261    0.59101       0.1", 
           " 0.52867   -0.79903    0.28646       0.1", 
           "-0.18297   -0.98140    0.05818       0.1", 
           "-0.86286    0.19578    0.46599       0.1", 
           " 0.05126    0.20181    0.97808       0.1", 
           "-0.66047    0.01890    0.75062       0.1", 
           " 0.65426   -0.18818    0.73249       0.1", 
           "-0.33052    0.11752    0.93646       0.1", 
           "-0.95407   -0.14465    0.26236       0.1", 
           "-0.14402   -0.86963    0.47224       0.1", 
           "-0.78028   -0.56460    0.26904       0.1", 
           "-0.75319   -0.31631    0.57675       0.1", 
           "-0.41392   -0.31357    0.85460       0.1", 
           " 0.56536    0.27275    0.77845       0.1", 
           "-0.73743    0.58153    0.34354       0.1", 
           " 0.29639   -0.44648    0.84428       0.1", 
           "-0.38689    0.75986    0.52243       0.1", 
           "-0.20081    0.52214    0.82888       0.1", 
           " 0.79424   -0.60284    0.07593       0.1", 
           " 0.98777    0.03678    0.15153       0.1", 
           " 0.29973    0.90469    0.30281       0.1", 
           "-0.95633    0.27860    0.08835       0.1", 
           "-0.49976   -0.83556    0.22819       0.1", 
           " 0.90574    0.41773    0.07170       0.1", 
           " 0.77284    0.45698    0.44032       0.1", 
           " 0.35339   -0.03796    0.93471       0.1", 
           " 0.84918    0.07921    0.52213       0.1", 
           "-0.50085    0.85590    0.12880       0.1", 
           " 0.47526    0.67317    0.56654       0.1", 
           " 0.87597   -0.32192    0.35923       0.1", 
           " 0.64492    0.75075    0.14303       0.1", 
           "-0.12687    0.96109    0.24536       0.1", 
           " 0.60869   -0.54852    0.57325       0.1", 
           "-0.03797   -0.21873    0.97505       0.1", 
           "-0.47780   -0.64010    0.60165       0.1", 
           "-0.57319    0.42677    0.69952       0.1", 
           " 0.18760   -0.96394    0.18874       0.1", 
           " 0.23342    0.50887    0.82859       0.1"}
          ];
        orderd = FindOrder[ConstantArray[gradd, Length[bvald]], bvald];
        ];
       
       (*stop gray background*)
       Pause[0.1];
       running = False;
       app = White;
       
       , Method -> "Queued", ImageSize -> {120, 23}],
      
      (*output buttens*)
      Button["to clipboard", CopyToClipboard[out], 
       ImageSize -> {80, 23}],
      Button["to file",
       file = SystemDialogInput["FileSave", "dti_vectors_input.txt"];
       If[! (file === $Canceled), Export[file, out, "Text"]], 
       ImageSize -> {80, 23}, Method -> "Queued"]}],
    
    (*disclaimer*)
    Delimiter,
    Row[{Style[
       "Made by Martijn Froeling, Phd \nm.froeling@umcutrecht.nl", \
{Small, Gray}]}],
    
    (*hidden dynamic local variables, using now control type*)
    {{points, {}}, ControlType -> None},
    {{pointspl, {}}, ControlType -> None},
    {{pointsc, {}}, ControlType -> None},
    {{ppoints, {}}, ControlType -> None},
    {{ppointspl, {}}, ControlType -> None},
    {{mpoints, {}}, ControlType -> None},
    {{gradd, ""}, ControlType -> None},
    {{gradm, ""}, ControlType -> None},
    {{grads, ""}, ControlType -> None},
    {{gradc, ""}, ControlType -> None},
    {{outd, ""}, ControlType -> None},
    {{outc, ""}, ControlType -> None},
    {{outm, ""}, ControlType -> None},
    {{outs, ""}, ControlType -> None},
    {{orderd, ""}, ControlType -> None},
    {{orderc, ""}, ControlType -> None},
    {{orderm, ""}, ControlType -> None},
    {{orders, ""}, ControlType -> None},
    {{out, ""}, ControlType -> None},
    {{file, ""}, ControlType -> None},
    
    {{show, {1, 2}}, ControlType -> None},
    {{showc, All}, ControlType -> None},
    {{weight, 0.5}, ControlType -> None},
    {{grid, 9}, ControlType -> None},
    {{gridf, False}, ControlType -> None},
    {{running, False}, ControlType -> None},
    {{app, White}, ControlType -> None},
    {{size, 430}, ControlType -> None},
    {{inter, True}, ControlType -> None},
    {{random, True}, ControlType -> None},
    {{sticks, False}, ControlType -> None},
    {{mirror, True}, ControlType -> None},
    {{proj, False}, ControlType -> None},
    
    {{vel, 1}, ControlType -> None},
    {{part, 1}, ControlType -> None},
    
    {{dirs1, 30}, ControlType -> None},
    {{dirs21, 15}, ControlType -> None},
    {{dirs22, 15}, ControlType -> None},
    {{dirs23, 15}, ControlType -> None},
    {{dirs24, 15}, ControlType -> None},
    {{dirs25, 15}, ControlType -> None},
    {{dirs26, 15}, ControlType -> None},
    
    {{bvald, {10, 20, 30, 40, 60, 80, 100, 200, 300, 500, 700, 1000}},
      ControlType -> None},
    {{bvall, {1000}}, ControlType -> None},
    {{bvals, Range[1000, 6000, 1000]}, ControlType -> None},
    {{bvalc, 9000}, ControlType -> None},
    
    {dirs2, ControlType -> None},
    {type, ControlType -> None},
    {typed, ControlType -> None},
    {sc, ControlType -> None},
    {scc, ControlType -> None},
    {shel, ControlType -> None},
    {{nshels, 2}, ControlType -> None},
    {len, ControlType -> None},
    {rlen, ControlType -> None},
    {{rlenc, {}}, ControlType -> None},
    {charge, ControlType -> None},
    
    {{names, {"Set_Name", "Shells_Name", "Grid_Name", "DWI_Name"}}, 
     ControlType -> None},
    
    {{vp, {1.3, -2.4, 2}}, ControlType -> None},
    {{va, 30. Degree}, ControlType -> None},
    {{vv, {0, 0, 1}}, ControlType -> None},
    
    {charts, ControlType -> None},
    {{viewvec, {0, 0}}, ControlType -> None},
    {{ctype, 1}, ControlType -> None},
    
    (*Manipulate settings*)
    ContentSize -> {450, 510},
    SaveDefinitions -> True,
    ControlPlacement -> Left,
    ContinuousAction -> True,
    Initialization :> {
      bvall = {1000},
      MakeChart[type_] := Module[{ranX, ranY, coors, lab},
        (*define grid*)
        ranX = Range[-180, 180, 30];
        ranY = Range[-90, 90, 10];
        (*define coordinate system*)
        coors = {
           
           N@Table[{j Cos[i Degree], 90 Sin[i Degree]}, {j, ranX}, {i,
               ranY}],
           N@Table[{j, i}, {j, ranX}, {i, ranY}],
           N@Table[{j, 90 Sin[i Degree]}, {j, ranX}, {i, ranY}]}[[
          type]];
        lab = {{"\[Phi] Cos[\[Theta]] (\[Degree])", 
            "\[Theta] Sin[\[Theta]] (\[Degree])"}, {"\[Phi] \
(\[Degree])", "\[Theta] (\[Degree])"}, {"\[Phi] (\[Degree])", 
            "\[Theta] Sin[\[Theta]] (\[Degree])"}}[[type]];
        Graphics[{
          {Lighter@Lighter@Lighter@Gray, 
           Polygon@Join[First@coors, Reverse@Last@coors]},
          {Lighter@Gray, Line[#]} & /@ coors,
          {Lighter@Gray, Line[#]} & /@ Transpose[coors]
          }, AspectRatio -> .8, ImageSize -> 400,
         PlotRange -> {{-185, 185}, {-95, 95}},
         LabelStyle -> 
          Directive[{Bold, Black, Medium, FontFamily -> "Helvetica"}],
         Frame -> True, FrameStyle -> Thick,
         Axes -> True, AxesStyle -> Thick,
         FrameTicks -> {{(Thread[{Round@coors[[1, All, 2]], ranY}][[
              1 ;; ;; 3]]), None}, {ranX, None}},
         FrameLabel -> lab
         ]
        ],
      charts = MakeChart /@ {1, 2, 3},
      SpherePlot[size_, op_] := If[size == 0 || size == 0.,
        Graphics3D[{},
         Lighting -> "Neutral", 
         PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}, {-1.1, 1.1}},
         ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
         ViewAngle -> Dynamic[va],
         SphericalRegion -> True]
        ,
        Graphics3D[{White, Opacity[op], Sphere[{0, 0, 0}, 0.95 size]},
         Lighting -> "Neutral",
         PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}, {-1.1, 1.1}},
         ViewPoint -> Dynamic[vp], ViewVertical -> Dynamic[vv], 
         ViewAngle -> Dynamic[va],
         SphericalRegion -> True]
        ],
      PlotChartPoints[grad_, {mirr_, col_}] := Block[{style},
        style = 
         If[ListQ[col], (Directive[#, PointSize[Large]] & /@ col), 
          Directive[col, PointSize[Large]]];
        Show[
         If[grad === {},
          (*if no gr then only chart*)
          Graphics[],
          Show[ListPlot[(*plot the gradients*)
            If[ArrayDepth[grad] == 2,
             CalcPolarPts[grad, ctype, viewvec],
             CalcPolarPts[#, ctype, viewvec] & /@ grad
             ], PlotStyle -> style],
           If[! mirr,(*plot the mirrored gradients*)
            Graphics[],
            ListPlot[If[ArrayDepth[grad] == 2,
              CalcPolarPts[-grad, ctype, viewvec],
              CalcPolarPts[-Flatten[grad, 1], ctype, viewvec]
              ], PlotStyle -> {Darker@Gray, PointSize[Large]}]
            ]]]]
        ]
      }
    , AppearanceElements -> None,
    AutorunSequencing -> {1}
    ];
  NotebookClose[gradwindow];
  gradwindow=CreateWindow[DialogNotebook[{CancelButton["Close",
  DialogReturn[]],pan},WindowSize->All,
  WindowTitle->"Generate gradients"]];
  ]


(* ::Subsubsection:: *)
(*CalcPolorPts*)


(*calculate polar coordinates from gradients*)
CalcPolarPts[grad_, type_] := CalcPolarPts[grad, type, {0., 0.}]
CalcPolarPts[grad_, type_, vec_] := Block[{sig, phi, gradp, rot},
  (*rotate viewpoint*)
  rot = (RotationMatrix[
      vec[[1]] 180 Degree, {0, 0, -1}].RotationMatrix[
      vec[[2]] 90 Degree, {0, 1, 0}]);
  (*calculate spherical coordinates*)
  {sig, phi} = 
   Transpose@(Quiet[({90, 
            0} + {-1, 1} ToSphericalCoordinates[rot.#][[2 ;;]]/
             Degree) & /@ grad] /. Indeterminate -> 0);
  (*transform to correct axes system*)
  {
    Transpose[{phi Cos[sig Degree], 90 Sin[sig Degree]}],
    Transpose[{phi , sig}],
    Transpose[{phi , 90 Sin[sig Degree]}]
    }[[type]]
  ]


(* ::Subsubsection:: *)
(*MakeChart*)


(*make grids*)
MakeChart[type_] := Module[{ranX, ranY, coors, lab},
   (*define grid*)
   ranX = Range[-180, 180, 30];
   ranY = Range[-90, 90, 10];
   (*define coordinate system*)
   coors = {
      N@Table[{j Cos[i Degree], 90 Sin[i Degree]}, {j, ranX}, {i, 
         ranY}],
      N@Table[{j, i}, {j, ranX}, {i, ranY}],
      N@Table[{j, 90 Sin[i Degree]}, {j, ranX}, {i, ranY}]}[[type]];
   lab = {{"\[Phi] Cos[\[Theta]] (\[Degree])", 
       "\[Theta] Sin[\[Theta]] (\[Degree])"}, {"\[Phi] (\[Degree])", 
       "\[Theta] (\[Degree])"}, {"\[Phi] (\[Degree])", 
       "\[Theta] Sin[\[Theta]] (\[Degree])"}}[[type]];
   Graphics[{
     {Lighter@Lighter@Lighter@Gray, 
      Polygon@Join[First@coors, Reverse@Last@coors]},
     {Lighter@Gray, Line[#]} & /@ coors,
     {Lighter@Gray, Line[#]} & /@ Transpose[coors]
     }, AspectRatio -> .8, ImageSize -> 400,
    PlotRange -> {{-185, 185}, {-95, 95}},
    LabelStyle -> 
     Directive[{Bold, Black, Medium, FontFamily -> "Helvetica"}],
    Frame -> True, FrameStyle -> Thick,
    Axes -> True, AxesStyle -> Thick,
    FrameTicks -> {{(Thread[{Round@coors[[1, All, 2]], ranY}][[
         1 ;; ;; 3]]), None}, {ranX, None}},
    FrameLabel -> lab
    ]
   ];


(* ::Subsubsection:: *)
(*Convert Grads*)


SyntaxInformation[ConvertGrads] = {"ArgumentsPattern" -> {_,_}};

(*convert gradient lists to txt output*)
ConvertGrads[gradi_, bv_] := Block[
  {depth, norm, gradu, grad, bval, bvalstr, gradstr, grad0str, list, 
   list0, part, listout, bvs, gr, name, nb},
  
  depth = ArrayDepth[gradi] /. 1 -> 3;
  norm = Map[Norm, gradi, {depth - 1}];
  
  gradu = DeleteDuplicates[Normalize /@ Flatten[gradi, depth - 2]];
  grad = Map[Normalize, gradi, {depth - 1}];
  bval = If[Length[bv] == 1, 
    bv[[1]]*(norm/Max[norm])^2, (norm/Max[norm])^2*bv];
  
  bvalstr = Flatten@Map[(
       bvs = 
        ToString[NumberForm[Round[Clip[#, {0, 35000}], 0.1], {7, 1}]];
       StringJoin[ConstantArray[" ", 10 - StringLength[bvs]]] <> bvs
       ) &, bval, {depth - 1}];
  gradstr = Flatten[Map[(
       gr = ToString[NumberForm[Round[#, 0.00001], {6, 5}]];
       If[StringTake[gr, 1] == "-", gr, " " <> gr]
       ) &, grad, {depth}], depth - 2];
  grad0str = Map[(
      gr = ToString[NumberForm[Round[#, 0.00001], {6, 5}]];
      If[StringTake[gr, 1] == "-", gr, " " <> gr]
      ) &, gradu, {2}];
  
  list = MapThread[
    StringJoin[Riffle[#1, "   "]] <> #2 &, {gradstr, bvalstr}];
  list0 = 
   Map[StringJoin[Riffle[#1, "   "]] <> "       0.1" &, grad0str];
  
  nb = ToString[Round[Max[Flatten[bval]]]];
  {list, list0, nb}
  ]


(* ::Subsubsection:: *)
(*FinalGrads*)


SyntaxInformation[FinalGrads] = {"ArgumentsPattern" -> {_,{_,_},{_,_}}};

FinalGrads[{listi_, list0_, nb_}, {inter_, int_}, {random_, ordr_}] :=
  Block[{part, listout, name, list},
  list = If[random, listi[[ordr]], listi];
  listout = DeleteDuplicates@Prepend[
     If[inter,
      part = Partition[list, int, int, 1, {}];
      Flatten@Riffle[part, list0[[;; Length[part]]]]
      , list
      ], " 0.00000    0.00000    1.00000       0.0"];
  name = " (" <> ToString[Length[listout]] <> ", " <> nb <> ")\n";
  name <> StringJoin[# <> "\n" & /@ listout]
  ]


(* ::Subsubsection:: *)
(*FindOrder*)


Options[FindOrder]={OrderSpan->"Auto"};

SyntaxInformation[FindOrder] = {"ArgumentsPattern" -> {_,_,OptionsPattern[]}};

FindOrder[grad_, bval_,OptionsPattern[]] := 
 Block[{local, minval, orderout, n, m, ord, val, temp, order,span},
  (*define the absvec and initial order*)
  temp = MakeAbsVec[grad, bval];
  orderout = order = Range[Length[temp]];
  
  span=OptionValue[OrderSpan];
  
  (*define initialization parameters*)
  local = Clip[If[span==="Auto",(*Round[Length[temp]/15],span]*)7,span],{5,Infinity}];
  (*local = 11;*)
  minval = Infinity;
  n = m = 0;
  
  (*brute force the solution, max 100000 itteration*)
  (*bo better solution after 25000 itteration probably already pritty good so stop*)
  While[n < 50000 && m < 150000,
   (*count itterations make new order and calculae load value*)
   m++;
   order = RandomSample[order];
   val = ValCalc[temp[[order]], local];
   
   (*if new optimum reset counter and update output*)
   If[val < minval,
    n = 0;
    minval = val;
    orderout = order;
    ,
    n++;
    ]];
  (*give optimal random order*)
  orderout
  ]

ValCalc =   Compile[{{vec, _Real, 2}, {local, _Integer, 0}}, Max[Mean[Transpose[Partition[vec, local, 1]]]]];
MakeAbsVec[grad_, bval_] := Partition[Flatten[Abs@(grad*Sqrt[bval])], 3];


(* ::Subsection::Closed:: *)
(*ConditionNumberCalc*)


SyntaxInformation[ConditionNumberCalc] = {"ArgumentsPattern" -> {_}};

ConditionNumberCalc[grad : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}] := Block[{gx, gy, gz},
  {gx, gy, gz} = Transpose[grad];
  (Max[#]/Min[#])&@SingularValueList[{gx^2, gy^2, gz^2, 2 gx gy, 2 gx gz, 2 gy gz}]
  ]

ConditionNumberCalc[mat_] := (Max[#]/Min[#])&@SingularValueList[mat];


(* ::Subsection::Closed:: *)
(*EnergyCalc*)


SyntaxInformation[EnergyCalc] = {"ArgumentsPattern" -> {_}};

EnergyCalc[grad_]:=Total@Flatten[Table[grad[[i]].grad[[j]],{i,2,Length[grad]},{j,1,i-1}]]


(* ::Subsection::Closed:: *)
(*OverPlus*)


SyntaxInformation[OverPlusCalc] = {"ArgumentsPattern" -> {_}};

OverPlusCalc[pts_]:=Module[{pt},pt=DeleteCases[pts//N,{0.,0.,0.}];Min[Norm/@((1/Max[pt])pt)]]


(* ::Subsection:: *)
(*Grad Rotation matrix*)


(* ::Subsubsection::Closed:: *)
(*GradRot*)


SyntaxInformation[GradRot] = {"ArgumentsPattern" -> {_,_}};

GradRot[Dfile_?StringQ,grad:{{_?NumberQ,_?NumberQ,_?NumberQ}..}]:=Module[{rot},
	rot=GradRotMat[Dfile];
	#*{1,1,-1} & /@ (grad.Inverse[rot])
]


(* ::Subsubsection::Closed:: *)
(*GradRotMat*)


SyntaxInformation[GradRotMat] = {"ArgumentsPattern" -> {_}};

GradRotMat[Dfile_String]:=
Module[{or},
	or = "ImageOrientation" /. Import[Dfile, "MetaInformation"];
	{or[[1;;3]],or[[4;;6]],Cross[or[[1;;3]],or[[4;;6]]]}
	]


(* ::Subsection::Closed:: *)
(*GradDirs*)


Options[GradDirs]={FullGrad->True};

SyntaxInformation[GradDirs] = {"ArgumentsPattern" -> {_,OptionsPattern[]}};

GradDirs[grad_?StringQ,OptionsPattern[]]:=
Module[{gradient},
	gradient=Switch[grad,
		(*Philips gradients*)
		"gr6",
		{{1,0,0},
		{0,1,0},
		{0,0,1},
		{0.7044,0.0882,0.7044},
		{0.7043,0.7044,0.0881},
		{0.0881,0.7044,0.7043}}
		,
		"gr6op",
		{{-0.3333,-0.6667,-0.6667},
		{-0.6667,-0.3333,0.6667},
		{0.6667,-0.6667,0.3333},
		{-0.7071,-0.7071,0},
		{0,-0.7071,0.7071},
		{0.7071,0,0.7071}}
		,
		"gr15",
		{{1,0,0},
		{0,1,0},
		{0,0,1},
		{-0.1789,-0.1113,-0.9776},
		{-0.0635,0.3767,-0.9242},
		{0.7108,0.0516,-0.7015},
		{0.6191,-0.4385,-0.6515},
		{0.2424,0.7843,-0.571},
		{-0.2589,-0.618,-0.7423},
		{-0.8169,0.1697,-0.5513},
		{-0.8438,0.5261,-0.106},
		{-0.2626,0.9548,-0.1389},
		{0.0001,0.9689,0.2476},
		{0.7453,0.6663,0.0242},
		{0.9726,0.2317,0.0209}}
		,
		"gr15op",
		{{-0.5000,-0.5000,-0.7071},
		{-0.5000,-0.5000,0.7071},
		{0.7071,-0.7071,0},
		{-0.1104,-0.7070,-0.6985},
		{0.2893,-0.6996,-0.6534},
		{0.6275,-0.3305,-0.7050},
		{0.6574,-0.2734,-0.7022},
		{-0.6725,-0.5421,-0.5038},
		{0.7038,-0.4911,0.5133},
		{-0.6930,-0.2531,0.6751},
		{-0.7066,-0.7071,0.02770},
		{-0.2821,-0.7070,0.6485},
		{0.2886,-0.7016,0.6514},
		{0.7058,-0.7063,0.05370},
		{0.7014,-0.2050,0.6827}}
		,
		"gr32",
		{{1,0,0},
		{0,1,0},
		{0,0,1},
		{-0.0425,-0.1146,-0.9925},
		{0.1748,-0.2005,-0.9639},
		{0.2323,-0.1626,-0.9589},
		{0.3674,0.0261,-0.9296},
		{0.1901,0.3745,-0.9075},
		{-0.1169,0.8334,-0.5402},
		{-0.2005,0.2527,-0.9466},
		{-0.4958,0.1345,-0.8579},
		{-0.0141,-0.6281,-0.7779},
		{-0.7445,-0.1477,-0.6511},
		{-0.7609,0.3204,-0.5643},
		{-0.1809,0.9247,-0.3351},
		{-0.6796,-0.4224,-0.5997},
		{0.7771,0.4707,-0.4178},
		{0.9242,-0.1035,-0.3676},
		{0.4685,-0.7674,-0.4378},
		{0.8817,-0.1893,-0.4322},
		{0.6904,0.7063,-0.1568},
		{0.239,0.7572,-0.6079},
		{-0.0579,0.9837,0.1703},
		{-0.5368,0.8361,-0.1135},
		{-0.9918,-0.1207,-0.0424},
		{-0.9968,0.0708,-0.0379},
		{-0.8724,0.478,-0.1014},
		{-0.2487,0.9335,0.2581},
		{0.1183,0.9919,-0.047},
		{0.3376,0.8415,0.4218},
		{0.5286,0.8409,0.1163},
		{0.9969,0.055,-0.057}}
		,
		"gr32op",
		{{-0.5,-0.5,-0.7071},
		{-0.5,-0.5,0.7071},
		{0.7071,-0.7071,0},
		{-0.6533,-0.2706,-0.7071},
		{-0.2087,-0.6756,-0.7071},
		{0.0197,-0.7068,-0.7071},
		{0.4212,-0.5679,-0.7071},
		{0.6899,-0.1549,-0.7071},
		{-0.6535,-0.2707,-0.7069},
		{-0.2929,-0.7071,-0.6436},
		{0.2945,-0.7064,-0.6436},
		{0.515,-0.4861,-0.7061},
		{0.7071,-0.2929,-0.6436},
		{-0.7071,-0.4725,-0.5261},
		{-0.4725,-0.7071,-0.5261},
		{0.5555,-0.6439,-0.5261},
		{0.7071,-0.4725,-0.5261},
		{-0.7071,-0.7071,-0.0002},
		{-0.7071,-0.4725,0.5261},
		{0.7071,-0.4725,0.5261},
		{0.4725,-0.7071,0.5261},
		{-0.7071,-0.7071,0.0078},
		{-0.6364,-0.4252,0.6436},
		{-0.706,-0.706,0.0547},
		{-0.2929,-0.7071,0.6436},
		{0.2929,-0.7071,0.6436},
		{0.7071,-0.7071,0.0078},
		{0.7071,-0.2929,0.6436},
		{-0.7063,-0.7063,0.0489},
		{0.0347,-0.7063,0.7071},
		{0.7071,-0.7071,0.0115},
		{0.7071,0,0.7071}}
		,
		(*Optimized gradients*)
		"opt6",
		{{0,0,1},
		{-0.8768,0.1769,0.4472},
		{0.8133,0.3722,0.4472},
		{0.6053,-0.6585,0.4472},
		{-0.1027,0.8885,0.4472},
		{-0.4392,-0.7792,0.4472}}
		,
		"opt6op",
		{{0.0011,-0.7071,0.7071},
		{0.7071,0.0061,0.7071},
		{0.7071,-0.707,0.0077},
		{-0.7071,0.0061,0.7071},
		{0.0011,0.7071,0.7071},
		{-0.707,-0.7071,0.0114}}
		,
		"opt10",
		{{1,0,0},
		{0,1,0},
		{0,0,1},
		{0.6861,-0.686,0.2422},
		{-0.5776,0.5773,0.5772},
		{-0.6861,-0.242,0.6861},
		{0.7071,0.7071,0.0002},
		{0.2414,0.6862,0.6862},
		{0,-0.7071,0.7071},
		{0.7071,0.0004,0.7071}}
		,
		"opt10op",
		{{0.6861,-0.1712,0.7071},
		{0.1727,0.6857,0.7071},
		{0.7071,0.4143,0.573},
		{0.7071,-0.7071,0.001},
		{-0.7071,-0.4068,0.5783},
		{-0.7071,0.1785,0.6842},
		{-0.409,0.7071,0.5769},
		{0.7071,0.7071,0},
		{0.412,-0.7071,0.5747},
		{-0.1768,-0.7071,0.6847}}
		,
		"opt15",
		{{1.,0.,0.},
		{0.,1.,0.},
		{0.,0.,1.},
		{-0.8088,-0.4165,0.4151},
		{0.5589,0.,0.8293},
		{0.8089,0.4153,0.4162},
		{0.2936,-0.6759,0.6759},
		{-0.2936,0.675,0.6768},
		{-0.5588,-0.0017,0.8293},
		{0.8088,-0.4163,0.4154},
		{0.2932,0.6761,0.6759},
		{-0.2931,-0.6768,0.6752},
		{0.5588,0.8293,0.0011},
		{-0.8089,0.415,0.4165},
		{-0.5589,0.8292,0.0013}}
		,
		"opt15op",
		{{-0.6165,-0.3464,0.7071},
		{0.5106,0.4891,0.7071},
		{0.1071,0.699,0.7071},
		{0.7071,-0.6583,0.2582},
		{0.7071,0.0995,0.7001},
		{0.7071,-0.325,0.628},
		{-0.7071,0.0576,0.7048},
		{-0.7071,-0.6978,0.1147},
		{-0.7071,0.4472,0.5477},
		{0.6327,0.7071,0.3159},
		{-0.3198,0.7071,0.6307},
		{-0.6856,0.7071,0.1732},
		{0.3546,-0.7071,0.6118},
		{-0.4934,-0.7071,0.5066},
		{-0.0724,-0.7071,0.7034}}
		,
		"opt24",
		{{1.,0,0},
		{0,1.,0},
		{0,0,1.},
		{0.8732,-0.4725,0.1199},
		{0.1275,-0.8588,0.4962},
		{0.4962,0.1275,0.8588},
		{-0.8627,-0.015,0.5055},
		{-0.4725,-0.1199,0.8732},
		{0.5566,0.8164,0.1539},
		{0.5809,-0.6408,0.502},
		{-0.383,-0.859,0.3397},
		{0.8164,-0.1539,0.5566},
		{0.859,0.3397,0.383},
		{0.502,0.5809,0.6408},
		{-0.1539,-0.5566,0.8164},
		{-0.6408,-0.502,0.5809},
		{-0.8588,-0.4962,0.1275},
		{-0.7917,0.4838,0.373},
		{0.5055,-0.8627,0.015},
		{-0.4838,0.373,0.7917},
		{0.015,0.5055,0.8627},
		{0.1199,0.8732,0.4725},
		{-0.373,0.7917,0.4838},
		{0.3397,-0.383,0.859}}
		,
		"opt24op",
		{{-0.7094,0.4985,0.4985},
		{-0.4985,0.7094,0.4985},
		{-0.6928,0.6928,0.2},
		{-0.4985,0.4985,0.7094},
		{-0.2,0.6928,0.6928},
		{-0.6928,0.2,0.6928},
		{0.2,0.6928,0.6928},
		{0.7094,0.4985,0.4985},
		{0.6928,0.6928,0.2},
		{0.4985,0.4985,0.7094},
		{0.4985,0.7094,0.4985},
		{0.6928,0.2,0.6928},
		{-0.7094,-0.4985,0.4985},
		{-0.6928,-0.2,0.6928},
		{-0.4985,-0.7094,0.4985},
		{-0.4985,-0.4985,0.7094},
		{-0.2,-0.6928,0.6928},
		{-0.6928,-0.6928,0.2},
		{0.4985,-0.4985,0.7094},
		{0.4985,-0.7094,0.4985},
		{0.6928,-0.6928,0.2},
		{0.6928,-0.2,0.6928},
		{0.7094,-0.4985,0.4985},
		{0.2,-0.6928,0.6928}}
		,
		"opt32",
		{{1.,0.,0.},
		{0.,1.,0.},
		{0.,0.,1.},
		{0.7045,0.,0.7097},
		{-0.9077,0.3751,0.1883},
		{0.3084,0.8999,0.3084},
		{0.5895,0.5432,0.5978},
		{0.4266,-0.9044,0.},
		{-0.4363,0.1203,0.8917},
		{-0.896,-0.1073,0.4309},
		{-0.744,0.279,0.6071},
		{0.9042,-0.2328,0.3581},
		{-0.6189,-0.2626,0.7403},
		{-0.7078,-0.5802,0.4029},
		{0.1679,0.7113,0.6825},
		{0.145,-0.6204,0.7708},
		{0.3538,-0.2368,0.9049},
		{-0.4406,-0.8807,0.174},
		{0.3821,0.2689,0.8841},
		{0.6209,-0.4731,0.625},
		{0.8814,0.2684,0.3887},
		{0.6847,0.7063,0.1798},
		{-0.9042,-0.4234,0.0559},
		{0.4036,-0.8077,0.4298},
		{-0.0617,0.4244,0.9034},
		{0.7631,-0.6243,0.1668},
		{-0.3494,-0.6855,0.6388},
		{-0.0326,-0.9037,0.4268},
		{-0.1807,0.8862,0.4265},
		{-0.2055,-0.3626,0.909},
		{-0.4043,0.5897,0.6991},
		{-0.6352,0.6971,0.3324}}
		,
		"opt32op",
		{{-0.681,-0.1904,0.7071},
		{0.5081,-0.4918,0.7071},
		{-0.6494,0.2797,0.7071},
		{-0.3606,-0.6083,0.7071},
		{0.5008,0.4992,0.7071},
		{0.6631,0.2455,0.7071},
		{-0.5655,-0.4245,0.7071},
		{-0.4842,0.5153,0.7071},
		{0.7071,0.6751,0.2105},
		{0.7071,-0.665,0.2403},
		{0.7071,0.4772,0.5218},
		{0.7071,-0.2428,0.6641},
		{0.7071,0.0012,0.7071},
		{0.7071,-0.488,0.5117},
		{-0.7071,-0.6542,0.2683},
		{-0.7071,-0.4929,0.5069},
		{-0.7071,0.6645,0.2417},
		{-0.7071,0.4826,0.5168},
		{-0.7071,0.0464,0.7056},
		{-0.7071,0.7071,0.0008},
		{0.5698,0.7071,0.4187},
		{-0.0933,0.7071,0.701},
		{-0.323,0.7071,0.629},
		{0.1392,0.7071,0.6933},
		{0.3661,0.7071,0.6049},
		{-0.5474,0.7071,0.4477},
		{-0.7065,-0.7071,0.0308},
		{0.0966,-0.7071,0.7005},
		{-0.512,-0.7071,0.4877},
		{0.3279,-0.7071,0.6265},
		{-0.1404,-0.7071,0.693},
		{0.5454,-0.7071,0.4501}}
		,
		"opt48",
		{{1.,0.,0.},
		{0.,1.,0.},
		{0.,0.,1.},
		{-0.5025,0.1031,0.8584},
		{0.674,-0.7165,0.1797},
		{0.3427,0.5997,0.7231},
		{-0.7936,0.5926,0.1376},
		{0.3885,-0.0601,0.9195},
		{0.116,-0.9154,0.3854},
		{-0.926,0.2704,0.2633},
		{0.1574,-0.3384,0.9277},
		{0.9175,-0.3773,0.1257},
		{0.9251,0.2712,0.2658},
		{-0.9162,-0.0819,0.3922},
		{-0.3321,0.7512,0.5704},
		{-0.6864,-0.1653,0.7082},
		{-0.465,0.4631,0.7545},
		{0.1689,-0.6875,0.7062},
		{0.1893,0.3176,0.9291},
		{-0.9323,-0.3458,0.1064},
		{-0.4875,-0.4884,0.7237},
		{-0.184,-0.7799,0.5982},
		{0.3645,-0.9254,0.1038},
		{-0.1966,0.2857,0.9379},
		{0.4924,0.7672,0.411},
		{0.6851,0.4864,0.5423},
		{-0.6783,0.5592,0.4766},
		{-0.5216,0.8244,0.2199},
		{-0.1795,0.9373,0.2988},
		{0.769,0.609,0.1941},
		{0.4582,-0.7524,0.4731},
		{-0.5054,-0.8629,0.},
		{0.4686,-0.4674,0.7496},
		{-0.1429,-0.5092,0.8487},
		{-0.251,-0.9343,0.2533},
		{0.2381,0.9406,0.2418},
		{-0.5149,-0.7365,0.4387},
		{-0.7739,-0.4296,0.4652},
		{0.7467,-0.472,0.4686},
		{-0.7674,0.2123,0.6049},
		{-0.3094,-0.1948,0.9308},
		{0.6804,-0.1765,0.7113},
		{-0.7407,-0.6506,0.1674},
		{0.5289,0.2571,0.8088},
		{0.7959,0.1401,0.589},
		{-0.06,0.5981,0.7992},
		{0.0854,0.8264,0.5566},
		{0.9193,-0.1233,0.3737}}
		,
		"opt48op",
		{{-0.6998,0.0438,0.7131},
		{-0.5762,0.5821,0.5738},
		{-0.7073,0.7069,0.0024},
		{-0.5612,0.4156,0.7158},
		{-0.6886,0.6918,0.2179},
		{-0.4075,0.5663,0.7165},
		{-0.6896,0.2491,0.6801},
		{-0.7176,0.5708,0.3993},
		{-0.41,0.7185,0.5621},
		{-0.2309,0.6891,0.687},
		{-0.5656,0.7192,0.4038},
		{-0.716,0.4219,0.5564},
		{0.6657,0.7192,0.1995},
		{0.6718,0.1732,0.7203},
		{0.3845,0.5695,0.7266},
		{0.5428,0.5926,0.5953},
		{0.198,0.6928,0.6936},
		{0.724,0.5966,0.3465},
		{0.5595,0.7138,0.4214},
		{0.5612,0.414,0.7168},
		{-0.0193,0.7068,0.7072},
		{0.693,0.5025,0.5171},
		{0.3832,0.7258,0.5714},
		{0.7232,0.3293,0.6072},
		{-0.5059,-0.5017,0.7018},
		{-0.5899,-0.5881,0.5534},
		{-0.7186,-0.1697,0.6745},
		{-0.724,-0.5653,0.3954},
		{-0.1766,-0.7183,0.6731},
		{-0.6871,-0.6924,0.2204},
		{-0.3356,-0.6045,0.7226},
		{-0.564,-0.7235,0.3982},
		{-0.6081,-0.331,0.7217},
		{-0.7167,-0.4113,0.5633},
		{-0.4123,-0.7154,0.5642},
		{-0.7141,-0.7,0.0095},
		{0.4092,-0.5662,0.7156},
		{0.0362,-0.6994,0.7139},
		{0.7179,-0.4151,0.5589},
		{0.2424,-0.6928,0.6793},
		{0.6917,-0.6899,0.2138},
		{0.5706,-0.7187,0.3975},
		{0.6829,-0.2455,0.6881},
		{0.5807,-0.5799,0.5716},
		{0.7131,-0.0394,0.7001},
		{0.561,-0.4175,0.715},
		{0.7192,-0.5681,0.4002},
		{0.4199,-0.7172,0.5563}}
		,
		"opt60",
		{{1.,0.,0.},
		{0.,1.,0.},
		{0.,0.,1.},
		{-0.9479,0.1784,0.264},
		{-0.2951,0.4896,0.8205},
		{-0.5432,0.5887,0.5986},
		{-0.0114,0.3426,0.9394},
		{0.9451,-0.2844,0.1611},
		{-0.6878,-0.7254,0.0289},
		{-0.5959,-0.0153,0.8029},
		{0.7874,0.5632,0.2507},
		{0.8152,0.2997,0.4956},
		{-0.7824,-0.2996,0.5459},
		{-0.9428,-0.1789,0.2812},
		{0.6047,-0.5814,0.5443},
		{0.2166,0.9404,0.2621},
		{-0.1561,0.9455,0.2856},
		{-0.7936,0.3481,0.4991},
		{0.587,0.3035,0.7505},
		{0.768,-0.0088,0.6404},
		{-0.3706,0.9286,0.0177},
		{-0.5827,0.2975,0.7563},
		{0.0444,0.8489,0.5267},
		{-0.2907,-0.174,0.9409},
		{0.7816,-0.5751,0.2416},
		{0.5229,-0.0026,0.8524},
		{0.0232,0.6313,0.7752},
		{0.9457,0.2661,0.1868},
		{-0.2761,-0.5166,0.8105},
		{0.0308,-0.6467,0.7622},
		{-0.4665,0.822,0.3266},
		{-0.1879,-0.9447,0.2688},
		{0.3432,0.7486,0.5673},
		{-0.7628,-0.5565,0.3293},
		{0.0208,-0.8592,0.5112},
		{-0.6623,0.749,0.0207},
		{-0.5441,-0.3238,0.774},
		{0.3246,0.5024,0.8014},
		{0.28,-0.1917,0.9406},
		{-0.8828,0.4602,0.094},
		{0.5198,-0.8049,0.2864},
		{0.5343,0.7967,0.2824},
		{0.3271,-0.7548,0.5685},
		{-0.8977,-0.4345,0.0726},
		{0.6127,0.5703,0.5471},
		{-0.5578,-0.5899,0.5839},
		{-0.3938,-0.9192,0.0105},
		{-0.5075,-0.8031,0.3122},
		{-0.3117,0.1623,0.9362},
		{0.2049,-0.9426,0.2636},
		{-0.0004,-0.3554,0.9347},
		{0.3281,-0.5011,0.8008},
		{-0.7226,0.6127,0.3201},
		{0.2798,0.194,0.9403},
		{-0.2615,0.7589,0.5964},
		{0.8155,-0.3202,0.4821},
		{-0.8268,0.0257,0.5619},
		{0.9338,-0.0156,0.3573},
		{0.5914,-0.3093,0.7447},
		{-0.2816,-0.7677,0.5756}}
		,
		"opt60op",
		{{-0.6078,0.7325,0.3066},
		{-0.7325,0.6078,0.3066},
		{-0.6204,0.6204,0.4799},
		{-0.1054,0.7032,0.7032},
		{-0.6204,0.4799,0.6204},
		{-0.4799,0.6204,0.6204},
		{-0.7325,0.3066,0.6078},
		{-0.4809,0.7332,0.4809},
		{-0.7332,0.4809,0.4809},
		{-0.3066,0.6078,0.7325},
		{-0.3066,0.7325,0.6078},
		{-0.6078,0.3066,0.7325},
		{-0.7032,0.7032,0.1054},
		{-0.7032,0.1054,0.7032},
		{-0.4809,0.4809,0.7332},
		{0.1054,0.7032,0.7032},
		{0.7032,0.7032,0.1054},
		{0.7325,0.3066,0.6078},
		{0.7332,0.4809,0.4809},
		{0.6204,0.6204,0.4799},
		{0.6078,0.7325,0.3066},
		{0.6204,0.4799,0.6204},
		{0.6078,0.3066,0.7325},
		{0.7325,0.6078,0.3066},
		{0.4809,0.7332,0.4809},
		{0.7032,0.1054,0.7032},
		{0.4809,0.4809,0.7332},
		{0.3066,0.7325,0.6078},
		{0.3066,0.6078,0.7325},
		{0.4799,0.6204,0.6204},
		{-0.3066,-0.7325,0.6078},
		{-0.4809,-0.4809,0.7332},
		{-0.7032,-0.7032,0.1054},
		{-0.7332,-0.4809,0.4809},
		{-0.6078,-0.7325,0.3066},
		{-0.6204,-0.4799,0.6204},
		{-0.4799,-0.6204,0.6204},
		{-0.4809,-0.7332,0.4809},
		{-0.7032,-0.1054,0.7032},
		{-0.6204,-0.6204,0.4799},
		{-0.1054,-0.7032,0.7032},
		{-0.3066,-0.6078,0.7325},
		{-0.6078,-0.3066,0.7325},
		{-0.7325,-0.6078,0.3066},
		{-0.7325,-0.3066,0.6078},
		{0.3066,-0.6078,0.7325},
		{0.1054,-0.7032,0.7032},
		{0.7032,-0.7032,0.1054},
		{0.6078,-0.3066,0.7325},
		{0.3066,-0.7325,0.6078},
		{0.4809,-0.7332,0.4809},
		{0.7325,-0.6078,0.3066},
		{0.6078,-0.7325,0.3066},
		{0.4799,-0.6204,0.6204},
		{0.4809,-0.4809,0.7332},
		{0.7032,-0.1054,0.7032},
		{0.7325,-0.3066,0.6078},
		{0.6204,-0.4799,0.6204},
		{0.6204,-0.6204,0.4799},
		{0.7332,-0.4809,0.4809}}
		,
		_,
		Return[Message[GradDirs::grad,grad]]
		];
	If[OptionValue[FullGrad],Prepend[gradient,{0,0,0}],gradient]
	]



(* ::Subsection:: *)
(*Bmatrix Functions*)


(* ::Subsubsection::Closed:: *)
(*B-matrix*)


Options[Bmatrix]={Method->"DTI"};

SyntaxInformation[Bmatrix] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

Bmatrix[{bvec_,grad_}, OptionsPattern[]]:=Switch[
	OptionValue[Method],
	"DTI",
	bvec*GradVecConv[grad,1],
	"DKI",
	If[(Length[grad]==Length[bvec])||(Length[bvec]==0),
		MapThread[Join[#1,#2]&,{bvec*GradVecConv[grad,1], bvec^2*GradVecConv[grad,2]}],
		MapThread[Join[#1,#2]&,{
			Flatten[#*GradVecConv[grad,1]& /@ bvec, 1], 
			Flatten[#^2*GradVecConv[grad,2]& /@ bvec, 1]
			}]
	]
]

Bmatrix[bvec_,grad_,OptionsPattern[]]:=Switch[
	OptionValue[Method],
	"DTI",
	Append[#,1]&/@(-bvec*GradVecConv[grad,1]),
	"DKI",
	If[(Length[grad]==Length[bvec])||(Length[bvec]==0),
		Append[#,1]&/@(MapThread[Join[#1,#2]&,{-bvec*GradVecConv[grad,1], (1/6)*bvec^2*GradVecConv[grad,2]}]),
		Append[#,1]&/@(MapThread[Join[#1,#2]&,{
			Flatten[-#*GradVecConv[grad,1]& /@ bvec, 1], 
			Flatten[(1/6)*#^2*GradVecConv[grad,2]& /@ bvec, 1]
			}])
	]
]

GradVecConv[grad_,type_]:=Block[{gx,gy,gz},
		{gx,gy,gz}=Transpose[grad];
		Transpose@Switch[type,
			1,{gx^2,gy^2,gz^2,2 gx gy,2 gx gz,2 gy gz},
			2,{gx^4, gy^4, gz^4, 4 gx^3 gy, 4 gx^3 gz, 6 gx^2 gy^2, 12 gx^2 gy gz, 6 gx^2 gz^2, 4 gx gy^3, 12 gx gy^2 gz, 12 gx gy gz^2, 4 gx gz^3, 4 gy^3 gz, 6 gy^2 gz^2, 4 gy gz^3}
	]
]



(* ::Subsubsection::Closed:: *)
(*BmatrixCalc*)


Options[BmatrixCalc] = {
   UseGrad -> {1, 1, {1, 1}, 1, 1},
   OutputType -> "Matrix",
   Method -> "Numerical",
   StepSizeI -> 0.05,
   UnitMulti -> 10^-3,
   PhaseEncoding -> "A",
   FlipAxes -> {{1, 1, 1}, {1, 1, 1}},
   SwitchAxes -> {{1, 2, 3}, {1, 2, 3}}};

SyntaxInformation[BmatrixCalc] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

BmatrixCalc[folder_, grads_, opts : OptionsPattern[]] := Module[
  {seq, Gt, hw, te, bmat, t},
  seq = ImportGradObj[folder];
  bmat = Map[(
      {Gt, hw, te} = 
       GradSeq[seq, t, #, UseGrad -> OptionValue[UseGrad], 
        UnitMulti -> OptionValue[UnitMulti], 
        FilterRules[{opts}, Options[GradSeq]]];
      Chop[GradBmatrix[
        Gt, hw, te, t, Method -> OptionValue[Method], StepSizeI -> OptionValue[StepSizeI], 
        FilterRules[{opts}, Options[GradBmatrix]]
        ]]
      ) &, grads];
  Switch[OptionValue[OutputType], "Matrix", bmat, "Gradient", 
   BmatrixInv[#] & /@ bmat]]


(* ::Subsubsection::Closed:: *)
(*BmatrixInv*)


SyntaxInformation[BmatrixInv] = {"ArgumentsPattern" -> {_, _.}};

BmatrixInv[bm_, bvi___] := Module[{bv, sigb, sign, gr},
  If[ArrayDepth[bm] == 2,
   
   Transpose[BmatrixInv[#,bvi]& /@ bm],
   
   bv = Total[bm[[1 ;; 3]]];
   sigb=Sign[bv];
   bv = If[bvi === Null, sigb bv, bvi];
   (*sign = Switch[
   	Sign[Sign[sigb bm[[4 ;; 6]]] + .00001], 
    {-1, -1, 1}, {-1, 1, 1},
    {-1, 1, -1}, {1, -1, 1}, 
    {1, -1, -1}, {-1, -1, 1},
    _, {1, 1, 1}
    ];*)
    
    sign = If[Sign[bm[[3]]] == 0,If[Sign[bm[[2]]] == 0,If[Sign[bm[[1]]] == 0,
     {1, 1, 1},Sign2[bm[[{1, 3, 4}]] sigb]],Sign2[bm[[{4, 2, 5}]] sigb]],Sign2[bm[[{5, 6, 3}]] sigb]];
    
   gr = sign*Sqrt[sigb bm[[1 ;; 3]]/(bv /. 0. -> Infinity)];
   
   {bv, gr}]
   ]

Sign2 = Sign[Sign[# + .00001]] &;



(* ::Subsubsection::Closed:: *)
(*BmatrixConv*)


SyntaxInformation[BmatrixConv] = {"ArgumentsPattern" -> {_}};

BmatrixConv[bmat_] := Module[{},
  If[Length[bmat[[1]]] == 6,
   Append[-#, 1] & /@ bmat,
   -bmat[[All, 1 ;; 6]]
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*BmatrixRot*)


SyntaxInformation[BmatrixRot] = {"ArgumentsPattern" -> {_, _}};

BmatrixRot[bmat_, rotmat_] := 
 Module[{sc1 = {1, 1, 1, 0.5, 0.5, 0.5}, sc2 = {1, 1, 1, 2, 2, 2}},
  If[Dimensions[bmat] == {6},
   sc2 TensVec[Transpose[rotmat].TensMat[sc1 bmat].rotmat],
   sc2 TensVec[Transpose[rotmat].TensMat[sc1 #].rotmat] & /@ bmat
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*BmatrixToggle*)


SyntaxInformation[BmatrixToggle] = {"ArgumentsPattern" -> {_, _, _}};

BmatrixToggle[bmat_, axes_, flip_] := 
 Block[{tmpa, tmpf, bmati, bmatn, outp, rule},
  (*error checking*)
  If[Sort[axes] != {"x", "y", "z"},
   Return[Message[BmatrixToggle::axes, axes]],
   If[!MemberQ[{{1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}}, flip],
    Return[Message[BmatrixToggle::flip, flip]],
    
    (*make bmat vecs of input*)
    tmpa = axes /. {x_, y_, z_} -> {x*x, y*y, z*z, x y, x z, y z};
    tmpf = flip /. {x_, y_, z_} -> {x*x, y*y, z*z, x y, x z, y z};
    
    (*check shape of bmat*)
    bmati = If[Length[bmat[[1]]] == 7, -bmat[[All, 1 ;; 6]], bmat];
    
    (*Toggle bmat*)
    bmatn = (
        rule = 
         Thread[{("x")^2, ("y")^2, ("z")^2, "x" "y", "x" "z", 
            "y" "z"} -> #];
        tmpf (tmpa /. rule)
        ) & /@ bmati;
    (*output bmat*)
    outp = If[Length[bmat[[1]]] == 7, Append[-#, 1] & /@ bmatn, bmatn]
    ]
   ]
  ]


(* ::Subsection::Closed:: *)
(*Import Grad object*)


SyntaxInformation[ImportGradObj] = {"ArgumentsPattern" -> {_}};

ImportGradObj[folder_] := 
 Module[{files, imp, obj}, 
 	files = FileNames["GR*.acq", folder];
 	(
 		imp = Import[#, "Lines"];
 		obj = (obj = StringSplit[StringDrop[StringTrim[#], -1], " = "];
 		obj[[1]] -> ToExpression[obj[[2]]] // N) & /@ imp[[2 ;;]];
 		StringReplace[imp[[1]], "$ modify_object " -> ""] -> obj
 		) & /@files
 	]

ImportGradObj[{base_, xbase_}] := 
 Module[{objectNames, objects, name, vals, props},
  objectNames = {"\"GR`blip\"", "\"GR`d_echo\"",
    "\"GR`diff[0]\"", "\"GR`diff[1]\"", "\"GR`diff[2]\"",
    "\"GR`diff_2nd[0]\"", "\"GR`diff_2nd[1]\"", "\"GR`diff_2nd[2]\"",
    "\"GR`r_diff[0]\"", "\"GR`r_diff[1]\"", "\"GR`r_diff[2]\"",
    "\"GR`diff_crush[0]\"", "\"GR`diff_crush[1]\"",
    "\"GR`m[0]\"", "\"GR`mc[0]\"", "\"GR`md\"",
    "\"GR`mf_base[0]\"", "\"GR`mf_base[1]\"", "\"GR`pf[0]\"", 
    "\"GR`pf[1]\"", "\"GR`sf_base[0]\"", "\"GR`sf_base[1]\"",
    "\"GR`py\"",
    "\"GR`r_echo\"", "\"GR`s_echo\"",
    "\"GR`r_ex\"", "\"GR`s_ex\"",
    "\"GR`s_ex1\"", "\"GR`s_ex2\"", "\"GR`s_diff_echo\"", 
    "\"GR`TM_crush\""
    };
  
  objects = Flatten[Partition[SplitBy[Import[#, "Lines"], (StringTake[#, 1] === "$" &)],2] & /@ {base, xbase}, 1];
  objects = Sort@DeleteCases[(
  	name = Last[StringSplit[First@First@#]];
  	If[MemberQ[objectNames, name],
  		vals = #[[2]];
  		vals = (
  			props = StringSplit[StringDrop[StringTrim[#], -1], " = "];
  			props[[1]] -> ToExpression[props[[2]]] // N
  		) & /@ vals;
  		name -> vals
  	]) & /@ objects, Null]
]


(* ::Subsection::Closed:: *)
(*Gradient sequence*)


Options[GradSeq] = {UseGrad -> {0, 1, {1, 0}, 1}, FlipGrad -> False, UnitMulti -> 1, 
	PhaseEncoding -> "A", FlipAxes->{{1,1,1},{1,1,1}},SwitchAxes->{{1,2,3},{1,2,3}}
	};

SyntaxInformation[GradSeq] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

GradSeq[pars_, t_, grad : {_, _, _}, OptionsPattern[]] := Module[{
   seq, rule, G, slope, leng, start, repn, repf, repi, dur, ori, flip,
    grepi1, grepi2, gr180, grdiff, grex, unit, func, name, str, seqt, 
   startRep, orVec, t180, te, hw, i, usegrad, grflow, g1, g2, g3, t1, 
   t2, t3, t4, t901, t902, t1801, t1802, flipg, AP},
  
  unit = OptionValue[UnitMulti];
  Clear[t];
  
  AP = If[OptionValue[PhaseEncoding] == "A"||OptionValue[PhaseEncoding] == "R", 1, -1];
  
  flipg = 
   If[MemberQ[pars[[All, 1]], "\"GR`s_diff_echo\""], 
   	{"\"GR`s_echo\"", "\"GR`d_echo\"", "\"GR`r_echo\""}, 
   	{"\"GR`s_echo\"", "\"GR`d_echo\"", "\"GR`r_echo\"", "\"GR`s_ex1\"", "\"GR`s_ex2\"", "\"GR`diff_crush[1]\"", "\"GR`diff_crush[0]\""}];
  
  grepi1 = {"\"GR`mc[0]\"", "\"GR`md\""};
  grepi2 = {"\"GR`blip\"", "\"GR`m[0]\"", "\"GR`py\""};
  
  gr180 = {"\"GR`d_echo\"", "\"GR`r_echo\"", "\"GR`s_echo\"", 
    "\"GR`diff_crush[0]\"", "\"GR`diff_crush[1]\"", 
    "\"GR`s_diff_echo\"", "\"GR`s_ex1\"", "\"GR`s_ex2\"", 
    "\"GR`TM_crush\""};
  grex = {"\"GR`r_ex\"", "\"GR`s_ex\""};
  
  grflow = {
  	"\"GR`sf_base[0]\"", "\"GR`sf_base[1]\"", 
    "\"GR`mf_base[0]\"", "\"GR`mf_base[1]\"",
     "\"GR`pf[0]\"", "\"GR`pf[1]\""
     };
  
  grdiff = {"\"GR`diff[2]\"", "\"GR`diff_2nd[2]\""};
    
  usegrad = OptionValue[UseGrad];
  (*{grex, gr180, {grepi1, grepi2}, grdiff, grflow}*)
  If[Length[usegrad] == 4, AppendTo[usegrad, 0]];
  seqt = Transpose[DeleteCases[(
         
         name = #[[1]];
         
         If[
          MemberQ[Flatten[
            usegrad {grex, gr180, {grepi1, grepi2}, grdiff, grflow}], 
           name], rule = #[[2]];
          
          G = (
          	str = If[("gr_str" /. rule) != 0,"gr_str" /. rule,("gr_str_step"*"gr_str_factor_max") /. rule];
          	If[("gr_lenc" /. rule) >= 0.,str,((str/"gr_slope")*("gr_slope" + "gr_lenc")) /. rule]
          	)unit;
          
          G = If[MemberQ[flipg, name] && OptionValue[FlipGrad], -G, G];
          G = If[MemberQ[Join[grepi1, grepi2], name], AP*G, G];
          
          (*{slope,leng,dur,start,repn,repf,repi,ori}={
          If[("gr_lenc"/.rule)\[GreaterEqual]0.,
          "gr_slope"/.rule,("gr_slope"+"gr_lenc")/.rule] unit,
          Abs[("gr_lenc"/.rule)] unit,
          ("gr_dur"/.rule) unit,
          (("gr_time"-"gr_ref")/.rule) unit,
          "gr_repetitions",
          "gr_rep_alt_factor",
          ("gr_interval"/.rule) unit,
          "gr_ori"
          }/.rule;*)
          
          slope = If[("gr_lenc" /. rule) >= 0.,"gr_slope" /. rule,("gr_slope" + "gr_lenc") /. rule] unit;
          leng = Abs[("gr_lenc" /. rule)] unit;
          dur = ("gr_dur" /. rule) unit;
          start = (("gr_time" - "gr_ref") /. rule) unit;
          repn = "gr_repetitions" /. rule;
          repf = "gr_rep_alt_factor" /. rule;
          repi = ("gr_interval" /. rule) unit;
          ori = "gr_ori" /. rule;
          
          orVec = RotateLeft[If[MemberQ[grdiff, name], grad // N, {0, 0, 1} // N],Abs[ori - 2]];
          
          func = Transpose[
            Table[flip = repf^(i + 1);
             startRep = start + (i - 1) (dur + repi);
             t1 = startRep;
             t2 = startRep + slope;
             t3 = startRep + slope + leng;
             t4 = startRep + 2 slope + leng;
             g1 = (G/slope) (t - t1);
             g2 = G;
             g3 = G - ((G/slope) (t - t3));
             
             orVec*Piecewise[{{flip g1, t1 <= t <= t2}, {flip g2, 
                 t2 <= t <= t3}, {flip g3, t3 <= t <= t4}}]
             , {i, 1, repn}]
            ];
          
          (*func = If[MemberQ[grdiff, name],
            {-1, 1, -1} func[[{2, 1, 3}]],
            Total[(#*func)] & /@ ({-1, 1, -1} slmat[[{2, 1, 3}]])
            ]*)
            
          func = If[MemberQ[grdiff, name],
            (OptionValue[FlipAxes][[1]])*func[[OptionValue[SwitchAxes][[1]]]],
            Total[#]&/@(OptionValue[FlipAxes][[2]]*func[[OptionValue[SwitchAxes][[2]]]])
            ]
          ]
         ) & /@ pars, Null]
     ];
     
  seq = Total[Flatten[#]] & /@ seqt;
  
  te = ("gr_time" /. Take[pars, Position[pars[[All, 1]], "\"GR`m[0]\""][[1]]][[1, 2]]) unit;
  
  If[!StringQ["\"GR`s_echo\"" /. pars] && StringQ["\"GR`s_diff_echo\"" /. pars],
   t180 = ("gr_time" /.Take[pars, Position[pars[[All, 1]], "\"GR`s_echo\""][[1]]][[1,2]]) unit;
   hw = Piecewise[{{1, 0 <= t <= t180}, {-1, t180 <= t (*<= te*)}}];];
  
  If[!StringQ["\"GR`s_echo\"" /. pars] && !StringQ["\"GR`s_diff_echo\"" /. pars],
   t1801 = ("gr_time" /. Take[pars, Position[pars[[All, 1]], "\"GR`s_echo\""][[1]]][[1,2]]) unit;
   t1802 = ("gr_time" /. Take[pars, Position[pars[[All, 1]], "\"GR`s_diff_echo\""][[1]]][[1, 2]]) unit;
   hw = Piecewise[{{1, 0 <= t <= t1801}, {-1, t1801 <= t <= t1802}, {1, t1802 <= t (*<= te*)}}];
   ];
  
  If[!StringQ["\"GR`s_ex1\"" /. pars] && !StringQ["\"GR`s_ex2\"" /. pars],
   t901 = ("gr_time" /. Take[pars, Position[pars[[All, 1]], "\"GR`s_ex1\""][[1]]][[1, 2]]) unit;
   t902 = ("gr_time" /. Take[pars, Position[pars[[All, 1]], "\"GR`s_ex2\""][[1]]][[1, 2]]) unit;
   hw = Piecewise[{{1, 0 <= t <= t901}, {-1, t902 <= t (*<= te*)}}];
   ];
  
  {PiecewiseExpand/@seq, PiecewiseExpand@hw, te}
  ]


(* ::Subsection::Closed:: *)
(*GradBmatrix*)


Options[GradBmatrix] = {OutputPlot -> False, Method -> "Analytical", StepSizeI -> 0.025};

SyntaxInformation[GradBmatrix] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

GradBmatrix[Gti_, hw_, te_, t_, OptionsPattern[]] := Module[{Ft, Ft2, Ft2i, s = 267.522 10^6, plot, bmat, Gtfn, Gt},
  (*Gt = {1, -1, -1} Gti[[{2, 1, 3}]];*)
  Gt=Gti;  
  Switch[OptionValue[Method],
   "Analytical",
   Ft = Chop[Integrate[s hw # 10^-3 // N, t], 10^-5] & /@ Gt;
   Ft2 = N[PiecewiseExpand[#]] & /@ {Ft[[1]] Ft[[1]], Ft[[2]] Ft[[2]], Ft[[3]] Ft[[3]], 2 Ft[[1]] Ft[[2]], 2 Ft[[1]] Ft[[3]], 2 Ft[[2]] Ft[[3]]};
   Ft2i = Map[Integrate[#, t] &, Ft2];
   bmat = ((Ft2i /. t -> te) - (Ft2i /. t -> 0));,
   "Numerical",
   Gtfn = Transpose[Table[s hw Gt 10^-3, {t, 0, te, OptionValue[StepSizeI]/1000}]];
   Ft = Integrate[ListInterpolation[#, {0, te}][t], t] & /@ Gtfn;
   Ft2 = {Ft[[1]] Ft[[1]], Ft[[2]] Ft[[2]], Ft[[3]] Ft[[3]], 2 Ft[[1]] Ft[[2]], 2 Ft[[1]] Ft[[3]], 2 Ft[[2]] Ft[[3]]};   
   bmat = Quiet[(NIntegrate[#, {t, 0, te}]) & /@ Ft2];
   ];
   
  If[OptionValue[OutputPlot],
   plot = GraphicsGrid[
       Partition[
        Plot[#1, {t, 0, te}, PlotRange -> {{-.1 te, 1.1 te}, Full}, 
           PlotPoints -> 500, Exclusions -> None,
           PlotRange -> Full, AspectRatio -> .2, 
           PlotStyle -> Directive[{Black, Thick}]] & /@ #, 3], 
       ImageSize -> 1000
       ] & /@ {Gt, Ft, Ft2};
       
   {bmat, plot},
   bmat]
   ]


(* ::Subsection:: *)
(*GetSliceNormal*)


(* ::Subsubsection::Closed:: *)
(*GetSliceNormal*)


SyntaxInformation[GetSliceNormal] = {"ArgumentsPattern" -> {_,_.}};

GetSliceNormal[folder_String,part_Integer] := Module[{or,files,grads,norm,gradRotmat},
	files=FileNames["*.dcm",folder];
	grads=-Round[If[part==1,GradRead[files[[1]]],GradRead/@files[[;;part]]],.00001];
	or = "ImageOrientation" /. Import[files[[1]], "MetaInformation"];
	gradRotmat=Transpose[{or[[1 ;; 3]], or[[4 ;; 6]],Cross[or[[1 ;; 3]], or[[4 ;; 6]]]}];
	norm={0, 0, 1}.gradRotmat;
	{norm,grads,gradRotmat}
   ]


(* ::Subsubsection::Closed:: *)
(*GetSliceNormalDir*)


SyntaxInformation[GetSliceNormalDir] = {"ArgumentsPattern" -> {_}};

GetSliceNormalDir[Dfile_String] := 
 Module[{meta, directions, slice, groups, orientation,grads,norm,gradRotmat},
  meta = Import[Dfile, "MetaInformation"];
  directions = "(2005,1415)" /. meta;
  slice = "(2001,1018)" /. meta;
  groups = If[("FrameCount"/slice /. meta) == directions, 
  	("(5200,9230)" /. meta)[[1 ;; directions]], 
  	("(5200,9230)" /. meta)[[1 ;; ("FrameCount"/"(2001,1018)") /. meta]]
  	];
  orientation = "ImageOrientation" /. ("(0020,9116)" /. groups[[1]]);
  
gradRotmat = {v1 = {-1, 1, -1} orientation[[{5, 4, 6}]],
 v2 = {1, -1, 1} orientation[[{2, 1, 3}]],
 -Cross[v1, v2]
 } // Transpose;
  
  (*gradRotmat = Transpose[{orientation[[1 ;; 3]], orientation[[4 ;; 6]], Cross[orientation[[1 ;; 3]], orientation[[4 ;; 6]]]}];*)
  grads = If[("(0018,9075)" /. #) == "NONE" || ("(0018,9075)" /. #) == "ISOTROPIC", 
   	{0, 0, 0},("(0018,9089)" /. ("(0018,9076)" /. #))
   	] & /@ (("(0018,9117)" /. #) & /@ groups);
  norm={0., 0, 1}.gradRotmat;
  {norm,grads,gradRotmat}
  ]


(* ::Subsubsection::Closed:: *)
(*CalculateMoments*)


SyntaxInformation[CalculateMoments] = {"ArgumentsPattern" -> {_, _}};

CalculateMoments[{Gt_, hw_, te_}, t_] := Module[{fun, M0, M1, M2, M3, vals},
	fun = N@PiecewiseExpand[#*hw] & /@ Gt;
	
	M0 = hw Integrate[PiecewiseExpand[# ]    , t, Assumptions -> t >= 0 && t <= te && t \[Element] Reals, GenerateConditions -> False] & /@ fun;
	M1 = hw Integrate[PiecewiseExpand[# t]   , t, Assumptions -> t >= 0 && t <= te && t \[Element] Reals, GenerateConditions -> False] & /@ fun;
	M2 = hw Integrate[PiecewiseExpand[# t^2 ], t, Assumptions -> t >= 0 && t <= te && t \[Element] Reals, GenerateConditions -> False] & /@ fun;
	M3 = hw Integrate[PiecewiseExpand[# t^3 ], t, Assumptions -> t >= 0 && t <= te && t \[Element] Reals, GenerateConditions -> False] & /@ fun;
	
	vals = {M0, M1, M2, M3} /. t -> te;
	
	{{PiecewiseExpand/@Gt, M0, M1, M2, M3}, vals}
  ]


(* ::Subsubsection::Closed:: *)
(*UniqueBvalPosition*)


SyntaxInformation[UniqueBvalPosition] = {"ArgumentsPattern" -> {_, _.}};

UniqueBvalPosition[bval_, num_: 1] := Block[{bvalU, pos},
  bvalU = Sort[DeleteDuplicates[bval]];
  pos = Flatten[Position[bval, #]] & /@ bvalU;
  Transpose[Select[Transpose[{bvalU, pos}], Length[#[[2]]] >= num &]]
  ]


(* ::Subsubsection::Closed:: *)
(*UniqueBvalPosition*)


(*Options[GetGradientScanOrder]={RoundBvalue->False}*)

SyntaxInformation[GetGradientScanOrder] = {"ArgumentsPattern" -> {_, _, _.(*,OptionsPattern[]*)}};

GetGradientScanOrder[grad_?ListQ, bval_?ListQ(*,opts:OptionsPattern[]*)] := Module[{file},
  file = FileSelect["FileOpen", {"*.txt"}, 
    WindowTitle -> "Select gradient input *.txt"];
  GetGradientScanOrder[file, grad, bval,opts]
  ]

GetGradientScanOrder[file_?StringQ, grd_?ListQ, bval_?ListQ,OptionsPattern[]] := 
 Module[{input, bvali, gradi, order1, p1, p2, int, order2,grad},
  input = Import[file, "Lines"];
  input = (ToExpression@
       DeleteCases[StringSplit[StringTrim[#]] // StringTrim, ""]) & /@
     input[[2 ;;]];
     input = input /. {_, _, _, 0.} -> {0., 0., 0., 0.};
(*  bvali = Append[input[[All, 4]], 0.1];
  gradi = Append[input[[All, 1 ;; 3]], {0., 0., 0.}];*)
  
  bvali=input[[All,4]];
  gradi=Normalize /@input[[All,1;;3]];
  
  (*bvali=If[OptionValue[RoundBvalue],Round[bvali],bvali];*)
  
  grad = Sign[Sign[grd[[All, 3]]] + 0.00001] grd;
  
  order1 = Flatten[MapThread[     (
       p1 = Position[Round[bval, .1], #1]; 
       p2 = Position[Round[grad, .1], #2]; 
       int = Intersection[p1, p2]
       ) &, {Round[bvali, .1], Round[gradi, .1]}
     ]
    ];
  order2 = Flatten[MapThread[(
  		p1 = Position[Round[bvali, .1], #1];
  		p2 = Position[Round[gradi, .1], #2]; 
  		int = Intersection[p1, p2]
       ) &, {Round[bval, .1], Round[grad, .1]}
     ]
    ];
  {order1, order2}
  ]


(* ::Subsection::Closed:: *)
(*BvecCreate*)


SyntaxInformation[BvecCreate] = {"ArgumentsPattern" -> {_,_}};

BvecCreate[bvalue_, gr_] := 
 Join[{0}, ConstantArray[bvalue, {Length[gr] - 1}]]



(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
