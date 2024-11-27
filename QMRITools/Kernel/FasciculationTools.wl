(* ::Package:: *)

(* ::Title:: *)
(*QMRITools FasciculationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`FasciculationTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`FasciculationTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


FindActivations::usage = 
"FindActivations[data] Finds the activation in MUMRI or DTI data after data normalization. 
FindActivations[data, mask] Finds the activation in MUMRI or DTI data after data normalization within the mask."

EvaluateActivation::usage =
"EvaluateActivation[out] allows to evaluate the activation detection using FindActivations, where out is the output of that function with the option ActivationOutput set to True.
EvaluateActivation[out, actS] The same with the extra analysis of the SelectActivations function output given as actS."

AnalyzeActivations::usage = 
"AnalyzeActivations[actMap, mask] Analysis of the activation map generated from the mask."

SelectActivations::usage = 
"SelectActivations[act] selects the activations above the given ActivationSize.
SelectActivations[act, vox] selects the activations above the given ActivationSize where the activation size is in mm^3.

SelectActivations[act, mask] selects the activations above the given ActivationSize within the given mask or masks. The mask can be 3D or 4D.
SelectActivations[act, {mask, back}] selects the activations above the given ActivationSize within the given mask or masks. All voxels outside the back are ignored.

Output is {actSelected, actTotal} is mask is 3D.
Output is {{actSelected, Total[actSelected]}, {actTotal, Total[actTotal]}} is mask is 4D where actSelected and actTotal are per mask."


(* ::Subsection::Closed:: *)
(*Options*)


ActivationThreshold::usage =
"ActivationThreshold is an option for FindActivations. Fist value is the number of standard deviations second is the percentage threshold."

ThresholdMethod::usage =
"ThresholdMethod is an option for FindActivations. Values can be \"StandardDeviation\", \"Fraction\" or \"Both\"."

IgnoreSlices::usage =
"IgnoreSlices is an option for FindActivations and SelectActivations. Determins how many slices of the start and end of the dataset are ignored."

ActivationBackground::usage = 
"ActivationBackground is an option for FindActivations. If all normalized signals, which range between 0-150, are below this value the algorithm does nothing."

ActivationIterations::usage = 
"ActivationIterations is an option for FindActivations. The maximum number of iteration that can be used for activation detection."

ActivationOutput::usage = 
"ActivationOutput is an option for ActivationOutput. If set to All also the mn and threshold values are returned."

ActivationSize::usage = 
"ActivationSize is an option for SelectActivations. Its the size of the activations selected defined in number of voxels if no voxel size is given. If a voxel size is given its the volume.";


(* ::Subsection::Closed:: *)
(*Error Messages*)


FindActivations::tresh = 
"Given thresholds are not valid. The sd should be >1 and is `1` and the fr should be < 1 and is `2`." ;

SelectActivations::dim = 
"The data and/or mask dimensions are wrong. The activation should be 4D, the mask 3/4D, and the background mask 3D. Current dimensions are `1`D, `2`D, and `3`D respectively.";

SelectActivations::size = 
"The data and/or mask size are wrong. The activation should be the same size as the mask and the background mask. Current sizes are `1`, `2`, and `3` respectively.";

AnalyzeActivations::size = 
"The mask and activation map must have the same Dimensions.";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*FindActivations*)


(* ::Subsubsection::Closed:: *)
(*FindActivations*)


Options[FindActivations] = Options[FindActivationsI] = {
	ActivationThreshold -> {3.0, 0.6}, 
	ThresholdMethod -> "Both", 
	ActivationOutput -> "Activation",
	MaskDilation -> 0, 
	IgnoreSlices -> {0, 0},
	ActivationBackground ->10,
	ActivationIterations-> 10
};

SyntaxInformation[FindActivations] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

FindActivations[data_, ops : OptionsPattern[]] := FindActivationsI[NormalizeData[data, NormalizeMethod -> "Volumes"],  ops]

FindActivations[data_, mask_, ops : OptionsPattern[]] := FindActivationsI[
		NormalizeData[data, DilateMask[mask, OptionValue[MaskDilation]], NormalizeMethod -> "Volumes"], ops]

FindActivationsI[data_, OptionsPattern[]] := Block[{met, sc, fr, start, stop, dat, act, mn ,tr, itt, back},
	
	(*Get and check options*)
	met = OptionValue[ThresholdMethod];
	{sc, fr} = OptionValue[ActivationThreshold];
	If[sc < 0 || fr > 1, Return[Message[FindActivations::tresh, sc, fr]]];
	{start, stop} = OptionValue[IgnoreSlices];
	itt = OptionValue[ActivationIterations];
	back = OptionValue[ActivationBackground];
	
	(*set the threshold*)
	{sc, fr} = Switch[met,
		"Both", {sc, fr},
		"Fraction", {0, fr},
		"StandardDeviation", {sc, 1}
	];
	
	(*perform the activation finding in the selected slices*)
	dat = RotateDimensionsLeft[Transpose[data[[start + 1 ;; -stop - 1]]]];	
	act = FindActC[dat, sc, fr, itt, back];

	(*create extra output if needed*)
	If[OptionValue[ActivationOutput]=!="Activation",
		{mn, tr, sc, fr} = RotateDimensionsRight[MeanThresh[dat, act, sc, fr, back]];
		mn = ToPackedArray@ArrayPad[mn, {{start, stop}, {0, 0}, {0, 0}}, 0.];
		tr = ToPackedArray@ArrayPad[Transpose[{tr, sc, fr}], {{start, stop}, {0, 0}, {0, 0}, {0, 0}}, 0.];
	];
		
	(*give output*)
	act = SparseArray[ArrayPad[Round[Transpose[RotateDimensionsRight[act]]], {{start, stop}, {0, 0}, {0, 0}, {0, 0}}]];

	If[OptionValue[ActivationOutput]==="Activation", {act, data}, {act, data, mn ,tr}] 
]


(* ::Subsubsection::Closed:: *)
(*FindActC*)


FindActC = Compile[{{t, _Real, 1}, {sc, _Real, 0}, {fr, _Real, 0}, {it, _Real, 0}, {back, _Real, 0}}, Block[
	{tSelect, ti, cont, err, i, mn, tr, out},
	
	If[Mean[t] <= back,
		(*if background do nothing*)
		0 t,
		
		(*find activation function*)
		tSelect = ti = t;
		cont = True;
		err = False;
		i = 0;
		
		(*keep find activation till convergence*)
		While[cont && i <= it, 
			i++;
			(*select on fr only*)
			mn = Mean[ti];
			tSelect = Select[t, # > fr mn &];
			If[Length[tSelect] <= 5, 
				err = True; 
				cont = False;
				,			
				(*perform selection on sd and fr*)
				tr = Max[{0.1, Min[{1 - (sc/mn) StandardDeviation[tSelect], fr}]}];
				tSelect = Select[t, # > tr mn &];
				If[Length[tSelect] <= 5, 
					err = True;
					cont = False;
					,
					cont = (tSelect =!= ti);
					ti = tSelect;
				];
			];
		];
		
		If[err,
			(*if while get into error do nothing*)
			0 t,
			(*based on data vector without dropouts find correct threshold*)
			mn = Mean[ti];
			tr = Min[{1 - (sc/mn) StandardDeviation[ti], fr}];
			UnitStep[-t + tr mn]
		]
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization->True]


(* ::Subsubsection::Closed:: *)
(*MeanThresh*)


MeanThresh = Compile[{{t, _Real, 1},{s, _Real, 1}, {sc, _Real, 0}, {fr, _Real, 0}, {back, _Real, 0}}, Block[
	{ti, mn, sd, tr},
	If[Mean[t] <= back,
		(*if background do nothing*)
		{0, 0, 0, 0}, 
		ti = Select[(1-s) t, #>0.&];
		mn = Mean[ti];
		sd = Ramp[1 - (sc/mn) StandardDeviation[ti]];
		tr = Min[{sd, fr}];
		mn {1, tr, sd, fr}
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization->True]


(* ::Subsection:: *)
(*SelectActivations*)


(* ::Subsubsection::Closed:: *)
(*SelectActivations*)


Options[SelectActivations] = {ActivationSize->4, IgnoreSlices->{0,0}};

SyntaxInformation[EvaluateActivation2]={"ArgumentsPattern"->{_,_.,_.,OptionsPattern[]}};

SelectActivations[act_?ArrayQ, ops:OptionsPattern[]]:=SelectActivations[act,{1,1},{1,1,1}]

SelectActivations[act_?ArrayQ, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, ops:OptionsPattern[]]:=SelectActivations[act,{1,1},vox]

SelectActivations[act_ ?ArrayQ, mask_?ArrayQ, ops:OptionsPattern[]]:=SelectActivations[act,{mask,1},{1,1,1}]

SelectActivations[act_ ?ArrayQ, {mask_?ArrayQ, bmask_?ArrayQ}, ops:OptionsPattern[]]:=SelectActivations[act,{mask,bmask},{1,1,1}]

SelectActivations[act_?ArrayQ, {mask_,bmask_}, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, ops:OptionsPattern[]]:=Block[{
	aDepth,mDepth,bDepth,size, aDim,nVol, mDim, bDim, start, stop,back, masks,out,sel
	},
	(*check data dimensions*)
	aDepth = ArrayDepth[act];
	mDepth = If[mask=!=1,ArrayDepth[mask],aDepth-1];
	bDepth = If[bmask=!=1,ArrayDepth[bmask],aDepth-1];
	If[aDepth=!=4||4<mDepth<3||bDepth=!=3, Return[Message[SelectActivations::dim,aDepth,mDepth,bDepth]]];
	
	(*check data sizes*)
	aDim = Dimensions@act[[All,1]];
	nVol = Length@act[[1]];
	mDim = If[mask=!=1,Dimensions@If[mDepth===4,mask[[All,1]],mask],aDim];
	bDim = If[bmask=!=1,Dimensions[bmask],aDim];
	If[aDim=!=mDim || aDim=!=bDim, Return[Message[SelectActivations::size,aDim,mDim,bDim]]];
	
	(*get detection size of fasc*)
	size=Round[OptionValue[ActivationSize]/(Times@@vox )];
	
	(*create the selection mask*)
	sel=If[mask===1, 0act[[All,1]]+1, mask ];
	
	(*create background mask*)
	{start,stop}=OptionValue[IgnoreSlices];
	back=If[bmask===1,If[mDepth===4,Total@Transpose@sel,sel],bmask ];
	back[[start+1;;-stop-1]]=back[[start+1;;-stop-1]]+1;
	
	(*find the masks for which to perform the segmentation*)
	masks=If[mDepth===4,Transpose[MaskData[sel, back]],{sel back}];
	
	(*make the activation masks per dataset*)
	out=SelectActivationI[act,masks,size];
	
	If[mDepth===4,
		{{out,Total@out},{Transpose[Total[Transpose[#]]&/@out],Total@Transpose[Total[out]]}},
		{First@out,Total@Transpose@First@out}
	]
]


(* ::Subsubsection::Closed:: *)
(*SelectActivationI*)


SelectActivationI[act_?ArrayQ,masks_?ArrayQ,size_?IntegerQ]:=Block[{msk,dat},
	SparseArray[(
	(*loop over masks for analysis*)
	msk=#;
	dat=SparseArray[Transpose[msk #&/@Transpose[act]]];
	SparseArray[Round[Map[SelectActivationI[#,size]&,dat,{2}]]]
)&/@masks]];

(*actual selection function*)
SelectActivationI[im_?MatrixQ,size_?IntegerQ]:=If[Total[Flatten[im]]<size,
	0im,
	ImageData[SelectComponents[Image[im,"Bit"], #Count>=size&]]
]; 


(* ::Subsection:: *)
(*AnalyzeActivations*)


(* ::Subsubsection::Closed:: *)
(*AnalyzeActivations*)


AnalyzeActivations[act_,msk_]:=AnalyzeActivations[act,msk,""]

AnalyzeActivations[act_,msk_,lab_]:=Block[{aDepth,aDim,mDepth,mDim,labs},
	aDepth=ArrayDepth[act];
	aDim=Dimensions[If[aDepth===4,act[[All,1]],act[[1,All,1]]]];
	mDepth=ArrayDepth[msk];
	mDim=Dimensions[If[mDepth==3,msk,msk[[All,1]]]];
	
	If[aDepth=!=(mDepth+1)||aDim=!=mDim, Return[Message[AnalyzeActivations::size]]];
	
	If[mDepth===3,
		AnalyzeActivationsI[act,msk,lab],
		labs=If[lab===""||Length[lab]=!=Length[act],"Vol_"<>StringPadLeft[ToString[#],3,"0"]&/@Range[Length[act]],lab];
		Association[MapThread[AnalyzeActivationsI[#1,#2,#3]&,{act,Transpose@msk,labs}]]
	]
]


(* ::Subsubsection::Closed:: *)
(*AnalyzeActivationsI*)


AnalyzeActivationsI[act_,msk_,lab_]:=Block[{sizes,nActs,mSize,mSizeT,nSlices,nVols,nObs,mSd,quants,chance,chanceO,chanceV,vals,out},
	sizes=N@Flatten[Map[If[Total[Flatten[#]]<=1,{},ComponentMeasurements[Image[#,"Bit"],"Count"][[All,2]]]&,act,{2}]];
	
	nActs=Length@sizes;
	
	mSize=Map[Total[Flatten[#]]&,msk];
	mSizeT=Total@mSize;
	
	nSlices=Total@Unitize@mSize;
	nVols=Length@act[[1]];
	nObs=nSlices nVols;
	
	(*mSd=If[nActs>2,{Mean[sizes],StandardDeviation[sizes]},{0,0}];
	quants=If[nActs>2,Quantile[sizes, {0.5,0.05,0.95}],{0,0,0}];*)
	
	mSd = Which[
		nActs > 2, {Mean[sizes], StandardDeviation[sizes]},
		0 < nActs <=2, {Mean[sizes], 0.},
		True, {0., 0.}];
	quants = If[0 < nActs, Quantile[sizes, {0.5, 0.05, 0.95}], {0. ,0. ,0.}];

	chance=100. nActs/nVols;
	chanceO=100. nActs/nObs;
	chanceV = 1000. 100. nActs/nObs / mSizeT;
	
	vals=Flatten@{mSizeT,nActs,nObs,chance,chanceO,chanceV,mSd,quants};
	
	out=Association[Thread[{"ROI vol","Amount","Observed","Chance/Vol","Chance/Obs","Chance/Vox","Mean Size","StDv Size","Median Size","5% Size","95% Size"}->vals]];
	If[lab==="",out, Association[lab->out]]
]


(* ::Subsection::Closed:: *)
(*EvaluateActivation*)


SyntaxInformation[EvaluateActivation]={"ArgumentsPattern"->{_,_.,_.,_.,_.,OptionsPattern[]}};

EvaluateActivation[act_,dat_,mn_,tr_]:=EvaluateActivation[act,dat,mn,tr,act]

EvaluateActivation[{act_,dat_,mn_,tr_}]:=EvaluateActivation[act,dat,mn,tr,act]

EvaluateActivation[{act_,dat_,mn_,tr_},actS_]:=EvaluateActivation[act,dat,mn,tr,actS]

EvaluateActivation[act_,dat_,mn_,tr_,actS_]:=Module[{datD,actD,actSD,mnD,trD,sc,dim,aim
	(*ddim, zz, dd, yy ,xx, sl, dyn*)},
	NotebookClose[plotWindow];
	
	actD=Normal[act];
	actSD=Normal[actS];
	datD=dat;mnD=mn;trD=tr;
	(*prep images*)
	sc = Quantile[Flatten[datD],{0.005,0.995}];
	dim = N@Clip[Rescale[datD,sc],{0,1}];
	aim = RotateDimensionsLeft[N@{actD-actSD, actSD, 0.actD, Clip[actD+actSD,{0,1}]}];
	
	ddim = Dimensions[datD];
	{zz, dd, yy, xx} = ddim;
	{sl,dyn} = Ceiling[{zz,dd}/2];
		
	PrintTemporary["Prepping Manipulate window"];
	pan = Manipulate[
		(*slice location*)
		sl = Clip[sl,{1,zz},{1,zz}];
		dyn = Clip[dyn,{1,dd},{1,dd}];
		
		l = Length[pos[[m]]];
		n = Clip[n,{1,l}];
		
		If[!sel,
			(*use selection bar*)
			{z,y,x} = pos[[m,n]];
			sl = z;
			c = Abs[{x,y}-{0,yy+1}];
			,
			(*use locator pane*)
			cor = Reverse[Abs[Ceiling[c]-{0,yy+1}]];
			{y,x} = cor;
			z = sl;
		];
		
		(*get the signals*)
		ddat = Range[dd];
		sig = datD[[z,All,y,x]];
		actt = actD[[z,All,y,x]];
		actSt = actSD[[z,All,y,x]];
		
		(*get the line values*)
		mean = mnD[[z,y,x]];
		tresh = trD[[z,All,y,x]];
		
		Row[{Show[
			ListLinePlot[Thread[{ddat,sig}],PlotStyle->Black,PlotMarkers->Automatic,PlotRange->{0,1.1 Max@sig},ImageSize->500,
				GridLines->If[grid,{
					{{dyn,Directive[Black,Thick]}},
					{{mean,Directive[Black,Thick]},
						{tresh[[1]],Directive[Red,Thick]},
						{tresh[[2]],Directive[Thick,RGBColor[0.812807, 0.518694, 0.303459], Dashed]},
						{tresh[[3]],Directive[Thick,RGBColor[0.253651, 0.344893, 0.558151], Dashed]}}},
						None
					]
			],
			Graphics[{
				RGBColor[0.812807, 0.518694, 0.303459],Line[{{0,tresh[[2]]},{Length@sig,tresh[[2]]}}],
				RGBColor[0.253651, 0.344893, 0.558151],Line[{{0,tresh[[2]]},{Length@sig,tresh[[2]]}}]
			}],
			ListPlot[Pick[Thread[{ddat,sig}],actt,1],PlotStyle->Red,PlotMarkers->{Automatic,10}],
			ListPlot[Pick[Thread[{ddat,sig}],actSt,1],PlotStyle->Green,PlotMarkers->{Automatic,10}]]
			,
			LocatorPane[Dynamic[c],ImageCompose[
				Image[dim[[sl,dyn]],ColorSpace->"Grayscale",ImageSize->400],{Image[aim[[sl,dyn]],ColorSpace->"RGB"], alpha}
			], Appearance->Style[If[crs,"+"," "], Blue,FontSize->40]
		]}]
		
		,{{sel,False,"Use Locator"},{True,False}}
		,Delimiter
		,{{m,1,"Number of act."},max,ControlType->SetterBar}
		,{{n,1,"Position number"}, 1, Dynamic[l], 1}
		,Delimiter
		,{{sl, 1, "Slice"}, 1, Dynamic[zz], 1}
		,{{dyn, 1, "Dynamic"}, 1, Dynamic[dd], 1}
		,Delimiter
		,{{alpha, 0.5, "Opacity"},0,1}
		,{{crs, True,"Show cross"},{True,False}}
		,Delimiter
		,{{grid,False},{True,False}}
		,Row[{Dynamic[{z,y,x}]}]
		,
		{zz,ControlType->None},{dd,ControlType->None},{xx,ControlType->None},{yy,ControlType->None},
		{max,ControlType->None},{pos,ControlType->None},{l,1,ControlType->None},
		
		{z,ControlType->None},{y,ControlType->None},{x,ControlType->None},
		{c,ControlType->None},{cor,ControlType->None},
		{ddat,ControlType->None},{sig,ControlType->None},
		{actt,ControlType->None},{actSt,ControlType->None},
		{mean,ControlType->None},{tresh,ControlType->None},
		
		Initialization:>{
			{zz, dd, yy, xx} = ddim,
			{sl,dyn} = Ceiling[{zz,dd}/2],
			max = Select[Sort[DeleteDuplicates[Flatten@Total@Transpose@actSD]],#>0&],
			pos = Position[Round@Total@Transpose@actSD,#]&/@max
		},
		
		SynchronousInitialization->True,
		SaveDefinitions->False,
		ControlPlacement->Right
	];
	
	plotWindow = CreateWindow[DialogNotebook[{CancelButton["Close",Clear[datD,actD,actSD,mnD,trD,dim,aim];DialogReturn[]],pan},WindowSize->All,WindowTitle->"Plot data window"]];
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
