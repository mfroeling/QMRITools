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


(* ::Subsection:: *)
(*Functions*)


FindActivations::usage = 
"FindActivations[data] Finds the activation in MUMRI or DTI data after data normalization. 
FindActivations[data, mask] Finds the activation in MUMRI or DTI data after data normalizeation within the mask."

EvaluateActivation::usage =
"EvaluateActivation[out] allows to evaluate the activation deterction using FindActivations, where out is the output of that function with the option Activationoutput set to True.
EvaluateActivation[out, actS] The same with the extra annalysis of the SelectActivations funcion output given as actS."


(* ::Subsection:: *)
(*Options*)


ActivationThreshold::usage =
"ActivationThresholdis an option for FindActivations. Fist value is the number of standard deviations second is the pecentage threshold."

ThresholdMethod::usage =
"ThresholdMethodis an option for FindActivations. Values can be \"StandardDeviation\", \"Fraction\" or \"Both\"."

IgnoreSlices::usage =
"IgnoreSlices is an option for FindActivations. Determins how many slices of the start and end of the dataset are ignored."

ActivationOutput::usage = 
"ActivationOutput is an option for ActivationOutput. If set to All aslo the mn and treshhold values are retured."


(* ::Subsection:: *)
(*Error Messages*)


FindActivations::tresh = "Given thresholds are not valid. The sd should be >1 and is `1` and the fr should be < 1 and is `2`." ;


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
    IgnoreSlices -> {0, 0}
};

SyntaxInformation[FindActivations] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

FindActivations[data_, ops : OptionsPattern[]] := FindActivationsI[NormalizeData[data, NormalizeMethod -> "Volumes"], ops]

FindActivations[data_, mask_, ops : OptionsPattern[]] := FindActivationsI[
		NormalizeData[data, DilateMask[mask, OptionValue[MaskDilation]], NormalizeMethod -> "Volumes"], ops]

FindActivationsI[data_, OptionsPattern[]] := Block[{met, sc, fr, start, stop, dat, act, mn ,tr},
	
	(*Get and check options*)
	met = OptionValue[ThresholdMethod];
	{sc, fr} = OptionValue[ActivationThreshold];
	If[sc < 1 || fr > 1, Return[Message[FindActivations::tresh, sc, fr]]];
	{start, stop} = OptionValue[IgnoreSlices];
	
	(*set the threshold*)
	{sc, fr} = Switch[met,
		"Both", {sc, fr},
		"Fraction", {0, fr},
		"StandardDeviation", {sc, 1}
	];
	
	(*perfomr the activation finding in the selected slices*)
	dat = RotateDimensionsLeft[Transpose[data[[start + 1 ;; -stop - 1]]]];
	act = FindActC[dat, sc, fr];
	
	(*create extra ouput if needed*)
	If[OptionValue[ActivationOutput]=!="Activation",
		{mn, tr, sc, fr} = RotateDimensionsRight[MeanTresh[dat, act, sc, fr]];
		mn = ToPackedArray@N@ArrayPad[mn, {{start, stop}, 0, 0}];
		tr = ToPackedArray@N@ArrayPad[Transpose[{tr, sc, fr}], {{start, stop}, 0, 0, 0}];
	];
	
	(*give outpu*)
	act = ToPackedArray@Round@ArrayPad[Transpose[RotateDimensionsRight[act]], {{start, stop}, 0, 0, 0}];
	If[OptionValue[ActivationOutput]==="Activation", {act, data}, {act, data, mn ,tr}] 
  ]


(* ::Subsubsection::Closed:: *)
(*FindActC*)


FindActC = Compile[{{t, _Real, 1}, {sc, _Real, 0}, {fr, _Real, 0}}, Block[{ts, ti, c, i, mn, tr},
	If[Total[t] <= 0.,
		(*if backgroud do nothing*)
		t, 
		(*find activation function*)
		ts = ti = t;
		c = True;
		i = 0;
		
		(*keep find activation till convergence*)
		While[c && i < 10, i++;
			mn = Mean[ti];
			ts = Select[t, # > fr mn &];
			tr = Max[{0., Min[{1 - sc StandardDeviation[ts/mn], fr}]}];
			ts = Select[t, # > tr mn &];
			c = (ts =!= ti);
			ti = ts;
		];
		
		(*based on data vector without dropouts find correct thresshold*)
		mn = Mean[ti];
		tr = Max[{0.1, Min[{1 - sc StandardDeviation[ti/mn], fr}]}];
		(*the activations*)
		UnitStep[-t + tr mn]
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*MeanTresh*)


MeanTresh = Compile[{{t, _Real, 1},{s, _Real, 1}, {sc, _Real, 0}, {fr, _Real, 0}}, Block[{ts, ti, c, i, mn, tr},
	If[Total[t] <= 0.,
		(*if backgroud do nothing*)
		{0, 0, 0, 0}, 
		ti = Pick[t, s, 0.];
		mn = Mean[ti];
		tr = Max[{0.1, Min[{1 - sc StandardDeviation[ti/mn], fr}]}];
		mn {1, tr, (1 - sc StandardDeviation[ti/mn]), fr}
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection:: *)
(*EvaluateActivation*)


SyntaxInformation[EvaluateActivation2]={"ArgumentsPattern"->{_,_.,_.,_.,_.,OptionsPattern[]}};

EvaluateActivation[act_,dat_,mn_,tr_]:=EvaluateActivation[act,dat,mn,tr,act]

EvaluateActivation[{act_,dat_,mn_,tr_}]:=EvaluateActivation[act,dat,mn,tr,act]

EvaluateActivation[{act_,dat_,mn_,tr_},actS_]:=EvaluateActivation[act,dat,mn,tr,actS]

EvaluateActivation[act_,dat_,mn_,tr_,actS_]:=Module[{datD,actD,actSD,mnD,trD,sc,dim,aim,
	ddim},
	NotebookClose[plotwindow];
	
	datD=dat;actD=act;actSD=actS;mnD=mn;trD=tr;
	(*prep images*)
	sc = Quantile[Flatten[datD],{0.005,0.995}];
	dim = N@Clip[Rescale[datD,sc],{0,1}];
	aim = RotateDimensionsLeft[N@{actSD,0.actD,actD-actSD,Clip[actD+actSD,{0,1}]}];
	
	ddim = Dimensions[datD];
		
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
				GridLines->{{{dyn,Directive[Black,Thick]}},{{mean,Directive[Black,Thick]},{tresh[[1]],Red},{tresh[[2]],Directive[Thick,Gray]},{tresh[[3]],Directive[Thick,Gray,Dashed]}}}],
			ListPlot[Pick[Thread[{ddat,sig}],actt,1],PlotStyle->Blue,PlotMarkers->{Automatic,10}],
			ListPlot[Pick[Thread[{ddat,sig}],actSt,1],PlotStyle->Red,PlotMarkers->{Automatic,10}]]
			,
			LocatorPane[Dynamic[c],ImageCompose[
				Image[dim[[sl,dyn]],ColorSpace->"Grayscale",ImageSize->400],{Image[aim[[sl,dyn]],ColorSpace->"RGB"], alpha}
			],Appearance->Style[If[crs,"+"," "],Green,FontSize->40]
		]}]
		
		,{{sel,False,"Use Locator"},{True,False}}
		,Delimiter
		,{{m,1,"Number of act."},max,ControlType->SetterBar}
		,{{n,1,"Position number"},1,Dynamic[l],1}
		,Delimiter
		,{{sl, Ceiling[zz/2], "Slice"},1,zz,1}
		,{{dyn, Ceiling[dd/2], "Dynamic"},1,dd,1}
		,Delimiter
		,{{alpha, 0.5, "Opacity"},0,1}
		,{{crs,True,"Show cross"},{True,False}}
		
		,
		{zz,ControlType->None},{dd,ControlType->None},{xx,ControlType->None},{yy,ControlType->None},
		{max,ControlType->None},{pos,ControlType->None},{l,1,ControlType->None},
		
		{z,ControlType->None},{y,ControlType->None},{x,ControlType->None},
		{c,ControlType->None},{cor,ControlType->None},
		{ddat,ControlType->None},{sig,ControlType->None},
		{actt,ControlType->None},{actSt,ControlType->None},
		{mean,ControlType->None},{tresh,ControlType->None},
		
		Initialization:>{
			{zz,dd,yy,xx} = ddim,
			{sl,dyn} = ddim[[1;;2]],
			max = Select[Sort[DeleteDuplicates[Flatten@Total@Transpose@actSD]],#>0&],
			pos = Position[Round@Total@Transpose@actSD,#]&/@max
		},
		
		SynchronousInitialization->False,
		SaveDefinitions->False,
		ControlPlacement->Right
	];
	
	plotwindow=CreateWindow[DialogNotebook[{CancelButton["Close",Clear[datD,actD,actSD,mnD,trD,dim,aim];DialogReturn[]],pan},WindowSize->All,WindowTitle->"Plot data window"]];
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
