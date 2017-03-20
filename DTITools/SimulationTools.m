(* ::Package:: *)

(* ::Title:: *)
(*DTITools SimulationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`SimulationTools`", {"Developer`"}];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`SimulationTools`*"];
ClearAll @@ Names["DTITools`SimulationTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


AddNoise::usage = 
"AddNoise[data, noise] ads rician noise to the data with a given sigma or SNR value."

Tensor::usage = 
"Tensor[{l1, l2, l3}] creates a diffuison tensor with vectors {{0,0,1},{0,1,0},{1,0,0}} and eigenvalues {l1, l2, l3}.
Tensor[{l1, l2, l3}, {e1, e2, e3}] creates a diffuison tensor with vectors {e1, e2, e3} and eigenvalues {l1, l2, l3}.
Tensor[{l1, l2, l3}, \"Random\"] creates a diffuison tensor with random eigenvectors and eigenvalues {l1, l2, l3}.
Tensor[{l1, l2, l3}, \"OrtRandom\"] creates a diffuison tensor with random Orthogonal eigenvectors {{1,0,0},{0,1,0},{0,0,1}} and eigenvalues {l1, l2, l3}."

Signal::usage = 
"Signal[par,TR,TE] calculates the MRI signal at a given TR and TE. par is defineds as {pd, T1, T2}."

CreateData::usage = 
"CreateData[sig, eig, bvec, gradients, dim] creates a DTI datasets of dimensions dim with sig as unweighted signal S0 and bvec and gradients. \
eig can be {l1, l2, l3}, {{l1, l2, l3}, {e1, e2, e3}}, {{l1, l2, l3}, \"Random\"} or {{l1, l2, l3}, \"OrtRandom\"}. Uses Tensor internally."

SimParameters::usage = 
"SimParameters[tens] caculates the diffusion parameters for tens. The output can be used in PlotSimulationHist and PlotSimulation."

PlotSimulationHist::usage = 
"PlotSimulationHist[pars, label, xdata, tr] plots the pars (output form Parameters). \
Using label as plotlabel and xdata as x axis label. tr are the true parameter values."

PlotSimulation::usage = 
"PlotSimulation[pars, xval, true, label, color] plots the pars (output form Parameters). Using label as PlotLabel and xval as x axis Thics.\
tr are the true parameter values. color are the color used for the plot."

SimAngleParameters::usage = 
"SimAngleParameters[tens,vec] caculates the diffusion eigenvectors for tens compared to the true values vec. \
The output can be used in PlotSimulationAngleHist and PlotSimulationAngle."

PlotSimulationAngleHist::usage = 
"PlotSimulationAngleHist[pars, label, xdata] plots pars (output from Anlge Parameters)."

PlotSimulationAngle::usage = 
"PlotSimulationAngle[par, xdata, label, col] plots pars (output from Anlge Parameters)."

PlotSimulationVec::usage =
"PlotSimulationVec[tens, xdata, label] plots the eigenvectors from simulated tensors."



(* ::Subsection:: *)
(*Options*)


NoiseSize::usage = 
"NoiseSize is an option for AddNoise. Values can be \"Sigma\", then the noise sigma is given or \"SNR\", then the SNR is given."

TensOutput::usage = 
"TensOutput is an option for Tensor. Values can be \"Vector\" or \"Matrix\"."

SortVecs::usage = 
"SortVecs is an option for PlotSimulationVec."



(* ::Subsection:: *)
(*Error Messages*)


Tensor::vec = "Eigenvectors must be a 3x3 matrix of \"Random\", not: `1`"

Tensor::val = "Eigenvalues mus be a vector of size 3 or a number, not:`1`"

CreateData::eig = "eigen system must be 3 eigenvalues {l1,l2,l2} in which cases fixed vectors wil be used 
	or eigenvalues with geven vectors {{l1,l2,l3},{e1,e2,e3}} 
	other possibility is 3 eigenvalues with a random or random ortogonal (fixed vectors with random sign)
		{{l1,l2,l3},Random}
		{{l1,l2,l3},OrtRandom}
	not : `1`"

AddNoise::opt = "AddNoise"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]





(* ::Subsection::Closed:: *)
(*Definitions*)


sizes={200,300,400,500,750,1000,1500,2000,2500,3000};
files={".pdf",".jpg",".gif",".tif",".png"};

labStyle = Directive[Bold,FontFamily->"Helvetica",14,Black];
unit = " [\!\(\*SuperscriptBox[\(10\), \(-3\)]\) \!\(\*SuperscriptBox[\(mm\), \(2\)]\)/s]";
lambda[x_] :="\!\(\*SubscriptBox[\"\[Lambda]\", \"" <> ToString[x] <> "\"]\)";
epsilon[x_] := "\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \"" <> ToString[x] <> "\"]\)"

RicianDistribution=Compile[{{Mu, _Real}, {Sigma, _Real}}, Sqrt[RandomReal[NormalDistribution[Mu, Sigma]]^2 + RandomReal[NormalDistribution[0, Sigma]]^2]];

Phi[x_]:=1/(E^(x^2/2)*Sqrt[2*Pi]);
CapitalPhi[x_]:=.5(1+Erf[(x)/Sqrt[2]]);
SkewNorm[x_,Omega_,Xi_,Alpha_]:=(2/Omega)Phi[(x-Xi)/Omega]CapitalPhi[Alpha (x-Xi)/Omega];

HalfNorm[x_,Theta_]:=PDF[HalfNormalDistribution[Theta],x]


(* ::Subsection:: *)
(*Tensor*)


(* ::Subsubsection::Closed:: *)
(*Tensor*)


Options[Tensor]={TensOutput->"Vector"}

SyntaxInformation[Tensor] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

Tensor[l_,OptionsPattern[]]:=Tensor[l,{{0,0,1},{0,1,0},{1,0,0}},TensOutput->OptionValue[TensOutput]]

Tensor[l_,vec_,OptionsPattern[]]:=
Module[{e,tens},
	e=If[vec==="Random",
		RandomMat[],
		If[vec==="OrtRandom",
			OrtRandomMat[],
			If[MatrixQ[vec]&&ArrayDepth[vec]==2,
				vec,
				Return[Message[Tensor::vec,vec]]
				]
			]
		];
	tens=Chop[If[NumberQ[l],
		Transpose[e].{{l,0,0},{0,l,0},{0,0,l}}.e,
		If[VectorQ[l]&&Length[l]==3,
			Transpose[e].{{l[[1]],0,0},{0,l[[2]],0},{0,0,l[[3]]}}.e,
			Return[Message[Tensor::val,l]]
			]
		]];
	Switch[OptionValue[TensOutput],"Vector",TensVec[tens],"Matrix",tens]
	]


(* ::Subsubsection::Closed:: *)
(*Tensor functions*)


OrtRandomMat[]:=RandomSample[{{0,0,1},{0,1,0},{1,0,0}},3]


RandomMat[]:=
Module[{l1,l2,l3},
	l1={1,0,0};
	l2=Normalize[{0,1,1}*RandomVec[]];
	l3=Cross[l1,l2];
	{l1,l2,l3}.RotationMatrix[{{1,0,0},RandomVec[]}]
	];


RandomVec[]:=Normalize[RandomReal[NormalDistribution[],3]]


(* ::Subsection::Closed:: *)
(*Signal*)


SyntaxInformation[Signal] = {"ArgumentsPattern" -> {_, _, _}};

Signal[par_,TR_,TE_]:=par[[1]](1-Exp[-TR/par[[2]]])Exp[-TE/par[[3]]]


(* ::Subsection:: *)
(*CreateData*)


(* ::Subsubsection::Closed:: *)
(*CreateData*)


SyntaxInformation[CreateData] = {"ArgumentsPattern" -> {_, _, _, _, _.}};

CreateData[S0_,eig_,bval_Real,grad_,dim_]:=CreateData[S0,eig,Prepend[ConstantArray[bval,Length[grad]],0],Prepend[grad,{0,0,0}],dim]

CreateData[S0_,eig_,bvec:{_?NumberQ..},grad_,dim_]:=CreateData[S0,eig,Bmatrix[bvec,grad],dim]

CreateData[S0_,eig_,bmat_?ArrayQ,dim_]:=
Module[{diff},
	
	diff=If[Dimensions[eig]=={3},
		ConstantArray[SignalTensor[S0,bmat,Tensor[eig]],dim],
		If[Dimensions[eig]=={2,3}&&Dimensions[eig[[2]]]=={3,3},
			ConstantArray[SignalTensor[S0,bmat,Tensor[eig[[1]],eig[[2]]]],dim],
			If[eig[[2]]==="Random"||eig[[2]]==="OrtRandom",
			Array[SignalTensor[S0,bmat,Tensor[eig[[1]],eig[[2]]]]&,dim],
			Return[Message[CreateData::eig,eig]]]
			]
		];
	
	Switch[Length[dim],
		1,Transpose[diff],
		2,Transpose[diff,{2,3,1}],
		3,Transpose[diff,{1,3,4,2}]
		]
	]


(* ::Subsubsection::Closed:: *)
(*SignalTensor*)


SignalTensor[S0_, bmat_, D_] := 
Module[{Dv},
	Dv=Append[If[Dimensions[D]=={3,3},TensVec[D],D],Log[S0]];
	Exp[bmat.Dv]
	]


(* ::Subsection::Closed:: *)
(*AddNoise*)


Options[AddNoise]={NoiseSize->"Sigma"};

SyntaxInformation[AddNoise] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

AddNoise[dat_,noise_,OptionsPattern[]]:=
Module[{sig,data=dat//N},
	sig=Switch[OptionValue[NoiseSize],
		"SNR",
		Switch[ArrayDepth[data],
		1,
		data[[1]]/noise,
		2,
		Mean[DeleteCases[data[[1,All]],0.]]/noise,
		3,
		Mean[DeleteCases[Flatten[data[[1,All,All]]],0.]]/noise,
		4,
		Mean[DeleteCases[Flatten[data[[All,1,All,All]]],0.]]/noise,
		_,
		Message[AddNoise::dat];Return[]
		],
		"Sigma",
		noise,
		_,
		Message[AddNoise::opt];Return[]
		];
	(*Print["Sigma Noise = ",Round[Sigma,0.001]];*)
	Map[RicianDistribution[#,sig]&,data,{ArrayDepth[data]}]
	]


(* ::Subsection::Closed:: *)
(*SimParameters*)


Options[SimParameters]={Reject->False}

SyntaxInformation[SimParameters] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SimParameters[tens_,OptionsPattern[]]:=
Module[{eig,ADC,FA,dataAll,rangy,bins,wbins,sol,fit,x,Omega,Xi,Alpha},
	(
		Off[NonlinearModelFit::"cvmit"];
		Off[NonlinearModelFit::"sszero"];
		eig=1000 EigenvalCalc[#,Reject->OptionValue[Reject],MonitorCalc->False];
		ADC=ADCCalc[eig];
		FA=FACalc[eig];
		dataAll=DeleteCases[DeleteCases[Flatten[#],0],0.]&/@{eig[[All,All,1]],eig[[All,All,2]],eig[[All,All,3]],ADC,FA};
		rangy={{0,3},{0,3},{0,3},{0,3},{0,1}};
		bins=MapThread[{Range[#2[[1]]+.5(#2[[2]]/100),#2[[2]],#2[[2]]/100],BinCounts[#1,{#2[[1]],#2[[2]],#2[[2]]/100}]}&,{dataAll,rangy}];
		wbins=Transpose/@MapThread[{#1[[1]],(#1[[2]]/Length[#3])/(#2[[2]]/100)}&,{bins,rangy,dataAll}];
		sol=NonlinearModelFit[FitData[#],SkewNorm[x,Omega,Xi,Alpha],{Omega,Xi,Alpha},x,Gradient->"FiniteDifference"]&/@dataAll;
		fit=Append[ParameterFit[#]&/@dataAll,{Length[dataAll[[1]]]}];
		{dataAll,wbins,sol,fit}
		)&/@tens
	]


(* ::Subsection::Closed:: *)
(*PlotSimulationHist*)


SyntaxInformation[PlotSimulationHist] = {"ArgumentsPattern" -> {_, _, _, _}};

PlotSimulationHist[pars_,label_,xdata_,tr_]:=
DynamicModule[{rangy,xlabel,exp},
	rangy={{0,3},{0,3},{0,3},{0,3},{0,1}};
	xlabel={lambda[1]<>unit, lambda[2]<>unit, lambda[3]<>unit, "MD"<>unit, "FA [-]"};
	Manipulate[
		If[!ListQ[pars]||!ListQ[rangy],
			Return[],
			exp=GraphicsRow[(Show[
				Histogram[Flatten[pars[[y,1,#]]],{rangy[[#,2]]/100},"ProbabilityDensity",PlotRange->{rangy[[#]],{0,1.1Max[pars[[y,2,#,All,2]]]}},
				PerformanceGoal->"Speed",AxesOrigin->{0,0},LabelStyle->labStyle,
				FrameLabel->{xlabel[[#]],"Probability density"},Axes->False,FrameStyle->Thick,Frame->{True,True,False,False},ChartBaseStyle->EdgeForm[{Thin,White}],ChartStyle->Gray],
				ListPlot[pars[[y,2,#]],Joined->True,PlotStyle->{Thick,Black},PlotRange->{rangy[[#]],{0,1.1Max[pars[[y,2,#,All,2]]]}}],
				Plot[pars[[y,3,#]][x],{x,rangy[[#,1]],rangy[[#,2]]},PlotStyle->{Thick,Red},PlotRange->{rangy[[#]],{0,1.1Max[pars[[y,2,#,All,2]]]}}],
				ListLinePlot[{{tr[[#]],0},{tr[[#]],1.1Max[pars[[y,2,#,All,2]]]}},PlotStyle->Directive[Thick,Black,Dashed]]
				]&/@{1,2,3,4,5})[[xx]],ImageSize->Length[xx]*400,PlotLabel->Style[label<>"  -  "<>ToString[xdata[[y]]],16],LabelStyle->labStyle]
			]
		,{{xx,{1,2,3},"Parameter"},{{1,2,3}->"eigenvalues",{4,5}->"MD\\FA"}},{{y,1,"Simulation Value"},1,Length[pars],1},
		Button["Export Plot To File",FileSave[Dynamic[exp],"jpg",2000],Method->"Queued"],
		{{size, 500, "Export Size"}, sizes}, {{file, ".jpg","File Type"}, files},
		SaveDefinitions->True
		]
	]


(* ::Subsection::Closed:: *)
(*PlotSimulation*)


Options[PlotSimulation]={PlotRange->{{0,3},{0,3},{0,3},{0,3},{0,1}}};

SyntaxInformation[PlotSimulation] = {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};

PlotSimulation[pars_,xval_,tr_,label_,color_,OptionsPattern[]]:=
Module[{pl,dat,err,rangx,rangy,ylabel,truey,off},

	off=0.025(Max[xval]-Min[xval]);
	rangx={Min[xval]-off,Max[xval]+off};
	rangy=OptionValue[PlotRange];
	ylabel={lambda[1]<>unit, lambda[2]<>unit, lambda[3]<>unit, "MD"<>unit, "FA [-]"};
	pl=(
	truey=ConstantArray[tr[[#]],{Length[xval]}];
	dat=pars[[All,4,#,1]];
	err=If[Length[pars[[1,4,#]]]==2,pars[[All,4,#,2]],0];
	ListLinePlot[
		{Transpose[{xval,truey}],Transpose[{xval,dat}],Transpose[{xval,dat-err}],Transpose[{xval,dat+err}]},
		Filling->{3->{4}},FillingStyle->Directive[Opacity[0.1],color],FrameTicks->{xval,Automatic},PlotMarkers->{"",{"\[FilledSmallCircle]",15},{"\[UpPointer]",15},{"\[DownPointer]",15}},
		PlotStyle->{Directive[Gray,Thick,Dashing[Medium]],Directive[color,Thick],Directive[color,Thin,Dashing[Medium]],Directive[color,Thin,Dashing[Medium]]},
		PlotRange->{rangx,rangy[[#]]},
		LabelStyle->labStyle,FrameLabel->{label,ylabel[[#]]},Axes->False,FrameStyle->Thick,Frame->{True,True,False,False}
		]
	)&/@{1,2,3,4,5};
	GraphicsGrid[Partition[pl, 3, 3, 1, {}], ImageSize -> 1200]
	]


(* ::Subsection::Closed:: *)
(*AngleParameters*)


SyntaxInformation[SimAngleParameters] = {"ArgumentsPattern" -> {_, _}};

SimAngleParameters[tens_,veci_]:=
Module[{bins,surface,wbins,fit,par,vec,ang,x,Theta},
	surface=N[Table[( Cos[(x-1) Degree]- Cos[x Degree]),{x,1,90,1}]];
	(
		vec=EigenvecCalc[#,MonitorCalc->False];
		ang=AngleCalc[vec[[All,All,#]],veci[[#]],Distribution->"0-90"]&/@{1,2,3};
		bins=Map[BinCounts[DeleteCases[Flatten[#],0.],{0,90,1}]&,ang];
		wbins=Map[{#/Total[#],#/surface/Total[#/surface]}&,bins];
		fit=Map[NonlinearModelFit[#[[2]],HalfNorm[x,Theta],{Theta},x,Gradient->"FiniteDifference"]&,wbins];
		par=Map[({Theta}/.#["BestFitParameters"])[[1]]&,fit];
		{ang,wbins,fit,par}
		)&/@tens
	]


(* ::Subsection::Closed:: *)
(*PlotSimulationAngleHist*)


SyntaxInformation[PlotSimulationAngleHist] = {"ArgumentsPattern" -> {_, _, _}};

PlotSimulationAngleHist[pars_,label_,xdata_]:=
Module[{exp},
	Manipulate[
		If[!ListQ[pars],
			Return[],
			exp=GraphicsRow[
				Show[
					Histogram[Flatten[pars[[y,1,#]]],{1},"ProbabilityDensity",PlotRange->{{0,90},{0,1.1Max[pars[[y,2,#]]]}},PerformanceGoal->"Speed",LabelStyle->labStyle,
					FrameLabel->{"Error "<>epsilon[#]<>" [\[Degree]]","Probability density"},
					Axes->False,FrameStyle->Thick,Frame->{True,True,False,False},ChartBaseStyle->EdgeForm[{Thin,White}],ChartStyle->Gray],
					ListPlot[pars[[y,2,#,1]],Joined->True,PlotStyle->{Thick,Darker[Gray]},PlotRange->{{0,90},{0,1.1Max[pars[[y,2,#]]]}}],
					ListPlot[pars[[y,2,#,2]],Joined->True,PlotStyle->{Thick,Black},PlotRange->{{0,90},{0,1.1Max[pars[[y,2,#]]]}}],
					Plot[pars[[y,3,#]][x],{x,1,90},PlotStyle->{Thick,Red},PlotRange->{{0,90},{0,1.1Max[pars[[y,2,#]]]}}]
					]&/@{1,2,3}
				,ImageSize->1200,LabelStyle->labStyle,PlotLabel->Style[label<>" - "<>ToString[xdata[[y]]],16]]
			]
		,{{y,1,"Simulation Value"},1,Length[pars],1},
		Button["Export Plot To File",FileSave[Dynamic[exp],"jpg",2000],Method->"Queued"],
		{{size, 500, "Export Size"}, sizes}, {{file, ".jpg","File Type"}, files},
		SaveDefinitions->True
		]
	]


(* ::Subsection::Closed:: *)
(*PlotSimulationAngle*)


Options[PlotSimulationAngle]={PlotRange->{0,90}}

SyntaxInformation[PlotSimulationAngleHist] = {"ArgumentsPattern" -> {_, _, _, _, _., OptionsPattern[]}};

PlotSimulationAngle[par_, xdata_, label_, col_, OptionsPattern[]]:=
PlotSimulationAngle[par, xdata, label, col, {.25,.5,.95}, PlotRange->OptionValue[PlotRange]]

PlotSimulationAngle[par_, xdata_, label_, col_, quantinp_, OptionsPattern[]] := 
Module[{quant,pars=Transpose[par[[All,4]]],e,pdat,xrange,off},
	off=0.025(Max[xdata]-Min[xdata]);
	xrange = {Min[xdata]-off,Max[xdata]+off};
	GraphicsRow[(
    e = #;
    pdat = (
        quant = #;
        Transpose[{xdata, 
          Quantile[HalfNormalDistribution[#], quant] & /@ pars[[e]]}]
        ) & /@ quantinp;
    ListLinePlot[pdat,
     PlotStyle -> {Directive[col, Thick], 
       Directive[col, Thick, Dashing[Large]], 
       Directive[col, Thick, Dashing[Medium]], 
       Directive[col, Thick, Dashing[Small]], 
       Directive[col, Thick, Dashing[Tiny]]},
     FrameTicks->{xdata,Automatic},
     LabelStyle -> labStyle,
     PlotMarkers -> {"\[FilledSmallCircle]", 14},
     PlotRange -> {xrange,OptionValue[PlotRange]},
     FrameLabel -> {label, "Error "<>epsilon[#]<>" [\[Degree]]"}, Axes -> False,
     FrameStyle -> Thick,
     FrameTicks -> {Automatic, {0, 15, 30, 45, 60, 75, 90}},
     Frame -> {True, True, False, False}]
    ) & /@ {1, 2, 3}
 , ImageSize -> 1000, Spacings -> 0]
]


(* ::Subsection:: *)
(*PlotSimulationAngle*)


(* ::Subsubsection::Closed:: *)
(*PlotSimulationAngle*)


Options[PlotSimulationVec]={SortVecs->True}

SyntaxInformation[PlotSimulationVec] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

PlotSimulationVec[tens_, xdata_, label_, OptionsPattern[]] := Module[
  {eigvm, eigvcm, arrows1, arrows2, arrows3, vp1, vp2, vp3, 
   va1, va2, va3, vv1, vv2, vv3, data, exp, xd},
  
  eigvm = Flatten[EigenvecCalc[#,MonitorCalc->False], 1] & /@ tens;
  
  eigvcm =If[OptionValue[SortVecs],
  Map[(
      {
       RandomReal[NormalDistribution[1.025, .025]] If[
         Negative[#[[1, 3]]], {1, -1, -1}, {1, -1, 1}]*#[[1]],
       RandomReal[NormalDistribution[1.025, .025]] If[
         Negative[#[[2, 2]]], {1, 1, 1}, {1, -1, 1}]*#[[2]],
       RandomReal[NormalDistribution[1.025, .025]] If[
         Negative[#[[3, 1]]], {-1, -1, 1}, {1, -1, 1}]*#[[3]]
       }
      ) &, eigvm, {2}]
  ,
  eigvm
  ];
  
  arrows1 = Graphics3D[
    {{Darker[Blue], Arrowheads[0.05], 
      Arrow[Tube[{{-1, 1, 0.5}, {-1, 1, 1}}, 0.025]]},
     {Darker[Red], Arrowheads[0.05], 
      Arrow[Tube[{{-1, 1, 0.5}, {-1, .5, 0.5}}, 0.025]]},
     {Darker[Green], Arrowheads[0.05], 
      Arrow[Tube[{{-1, 1, 0.5}, {-.5, 1, 0.5}}, 0.025]]
      }}];
  arrows2 = Graphics3D[{
     {Darker[Blue], Arrowheads[0.05], 
      Arrow[Tube[{{-1, -.5, -1}, {-1, -.5, -.5}}, 0.025]]},
     {Darker[Red], Arrowheads[0.05], 
      Arrow[Tube[{{-1, -.5, -1}, {-1, -1, -1}}, 0.025]]},
     {Darker[Green], Arrowheads[0.05], 
      Arrow[Tube[{{-1, -.5, -1}, {-.5, -.5, -1}}, 0.025]]
      }}];
  arrows3 = Graphics3D[{
     {Darker[Blue], Arrowheads[0.05], 
      Arrow[Tube[{{.5, 1, -1}, {.5, 1, -.5}}, 0.025]]},
     {Darker[Red], Arrowheads[0.05], 
      Arrow[Tube[{{.5, 1, -1}, {.5, .5, -1}}, 0.025]]},
     {Darker[Green], Arrowheads[0.05], 
      Arrow[Tube[{{.5, 1, -1}, {1, 1, -1}}, 0.025]]
      }}];
    
  vp1 = {0.426945, -0.474858, 3.32298}; 
  vv1 = {0.0421748, -0.881393, 0.470498}; va1 = 25 Degree;
  (*{Dynamic[vp1],Dynamic[vv1],Dynamic[va1]}*)
  vp2 = {0.429321, -3.30821, 0.56695}; 
  vv2 = {0.058434, -0.443718, 0.894259}; va2 = 25 Degree;
  (*{Dynamic[vp2],Dynamic[vv2],Dynamic[va2]}*)
  vp3 = {3.1877, 0.899157, 0.692886}; 
  vv3 = {0.416537, 0.112532, 0.902127}; va3 = 25 Degree;
  (*{Dynamic[vp3],Dynamic[vv3],Dynamic[va3]}*)

  Manipulate[
  	If[!ListQ[eigvcm],Return[],
	   data = eigvcm[[set]] // Transpose;
	   xd = xdata[[set]];
	      
	   exp = GraphicsRow[{
	      EigPlot[data, vp1, vv1, va1, arrows1, 1], 
	      EigPlot[data, vp2, vv2, va2, arrows2, 2], 
	      EigPlot[data, vp3, vv3, va3, arrows3, 3]
	      }, 
	      ImageSize -> {1000}, Spacings -> 0,
	      PlotLabel-> Style[label <> " - " <>ToString[xdata[[set]]],18],
	      LabelStyle -> labStyle
      ]
  	]
   ,
   {{set, 1, "Simulation Value"}, 1, Length[xdata], 1},
   Button["Export Plot To File", FileSave[exp, "jpg", 2000], Method -> "Queued"],
   {{size, 500, "Export Size"}, sizes}, {{file, ".jpg","File Type"}, files},
   SaveDefinitions->True
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*PlotSimulationAngle*)


EigPlot[data_, vp_, vv_, va_, arrows_, val_] := Module[{sphere,line},
	 sphere = SphericalPlot3D[.975, {Theta, 0, Pi}, {Phi, 0, 2 Pi}, Lighting -> "Neutral"];
	 line = ParametricPlot3D[{{Sin[u], 0, Cos[u]}, -{0, Sin[u], Cos[u]}}, {u, -Pi, Pi},
	 	PlotStyle -> {Directive[Dashed, Thickness[.0125], Darker[Green]], Directive[Dashed, Thickness[.0125], Darker[Red]]}
	 	];
	 	Show[
	 		ListPointPlot3D[data, 
	 			ViewPoint -> vp, ViewVertical -> vv, ViewAngle -> va,
	 			PlotLabel -> {"First", "Second", "Third"}[[val]]<>" eigenvector",
	 			PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}, {-1.1, 1.1}}, 
			   BoxRatios -> 1, PlotStyle -> {Darker[Blue], Darker[Red], Darker[Green]}, 
			   Lighting -> "Neutral", Axes -> False, Boxed -> False, 
			   SphericalRegion -> True, 
			   LabelStyle -> labStyle
	 			], arrows, sphere, line
	 			]
    		];


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`SimulationTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]& /@ Names["DTITools`SimulationTools`*"];

EndPackage[]
