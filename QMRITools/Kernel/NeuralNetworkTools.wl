(* ::Package:: *)

(* ::Title:: *)
(*QMRITools NeuralNetworkTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`NeuralNetworkTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`NeuralNetworkTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


MakeUnet::usage = 
"MakeUnet[nClasses, dimIn] Generates a UNET with one channel as input and nClasses as output. 
MakeUnet[nChannels, nClasses, dimIn] Generates a UNET with nChannels as input and nClasses as output. 
he number of parameter of the first convolution layer can be set with dep. The data dimensions can be 2D or 3D and each 
of the dimensions should be 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240 or 256. However dimensions can be different
based on the network depth and the block type. The implemented block types are \"Conv\", \"UNet\", \"ResNet\", \"DenseNet\", \"Inception\", or \"U2Net\"."

MakeNode::usage = 
"MakeNode[scale, conn, blockConfig] makes a node for a UNET. 
The input scale defines the input and output scaling, is either an integer or a vector of length dim. 
The input conn defines the connections, is a list of two integer values defining the number of input and output ports.
The blockConfig is defined as {{blockType, settings}, {features,..}, {act, dim}}."

ActivationLayer::usage=
"ActivationLayer[type] is a wrapper around ElementwiseLayer needed in MakeNode."

MakeClassifyNetwork::usage = 
"MakeClassifyNetwork[classes] makes a classify network with three convolution layers and 3 fully connected layers. 
The input classes should be a list of strings. The input image dimensions should not be smaller than 64x64."

MakeClassifyImage::usage =
"MakeClassifyImage[data] makes a image of the input data. The data is automatically cropped to remove the background and normalized.
If the input data is 3D a list of images is returned."


NetSummary::usage = 
"NetSummary[net] gives a short summary of the convolution kernels and array elements in the network.
NetSummary[net, what] does the same but what can be \"Full\" which also includes net and node images or \"Mem\" which only reports the memory."

NetDimensions::usage = 
"NetDimensions[net] extracts the input channels, output classes, the input patch dimension, and the number of input filters."

ChangeNetDimensions::usage =
"ChangeNetDimensions[netIn] changes input channels, output classes, the input patch dimension of the input network netIn."

GetNetNodes::usage = 
"GetNetNodes[net] returns a list of all the nodes in the network net."


AddLossLayer::usage = 
"AddLossLayer[net] adds loss layers to a NetGraph. The DiceLossLayer, JaccardLossLayer, TverskyLossLayer, MSDLossLayer, TopK, and CELossLayer are added."

DiceLossLayer::usage = 
"DiceLossLayer[] represents a net layer that computes the Dice loss by comparing input class probability vectors with the target class vector.
DiceLossLayer[n] does the same but n defines the power of the denominator, with n=2 the squared dice score, is calculated."

JaccardLossLayer::usage =
"JaccardLossLayer[] represents a net layer that computes the Jaccard loss by comparing input class probability vectors with the target class vector.
JaccardLossLayer[n] does the same but n defines the power of the denominator, with n=2 the squared Jaccard score is calculated."

TverskyLossLayer::usage =
"TverskyLossLayer[] represents a net layer that computes the Tversky loss by comparing input class probability vectors with the target class vector.
TverskyLossLayer[b] does the same but b defines the Tversky beta factor. With beta = 0.5 its is the Dice coefficient. Here alpha + beta = 1."

OverlapLossFunction::usage =
"OverlapLossFunction[] is a generalization of TverskyLossLayer that can also generate Jaccard and Dice."

FocalLossLayer::usage =
"FocalLossLayer[] represents a net layer that computes the Focal loss by comparing input class probability vectors with the target class vector.
FocalLossLayer[g] does the same but uses g as the tunable focusing parameter gamma which needs to be larger than one.
FocalLossLayer[g, a] does the same but uses as the balancing factor alpha."

TopKLossLayer::usage =
"TopKLossLayer[net] represents a net layer that computes the topK 10% loss.
TopKLossLayer[net, k] does the same but k defines the topK between 0 and 1."


ClassEncoder::usage = 
"ClassEncoder[label] encodes Integer label data of 0 to max value of label into a nClass + 1 vector of 1 and 0 as the last dimension.
ClassEncoder[label, nClass] encodes Integer label data of 0 to nClass into a nClass + 1 vector of 1 and 0 as the last dimension."

ClassDecoder::usage = 
"ClassDecoder[probability] decodes a probability vector of 1 and 0 into Integers of 0 to the value of the last dimension of probability minus one.
ClassDecoder[probability, nClass] decodes a probability vector of 1 and 0 into Integers of 0 to nClass - 1."


AnalyzeNetworkFeatures::usage = 
"AnalyzeNetworkFeatures[net, data] gives overview of the information density of the network features by analyzing them with SVD."


(* ::Subsection::Closed:: *)
(*Options*)


BlockType::usage = 
"BlockType is an option for MakeUnet. It specifies the type of block used in the network. It can be \"Conv\", \"UNet\", \"ResNet\", \"DenseNet\", \"Inception\", or \"U2Net\"."

DropoutRate::usage = 
"DropoutRate is an option for MakeUnet. It specifies how much dropout is used after each block. It is a value between 0 and 1, default is .2."

RescaleMethod::usage =
"RescaleMethod is an option for MakeUnet. It specifies how the network rescales. It can be \"Conv\" or \"Pool\"."

CatenateMethod::usage =
"CatenateMethod is an option for MakeUnet. It specifies how the network catenates the skip and upscale channel. It can be \"Tot\", or \"Cat\"."

NetworkDepth::usage = 
"NetworkDepth is an option for MakeUnet. It specifies how deep the UNET will be."

DownsampleSchedule::usage = 
"DownsampleSchedule is an option for MakeUnet. It defines how the data is downsampled for each of the deeper layers of the Unet. 
By default is is a factor two for each layer. A custom schedule for a 5 layer 3D Unet could be {{2,2,2},{1,2,2},{2,2,2},{1,2,2}, 1}.
The deepest layer is always downsampled by 1 and therefore not needed to be specified."

SettingSchedule::usage =
"SettingSchedule is an option for MakeUnet. It defines the settings for the Unet blocks. If one setting is given it applied to all layers.
If a list of settings is given the settings can be different per layer. The following settings are the default settings. 
\"Unet\": convolution block repetitions, 2, \"ResNet\" -> convolution block repetitions, 2, \"DenseNet\" -> {dense depth, block repetitions}, {4,2},
\"Inception\" -> {inception width, block repetitions}, {4,2}, \"U2Net\"-> {Unet depth, downscale}, {5, True}."

FeatureSchedule::usage = 
"FeatureSchedule is an option for MakeUnet. It defines how the number of features is up-sampled for each of the deeper layers of the Unet.
By default it increases the number of features by a factor 2 each layer, i.e. {1, 2, 4, 8, 16}."

NetworkArchitecture::usage = 
"NetworkArchitecture is an option for MakeUnet. It defines the architecture of the network. It can be \"UNet\", \"UNet+\", or \"UNet++\".
For \"UNet+\" or \"UNet++\" it can also be {arch, i} where i specifies how many of the top layers are connected to the mapping layer."

ActivationType::usage = 
"ActivationType is an option for MakeUnet. It specifies which activation layer is used in the network. It can be \"LeakyRELU\" or any type allowed 
by a \"name\" definition in ElementwiseLayer."

NormalizationType::usage = 
"NormalizationType is an option for MakeUnet. It specifies which normalization layer is used in the network. It can be \"Batch\" or \"Instance\" for 
either BatchNormalization or NormalizationLayer which is instance normalization."


(* ::Subsection::Closed:: *)
(*Error Messages*)


MakeUnet::scale = "The scaling input is not valid. It should be Automatic. It can also be a integer or a list of integers that will be applied to the Layers. 
It can also be a vector of integers per layer where the length of the vector should be equal to the depth of the network.";

MakeUnet::sett = "The setting input is not valid. It can be a number or a list of numbers that will be applied to the Layers.";

MakeUnet::drop = "The dropout input is not valid. It can be a number or a list of numbers that will be applied to the Layers.";

MakeUnet::feat = "The feature input is not valid. It can be a number or a list of numbers that will be applied to the Layers.";

MakeUnet::arch = "The architecture input is not valid. It can be \"UNet\", \"UNet+\", or \"UNet++\".";

MakeUnet::block = "The block type input is not valid. It can be \"Conv\", \"UNet\", \"ResNet\", \"DenseNet\", \"Inception\", or \"U2Net\".";


ActivationLayer::type = "Not a correct activation layer `1`";


SurfaceDistance::met = "Method `1` not recognized";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


debugUnet[x___] := If[$debugUnet, Print[x]];


(* ::Subsection:: *)
(*MakeUnet*)


(* ::Subsubsection::Closed:: *)
(*Settings*)


netDefaults={
	"Conv"->1, 
	"UNet"->2, 
	"ResNet"->2, 
	"ResNetL"->2, 
	"DenseNet"->{4,2}, 
	"Inception"->{4,2}, 
	"U2Net"->{3,True}, 
	"Map"->1
};


(* ::Subsubsection::Closed:: *)
(*MakeUnet*)


Options[MakeUnet] = {
	NetworkArchitecture -> "UNet",
	BlockType -> "ResNet", 
	ActivationType -> "GELU",
	NormalizationType -> "Instance",
	RescaleMethod -> "Conv",
	CatenateMethod -> "Cat",
	
	NetworkDepth -> 5,
	DownsampleSchedule -> Automatic,
	SettingSchedule -> Automatic,
	FeatureSchedule -> 32,

	DropoutRate -> 0.2, 
	MonitorCalc -> False
};

SyntaxInformation[MakeUnet] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeUnet[nClass_?IntegerQ, dimIn_, opts : OptionsPattern[]] := MakeUnet[1, nClass, dimIn, opts]

MakeUnet[nChan_?IntegerQ, nClass_?IntegerQ, dimIn_, OptionsPattern[]] := Block[{
		architecture, blockType, actType, normType, metSc, metCat,
		drop, depth, scaling, feature, setting, mon,
		ndim, nam, conf, net, boolSc, cons, depthj, fStart, mapCon
	},

	{architecture, blockType, actType, normType, metSc, metCat} = OptionValue[{
		NetworkArchitecture, BlockType, ActivationType, NormalizationType, RescaleMethod, CatenateMethod}];
	{drop, depth, scaling, feature, setting, mon} = 
		OptionValue[{DropoutRate, NetworkDepth, DownsampleSchedule, FeatureSchedule, SettingSchedule, MonitorCalc}];
	mon = If[mon, MonitorFunction, List];

	(*check the input*)
	{architecture, mapCon} = If[StringQ[architecture], {architecture, Automatic}, 
		If[Length[architecture]===2, architecture, Return[Message[MakeUnet::arch]; $Failed]]];
	If[!MemberQ[{"UNet", "UNet+", "UNet++"}, architecture], 
		Return[Message[MakeUnet::arch]; $Failed]];
	If[!MemberQ[{"Conv", "UNet", "ResNet", "ResNetL", "DenseNet", "Inception", "U2Net"}, blockType], 
		Return[Message[MakeUnet::block]; $Failed]];
	conf = {actType, normType};

	(*is the network 2D or 3D*)
	ndim = Length@dimIn;

	mon[{architecture, blockType, actType, normType}, "Network block type: "];
	mon[{dimIn, ndim}, "Network dimension order: "];

	drop = Which[
		drop === False || drop === None, ConstantArray[0., depth],
		NumberQ[drop], ConstantArray[drop, depth],
		ListQ[drop],
			If[Length@drop =!= depth, Message[MakeUnet::drop]];
			PadRight[drop, depth, Last[drop]]
	];

	(*define the scaling, automatic is 2 for all layers, integer is for all layers, list is padded with last value*)
	scaling = Which[
		scaling === Automatic, ConstantArray[2, depth],
		IntegerQ[scaling], ConstantArray[scaling, depth],
		ListQ[scaling],	
			If[Length@scaling =!= depth, Message[MakeUnet::scale]];
			PadRight[scaling, depth, 1],
		True, Return[Message[MakeUnet::scale]; $Failed]
	];

	boolSc = (Times @@ If[IntegerQ[#], ConstantArray[#, ndim], #])=!=1&/@scaling;
	mon[scaling, "Network scaling schedule: "];

	(*define the setting, can be Automatic, number or {set}, or list of settings*)
	setting = Which[
		setting === Automatic, ConstantArray[blockType /. netDefaults, depth],
		IntegerQ[setting], ConstantArray[setting, depth],
		ListQ[setting], Which[
			Length[setting] == 2 && VectorQ[setting],
				ConstantArray[setting, depth],
			ListQ[setting],
				If[Length@setting =!= depth, Message[MakeUnet::sett]];
				PadRight[setting, depth, Last[setting]],
			True, Return[Message[MakeUnet::sett]; $Failed]
		]
	];
	mon[setting, "Network setting schedule: "];

	(*define the number of features per layer*)
	feature = Round[Which[
		IntegerQ[feature], feature Switch[blockType,
			"UNet" | "ResNet" | "ResNetL", 2^Range[0, depth - 1],
			"DenseNet" | "Inception", 0.5 Range[2, depth + 1],
			"U2Net", ConstantArray[1, depth]
		],
		ListQ[feature], If[Length[feature] == 2 && VectorQ[feature],
			ConstantArray[feature, depth], 
			If[Length@feature =!= depth, Message[MakeUnet::feat]];
			PadRight[feature, depth, Last[feature]]
		],
		True, Return[Message[MakeUnet::feat]; $Failed]
	]];
	If[VectorQ[feature], feature = Round[Transpose[{feature, Prepend[Most[feature], First[feature]/2]}]]];
	fStart = Last@Flatten@{First@feature};
	mon[feature, "Network feature schedule: "];

	If[metSc === "ConvS", blockType = blockType /. {"ResNet"->"ResNetS"}];

	(*make the network*)
	nam = If[StringQ[#1], #1<>ToString[#2], "node_"<>ToString[#1]<>"_"<>ToString[#2]]&;
	Switch[architecture,
		"UNet",
		net = NetGraph[
			(*define the blocks*)
			Sort@Join[
				(*encoding layers*)
				Table[debugUnet[nam["enc_", i]]; nam["enc_", i] -> MakeNode[
					(*scale up -> never, scale down*)
					{If[i==depth, scaling[[i]], False], scaling[[i]]}, 
					(*skip in -> only input from above, skip out -> always for enc*)
					{0, i=!=depth}, 
					(*config*)
					{{blockType, setting[[i]]}, feature[[i]], {conf, ndim}, scaling[[i]]},
					DropoutRate -> drop[[i]], RescaleMethod -> metSc, CatenateMethod -> metCat
				], {i, 1, depth}],
				(*decoding layers*)
				Table[debugUnet[nam["dec_", i]]; nam["dec_", i] -> MakeNode[
					(*scale up -> always, scale down -> never*)
					{If[ i ==1, False, scaling[[i]]], False}, 
					(*skip in -> accepts one skip, skip out -> never for dec*)
					{1, False}, 
					(*config*)
					{{blockType, setting[[i]]}, feature[[i]], {conf, ndim}}, 
					DropoutRate -> drop[[i]], RescaleMethod -> metSc, CatenateMethod -> metCat
				], {i, 1, depth - 1}],
				(*start and mapping*)
				{"start" -> UNetStart[nChan, fStart, dimIn, conf],
				"map" -> UNetMap[ndim, nClass]}
			],

			(*define the connections*)
			cons = Join[
				Flatten@Table[{
					(*connect encoding layers*)
					NetPort[nam["enc_", i], "Skip"] -> NetPort[nam["enc_", i + 1], "Down"],
					(*connect encoding to decoding layers*)
					NetPort[nam["enc_", i], "Skip"] -> NetPort[nam["dec_", i], "Skip"],
					(*connect decoding layers*)
					NetPort[If[i === depth - 1, nam["enc_", i + 1], nam["dec_", i + 1]], "Up"] -> NetPort[nam["dec_", i], "Scale"]
					
				}, {i, 1, depth - 1}],
				(*attach the start and map layers*)
				{"start" -> "enc_1", "dec_1" -> "map"}
			];
			debugUnet["The node connection list ", cons]; cons, 

			(*network settings and options*)
			"Input" -> Join[{nChan}, dimIn]
		];

		,
		"UNet+"|"UNet++",
		mapCon = If[mapCon === Automatic, depth -1, If[IntegerQ[mapCon] && mapCon<depth, mapCon, 1]];
		depthj = depth + 1;
		net = NetGraph[
			Join[
				(*Make all the nodes*)
				Flatten@Table[debugUnet[nam[i, j]]; nam[i, j] -> MakeNode[
					(*upscale for all nodes accept backbone -1, downscale only for backbone*)
					{If[j > 1, scaling[[i]], 1], If[j =!= 1 || i == depth, 1, scaling[[i]]]}, 
					(*skip in for all accept backbone, for UNET++ name the skips, skip out for all except right most upscale*)
					{If[j > 1, If[architecture==="UNet++", j - 1, 1], 0], If[j < depthj - i, True, False]}, 
					(*config*)
					{{blockType, setting[[i]]}, feature[[i]], {conf, ndim}}, 
					DropoutRate -> drop[[i]], RescaleMethod -> metSc, CatenateMethod -> metCat
				], {i, 1, depth}, {j, 1, depthj - i}],
				(*start and mapping*)
				{"start" -> UNetStart[nChan, fStart, dimIn, conf],
				"map" -> UNetMap[ndim, nClass, mapCon]}
			],
			cons = Join[
				Flatten@Table[{
					(*connect the backbone, the downscaling*)
					If[j === 1 && i=!= depth, NetPort[nam[i, j], If[boolSc[[i]], "Down", "Skip"]] -> nam[i + 1, j], Nothing],
					(*connect the nodes with up scaling*)
					If[1 < i <= depth, If[i=== depth, nam[i, j], NetPort[nam[i, j], If[j==depthj-i, "Up", "Skip"]]] -> NetPort[nam[i-1, j+1], "Scale"], Nothing],
					(*connect the node skip connection, for UNet++ its a dense connection.*)
					Switch[architecture,
						"UNet+",
						If[j < depthj - i && i=!= depth, NetPort[nam[i, j], "Skip"] -> NetPort[nam[i, j + 1], "Skip"], Nothing],
						"UNet++",
						If[j < depthj - i && i=!= depth, 
							NetPort[nam[i, j], "Skip"] -> Table[NetPort[nam[i, ji], If[ji===2, "Skip", nam["Skip",j]]], {ji, j + 1, depthj - i}]
						, Nothing]
					]
				}, {i, 1, depth}, {j, 1, depthj - i}],
				(*attach the start and map layers*)
				{"start" -> "node_1_1", Table[NetPort["node_1_"<>ToString[n], If[n=!=depth, "Skip", "Up"]], {n, depthj - mapCon, depth}] -> "map"}
			];
			debugUnet["The node connection list: ", cons]; 
			cons, 

			(*network settings and options*)
			"Input" -> Join[{nChan}, dimIn]		
		]
	];
	(*Monitor network properties*)
	mon[NetSummary[net], "Network description: "];

	(*return network*)
	net
]


(* ::Subsubsection::Closed:: *)
(*UNetStart*)


UNetStart[nChan_, feat_, dimIn_, conf_] := NetGraph[
	NetChain[Conv[feat , {Length@dimIn, 1}, (*conf*){"None", "None"}], 
		"Input" -> Prepend[dimIn, nChan]]
]


(* ::Subsubsection::Closed:: *)
(*ClassMap*)


UNetMap[dim_, nClass_] := UNetMap[dim, nClass, 1]

UNetMap[dim_, nClass_, n_] := Block[{map},
	map = Flatten[{Conv[nClass, {dim, 1}, {"None", "None"}], 
		If[nClass > 1, 
			{TransposeLayer[Which[dim === 2, {3, 1, 2}, dim === 3, {4, 1, 2, 3}]], SoftmaxLayer[]}, 
			{LogisticSigmoid, FlattenLayer[1]}
		]}];
	map = If[n>1, Prepend[map, CatenateLayer[]], map];
	NetGraph@NetChain@map
]


(* ::Subsubsection::Closed:: *)
(*MakeNode*)


Options[MakeNode] = {
	DropoutRate -> 0,
	RescaleMethod -> "Conv",
	CatenateMethod -> "Cat"
}

SyntaxInformation[MakeNode] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeNode[{scUp_ ,scDown_}, {skIn_, skOut_}, blockConfig_, OptionsPattern[]] := Block[{
		drop, metSc, metCat, chan, dim, block, node, boolCat, cat, catIn, conIn,
		boolScUp, scaleUp, blockUp, boolScDown, scaleDown, blockDown, boolDr, conOut, skip
	},
	(*get the node parameters*)
	{drop, metSc, metCat} = OptionValue[{DropoutRate, RescaleMethod, CatenateMethod}];
	chan = Last@Flatten@{blockConfig[[2]]};
	dim = blockConfig[[3, 2]];

	debugUnet[Column@{{scUp ,scDown}, {skIn, skOut}, {drop, metSc}, blockConfig}];

	(*define the block*)
	block = "block" -> ConvBlock @@ blockConfig;

	(*define the up scaling*)
	{boolScUp, scaleUp} = Which[
		scUp === False, {False, {1,1,1}},
		IntegerQ[scUp], {True, ConstantArray[scUp, dim]},
		ListQ[scUp], {True, PadRight[scUp, dim, 1]}
	];
	blockUp = Times @@ scaleUp =!= 1;
	scaleUp = If[boolScUp, "scaleU" -> ConvScale[{"Decode", metSc}, scaleUp, chan], Nothing];
	(*scale = If[boolScUp, NetPort["Scale"] -> "scaleU", NetPort["Scale"]];*)

	(*define the downscaling*)
	{boolScDown, scaleDown} = Which[
		scDown === False, {False, {1,1,1}},
		IntegerQ[scDown], {True, ConstantArray[scDown, dim]},
		ListQ[scDown], {True, PadRight[scDown, dim, 1]}
	];
	blockDown = (Times @@ scaleDown =!= 1) && metSc =!= "ConvS";
	scaleDown = If[boolScDown && blockDown, 
		"scaleD" -> ConvScale[{"Encode", metSc}, scaleDown, chan], 
		Nothing
	];

	(*define the dropout*)
	boolDr = 0 < drop;
	drop = If[boolDr, "drop" -> DropoutLayer[drop], Nothing];
	conOut = If[boolDr, "drop", "block"];

	(*define the skip connections*)
	skip = If[skOut=!=False, 
		conOut -> NetPort["Skip"<>If[IntegerQ[skOut], ToString[skOut], ""]], 
		Nothing
	];

	(*define the catenation*)
	boolCat = skIn > 0 (*&& boolScUp*);
	cat = If[boolCat, 
		"cat" -> If[(metSc === "Conv" || metSc === "ConvS") && metCat === "Tot", TotalLayer[], CatenateLayer[]], 
		Nothing
	];
	catIn = If[boolCat, 
		Append[If[skIn>1, Table[NetPort["Skip"<>ToString[i]], {i, skIn}], {NetPort["Skip"]}], 
			NetPort["Scale"]], 
		Nothing
	];
	conIn = If[boolCat, "cat", "block"];

	(*return the node*)
	node = NetGraph@NetFlatten[NetGraph[{block, drop, scaleUp, scaleDown, cat
		(*cat, block, scaleUp, scaleDown, drop*)}, {
			If[boolDr, "block" -> "drop", Nothing],
			If[boolScUp, If[blockUp, conOut -> "scaleU" -> NetPort["Up"], conOut -> NetPort["Up"]], Nothing], 
			If[boolScDown, If[blockDown, NetPort["Down"] -> "scaleD" -> conIn, NetPort["Down"] -> conIn], Nothing],
			If[boolCat, catIn -> "cat" -> "block", Nothing],
			skip
	}],1];

	debugUnet[node];
	node
]


(* ::Subsubsection::Closed:: *)
(*ConvScale*)


ConvScale[type_, scaleVec_, chan_] := Switch[type,
	{"Encode", "Pool"}, 
	{PoolingLayer[scaleVec, scaleVec]},
	{"Encode", "Conv"}, 
	{ConvolutionLayer[chan, 0 scaleVec + 3, "Stride" -> scaleVec,
		"PaddingSize"-> 1, "Dilation"-> 1]},
	{"Decode", "Pool"}, 
	{ResizeLayer[Scaled /@ scaleVec, Resampling -> "Nearest"]},
	{"Decode", "Conv"}, 
	{
		ResizeLayer[Scaled /@ scaleVec, Resampling -> "Nearest"], 
		ConvolutionLayer[chan, 0 scaleVec + 3, "Stride" -> 1, 
			"PaddingSize" -> 1, "Dilation" -> 1]
	},
	{"Decode", "ConvS"}, 
	{
		ConvolutionLayer[chan, 0 scaleVec + 1, "Stride" -> 1, 
			"PaddingSize" -> 0, "Dilation" -> 1],
		ResizeLayer[Scaled /@ scaleVec, Resampling -> "Linear"]
	}
]


(* ::Subsubsection::Closed:: *)
(*ConvBlock*)


ConvBlock[block_, feat_?IntegerQ, dim_?IntegerQ] := ConvBlock[block, {feat, Round[feat/2]}, {{"None", "Instance"}, dim}]

ConvBlock[block_, feat_?IntegerQ, {act_, dim_}] := ConvBlock[block, {feat, Round[feat/2]}, {act, dim}]

ConvBlock[block_, {featOut_, featInt_}, {act_, dim_}, scale_:1] := Block[{
		blockType, type, blockSet, repBlock, dep, rep, nam, scaleF, cons
	},

	(*short notation for naming layers*)
	nam = #1 <> ToString[#2] &;
	(*get the block settings if not defined*)
	{blockType, blockSet} = If[Length[block] === 2, block, {block, block /. netDefaults}];

	(*switch between the different block types*)
	Switch[blockType,
		"Conv", 
		dep = blockSet;
		ConvBlock[{"UNet", blockSet}, {featOut, featInt}, {act, dim}],

		"UNet",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"conv" -> Flatten[Table[Conv[featOut, {dim, 3}, act], {i, 1, dep}]]}, {}
		],

		"ResNet",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"con" -> Flatten[Table[
				If[i =!= dep, Conv[Round[featOut/2], {dim, 3}, act], Conv[featOut, {dim, 3}, act]], {i, 1, dep}]],
			"skip" -> Conv[featOut, {dim, 1}, {"None", Last@Flatten[{act}]}],
			"tot" -> {TotalLayer[]}
			}, {{"con", "skip"} -> "tot"}
		],

		"ResNetS",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"con" -> Flatten[Table[Conv[If[i =!= dep, Round[featOut/2], featOut], {dim, 3}, act, "Stride"-> If[i===1, scale, 1]], {i, 1, dep}]],
			"skip" -> Conv[featOut, {dim, scale}, {"None", Last@Flatten[{act}]}, "Stride" -> scale],
			"tot" -> {TotalLayer[]}
			}, {{"con", "skip"} -> "tot"}
		],

		"ResNetL",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"con" -> Flatten[Table[If[i =!= dep, Conv[featInt, {dim, 1}, act], Conv[featOut, {dim, 3}, {"None", Last@Flatten[{act}]}]], {i, 1, dep}]],
			"skip" -> Conv[featOut, {dim, 1}, {"None", Last@Flatten[{act}]}],
			"tot" -> {TotalLayer[], ActivationLayer[First[Flatten[{act}]]]}
			}, {{"con", "skip"} -> "tot"}
		],

		"DenseNet",
		{dep, rep} = If[IntegerQ[blockSet], {blockSet, 2}, blockSet];
		repBlock = NetGraph[
			Table[nam["bl", i] -> (Prepend[Conv[If[i === dep, featOut, featInt], {dim, 3}, act], If[i == 1, Nothing, CatenateLayer[]]]), {i, 1, dep}],
			Table[Table[If[j == 0, NetPort["Input"], nam["bl", j]], {j, 0, i - 1}] -> nam["bl", i], {i, 1, dep}]
		];
		NetGraph@NetFlatten[NetGraph[
			Table[nam["rep", i] -> repBlock, {i, rep}],
			Table[nam["rep", i] -> nam["rep", i + 1], {i, rep - 1}]
		], 2],

		"Inception",
		{dep, rep} = If[IntegerQ[blockSet], {blockSet, 2}, blockSet];
		repBlock = NetGraph[Join[
				Table[nam["dil", 2 (i - 1) + 1] -> Conv[featInt, {dim, 3}, act, "Dilation"->2 (i - 1) + 1], {i, 1, dep}],
				{"cat" -> CatenateLayer[]}
			], {NetPort["Input"] -> Table[nam["dil", 2 (i - 1) + 1], {i, 1, dep}] -> "cat"}
		];
		NetGraph@NetFlatten[NetGraph[
			Append[Table[nam["rep", i] -> repBlock, {i, rep}], "out" -> Conv[featOut, {dim, 1}, act]],
			Append[Table[nam["rep", i] -> nam["rep", i + 1], {i, rep - 1}], 
			nam["rep", rep] -> "out"]
		], 2],

		"U2Net",
		{dep, type} = If[IntegerQ[blockSet], {blockSet, True}, blockSet];
		scaleF = If[! #1, 1, 2] &;
		NetGraph[
			Join[
				Table[debugUnet[nam["U2enc_", i]]; nam["U2enc_", i] -> MakeNode[
					(*scale up -> never, scale down*)
					{1, If[i==dep, 1, scaleF[type]]}, 
					(*skip in -> only input from above, skip out -> always for enc*)
					{0, True}, 
					(*config*)
					{"Conv", featInt, {act, dim}}
				], {i, 1, dep}],
				Table[debugUnet[nam["U2dec_", i]]; nam["U2dec_", i] -> MakeNode[
					(*scale up -> always, scale down -> never*)
					{scaleF[type], 1}, 
					(*skip in -> accepts one skip, skip out -> never for dec*)
					{1, False}, 
					(*config*)
					{"Conv", If[i === 1, featOut, featInt], {act, dim}}
				], {i, 1, dep - 1}],
				{"start" -> Conv[featOut, {dim, 1}, act], "add" -> TotalLayer[]}
			], 
			cons = Join[
				Flatten@Table[{
					NetPort[nam["U2enc_", i], If[type, "Down", "Skip"]] -> nam["U2enc_", i + 1], 
					NetPort[nam["U2enc_", i], "Skip"] -> NetPort[nam["U2dec_", i], "Skip"],
					If[i === dep - 1, nam["U2enc_", i + 1], nam["U2dec_", i + 1]] -> NetPort[nam["U2dec_", i], "Scale"]}
				, {i, 1, dep - 1}],
				{"start" -> "U2enc_1", {"start", "U2dec_1"} -> "add"}
			];
			debugUnet[cons]; cons
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*Conv*)

Options[Conv] = {"Dilation" -> 1, "Stride"->1, "ChannelGroups"->1}

Conv[featOut_?IntegerQ, {dim_, kern_}, opts:OptionsPattern[]] := Conv[featOut, {dim, kern}, {"None", "Instance", opts}]

Conv[featOut_?IntegerQ, {dim_, kern_}, act_?StringQ, opts:OptionsPattern[]] := Conv[featOut, {dim, kern}, {act, "Instance"}, opts]

Conv[featOut_?IntegerQ, {dim_, kern_}, {act_?StringQ, bat_?StringQ}, OptionsPattern[]] := Block[{
		dil, str
	},

	dil = OptionValue["Dilation"];
	str = OptionValue["Stride"];
	gr = OptionValue["ChannelGroups"];
	ker = If[Length[kern]===dim, kern, ConstantArray[kern, dim]];

	{
		(*The most basic conv block CON>BN>ACT used in each of the advanced conv blocks*)
		ConvolutionLayer[featOut, ker, 
			"PaddingSize" -> (Ceiling[ker/2] - 1) dil, 
			"Stride" -> str, "Dilation" -> dil, "ChannelGroups" -> gr],	
		Switch[bat, 
			"Group", NormalizationLayer[All, 1, "GroupNumber" -> {8, 1, 1, 1}], 
			"Instance", NormalizationLayer[],
			"None", Nothing,
			_, BatchNormalizationLayer[]
		], 
		ActivationLayer[act]
	}
]


(* ::Subsubsection::Closed:: *)
(*ActivationLayer*)


ActivationLayer[] := ActivationLayer["GELU"]

ActivationLayer[actType_] := If[Head[actType]===ParametricRampLayer||Head[actType]===ElementwiseLayer,
	actType,
	Switch[actType, 
		"LeakyRELU", ParametricRampLayer[], 
		"None" | "", Nothing,
		_, If[StringQ[actType],ElementwiseLayer[actType], Message[ActivationLayer::type, actType];Nothing]
	]
]


(* ::Subsection:: *)
(*NetDimensions*)


(* ::Subsubsection::Closed:: *)
(*NetDimensions*)


NetDimensions[net_] := NetDimensions[net, ""]

NetDimensions[net_, port_] := Block[{block, neti},Switch[port,
	"Input", 
	Information[net,"InputPorts"]["Input"],

	"Output", 
	Information[net,"OutputPorts"]["Output"],

	"FirstEncodingIn",
	Last@Values@Information[
		NetTake[net, {block = First[Select[Keys[net[[All, 1]]], StringContainsQ[#, "enc_"| ("node_" ~~ __ ~~ "_1")] &]], 
			block}], "InputPorts"],

	"FirstEncodingOut", 
	Values@Information[
		NetTake[net, {block = First[Select[Keys[net[[All, 1]]], StringContainsQ[#, "enc_"| ("node_" ~~ __ ~~ "_1")] &]], 
			block}], "OutputPorts"],

	"LastEncodingIn", 
	Last@Values@Information[
		NetTake[net, {block = Last[Select[Keys[net[[All, 1]]], StringContainsQ[#, "enc_"| ("node_" ~~ __ ~~ "_1")] &]], 
			block}], "InputPorts"],

	"LastEncodingOut", 
	Last@Values@Information[
		NetTake[net, {block = Last[Select[Keys[net[[All, 1]]], StringContainsQ[#, "enc_"| ("node_" ~~ __ ~~ "_1")] &]], 
			block}], "OutputPorts"],

	"MinEncodingOut",
	Min /@ Transpose[Flatten[Values[Information[#, "OutputPorts"]] & /@ Information[net, "LayersList"], 1]],

	"AllEncodingOut",
	Max[Values[Information[NetTake[net, {#}], "OutputPorts"]][[All, 1]]] & /@ Keys[net[[All, 1]]],

	"AllMaxChannels",
	Max[Values[Information[#, "OutputPorts"] & /@ Select[Information[NetTake[net, {#, #}], "LayersList"], Head[#] === ConvolutionLayer &]][[All, 1, 1]]] & /@ Keys[net[[All, 1]]],

	"U2Encoding",
	block = Select[Keys[net[[All, 1]]], StringContainsQ[#, "enc_"|("node_" ~~ __ ~~ "_1")] &];
	If[block==={}, Return[$Failed],
		block = First[block];
		neti = NetFlatten[NetTake[net, {block, block}], 1];
		block = Select[Keys[Normal[neti]], StringContainsQ[#, "block/U2enc_"] &];
		If[block === {}, Return[$Failed],
			block = Last[block];
			Last@Values@Information[NetTake[neti, {block, block}], "OutputPorts"]
	]],

	_, 
	{First@NetDimensions[net, "Input"], Last@NetDimensions[net, "Output"], 
		Rest@NetDimensions[net, "Input"], First@NetDimensions[net, "FirstEncodingIn"]}
]]


(* ::Subsubsection::Closed:: *)
(*ChangeNetDimensions*)


Options[ChangeNetDimensions] = {
	"Dimensions" -> None,
	"Channels" -> None,
	"Classes" -> None
};

SyntaxInformation[ChangeNetDimensions] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ChangeNetDimensions[netIn_, OptionsPattern[]] := Block[{
		dimIn, nChanIn, nClassIn,
		dimOut, nChanOut, nClassOut, filt, netOut, start
	},

	{dimOut, nChanOut, nClassOut} = OptionValue[{"Dimensions", "Channels", "Classes"}];

	(*Figure out all the dimensions*)
	netOut = netIn;
	dimIn = NetDimensions[netIn, "Input"];
	nChanIn = First@dimIn;
	dimIn = Rest@dimIn;
	nClassIn = Last@NetDimensions[netIn, "Output"];

	(*Change network Dimensions if needed*)
	If[dimOut =!= None && dimOut =!= dimIn, 
		netOut = NetReplacePart[netOut, "Input" -> Round[Prepend[dimOut, nChanIn]]],
		dimOut = dimIn
	];

	(*Change input Channels if needed*)
	If[nChanOut =!= None && nChanOut=!=nChanIn,		
		netOut = NetReplacePart[netOut, {
			"Input" -> Prepend[Round@dimOut, Round@nChanOut],
			"start" -> UNetStart[Round@nChanOut, First@NetDimensions[netOut, "FirstEncodingIn"], 
				Round@dimOut, NetFlatten[NetTake[netOut, "start"]][[-1]]]
		}]
	];

	(*Change output Classes if needed*)
	If[nClassOut =!= None && nClassOut=!=nClassIn,
		netOut = NetReplacePart[netOut, {
			"map" -> UNetMap[Length@dimIn, Round@nClassOut], 
			"Output" -> Round@Append[dimOut, nClassOut]
		}]
	];

	netOut
]


(* ::Subsubsection::Closed:: *)
(*GetNetNodes*)


SyntaxInformation[GetNetNodes] = {"ArgumentsPattern" -> {_}};

GetNetNodes[net_] := DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]]


(* ::Subsection:: *)
(*LossLayers*)


(* ::Subsubsection::Closed:: *)
(*AddLossLayer*)


SyntaxInformation[AddLossLayer] = {"ArgumentsPattern" -> {_}};

(*http://arxiv.org/abs/2312.05391*)
AddLossLayer[net_] := NetGraph[<|
	"net"->net,
	"Dice" -> NetFlatten@NetGraph@NetChain@{DiceLossLayer[2]}, 
	"Jaccard" -> NetFlatten@NetGraph@NetChain@{JaccardLossLayer[2]},
	"Tversky" -> NetFlatten@NetGraph@NetChain@{TverskyLossLayer[0.7]}, 
	"MSD" -> NetFlatten@NetGraph@NetChain@{MeanSquaredLossLayer[], ElementwiseLayer[50 #&]},
	"Focal" -> NetFlatten@NetGraph@NetChain@{FocalLossLayer[2, 0.25], ElementwiseLayer[5 #&]},
	"TopK" -> NetFlatten@NetGraph@NetChain@{TopKLossLayer[net, 0.1]},
	"CE" -> NetFlatten@NetGraph@NetChain@{CrossEntropyLossLayer["Probabilities"]}
|>, {
	{"net", NetPort["Target"]} -> "Dice" -> NetPort["Dice"], (*using squared dice, F1score*)
	{"net", NetPort["Target"]} -> "Jaccard" -> NetPort["Jaccard"],(*using squared Intersection over union*)
	{"net", NetPort["Target"]} -> "Tversky" -> NetPort["Tversky"], (*recall more than precision, 0.5 is dice*)
	{"net", NetPort["Target"]} -> "MSD" -> NetPort["MSD"], (*Brier Score*)
	{"net", NetPort["Target"]} -> "Focal" -> NetPort["Focal"], (*scaled cross entropy for hard examples*)
	{"net", NetPort["Target"]} -> "TopK" -> NetPort["TopK"], (*scaled cross entropy for hard examples*)
	{"net", NetPort["Target"]} -> "CE" -> NetPort["CE"]
}]


(* ::Subsubsection::Closed:: *)
(*DiceLossLayer*)


SyntaxInformation[DiceLossLayer] = {"ArgumentsPattern" -> {_., _.}};

DiceLossLayer[] := DiceLossLayer[]

DiceLossLayer[n_?IntegerQ] := NetFlatten[
	(*10.48550/arXiv.1911.02855 and https://doi.org/10.48550/arXiv.1606.04797 for squared dice loss look at v-net*)
	(*https://arxiv.org/abs/1707.03237 for the generalized with class weighting*)
	NetGraph[<|
		(*flatten input and target; function layer allows to switch to L2 norm if n = 2*)
		"input" -> {If[n > 1, FunctionLayer[#^n &], Nothing], AggregationLayer[Total, ;; -2]},
		"target" -> {AggregationLayer[Total, ;; -2]},
		(*intersection or TP*)
		"intersection" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; -2]},
		(*2*intersection / (input + target), +1 is for numerical stability*)
		"dice" -> {ThreadingLayer[1. - (2. #1 + 1.) / (#2 + #3 + 1.) &], AggregationLayer[Mean, 1]}
	|>, {
		{NetPort["Target"], NetPort["Input"]} -> "intersection",
		NetPort["Input"] -> "input",
		NetPort["Target"] -> "target",
		{"intersection", "target", "input"} -> "dice" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*JaccardLossLayer*)


SyntaxInformation[JaccardLossLayer] = {"ArgumentsPattern" -> {_., _.}};

JaccardLossLayer[] := JaccardLossLayer[2]

JaccardLossLayer[n_?IntegerQ] := NetFlatten[
	(*https://arxiv.org/pdf/1906.11600 for the definition of squared jaccard loss*)
	(*https://arxiv.org/html/2302.05666v5 soft jaccard*)
	NetGraph[<|
		(*flatten input and target; function layer allows to switch to L2 norm if n = 1*)
		"input" -> {If[n > 1, FunctionLayer[#^n &], Nothing], AggregationLayer[Total, ;; -2]},
		"target" -> {AggregationLayer[Total, ;; -2]},
		(*intersection or TP*)
		"intersection" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; -2]},
		(*intersection / union with union = (input + target - intersection), +1 is for numerical stability*)
		"Jaccard" -> {ThreadingLayer[1. - (#1 + 1.) / (#3 + #2 - #1 + 1.) &], AggregationLayer[Mean, 1]}
	|>, {
		NetPort["Input"] -> "input",
		NetPort["Target"] -> "target",
		{NetPort["Target"], NetPort["Input"]} -> "intersection",
		{"intersection", "target", "input"} -> "Jaccard" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*TverskyLossLayer*)


SyntaxInformation[TverskyLossLayer] = {"ArgumentsPattern" -> {_.}};

TverskyLossLayer[] := TverskyLossLayer[0.7]

TverskyLossLayer[beta_?NumberQ] := TverskyLossLayer[1- beta, beta]

TverskyLossLayer[alpha_?NumberQ, beta_?NumberQ] := NetFlatten[
	(*https://doi.org/10.48550/arXiv.1706.05721 Tversky loss function for 3D segmentation*)
	(*generalization of dice and jaccard, alpha = beta = 0.5 means dice, alpha = beta = 1 means jaccard*)
	NetGraph[<|
		(*intersection or TP, first input is Target, second is Input*)
		"truePos" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; -2]},
		"falsePos" -> {ThreadingLayer[(1 - #1) #2 &], AggregationLayer[Total, ;; -2]},
		"falseNeg" -> {ThreadingLayer[#1 (1 - #2) &], AggregationLayer[Total, ;; -2]},
		(*TP / (TP + a FP + b FN), +1 is for numerical stability, alpha = 1 - beta *)
		"Tversky" -> {ThreadingLayer[1. - (#1 + 1) / (#1 + alpha #2 + beta #3 + 1) &], 
			AggregationLayer[Mean, 1]}
	|>, {
		{NetPort["Target"], NetPort["Input"]} -> "truePos",
		{NetPort["Target"], NetPort["Input"]} -> "falsePos",
		{NetPort["Target"], NetPort["Input"]} -> "falseNeg",
		{"truePos", "falsePos", "falseNeg"} -> "Tversky" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


SyntaxInformation[
OverlapLossFunction] = {"ArgumentsPattern" -> {_., _.}};

OverlapLossFunction[] := OverlapLossFunction[0.7, 1]

OverlapLossFunction[beta_?NumberQ, 1] := OverlapLossFunction[{1 - beta, beta}, 1]

OverlapLossFunction[beta_?NumberQ, n_?IntegerQ] := OverlapLossFunction[{1 - beta, beta}, n]

OverlapLossFunction[{alpha_?NumberQ, beta_?NumberQ}] := OverlapLossFunction[{alpha, beta}, 1]

OverlapLossFunction[{alpha_?NumberQ, beta_?NumberQ}, n_?IntegerQ] := NetFlatten[
(*https://doi.org/10.48550/ arXiv.1706.05721 Tversky loss function for 3D segmentation*)
(*https://arxiv.org/pdf/1906.11600 for the definition of squared jaccard loss*)
(*https://arxiv.org/abs/ 1707.03237 for the generalized with class weighting*)
(*10.48550/arXiv.1911.02855 and 10.48550/arXiv.1606.04797 for squared dice loss look at v-net*)
(*https://arxiv.org/abs/1707.03237 for the generalized with class weighting*)
(*tversky generalization of dice and jaccard,alpha=beta=0.5 means dice,alpha=beta=1 means jaccard*)
	NetGraph[<|
	If[n > 1, "input" -> FunctionLayer[#^n &], Nothing],
	(*intersection or TP,first input is Target,second is Input*)
	"TP" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; -2]},
	"TPn" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; -2]},
	"FPn" -> {ThreadingLayer[(1. - #1) #2 &], AggregationLayer[Total, ;; -2]},
	"FNn" -> {ThreadingLayer[#1 (1. - #2) &], AggregationLayer[Total, ;; -2]},
	(*the loss function TP / (TP + a FP + b FN), + 1 is for numerical stability, form of laplace smoothing*)
	"loss" -> {ThreadingLayer[1. - (#1 + 1)/(#2 + alpha #3 + beta #4 + 1) &], AggregationLayer[Mean, 1]}
	(*"weight"->If[w>0.,{FunctionLayer[If[#<1.,#,0.]&/@(1./(#^w+1.))&],FunctionLayer[#/Mean[#]&]},{FunctionLayer[0. #+1.&]}]*)
	|>, {
		{NetPort["Target"], NetPort["Input"]} -> "TP",
		If[n > 1, NetPort["Input"] -> "input", Nothing], 
		{NetPort["Target"], If[n > 1, "input", NetPort["Input"]]} -> {"TPn", "FPn", "FNn"},
		{"TP", "TPn", "FPn", "FNn"} -> "loss" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*FocalLossLayer*)


SyntaxInformation[FocalLossLayer] = {"ArgumentsPattern" -> {_., _.}};

FocalLossLayer[] := FocalLossLayer[2, 0.25]

FocalLossLayer[g_] := FocalLossLayer[g, 0.25]

FocalLossLayer[g_, a_] := NetFlatten[
	(*https://arxiv.org/abs/1708.02002v2 for definition of focal loss, 10^-15 is for numerical stability*)
	NetGraph[{
		"flatPr" -> {ThreadingLayer[#1 #2 &], AggregationLayer[Total, {-1}], FlattenLayer[]},
		"focal" -> {ThreadingLayer[-a Log[#1 + 10^-15](1 - #1)^g &], AggregationLayer[Mean, 1]}
		}, {
			{NetPort["Input"], NetPort["Target"]} -> "flatPr" -> "focal" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*TopKLossLayer*)


SyntaxInformation[TopKLossLayer] = {"ArgumentsPattern" -> {_.}};

TopKLossLayer[net_] := TopKLossLayer[net, 0.1]

TopKLossLayer[net_, k_] := Block[{
		i = Round[k Times @@ NetDimensions[net, "Input"]]
	},
	NetFlatten[
		(*10.48550/arXiv.1512.00486 for definition of top k losses, 10^-15 is for numerical stability*)
		NetGraph[{
			"flatPr" -> {ThreadingLayer[#1 #2 &], AggregationLayer[Total, {-1}], FlattenLayer[]},
			"topK" -> {NetChain[{FunctionLayer[Sort[-Log[# + 10^-15]] &],PartLayer[-i;;-1], AggregationLayer[Mean, 1]}]}
			}, {
				{NetPort["Input"], NetPort["Target"]} -> "flatPr" -> "topK" -> NetPort["Loss"]
		}, "Loss" -> "Real"]
	]
]


(* ::Subsection:: *)
(*Encoders*)


(* ::Subsubsection::Closed:: *)
(*ClassEncoder*)


SyntaxInformation[ClassEncoder] = {"ArgumentsPattern" -> {_, _.}};

ClassEncoder[data_] := ClassEncoder[data, Round[Max@data]+1]

ClassEncoder[data_, nClass_] := ToPackedArray@Round@If[nClass === 1, data, ClassEncoderC[data, UnitVector[nClass, 1]]]

ClassEncoderC = Compile[{{class, _Integer, 0}, {vec, _Integer, 1}}, 
	RotateRight[vec, class], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*ClassDecoder*)


SyntaxInformation[ClassDecoder] = {"ArgumentsPattern" -> {_, _.}};

ClassDecoder[data_] := ClassDecoder[data, Last@Dimensions@data]

ClassDecoder[data_, nClass_] := ToPackedArray@Round@ClassDecoderC[data, Range[nClass]]

ClassDecoderC = Compile[{{prob, _Real, 1}, {classes, _Integer, 1}}, 
	Max[classes (1 - Unitize[Chop[(prob/Max[prob]) - 1]])] - 1, 
RuntimeAttributes -> {Listable}]


(* ::Subsection::Closed:: *)
(*NetSummary*)


NetSummary[net_] := NetSummary[net, ""]

NetSummary[net_?ListQ, rep_?StringQ] := NetSummary[#, rep]&/@net

NetSummary[net_, rep_?StringQ] := Block[{
		toK, st, quantStr, quantStrG, lays, convs, kerns, count, nKern, kernWeights,
		norm, normWeights, nelem, elems, elemSize, arrSize, netSize, table, nodes, netIm, nodeIm
	},

	toK = Which[
		# > 1000000, ToString[NumberForm[#/1000000., {Infinity, 1}]] <> " M",
		# > 1000, ToString[NumberForm[Round[#/1000.], {Infinity, 0}]] <> " K",
		True, #] &;
	st = Style[#1, Bold] & ;
	quantStr = ToString[Round[QuantityMagnitude[#], .1]] <> " " <> (QuantityUnit[#] /. {"Megabytes" -> "MB", "Gigabytes" -> "GB"}) &;
	quantStrG = ToString[Round[QuantityMagnitude[UnitConvert[#, "Gigabytes"]], .1]] <> " GB" &;

	lays = Information[net, "LayersList"];
	convs = Select[lays, Head[#] === ConvolutionLayer &];
	kerns = Information[#, "ArraysDimensions"][{"Weights"}] & /@ convs;
	kerns = GatherBy[{#[[3 ;;]], Times @@ #[[1 ;; 2]], Times @@ #[[1 ;;]]} & /@	kerns, First];

	count = Sort[{Length[#], #[[1, 1]], Total[#[[All, 2]]], Total[#[[All, 3]]]} & /@ kerns];
	nKern = Total[count[[All, 3]]];
	kernWeights = Total[count[[All, 4]]];
	count[[All, 3]] = toK /@ count[[All, 3]];
	count[[All, 4]] = toK /@ count[[All, 4]];

	norm = Select[lays, (Head[#] === BatchNormalizationLayer || Head[#] === NormalizationLayer )&];
	normWeights = If[norm==={}, 0, First@Total[Total[Information[#, "ArraysDimensions"] /@ Keys[Information[norm[[1]], "ArraysDimensions"]]] & /@ norm]];

	nelem = Information[net, "ArraysTotalElementCount"];
	elems = Round[nelem/1000000, .01];
	elemSize = UnitConvert[Quantity[32. nelem, "Bits"], "MB"];
	arrSize = UnitConvert[Quantity[32. Total[Times @@@ Values[(Information[#, "OutputPorts"]["Output"]) & /@ Information[net, "Layers"]]], "Bits"], "MB"];
	netSize = UnitConvert[elemSize + arrSize, "GB"];

	table = Grid[{{ Grid[{
			{st@"Number of Norm. Layers: ", st@Length@norm},
			{st@" - Number of Weights: ", toK@normWeights},
			{st@"Number of Conv. Layers: ", st@Length@convs},
			{st@" - Number of Kernels: ", toK@nKern},
			{st@" - Number of Weights: ", toK@kernWeights},
			{""},
			{st@"Convolution Kernel Distribution:", SpanFromLeft},
			{Item[Grid[Join[{Style[#, Bold] & /@ {"Count", "Size", "Kernels", "Weights"}}, count], 
				Alignment -> Right, Spacings -> {1.5, 1}], Alignment->Right], SpanFromLeft},
			{""},
			{st@"Estimated Memory:", SpanFromLeft},
			{st@" - Total Weight Memory", quantStr@elemSize},
			{st@" - Total Network Memory", quantStr@netSize},
			{st@" - Estimaged Train Memory", "n * " <> quantStrG[2 * (arrSize + elemSize)]}
		}, Alignment -> {{Left, Right}}, Spacings -> {1, 1}(*, Background -> LightDarkV[GrayLevel[.95], GrayLevel[.2]]]*)]
	}}, Frame -> All, Spacings -> {2, 2}, Background -> LightDarkV[GrayLevel[.95], GrayLevel[.2]]];

	Switch[rep,
		"Full",
		makeNetIm = With[{im = Information[#, "SummaryGraphic"]}, 
			Show[im, AspectRatio -> 0.5, ImageSize -> Max[AbsoluteOptions[im, ImageSize][[1, 2]]]]
		]&;
		nodes = DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]];
		netIm = makeNetIm @ net;
		nodeIm = Information[NetTake[net, {#, #}][[1]], "SummaryGraphic"] & /@ nodes;
		TabView[Join[{"summary" -> table, "net" -> netIm}, Thread[nodes -> nodeIm]],
			ControlPlacement -> Left, Alignment -> Center],

		"Mem",
		Grid[{
			{st["Weight mem: "], quantStr[elemSize]}, 
			{st["Net mem: "], quantStr[netSize]},
			{st["Train mem: "], "n * " <> quantStrG[2 * (arrSize + elemSize)]}
		}, Alignment -> Left],
		"MemVal",
		netSize,
		_, 
		table
	]
]


(* ::Subsection::Closed:: *)
(*AnalyzeNetworkFeatures*)


AnalyzeNetworkFeatures[net_, datI_] := AnalyzeNetworkFeatures[net, datI, ""]

AnalyzeNetworkFeatures[net_, datI_, met_] := Block[{
		data, dim, dataP, netP, nodes, vals, cutoff, table, plot, feat, nfeat, ttt, n, col
	},

	data = datI;
	If[ArrayDepth[data] === 4, data = data[[All, 1]]];

	(*find the patch dimensions and adjust data and network*)
	dim = QMRITools`SegmentationTools`Private`FindPatchDim[net, Dimensions@data];
	dataP = NormalizeData[PadToDimensions[data, dim], NormalizeMethod -> "Uniform"];
	netP = ChangeNetDimensions[net, "Dimensions" -> dim];

	(*extract the network nodes*)
	nodes = GetNetNodes[netP];
	nodes = nodes[[2;;-2]];

	col = ColorData["SolarColors"] /@Rescale[Range[Ceiling[Length[nodes]/2]]];
	col = Join[col , Reverse[col][[2;;]]];

	(*calculate the singular values for plotting and reporting*)
	{nfeat, table, plot} = Transpose[(
		(*get the features*)
		feat = NetTake[netP, #][{dataP}, TargetDevice -> "CPU"];
		If[Head[feat] === Association, feat = Last@feat];
		feat = Map[Flatten, feat];

		(*calculate the singular values of the features*)
		nfeat = Length@feat;
		vals = Diagonal[SingularValueDecomposition[feat, UpTo[nfeat]][[2]]];
		vals = 100 Rescale[Accumulate[vals]];

		(*find cuoff index and percentage*)
		cutoff = First[Position[UnitStep[vals - 99], 1]] - 1;

		(*give the output*)
		{
			nfeat, 
			Flatten[{Style[#, 14, Bold], cutoff, Round[100 cutoff/nfeat, .1]}], 
			Transpose[{100 Rescale[Range[1., nfeat]], vals}]
		}
	) & /@ nodes];

	(*find the nodes with the highest singular values*)
	ttt = table[[All, 2]];
	n = Ceiling[(Length@ttt)/2];
	n = Max /@ Thread[{ttt[[ ;; n]], Reverse[ttt[[n ;; ]]]}];

	(*output based on method*)
	If[met === "",
		MonitorFunction[Thread[{nodes, nfeat}]];
		MonitorFunction[n];

		(*dynamic plot output*)
		DynamicModule[{cols, tab = table, pl = plot, nods = nodes, clist = col, ln, pcol},

			(*define colors for plotting*)
			cols = Table[Directive[{clist[[i]], Dashed, Thick}], {i, Length@nodes}];
			ln = Range@Length@nods;

			(*define the plots within a manipulate that allows to select the nodes of the network*)
			Manipulate[
				pcol = cols;
				pcol[[k]] = Directive[{(*Black*)clist[[k]], Dashing[None], Thickness[.01]}];

				Column[{
				Grid[Transpose@tab, Frame -> All, Background -> {{k -> Gray}, None, 
					Thread[Thread[{1, ln}] -> (Lighter/@clist)]}, Spacings -> {1.2, 1.2}],
				Show[
					ListLinePlot[pl, PlotStyle -> pcol(*RotateRight[cols, k - 1]*), GridLines -> {{tab[[k, 3]]}, {99}}, 
						ImageSize -> 500, AspectRatio -> 1, $plotOptions
					],
					Plot[x, {x, 0, 100}, PlotStyle -> Directive[{Thick, Gray, Dotted}]]
				]
			}, Alignment -> Center],
			{{k, 1, ""}, Thread[ln -> (Style[#, LightDarkV[Black, White], 14, Bold] & /@ nods)], ControlType -> SetterBar}
			]
		],

		(*value per node*)
		table[[All, 2]]
	]
]


(* ::Subsection:: *)
(*MakeClassify*)


(* ::Subsubsection::Closed:: *)
(*MakeClassifyNetwork*)


Options[MakeClassifyNetwork] = {
	ImageSize -> {128, 128}
};

SyntaxInformation[MakeClassifyNetwork] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

MakeClassifyNetwork[classes_, OptionsPattern[]] := Block[{imSize, enc, dec, conv, head, heads, nodes, connection},
	
	(*make the image enoder*)
	imSize = OptionValue[ImageSize];
	enc = NetEncoder[{"Image", imSize, ColorSpace -> "Grayscale"}];

	(*make the decoders: switch between single list of classes or named association of classes*)
	dec = If[ListQ[classes],
		Thread[{"Output"} -> {NetDecoder[{"Class", classes}]}]		,
		Thread[Keys[classes] -> (NetDecoder[{"Class", #}] & /@ Values[classes])]
	];

	(*general convolution layer*)
	conv = NetChain[{
		ConvolutionLayer[16, 7, "Stride" -> 1, PaddingSize -> 3], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[2, 2],
		ConvolutionLayer[32, 5, "Stride" -> 1, PaddingSize -> 2], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[2, 2],
		ConvolutionLayer[32, 5, "Stride" -> 1, PaddingSize -> 2], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[2, 2],
		ConvolutionLayer[64, 3, "Stride" -> 1, PaddingSize -> 1], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[2, 2],
		ConvolutionLayer[64, 3, "Stride" -> 1, PaddingSize -> 1], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[4, 4], 
		FlattenLayer[]
	}];

	(*class specific head function*)
	head = NetChain[{
		LinearLayer[256], BatchNormalizationLayer[], ElementwiseLayer["GELU"],
		LinearLayer[128], BatchNormalizationLayer[], ElementwiseLayer["GELU"],
		LinearLayer[64], BatchNormalizationLayer[], ElementwiseLayer["GELU"],
		LinearLayer[32], BatchNormalizationLayer[], ElementwiseLayer["GELU"],
		LinearLayer[Length[#]], SoftmaxLayer[]}
	] &;

	(*make the heads and the convolution layers and the connections*)
	;
	
	(*make the network*)
	NetGraph[
		Association[Join[{"Conv" -> conv},  Thread[Keys[classes] -> (head /@ Values[classes])]]], 
		("Conv" -> # -> NetPort[#] & /@ Keys[classes]),
		"Input" -> enc, ##
	] & @@ dec
]


(* ::Subsubsection::Closed:: *)
(*MakeClassifyImage*)


Options[MakeClassifyImage]={
	ImageSize -> {128, 128}
};

SyntaxInformation[MakeClassifyImage] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

MakeClassifyImage[dat_, opts:OptionsPattern[]] := Switch[ArrayDepth[dat],
	2, MakeClassifyImage[dat, opts],
	3, MakeClassifyImage[#, opts]&/@First[AutoCropData[dat]],
	4, MakeClassifyImage[#, opts]&/@First[AutoCropData[dat[[All, 1]]]],
	_, $Failed
]

MakeClassifyImage[dat_?MatrixQ, OptionsPattern[]] := Block[{imSize},
	imSize = OptionValue[ImageSize];
	If[Total[Flatten[dat]]<10,
		Image@ConstantArray[0., imSize],
		ImageResize[Image[Rescale[First[AutoCropData[{dat}, CropPadding -> 5][[1]]]]], imSize]
	]
];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
