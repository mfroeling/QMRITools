(* ::Package:: *)

(* ::Title:: *)
(*QMRITools SegmentationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`SegmentationTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`SegmentationTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


MakeUnet::usage = 
"MakeUnet[nChannels, nClasses, dep, dimIn] Generates a UNET with nChannels as input and nClasses as output. 
he number of parameter of the first convolution layer can be set with dep.\n
The data dimensions can be 2D or 3D and each of the dimensions should be 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240 or 256."

AddLossLayer::usage = 
"AddLossLayer[net] adds three loss layers to a NetGraph, a SoftDiceLossLayer, BrierLossLayer and a CrossEntropyLossLayer."

SoftDiceLossLayer::usage = 
"SoftDiceLossLayer[dim] represents a net layer that computes the SoftDice loss by comparing input class probability vectors with the target class vector."

BrierLossLayer::usage = 
"BrierLossLayer[dim] represents a net layer that computes the Brier loss by comparing input class probability vectors with the target class vector."


ClassEncoder::usage = 
"ClassEncoder[label] encodes Integer label data of 1 to max value of label into a nClass vector of 1 and 0 as the last dimension.
ClassEncoder[label, nClass] encodes Integer label data of 1 to nCalss into a nClass vector of 1 and 0 as the last dimension."

ClassDecoder::usage = 
"ClassDecoder[probability] decodes a probability vector of 1 and 0 into Integers of 1 to the value of the last dimension of probability.
ClassDecoder[probability, nClass] decodes a probability vector of 1 and 0 into Integers of 1 to nClass."


DiceSimilarity::usage = 
"DiceSimilarity[ref, pred] gives the Dice Similarity between 1 and 0 of segmentations ref and pred for class equals 1.
DiceSimilarity[x, y, class] gives the Dice Similarity of segmentations ref and pred for class.
DiceSimilarity[x, y, {class, ..}] gives the Dice Similarity of segmentations ref and pred for the list of gives classes."

MeanSurfaceDistance::usage = 
"MeanSurfaceDistance[ref, pred] gives the mean surface distance of segmentations ref and pred for class equals 1 in voxels.
MeanSurfaceDistance[x, y, class] gives the mean surface distance of segmentations ref and pred for class in voxels.
MeanSurfaceDistance[x, y, {class, ..}] gives the mean surface distance of segmentations ref and pred for the list of gives classes in voxels.
MeanSurfaceDistance[x, y, class , vox] gives the mean surface distance of segmentations ref and pred for class in milimeter.
MeanSurfaceDistance[x, y, {class, ..}, vox] gives the mean surface distance of segmentations ref and pred for the list of gives classes in milimeters."


SegmentData::usage = 
"SegmentData[data, net] segements the data using the pretrained net."


PatchesToData::usage = 
"PatchesToData[patches, ran] creates a continous dataset from the patches. For each patch the range in the data nees to be specified in ran.
The patches are have dimensions {x, y, z} each and ran is speciefied as {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}}.
PatchesToData[patches, ran, dim] creates a continous dataset from the patches with dimensions dim."

PatchesToSegmentation::usage = 
"PatchesToSegmentation[patches, ran]
PatchesToSegmentation[patches, ran, dim] ...
"

DataToPatches::usage =
"DataToPatches[data, patchSize] creates the maximal number of patches with patchSize from data, where the patches have minimal overlap.
DataToPatches[data, patchSize, n] gives n random patches from the maximal number of patches with patchSize from data, where the patches have minimal overlap."


GetTrainData::usage =
"GetTrainData[data, batchsize, patch] creates a training batch of size batchsize with patchsize patch. 
The input data can be out of memory in the form of a list of \"*wxf\" files that contain the data, segmentation and voxel size or a list of \"*.nii\" files in the form
{{\"data.nii\", \"segmentation.nii\"}..}. The input data can be in memory in a list in the form {{data, segmentation, vox}..}
GetTrainData[data, batchsize, patch, nClass] If nClass is set to an value n > 0 the segmentations are decoded in n classes."


AugmentTrainingData::usage = 
"AugmentTrainingData[{data, segmentation}, vox] augments the data and segmentation in the same way.
AugmentTrainingData[{data, segmentation}, vox, aug] by setting aug to True or False the autmentation can be turend on or off."


MakeChannelClassImage::usage = 
"MakeChannelClassImage[data, label]\[IndentingNewLine]MakeChannelClassImage[data, label, {off,max}]\[IndentingNewLine]MakeChannelClassImage[data, label, vox]
MakeChannelClassImage[data, label, {off,max}, vox]
..."


MakeChannelImage::usage = 
"MakeChannelImage[label]\[IndentingNewLine]MakeChannelImage[label, {off,max}]\[IndentingNewLine]MakeChannelImage[label, vox]
MakeChannelImage[label, {off,max}, vox]
..."


MakeClassImage::usage = 
"MakeClassImage[label]\[IndentingNewLine]MakeClassImage[data, vox]
..."


CropDatSeg::usage=
"CropDatSeg..."


(* ::Subsection:: *)
(*Options*)


BlockType::usage = 
"BlockType is an option for MakeUnet. It specifies which block are used to build the network. 
Values can be \"UNET\", \"ResNet\", \"UResNet\", \"DenseNet\" or \"UDenseNet\"."

DropoutRate::usage = 
"DropoutRate is an option for MakeUnet. It specifies how musch dropout is used after each block. It is a value between 0 and 1, default is .2."

NetworkDepth::usage = 
"NetworkDepth is an option for MakeUnet. It specifief how deep the UNET will be."

DownsampleSchedule::usage = 
"DownsampleSchedule is an option for MakeUnet. It defines how the data is downsampled for each of the deeper layers of the Unet. 
By default is is a factor two for each layer. A custum schedual for a 4 layer 3D Unet could be {{2,2,2},{1,2,2},{2,2,2},{1,2,2}}."

InputFilters::usage = 
"InputFilters is an option for MakeUnet. It defines the amount of convolutional filters of the the first UNET block."

ActivationType::usage = 
"InputFilters is an option for MakeUnet. It sepecifies which activation layer is used in the network. It can be \"LeakyRELU\" or any type allowed 
by a \"name\" definition in ElementwiseLayer."


PatchesPerSet::usage =
"PatchesPerSet is an option for GetTrainData. Defines how many random patches per dataset are created within the batch."

AugmentData::usage = 
"AugmentData is an option for GetTrainData. If set True the trainingdata is augemnted"


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*UNET*)


(* ::Subsubsection::Closed:: *)
(*MakeUnet*)


Options[MakeUnet] = {
	BlockType -> "ResNet", 
	DropoutRate -> 0.2, 
	NetworkDepth -> 5, 
	DownsampleSchedule -> Automatic, 
	InputFilters -> 32, 
	ActivationType -> "GELU"
}

SyntaxInformation[MakeUnet] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

MakeUnet[nChan_, nClass_, dimIn_, OptionsPattern[]] := Block[{
		dep, drop, type, dim, filt, enc, dec, stride, filtIn, actType
	},

	enc ="enc_" <> ToString[#]&;
	dec ="dec_" <> ToString[#]&;
	{dep, drop, type, stride, filtIn, actType} = OptionValue[
		{NetworkDepth, DropoutRate, BlockType, DownsampleSchedule, InputFilters, ActivationType}
	];
	
	dim = Switch[Length@dimIn, 2, "2D", 3, "3D"];
	filt = Switch[type, 
		"DenseNet" | "UDenseNet", Table[{filtIn, 1 + i}, {i, {1, 2, 4, 6, 8}}],
		_, filtIn {1, 2, 4, 8, 16}
	];
	stride = Prepend[If[stride===Automatic, ConstantArray[2, {dep-1, 3}], stride], {1, 1, 1}];

	NetGraph[
		Association@Join[
			Table[
				enc[i] -> ConvNode[filt[[i]], 
					"Dropout" -> drop, "Dimensions" -> dim, "Stride" -> stride[[i]],
					"ConvType" -> type, "NodeType" -> "Encode", "ActivationType" -> actType
				]
			, {i, 1, dep}],
			Table[
				dec[i] -> ConvNode[filt[[i]], 
					"Dropout" -> drop, "Dimensions" -> dim, "Stride" -> stride[[i+1]],
					"ConvType" -> type, "NodeType" -> "Decode", "ActivationType" -> actType
				]
			, {i, 1, dep - 1}],
			{"map" -> ClassMap[dim, nClass]}
		],
		Join[
			Table[If[i === 1, NetPort["Input"] -> enc[i], enc[i - 1] -> enc[i]], {i, 1, dep}],
			Table[If[i === dep - 1,	{enc[i + 1], enc[i]} -> dec[i],	{dec[i + 1], enc[i]} -> dec[i] ], {i, 1, dep - 1}],
			{"dec_1" -> "map"}],
		"Input" -> Prepend[dimIn, nChan]
	]
]


(* ::Subsubsection::Closed:: *)
(*ClassMap*)


ClassMap[dim_, nClass_] :=  Flatten[{
	ConvolutionLayer[nClass, 1], If[nClass > 1,	
		{TransposeLayer[If[dim === "2D", {3, 1, 2}, {4, 1, 2, 3}]], SoftmaxLayer[]},
		{LogisticSigmoid, FlattenLayer[1]}
	]
}]


(* ::Subsubsection::Closed:: *)
(*ConvNode*)


Options[ConvNode] = {
	"Dimensions" -> "3D",
	"ActivationType" -> "GELU",
	"Dropout" -> 0.2,
	"ConvType" -> "ResNet",
	"NodeType" -> "Encode",(*encode, decode, start*)
	"Stride" -> Automatic
};

ConvNode[chan_, OptionsPattern[]] := Block[{
		convType, nodeType, actType, mode, node, drop, dim, stride
	},
	{convType, nodeType, actType, drop, dim, stride} = OptionValue[{"ConvType", "NodeType", "ActivationType", "Dropout", "Dimensions", "Stride"}];

	(*mode is encoding or decoding, decoding is solved later and treated as normal here*)
	mode = If[nodeType === "Encode", "down", "normal"];

	(*make convblocks for various convolution types*)
	node = Switch[convType,	
		"UResNet", 
		Flatten[{
			ConvBlock[chan/2, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], 
			ConvBlock[chan, "ActivationType" -> actType, "Stride"->stride]
		}],

		"ResNet",
		{<|
			"con" -> Join[
				ConvBlock[chan/2, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], 
				ConvBlock[chan, "ActivationType" -> "None"]], 
			"skip" -> ConvBlock[chan, "ConvMode" -> mode<>"S", "ActivationType" -> "None", "Stride"->stride],
			"tot" -> {TotalLayer[], ActivationLayer[actType]}
		|>, {
			{"con", "skip"} -> "tot"
		}},

		"DenseNet",
		With[{n = chan[[1]], dep = chan[[2]], layName = "lay_" <> ToString[#] &},{
			Join[
				<|If[mode === "down", "down"->ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], Nothing]|>,
				Association@Table[If[rep==dep, "lay_end", layName[rep]] -> ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> "catenate"], {rep, 1, dep}]
			],
			Table[Table[If[rr == 0, If[mode==="down", "down", NetPort["Input"]], layName[rr]], {rr, 0, rep - 1}] -> If[rep==dep, "lay_end", layName[rep]], {rep, 1, dep}]
		}],

		"UDenseNet", 
		Flatten[{If[mode === "down", ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], Nothing], ConstantArray[ConvBlock[chan[[1]], "ActivationType" -> actType], chan[[2]]]}],

		_,
		Flatten[{ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> mode, "Stride"->stride], ConvBlock[chan, "ActivationType" -> actType]}]
	];


	(*Add dropout and upconv for deconding block*)
	If[nodeType === "Decode",

		(*convert to decoding block and add dropout*)
		NetGraph[<|
			"upconv" -> ConvBlock[chan, "ActivationType" -> actType, "ConvMode" -> "up", Dimensions -> dim, "Stride"->stride],
			"conv" -> If[convType==="ResNet"||convType==="DenseNet",
				NetGraph[
					Join[node[[1]], <|"cat"->CatenateLayer[],"drop"->DropoutLayer[drop]|>],
					Switch[convType,
						"ResNet", Join[node[[2]], {"cat"->{"con","skip"}, "tot"->"drop"}],
						"DenseNet", Join[node[[2]] /. NetPort["Input"]->"cat", {"lay_end"->"drop"}]
					]
				],
				Flatten[{CatenateLayer[], node, DropoutLayer[drop]}]
			]
		|>, {
			{NetPort["Input2"] -> "upconv", NetPort["Input1"]} -> "conv"
		}]
		,

		(*add dropout to encoding block*)
		If[convType==="ResNet"||convType==="DenseNet",
				NetGraph[
					Join[node[[1]], <|"drop"->DropoutLayer[drop]|>],
					Join[node[[2]], {Switch[convType,"ResNet", "tot", "DenseNet", "lay_end"]->"drop"}]
				],
				NetChain[Flatten@{node, DropoutLayer[drop]}]
			]
		
	]
]


(* ::Subsubsection::Closed:: *)
(*ConvBlock*)


Options[ConvBlock] = {
	"Dimensions" -> "3D",
	"ActivationType" -> "GELU",
	"ConvMode" -> "normal"(*normal, up, down, catenate*),
	"Stride" -> 2
};

ConvBlock[channels_, OptionsPattern[]] := Block[{chan, kern,  actType, pad, actLayer, convMode, dim, str},
	{actType, convMode, dim, str} = OptionValue[{"ActivationType", "ConvMode", "Dimensions", "Stride"}];
	chan = Round@First@Flatten@{channels};
	
	Switch[convMode,
		"up", 
		{ResizeLayer[Scaled/@str, Resampling -> "Nearest"], ConvolutionLayer[chan, 2, "PaddingSize" -> ConstantArray[{0,1},Length[str]], "Stride" -> 1]},
		"down"|"downS", 
		{ConvolutionLayer[chan, str, "PaddingSize" -> 0, "Stride" -> str], BatchNormalizationLayer[],  ActivationLayer[actType]},
		"normal", 
		{ConvolutionLayer[chan, 3, "PaddingSize" -> 1, "Stride" -> 1], BatchNormalizationLayer[],  ActivationLayer[actType]},
		"normalS", 
		{ConvolutionLayer[chan, 1, "PaddingSize" -> 0, "Stride" -> 1], BatchNormalizationLayer[],  ActivationLayer[actType]},
		"catenate", 
		{CatenateLayer[], ConvolutionLayer[chan, 3, "PaddingSize" -> 1, "Stride" -> 1], BatchNormalizationLayer[],  ActivationLayer[actType]}
	]
]


(* ::Subsubsection::Closed:: *)
(*ActivationLayer*)


ActivationLayer[actType_] := Switch[actType, "LeakyRELU", ParametricRampLayer[], "None", Nothing, _, ElementwiseLayer[actType]]


(* ::Subsection:: *)
(*LossLayers*)


(* ::Subsubsection::Closed:: *)
(*AddLossLayer*)


SyntaxInformation[AddLossLayer] = {"ArgumentsPattern" -> {_}};

AddLossLayer[net_]:=Block[{dim},
	dim = Length[Information[net,"OutputPorts"][[1]]]-1;
	NetGraph[<|
		"net"->net,
		"SoftDice" -> SoftDiceLossLayer[dim],
		"CrossEntropy" -> CrossEntropyLossLayer["Probabilities"],
		"Brier" -> {BrierLossLayer[dim], FunctionLayer[10 #&]}
		|>,{
		NetPort["Input"]->"net"->NetPort["Output"],
		{"net",NetPort["Target"]}->"SoftDice"->NetPort["SoftDice"],
		{"net",NetPort["Target"]}->"CrossEntropy"->NetPort["CrossEntropy"],
		{"net",NetPort["Target"]}->"Brier"->NetPort["Brier"]
	}]
]


(* ::Subsubsection::Closed:: *)
(*SoftDiceLossLayer*)


SyntaxInformation[SoftDiceLossLayer] = {"ArgumentsPattern" -> {_}};

SoftDiceLossLayer[dim_]:=NetGraph[<|
	"times" -> ThreadingLayer[Times],
	"flattot1" -> FlatTotLayer[dim - 1],
	"flattot2" -> FlatTotLayer[dim - 1],
	"flattot3" -> FlatTotLayer[dim - 1],
	"total" -> TotalLayer[],
	"devide" -> {ThreadingLayer[Divide], AggregationLayer[Mean, 1], ElementwiseLayer[1 - 2 # &]},
	"weight" -> ElementwiseLayer[1/(# + 1) &],
	"times1" -> ThreadingLayer[Times],
	"times2" -> ThreadingLayer[Times]
|>, {
	{NetPort["Input"], NetPort["Target"]} -> "times" -> "flattot1",
	NetPort["Input"] -> "flattot2",
	NetPort["Target"] -> "flattot3",
	{"flattot2", "flattot3"} -> "total", "flattot3" -> "weight",
	{"flattot1", "weight"} -> "times1",
	{"total", "weight"} -> "times2",
	{"times1", "times2"} -> "devide" -> NetPort["Loss"]
}, "Loss" -> "Real"]


(* ::Subsubsection::Closed:: *)
(*BrierLossLayer*)


SyntaxInformation[BrierLossLayer] = {"ArgumentsPattern" -> {_}};
   
BrierLossLayer[dim_] := NetGraph[<|
	"sub" -> ThreadingLayer[Subtract],
	"SqMn" -> {ElementwiseLayer[#^2 &], FlattenLayer[dim - 1], TransposeLayer[], AggregationLayer[Mean]},
	"weigth" -> {FlatTotLayer[dim - 1], ElementwiseLayer[1/(# + 1) &]},
	"times" -> ThreadingLayer[Times],
	"tot1" -> AggregationLayer[Total, 1],
	"tot2" -> AggregationLayer[Total, 1],
	"devide" -> ThreadingLayer[Divide]
|>, {
	{NetPort["Input"], NetPort["Target"]} -> "sub" -> "SqMn",
	NetPort["Target"] -> "weigth",
	{"weigth", "SqMn"} -> "times" -> "tot1",
	"weigth" -> "tot2",
	{"tot1", "tot2"} -> "devide" -> NetPort["Loss"]
}, "Loss" -> "Real"]


(* ::Subsubsection::Closed:: *)
(*FlatTotLayer*)


FlatTotLayer[lev_]:=NetChain[{FlattenLayer[lev],AggregationLayer[Total,1]}];


(* ::Subsection:: *)
(*Encoders*)


(* ::Subsubsection::Closed:: *)
(*ClassEndocer*)


SyntaxInformation[ClassEncoder] = {"ArgumentsPattern" -> {_, _.}};

ClassEncoder[data_]:= If[nClass === 1, data, ClassEncoderC[data, Max@data]]

ClassEncoder[data_, nClass_]:= If[nClass === 1, data, ClassEncoderC[data, nClass]]

ClassEncoderC = Compile[{{data, _Integer, 2}, {n, _Integer, 0}},
	Transpose[1 - Unitize[ConstantArray[data, n] - Range[n]], {3, 1, 2}]
,RuntimeAttributes -> {Listable}]


(* ::Subsubsection::Closed:: *)
(*ClassDecoder*)


SyntaxInformation[ClassDecoder] = {"ArgumentsPattern" -> {_, _.}};

ClassDecoder[data_]:=ClassDecoderC[data, Last@Dimensions@data]

ClassDecoder[data_, nClass_]:=ClassDecoderC[data, nClass]

ClassDecoderC = Compile[{{data, _Real, 1}, {n, _Integer, 0}},  
	Total[Range[n] (1 - Unitize[Chop[(data/Max[data]) - 1]])]
, RuntimeAttributes -> {Listable}]


(* ::Subsection:: *)
(*Distance measures*)


(* ::Subsubsection::Closed:: *)
(*DiceSimilarity*)


SyntaxInformation[DiceSimilarity] = {"ArgumentsPattern" -> {_, _, _}};

DiceSimilarity[ref_, pred_, nClasses_?ListQ] := Table[DiceSimilarity[ref, pred, c], {c, nClasses}]

DiceSimilarity[ref_, pred_] := DiceSimilarityC[Flatten[ref], Flatten[pred], 1]

DiceSimilarity[ref_, pred_, c_?IntegerQ] := DiceSimilarityC[Flatten[ref], Flatten[pred], c]


DiceSimilarityC = Compile[{{ref, _Integer, 1}, {pred, _Integer, 1}, {class, _Integer, 0}}, Block[{refv, predv, denom},
	refv = Flatten[1 - Unitize[ref - class]];
	predv = Flatten[1 - Unitize[pred - class]];
	denom = (Total[refv] + Total[predv]);
	If[denom === 0., 1., N[2 Total[refv predv] / denom]]
 ], RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MeanSurfaceDistance*)


SyntaxInformation[MeanSurfaceDistance] = {"ArgumentsPattern" -> {_, _, _, _.}};

MeanSurfaceDistance[ref_, pred_] := MeanSurfaceDistance[ref, pred, 1, {1, 1, 1}]

MeanSurfaceDistance[ref_, pred_, class_?IntegerQ] := MeanSurfaceDistance[ref, pred, class, {1, 1, 1}]

MeanSurfaceDistance[ref_, pred_, nClasses_?ListQ] := MeanSurfaceDistance[ref, pred, nClasses, {1, 1, 1}]

MeanSurfaceDistance[ref_, pred_, nClasses_?ListQ, vox_] := Table[MeanSurfaceDistance[ref, pred, class, vox], {class, nClasses}]

MeanSurfaceDistance[ref_, pred_, class_?IntegerQ, vox_] := Block[{coorRef, coorPred, fun},
	coorRef = Transpose[vox Transpose[GetEdge[ref, class]["ExplicitPositions"]]];
	coorPred = Transpose[vox Transpose[GetEdge[pred, class]["ExplicitPositions"]]];
	If[coorRef==={}||coorPred==={},
		"noSeg",
		fun = Nearest[coorRef];
		Mean@Sqrt@Total[(fun[coorPred,1][[All,1]]-coorPred)^2,{2}]
	]
]


(* ::Subsubsection::Closed:: *)
(*GetEdge*)


GetEdge[lab_, class_] := Block[{seg, n, per, pts},
	seg = 1 - Unitize[lab - class];
	n = Length[seg];
	per = Flatten[Table[
		pts = Ceiling[Flatten[ComponentMeasurements[Image[seg[[i]]], "PerimeterPositions"][[All, 2]], 2]];
		If[pts === {}, Nothing, Join[ConstantArray[i, {Length[pts], 1}], pts[[All, {2, 1}]], 2]]
	, {i, 1, n}], 1];
	Reverse[SparseArray[per -> 1, Dimensions[seg]], 2]
]


(* ::Subsection:: *)
(*SegmentData*)


(* ::Subsubsection:: *)
(*SegmentData*)


Options[SegmentData]={TargetDevice->"GPU"}

SyntaxInformation[SegmentData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

SegmentData[data_, net_, OptionsPattern[]]:=Block[{dev, inp, out, dim, patch, pts, dimN, seg},
	dev = OptionValue[TargetDevice];

	inp = NetDimensions[net,"Input"];
	out = Rest[NetDimensions[net, "LastEncoding"]];
	dim = Dimensions[data];
	dimN = Rest[inp]/out;

	(*if patching padding is needed for bad edge behaviour of CNN*)
	lim = 112;
	pad = 8;

	If[CubeRoot[N[Times @@ dim]] > lim && dev=!="GPU",
		(*run data on small CPU*)
		dim = (dim+2 pad);
		patch = Min[{lim, #}] & /@ (out Ceiling[dim/out]);
		patch = dimN Ceiling[patch/dimN];
		Echo[patch, "Data is probably to big for memmory. Patching with patchsize:"];

		(*patch padded data and run each patch seperate*)	
		{patch, pts} = DataToPatches[ArrayPad[data, pad, 0.], patch, PatchNumber -> 0];
		patch = SegmentData[#, net, TargetDevice->dev] & /@ patch;

		(*contract all the patches and thus the pts, keep expanded dim for depatching, then contract*)
		ArrayPad[PatchesToSegmentation[ArrayPad[#, -pad]&/@patch, Map[# + {pad, -pad} &, pts, {2}], dim], -pad]
		,
		dimN = dimN Ceiling[dim/dimN];
		(*run full data on big GPU*)
		seg = {PadToDimensions[NormDat[data], dimN, PadDirection->"Right"]};
		(*apply network*)
		seg = ClassDecoder@NetReplacePart[net, "Input" -> Prepend[dimN, First@inp]][seg, TargetDevice->dev];
		(*decode output and crop to original dimensions*)
		Round[PadToDimensions[seg, dim, PadDirection->"Right"] - 1]
	]
]


(* ::Subsubsection:: *)
(*NetDimensions*)


NetDimensions[net_, port_]:=Switch[port,
	"Input",Information[net,"InputPorts"]["Input"],
	"Output",Information[net,"OutputPorts"]["Output"],
	"LastEncoding",Information[NetTake[net,Last[Select[Keys[net[[All,1]]],StringContainsQ[#,"enc_"]&]]],"OutputPorts"]["Output"]
] 


(* ::Subsection:: *)
(*Prepare Data*)


CropDatSeg[dat_,seg_,labs_]:=Block[{datO, segO, lab, dim,cr},
	dim=Dimensions@dat;
	cr=FindCrop[dat Mask[NormalizeData[dat],5,MaskDilation->1]];
	datO=NormDat[ApplyCrop[dat,cr]];
	{segO,lab}=NormSeg[ApplyCrop[seg,cr],labs];
	
	{{datO,segO},{lab,dim,Dimensions@datO,cr}}
]


NormDat = Compile[{{dat, _Real, 3}}, Block[{data},
	If[Min[dat]=!=Max[dat],
		data = Flatten[dat];
		0.5 dat/Median[Pick[data, Unitize[data], 1]],
		dat
	]
], RuntimeOptions -> "Speed", RuntimeAttributes->Listable];


NormSeg[seg_, labs_]:=Block[{zero, segT, labT, sel},
	zero=0 seg;
	{segT,labT}=SplitSegmentations[seg];
	
	sel=Flatten[(Flatten[Position[labT,#]]&/@labs)/.{}->{0}];
	segT=Transpose[If[#===0,zero,segT[[All,#]]]&/@sel];
	
	{MergeSegmentations[segT,Range[Length[sel]]],Unitize[sel]labs}
]


(* ::Subsection:: *)
(*Data Patching*)


(* ::Subsubsection::Closed:: *)
(*PatchesToData*)


SyntaxInformation[PatchesToData] = {"ArgumentsPattern" -> {_, _, _.}};

PatchesToData[patches_, ran_]:=PatchesToData[patches, ran,  Max/@Transpose[ran[[All,All,2]]]]

PatchesToData[patches_, ran_, dim_]:=Block[{sel,zero,a1,a2,b1,b2,c1,c2,out},
	zero = ConstantArray[0., dim];
	out = MapThread[(
		out = zero;
		{{a1,a2},{b1,b2},{c1,c2}} = #2;
		out[[a1;;a2, b1;;b2, c1;;c2]] = #1;
		out
	)&,{patches, ran}];

	DevideNoZero[N[Total[out]], N[Total[Unitize /@ out]]]
]


(* ::Subsubsection::Closed:: *)
(*PatchesToSegmentation*)


SyntaxInformation[PatchesToSegmentation] = {"ArgumentsPattern" -> {_, _, _.}};

PatchesToSegmentation[patches_, ran_] := PatchesToSegmentation[patches, ran, Max /@ Transpose[ran[[All, All, 2]]]]

PatchesToSegmentation[patches_, ran_, dim_] := Block[{segs, labs, allLab, pos, roi, seg},
	{segs, labs} = Transpose[SplitSegmentations /@ patches];
	segs = Transpose /@ segs; (*first index of segs is patches and second is labs*)
	
	allLab = Sort[DeleteDuplicates[Flatten[labs]]];
	pos = Position[labs, #] & /@ allLab;
	
	segs = PatchesToData[N[Normal[segs[[#[[1]], #[[2]]]] + 1 & /@ #]], ran[[#[[All,1]]]], dim] & /@ pos;
	MergeSegmentations[ Transpose[SparseArray[Round[Ramp[# - 1]]] & /@ segs], allLab]
]


(* ::Subsubsection::Closed:: *)
(*DataToPatches*)


Options[DataToPatches] = {PatchNumber->2}

SyntaxInformation[DataToPatches] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, opts:OptionsPattern[]]:=DataToPatches[dat, patch, "All", opts] 

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}]:={GetPatch[dat, patch, pts], pts}

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, nPatch_, OptionsPattern[]]:=Block[{ptch, pts},
	pts = GetPatchRanges[dat, patch, If[IntegerQ[nPatch], nPatch, "All"], OptionValue[PatchNumber]];
	ptch = GetPatch[dat, patch, pts];
	{ptch, pts}
] 


(* ::Subsubsection::Closed:: *)
(*GetPatch*)


GetPatch[dat_, patch:{_,_,_}, pts:{{{_,_},{_,_},{_,_}}..}]:=GetPatch[dat, patch, #]&/@pts

GetPatch[dat_, patch:{_,_,_}, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}]:=PadRight[dat[[i1;;i2,j1;;j2,k1;;k2]], patch, 0.]


(* ::Subsubsection::Closed:: *)
(*GetPatchRanges*)


GetPatchRanges[dat_, patch_, nPatch_, nRan_]:=Block[{pts},
	pts = GetPatchRangeI[Dimensions[dat], patch, nRan];
	If[nPatch==="All", pts, RandomSample[pts, Min[{Length@pts, nPatch}]]]
]


GetPatchRangeI[datDim_?ListQ, patchDim_?ListQ, n_]:=Tuples@MapThread[GetPatchRangeI[#1,#2, n]&, {datDim, patchDim}]

GetPatchRangeI[d_?IntegerQ, p_?IntegerQ, n_]:=Block[{i,st},
	i = Ceiling[d/p]+n;
	If[i>1,
		st=Round[Range[0, 1, 1./(i - 1)](d-p)]+1;
		Thread[{st, st+p-1}],
		{{1,d}}
	]
]


(* ::Subsection:: *)
(*Get Train Data*)


(* ::Subsubsection::Closed:: *)
(*GetTrainData*)


Options[GetTrainData] = {
	PatchesPerSet -> 1, 
	AugmentData -> True
};

GetTrainData[datas_, nBatch_, patch_, opts:OptionsPattern[]]:=GetTrainData[datas, nBatch, patch, False, opts]

GetTrainData[datas_, nBatch_, patch_, nClass_, OptionsPattern[]] := Block[{
		itt, i, datO, segO, dat, seg, vox, dim, aug, ran, nSet
	},

	itt = 0;
	datO = segO = {};

	{aug, nSet} = OptionValue[{AugmentData, PatchesPerSet}];
	
	aug = If[BooleanQ[aug], # && aug & /@ {True, True, True, True, False, False, False}, PadRight[aug, 7, False]];

	itt = Ceiling[nBatch/nSet];

	Do[
		dat = RandomChoice[datas];
		
		If[StringQ[dat], 
			(*data is wxf file format*)
			{dat, seg, vox} = Import[dat];
			,
			If[Length[dat]===2, 
				(*datas is list of nii files {dat.nii, seg.nii}*)
				{dat, vox} = ImportNii[dat[[1]]];
				{seg, vox} = ImportNii[dat[[2]]];
				,
				(*data is in memory*)
				{dat, seg, vox} = dat;
			]
		];

		dim = Max /@ Transpose[{Dimensions@dat, patch}];
		dat = PadToDimensions[dat, dim];
		seg = PadToDimensions[seg, dim];

		ran = RandomChoice[{True, True, True, False}];
		{dat, seg} = AugmentTrainingData[{dat, seg}, vox, ran && # & /@ aug];
		{dat, seg} = PatchTrainingData[{dat, seg}, patch, nSet];

		datO = Join[datO, dat];
		segO = Join[segO, seg];
	, itt];

	If[IntegerQ[nClass],
		Thread[Transpose[{datO[[;; nBatch]]}] -> ClassEncoder[segO[[;; nBatch]] + 1, nClass]],
		Thread[Transpose[{datO[[;; nBatch]]}] -> segO[[;; nBatch]] + 1]
	]
];


(* ::Subsubsection::Closed:: *)
(*AugmentTrainingData*)


AugmentTrainingData[{dat_, seg_}, vox_] := AugmentTrainingData[{dat, seg}, vox, {True, True, True, True, False, False, False}]

AugmentTrainingData[{dat_, seg_}, vox_, aug_?BooleanQ] := AugmentTrainingData[{dat, seg}, vox, {aug, aug, aug, aug, False, False, False}]

AugmentTrainingData[{dat_, seg_}, vox_, {flip_, rot_, trans_, scale_, noise_, blur_, bright_}] := Block[{datT, segT, w, r, t, s},
	datT = ToPackedArray[N[dat]];
	segT = ToPackedArray[N[seg]];
	
	(*Augmentation of mirroring*)
	If[flip && RandomChoice[{True, False}], {datT, segT} = ReverseC[{datT, segT}]];
	
	(*Augmentation of orientation and scale*)
	If[rot || trans || scale,
		r = {If[rot, RandomReal[{-90, 90}] , 0.], 0., 0.};
		t = If[trans, (Dimensions[dat] vox) {0, 1, 1} RandomReal[{-0.1, 0.1}, 3], {0., 0., 0.}];
		s = If[scale, Append[RandomReal[{0.5, 1.5}, 2], 1.], {1., 1., 1.}];
		w = Join[r, t, s, {0., 0., 0.}];
		datT = DataTransformation[datT, vox, w, InterpolationOrder -> 0];
		segT = DataTransformation[segT, vox, w, InterpolationOrder -> 0];
	];
	
	(*Augmentations of sharpness intensity and noise*)
	If[bright, datT = RandomChoice[{RandomReal[{1, 1.5}], 1/RandomReal[{1, 1.5}]}] datT];
	If[blur, datT = GaussianFilter[datT, RandomReal[{0.1, 1.5}]]];
	If[noise && RandomChoice[{True, False, False}], datT = AddNoise[datT, 0.5/RandomReal[{10, 50}]]];
	
	{ToPackedArray[N[datT]], ToPackedArray[Round[segT]]}
]


ReverseC = Compile[{{dat, _Real, 1}}, Reverse[dat], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection:: *)
(*PatchTrainingData*)


PatchTrainingData[{dat_,seg_}, patch_, n_]:=Block[{pts,datP,segP},
	{datP, pts}=DataToPatches[dat, patch, n];
	segP=DataToPatches[seg, patch, pts][[1]];
	{ToPackedArray[N[#]]&/@datP,ToPackedArray[Round[#]]&/@segP}
]


(* ::Subsubsection:: *)
(*PatchTrainingData*)


SyntaxInformation[MakeChannelClassImage]={"ArgumentsPattern"->{_, _, _., _.}};

MakeChannelClassImage[data_,label_]:=MakeChannelClassImage[data,label,MinMax[label],{1,1,1}]

MakeChannelClassImage[data_,label_,{off_,max_}]:=MakeChannelClassImage[data,label,{off,max},{1,1,1}]

MakeChannelClassImage[data_,label_,vox_]:=MakeChannelClassImage[data,label,MinMax[label],vox]

MakeChannelClassImage[data_,label_,{off_,max_},vox_]:=Block[{i1,i2},
	i1=MakeClassImage[label,vox];
	i2=MakeChannelImage[data,vox];
	ImageCollage[ImageCompose[#,SetAlphaChannel[i1,0.5AlphaChannel[i1]]]&/@i2]
]


SyntaxInformation[MakeClassImage]={"ArgumentsPattern"->{_, _., _.}};

MakeClassImage[label_]:=MakeClassImage[label,MinMax[label],{1,1,1}]

MakeClassImage[label_,{off_,max_}]:=MakeClassImage[label,{off,max},{1,1,1}]

MakeClassImage[label_,vox_]:=MakeClassImage[label,MinMax[label],vox]

MakeClassImage[label_,{off_,max_}, vox_]:=Block[{cols, im, rat},
	cols=Prepend[ColorData["DarkRainbow"][#]&/@Rescale[Range[off+1,max]],Transparent];
	
	im=If[ArrayDepth[label]===3,label[[Round[Length@label/2]]],label]-off+1;
	rat=vox[[{2,3}]]/Min[vox[[{2,3}]]];
	
	ImageResize[Image[cols[[#]]&/@im], Round@Reverse[rat Dimensions[im]], Resampling->"Nearest"]
]


SyntaxInformation[MakeChannelImage]={"ArgumentsPattern"->{_, _., _.}};

MakeChannelImage[data_]:=MakeChannelImage[data,{1,1,1}]

MakeChannelImage[data_,vox_]:=Block[{dat, im, rat},
	dat = NormDat[data];
	rat=vox[[{2,3}]]/Min[vox[[{2,3}]]];
	(
		im=#;
		im=If[ArrayDepth[#]===3, im[[Round[Length@im/2]]], im];
		ImageResize[Image[Clip[im,{0,1}]],Round@Reverse[rat Dimensions[im]],Resampling->"Nearest"]
	)&/@dat
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
