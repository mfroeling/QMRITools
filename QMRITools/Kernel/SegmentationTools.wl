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


GetNeuralNet::usage = 
"GetNeuralNet[name] loads a pretrained neural net that come with the toolbox. Current named nets 
are \"LegSide\", \"LegSide\", \"SegThighMuscle\", \"SegLegMuscle\", and \"SegLegBones\". The loading is cashed within a session."


MakeUnet::usage = 
"MakeUnet[nClasses, dimIn] Generates a UNET with one channel as input and nClasses as output. 
MakeUnet[nChannels, nClasses, dimIn] Generates a UNET with nChannels as input and nClasses as output. 
he number of parameter of the first convolution layer can be set with dep.\n
The data dimensions can be 2D or 3D and each of the dimensions should be 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240 or 256.
..."

MakeNode::usage = 
"MakeNode[] ..."

MakeClassifyNetwork::usage = 
"MakeClassifyNetwork[classes] makes a classify network with three convolusion layers and 3 fully connected layers. 
The input classes should be a list of strings. The imput image dimensions should not be smaller thand 64x64."

MakeClassifyImage::usage =
"MakeClassifyImage[data] makes a image of the input data. The data is automatically cropped to remove the background and normalized.
If the input data is 3D a list of images is returned."


PrintKernels::usage = 
"PrintKernels[net] gives a short summary of the convolution kernels and array elements in the network."

NetDimensions::usage = 
"NetDimensions[net] extracts the input channels, output classes, the input patch dimension, and the number of input filters."

ChangeNetDimensions::usage =
"ChangeNetDimensions[netIn] changes input channels, output classes, the input patch dimension of the input network netIn."


AddLossLayer::usage = 
"AddLossLayer[net] adds three loss layers to a NetGraph, a SoftDiceLossLayer, BrierLossLayer and a CrossEntropyLossLayer."

SoftDiceLossLayer::usage = 
"SoftDiceLossLayer[dim] represents a net layer that computes the SoftDice loss by comparing input class probability vectors with the target class vector."


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
"SegmentData[data, what] segements the data. The what specifies the segmentation to be done.
It currently allows for \"LegBones\" for the bones or \"Legs\" for the muscles."

ApplySegmentationNetwork::usage = 
"ApplySegmentationNetwork[data, net] segements the data using the pretrained net."

ClassifyData::usage = 
"ClassifyData[data, method] classifies the input data using the given method. The data is converted to images using MakeClassifyImages.
The input method can be a filename of a classify network or a classify network. 
Additionally the input method can be one of the predefined methods \"LegPosition\" or \"LegSide\"."


ShowTrainLog::usage =
"ShowTrainLog[fol] shows the training log of the network training stored in fol.
ShowTrainLog[fol, min] shows the training log of the network training stored in fol if the log has a minimum of min rounds."

TrainSegmentationNetwork::usage =
"TrainSegmentationNetwork[{inFol, outFol}] trains a segmentation network. The correctly prepared training data should be stored in inFol. The progress each round will be saved in outFol.
TrainSegmentationNetwork[{inFol, outFol}, netCont] does the same but defines how to continue with netCont. If netCont is \"Start\" training will be restarted.
If netCont is a initialized network or network file (wlnet) this will be used. If netCont is a a outFol the last saved network will be used."

GetTrainData::usage =
"GetTrainData[data, batchsize, patch] creates a training batch of size batchsize with patchsize patch. 
The input data can be out of memory in the form of a list of \"*wxf\" files that contain the data, segmentation and voxel size or a list of \"*.nii\" files in the form
{{\"data.nii\", \"segmentation.nii\"}..}. The input data can be in memory in a list in the form {{data, segmentation, vox}..}
GetTrainData[data, batchsize, patch, nClass] If nClass is set to an value n > 0 the segmentations are decoded in n classes."

PrepTrainData::usage=
"PrepTrainData[data, segmentation] crops and normalizes the data and segementation such that it is optimal for training CCN for segmentation.
PrepTrainData[data, segmentation, labin] does the same but only selects the labin from the segmentation.
PrepTrainData[data, segmentation, {labin, labout}] does the same but only selects the labin from the segmentation and replaces it with labout."


DataToPatches::usage =
"DataToPatches[data, patchSize] creates the maximal number of patches with patchSize from data, where the patches have minimal overlap.
DataToPatches[data, patchSize, n] gives n random patches from the maximal number of patches with patchSize from data, where the patches have minimal overlap."

PatchesToData::usage = 
"PatchesToData[patches, ran] creates a continous dataset from the patches. For each patch the range in the data nees to be specified in ran.
The patches are have dimensions {x, y, z} each and ran is speciefied as {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}}.
PatchesToData[patches, ran, dim] creates a continous dataset from the patches with dimensions dim."


AugmentTrainingData::usage = 
"AugmentTrainingData[{data, segmentation}, vox] augments the data and segmentation in the same way.
AugmentTrainingData[{data, segmentation}, vox, aug] by setting aug to True or False the autmentation can be turend on or off.
The value aug can also be a list of boolean values contoling various augentation parameters {flip, rotate, translate, scale, noise, blur, brightness}.
The defualt settings are {True, True, True, True, False, False, False}."

AugmentImageData::usage = 
"AugmentImageData[image, {rotate, flip}] augments the input image by rotating between -180 and 180 degrees and flipping. The inputs rotate and flip
can be set to True or False.
AugmentImageData[{image, ..}, {rotate, flip}] same but for a list of images."


MakeChannelImage::usage = 
"MakeChannelImage[data] makes a crossectional image of the channels data of a training dataset generated by GetTrainData.
MakeChannelImage[data, vox] same but with the aspect ratio determined by vox."

MakeClassImage::usage = 
"MakeClassImage[label ] makes a crossectional image of the classes label of a training dataset generated by GetTrainData
MakeChannelImage[label, {b, n}] same but with explicit definition of background value b and number of classes n. 
MakeClassImage[data, vox] same but with the aspect ratio determined by vox.
MakeChannelImage[label, {b, n}, vox] same with explicit definition and aspect ratio definition."

MakeChannelClassImage::usage = 
"MakeChannelClassImage[data, label] makes a crossectional image of the channels data overlaid with a crossectional image of the classes label of a training dataset generated
MakeChannelClassImage[data, label, {off,max}] same but with explicit definition of background value b and number of classes n. 
MakeChannelClassImage[data, label, vox] same but with the aspect ratio determined by vox.
MakeChannelClassImage[data, label, {off,max}, vox] same with explicit definition and aspect ratio definition."


SplitDataForSegementation::usage = 
"SplitDataForSegementation[data] is a specific function for leg data to prepare data for segmentation. It detects the side and location and will split and label the data accordingly.
SplitDataForSegementation[data ,seg] does the same but is rather used when preparing training data. Here the seg is split in exaclty the same way as the data."


MuscleLabelToName::usage =
"MuscleLabelToName[{lab, ..}] converts list of lab, which need to be integers to names using the file GetAssetLocation[\"MusclesLegLabels\"].
MuscleLabelToName[{lab, ..}, file] does the same but uses a user defined ITKSnap label definition file."

MuscleNameToLabel::usage = 
"MuscleNameToLabel[{name, ..}] converts list of muscle names to integer labels using the file GetAssetLocation[\"MusclesLegLabels\"]
MuscleNameToLabel[{name, ..}, file] does the same but uses a user defined ITKSnap label definition file."

ImportITKLabels::usage = 
"ImportITKLabels[file] imports the ITKSnap label file."

SegmentDataGUI::usage = 
"SegmentDataGUI[] is a function that creates a graphical user interface (GUI) for segmenting data. 
It prompts the user to enter the paths for the input and output files, and allows them to select the segmentation type." 


(* ::Subsection::Closed:: *)
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

SettingSchedule::usage =
"SettingSchedule is an option for MakeUnet. It defines the settings for the Unet blocks. If one setting is given it applied to all layers.
If a list of settings is given the settings can be different per layer. The following settings are the default settings. 
\"Unet\": convblock repetitions, 2, \"ResNet\" -> convblock repetitions, 2, \"DenseNet\" -> {dense depth, block repetitions}, {4,2},
\"Inception\" -> {inception width, block repetitions}, {4,2}, \"U2Net\"-> {Unet depth, downscale}, {5, True}."

FeatureSchedule::usage = 
"FeatureSchedule is an option for MakeUnet. It defines how the number of features is upsampled for each of the deeper layers of the Unet.
By default it increases the number of features by a factor 2 each layer, i.e. {1, 2, 4, 8, 16}."

InputFilters::usage = 
"InputFilters is an option for MakeUnet. It defines the amount of convolutional filters of the the first UNET block."

ActivationType::usage = 
"ActivationType is an option for MakeUnet. It sepecifies which activation layer is used in the network. It can be \"LeakyRELU\" or any type allowed 
by a \"name\" definition in ElementwiseLayer."


PatchesPerSet::usage =
"PatchesPerSet is an option for GetTrainData. Defines how many random patches per dataset are created within the batch."

AugmentData::usage = 
"AugmentData is an option for GetTrainData and TrainSegmentationNetwork. If set True the trainingdata is augmented."

PatchSize::usage =
"PatchSize is an option for TrainSegmentationNetwork. Defines the patch size used in the network training."

RoundLength::usage = 
"RoundLength is an option for TrainSegmentationNetwork. Defines how many batches will be seen during eacht training round."


MaxPatchSize::usage = 
"MaxPatchSize is an option for SegmentData and ApplySegmentationNetwork. Defines the patch size used when segmenting data. Bigger patches are better."

DataPadding::usage = 
"DataPadding is an option for ApplySegmentationNetwork. Defines how much to pad the data patches in all directions."


PatchNumber::usage = 
"PatchNumber is an option for DataToPatches. Can be an integer value >= 0. The larger the number the more overlap the patches have.
The minimal number of patches in each direction is calculated, and then for each dimension the given number is added."

PatchPadding::usage = 
"PatchPadding is an option for DataToPatches. Can be an integer value >= 0. It padds the chosen patch size with the given number."


(* ::Subsection:: *)
(*Error Messages*)


TrainSegmentationNetwork::net = "The net input is not \"Start\", a network, a network file, or a previous train folder."

TrainSegmentationNetwork::cont = "Could not find a previous network in the specified folder."

TrainSegmentationNetwork::inp = "The string input given is not a network file or a directory."

TrainSegmentationNetwork::itt = "Not enough itterations specified for training. Remaining itterations are less than 5."


MakeUnet::scale = "";

MakeUnet::sett = "";

MakeUnet::feat = "";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*MakeUnet - new*)


(* ::Subsubsection:: *)
(*Settings*)


netDefaults={"Conv"->1,"UNet"->2,"ResNet"->2,"DenseNet"->{4,2},"Inception"->{4,2},"U2Net"->{3,True},"Map"->1};


(* ::Subsubsection::Closed:: *)
(*MakeUnet*)


(*"
features -> {layer input and output features, internal layer features}

Unet -> repetitions
ResNet -> repetitions
DenseNet -> {depth, repetitions}
Inception -> {width, repetitions}
U2Net-> {depth, downscale}

Make convblocks
- {"Conv", rep} all use featOut
- {"UNet", rep} all use featOut
- "ResNet" skip and last use featOut, first resblock uses featInt
- {"DenseNet",{dep, rep}} the dense layers use featInt, last one uses \
featOut. contains rep dense blocks of depht dep
- {"Inception", {dep, rep}} the inceptenion layers use featInt, last \
one uses featOut. contains rep inception blocks with dep number of \
layers with increasing dilation (1,3,5,7,9)
- {"U2Net", {dep, bool}} an U2Net block that is dep layers deep. If \
bool is set True all layers are dilation instead of scaling. 
"*)


Options[MakeUnet] = {
	DropoutRate -> 0.2, 
	BlockType -> "ResNet", 
	ActivationType -> "GELU",
	NetworkDepth -> 5,

	DownsampleSchedule -> Automatic,
	SettingSchedule -> Automatic,
	FeatureSchedule -> 32,

	MonitorCalc -> False
};

MakeUnet[nClass_, dimIn_, opts : OptionsPattern[]] := MakeUnet[1, nClass, dimIn, opts]

MakeUnet[nChan_, nClass_, dimIn_, OptionsPattern[]] := Block[{
		ndim, depth, nam, drop, blockType, actType, scaling, feature, setting, mon, net
	},

	{drop, blockType, actType, depth, scaling, feature, setting, mon} = 
		OptionValue[{DropoutRate, BlockType, ActivationType, NetworkDepth, DownsampleSchedule, FeatureSchedule, SettingSchedule, MonitorCalc}];

	If[mon, Echo[{blockType, actType}, "Network block type: "]];

	(*is the network 2D or 3D*)
	ndim = Length@dimIn;
	If[mon, Echo[{ndim, dimIn}, "Network dimension order: "]];

	(*define the scaling, automatic is 2 for all layers, integer is for all layers, list is padded with last value*)
	scaling = Which[
		scaling === Automatic, ConstantArray[2, depth],
		IntegerQ[scaling], ConstantArray[scaling, depth],
		ListQ[scaling],	
			If[Length@scaling =!= depth, Message[MakeUnet::scale]];
			PadRight[scaling, depth, Last[scaling]]
	];
	If[mon, Echo[scaling, "Network scaling shedual: "]];

	(*define the setting, can be Automatic, number or {set}, or list of settings*)
	setting = Which[
		setting === Automatic, ConstantArray[blockType /. netDefaults, depth],
		IntegerQ[setting], ConstantArray[setting, depth],
		ListQ[setting], If[Length[setting] == 2 && VectorQ[setting],
			ConstantArray[setting, depth],
			If[Length@setting =!= depth, Message[MakeUnet::sett]];
			PadRight[setting, depth, Last[setting]]
		]
	];
	If[mon, Echo[setting, "Network setting shedual: "]];

	(*define the number of features per layer*)
	feature = Round[Which[
		IntegerQ[feature], feature  Switch[blockType,
			"UNet" | "ResNet", 2^Range[0, depth - 1],
			"DenseNet" | "Inception", 0.5 Range[2, depth + 1],
			"U2Net", ConstantArray[1, depth]
		],
		ListQ[feature], If[Length[feature] == 2 && VectorQ[feature],
			ConstantArray[feature, depth], 
			If[Length@feature =!= depth, Message[MakeUnet::feat]];
			PadRight[feature, depth, Last[feature]]
		]
	]];
	If[VectorQ[feature], feature = Round[Transpose[{feature, feature/2}]]];
	If[mon, Echo[feature, "Network feature shedual: "]];

	(*make the network*)
	nam = #1 <> ToString[#2] &;
	net = NetGraph[
		(*define the blocks*)
		Join[
			(*encoding layers*)
			Table[nam["enc_", i] -> MakeNode[{"Encode", If[i === 1, 1, scaling[[i]]], drop}, {{blockType, setting[[i]]}, feature[[i]], {actType, ndim}}, If[i===1, 0, feature[[i-1, 1]]]], {i, 1, depth}],
			(*deconding layers*)
			Table[nam["dec_", i] -> MakeNode[{"Decode", scaling[[i + 1]], If[i === 1, 0, drop]}, {{blockType, setting[[i]]}, feature[[i]], {actType, ndim}}, feature[[i,1]]], {i, 1, depth - 1}],
			(*start and mapping*)
			{"start" -> NetGraph@NetChain@Conv[Last@Flatten@{First@feature}, {ndim, 1}, actType],
			"map" -> NetGraph@NetChain@ClassMap[ndim, nClass]}
		],

		(*define the connections*)
		Join[
			Flatten@Table[{
				nam["enc_", i] -> nam["enc_", i + 1],
				nam["enc_", i] -> nam["dec_", i],
				If[i === depth - 1, nam["enc_", i + 1], nam["dec_", i + 1]] -> nam["dec_", i]
			}, {i, 1, depth - 1}],
			{"start" -> "enc_1", "dec_1" -> "map"}
		], 
		
		(*network settings and options*)
		"Input" -> Join[{nChan}, dimIn]
	];

	(*Monitor network properties*)
	If[mon, Echo[MBCount[NetInitialize@net], "Network size: "]];
	If[mon, Echo[PrintKernels[net], "Network discription: "]];

	(*return network*)
	net
]


(* ::Subsubsection::Closed:: *)
(*MakeNode*)


MakeNode[nodeType_?StringQ, blockConfig_] := MakeNode[{nodeType, 1, 0}, blockConfig, 0]

MakeNode[nodeType_?StringQ, blockConfig_, chan_?IntegerQ] := MakeNode[{nodeType, 1, 0}, blockConfig, chan]

MakeNode[{nodeType_, scale_}, blockConfig_] := MakeNode[{nodeType, scale, 0}, blockConfig, 0]

MakeNode[{nodeType_, scale_, drop_}, blockConfig_]:= MakeNode[{nodeType, scale, drop}, blockConfig, 0]

MakeNode[{nodeType_, scale_}, blockConfig_, chan__?IntegerQ] := MakeNode[{nodeType, scale, 0}, blockConfig, chan]

MakeNode[{nodeType_, scale_, drop_}, blockConfig_, chan__?IntegerQ] := Block[{
		block, scaleVec, dropL, dropC, scaleL, scaleC, down, up
	},

	(*make the block*)
	block = ConvBlock @@ blockConfig;
	(*chan = Last@Flatten@{blockConfig[[2]]};*)

	(*make correct scale vector*)
	scaleVec = If[IntegerQ[scale], ConstantArray[scale, blockConfig[[3, 2]](*dim*)], scale];

	(*figure out dropout*)
	dropL = If[0 < drop, "drop" -> DropoutLayer[drop], Nothing];
	dropC = If[0 < drop, "block" -> "drop", Nothing];

	NetGraph@NetFlatten[Switch[nodeType, 
		"Encode", 
		scaleL = If[Times @@ scaleVec =!= 1, "scale" -> ConvScale["Encode", scaleVec, chan], Nothing];
		scaleC = If[Times @@ scaleVec =!= 1, "scale" -> "block", Nothing];
		NetGraph[{"block" -> block, scaleL, dropL}, {scaleC, dropC}],

		"Decode",
		scaleL = If[Times @@ scaleVec =!= 1, "scale" -> ConvScale["Decode", scaleVec, chan], Nothing];
		scaleC = If[Times @@ scaleVec =!= 1, NetPort["Input2"] -> "scale", NetPort["Input2"]];
		NetGraph[{"cat" -> CatenateLayer[], "block" -> block, scaleL, dropL}, {{NetPort["Input1"], scaleC} -> "cat" -> "block", dropC}]
	], 1]
]


(* ::Subsubsection::Closed:: *)
(*ConvBlock*)


ConvBlock[block_, feat_?IntegerQ, dim_?IntegerQ] := ConvBlock[block, {feat, Round[feat/2]}, {"None", dim, 1}]

ConvBlock[block_, feat_?IntegerQ, {act_, dim_}] := ConvBlock[block, {feat, Round[feat/2]}, {act, dim, 1}]

ConvBlock[block_, feat_?IntegerQ, {act_, dim_, dil_}] := ConvBlock[block, {feat, Round[feat/2]}, {act, dim, dil}]

ConvBlock[block_, {featOut_, featInt_}, {act_, dim_}] := ConvBlock[block, {featOut, featInt}, {act, dim, 1}]

ConvBlock[block_, {featOut_, featInt_}, {act_, dim_, dil_}] := Block[{
		blockType, type, blockSet, repBlock, dep, rep, nam, dilf, sclf
	},
(*short notation for naming layers*)
nam = #1 <> ToString[#2] &;
(*get the block settings if not defined*)
{blockType, blockSet} = If[Length[block] === 2, block, {block, block /. netDefaults}];

(*swtich between the different block types*)
(*NetGraph@NetFlatten@*)
	Switch[blockType,
		"Conv", 
		dep = blockSet;
		ConvBlock[{"UNet", blockSet}, {featOut, featInt}, {act, dim, dil}],

		"UNet",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"conv" -> Flatten[Table[Conv[If[i === dep, featOut, featOut], {dim, 3, dil}, act], {i, 1, dep}]]}, {}
		],

		"ResNet",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"con" -> Flatten[Table[If[i =!= dep, Conv[featInt, {dim, 3}, act], Conv[featOut, {dim, 3}]], {i, 1, dep}]],
			"skip" -> Conv[featOut, {dim, 1}],
			"tot" -> {TotalLayer[], ActivationLayer[act]}
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
		], 1],

		"Inception",
		{dep, rep} = If[IntegerQ[blockSet], {blockSet, 2}, blockSet];
		repBlock = NetGraph[
			Join[
				Table[nam["dil", 2 (i - 1) + 1] -> Conv[featInt, {dim, 3, 2 (i - 1) + 1}, act], {i, 1, dep}],
				{"cat" -> CatenateLayer[]}
			],
			{NetPort["Input"] -> Table[nam["dil", 2 (i - 1) + 1], {i, 1, dep}] -> "cat"}
		];
		NetGraph@NetFlatten[NetGraph[
			Append[Table[nam["rep", i] -> repBlock, {i, rep}], "out" -> Conv[featOut, {dim, 1}, "None"]],
			Append[Table[nam["rep", i] -> nam["rep", i + 1], {i, rep - 1}], 
			nam["rep", rep] -> "out"]
		], 1],

		"U2Net",
		{dep, type} = If[IntegerQ[blockSet], {blockSet, True}, blockSet];
		dilf = If[! #1, 2^(#2 - 1), #3] &;
		sclf = If[! #1, 1, 2] &;
		NetGraph[
			Join[
				Table[nam["enc_", i] -> MakeNode[{"Encode", If[i == 1, 1, sclf[type]]}, {"Conv", featInt, {act, 3, dilf[type, i, 1]}}, featInt], {i, 1, dep}],
				Table[nam["dec_", i] -> MakeNode[{"Decode", sclf[type]}, {"Conv", If[i === 1, featOut, featInt], {act, 3, dilf[type, i, 1]}}, featInt], {i, 1, dep - 1}],
				{"start" -> Conv[featOut, {dim, 1}, act], "add" -> TotalLayer[]}
			], 
			Join[
				Flatten@Table[{
					nam["enc_", i] -> nam["enc_", i + 1], 
					nam["enc_", i] -> nam["dec_", i],
					If[i === dep - 1, nam["enc_", i + 1], nam["dec_", i + 1]] -> nam["dec_", i]}
				, {i, 1, dep - 1}],
				{"start" -> "enc_1", {"start", "dec_1"} -> "add"}
			]
		]
	]
]



(* ::Subsubsection::Closed:: *)
(*ConvScale*)


ConvScale["Encode", scaleVec_, chan_] := {ConvolutionLayer[chan, scaleVec, "Stride" -> scaleVec]}

ConvScale["Decode", scaleVec_, chan_] := {
	ResizeLayer[Scaled /@ scaleVec, Resampling -> "Nearest"],
	ConvolutionLayer[chan, scaleVec, "Stride" -> 0 scaleVec + 1, 
	"PaddingSize" -> (If[OddQ[#], {(Ceiling[#/2] - 1), (Ceiling[#/2] - 1)}, {0, #/2}] & /@ scaleVec)]
}


(* ::Subsubsection::Closed:: *)
(*Conv*)


Conv[featOut_, {dim_, kern_}] := Conv[featOut, {dim, kern, 1}, "None"]

Conv[featOut_, {dim_, kern_}, act_] := Conv[featOut, {dim, kern, 1}, act]

Conv[featOut_, {dim_, kern_, dil_}] := Conv[featOut, {dim, kern, dil}, "None"]

Conv[featOut_, {dim_, kern_, dil_}, act_] := {
	(*The most basic conv block CON>BN>ACT used in each of the advanced conv blcoks*)
	ConvolutionLayer[featOut, ConstantArray[kern, dim], 
		"PaddingSize" -> (Ceiling[kern/2] - 1) dil, 
		"Stride" -> 1, 
		"Dilation" -> dil],	
	BatchNormalizationLayer[], 
	ActivationLayer[act]
}


(* ::Subsubsection::Closed:: *)
(*ClassMap*)


ClassMap[dim_, nClass_] := Flatten[{ConvolutionLayer[nClass, 1], 
	If[nClass > 1, 
		{TransposeLayer[Which[dim === 2, {3, 1, 2}, dim === 3, {4, 1, 2, 3}]], SoftmaxLayer[]}, 
		{LogisticSigmoid, FlattenLayer[1]}
	]
}]


(* ::Subsubsection::Closed:: *)
(*ActivationLayer*)


ActivationLayer[] := ActivationLayer["GELU"]

ActivationLayer[actType_] := Switch[actType, 
	"LeakyRELU", ParametricRampLayer[], 
	"None" | "", Nothing, 
	_, ElementwiseLayer[actType]
]


(* ::Subsection::Closed:: *)
(*PrintKernels*)


PrintKernels[net_] := Block[{convs, kerns, count, pars},
	convs = Select[Information[net, "LayersList"], Head[#] === ConvolutionLayer &];
	kerns = Information[#, "ArraysDimensions"][{"Weights"}] & /@ convs;
	count = Sort[{Length[#], #[[1, 1]], Total[#[[All, 2]]], Total[#[[All, 3]]]} & /@ 
		GatherBy[{#[[3 ;;]], Times @@ #[[1 ;; 2]], Times @@ #[[1 ;;]]} & /@ 
			kerns, First]];
	pars = Information[net, "ArraysTotalElementCount"];
	Column[{Grid[
		Transpose[{Style[#, Bold] & /@ {"Convolution Layers", "Total Kernels", "Total Weighths"},
		{Length[convs], Total[count[[All, 3]]], 
		Total[count[[All, 4]]]}}], Alignment -> Left, Spacings -> {2, 1}],
		"",
		Grid[Join[{Style[#, Bold] & /@ {"Count", "Size", "Kernels", "Weights"}}, 
		count], Alignment -> Left, Spacings -> {2, 1}]}, 
	Alignment -> Center]
]


(* ::Subsection:: *)
(*LossLayers*)


(* ::Subsubsection::Closed:: *)
(*AddLossLayer*)


SyntaxInformation[AddLossLayer] = {"ArgumentsPattern" -> {_}};

AddLossLayer[net_]:=Block[{dim},
	(*http://arxiv.org/abs/2312.05391*)
	dim = Length[Information[net,"OutputPorts"][[1]]]-1;
	NetGraph[<|
		"net"->net,
		"SoftDice" -> SoftDiceLossLayer[dim, 2],
		"Jaccard" -> SofJaccardLossLayer[dim],
		"Tversky" -> TverskyLossLayer[dim, 0.7],
		"SquaredDiff" -> {MeanSquaredLossLayer[], ElementwiseLayer[50 #&]},
		"CrossEntropy" -> {CrossEntropyLossLayer["Probabilities"]}
	|>,{
		{"net", NetPort["Target"]}->"SoftDice"->NetPort["SoftDice"],(*using squared dice*)
		{"net", NetPort["Target"]}->"Jaccard"->NetPort["Jaccard"],
		{"net", NetPort["Target"]}->"Tversky"->NetPort["Tversky"],
		{"net", NetPort["Target"]}->"SquaredDiff"->NetPort["SquaredDiff"],(*Brier Score*)
		{"net", NetPort["Target"]}->"CrossEntropy"->NetPort["CrossEntropy"]
	}]
]


(* ::Subsubsection::Closed:: *)
(*SoftDiceLossLayer*)


SyntaxInformation[SoftDiceLossLayer] = {"ArgumentsPattern" -> {_, _.}};

SoftDiceLossLayer[dim_ ] := SoftDiceLossLayer[dim, 1]

SoftDiceLossLayer[dim_, n_] := Block[{smooth},
	(*10.48550/arXiv.1911.02855 and 10.48550/arXiv.1606.04797 for scquared dice loss look at v-net*)
	smooth =1;
	NetGraph[<|
		(*flatten input and target; function layer allows to switch to L2 norm if #^2*)
		"input" -> {FunctionLayer[#^n &], AggregationLayer[Total, ;; dim]},
		"target" -> {FunctionLayer[#^n &], AggregationLayer[Total, ;; dim]},
		(*intersection or TP*)
		"intersection" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; dim]},
		(*the loss function 2*intersection / (input + target)*)
		"dice" -> {ThreadingLayer[1. - ((2. #1 + smooth) / (#2 + #3 + smooth)) &], AggregationLayer[Mean, 1]}
	|>, {
		NetPort["Input"] -> "input",
		NetPort["Target"] -> "target",
		{NetPort["Target"], NetPort["Input"]} -> "intersection",
		{"intersection",  "target", "input"} -> "dice" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*SofJaccardLossLayer*)


SyntaxInformation[SofJaccardLossLayer] = {"ArgumentsPattern" -> {_}};

SofJaccardLossLayer[dim_ ] := SofJaccardLossLayer[dim, 1]

SofJaccardLossLayer[dim_, n_]:= Block[{smooth},
	smooth = 1;
	NetGraph[<|
		(*flatten input and target; function layer allows to switch to L2 norm if #^2*)
		"input" -> {FunctionLayer[#^n &], AggregationLayer[Total, ;; dim]},
		"target" -> {FunctionLayer[#^n &], AggregationLayer[Total, ;; dim]},
		(*intersection or TP*)
		"intersection" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; dim]},
		(*the loss function intersection / union with union = (input + target - intersection)*)
		"Jaccard" -> {ThreadingLayer[1. - ((#1 + smooth) / ((#2 + #3) - #1 + smooth)) &], AggregationLayer[Mean, 1]}
	|>, {
		NetPort["Input"] -> "input",
		NetPort["Target"] -> "target",
		{NetPort["Target"], NetPort["Input"]} -> "intersection",
		{"intersection", "target", "input"} -> "Jaccard" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*TwerskyLossLayer*)


SyntaxInformation[TwerskyLossLayer] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

TverskyLossLayer[dim_] := TverskyLossLayer[dim, 0.7]

TverskyLossLayer[dim_, beta_?NumberQ] := Block[{smooth, alpha},
	smooth = 1;
	alpha = 1- beta;
	(* https://doi.org/10.48550/arXiv.1706.05721 *)
	NetGraph[<|
		(*intersection or TP*)
		"truePos" -> {ThreadingLayer[Times], AggregationLayer[Total, ;; dim]},
		"falsePos" -> {ThreadingLayer[(1 - #1) #2 &], AggregationLayer[Total, ;; dim]},
		"falseNeg" -> {ThreadingLayer[#1 (1 - #2) &], AggregationLayer[Total, ;; dim]},
		(*the loss function TP / (TP + a FP + b FN)*)
		"Twersky" -> {ThreadingLayer[1. - (#1 + smooth) / (#1 + alpha #2 + beta #3 + smooth) &], AggregationLayer[Mean, 1]}
	|>, {
		{NetPort["Target"], NetPort["Input"]} -> "truePos",
		{NetPort["Target"], NetPort["Input"]} -> "falsePos",
		{NetPort["Target"], NetPort["Input"]} -> "falseNeg",
		{"truePos", "falsePos", "falseNeg"} -> "Twersky" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsection:: *)
(*Encoders*)


(* ::Subsubsection::Closed:: *)
(*ClassEndocer*)


SyntaxInformation[ClassEncoder] = {"ArgumentsPattern" -> {_, _.}};

ClassEncoder[data_]:= ClassEncoder[data, Round[Max@data]]

ClassEncoder[data_, nClass_]:= ToPackedArray@Round@If[nClass === 1, data, ClassEncoderC[data, nClass]]

ClassEncoderC = Compile[{{data, _Integer, 2}, {n, _Integer, 0}},
	Transpose[1 - Unitize[ConstantArray[data, n] - Range[n]], {3, 1, 2}]
, RuntimeAttributes -> {Listable}]


(* ::Subsubsection::Closed:: *)
(*ClassDecoder*)


SyntaxInformation[ClassDecoder] = {"ArgumentsPattern" -> {_, _.}};

ClassDecoder[data_]:= ClassDecoder[data, Last@Dimensions@data]

ClassDecoder[data_, nClass_]:=ToPackedArray@Round@ClassDecoderC[data, Range[nClass]]

ClassDecoderC = Compile[{{data, _Real, 1}, {nl, _Integer, 1}}, Block[{cl},
	cl = (1 - Unitize[Chop[(data/Max[data]) - 1]]);
	If[Total[cl] > 1, 1, Total[nl cl]]
], RuntimeAttributes -> {Listable}]


(* ::Subsection:: *)
(*MakeClassify*)


(* ::Subsubsection::Closed:: *)
(*MakeClassifyNetwork*)


Options[MakeClassifyNetwork]={ImageSize->{128,128}};

MakeClassifyNetwork[classes_,OptionsPattern[]]:=Block[{enc, dec, net,imSize},
	imSize = OptionValue[ImageSize];
	enc = NetEncoder[{"Class",classes,"IndicatorVector"}];
	dec = NetDecoder[{"Class",classes}];
	net = NetChain[{
		ConvolutionLayer[16,7,"Stride"->1,PaddingSize->3],BatchNormalizationLayer[],ElementwiseLayer["GELU"],PoolingLayer[4,4],
		ConvolutionLayer[32,5,"Stride"->1,PaddingSize->2],BatchNormalizationLayer[],ElementwiseLayer["GELU"],PoolingLayer[4,4],
		ConvolutionLayer[64,3,"Stride"->1,PaddingSize->1],BatchNormalizationLayer[],ElementwiseLayer["GELU"],PoolingLayer[4,4],
		FlattenLayer[],LinearLayer[128],BatchNormalizationLayer[],ElementwiseLayer["GELU"],LinearLayer[64],
		BatchNormalizationLayer[],ElementwiseLayer["GELU"],LinearLayer[Length@classes],SoftmaxLayer[]
		},"Input"->Prepend[imSize,1]
	];
	NetFlatten@NetChain[{net},"Input"->NetEncoder[{"Image",imSize,ColorSpace->"Grayscale"}],"Output"->dec]
]


(* ::Subsubsection::Closed:: *)
(*MakeClassifyImage*)


Options[MakeClassifyImage]={ImageSize->{128,128}};

MakeClassifyImage[dat_, opts:OptionsPattern[]]:=Switch[ArrayDepth[dat],
	2, MakeClassifyImage[dat, opts],
	3, MakeClassifyImage[#, opts]&/@dat,
	_, $Failed
]

MakeClassifyImage[dat_?MatrixQ, OptionsPattern[]]:=Block[{imSize},
	imSize=OptionValue[ImageSize];
	If[Total[Flatten[dat]]<10,
		Image@ConstantArray[0.,imSize],
		ImageResize[Image[Rescale[First[AutoCropData[{dat},CropPadding->0][[1]]]]],imSize]
	]
];


(* ::Subsection:: *)
(*ClassifyData*)


(* ::Subsubsection::Closed:: *)
(*ClassifyData*)


ClassifyData[dat_, met_]:=Block[{data,len,ran,datF,kneeStart,kneeEnd,pos},
	data = MakeClassifyImage[dat];
	Which[
		FileExistsQ[met]&&FileExtension[met]==="wlnet", Import[met][data],
		StringQ[met],Switch[met,
			"LegSide",FindLegSide[data],
			"LegPosition",FindLegPos[data]],
		Head[met]===NetChain||Head[met]===NetGraph,met[data]
	]
]


(* ::Subsubsection::Closed:: *)
(*FindLegSide*)


FindLegSide[data_]:=Block[{net,imSize},
	net=GetNeuralNet["LegSide"];
	If[net===$Failed,$Failed,
		imSize=NetDimensions[NetReplacePart[net,"Input"->None],"Input"][[2;;]];
		If[!AllTrue[ImageDimensions/@data,#===imSize&],$Failed,
			Last@Keys@Sort@Counts[net[data]]
]]]


(* ::Subsubsection::Closed:: *)
(*FindLegPos*)


FindLegPos[data_]:=Block[{net,len,ran,datF,kneeStart,kneeEnd,pos,imSize},
	net=GetNeuralNet["LegPosition"];
	If[net===$Failed,$Failed,
		imSize=NetDimensions[NetReplacePart[net,"Input"->None],"Input"][[2;;]];
		If[!AllTrue[ImageDimensions/@data,#===imSize&],$Failed,
			len=Length[data];
			ran=Range[1,len];
			(*find loc per slice*)
			datF=MedianFilter[net[data]/.Thread[{"Lower","Knee","Upper"}->{1.,2.,3.}],1];
			{kneeStart,kneeEnd}=First[SortBy[Flatten[Table[{a,b,PosFunc[a,b,len,ran,datF]},{a,0,len},{b,a+1,len}],1],Last]][[1;;2]];
			pos=Which[kneeStart==0.&&kneeEnd==len,"Knee",kneeStart>0&&kneeEnd>=len,"Lower",kneeStart==0.&&kneeEnd=!=len,"Upper",kneeStart=!=0.&&kneeEnd=!=len,"Both"];
			{pos,{kneeStart+1,kneeEnd}}
]]]


PosFunc=Compile[{{a,_Integer,0},{b,_Integer,0},{l,_Integer,0},{x,_Integer,1},{d,_Real,1}},Total[((Which[1<=#<=a,1.,a<=#<=b,2.,b<=#<=l,3.,True,0]&/@x)-d)^2]];


(* ::Subsection::Closed:: *)
(*SegmentDataGUI*)


SegmentDataGUI[] := DynamicModule[{inputFile, outputFile}, Block[{dat, vox, seg, status, diag, option},
	NotebookClose[segwindow];

	option = "Legs";
	
	diag = DialogNotebook[
		status = TextCell@"";
		{
			TextCell["Please enter the paths for the input and output files:"],

			Grid[{{
				TextCell["Status: "], Dynamic[status]
			}, {
				TextCell["Input File: "],
				InputField[Dynamic[inputFile], String, 
				FieldHint -> "Enter input file path", FieldSize -> {25, 1}],
				Button["Browse", inputFile = SystemDialogInput["FileOpen"], Method -> "Queued"]
			}, {
				TextCell["Output File: "], 
				InputField[Dynamic[outputFile], String, 
				FieldHint -> "Enter output file path", 
				FieldSize -> {25, 1}],
				Button["Browse", outputFile = SystemDialogInput["FileSave"], Method -> "Queued"]
			}, {
				TextCell["Segmentation type"], 
				PopupMenu[Dynamic[option], {"Legs", "LegBones"}]
			},{
				Button["Start Segmentation, please be patient", 
					If[! NiiFileExistQ[inputFile], 
						MessageDialog["Input file could not be foud."]
						,
						status = TextCell@"Importing";
						{dat, vox} = ImportNii[inputFile];
						status = TextCell@"Segmenting Data";
						seg = SegmentData[dat, option, TargetDevice -> "CPU"];
						status = TextCell@"Exporting";

						CopyFile[GetAssetLocation["MusclesLegLabels"], 
						ConvertExtension[outputFile, ".txt"], 
						OverwriteTarget -> True];
						ExportNii[seg, vox, outputFile];
						status = Button["Go to " <> FileBaseName@outputFile, 
							SystemOpen[DirectoryName@outputFile]];
					],
				Method -> "Queued"]
			}}, Alignment -> Left],
			Row[{
				DefaultButton[],
				CancelButton[]
			}]
		}
	];

	segwindow = CreateWindow[diag, WindowTitle -> "Muscle segmentation", WindowSize -> All];
];];


(* ::Subsection:: *)
(*SegmentData*)


(* ::Subsubsection::Closed:: *)
(*SegmentData*)


Options[SegmentData] = {TargetDevice -> "GPU", MaxPatchSize->Automatic, Monitor->False};

SegmentData[data_, what_, OptionsPattern[]] := Block[{
		dev, max, mon, patch, pts, dim ,loc, set, net, type, segs, all
	},

	{dev, max, mon} = OptionValue[{TargetDevice, MaxPatchSize, Monitor}];

	(*split the data in upper and lower legs and left and right*)
	If[mon, Echo[Dimensions@data, "Analyzing the data with dimensions:"]];
	{{patch, pts, dim}, loc, set} = SplitDataForSegementation[data];
	If[mon, Echo[Thread[{loc,Dimensions/@ patch}], "Segmenting "<>what<>" locations with dimenisons:"]];

	(*get the network name and data type*)
	{net, type} = Switch[what,
		"LegBones", {"SegLegBones"&, "Bones"},
		"Legs",	{(#[[1]] /. {"Upper" -> "SegThighMuscle", "Lower" -> "SegLegMuscle"})&, "Muscle"},
		_, Return[]
	];

	(*Perform the segmentation*)
	segs = MapThread[(
		If[mon, Echo[{#2, net[#2]}, "Performing segmentation for: "]];
		segs = ApplySegmentationNetwork[#1, net[#2], TargetDevice -> dev, MaxPatchSize->max, Monitor->mon];
		ReplaceLabels[segs, #2, type]
	) &, {patch, loc}];

	(*Merge all segmentations for all expected labels*)
	all = Select[DeleteDuplicates[Sort[Flatten[GetSegmentationLabels/@segs]]], IntegerQ];
	If[mon, Echo[all, "Putting togeteher the segmenations with lables"]];
	(*after this only one cluster per label remains*)
	PatchesToData[segs, pts, dim, all]
]


(* ::Subsubsection::Closed:: *)
(*ReplaceLabels*)


ReplaceLabels[seg_, loc_, type_] := Block[{what, side, labNam, labIn, labOut, file},
	{what, side} = loc;
	labIn = GetSegmentationLabels[seg];

	(*for now overwrite the labIn with custom values since some muscles are not segemnted *)
	file = GetAssetLocation@Switch[type,
		"Muscle", Switch[what, 
			"Upper", "MusclesLegUpperLabels",
			"Lower", "MusclesLegLowerLabels"],
		"Bones", "BonesLegLabels"];
	labNam = # <> "_" <> side & /@ MuscleLabelToName[labIn, file];
	labOut = MuscleNameToLabel[labNam, GetAssetLocation["MusclesLegLabels"]];

	ReplaceSegmentations[seg, labIn, labOut]
]


(* ::Subsubsection::Closed:: *)
(*ApplySegmentationNetwork*)


Options[ApplySegmentationNetwork]={TargetDevice->"GPU", DataPadding->0, MaxPatchSize->Automatic, Monitor->False}

SyntaxInformation[ApplySegmentationNetwork] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ApplySegmentationNetwork[dat_, netI_, OptionsPattern[]]:=Block[{
		dev, pad , lim, data, crp, net, dim, ptch, 
		patch, pts, seg, mon, dimi, dimc, lab, class
	},
	
	{dev, pad, lim, mon} = OptionValue[{TargetDevice, DataPadding, MaxPatchSize, Monitor}];
	If[lim === Automatic, lim = If[dev==="GPU", 224, 224]];

	data=If[First@Dimensions@dat===1, First@dat, dat];

	dimi = Dimensions@data;
	{data, crp} = AutoCropData[data, CropPadding->0];
	dimc = Dimensions@data;
	data = N@ArrayPad[data, pad, 0.];
	dim = Dimensions[data];

	If[mon, Echo[{dimi, dimc, dim}, "Data dimensions before and after cropping and padding are: "]];

	net = If[StringQ[netI],	GetNeuralNet[netI],	netI];

	If[net===$Failed, $Failed,
		(*calculate the patch size for the data *)
		ptch = FindPatchDim[net, dim, lim];

		(*create the patches*)
		{patch, pts} = DataToPatches[data, ptch, PatchNumber -> 0, PatchPadding->pad];
		If[mon, Echo[{ptch, Length@patch}, "Patch size and created number of patches is:"]];

		(*actualy perform the segmentation with the NN*)
		net = ChangeNetDimensions[net, "Dimensions" ->ptch];
		seg = ToPackedArray[Round[ClassDecoder[net[{NormDatH[#]}, TargetDevice->dev]]&/@patch]];
		
		(*reverse all the padding and cropping and merged the patches if needed*)
		class = NetDimensions[net,"Output"][[-1]];
		If[mon, Echo[{Dimensions[seg], Sort@Round@DeleteDuplicates[Flatten[seg]]}, 
			"Segmentations dimensions and labels:"]];
		seg = ArrayPad[PatchesToData[ArrayPad[#, -pad] & /@ seg, Map[# + {pad, -pad} &, pts, {2}], 
			dim, Range[class]], -pad];
		seg = Ramp[seg - 1]; (*set background to zero*)
		seg = ToPackedArray@Round@ReverseCrop[seg, dimi, crp];
		If[mon, Echo[{Dimensions[seg], Sort@Round@DeleteDuplicates[Flatten[seg]]}, 
			"Output segmentations dimensions and labels:"]];
		
		(*give the output*)
		seg
	]
]


(* ::Subsubsection::Closed:: *)
(*NormDat*)


NormDat[dat_] := Block[{q = Quantile[Flatten[dat], 0.9], m = Max[dat]}, If[q <= 0.5 m, If[m===0., dat, dat/m], If[q===0., dat, 0.75 dat/q]]]


NormDatH = Compile[{{dat, _Real, 3}}, Block[{fl, flp, min, max, bins, cdf, n},
	n = 512;
	fl = Flatten[dat];
	flp = Pick[fl, Unitize[fl], 1];
	{min, max} = MinMax[flp];
	bins = BinCounts[flp, {min, max, (max - min)/n}];
	cdf = Prepend[N[Accumulate[bins]/Total[bins]], 0.];
	Map[cdf[[# + 1]] &, Clip[Floor[(dat - min)/(max - min)  n], {0, n}, {0, n}], {-2}]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True]


(* ::Subsubsection::Closed:: *)
(*FindPatchDim*)


FindPatchDim[net_, dims_] := FindPatchDim[net, dims, 1000]

FindPatchDim[net_, dims_, lim_] := Block[{
		dim, inp, out, class, sc, ptch, cont, u, dimM, dimN
	},

	(*if dim is data array then get the dimensions, 
	if it is vector it is already dim*)
	dim = If[! VectorQ[dims], Dimensions@dims, dims];

	(*get net properties*)
	inp = Rest[NetDimensions[net, "Input"]];
	out = Rest[NetDimensions[net, "LastEncoding"]];
	class = NetDimensions[net, "Output"][[-1]];
	sc = inp/out;

	(*initialize dim*)
	dimM = sc  Ceiling[dim/sc];

	(*if memeory limit is given find patch that fits*)
	If[lim =!= 1000,
		If[CubeRoot[N[Times @@ dimM]] < lim, dimN = dimM, u = 1;
		cont = True;
		While[cont, u++;
			dimN = u  sc;
			dimN = Min /@ Transpose[{dimN, dimM}];
			cont = CubeRoot[Times @@ dimN] < lim && Min [dimN] < Min[dim]]],
		dimN = dimM
	];

	(*output the patch dim*)
	Max /@ Thread[{dimN, inp}]
]


(* ::Subsubsection::Closed:: *)
(*NetDimensions*)


NetDimensions[net_]:=NetDimensions[net, ""]

NetDimensions[net_, port_]:=Switch[port,
	"Input", Information[net,"InputPorts"]["Input"],
	"Output", Information[net,"OutputPorts"]["Output"],
	"FirstEncoding", Information[NetTake[net, First[Select[Keys[net[[All, 1]]], StringContainsQ[#, "enc_"] &]]], "OutputPorts"]["Output"],
	"LastEncoding", Information[NetTake[net,Last[Select[Keys[net[[All,1]]],StringContainsQ[#,"enc_"]&]]],"OutputPorts"]["Output"],
	_, {First@NetDimensions[net, "Input"], Last@NetDimensions[net, "Output"], Rest@NetDimensions[net, "Input"], First@NetDimensions[net, "FirstEncoding"]}
] 


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
		filt = First@NetDimensions[netOut, "FirstEncoding"];
		start = NetFlatten[NetTake[netOut, "start"]];
		netOut = NetReplacePart[netOut, {"Input" -> Prepend[Round@dimOut, Round@nChanOut], "start" -> UNetStart[filt, Round@nChanOut, Round@dimOut, start[[-1]]]}]
	];

	(*Change output Classes if needed*)
	If[nClassOut =!= None && nClassOut=!=nClassIn,
		netOut = NetReplacePart[netOut, {"map" -> UNetMap[Round@dimOut, Round@nClassOut], "Output" -> Round@Append[dimOut, nClassOut]}]
	];
	netOut
]


(* ::Subsubsection::Closed:: *)
(*GetNeuralNet*)


SyntaxInformation[GetNeuralNet] = {"ArgumentsPattern" -> {_}};

GetNeuralNet[name_?StringQ]:= GetNeuralNetI[name]

GetNeuralNetI[name_]:=GetNeuralNetI[name]=Which[
	FileExistsQ[name], Import[name],
	FileExistsQ[GetAssetLocation[name]], Import[GetAssetLocation[name]],
	True, $Failed
]


(* ::Subsection:: *)
(*Data Patching*)


(* ::Subsubsection::Closed:: *)
(*PatchesToData*)


SyntaxInformation[PatchesToData] = {"ArgumentsPattern" -> {_, _, _., _.}};

PatchesToData[patches_, ran_] := PatchesToData[patches, ran, Max /@ Transpose[ran[[All, All, 2]]], {}]

PatchesToData[patches_, ran_, dim : {_?IntegerQ, _?IntegerQ, _?IntegerQ}] := PatchesToData[patches, ran, dim, {}]

PatchesToData[patches_, ran_, dim : {_?IntegerQ, _?IntegerQ, _?IntegerQ}, labs_?ListQ] := Block[{
		sel, zero, a1, a2, b1, b2, c1, c2, dat, wt,
		seg, lab, pos, si, pi, overlap
	},
	zero = SparseArray[{}, dim];

	If[labs === {},
		(*no labs given, assuming normal data*)
		{dat, wt} = Total /@ Transpose@MapThread[(
			dat = zero;
			{{a1, a2}, {b1, b2}, {c1, c2}} = #2;
			dat[[a1 ;; a2, b1 ;; b2, c1 ;; c2]] = #1[[1 ;; a2 - a1 + 1, 1 ;; b2 - b1 + 1, 1 ;; c2 - c1 + 1]];
			{dat, SparseArray[Unitize[dat]]}
		) &, {patches, ran}];
		SparseArray[dat["ExplicitPositions"] -> dat["ExplicitValues"] / wt["ExplicitValues"], dim]
		,
		(*labs given, assuming segmentations*)
		{seg, lab} = Transpose[SplitSegmentations /@ patches];
		seg = Transpose /@ seg;

		seg = (
			pos = Position[lab, #];
			If[pos === {},
				zero,
				si = seg[[##]] & @@ # & /@ pos;
				pi = ran[[#[[1]]]] & /@ pos;
				Ceiling[PatchesToData[si, pi, dim]]
			]
		) & /@ labs;

		(*set overlapping to zero and then add to background label*)
		seg = TakeLargestComponent[#, True] &/@ seg;
		overlap = SparseArray[1 - UnitStep[Total[seg] - 2]];
		If[Min[overlap]===1, seg = TakeLargestComponent[overlap #] &/@ seg];
		MergeSegmentations[Transpose[seg], labs]
	]
]


TakeLargestComponent[seg_]:=TakeLargestComponent[seg, False]

TakeLargestComponent[seg_, err_] := Block[{dim, segc, cr},
	If[Max[seg] < 1,
		seg,
		dim = Dimensions[seg];
		{segc, cr} = AutoCropData[seg];
		segc = ToPackedArray@N@If[err,
			ImageData[Dilation[SelectComponents[Erosion[Image3D[NumericArray[segc, "Integer8"]], CrossMatrix[{1, 1, 1}]], "Count", -1, CornerNeighbors -> False], CrossMatrix[{1, 1, 1}]]],
			ImageData[SelectComponents[Image3D[NumericArray[segc, "Integer8"]], "Count", -1, CornerNeighbors -> False]]
		];
		Round@SparseArray@ReverseCrop[segc, dim, cr]
	]
]



(* ::Subsubsection::Closed:: *)
(*DataToPatches*)


Options[DataToPatches] = {PatchNumber->0, PatchPadding->0}

SyntaxInformation[DataToPatches] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, opts:OptionsPattern[]]:=DataToPatches[dat, patch, "All", opts] 

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}]:={GetPatch[dat, patch, pts], pts}

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, nPatch_, OptionsPattern[]]:=Block[{
		ptch, pts, nRan, pad
	},
	{nRan, pad} = OptionValue[{PatchNumber, PatchPadding}];
	If[Or @@ (#<=2 pad &/@ patch),
		$Failed,
		pts = GetPatchRanges[dat, patch, If[IntegerQ[nPatch], nPatch, "All"], {nRan, pad}];
		ptch = GetPatch[dat, patch, pts];
		{ptch, pts}
	]
] 


(* ::Subsubsection::Closed:: *)
(*GetPatch*)


GetPatch[dat_, pts:{{{_,_},{_,_},{_,_}}..}]:=GetPatch[dat, #]&/@pts

GetPatch[dat_, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}]:=dat[[i1;;i2,j1;;j2,k1;;k2]]

GetPatch[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}]:=GetPatch[dat, patch, #]&/@pts

GetPatch[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}]:=PadRight[dat[[i1;;i2,j1;;j2,k1;;k2]], patch, 0.]


(* ::Subsubsection::Closed:: *)
(*GetPatchRanges*)


GetPatchRanges[dat_, patch_, nPatch_, {nRan_, pad_}]:=Block[{pts},
	pts = GetPatchRangeI[Dimensions[dat], patch, {nRan, pad}];
	If[nPatch==="All", pts, RandomSample[pts, Min[{Length@pts, nPatch}]]]
]


GetPatchRangeI[datDim_?ListQ, patchDim_?ListQ, {nr_, pad_}]:=Tuples@MapThread[GetPatchRangeI[#1,#2, {nr, pad}]&, {datDim, patchDim}]

GetPatchRangeI[dim_?IntegerQ, patch_?IntegerQ, {nr_, pad_}]:=Block[{i,st},
	i = Ceiling[(dim - 2 pad)/(patch - 2 pad)] + nr;
	If[!(dim > patch && i > 1),
		{{1,dim}},
		st = Round[Range[0, 1, 1./(i - 1)](dim - patch)];
		Thread[{st + 1, st + patch}]
	]
]


(* ::Subsection:: *)
(*TrainSegmentationNetwork*)


(* ::Subsubsection::Closed:: *)
(*TrainSegmentationNetwork*)


Options[TrainSegmentationNetwork] = {
	PatchSize -> {32, 96, 96},
	DownsampleSchedule -> {{1, 1, 1},{1, 2, 2}, {2, 2, 2} , {1, 2, 2}, {2, 2, 2}},
	SettingSchedule -> Automatic,
	FeatureSchedule -> 32,
	NetworkDepth -> 5,
	DropoutRate -> 0.2,
	BlockType -> "ResNet",

	BatchSize -> 4,
	RoundLength -> 256,
	MaxTrainingRounds -> 500,
	AugmentData -> True,
	PatchesPerSet -> 1,
	LoadTrainingData -> True,
	LossFunction -> All
};

TrainSegmentationNetwork[{inFol_?StringQ, outFol_?StringQ}, opts : OptionsPattern[]] := TrainSegmentationNetwork[{inFol, outFol}, "Start", opts]

TrainSegmentationNetwork[{inFol_?StringQ, outFol_?StringQ}, netCont_, opts : OptionsPattern[]] := Block[{
		netOpts, batch, roundLength, rounds, data, dDim, nChan, nClass, outName, ittName, makeIm,
		patch, augment, netIn, ittTrain, testData, testVox, testSeg, im, patches,
		netName, monitorFunction, netMon, netOut, trained, 
		validation, files, loss
	},

	(*------------ Get all the configuration struff -----------------*)

	(*getting all the options*)
	netOpts = Join[FilterRules[{opts}, Options@MakeUnet], FilterRules[Options@TrainSegmentationNetwork, Options@MakeUnet]];
	{batch, roundLength, rounds, augment, patches, loss} = OptionValue[{BatchSize, RoundLength, MaxTrainingRounds, AugmentData, PatchesPerSet, LossFunction}];
	
	(*import all the train data*)
	files = FileNames["*.wxf", inFol];

	(*figure out network properties from train data*)
	(*background is 0 but for network its 1 so class +1*)
	data = Import@First@files;
	dDim = Dimensions@data[[1]][[1]];
	nChan = If[Length@dDim === 2, 1, dDim[[2]]];
	nClass = Round[Max@data[[2]] + 1];
	patch = OptionValue[PatchSize];

	(*define the training loss funciton*)
	loss = Which[
		loss === All, {"SoftDice", "SquaredDiff", "Tversky" , "CrossEntropy", "Jaccard"},
		StringQ[loss], {loss}
	];
	If[! And @@ (MemberQ[{"SoftDice", "SquaredDiff", "Tversky", "CrossEntropy", "Jaccard"}, #] & /@ loss), 
		Return[Message[TrainSegmentationNetwork::loss]; $Failed]];

	(*------------ Define the network -----------------*)

	(*make or import network, netCont can be a network or a previous train folder*)
	{netIn, ittTrain} = Which[
		(*start with a clean network*)
		netCont === "Start",
		{MakeUnet[nChan, nClass, patch, Sequence@netOpts], 0}
		,
		(*continue with given network*)
		Head@netCont === NetGraph,
		{netCont , 0}
		,
		(*string can be different things*)
		StringQ[netCont],
		Which[
			(*is wlnet import and start*)
			FileExtension[netCont] === "wlnet",
			If[FileExistsQ[netCont], 
				{Import[netCont], 0}, 
				Return[$Failed]
			]
			,
			(*is an output directory from previous training look for net*)
			DirectoryQ[netCont],
			netIn = FileNames["*_itt_*.wlnet", netCont];
			If[Length[netIn] > 0, 
				{Import@Last@netIn, ToExpression@Last@StringSplit[FileBaseName@Last@netIn, "_"]},
				Return[
					Message[TrainSegmentationNetwork::net];
					Message[TrainSegmentationNetwork::cont]; 
					$Failed]
			]
			,
			True, Return[
				Message[TrainSegmentationNetwork::net]; 
				Message[TrainSegmentationNetwork::inp]; 
				$Failed]
		],
		True, Return[Message[TrainSegmentationNetwork::net]; $Failed]
	];


	(*------------ Define the network -----------------*)


	If[rounds - ittTrain < 5,
		Return[Message[TrainSegmentationNetwork::itt]; $Failed]
		,
		
		(*if the network already exists make the dimensions, classes en channels match the input*)
		netIn = NetInitialize@ChangeNetDimensions[netIn, "Dimensions" -> patch, "Channels" -> nChan, "Classes" -> nClass];
		Echo[netIn, "Network for training"];

		(*---------- Stuff for monitoring ----------------*)
		
		(*Functions for consistant file names*)
		outName = FileNameJoin[{outFol, Last[FileNameSplit[outFol]] <> "_" <> #}]&;
		netName = outName["itt_" <> StringPadLeft[ToString[#], 4, "0"] <> ".wlnet"]&;
		ittName = FileNameJoin[{outFol, "itt_" <> StringPadLeft[ToString[#], 4, "0"] <> ".png"}]&;
		
		(*Make and export test data*)
		{testData, testVox} = MakeTestData[data, 2, patch];
		ExportNii[First@testData, testVox, outName["testSet.nii"]];

		(*segment test data *)
		testSeg = Ramp[ClassDecoder[netIn[testData, TargetDevice -> "GPU"]] - 1];
		ExportNii[testSeg, testVox, outName["testSeg.nii"]];
		
		(*make and export image*)
		makeIm[data_, label_] := ImageAssemble[Partition[
			MakeChannelClassImage[data[[All, #]], label[[#]], {0, nClass - 1}
				] & /@ (Round[Range[2, Length[label] - 1, (Length[label] - 2) / 9.]])
		, 3]];
		im = makeIm[testData, testSeg];
		Export[ittName[ittTrain], im];

		(*Print progress function*)
		Echo[Dynamic[Column[
			{Style["Training Round: " <> ToString[ittTrain], Bold, Large], Image[im, ImageSize->400]}
		, Alignment -> Center]], "Progress"];

		(*define the monitor function, exports image and last net and Nii of result*)
		monitorFunction = (
			ittTrain++;
			(*perform segmentation and export*)
			netMon = NetExtract[#Net, "net"];
			testSeg = Ramp[ClassDecoder[netMon[testData, TargetDevice -> "GPU"]] - 1];
			ExportNii[testSeg, testVox, outName["testSeg.nii"]];
			(*make and export test image*)
			im = makeIm[testData, testSeg];
			Export[ittName[ittTrain], im];
			(*export the network and delete the one from the last itteration*)
			Export[netName[ittTrain], netMon];
			Quiet@DeleteFile[netName[ittTrain - 1]];
		)&;


		(*---------- Train the network ----------------*)

		Echo[DateString[], "Starting training"];

		(*import all train data or train out of memory*)
		data = If[OptionValue[LoadTrainingData]===True, Import/@files, files];
		
		(*prepare a validation set*)
		validation = GetTrainData[data, Round[0.2 roundLength], patch, nClass, AugmentData->augment];
		
		(*train the network*)
		trained = NetTrain[
			AddLossLayer@netIn, 
			{GetTrainData[data, #BatchSize, patch, nClass, AugmentData -> augment, PatchesPerSet->patches] &, 
				"RoundLength" -> roundLength}, 
			All, ValidationSet -> validation,
			
			LossFunction -> {"SoftDice", "SquaredDiff", "Tversky" , "CrossEntropy", "Jaccard"}, 

			MaxTrainingRounds -> rounds - ittTrain, BatchSize -> batch,
			TargetDevice -> "GPU", WorkingPrecision -> "Mixed",
			LearningRate -> 0.025, Method -> {"ADAM", "Beta1" -> 0.99, "Beta2"->0.999, "Epsilon"->10^-5},

			TrainingProgressFunction -> {monitorFunction, "Interval" -> Quantity[1, "Rounds"]},
			TrainingProgressReporting -> File[outName[StringReplace[DateString["ISODateTime"], 
				":" | "-" -> ""] <> ".json"]]
		];


		(*---------- Export the network ----------------*)

		netOut = NetExtract[trained["TrainedNet"], "net"];
		Export[outName["trained" <> ".wxf"], trained];
		Export[outName["final" <> ".wlnet"], netOut];
		Export[outName["final" <> ".onnx"], netOut];
	]
]


(* ::Subsubsection::Closed:: *)
(*MakeTestData*)


MakeTestData[data_, n_, patch_] := Block[{testData, len, sel, testDat},
	testData = data[[1]];
	len = Length@testData;
	If[len > First@patch,
		sel = Range @@ Clip[Round[(Clip[Round[len/3 - (0.5 n) First@patch], {0, Infinity}] + {1, n First@patch})], {1, len}, {1, len}];
		testData = First@AutoCropData[testData[[sel]]]
	];
	{{NormDatH@PadToDimensions[testData, patch]}, data[[3]]}
];


(* ::Subsubsection::Closed:: *)
(*ShowTrainLog*)


ShowTrainLog[fol_] := ShowTrainLog[fol, 5]

ShowTrainLog[fol_, max_] := Block[{files, log, keys, plots},
	{keys, log} = LoadLog[fol, max];

	(* Create a dynamic module to display the interactive plot *)
	DynamicModule[{pdat = log, klist = keys, plot, folder=fol},
		Manipulate[
			plot = Transpose[Values /@ Normal[pdat[All, key]][[All, All]]];
			ymaxv = Max[{1.1, 1.1 If[plot==={}, 1, Max[Flatten@plot]]}];
			ymax = Min[{ymax, ymaxv}];
			(* Plot the selected metrics *)
			ListLinePlot[If[key === {}, {}, plot], 
				PlotLegends -> Placed[key, Right], ImageSize -> 600, PlotRange->{{xmin,xmax} ,{0,ymax}}],
			(*the controls*)
			Row[{
				InputField[Dynamic[folder], String, Enabled -> False, FieldSize -> 50], 
				Button["Browse", 
					temp = SystemDialogInput["Directory", folder];
					If[StringQ[temp], folder = temp; {klist, pdat} = LoadLog[folder, max];xmax = Length[pdat];];
					, ImageSize -> {60, Automatic}, Method->"Queued"]}
			],
			Button["Reload", {klist, pdat} = LoadLog[folder, max];],
			Control[{{key, {}, ""}, klist, ControlType -> TogglerBar, 
				Appearance -> "Vertical" -> {Automatic, 4}, BaseStyle -> Medium}],
			Row[{
				Control[{{xmin, 1,"x min"},1, Dynamic[xmax-1], 1}], "   ", 
				Control[{{xmax,Length[pdat],"x max"}, Dynamic[xmin+1], Dynamic[Length[pdat]], 1}]}],
			Row[{
				Control[{{ymax, 1, "y max"}, 0.01, Dynamic[ymaxv]}],
				Button["Autoscale", ymax = Max[{1.1, 1.1 If[plot==={}, 1, Max[Flatten@plot]]}]]
			}],
			{{key, {}}, ControlType -> None},
			Initialization :> (
				key = {};
				plot = Transpose[Values /@ Normal[pdat[All, key]][[All, All]]];
				xmin = 1; xmax = Length@pdat; ymax = ymaxv = 1;
			)
		]
	]
]

LoadLog[fol_, max_]:=Block[{files, keys, log},
	(* Get a list of log files in the specified folder *)
	files = Sort[FileNames["*.json", fol]];
	
	(* Read the log files and extract the relevant information *)
	log = Flatten[Select[(Select[Import[#, "Lines"], StringContainsQ[#, "ProgressFraction"] &] & /@ files), Length[#] > max &], 1];
	
	(* Convert the log data into a dataset *)
	log = "[\n" <> StringDrop[StringRiffle[If[StringTake[#, -1] === "}", # <> ",", #] & /@ log, "\n"], -1] <> "\n]";
	log = Dataset[Association /@ Import[Export[FileNameJoin[{$TemporaryDirectory, "log.json"}], log, "text"]]];
	
	(* Get the unique keys (metrics) in the log data *)
	keys = Sort@DeleteDuplicates[Flatten[Normal@log[All, Keys]]];

	{keys, log}
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
		itt, i, datO, segO, dat, seg, vox, dim, augI, aug, nSet
	},

	itt = 0;
	datO = segO = {};

	{augI, nSet} = OptionValue[{AugmentData, PatchesPerSet}];

	itt = Ceiling[nBatch/nSet];

	Do[
		dat =RandomChoice[datas];
		
		If[StringQ[dat], 
			(*data is wxf file format*)
			{dat, seg, vox} = Import[dat];
			,
			If[Length[dat]===2, 
				(*datas is list of nii files {dat.nii, seg.nii}*)
				{seg, vox} = ImportNii[dat[[2]]];
				{dat, vox} = ImportNii[dat[[1]]];
				,
				(*data is in memory*)
				{dat, seg, vox} = dat;
			]
		];

		(*normalize the data and segmentations and pad to at least patch dimensions*)
		dim = Max /@ Transpose[{Dimensions@dat, patch}];
		dat = NormDatH@PadToDimensions[dat, dim];
		seg = PadToDimensions[seg, dim];

		(*check if augmentation is a boolean or a list*)
		aug = Which[
			augI === 1, True,
			augI === 0, False,
			0 < augI < 1, RandomChoice[{augI, 1 - augI} -> {True, False}],
			BooleanQ[augI], augI,
			True, True];
		(* {flip, rot, trans, scale, noise, blur, bright} *)
		aug = (# && aug) & /@ {True, True, True, True, True, False, False};

		(*perform augmentation on full data and get the defined number of patches*)
		{dat, seg} = AugmentTrainingData[{dat, seg}, vox, aug];
		{dat, seg} = PatchTrainingData[{dat, seg}, patch, nSet];

		datO = Join[datO, dat];
		segO = Join[segO, seg];
	, itt];

	datO = datO[[;; nBatch]];
	segO = If[IntegerQ[nClass], ClassEncoder[segO[[;; nBatch]] + 1, nClass], segO[[;; nBatch]] + 1];
	Thread[Transpose[{datO}] -> segO]
];


(* ::Subsubsection::Closed:: *)
(*AugmentImageData*)


AugmentImageData[im_?ListQ, {rot_, flip_}]:=AugmentImageData[#, {rot, flip}]&/@im

AugmentImageData[im_, {rot_, flip_}]:=Block[{rt, fl, tr},
	rt = If[rot, RotationTransform[RandomReal[{-90, 90}]Degree], TranslationTransform[{0, 0}]];
	fl = If[flip&&RandomChoice[{True, False}], ReflectionTransform[{1, 0}], TranslationTransform[{0, 0}]];
	tr = rt . fl;
	If[Head[im]===Rule,
		ImageTransformation[im[[1]], tr, DataRange->{{-.5, .5}, {-.5, .5}}]->im[[2]],
		ImageTransformation[im, tr, DataRange->{{-.5, .5}, {-.5, .5}}]
	]]


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
		w = Join[
			If[rot, {1., 1., 1.} {RandomReal[{-180, 180}], RandomReal[{-10, 10}], RandomReal[{-10, 10}]}, {0., 0., 0.}],
			If[trans, {1., 1., 1.} (Dimensions[dat] vox)  RandomReal[{-0.05, 0.05}, 3], {0., 0., 0.}],
			If[scale, {1., 1., 1.} RandomReal[{0.5, 2}, 3], {1., 1., 1.}],
			{0., 0., 0.}
		];
		
		{datT, segT} = DataTransformation[#, vox, w, InterpolationOrder -> 0, PadOutputDimensions -> False]&/@{datT, segT};
	];
	
	(*Augmentations of sharpness intensity and noise*)
	If[bright, datT = RandomChoice[{RandomReal[{1, 1.5}], 1/RandomReal[{1, 1.5}]}] datT];
	If[blur, datT = GaussianFilter[datT, RandomReal[{0.1, 1.5}]]];
	If[(noise && RandomChoice[{0.3, 0.7} -> {True, False}]), datT = addNoise[datT, Mean[Flatten[datT]]/RandomReal[{10, 150}], RandomChoice[{0.8, 0.2} -> {0, 1}] RandomReal[{0, 0.5}]/100]];
	
	{ToPackedArray[N[datT]], ToPackedArray[Round[segT]]}
]


ReverseC = Compile[{{dat, _Real, 1}}, Reverse[dat], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


addNoise[data_, sigma_, p_] := Block[{dims, g1, g2, num, coors},
	dims = Dimensions[data];
	(*random noise*)
	{g1, g2} = RandomReal[NormalDistribution[0., sigma], Prepend[dims, 2]];
	(*salt and pepper noise*)
	num = Max[{1, Round[p*Times @@ dims/2]}];
	coors = Transpose[RandomInteger[{1, #}, 2 num] & /@ dims];
	(*add the noise*)
	SaltAndRiceC[data, coors, num, g1, g2]
];


SaltAndRiceC = Compile[{{data, _Real, 3}, {coors, _Integer, 2}, {num, _Integer, 0}, {g1, _Real, 3}, {g2, _Real, 3}}, Block[{newData},
	newData = Sqrt[(data + g1)^2 + g2^2];
	(newData[[#[[1]], #[[2]], #[[3]]]] = 1) & /@ coors[[1 ;; num]];
	(newData[[#[[1]], #[[2]], #[[3]]]] = 0) & /@ coors[[num + 1 ;; -1]];
	newData
], RuntimeAttributes -> {Listable}, Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*PatchTrainingData*)


PatchTrainingData[{dat_,seg_}, patch_, n_]:=Block[{pts,datP,segP},
	{datP, pts} = DataToPatches[dat, patch, n, PatchNumber->2];
	segP = DataToPatches[seg, patch, pts][[1]];
	{ToPackedArray[N[#]]&/@datP, ToPackedArray[Round[#]]&/@segP}
]



(* ::Subsubsection::Closed:: *)
(*PrepTrainData*)


PrepTrainData[dat_, seg_]:= PrepTrainData[dat, seg, {0}]

PrepTrainData[dat_, seg_, labi_?VectorQ]:= PrepTrainData[dat, seg, {labi, labi}]

PrepTrainData[dat_, seg_, {labi_?VectorQ, labo_?VectorQ}] := Block[{cr},
	cr = FindCrop[dat  Mask[NormalizeData[dat], 5, MaskDilation -> 1]];
	{
		NormDat[ApplyCrop[dat, cr]], 
		If[labi==={0},
			ApplyCrop[seg,cr],
			ReplaceSegmentations[ApplyCrop[seg, cr], labi, labo]
		]
	}
]


(* ::Subsubsection::Closed:: *)
(*SplitDataForSegmentation*)


SyntaxInformation[SplitDataForSegementation] = {"ArgumentsPattern" -> {_, _.}};

SplitDataForSegementation[data_, seg_]:=Block[{dat,pts,dim,loc,set, segp},
	{{dat, pts, dim}, loc, set} = SplitDataForSegementation[data];
	segp = GetPatch[seg, pts];
	{{dat, pts, dim}, {segp, pts, dim}, loc,set}
]

SplitDataForSegementation[data_]:=Block[{dim,whatSide,side,whatPos,pos,dat,right,left,cut,pts,loc},
	dim = Dimensions[data];

	(*find which side using NN*)
	whatSide = ClassifyData[data, "LegSide"];

	(*based on side cut data or propagate*)
	dat=If[whatSide==="Both",
		{right, left, cut} = CutData[data];
		{{right, {"Right", {1, cut}}}, {left, {"Left", {cut+1, dim[[3]]}}}},
		{{data, cut=0; {whatSide, {1, dim[[3]]}}}}
	];

	(*loop over data to find upper or lower*)
	dat = Flatten[(
		{dat, side} = #;
		{whatPos, pos} = ClassifyData[dat, "LegPosition"];

		Switch[whatPos,
			(*if upper and lower split upper and lower*)
			"Both",{{dat[[pos[[1]];;]],{"Upper",{pos[[1]],dim[[1]]}},side},{dat[[;;pos[[2]]]],{"Lower",{1,pos[[2]]}},side}},
			(*if only knee data duplicate for both networks*)
			"Knee",{{dat,{"Upper",{1,dim[[1]]}},side},{dat,{"Lower",{1,dim[[1]]}},side}},
			(*if only upper or only lower return what it is*)
			_,{{dat, {whatPos, {1, dim[[1]]}}, side}}
		]
	)&/@dat, 1];

	{dat, pts, loc} = Transpose[CropPart/@dat];

	{{dat,pts,dim}, loc, {{whatSide, cut}, {whatPos, pos}}}
]


(* ::Subsubsection::Closed:: *)
(*CropPart*)


CropPart[data_]:=Block[{dat,up,sid,upst,upend,sidst,sidend,crp},
	{dat, {up, {upst, upend}}, {sid, {sidst, sidend}}} = data;

	{dat, crp} = AutoCropData[Dilation[Normal[TakeLargestComponent[Mask[NormalizeData[dat],10]]],1] dat, CropPadding->0];
	{dat, Partition[crp,2]+{upst-1,0,sidst-1}, {up,sid}}
]


(* ::Subsection:: *)
(*Make evaluation images*)


(* ::Subsubsection::Closed:: *)
(*MakeChannelClassImage*)


SyntaxInformation[MakeChannelClassImage]={"ArgumentsPattern"->{_, _, _., _.}};

MakeChannelClassImage[data_, label_]:=MakeChannelClassImage[data,label,MinMax[label], {1,1,1}]

MakeChannelClassImage[data_,label_, {off_, max_}]:=MakeChannelClassImage[data, label,{off,max}, {1,1,1}]

MakeChannelClassImage[data_, label_, vox_]:=MakeChannelClassImage[data, label, MinMax[label], vox]

MakeChannelClassImage[data_, label_, {off_, max_}, vox_]:=Block[{i1, i2},
	i1 = MakeClassImage[label, {off, max}, vox];
	i2 = MakeChannelImage[data, vox];
	ImageCollage[ImageCompose[#, SetAlphaChannel[i1, 0.4 AlphaChannel[i1]]]& /@ i2]
]


(* ::Subsubsection::Closed:: *)
(*MakeClassImage*)


SyntaxInformation[MakeClassImage]={"ArgumentsPattern"->{_, _., _.}};

MakeClassImage[label_]:=MakeClassImage[label, MinMax[label], {1,1,1}]
 
MakeClassImage[label_, {off_, max_}]:=MakeClassImage[label, {off, max}, {1,1,1}]

MakeClassImage[label_, vox_]:=MakeClassImage[label, MinMax[label], vox]

MakeClassImage[label_,{off_, maxI_}, vox_]:=Block[{max, cols, imlab, rat},
	(*SeedRandom[1345];
		cols = Prepend[ColorData["DarkRainbow"][#]&/@RandomSample[Rescale[Range[off+1, max]]],Transparent];
		cols = Prepend[ColorData["RomaO"][#]&/@Rescale[Range[off+1, max]],Transparent];
	*)
	max = Max[{Max[label], maxI}];
 	cols = Prepend[ColorData["RomaO"][#]&/@Rescale[Join[Select[Range[off + 1, max], EvenQ], Select[Range[off + 1, max], OddQ]]],Transparent];
	imlab = Round@Clip[If[ArrayDepth[label] === 3, label[[Round[Length@label/2]]], label] - off + 1, {1, max + 1}, {1, 1}];
	rat = vox[[{2,3}]]/Min[vox[[{2,3}]]];
	ImageResize[Image[cols[[#]]&/@imlab], Round@Reverse[rat Dimensions[imlab]], Resampling->"Nearest"]
]


(* ::Subsubsection::Closed:: *)
(*MakeChannelImage*)


SyntaxInformation[MakeChannelImage]={"ArgumentsPattern"->{_, _., _.}};

MakeChannelImage[data_]:=MakeChannelImage[data, {1, 1, 1}]

MakeChannelImage[data_, vox_]:=Block[{dat, imdat, rat},
	dat = Rescale[data];
	rat = vox[[{2, 3}]] / Min[vox[[{2, 3}]]];
	(
		imdat = #;
		imdat = If[ArrayDepth[#]===3, imdat[[Round[Length@imdat/2]]], imdat];
		ImageResize[Image[Clip[imdat, {0,1}]], Round@Reverse[rat Dimensions[imdat]], Resampling->"Nearest"]
	) &/@ dat
]


(* ::Subsection:: *)
(*Muscle Names*)


(* ::Subsubsection::Closed:: *)
(*ImportITKLabels*)


ImportITKLabels[file_] := Block[{lines, muscleNames, muscleLabels},
	(*import*)
	lines = Select[Import[file, "Lines"], StringTake[#, 1] =!= "#" &];
	(*extract names and numbers*)
	muscleNames = StringRiffle[Capitalize[ToLowerCase[Select[#, ! IntegerQ[ToExpression[#]] &]]], "_"] & /@ 
		StringSplit[(Select[StringTrim[#], (# =!= "\t" && # =!= "") &] & /@ 
		StringSplit[lines, "\""])[[All, -1]]];
	muscleLabels = ToExpression[StringSplit[#, " "][[1]]] & /@ lines;
	
	(*output*)
	{muscleNames, muscleLabels}
]


(* ::Subsubsection::Closed:: *)
(*MuscleLabelToName*)


MuscleLabelToName[num_] := num /. Thread[#[[2]] -> #[[1]]] &[ImportITKLabels[GetAssetLocation["MusclesLegLabels"]]]

MuscleLabelToName[num_, file_] := Block[{muscleNames, muscleLabels},
	{muscleNames, muscleLabels} = ImportITKLabels[file];
	num /. Thread[muscleLabels -> muscleNames]
]


(* ::Subsubsection::Closed:: *)
(*MuscleLabelToName*)


MuscleNameToLabel[num_] := num /. Thread[#[[1]] -> #[[2]]] &[ImportITKLabels[GetAssetLocation["MusclesLegLabels"]]]

MuscleNameToLabel[num_, file_] := Block[{muscleNames, muscleLabels},
	{muscleNames, muscleLabels} = ImportITKLabels[file];
	num /. Thread[muscleNames -> muscleLabels]
]


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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
