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


(* ::Subsection::Closed:: *)
(*Functions*)


GetNeuralNet::usage = 
"GetNeuralNet[name] loads a pre trained neural net that come with the toolbox. Current named nets 
are \"LegSide\", \"LegSide\", \"SegThighMuscle\", \"SegLegMuscle\", and \"SegLegBones\". The loading is cashed within a session."


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
"AddLossLayer[net] adds three loss layers to a NetGraph, a DiceLossLayer, JaccardLossLayer, TverskyLossLayer, MeanSquaredLossLayer and a CrossEntropyLossLayer are added."

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


ClassEncoder::usage = 
"ClassEncoder[label] encodes Integer label data of 0 to max value of label into a nClass + 1 vector of 1 and 0 as the last dimension.
ClassEncoder[label, nClass] encodes Integer label data of 0 to nClass into a nClass + 1 vector of 1 and 0 as the last dimension."

ClassDecoder::usage = 
"ClassDecoder[probability] decodes a probability vector of 1 and 0 into Integers of 0 to the value of the last dimension of probability minus one.
ClassDecoder[probability, nClass] decodes a probability vector of 1 and 0 into Integers of 0 to nClass - 1."


DiceSimilarity::usage = 
"DiceSimilarity[ref, pred] gives the Dice Similarity between 1 and 0 of segmentations ref and pred for class equals 1.
DiceSimilarity[x, y, class] gives the Dice Similarity of segmentations ref and pred for class.
DiceSimilarity[x, y, {class, ..}] gives the Dice Similarity of segmentations ref and pred for the list of gives classes."

JaccardSimilarity::usage = 
"JaccardSimilarity[ref, pred] gives the Jaccard Similarity between 1 and 0 of segmentations ref and pred for class equals 1.
JaccardSimilarity[x, y, class] gives the Jaccard Similarity of segmentations ref and pred for class.
JaccardSimilarity[x, y, {class, ..}] gives the Jaccard Similarity of segmentations ref and pred for the list of gives classes."

SurfaceDistance::usage = 
"SurfaceDistance[ref, pred] gives the mean surface distance of segmentations ref and pred for class equals 1 in voxels.
SurfaceDistance[x, y, class] gives the mean surface distance of segmentations ref and pred for class in voxels.
SurfaceDistance[x, y, {class, ..}] gives the mean surface distance of segmentations ref and pred for the list of gives classes in voxels.
SurfaceDistance[x, y, class , vox] gives the mean surface distance of segmentations ref and pred for class in millimeter.
SurfaceDistance[x, y, {class, ..}, vox] gives the mean surface distance of segmentations ref and pred for the list of gives classes in millimeters."

MakeDistanceMap::usage = 
"MakeDistanceMap[mask] makes a distance map of the given mask in voxels. The distance map is negative inside the mask and positive outside the mask.
MakeDistanceMap[mask, vox] makes a distance map of the given mask in the same unit as vox. The distance map is negative inside the mask and positive outside the mask."


SegmentData::usage = 
"SegmentData[data, what] segments the data. The what specifies the segmentation to be done.
It currently allows for \"LegBones\" for the bones or \"Legs\" for the muscles."

ApplySegmentationNetwork::usage = 
"ApplySegmentationNetwork[data, net] segments the data using the pre trained net.
ApplySegmentationNetwork[data, net, node] segments the data using the pre trained net but only use the network upto node."

ClassifyData::usage = 
"ClassifyData[data, method] classifies the input data using the given method. The data is converted to images using MakeClassifyImages.
The input method can be a filename of a classify network or a classify network. 
Additionally the input method can be one of the predefined methods \"LegPosition\" or \"LegSide\"."


ShowTrainLog::usage =
"ShowTrainLog[log] shows the training log of a network training."


TrainSegmentationNetwork::usage =
"TrainSegmentationNetwork[{inFol, outFol}] trains a segmentation network. The correctly prepared training data should be stored in inFol. The progress each round will be saved in outFol.
TrainSegmentationNetwork[{inFol, outFol}, netCont] does the same but defines how to continue with netCont. If netCont is \"Start\" training will be restarted.
If netCont is a initialized network or network file (wlnet) this will be used. If netCont is a a outFol the last saved network will be used.
Possible loss functions are {\"SoftDice\", \"SquaredDiff\", \"Tversky\" , \"CrossEntropy\", \"Jaccard\"}."

GetTrainData::usage =
"GetTrainData[data, batch size, patch] creates a training batch of size batch size with patch size patch. 
The input data can be out of memory in the form of a list of \"*wxf\" files that contain the data, segmentation and voxel size or a list of \"*.nii\" files in the form
{{\"data.nii\", \"segmentation.nii\"}..}. The input data can be in memory in a list in the form {{data, segmentation, vox}..}
GetTrainData[data, batch size, patch, nClass] If nClass is set to an value n > 0 the segmentations are decoded in n classes."

PrepareTrainingData::usage = 
"PrepareTrainingData[inFolder, outFolder] prepares the data in de inFolder for training a neural network for segmentation and outputs in outFolder.
PrepareTrainingData[{labFolder, datFolder}, outFolder] does the same but the labels are stored in labFolder and data is stored in datFolder."


CheckSegmentation::usage=
"CheckSegmentation[seg] checks the segmentation for errors and returns a vector of two numbers, the first indicates if the segmentation has more than one region, the second indicates if it has holes."

DataToPatches::usage =
"DataToPatches[data, patchSize] creates the maximal number of patches with patchSize from data, where the patches have minimal overlap.
DataToPatches[data, patchSize, n] gives n random patches from the maximal number of patches with patchSize from data, where the patches have minimal overlap."

PatchesToData::usage = 
"PatchesToData[patches, ran] creates a continuous dataset from the patches. For each patch the range in the data needs to be specified in ran.
The patches are have dimensions {x, y, z} each and ran is specified as {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}}.
PatchesToData[patches, ran, dim] creates a continuous dataset from the patches with dimensions dim."


AugmentTrainingData::usage = 
"AugmentTrainingData[{data, segmentation}, vox] augments the data and segmentation in the same way.
AugmentTrainingData[{data, segmentation}, vox, aug] by setting aug to True or False the augmentation can be turned on or off.
The value aug can also be a list of boolean values controlling various augmentation parameters {flip, rotate, translate, scale, noise, blur, brightness}.
The default settings are {True, True, True, True, False, False, False}."

AugmentImageData::usage = 
"AugmentImageData[image, {rotate, flip}] augments the input image by rotating between -180 and 180 degrees and flipping. The inputs rotate and flip
can be set to True or False.
AugmentImageData[{image, ..}, {rotate, flip}] same but for a list of images."


MakeChannelImage::usage = 
"MakeChannelImage[data] makes a cross-sectional image of the channels data of a training dataset generated by GetTrainData.
MakeChannelImage[data, vox] same but with the aspect ratio determined by vox."

MakeClassImage::usage = 
"MakeClassImage[label ] makes a cross-sectional image of the classes label of a training dataset generated by GetTrainData
MakeChannelImage[label, {b, n}] same but with explicit definition of background value b and number of classes n. 
MakeClassImage[data, vox] same but with the aspect ratio determined by vox.
MakeChannelImage[label, {b, n}, vox] same with explicit definition and aspect ratio definition."

MakeChannelClassImage::usage = 
"MakeChannelClassImage[data, label] makes a cross-sectional image of the channels data overlaid with a cross-sectional image of the classes label of a training dataset generated
MakeChannelClassImage[data, label, {off,max}] same but with explicit definition of background value b and number of classes n. 
MakeChannelClassImage[data, label, vox] same but with the aspect ratio determined by vox.
MakeChannelClassImage[data, label, {off,max}, vox] same with explicit definition and aspect ratio definition."

MakeChannelClassGrid::usage =
"MakeChannelClassGrid[data, label] makes a 3 x 3 grid of cross-sectional images of the channels data overlaid with a cross-sectional image of the classes label of a training dataset generated
MakeChannelClassGrid[data, label, n] makes a n x n.
MakeChannelClassGrid[data, label, {n, m}] makes a n x m."


SplitDataForSegmentation::usage = 
"SplitDataForSegmentation[data] is a specific function for leg data to prepare data for segmentation. It detects the side and location and will split and label the data accordingly.
SplitDataForSegmentation[data ,seg] does the same but is rather used when preparing training data. Here the seg is split in exactly the same way as the data."


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

LoadTrainingData::usage =
"LoadTrainingData is an option for TrainSegmentationNetwork. If set to True the training data is loaded from the disk."

MonitorInterval::usage =
"MonitorInterval is an option for TrainSegmentationNetwork. It defines how often the training is monitored."

L2Regularization::usage =
"L2Regularization is an option for TrainSegmentationNetwork. It defines the L2 regularization factor."

MultiChannel::usage = 
"MultiChannel ..."


PatchesPerSet::usage =
"PatchesPerSet is an option for GetTrainData. Defines how many random patches per dataset are created within the batch."

AugmentData::usage = 
"AugmentData is an option for GetTrainData and TrainSegmentationNetwork. If set True the training data is augmented."

PadData::usage =
"PadData is an option for GetTrainData and TrainSegmentationNetwork.. If set to an integers the that number of slices on the top and bottom of the 
data can be made 0. This is done to learn cut of datasets."

PatchSize::usage =
"PatchSize is an option for TrainSegmentationNetwork. Defines the patch size used in the network training."

RoundLength::usage = 
"RoundLength is an option for TrainSegmentationNetwork. Defines how many batches will be seen during each training round."


MaxPatchSize::usage = 
"MaxPatchSize is an option for SegmentData and ApplySegmentationNetwork. Defines the patch size used when segmenting data. Bigger patches are better."

DataPadding::usage = 
"DataPadding is an option for ApplySegmentationNetwork. Defines how much to pad the data patches in all directions."


PatchNumber::usage = 
"PatchNumber is an option for DataToPatches. Can be an integer value >= 0. The larger the number the more overlap the patches have.
The minimal number of patches in each direction is calculated, and then for each dimension the given number is added."

PatchPadding::usage = 
"PatchPadding is an option for DataToPatches. Can be an integer value >= 0. It pads the chosen patch size with the given number."

LabelTag::usage = 
"LabelTag is an option for PrepareTrainingData. It defines the tag used in the filenames of the label data."

DataTag::usage = 
"DataTag is an option for PrepareTrainingData. It defines the tag used in the filenames of the data."

InputLabels::usage = 
"InputLabels is an option for PrepareTrainingData. Can be set to a list of integers corresponding to the labels to be used from the given segmentation."

OutputLabels::usage = 
"OutputLabels is an option for PrepareTrainingData. Can be set to a list of integers. The used label number will be replaced by these numbers."

TestRun::usage = 
"TestRun is an option for PrepareTrainingData. If set to True the data is not saved only analyzed."

CleanUpSegmentations::usage = 
"CleanUpSegmentations is an option for PrepareTrainingData. If set to True the segmentations are cleaned up by removing holes reducing to one volume and smoothing."

TrainVoxelSize::usage =
"TrainVoxelSize is an option for PrepareTrainingData. It defines the voxel size of the training data. When set to Automatic the voxel size is that of the data."

DistanceRange::usage =
"DistanceRange is an option for MakeDistanceMap. It defines the range of the distance map outside the segmentation in voxels.
Values can be Automatic, All, or a integer value. If All the distance map is calculated for the whole image. If 0 the distance map is only calculated inside the segmentation."



(* ::Subsection::Closed:: *)
(*Error Messages*)


TrainSegmentationNetwork::net = "The net input is not \"Start\", a network, a network file, or a previous train folder."

TrainSegmentationNetwork::cont = "Could not find a previous network in the specified folder."

TrainSegmentationNetwork::inp = "The string input given is not a network file or a directory."

TrainSegmentationNetwork::itt = "Not enough iterations specified for training. Remaining iterations are less than 5."


GetTrainData::aug = "The augmentation input is not a number or a boolean value. Using False by default."


MakeUnet::scale = "The scaling input is not valid. It should be Automatic. It can also be a integer or a list of integers that will be applied to the Layers. 
It can also be a vector of integers per layer where the length of the vector should be equal to the depth of the network.";

MakeUnet::sett = "The setting input is not valid. It can be a number or a list of numbers that will be applied to the Layers.";

MakeUnet::drop = "The dropout input is not valid. It can be a number or a list of numbers that will be applied to the Layers.";

MakeUnet::feat = "The feature input is not valid. It can be a number or a list of numbers that will be applied to the Layers.";

MakeUnet::arch = "The architecture input is not valid. It can be \"UNet\", \"UNet+\", or \"UNet++\".";

MakeUnet::block = "The block type input is not valid. It can be \"Conv\", \"UNet\", \"ResNet\", \"DenseNet\", \"Inception\", or \"U2Net\".";


SurfaceDistance::met = "Method `1` not recognized";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


debugUnet[x___] := If[$debugUnet, Print[x]];


(* ::Subsection::Closed:: *)
(*GetNeuralNet*)


SyntaxInformation[GetNeuralNet] = {"ArgumentsPattern" -> {_}};

(*allows for caching and clearing while GetNeuralNet remains protected*)
GetNeuralNet["Clear"] := (
	ClearAll[GetNeuralNetI];
	GetNeuralNetI[name_] := GetNeuralNetI[name] = NeuralNetFunc[name]
)

GetNeuralNet[name_?StringQ] := GetNeuralNetI[name]

GetNeuralNetI[name_] := GetNeuralNetI[name] = NeuralNetFunc[name]

NeuralNetFunc[name_] := GetNeuralNetI[name] = Which[
	FileExistsQ[name], Import[name],
	FileExistsQ[GetAssetLocation[name]], Import[GetAssetLocation[name]],
	True, $Failed
]


(* ::Subsection:: *)
(*MakeUnet*)


(* ::Subsubsection::Closed:: *)
(*Settings*)


netDefaults={"Conv"->1,"UNet"->2,"ResNet"->2,"DenseNet"->{4,2},"Inception"->{4,2},"U2Net"->{3,True},"Map"->1};


(* ::Subsubsection::Closed:: *)
(*MakeUnet*)


Options[MakeUnet] = {
	NetworkArchitecture -> "UNet",
	BlockType -> "ResNet", 
	ActivationType -> "GELU",
	RescaleMethod -> "Conv",
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
		ndim, depth, nam, drop, blockType, actType, scaling, feature, setting, mon, net, boolSc, metSc,
		architecture, cons, depthj, fStart, mapCon
	},

	{drop, blockType, actType, depth, scaling, feature, setting, mon, metSc, architecture} = 
		OptionValue[{DropoutRate, BlockType, ActivationType, NetworkDepth, 
			DownsampleSchedule, FeatureSchedule, SettingSchedule, MonitorCalc, 
			RescaleMethod, NetworkArchitecture}];

	(*check the input*)
	{architecture, mapCon} = If[StringQ[architecture], {architecture, Automatic}, 
		If[Length[architecture]===2, architecture, Return[Message[MakeUnet::arch]; $Failed]]];
	If[!MemberQ[{"UNet", "UNet+", "UNet++"}, architecture], 
		Return[Message[MakeUnet::arch]; $Failed]];
	If[!MemberQ[{"Conv", "UNet", "ResNet", "DenseNet", "Inception", "U2Net"}, blockType], 
		Return[Message[MakeUnet::block]; $Failed]];

	(*is the network 2D or 3D*)
	ndim = Length@dimIn;

	If[mon, 
		Echo[{architecture, blockType, actType}, "Network block type: "];
		Echo[{dimIn, ndim}, "Network dimension order: "];
	];

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
	If[mon, Echo[scaling, "Network scaling schedule: "]];

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
	If[mon, Echo[setting, "Network setting schedule: "]];

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
		],
		True, Return[Message[MakeUnet::feat]; $Failed]
	]];
	If[VectorQ[feature], feature = Round[Transpose[{feature, feature/2}]]];
	fStart = Last@Flatten@{First@feature};
	If[mon, Echo[feature, "Network feature schedule: "]];

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
					{1, If[i==depth, 1, scaling[[i]]]}, 
					(*skip in -> only input from above, skip out -> always for enc*)
					{0, True}, 
					(*config*)
					{{blockType, setting[[i]]}, feature[[i]], {actType, ndim}},	DropoutRate -> drop[[i]], RescaleMethod -> metSc
				], {i, 1, depth}],
				(*decoding layers*)
				Table[debugUnet[nam["dec_", i]]; nam["dec_", i] -> MakeNode[
					(*scale up -> always, scale down -> never*)
					{scaling[[i]], 1}, 
					(*skip in -> accepts one skip, skip out -> never for dec*)
					{1, False}, 
					(*config*)
					{{blockType, setting[[i]]}, feature[[i]], {actType, ndim}}, DropoutRate -> drop[[i]], RescaleMethod -> metSc
				], {i, 1, depth - 1}],
				(*start and mapping*)
				{"start" -> UNetStart[nChan, fStart, dimIn, actType],
				"map" -> UNetMap[ndim, nClass]}
			],

			(*define the connections*)
			cons = Join[
				Flatten@Table[{
					(*connect encoding layers*)
					NetPort[nam["enc_", i], If[boolSc[[i]], "Down", "Skip"]] -> nam["enc_", i + 1],
					(*connect encoding to decoding layers*)
					NetPort[nam["enc_", i], "Skip"] -> NetPort[nam["dec_", i], "Skip"],
					(*connect decoding layers*)
					If[i === depth - 1, NetPort[nam["enc_", i + 1], "Skip"], nam["dec_", i + 1]] -> NetPort[nam["dec_", i], "Scale"]
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
				Flatten@Table[debugUnet[nam[i, j]];nam[i, j] -> MakeNode[
					(*upscale for all nodes accept backbone -1, downscale only for backbone*)
					{If[j > 1, scaling[[i]], 1], If[j =!= 1 || i == depth, 1, scaling[[i]]]}, 
					(*skip in for all accept backbone, for UNET++ name the skips, skip out for all except right most upscale*)
					{If[j > 1, If[architecture==="UNet++", j - 1, 1], 0], If[j < depthj - i, True, False]}, 
					(*config*)
					{{blockType, setting[[i]]}, feature[[i]], {actType, ndim}},	DropoutRate -> drop, RescaleMethod -> metSc
				], {i, 1, 5}, {j, 1, depthj - i}],
				(*start and mapping*)
				{"start" -> UNetStart[nChan, fStart, dimIn, actType],
				"map" -> UNetMap[ndim, nClass, mapCon]}
			],
			cons = Join[
				Flatten@Table[{
					(*connect the backbone, the downscaling*)
					If[j === 1 && i=!= depth, NetPort[nam[i, j], If[boolSc[[i]], "Down", "Skip"]] -> nam[i + 1, j], Nothing],
					(*connect the nodes with up scaling*)
					If[1 < i <= depth, If[i=== depth, nam[i, j], NetPort[nam[i, j], If[j==depthj-i, "Up", "Skip"]]] -> NetPort[nam[i-1, j+1], "Scale"],  Nothing],
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
			debugUnet["The node connection list: ", cons]; cons, 

			(*network settings and options*)
			"Input" -> Join[{nChan}, dimIn]		
		]
	];
	(*Monitor network properties*)
	If[mon, Echo[NetSummary[net], "Network description: "]];

	(*return network*)
	net
]


(* ::Subsubsection::Closed:: *)
(*UNetStart*)


UNetStart[nChan_, feat_, dimIn_, actType_] := NetGraph[
	NetChain[Conv[feat , {Length@dimIn, 1}, actType], 
		"Input" -> Prepend[dimIn, nChan]]
]


(* ::Subsubsection::Closed:: *)
(*ClassMap*)


UNetMap[dim_, nClass_]:=UNetMap[dim, nClass, 1]

UNetMap[dim_, nClass_, n_] := Block[{map},
	map = Flatten[{ConvolutionLayer[nClass, 1], 
		If[nClass > 1, 
			{TransposeLayer[Which[dim === 2, {3, 1, 2}, dim === 3, {4, 1, 2, 3}]], SoftmaxLayer[]}, 
			{LogisticSigmoid, FlattenLayer[1]}
		]}];
	map = If[n>1,Prepend[map, CatenateLayer[]], map];
	NetGraph@NetChain@map
]


(* ::Subsubsection::Closed:: *)
(*MakeNode*)


Options[MakeNode] = {
	DropoutRate -> 0,
	RescaleMethod -> "Pool"
}

SyntaxInformation[MakeNode] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

MakeNode[{scUp_ ,scDown_}, {skIn_, skOut_}, blockConfig_, OptionsPattern[]] := Block[{
		drop, metSc, chan, dim, block, scaleUp, boolScUp, scale, scaleDown, boolScDown, 
		boolDr, con, boolCat, cat, catIn, skip, node
	},
	(*get the node parameters*)
	{drop, metSc} = OptionValue[{DropoutRate, RescaleMethod}];
	chan = First@Flatten@{blockConfig[[2]]};
	dim = blockConfig[[3, 2]];

	debugUnet[Column@{{scUp ,scDown}, {skIn, skOut}, {drop, metSc}, blockConfig}];

	(*define the block*)
	block = "block" -> ConvBlock @@ blockConfig;

	(*define the up scaling*)
	scaleUp = scUp /. False -> 1;
	scaleUp = If[IntegerQ[scaleUp], ConstantArray[scaleUp, dim], scaleUp];
	boolScUp = Times @@ scaleUp =!= 1;
	scaleUp = If[boolScUp, "scaleU" -> ConvScale[{"Decode", metSc}, scaleUp, chan], Nothing];
	scale = If[boolScUp, NetPort["Scale"] -> "scaleU", NetPort["Scale"]];

	(*define the downscaling*)
	scaleDown = scDown /. False -> 1;
	scaleDown = If[IntegerQ[scaleDown], ConstantArray[scaleDown, dim], scaleDown];
	boolScDown = Times @@ scaleDown =!= 1;
	scaleDown = If[boolScDown, "scaleD" -> ConvScale[{"Encode", metSc}, scaleDown, chan], Nothing];

	(*define the dropout*)
	boolDr = 0 < drop;
	drop = If[boolDr, "drop" -> DropoutLayer[drop], Nothing];
	con = If[boolDr, "drop", "block"];

	(*define the skip connections*)
	boolCat = skIn > 0 (*&& boolScUp*);
	cat = If[boolCat, "cat" -> CatenateLayer[], Nothing];
	catIn = If[boolCat, Append[If[skIn>1, Table[NetPort["Skip"<>ToString[i]], {i, skIn}], {NetPort["Skip"]}], scale], Nothing];
	skip = If[skOut=!=False, con -> NetPort["Skip"<>If[IntegerQ[skOut], ToString[skOut], ""]], Nothing];

	(*return the node*)
	node = NetGraph@NetFlatten[NetGraph[{block, drop, scaleUp, scaleDown, cat
		(*cat, block, scaleUp, scaleDown, drop*)}, {
			If[boolDr, "block" -> "drop", Nothing],
			If[boolScUp && skOut===False, con -> NetPort["Up"], Nothing],
			If[boolScDown, con -> "scaleD" -> NetPort["Down"], Nothing],
			If[boolCat, catIn -> "cat" -> "block", Nothing],
			skip
	}],1];

	debugUnet[node];
	node
]


(* ::Subsubsection::Closed:: *)
(*ConvScale*)


ConvScale[type_, scaleVec_, chan_] := Switch[type,
	{"Encode","Pool"}, 
	{PoolingLayer[scaleVec, scaleVec]},
	{"Encode","Conv"}, 
	{ConvolutionLayer[chan, scaleVec, "Stride" -> scaleVec]},
	{"Decode","Pool"}, 
	{ResizeLayer[Scaled /@ scaleVec, Resampling -> "Nearest"]},
	{"Decode","Conv"}, 
	{ResizeLayer[Scaled /@ scaleVec, Resampling -> "Nearest"], ConvolutionLayer[chan, scaleVec, 
		"Stride" -> 0 scaleVec + 1, 
		"PaddingSize" -> (If[OddQ[#], {(Ceiling[#/2] - 1), (Ceiling[#/2] - 1)}, {0, #/2}] & /@ scaleVec)]}
]


(* ::Subsubsection::Closed:: *)
(*ConvBlock*)


ConvBlock[block_, feat_?IntegerQ, dim_?IntegerQ] := ConvBlock[block, {feat, Round[feat/2]}, {"None", dim, 1}]

ConvBlock[block_, feat_?IntegerQ, {act_, dim_}] := ConvBlock[block, {feat, Round[feat/2]}, {act, dim, 1}]

ConvBlock[block_, feat_?IntegerQ, {act_, dim_, dil_}] := ConvBlock[block, {feat, Round[feat/2]}, {act, dim, dil}]

ConvBlock[block_, {featOut_, featInt_}, {act_, dim_}] := ConvBlock[block, {featOut, featInt}, {act, dim, 1}]

ConvBlock[block_, {featOut_, featInt_}, {act_, dim_, dil_}] := Block[{
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
		ConvBlock[{"UNet", blockSet}, {featOut, featInt}, {act, dim, dil}],

		"UNet",
		dep = blockSet;
		NetGraph@NetFlatten@NetGraph[{
			"conv" -> Flatten[Table[Conv[featOut, {dim, 3, dil}, act], {i, 1, dep}]]}, {}
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
		], 2],

		"Inception",
		{dep, rep} = If[IntegerQ[blockSet], {blockSet, 2}, blockSet];
		repBlock = NetGraph[Join[
				Table[nam["dil", 2 (i - 1) + 1] -> Conv[featInt, {dim, 3, 2 (i - 1) + 1}, act], {i, 1, dep}],
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


Conv[featOut_, {dim_, kern_}] := Conv[featOut, {dim, kern, 1}, "None"]

Conv[featOut_, {dim_, kern_}, act_] := Conv[featOut, {dim, kern, 1}, act]

Conv[featOut_, {dim_, kern_, dil_}] := Conv[featOut, {dim, kern, dil}, "None"]

Conv[featOut_, {dim_, kern_, dil_}, act_] := {
	(*The most basic conv block CON>BN>ACT used in each of the advanced conv blocks*)
	ConvolutionLayer[featOut, ConstantArray[kern, dim], 
		"PaddingSize" -> (Ceiling[kern/2] - 1) dil, 
		"Stride" -> 1, 
		"Dilation" -> dil],	
	BatchNormalizationLayer[], 
	ActivationLayer[act]
}


(* ::Subsubsection::Closed:: *)
(*ActivationLayer*)


ActivationLayer[] := ActivationLayer["GELU"]

ActivationLayer[actType_] := If[Head[actType]===ParametricRampLayer||Head[actType]===ElementwiseLayer,actType,
	Switch[actType, 
	"LeakyRELU", ParametricRampLayer[], 
	"None" | "", Nothing,
	_, If[StringQ[actType],ElementwiseLayer[actType], Echo["not a correct activation"];Nothing]
]]


(* ::Subsection:: *)
(*NetDimensions*)


(* ::Subsubsection::Closed:: *)
(*NetDimensions*)


NetDimensions[net_]:=NetDimensions[net, ""]

NetDimensions[net_, port_]:=Block[{block, neti},Switch[port,
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


(* ::Subsection:: *)
(*LossLayers*)


SyntaxInformation[GetNetNodes] = {"ArgumentsPattern" -> {_}};

GetNetNodes[net_]:=DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]]


(* ::Subsection:: *)
(*LossLayers*)


(* ::Subsubsection::Closed:: *)
(*AddLossLayer*)


SyntaxInformation[AddLossLayer] = {"ArgumentsPattern" -> {_}};

(*http://arxiv.org/abs/2312.05391*)
AddLossLayer[net_]:=NetGraph[<|
	"net"->net,
	"Dice" -> DiceLossLayer[2], 
	"Jaccard" -> JaccardLossLayer[2],
	"Tversky" -> TverskyLossLayer[0.7], 
	"Focal" ->  NetGraph[NetChain[{FocalLossLayer[2, 0.25], ElementwiseLayer[5 #&]}]],
	"SquaredDiff" -> NetGraph[NetChain[{MeanSquaredLossLayer[], ElementwiseLayer[50 #&]}]],
	"CrossEntropy" -> NetGraph[NetChain[{CrossEntropyLossLayer["Probabilities"]}]]
|>,{
	{"net", NetPort["Target"]} -> "Dice" -> NetPort["Dice"], (*using squared dice, F1score*)
	{"net", NetPort["Target"]} -> "Jaccard" -> NetPort["Jaccard"],(*using squared Intersection over union*)
	{"net", NetPort["Target"]} -> "Tversky" -> NetPort["Tversky"], (*recall more than precision, 0.5 is dice*)
	{"net", NetPort["Target"]} -> "SquaredDiff" -> NetPort["SquaredDiff"], (*Brier Score*)
	{"net", NetPort["Target"]} -> "Focal" -> NetPort["Focal"], (*scaled cross entropy for hard examples*)
	{"net", NetPort["Target"]} -> "CrossEntropy" -> NetPort["CrossEntropy"]
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
		{"intersection",  "target", "input"} -> "dice" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsubsection::Closed:: *)
(*JaccardLossLayer*)


SyntaxInformation[JaccardLossLayer] = {"ArgumentsPattern" -> {_., _.}};

JaccardLossLayer[] := JaccardLossLayer[2]

JaccardLossLayer[n_?IntegerQ]:= NetFlatten[
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
	(* https://doi.org/10.48550/arXiv.1706.05721 Tversky loss function for 3D segmentation*)
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
		"flatPr" -> {ThreadingLayer[#1  #2 &], AggregationLayer[Total, {-1}], FlattenLayer[]},
		"focal" -> {ThreadingLayer[-a  Log[#1 + 10^-15](1 - #1)^g &], AggregationLayer[Mean, 1]}
		}, {
			{NetPort["Input"], NetPort["Target"]} -> "flatPr" -> "focal" -> NetPort["Loss"]
	}, "Loss" -> "Real"]
]


(* ::Subsection:: *)
(*Encoders*)


(* ::Subsubsection::Closed:: *)
(*ClassEncoder*)


SyntaxInformation[ClassEncoder] = {"ArgumentsPattern" -> {_, _.}};

ClassEncoder[data_]:= ClassEncoder[data, Round[Max@data]+1]

ClassEncoder[data_, nClass_]:= ToPackedArray@Round@If[nClass === 1, data, ClassEncoderC[data, UnitVector[nClass, 1]]]

ClassEncoderC = Compile[{{class, _Integer, 0}, {vec, _Integer, 1}}, 
	RotateRight[vec, class], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*ClassDecoder*)


SyntaxInformation[ClassDecoder] = {"ArgumentsPattern" -> {_, _.}};

ClassDecoder[data_]:= ClassDecoder[data, Last@Dimensions@data]

ClassDecoder[data_, nClass_]:=ToPackedArray@Round@ClassDecoderC[data, Range[nClass]]

ClassDecoderC = Compile[{{prob, _Real, 1}, {classes, _Integer, 1}}, 
	Max[classes (1 - Unitize[Chop[(prob/Max[prob]) - 1]])] - 1, 
RuntimeAttributes -> {Listable}]


(* ::Subsection:: *)
(*Data Patching*)


(* ::Subsubsection::Closed:: *)
(*PatchesToData*)


SyntaxInformation[PatchesToData] = {"ArgumentsPattern" -> {_, _, _., _.}};

PatchesToData[patches_, location_] := PatchesToData[patches, location, Max /@ Transpose[location[[All, All, 2]]], {}]

PatchesToData[patches_, location_, dim : {_?IntegerQ, _?IntegerQ, _?IntegerQ}] := PatchesToData[patches, location, dim, {}]

PatchesToData[patches_, location_, dim : {_?IntegerQ, _?IntegerQ, _?IntegerQ}, labs_?ListQ] := Block[{
		sel, zero, a1, a2, b1, b2, c1, c2, dat, wt,
		seg, lab, pos, si, pi, overlap
	},

	If[labs === {},
		(*no labs given, assuming normal data*)
		PatchesToDataI[patches, location, dim]
		,
		(*labs given, assuming segmentations*)
		{seg, lab} = Transpose[SplitSegmentations /@ patches];
		seg = Transpose /@ seg;
		zero = SparseArray[{}, dim];

		(*get the positions for each expected label*)
		pos = Position[lab, #]& /@ labs;

		(*make the segmentation for each expected label*)
		seg = Table[If[p === {}, zero,
			si = seg[[##]] & @@ # & /@ p;
			pi = location[[#]] & /@ p[[All,1]];
			Unitize[PatchesToDataI[si, pi, dim]]
		], {p, pos}];

		(*only keep the largest connected segmentation remove the overlap and merge*)
		seg = SmoothMask[#, MaskComponents -> 1, MaskClosing -> False, SmoothIterations -> 0] &/@ seg;
		MergeSegmentations[RemoveMaskOverlaps@Transpose[seg], labs]
	]
]


PatchesToDataI[data_, loc_?MatrixQ, dim_] := Block[{dat}, 
	dat = SparseArray[{}, dim];
	(*with patch creation data can be padded with 0s on right side therefore clip here to bounds*)
	dat[[loc[[1, 1]] ;; loc[[1, 2]], loc[[2, 1]] ;; loc[[2, 2]], loc[[3, 1]] ;; loc[[3, 2]]]] = 
		data[[
			1 ;; loc[[1, 2]] - loc[[1, 1]] + 1, 
			1 ;; loc[[2, 2]] - loc[[2, 1]] + 1, 
			1 ;; loc[[3, 2]] - loc[[3, 1]] + 1
		]]; 
	dat
]

PatchesToDataI[data_, loc_?ArrayQ, dim_] := Block[{dat, wt},
	If[Length@data === 1,
		PatchesToDataI[data[[1]], loc[[1]], dim],
		dat = MapThread[PatchesToDataI[#1, #2, dim] &, {data, loc}];
		wt = N@Total@Unitize@dat;
		dat = N@Total@dat;
		SparseArray[dat["ExplicitPositions"] -> dat["ExplicitValues"]/wt["ExplicitValues"], dim]
	]
]


(* ::Subsubsection::Closed:: *)
(*DataToPatches*)


Options[DataToPatches] = {
	PatchNumber->0, 
	PatchPadding->0
}

SyntaxInformation[DataToPatches] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, opts:OptionsPattern[]]:=DataToPatches[dat, patch, "All", opts] 

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}]:={GetPatch[dat, patch, pts], pts}

DataToPatches[dat_, pts:{{{_,_},{_,_},{_,_}}..}]:={GetPatch[dat, pts], pts}

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

GetPatch[dat_, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}]:=ToPackedArray@dat[[i1;;i2,j1;;j2,k1;;k2]]

GetPatch[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}]:=GetPatch[dat, patch, #]&/@pts

GetPatch[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}]:=ToPackedArray@PadRight[dat[[i1;;i2,j1;;j2,k1;;k2]], patch, 0.]


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
(*MakeClassify*)


(* ::Subsubsection::Closed:: *)
(*MakeClassifyNetwork*)


Options[MakeClassifyNetwork] = {
	ImageSize -> {128, 128}
};

SyntaxInformation[MakeClassifyNetwork] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

MakeClassifyNetwork[classes_, OptionsPattern[]]:=Block[{enc, dec, net,imSize},
	imSize = OptionValue[ImageSize];

	enc = NetEncoder[{"Class", classes, "IndicatorVector"}];
	dec = NetDecoder[{"Class", classes}];

	net = NetChain[{
		ConvolutionLayer[16, 7, "Stride"->1, PaddingSize->3], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[4,4],
		ConvolutionLayer[32, 5, "Stride"->1, PaddingSize->2], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[4,4],
		ConvolutionLayer[64, 3, "Stride"->1, PaddingSize->1], BatchNormalizationLayer[], ElementwiseLayer["GELU"], PoolingLayer[4,4],
		FlattenLayer[], LinearLayer[128], BatchNormalizationLayer[], ElementwiseLayer["GELU"], LinearLayer[64],
		BatchNormalizationLayer[], ElementwiseLayer["GELU"], LinearLayer[Length@classes], SoftmaxLayer[]
	}, "Input" -> Prepend[imSize, 1]];

	NetFlatten@NetChain[{net}, "Input"->NetEncoder[{"Image", imSize, ColorSpace->"Grayscale"}], "Output"->dec]
]


(* ::Subsubsection::Closed:: *)
(*MakeClassifyImage*)


Options[MakeClassifyImage]={
	ImageSize->{128, 128}
};

SyntaxInformation[MakeClassifyImage] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

MakeClassifyImage[dat_, opts:OptionsPattern[]]:=Switch[ArrayDepth[dat],
	2, MakeClassifyImage[dat, opts],
	3, MakeClassifyImage[#, opts]&/@dat,
	_, $Failed
]

MakeClassifyImage[dat_?MatrixQ, OptionsPattern[]]:=Block[{imSize},
	imSize = OptionValue[ImageSize];
	If[Total[Flatten[dat]]<10,
		Image@ConstantArray[0., imSize],
		ImageResize[Image[Rescale[First[AutoCropData[{dat}, CropPadding -> 0][[1]]]]], imSize]
	]
];


(* ::Subsection:: *)
(*ClassifyData*)


(* ::Subsubsection::Closed:: *)
(*ClassifyData*)


Options[ClassifyData] = {
	TargetDevice -> "GPU"
};

SyntaxInformation[ClassifyData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ClassifyData[dat_, met_, OptionsPattern[]]:=Block[{
		data, dev, len, ran, datF, kneeStart, kneeEnd, pos
	},
	data = MakeClassifyImage[dat];
	dev = OptionValue[TargetDevice];
	Which[
		FileExistsQ[met] && FileExtension[met]==="wlnet", Import[met][data],
		StringQ[met], Switch[met,
			"LegSide", FindLegSide[data, dev],
			"LegPosition", FindLegPos[data ,dev]],
		Head[met]===NetChain || Head[met]===NetGraph, met[data]
	]
]


(* ::Subsubsection::Closed:: *)
(*FindLegSide*)


FindLegSide[data_, dev_]:=Block[{net, imSize},
	net = GetNeuralNet["LegSide"];
	If[net === $Failed, $Failed,
		imSize = NetDimensions[NetReplacePart[net, "Input"->None], "Input"][[2;;]];
		If[!AllTrue[ImageDimensions/@data, #===imSize&], $Failed,
			Last@Keys@Sort@Counts[net[data, TargetDevice -> dev]]
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*FindLegPos*)


FindLegPos[data_, dev_]:=Block[{
		net, len, ran, datF, kneeStart, kneeEnd, pos, imSize
	},
	net = GetNeuralNet["LegPosition"];
	If[net === $Failed, $Failed,
		imSize = NetDimensions[NetReplacePart[net,"Input"->None],"Input"][[2;;]];
		If[!AllTrue[ImageDimensions/@data,#===imSize&],	$Failed,
			len = Length[data];
			ran = Range[1, len];
			(*find loc per slice*)
			datF = MedianFilter[net[data, TargetDevice -> dev]/.Thread[{"Lower","Knee","Upper"}->{1.,2.,3.}], 1];
			{kneeStart, kneeEnd} = First[SortBy[Flatten[
					Table[{a, b, PosFunc[a, b, len, ran, datF]}, {a, 0, len}, {b, a+1, len}]
				, 1],Last]][[1;;2]];
			pos = Which[
				kneeStart == 0. && kneeEnd == len, "Knee",
				kneeStart > 0 && kneeEnd >= len, "Lower",
				kneeStart == 0. && kneeEnd =!= len, "Upper",
				kneeStart =!= 0. && kneeEnd =!= len, "Both"
			];
			{pos, {kneeStart + 1, kneeEnd}}
		]
	]
]


PosFunc = Compile[{{a, _Integer, 0},{b, _Integer, 0},{l, _Integer, 0},{x, _Integer, 1},{d, _Real, 1}},
	Total[((Which[1<=#<=a, 1., a<=#<=b, 2., b<=#<=l, 3., True, 0] &/@ x) - d)^2]];


(* ::Subsection:: *)
(*SegmentData*)


(* ::Subsubsection::Closed:: *)
(*SegmentData*)


Options[SegmentData] = {
	TargetDevice -> "GPU", 
	MaxPatchSize->Automatic, 
	Monitor->False
};

SyntaxInformation[SegmentData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SegmentData[data_, opts:OptionsPattern[]]:=SegmentData[data, "Legs", opts]

SegmentData[data_, what_, OptionsPattern[]] := Block[{
		dev, max, mon, patch, pts, dim ,loc, set, net, type, segs, all, time, timeAll
	},

	timeAll = First@AbsoluteTiming[
		{dev, max, mon} = OptionValue[{TargetDevice, MaxPatchSize, Monitor}];

		(*split the data in upper and lower legs and left and right*)
		If[mon, Echo[Dimensions@data, "Analyzing the data with dimensions:"]];
		time = First@AbsoluteTiming[
			{{patch, pts, dim}, loc, set} = SplitDataForSegmentation[data, what, Monitor -> mon, TargetDevice -> dev]
		];
		If[mon, Echo[Round[time, .1], "Total time for analysis [s]: "]];
		If[mon, Echo[Column@Thread[{loc,Dimensions/@ patch}], "Segmenting \""<>what<>"\" locations with dimensions:"]];

		
		(*get the network name and data type*)
		{net, type} = Switch[what,
			"LegBones", {"SegLegBones"&, "Bones"},
			"Legs"|"UpperLegs"|"LowerLegs",	{
				(#[[1]] /. {"Upper" -> "SegThighMuscle", "Lower" -> "SegLegMuscle"})&, "Muscle"
			},
			_, Return[]
		];

		(*Perform the segmentation*)
		time = First@AbsoluteTiming[segs = MapThread[(
			If[mon, Echo[{#2, net[#2]}, "Performing segmentation for: "]];
			segs = ApplySegmentationNetwork[#1, net[#2], TargetDevice -> dev, MaxPatchSize->max, Monitor->mon];
			ReplaceLabels[segs, #2, type]
		) &, {patch, loc}]];
		If[mon, Echo[Round[time, .1], "Total time for segmentations [s]: "]];

		(*Merge all segmentations for all expected labels*)
		all = Select[DeleteDuplicates[Sort[Flatten[GetSegmentationLabels/@segs]]], IntegerQ];
		If[mon, Echo[all, "Putting together the segmentations with labels"]];

		(*after this only one cluster per label remains*)
		time = First@AbsoluteTiming[segs = PatchesToData[segs, pts, dim, all]];
		If[mon, Echo[Round[time, .1], "Total time for final evaluation [s]: "]];
	];
	If[mon, Echo[Round[timeAll, .1], "Total evaluation time [s]: "]];

	segs
]


(* ::Subsubsection::Closed:: *)
(*ReplaceLabels*)


ReplaceLabels[seg_, loc_, type_] := Block[{what, side, labNam, labIn, labOut, file},
	{what, side} = loc;
	labIn = GetSegmentationLabels[seg];

	(*for now overwrite the labIn with custom values since some muscles are not segmented *)
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
(*SplitDataForSegmentation*)


Options[SplitDataForSegmentation] = {
	Monitor -> False,
	TargetDevice -> "GPU"
};

SyntaxInformation[SplitDataForSegmentation] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};


SplitDataForSegmentation[data_?ArrayQ, seg_?ArrayQ, opt:OptionsPattern[]]:=SplitDataForSegmentation[data, seg, "Legs", opt]

SplitDataForSegmentation[data_?ArrayQ, seg_?ArrayQ, what_?StringQ, opt:OptionsPattern[]]:=Block[{dat,pts,dim,loc,set, segp},
	{{dat, pts, dim}, loc, set} = SplitDataForSegmentation[data, what, opt];
	segp = GetPatch[seg, pts];
	{{dat, pts, dim}, {segp, pts, dim}, loc, set}
]


SplitDataForSegmentation[data_?ArrayQ, opt:OptionsPattern[]]:=SplitDataForSegmentation[data, "Legs", opt]

SplitDataForSegmentation[data_?ArrayQ, what_?StringQ, opt:OptionsPattern[]]:=Block[{
		dim, whatSide, side, whatPos, pos, dat, right, left, cut, pts, loc, time, mon, dev
	},
	dim = Dimensions[data];
	{mon, dev} = OptionValue[{Monitor, TargetDevice}];

	Switch[what,
		"Legs"|"UpperLegs"|"LowerLegs",
		(*split the data in upper and lower legs and left and right*)

		(*find which side using NN*)
		time= First@AbsoluteTiming[whatSide = ClassifyData[data, "LegSide", TargetDevice -> dev]];
		If[mon, Echo[whatSide, "Data contains sides: "]];
		If[mon, Echo[Round[time, .1], "Time for side estimation [s]:"]];

		(*based on side cut data or propagate*)
		dat = Switch[whatSide,
			(*both sides which need to be split*)
			"Both",
			{right, left, cut} = CutData[data];
			{
				{right, {"Right", {1, cut}}}, 
				{left, {"Left", {cut+1, dim[[3]]}}}
			},
			_,
			(*only one side, no split*)
			cut=0; 
			{
				{data, {whatSide, {1, dim[[3]]}}}
			}
		];

		(*loop over data to find upper or lower*)
		time= First@AbsoluteTiming[dat = Flatten[(
			{dat, side} = #;
			{whatPos, pos} = Switch[what,
				"Legs", ClassifyData[dat, "LegPosition", TargetDevice -> dev],
				"UpperLegs", {"Upper", dim[[1]]},
				"LowerLegs", {"Lower", dim[[1]]}
			];

			Switch[whatPos,
				(*if upper and lower split upper and lower*)
				"Both", {
					{dat[[pos[[1]];;]], {"Upper", {pos[[1]], dim[[1]]}}, side}, 
					{dat[[;;pos[[2]]]], {"Lower", {1, pos[[2]]}}, side}
				},
				(*if only knee data duplicate for both networks*)
				"Knee", {
					{dat, {"Upper", {1, dim[[1]]}}, side}, 
					{dat, {"Lower", {1, dim[[1]]}}, side}
				},
				(*if only upper or only lower return what it is*)
				_, {
					{dat, {whatPos, {1, dim[[1]]}}, side}
				}
			]
		)&/@dat, 1]];

		If[mon, Echo[whatPos, "Data contains positions: "]];
		If[mon, Echo[Round[time, .1], "Time for position estimation [s]:"]];

		(*output the selected data with the correct label and coordinates*)
		{dat, pts, loc} = Transpose[CropPart/@dat];
		{{dat, pts, dim}, loc, {{whatSide, cut}, {whatPos, pos}}}

		,
		_,
		$Failed		
	]
]


(* ::Subsubsection::Closed:: *)
(*CropPart*)


CropPart[{dat_, {up_, {upStart_, upEnd_}}, {side_, {sideStart_, sideEnd_}}}]:= Block[{data, crp},
	{data, crp} = AutoCropData[dat, CropPadding->0];
	{data, Partition[crp, 2] + {upStart-1, 0, sideStart-1}, {up, side}}
]


(* ::Subsubsection::Closed:: *)
(*ApplySegmentationNetwork*)


Options[ApplySegmentationNetwork] = {
	TargetDevice->"GPU", 
	DataPadding->0, 
	MaxPatchSize->Automatic, 
	Monitor->False
}

SyntaxInformation[ApplySegmentationNetwork] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};


(*Apply segmentation network on folder of datasets*)
ApplySegmentationNetwork[inp_?(!(TensorQ[#, NumericQ] || StringQ[#])&), rest___] := Block[{
		dat, datFol, outFol, inLab, outLab, i, rule, files, netI, node, seg, vox, im
	},
	(*Figure out what the input is *)
	{{datFol, outFol}, {inLab, outLab}, i} = Which[
		Length[inp] == 2, Which[
			VectorQ[inp, StringQ], {inp, {"data", "label_NN"}, 1},
			VectorQ[First@inp, StringQ] && VectorQ[Last@inp, StringQ], Append[inp, 1],
			VectorQ[First@inp, StringQ] && IntegerQ[Last@inp], {First@inp, {"data", "label_NN"}, Last@inp},
			True, Return[$Failed]
		],
		Length[inp] === 3, inp,
		True, Return[$Failed]
	];

	(*get the files from the input folder and define rule for output folder*)
	rule = {"_" <> inLab -> "_" <> outLab, datFol -> outFol};
	files = FileNames["*" <> inLab <> ".nii.gz", datFol][[i ;;]];

	(*loop over actual apply segment function for all files*)
	Table[EchoTiming[
		{dat, vox} = ImportNii[f];
		seg = ApplySegmentationNetwork[dat, rest];
		ExportNii[seg, vox, StringReplace[f, rule]];
		im = ImageResize[MakeChannelClassGrid[{dat}, seg, 5], 1200];
		Export[StringReplace[f, Flatten@{rule, ".nii.gz" -> ".png"}], im, "ColorMapLength" -> 256];
	, f], {f, files}]
]


(*Apply segmentation network on single dataset*)
ApplySegmentationNetwork[dat_, netI_, opt:OptionsPattern[]]:=ApplySegmentationNetwork[dat, netI, "", opt]

ApplySegmentationNetwork[dat_, netI_, node_, OptionsPattern[]]:=Block[{
		dev, pad , lim, data, crp, net, dim, ptch, time, nodes,
		patch, pts, seg, mon, dimi, dimc, lab, nClass, prec, normF
	},

	{dev, pad, lim, mon} = OptionValue[{TargetDevice, DataPadding, MaxPatchSize, Monitor}];
	If[lim === Automatic, lim = If[dev==="GPU", 200, 224]];
	prec = If[(dev==="GPU") && ($OperatingSystem === "Windows"), "Mixed", "Real32"];

	data = Which[
		StringQ[dat], First@ImportNii[dat], 
		TensorQ[{dat}, NumericQ], If[First@Dimensions@dat===1, First@dat, dat],
		True, Return[$Failed]
	];

	dimi = Dimensions@data;
	{data, crp} = AutoCropData[data, CropPadding->0];
	dimc = Dimensions@data;
	data = N@ArrayPad[data, pad, 0.];
	dim = Dimensions[data];

	If[mon, Echo[{"in", dimi, "crop", dimc, "pad", dim}, "Data dimensions before and after cropping and padding are: "]];

	net = If[StringQ[netI], GetNeuralNet[netI], netI];

	If[net===$Failed, $Failed,
		(*calculate the patch size for the data *)
		ptch = FindPatchDim[net, dim, lim];

		(*create the patches*)
		{patch, pts} = DataToPatches[data, ptch, PatchNumber -> 0, PatchPadding->pad];
		If[mon, Echo[{ptch, Length@patch}, "Patch size and created number of patches is:"]];

		(*create the network*)
		net = ChangeNetDimensions[net, "Dimensions" ->ptch];
		nClass = NetDimensions[net,"Output"][[-1]];

		(*dataNormalization function*)
		normF = NormalizeData[#, NormalizeMethod -> "Uniform"]&;

		(*perform the segmentation*)
		If[node==="",
			time = First@AbsoluteTiming[
				(*actually perform the segmentation with the NN*)
				seg = ClassDecoder[net[{normF[#]}, TargetDevice->dev, WorkingPrecision ->prec]]&/@patch;
				(*reverse all the padding and cropping and merged the patches if needed*)
				seg = ReverseCrop[ArrayPad[
						PatchesToData[ArrayPad[#, -pad] & /@ seg, Map[# + {pad, -pad} &, pts, {2}], 
							dim, Range[nClass]]
					, -pad], dimi, crp];
			];
			If[mon, 
				Echo[{Dimensions[seg], Sort@Round@DeleteDuplicates[Flatten[seg]]}, "Output segmentations dimensions and labels:"];
				Echo[Round[time, .1], "Time for segmentation [s]: "]
			];

			,
			(*perform the segmentation on a specific node*)
			(*check if node is part of the network*)
			nodes = DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]];
			If[!MemberQ[nodes, node],
				Echo["The node "<>node<>" is not part of the network"];
				Echo[nodes];
				Return[$Failed]
				,
				seg = NetTake[net, node][{normF[First@patch]}, TargetDevice->dev, WorkingPrecision ->prec];
				If[Head[seg] === Association, seg = Last@seg];
				If[node == nodes[[-1]], seg = RotateDimensionsRight[seg]];
			];
		];

		(*give the output*)
		seg
	]
]


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
	class = NetDimensions[net, "Output"][[-1]];

	(*check smallest net dimensions output*)
	out = Rest[NetDimensions[net, "MinEncodingOut"]];

	(*needed scaling*)
	sc = inp/out;
	dimM = sc  Ceiling[dim/sc];
	(*dimM = sc ({Floor[#[[1]]], Ceiling[#[[2]]], Ceiling[#[[3]]]} & [dim/sc]);*)

	(*if memory limit is given find patch that fits*)
	If[!(lim === 1000 || lim === Automatic),
		If[CubeRoot[N[Times @@ dimM]] < lim, 
			dimN = dimM, 
			u = 1;
			cont = True;
			While[cont, u++;
				dimN = u  sc;
				dimN = Min /@ Transpose[{dimN, dimM}];
				cont = CubeRoot[Times @@ dimN] < lim && Min [dimN] < Min[dim]]
		],
		dimN = dimM
	];

	(*output the patch dim*)
	{Min@#[[1]], Max@#[[2]], Max@#[[3]]} & [Thread[{dimN, {2, 1, 1} inp}]]
]


(* ::Subsection:: *)
(*TrainSegmentationNetwork*)


(* ::Subsubsection::Closed:: *)
(*TrainSegmentationNetwork*)


Options[TrainSegmentationNetwork] = {
	LoadTrainingData -> True,
	MonitorInterval -> 1,

	PatchSize -> {32, 112, 112},
	PatchesPerSet -> 1,
	BatchSize -> 4,
	RoundLength -> 512,
	MaxTrainingRounds -> 150,

	BlockType -> "ResNet",
	NetworkArchitecture -> "UNet",
	ActivationType -> "GELU",
	DownsampleSchedule -> 2,
	SettingSchedule -> Automatic,
	FeatureSchedule -> 32,
	NetworkDepth -> 5,

	MultiChannel -> False,

	AugmentData -> True,
	PadData-> False,

	LossFunction -> All,
	DropoutRate -> 0.2,
	LearningRate -> 0.001,
	L2Regularization -> 0.0001,

	MonitorCalc -> False
}


SyntaxInformation[TrainSegmentationNetwork] = {"ArgumentsPattern" -> {{_, _}, _., OptionsPattern[]}};

TrainSegmentationNetwork[{inFol_?StringQ, outFol_?StringQ}, opts : OptionsPattern[]] := TrainSegmentationNetwork[{inFol, outFol}, "Start", opts]

TrainSegmentationNetwork[{inFol_?StringQ, outFol_?StringQ}, netCont_, opts : OptionsPattern[]] := Block[{
		netOpts, batch, roundLength, rounds, data, depth, nChan, nClass, outName, ittString, multi,
		patch, augment, netIn, ittTrain, testData, testVox, testSeg, im, patches,
		monitorFunction, netMon, netOut, trained, l2reg, pad, batchFunction,
		validation, files, loss, rep, learningRate, schedule
	},


	(*------------ Get all the configuration stuff -----------------*)

	(*getting all the options*)
	netOpts = Join[FilterRules[{opts}, Options@MakeUnet], FilterRules[Options@TrainSegmentationNetwork, 
		Options@MakeUnet]];

	{batch, roundLength, rounds, augment, pad, patch, patches, 
		loss, rep, learningRate, l2reg, multi} = OptionValue[
		{BatchSize, RoundLength, MaxTrainingRounds, AugmentData, PadData, PatchSize, PatchesPerSet, 
			LossFunction, MonitorInterval, LearningRate, L2Regularization, MultiChannel}];
	pad = If[NumberQ[pad], Round[pad], False];

	(*get the train data files*)
	files = FileNames["*.wxf", inFol];

	(*figure out network properties from train data*)
	testData = Import@First@files;
	(*figure out how to treat multi channel data*)	
	depth = ArrayDepth@testData;
	nChan = If[depth === 4 && multi, Length@First@testData, 1];
	(*background is 0 but for network its 1 so class +1*)
	nClass = Round[Max@testData[[2]] + 1];


	(*------------ Define the network -----------------*)

	Echo[DateString[], "Preparing the network"];

	(*make or import network, netCont can be a network or a previous train folder*)
	{netIn, ittTrain} = Which[
		(*start with a clean network*)
		netCont === "Start",
		{MakeUnet[nChan, nClass, patch, MonitorCalc->OptionValue[MonitorCalc], Sequence@netOpts], 0}
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

	(*quit if there are not enough rounds*)
	If[rounds - ittTrain < 5, Return[Message[TrainSegmentationNetwork::itt]; $Failed]];

	(*define and check the training loss function*)
	loss = Which[
		loss === All, {"Dice", "SquaredDiff", "Tversky" , "CrossEntropy", "Jaccard", "Focal"},
		StringQ[loss], {loss},
		True, loss];
	If[!And @@ (MemberQ[{"Dice", "SquaredDiff", "Tversky", "CrossEntropy", "Jaccard", "Focal"}, #] & /@ loss), 
		Return[Message[TrainSegmentationNetwork::loss]; $Failed]];

	(*if the network already exists make the dimensions, classes en channels match the input*)
	netIn = NetInitialize@ChangeNetDimensions[netIn, "Dimensions" -> patch, "Channels" -> nChan, "Classes" -> nClass];
	Echo[NetSummary[netIn, "Mem"], "Network summary: "];

	(*define the network for training*)
	netIn = AddLossLayer@netIn;
	Echo[netIn, "Network is ready"];


	(*---------- Training functions ----------------*)

	(*Local functions*)
	outName = FileNameJoin[{outFol, Last[FileNameSplit[outFol]] <> "_" <> #}]&;
	ittString = "itt_" <> StringPadLeft[ToString[#], 4, "0"]&;

	(*Monitor function*)
	monitorFunction = (
		ittTrain++;
		(*perform segmentation and export*)
		netMon = NetExtract[#Net, "net"];
		testSeg = Ramp[ClassDecoder[netMon[testData, TargetDevice -> "GPU", WorkingPrecision -> "Mixed"]]];
		ExportNii[testSeg, testVox, outName[ittString[ittTrain]<>".nii"]];
		(*make and export test image*)
		im = MakeChannelClassGrid[testData, {testSeg, {0, nClass-1}}, 3];
		Export[outName[ittString[ittTrain] <> ".png"], im , "ColorMapLength" -> 256];
		(*export the network and delete the one from the last iteration*)
		Export[outName[ittString[ittTrain] <> ".wlnet"], netMon];
		Quiet@DeleteFile[outName[ittString[ittTrain - 1] <> ".wlnet"]];
	)&;

	(*batch function*)
	batchFunction = GetTrainData[data, #BatchSize, patch, nClass, 
		PatchesPerSet -> patches, AugmentData -> augment, PadData -> pad] &;

	(*oneCycle learning rate schedual*)
	schedule = With[{n1 = 0.2 #2, n2 = 0.25 #2, n3 = 0.9 #2, sc = 0.3},	((1 - sc) Which[
			#1 < n1, Cos[Pi (#1 - n1)/n1], 
			#1 < n2, 1, 
			#1 < n3, Cos[Pi (#1 - n2)/(n3 - n2)], 
			True, -1] + 1 + sc
		)/2 ] &;	


	(*---------- Prepare the data ----------------*)
	
	(*Make and export test data*)
	{testData, testVox} = MakeTestData[testData, 2, patch];
	ExportNii[First@testData, testVox, outName["testSet.nii"]];

	(*export first itt*)
	ittTrain--;
	monitorFunction[<|"Net"->netIn|>];
	
	Echo[DateString[], "Preparing the data"];
	(*import all train data or train out of memory*)
	data = If[OptionValue[LoadTrainingData] === True, Import /@ files, files];

	(*prepare a validation set which is 10% of round*)
	validation = batchFunction[<|"BatchSize" -> Round[0.1 roundLength]|>];
	Echo[{Length@data, Length@validation}, "data / validation: "];


	(*---------- Train the network ----------------*)

	(*Print progress function*)
	Echo[{DateString[], loss}, "Starting training"];
	PrintTemporary[Dynamic[Column[{
		Style["Training Round: " <> ToString[ittTrain], Bold, Large], 
		Image[im, ImageSize->400]
	}, Alignment -> Center]]];

	(*train the network*)
	trained = NetTrain[netIn, {batchFunction, "RoundLength" -> roundLength}, All, 

		ValidationSet -> validation,
		LossFunction -> loss,
		TargetDevice -> "GPU", WorkingPrecision -> "Mixed",

		MaxTrainingRounds -> rounds - ittTrain, BatchSize -> batch, LearningRate -> learningRate, 
		Method -> {"ADAM", "Beta1" -> 0.9, "Beta2" -> 0.999, "Epsilon" -> 10^-5, 
			"L2Regularization" -> l2reg, "LearningRateSchedule" -> schedule},

		TrainingProgressFunction -> {monitorFunction, "Interval" -> Quantity[rep, "Rounds"]},
		TrainingProgressReporting -> File[outName[StringReplace[DateString["ISODateTime"], ":" | "-" -> ""] <> ".json"]]
	];

	(*---------- Export the network ----------------*)

	netOut = NetExtract[trained["TrainedNet"], "net"];
	Export[outName["trained" <> ".wxf"], trained];
	Export[outName["final" <> ".wlnet"], netOut];
	Export[outName["final" <> ".onnx"], netOut];
]


(* ::Subsubsection::Closed:: *)
(*MakeTestData*)


MakeTestData[data_, n_, patch_] := Block[{testData, len, sel, testDat},
	testData = data[[1]];
	(*figure out how to treat multi channel data - for now just take the first volume*)
	If[ArrayDepth@testData === 4, testData = testData[[All, 1]]];

	(*crop the test data*)
	len = Length@testData;
	If[len > First@patch,
		sel = Range @@ Clip[Round[(Clip[Round[len/3 - (0.5 n) First@patch], {0, Infinity}] + {1, n First@patch})], 
			{1, len}, {1, len}];
		testData = First@AutoCropData[testData[[sel]]]
	];

	{{NormalizeData[PadToDimensions[testData, patch], NormalizeMethod -> "Uniform"]}, data[[3]]}
];


(* ::Subsection:: *)
(*AugmentTrainingData*)


(* ::Subsubsection::Closed:: *)
(*AugmentTrainingData*)


SyntaxInformation[AugmentTrainingData] = {"ArgumentsPattern" -> {_, _, _.}};

AugmentTrainingData[{dat_?ArrayQ, seg_?ArrayQ}, vox_] := AugmentTrainingDataI[{dat, seg}, vox, {True, True, True, True, True}]

AugmentTrainingData[{dat_?ArrayQ, seg_?ArrayQ}, vox_, aug_?BooleanQ] := AugmentTrainingDataI[{dat, seg}, vox, {aug, aug, aug, aug, aug}]

AugmentTrainingData[{dat_?ArrayQ, seg_?ArrayQ}, vox_, aug_?ListQ] := AugmentTrainingDataI[{dat, seg}, vox, aug]

AugmentTrainingData[dat_?ArrayQ, vox_] := First@AugmentTrainingDataI[{dat, dat}, vox, {True, True, True, True, True}]

AugmentTrainingData[dat_?ArrayQ, vox_, aug_?BooleanQ] := First@AugmentTrainingDataI[{dat, dat}, vox, {aug, aug, aug, aug, aug}]

AugmentTrainingData[dat_?ArrayQ, vox_, aug_?ListQ] := First@AugmentTrainingDataI[{dat, dat}, vox, aug]


AugmentTrainingDataI[{dat_?ArrayQ, seg_?ArrayQ}, vox_, aug_?ListQ] := Block[{
	datT, segT, cr, w, r, t, s, flip, rot, trans, scale, noise, blur},

	(*prep data*)
	{flip, rot, scale, noise, blur} = aug;
	datT = ToPackedArray[N[dat]];
	segT = ToPackedArray[N[seg]];

	(*Augmentations sharpness*)
	If[blur && Coin[], datT = GaussianFilter[datT, RandomReal[{0, 2}]]];
	(*Augmentation of noise*)
	If[noise && Coin[], datT = AddSaltAndRice[datT, RandomReal[{5, 50}], CoinN[] RandomReal[{0.001, 0.01}]]];
	(*Augmentation of mirroring*)
	If[flip && Coin[], {datT, segT} = ReverseC[{datT, segT}]];
	(*Augmentation of orientation and scale*)
	If[rot || scale,
		w = {
			(*rotation around z, y and x axis, z is in-plane for axial*)	
			If[rot && Coin[], RandomReal[{-30, 30}], 0.],
			If[rot && Coin[], RandomReal[{-15, 15}], 0.], 
			If[rot && Coin[], RandomReal[{-15, 15}], 0.],
			(*no Translation, never needed since how data is padded and cropped its always centered*)
			0., 0., 0.,
			(*Scaling*)
			If[scale && Coin[], RandomReal[{0.6, 1.6}], 1.],
			If[scale && Coin[], RandomReal[{0.6, 1.6}], 1.],
			If[scale && Coin[], RandomReal[{0.6, 1.6}], 1.],
			(*no skewing*)
			0., 0., 0.
		};
		(*actual transform the data*)
		{datT, segT} = DataTransformation[#, vox, w, InterpolationOrder -> 0, 
			PadOutputDimensions -> True]& /@ {datT, segT};
		(*make sure there is no unneeded background before rest of processing*)
		cr = FindCrop[datT, CropPadding -> 0];
		{datT, segT} = ApplyCrop[#, cr]& /@ {datT, segT};
	];
	
	(*output augmented data*)
	{ToPackedArray[N[datT]], ToPackedArray[Round[segT]]}
]


(* ::Subsubsection::Closed:: *)
(*Coin*)


Coin[] := Coin[0.5];
Coin[t_] := RandomChoice[{t, 1 - t} -> {True, False}]

CoinN[] := CoinN[0.5];
CoinN[t_] := RandomChoice[{t, 1 - t} -> {1, 0}]


(* ::Subsubsection::Closed:: *)
(*ReverseC*)


ReverseC = Compile[{{dat, _Real, 1}}, Reverse[dat], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*AddSaltAndRice*)


AddSaltAndRice[data_, snr_, p_] := Block[{dims, sp, sigma, noise, coors},
	dims = Dimensions[data];
	sp = SparseArray[N@data];
	(*random noise for rice distribution*)
	sigma = Mean[sp["ExplicitValues"]] / (snr + 1);
	noise = RandomReal[NormalDistribution[0., sigma], dims];
	(*salt and pepper noise locations*)
	coors = sp["ExplicitPositions"];
	coors = RandomSample[coors, 2 Max[{1, Round[0.5 p Length[coors]]}]];
	(*add noise*)
	SaltAndRiceC[data, noise, coors]
];


SaltAndRiceC = Compile[{{data, _Real, 3}, {noise, _Real, 3}, {coors, _Integer, 2}}, Block[{newData, num},
    newData = Unitize[data] Sqrt[(data + noise)^2. + RandomSample[noise]^2.];
    num = Round[Length[coors]/2];
    Do[
		newData[[coors[[i, 1]], coors[[i, 2]], coors[[i, 3]]]] = 1.;
		newData[[coors[[i + num, 1]], coors[[i + num, 2]], coors[[i + num, 3]]]] = 0.;
	, {i, 1, num}];
	newData
], RuntimeAttributes -> {Listable}, Parallelization -> True];


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


(* ::Subsection:: *)
(*Get Train Data*)


(* ::Subsubsection::Closed:: *)
(*GetTrainData*)


Options[GetTrainData] = {
	PatchesPerSet -> 1, 
	AugmentData -> True,
	PadData -> False
};

SyntaxInformation[GetTrainData] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

GetTrainData[dataSets_, nBatch_, patch_, opts:OptionsPattern[]]:=GetTrainData[dataSets, nBatch, patch, False, opts]

GetTrainData[dataSets_, nBatch_, patch_, nClass_, OptionsPattern[]] := Block[{
		itt, datO, segO, dat, seg, vox, augI, aug, nSet, pad
	},

	itt = 0;
	datO = segO = {};

	{augI, nSet, pad} = OptionValue[{AugmentData, PatchesPerSet, PadData}];
	aug = If[BooleanQ[augI], augI, True];

	(*get the correct number of sets*)
	itt = Ceiling[nBatch/nSet];
	Do[
		(*get random dataset*)
		dat = RandomChoice[dataSets];
		{dat, seg, vox} = Which[
			(*data is wxf file format*)
			StringQ[dat], Import[dat],
			(*dataSets is list of nii files {dat.nii, seg.nii}*)
			Length[dat]===2, Join[ImportNii[dat[[1]]][[{1}]], ImportNii[dat[[2]]]],
			(*data is in memory*)
			True, dat
		];

		(*figure out how to treat multichannel - for now take a random volume*)
		dat = If[ArrayDepth[dat] === 4, RandomChoice@Transpose@dat, dat];

		(*perform augmentation on full data and get the defined number of patches*)
		{dat, seg} = AugmentTrainingData[{dat, seg}, vox, aug];
		{dat, seg} = PatchTrainingData[{dat, seg}, patch, nSet];
		datO = Join[datO, dat];
		segO = Join[segO, seg];
	, itt];

	(*make the output*)
	datO = datO[[;; nBatch]];
	segO = segO[[;; nBatch]];
	If[IntegerQ[pad], {datO, segO} = AddPadding[datO, segO, pad]];
	datO = NormalizeData[#, NormalizeMethod -> "Uniform"]& /@ datO;
	segO = If[IntegerQ[nClass], ClassEncoder[segO, nClass], segO + 1];
	Thread[Transpose[{ToPackedArray@N@datO}] -> ToPackedArray@Round@segO]
];


AddPadding[dat_, seg_, p_]:=Block[{datp, segp, pad},
	Transpose@MapThread[(
		datp = #1;
		segp = #2;
		If[RandomChoice[{0.3, 0.7} -> {True, False}], 
			pad = RandomInteger[{-p, p}];
			Which[
				pad < 0, datp[[pad ;;, All, All]] = segp[[pad ;;, All, All]] = 0.,
				pad > 0, datp[[;; pad, All, All]] = segp[[;; pad, All, All]] = 0.
			]
		];
		{datp, segp}
	)& ,{dat, seg}]
]


(* ::Subsubsection::Closed:: *)
(*PatchTrainingData*)


PatchTrainingData[{dat_, seg_}, patch_, n_]:=Block[{pts,datP,segP},
	(*by overlapping patches with PatchNumber more random patches per data are created*)
	{datP, pts} = DataToPatches[dat, patch, n, PatchNumber->2];
	segP = DataToPatches[seg, patch, pts][[1]];
	{ToPackedArray[N[#]]&/@datP, ToPackedArray[Round[#]]&/@segP}
]


(* ::Subsection:: *)
(*PrepareTrainingData*)


(* ::Subsubsection::Closed:: *)
(*PrepareTrainingData*)


Options[PrepareTrainingData] = {
	LabelTag -> "label",
	DataTag -> "data",
	InputLabels -> Automatic,
	OutputLabels -> Automatic,
	TrainVoxelSize -> Automatic,
	CleanUpSegmentations -> True,
	TestRun -> False
}

SyntaxInformation[PrepareTrainingData] = {"ArgumentsPattern" -> {_, _,OptionsPattern[]}};

PrepareTrainingData[labFol_?StringQ, outFol_?StringQ, opt:OptionsPattern[]]:=PrepareTrainingData[{labFol, labFol}, outFol, opt]

PrepareTrainingData[{labFol_?StringQ, datFol_?StringQ}, outFol_?StringQ, OptionsPattern[]] := Block[{
		labT, datT, inLab, outLab, test, segFiles, datFiles, name, i, df, voxOut, dimSeg, dimDat,
		seg, err, vox, voxd, dat, im, nl, outf, gr, clean, legend, head, out
	},

	{labT, datT, inLab, outLab, test, clean, voxOut} = OptionValue[{LabelTag, DataTag, InputLabels, OutputLabels, 
		TestRun, CleanUpSegmentations, TrainVoxelSize}];
	{inLab, outLab} = {inLab, outLab} /. Automatic -> {0};

	(*look for the files in the given folder*)
	segFiles = FileNames["*" <> labT <> ".nii.gz", labFol];
	datFiles = FileNames["*" <> datT <> ".nii.gz", datFol];

	(*prepare stuff for monitoring*)
	i = 1; 
	out = {i, "", ""};
	im = Image[{{0}}];
	
	PrintTemporary["Number of segmentation files: ", Length@segFiles];
	PrintTemporary[Dynamic[out]];
	If[!test, PrintTemporary[Dynamic[Show[im, ImageSize -> 300]]]];

	(*loop over segmentation files check for data and validate*)
	out = Table[
		name = StringTrim[StringReplace[FileBaseName@FileBaseName[sf], {labT -> ""}], "_" ...];
		df = Select[datFiles, StringContainsQ[#, StringReplace[Last@FileNameSplit[sf], labT -> datT]] &];

		(*check if data file exist*)
		If[df === {},
			out = {i++, name, "Data file does not exist"},

			(*import data and label*)
			{seg, vox} = ImportNii@sf;
			{dat, voxd} = ImportNii@First@df;
			If[voxOut === Automatic, voxOut = vox];

			dimSeg = Dimensions[seg]; 
			dimDat = If[ArrayDepth[dat] === 4, Dimensions[dat[[All,1]]], Dimensions[dat]];

			(*check if data file exist*)

			(*check dimensions and voxel size*)
			If[vox =!= voxd,
				out = {i++, name, "Data and segmentation have different voxel size."},
				If[dimSeg =!= dimDat,
					out = {i++, name, "Data and segmentation have different dimensions size."},

					(*Prepare and analyze the training data and segmentation*)
					{dat, seg} = PrepTrainData[{dat, seg}, {inLab, outLab}, {vox, voxOut}];

					(*output label check*)
					err = CheckSegmentation[seg];
					out = {i++, name, err};

					(*Cleanup if needed*)
					If[clean && (!test), 
						seg = SmoothSegmentation[seg, MaskComponents -> 1, MaskClosing -> True, SmoothIterations -> 1]
					];

					(*export*)
					If[!test,
						im = MakeChannelClassGrid[{dat}, seg, 5];
						outf = FileNameJoin[{outFol, name}];
						ExportNii[dat, voxOut, outf <> "_data.nii"];
						ExportNii[seg, voxOut, outf <> "_label.nii"];
						Export[outf <> ".png", im, "ColorMapLength" -> 256];
						Export[outf <> ".wxf", {dat, seg, voxOut}, PerformanceGoal -> "Size", Method -> {"PackedArrayRealType" -> "Real32"}];
					];

					out
				]
			]
		]
	, {sf, segFiles}];

	(*export the overview of what has happened*)
	legend = Grid[{{}, Join[{""}, Item[Style[#[[1]], White, Bold], Background -> #[[2]]] & /@ {{"hole & n > 1", Red}, {"n > 1", Purple}, {"hole", Blue}}, {""}], {}}, Spacings -> {1, 0.5}];
	head = Style[#, Bold] & /@ {"#", "Name", "Labels"};

	out = Grid[Append[Prepend[out, head], {"", legend, SpanFromLeft}], 
		Spacings -> {1, 1}, Background -> {None, {{White, Lighter@LightGray}}}, Alignment -> Left];
	Export[FileNameJoin[{outFol, "summary.png"}], ImagePad[Rasterize[out], 6, White]];

	out
]


(* ::Subsubsection::Closed:: *)
(*PrepTrainData*)


SyntaxInformation[PrepTrainData] = {"ArgumentsPattern" -> {_, _, _.}};

PrepTrainData[{dat_?ArrayQ, seg_?ArrayQ}] := PrepTrainData[{dat, seg}, {{0}, {0}}, {{1,1,1}, {1,1,1}}]

PrepTrainData[{dat_?ArrayQ, seg_?ArrayQ}, {labi_?VectorQ, labo_?VectorQ}] := PrepTrainData[{dat, seg}, {labi, labo}, {{1,1,1}, {1,1,1}}]

PrepTrainData[{daI_?ArrayQ, segI_?ArrayQ}, {labi_?VectorQ, labo_?VectorQ}, {voxi_?VectorQ, voxo_?VectorQ}] := Block[{
		cr, dat, seg
	},
	(*rescale if needed*)
	{dat, seg} = If[voxi===voxo, 
		{daI, segI}, 
		{RescaleData[daI, {voxi, voxo}], RescaleSegmentation[segI, {voxi, voxo}]}
	];

	(*remove background and normalize data and figure out what to do with multi channel data*)
	cr = FindCrop[Mask[If[ArrayDepth[dat] === 3, NoramlizeData, NormalizeMeanData][dat], 5, MaskDilation -> 1]];
	{
		ApplyCrop[dat, cr],
		If[labi === {0},
			ApplyCrop[seg, cr],
			ReplaceSegmentations[ApplyCrop[seg, cr], labi, labo]
		]
	}
]


(* ::Subsubsection::Closed:: *)
(*SplitSegmentations*)


SyntaxInformation[CheckSegmentation] = {"ArgumentsPattern" -> {_, _.}};

CheckSegmentation[seg_]:=CheckSegmentation[seg, "label"]

CheckSegmentation[seg_, out_?StringQ] := Block[{arrDep, segs, lab, err},
	arrDep = ArrayDepth[seg];
	If[arrDep === 3, {segs, lab} = SplitSegmentations[seg], lab = Range@Length@First@seg; segs = seg];
	segs = Transpose[segs];
	err = Thread[{CheckSegmentation[segs, 0], CheckSegmentation[segs, 1]}];
	lab = Grid[{Switch[#[[2]],
			{0, 0}, Style[#[[1]], Black, Bold],
			{1, 1}, Item[Style[#[[1]], White, Bold], Background -> Red],
			{0, 1}, Item[Style[#[[1]], White, Bold], Background -> Blue],
			{1, 0}, Item[Style[#[[1]], White, Bold], Background -> Purple]
		] & /@ Thread[{lab, err}]
	}];

	If[out==="label", lab, err]
]

CheckSegmentation[seg_, i_?IntegerQ] := Unitize[Max[MorphologicalComponents[Image3D[NumericArray[Abs[i - First@AutoCropData@#], "Integer8"]], CornerNeighbors -> False]] - 1 & /@ seg]


(* ::Subsection:: *)
(*Make evaluation images*)


(* ::Subsubsection::Closed:: *)
(*MakeChannelClassGrid*)


SyntaxInformation[MakeChannelClassGrid] = {"ArgumentsPattern"->{_, _, _.}};

MakeChannelClassGrid[dat_, lab_] := MakeChannelClassGrid[dat, lab, 3]

MakeChannelClassGrid[dat_, lab_, ni_] := Block[{len, n1, n2, labp, ran},
	len = Length@First@dat;
	If[IntegerQ[ni],
		n1 = n2 = Min[{Floor[Sqrt[len]], ni}],
		{n1, n2} = ni;
		While[n1 n2 > l, n1--; n2--;]
	];

	{labp, ran} = If[TensorQ[lab], {lab, MinMax@lab}, lab];
	If[IntegerQ[Round[ran]], ran = {0, ran}];

	RemoveAlphaChannel@ImageAssemble@Partition[
		ImagePad[MakeChannelClassImage[dat[[{1}, #]], labp[[#]], ran], 4, White
	] & /@ (Round[Range[1., len, (len - 1)/(n1 n2 - 1)]]), n1]
]


(* ::Subsubsection::Closed:: *)
(*MakeChannelClassImage*)


SyntaxInformation[MakeChannelClassImage]={"ArgumentsPattern"->{_, _, _., _.}};

MakeChannelClassImage[data_, label_]:=MakeChannelClassImage[data, label, MinMax[label], {1,1,1}]

MakeChannelClassImage[data_, label_, {off_, max_}]:=MakeChannelClassImage[data, label, {off,max}, {1,1,1}]

MakeChannelClassImage[data_, label_, vox_]:=MakeChannelClassImage[data, label, MinMax[label], vox]

MakeChannelClassImage[data_, label_, {off_, max_}, vox_]:=Block[{i1, i2},
	i1 = MakeClassImage[label, {off, max}, vox];
	i2 = MakeChannelImage[data, vox];
	ImageCollage[ImageCompose[#, SetAlphaChannel[i1, 0.4 AlphaChannel[i1]]]& /@ i2]
]


(* ::Subsubsection::Closed:: *)
(*MakeClassImage*)


SyntaxInformation[MakeClassImage]={"ArgumentsPattern"->{_, _., _.}};

MakeClassImage[label_]:=MakeClassImage[label, Round@MinMax[label], {1,1,1}]

MakeClassImage[label_, {off_?NumberQ, max_?NumberQ}]:=MakeClassImage[label, {off, max}, {1,1,1}]

MakeClassImage[label_, vox_?VectorQ]:=MakeClassImage[label, Round@MinMax[label], vox]

MakeClassImage[labelI_,{offI_?NumberQ, maxI_?NumberQ}, vox_?VectorQ]:=Block[{max, cols, imlab, rat, label, off},
	(*SeedRandom[1345];
		cols = Prepend[ColorData["DarkRainbow"][#]&/@RandomSample[Rescale[Range[off+1, max]]],Transparent];
		cols = Prepend[ColorData["RomaO"][#]&/@Rescale[Range[off+1, max]],Transparent];
	*)
	{label, off, max} = Round[{labelI, offI, maxI}];
	max = Max[{Max[label], max}];
	cols = Prepend[ColorData["RomaO"][#]&/@Rescale[
			Join[Select[Range[off + 1, max], EvenQ], Select[Range[off + 1, max], OddQ]]
		],Transparent];
	imlab = Round@Clip[If[ArrayDepth[label] === 3, label[[Round[Length@label/2]]], label] - off + 1, {1, max + 1}, {1, 1}];
	rat = vox[[{2,3}]]/Min[vox[[{2,3}]]];
	ImageResize[Image[cols[[#]]&/@imlab], Round@Reverse[rat Dimensions[imlab]], Resampling->"Nearest"]
]


(* ::Subsubsection::Closed:: *)
(*MakeChannelImage*)


SyntaxInformation[MakeChannelImage]={"ArgumentsPattern"->{_, _., _.}};

MakeChannelImage[data_]:=MakeChannelImage[data, {1, 1, 1}]

MakeChannelImage[data_, vox_]:=Block[{dat, imdat, rat},
	dat=Clip[Rescale[data, Quantile[Flatten[data], {0.01, 0.99}]], {0., 1.}];
	(*dat = Rescale[data];*)
	rat = vox[[{2, 3}]] / Min[vox[[{2, 3}]]];
	(
		imdat = #;
		imdat = If[ArrayDepth[#]===3, imdat[[Round[Length@imdat/2]]], imdat];
		ImageResize[Image[imdat], Round@Reverse[rat Dimensions[imdat]], Resampling->"Nearest"]
	) &/@ dat
]


(* ::Subsection:: *)
(*Muscle Names*)


(* ::Subsubsection::Closed:: *)
(*ImportITKLabels*)


SyntaxInformation[ImportITKLabels] = {"ArgumentsPattern"->{_.}};

ImportITKLabels[] := ImportITKLabels[GetAssetLocation["MusclesLegLabels"]];

ImportITKLabels[file_] := Block[{fileL, lines, muscleNames, muscleLabels},
	fileL = If[FileExistsQ[file], file, GetAssetLocation[file]]; 
	If[fileL === $Failed, Return["specified name is not file or asset"]];
	(*import*)
	lines = Select[Import[fileL, "Lines"], StringTake[#, 1] =!= "#" &];
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


SyntaxInformation[MuscleLabelToName] = {"ArgumentsPattern"->{_, _.}};

MuscleLabelToName[num_] := num /. Thread[#[[2]] -> #[[1]]] &[ImportITKLabels[GetAssetLocation["MusclesLegLabels"]]]

MuscleLabelToName[num_, file_] := Block[{muscleNames, muscleLabels},
	{muscleNames, muscleLabels} = ImportITKLabels[file];
	num /. Thread[muscleLabels -> muscleNames]
]


(* ::Subsubsection::Closed:: *)
(*MuscleLabelToName*)


SyntaxInformation[MuscleNameToLabel] = {"ArgumentsPattern"->{_, _.}};

MuscleNameToLabel[num_] := num /. Thread[#[[1]] -> #[[2]]] &[ImportITKLabels[GetAssetLocation["MusclesLegLabels"]]]

MuscleNameToLabel[num_, file_] := Block[{muscleNames, muscleLabels},
	{muscleNames, muscleLabels} = ImportITKLabels[file];
	num /. Thread[muscleNames -> muscleLabels]
]


(* ::Subsection:: *)
(*Distance measures*)


(* ::Subsubsection::Closed:: *)
(*DiceSimilarity*)


SyntaxInformation[DiceSimilarity] = {"ArgumentsPattern" -> {_, _, _.}};

DiceSimilarity[ref_, pred_, nClasses_?ListQ] := Block[{refF, predF},
	refF = flatRound@ref;
	predF = flatRound@pred;
	Table[DiceSimilarityC[refF, predF, class], {class, nClasses}]
]

DiceSimilarity[ref_, pred_] := DiceSimilarityC[flatRound@ref, flatRound@pred, 1]

DiceSimilarity[ref_, pred_, class_?IntegerQ] := DiceSimilarityC[flatRound@ref, flatRound@pred, class]


DiceSimilarityC = Compile[{{ref, _Integer, 1}, {pred, _Integer, 1}, {class, _Integer, 0}}, Block[{refv, predv, inter},
	refv = 1 - Unitize[ref - class];
	predv = 1 - Unitize[pred - class];
	inter = Total[refv predv];
	N[(2 inter + 1) / (Total[refv] + Total[predv] + 1)]]
, RuntimeOptions -> "Speed", Parallelization -> True];


flatRound=Flatten@Round@ToPackedArray@#&


(* ::Subsubsection::Closed:: *)
(*JaccardSimilarity*)


SyntaxInformation[JaccardSimilarity] = {"ArgumentsPattern" -> {_, _, _.}};

JaccardSimilarity[ref_, pred_, nClasses_?ListQ] := Block[{refF, predF},
	refF = flatRound@ref;
	predF = flatRound@pred;
	Table[JaccardSimilarityC[refF, predF, class], {class, nClasses}]
]

JaccardSimilarity[ref_, pred_] := JaccardSimilarityC[flatRound@ref, flatRound@pred, 1]

JaccardSimilarity[ref_, pred_, class_?IntegerQ] := JaccardSimilarityC[flatRound@ref, flatRound@pred, class]


JaccardSimilarityC = Compile[{{ref, _Integer, 1}, {pred, _Integer, 1}, {class, _Integer, 0}}, Block[{refv, predv, inter},
	refv = 1 - Unitize[ref - class];
	predv = 1 - Unitize[pred - class];
	inter = Total[refv predv];
	N[(inter + 1) / (Total[refv] + Total[predv] - inter + 1)]]
, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SurfaceDistance*)


Options[SurfaceDistance] = {
	Method->"HD95"
};

SyntaxInformation[SurfaceDistance] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

SurfaceDistance[ref_, pred_, opts : OptionsPattern[]] := SurfaceDistance[Round@ToPackedArray@ref, Round@ToPackedArray@pred, 1, {1, 1, 1}, opts]

SurfaceDistance[ref_, pred_, class_?IntegerQ, opts : OptionsPattern[]] := SurfaceDistance[Round@ToPackedArray@ref, Round@ToPackedArray@pred, class, {1, 1, 1}, opts]

SurfaceDistance[ref_, pred_, nClasses_?ListQ, opts : OptionsPattern[]] := SurfaceDistance[ref, pred, nClasses, {1, 1, 1}, opts]

SurfaceDistance[ref_, pred_, nClasses_?ListQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, opts : OptionsPattern[]] := Block[{refF, predF},
	refF = Round@ToPackedArray@ref;
	predF = Round@ToPackedArray@pred;
	Table[SurfaceDistance[refF, predF, class, vox, opts], {class, nClasses}]
]

SurfaceDistance[ref_, pred_, class_?IntegerQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, opts :OptionsPattern[]] := Block[{
		coorRef, coorPred, funRef, funPred, met, dist
	},

	coorRef = GetEdge[ref, class, vox];
	coorPred = GetEdge[pred, class, vox];
	If[coorRef==={}||coorPred==={},
		"noSeg",
		met = OptionValue[Method];

		funRef = Nearest[coorRef, DistanceFunction -> EuclideanDistance];
		funPred = Nearest[coorPred, DistanceFunction -> EuclideanDistance];
		dist = Sqrt@Total[Join[
			funRef[coorPred, 1][[All,1]] - coorPred,
			funPred[coorRef, 1][[All,1]] - coorRef
		]^2, {2}];

		If[ListQ[met],
			SufDistFunc[dist, #]&/@met,
			SufDistFunc[dist, met]
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*SufDistFunc*)


SufDistFunc[dist_, met_] := Switch[met,
	"Mean", Mean@dist,
	"Median", Median@dist,
	"RootMeanSquare"|"RMS", Sqrt[Mean[dist^2]],
	"Max"|"Hausdorff"|"HD", Max@dist,
	"Hausdorff95"|"HD95", Quantile[dist,.95],
	"Std"|"StandardDeviation", StandardDeviation@dist,
	_ , Message[SurfaceDistance::met, met]; $Failed
]


(* ::Subsubsection::Closed:: *)
(*GetEdge*)


GetEdge[lab_] := GetEdge[lab, -1, {1, 1, 1}]

GetEdge[lab_, vox_?VectorQ] := GetEdge[lab, -1, vox]

GetEdge[lab_, class_?IntegerQ] := GetEdge[lab, class, {1, 1, 1}]

GetEdge[lab_, class_?IntegerQ, vox_?VectorQ] := Block[{edge, im},
	im = Image3D[If[class === -1, lab, 1 - Unitize[lab - class]], "Bit"];
	edge = SparseArray[ImageData[MorphologicalPerimeter[im, CornerNeighbors -> False, Padding -> 0], "Bit"]]["ExplicitPositions"];
	If[edge =!= {}, Transpose[vox Transpose[edge]], edge]
]


(* ::Subsection::Closed:: *)
(*MakeDistanceMap*)


Options[MakeDistanceMap] = {
	DistanceRange -> Automatic
};

SyntaxInformation[MakeDistanceMap] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

MakeDistanceMap[data_, opt:OptionsPattern[]] := MakeDistanceMap[data, {1., 1., 1.}, opt]

MakeDistanceMap[data_, vox_, OptionsPattern[]] := Block[{
		dat, dim, datDil, datc, dimC, edge, inner, outer, nearFun, din, dout, cr, outRange
	},
	outRange = OptionValue[DistanceRange];

	dat = Unitize@Round@Normal@data;
	dim = Dimensions[dat];
	If[outRange === Automatic, outRange = Round[Max[1 /. DeleteCases[ComponentMeasurements[#, "Length"] & /@ dat, {}]]/3]];

	datDil = Switch[outRange,
		0, SparseArray@dat,
		All, SparseArray[1 - dat],
		_, SparseArray@Dilation[dat, outRange]
	];
	{datDil, cr} = AutoCropData[datDil, CropPadding -> 0];

	datc = ApplyCrop[SparseArray@dat, cr];
	dimC = Dimensions[datc];

	edge = GetEdge[datc, vox];
	inner = datc["ExplicitPositions"];
	outer = (datDil - datc)["ExplicitPositions"];

	nearFun = Nearest[edge];
	din = DistFun[nearFun, vox # & /@ inner];
	If[outer === {}, {}, dout = -DistFun[nearFun, vox # & /@ outer]];

	ReverseCrop[SparseArray[Join[Thread[inner -> din], Thread[outer -> dout]], dimC,0.], dim, cr]
]

DistFun[fun_, pts_] := Sqrt[Total[(Flatten[fun[#, 1] & /@ pts, 1] - pts)^2, {2}]];


(* ::Subsection::Closed:: *)
(*NetSummary*)


NetSummary[net_] := NetSummary[net, ""]

NetSummary[net_?ListQ, rep_?StringQ]:=NetSummary[#, rep]&/@net

NetSummary[net_, rep_?StringQ] := Block[{
		toK, st, quantStr, lays, convs, kerns, count, nKern, kernWeights,
		norm, normWeights, nelem, elems, elemSize, arrSize, netSize, table, nodes, netIm, nodeIm
	},

	toK = Which[
		# > 1000000, ToString[NumberForm[#/1000000., {Infinity, 2}]] <> " M",
		# > 1000, ToString[NumberForm[#/1000., {Infinity, 2}]] <> " K",
		True, #] &;
	st = Style[#1, Bold] & ;
	quantStr = ToString[Round[QuantityMagnitude[#], .01]] <> " " <> (QuantityUnit[#] /. {"Megabytes" -> "MB", "Gigabytes" -> "GB"}) &;

	lays = Information[net, "LayersList"];
	convs = Select[lays, Head[#] === ConvolutionLayer &];
	kerns = Information[#, "ArraysDimensions"][{"Weights"}] & /@ convs;
	kerns = GatherBy[{#[[3 ;;]], Times @@ #[[1 ;; 2]], Times @@ #[[1 ;;]]} & /@	kerns, First];

	count = Sort[{Length[#], #[[1, 1]], Total[#[[All, 2]]], Total[#[[All, 3]]]} & /@ kerns];
	nKern = Total[count[[All, 3]]];
	kernWeights = Total[count[[All, 4]]];
	count[[All, 3]] = toK /@ count[[All, 3]];
	count[[All, 4]] = toK /@ count[[All, 4]];

	norm = Select[lays, Head[#] === BatchNormalizationLayer &];
	normWeights = First@Total[Total[Information[#, "ArraysDimensions"] /@ Keys[Information[norm[[1]], "ArraysDimensions"]]] & /@ norm];

	nelem = Information[net, "ArraysTotalElementCount"];
	elems = Round[nelem/1000000, .01];
	elemSize = UnitConvert[Quantity[32. nelem, "Bits"], "MB"];
	arrSize = UnitConvert[Quantity[32. Total[Times @@@ Values[(Information[#, "OutputPorts"]["Output"]) & /@ Information[net, "Layers"]]], "Bits"], "MB"];
	netSize = UnitConvert[elemSize + arrSize, "GB"];

	table = Grid[{
		{st@"Number of batch norm. Layers: ", st@Length@norm},
		{st@" - Number of Weights: ", toK@normWeights},
		{st@"Number of convolution Layers: ", st@Length@convs},
		{st@" - Number of Kernels: ", toK@nKern},
		{st@" - Number of Weighths: ", toK@kernWeights},
		{""},
		{st@"Convolution Kernel Distribution:", SpanFromLeft},
		{Item[Grid[Join[{Style[#, Bold] & /@ {"Count", "Size", "Kernels", "Weights"}}, count], 
			Alignment -> Right, Spacings -> {1.5, 1}], Alignment->Right], SpanFromLeft},
		{""},
		{st@"Total Weight Memory", quantStr@elemSize},
		{st@"Total Network Memory", quantStr@netSize}
	}, Alignment -> {{Left, Right}}, Spacings -> {1, 1}, Background -> GrayLevel[.95]];

	Switch[rep,
		"Full",
		makeNetIm = With[{im = Information[#, "SummaryGraphic"]}, Show[im, AspectRatio -> 0.5, ImageSize -> Max[AbsoluteOptions[im, ImageSize][[1, 2]]]]] &;
		nodes = DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]];
		netIm = makeNetIm @ net;
		nodeIm = Information[NetTake[net, {#, #}][[1]], "SummaryGraphic"] & /@ nodes;
		TabView[Join[{"summary" -> table, "net" -> netIm}, Thread[nodes -> nodeIm]],
			ControlPlacement -> Left, Alignment -> Center],

		"Mem",
		Grid[{{st["Weight mem: "], quantStr[elemSize]}, {st["Net mem: "], quantStr[netSize]}}, Alignment -> Left],

		_, 
		table
	]
]


(* ::Subsection::Closed:: *)
(*AnalyzeNetworkFeatures*)


AnalyzeNetworkFeatures[net_, data_] := AnalyzeNetworkFeatures[net, data, ""]

AnalyzeNetworkFeatures[net_, data_, met_] := Block[{
		dim, dataP, netP, nodes, vals, cutoff, table, plot, feat, nfeat, ttt, n, col
	},

	(*find the patch dimensions and adjust data and network*)
	dim = FindPatchDim[net, Dimensions@data];
	dataP = NormalizeData[PadToDimensions[data, dim], NormalizeMethod -> "Uniform"];
	netP = ChangeNetDimensions[net, "Dimensions" -> dim];

	(*extract the network nodes*)
	nodes = GetNetworkNodes[netP];
	nodes = nodes[[2;;-2]];

	col = ColorData["SolarColors"] /@Rescale[Range[Ceiling[Length[nodes]/2]]];
	col = Join[col , Reverse[col][[2;;]]];

	(*calculate the singular values for plotting and reporting*)
	{nfeat, table, plot} = Transpose[(
		(*get the features*)
		feat = NetTake[netP, #][{dataP}, TargetDevice -> "GPU", WorkingPrecision -> "Mixed"];
		If[Head[feat] === Association, feat = Last@feat];
		(*If[# == nodes[[-1]], feat = RotateDimensionsRight[feat]];*)(*only for map*)
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
		Echo[Thread[{nodes, nfeat}]];
		Echo[n];

		(*dynamic plot output*)
		DynamicModule[{cols, tab = table, pl = plot, nods =  nodes, clist = col, ln, pcol},

			(*define colors for plotting*)
			cols = Table[Directive[{(*GrayLevel[.5 + i/100]*)clist[[i]], Dashed, Thick}], {i, Length@nodes}];
			(*cols[[1]] = Directive[{Red, Dashing[None], Thick}];*)

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
						ImageSize -> 500, AxesStyle -> Directive[{Black, Thick}], AspectRatio -> 1, 
						LabelStyle -> Directive[{Black, 14, Bold}]
					],
					Plot[x, {x, 0, 100}, PlotStyle -> Directive[{Thick, Gray, Dotted}]]
				]
			}, Alignment -> Center],
			{{k, 1, ""}, Thread[ln -> (Style[#, Black, 14, Bold] & /@ nods)], ControlType -> SetterBar}
			]
		],

		(*value per node*)
		table[[All, 2]]
	]
]


(* ::Subsection:: *)
(*ShowTrainLog*)


(* ::Subsubsection::Closed:: *)
(*ShowTrainLog*)


ShowTrainLog[fol_] := ShowTrainLog[fol, 5]

ShowTrainLog[fol_, max_] := Block[{files, log, keys, leng, plots},
	{keys, log, leng} = LoadLog[fol, max];
	log = log[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];

	(* Create a dynamic module to display the interactive plot *)
	DynamicModule[{pdat = log, klist = keys, folder = fol, len = leng, plot, ymaxv, xmin, xmax, key, ymin, ymax, temp, plots},
		key1 = Select[klist, ! StringContainsQ[#, "Current"] &];
		key2 = Select[Select[key1, StringContainsQ[#, "Loss"] &], # =!= "RoundLoss" && # =!= "ValidationLoss" &];
		key1 = Select[Complement[key1, key2], # =!= "RoundLoss" && # =!= "ValidationLoss" &];
		key0 = {"RoundLoss", "ValidationLoss"};

		Manipulate[
			plot = Transpose[Values /@ Normal[pdat[All, key]][[All, All]]];
			plotf = If[filt, GaussianFilter[#, fsize]&/@plot, plot];

			ymaxv = Max[{1.1, 1.1 If[plot==={}, 1, Max[Select[Flatten@plot,NumberQ]]]}];
			ymax = Min[{ymax, ymaxv}];

			(* Plot the selected metrics *)
			If[logp, ListLogPlot, ListLinePlot][If[key === {}, {}, plotf], Joined -> True, 
				PlotLegends -> Placed[key, Right], ImageSize -> 600, PlotRange->{{xmin,xmax} ,{ymin,ymax}},
				If[grid, GridLines -> {len, Automatic}, GridLines -> {len, None}], PlotHighlighting -> "YSlice"],

			(*the controls*)
			Row[{
				InputField[Dynamic[folder], String, Enabled -> True, FieldSize -> 50], 
				Button["Browse", 
					temp = SystemDialogInput["Directory", folder];
					If[StringQ[temp], folder = temp; 
						{klist, pdat, len} = LoadLog[folder, max];
						pdat = pdat[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];
						xmax = Length[pdat];
					];
					, ImageSize -> {60, Automatic}, Method->"Queued"]}
			],
			Button["Reload", 
				{klist, pdat, len} = LoadLog[folder, max]; xmax = Length[pdat];
				pdat = pdat[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];
			],

			Delimiter,
			{{filt, False, "Filter"}, {True, False}},
			{{fsize, 5, "FilterSize"}, 1, 10, 1},
			{{grid, False, "Grid"}, {True, False}},
			{{logp, False, "Log"}, {True, False}},

			Delimiter,
			(*Control[{{key, {}, ""}, klist, ControlType -> TogglerBar, Appearance -> "Vertical" -> {Automatic, 4}, BaseStyle -> Medium}],*)
			Control[{{key, {}, ""}, key0, ControlType -> TogglerBar, 
			Appearance -> "Vertical" -> {Automatic, 4}, BaseStyle -> Medium}],
			Delimiter,
			Control[{{key, {}, ""}, key2, ControlType -> TogglerBar, 
			Appearance -> "Vertical" -> {Automatic, 4}, BaseStyle -> Medium}],
			Delimiter,
			Control[{{key, {}, ""}, key1, ControlType -> TogglerBar, 
			Appearance -> "Vertical" -> {2, Automatic}, BaseStyle -> Medium}],
			Delimiter,

			Row[{
				Control[{{xmin, 1,"X min"},1, Dynamic[xmax-1], 1}], "   ", 
				Control[{{xmax,Length[pdat],"X max"}, Dynamic[xmin+1], Dynamic[Length[pdat]], 1}]
			}],
			Row[{
				Control[{{ymin, 0, "Y min"}, 0, Dynamic[ymax-0.01]}], "   ",
				Control[{{ymax, 1, "Y max"}, Dynamic[ymin+0.01], Dynamic[ymaxv]}]
			}],
			Row[{
				Button["Autoscale X", {xmax, xmax} = {1, Length[pdat]}], "   ",
				Button["Autoscale Y", {ymin, ymax} = {0, Max[{1.1, 1.1 If[plot==={}, 1, Max[Select[Flatten@plot,NumberQ]]]}]}]
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


(* ::Subsubsection::Closed:: *)
(*LoadLog*)


LoadLog[fol_, max_]:=Block[{files, keys, log, leng},
	(* Get a list of log files in the specified folder *)
	files = Sort[FileNames["*.json", fol]];

	(* Read the log files and extract the relevant information *)
	log = Select[(Select[Import[#, "Lines"], StringContainsQ[#, "ProgressFraction"] &] & /@ files), Length[#] > max &];
	leng = Accumulate[Length /@ log];

	(* Convert the log data into a dataset *)
	log = "[\n" <> StringDrop[StringRiffle[If[StringTake[#, -1] === "}", # <> ",", #] & /@ Flatten[log, 1], "\n"], -1] <> "\n]";
	log = Dataset[Association /@ Import[Export[FileNameJoin[{$TemporaryDirectory, "log.json"}], log, "text"]]];

	(* Get the unique keys (metrics) in the log data *)
	keys = Sort@DeleteDuplicates[Flatten[Normal@log[All, Keys]]];

	{keys, log, leng}
]


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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
