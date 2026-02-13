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
Possible loss functions are {\"SoftDice\", \"MSD\", \"Tversky\" , \"CE\", \"Jaccard\"}."

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


(* ::Subsection::Closed:: *)
(*Options*)


LoadTrainingData::usage =
"LoadTrainingData is an option for TrainSegmentationNetwork. If set to True the training data is loaded from the disk."

MonitorInterval::usage =
"MonitorInterval is an option for TrainSegmentationNetwork. It defines how often the training is monitored."

L2Regularization::usage =
"L2Regularization is an option for TrainSegmentationNetwork. It defines the L2 regularization factor."

MultiChannel::usage = 
"MultiChannel is an option for TrainSegmentationNetwork, If set to True it will train on multi channel input data. If set to False it will select a random channel."


PatchesPerSet::usage =
"PatchesPerSet is an option for GetTrainData. Defines how many random patches per dataset are created within the batch."

AugmentData::usage = 
"AugmentData is an option for GetTrainData and TrainSegmentationNetwork. If set True the training data is augmented.
It can also be set to \"2D\" or \"3D\" to control if augmentation is done trought plane or only inplane."

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


SurfaceDistance::met = "Method `1` not recognized";


ApplySegmentationNetwork::node = "The node ``` is not part of the network"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


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

MuscleNameToLabel[name_] := name /. Thread[#[[1]] -> #[[2]]] &[ImportITKLabels[GetAssetLocation["MusclesLegLabels"]]]

MuscleNameToLabel[name_, file_] := Block[{muscleNames, muscleLabels},
	{muscleNames, muscleLabels} = ImportITKLabels[file];
	name /. Thread[muscleNames -> muscleLabels]
]


(* ::Subsection:: *)
(*ClassifyData*)


(* ::Subsubsection::Closed:: *)
(*ClassifyData*)


Options[ClassifyData] = {
	TargetDevice -> "CPU"
};

SyntaxInformation[ClassifyData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ClassifyData[dat_, met_, OptionsPattern[]] := Block[{
		dev, net, imSize, ims
	},

	dev = OptionValue[TargetDevice];

	(*get the network*)
	net = Which[
		FileExistsQ[met] && FileExtension[met]==="wlnet", Import[met][data],
		StringQ[met], net = GetNeuralNet[met],
		Head[met]===NetChain || Head[met]===NetGraph, met,
		True, $Failed
	];
	If[net === $Failed, Return[$Failed]];

	(*convert data *)
	imSize = NetDimensions[NetReplacePart[net, "Input"->None], "Input"][[2;;]];
	ims = MakeClassifyImage[dat, ImageSize -> imSize];

	Switch[met,
		"LegSide"|"ShoulderSide", FindSideClass[ims, net, dev],
		"LegPosition", FindLegPos[ims, net, dev],
		"ShoulderPosition", FindShoulderPos[dat, ims, net, dev]
	]
]


(* ::Subsubsection::Closed:: *)
(*FindSideClass*)


FindSideClass[data_, net_, dev_] := Last@Keys@Sort@Counts@net[data, TargetDevice -> dev]


(* ::Subsubsection::Closed:: *)
(*FindSidePos*)


FindSidePos[data_, net_, dev_] := First@Mean@net[data, TargetDevice -> dev]


(* ::Subsubsection::Closed:: *)
(*FindShoulderPos*)


FindShoulderPos[data_, ims_, net_, dev_] := Block[{pos},
	pos = FindSidePos[ims, net, dev];
	Round[Dimensions[data][[-2]]{pos, 0.1}]
]


(* ::Subsubsection::Closed:: *)
(*FindLegPos*)


FindLegPos[data_, net_, dev_] := Block[{
		len, ran, datF, kneeStart, kneeEnd, pos, imSize
	},

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


PosFunc = Compile[{{a, _Integer, 0}, {b, _Integer, 0}, {l, _Integer, 0}, {x, _Integer, 1}, {d, _Real, 1}},
	Total[((Which[1<=#<=a, 1., a<=#<=b, 2., b<=#<=l, 3., True, 0] &/@ x) - d)^2]
];


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

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, opts:OptionsPattern[]] := DataToPatches[dat, patch, "All", opts] 

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}] := {GetPatch[dat, patch, pts], pts}

DataToPatches[dat_, pts:{{{_,_},{_,_},{_,_}}..}] := {GetPatch[dat, pts], pts}

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, nPatch_, OptionsPattern[]] := Block[{
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


GetPatch[dat_, pts:{{{_,_},{_,_},{_,_}}..}] := GetPatch[dat, #]&/@pts

GetPatch[dat_, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}] := ToPackedArray@dat[[i1;;i2,j1;;j2,k1;;k2]]

GetPatch[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}] := GetPatch[dat, patch, #]&/@pts

GetPatch[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, {{i1_,i2_}, {j1_,j2_}, {k1_,k2_}}] := ToPackedArray@PadRight[dat[[i1;;i2,j1;;j2,k1;;k2]], patch, 0.]


(* ::Subsubsection::Closed:: *)
(*GetPatchRanges*)


GetPatchRanges[dat_, patch_, nPatch_, {nRan_, pad_}] := Block[{pts},
	pts = GetPatchRangeI[Dimensions[dat], patch, {nRan, pad}];
	If[nPatch==="All", pts, RandomSample[pts, Min[{Length@pts, nPatch}]]]
]


GetPatchRangeI[datDim_?ListQ, patchDim_?ListQ, {nr_, pad_}] := Tuples@MapThread[GetPatchRangeI[#1,#2, {nr, pad}]&, {datDim, patchDim}]

GetPatchRangeI[dim_?IntegerQ, patch_?IntegerQ, {nr_, pad_}] := Block[{i,st},
	i = Ceiling[(dim - 2 pad)/(patch - 2 pad)] + nr;
	If[!(dim > patch && i > 1),
		{{1,dim}},
		st = Round[Range[0, 1, 1./(i - 1)](dim - patch)];
		Thread[{st + 1, st + patch}]
	]
]


(* ::Subsection:: *)
(*SegmentData*)


(* ::Subsubsection::Closed:: *)
(*SegmentData*)


Options[SegmentData] = {
	TargetDevice -> "CPU", 
	MaxPatchSize->Automatic, 
	Monitor->False
};

SyntaxInformation[SegmentData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SegmentData[datI_, opts:OptionsPattern[]] := SegmentData[datI, "Legs", opts]

SegmentData[datI_, whati_, OptionsPattern[]] := Block[{
		dev, max, mon, patch, pts, dim ,loc, set, net, type, segs, all, data,
		time, timeAll, what, netFile, custom, monO
	},

	timeAll = First@AbsoluteTiming[
		{dev, max, mon} = OptionValue[{TargetDevice, MaxPatchSize, Monitor}];
		monO = mon;
		mon = If[mon, MonitorFunction, List];

		custom = If[ListQ[whati], 
			{what, netFile} = whati; True, 
			what = whati; False
		];

		(*split the data in upper and lower legs and left and right*)
		mon[Dimensions@datI, "Analyzing the data with dimensions:"];
		data = Switch[ArrayDepth[datI], 4, datI[[All,1]], _, datI];
		time = First@AbsoluteTiming[
			{{patch, pts, dim}, loc, set} = SplitDataForSegmentation[data, what, Monitor -> monO, TargetDevice -> dev]
		];
		mon[Round[time, .1], "Total time for analysis [s]: "];
		mon[set, "Settings: "];
		mon[Column@Thread[{loc, Dimensions/@ patch}], "Segmenting \""<>what<>"\" locations with dimensions:"];

		(*get the network name and data type*)
		{net, type} = Switch[what,
			"LegBones", {"SegLegBones"&, "Bones"},
			"Legs"|"UpperLegs"|"LowerLegs",	{
				(#[[1]] /. {"Upper" -> "SegThighMuscle", "Lower" -> "SegLegMuscle"})&, "Muscle"
			},
			"Shoulder", {"SegShoulderMuscle"&, "Muscle"},
			"Back", {"SegBackMuscle"&, "Muscle"},
			_, Return[]
		];
		If[custom, net = If[ListQ[netFile],
			(#[[1]] /. {"Upper" -> netFile[[1]], "Lower" -> netFile[[2]]})&,
			netFile&]
		];

		(*Perform the segmentation*)
		time = First@AbsoluteTiming[
			segs = MapThread[(
				mon[{#2, net[#2]}, "Performing segmentation for: "];
				segs = ApplySegmentationNetwork[#1, net[#2], TargetDevice -> dev, MaxPatchSize->max, Monitor->monO];
				If[custom, segs, ReplaceLabels[segs, #2, type]]
		) &, {patch, loc}]];
		mon[Round[time, .1], "Total time for segmentations [s]: "];

		(*Merge all segmentations for all expected labels*)
		all = Select[DeleteDuplicates[Sort[Flatten[GetSegmentationLabels/@segs]]], IntegerQ];
		mon[all, "Putting together the segmentations with labels"];

		(*after this only one cluster per label remains*)
		time = First@AbsoluteTiming[segs = PatchesToData[segs, pts, dim, all]];
		mon[Round[time, .1], "Total time for final evaluation [s]: "];
	];
	mon[Round[timeAll, .1], "Total evaluation time [s]: "];

	segs
]


(* ::Subsubsection::Closed:: *)
(*ReplaceLabels*)


ReplaceLabels[seg_, loc_, type_] := Block[{what, side, labNam, labNamS, labIn, labOut, labOutS,file, fileOut},
	{what, side} = loc;
	labIn = GetSegmentationLabels[seg];

	(*for now overwrite the labIn with custom values since some muscles are not segmented *)
	file = GetAssetLocation@Switch[type,
		"Muscle", Switch[what, 
			"Upper", "MusclesLegUpperLabels",
			"Lower", "MusclesLegLowerLabels",
			"Shoulder", "MusclesShoulderLabels"],
		"Bones", "BonesLegLabels"];

	labNam = MuscleLabelToName[labIn, file];
	labNamS = # <> "_" <> side & /@ labNam;
	
	fileOut = GetAssetLocation@Switch[type,
		"Muscle", Switch[what, 
			"Upper", "MusclesLegLabels",
			"Lower", "MusclesLegLabels",
			"Shoulder", "MusclesShoulderAllLabels"],
		"Bones", "MusclesLegLabels"];

	labOut = MuscleNameToLabel[labNam, fileOut];
	labOutS = MuscleNameToLabel[labNamS, fileOut];

	labOut = Cases[Transpose[{labOut, labOutS}], _Integer, 2];

	ReplaceSegmentations[seg, labIn, labOut]
]


(* ::Subsubsection::Closed:: *)
(*SplitDataForSegmentation*)


Options[SplitDataForSegmentation] = {
	Monitor -> False,
	TargetDevice -> "CPU"
};

SyntaxInformation[SplitDataForSegmentation] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};


SplitDataForSegmentation[data_?ArrayQ, seg_?ArrayQ, opt:OptionsPattern[]] := SplitDataForSegmentation[data, seg, "Legs", opt]

SplitDataForSegmentation[data_?ArrayQ, seg_?ArrayQ, what_?StringQ, opt:OptionsPattern[]] := Block[{dat,pts,dim,loc,set, segp},
	{{dat, pts, dim}, loc, set} = SplitDataForSegmentation[data, what, opt];
	segp = GetPatch[seg, pts];
	{{dat, pts, dim}, {segp, pts, dim}, loc, set}
]


SplitDataForSegmentation[data_?ArrayQ, opt:OptionsPattern[]] := SplitDataForSegmentation[data, "Legs", opt]

SplitDataForSegmentation[data_?ArrayQ, what_?StringQ, opt:OptionsPattern[]] := Block[{
		dim, whatSide, side, whatPos, pos, dat, right, left, cut, pts, loc, time, mon, dev, over
	},
	dim = Dimensions[data];
	{mon, dev} = OptionValue[{Monitor, TargetDevice}];
	mon = If[mon, MonitorFunction, List];

	Switch[what,
		"Legs"|"UpperLegs"|"LowerLegs",
		(*split the data in upper and lower legs and left and right*)

		(*find which side using NN*)
		time= First@AbsoluteTiming[whatSide = ClassifyData[data, "LegSide", TargetDevice -> dev]];
		mon[whatSide, "Data contains sides: "];
		mon[Round[time, .1], "Time for side estimation [s]:"];

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

		mon[whatPos, "Data contains positions: "];
		mon[Round[time, .1], "Time for position estimation [s]:"];

		(*output the selected data with the correct label and coordinates*)
		{dat, pts, loc} = Transpose[CropPart/@dat];
		{{dat, pts, dim}, loc, {{whatSide, cut}, {whatPos, pos}}}

		,
		"Shoulder",
		(*find which side using NN*)
		time= First@AbsoluteTiming[
			whatSide = ClassifyData[data, "ShoulderSide", TargetDevice -> dev];

			(*based on side cut data or propagate*)
			dat = Switch[whatSide,
				(*both sides which need to be split*)
				"Both",
				(*{cut, over} = ClassifyData[data, "ShoulderPosition"];*)
				{cut, over} = Round[{0.5, 0.1} Last@Dimensions[data]];
				{right, left, cut} = CutData[data, {cut, over}];
				{
					{right, {"Shoulder", {1, dim[[1]]}}, {"Right", {1, cut + over}}}, 
					{left, {"Shoulder", {1, dim[[1]]}}, {"Left", {cut - over + 1, dim[[3]]}}}
				},
				_,
				(*only one side, no split*)
				cut = over = 0;
				{
					{data, {"Shoulder", {1, dim[[1]]}}, {whatSide, {1, dim[[3]]}}}
				}
			];
		];

		mon[whatSide, "Data contains sides: "];
		mon[Round[time, .1], "Time for side estimation [s]:"];

		{dat, pts, loc} = Transpose[CropPart/@dat];
		{{dat, pts, dim}, loc, {{whatSide, {cut, over}}, {"Shoulder", 0}}}

		,
		"Back",
		dat = {{data, {"Both", {1, dim[[1]]}}, {"Both", {1, dim[[3]]}}}};
		{dat, pts, loc} = Transpose[CropPart/@dat];
		{{dat, pts, dim}, loc, {{whatSide, 0}, {whatPos, 0}}}

		,
		_,
		$Failed		
	]
]


(* ::Subsubsection::Closed:: *)
(*CropPart*)


CropPart[{dat_, {up_, {upStart_, upEnd_}}, {side_, {sideStart_, sideEnd_}}}] := Block[{data, crp},
	{data, crp} = AutoCropData[dat, CropPadding->0];
	{data, Partition[crp, 2] + {upStart-1, 0, sideStart-1}, {up, side}}
]


(* ::Subsubsection::Closed:: *)
(*ApplySegmentationNetwork*)


Options[ApplySegmentationNetwork] = {
	TargetDevice->"CPU", 
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
ApplySegmentationNetwork[dat_, netI_, opt:OptionsPattern[]] := ApplySegmentationNetwork[dat, netI, "", opt]

ApplySegmentationNetwork[dat_, netI_, node_, OptionsPattern[]] := Block[{
		dev, pad , lim, data, crp, net, dim, ptch, time, nodes,
		patch, pts, seg, mon, dimi, dimc, lab, nClass, prec, normF
	},

	{dev, pad, lim, mon} = OptionValue[{TargetDevice, DataPadding, MaxPatchSize, Monitor}];
	(*If[lim === Automatic, lim = If[dev==="GPU", 600, 150]];*)
	mon = If[mon, MonitorFunction, List];

	prec = If[(dev==="GPU") && ($OperatingSystem === "Windows"), "Mixed", "Real32"];

	data = Which[
		StringQ[dat], First@ImportNii[dat], 
		TensorQ[{dat}, NumericQ], If[First@Dimensions@dat===1, First@dat, dat],
		True, Return[$Failed]
	];
	If[ArrayDepth[data] === 4, data = data[[All, 1]]];

	dimi = Dimensions@data;
	{data, crp} = AutoCropData[data, CropPadding->0];
	dimc = Dimensions@data;
	data = N@ArrayPad[data, pad, 0.];
	dim = Dimensions[data];

	mon[{"in", dimi, "crop", dimc, "pad", dim}, "Data dimensions before and after cropping and padding are: "];

	net = If[StringQ[netI], GetNeuralNet[netI], netI];

	If[net===$Failed, $Failed,
		(*calculate the patch size for the data *)
		ptch = FindPatchDim[net, dim, lim];

		(*create the patches*)
		{patch, pts} = DataToPatches[data, ptch, PatchNumber -> 0, PatchPadding->pad];
		mon[{ptch, Length@patch}, "Patch size and created number of patches is:"];

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
			mon[{Dimensions[seg], Sort@Round@DeleteDuplicates[Flatten[seg]]}, "Output segmentations dimensions and labels:"];
			mon[Round[time, .1], "Time for segmentation [s]: "];

			,
			(*perform the segmentation on a specific node*)
			(*check if node is part of the network*)
			nodes = DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]];
			If[!MemberQ[nodes, node],
				Message[ApplySegmentationNetwork::node, node];
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
	dimM = sc Ceiling[dim/sc];
	(*dimM = sc ({Floor[#[[1]]], Ceiling[#[[2]]], Ceiling[#[[3]]]} & [dim/sc]);*)

	(*if memory limit is given find patch that fits*)
	If[!(lim === 1000 || lim === Automatic),
		If[CubeRoot[N[Times @@ dimM]] < lim,
			(*N[Sqrt[Times@@Rest[dimM]]]<lim,*)
			dimN = dimM
			,
			u = 1;
			cont = True;

			While[cont, u++;
				dimN = u sc;
				dimN = Min /@ Transpose[{dimN, dimM}];
				cont = ((CubeRoot[N[Times @@ dimN]] < lim && dimN =!= dimM) && u < 20)
			]
		],
		dimN = dimM
	];

	(*if no memory limit is given use the smallest patch that fits the network*)

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
	NormalizationType -> "Batch",
	RescaleMethod -> "Conv",
	CatenateMethod -> "Cat",

	NetworkDepth -> 5,
	DownsampleSchedule -> 2,
	SettingSchedule -> Automatic,
	FeatureSchedule -> 32,

	MultiChannel -> False,

	AugmentData -> True,
	PadData-> False,

	LossFunction -> All,
	DropoutRate -> 0.2,
	LearningRate -> 0.001,
	L2Regularization -> 0.0001,

	MonitorCalc -> False,
	TargetDevice -> "GPU"
}


SyntaxInformation[TrainSegmentationNetwork] = {"ArgumentsPattern" -> {{_, _}, _., OptionsPattern[]}};

TrainSegmentationNetwork[{inFol_?StringQ, outFol_?StringQ}, opts : OptionsPattern[]] := TrainSegmentationNetwork[{inFol, outFol}, "Start", opts]

TrainSegmentationNetwork[{inFol_?StringQ, outFol_?StringQ}, netCont_, opts : OptionsPattern[]] := Block[{
		netOpts, batch, roundLength, rounds, data, depth, nChan, nClass, outName, ittString, multi,
		patch, augment, netIn, ittTrain, testData, testVox, testSeg, im, patches, pLen, is2D,
		monitorFunction, netMon, netOut, trained, l2reg, pad, batchFunction, n, it, ti, br,
		validation, files, loss, rep, learningRate, schedule, dims, tar
	},

	(*------------ Get all the configuration stuff -----------------*)

	(*getting all the options*)
	netOpts = Join[FilterRules[{opts}, Options@MakeUnet], FilterRules[Options@TrainSegmentationNetwork, 
		Options@MakeUnet]];

	{batch, roundLength, rounds, augment, pad, patch, patches, 
		loss, rep, learningRate, l2reg, multi, tar} = OptionValue[
		{BatchSize, RoundLength, MaxTrainingRounds, AugmentData, PadData, PatchSize, PatchesPerSet, 
			LossFunction, MonitorInterval, LearningRate, L2Regularization, MultiChannel, TargetDevice}];
	pLen = Length@patch;
	is2D = pLen===2;
	pad = If[NumberQ[pad] && !is2D, Round[pad], False];

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

	MonitorFunction[DateString[], "Preparing the network"];

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
		loss === All, {"Dice", "MSD", "Tversky" , "CE", "Jaccard", "Focal"},
		StringQ[loss], {loss},
		True, loss];
	If[!And @@ (MemberQ[{"Dice", "MSD", "Tversky", "CE", "Jaccard", "Focal", "TopK"}, #] & /@ loss), 
		Return[Message[TrainSegmentationNetwork::loss]; $Failed]];

	(*if the network already exists make the dimensions, classes en channels match the input*)
	netIn = NetInitialize[ChangeNetDimensions[netIn, "Dimensions" -> patch, "Channels" -> nChan, "Classes" -> nClass],
		Method -> {"Kaiming", "Distribution" -> "Normal"}];
	MonitorFunction[NetSummary[netIn, "Mem"], "Network summary: "];

	(*define the network for training*)
	netIn = AddLossLayer@netIn;
	MonitorFunction[netIn, "Network is ready"];

	(*---------- Training functions ----------------*)

	(*Local functions*)
	outName = FileNameJoin[{outFol, Last[FileNameSplit[outFol]] <> "_" <> #}]&;
	ittString = "itt_" <> StringPadLeft[ToString[#], 4, "0"]&;

	(*Monitor function*)
	monitorFunction = (
		ittTrain++;
		(*perform segmentation and export*)
		netMon = NetExtract[#Net, "net"];
		testSeg = Ramp[ClassDecoder[netMon[testData, TargetDevice -> "CPU"]]];
		ExportNii[testSeg, testVox, outName[ittString[ittTrain]<>".nii"]];
		(*make and export test image*)
		im = MakeChannelClassGrid[If[is2D, Transpose@testData, testData], {testSeg, {0, nClass-1}}, 3];
		Export[outName[ittString[ittTrain] <> ".png"], im , "ColorMapLength" -> 256];
		(*export the network and delete the one from the last iteration*)
		Export[outName[ittString[ittTrain] <> ".wlnet"], netMon];
		Quiet@DeleteFile[outName[ittString[ittTrain - 1] <> ".wlnet"]];
	)&;

	(*batch function*)
	batchFunction = GetTrainData[data, #BatchSize, patch, nClass, 
		PatchesPerSet -> patches, AugmentData -> augment, PadData -> pad] &;

	(*oneCycle learning rate schedual function*)
	br = roundLength / batch;
	n = {0.25, 0.3, 0.80} rounds br;
	it = ittTrain br;
	schedule = (ti = #1 + it; Which[
		ti < n[[1]], Rescale[Cos[Pi ti / n[[1]]], {1, -1}, {1/5, 1}],
		ti < n[[2]], 1.,
		ti < n[[3]],Rescale[Cos[Pi (ti - n[[2]]) / (n[[3]] - n[[2]])], {1, -1}, {1/10, 1}],
		True, 1./10
	])&;

	(*---------- Prepare the data ----------------*)

	MonitorFunction[DateString[], "Preparing the data"];

	(*import all train data or train out of memory*)
	data = If[OptionValue[LoadTrainingData] === True, Import /@ files, files];
	dims = MeanRange[#, 0] & /@ Transpose[If[ArrayDepth[#] === 3, 
		Dimensions[Transpose[{#}]], Dimensions[#]] & /@ data[[All, 1]]];
	MonitorFunction[dims, "Data Dimensions: "];

	(*Make and export test data*)
	{testData, testVox} = MakeTestData[testData, 2, patch];
	ExportNii[If[is2D, testData[[All,1]], First@testData], testVox, outName["testSet.nii"]];

	(*prepare a validation set which is 10% of round*)
	validation = batchFunction[<|"BatchSize" -> Round[0.1 roundLength]|>];
	Export[outName["validation.wxf"], validation];
	MonitorFunction[{Length@data, Length@validation}, "data / validation: "];

	(*---------- Train the network ----------------*)

	MonitorFunction[{DateString[], loss}, "Starting training"];

	(*export first itt*)
	ittTrain--;
	monitorFunction[<|"Net"->netIn|>];
	logFile = File[outName[StringReplace[DateString["ISODateTime"], ":" | "-" -> ""] <> ".json"]];

	(*Print progress function*)
	PrintTemporary[Dynamic[Column[{
		Style["Training Round: " <> ToString[ittTrain], Bold, Large], 
		Image[im, ImageSize->400]
	}, Alignment -> Center]]];

	(*train the network*)
	trained = NetTrain[netIn, {batchFunction, "RoundLength" -> roundLength}, All, 
		ValidationSet -> validation, LossFunction -> loss,
		TargetDevice -> tar, WorkingPrecision -> "Mixed",
		MaxTrainingRounds -> rounds - ittTrain, BatchSize -> batch, 
		LearningRate -> learningRate, Method -> {"ADAM", "L2Regularization" -> l2reg, 
			"LearningRateSchedule" -> schedule, "Beta1" -> 0.9, "Beta2" -> 0.99, 
			"Epsilon" -> 10^-5, "GradientClipping" -> 1},
		TrainingProgressFunction -> {monitorFunction, "Interval" -> Quantity[rep, "Rounds"]},
		TrainingProgressReporting -> logFile
	];

	(*---------- Export the network ----------------*)

	netOut = NetExtract[trained["TrainedNet"], "net"];
	Export[outName["trained.wxf"], trained];
	Export[outName["final.wlnet"], netOut];
	Export[outName["final.onnx"], netOut];
]


(* ::Subsubsection::Closed:: *)
(*MakeTestData*)


MakeTestData[data_, n_, patch_] := Block[{testData, len, sel, testDat},
	testData = data[[1]];
	(*figure out how to treat multi channel data - for now just take the first volume*)
	If[ArrayDepth@testData === 4, testData = testData[[All, 1]]];

	(*crop the test data*)
	len = Length@testData;
	If[len > First@patch && Length@patch===3,
		sel = Range @@ Clip[Round[(Clip[Round[len/3 - (0.5 n) First@patch], {0, Infinity}] + {1, n First@patch})], 
			{1, len}, {1, len}];
		testData = First@AutoCropData[testData[[sel]]]
	];

	testData = If[Length@patch===2,
		{NormalizeData[PadToDimensions[#, patch], NormalizeMethod -> "Uniform"]}&/@testData[[;;Min[{len, 9}]]],
		{NormalizeData[PadToDimensions[testData, patch], NormalizeMethod -> "Uniform"]}
	];

	{testData, data[[3]]}
];


(* ::Subsection:: *)
(*AugmentTrainingData*)


(* ::Subsubsection::Closed:: *)
(*AugmentTrainingData*)

Options[AugmentTrainingData] = Options[AugmentTrainingDataI] ={
	"Augment2D" -> False
}

SyntaxInformation[AugmentTrainingData] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

AugmentTrainingData[{dat_?ArrayQ, seg_?ArrayQ}, vox_, opts:OptionsPattern[]] := 
	AugmentTrainingDataI[{dat, seg}, vox, {True, True, True, True, True}, opts]

AugmentTrainingData[{dat_?ArrayQ, seg_?ArrayQ}, vox_, aug_?BooleanQ, opts:OptionsPattern[]] := 
	AugmentTrainingDataI[{dat, seg}, vox, {aug, aug, aug, aug, aug}, opts]

AugmentTrainingData[{dat_?ArrayQ, seg_?ArrayQ}, vox_, aug_?ListQ, opts:OptionsPattern[]] := 
	AugmentTrainingDataI[{dat, seg}, vox, aug, opts]

AugmentTrainingData[dat_?ArrayQ, vox_, opts:OptionsPattern[]] := 
	First@AugmentTrainingDataI[{dat, dat}, vox, {True, True, True, True, True}, opts]

AugmentTrainingData[dat_?ArrayQ, vox_, aug_?BooleanQ, opts:OptionsPattern[]] := 
	First@AugmentTrainingDataI[{dat, dat}, vox, {aug, aug, aug, aug, aug}, opts]

AugmentTrainingData[dat_?ArrayQ, vox_, aug_?ListQ, opts:OptionsPattern[]] := 
	First@AugmentTrainingDataI[{dat, dat}, vox, aug, opts]


AugmentTrainingDataI[{dat_?ArrayQ, seg_?ArrayQ}, vox_, aug_?ListQ, OptionsPattern[]] := Block[{
	datT, segT, cr, w, r, t, s, flip, rot, trans, scale, noise, blur, isNot2D},

	(*prep data*)
	{flip, rot, scale, noise, blur} = aug;
	datT = ToPackedArray[N[dat]];
	segT = ToPackedArray[N[seg]];
	isNot2D = !OptionValue["Augment2D"];

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
			If[rot && Coin[] && isNot2D, RandomReal[{-15, 15}], 0.], 
			If[rot && Coin[] && isNot2D, RandomReal[{-15, 15}], 0.],
			(*no Translation, never needed since how data is padded and cropped its always centered*)
			0., 0., 0.,
			(*Scaling*)
			If[scale && Coin[] && isNot2D, RandomReal[{0.6, 1.6}], 1.],
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


AugmentImageData[im_?ListQ, {rot_, flip_}] := AugmentImageData[#, {rot, flip}]&/@im

AugmentImageData[im_, {rot_, flip_}] := Block[{rt, fl, tr},
	rt = If[rot, RotationTransform[RandomReal[{-90, 90}]Degree], TranslationTransform[{0, 0}]];
	fl = If[flip && RandomChoice[{True, False}], ReflectionTransform[{1, 0}], TranslationTransform[{0, 0}]];
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

GetTrainData[dataSets_, nBatch_, patch_, opts:OptionsPattern[]] := GetTrainData[dataSets, nBatch, patch, False, opts]

GetTrainData[dataSets_, nBatch_, patch_, nClass_, OptionsPattern[]] := Block[{
		itt, datO, segO, dat, seg, vox, augI, aug, nSet, pad, sel, is2D
	},

	itt = 0;
	datO = segO = {};

	(*figure out how to augment the data*)
	is2D = Length[patch]===2;
	{augI, nSet, pad} = OptionValue[{AugmentData, PatchesPerSet, PadData}];
	{aug, is2D} = Which[
		BooleanQ[augI], {augI, is2D}, 
		augI==="2D", {True, True}, 
		augI==="3D", {True, False}, 
		True, {True, is2D}
	];

	(*generate sets untill batch size is reached*)
	While[Length@datO < nBatch,
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
		{dat, seg} = AugmentTrainingData[{dat, seg}, vox, aug, "Augment2D" -> is2D];
		{dat, seg} = PatchTrainingData[{dat, seg}, patch, nSet];
		datO = Join[datO, dat];
		segO = Join[segO, seg];
	];

	(*randomly sample the batch to correct length*)
	sel = RandomSample[Range[Length@datO]][[;; nBatch]];
	datO = datO[[sel]];
	segO = segO[[sel]];

	(*pad the data with extra background if needed*)
	If[IntegerQ[pad], {datO, segO} = AddPadding[datO, segO, pad]];
	datO = NormalizeData[#, NormalizeMethod -> "Uniform"]& /@ datO;
	segO = If[IntegerQ[nClass], ClassEncoder[segO, nClass], segO + 1];
	Thread[Transpose[{ToPackedArray@N@datO}] -> ToPackedArray@Round@segO]
];


AddPadding[dat_, seg_, p_] := Block[{datp, segp, pad},
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


PatchTrainingData[{dat_, seg_}, patch_, n_] := Block[{
		pLen, dLen, patchI, nI, pts, datP, segP, slice, sel
	},
	pLen = Length@patch;
	dLen = Length@dat;
	{patchI, nI} = If[pLen === 2, {Prepend[patch, dLen], 1}, {patch, n}];

	(*by overlapping patches with PatchNumber more random patches per data are created*)
	{datP, pts} = DataToPatches[dat, patchI, nI, PatchNumber -> 2];
	segP = First@DataToPatches[seg, patchI, pts];

	(*for 2D patches randomly select slices from single patch*)
	{datP, segP} = Switch[pLen,
		2,
		slice = RandomSample[Range[1, dLen], Min[{n, dLen}]];
		{ToPackedArray[N[#]] & /@ datP[[1, slice]], ToPackedArray[Round[#]] & /@ segP[[1, slice]]},
		3,
		{ToPackedArray[N[#]] & /@ datP, ToPackedArray[Round[#]] & /@ segP}
	];

	(*only select patches where the data is not 0.*)
	sel = Unitize[Mean[Flatten[#]]& /@ datP];
	{Pick[datP, sel, 1], Pick[segP, sel, 1]}
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

PrepareTrainingData[labFol_?StringQ, outFol_?StringQ, opt:OptionsPattern[]] := PrepareTrainingData[{labFol, labFol}, outFol, opt]

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
	cr = FindCrop[Mask[If[ArrayDepth[dat] === 3, NormalizeData, NormalizeMeanData][dat], 5, MaskDilation -> 1]];
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

CheckSegmentation[seg_] := CheckSegmentation[seg, "label"]

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

MakeChannelClassImage[data_, label_] := MakeChannelClassImage[data, label, MinMax[label], {1,1,1}]

MakeChannelClassImage[data_, label_, {off_, max_}] := MakeChannelClassImage[data, label, {off,max}, {1,1,1}]

MakeChannelClassImage[data_, label_, vox_] := MakeChannelClassImage[data, label, MinMax[label], vox]

MakeChannelClassImage[data_, label_, {off_, max_}, vox_] := Block[{i1, i2},
	i1 = MakeClassImage[label, {off, max}, vox];
	i2 = MakeChannelImage[data, vox];
	ImageCollage[ImageCompose[#, SetAlphaChannel[i1, 0.4 AlphaChannel[i1]]]& /@ i2]
]


(* ::Subsubsection::Closed:: *)
(*MakeClassImage*)


SyntaxInformation[MakeClassImage]={"ArgumentsPattern"->{_, _., _.}};

MakeClassImage[label_] := MakeClassImage[label, Round@MinMax[label], {1,1,1}]

MakeClassImage[label_, {off_?NumberQ, max_?NumberQ}] := MakeClassImage[label, {off, max}, {1,1,1}]

MakeClassImage[label_, vox_?VectorQ] := MakeClassImage[label, Round@MinMax[label], vox]

MakeClassImage[labelI_,{offI_?NumberQ, maxI_?NumberQ}, vox_?VectorQ] := Block[{max, cols, imlab, rat, label, off},
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

MakeChannelImage[data_] := MakeChannelImage[data, {1, 1, 1}]

MakeChannelImage[data_, vox_] := Block[{dat, imdat, rat},
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
	If[outRange === Automatic, 
		outRange = Round[Max[1 /. DeleteCases[ComponentMeasurements[#, "Length"] & /@ dat, {}]]/3]
	];

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
	dout = If[outer === {}, {}, -DistFun[nearFun, vox # & /@ outer]];

	ReverseCrop[SparseArray[Join[Thread[inner -> din], Thread[outer -> dout]], dimC, 0.], dim, cr]
]

DistFun[fun_, pts_] := Sqrt[Total[(Flatten[fun[#, 1] & /@ pts, 1] - pts)^2, {2}]];


(* ::Subsection:: *)
(*ShowTrainLog*)


(* ::Subsubsection::Closed:: *)
(*ShowTrainLog*)


ShowTrainLog[fol_] := ShowTrainLog[fol, 5]

ShowTrainLog[fol_, max_] := DynamicModule[{
		pdat, klist, folder = fol, len, plot, plotf, ymaxv,
		xmin, xmax, ymin, ymax, temp, key, key0, key1, key2, filt, fsize, 
		grid, logp
	},

	{klist, pdat, len} =LoadLog[fol, max];
	pdat = pdat[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];

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
			If[grid, GridLines -> {len, Automatic}, GridLines -> {len, None}], 
			PlotHighlighting -> "Dropline"],

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
		{{grid, True, "Grid"}, {True, False}},
		{{logp, True, "Log"}, {True, False}},

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
			Control[{{xmin, 1,"X min"},1, Dynamic[xmax-1], 1}], "  ", 
			Control[{{xmax,Length[pdat],"X max"}, Dynamic[xmin+1], Dynamic[Length[pdat]], 1}]
		}],
		Row[{
			Control[{{ymin, 0.01, "Y min"}, 0, Dynamic[ymax-0.01]}], "  ",
			Control[{{ymax, 1.1, "Y max"}, Dynamic[ymin+0.01], Dynamic[ymaxv]}]
		}],
		Row[{
			Button["Autoscale X", {xmax, xmax} = {1, Length[pdat]}], "  ",
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


(* ::Subsubsection::Closed:: *)
(*LoadLog*)


LoadLog[fol_, max_] := Block[{files, keys, log, leng},
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
				Button["Browse", inputFile = SystemDialogInput["FileOpen"], 
					Method -> "Queued"]
			}, {
				TextCell["Output File: "], 
				InputField[Dynamic[outputFile], String, 
				FieldHint -> "Enter output file path", 
				FieldSize -> {25, 1}],
				Button["Browse", outputFile = SystemDialogInput["FileSave"], 
					Method -> "Queued"]
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
