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

CopyTrainedNetwork::usage = 
"CopyTrainedNetwork[fileIn, loc, dim] copy a trained network to final location."


GetNeuralNet::usage = 
"GetNeuralNet[name] loads a pre-trained neural net by asset name. The net is cached within a session.
GetNeuralNet[file] loads a net from a file path directly.
GetNeuralNet[net] passes through an existing NetGraph or NetChain unchanged.
GetNeuralNet[\"Clear\"] clears the session cache."


DiceSimilarity::usage = 
"DiceSimilarity[ref, pred] gives the Dice Similarity between 1 and 0 of segmentations ref and pred for class equals 1.
DiceSimilarity[x, y, class] gives the Dice Similarity of segmentations ref and pred for class.
DiceSimilarity[x, y, {class, ..}] gives the Dice Similarity of segmentations ref and pred for the list of gives classes."

JaccardSimilarity::usage = 
"JaccardSimilarity[ref, pred] gives the Jaccard Similarity between 1 and 0 of segmentations ref and pred for class equals 1.
JaccardSimilarity[x, y, class] gives the Jaccard Similarity of segmentations ref and pred for class.
JaccardSimilarity[x, y, {class, ..}] gives the Jaccard Similarity of segmentations ref and pred for the list of gives classes."

SurfaceDistance::usage = 
"SurfaceDistance[ref, pred] gives the surface distance of segmentations ref and pred for class 1, assuming isotropic unit voxels. Use the vox argument for correct physical distances.
SurfaceDistance[ref, pred, class] gives the surface distance for the given integer class.
SurfaceDistance[ref, pred, {class, ..}] gives the surface distance for each class in the list.
SurfaceDistance[ref, pred, class, vox] gives the surface distance in millimeters using voxel size vox.
SurfaceDistance[ref, pred, {class, ..}, vox] does the same for a list of classes.
The Method option controls the distance metric: \"Mean\", \"Median\", \"RMS\", \"HD\", \"HD95\" (default), or \"Std\". A list of methods can be given to return multiple metrics at once."

MakeDistanceMap::usage = 
"MakeDistanceMap[mask] makes a distance map of the given mask in voxels. The distance map is negative inside the mask and positive outside the mask.
MakeDistanceMap[mask, vox] makes a distance map of the given mask in the same unit as vox. The distance map is negative inside the mask and positive outside the mask."


SegmentData::usage = 
"SegmentData[data] segments the data using the default \"Legs\" method.
SegmentData[data, what] segments using the specified anatomical region. What can be \"Legs\", \"LegsHip\", \"UpperLegs\", \"LowerLegs\", \"Shoulder\", \"Hip\", or \"Body\".
SegmentData[data, {what, netFile}] uses a custom network file instead of the built-in net."

ApplySegmentationNetwork::usage = 
"ApplySegmentationNetwork[data, net] segments data using net. Data can be an array or a nii file path. Net can be a network name, file, or NetGraph.
ApplySegmentationNetwork[data, net, node] returns the network output at the specified intermediate node rather than the final segmentation.
ApplySegmentationNetwork[{datFol, outFol}, net] processes all nii files in datFol and saves results to outFol.
ApplySegmentationNetwork[{{datFol, outFol}, {inTag, outTag}}, net] uses custom file tags to find input and name output files."

ClassifyData::usage = 
"ClassifyData[data, method] classifies the input data using the Body classification network. 
Method \"Body\" returns {side, positions}. Any method containing \"Side\" returns the detected side only. 
Any other method returns the raw network output."


ShowTrainLog::usage =
"ShowTrainLog[folder] shows an interactive plot of the training log stored in folder as JSON files.
ShowTrainLog[folder, min] filters rounds with fewer than min progress entries. Default is 5.
The interface allows selecting metrics, toggling log scale, Gaussian filtering, and live folder reloading."


TrainSegmentationNetwork::usage =
"TrainSegmentationNetwork[{inFol, outFol}] trains a segmentation network from scratch. Training data should be wxf files in inFol. inFol can also be a list of folders, in which case wxf files from all folders are used. Progress is saved each round to outFol.
TrainSegmentationNetwork[{inFol, outFol}, netCont] continues training. netCont can be \"Start\" to restart, a NetGraph, a wlnet file, or a previous outFol to continue from the last saved network.
Loss functions can be All or a subset of {\"Dice\", \"Jaccard\", \"Tversky\", \"MSD\", \"CE\", \"Focal\", \"TopK\"}."

FreezeEncoderLayers::usage =
"FreezeEncoderLayers[net] freezes the \"start\" stem node and the first 3 encoder blocks of net (\"enc_1\" .. \"enc_3\") for transfer learning,
returning a LearningRateMultipliers spec to pass to NetTrain.
FreezeEncoderLayers[net, n] freezes \"start\" plus the first n encoder blocks instead.
FreezeEncoderLayers[net, n, includeStart] controls whether \"start\" is included, default True \[Dash]
set to False if the input channel count differs from the pretrained net, since \"start\" is then freshly initialized and must stay trainable."

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
"DataToPatches[data, patchSize] creates non-overlapping patches covering data with minimal overlap.
DataToPatches[data, patchSize, n] creates patches with n additional overlap steps per dimension.
DataToPatches[data, patchSize, pts] extracts patches at the given pre-computed ranges pts.
DataToPatches[data, pts] extracts patches at ranges pts without resizing.
Output is {patches, ranges} where ranges can be passed to PatchesToData."

PatchesToData::usage = 
"PatchesToData[patches, ranges] reconstructs data from patches at the given ranges.
PatchesToData[patches, ranges, dim] reconstructs into an array of dimensions dim. Overlapping patches are averaged.
PatchesToData[patches, ranges, dim, labels] reconstructs segmentation patches. For each label only the largest connected component is kept and overlapping labels are resolved. Returns a merged integer segmentation."


AugmentTrainingData::usage = 
"AugmentTrainingData[{data, seg}, vox] augments data and segmentation consistently. 
AugmentTrainingData[{data, seg}, vox, aug] controls augmentation. aug can be True/False for all, or a list {flip, rotate, scale, noise, blur}. Default is all True.
AugmentTrainingData[data, vox] augments data only, returning the data without segmentation."

AugmentImageData::usage = 
"AugmentImageData[image, {rotate, flip}] augments the input image by rotating between -180 and 180 degrees and flipping. The inputs rotate and flip
can be set to True or False.
AugmentImageData[{image, ..}, {rotate, flip}] same but for a list of images."


MakeChannelImage::usage = 
"MakeChannelImage[data] makes a cross-sectional image of the channels data of a training dataset generated by GetTrainData.
MakeChannelImage[data, vox] same but with the aspect ratio determined by vox."

MakeClassImage::usage = 
"MakeClassImage[label] makes a cross-sectional image of the classes label of a training dataset generated by GetTrainData
MakeClassImage[label, {b, n}] same but with explicit definition of background value b and number of classes n. 
MakeClassImage[data, vox] same but with the aspect ratio determined by vox.
MakeClassImage[label, {b, n}, vox] same with explicit definition and aspect ratio definition."

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
"SplitDataForSegmentation[data] splits data for \"Legs\" segmentation, detecting left/right side and upper/lower leg position automatically using a classification network.
SplitDataForSegmentation[data, what] splits for the specified region. What can be \"Legs\", \"LegsHip\", \"UpperLegs\", \"LowerLegs\", \"Shoulder\", \"Hip\", or \"Body\".
SplitDataForSegmentation[data, seg] splits both data and segmentation identically for \"Legs\".
SplitDataForSegmentation[data, seg, what] does the same for the specified region.
Output is {{patches, ranges, dim}, locations}."


MuscleLabelToName::usage =
"MuscleLabelToName[{lab, ..}] converts list of lab, which need to be integers to names using the file GetAssetLocation[\"MusclesLegLabels\"].
MuscleLabelToName[{lab, ..}, file] does the same but uses a user defined ITKSnap label definition file."

MuscleNameToLabel::usage = 
"MuscleNameToLabel[{name, ..}] converts list of muscle names to integer labels using the file GetAssetLocation[\"MusclesLegLabels\"]
MuscleNameToLabel[{name, ..}, file] does the same but uses a user defined ITKSnap label definition file."

ImportITKLabels::usage = 
"ImportITKLabels[] loads the default MusclesLegLabels asset.
ImportITKLabels[file] imports an ITKSnap label file by path or asset name.
ImportITKLabels[file, \"Labels\"] returns label integer -> name replacement rules.
ImportITKLabels[file, \"Names\"] returns name -> label integer replacement rules.
ImportITKLabels[file, \"List\"] returns {names, labels} as separate lists."

SegmentDataGUI::usage = 
"SegmentDataGUI[] is a function that creates a graphical user interface (GUI) for segmenting data. 
It prompts the user to enter the paths for the input and output files, and allows them to select the segmentation type." 


RunMuscleMap::usage = 
"RunMuscleMap[file] will run MuscleMap segmentation on the selected nii file. 
RunMuscleMap[{data, vox}] will run MuscleMap segmentation of the data."


(* ::Subsection::Closed:: *)
(*Options*)


LoadTrainingData::usage =
"LoadTrainingData is an option for TrainSegmentationNetwork. If set to True the training data is loaded from the disk."

UseParallelKernels::usage =
"UseParallelKernels is an option for TrainSegmentationNetwork. If set to True, data loading and batch augmentation run on separate parallel kernels while NetTrain runs concurrently on its own kernel, which can substantially speed up training when augmentation is the bottleneck. If set to {True, n} at most n producer kernels are used, useful for debugging on a small scale. If set to False (default) batches are generated in-process as before."

MonitorInterval::usage =
"MonitorInterval is an option for TrainSegmentationNetwork. It defines how often the training is monitored."

L2Regularization::usage =
"L2Regularization is an option for TrainSegmentationNetwork. It defines the L2 regularization factor."

MultiChannel::usage =
"MultiChannel is an option for TrainSegmentationNetwork, If set to True it will train on multi channel input data. If set to False it will select a random channel."

FreezeEncoderDepth::usage =
"FreezeEncoderDepth is an option for TrainSegmentationNetwork. If set to an integer n, freezes the first n encoder blocks of a pretrained input network
so their weights are not updated during training, useful for transfer learning from a pretrained encoder. Only applies when netCont is not \"Start\".
Default is None, which trains all layers."


PatchesPerSet::usage =
"PatchesPerSet is an option for GetTrainData. Defines how many random patches per dataset are created within the batch."

AugmentData::usage = 
"AugmentData is an option for GetTrainData and TrainSegmentationNetwork. If set True the training data is augmented.
It can also be set to \"2D\" or \"3D\" to control if augmentation is done through plane or only in-plane."

PadData::usage =
"PadData is an option for GetTrainData and TrainSegmentationNetwork.. If set to an integers the that number of slices on the top and bottom of the 
data can be made 0. This is done to learn cut of datasets."

PatchSize::usage =
"PatchSize is an option for TrainSegmentationNetwork. Defines the patch size used in the network training."

RoundLength::usage = 
"RoundLength is an option for TrainSegmentationNetwork. Defines how many batches will be seen during each training round."


MaxMemorySize::usage = 
"MaxMemorySize is an option for SegmentData and ApplySegmentationNetwork. Defines the maximum memory a patch can use in Gb. Higher is better."

SegmentationDimension::usage = 
"SegmentationDimension is an option for SegmentData. The default value is \"3D\" but can also be set to \"2D\" to used 2D slice by slice segmentation."

SegmentationResolution::usage=
"SegmentationResolution is an option for SegmentData. If specified it allows to rescale the data before its processed by the segmentation network. It can be 
{voxData, voxSegment} or {dimension} which are the same inputs as for the RescaleData function which is used internally."

DataPadding::usage = 
"DataPadding is an option for ApplySegmentationNetwork. Defines how much to pad the data patches in all directions."


SplitOverlap::usage = 
"SplitOverlap is an option of SplitDataForSegmentation. Has to be a value between 0 and 0.2 and defines how much the left
and right side overlap."


PatchNumber::usage = 
"PatchNumber is an option for DataToPatches. Can be an integer value >= 0. The larger the number the more overlap the patches have.
The minimal number of patches in each direction is calculated, and then for each dimension the given number is added."

PatchPadding::usage = 
"PatchPadding is an option for DataToPatches. Can be an integer value >= 0. It pads the chosen patch size with the given number."


LabelTag::usage = 
"LabelTag is an option for PrepareTrainingData. It defines the tag used in the filenames of the label data."

DataTag::usage = 
"DataTag is an option for PrepareTrainingData. It defines the tag used in the filenames of the data."

OutputTag::usage = 
"OutputTag is an option for PrepareTrainingData. It defines the tag that will be added to the output training data."

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

SegmentationsPerSlice::usage = 
"SegmentationsPerSlice is an option for PrepareTrainingData. It defines how many segmentations should be present in a slice 
for it to be used in the training data."


DistanceRange::usage =
"DistanceRange is an option for MakeDistanceMap. It defines the range of the distance map outside the segmentation in voxels.
Values can be Automatic, All, or a integer value. If All the distance map is calculated for the whole image. If 0 the distance map is only calculated inside the segmentation."


MuscleMapPythonEnvironment::usage =
"MuscleMapPythonEnvironment is the path to the conda environment of the MuscleMap installation. Fhe folder should contain the python.exe file."

MuscleMapPath::usage = 
"MuscleMapPath is the path to the MuscleMap code that is cloned from github and installed according to the muscle map readme."


(* ::Subsection::Closed:: *)
(*Error Messages*)


TrainSegmentationNetwork::net = "The net input is not \"Start\", a network, a network file, or a previous train folder."

TrainSegmentationNetwork::cont = "Could not find a previous network in the specified folder."

TrainSegmentationNetwork::inp = "The string input given is not a network file or a directory."

TrainSegmentationNetwork::itt = "Not enough iterations specified for training. Remaining iterations are less than 5."

TrainSegmentationNetwork::loss = "Unknown loss function should be one of {\"Dice\", \"MSD\", \"Tversky\", \"CE\", \"Jaccard\", \"Focal\", \"TopK\"."


SurfaceDistance::met = "Method `1` not recognized";


ApplySegmentationNetwork::node = "The node ``` is not part of the network"


RunMuscleMap::noEnv = "Could not automatically locate the MuscleMap conda environment in the usual conda locations. Specify \"PythonEnv\" explicitly.";

RunMuscleMap::noScript = "Could not automatically locate mm_segment.py via pip. Specify \"ScriptPath\" explicitly.";

RunMuscleMap::badEnv = "Python environment `1` does not exist.";

RunMuscleMap::badScript = "Script path `1` does not exist.";

RunMuscleMap::noOutFol = "Could not create output folder `1`.";

RunMuscleMap::runFail = "MuscleMap process failed with exit code `1`. See log: `2`.";

RunMuscleMap::noOutput = "MuscleMap ran successfully but expected output file `1` was not found.";

(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*Locations and grouping*)


$BodyPositionClasses = {"LowerLegs", "Knee", "UpperLegs", "Hip", "Torso", "Shoulder", "HeadNeck"};


$SegmentationLocations = <|
	"LowerLegs" -> <|"Net2D" -> "SegLegMuscle2D", "Net3D" -> "SegLegMuscle3D",
		"TrainLabels" -> "LegLowerTrainLabels", "PositionClasses" -> {"LowerLegs", "Knee"}, 
		"Offset" -> {0, 0}|>,
	"UpperLegs" -> <|"Net2D" -> "SegThighMuscle2D", "Net3D" -> "SegThighMuscle3D",
		"TrainLabels" -> "LegUpperTrainLabels", "PositionClasses" -> {"Knee", "UpperLegs", "Hip"}, 
		"Offset" -> {0, 0}|>,
	"Hip" -> <|"Net2D" -> "SegHipMuscle2D", "Net3D" -> "SegHipMuscle3D",
		"TrainLabels" -> "HipTrainLabels", "PositionClasses" -> {"UpperLegs", "Hip", "Torso"}, 
		"Offset" -> {-5, 5}|>,
	"Torso" -> <|"Net2D" -> Missing["NotImplemented"], "Net3D" -> Missing["NotImplemented"],
		"TrainLabels" -> "TorsoTrainLabels", "PositionClasses" -> {"Hip", "Torso", "Shoulder"}, 
		"Offset" -> {-5, 0}|>,
	"Shoulder" -> <|"Net2D" -> "SegShoulderMuscle2D", "Net3D" -> "SegShoulderMuscle3D",
		"TrainLabels" -> "ShoulderTrainLabels", "PositionClasses" ->{"Torso", "Shoulder", "HeadNeck"}, 
		"Offset" -> {-5, 0}|>,
	"HeadNeck" -> <|"Net2D" -> Missing["NotImplemented"], "Net3D" -> Missing["NotImplemented"],
		"TrainLabels" -> "HeadNeckTrainLabels", "PositionClasses" -> {"Shoulder", "HeadNeck"}, 
		"Offset" -> {-5, 0}|>,
	"Arm" -> <|"Net2D" -> Missing["NotImplemented"], "Net3D" -> "SegArmMuscle3D",
		"TrainLabels" -> "ArmTrainLabels" (*no PositionClasses/Offset*)|>
|>;


$SegmentationGroups = <|
	(*output to full body labels*)
	"Body" -> <|"Locations" -> {"LowerLegs", "UpperLegs", "Hip", "Shoulder"}, "Split" -> "Auto",
		"Classify" -> "Position", "OutputLabels" -> "MuscleLabels"|>,
	"LegsBody" -> <|"Locations" -> {"LowerLegs", "UpperLegs"}, "Split" -> "Find",
		"Classify" -> "Position", "OutputLabels" -> "MuscleLabels"|>,
	"LegsHipBody" -> <|"Locations" -> {"LowerLegs", "UpperLegs", "Hip"}, "Split" -> "Find",
		"Classify" -> "Position", "OutputLabels" -> "MuscleLabels"|>,
	"HipBody" -> <|"Locations" -> {"UpperLegs", "Hip"}, "Split" -> "Auto",
		"Classify" -> "Position", "OutputLabels" -> "MuscleLabels"|>,
	"UpperLegsBody" -> <|"Locations" -> {"UpperLegs"}, "Split" -> "Find",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleLabels"|>,
	"LowerLegsBody" -> <|"Locations" -> {"LowerLegs"}, "Split" -> "Find",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleLabels"|>,
	"ShoulderBody" -> <|"Locations" -> {"Shoulder"}, "Split" -> "Auto",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleLabels"|>,

	(*Output to legacy leg labels*)
	"Legs" -> <|"Locations" -> {"LowerLegs","UpperLegs"}, "Split" -> "Find",
		"Classify" -> "Position", "OutputLabels" -> "MuscleLegLabels"|>,
	"LegsHip" -> <|"Locations" -> {"LowerLegs", "UpperLegs", "Hip"}, "Split" -> "Find",
		"Classify" -> "Position", "OutputLabels" -> "MuscleLegLabels"|>,
	"UpperLegs" -> <|"Locations" -> {"UpperLegs"}, "Split" -> "Find",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleLegLabels"|>,
	"LowerLegs" -> <|"Locations" -> {"LowerLegs"}, "Split" -> "Find",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleLegLabels"|>,
	"Hip" -> <|"Locations" -> {"Hip"}, "Split" -> "Auto",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleLegLabels"|>,

	(*Output to legacy shoulder labels*)
	"Shoulder" -> <|"Locations" -> {"Shoulder"}, "Split" -> "Auto",
		"Classify" -> "Side",  "OutputLabels" -> "MuscleShoulderLabels"|>,

	(*Output to arm labels*)
	"Arm" -> <|"Locations" -> {"Arm"}, "Split" -> "Auto",
		"Classify" -> "None", "OutputLabels" -> "MuscleArmLabels"|>
|>;


(* ::Subsection::Closed:: *)
(*CopyTrainedNetwork*)


CopyTrainedNetwork[file_, loc_, dim_] := Block[{fileIn, fileOut},
	fileOut = FileNameJoin[{First[PacletFind["QMRITools"]]["Location"], "NeuralNetworks", Switch[loc,
		"HeadNeck", "N1_HeadNeck_",
		"Shoulder", "N2_Shoulder_",
		"Torso", "N3_Torso_",
		"Hip", "N4_Hip_",
		"UpperLeg", "N5_UpperLeg_",
		"LowerLeg", "N6_LowerLeg_",
		"Arm", "N7_Arm_",
		_, Return[$Failed]
	] <> dim <> ".wlnet"}];
	fileIn = If[DirectoryQ[file], Last[SortBy[FileNames["*.wlnet",file], FileDate]], If[
		FileExtension[file]==="wlnet", file, Return[$Failed]]];
	CopyFile[fileIn, fileOut, OverwriteTarget -> True]
]


(* ::Subsection::Closed:: *)
(*GetNeuralNet*)


SyntaxInformation[GetNeuralNet] = {"ArgumentsPattern" -> {_}};

(*allows for caching and clearing while GetNeuralNet remains protected*)
GetNeuralNet["Clear"] := (
	ClearAll[GetNeuralNetI];
	GetNeuralNetI[name_] := GetNeuralNetI[name] = NeuralNetFunc[name]
)

GetNeuralNet[net : (_NetChain | _NetGraph)] := net

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


SyntaxInformation[ImportITKLabels] = {"ArgumentsPattern"->{_., _.}};

ImportITKLabels[x___]:=ImportITKLabelsI[x]

ImportITKLabelsI[] := ImportITKLabelsI["MuscleLegLabels"];

ImportITKLabelsI[file_]:=ImportITKLabelsI[file, "List"]

ImportITKLabelsI[file_, outType_] := (*ImportITKLabelsI[file, outType] =*) Block[{fileL, lines, muscleNames, muscleLabels},
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
	Switch[outType,
		"Labels", Thread[muscleLabels->muscleNames],
		"Names", Thread[muscleNames->muscleLabels],
		_, {muscleNames, muscleLabels}
	]
]


(* ::Subsubsection::Closed:: *)
(*MuscleLabelToName*)


SyntaxInformation[MuscleLabelToName] = {"ArgumentsPattern"->{_, _.}};

MuscleLabelToName[num_] := MuscleLabelToName[num, "MuscleLegLabels"]

MuscleLabelToName[num_, file_] := num /. ImportITKLabels[file, "Labels"]


(* ::Subsubsection::Closed:: *)
(*MuscleNameToLabel*)


SyntaxInformation[MuscleNameToLabel] = {"ArgumentsPattern"->{_, _.}};

MuscleNameToLabel[name_] := MuscleNameToLabel[name, "MuscleLegLabels"]

MuscleNameToLabel[name_, file_] := name /. ImportITKLabels[file, "Names"]


(* ::Subsection:: *)
(*ClassifyData*)


(* ::Subsubsection::Closed:: *)
(*ClassifyData*)


Options[ClassifyData] = {
	TargetDevice -> "CPU",
	Monitor -> False
};


SyntaxInformation[ClassifyData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ClassifyData[dat_, met_, OptionsPattern[]] := Block[{
		dev, mask, mon, net, imSize, ims, class, side
	},

	dev = OptionValue[TargetDevice];
	mon = OptionValue[Monitor];

	(*get the network*)
	net = GetNeuralNet["Body"];
	If[net === $Failed, Return[$Failed]];

	(*convert data *)
	mask = Mask[NormalizeData[dat], 10, MaskSmoothing -> True, MaskClosing -> 5];
	imSize = NetDimensions[NetReplacePart[net, "Input"->None], "Input"][[2;;]];
	class = net[MakeClassifyImage[MaskData[dat, mask], ImageSize -> imSize], TargetDevice -> dev];
	side = Last@Keys@Sort@Counts@class["Side"];

	Switch[met,
		"Side", side,
		"Body", {side, FindBodyPos[class["Position"], mon]},
		_, class
	]
]


(* ::Subsubsection::Closed:: *)
(*FindLegPos*)


FindBodyPos[class_] := FindBodyPos[class, False, False]

FindBodyPos[class_, mon_]:=FindBodyPos[class, mon, False]

FindBodyPos[class_, mon_, debug_] := Block[{selection, locations, locationsR, len, classI, classN, n, xVars, eVars, dVars,
	pad, x, e, d, cons, sol, classF, offset, what, lab, pos, found},

	(*pull per-location classifier ranges and offsets from the central location table*)
	locations = Thread[Range[Length[$BodyPositionClasses]] -> $BodyPositionClasses];
	locationsR = Reverse /@ locations;

	len = Length@class;
	pad = Max[{5, Round[0.05 len]}];

	(*create a smoothed and padded list of label numbers*)
	classI = Round@MedianFilter[class /. locationsR, 1];
	classN = Round@MedianFilter[ArrayPad[ArrayPad[classI, {pad, 0}, Min[classI]], {0, pad}, Max[classI]], 1];
	n = Length[classN];

	(*define the model parameters*)
	xVars = Table[x[i], {i, n}];(*solution at each position*)
	dVars = Table[d[i], {i, n - 1}]; (*jump indicators*)
	eVars = Table[e[i], {i, n}];(*error between data and solution*)

	(*define the fit constrains*)
	cons = Join[
		(*define start and end and keep all x larger than 0*)
		{x[1] == Min[classN], x[n] == Max[classN]},
		Table[x[i] >= 0, {i, n - 1}],
		
		(*define jumps and force them between 0 and 1, and no jump can happen within 6 slices*)
		Table[x[i + 1] - x[i] == d[i], {i, n - 1}],
		Table[0 <= d[i] <= 1, {i, n - 1}],
		Table[Total[dVars[[i ;; i + 3]]] <= 1, {i, n - 1 - 3}],

		(*define the minimization error (L1) at each point*)
		Table[x[i] - classN[[i]] <= e[i], {i, n}],
		Table[classN[[i]] - x[i] <= e[i], {i, n}]
	];

	(*perform the fitting where everything should be integers*)
	sol = LinearOptimization[Total[eVars], cons, (Join[xVars, eVars, dVars] \[Element] Integers)];
	classF = ArrayPad[xVars /. sol, {-pad, -pad}];

	If[debug, 
		plot = Show[ListPlot[classI, PlotStyle->PointSize[Large]],ListLinePlot[classF, PlotStyle->Red]];
		MonitorFunction[plot, "Fit"];
	];

	(*figure out which slices belong to which body pos*)
	selection = Select[{#, $SegmentationLocations[#, "PositionClasses"]}& /@ Keys[$SegmentationLocations], ListQ[Last[#]] &];
	offset = $SegmentationLocations[#, "Offset"]&;

	(*classes actually found*)
	found = DeleteDuplicates[classF] /. locations;
	what = Select[(
		{lab, pos} = #;
		pos = pos /. locationsR;
		pos = Flatten[Position[classF /. Append[Thread[pos -> 1], _Integer -> 0], 1]];
		pos = If[pos =!= {} && MemberQ[found, lab], Clip[pos[[{1, -1}]] + offset[lab], {1, len}], {}];
		{lab, pos}
	) & /@ selection, #[[2]] =!= {} &];

	If[mon, MonitorFunction[Row[{found, Column@what}, " | "], "Found locations:"]];

	what
]


(* ::Subsection:: *)
(*Data Patching*)


(* ::Subsubsection::Closed:: *)
(*PatchesToData*)


SyntaxInformation[PatchesToData] = {"ArgumentsPattern" -> {_, _, _., _.}};

PatchesToData[patches_, location_] := PatchesToData[patches, location, Max /@ Transpose[location[[All, All, 2]]], {}]

PatchesToData[patches_, location_, dim : {_?IntegerQ, _?IntegerQ, _?IntegerQ}] := PatchesToData[patches, location, dim, {}]

PatchesToData[patches_, location_, dim : {_?IntegerQ, _?IntegerQ, _?IntegerQ}, labs_?ListQ] := Block[{
		sel, zero, a1, a2, b1, b2, c1, c2, dat, wt,
		seg, lab, pos, si, pi, over
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
			pi = location[[#]] & /@ p[[All, 1]];
			Unitize[PatchesToDataI[si, pi, dim]]
		], {p, pos}];

		(*only keep the largest connected segmentation*) 
		seg = Transpose[SmoothMask[#, MaskComponents -> 1, 
			MaskClosing -> False, SmoothIterations -> 0] &/@ seg];

		(*if needed remove the overlap*)
		If[MinMax[seg] =!= {0, 0}, 
			seg = RemoveMaskOverlaps@seg;
			(*over = 1 - Unitize[Ramp[Total[Transpose@seg] - 1]];
			seg = MaskData[seg, over];*)
		];

		(*merge the segmentations*)	
		MergeSegmentations[seg, labs]
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
	PatchNumber -> 0, 
	PatchPadding -> 0
}


SyntaxInformation[DataToPatches] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, opts:OptionsPattern[]] := DataToPatches[dat, patch, "All", opts] 

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, pts:{{{_,_},{_,_},{_,_}}..}] := {GetPatch[dat, patch, pts], pts}

DataToPatches[dat_, pts:{{{_,_},{_,_},{_,_}}..}] := {GetPatch[dat, pts], pts}

DataToPatches[dat_, patch:{_?IntegerQ, _?IntegerQ, _?IntegerQ}, nPatch_, OptionsPattern[]] := Block[{
		patchOut, pts, nRan, pad
	},
	{nRan, pad} = OptionValue[{PatchNumber, PatchPadding}];
	If[Or @@ (#<=2 pad &/@ patch),
		$Failed,
		pts = GetPatchRanges[dat, patch, If[IntegerQ[nPatch], nPatch, "All"], {nRan, pad}];
		patchOut = GetPatch[dat, patch, pts];
		{patchOut, pts}
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
(*SetMXenvironment*)

SetMXenvironment[what_]:=SetEnvironment[Switch[what,
	"StartSegment",
	{
		"MXNET_CUDNN_AUTOTUNE_DEFAULT" -> "0",
		"MXNET_GPU_MEM_POOL_TYPE" -> "Unpooled"
	},
	"StartTrain",
	{
		"MXNET_CUDNN_AUTOTUNE_DEFAULT" -> "1",
		"MXNET_GPU_MEM_POOL_TYPE" -> "Round",
		"MXNET_CUDA_ALLOW_TENSOR_CORE" -> "1",
		"MXNET_CUDA_TENSOR_OP_MATH_ALLOW_CONVERSION" -> "1"
	},
	"Reset",
	{
		"MXNET_CUDNN_AUTOTUNE_DEFAULT" -> "1",
		"MXNET_GPU_MEM_POOL_TYPE" -> "Round",
		"MXNET_CUDA_ALLOW_TENSOR_CORE" -> "1",
		"MXNET_CUDA_TENSOR_OP_MATH_ALLOW_CONVERSION" -> "0"
	}
]];


(* ::Subsubsection::Closed:: *)
(*SegmentData*)


Options[SegmentData] = {
	TargetDevice -> "CPU",
	MaxMemorySize->Automatic,
	SegmentationDimension -> "3D",
	SegmentationResolution -> Automatic,
	Monitor->False
};


SyntaxInformation[SegmentData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SegmentData[datI_, opts:OptionsPattern[]] := SegmentData[datI, "Body", opts]

SegmentData[datI_, what_, OptionsPattern[]] := Block[{
		dev, max, mon, patch, pts, dim ,loc, net, seg, all, data, mask, 
		time, timeAll, netFile, monO, sDim, dimI, rescale, labs
	},

	SetMXenvironment["StartSegment"];

	timeAll = First@AbsoluteTiming[
		{dev, max, mon, sDim, rescale} = OptionValue[{TargetDevice, MaxMemorySize, Monitor, 
			SegmentationDimension, SegmentationResolution}];
		sDim = If[MemberQ[{"2D", "3D"}, sDim], sDim, "2D"];
		monO = mon;
		mon = If[mon, MonitorFunction, List];

		mon["--------------------"];
		mon[what, "Segmenting location:"];
		mon[Dimensions@datI, "Analyzing the data with dimensions:"];
		data = Switch[ArrayDepth[datI], 4, datI[[All,1]], _, datI];
		dimI = Dimensions@data;

		(*figure out if the data needs rescaling*)
		If[rescale =!= Automatic, If[rescale[[1]]=!=rescale[[2]],
			data = RescaleData[data, rescale, InterpolationOrder -> 1];
			mon[Dimensions@data, "Data is rescaled to:"];
		]];

		mask = Mask[NormalizeData[data], 10, MaskSmoothing -> True, MaskClosing -> 5];
		data = MaskData[data, mask];

		(*split the data in anatomical based patches for segmentation*)
		time = First@AbsoluteTiming[
			{{patch, pts, dim}, loc} = SplitDataForSegmentation[data, what, 
				Monitor -> monO, TargetDevice -> dev];
		];

		mon["--------------------"];
		mon[Round[time, .1], "Total time for analysis [s]: "];
		mon["--------------------"];
		mon[Column@Thread[{loc, Dimensions/@ patch}], "Segmenting \""<>what<>"\" locations with dimensions:"];

		(*check if the network is correct*)
		GetNetwork[loc_, sDim_] := Lookup[$SegmentationLocations[loc], "Net"<>sDim, $Failed];
		If[MemberQ[GetNetwork[#[[1]], sDim]& /@ loc, $Failed], Return[$Failed]];

		(*Perform the segmentation*)
		time = First@AbsoluteTiming[
			seg = MapThread[(
				net = GetNetwork[#2[[1]], sDim];
				mon["--------------------"];
				mon[{#2, net}, "Performing segmentation for: "];
				seg = ApplySegmentationNetwork[#1, net, TargetDevice -> dev, 
					MaxMemorySize -> max, Monitor -> monO];
				ReplaceLabels[seg, #2, what]
			) &, {patch, loc}]];
		mon["--------------------"];
		mon[Round[time, .1], "Total time for segmentations [s]: "];
		mon["--------------------"];

		(*Merge all segmentations for all expected labels*)
		labs = GetSegmentationLabels /@ seg;
		all = Select[Sort[DeleteDuplicates[Flatten[labs]]], IntegerQ];
		mon[Column[labs], "Putting together the segmentations with labels: "];

		(*after this only one cluster per label remains*)
		time = First@AbsoluteTiming[seg = PatchesToData[seg, pts, dim, all]];
		mon[Round[time, .1], "Total time for final evaluation [s]: "];
		mon["--------------------"];

		time = First@AbsoluteTiming[If[rescale =!= Automatic, If[rescale[[1]]=!=rescale[[2]], 
			seg = RescaleSegmentation[seg, dimI];
			mon[Round[time, .1], "Total time for rescaling [s]: "];
			mon["--------------------"];
		]]];
	];

	SetMXenvironment["Reset"];
	mon[Round[timeAll, .1], "Total evaluation time [s]: "];
	mon["--------------------"];

	seg
]


(* ::Subsubsection::Closed:: *)
(*ReplaceLabels*)


ReplaceLabels[seg_, locI_, what_] := Block[{
		loc, side, labIn, fIn, fOut, labNam, labOut, labOutS
	},

	{loc, side} = locI;

	If[KeyExistsQ[$SegmentationLocations, loc] && Max[seg] > 1,

		(*which labels were used for training and output*)
		labIn = GetSegmentationLabels[seg];
		fIn = $SegmentationLocations[loc, "TrainLabels"];
		fOut = $SegmentationGroups[what, "OutputLabels"];

		(*some labels have side encoding some don't that is figured out here*)
		labNam = MuscleLabelToName[labIn, fIn];
		labOut = MuscleNameToLabel[labNam, fOut];
		labOutS = MuscleNameToLabel[(# <> "_" <> side & /@ labNam), fOut];
		labOut = Cases[Transpose[{labOut, labOutS}], _Integer, 2];

		(*Replace train labels with final labels*)
		ReplaceSegmentations[seg, labIn, labOut]
		,
		seg
	]
]


(* ::Subsubsection::Closed:: *)
(*SplitDataForSegmentation*)


Options[SplitDataForSegmentation] = {
	Monitor -> False,
	TargetDevice -> "CPU",
	SplitOverlap -> 0.05
};


SyntaxInformation[SplitDataForSegmentation] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};


SplitDataForSegmentation[data_?ArrayQ, seg_?ArrayQ, opt:OptionsPattern[]] := SplitDataForSegmentation[data, seg, "Legs", opt]

SplitDataForSegmentation[data_?ArrayQ, seg_?ArrayQ, what_?StringQ, opt:OptionsPattern[]] := Block[{dat,pts,dim,loc, segPatch},
	{{dat, pts, dim}, loc} = SplitDataForSegmentation[data, what, opt];
	segPatch = GetPatch[seg, pts];
	{{dat, pts, dim}, {segPatch, pts, dim}, loc}
]


SplitDataForSegmentation[data_?ArrayQ, opt:OptionsPattern[]] := SplitDataForSegmentation[data, "Legs", opt]

SplitDataForSegmentation[data_?ArrayQ, what_?StringQ, opt:OptionsPattern[]] := Block[{
		dim, whatSide, side, whatPos, pos, dat, right, left, cut, pts, loc, 
		time, monO, mon, overP, dev, over, locs
	},
	dim = Dimensions[data];
	{monO, dev, overP} = OptionValue[{Monitor, TargetDevice, SplitOverlap}];
	mon = If[monO, MonitorFunction, List];

	(*Based on the body location tag decide what to do with the classification*)
	{whatSide, whatPos} = Switch[$SegmentationGroups[what, "Classify"],
		"Position", ClassifyData[data, "Body", TargetDevice -> dev, Monitor -> monO],
		"Side", {ClassifyData[data, "Side", TargetDevice -> dev, Monitor -> monO], 
			{#, {1, dim[[1]]}}& /@ $SegmentationGroups[what, "Locations"]},
		"None", {"Both", {#, {1, dim[[1]]}}& /@ $SegmentationGroups[what, "Locations"]},
		_, Return[$Failed];
	];
	mon[whatSide, "Data contains sides: "];

	(*based on side cut data or propagate*)
	dat = Switch[whatSide,
		(*both sides which need to be split*)
		"Both",
		{cut, over} = Switch[$SegmentationGroups[what, "Split"],
			"Find", {FindMiddle[data], Round[0.5 overP Last@dim]},
			"Auto", Round[{0.5, overP} Last@dim],
			_, {$Failed, $Failed}
		];
		(*cut the data*)
		{right, left, cut} = CutData[data, {cut, over}];
		{{right, {"Right", {1, cut + over}}}, {left, {"Left", {cut + 1 - over, dim[[3]]}}}},

		(*only one side, no split*)
		_, {{data, {whatSide, {1, dim[[3]]}}}}
	];

	(*Select only the locations in that region and the correct position to be segmented*)
	locs = $SegmentationGroups[what, "Locations"];
	whatPos = Select[whatPos, MemberQ[locs, #[[1]]]&& #[[2,1]]=!=#[[2,2]]&];
	mon[whatPos[[All, 1]], "Selected positions \""<>what<>"\": "];

	(*loop over the locations and select the correct data*)
	dat = Flatten[Transpose[(
		{dat, side} = #;
		{dat[[#[[2,1]];;#[[2,2]]]], #, side} & /@ whatPos
	) & /@ dat], 1];

	(*output the selected data with the correct label and coordinates*)
	{dat, pts, loc} = Transpose[CropPart /@ dat];
	{{dat, pts, dim}, loc}
]


(* ::Subsubsection::Closed:: *)
(*CropPart*)


CropPart[{dat_, {loc_, {locStart_, locEnd_}}, {side_, {sideStart_, sideEnd_}}}] := Block[{data, crp},
	{data, crp} = AutoCropData[dat, CropPadding->0];
	{data, Partition[crp, 2] + {locStart-1, 0, sideStart-1}, {loc, side}}
]


(* ::Subsubsection::Closed:: *)
(*ApplySegmentationNetwork*)


Options[ApplySegmentationNetwork] = {
	TargetDevice->"CPU", 
	DataPadding->0, 
	MaxMemorySize->Automatic, 
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
		dev, pad , lim, data, crp, net, dim, time, nodes, mem, dimO, size,
		patch, pts, seg, mon, lab, nClass, precision, normF, is2D
	},

	{dev, pad, lim, mon} = OptionValue[{TargetDevice, DataPadding, MaxMemorySize, Monitor}];
	If[lim === Automatic, lim = If[dev==="CPU", 32, 8]];(*memory limit in GB*)
	mon = If[mon, MonitorFunction, List];

	precision = If[(dev=!="CPU") && ($OperatingSystem === "Windows"), "Mixed", "Real32"];

	data = Which[
		StringQ[dat], First@ImportNii[dat], 
		TensorQ[{dat}, NumericQ], If[First@Dimensions@dat===1, First@dat, dat],
		True, Return[$Failed]
	];
	If[ArrayDepth[data] === 4, data = data[[All, 1]]];

	dimO = Dimensions@data;
	{data, crp} = AutoCropData[data, CropPadding->0];
	data = N@ArrayPad[data, pad, 0.];
	dim = Dimensions[data];

	net = GetNeuralNet[netI];

	If[net===$Failed, $Failed,
		(*calculate the patch size for the data *)
		is2D = Length[Rest[NetDimensions[net, "Input"]]] === 2;
		{mem, size} = FindPatchDim[net, dim, lim];
		mem = If[mem < 1, ToString[Round[1000 mem]] <> " MB", ToString[Round[mem, .1]] <> " GB"];

		(*create the network*)
		net = ChangeNetDimensions[net, "Dimensions" ->If[is2D, Rest@size, size]];
		nClass = NetDimensions[net,"Output"][[-1]];

		(*dataNormalization function*)
		normF = NormalizeData[#, NormalizeMethod -> "Uniform"]&;

		(*create the patches*)
		{patch, pts} = DataToPatches[data, size, PatchNumber -> 0, PatchPadding->pad];
		mon[{size, Length@patch}, "Patch size and created number of patches is:"];
		mon[mem, "Estimated memory need is:"];
		patch = If[is2D, ({#}& /@ normF[#])& /@ patch, {normF[#]}&/@patch];

		(*perform the segmentation*)
		If[node==="",
			time = First@AbsoluteTiming[
				(*actually perform the segmentation with the NN*)
				seg = net[#, TargetDevice->dev, WorkingPrecision ->precision]&/@patch;
				seg = ClassDecoder /@ seg;
				(*reverse all the padding and cropping and merged the patches if needed*)
				seg = PatchesToData[ArrayPad[#, -pad] & /@ seg, Map[# + {pad, -pad} &, pts, {2}], dim, Range[nClass]];
				seg = ReverseCrop[ArrayPad[seg, -pad], dimO, crp];
			];
			mon[{Dimensions[seg], If[Max[seg]<1,{},MinMax[GetSegmentationLabels[seg]]]}, "Output segmentations dimensions and labels: "];
			mon[Round[time, .1], "Time for segmentation [s]: "];

			,
			(*perform the segmentation on a specific node*)
			(*check if node is part of the network*)
			nodes = DeleteDuplicates[Keys[Information[net, "Layers"]][[All, 1]]];
			If[!MemberQ[nodes, node],
				Message[ApplySegmentationNetwork::node, node];
				Return[$Failed]
				,
				seg = NetTake[net, node][{normF[First@patch]}, TargetDevice->dev, WorkingPrecision ->precision];
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


FindPatchDim[net_, dim_] := FindPatchDim[net, dim, 8]

FindPatchDim[net_, dim_, lim_] := Block[{
		dz, dy, dx, inp, sz, sy, sx, base, max, out, rat, netMem, x, y, z,
		budget, zMin, zC, yMin, yC, xC, zMax, yMax, xMax
	},

	(*network memory and dimensions*)
	netMem = QuantityMagnitude[UnitConvert[Quantity[16. (
			Information[#, "ArraysTotalElementCount"] + 
			Total[Times @@@ Values[(Information[#, "OutputPorts"]["Output"]) & /@ Information[#, "Layers"]]]
		), "Bits"], "GB"]] &;
	inp = Rest[NetDimensions[net, "Input"]];

	If[Length[inp]===2,
		(*2D network*)
		inp = Rest[NetDimensions[net, "Input"]];
		{sy, sx} = inp / Rest[NetDimensions[net, "MinEncodingOut"]];

		(*make the images as big as possible*)
		out = {dim[[1]], Ceiling[dim[[2]], sy], Ceiling[dim[[3]], sx]};
		{Round[netMem[ChangeNetDimensions[net, "Dimensions" -> Rest@out]], .01], out}

		,
		(*3D network*)
		(*figure out net dimensions and allowed steps*)
		{dz, dy, dx} = inp;
		{sz, sy, sx} = inp / Rest[NetDimensions[net, "MinEncodingOut"]];
		out = {Ceiling[dim[[1]], sz], Ceiling[dim[[2]], sy], Ceiling[dim[[3]], sx]};
		rat = N[out[[2]] / out[[3]]];
		
		(*figure out net memory and memory steps*)
		base = netMem[net];
		budget = 0.9 lim;
		z = zMax = out[[1]];
		(*make sure largest of yx dims is used first*)
		{y, x} = {yMax, xMax} = If[rat > 1, out[[{2, 3}]], 
			{dy, dx} = {dx, dy}; {sy, sx} = {sx, sy}; out[[{3, 2}]]];

		(*first shrink z then y then x until fit in memory budget*)
		If[base (z y x) / (dz dy dx) > budget,
			(* stage 1: shrink z until it hits the thin limit *)
			zMin = Mean[{y, x}] / 4;
			zC = budget (dz dy dx) / (base y x);
			z = Clip[dz + sz Floor[(If[zC > zMin, Min[zC, zMax], zMin] - dz) / sz], {dz, zMax}];
			(* stage 2: if z was not enough shrink y until y / x hits 75% of initial ratio*)
			If[zC <= zMin,
				yMin = 0.75 x;
				yC = budget dz dy dx/(base z x);
				y = Clip[dy + sy Floor[(If[yC > yMin, Min[yC, yMax], yMin] - dy) / sy], {dy, yMax}];
				(* stage 3: shrink x as a last resort *)
				If[yC <= yMin,
					xC = budget dz dy dx/(base z y);
					x = Clip[dx + sx Floor[(xC - dx) / sx], {dx, xMax}];
				];
			];
		];

		out = If[rat > 1, {z, y, x}, {z, x, y}];
		{Round[netMem[ChangeNetDimensions[net, "Dimensions" -> out]], .1], out}
	]
]


(* ::Subsection:: *)
(*TrainSegmentationNetwork*)


(* ::Subsubsection::Closed:: *)
(*TrainSegmentationNetwork*)


Options[TrainSegmentationNetwork] = {
	LoadTrainingData -> True,
	UseParallelKernels -> False,
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
	FreezeEncoderDepth -> None,

	AugmentData -> True,
	PadData-> False,

	LossFunction -> All,
	DropoutRate -> 0.2,
	LearningRate -> 0.001,
	L2Regularization -> 0.0001,

	Monitor -> False,
	TargetDevice -> "GPU"
}


SyntaxInformation[TrainSegmentationNetwork] = {"ArgumentsPattern" -> {{_, _}, _., OptionsPattern[]}};

TrainSegmentationNetwork[{inFol : (_?StringQ | {__?StringQ}), outFol_?StringQ}, opts : OptionsPattern[]] := TrainSegmentationNetwork[{inFol, outFol}, "Start", opts]

TrainSegmentationNetwork[{inFol : (_?StringQ | {__?StringQ}), outFol_?StringQ}, netCont_, opts : OptionsPattern[]] := Block[{
		netOpts, batch, roundLength, rounds, data, depth, nChan, nClass, outName, ittString, multi,
		patch, augment, netIn, ittTrain, testData, testVox, testSeg, im, patches, pLen, is2D,
		monitorFunction, netMon, netOut, trained, l2reg, pad, batchFunction, trainFunc, trainOpts, base,
		validation, files, loss, rep, learningRate, schedule, dims, tar, logFile, allOpts, parallel,
		nProducers, loadData, queueVars, produced, used, ready, activeProducers, trainingDone, index,
		nVal, makeVal, maxProducers, producerStatus, roundImage, freezeDepth, chanIn, lrMult
	},

	SetMXenvironment["StartTrain"];

	(*------------ Get all the configuration stuff -----------------*)

	(*getting all the options*)
	netOpts = Join[FilterRules[{opts}, Options@MakeUnet],
		FilterRules[Options@TrainSegmentationNetwork, Options@MakeUnet]];

	{batch, roundLength, rounds, augment, pad, patch, patches,
		loss, rep, learningRate, l2reg, multi, tar, freezeDepth} = OptionValue[
		{BatchSize, RoundLength, MaxTrainingRounds, AugmentData, PadData, PatchSize,
			PatchesPerSet, LossFunction, MonitorInterval, LearningRate, L2Regularization,
			MultiChannel, TargetDevice, FreezeEncoderDepth}];

	(*False, True, or {True,n} to cap producer count*)
	parallel = OptionValue[UseParallelKernels];
	{parallel, maxProducers} = Which[
		parallel === True, {True, Automatic},
		MatchQ[parallel, {True, _Integer?Positive}], {True, Last[parallel]},
		True, {False, Automatic}
	];

	(*export the training settings*)
	allOpts = Options[TrainSegmentationNetwork][[All, 1]];
	allOpts = Association[Thread[allOpts->OptionValue[allOpts]]];
	Export[FileNameJoin[{outFol,FileBaseName[outFol]<>"_settings.txt"}], 
		KeyValueMap[ToString[#1] <> " -> " <> ToString[#2] &, allOpts], "List"];

	pLen = Length@patch;
	is2D = pLen===2;
	pad = If[NumberQ[pad] && !is2D, Round[pad], False];

	(*get the train data files*)
	files = Flatten[FileNames["*.wxf", #] & /@ Flatten[{inFol}]];

	(*figure out network properties from train data*)
	testData = Normal/@Import[First@files];
	(*figure out how to treat multi channel data*)	
	depth = ArrayDepth@testData;
	nChan = If[depth === 4 && multi, Length@First@testData, 1];
	(*background is 0 but for network its 1 so class +1*)
	nClass = Round[Max@testData[[2]] + 1];
	{testData, testVox} = MakeTestData[testData, 2, patch];

	(*------------ Define the network -----------------*)

	MonitorFunction[DateString[], "Preparing the network: "];

	(*netCont: a network, or a previous train folder*)
	{netIn, ittTrain} = Which[
		(*start with a clean network*)
		netCont === "Start",
		{MakeUnet[nChan, nClass, patch, Monitor->OptionValue[Monitor], Sequence@netOpts], 0}
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

	(*match dimensions, classes, and channels to the input*)
	chanIn = First@NetDimensions[netIn, "Input"];
	netIn = NetInitialize[ChangeNetDimensions[netIn,
		"Dimensions" -> patch, "Channels" -> nChan, "Classes" -> nClass],
		Method -> {"Kaiming", "Distribution" -> "Normal"}];

	(*freeze the pretrained encoder for transfer learning, only applies when continuing from an existing network*)
	lrMult = If[netCont =!= "Start" && freezeDepth =!= None,
		FreezeEncoderLayers[netIn, freezeDepth, chanIn === nChan],
		{}
	];
	MonitorFunction[NetSummary[netIn, "Mem"], "Network summary: "];

	(*define the network for training*)
	netIn = AddLossLayer@netIn;
	MonitorFunction[netIn, "Network is ready"];
	MonitorFunction["--------------------"];

	(*---------- Training functions ----------------*)

	(*Local functions*)
	outName = FileNameJoin[{outFol, Last[FileNameSplit[outFol]] <> "_" <> #}]&;
	ittString = "itt_" <> StringPadLeft[ToString[#], 4, "0"]&;
	roundImage[] := Column[{
		Style["Training Round: " <> ToString[ittTrain], Bold, Large],
		Image[im, ImageSize->400]
	}, Alignment -> Center];

	ExportNii[If[is2D, testData[[All, 1]], First@testData], testVox, outName["testSet.nii"]];
	(*Monitor function*)
	With[{mon = <|"outName" -> outName, "ittString" -> ittString, "testData" -> testData,
			"testVox" -> testVox, "nClass" -> nClass, "is2D" -> is2D|>}, 
		monitorFunction = (
			ittTrain++;
			base = mon["outName"][mon["ittString"][#1]<>#2]&;
			(*perform segmentation and export*)
			netMon = NetExtract[#Net, "net"];
			testSeg = Ramp[ClassDecoder[netMon[mon["testData"], TargetDevice -> "CPU"]]];
			ExportNii[testSeg, mon["testVox"], base[ittTrain,".nii"]];
			(*make and export test image*)
			im = MakeChannelClassGrid[If[mon["is2D"], Transpose@mon["testData"], mon["testData"]],
				{testSeg, {0, mon["nClass"]-1}}, 3];
			Export[base[ittTrain,".png"], im, "ColorMapLength" -> 256];
			(*export network, delete the previous itt*)
			Export[base[ittTrain,".wlnet"], netMon];
			Quiet@DeleteFile[base[ittTrain-2,".wlnet"]];
		)&;
	];

	(*batch function*)
	With[{
			patch = {patch, nClass}, 
			batchOpts = Sequence[PatchesPerSet -> patches, AugmentData -> augment, PadData -> pad]
		},
		batchFunction[n_Integer] := batchFunction[<|"BatchSize" -> n|>];
		batchFunction[assoc_?AssociationQ] := GetTrainData[data, assoc["BatchSize"], patch, batchOpts];
	];

	(*oneCycle learning rate schedule function*)
	schedule = OneCycleSchedule[roundLength / batch, rounds, ittTrain];

	(*---------- Load the data ----------------*)

	MonitorFunction[DateString[], "Load the data and make validation: "];
	loadData = OptionValue[LoadTrainingData] === True;

	(*check if validation can be loaded if not do later*)
	nVal = Min[{50, Round[0.2 roundLength]}];
	makeVal = True;
	If[ittTrain > 0 && FileExistsQ[outName["validation.wxf"]],
		validation = Import[outName["validation.wxf"]];
		makeVal = Dimensions[validation[[1, 1, 1]]] =!= patch;
	];

	(*import all train data or train out of memory and create 20% of round as validation*)
	If[parallel,
		nProducers = LaunchTrainingKernels[maxProducers];
		dims = LoadProducerData[files, {nProducers, nVal}, {loadData, makeVal, batchFunction}];
		If[makeVal, {dims, validation} = dims, dims = First@dims];
		,
		data = If[loadData, Import /@ files, files];
		dims = If[loadData, If[ArrayDepth[#] === 3,
			Dimensions[Transpose[{#}]], Dimensions[#]] & /@ data[[All, 1]], {}];
		If[makeVal, validation = batchFunction[nVal]];
	];

	(*export the validation set*)
	If[makeVal, Export[outName["validation.wxf"], validation]];

	(*data logging*)
	dims = If[loadData, MeanRange[#, 0] & /@ Transpose[dims], Missing["NotLoaded"]];
	MonitorFunction[Row[{Length@files, Length@validation, Column[dims]}, " / "], "Data / Validation / Dimensions: "];
	MonitorFunction["--------------------"];

	(*---------- Train the network ----------------*)

	MonitorFunction[DateString[], "Starting training: "];
	MonitorFunction[loss, "Using loss functions: "];

	(*make the train function*)
	logFile = File[outName[StringReplace[DateString["ISODateTime"], ":" | "-" -> ""] <> ".json"]];
	trainOpts = Sequence[
		ValidationSet -> validation, LossFunction -> loss,
		TargetDevice -> tar, WorkingPrecision -> "Mixed",
		MaxTrainingRounds -> rounds - ittTrain, BatchSize -> batch,
		LearningRateMultipliers -> lrMult,
		LearningRate -> learningRate, Method -> {"ADAM", "L2Regularization" -> l2reg,
			"LearningRateSchedule" -> schedule, "Beta1" -> 0.9, "Beta2" -> 0.99,
			"Epsilon" -> 10^-5, "GradientClipping" -> 1},
		TrainingProgressFunction -> {monitorFunction, "Interval" -> Quantity[rep, "Rounds"]},
		TrainingProgressReporting -> logFile];
	trainFunc = With[{netIn = netIn, roundLength = roundLength, trainOpts = trainOpts},
		NetTrain[netIn, {#, "RoundLength" -> roundLength}, All, trainOpts]&];

	(*export first itt*)
	ittTrain--; monitorFunction[<|"Net"->netIn|>];
	(*start training*)
	If[!parallel,
		(*sequential monitor and train*)
		PrintTemporary[Dynamic[roundImage[]]];
		trained = trainFunc[batchFunction];
		,
		(*branch with parallel evaluation, first set all shared variables*)
		queueVars = Table["queue" <> ToString[i], {i, nProducers}];
		Clear /@ queueVars;
		Table[ToExpression[queueVars[[i]] <> " = {}"], {i, nProducers}];
		produced = used = Table[0, nProducers];
		ready = Table[False, nProducers];
		activeProducers = nProducers;
		trainingDone = False;
		index = 0;

		ToExpression["SetSharedVariable[" <> StringRiffle[queueVars, ","] <> "]"];
		SetSharedVariable[produced, used, ready, activeProducers, 
			trainingDone, index, ittTrain, im];

		(*Monitor training and producing*)
		producerStatus[] := Grid[{
			{Style["Active producers: " <> ToString[activeProducers], Bold], SpanFromLeft},
			Join[{Style["Produced:", Bold], Total[produced], " - "}, produced],
			Join[{Style["Used:", Bold], Total[used], " - "}, used],
			Join[{Style["Ready:", Bold], Count[ready, True], " - "}, Boole[ready]]
		}, Alignment -> Center];
		PrintTemporary[Dynamic[producerStatus[]]];
		PrintTemporary[Dynamic[roundImage[]]];

		(*launch and activate the producer and trainer*)
		trained = Last[WaitAll[Append[
			RunProducerKernel[#, batch, batchFunction]& /@ Range[nProducers],
			RunTrainerKernel[trainFunc, nProducers]
		], ProgressReporting -> False]];

		(*Final producer feedback*)
		MonitorFunction[producerStatus[]];
	];

	(*Finalize monitoring and training*)
	MonitorFunction[roundImage[]];
	MonitorFunction[DateString[], "Done training: "];
	SetMXenvironment["Reset"];

	(*---------- Export the network ----------------*)

	netOut = NetExtract[trained["TrainedNet"], "net"];
	Export[outName["trained.wxf"], trained];
	Export[outName["final.wlnet"], netOut];
	Export[outName["final.onnx"], netOut];
]


(* ::Subsubsection::Closed:: *)
(*FreezeEncoderLayers*)


SyntaxInformation[FreezeEncoderLayers] = {"ArgumentsPattern" -> {_, _., _.}};

FreezeEncoderLayers[net_] := FreezeEncoderLayers[net, 3]

FreezeEncoderLayers[net_, n_Integer] := FreezeEncoderLayers[net, n, True]

FreezeEncoderLayers[net_, n_Integer, includeStart_?BooleanQ] := Block[{nodes},
	nodes = Select[Keys[net[[All, 1]]], StringMatchQ[#, "enc_" ~~ DigitCharacter ..] &];
	nodes = Select[nodes, ToExpression[StringDelete[#, "enc_"]] <= n &];
	If[includeStart, nodes = Prepend[nodes, "start"]];
	Thread[nodes -> 0]
]


(* ::Subsubsection::Closed:: *)
(*OneCycleSchedule*)


OneCycleSchedule[br_, rounds_, ittTrain_] := With[{
		n = {0.25, 0.4, 0.95} rounds br, it = ittTrain br
	}, (
	ti = #1 + it; 
	Which[
		ti < n[[1]], Rescale[Cos[Pi ti / n[[1]]], {1, -1}, {1./5, 1.}],
		ti < n[[2]], 1.,
		ti < n[[3]], Rescale[Cos[Pi (ti - n[[2]]) / (n[[3]] - n[[2]])], {1, -1}, {1., 1./10}],
		True, 1./10
	]
)& ]


(* ::Subsubsection::Closed:: *)
(*LaunchTrainingKernels*)


LaunchTrainingKernels[maxProducers_:Automatic] := With[{
		setMX = SetMXenvironment,
		load = If[StringContainsQ[First[PacletFind["QMRITools"]]["Location"], "workspace"],
			"QMRIToolsDev`", "QMRITools`"]
	}, Block[{nKernels},
	CloseKernels[];
	nKernels = Length[If[maxProducers === Automatic, 
		LaunchKernels[], LaunchKernels[maxProducers + 1]]];
	MonitorFunction[nKernels, "Starting parallel kernels: "];
	Quiet@ParallelEvaluate[
		Get[load];
		setMX["StartTrain"];
		Quiet@System`SetSystemOptions["ParallelOptions" -> 
			{"MKLThreadNumber" -> 1, "ParallelThreadNumber" -> 1}]
	, ProgressReporting -> False];
	nKernels - 1
]]


(* ::Subsubsection::Closed:: *)
(*LoadProducerData*)


LoadProducerData[files_, {nP_,nV_}, {loadData_, makeVal_, batchFunction_}] := Block[{
		pFiles, jobs, dims, validation
	},
	(*repeat the files so every worker gets at least 3 and then partition*)
	pFiles = RandomSample[Flatten@Table[files, Ceiling[3 nP/Length[files]]]];
	pFiles = Table[pFiles[[i ;; ;; nP]], {i, nP}];
	MonitorFunction[Length /@ pFiles, "Files per kernel: "];
	(*one submit per worker that loads its own part, plus one empty submit for the fitter kernel*)
	DistributeDefinitions[pFiles, nP, nV, batchFunction];
	jobs = ParallelSubmit[
		data = If[loadData, Import /@ #, #];
		dims = If[loadData, If[ArrayDepth[#] === 3, 
			Dimensions[Transpose[{#}]], Dimensions[#]] & /@ data[[All, 1]], {}];
		validation = If[makeVal, batchFunction[Round[nV / nP]], {}];
		{dims, validation}
	]& /@ pFiles;
	{dims, validation} = Transpose[Most[WaitAll[Append[jobs, ParallelSubmit[{}]], ProgressReporting -> False]]];
	{Join @@ dims, Join @@ validation}
]


(* ::Subsubsection::Closed:: *)
(*RunProducerKernel*)


RunProducerKernel[qi_, batch_, batchFunction_] := Block[{chunk},
	ParallelSubmit[
		While[!trainingDone, If[ready[[qi]], Pause[0.1],
			chunk = batchFunction[batch];
			ToExpression["queue" <> ToString[qi] <> " = QMRITools`SegmentationTools`Private`chunk"];
			produced[[qi]]++;
	(*If[produced[[qi]]===1, Print["run: ",Dimensions/@ #& /@ ToExpression["queue" <> ToString[qi]]]];*)
			ready[[qi]] = True
		]];
		activeProducers--;
	]
]


(* ::Subsubsection::Closed:: *)
(*RunTrainerKernel*)


RunTrainerKernel[trainFunc_, nProducers_] := ParallelSubmit[CheckAbort[
	First@{trainFunc[GetFromBatchQueue[nProducers]&], trainingDone = True}, 
	trainingDone = True; $Aborted
]]


(* ::Subsubsection::Closed:: *)
(*PopBatchQueue*)


GetFromBatchQueue[nProducers_] := Block[{startTime, chunk},
	startTime = AbsoluteTime[];
	Catch[While[True,
		Do[
			index = Mod[index, nProducers] + 1;
			If[ready[[index]],
					used[[index]]++;
	(*If[used[[index]]===1, Print["Get: ",index," ",Dimensions/@ #& /@ ToExpression["queue" <> ToString[index]]]];*)
					ready[[index]] = False;
					Throw[ToExpression["queue" <> ToString[index]]]
			], {nProducers}
		];
		(*no producers left, or none delivered enough new data in time, give up*)
		If[activeProducers <= 0 || AbsoluteTime[] - startTime > 20, Abort[]];
		Pause[0.1]
	]]
]


(* ::Subsubsection::Closed:: *)
(*MakeTestData*)


MakeTestData[data_, n_, patch_] := Block[{testData, len, sel, testDat},
	testData = data[[1]];
	(*figure out how to treat multi channel data - for now just take the first volume*)
	If[ArrayDepth@testData === 4, testData = testData[[All, 1]]];

	(*crop the test data*)
	len = Length@testData;
	If[Length@patch===2,
		sel = If[len < 9, Range[len], Round[Range[1, len, len / 11]][[2 ;; -2]]];
		testData = First@AutoCropData[testData[[sel]]]
		,
		If[len > First@patch && Length@patch===3,
			sel = Range @@ Clip[Round[(Clip[Round[len/3 - (0.5 n) First@patch], {0, Infinity}] + {1, n First@patch})], 
				{1, len}, {1, len}];
			testData = First@AutoCropData[testData[[sel]]]
		]
	];

	testData = If[Length@patch===2,
		{NormalizeData[PadToDimensions[#, patch], NormalizeMethod -> "Uniform"]}&/@testData,
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
		datT = DataTransformation[datT, vox, w, InterpolationOrder -> 1, PadOutputDimensions -> True];
		segT = DataTransformation[segT, vox, w, InterpolationOrder -> 0, PadOutputDimensions -> True];
		(*make sure there is no unneeded background before rest of processing*)
		cr = FindCrop[datT, CropPadding -> 0];
		{datT, segT} = ApplyCrop[#, cr]& /@ {datT, segT};
	];
	(*Augmentations sharpness*)
	If[blur && Coin[], datT = GaussianFilter[datT, RandomReal[{0, 2}]]];
	(*Augmentation of noise*)
	If[noise && Coin[], datT = AddSaltAndRice[datT, RandomReal[{5, 50}], CoinN[] RandomReal[{0.001, 0.01}]]];

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


ReverseC = Compile[{{dat, _Real, 1}}, 
	Reverse[dat]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


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

AugmentImageData[im_, {rotI_, flip_}] := Block[{rot, ang, rt, fl, tr},
	{rot, ang} = If[NumberQ[rotI], {True, rotI}, {rotI, 90}];
	rt = If[rot, RotationTransform[RandomReal[{-ang, ang}]Degree], TranslationTransform[{0, 0}]];
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

GetTrainData[dataSets_, nBatch_, patch_, nClass_, opts:OptionsPattern[]]:= GetTrainData[dataSets, nBatch, {patch, nClass}, opts]

GetTrainData[dataSets_, nBatch_, {patch_, nClass_}, OptionsPattern[]] := Block[{
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

	(*generate sets until batch size is reached*)
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
		dat = If[ArrayDepth[dat] === 4, RandomChoice@Normal@Transpose@dat, Normal@dat];
		seg = Normal@seg;

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

	Thread[
		(NumericArray[{#}, "Real32"] & /@ N[datO]) -> 
		(NumericArray[#, "Byte"] & /@ Round[segO])
	]
];


AddPadding[dat_, seg_, p_] := Block[{datPatch, segPatch, pad},
	Transpose@MapThread[(
		datPatch = #1;
		segPatch = #2;
		If[RandomChoice[{0.3, 0.7} -> {True, False}], 
			pad = RandomInteger[{-p, p}];
			Which[
				pad < 0, datPatch[[pad ;;, All, All]] = segPatch[[pad ;;, All, All]] = 0.,
				pad > 0, datPatch[[;; pad, All, All]] = segPatch[[;; pad, All, All]] = 0.
			]
		];
		{datPatch, segPatch}
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
	OutputTag -> "",
	InputLabels -> Automatic,
	OutputLabels -> Automatic,
	TrainVoxelSize -> Automatic,
	CleanUpSegmentations -> True,
	TestRun -> False,
	SegmentationsPerSlice -> All
}

SyntaxInformation[PrepareTrainingData] = {"ArgumentsPattern" -> {_, _,OptionsPattern[]}};


PrepareTrainingData[labFol_?StringQ, outFol_?StringQ, opt:OptionsPattern[]] := PrepareTrainingData[{labFol, labFol}, outFol, opt]

PrepareTrainingData[{labFol_?StringQ, datFol_?StringQ}, outFol_?StringQ, OptionsPattern[]] := Block[{
		labT, datT, inLab, outLab, test, segFiles, datFiles, name, i, df, voxOut, dimSeg, dimDat,
		seg, err, vox, voxD, dat, im, nl, outFile, gr, clean, legend, head, out, segSl, outT
	},

	{labT, datT, outT, inLab, outLab, test, clean, voxOut, segSl} = OptionValue[{LabelTag, DataTag, OutputTag, 
		InputLabels, OutputLabels, TestRun, CleanUpSegmentations, TrainVoxelSize, SegmentationsPerSlice}];
	{inLab, outLab} = {inLab, outLab} /. Automatic -> {0};
	segSl = segSl /. All -> 0;

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
			{dat, voxD} = ImportNii@First@df;
			If[voxOut === Automatic, voxOut = vox];
			(*Print[{vox, voxD, voxOut}];*)

			dimSeg = Dimensions[seg]; 
			dimDat = If[ArrayDepth[dat] === 4, Dimensions[dat[[All,1]]], Dimensions[dat]];

			(*check dimensions and voxel size*)
			If[vox =!= voxD,
				out = {i++, name, "Data and segmentation have different voxel size."},
				If[dimSeg =!= dimDat,
					out = {i++, name, "Data and segmentation have different dimensions size."},

					(*Prepare and analyze the training data and segmentation*)
					{dat, seg} = PrepTrainData[{dat, seg}, {inLab, outLab}, {vox, voxOut}];

					(*Select slice range*)
					{dat, seg} = SelectTrainData[{dat, seg}, segSl];

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
						outFile = FileNameJoin[{outFol, name<>If[outT === "", "", "_"<>outT]}];
						ExportNii[dat, voxOut, outFile <> "_data.nii"];
						ExportNii[seg, voxOut, outFile <> "_label.nii"];
						Export[outFile <> ".png", im, "ColorMapLength" -> 256];
						Export[outFile <> ".wxf", {NumericArray[dat,"Real32"], NumericArray[seg,"Integer16"], voxOut}, 
							PerformanceGoal -> "Size"(*, Method -> {"PackedArrayRealType" -> "Real32"}*)];
					];

					out
				]
			]
		]
	, {sf, segFiles}];

	(*export the overview of what has happened*)
	legend = Grid[{{}, Join[{""}, Item[Style[#[[1]], White, Bold], Background -> #[[2]]] & /@ {{"hole & n > 1", Red}, {"n > 1", Purple}, {"hole", Blue}}, {""}], {}}, Spacings -> {1, 0.5}];
	head = Style[#, Bold] & /@ {" # ", "Name", "Labels"};

	out = Grid[{{Grid[Append[Prepend[out, head], {"", legend, SpanFromLeft}], 
		Spacings -> {2, 1}, Background -> {None, {{None, LightDarkV[Lighter@LightGray, Darker@Gray]}}}, Alignment -> Center]
		}}, Spacings -> {2, 2}, Background-> LightDarkV[White,GrayLevel[0.1]]];
	Export[FileNameJoin[{outFol, "summary.png"}], ImagePad[Rasterize[out], 6, White]];

	out
]


SelectTrainData[{dat_, seg_}, n_]:=Block[{segPerSlice, min ,max},
	segPerSlice = Total[Max[#] & /@ #] & /@ First[SplitSegmentations[seg]];
	{min, max} = MinMax[Position[UnitStep[segPerSlice - n], 1]];
	{dat[[min;;max]], seg[[min;;max]]}
]


(* ::Subsubsection::Closed:: *)
(*PrepTrainData*)


SyntaxInformation[PrepTrainData] = {"ArgumentsPattern" -> {_, _, _.}};

PrepTrainData[{dat_?ArrayQ, seg_?ArrayQ}] := PrepTrainData[{dat, seg}, {{0}, {0}}, {{1,1,1}, {1,1,1}}]

PrepTrainData[{dat_?ArrayQ, seg_?ArrayQ}, {labI_?VectorQ, labO_?VectorQ}] := PrepTrainData[{dat, seg}, {labI, labO}, {{1,1,1}, {1,1,1}}]

PrepTrainData[{daI_?ArrayQ, segI_?ArrayQ}, {labI_?VectorQ, labO_?VectorQ}, {voxI_?VectorQ, voxO_?VectorQ}] := Block[{
		cr, dat, seg
	},
	(*rescale if needed*)
	{dat, seg} = If[voxI===voxO, 
		{daI, segI}, 
		{RescaleData[daI, {voxI, voxO}], RescaleSegmentation[segI, {voxI, voxO}]}
	];

	(*remove background and normalize data and figure out what to do with multi channel data*)
	cr = FindCrop[Mask[If[ArrayDepth[dat] === 3, NormalizeData, NormalizeMeanData][dat], 5, MaskDilation -> 1]];
	{
		ApplyCrop[dat, cr],
		If[labI === {0},
			ApplyCrop[seg, cr],
			ReplaceSegmentations[ApplyCrop[seg, cr], labI, labO]
		]
	}
]


(* ::Subsubsection::Closed:: *)
(*CheckSegmentation*)


SyntaxInformation[CheckSegmentation] = {"ArgumentsPattern" -> {_, _.}};

CheckSegmentation[seg_] := CheckSegmentation[seg, "label"]

CheckSegmentation[seg_, out_?StringQ] := Block[{arrDep, segs, lab, err},
	arrDep = ArrayDepth[seg];
	If[arrDep === 3, {segs, lab} = SplitSegmentations[seg], lab = Range@Length@First@seg; segs = seg];
	segs = Transpose[segs];
	err = Thread[{CheckSegmentation[segs, 0], CheckSegmentation[segs, 1]}];
	lab = Grid[{Switch[#[[2]],
			{0, 0}, Style[#[[1]], LightDarkV[Black,White], Bold],
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

MakeChannelClassGrid[dat_, lab_, ni_] := Block[{len, n1, n2, labP, ran},
	len = Length@First@dat;
	If[IntegerQ[ni],
		n1 = n2 = Min[{Floor[Sqrt[len]], ni}],
		{n1, n2} = ni;
		While[n1 n2 > len, n1--; n2--;]
	];

	{labP, ran} = If[TensorQ[lab], {lab, MinMax@lab}, lab];
	If[IntegerQ[Round[ran]], ran = {0, ran}];

	RemoveAlphaChannel@ImageAssemble@Partition[
		ImagePad[MakeChannelClassImage[dat[[{1}, #]], labP[[#]], ran], 4, White
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

MakeClassImage[labelI_,{offI_?NumberQ, maxI_?NumberQ}, vox_?VectorQ] := Block[{max, cols, imLab, rat, label, off},
	(*SeedRandom[1345];
		cols = Prepend[ColorData["DarkRainbow"][#]&/@RandomSample[Rescale[Range[off+1, max]]],Transparent];
		cols = Prepend[ColorData["RomaO"][#]&/@Rescale[Range[off+1, max]],Transparent];
	*)
	{label, off, max} = Round[{labelI, offI, maxI}];
	max = Max[{Max[label], max}];
	cols = Prepend[ColorData["RomaO"][#]&/@Rescale[
			Join[Select[Range[off + 1, max], EvenQ], Select[Range[off + 1, max], OddQ]]
		],Transparent];
	imLab = Round@Clip[If[ArrayDepth[label] === 3, label[[Round[Length@label/2]]], label] - off + 1, {1, max + 1}, {1, 1}];
	rat = vox[[{2,3}]]/Min[vox[[{2,3}]]];
	ImageResize[Image[cols[[#]]&/@imLab], Round@Reverse[rat Dimensions[imLab]], Resampling->"Nearest"]
]


(* ::Subsubsection::Closed:: *)
(*MakeChannelImage*)


SyntaxInformation[MakeChannelImage]={"ArgumentsPattern"->{_, _., _.}};

MakeChannelImage[data_] := MakeChannelImage[data, {1, 1, 1}]

MakeChannelImage[data_, vox_] := Block[{dat, imDat, rat},
	dat=Clip[Rescale[data, Quantile[Flatten[data], {0.01, 0.99}]], {0., 1.}];
	(*dat = Rescale[data];*)
	rat = vox[[{2, 3}]] / Min[vox[[{2, 3}]]];
	(
		imDat = #;
		imDat = If[ArrayDepth[#]===3, imDat[[Round[Length@imDat/2]]], imDat];
		ImageResize[Image[imDat], Round@Reverse[rat Dimensions[imDat]], Resampling->"Nearest"]
	) &/@ dat
]



(* ::Subsection:: *)
(*Train data preparation*)


(* ::Subsubsection::Closed:: *)
(*MakeTrainData*)


NormDat[dat_] :=  Block[{q = Quantile[Flatten[dat], 0.9], m = Max[dat]}, 
	If[q <= 0.5 m, If[m === 0., dat, dat/m], If[q === 0., dat, 0.75 dat/q]]];

MakeTrainData[{dat1_?ArrayQ, dat2_?ArrayQ, i_?IntegerQ}] := MakeTrainData[{{dat1, dat2, i}}]

MakeTrainData[datI : {{_?ArrayQ, _?ArrayQ, _?IntegerQ} ..}] := Block[{d1, d2, i, dat, mask},
	dat =  Join @@ (({d1, d2, i} = #;	NormDat[(1 - #) d1 + # d2] & /@ Range[0, 1, 1/(i - 1)]) & /@ datI);
	mask = Mask[NormalizeData[Total@dat], 10, MaskSmoothing -> True, 
		MaskComponents -> 1, MaskClosing -> 50, MaskDilation -> 3];
	MaskData[Transpose@dat, mask]
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


DiceSimilarityC = Compile[{{ref, _Integer, 1}, {pred, _Integer, 1}, {class, _Integer, 0}}, Block[{referenceVector, predictionVector, inter},
	referenceVector = 1 - Unitize[ref - class];
	predictionVector = 1 - Unitize[pred - class];
	inter = Total[referenceVector predictionVector];
	N[(2 inter + 1) / (Total[referenceVector] + Total[predictionVector] + 1)]]
, RuntimeOptions -> {"Speed", "WarningMessages" -> False}, Parallelization -> True];


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


JaccardSimilarityC = Compile[{{ref, _Integer, 1}, {pred, _Integer, 1}, {class, _Integer, 0}}, Block[{referenceVector, predictionVector, inter},
	referenceVector = 1 - Unitize[ref - class];
	predictionVector = 1 - Unitize[pred - class];
	inter = Total[referenceVector predictionVector];
	N[(inter + 1) / (Total[referenceVector] + Total[predictionVector] - inter + 1)]]
, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


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
	"RootMeanSquare" | "RMS", Sqrt[Mean[dist^2]],
	"Max" | "Hausdorff" | "HD", Max@dist,
	"Hausdorff95" | "HD95", Quantile[dist,.95],
	"Std" | "StandardDeviation", StandardDeviation@dist,
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
		dat, dim, datDil, dataCrop, dimC, edge, inner, outer, nearFun, din, dOut, cr, outRange
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

	dataCrop = ApplyCrop[SparseArray@dat, cr];
	dimC = Dimensions[dataCrop];

	edge = GetEdge[dataCrop, vox];
	inner = dataCrop["ExplicitPositions"];
	outer = (datDil - dataCrop)["ExplicitPositions"];

	nearFun = Nearest[edge];
	din = DistFun[nearFun, vox # & /@ inner];
	dOut = If[outer === {}, {}, -DistFun[nearFun, vox # & /@ outer]];

	ReverseCrop[SparseArray[Join[Thread[inner -> din], Thread[outer -> dOut]], dimC, 0.], dim, cr]
]


DistFun[fun_, pts_] := Sqrt[Total[(Flatten[fun[#, 1] & /@ pts, 1] - pts)^2, {2}]];


(* ::Subsection:: *)
(*ShowTrainLog*)


(* ::Subsubsection::Closed:: *)
(*ShowTrainLog*)


SyntaxInformation[ShowTrainLog] = {"ArgumentsPattern" -> {_, _.}};

ShowTrainLog[fol_] := ShowTrainLog[fol, 5]

ShowTrainLog[fol_, max_] := DynamicModule[{
		plotDat, keyList, folder = fol, len, plot, plotFilter, ymaxMax,
		xmin, xmax, ymin, ymax, temp, key, key0, key1, key2, filt, filtSize, 
		grid, logFunc
	},

	{keyList, plotDat, len} = LoadLog[fol, max];
	plotDat = plotDat[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];

	key1 = Select[keyList, ! StringContainsQ[#, "Current"] &];
	key2 = Select[Select[key1, StringContainsQ[#, "Loss"] &], # =!= "RoundLoss" && # =!= "ValidationLoss" &];
	key1 = Select[Complement[key1, key2], # =!= "RoundLoss" && # =!= "ValidationLoss" &];
	key0 = {"RoundLoss", "ValidationLoss"};

	Manipulate[
		plot = Transpose[Values /@ Normal[plotDat[All, key]][[All, All]]];
		plotFilter = If[filt, GaussianFilter[#, filtSize]&/@plot, plot];

		ymaxMax = Max[{1.1, 1.1 If[plot==={}, 5, Max[Select[Flatten@plot,NumberQ]]]}];
		ymax = Min[{ymax, ymaxMax}];

		(* Plot the selected metrics *)
		If[logFunc, ListLogPlot, ListLinePlot][If[key === {}, {}, plotFilter], Joined -> True, 
			PlotLegends -> Placed[key, Right], ImageSize -> 600, PlotRange->{{xmin,xmax} ,{ymin,ymax}},
			If[grid, GridLines -> {len, Automatic}, GridLines -> {len, None}], 
			PlotHighlighting -> "Dropline"],

		(*the controls*)
		{{filt, False, "Filter"}, {True, False}},
		{{filtSize, 5, "FilterSize"}, 1, 10, 1},
		{{grid, True, "Grid"}, {True, False}},
		{{logFunc, True, "Log"}, {True, False}},

		Delimiter,
		(*Control[{{key, {}, ""}, keyList, ControlType -> TogglerBar, Appearance -> "Vertical" -> {Automatic, 4}, BaseStyle -> Medium}],*)
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
			Control[{{xmin, 1, "X min"},1, Dynamic[xmax-1], 1}], "  ", 
			Control[{{xmax, Length[plotDat], "X max"}, Dynamic[xmin+1], Dynamic[Length[plotDat]], 1}]
		}],
		Row[{
			Control[{{ymin, 0.05, "Y min"}, 0, Dynamic[ymax-0.01]}], "  ",
			Control[{{ymax, 5.1, "Y max"}, Dynamic[ymin+0.01], Dynamic[ymaxMax]}]
		}],
		Row[{
			Button["Autoscale X", {xmax, xmax} = {1, Length[plotDat]}], "  ",
			Button["Autoscale Y", {ymin, ymax} = {0, Max[{1.1, 1.1 If[plot==={}, 1, Max[Select[Flatten@plot,NumberQ]]]}]}]
		}],

		Delimiter,
		Row[{
			InputField[Dynamic[folder], String, Enabled -> True, FieldSize -> 50], 
			Button["Browse", 
				temp = SystemDialogInput["Directory", folder];
				If[StringQ[temp], folder = temp; 
					{keyList, plotDat, len} = LoadLog[folder, max];
					plotDat = plotDat[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];
					xmin = 1;
				];
				, ImageSize -> {60, Automatic}, Method->"Queued"]}
		],
		Button["Reload", 
			{keyList, plotDat, len} = LoadLog[folder, max];
			plotDat = plotDat[All, <|#, "LearningRate" -> #["LearningRate"]*1000|> &];
		, ImageSize -> {60, Automatic}, Method->"Queued"],

		{{key, {}}, ControlType -> None}
	]
]


(* ::Subsubsection::Closed:: *)
(*LoadLog*)


LoadLog[fol_, max_] := Block[{files, keys, log, length},
	(* Get a list of log files in the specified folder *)
	files = Sort[FileNames["*.json", fol]];

	(* Read the log files and extract the relevant information *)
	log = Select[(Select[Import[#, "Lines"], StringContainsQ[#, "ProgressFraction"] &] & /@ files), Length[#] > max &];
	length = Accumulate[Length /@ log];

	(* Convert the log data into a dataset *)
	log = "[\n" <> StringDrop[StringRiffle[If[StringTake[#, -1] === "}", # <> ",", #] & /@ Flatten[log, 1], "\n"], -1] <> "\n]";
	log = Dataset[Association /@ Import[Export[FileNameJoin[{$TemporaryDirectory, "log.json"}], log, "text"]]];

	(* Get the unique keys (metrics) in the log data *)
	keys = Sort@DeleteDuplicates[Flatten[Normal@log[All, Keys]]];

	{keys, log, length}
]


(* ::Subsection::Closed:: *)
(*SegmentDataGUI*)


SegmentDataGUI[] := DynamicModule[{inputFile, outputFile}, Block[{dat, vox, seg, status, diag, what},
	NotebookClose[segmentWindow];

	what = "Legs";

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
				PopupMenu[Dynamic[what], {"Legs", "UpperLegs", "LowerLegs", "Shoulder"}]
			},{
				Button["Start Segmentation, please be patient", 
					If[! NiiFileExistQ[inputFile], 
						MessageDialog["Input file could not be found."]
						,
						status = TextCell@"Importing";
						{dat, vox} = ImportNii[inputFile];
						status = TextCell@"Segmenting Data";
						seg = SegmentData[dat, what, TargetDevice -> "CPU"];
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

	segmentWindow = CreateWindow[diag, WindowTitle -> "Muscle segmentation", WindowSize -> All];
];];


(* ::Subsection::Closed:: *)
(*RunMuscleMap*)


Options[RunMuscleMap] = {
	MuscleMapPythonEnvironment -> Automatic,
	MuscleMapPath -> Automatic,
	MuscleMapLabels -> False
};

SyntaxInformation[RunMuscleMap] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

RunMuscleMap[input_, OptionsPattern[]] := Block[{
		pyEnvironment, scriptMuscleMap, inFile, outFolder, logFile,
		args, result, out, outMM, suffix, version, gpu, json,
		namMM, numMM, nam, num, lab
	},

	(*resolve python environment*)
	pyEnvironment = OptionValue[MuscleMapPythonEnvironment];
	If[pyEnvironment === Automatic, pyEnvironment = FindMuscleMapEnv[]];
	If[pyEnvironment === $Failed, Message[RunMuscleMap::noEnv];
	Return[$Failed]];
	If[! FileExistsQ[pyEnvironment], 
	Message[RunMuscleMap::badEnv, pyEnvironment];
	Return[$Failed]];

	(*resolve mm_segment.py script location*)
	scriptMuscleMap = OptionValue[MuscleMapPath];
	If[scriptMuscleMap === Automatic, 
	scriptMuscleMap = FindMuscleMap[pyEnvironment]];
	If[scriptMuscleMap === $Failed, Message[RunMuscleMap::noScript];
	Return[$Failed]];
	If[! FileExistsQ[scriptMuscleMap], 
	Message[RunMuscleMap::badScript, scriptMuscleMap];
	Return[$Failed]];

	(*fixed settings for this first version*)
	suffix = "MMap";
	version = "1.4";
	gpu = "Y";

	(*output folder+log file setup*)
	outFolder = FileNameJoin[{$TemporaryDirectory, "QMRIToolsMM"}];
	If[! DirectoryQ[outFolder], CreateDirectory[outFolder]];
	If[! DirectoryQ[outFolder], Message[RunMuscleMap::noOutFol, outFolder]; 
	Return[$Failed]];
	logFile = FileNameJoin[{outFolder, "MuscleMap.log"}];

	(*prepare input file:
	either use path directly or export from Mathematica data*)
	If[StringQ[input],
		inFile = input,
		inFile = FileNameJoin[{outFolder, "MuscleMap.nii.gz"}];
		ExportNii[input[[1]], input[[2]], inFile];
	];

	(*expected output paths based on MuscleMap's naming convention*)
	outMM = FileNameJoin[{outFolder, StringReplace[FileNameTake[inFile], ".nii.gz" -> "_" <> suffix <> ".nii.gz"]}];
	(*build and run the command*)
	args = {pyEnvironment, scriptMuscleMap, "-i", inFile, "-o", outFolder, 
		"-g", gpu, "-v", version, "-x", suffix};
	result = RunProcess[args, All];
	
	(*always write the log,regardless of success/failure*)
	Export[logFile,
	StringJoin["=== ", DateString[], " ===\n",
		"Command: ", StringRiffle[args, " "], "\n\n",
		"--- STDOUT ---\n", result["StandardOutput"], "\n\n",
		"--- STDERR ---\n", result["StandardError"], "\n"],
		"Text", "OverwriteTarget" -> True];

	(*check exit code before trusting any output files exist*)
	If[result["ExitCode"] != 0, 
	Message[RunMuscleMap::runFail, result["ExitCode"], logFile]; 
	Return[$Failed]];

	(*sanity check the expected output actually exists before importing*)
	If[! FileExistsQ[outMM], Message[RunMuscleMap::noOutput, outMM]; 
	Return[$Failed]];

	(*import the nii and json give the output*)
	out = First@ImportNii[outMM];

	If[!OptionValue[MuscleMapLabels],
		{namMM, numMM} = ImportMMLabels[ Last@FileNames["*contrast_agnostic_wholebody_model.json", DirectoryName[scriptMuscleMap], Infinity]];
		{nam, num} = ImportITKLabels["MuscleLabels"];
		lab = GetSegmentationLabels[out];
		out = ReplaceSegmentations[out, lab, (lab /. Thread[numMM -> namMM]) /. Thread[nam -> num]];
	];
	json = Import[StringReplace[outMM, ".nii.gz" -> ".json"], "RawJSON"];
	{out, json}
]

(*locate mm_segment.py via pip's editable-install metadata*)
FindMuscleMap[pyEnvironment_] := Block[{pipShow, editableLine, repoRoot},
	pipShow = RunProcess[{pyEnvironment, "-m", "pip", "show", "scripts"}, "StandardOutput"];
	editableLine = SelectFirst[StringSplit[pipShow, "\n"], StringStartsQ[#, "Editable project location:"] &];

	If[MissingQ[editableLine] || editableLine === Missing["NotFound"], Return[$Failed]];

	repoRoot = StringTrim[StringDrop[editableLine, StringLength["Editable project location:"]]];
	FileNameJoin[{repoRoot, "scripts", "mm_segment.py"}]
]

(*locate the MuscleMap conda environment in common install locations*)
FindMuscleMapEnv[] := Block[{possibleCondaRoots, envDir},
	possibleCondaRoots = {
		FileNameJoin[{$HomeDirectory, ".conda", "envs"}],
		FileNameJoin[{$HomeDirectory, "miniconda3", "envs"}],
		FileNameJoin[{$HomeDirectory, "anaconda3", "envs"}]
	};
	envDir = SelectFirst[FileNameJoin[{#, "MuscleMap"}] & /@ possibleCondaRoots, DirectoryQ];

	If[MissingQ[envDir], Return[$Failed]];
	FileNameJoin[{envDir, "python.exe"}]
]


ImportMMLabels[file_] := Block[{alt, labMM, anatMM, sideMM, numMM, nameMM},
	alt = {
		"gemelli and quadratus femoris" -> "quadratus femoris",
		"biceps femoris long head" -> "biceps femoris long",
		"biceps femoris short head" -> "biceps femoris short",
		"thoracolumbar multifidus" -> "transversospinalis",
		"extensor digitorum / hallucis longus" -> 
		"extensor digitorum longus",
		"peroneus longus" -> "fibularis longus",
		"semispinalis cervicis and multifidus" -> "semispinalis cervicis"
		};

	labMM = Tabular[Import[file, "RawJSON"]["labels"]];
	anatMM = Normal[labMM[[All, "anatomy"]]] /. alt;
	sideMM = Normal[labMM[[All, "side"]]];
	numMM = Normal[labMM[[All, "value"]]];
	nameMM = (StringReplace[
		StringReplace[StringRiffle[#, " "], " " -> "_"], 
		"_No_side" -> ""] & /@ 
		Transpose[{Capitalize[anatMM, "AllWords"], Capitalize[sideMM]}]);

	{nameMM, numMM}
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
