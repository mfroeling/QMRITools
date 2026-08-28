(* ::Package:: *)

PacletObject[<|
	"Name" -> "QMRITools",
	"Version" -> "4.10.1",
	"WolframVersion" -> "15.0+",
	"SystemID" -> All,
	"Description" -> "Toolbox for Quantitative MRI.",
	"Creator" -> "Martijn Froeling <m.froeling@gmail.com>",
	"Support" -> "https://github.com/mfroeling/QMRITools",
	"Icon" -> "Resources/icon.png",
	"URL" -> "https://www.qmritools.com",

	"Extensions" ->	{
		(*context and documentation*)
		{"Kernel", "Root" -> "Kernel", "Context" -> "QMRITools`"},
		{"Documentation", "Language" -> "English", "MainPage" -> "Guides/QMRITools"},
		(* ---- OS independent assets ---- *)
		(*files that need to be included in the build*)
		{"Asset", "Root" -> "Resources", "Assets" -> {
			{"Logo", "icon.png"},
			{"Functions", "All-Functions.nb"},
			{"Demo", "Demonstrations.nb"},
			{"Demo_Unet", "Demo_UNet.nb"},
			{"DemoData", "DemoData.zip"},
			{"ColorData", "SCMv8txt.zip"},
			{"GradientTool", "GradientGUI-v14.cdf"}
		}},
		(*Neural Networks classification and segmentation*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {
			{"Body", "Body_Pos_Side.wlnet"},
			{"SegShoulderMuscle3D", "N2_Shoulder_3D.wlnet"},
			{"SegShoulderMuscle2D", "N2_Shoulder_2D.wlnet"},
			{"SegHipMuscle3D", "N4_Hip_3D.wlnet"},
			{"SegHipMuscle2D", "N4_Hip_2D.wlnet"},
			{"SegThighMuscle3D", "N5_UpperLeg_3D.wlnet"},
			{"SegThighMuscle2D", "N5_UpperLeg_2D.wlnet"},
			{"SegLegMuscle3D", "N6_LowerLeg_3D.wlnet"},
			{"SegLegMuscle2D", "N6_LowerLeg_2D.wlnet"},
			{"SegArmMuscle3D", "N7_Arm_3D.wlnet"},
			(*Segmentation train labels*)
			{"HeadNeckTrainLabels", "N1_HeadNeck.txt"},
			{"ShoulderTrainLabels", "N2_Shoulder.txt"},
			{"TorsoTrainLabels", "N3_Torso.txt"},
			{"HipTrainLabels", "N4_Hip.txt"},
			{"LegUpperTrainLabels", "N5_UpperLeg.txt"},
			{"LegLowerTrainLabels", "N6_LowerLeg.txt"},
			{"ArmTrainLabels", "N7_Arm_Up.txt"},
			(*Segmentation output labels*)
			{"MuscleLabels", "S1_Muscles_body_side.txt"},
			{"MuscleLegLabels", "S2_Muscles_leg_side.txt"},
			{"MuscleShoulderLabels", "S3_Muscles_shoulder_side.txt"},
			{"MuscleArmLabels", "S4_Muscles_arm_upper_side.txt"},
			{"MusclesAllLabels", "T1_Muscles_body.txt"},
			{"MusclesLegAllLabels", "T2_Muscles_leg.txt"}
		}},
		(* ---- OS dependant assets ---- *)
		(*Windows*)
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {
			{"Elastix", "elastix.exe"},
			{"Transformix", "transformix.exe"},
			{"ElastixLib", "elxANNlib.dll"},
			{"DcmToNii", "dcm2niix-20260416.exe"},
			{"pigz", "pigz.exe"},
			(*windows old dcm2nii versions*)
			{"DcmToNii-25", "dcm2niix-20250506.exe"},
			{"DcmToNii-24", "dcm2niix-20240202.exe"},
			{"DcmToNii-23", "dcm2niix-20230411.exe"},
			{"DcmToNii-21", "dcm2niix-20210317.exe"},
			{"DcmToNii-20", "dcm2niix-20201102.exe"},
			{"DcmToNii-19", "dcm2niix-20190902.exe"},
			{"DcmToNii-17", "dcm2niix-20171204.exe"}
		}},
		(*Mac-x86 and Mac-ARM (both served from the same universal binaries)*)
		{"Asset", "SystemID" -> "MacOSX-x86-64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {
			{"Elastix", "bin/elastix"},
			{"Transformix", "bin/transformix"},
			{"ElastixLib", "lib/libelxANNlib.dylib"},
			{"ElastixLib", "lib/libelxANNlib.1.dylib"},
			{"DcmToNii", "bin/dcm2niix"}
		}},
		{"Asset", "SystemID" -> "MacOSX-ARM64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {
			{"Elastix", "bin/elastix"},
			{"Transformix", "bin/transformix"},
			{"ElastixLib", "lib/libelxANNlib.dylib"},
			{"ElastixLib", "lib/libelxANNlib.1.dylib"},
			{"DcmToNii", "bin/dcm2niix"}
		}}
	}
|>]
