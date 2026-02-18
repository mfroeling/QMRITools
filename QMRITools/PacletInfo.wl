(* ::Package:: *)

PacletObject[<|
	"Name" -> "QMRITools",
	"Version" -> "4.7.0",
	"WolframVersion" -> "14.0+",
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
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Logo", "icon.png"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Functions", "All-Functions.nb"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Demo", "Demonstrations.nb"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Demo_Unet", "Demo_UNet.nb"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"DemoData", "DemoData.zip"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"ColorData", "SCMv8txt.zip"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"GradientTool", "GradientGUI-v14.cdf"}}},

		(*Neural Networks classification*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"LegSide", "Side_Leg.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"LegPosition", "Pos_Leg.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"Body", "Body.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"ShoulderSide", "Side_Shoulder.wlnet"}}},

		(*Neural Networks segmentation*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegShoulderMuscle", "N2_Shoulder.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegShoulderMuscle2D", "N2_Shoulder_2D.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegThighMuscle", "N5_UpperLeg.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegThighMuscle2D", "N5_UpperLeg_2D.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegLegMuscle", "N6_LowerLeg.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegLegMuscle2D", "N6_LowerLeg_2D.wlnet"}}},

		(*Segmentation train labels*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"HeadNeckTrainLabels", "N1_HeadNeck.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"ShoulderTrainLabels", "N2_Shoulder.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"TorsoTrainLabels", "N3_Torso.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"HipTrainLabels", "N4_Hip.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"LegUpperTrainLabels", "N5_UpperLeg.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"LegLowerTrainLabels", "N6_LowerLeg.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"ArmTrainLabels", "N7_Arm.txt"}}},

		(*Segmentation output labels*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MuscleLabels", "N0_Body.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MuscleLegLabels", "Muscles_leg.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MuscleShoulderLabels", "Muscles_shoulder.txt"}}},

		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MusclesLegAllLabels", "Muscles_leg_all.txt"}}},

		(* ---- OS dependant assets ---- *)

		(*Windows*)
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"Elastix", "elastix.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"Transformix", "transformix.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"ElastixLib", "ANNlib-5.2.dll"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii", "dcm2niix-20241211.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"pigz", "pigz.exe"}}},

		(*windows olf dcm2nii versions*)
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-24", "dcm2niix-20240202.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-23", "dcm2niix-20230411.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-21", "dcm2niix-20210317.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-20", "dcm2niix-20201102.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-19", "dcm2niix-20190902.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-17", "dcm2niix-20171204.exe"}}},

		(*Mac-x86*)
		{"Asset", "SystemID" -> "MacOSX-x86-64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"Elastix", "bin/elastix"}}},
		{"Asset", "SystemID" -> "MacOSX-x86-64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"Transformix", "bin/transformix"}}},
		{"Asset", "SystemID" -> "MacOSX-x86-64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-5.2.dylib"}}},
		{"Asset", "SystemID" -> "MacOSX-x86-64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-5.2.1.dylib"}}},
		{"Asset", "SystemID" -> "MacOSX-x86-64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"DcmToNii", "bin/dcm2niix"}}},

		(*Mac-x86*)
		{"Asset", "SystemID" -> "MacOSX-ARM64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"Elastix", "bin/elastix"}}},
		{"Asset", "SystemID" -> "MacOSX-ARM64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"Transformix", "bin/transformix"}}},
		{"Asset", "SystemID" -> "MacOSX-ARM64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-5.2.dylib"}}},
		{"Asset", "SystemID" -> "MacOSX-ARM64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-5.2.1.dylib"}}},
		{"Asset", "SystemID" -> "MacOSX-ARM64", "Root" -> "Applications/MacOSX-x86-64", "Assets" -> {{"DcmToNii", "bin/dcm2niix"}}},

		(*Linux-x86*)
		{"Asset", "SystemID" -> "Linux-x86-64", "Root" -> "Applications/Linux-x86-64", "Assets" -> {{"Elastix", "bin/elastix"}}},
		{"Asset", "SystemID" -> "Linux-x86-64", "Root" -> "Applications/Linux-x86-64", "Assets" -> {{"Transformix", "bin/transformix"}}},
		{"Asset", "SystemID" -> "Linux-x86-64", "Root" -> "Applications/Linux-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-5.2.so"}}},
		{"Asset", "SystemID" -> "Linux-x86-64", "Root" -> "Applications/Linux-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-5.2.so.1"}}},
		{"Asset", "SystemID" -> "Linux-x86-64", "Root" -> "Applications/Linux-x86-64", "Assets" -> {{"DcmToNii", "bin/dcm2niix"}}}
	}
|>]
