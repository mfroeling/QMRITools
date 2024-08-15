(* ::Package:: *)

PacletObject[<|
	"Name" -> "QMRITools",
	"Version" -> "3.31.0",
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

		(* ---- OS independent assest ---- *)

		(*files that need to be included in the build*)
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Logo", "icon.png"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Functions", "All-Functions.nb"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"Demo", "Demonstrations.nb"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"DemoData", "DemoData.zip"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"ColorData", "SCMv8txt.zip"}}},
		{"Asset", "Root" -> "Resources", "Assets" -> {{"GradientTool", "GradientGUI-v14.cdf"}}},
		
		(*Neural Networks*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegThighMuscle", "Muscles_leg_upper.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegLegMuscle", "Muscles_leg_lower.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"SegLegBones", "Bones_leg_full.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"LegSide", "Leg_side_full.wlnet"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"LegPosition", "Leg_pos_full.wlnet"}}},

		(*muscle segmentation labels*)
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"BonesLegLabels", "Bones_leg.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MusclesLegLabels", "Muscles_leg.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MusclesLegAllLabels", "Muscles_leg_all.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MusclesLegHipLabels", "Muscles_leg_hip.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MusclesLegLowerLabels", "Muscles_leg_lower.txt"}}},
		{"Asset", "Root" -> "NeuralNetworks", "Assets" -> {{"MusclesLegUpperLabels", "Muscles_leg_upper.txt"}}},

		(* ---- OS dependant assest ---- *)

		(*Windows*)
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"Elastix", "elastix.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"Transformix", "transformix.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"ElastixLib", "ANNlib-5.2.dll"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii", "dcm2niix.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"pigz", "pigz.exe"}}},

		(*windows olf dcm2nii versions*)
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-23", "dcm2niix-20230411.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-21", "dcm2niix-20210317.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-20", "dcm2niix-20201102.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-19", "dcm2niix-20190902.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-17", "dcm2niix-20171204.exe"}}},
		{"Asset", "SystemID" -> "Windows-x86-64", "Root" -> "Applications/Windows-x86-64", "Assets" -> {{"DcmToNii-0", "dcm2nii.exe"}}},

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
