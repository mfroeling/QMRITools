(* ::Package:: *)

(* Paclet Info File *)

(* created 2022/01/21*)

Paclet[
    Name -> "QMRITools",
    Version -> "2.8.1",
    WolframVersion -> "13.0+",
    SystemID -> All, 
    Description -> "Toolbox for Quantitative MRI.",
    Creator -> "Martijn Froeling <m.froeling@gmail.com>",
    Support -> "https://github.com/mfroeling/QMRITools",
    Icon -> "Resources/icon.png",
    URL -> "https://www.qmritools.com",
    Extensions -> 
        {
            {"Kernel", Root -> "Kernel", Context -> "QMRITools`"}, 
            {"Documentation", Language -> "English", MainPage -> "Guides/QMRITools"},
            
            {"Resource", Root -> "Resources", Resources -> {{"Logo", "icon.png"}}},
            {"Resource", Root -> "Resources", Resources -> {{"Functions", "All-Functions.nb"}}},
            
            {"Asset", "Root" -> "Applications/Windows-x86-64", "SystemID" -> "Windows-x86-64", "Assets" -> {{"Elastix", "elastix.exe"}}},
            {"Asset", "Root" -> "Applications/MacOSX-x86-64", "SystemID" -> "MacOSX-x86-64", "Assets" -> {{"Elastix", "bin/elastix"}}},
            {"Asset", "Root" -> "Applications/Linux-x86-64", "SystemID" -> "Linux-x86-64", "Assets" -> {{"Elastix", "bin/elastix"}}},
            
            {"Asset", "Root" -> "Applications/Windows-x86-64", "SystemID" -> "Windows-x86-64", "Assets" -> {{"Transformix", "transformix.exe"}}},
            {"Asset", "Root" -> "Applications/MacOSX-x86-64", "SystemID" -> "MacOSX-x86-64", "Assets" -> {{"Transformix", "bin/transformix"}}},
            {"Asset", "Root" -> "Applications/Linux-x86-64", "SystemID" -> "Linux-x86-64", "Assets" -> {{"Transformix", "bin/transformix"}}},
            
            {"Asset", "Root" -> "Applications/Windows-x86-64", "SystemID" -> "Windows-x86-64", "Assets" -> {{"ElastixLib", "ANNlib-5.0.dll"}}},
            {"Asset", "Root" -> "Applications/MacOSX-x86-64", "SystemID" -> "MacOSX-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-4.9.1.dylib"}}},
            {"Asset", "Root" -> "Applications/Linux-x86-64", "SystemID" -> "Linux-x86-64", "Assets" -> {{"ElastixLib", "lib/libANNlib-4.9.so"}}},
            
            {"Asset", "Root" -> "Applications/Windows-x86-64", "SystemID" -> "Windows-x86-64", "Assets" -> {{"DcmToNii", "dcm2niix.exe"}}},
            {"Asset", "Root" -> "Applications/Windows-x86-64", "SystemID" -> "Windows-x86-64", "Assets" -> {{"DcmToNii-2", "dcm2niix-20210317.exe"}}},
            {"Asset", "Root" -> "Applications/MacOSX-x86-64", "SystemID" -> "MacOSX-x86-64", "Assets" -> {{"DcmToNii", "bin/dcm2niix"}}},
            {"Asset", "Root" -> "Applications/Linux-x86-64", "SystemID" -> "Linux-x86-64", "Assets" -> {{"DcmToNii", "bin/dcm2niix"}}}
        }
]


