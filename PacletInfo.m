(* Paclet Info File *)

(* created 2022/03/04*)

Paclet[
    Name -> "QMRITools",
    Version -> "2.6.6",
    WolframVersion -> "13.0+",
    Description -> "Toolbox for Quantitative MRI.",
    Creator -> "Martijn Froeling <m.froeling@gmail.com>",
    Support -> "https://github.com/mfroeling/QMRITools",
    Icon -> "icon.png",
    Extensions -> 
        {
            {"Kernel", Root -> "Kernel", Context -> "QMRITools`"}, 
            {"Documentation", Language -> "English", MainPage -> "Guides/QMRITools"}, 
            {"Resource", Root -> "Resources", Resources -> 
                {
                    {"Logo", "icon.png"}
                }}, 
            {"Resource", Root -> "Resources", Resources -> 
                {
                    {"Functions", "All-Functions.nb"}
                }}, 
            {"Asset", "Root" -> "Windows-x86-64", "SystemID" -> "Windows-x86-64", "Assets" -> 
                {
                    {"Elastix", "elastix.exe"}
                }}, 
            {"Asset", "Root" -> "MacOSX-x86-64", "SystemID" -> "MacOSX-x86-64", "Assets" -> 
                {
                    {"Elastix", "elastix"}
                }}, 
            {"Asset", "Root" -> "Linux-x86-64", "SystemID" -> "Linux-x86-64", "Assets" -> 
                {
                    {"Elastix", "elastix"}
                }}
        }
]


