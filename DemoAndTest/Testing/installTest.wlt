BeginTestSection["InstallTest"]

VerificationTest[(* 1 *)
	CompoundExpression[Set[fol1, FileNameJoin[List[$UserBaseDirectory, "Applications", "QMRITools"]]], Set[fol2, FileNameJoin[List[$BaseDirectory, "Applications", "QMRITools"]]], Set[fol1Q, DirectoryQ[fol1]], Set[fol2Q, DirectoryQ[fol2]], Set[folT, Or[fol1Q, And[fol2Q, Not[Equal[fol1Q, fol2Q]]]]]]
	,
	True	
]

VerificationTest[(* 2 *)
	CompoundExpression[Set[pac, PacletInformation[First[PacletFind["QMRITools"]]]], Set[pacT, SameQ[pac, List[Rule["Name", "QMRITools"], Rule["Version", "2.0.1"], Rule["BuildNumber", ""], Rule["Qualifier", ""], Rule["WolframVersion", "11.0+"], Rule["SystemID", All], Rule["Description", "Toolbox for Qunatitative MRI."], Rule["Category", ""], Rule["Creator", "Martijn Froeling <m.froeling@gmail.com>"], Rule["Publisher", ""], Rule["Support", ""], Rule["Internal", False], Rule["Location", "C:\\ProgramData\\Mathematica\\Applications\\QMRITools"], Rule["Context", List["QMRITools`"]], Rule["Enabled", True], Rule["Loading", Manual]]]]]
	,
	True	
]

VerificationTest[(* 3 *)
	CompoundExpression[Set[fol, Part[Select[List[fol1, fol2], DirectoryQ], 1]], Set[files, Sort[Map[Function[FileBaseName[Last[FileNameSplit[Slot[1]]]]], FileNames["*.wl", fol]]]], Set[fileT, SameQ[files, Sort[List["CardiacTools", "CoilTools", "DenoiseTools", "DixonTools", "ElastixTools", "GeneralTools", "GradientTools", "ImportTools", "IVIMTools", "JcouplingTools", "MaskingTools", "NiftiTools", "PhysiologyTools", "PlottingTools", "ProcessingTools", "RelaxometryTools", "SimulationTools", "TensorTools", "VisteTools"]]]]]
	,
	True	
]

VerificationTest[(* 4 *)
	CompoundExpression[Set[folApp, FileNameJoin[List[fol, "Applications"]]], Set[apps, Sort[Map[Function[Part[FileNameSplit[Slot[1]], -1]], Join[FileNames["*.exe", folApp], FileNames["*.dll", folApp]]]]], Set[appsT, SameQ[apps, Sort[List["ANNlib-4.9.dll", "dcm2niix.exe", "elastix.exe", "transformix.exe"]]]]]
	,
	True	
]

EndTestSection[]
