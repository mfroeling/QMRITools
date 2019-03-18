BeginTestSection["InstallTest"]

VerificationTest[(* 1 *)
	CompoundExpression[Set[fol1Q, DirectoryQ[FileNameJoin[List[$UserBaseDirectory, "Applications", "QMRITools"]]]], Set[fol2Q, DirectoryQ[FileNameJoin[List[$BaseDirectory, "Applications", "QMRITools"]]]], Set[folT, Or[fol1Q, And[fol2Q, Not[Equal[fol1Q, fol2Q]]]]]]
	,
	True	
]

VerificationTest[(* 2 *)
	CompoundExpression[Set[pac, Sort[PacletInformation[First[PacletFind["QMRITools"]]]]], Set[pacI, Map[Function[Sort[List[Rule["Name", "QMRITools"], Rule["Version", "2.1.0"], Rule["BuildNumber", ""], Rule["Qualifier", ""], Rule["WolframVersion", "11.0+"], Rule["SystemID", All], Rule["Description", "Toolbox for Quantitative MRI."], Rule["Category", ""], Rule["Creator", "Martijn Froeling <m.froeling@gmail.com>"], Rule["Publisher", ""], Rule["Support", ""], Rule["Internal", False], Rule["Location", FileNameJoin[List[Slot[1], "Applications", "QMRITools"]]], Rule["Context", List["QMRITools`"]], Rule["Enabled", True], Rule["Loading", Manual]]]], List[$UserBaseDirectory, $BaseDirectory]]], Set[pacT, AnyTrue[pacI, Function[SameQ[pac, Slot[1]]]]]]
	,
	True	
]

VerificationTest[(* 3 *)
	CompoundExpression[Set[fol, ReplaceAll["Location", pac]], Set[files, Sort[Map[Function[FileBaseName[Last[FileNameSplit[Slot[1]]]]], FileNames["*.wl", fol]]]], Set[filesI, Sort[List["CardiacTools", "CoilTools", "DenoiseTools", "DixonTools", "ElastixTools", "GeneralTools", "GradientTools", "ImportTools", "IVIMTools", "JcouplingTools", "MaskingTools", "NiftiTools", "PhysiologyTools", "PlottingTools", "ProcessingTools", "RelaxometryTools", "SimulationTools", "TaggingTools", "TensorTools", "VisteTools"]]], Set[fileT, SameQ[files, filesI]]]
	,
	True	
]

VerificationTest[(* 4 *)
	CompoundExpression[Set[List[sys, appsI], Switch[$OperatingSystem, "Windows", List["Windows-x86-64", List["ANNlib-4.9.dll", "dcm2niix.exe", "elastix.exe", "transformix.exe"]], "MacOSX", List["MacOSX-x86-64", List["bin", "dcm2niix", "elastix", "transformix", "lib", "libANNlib-4.9.1.dylib"]], "Unix", List["Linux-x86-64", List["bin", "dcm2niix", "elastix", "transformix", "lib", "libANNlib-4.9.so"]]]], Set[apps, Map[Function[Part[FileNameSplit[Slot[1]], -1]], FileNames["*", FileNameJoin[List[fol, "Applications", sys]], 3]]], Set[appsT, SameQ[apps, appsI]]]
	,
	True	
]

EndTestSection[]
