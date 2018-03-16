(* ::Package:: *)

(* ::Title:: *)
(*DTITools ExportTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`ExportTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DatWrite::usage = 
"DatWrite[file, data] exports data to *.dat file as binary data using Real32 format."

DTItoolExpInd::usage = 
"DTItoolExpInd[data, file] exports a 3D array data to the file filename DTItool format (*.dat) using DatWrite.
DTItoolExpInd[data, file ,folder] exports data to given file and folder.
DTItoolExpInd[data, file ,folder, add] exports data to given file and folder and adds -add to the filename.";

DTItoolExpTens::usage =
"DTItoolExpTens[tensor] exports a diffustion tensor array to the DTItool format (*.dat).
DTItoolExpTens[tensor, add] exports tensor and adds - add to the filenames.
DTItoolExpTens[tensor, add, folder] exports tensor to the given folder and adds - add to the filenames."

DTItoolExpFile::usage = 
"DTItoolExpFile[file, background, add, voxsize] exports a *.dti text file."

DTItoolExp::usage =
"DTItoolExp[tensor, voxsize] exports tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} and uses XX.dat as background and generates corresponding *dti files.
DTItoolExp[tensor, voxsize, folder] exports tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and uses XX.dat as background and generates corresponding *dti files. 
DTItoolExp[tensor, voxsize, folder, add] exports tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and uses XX.dat as background and generates corresponding *dti files adds - add to the filenames.
DTItoolExp[back, tensor, voxsize] exports background to back.dat and tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} and generates corresponding *dti files. 
DTItoolExp[back, tensor, voxsize, folder] exports background to back.dat and tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and generates corresponding *dti files.
DTItoolExp[back, tensor, voxsize, folder, add] exports background to back.dat and tensor to {XX.dat, YY.dat, ZZ.dat, XY.dat, XZ.dat, YZ.dat} to the given folder and generates corresponding *dti files and adds - add to the filenames."

ExportVol::usage =
"ExportVol[filename, data, voxsize] exports a .vol and .raw file which can be loaded in DTItool 3.0."

ExportMhdRaw::usage = 
"ExportMhdRaw[filename, data, voxsize] exports a .mhd and .raw file which can be loaded in elastix."

ExpHistInd::usage = 
"ExpHistInd[dat, name, text] exports 5 histograms to individual files named name with the addition of text using the function Hist.
dat must be {l1, l2, l3, MD, FA}."

ExpHistAll::usage = 
"ExpHistAll[dat, name, text] exports 5 histograms to one file named name with the addition of text using the function Hist.
dat must be {l1, l2, l3, MD, FA}."

ExpPlotsInd::usage = 
"ExpPlotsInd[dat, name, text, text2] exports 5 error plots to individual files named name with the addition of text using the function ErrorPlot.
text2 is used to label the individual plots.
dat must be {l1, l2, l3, MD, FA}."

ExpPlotsAll::usage = 
"ExpPlotsAll[dat, name, text, text2] exports 5 error plots to one file named name with the addition of text using the function ErrorPlot.
text2 is used to label the individual plots.\
dat must be {l1, l2, l3, MD, FA}."

SaveImage::usage = 
"SaveImage[image] exports graph to image, ImageSize, FileType and ImageResolution can be given as options.
SaveImage[image, \"filename\"] exports graph to image with \"filname\", ImageSize, FileType and ImageResolution can be given as options."


(* ::Subsection:: *)
(*Options*)


BinaryType::usage = 
"BinaryType is an option for ExportVol and must be \"Integer16\" for an integer array and \"Real32\" for a Double array."

ExportFile::usage = 
"ExportFile is an option for ExpHistInd, ExpHistAll, ExpPlotsInd and ExpPlotsAll.\
Default value is jpg. Can be any image file type extention."


(* ::Subsection:: *)
(*Error Messages*)


DTItoolExpInd::dim = "data is a `1` dimensional array and must be 3 dimensional {slices, x, y}."

DTItoolExpTens::dim = "data is a `1` dimensional array and must be 4 dimensional {6 tensor elements, slices, x, y}."

DatWrite::dim = "Dat write is used for exporting data. The data is a `1` dimensional array and must be 3 dimensional {slices, x, y}."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*SaveImage*)


Options[SaveImage] = {ImageSize -> 6000, FileType -> ".jpg", ImageResolution -> 300};

SyntaxInformation[SaveImage] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SaveImage[exp_, opts : OptionsPattern[]] := Module[{input, type},
 	type = OptionValue[FileType];
 	type = If[StringTake[type, 2] === "*.", StringDrop[type,1], If[StringTake[type, 1] === ".",type, "."<>type]];
 	
 	input = FileSelect["FileSave", {"*"<>type}];
 	
 	If[input != "Canceled!",
 		input=If[StringTake[input,-4]===type,input,input<>type];
 		SaveImage[exp, input, opts]
  	];
  ]

SaveImage[exp_, filei_String, OptionsPattern[]] := Module[{file,imsize,res,type},
	type = OptionValue[FileType];
 	type = If[StringTake[type, 2] === "*.", StringDrop[type,1], If[StringTake[type, 1] === ".",type, "."<>type]];
	
	file=If[StringTake[filei,-4]===type||StringTake[filei,-5]===type,filei,filei<>type];
	
	imsize=OptionValue[ImageSize];
	res=OptionValue[ImageResolution];
	
	If[OptionValue[FileType]===".tiff"||OptionValue[FileType]===".tif",
  	  Export[file, Rasterize[exp, ImageResolution -> 2*res],(*RasterSize -> imsize,*) ImageResolution -> res,"ImageEncoding"->"LZW"],
  	  Export[file, Rasterize[exp, ImageResolution -> 2*res],(*RasterSize -> imsize,*) ImageResolution -> res]
	];
	
	Print["File was saved to: " <> file];
  ]


(* ::Subsection:: *)
(*DTItool Export*)


(* ::Subsubsection::Closed:: *)
(*DatWrite*)


SyntaxInformation[DatWrite] = {"ArgumentsPattern" -> {_, _}};

DatWrite[file_?StringQ,data_?ArrayQ]:=
Module[{strm},
	If[!ArrayQ[data,3],
		Message[DatWrite::dim,ArrayDepth[data]],
		strm=OpenWrite[file,BinaryFormat->True];
		BinaryWrite[file,Reverse[Dimensions[data]],"UnsignedInteger16",ByteOrdering->-1];
		BinaryWrite[file,N[Flatten[data]],"Real32",ByteOrdering->-1];
		Close[strm];]
		]


(* ::Subsubsection::Closed:: *)
(*DTItoolExpInd*)


SyntaxInformation[DTItoolExpInd] = {"ArgumentsPattern" -> {_, _, _., _.}};

DTItoolExpInd[data_?ArrayQ,filename_String,folder_String:"",add_String:""]:=
Module[{folderi,addi,file},
	If[!ArrayQ[data,3],
		Message[DTItoolExpInd::dim,ArrayDepth[data]],
		addi=If[add!="","-"<>add,""];
		folderi=If[folder!="",folder<>"\\",""];
		file=folder<>filename<>add<>".dat";
		DatWrite[file,N[data]];
		Print["Exported: ",file]]
		];


(* ::Subsubsection::Closed:: *)
(*DTItoolExpTens*)


SyntaxInformation[DTItoolExpTens] = {"ArgumentsPattern" -> {_, _., _.}};

DTItoolExpTens[data_?ArrayQ,add_String:"",folder_String:""]:=
Module[{dir,addi,folderi,file,files},
	If[!ArrayQ[data,4] || Length[data]!=6,
		Message[DTItoolExpTens::dim,ArrayDepth[data]],
		dir={"XX","YY","ZZ","XY","XZ","YZ"};
		addi=If[add!="","-"<>add,""];
		folderi=If[folder!="",folder<>"\\",""];
		files=MapThread[(file=folder<>#1<>add<>".dat";DatWrite[file,N[#2]];file)&,{dir,data}];
		Print["Exported: ",files];
		]];


(* ::Subsubsection::Closed:: *)
(*DTItoolExp*)


SyntaxInformation[DTItoolExp] = {"ArgumentsPattern" -> {_, _, _., _., _.}};

DTItoolExp[tens_?ArrayQ,vox_?ListQ,folder_String:"",add_String:""]:=DTItoolExpi[{},tens,vox,folder,add];

DTItoolExp[back:{_?ArrayQ,_?StringQ},tens_?ArrayQ,vox_?ListQ,folder_String:"",add_String:""]:=DTItoolExpi[{back},tens,vox,folder,add];

DTItoolExp[back:{{_?ArrayQ,_?StringQ}..},tens_?ArrayQ,vox_?ListQ,folder_String:"",add_String:""]:=DTItoolExpi[back,tens,vox,folder,add];


(* ::Subsubsection::Closed:: *)
(*DTItoolExpi*)


DTItoolExpi[back_,tens_,vox_,folder_,add_]:=
Module[{addi,folderi,file,files,fileb,filesb,print},
	addi=If[add!="","-"<>add,""];
	folderi=If[folder!="",folder<>"\\",""];
	
	print=Reap[
		If[back=={},
			fileb=folderi<>"tens"<>addi;
			DTItoolExpFile[fileb<>".dti","XX",addi,vox];
			Sow["No background given, using XX as background.\nExported: "<>fileb<>".dti\n"];
			,
			filesb=Map[(
				fileb=folderi<>#[[2]]<>addi;
				DatWrite[fileb<>".dat",N[#[[1]]]];
				DTItoolExpFile[fileb<>".dti",#[[2]],addi,vox];
				fileb
				)&,back];
			Sow[StringJoin[
				{"Exported: ",Riffle[(#<>".dat"&/@filesb),", "],"as background.\nExported: ",
				Riffle[(#<>".dti"&/@filesb),", "],"\n"
				}]];
			];
		files=MapThread[(
			file=folderi<>#1<>addi<>".dat";
			DatWrite[file,N[#2]];file)&,{{"XX","YY","ZZ","XY","XZ","YZ"},tens}];
			Sow[StringJoin[{"Exported: ",Riffle[files,", "]}]];
	][[2,1]];
	StringJoin[print]
]


(* ::Subsubsection::Closed:: *)
(*DTItoolExpFile*)


SyntaxInformation[DTItoolExpFile] = {"ArgumentsPattern" -> {_, _, _, _}};

DTItoolExpFile[file_String,back_String,add_String,vox_List]:= 
Module[{},
	Export[file,
"/* DTI BMT format */
T
float
XX XX"<>add<>".dat
YY YY"<>add<>".dat
ZZ ZZ"<>add<>".dat
XY XY"<>add<>".dat
XZ XZ"<>add<>".dat
YZ YZ"<>add<>".dat
I "<>back<>add<>".dat
"<>ToString[vox[[2]]]<>" "<>ToString[vox[[3]]]<>" "<>ToString[vox[[1]]]
		,"Text"]]


(* ::Subsection:: *)
(*Volume files *)


(* ::Subsubsection::Closed:: *)
(*ExportVol*)


Options[ExportVol]={BinaryType->"Integer16"};

SyntaxInformation[ExportVol] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ExportVol[file_String,dat_?ArrayQ,voxsize_List,OptionsPattern[]]:=
Module[{data,type},
	{data,type}=If[Element[Flatten[dat],Integers],
		{dat,
			If[Max[dat] > 32768 || Min[dat] < -32768,
				"Integer32",
				"Integer16"
				]},
		{N[dat],"Real32"}
		];
	ExportVolRaw[data,file,type];
	ExportVolFile[file,Dimensions[data],voxsize,type];
	]


(* ::Subsubsection::Closed:: *)
(*ExportVolRaw*)


ExportVolRaw[data_,file_,bit_]:=
Module[{strm},
	strm=OpenWrite[file<>".raw",BinaryFormat->True];
	BinaryWrite[strm,Flatten[data],bit];
	Close[strm];
	]


(* ::Subsubsection::Closed:: *)
(*ExportVolFile*)


ExportVolFile[file_,dim_,vox_,bit_]:=
Module[{},
	Export[file<>".vol",
"Data.FileName     = "<>Last[StringSplit[file,"\\"]]<>".raw
Data.Type         = raw
Data.Dimensions   = "<>StringReplace[StringTrim[ToString[Reverse[dim]],("{"|"}")...],","->""]<>"
Data.PixelSpacing = "<>StringReplace[StringTrim[ToString[Reverse[vox]],("{"|"}")...],","->""]<>"
Data.NrBits       = "<>StringCases[bit,NumberString]<>"
Data.NrComponents = 1
"
,"Text"]
]


(* ::Subsection:: *)
(*MhdRaw files *)


(* ::Subsubsection::Closed:: *)
(*ExportMhdRaw*)


SyntaxInformation[ExportMhdRaw] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

ExportMhdRaw[file_String, dat_?ArrayQ, voxsize_List, OptionsPattern[]] := 
 Module[{data, type1, type2}, {data, {type1, type2}} = 
   If[Element[Flatten[dat], Integers],
    {dat, If[Max[dat] > 32768 || Min[dat] < -32768,
      {"Integer32", "MET_INT"},
      {"Integer16", "MET_SHORT"}]},
    {N[dat], {"Real32", "MET_DOUBLE"}}
    ];
  ExportMhdRawdat[data, file, type1];
  ExportMhdRawFile[file, Dimensions[data], voxsize, type2];
]


(* ::Subsubsection::Closed:: *)
(*ExportMhdRawdat*)


ExportMhdRawdat[data_, file_, bit_] := 
 Module[{strm}, strm = OpenWrite[file <> ".raw", BinaryFormat -> True];
  BinaryWrite[strm, Flatten[data], bit];
  Close[strm];]


(* ::Subsubsection::Closed:: *)
(*ExportMhdRawFile*)


ExportMhdRawFile[file_, dim_, vox_, bit_] := Module[{},
  Export[file <> ".mhd", "
    ObjectType = Image
    NDims = " <> ToString[Length[dim]] <> "
    BinaryData = True
    BinaryDataByteOrderMSB = False
    ElementSpacing =" <> 
    StringReplace[StringTrim[ToString[Reverse[vox]], ("{" | "}") ...],
      "," -> ""] <> "
    DimSize = " <> 
    StringReplace[StringTrim[ToString[Reverse[dim]], ("{" | "}") ...],
      "," -> ""] <> "
    ElementType = " <> bit <> "
    ElementDataFile = " <> Last[StringSplit[file, "\\"]] <> ".raw
    ", "Text"]]


(* ::Subsection:: *)
(*Batch export error and hist plots*)


(* ::Subsubsection::Closed:: *)
(*ExpHistAll*)


Options[ExpHistAll]={ExportFile->"jpg"}

SyntaxInformation[ExpHistAll] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

ExpHistAll[dat_, name_, text_, OptionsPattern[]]:=
Module[{data,file},
	file=OptionValue[ExportFile];
	data=DeleteCases[dat,Null,{3}];
	data=DeleteCases[dat,Null,{4}];
	Export[name<>"-hist-all-"<>text<>"."<>file,
		GraphicsGrid[{{
			Hist[data[[1]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\) [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]","A) \!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\)"],
			Hist[data[[2]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\) [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]","B) \!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\)"],
			Hist[data[[3]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\) [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]","C) \!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\)"]},{
			Hist[data[[4]],{0,3},"MD [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]","D) MD"],
			Hist[data[[5]],{0,.8},"FA [-]","E) FA"]
			}},Frame-> None,ImageSize->900],ImageResolution->300];
			]


(* ::Subsubsection::Closed:: *)
(*ExpHistInd*)


Options[ExpHistInd]={ExportFile->"jpg"}

SyntaxInformation[ExpHistInd] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

ExpHistInd[dat_,name_,text_,OptionsPattern[]]:=
Module[{data,file},
	file=OptionValue[ExportFile];
	data=DeleteCases[dat,Null,{3}];
	data=DeleteCases[dat,Null,{4}];
	
	Export[name<>"-first-"<>text<>"."<>file,
		Hist[data[[1]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\)","\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\) [\!\(\*SuperscriptBox[\"10\",
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"],ImageSize->300,ImageResolution->300
		];
		
	Export[name<>"-second-"<>text<>"."<>file,
		Hist[data[[2]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\)","\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\) [\!\(\*SuperscriptBox[\"10\",
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"],ImageSize->300,ImageResolution->300
		];
				
	Export[name<>"-third-"<>text<>"."<>file,
		Hist[data[[3]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\)","\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\) [\!\(\*SuperscriptBox[\"10\",
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"],ImageSize->300,ImageResolution->300
		];

	Export[name<>"-ADC-"<>text<>"."<>file,
		Hist[data[[4]],{0,3},"MD","MD [\!\(\*SuperscriptBox[\"10\",
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"],ImageSize->300,ImageResolution->300
		];
			
	Export[name<>"-FA-"<>text<>"."<>file,
		Hist[data[[5]],{0,.8},"FA","FA [-]"],ImageSize->300,ImageResolution->300
		];
		]


(* ::Subsubsection::Closed:: *)
(*ExpPlotsInd*)


Options[ExpPlotsInd]={ExportFile->"jpg"}

SyntaxInformation[ExpPlotsInd] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

ExpPlotsInd[dat_,name_,text_,text2_,OptionsPattern[]]:=
Module[{data,file},
	file=OptionValue[ExportFile];
	data=DeleteCases[dat,Null,{3}];
	data=DeleteCases[dat,Null,{4}];
	
	Export[name<>"-first-"<>text<>"."<>file,
		ErrorPlot[data[[1]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\)",{text2,"\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\) [\!\(\*SuperscriptBox[\"10\",
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],ImageSize->300,ImageResolution->300
		];
	
	Export[name<>"-second-"<>text<>"."<>file,
		ErrorPlot[data[[2]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\)",{text2,"\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\) [\!\(\*SuperscriptBox[\"10\", 
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],ImageSize->300,ImageResolution->300
		];
	
	Export[name<>"-third-"<>text<>"."<>file,
		ErrorPlot[data[[3]],{0,3},"\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\)",{text2,"\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\) [\!\(\*SuperscriptBox[\"10\", 
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],ImageSize->300,ImageResolution->300
		];
	
	Export[name<>"-ADC-"<>text<>"."<>file,
		ErrorPlot[data[[4]],{0,3},"MD",{text2,"MD [\!\(\*SuperscriptBox[\"10\", 
		RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],ImageSize->300,ImageResolution->300
		];
	
	Export[name<>"-FA-"<>text<>"."<>file,
		ErrorPlot[data[[5]],{0,.8},"FA",{text2,"FA [-]"}],ImageSize->300,ImageResolution->300
		];
		];


(* ::Subsubsection::Closed:: *)
(*ExpPlotsAll*)


Options[ExpPlotsAll]={ExportFile->"jpg"}

SyntaxInformation[ExpPlotsAll] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

ExpPlotsAll[dat_,name_,text_,text2_,sOptionsPattern[]]:=
Module[{data,file},
	file=OptionValue[ExportFile];
	data=DeleteCases[dat,Null,{3}];
	data=DeleteCases[dat,Null,{4}];
	
	Export[name<>"-"<>text<>"."<>file,
		GraphicsGrid[{{
			ErrorPlot[data[[1]],{0,3},"A) \!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\)",{text2,"\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\) [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],
			ErrorPlot[data[[2]],{0,3},"B) \!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\)",{text2,"\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\) [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],
			ErrorPlot[data[[3]],{0,3},"C) \!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\)",{text2,"\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\) [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}]},{
			ErrorPlot[data[[4]],{0,3},"D) MD",{text2,"MD [\!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"3\"}]]\) \!\(\*SuperscriptBox[\"mm\", \"2\"]\)/s]"}],
			ErrorPlot[data[[5]],{0,.8},"E) FA",{text2,"FA [-]"}]
			}},Frame-> None,ImageSize->900],ImageResolution->300
			];
			]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
