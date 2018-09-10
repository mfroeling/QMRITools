(* ::Package:: *)

(* ::Title:: *)
(*DTITools ElastixTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`ElastixTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


RegisterData::usage =
"RegisterData[data] registers the data series. If data is 3D it performs multiple 2D registration, if data is 4D it performs multipe 3D registration. The target is the first image orvolume in the series.
RegisterData[{data, vox}] registers the data series using the given voxel size.
RegisterData[{data, mask}] registers the data series only using data whithin the mask.
RegisterData[{data, mask, vox}] registers the data series using the given voxel size only using data within the mask.

RegisterData[target, moving] registers the moving data to the target data. target can be 2D or 3D. moving can be the same of one dimension higher than the target.
RegisterData[{target, mask, vox},{moving, mask, vox}] registers the data using the given voxel size only using data within the mask.
RegisterData[{target, vox}, moving] registers the data using the given voxel size.
RegisterData[target, {moving, vox}] registers the data using the given voxel size.
RegisterData[{target, vox}, {moving, vox}] registers the data using the given voxel size.

RegisterData[{target, mask}, moving] registers the data series only using data whithin the mask.
RegisterData[target, {moving, mask}] registers the data series only using data whithin the mask.
RegisterData[{target, mask}, moving] registers the data series only using data whithin the mask.
RegisterData[{target, mask}, {moving, mask}] registers the data series only using data whithin the mask.

RegisterData[target, {moving, mask, vox}] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, mask}, {moving, mask, vox}] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, vox}, {moving, mask, vox}] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, mask, vox}, moving] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, mask, vox}, {moving, mask}] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, mask, vox}, {moving, vox}] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, mask}, {moving, vox}] registers the data series using the given voxel size only using data within the mask.
RegisterData[{target, vox}, {moving, mask}] registers the data series using the given voxel size only using data within the mask.

Output is the registered data with the dimensions of the moving data. 
If OutputTransformation is True it also outputs the translation, rotation scale and skew of all images or volumes."

RegisterCardiacData::usage =
"RegisterCardiacData[data] registers the data using a 2D algorithm. data can be 3D or 4D.  
RegisterCardiacData[{data,vox}] registers the data series using the given voxel size.
RegisterCardiacData[{data,mask}] registers the data series only using data whithin the mask.
RegisterCardiacData[{data,mask,vox}] registers the data series using the given voxel size only using data within the mask.

Output is the registered data."

RegisterDiffusionData::usage =
"RegisterDiffusionData[{dtidata, vox}] registers a diffusion dataset. dtidata should be 4D {slice, diff, x, y}. vox is the voxelsize of the data.
RegisterDiffusionData[{dtidata, dtimask, vox}] registers the data series using the given voxel size only using data within the mask.
RegisterDiffusionData[{dtidata ,vox}, {anatdata, voxa}] registers a diffusion dataset. The diffusion data is also registered to the anatdata.
RegisterDiffusionData[{dtidata, dtimask, vox}, {anatdata, voxa}] registers the data series using the given voxel size only using data within the mask.
RegisterDiffusionData[{dtidata,vox}, {anatdata, anatmask, voxa}] registers the data series using the given voxel size only using data within the mask.
RegisterDiffusionData[{dtidata, dtimask, vox}, {anatdata, anatmask, voxa}] registers the data series using the given voxel size only using data within the mask.

Output is the registered dtidata and, if anatdata is given, the registered dtidata in anatomical space. If OutputTransformation is True it also outputs the translation, rotation scale and skew of all images or volumes."

RegisterDiffusionDataSplit::usage = 
"RegisterDiffusionDataSplit[dtidata, vox] is identical to Register diffusion data however left and right side of the data are registered seperately.
RegisterDiffusionDataSplit[{dtidata, vox}, {anatdata, voxa}] is identical to Register diffusion data however left and right side of the data are registered seperately.
RegisterDiffusionDataSplit[{dtidata, dtimask, vox}, {anatdata, anatmask, voxa}] is identical to Register diffusion data however left and right side of the data are registered seperately.

Splitting the data is done using the function CutData and merged wit Stich data.
Output is the registered data."

RegisterDataSplit::usage = 
"RegisterDataSplit[target, moving] is identical to RegisterData data however left and right side of the data are registered seperately.

Splitting the data is done using the function CutData and merged wit Stich data.
Output is the registered data."

ReadTransformParameters::usage = 
"ReadTransformParameters[directory] reads the tranfomation parameters generated by RegisterData. The directory should be the TempDirectory were the registration is stored. DeleteTempDirectory should be False.

Output is the affine transformation vector per volume."

CorrectGradients::usage = 
"CorrectGradients[grad, transformation] corrects the gradient directions grad with the tranformation parameters from RegisterData or RegisterDiffusionData.

Output is the corrected gradient vector."

CorrectBmatrix::usage = 
"CorrectBmatrix[bmat, transformation] corrects the bmatrix bmat with the tranformation parameters from RegisterData or RegisterDiffusionData.

Output is the corrected bmatrix."

TransformData::usage = 
"TransformData[{data,vox}] deforms the data according to the last output of register data.
The directory should be the TempDirectory were the registration is stored. DeleteTempDirectory should be False."

RegisterDataTransform::usage = 
"RegisterDataTransform[target, moving, {moving2nd, vox}] performs the registration exactly as RegisterData. target and moving are the inputs for Registerdata, which can be {data,mask,vox}.
After the registeration is done the moving2nd data is deformed acording to the output of the registrtion of moving.

moving2nd can have the same dimensions of moving or one dimension higher (e.g. 3D and 3D or 3D and 4D). 

Output is {registered moving, deformed moving2nd}."

RegisterDataTransformSplit::usage = 
"RegisterDataTransformSplit[target, moving, {moving2nd, vox}] is idenditcal to RegisterDataTransform with the same functionality as RegisterDataSplit.
This means the data is split in two using the function CutData and merged wit Stich data.

Output is {registered moving, deformed moving2nd}."


(* ::Subsection::Closed:: *)
(*Options*)


Iterations::usage =
"Iterations is an options for RegisterData, RegisterDiffusionData, and RegisterDataTransform. 
It specifies the number of iterations used by the registration functions."

Resolutions::usage =
"Resolutions is an options for RegisterData, RegisterDiffusionData, and RegisterDataTransform. 
It specifies the number of scale space resolutions used by the registration functions."

HistogramBins::usage =
"HistogramBins is an options for RegisterData, RegisterDiffusionData, and RegisterDataTransform. 
It specifies the number of bins of the joined histogram used by the registration functions."

NumberSamples::usage =
"NumberSamples is an options for RegisterData, RegisterDiffusionData, and RegisterDataTransform. 
It specifies the number of random samples that are taken each iteration used by the registration functions."

OutputImage::usage =
"OutputImage is an options for RegisterData, RegisterDiffusionData, and RegisterDataTransform. 
It specifies if the result image should be writen in the TempDirectory as nii file."

InterpolationOrderReg::usage =
"InterpolationOrderReg is an options for RegisterData, RegisterDiffusionData, and RegisterDataTransform. 
It specifies the interpolation order used in the registration functions."

MethodReg::usage = 
"MethodReg is an options for RegisterData, RegisterDiffusionData, RegisterCardiacData and RegisterDataTransform. 
It spefifies which registration method to use. Mehtods can be be \"rigid\",\"affine\", \"bspline\" or \"cyclyc\"."

BsplineSpacing::usage =
"BsplineSpacing is an options for RegisterData, RegisterDiffusionData, RegisterCardiacData and RegisterDataTransform. 
It specifies the spacing of the bsplines if the method is \"bspline\"."

TempDirectory::usage = 
"TempDirectory is an options for RegisterData, RegisterDiffusionData, RegisterCardiacData and RegisterDataTransform. 
It specifies the temprary directory used to perform and output the registration."

DeleteTempDirectory::usage =
"DeleteTempDirectory an options for RegisterData, RegisterDiffusionData, RegisterCardiacData and RegisterDataTransform. 
It specifies if the temp directory should be deleted after the registration is finisched."

PrintTempDirectory::usage = 
"PrintTempDirectory is an options for RegisterData, RegisterDiffusionData, RegisterCardiacData and RegisterDataTransform. 
It spefifies if the location of the temp directory should be deplayed."

RegistrationTarget::usage = 
"RegistrationTarget is an option for RegisterDiffusionData and RegisterCardiacData. Specifies which target to uses for registration if using \"rigid\", \"affine\" or \"bspline\" as MethodReg.
If the MethodReg is \"cyclyc\" or \"PCA\" it does not need a target and this options does nothing. 
Values can be \"First\", \"Mean\" or \"Median\"."

BsplineDirections::usage = 
"BsplineDirections is an option for RegisterData ad RegisterDiffusionData. 
It gives the direction in which the bsplines are allowed to move when registering diffusion data to anatomical space."

AffineDirections::usage = 
"AffineDirections is an option for RegisterData ad RegisterDiffusionData. 
It gives the directions in which data can be moved when registering diffusion data to anatomical space."

OutputTransformation::usage =
"OutputTransformation is an option for RegisterData ad RegisterDiffusionData.
It specifies if the tranformation paramters (translation, rotation, scale and skew) should be given as output in the registration functions."

IterationsA::usage = 
"IterationsA is an option for RegisterDiffusionData. 
It specifies the number of iterations used when registering diffusion data to anatomical space."

ResolutionsA::usage =
"ResolutionsA is an option for RegisterDiffusionData.
It specifies the number of scale space resolutions used when registering diffusion data to anatomical space."

HistogramBinsA::usage =
"HistogramBinsA is an option for RegisterDiffusionData.
It specifies the number of bins of the joined histogram used when registering diffusion data to anatomical space."

NumberSamplesA::usage =
"NumberSamplesA is an option for RegisterDiffusionData.
It specifies the number of random samples that are taken each iteration when registering diffusion data to anatomical space."

InterpolationOrderRegA::usage =
"InterpolationOrderRegA is an option for RegisterDiffusionData.
It specifies the interpolation order used in the registration functions when registering diffusion data to anatomical space."
 
MethodRegA::usage =
"MethodRegA is an option for RegisterDiffusionData.
It spefifies which registration method to use when registering diffusion data to anatomical space. Mehtods can be be \"rigid\",\"affine\" or \"bspline\"."

UseGPU::usage = 
"UseGPU is an option for RegisterData. The value is {bool, gpu} where bool is True or False, and gpu is the gpu ID which is an integer or Automatic."

PCAComponents::usage = 
"PCAComponents is an option for RegisterData. It speciefies how many PCA components are used if method is set to \"PCA\""

FindTransform::usage = 
"FindTransform is an option for TransformData and RegisterTransformData. It specifies where to find the transformfile."

SplitMethod::usage = 
"SplitMethod is an option for RegisterDataSplit and RegisterDataTransformSplit. values can be \"mean\", \"moving\", \"target\""


(* ::Subsection::Closed:: *)
(*Error Messages*)


RegisterData::vol="The `1`D datasets should have 2 or more volumes, it has `2` volumes."

RegisterData::dim="Datasets should both be 2D or 3D, or 2D and 3D, or 3D and 4D, current datasets are `1`D and `2`D."

RegisterData::dims="Dataset should be 3D or 4D, current dataset is `1`D."

RegisterData::vox="voxel size should be {z,x,y} and numeric, current sizes are `1` and `2`."

RegisterData::voxs="voxel size should be {z,x,y} and numeric, current size is `1`."

RegisterData::met="MethodReg should be \"rigid\",\"affine\", \"bspline\" or \"cyclyc\", current method is `1`."

RegisterData::metc="If the MethodReg is \"cyclyc\" no target can be given."

RegisterData::mask="The mask dimensions `1` should be equal to the data dimensions `2`."

RegisterData::dir="Temporary directory not created."

RegisterData::elastix="Elastix not found, check if DTITools is installed in the $BaseDirectory or $UserBaseDirectory."

RegisterData::par="`1` should be a number or a list of numbers with length `2`."

RegisterData::fatal="Fatal error encountered."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Support Functions*)


(* ::Subsubsection::Closed:: *)
(*ParString*)


SchedulePar[res_]:=ListToString[ToString[2^#1]<>" "<>ToString[2^#]<>" 0" &/@ Reverse[Range[res]-1]];

ListToString[list_]:=StringJoin[Riffle[ToString/@list," "]]

ParString[{itterations_,resolutions_,bins_,samples_,intOrder_},{type_,output_},{grid_, derscB_, derscA_,pca_},{openCL_,gpu_}]:="// *********************
// * "<>type<>"
// *********************

// *********************
// * ImageTypes
// *********************
(FixedInternalImagePixelType \"float\")
(MovingInternalImagePixelType \"float\")
(UseDirectionCosines \"true\")

// *********************
// * Components
// *********************
(FixedImagePyramid \"FixedSmoothingImagePyramid\")
(MovingImagePyramid \"MovingSmoothingImagePyramid\")
(Registration \"MultiResolutionRegistration\")
"<>Switch[type,
"cyclyc",
"(Interpolator \"ReducedDimensionBSplineInterpolator\")
(ResampleInterpolator \"FinalReducedDimensionBSplineInterpolator\")
(Metric \"VarianceOverLastDimensionMetric\")",
"PCA",
"(Interpolator \"ReducedDimensionBSplineInterpolator\")
(ResampleInterpolator \"FinalReducedDimensionBSplineInterpolator\")
(Metric \"PCAMetric2\")
(NumEigenValues "<>ToString[pca]<>")",
_,
"(Interpolator \"BSplineInterpolator\")
(ResampleInterpolator \"FinalBSplineInterpolator\")
(Metric \"AdvancedMattesMutualInformation\")"
]<>"
(BSplineInterpolationOrder "<>ToString[intOrder]<>")
"<>If[openCL,
"(OpenCLResamplerUseOpenCL \"true\")
(OpenCLDeviceID \""<>ToString[gpu]<>"\")
(Resampler \"OpenCLResampler\")"
,
"(Resampler \"DefaultResampler\")"
]<>"
(Optimizer \"AdaptiveStochasticGradientDescent\")
"<>Switch[type,
"translation",
"(Transform \"TranslationTransform\")
(MovingImageDerivativeScales "<>ToString[Clip[derscA[[3]]]] <> " " <>ToString[Clip[derscA[[2]]]] <> " " <> ToString[Clip[derscA[[1]]]]<>")",
"rigid",
"(Transform \"EulerTransform\")",
"rigidDTI",
"(Transform \"AffineDTITransform\")",
"affine",
"(Transform \"AffineTransform\")",
"affineDTI",
"(Transform \"AffineDTITransform\")
(MovingImageDerivativeScales "<>ToString[Clip[derscA[[3]]]] <> " " <>ToString[Clip[derscA[[2]]]] <> " " <> ToString[Clip[derscA[[1]]]]<>")",
"bspline",
"(Transform \"RecursiveBSplineTransform\")
(FinalGridSpacingInPhysicalUnits "<>ToString[grid[[3]]]<>" "<>ToString[grid[[2]]]<>" "<>ToString[grid[[1]]]<>")
(MovingImageDerivativeScales "<>ToString[Clip[derscB[[3]]]] <> " " <>ToString[Clip[derscB[[2]]]] <> " " <> ToString[Clip[derscB[[1]]]]<>")",
"PCA",
"(Transform \"BSplineStackTransform\")
(FinalGridSpacingInPhysicalUnits "<>ToString[grid[[3]]]<>" "<>ToString[grid[[2]]]<>" 0)

// *********************
// * Metric settings
// *********************
(MovingImageDerivativeScales 1.0 1.0 0.0)
(NumberOfSamplesForExactGradient 50000)
(SubtractMean \"true\")",
"cyclyc",
"(Transform \"BSplineStackTransform\")
(FinalGridSpacingInPhysicalUnits "<>ToString[grid[[3]]]<>" "<>ToString[grid[[2]]]<>" 0)

// *********************
// * Metric settings
// *********************
(MovingImageDerivativeScales 1.0 1.0 0.0)
(SubtractMean \"true\")"
]<>"

// *********************
// * Mask settings
// *********************
(ErodeMask \"false\")
(ErodeFixedMask \"false\")

// *********************
// * Optimizer settings
// *********************
(NumberOfResolutions "<>ToString[resolutions]<>")
(MaximumNumberOfIterations "<>ToString[itterations]<>")
(ASGDParameterEstimationMethod \"Original\")
(AutomaticParameterEstimation \"true\")
"<>Switch[type,
"PCA",
"",
"cyclyc",
"",
_,
"(AutomaticTransformInitialization \"true\")
(AutomaticScalesEstimation \"true\")"
]<>"

// *********************
// * Transform settings
// *********************
(HowToCombineTransforms \"Compose\")

// *********************
// * Pyramid settings
// *********************
(NumberOfHistogramBins "<>ToString[bins]<>")
"<>Switch[type,
"affineDTI",""(*"(Scales -1.000000e+00 -1.000000e+00 -1.000000e+00  1.000000e+06  1.000000e+06  1.000000e+06  1.000000e+06  1.000000e+06  1.000000e+06 -1.000000e+00 -1.000000e+00 -1.000000e+00)"*),
"rigidDTI","(Scales -1.000000e+00 -1.000000e+00 -1.000000e+00  3.000000e+38  3.000000e+38  3.000000e+38  3.000000e+38  3.000000e+38  3.000000e+38 -1.000000e+00 -1.000000e+00 -1.000000e+00)",
"PCA",
"(ImagePyramidSchedule "<>SchedulePar[resolutions]<>")",
"cyclyc",
"(ImagePyramidSchedule "<>SchedulePar[resolutions]<>")",
_,""
]<>"

// *********************
// * Sampler parameters
// *********************
(NumberOfSpatialSamples "<>ToString[samples]<>")
"<>Switch[type,
"cyclyc",
"(ImageSampler \"Random\")",
_,
"(ImageSampler \"RandomCoordinate\")"
]<>"
(CheckNumberOfSamples \"false\")
(NewSamplesEveryIteration \"true\")
(MaximumNumberOfSamplingAttempts 5)
(FinalBSplineInterpolationOrder "<>ToString[intOrder]<>")

// *********************
// * Output settings
// *********************
(DefaultPixelValue 0)
(WriteTransformParametersEachIteration \"false\")
(WriteResultImage  \""<>output<>"\")
(ResultImageFormat \"nii.gz\")
(ResultImagePixelType \"float\")
"


(* ::Subsubsection::Closed:: *)
(*FindElastix*)


FindElastix[]:=Module[{fil1,fil2},
	Switch[$OperatingSystem,
		"Windows",
		fil1=$UserBaseDirectory<>"\\Applications\\DTITools\\Applications\\elastix.exe";
		fil2=$BaseDirectory<>"\\Applications\\DTITools\\Applications\\elastix.exe";
		,
		"MacOSX",
		fil1=$UserBaseDirectory<>"/Applications/DTITools/Applications/bin/elastix";
		fil2=$BaseDirectory<>"/Applications/DTITools/Applications/bin/elastix";
		];
	If[FileExistsQ[fil1],fil1,If[FileExistsQ[fil2],fil2,"error"]]
]


(* ::Subsubsection::Closed:: *)
(*FindTransformix*)


FindTransformix[]:=Module[{fil1,fil2},
	Switch[$OperatingSystem,
		"Windows",
		fil1=$UserBaseDirectory<>"\\Applications\\DTITools\\Applications\\transformix.exe";
		fil2=$BaseDirectory<>"\\Applications\\DTITools\\Applications\\transformix.exe";
		,
		"MacOSX",
		fil1=$UserBaseDirectory<>"/Applications/DTITools/Applications/bin/transformix";
		fil2=$BaseDirectory<>"Applications/DTITools/Applications/bin/transformix";
	];
	If[FileExistsQ[fil1],fil1,If[FileExistsQ[fil2],fil2,"error"]]
]


(* ::Subsubsection::Closed:: *)
(*RunElastix*)


RunElastix[elastix_,tempdir_,parfile_,{inpfol_,movfol_,outfol_},{fixedi_,movingi_,out_},{maskfi_,maskmi_}]:=Module[
	{fixed,moving, maskf, maskm, command,inpfold,outfold,movfold,parfiles,copy,maskfFile,maskmFile,elastixFol},
	
	(*make files into gz where needed*)
	fixed = If[FileExtension[fixedi] == "nii", fixedi<>".gz", fixedi];
	moving = If[FileExtension[movingi] == "nii", movingi<>".gz", fixedi];
	maskf = If[FileExtension[maskfi] == "nii", maskfi<>".gz", maskfi];
	maskm = If[FileExtension[maskmi] == "nii", maskmi<>".gz", maskmi];
	
	(*make elastix command based on operating system*)
	Switch[$OperatingSystem,
		"Windows",
		inpfold=If[inpfol=="",tempdir,tempdir<>inpfol<>"\\"];
		movfold=If[movfol=="",tempdir,tempdir<>movfol<>"\\"];
		outfold=If[outfol=="",StringDrop[tempdir,-1],tempdir<>outfol];
		
		maskfFile=If[maskf==="",""," -fMask \""<>tempdir<>maskf<>"\""];
		maskmFile=If[maskm==="",""," -mMask \""<>tempdir<>maskm<>"\""];
		
		parfiles=StringJoin[" -p \""<>tempdir<>#<>"\""&/@parfile];
		copy=If[out=="","","copy \""<>tempdir<>outfol<>"\\result."<>ToString[Length[parfile]-1]<>".nii.gz\" \""<>tempdir<>outfol<>"\\"<>out<>"\""];
		
		command=
		"\""<>elastix<>
		"\" -f \""<>inpfold<>fixed<>
		"\" -m \""<>movfold<>moving<>
		"\" -out \""<>outfold<>"\""<>
		maskfFile<>
		maskmFile<>
		(*"\" -p \""<>tempdir<>parfile<>*)
		parfiles<>" > "<>movfold<>"\"output.txt\" \n"<>
		copy<>" \n"<>
		"exit \n";

		,
		"MacOSX",	
		inpfold=If[inpfol=="",tempdir,tempdir<>inpfol<>"/"];
		movfold=If[movfol=="",tempdir,tempdir<>movfol<>"/"];
		outfold=If[outfol=="",StringDrop[tempdir,-1],tempdir<>outfol];
		
		maskfFile=If[maskf==="",""," -fMask "<>tempdir<>maskf<>" "];
		maskmFile=If[maskm==="",""," -mMask "<>tempdir<>maskm<>" "];
		
		(*parfiles=StringJoin["-p "<>tempdir<>#<>""&/@parfile];*)
		parfiles=StringJoin[" -p "<>tempdir<>#&/@parfile];
		copy=If[out=="",""," cp "<>tempdir<>outfol<>"result."<>ToString[Length[parfile]-1]<>".nii "<>tempdir<>outfol<>""<>out];
		
		elastixFol=StringDrop[DirectoryName[elastix, 2], -1];
		
		command=
		"export PATH="<>elastixFol<>"/bin:$PATH 
		export DYLD_LIBRARY_PATH="<>elastixFol<>"/lib:$DYLD_LIBRARY_PATH \n"<>
		elastix<>
		" -f "<>inpfold<>fixed<>
		" -m "<>movfold<>moving<>
		" -out "<>outfold<>
		maskfFile<>
		maskmFile<>
		parfiles<>" > "<>movfold<>"output.txt"<>
		"\n"<>copy<>" \n"<>
		"exit";
	];
	(*Print[command];*)
	(*perform elastix on system shell*)	
	RunProcess[$SystemShell,"StandardOutput",command];
]


(* ::Subsubsection::Closed:: *)
(*ElastixCommand*)


ElastixCommand[elastix_,tempdir_,parfile_,{inpfol_,movfol_,outfol_},{fixed_,moving_,out_},{maskf_,maskm_}]:=Block[
	{command,inpfold,outfold,movfold,maskfFile,maskmFile,outfile,parfiles,copy, elastixFol},
	
	(*make elastix command based on operating system*)
	Switch[$OperatingSystem,
		"Windows",
		inpfold=If[inpfol=="",tempdir,tempdir<>inpfol<>"\\"];
		movfold=If[movfol=="",tempdir,tempdir<>movfol<>"\\"];
		outfold=If[outfol=="",StringDrop[tempdir,-1],tempdir<>outfol];
		outfile=outfold<>"\\"<>out;
		
		maskfFile=If[maskf==="",""," -fMask \""<>tempdir<>maskf<>"\""];
		maskmFile=If[maskm==="",""," -mMask \""<>tempdir<>maskm<>"\""];
		
		parfiles=StringJoin[" -p \""<>tempdir<>#<>"\""&/@parfile];
		copy=If[out=="","","@ copy \""<>outfold<>"\\result."<>ToString[Length[parfile]-1]<>".nii.gz\" \""<>outfile<>"\""];
		
		command=
		"@ \""<>elastix<>
		"\" -f \""<>inpfold<>fixed<>
		"\" -m \""<>movfold<>moving<>
		"\" -out \""<>outfold<>"\""<>
		maskfFile<>
		maskmFile<>
		parfiles<>" > \""<>movfold<>"output.txt\" \n"<>
		copy<>" \n";
		,
		"MacOSX",
		inpfold=If[inpfol=="",tempdir,tempdir<>inpfol<>"/"];
		movfold=If[movfol=="",tempdir,tempdir<>movfol<>"/"];
		outfold=If[outfol=="",StringDrop[tempdir,-1],tempdir<>outfol];
		outfile=outfold<>"/"<>out;
		
		maskfFile=If[maskf==="",""," -fMask "<>tempdir<>maskf<>" "];
		maskmFile=If[maskm==="",""," -mMask "<>tempdir<>maskm<>" "];
		
		parfiles=StringJoin[" -p "<>tempdir<>#&/@parfile];
		copy=If[out=="",""," cp "<>outfold<>"/result."<>ToString[Length[parfile]-1]<>".nii "<>outfile];
		
		elastixFol = StringDrop[DirectoryName[elastix, 2], -1];
		
		command=
		"export PATH="<>elastixFol<>"/bin:$PATH 
		export DYLD_LIBRARY_PATH="<>elastixFol<>"/lib:$DYLD_LIBRARY_PATH \n"<>
		elastix<>
		" -f "<>inpfold<>fixed<>
		" -m "<>movfold<>moving<>
		" -out "<>outfold<>" "<>
		maskfFile<>
		maskmFile<>
		parfiles<>" > "<>movfold<>"output.txt \n"<>
		copy<>" \n";	
	];
	
	{command,outfile}
]


(* ::Subsubsection::Closed:: *)
(*RunBatfile*)


RunBatfile[tempdir_,command_]:=Block[{batfile,com},
	(*make elastix sh/bat based on operating system*)
	Switch[$OperatingSystem,
		"Windows",
		batfile = tempdir<>"elastix-batch.bat";
		Export[batfile,StringJoin[command],"TEXT"];
		com = "\"" <> batfile <> "\"\n exit \n";
		,
		"MacOSX",
		batfile = tempdir<>"/elastix-bash.sh";
		Export[batfile,StringJoin[command],"TEXT"];
		com= "chmod 700 "<>batfile<>"\n"<>batfile<> "\n exit \n";
	];
	
	(*perform sh/bat on system shell*)	
	RunProcess[$SystemShell, "StandardOutput", com];
]


(* ::Subsubsection::Closed:: *)
(*StringPad*)


StringPad[x_] := 
 StringJoin[
  PadLeft[{ToString[x]}, 5 - StringLength[ToString[x]], "0"]
]


(* ::Subsubsection::Closed:: *)
(*ConcatenateTransformFiles*)


ConcatenateTransformFiles[files_, outDir_] := Block[{len, filesi, tfile,slash},
  (*import the transform files*)
  len = Range[Length[files]];
  filesi = Import[#, "Lines"] & /@ files;
  slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
  
  (*concatenate the transform files*)
  (
  	tfile = If[# == 1, "NoInitialTransform", outDir <> slash <> "FinalTransform." <> ToString[# - 2] <> ".txt"];
  	filesi[[#, 4]] = "(InitialTransformParametersFileName \"" <> tfile <> "\")";
  	Export[outDir <> slash <> "FinalTransform." <> ToString[# - 1] <> ".txt", filesi[[#]]];
  ) & /@ len;
  ]


(* ::Subsubsection::Closed:: *)
(*RunBatfileT*)


RunBatfileT[tempdir_, command_] := Block[{batfile, com},
	Switch[$OperatingSystem,
		"Windows",
		batfile = tempdir <> "\\transformix-batch.bat";
		Export[batfile, StringJoin[command], "TEXT"];
		com = batfile <> "\n exit \n";
		,
		"MacOSX",
		batfile = tempdir <> "/transformix-bash.sh";
		Export[batfile, StringJoin[command], "TEXT"];
		com = "chmod 700 "<>batfile<>"\n"<>batfile<> "\n exit \n";
	];
	
	RunProcess[$SystemShell, "StandardOutput", com];
]


(* ::Subsubsection::Closed:: *)
(*TransformixCommand*)


TransformixCommand[tempDir_] := Block[{volDirs, transformix, transFol},
  transformix = FindTransformix[];
  transFol = StringDrop[DirectoryName[transformix, 2], -1];
  
  volDirs = FileNames["vol*", tempDir, 1];
  
  Switch[$OperatingSystem,
  	"Windows",
  	(
  		"@ \"" <> transformix <>
  		"\" -in \"" <> First[FileNames["moving*", #]] <>
  		"\" -out \"" <> # <>
  		"\" -tp \"" <> Last[FileNames["FinalTransform*", #]] <> "\"" <>
  		" > \"" <> # <> "\\outputa.txt\" \n" <>
  		"@ rename \"" <> # <> "\\result.nii.gz\" resultA-3D.nii.gz \n"
  	) & /@ volDirs
  	,
  	"MacOSX",
  	(
  		"export PATH="<>transFol<>"/bin:$PATH 
		export DYLD_LIBRARY_PATH="<>transFol<>"/lib:$DYLD_LIBRARY_PATH \n"<>
		transformix <>
		" -in " <> First[FileNames["moving*", #]] <>
  		" -out " <> # <>
  		" -tp " <> Last[FileNames["FinalTransform*", #]] <> 
	    " > " <> # <> "/outputa.txt \n" <>
  		" mv " <> # <> "/result.nii.gz "<> # <> "/resultA-3D.nii.gz \n"
  	) & /@ volDirs
  ]
]


(* ::Subsection:: *)
(*RegisterData/Split*)


(* ::Subsubsection::Closed:: *)
(*RegisterData*)


Options[RegisterData]={
Iterations->1000,
Resolutions->1,
HistogramBins->64,
NumberSamples->2000,
InterpolationOrderReg->3,
BsplineSpacing->30,
BsplineDirections->{1,1,1},
AffineDirections->{1,1,1},
MethodReg->"affine",
OutputImage->True,
TempDirectory->"Default",
DeleteTempDirectory->True,
PrintTempDirectory->True,
OutputTransformation->False,
UseGPU->{False,Automatic},
PCAComponents->3
};

SyntaxInformation[RegisterData]={"ArgumentsPattern"->{_,_.,OptionsPattern[]}};


(* ::Subsubsection::Closed:: *)
(*RegisterData Series*)


(*series have no target defninition*)
(*register series of data sets, no vox definition, no mask definition*)
RegisterData[series_?ArrayQ,opts:OptionsPattern[]]:=RegisterData[{series,{1,1,1}},opts]
(*register series of data sets, vox definition, no mask definition*)
RegisterData[{series_?ArrayQ,vox:{_?NumberQ,_?NumberQ,_?NumberQ}},opts:OptionsPattern[]]:=RegisterData[{series,{1},vox},opts]
(*register series of data sets, no vox definition, mask definition*)
RegisterData[{series_?ArrayQ,mask_?ArrayQ},opts:OptionsPattern[]]:=RegisterData[{series,mask,{1,1,1}},opts]

(*register series of data sets, vox definition, mask definition, no target definition*)
RegisterData[{series_?ArrayQ,mask_?ArrayQ,vox:{_?NumberQ,_?NumberQ,_?NumberQ}},opts:OptionsPattern[]]:=Module[
{depthS,error,dim,dimm,dimL,target,moving,dataout,
voxL,output,cyclyc,maskf,maskm},

(*set error*)
error=False;

(*get data properties*)
depthS=ArrayDepth[series];
dim=Dimensions[series];
dimL=If[depthS==3,dim[[1]],dim[[2]]];
dimm=Dimensions[mask];
voxL=Length[vox];

cyclyc=OptionValue[MethodReg]==="cyclyc"||OptionValue[MethodReg]==="PCA";

(*check dimensions*)
(*series must be 3 of 4D*)
If[!(depthS==3||depthS==4),Message[RegisterData::dims,depthS];Return[Message[DiffusionReg::fatal]]];
(*sereis must have 2 or more volumes*)
If[!dimL>=2,Message[RegisterData::vol,depthS,dim[[1]]];Return[Message[DiffusionReg::fatal]]];

(*check voxel sizes*)
If[voxL!=3||!(NumberQ@Total@vox),Message[RegisterData::voxs,vox];Return[Message[DiffusionReg::fatal]]];

(*check mask*)
If[mask!={1},If[cyclyc,
(*cyclyc mask needs to be same dimensions as moving data*)
If[dim!=dimm,Message[RegisterData::mask,dimm,dim];Return[Message[DiffusionReg::fatal]]],
If[depthS==3,
(*normal mask, one mask for all or one mask per volume*)
If[!(dim[[2;;3]]==dimm||dim==dimm),Message[RegisterData::mask,dimm,dim];Return[Message[DiffusionReg::fatal]]],
If[!(dim[[{1,3,4}]]==dimm||dim==dimm),Message[RegisterData::mask,dimm,dim];Return[Message[DiffusionReg::fatal]]]
]
]];

(*check if method is cyclyc*)
If[cyclyc,
(*cyclyc series*)
(*define moving and target voluems*)
target=moving=series;

(*go to registration function*)
output=RegisterDatai[{target,mask,vox},{moving,mask,vox},OptionValue[MethodReg],opts];
output
,
(*normal series*)
(*define moving and target voluems*)
{target,moving}=If[depthS==3,
{series[[1]],series[[2;;]]},
{series[[All,1]],Transpose@series[[All,2;;]]}
];

{maskf,maskm}=If[dimm==dim,If[depthS==3,{mask[[1]],mask[[2;;]]},{mask[[All,1]],Transpose@mask[[All,2;;]]}],{mask,mask}];

(*go to registration function*)
output=RegisterDatai[{target,maskf,vox},{moving,maskm,vox},"series",opts];

If[OptionValue[OutputTransformation],
		(*output data with tranformation parameters*)
		dataout=Prepend[output[[1]],target];
		{If[depthS==4,Transpose@dataout,dataout],Prepend[output[[2]],{0,0,0,0,0,0,1,1,1,0,0,0}]}
		,
		(*output dat without transformation parameters*)
		dataout=Prepend[output,target];
		If[depthS==4,Transpose@dataout,dataout]
]
]
]


(* ::Subsubsection::Closed:: *)
(*RegisterData Volumes*)


(*Volumes do have target defninition*)

(*register two data sets, vox definition and mask definition*)
RegisterData[
{target_?ArrayQ,maskt:{_?ListQ..},voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
{moving_?ArrayQ,maskm:{_?ListQ..}}
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,voxt},{moving,maskm,{1,1,1}},opts];
RegisterData[
{target_?ArrayQ,maskt:{_?ListQ..}},
{moving_?ArrayQ,maskm:{_?ListQ..},voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,{1,1,1}},{moving,maskm,voxm},opts];

(*register two data sets, vox definition and maskt definition*)
RegisterData[
{target_?ArrayQ,maskt:{_?ListQ..},voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
moving_?ArrayQ
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,voxt},{moving,{1},{1,1,1}},opts];
RegisterData[
{target_?ArrayQ,maskt:{_?ListQ..}},
{moving_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,{1,1,1}},{moving,{1},voxm},opts];
RegisterData[
{target_?ArrayQ,maskt:{_?ListQ..},voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
{moving_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,voxt},{moving,{1},voxm},opts];

(*register two data sets, vox definition and maskm definition*)
RegisterData[
{target_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
{moving_?ArrayQ,maskm:{_?ListQ..}}
,opts:OptionsPattern[]]:=RegisterData[{target,{1},voxt},{moving,maskm,{1,1,1}},opts];
RegisterData[
target_?ArrayQ,
{moving_?ArrayQ,maskm:{_?ListQ..},voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,{1},{1,1,1}},{moving,maskm,voxm},opts];
RegisterData[
{target_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
{moving_?ArrayQ,maskm:{_?ListQ..},voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,{1},voxt},{moving,maskm,voxm},opts];

(*register two data sets, no vox definition and mask definition*)
RegisterData[
{target_?ArrayQ,maskt_:{_?ListQ..}},
moving_?ArrayQ
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,{1,1,1}},{moving,{1},{1,1,1}},opts];
RegisterData[
target_?ArrayQ,
{moving_?ArrayQ,maskm:{_?ListQ..}}
,opts:OptionsPattern[]]:=RegisterData[{target,{1},{1,1,1}},{moving,maskm,{1,1,1}},opts];
RegisterData[
{target_?ArrayQ,maskt:{_?ListQ..}},
{moving_?ArrayQ,maskm:{_?ListQ..}}
,opts:OptionsPattern[]]:=RegisterData[{target,maskt,{1,1,1}},{moving,maskm,{1,1,1}},opts];

(*register two data sets, vox definition and no mask definition*)
RegisterData[
{target_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
moving_?ArrayQ
,opts:OptionsPattern[]]:=RegisterData[{target,{1},voxt},{moving,{1},{1,1,1}},opts]
RegisterData[
target_?ArrayQ,
{moving_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,{1},{1,1,1}},{moving,{1},voxm},opts]
RegisterData[
{target_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
{moving_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}}
,opts:OptionsPattern[]]:=RegisterData[{target,{1},voxt},{moving,{1},voxm},opts]

(*register two data sets, no vox definition and no mask definition*)
RegisterData[
target_?ArrayQ,
moving_?ArrayQ
,opts:OptionsPattern[]]:=RegisterData[{target,{1},{1,1,1}},{moving,{1},{1,1,1}},opts]

(*register two data sets, mask and vox definition*)
RegisterData[
{target_?ArrayQ,maskt_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},
{moving_?ArrayQ,maskm_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}},opts:OptionsPattern[]]:=Module[
{depthT,depthM,voxtL,voxmL,error,dim,type,mov,output},

(*set error*)
error=False;

(*Check Method, cyclyc and PCA only possible for series*)
If[OptionValue[MethodReg]==="cyclyc"||OptionValue[MethodReg]==="PCA",error=True;Message[RegisterData::metc];];

(*get data properties*)
depthT=ArrayDepth[target];
depthM=ArrayDepth[moving];
dim=Dimensions[moving];
voxtL=Length[voxt];
voxmL=Length[voxm];

(*check dimensions and determine type*)
(*2D-2D, 3D-3D*)
type=If[depthT==depthM,
"vol",
(*2D-3D, 3D-4D*)
If[(depthT==2||depthT==3)&&depthM==depthT+1,
"series",
(*error*)
error=True;Message[RegisterData::dim,depthT,depthM];
]
];

(*check voxel sies*)
If[voxtL!=3||voxmL!=3||!(NumberQ@Total@voxt)||!(NumberQ@Total@voxm),Message[RegisterData::vox,voxt,voxm];Return[Message[DiffusionReg::fatal]]];

(*if error found quit*)
If[error,Return[Message[DiffusionReg::fatal]]];

(*define moving voluems*)
mov=If[depthM==4,Transpose@moving,moving];

(*go to registration function*)
output=RegisterDatai[{target,maskt,voxt},{mov,maskm,voxm},type,opts];

If[OptionValue[OutputTransformation],
		{If[depthM==4,Transpose@output[[1]],output[[1]]],Prepend[output[[2]],{0,0,0,0,0,0,1,1,1,0,0,0}]},
		If[depthM==4,Transpose@output,output]
]
]


(* ::Subsubsection::Closed:: *)
(*RegisterDatai*)


Options[RegisterDatai]=Options[RegisterData];

RegisterDatai[{target_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},{moving_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}},type_,opts:OptionsPattern[]]:=
RegisterDatai[{target,{1},voxt},{moving,{1},voxm},type,opts]

RegisterDatai[{target_?ArrayQ,maskt_,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},{moving_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}},type_,opts:OptionsPattern[]]:=
RegisterDatai[{target,maskt,voxt},{moving,{1},voxm},type,opts]

RegisterDatai[{target_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},{moving_?ArrayQ,maskm_,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}},type_,opts:OptionsPattern[]]:=
RegisterDatai[{target,{1},voxt},{moving,maskm,voxm},type,opts]

RegisterDatai[
{target_?ArrayQ,maskt_?ArrayQ,voxt:{_?NumberQ,_?NumberQ,_?NumberQ}},{moving_?ArrayQ,maskm_?ArrayQ,voxm:{_?NumberQ,_?NumberQ,_?NumberQ}},
type_,OptionsPattern[]]:=Module[{
	tdir, tempdir, elastix, targetFile, parstring, outputImg, iterations, resolutions,
	histogramBins, numberSamples, derivativeScaleA, derivativeScaleB, interpolationOrder,
	method, bsplineSpacing, data, vox, dimmov, dimtar, dimmovm, dimtarm, inpfol, movfol, outfol, 
	fixedF, movingF, outF, parF, depth, index, error, regpars, lenMeth, command, outfile, 
	fmaskF, mmaskF, maske, maske2, w, openCL, gpu, pca, slash},
	
	w={{0,0,0,0,0,0,1,1,1,0,0,0}};
	slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
	
	(*set error*)
	error=False;
	maske=False;
	
	(*get option values*)
	tdir=OptionValue[TempDirectory];
	outputImg=ToLowerCase[ToString[OptionValue[OutputImage]]];
	
	method=OptionValue[MethodReg];
	(*Print[method];*)
	
	bsplineSpacing=OptionValue[BsplineSpacing];
	bsplineSpacing=If[!ListQ[bsplineSpacing],ConstantArray[bsplineSpacing,3],bsplineSpacing];
	derivativeScaleB=OptionValue[BsplineDirections];
	derivativeScaleA=OptionValue[AffineDirections];
	
	{openCL,gpu}=OptionValue[UseGPU];
	gpu=If[gpu===Automatic,0,gpu];
	pca=OptionValue[PCAComponents];
	
	(*Print[{derivativeScaleA,derivativeScaleB}];*)
	
	iterations=OptionValue[Iterations];
	resolutions=OptionValue[Resolutions];
	histogramBins=OptionValue[HistogramBins];
	numberSamples=OptionValue[NumberSamples];
	interpolationOrder=OptionValue[InterpolationOrderReg];
	regpars={iterations,resolutions,histogramBins,numberSamples,interpolationOrder};
	
	dimmov=Dimensions[moving];
	dimtar=Dimensions[target];
	dimmovm=Dimensions[maskm];
	dimtarm=Dimensions[maskt];
	
	(*find the elastix program*)
	elastix=FindElastix[];
	If[elastix=="error",error=True;Message[RegisterData::elastix];];
	
	(*create temp directory*)
	tdir=(If[StringQ[tdir],tdir,"Default"]/. {"Default"->$TemporaryDirectory});
	
	tdir=If[Last[FileNameSplit[tdir]] === "DTItoolsReg" || Last[FileNameSplit[tdir]] === "anat",
		tdir,tdir<>slash<>"DTItoolsReg"
	];
	
	If[DirectoryQ[tdir],DeleteDirectory[tdir,DeleteContents->True]];
	tempdir=CreateDirectory[tdir]<>slash;
	If[!DirectoryQ[tempdir],Message[RegisterData::dir];Return[Message[DiffusionReg::fatal]]];
	
	(*check registration method*)
	method=If[StringQ[method],{method},method];
	If[!MemberQ[{"rigid","affine","rigidDTI","affineDTI","bspline","cyclyc","translation","PCA"},#],
		Message[RegisterData::met,#];
		Return[Message[DiffusionReg::fatal],Module]
		]&/@method; 
	lenMeth=Length[method];
	
	(*only cyclyc is possible*)
	If[MemberQ[method,"cyclyc"]&&lenMeth>1,error=True];
	
	(*create parameter list*)
	regpars=If[NumberQ[#],ConstantArray[#,lenMeth],
	If[Length[#]==lenMeth,#,Message[RegisterData::par,#,lenMeth];Return[Message[DiffusionReg::fatal]];
	]]&/@regpars;
	
	(*check mask dimensions*)
	
	(*if error quit*)
	If[error,Return[Message[RegisterData::fatal]]];
	If[OptionValue[PrintTempDirectory],PrintTemporary["using as temp directory: "<>tdir]];
	
	(*create parameter files*)
	parF=MapThread[(
	parstring=ParString[#2,{#1,outputImg},{bsplineSpacing,derivativeScaleB,derivativeScaleA,pca},{openCL,gpu}];
	parF="parameters-"<>#1<>".txt";
	Export[tempdir<>parF,parstring];
	parF
	)&,{method,Transpose[regpars]}];
	
	(*create target file*)
	depth=If[type==="cyclyc"||type==="PCA",ToString[ArrayDepth[target]-1]<>"D-t",ToString[ArrayDepth[target]]<>"D"];
	fixedF="target-"<>depth<>".nii";
	targetFile=tempdir<>fixedF;
	
	ExportNii[target,voxt,targetFile];
	
	(*perform registration*)
	Switch[type,
	
	"vol",
	{inpfol,movfol,outfol}={"","",""};
	{movingF,outF}={"moving-"<>depth<>".nii","result-"<>depth<>".nii.gz"};
	ExportNii[moving,voxm,tempdir<>movingF];
	
	{fmaskF,mmaskF}={"",""};
	
	(*check if target mask is needed*)
	If[dimtarm == dimtar && maskt!={1},
		fmaskF="targetMask.nii";
		ExportNii[maskt,voxm,tempdir<>fmaskF]];
	
	(*check if moving mask is needed*)
	If[(dimmovm == dimmov && maskm!={1}),
		mmaskF="moveMask.nii";
		ExportNii[maskm,voxm,tempdir<>mmaskF]];
	
	RunElastix[elastix,tempdir,parF,{inpfol,movfol,outfol},{fixedF,movingF,outF},{fmaskF,mmaskF}];
	{data,vox}=ImportNii[tempdir<>outfol<>outF];
	,
	
	"cyclyc",
	{inpfol,movfol,outfol}={"","",""};
	{fmaskF,mmaskF}={"",""};
	{movingF,outF}={"moving-"<>depth<>".nii","result-"<>depth<>".nii.gz"};
	ExportNii[moving,voxm,tempdir<>movingF];
	If[maskm!={1},mmaskF="moveMask.nii";ExportNii[maskm,voxm,tempdir<>mmaskF]];
	If[maskt!={1},fmaskF="targetMask.nii";ExportNii[maskt,voxm,tempdir<>fmaskF]];
	RunElastix[elastix,tempdir,parF,{inpfol,movfol,outfol},{fixedF,movingF,outF},{fmaskF,mmaskF}];
	{data,vox}=ImportNii[tempdir<>outfol<>outF];
	data=ToPackedArray[data];
	,
	
	"PCA",
	{inpfol,movfol,outfol}={"","",""};
	{fmaskF,mmaskF}={"",""};
	{movingF,outF}={"moving-"<>depth<>".nii","result-"<>depth<>".nii.gz"};
	ExportNii[moving,voxm,tempdir<>movingF];
	If[maskm!={1},mmaskF="moveMask.nii";ExportNii[maskm,voxm,tempdir<>mmaskF]];
	If[maskt!={1},fmaskF="targetMask.nii";ExportNii[maskt,voxm,tempdir<>fmaskF]];
	RunElastix[elastix,tempdir,parF,{inpfol,movfol,outfol},{fixedF,movingF,outF},{fmaskF,mmaskF}];
	{data,vox}=ImportNii[tempdir<>outfol<>outF];
	data=ToPackedArray[data];
	,
	
	"series",
	
	inpfol="";
	{fmaskF,mmaskF}={"",""};
	{movingF,outF}={"moving-"<>depth<>".nii","result-"<>depth<>".nii.gz"};
	
	(*export one mask for every volume in the series*)
	If[dimtarm == dimtar && maskt!={1},
	fmaskF="targetMask.nii";
	ExportNii[maskt,voxm,tempdir<>fmaskF]];
	
	(*check if mask needs to be exported for each volume*)
	maske=(dimmovm == dimmov && maskm!={1});
	(*check if same mask for all volumes*)
	maske2=(dimmovm == Drop[dimmov,1] && maskm!={1});
	
	(*export data*)
	{command,outfile}=Transpose@(
	(
		index=StringPad[#];
		movfol=outfol="vol"<>index;
		CreateDirectory[tempdir<>outfol];
		ExportNii[moving[[#]],voxm,tempdir<>movfol<>slash<>movingF];
		
		(*export mask*)
		If[maske,
		mmaskF=movfol<>slash<>"moveMask.nii";
		ExportNii[maskm[[#]],voxm,tempdir<>mmaskF];
		];
		
		If[maske2,
		mmaskF=movfol<>slash<>"moveMask.nii";
		ExportNii[maskm,voxm,tempdir<>mmaskF];
		];
		
		ElastixCommand[elastix,tempdir,parF,{inpfol,movfol,outfol},{fixedF,movingF,outF},{fmaskF,mmaskF}]
	)&/@Range[Length[moving]]);
	(*create and run batch*)
	RunBatfile[tempdir,command];
	
	(*Import data*)
	data=(First@ImportNii[#])&/@outfile;
	
	If[OptionValue[OutputTransformation], w = ReadTransformParameters[tempdir]];

	];
	
	data=ToPackedArray[Chop[Clip[data,MinMax[moving]]]];
	
	If[OptionValue[DeleteTempDirectory],DeleteDirectory[tempdir,DeleteContents->True]];
	If[OptionValue[OutputTransformation], {data,w},	data]
	
]


(* ::Subsubsection::Closed:: *)
(*RegisterDataSplit*)


Options[RegisterDataSplit] = Join[Options[RegisterData],{SplitMethod->"Mean"}];

SyntaxInformation[RegisterDataSplit] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

RegisterDataSplit[targeti_, movingi_, opts : OptionsPattern[]] := Block[{
	reg, mov,
	target ,maskT, voxT,
	moving, maskM, voxM,
	targetl, targetr, maskTl, maskTr, cut1,
	movingl, movingr, maskMl, maskMr, cut2,
	regl, regr, movl, movr
	},
	
	(*prepare the input*)
	{target ,maskT, voxT}=SplitInput[targeti];
	{moving, maskM, voxM}=SplitInput[movingi];
	
	(*find the common split*)	
	{targetl, targetr, cut1}=CutData[target];
	{movingl, movingr, cut2}=CutData[moving];
	
	{cut1, cut2} = Switch[OptionValue[SplitMethod],
		"target", Round[{cut1, (cut1 voxT[[2]])/voxM[[2]]}],
		"moving", Round[{(cut2 voxM[[2]])/voxT[[2]], cut2}],
		"nearest", Round[First@Nearest[{cut1 Last@voxT, cut2 Last@voxM}, (Last@Dimensions[target]/2) Last[voxT]]/{Last@voxT, Last@voxM}],
		_, Round[Mean[{cut1 voxT[[2]], cut2 voxM[[2]]}]/{voxT[[2]], voxM[[2]]}]
		];
	
	(*cut data*)
	{targetl, targetr, cut1}=CutData[target,cut1];
	{movingl, movingr, cut2}=CutData[moving,cut2];
	(*cut masks*)
	{maskTl, maskTr}=If[maskT==={1},{{1},{1}},CutData[maskT,cut1][[;;-2]]];
	{maskMl, maskMr}=If[maskM==={1},{{1},{1}},CutData[maskM,cut2][[;;-2]]];
	
	(*register left part*)
	regl = RegisterData[{targetl, maskTl, voxT}, {movingl, maskMl, voxM},  Sequence@@FilterRules[{opts}, Options[RegisterData]]];
	(*register right part*)
	regr = RegisterData[{targetr, maskTr, voxT}, {movingr, maskMr, voxM},  Sequence@@FilterRules[{opts}, Options[RegisterData]]];
	
	StichData[regl,regr]
	
  ]


(* ::Subsubsection::Closed:: *)
(*SplitInput*)


SplitInput[input_]:=Module[{data,mask,vox},
	(*split the target data*)
	If[ArrayQ[input],
		(*data*)
		data=input;mask={1};vox={1,1,1};
		,
		If[Length[input]==2 && ArrayQ[input[[1]]] && Length[input[[2]]]==3,
			(*data and vox*)
			data=input[[1]];mask={1};vox=input[[2]];
			,
			(*data, mask and vox*)
			data=input[[1]];mask=input[[2]];vox=input[[3]];
		]
	];
	{data,mask,vox}
]


(* ::Subsection:: *)
(*TransformData*)


(* ::Subsubsection::Closed:: *)
(*TransformData*)


Options[TransformData] = {TempDirectory -> "Default", FindTransform -> "Auto", DeleteTempDirectory -> "All",PrintTempDirectory->True}

TransformData[{data_, vox_}, OptionsPattern[]] := Module[{tdir, command, output,slash},
	slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
	
	(*define the directory*)
	tdir = OptionValue[TempDirectory];
	tdir = (If[StringQ[tdir], tdir, "Default"] /. {"Default" -> $TemporaryDirectory}) <>slash<>"DTItoolsReg"<>slash<>"transform";
	
	(*create and print the directory*)
	If[OptionValue[PrintTempDirectory],PrintTemporary[tdir]];
	If[DirectoryQ[tdir],DeleteDirectory[tdir,DeleteContents->True]];
	CreateDirectory[tdir];
	
	(*Export and transform*)
	ExportNii[data, vox, tdir <> slash <> "trans.nii"];
	command = TransformixCommandInd[tdir];

	RunProcess[$SystemShell, "StandardOutput", command];
	output = ToPackedArray[ImportNii[tdir <> slash <> "result.nii"][[1]]];
	
	(*Delete temp directory*)
	Switch[OptionValue[DeleteTempDirectory],
		"All", DeleteDirectory[FileNameTake[tdir, {1, -2}],  DeleteContents -> True],
		"Trans", DeleteDirectory[tdir, DeleteContents -> True],
		_, Null];
		
	(*give the output*)
	Chop[Clip[output,MinMax[data]],10^-6]
]


(* ::Subsubsection::Closed:: *)
(*TransformixCommandInd*)


TransformixCommandInd[tempDir_] := Block[{transformix, transfile,transFol},
	transformix = FindTransformix[];
	transFol = StringDrop[DirectoryName[transformix, 2], -1];
	transfile = Last[SortBy[
		FileNames["TransformParameters*", FileNameTake[tempDir, {1, -2}]],
		FileDate[#, "Creation"] &]];
	
	Switch[$OperatingSystem,
		"Windows",
		"@ \"" <> transformix <>
		"\" -in \"" <> First[FileNames["trans*", tempDir]] <>
		"\" -out \"" <> tempDir <>
		"\" -tp \"" <> transfile <>
		"\" -def all" <>
		" > \"" <> tempDir <> "\\outputT.txt\" \n exit \n"
		,
		"MacOSX",
		"export PATH="<>transFol<>"/bin:$PATH 
		export DYLD_LIBRARY_PATH="<>transFol<>"/lib:$DYLD_LIBRARY_PATH \n"<>
		transformix <>
		" -in " <> First[FileNames["trans*", tempDir]] <>
		" -out " <> tempDir <>
		" -tp " <> transfile <>
		" > " <> tempDir <> "/outputT.txt \n"
	]
  ]


(* ::Subsubsection::Closed:: *)
(*ParametersToTransform*)


ParametersToTransform[w_, opt_] := 
 Block[{tx, ty, tz, rx, ry, rz, sx, sy, sz, gx, gy, gz, T, R, G, S, 
   Rx, Ry, Rz, Gx, Gy, Gz, mat, mats, matL},
  {rx, ry, rz, tx, ty, tz, sx, sy, sz, gx, gy, gz} = w;
  {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  rx = rx Degree; ry = ry Degree; rz = rz Degree;
  T = {
    {1, 0, 0, tx},
    {0, 1, 0, ty},
    {0, 0, 1, tz},
    {0, 0, 0, 1}};
  Rx = {
    {1, 0, 0, 0},
    {0, Cos[rx], Sin[rx], 0},
    {0, -Sin[rx], Cos[rx], 0},
    {0, 0, 0, 1}};
  Ry = {
    {Cos[ry], 0, -Sin[ry], 0},
    {0, 1, 0, 0},
    {Sin[ry], 0, Cos[ry], 0},
    {0, 0, 0, 1}};
  Rz = {
    {Cos[rz], Sin[rz], 0, 0},
    {-Sin[rz], Cos[rz], 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}};
  R = Rx.Ry.Rz;
  Gx = {
    {1, 0, gx, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}};
  Gy = {
    {1, 0, 0, 0},
    {gy, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}};
  Gz = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, gz, 1, 0},
    {0, 0, 0, 1}};
  G = Gx.Gy.Gz;
  S = {
    {sx, 0, 0, 0},
    {0, sy, 0, 0},
    {0, 0, sz, 0},
    {0, 0, 0, 1}};
  
  mat = Switch[opt,
    "Full", T.R.G.S,
    "Rotation", R,
    _, R
    ];
  
  mats = mat[[1 ;; 3, 1 ;; 3]];
  
  (MatrixPower[mats.Transpose[mats], -(1/2)].mats)
  ]



(* ::Subsection:: *)
(*RegisterDiffusionData/Split*)


(* ::Subsubsection::Closed:: *)
(*RegisterDiffusionData*)


Options[RegisterDiffusionData] = 
  Join[Options[RegisterData] /. {{1, 1, 1} -> {0, 1, 1}, "affine" -> "affineDTI", "rigid" -> "rigidDTI"},
   {IterationsA -> 1000, 
   	ResolutionsA -> 1, 
   	HistogramBinsA -> 64, 
    NumberSamplesA -> 20000, 
    InterpolationOrderRegA -> 1, 
    MethodRegA -> {"rigidDTI", "bspline"},
    RegistrationTarget->"Fist"
    }];

SyntaxInformation[RegisterDiffusionData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

(*No anatomical data, goto Registerdata*)
RegisterDiffusionData[
	{dtidata_?ArrayQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}},
	opts : OptionsPattern[]] := RegisterDiffusionData[{dtidata, {1}, vox},opts]

RegisterDiffusionData[
	{dtidata_?ArrayQ, dtimask_?ArrayQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}},
	opts:OptionsPattern[]] := (
	RegisterData[{dtidata, dtimask, vox},(*OutputTransformation->True,*) 
		MethodReg-> (OptionValue[MethodReg] /. {"affine" -> "affineDTI", "rigid" -> "rigidDTI"}),
		AffineDirections -> {1, 1, 1},
		FilterRules[{opts}, Options[RegisterData]]]
		)

(**)
RegisterDiffusionData[
  {dtidata_?ArrayQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}},
  {anatdata_?ArrayQ, voxa : {_?NumberQ, _?NumberQ, _?NumberQ}},
  opts : OptionsPattern[]
  ] := RegisterDiffusionData[{dtidata, {1}, vox}, {anatdata, {1},voxa}, opts]

RegisterDiffusionData[
  {dtidata_?ArrayQ, dtimask_?ArrayQ, 
   vox : {_?NumberQ, _?NumberQ, _?NumberQ}},
  {anatdata_?ArrayQ, voxa : {_?NumberQ, _?NumberQ, _?NumberQ}},
  opts : OptionsPattern[]
  ] := RegisterDiffusionData[{dtidata, dtimask, vox}, {anatdata, {1}, voxa}, opts]

RegisterDiffusionData[
  {dtidata_?ArrayQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}},
  {anatdata_?ArrayQ, anatmask_?ArrayQ, 
   voxa : {_?NumberQ, _?NumberQ, _?NumberQ}},
  vox_, opts : OptionsPattern[]
  ] := RegisterDiffusionData[{dtidata, {1}, vox}, {anatdata, anatmask, voxa}, opts]

RegisterDiffusionData[
  {dtidata_?ArrayQ, dtimask_?ArrayQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}},
  {anatdata_?ArrayQ, anatmask_?ArrayQ, voxa : {_?NumberQ, _?NumberQ, _?NumberQ}},
  opts : OptionsPattern[]] := Module[{
  	dtidatar, tempDir, tempDira, volDirs, w,tFilesA, tFilesD, dtidatarA, cmd, target, movingdata, slash
  	},
  (*Print["RegisterDiffusionData"];*)
  
  slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
  (*get the current temp dir and define the anat tempdir*)
  tempDir = OptionValue[TempDirectory];
  tempDir = (If[StringQ[tempDir], tempDir, "Default"]/. {"Default"->$TemporaryDirectory})<>slash<>"DTItoolsReg";
  tempDira = tempDir <> slash <> "anat";

  (*perform DTI registration*)
  dtidatar = RegisterData[{dtidata, dtimask, vox},
    TempDirectory -> tempDir, 
    DeleteTempDirectory -> False, 
    OutputTransformation->OptionValue[OutputTransformation], 
    MethodReg-> (OptionValue[MethodReg] /. {"affine" -> "affineDTI", "rigid" -> "rigidDTI"}),
    (*AffineDirections -> {1, 1, 1},*)
    FilterRules[{opts} , Options[RegisterData]]];

  If[OptionValue[OutputTransformation],{dtidatar,w}=dtidatar];
    
  target = OptionValue[RegistrationTarget];
  movingdata=If[ListQ[target] && AllTrue[target, IntegerQ] && Min[target] > 0 && Max[target] <= Length[dtidatar[[1]]],
  	Median /@ dtidatar[[All, DeleteDuplicates[target]]],
  	Switch[target,
  		"Median", Median /@ dtidatar,
  		"First", dtidatar[[All, 1]],
  		_, Mean /@ dtidatar
  		]];
  
  (*perform anat registration*)
  RegisterData[{anatdata, anatmask, voxa}, {movingdata, dtimask, vox},
   TempDirectory -> tempDira, 
   DeleteTempDirectory -> False,
   Iterations -> OptionValue[IterationsA], 
   Resolutions -> OptionValue[ResolutionsA],
   HistogramBins -> OptionValue[HistogramBinsA], 
   NumberSamples -> OptionValue[NumberSamplesA],
   InterpolationOrderReg -> OptionValue[InterpolationOrderRegA],
   BsplineSpacing -> OptionValue[BsplineSpacing], 
   BsplineDirections -> OptionValue[BsplineDirections],
   AffineDirections -> OptionValue[AffineDirections],
   MethodReg -> OptionValue[MethodRegA], 
   FilterRules[{opts}, Options[RegisterData]]
   ];
  
  (*transform all diffusion files to anatomy*)
  
  (*export diffusion reg target*)
  CreateDirectory[tempDir<>slash<>"vol0000"];
  ExportNii[dtidatar[[All,1]],vox,tempDir<>slash<>"vol0000"<>slash<>"moving-3D.nii"];
  
  (*get vol folders and anat transform files*)
  volDirs = FileNames["vol*", tempDir, 1];
  tFilesA = FileNames["TransformParameters*", tempDira];
  
  (*create Final Transform files*)
  (
     tFilesD = FileNames["TransformParameters*", #];
     ConcatenateTransformFiles[Join[tFilesD, tFilesA], #]
     ) & /@ volDirs;
  
  (*call transformix*)
  cmd = TransformixCommand[tempDir];
  PrintTemporary["Combining transformations"];
  RunBatfileT[tempDir, cmd];
  
  (*import dti data in anat space*)
  dtidatarA = Transpose[ImportNii[#][[1]] & /@ FileNames["resultA*", tempDir, 2]];
  
  (*finalize by deleting temp director*)
  If[OptionValue[DeleteTempDirectory],DeleteDirectory[tempDir,DeleteContents->True]];
  
  (*output data*)
  If[OptionValue[OutputTransformation],
  	{dtidatar, dtidatarA, w},
  	{dtidatar, dtidatarA}
  ]
]


(* ::Subsubsection::Closed:: *)
(*RegisterDiffusionDataSplit*)


Options[RegisterDiffusionDataSplit] := Options[RegisterDiffusionData];

SyntaxInformation[RegisterDiffusionDataSplit] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

RegisterDiffusionDataSplit[data_, vox: {_?NumberQ, _?NumberQ, _?NumberQ}, opts : OptionsPattern[]] := 
  Block[{datal, datar, cut},
   {datal, datar, cut} = CutData[data];
   datal = RegisterDiffusionData[{datal, vox}, opts];
   datar = RegisterDiffusionData[{datar, vox}, opts];
   StichData[datal, datar]
   ];

RegisterDiffusionDataSplit[{data_, mask_, vox: {_?NumberQ, _?NumberQ, _?NumberQ}}, opts : OptionsPattern[]] := Block[
	{datal, datar, cut,maskr,maskl},
	
	{datal, datar, cut} = CutData[data];
	{maskl, maskr, cut} = CutData[mask,cut];
	datal = RegisterDiffusionData[{datal, maskl, vox}, opts];
	datar = RegisterDiffusionData[{datar, maskr, vox}, opts];
	StichData[datal, datar]
   ];

RegisterDiffusionDataSplit[{data_, vox: {_?NumberQ, _?NumberQ, _?NumberQ}}, {dataa_, voxa: {_?NumberQ, _?NumberQ, _?NumberQ}}, opts : OptionsPattern[]] := Block[
   	{datal, datar, dataal, dataar, cut1, cut2},
   	
   (*find cuts*)
   {datal, datar, cut1} = CutData[data];
   {dataal, dataar, cut2} = CutData[dataa];
   (*align cuts*)
   {cut1,cut2}=Round[First@Nearest[{cut1 Last@vox, cut2 Last@voxa}, Round[Last@Dimensions[data]/2] Last@vox] / {Last@vox, Last@voxa}];
   (*{cut1,cut2}=Round[Mean[{cut1 vox[[2]], cut2 voxa[[2]]}]/{vox[[2]],voxa[[2]]}];*)
   
   (*cut with the aligned cuts*)
   {datal, datar, cut1} = CutData[data, cut1];
   {dataal, dataar, cut2} = CutData[dataa, cut2];
   
   datal = RegisterDiffusionData[{datal, vox}, {dataal, voxa}, opts][[2]];
   datar = RegisterDiffusionData[{datar, vox}, {dataar, voxa}, opts][[2]];
   StichData[datal, datar]
   ];

RegisterDiffusionDataSplit[{data_, mask_, vox: {_?NumberQ, _?NumberQ, _?NumberQ}}, {dataa_, maska_, voxa: {_?NumberQ, _?NumberQ, _?NumberQ}}, opts : OptionsPattern[]] := Block[
	{datal, datar, dataal, dataar, maskl, maskr, maskal, maskar,cut1,cut2},
	
	(*find cuts*)
   {datal, datar, cut1} = CutData[data];
   {dataal, dataar, cut2} = CutData[dataa];
   
   (*align cuts*)
   {cut1,cut2}=Round[First@Nearest[{cut1 Last@vox, cut2 Last@voxa}, Round[Last@Dimensions[data]/2] Last@vox] / {Last@vox, Last@voxa}];
   (*{cut1,cut2}=Round[Mean[{cut1 vox[[2]], cut2 voxa[[2]]}]/{vox[[2]],voxa[[2]]}];*)
   
   (*cut with the aligned cuts*)   
   {datal, datar, cut1} = CutData[data,cut1];
   {maskl, maskr, cut1} = CutData[mask,cut1];
   {dataal, dataar, cut2} = CutData[dataa,cut2];
   {maskal, maskar, cut2} = CutData[maska,cut2];
  
   
   datal = RegisterDiffusionData[{datal, maskl, vox}, {dataal, maskal, voxa}, opts][[2]];
   datar = RegisterDiffusionData[{datar, maskr, vox}, {dataar, maskar, voxa}, opts][[2]];
   
   StichData[datal, datar]
   ];


(* ::Subsection:: *)
(*RegisterDataTransform/Split*)


(* ::Subsubsection::Closed:: *)
(*RegisterDataTransform*)


Options[RegisterDataTransform] = Options[RegisterData];

SyntaxInformation[RegisterDataTransform] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

RegisterDataTransform[target_, moving_, {moving2_, vox_}, opts : OptionsPattern[]] := Block[{reg, mov,tdir, slash},
	reg = RegisterData[target, moving, DeleteTempDirectory -> False, opts];
	
	mov = If[
		ArrayDepth[moving2]==4 && ArrayDepth[reg]==3 ,
		Transpose[TransformData[{#, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False] & /@ Transpose[moving2]] ,
		If[ArrayDepth[moving2]==3 && ArrayDepth[reg]==2 ,
			TransformData[{#, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False] & /@ moving2,
			TransformData[{moving2, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False]
			]
		];
	
	slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
	
	tdir=OptionValue[TempDirectory];
	tdir=(If[StringQ[tdir],tdir,"Default"]/. {"Default"->$TemporaryDirectory})<>slash<>"DTItoolsReg";
	
	If[OptionValue[DeleteTempDirectory],DeleteDirectory[tdir,DeleteContents->True]];		
		
	{reg, mov}
  ]


(* ::Subsubsection::Closed:: *)
(*RegisterDataTransformSplit*)


Options[RegisterDataTransformSplit] = Join[Options[RegisterData],{SplitMethod->"Mean"}];

SyntaxInformation[RegisterDataTransformSplit] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

RegisterDataTransformSplit[targeti_, movingi_, {moving2_, vox_}, opts : OptionsPattern[]] := Block[{reg, mov,
	target ,maskT, voxT,
	moving, maskM, voxM,
	targetl, targetr, maskTl, maskTr, cut1,
	movingl, movingr, maskMl, maskMr, cut2,
	moving2l, moving2r, tdir,
	regl, regr, movl, movr, slash
	},
	
	(*prepare the input*)
	{target ,maskT, voxT}=SplitInput[targeti];
	{moving, maskM, voxM}=SplitInput[movingi];
	
	(*find the common split*)	
	{targetl, targetr, cut1}=CutData[target];
	{movingl, movingr, cut2}=CutData[moving];
	
	{cut1, cut2} = Switch[OptionValue[SplitMethod],
		"Target", Round[{cut1, (cut1 voxT[[2]])/voxM[[2]]}],
		"Moving", Round[{(cut2 voxM[[2]])/voxT[[2]], cut2}],
		"Nearest", Round[First@Nearest[{cut1 Last@voxT, cut2 Last@voxM}, (Last@Dimensions[target]/2) Last[voxT]]/{Last@voxT, Last@voxM}],
		_, Round[Mean[{cut1 voxT[[2]], cut2 voxM[[2]]}]/{voxT[[2]], voxM[[2]]}]
		];
	
	(*cut data*)
	{targetl, targetr, cut1}=CutData[target,cut1];
	{movingl, movingr, cut2}=CutData[moving,cut2];
	(*cut masks*)
	{maskTl, maskTr}=If[maskT==={1},{{1},{1}},CutData[maskT,cut1][[;;-2]]];
	{maskMl, maskMr}=If[maskM==={1},{{1},{1}},CutData[maskM,cut2][[;;-2]]];
	
	(*split the moving2 data*)
	{moving2l, moving2r, cut2} = CutData[moving2, cut2];
	
	(*register left part*)
	regl = RegisterData[{targetl, maskTl, voxT}, {movingl, maskMl, voxM}, DeleteTempDirectory -> False, Sequence@@FilterRules[{opts}, Options[RegisterData]]];
	(*transform the left part*)
	movl = If[ArrayDepth[moving2l] == 4,
		Transpose[TransformData[{#, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False] & /@ Transpose[moving2l]],
		TransformData[{moving2l, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False]
		];
		
	(*register right part*)
	regr = RegisterData[{targetr, maskTr, voxT}, {movingr, maskMr, voxM}, DeleteTempDirectory -> False, Sequence@@FilterRules[{opts}, Options[RegisterData]]];
	(*transform the right part*)
	movr = If[ArrayDepth[moving2r] == 4,
		Transpose[TransformData[{#, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False] & /@ Transpose[moving2r]],
		TransformData[{moving2r, vox}, DeleteTempDirectory -> False, PrintTempDirectory -> False]
		];
	
	slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
	tdir=OptionValue[TempDirectory];
	tdir=(If[StringQ[tdir],tdir,"Default"]/. {"Default"->$TemporaryDirectory})<>slash<>"DTItoolsReg";
	
	If[OptionValue[DeleteTempDirectory],DeleteDirectory[tdir,DeleteContents->True]];	
	
	{StichData[regl,regr],StichData[movl,movr]}
	
  ]



(* ::Subsection::Closed:: *)
(*RegisterCardiacData*)


Options[RegisterCardiacData]=Join[{RegistrationTarget->"Mean"},Options[RegisterData]];

SyntaxInformation[RegisterCardiacData]={"ArgumentsPattern"->{_,OptionsPattern[]}};

(*data only*)
RegisterCardiacData[data_?ArrayQ,opts:OptionsPattern[]]:=RegisterCardiacData[{data,{1},{1,1,1}},opts]
(*data with voxel*)
RegisterCardiacData[{data_?ArrayQ,vox:{_?NumberQ,_?NumberQ,_?NumberQ}},opts:OptionsPattern[]]:=RegisterCardiacData[{data,{1},vox},opts]
(*data with mask*)
RegisterCardiacData[{data_?ArrayQ,mask_?ArrayQ},opts:OptionsPattern[]]:=RegisterCardiacData[{data,mask,{1,1,1}},opts]
(*data with mask and voxel*)
RegisterCardiacData[{data_?ArrayQ,mask_?ArrayQ,vox:{_?NumberQ,_?NumberQ,_?NumberQ}},opts:OptionsPattern[]]:=Block[
{tdir, datar, slices, maskr, i, size, target, slash},

slash = Switch[$OperatingSystem, "Windows", "\\", "MacOSX", "/"];
tdir=OptionValue[TempDirectory];
tdir=(If[StringQ[tdir],tdir,"Default"]/. {"Default"->$TemporaryDirectory})<>slash<>"DTItoolsReg";

If[OptionValue[PrintTempDirectory],PrintTemporary["using as temp directory: "<>tdir]];

slices=Range[Length[data]];
size=Length[data[[1]]];
maskr=If[mask=={1},ConstantArray[1,Dimensions[data[[All,1]]]],mask];

target=If[OptionValue[MethodReg]==="PCA"||OptionValue[MethodReg]==="cyclyc","stack",OptionValue[RegistrationTarget]];

(*monitro over slices*)
Monitor[
	i=0;
	datar=Switch[
	target,
	"Mean",
	(i++;RegisterData[{N[Mean@data[[#]]],maskr[[#]],vox},{data[[#]],maskr[[#]],vox},
		OutputTransformation->False, PrintTempDirectory->False,FilterRules[{opts},Options[RegisterData]]])&/@slices,
	"Median",
	(i++;RegisterData[{N[Median@data[[#]]],maskr[[#]],vox},{data[[#]],maskr[[#]],vox},
		OutputTransformation->False, PrintTempDirectory->False,FilterRules[{opts},Options[RegisterData]]])&/@slices,
	"First",
	(i++;RegisterData[{data[[#,1]],maskr[[#]],vox},{data[[#]],maskr[[#]],vox},
		OutputTransformation->False, PrintTempDirectory->False,FilterRules[{opts},Options[RegisterData]]])&/@slices,
	"stack",
	(i++;RegisterData[{data[[#]],ConstantArray[maskr[[#]],size],vox},
		OutputTransformation->False, PrintTempDirectory->False,FilterRules[{opts},Options[RegisterData]]])&/@slices
	]
	,ProgressIndicator[i,{0,Length[data]}]
];
datar
]


(* ::Subsection:: *)
(*Correct Gradients*)


(* ::Subsubsection::Closed:: *)
(*ReadTransformParameters*)


SyntaxInformation[ReadTransformParameters]={"ArgumentsPattern"->{_}};

ReadTransformParameters[dir_] := Block[{files,filenum,cor},
  files = FileNames["TransformParameters*", dir, 3];
  filenum = If[Length[files] == 1,
  	{1},
  	ToExpression[First[StringCases[FileNameSplit[#][[-2]],DigitCharacter ..]]] & /@ files
  	];
  files = files[[Ordering[filenum]]];
  cor =  
    Partition[
        ToExpression[
         StringSplit[StringTake[Import[#, "Lines"][[3]], {2, -2}]][[
          2 ;;]]], 3][[{1, 4, 3, 2}, {2, 1, 3}]] & /@ files;
  cor[[All, 1]] = cor[[All, 1]]/Degree;
  Flatten /@ cor
  ]


(* ::Subsubsection::Closed:: *)
(*CorrectBmatrix*)


Options[CorrectBmatrix] = {MethodReg -> "Full"}

SyntaxInformation[CorrectBmatrix] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
 
 CorrectBmatrix[bmati_, w_,OptionsPattern[]] := 
  Block[{bmat, trans, bmi, rot, bm, bmnew, bminew},
   bmat = If[Length[First[bmati]] == 7, BmatrixConv[bmati], bmati];
   MapThread[(
      trans = #1;
      bmi = #2;
      rot = ParametersToTransform[trans, OptionValue[MethodReg]];
      bm = TensMat[(bmi/{1, 1, 1, 2, 2, 2})[[{2, 1, 3, 4, 6, 5}]]];
      bmnew = rot.bm.Transpose[rot];
      
      (*Print[MatrixForm/@Round[{bm,bmnew,rot},.00001]];*)
      bminew = ({1, 1, 1, 2, 2, 2} TensVec[bmnew]
      )[[{2, 1, 3, 4, 6, 5}]]
      ) &, {w, bmat}]
   ]


(* ::Subsubsection::Closed:: *)
(*CorrectGradients*)


Options[CorrectGradients] = {MethodReg -> "Rotation"}

SyntaxInformation[CorrectGradients] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
 
 CorrectGradients[grad_, w_,OptionsPattern[]] := 
  Block[{gr,grnew,trans,rot},
   MapThread[(
      ParametersToTransform[#1, OptionValue[MethodReg]].#2
      ) &, {w, grad}]
   ]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
