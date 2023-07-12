(* ::Package:: *)

(* ::Title:: *)
(*QMRITools TaggingTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`TaggingTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`TaggingTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


AnnalyzeTagging::usage = 
"AnnalyzeTagging[gridC] ...";

CalculateDispacementParameters::usage = 
"CalculateDispacementParameters[{motx, moty}, mask] ...";


(* ::Subsection::Closed:: *)
(*Options*)


HistoryWeighting::usage = 
"HistoryWeighting is an options for AnnalyzeTagging."

MonitorTagging::usage = 
"MonitorTagging is an options for AnnalyzeTagging."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*CalculateWaveVector*)


CalculateWaveVector[dat_]:=Block[{all,dim,cent,pts,vecs,phi,rw,vals,vs,ang,n,npts,alli},
	(*make power Spectrum*)
	alli=Rescale@GaussianFilter[Log@Abs@Shift@FFT@Mean[dat[[All,1]]],Round[Max[Dimensions[dat[[1,1]]]]/50]];
	cent=(dim=Dimensions[alli])/2;
		
	(*If less then 5 points retry*)
	n=0.75;npts=1;
	While[npts<5 && n>0.5,
		n-=0.05;
		all=Mask[alli,n];
		(*find the peaks in the power spectrum, one in center and four arround*)
		pts=ComponentMeasurements[Image@all,"Centroid"][[All,2]];
		pts=Nearest[pts,cent,5];
		npts=Length[pts];
	];
		
	(*Print[Image[all]];*)
	(*define the vectors of the points*)
	vecs=(#-pts[[1]])&/@pts[[2;;]];
	vecs=Sign[Sign[#[[2]]+0.00001]]#&/@vecs;
	
	(*calculate the angles and size*)
	phi=Mod[ArcTan[#[[1]],#[[2]]]+Pi/2,Pi]-Pi/2&/@vecs;
	rw=Norm/@vecs;
	
	(*sort for angle, first two are closest to horizontal*)
	vals=Sort@Transpose@{phi,rw};
	
	(*find the mean vectors*)
	vs={Mean[vals[[1;;2]]],Mean[vals[[3;;]]]};
	vs=If[-Pi/4<#[[1]]<Pi/4,dim[[1]]/#[[2]],dim[[2]]/#[[2]]]{Cos[#[[1]]],Sin[#[[1]]]}&/@vs;
	
	(*make the vectors orthogonal*)
	ang=(VectorAngle[vs[[1]],vs[[2]]]-Pi/2);
	
	Chop[{RotationMatrix[ang] . vs[[1]],RotationMatrix[-ang] . vs[[2]]}]
]


(* ::Subsection::Closed:: *)
(*AnnalyzeTagging*)


Options[AnnalyzeTagging] = {HistoryWeighting -> 0.7, MonitorTagging -> True}

AnnalyzeTagging[gridC_, OptionsPattern[]] := Block[{
   waveVecs, im0, im1, im1p, dux, duy, ux, uy, Uall, W0, W1, imi, uxF, uyF, im, w, w0, uxAll, uyAll, s, f, smax, fmax, alpha
   },
  (*get the wave vec is needed*)
  alpha = OptionValue[HistoryWeighting];
  waveVecs = CalculateWaveVector[gridC];
  (*Dynamic monitoring*)
  im0 = im1 = im1p = dux = duy = ux = uy = 0. gridC[[1, 1]];
  s = f = 0;
  smax = Length[gridC];
  fmax = Length[gridC[[1]]];
  
  PrintTemporary[Dynamic[If[OptionValue[MonitorTagging],
     Column[{
       "Slice: " <> ToString[s] <> "/" <> ToString[smax] <> " - Frame: " <> ToString[f] <> "/" <> ToString[fmax],
       Grid@{{MakeImage[im0], MakeImage[im1p]}, {MakeImage[dux], MakeImage[duy]}, {MakeImage[ux], MakeImage[uy]}}
       }, Alignment -> Center]
     ,
     "Slice: " <> ToString[s] <> "/" <> ToString[smax] <> " - Frame: " <> ToString[f] <> "/" <> ToString[fmax]
     ]]];
  
  (*perform analisys*)
  Uall = Table[
    Table[
     {s, f} = {slice, frame};
     If[frame == 1,
      {W0, im0} = ImageToWave[gridC[[slice, frame]], waveVecs];
      ux = uy = dux = duy = 0 im0;
      im1 = imi = im0;
      W1 = W0;
      ,
      (*get the next frame*)
      im1 = gridC[[slice, frame]];
      
      (*extend the motion field*)
      {uxF, uyF} = WExtend[ux, uy, W0, waveVecs];
      
      (*Distort the frame according to the accumulated motion up to the previous frame*)
      im1 = MotionToImage[im1, uxF, uyF];
      {W1, im1} = ImageToWave[im1, waveVecs];
      imi = im1;
      
      (*calculate the displacement between the two frames*)
      {dux, duy} = ImageToMotion[im0, im1, waveVecs];
      
      (*displace the image and increment displacement*)
      im1 = im1p = MotionToImage[im1, dux, duy];
      W1 = MotionToImage[W1, dux, duy];
      ux = uxF + dux;
      uy = uyF + duy;
      
      (*update refernece images*)
      im0 = (1-alpha) im0 + alpha im1;
      W0 = (1-alpha) W0 + alpha W1;
      ];
     {im0, im1, W0, W1, ux, uy}
     , {frame, 1, fmax, 1}]
    , {slice, 1, smax, 1}];
  
  (*create output*)
  {im0, im, w0, w, uxAll, uyAll} = Transpose[Uall, {2, 3, 1, 4, 5}];
  {{uxAll, uyAll}, {im0, im, w0, w}}
  ]

(*plot image*)
MakeImage[im_] := Image[Rescale@N@im, ImageSize -> 200]


(* ::Subsection::Closed:: *)
(*ImageToMotion*)


ImageToMotion[im0_, im1_, waveVecs_] := Block[{
   imRef, imDefi, imDef, dim, row, col, band, uxi, uyi, x, y, dGrid, ker, pad, iwRef, iwDef, win,
   mask, bf, wrg, xmin, xmax, ymin, ymax, iwAux, jm1, jm2, jmD1, jmD2, f1, f2, g1, g2, fmap, gmap, wmap, thrW, bf2, dUr,
   weightL, weightH, phase, tmap, phaseEst, nmap, dU, dUx, dUy
   },
  
  (*initialize parameters and band filters*)
  imRef = im0;
  imDefi = imDef = im1;
  dim = {row, col} = Dimensions[imRef];
  uxi = uyi = 0. imRef;
  (*bandfilters are shifted to center*)
  band = BandFilter[dim, waveVecs, True];
  
  (*loop over waveVectors*)
  Table[
   (*get the filters*)
   {mask, bf, wrg} = band[[All, i]];
   
   (*get wave vecotr properties*)
   {x, y} = Normalize[waveVecs[[i]]];
   dGrid = Norm[waveVecs[[i]]];
   {ker, pad} = GetKernel[waveVecs[[i]], 1];
   win = GetWin[ker, pad, dim];
   
   (*fourier transform of images with shift to center*)
   iwRef = Shift@FFT[win imRef];
   iwDef = Shift@FFT[win imDef];
   
   (*crop to frequncey area*)
   {{xmin, xmax}, {ymin, ymax}} = (MinMax[#] + {-1, 1}) & /@ Transpose[Position[mask, 1]];
   {iwRef, iwDef, mask, bf, wrg} = #[[xmin ;; xmax, ymin ;; ymax]] & /@ {iwRef, iwDef, mask, bf, wrg};
   bf2 =mask( bf/( Sqrt[wrg]+10^-10));
   
   (*weigthed background filter to calculate motion and derivatives*)
   (*reference image*)
   iwAux = iwRef bf2; jm1 = IFFT[iwAux];
   iwAux = wrg iwAux; jmD1 = IFFT[iwAux];
   (*deformed images*)
   iwAux = iwDef bf2; jm2 = IFFT[iwAux];
   iwAux = wrg iwAux; jmD2 = IFFT[iwAux];
   
   (*make abs and calculated needed maps*)
   {f1, f2, g1, g2} = Abs[{jm1, jm2, jmD1, jmD2}];
   fmap = Sqrt[f1^2 + f2^2];
   gmap = Sqrt[g1^2 + g2^2];
   wmap = f1*f2 + g1*g2;
   
   (*find Threshold for weighting map and calculate wheight*)
   thrW = 0.3 Mean[Flatten[wmap]];
   weightL = DevideNoZero[thrW, (wmap + thrW)];
   weightH = 1 - weightL;
   
   (*calculate phases and anti aliasing*)
   phase = Arg[jm1 Conjugate[jm2] + jmD1 Conjugate[jmD2]];
   tmap = fmap phase;
   phaseEst = ListConvolve[ker, tmap, pad, 0.]/ListConvolve[ker, fmap, pad, 0.];
   tmap = fmap (Mod[phase + Pi - phaseEst, 2 Pi] + phaseEst - Pi);
   nmap = (2 Pi/dGrid) gmap;
   
   (*estimate dispacement map U with low spatial smaling rate, in low regions smooth extra*)
   dUr = (tmap weightH + ListConvolve[ker, tmap, pad, 0.] weightL)/(nmap weightH + ListConvolve[ker, nmap, pad, 0.] weightL );
   
   dU = RescaleData[dUr, dim, InterpolationOrder -> 1];
   
   (*update motion and deform image*)
   uxi += x dU; uyi += y dU;
   (*imDef=MotionToImage[imDefi,uxi,uyi]*)
   
   , {i, 1, Length[waveVecs]}];
  
  (*output*)
  {uxi, uyi}
  ]


(* ::Subsection::Closed:: *)
(*ImageToWave*)


ImageToWave[data_, waveVecs_] := Block[{
	dim, ker, pad, win, iw, bf, wrg, back, grid, iwX, fiwX, mask
	},
	
	(*get filter kernel and make band filters*)
	dim = Dimensions[data];
	{ker, pad} = GetKernel[waveVecs, 1];
	{mask, bf, wrg} = BandFilter[dim, waveVecs];
	win = GetWin[ker, pad, dim];
	
	(*perform FFT and aply band pass filters in both waveVec direcions*)
	iw = FFT[data];
	iwX = iw # & /@ bf;
	fiwX = IFFT /@ iwX;
	
	(*homoginize background and filter*)
	back = Total[Abs[fiwX]^2];
	back = back/(back + Median[Flatten[back]]);
	back = win ListConvolve[ker, back, pad, 0.];
	
	(*weigthing of tag grid for background*)
	grid = Total[Re[fiwX]];
	grid = grid back;
	
	(*output with removed DC offset from grid*)
	{back, grid - Mean@Flatten@grid}
]


(* ::Subsection::Closed:: *)
(*MotionToImage*)


MotionToImage[data_, ux_, uy_] := Block[{dim, xcor, ycor, intdata, intFun, zx, zy},
	(*get the cordinate system*)
	dim = Dimensions[data];
	xcor = Transpose@ConstantArray[Range[1, dim[[1]]], dim[[2]]];
	ycor = ConstantArray[Range[1, dim[[2]]], dim[[1]]];
	
	(*define the linear interpolation function*)
	intFun = ListInterpolation[data, InterpolationOrder -> 1];
	
	(*get the deformed grid*)
	{zx, zy} = {Clip[ux + xcor, {1, dim[[1]]}], Clip[uy + ycor, {1, dim[[2]]}]};
	
	(*interpolate the distorted image*)
	intFun[zx, zy]
]


(* ::Subsection::Closed:: *)
(*MaskToCoordinates*)


MaskToCoordinates[mask_]:=Block[{x,y,z2,z1,p,xm,ym,r,phi,rad, dim, xcor,ycor},
	(*get coordiantes*)
	{x,y}=Transpose@Position[mask,1];
	
	(*fit circle*)
	z2=N[x^2+y^2];
	z1=N[Transpose@{x,y,0x+1}];
	p=(Inverse[(Transpose[z1] . z1)] . (Transpose[z1] . z2));
	
	(*find center and radius*)
	{xm,ym}=p[[1;;2]]/2;
	rad=Sqrt[xm^2+ym^2+p[[3]]];
	
	(*get the polar coordiantes*)
	dim=Dimensions[mask];
	xcor=Transpose@ConstantArray[Range[1,dim[[1]]],dim[[2]]]-xm;
	ycor=ConstantArray[Range[1,dim[[2]]],dim[[1]]]-ym;
	
	{rad,phi}=RotateDimensionsRight[CoordinateTransform["Cartesian"->"Polar",RotateDimensionsLeft[{xcor,ycor}]]];
	{xcor,ycor,rad,phi}
];


(* ::Subsection:: *)
(*Filters*)


(* ::Subsubsection::Closed:: *)
(*Wextend*)


WExtend[ux_, uy_, back_, waveVecs_] := Block[{ker, pad, backF},
	{ker, pad} = GetKernel[waveVecs, 0.5];
	(*filter background*)
	backF=ListConvolve[ker, back, pad, 0.]+10^-5;
	(*filter motion fields weigthed for background*)
	{ListConvolve[ker, ux back, pad, 0.]/ backF, ListConvolve[ker, uy back, pad, 0.]/ backF}
]


(* ::Subsubsection::Closed:: *)
(*GetKernel*)


GetKernel[waveVecs_,sc_]:=GetKernel[waveVecs,sc]=Block[{dGrid,v,pad,ker,win},
	(*define the kernel width*)
	v=Round[sc Norm[waveVecs]]-1;
	pad={v+1,-v-1};
	(*make and normalize the kernel*)
	ker=Table[HannWindow[(x/(2v+2))],{x,-v,v,1.}];
	ker=ConstantArray[ker/Total[ker],Length[ker]];
	(*ouput*)
	{ker*Transpose[ker],pad}
]


(* ::Subsubsection::Closed:: *)
(*GetWin*)


(*make edge fileter*)
GetWin[ker_,pad_,dim_]:=GetWin[ker,pad,dim]=Clip[2ListConvolve[ker,ConstantArray[1.,dim],pad,0.]-1,{0,1}]^2;


(* ::Subsubsection::Closed:: *)
(*BandFilter*)


BandFilter[dim_,waveVecs_,True]:=Map[Shift,BandFilter[dim,waveVecs],{2}]

BandFilter[dim_,waveVecs_]:=BandFilter[dim,waveVecs]=Block[{
	wx,wy,wxMat,wyMat,rMat,phiMat,dgrid,phi,rw,lnrw,lnrW,phiW,rBf,mask, BfRg,wrgMat,waveVec
	},
	
	(*define coordinate system in shifted kspace*)
	{wx,wy}=RotateLeft[N[Range[0,#-1]/#]-0.5,Floor[#/2]]&/@dim;
	{wxMat,wyMat}={Transpose@ConstantArray[wx,dim[[2]]],ConstantArray[wy,dim[[1]]]};
	
	(*get polar cordiantes*)
	{rMat,phiMat}=RotateDimensionsRight[Map[If[#=={0.,0.},#,ToPolarCoordinates[#]]&,N@RotateDimensionsLeft[{wxMat,wyMat}],{2}]];
	rMat[[1,1]]=10^-10.;
	
	(*make fileters*)
	Transpose[(
		waveVec=#;
		dgrid=Sqrt[Total[waveVec^2]];
		phi=ArcTan@@waveVec;
		rw=1/dgrid;
		lnrw=Log[rw];
		
		lnrW=Log[rMat]-lnrw;
		phiW=Mod[phiMat-phi+Pi,2 Pi]-Pi;
		
		rBf=Sqrt[lnrW^2+phiW^2];
		mask=UnitStep[1-rBf];
		
		BfRg=mask Cos[Pi/2*rBf]^2;
		wrgMat=mask Total[waveVec{wxMat,wyMat}];
				
		{mask,BfRg,wrgMat}
	)&/@waveVecs]
]


(* ::Subsection:: *)
(*Fourier*)


(* ::Subsubsection::Closed:: *)
(*Fourier*)


(*perfomr FFT and inverse FFT*)
FFT[im_]:=Fourier[im,FourierParameters->{1,-1}];
IFFT[im_]:=InverseFourier[im,FourierParameters->{1,-1}]


(* ::Subsubsection::Closed:: *)
(*Shift*)


Shift[img_]:=Block[{dx,dy},
	(*shift image by half field of view*)
	{dx,dy}=Round[Dimensions[img]/2];
	Transpose[RotateRight[Transpose[RotateRight[img,dx]],dy]]
]


(* ::Subsection::Closed:: *)
(*CalculateDispacementParameters*)


CalculateDispacementParameters[{motx_,moty_},mask_]:=Block[{
	v,kx,ky,xcor,ycor,rad,phi,sphi,cphi,ss,cc,cs2,
	du,ux,uy,v1,v2,nv1,nv2,v2n,v1n,dot,crs,rot,
	Fxx,Fxy,Fyx,Fyy,Exx,Eyy,Exy,Ecc,Ecr,out
	},
	(*denife kernel for derivative*)
	v=6;
	(*size of kernel*)
	kx=GaussianMatrix[{v},{1,0}];
	ky=GaussianMatrix[{v},{0,1}];
	
	(*calculate strains*)
	out=Table[
		{xcor,ycor,rad,phi}=MaskToCoordinates[mask[[i]]];
		
		(*define base vector and norm*)
		v1={xcor,ycor};
		nv1=rad;(*Sqrt[v1[[1]]^2+v1[[2]]^2];*)
		v1n=DevideNoZero[#,nv1]&/@v1;
		
		(*define angles for calculating Ecc and Ecr*)
		sphi = Sin[phi];cphi = Cos[phi];
		ss=sphi^2; cc=cphi^2; cs2 = 2 sphi cphi;
		
		Table[
			
			(*displacement and posision vectors*)
			du={ux,uy}={motx[[i,j]],moty[[i,j]]};
			v2=v1+du;
			
			(*calculate the norm and normalized displaced vector*)
			nv2=Sqrt[v2[[1]]^2+v2[[2]]^2];
			v2n=DevideNoZero[#,nv2]&/@v2;
			
			(*angle between two vectors = AcrTan[Dot[A,B],Cross[A,B]]*)
			dot=v1n[[1]]v2n[[1]]+v1n[[2]]v2n[[2]];(*dot*)
			crs=v1n[[1]]v2n[[2]]-v1n[[2]]v2n[[1]];(*cross*)
			
			rot=N@ArcTan[dot,crs]/Degree;
			
			(*calculate strain*)
			(*Displacement tensor G = F-I where F = I + (grad u)*)
			Fxx=ListConvolve[kx,ux,{v+1,-v-1},0];
			Fxy=ListConvolve[ky,ux,{v+1,-v-1},0];
			Fyx=ListConvolve[kx,uy,{v+1,-v-1},0];
			Fyy=ListConvolve[ky,uy,{v+1,-v-1},0];
			
			(*strain tensor E = 0.5 F`*F-I*)
			Exx = Fxx + .5(Fxx^2 + Fyx^2);
			Eyy = Fyy + .5(Fyy^2 + Fxy^2);
			Exy =0.5( Fxy + Fyx + Fxx Fxy + Fyx Fyy);
			
			(*get the circular and radial strain components R`.E.R*)
			Ecc=ss Exx + cc Eyy - cs2 Exy;
			Ecr=cc Exx + ss Eyy + cs2 Exy;
			
			(*output*)
			{Exx,Eyy,Exy,Ecc,Ecr,rot}
		,{j,1,Length[motx[[1]]],1}]
	,{i,1,Length[motx],1}];
	
	(*make parameters first dimension*)
	Transpose[out,{2,3,1,4,5}]
]


(* ::Subsection::Closed:: *)
(*TaggingParPlot*)


col[n_]:=Blend[{Darker[Darker[Red]],Blend[{Yellow,Darker[Yellow]},.5]},#]&/@(Range[0,n]/(n-1))

TaggingParPlot[dati_,lab_,rani_:0]:=Block[{dat,l,d,style,leg,r,ran},
	dat=Transpose@dati;
	{l,d}=Dimensions@dat;
	
	style=Directive[{Thickness[.02],#}]&/@col[l];
	leg="Slice "<>ToString[#]&/@Range[l];
	
	ran=If[rani===0,r=1.2Max@Abs[dat];{-r,r},rani];
	
	ListLinePlot[dat,PlotStyle->style,PlotRange->{{-3,d+2},ran},AxesOrigin->{1,0},Ticks->False,
		Frame->{{True,False},{True,False}},FrameStyle->Directive[{Thick,Black}],
		AspectRatio->0.5,PlotLegends->leg, ImageSize->300,
		PlotLabel->Style[lab,Bold,12,Black]
	]
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
