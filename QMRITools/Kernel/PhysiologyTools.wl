(* ::Package:: *)

(* ::Title:: *)
(*QMRITools PhysiologyTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`PhysiologyTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`PhysiologyTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


AlignRespLog::usage =
"AlignRespLog[physLog, respirect, scanTime] aligns respirect and physlog data. physLog is output from ImportPhyslog.
resirect is the first output from ImportRespirect.";

ImportPhyslog::usage =
"ImportPhyslog[] imports all physlog files from the folder selcted.
ImportPhyslog[\"forder\"] imports all physlog files from \"folder\" selcted."

ImportRespirect::usage = 
"ImportRespirect[] impors all the respirect log files from the folder selcted.
ImportRespirect[\"folder\"] impors all the respirect log files from the \"folder\" selcted."

PlotPhyslog::usage =
"PlotPhyslog[{time, resp}, {start, stop}] plots the physlog from ImportPhyslog.
PlotPhyslog[{time, resp}, {start, stop}, scanTime] plots the physlog from ImportPhyslog."

PlotRespiract::usage =
"PlotRespiract[data, dataP, scantimes] plots the respirect data to correct peaks. data and dataP are the first outputs of ImportResirect. scantimes is the output from AlignRespLog. 
PlotRespiract[data, dataP, scantimes, steps]."


(* ::Subsection::Closed:: *)
(*Options*)


OutputMethod::usage = "OutputMethod can be \"val\" or \"plot\"."

SampleStep::usage= "SampleStep is an option for AlignRespiract."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] (* Begin Private Context *) 


(* ::Subsection::Closed:: *)
(*AlignRespLog*)


Options[AlignRespLog] = Options[AlignRespLogi] = {OutputMethod -> "val", SampleStep -> 0.005};

SyntaxInformation[AlignRespLog] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

AlignRespLog[physLog_, respirect_, scanTime_, opts:OptionsPattern[]]:=AlignRespLog[physLog, respirect, scanTime, Range[Length[respirect]], opts]

AlignRespLog[physLog_, respirect_, scanTime_, order_, opts:OptionsPattern[]]:= Block[{n},
	Switch[OptionValue[OutputMethod],
	"val",
	Transpose[AlignRespLogi[physLog, respirect[[#]], scanTime[[#]], opts] & /@ order],
	"plot",
	Manipulate[AlignRespLogi[physLog, respirect[[order[[n]]]], scanTime[[order[[n]]]], opts] ,{{n,1,"dataset"},1,Length[order],1}]
]
];

AlignRespLogi[physLog_, respirect_, scanTime_, OptionsPattern[]] := 
 Block[{ptime, pdat, start, stop, rtime, rdat, rstart, rend, pstart, co2dataSel,
   pend, pran, samp, rtimei, rint, ptimei, pint, corr, poff, len,corrAll,n,nsel,
   startscan, stopscan},
  
  co2dataSel=respirect[[1]];
  
  {rtime, rdat} = co2dataSel[[All, {1, 4}]] // Transpose;
  samp = OptionValue[SampleStep];
  len = Length[physLog];
  
  corrAll = {corr, rstart, samp, pint, ptime, pdat, stop} = Transpose[
     (
        n = #;
        {ptime, pdat} = physLog[[n, 1]];
        {start, stop} = physLog[[n, 2]];
        
        {rstart, rend} = StartEnd[rtime];
        {pstart, pend} = StartEnd[ptime];
        
        pran = pend - pstart;
        
        rtimei = Range[rstart, rend, samp];
        rint = Interpolation[{rtime, rdat} // Transpose][rtimei];
        ptimei = Range[pstart, pend, samp] + rstart;
        pint = 
         Interpolation[{ptime + rstart, pdat} // Transpose][ptimei];
        corr = Abs[ListCorrelate[pint, rint, {-1, 1}, 0]];
        
        {corr, rstart, samp, pint, ptime, pdat, stop}
        ) & /@ Range[len]];
  
  nsel = First[Flatten[Position[corrAll[[1]], Max[corrAll[[1]]]]]];
  {corr, rstart, samp, pint, ptime, pdat, stop} = corrAll[[All, nsel]];
  poff = rstart + samp (First@Flatten[Position[corr, Max[corr]]] - Length[pint]);
  {startscan, stopscan} = {stop - scanTime, stop} + poff;
  
  Switch[OptionValue[OutputMethod],
   "val", {{startscan, stopscan},physLog[[nsel]],respirect},
   "plot", ListLinePlot[
    {{rtime, 1 - (rdat - Min[rdat])/(Max[rdat] - Min[rdat])} // Transpose,
     {ptime + poff, (pdat - Min[pdat])/(Max[pdat] - Min[pdat])} // Transpose
     },
    AspectRatio -> 0.2, PlotRange -> {Full, Full}, ImageSize -> 1000, 
    PlotStyle -> {Red, Black}, 
    GridLines -> {{{startscan, 
        Directive[{Thickness[.005], Green}]}, {stopscan, 
        Directive[{Thickness[.005], Red}]}}, None}]]
  ]

StartEnd = {First[#], Last[#]} &;


(* ::Subsection::Closed:: *)
(*AlignRespLog*)


SyntaxInformation[ImportPhyslog] = {"ArgumentsPattern" -> {_.}};

ImportPhyslog[] := ImportPhyslog[FileSelect["Directory" ,WindowTitle->"Select directory that contains the physlogs"]]
ImportPhyslog[folder_] := Block[
	{files, file, resp, mark, time, start, stop, sel},
  files = Sort[FileNames["*.log", {folder}]];
  (
     {resp, mark} = 
      Drop[ToExpression[
            "{" <> StringReplace[#, " " .. -> ","] <> "}"] & /@ 
          DeleteCases[Import[files[[#]], "Lines"][[6 ;;]], 
           "#"], -1][[All, {6, 10}]] // Transpose;
     time = Range[0, Length[resp] - 1] (1/500/60.);
     
     start = time[[First[Flatten[Position[mark, 10]]]]];
     
     stop = Flatten[Position[mark, 20]];
     stop = time[[If[stop == {}, Length[mark], Last[stop]]]];
     
     sel = Range[1, Length[time], 50];
     
     {{time, resp}[[All, sel]], {start, stop}}
     
     ) & /@ Range[Length[files]]
  ]


(* ::Subsection::Closed:: *)
(*PlotPhyslog*)


SyntaxInformation[PlotPhyslog] = {"ArgumentsPattern" -> {_, _, _.}};

PlotPhyslog[{time_, resp_}, {start_, stop_}, scanTime__: 0] := 
 Block[{start2},
  start2 = If[scanTime == 0, start, stop - scanTime];
  ListLinePlot[Transpose[{time, resp}], GridLines -> {
     {{start, Directive[{Dashed, Thickness[0.005], Green}]}, {start2, 
       Directive[{Thickness[0.005], Green}]}, {stop, 
       Directive[{Thickness[0.005], Red}]}}
     , None}, PlotStyle -> Black, AspectRatio -> 0.1,
     Axes->False,Frame->{{True,False},{True,False}}, 
   ImageSize -> 1000, PlotLabel -> {stop - start, stop - start2}]
  ]


(* ::Subsection::Closed:: *)
(*ImportRespirect*)


SyntaxInformation[ImportRespirect] = {"ArgumentsPattern" -> {_.}};

ImportRespirect[] := ImportRespirect[FileSelect["Directory",WindowTitle->"Select directory that contains the respiract data"]]
ImportRespirect[folder_] := 
 Module[{co2data, co2dataP, events, events2, checks, xsel, timesel, sel, start, end, a, b},

  If[folder === Null, 
  	
  	  Return[],
  
	  co2dataP = (Import[FileNames[{"bbb_*.txt"}, {folder}][[1]], "Data"])[[2 ;;, 1 ;; 3]];
	  co2data = (DeleteCases[Import[FileNames[{"raw_*.txt"}, {folder}][[1]],"Data"], {""}])[[3 ;;, {1, 4, 5, 2}]];
	  events = (Import[FileNames[{"events_*.txt"}, {folder}][[1]], "Data"])[[2 ;;]];
	  
	  events2 = Delete[events, Position[events, "END sequence"][[All, {1}]]];
	  
	  DynamicModule[{
	    x = ConstantArray[False, Length[events2]],
	    time = 6},
	   
	   checks = {Checkbox[Dynamic[x[[#]]]]} & /@ Range[Length[events2]];
	   
	   {xsel, timesel} = DialogInput[
	     {
	      TextCell["Choose events: "],
	      Row[{Thread[{checks, events2[[All, 2]]}] // TableForm}],
	      TextCell["Enter the respirect experiment duration: "],
	      InputField[Dynamic[time], Number],
	      DefaultButton[DialogReturn[{x, time}]]
	      }, Modal -> True];
	   ];
	  
	  If[Total[Boole[xsel]] == 0,
	   
	   Return[Print["Error, more than one event need to be selected"]],
	   
	   sel = Flatten[Position[xsel, True]];
	   
	   ((
	        start = #;
	        end = start + timesel;
	        a = Select[co2data, end > #[[1]] > start &];
	        b = Select[co2dataP, end > #[[1]] > start &];
	        {a, b}
	        ) & /@ events2[[sel, 1]]) 
	   ]
   ]
  ]


(* ::Subsection:: *)
(*ImportRespirect*)


(* ::Subsubsection::Closed:: *)
(*ImportRespirect*)


SyntaxInformation[PlotRespiract] = {"ArgumentsPattern" -> {_, _, _.}};

PlotRespiract[dataAll_, scantimes_, steps__: 10] := PlotRespiracti[dataAll, scantimes, steps][[2]]

PlotRespiracti[dataAll_, scantimes_, steps_] := DynamicModule[
  {datas, samplesO2, allx, samplesCO2, pos, CO2, O2, pt, pd, peak,data,dataP,
   allCO2, allO2,  max, min, scant, r, but, temp,range, allPmouth, samplesx,
    rr,  CO2tot, O2tot, dataout, len, shift, n, pointsall, x,cent, span, pmin, pmax, pran},
  
  {data, dataP}=Transpose[dataAll];
  
  range = Range[1, Length[data[[1, All, 1]]], steps];
  len = Length[data];
  shift = {0, 1.2, 0, 0};
  
  CO2tot = O2tot = ConstantArray[0, len];
  dataout = Transpose /@ dataP;
  
  datas = Transpose /@ data[[All, range]];
  
  {min, max} = 
   Transpose[Transpose[{Min[#], Max[#]} & /@ #] & /@ datas];
  rr = max - min; rr[[All, 1]] = 1;
  
  {allx, allCO2, allO2, allPmouth} = Transpose[# + shift & /@ ((datas - min)/rr)];
  
  {samplesx, samplesCO2, samplesO2} = 
   Transpose[# + shift[[1 ;; 3]] & /@ ((dataout - min[[All, 1 ;; 3]])/
       rr[[All, 1 ;; 3]])];
  
  scant = scantimes - min[[All, 1]];
  pointsall = MapThread[Transpose[{#1, #2}] &, {samplesx, samplesCO2}];
  temp = pointsall = MapThread[(r = #2; Select[#1, r[[1]] < #[[1]] < r[[2]] &]) &, {pointsall, scant}];
  n = 1; but = False;
  
  {pmin,pmax}={Min[allx],Max[allx]};
  pran=pmax-pmin;
  
  DialogInput[{
     Manipulate[
     	cent = Clip[cent,{pmin+span,pmax-span}];
	    LocatorPane[
	     Dynamic[pointsall[[n]]],
	     Dynamic[
	      n = Clip[n, {1, len}, {1, len}];
	      
	      pos = 
	       Position[allx[[n]], First[Nearest[allx[[n]], #, 1]]][[1, 
	           1]] & /@ pointsall[[n, All, 1]];
	      
	      x = allx[[n, pos]];
	      CO2 = pointsall[[n, All, 2]] = allCO2[[n, pos]];
	      O2 = allO2[[n, pos]];
	      
	      dataout[[n]] = ({x, CO2 - shift[[2]], O2} rr[[n, 1 ;; 3]]) + min[[n, 1 ;; 3]];
	      
	      CO2tot[[n]] = Round[Mean[dataout[[n, 2]]], .1];
	      O2tot[[n]] = Round[Mean[dataout[[n, 3]]], .1];
	      
	     
	      	Graphics[
		       {
		        Thickness[.0015],
		        Darker[Darker[Red]], Line[Transpose[{allx[[n]], allCO2[[n]]}]],
		        Darker[Darker[Blue]], Line[Transpose[{allx[[n]], allO2[[n]]}]]
		        ,
		        PointSize[Large],
		        Red, Point[Transpose[{x, CO2}]],
		        Blue, Point[Transpose[{x, O2}]]
		        ,
		        Thickness[.001],
		        Red, Line[SortBy[Transpose[{x, CO2}], First]],
		        Blue, Line[SortBy[Transpose[{x, O2}], First]],
		        
		        Opacity[.5], Thickness[.0025], Green,
		        Line[{{scant[[n, 1]], -0.1}, {scant[[n, 1]], 2.3}}],
		        Red,
		        Line[{{scant[[n, 2]], -0.1}, {scant[[n, 2]], 2.3}}]
		        },
		        
		       Background -> White,
		       PlotRange->{{cent-span,cent+span},{-0.1,2.4}},
		       ImageSize -> 1200, AspectRatio -> 0.4, Frame -> True, 
		       FrameTicks -> {True, False}, LabelStyle -> Large,
		       GridLines->{x,None},
		       LabelStyle->Black,
		       PlotLabel -> Column[{
		       	Style["Dataset "<>ToString[n],Bold,Large],
		       	Style[LabFun[x, CO2tot[[n]], O2tot[[n]]], Bold, Medium], 
		        TableForm[{CO2tot, O2tot}, TableHeadings -> {{"CO2", "O2"}, Range[len]}]}, Alignment -> Center]
		       ]
		   
		   ], LocatorAutoCreate -> True, Appearance -> None]
	   
      ,{{cent,pran/2,"Center"},Dynamic[pmin+span],Dynamic[pmax-span],.1}
	  ,{{span,pran/2,"Span"},1,pran/2,.1}
	  ]
    ,
    Dynamic[Grid[{
       {
        Button["   Back   ", If[n > 1, n--]],
        Button["   Next   ", If[n < len, n++]],
        CancelButton[]
        
        },
       {
        Button["Peak Detect", (
          temp[[n]] = pointsall[[n]];
          {pt, pd} = 
           Transpose[
            Select[Transpose[{allx[[n]], allCO2[[n]]}], 
             scant[[n, 1]] < #1[[1]] < scant[[n, 2]] &]];
          
          peak = (Partition[
              PeakDetect[pd/Mean[pd], 20/steps, 0, 0.5, Padding -> 10], 
              4, 1, 4] /. {{0, 0, 0, 1} -> 1, {_, _, _, _} -> 0});
          
          pos = Position[allx[[n]], #][[1, 1]] & /@ (pt = 
              DeleteCases[(peak pt), 0.]);
          pointsall[[n]] = Transpose[{pt, allCO2[[n, pos]]}];
          )]
        ,
        Button["    Undo    ", pointsall[[n]] = temp[[n]]],
        If[n == len || but, but = True; 
         DefaultButton[
          DialogReturn[{CO2tot, O2tot, Transpose[{data, SortBy[Transpose[#],1]& /@ dataout}]}]]]
        }
       }]]
    }, Modal -> True, WindowFloating -> True]
  ]


(* ::Subsubsection::Closed:: *)
(*LabFun*)


LabFun[x_, co2_, o2_] := Block[{st, en},
  StringJoin@(ToString /@ {
      "CO2:  ", co2, "   ",
      "O2:  ", o2, "\n",
      "start:  ", st = Round[Min[x], .1], "  ",
      "end:  ", en = Round[Max[x], .1], "  ",
      "span:  ", en - st, "  "
      })]


(* ::Section:: *)
(*End Package*)


End[] (* End Private Context *)

EndPackage[]
