#!/usr/bin/env wolframscript
(* Compares each documentation guide with the public functions of its toolbox.
   Run:  wolframscript -file agents/guidecheck.wls [Toolbox ..]
   Per guide it prints functions missing from the guide, links to symbols that no longer exist (stale, or renamed),
   duplicate links and functions without a reference page. Options are ignored. Loads the working tree in a fresh kernel. *)

root = ParentDirectory[DirectoryName[$InputFileName]];
PacletDirectoryLoad[FileNameJoin[{root, "QMRITools"}]];
Get["QMRITools`"];

dir = FileNameJoin[{root, "QMRITools", "Documentation", "English"}];
guides = Complement[FileBaseName /@ FileNames["*.nb", FileNameJoin[{dir, "Guides"}]], {"QMRITools", "ScientificColorData"}];
If[Length[Rest[$ScriptCommandLine]] > 0, guides = Intersection[guides, Rest[$ScriptCommandLine]]];

(*options that are not listed by Options[..] are recognised by their usage text*)
optionQ[name_String, context_String] := StringMatchQ[
	Quiet@ToString@ToExpression[context <> name <> "::usage"], __ ~~ (" is an option" | " is an options") ~~ __];

Do[
	r = Quiet@QMRIToolsFunctions[g];
	If[! MatchQ[r, {{_, _List, _List}}], Print[g, ": no toolbox"]; Continue[]];
	funcs = Select[StringReplace[r[[1, 2]], "QMRITools`" ~~ __ ~~ "`" -> ""], ! optionQ[#, "QMRITools`" <> g <> "`"] &];
	txt = Import[FileNameJoin[{dir, "Guides", g <> ".nb"}], "Text"];
	links = StringCases[txt, "paclet:QMRITools/ref/" ~~ n : (("$" | LetterCharacter) ~~ WordCharacter ...) :> n];
	dups = Keys@Select[Counts[links], # > 1 &];
	noRef = Select[funcs, ! FileExistsQ[FileNameJoin[{dir, "ReferencePages", "Symbols", # <> ".nb"}]] &];
	Print["== ", g, " (", Length[funcs], " functions, ", Length@Union[links], " links)"];
	Print["  missing: ", Complement[funcs, links]];
	Print["  stale:   ", Complement[Union[links], funcs]];
	If[dups =!= {}, Print["  dups:    ", dups]];
	If[noRef =!= {}, Print["  no ref page: ", noRef]],
	{g, guides}];
