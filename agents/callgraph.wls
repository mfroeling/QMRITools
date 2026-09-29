#!/usr/bin/env wolframscript
(* Call graph of QMRITools/Kernel via the CodeGraph paclet, cloned next to this repo or installed.
   Build:  wolframscript -file agents/callgraph.wls
           writes agents/callgraph/overview.md, defs.tsv and edges.tsv and prints a summary
   Query:  wolframscript -file agents/callgraph.wls callers|callees <symbol> [depth] *)

root = DirectoryName[$InputFileName];
clone = FileNameJoin[{ParentDirectory[root, 2], "CodeGraph", "CodeGraph", "Kernel", "CodeGraph.wl"}];
(*a clone next to this repo is the development version and wins over an installed release*)
If[FileExistsQ[clone], Get[clone], Needs["CodeGraph`"]];

out = FileNameJoin[{root, "callgraph"}];
args = Rest[$ScriptCommandLine];

If[Length[args] >= 2 && MemberQ[{"callers", "callees"}, args[[1]]],
	Print[CodeGraphQuery[ImportCodeGraph[out], args[[1]], args[[2]], If[Length[args] > 2, FromDigits[args[[3]]], 1]]],

	(*Legacy is unmaintained and ScientificColorData only loads color tables*)
	cg = BuildCodeGraph[FileNameJoin[{ParentDirectory[root], "QMRITools", "Kernel"}], "Exclude" -> {"Legacy", "ScientificColorData"}];
	ExportCodeGraph[out, cg];
	KeyValueMap[Print[#1, ": ", #2] &, CodeGraphSummary[cg]]
];
