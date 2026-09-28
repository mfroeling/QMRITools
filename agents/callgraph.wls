#!/usr/bin/env wolframscript
(* Call graph of QMRITools/Kernel via the CodeGraph paclet, installed or cloned next to this repo.
   Build:  wolframscript -file agents/callgraph.wls
           writes agents/callgraph/overview.md, defs.tsv and edges.tsv and prints a summary
   Query:  wolframscript -file agents/callgraph.wls callers|callees <symbol> [depth] *)

root = DirectoryName[$InputFileName];
If[PacletFind["CodeGraph"] === {}, PacletDirectoryLoad[FileNameJoin[{ParentDirectory[root, 2], "CodeGraph", "CodeGraph"}]]];
Needs["CodeGraph`"];

out = FileNameJoin[{root, "callgraph"}];
args = Rest[$ScriptCommandLine];

If[Length[args] >= 2 && MemberQ[{"callers", "callees"}, args[[1]]],
	Print[CodeGraphQuery[ImportCodeGraph[out], args[[1]], args[[2]], If[Length[args] > 2, FromDigits[args[[3]]], 1]]],

	cg = BuildCodeGraph[FileNameJoin[{ParentDirectory[root], "QMRITools", "Kernel"}]];
	ExportCodeGraph[out, cg];
	KeyValueMap[Print[#1, ": ", #2] &, CodeGraphSummary[cg]]
];
