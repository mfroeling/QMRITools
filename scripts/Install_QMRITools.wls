#!/usr/bin/env wolframscript
(* ::Package:: *)

inst = PacletInstall[Last@$ScriptCommandLine];

Print[inst];

If[inst === $Failed,
	Print["Installation Failed"]
	,
	load = Get["QMRITools`"];
	pack = Last@PacletFind["QMRITools"];
	Print["Installed ",pack["Location"]];
]
