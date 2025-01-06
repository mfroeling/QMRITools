#!/usr/bin/env wolframscript
(* ::Package:: *)

inst = If[Length[$ScriptCommandLine]===1,
	ResourceFunction["GitHubInstall"]["mfroeling", "QMRITools"],
	PacletInstall[Last@$ScriptCommandLine]
];

If[inst === $Failed,
	Print["Installation Failed"]
	,
	load = Get["QMRITools`"];
	pack = Last@PacletFind["QMRITools"];
	Print["Installed ",pack["Location"]];
]
