#!/usr/bin/env wolframscript
(* ::Package:: *)

inst = If[Length[$ScriptCommandLine]===1,
	ResourceFunction["GitHubInstall"]["mfroeling", "QMRITools"],
	PacletInstall[Last@$ScriptCommandLine]
];

If[!PacletObjectQ[inst], Print["Installation Failed"]; Exit[1]];

Print["Installed versions: ", StringRiffle[#["Version"]& /@ PacletFind["QMRITools"], ", "]];
Print["Just installed: ", inst["Version"], " at ", inst["Location"]];

If[Get["QMRITools`"]===$Failed, Print["Installed but failed to load"]; Exit[1]];
