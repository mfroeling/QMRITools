#!/usr/bin/env wolframscript
(* ::Package:: *)

(*load QMRITools, the dev version when it is installed*)
Print["\nLoading QMRITools"];
dev = PacletFind["QMRIToolsDev"] =!= {};
t1 = First@AbsoluteTiming[load = Get[If[dev, "QMRIToolsDev`", "QMRITools`"]]];
Print[">> Loading time [s]: ", ToString[Round[t1, .1]]];

If[load === $Failed,
	Print["\nCould not load QMRITools"];
	Quit[1]
];

pack = SelectFirst[PacletFind["QMRITools"], StringStartsQ[FindFile["QMRITools`"], #["Location"], IgnoreCase->True]&];
Print[">> ", "Using QMRITools " <> pack["Version"] <> " at " <> pack["Location"] <> If[dev, " (dev)", ""]];


(*steps in pipeline order and the report steps that check them*)
steps = <|
	"DcmToNii" -> BidsDcmToNii,
	"Convert" -> MuscleBidsConvert,
	"Process" -> MuscleBidsProcess,
	"Merge" -> MuscleBidsMerge,
	"Segment" -> MuscleBidsSegment,
	"Tractography" -> MuscleBidsTractography,
	"Analysis" -> MuscleBidsAnalysis
|>;

reportSteps = <|
	"DcmToNii" -> {"Protocol names", "Config labels"},
	"Convert" -> {"Convert", "Custom configs"},
	"Process" -> {"Process"},
	"Merge" -> {"Merge"},
	"Segment" -> {"Segment"},
	"Tractography" -> {"Tractography"},
	"Analysis" -> {"Analysis"}
|>;

(*report steps of the check files, there missing after a run means a failed dataset*)
checkSteps = {"Process", "Merge", "Segment", "Tractography", "Analysis"};

(*values of the options with a fixed list*)
methods = {"Full", "Tractography", "Segmentation"};
images = {"All", "Quantitative", "Segmentation", "Tractography", "None"};

(*flags with a value, flags without a value and flags with an optional value*)
valueFlags = {"--dir", "--step", "--subs", "--method", "--images", "--report"};
switchFlags = {"--versioncheck", "--keepraw", "--check", "--h", "--help"};
optionalFlags = {"--debug"};


(*the help text*)
help[] := Print[StringRiffle[{
	"",
	"Usage:",
	"  wolframscript -script MuscleBids.wls --dir <study> --step <steps> [options]",
	"  wolframscript -script MuscleBids.wls --dir <study> --check | --report <file>",
	"",
	"Flags:",
	"  --dir <folder>       Study folder with the config.json (required)",
	"  --step <steps>       Steps to run, space or comma separated, or All:",
	"                       " <> StringRiffle[Keys[steps], " "],
	"                       Steps always run in pipeline order",
	"  --subs <names>       Subjects to process, default All. Use the dicom folder name (001_1)",
	"                       or sub-001_ses-1. A name without session selects session 001",
	"  --versioncheck       Redo data that was processed with an older QMRITools version",
	"  --keepraw            Keep the raw nii of dcm2niix after Convert. Convert has no check file,",
	"                       so with --keepraw every run converts these files again",
	"  --method <method>    Tractography steps: " <> StringRiffle[methods, ", "] <> " (default Full)",
	"  --images <images>    Analysis images: " <> StringRiffle[images, ", "] <> " (default All)",
	"  --check              Only show the problem summary of all steps, nothing is run",
	"  --report <file>      Write the problem report of all steps as .md or .pdf, after the steps or alone",
	"  --debug [tags]       Print debug output, default tag Bids",
	"  --h, --help          Show this help message",
	"",
	"Rerun after a crash:",
	"  Process to Analysis only redo the datasets that did not finish.",
	"  DcmToNii skips a session when its raw folder is not empty, empty that folder to convert it again.",
	"",
	"Exit codes:",
	"  0  done without problems",
	"  1  invalid input, nothing was run",
	"  2  a step failed, the later steps were not run",
	"  3  done, but with messages or with problems in the steps that ran"
}, "\n"]]

(*print an error and stop*)
fail[msg__] := (
	Print["\nError: ", msg];
	Quit[1]
)

(*the item of list that equals val, ignoring case*)
match[val_, list_] := SelectFirst[list, ToLowerCase[#] === ToLowerCase[val]&, Missing[]]

(*split a flag value on spaces and commas*)
split[val_] := StringSplit[val, {" ", ","}]

(*summary line per report step*)
summaryLine[sum_, st_] := Block[{counts, text},
	counts = Lookup[sum["Counts"], st, <||>];
	text = Which[
		!MemberQ[sum["Checked"], st], "not checked",
		counts === <||>, "ok",
		True, StringRiffle[KeyValueMap[ToString[#2] <> " " <> #1&, counts], ", "]
	];
	Print["  ", StringPadRight[st, 16], text]
]


(*-------------- Arguments --------------*)

args = ParseCommandLine[Rest@$ScriptCommandLine];

If[KeyExistsQ[args, "--h"] || KeyExistsQ[args, "--help"],
	help[];
	Quit[0]
];

(*the new public functions of the summary must exist*)
If[Names["QMRITools`*`MuscleBidsSummary"] === {} || Names["QMRITools`*`CheckBidsSubjects"] === {},
	fail["QMRITools ", pack["Version"], " is too old for this script, use version 5.3.1 or newer."]
];

(*unknown flags and flags with a wrong value*)
bad = Complement[Keys[args], valueFlags, switchFlags, optionalFlags];
If[bad =!= {},
	fail["Unknown flag ", StringRiffle[bad, ", "], ", see --help."]
];
bad = Select[valueFlags, KeyExistsQ[args, #] && !StringQ[args[#]]&];
If[bad =!= {},
	fail["The flag ", StringRiffle[bad, ", "], " needs a value."]
];
bad = Select[switchFlags, KeyExistsQ[args, #] && args[#] =!= True&];
If[bad =!= {},
	fail["The flag ", StringRiffle[bad, ", "], " takes no value."]
];

(*there must be something to do*)
If[!KeyExistsQ[args, "--step"] && !KeyExistsQ[args, "--check"] && !KeyExistsQ[args, "--report"],
	help[];
	Quit[1]
];

(*study folder with a readable config*)
If[!KeyExistsQ[args, "--dir"],
	fail["No study folder given, use --dir."]
];
dir = ExpandFileName[args["--dir"]];
If[!DirectoryQ[dir],
	fail["The study folder ", dir, " does not exist."]
];
If[!AssociationQ[Quiet@GetConfig[dir]],
	fail["No readable config.json in ", dir, "."]
];

(*steps to run in pipeline order*)
run = If[KeyExistsQ[args, "--step"],
	split[args["--step"]],
	{}
];
bad = Select[run, MissingQ[match[#, Append[Keys[steps], "All"]]]&];
If[bad =!= {},
	fail["Unknown step ", StringRiffle[bad, ", "], ", use ", StringRiffle[Keys[steps], ", "], " or All."]
];
run = If[MemberQ[ToLowerCase[run], "all"],
	Keys[steps],
	Select[Keys[steps], MemberQ[ToLowerCase[run], ToLowerCase[#]]&]
];

(*options with a fixed list of values*)
method = match[Lookup[args, "--method", "Full"], methods];
If[MissingQ[method],
	fail["Unknown tractography method ", args["--method"], ", use ", StringRiffle[methods, ", "], "."]
];
imOut = match[Lookup[args, "--images", "All"], images];
If[MissingQ[imOut],
	fail["Unknown analysis images ", args["--images"], ", use ", StringRiffle[images, ", "], "."]
];

(*report file*)
report = Lookup[args, "--report", None];
If[StringQ[report] && !MemberQ[{"md", "pdf"}, ToLowerCase[FileExtension[report]]],
	fail["The report file ", report, " must be .md or .pdf."]
];

(*subjects, every name must select a session*)
subs = If[!KeyExistsQ[args, "--subs"] || ToLowerCase[args["--subs"]] === "all",
	All,
	split[args["--subs"]]
];
sessions = CheckBidsSubjects[dir, subs];
bad = Keys@Select[sessions, # === {}&];
If[subs =!= All && bad =!= {},
	fail["No session found for ", StringRiffle[bad, ", "], ". Use the dicom folder name (001_1) or sub-001_ses-1, a name without session selects session 001."]
];
If[Union @@ Values[sessions] === {},
	fail["No sessions found in the dicom or raw data folder of ", dir, "."]
];

(*debug output*)
If[KeyExistsQ[args, "--debug"],
	QMRITools`$Debug = If[StringQ[args["--debug"]],
		split[args["--debug"]],
		"Bids"
	]
];

(*options of the steps, each step gets only its own*)
opts = {
	ProcessSubjects -> subs,
	VersionCheck -> KeyExistsQ[args, "--versioncheck"],
	DeleteAfterConversion -> !KeyExistsQ[args, "--keepraw"],
	BidsTractographyMethod -> method,
	BidsOutputImages -> imOut
};

Print["\nStudy folder:\n>> ", dir];
Print["Sessions:\n>> ", Length[Union @@ Values[sessions]]];


(*-------------- Run the steps --------------*)

status = <||>;
Do[
	Print["\n==================== ", st, " ===================="];
	t0 = AbsoluteTime[];
	msgs = False;
	(*run the step, messages are flagged and an abort fails the step*)
	res = CheckAbort[Check[
		Quiet[steps[st][dir, Sequence @@ FilterRules[opts, Options[steps[st]]]], FrontEndObject::notavail],
		msgs = True
	], $Aborted];
	(*status of the step*)
	status[st] = Which[
		res === $Aborted || res === $Failed, "failed",
		msgs, "done with messages",
		True, "done"
	];
	Print["\n>> ", st, ": ", status[st], " (", Round[AbsoluteTime[] - t0], " s)"];
	(*later steps read the output of this step*)
	If[status[st] === "failed",
		Break[]
	];
, {st, run}];


(*-------------- Summary and report --------------*)

code = 0;
If[run =!= {} || KeyExistsQ[args, "--check"],
	Print["\nProblem summary:"];
	sum = MuscleBidsSummary[dir, "All", ProcessSubjects -> subs, BidsOutputImages -> imOut];
	(*the report steps of the steps that ran, or all with --check*)
	shown = If[run === {},
		Flatten[Values[reportSteps]],
		Flatten[Lookup[reportSteps, Keys[status]]]
	];
	summaryLine[sum, #]& /@ shown;
	(*problems in the check files of the steps that ran, outdated only counts with --versioncheck*)
	classes = If[KeyExistsQ[args, "--versioncheck"],
		{"error", "unreadable", "missing", "outdated"},
		{"error", "unreadable", "missing"}
	];
	problems = If[run === {},
		0,
		Total[Flatten[Lookup[Lookup[sum["Counts"], #, <||>], classes, 0]& /@ Intersection[shown, checkSteps]]]
	];
	If[problems > 0,
		Print["Problems found, details with --report <file.md>"]
	];
	code = Which[
		MemberQ[Values[status], "failed"], 2,
		MemberQ[Values[status], "done with messages"] || problems > 0, 3,
		True, 0
	];
];

If[StringQ[report],
	Print["\nWriting the report"];
	file = MuscleBidsReport[dir, "All", report, ProcessSubjects -> subs, BidsOutputImages -> imOut];
	If[file === $Failed,
		Print[">> The report could not be written"],
		Print[">> ", file]
	];
];

Quit[code]
