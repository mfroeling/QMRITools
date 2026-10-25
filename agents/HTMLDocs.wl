(*HTML documentation build for the paclet, run from a notebook or script:
    Get["<repo>/agents/HTMLDocs.wl"];
    BuildHTMLDocs[location, builDirhtml];   (*only new or changed notebooks, continues after a crash*)
    FixHTMLDocs[builDirhtml];               (*repairs the pages, can be run again*)
  Written for Wolfram 15.0, where DocumentationBuildHTML links assets that its own bundle does not contain. Use forward
  slashes in the paths. The online assets are downloaded once into a folder next to the build folder, not into the repo.*)

Needs["DocumentationBuild`"];

Options[BuildHTMLDocs] = {"Batch" -> 25, "Force" -> False, "Only" -> All};

(*Builds the HTML of every notebook in paclet/Documentation/English that has no up to date html in out/html.
  Rerunning after a crash or a stop continues where it left off, because done is judged from the files on disk.*)
BuildHTMLDocs[paclet_String, out_String, OptionsPattern[]] := Module[{nbs, rel, htmlOf, todo, failed = {}, n = 0, res},
	nbs = FileNames["*.nb", FileNameJoin[{paclet, "Documentation", "English"}], Infinity];
	If[OptionValue["Only"] =!= All, nbs = Select[nbs, MemberQ[OptionValue["Only"], FileBaseName[#]] &]];

	(*html location in the build output: Guides -> guide, ReferencePages/Symbols -> ref*)
	htmlOf[nb_] := FileNameJoin[{out, "html", If[StringContainsQ[nb, "Guides"], "guide", "ref"], FileBaseName[nb] <> ".html"}];

	todo = If[OptionValue["Force"], nbs,
		Select[nbs, ! FileExistsQ[htmlOf[#]] || FileDate[htmlOf[#]] < FileDate[#] &]];
	Print["HTML files to build: ", Length[todo], " of ", Length[nbs]];

	Do[
		(*the filter gets paths relative to the paclet, so match on exactly those*)
		rel = StringDrop[#, StringLength[paclet] + 1] & /@ batch;
		res = Quiet@DocumentationBuild`DocumentationBuildHTML[paclet, out, MemberQ[rel, #] &,
			"UseLanguageExtension" -> False, "CompleteHTMLQ" -> True];
		failed = Join[failed, If[AssociationQ[res["ProcessedFiles"]], Keys@Select[res["ProcessedFiles"], FailureQ], batch]];
		n += Length[batch];
		Print[n, " / ", Length[todo], "  failed so far: ", Length[failed]],
		{batch, Partition[todo, UpTo[OptionValue["Batch"]]]}];

	If[failed =!= {}, Print["Failed (rerun to retry): ", failed]];
	failed
]


(*Repairs the html written by DocumentationBuildHTML (Wolfram 15.0): the generated pages link language-assets-v2 and use a
  markup that the language-assets bundle shipped with Wolfram has no css or js for. The v2 files are taken from the online
  documentation, cached next to the build folder, and copied over the shipped ones. Safe to run again.*)

$assetSite = "https://reference.wolfram.com/language-assets-v2/";

$cssFiles = {"entities", "footer", "global", "guide", "header-language", "iconography", "icons", "in-out", "language-specific", "language",
	"man-page", "related-links", "root-guide-language", "special-styles", "tables", "technotes", "tooltips", "workflows",
	"sidebar-nav", "stickynavigation", "header-v2"};
$jsFiles = {"dropdowns.js", "global.js", "language-specific.js.en", "language.js", "temporary.js", "toggles.js", "tooltips.js", "workflows.js", "stickynavigation.js"};

(*download a file once into the cache, give back its path or $Failed*)
cached[cache_String, rel_String] := Module[{f = FileNameJoin[{cache, Sequence @@ StringSplit[rel, "/"]}]},
	If[! FileExistsQ[f],
		If[! DirectoryQ[DirectoryName[f]], CreateDirectory[DirectoryName[f], CreateIntermediateDirectories -> True]];
		Quiet@URLDownload[$assetSite <> rel, f]];
	If[FileExistsQ[f], f, $Failed]
];

copyTo[out_String, src_String, rel_String] := Module[{dest = FileNameJoin[{out, "language-assets", Sequence @@ StringSplit[rel, "/"]}]},
	If[! DirectoryQ[DirectoryName[dest]], CreateDirectory[DirectoryName[dest], CreateIntermediateDirectories -> True]];
	CopyFile[src, dest, OverwriteTarget -> True];
	dest
];

fixAssets[out_String] := Module[{cache = out <> "-assets", css, imgs, f, t},
	Scan[If[(f = cached[cache, "css/" <> # <> ".css"]) =!= $Failed, copyTo[out, f, "css/" <> # <> ".css"]] &, $cssFiles];
	Scan[If[(f = cached[cache, "js/" <> #]) =!= $Failed, copyTo[out, f, "js/" <> #]] &, $jsFiles];

	(*the jquery shipped with wolfram 15.0 ends in a build log line, a syntax error that stops every script on the page*)
	f = FileNameJoin[{out, "common", "js", "jquery", "3.7.1", "jquery.min.js"}];
	If[FileExistsQ[f] && StringContainsQ[t = ByteArrayToString[ReadByteArray[f], "UTF-8"], "dependencycheck"],
		Export[f, StringTake[t, Last[StringPosition[t, "ce});"]][[2]]] <> "
", "Text", CharacterEncoding -> "UTF-8"]];

	(*root relative urls in the css cannot work from a local folder*)
	css = FileNames["*.css", FileNameJoin[{out, "language-assets", "css"}]];
	Scan[(Export[#, StringReplace[ByteArrayToString[ReadByteArray[#], "UTF-8"], {
		"url(/language-assets/img/" -> "url(../img/", "url('/language-assets/img/" -> "url('../img/",
		RegularExpression["url\\(['\"]?https://[^)]*/img/"] -> "url(../img/"}], "Text", CharacterEncoding -> "UTF-8"]) &, css];

	(*images the css and the pages ask for, not in the shipped bundle*)
	imgs = Union@Flatten[StringCases[ByteArrayToString[ReadByteArray[#], "UTF-8"], RegularExpression["\\.\\./img/([^'\")?#\\s]+)"] -> "$1"] & /@ css];
	imgs = Join[imgs, {"DetailsOptions-gray.png", "DetailsOptions-hover.png"}];
	Scan[If[! FileExistsQ[FileNameJoin[{out, "language-assets", "img", #}]] && (f = cached[cache, "img/" <> #]) =!= $Failed, copyTo[out, f, "img/" <> #]] &, imgs];
]

(*extra rules: no wolfram.com header or footer placeholders, no search (it posts to wolfram.com)*)
fixCss = "html { border-top: 0 !important; }
html, body { background: #fff !important; }
#_header, #_footer, #_footer-offset, .input-container, .search-content, #sticky-search-toggle { display: none !important; }
";

(*first the folder name, so that the root relative paths it produces are caught by the next rule*)
fixPage[t_String] := StringReplace[StringReplace[t, "language-assets-v2/" -> "language-assets/"], {
	"hoverintent.min.j2" -> "hoverintent.min.js",
	(*code-clipboard.js is an es module*)
	"<script src=\"../../common/js/clipboard/2.0/code-clipboard.js\">" -> "<script type=\"module\" src=\"../../common/js/clipboard/2.0/code-clipboard.js\">",
	(*the badge letter spaces the capitals of the package name*)
	"Q M R I TOOLS " -> "QMRITOOLS ",
	(*root relative asset urls, used by the sidebar images, only work from a web root*)
	"\"/language-assets/img/" -> "\"../../language-assets/img/",
	(*the title was utf-8 decoded twice, an em dash comes out as three or four characters*)
	"&#226;&#128;&#148;" | "&#226;&#8364;&#8221;" -> "&#8212;",
	(*context joined to the url without a slash: ../QMRIToolsSegmentationTools/ref/x.html*)
	RegularExpression["href=\"\\.\\./QMRITools\\w+/(ref|guide)/"] :> "href=\"../$1/"}];

(*the sidebar, sticky header and fix stylesheets and the sticky header script are not in the generated pages*)
addAssets[t_String] := StringReplace[
	StringReplace[t, If[StringContainsQ[t, "qmritools-fix.css"], {}, {"</head>" :> "<link rel=\"stylesheet\" href=\"../../language-assets/css/sidebar-nav.css\"><link rel=\"stylesheet\" href=\"../../language-assets/css/stickynavigation.css\"><link rel=\"stylesheet\" href=\"../../language-assets/css/header-v2.css\"><link rel=\"stylesheet\" href=\"../../language-assets/css/qmritools-fix.css\">
</head>"}], 1],
	If[StringContainsQ[t, "stickynavigation.js"], {}, {"<script src=\"../../language-assets/js/dropdowns.js\">" :> "<script src=\"../../language-assets/js/stickynavigation.js\"></script><script src=\"../../language-assets/js/dropdowns.js\">"}], 1]

(*every reference page has an Examples section holding only the rendered input Needs["QMRITools`"], one png per page.
  That section and its sidebar entry are removed, a section with any other example is left alone.*)
fixExamples[t_String] := If[
	StringCount[First[StringCases[t, RegularExpression["(?s)<section id=\"Examples\">.*?</section>"], 1], ""], "<figure"] === 1 &&
		StringContainsQ[First[StringCases[t, RegularExpression["(?s)<section id=\"Examples\">.*?</section>"], 1], ""], "Needs[&quot;QMRITools`&quot;]"],
	StringReplace[t, {
		RegularExpression["(?s)<section id=\"Examples\">.*?</section>"] -> "",
		RegularExpression["(?s)<div class=\"examples-container\">\\s*<div class=\"primaryexamplesection menu-item\">.*?</div>\\s*</div>\\s*</div>"] -> ""}],
	t];

FixHTMLDocs[out_String] := Module[{pages, n = 0, t, new},
	fixAssets[out];
	Export[FileNameJoin[{out, "language-assets", "css", "qmritools-fix.css"}], fixCss, "Text"];
	pages = FileNames["*.html", FileNameJoin[{out, "html"}], Infinity];
	Scan[(
		t = ByteArrayToString[ReadByteArray[#], "UTF-8"];
		new = addAssets[fixExamples[fixPage[t]]];
		If[new =!= t, n++; Export[#, new, "Text", CharacterEncoding -> "UTF-8"]]
	) &, pages];
	Print["html pages fixed: ", n, " of ", Length[pages]];
	(*the Files folders hold the example images, delete them when no page links one*)
	If[! AnyTrue[pages, StringContainsQ[ByteArrayToString[ReadByteArray[#], "UTF-8"], "data-src=\"Files/"] &],
		Scan[If[DirectoryQ[#], DeleteDirectory[#, DeleteContents -> True]] &, FileNameJoin[{out, "html", #, "Files"}] & /@ {"ref", "guide"}]];
]
