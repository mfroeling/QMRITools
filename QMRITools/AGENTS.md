# QMRITools: agent entry point

QMRITools is a Wolfram Language paclet for quantitative (muscle) MRI. It is written and maintained by Martijn
Froeling. The source is `QMRITools/Kernel/*.wl` (~30 subpackages, one context each, `QMRITools`<Name>``).

This file is the index for all agent documentation. The guides live in the `agents/` folder at the repository root
(one level above this file). Read the guide for the file you are about to change **before** editing it.

## Documentation

| Guide | Read when | Covers |
| --- | --- | --- |
| [CodeStyle.md](../agents/CodeStyle.md) | **always**, before any edit | house style (comments, `Block`/`With`, overloads, no duplication), file anatomy, logging, error handling, intentional designs that look like bugs, how to work with the user, pre-edit checklist |
| [DataConventions.md](../agents/DataConventions.md) | passing data between functions, writing a fit, NIfTI or gradient I/O, anything with units | array layouts (`{z, t, y, x}`, tensor `{6, z, y, x}`), `vox = {z, y, x}`, NIfTI import/export defaults, diffusion/IVIM/relaxometry/Dixon units (including the Dixon seconds-vs-ms discrepancy), masks and labels, normalisation, sanity-check ranges |
| [ToolboxStructure.md](../agents/ToolboxStructure.md) | changing the loader (`Kernel/QMRITools.wl`), `PacletInfo.wl`, the package list, assets or docs; adding a package or public function; "works only on second load" problems | paclet layout, assets per OS, the loader step by step, the depth-first cascade and required second `Get` pass (shadow symbols), dev reload caveats, package map with measured dependencies and cycles, documentation and build, doc gaps, checklists |
| [MuscleBidsTools.md](../agents/MuscleBidsTools.md) | changing `MuscleBidsTools.wl`, or reading or editing a study `config.json` | the 7-step Muscle-BIDS pipeline: call chain, folder layout, file naming, every config key the code reads, check files, differences from `BIDS-config.docx`, quirk list, recipes |
| [SegmentationTools.md](../agents/SegmentationTools.md) | changing `SegmentationTools.wl`, or working on segmentation inference or network training | anatomy tables, the `SegmentData` pipeline, patching, `TrainSegmentationNetwork` and producer kernels, augmentation and masked pretraining, data prep, metrics, quirk list, recipes |

### Cross-references between the guides

Each guide's header links to the related guides. The shared topics are:

- **Structure → every guide.** ToolboxStructure.md explains how a single subpackage (CodeStyle.md §2) fits into the
  loader and the dependency graph. Where a helper lives (for example `JoinSets` in ProcessingTools) is listed there.
- **Style → every guide.** Both code maps assume CodeStyle.md. Its §11 ("intentional designs") lists behaviour in both
  files that must not be "fixed".
- **BIDS ↔ Segmentation.** `MuscleBidsSegment` calls `SegmentData` (MuscleBidsTools.md §6.5 ↔ SegmentationTools.md
  §2–3). A BIDS config's `Segment.Location` is a `what` group from `$SegmentationGroups`. `MuscleBidsAnalysis` uses the
  label files and `MuscleLabelToName` (SegmentationTools.md §9).
- **Quirk lists.** Each code map ends with a numbered quirk table (MuscleBidsTools.md §10: Q1–Q14,
  SegmentationTools.md §10: Q1–Q13). Each entry is marked fixed, intentional, or open. **Ask the user before fixing an
  open one.**

## Essentials

- **Load the dev version**: `<< QMRIToolsDev`` in a notebook. In a script:
  `PacletDirectoryLoad["D:\\werk\\workspace\\QMRITools\\QMRITools"]; Get["QMRITools`"]`.
- **Verify semantics empirically** with wolframscript (`C:/Program Files/Wolfram Research/WolframScript/wolframscript.exe
  -file test.wls`), keeping test scripts outside the repo. Private symbols need their full context,
  e.g. `` QMRITools`SegmentationTools`Private`NormDat ``.
- **Public symbols** need a `::usage` in the front "Usage Notes" section, and `Options`/`SyntaxInformation`.
  Otherwise the function stays private.
- **Research notebooks** that call the package are outside the repo, under `D:\Werk\Research\` (for example
  `Segmentation`, `neural segment`, `MOTOR_shoulder_2`, and study folders with a `config.json`). Check them before
  calling a function unused.
- Commit or push only when asked. The main branch is `master`.

## Keeping the guides current

When you change behaviour, naming, config keys, or a quirk's status in a documented file, update the matching guide
in the same change (section, quirk table, function index). Line numbers in the guides are approximate. Grep for the
function name instead of trusting them. New guides go in `agents/` and get a row in the table above.
