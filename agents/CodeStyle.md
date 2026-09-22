# QMRITools code style and working agreements (for agents)

These are the house conventions for the Wolfram Language code in `QMRITools/Kernel/*.wl` (~30 files, ~45k lines),
and for how to work with the author (Martijn Froeling). They come from repeated corrections across many sessions and
from a survey of the whole Kernel folder. **Match the surrounding code, and follow these rules even when your own
default differs.**

Related: [MuscleBidsTools.md](MuscleBidsTools.md) (map of the BIDS pipeline).

---

## 1. The five rules that get corrected most often

1. **Comments: at most one short line, stating *what* the code does, never *why*.** No comment blocks above a
   function, no multi-line comments, and no rationale in an edit unless the user asks for it. The default number of
   new comments is zero.
2. **`Block[{a, b, c}, ...]` with bare symbols only.** No initialisers in the header. Write
   `Block[{t0}, t0 = AbsoluteTime[]; ...]`, not `Block[{t0 = AbsoluteTime[]}, ...]`.
3. **`Block` is the default scoping construct. Use `With` only when a value must survive scope exit** (a returned
   closure, a value embedded in a data structure, code sent to another kernel). Avoid `Module` in new code.
4. **No duplicated logic.** Before adding a case that resembles existing ones, restructure the dispatch: a `Switch`
   becomes gated `If`s, a pattern alternation gets broader, a helper gets extracted. Never copy a body and tweak it.
5. **Compact code.** Pull setup or teardown shared by both `If` branches out so it runs once. Test a condition once
   instead of three times nearby. Pure functions and `/@` beat `Do`/`For`/`Table[With[...]]`.

---

## 2. File anatomy (every Kernel file)

```wolfram
(* ::Package:: *)

(* ::Title:: *)
(*QMRITools XxxTools*)

(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)

(* ::Section:: *)
(*Begin Package*)

BeginPackage["QMRITools`XxxTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`XxxTools`"}]]];

(* ::Section:: *)
(*Usage Notes*)

(* ::Subsection::Closed:: *)
(*Functions*)
PublicFunction::usage = "PublicFunction[x] does ...
PublicFunction[x, y] same but ...";

(* ::Subsection::Closed:: *)
(*Options*)
SomeOption::usage = "SomeOption is an option for PublicFunction. ...";

(* ::Subsection::Closed:: *)
(*Error Messages*)
PublicFunction::tag = "Message text with `1`.";

(* ::Section:: *)
(*Functions*)

Begin["`Private`"]

(* ::Subsection:: *)
(*Group name*)

(* ::Subsubsection::Closed:: *)
(*PublicFunction*)

Options[PublicFunction] = {...};

SyntaxInformation[PublicFunction] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

PublicFunction[...] := ...

(* ::Section:: *)
(*End Package*)

End[]

EndPackage[]
```

- **All `::usage` strings go in the front "Usage Notes" section**, never next to the implementation. Options get their
  own subsection.
- **All `f::tag` message declarations go in the front "Error Messages" subsection**, one block per function.
- Each function sits in its own `(* ::Subsubsection::Closed:: *)` cell. These cell markers are the navigation
  structure and must be preserved exactly (the files are also edited as notebooks).
- A private helper has no `::usage` and lives after `Begin["`Private`"]`.
- Every subpackage sees every other subpackage's public symbols through `QMRITools`$Contexts`.
- Leave two blank lines between cells and one blank line between overloads of the same function.

## 3. Formatting

- **Tabs** for indentation.
- A long `Block` locals list wraps like this:

  ```wolfram
  F[x_] := Block[{
  		a, b, c, d,
  		e, f
  	},
  	body
  ]
  ```

- Closing brackets of long `Table`/`Switch`/`If` bodies get a trailing comment:
  `, {set, sets}](*close loop over sets*)`, `];(*close method switch*)`.
- Visual separators inside long functions: `(*-------------------------------------------*)` banners around
  each `Switch` case, and `(* -------------- Logging -------------- *)` style dividers.
- Ordinary tab indentation. No alignment gymnastics.

## 4. Definitions and dispatch

- Public function scaffold: `Options[f]` → `SyntaxInformation[f]` → the definitions.
- **Optional or polymorphic arguments use overloads, not internal `If[Head[x]===...]`:**
  - *Arity overloading*: the short form fills defaults and tail-calls the full form:
    `Mask[dat_?ArrayQ, opts:OptionsPattern[]] := Mask[dat, {0,0}, opts]`.
  - *Shape overloading*: dispatch on pattern tests (`_?StringQ`, `_?AssociationQ`, `_?ListQ`, `list:{_?AssociationQ ..}`).
    Lists usually map onto the scalar form: `GenerateBidsName[list:{_?AssociationQ..}] := GenerateBidsName /@ list`.
  - Pipeline entry points use a 3-tier chain: `f[dir]` → `f[dir, config]` → `f[in, out, desc]`
    (see MuscleBidsTools.md §3).
- **`Switch` on string option values** for mutually exclusive dispatch. Use `Which` for condition chains. When a new
  variant shares behaviour with several existing cases, **convert to independent `MemberQ`-gated `If`s** so the
  bodies are reused. The precedent is `Dixon-A` in `MuscleBidsProcessI`:

  ```wolfram
  If[MemberQ[{"Dixon-S", "Dixon-P", "Dixon-A"}, met], (*branch 1*) ...];
  If[MemberQ[{"Dixon", "Dixon-B", "Dixon-A"}, met], (*branch 2*) ...];
  ```

  Each branch has its own "files not found → skip" guard, so no branch needs to know about the others.
- Options: read them all at once, `{a, b, c} = OptionValue[{OptA, OptB, OptC}];`.
- Pass options through with `opts` / `FilterRules[{opts}, Options[g]]` where needed.

## 5. Scoping in detail

- `Block` everywhere, with dynamic scope. `Module` is legacy (mostly in `Legacy.wl` and a few older files), so do not
  introduce it.
- `With` is correct and required when:
  - a closure or pure function is **returned** or **stored** (for example in an options list), which would otherwise
    capture a `Block` local that vanishes;
  - code is sent to a remote kernel with `LinkWrite[link, Unevaluated[...]]`. Master-side values must be spliced in:
    `With[{l = load}, LinkWrite[link, Unevaluated[Get[l]]]]`. Plain function parameters (`f[link_, x_] := ...`,
    `#1`) are already substituted and **do not** need `With`;
  - `ParallelSubmit`/`DistributeDefinitions` closures bake in runtime-local values.
- When a nested `With` exists only to share one computed value, promote that value to a **function parameter**
  instead (`OneCycleSchedule[br_, rounds_, itt_]`).
- Prefer `f[#] & /@ Range[n]` over `Table[With[{qi = i}, f[qi]], {i, n}]` when `f` is not HoldAll.
- Private symbols are not on `$ContextPath` outside the package. In test scripts use the full context,
  `` QMRITools`SegmentationTools`Private`SetMXenvironment `` or `ToExpression["QMRITools`MuscleBidsTools`Private`" <> "SubNameToBids"]`.

## 6. Naming

- Public functions are PascalCase and verb-first: `MaskData`, `TensorCalc`, `MuscleBidsMerge`.
- Internal worker of a public function: suffix `I` (`MuscleBidsProcessI`, `BidsDcmToNiiI`).
- Compiled variants: suffix `C` (`DotC`, `SaltAndRiceC`). Vectorised or interpreted siblings: `V` or `I`.
- Locals are short lowercase abbreviations, reused consistently: `dat`/`data`, `vox`, `dim`, `mask`, `seg`, `grad`,
  `val` (bvals), `mon`, `fol`, `out`, `met`, `opts`, `pos`, `n`.
- Package flags get a `$` prefix: `$debugBids`. Message tags are short fragments: `Mask::tresh`,
  `CheckDataDescription::man`.
- MuscleBids outputs: **the local variable name equals the output file suffix**. Files are exported with
  `ToExpression[Context[con] <> name]`, so renaming a variable renames a file.

## 7. Comments: examples

Good (existing style):

```wolfram
(*get random dataset*)
(*check if files are already done*)
(*returns {links, nProducers} - links is {} for the "Parallel" path*)
```

Bad (these were rejected):

```wolfram
(* This helper launches the training kernels. It returns a pair where the first
   element is the list of links and the second is the number of producers. *)
F[x_] := ...

(*use 0.2 because 0.1 made the validation set too small for the smaller cohorts*)
```

- Keep comments **accurate**. After changing a value or reordering lines, re-check nearby comments. A stale
  "10%" after a change to `0.2` was caught, and so was a comment left above the wrong line.
- **Keep commented-out debug scaffolding** (`(*If[...,Print[...]];*)`, `(*echo = Lookup[...];*)`). The user reuses it.
  Do not clean it up unless asked.
- `(*TODO ...*)` blocks the user wrote stay as they are.

## 8. Logging and debugging idioms

- BIDS code: `(*-----*)AddToLog[{"msg", value}, level, True]`. The `(*-----*)` prefix marks log lines visually.
  `level` is indentation (0–5). `True` adds a timestamp (it can go before or after the level).
- `debugBids[...]` prints only when `$debugBids` is True. Sprinkle these generously in pipeline code (paths,
  dimensions). They are cheap and the user relies on them.
- Other files: a `MonitorCalc`/`Monitor` option controls `PrintTemporary`/progress output. Pass `MonitorCalc->False`
  when calling from a pipeline.

## 9. Error handling

- Targeted, not exhaustive. Guard the 2–3 failure modes that are known to happen near the top of a function:
  `Return[Message[f::tag, x]; $Failed]`. Nested guards use `Return[..., Block]`, the deliberate two-argument form.
- Much validation is implicit: a pattern test fails to match and the call stays unevaluated. **Do not add broad
  defensive checks.**
- Pipeline code logs and skips (`AddToLog[...]; Return[]`) instead of throwing, so one bad subject does not stop a
  batch.
- `Quiet@` is used around known-noisy calls (`Quiet[CreateDirectory[...]]`, `Quiet@GetConfig`). Do not spread it
  further.

## 10. Performance idioms

- Numeric inner loops use typed `Compile[..., RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]`
  (sometimes `{"Speed", "WarningMessages" -> False}`). Copy the options of the neighbouring compiled functions.
- Keep arrays packed: `ToPackedArray@N@...`, and avoid mixing Integer and Real.
- The standard parallel idiom is `fun = If[par, DistributeDefinitions[...]; ParallelMap, Map]`. The link-based
  producer/consumer pool in `TrainSegmentationNetwork` is a bespoke exception, not a template.

## 11. Intentional designs, not bugs (do not "fix")

- `AugmentTrainingData` (SegmentationTools.wl): salt/pepper hard-coded to `1.`/`0.` (data is pre-scaled so
  Q99 = 1); only the last axis is flipped (left-right mirroring); scale range `{0.6, 1.6}` per axis. (Data now uses
  interpolation order 1 and the seg order 0.)
- `FiberTractography`: the serial branch builds one multi-valued interpolator and the parallel branch builds N depth-3
  ones (cheaper to broadcast). The commented-out `vecInt = MakeInt[...]` in the parallel branch documents this on
  purpose.
- Pretraining (`AugmentMask`): `input = data(1-mask) + fill·mask`, `target = data·mask` (not dilated, which was
  tried and reverted); the `MakeBlockMask` `2*RandomReal[{0.1,0.5}]` factor; `FreezeEncoderDepth` is a contiguous
  cutoff.
- MuscleBids: the Dixon gated-`If` branches (§4). Also see the quirk list in MuscleBidsTools.md §10: ask before
  changing any of it.

## 12. Working with the user

- **Treat proposals as drafts.** The user iterates, often redirecting a fix toward a cleaner shape (a parameter
  instead of `With`, gated `If`s instead of a new case). Offer the diff and expect a redirect. Do not defend the first
  version.
- **Verify semantics empirically before asserting them**, especially holding, evaluation order, `With`/`Block`
  capture, and what a pattern matches. The user pushes back with "are you sure". Answer with a small test:

  ```powershell
  & "C:/Program Files/Wolfram Research/WolframScript/wolframscript.exe" -file "<scratchpad>\test.wls"
  ```

  To load the dev package in a script (this is what `<< QMRIToolsDev`` does):

  ```wolfram
  PacletDirectoryLoad["D:\\werk\\workspace\\QMRITools\\QMRITools"];
  Quiet@Get["QMRITools`"];
  ```

  The first load takes about a minute. Keep test scripts in the session scratchpad, not in the repo.
- The user debugs from live notebook output they paste in (including parallel-kernel prints). Trace the actual
  evaluation to find the mismatch instead of guessing from the error text.
- When fixing ordering or capture bugs, **trace every affected call site**. Do not patch only the reported symptom.
- Scope discipline: do what was asked. Mention nearby issues you notice, but do not fix them unasked (especially in
  MuscleBidsTools, where file names on disk depend on current behaviour).
- Before `LinkLaunch`-based development, the front-end setting "Launch parallel kernels: At startup" must be off.
  Otherwise every launched sub-kernel spawns its own pool, and hundreds of processes appear.
- Commit or push only when asked. The main branch is `master`.

## 13. Pre-flight checklist for any edit

- [ ] New public symbol? Usage string in the front section, and `Options`/`SyntaxInformation` present.
- [ ] New message? Declared in the front Error Messages block.
- [ ] `Block` header has bare symbols only, and every new local is declared (no leaked globals).
- [ ] `With` only where a value must escape the scope.
- [ ] No copied branch bodies. Dispatch restructured instead.
- [ ] Comments: zero or one short line each, no rationale, nearby comments still accurate.
- [ ] Debug scaffolding and `(*close ...*)` markers preserved.
- [ ] Semantics you relied on were verified with wolframscript when not obvious.
- [ ] MuscleBids: output names unchanged, or every downstream stage updated (see MuscleBidsTools.md §4, §11).
