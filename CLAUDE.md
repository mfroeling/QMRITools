# QMRITools: agent notes

Wolfram Language paclet for quantitative muscle MRI. Source lives in `QMRITools/Kernel/*.wl`. Load the dev version
with `<< QMRIToolsDev`` (it points `PacletDirectoryLoad` at `QMRITools/QMRITools`).

Read before editing:

- [agents/CodeStyle.md](agents/CodeStyle.md): house style and working agreements (comments, `Block`/`With`, no
  duplication, verifying with wolframscript). Always applies.
- [agents/MuscleBidsTools.md](agents/MuscleBidsTools.md): map of the config-driven Muscle-BIDS pipeline
  (`MuscleBidsTools.wl`): stages, naming, config keys, check files, known quirks. Read it before touching that file
  or a study `config.json`.
- [agents/SegmentationTools.md](agents/SegmentationTools.md): map of `SegmentationTools.wl`: inference
  (`SegmentData`), anatomy tables, patching, training loop and producer kernels, augmentation and pretraining, data
  prep, metrics, known quirks.
