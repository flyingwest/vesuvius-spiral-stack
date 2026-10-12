# P0-T10 bounded positive-control plan

Fixed before CT discrimination/inference; label-only scouting is disclosed.

Hypothesis: the unchanged T9 seed42 step-075000 checkpoint is ink-sensitive on
an official producer-rendered, registered surface. This is a training-corpus
positive control, not held-out performance and not a replacement for the
Direct Thread-Stack method. Preserve all T1–T9 artifacts and source.

Memory retry: execute the sealed T9 smoke function with only its output path
redirected to T10. No full T9 inference. A pass establishes current execution
compatibility; it does not prove Nemotron was the sole historical OOM cause.

Control: phercparis4-w02, public segment 20231031143852, official 2.4 µm
78 keV surface volume, XY level2. Smallest rendered area among the eight
Paris4 cases in the released ink_9um source table, not a global minimum.
Label-only first qualifying center-out 512×512 crop: Y[4096,4608),
X[5248,5760). Initial crop has 145587 supervised pixels, including 43032 ink.
The selection rule and alternative fixed candidate order are saved. No
control CT/model outputs were inspected to choose this crop.

Acquire 16 original CT chunks plus producer geometry coordinates; HEAD all
objects first and require total <100,000,000 bytes. Preserve exact raw payloads.
Use official prepare_isotropic_input: source depth[13,97), rounded means of
four planes, 21 pooled planes. Use official forward center crop[2,19), depth17,
128×128 XY, 64pixel inference stride, official Hann overlap blending, batch8,
float16 autocast, official checkpoint-derived robust preprocessing. No TTA,
normal reversal, depth search, fitting or change after evidence. Compare raw
floating output with the complete official flat-inference API (encoded error
at most one uint8 level). The official rendered surface is consumed, not
reconstructed or modified.

Evaluate registered producer labels only inside supervision, excluding64pixel
crop borders and pixels within2pixels of the opposite label class. Require at
least1024pixels per class. Report ROC AUC, PR AUC as average precision, class
probability quantiles/means, standardized mean difference, and IoU/Dice at
fixed p>0.5. Bootstrap1000resamples of64×64spatial blocks, seed2910; no iid-pixel
confidence claim. Positive-control PASS requires all: ROC AUC≥0.75, AP≥class
prevalence+0.15, effect d>0.5, and spatial-bootstrap95%AUC lower bound>0.5.
These are engineering positive-control gates, not clinical/generalization
claims. A degenerate single-class bootstrap draw makes the result INCONCLUSIVE;
constant probabilities use d=0 rather than undefined division. No gates change.
Insufficient data or registration/integrity failure is INCONCLUSIVE,
not proof of invalid model integration.

If control passes, compare fixed control patches at local origins(128,128),
(192,192),(256,256) with the already saved T9 smoke patches Z128/152/177.
Use identical uint8 display range0..255, central depth and orthogonal cuts.
Do not identify CT darkness as ink. Audit measured surface frame and placement
differences; distinguish implicating geometry from proving it is the sole cause.
Other confounds include scan energy/native resolution and sampling density.

Decision: A if control passes and relevant geometry/placement differences
are demonstrated; D if it passes but investigated differences do not explain
failure; B if a suitable registered control fails its frozen model gates;
C if no suitable reliable control can be established. Do not alter T9 FAIL.

Replay all49control patches from saved source, requiring exact raw-patch hashes,
independent exact pooling and probability max error≤0.002. Verify predecessor
preservation and checksum raw/lossless evidence. Commit/push final T10 only,
without force. No local LLM; no training, OCR, language priors or scale-up.
