# P0-T13 — local cross-slice ribbon: execution blocked

P0_T13_COMPLETE = NO. No scientific A/B/C/D result is assigned.
Correspondence rate is unmeasured; local ribbon and geometry gate are not evaluated.
Smoke and full inference are NOT_RUN; ink structure is NOT_ASSESSED.
READY_FOR_NEXT_STEP = NO pending the required control inputs or a scope decision.

## Hypothesis and frozen scope

Neighboring independent frozen 2-D threads within Z−2..Z+2 may provide enough
local geometric evidence to construct a bounded 3-D ribbon, without changing
any thread coordinates, bridge decisions, masks, provenance, arc length or row
ordering. The intended inputs remain T5 Z128..177 and the exact 72 T11/T12
locations at Z128..130,151..153,175..177. No experiment has tested this T13
hypothesis yet. T12 remains complete, result D, with failed local-sheet and
patch-geometry gates.

## Required input blocker

The existing sealed `outputs/p0_t13/CHECKPOINT.md` and `preflight_blocker.json`,
committed at `37e0baf30ca52ebe00e786ed869c9c6a5a3be43b`, already document the
missing prerequisite. The checkpoint states:

> Resume only after receiving existing independent control threads with registration, or authorization for a separate bounded control-only thread reconstruction.

T5/T9/T11 target coordinates use `54keV_7.91um_Scroll1A.zarr`, native 7.91 µm.
The trusted T10 segment `20231031143852` uses the 2.4 µm/78 keV scan
`2.4um-0.22m-78keV-volume-20260411134726.zarr`. The sealed T10 report explicitly
identifies different locations and no target-wise reference normals. No registered
overlap or independently traced control threads have been supplied for this resume.
Scan names alone do not establish different physical specimens; the missing
registration is the relevant limitation.

Consequently, surface-normal error, longitudinal-direction error, center-placement
error, ribbon-to-trusted-surface distance and correspondence correctness cannot
be measured. Similar median tilt distributions (T11 approximately 24.73°, T10
approximately 25°) are not paired frame-error measurements. T10's producer patches
also cannot serve as identical-source, identical-location D patches for these targets.
Extracting input threads from trusted geometry would violate its evaluation-only
role. Constructing new independent control threads would exceed this frozen scope.

Gate condition 2, material improvement in trusted-control frame error over T11,
therefore cannot be evaluated. Missing inputs are not a scientific FAIL or result D.
No gate was weakened or replaced, no correspondence parameters were preregistered,
and no model was loaded. No inference or geometric success is claimed.

## Work performed and preservation

README.md, AGENTS.md and the T1–T12 reports were reviewed, including T11's support
amendment and T12's authoritative failed result. Prior implementations and input
schemas were inspected read-only. An unexecuted implementation/plan draft was
prepared; its setup refused the existing T13 evidence directory before writing
anything. Inspection then identified the sealed preflight checkpoint. Both drafts
were removed without executing correspondence, ribbon construction or rendering.
They are not a registered scientific method or negative experimental result.

A fresh read-only SHA-256 audit verified **3,999 predecessor/checkpoint entries**,
with **zero mismatches**. This includes the T13 checkpoint and its original seal;
P0-T1 through P0-T12 remain unchanged. Audit receipt:
`outputs/p0_t13/resume_audit/verification.json`.
The original `outputs/p0_t13/SHA256SUMS` remains unchanged. A separate
`outputs/p0_t13/resume_audit/SHA256SUMS` seals this report and the new receipt.
This is preservation verification, not correspondence/ribbon or model replay.

No correspondence metrics, candidate scores, ribbon frames, frame-error plots,
rendered-patch comparison, patch-quality metrics or smoke metrics were generated.
Their absence is intentional: fabricated zero values or placeholder scientific
results would misrepresent an experiment that did not run.

## Required next input and stopping point

Supply the existing independent frozen control threads and their registration to
T10's trusted surface, or separately authorize and specify a bounded control-only
thread reconstruction. The latter is a research-scope decision for the user and
ChatGPT, not an implementation choice made here. No such reconstruction was started.
Once the prerequisite is resolved, preregister the bounded correspondence/ribbon
method and numerical gates before collecting T13 evidence.

This blocked-execution report and preservation receipt are published as a checkpoint,
not a sealed completed P0-T13 result. Commit and remote identity are reported in the
execution response to avoid embedding a commit's own hash in its content. No force
push is used. Stop here; no next experiment executed.
