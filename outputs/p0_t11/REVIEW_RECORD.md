# P0-T11 pre-evidence review

Reviewer: online independent-context agent `/root/t11_preflight_review`,
same inherited Codex Astra model/provider, not an independent provider.
No local LLM, target CT evidence, or model inference used by reviewer.

Final VERDICT: APPROVE-FOR-CT-FRAME-EVIDENCE.

MAJOR full replay missing: FIXED before evidence. Distinct full replay branch
verifies50row coordinate/index/stack mappings and recomputes one deterministic
frame/prediction per row; tolerance .002. Static review, not execution proof.

MAJOR incomplete missing-data handling: FIXED before evidence. Either old or
corrected missing CT is recorded as DATA_BLOCKER and terminates before
scientific classification, CSV/plots or inference. All72paired patches required.

MINOR rejected smoke frames differ from full policy: FIXED. Both retain the
original T9 frame for rejected estimates; all records preserved.

MINOR model-plane labels incorrectly native XY: FIXED, now model X-Y.

MINOR sparse-station angles are not local smoothness: DISCLOSED in plan/report.
Perturbation test is local +/-2native Z; no longitudinal continuity claim.

Synthetic regression PASS: six tilted sheets, blank-volume rejection,
handedness/depth sign, affine-coordinate axis sampling. Initial development
einsum covariance expression failed this test and was replaced by the explicit
weighted matrix product before any target evidence. No gates were relaxed.

Scientific limitations: constrained local planar patch, not curved surface
following; CT quality proxies are not true sheet identity or ink evidence.
Conditional inference requires frozen geometry/smoke gates. Final evaluator
fidelity, replay, preservation, and Git push require separate verification.
