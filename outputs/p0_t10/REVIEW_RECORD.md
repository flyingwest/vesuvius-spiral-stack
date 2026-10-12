# Pre-evidence review

Reviewer: separate online agent `/root/t10_preflight_review`, same inherited
model/provider as author (Codex Astra). Independent context, not independent
provider/model. Read-only; no local LLM, inference, or file changes by reviewer.

VERDICT: APPROVE-FOR-EVIDENCE, not approval of a final A/B/C/D claim.
BLOCKER: none.

MAJOR: Final geometry conclusion needs measured evidence. The plan says
“A if control passes and relevant geometry/placement differences are
demonstrated”. Paired images alone do not quantify trusted frames, sign,
depth origin or papyrus placement. Distinct methods alone are insufficient;
scan resolution/energy and output sampling remain confounders.
Disposition: open until measured audit; no automatic A from control PASS.

MINOR: Official API comparison reuses `configured_model=configured` and
official helpers. It verifies scheduling/blending conformance, not an
independent model implementation. Label discrimination supplies sensitivity.
Disposition: disclosed; no independent-model claim.

MINOR: Single-class bootstrap draws and zero within-class SD could abort or
produce undefined metrics. Original code asserted all bootstrap draws finite.
Disposition: corrected before model evidence. Any degenerate bootstrap yields
INCONCLUSIVE; zero pooled SD uses effect0. No threshold was changed.

Verified by reviewer: source/label URL and level match; source depth[13,97),
four-plane means, annotation center10 retained as model depth8; forward
layers[2,19); source axes and chunk conventions; official/T9 normalization
equivalence; frozen stride, blending and training-corpus scope.

Pre-evidence source/plan commit: `9b8c7e2`.

Follow-up: reviewer confirmed the degenerate-metric fixes and unchanged gates
at `9b8c7e25ce2542715b1900dec047c19f3030d818`, still APPROVE-FOR-EVIDENCE.
Residual MINOR fixed before inference: conditional patch audit now checks
`result['status'] == 'PASS'`, not the numerical-gates-only `passed` variable,
so an INCONCLUSIVE bootstrap cannot enter the pass-only branch.

Post-evidence review: APPROVE scientific conclusion A with causal limits.
MAJOR CLOSED: reviewer independently recomputed median |normal Z|=0.349152
over10,712 trusted-grid points, checked all finite-difference neighbors and
nondegenerate cross products, and confirmed median V tilt25.03degrees versus
T9's enforced scan-Z direction and nz=0. Identically displayed patches show
repeated cross-cut bands/empty regions in T9. Control passes all frozen gates;
all49patches and pooling replay exactly, probability maximum error0.

Limits retained: different locations/scans; not true normals at T9 centers;
producer sign/interpolation flags unavailable; intensity maxima not an oracle
for ink or seating; training-corpus control; API comparison shares helpers.
Geometry is implicated, not established as the sole cause. No new scientific
blocker. Preservation/checksum/Git checks are separate from review approval.
