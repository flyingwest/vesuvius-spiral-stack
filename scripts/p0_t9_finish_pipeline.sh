#!/usr/bin/env bash
set -euo pipefail
cd /home/wpotech/deployments/vesuvius-spiral-stack
# Await only the already active acquisition; never restart an active downloader.
while kill -0 3635784 2>/dev/null; do sleep 20; done
test -f outputs/p0_t9/acquisition.json
image='nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a'
docker run --rm --gpus all -v "$PWD:/work" -w /work -e PYTHONPATH=/work/data/cache/p0_t9/python_site:/work/external/villa/vesuvius/src "$image" python -B scripts/p0_t9_smoke.py > outputs/p0_t9/real_ct_smoke.log 2>&1
.venv/bin/python -c "import json; assert json.load(open('outputs/p0_t9/smoke_test.json'))['status']=='PASS'"
docker run --rm --gpus all -v "$PWD:/work" -w /work -e PYTHONPATH=/work/data/cache/p0_t9/python_site:/work/external/villa/vesuvius/src "$image" python -B scripts/p0_t9_infer.py > outputs/p0_t9/full_inference.log 2>&1
MPLCONFIGDIR=/home/wpotech/.hermes/cache/scratch/mpl-t9 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B scripts/p0_t9_evaluate.py > outputs/p0_t9/evaluation.log 2>&1
