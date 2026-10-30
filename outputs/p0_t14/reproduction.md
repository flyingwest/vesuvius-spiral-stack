# Execution record / reproduction

Use repository Python `.venv/bin/python -B`. The analysis script modes are:
`preserve`, `geometry <fresh-subdirectory>`, and `stack <fresh-subdirectory>`.
Primary destinations were `geometry`, `stack`; replays `geometry_replay`,
`stack_replay`. Stack consumes `geometry/threads.npz` and freshly inferred
`inference/official_control_probability.npy`, never the old T10 prediction.

Inference was executed twice, second with `--replay` appended:

```bash
docker run --rm --gpus all --ipc=host \
  -v "$PWD:/work:ro" -v "$PWD/outputs/p0_t14:/work/outputs/p0_t14" \
  -w /work \
  -e PYTHONPATH=/work/data/cache/p0_t9/python_site:/work/external/villa/vesuvius/src \
  -e TMPDIR=/work/outputs/p0_t14/runtime_tmp --entrypoint python \
  nvcr.io/nvidia/pytorch@sha256:3becd068f49bd2ad38f90db5f9a4803019a76933a24e63d821376c44e7a9200a \
  -B scripts/p0_t14_infer.py
```

Only T14 output was writable in this container. No downloads or model updates.
Primary stdout/stderr and replay logs are retained. Old patches/source arrays were
read-only. Analysis used NumPy/SciPy/tifffile/Pillow/matplotlib in the existing venv.

Verification was `scripts/p0_t14_verify.py`; sealing was `scripts/p0_t14_seal.py`.
Do not rerun write stages in the sealed directory. Clone/use isolated destinations
for future replay. Read-only verification is `sha256sum -c outputs/p0_t14/SHA256SUMS`.
