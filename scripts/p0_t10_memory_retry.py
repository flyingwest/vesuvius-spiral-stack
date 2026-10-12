"""Run the exact sealed T9 smoke with only its output destination redirected."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import torch
import p0_t9_smoke as smoke

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/p0_t10'

def main():
    retry = OUT / 'memory_retry'
    retry.mkdir(exist_ok=False)
    shutil.copyfile(ROOT / 'outputs/p0_t9/acquisition.json', retry / 'acquisition.json')
    free, total = torch.cuda.mem_get_info()
    record = {'status': 'RUNNING', 'checkpoint': 'scrollprize/ink_9um/hybrid_3d2d-seed42/step-075000.pth', 'cuda_free_before_load_bytes': free, 'cuda_total_bytes': total, 'device': torch.cuda.get_device_name(), 'cuda_version': torch.version.cuda, 'system_meminfo_before': Path('/proc/meminfo').read_text(), 'gpu_before': subprocess.check_output(['nvidia-smi'], text=True), 'sealed_smoke_sha256': hashlib.sha256((ROOT / 'scripts/p0_t9_smoke.py').read_bytes()).hexdigest(), 'output_redirect_only': True, 'full_T9_rerun': False, 'local_llm_dependency': 'NONE'}
    smoke.OUT = retry
    try:
        smoke.main()
        record['status'] = 'PASS'
        record['smoke'] = json.loads((retry / 'smoke_test.json').read_text())
    except Exception as exc:
        record['status'] = 'FAIL'
        record['error'] = repr(exc)
        raise
    finally:
        record['gpu_after'] = subprocess.check_output(['nvidia-smi'], text=True)
        record['system_meminfo_after'] = Path('/proc/meminfo').read_text()
        record['causal_scope'] = 'Absence of the old OOM under current resources, not proof that Nemotron was the sole historical cause.'
        (OUT / 'gpu_memory_retry.json').write_text(json.dumps(record, indent=2) + '\n')

if __name__ == '__main__':
    main()
