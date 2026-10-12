"""Official architecture/weight compatibility and resource smoke check, not CT evidence."""
import json,time,hashlib,traceback
from pathlib import Path
from types import SimpleNamespace
import torch
OUT=Path('/work/outputs/p0_t9');checkpoint=Path('/work/data/raw/p0_t9/step-075000.pth')
from vesuvius.ink_detection.inference.infer import configure_model,logits_to_probabilities
from vesuvius.ink_detection.models.checkpoint import load_checkpoint,select_inference_weights
from vesuvius.ink_detection.config import InkConfig
from vesuvius.ink_detection.models.model import make_model
payload=load_checkpoint(checkpoint);config=InkConfig.from_mapping(payload['config']);which,state=select_inference_weights(payload,source=checkpoint)
model=make_model(config);keys=model.load_state_dict(state,strict=True);del model
meta={'model':'scrollprize/ink_9um/hybrid_3d2d-seed42/step-075000.pth','revision':'7109667e2607db1b90c37c8b09cb876ea7fe7bb1','sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),'selected_weights':which,'strict_weights_match':True,'checkpoint_config':payload['config'],'torch':torch.__version__}
(OUT/'model_metadata.json').write_text(json.dumps(meta,indent=2,default=str)+'\n')
c=configure_model(SimpleNamespace(checkpoint=checkpoint,amp_dtype='auto'))
r={'official_import':'PASS','strict_state_dict':'PASS','patch_size':c.patch_size,'input_depth':c.input_depth,'preprocessing':c.preprocessing,'smoke_input':'zeros only for runtime test; never used as scientific CT evidence'}
try:
 c.model.cuda();torch.cuda.synchronize();t=time.perf_counter()
 with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):
  for _ in range(3):out=logits_to_probabilities(c.model(torch.zeros((1,1,c.input_depth,c.patch_size,c.patch_size),device='cuda')),image_hw=(c.patch_size,c.patch_size))
 torch.cuda.synchronize();r.update(cuda_forward='PASS',three_forward_seconds=time.perf_counter()-t,output_shape=list(out.shape),peak_allocated_bytes=torch.cuda.max_memory_allocated())
except Exception as e:r.update(cuda_forward='FAIL',error=repr(e),traceback=traceback.format_exc())
(OUT/'runtime_check.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)
