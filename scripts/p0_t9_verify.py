"""Verify finalized stacks and replay fixed probes from frozen coordinates."""
from p0_t9_infer import *
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 torch.set_num_threads(4)
 prior=json.loads((OUT/'preservation.json').read_text());assert all(sha(ROOT/p)==h for p,h in prior.items())
 params=json.loads((OUT/'inference_parameters.json').read_text());assert params['model_spacing_um']==9.6 and params['batch_size']==8
 mask=np.load(OUT/'mask.npy');stack=np.load(OUT/'ink_probability.npy');assert np.array_equal(np.isnan(stack),mask)
 for name in ['mask','provenance','s_um','z_indices','z_relative_um']:assert sha(OUT/f'{name}.npy')==sha(ROOT/f'outputs/p0_t5/{name}.npy')
 rows=0;centers=0
 for z in range(128,178):
  row=np.load(OUT/f'slice_{z}_predictions.npz');original=np.load(ROOT/f'outputs/p0_t5/slice_{z}/path_and_samples.npz');frame=np.load(ROOT/f'outputs/p0_t6/slice_{z}_profiles.npz');ix=row['sample_indices']
  assert np.array_equal(ix,np.flatnonzero(~original['mask']));assert np.array_equal(row['xy'],original['xy'][ix]);assert np.array_equal(row['xy'],frame['xy']);assert np.array_equal(row['tangent'],frame['tangent']);assert np.array_equal(row['normal'],frame['normal'])
  assert np.array_equal(stack[z-128,ix],row['probability']);rows+=1;centers+=len(ix)
 c=configure_model(SimpleNamespace(checkpoint=ROOT/'data/raw/p0_t9/step-075000.pth',amp_dtype='auto'));c.model.cuda().eval();replay=[]
 for file in sorted(OUT.glob('probe_*.npz')):
  z=int(file.stem.split('_')[1]);s=np.load(file);raw=patch(z,s['xy'],s['tangent'],s['normal']);assert np.array_equal(raw,s['raw_patch'])
  cases=[raw,raw[::-1].copy(),patch(z,s['xy'],s['tangent'],s['normal'],-2),patch(z,s['xy'],s['tangent'],s['normal'],2)]
  x=np.stack([normalize_flat_patch(a,c.preprocessing) for a in cases])[:,None]
  with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):q=logits_to_probabilities(c.model(torch.from_numpy(x).cuda()),image_hw=(128,128))[:,0,63:65,63:65].mean((1,2)).cpu().numpy()
  err=float(np.max(abs(q-s['probabilities'])));assert err<=.002,(file,err);replay.append({'file':file.name,'raw_patch_exact':True,'maximum_probability_error':err})
 assert len(replay)==96
 result={'status':'PASS','rows_verified':rows,'all_valid_coordinate_mappings_verified':centers,'mask_provenance_physical_coordinates_byte_identical':True,'prior_files_preserved':len(prior),'probe_raw_CT_patches_replayed_exact':len(replay),'model_predictions_replayed':len(replay)*4,'probability_replay_tolerance':.002,'maximum_probability_replay_error':max(a['maximum_probability_error'] for a in replay),'replay_scope':'all50row assemblies and130629coordinate mappings;96fixed real-CT probes and their3orientation/depth variants, not a second full130629center inference','real_ct_smoke':json.loads((OUT/'smoke_test.json').read_text()),'replay':replay}
 (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['replay','real_ct_smoke']},indent=2),flush=True)
if __name__=='__main__':main()
