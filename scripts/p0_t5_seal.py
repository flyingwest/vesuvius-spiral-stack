"""Seal complete T5 evidence and immediately verify every checksum."""
import json,subprocess
from pathlib import Path
from p0_t5_consecutive import ROOT,OUT,sha,dump
from p0_t5_evaluate import barrier

def main():
    manifest=OUT/'SHA256SUMS';assert not manifest.exists(),'Already sealed; use sha256sum -c for verification'
    barrier()
    for name in ['verification.json','geometry_replay_validation.json','diagnostic_validation.json']:
        assert json.loads((OUT/name).read_text())['status']=='PASS'
    for name,h in json.loads((OUT/'preservation.json').read_text()).items():assert sha(ROOT/name)==h,name
    required=['raw_stack.png','raw_stack_contrast.png','physical_aspect_stack.png','mask_stack.png','provenance_stack.png','coherence_diagnostic.png','lag_diagnostic.png','local_window_coherence.png','adjacent_correlations.csv','lag_correlations.csv','path_metrics.csv','row_extents.csv','stack.npy','mask.npy','provenance.npy']
    assert all((OUT/f).is_file() for f in required)
    paths=sorted(p for p in OUT.rglob('*') if p.is_file())
    rel=[str(p.relative_to(ROOT)) for p in paths]
    probe=subprocess.run(['git','check-ignore','--stdin'],input='\n'.join(rel)+'\n',capture_output=True,text=True,cwd=ROOT)
    assert probe.returncode in (0,1);ignored=set(probe.stdout.splitlines())
    inventory=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'git':'excluded; retained locally' if str(p.relative_to(ROOT)) in ignored else 'milestone inclusion'} for p in paths]
    dump(OUT/'git_artifact_inventory.json',{'scope':'existing T5 files before this inventory/receipt/manifest; those three additional metadata files are included in Git','files':inventory,
        'included_bytes':sum(r['bytes'] for r in inventory if r['git']=='milestone inclusion'),'excluded_local_bytes':sum(r['bytes'] for r in inventory if r['git']!='milestone inclusion'),
        'note':'Git clone is not full replay package; all excluded artifacts remain local and are sealed by SHA256SUMS'})
    dump(OUT/'seal_receipt.json',{'P0_T5_COMPLETE':'YES','P0_T5_RESULT':'PARTIAL','READY_FOR_NEXT_STEP':'NO','ready_meaning':'not recommended for unchanged larger-Z scale-up',
        'all50_finalized_hashes':'PASS','prior_T1_through_T4_and_protected_files':'PASS','numeric_and_saved_geometry_replay':'PASS','lag_diagnostic_validation':'PASS','no_reconstruction_rerun_during_finalization':True,
        'checksum_verification':'seal script requires every entry to pass immediately after writing; command output is final confirmation'})
    files=sorted([p for p in OUT.rglob('*') if p.is_file() and p!=manifest]+list((ROOT/'scripts').glob('p0_t5_*.py'))+[ROOT/'docs/p0_t5_consecutive_thread_stack.md'])
    manifest.write_text(''.join(sha(p)+'  '+str(p.relative_to(ROOT))+'\n' for p in files))
    result=subprocess.run(['sha256sum','-c',str(manifest)],cwd=ROOT,capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    print('PASS: every checksum verified:',len(files),'entries; prior preservation; finalized rows; required outputs; all validation receipts')
    print('Git payload bytes',sum(r['bytes'] for r in inventory if r['git']=='milestone inclusion'),'excluded local bytes',sum(r['bytes'] for r in inventory if r['git']!='milestone inclusion'))
if __name__=='__main__':main()
