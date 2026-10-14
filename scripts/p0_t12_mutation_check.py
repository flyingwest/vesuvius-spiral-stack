"""Check that synthetic tests actually catch three critical guard regressions."""
from pathlib import Path
import tempfile,subprocess,os,json
ROOT=Path(__file__).resolve().parents[1]
def main():
    original=(ROOT/'scripts/p0_t12_core.py').read_text();tests=(ROOT/'scripts/p0_t12_synthetic.py').read_text();results=[]
    mutations=[('displacement_guard','if distance>limit:', 'if False:'),('competitor_guard','if len(ordered)>1 and ', 'if False and len(ordered)>1 and '),('quadratic_branch',"quad['accepted']=bool(all(quad['checks'].values()))","quad['accepted']=False")]
    for name,old,new in mutations:
        assert original.count(old)==1
        with tempfile.TemporaryDirectory(prefix='t12-mutant-',dir='/home/wpotech/.hermes/cache/scratch') as directory:
            p=Path(directory);(p/'p0_t12_core.py').write_text(original.replace(old,new));(p/'check.py').write_text(tests)
            env={**os.environ,'PYTHONPATH':str(ROOT/'scripts'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
            r=subprocess.run([str(ROOT/'.venv/bin/python'),'-B',str(p/'check.py')],env=env,capture_output=True,text=True)
            assert r.returncode!=0 and 'AssertionError' in r.stderr,(name,r.stdout,r.stderr)
            results.append({'mutant':name,'killed':True})
    (ROOT/'outputs/p0_t12/synthetic_checks.json').write_text(json.dumps({'synthetic_baseline':'PASS','guard_mutants':results},indent=2)+'\n');print('PASS: all3guard mutants killed')
if __name__=='__main__':main()
