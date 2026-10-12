"""Acquire only preflighted native chunks, with resumable size/hash checks."""
from p0_t9_preflight import *
from numcodecs import Blosc
DATA=ROOT/'data/raw/p0_t9/native_chunks'
def main():
 budget=json.loads((OUT/'native_chunk_budget.json').read_text());assert budget['errors']==0 and budget['compressed_bytes']<4500000000
 assert json.loads((OUT/'runtime_check.json').read_text())['cuda_forward']=='PASS'
 DATA.mkdir(parents=True,exist_ok=True)
 dump('acquisition_plan.json',{'maximum_native_payload_bytes':4500000000,'measured_native_payload_bytes':budget['compressed_bytes'],'checkpoint_bytes':138360039,'no_other_CT_regions_authorized_in_this_run':True})
 def get(rec):
  key='/'.join(map(str,rec['key']));p=DATA/key.replace('/','_')
  if not p.exists():
   last=None
   for attempt in range(3):
    try:
     with urllib.request.urlopen(BASE+'0/'+key,timeout=120) as r:b=r.read()
     assert len(b)==rec['bytes'];assert len(Blosc().decode(b))==128**3
     q=p.with_suffix('.partial');q.write_bytes(b);q.rename(p);break
    except Exception as e:last=e
   else:raise last
  assert p.stat().st_size==rec['bytes']
  return {'key':rec['key'],'source':BASE+'0/'+key,'filename':str(p.relative_to(ROOT)),
          'expected_bytes':rec['bytes'],'actual_bytes':p.stat().st_size,
          'bytes':p.stat().st_size,'sha256':sha(p),'completion_status':'COMPLETE'}
 records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=48) as ex:
  for r in ex.map(get,budget['chunks']):
   records.append(r)
   if len(records)%250==0:print('acquired',len(records),flush=True)
 dump('acquisition.json',{'source':BASE,'chunks':records,'payload_bytes':sum(r['bytes'] for r in records)})
 print('Acquisition complete',flush=True)
if __name__=='__main__':main()
