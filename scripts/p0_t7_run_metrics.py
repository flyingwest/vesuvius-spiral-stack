"""Run the identical T6 metrics on a frozen T7 scalar or combination set."""
import argparse,json
import p0_t7_evaluator as evaluator
from p0_t7_family import OUT,FEATURES
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['scalars','combinations']);args=parser.parse_args()
 evaluator.OUT=OUT/args.phase
 if args.phase=='scalars':evaluator.FEATURES=FEATURES
 else:
  specs=json.loads((OUT/'combination_plan.json').read_text())['specifications']
  evaluator.FEATURES={k:(v['definition'],v['display_range']) for k,v in specs.items()}
 evaluator.main()
