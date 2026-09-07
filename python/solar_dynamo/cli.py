import argparse
from pathlib import Path
from . import Model,Parameters

def main():
    parser=argparse.ArgumentParser(description="Exact C++ port of surya_dynamo_v3.f")
    parser.add_argument("--init",type=Path,default=Path("init.dat"));parser.add_argument("--output-dir",type=Path,default=Path("output"))
    parser.add_argument("--steps",type=int,default=None);parser.add_argument("--tmax",type=float,default=10.0)
    parser.add_argument("--analytic-init",action="store_true",help="Use the original irelax=0 branch")
    args=parser.parse_args();p=Parameters();p.tmax=args.tmax;p.relaxed_initial_state=not args.analytic_init
    model=Model(p);model.initialize(args.init)
    if args.steps is None:model.run(args.output_dir)
    else:model.step_with_output(args.steps,args.output_dir);model.write_final_outputs(args.output_dir)
    print(f"completed {model.step_number} steps; t={model.time:.7f}; outputs in {args.output_dir}")
if __name__=="__main__":main()
