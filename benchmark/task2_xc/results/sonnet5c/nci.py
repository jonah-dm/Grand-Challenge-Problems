import sys, os, json
sys.path.insert(0,'/home/user/Grand-Challenge-Problems/benchmark/task2_xc'); sys.path.insert(0,'/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit')
import score
from pyscf import lib
lib.num_threads(int(os.environ.get('NTHREADS','1')))
A=[('O',(-1.551007,-0.114520,0.0)),('H',(-1.934259,0.762503,0.0)),('H',(-0.599677,0.040712,0.0))]
B=[('O',(1.350625,0.111469,0.0)),('H',(1.680398,-0.373741,-0.758561)),('H',(1.680398,-0.373741,0.758561))]
sub=score.load_submission('/home/user/Grand-Challenge-Problems/benchmark/task2_xc/results/sonnet5c/submission.py')
basis=sys.argv[1] if len(sys.argv)>1 else 'def2-qzvp'
E={}
for nm,at in (('dimer',A+B),('A',A),('B',B)):
    E[nm]=score.run_species(at,0,0,basis,None,sub)[0]; print(nm,E[nm],flush=True)
print('De(kcal/mol)= %.2f   (CCSD(T)/CBS ref -4.99)'%((E['dimer']-E['A']-E['B'])*627.509))
