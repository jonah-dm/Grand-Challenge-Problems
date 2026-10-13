import sys, time
sys.path.insert(0,'.'); sys.path.insert(0,'kit')
import score
sub=score.load_submission('results/opus5t/submission.py')
atoms=[('H',(0,0,-0.93)),('H',(0,0,0)),('H',(0,0,0.93))]
for bas in sys.argv[1:]:
    e,conv,dt=score.run_species(atoms,0,1,bas,None,sub); print(bas,e,conv,'%.0fs'%dt,flush=True)
