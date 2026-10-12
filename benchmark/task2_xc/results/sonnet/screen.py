import sys, os, json
sys.path.insert(0,'kit'); sys.path.insert(0,'.')
import dev_set as S, score
from pyscf import lib
lib.num_threads(1)
xc=sys.argv[1]; basis=sys.argv[2]
ref=json.load(open('kit/dev_reference.json'))['items']
E={}
for n,(a,c,s) in S.SPECIES.items():
    E[n]=score.run_species(a,c,s,basis,xc,None)[0]
ae=[];ip=[]
for k,(cls,ex) in S.ITEMS.items():
    v=sum(c*E[s] for s,c in ex.items())*627.5094740631
    (ae if cls=='AE' else ip).append((k,v-ref[k]['kcal'],abs(ref[k]['kcal'])))
def sc(l): return sum(abs(e) for _,e,_ in l)/len(l)/(sum(r for *_,r in l)/len(l))*100
w=(sc(ae)*len(ae)+sc(ip)*len(ip))/(len(ae)+len(ip))
print(xc,basis,'AEmae %.2f IPmae %.2f WTMAD %.2f'%(sum(abs(e) for _,e,_ in ae)/len(ae),sum(abs(e) for _,e,_ in ip)/len(ip),w),' IPerr',[round(e,1) for _,e,_ in ip],flush=True)
json.dump(E,open('results/sonnet/E_%s_%s.json'%(xc,basis),'w'))
