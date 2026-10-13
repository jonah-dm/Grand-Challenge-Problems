import numpy as np, pickle, sys
from collections import Counter
sys.path.insert(0,'/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit')
sys.path.insert(0,'.')
import diag_x
Dg=pickle.load(open('diag_x_def2-qzvp.pkl','rb'))
KAP=0.804; MU=0.2195149727645171; CX=-0.75*(3/np.pi)**(1/3)
def zsat(z,k=40.0): return 1-np.log(1+np.exp(k*(1-z)))/k
def ex_chan(r,sg,t,mode,p1,ztot=None,pw=2):
    m=r>1e-12
    out=np.zeros_like(r)
    r_,sg_,t_=r[m],sg[m],t[m]
    n2=2*r_; kf=(3*np.pi**2*n2)**(1/3); s2=4*sg_/(4*kf**2*n2**2)
    hg=KAP*(1-1/(1+MU*s2/KAP))
    if mode=='pbe': h=hg
    else:
        z=zsat(sg_/(8*r_*np.maximum(t_,1e-20))) if ztot is None else ztot[m]
        w=z**pw; a,b=p1; h1=a*s2/(1+b*s2*s2); h=(1-w)*hg+w*h1
    out[m]=0.5*CX*n2**(4/3)*(1+h)
    return out
def Ex_sl(d,mode,p1,zmode='chan',pw=2):
    ra,rb=d['ra'],d['rb']
    ga=ra[1:4]; gb=rb[1:4]
    saa=np.sum(ga**2,0); sbb=np.sum(gb**2,0); ta,tb=ra[-1],rb[-1]
    zt=None
    if zmode=='tot':
        n=ra[0]+rb[0]; sig=np.sum((ga+gb)**2,0); zt=zsat(sig/(8*np.maximum(n,1e-14)*np.maximum(ta+tb,1e-20)))
    e=ex_chan(ra[0],saa,ta,mode,p1,zt,pw)+ex_chan(rb[0],sbb,tb,mode,p1,zt,pw)
    return np.sum(d['w']*e)
forms={'PBE-form':('pbe',None,'chan'),'H-exact chan':('b',(MU,0.037711),'chan'),'H-exact tot':('b',(MU,0.037711),'tot'),'HHe chan':('b',(1.9843,16.6604),'chan'),'HHe tot':('b',(1.9843,16.6604),'tot')}
mols=[k for k in Dg if k not in diag_x.SP or True]
atoms={k:Dg[k] for k in ['H','Li','B','C','N','O','F','He','Be','Ne','Mg','Ar','Na','Si','P','S','Cl']}
comp={'H2':'HH','CH4':'CHHHH','H2O':'OHH','N2':'NN','HCl':'ClH','H2S':'SHH','PH3':'PHHH','SiH4':'SiHHHH','LiF':'LiF','BeH2':'BeHH','NaCl':'NaCl','Li2':'LiLi','BH':'BH'}
import re
def parse(s): return Counter(re.findall('[A-Z][a-z]?',s))
res={}
print('%-8s %9s'%('mol','exactAEx'),*['%13s'%f for f in forms])
tot={f:[] for f in forms}
for m,c in comp.items():
    cnt=parse(c)
    ae_ex=sum(n*Dg[a]['Ex'] for a,n in cnt.items())-Dg[m]['Ex']
    row=[]
    for f,(mode,p1,zm) in forms.items():
        ae=sum(n*Ex_sl(Dg[a],mode,p1,zm) for a,n in cnt.items())-Ex_sl(Dg[m],mode,p1,zm)
        row.append((ae-ae_ex)*627.509); tot[f].append(abs(row[-1]))
    print('%-8s %9.1f'%(m,ae_ex*627.509),*['%+13.1f'%x for x in row])
print('MAE',*['%13.1f'%np.mean(v) for v in tot.values()])
# atoms absolute
print('atoms: (Ex_sl-Ex_exact) mEh')
for a in atoms:
    print(a,*['%+7.1f'%(1000*(Ex_sl(atoms[a],*forms[f])-atoms[a]['Ex'])) for f in forms])
