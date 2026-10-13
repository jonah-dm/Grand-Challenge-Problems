import numpy as np, sys
sys.argv=['x']
exec(open('ex_eval.py').read().split("forms={")[0])
atoms=['H','Li','B','C','N','O','F','He','Be','Ne','Mg','Ar','Na','Si','P','S','Cl']
print('atom exact-Ex, then (Ex_sl-Ex_exact) in mEh for forms; p = weight exponent')
cfg=[('pbe',None,'chan',2)]+[('b',(MU,0.037711),zm,p) for zm in ('chan','tot') for p in (1,2,3,4,6,8)]
print('%-5s'%'',*['%s-%s'%(('pbe' if c[0]=='pbe' else 'H'),(c[2] if c[0]!='pbe' else '')+str(c[3])) for c in cfg])
tot=np.zeros(len(cfg))
for a in atoms:
    d=Dg[a]; row=[1000*(Ex_sl(d,c[0],c[1],c[2],c[3])-d['Ex']) for c in cfg]
    tot+=np.abs(row); print('%-5s'%a,*['%7.1f'%x for x in row])
print('MAE ',*['%7.1f'%x for x in tot/len(atoms)])
