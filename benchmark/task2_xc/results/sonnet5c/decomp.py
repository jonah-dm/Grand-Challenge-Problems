import numpy as np, pickle, sys, json, os, importlib
sys.path.insert(0,'/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit'); sys.path.insert(0,'.')
import dev_set
Dg=pickle.load(open('diag_all_def2-qzvp.pkl','rb'))
ref=json.load(open('/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit/dev_reference.json'))['items']
H2K=627.5094740631
def load(env):
    for k in list(os.environ):
        if k.startswith('S5_'): del os.environ[k]
    os.environ.update(env)
    import submission as s; importlib.reload(s); return s
def corr_E(s,d):
    ra,rb=d['ra'],d['rb']; ga=ra[1:4]; gb=rb[1:4]
    args=[np.maximum(ra[0],1e-14),np.maximum(rb[0],1e-14),np.sum(ga**2,0),np.sum(ga*gb,0),np.sum(gb**2,0),np.maximum(ra[-1],1e-20),np.maximum(rb[-1],1e-20)]
    m=(args[0]+args[1])>1e-11
    old=s.HYB; s.HYB=1.0
    e=np.zeros_like(args[0]); e[m]=s._E(*[a[m] for a in args]); s.HYB=old
    return np.sum(d['w']*e)
def exch_E(s,d):
    ra,rb=d['ra'],d['rb']; ga=ra[1:4]; gb=rb[1:4]
    a=[np.maximum(ra[0],1e-14),np.maximum(rb[0],1e-14),np.sum(ga**2,0),np.sum(ga*gb,0),np.sum(gb**2,0),np.maximum(ra[-1],1e-20),np.maximum(rb[-1],1e-20)]
    m=(a[0]+a[1])>1e-11
    n=a[0]+a[1]; sig=a[2]+2*a[3]+a[4]
    zt=s._zsat(sig/(8*n*(a[5]+a[6]))) if s.ZX=='tot' else None
    e=np.zeros_like(a[0])
    e[m]=s._ex_chan(a[0][m],a[2][m],a[5][m],None if zt is None else zt[m])+s._ex_chan(a[1][m],a[4][m],a[6][m],None if zt is None else zt[m])
    return np.sum(d['w']*e)
def run(label,env):
    s=load(env)
    cache={}
    def E(name):
        if name not in cache:
            d=Dg[name]; cache[name]=(exch_E(s,d),corr_E(s,d))
        return cache[name]
    rows=[]
    for k,(cls,expr) in dev_set.ITEMS.items():
        if cls!='AE': continue
        mol=[x for x in expr if not x.endswith('_atom')][0]
        ex_err=0;c_sl=0;ex_hf=0;hf=0
        for sp,c in expr.items():
            exs,cs=E(sp); ex_err+=c*(0.75*exs+0.25*Dg[sp]['Ex']-Dg[sp]['Ex']); c_sl+=c*cs; hf+=c*Dg[sp]['E']
        ae_ref=ref[k]['kcal']; ae_hf=hf*H2K
        ec_ref=ae_ref-ae_hf       # exact AE correlation part (CCSD(T) minus HF)
        ec_err=c_sl*H2K-ec_ref
        rows.append((k,ex_err*H2K,ec_err,ex_err*H2K+ec_err))
    print(label)
    for r in rows: print('  %-14s xerr %+7.1f  cerr %+7.1f  total %+7.1f'%r)
    print('  MAE total %.2f  xerr %.2f cerr %.2f'%(np.mean([abs(r[3]) for r in rows]),np.mean([abs(r[1]) for r in rows]),np.mean([abs(r[2]) for r in rows])))
if __name__=='__main__':
    run('c1 none/pair',{'S5_X1':'none','S5_CORR':'pair'})
    run('none/pbe-corr',{'S5_X1':'none','S5_CORR':'pbe'})
    run('H tot/pair',{'S5_X1':'H','S5_ZX':'tot'})
