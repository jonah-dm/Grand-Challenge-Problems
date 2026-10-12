import json,sys
fs=sys.argv[1:]
ds=[json.load(open(f)) for f in fs]
keys=list(ds[0]['items'])
print('%-16s'%'item'+''.join('%10s'%f.split('/')[-1][:9] for f in fs))
for k in keys:
    print('%-16s'%k+''.join('%10.1f'%d['items'][k]['error'] if k in d['items'] else '%10s'%'-' for d in ds))
for c in ['AE','IP']:
    print('%-16s'%(c+' MAE')+''.join('%10.2f'%d['summary'][c]['MAE'] for d in ds))
print('%-16s'%'WTMAD'+''.join('%10.2f'%d['WTMAD'] for d in ds))
print('%-16s'%'conv'+''.join('%10s'%d['all_converged'] for d in ds))
