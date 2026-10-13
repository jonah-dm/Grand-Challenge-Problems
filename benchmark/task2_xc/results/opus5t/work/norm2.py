import numpy as np
exec(open('norm.py').read().split("c=0.3125")[0])
MU=0.2195149727645171
Fq=lambda s,c: (1+MU*s*s)/(1+MU*s**4/c)
g=lambda c: Ex(lambda s: Fq(s,c))+0.3125
for c in (0.5,1,2,5,10,50): print(c, g(c))
