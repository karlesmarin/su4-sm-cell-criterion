# Confirm: (1) J=n1+n2-n3-n4 changes ONLY under alpha2 (middle root).
#          (2) spread of (2q8+q15) over RH singlets == 0  <=>  b==0, and grows with b.
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2): rem[m2]=rem.get(m2,0)-m
        out.extend([top+1]*m)
    return out
def singlet_J_spread(L):
    boxes=sum((i+1)*L[i] for i in range(3)); buckets={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        if (-1)**n4!=-1 or (-1)**(n3+n4)!=-1: continue   # RH
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        bb=buckets.setdefault((q8,q15),{}); bb[n1-n2]=bb.get(n1-n2,0)+mult
    Jvals=[]
    for (q8,q15),bym in buckets.items():
        for d in su2_decompose(bym):
            if d==1: Jvals.append(2*q8+q15)
    if not Jvals: return None
    return max(Jvals)-min(Jvals)

# (1) alpha-move check on raw weights: does 2q8+q15 shift by 0,+-6 under alpha1,alpha3 and by +-6 under alpha2?
# (2q8+q15)=3J; alpha2 shifts J by -2 => shifts (2q8+q15) by -6. alpha1,alpha3 shift by 0.
print("=== (1) root-shift of (2q8+q15)=3J: e_i-e_{i+1} lowering ===")
def qc(n): 
    n1,n2,n3,n4=n; return 2*(n1+n2-2*n3)+(n1+n2+n3-3*n4)
base=(5,3,2,1)
print(f"  base {base}: 2q8+q15={qc(base)}")
for name,d in [("alpha1 (e1-e2)",(-1,1,0,0)),("alpha2 (e2-e3)",(0,-1,1,0)),("alpha3 (e3-e4)",(0,0,-1,1))]:
    nn=tuple(base[i]+d[i] for i in range(4))
    print(f"    after {name}: {nn}  2q8+q15={qc(nn)}  delta={qc(nn)-qc(base)}")

# (2) spread of (2q8+q15) over RH singlets vs b
print("\n=== (2) spread of (2q8+q15) over RH singlets, grouped by b (odd n-ality, dim<=320) ===")
import collections
byb=collections.defaultdict(set)
zerob_nonzero_spread=0; nonzero_b_zero_spread=[]
for a in range(11):
 for b in range(11):
  for c in range(11):
   L=(a,b,c)
   if L==(0,0,0): continue
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>320: continue
   sp=singlet_J_spread(L)
   if sp is None: continue
   byb[b].add(sp)
   if b==0 and sp!=0: zerob_nonzero_spread+=1
   if b>=1 and sp==0: nonzero_b_zero_spread.append((L,sum(L)))
for b in sorted(byb):
    print(f"  b={b}: spreads observed = {sorted(byb[b])}")
print(f"\n  b=0 with NONZERO spread (should be 0): {zerob_nonzero_spread}")
print(f"  b>=1 with ZERO spread (these are the a+b+c=2 rank-1 cases): {nonzero_b_zero_spread}")
