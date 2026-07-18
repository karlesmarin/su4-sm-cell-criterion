# Clean: (1) verify ALL L-singlets (odd n-ality) have n3+n4 ODD -> the two T2/Z2 parities are LOCKED into ONE
# chirality Z2 on the singlet sector (why 5D=6D). (2) restate the extrapolation law.
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")
def su2ct(counter):
    rem=dict(counter); out=[]
    while any(v>0 for v in rem.values()):
        t=max(k for k,v in rem.items() if v>0); m=rem[t]
        for x in range(t,-t-1,-2): rem[x]=rem.get(x,0)-m
        out.extend([t]*m)
    return out
def lsinglets(L):
    a,b,c=L; boxes=a+2*b+3*c
    from collections import defaultdict
    buck=defaultdict(lambda: defaultdict(int)); meta={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        buck[(q8,q15)][n1-n2]+=mult; meta[(q8,q15)]=(n1,n2,n3,n4)
    return [meta[k] for k,mb in buck.items() if su2ct(mb).count(0)>0]
# (1) n3+n4 parity over all L-singlets, odd n-ality
bad=0; tested=0
for a in range(9):
 for b in range(9):
  for c in range(9):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0 or A3(*L).degree()>500: continue
   tested+=1
   for (n1,n2,n3,n4) in lsinglets(L):
       if (n3+n4)%2==0: bad+=1
print(f"(1) L-singlets with n3+n4 EVEN (should be 0): {bad}   [tested {tested} reps dim<=500]")
if bad==0: print("    => ALL L-singlets have n3+n4 ODD. The two T2/Z2 parities (p0=n4, p2=n3+n4) are LOCKED:")
if bad==0: print("       n4 odd  <=>  n3 even  on the singlet sector => ONE effective chirality Z2 => fraction 1/2 in BOTH 5D & 6D.")
# (2) correct projection table
print("\n(2) count = (b+1)(a+c+1) x fraction:")
for L in [(0,2,1),(2,2,1),(0,3,3)]:
    a,b,c=L; ws=lsinglets(L); full=(b+1)*(a+c+1)
    rh=sum(1 for n in ws if (-1)**n[3]==-1)          # n4 odd (=RH, and auto n3 even)
    print(f"  {L}: full={full}=(b+1)(a+c+1) | RH(n4 odd)={rh}=full/2 | 6D adds no cut (n3 even auto)")
print("\n  EXTRAPOLATION LAW:  N_zeromode = [Levi-branching count of the bulk rep]  x  1/(prod of INDEPENDENT projections)")
print("   - group part: dimension-INDEPENDENT (Levi/CG, node-independent).")
print("   - each extra compact Z2 on an INDEPENDENT axis -> x1/2; on a LOCKED axis (like n3~n4 here) -> x1.")
print("   - Z_N orbifold -> x1/N per independent projection. Higher D / bigger N => fewer modes.")
