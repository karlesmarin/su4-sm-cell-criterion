# C4 test: is SU(4) the UNIQUE rank-3 simple group hosting the SM cell? (center Z4 vs Z2 the hinge?)
# Core group-theory question (no orbifold): does an irrep contain an SU(2)_L-doublet at Y=1/6 and
# SU(2)_L-singlets at Y=2/3,-1/3 under ONE free SU(2)_L-invariant Cartan functional Y (Y.alpha_L=0)?
# rank 3 => Y has 2 free params, 3 constraints => over-determined (rigid) for ALL three; center decides satisfiability.
# author: Carles Marin (Claude as AI assistant).
from collections import defaultdict
def cell_exists(ct, node, dimcap, maxlabel):
    R=WeylCharacterRing(ct, style="coroots"); amb=R.space()
    idx=list(amb.index_set()); aL=amb.simple_roots()[idx[node]]; aLv=amb.simple_coroots()[idx[node]]
    import itertools, math
    def strings(rep):
        # group weights by projection orthogonal to aL; within each, the aL-weights form SU(2) strings.
        wm=rep.weight_multiplicities()
        by=defaultdict(dict)
        for w,mult in wm.items():
            k=ZZ(w.scalar(aLv))                       # aL-weight (2*sz)
            key=tuple(w-QQ(k)/2*aL) if False else tuple((w).to_vector())  # keep full w; group by (w mod aL)
            # projection modulo aL: subtract (k/<aL,aLv>)*aL -> component along aL removed
            proj=w - QQ(w.scalar(aLv))/ZZ(aL.scalar(aLv))*aL
            by[tuple(proj.to_vector())][k]=by[tuple(proj.to_vector())].get(k,0)+mult
        D=[]; S=[]   # doublets (represented by proj vector), singlets
        for pj,km in by.items():
            # SU(2) decompose km (weights k, step 2)
            rem=dict(km)
            while any(v>0 for v in rem.values()):
                top=max(x for x,v in rem.items() if v>0); m=rem[top]
                for x in range(top,-top-1,-2): rem[x]=rem.get(x,0)-m
                if top==1:  D.append(vector(QQ,pj))
                if top==0:  S.append(vector(QQ,pj))
        return D,S
    aLvec=vector(QQ,aL.to_vector())
    for lab in itertools.product(range(maxlabel), repeat=len(idx)):
        if sum(lab)==0: continue
        rep=R(*lab)
        if rep.degree()>dimcap: continue
        D,S=strings(rep)
        for Q in D:
            for u in S:
                for d in S:
                    if tuple(u)==tuple(d): continue
                    M=Matrix(QQ,[list(Q),list(u),list(d),list(aLvec)])
                    b=vector(QQ,[QQ(1)/6,QQ(2)/3,QQ(-1)/3,0])
                    if M.rank()==M.augment(b).rank() and M.rank()>=3:
                        return (lab,tuple(Q),tuple(u),tuple(d))
    return None
print("=== C4: does the SM cell exist (free Y, no orbifold) in rank-3 simple groups? ===")
for ct in ["A3","B3","C3"]:
    R=WeylCharacterRing(ct); nnodes=R.cartan_type().rank()
    found=None
    for node in range(nnodes):
        r=cell_exists(ct, node, 200, 5)
        if r: found=(node,r); break
    print(f"  {ct}: cell found = {'YES node '+str(found[0])+' rep '+str(found[1][0]) if found else 'NO (dim<=200,all nodes)'}")
print("\n  A3=SU(4) center Z4 ; B3=SO(7),C3=Sp(6) center Z2. If only A3 -> center (not rank) is the selector (C4 survives).")
