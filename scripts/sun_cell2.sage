# Fixed: hypercharge Y must be TRACELESS: Y . tau = 0, tau=(2,1,...,1) in c-coords (n=(1,..,1)).
# Cell closes <=> exists doublet Q, singlets u!=d, traceless Y with Y(Q)=1/6,Y(u)=2/3,Y(d)=-1/3.
# Tag: free-Y params = N-2 (=rank-1); 3 constraints; symmetric cell makes it effectively 2. Rigid iff N-2<=2 i.e. N<=4.
# author: Carles Marin (Claude as AI assistant).
def make(N): return WeylCharacterRing(f"A{N-1}", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2): rem[m2]=rem.get(m2,0)-m
        out.extend([top+1]*m)
    return out
def modes(ring,N,L,parity):
    boxes=sum((i+1)*L[i] for i in range(N-1)); buckets={}
    for w,mult in ring(*L).weight_multiplicities().items():
        n=[ZZ(w[i]+QQ(boxes)/N) for i in range(N)]
        p0=(-1)**n[N-1]; p2=(-1)**(n[N-2]+n[N-1])
        if (p0,p2)!=parity: continue
        c=tuple([n[0]+n[1]]+[n[i] for i in range(2,N)])
        buckets.setdefault(c,{})[n[0]-n[1]]=buckets.setdefault(c,{}).get(n[0]-n[1],0)+mult
    out=[]
    for c,bym in buckets.items():
        for d in su2_decompose(bym): out.append((d,c))
    return out
def cell(ring,N,L):
    D=list(set(c for (d,c) in modes(ring,N,L,(1,1)) if d==2))
    S=list(set(c for (d,c) in modes(ring,N,L,(-1,-1)) if d==1))
    tau=vector(QQ,[2]+[1]*(N-2))
    for Q in D:
        for u in S:
            for d in S:
                if d==u: continue
                A=Matrix(QQ,[list(Q),list(u),list(d),list(tau)]); bb=vector(QQ,[QQ(1)/6,QQ(2)/3,QQ(-1)/3,0])
                Aug=A.augment(bb)
                if A.rank()==Aug.rank() and A.rank()>=3:   # solvable, Y nontrivially determined
                    return (Q,u,d)
    return None
def nality(N,L): return sum((i+1)*L[i] for i in range(N-1))%N

# VALIDATE SU(4)
r4=make(4); bad=0;n=0
for a in range(6):
 for b in range(6):
  for c in range(6):
   L=(a,b,c)
   if L==(0,0,0) or r4(*L).degree()>200: continue
   n+=1
   gate=((a+2*b+3*c)%2==1) and b>=1 and (a+b+c)>=3
   if (cell(r4,4,L) is not None)!=gate: bad+=1
print(f"SU(4) validation: tested {n}, mismatches vs 3-gate = {bad}  {'VALID' if bad==0 else 'BROKEN'}")

# SU(5): scan, report which reps admit + their center charge (mod 5) + size
print("\n=== SU(5): which irreps admit the cell? (dim<=200) ===")
r5=make(5); adm5=[]
for L in Tuples(range(5),4):
   L=tuple(L)
   if sum(L)==0: continue
   if r5(*L).degree()>200: continue
   if cell(r5,5,L) is not None: adm5.append((r5(*L).degree(),L,nality(5,L),sum(L)))
for dim,L,na,s in sorted(adm5)[:25]:
    print(f"  {L} dim={dim} nality(mod5)={na} sum={s}")
print(f"  total admitting SU(5) reps (dim<=200): {len(adm5)}")
print(f"  center charges (mod5) of admitting: {sorted(set(x[2] for x in adm5))}")

# SU(5): full stats admit vs not
print("\n=== SU(5) full stats (dim<=160) ===")
tot5=0;a5=0
for L in Tuples(range(5),4):
   L=tuple(L)
   if sum(L)==0 or r5(*L).degree()>160: continue
   tot5+=1
   if cell(r5,5,L) is not None: a5+=1
print(f"  SU(5): {a5}/{tot5} irreps admit  ({round(100*a5/tot5)}%)")

# SU(6)
print("\n=== SU(6) stats (dim<=252) ===")
r6=make(6); tot6=0;a6=0;cent=set()
for L in Tuples(range(4),5):
   L=tuple(L)
   if sum(L)==0 or r6(*L).degree()>252: continue
   tot6+=1
   if cell(r6,6,L) is not None:
       a6+=1; cent.add(sum((i+1)*L[i] for i in range(5))%6)
print(f"  SU(6): {a6}/{tot6} irreps admit  ({round(100*a6/tot6) if tot6 else 0}%)")
print(f"  center charges (mod6) of admitting: {sorted(cent)}")

# The mechanism, stated: free-Y params = N-2 ; effective cell constraints = 2 (symmetric cell) ; rigid iff N-2<=2 <=> N<=4
print("\n=== freedom-counting mechanism ===")
for N in [4,5,6]:
    print(f"  SU({N}): rank={N-1}, hypercharge-free-params=N-2={N-2}, cell-constraints(eff)=2 -> {'RIGID gate (SU(4) unique)' if N-2<=2 else 'DISSOLVED (slack '+str(N-2-2)+')'}")

print("\n\n=== MECHANISM TEST 1: restrict SU(5) hypercharge to 2 free params -> does a rigid gate reappear? ===")
# Force Y to lie in a 2-dim subspace of c-functionals (drop the extra U(1)): require Y . e_extra = 0 for one extra direction.
# Implement: add EXTRA constraints Y.k=0 for (rank-1-2) directions k, shrinking Y to 2 params, then test admission.
def cell_restricted(ring,N,L,nparams):
    D=list(set(c for (d,c) in modes(ring,N,L,(1,1)) if d==2))
    S=list(set(c for (d,c) in modes(ring,N,L,(-1,-1)) if d==1))
    tau=vector(QQ,[2]+[1]*(N-2))
    # extra kill-directions to shrink Y to nparams: use standard basis e_3.. as needed
    extra=[vector(QQ,[0]*(i)+[1]+[0]*(N-1-i-1)) for i in range(2, 2+((N-2)-nparams))]  # kill (N-2-nparams) dirs
    for Q in D:
        for u in S:
            for d in S:
                if d==u: continue
                rows=[list(Q),list(u),list(d),list(tau)]+[list(e) for e in extra]
                rhs=[QQ(1)/6,QQ(2)/3,QQ(-1)/3,0]+[0]*len(extra)
                A=Matrix(QQ,rows); bb=vector(QQ,rhs)
                if A.rank()==A.augment(bb).rank() and A.rank()>=3:
                    return True
    return False
r5=make(5)
tot=0; adm_full=0; adm_2p=0; cent2=set()
for L in Tuples(range(5),4):
   L=tuple(L)
   if sum(L)==0 or r5(*L).degree()>160: continue
   tot+=1
   if cell(r5,5,L) is not None: adm_full+=1
   if cell_restricted(r5,5,L,2):
       adm_2p+=1; cent2.add(sum((i+1)*L[i] for i in range(4))%5)
print(f"  SU(5) full (3 params): {adm_full}/{tot} admit")
print(f"  SU(5) restricted to 2 params: {adm_2p}/{tot} admit   centers={sorted(cent2)}")
print(f"  -> if restricted is MUCH sparser & center-selective, the mechanism is the PARAM COUNT.")

print("\n=== MECHANISM TEST 2: do SU(5) NON-admitting reps simply lack a doublet or 2 indep singlets? ===")
nogate=0; realgap=[]
for L in Tuples(range(5),4):
   L=tuple(L)
   if sum(L)==0 or r5(*L).degree()>160: continue
   if cell(r5,5,L) is not None: continue
   D=set(c for (d,c) in modes(r5,5,L,(1,1)) if d==2)
   S=list(set(c for (d,c) in modes(r5,5,L,(-1,-1)) if d==1))
   # enough modes? need a doublet and 2 independent singlets
   indep2 = any(Matrix(QQ,[list(S[i]),list(S[j])]).rank()==2 for i in range(len(S)) for j in range(i+1,len(S)))
   enough = (len(D)>=1 and indep2)
   if not enough: nogate+=1
   else: realgap.append((L,len(D),len(S)))
print(f"  SU(5) non-admitting that simply LACK modes (doublet or 2 indep singlets): {nogate}")
print(f"  SU(5) non-admitting DESPITE having modes (a real residual gate?): {len(realgap)}")
for x in realgap[:10]: print("     residual:",x)
