# Cell theorem: reformulate SM-admissibility as an exact affine/rank condition on the weight diagram.
# author: Carles Marin (Claude as AI assistant). Attack on gates 2-3 = the isospin-cell theorem.
A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2): rem[m2]=rem.get(m2,0)-m
        out.extend([top+1]*m)
    return out
def modes_of(labels,parity):
    boxes=sum((i+1)*labels[i] for i in range(3)); buckets={}
    for w,mult in A3(*labels).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        p0=(-1)**n4; p2=(-1)**(n3+n4)
        if (p0,p2)!=parity: continue
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        b=buckets.setdefault((q8,q15),{}); b[n1-n2]=b.get(n1-n2,0)+mult
    out=[]
    for (q8,q15),bym in buckets.items():
        for d in su2_decompose(bym): out.append((d,q8,q15))
    return out

def doublets_singlets(L):
    D=sorted(set((m[1],m[2]) for m in modes_of(L,(1,1)) if m[0]==2))
    S=sorted(set((m[1],m[2]) for m in modes_of(L,(-1,-1)) if m[0]==1))
    return D,S

# --- Reformulation: cell closes iff exists doublet Q, singlets u!=d with a single Y=a*q8+b*q15
#     s.t. Y(Q)=1/6, Y(u)=2/3, Y(d)=-1/3.  <=>  3x3 det = 0 with top 2x2 minor != 0.
def cell_by_det(L):
    D,S=doublets_singlets(L)
    for Q in D:
        for u in S:
            M2=Matrix(QQ,[[Q[0],Q[1]],[u[0],u[1]]])
            if M2.det()==0: continue
            for d in S:
                if d==u: continue
                M3=Matrix(QQ,[[Q[0],Q[1],QQ(1)/6],[u[0],u[1],QQ(2)/3],[d[0],d[1],QQ(-1)/3]])
                if M3.det()==0:
                    return (Q,u,d)
    return None

def admits(L): return cell_by_det(L) is not None

# 1) VERIFY the det reformulation matches the original admits over a big range
print("=== 1) det=0 reformulation vs baseline (recomputed) — sanity ===")
def admits_base(L):
    D,S=doublets_singlets(L)
    for Q in D:
        for u in S:
            M=Matrix(QQ,[[Q[0],Q[1]],[u[0],u[1]]])
            if M.det()==0: continue
            al,be=M.solve_right(vector(QQ,[QQ(1)/6,QQ(2)/3]))
            for d in S:
                if d!=u and al*d[0]+be*d[1]==QQ(-1)/3: return True
    return False
mis=0;tot=0
for a in range(8):
 for b in range(8):
  for c in range(8):
   L=(a,b,c)
   if L==(0,0,0):continue
   if A3(*L).degree()>200:continue
   tot+=1
   if admits(L)!=admits_base(L): mis+=1; print("  MISMATCH",L)
print(f"  tested {tot}, mismatches {mis}  -> det=0 reformulation is EXACT" if mis==0 else "  BROKEN")

# 2) GATE 2 dissection: for b=0 reps WITH doublets, WHY does det=0 never hold?
print("\n=== 2) GATE 2 (b=0) obstruction geometry ===")
for L in [(2,0,1),(1,0,2),(4,0,1),(2,0,3),(1,0,4)]:
    if (sum((i+1)*L[i] for i in range(3)))%2==0: continue
    D,S=doublets_singlets(L)
    print(f"\n {L} dim={A3(*L).degree()}: #dbl={len(D)} #sng={len(S)}")
    print(f"   doublets(q8,q15): {D}")
    print(f"   singlets(q8,q15): {S}")
    # key geometric question: are all doublets & singlets collinear through origin (rank of the point set)?
    pts=[vector(QQ,p) for p in D+S]
    R=Matrix(QQ,pts).rank() if pts else 0
    print(f"   rank of {{doublets ∪ singlets}} as vectors: {R}  (rank 1 => all on a line through 0 => no Y separates)")

# 3) GATE 3: among b>=1 odd-nality, the singlet point-set — why >=3 needed
print("\n=== 3) GATE 3 geometry: singlet set at a+b+c=2 (fail) vs =3 (works) ===")
for L in [(0,1,1),(1,1,0),(0,2,1),(1,2,0),(0,1,3),(2,1,1)]:
    nb=(sum((i+1)*L[i] for i in range(3)))%2
    if nb==0: continue
    D,S=doublets_singlets(L)
    cell=cell_by_det(L)
    print(f"\n {L} a+b+c={sum(L)} dim={A3(*L).degree()}: #dbl={len(D)} #sng={len(S)} admit={'YES' if cell else 'no'}")
    print(f"   singlets: {S}")
    if cell: print(f"   CELL Q={cell[0]} u={cell[1]} d={cell[2]}")
