# Doublet-cooperation step: show every admitting rep has a SYMMETRIC cell Q=(u+d)/2 with Q a LH doublet.
# If so, Y(Q)=(Y(u)+Y(d))/2=1/6 is automatic and the doublet cooperation is a clean lattice fact.
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2): rem[m2]=rem.get(m2,0)-m
        out.extend([top+1]*m)
    return out
def DS(L):
    boxes=sum((i+1)*L[i] for i in range(3)); Dc={}; Sc={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        p0=(-1)**n4; p2=(-1)**(n3+n4); q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        if (p0,p2)==(1,1):  Dc.setdefault((q8,q15),{})[n1-n2]=Dc.setdefault((q8,q15),{}).get(n1-n2,0)+mult
        if (p0,p2)==(-1,-1):Sc.setdefault((q8,q15),{})[n1-n2]=Sc.setdefault((q8,q15),{}).get(n1-n2,0)+mult
    D=set(); S=set()
    for k,bym in Dc.items():
        for d in su2_decompose(bym):
            if d==2: D.add(k)
    for k,bym in Sc.items():
        for d in su2_decompose(bym):
            if d==1: S.add(k)
    return D,S
def admits(L):
    D,S=DS(L)
    for Q in D:
        for u in S:
            if Matrix(QQ,[[Q[0],Q[1]],[u[0],u[1]]]).det()==0: continue
            for d in S:
                if d!=u and Matrix(QQ,[[Q[0],Q[1],QQ(1)/6],[u[0],u[1],QQ(2)/3],[d[0],d[1],QQ(-1)/3]]).det()==0:
                    return True
    return False
def has_symmetric_cell(L):
    D,S=DS(L); Sset=set(S)
    for u in S:
        for d in S:
            if u==d: continue
            mid=((u[0]+d[0]),(u[1]+d[1]))
            if mid[0]%2 or mid[1]%2: continue
            Q=(mid[0]//2,mid[1]//2)
            if Q in D:
                # also need u,d independent (Y well-defined) and Y(u)-Y(d)=1 solvable: it's automatic since 2x2 nonsing
                if Matrix(QQ,[[u[0],u[1]],[d[0],d[1]]]).det()!=0:
                    return (Q,u,d)
    return None

print("=== does EVERY admitting rep have a symmetric cell Q=(u+d)/2 ? (odd n-ality, dim<=900) ===")
nadm=0; nosym=[]
for a in range(16):
 for b in range(16):
  for c in range(16):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>900: continue
   if not admits(L): continue
   nadm+=1
   if has_symmetric_cell(L) is None: nosym.append((L,A3(*L).degree()))
print(f"admitting reps tested: {nadm}")
print(f"admitting reps WITHOUT a symmetric cell: {len(nosym)}")
if not nosym:
    print("*** every admitting rep has a SYMMETRIC cell Q=(u+d)/2, Q a LH doublet ***")
    print("    => Y(Q)=1/6 automatic; doublet cooperation reduces to: exists singlet pair whose midpoint is a doublet.")
else:
    for x in nosym[:20]: print("   no sym cell:",x)
