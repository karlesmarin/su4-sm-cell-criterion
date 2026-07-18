# Harden the theorem to higher dim: admit <=> (odd nality & b>=1 & a+b+c>=3) <=> singlets span 2D.
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
        p0=(-1)**n4; p2=(-1)**(n3+n4)
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        if (p0,p2)==(1,1):  Dc.setdefault((q8,q15),{})[n1-n2]=Dc.setdefault((q8,q15),{}).get(n1-n2,0)+mult
        if (p0,p2)==(-1,-1):Sc.setdefault((q8,q15),{})[n1-n2]=Sc.setdefault((q8,q15),{}).get(n1-n2,0)+mult
    D=set(); S=set()
    for k,bym in Dc.items():
        for d in su2_decompose(bym):
            if d==2: D.add(k)
    for k,bym in Sc.items():
        for d in su2_decompose(bym):
            if d==1: S.add(k)
    return sorted(D),sorted(S)
def admits(L):
    D,S=DS(L)
    for Q in D:
        for u in S:
            if Matrix(QQ,[[Q[0],Q[1]],[u[0],u[1]]]).det()==0: continue
            for d in S:
                if d!=u and Matrix(QQ,[[Q[0],Q[1],QQ(1)/6],[u[0],u[1],QQ(2)/3],[d[0],d[1],QQ(-1)/3]]).det()==0:
                    return True
    return False
def singlets_2D(S):
    S=list(S); n=len(S)
    for i in range(n):
     for j in range(i+1,n):
      for k in range(j+1,n):
        if Matrix(QQ,[vector(QQ,S[j])-vector(QQ,S[i]),vector(QQ,S[k])-vector(QQ,S[i])]).det()!=0: return True
    return False

bad1=0;bad2=0;n=0;adm=0
for a in range(16):
 for b in range(16):
  for c in range(16):
   L=(a,b,c)
   if L==(0,0,0): continue
   if A3(*L).degree()>900: continue
   n+=1
   gates = ((a+2*b+3*c)%2==1) and (b>=1) and (a+b+c>=3)
   A=admits(L)
   if A: adm+=1
   if A!=gates: bad1+=1; print(f"  GATE MISMATCH {L} dim={A3(*L).degree()} admit={A} gates={gates}")
   if (a+2*b+3*c)%2==1:
       G2=singlets_2D(DS(L)[1])
       if A!=G2: bad2+=1; print(f"  2D MISMATCH {L} admit={A} singlets2D={G2}")
print(f"\ntested {n} reps (dim<=900), admit={adm}")
print(f"admit<=>(odd&b>=1&a+b+c>=3): mismatches={bad1}")
print(f"admit<=>singlets-2D (odd nality): mismatches={bad2}")
if bad1==0 and bad2==0: print("*** THEOREM holds to dim 900: gates AND geometric form both exact ***")
