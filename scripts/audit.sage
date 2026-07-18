# Authoritative formula audit under Sage preparser (integer/rational arithmetic correct).
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")
def su2ct(counter):
    rem=dict(counter); out=[]
    while any(v>0 for v in rem.values()):
        t=max(k for k,v in rem.items() if v>0); m=rem[t]
        for x in range(t,-t-1,-2): rem[x]=rem.get(x,0)-m
        out.extend([t+1]*m)
    return out
def modes(L,par):
    boxes=sum((i+1)*L[i] for i in range(3)); bk={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        if ((-1)**n4,(-1)**(n3+n4))!=par: continue
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        bk.setdefault((q8,q15),{})[n1-n2]=bk.setdefault((q8,q15),{}).get(n1-n2,0)+mult
    return [(d,k[0],k[1]) for k,bym in bk.items() for d in su2ct(bym)]
def admits(L):
    D=list(set((m[1],m[2]) for m in modes(L,(1,1)) if m[0]==2))
    S=list(set((m[1],m[2]) for m in modes(L,(-1,-1)) if m[0]==1))
    for Q in D:
        for u in S:
            for d in S:
                if u==d: continue
                if Matrix(QQ,[[Q[0],Q[1],QQ(1)/6],[u[0],u[1],QQ(2)/3],[d[0],d[1],QQ(-1)/3]]).det()==0 \
                   and Matrix(QQ,[[Q[0],Q[1]],[u[0],u[1]]]).det()!=0: return True
    return False
cat=[]
for a in range(9):
 for b in range(9):
  for c in range(9):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0: continue
   dim=int(A3(*L).degree())
   if dim>400: continue
   if admits(L): cat.append((dim,L))
print("AUDIT1 catalog dims<=400:", sorted(d for d,_ in cat))
print("AUDIT1 min:", min(cat))
mism=0
for a in range(9):
 for b in range(9):
  for c in range(9):
   L=(a,b,c)
   if L==(0,0,0) or A3(*L).degree()>300: continue
   g=((a+2*b+3*c)%2==1) and b>=1 and (a+b+c)>=3
   if g!=admits(L): mism+=1
print("AUDIT2 criterion-vs-det mismatches (dim<=300):", mism)
badN=0
for a in range(9):
 for b in range(9):
  for c in range(9):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0 or A3(*L).degree()>400: continue
   if len(set((m[1],m[2]) for m in modes(L,(-1,-1)) if m[0]==1))!=(b+1)*(a+c+1)//2: badN+=1
print("AUDIT3 N=(b+1)(a+c+1)/2 mismatches (dim<=400):", badN)
# extent 12b check
badE=0
for a in range(7):
 for b in range(7):
  for c in range(7):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0 or A3(*L).degree()>300: continue
   ms=[2*m[1]+m[2] for m in modes(L,(-1,-1)) if m[0]==1]
   if ms and (max(ms)-min(ms))!=12*b: badE+=1
print("AUDIT4 extent=12b mismatches (dim<=300):", badE)
print("AUDIT5 Y check: Q,u,d =", [QQ(-7*q8+4*q15)/18 for q8,q15 in [(-1,-1),(0,3),(-2,-5)]])
