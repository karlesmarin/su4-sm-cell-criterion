# Test: does admit (given odd n-ality) <=> RH color-singlet set contains 3 AFFINELY-INDEPENDENT points?
# i.e. gates 2 AND 3 collapse to ONE condition: "the singlets are 2-dimensional (not collinear)".
# author: Carles Marin (Claude as AI assistant).
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
        bb=buckets.setdefault((q8,q15),{}); bb[n1-n2]=bb.get(n1-n2,0)+mult
    out=[]
    for (q8,q15),bym in buckets.items():
        for d in su2_decompose(bym): out.append((d,q8,q15))
    return out
def DS(L):
    D=sorted(set((m[1],m[2]) for m in modes_of(L,(1,1)) if m[0]==2))
    S=sorted(set((m[1],m[2]) for m in modes_of(L,(-1,-1)) if m[0]==1))
    return D,S
def admits(L):
    D,S=DS(L)
    for Q in D:
        for u in S:
            if Matrix(QQ,[[Q[0],Q[1]],[u[0],u[1]]]).det()==0: continue
            for d in S:
                if d==u: continue
                if Matrix(QQ,[[Q[0],Q[1],QQ(1)/6],[u[0],u[1],QQ(2)/3],[d[0],d[1],QQ(-1)/3]]).det()==0:
                    return True
    return False
def singlets_2D(S):   # exists 3 affinely-independent singlets?
    n=len(S)
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                v1=vector(QQ,S[j])-vector(QQ,S[i]); v2=vector(QQ,S[k])-vector(QQ,S[i])
                if Matrix(QQ,[v1,v2]).det()!=0: return True
    return False

print("=== admit  vs  (RH singlets contain 3 affinely-independent pts) — odd n-ality, dim<=300 ===")
mism=[]; n=0; both=0
for a in range(12):
 for b in range(12):
  for c in range(12):
   L=(a,b,c)
   if L==(0,0,0): continue
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>300: continue
   n+=1
   A=admits(L); G=singlets_2D(DS(L)[1])
   if A: both+=1
   if A!=G: mism.append((L,sum(L),A,G))
print(f"tested {n} odd-nality reps (dim<=300), admit={both}")
if not mism:
    print("*** THEOREM (empirical): among odd n-ality, admit <=> RH color-singlets span 2D (3 non-collinear) ***")
    print("    => gates 2 (b>=1) and 3 (a+b+c>=3) UNIFY into one geometric fact: the singlet diagram is 2-dimensional.")
else:
    print("counterexamples:")
    for m in mism[:20]: print("   ",m)

# Now: WHY is 'singlets 2D' equivalent to the two arithmetic gates? Show b=0 => singlets collinear; a+b+c=2 => <=2 singlets.
print("\n=== structural: singlet affine-rank as a function of (b, a+b+c) ===")
def affrank(S):
    if len(S)<=1: return len(S)-1 if S else -1
    base=vector(QQ,S[0]); M=Matrix(QQ,[vector(QQ,s)-base for s in S[1:]])
    return M.rank()
import collections
tab=collections.defaultdict(list)
for a in range(10):
 for b in range(10):
  for c in range(10):
   L=(a,b,c)
   if L==(0,0,0): continue
   if (a+2*b+3*c)%2: 
     if A3(*L).degree()>300: continue
     S=DS(L)[1]
     tab[(min(b,1),'b0' if b==0 else 'b>=1')].append((sum(L),affrank(S),len(S)))
for key in sorted(tab):
    vals=tab[key]; ranks=sorted(set(r for _,r,_ in vals))
    print(f"  {key[1]:>5}: singlet affine-ranks observed = {ranks}")
