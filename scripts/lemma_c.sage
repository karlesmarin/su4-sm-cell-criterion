# Lemma (c): closed form for #RH-color-singlets(a,b,c); show (with b>=1, odd n-ality) that #>=3 <=> a+b+c>=3.
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2): rem[m2]=rem.get(m2,0)-m
        out.extend([top+1]*m)
    return out
def nsinglets(L):
    boxes=sum((i+1)*L[i] for i in range(3)); buckets={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        if (-1)**n4!=-1 or (-1)**(n3+n4)!=-1: continue
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        bb=buckets.setdefault((q8,q15),{}); bb[n1-n2]=bb.get(n1-n2,0)+mult
    S=set()
    for (q8,q15),bym in buckets.items():
        for d in su2_decompose(bym):
            if d==1: S.add((q8,q15))
    return len(S)

# tabulate n_singlets(a,b,c) for small labels, odd n-ality
print("=== #RH singlets vs (a,b,c) [odd n-ality] ===")
data=[]
for a in range(7):
 for b in range(1,6):
  for c in range(7):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>500: continue
   ns=nsinglets(L)
   data.append((L,sum(L),ns))
# guess: does ns depend only on some simple combo? print grouped
for L,s,ns in sorted(data,key=lambda x:(x[1],x[0])):
    print(f"  {L} a+b+c={s:>2}  #singlets={ns}")

print("\n=== VERIFY closed form  N_singlets = (b+1)*(a+c+1)/2  [odd n-ality] ===")
bad=0; n=0
for a in range(9):
 for b in range(9):
  for c in range(9):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>600: continue
   n+=1
   pred=(b+1)*(a+c+1)//2
   act=nsinglets(L)
   if pred!=act:
       bad+=1
       if bad<=15: print(f"  MISMATCH {L}: pred={pred} act={act}")
print(f"tested {n} odd-nality reps (dim<=600): mismatches={bad}")
if bad==0:
    print("*** N_singlets(a,b,c) = (b+1)*(a+c+1)/2  exactly (odd n-ality) ***")
    print("  Gate 3 derivation: with b>=1 and a+c>=1 odd,")
    print("    N=2  <=>  (b+1)(a+c+1)=4  <=>  b=1 & a+c=1  <=>  a+b+c=2   (the dim-20 minimal reps)")
    print("    N>=3 for ALL other (b>=1) reps  <=>  a+b+c>=3.   Gate 3 = 'not the b=1,a+c=1 corner'.")
    print("  And N>=3 (b>=1) <=> singlets span 2D <=> cell closes.  FULL THEOREM assembled.")
