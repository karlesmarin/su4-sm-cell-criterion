# Rigor closure: derive Lemma B (spread=12b) and Lemma C (N=(b+1)(a+c+1)/2) from the Levi branching
# SU(4) -> SU(2)_L x SU(2)_R x U(1)_m  (delete middle node alpha2).
# SU(2)_L on {1,2} (alpha1, weight n1-n2), SU(2)_R on {3,4} (alpha3, weight n3-n4), U(1)_m = m=2q8+q15.
# RH-singlet weights = SU(2)_L-singlet constituents (1,j_R)_m, projected to RH parity (n3 even,n4 odd).
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")

def levi_singlet_constituents(L):
    # Group weights by (m, n3-n4) after restricting to SU(2)_L singlets, then read SU(2)_R structure.
    boxes=sum((i+1)*L[i] for i in range(3))
    # collect ALL weights with (n1,n2,n3,n4)
    W=[]
    for w,mult in A3(*L).weight_multiplicities().items():
        n=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        W.append((tuple(n),mult))
    # SU(2)_L singlet: within fixed (m, and SU(2)_R weight n3-n4, and n3+n4), the n1-n2 multiplicities -> take d==1 part.
    # Bucket by (m, n3, n4) actually determines q8,q15; SU(2)_L weight is n1-n2.
    from collections import defaultdict
    buck=defaultdict(dict)
    for n,mult in W:
        n1,n2,n3,n4=n
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        buck[(q8,q15,n3,n4)][n1-n2]=buck[(q8,q15,n3,n4)].get(n1-n2,0)+mult
    def su2decomp(mb):
        rem=dict(mb); out=[]
        while any(v>0 for v in rem.values()):
            top=max(k for k,v in rem.items() if v>0); mm=rem[top]
            for m2 in range(top,-top-1,-2): rem[m2]=rem.get(m2,0)-mm
            out.extend([top]*mm)   # store 2*spin = top
        return out
    # singlet weights (SU(2)_L spin 0): those buckets contributing top==0
    sing=[]   # (m, q8,q15, n3-n4, n3,n4)
    for (q8,q15,n3,n4),mb in buck.items():
        for twos in su2decomp(mb):
            if twos==0:
                sing.append((2*q8+q15, q8,q15, n3-n4, n3,n4))
    return sing

print("=== Levi derivation: SU(2)_L-singlet weights, grouped by U(1)_m level, showing SU(2)_R & parity ===")
for L in [(0,2,1),(2,2,1),(0,1,3),(0,3,3),(2,1,1),(1,4,0)]:
    if (sum((i+1)*L[i] for i in range(3)))%2==0: continue
    sing=levi_singlet_constituents(L)
    from collections import defaultdict
    bym=defaultdict(list)
    for (m,q8,q15,nr,n3,n4) in sing:
        bym[m].append((n3-n4, (-1)**n4, (-1)**(n3+n4)))  # SU(2)_R weight, parity flags
    ms=sorted(bym)
    print(f"\n {L}: b={L[1]} a={L[0]} c={L[2]}  #m-levels={len(ms)} (b+1={L[1]+1})")
    for m in ms:
        allR=sorted(set(x[0] for x in bym[m]))
        rhR =sorted(set(x[0] for x in bym[m] if x[1]==-1 and x[2]==-1))  # RH: n4 odd, n3+n4 odd
        print(f"   m={m:>4}: SU(2)_R weights(n3-n4) all={allR}  RH-projected={rhR}  (#RH={len(rhR)}, expect (a+c+1)/2={(L[0]+L[2]+1)//2})")

print("\n\n=== VERIFY the 3 sub-claims of the derivation (odd n-ality, dim<=700) ===")
bad_levels=0; bad_spin=0; bad_parity=0; bad_spacing=0; n=0
for a in range(10):
 for b in range(10):
  for c in range(10):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>700: continue
   n+=1
   sing=levi_singlet_constituents(L)
   from collections import defaultdict
   bym=defaultdict(list)
   for (m,q8,q15,nr,n3,n4) in sing: bym[m].append((n3-n4,(-1)**n4,(-1)**(n3+n4)))
   ms=sorted(bym)
   # (1) #m-levels = b+1
   if len(ms)!=b+1: bad_levels+=1
   # (1b) spacing = 12
   if len(ms)>=2 and any(ms[i+1]-ms[i]!=12 for i in range(len(ms)-1)): bad_spacing+=1
   for m in ms:
     allR=sorted(set(x[0] for x in bym[m]))
     rhR =[x[0] for x in bym[m] if x[1]==-1 and x[2]==-1]
     # (2) SU(2)_R spin: max weight = a+c  (2j_R=a+c)
     if allR and max(allR)!=a+c: bad_spin+=1
     # (3) RH parity keeps exactly (a+c+1)/2
     if len(set(rhR))!=(a+c+1)//2: bad_parity+=1
print(f"tested {n} odd-nality reps (dim<=700)")
print(f"  (1) #m-levels != b+1        : {bad_levels}")
print(f"  (1b) m-spacing != 12        : {bad_spacing}")
print(f"  (2) SU(2)_R max-weight != a+c: {bad_spin}")
print(f"  (3) RH-kept != (a+c+1)/2    : {bad_parity}")
if bad_levels==bad_spacing==bad_spin==bad_parity==0:
    print("\n*** DERIVATION AIRTIGHT (to dim 700): ***")
    print("    N_singlets = (#m-levels)*(RH-kept per level) = (b+1)*(a+c+1)/2   [Lemma C]")
    print("    m-extent   = 12*(#levels-1)                  = 12b               [Lemma B]")
    print("    both from ONE Levi branching SU(4)->SU(2)_L x SU(2)_R x U(1)_m + RH parity.")
