# CORRECT F1 derivation. Levi SU(4)->SU(2)_L x SU(2)_R x U(1)_m (delete middle node).
# Coordinates: w_L=n1-n2, w_R=n3-n4, m=3(n1+n2-n3-n4)=2q8+q15.  (w_L,w_R,m)+boxes determine (n1..n4).
# SU(2)_L-singlets: bucket by (w_R,m), decompose w_L-multiplicities, take spin-0 count.
# CLAIM (cartographer): at each m-level the SU(2)_R content of the L-singlets = CG tower [a/2](x)[c/2].
#   #m-levels=b+1; distinct (q8,q15) after RH parity=(a+c+1)/2 (extent=top spin); N=(b+1)(a+c+1)/2.
# author: Carles Marin (Claude as AI assistant).
A3 = WeylCharacterRing("A3", style="coroots")
def su2_multiset(counter):   # counter: w->mult ; return sorted list of 2j (dim-1) present
    rem=dict(counter); out=[]
    while any(v>0 for v in rem.values()):
        t=max(k for k,v in rem.items() if v>0); m=rem[t]
        for x in range(t,-t-1,-2): rem[x]=rem.get(x,0)-m
        out.extend([t]*m)
    return sorted(out)
def cg_tower(a,c):
    return sorted(range(abs(a-c),a+c+1,2))   # 2j values of [a/2](x)[c/2]
def analyze(L):
    a,b,c=L; boxes=a+2*b+3*c
    from collections import defaultdict
    buck=defaultdict(lambda: defaultdict(int))      # (w_R,m) -> {w_L: mult}
    for w,mult in A3(*L).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        wL=n1-n2; wR=n3-n4; m=3*(n1+n2-n3-n4)
        buck[(wR,m)][wL]+=mult
    # singlet (spin-0) count per (w_R,m)
    sing=defaultdict(int)
    for (wR,m),mb in buck.items():
        twos=su2_multiset(mb)
        sing[(wR,m)] += twos.count(0)
    # per m-level: SU(2)_R content = decompose w_R multiplicities (each wR with weight = sing count)
    bym=defaultdict(lambda: defaultdict(int))
    rh_distinct=defaultdict(set)
    for (wR,m),k in sing.items():
        if k>0:
            bym[m][wR]+=k
            # parity: recover n3,n4. t=n3+n4=(boxes - m/3)/2 ; n3=(t+wR)/2,n4=(t-wR)/2
            t=(boxes-m/3)/2
            if t in ZZ and (t+wR)%2==0:
                n3=(t+wR)//2; n4=(t-wR)//2
                if (-1)**n4==-1 and (-1)**(n3+n4)==-1:   # RH sector
                    q8=n1+n2-2*n3 if False else (m/3 - wR)  # q8=m/3 - wR
                    q15=m/3 + 2*wR
                    rh_distinct[m].add((q8,q15))
    tower=cg_tower(a,c)
    tower_ok=all(su2_multiset(bym[m])==tower for m in bym)
    nlevels=len(bym)
    rh_ok=all(len(rh_distinct[m])==(a+c+1)//2 for m in bym if m in rh_distinct)
    Npred=(b+1)*(a+c+1)//2
    Nrh=sum(len(rh_distinct[m]) for m in rh_distinct)
    return tower,tower_ok,nlevels,(b+1),rh_ok,Nrh,Npred

print("=== F1 BLINDAJE: L-singlet SU(2)_R = CG tower? #levels=b+1? N=(b+1)(a+c+1)/2? ===")
bad=0; tested=0
import itertools
for a in range(9):
 for b in range(9):
  for c in range(9):
   L=(a,b,c)
   if (a+2*b+3*c)%2==0: continue
   if A3(*L).degree()>600: continue
   tested+=1
   tower,tok,nl,b1,rhok,Nrh,Npred=analyze(L)
   if not(tok and nl==b1 and rhok and Nrh==Npred):
       bad+=1
       if bad<=12: print(f"  FAIL {L}: tower_ok={tok} levels={nl}/{b1} rh_ok={rhok} N={Nrh}/{Npred}")
print(f"tested {tested} odd-nality reps (dim<=600): failures={bad}")
if bad==0:
    print("*** F1 DERIVATION AIRTIGHT: L-singlet sector = (b+1) copies of CG tower [a/2](x)[c/2];")
    print("    the tower's EXTENT (top spin (a+c)/2) sets distinct (q8,q15); RH parity halves -> (a+c+1)/2/level;")
    print("    N=(b+1)(a+c+1)/2. Prose F3 corrected: it is the TOWER, extent gives the count.")
# show one explicit tower vs single-spin contrast
print("\n  explicit (2,2,1): CG tower 2j =", cg_tower(2,1), " (i.e. [3/2](+)[1/2], NOT a single spin 3/2)")
