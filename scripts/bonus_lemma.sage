# Bonus lemma: dim of single-simple-root-Levi-SU(2) invariants of SU(n) irrep lambda = prod_j(lambda_j+1),
# node-independent. Verify for SU(5)=A4 and SU(6)=A5.
# For simple root i, the SU(2)-weight of an ambient weight w is w[i]-w[i+1]; #singlets = N(0)-N(2).
# author: Carles Marin (Claude as AI assistant).
def levi_singlets(ring, n, L, node):
    from collections import defaultdict
    N=defaultdict(int)
    for w,mult in ring(*L).weight_multiplicities().items():
        v = w[node]-w[node+1]
        N[v]+=mult
    return N[0]-N[2]     # multiplicity of spin-0 in the node-i SU(2)

def test(n, dimcap, maxlabel):
    ring = WeylCharacterRing(f"A{n-1}", style="coroots")
    bad_prod=0; bad_nodeindep=0; tested=0
    import itertools
    for L in itertools.product(range(maxlabel), repeat=n-1):
        if sum(L)==0: continue
        if ring(*L).degree()>dimcap: continue
        tested+=1
        prod=1
        for x in L: prod*=(x+1)
        vals=[levi_singlets(ring,n,L,i) for i in range(n-1)]
        if any(v!=prod for v in vals): 
            bad_prod+=1
            if bad_prod<=6: print(f"   SU({n}) {L}: nodes={vals} prod={prod}")
        if len(set(vals))!=1: bad_nodeindep+=1
    print(f"SU({n}): tested {tested} irreps (dim<={dimcap}) | prod-formula fails={bad_prod} | node-indep fails={bad_nodeindep}")
    return bad_prod,bad_nodeindep

print("=== BONUS LEMMA: dim(single-node Levi-SU(2) invariants) = prod_j(lambda_j+1), node-independent ===\n")
r4=test(4, 400, 7)
r5=test(5, 400, 6)
r6=test(6, 600, 5)
print()
if r4==(0,0) and r5==(0,0) and r6==(0,0):
    print("*** LEMMA HOLDS for SU(4),SU(5),SU(6): dim of ANY single-simple-root Levi-SU(2) invariants")
    print("    of SU(n) irrep (lambda_1..lambda_{n-1}) = prod_j (lambda_j+1), independent of which node. ***")
    print("    => N=(b+1)(a+c+1)/2 for SU(4) is a corollary: full L-singlet space=(a+1)(b+1)(c+1);")
    print("       the SM-relevant count is the U(1)-graded EXTENT halved, but the invariant-DIM is this product.")
