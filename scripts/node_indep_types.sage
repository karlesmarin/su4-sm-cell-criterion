# Is node-independence (dim of single-simple-root SU(2)-invariants = same for all nodes) A-type-special,
# or does it hold for B,C,D,G2 (non-simply-laced included, where the node-SU(2)s have different root lengths)?
# author: Carles Marin (Claude as AI assistant).
from collections import defaultdict
def test_type(ct, dimcap, maxlabel):
    R = WeylCharacterRing(ct, style="coroots")
    amb = R.space()
    coroots = [amb.simple_coroots()[i] for i in amb.index_set()]
    rank = R.cartan_type().rank()
    import itertools
    bad=0; tested=0; examples=[]
    for lab in itertools.product(range(maxlabel), repeat=rank):
        if sum(lab)==0: continue
        rep=R(*lab)
        if rep.degree()>dimcap: continue
        tested+=1
        wm=rep.weight_multiplicities()
        vals=[]
        for cr in coroots:
            N=defaultdict(int)
            for w,mult in wm.items():
                N[ZZ(w.scalar(cr))]+=mult
            vals.append(N[0]-N[2])
        if len(set(vals))!=1:
            bad+=1
            if len(examples)<4: examples.append((lab,vals))
    print(f"  {ct}: tested {tested} irreps (dim<={dimcap}) | node-indep FAILS = {bad}")
    for e in examples: print(f"       {e[0]}: node-values={e[1]}")
    return bad

print("=== Node-independence across Lie types (does dim of single-node SU(2)-invariants depend on node?) ===")
for ct in ["A2","A3","A4","B2","B3","C3","D4","G2"]:
    test_type(ct, 200, 5)
