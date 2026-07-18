# C2: closed form for the node-independent invariant dim of SU(n). Conjecture (prospector): it's the
# SU(2)-trivial part of the Levi branching SU(n)->SU(2)xSU(n-2)xU(1)^2, an LR/Kostka sum.
# Test candidate: dim Inv_{sl2(a1)}(V_lambda) = sum over SU(n-2) content = a clean combinatorial number?
# For SU(5) (n-2=3): compare to explicit branching to SU(2)xSU(3)xU(1), take SU(2)-singlet dim.
# author: Carles Marin (Claude as AI assistant).
from collections import defaultdict
def node_inv(ct, lab, node=0):
    R=WeylCharacterRing(ct, style="coroots"); amb=R.space()
    cr=amb.simple_coroots()[amb.index_set()[node]]
    N=defaultdict(int)
    for w,mult in R(*lab).weight_multiplicities().items(): N[ZZ(w.scalar(cr))]+=mult
    return N[0]-N[2]
def branch_su2_singlet_dim(ct, lab):
    # branch to SU(2)xSU(n-2) via 'levi' deleting node 0 (removes alpha_1) -> A1 x A_{n-3}; sum dims of SU(2)-trivial pieces
    R=WeylCharacterRing(ct, style="coroots")
    n=R.cartan_type().rank()+1  # SU(n)
    try:
        S=WeylCharacterRing(f"A1xA{n-3}", style="coroots") if n-3>=1 else WeylCharacterRing("A1", style="coroots")
        br=R(*lab).branch(S, rule="levi")
    except Exception as e:
        return None
    # sum multiplicities*dim of components whose A1 part is trivial (highest weight 0 on first factor)
    tot=0
    for comp,mult in br.monomial_coefficients().items():
        # comp is a tuple of labels for A1 x A_{n-3}; A1 label = comp[0]
        lab_c=comp
        if lab_c[0]==0:  # SU(2)-trivial
            # dim of the remaining A_{n-3} irrep
            if n-3>=1:
                rest=WeylCharacterRing(f"A{n-3}", style="coroots"); tot+=mult*rest(*lab_c[1:]).degree()
            else:
                tot+=mult
    return tot
print("=== C2: SU(5) node-invariant values + test as SU(2)-singlet Levi-branching dim ===")
import itertools
R5=WeylCharacterRing("A4", style="coroots")
rows=[]
for lab in itertools.product(range(4),repeat=4):
    if sum(lab)==0: continue
    if R5(*lab).degree()>150: continue
    ni=node_inv("A4",lab)
    rows.append((R5(*lab).degree(),lab,ni))
for dim,lab,ni in sorted(rows)[:22]:
    print(f"  {lab} dim={dim:>3} node_inv={ni}")
# try candidate closed forms for SU(5): does node_inv match dim of SU(3) rep (lab[1],lab[2]) times something? print a few to eyeball
print("\n  eyeball vs SU(3)=(b,c) dim and simple products:")
r3=WeylCharacterRing("A2",style="coroots")
for lab in [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,0,0),(0,1,1,0),(2,0,0,0),(1,0,0,1),(0,2,0,0)]:
    ni=node_inv("A4",lab)
    print(f"    {lab}: node_inv={ni}")
