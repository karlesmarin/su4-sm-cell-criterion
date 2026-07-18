from sage.all import *
# Data-driven publication figures for Part II (Standard-Model cell criterion).
# author: Carles Marin (Claude as AI assistant).
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.size':11,
                     'axes.linewidth':0.8,'figure.dpi':150})
RB='#1F4E79'; GN='#2E8C36'; BR='#B5530F'; GY='#93a1b0'
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
    out=[]
    for k,bym in bk.items():
        for d in su2ct(bym): out.append((d,k[0],k[1]))
    return out

# ---------- FIG A: admissibility landscape ----------
smax, bmax = 11, 8
S=[s for s in range(1,smax+1) if s%2==1]
Bv=list(range(0,bmax+1))
Nmat=np.full((len(Bv),len(S)),np.nan)
adm=np.zeros((len(Bv),len(S)),dtype=bool)
for i,b in enumerate(Bv):
    for j,s in enumerate(S):
        Nmat[i,j]=(b+1)*(s+1)//2
        adm[i,j]= (b>=1) and (b+s>=3)
fig,ax=plt.subplots(figsize=(6.6,4.3))
cmap=LinearSegmentedColormap.from_list('cell',['#eef4fb','#9ec6e8','#3d7cc0','#1F4E79'])
im=ax.imshow(Nmat,origin='lower',cmap=cmap,aspect='auto',extent=[S[0]-1,S[-1]+1,-0.5,bmax+0.5])
cb=fig.colorbar(im,ax=ax,pad=0.02); cb.set_label(r'zero-mode count $N=(b{+}1)\frac{a+c+1}{2}$')
for i,b in enumerate(Bv):
    for j,s in enumerate(S):
        if adm[i,j]:
            ax.scatter(s,b,s=118,facecolors='none',edgecolors=GN,linewidths=1.8,zorder=3)
        else:
            ax.scatter(s,b,marker='x',s=34,color=BR,alpha=0.85,zorder=3)
ax.axhline(0.5,color=BR,lw=1.6,ls='--')
ax.text(smax-0.2,0.10,r'$b\geq1$  (middle node excited)',color=BR,ha='right',fontsize=9)
ax.text(1.05,2.5,r'$a{+}b{+}c\geq3$',color=BR,fontsize=9)
ax.annotate(r'$(0,2,1)=\mathbf{60}$  (Part I minimum)',xy=(1,2),xytext=(3.1,3.6),fontsize=9,color=GN,
            arrowprops=dict(arrowstyle='->',color=GN,lw=1.1))
ax.set_xlabel(r'$a+c$   (odd $\Leftrightarrow$ centre charge odd)'); ax.set_ylabel(r'middle label $b$')
ax.set_xticks(S); ax.set_yticks(Bv)
ax.set_title(r'The admissibility landscape: green rings pass all three gates',fontsize=10.5)
fig.tight_layout(); fig.savefig('fig_cell_landscape.pdf'); plt.close(fig)
print("fig_cell_landscape.pdf done")

# ---------- FIG B: dissolution with rank ----------
frac=[0.06,0.73,0.85]
cols=[BR,GN,RB]
status=['rigid gate\n(odd centre only)','dissolved\n(all centres)','dissolved\n(all centres)']
fig,ax=plt.subplots(figsize=(6.4,4.4))
ax.bar(range(3),[f*100 for f in frac],color=cols,width=0.60,zorder=3,edgecolor='white')
for i,f in enumerate(frac):
    ax.text(i,f*100+2.6,str(int(round(f*100)))+'%',ha='center',fontsize=13,color=cols[i],fontweight='bold')
    if f>0.3:
        ax.text(i,f*100-12,status[i],ha='center',fontsize=8.8,color='white')
    else:
        ax.text(i,f*100+12,status[i],ha='center',fontsize=8.8,color=cols[i])
ax.set_xticks(range(3))
ax.set_xticklabels(['$SU(4)$\nrank 3','$SU(5)$\nrank 4','$SU(6)$\nrank 5'],fontsize=11)
ax.set_ylabel('representations admitting the SM cell  (%)')
ax.set_ylim(0,100); ax.grid(axis='y',alpha=0.22,zorder=0)
ax.axvline(0.5,color='0.6',lw=1.0,ls=':')
ax.text(0.5,48,'critical rank',rotation=90,va='center',ha='center',fontsize=8.8,color='0.45',
        bbox=dict(fc='white',ec='none',pad=1.5))
ax.text(-0.44,98,'free hypercharge:  (rank$-$1) parameters vs 3 constraints',fontsize=8.8,color='0.32',va='top')
ax.set_title('The rigid cell gate dissolves above the critical rank',fontsize=11)
fig.tight_layout(); fig.savefig('fig_cell_dissolution.pdf'); plt.close(fig)
print("fig_cell_dissolution.pdf done")

# ---------- FIG C: real sheared singlet lattice for (0,3,3) ----------
L=(0,3,3)
Sng=sorted(set((int(m[1]),int(m[2])) for m in modes(L,(-1,-1)) if m[0]==1))
Dbl=sorted(set((int(m[1]),int(m[2])) for m in modes(L,(1,1)) if m[0]==2))
Dset=set(Dbl); cell=None
for u in Sng:
    for d in Sng:
        if u==d: continue
        mid=((u[0]+d[0]),(u[1]+d[1]))
        if mid[0]%2 or mid[1]%2: continue
        Q=(mid[0]//2,mid[1]//2)
        if Q in Dset and Matrix(QQ,[[u[0],u[1]],[d[0],d[1]]]).det()!=0:
            cell=(Q,u,d); break
    if cell: break
from collections import defaultdict
mvals=sorted(set(2*p[0]+p[1] for p in Sng))
cmapL=plt.get_cmap('viridis')
fig,ax=plt.subplots(figsize=(6.4,4.6))
bym=defaultdict(list)
for p in Sng: bym[2*p[0]+p[1]].append(p)
for m,ps in bym.items():
    ps=sorted(ps); c=cmapL(float(mvals.index(m))/max(1,len(mvals)-1))
    ax.plot([q[0] for q in ps],[q[1] for q in ps],color=c,lw=1.1,alpha=0.7,zorder=2)
for p in Dbl: ax.scatter(p[0],p[1],marker='s',s=26,color=GY,alpha=0.5,zorder=1)
for p in Sng:
    m=2*p[0]+p[1]; c=cmapL(float(mvals.index(m))/max(1,len(mvals)-1))
    ax.scatter(p[0],p[1],s=95,color=c,edgecolor='white',linewidth=0.7,zorder=3)
if cell:
    Q,u,d=cell
    ax.plot([u[0],d[0]],[u[1],d[1]],color=BR,lw=1.9,ls='--',zorder=4)
    ax.scatter(*Q,marker='*',s=360,color=BR,edgecolor='white',linewidth=0.8,zorder=5)
    for pt,lab,dx,dy in [(Q,r'$Q\,(\frac{1}{6})$',0.35,0.6),(u,r'$u\,(\frac{2}{3})$',0.3,-1.2),(d,r'$d\,(-\frac{1}{3})$',0.35,0.6)]:
        ax.annotate(lab,xy=pt,xytext=(pt[0]+dx,pt[1]+dy),fontsize=10.5,color=BR)
ax.set_xlabel(r'$q_8$'); ax.set_ylabel(r'$q_{15}$')
ax.set_title(r'Real singlet lattice of $(0,3,3)$: four sheared middle-node levels',fontsize=10.3)
sm=plt.cm.ScalarMappable(cmap=cmapL); sm.set_array([0,1])
cb=fig.colorbar(sm,ax=ax,pad=0.02); cb.set_label(r'middle-node charge $m=2q_8+q_{15}$')
cb.set_ticks(list(np.linspace(0,1,len(mvals)))); cb.set_ticklabels([str(v) for v in mvals])
ax.grid(alpha=0.18); fig.tight_layout(); fig.savefig('fig_cell_lattice.pdf'); plt.close(fig)
print("fig_cell_lattice.pdf done; cell =",cell)
