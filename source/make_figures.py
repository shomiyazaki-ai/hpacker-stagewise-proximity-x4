from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator,LogFormatterMathtext
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'code'))
from verify_e1_from_pdb import parse_pdb
D=R/'data/authoritative';F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'pdf.fonttype':42})
e1=json.loads((R/'data/expected_e1.json').read_text());cols=['#176b9e','#ce7624']
fig,axes=plt.subplots(1,2,figsize=(7.0,2.85),layout='constrained')
for ax,vals,col,label in zip(axes,[e1['general_stage_counts'],e1['phe_stage_counts']],cols,['A  General proximity','B  PHE6-PHE55 proximity']):
 ax.bar(range(6),vals,color=col,width=.65);ax.set(ylim=(0,200),yticks=[0,50,100,150,200],xticks=range(6),xlabel='HPacker stage',ylabel='Source frames / 200',title=label)
 for x,y in enumerate(vals):ax.text(x,y+4,str(y),ha='center',fontsize=9)
 ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
fig.savefig(F/'Figure_1.png');fig.savefig(F/'Figure_1.pdf');plt.close(fig)
# Six 3D views of the exact stored PHE side chains; fixed limits, camera, no alignment per stage.
data=[]
for s in range(6):
 k,x=parse_pdb(R/'data/pdb/WT_source000'/f'stage_{s:02d}.pdb');data.append({(kk[1],kk[-1]):v for kk,v in zip(k,x) if kk[1] in [6,55]})
show={'CA','CB','CG','CD1','CD2','CE1','CE2','CZ'}
coords=np.array([v for d in data for (r,a),v in d.items() if a in show]);center=(coords.max(0)+coords.min(0))/2;span=max(np.ptp(coords,axis=0))*.59
bonds=[('CA','CB'),('CB','CG'),('CG','CD1'),('CG','CD2'),('CD1','CE1'),('CD2','CE2'),('CE1','CZ'),('CE2','CZ')]
f0=list(csv.DictReader((D/'Figure_2_Source_0.csv').open()))
fig=plt.figure(figsize=(7.2,4.85))
for s,d in enumerate(data):
 ax=fig.add_subplot(2,3,s+1,projection='3d');ax.set_proj_type('ortho');ax.view_init(elev=22,azim=-63);ax.set_box_aspect((1,1,1))
 for r,c in [(6,cols[0]),(55,cols[1])]:
  for a,b in bonds:
   pts=np.array([d[r,a],d[r,b]]);ax.plot(*pts.T,color=c,lw=2.1)
  pts=np.array([v for (rr,aa),v in d.items() if rr==r and aa in show]);ax.scatter(*pts.T,s=13,color=c,depthshade=False)
 row=f0[s];p=np.array([d[6,row['phe655_atom6']],d[55,row['phe655_atom55']]]);ax.plot(*p.T,ls='--',color='#444444',lw=1.0)
 ax.set(xlim=(center[0]-span,center[0]+span),ylim=(center[1]-span,center[1]+span),zlim=(center[2]-span,center[2]+span));ax.set_axis_off()
 ax.set_title(f"Stage {s}  |  {float(row['phe655_min_distance_A']):.3f} Å",fontsize=10,pad=-7)
 ax.text2D(.5,-.02,row['phe655_atom6']+' - '+row['phe655_atom55'],transform=ax.transAxes,ha='center',fontsize=8,color='#444444')
fig.text(.5,.015,'PHE6 (blue)   /   PHE55 (orange)   |   Identical camera and coordinate limits',ha='center',fontsize=9)
fig.subplots_adjust(left=.01,right=.99,top=.94,bottom=.085,wspace=.0,hspace=.15);fig.savefig(F/'Figure_2.png');fig.savefig(F/'Figure_2.pdf');plt.close(fig)
# Exact minimum-distance trace (not projected lengths).
ys=[float(x['phe655_min_distance_A']) for x in f0];fig,ax=plt.subplots(figsize=(6.5,2.8),layout='constrained');ax.plot(range(6),ys,'o-',color=cols[0]);ax.axhline(1.6,color='#bb493e',ls='--',lw=1,label='1.6 Å threshold')
for i,y in enumerate(ys):ax.annotate(f'{y:.3f}',(i,y),xytext=(0,9),textcoords='offset points',ha='center',fontsize=9)
ax.set(xticks=range(6),xlabel='HPacker stage',ylabel='PHE6-PHE55 minimum\ndistance (Å)',ylim=(0,5.1));ax.legend(frameon=False,loc='upper right');ax.grid(axis='y',alpha=.17);fig.savefig(F/'Figure_S1.png');fig.savefig(F/'Figure_S1.pdf');plt.close(fig)
lj=list(csv.DictReader((R/'verification/current/e4_lj.csv').open()));fig,ax=plt.subplots(figsize=(6.5,3.1),layout='constrained')
for s,c in [(0,cols[0]),(1,cols[1])]:
 y=[float(r['archived_summary_kj_mol']) for r in lj if int(r['source'])==s];xx=np.array([0,1])+(-.05 if s==0 else .05);ax.plot(xx,y,'o-',color=c,label=f'Source {s}',ms=6)
 for x,v in zip(xx,y):ax.annotate(f'{v:.3e}',(x,v),xytext=((-10 if s==0 else 12),(12 if s==0 else -16) if x<.5 else (24 if s==0 else 20)),textcoords='offset points',ha='right' if s==0 else 'left',fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.3))
ax.set_yscale('log');ax.set(xticks=[0,1],xticklabels=['Original','Locally perturbed'],xlim=(-.55,1.62),ylim=(1e4,1e14),ylabel='Net heavy-atom LJ proxy\n(kJ/mol; logarithmic axis)');ax.yaxis.set_major_locator(LogLocator(base=10));ax.yaxis.set_major_formatter(LogFormatterMathtext());ax.grid(axis='y',alpha=.18);ax.legend(frameon=False,loc='upper right');fig.savefig(F/'Figure_S2.png');fig.savefig(F/'Figure_S2.pdf');plt.close(fig)
print('Rendered four data-based figures (PNG + vector PDF)')
