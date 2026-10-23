"""Reproduce six main figures and Figure S11 from frozen numerical inputs."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR","/tmp/aim1-freeze-mpl")
import runpy,json,hashlib,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
H=Path(__file__).resolve().parent
(H/'supporting_figures').mkdir(exist_ok=True)
runpy.run_path(str(H/'reproduce_legacy_plots.py'),run_name='__main__')
F=H/'reproduced_figures';T=H/'reproduced_tables'
for ext in ['png','pdf','svg']:
 shutil.copy2(F/('Figure_4_conditional_PMF.'+ext),H/'supporting_figures'/('S11.'+ext))
 shutil.copy2(F/('Figure_6_structural_response.'+ext),F/('Figure_4_structural_response.'+ext))
# Assemble existing molecular renders with exact quantitative tables; no synthetic structures.
fig=plt.figure(figsize=(7.1,5.7));gs=fig.add_gridspec(2,3,left=.08,right=.93,bottom=.12,top=.90,height_ratios=[1.3,1],hspace=.62,wspace=.32)
for j,(stem,label) in enumerate([('patch_a_nx-c_medoid','Family A (NX-C)'),('patch_b_py-c_medoid','Family B (PY-C)'),('patch_c_pz-c_medoid','Family C (PZ-C)')]):
 ax=fig.add_subplot(gs[0,j]);im=plt.imread(H/'molecular_inputs'/(stem+'.png'));rgb=im[:,:,:3];mask=np.any(rgb<.985,axis=2);rows,cols=np.where(mask);im=im[max(0,rows.min()-20):min(im.shape[0],rows.max()+21),max(0,cols.min()-20):min(im.shape[1],cols.max()+21)];ax.imshow(im);ax.axis('off');ax.set_title(label,fontsize=10,pad=8)
fig.text(.08,.96,'A  Several protein interfaces on the same polymer model',fontsize=11,fontweight='bold')
fig.text(.08,.505,'Gold: polystyrene   Gray: ubiquitin/contact participants   Red: recurrent contact residues',fontsize=8)
rel=pd.read_csv(T/'completed_release_summary.csv');faces=['px','nx','py','ny','pz','nz'];key='starting_face' if 'starting_face' in rel else 'face';counts=rel[rel.outcome=='retained'].groupby(key).size().reindex(faces).fillna(0)
a=fig.add_subplot(gs[1,0]);y=np.arange(6);a.barh(y,[3]*6,facecolor='white',edgecolor='.55',height=.62);a.barh(y,counts,color='.25',height=.62);a.set(yticks=y,yticklabels=['+X','−X','+Y','−Y','+Z','−Z'],xticks=[0,1,2,3],xlim=(0,3.25),xlabel='Retained trials / 3');a.invert_yaxis();a.set_title('B  Retention from every face',loc='left',fontsize=9,fontweight='bold',pad=12)
b=fig.add_subplot(gs[1,1:]);atlas=pd.read_csv(T/'F1_family_mean_occupancy.csv').set_index('patch_family');im=b.imshow(atlas.values,aspect='auto',cmap='Greys',vmin=0,vmax=1,extent=[.5,76.5,2.5,-.5]);b.set(yticks=[0,1,2],yticklabels=['A (7)','B (3)','C (3)'],xticks=[1,20,40,60,76],xlabel='Ubiquitin residue');b.set_title('C  Distinct residue-contact patterns',loc='left',fontsize=9,fontweight='bold',pad=12);cb=fig.colorbar(im,ax=b,fraction=.035,pad=.03);cb.set_ticks([0,1]);cb.set_label('Contact occupancy',fontsize=8)
fig.text(.08,.025,'18 prepared encounters → 13 retained releases → 3 descriptive families; no equilibrium populations.',fontsize=8)
def save(fig,name):
 fig.set_size_inches(7.0,fig.get_size_inches()[1]*7.0/fig.get_size_inches()[0],forward=True)
 for ext in ['pdf','svg','png']:fig.savefig(F/(name+'.'+ext),dpi=600)
 fig.savefig(H/'review'/(name+'_final_size.png'),dpi=180);plt.close(fig)
save(fig,'Figure_1_interfaces')
# Evidence synthesis: arrows connect observations to interpretation, not a reaction trajectory.
fig,ax=plt.subplots(figsize=(7.1,3.8));fig.subplots_adjust(left=.02,right=.98,bottom=.04,top=.96);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax.text(.5,.96,'An irregular polystyrene surface accommodates multiple interfaces',ha='center',va='top',fontsize=11,fontweight='bold')
items=[(.02,'Interface diversity','Different protein faces\nretain distinct contacts','Orientation releases\n+ exposure-normalized chemistry'),(.35,'Direct attraction','Larger buried area\naccompanies stronger LJ','Matched-distance energies\n+ atom-group localization'),(.68,'Retained hydration','Blocked space loses water;\naccessible space stays hydrated','Spatial masks + hydration\nenergy and entropy controls')]
for x,title,obs,evidence in items:
 ax.add_patch(FancyBboxPatch((x,.43),.29,.37,boxstyle='round,pad=.009',facecolor='#f4f7f8',edgecolor='#647681',linewidth=.8));ax.text(x+.145,.75,title,ha='center',fontsize=10,fontweight='bold');ax.text(x+.145,.65,obs,ha='center',va='center',fontsize=9,linespacing=1.35);ax.text(x+.145,.49,evidence,ha='center',va='center',fontsize=8,color='.3',linespacing=1.3);ax.annotate('',xy=(.5,.31),xytext=(x+.145,.42),arrowprops=dict(arrowstyle='->',lw=.9,color='#647681'))
ax.add_patch(FancyBboxPatch((.16,.17),.68,.14,boxstyle='round,pad=.009',facecolor='#e8f1f5',edgecolor='#0072b2',linewidth=1));ax.text(.5,.24,'Chemically distinct, attractive, partially hydrated contacts\nwith much of the sampled ubiquitin core retained',ha='center',va='center',fontsize=10,linespacing=1.3)
ax.text(.5,.065,'Evidence synthesis, not equilibrium state weights or a free-energy decomposition.',ha='center',fontsize=8)
save(fig,'Figure_6_mechanism')
# Keep output directory unambiguous: only the six current main figures.
for stem in ['Figure_4_conditional_PMF','Figure_6_structural_response']:
 for ext in ['png','pdf','svg']:(F/(stem+'.'+ext)).unlink()
manifest={str(p.relative_to(H)):hashlib.sha256(p.read_bytes()).hexdigest() for d in [F,H/'supporting_figures'] for p in d.iterdir() if p.is_file()}
(H/'REPRODUCED_FIGURES.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Six main figures and supporting Figure S11 have been rebuilt.')
