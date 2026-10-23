from pathlib import Path
import os,json,hashlib,shutil
os.environ['MPLCONFIGDIR']='/tmp/aim1-polish-mpl'
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle,Ellipse
H=Path(__file__).resolve().parent;F=H/'reproduced_figures';T=H/'reproduced_tables';pins={};exports={};facts={}
for directory in [F,T,H/'review']:directory.mkdir(exist_ok=True)
INPUT_MAP=json.loads((H/'INPUT_MAP.json').read_text())
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Nimbus Roman','Liberation Serif','DejaVu Serif'],'mathtext.fontset':'stix','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'lines.linewidth':1.2,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
C={'A':'#0072B2','B':'#D55E00','C':'#009E73'};markers={'A':'o','B':'s','C':'^'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path,name=None):
 p=H/'inputs'/INPUT_MAP[path];pins[str(p.relative_to(H))]=sha(p);q=T/(name or Path(path).name);shutil.copy2(p,q);return pd.read_csv(q)
def data(d,n):d.to_csv(T/(n+'.csv'),index=False)
def title(ax,label,text):ax.set_title(label+'  '+text,loc='left',pad=10,fontweight='bold')
def save(fig,n):
 fig.set_size_inches(7.0,fig.get_size_inches()[1]*7.0/fig.get_size_inches()[0],forward=True)
 fig.savefig(F/(n+'.pdf'));fig.savefig(F/(n+'.svg'));fig.savefig(F/(n+'.png'),dpi=600)
 fig.savefig(H/'review'/(n+'_final_size.png'),dpi=150)
 exports[n]={e:sha(F/(n+'.'+e)) for e in ['pdf','svg','png']};plt.close(fig)
S='manuscript_execution_20260920/aim1/'
rel=read(S+'corrected_orientation/completed_release_summary.csv');fam=read(S+'corrected_patch_atlas/patch_family_assignments.csv');occ=read(S+'corrected_orientation/completed_release_residue_occupancy.csv')
occ=occ.merge(fam[['run_name','patch_family']],on='run_name');occ['resid']=occ.resid-864
atlas=occ.groupby(['patch_family','resid']).contact_occupancy_1_to_5_ns.mean().unstack().reindex(columns=range(1,77)).fillna(0)
data(atlas.reset_index(),'F1_family_mean_occupancy')
faces=['px','nx','py','ny','pz','nz'];counts=rel[rel.outcome=='retained'].groupby('face').size().reindex(faces) if 'face' in rel else None
# The release source uses starting_face, inspect its original schema rather than infer row ordering.
if counts is None:counts=rel[rel.outcome=='retained'].groupby('starting_face').size().reindex(faces)
fig=plt.figure(figsize=(7.1,3.1));gs=fig.add_gridspec(1,2,width_ratios=[1,2.25],left=.09,right=.975,bottom=.22,top=.80,wspace=.40)
a=fig.add_subplot(gs[0]);b=fig.add_subplot(gs[1]);y=np.arange(6)
a.barh(y,[3]*6,facecolor='white',edgecolor='.8',height=.65);a.barh(y,counts.values,color='.25',height=.65);a.set(yticks=y,yticklabels=['+X','−X','+Y','−Y','+Z','−Z'],xticks=[0,1,2,3],xlabel='Retained trials / 3',ylabel='Starting protein face');a.invert_yaxis();title(a,'A','Contact retention')
for i,v in enumerate(counts):a.text(v+.07,i,str(v),va='center',fontsize=8)
im=b.imshow(atlas.values,aspect='auto',cmap='Greys',vmin=0,vmax=1,extent=[.5,76.5,2.5,-.5]);b.set(yticks=[0,1,2],yticklabels=['Family A (7)','Family B (3)','Family C (3)'],xticks=[1,20,40,60,76],xlabel='Ubiquitin residue');title(b,'B','Contact location by family');cb=fig.colorbar(im,ax=b,orientation='horizontal',fraction=.075,pad=.33,aspect=32,ticks=[0,.5,1]);cb.set_label('Mean fraction of frames in contact',fontsize=8)
fig.text(.09,.94,'18 encounters → 5 ns unbiased releases → 13 retained trials',fontsize=10)
save(fig,'Figure_1_interfaces');facts['retained_trials']=int(counts.sum());facts['total_trials']=len(rel);facts['family_sizes']=fam.groupby('patch_family').size().to_dict();data(counts.rename('retained').reset_index(),'F1_face_counts')
chem=read('aim1_mechanism_completion_20260930/chemistry/RELEASE_FAMILY_CHEMISTRY.csv');native=read('aim1_mechanism_completion_20260930/chemistry/NATIVE_MATCHED_CHEMISTRY.csv');raw=read('aim1_mechanism_completion_20260930/chemistry/RELEASE_CHEMISTRY.csv')
classes=['Aliphatic','Polar','Lys/Arg','Acidic','Aromatic','His','Gly'];labels=['Aliphatic','Polar, uncharged','Lys / Arg','Asp / Glu','Aromatic','Histidine','Glycine'];n=native[(native.phase=='common_2_5ns')&np.isclose(native.lo_nm,4.3)&np.isclose(native.hi_nm,4.5)]
fig,ax=plt.subplots(1,2,figsize=(7.1,3.7),sharey=True);fig.subplots_adjust(left=.17,right=.975,bottom=.20,top=.83,wspace=.16)
for i,arm in enumerate('ABC'):
 q=chem[chem.family=='Patch '+arm].set_index('chemical_class').loc[classes]; yy=np.arange(7)+(i-1)*.20
 # Individual trajectory values make the small sample visible; black markers denote equal-trial means.
 rr=raw[raw.family=='Patch '+arm]
 for j,cl in enumerate(classes):
  vals=rr[rr.chemical_class==cl].own_conformation_enrichment.values;ax[0].scatter(vals,np.full(len(vals),yy[j]),s=7,marker=markers[arm],c='.75',alpha=.7,zorder=1)
 ax[0].scatter(q.own_enrichment,yy,s=23,marker=markers[arm],c='black',label='Family '+arm,zorder=3)
 q=n[n.arm==arm].set_index('chemical_class').loc[classes];ax[1].scatter(q.enrichment,yy,s=23,marker=markers[arm],c=C[arm],label='Preparation '+arm)
for a in ax:a.axvline(1,c='.55',ls='--',lw=.8);a.set(xlim=(-.15,7.5),xticks=[0,1,3,5,7],yticks=np.arange(7),yticklabels=labels,xlabel='Contact share / exposed-area share');a.invert_yaxis();a.legend(frameon=False,fontsize=7,loc='lower left',bbox_to_anchor=(-.03,1.13),ncol=3,handletextpad=.2,columnspacing=.5)
# Shared inverted axis should be inverted once.
ax[0].set_ylim(6.7,-.7)
title(ax[0],'A','Release families');title(ax[1],'B','Native contact preparations');fig.text(.17,.055,'Dashed line: contact composition equals exposed composition. Native comparison: 4.3–4.5 nm.',fontsize=8)
save(fig,'Figure_2_chemistry');data(n,'F2_native_contact_chemistry')
e=read('aim1_contact_mechanism_20260930/analysis/ABC_MATCHED_ENERGIES.csv');s=read('aim1_mechanism_completion_20260930/structure/MATCHED_STRUCTURE.csv');q=e[(e.phase=='common_2_5ns')&(e.band=='contact')&(e.matching=='distance_only')].set_index('arm').loc[list('ABC')]
q2=s[(s.phase=='common_2_5ns')&np.isclose(s.lo_nm,4.3)&np.isclose(s.hi_nm,4.5)].set_index('arm').loc[list('ABC')];assert np.allclose(q.buried_protein_SASA_nm2,q2.buried_protein_SASA_nm2)
data(q.reset_index(),'F3_matched_contact');fig,axs=plt.subplots(1,3,figsize=(7.1,3.0));fig.subplots_adjust(left=.095,right=.985,bottom=.24,top=.83,wspace=.62)
keys=['buried_protein_SASA_nm2','LJ_SR','LJ_per_buried_nm2_ratio_of_means'];titles=['Buried protein area','Direct LJ attraction','Attraction per area'];ys=['Area (nm²)','Mean LJ energy (kJ/mol)',r'LJ / area (kJ mol$^{-1}$ nm$^{-2}$)']
for a,key,tt,yy,letter in zip(axs,keys,titles,ys,'ABC'):
 vals=q[key].values;a.bar(range(3),vals,color=[C[v] for v in 'ABC'],width=.55);a.set(xticks=range(3),xticklabels=list('ABC'),xlabel='Native preparation',ylabel=yy);a.axhline(0,c='.4',lw=.8);title(a,letter,tt)
 for j,v in enumerate(vals):a.annotate(f'{v:.2f}' if key!= 'LJ_SR' else f'{v:.1f}',(j,v),xytext=(0,5 if v>=0 else -12),textcoords='offset points',ha='center',fontsize=8)
 if np.min(vals)<0:a.set_ylim(np.min(vals)*1.2,0)
 else:a.set_ylim(0,np.max(vals)*1.22)
fig.text(.095,.055,'Same separation (4.3–4.5 nm) and time (2–5 ns); descriptive means, not independent replicas.',fontsize=8);save(fig,'Figure_3_area_attraction')
facts['matched_contact']={a:{k:float(q.loc[a,k]) for k in keys} for a in 'ABC'}
pmf=read('bridge_endpoint_analysis_20260929/B_baseline/primary/profile.csv','B_conditional_PMF.csv');bands=read('aim1_contact_mechanism_20260930/analysis/B_BANDS.csv');groups=read('aim1_contact_mechanism_20260930/analysis/B_ATOM_GROUPS.csv');atoms=read('aim1_contact_mechanism_20260930/analysis/B_HEAVY_ATOM_CONTRIBUTORS.csv');shape=read('aim1_contact_mechanism_20260930/analysis/B_DISTANCE_PROFILE.csv')
fig,axs=plt.subplots(1,2,figsize=(7.1,3.7));fig.subplots_adjust(left=.10,right=.975,bottom=.30,top=.84,wspace=.48)
a,b=axs;a.plot(pmf.r_nm,pmf.W_wham,c='.15');a.axhline(0,c='.7',lw=.8);a.axvspan(4.3,4.5,color='.75',alpha=.5);a.axvspan(4.7,4.9,color='#E69F00',alpha=.25);a.axvspan(7,7.2,color='.85',alpha=.7);a.set(xlim=(4.25,7.25),ylim=(-14,37),xlabel='COM separation (nm)',ylabel='Conditional W(R) (kJ/mol)');title(a,'A','B distance profile')
a.annotate('Profile minimum',xy=(4.79,float(pmf.loc[(pmf.r_nm-4.79).abs().idxmin(),'W_wham'])),xytext=(5.15,13),arrowprops={'arrowstyle':'->','lw':.8},fontsize=8)
a.text(7.0,30,'Zero\nreference',ha='center',fontsize=8)
for i,(band,label,col,m) in enumerate([('contact','Closer contact\n4.3–4.5 nm','.25','o'),('PMF_minimum_band','Profile-minimum band\n4.7–4.9 nm','#B67D00','s')]):
 g=bands[bands.band==band].set_index('metric')
 if len(g)==0:raise ValueError(bands.band.unique())
 for j,key in enumerate(['LJ_SR','Coul_PME']):
  r=g.loc[key];b.errorbar(r['mean'],j+(i-.5)*.24,xerr=[[r['mean']-r.low],[r.high-r['mean']]],fmt=m,color=col,ms=5,capsize=3,label=label if j==0 else None)
  b.text(r['mean'],j+(i-.5)*.24-.12,f'{r["mean"]:.1f}',ha='center',fontsize=8,color=col)
 b.set(yticks=[0,1],yticklabels=['Lennard–Jones','Electrostatic\n(full PME)'],xlabel='Mean direct energy (kJ/mol)',xlim=(-155,10),ylim=(1.45,-.6));b.axvline(0,c='.6',lw=.8)
fig.legend(*b.get_legend_handles_labels(),frameon=False,fontsize=8,loc='center',bbox_to_anchor=(.72,.12),ncol=2,handletextpad=.35,columnspacing=.7);title(b,'B','B interaction energies')
fig.text(.10,.04,'A: illustrative profile; equilibrium uncertainty unqualified.  B: 95% conditional block intervals.',fontsize=8)
save(fig,'Figure_4_conditional_PMF');facts['B_band_means']=bands[['band','metric','mean','low','high']].to_dict('records')
g=groups[groups.band=='contact'];total=g.net.sum();facts['phenyl_percent']=100*g[g.group.str.contains('PS_phenyl')].net.sum()/total;facts['hetero_percent']=100*g[g.group.str.contains('hetero')].net.sum()/total
hyd=read('gist_B_progress_20261005/integration_review/PRIMARY_CONTRASTS.csv');density=read('gist_B_progress_20261005/integration_review/ABC_COMMON_ACCESSIBLE_DENSITY.csv');data(hyd,'F5_hydration_contrasts')
fig=plt.figure(figsize=(7.1,5.5));gs=fig.add_gridspec(2,2,left=.15,right=.97,bottom=.14,top=.93,wspace=.69,hspace=.77);a,b,c,d=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
# Diagram of the mask definition only, not a molecular snapshot or inferred water count.
for x in [.27,.73]:a.add_patch(Circle((x,.36),.20,color='#B99955',alpha=.7));a.text(x,.11,'PS',ha='center',fontsize=8)
a.add_patch(Ellipse((.27,.62),.20,.22,facecolor='.65',edgecolor='.3'));a.add_patch(Ellipse((.73,.83),.20,.22,facecolor='.65',edgecolor='.3'));a.plot([.15,.39],[.52,.52],color='black',lw=3);a.text(.27,.97,'Contact',ha='center',fontsize=8);a.text(.73,.97,'Separated',ha='center',fontsize=8);a.annotate('Blocked space',xy=(.27,.52),xytext=(.05,.79),fontsize=7,arrowprops={'arrowstyle':'->','lw':.7});a.text(.73,.63,'Ubiquitin',ha='center',fontsize=7);a.text(.5,-.10,'Same region; solute geometry defines access',ha='center',fontsize=8);a.set(xlim=(0,1),ylim=(0,1.1));a.axis('off');title(a,'A','Accessibility control (schematic)')
regions=['full_interface','blocked_interface','common_interface'];labels=['Whole interface','Space blocked\nat contact','Space accessible\nin both states']
for j,(region,label) in enumerate(zip(regions,labels)):
 r=hyd[(hyd.region==region)&(hyd.metric=='N')].iloc[0];b.errorbar(r.delta,j,xerr=[[r.delta-r.low],[r.high-r.delta]],fmt='o',color='black',capsize=3);b.text(r.delta,j-.18,f'{r.delta:+.2f}',ha='center',fontsize=8)
b.axvline(0,c='.6',lw=.8);b.set(yticks=range(3),yticklabels=labels,ylim=(2.5,-.6),xlim=(-21,5),xlabel='Water-count change');title(b,'B','Water displacement and retention')
metrics=['Esw','Eww_excess','TS_trans','TS_orient'];labs=['Solute–water','Water–water\n(relative to bulk)','−TΔS, translation','−TΔS, orientation'];z=hyd[hyd.region=='common_interface'].set_index('metric')
for j,key in enumerate(metrics):
 r=z.loc[key];sgn=-1 if key.startswith('TS') else 1;v=r.delta*sgn;lo,hi=sorted([r.low*sgn,r.high*sgn]);c.errorbar(v,j,xerr=[[v-lo],[hi-v]],fmt='o',color='black',capsize=3)
c.set(yticks=range(4),yticklabels=labs,ylim=(3.5,-.5),xlim=(-165,100),xlabel='Regional change (kJ/mol)');c.axvline(0,c='.6',lw=.8);title(c,'C','Accessible-space solvent terms')
mask=H/'inputs'/INPUT_MAP['gist_B_progress_20261005/MASK_SENSITIVITY_REVIEW.json'];pins[str(mask.relative_to(H))]=sha(mask);shutil.copy2(mask,T/mask.name);z=pd.read_json(mask);z=z[(z.metric=='proxy')&(z.region=='common_interface')].copy();z[['probe_A','delta','low','high']]=z[['probe_A','delta','low','high']].astype(float);z=z.sort_values('probe_A');d.errorbar(z.probe_A,z.delta,yerr=[z.delta-z.low,z.high-z.delta],fmt='o',color='black',capsize=3);d.axhline(0,c='.6',lw=.8);d.set(xticks=[1.2,1.4,1.6],xlim=(1.1,1.7),xlabel='Accessibility probe radius (Å)',ylabel='Regional proxy change (kJ/mol)');title(d,'D','Mask sensitivity')
fig.text(.15,.035,'Contact − separated; 95% velocity-history intervals. Regional solvent terms are not binding free energies.',fontsize=8);save(fig,'Figure_5_hydration')
facts['hydration_primary']=hyd[['region','metric','delta','low','high']].to_dict('records')
w=read('aim1_mechanism_completion_20260930/structure/WINDOW_STRUCTURE.csv');w=w[w.phase=='common_2_5ns'];data(w,'F6_window_structure')
fig,axs=plt.subplots(1,3,figsize=(7.1,3.1));fig.subplots_adjust(left=.09,right=.98,bottom=.24,top=.81,wspace=.66)
for arm in 'ABC':
 q=w[w.arm==arm].sort_values('center_nm')
 for a,key in zip(axs,['reference_Q_1.2','core_CA_RMSD_nm','ps_rg_nm']):a.plot(q.center_nm,q[key],marker=markers[arm],color=C[arm],ms=3,lw=.9,label='Preparation '+arm)
for a,tt,yy,letter in zip(axs,['Protein contact retention','Protein core deviation','Polymer size'],['Reference contact fraction Q','Cα RMSD (nm)','Radius of gyration (nm)'],'ABC'):title(a,letter,tt);a.set(xlabel='Umbrella center (nm)',ylabel=yy,xticks=[4,5,6,7])
axs[0].set_ylim(0,1.03);axs[1].set_ylim(0,.15);axs[2].set_ylim(2.60,2.73);fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,fontsize=8,bbox_to_anchor=(.55,1.0));fig.text(.09,.055,'Each point: one biased-window mean over 2–5 ns. Polymer axis is expanded; histories share ancestry.',fontsize=8)
save(fig,'Figure_6_structural_response');facts['structure_ranges']={key:[float(w[key].min()),float(w[key].max())] for key in ['reference_Q_1.2','core_CA_RMSD_nm','ps_rg_nm']}
(H/'FACTS.json').write_text(json.dumps(facts,indent=2)+'\n');(H/'SOURCE_PINS.json').write_text(json.dumps(pins,indent=2)+'\n');(H/'FIGURE_MANIFEST.json').write_text(json.dumps({'script_sha256':sha(Path(__file__)),'width_inches':7.1,'png_dpi':600,'sources':pins,'outputs':exports},indent=2)+'\n')
print('Built six figures from pinned data; no new MD, resampling or invented error bars.')
