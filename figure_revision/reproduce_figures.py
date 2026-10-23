#!/usr/bin/env python3
"""Recreate repaired SI Figures S10 and S12 from supplied, hash-pinned data.
No trajectory extraction, resampling, or scientific selection is performed.
"""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cut_tree
H=Path(__file__).resolve().parent; I=H/'inputs'; O=H/'figures'; O.mkdir(exist_ok=True)
def read(name):
 with (I/name).open(newline='') as f:return list(csv.DictReader(f))
def save(fig,name):
 for ext in ['pdf','svg','png']:
  fig.savefig(O/f'{name}.{ext}',dpi=400 if ext=='png' else None,bbox_inches='tight')
 plt.close(fig)
for p in json.loads((H/'SOURCE_PINS.json').read_text()):
 assert hashlib.sha256((H/p['copied']).read_bytes()).hexdigest()==p['sha256']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.titlesize':11,'axes.linewidth':.7})
original_assignments=read('patch_family_assignments.csv'); a=original_assignments.copy(); a.sort(key=lambda x:(x['patch_family'],x['is_family_medoid']!='True',x['display_id']))
assert len(a)==13 and [sum(x['patch_family']==f'Patch {f}' for x in a) for f in 'ABC']==[7,3,3]
occ={(x['run_name'],int(x['resid'])-864):float(x['contact_occupancy_1_to_5_ns']) for x in read('completed_release_residue_occupancy.csv')}
values=np.array([[occ[(x['run_name'],i)] for i in range(1,77)] for x in a]); assert values.shape==(13,76) and np.all((values>=0)&(values<=1))
dist=squareform(pdist(values>=.2,metric='jaccard')); sim=1-dist
primary=np.array([x['patch_family'][-1] for x in original_assignments]); source_values=np.array([[occ[(x['run_name'],i)] for i in range(1,77)] for x in original_assignments]); pairs=np.triu_indices(13,1); primary_same=(primary[:,None]==primary)[pairs]
results=[]; detail=[]
for r in read('patch_family_sensitivity.csv'):
 d=pdist(source_values>=float(r['persistent_contact_threshold']),metric='jaccard'); labels=cut_tree(linkage(d,method=r['linkage']),n_clusters=3).ravel(); same=(labels[:,None]==labels)[pairs]; agreed=int(np.sum(same==primary_same)); fraction=agreed/78
 assert abs(fraction-float(r['pairwise_agreement_with_primary']))<1e-14
 results.append(dict(threshold=float(r['persistent_contact_threshold']),linkage=r['linkage'],agreed_pairs=agreed,total_pairs=78,agreement=fraction))
 detail.extend(dict(threshold=float(r['persistent_contact_threshold']),linkage=r['linkage'],display_id=x['display_id'],primary_family=x['patch_family'][-1],sensitivity_cluster=int(l)) for x,l in zip(original_assignments,labels))
geometry=read('contact_patch_geometry_residue_motifs.csv'); verified=[]
for r in read('patch_family_geometry.csv'):
 runs=[x['run_name'] for x in a if x['patch_family']==r['patch_family']]
 v=[max([float(g['motif_frame_fraction']) for g in geometry if g['run_name']==run and g['motif']==r['motif']] or [0]) for run in runs]
 assert abs(np.mean(v)-float(r['mean_trajectory_motif_fraction']))<1e-14
 assert np.count_nonzero(v)==int(r['trajectories_with_geometry'])
 verified.append(dict(family=r['patch_family'],motif=r['motif'],mean_max_residue_occupancy=float(np.mean(v)),trajectories_with_geometry=int(np.count_nonzero(v)),trajectories=len(v)))
for name,rows in [('CLUSTER_SENSITIVITY_RECOMPUTED.csv',results),('CLUSTER_MEMBERSHIP_RECOMPUTED.csv',detail),('RING_GEOMETRY_RECOMPUTED.csv',verified)]:
 with (H/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
labels=[x['display_id'] for x in a]
fig=plt.figure(figsize=(8,8.15),layout='constrained'); outer=fig.add_gridspec(2,1,height_ratios=[1,1.05]); top=outer[0].subgridspec(1,2,width_ratios=[1.5,.75],wspace=.28)
ax=fig.add_subplot(top[0]); image=ax.pcolormesh(np.arange(14)-.5,np.arange(14)-.5,sim,cmap='Greys',vmin=0,vmax=1,shading='flat',rasterized=False);ax.set_ylim(12.5,-.5);ax.set_aspect('equal');ax.set_xticks(range(13),labels,rotation=45,ha='right',fontsize=8);ax.set_yticks(range(13),labels,fontsize=8);ax.set_title('A  Shared persistent residues',loc='left',pad=36,fontweight='bold')
for lo,n in [(0,7),(7,3),(10,3)]:ax.add_patch(plt.Rectangle((lo-.5,lo-.5),n,n,fill=False,edgecolor='#222222',lw=1.1))
sax=ax.secondary_xaxis('top');sax.set_xticks([3,8,11],['Family A\n(n=7)','Family B\n(n=3)','Family C\n(n=3)']);sax.tick_params(length=0,labelsize=8);sax.spines['top'].set_visible(False)
fig.colorbar(image,ax=ax,fraction=.035,pad=.025).set_label('Jaccard similarity')
ax=fig.add_subplot(top[1]);faces=['PX','NX','PY','NY','PZ','NZ'];counts=np.array([[sum(x['starting_face']==face and x['patch_family']==f'Patch {f}' for x in a) for f in 'ABC'] for face in faces]); assert counts.sum()==13
# Vector cells replace the malformed narrow embedded raster of the former panel B.
ax.pcolormesh(np.arange(4)-.5,np.arange(7)-.5,counts,cmap='Greys',vmin=0,vmax=3,shading='flat',edgecolors='#bbbbbb',linewidth=.3,rasterized=False);ax.set_ylim(5.5,-.5);ax.set_aspect('equal');ax.set_xticks(range(3),['A','B','C']);ax.set_yticks(range(6),faces);ax.set_xlabel('Release family');ax.set_ylabel('Starting face');ax.set_title('B  Starting-face outcomes',loc='left',pad=11,fontweight='bold')
for i in range(6):
 for j in range(3):ax.text(j,i,str(counts[i,j]),ha='center',va='center',color='white' if counts[i,j]>=2 else '#222222',weight='bold')
ax=fig.add_subplot(outer[1]); image=ax.pcolormesh(np.arange(77)+.5,np.arange(14)-.5,values,cmap='Greys',vmin=0,vmax=1,shading='flat',rasterized=False);ax.set_ylim(12.5,-.5);ax.set_yticks(range(13),[f'{x["display_id"]} [{x["patch_family"][-1]}]' for x in a],fontsize=8);ax.set_xticks([1,10,20,30,40,50,60,70,76]);ax.set_xlabel('Ubiquitin residue number');ax.set_title('C  Residue contact occupancy',loc='left',pad=10,fontweight='bold')
for b in [6.5,9.5]:ax.axhline(b,color='#222222',lw=.9)
fig.colorbar(image,ax=ax,fraction=.025,pad=.015).set_label('Fraction of analyzed frames in contact')
save(fig,'figure_s10')
# S12 uses exactly the existing plotted points and outlines. Only metric naming changes.
colors={'A':'#0072B2','B':'#D55E00','C':'#009E73'}
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(7.35,6.2));fig.subplots_adjust(wspace=.34,hspace=.48,left=.10,right=.98,bottom=.095,top=.93)
for panel,ax,plane in [('a',axes[0,0],'xy'),('b',axes[0,1],'xz')]:
 for arm in 'ABC':
  z=np.load(I/f'{arm}_4.8_all_external.npz');v=z['occupancy'];o=z['origin_A'];s=z['spacing_A'];xyz=[(o[i]+np.arange(v.shape[i])*s[i])/10 for i in range(3)]
  if plane=='xy':data=v[:,:,np.argmin(abs(xyz[2]))];vertical=xyz[1]
  else:data=v[:,np.argmin(abs(xyz[1])),:];vertical=xyz[2]
  ax.contour(xyz[0],vertical,data.T,levels=[.5],colors=[colors[arm]],linewidths=1.35)
 ax.set_aspect('equal');ax.set_xlim(-4.5,4.5);ax.set_ylim(-4.5,4.5);ax.set_xticks([-4,-2,0,2,4]);ax.set_yticks([-4,-2,0,2,4]);ax.set_xlabel('Particle-fixed x (nm)');ax.set_ylabel('Particle-fixed '+plane[1]+' (nm)');ax.set_title('Local surface outlines: '+plane+' section',loc='left',pad=8);ax.text(-.19,1.06,panel,transform=ax.transAxes,fontweight='bold',fontsize=13)
fig.legend([Line2D([0],[0],color=colors[a],lw=2) for a in 'ABC'],['A','B','C'],loc='upper center',bbox_to_anchor=(.53,.995),frameon=False,ncol=3,handlelength=1.4,columnspacing=2,borderpad=0)
r=read('SURFACE_COMPARISONS.csv');ax=axes[1,0]
for j,center in enumerate((4.4,4.8,7.2)):
 within=[float(x['redistributed_fraction'])*100 for x in r if float(x['center_nm'])==center and x['comparison']=='within_preparation_early_late'];between=[float(x['redistributed_fraction'])*100 for x in r if float(x['center_nm'])==center and x['comparison']=='between_preparation_all'];assert len(within)==len(between)==3
 ax.scatter(j-.13+np.array([-.035,0,.035]),within,color='#777777',s=31,marker='o',zorder=3);ax.scatter(j+.13+np.array([-.035,0,.035]),between,color='#5A3D8A',s=35,marker='D',zorder=3)
ax.set_xticks(range(3),['4.4','4.8','7.2\n(separated)']);ax.set_ylabel('Normalized envelope-map\ndifference (%)');ax.set_xlabel('Window center (nm)');ax.set_ylim(0,10.8);ax.grid(axis='y',color='.9',zorder=0);ax.set_title('Between preparations > within one',loc='left',pad=10);ax.text(-.19,1.06,'c',transform=ax.transAxes,fontweight='bold',fontsize=13)
ax.legend([Line2D([0],[0],marker='o',linestyle='none',color='#777777'),Line2D([0],[0],marker='D',linestyle='none',color='#5A3D8A')],['Early vs late (same trajectory)','A/B, A/C and B/C'],frameon=False,fontsize=8,loc='upper left',borderpad=0)
r=read('SASA_CONTROL.csv');ax=axes[1,1]
for ai,arm in enumerate('ABC'):
 means=[]
 for j,center in enumerate((4.4,4.8,7.2)):
  v=[float(x['SASA_nm2']) for x in r if x['key']==f'{arm}_{center:.1f}' and int(x['sphere_points'])==960];assert len(v)==3;means.append(np.mean(v));ax.scatter(j+(ai-1)*.13+np.array([-.025,0,.025]),v,color=colors[arm],s=17,alpha=.8,zorder=3)
 ax.plot(np.arange(3)+(ai-1)*.13,means,color=colors[arm],lw=1,label=arm)
ax.set_xticks(range(3),['4.4','4.8','7.2']);ax.set_ylim(250,315);ax.set_ylabel('Polymer-only accessible area (nm²)');ax.set_xlabel('Window center (nm)');ax.grid(axis='y',color='.9');ax.set_title('Surface area also differs',loc='left',pad=10);ax.text(-.19,1.06,'d',transform=ax.transAxes,fontweight='bold',fontsize=13);ax.legend(frameon=False,ncol=3,loc='upper right',handlelength=1,columnspacing=.8)
save(fig,'figure_s12')
(H/'FIGURE_CHECKS.json').write_text(json.dumps(dict(status='PASS_NUMERICAL_SOURCE_CHECKS',S10_retained_trajectories=13,S10_residue_observations=988,S10_family_sizes=[7,3,3],S10_face_counts=counts.tolist(),S10_vector_heatmaps=True,clustering_six_sensitivity_rows_recomputed=True,geometry_six_max_residue_summaries_recomputed=True,S12_comparison_points=18,S12_SASA_points=27,S12_original_inputs_preserved=True,claim_limit='Revised presentation of existing data; no new trajectory or clustering-selection rule.'),indent=2)+'\n')
print('Repaired S10 and S12; numerical source checks PASS.')
