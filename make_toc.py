from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/aim1-freeze-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
H=Path(__file__).resolve().parent
fig=plt.figure(figsize=(3.25,1.75));fig.text(.5,.94,'One disordered polystyrene surface',ha='center',fontsize=8,fontweight='bold')
for j,stem in enumerate(['patch_a_nx-c_medoid','patch_b_py-c_medoid','patch_c_pz-c_medoid']):
 ax=fig.add_axes([.02+j*.325,.22,.31,.59]);im=plt.imread(H/'molecular_inputs'/(stem+'.png'));rows,cols=np.where(np.any(im[:,:,:3]<.985,axis=2));ax.imshow(im[max(0,rows.min()-10):rows.max()+11,max(0,cols.min()-10):cols.max()+11]);ax.axis('off')
fig.text(.5,.12,'Multiple ubiquitin interfaces',ha='center',fontsize=8)
fig.text(.5,.025,'Attractive contacts with retained hydration',ha='center',fontsize=7)
for ext in ['pdf','svg','png']:fig.savefig(H/('TOC_graphic.'+ext),dpi=600)
