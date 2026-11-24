# Both figures are created first and shown together by a single plt.show() at the end; they are also saved to results/
import numpy as np, pandas as pd, matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LightSource
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter, uniform_filter, distance_transform_edt
from scipy import stats
import os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA=os.path.join(ROOT,'data')+os.sep; RES=os.path.join(ROOT,'results')+os.sep
os.makedirs(DATA,exist_ok=True); os.makedirs(RES,exist_ok=True)
N=100; cs=0.1  # 100 m cells -> 10 km x 10 km = 100 km2
yy,xx=np.mgrid[0:N,0:N]
def sn(sig,seed):
    g=gaussian_filter(np.random.default_rng(seed).normal(size=(N,N)),sig,mode='reflect'); return (g-g.mean())/g.std()
z=lambda a:(a-a.mean())/a.std()
nm=lambda a:(a-a.min())/(a.max()-a.min())
# ---- synthetic DEM, slope, D8 flow accumulation
dem=900-250*(xx/N*0.6+yy/N*0.4)+120*np.exp(-((xx-15)**2+(yy-15)**2)/(2*20**2))+25*sn(6,1)+10*sn(3,2)
gy,gx=np.gradient(dem,cs*1000); slope=np.hypot(gx,gy)*100
acc=np.ones((N,N)); nb=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
for idx in np.argsort(-dem.ravel()):
    r,c=divmod(idx,N); best=None; bd=0
    for dr,dc in nb:
        rr,cc=r+dr,c+dc
        if 0<=rr<N and 0<=cc<N:
            d=(dem[r,c]-dem[rr,cc])/np.hypot(dr,dc)
            if d>bd: bd=d; best=(rr,cc)
    if best: acc[best]+=acc[r,c]
stream=acc>=50
dd=uniform_filter(stream.astype(float),11)*10   # km/km2
# ---- other synthetic layers
rain=1100-450*(xx/N)+40*sn(8,7)                  # mm, west->east rain-shadow gradient
n1=sn(5,3); n4=sn(3,6)
lulc=np.full((N,N),5)                             # 1 built-up,2 barren,3 scrub,4 forest,5 agri
lulc[(slope>6)&(n1>0)]=3; lulc[(dem>np.quantile(dem,.85))]=4
lulc[(slope>9)&(sn(4,8)>0.3)]=2; lulc[n4>1.7]=1
soil_raw=sn(7,4)-0.3*z(slope); q=np.quantile(soil_raw,[.2,.5,.8])
soil=5-np.digitize(soil_raw,q)                    # 5 sandy loam ... 2 clay
gm_raw=-0.7*z(slope)+0.5*z(np.log(acc))+0.5*sn(6,5); q=np.quantile(gm_raw,[.15,.35,.65,.85])
geom=1+np.digitize(gm_raw,q)                      # 1 hill/ridge ... 5 valley fill
rng=np.random.default_rng(11); lin=np.zeros((N,N))
for _ in range(45):
    r0,c0=rng.uniform(0,N,2); a=np.deg2rad(rng.normal(135,25)); L=rng.uniform(10,30)
    for t in np.linspace(0,L,int(L*3)):
        r,c=int(r0+t*np.sin(a)),int(c0+t*np.cos(a))
        if 0<=r<N and 0<=c<N: lin[r,c]=1
ld=uniform_filter(lin,11)*10
# ---- reclassification to scores 1-5 (5 = high recharge)
qs=lambda a,rev=False:(np.digitize(a,np.quantile(a,[.2,.4,.6,.8]))+1 if not rev else 5-np.digitize(a,np.quantile(a,[.2,.4,.6,.8])))
S={'Geology/geomorphology':geom,
   'Slope':np.select([slope<2,slope<5,slope<10,slope<15],[5,4,3,2],1),
   'Soil texture':soil,'Lineament density':qs(ld),
   'Rainfall':np.digitize(rain,[650+90*i for i in range(1,5)])+1,
   'Drainage density':qs(dd,True),
   'LULC':np.select([lulc==5,lulc==4,lulc==3,lulc==2,lulc==1],[4,4,3,2,1])}
# ---- AHP (Saaty 1-9; importance scores -> pairwise matrix)
names=list(S); imp=np.array([8,7,6,5,5,4,3])
A=np.array([[1+(i-j) if i>=j else 1/(1+(j-i)) for j in imp] for i in imp],float)
ev,evec=np.linalg.eig(A); k=np.argmax(ev.real); w=np.abs(evec[:,k].real); w/=w.sum()
lam=ev.real[k]; CI=(lam-7)/6; CR=CI/1.32
print('AHP weights:',dict(zip(names,np.round(w*100,1)))); print('lambda_max=%.3f CI=%.3f CR=%.3f'%(lam,CI,CR))
gwrpi=sum(w[i]*S[n] for i,n in enumerate(names))
edges=np.linspace(gwrpi.min(),gwrpi.max(),6); zone=np.clip(np.digitize(gwrpi,edges[1:-1])+1,1,5)
zn=['Very low','Low','Moderate','High','Very high']
print('GWRPI range %.2f-%.2f'%(gwrpi.min(),gwrpi.max()))
for i in range(1,6): print(zn[i-1],'%.1f km2 (%.1f%%)'%((zone==i).sum()*cs*cs,(zone==i).mean()*100))
# ---- synthetic observation wells (independent 'truth' = perturbed weights + local effects)
cand=np.argwhere(lulc!=1); sel=[]
for i in rng.permutation(len(cand)):
    p=cand[i]
    if all(np.hypot(*(p-q))>12 for q in sel): sel.append(p)
    if len(sel)==28: break
sel=np.array(sel); wr,wc=sel[:,0],sel[:,1]
w2=rng.dirichlet(w*60); truth=sum(w2[i]*S[n] for i,n in enumerate(names))+0.35*sn(4,9)
tn=(truth-truth.min())/(truth.max()-truth.min())
fluct=np.clip(1.5+5.5*tn[wr,wc]+rng.normal(0,.9,28),.2,None)      # m, post- minus pre-monsoon rise
depth=6+14*(1-tn[wr,wc])+rng.normal(0,2,28)                         # m bgl pre-monsoon
rho,pv=stats.spearmanr(gwrpi[wr,wc],fluct); print('Spearman rho=%.2f p=%.4f n=28'%(rho,pv))
wells=pd.DataFrame(dict(row=wr,col=wc,x_km=wc*cs,y_km=(N-1-wr)*cs,GWRPI=gwrpi[wr,wc].round(2),zone=[zn[z_-1] for z_ in zone[wr,wc]],fluct_m=fluct.round(2),pre_monsoon_depth_m=depth.round(1)))
wells.to_csv(DATA+'synthetic_wells.csv',index=False)
print(wells.groupby('zone').fluct_m.agg(['count','mean']).round(2))
# ---- need (IDW of depth) and priority score
d2=(yy[...,None]-wr)**2+(xx[...,None]-wc)**2+1e-6; wt=1/d2; need=(wt*depth).sum(-1)/wt.sum(-1)
score=0.5*nm(gwrpi)+0.3*nm(np.log(acc))+0.2*nm(need)
dist=distance_transform_edt(~stream)
masks={'Check dam':stream&(acc>=80)&(acc<=1000)&(slope<10)&(zone>=3)&(lulc!=1),
 'Percolation tank':(dist<=3)&(slope<5)&(zone>=4)&(soil>=3)&(geom>=3)&(lulc!=1)&(acc>=30),
 'Farm pond':(lulc==5)&(slope>=1)&(slope<=6)&(soil<=3)&(zone>=2)&(acc>=5)&(acc<=40)&(dist>3)}
quota={'Check dam':5,'Percolation tank':4,'Farm pond':6}; sites=[]; allchosen=[]
edge=np.zeros((N,N),bool); edge[6:-6,6:-6]=True
for t in masks: masks[t]=masks[t]&edge
for t,m in masks.items():
    print(t,'candidate cells:',m.sum()); idx=np.argwhere(m); order=idx[np.argsort(-score[m])] if False else None
    cells=sorted(map(tuple,idx),key=lambda p:-score[p]); chosen=[]
    for p in cells:
        if all(np.hypot(p[0]-q[0],p[1]-q[1])>8 for q in chosen) and all(np.hypot(p[0]-q[0],p[1]-q[1])>6 for q in allchosen): chosen.append(p)
        if len(chosen)==quota[t]: break
    allchosen+=chosen
    for p in chosen: sites.append(dict(type=t,row=p[0],col=p[1],x_km=round(p[1]*cs,1),y_km=round((N-1-p[0])*cs,1),
        GWRPI=round(gwrpi[p],2),zone=zn[zone[p]-1],slope_pct=round(slope[p],1),catchment_km2=round(acc[p]*cs*cs,2),priority_score=round(score[p],3)))
sd=pd.DataFrame(sites).sort_values('priority_score',ascending=False).reset_index(drop=True); sd.insert(0,'rank',sd.index+1)
sd['class']=np.where(sd['rank']<=5,'Priority 1',np.where(sd['rank']<=10,'Priority 2','Priority 3'))
sd.to_csv(RES+'ranked_recharge_sites.csv',index=False); print(sd.to_string())
# ---- FIGURE 1
cols=['#a6611a','#dfc27d','#f5f0c8','#80cdc1','#018571']
fig,ax=plt.subplots(figsize=(10.5,8.2)); ext=[0,10,0,10]
hs=LightSource(315,45).hillshade(dem,vert_exag=3)
ax.imshow(hs,cmap='gray',extent=ext,origin='upper')
ax.imshow(zone,cmap=ListedColormap(cols),extent=ext,origin='upper',alpha=.72,vmin=1,vmax=5)
sm=np.ma.masked_where(~(acc>=150),np.ones((N,N))); ax.imshow(sm,cmap=ListedColormap(['#08306b']),extent=ext,origin='upper')
ax.scatter(wc*cs+.05,(N-1-wr)*cs+.05,marker='^',c='k',s=45,edgecolor='w',lw=.6,zorder=5)
mk={'Check dam':('s','#e41a1c'),'Percolation tank':('*','#ff7f00'),'Farm pond':('o','#984ea3')}
for _,r in sd.iterrows():
    m,c=mk[r['type']]; ax.scatter(r.x_km+.05,r.y_km+.05,marker=m,s=330 if m=='*' else 190,c=c,edgecolor='k',lw=1,zorder=6)
    ax.text(r.x_km+.05,r.y_km+.05,str(r['rank']),ha='center',va='center',fontsize=7.5,color='w',weight='bold',zorder=7)
ax.set_xlabel('Easting within sub-watershed (km)',fontsize=12); ax.set_ylabel('Northing within sub-watershed (km)',fontsize=12)
ax.set_title('Fig. 1  Groundwater recharge-potential zones and ranked recharge structures\n(synthetic Upper Bhima-type sub-watershed, 100 km², 100 m cells)',fontsize=12)
h=[plt.Rectangle((0,0),1,1,fc=c) for c in cols]+[Line2D([],[],color='#08306b',lw=2),Line2D([],[],marker='^',color='k',ls='',mfc='k',ms=8)]+[Line2D([],[],marker=m,color=c,ls='',mec='k',ms=10) for m,c in mk.values()]
l=[f'{n} ({(zone==i+1).mean()*100:.0f}%)' for i,n in enumerate(zn)]+['Streams','Observation wells']+[f'{k} (no. = rank)' for k in mk]
ax.legend(h,l,loc='upper left',bbox_to_anchor=(1.01,1),title='GWRPI zone (% area)',fontsize=10,title_fontsize=11,frameon=False)
plt.tight_layout(); plt.savefig(RES+'Fig1_recharge_potential_map.png',dpi=200,bbox_inches='tight')
# ---- FIGURE 2
fig,(a,b)=plt.subplots(1,2,figsize=(12,5.2))
a.scatter(wells.GWRPI,wells.fluct_m,c='#2c7fb8',s=55,edgecolor='k')
m_,c_=np.polyfit(wells.GWRPI,wells.fluct_m,1); xs=np.linspace(wells.GWRPI.min(),wells.GWRPI.max(),10); a.plot(xs,m_*xs+c_,'r--')
a.set_xlabel('Groundwater recharge-potential index (GWRPI, 1-5)',fontsize=12); a.set_ylabel('Water-level rise, pre- to post-monsoon (m)',fontsize=12)
a.set_title(f'(a) Index vs. well response (n=28)\nSpearman ρ = {rho:.2f}, '+('p < 0.001' if pv<0.001 else f'p = {pv:.3f}')+f'',fontsize=12); a.tick_params(labelsize=11)
order=[z_ for z_ in zn if z_ in set(wells.zone)]; b.boxplot([wells[wells.zone==z_].fluct_m for z_ in order],patch_artist=True,boxprops=dict(facecolor='#c7e9c0'))
b.set_xticks(range(1,len(order)+1)); b.set_xticklabels([f'{z_}\n(n={(wells.zone==z_).sum()})' for z_ in order])
b.set_xlabel('Recharge-potential zone',fontsize=12); b.set_ylabel('Water-level rise (m)',fontsize=12); b.set_title('(b) Well response by zone',fontsize=12); b.tick_params(labelsize=11)
fig.suptitle('Fig. 2  Consistency check of recharge zones against (synthetic) observation wells',fontsize=13,y=1.0)
plt.tight_layout(); plt.savefig(RES+'Fig2_well_validation.png',dpi=200,bbox_inches='tight')
print('Figures created:',plt.get_fignums(),'-> saved in',RES)
plt.show()   # one show() at the end displays BOTH figures (Fig 1 and Fig 2)
