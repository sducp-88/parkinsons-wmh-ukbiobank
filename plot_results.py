"""Publication plots from aggregate-only pd_wmh.py outputs. No individual rows."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'xtick.labelsize':8,
    'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white'})
INK='#26343D';BLUE='#3E5568';GREY='#8A959D'
def save(fig,out,name):
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/(name+'.png'),dpi=220)
    fig.savefig(out/(name+'.pdf'))
    fig.savefig(out/(name+'.svg'))
    # Raster line art at 1200 dpi, LZW compressed TIFF; vector masters retained.
    fig.savefig(out/(name+'.tif'),dpi=1200,pil_kwargs={'compression':'tiff_lzw'})
    tif=out/(name+'.tif')
    with Image.open(tif) as im:
        rgb=im.convert('RGB')
    rgb.save(tif,compression='tiff_lzw',dpi=(1200,1200));rgb.close()
    plt.close(fig)
def pvalue(p):return '<0.001' if p<.001 else f'{p:.3f}'
def forest(df,out,name,cognitive=False):
    fig,ax=plt.subplots(figsize=(6.85,3.1 if cognitive else 2.8))
    fig.subplots_adjust(left=.30,right=.70,bottom=.26,top=.83)
    y=np.arange(len(df))[::-1]
    v='beta' if cognitive else 'pct';lo='ci_low' if cognitive else 'pct_low';hi='ci_high' if cognitive else 'pct_high'
    ax.axvline(0,color=GREY,lw=.8,ls='--',zorder=0)
    ax.errorbar(df[v],y,xerr=np.array([df[v]-df[lo],df[hi]-df[v]]),fmt='s',ms=4.7,lw=1.1,
                capsize=3,color=BLUE,ecolor=BLUE)
    ax.set_yticks(y,df.outcome);ax.tick_params(axis='y',length=0,pad=8)
    ax.spines['left'].set_visible(False);ax.set_ylim(-.6,len(df)-.3)
    if cognitive:ax.set_xlim(-.60,.36);ax.set_xlabel('Adjusted difference in cognitive score (SD)')
    else:ax.set_xlim(-6,88);ax.set_xlabel('Difference in geometric mean of 1 + WMH (%)')
    fmt=lambda x:f'{x:.3f}' if cognitive else f'{x:.1f}'
    for (_,r),yy in zip(df.iterrows(),y):
        ax.text(1.045,yy,f'{fmt(r[v])} ({fmt(r[lo])}, {fmt(r[hi])})',transform=ax.get_yaxis_transform(),va='center',fontsize=8.4)
        ax.text(1.56,yy,pvalue(r.p),transform=ax.get_yaxis_transform(),va='center',fontsize=8.4)
    ax.text(-.70,1.10,'Cognitive test' if cognitive else 'WMH phenotype',transform=ax.transAxes,weight='bold',fontsize=9)
    ax.text(1.045,1.10,'Difference (95% CI)',transform=ax.transAxes,weight='bold',fontsize=8.4)
    ax.text(1.56,1.10,'P',transform=ax.transAxes,weight='bold',fontsize=8.4)
    save(fig,out,name)
def balance_plot(result,out):
    a=pd.read_csv(result/'balance_full.csv');b=pd.read_csv(result/'balance_restricted.csv')
    labels={'Age_at_Instance_2':'Age at MRI','Townsend_deprivation_index_at_recruitment':'Townsend index',
      'Body_mass_index_BMI_Instance_0':'Body mass index','Alcohol_intake_frequency_ordinal':'Alcohol frequency',
      'C(Sex)[T.1]':'Male sex','C(Genetic_ethnic_grouping)[T.1]':'Genetically White British',
      'C(Smoking_Ever)[T.1.0]':'Ever smoking','C(has_degree)[T.1.0]':'Degree qualification',
      'C(CMC_score_cat)[T.1]':'One recorded comorbidity','C(CMC_score_cat)[T.2]':'Two or more recorded comorbidities'}
    fig,axes=plt.subplots(1,2,figsize=(6.85,4.2),sharey=True)
    fig.subplots_adjust(left=.37,right=.98,bottom=.19,top=.87,wspace=.15)
    for i,(ax,d,name) in enumerate(zip(axes,[a,b],['Full cohort','Restricted cohort'])):
        y=np.arange(len(d))[::-1]
        ax.scatter(d.smd_before.abs(),y,s=23,facecolors='none',edgecolors=GREY,marker='o',label='Before weighting')
        ax.scatter(d.smd_after.abs(),y,s=21,color=BLUE,marker='s',label='After weighting')
        ax.axvline(.1,color=GREY,lw=.8,ls='--');ax.set_xlim(-.025,.75)
        ax.set_yticks(y,[labels.get(x,x) for x in d.covariate]);ax.tick_params(axis='y',length=0)
        ax.spines['left'].set_visible(False);ax.set_xlabel('Absolute SMD')
        ax.text(0,1.04,chr(65+i)+'  '+name,transform=ax.transAxes,weight='bold',fontsize=9)
    axes[1].legend(loc='upper center',bbox_to_anchor=(-.05,-.19),ncol=2,frameon=False,fontsize=8)
    save(fig,out,'Figure_S1_Covariate_balance')
def pathways(result,out):
    d=pd.read_csv(result/'path_products.csv');names=list(dict.fromkeys(d.outcome))
    fig,axes=plt.subplots(2,2,figsize=(6.85,4.9),sharex=True,sharey=True)
    fig.subplots_adjust(left=.27,right=.98,bottom=.14,top=.92,hspace=.48,wspace=.12)
    for i,(ax,label) in enumerate(zip(axes.flat,names)):
        q=d.loc[d.outcome==label];y=np.arange(len(q))[::-1]
        ax.axvline(0,color=GREY,lw=.8,ls='--')
        ax.errorbar(q.ab,y,xerr=np.array([q.ab-q.ci_low,q.ci_high-q.ab]),fmt='s',ms=4,lw=1.1,capsize=3,color=BLUE)
        ax.set_yticks(y,q.mediator);ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
        ax.set_xlim(-.034,.004);ax.set_xticks([-.03,-.02,-.01,0]);ax.set_ylim(-.5,2.6)
        ax.text(0,1.06,chr(65+i)+'  '+label,transform=ax.transAxes,weight='bold',fontsize=9)
        if i>=2:ax.set_xlabel('Statistical path product (SD)')
    save(fig,out,'Figure_S2_Exploratory_path_products')
def flowchart(flow,out):
    d=json.loads(flow.read_text(encoding='utf-8'))
    fig,ax=plt.subplots(figsize=(6.85,5.2));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    def box(x,y,w,h,text,fontsize=9):
        ax.add_patch(Rectangle((x-w/2,y-h/2),w,h,facecolor='white',edgecolor=INK,lw=.8))
        ax.text(x,y,text,ha='center',va='center',fontsize=fontsize,color=INK)
    def arrow(x1,y1,x2,y2):ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'->','lw':.9,'color':INK})
    box(.28,.92,.50,.10,f"UK Biobank extract\nn = {d['extract_n']:,}")
    box(.80,.80,.36,.11,f"Total WMH unavailable\nn = {d['wmh_missing_n']:,}")
    arrow(.28,.87,.28,.74);arrow(.28,.80,.62,.80)
    box(.28,.69,.50,.10,f"Total WMH available\nn = {d['wmh_available_n']:,}")
    box(.80,.58,.36,.10,f"PD/control definition not met\nn = {d['definition_excluded_n']:,}",8.3)
    arrow(.28,.64,.28,.52);arrow(.28,.58,.62,.58)
    box(.28,.46,.50,.12,f"Full analysis cohort  n = {d['full_n']:,}\nRecorded PD  n = {d['pd_n']:,}\nControls  n = {d['control_n']:,}")
    arrow(.28,.40,.28,.30)
    box(.28,.24,.50,.12,f"Restricted cohort  n = {d['restricted_n']:,}\nRecorded PD  n = {d['restricted_pd_n']:,}\nControls  n = {d['restricted_control_n']:,}")
    box(.80,.35,.36,.12,f"Recorded neurologic codes\nbefore or at MRI\nExcluded  n = {d['neurologic_excluded_n']:,}",8.5)
    arrow(.28,.35,.62,.35)
    ax.text(.03,.07,f"Regional WMH available\nFull cohort: {d['regional_full_n']:,} ({d['pd_n']} PD)\nRestricted cohort: {d['regional_restricted_n']:,} ({d['restricted_pd_n']} PD)",fontsize=8.5,va='center')
    save(fig,out,'Figure_1_Cohort_flow')
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--flow-json',type=Path)
    a=p.parse_args();d=pd.read_csv(a.results/'wmh_results.csv')
    forest(d.loc[d.cohort=='Full overlap weighted'],a.output,'Figure_2_PD_WMH')
    c=pd.read_csv(a.results/'cognitive_results.csv')
    forest(c.loc[~c.weighted],a.output,'Figure_3_PD_Cognition',True)
    balance_plot(a.results,a.output);pathways(a.results,a.output)
    if a.flow_json:flowchart(a.flow_json,a.output)
if __name__=='__main__':main()
