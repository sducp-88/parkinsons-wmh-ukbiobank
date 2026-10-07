"""Aggregate-only analysis of Parkinson disease and WMH in an authorised UKB workspace.

Input is an author-prepared, consent-filtered canonical CSV; no participant rows
are exported. This release repairs the archived propensity model and diagnosis
column aliases. It does not validate the upstream consent or access approval.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sm
from patsy import dmatrices, dmatrix
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

CORE=['Sex','Age_at_Instance_2','Townsend_deprivation_index_at_recruitment',
      'Body_mass_index_BMI_Instance_0','Genetic_ethnic_grouping','Smoking_Ever',
      'Alcohol_intake_frequency_ordinal','has_degree','CMC_score_cat']
CAT={'Sex','Genetic_ethnic_grouping','Smoking_Ever','has_degree','CMC_score_cat'}
SCALE='Volumetric_scaling_from_T1_head_image_to_standard_space_Instance_2'
MRI='Date_of_attending_assessment_centre_Instance_2'
PDDATE='Date_of_parkinson_s_disease_report'
ALLDATE='Date_of_all_cause_parkinsonism_report'
RAW={
    'Total WMH':'Total_volume_of_white_matter_hyperintensities_from_T1_and_T2_FLAIR_images_Instance_2',
    'Periventricular WMH':'Total_volume_of_peri_ventricular_white_matter_hyperintensities_Instance_2',
    'Deep WMH':'Total_volume_of_deep_white_matter_hyperintensities_Instance_2',
}
COG={
    'Reaction time':'Mean_time_to_correctly_identify_matches_Instance_2',
    'Trail Making Test B':'Duration_to_complete_alphanumeric_path_trail_2_Instance_2',
    'Fluid intelligence':'Fluid_intelligence_score_Instance_2',
    'Digit span':'Maximum_digits_remembered_correctly_Instance_2',
}
NEURO=('I63','I64','G45','G35','I67')
CMC={
    'CMC_hypertension':('I10','I11','I12','I13','I15'),
    'CMC_hyperlipidemia':('E78',),
    'CMC_arrhythmia':('I47','I48','I49'),
    'CMC_CAD':('I20','I21','I22','I23','I24','I25'),
    'CMC_heart_failure':('I50',),
    'CMC_diabetes':('E10','E11','E12','E13','E14'),
    'CMC_stroke':('I60','I61','I62','I63','I64'),
}

def require_columns(df,cols):
    absent=sorted(set(cols)-set(df.columns))
    if absent: raise ValueError('Missing required schema columns: '+', '.join(absent))

def prepare(df,repair_cmc=True):
    df=df.copy()
    require_columns(df,CORE+['e4_count',SCALE,MRI,PDDATE,ALLDATE]+list(RAW.values())+list(COG.values()))
    if 'treatment_var' not in df:
        require_columns(df,['group'])
        df['treatment_var']=df['group'].map({'Control':0,'Study':1})
    if df['treatment_var'].isna().any() or not df['treatment_var'].isin([0,1]).all():
        raise ValueError('Exposure must have only 0/1 values.')
    if 'Participant_ID' in df and df.Participant_ID.duplicated().any():
        raise ValueError('Participant identifiers are not unique.')
    t=df.treatment_var.astype(int)
    dates={c:pd.to_datetime(df[c],errors='coerce') for c in [MRI,PDDATE,ALLDATE]}
    unknown=(dates[PDDATE]==pd.Timestamp('1900-01-01'))|(dates[ALLDATE]==pd.Timestamp('1900-01-01'))
    valid_case=(dates[PDDATE].notna()&dates[ALLDATE].notna()&dates[MRI].notna()&(dates[PDDATE]<=dates[MRI])&~unknown)
    valid_control=dates[PDDATE].isna()&dates[ALLDATE].isna()
    if ((t==1)&~valid_case).any() or ((t==0)&~valid_control).any():
        raise ValueError('Archived exposure fails date-anchored PD/control definition; update upstream cohort locally.')
    for c in CORE+['e4_count',SCALE]+list(RAW.values())+list(COG.values()):
        df[c]=pd.to_numeric(df[c],errors='coerce')
    imputation=[]
    for c in CORE+['e4_count']:
        n=int(df[c].isna().sum())
        if not n:continue
        if df[c].notna().sum()==0:raise ValueError('Entirely missing covariate '+c)
        strategy='mode' if c in CAT or c=='e4_count' else 'median'
        fill=float(df[c].mode().iloc[0]) if strategy=='mode' else float(df[c].median())
        df[c]=df[c].fillna(fill)
        imputation.append({'variable':c,'missing_n':n,'method':strategy,'fill':fill})
    if (~np.isfinite(df[SCALE])|(df[SCALE]<=0)).any():raise ValueError('Invalid scaling factor.')
    cmc_audit={'repaired':repair_cmc,'available_codes':{},'absent_codes':{},'changed_count':0,'additional_imputation':imputation}
    old=df.CMC_score_cat.copy()
    if repair_cmc:
        for name,codes in CMC.items():
            matched={code:[c for c in df if c.startswith('Date_'+code+'_first_reported_')] for code in codes}
            cmc_audit['available_codes'][name]=[x for x,cols in matched.items() if cols]
            cmc_audit['absent_codes'][name]=[x for x,cols in matched.items() if not cols]
            if not any(matched.values()):raise ValueError('No available diagnosis fields for '+name)
            mask=pd.Series(False,index=df.index)
            for cols in matched.values():
                for c in cols:
                    z=pd.to_datetime(df[c],errors='coerce')
                    mask |= z.notna()&(z>pd.Timestamp('1900-01-01'))&(z<=dates[MRI])
            df[name]=mask.astype(int)
        df['CMC_score_raw']=df[list(CMC)].sum(axis=1)
        df['CMC_score_cat']=df.CMC_score_raw.clip(upper=2)
        cmc_audit['changed_count']=int((old!=df.CMC_score_cat).sum())
    for j,(label,col) in enumerate(RAW.items()):
        z=df[col]
        if (z.dropna()<0).any():raise ValueError('Negative WMH volume in '+label)
        df['wmh_'+str(j)]=np.log1p(z*df[SCALE])
    for j,(label,col) in enumerate(COG.items()):
        z=df[col].where(df[col]>=0)
        if j<2:z=-np.log(z.where(z>0))
        sd=z.std(ddof=1)
        if not sd>0:raise ValueError('No variation in cognitive endpoint '+label)
        df['cog_'+str(j)]=(z-z.mean())/sd
    return df,cmc_audit

def restricted_mask(df):
    mri=pd.to_datetime(df[MRI],errors='coerce');excluded=pd.Series(False,index=df.index)
    for code in NEURO:
        cols=[c for c in df if c.startswith('Date_'+code+'_first_reported_')]
        if not cols:raise ValueError('Missing restriction code '+code)
        for c in cols:
            z=pd.to_datetime(df[c],errors='coerce')
            excluded |= z.notna()&(z>pd.Timestamp('1900-01-01'))&(z<=mri)
    return ~excluded

def terms(covs):return ' + '.join('C('+c+')' if c in CAT else c for c in covs)

def estimate_weights(df):
    x=dmatrix(terms(CORE),df,return_type='dataframe').drop(columns='Intercept')
    y=df.treatment_var.to_numpy(dtype=int)
    xs=StandardScaler().fit_transform(x)
    # C=inf requests no penalty across recent sklearn versions without deprecated penalty=None.
    m=LogisticRegression(C=np.inf,solver='lbfgs',max_iter=5000,tol=1e-10)
    with warnings.catch_warnings(record=True) as seen:
        m.fit(xs,y)
    if any('converge' in str(w.message).lower() for w in seen):raise RuntimeError('PS did not converge.')
    ps=m.predict_proba(xs)[:,1]
    if not np.isfinite(ps).all() or ((ps<=0)|(ps>=1)).any():raise ValueError('PS violates numerical positivity.')
    w=np.where(y==1,1-ps,ps)
    rows=[]
    for j,c in enumerate(x.columns):
        a=x[c].to_numpy(dtype=float);s=np.sqrt((a[y==1].var(ddof=1)+a[y==0].var(ddof=1))/2)
        before=(a[y==1].mean()-a[y==0].mean())/s if s else 0
        after=(np.average(a[y==1],weights=w[y==1])-np.average(a[y==0],weights=w[y==0]))/s if s else 0
        rows.append({'covariate':c,'smd_before':before,'smd_after':after})
    balance=pd.DataFrame(rows)
    if balance.smd_after.abs().max()>.1:raise RuntimeError('Covariate balance gate failed.')
    ess=lambda q:float(q.sum()**2/(q@q))
    diagnostics={'n':len(df),'pd_n':int(y.sum()),'control_n':int((1-y).sum()),
                 'max_abs_smd_before':float(balance.smd_before.abs().max()),
                 'max_abs_smd_after':float(balance.smd_after.abs().max()),
                 'ess_all':ess(w),'ess_pd':ess(w[y==1]),'ess_control':ess(w[y==0]),
                 'iterations':m.n_iter_.tolist()}
    df=df.copy();df['overlap_weight']=w
    return df,balance,diagnostics

def model(df,outcome,covs=CORE,weighted=False,extra=''):
    f=outcome+' ~ treatment_var + '+terms(covs)+extra
    yy,xx=dmatrices(f,df,return_type='dataframe')
    fit=(sm.WLS(yy,xx,weights=df.loc[yy.index,'overlap_weight']) if weighted else sm.OLS(yy,xx)).fit(cov_type='HC3')
    b=fit.params['treatment_var'];lo,hi=fit.conf_int().loc['treatment_var']
    z={'beta':float(b),'ci_low':float(lo),'ci_high':float(hi),'p':float(fit.pvalues['treatment_var']),
       'n':int(fit.nobs),'pd_n':int(df.loc[yy.index,'treatment_var'].sum()),'weighted':weighted,
       'pct':float(np.expm1(b)*100),'pct_low':float(np.expm1(lo)*100),'pct_high':float(np.expm1(hi)*100)}
    return z

def primary(df,label,weighted):
    rows=[]
    for j,o in enumerate(RAW):
        z=model(df,'wmh_'+str(j),weighted=weighted);z.update(cohort=label,outcome=o);rows.append(z)
    q=multipletests([rows[1]['p'],rows[2]['p']],method='fdr_bh')[1]
    rows[0]['q']=None;rows[1]['q']=float(q[0]);rows[2]['q']=float(q[1])
    return rows

def baseline(df):
    rows=[]
    vars=CORE+list(CMC)+['e4_count']
    for c in vars:
        for t in [0,1]:
            sub=df.loc[df.treatment_var==t,c]
            w=df.loc[df.treatment_var==t,'overlap_weight']
            rows.append({'variable':c,'group':t,'n':len(sub),'mean':float(sub.mean()),'sd':float(sub.std(ddof=1)),
                         'weighted_mean':float(np.average(sub,weights=w)),
                         'counts':json.dumps({str(k):int(v) for k,v in sub.value_counts().sort_index().items()}) if c in CAT or c in CMC or c=='e4_count' else '',
                         'weighted_proportions':json.dumps({str(k):float(w[sub==k].sum()/w.sum()) for k in sorted(sub.unique())}) if c in CAT or c in CMC or c=='e4_count' else ''})
    return rows

def hypertension_sensitivity(df):
    """Exploratory adjustment for baseline *reported* hypertension, not true absence.

    Uses available Field 20002 Instance 0 arrays, coding 6 values 1065 and 1072.
    Empty illness arrays are treated as no reported hypertension. This does not
    replace complete BP/history ascertainment or repair missing ICD I10-I12.
    """
    cols=[c for c in df if c.startswith('Non_cancer_illness_code_self_reported_Instance_0_Array_')]
    if not cols:return [],{'available':False}
    date0='Date_of_attending_assessment_centre_Instance_0'
    require_columns(df,[date0])
    baseline_date=pd.to_datetime(df[date0],errors='coerce')
    mri_date=pd.to_datetime(df[MRI],errors='coerce')
    if baseline_date.isna().any() or mri_date.isna().any() or (baseline_date>mri_date).any():
        raise ValueError('Baseline self-report is not confirmed to precede MRI.')
    codes=df[cols].apply(pd.to_numeric,errors='coerce')
    d=df.copy();d['reported_hypertension']=codes.isin([1065,1072]).any(axis=1).astype(int)
    rows=[]
    for j,o in enumerate(RAW):
        z=model(d,'wmh_'+str(j),CORE+['reported_hypertension'],True)
        z.update(outcome=o,analysis='Full weights plus reported baseline hypertension');rows.append(z)
    report={'available':True,'codes':[1065,1072],'arrays_used':len(cols),
            'pd_reported_n':int(d.loc[d.treatment_var==1,'reported_hypertension'].sum()),
            'control_reported_n':int(d.loc[d.treatment_var==0,'reported_hypertension'].sum()),
            'empty_arrays_n':int(codes.isna().all(axis=1).sum()),
            'note':'No reported hypertension is not clinically confirmed absence; full-cohort weights retained.'}
    return rows,report

def path_decomposition(df,n_boot,seed,outdir):
    """Unweighted OLS path products; resample participants, no causal interpretation.

    Pair bootstrap uses multinomial row multiplicities and exact sufficient
    statistics. The same complete cases are used for a, b and c per comparison.
    Covariates and cognitive standardisation are fixed as in the observed model.
    """
    rows=[]
    for k,outcome in enumerate(COG):
        for j,mediator in enumerate(RAW):
            names=['treatment_var']+CORE+['e4_count','wmh_'+str(j),'cog_'+str(k)]
            d=df[names].dropna().copy()
            x=dmatrix('treatment_var + '+terms(CORE+['e4_count']),d,return_type='dataframe')
            # Stabilise linear algebra without changing treatment coefficients.
            for c in x:
                if c not in {'Intercept','treatment_var'} and x[c].std()>0:
                    x[c]=(x[c]-x[c].mean())/x[c].std()
            xx=x.to_numpy();m=d['wmh_'+str(j)].to_numpy();y=d['cog_'+str(k)].to_numpy()
            ai=list(x.columns).index('treatment_var')
            am=sm.OLS(m,xx).fit(cov_type='HC3')
            xy=np.column_stack([xx,m]);bm=sm.OLS(y,xy).fit(cov_type='HC3')
            cm=sm.OLS(y,xx).fit(cov_type='HC3')
            a=float(am.params[ai]);b=float(bm.params[-1]);ab=a*b
            vals=[];rng=np.random.default_rng(seed+100*k+j)
            if n_boot:
                # Crossproducts of [X,M,Y] cover both path regressions.
                full=np.column_stack([xx,m,y]);p=full.shape[1];ii,jj=np.triu_indices(p)
                features=np.ascontiguousarray(full[:,ii]*full[:,jj])
                prob=np.full(len(d),1/len(d));batch=16
                for start in range(0,n_boot,batch):
                    count=min(batch,n_boot-start)
                    weights=rng.multinomial(len(d),prob,size=count).astype(float)
                    sums=weights@features
                    for v in sums:
                        gram=np.zeros((p,p));gram[ii,jj]=v;gram[jj,ii]=v
                        q=xx.shape[1]
                        try:
                            aa=np.linalg.solve(gram[:q,:q],gram[:q,q])[ai]
                            bb=np.linalg.solve(gram[:q+1,:q+1],gram[:q+1,q+1])[-1]
                            if np.isfinite(aa*bb):vals.append(aa*bb)
                        except np.linalg.LinAlgError:pass
                    if (start+count)%256==0:print('bootstrap',outcome,mediator,start+count,flush=True)
                if len(vals)<.95*n_boot:raise RuntimeError('Too many nonestimable bootstrap draws.')
            lo,hi=np.quantile(vals,[.025,.975]) if vals else [np.nan,np.nan]
            rows.append({'outcome':outcome,'mediator':mediator,'n':len(d),'pd_n':int(d.treatment_var.sum()),
                         'a':a,'b':b,'ab':ab,'ab_boot_mean':float(np.mean(vals)) if vals else None,
                         'ci_low':float(lo),'ci_high':float(hi),'total_beta':float(cm.params[ai]),
                         'direct_beta':float(bm.params[ai]),'bootstrap_success':len(vals),'seed':seed+100*k+j})
            pd.DataFrame(rows).to_csv(outdir/'path_products.csv',index=False)
            print('path complete',outcome,mediator,len(d),len(vals),flush=True)
    return rows

def run(args):
    src=Path(args.input);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(src,low_memory=False)
    withdrawals={'provided':bool(args.withdrawn_list),'excluded_n':None}
    if args.withdrawn_list:
        require_columns(d,['Participant_ID'])
        withdrawn=pd.read_csv(args.withdrawn_list,header=None,usecols=[0],dtype=str).iloc[:,0].str.strip()
        ids=d.Participant_ID.astype(str).str.replace(r'\.0$','',regex=True).str.strip()
        keep=~ids.isin(set(withdrawn));withdrawals['excluded_n']=int((~keep).sum())
        d=d.loc[keep].copy()
    d,cmc_audit=prepare(d,repair_cmc=not args.keep_legacy_cmc)
    d,balance,diag=estimate_weights(d);balance.to_csv(out/'balance_full.csv',index=False)
    r=d.loc[restricted_mask(d)].copy();r,rb,rd=estimate_weights(r);rb.to_csv(out/'balance_restricted.csv',index=False)
    rows=primary(d,'Full overlap weighted',True)+primary(r,'Restricted overlap weighted',True)+primary(d,'Unweighted adjusted',False)
    pd.DataFrame(rows).to_csv(out/'wmh_results.csv',index=False)
    cog=[]
    for weighted in [False,True]:
        for j,o in enumerate(COG):
            z=model(d,'cog_'+str(j),CORE+['e4_count'],weighted);z.update(outcome=o);cog.append(z)
    pd.DataFrame(cog).to_csv(out/'cognitive_results.csv',index=False)
    pd.DataFrame(baseline(d)).to_csv(out/'baseline_aggregate.csv',index=False)
    vascular,vascular_audit=hypertension_sensitivity(d)
    pd.DataFrame(vascular).to_csv(out/'hypertension_sensitivity.csv',index=False)
    subgroup=[]
    duration=(pd.to_datetime(d[MRI])-pd.to_datetime(d[PDDATE],errors='coerce')).dt.days/365.25
    for label,mask in [('Male',d.Sex==1),('Female',d.Sex==0),('PD interval <=5 years',(d.treatment_var==0)|(duration<=5)),('PD interval >5 years',(d.treatment_var==0)|(duration>5))]:
        sub=d.loc[mask].copy()
        covs=[c for c in CORE if c!='Sex'] if label in ['Male','Female'] else CORE
        for j,o in enumerate(RAW):
            z=model(sub,'wmh_'+str(j),covs,True);z.update(stratum=label,outcome=o);subgroup.append(z)
    pd.DataFrame(subgroup).to_csv(out/'subgroup_results.csv',index=False)
    interaction=[]
    for j,o in enumerate(RAW):
        yy,xx=dmatrices('wmh_'+str(j)+' ~ treatment_var * C(Sex) + '+terms([c for c in CORE if c!='Sex']),d,return_type='dataframe')
        m=sm.WLS(yy,xx,weights=d.loc[yy.index,'overlap_weight']).fit(cov_type='HC3')
        cols=[c for c in xx if 'treatment_var:C(Sex)' in c]
        rr=np.eye(len(xx.columns))[[list(xx).index(c) for c in cols]]
        interaction.append({'outcome':o,'sex_interaction_p':float(m.wald_test(rr,scalar=True).pvalue)})
    pd.DataFrame(interaction).to_csv(out/'sex_interactions.csv',index=False)
    pathways=path_decomposition(d,args.bootstrap,args.seed,out) if args.paths else []
    report={'version':'1.0.0','python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
            'scipy':scipy.__version__,'statsmodels':sm.__version__,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
            'cmc_audit':cmc_audit,'full':diag,'restricted':rd,'restriction_excluded':len(d)-len(r),
            'hypertension_sensitivity':vascular_audit,
            'bootstrap_requested':args.bootstrap if args.paths else 0,'seed':args.seed,
            'participant_rows_exported':False,'withdrawals':withdrawals,
            'consent_list_currentness_verified_by_this_program':False,
            'analysis_note':'Method corrections documented on 2026-10-07; analyses are cross-sectional.'}
    (out/'run_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--withdrawn-list',help='Private headerless current withdrawal IDs; never redistribute this file.')
    p.add_argument('--paths',action='store_true');p.add_argument('--bootstrap',type=int,default=2000)
    p.add_argument('--seed',type=int,default=20261007);p.add_argument('--keep-legacy-cmc',action='store_true')
    a=p.parse_args()
    if a.bootstrap<0:p.error('bootstrap must be nonnegative')
    with threadpool_limits(limits=4):run(a)
