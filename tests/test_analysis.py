"""Synthetic validation only; no study participant data or result fixtures."""
import sys,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pd_wmh as p

def synthetic(n=300):
    rng=np.random.default_rng(912)
    d=pd.DataFrame({'Participant_ID':np.arange(n)+900000,'treatment_var':np.r_[np.ones(60),np.zeros(n-60)],
      'Sex':rng.binomial(1,.5,n),'Age_at_Instance_2':rng.normal(65,7,n),
      'Townsend_deprivation_index_at_recruitment':rng.normal(0,2,n),'Body_mass_index_BMI_Instance_0':rng.normal(26,3,n),
      'Genetic_ethnic_grouping':rng.binomial(1,.7,n),'Smoking_Ever':rng.binomial(1,.3,n),
      'Alcohol_intake_frequency_ordinal':rng.integers(0,6,n),'has_degree':rng.binomial(1,.5,n),
      'CMC_score_cat':rng.integers(0,3,n),'e4_count':rng.binomial(2,.15,n),p.SCALE:np.ones(n),
      p.MRI:'2020-01-01',p.PDDATE:np.where(np.arange(n)<60,'2017-01-01',None),
      p.ALLDATE:np.where(np.arange(n)<60,'2017-01-01',None)})
    d.loc[d.treatment_var==1,'Age_at_Instance_2']+=4
    for codes in p.CMC.values():
        for code in codes:d['Date_'+code+'_first_reported_synthetic']=np.where(rng.random(n)<.04,'2015-01-01',None)
    for code in p.NEURO:
        if not any(x.startswith('Date_'+code+'_first_reported_') for x in d):d['Date_'+code+'_first_reported_synthetic']=None
    for col in p.RAW.values():d[col]=rng.lognormal(5,.8,n)
    for j,col in enumerate(p.COG.values()):d[col]=rng.uniform(20,80,n) if j<2 else rng.integers(2,12,n)
    return d

class AnalysisValidation(unittest.TestCase):
    def setUp(self):self.d=synthetic()
    def test_date_anchor_and_unknown_date_rejected(self):
        for date in ['2022-01-01','1900-01-01']:
            d=self.d.copy();d.loc[0,p.PDDATE]=date
            with self.assertRaisesRegex(ValueError,'date-anchored'):p.prepare(d)
    def test_schema_duplicates_and_scaling_errors(self):
        with self.assertRaisesRegex(ValueError,'Missing required'):p.prepare(self.d.drop(columns=p.SCALE))
        d=self.d.copy();d.loc[1,'Participant_ID']=d.loc[0,'Participant_ID']
        with self.assertRaisesRegex(ValueError,'not unique'):p.prepare(d)
        d=self.d.copy();d.loc[0,p.SCALE]=0
        with self.assertRaisesRegex(ValueError,'scaling'):p.prepare(d)
    def test_prefix_alias_recovers_diabetes(self):
        d=self.d.copy()
        for col in [x for x in d if x.startswith('Date_E1')]:d[col]=None
        d['Date_E11_first_reported_unexpected_label']=None;d.loc[0,'Date_E11_first_reported_unexpected_label']='2016-01-01'
        out,report=p.prepare(d)
        self.assertEqual(out.loc[0,'CMC_diabetes'],1);self.assertIn('E11',report['available_codes']['CMC_diabetes'])
        d.loc[0,'Date_E11_first_reported_unexpected_label']='2021-01-01'
        out,_=p.prepare(d);self.assertEqual(out.loc[0,'CMC_diabetes'],0)
    def test_unpenalized_overlap_balance(self):
        d,_=p.prepare(self.d)
        with threadpool_limits(4):out,b,report=p.estimate_weights(d)
        self.assertLess(b.smd_after.abs().max(),1e-4)
        self.assertAlmostEqual(out.loc[out.treatment_var==1,'overlap_weight'].sum(),out.loc[out.treatment_var==0,'overlap_weight'].sum(),places=4)
        self.assertLessEqual(report['ess_pd'],60.000001)
    def test_cognitive_direction_and_no_outcome_imputation(self):
        d=self.d.copy();col=list(p.COG.values())[0];d.loc[0,col]=25;d.loc[1,col]=75;d.loc[2,col]=np.nan
        out,_=p.prepare(d)
        self.assertGreater(out.loc[0,'cog_0'],out.loc[1,'cog_0']);self.assertTrue(pd.isna(out.loc[2,'cog_0']))
    def test_missing_covariate_and_missing_component_fail(self):
        d=self.d.copy();d['Sex']=np.nan
        with self.assertRaisesRegex(ValueError,'Entirely missing'):p.prepare(d)
        d=self.d.drop(columns=[x for x in self.d if x.startswith('Date_E78_')])
        with self.assertRaisesRegex(ValueError,'No available diagnosis'):p.prepare(d)
    def test_path_identity_bootstrap_and_no_participant_export(self):
        d,_=p.prepare(self.d)
        test_root=Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(dir=test_root) as td,threadpool_limits(4):
            self.assertTrue(Path(td).resolve().is_relative_to(test_root))
            rows=p.path_decomposition(d,64,11,Path(td))
            self.assertEqual(len(rows),12)
            for row in rows:
                self.assertAlmostEqual(row['total_beta'],row['direct_beta']+row['ab'],places=9)
                self.assertGreaterEqual(row['bootstrap_success'],61)
                self.assertLessEqual(row['ci_low'],row['ci_high'])
            f=Path(td)/'path_products.csv';header=f.read_text().splitlines()[0]
            self.assertNotIn('Participant_ID',header);self.assertNotIn('Date_',header)
            self.assertEqual({x.name for x in Path(td).iterdir()},{'path_products.csv'})
if __name__=='__main__':unittest.main()
