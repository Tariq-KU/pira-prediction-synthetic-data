"""Generate wholly synthetic RRMS landmark data. No patient records are used.

Run: python generate_data.py --out /path/to/new/output
Dependencies: Python 3, numpy, pandas. The seed is fixed; no fitted models are used.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import numpy as np
import pandas as pd

SEED = 20261022
N_DEV, N_VAL = 1000, 500
SHAPE = 1.35
TARGET_RISK = 0.15
DMT = np.array(['anti_cd20', 'other_high_efficacy', 'moderate_oral', 'injectable'])
CLINICAL = ['age_m12_years', 'sex', 'disease_duration_m12_years', 'edss_m12',
            'relapses_m0_m12', 'dmt_class_m12', 'dmt_duration_m12_months',
            'bmi_m12_kg_m2', 'egfr_m12_ml_min_1_73m2']
MRI = ['brain_parenchymal_fraction_m12', 'brain_volume_change_pct_m0_m12',
       't2_lesion_volume_m12_ml', 'cord_area_m12_mm2', 'cord_area_change_pct_m0_m12']
BLOOD = ['snfl_m12_pg_ml', 'sgfap_m12_pg_ml', 'snfl_log2_change_m0_m12',
         'sgfap_log2_change_m0_m12']
SOURCES = [
 ['S1', 'Morris et al. (2019)', 'Using simulation studies to evaluate statistical methods',
  'https://doi.org/10.1002/sim.8086', 'Simulation design, reproducibility and explicit generating assumptions; not the numerical coefficients.'],
 ['S2', 'Müller et al. (2025)', 'Standardized Definition of PIRA in Relapsing-Remitting Multiple Sclerosis',
  'https://doi.org/10.1001/jamaneurol.2025.0495', 'Intended EDSS and relapse-based clinical endpoint; this dataset does not implement full adjudication.'],
 ['S3', 'Meier et al. (2023)', 'Serum GFAP Compared With NfL as a Biomarker for Disease Progression in MS',
  'https://doi.org/10.1001/jamaneurol.2022.5250', 'Broad assay concentration scale, differing inflammatory associations and repeated measurements; not population calibration.'],
 ['S4', 'Benkert et al. (2025)', 'Serum GFAP and NfL Reflect Different Mechanisms of Disease Progression under B-Cell Depleting Treatment in MS',
  'https://doi.org/10.1002/ana.27096', 'Rationale for examining both levels and changes; does not establish the assumed incremental effects in RRMS.'],
 ['S5', 'Riley et al. (2020)', 'Calculating the sample size required for developing a clinical prediction model',
  'https://doi.org/10.1136/bmj.m441', 'Why eventual clinical recruitment needs separate justification; 1000/500 are illustrative, not calculated requirements.'],
]

def clamp(x, lo, hi):
    return np.clip(x, lo, hi)

def make_people(n, validation, seed):
    """All first-year predictors are generated before any event time."""
    r = np.random.default_rng(seed)
    nc = 3 if validation else 6
    centre = r.integers(1, nc+1, n)
    centre_z = (centre - (nc+1)/2) / nc
    age0 = clamp(r.normal(41 + 2*validation, 10, n) + 2*centre_z, 18, 70)
    female = r.random(n) < (0.69 if validation else 0.71)
    dur0 = np.minimum(clamp(r.gamma(2.2, 4.0, n), 1.5, 35), age0-12)
    bmi = clamp(r.normal(26+0.4*validation, 4.6, n), 17, 43)
    egfr = clamp(108-0.65*(age0-30)+r.normal(0,12,n), 40, 135)
    burden, inflammation = r.multivariate_normal([0,0], [[1,.30],[.30,1]], n).T
    edss0 = clamp(np.round((1.8 + .045*dur0 + .015*(age0-40) + .70*burden + r.normal(0,.6,n))*2)/2,0,6.5)
    edss0[edss0 == .5] = 1.0
    # All simulated individuals are already eligible at month 12.
    # Stable EDSS over the predictor year is an explicit simplifying assumption.
    probs = np.array([.42,.23,.23,.12]) if not validation else np.array([.48,.22,.20,.10])
    dmt_ix = r.choice(4,n,p=probs)
    dmt0 = np.minimum(12+r.gamma(2,15,n),dur0*12)
    relapse_lambda = np.exp(-1.35+.40*inflammation + np.array([-.35,-.25,0,.15])[dmt_ix])
    relapses = np.minimum(r.poisson(relapse_lambda),4)
    # Latent shared factors induce realistic-direction correlations; magnitudes are assumptions.
    icv = clamp(r.normal(1450+120*(~female),90,n),1150,1900)
    bpf0 = clamp(.815-.0008*(age0-40)-.0012*dur0-.010*burden+r.normal(0,.015,n),.68,.88)
    brain0 = icv*bpf0
    bpct = clamp(-.40-.12*burden-.004*(age0-40)-.045*relapses+r.normal(0,.17,n),-1.8,.25)
    brain12 = brain0*(1+bpct/100)
    lesion0 = clamp(np.exp(np.log(4.0)+.045*dur0+.45*burden+.40*inflammation+r.normal(0,.55,n)),.05,65)
    lesion12 = clamp(lesion0*np.exp(.03+.08*relapses+r.normal(0,.12,n)),.05,80)
    cord0 = clamp(78-.15*dur0-.07*(age0-40)-2.5*burden+r.normal(0,5,n),45,98)
    cpct = clamp(-.65-.20*burden-.008*(age0-40)+r.normal(0,.35,n),-3.0,.5)
    cord12 = cord0*(1+cpct/100)
    ln0 = np.log(8.0)+.020*(age0-40)-.020*(bmi-26)+.004*(95-egfr)+.14*burden+.22*inflammation+.09*relapses+r.normal(0,.22,n)
    lg0 = np.log(65.0)+.018*(age0-40)+.08*female-.012*(bmi-26)+.003*(95-egfr)+.24*burden+r.normal(0,.24,n)
    nslope = .015+.045*burden+r.normal(0,.085,n)
    gslope = .015+.035*burden+r.normal(0,.06,n)
    blood = []
    nf, gf = {}, {}
    for month in [0,6,12]:
        # A sporadic transient component affects NfL more strongly than GFAP.
        pulse = r.binomial(1,clamp(.025+.10*relapses,0,.35)) * r.uniform(.15,.45,n)
        nf[month] = clamp(np.exp(ln0+nslope*month/12+pulse+r.normal(0,.10,n)),2,120).round(4)
        gf[month] = clamp(np.exp(lg0+gslope*month/12+r.normal(0,.07,n)),15,600).round(4)
    prefix = 'SYN_VAL' if validation else 'SYN_DEV'
    ids = [f'{prefix}_{i+1:04d}' for i in range(n)]
    d = pd.DataFrame({
      'participant_id':ids, 'centre_id':[f'{"V" if validation else "D"}{c:02d}' for c in centre],
      'age_m12_years':np.round(age0+1,2), 'sex':np.where(female,'female','male'),
      'disease_duration_m12_years':np.round(dur0+1,2), 'edss_m12':edss0,
      'relapses_m0_m12':relapses, 'dmt_class_m12':DMT[dmt_ix],
      'dmt_duration_m12_months':np.round(dmt0+12,1), 'bmi_m12_kg_m2':np.round(bmi,2),
      'egfr_m12_ml_min_1_73m2':np.round(egfr,2),
      'brain_parenchymal_fraction_m12':np.round(brain12/icv,6),
      'brain_volume_change_pct_m0_m12':np.round(bpct,4),
      't2_lesion_volume_m12_ml':np.round(lesion12,4),
      'cord_area_m12_mm2':np.round(cord12,4), 'cord_area_change_pct_m0_m12':np.round(cpct,4),
      'snfl_m12_pg_ml':nf[12], 'sgfap_m12_pg_ml':gf[12],
      'snfl_log2_change_m0_m12':np.round(np.log2(nf[12]/nf[0]),6),
      'sgfap_log2_change_m0_m12':np.round(np.log2(gf[12]/gf[0]),6)
    })
    for m in [0,6,12]:
        blood.append(pd.DataFrame({'participant_id':ids,'study_month':m,'snfl_pg_ml':nf[m],'sgfap_pg_ml':gf[m]}))
    imaging=[]
    for m,b,l,c in [(0,brain0,lesion0,cord0),(12,brain12,lesion12,cord12)]:
        imaging.append(pd.DataFrame({'participant_id':ids,'study_month':m,
        'intracranial_volume_ml':icv.round(4),'brain_volume_ml':b.round(4),
        't2_lesion_volume_ml':l.round(4),'cord_area_mm2':c.round(4)}))
    # Information for eligibility audit, not fitted predictors or real case-report forms.
    eligibility = pd.DataFrame({'participant_id':ids,'rrms_at_enrolment':1,
      'age_m0_years':age0.round(2),'disease_duration_m0_years':dur0.round(2),
      'same_dmt_duration_at_enrolment_months':dmt0.round(1),'edss_m0':edss0,
      'edss_m6':edss0,'edss_m12':edss0,'unresolved_worsening_at_m12':0,
      'dmt_switch_m0_m12':0})
    return d,pd.concat(blood,ignore_index=True),pd.concat(imaging,ignore_index=True),eligibility

def log_hazard(d, scenario):
    """All coefficients are declared assumptions, not literature estimates."""
    dmt = d.dmt_class_m12.map(dict(zip(DMT,[0,.04,.06,.08]))).to_numpy()
    z = (.28*(d.age_m12_years-42)/10 + .12*(d.sex=='male')
         + .12*(d.disease_duration_m12_years-9)/5 + .28*(d.edss_m12-2.5)/1.5
         + .08*d.relapses_m0_m12 + dmt + .05*(d.bmi_m12_kg_m2-26)/5
         + .05*(95-d.egfr_m12_ml_min_1_73m2)/15
         + .22*(.80-d.brain_parenchymal_fraction_m12)/.035
         + .25*(-d.brain_volume_change_pct_m0_m12-.45)/.35
         + .08*np.log2((d.t2_lesion_volume_m12_ml+1)/6)
         + .17*(75-d.cord_area_m12_mm2)/8
         + .12*(-d.cord_area_change_pct_m0_m12-.7)/.6)
    if scenario == 'modest_added_signal':
        z += (.25*np.log2(d.snfl_m12_pg_ml/8) + .35*np.log2(d.sgfap_m12_pg_ml/65)
              + .10*d.snfl_log2_change_m0_m12 + .16*d.sgfap_log2_change_m0_m12)
    return np.asarray(z,float)

def reference_hazard(scenario):
    # Calibrate expected risk in a separate synthetic reference population.
    # No outcome from the delivered development or validation cohorts is used.
    ref,*_ = make_people(40000,False,SEED+700)
    hr=np.exp(log_hazard(ref,scenario))
    lo,hi=1e-6,3.0
    for _ in range(70):
        mid=(lo+hi)/2
        if np.mean(1-np.exp(-mid*hr))>TARGET_RISK: hi=mid
        else: lo=mid
    return (lo+hi)/2

def add_outcome(d,scenario,seed,h0):
    r=np.random.default_rng(seed)
    lp=log_hazard(d,scenario)
    u=r.uniform(1e-10,1-1e-10,len(d))
    continuous=24*(-np.log(u)/(h0*np.exp(lp)))**(1/SHAPE)
    detected=np.ceil(continuous/6)*6
    # Independent analytical censoring with a 12-month confirmation tail.
    early=r.random(len(d))<.12
    c=np.where(early,r.choice([6.,12.,18.],len(d),p=[.25,.35,.40]),24.)
    event=(detected<=c).astype(int)
    time=np.where(event==1,detected,c)
    out=d.copy()
    out['followup_months']=time
    out['pira_event']=event
    out['pira_by_24m']=np.where(event==1,1,np.where(c==24,0,np.nan))
    out['pira_onset_study_month']=np.where(event==1,12+detected,np.nan)
    out['pira_confirmation_study_month']=np.where(event==1,24+detected,np.nan)
    out['adjudicable_through_study_month']=12+c
    out['last_followup_study_month']=24+c
    out['outcome_status']=np.where(event==1,'confirmed_pira',np.where(c==24,'no_pira_through_month36','censored_before_month36'))
    truth=pd.DataFrame({'participant_id':d.participant_id,'scenario':scenario,
      'generating_log_hazard':lp,'true_risk_by_24m':1-np.exp(-h0*np.exp(lp)),
      'latent_continuous_time_months':continuous,'true_detection_time_months':detected,
      'true_pira_by_24m':(detected<=24).astype(int),'analytical_censor_months':c})
    return out,truth

def mask_predictors(d,blood,mri,seed):
    """MCAR missingness, independent of the event outcome. No imputation here."""
    r=np.random.default_rng(seed)
    observed=d.copy(); b=blood.copy(); m=mri.copy()
    b['snfl_unavailable']=0; b['sgfap_unavailable']=0
    for col in ['snfl','sgfap']:
        miss=r.random(len(b))<.025
        b.loc[miss,col+'_pg_ml']=np.nan
        b.loc[miss,col+'_unavailable']=1
    m['scan_unusable']=0
    unusable=(m.study_month==12)&(r.random(len(m))<.04)
    m.loc[unusable,['brain_volume_ml','t2_lesion_volume_ml','cord_area_mm2']]=np.nan
    m.loc[unusable,'scan_unusable']=1
    for col in ['snfl','sgfap']:
        p=b.pivot(index='participant_id',columns='study_month',values=col+'_pg_ml').reindex(d.participant_id)
        observed[col+'_m12_pg_ml']=p[12].to_numpy()
        observed[col+'_log2_change_m0_m12']=np.round(np.log2(p[12]/p[0]).to_numpy(),6)
    ids=set(m.loc[unusable,'participant_id'])
    observed.loc[observed.participant_id.isin(ids),MRI]=np.nan
    return observed,b,m

def dictionary():
    rows=[]
    def add(name,role,unit,description,transform='None'):
        rows.append([name,role,unit,description,transform])
    add('participant_id','Identifier; exclude from predictors','text','Wholly fictional participant identifier. Never join scenarios as separate people.')
    add('centre_id','Grouping; exclude from primary predictors','text','Fictional development D01-D06 or validation V01-V03 centre. No overlapping centres.')
    desc={
      'age_m12_years':('years','Age at the month-12 prediction landmark.'),
      'sex':('category','female/male, a simplified recorded-sex variable; not gender identity.'),
      'disease_duration_m12_years':('years','Time since first clinical demyelinating event at month 12.'),
      'edss_m12':('EDSS points','Month-12 EDSS; valid half-point values except 0.5; restricted to 0-6.5.'),
      'relapses_m0_m12':('count','Count of first-year relapses. No relapse dates or post-landmark histories are simulated.'),
      'dmt_class_m12':('category','anti_cd20 / other_high_efficacy / moderate_oral / injectable. Broad teaching categories, not complete drug taxonomy.'),
      'dmt_duration_m12_months':('months','Duration of the same DMT at month 12. At least 24 months, including >=12 months before enrolment.'),
      'bmi_m12_kg_m2':('kg/m2','Body mass index at month 12.'),
      'egfr_m12_ml_min_1_73m2':('mL/min/1.73m2','Estimated kidney function; included equally in all four models.'),
      'brain_parenchymal_fraction_m12':('fraction 0-1','Brain volume divided by intracranial volume at month 12.'),
      'brain_volume_change_pct_m0_m12':('percentage points','100*(brain volume month12 / baseline - 1). Negative values indicate volume loss. Numeric proxy, not an actual SIENA output.'),
      't2_lesion_volume_m12_ml':('mL','T2 lesion volume at month 12; positive and right-skewed.'),
      'cord_area_m12_mm2':('mm2','Cervical cord cross-sectional area at month 12. Numeric proxy, not actual image segmentation.'),
      'cord_area_change_pct_m0_m12':('percentage points','100*(cord area month12 / baseline - 1). Negative values indicate area loss.'),
      'snfl_m12_pg_ml':('pg/mL','Month-12 serum neurofilament light concentration.'),
      'sgfap_m12_pg_ml':('pg/mL','Month-12 serum GFAP concentration.'),
      'snfl_log2_change_m0_m12':('log2 ratio','log2(month12 / baseline). 0=no change; 1=doubling; -1=halving.'),
      'sgfap_log2_change_m0_m12':('log2 ratio','log2(month12 / baseline). 0=no change; 1=doubling; -1=halving.')}
    for name in CLINICAL+MRI+BLOOD:
        role='Clinical predictor' if name in CLINICAL else ('MRI predictor' if name in MRI else 'Blood predictor')
        transform='Categorical encoding' if name in ['sex','dmt_class_m12'] else 'Scale using training data only'
        if name in ['snfl_m12_pg_ml','sgfap_m12_pg_ml']:transform='log2 concentration, then scale within training data'
        if name=='t2_lesion_volume_m12_ml':transform='log2(1 + volume), then scale within training data'
        add(name,role,*desc[name],transform)
    add('followup_months','Survival outcome; exclude from predictors','months after month12','Time to confirmed-event onset, otherwise effective censoring time. Use with pira_event, not alone.')
    add('pira_event','Survival outcome; exclude from predictors','0/1','1=qualifying event observed and confirmed. 0=no observed confirmed event before censoring; may have unknown 24-month status.')
    add('pira_by_24m','Horizon outcome; exclude from predictors','0/1/blank','1=event by month36; 0=known no event through month36; blank=unknown horizon status. Never replace blanks with 0.')
    add('pira_onset_study_month','Outcome audit; exclude from predictors','study month','18/24/30/36 for observed events. Blank otherwise. First qualifying detection, not biological onset.')
    add('pira_confirmation_study_month','Outcome audit; exclude from predictors','study month','Observed onset +12 months. Confirmation is assumed successful; no full EDSS/relapse adjudication is simulated.')
    add('adjudicable_through_study_month','Outcome audit; exclude from predictors','study month','Latest candidate onset date with a full 12-month confirmation tail available; last_followup_study_month minus12.')
    add('last_followup_study_month','Outcome audit; exclude from predictors','study month','30/36/42/48; the last available scheduled clinical visit. This is not the Cox time variable.')
    add('outcome_status','Outcome audit; exclude from predictors','category','confirmed_pira / no_pira_through_month36 / censored_before_month36.')
    return pd.DataFrame(rows,columns=['variable','role','unit_or_coding','definition','planned_preprocessing'])

def assumptions():
    return pd.DataFrame([
      ['Purpose','One synthetic learning dataset plus a no-added-blood-signal companion','Design choice','No real participants, patient-level source records, lab tests or scans. S1'],
      ['Cohorts','1000 development; 500 validation; six D centres and three V centres','Assumption','Independent draws, not a random split of one delivered cohort. Not a Riley sample-size result. S5'],
      ['Population','Adults with RRMS; stable treatment >=12 months before entry; eligible at month12','Protocol alignment / simplification','No progressive MS at entry. No first-year DMT switches or unresolved EDSS worsening are simulated.'],
      ['Clinical profiles','Age at enrolment: clipped normal mean41 SD10 (validation mean43), range18-70; ~71% female (~69% validation)','Assumption','Illustrative treated RRMS case mix; not calibrated to a country or one published cohort.'],
      ['Measurement schedule','Blood: months0,6,12; MRI: months0,12; outcome detection on 6-month visit grid','Protocol alignment','MRI/blood are numerical simulations. No raw image or biological sample exists.'],
      ['Blood distributions','Positive correlated log-normal measurements; reference levels8 pg/mL NfL and65 pg/mL GFAP','Literature-informed scale; assumed distribution','Meier reports remission/stable/active values on similar scales; no diagnostic or normative thresholds are implied. S3'],
      ['Blood associations','Age, BMI, kidney function and shared latent disease factors influence concentrations; stronger transient component for NfL','Qualitatively informed; coefficients assumed','Exact joint distribution, slopes, correlations and error levels were not extracted from patient records. S3-S4'],
      ['MRI','Correlated brain fraction, lesion volume and cervical cord area; small noisy annual volume/area changes','Assumption','Generated numerical proxies. SIENA and Spinal Cord Toolbox were not run.'],
      ['Expected outcome frequency','15% cumulative qualifying onset risk over24 months in a separate development-like reference population','Assumption','Expected risk before censoring; realised event counts can differ. Not an estimate of RRMS PIRA incidence.'],
      ['Outcome mechanism','Weibull PH event times, shape1.35; detection rounded up to the next6-month assessment','Assumption','Compatible with Cox as an initial teaching case. It cannot establish Cox superiority over other architectures.'],
      ['Added signal scenario','Blood log-hazard coefficients0.25 NfL level,0.35 GFAP level,0.10 NfL log2 change,0.16 GFAP log2 change','Assumption','Level coefficients apply per doubling. They are not reported treatment or biomarker effect estimates.'],
      ['Null companion','All four conditional blood coefficients set to zero; clinical/MRI effects and predictors retained','Assumption','No extra blood information given complete clinical/MRI predictors. With missing clinical/MRI values, blood may still act as a proxy.'],
      ['Confirmation','Event represents an already-qualified PIRA endpoint with a12-month persistence tail','Simplification','Full EDSS trajectories, relapse-exclusion logic and confirmation failures are not generated. S2'],
      ['Censoring','12% chance of an early follow-up end; effective onset censoring at6,12,18 months with probabilities0.25,0.35,0.40','Assumption','Independent of predictors and event times. Physical follow-up extends another12 months to support confirmation. No very early dropouts or death.'],
      ['Missing measurements','2.5% per biomarker measurement;4% probability of unusable month12 MRI','Assumption','MCAR: random missingness unrelated to values or outcomes. All selected MRI summaries missing when month12 scan unusable.'],
      ['Validation differences','Slightly older case mix, different DMT proportions and independent fictional centres','Assumption','Same biological/measurement generating rules. This is synthetic validation, not external clinical validation.'],
      ['Scope limits','No first-year treatment changes, death, SPMS transition, T25FW/9HPT, informative dropout or alternate clinical PIRA definitions','Initial scope','Those protocol features need later explicit simulation extensions. Do not claim they were evaluated.'],
      ['Reproducibility','Master seed20261022; documented child seeds; exact generating code supplied','Implementation','No model was fitted or dataset selected to obtain favourable performance. One realised dataset per scenario is not a repeated simulation study.'],
    ],columns=['topic','specification','basis','interpretation'])

def validate_data(obs,complete,blood,mri,elig):
    assert obs.participant_id.is_unique
    assert (elig.age_m0_years>=18).all()
    assert (elig.same_dmt_duration_at_enrolment_months>=12).all()
    assert (obs.dmt_duration_m12_months>=24).all()
    assert (obs.disease_duration_m12_years<=obs.age_m12_years-12+.02).all()
    assert not (obs.edss_m12==.5).any()
    assert set(obs.followup_months).issubset({6.,12.,18.,24.})
    ev=obs.pira_event==1
    assert np.allclose(obs.loc[ev,'pira_onset_study_month']-12,obs.loc[ev,'followup_months'])
    assert np.allclose(obs.loc[ev,'pira_confirmation_study_month']-obs.loc[ev,'pira_onset_study_month'],12)
    assert (obs.loc[ev,'pira_confirmation_study_month']<=obs.loc[ev,'last_followup_study_month']).all()
    assert (obs.loc[~ev,'pira_onset_study_month'].isna()).all()
    assert obs.pira_by_24m.isna().equals((obs.pira_event==0)&(obs.adjudicable_through_study_month<36))
    assert not blood.duplicated(['participant_id','study_month']).any()
    assert not mri.duplicated(['participant_id','study_month']).any()
    assert len(blood)==3*len(obs) and len(mri)==2*len(obs)
    assert set(obs.columns)==set(dictionary().variable)
    for marker in ['snfl','sgfap']:
        p=blood.pivot(index='participant_id',columns='study_month',values=marker+'_pg_ml').reindex(obs.participant_id)
        assert np.allclose(obs[marker+'_m12_pg_ml'],p[12],equal_nan=True)
        assert np.allclose(obs[marker+'_log2_change_m0_m12'],np.log2(p[12]/p[0]),atol=6e-7,equal_nan=True)
    brain=mri.pivot(index='participant_id',columns='study_month',values='brain_volume_ml').reindex(obs.participant_id)
    cord=mri.pivot(index='participant_id',columns='study_month',values='cord_area_mm2').reindex(obs.participant_id)
    assert np.allclose(obs.brain_volume_change_pct_m0_m12,100*(brain[12]/brain[0]-1),atol=6e-5,equal_nan=True)
    assert np.allclose(obs.cord_area_change_pct_m0_m12,100*(cord[12]/cord[0]-1),atol=3e-4,equal_nan=True)
    return True

def write_csv(df,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(path,index=False,na_rep='',float_format='%.8g')

def generate(out):
    out=Path(out); (out/'data').mkdir(parents=True,exist_ok=True)
    (out/'documentation').mkdir(exist_ok=True); (out/'reproducibility').mkdir(exist_ok=True)
    all_summary=[]; all_b=[]; all_m=[]; all_e=[]; truths=[]; all_complete=[]; work={}
    hs={s:reference_hazard(s) for s in ['modest_added_signal','no_added_signal']}
    for k,(cohort,n,isval) in enumerate([('development',N_DEV,False),('validation',N_VAL,True)]):
        base,b,m,e=make_people(n,isval,SEED+100+k)
        masked,bobs,mobs=mask_predictors(base,b,m,SEED+300+k)
        all_b.append(bobs); all_m.append(mobs); all_e.append(e)
        for s in hs:
            full,tr=add_outcome(base,s,SEED+200+k,hs[s])
            obs=full.copy(); obs[CLINICAL+MRI+BLOOD]=masked[CLINICAL+MRI+BLOOD]
            assert validate_data(obs,full,bobs,mobs,e)
            suffix='' if s=='modest_added_signal' else 'null_'
            write_csv(obs,out/'data'/f'{suffix}{cohort}.csv')
            full['scenario']=s; full['cohort']=cohort; all_complete.append(full)
            tr['cohort']=cohort; truths.append(tr)
            all_summary.append({'scenario':s,'cohort':cohort,'n':n,
              'confirmed_events':int(obs.pira_event.sum()),
              'known_non_events_by24':int((obs.pira_by_24m==0).sum()),
              'unknown_24m_outcomes':int(obs.pira_by_24m.isna().sum()),
              'early_followup_end':int((obs.last_followup_study_month<48).sum()),
              'any_missing_predictor':int(obs[CLINICAL+MRI+BLOOD].isna().any(axis=1).sum()),
              'event_fraction_observed':float(obs.pira_event.mean()),
              'median_age_m12':float(obs.age_m12_years.median()),
              'median_edss':float(obs.edss_m12.median()),
              'median_snfl':float(obs.snfl_m12_pg_ml.median()),
              'median_sgfap':float(obs.sgfap_m12_pg_ml.median())})
            if not suffix:work[cohort]=obs
    assert not set(work['development'].participant_id)&set(work['validation'].participant_id)
    assert not set(work['development'].centre_id)&set(work['validation'].centre_id)
    # Pure determinism check, without choosing a favourable seed.
    assert make_people(20,False,SEED+100)[0].equals(make_people(20,False,SEED+100)[0])
    dd=dictionary(); aa=assumptions()
    write_csv(pd.concat(all_b,ignore_index=True),out/'data'/'blood_measurements_long.csv')
    write_csv(pd.concat(all_m,ignore_index=True),out/'data'/'mri_measurements_long.csv')
    write_csv(pd.concat(all_e,ignore_index=True),out/'data'/'eligibility_context.csv')
    write_csv(dd,out/'documentation'/'data_dictionary.csv')
    write_csv(aa,out/'documentation'/'assumptions.csv')
    write_csv(pd.DataFrame(SOURCES,columns=['id','citation','title','url','role']),out/'documentation'/'sources.csv')
    write_csv(pd.DataFrame(all_summary),out/'documentation'/'dataset_summary.csv')
    # These files are supplied for audit only and are never model input.
    audit=out/'reproducibility'/'DO_NOT_USE_AS_PREDICTORS'
    write_csv(pd.concat(truths,ignore_index=True),audit/'generating_truth.csv')
    write_csv(pd.concat(all_complete,ignore_index=True),audit/'complete_values_before_missingness.csv')
    config={'master_seed':SEED,'n_development':N_DEV,'n_validation':N_VAL,
      'expected_reference_risk':TARGET_RISK,'weibull_shape':SHAPE,'baseline_cumulative_hazard_at24':hs,
      'seed_offsets':{'development_predictors':100,'validation_predictors':101,'development_outcome_and_censor':200,'validation_outcome_and_censor':201,'development_missingness':300,'validation_missingness':301,'independent_reference_population':700},
      'clinical_predictors':CLINICAL,'mri_predictors':MRI,'blood_predictors':BLOOD,
      'primary_comparison':['clinical_plus_mri','clinical_plus_mri_plus_blood'],
      'clinical_parameter_count':11,'mri_parameter_count':5,'blood_parameter_count':4,
      'total_candidate_parameters_excluding_baseline_hazard':20,
      'fit_model_now':False,'python_version':platform.python_version(),
      'numpy_version':np.__version__,'pandas_version':pd.__version__,
      'missing_value_encoding':'Blank; never zero',
      'analysis_outcome':['followup_months','pira_event']}
    (out/'reproducibility'/'simulation_config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'documentation'/'quality_checks.json').write_text(json.dumps({'status':'passed','checks':['unique identifiers','non-overlapping cohorts and centres','age and DMT eligibility','valid EDSS levels','event/censoring time bounds','12-month confirmation timing','unknown horizon outcomes preserved','biomarker longitudinal summaries reconcile','MRI change summaries reconcile within export rounding','complete dictionary coverage','deterministic seed']},indent=2)+'\n')
    payload={'development':work['development'].to_dict('records'),'validation':work['validation'].to_dict('records'),
      'dictionary':dd.to_dict('records'),'assumptions':aa.to_dict('records'),'sources':SOURCES,'summary':all_summary}
    # JSON nulls preserve blank values for the workbook builder.
    (out/'reproducibility'/'workbook_payload.json').write_text(pd.Series(payload).to_json(indent=2))
    manifest={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*.csv'))}
    (out/'reproducibility'/'csv_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'summary':all_summary,'checks':'passed','hazards':hs},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--out',type=Path,required=True)
    generate(parser.parse_args().out)
