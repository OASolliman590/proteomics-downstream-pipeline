"""Independent arithmetic and archive checks; does not refit limma or run enrichment."""
import argparse
import hashlib
import itertools
import json
import pathlib
import posixpath
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--original', type=pathlib.Path, default=pathlib.Path(__file__).parent/'original')
ap.add_argument('--output', type=pathlib.Path, default=pathlib.Path(__file__).parent/'data_checks.json')
args = ap.parse_args()
root = args.original.resolve()
a = root/'01_METHOD_A_SOURCE_PRECOMPUTED_LIMMA_INTAKE'
b = root/'02_METHOD_B_ABUNDANCE_REFIT_FULL_PIPELINE'/'01_FULL_VERIFIED_EXPORT'
def table(path):
    return pd.read_csv(path, sep='\t' if path.suffix == '.tsv' else ',')
def locate(name, base=b):
    paths = list(base.rglob(name))
    if len(paths)!=1: raise ValueError((name, paths))
    return paths[0]
def bh(values):
    p = np.array(values, dtype=float)
    valid = np.isfinite(p)
    order = np.argsort(p[valid])
    # Match R p.adjust's default n, evaluated after removing NA values.
    q = np.minimum.accumulate((p[valid][order]*valid.sum()/np.arange(1,valid.sum()+1))[::-1])[::-1]
    result = np.full(len(p), np.nan)
    target = np.flatnonzero(valid)[order]
    result[target] = np.minimum(q,1)
    return result
report = {'scope':'Independent archive, workbook, arithmetic and stored-result checks. No new limma or fgsea execution.'}
report['manifests'] = []
for manifest in [root/'00_INDEX/MANIFEST_SHA256.csv',b/'00_INDEX/MANIFEST_SHA256.csv']:
    rows = table(manifest)
    missing=[]; mismatch=[]
    base=manifest.parent.parent
    for row in rows.to_dict('records'):
        p = base/row['relative_path']
        if not p.is_file(): missing.append(row['relative_path']); continue
        if p.stat().st_size!=row['size_bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            mismatch.append(row['relative_path'])
    report['manifests'].append({'path':str(manifest.relative_to(root)), 'entries':len(rows),'missing':missing,'mismatch':mismatch})

meta = table(locate('sample_meta.csv', b/'01_INPUT_DATA'))
matrix = table(locate('abundance_input_log2.csv')).set_index('Accession')
analysis = table(locate('abundance_analysis_log2.csv')).set_index('Accession')
annot = table(locate('protein_annot.csv'))
ids = meta.sample_id.tolist()
x = matrix[ids]
coverage = {g:x[meta.loc[meta.group==g,'sample_id']].notna().sum(axis=1) for g in meta.group.unique()}
joint_keep = pd.DataFrame(coverage).max(axis=1)>=3
axis_keep = {axis:pd.DataFrame(coverage)[['CTRL',axis,axis+'_Dapa']].min(axis=1)>=3 for axis in ['PreDM','DM']}
amat = table(locate('abundance_input_log2.tsv',a)).set_index('Accession')
report['input']={'proteins':len(x),'samples':len(ids),'group_counts':meta.group.value_counts().to_dict(),'missing_cells':int(x.isna().sum().sum()),'missing_pct':float(x.isna().mean().mean()*100),'sample_missing_pct':(x.isna().mean()*100).to_dict(),'range':[float(x.min().min()),float(x.max().max())],'duplicate_accessions':int(matrix.index.duplicated().sum()),'duplicate_samples':int(meta.sample_id.duplicated().sum()),'duplicate_rat_ids':int(meta.rat_id.duplicated().sum()),'duplicate_annotation_accessions':int(annot.Accession.duplicated().sum()),'blank_gene_symbols':int(annot.gene_symbol.isna().sum()),'repeated_gene_rows':int(annot.gene_symbol.dropna().duplicated().sum()),'joint_retained':int(joint_keep.sum()),'axis_retained':{k:int(v.sum()) for k,v in axis_keep.items()},'method_A_B_matrix_equal':bool(np.allclose(x,amat.reindex(index=x.index,columns=ids),equal_nan=True)),'filtered_matrix_equal':bool(np.allclose(analysis.reindex(index=x.index[joint_keep],columns=ids), x.loc[joint_keep],equal_nan=True)),'metadata_columns':meta.columns.tolist()}

# Independently parse OOXML cells, ignoring drawings but checking workbook relationships.
book=locate('limma_results.xlsx',a)
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(book) as z:
    strings=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))]
    rels={r.attrib['Id']:posixpath.normpath('xl/'+r.attrib['Target']) for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    sheets={}
    for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets',ns):
        relid=s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        raw=[]
        for row in ET.fromstring(z.read(rels[relid])).findall('.//m:sheetData/m:row',ns):
            vals={}
            for c in row:
                ref=''.join(t for t in c.attrib['r'] if t.isalpha())
                v=c.find('m:v',ns)
                if v is None: vals[ref]=None
                elif c.attrib.get('t')=='s': vals[ref]=strings[int(v.text)]
                else:
                    try: vals[ref]=float(v.text)
                    except (TypeError,ValueError): vals[ref]=v.text
            raw.append(vals)
        header=raw[0]
        sheets[s.attrib['name']]=pd.DataFrame([{header[k]:r.get(k) for k in header} for r in raw[1:]])
    missing_relationships=[]
    for name in z.namelist():
        if not name.endswith('.rels'):continue
        for rel in ET.fromstring(z.read(name)):
            if rel.attrib.get('TargetMode')=='External':continue
            parent=posixpath.dirname(posixpath.dirname(name))
            target=posixpath.normpath(posixpath.join(parent,rel.attrib['Target'])).lstrip('/')
            if target not in z.namelist():missing_relationships.append({'rels':name,'target':target})
report['workbook']={'sha256':hashlib.sha256(book.read_bytes()).hexdigest(),'sheets':{},'missing_relationships':missing_relationships}
for name, df in sheets.items():
    d=df.set_index('Accession')
    discrepancies=[]
    for m in meta.to_dict('records'):
        if m['source_sample_label'] in d:
            if not np.allclose(d.loc[x.index,m['source_sample_label']].astype(float),x[m['sample_id']],equal_nan=True):discrepancies.append(m['sample_id'])
    report['workbook']['sheets'][name]={'rows':len(df),'sample_mismatch':discrepancies,'bh_max_abs_error':float(np.max(np.abs(bh(df.pvalue)-df['p.adj']))),'min_fdr':float(df['p.adj'].min())}

groups={g:x[meta.loc[meta.group==g,'sample_id']].mean(axis=1) for g in meta.group.unique()}
means={
'PreDM_disease_vs_CTRL':groups['PreDM']-groups['CTRL'],'DM_disease_vs_CTRL':groups['DM']-groups['CTRL'],
'PreDM_Dapa_vs_PreDM':groups['PreDM_Dapa']-groups['PreDM'],'DM_Dapa_vs_DM':groups['DM_Dapa']-groups['DM'],
'PreDM_Dapa_vs_CTRL':groups['PreDM_Dapa']-groups['CTRL'],'DM_Dapa_vs_CTRL':groups['DM_Dapa']-groups['CTRL']}
mapping={'D1_PreDM_vs_CTRL':'PreDM_disease_vs_CTRL','D2_DM_vs_CTRL':'DM_disease_vs_CTRL','T1_PreDM_Dapa_vs_PreDM':'PreDM_Dapa_vs_PreDM','T2_DM_Dapa_vs_DM':'DM_Dapa_vs_DM','R1_PreDM_Dapa_vs_CTRL':'PreDM_Dapa_vs_CTRL','R2_DM_Dapa_vs_CTRL':'DM_Dapa_vs_CTRL'}
means['S1_DM_vs_PreDM']=groups['DM']-groups['PreDM']
means['S2_DM_Dapa_vs_PreDM_Dapa']=groups['DM_Dapa']-groups['PreDM_Dapa']
means['I1_stage_by_dapa']=means['DM_Dapa_vs_DM']-means['PreDM_Dapa_vs_PreDM']
axis_de=table(locate('axis_dea_all.csv'))
joint_files=list((b/'04_JOINT_MODEL_SENSITIVITY'/'01_differential_abundance').glob('dea_*.csv'))
de_tables=[('axis',axis_de)]+[('joint',table(p)) for p in joint_files if p.name not in ['dea_master.csv','dea_summary.csv']]
report['dea']=[]
for family,df in de_tables:
    if 'P.Value' not in df:continue
    for cn,d in df.groupby('contrast'):
        expected=means[mapping.get(cn,cn)].reindex(d.Accession).to_numpy()
        q=bh(d['P.Value'])
        report['dea'].append({'family':family,'contrast':cn,'n':len(d),'logfc_max_abs_error':float(np.nanmax(abs(expected-d.logFC.to_numpy()))),'bh_max_abs_error':float(np.nanmax(abs(q-d['adj.P.Val'].to_numpy()))),'raw_p_lt_005':int((d['P.Value']<.05).sum()),'nominal_effect':int(((d['P.Value']<.05)&(abs(d.logFC)>=.58)).sum()),'fdr_lt_005':int((d['adj.P.Val']<.05).sum()),'fdr_lt_010':int((d['adj.P.Val']<.10).sum()),'min_fdr':float(d['adj.P.Val'].min()),'flag_mismatch':int((d.formal_DEP!=((d['adj.P.Val']<.05)&(abs(d.logFC)>=.25))).sum())})

rev=table(locate('axis_reversal_all.csv'))
report['reversal']=[]
for axis,d in rev.groupby('axis'):
    c=d[d.candidate_disease_protein]
    report['reversal'].append({'axis':axis,'candidates':len(c),'classes':c.reversal_class.value_counts().to_dict(),'formal_rescues':int(d.formal_statistical_rescue.sum()),'ri_error':float(np.nanmax(abs(d.reversal_index+d.treatment_logFC/d.disease_logFC))),'full_with_RI_gt_1_2':int(((c.reversal_class=='Full reversal')&(c.reversal_index>1.2)).sum()),'full_with_RI_gt_2':int(((c.reversal_class=='Full reversal')&(c.reversal_index>2)).sum()),'max_candidate_RI':float(c.reversal_index.max())})
scores=table(locate('axis_disease_axis_sample_scores.csv'))
tests=table(locate('axis_disease_axis_treatment_tests.csv'))
report['permutation']=[]
for axis,d in scores.groupby('axis'):
    u=d[d.group==axis].disease_axis_score.to_numpy(); t=d[d.group==axis+'_Dapa'].disease_axis_score.to_numpy()
    pool=np.r_[u,t]; obs=t.mean()-u.mean(); diffs=[]
    for idx in itertools.combinations(range(len(pool)),len(u)):
        mask=np.zeros(len(pool),dtype=bool);mask[list(idx)]=True
        diffs.append(pool[~mask].mean()-pool[mask].mean())
    k=int(np.sum(np.abs(diffs)>=abs(obs)-1e-12)); n=len(diffs)
    report['permutation'].append({'axis':axis,'assignments':n,'extreme':k,'stored_p':float(tests.loc[tests.axis==axis,'exact_permutation_p'].iloc[0]),'exhaustive_fixed_score_p':k/n,'plus_one_p':(k+1)/(n+1),'caveat':'These are conditional fixed-score calculations, not selection-adjusted valid treatment tests.'})

report['pathways']=[]
for family,name in [('axis','axis_gsea_all.csv'),('joint','gsea_all.csv')]:
    df=table(locate(name));globq=bh(df.pval)
    for (cn,coll),d in df.groupby(['contrast','collection']):
        report['pathways'].append({'family':family,'contrast':cn,'collection':coll,'tested':len(d),'nonfinite_p':int((~np.isfinite(d.pval)).sum()),'fdr_lt_005':int((d.padj<.05).sum()),'fdr_lt_010':int((d.padj<.1).sum()),'bh_within_family_error':float(np.nanmax(abs(bh(d.pval)-d.padj))),'fdr_lt_005_if_all_contrasts_collections_pooled':int((globq[d.index]<.05).sum())})
    report[family+'_pathway_totals']={'rows':len(df),'significant_within_family':int((df.padj<.05).sum()),'significant_pooled_diagnostic':int((globq<.05).sum())}
report['qc_watchlist']=table(locate('sample_qc.csv')).query('technical_watch').sample_id.tolist()
report['native_rds_present']=[str(p.relative_to(root)) for p in root.rglob('*.rds')]
args.output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ['pathways','workbook']},indent=2))
