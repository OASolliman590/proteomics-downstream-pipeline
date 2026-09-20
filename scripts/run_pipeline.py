"""Portable preparation and execution of the recovered historical analysis.

Only preparation/preflight were validated during recovery. Statistical limitations
are retained and documented in docs/audit/SCIENTIFIC_AUDIT.md.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
GROUPS = {'CTRL', 'PreDM', 'DM', 'PreDM_Dapa', 'DM_Dapa'}
STEPS = ['01_qc_dea.R', '02_reversal.R', '03_enrichment.R',
         '05_axis_specific_limma.R', '06_axis_specific_enrichment.R',
         '07_reference_reversal_figure.R']
REQUIRED_R = ['limma','readr','dplyr','tidyr','purrr','ggplot2','ggrepel',
              'pheatmap','impute','matrixStats','fgsea','msigdbr','patchwork','scales','statmod']

def read_table(path, delimiter):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        names = reader.fieldnames
        if not names or len(names) != len(set(names)):
            raise ValueError(f'Missing or duplicate column names: {path.name}')
        rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError(f'Empty or nonrectangular table: {path.name}')
    return names, rows

def validate_inputs(directory):
    _, meta = read_table(directory/'sample_meta.csv', ',')
    fields, abundance = read_table(directory/'abundance_input_log2.tsv', '\t')
    afields, annotation = read_table(directory/'protein_annot.tsv', '\t')
    required_meta = {'sample_id','group','rat_id','biological_unit','dapagliflozin'}
    if not required_meta <= set(meta[0]): raise ValueError('Missing required sample metadata columns')
    ids = [r['sample_id'] for r in meta]
    if len(set(ids)) != len(ids) or any(not x for x in ids): raise ValueError('Duplicate or empty sample IDs')
    rats = [r['rat_id'] for r in meta]
    if len(set(rats)) != len(rats) or any(not x for x in rats): raise ValueError('Duplicate or empty rat IDs')
    counts = Counter(r['group'] for r in meta)
    if set(counts) != GROUPS or any(n != 4 for n in counts.values()):
        raise ValueError('Historical pipeline requires the documented five groups with four samples per group')
    if set(fields) != {'Accession', *ids}: raise ValueError('Abundance sample columns must exactly match metadata')
    for r in meta:
        if r['biological_unit'] != 'rat': raise ValueError('Historical study requires rat biological units')
        expected = 'yes' if r['group'].endswith('_Dapa') else 'no'
        if r['dapagliflozin'] != expected: raise ValueError('Treatment metadata disagrees with group')
    accessions = [r['Accession'] for r in abundance]
    if any(not x for x in accessions) or len(set(accessions)) != len(accessions):
        raise ValueError('Duplicate or empty abundance accessions')
    if not {'Accession','gene_symbol','description'} <= set(afields): raise ValueError('Missing annotation columns')
    aids=[r['Accession'] for r in annotation]
    if len(set(aids))!=len(aids) or set(aids)!=set(accessions): raise ValueError('Annotation accessions must match abundance exactly')
    missing=0
    for row in abundance:
        observed=0
        for sample in ids:
            value=row[sample]
            if value in ('','NA','NaN'): missing+=1; continue
            if not math.isfinite(float(value)): raise ValueError('Nonfinite abundance value; use NA for missing')
            observed+=1
        if not observed: raise ValueError('Entirely missing protein row')
    return {'proteins':len(abundance),'samples':len(meta),'group_counts':dict(counts),'missing_cells':missing,
            'scale':'processed log2, asserted by historical source; not inferred from value range'}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare(export, run_dir):
    export=export.resolve();run_dir=run_dir.resolve()
    if run_dir.exists(): raise ValueError(f'Refusing to overwrite existing run directory: {run_dir}')
    if run_dir.is_relative_to(export) or export.is_relative_to(run_dir): raise ValueError('Run directory and historical export must be separate')
    src=export/'01_INPUT_DATA'/'02_processed_rat_level'
    inputs={name:src/name for name in ['sample_meta.csv','abundance_input_log2.csv','protein_annot.csv']}
    for p in inputs.values():
        if not p.is_file(): raise ValueError(f'Missing source input: {p}')
    # Read every source first; quoted CSV is parsed, not split on commas.
    loaded={name:read_table(p,',') for name,p in inputs.items()}
    proc=run_dir/'data'/'processed';proc.mkdir(parents=True)
    for name,(fields,rows) in loaded.items():
        target=proc/(name if name=='sample_meta.csv' else name.replace('.csv','.tsv'))
        with target.open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter=',' if target.suffix=='.csv' else '\t')
            w.writeheader();w.writerows(rows)
    validation=validate_inputs(proc)
    provenance={'source':'archived Method B export; metadata claims inherited, not independently confirmed',
                'created_utc':datetime.now(timezone.utc).isoformat(),
                'inputs':[{'source_name':n,'sha256':digest(p)} for n,p in inputs.items()],
                'prepared':[{'name':p.name,'sha256':digest(p)} for p in sorted(proc.iterdir())],
                'validation':validation}
    (run_dir/'input_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    return provenance

def preflight(run_dir, rscript):
    validation=validate_inputs(run_dir/'data'/'processed')
    executable=shutil.which(rscript)
    result={'input_validation':validation,'rscript':executable,'r_environment':'unavailable'}
    if not executable: return result,False
    expression='p <- c('+','.join(json.dumps(p) for p in REQUIRED_R)+'); missing <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop(paste("Missing R packages:",paste(missing,collapse=", "))); cat(R.version.string,"\\n"); for(x in p) cat(x,as.character(packageVersion(x)),"\\n")'
    completed=subprocess.run([executable,'--vanilla','-e',expression],capture_output=True,text=True)
    result['r_environment']=completed.stdout+completed.stderr
    return result,completed.returncode==0

def run(run_dir, rscript):
    status,ok=preflight(run_dir,rscript)
    if not ok: raise ValueError(json.dumps(status,indent=2))
    # A run is immutable after execution starts; use a new directory for retries.
    if (run_dir/'execution.json').exists() or (run_dir/'results').exists():
        raise ValueError('Existing execution/results found; prepare a new run directory')
    for folder in ['logs','report','reproducibility']: (run_dir/folder).mkdir(exist_ok=True)
    env=os.environ.copy();env['DAPA_ANALYSIS_ROOT']=str(run_dir.resolve())
    record={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'RUNNING','preflight':status,'steps':[],
            'limitation':'Historical statistical logic retained; execution does not resolve scientific audit findings.'}
    def save(): (run_dir/'execution.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    save()
    for name in STEPS:
        script=REPO/'pipeline'/name
        cmd=[status['rscript'],'--vanilla',str(REPO/'scripts'/'run_stage.R'),str(script),str(run_dir/'reproducibility'/(name+'.session.txt'))]
        print(f'Running {name}',flush=True)
        with (run_dir/'logs'/(name+'.log')).open('w',encoding='utf-8') as log:
            completed=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,text=True)
        record['steps'].append({'script':name,'sha256':digest(script),'exit_code':completed.returncode})
        save()
        if completed.returncode:
            record['status']='FAILED';save();raise ValueError(f'{name} failed; inspect its log')
    record['status']='COMPLETED_WITH_KNOWN_SCIENTIFIC_LIMITATIONS'
    record['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    (run_dir/'report'/'execution_summary.md').write_text('# Historical pipeline execution\n\nSix recovered R stages completed. This verifies execution, not scientific validity.\n\nSee execution.json for commands, versions and status, and repository/docs/audit/SCIENTIFIC_AUDIT.md for unresolved interpretation issues. The historical finalizer with hardcoded PASS claims was not executed.\n',encoding='utf-8')
    manifest=[{'path':p.relative_to(run_dir).as_posix(),'sha256':digest(p)} for p in sorted(run_dir.rglob('*')) if p.is_file() and p.name!='manifest.json']
    (run_dir/'reproducibility'/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return record

def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare');p.add_argument('--export',type=Path,required=True);p.add_argument('--run-dir',type=Path,required=True)
    for name in ['preflight','run']:
        p=sub.add_parser(name);p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--rscript',default='Rscript')
    args=ap.parse_args()
    try:
        if args.command=='prepare': result=prepare(args.export,args.run_dir)
        elif args.command=='preflight':
            result,ok=preflight(args.run_dir,args.rscript);print(json.dumps(result,indent=2));return 0 if ok else 2
        else: result=run(args.run_dir,args.rscript)
        print(json.dumps(result,indent=2));return 0
    except (ValueError,OSError,KeyError) as error:
        print(f'ERROR: {error}',file=sys.stderr);return 2
if __name__=='__main__': raise SystemExit(main())
