"""Check Spec Kit traceability and local documentation links without running analysis."""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
EPIC=ROOT/'specs/001-downstream-proteomics'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args()
    trace=json.loads((EPIC/'traceability.json').read_text(encoding='utf-8'))['requirements']
    roadmap=json.loads((EPIC/'roadmap.json').read_text(encoding='utf-8'))['slices']
    assert len(trace)==120 and len(roadmap)==12
    for key in ['requirement','task','acceptance']:
        assert len({x[key] for x in trace})==120,f'Duplicate {key}'
    packet_ids={r['id'] for r in roadmap};seen=set()
    for row in roadmap:
        assert set(row['depends_on'])<=packet_ids,f'Unknown dependency: {row}'
        assert set(row['depends_on'])<=seen,f'Dependency order invalid: {row}'
        folder=ROOT/'specs'/row['slice']
        for filename in ['spec.md','plan.md','tasks.md']:assert (folder/filename).is_file()
        assigned=[t for t in trace if t['slice']==row['slice']]
        assert len(assigned)==10
        for t in assigned:
            assert t['requirement'] in (folder/'spec.md').read_text(encoding='utf-8')
            assert re.search(r'^- \[[ x]\] '+t['task']+r'\b',(folder/'tasks.md').read_text(encoding='utf-8'),re.M)
            assert t['acceptance'] in (folder/'spec.md').read_text(encoding='utf-8')
        seen.add(row['id'])
    broken=[]
    docs=list((ROOT/'specs').rglob('*.md'))+list((ROOT/'.specify').rglob('*.md'))
    for p in docs:
        txt=p.read_text(encoding='utf-8')
        assert '[NEEDS CLARIFICATION' not in txt and '[FEATURE NAME]' not in txt,p
        for link in re.findall(r'\]\(([^)]+)\)',txt):
            if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:',link) or link.startswith('#'):continue
            path=link.split('#',1)[0]
            if not (p.parent/path).exists():broken.append({'file':p.relative_to(ROOT).as_posix(),'link':link})
    assert not broken,broken
    for p in (ROOT/'specs').rglob('*.json'):json.loads(p.read_text(encoding='utf-8'))
    source_registry=json.loads((ROOT/'docs/audit/source_provenance.json').read_text(encoding='utf-8'))
    for row in source_registry:
        assert hashlib.sha256((ROOT/row['destination']).read_bytes()).hexdigest()==row['sha256_packaged'],row['destination']
    record={'status':'PASS','scope':'Specification structure/traceability/links and preserved-source hashes; not implementation validation',
            'slices':len(roadmap),'requirements':len(trace),'tasks':len(trace),'acceptance_cases':len(trace),'markdown_files_checked':len(docs),
            'implementation_status':'NOT_STARTED','schema_note':'Draft 2020-12 schema and three examples separately validated with jsonschema 4.25.1; semantic/runtime validation remains an implementation task.'}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))
if __name__=='__main__':main()
