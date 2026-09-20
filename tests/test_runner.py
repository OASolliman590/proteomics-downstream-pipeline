import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('runner',Path(__file__).resolve().parents[1]/'scripts'/'run_pipeline.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)

class InputValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.meta=[{'sample_id':f'{g}_{n}','rat_id':f'{g}_{n}','group':g,'biological_unit':'rat','dapagliflozin':'yes' if g.endswith('_Dapa') else 'no'} for g in sorted(runner.GROUPS) for n in range(1,5)]
        self.ab=[{'Accession':'P1',**{m['sample_id']:'24.0' for m in self.meta}}]
        self.ann=[{'Accession':'P1','gene_symbol':'Gene1','description':'Synthetic protein, for tests'}]
    def write(self,name,rows,delimiter):
        with (self.root/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=rows[0],delimiter=delimiter);w.writeheader();w.writerows(rows)
    def validate(self):
        self.write('sample_meta.csv',self.meta,',');self.write('abundance_input_log2.tsv',self.ab,'\t');self.write('protein_annot.tsv',self.ann,'\t')
        return runner.validate_inputs(self.root)
    def test_valid_missing_preserved(self):
        self.ab[0][self.meta[0]['sample_id']]='NA'
        self.assertEqual(self.validate()['missing_cells'],1)
    def test_duplicate_rat_rejected(self):
        self.meta[1]['rat_id']=self.meta[0]['rat_id']
        with self.assertRaisesRegex(ValueError,'rat IDs'):self.validate()
    def test_unknown_group_rejected(self):
        self.meta[0]['group']='PDM'
        with self.assertRaisesRegex(ValueError,'five groups'):self.validate()
    def test_treatment_disagreement_rejected(self):
        self.meta[0]['dapagliflozin']='yes'
        with self.assertRaisesRegex(ValueError,'Treatment'):self.validate()
    def test_sample_mismatch_rejected(self):
        self.ab[0].pop(self.meta[0]['sample_id'])
        with self.assertRaisesRegex(ValueError,'sample columns'):self.validate()
    def test_infinity_rejected(self):
        self.ab[0][self.meta[0]['sample_id']]='inf'
        with self.assertRaisesRegex(ValueError,'Nonfinite'):self.validate()
    def test_duplicate_protein_rejected(self):
        self.ab.append(self.ab[0].copy())
        with self.assertRaisesRegex(ValueError,'accessions'):self.validate()
    def test_annotation_mismatch_rejected(self):
        self.ann[0]['Accession']='P2'
        with self.assertRaisesRegex(ValueError,'Annotation'):self.validate()
    def test_existing_run_refused(self):
        with self.assertRaisesRegex(ValueError,'overwrite'):runner.prepare(self.root/'export',self.root)
    def test_missing_r_does_not_create_results(self):
        self.validate();proc=self.root/'data'/'processed';proc.mkdir(parents=True)
        for name in ['sample_meta.csv','abundance_input_log2.tsv','protein_annot.tsv']:(self.root/name).rename(proc/name)
        status,ok=runner.preflight(self.root,'nonexistent-rscript-for-unit-test')
        self.assertFalse(ok);self.assertIsNone(status['rscript']);self.assertFalse((self.root/'results').exists())

if __name__=='__main__':unittest.main()
