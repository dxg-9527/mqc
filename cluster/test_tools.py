"""Runner mechanics only; synthetic fixtures are not biological validation."""
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock
import tarfile

from common import load_inputs, sha, under, valid_completed, write_csv, write_json
from submit import assignments
from worker import compress_evidence, run_one
from check_environment import inspect


class ClusterTests(unittest.TestCase):
    def test_environment_rejects_unvalidated_default_solver(self):
        cobra = SimpleNamespace(Configuration=lambda: SimpleNamespace(
            solver=SimpleNamespace(__name__='optlang.gurobi_interface')))
        with mock.patch.dict('sys.modules',{'cobra':cobra,'cplex':SimpleNamespace()}):
            result=inspect(probe=True)
        self.assertEqual(result['cplex_probe'],'failed')
        self.assertIn('validated environment uses CPLEX',result['probe_error'])

    def test_assignment_is_complete_unique_balanced(self):
        rows=[dict(job_id=str(i),reactions=str(i+10)) for i in range(6543)]
        groups=assignments(rows,1000)
        flattened=[x for g in groups for x in g]
        self.assertEqual(len(groups),1000)
        self.assertEqual(len(flattened),6543)
        self.assertEqual(set(flattened),{r['job_id'] for r in rows})
        self.assertEqual(assignments(rows,1000),groups)

    def test_escape_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(ValueError): under(folder,'../outside')

    def test_atomic_json_normalizes_nonfinite(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'status.json'
            write_json(path,dict(value=float('nan')))
            self.assertIsNone(json.loads(path.read_text())['value'])

    def test_modified_input_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            rows=[dict(job_id='test_1',input_relative_path='models/test.xml.gz')]
            write_csv(root/'manifest.csv',rows,list(rows[0]))
            write_json(root/'bundle.json',dict(model_count=1,manifest_sha256=sha(root/'manifest.csv')))
            self.assertEqual(len(load_inputs(root)[0]),1)
            (root/'manifest.csv').write_text('bad')
            with self.assertRaises(ValueError): load_inputs(root)

    def test_archive_excludes_regenerable_input_and_keeps_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); work=root/'work'; attempt=root/'attempt'
            work.mkdir(); attempt.mkdir(); (work/'input').mkdir()
            (work/'input/model.xml').write_text('original')
            (work/'pipeline.log').write_text('test failure evidence')
            digest=compress_evidence(work,attempt)
            self.assertEqual(digest,sha(attempt/'evidence.tar.gz'))
            with tarfile.open(attempt/'evidence.tar.gz') as tar:
                self.assertEqual(tar.getnames(),['pipeline.log'])
            self.assertTrue((work/'input/model.xml').exists())

    def test_failed_attempt_is_archived_and_original_is_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); inputs=root/'inputs'; inputs.mkdir()
            model=inputs/'model.xml.gz'
            xml=b'<sbml/>'
            with gzip.open(model,'wb') as stream: stream.write(xml)
            row=dict(job_id='test_1',database='BiGG',stable_model_id='test',input_filename='test.xml',
                     input_relative_path='model.xml.gz',gzip_sha256=sha(model),
                     xml_sha256=hashlib.sha256(xml).hexdigest())
            plan=dict(results=str(root/'results'),inputs=str(inputs),release_lock_sha256='release',
                      input_manifest_sha256='manifest',model_timeout_seconds=10)
            class FakeChild:
                def __init__(self,command,**kw):
                    work=Path(command[command.index('--work')+1])
                    write_json(work/'outcome.json',dict(status='failed',error='fixture_failure'))
                def wait(self,timeout=None): return 1
            with mock.patch('worker.subprocess.Popen',FakeChild),mock.patch.dict('os.environ',{},clear=True):
                self.assertFalse(run_one(row,plan,{}))
            attempts=list((root/'results/tasks/test_1').glob('attempt_*'))
            self.assertEqual(len(attempts),1)
            record=json.loads((attempts[0]/'status.json').read_text())
            self.assertEqual(record['status'],'failed')
            self.assertTrue((attempts[0]/'evidence.tar.gz').is_file())
            self.assertEqual(sha(model),row['gzip_sha256'])
            self.assertFalse(Path(record['scratch']).exists())

    def test_completed_requires_identity_and_archive_hash(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder); (p/'evidence.tar.gz').write_bytes(b'test')
            record=dict(status='completed',release_lock_sha256='r',input_manifest_sha256='m',
                        input_xml_sha256='i',archive_sha256=sha(p/'evidence.tar.gz'))
            write_json(p/'status.json',record)
            self.assertTrue(valid_completed(p,'r','m',dict(xml_sha256='i')))
            self.assertFalse(valid_completed(p,'wrong','m',dict(xml_sha256='i')))
            (p/'evidence.tar.gz').write_bytes(b'changed')
            self.assertFalse(valid_completed(p,'r','m',dict(xml_sha256='i')))


if __name__=='__main__': unittest.main()
