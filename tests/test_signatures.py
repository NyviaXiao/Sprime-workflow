import ast
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def copied_project():
    """Copy only signature inputs to an isolated temporary project root."""
    tmp = Path(tempfile.mkdtemp())
    scripts = tmp / 'workflow/scripts'
    scripts.mkdir(parents=True)
    for path in (ROOT / 'workflow/scripts').iterdir():
        if path.is_file():
            shutil.copy2(path, scripts / path.name)
    return tmp


def signatures(project_root):
    """Evaluate the Snakefile preamble against a temporary script directory."""
    from run import resolve_config

    config = resolve_config(ROOT / 'config/config.yaml', {})
    config['project_root'] = str(project_root)
    namespace = {'config': config}
    source = (ROOT / 'workflow/Snakefile').read_text(encoding='utf-8')
    exec(source.split('wildcard_constraints:')[0], namespace)
    return namespace['SIG']


def mutate_task_branch(path, stage):
    """Insert a harmless assignment into exactly one main() elif body."""
    tree = ast.parse(path.read_text(encoding='utf-8'))
    main = next(node for node in tree.body
                if isinstance(node, ast.FunctionDef) and node.name == 'main')
    for node in main.body:
        if not isinstance(node, ast.If):
            continue
        branch = node
        while isinstance(branch, ast.If):
            if any(isinstance(x, ast.Constant) and x.value == stage
                   for x in ast.walk(branch.test)):
                branch.body.insert(0, ast.Assign(
                    targets=[ast.Name(id='_signature_test_marker', ctx=ast.Store())],
                    value=ast.Constant(value=1)))
                ast.fix_missing_locations(tree)
                path.write_text(ast.unparse(tree) + '\n', encoding='utf-8')
                return
            branch = branch.orelse[0] if len(branch.orelse) == 1 else None
    raise AssertionError(f'No task branch for {stage}')


def mutate_core_function(path, name):
    """Insert a harmless assignment into one selected core helper."""
    tree = ast.parse(path.read_text(encoding='utf-8'))
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == name)
    function.body.insert(0, ast.Assign(
        targets=[ast.Name(id='_signature_test_marker', ctx=ast.Store())],
        value=ast.Constant(value=1)))
    ast.fix_missing_locations(tree)
    path.write_text(ast.unparse(tree) + '\n', encoding='utf-8')


class SignatureTests(unittest.TestCase):
    def test_branch_mutations_are_independent(self):
        for stage, changed in [('individual', 'individual'),
                               ('gmm', 'gmm'), ('classify', 'classify')]:
            with self.subTest(stage=stage):
                project = copied_project()
                try:
                    before = signatures(project)
                    mutate_task_branch(project / 'workflow/scripts/task.py', stage)
                    after = signatures(project)
                    self.assertNotEqual(after[changed], before[changed])
                    for other in ('gmm', 'classify', 'individual'):
                        if other != changed:
                            self.assertEqual(after[other], before[other])
                finally:
                    shutil.rmtree(project)

    def test_unrelated_task_import_does_not_invalidate_stages(self):
        project = copied_project()
        try:
            before = signatures(project)
            path = project / 'workflow/scripts/task.py'
            tree = ast.parse(path.read_text(encoding='utf-8'))
            tree.body.insert(0, ast.Import(names=[ast.alias(name='math', asname='signature_test')]))
            ast.fix_missing_locations(tree)
            path.write_text(ast.unparse(tree) + '\n', encoding='utf-8')
            after = signatures(project)
            self.assertEqual(after, before)
        finally:
            shutil.rmtree(project)

    def test_shared_helpers_change_only_their_consumers(self):
        project = copied_project()
        try:
            before = signatures(project)
            mutate_core_function(project / 'workflow/scripts/core.py', 'classify')
            after = signatures(project)
            self.assertNotEqual(after['classify'], before['classify'])
            self.assertEqual(after['gmm'], before['gmm'])
            self.assertEqual(after['individual'], before['individual'])
        finally:
            shutil.rmtree(project)

        project = copied_project()
        try:
            before = signatures(project)
            mutate_core_function(project / 'workflow/scripts/core.py', 'gmm_pass')
            after = signatures(project)
            self.assertNotEqual(after['gmm_input'], before['gmm_input'])
            self.assertEqual(after['gmm'], before['gmm'])
            self.assertEqual(after['classify'], before['classify'])
            self.assertEqual(after['individual'], before['individual'])
        finally:
            shutil.rmtree(project)

    def test_stage_signature_inventory(self):
        project = copied_project()
        try:
            sig = signatures(project)
            self.assertEqual(len(sig), 20)
        finally:
            shutil.rmtree(project)


if __name__ == '__main__':
    unittest.main()
