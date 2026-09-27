"""Standard-library preflight for the reviewed s207 replay continuation.

Runs on the login node without importing Torch or starting any training.
"""
import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA256 = '886cba9805a2c83e1a7acc35cb14b5e8a0f7d3daa04206706bd862e7bc0ffd11'
ACS_SHA256 = '687b3164cf8aa91a202333ba8b5c4aeb0f5aa936f19d629c846990df6e2a0e19'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify(source_dir, base_dir):
    source_dir, base_dir = Path(source_dir), Path(base_dir)
    result = json.loads((source_dir / 'result.json').read_text())
    diagnostic = json.loads((source_dir / 'diagnostic/report.json').read_text())
    audit_path = base_dir / 'final.audit.full.seed900001.json'
    audit = json.loads(audit_path.read_text())
    source_hash = digest(source_dir / 'final.pt')
    points_hash = digest(source_dir / 'diagnostic/diagnostic_points.pt')
    base_hash, audit_hash = digest(base_dir / 'final.pt'), digest(audit_path)
    expected = {'case': 's4xs1xs1', 'task': 'cs', 'mode': 'replay', 'seed': 207,
                'steps_done': 2000, 'steps_planned': 2000, 'selected_step': 2000,
                'status': 'completed_pending_independent_audit'}
    if any(result.get(k) != v for k, v in expected.items()):
        raise ValueError('Expected the completed, selected step-2000 s207 replay run')
    if result.get('source', {}).get('sha256') != PARENT_SHA256:
        raise ValueError('Replay run does not descend from the reviewed s107 checkpoint')
    if (diagnostic.get('status') != 'completed_diagnostic'
            or diagnostic.get('sha256') != source_hash
            or diagnostic.get('points_sha256') != points_hash):
        raise ValueError('Replay diagnostic/checkpoint/witness checksum mismatch')
    if (base_hash != ACS_SHA256 or audit.get('sha256') != base_hash
            or audit.get('case') != 's4xs1xs1' or audit.get('audit_tier') != 'full'
            or audit.get('numerically_validated') is not True
            or audit.get('screening_targets_passed') is not True):
        raise ValueError('Expected the previously audited ACS base')
    for report in (result, diagnostic):
        evidence = report.get('audited_base', {})
        if evidence.get('sha256') != base_hash or evidence.get('audit_sha256') != audit_hash:
            raise ValueError('Source/diagnostic ACS provenance mismatch')
    return source_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', required=True)
    parser.add_argument('--base-dir', required=True)
    args = parser.parse_args()
    try:
        print(verify(args.source_dir, args.base_dir))
    except (OSError, ValueError) as exc:
        parser.exit(2, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
