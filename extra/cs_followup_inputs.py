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


def verify(source_dir, base_dir, stage='s207'):
    source_dir, base_dir = Path(source_dir), Path(base_dir)
    result = json.loads((source_dir / 'result.json').read_text())
    diagnostic = json.loads((source_dir / 'diagnostic/report.json').read_text())
    audit_path = base_dir / 'final.audit.full.seed900001.json'
    audit = json.loads(audit_path.read_text())
    source_hash = digest(source_dir / 'final.pt')
    points_hash = digest(source_dir / 'diagnostic/diagnostic_points.pt')
    base_hash, audit_hash = digest(base_dir / 'final.pt'), digest(audit_path)
    if stage not in ('s207', 's307'):
        raise ValueError('Unknown reviewed source stage')
    seed, steps, selected, parent = ((207, 2000, 2000, PARENT_SHA256) if stage == 's207' else
                                    (307, 5000, 4750, '18abf4ce5e70a3bffe3d9264016cda7f4ff3141de84f38124f2e4e2c01bdb4d4'))
    expected = {'case': 's4xs1xs1', 'task': 'cs', 'mode': 'replay', 'seed': seed,
                'steps_done': steps, 'steps_planned': steps, 'selected_step': selected,
                'status': 'completed_pending_independent_audit'}
    if any(result.get(k) != v for k, v in expected.items()):
        raise ValueError(f'Expected the completed, selected step-{selected} {stage} replay run')
    if result.get('source', {}).get('sha256') != parent:
        raise ValueError('Replay run does not descend from the reviewed parent checkpoint')
    if stage == 's307' and result.get('config', {}).get('matched', {}).get('sampling', {}).get('target_count') != 4:
        raise ValueError('Expected the selected 12.5-percent replay arm')
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
    parser.add_argument('--stage', choices=['s207', 's307'], default='s207')
    args = parser.parse_args()
    try:
        print(verify(args.source_dir, args.base_dir, args.stage))
    except (OSError, ValueError) as exc:
        parser.exit(2, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
