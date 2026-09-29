"""Standard-library preflight for an audited precision campaign (no Torch)."""
import argparse
import json
import math
from pathlib import Path

try:
    from .cs_followup_inputs import digest
except ImportError:
    from cs_followup_inputs import digest


def verify(source_root):
    root = Path(source_root).resolve()
    files = {}

    def record(path, expected=None):
        path = Path(path).resolve()
        fingerprint = digest(path)
        if expected is not None and fingerprint != expected:
            raise ValueError(f'Input hash mismatch: {path}')
        files[str(path)] = fingerprint
        return fingerprint

    def read(name):
        path = root/name
        record(path)
        return json.loads(path.read_text())

    comparison = read('comparison.json')
    if (comparison.get('all_jobs_valid') is not True
            or comparison.get('full_audit_protocol') is not True
            or comparison.get('errors')):
        raise ValueError('Source must have a complete, valid full precision audit')
    selected = comparison['selected']
    checkpoint = Path(selected['checkpoint']).resolve()
    fingerprint = record(checkpoint)
    prepared = read('prepare.json')
    if prepared.get('status') != 'complete':
        raise ValueError('Source preparation is incomplete')
    for path, expected in prepared['inputs']['files'].items():
        record(path, expected)
    for path, expected in selected['artifacts'].items():
        record(path, expected)
    clouds = []
    for part in range(4):
        audit = read(f'audit_{part}.json')
        if (audit != comparison['audits'][part] or audit.get('status') != 'complete'
                or audit.get('sha256') != fingerprint
                or Path(audit['checkpoint']).resolve() != checkpoint
                or not audit['crosscheck']['passed'] or not audit['growth_passed']
                or audit['combined']['acs_max'] > 1e-10
                or audit['combined']['tangency_max'] > 1e-10):
            raise ValueError(f'Source audit {part} is inconsistent or numerically invalid')
        if not all(math.isfinite(audit['combined'][key]) and audit['combined'][key] >= 0
                   for key in ('acs_max', 'tangency_max')):
            raise ValueError('Nonfinite/negative source algebraic residual')
        if part == 0:
            if audit['groups']['uniform']['count'] < 100000 or audit['groups']['structured']['count'] < 20480:
                raise ValueError('Source uniform/structured audit is incomplete')
        elif audit['combined']['count'] < 16896:
            raise ValueError('Source adversarial audit is incomplete')
        for path, expected in audit['artifacts'].items():
            record(path, expected)
        cloud = root/f'audit_{part}.points.pt'
        if str(cloud) not in audit['artifacts']:
            raise ValueError('Source audit must bind its point cloud')
        clouds.append(str(cloud))
    result = dict(prepared['inputs'])
    result.update(source_root=str(root), checkpoint=str(checkpoint), checkpoint_sha256=fingerprint,
                  clouds=clouds, files=files)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.source_root), sort_keys=True))
    except (OSError, KeyError, ValueError, TypeError) as exc:
        parser.exit(2, f'Preflight failed: {exc}\n')
