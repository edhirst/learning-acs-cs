"""Standard-library-only preflight for the completed s407 long campaign."""
import argparse
import json
from pathlib import Path

try:
    from .cs_followup_inputs import ACS_SHA256, digest
except ImportError:  # direct login-node script invocation
    from cs_followup_inputs import ACS_SHA256, digest


def verify(source_run, base_dir):
    root, base = Path(source_run).resolve(), Path(base_dir).resolve()
    result_path = root / 'result.json'
    result = json.loads(result_path.read_text())
    audit = base / 'final.audit.full.seed900001.json'
    a = json.loads(audit.read_text())
    if (digest(base / 'final.pt') != ACS_SHA256 or a.get('sha256') != ACS_SHA256
            or a.get('audit_tier') != 'full' or a.get('case') != 's4xs1xs1'
            or a.get('numerically_validated') is not True or a.get('screening_targets_passed') is not True):
        raise ValueError('Expected the reviewed full-audited ACS base')
    if (result.get('status') != 'completed_pending_independent_audit'
            or result.get('steps_done') != 50000 or result.get('ready_for_diagnostics') is not True):
        raise ValueError('Expected the completed 50000-update long campaign')
    paths = [base / 'final.pt', audit, result_path]
    for name in ('average', 'tail'):
        ck = root / f'final_{name}.pt'
        report_path = root / f'diagnostic_{name}/report.json'
        points = root / f'diagnostic_{name}/diagnostic_points.pt'
        report = json.loads(report_path.read_text())
        if (report.get('status') != 'completed_diagnostic' or report.get('sha256') != digest(ck)
                or report.get('points_sha256') != digest(points)
                or report.get('audited_base', {}).get('sha256') != ACS_SHA256
                or report.get('audited_base', {}).get('audit_sha256') != digest(audit)):
            raise ValueError(f'{name}: diagnostic/checkpoint/base provenance mismatch')
        paths += [ck, report_path, points]
    parent = Path(result['source']['path'])
    if digest(parent) != result['source']['sha256']:
        raise ValueError('Previous source checkpoint mismatch')
    paths.append(parent)
    return {'source_run': str(root), 'base_dir': str(base), 'base_audit': str(audit),
            'previous': str(parent), 'files': {str(p): digest(p) for p in paths}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', required=True)
    parser.add_argument('--base-dir', required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.source_run, args.base_dir), sort_keys=True))
    except (OSError, KeyError, ValueError) as exc:
        parser.exit(2, f'Preflight failed: {exc}\n')
