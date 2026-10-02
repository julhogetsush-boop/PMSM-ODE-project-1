"""Tutorial 4: reproduce a reported number, review formulas, compare cost.

python code/week4_reproduction.py
python code/week4_reproduction.py --require-clean-commit --expected-commit SHA
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp

from methods import solve_fixed
from model import load_parameters, matrices, make_model, stability_limit
from week3_review import run_week3
from checks import run_checks

ROOT = Path(__file__).resolve().parents[1]


def git_evidence(require_clean=False, expected_commit=None):
    git = shutil.which('git')
    if git is None:
        if require_clean or expected_commit:
            raise RuntimeError('Git is required for exact-commit verification.')
        return {'status': 'PENDING_REAL_COMMIT', 'commit': None}
    def call(*args):
        return subprocess.run([git, '-C', str(ROOT), *args], text=True,
                              capture_output=True, check=False)
    head = call('rev-parse', '--verify', 'HEAD')
    if head.returncode:
        if require_clean or expected_commit:
            raise RuntimeError('No real Git commit found. Commit the reviewed files first.')
        return {'status': 'PENDING_REAL_COMMIT', 'commit': None}
    commit = head.stdout.strip()
    if expected_commit:
        wanted = call('rev-parse', '--verify', expected_commit + '^{commit}')
        if wanted.returncode or wanted.stdout.strip() != commit:
            raise RuntimeError('HEAD does not match --expected-commit.')
    source_files = list((ROOT / 'code').glob('*.py')) + list((ROOT / 'report').glob('*.tex')) + [ROOT / 'parameters.json',
                    ROOT / 'requirements.txt', ROOT / 'report/report_claims.json']
    missing = [p.relative_to(ROOT).as_posix() for p in source_files
               if call('ls-files', '--error-unmatch', p.relative_to(ROOT).as_posix()).returncode]
    dirty = call('status', '--porcelain', '--untracked-files=all', '--', '.').stdout.splitlines()
    clean = not missing and not dirty
    if require_clean and not clean:
        raise RuntimeError('Commit verification needs tracked source files and a clean folder. '
                           f'Untracked sources: {missing}; changes: {dirty}')
    return {'status': 'EXACT_CLEAN_COMMIT' if clean else 'COMMIT_WITH_LOCAL_CHANGES',
            'commit': commit, 'untracked_sources': missing, 'changes_before_run': dirty}


def source_fingerprint():
    paths = list((ROOT / 'code').glob('*.py')) + list((ROOT / 'report').glob('*.tex')) + [ROOT / 'parameters.json',
                    ROOT / 'requirements.txt', ROOT / 'report/report_claims.json']
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(paths)}
    combined = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    return {'sha256': combined, 'files': hashes}


def mathematical_review(p):
    A, b = matrices(p)
    D = np.diag([p['Ld'], p['Lq']])
    P = D @ D
    B = D @ A @ np.linalg.inv(D)
    expected_B = np.array([[-p['R']/p['Ld'], p['omega_e']],
                           [-p['omega_e'], -p['R']/p['Lq']]])
    lyapunov = A.T @ P + P @ A
    target = np.diag([-2*p['R']*p['Ld'], -2*p['R']*p['Lq']])
    assert np.allclose(B, expected_B, rtol=1e-13, atol=1e-13)
    assert np.allclose(lyapunov, target, rtol=1e-13, atol=1e-13)
    assert np.all(np.linalg.eigvalsh(lyapunov) < 0)
    det_expected = p['R']**2/(p['Ld']*p['Lq']) + p['omega_e']**2
    assert np.isclose(np.linalg.det(A), det_expected, rtol=1e-12)
    return {'status': 'PASS', 'jacobian_per_s': A.tolist(), 'forcing_A_per_s': b.tolist(),
            'transformed_jacobian_per_s': B.tolist(), 'lyapunov_matrix': lyapunov.tolist(),
            'lyapunov_identity_max_defect': float(np.max(abs(lyapunov-target))),
            'determinant_per_s2': float(np.linalg.det(A)),
            'scope': 'Global constant-coefficient perturbation result; excludes variable speed and saturation.'}


def cost_at_matched_accuracy(p):
    f, jac = make_model(p)
    ref = solve_ivp(f, (0, p['T']), p['y0'], method='Radau', jac=jac,
                    rtol=1e-12, atol=1e-14, dense_output=True)
    if not ref.success:
        raise RuntimeError(ref.message)
    targets = [1e-2, 5e-3]
    sweep, matched = [], []
    for method in ['euler', 'rk4', 'implicit']:
        remaining = set(targets)
        for N in [30*2**k for k in range(13)]:
            calls = {'rhs': 0, 'jacobian': 0}
            def rhs(t, y):
                calls['rhs'] += 1
                return f(t, y)
            def counted_jac(t, y):
                calls['jacobian'] += 1
                return jac(t, y)
            t, Y, counts = solve_fixed(rhs, p['y0'], 0, p['T'], p['T']/N,
                                       method=method, jac=counted_jac)
            error = float(np.max(np.abs(Y-ref.sol(t).T)))
            row = {'method': method, 'N': N, 'h_s': p['T']/N,
                   'max_grid_component_error_A': error,
                   'rhs_evaluations': calls['rhs'], 'jacobian_evaluations': calls['jacobian'],
                   'newton_updates': counts['newton_updates'], 'steps': counts['steps']}
            sweep.append(row)
            for target in sorted(remaining, reverse=True):
                if error <= target:
                    matched.append(dict(target_error_A=target, **row))
                    remaining.remove(target)
            if not remaining:
                break
        if remaining:
            raise RuntimeError(f'{method} did not attain targets {remaining}; extend refinement grid.')
    for name, rows in [('week4_cost_sweep.csv', sweep), ('week4_matched_accuracy.csv', matched)]:
        with (ROOT / 'results' / name).open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for method, colour, marker in [('euler', '#0072B2', 'o'), ('rk4', '#009E73', 's'),
                                    ('implicit', '#D55E00', '^')]:
        # Keep very large coarse-grid errors in CSV; the figure focuses on target crossing.
        rows = [r for r in sweep if r['method'] == method
                and r['max_grid_component_error_A'] <= 100]
        ax.loglog([r['rhs_evaluations'] for r in rows],
                  [r['max_grid_component_error_A'] for r in rows],
                  marker=marker, color=colour, label=method)
    for target in targets:
        ax.axhline(target, color='0.4', linestyle='--', linewidth=0.9,
                   label=f'shared error target {target:g} A')
    ax.set(xlabel='RHS evaluations [count]', ylabel='Max grid/component error [A]',
           title='RK4 needs fewer RHS evaluations at the same error target')
    ax.grid(True, which='both', linestyle=':', alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/week4_cost_accuracy.png', dpi=180)
    plt.close(fig)
    best = {str(target): min([r for r in matched if r['target_error_A']==target],
                             key=lambda r: r['rhs_evaluations'])['method'] for target in targets}
    # Guard the figure's stated conclusion if parameters are changed.
    if not all(method == 'rk4' for method in best.values()):
        raise RuntimeError('The cost figure conclusion needs revision for these parameters.')
    return {'metric': 'RHS evaluations including Newton residual evaluations',
            'error_metric': 'maximum component error on each solver grid versus tight Radau',
            'targets_A': targets, 'matched': matched, 'minimum_rhs_method': best,
            'selection': 'first qualifying grid in a powers-of-two sweep; shared accuracy ceilings, not equal errors',
            'excluded_cost': 'Jacobian assembly, linear solves, setup, reference solve, plotting; no CPU-time claim'}


def write_report(cost, math_review, reproduced):
    numbers = {row['id']: row for row in reproduced}
    rows = ['# Week 4 · 实际复现结果', '', '本页由 `python code/run_all.py` 或 '
            '`python code/week4_reproduction.py` 生成。参数为教学示例。', '',
            '## 与报告逐项比对', '', '| 数字 | 报告值（3位有效数字） | 本次值 | 单位 | 状态 |',
            '|---|---:|---:|---|---|']
    for r in reproduced:
        rows.append(f"| {r['label']} | {r['expected_3sf']} | {r['actual_3sf']} | {r['units']} | PASS |")
    rows += ['', '## 匹配精度的成本比较', '',
             '| 共同误差上限 (A) | 方法 | N | 实测网格误差 (A) | RHS | Jacobian | Newton 修正 |',
             '|---:|---|---:|---:|---:|---:|---:|']
    for r in cost['matched']:
        rows.append(f"| {r['target_error_A']:g} | {r['method']} | {r['N']} | "
                    f"{r['max_grid_component_error_A']:.6g} | {r['rhs_evaluations']} | "
                    f"{r['jacobian_evaluations']} | {r['newton_updates']} |")
    rows += ['', '图展示什么：三种方法从各自细化网格得到的误差与 RHS 次数。'
             '为什么重要：统一误差上限之后才能比较计算成本。这里 RK4 的 RHS 次数最少；'
             '统计不含 Jacobian/线性求解成本，不能改写为运行时间排名。', '',
             '![Cost at shared accuracy targets](../figures/week4_cost_accuracy.png)', '',
             f"Lyapunov 恒等式数值缺陷：{math_review['lyapunov_identity_max_defect']:.3e}。",
             '', 'Git commit、运行版本、源码哈希和命令见 `results/week4_reproduction.json`。'
             '未建立真实提交时，Git 状态会明确显示 PENDING_REAL_COMMIT。']
    (ROOT / 'notes/04_Week4_实际运行结果.md').write_text('\n'.join(rows)+'\n', encoding='utf-8')
    english = ['# Week 4 reproduction and cost evidence', '',
        f"For the documented illustrative PMSM case, the RK4 endpoint-order estimate was "
        f"recomputed as {numbers['rk4_endpoint_order']['actual']:.16g}, matching the report "
        f"to three significant digits ({numbers['rk4_endpoint_order']['actual_3sf']}). "
        f"The Euler stability boundary is {numbers['euler_stability_ms']['actual']:.16g} ms "
        f"({numbers['euler_stability_ms']['actual_3sf']} ms to three significant digits). "
        'Reproduction uses the exact parameter file and the scripts shipped with this folder. '
        'The source SHA256 and environment are recorded with the run; a genuine Git commit '
        'must additionally be recorded before claiming the exact-commit Tutorial 4 requirement.', '',
        '## Matched accuracy study', '',
        'Cost is counted in right-hand-side evaluations, including implicit Newton residual '
        'evaluations. The same error ceiling is imposed on all methods; the first qualifying '
        'grid in the recorded powers-of-two sweep is selected. Errors are maximum component '
        'errors on each method\'s own time grid versus a Radau reference '
        '(rtol=1e-12, atol=1e-14). This is sampled trajectory error, not a proven supremum '
        'between nodes. Jacobian counts and Newton updates are reported separately. '
        'Jacobian assembly, linear solves, reference construction, setup and plotting are '
        'excluded from the primary count, so this experiment does not rank wall-clock time.', '',
        '| Error ceiling [A] | Method | N | Measured error [A] | RHS calls |',
        '|---:|---|---:|---:|---:|']
    for r in cost['matched']:
        english.append(f"| {r['target_error_A']:g} | {r['method']} | {r['N']} | "
                       f"{r['max_grid_component_error_A']:.6g} | {r['rhs_evaluations']} |")
    english += ['', '**Figure caption.** Cost and accuracy for the same constant-speed PMSM '
        'IVP. Each curve shows the measured error against Radau versus counted RHS evaluations. '
        'Dashed lines indicate the shared accuracy ceilings. RK4 reaches both ceilings with '
        'fewer RHS evaluations than either first-order method in this experiment, which matters '
        'because an efficiency comparison must control the attained accuracy.', '',
        '![Cost and accuracy](../figures/week4_cost_accuracy.png)', '',
        'The model has a pair of equally damped complex eigenvalues; its explicit-step restriction '
        'is primarily oscillatory, rather than evidence of widely separated decay scales. '
        'The cost result is specific to this affine teaching case.']
    (ROOT / 'report/week4_verified_results.md').write_text('\n'.join(english)+'\n', encoding='utf-8')


def run_week4(week3=None, checks=None, require_clean=False, expected_commit=None,
              command='python code/week4_reproduction.py'):
    for part in ['results', 'figures', 'notes', 'report']:
        (ROOT / part).mkdir(exist_ok=True)
    # Inspect source provenance before this run rewrites generated outputs.
    provenance = git_evidence(require_clean, expected_commit)
    fingerprint = source_fingerprint()
    p = load_parameters()
    checks = run_checks() if checks is None else checks
    week3 = run_week3() if week3 is None else week3
    actual = {'rk4_endpoint_order': week3['slope'],
              'euler_stability_ms': 1000*stability_limit('euler', np.linalg.eigvals(matrices(p)[0]))}
    claims = json.loads((ROOT / 'report/report_claims.json').read_text(encoding='utf-8'))['claims']
    reproduced = []
    lines = []
    for claim in claims:
        expected_text = format(claim['expected'], '.3g')
        actual_text = format(actual[claim['id']], '.3g')
        if expected_text != actual_text:
            raise RuntimeError(f"Report mismatch for {claim['id']}: report={expected_text}, "
                               f"run={actual_text}. Review and fix the report before rerunning.")
        reproduced.append(dict(claim, actual=actual[claim['id']], expected_3sf=expected_text,
                               actual_3sf=actual_text, status='PASS'))
        line = f"REPRODUCED {claim['id']}: report={expected_text}; run={actual_text}; units={claim['units']}; PASS"
        lines.append(line)
        print(line, flush=True)
    math_review = mathematical_review(p)
    cost = cost_at_matched_accuracy(p)
    write_report(cost, math_review, reproduced)
    if require_clean:
        command += ' --require-clean-commit'
    if expected_commit:
        command += ' --expected-commit ' + expected_commit
    result = {'status': 'PASS', 'timestamp_utc': datetime.now(timezone.utc).isoformat(),
              'command': command, 'git': provenance, 'source_snapshot': fingerprint,
              'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                           'scipy': scipy.__version__, 'matplotlib': matplotlib.__version__},
              'reproduced_claims': reproduced, 'independent_checks': checks,
              'mathematical_review': math_review, 'cost_accuracy': cost,
              'course_submission_status': 'Pending real team report/slides, personal ICS, peer review and PM posts.'}
    (ROOT / 'results/week4_reproduction.json').write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    (ROOT / 'results/week4_reproduction_line.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print('Week 4 numerical reproduction and mathematical checks PASS.', flush=True)
    print('Exact-commit status:', provenance['status'], flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-clean-commit', action='store_true')
    parser.add_argument('--expected-commit')
    args = parser.parse_args()
    run_week4(require_clean=args.require_clean_commit, expected_commit=args.expected_commit)
