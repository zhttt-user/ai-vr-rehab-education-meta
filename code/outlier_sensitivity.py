# -*- coding: utf-8 -*-
"""离群点敏感性分析：剔除95%CI与合并CI不重叠的研究，看重合并是否逆转/大幅波动"""
import sys, io, math
from pathlib import Path
sys.stdout = io.StringIO()
import importlib.util
spec = importlib.util.spec_from_file_location(
    'ma', str(Path(__file__).resolve().parent / 'meta_analysis.py'))
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
sys.stdout = sys.__stdout__

OUT = Path(__file__).resolve().parent.parent / 'results'
report = ['\n## 补充验证2：离群点敏感性分析（漏斗图CI不重叠法）\n']
report.append('离群点定义：该研究95%CI与本亚组合并效应95%CI完全不重叠（Harrer/dmetar惯例）。\n')
report.append('| 池 | 亚组 | 原合并g [95%CI] | 离群研究 | 剔除后k | 剔除后g [95%CI] | 变化 | 逆转？ |')
report.append('|---|---|---|---|---|---|---|---|')

for pname, data in ma.pools.items():
    eff = ma.build_effects(data)
    for grp in ('AI', 'VR'):
        sub = {k: v for k, v in eff.items() if v['grp'] == grp}
        if len(sub) < 3:
            continue
        es = [(v['g'], v['var']) for v in sub.values()]
        pool = ma.dl_pool(es)
        plo, phi = pool['lo'], pool['hi']
        outliers = []
        for k, v in sub.items():
            se = math.sqrt(v['var'])
            lo, hi = v['g']-1.96*se, v['g']+1.96*se
            if hi < plo or lo > phi:
                outliers.append(k)
        if not outliers:
            report.append(f'| {pname} | {grp} | {pool["g"]:.2f} [{plo:.2f},{phi:.2f}] | 无 | — | — | — | — |')
            continue
        es2 = [(v['g'], v['var']) for k, v in sub.items() if k not in outliers]
        pool2 = ma.dl_pool(es2)
        chg = pool2['g'] - pool['g']
        flip = '方向逆转!' if (pool['g'] > 0) != (pool2['g'] > 0) else ('显著性消失' if pool['p'] < 0.05 <= pool2['p'] else '否')
        report.append(f'| {pname} | {grp} | {pool["g"]:.2f} [{plo:.2f},{phi:.2f}] | {"、".join(o.split(" ",1)[-1] if " " in o else o for o in outliers)} | '
                      f'{pool2["k"]} | {pool2["g"]:.2f} [{pool2["lo"]:.2f},{pool2["hi"]:.2f}] | {chg:+.2f} | {flip} |')

report.append('\n注：高异质性背景下，CI不重叠法会同时标记正向大效应研究（如朱丹、王树本、张通），')
report.append('剔除目的是诊断合并结果的依赖度，而非主张这些研究应被排除。')
txt = '\n'.join(report)
with open(OUT/'稳健性交叉验证报告.md', 'a', encoding='utf-8') as f:
    f.write('\n' + txt)
print(txt)
