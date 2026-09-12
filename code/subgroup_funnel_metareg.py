# -*- coding: utf-8 -*-
"""亚组漏斗图(AI/VR分开) + Meta回归(干预类型作为调节变量)
回应审稿建议：混合漏斗图的不对称可能源于亚组效应差异(混淆)而非发表偏倚
"""
import sys, io, math
from pathlib import Path
sys.stdout = io.StringIO()
import importlib.util
spec = importlib.util.spec_from_file_location(
    'ma', str(Path(__file__).resolve().parent / 'meta_analysis.py'))
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
sys.stdout = sys.__stdout__

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

OUT = Path(__file__).resolve().parent.parent / 'results'

def betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
    qab, qap, qam = a+b, a+1.0, a-1.0
    c, d = 1.0, 1.0-qab*x/qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0/d; h = d
    for m in range(1, MAXIT+1):
        m2 = 2*m
        aa = m*(b-m)*x/((qam+m2)*(a+m2))
        d = 1.0+aa*d;  d = FPMIN if abs(d)<FPMIN else d
        c = 1.0+aa/c;  c = FPMIN if abs(c)<FPMIN else c
        d = 1.0/d; h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d = 1.0+aa*d;  d = FPMIN if abs(d)<FPMIN else d
        c = 1.0+aa/c;  c = FPMIN if abs(c)<FPMIN else c
        d = 1.0/d; de = d*c; h *= de
        if abs(de-1.0) < EPS: break
    return h
def betai(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log(1-x))
    if x < (a+1)/(a+b+2): return bt*betacf(a,b,x)/a
    return 1.0 - bt*betacf(b,a,1-x)/b
def t_p2(t, df):
    x = df/(df+t*t); ib = betai(df/2, 0.5, x)
    return ib  # I_{df/(df+t²)}(df/2,1/2) 即双侧p，与t符号无关

report = ['\n## 补充验证：亚组漏斗图 + Meta回归（干预类型=异质性来源？）\n']
report.append('| 结局池 | k | Meta回归系数β(VR-AI) | SE | t | p(df=k-2) | 对应Qb检验p | 结论一致？ |')
report.append('|---|---|---|---|---|---|---|---|')

for pname, data in ma.pools.items():
    eff = ma.build_effects(data)
    items = [(v['grp'], v['g'], v['var']) for v in eff.values()]
    k = len(items)
    if k < 7:  # 满意度池k=5不画亚组漏斗、不做回归
        report.append(f'| {pname} | {k} | — | — | — | — | — | k<7，研究数不足，不做 |')
        continue
    g = np.array([i[1] for i in items]); v = np.array([i[2] for i in items])
    x = np.array([1.0 if i[0]=='VR' else 0.0 for i in items])
    # DL τ²（全池）
    w0 = 1/v
    gm0 = np.sum(w0*g)/np.sum(w0)
    Q = float(np.sum(w0*(g-gm0)**2))
    C = np.sum(w0) - np.sum(w0**2)/np.sum(w0)
    tau2 = max(0.0, (Q-(k-1))/C)
    # 混合效应Meta回归：g = b0 + b1*VR, 权重1/(v+τ²)
    w = 1/(v+tau2)
    X = np.column_stack([np.ones(k), x])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    beta = np.linalg.solve(XtWX, X.T @ W @ g)
    # K-H式方差校正
    resid = g - X@beta
    s2 = float(resid @ W @ resid/(k-2))
    covb = s2*np.linalg.inv(XtWX)
    se_b1 = math.sqrt(covb[1,1])
    t = beta[1]/se_b1
    p = t_p2(t, k-2)
    # 对应Qb
    ai = [(i[1],i[2]) for i in items if i[0]=='AI']; vr = [(i[1],i[2]) for i in items if i[0]=='VR']
    r1, r2 = ma.dl_pool(ai), ma.dl_pool(vr)
    qb, qp = ma.q_between(r1, r2, ai, vr)
    agree = '一致' if (p<0.05) == (qp<0.05) else '不一致'
    report.append(f'| {pname} | {k} | {beta[1]:.2f} | {se_b1:.2f} | {t:.2f} | {p:.4f} | {qp:.4f} | {agree} |')

    # 亚组漏斗图：左右两面板
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
    for ax, grp, color in [(axes[0], 'AI', '#1f77b4'), (axes[1], 'VR', '#d62728')]:
        sub = [(i[1], i[2]) for i in items if i[0] == grp]
        pool = ma.dl_pool(sub)
        for gi, vi in sub:
            ax.scatter(gi, math.sqrt(vi), s=45, color=color, zorder=3)
        se_max = max(math.sqrt(vi) for _, vi in sub)*1.15
        ses = np.linspace(0.001, se_max, 100)
        ax.plot(pool['g']+1.96*ses, ses, '--', color='#888', lw=1)
        ax.plot(pool['g']-1.96*ses, ses, '--', color='#888', lw=1)
        ax.axvline(pool['g'], color='k', lw=1)
        ax.set_ylim(se_max, 0)
        note = '（k≥10可参考）' if len(sub) >= 10 else '（k<10，仅描述性）'
        ax.set_title(f'{grp}亚组 k={len(sub)} {note}', fontsize=11)
        ax.set_xlabel("Hedges' g")
    axes[0].set_ylabel('标准误 SE')
    fig.suptitle(f'{pname}：AI/VR亚组漏斗图（分开评估，排除亚组混淆）', fontsize=12)
    plt.tight_layout()
    fig.savefig(OUT/f'亚组漏斗图_{pname}.png', bbox_inches='tight', dpi=150)
    plt.close(fig)

report.append('\n说明：单一二分类调节变量的Meta回归与Qb亚组检验数学等价，两者p值互为印证；')
report.append('k<10的亚组漏斗图/Egger检验功效不足，仅作描述性展示，论文中如实标注。')
txt = '\n'.join(report)
with open(OUT/'稳健性交叉验证报告.md', 'a', encoding='utf-8') as f:
    f.write('\n' + txt)
print(txt)
