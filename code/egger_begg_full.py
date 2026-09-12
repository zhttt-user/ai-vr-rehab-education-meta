# -*- coding: utf-8 -*-
"""全池Egger's + Begg's检验（回应建议：亚组k<10时用全池检验发表偏倚）"""
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
def norm_p2(z):
    from math import erf
    return 2*(1-0.5*(1+erf(abs(z)/math.sqrt(2))))

def eggers(effects):
    k = len(effects)
    g = np.array([e[0] for e in effects]); v = np.array([e[1] for e in effects])
    se = np.sqrt(v); z = g/se; prec = 1/se
    X = np.column_stack([np.ones(k), prec])
    beta = np.linalg.lstsq(X, z, rcond=None)[0]
    resid = z - X@beta
    s2 = float(resid@resid/(k-2))
    se_a = math.sqrt(s2*np.linalg.inv(X.T@X)[0,0])
    return beta[0], beta[0]/se_a, t_p2(beta[0]/se_a, k-2)

def begg(effects):
    """Begg-Mazumdar秩相关检验（Kendall's tau，效应量 vs 方差）"""
    k = len(effects)
    g = [e[0] for e in effects]; v = [e[1] for e in effects]
    S = 0
    for i in range(k):
        for j in range(i+1, k):
            sg = int(g[i] > g[j]) - int(g[i] < g[j])
            sv = int(v[i] > v[j]) - int(v[i] < v[j])
            S += sg*sv
    var_s = k*(k-1)*(2*k+5)/18
    z = (abs(S)-1)/math.sqrt(var_s) * (1 if S > 0 else -1 if S < 0 else 0)
    tau = S/(k*(k-1)/2)
    return tau, z, norm_p2(z)

report = ['\n## 补充验证3：全池Egger\'s与Begg\'s检验（回应"亚组k<10时用全池"建议）\n']
report.append('| 结局池 | k | Egger截距 | Egger p | Begg tau | Begg z | Begg p | 结论 |')
report.append('|---|---|---|---|---|---|---|---|')
results = {}
for pname, data in ma.pools.items():
    eff = ma.build_effects(data)
    es = [(v['g'], v['var']) for v in eff.values()]
    if len(es) < 10:
        report.append(f'| {pname} | {len(es)} | — | — | — | — | — | k<10，检验功效不足，不做 |')
        continue
    a, t, pe = eggers(es)
    tau, z, pb = begg(es)
    concl = '无明显发表偏倚证据' if pe > 0.1 and pb > 0.1 else '存在不对称，需结合亚组混淆解释'
    report.append(f'| {pname} | {len(es)} | {a:.2f} | {pe:.3f} | {tau:.2f} | {z:.2f} | {pb:.3f} | {concl} |')
    results[pname] = (a, pe, tau, z, pb)
report.append('\n注：全池检验将AI/VR合并，若理论池检出不对称，应结合Meta回归（干预类型β=0.58, P=0.049）')
report.append('判断其为亚组效应差异所致的"假性不对称"，而非发表偏倚。亚组内Egger结果见前文表7。')
txt = '\n'.join(report)
with open(OUT/'稳健性交叉验证报告.md', 'a', encoding='utf-8') as f:
    f.write('\n' + txt)
print(txt)
