# -*- coding: utf-8 -*-
"""稳健性交叉验证：回应5条方法学质疑
1. Knapp-Hartung校正 vs D-L（小研究数τ²不稳）
2. 零格子核查 + Peto OR（0.5校正质疑）
3. Borenstein r 谱系分析 r∈{0,0.25,0.5,0.75,0.9}
4. Qb多重比较核查（仅2亚组+Bonferroni）
5. Egger's检验 + 漏斗图 + 剪补法(L0)
"""
import sys, io, math
from pathlib import Path
sys.stdout = io.StringIO()  # 静默导入主脚本（复用其数据与函数）
import importlib.util
spec = importlib.util.spec_from_file_location(
    'ma', str(Path(__file__).resolve().parent / 'meta_analysis.py'))
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
sys.stdout = sys.__stdout__

import numpy as np
OUT = Path(__file__).resolve().parent.parent / 'results'

# ---------- t分布（无scipy，用不完全贝塔函数） ----------
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
def t_sf(t, df):
    """单侧 P(T>t)"""
    x = df/(df+t*t)
    ib = betai(df/2, 0.5, x)
    return ib/2 if t >= 0 else 1-ib/2
def t_p2(t, df):  # 双侧
    return 2*t_sf(abs(t), df)

# ---------- Knapp-Hartung ----------
def kh_pool(effects):
    """DL τ² + KH方差校正，返回 dict"""
    k = len(effects)
    g = np.array([e[0] for e in effects]); v = np.array([e[1] for e in effects])
    r = ma.dl_pool(effects)
    tau2 = r['tau2']
    w = 1/(v+tau2)
    gm = float(np.sum(w*g)/np.sum(w))
    q_resid = float(np.sum(w*(g-gm)**2))
    var_kh = q_resid/((k-1)*np.sum(w))
    se = math.sqrt(var_kh)
    t = gm/se
    p = t_p2(t, k-1)
    # t分布CI
    from math import sqrt
    # 近似t临界值用数值反解
    lo_t, hi_t = 0.0, 100.0
    for _ in range(100):
        mid = (lo_t+hi_t)/2
        if t_p2(mid, k-1) > 0.05: lo_t = mid
        else: hi_t = mid
    tcrit = (lo_t+hi_t)/2
    return dict(k=k, g=gm, se=se, lo=gm-tcrit*se, hi=gm+tcrit*se, t=t, p=p, tau2=tau2)

report = ['# 稳健性交叉验证报告（回应5条方法学质疑）\n']

# 各池效应量
pools_eff = {}
for pname, data in ma.pools.items():
    eff = ma.build_effects(data)
    pools_eff[pname] = eff

report.append('\n## 漏洞1验证：D-L vs Knapp-Hartung（各池各亚组）\n')
report.append('| 结局池 | 亚组 | k | D-L: g [95%CI], p | K-H: g [95%CI], p | 结论是否翻转 |')
report.append('|---|---|---|---|---|---|')
kh_rows = []
for pname, eff in pools_eff.items():
    for grp in ('AI','VR'):
        es = [(v['g'],v['var']) for v in eff.values() if v['grp']==grp]
        if len(es) < 2: continue
        dl = ma.dl_pool(es); kh = kh_pool(es)
        flip = ('显著→不显著' if (dl['p']<0.05 and kh['p']>=0.05)
                else '不显著→显著' if (dl['p']>=0.05 and kh['p']<0.05) else '一致')
        report.append(f"| {pname} | {grp} | {len(es)} | {dl['g']:.2f} [{dl['lo']:.2f},{dl['hi']:.2f}], p={dl['p']:.4f} | "
                      f"{kh['g']:.2f} [{kh['lo']:.2f},{kh['hi']:.2f}], p={kh['p']:.4f} | {flip} |")
        kh_rows.append((pname, grp, len(es), dl, kh, flip))

# 满意率OR的KH
es_or = [(e[2], e[3]) for e in ma.or_pool(ma.BINARY_SAT) if '敏感性替代' not in e[4]]
dl_or = ma.dl_pool(es_or); kh_or = kh_pool(es_or)
report.append(f"| 满意率(OR,log尺度) | 合并 | 3 | OR={math.exp(dl_or['g']):.2f}, p={dl_or['p']:.4f} | "
              f"OR={math.exp(kh_or['g']):.2f} [{math.exp(kh_or['lo']):.2f},{math.exp(kh_or['hi']):.2f}], p={kh_or['p']:.4f} | "
              f"{'显著→不显著' if kh_or['p']>=0.05 else '一致'} |")

# ---------- 漏洞2：零格子核查 + Peto ----------
report.append('\n## 漏洞2验证：零格子核查 + Peto OR交叉验证\n')
report.append('| 研究 | a | b | c | d | 最小格子 | 触发0.5校正？ | 对照组满意率 | Peto OR [95%CI] |')
report.append('|---|---|---|---|---|---|---|---|---|')
peto_oe, peto_v = 0.0, 0.0
for (name, grp, a, b, c, d, rob, tag) in ma.BINARY_SAT + ma.BINARY_EFFECT:
    N = a+b+c+d
    O_E = a - (a+b)*(a+c)/N
    V = (a+b)*(c+d)*(a+c)*(b+d)/(N*N*(N-1))
    lo = math.exp((O_E-1.96*math.sqrt(V))/V); hi = math.exp((O_E+1.96*math.sqrt(V))/V)
    peto = math.exp(O_E/V)
    zero = '是' if min(a,b,c,d)==0 else '否'
    if '敏感性替代' not in tag and name != '中文VR04 王婷':
        peto_oe += O_E; peto_v += V
    report.append(f"| {name} | {a} | {b} | {c} | {d} | {min(a,b,c,d)} | {zero} | {c/(c+d)*100:.0f}% | {peto:.2f} [{lo:.2f},{hi:.2f}] |")
peto_pool = math.exp(peto_oe/peto_v)
plo = math.exp((peto_oe-1.96*math.sqrt(peto_v))/peto_v); phi = math.exp((peto_oe+1.96*math.sqrt(peto_v))/peto_v)
report.append(f"\n满意率Peto合并: OR={peto_pool:.2f} [{plo:.2f},{phi:.2f}]（对比0.5校正D-L合并 OR=2.81 [1.28,6.17]）")

# ---------- 漏洞3：r谱系 ----------
report.append('\n## 漏洞3验证：Borenstein r 谱系分析（中文AI01高阶，对高阶AI池的影响）\n')
report.append('| r | AI01合成g | AI01合成var | 高阶AI池合并g [95%CI] | p | 结论 |')
report.append('|---|---|---|---|---|---|')
es3 = [ma.hedges_g(*x[:3], *x[3:]) for x in ma.AI01_HIGHER]
gs = np.array([e[0] for e in es3]); vs = np.array([e[1] for e in es3])
others = [(v['g'],v['var']) for k2,v in ma.build_effects(ma.HIGHER).items()
          if v['grp']=='AI' and k2 != '中文AI01 饶柳']
for r_ in (0, 0.25, 0.5, 0.75, 0.9):
    var_c = (vs.sum() + 2*r_*(math.sqrt(vs[0]*vs[1])+math.sqrt(vs[0]*vs[2])+math.sqrt(vs[1]*vs[2])))/9
    pool = ma.dl_pool(others + [(float(gs.mean()), var_c)])
    report.append(f"| {r_} | {gs.mean():.3f} | {var_c:.4f} | {pool['g']:.2f} [{pool['lo']:.2f},{pool['hi']:.2f}] | {pool['p']:.5f} | {'显著' if pool['p']<0.05 else '不显著'} |")

# ---------- 漏洞4：Qb多重比较 ----------
report.append('\n## 漏洞4核查：Qb多重比较\n')
report.append('本研究每个结局池仅2个亚组（AI vs VR），Qb为df=1的单一检验，不存在3组以上的多重比较问题。')
report.append('4个结局池共做4次Qb检验，Bonferroni校正后α=0.0125；理论知识Qb p=0.0246，未校正时显著，校正后不再显著，亚组差异结论应作探索性解读。')
report.append('按建议将亚组分析定位为"探索性"写入方法学限制。')

# ---------- 漏洞5：Egger + 漏斗图 + 剪补法 ----------
def eggers(effects):
    k = len(effects)
    g = np.array([e[0] for e in effects]); v = np.array([e[1] for e in effects])
    se = np.sqrt(v); z = g/se; prec = 1/se
    X = np.column_stack([np.ones(k), prec])
    beta = np.linalg.lstsq(X, z, rcond=None)[0]
    resid = z - X@beta
    s2 = float(resid@resid/(k-2))
    XtXinv = np.linalg.inv(X.T@X)
    se_a = math.sqrt(s2*XtXinv[0,0])
    t = beta[0]/se_a
    return dict(a=beta[0], se=se_a, t=t, p=t_p2(t, k-2), k=k)

def trimfill_L0(effects, side='right', maxit=100):
    """Duval-Tweedie 剪补法 L0估计量，返回(补后DL合并, k0)"""
    g = np.array([e[0] for e in effects]); v = np.array([e[1] for e in effects])
    k = len(g)
    def fe_pool(gg, vv):
        w = 1/vv; return float(np.sum(w*gg)/np.sum(w))
    theta = fe_pool(g, v)
    k0_prev = -1
    for _ in range(maxit):
        dev = g - theta
        ranks = np.empty(k); order = np.argsort(np.abs(dev)); ranks[order] = np.arange(1, k+1)
        if side == 'right':
            idx = np.where(dev > 0)[0]
        else:
            idx = np.where(dev < 0)[0]
        if len(idx) == 0: k0 = 0
        else:
            Tn = ranks[idx].sum()
            num = 4*np.sum(ranks[idx]**2) - k*Tn
            den = 2*Tn - k
            k0 = int(round(num/den)) if den > 0 else 0
            k0 = max(0, min(k0, k-2))
        if k0 == k0_prev: break
        k0_prev = k0
        if k0 > 0:
            trim_i = np.argsort(np.abs(dev))[::-1][:k0]
            keep = np.array([i for i in range(k) if i not in trim_i])
            theta = fe_pool(g[keep], v[keep])
        else:
            theta = fe_pool(g, v)
    if k0 == 0:
        return ma.dl_pool(effects), 0, theta
    dev = g - theta
    trim_i = np.argsort(np.abs(dev))[::-1][:k0]
    g_aug = list(g) + [2*theta - g[i] for i in trim_i]
    v_aug = list(v) + [v[i] for i in trim_i]
    return ma.dl_pool(list(zip(g_aug, v_aug))), k0, theta

report.append('\n## 漏洞5验证：Egger检验 + 剪补法\n')
report.append('| 池 | k | Egger截距 | t | p | 剪补法(右侧)k0 | 补后g [95%CI] | 原g | 变化 |')
report.append('|---|---|---|---|---|---|---|---|---|')
tf_results = {}
for pname, eff in pools_eff.items():
    for grp in ('AI','VR'):
        es = [(v['g'],v['var']) for v in eff.values() if v['grp']==grp]
        if len(es) < 5:
            report.append(f'| {pname}-{grp} | {len(es)} | — | — | — | — | — | — | k<5不适用 |')
            continue
        eg = eggers(es)
        orig = ma.dl_pool(es)
        # 两侧都试，取k0>0的一侧
        pool_r, k0r, _ = trimfill_L0(es, 'right')
        pool_l, k0l, _ = trimfill_L0(es, 'left')
        if k0r > 0:
            pool_tf, k0, sd = pool_r, k0r, '右'
        elif k0l > 0:
            pool_tf, k0, sd = pool_l, k0l, '左'
        else:
            pool_tf, k0, sd = orig, 0, '-'
        chg = abs(pool_tf['g']-orig['g'])/abs(orig['g'])*100 if orig['g'] != 0 else 0
        report.append(f"| {pname}-{grp} | {len(es)} | {eg['a']:.2f} | {eg['t']:.2f} | {eg['p']:.3f} | "
                      f"{k0}({sd}) | {pool_tf['g']:.2f} [{pool_tf['lo']:.2f},{pool_tf['hi']:.2f}] | {orig['g']:.2f} | {chg:.0f}% |")
        tf_results[(pname,grp)] = (eg, orig, pool_tf, k0)

# ---------- 漏斗图（k≥7的池） ----------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
for pname, eff in pools_eff.items():
    es_all = [(k2, v['grp'], v['g'], v['var']) for k2, v in eff.items()]
    if len(es_all) < 7: continue
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for name, grp, g, var in es_all:
        ax.scatter(g, math.sqrt(var), s=45, color='#1f77b4' if grp=='AI' else '#d62728', zorder=3)
    all_e = [(g, var) for _,_,g,var in es_all]
    pool = ma.dl_pool(all_e)
    se_max = max(math.sqrt(v) for _,_,_,v in es_all)*1.1
    ses = np.linspace(0.001, se_max, 100)
    ax.plot(pool['g']+1.96*ses, ses, '--', color='#888', lw=1)
    ax.plot(pool['g']-1.96*ses, ses, '--', color='#888', lw=1)
    ax.axvline(pool['g'], color='k', lw=1)
    ax.set_ylim(se_max, 0)
    ax.set_xlabel("Hedges' g"); ax.set_ylabel('标准误 SE')
    ax.set_title(f'{pname}漏斗图（蓝=AI，红=VR）')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color='#1f77b4',label='AI'),
                       Line2D([],[],marker='o',ls='',color='#d62728',label='VR')])
    plt.tight_layout()
    fig.savefig(OUT/f'漏斗图_{pname}.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
report.append('\n漏斗图已输出：漏斗图_理论知识.png / 漏斗图_技能操作.png / 漏斗图_高阶能力.png')

# ---------- 样本量汇总 ----------
report.append('\n## 各池研究数与总样本量\n')
report.append('| 池 | AI k | AI总n | VR k | VR总n | 合计k | 合计n |')
report.append('|---|---|---|---|---|---|---|')
for pname, eff in pools_eff.items():
    ai = [v for v in eff.values() if v['grp']=='AI']; vr = [v for v in eff.values() if v['grp']=='VR']
    report.append(f"| {pname} | {len(ai)} | {sum(v['n'] for v in ai)} | {len(vr)} | {sum(v['n'] for v in vr)} | {len(ai)+len(vr)} | {sum(v['n'] for v in ai)+sum(v['n'] for v in vr)} |")

txt = '\n'.join(report)
(OUT/'稳健性交叉验证报告.md').write_text(txt, encoding='utf-8')
print(txt)
