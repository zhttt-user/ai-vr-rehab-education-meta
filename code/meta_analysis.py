# -*- coding: utf-8 -*-
"""AI/VR赋能康复教育Meta分析：四大结局池 + 二分类满意率/教学有效率
方法：Hedges' g (SMD) + DerSimonian-Laird随机效应；OR(二分类) + DL随机效应
"""
import sys, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

OUT = Path(__file__).resolve().parent.parent / 'results'
OUT.mkdir(exist_ok=True)

# ============ 数据（来源：Meta数据提取表_V2_23篇.xlsx 定稿版）============
# 字段: 研究, 组(AI/VR), ne, me, sde, nc, mc, sdc, RoB总体, 设计, 标记
THEORY = [
 ('中文AI01 饶柳','AI',60,83.42,3.54,60,81.32,3.27,'有一些担忧','RCT',''),
 ('中文AI02 万丽丽','AI',46,17.63,1.85,45,15.97,2.46,'有一些担忧','RCT','合成'),
 ('中文AI03 李梦晓','AI',30,77.37,7.66,30,74.17,7.62,'有一些担忧','RCT',''),
 ('中文AI04 朱丹','AI',30,91.45,3.20,30,84.23,4.87,'有一些担忧','RCT','AI+BOPPPS'),
 ('中文AI05 张璐','AI',36,81.74,3.32,36,78.75,4.11,'有一些担忧','RCT',''),
 ('中文AI06 张佳钰','AI',25,79.40,6.60,87,73.80,8.90,'严重风险','准实验','S1剔除'),
 ('中文AI07 陈鹏鑫','AI',35,93.97,4.53,36,93.22,4.04,'高风险','准实验','S1剔除/阴性'),
 ('英文AI03 Sigirtmac','AI',11,4.91,1.28,14,5.54,1.22,'有一些担忧','RCT','阴性'),
 ('英文AI04 Ergezen','AI',19,24.20,11.50,16,22.60,11.50,'有一些担忧','RCT','差值/SD近似/S2剔除'),
 ('英文AI06 Okuno','AI',92,107.67,45.93,77,80.00,44.44,'中等风险','准实验','Wan换算/S2剔除'),
 ('中文VR01 李鹏','VR',33,89.76,2.96,33,86.61,3.64,'有一些担忧','RCT',''),
 ('中文VR02 辜茜娟','VR',30,86.80,4.18,30,80.63,5.06,'有一些担忧','RCT',''),
 ('中文VR03 王树本','VR',69,16.57,2.13,69,10.87,3.32,'有一些担忧','RCT',''),
 ('中文VR04 王婷','VR',34,9.12,1.34,34,6.98,1.23,'有一些担忧','RCT',''),
 ('中文VR05 余美华','VR',50,85.71,6.14,50,80.76,8.63,'有一些担忧','RCT',''),
 ('中文VR06 张通','VR',14,89.07,1.05,12,87.13,0.87,'有一些担忧','RCT','阴性'),
 ('英文VR01 García','VR',34,11.40,2.60,33,6.50,2.60,'有一些担忧','RCT',''),
 ('英文VR02 He','VR',20,88.00,5.62,20,82.60,7.49,'有一些担忧','RCT','A vs D'),
 ('英文VR03 Kurul','VR',36,33.26,22.86,36,10.33,10.13,'有一些担忧','RCT','差值'),
 ('英文VR05 Twose','VR',18,22.17,2.96,20,22.00,1.48,'有一些担忧','RCT','Wan换算/非劣效试点/S2剔除'),
 ('英文VR06 Stam','VR',9,72.50,16.67,9,75.83,12.96,'有一些担忧','RCT','Wan换算/交叉第一阶段/S2剔除'),
]
SKILL = [
 ('中文AI02 万丽丽','AI',46,86.90,4.97,45,80.83,5.50,'有一些担忧','RCT','合成'),
 ('中文AI03 李梦晓','AI',30,76.91,2.29,30,73.18,2.40,'有一些担忧','RCT',''),
 ('中文AI04 朱丹','AI',30,93.27,2.85,30,85.17,4.95,'有一些担忧','RCT','OSCE'),
 ('中文AI05 张璐','AI',36,89.16,3.45,36,84.42,3.26,'有一些担忧','RCT',''),
 ('中文AI07 陈鹏鑫','AI',35,87.65,4.18,36,88.22,4.70,'高风险','准实验','S1剔除/阴性'),
 ('英文AI03 Sigirtmac','AI',11,5.57,1.51,14,6.23,1.08,'有一些担忧','RCT','阴性'),
 ('英文AI04 Ergezen','AI',19,38.00,13.10,16,31.30,18.50,'有一些担忧','RCT','Mini-CEX阴性'),
 ('中文VR01 李鹏','VR',33,90.42,2.98,33,87.58,2.92,'有一些担忧','RCT',''),
 ('中文VR02 辜茜娟','VR',30,81.80,5.24,30,75.37,6.63,'有一些担忧','RCT',''),
 ('中文VR03 王树本','VR',69,53.24,2.36,69,44.02,4.16,'有一些担忧','RCT',''),
 ('中文VR04 王婷','VR',34,8.97,1.25,34,7.23,1.33,'有一些担忧','RCT',''),
 ('中文VR05 余美华','VR',50,90.79,6.20,50,87.45,6.43,'有一些担忧','RCT',''),
 ('中文VR06 张通','VR',14,91.32,0.67,12,89.45,0.12,'有一些担忧','RCT','SD极小'),
 ('英文VR02 He','VR',20,82.30,5.92,20,77.70,6.03,'有一些担忧','RCT','混合分/S3剔除'),
 ('英文VR04 Bonnin','VR',23,9.17,1.85,25,11.33,3.70,'有一些担忧','RCT','Wan/阴性/S2剔除'),
]
HIGHER = [
 ('中文AI01 饶柳','AI',60,None,None,60,None,None,'有一些担忧','RCT','Borenstein合成'),
 ('中文AI02 万丽丽','AI',46,22.83,2.20,45,20.59,2.46,'有一些担忧','RCT','合成'),
 ('中文AI04 朱丹','AI',30,86.88,5.15,30,75.77,8.75,'有一些担忧','RCT',''),
 ('中文AI06 张佳钰','AI',25,25.20,2.00,87,22.70,3.70,'严重风险','准实验','S1剔除'),
 ('中文AI07 陈鹏鑫','AI',35,253.00,22.89,36,237.83,23.60,'高风险','准实验','S1剔除'),
 ('英文AI02 Yildiz','AI',30,6.70,4.62,30,4.33,2.43,'有一些担忧','RCT','差值/阴性'),
 ('英文AI03 Sigirtmac','AI',11,15.52,3.42,14,17.04,2.78,'有一些担忧','RCT','阴性'),
 ('英文AI04 Ergezen','AI',19,19.80,11.70,16,-3.90,15.00,'有一些担忧','RCT','差值'),
 ('中文VR01 李鹏','VR',33,86.82,2.10,33,85.30,2.30,'有一些担忧','RCT',''),
 ('中文VR02 辜茜娟','VR',30,86.13,3.82,30,81.79,4.64,'有一些担忧','RCT','合成'),
 ('中文VR03 王树本','VR',69,9.15,0.30,69,8.18,0.54,'有一些担忧','RCT','合成'),
]
SATIS = [
 ('中文AI05 张璐','AI',36,8.26,0.73,36,7.32,0.91,'有一些担忧','RCT','合成'),
 ('中文AI07 陈鹏鑫','AI',35,76.31,7.51,36,73.19,9.99,'高风险','准实验','S1剔除/阴性'),
 ('中文VR01 李鹏','VR',33,16.92,0.94,33,15.73,1.01,'有一些担忧','RCT','合成'),
 ('中文VR05 余美华','VR',50,3.76,0.77,50,3.28,1.04,'有一些担忧','RCT','合成'),
 ('英文VR01 García','VR',34,9.60,0.70,33,3.40,2.50,'有一些担忧','RCT',''),
]
# 中文AI01 高阶 3量表（Borenstein合成）
AI01_HIGHER = [
 (60,93.45,4.27,60,91.00,5.15),   # 深度学习
 (60,264.47,15.10,60,255.60,11.15), # 批判性思维CTDI
 (60,76.20,3.19,60,74.12,2.42),   # 积极主动性
]
# 二分类：满意率（a=满意 b=不满意 | c=满意 d=不满意）
BINARY_SAT = [
 ('中文AI03 李梦晓(条目合计)','AI',256,14,246,24,'有一些担忧','n放大注意'),
 ('中文AI03 李梦晓(授课方式条目)','AI',29,1,28,2,'有一些担忧','敏感性替代'),
 ('中文VR03 王树本','VR',65,4,56,13,'有一些担忧',''),
 ('中文VR06 张通','VR',12,2,5,7,'有一些担忧','喜欢该学习方法'),
]
BINARY_EFFECT = [  # 教学有效率（高阶能力二分类）
 ('中文VR04 王婷','VR',32,2,25,9,'有一些担忧','总有效率'),
]

# ============ 效应量计算 ============
def hedges_g(ne, me, sde, nc, mc, sdc):
    sp = math.sqrt(((ne-1)*sde**2 + (nc-1)*sdc**2) / (ne+nc-2))
    d = (me - mc) / sp
    df = ne + nc - 2
    J = 1 - 3/(4*df - 1)
    g = J * d
    var = (ne+nc)/(ne*nc) + g**2/(2*(ne+nc-2))
    return g, var

def dl_pool(effects):
    """effects: list of (g, var) → DL随机效应合并"""
    k = len(effects)
    g = np.array([e[0] for e in effects]); v = np.array([e[1] for e in effects])
    w = 1/v
    g_fe = np.sum(w*g)/np.sum(w)
    Q = float(np.sum(w*(g-g_fe)**2))
    df = k-1
    C = np.sum(w) - np.sum(w**2)/np.sum(w)
    tau2 = max(0.0, (Q-df)/C) if C>0 else 0.0
    I2 = max(0.0, (Q-df)/Q)*100 if Q>0 else 0.0
    wr = 1/(v+tau2)
    g_re = float(np.sum(wr*g)/np.sum(wr))
    se = math.sqrt(1/np.sum(wr))
    z = g_re/se; 
    from math import erf
    p = 2*(1-0.5*(1+erf(abs(z)/math.sqrt(2))))
    return dict(k=k, g=g_re, se=se, lo=g_re-1.96*se, hi=g_re+1.96*se,
                Q=Q, df=df, I2=I2, tau2=tau2, z=z, p=p)

def q_between(res1, res2, eff1, eff2):
    """亚组间差异检验(Qb)：基于两亚组各自DL合并"""
    allw = None
    g1, g2 = res1['g'], res2['g']
    se1, se2 = res1['se'], res2['se']
    # 用固定效应框架的Qb近似
    g_all = np.array([e[0] for e in eff1]+[e[0] for e in eff2])
    v_all = np.array([e[1] for e in eff1]+[e[1] for e in eff2])
    tau_all = (res1['tau2']*res1['df']+res2['tau2']*res2['df'])/max(1,(res1['df']+res2['df']))
    w = 1/(v_all+tau_all)
    gm = np.sum(w*g_all)/np.sum(w)
    Qb = float(w[:len(eff1)].sum()*(g1-gm)**2 + w[len(eff1):].sum()*(g2-gm)**2)
    from math import erf
    # chi2 df=1 的p值
    p = 1-0.5*(1+erf(math.sqrt(Qb/2)))*2/2  # 简化: chi2_1 sf = 2*(1-Phi(sqrt(Qb)))
    p = 2*(1-0.5*(1+erf(math.sqrt(Qb)/math.sqrt(2))))
    return Qb, p

# ============ 主分析 ============
def build_effects(data):
    out = {}
    for (name, grp, ne, me, sde, nc, mc, sdc, rob, design, tag) in data:
        if me is None:  # 中文AI01 Borenstein
            es = [hedges_g(*x[:3], *x[3:]) for x in AI01_HIGHER]
            gs = np.array([e[0] for e in es]); vs = np.array([e[1] for e in es])
            r = 0.5
            g_c = gs.mean()
            var_c = (vs.sum() + 2*r*(math.sqrt(vs[0]*vs[1])+math.sqrt(vs[0]*vs[2])+math.sqrt(vs[1]*vs[2])))/9
            out[name] = dict(grp=grp, g=g_c, var=var_c, n=120, rob=rob, design=design, tag=tag)
        else:
            g, var = hedges_g(ne, me, sde, nc, mc, sdc)
            out[name] = dict(grp=grp, g=g, var=var, n=ne+nc, rob=rob, design=design, tag=tag)
    return out

def analyze_pool(pname, data, exclusions=()):
    eff = build_effects(data)
    eff = {k:v for k,v in eff.items() if k not in exclusions}
    ai = [(v['g'],v['var']) for v in eff.values() if v['grp']=='AI']
    vr = [(v['g'],v['var']) for v in eff.values() if v['grp']=='VR']
    res = {}
    res['AI'] = dl_pool(ai) if ai else None
    res['VR'] = dl_pool(vr) if vr else None
    res['ALL'] = dl_pool(ai + vr) if (ai or vr) else None  # 总合并（不分技术类型）
    if ai and vr and len(ai)>1 and len(vr)>1:
        res['Qb'] = q_between(res['AI'], res['VR'], ai, vr)
    return eff, res

def fmt(res):
    if not res: return '—'
    sig = '*' if res['p']<0.05 else ''
    return (f"k={res['k']}, g={res['g']:.2f} [{res['lo']:.2f},{res['hi']:.2f}]{sig}, "
            f"I²={res['I2']:.0f}%, τ²={res['tau2']:.2f}, p={res['p']:.4f}")

pools = {'理论知识':THEORY, '技能操作':SKILL, '高阶能力':HIGHER, '满意度(连续)':SATIS}
report = []
all_eff = {}
for pname, data in pools.items():
    eff, res = analyze_pool(pname, data)
    all_eff[pname] = eff
    report.append(f"\n### {pname}（主分析）")
    report.append(f"  总合并(不分技术): {fmt(res['ALL'])}")
    report.append(f"  AI组: {fmt(res['AI'])}")
    report.append(f"  VR组: {fmt(res['VR'])}")
    if 'Qb' in res:
        report.append(f"  亚组差异 Qb={res['Qb'][0]:.2f}, p={res['Qb'][1]:.4f}")

# ============ 敏感性分析 ============
S1 = ['中文AI06 张佳钰','中文AI07 陈鹏鑫']  # 严重+高风险
S2 = ['英文AI04 Ergezen','英文AI06 Okuno','英文VR04 Bonnin','英文VR05 Twose','英文VR06 Stam']  # 数据形式近似/Wan
S3 = ['英文VR02 He']  # 混合构念
report.append('\n\n### 敏感性分析')
for pname, data in pools.items():
    for sname, excl in [('S1(剔除严重/高RoB)', S1), ('S2(剔除换算/近似)', S2), ('S3(剔除混合构念)', S3)]:
        excl_use = [e for e in excl if any(d[0]==e for d in data)]
        if not excl_use: continue
        eff2 = build_effects(data)
        eff2 = {k:v for k,v in eff2.items() if k not in excl_use}
        ai = [(v['g'],v['var']) for v in eff2.values() if v['grp']=='AI']
        vr = [(v['g'],v['var']) for v in eff2.values() if v['grp']=='VR']
        line = f"\n  {pname}-{sname}（剔除{len(excl_use)}篇）:"
        if ai: line += f"\n    AI组: {fmt(dl_pool(ai))}"
        if vr: line += f"\n    VR组: {fmt(dl_pool(vr))}"
        report.append(line)

# ============ 二分类 OR ============
def or_pool(data):
    effects = []
    for (name, grp, a, b, c, d, rob, tag) in data:
        aa,bb,cc,dd = [x+0.5 if min(a,b,c,d)==0 else x for x in (a,b,c,d)] if min(a,b,c,d)==0 else (a,b,c,d)
        lor = math.log((aa*dd)/(bb*cc))
        var = 1/aa+1/bb+1/cc+1/dd
        effects.append((name, grp, lor, var, tag))
    return effects

def dl_pool_or(effects):
    es = [(e[2], e[3]) for e in effects]
    r = dl_pool(es)
    r['OR'], r['ORlo'], r['ORhi'] = math.exp(r['g']), math.exp(r['lo']), math.exp(r['hi'])
    return r

report.append('\n\n### 满意率（二分类，OR）')
eff_sat = or_pool(BINARY_SAT)
for (name, grp, lor, var, tag) in eff_sat:
    report.append(f"  {name}: OR={math.exp(lor):.2f} [{math.exp(lor-1.96*math.sqrt(var)):.2f},{math.exp(lor+1.96*math.sqrt(var)):.2f}] {tag}")
main_sat = [e for e in eff_sat if '敏感性替代' not in e[4]]
r_sat = dl_pool_or(main_sat)
report.append(f"  合并(3篇): OR={r_sat['OR']:.2f} [{r_sat['ORlo']:.2f},{r_sat['ORhi']:.2f}], I²={r_sat['I2']:.0f}%, p={r_sat['p']:.4f}")
sens_sat = [e for e in eff_sat if '条目合计' not in e[4]]
r_sat2 = dl_pool_or(sens_sat)
report.append(f"  敏感性(李梦晓换单条目): OR={r_sat2['OR']:.2f} [{r_sat2['ORlo']:.2f},{r_sat2['ORhi']:.2f}], I²={r_sat2['I2']:.0f}%, p={r_sat2['p']:.4f}")

report.append('\n### 教学有效率（二分类，OR，单篇）')
for (name, grp, a, b, c, d, rob, tag) in BINARY_EFFECT:
    lor = math.log((a*d)/(b*c)); var = 1/a+1/b+1/c+1/d
    report.append(f"  {name}: OR={math.exp(lor):.2f} [{math.exp(lor-1.96*math.sqrt(var)):.2f},{math.exp(lor+1.96*math.sqrt(var)):.2f}] {tag}")

# ============ 森林图 ============
def forest(pname, eff, res, fname):
    studies = list(eff.items())
    studies.sort(key=lambda kv: (kv[1]['grp']!='AI', -kv[1]['g']))
    n = len(studies)
    fig, ax = plt.subplots(figsize=(10, 0.55*n+2.6))
    y = 0; ylabels = []; ypos = []
    for name, e in studies:
        se = math.sqrt(e['var'])
        lo, hi = e['g']-1.96*se, e['g']+1.96*se
        ax.plot([lo,hi],[y,y],'-',color='#555',lw=1.2)
        ax.plot(e['g'],y,'s',color='#1f77b4' if e['grp']=='AI' else '#d62728', ms=7)
        tag = f" [{e['tag']}]" if e['tag'] else ''
        ylabels.append(f"{name} ({e['grp']}){tag}")
        ypos.append(y); y += 1
    # 亚组合并菱形
    for grp, color in [('AI','#1f77b4'),('VR','#d62728')]:
        r = res[grp]
        if not r: continue
        y += 0.4
        cx = r['g']; half_lo = cx-r['lo']; half_hi = r['hi']-cx
        ax.fill([r['lo'],cx,r['hi'],cx],[y,y+0.28,y,y-0.28],color=color,alpha=0.75)
        ylabels.append(f"◆ {grp}组合并 (I²={r['I2']:.0f}%, p={r['p']:.3f})")
        ypos.append(y); y += 1.2
    ax.axvline(0,color='k',lw=0.8,ls='--')
    ax.set_yticks(ypos); ax.set_yticklabels(ylabels,fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Hedges' g（SMD，>0 偏向技术组）")
    ax.set_title(f"{pname}：AI组 vs VR组 随机效应合并")
    plt.tight_layout()
    fig.savefig(OUT/fname, bbox_inches='tight', dpi=150)
    plt.close(fig)

for pname, data in pools.items():
    eff = all_eff[pname]
    ai = [(v['g'],v['var']) for v in eff.values() if v['grp']=='AI']
    vr = [(v['g'],v['var']) for v in eff.values() if v['grp']=='VR']
    res = {'AI': dl_pool(ai) if ai else None, 'VR': dl_pool(vr) if vr else None}
    forest(pname, eff, res, f"森林图_{pname.replace('(','').replace(')','')}.png")

# 满意率森林图
fig, ax = plt.subplots(figsize=(9,3.6))
ypos=[]; ylabels=[]
for i,(name, grp, lor, var, tag) in enumerate(eff_sat):
    lo,hi = lor-1.96*math.sqrt(var), lor+1.96*math.sqrt(var)
    ax.plot([math.exp(lo),math.exp(hi)],[i,i],'-',color='#555',lw=1.2)
    ax.plot(math.exp(lor),i,'s',color='#2ca02c',ms=7)
    ypos.append(i); ylabels.append(f"{name} ({grp}) {tag}")
i=len(eff_sat)+0.5
ax.fill([r_sat['ORlo'],r_sat['OR'],r_sat['ORhi'],r_sat['OR']],[i,i+0.28,i,i-0.28],color='#2ca02c',alpha=0.75)
ypos.append(i); ylabels.append(f"◆ 合并 (I²={r_sat['I2']:.0f}%, p={r_sat['p']:.3f})")
ax.axvline(1,color='k',lw=0.8,ls='--')
ax.set_yticks(ypos); ax.set_yticklabels(ylabels,fontsize=9); ax.invert_yaxis()
ax.set_xscale('log'); ax.set_xlabel('OR（>1 偏向技术组，log刻度）')
ax.set_title('教学满意率（二分类）随机效应合并')
plt.tight_layout(); fig.savefig(OUT/'森林图_满意率OR.png', bbox_inches='tight', dpi=150); plt.close(fig)

# 逐研究效应量表（含样本量n与各亚组内随机效应权重%）
rows_out = []
for pname, eff in all_eff.items():
    # 先按亚组算DL合并，得到该亚组τ²，再求每篇权重
    wt = {}
    for grp in ('AI','VR'):
        sub = {k:v for k,v in eff.items() if v['grp']==grp}
        if not sub: continue
        r = dl_pool([(v['g'],v['var']) for v in sub.values()])
        wr = {k: 1/(v['var']+r['tau2']) for k,v in sub.items()}
        tot = sum(wr.values())
        for k,w in wr.items():
            wt[k] = w/tot*100
    for name, e in eff.items():
        se = math.sqrt(e['var'])
        rows_out.append(dict(结局池=pname, 研究=name, 组=e['grp'], n=e['n'], g=round(e['g'],3),
                             SE=round(se,3), CI_low=round(e['g']-1.96*se,3), CI_high=round(e['g']+1.96*se,3),
                             权重pct=round(wt.get(name,0),1),
                             RoB=e['rob'], 设计=e['design'], 标记=e['tag']))
try:
    pd.DataFrame(rows_out).to_csv(OUT/'各研究效应量_SMD.csv', index=False, encoding='utf-8-sig')
except PermissionError:
    pd.DataFrame(rows_out).to_csv(OUT/'各研究效应量_SMD_含权重.csv', index=False, encoding='utf-8-sig')
    print('原CSV被占用，已另存为 各研究效应量_SMD_含权重.csv')

(OUT/'合并结果报告.md').write_text('\n'.join(report), encoding='utf-8')
print('\n'.join(report))
print('\n输出文件:', [f.name for f in OUT.iterdir()])
