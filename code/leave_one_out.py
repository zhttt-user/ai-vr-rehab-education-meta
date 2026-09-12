# -*- coding: utf-8 -*-
"""逐篇排除(leave-one-out)异质性分解分析
目的:对I²>50%的各结局亚组池,逐篇剔除后重合并,定位异质性来源是
"单篇驱动"还是"弥漫性"。数据复用 meta_analysis.py 顶部唯一数据源。
输出: 逐篇排除_LOO结果.csv + 逐篇排除_LOO报告.md
"""
import sys, io
from pathlib import Path
sys.stdout = io.StringIO()  # 静默导入主脚本
import importlib.util
spec = importlib.util.spec_from_file_location(
    'ma', str(Path(__file__).resolve().parent / 'meta_analysis.py'))
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
sys.stdout = sys.__stdout__

import csv
OUT = Path(__file__).resolve().parent.parent / 'results'

pools = {'理论知识': ma.THEORY, '技能操作': ma.SKILL, '高阶能力': ma.HIGHER, '满意度(连续)': ma.SATIS}

rows = []
for pname, data in pools.items():
    eff = ma.build_effects(data)
    for grp in ['AI', 'VR']:
        members = {k: v for k, v in eff.items() if v['grp'] == grp}
        if len(members) < 3:
            continue
        full = ma.dl_pool([(v['g'], v['var']) for v in members.values()])
        # 逐篇剔除
        loo = []
        for name, v in members.items():
            rest = [(x['g'], x['var']) for k2, x in members.items() if k2 != name]
            r = ma.dl_pool(rest)
            loo.append((name, r['I2'], r['g'], full['I2'] - r['I2']))
        loo.sort(key=lambda x: -x[3])
        # 判定性质:最大单篇剔除后I²仍>50% → 弥漫性;否则单篇驱动
        nature = '单篇驱动' if loo[0][1] <= 50 else '弥漫性'
        for rank, (name, i2_new, g_new, drop) in enumerate(loo, 1):
            rows.append([pname, grp, len(members), f"{full['I2']:.1f}", f"{full['g']:.2f}",
                         rank, name, f"{i2_new:.1f}", f"{drop:+.1f}", f"{g_new:.2f}", nature])

with open(OUT / '逐篇排除_LOO结果.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['结局池', '亚组', 'k', '全池I2%', '全池g', '剔除排序', '剔除研究', '剔除后I2%', 'I2变化(点)', '剔除后g', '异质性性质'])
    w.writerows(rows)

# 汇总报告(每池只列前3名)
lines = ['# 逐篇排除(Leave-one-out)异质性分解报告', '',
         '方法:对k≥3的各结局亚组池,逐篇剔除后按D-L随机效应模型重合并,记录I²与g的变化。',
         '判定:最大单篇剔除后I²≤50%为"单篇驱动",仍>50%为"弥漫性"。', '']
cur = None
for r in rows:
    key = (r[0], r[1])
    if key != cur:
        cur = key
        lines.append(f"## {r[0]}-{r[1]} (k={r[2]}, 全池I²={r[3]}%, g={r[4]}) → {r[10]}")
        lines.append('')
        lines.append('| 排序 | 剔除研究 | 剔除后I²% | I²变化 | 剔除后g |')
        lines.append('|---|---|---|---|---|')
    if r[5] <= 3:
        lines.append(f'| {r[5]} | {r[6]} | {r[7]} | {r[8]} | {r[9]} |')
    if r[5] == 3:
        lines.append('')
(OUT / '逐篇排除_LOO报告.md').write_text('\n'.join(lines), encoding='utf-8')

# 控制台摘要
print('逐篇排除LOO完成,各池性质判定:')
seen = set()
for r in rows:
    key = (r[0], r[1])
    if key not in seen:
        seen.add(key)
        top = [x for x in rows if (x[0], x[1]) == key][0]
        print(f"  {r[0]}-{r[1]}: 全池I²={r[3]}% → 剔「{top[6]}」后I²={top[7]}% ({top[8]}点) → {r[10]}")
