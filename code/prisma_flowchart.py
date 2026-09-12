# -*- coding: utf-8 -*-
"""PRISMA 文献筛选流程图（双支路：主检索 + Embase 补充检索）
2026-09-12 更新：右侧新增 Embase 补检支路（n=2077 → 纳入2篇）"""
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

OUT = Path(__file__).resolve().parent.parent / 'results' / 'PRISMA流程图.png'

fig, ax = plt.subplots(figsize=(11.5, 11))
ax.set_xlim(0, 14)
ax.set_ylim(0, 13)
ax.axis("off")

LX0, LX1 = 0.3, 4.5      # 左支路（主检索）
RX0, RX1 = 9.5, 13.7     # 右支路（Embase补检）
MX0, MX1 = 4.7, 9.3      # 中间排除框列
LC, RC, MC = (LX0+LX1)/2, (RX0+RX1)/2, (MX0+MX1)/2

def box(x0, x1, y0, y1, lines):
    ax.add_patch(plt.Rectangle((x0, y0), x1-x0, y1-y0, fill=False, lw=1.6, ec="black"))
    n = len(lines)
    ys = [y0 + (y1-y0)*(n-i)/(n+1) for i in range(n)]
    for (txt, bold, size), y in zip(lines, ys):
        ax.text((x0+x1)/2, y, txt, ha="center", va="center", fontsize=size,
                fontweight="bold" if bold else "normal")

def excl_box(x0, x1, y0, y1, title, items=None):
    ax.add_patch(plt.Rectangle((x0, y0), x1-x0, y1-y0, fill=False, lw=1.6, ec="black"))
    if items:
        ax.text((x0+x1)/2, y1-0.42, title, ha="center", va="center", fontsize=12)
        n = len(items)
        ys = [y1-0.95 - i*((y1-y0-1.25)/max(1, n-1)) for i in range(n)]
        for it, y in zip(items, ys):
            ax.text(x0+0.18, y, it, ha="left", va="center", fontsize=9.5)
    else:
        ax.text((x0+x1)/2, (y0+y1)/2, title, ha="center", va="center", fontsize=12)

def varrow(x, y_from, y_to):
    ax.add_patch(FancyArrowPatch((x, y_from), (x, y_to), arrowstyle="-|>",
                                 mutation_scale=22, lw=1.6, color="black"))

def elbow_to_excl(x_from, y_box, x_to):
    """支路框右边/左边 → 水平进入中间排除框"""
    ax.add_patch(FancyArrowPatch((x_from, y_box), (x_to, y_box), arrowstyle="-|>",
                                 mutation_scale=20, lw=1.5, color="black"))

# ---- 左支路：主检索 ----
box(LX0, LX1, 11.1, 12.7, [
    ("数据库检索获得相关文献", True, 13.5),
    ("(n=1 597)", True, 12.5),
    ("WoS核心合集841、PubMed 601", False, 10),
    ("中文库155（知网/万方/维普）", False, 10),
])
box(LX0, LX1, 8.5, 9.9, [
    ("合并去重，阅读文题和摘要初筛", True, 12),
    ("(n=1 597)", True, 12),
])
box(LX0, LX1, 5.7, 6.9, [
    ("阅读全文复筛", True, 12.5),
    ("(n=68)", True, 12),
])
# ---- 右支路：Embase 补检 ----
box(RX0, RX1, 11.1, 12.7, [
    ("补充检索Embase获得文献", True, 13.5),
    ("(n=2 077)", True, 12.5),
    ("（2026年9月补充检索）", False, 10),
])
box(RX0, RX1, 8.5, 9.9, [
    ("阅读文题和摘要初筛", True, 12),
    ("(n=2 077)", True, 12),
])
box(RX0, RX1, 5.7, 6.9, [
    ("阅读全文复筛", True, 12.5),
    ("(n=2)", True, 12),
])
# ---- 中间排除框 ----
excl_box(MX0, 6.9, 7.5, 8.1, "排除（n=1 529）")
excl_box(7.3, MX1, 7.5, 8.1, "排除（n=2 075）")
excl_box(MX0, 6.9, 2.75, 4.95, "排除（n=45）", [
    "● 不符合纳入排除",
    "　标准（n=34）",
    "● 研究对象或研究",
    "　类型不符（n=11）",
])
excl_box(7.3, MX1, 3.35, 4.45, "排除（n=0）")
# ---- 底部合并框 ----
box(LX0, RX1, 0.5, 2.1, [
    ("纳入Meta分析的文献 (n=25)", True, 14),
    ("中文13篇、英文12篇；AI组13篇、VR组12篇", False, 12),
])

# ---- 箭头 ----
varrow(LC, 11.1, 9.9)
varrow(LC, 8.5, 6.9)
varrow(RC, 11.1, 9.9)
varrow(RC, 8.5, 6.9)
# 支路 → 中间排除框
elbow_to_excl(LX1, 7.8, MX0)       # 左支路初筛排除
elbow_to_excl(RX0, 7.8, MX1)       # 右支路初筛排除
elbow_to_excl(LX1, 3.9, MX0)       # 左支路全文排除
elbow_to_excl(RX0, 3.9, MX1)       # 右支路全文排除
# 汇入底部合并框
varrow(LC, 5.7, 2.1)
varrow(RC, 5.7, 2.1)

fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("saved:", OUT)
