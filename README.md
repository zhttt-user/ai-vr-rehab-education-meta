# 人工智能技术在康复教育中的应用效果与虚拟现实技术比较的Meta分析：分析代码与提取数据

本数据集为同名论文的关联数据，支撑论文全部定量结果的可重复验证。

更新记录：2026-09-12 经 Embase 补充检索新纳入 2 篇研究（Twose 2024、Stam 2023），纳入研究总数由 23 篇更新为 25 篇（AI 13 篇、VR 12 篇），理论知识池合并结果相应更新（总合并 k=21，g=0.86，95%CI：0.59～1.12），全部代码、数据与参考输出已同步更新。

## 目录结构

```
code/           分析代码（Python）
data/           提取数据与效应量结果表（CSV）
results/        参考输出：三份结果报告（md）与全部森林图、漏斗图、PRISMA 流程图（png）
```

## 运行环境

- Python ≥ 3.10
- 依赖见 requirements.txt：`pip install -r requirements.txt`

## 复现步骤

在 `code/` 目录下按下表顺序运行：

| 顺序 | 脚本 | 内容 |
|---|---|---|
| 1 | meta_analysis.py | 四个指标的总合并（不分技术类型）与 AI/VR 亚组合并、亚组差异检验（Qb）、敏感性分析 S1～S3、二分类满意率 OR 合并、森林图、各研究效应量表 |
| 2 | robustness_checks.py | Knapp-Hartung 校正、Peto 法交叉验证、Borenstein 相关谱系分析、Egger's 检验与剪补法、全池漏斗图 |
| 3 | outlier_sensitivity.py | 离群点识别与剔除敏感性分析（结果追加写入稳健性交叉验证报告） |
| 4 | leave_one_out.py | 逐篇剔除（leave-one-out）异质性分解 |
| 5 | egger_begg_full.py | 全池与亚组 Egger's / Begg's 检验汇总 |
| 6 | subgroup_funnel_metareg.py | 亚组漏斗图与 Meta 回归 |
| 7 | prisma_flowchart.py | PRISMA 文献筛选流程图 |

运行后全部输出写入 `results/` 目录（报告、图表及重新生成的 CSV），可与随包提供的参考输出逐一比对；`data/` 中的两份 CSV 即为论文附表所用的效应量数据。效应量算法与 R 语言 metafor 程序包一致。

## 许可

CC-BY 4.0（见 LICENSE.txt）。使用本数据集请注明出处。

## 引用格式

作者. 人工智能技术在康复教育中的应用效果与虚拟现实技术比较的Meta分析：分析代码与提取数据[DS/OL]. 科学数据银行（ScienceDB）, 2026. DOI: 10.57760/sciencedb.xxxxx（发布时填写）.
