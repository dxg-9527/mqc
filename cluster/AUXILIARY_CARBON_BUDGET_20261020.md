# 辅助供给修复和集群补跑

本次修复以 `d7a49ce319edc9e63d50b693d4d88737f7f148e8` 为基线。原生物量回补开放的辅助交换和输入型 sink，纳入 MQC 已有的 pFBA 比例碳预算。限制写入反应边界，导出 XML 重读后再验证；循环罚分、回补选路和零输入能量通过标准不变。

没有使用实验生长率设定边界。这里调整的是 MQC 的标准化供给条件，不是宣称重建了真实培养基，也不是给生长率设置人为上限。

## 需要补跑哪些模型

已读取上次 6543 个模型队列中所有已完成模型的原始辅助供给记录，结合最终交换方向确定范围。不依赖模型名称或生长高低筛选。

| 方案 | ModelSEED | VMH | BiGG | 总数 |
|---|---:|---:|---:|---:|
| 命中辅助供给检查的已完成模型 | 74 | 328 | 1 | 403 |
| 两个数据库整体补跑并补上 BiGG 例外 | 79 | 806 | 1 | 886 |

默认建议按 403 个的清单补跑，可以省去未命中此次修复分支的模型。若希望两个数据库的全部记录都重新生成，则采用 886 个方案；其中包括上次未完成的 20 个模型，不保证此次供给修复能解决其原有失败原因。

BiGG 例外是 `bigg_005614 / iIS312_Amastigote`。局部验证中生长值未变化，但供给边界有所收紧，因此纳入清单。5580 个已完成的 CarveMe 模型未启用此项辅助含碳供给，无需因本次修复再批量重跑；其原有 5 个失败也不因本次局部修复自动变成通过。

两个清单已经包含在仓库：

- `cluster/rerun_auxiliary_affected_20261020.csv`：403 个，默认推荐。
- `cluster/rerun_auxiliary_databases_20261020.csv`：886 个，可选的整库方案。

只选一个，不要分别提交后重复运行交集。继续使用原 `mqc_inputs_6543` 模型包和已验证环境，不需要重新传输输入模型。

## 从 GitHub 更新并提交

在集群源码目录确认没有需要保留的未提交改动，然后更新优化分支；若工作区有本地修改，先保留它们，不使用强制覆盖命令。

```bash
git fetch origin
git switch optimized-20261018
git pull --ff-only origin optimized-20261018
git rev-parse HEAD
git status --short
sha256sum cluster/release_lock.json
```

以下 `INPUTS_ABS`、`NEW_RESULTS_ABS`、`HANDOFF_ABS` 均须替换成集群实际绝对路径。**新建结果目录，不与上一版结果混写。**

```bash
python cluster/check_environment.py --output HANDOFF_ABS/environment_before_auxiliary.json
python cluster/verify_inputs.py --inputs INPUTS_ABS

python cluster/submit.py \
  --inputs INPUTS_ABS --results NEW_RESULTS_ABS \
  --ids-file cluster/rerun_auxiliary_affected_20261020.csv \
  --partition qcpu_23i --parallel 32 --cpus 10 --mem 8G \
  --time 24:00:00 --model-timeout 7200
```

先检查计划显示 `Models: 403`，再用同一命令追加 `--submit`。若选择整库方案，仅替换 `--ids-file` 为另一个文件，计划应显示 `Models: 886`。分区、账户、并发额度沿用集群已跑通的配置；整批运行期间不要更新源码或环境。

## 验收和回传

除了原有输入身份、XML 重读、生长、pFBA 收支和零输入能量检查，还要保留：

- `control_result.json`、`model_info.json` 中的 `auxiliary_supply_budget`。
- 每条被收紧输入的原边界、新边界及 pFBA 分配依据。
- 分子式覆盖标记。带 `R/X` 的公式只能按明确列出的碳原子计算名义预算，未知基团不能解释为已知完整组成。
- 人工输入型 sink 的依赖，不能把人工供给称为已经补通的生物合成途径。

`partial_formula_inputs` 列出带未定基团但有明确碳数的输入；`unclassified_existing_inputs` 列出原有记账或未分类边界，不能把它们默认为已知的零碳化合物。本轮不擅自删除原有 `X` 记账边界；新开放的辅助入口若连已知碳数都无法确定，则不能作为已完成的预算分配放行。`carbon_accounting` 明确区分完整输入公式与仅按明确原子计算的预算，后者不代表完整元素收支已获证明。

所有任务停止后：

```bash
python cluster/collect.py --inputs INPUTS_ABS --results NEW_RESULTS_ABS --pack
```

收集器仍以 6543 个输入汇总，所以没有提交的模型会显示未运行。**这不是本轮失败清单：403 个方案有 6140 个未提交，886 个方案有 5657 个未提交。不要不经筛选就提交自动生成的全部 `rerun.csv`。**

回传完整的 `return_packages/日期时间/`，包括报告压缩包、ModelSEED、VMH、BiGG 证据包、校验和及清单。默认打包若生成了空的 CarveMe 归档，也保留在整包中便于校验。另附本次选用的补跑清单和独立验收报告；成功与失败证据都保留。

后续合并按 `job_id` 替换本次计划运行的模型，并保留每个结果的源码提交、发布锁、输入及输出哈希。未受影响而沿用的旧结果必须标记为保留结果，不能伪称都用新提交重新跑过。流程完成、有限供给及能量验收通过，均不等于实验预测准确率已获证明。
