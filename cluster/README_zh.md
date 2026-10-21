# 优化版 MQC 集群运行指南

本指南以论文的6,543个原始模型为例，只运行优化版修正流程，并与论文历史质控结果配对比较，不重跑旧版。优化版已于2026-10-21从 `optimized-20261018` 合入 `main`，科学代码和运行器仍为已验收的 `c9fae34` 版本；AWS网站已于2026-10-20部署同一科学计算版本。新的批次可直接使用 `main`，运行时仍需固定提交号与源码锁。

原批次及定向补跑最终接纳6,516项，其中6,113项沿用 `d7a49ce`，403项使用 `c9fae34`；另有19项失败、6项超时、2项不适用。不要把这一集合标成所有模型都在同一个提交上重跑，也不要把完成数量当成实验预测准确率。

## 需要传到集群的文件

传入 `mqc_inputs_6543.tar` 和 `mqc_inputs_6543.tar.sha256`。输入包约 1.05 GB，包含 CarveMe 5,585、BiGG 73、ModelSEED 79、VMH 806 个原始模型及身份清单。不要改用本机 `tmp/*_model` 的历史质控后模型。

模型不上传 GitHub。GitHub 只提供源码、原有参考库、运行工具和说明。批量程序不读取六株实验数值，也不接收菌株专用补全任务。

以下以旧脚本所在集群的 `/hpcfs/fhome/dengxg/data` 和 `qcpu_23i` 分区为例。如果位置或分区已变，替换相应参数即可。不得覆盖之前的集群结果目录。

```bash
mkdir -p /hpcfs/fhome/dengxg/data/mqc_batch_20261018
cd /hpcfs/fhome/dengxg/data/mqc_batch_20261018
# 将两个输入文件上传到这里后执行
sha256sum -c mqc_inputs_6543.tar.sha256
tar -xf mqc_inputs_6543.tar
git clone --branch main --single-branch https://github.com/dxg-9527/mqc.git mqc-optimized
cd mqc-optimized
git rev-parse HEAD
```

保留 `git rev-parse HEAD` 的输出。整批运行期间不要 `git pull`、改源码、改参考库或更换环境；运行器会校验源码和输入哈希。

## 准备 Python 和 CPLEX

已有集群环境可先 `conda activate dxbase` 检查。但旧环境不一定与本轮验证环境一致；建议单独建环境，以免影响旧 MQC。

```bash
conda create -n mqcopt20261018 python=3.8.16 -y
conda activate mqcopt20261018
python -m pip install -r cluster/requirements-validated.txt
python -m pip install --no-deps -e .
python cluster/check_environment.py --output ../environment_login.json
python cluster/verify_inputs.py --inputs ../mqc_inputs_6543
```

`check_environment.py` 默认要求关键依赖版本与已验证环境一致。每个计算任务还会执行一个超过 1,000 个变量的 CPLEX 小问题，检查求解器是否能处理模型规模。**安装 Python 的 cplex 包不代表已经取得足够规模的求解授权**；若失败，应使用集群已有的合规 CPLEX 安装或联系管理员配置，不跳过测试、不静默换求解器。

还会检查 COBRA 默认求解器是否为本轮验证使用的 CPLEX，避免旧环境里额外安装的求解器被自动优先选中。若被拦截，使用上述独立环境，不静默改变求解器继续计算。

如果只能使用不同版本，先记录差异并小样复核。`submit.py --allow-environment-drift` 是显式允许版本差异的选项，不允许缺少依赖、默认求解器不一致或 CPLEX 测试失败。不同环境应使用不同结果目录，不能混成同一批验证。

不要求集群装 Docker，不上传 AWS 镜像、凭据或 CPLEX 许可证。所有模型计算都在计算节点执行，不在登录节点运行整套 MQC。

## 先提交 24 个小样

小样已按四个数据库和三个模型规模层预选，共 24 个，每库 6 个。先验证集群环境、实际耗时和内存，再扩大并发。

```bash
python cluster/submit.py \
  --inputs ../mqc_inputs_6543 --results ../results_optimized \
  --partition qcpu_23i --pilot --parallel 4 --cpus 10 --mem 8G \
  --time 04:00:00 --model-timeout 7200 --submit
```

如需指定账户，加 `--account 你的账户`。`--python /绝对路径/bin/python` 可以明确选择环境；默认使用当前 Python。提交时继承激活环境，但不会把环境变量或许可证内容写进回传文件。

去掉 `--submit` 只生成计划，不向 Slurm 提交。计划、Job ID 和所有资源参数保存在结果目录的 `plans/`。

```bash
squeue -u "$USER"
# 用提交输出的 Job ID 替换下面的 123456
sacct -j 123456 -P --format=JobID,State,ExitCode,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

注意 `MaxRSS` 受集群记账方式影响，不能未经检查就当作所有子进程的总峰值。退出码 0 也不是验收的唯一条件；运行器还检查六个阶段、导出可加载、末次约束检查和归档完整性。

小样结束后执行：

```bash
python cluster/collect.py --inputs ../mqc_inputs_6543 --results ../results_optimized
```

汇总仍以完整 6,543 个模型为分母。因此小样成功 24 个时，其余 6,519 个显示未完成是正常的。先检查这 24 个的日志、版本、耗时和失败情况，再提交全量；不要把“成功数 24”误读成全量成功。

## 提交全量优化版

```bash
python cluster/submit.py \
  --inputs ../mqc_inputs_6543 --results ../results_optimized \
  --partition qcpu_23i --parallel 32 --cpus 10 --mem 8G \
  --time 24:00:00 --model-timeout 7200 --submit
```

已成功且归档校验通过的小样会自动跳过。失败或不完整的尝试不会覆盖，而是建立新的 attempt 目录。

等上一批任务全部结束再提交下一批，不要同时向同一个结果目录重复提交。并发上限 `%N` 只约束一组数组；同时提交两组会累加资源使用，也可能重复计算尚未完成的模型。

这里始终只使用 **一个 `array.sbatch` 脚本**，不是生成 6,543 份脚本。运行器读取集群 `MaxArraySize`；默认最多 1,000 个数组索引。全量时通常每个数组任务依次处理约 7 个模型，具体分配按反应规模平衡，计划里会列出最大个数。`--parallel 32` 限制整组数组最多同时运行 32 个模型，而不是每个数据库都各跑 32 个。

若集群限制提交任务数，可加 `--array-slots 256`，减少数组任务数量、增加每个任务串行处理的模型数，并相应留足 `--time`。某个模型超时后会终止它的子进程、保存失败证据，再处理该数组任务的下一个模型；Slurm 整个数组任务超时造成的未启动模型会在汇总中显示未完成，可重新分配补跑。

旧脚本只申请 1 个 CPU，而 MQC 内部两个检查阶段会开 10 个子进程。新默认申请 10 个 CPU，不改原内部进程数和算法。批量吞吐还可以通过小样比较 `--cpus 4` 与 `--cpus 10`：4 核时内部进程共享已分配 CPU，必须依照集群政策和实际吞吐选择，不能宣称必然更快。所有任务都设置 BLAS 单线程，避免再嵌套产生大量计算线程。

默认 8 GiB 是起始预算，并非所有模型的内存上限。OOM 模型可单独以 16 GiB 补跑。每模型 2 小时超时是运行保护，不代表超时模型完成质控；必要时补跑时增加 `--model-timeout` 和 `--time`。

## 时间和存储预算

以平均每模型 10–20 分钟、额外 25% 调度余量计算，**仅优化版**理论墙钟时间如下；实际排队、配额、复杂模型长尾另计：

| 同时运行模型数 | 申请 CPU 总数（每任务 10 核） | 申请内存总量（每任务 8 GiB） | 时间估计 |
|---|---:|---:|---:|
| 32 | 320 | 256 GiB | 1.8–3.6 天 |
| 64 | 640 | 512 GiB | 0.9–1.8 天 |
| 128 | 1,280 | 1,024 GiB | 0.4–0.9 天 |

这是持续拿到相应资源的估计，不是实际集群完成承诺。增加 `--parallel` 不能突破用户资源配额。小样结束后用实测值更新预算。

优先使用集群设置的 `SLURM_TMPDIR` 保存单模型临时文件；没有该变量时使用结果目录的 `scratch/`。每个模型完成后，将报告、模型和通量压成一个归档，验证归档后只清理该次新建的临时目录。不删除输入、旧结果或失败尝试。

建议结果存储区至少准备 150 GiB；如果最后还要在同一磁盘生成整套回传包，建议 200 GiB 或以上，具体按小样实际输出量调整。每个模型的 `evidence.tar.gz` 已压缩，整库打包用 `.tar`，不会反复压缩模型。

## 汇总并补跑失败模型

先确认本批 `squeue` 中的任务均已结束，再汇总：

```bash
python cluster/collect.py --inputs ../mqc_inputs_6543 --results ../results_optimized
```

工具会显示本次报告目录 `results_optimized/reports/日期时间/`，包含：

- `summary.csv`：每个模型一行，全部 6,543 个模型都在。
- `counts.json`：各库成功/未完成数。
- `attempts.csv`：全部尝试，包括失败、超时及重复运行。
- `rerun.csv`：需要补跑的模型列表。
- `sacct.txt`：Slurm 状态、退出码、耗时和内存记账；记账不可用时保留错误说明。

有未完成模型时，复制最新 `rerun.csv` 的实际路径：

```bash
python cluster/submit.py \
  --inputs ../mqc_inputs_6543 --results ../results_optimized \
  --ids-file ../results_optimized/reports/实际日期时间/rerun.csv \
  --partition qcpu_23i --parallel 8 --cpus 10 --mem 16G \
  --time 24:00:00 --model-timeout 14400 --submit
```

内存不足和时间不足应按失败原因分别调资源，不必把所有任务都加倍。空 `rerun.csv` 不需要提交。模型自身无法完成时，保留失败日志并回传，不删掉该模型，也不要人工加反应让它通过。

## 最后带回什么

所有任务停止后执行：

```bash
python cluster/collect.py --inputs ../mqc_inputs_6543 --results ../results_optimized --pack
```

将打印的 `return_packages/日期时间/` **整个目录**带回来，至少包括：

1. `MQC_optimized_reports.tar.gz`：汇总、失败清单、Slurm 日志、环境版本、运行计划、模型清单。
2. `MQC_optimized_CarveMe_evidence.tar`。
3. `MQC_optimized_BiGG_evidence.tar`。
4. `MQC_optimized_ModelSEED_evidence.tar`。
5. `MQC_optimized_VMH_evidence.tar`。
6. `SHA256SUMS.txt` 和 `return_manifest.csv`。

每个数据库包包含每个模型的成功与失败尝试。单模型证据中有最终模型、`control_result.json`、`model_info.json`、六个阶段的模型快照、原流程实际产生的中间通量 TSV、最终模型重新加载后的 `final_pfba.tsv`、日志及校验信息。它不声称记录了原流程每一次内部 pFBA 调用；保留的阶段模型可用于必要的进一步分析。

`final_objective_value` 是模型保存环境下的目标通量，不自动等于已归一化的实验生长速率。不要只回传这个字段或 Excel 汇总；缺少最终模型和边界信息就无法做可信的共同条件复检。

回到本机后先在回传目录运行 `sha256sum -c SHA256SUMS.txt`，再告诉我目录位置。我会核对 6,543 个配对、版本与完成状态，并用同一评价标准比较优化版与论文历史模型。改善、退步和失败都保留；这轮不重跑旧版。历史集群环境若无法完整恢复，比较中会将这个来源差异与算法变化分开说明。

不要在任务仍运行时打最终回传包；也不要只传成功模型。可先回传 reports 包查看进度，但完整评价需要四个数据库证据包。

## Slurm 依据

数组索引、并发限制 `%N` 和数组大小上限见 [Slurm Job Array Support](https://slurm.schedmd.com/job_array.html)。运行记账字段和资源统计的适用范围见 [Slurm sacct](https://slurm.schedmd.com/sacct.html)。本机没有目标集群连接，因此分区、账户、配额与排队时间需要在集群端确认。
