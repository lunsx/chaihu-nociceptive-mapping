# UniProt 小项补核与 TrEMBL 裁决报告（第 6 轮）

日期：2026-11-14。脚本 `uniprot_round2_check.py`、`uniprot_trembl_canonical.py`、`uniprot_trembl_adjudicate.py`、`uniprot_name_resolve.py`。产物 `uniprot_round2_check.csv`、`uniprot_trembl_canonical.csv`、`uniprot_trembl_adjudicated.csv`、`uniprot_name_resolved.csv`。

## 一、两个"未解决"项已定案

| 项目 | 原来状态 | 现在 |
|---|---|---|
| **NO66 / URLC2**（赖氨酸去甲基化酶 NO66） | UNRESOLVED | **人类确无 reviewed（Swiss-Prot）条目**；唯一人类记录是 TrEMBL `A0A0A8K8C0`（"双功能赖氨酸去甲基化酶/组氨酸羟化酶"）。结论：保留在"未映射"里，若必须引用则注明"仅 TrEMBL、未审阅" |
| **DHRS4**（脱氢酶/还原酶 SDR 家族 4） | UNRESOLVED | **已定案 Q9BTZ2**（DHRS4_HUMAN，reviewed）。同时确认 `Q9P1J3` 是反义链 `DHRS4-AS1`，**不能当酶用**；另有 3 条 TrEMBL 片段记录（A0ACP2FI66/HAF0/HAF6）应弃 |

## 二、18 条 TrEMBL（未审阅）映射的逐条裁决

先把自动比对的坑说清楚：按"靶点名与 reviewed 蛋白名词语重叠"自动判，会把 **PFKP（血小板型磷酸果糖激酶）和 UMPS（乳清酸核苷酸脱羧酶）** 误判成错配，又会把 **IRP1、MAOA** 这样的"撞名不同蛋白"放过。所以最终表是**逐条人工读 UniProt 条目后裁决**的。

| 结论 | 条数 | 明细 |
|---|---|---|
| **REPLACE**（换成 reviewed 同蛋白） | 2 | PFKP → **Q01813**（ATP-dependent 6-phosphofructokinase, platelet type）；UMPS → **P11172**（UMP synthase，含乳清酸核苷酸脱羧酶活性） |
| **AMBIGUOUS**（同名同工酶，需指定） | 2 | L-岩藻糖异构酶 → FUOM **A2VDF0**（岩藻糖变旋酶，活性不同，历史注释混用）；Uridine phosphorylase → UPP2 **O95045**（与 UPP1 并存） |
| **REJECT_WRONG**（基因符号指向别的蛋白） | 6 | IRP1（顺乌头酸酶/铁反应元件结合蛋白，不是异丙基苹果酸合酶）；CYTB（细胞色素 b ≠ 细胞色素 c4）；Endoglucanase F（微生物酶，人类无）；Gag-Pol（病毒多聚蛋白）；MAOA（FAD 单胺氧化酶 ≠ 膜结合原胺氧化酶 AOC3/VAP-1）；WARS1（靶点写的是 synthetase **2**，应为线粒体 **WARS2 Q9UGM6**） |
| **KEEP_TR**（确无可换） | 8 | 无基因符号的片段记录为主；其中 URLC2 人类无 reviewed |

## 三、按蛋白名补定的 4 条（无基因符号那批）

| 靶点名 | 定到的 reviewed 条目 | 备注 |
|---|---|---|
| Adenylsuccinate synthetase | **ADSS2 P30520** / ADSS1 Q8N142 | 注意：`ASS1 P00966` 是**精氨琥珀酸合成酶**，不是腺苷琥珀酸合成酶——这是很容易顺手写错的一对，已避开 |
| Methylmalonyl-CoA mutase | **MMUT P22033** | 线粒体 |
| Branched-chain aminotransferase | **BCAT1 P54687** / **BCAT2 O15382** | 胞质 / 线粒体两个同工酶 |
| Alkaline phosphatase | **ALPL P05186 / ALPP P05187 / ALPG P10696 / ALPI P09923** | 四个同工酶，泛名无法唯一化，保留歧义 |
| Serine hydroxymethyltransferase（1/2） | **SHMT1 P34896 / SHMT2 P34897** | 已在 `uniprot_trembl_adjudicated.csv` 内标注 |

## 四、对论文数字的影响

- 691 行已映射 / 255 未映射、662 唯一人类基因这套口径**不受影响**：这 18 条本来就在"已回捞"的 74 条里，换成 reviewed 的对照号只是提高可信度，不改变基因数量。
- 但**如果审稿人查补充表里的 UniProt 号**，会有 2 条（PFKP、UMPS）显示"unreviewed"，现在可以直接给 reviewed 号；另有 6 条属于"撞名错配"，如果这些靶点进了正文的靶点清单，应该标注为非人源/错配，而不是当作人类靶点。

## 五、还没做

- 没有回头把这 18 条的**更正号回写进 `target_list_enriched.json` 与 `靶点总清单_证据分级.csv`**（这两份是论文与补充表的源头）。回写属于"改数据"，应在末轮统一改版时一次做完并同步 CHECKSUMS。
- 没有对 946 条里**其余 691 条 reviewed 记录**做同样的"名称-基因一致性"抽检（本轮只查了有疑问的 TrEMBL 段）。
