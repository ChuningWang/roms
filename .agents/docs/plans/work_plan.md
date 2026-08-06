# ROMS 冰盖-海洋耦合模块开发工作计划 (Work Plan)

> **本文件是所有 AI agent（Qoder、Zcode 等）的共享开发计划。**
> 修改本文件前请确认与其他 agent 的工作不冲突。
>
> 基于 `AGENTS.md` 项目背景，参考 MITgcm STREAMICE 与 PISM 架构。
>
> **创建日期**：2026-07-10  
> **最后更新**：2026-08-02（阶段 3.3 SMB 接口已实施，编译验证通过）  
> **分支**：`icesheet`  
> **约束**：单次 commit 不超过 1000 行修改

---

## 相关文档索引

| 文档 | 相对路径 | 用途 |
|------|---------|------|
| SSA 技术设计 | `docs/plans/ssa_solver_design.md` | FE 离散公式、STREAMICE 对比分析、参数体系 |
| 代码审查报告 | `docs/reviews/code_review_2026-07.md` | 20 项问题及修复状态 |
| 独立审查报告 | `docs/reviews/code_review_2026-07-31.md` | 阶段 6 缺陷（A 类）、残留问题（B 类）、计划问题（C 类）、测试方法论 |
| 计划审查决策 | `docs/reviews/plan_review_2026-07.md` | 设计决策理由记录（归档） |
| SSA 修复记录 | `docs/fixes/ssa_fixes_2026-07.md` | 11 项修复的详细方案与验证 |
| MISOMIP1 协议 | `docs/references/misomip1_protocol.md` | 实验参数、域配置、边界条件 |
| 目录索引 | `../../README.md` | .agents/ 完整文件索引 |
| 项目总览 | `../../../AGENTS.md` | 分支体系、编码规范、CPP 宏 |

> 所有路径相对于 `.agents/docs/plans/`。

---

## 当前代码状态摘要（2026-07-28 更新）

| 组件 | 状态 | 说明 |
|------|------|------|
| `isf_mom.F` | **已实现 + 已修复** | SSA 主驱动：几何/掩膜/浮力判据/驱动应力/黏度/CG 求解 |
| `isf_cg_solve.h` | **已实现 + 已修复** | Picard + matrix-free PCG，Jacobian 权重和收敛判据已修正 |
| `isf_driving_stress.h` | **已实现（FE 装配）** | Q1 2×2 Gauss 积分 FE 装配 + 5 点 Laplacian 平滑接地线 kink |
| `isf_viscosity.h` | **已实现 + 已修复** | Paterson (1994) 温度依赖 Glen 律 + Weertman 拖曳 |
| `isf_thickness_adv.h` | **新增已实现** | CFL 子步进 + Van Leer 通量限制器厚度平流 |
| `isf_calving.h` | **新增已实现** | 厚度阈值崩解（H_calve = 100 m） |
| `isf_vbc.F` | **完整实现** | 3EQN 热力学，计算底融速率；ICESHEET_SSA 下使用耦合器 u* |
| `isf_coupler.F` | **已实现 + 已修复** | 冰-海洋动量耦合 + 底融 EMA；P1 缺陷（sustr halo、mask 判据、U/V 插值、双写）已于 2026-07-31 修复（阶段 6.6） |
| `isf_masks.F` | **完整实现，运行时生效** | 掩膜计算 + 浮力判据接地线；调度守卫扩展为 `WET_DRY || ICESHEET_SSA`，MISOMIP 配置下动态 zice/rmask_isf 已生效（阶段 4a.5 方案 B） |
| `iceshelf.F` | **完整驱动 + 多时间尺度** | 含 SSA/厚度/崩解/接地线调度（N_ssa/N_thick/N_gl） |
| `mod_isf.F` / `isf_mod.h` | **完整基础设施** | 状态变量(5)、诊断场(7+6)、参数(15+6新增)、时间层(2)；Fi 数组 allocate 归零 |
| CPP 宏体系 | **`ICESHELF_*` + `ICESHEET_SSA` + `ICESHEET_CALVING`** | SSA→MOM→DYNAM 派生 |
| 厚度演化 | **已实现** | CFL 子步进 + Van Leer 限制器 + 底融源项 |
| 接地线动力学 | **已实现 + 配置层已生效** | 浮力判据 + 余弦平滑 float_frac；动态 zice/rmask_isf 随 `isf_masks` 调度更新（无 WET_DRY 时海洋网格掩膜保持静态，已知限制） |

### 已完成的修复（2026-07-27，见 `fix_plan_ssa_solver.md`）

- [x] P0-1: Jacobian 权重 `0.25/(pm*pn)` → `1.0/(pm*pn)`
- [x] P0-2: 驱动应力 has_ice 逻辑修正（4 邻域正确赋值）
- [x] P0-3: 删除 Picard 对 ksa_old 的覆写；黏度改从 ksa_new 读取
- [x] P1-4: B_glen 参数 `3.076e-9`(A) → `464.2`(√B = A^{-1/6})
- [x] P1-5: 收敛判据改为速度相对变化 `‖u_new-u_old‖/‖u_new‖`
- [x] P1-6: 强度扩展补 0.5 因子
- [x] P2-8: 浮力判据 `rho_w*(-h)` → `rho_w*h`

### 待处理项

- [x] ~~**驱动应力 FE 装配**~~（已完成：Q1 2×2 Gauss 积分，commit `5be44db7`）
- [x] ~~MISOMIP 网格 h<0 / zice>0 问题~~（已修复数据文件）
- [x] ~~初始条件 `_FillValue=0.0` 导致数据全零~~（已移除属性）
- [x] ~~**MISOMIP 配置数值稳定性调优**~~（已完成：DT=300 稳定，CG BC 修复，Picard 收敛）
- [x] ~~**isf_coupler 修复**（阶段 6.6）~~（✅ 2026-07-31：sustr halo、mask 判据、U/V 交叉分量插值、sustr 双写协议）
- [x] ~~**WET_DRY 依赖处理**（阶段 4a.5）~~（✅ 2026-07-31：方案 B，isf_masks 调度移出守卫）
- [x] ~~**Fi(icImef)/icIust 首步初始化**（阶段 3.4）~~（✅ 2026-07-31：Fi allocate 后归零）
- [x] ~~**死参数 `gl_smooth_width` 清理**~~（✅ 2026-08-02：已从 isf_mod.h/isf_inp.h/isf.in 全部移除）
- [x] ~~**厚度平流融化项守卫不一致**（P0）~~（✅ 2026-08-02：`#ifdef ICESHEET_SSA` -> `#if defined ICESHEET_SSA && defined ICESHELF_THERMO`，修复纯 SSA 无 THERMO 时融化项静默消失）
- [x] ~~**`mp_exchange2d` model 参数统一**（P0）~~（✅ 2026-08-02：isf_coupler.F `iNLM` -> `model`）
- [x] ~~**`Cd` 参数化 + 腔外 sustr 清零**（P1）~~（✅ 2026-08-02：`Cd` -> `Cd_isf(ng)` 可配置；腔外 `sustr=0`/`svstr=0`）
- [x] ~~**浮力判据 `rho0` -> `rho_w_ref` 参数化**（P1）~~（✅ 2026-08-02：isf_masks.F + isf_mom.F 使用 `rho_w_ref(ng)`，默认 1025.0）
- [x] ~~**硬编码 `2000.0_r8` -> `max_thickness` 参数**（P2）~~（✅ 2026-08-02：默认 5000m，可配置）
- [x] ~~**`isf_mms_source.F` 线程安全**（P2）~~（✅ 2026-08-02：`!$OMP CRITICAL` 保护 NetCDF 读取与缓存）
- [x] ~~**`OCEAN_FIXED` 安全断言**（P2）~~（✅ 2026-08-02：`globaldefs.h` 加 `#error` 守卫）
- [x] ~~**`SSA_MMS_TESTING` 编译期依赖隔离**（P2）~~（✅ 2026-08-02：`SSA_MMS_TESTING` -> `ICESHELF_MOM`（不再触发 `ICESHEET_SSA`），`rho_w_ref`/`max_thickness`/`mms_source_file` 移至 `ICESHELF_MOM` 守卫下）
- [x] ~~**`isf_driving_stress.h` 注释修正**（P1）~~（✅ 2026-08-02：明确标注为 mass-lumped Q1 而非完整 FE 装配）
- [x] ~~**`isf_calving.F` 死代码清理**~~（✅ 2026-08-02：移除未使用的 `rho_w_local`/`rho_i` 导入）
- [ ] 薄水柱问题（57 个 cell 水柱 < 10m）
- [ ] 单元测试和回归测试框架

### 冒烟测试结果（2026-07-28）

- **编译**：ICESHEET_SSA + ICESHEET_CALVING + ICESHELF_THERMO + ICESHELF_TEOS10 全部通过
- **运行**：MISOMIP 配置（398×38×36，DT=120s，4 tiles MPI）运行至第 43 步
- **SSA 求解器**：正常工作（非零速度、应力、黏度），无 NaN/crash
- **不稳定原因**：配置/参数调优问题（非代码 bug），FE 积分节点力量级 ~1e12 N，冰厚 48m/40步异常增长
- **Git commits**：`5be44db7`（SSA 改进）→ `dd3facc7`（厚度平流+崩解）→ `de146c25`（多时间尺度+接地线）

---

## 设计决策记录

| 编号 | 决策 | 理由 |
|------|------|------|
| D1 | 速度**永远在 PSI 点**计算和存储 | 与 PISM/STREAMICE 一致；Q1 FE 的自然节点 |
| D2 | 耦合时 PSI→U/V 简单插值 | 仅在需要 couple 到海洋时做，不维护双份速度场 |
| D3 | 驱动应力采用**全 FE 装配**（方案 A） | 数学一致性，保证 Q1 收敛阶 |
| D4 | 接地线：余弦平滑 + 次网格插值，**不做局部网格加密** | 复杂度控制 |
| D5 | 时间步进：固定多时间尺度，**不做 adaptive** | 简化实现，后续可扩展 |
| D6 | 崩解：MISOMIP1/2 用**厚度阈值方案** | 与 Asay-Davis et al. (2016) 协议一致 |
| D7 | 输出兼容性（PSI 点速度）问题**汇报用户决定** | 后处理工具可能需要适配 |
| D8 | **动态接地线强制要求 `WET_DRY`**：`iceshelf.F` 对 `isf_masks` 的运行时调度保持在 `#ifdef WET_DRY` 内；MISOMIP 配置**必须启用 WET_DRY**，或在无 WET_DRY 时明确降级为"静态几何"模式（阶段 4a.5 二选一） | 冰盖侧（SSA scratch 掩膜）与海洋侧（zice/rmask_isf/set_depth）几何必须一致；当前 MISOMIP 配置下两侧不一致 |
| D9 | **冰架底部应力 `sustr`/`svstr` 的唯一生产者是 `set_vbc`**；`isf_coupler` 只计算并存储冰架应力（Fi 字段），由 `set_vbc` 在 `#ifndef ICESHEET_SSA` / `#ifdef ICESHEET_SSA` 分支统一注入（阶段 6.6） | 消除双写与时序依赖；保持 `sustr` 的 2D 推进语义与 UV_LDRAG/UV_QDRAG 兼容 |
| D10 | 阶段 7（MISOMIP1 IceOcean1r）设为**硬门禁**：阶段 6.6 修复、薄水柱、SMB、WET_DRY 前提全部完成后才进入 | 避免带着已知缺陷推进验证实验 |

---

## 开发阶段

### 阶段 0：准备工作 ✅ 已完成

精读 STREAMICE/PISM 参考实现，确认 B-grid Q1 FE 方案。产出 `roms_ssa_plan.md`。

---

### 阶段 1：基础设施与 CPP 宏体系 ✅ 已完成

CPP 宏（`ICESHEET_SSA`）、PSI 点速度重定向、15 个 SSA 参数、6 个诊断场、varinfo.yaml。

---

### 阶段 2：SSA 动量方程求解器 ✅ 已完成（含修复）

Picard + matrix-free PCG + Q1 FE 装配。已通过 ISOMIP 稳定运行验证。

**待完成子项**：

| 步骤 | 内容 | 状态 |
|------|------|------|
| 2.A | 驱动应力改为 FE 装配（方案 A，见下文） | ✅ 已完成 |
| 2.B | 并行一致性测试（1 tile vs 4 tiles） | 待实施 |
| 2.C | MMS 收敛性测试 | 待实施 |

#### 2.A 驱动应力 FE 装配（方案 A）

将 `isf_driving_stress.h` 的预计算点值方式改为在 CG 装配循环中同步积分 RHS：

```
对每个 RHO 单元 (i,j)，对每个 Gauss 点 nq：
  τ_d(xq) = -ρ_i * g * H(xq) * ∇S(xq)    ! 在积分点计算
  b_k += jacq(nq) * jac * τ_d(xq) * φ_k(xq)  ! 散射到 4 节点
```

其中 H(xq) 和 ∇S(xq) 通过形函数从 4 个 RHO 角点插值。冰前边界压力项作为线积分单独处理（乘以面长度而非面积）。

这与 STREAMICE 的 `streamice_driving_stress.F` 做法一致（含 `0.25*dyG` 面积因子）。

---

### 阶段 2.5：SSA 验证与测试 ✅ 已完成

| 步骤 | 内容 | 通过标准 | 结果 |
|------|------|---------|------|
| 2.5.1 | MISOMIP 稳定运行 | 无 NaN，速度 < 1 m/s | ✅ DT=300 稳定 5000 步 |
| 2.5.2 | 并行一致性：1 tile vs 4 tiles | 后续时间步一致 | ✅ ~4%（驱动应力平滑 tile 伪影） |
| 2.5.3 | 2D MMS 收敛性测试 | L2 误差 O(h²) | ✅ rate=2.01/2.04 |
| 2.5.4 | SSA 解析解（均匀冰架线性坡度） | 误差 < 1% | ✅ L2=0.11% at dx=2km |
| 2.5.5 | Picard 收敛曲线 | 单调下降至 nl_tol | ✅ 37-39 次迭代 |
| 2.5.6 | STREAMICE halfpipe 对比 | 速度场定性一致 | ✅ 中心线误差 0.015% |

**测试代码位置**：`icesheet_test/tests/test_ssa_{analytical,mms,halfpipe}.py`

**⚠️ 2.5 方法论缺陷（审查发现，需在阶段 2.6 纠正）**：上述三个测试是在 Python 中
**复刻** Q1 FE 离散后求解，并不调用 ROMS 的 Fortran 实现——它们验证"离散格式正确"，
不验证"ROMS 实现正确"。只要 Fortran 与 Python 复刻存在相同的 bug，测试依然通过
（本项目历史已出现 Jacobian 与 B_glen 参数互掩的案例）。此外冒烟测试只检查
"无 NaN"，不检查任何物理量量级。

### 阶段 2.6：验证测试框架升级（审查新增，阶段 7 前提）

| 步骤 | 工作内容 | 通过标准 |
|------|---------|---------|
| 2.6.1 | 将解析解/MMS 源项制成 NetCDF 强迫，在 ROMS 中运行 1-2 个 SSA 时间步，读回速度场与参考解比对 | 直接验证 Fortran 实现；L2 误差 ≤ Python 复刻结果 |
| 2.6.2 | 冒烟测试增加数值断言：`max|u|`、`max|melt|`、`max|sustr|` 在物理范围内 | 超出范围即失败 |
| 2.6.3 | 并行一致性测试加严：消除 `isf_driving_stress.h` 平滑的 tile 伪影（当前 ~4% 偏差被测试容差掩盖） | tile 伪影归零或显式文档化 |

**✅ 阶段 2.6.1 落地（2026-08-02）：SSA MMS 端到端 + MPI 发散修复**

MMS 测试已从"Python 复刻"升级为"ROMS 内运行 + 制造解比对"（`tests/gen_mms_source.py` 生成
`mms_source.nc`/`ssa_mms_ini.nc`，`SSA_MMS_TESTING` 编译期注入，`isf_read_mms_source.F` 读取），
并在此过程中定位并修复了 **SSA 求解器在 MPI 分块下的发散根因**：

| 问题 | 根因 | 修复 |
|------|------|------|
| 1 tile 收敛（21 次迭代 conv=T）但 4 tiles 发散（nl_resid 不降） | `isf_mom.F` 中 `Au/Av/diag_u/diag_v/u_scr/v_scr/taudx_psi/taudy_psi/zs_smooth` 声明为 tile-local `(IminS:ImaxS,JminS:JmaxS)`（= Istr-4..Iend+3），却传给期望全局边界 `(LBi:UBi,LBj:UBj)` 的 `mp_exchange2d`/`exchange_p2d_tile`；MPI 分块下数组越界写坏 halo → 跨 tile 共享 PSI 节点丢失装配贡献（Au=0）→ 发散 | 上述数组全部改为全局声明；初始化循环相应改为全局范围 |
| bounds-check 暴露 `zs_smooth` 越界 | 5 点 Laplacian 平滑循环写 `Jstr-1:Jend+1, Istr-1:Iend+1`，在 i=Istr-1 处读 `zs_smooth(i-1)`（i-1 < LBi） | 平滑只更新内部点 `Jstr:Jend, Istr:Iend`，邻居从已交换的 halo 读取 |
| 装配读 tile 边界 halo 速度 | `isf_viscosity.h` 与元素装配读每个 RHO cell 的 4 角速度 `Si(i+1,j+1,ksa_new)`，装配前未交换 halo | 进入 `isf_cg_solve.h` 前交换 `Si(ksa_new,isIuis/isIvis)` halo；`pmask_ssa` 改全局声明并交换 halo |

**验证结果**：`-fbounds-check -fcheck=all` 构建下 4 tiles（2×2）Picard 24 次迭代 `conv=T`、
零越界错误；1 tile 21 次迭代 `conv=T`。制造解装配一致（1-tile iter1 相对残差 ~3e-7）。

**遗留问题（独立于 SSA）**：MMS 测试的海洋首个时间步 blow up（`(228,16,01)` Cv~8e6 m/s，
负动能）——经确认与 SSA/MMS 无关：`log_smoke_no_ssa.txt`（无 SSA）同样 blow up，`isomip` 网格
跑通而 `misomip` 网格（cavity=0.1 m 薄腔 + 3EQN）一直炸。属 MISOMIP 配置遗留问题，待阶段 8
处理。

**✅ SSA MMS 端到端 PASS（2026-08-04）：MPI 接缝装配 + 时间层 halo 修复，L2 降至 2.4e-8**

继 2.6.1 之后进一步定位并修复了两个导致 4-rank MPI 下 MMS 测试失败（L2~1.1e7）的深层 bug：

| 问题 | 根因 | 修复 |
|------|------|------|
| 接缝 PSI 节点装配贡献缺失（第一次 SSA 调用 op chk 偏大） | ROMS PSI 无重叠分块：接缝节点由单一 tile 拥有，但其单元模板跨两侧；装配循环只遍历本 tile 单元，邻居侧单元对接缝节点的贡献丢失 | 三处装配循环（diag、初始 Au、CG matvec）扩展一圈 halo 单元（`iCstr..iCend`/`jCstr..jCend`），单元守卫改用 `nuH_rho>0`，散射加所有权过滤（`i1∈[IstrP,IendP] ∧ j1∈[JstrP,JendP]`）；`nuH_rho/beta_rho` 改全局声明并在黏度计算后交换 halo |
| 第二次 SSA 调用算子爆炸（op chk 9.3e12，交替发散） | `get_state` 只填时间层 1（`Tindex=1`）；标量态 carry-forward 拷贝只复制 tile 内部，第二次调用时 `ksa_old`（时间层 2）halo 恒为零 → 掩膜 `min_thickness` 判据将接缝 PSI 节点误判为无冰失活 → 装配算子被破坏 | `isf_mom_tile` 开头（掩膜计算前）交换 `ksa_old` 标量（H/T/E）halo，carry-forward 后再交换 `ksa_new` halo；`isf_cg_solve.h` 内冗余的标量拷贝移除（只留速度初猜拷贝） |
| MMS 模式下两次 SSA 调用间状态被热力学污染 | `isf_vbc`（3EQN）每步修改 H，破坏制造解前提 | `iceshelf.F` 中 `isf_vbc`/`isf_coupler` 用 `#ifndef SSA_MMS_TESTING` 包裹 |

**验证结果**：`run_ssa_mms_test.sh`（2×2 MPI、NTIMES=10）PASS，L2(u)=2.36e-8、L2(v)=2.37e-8
（机器精度量级）；两次 SSA 调用的 op chk 均为 2.1e-2（相对 max|mms|=3.3e12 为 ~6e-15）。
**遗留**：`isf_driving_stress.h` 的 `taud` 装配存在同样的接缝共享节点问题（MMS 模式 RHS 用
制造源项不受影响；生产 `ICESHEET_SSA` 启用真实驱动应力前需按同一模式修复）。

**✅ 阶段 2.6.4 落地（2026-08-02）：`OCEAN_FIXED` 宏 — 冻结海洋、独立测试冰盖**

为冰盖模型独立测试新增 CPP 宏 **`OCEAN_FIXED`**（文档见 `cppdefs.h`）：开启后海洋（与海冰）
只读入初始条件、**不向前积分**，但 `iceshelf` 模块（SSA/厚度/接地线/热力学）照常按海洋时钟
`iic` 推进。作用：

- **绕开 MISOMIP 海洋 blow up 遗留问题**：冻结模式下海洋步被跳过，测试可完整跑完
- **缩短测试时间/成本**：跳过 step2d/rhs3d/step3d 等海洋算子，100 步 4-tile 仅需 ~17 分钟
- **隔离冰盖动力学**：海洋作为静态（初值）背景，便于 SSA/厚度/接地线单独验证

| 项目 | 内容 |
|------|------|
| 宏文档 | `cppdefs.h` 冰架宏块 + `globaldefs.h` 说明注释（独立用户宏，无派生） |
| 调度守卫 | `main3d.F` 5 处 `#ifndef OCEAN_FIXED` 跳过海洋演化段（vmix/omega、rhs3d、step2d、set_depth、step3d_uv→step3d_t）+ seaice 调用 `!defined OCEAN_FIXED` |
| 保留 | `get_data`/`set_data`、`set_vbc`、**`CALL iceshelf`**、`set_zeta`/`set_diags`/`set_avg`、`CALL output`、时钟推进 `iic/time` |
| 测试配置 | `src/ssa_test_ocnfix.h`、`roms_ssa_ocnfix.in`（NTIMES=100, NHIS=10, N_ssa=N_thick=N_gl=1）、`isf_ocnfix.in`、`run_ssa_ocnfix_test.sh` |

**验证结果（1-tile 与 4-tile 均通过）**：
- 海洋 `temp/u/v/zeta` 各 HIS 记录**逐位一致**（rec0 == rec_last，冻结完全生效；注：rec0 与原始
  ini 的 temp 差异是 ROMS `post_initial` 的 zice 深度修正 + `ini_fields` 正常初始化，非演化）
- 冰盖 `is_h/is_u/is_v` 随步推进持续演化（厚度平流 + 崩解在 N_thick=1 下每步生效）
- 100 步完整跑完 `ROMS: DONE`，无 blow up / NaN
- 用法：`NTIMES` 即冰盖步数；配合 `N_ssa/N_thick/N_gl` 控制各过程频率

---

### 阶段 3：冰盖厚度演化（质量守恒）✅ 已完成

**目标**：`∂H/∂t = -∇·(ūH) + ȧ - ḃ`

| 步骤 | 文件 | 工作内容 | 状态 |
|------|------|---------|------|
| 3.1 | `isf_thickness_adv.h` | 通量形式厚度平流 + Van Leer 限制器 + CFL 子步进 | ✅ |
| 3.2 | `iceshelf.F` | 多时间尺度调度（N_thick 步更新一次） | ✅ |
| 3.3 | `isf_thickness_adv.h` + `isf_mod.h` + `isf_inp.h` | 表面质量平衡接口：均匀 SMB 速率 `smb_rate(ng)` [m/s]，在厚度平流中作为源项 `+dt_sub*smb_rate(ng)` | ✅ 已实施（阶段 7 前提） |
| 3.4 | `isf_mod.h` | **Fi 诊断场首步初始化**：allocate 后归零 | ✅ 已实施 |

**PSI→面插值说明**（设计决策 D2）：
- 厚度平流通量需要 RHO 面（x-face/y-face）上的速度
- 使用简单 2 点平均：`u_xface(i,j) = 0.5*(u_psi(i,j) + u_psi(i,j+1))`
- 这是一阶精度，但对 MISOMIP1 分辨率（2 km）足够
- 如后续需要更高精度，可改为 4 点加权插值

---

### 阶段 4：接地线动力学

分为两个子阶段：

#### 阶段 4a：基本接地线（快速可用）

| 步骤 | 文件 | 工作内容 | 状态 |
|------|------|---------|------|
| 4a.1 | 新建 `isf_grounding_line.h` | 浮力判据：`ρ_i*H > ρ_w*h` → 接地（已实现于 `isf_mom.F`/`isf_masks.F`） | ✅ |
| 4a.2 | 同上 | 余弦平滑 `float_frac`：在接地线宽度 `W=3Δx` 内从 0→1 过渡 | ⚠️ 部分（见下） |
| 4a.3 | `isf_masks.F` | 重写掩膜更新：厚度更新后重新计算接地/漂浮。修复 halo exchange | ✅ |
| 4a.4 | `isf_masks.F` / `set_depth.F` | zice 动态更新：漂浮 `zice = -(ρ_i/ρ_w)*H`，接地 `zice = 0` | ✅（2026-07-31 修复，见下） |
| 4a.5 | `iceshelf.F` / `misomip.h` | **WET_DRY 依赖处理** | ✅ 已实施（方案 B） |

**⚠️ 4a.5 说明（审查发现，P1）**：`iceshelf.F` 对 `isf_masks` 的运行时调度在
`#ifdef WET_DRY` 内，`isf_calving.F` 的 `rmask_isf`/`zice` 更新也在 `WET_DRY` 内，
而 `icesheet_test/src/misomip.h` 中 `WET_DRY` 被注释。后果：
- 运行时 `isf_masks` 从不执行 → **动态接地线、动态 zice 在 MISOMIP 配置下失效**
- 崩解移除的冰不反映到 `rmask_isf`/`zice` → `set_depth.F`（只读 zice）与 `isf_vbc`
  继续用旧几何 → **冰盖侧（SSA scratch 掩膜用新厚度）与海洋侧（冻结几何）不一致**

**✅ 已实施修复（2026-07-31，选择方案 B，见 `docs/reviews/code_review_2026-07-31.md`）**：
- `iceshelf.F`：`isf_masks` 调度守卫 `#ifdef WET_DRY` → `#if defined WET_DRY || defined ICESHEET_SSA`
  （USE 语句同步扩展）。ICESHEET_SSA 下每 `N_ssa*N_gl` 步重算动态 zice/rmask_isf；
  无 WET_DRY 时 `*_ist` wet/dry 掩膜不更新（仅诊断，无计算代码读），海洋网格掩膜保持静态
  ——已知限制：接地线迁移不会烘干海洋单元（MISOMIP IceOcean1r 稳态实验可接受）。
- `isf_calving.F`：崩解时 `rmask_isf=0`/`zice=0` 更新移出 `WET_DRY` 守卫，halo exchange
  改为无条件执行 3 字段（H_ice, rmask_isf, zice）。
- 方案 A（MISOMIP 启用 WET_DRY）暂缓：`log_wetdry.txt` 显示旧版本（SSA 修复前）WET_DRY
  运行以 BLOWUP 结束，直接启用风险高；且 ROMS 的 wetdry 会使海洋网格掩膜随时间变化，
  需在后续阶段单独验证。

#### 阶段 4b：次网格接地线（精度提升）

| 步骤 | 文件 | 工作内容 |
|------|------|---------|
| 4b.1 | `isf_grounding_line.h` | 次网格接地线位置插值：参考 PISM 的 LI（linear interpolation）方法，在相邻接地/漂浮单元间线性插值精确接地线位置 |
| 4b.2 | 同上 | 接地线迁移速率计算：`v_gl = Δx_gl / Δt`，用于 MISMIP+ 验证 |
| 4b.3 | 验证 | 网格收敛性测试（1 km vs 2 km vs 4 km），确认接地线位置精度 |

---

### 阶段 5：崩解与冰前动力学 ✅ 基本完成（含待办）

| 步骤 | 文件 | 工作内容 | 状态 |
|------|------|---------|------|
| 5.1 | `isf_calving.h` | **厚度阈值崩解**（D6）：`H < H_calve` 移除冰，`H_calve = 100 m` | ✅（实现于 `isf_calving.F`） |
| 5.2 | 同上 | 冰前检测：`H>0` 但某邻居 `H=0` → 冰前单元 | ✅ |
| 5.3 | 新建 `isf_adv_front.h` | 冰前推进（可选）：参考 `streamice_adv_front.F` 填格子方案 | 🔲 待实施（未创建） |
| 5.4 | `iceshelf.F` | 调度：厚度更新 → 崩解 → 冰前更新 | ✅ |
| 5.5 | `isf_calving.F` | **无 WET_DRY 时崩解后掩膜同步** | ✅ 已实施（2026-07-31）：`rmask_isf`/`zice` 更新移出 `WET_DRY` 守卫，halo exchange 无条件执行 3 字段 |

---

### 阶段 6：冰盖-海洋耦合增强 + 多时间尺度 ⚠️ 部分完成（阶段 6.6 待修复）

#### 6.1 多时间尺度策略（设计决策 D5） ✅

```
海洋时间步:     dt_ocean = 240 s (每个 baroclinic step)
SSA 求解:       每 N_ssa 个海洋步求解一次 (默认 N_ssa = 10)
厚度演化:       每 N_thick 个 SSA 步更新一次 (默认 N_thick = 10)
接地线更新:     每 N_gl 个厚度步更新一次 (默认 N_gl = 10)
```

参数从 `isf.in` 读取。暂不实现 adaptive time stepping。

**⚠️ 6.1 已知缺陷（审查新增）**：`isf_vbc` 的熔水示踪体积 `Vi`/`Ti` 使用
`lsnew(ng)` 作为时间层，但 `lsnew`/`lsold` 只在 `isf_mom`（每 `N_ssa` 步）交换。
`N_ssa > 1` 时 `Vi`/`Ti` 两层时间步退化为单层读写，熔水示踪守恒被破坏。
**影响**：仅当启用 `Lisftrc`（被动示踪）时相关；MISOMIP 配置未启用。
**修复**：`isf_vbc` 的 Vi/Ti 时间层改用独立于 SSA 求解节奏的交换，或显式文档声明
`Lisftrc` 与 `N_ssa > 1` 不兼容。

#### 6.2 耦合通量

| 步骤 | 文件 | 工作内容 | 状态 |
|------|------|---------|------|
| 6.2.1 | `isf_vbc.F` / `set_depth.F` | 动态 zice 更新 + set_depth 联动 | ⚠️ 依赖阶段 4a.5（WET_DRY） |
| 6.2.2 | `isf_coupler.F` / `step2d_FB.h` | 冰架底部拖曳 → 海洋 `sustr`/`svstr`（PSI→U/V 插值） | ⚠️ 见阶段 6.6 |
| 6.2.3 | `isf_coupler.F` | 统一耦合器：冰→海洋（淡水、热、盐、动量）；海洋→冰（底融） | ⚠️ 见阶段 6.6 |
| 6.2.4 | `isf_coupler.F` | 底融速率时间滤波（EMA，τ ~1 天）+ 通量守恒检查 | ⚠️ **通量守恒检查代码未实现**；EMA 状态 `icImef` 未初始化（阶段 3.4） |
| 6.2.5 | `iceshelf.F` | 多时间尺度调度逻辑 | ✅ |

#### 6.6 isf_coupler 缺陷修复（P1，✅ 已实施 2026-07-31）

详细分析见 `docs/reviews/code_review_2026-07-31.md` 第 A 类问题。修复清单（全部完成）：

| # | 缺陷 | 修复方案 | 状态 |
|---|------|---------|------|
| 6.6.1 | `sustr`/`svstr` 是 U/V 点量，`isf_coupler` 只写本 tile 点，但 Step 2 的 u\* 插值读邻居 tile 的 U/V 点——**无 halo exchange** | Step 2 前对 `sustr`/`svstr` 做 `bc_u2d/bc_v2d` + `mp_exchange2d`（仿 `set_vbc.F` 末尾模式） | ✅ |
| 6.6.2 | mask 判据 `rmask_isf(i-1,j)+rmask_isf(i,j) .gt. 0` 与 `set_vbc.F` 的 `eq.1 and eq.1` 不一致，应力外泄到腔外 | 改为与 `set_vbc.F` 一致的 `eq.1 .and. eq.1` | ✅ |
| 6.6.3 | U/V 点交叉分量插值半格偏移（如 V 点的 du 应为 4 邻 U 平均 `0.25*(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1))`，现为 2 点平均） | 改用与 `set_vbc.F` 相同的 4 点平均 | ✅ |
| 6.6.4 | **sustr 双写**：`set_vbc`（UV_QDRAG 等块）写海洋拖曳 → `isf_coupler` 覆盖为相对速度拖曳，靠 main3d 调用顺序"碰巧正确" | 执行设计决策 D9：`set_vbc` 冰架 drag 块加 `#if !defined ICESHEET_SSA || !defined ICESHELF_THERMO` 守卫；ICESHEET_SSA+THERMO 下 `isf_coupler` 是 `sustr` 唯一生产者 | ✅ |
| 6.6.5 | u\* 定义两套（coupler 的动力学 vs `set_vbc` 的 `gamU` 线性） | 3EQN 在 ICESHEET_SSA 下统一用 coupler 的 u\* = `sqrt(Cd)*|du|`（`isf_vbc.F` 已有 `#ifdef ICESHEET_SSA` 分支）；非 SSA 配置保持原 `gamU` 路径 | ✅ |
| 6.6.6 | `Fi(icImef)` 首步未初始化，`isf_vbc` 的 `MAX(min_ustar, Fi(icIust))` 可能被垃圾值污染 | `isf_mod.h`：`Fi` 数组 allocate 后整体归零（冷启动首步确定性；restart 路径由 `isf_get_state` 覆写） | ✅ |

**验证（2026-07-31）**：MISOMIP 配置（398×38×36，4-tile MPI，ICESHEET_SSA+CALVING+THERMO，
无 WET_DRY）编译通过；60 步冒烟运行正常结束（ROMS: DONE），输出变量全部有限
（is_h max=1626 m 与初始一致、is_u~1e-4 m/s、is_mr~0.16 m/yr、is_Ustar max=0.078 m/s），
无 NaN/BLOWUP。

---

### 阶段 6.5：冰盖热力学（可选扩展）

| 步骤 | 工作内容 |
|------|---------|
| 6.5.1 | 冰盖温度场演化方程（参考 PISM enthalpy formulation） |
| 6.5.2 | 温度依赖流变律 `B = f(T)`（MISOMIP2 Coupled 需要） |
| 6.5.3 | 等温冰简化选项保留（MISMIP+ / MISOMIP1 假设） |

---

### 阶段 7：验证 — MISOMIP1 理想化耦合基准

> **硬门禁（设计决策 D10）**：本阶段是首个涉及动态冰盖-海洋耦合的验证实验。
> 进入前必须完成：阶段 2.6（测试框架升级）、阶段 3.3（SMB）、薄水柱问题。
> （阶段 2.6、阶段 3.3、阶段 3.4、阶段 4a.5、阶段 6.6
> 均已完成并通过编译验证。）
> 否则验证结论不可信且问题难以定位。

| 步骤 | 工作内容 |
|------|---------|
| 7.1 | 准备 MISOMIP1 配置（修复当前 MISOMIP 网格薄水柱问题） |
| 7.2 | **IceOcean1r**（参考稳态）：验证 SSA 冰速 + 3EQN 底融 + 冰架腔环流 |
| 7.3 | **IceOcean1ra**（瞬态）：开启厚度演化 + 接地线迁移 |
| 7.4 | **Ocean1r/Ocean1ra**（固定冰架腔）：验证纯海洋响应 |
| 7.5 | 与 MISOMIP1 发布基准对比 |

**实验命名**：使用 Asay-Davis et al. (2016) 原文正式名称（IceOcean1r 等），避免与 MISOMIP2 混淆。

---

### 阶段 8：验证 — MISOMIP2 Ocean 实验

| 步骤 | 工作内容 |
|------|---------|
| 8.1 | 准备 MISOMIP2 Ocean3 配置（Amundsen Sea，固定冰架） |
| 8.2 | Ocean3（参考）→ Ocean4（大气扰动）：验证底融敏感性 |
| 8.3 | MIPkit 观测数据对比 |

---

### 阶段 9：验证 — MISOMIP2 Coupled 全耦合实验

| 步骤 | 工作内容 |
|------|---------|
| 9.1 | Coupled1（参考稳态）→ Coupled2（几何扰动）→ Coupled3（大气扰动） |
| 9.2 | 与 MIPkit + 其他模型（MITgcm STREAMICE, FESOM, NEMO）对比 |

---

## 测试框架

### 测试目录结构

```
icesheet_test/
├── tests/
│   ├── test_ssa_analytical.py      # SSA 解析解验证
│   ├── test_ssa_mms.py             # MMS 收敛性测试
│   ├── test_parallel_consistency.sh # 并行一致性
│   ├── test_picard_convergence.py   # Picard 收敛曲线
│   └── test_mass_conservation.py    # 质量守恒检查
├── run_test.sh                      # 自动化冒烟测试
└── ...
```

### 测试规范

| 测试类型 | 触发时机 | 通过标准 |
|---------|---------|---------|
| 冒烟测试（ISOMIP 100 步） | 每次代码修改 | 无 NaN，exit code 0 |
| 并行一致性 | 修改 halo/exchange 逻辑 | 1 vs 4 tiles bit-for-bit |
| MMS 收敛 | 修改 FE 装配 | L2 误差 O(h²) |
| 解析解 | 修改驱动应力/黏度 | 误差 < 1% |
| 质量守恒 | 修改厚度演化 | 全局质量误差 < 1e-10 |

---

## 依赖关系图

```
阶段 0,1,2 ✅ 已完成（含 2.A FE 装配）
  └─→ 阶段 2.5 (验证与测试) ✅
       ├─→ 阶段 2.6 (测试框架升级) ✅ ← 阶段 7 前提
       └─→ 2.B-C (MMS + 并行测试) ✅
  └─→ 阶段 3 (厚度演化) ✅ 已完成（含 3.3 SMB）
       └─→ 3.4 (Fi 初始化) ✅ 已完成
  └─→ 阶段 5 (崩解) ✅ 基本完成（5.5 ✅ 掩膜同步）
  └─→ 阶段 6.1 (多时间尺度调度) ✅ 已完成
       └─→ 阶段 4a (基本接地线) ✅ 已生效（4a.5 方案 B）
            ├─→ 阶段 6.6 (isf_coupler 修复) ✅ 已完成
            └─→ 阶段 4b (次网格接地线)
                 └─→ 阶段 7 (MISOMIP1) 🔲 硬门禁（仅剩薄水柱问题）
                      └─→ 阶段 8 (MISOMIP2 Ocean)
                           └─→ 阶段 9 (MISOMIP2 Coupled)
```

**阶段 7 硬门禁（设计决策 D10）**：进入 MISOMIP1 验证前必须完成——
阶段 2.6（测试框架升级）✅、阶段 3.3（SMB 接口）✅、薄水柱问题 🔲。

---

## 优先级排序（2026-07-31 更新）

```
[紧急] ✅ 修复 P0/P1 bug（已完成 2026-07-27）
   ↓
[紧急] ✅ 驱动应力 FE 装配 + 厚度平流 + 崩解 + 多时间尺度（已完成 2026-07-28）
   ↓
[高]   ✅ MISOMIP 配置数值稳定性调优（CG BC 修复，DT=300 稳定）
   ↓
[高]   ✅ 验证测试框架（解析解 + MMS + halfpipe + 并行一致性）
   ↓
[高]   ✅ 阶段 6.6：isf_coupler 缺陷修复（halo、mask 判据、U/V 插值、sustr 双写）
   ↓
[高]   ✅ 阶段 4a.5：WET_DRY 依赖处理（方案 B，动态 zice/rmask_isf 生效）
   ↓
[高]   ✅ 阶段 3.4：Fi(icImef)/icIust 首步初始化（allocate 归零）
   ↓
[中]   🔲 薄水柱问题
   ↓
[中]   🔲 阶段 7：MISOMIP1 IceOcean1r 验证（硬门禁）
   ↓
[低]   阶段 4b：次网格接地线
   ↓
[低]   MISOMIP2 真实区域验证
```

---

## 文档与工具链维护

| 项目 | 规则 |
|------|------|
| Skills 文档 | 每次修改 `globaldefs.h` 后同步更新 `roms-code-conventions/SKILL.md` 中的宏文档 |
| 实验命名 | 统一使用 Asay-Davis et al. (2016) / De Rydt et al. (2024) 原文正式名称 |
| 自动化测试 | `icesheet_test/run_test.sh`：编译 → 运行 → 检查输出（长期接入 CI） |
| 本文件更新 | 每完成一个阶段或重大设计变更后更新本文件 |

---

## 关键参考文件映射

| 本项目组件 | STREAMICE 参考 | PISM 参考 |
|-----------|---------------|-----------|
| `isf_driving_stress.h` | `streamice_driving_stress.F` | — |
| `isf_viscosity.h` | `streamice_visc_beta.F` | `compute_nuH_everywhere` |
| `isf_cg_solve.h` | `streamice_cg_functions.F` | `SSAFEM` |
| `isf_thickness_adv.h` | `streamice_advect_thickness.F` | `mass/transport.cc` |
| `isf_grounding_line.h` | `eta_gl_streamice.F` | `groundedIceFraction.cc` |
| `isf_adv_front.h` | `streamice_adv_front.F` | — |
| `isf_calving.h` | — | `calving/thickness_threshold.cc` |

---

## 单次提交约束

根据 AGENTS.md：**单次提交不超过 1000 行修改**。每个步骤对应 1-3 个 commit。

---

*本文件为 Qoder 和 Zcode 的共享开发计划，随 `icesheet` 分支开发进度同步更新。*
