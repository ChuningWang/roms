# ROMS 冰盖-海洋耦合模块开发工作计划 (Work Plan)

> **本文件是所有 AI agent（Qoder、Zcode 等）的共享开发计划。**
> 修改本文件前请确认与其他 agent 的工作不冲突。
>
> 基于 `AGENTS.md` 项目背景，参考 MITgcm STREAMICE 与 PISM 架构。
>
> **创建日期**：2026-07-10  
> **最后更新**：2026-07-28（阶段 2.5 全部完成：CG BC 修复 + 解析解/MMS/halfpipe 验证）  
> **分支**：`icesheet`  
> **约束**：单次 commit 不超过 1000 行修改

---

## 相关文档索引

| 文档 | 相对路径 | 用途 |
|------|---------|------|
| SSA 技术设计 | `docs/plans/ssa_solver_design.md` | FE 离散公式、STREAMICE 对比分析、参数体系 |
| 代码审查报告 | `docs/reviews/code_review_2026-07.md` | 20 项问题及修复状态 |
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
| `isf_vbc.F` | **完整实现** | 3EQN 热力学，计算底融速率 |
| `isf_masks.F` | **部分实现** | 掩膜计算可用，halo exchange 待修复 |
| `iceshelf.F` | **完整驱动 + 多时间尺度** | 含 SSA/厚度/崩解/接地线调度（N_ssa/N_thick/N_gl） |
| `mod_isf.F` / `isf_mod.h` | **完整基础设施** | 状态变量(5)、诊断场(7+6)、参数(15+6新增)、时间层(2) |
| CPP 宏体系 | **`ICESHELF_*` + `ICESHEET_SSA` + `ICESHEET_CALVING`** | SSA→MOM→DYNAM 派生 |
| 厚度演化 | **已实现** | CFL 子步进 + Van Leer 限制器 + 底融源项 |
| 接地线动力学 | **基本浮力判据 + 多时间尺度** | 余弦平滑 float_frac，N_gl 步更新一次 |

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
- [ ] MISOMIP 网格薄水柱问题（57 个 cell 水柱 < 10m）
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

---

### 阶段 3：冰盖厚度演化（质量守恒）✅ 已完成

**目标**：`∂H/∂t = -∇·(ūH) + ȧ - ḃ`

| 步骤 | 文件 | 工作内容 | 状态 |
|------|------|---------|------|
| 3.1 | `isf_thickness_adv.h` | 通量形式厚度平流 + Van Leer 限制器 + CFL 子步进 | ✅ |
| 3.2 | `iceshelf.F` | 多时间尺度调度（N_thick 步更新一次） | ✅ |
| 3.3 | 新建 `isf_smb.h` | 表面质量平衡接口（MISOMIP1 用固定 SMB 强迫场） | 待实施 |

**PSI→面插值说明**（设计决策 D2）：
- 厚度平流通量需要 RHO 面（x-face/y-face）上的速度
- 使用简单 2 点平均：`u_xface(i,j) = 0.5*(u_psi(i,j) + u_psi(i,j+1))`
- 这是一阶精度，但对 MISOMIP1 分辨率（2 km）足够
- 如后续需要更高精度，可改为 4 点加权插值

---

### 阶段 4：接地线动力学

分为两个子阶段：

#### 阶段 4a：基本接地线（快速可用）

| 步骤 | 文件 | 工作内容 |
|------|------|---------|
| 4a.1 | 新建 `isf_grounding_line.h` | 浮力判据：`ρ_i*H > ρ_w*h` → 接地（已实现，需从 isf_mom.F 提取） |
| 4a.2 | 同上 | 余弦平滑 `float_frac`：在接地线宽度 `W=3Δx` 内从 0→1 过渡 |
| 4a.3 | `isf_masks.F` | 重写掩膜更新：厚度更新后重新计算接地/漂浮。修复 halo exchange |
| 4a.4 | `isf_masks.F` / `set_depth.F` | zice 动态更新：漂浮 `zice = -(ρ_i/ρ_w)*H`，接地 `zice = 0` |

#### 阶段 4b：次网格接地线（精度提升）

| 步骤 | 文件 | 工作内容 |
|------|------|---------|
| 4b.1 | `isf_grounding_line.h` | 次网格接地线位置插值：参考 PISM 的 LI（linear interpolation）方法，在相邻接地/漂浮单元间线性插值精确接地线位置 |
| 4b.2 | 同上 | 接地线迁移速率计算：`v_gl = Δx_gl / Δt`，用于 MISMIP+ 验证 |
| 4b.3 | 验证 | 网格收敛性测试（1 km vs 2 km vs 4 km），确认接地线位置精度 |

---

### 阶段 5：崩解与冰前动力学 ✅ 基本完成

| 步骤 | 文件 | 工作内容 | 状态 |
|------|------|---------|------|
| 5.1 | `isf_calving.h` | **厚度阈值崩解**（D6）：`H < H_calve` 移除冰，`H_calve = 100 m` | ✅ |
| 5.2 | 同上 | 冰前检测：`H>0` 但某邻居 `H=0` → 冰前单元 | ✅ |
| 5.3 | 新建 `isf_adv_front.h` | 冰前推进（可选）：参考 `streamice_adv_front.F` 填格子方案 | 待实施 |
| 5.4 | `iceshelf.F` | 调度：厚度更新 → 崩解 → 冰前更新 | ✅ |

---

### 阶段 6：冰盖-海洋耦合增强 + 多时间尺度

#### 6.1 多时间尺度策略（设计决策 D5）

```
海洋时间步:     dt_ocean = 240 s (每个 baroclinic step)
SSA 求解:       每 N_ssa 个海洋步求解一次 (默认 N_ssa = 10)
厚度演化:       每 N_thick 个 SSA 步更新一次 (默认 N_thick = 10)
接地线更新:     每 N_gl 个厚度步更新一次 (默认 N_gl = 10)
```

参数从 `isf.in` 读取。暂不实现 adaptive time stepping。

#### 6.2 耦合通量

| 步骤 | 文件 | 工作内容 |
|------|------|---------|
| 6.2.1 | `isf_vbc.F` / `set_depth.F` | 动态 zice 更新 + set_depth 联动 |
| 6.2.2 | `isf_mom.F` / `step2d_FB.h` | 冰架底部拖曳 → 海洋 `bustr`/`bvstr`（PSI→U/V 插值） |
| 6.2.3 | 新建 `isf_coupler.h` | 统一耦合器：冰→海洋（淡水、热、盐、动量）；海洋→冰（底融） |
| 6.2.4 | `isf_coupler.h` | 底融速率时间滤波（指数移动平均，τ ~1 天）+ 通量守恒检查 |
| 6.2.5 | `iceshelf.F` | 多时间尺度调度逻辑 |

---

### 阶段 6.5：冰盖热力学（可选扩展）

| 步骤 | 工作内容 |
|------|---------|
| 6.5.1 | 冰盖温度场演化方程（参考 PISM enthalpy formulation） |
| 6.5.2 | 温度依赖流变律 `B = f(T)`（MISOMIP2 Coupled 需要） |
| 6.5.3 | 等温冰简化选项保留（MISMIP+ / MISOMIP1 假设） |

---

### 阶段 7：验证 — MISOMIP1 理想化耦合基准

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
  └─→ 阶段 2.5 (验证与测试) ✅ 已完成
       └─→ 2.B-C (MMS + 并行测试) ✅
  └─→ 阶段 3 (厚度演化) ✅ 已完成
       └─→ 阶段 5 (崩解) ✅ 基本完成
  └─→ 阶段 6.1 (多时间尺度调度) ✅ 已完成
       └─→ 阶段 4a (基本接地线) 🔄 部分完成
            ├─→ 阶段 6.2-6.5 (耦合增强)
            └─→ 阶段 4b (次网格接地线)
                 └─→ 阶段 7 (MISOMIP1)
                      └─→ 阶段 8 (MISOMIP2 Ocean)
                           └─→ 阶段 9 (MISOMIP2 Coupled)
```

---

## 优先级排序（2026-07-28 更新）

```
[紧急] ✅ 修复 P0/P1 bug（已完成 2026-07-27）
   ↓
[紧急] ✅ 驱动应力 FE 装配 + 厚度平流 + 崩解 + 多时间尺度（已完成 2026-07-28）
   ↓
[高]   ✅ MISOMIP 配置数值稳定性调优（CG BC 修复，DT=300 稳定）
   ↓
[高]   ✅ 验证测试框架（解析解 + MMS + halfpipe + 并行一致性）
   ↓
[中]   阶段 4a：接地线动力学完善（halo exchange 修复）
   ↓
[中]   MISOMIP1 IceOcean1r 验证
   ↓
[低]   阶段 4b：次网格接地线
   ↓
[低]   阶段 6：耦合增强（isf_coupler.h、动量交换）
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
