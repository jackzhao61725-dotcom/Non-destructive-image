# Dual-port Faraday imaging：实验室测量执行计划

**版本：** 2026-07-19
**用途：** 从尚未定型的光路开始，逐步完成 dual-port Faraday 成像、仪器标定、单帧测量、重复曝光测量和模型反馈。
**适用边界：** 本文件是实验执行清单，不是 dissertation 正文。所有未测参数均保持 `TBD`，不得用 Chapter 4–5 的 simulation reference values 冒充实验值。

---

## 0. 当前已确定与尚未确定的事项

### 0.1 已确定的设计方向

- principal readout 是 **dual-port Faraday imaging**；
- analyser 使用 **half-wave plate (HWP) + Wollaston prism (WP)**，不是 conventional PBS；
- WP 后的两个偏振输出由同一 tube lens 成像到同一 camera 的两个独立 ROI；
- camera 在 resonant absorption imaging (RAI) 与 Faraday imaging 之间原则上不移动；
- resonant 与 far-detuned probe 由不同 fibre input 提供，并在进入共同成像路径前切换；
- 两端口原始图像 \(H\) 和 \(V\) 必须分别保存，不能只保存处理后的 \(S\)；
- Chapter 5 的 \(|\Delta|/2\pi=1.5\ \mathrm{GHz}\) 和 fluence points 只作为首轮实验的参考，不是硬件规格；
- repeated-exposure 测量首先约束 **net disturbance**，不试图从同一条 depletion curve 中分别反演 recoil heating、reabsorption 和 \(\eta_{\mathrm{coll}}\)。

### 0.2 尚未确定，必须通过本计划解决的事项

| 项目 | 当前状态 | 在何处确定 |
|---|---|---|
| 实际 magnification | `TBD` | T1、T7 |
| installed effective NA / PSF | `TBD` | T7 |
| Wollaston 型号、分束角和安装位置 | `TBD` | T2 |
| 两端口在 camera 上的实际间距与 ROI 尺寸 | `TBD` | T2 |
| H/V 端口命名和 \(S\) 的正号约定 | `TBD` | T4、T10 |
| resonant 与 detuned beams 在原子处的重合度 | `TBD` | T3 |
| WP 保留时能否通过 \(H+V\) 完成 RAI | `TBD` | T3 |
| 同一 realization 中 Faraday 后接 terminal RAI 是否可行 | `TBD` | T3 |
| 实际 line centre 与 detuning uncertainty | `TBD` | T5 |
| 原子处功率、beam radii 和 pulse energy | `TBD` | T5、T6 |
| camera offset、read noise、gain 和 linear range | `TBD` | T8 |
| input polarisation、port balance 与背景漂移 | `TBD` | T4、T9 |
| reference condensate distribution | `TBD` | T11 |
| 有效 Faraday coefficient \(\kappa_F\) | `TBD` | T13 |
| net repeated-exposure disturbance | `TBD` | T15 |

---

## 1. 总体执行顺序与硬性 gate

实验按依赖关系推进。只有前一 gate 通过，才进入下一阶段的大规模数据采集。

| 阶段 | 任务 | Gate |
|---|---|---|
| **A. 光路定义** | T1–T4 | 两端口可同时、完整、稳定地落在同一 camera；resonant/detuned 路径兼容；背景可控制 |
| **B. Probe 与 timing** | T5–T6 | 实际 \(\Delta\)、\(P_{\mathrm{atoms}}\)、\(\tau_{\mathrm{eff}}\) 和 pulse energy 可追溯且稳定 |
| **C. Detector 与空间响应** | T7–T9 | counts、noise、port registration 和空间采样能够定量解释 |
| **D. Atomic baseline 与首图** | T10–T12 | 获得可重复、与 condensate 空间位置一致的 atom-dependent dual-port signal |
| **E. 定量标定与 operating region** | T13–T14 | 得到有效 \(\kappa_F\) 和实验 single-frame SNR–fluence relation |
| **F. Repeated exposure 与模型评估** | T15–T16 | 得到 net disturbance，并完成至少一组 frozen-parameter held-out comparison |

> **重要：** “先得到首图”不要求提前完成所有高精度标定。T7–T9 可以先完成足以支持首图的 provisional version，随后再做精确版本；但任何进入 dissertation 的绝对数值必须来自精确版本。

---

# Stage A — 光路定义

## T1. Survey the existing imaging path

### 目标

把现有 RAI 成像路径转换成一张可测量、可复现的光路图，确认 WP、HWP 和双端口传播可以实际放置的位置。此任务不依赖原子。

### 所需材料

- 现有 optical table 和 imaging path；
- 元件型号记录或照片；
- ruler/caliper、beam card、low-power alignment beam；
- camera live view；
- 当前 optical layout SVG 的打印件或电子版；
- lab laser-safety SOP。

### 操作步骤

1. 从 condensate/object plane 开始，沿传播方向逐个记录：
   - viewport；
   - objective；
   - relay lenses；
   - mirrors；
   - filters；
   - aperture/iris；
   - camera window 和 sensor plane。
2. 对每个元件记录：
   - 型号；
   - nominal focal length；
   - clear aperture；
   - coating/wavelength range；
   - 当前位置 \(z\)；
   - 是否可移动；
   - 是否与其他实验共用。
3. 标出所有 approximately collimated sections、pupil planes 和 image planes。不要仅凭图纸判断，使用 beam size 随 \(z\) 的变化或现有成像关系确认。
4. 测量 camera 前可用直线传播距离、横向安装空间和 beam-height constraints。
5. 在 camera 上记录当前单束图像的：
   - centroid；
   - full illuminated area；
   - condensate expected ROI；
   - 可用于 background estimation 的周围区域；
   - 到 sensor 四边的余量。
6. 记录当前 magnification 的来源。如果只是 lens ratio 或旧记录，标为 nominal，不写成 measured。
7. 更新光路图，严格区分：
   - installed；
   - available but not installed；
   - proposed；
   - position not yet fixed。

### 必须记录

| 字段 | 记录值 |
|---|---|
| Objective model / focal length / NA | `TBD` |
| Relay-lens models and positions | `TBD` |
| Camera model and active sensor size | `TBD` |
| Current nominal magnification | `TBD` |
| Available collimated length for WP | `TBD` |
| Distance from candidate WP position to camera | `TBD` |
| Maximum usable horizontal/vertical sensor span | `TBD` |
| Shared components that cannot be moved | `TBD` |

### 通过标准

- 能够从 object plane 一直追踪到 camera，不存在未识别的 relay stage；
- 至少确定一个 physically accessible WP candidate position；
- 已知 camera 上可容纳双端口图像的最大范围；
- 没有把 nominal magnification 或 simulated NA 当成 measured value。

### 交付物

- `optical_path_as_installed.svg/pdf`；
- `optical_path_inventory.csv`；
- 有比例尺的 table photograph；
- T1 run note。

---

## T2. Determine the Wollaston dual-port geometry

### 目标

确认所选 WP 能使 \(H\) 和 \(V\) 两幅像完整落在同一 camera 上，且不重叠、不 clipping、焦点可同时接受。

### 前置条件

- T1 完成；
- WP 型号或至少候选型号已知；
- low-power alignment beam 可进入共同成像路径。

### 预计算

从 WP datasheet 记录 angular separation \(\alpha_{\mathrm{WP}}\)。若 WP 位于近似 collimated section，camera-plane separation 的初始估计为

\[
d_{\mathrm{est}}\simeq f_{\mathrm{eff}}\alpha_{\mathrm{WP}},
\]

其中 \(f_{\mathrm{eff}}\) 是 WP 后将角度转换为位移的有效焦距。这个关系只用于选初始位置，最终值必须实测。

### 操作步骤

1. 在低功率、无原子条件下安装 HWP 和 WP。
2. 调整 HWP，使两端口功率大致相等；此时不要追求最终 balance calibration。
3. 在 camera 上识别两个输出，并记录：
   - centroid separation；
   - horizontal/vertical orientation；
   - beam/image size；
   - relative focus；
   - sensor-edge clearance。
4. 平移 WP 或调整 relay distance，找到同时满足以下条件的位置：
   - 两个 atom ROI 不重叠；
   - 每个 ROI 周围保留 background area；
   - 无 visible aperture clipping；
   - 两端口都处于可接受焦点；
   - expected pointing drift 不会把图像推到 sensor 边缘。
5. 扫描 camera focus 或 tube-lens position，分别得到两个端口的 sharpness-versus-position curve。
6. 轻微改变 input pointing，测量两个端口的共同移动和相对移动，确认共路传播是否保持 port registration。
7. 暂时固定 WP、HWP、tube lens 和 camera。记录 mount positions，避免后续每次重调改变几何关系。

### 建议 acquisition

- 每个候选 WP position：至少 10 张无原子双端口图；
- focus scan：覆盖共同 best focus 两侧，每点至少 5 张；
- pointing repeatability：至少 20 次重复或连续记录。

### 分析

对每张图拟合两端口 centroid 和 second-moment size，输出：

- measured separation \(d\)；
- port sizes；
- focus difference；
- centroid drift；
- sensor margins；
- clipping indicator，例如 total power 对小幅 pointing change 是否非对称变化。

### 通过标准

- \(H/V\) 两幅像在完整 expected cloud extent 内不重叠；
- 每幅图像周围均有可用于背景估计的像素；
- 两端口无明显 clipping；
- 两端口可使用同一 camera focus；
- port separation 和 centroid drift 小于预留 margin；
- 选定几何被记录为 measured configuration，而不是仅靠 SVG 尺寸。

### 失败时的处理顺序

1. 改变 WP 到 camera 的有效传播几何；
2. 改变 tube-lens focal length 或 relay arrangement；
3. 调整图像 magnification；
4. 更换不同 separation angle 的 WP；
5. 最后才考虑双 camera。不要在前四项未排除前默认需要第二台 camera。

---

## T3. Test resonant/detuned path compatibility and the RAI route

### 目标

确定 resonant 与 far-detuned fibre inputs 是否在原子处共享同一空间模式，并决定 terminal RAI 的实际实现方式。

### 必须回答的核心问题

1. 两个 fibre outputs 切换后，beam centroid、angle、waist 和 polarisation 是否保持一致？
2. 是否可以在 WP 保留的情况下，用两个端口之和完成 RAI？
3. Faraday pulse 后能否在同一 realization 中自动切换到 resonant pulse 并记录 terminal RAI？

### 操作步骤 A：两束光的空间重合

1. 依次选择 resonant 和 detuned fibre input，不移动 downstream imaging path。
2. 在至少两个相距较远的诊断平面记录 beam centroid 和 radius。
3. 由两个平面的 centroid difference 分离 lateral offset 与 angular mismatch。
4. 在 atom-equivalent plane 或通过 relay 映射到该平面，记录两束光的 beam profile。
5. 反复切换输入，测量 switching repeatability。
6. 调整 fibre output coupler 或 combining optics，直到两束光的 offset、angle 和 waist difference 满足 condensate illumination 要求。

### 操作步骤 B：验证 WP-on absorption readout

1. 保持 HWP、WP、tube lens 和 camera 不动。
2. 对 resonant path 采集两端口 probe frames 和 dark frames。
3. 做 port gain correction 后构造

\[
C_{\Sigma}=C_H+C_V.
\]

4. 扫描 HWP 小角度范围，检查 \(C_{\Sigma}\) 是否独立于 H/V 功率重新分配。
5. 改变输入总功率，检查 \(C_{\Sigma}\) 的线性和稳定性。
6. 如果可获得原子，分别比较：
   - WP 保留时由 \(H+V\) 得到的 absorption image；
   - 当前标准 RAI configuration 得到的 absorption image。
7. 比较 atom number、radii、column-density profile 和 background residuals。

### Decision A：WP-on RAI 通过

若 \(H+V\) 在 gain correction 后能够稳定复现 RAI，则保留 WP，并设计同一 realization 的时序：

1. detuned Faraday pulse；
2. 保存 \(H/V\) frame；
3. 切换到 resonant fibre input；
4. terminal absorption pulse；
5. 保存 resonant \(H/V\) atom image；
6. 取得相应 probe 和 dark references。

### Decision B：WP-on RAI 不通过

若必须物理移除 WP，而移除过程不能在单次实验序列中自动完成，则：

- 不得写“same-realisation paired Faraday–RAI calibration”；
- 使用 matched ensembles：Faraday shots 与 RAI-only shots 采用相同 preparation window；
- 将 condensate shot-to-shot variation 明确计入 \(\kappa_F\) uncertainty；
- 优先评估可重复翻转/平移 mount 或其他不移动 camera 的切换方案。

### 通过标准

- 两 fibre paths 的重合度和切换重复性已测量；
- RAI 的实际 downstream configuration 已决定；
- “same-realisation”或“matched-ensemble”路线已明确，不能在后续分析中混用；
- 对应 sequence timing 在 control system 中可实现。

---

## T4. Polarisation preparation, port convention and background

### 目标

建立 reproducible linear input polarisation，确定 HWP working point、H/V port naming、\(S\) 的符号约定和无原子 polarisation background。

### 操作步骤

1. 在进入 atoms 前的可测位置检查 fibre output polarisation。
2. 安装 input polarisation preparation optics，并记录 polariser/HWP 的机械零点。
3. 无原子条件下扫描 HWP angle \(\psi\)，同时记录两个端口总 counts。
4. 对每个角度采集多帧，拟合两端口功率随 \(\psi\) 的互补变化。
5. 选择 \(H/V\) 平衡点，使无旋转时两个端口的 integrated counts 在 gain correction 后相等。
6. 在平衡点记录长时间无原子序列，得到：
   - mean \(S_0(y,z)\) map；
   - frame-to-frame rms；
   - spatial background structure；
   - integrated imbalance；
   - drift versus time。
7. 改变 bias-field direction 或使用一个已知 polarisation rotation，确定哪一个端口定义为 H、哪一个定义为 V，以及 positive \(\theta_F\) 对应 \(S\) 的正负号。
8. 在 optical path 不变的情况下重复上述步骤，判断 viewport/relay birefringence 是否产生稳定 offset 或空间 pattern。

### 建议 acquisition

- HWP scan：覆盖至少一个完整 balance crossing，每点 10–20 帧；
- balance-point stability：至少覆盖一次典型实验 run 的持续时间；
- run 开始和结束各采集一组 no-atom H/V references。

### 端口校正定义

暗计数扣除后，定义

\[
C_H=H-D_H,\qquad C_V=\beta\,(V-D_V),
\]

其中 \(\beta\) 由无原子平衡数据确定。随后计算

\[
S=\frac{C_H-C_V}{C_H+C_V},\qquad \Delta S=S_{\mathrm{atoms}}-S_0.
\]

不要在每张 atom frame 上重新拟合 \(\beta\)，否则会吸收真实 atom-dependent signal。\(\beta\) 应由独立 calibration block 冻结。

### 通过标准

- H/V port labels 和 \(S\) sign convention 已写入配置；
- balance angle 可重复找到；
- no-atom background map 和技术噪声已量化；
- background drift 在首个候选 atom ROI 内小于预期或已观测 atom-dependent signal；
- 若背景不稳定，先处理 polarisation/path drift，不进入 \(\kappa_F\) 标定。

---

# Stage B — Probe generation and timing

## T5. Reconstruct and measure the optical frequency chain

### 目标

用实际 laser lock point、AOM frequencies、diffraction orders 和 pass numbers 确定原子处的 probe detuning，而不是直接把 commanded value 当成 \(|\Delta|/2\pi=1.5\ \mathrm{GHz}\)。

### 频率账本

对每个 frequency-shifting element 记录：

| Element | RF frequency | Diffraction order \(m_j\) | Pass number \(n_j\) | Shift sign | Optical shift \(m_jn_jf_j\) |
|---|---:|---:|---:|---:|---:|
| Laser lock / offset | `TBD` | — | — | `TBD` | `TBD` |
| AOM 1 | `TBD` | `TBD` | 1 or 2 | `TBD` | `TBD` |
| AOM 2 | `TBD` | `TBD` | 1 or 2 | `TBD` | `TBD` |
| Additional offset | `TBD` | — | — | `TBD` | `TBD` |

计算

\[
\nu_{\mathrm{probe}}=\nu_{\mathrm{lock}}+\sum_j m_jn_jf_j,
\qquad
\Delta=\nu_{\mathrm{probe}}-\nu_0.
\]

### 操作步骤

1. 从 laser reference/lock 开始逐级追踪 frequency chain。
2. 记录每个 AOM 的 RF source、frequency readback、order sign 和 single/double-pass geometry。
3. 记录 fibre input 切换前后的 frequency path，避免把 resonant 与 detuned branch 的 offset 混合。
4. 使用 lab 可用的 spectroscopy、beat-note、wavemeter 或 RF-referenced method 确定 line centre 和 optical offset。选择的方法及其 uncertainty 必须记录。
5. 在一次典型 run 内连续记录 frequency reference，测量 short-term drift。
6. 在多天重复 reference measurement，区分 within-run 与 day-to-day drift。

### 输出

- measured \(\Delta/2\pi\)；
- statistical 和 systematic uncertainty；
- stable operating range；
- run-to-run recalibration rule；
- 完整 frequency-chain diagram。

### 通过标准

- 任何 shot 的 detuning 可以由 metadata 重建；
- 所有 AOM order signs 明确；
- 1.5 GHz 是 measured setting with uncertainty，而不是只有 nominal command；
- detuning drift 相对于 planned scan spacing 足够小。

---

## T6. Measure delivered power, beam profile and pulse timing

### 目标

确定每个 imaging pulse 在原子处实际提供的 peak intensity 和 integrated energy，并确认 camera exposure 与 optical pulse 的时序关系。

### 操作步骤 A：功率链路

1. 在原子前最后一个可访问位置测量 optical power。
2. 测量从该位置到 atom plane 的 transmission；若不能直接测 atom plane，记录每个 downstream loss 和总 transmission uncertainty。
3. 在所有候选 AOM RF amplitudes 下测量 delivered power。
4. 重复 power ramp 上行和下行，检查 hysteresis 或 thermal drift。
5. 在典型 run 时间内记录 power drift。

### 操作步骤 B：beam profile

1. 在 atom-equivalent plane 测量 2D beam profile。
2. 拟合 beam centre、\(w_y\)、\(w_z\)、ellipticity 和 local intensity gradient。
3. 将 expected condensate position 映射到 beam profile，得到 atoms-at-centre intensity。
4. 检查不同 fibre input 和不同 AOM frequencies 是否改变 pointing 或 beam size。

### 操作步骤 C：pulse timing

1. 使用 fast photodiode 同时记录 AOM gate、实际 optical pulse 和 camera trigger/exposure marker。
2. 对每个候选 commanded duration 测量：
   - trigger delay；
   - rise time；
   - fall time；
   - ringing；
   - pulse-to-pulse jitter；
   - integrated pulse area。
3. 定义实际有效时长和 pulse energy：

\[
E_{\mathrm{pulse}}=\int P(t)\,dt,
\qquad
F_{\mathrm{eff}}\propto E_{\mathrm{pulse}}
\]

其中 simulation coordinate \(F=P\tau\) 只有在矩形 pulse 近似成立时才等于 measured pulse area。
4. 确认 optical pulse 完整落在 camera exposure 内，并记录 timing margin。
5. 测量 shortest reliable pulse，不以 command generator 的最短 setting 代替实测结果。

### 建议 acquisition

- 每个 candidate duration：至少 50 个 photodiode traces；
- 每个 candidate power：run 开始、中央和结束各测一次；
- beam profile：至少两个不同日期重复。

### 输出

- \(P_{\mathrm{atoms}}\) 与 uncertainty；
- \(w_y,w_z\)、beam centre 与 local intensity；
- \(\tau_{\mathrm{eff}}\)、pulse area 和 jitter；
- camera–probe timing diagram；
- 可实现的 \(P\)-\(\tau\) region。

### 通过标准

- 每个 shot 的 \(P_{\mathrm{atoms}}\)、\(\tau_{\mathrm{eff}}\) 和 pulse area 可追溯；
- 候选 pulse 不存在严重 clipping、ringing 或 exposure-window loss；
- 使用 measured pulse area，而不是只使用 programmed \(P\tau\)。

---

# Stage C — Detector and spatial response

## T7. Measure magnification, port registration and effective spatial response

### 目标

用 installed system 的实测值替换 nominal magnification 和 simulated ideal aperture，并建立两个端口之间固定的几何映射。

### Magnification measurement

优先使用以下方法之一，选择实际可实现且 uncertainty 最清楚的方法：

1. object-equivalent plane 的 calibrated target；
2. 已知距离的 object/stage translation；
3. 经过独立标定的 atomic displacement。

对两个方向分别计算

\[
M_y=\frac{\Delta y_{\mathrm{camera}}}{\Delta y_{\mathrm{object}}},
\qquad
M_z=\frac{\Delta z_{\mathrm{camera}}}{\Delta z_{\mathrm{object}}}.
\]

不要强制令 \(M_y=M_z\)。记录 rotation、skew 和 uncertainty。

### Port registration

1. 使用相同 target 或稳定 spatial pattern 同时记录 H/V。
2. 拟合从 V 到 H 的 affine transform：translation、rotation、scale；只有数据确实要求时才加入 shear。
3. 使用独立图像检查 registration residual。
4. 将 transform 固定在 calibration file 中，不允许逐张 atom image 自由拟合。

### Spatial response

1. 若有 calibrated edge，测 edge-spread function，并求 line-spread function/MTF。
2. 若有 point-like or narrow feature，直接测 PSF。
3. 若只能使用已知 atomic feature，应明确这是 effective response estimate，而不是独立 optical PSF。
4. 分别检查 H/V 两端口的 response；若一致，可使用共同 response；若不一致，模型中分别保留。

### 输出

- \(M_y,M_z\) 与 object-plane pixel size；
- port registration transform；
- measured effective PSF/MTF；
- focus position 和 uncertainty；
- 对 Chapter 4 ideal circular pupil approximation 的适用范围判断。

### 通过标准

- 同一结构在 H/V registration 后的 residual 小于需要解析的最小 condensate feature；
- magnification 有实测 uncertainty；
- 后续 \(\kappa_F\) fitting 使用 measured spatial response 或明确标记 remaining approximation。

---

## T8. Characterise camera offset, read noise, gain and linearity

### 目标

建立从 raw ADU 到 detector noise 的实验转换，替换 simulation 中 provisional \(\mathrm{QE}\) 和 \(\sigma_r\)。

### 固定 camera configuration

在开始前记录并锁定：

- exposure time；
- trigger mode；
- bit depth；
- analogue/digital gain；
- readout mode；
- ROI；
- temperature/cooling state；
- binning；
- frame rate。

任何一项变化都需要新的 camera calibration identifier。

### Dark frames

1. 遮光并使用与实验相同的 exposure 和 trigger。
2. 采集至少 200 张 dark frames 作为初始基线。
3. 对每个 pixel 计算 mean offset 和 temporal variance。
4. 生成 hot/dead/unstable pixel mask。
5. 分别报告 temporal read noise 和 fixed-pattern offset；不要把未扣除的固定图样当成 read noise。

### Photon-transfer measurement

1. 使用 spatially uniform、stable illumination，覆盖从低 counts 到接近但不进入 saturation 的约 8–12 个 level。
2. 每个 level 采集至少 50 对 frames。
3. 对每对图像作差，以去除固定 illumination pattern；差分 variance 除以 2 得单帧 temporal variance。
4. 拟合 linear region 中 variance 对 mean 的关系。
5. 明确 gain convention：若 \(g\) 的单位为 electron/ADU，则

\[
\operatorname{var}_{\mathrm{ADU}}=\frac{\overline C_{\mathrm{ADU}}}{g}+\sigma_{r,\mathrm{ADU}}^2.
\]

6. 检查 H/V 两个 sensor regions 的 gain、offset 和 read noise 是否一致。
7. 扫描到接近 full-well，但在出现 clipping 前停止，确定 calibrated linear range。

### 关于 QE

- photon-transfer curve 给出 conversion gain，不自动给出 QE；
- 只有当入射到 sensor 的 photon flux 被独立标定时，才能报告 effective QE；
- 若做不到，模型使用 measured effective count scale，不把 catalogue QE 写成实验测量。

### 输出

- offset maps；
- \(\sigma_r\) in ADU and electrons；
- conversion gain；
- linear range；
- pixel mask；
- H/V ROI consistency；
- camera calibration ID。

### 通过标准

- 所有 atom/reference frames 位于 calibrated linear range；
- 两端口 detector differences 已被量化；
- Chapter 5 的 \(3e^-\) provisional value 不再被当成 installed-camera measurement。

---

## T9. Quantify no-atom dual-port noise and common-mode rejection

### 目标

确定实际 dual-port subtraction 是否抑制 probe technical noise，并建立 atom signal 之前的 camera-level noise baseline。

### 操作步骤

1. 使用最终 optical geometry、detuned probe 和 camera configuration。
2. 采集同步 no-atom H/V frames，覆盖典型 run duration。
3. 同时记录 power monitor、detuning monitor 和 relevant trigger metadata。
4. 用 T4 冻结的 \(\beta\)、background map 和 registration transform 计算 \(S_0(y,z)\)。
5. 比较：
   - 单端口 fractional noise；
   - \(H-V\) difference noise；
   - \(S\) noise；
   - \(H+V\) total-intensity noise；
   - H/V covariance。
6. 检查噪声是否随 ROI position、total counts、time 或 HWP angle 变化。
7. 将 measured variance 与 Poisson + read-noise prediction 比较。额外部分记为 technical-noise floor，不强行归入 read noise。

### 输出

- no-atom \(S_0\) mean map；
- variance map；
- port covariance；
- common-mode rejection factor；
- technical-noise floor；
- reference update interval。

### 通过标准

- \(S_0\) 在 atom ROI 内稳定；
- background references 的更新周期已定义；
- 若 technical noise 明显超过 photon + read noise，先定位 source，不把它隐藏在 empirical SNR 中。

---

# Stage D — Atomic baseline and first Faraday image

## T10. Freeze the analysis convention before atom data

### 目标

在看到 atom signal 前固定最小分析流程，防止逐张调参制造信号。

### 必须冻结

- H/V port labels；
- positive \(S\) convention；
- dark subtraction；
- bad-pixel mask；
- port gain factor \(\beta\)；
- registration transform；
- background/reference construction；
- atom ROI 和 background ROI；
- central \(3\times3\) SNR estimator；
- exclusion criteria；
- configuration ID 和 code commit。

### Single-frame processing order

1. 读取 raw H/V；
2. 扣除对应 dark maps；
3. 应用固定 pixel mask；
4. 将 V registration 到 H；
5. 应用固定 port-gain correction；
6. 计算 \(S\)；
7. 减去独立 no-atom \(S_0\)；
8. 在预先定义 ROI 内计算 signal 和 SNR；
9. 保存 processed maps，但不覆盖 raw files。

### 通过标准

- 使用 synthetic/no-atom data 时，pipeline 不产生 condensate-shaped artefact；
- 处理顺序和参数写入 machine-readable config；
- 不允许在每个 atom shot 上重新选择 HWP balance、registration 或 ROI。

---

## T11. Establish the RAI and reference-condensate baseline

### 目标

确定用于 Faraday calibration 的 reference condensate 不是单个假定的 \(N_0\)，而是一个带有 shot-to-shot distribution 和 acceptance window 的实验对象。

### RAI baseline

1. 使用现有 validated RAI pipeline 采集 atom/probe/dark triplets。
2. 在多个 resonant probe intensities 下重复测量，检查 saturation correction 后的 atom number stability。
3. 记录使用的 effective cross-section 或 saturation correction convention。
4. 若 WP-on RAI 已通过 T3，确认 summed-port RAI 与标准 RAI 的一致性。

### Reference-condensate acquisition

1. 使用固定 preparation sequence，连续采集至少 20–30 个 unprobed reference condensates 作为首轮分布；正式统计可扩展到多个日期。
2. 每个 realization 提取：
   - \(N_0\)；
   - \(N_{\mathrm{tot}}\)；
   - temperature；
   - condensate fraction；
   - radii；
   - centre；
   - column-density profile；
   - fit residual。
3. 记录 bias field、optical pumping/spin preparation 和 trap timing identifiers。
4. 根据分布预先定义 calibration acceptance window。不得只保留“看起来好”的 shots 而不记录 rejection reason。

### 输出

- reference-condensate parameter distributions；
- RAI scale uncertainty；
- accepted preparation window；
- shot rejection criteria；
- spin-preparation stability proxy。

### 通过标准

- corrected RAI atom number 在 accepted probe-intensity range 内稳定；
- reference cloud distribution 足以区分 imaging disturbance 与 preparation noise；
- \(\kappa_F\) 标定只使用预先定义的 accepted preparation range。

---

## T12. Obtain the first repeatable dual-port Faraday signal

### 目标

先证明 atom-dependent signal 存在、可重复且空间上与 condensate 一致，再进行大规模 scan。

### 首轮 operating conditions

- detuning 以 measured 1.5 GHz reference 为起点；
- fluence 从实验允许的低-disturbance level 开始；
- Chapter 5 的 \(F=30,50,90,150\ \mathrm{mW\,\mu s}\) 只用作 bracketing reference；
- 不得在 T6 未确认 pulse area 前直接按 commanded \(P\tau\) 使用这些值。

### 每个 condition 的 acquisition block

1. camera dark block；
2. no-atom H/V reference block；
3. 约 10–20 个 accepted atom preparations 的 raw H/V frames；
4. 对应 RAI reference：优先 same-realisation terminal RAI；若 T3 不支持，则使用 interleaved matched RAI shots；
5. no-atom H/V reference block；
6. camera dark block。

### Escalation rule

1. 从最低 planned pulse area 开始；
2. 若没有 repeatable signal，先检查 background、sign、port registration、frequency 和 atom overlap；
3. 只有上述检查通过后才提高 pulse area；
4. 每次提升后检查 terminal RAI 或 matched controls，避免用“更强信号”掩盖明显 disturbance；
5. 任何 saturation、clipping、unstable pulse 或 preparation drift 出现时停止 escalation。

### 首图判据

必须同时满足：

- mean atom-on minus no-atom map 在 condensate expected position 出现结构；
- 结构在重复 shots 中保持相同 sign 和大致 profile；
- signal 不随 camera ROI 或 registration 微调而消失；
- no-atom blocks 中不出现同样结构；
- empirical SNR 高于 background distribution，但此时不强行规定 universal threshold；
- terminal RAI 或 matched RAI 确认该 shot block 中存在 accepted condensate。

### 交付物

- 第一组 repeatable raw H/V atom images；
- corresponding \(S\) and \(\Delta S\) maps；
- no-atom distribution；
- empirical single-frame SNR；
- complete metadata；
- first-signal run report。

---

# Stage E — Quantitative calibration and operating region

## T13. Determine the effective Faraday coefficient \(\kappa_F\)

### 目标

确定在选定 transition、detuning、spin preparation 和 optical configuration 下，将 absorption-derived scalar phase response 映射到 measured Faraday signal 的 **effective coefficient**。该值不是 microscopic erbium calculation。

### Route A：same-realisation paired measurement

仅当 T3 已证明 automatic detuned-to-resonant switching 和 WP-on RAI 可行时使用。

每个 realization：

1. 准备 accepted condensate；
2. 先拍一张低-disturbance Faraday H/V frame；
3. 随后拍 terminal resonant absorption H/V frame；
4. 采集对应 probe and dark references；
5. 保存完整 timing 和 source-switch metadata。

### Route B：matched-ensemble measurement

若不能在同一 realization 做 terminal RAI：

1. 在同一 preparation block 中随机交错 Faraday shots 和 RAI-only shots；
2. 使用 T11 的 accepted preparation window；
3. 比较 ensemble means 和 distributions；
4. 将 shot-to-shot condensate variation 纳入 \(\kappa_F\) uncertainty；
5. 不把它称为 one-to-one paired calibration。

### Acquisition range

- 至少覆盖多个 accepted column densities；
- 首轮保持 detuning 和 spin preparation 不变；
- fluence 保持足够低，使 Faraday pulse 对 terminal RAI 的影响小于所需 calibration precision；
- 若需要检查 linearly, 使用多个 fluence points，但不要同时自由改变 density、detuning、spin state 和 fluence。

### 分析步骤

1. 使用 RAI 得到 \(n_{\mathrm{col}}(y,z)\) 及其 uncertainty。
2. 使用 Chapter 3 的 scalar phase relation 构造 \(\varphi(y,z)\)，并保留 RAI effective-cross-section uncertainty。
3. 使用 T7 measured spatial response 和 T10 fixed dual-port pipeline 生成 forward-model \(S(y,z;\kappa_F)\)。
4. 直接拟合 \(\kappa_F\) 到 measured \(S\) map；小角度关系 \(S\simeq2\theta_F\) 只作为 diagnostic，不替代 finite-aperture fitting。
5. 同时检查 amplitude residual、spatial residual 和 density dependence。
6. 对每个 shot/block 拟合后，判断 \(\kappa_F\) 是否在 accepted density/time range 内稳定。

### 输出

- effective \(\kappa_F\)；
- statistical 和 systematic uncertainty；
- applicable detuning、spin state、density 和 time range；
- spatial residual map；
- calibration route A or B；
- parameter-freeze interval。

### 通过标准

- \(\kappa_F\) 在定义的 operating range 内无显著 density/time drift，或其 dependence 已被显式建模；
- residual 不表现出明显未处理的 registration/PSF/background artefact；
- 不引用 rubidium scale factor 作为 erbium \(\kappa_F\)。

---

## T14. Measure single-frame performance versus fluence

### 目标

用实验数据确定最低能够提供所需 spatial information 的 pulse area，并检验 simulation 中 \(F=P\tau\) 作为主要 scan coordinate 的适用程度。

### 初始 scan design

1. 固定 measured \(|\Delta|/2\pi\approx1.5\ \mathrm{GHz}\)。
2. 选择能够 bracket transition 的 fluence/pulse-area points。优先参考：
   - low：约 \(30\ \mathrm{mW\,\mu s}\)；
   - intermediate：约 \(50\ \mathrm{mW\,\mu s}\)；
   - reference：约 \(90\ \mathrm{mW\,\mu s}\)；
   - high：约 \(150\ \mathrm{mW\,\mu s}\)。
3. 这些数值必须按 T6 measured pulse area 转换；若硬件范围不同，应选择新的 low/intermediate/reference/high points，而不是强行达到模拟值。
4. 每点以约 15–20 个 accepted condensates 作为 pilot；根据 observed variance 决定正式样本量。
5. 随机化 condition 顺序，避免把 run drift 误认为 fluence trend。

### Fixed-fluence power–duration split

在至少一个中间 fluence 上选择两组可实现的 \(P,\tau\)：

- 一组较高 \(P\)、较短 \(\tau\)；
- 一组较低 \(P\)、较长 \(\tau\)；
- 两组 measured pulse area 尽量相同。

比较 signal、SNR、technical noise、pulse reproducibility 和 disturbance。其目的不是重新证明公式，而是检查短曝光是否受到 AOM transient 或 camera timing 限制。

### 每点记录

- raw H/V；
- dark/no-atom references；
- measured \(\Delta,P(t),E_{\mathrm{pulse}}\)；
- \(S\) and \(\Delta S\)；
- central \(3\times3\) SNR；
- spatial profile/edge visibility；
- corresponding condensate reference；
- rejection reason。

### 分析

1. 绘制 empirical single-frame \(\mathrm{SNR}_{3\times3}\) versus measured pulse area。
2. 与 photon + read-noise model 和 measured technical-noise floor 比较。
3. 显示代表性 raw/reconstructed camera images，而不是只报告 SNR 数字。
4. 确定最低满足实际 measurement task 的 pulse area：
   - only detect condensate；
   - resolve centre/elongation；
   - track boundary or shape change。
5. 不声明 universal SNR threshold；记录本实验任务采用的 working image-quality requirement。

### 输出

- experimental SNR–fluence curve；
- representative camera frames；
- working single-frame image-quality criterion；
- initial experimental fluence region；
- fixed-fluence split result。

---

# Stage F — Repeated exposure and model assessment

## T15. Measure net repeated-exposure disturbance

### 目标

测量 repeated Faraday pulses 对 condensate state 和 framewise image quality 的联合影响，得到实际 sequence length。首先约束 net effect，不分别拟合 heating、reabsorption 和 collisional depletion。

### Experimental design

每个 pulse count \(q\) 使用独立 condensate realization，因为 terminal RAI 是 destructive 的。建议首轮从

\[
q=0,1,3,5,8,10,12,15,20,25
\]

中选择能够覆盖无扰动区、transition 和明显 depletion 区的点；pilot 后可调整。不要机械保留全部点。

### Probed sequence

对每个 \(q\)：

1. 准备 accepted condensate；
2. 施加 \(q\) 个相同 Faraday pulses；
3. 保存每一帧 raw H/V；
4. 保持与 planned dynamical sequence 相同的 frame interval；
5. 结束后取得 terminal RAI；
6. 保存 dark、no-atom reference 和完整 metadata。

### Time-matched no-probe control

每个 \(q\) 必须有 control：

1. 使用相同 condensate preparation；
2. 保持相同总 hold time；
3. 发出相同 triggers 和 camera exposures；
4. 不让 probe light 到达 atoms；
5. 取得 terminal RAI。

control 用于扣除自然演化、trap loss、hold-time heating 和 preparation drift。没有 matched control 的 depletion curve 不用于定量 disturbance fitting。

### Acquisition ordering

- 随机化 \(q\) 的顺序；
- probed 与 control blocks 交错；
- run 开始、中间、结束重复 reference condensate 和 no-atom optical blocks；
- 若 background、power、detuning 或 preparation 离开 acceptance window，暂停并重新建立 baseline。

### 每个序列提取

从 Faraday frames：

- framewise signal；
- framewise \(\mathrm{SNR}_{3\times3}\)；
- profile width/centre；
- background drift；
- total H+V intensity；
- port imbalance。

从 terminal RAI：

- \(N_0(q)\)；
- \(N_{\mathrm{tot}}(q)\)；
- \(T(q)\)；
- condensate fraction；
- radii；
- profile residual。

### 结果定义

分别报告两个 boundary：

1. **image-quality boundary**：framewise signal 不再满足 T14 定义的 working requirement；
2. **condensate-state boundary**：condensate depletion/shape change 超过实验允许范围。

实际可用 sequence 截止于两者中先发生的一项。不要只用固定 30% loss，也不要只用 SNR threshold。

### Net-disturbance fitting

1. 先比较 probed 与 matched-control distributions。
2. 将 measured net change 拟合到 forward model 的 effective disturbance stage。
3. 可以更新 effective energy/depletion mapping，但不声称从同一曲线分别测出了：
   - one- versus two-recoil heating；
   - reabsorption fraction；
   - \(\eta_{\mathrm{coll}}\)。
4. 只有额外独立 observable 能分离这些机制时才增加参数。

### 输出

- \(N_0,T,\) fraction and radii versus \(q\) relative to controls；
- framewise SNR versus frame number；
- net-disturbance curve；
- image-quality and condensate-state boundaries；
- experimentally available sequence length；
- refined effective disturbance parameter(s)。

---

## T16. Frozen-parameter held-out comparison

### 目标

检验经过 calibration 的 forward model 是否能描述未用于调参的新数据。只有这一阶段能够支持 experimental validation language。

### 在 calibration 前预留

至少预先指定一种 held-out variation：

- 一个未用于 \(\kappa_F\) fitting 的 column-density range；
- 一个未用于 fluence fitting 的 pulse area；
- 若实验允许，一个邻近 detuning；
- 若干未用于 disturbance fitting 的 \(q\) values；
- 一个独立 acquisition day/block。

不要从同一 sequence 随机抽几帧作为 held-out data，因为它们共享 preparation、background 和 drift。

### Protocol

1. 完成 T7–T15 calibration；
2. 冻结：camera calibration、registration、background procedure、\(\kappa_F\)、spatial response 和 net-disturbance parameters；
3. 用冻结参数为 held-out conditions 生成 prediction；
4. 保存 prediction 和 config hash；
5. 之后采集 held-out data；
6. 按冻结 pipeline 处理，不重新拟合；
7. 比较 amplitude、profile、SNR、state evolution 和 sequence boundary；
8. 报告 residuals 和 uncertainty。

### 语言规则

- 用于确定参数的数据：**calibration**；
- 调参过程中查看的数据：**model–data comparison**；
- 参数冻结后、条件未参与 calibration 的数据：才可称 **held-out validation/assessment**；
- 若 held-out data 不足，dissertation 中写 assessment 或 validation plan，不夸大。

### 输出

- frozen model release/config；
- held-out residual plots；
- passed/failed aspects；
- 需要 refinement 的具体 model stage；
- final experimentally informed operating region。

---

# 2. 每次实验 run 的标准执行模板

## 2.1 Run 开始前

- [ ] Laser safety 和实验室 interlock 按本组 SOP 检查；
- [ ] optical configuration ID 与实际安装一致；
- [ ] HWP、WP、tube lens、camera positions 已记录；
- [ ] resonant/detuned source selection verified；
- [ ] frequency chain 和 AOM settings recorded；
- [ ] power meter/photodiode calibration ID recorded；
- [ ] camera configuration ID recorded；
- [ ] code commit/config hash recorded；
- [ ] raw-data storage location 已建立；
- [ ] 当天 dataset roles 已预先指定：pilot/calibration/control/held-out。

## 2.2 每个 acquisition block 的顺序

1. dark frames；
2. no-atom H/V references；
3. reference condensate/RAI checks；
4. atom acquisition blocks；
5. interleaved controls；
6. no-atom H/V references；
7. dark frames；
8. power/frequency/timing end checks。

## 2.3 Run 中止或暂停条件

- camera saturation/clipping；
- optical pulse 不完整落在 exposure window；
- pulse area drift 超过 run acceptance；
- H/V centroid 或 separation 离开预定 ROI；
- background \(S_0\) 显著变化；
- resonant/detuned beam overlap 漂移；
- condensate preparation 离开 T11 acceptance window；
- frequency reference 失锁或 uncertainty 不可追溯；
- raw metadata 缺失。

发生上述情况时，保留所有 raw data 并标注 excluded reason，不删除文件。

---

# 3. Raw-data 与 metadata 规范

## 3.1 建议目录结构

```text
YYYYMMDD_runNN/
  run_manifest.yaml
  raw/
    dark/
    no_atom_HV/
    faraday_HV/
    rai_HV_or_triplets/
    controls/
    diagnostics/
  calibration/
  processed/
  figures/
  notes/
```

raw 目录只读保留；processed output 不得覆盖 raw frames。

## 3.2 每个 shot 必须记录的 metadata

| 类别 | 字段 |
|---|---|
| Identification | run ID、shot ID、timestamp、operator |
| Dataset role | pilot、calibration、control、held-out、excluded |
| Probe | resonant/detuned、line reference、physical detuning、AOM frequencies/orders |
| Pulse | commanded duration、measured \(\tau_{\mathrm{eff}}\)、\(P_{\mathrm{atoms}}\)、pulse energy、pulse count \(q\) |
| Beam | beam-profile calibration ID、beam centre/radii |
| Camera | camera ID、ROI、exposure、gain mode、trigger mode、calibration ID |
| Polarisation | HWP angle、WP configuration、H/V sign convention、background calibration ID |
| Atomic state | preparation ID、bias-field setting、spin-preparation setting、hold time |
| Files | H path、V path、dark paths、no-atom reference paths、RAI paths |
| Analysis | code commit、config hash、processing version |
| Quality | accepted/excluded、reason、operator note |

## 3.3 每个 calibration parameter 的记录

- central value；
- statistical uncertainty；
- systematic uncertainty；
- units；
- number of realisations；
- date range；
- fitting method；
- raw-data path；
- code/config version；
- applicable range；
- whether provisional, calibrated or frozen。

---

# 4. 建议的首周工作顺序

这不是固定日程，而是最短依赖路径。

## Session 1 — 不开原子的大光路检查

- 完成 T1；
- 确认 camera active area 和 current image footprint；
- 填完 WP candidate geometry；
- 输出 as-installed layout。

## Session 2 — Wollaston 双端口

- 完成 T2 初测；
- 确定两个 ROI 是否可同时放入同一 camera；
- 选择 provisional HWP balance angle；
- 冻结 provisional mechanical positions。

## Session 3 — Resonant/detuned 与 RAI 决策

- 完成 T3 beam-overlap measurement；
- 测试 WP-on \(H+V\) readout；
- 决定 same-realisation 或 matched-ensemble calibration route。

## Session 4 — Frequency, power and timing

- 完成 T5 frequency-chain table；
- 完成 T6 photodiode timing；
- 得到首轮可用 pulse-area range。

## Session 5 — Camera 与 polarisation baseline

- 完成 T4 provisional background；
- 完成 T8 dark/read-noise/linearity；
- 完成 T9 no-atom common-mode data。

## Session 6 — Reference condensate 与首张 Faraday 图

- 完成 T11 initial reference block；
- 按 T12 从低 pulse area 逐步寻找 repeatable signal；
- 不在首日同时开展完整 fluence scan。

完成首图以后，再安排 T7 精确空间标定、T13 \(\kappa_F\)、T14 fluence scan 和 T15 repeated-exposure sequence。

---

# 5. 最终实验交付清单

- [ ] as-installed optical layout；
- [ ] measured dual-port geometry；
- [ ] confirmed RAI route；
- [ ] measured frequency chain and detuning；
- [ ] measured power, beam profile and pulse timing；
- [ ] camera calibration；
- [ ] measured magnification and spatial response；
- [ ] frozen port registration/background procedure；
- [ ] reference-condensate distribution；
- [ ] first repeatable dual-port Faraday image；
- [ ] effective \(\kappa_F\) with uncertainty and applicability range；
- [ ] experimental single-frame SNR–fluence relation；
- [ ] fixed-fluence \(P/\tau\) comparison；
- [ ] net repeated-exposure disturbance with matched controls；
- [ ] image-quality and condensate-state sequence boundaries；
- [ ] frozen-parameter held-out comparison；
- [ ] updated Chapter 5 figures from measured parameters；
- [ ] final experimentally informed operating region。
