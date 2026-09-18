# Week 1 · 从课程 IVP 到 PMSM 电流模型

目标：能够从电压关系推导出 `f(t,y)`，说明每个符号、初值和假设，独立求平衡点与 Jacobian，再用 Euler 和参考解检验。对应分工图片第一栏以及 Lecture 1 的“模型—方法—验证”路径。

## 1. 课程概念和本项目的对应

| 课程概念 | 课程例子 | PMSM 中的具体对象 |
|---|---|---|
| 状态向量与 IVP | Logistic 的 y；Van der Pol 的 (位置,速度) | y=(i_d,i_q)^T，给定 t=0 的两轴电流 |
| 向量场 f(t,y) | 生长率或振子的加速度 | 电压差除以电感得到电流变化率 |
| 适定性 | 局部 Lipschitz 与 Picard–Lindelöf | f=Ay+b 对 y 全局 Lipschitz；常系数下存在全局唯一解 |
| Euler | 左端点矩形积分 | 用当前两轴电流同时更新下一时刻的两轴电流 |
| 平衡点、Jacobian | Logistic 的 0、K；Van der Pol 的导数矩阵 | 求电气稳态电流与扰动模式 |
| 参考解 | Logistic 解析解；独立 SciPy 解 | 常系数矩阵指数解，另以 Radau 与零转速解析特例交叉验证 |
| 收敛与稳定 | 误差 log–log 图；y'=-y 的 h=2 边界 | 全轨迹电流误差；复特征值对 Euler 的更严格步长限制 |

## 2. 明确建模边界

研究永磁同步电机（PMSM）在转子同步旋转 d-q 坐标系下的**电气电流暂态**。取 d 轴沿永磁体磁链方向，q 轴采用 Park 变换的 `-sin` 行约定；三相平衡，零序省略。

前两周基准假设为：R、L_d、L_q、ψ_f 均为常数；电角速度 ω_e 由外部给定且恒定；u_d、u_q 在本算例区间内恒定；忽略饱和、温升、铁耗、PWM 开关纹波和机械反馈。数值模型是平均电气模型，不是完整电机启动或闭环控制模型。

**恒速不等于零转速。** 基准 ω_e=1000 rad/s。若给的是机械转速，必须先换算 ω_e=pω_m，其中 p 是极对数；若给 rpm，ω_m=2π·rpm/60。代码直接接收电角速度，不再乘一次 p。

| 符号 | 含义 | SI 单位 | 教学基准值 |
|---|---|---|---:|
| t | 时间 | s | 0≤t≤0.03 |
| i_d,i_q | 两轴电流，状态顺序固定为 d 后 q | A | 初值 (0,0) |
| u_d,u_q | 两轴施加电压 | V | (0,60) |
| R | 相电阻，沿用同一变换约定 | Ω | 0.5 |
| L_d,L_q | 两轴电感 | H | (0.003,0.004) |
| ψ_f | 永磁体磁链 | Wb | 0.05 |
| ω_e | 电角速度 | rad/s | 1000 |

这些数值用于可复现的教学验证，不来自实验或报告。项目若另有定稿参数，统一修改 `parameters.json` 并重新计算所有图表。

## 3. 从磁链与电压得到 ODE

电压在旋转坐标系中的形式为

$$
\begin{pmatrix}u_d\\u_q\end{pmatrix}
=R\begin{pmatrix}i_d\\i_q\end{pmatrix}
+\frac{d}{dt}\begin{pmatrix}\psi_d\\\psi_q\end{pmatrix}
+\omega_e\begin{pmatrix}0&-1\\1&0\end{pmatrix}
\begin{pmatrix}\psi_d\\\psi_q\end{pmatrix}.
$$

旋转项体现坐标轴本身转动带来的电压。取

$$\psi_d=L_di_d+\psi_f,\qquad \psi_q=L_qi_q,$$

且电感、永磁磁链为常数，得到

$$u_d=Ri_d+L_d\dot i_d-\omega_eL_qi_q,$$
$$u_q=Ri_q+L_q\dot i_q+\omega_e(L_di_d+\psi_f).$$

这一符号约定与 [MathWorks PMSM 电气方程](https://www.mathworks.com/help/sps/ref/pmsm.html) 一致。若另一资料使用相反的 q 轴方向，必须连同变换、电流和电压一起转换，不能只挑一个符号改动。

把导数移到左边：

$$
\dot i_d=\frac{u_d-Ri_d+\omega_eL_qi_q}{L_d},\qquad
\dot i_q=\frac{u_q-Ri_q-\omega_eL_di_d-\omega_e\psi_f}{L_q}.
$$

因此完整 IVP 是

$$
y'=f(t,y)=Ay+b,\quad y(0)=\begin{pmatrix}0\\0\end{pmatrix},\quad 0\le t\le0.03,
$$

$$
A=\begin{pmatrix}
-R/L_d&\omega_eL_q/L_d\\
-\omega_eL_d/L_q&-R/L_q
\end{pmatrix},\qquad
b=\begin{pmatrix}u_d/L_d\\(u_q-\omega_e\psi_f)/L_q\end{pmatrix}.
$$

量纲核查：V/H=A/s；ω_eL_qi_q 的单位是 V；ω_eψ_f 的单位也是 V。A 的元素单位为 1/s，b 的单位为 A/s。弧度在量纲上视为无量纲。

代码 `model.py:matrices` 的四个 A 元素与两个 b 元素应逐项核对，尤其是 A_12 为正、A_21 为负，以及 q 轴 b 中的负反电动势。

## 4. 为什么这个初值问题有唯一解

对任意 y、z，

$$\|f(t,y)-f(t,z)\|=\|A(y-z)\|\le\|A\|\|y-z\|.$$

因此 f 对 y 全局 Lipschitz。常系数线性系统不发生有限时刻爆破，解在任意有限时间区间唯一存在。Picard 迭代的积分映射在长度 Δ 满足 Δ‖A‖<1 的短区间上是压缩映射，再逐段延拓。

对初始扰动 δy_0，Gronwall 给出一般估计

$$\|\delta y(t)\|\le e^{\|A\|t}\|\delta y_0\|.$$

该上界保证连续依赖，但通常很保守；本模型实际有更强的衰减性质。不要把“上界里有指数增长”误读成真实解一定增长。

## 5. 平衡点：先算出电气稳态

令导数为零，而不是令两个电流都为零：

$$
\begin{pmatrix}R&-\omega_eL_q\\\omega_eL_d&R\end{pmatrix}
\begin{pmatrix}i_d^*\\i_q^*\end{pmatrix}
=\begin{pmatrix}u_d\\u_q-\omega_e\psi_f\end{pmatrix}.
$$

记 D=R²+ω_e²L_dL_q>0，则

$$i_d^*=\frac{Ru_d+\omega_eL_q(u_q-\omega_e\psi_f)}{D},$$
$$i_q^*=\frac{R(u_q-\omega_e\psi_f)-\omega_eL_du_d}{D}.$$

因此本基准 D=12.25，

$$i_d^*=40/12.25\approx3.26530612\ \mathrm A,\qquad
i_q^*=5/12.25\approx0.40816327\ \mathrm A.$$

检查方法：代回两条电压平衡关系，并在代码中检查 ‖f(0,y*)‖ 接近浮点误差。本模型中“初值是零”和“平衡点是零”是两件事。

## 6. Jacobian 与连续系统稳定性

逐项对电流求偏导：

$$J_f=\frac{\partial(f_d,f_q)}{\partial(i_d,i_q)}=A.$$

这里 J_f 在任何时间、任何状态都相同；它既是局部线性化，也是整个扰动方程的精确矩阵。电感不相等会让两个耦合系数不同，不代表 f 对电流非线性。

从二阶特征方程

$$\lambda^2+R\left(\frac1{L_d}+\frac1{L_q}\right)\lambda+
\left(\frac{R^2}{L_dL_q}+\omega_e^2\right)=0$$

得到

$$
\lambda_{\pm}=-\frac R2\left(\frac1{L_d}+\frac1{L_q}\right)
\pm\sqrt{\frac{R^2}{4}\left(\frac1{L_d}-\frac1{L_q}\right)^2-\omega_e^2}.
$$

R>0、L_d,L_q>0 时迹为负、行列式为正，因此两个特征值实部都为负。若根号内为负，扰动呈衰减振荡。基准约为 -145.8333±999.7830i，单位 1/s；衰减时间约 6.86 ms，振荡周期约 6.28 ms。

还可以作不依赖特征向量条件数的核查。记 x=y-y*，取

$$W(x)=\tfrac12(L_d^2x_d^2+L_q^2x_q^2).$$

沿 x'=Ax 有

$$\dot W=-RL_dx_d^2-RL_qx_q^2
\le-\frac{2R}{\max(L_d,L_q)}W.$$

两项交叉耦合恰好抵消，因此加权扰动范数衰减。这是数学 Lyapunov 函数，**不应叫作电机的总物理能量**。带电源、永磁体和耗散的本模型也不具备 Robertson 式“电流和守恒”。

## 7. 自建参考解，而不是借用报告数字

令 x=y-y*，则 x'=Ax，因此

$$y(t)=y^*+e^{At}(y_0-y^*).$$

可验证初值和导数都正确。`exact_solution` 用 SciPy 的 `expm` 计算矩阵指数，公式在数学上精确，计算机结果仍有浮点误差。

二维情形也能看出解析结构：设 a=tr(A)/2、B=A-aI、κ²=a²-det(A)，则 B²=κ²I，由幂级数分组

$$e^{At}=e^{at}\left[\cosh(\kappa t)I+
\frac{\sinh(\kappa t)}{\kappa}B\right].$$

κ=0 时取极限 I+tB；κ=iβ 时换成 cos(βt) 和 sin(βt)/β。它解释了衰减包络与旋转频率。

独立检查有两条：第一，Radau 在 rtol=10^-12、atol=10^-14 下和矩阵指数解比较；第二，ω_e=0 时两轴解耦，分别具有

$$i_j(t)=u_j/R+(i_j(0)-u_j/R)e^{-Rt/L_j},\quad j=d,q,$$

可直接检查，不依赖 PMSM 耦合矩阵的特征分解。

## 8. Euler 从积分出发

积分恒等式是

$$y(t_{n+1})=y(t_n)+\int_{t_n}^{t_n+h} f(s,y(s))\,ds.$$

把区间内斜率近似为左端点 f(t_n,y_n)，得到

$$y_{n+1}=y_n+h f(t_n,y_n).$$

PMSM 中对应两条式子：

$$i_{d,n+1}=i_{d,n}+\frac h{L_d}(u_d-Ri_{d,n}+\omega_eL_qi_{q,n}),$$
$$i_{q,n+1}=i_{q,n}+\frac h{L_q}(u_q-Ri_{q,n}-\omega_eL_di_{d,n}-\omega_e\psi_f).$$

必须同时使用旧的 i_d,n 与 i_q,n。先改 i_d 再拿新 i_d 更新 i_q 会变成另一种算法。NumPy 一次返回向量可以避免这种错误。

## 9. 为什么局部二阶、全局一阶

从精确初值开始的一步缺陷定义为

$$d_{n+1}=y(t_n+h)-y(t_n)-hf(t_n,y(t_n))
=\tfrac12h^2y''(t_n)+O(h^3).$$

因此 d=O(h²)，若使用归一化 LTE 则 τ=d/h=O(h)。令全局误差 e_n=y_n-y(t_n)，有

$$\|e_{n+1}\|\le(1+hL)\|e_n\|+Ch^2.$$

e_0=0 时求几何级数：

$$\|e_n\|\le Ch^2\sum_{j=0}^{n-1}(1+hL)^j
\le\frac C L(e^{LT}-1)h.$$

L=0 时用极限 CT h。结论依赖固定 T、足够光滑、传播受控和足够小的 h；不是任何大步长都保证一阶精度。PMSM 更可直接检查 d≈h²A²(y-y*)/2。

## 10. 第一周应能交出的证据

运行 `python code/run_all.py`，检查 Logistic 两步手算、标量 y'=-y 的边界试验、Van der Pol 指定个体练习，以及 PMSM 的电流和误差。PMSM 报告误差取

$$E(h)=\max_n\max_{j=d,q}|y_{n,j}-y_j(t_n)|,$$

同时保留末点误差。固定同一 T、初值、参数、单位和范数，只改变 h；观察阶

$$p_{\mathrm{obs}}=\log(E(h)/E(h/2))/\log 2.$$

本包的实际数值在 [运行结果](05_运行结果.md)。图像相似只是起点；解析公式、独立参考、收敛阶和稳定性预测共同构成可信证据。

依据：Lecture 1 第 8–13、18–19 页；Week1 Supplementary 第 15–40、48–51 页；Tutorial 1 第 4、9–10 页；Griffiths–Higham 第 2、6–7 章；PMSM 电气方程的来源见第 3 节。
