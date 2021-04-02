
# <center>时间序列中的平稳化方法</center>

## 1.不同的平稳化方法的使用顺序对平稳化序列的影响

首先将一个时间序列表示为

$$
{X_t} = {T_t} + {S_t} + {\eta _t}
$$


## 1.1 一般差分 + 季节差分

设经过一般$k$阶差分和$T$期季节差分后得到平稳序列

$$
{\nabla ^k}{\nabla _T}{X_t} = \left\{ {{\eta _t}} \right\}
$$


由于
$$
{\nabla ^k}{\nabla _T}{X_t} = {\left( {1 - \Im } \right)^k}\left( {1 - {\Im ^T}} \right){X_t} = \left( {1 - {\Im ^T}} \right){\left( {1 - \Im } \right)^k}{X_t} = {\nabla _T}{\nabla ^k}{X_t}
$$


故一般的差分与季节差分的使用顺序对平稳序列无影响.

## 1.2 一般差分 + 季节指数

先做$k$阶差分(注意，这里差分过后数据量将会减少$k$)

$$
\begin{array}{l}
{\nabla ^k}{X_t} = {S_t} + {\eta _t}\\
X_t^s = {\nabla ^k}{X_t}\\
{{\hat S}_t} = \bar X_t^s \times {I_j}\\
{\eta _t}{\rm{ = }}{\nabla ^k}{X_t}{\rm{ - }}\bar X_t^s \times {I_j}
\end{array}
$$


先做季节指数

$$
\begin{array}{l}
{{\hat S}_t} = {{\bar X}_t} \times {I_j}\\
{X_t} - {{\bar X}_t} \times {I_j} = {T_t} + {\eta _t}\\
{\eta _t} = {\nabla ^k}\left( {{X_t} - {{\bar X}_t} \times {I_j}} \right)
\end{array}
$$


故一般差分与季节指数的使用顺序对最终的平稳序列有影响，这是因为先做了差分会导致数据量减少，进而影响季节指数相关均值的计算.

## 1.3 函数回归 + 季节差分

### 1.3.1 先函数回归, 再季节差分

由最小二乘法,先对数据进行函数回归有

$$
\begin{array}{l}
\xi  = {\left( {{A^T}A} \right)^{ - 1}}{A^T}{X_t}\\
{{\hat T}_t} = A\xi  = A{\left( {{A^T}A} \right)^{ - 1}}{A^T}{X_t}\\
{X_t} = A\xi  + {S_t} + {\eta _t}\\
{S_t} - {\eta _t^1} = \left[ {E - A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}} \right]{X_t}
\end{array}
$$


其中
$$
A = \left[ {\begin{array}{*{20}{c}}
1&{{x_1}}& \cdots &{x_1^k}\\
1&{{x_2}}& \cdots &{x_2^k}\\
 \vdots & \vdots & \ddots & \vdots \\
1&{{x_n}}& \cdots &{x_n^k}
\end{array}} \right]
$$


再由季节差分，可得

$$
\begin{array}{l}
{\nabla _T} = \left( {1 - {\Im ^T}} \right)\\
X_t^R = \left[ {E - A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}} \right]{X_t}\\
{\nabla _T}X_t^R = {\eta _t^1}
\end{array}
$$


最终得到
$$
{\eta _t^1} = \left( {1 - {\Im ^T}} \right)\left[ {E - A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}} \right]{X_t}
$$


### 1.3.2 先季节差分,再函数回归

  先季节差分得
$$
{\nabla _T}{X_t} = \left( {1 - {\Im ^T}} \right){X_t} = {T_t} + {\eta _t^2}
$$


再函数回归有

$$
\begin{array}{l}
X_t^D = {\nabla _T}{X_t}\\
{A_d} = {A_{n - T}}
\end{array}
$$


*注意,此时因为先做了季节差分,所以数据少了一个周期T的量*

$$
{{\hat T}_t} = {A_d}\xi  = {A_d}{\left( {{A_d}^T{A_d}} \right)^{ - 1}}{A_d}^TX_t^D
$$


最终得到
$$
{\eta _t^2} = \left[ {E - {A_d}{{\left( {{A_d}^T{A_d}} \right)}^{ - 1}}{A_d}^T} \right]\left( {1 - {\Im ^T}} \right){X_t}
$$


故$$\eta _t^1 \ne \eta _t^2$$

## 1.4 函数回归 + 季节指数

### 1.4.1 先函数回归

$$
\begin{array}{l}
\xi  = {\left( {{A^T}A} \right)^{ - 1}}{A^T}{X_t}\\
{{\hat T}_t} = A\xi  = A{\left( {{A^T}A} \right)^{ - 1}}{A^T}{X_t}\\
{X_t} = A\xi  + {S_t} + {\eta _t}\\
{S_t} - {\eta _t^1} = \left[ {E - A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}} \right]{X_t}
\end{array}
$$



记$${Y_t} = \left[ {E - A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}} \right]{X_t}$$

则
$$
\eta _t^1 = {{\bar Y}_t} \times {I_j} - {Y_t}
$$


### 1.4.2 先季节指数

$$
\begin{array}{l}
{{\hat S}_t} = {{\bar X}_t} \times {I_j}\\
{X_t} - {{\bar X}_t} \times {I_j} = {T_t} + {\eta _t^2}
\end{array}
$$



记$${Y_t} = {X_t} - {{\bar X}_t} \times {I_j}$$

则
$$
\left\{ {\begin{array}{*{20}{l}}
{\xi  = {{\left( {{A^T}A} \right)}^{ - 1}}{A^T}{Y_t}}\\
{{{\hat T}_t} = A\xi  = A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}{Y_t}}\\
\begin{array}{l}
{Y_t} = {{\hat T}_t} + \eta _t^2\\
\eta _t^2 = {Y_t} - {{\hat T}_t}
\end{array}\\
{\eta _t^2 = \left( {{X_t} - {{\bar X}_t} \times {I_j}} \right) - \left[ {A{{\left( {{A^T}A} \right)}^{ - 1}}{A^T}\left( {{X_t} - {{\bar X}_t} \times {I_j}} \right)} \right]}
\end{array}} \right.
$$


# 2. 测试


```python
import numpy as np
import matplotlib.pyplot as plt

# typedef
plot = plt.plot
show = plt.show
subplot = plt.subplot
legend = plt.legend
title = plt.title

import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['SimHei']   # 设置简黑字体
mpl.rcParams['axes.unicode_minus'] = False     # 解决 '-' bug
```

## 2.1 生成数据


```python
%matplotlib notebook
x = np.linspace(0, 10, 100)

np.random.seed(1)
e = 10 + np.exp(x/3)
y = x*np.sin(np.linspace(0, 10*np.pi, 100))
z = np.random.randn(1, 100)

Y = (e + y +z)[0]

xx, yy = x, Y

plot(yy)
show()
```

<img src="/img/1.png" width="720">


## 2.2 几种方法


```python
%matplotlib notebook
# 季节差分[去除周期] <-- seq: 原始序列 T:周期间隔
def dertSeason(seq, T):
    n = seq.shape[0]
    dertS = np.zeros((1, n-T))
    for i in range(n-T):
        dertS[0, i] = seq[T+i] - seq[i]
    return dertS[0]


# 季节指数[去除周期] <-- seq: 原始序列 T:周期间隔
def indexSeason(seq, T):
    n = seq.shape[0]
    cow = int(n/T)
    # ary = seq.reshape(T, cow)
    ary = np.reshape(seq, (T, cow))
    _Xmean = np.mean(ary)
    _Xkmean = np.mean(ary, 1)
    _Ij = _Xkmean / _Xmean
    _Ij = _Ij.reshape(T, 1)
    multA = _Xmean * _Ij
    return (ary - multA).reshape(1, n)[0]


# 回归[去除趋势] <-- seq: 原始序列 func: 回归函数
# 针对本实例, 选用 指数函数 回归 <-- func = beta * exp(ct)
# 最小二乘法 <-- Least square method
# 为了避免对数无法对负数取对数, 改用二次函数: func = a0 + a1*x + a2*x^2
def LS(xx2, seq, func='exp'):
    n = seq.shape[0]
    # _Lseq = np.log(seq).reshape(n, 1)
    _Lseq = seq.reshape(n, 1)
    A = np.ones((n, 3))
    A[:, 1] = xx2
    A[:, 2] = xx2 ** 2
    # [beta[0], c[0]] = np.linalg.inv(A.T @ A) @ A.T @ _Lseq
    [a, b, c] = np.linalg.inv(A.T @ A) @ A.T @ _Lseq
    return [a, b, c]


def funcRegression(xx, seq, func='exp'):
    [a, b, c] = LS(xx, seq)
    # _Rseq = seq - beta * np.exp(c*xx)
    _Rseq = a + b*xx + c*xx**2
    return seq - _Rseq

plot(yy)
plot(yy-funcRegression(xx, yy))
title('二次回归')
show()

# 作K阶差分[去除趋势] <-- seq: 原始序列 K: 差分阶数
def dertK(seq, k):
    n = seq.shape[0]
    dertK = np.zeros((1, n-k))
    for i in range(n-k):
        dertK[0, i] = seq[i+k] - seq[i]
    return dertK[0]
```

<img src="/img/2.png" width="720">




## 2.3 比较不同顺序

### 2.3.1 一般差分 + 季节差分


```python
%matplotlib notebook
# 先季节差分, 再做一般k阶差分
seaSeq1 = dertSeason(yy, 20)
subplot(3, 2, 1)
plot(seaSeq1)
title('1.季节差分')

kSeq1 = dertK(seaSeq1, 1)
subplot(3, 2, 2)
plot(kSeq1)
title('2.一般差分')

# 先k阶差分, 再季节差分
kSeq2 = dertK(yy, 1)
subplot(3, 2, 3)
plot(kSeq2)
title('1.一般差分')

seaSeq2 = dertSeason(kSeq2, 20)
subplot(3, 2, 4)
plot(seaSeq2)
title('2.季节差分')

subplot(3, 1, 3)
plot(kSeq1, 'r--')
plot(seaSeq2, 'g--')
legend(['季节差分->一般差分', '一般差分->季节差分'])
show()
print('两条曲线的距离差:', np.mean(np.abs(kSeq1 - seaSeq2)))
```

<img src="/img/3.png" width="1080">


    两条曲线的距离差: 0.0


### 2.3.2 函数回归 + 季节差分


```python
%matplotlib notebook
RS1 = funcRegression(xx, yy)
subplot(3, 2, 1)
plot(RS1)
title('1.二次回归')

dertSeq1 = dertSeason(RS1, 20)
subplot(3, 2, 2)
plot(dertSeq1)
title('2.季节差分')

dertSeq2 = dertSeason(yy, 20)
subplot(3,2, 3)
plot(dertSeq2)
title('1.季节差分')

RS2 = funcRegression(xx[20:], dertSeq2)
subplot(3, 2, 4)
plot(RS2)
title('2.二次回归')

subplot(3, 1, 3)
plot(dertSeq1, 'r--')
plot(RS2, 'g--')
legend(['二次回归->季节差分', '季节差分->二次回归'])
show()
print('两条曲线的距离差:', np.mean(np.abs(dertSeq1 - RS2)))
```

<img src="/img/4.png" width="1080">


    两条曲线的距离差: 0.965382053089207


### 2.3.3 一般差分 + 季节指数


```python
%matplotlib notebook
dertSeq31 = dertK(yy, 1)
subplot(3, 2, 1)
plot(dertSeq31)
title('1.一般差分')

seaIdx31 = indexSeason(np.hstack((dertSeq31, np.array([0]))), 20)
subplot(3, 2, 2)
plot(seaIdx31)
title('2.季节指数')

seaIdx32 = indexSeason(yy, 20)
subplot(3, 2, 3)
plot(seaIdx32)
title('1.季节指数')

dertSeq32 = dertK(seaIdx32, 1)
subplot(3, 2, 4)
plot(dertSeq32)
title('2.一般差分')

subplot(3, 1, 3)
plot(seaIdx31, 'r--')
plot(dertSeq32, 'g--')
legend(['一般差分->季节指数', '季节指数->一般差分'])
show()
print('两条曲线的距离差:', np.mean(np.abs(seaIdx31[:-1] - dertSeq32)))
```

<img src="/img/5.png" width="1080">


    两条曲线的距离差: 1.537409833561132


### 2.3.4 函数回归 + 季节指数


```python
%matplotlib notebook
RS41 = funcRegression(xx, yy)
subplot(3, 2, 1)
plot(RS41)
title('1.二次回归')

seaIdx41 = indexSeason(RS41, 20)
subplot(3, 2, 2)
plot(seaIdx41)
title('2.季节指数')

seaIdx42 = indexSeason(yy, 20)
subplot(3, 2, 3)
plot(seaIdx42)
title('2.季节指数')

RS42 = funcRegression(xx, seaIdx42)
subplot(3, 2, 4)
plot(RS42)
title('2.二次回归')

subplot(3, 1, 3)
plot(seaIdx41, 'r--')
plot(RS42, 'g--')
legend(['二次回归->季节差分', '季节差分->二次回归'])
show()
print('两条曲线的距离差:', np.mean(np.abs(seaIdx41 - RS42)))
```

<img src="/img/6.png" width="1080">


    两条曲线的距离差: 0.27429634614795284

