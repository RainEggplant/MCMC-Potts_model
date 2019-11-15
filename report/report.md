# 马氏链蒙特卡洛 Project 报告

（仅作草稿用）



## 一、二维高斯分布的估计





## 二、Potts 模型

### （一）数学推导

$$
p(x)=\frac{1}{Z(T)} \mathrm{exp}[-\frac{u(x)}{T}]
$$

$$
u(x)=-\sum_{i \leftrightarrow j} 1(x_i=x_j)
$$

$$
Z(T)=\sum_x \mathrm{exp}[-\frac{u(x)}{T}]
$$



> Zhiqiang Tan. 2015. Optimally adjusted mixture sampling and locally weighted
> histogram. In Technical Report, Department of Statistics, Rutgers University. 

统计上，Potts 模型的分布属于指数分布族，其充分统计量为 $-u(x)$, 自然参数 $\beta=T^{-1}$。令 $U=E\{u(x)\}$, 则有 $U=-\frac{\mathrm{d}}{\mathrm{d}\beta} \mathrm{ln}Z(T)$.

因为当 $\beta=0$ 时，由 (3) 式可得 $Z(\infty)=\sum_x 1=q^{K^2}=10^{400}$, 从而 $\mathrm{ln}Z(T)=-\int_0^{1/T} U \mathrm{d}\beta=400 \mathrm{ln}10-\frac{U}{T}$ .

