## 如何分析RDS数据
[English Version](https://github.com/Kin-Wang/EPI_Shelves/blob/main/RDS%20Analyst.pdf) and [R Code](https://github.com/Kin-Wang/EPI_Shelves/blob/main/RDS.R)

### 1. RDS是什么
- RDS，全名Respondent Driven Sampling（同伴推动抽样）。可以理解为“滚雪球抽样”的进阶版。适用于Hidden Population，也就是没有抽样框进行概率抽样的人群。例如，LGBTQ、性工作者、非法移民。
- “滚雪球抽样“是研究者通过一个一个的被访者介绍下一个被访者。这样收集到的数据很大可能是同质化的，也就是很可能这些人都差不多，甚至大部分都认识。这不是概率抽样，所以没有办法推断总体。
- RDS则是：研究者选择一定数量的种子，这些种子可以尽可能的多元化（例如，年龄，种族）。然后，这些种子去自行推荐熟人来参与调查，那些熟人又可以再推荐他们的熟人来参加。
- RDS的设计和操作并不是这篇文章的重点。如果想要了解如何使用RDS进行调查，那可以去自行搜索一下。

### 2. RDS数据集需要有哪些变量
- ID：每一个被调查者的独特编号
- Recruiter ID：每一个被调查者是被谁直接介绍来参加调查的；如果被调查者刚好是种子，那这个Recruiter ID就是“seed“；如果被调查者不是种子，那Recruiter ID就是介绍这位被调查者来参加调查的人的ID
- Network Size：如果被调查者可以无限介绍熟人来参加，理想状态下，他能介绍多少人来；也就是他到底认识多少人；当然必须要限制在你所关注的总体中（例如，男同性恋人群）
- Outcome Variables：你所关心的变量，例如年龄、种族、疾病状态等。
- 如果在你的数据集中没有Recruiter ID：那就必须有关于coupons分配的变量，包括：1）每一个被访者代表自己身份的那个独一无二的Coupon ID；2）每一个被访者收到的那些用来分配给自己的熟人们的Coupon IDs（例如，你每个人给五张coupon，让每个人最多只能招募五个熟人来参加调查，那你就要有五个变量来代表这些coupons，分别为Coupon.1 ID，Coupon.2 ID，直到Coupon.5 ID）

### 3. RDS诊断可能需要的一些问题
- 你认识多少人（需要具体，如果你研究性工作者，你就要问，你认识多少性工作者？）
- 你认识的那些人当中，有多少人住在XXX（地点）
- 最近两周，你认识的那些人当中，你见过多少个？
- 如果你需要多少coupon就给你多少，你觉得一天之内你最多能给多少人？
- 如果你需要多少coupon就给你多少，你觉得两周之内你最多能给多少人？
- 有多人你试图给他们coupon，但他们已经参加了这项调查？
- 除了那个给你coupon的人，你还知道有多少人参加了这项调查？
- 你给了几个人coupon？
- 有多少人拒绝了你的coupon？

### 4. 目前用来分析RDS数据的主要方法
- R语言当中的一个包（“RDS”）
- 基于R语言并结合Java编写的一个易于用户使用的互动性分析工具（RDS Analyst），使用起来的感觉有点像SPSS

### 5. 承载RDS数据的文件应该什么格式
- 推荐：.csv （使用Excel就可以将文件保存成这个格式）
- 其他：如果你使用R语言直接写代码分析数据，那么什么文件类型都是没有问题的，R语言里有各种匹配的包用来导入数据
- 其他：如果你使用RDS Analyst分析数据，那可接受的文件类型还包括：*.txt, *.sav, *.xpt, *.bdf, *.dta, *.sys, *.syd, *.arff, *.rec, *.mtp, *.s3

### 6. 关键概念解析
- sample fraction，样本比例，指的是样本量占总体的比例
- homophily，同质性，指的是具有某一特征的种子基本只招募同一特征的人，最后导致不同特征的人彼此隔绝，例如HIV感染者都只招募感染者，HIV非感染者都只招募非感染者；同质性越接近1程度越低。
- differential activity，差异性活动，指的是不同特征的人在网络数量（平均认识多少熟人）上的差别，例如HIV感染者认识的熟人整体很低，而HIV非感染者认识的熟人整体很高；越接近1差异程度越低。
- seed dependency，种子依赖性，种子的选择不随机，明显偏离真实群体情况，例如10个种子全部都是HIV感染者而没有一个非感染者，按照HIV流行率来看就很明显出现偏差。

### 7. 如何选择合适的estimator（estimator是用来利用样本情况估计总体情况的函数）
条件 |  RDS-I  |  RDS-II  |  Gile‘s SS  |  HCG
---  |---      |---      |---           |---
不知道总体数量 |  可以 |  可以 | 不可以 | 不可以
很大的样本比例 |  会偏倚 | 会偏倚 | 可以 | 可以
很短的波浪waves  |  会偏倚 | 会偏倚 | 可以 | 可以
不知道招募时间 |  可以 |  可以 | 可以 | 不可以
很高的同质性  |  有些偏倚 | 偏倚 | 偏倚 | 可以
种子依赖性   |  有些偏倚 | 偏倚 | 有些偏倚 | 可以
很高的差异性活动 |  不可以 |  不可以 | 可以 | 可以
连续性变量 |  不可以 |  可以 | 可以 | 可以

### 8. 几个estimators的关键性区别
- RDS-I和RDS-II全部基于放回抽样假设，而SS和HCG全部基于不放回抽样假设
  - 在实际调查中，放回抽样假设是不现实的，不会有研究者让被访者可以随机多次参与调查
- RDS-I基于马科夫链条假设，RDS-II利用Random Walk，而SS则是Successive Sampling加Configuration Graph Network Model，HCG基于Configuration Graph Network Model并且将同质性单独加入

### 9. 使用R语言的RDS包来分析数据
``` R
# 以一个格式为.csv的文件为例导入数据到R
dat1 <- read.csv("..path../name.csv")

# 检查你的数据
## 1. 重要的RDS变量有没有缺失值，以ID为例
anyNA(dat1$id)
## 2. 重要的RDS变量有没有重复
anyDuplicated(dat1$id)
## 3. 重要的RDS变量有没有超出正常范围
summary(dat1$id)
### 有时候，人为输入数据的时候，可能会把种子的coupon id输入成-1，需要将其改成0或者NA

# 如果没有Recruiter ID，需要利用Coupon IDs创建Recruiter ID
# 例如，你的数据中有这几个关于coupon的变量：own.coupon, coupon.1, coupon.2, coupon.3, coupon.4, coupon.5
dat1$recruiter.id <- rid.from.coupons(dat1,subject.coupon = "own.coupon",
                                     coupon.variables = paste0("coupon.",1:5),
                                     subject.id = "id")

# 将r.data.frame转换成rds.data.frame
# 数据中所需的几个变量：id，recruiter.id，network.size
install.packages("RDS")
library(RDS)
data_rds <- as.rds.data.frame(dat1, id = "id", recruiter.id = "recruiter.id", 
                              network.size = "network.size")
assert.valid.rds.data.frame(dat1)
```
``` R
# 使用图去了解样本的招募状况以及不同阶段的network size
# 以RDS包本身的数据为例
# 数据集为fauxmadrona：总体数量为1000，样本量为500，样本比例（sample fraction）为50%，同质性（Homophily）是5，差异性活动（Differential Activity）是1.8，不存在种子依赖性。
# disease是其中一个变量，是否得病，这样可以展示出每一个样本是否得病
data(fauxmadrona)
plot(fauxmadrona, plot.type = "Recruitment tree", stratify.by = "disease")
plot(fauxmadrona, plot.type = "Network size by wave", stratify.by = "disease")
```
1. 招募状况树（种子当中，患病者占比较低，符合真实状况）
![1](https://user-images.githubusercontent.com/60868837/103435055-9c769100-4bd7-11eb-94b0-29dd035ba7f9.png)
2. 网络数量图（患病人群和未患病人群的网络数量存在差异，也就是差异性活动比较大）
![2](https://user-images.githubusercontent.com/60868837/103435056-a13b4500-4bd7-11eb-8b4c-83b9310cac43.png)

``` R
# 使用estimators去估计总体，以fauxmadrona
# fauxmadrona当中没有“招募时间”，所以无法使用HCG
# RDS-I
RDS.I.estimates(fauxmadrona, outcome.variable = "disease")
# RDS-II
RDS.II.estimates(fauxmadrona, outcome.variable = "disease")
# 计算CI假设总体数量为1000; 对于单独的流行率结果来说，RDS-I和RDS-II都不依赖于总体数量
# SS
RDS.SS.estimates(fauxmadrona, outcome.variable = "disease", N=1000)
# 真实的总体疾病流行率应该是0.259; 
# RDS-I的结果是0.1592 (95% CI: 0.1396, 0.1788)
# RDS-II的结果是0.1641 (95% CI: 0.1367, 0.1916)
# Gile's SS的结果是0.1943 (95% CI: 0.1674, 0.2212)
```
``` R
# 诊断，以fauxmadrona为例
# 使用不同的estimators去诊断差异性活动
# RDS-I
differential.activity.estimates(fauxmadrona, outcome.variable = "disease", 
                                weight.type = "RDS-I")
# RDS-II
differential.activity.estimates(fauxmadrona, outcome.variable = "disease", 
                                weight.type = "RDS-II")
# Gile SS
differential.activity.estimates(fauxmadrona, outcome.variable = "disease", 
                                weight.type = "Gile's SS")
# 真正的差异性活动是1.8，也就是得病人平均网络数量是没得病人的1.8倍
# RDS-I的结果是1.66
# RDS-II的结果是1.77
# Gile's SS的结果是1.76

# 使用不同的estimators去诊断总体的同质性
# RDS-I
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "RDS-I")
# RDS-II
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "RDS-II")
# Gile's SS
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "Gile's SS")
# 真正的总体同质性是5
# RDS-I的结果是1.58
# RDS-II的结果是1.60
# Gile's SS的结果是1.66
# 三个都低于真正的总体同质性
```
``` R
# bottleneck plot
# 用来诊断同质性
# RDS-I
bottleneck.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.I.estimates)
# RDS-II
bottleneck.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.II.estimates)
# Gile's SS
bottleneck.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.SS.estimates)
```
3. RDS-I（没有明显的分支，说明同质性并没有很高）
![3](https://user-images.githubusercontent.com/60868837/103435057-a39d9f00-4bd7-11eb-80f9-8d551789ecf5.png)
4. RDS-II（没有明显的分支，说明同质性并没有很高）
![4](https://user-images.githubusercontent.com/60868837/103435058-a5fff900-4bd7-11eb-8377-6c3c6549f0e3.png)
5. SS（没有明显的分支，说明同质性并没有很高）
![5](https://user-images.githubusercontent.com/60868837/103435059-a8625300-4bd7-11eb-8455-4dca15ab4a49.png)

``` R
# convergence plot
# 用来诊断样本是否足够
# RDS-I
convergence.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.I.estimates)
# RDS-II
convergence.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.II.estimates)
# Gile's SS
convergence.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.SS.estimates)
```
6. RDS-I（随着样本量增加，出现聚合，说明样本量足够）
![6](https://user-images.githubusercontent.com/60868837/103435060-aa2c1680-4bd7-11eb-806c-458afccfcada.png)
7. RDS-II
![7](https://user-images.githubusercontent.com/60868837/103435062-adbf9d80-4bd7-11eb-8cd3-fdd67e3fd345.png)
8. SS
![8](https://user-images.githubusercontent.com/60868837/103435063-b021f780-4bd7-11eb-84de-6c8dd4c3d3c6.png)

``` R
# 使用fauxtime数据集来看看HCG的情况
# fauxtime有“招募时间”这个变量
# 如果你自己的数据中也有“招募时间”这个变量，别忘了在转换rds.data.frame的时候多加一个argument
data_time <- as.rds.data.frame(dat2, id = "id", recruiter.id = "recruiter.id", 
                              network.size = "network.size", time = "date")
# 估计总体
RDS.HCG.estimates(fauxtime, outcome.variable = "var1", N=1000)
# 诊断差异性活动
differential.activity.estimates(fauxtime, outcome.variable = "var1", 
                                weight.type = "HCG")
# 诊断同质性
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "HCG")
# bottleneck plot
bottleneck.plot(fauxtime, outcome.variable = "var1", est.func = RDS.HCG.estimates)
# convergence plot
convergence.plot(fauxtime, outcome.variable = "var1", est.func = RDS.HCG.estimates)
```
### 10. 如何用RDS Analyst分析数据
- 软件下载地址：http://wiki.stat.ucla.edu/hpmrg/index.php/RDS_Analyst_Install
- Windows和Mac平台都可以安装和运行

#### 1）将格式为.csv的数据导入到RDS Analyst当中
- 打开RDS Analyst，点击Open Data，选择要导入的数据
- 导入后，自动跳出以下界面（选择两种导入模式：一种需要有Recruiter ID，另一种只需要Coupon IDs）
- Coupon ID模式（请确保你的数据没有缺失值、非正常值等问题）

![Picture1](https://user-images.githubusercontent.com/60868837/103435543-190c6e00-4bde-11eb-831b-053ed6b58c4f.png)

- Recruiter ID模式（请确保你的数据没有缺失值、非正常值等问题）

![Picture2](https://user-images.githubusercontent.com/60868837/103435545-1b6ec800-4bde-11eb-8006-cf523bbadabf.png)

#### 2）可视化招募状况和网络数量
- 点击RDS Sample，再点击Diagnostic Plots
- 将想要分层的变量导入
- 选择Recruitment tree和Network size by wave

![Screen Shot 2021-01-01 at 02 55 37](https://user-images.githubusercontent.com/60868837/103435551-2c1f3e00-4bde-11eb-8ced-f7d72df05999.png)

#### 3) 估算总体并产生Bottleneck Plot和Convergence Plot
- 点击RDS Population，再点击Frequency Estimates
- 将想要分析的变量导入，例如疾病状态
- 选择合适的estimator
- 勾选bottleneck plot和convergence plot

![Screen Shot 2021-01-01 at 02 47 55](https://user-images.githubusercontent.com/60868837/103435554-30e3f200-4bde-11eb-9a78-d834eba01247.png)

#### 4) 诊断同质性和差异性活动
- 点击RDS Population，再点击Population Homophily
- 选择想要分析的变量导入，填写估计的总体数量
- 点击RDS Population，再点击Differential Activity
- 选择想要分析的变量导入，选择合适的estimator

#### 5）总结
- RDS Analyst还有很多功能，例如利用估算总体数量（后验），时间趋势数据的检验，等
- 上面只是最简单最基本的分析

### 11.参考文献
- Salganik MJ, Heckathorn DD. Sampling and estimation in hidden populations using respondent-driven sampling. Sociol Methodol. 2004;34(1):193-240.
- Volz E, Heckathorn DD. Probability based estimation theory for respondent driven sampling. J Off Stat. 2008;24(1):79-97.
- Gile KJ, Handcock MS. Respondent-driven sampling: an assessment of current methodology. Sociol Methodol. 2010;40(1):285-327.
- Gile KJ. Improved inference for respondent-driven sampling data with application to HIV prevalence estimation. J Am Stat Assoc. 2011;106(493):135-146.
- Gile KJ, Handcock MS. Network model-assisted inference from respondent-driven sampling data. J Royal Stat Soc Ser A (Stat Soc). 2015;178(3):619-639.
- Fellows IE. Respondent-Driven Sampling and the Homophily Configuration Graph Statistics in Medicine. 2019;38:131–150. 



