## 如何区分混杂作用（Confounding Effect）和交互作用（Interaction)


### 写在前面：
 - 不是为了深入讨论，只是想尽量简单清晰地解释二者，以形成对两个概念大概的判断
 - 本文立足于流行病学领域，很可能其他领域有略微不同的理解和讨论
 - 文章里使用的“因果关系”是简化版的因果关系，为了方便理解，真实的因果关系推断很复杂
 - 这里不是讨论到底是不是因果。这里的假设就是，以下均认为存在暴露和结果之间的因果关系

#### 几个关键术语先解释一下：
```markdown

**偏移(Bias)、混杂(Confounding)、效应修饰(Effect modification)、交互作用(Interaction)**

1. 偏倚(Bias)：又叫做系统误差(Systematic Error)，基本上就是你的研究从设计上到操作上再到分析上有问题
  - 简单说就是，这种误差是你研究的问题，你得尽量想办法把这些误差减少，甚至理想状态下把它们消除掉

2. 混杂(Confounding)：自然存在的，它们不是你研究本身所造成的问题
  - 你不可能用任何方式把它们干掉，你只能想办法识别它们并控制它们，以限制它们对你的因果推断造成歪曲的影响

3. 效应修饰(Effect modification)：准确来说，流行病学所说的效应修饰是生物意义上的或者在因果理论意义上的第三方，影响我们的因果关系
  - 如果单纯放到统计学意义上，那我们会用交互作用(Interaction)来指代
  - 所以，很多时候，二者都被当成一个东西
  - 但一定要明白，交互作用必须是三个变量存在
```

#### 然后，我们来看一下怎么确定混杂和交互：
```markdown
# 混杂（Confounding）：判断一个第三方因素到底是不是混杂因素（Confounder），我们基本都知道有这几个判断标准
1. 混杂因素要和结果（Outcome；因变量）有因果关系
2. 混杂因素要和暴露（Exposure；自变量）有关系，但因果关系与否无所谓
3. 混杂因素绝对不能在暴露和结果的因果关系链条中间，也就是说，绝对不可能出现”暴露导致混杂因素进而导致结果“
- 满足这三个条件，那这个第三方因素基本就是混杂因素了，可是有时候我们用这三个条件去判断混杂与否并不容易
- 一个很简单很重要的特性是，混杂因素一定会歪曲你本来的因果关系。那我们到底如何知道，因果关系是不是被歪曲了呢？
- 我们可以用关联性测量来看(OR，RR，RD等)

**假设：暴露（自变量）是X，结果（因变量）是Y，第三方因素是Z。你原本的X和Y之间的Risk Ratio（RR）是3（Crude；Marginal)**

## 你可以用分层的方法来看RR有没有变化
 - 当然，前提是Z本身没有很多类别，最好就只有两个类别，有/无，男/女；如果类别很多，那你就要靠统计软件了

1. 第一种情况：当Z存在：RRz+ = 2；当Z不存在：RRz- = 2.1；

- 你发现分层后的两个RR在程度上基本一样，而且它们都和原本的RR不同
- RRz+和RRz-都和原本在RR（Crude）在一个方向（在RR>1的方向），但是RRz+和RRz-都更靠近1了
- 这种情况，我们说存在正向混杂作用（Positive Confounding Effect）
- 如果我们去算校正RR（Adjusted RR），就会发现这样的特点：校正RR比未校正RR更靠近1，而二者的方向都一样（要么都大于1，要么都小于1）

2. 第二种情况：当Z存在：RRz+ = 5；当Z不存在：RRz- = 4.8；

- 你发现分层后的两个RR在程度上基本一样，而且它们都和原本的RR不同
- RRz+和RRz-都和原本在RR（Crude）在一个方向（在RR>1的方向），但是RRz+和RRz-都更远离1了
- 这种情况，我们说存在负向混杂作用（Negative Confounding Effect）
- 如果我们去算校正RR（Adjusted RR），就会发现这样的特点：校正RR比未校正RR更远离1，而二者的方向都一样（要么都大于1，要么都小于1）

3. 第三种情况：当Z存在：RRz+ = 0.5；当Z不存在：RRz- = 0.7；

- 你发现分层后的两个RR在程度上基本一样，而且它们都和原本的RR不同
- RRz+和RRz-都和原本在RR（Crude）在不同方向
- 这种情况，我们说存在质性混杂作用（Qualitative Confounding Effect）
- 因为方向完全反了，性质完全变了
- 如果我们去算校正RR（Adjusted RR），就会发现这样的特点：校正RR和未校正RR在1的两边，方向不一样，一个大于1，另一个就要小于1

[Link](https://www.ctspedia.org/do/view/CTSpedia/BiasConfoundTypes/)

## 总结：到底怎么确定存在混杂作用呢？
- 首先，根据这个第三方因素分层以后的关联测量数（例如RRz+和RRz-），它们两个必须都在同一个方向
- 如果>1，两个都必须>1；如果<1，两个都必须<1；如果大于原本的Crude RR，那都要大于原本的Crude RR，绝对不能一个大一个小；要共进退
- 它们的程度要差不多(RRz+和RRz-差不多，不一定要相同，现实状况下，很可能都不同)
- 这样，我们就可以说，第三方因素和混杂因素，而且对我们的因果推断产生了混杂作用，因为它歪曲了我们的因果关系
- 它存在，所以我们得到了我们原本的Crude RR，可如果我们把它控制住看真实的RR就发现它和Crude RR不同


# 然后，我们来看交互作用。如果理解了上面的混杂，那交互就不难了。
- 交互作用有一个非常非常重要的前提，无论它的存在还是消失，你本来的因果关系都要是真实存在的因果关系
- 只是因为这个第三方因素的存在或消失，你的因果关系程度存在差异

**假设：暴露是X，结果是Y，第三方因素是Z。你原本的X和Y之间的RR是3(Crude;marginal)**

## 你可以用分层的方法来看RR有没有变化
  - 当然，前提是Z本身没有很多类别，最好就只有两个类别，有/无，男/女；如果类别很多，那你就要靠统计软件了
  - 如果是连续性变量，可以用回归分析等方法去看有没有交互作用）

1. 第一种情况：当Z存在：RRz+ = 5；当Z不存在：RRz- = 2；

- 你发现分层后的RR在程度上不同，RRz+ > Crude RR > RRz-
- 但是，三个RR的方向都一样，它们都>1，这样你可以说存在量化交互作用（Quantitative Interaction）
- 注意注意，Crude RR一定要在二者中间啊，如果都大于Crude RR或都小于Crude RR，基本就要往混杂作用上去看了

2. 第二种情况：当Z存在：RRz+ = 5；当Z不存在：RRz- = 1.01 或者 RRz- = 0.6；

- 你发现分层后的RR在程度上不同，RRz+ > Crude RR > RRz-
- 而且，有一个分层RR（RRz-）要么差不多到1了要么方向都变了，直接跨越1到另一边了
- 但是，另一个分层RR（RRz+）和Crude RR在同一边，而且要比Crude RR程度强
- 这种情况，我们说存在质性交互作用（Qualitative Interaction）
- 注意注意，Crude RR一定要在二者中间啊，如果都大于Crude RR或都小于Crude RR，基本就要往混杂作用上去看了

## 总结：
- 分层后的两个RR程度要不同，方向可以相同，也可以不同
- 但是必须有一个分层RR要比原来的RR（Crude RR）程度要强（远离1），而另一个要接近1甚至可以跨越1
- 简单理解就是，Crude RR要在两个分层RR之间，类似于你去用分层变量算出来一个加权后的RR，而这个加权后的RR也是真实的
- 只不过，分层以后，某一层的RR程度更强，另一层的RR程度变弱
- 这个第三方因素并没有扭曲你原本的因果关系

[Link](https://www.ctspedia.org/do/view/CTSpedia/InteractionGraph/)
```

