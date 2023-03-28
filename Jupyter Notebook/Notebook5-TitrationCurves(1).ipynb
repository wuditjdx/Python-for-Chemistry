{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "本文将介绍计算 51.3毫升 0.131mol/L 乙酸与 0.0953mol/L NaOH溶液的滴定曲线的方法。首先，导入库。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import scipy.optimize as opt"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "单质子酸和碱的滴定曲线分为三部分。对于弱酸与强碱的滴定来说，这三个部分分别是添加NaOH到半中和点（HA等于A$^{-}$），从半中和点到中和点（A$^{-}$大于HA），以及酸碱中和之后的部分。记录三个过程中所消耗的NaOH的量很重要。因此，在开始编程之前，让我们先输入达到中和点和半中和点时所添加的NaOH的体积吧！（中和点的时候弱酸的摩尔数*正好*等于强碱的摩尔数）"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "滴定前乙酸有 0.00672 mol。\n",
      "到达中和点时消耗了 0.0705 L 0.0953 mol/L NaOH。\n",
      "到达半中和点时消耗了 0.0353 L 0.0953 mol/L NaOH。\n"
     ]
    }
   ],
   "source": [
    "acid_VolI = 0.0513\n",
    "acid_ConI = 0.131\n",
    "\n",
    "acid_MolI = acid_VolI * acid_ConI\n",
    "print(\"滴定前乙酸有\", \"{:.3}\".format(acid_MolI), \"mol。\")\n",
    "\n",
    "base_Con = 0.0953\n",
    "eq_pt = acid_MolI/base_Con\n",
    "\n",
    "print(\"到达中和点时消耗了\", \"{:.3}\".format(eq_pt), \"L\", base_Con, \"mol/L NaOH。\")\n",
    "\n",
    "half_eq_pt = eq_pt/2\n",
    "print(\"到达半中和点时消耗了\", \"{:.3}\".format(half_eq_pt), \"L\", base_Con, \"mol/L NaOH。\")\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "现在我们得到了到达半中和点时所消耗的NaOH的量。接下来，我们需要生成一个数组，其中包含了一系列从初始到半中和点添加的NaOH的体积。同理，再生成一个数组，其中包含了一系列从半中和点到中和点添加的NaOH的体积。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "datapts = 100\n",
    "\n",
    "xrange = np.linspace(0,half_eq_pt,num=datapts,endpoint=False,retstep=True)\n",
    "PreHalfEqVols = xrange[0]\n",
    "VolStep = xrange[1]\n",
    "\n",
    "xrange = np.linspace(half_eq_pt,eq_pt,num=datapts,endpoint=False,retstep=True)\n",
    "PostHalfEqVols = xrange[0]\n",
    "VolStep2 = xrange[1]\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "现在我们有了两个数组，它将作为滴定曲线的x轴，我们需要对每个点进行平衡计算。从滴定开始到半中和点（直到加入的乙酸比乙酸盐多），下面是该平衡：\n",
    "\n",
    "$$\\textrm{CH}_{3}\\textrm{COOH}_{(aq)} + \\textrm{H}_{2}\\textrm{O} \\rightarrow \\textrm{H}_{3}\\textrm{O}^{+}_{(aq)} + \\textrm{CH}_{3}\\textrm{COO}^{-}_{(aq)} \\qquad \\textrm{K}_{a} = 1.76 \\textrm{x} 10^{-5}$$\n",
    "\n",
    "对于每个点的平衡计算，三段式中乙酸初始摩尔数将是乙酸的初始摩尔数*减去*添加的NaOH的摩尔数，乙酸盐的初始摩尔数将*等于*添加的NaOH的摩尔数，总体积将是酸溶液的初始体积*加上*添加的NaOH的体积。为了求解平衡，我们需要求解一个二次方程，所以让我们首先从 AgCl 溶解度notebook中复制二次方程求解器，并对其进行测试。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "该溶液的初始pH值为： 2.82\n"
     ]
    }
   ],
   "source": [
    "def quad_solve(a, b, c):\n",
    "    discrim = np.sqrt(b**2 - 4 * a * c)\n",
    "    sol1 = (-b + discrim)/(2 * a)\n",
    "    sol2 = (-b - discrim)/(2 * a)\n",
    "    return [sol1,sol2]\n",
    "\n",
    "solutions = quad_solve(1,1.76E-5,-0.131*1.76E-5)\n",
    "conc_H3O = solutions[0]\n",
    "pH = -np.log10(conc_H3O)\n",
    "\n",
    "print(\"该溶液的初始pH值为：\", \"{:.3}\".format(pH))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "好像运行成功了（当然，你可以自己算一下看看是否正确）。\n",
    "现在让我们计算这个过程中的pH值。\n",
    "注意浓度的计算需要考虑到总体积的变化！"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "equilibConst = 1.76E-5\n",
    "H3O_water = 1E-7\n",
    "init_baseConc = (PreHalfEqVols*base_Con)/(PreHalfEqVols+acid_VolI)\n",
    "init_acidConc = (acid_ConI*acid_VolI-PreHalfEqVols*base_Con)/(PreHalfEqVols+acid_VolI)\n",
    "\n",
    "pre_solutions = quad_solve(1,init_baseConc+equilibConst+H3O_water,init_baseConc*H3O_water-init_acidConc*equilibConst)\n",
    "preH3Oconc = pre_solutions[0]\n",
    "pre_pH = -np.log10(preH3Oconc)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "让我们检查一下计算出来的pH值是否正确。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.plot(PreHalfEqVols, pre_pH, 'ro')\n",
    "plt.rcParams['font.sans-serif'] = ['SimHei']\n",
    "plt.title(\"图1 - 用氢氧化钠滴定的乙酸\")\n",
    "plt.xlabel('添加的NaOH的体积')\n",
    "plt.ylabel('溶液的PH')\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {
    "collapsed": true
   },
   "source": [
    "图表具有我们预期的形状，并且pH的范围是我们期望的，看起来我们做对了！\n",
    "\n",
    "现在，按照上面的方法计算溶液从半中和点到当量点的pH值。注意！因为溶液中乙酸盐的含量比乙酸更高，所以主要反应是：\n",
    "\n",
    "$$\\textrm{CH}_{3}\\textrm{COO}^{-}_{(aq)} + \\textrm{H}_{2}\\textrm{O} \\rightarrow \\textrm{OH}^{-}_{(aq)} + \\textrm{CH}_{3}\\textrm{COOH}_{(aq)} \\qquad \\textrm{K}_{b} = \\frac{\\textrm{K}_{w}}{\\textrm{K}_{a}} = 5.67 \\textrm{x} 10^{-10}$$\n",
    "\n",
    "请注意，我已经创建了表示OH$^{-}$ 含量的列表（名为 PostHalfEqVols ）。同时请注意，随着乙酸的浓度变小，来自水电离的OH$^{-}$将变得很重要，因此您需要将其考虑在在计算中（我在上面并没有真正考虑它）。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "现在你已经计算了从半中和点到中和点溶液的pH值，请将所得的数据点与上面的数据点相结合，并用图标画出来，x轴为添加的NaOH的体积，纵轴为溶液的PH值。其中半中和点前画成红色，半中和点画成蓝色。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "最后，计算从中和点到添加0.15LNaOH 的pH值。这应该比前文更为直观，因为只需要考虑NaOH的OH$^{-}$。将得到的数据点与上面的数据点相结合，绘制整个滴定曲线，红色显示到半中和点的数据，蓝色显示从半中和点到中和点的数据，中和点之后的点用绿色来表示。"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
