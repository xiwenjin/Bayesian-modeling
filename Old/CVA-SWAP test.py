{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "from  QuantLib import *\n",
    "import numpy as np\n",
    "from math import *\n",
    "from pylab import *\n",
    "\n",
    "def A(t,T):\n",
    "    forward = crvToday.forwardRate(t, t,Continuous, NoFrequency).rate()\n",
    "    value = B(t,T)*forward - 0.25*sigma*B(t,T)*sigma*B(t,T)*B(0.0,2.0*t);\n",
    "    return exp(value)*crvToday.discount(T)/crvToday.discount(t);\n",
    "\n",
    "def B(t,T):\n",
    "    return (1.0-exp(-a*(T-t)))/a;\n",
    "\n",
    "def gamma(t):\n",
    "        forwardRate =crvToday.forwardRate(t, t, Continuous, NoFrequency).rate()\n",
    "        temp = sigma*(1.0 - exp(-a*t))/a\n",
    "        return (forwardRate + 0.5*temp*temp)\n",
    "\n",
    "def gamma_v(t):\n",
    "    res=np.zeros(len(t))\n",
    "    for i in xrange(len(t)):\n",
    "        res[i]=gamma(t[i])\n",
    "    return res"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "Nsim=100\n",
    "a=0.1\n",
    "sigma=0.02\n",
    "\n",
    "todaysDate=Date(30,12,2013);\n",
    "Settings.instance().evaluationDate=todaysDate;\n",
    "crvToday=FlatForward(todaysDate,0.02,Actual360())\n",
    "\n",
    "r0=forwardRate =crvToday.forwardRate(0,0, Continuous, NoFrequency).rate()\n",
    "months=range(1,12*5,1)\n",
    "sPeriods=[str(month)+\"m\" for month in months]\n",
    "Dates=[todaysDate]+[todaysDate+Period(s) for s in sPeriods]\n",
    "T=[0]+[Actual360().yearFraction(todaysDate,Dates[i]) for i in xrange(1,len(Dates))]\n",
    "T=np.array(T)\n",
    "rmean=r0*np.exp(-a*T)+ gamma_v(T) -gamma(0)*np.exp(-a*T)\n",
    "\n",
    "np.random.seed(1)\n",
    "stdnorm = np.random.standard_normal(size=(Nsim,len(T)-1))\n",
    "\n",
    "rmat=np.zeros(shape=(Nsim,len(T)))\n",
    "rmat[:,0]=r0\n",
    "\n",
    "for iSim in xrange(Nsim):\n",
    "    for iT in xrange(1,len(T)):\n",
    "        mean=rmat[iSim,iT-1]*exp(-a*(T[iT]-T[iT-1]))+gamma(T[iT])-gamma(T[iT-1])*exp(-a*(T[iT]-T[iT-1]))\n",
    "        var=0.5*sigma*sigma/a*(1-exp(-2*a*(T[iT]-T[iT-1])))\n",
    "        rmat[iSim,iT]=mean+stdnorm[iSim,iT-1]*sqrt(var)\n",
    "\n",
    "startDate=Date(30,12,2013);\n",
    "\n",
    "crvMat= [ [ 0 for i in xrange(len(T)) ] for iSim in range(Nsim) ]\n",
    "npvMat= [ [ 0 for i in xrange(len(T)) ] for iSim in range(Nsim) ]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "for row in crvMat:\n",
    "    row[0]=crvToday\n",
    "\n",
    "for iT in xrange(1,len(T)):\n",
    "    for iSim in xrange(Nsim):\n",
    "        crvDate=Dates[iT];\n",
    "        crvDates=[crvDate]+[crvDate+Period(k,Years) for k in xrange(1,21)]\n",
    "        crvDiscounts=[1.0]+[A(T[iT],T[iT]+k)*exp(-B(T[iT],T[iT]+k)*rmat[iSim,iT]) for k in xrange(1,21)]\n",
    "        crvMat[iSim][iT]=DiscountCurve(crvDates,crvDiscounts,Actual360(),TARGET())\n",
    "\n",
    "#indexes definitions\n",
    "forecastTermStructure = RelinkableYieldTermStructureHandle()\n",
    "index = Euribor(Period(\"6m\"),forecastTermStructure)\n",
    "\n",
    "#swap 1 definition\n",
    "maturity = Date(30,12,2018);\n",
    "fixedSchedule = Schedule(startDate, maturity,Period(\"6m\"), TARGET(),ModifiedFollowing,ModifiedFollowing,DateGeneration.Forward, False)\n",
    "floatingSchedule = Schedule(startDate, maturity,Period(\"6m\"),TARGET() ,ModifiedFollowing,ModifiedFollowing,DateGeneration.Forward, False)\n",
    "swap1 = VanillaSwap(VanillaSwap.Receiver, 1000000,fixedSchedule,0.05 , Actual360(),floatingSchedule, index, 0,Actual360())  #0.01215\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<QuantLib.QuantLib.Euribor; proxy of <Swig Object of type 'EuriborPtr *' at 0x107366bd0> >"
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "index"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "for iT in xrange(len(T)):\n",
    "    Settings.instance().evaluationDate=Dates[iT]\n",
    "    allDates= list(floatingSchedule)\n",
    "    fixingdates=[index.fixingDate(floatingSchedule[iDate]) for iDate in xrange(len(allDates)) if index.fixingDate(floatingSchedule[iDate])<=Dates[iT]]\n",
    "    if fixingdates:\n",
    "        for date in fixingdates[:-1]:\n",
    "            try:index.addFixing(date,0.0)\n",
    "            except:pass\n",
    "        try:index.addFixing(fixingdates[-1],rmean[iT])\n",
    "        except:pass\n",
    "    discountTermStructure = RelinkableYieldTermStructureHandle()\n",
    "    swapEngine = DiscountingSwapEngine(discountTermStructure)\n",
    "    swap1.setPricingEngine(swapEngine)\n",
    "\n",
    "    for iSim in xrange(Nsim):\n",
    "        crv=crvMat[iSim][iT]\n",
    "        discountTermStructure.linkTo(crv)\n",
    "        forecastTermStructure.linkTo(crv)\n",
    "        npvMat[iSim][iT]=swap1.NPV()\n",
    "\n",
    "npvMat=np.array(npvMat)\n",
    "npvMat[npvMat<0]=0\n",
    "EE=np.mean(npvMat,axis=0)\n",
    "\n",
    "S=0.05 #constant CDS spread\n",
    "R=0.4  #Recovery rate 40%\n",
    "\n",
    "sum=0"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "CVA= 17008.5902689\n",
      "\n",
      "EE\n",
      "[[ 143606.52236152  110892.25913033  123290.20415129 ...,   13967.58057333\n",
      "    13977.69516598   14005.16151835]\n",
      " [ 143606.52236152  126227.03001656  141763.93414938 ...,   13899.14723321\n",
      "    13944.45211361   14001.20809908]\n",
      " [ 143606.52236152  121295.88270064  113094.03325023 ...,   13806.22098036\n",
      "    13884.49671652   13951.8993805 ]\n",
      " ..., \n",
      " [ 143606.52236152  122839.24885826  100364.02025741 ...,   13945.0936529\n",
      "    13977.09778729   14012.67541012]\n",
      " [ 143606.52236152  132998.08151752  166252.5282214  ...,   14097.25299989\n",
      "    14095.04320873   14067.92423857]\n",
      " [ 143606.52236152  177737.59151528  194720.39115905 ...,   14071.24723422\n",
      "    14064.61259999   14040.44191233]]\n"
     ]
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAY0AAAD8CAYAAACLrvgBAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADl0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uIDIuMi4yLCBodHRwOi8vbWF0cGxvdGxpYi5vcmcvhp/UCwAAIABJREFUeJzt3Xt0nPV95/H3d0Z3WXfJV0mWfDd3jCzb0ISEq0mTQNN2Q9oGNyFLN4EmLdkksN1dctmkSbMnyZKmtARoTEpDaC7FJ4E4jiEhJDZYNtjG2MayjW35KluSJduSdZnv/jGPYDCyNZJmNNbo8zpnzsz8nt/zzPc5YH1mfr/nYu6OiIhIPEKpLkBERMYOhYaIiMRNoSEiInFTaIiISNwUGiIiEjeFhoiIxE2hISIicVNoiIhI3BQaIiISt4xUF5Bo5eXlXlNTk+oyRETGlPXr1x9194rB+qVdaNTU1NDQ0JDqMkRExhQz2xNPPw1PiYhI3BQaIiISN4WGiIjETaEhIiJxU2iIiEjcFBoiIhI3hYaIiMQt7c7TOB9190bYc+wkO5tPsLflFO+5eAqVJXmpLktEZMgUGiPk7vxHQxP//JudOJCdESI7M0x2RojMsHGgrYu9Lafoi7x5L/bXj53iK390ceqKFhEZJoXGCJw83cv/+s9X+MlL+7msqpiq0jxO9/RxujfC6d4+Orv7uGBKIe+7ZAozJ05gZsUE/mHldtbuPJbq0kVEhkWhMUzbD3XwicfWs+voSf72ujncdc0swiEbdL13zCrny09t5XB7F5MKc0ahUhGRxNFE+BC5Oz9ct5ebv/M87V29PPaxRXzqutlxBQbAkpllAKzdpV8bIjL2KDSGwN35ylNb+dyPN3PF9BKe+uQ7uHJm+ZC2MX9KIYU5GazREJWIjEEanopTf2B897e7uW3JdO5734Vx/7qIFQ4Z9bVlrNEvDREZg/RLIw7uzt8/vY3v/nY3y5ZM5wvvH15g9Fsys4w9x05xoK0zgVWKiCSfQmMQ/YHx4HO7uG3JdD7//gsxG35gACyZoXkNERmbFBrnEIk4X40JjC8kIDAA5k0uoDgvU/MaIjLmaE7jDJ3dffyu8Sirtx1m9dYjHOk4zYcXJy4wAEIhY1FtqeY1RGTMUWgEVm45xH807OP5xqN09USYkJ3BO+eUc+OFk3n/pVMTFhj9lswoY+WWw+xrOUVVqS4pIiJjg0IjsHFfG1sPdnDrwmqunT+RRbVlZGUkb/Ruccz5GgoNERkrFBqBT103m8/cODfhvyjOZs7EAkrzs1iz6xh/Wlc1Kp8pIjJSg36VNrNHzOyImb0ywLL/bmZuZuXBezOz+82s0cw2mdmCmL7LzGxH8FgW036FmW0O1rnfgr/aZlZqZquC/qvMrCQxuzyw7IzwqAUGROc1Fs8o5YVdLbj74CuIiJwH4hl/+R6w9MxGM6sCrgf2xjTfBMwOHncADwR9S4H7gEVAPXBfTAg8EPTtX6//s+4BVrv7bGB18D6tLJ5Rxv62Tva16HwNERkbBg0Nd38OaBlg0TeBzwKxX5NvBh71qLVAsZlNAW4EVrl7i7u3AquApcGyQndf49Gv248Ct8Rsa3nwenlMe9roP19jza6jKa5ERCQ+w5rpNbP3A/vdfeMZi6YB+2LeNwVt52pvGqAdYJK7HwQInieeo547zKzBzBqam5uHsUepMWviBMonZOl8DREZM4YcGmaWB/wd8L8HWjxAmw+jfUjc/UF3r3P3uoqKiqGunjJmxqIZZazVvIaIjBHD+aUxE6gFNprZ60AlsMHMJhP9pRB7KFAlcGCQ9soB2gEOB8NXBM9HhlHreW/JjDIOtXfx+rFTqS5FRGRQQw4Nd9/s7hPdvcbda4j+4V/g7oeAFcBtwVFUi4HjwdDSSuAGMysJJsBvAFYGyzrMbHFw1NRtwJPBR60A+o+yWhbTnlb676+hISoRGQviOeT2B8AaYK6ZNZnZ7efo/hSwC2gEvgt8AsDdW4AvAeuCxxeDNoCPAw8F6+wEng7avwpcb2Y7iB6l9dWh7drYMKM8n4qCbJ5+5SDdvZFUlyMick6WbmPpdXV13tDQkOoyhuQ7zzby9ZXbuWhaId/64OXMmjgh1SWJyDhjZuvdvW6wfrrK7XngznfP4l8+fAVNrZ2899u/5d9f2KuJcRE5L+kyIueJGy+czGVVxdz9xMv8j59u5jevHeHvP3AJpflZcW9jz7GT/Oa1Zp57rZkXdrVQmJvJnEkTmDOpgNmTCpg7qYB5UwrIDOu7gogMj4anzjORiPPw87v5h5Xb6I04tWX5zJ1cwLzJhcydXEBVaS7tnb20nOym5VQ3LSe6OdTeye93HmNPcARWdWkeV80q4+TpPl473MGu5pN090XnSz5w+TS+8cHLUrmLInIeind4Sr80zjOhkPFf3zmDd8wp5+nNh9h2qJ2tB9v5xZZDnC3fi3IzuWJ6CR+9qpar51RQU57/luW9fRFeP3aKL/3sVZ5vPIq7j+p1tkQkfSg0zlPzJhcyb3LhG+9Pdffy2uETHGjrpDg3k5L8LMrysyjJzxp0uCkjHGLWxAlcO38iv3mtmf1tnVSW6HLsIjJ0Co0xIi8rg8uqirmsqnjY21hQHb1G5Ia9bQoNERkWzYiOI/MmF5CbGWbDntZUlyIiY5RCYxzJCIe4pLKIDXsVGiIyPAqNcWbB9BJePdBOV09fqksRkTFIoTHOLKguoTfibGo6nupSRGQMUmiMM5dXRyfSNUQlIsOh0BhnyidkM70sT5PhIjIsCo1x6IrqEjbsbdP1rURkyBQa49Dl00s4euI0Ta2dqS5FRMYYhcY4tEDzGiIyTAqNcWjupALysnSSn4gMnUJjHMoIh7i0spgNe9tSXYqIjDEKjXFqwfRith5sp7N7dE7y23vsFH0RTbyLjHUKjXHqzZP8kv9r47EX9vDOrz/LzzYdSPpniUhyDRoaZvaImR0xs1di2r5uZtvMbJOZ/dTMimOW3WtmjWa23cxujGlfGrQ1mtk9Me21ZvaCme0wsx+aWVbQnh28bwyW1yRqpwUuj7nibTI9+fJ+/ud/Rv/XWfd6S1I/S0SSL55fGt8Dlp7Rtgq4yN0vAV4D7gUwswuAW4ELg3X+yczCZhYGvgPcBFwAfCjoC/A14JvuPhtoBW4P2m8HWt19FvDNoJ8kSGl+FrXl+Uk9guqZbYf59BMbWVhTyhXTS9i4T5cuERnrBg0Nd38OaDmj7Zfu3hu8XQtUBq9vBh5399PuvhtoBOqDR6O773L3buBx4GaL3j7uGuBHwfrLgVtitrU8eP0j4FrT7eYS6vLqYl7a25qUk/zW7jrGx/9tA/OnFPLwsjoW1pSy7ZAulCgy1iViTuOjwNPB62nAvphlTUHb2drLgLaYAOpvf8u2guXHg/6SIAuqSzh6opt9LYk9yW9TUxsfW95AVWkeyz9aT0FOJpdVFdHT52w92J7QzxKR0TWiO/eZ2d8BvcBj/U0DdHMGDic/R/9zbWugOu4A7gCorq4+R8USq/9Ofuv3tlBdNrQ7+Z083cvvGo+yqek4HV09dHT10nG6l46uHrYcaKcoN5Pv315PaX4WAJcGdxzcuK/tjfkUERl7hh0aZrYMeC9wrb85vtEEVMV0qwT6D5kZqP0oUGxmGcGvidj+/dtqMrMMoIgzhsn6ufuDwIMAdXV1Oq4zTnMnF5CfFWbVq4eZXpZPQXYGBTmZTMjJIC8zzJmDgftaOnlm22Ge2d7M2p3H6O6LEDIoyMmkICeDCdkZFOZk8o7Z5Xz2xnlMKcp9Y93JhTlUFGTrkuwiY9ywQsPMlgKfA65291Mxi1YA/25m3wCmArOBF4n+aphtZrXAfqKT5X/m7m5mzwJ/QnSeYxnwZMy2lgFrguXPuK6wl1DhkFFfW8pTmw/x1OZDca83ozyf25ZM55p5E6mrKSUrY/BRTjPj0spiXh6FQ3xFJHkGDQ0z+wHwLqDczJqA+4geLZUNrArmpte6+39z9y1m9gTwKtFhqzvdvS/Yzl3ASiAMPOLuW4KP+BzwuJn9H+Al4OGg/WHg+2bWSPQXxq0J2F85w/0fupzXDnfQ3tXLia5eTgRDTKcGOOmvJC+Lq+dUUFOeP6zPuqyqiF9tPczxzh6KcjNHWrqIpICl25f3uro6b2hoSHUZMoDf7mjmww+/yGMfW8RVs8pTXY6IxDCz9e5eN1g/nREuo+aSadHJ8Jf3aYhKZKxSaMioKcrLpLY8n40KDZExS6Eho+rSyiIdQSUyhik0ZFRdUlnMofYuDh3vSnUpIjIMCg0ZVW+c5KdDb0XGJIWGjKoLpxaSEbJRuSS7iCSeQkNGVU5mmLmTC3TFW5ExSqEho+7SqmI2NrUR0Z38RMYchYaMussqi+no6uX1YydTXYqIDJFCQ0bdJVVFgCbDRcYihYaMutkTC8jLCmteQ2QMUmjIqAuHjIumFemXhsgYpNCQlLi0sogtB9rp7o2kuhQRGYIR3blPZLgurSqm+7e72X6og4sri4a8/pYDx3n4t7vp7Okj4k7EwR0yw8bNl03jhgsmEQrplvIiiabQkJS4tDJ6ZvhzO5qpKs2lMCczrj/ykYjz0PO7+PrK7eRmhplclINhmEHIjJaT3Tz9yiHmTS7grmtmcdNFUwgrPEQSRvfTkJRwd+q/sprmjtMAmEFRbiYleVlcMKWQ9182lXfNrSA7I/zGOgfaOvn0ExtZs+sYN144ib//wCVv3IO8X29fhJ9tOsi3n9nBzuaTzKzI56+vmc37Lp2q8BA5h3jvp6HQkJR57XAHm5uO09bZQ9upblpPddN6soe1u45x7GQ3BTkZvOeiKdx82VSOnuzmf/50M70R5/Pvu5A/ravEzryJeYy+iPP0Kwf59upGth/u4DM3zuXOd88axb0TGVsUGjJm9fZF+N3OYzz58n5WvnKIk8GtZy+vLuZbH7yM6WXx3242EnE+8MDvCRn85BNXJatkkTEv3tDQnIacdzLCIa6eU8HVcyrovKWP1dsOc/J0L3+8oJKM8NAO+AuFjMUzynj4+V10dveRmxUefCUROSsdcivntdysMO+9ZCofXFg95MDoV19bQk+f89K+1gRXJzL+DPqv0MweMbMjZvZKTFupma0ysx3Bc0nQbmZ2v5k1mtkmM1sQs86yoP8OM1sW036FmW0O1rnfgoHqs32GyFBdMb0UM1i3W6EhMlLxfHX7HrD0jLZ7gNXuPhtYHbwHuAmYHTzuAB6AaAAA9wGLgHrgvpgQeCDo27/e0kE+Q2RIinIzmTe5kHWvt6S6FJExb9DQcPfngDP/td0MLA9eLwduiWl/1KPWAsVmNgW4EVjl7i3u3gqsApYGywrdfY1HZ+QfPWNbA32GyJDV15Swfk8rPX06A11kJIY7pzHJ3Q8CBM8Tg/ZpwL6Yfk1B27namwZoP9dniAxZfW0ZnT19bDnQnupSRMa0RE+ED3TgvA+jfWgfanaHmTWYWUNzc/NQV5dxYGFtdDR03W4NUYmMxHBD43AwtETwfCRobwKqYvpVAgcGaa8coP1cn/E27v6gu9e5e11FRcUwd0nS2cSCHGrK8nhBoSEyIsMNjRVA/xFQy4AnY9pvC46iWgwcD4aWVgI3mFlJMAF+A7AyWNZhZouDo6ZuO2NbA32GyLAsrCmlYU+LbjMrMgLxHHL7A2ANMNfMmszsduCrwPVmtgO4PngP8BSwC2gEvgt8AsDdW4AvAeuCxxeDNoCPAw8F6+wEng7az/YZIsNSX1tK26keGptPpLoUkTFr0DPC3f1DZ1l07QB9HbjzLNt5BHhkgPYG4KIB2o8N9Bkiw1VfWwrAC7tbmDOpIMXViIxNOiNcxo3q0jwmFmRrMlxkBBQaMm6YGQtrS3lxdwvpdqFOkdGi0JBxZVFtKYfau2hq7Ux1KSJjkkJDxpWFNdF5jRc1RCUyLAoNGVfmTiqgMCdD16ESGSaFhowroZCxsKZUvzREhkmhIePOwtpSdh09+cb9yUUkfgoNGXf65zUaNEQlMmQKDRl3Lp5WRE5miO+v3cOxE/q1ITIUCg0Zd7IyQnzmxnm8uLuF677xG368vinh523sPXaKf/nNTj7yry+y/VBHQrctkkqDXkZEJB3d/ge1/MGscu79ySY+/R8b+clLTXz5loupKc9/o09PX4SDbV0AVJflDbrNxiMdPL35EE+/cohXD7553445kwu496b5id8JkRSwdDsztq6uzhsaGlJdhowRkYjz2It7+drT2+jpi3DDhZNp7uhiX0snh9q76AuuiHvh1EL+6PJpvO/SqUwqzAHA3dlx5AQ/33SQpzYfZMeR6IUQF1QXc9NFU1h60WT+5ocvE3Hnp5+4KmX7KBIPM1vv7nWD9lNoiMCh41186Wevsn5PK9NKcqkqyaWqNI/Kklw6unpZsfEAm5qOYwZXzizjoqlF/GrrYXY2n8QM6mtKuemiySy9aAqTi3Le2O7XfrGN7z63i02fv4G8LP2wl/NXvKGh/4tFgMlFOXznzxecdfnH3jGDnc0nePLlA/znS/tZs/MYi2rL+MurarnxwklMLMgZcL362lIe+PVOXt7bxpWzypNVvsioUWiIxGlmxQTuvn4Of3vdbLp6IuRmhQdd54rpJZhFL8eu0JB0oKOnRIbIzOIKDIDCnEwumFKoy5ZI2lBoiCTZwppSNuxtpbs3kupSREZMoSGSZItqS+nqifDKgeOpLkVkxBQaIklWp8uxSxpRaIgkWUVBNjMq8nWbWUkLIwoNM/tbM9tiZq+Y2Q/MLMfMas3sBTPbYWY/NLOsoG928L4xWF4Ts517g/btZnZjTPvSoK3RzO4ZSa0iqbSotpR1r7cQiaTXeVEy/gw7NMxsGvBJoM7dLwLCwK3A14BvuvtsoBW4PVjldqDV3WcB3wz6YWYXBOtdCCwF/snMwmYWBr4D3ARcAHwo6Csy5iysKaW9q5fth3UdKhnbRjo8lQHkmlkGkAccBK4BfhQsXw7cEry+OXhPsPxaM7Og/XF3P+3uu4FGoD54NLr7LnfvBh4P+oqMOfW1mteQ9DDs0HD3/cD/BfYSDYvjwHqgzd17g25NwLTg9TRgX7Bub9C/LLb9jHXO1i4y5lSW5DG1KIcXdb6GjHEjGZ4qIfrNvxaYCuQTHUo6U/8grp1l2VDbB6rlDjNrMLOG5ubmwUoXSYn62uhtZtPtem8yvoxkeOo6YLe7N7t7D/AT4EqgOBiuAqgEDgSvm4AqgGB5EdAS237GOmdrfxt3f9Dd69y9rqKiYgS7JJI8C2tLae44zZ5jp1JdisiwjSQ09gKLzSwvmJu4FngVeBb4k6DPMuDJ4PWK4D3B8mc8+pVrBXBrcHRVLTAbeBFYB8wOjsbKIjpZvmIE9Yqk1CLNa0gaGMmcxgtEJ7Q3AJuDbT0IfA6428waic5ZPBys8jBQFrTfDdwTbGcL8ATRwPkFcKe79wXzHncBK4GtwBNBX5ExaWbFBErzszSvIWOa7qchMor+6vsNbD3YwXOffXeqSxF5i3jvp6EzwkVGUX1tGXtbTnHoeFeqSxEZFoWGyCiq778OVYKGqNydZ7cd4c5/38D9q3ew9WC7js6SpNLwlMgo6u2LcNkXV9HdFyEzZPS5E4lAnzuZYWNCdiaFORkU5GRQkJPJ5KIc3jG7nKvnVFCcl/XGdtyd53Yc5ZurXuPlfW0U52VyvLMHd5hWnMv1F0ziuvmTWDKzjHBooKPXRd5K9wgXOU+t2HiAl/a2EjYjHDJCISNsRndfhI6uHtq7euno6qWjq4fXj56k9VQPIYMF1SW8e95Easvzefj53dH7mRfnctc1s/jjBZW0dXbz7LYjrHr1CM83NtPVE+Fvr5vDp66bnepdljFAoSGSBvoizsamNn697QjPbm9m8/7oPTkmF+Zw5zWz+C91lWRnvP0ugp3dfXzou2tx4Mk7rxrlqmUsijc0dI9wkfNYOGQsqC5hQXUJd98wlyMdXWw72EF9bSk5mWe/5WxuVphFM0p55PnddPX0nbOvyFBoIlxkDJlYkMM751TEFQILqkvo6XO26I6BkkAKDZE0taC6BIANe9pSXImkE4WGSJqqKMimqjSXDXtbU12KpBGFhkgau7yqhA17W3XuhiSMQkMkjS2oLuZw+2kO6gx0SRCFhkgaWzA9mNfQEJUkiEJDJI3Nm1xIdkZIk+GSMAoNkTSWlRHiksoiXtqnXxqSGAoNkTS3oLqELfvbOd3bl+pSJA0oNETS3OXVxXT3RXhlf3uqS5E0oNAQSXP9J/m9pMlwSQCFhkiam1iYw7TiXF7aq8lwGTmFhsg4sGB6iQ67lYRQaIiMA5dXFXPweBcHj3emuhQZ40YUGmZWbGY/MrNtZrbVzJaYWamZrTKzHcFzSdDXzOx+M2s0s01mtiBmO8uC/jvMbFlM+xVmtjlY534z0y3IRIah/yQ/DVHJSI30l8b/A37h7vOAS4GtwD3AanefDawO3gPcBMwOHncADwCYWSlwH7AIqAfu6w+aoM8dMestHWG9IuPSBVMKycoIsWGPhqhkZIYdGmZWCLwTeBjA3bvdvQ24GVgedFsO3BK8vhl41KPWAsVmNgW4EVjl7i3u3gqsApYGywrdfY1Hr7b2aMy2RGQIsjJCXDytSPMaMmIj+aUxA2gG/tXMXjKzh8wsH5jk7gcBgueJQf9pwL6Y9ZuCtnO1Nw3Q/jZmdoeZNZhZQ3Nz8wh2SSR9Lagu5pUDOslPRmYkoZEBLAAecPfLgZO8ORQ1kIHmI3wY7W9vdH/Q3evcva6iouLcVYuMUwuqS+jujfDqAZ3kJ8M3knuENwFN7v5C8P5HREPjsJlNcfeDwRDTkZj+VTHrVwIHgvZ3ndH+66C9coD+IjIMlwcn+f1800GOd/bQ3tVLe2cP7V09tHf20tH1ZltHVw+VJXncff0casrzU1y5nE+GHRrufsjM9pnZXHffDlwLvBo8lgFfDZ6fDFZZAdxlZo8TnfQ+HgTLSuArMZPfNwD3unuLmXWY2WLgBeA24NvDrVdkvJtclENVaS4PPb+bh57f/ZZlmWGjKDeTwpxMCnIyKMjJZPXWwzz9ykFuW1LDJ6+ZTVFeZooql/PJSH5pAPw18JiZZQG7gI8QHfJ6wsxuB/YCfxr0fQp4D9AInAr6EoTDl4B1Qb8vuntL8PrjwPeAXODp4CEiw7T8I/XsbTlFYRAQhbkZFOZkkp0R4swj2o90dPGNX77GI7/bzY83NPGpa2fzF4unkxnW6V3jmaXbbSDr6uq8oaEh1WWIpI1XD7Tz5ade5XeNxyjOy2R6WT5VJblUluRRVZpLdWkeMyomMLUo523BI2OHma1397pB+yk0RGQw7s6vtzfzy1cPsa+lk6bWU+xv66Sn782/H3lZYWZWTGDWxAnMn1LAX15ZS1aGfpWMFfGGxkiHp0RkHDAz3j1vIu+eN/GNtr6Ic7i9iz3HTrGz+QSNR06ws/kEa3Ye46cv7aeqJI+bLp6SwqolGRQaIjIs4ZAxtTiXqcW5LJlZ9kZ7T1+ES7/wS9bsOqbQSEP67SgiCZUZDlFXU8rvdx5LdSmSBAoNEUm4JTPKaDxygiMdXakuRRJMoSEiCXdlMFy1dlfLID1lrFFoiEjCXTi1kILsDNZoiCrtKDREJOEywiHqa0tZu0uhkW4UGiKSFEtmlrH76EndLTDNKDREJCn6D8PVEFV6UWiISFLMn1xIcV6mQiPNKDREJClCIWNRbSlrNK+RVhQaIpI0S2aU0dTayb6WU6kuRRJEoSEiSbNkZjmgeY10otAQkaSZM2kCZflZGqJKIwoNEUkaM2PxzDLW7DxGut2GYbxSaIhIUi2ZUcah9i5eP6Z5jXSg0BCRpOo/X+P3O4+muBJJBIWGiCTVjPJ8JhZkazI8TSg0RCSpzIwrZ5axdleL5jXSwIhDw8zCZvaSmf0seF9rZi+Y2Q4z+6GZZQXt2cH7xmB5Tcw27g3at5vZjTHtS4O2RjO7Z6S1ikhqLJlZxtETp2k8ciIp2+/pi/Da4Q6efHk/96/eweF23ccjWRJxu9dPAVuBwuD914BvuvvjZvbPwO3AA8Fzq7vPMrNbg34fNLMLgFuBC4GpwK/MbE6wre8A1wNNwDozW+HuryagZhEZRUtmRM/X+Kvvr2dGRT4VBdlUTMimvCCbgpwMcjPDZGeGyQ0eGWHDsLdso6cvQsupblpPdtNysptjJ7s53N7F9kMd7Dh8gu6+yBt9T/f28Zkb543qPo4XIwoNM6sE/hD4MnC3mRlwDfBnQZflwOeJhsbNwWuAHwH/GPS/GXjc3U8Du82sEagP+jW6+67gsx4P+io0RMaYqtJc/vqaWWxsOs7+ti5e3neclpOniYxgtCojZJRPyGb2pAl85Koa5k0pYP6UQj77o02s292auOLlLUb6S+NbwGeBguB9GdDm7r3B+yZgWvB6GrAPwN17zex40H8asDZmm7Hr7DujfdFARZjZHcAdANXV1SPYHRFJBjPj0zfMfUtbX8RpOdnNidO9dPX00dnTR1d3H129fXT3vj1NMkJGSX4WpcGjMCeD6PfOt1o8o4zv/e51unr6yMkMJ22fxqthh4aZvRc44u7rzexd/c0DdPVBlp2tfaD5lgG/l7j7g8CDAHV1dZppExkDwiGLDlMVZCd0uwtrSnnwuV1sajpOfW1pQrctI5sIvwp4v5m9DjxOdFjqW0CxmfWHUSVwIHjdBFQBBMuLgJbY9jPWOVu7iMhZLawpAeDF3TrENxmGHRrufq+7V7p7DdGJ7Gfc/c+BZ4E/CbotA54MXq8I3hMsf8ajx9+tAG4Njq6qBWYDLwLrgNnB0VhZwWesGG69IjI+FOdlMXdSAS/sbkl1KWkpGedpfI7opHgj0TmLh4P2h4GyoP1u4B4Ad98CPEF0gvsXwJ3u3hfMi9wFrCR6dNYTQV8RkXOqry1lw55WemOOqJLEsHQ72aaurs4bGhpSXYaIpNCKjQf45A9eYsVdV3FJZXGqyxkTzGy9u9cN1k9nhItI2qmviU6Av6ghqoRTaIhI2plclEN1aZ5CIwkUGiKSluprS1n3uq53lWgKDRFJS/U1pbSe6kna9a7GK4WGiKSl/hP7XnxdQ1SJpNCfkUTKAAAH/UlEQVQQkbQ0vSyPioJszWskmEJDRNKSmVFfW8qLuzWvkUgKDRFJW4tqSzl4vIum1s5Ul5I2FBoikrYWBudrrNO8RsIoNEQkbc2dVEBhTobmNRJIoSEiaSsUMhbWlCo0EkihISJprb62lF1HT9LccTrVpaQFhYaIpLWFwfka31+7h01NbRw9cVpHU43ASG/3KiJyXrt4WhHlE7K4f/UO7l+9A4DsjBDTinOZWpz75nNJLlOLc8jJDNPVHb397KnguSt4faq7/3UvfREnPyuDCTkZTMjOoCAng6qSPK6cVZ7iPU4uhYaIpLXMcIjffObd7D56kgNtnexv63zjeX9bF89sPzKkoauczBB5WRmEQ8bJ072c6u57y/Jf3X01syZOSPRunDcUGiKS9vKzM7hoWhEXTSsacHlXTx+Hjnexv62Tnr4IuZlhcrPC5GaGyQle52WFyckIEwrZW9btizgnTvey99gp3vePz7N662GFhohIOsvJDFNTnk9Nef6Q1w2HjKLcTC6uLGL+lEJWbz3CX109MwlVnh80ES4ikiDXzZ9Iw54WWk92p7qUpFFoiIgkyLXzJxFx+PVrR1JdStIMOzTMrMrMnjWzrWa2xcw+FbSXmtkqM9sRPJcE7WZm95tZo5ltMrMFMdtaFvTfYWbLYtqvMLPNwTr3m5m9vRIRkfPDJdOKqCjI5ldbFRoD6QU+7e7zgcXAnWZ2AXAPsNrdZwOrg/cANwGzg8cdwAMQDRngPmARUA/c1x80QZ87YtZbOoJ6RUSSKhQyrpk7kee2N9PdG0l1OUkx7NBw94PuviF43QFsBaYBNwPLg27LgVuC1zcDj3rUWqDYzKYANwKr3L3F3VuBVcDSYFmhu6/x6Jk4j8ZsS0TkvHTt/Il0nO6lIU0vkpiQOQ0zqwEuB14AJrn7QYgGCzAx6DYN2BezWlPQdq72pgHaRUTOW38wu5ysjFDaDlGNODTMbALwY+Bv3L39XF0HaPNhtA9Uwx1m1mBmDc3NzYOVLCKSNHlZGVw5s4zV2w6n5eVKRhQaZpZJNDAec/efBM2Hg6Elguf+uG0CqmJWrwQODNJeOUD727j7g+5e5+51FRUVI9klEZERu3b+JPYcO8XO5hOpLiXhRnL0lAEPA1vd/Rsxi1YA/UdALQOejGm/LTiKajFwPBi+WgncYGYlwQT4DcDKYFmHmS0OPuu2mG2JiJy3rp0XHZVPxyGqkfzSuAr4MHCNmb0cPN4DfBW43sx2ANcH7wGeAnYBjcB3gU8AuHsL8CVgXfD4YtAG8HHgoWCdncDTI6hXRGRUTC3O5YIphazeejjVpSTcsC8j4u7PM/C8A8C1A/R34M6zbOsR4JEB2huAi4Zbo4hIqlw3fyL/+GwjrSe7KcnPSnU5CaMzwkVEkiBdzw5XaIiIJMHFaXp2uK5yKyKSBKGQce28ifx800G6eyNkZYz8O3ok4kTciTjBc/R1X8Rxd/KyMhLyOeei0BARSZJr5k3k8XX7uPYbvyYjFHrzD33krX/03Z2+SBAGQTD0DbBsMMs/Ws/Vc5J72oFCQ0QkSa6eW8GfLarmeGcPYTNCBiEzzIxw6M3X/e3hkGFGtG/IguWx76P9QsYby8NBn5AZM4ZxP5ChUmiIiCRJdkaYr/zRxakuI6E0ES4iInFTaIiISNwUGiIiEjeFhoiIxE2hISIicVNoiIhI3BQaIiISN4WGiIjEzdLtdoRm1gzsGebq5cDRBJYzFmifxwft8/gwkn2e7u6DXoMk7UJjJMyswd3rUl3HaNI+jw/a5/FhNPZZw1MiIhI3hYaIiMRNofFWD6a6gBTQPo8P2ufxIen7rDkNERGJm35piIhI3BQaATNbambbzazRzO5JdT3JZmaPmNkRM3sl1bWMFjOrMrNnzWyrmW0xs0+luqZkM7McM3vRzDYG+/yFVNc0WswsbGYvmdnPUl3LaDCz181ss5m9bGYNSfscDU9F/+cCXgOuB5qAdcCH3P3VlBaWRGb2TuAE8Ki7X5TqekaDmU0Bprj7BjMrANYDt6T5f2cD8t39hJllAs8Dn3L3tSkuLenM7G6gDih09/emup5kM7PXgTp3T+q5KfqlEVUPNLr7LnfvBh4Hbk5xTUnl7s8BLamuYzS5+0F33xC87gC2AtNSW1VyedSJ4G1m8Ej7b4pmVgn8IfBQqmtJNwqNqGnAvpj3TaT5H5PxzsxqgMuBF1JbSfIFwzQvA0eAVe6e9vsMfAv4LBBJdSGjyIFfmtl6M7sjWR+i0IiyAdrS/tvYeGVmE4AfA3/j7u2prifZ3L3P3S8DKoF6M0vr4Ugzey9wxN3Xp7qWUXaVuy8AbgLuDIagE06hEdUEVMW8rwQOpKgWSaJgXP/HwGPu/pNU1zOa3L0N+DWwNMWlJNtVwPuDMf7HgWvM7N9SW1LyufuB4PkI8FOiw+4Jp9CIWgfMNrNaM8sCbgVWpLgmSbBgUvhhYKu7fyPV9YwGM6sws+LgdS5wHbAttVUll7vf6+6V7l5D9N/yM+7+FykuK6nMLD84uAMzywduAJJyZKRCA3D3XuAuYCXRydEn3H1LaqtKLjP7AbAGmGtmTWZ2e6prGgVXAR8m+s3z5eDxnlQXlWRTgGfNbBPRL0er3H1cHII6zkwCnjezjcCLwM/d/RfJ+CAdcisiInHTLw0REYmbQkNEROKm0BARkbgpNEREJG4KDRERiZtCQ0RE4qbQEBGRuCk0REQkbv8fKY70lV15W1MAAAAASUVORK5CYII=\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "#calculate CVA\n",
    "for i in xrange(len(T)-1):\n",
    "    sum=sum+0.5*(EE[i]*crvToday.discount(T[i])+EE[i+1]*crvToday.discount(T[i+1]))*(exp(-S*T[i]/(1.0-R))-exp(-S*T[i+1]/(1.0-R)))\n",
    "CVA=(1.0-R)*sum\n",
    "\n",
    "print \"\\nCVA=\",CVA\n",
    "\n",
    "plot(T,EE)\n",
    "print \"\\nEE\\n\",npvMat\n",
    "show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 2",
   "language": "python",
   "name": "python2"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 2
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython2",
   "version": "2.7.15"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
