
# coding: utf-8

# American Monte Carlo for Exposure Simulation
# Copyright (c) 2016 Matthias Groncki
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# - Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# This disclaimer is taken from the QuantLib license

# In[1]:


# import the used libraries
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import QuantLib as ql
get_ipython().magic(u'matplotlib inline')


# In[2]:


# Setting evaluation date
today = ql.Date(7,4,2015)
ql.Settings.instance().setEvaluationDate(today)


# In[3]:


# Setup Marketdata
rate = ql.SimpleQuote(0.03)
rate_handle = ql.QuoteHandle(rate)
dc = ql.Actual365Fixed()
yts = ql.FlatForward(today, rate_handle, dc)
yts.enableExtrapolation()
hyts = ql.RelinkableYieldTermStructureHandle(yts)
t0_curve = ql.YieldTermStructureHandle(yts)
euribor6m = ql.Euribor6M(hyts)
cal = ql.TARGET()


# In[4]:


# Setup a dummy portfolio with two Swaps
def makeSwap(start, maturity, nominal, fixedRate, index, typ=ql.VanillaSwap.Payer):
    """
    creates a plain vanilla swap with fixedLegTenor 1Y
    
    parameter:
        
        start (ql.Date) : Start Date
        
        maturity (ql.Period) : SwapTenor
        
        nominal (float) : Nominal
        
        fixedRate (float) : rate paid on fixed leg
        
        index (ql.IborIndex) : Index
        
    return: tuple(ql.Swap, list<Dates>) Swap and all fixing dates
    
        
    """
    end = ql.TARGET().advance(start, maturity)
    fixedLegTenor = ql.Period("1y")
    fixedLegBDC = ql.ModifiedFollowing
    fixedLegDC = ql.Thirty360(ql.Thirty360.BondBasis)
    spread = 0.0
    fixedSchedule = ql.Schedule(start,
                                end, 
                                fixedLegTenor, 
                                index.fixingCalendar(), 
                                fixedLegBDC,
                                fixedLegBDC, 
                                ql.DateGeneration.Backward,
                                False)
    floatSchedule = ql.Schedule(start,
                                end,
                                index.tenor(),
                                index.fixingCalendar(),
                                index.businessDayConvention(),
                                index.businessDayConvention(),
                                ql.DateGeneration.Backward,
                                False)
    swap = ql.VanillaSwap(typ, 
                          nominal,
                          fixedSchedule,
                          fixedRate,
                          fixedLegDC,
                          floatSchedule,
                          index,
                          spread,
                          index.dayCounter())
    return swap, [index.fixingDate(x) for x in floatSchedule][:-1]


def makeSwaption(swap, callDates, settlement):
    """
    Creates the swaption of the underlying swap.
    
    If there is only one callDate in the list of call dates it will be a European swaption
    otherwise a Bermudan.
    """
    if len(callDates) == 1:
        exercise = ql.EuropeanExercise(callDates[0])
    else:
        exercise = ql.BermudanExercise(callDates)
    return ql.Swaption(swap, exercise, settlement)


# Setup the Bermudan Swaption
# 1) Create a ATM plain vanilla swap using the helper function above and
# 
# 2) create a bermudan swaption with yearly exercise dates.

# In[257]:


settlementDate = today + ql.Period("2D")

swaps = [makeSwap(settlementDate,
                  ql.Period("5Y"),
                  1e6,
                  0.1,
                  euribor6m)
        ]

calldates = [ql.Date(7,4,2016), ql.Date(6,4,2017), ql.Date(5,4, 2018), ql.Date(5,4,2019)]
 
swaptions = [makeSwaption(swap, 
                          calldates, 
                          ql.Settlement.Physical) 
             for swap, fd in swaps]

calldates = [euribor6m.valueDate(d) for d in calldates]


# In[258]:


#%%timeit
# Setup pricing engine and calculate the npv of the underlying swap
engine = ql.DiscountingSwapEngine(hyts)
for swap, fixingDates in swaps:
    swap.setPricingEngine(engine)
    print("Swap NPV at time 0: %.4f" % swap.NPV())


# Setup the Gaussian Shortrate model (a.k.a Hull White model)
# Don't worry about calibration, assume we know the calbriated model parameter

# In[259]:


# Assume the model is already calibrated either historical or market implied
volas = [ql.QuoteHandle(ql.SimpleQuote(0.0075)),
         ql.QuoteHandle(ql.SimpleQuote(0.0075))]
meanRev = [ql.QuoteHandle(ql.SimpleQuote(0.002))]
model = ql.Gsr(t0_curve, [today+100], volas, meanRev, 16.)
process = model.stateProcess()


# Calculate the swaption price using an integral pricing engine

# In[260]:


swaptionEngine = ql.Gaussian1dSwaptionEngine(model)
for swaption in swaptions:
    swaption.setPricingEngine(swaptionEngine)
    print("Swaption NPV : %.2f" % swaption.NPV())


# Pricing with an Monte Carlo method
# Create a swap path pricer in Python
# Convert all Dates into times in years (using the same DayCounter as in the yieldTermStructure and store all fix cashflows in a numpy array.

# In[261]:


mcDC = yts.dayCounter()

def timeFromReferenceFactory(daycounter, ref):
    """
    returns a function, that calculate the time in years
    from a the reference date *ref* to date *dat* 
    with respect to the given DayCountConvention *daycounter*
    
    Parameter:
        dayCounter (ql.DayCounter)
        ref (ql.Date)
        
    Return:
    
        f(np.array(ql.Date)) -> np.array(float)
    """
    def impl(dat):
        return daycounter.yearFraction(ref, dat)
    return np.vectorize(impl)

timeFromReference = timeFromReferenceFactory(mcDC, today)


# In[262]:


def getFixedLeg(swap, t):
    """
    returns all future payment times and amounts of the fixed leg of the underlying swap
    
    Parameter:
        swap (ql.Swap)
        t (float) 
        
    Return:
        (np.array, np.array) (times, amounts)

    """
    fixed_leg = swap.leg(0)
    n = len(fixed_leg)
    fixed_times=[]
    fixed_amounts=[]
    npv = 0
    for i in range(n):
        cf = fixed_leg[i]
        t_i = timeFromReference(cf.date())
        if t_i > t:
            fixed_times.append(t_i)
            fixed_amounts.append(cf.amount())
    return np.array(fixed_times), np.array(fixed_amounts)


def getFloatingLeg(swap, t):
    """
    returns all future payment, fixing and accrual start and end times and amounts and nominals for all non fixed 
    periods of the floatiing leg
    
    Parameter:
        swap (ql.Swap)
        t (float) 
        
    Return:
        (np.array, np.array, np.array, np.array, np.array) (payment_times, accrual_period, accrual_start_time, accrual_end_time, nominals)

    """
    float_leg = swap.leg(1)
    n = len(float_leg)
    float_times = []
    float_dcf = []
    accrual_start_time = []
    accrual_end_time = []
    nominals = []
    for i in range(n):
        # convert base classiborstart_idx Cashflow to
        # FloatingRateCoupon
        cf = ql.as_floating_rate_coupon(float_leg[i])
        value_date = cf.referencePeriodStart()
        t_fix_i = timeFromReference(value_date)
        t_i = timeFromReference(cf.date()) 
        if t_fix_i >= t:
            iborIndex = cf.index()
            index_mat = cf.referencePeriodEnd()
            # year fraction
            float_dcf.append(cf.accrualPeriod())
            # calculate the start and end time
            accrual_start_time.append(t_fix_i)
            accrual_end_time.append(timeFromReference(index_mat))
            # payment time
            float_times.append(t_i)
            # nominals 
            nominals.append(cf.nominal())
    return np.array(float_times), np.array(float_dcf), np.array(accrual_start_time), np.array(accrual_end_time), np.array(nominals)

def getFixedFloatingPeriod(swap, t):
    """
    
    """
    float_leg = swap.leg(1)
    n = len(float_leg)
    for i in range(n):
        cf = ql.as_floating_rate_coupon(float_leg[i])
        value_date = cf.referencePeriodStart()
        t_fix_i = timeFromReference(value_date)
        t_i = timeFromReference(cf.date()) 
        if t_fix_i < t and t < t_i:
            iborIndex = cf.index()         
            index_mat = cf.referencePeriodEnd()
            # year fraction
            float_dcf = cf.accrualPeriod()
            # calculate the start and end time
            accrual_start_time = t_fix_i
            accrual_end_time = timeFromReference(index_mat)
            # payment time
            float_times = t_i
            # nominals 
            nominals = cf.nominal()
            return (float(float_times), float(float_dcf), float(accrual_start_time), float(accrual_end_time), float(nominals))
    return (float(t), 1., float(t), float(t), 0.)

def swapPathNPV(swap, t, timegrid):
    """
    Generate a path pricer. 
    
    The path pricer calculate the npv of the swap conditional 
    at time t on a given path of the short rate.
    """
    fixed_times, fixed_amounts = getFixedLeg(swap, t)
    float_times, float_dcf, accrual_start_time, accrual_end_time, nominals = getFloatingLeg(swap, t)
    df_times = np.concatenate([fixed_times, 
                           accrual_start_time, 
                           accrual_end_time, 
                           float_times])
    df_times = np.unique(df_times)
    # Store indices of fix leg payment times in 
    # the df_times array
    fix_idx = np.in1d(df_times, fixed_times, True)
    fix_idx = fix_idx.nonzero()
    # Indices of the floating leg payment times 
    # in the df_times array
    float_idx = np.in1d(df_times, float_times, True)
    float_idx = float_idx.nonzero()
    # Indices of the accrual start and end time
    # in the df_times array
    accrual_start_idx = np.in1d(df_times, accrual_start_time, True)
    accrual_start_idx = accrual_start_idx.nonzero()
    accrual_end_idx = np.in1d(df_times, accrual_end_time, True)
    accrual_end_idx = accrual_end_idx.nonzero()
    paytime_ffp, float_dcf_ffp, accrual_start_time_ffp, accrual_end_time_ffp, nominals_ffp = getFixedFloatingPeriod(swap, t)
    # Calculate NPV
    def calc(path):
        """
        Calculate the npv conditional on the given path 
        of the short rate.
        """
        if len(df_times)==0:
            return 0
        i = np.where(timegrid == t)[0][0]
        x_t = path[i]
        discount = np.vectorize(lambda T: model.zerobond(T, t, x_t))
        dfs = discount(df_times)
        # Calculate fixed leg npv
        fix_leg_npv = np.sum(fixed_amounts * dfs[fix_idx])
        # Estimate the index fixings
        index_fixings = (dfs[accrual_start_idx] / dfs[accrual_end_idx] - 1) 
        index_fixings /= float_dcf
        # Calculate the floating leg npv
        float_leg_npv = np.sum(nominals * index_fixings * float_dcf * dfs[float_idx])
        # Calculate the already fixed accrual period of the floating leg
        t_f = accrual_start_time_ffp
        i = np.where(timegrid == t_f)[0][0]
        x_f = path[i]
        df_e = model.zerobond(accrual_end_time_ffp, t_f, x_f)
        npv_accrualperiod = (1. / df_e - 1) * nominals_ffp * model.zerobond(paytime_ffp, t, x_t)
        # Calculate swap npv
        npv = float_leg_npv + npv_accrualperiod - fix_leg_npv
        return npv
    return calc


# In[263]:


# Convert call date to time
callTimes = timeFromReference(calldates)
callTimes


# In[264]:


swap = swaps[0][0]
swaption = swaptions[0]


# In[265]:


npv = swapPathNPV(swap, 0., np.array([0.]))(np.array([0.]))
print("Swap NPV at time 0: %.4f" % npv)
print("Error : %.8f" % (npv - swap.NPV()))


# Monte Carlo Simulation
# Generate time grid and paths

# In[266]:


def fixingdates(swap):
    leg = swap.leg(1)
    n = len(leg)
    fixing_dates = []
    for i in range(0, n):
        cf = ql.as_floating_rate_coupon(leg[i])
        value_date = cf.referencePeriodStart()
        fixing_dates.append(value_date)
    return fixing_dates

# Define evaluation grid
fixing_dates = fixingdates(swap)
fixing_times = timeFromReference(fixing_dates )

date_grid = [today + ql.Period(i, ql.Months) for i in range(0,66)] + calldates + fixing_dates

date_grid = np.unique(np.sort(date_grid))
time_grid = np.vectorize(lambda x: ql.Actual365Fixed().yearFraction(today, x))(date_grid)
time_grid = np.unique(time_grid)
dt = time_grid[1:] - time_grid[:-1]


# In[267]:


# Random number generator
seed = 1
urng = ql.MersenneTwisterUniformRng(seed)
usrg = ql.MersenneTwisterUniformRsg(len(time_grid)-1,urng)
generator = ql.InvCumulativeMersenneTwisterGaussianRsg(usrg)


# In[268]:


#%%timeit
# Generate N paths
M = 10000
m = len(time_grid)
x = np.zeros((M, m))
y = np.zeros((M, m))
numeraires = np.zeros((M, m))
                      
for n in range(0, M):
    numeraires[n, 0] = model.numeraire(0, 0)
    
for n in range(0,M):
    dWs = generator.nextSequence().value()
    j = 1
    for i in range(1, len(time_grid)):
        t0 = time_grid[i-1]
        t1 = time_grid[i]
        e = process.expectation(t0, 
                                x[n,i-1], 
                                dt[i-1])
        std = process.stdDeviation(t0,
                                   x[n,i-1],
                                   dt[i-1])
        x[n,i] = e + dWs[i-1] * std 
        e_0_0 = process.expectation(0,0,t1)
        std_0_0 = process.stdDeviation(0,0,t1)
        y[n,i] = (x[n,i] - e_0_0) / std_0_0
        numeraires[n, i] = model.numeraire(t1, y[n,i])
        #df_times_temp = df_times.copy()
        #df_times_temp[df_times_temp <= t1] = t1
        #dfs[n,i] = np.vectorize(lambda T: model.zerobond(T, t1, y[n,i]))(df_times_temp)


# In[269]:


swap_npvs = np.zeros((M, m))
swaption_npvs = np.zeros((M, m))
cont_value = np.zeros(numeraires[:,i].shape)
for i in range(m-1, 0, -1):
    t = time_grid[i]
    print(t)
    pricer = swapPathNPV(swap, t, time_grid)
    swap_npvs[:, i] = np.apply_along_axis(pricer, 1, y) / numeraires[:, 0]
    exercise_values = np.zeros(numeraires[:,i].shape)
    if t in callTimes:
        exercise_values = swap_npvs[:, i].copy()
        exercise_values[exercise_values < 0] = 0
    states = y[:, i]
    Y = np.column_stack((states, states**2, states**3, states**4))
    Y = sm.add_constant(Y)
    ols = sm.OLS(cont_value, Y)
    ols_result = ols.fit()
    cont_value_hat = np.sum(ols_result.params * Y, axis=1)
    swaption_npvs[:,i] = np.maximum(cont_value_hat, exercise_values)
    if t in callTimes:
        print("Update")
        cont_value = np.maximum(cont_value_hat, exercise_values)
        swaption_npvs[cont_value_hat < exercise_values, i:] = swap_npvs[cont_value_hat < exercise_values, i:].copy()
swaption_npvs[:,0] = np.mean(cont_value)
swap_npvs[:, 0] = np.apply_along_axis(swapPathNPV(swap, 0, time_grid), 1, y) / numeraires[:, 0]


# In[270]:


swaption_npvs *= numeraires[0,0]
swap_npvs *= numeraires[0,0]


# In[271]:


plt.figure(figsize=(10,5))
for i in range(1, 100):
    plt.plot(time_grid, swaption_npvs.T[:,i])
    plt.title("Simulted Swaption Exposures")
    plt.xlabel("Time")
    plt.ylabel("NPV")


# In[272]:


plt.figure(figsize=(10,5))
for i in range(1, 100):
    plt.plot(time_grid, swap_npvs.T[:,i])
    plt.title("Simulted Swap Exposures")
    plt.xlabel("Time")
    plt.ylabel("NPV")


# In[273]:


# Expected Exposure
swap_npvs[swap_npvs<0] = 0
swaption_npvs[swaption_npvs<0] = 0
EE_swaption = np.mean(swaption_npvs, axis=0)
EE_swap = np.mean(swap_npvs, axis=0)


# In[274]:


plt.figure(figsize=(10,5))
plt.title("Positive Expected Exposures")
plt.plot(time_grid, EE_swaption)
plt.plot(time_grid, EE_swap)
plt.xlabel("Time")
plt.ylabel("NPV")
plt.legend(["Swaption", "Underlying Swap"])


# CVA

# In[275]:


# Setup Default Curve 
pd_dates =  [today + ql.Period(i, ql.Years) for i in range(11)]
hzrates = [0.02 * i for i in range(11)]
pd_curve = ql.HazardRateCurve(pd_dates,hzrates,ql.Actual365Fixed())
pd_curve.enableExtrapolation()


# In[276]:


# Plot curve
# Calculate default probs on grid *times*
times = np.linspace(0,30,100)
dp = np.vectorize(pd_curve.defaultProbability)(times)
sp = np.vectorize(pd_curve.survivalProbability)(times)
dd = np.vectorize(pd_curve.defaultDensity)(times)
hr = np.vectorize(pd_curve.hazardRate)(times)
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
ax1.plot(times, dp)
ax2.plot(times, sp)
ax3.plot(times, dd)
ax4.plot(times, hr)
ax1.set_xlabel("Time in years")
ax2.set_xlabel("Time in years")
ax3.set_xlabel("Time in years")
ax4.set_xlabel("Time in years")
ax1.set_ylabel("Probability")
ax2.set_ylabel("Probability")
ax3.set_ylabel("Density")
ax4.set_ylabel("HazardRate")
ax1.set_title("Default Probability")
ax2.set_title("Survival Probability")
ax3.set_title("Default density")
ax4.set_title("Harzard rate")


# In[277]:


# Calculation of the default probs
defaultProb_vec = np.vectorize(pd_curve.defaultProbability)
dPD = defaultProb_vec(time_grid[:-1], time_grid[1:])


# In[278]:


# Calculation of the CVA
recovery = 0.4
CVA = (1-recovery) * np.sum(EE_swaption[1:] * dPD)
CVA


# In[279]:


# Calculation of the CVA
recovery = 0.4
CVA = (1-recovery) * np.sum(EE_swap[1:] * dPD)
CVA


# 
# 
# 
# fit GP

# In[31]:


import os
os.environ['TF_C_API_GRAPH_CONSTRUCTION']='0'
import tensorflow as tf
import edward as ed
import matplotlib.pyplot as plt
import numpy as np
from edward.models import Normal
import seaborn as sns
plt.style.use('ggplot')
from edward.models import MultivariateNormalTriL
from edward.util import rbf


# In[ ]:





# In[242]:


# CVA vs different fix rate of swap, change from 0.01 to 0.1 
CVA_fix = np.array([5901.3797828005836,2883.5699104367236,913.4,204.82055378641272,38.164562723122458,5.8741266072060876,0.78326262396065416,0.057947311871232809,0,0])
plt.plot(CVA_fix) 


# In[296]:


# Polynomial regression with 5 datasets
x_poly = np.arange(0.01, 0.1, 0.02)
y_poly = np.array([5901.3797828005836, 913.4, 38.164562723122458, 0.78326262396065416, 0])    # manually precalculated 
reg = np.polyfit(x_poly, y_poly, 4)
CVA_reg = np.polyval(reg, x_poly) 
plt.plot(CVA_reg)


# In[302]:


#Fit GP
x = np.array([ 0.01 ,  0.011,  0.012,  0.013,  0.014,  0.015,  0.016,  0.017, 0.018,  0.019,  0.02 ,  0.021,  0.022,  0.023,  0.024,  0.025, 0.026,  0.027,  0.028,  0.029,  0.03 ,  0.031,  0.032,  0.033,  0.034,  0.035,  0.036,  0.037,  0.038,  0.039,  0.04 ,  0.041, 0.042,  0.043,  0.044,  0.045,  0.046,  0.047,  0.048,  0.049, 0.05 ,  0.051,  0.052,  0.053,  0.054,  0.055,  0.056,  0.057,  0.058,  0.059,  0.06 ,  0.061,  0.062,  0.063,  0.064,  0.065, 0.066,  0.067,  0.068,  0.069,  0.07 ,  0.071,  0.072,  0.073,  0.074,  0.075,  0.076,  0.077,  0.078,  0.079,  0.08 ,  0.081, 0.082,  0.083,  0.084,  0.085,  0.086,  0.087,  0.088,  0.089, 0.09 ,  0.091,  0.092,  0.093,  0.094,  0.095,  0.096,  0.097, 0.098,  0.099,  0.1  ])
y = np.array([ 5901.3797828005836, 5569.8116458886179, 5243.3097029338842, 4922.3979386479214, 4607.5602210838551, 4299.5616852465564, 3999.1228967684647, 3706.7000907939928, 3422.8524275438476, 3148.2089286515843, 2883.5699104367236, 2629.537346227004, 2386.823136855875, 2156.014865360778, 1937.5984519417589, 1732.0833013384681, 1539.9393148366837,1361.5937012019269, 1197.5548932538395, 1048.1133118295052, 913.40656293382506, 793.4222251934309, 687.46147354480217, 594.46771084750969, 513.25449732168795, 442.40506161634249, 380.71257781868701, 326.98840909978213, 280.34691845840149, 239.88019466937141, 204.82055378641272, 174.55792656179881, 148.45235396283485, 125.97053659837852, 106.71776671789708, 90.278915352756698, 76.251110323898402, 64.300845261507632, 54.138179894453842, 45.488719497591319, 38.164562723122458, 31.98223527952819, 26.799747864388653, 22.438689471651447, 18.702216397880179, 15.509269810375406, 12.821960656286787, 10.56740847499926, 8.6862115384909551, 7.1410532666733078, 5.8741266072061009, 4.8329702671515413, 3.9738436699678594, 3.2591733644758158, 2.6681698442667749, 2.1824912932366822, 1.7861856622619241, 1.4620874326173956, 1.1909415126145204, 0.96550097869537277, 0.7832626239606888, 0.63643603110208258, 0.5161794857916463, 0.41787126922188944, 0.33462414435446747, 0.26336250993779964, 0.20478039696807784, 0.15558621751956891, 0.11537030948219035, 0.082752956812375225, 0.057947311871237493, 0.040262272662531871, 0.028471943469182071, 0.019720930732686406, 0.013142989767516253, 0.0075443601354315479, 0.0030754044425918806, 0.0007066804798661916, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
plt.plot(x,y)  # precalculated data


# In[430]:


sess = ed.get_session()
training_number = 17
idx = np.linspace(0, len(x)-1, training_number, dtype = 'int')

x_train = np.array(x[idx], dtype='float32').reshape(training_number, 1)
y_train = np.array(y[idx], dtype='float32').reshape(training_number, 1)
x_test = np.array(x, dtype='float32').reshape(len(x), 1)
y_test = np.array(y, dtype='float32').reshape(len(x), 1)


# In[431]:


# mean
Kernel = rbf(x_train).eval()
K_noise = Kernel + np.eye(training_number) * 0.01  # without noise, the cov band converge to 0 at the training points 
k_s = rbf(x_test, x_train).eval()                   
L = np.linalg.cholesky(K_noise)
alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train))
predict_mean = np.dot(k_s, alpha)

# cov
v = np.linalg.solve(L, k_s.T)
var = rbf(x_test).eval() - np.dot(v.T, v)


# In[433]:


# plot with cov band
up = predict_mean.reshape(len(x),) - 2 * (np.sqrt(np.diag(var)))   # 95% confident interval 
down = predict_mean.reshape(len(x),) + 2 * (np.sqrt(np.diag(var)))

plt.figure(facecolor='white', edgecolor='black')
plt.plot(x_test, predict_mean, color = 'red', label = 'GP Prediction')
plt.plot(x_test, y_test, color = 'black', label = 'Analytical Model')
plt.scatter(x_train, y_train, s = 150, color = 'black', marker = "+")
plt.fill_between(x, up, down, color = 'blue', alpha=0.9)
plt.grid(True)
#plt.ylim(-1000, 6100)
plt.xlabel('Fixed rate')
plt.ylabel('CVA')
plt.legend(loc = 'best', prop={'size':10})


# In[434]:


Kernel


# In[439]:


plt.fill_between(x, up, down, color = 'blue', alpha=0.9)


# In[442]:


rbf([[0.2, 0.01]], [[0.04, 0.11]]).eval()

