# -*- coding: utf-8 -*-
"""
simulates the pressure and flow pattern during mechanical ventilation
and computes iPEEP and VT, depending on tube resistance and patient's lung mechanics'
"""

import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
# from scipy.optimize import fsolve
# faster but without errorchecks
def fsolve(func, x0, tol=1e-6, max_iter=100):
    x1 = x0
    x2 = x0 * 1.1 if x0 != 0 else 1e-3
    for _ in range(max_iter):
        f1 = func(x1)
        f2 = func(x2)
        denom = f2 - f1
        if abs(denom) < tol:
            break
        x_new = x2 - f2 * (x2 - x1) / denom
        if abs(x_new - x2) < tol:
            x2 = x_new
            break
        x1, x2 = x2, x_new
    return x2
from scipy.signal import find_peaks
plt.close('all')

dt = 0.001  # delta t for numerical integration in s

# %% Patient parameters
C = 50  # patient's total (static) respiratory compliance in ml/mbar [5 ... 500]
R_aw = 1  # patient's airway resistance in mbar/(L/s) [1 ... 100]
ti = 1  # inspiratory time (in sec) [0.1 ... 10]
te = 1  # expiratory time (in sec) [0.1 ... 10]
PMusMax = 5.0  # maxumum inspiratory muscular pressure of the patient during spontaneous breathing in mbar [0 ... 100]
mus_relax_time = 0.4  # time for relaxing the muscles at the end of the inspiration [0 ... te]
if mus_relax_time > te:
    mus_relax_time = te  # muscles must be relaxed at the end of the expiratory phase
NoB = 5  # number of breaths to be simulated [1 ... 50]

# Ventilator parameters
VC = False  # volume controlled (VC = True) or pressure controlled (VC = False) ventilation
SP = True  # set SP (spontaneously breathing) to True if the timing of inspiration is determined by the airway pressure falling below the trigger.
if VC:
    SP = False  # when VC is true, SP must be false
    PMusMax = 0.0  # and the pateint is fully sedated

IPS = 15  # inspiratory pressure support (above PEEP) in mbar [0 ... 50]
PEEP = 5  # PEEP in mbar [0 ... 50]
VT = 500  # tidal volume in ml [0 ... 2500]
PTrigger = 0.5  # Trigger pressure below PEEP for initiating an inspiratory pressure support cycle [0 ... 10]
if IPS > 0.5:
    FlowTermination = 0.25  # expiration is initiated during IPS when the flow falls below FlowTermination*V'max  [0.01 .. 0.99]
else:
    FlowTermination = 0.01
t_rise = 0.2  # pressure rise time of the pressure support in s [0 .. 10]

# resistance and demand flow control parameters of the ventilator
R_PEEPvalve = 2  # resistance of the expiratory (PEEP) valve [0 .. 20]
Ideal = False  # set to True for an ideal ventilator where the airway pressure is delivered exactly as intended
# set to False to model the demand-flow characteristics of a real ventilator, which is
# determined by the Integral control parameter and the cutoff frequency
# this setting has no effect in VC mode
Integral = 2  # Integral control parameter for the demand flow controler [0.01 .. 10]
f_cutoff = 6  # cutoff frequency of the demand flow controler [0.1 .. 50]

# ETT - tube resistance parameters
ETT_ID = 7.0  # inner diameter of the endotracheal tube [6 .. 12]
Cath_OD = 0  # outher diameter of a bronchoscopy tube [0,3.8 .. 6]
if Cath_OD >= ETT_ID:
    Cath_OD = ETT_ID - 0.5  # outer diameter of the bronchoscopy catheter must be 0.5 mm smaller than the inner tube diameter

# %% initialization
k1 = 0.72
if Cath_OD > 0:
    beta = -5 #obstructed tube
else:
    beta = -3.6 #unobstructed tube


deff = np.sqrt(ETT_ID ** 2 - Cath_OD ** 2)
d0 = 10 #mm reference brochoscope

k1_in = k1 * (deff/d0) ** beta
k2_in = 4 * k1_in #changed from 4.3 to 4
k1_ex = k1 * (deff/d0) ** beta
k2_ex = 4 * k1_ex #changed from 4.3 to 4

i_e_r = ti / te  # i:e ratio
rr = 60 / (ti + te)  # respiratory rate (breaths/min)
T = 60 / rr  # period time of a cycle
flow_i = VT / ti / 1000.0  # inspiratory flow during VC ventilation in l/s
smoothing = np.exp(-2 * np.pi * f_cutoff * dt)


# %%
def iterate_flow_in(flow):
    R = R_aw + (k1_in * np.abs(flow) + k2_in * (flow ** 2)) / (np.abs(flow) + 1e-6)
    predicted_flow = (PS - ((V[-1] + flow * dt) * 1000.0 / C - P_mus[-1])) / R
    return flow - predicted_flow


def iterate_flow_ex(flow):
    R = R_aw + R_PEEPvalve + (k1_ex * np.abs(flow) + k2_ex * (flow ** 2)) / (np.abs(flow) + 1e-6)
    predicted_flow = (P_mus[-1] - (V[-1] + flow * dt) * 1000.0 / C) / R
    return flow - predicted_flow


# %% initialization

F = [0]
V = [0]
P_aw = [PEEP]
P_trach = [PEEP]
P_alv = [PEEP]
P_mus = [0]
t = [0]
support_start = []
support_end = []

runtime = 0
P = PEEP
R = R_aw
flow = 0.001
FlowPeakIns = 0
isInspiration = False
ErrorIntegral = 0

# %% model simulation
overtime = 0.5  # continue to simulate beyond the last breath by some runtime
while runtime < (NoB * T + overtime):
    t_cycle = runtime % T  # determine cycle time
    if t_cycle < dt:
        # flow = 0.001
        print('breath # ', int(runtime / T))
        # IPS = IPS + 1

    # compute muscle effort
    if t_cycle <= ti:
        # P_mus.append(PMusMax * 1 * (0 + np.sin((t_cycle) / ti * np.pi / 2)))
        P_mus.append(PMusMax * 0.5 * (1 - np.cos((t_cycle) / ti * np.pi / 1)))
    else:
        if t_cycle < ti + mus_relax_time:
            P_mus.append(PMusMax * 0.5 * (1 + np.cos((t_cycle - ti) / mus_relax_time * np.pi / 1)))
        else:
            P_mus.append(0)

    # check if inspiration or expiration
    if (SP == False):  # pressure or volume controlled ventilation, but no spontaneous breathing
        if t_cycle <= ti:  # inspiration
            isInspiration = True
            TriggertimeIns = runtime - t_cycle
        else:
            isInspiration = False
            ErrorIntegral = 0

    if isInspiration:
        if VC:  # volume controlled ventilation
            flow = flow_i
        else:  # Pressure controlled
            if (SP == False):  # Pressure controlled without spontaneous breathing
                if t_cycle <= t_rise:
                    PS = (IPS * t_cycle / t_rise)
                else:
                    PS = IPS


            else:  # Pressure controlled with spontaneous breathing (IPS)
                if runtime <= TriggertimeIns + t_rise:
                    PS = (IPS * (runtime - TriggertimeIns) / t_rise)
                else:
                    PS = IPS

            if Ideal:
                # compute ideal flow
                flow = fsolve(iterate_flow_in, 1)
                if flow < 0.0:
                    flow = 0.0
            else:
                # compute demand flow
                Pmeas = P_aw[-1]
                Error = (PEEP + PS) - Pmeas
                ErrorIntegral = ErrorIntegral + Error * Integral * dt
                flow = smoothing * F[-1] + (1 - smoothing) * ErrorIntegral

            if ((flow > FlowPeakIns) and (flow > 0)):
                FlowPeakIns = flow

            # check the termination criterion for IPS mode
            if (flow < FlowTermination * FlowPeakIns) and SP:
                # R = R_aw + (k1_ex * abs(0.001)**k2_ex)/abs(0.001)   #secant ETT resistance
                # die()
                isInspiration = False
                support_end.append(int(runtime / dt))
                ErrorIntegral = 0
                R = R_aw + R_PEEPvalve + k1_ex
                FlowPeakIns = -2;

            else:
                R = R_aw + (k1_in * np.abs(flow) ** k2_in) / (
                            np.abs(flow) + 1e-6)  # secant ETT resistance


    else:  # if not(isInspiration):
        PS = 0
        flow = fsolve(iterate_flow_ex, -1)
        if flow > 0.0:
            flow = 0.0
    F.append(flow)
    V.append(V[-1] + flow * dt)
    P_alv.append(PEEP + V[-1] * 1000.0 / C - P_mus[-1])
    P_trach.append(P_alv[-1] + flow * R_aw)
    if flow > 0:
        dP_ETT = k1_in * flow + k2_in * (flow ** 2)
    else:
        dP_ETT = k1_ex * flow - k2_ex * (flow ** 2)
    P_aw.append(P_trach[-1] + dP_ETT)
    runtime = runtime + dt
    t.append(runtime)

    if ((VC == False) and (SP == True)):  # IPS
        if (isInspiration == False):  # check if the ventilator recognizes an inspiratory effort
            if (P_aw[-1] < (PEEP - PTrigger)):  # trigger criterion for inspiration
                TriggertimeIns = runtime
                isInspiration = True
                support_start.append(int(runtime / dt))
        else:  # check if the ventilator recognizes an expiratory effort
            if (flow > 0) and (F[-1] < (FlowTermination * FlowPeakIns)):
                isInspiration = False
                support_end.append(int(runtime / dt))
                ErrorIntegral = 0
                R = R_aw + R_PEEPvalve + k1_ex
                FlowPeakIns = -2;

t = np.asarray(t)

# %% plotting of the results
# compute some parameters, based on the last breath
in_start = int((NoB - 1) * T / dt)
in_end = int(((NoB - 1) * T + ti) / dt)
ex_start = int(((NoB - 1) * T + ti + dt) / dt)
ex_end = int((NoB * T - dt) / dt)

plt.close('all')

fig, ax = plt.subplots(3, sharex='all', figsize=(18, 15))
plt.subplots_adjust(left=0.22)
plt.subplots_adjust(top=0.95)
plt.subplots_adjust(bottom=0.12)

Tpast = NoB  # breaths
Future = overtime  # s
a = int((-Tpast * T - Future) / dt)
b = int(-Future / dt)

C0 = '#1f77b4'
C1 = '#ff7f0e'
C2 = '#2ca02c'
C3 = '#d62728'
font = {'family': 'sans-serif',
        'sans-serif': ['Arial'],
        'weight': 'normal',
        'size': 16}
plt.rc('font', **font)
plt.rc('legend', fontsize=10)
plt.rc('axes', titlesize=18)
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

# t = t - t[int((-Tpast*T-overtime)/dt)]
P_mus = 0.0 - np.asarray(P_mus)
# (P_alv_line, ) = ax[0].plot(t[int((-8*T-2)/dt):], P_alv[int((-8*T-2)/dt):], color = C0, linewidth = 2)
ax[0].plot(t[a:b], t[a:b] * 0 + PEEP, color='gray', linewidth=1, linestyle='--')
(P_aw_line,) = ax[0].plot(t[a:b], P_aw[a:b], color=C3, linewidth=2)
(P_tr_line,) = ax[0].plot(t[a:b], P_trach[a:b], color=C0, linewidth=2)
(P_mus_line,) = ax[0].plot(t[a:b], P_mus[a:b], color='gray', linewidth=2)

(flow_line,) = ax[1].plot(t[a:b], F[a:b], color=C1, linewidth=2)
(vol_line,) = ax[2].plot(t[a:b], V[a:b], color=C2, linewidth=2)

ax[0].set_ylim(np.min(P_mus) - 1, np.max([np.max(P_aw), np.max(P_trach)]) + 1)
ax[1].set_ylim(np.min([-1, np.min(F) - 0.1]), np.max([1, np.max(F) + 0.1]))
ax[2].set_ylim(-0.02, np.max(V) + 0.02)

ax[0].set_ylabel('pressure (mbar)')
ax[1].set_ylabel('flow (L/s)')
ax[2].set_ylabel('volume (L)')
ax[2].set_xlabel('time (s)')

ax[0].grid(True, axis='y')
ax[1].grid(True, axis='y')
ax[2].grid(True, axis='y')

print('iPEEP    = %3.12f  mbar' % P_alv[in_start])
print('Palv_max = %3.12f  mbar' % P_alv[ex_start - 1])
print('Paw_max  = %3.12f  mbar' % P_aw[ex_start - 1])
V = np.asarray(V)
support_start = np.asarray(support_start)
support_end = np.asarray(support_end)
spontan_start = []
spontan_end = []
for i in range(Tpast):
    spontan_start.append(i * T / dt)
    spontan_end.append((i * T + ti) / dt)
spontan_start = np.asarray(spontan_start)
spontan_end = np.asarray(spontan_end)

# Loop over the start and end times
for start, end in zip(support_start, support_end):
    start = start * dt
    end = end * dt
    width = end - start  # Calculate the width of each rectangle
    height = ax[0].get_ylim()[1] - PEEP  # Dynamically covers the y-axis range

    # Create a rectangle patch for each time range
    rect0 = patches.Rectangle((start, PEEP), width, height, linewidth=1, edgecolor='black', facecolor='gray', alpha=0.1,
                              hatch='////')
    rect1 = patches.Rectangle((start, ax[1].get_ylim()[0]), width, height, linewidth=1, edgecolor='black',
                              facecolor='gray', alpha=0.1, hatch='////')
    rect2 = patches.Rectangle((start, ax[2].get_ylim()[0]), width, height, linewidth=1, edgecolor='black',
                              facecolor='gray', alpha=0.1, hatch='////')
    ax[0].add_patch(rect0)  # Add the rectangle to ax[0]

br = 0
for start, end in zip(spontan_start, spontan_end):
    start = start * dt
    end = end * dt
    width = end - start  # Calculate the width of each rectangle
    # height = ax[0].get_ylim()[1] - ax[0].get_ylim()[0]  # Dynamically covers the y-axis range
    height = PEEP - ax[0].get_ylim()[0]  # Dynamically covers the y-axis range

    # Create a rectangle patch for each time range
    rect0 = patches.Rectangle((start, ax[0].get_ylim()[0]), width, height, linewidth=1, edgecolor='black',
                              facecolor='gray', alpha=0.1, hatch='\\\\\\\\')
    rect1 = patches.Rectangle((start, ax[1].get_ylim()[0]), width, height, linewidth=1, edgecolor='black',
                              facecolor='gray', alpha=0.1, hatch='\\\\\\\\')
    rect2 = patches.Rectangle((start, ax[2].get_ylim()[0]), width, height, linewidth=1, edgecolor='black',
                              facecolor='gray', alpha=0.1, hatch='\\\\\\\\')
    ax[0].add_patch(rect0)  # Add the rectangle to ax[0]

    br = br + 1

ax[0].text(t[-1], PEEP - 3, 'P_mus', bbox=dict(facecolor='white', edgecolor='none'), color='gray', fontsize=14,
           zorder=10)
ax[0].text(t[-1], PEEP + 0, 'P_trach', bbox=dict(facecolor='white', edgecolor='none'), color=C0, fontsize=14, zorder=10)
ax[0].text(t[-1], PEEP + 6, 'P_aw', bbox=dict(facecolor='white', edgecolor='none'), color=C3, fontsize=14, zorder=10)

# compute minute ventilation
if SP == True:
    V = np.array(V)

    maxima, _ = find_peaks(V,distance=(ti + te - 0.1)/dt)
    minima, _ = find_peaks(-V, distance=(ti + te - 0.1) / dt)
    last_minimum = np.max(minima[minima < maxima[-1]])
    second_last_maximum = np.max(maxima[maxima < last_minimum])
    breath_volume = 0.5*(V[second_last_maximum] - V[last_minimum] + V[maxima[-1]] - V[last_minimum])
    Vmin = breath_volume * 60 / ((maxima[-1]-second_last_maximum) * dt)
    breath_duration = np.diff(support_start) * dt

else:
    breath_volume = 0.5 * (V[in_end] - V[in_start] + V[ex_start] - V[ex_end])
    breath_duration = T
    Vmin = breath_volume * 60 / T
print('rr       = %3.2f  br/min' % (60 / np.mean(breath_duration)))
print('VT       = %3.0f  ml' % (1000 * np.mean(breath_volume)))
print('Vmin     = %3.2f  L/min' % Vmin)

plt.show()
