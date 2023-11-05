import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
import os.path
import xlwt
import csv

########### Transition Pathways  ###########
GET_case = ['Baseline', 'Linear', 'Delayed', 'Fast']
cs = 3
########### PARAMETERS  ###########
T = 50
tstep  = 5                                          #Time Step

## Preferences
elasmu = 1.45                                       # Elasticity of marginal utility of consumption
prstp = .015                                        # Initial rate of social time preference per year


alphab0 = 0.80                                      # Initial share of Brown good in the consumption composite (https://www.iea.org/reports/renewable-electricity)
alphag0 = 1 - alphab0                               # Initial share of Green good in the consumption composite

## Population and technology

dk = 0.100                                          # Depreciation rate on capital (per year)
q0 = 105.5                                          # Initial world gross output 2015 (trill 2010 USD)
k0 = 223                                            # Initial capital value 2015 (trill 2010 USD)

ab0 = 5                                             # Initial level of capital productivity in Brown sector
ag0 = 5                                             # Initial level of capital productivity in Green sector

bb0 = 5                                             # Initial level of labor productivity in Brown sector
bg0 = 5                                             # Initial level of labor productivity in Green sector

## Carbon cycle
# Initial Conditions
mat0 = 851                                          # Initial Concentration in atmosphere 2015 (GtC)
mu0 = 460                                           # Initial Concentration in upper strata 2015 (GtC)
ml0 = 1740                                          # Initial Concentration in lower strata 2015 (GtC)
mateq = 588                                         # Equilibrium concentration atmosphere  (GtC)
mueq = 360                                          # Equilibrium concentration in upper strata (GtC)
mleq = 1720                                         # Equilibrium concentration in lower strata (GtC)
# Flow paramaters
b12 = .12                                           # Carbon cycle transition matrix
b23 = 0.007                                         # Carbon cycle transition matrix
## Climate model parameters
t2xco2 = 3.1                                        # Equilibrium temp impact (oC per doubling CO2)
fex0 = 0.5                                          # 2015 forcings of non-CO2 GHG (Wm-2)
fex1 = 1.0                                          # 2100 forcings of non-CO2 GHG (Wm-2)
tocean0 = .0068                                     # Initial lower stratum temp change (C from 1900)
tatm0 = 0.85                                        # Initial atmospheric temp change (C from 1900)
c1 = 0.1005                                         # Climate equation coefficient for upper level
c3 = 0.088                                          # Transfer coefficient upper to lower stratum
c4 = 0.025                                          # Transfer coefficient for lower level 
fco22x = 3.6813                                     # Forcings of equilibrium CO2 doubling (Wm-2)
## Climate damage parameters
ab2 = 0.00236                                       # Damage quadratic term in Brown sector
ag2 = 0.00236                                       # Damage quadratic term in Green sector

ab3 = 2.00                                          # Damage exponent in Brown sector
ag3 = 2.00                                          # Damage exponent in Green sector
## Abatement cost
# Learning by doing
gammabar = 2.5
gamma0 = 0.1
nu =  1.5#0.85

# Program control variables

# PARAMETERS
# population
pop0 = 7403                                         # Initial world population 2015 (millions)
popadj = 0.134                                      # Growth rate to calibrate to 2050 pop projection
popasym = 11500                                     # Asymptotic population (millions)
l = [pop0] * T                                      # Level of population and labor
lr = [0] * T                                        # Population growth rate

lb0 = alphab0 * pop0                                # Initial capital in Brown sector
lg0 = alphag0 * pop0                                # Initial capital in Green sector

kb0 = alphab0 * k0                                  # Initial capital in Brown sector
kg0 = alphag0 * k0                                  # Initial capital in Green sector

sigmag = [0] * T                                    # CO2-equivalent-emissions output ratio
sigmab = [0] * T                                    # CO2-equivalent-emissions output ratio

gsigg = [0] * T                                     # Change in sigma (cumulative improvement of energy efficiency) in Green sector
gsigb = [0] * T                                     # Change in sigma (cumulative improvement of energy efficiency) in brown sector

rr = [1] * T                                        # Average utility social discount rate

forcoth = [fex0] * T                                # Exogenous forcing for other greenhouse gases

etree = [0] * T                                     # Emissions from deforestation
cumetree = [0] * T                                  # Cumulative from land

# Parameters for long-run consistency of carbon cycle
b11 = 1 - b12
b21 = b12 * mateq/mueq
b22 = 1 - b21 - b23
b32 = b23 * mueq/mleq
b33 = 1 - b32

## Emissions parameters
eland0 = 2.6                                        # Carbon emissions from land 2015 (GtCO2 per year)
deland = .115                                       # Decline rate of land emissions (per period)
e0 = 35.85                                          # Industrial emissions 2015 (GtCO2 per year) 
miu0 = .03                                          # Initial emissions control rate for base case 2015

gsigmab1 = -0.0152                                  # Initial growth of sigma (per year)
gsigmag1 = -0.0152                                  # Initial growth of sigma (per year)

dsigb = -0.001                                      # Decline rate of decarbonization (per period) *5
dsigg = -0.001                                      # Decline rate of decarbonization (per period)

sigg0 = (e0/q0) * alphag0                           # Initial emission intensity of output in Green sector
sigb0 = (e0/q0)                                     # Initial emission intensity of output in Brown sector

gsigb[0] = gsigmab1                                 # Initial growth rate of emission intensity of output in Brown sector
gsigg[0] = gsigmag1                                 # Initial  growth rate of emission intensity of output in Green sector

sigmab[0] = sigb0                                   # Initial emission intensity of output in Brown sector
sigmag[0] = sigg0                                   # Initial emission intensity of output in Green sector

etree[0] = eland0                                   # Exogenous emissions from land-use change
cumetree[0] = 100                                   # Initial levbel of exogenous emissions from land-use change

for t in range(T - 1):
    l[t + 1] = l[t] * (popasym/l[t])**popadj
    lr[t] = l[t + 1]/l[t]
    
    gsigg[t + 1] = gsigg[t] * (1 + dsigg)**tstep
    sigmag[t + 1] = sigmag[t] * np.exp(gsigg[t] * tstep)
    
    gsigb[t + 1] = gsigb[t] * (1 + dsigb)**tstep
    sigmab[t + 1] = sigmab[t] * np.exp(gsigb[t] * tstep)
    
    etree[t + 1] = eland0 * (1 - deland)**(t + 1)
    cumetree[t + 1] = cumetree[t] + etree[t] * (5/3.666)

    rr[t + 1] = 1/((1 + prstp)**(tstep*(t + 1)))
    if t<17:
        forcoth[t + 1] = fex0 + (1/17) * (fex1 - fex0) * (t + 1)
    else:
        forcoth[t + 1] = fex1

# ============================================ DICE Variables ============================================ #

#  State Variables 
GET_Agt = [ag0] * T                        # Level of Capital productivity in Green sector
GET_Abt = [ab0] * T                        # Level of Capital productivity in Brown sector

GET_Bgt = [bg0] * T                        # Level of Labor productivity in Green sector
GET_Bbt = [bb0] * T                        # Level of Labor productivity in Brown sector

GET_Kgt = [kg0]*T                          # Capital ($trill, 2005$) in Green sector
GET_Kbt = [kb0]*T                          # Capital ($trill, 2005$) in Brown sector

GET_Lgt = [lg0]*T                          # Labor in Green sector (millions)
GET_Lbt = [lb0]*T                          # Labor in Brown sector (millions)

GET_MAT = [mat0]*T                         # Atmospheric concentration of carbon (GTC)
GET_ML = [ml0]*T                           # Concentration in biosphere and upper oceans (GTC)
GET_MU = [mu0]*T                           # Concentration in deep oceans (GTC)
GET_TATm = [tatm0]*T                       # Atmospheric temperature (degrees Celsius above preindustrial)
GET_TOCEAN = [tocean0]*T                   # Lower ocean temperature (degrees Celsius above preindustrial)
GET_FORC = [0]*T                           # Total increase in radiative forcing since preindustrial (Watts per square meter)
GET_CCA = [400]*T                          # Cumulative industrial carbon emissions (GTC)

 # == Other variables == #
GET_Yt = [0]*T                             # Output net of abatement cost and climate damage ($trill)
GET_Ggt = [0]*T                            # Gross output in Green sector
GET_Gbt = [0]*T                            # Gross output in Brown sector
GET_Ygt = [0]*T                            # Net output in Green sector
GET_Ybt = [0]*T                            # Net output in Brown sector

GET_CTLbgt = [0]*T                         # Labor transition cost from Brown to Green (trill 2010 USD)
GET_CTLgbt = [0]*T                         # Labor transition cost from Green to Brown (trill 2010 USD)
GET_CTKbgt = [0]*T                         # Capital transition cost from Brown to Green (trill 2010 USD)
GET_CTKgbt = [0]*T                         # Capital transition cost from Green to Brown (trill 2010 USD)

GET_Et = [0]*T                             # Total carbon emissions (GTCO2 per year)
GET_Egt = [0]*T                            # Carbon emissions in Green sector (GTCO2 per year)
GET_Ebt = [0]*T                            # Carbon emissions in Brown sector (GTCO2 per year)

GET_Eind = [0]*T                           # Industrial emissions (GTCO2 per year)
GET_Ccatot = [0]*T                         # Industrial emissions (GTCO2 per year)

GET_ABTFRAcb = [0]*T                       # Abatement cost in Brown sector (fraction of gross output)
GET_ABTFRAcg = [0]*T                       # Abatement cost in Green sector (fraction of gross output)

GET_ABTCOSb = [0]*T                        # Abatement cost in Brown sector ($ trillion)
GET_ABTCOSg = [0]*T                        # Abatement cost in Green sector ($ trillion)

GET_DAMFRACb = [0]*T                       # Total damage in Brown sector (fraction of gross output)
GET_DAMFRACg = [0]*T                       # Total damage in Green sector (fraction of gross output)

GET_DAMAGEb = [0]*T                        # Total damage in Brown sector ($ trillion)
GET_DAMAGEg = [0]*T                        # Total damage in Green sector ($ trillion)

GET_Igt = [0]*T                            # Saving in Green sector ($trill per year)
GET_Ibt = [0]*T                            # Saving in Brown sector ($trill per year)

GET_Ct = [0]*T                             # Consumption ($trill per year)
GET_Cgt = [0]*T                            # Consumption in Green sector ($trill per year)
GET_Cbt = [0]*T                            # Consumption in Brown sector ($trill per year)

GET_CPC = [0]*T                            # Consumption per capita ($thous per year)
GET_CEMUTOTPER = [0]*T                     # Utility of p. c. consumption
GET_PERIODU = [0]*T                        # Total period utility
GET_CPRICE = [0]*T                         # Carbon price (2005$ per ton of CO2)

GET_TKbgt = [0]*T                          # Capital transfer ($trill, 2005$) from Brown to Green sector
GET_TKgbt = [0]*T                          # Capital transfer ($trill, 2005$) from Green to Brown sector

GET_TLbgt = [0]*T                          # Labor transfer (millions) from Brown to Green sector
GET_TLgbt = [0]*T                          # Labor transfer (millions) from Green to Brown sector

GET_srg = [0]*T                            # Optimal saving rate in Green sector
GET_srb = [0]*T                            # Optimal saving rate in Brown sector

GET_ktr = [0]*T                            # Optimal share of capital transfer between sectors 
GET_ltr = [0]*T                            # Optimal share of labor transfer between sectors

GET_lgb = [0]*T                            # Optimal share of Green labor transferred to Brown sector
GET_lbg = [0]*T                            # Optimal share of Brown labor transferred to Green sector

GET_rd = [0]*T                             # Optimal level of R&D allocated to improving productivity
GET_rdk = [0]*T                            # Optimal share of R&D allocated to improving capital productivity
GET_rdl = [0]*T                            # Optimal share of R&D allocated to improving labor productivity

GET_kgb = [0]*T                            # Optimal share of Green capital transferred to Brown sector
GET_kbg = [0]*T                            # Optimal share of Brown capital transferred to Green sector

GET_rb = [0]*T                             # Optimal R&D share in Brown sector
GET_rg = [0]*T                             # Optimal R&D share in Green sector

GET_eb = [0]*T                             # Optimal Education share in Brown sector
GET_eg = [0]*T                             # Optimal Education share in Green sector

GET_miub = [0]*T                           # Optimal Abatement rate in Brown sector
GET_miug = [0]*T                           # Optimal Abatement rate in Green sector

GET_gamma = [0]*T                          # Abatement cost coefficient
GET_alphab = [0]*T                         # Brown sector share
GET_alphag = [0]*T                         # Green sector share
# ============================================ Model Parameteres ============================================ #
thetab = 0.30                              # Capital share in Brown sector 
thetag = 0.30                              # Capital share in Green sector 

epsilon = 0.5                              # Elasticity of substitution between Brown and Green outputs

wrb0 = 0.5/alphab0                         # Coefficient of capital productivity growth in Brown sector 
wrg0 = 0.5/alphab0                         # Coefficient of capital productivity growth in Green sector 

web0 = 0.5/alphab0                         # Coefficient of human productivity growth in Brown sector 
weg0 = 0.5/alphab0                         # Coefficient of human productivity growth in Green sector 

pab = 0.5                           # Power of capital productivity growth in Brown sector 
pag = 0.5                           # Power of capital productivity growth in Green sector 

pbb = 0.2                           # Power of human productivity growth in Brown sector 
pbg = 0.2                           # Power of human productivity growth in Green sector 

rdshare = 0.02                             # Share of labour and capital devoted to R&D/Education

scale1 = 0.0302455265681763                # Multiplicative scaling coefficient 
scale2 = -10993.704                        # Additive scaling coefficient

tclabg_l = 0 #1e-2                         # Linear transfer cost parameters
tclabb_l = 0 #1e-2
tccapg_l = 0 #1e-1
tccapb_l = 0 #1e-1
tclabg_q = 1e-3                            # Quadratic transfer cost parameters
tclabb_q = 1e-3
tccapg_q = 1e-3
tccapb_q = 1e-3

alphagLin = [0.20 + (0.9 - 0.2)/85 * 5 * i for i in range(18)] + 82*[0.90]
alphagDel = 7*[0.20] + [0.20 + (0.9 - 0.2)/50 * 5 * i for i in range(11)] + 82*[0.90]
alphagFst = [0.20 + (0.9 - 0.2)/35 * 5 * i for i in range(8)] + 10*[0.90]  + 82*[0.9]
# ========================================= Transition Function ========================================= #

def state( ST1, EX1, Act1 ):
    
    [Ab1, Ag1, Bb1, Bg1, Kb1, Kg1, Lb1, Lg1, MAT1, MU1, ML1, TATM1, TOCEAN1] = ST1
    [sigb1, sigg1, Etree1, pop1, popr1, rr1, forcoth2, i1] = EX1
    [mub1, mug1, ktbg1, ktgb1, ltbg1, ltgb1, rdk1, rb1, eb1, srb1, srg1] = Act1

    # Green nudge
    if cs == 0:
        # Baseline case
        alphag1 = alphag0

    if cs == 1:
        # Linear case
        alphag1 = alphagLin[i1]
    
    if cs == 2:
        # Delayed case
        alphag1 = alphagDel[i1]
    
    if cs == 3:
        # Fast case
        alphag1 = alphagFst[i1]

    alphab1 = 1 - alphag1

    Lby1 = (1 - rdshare) * Lb1
    Lgy1 = (1 - rdshare) * Lg1

    # Labor transfer between two sectors     
    TLbg1 = ltbg1 * Lby1
    TLgb1 = ltgb1 * Lgy1
    
    CTLbg1 = TLbg1 * tclabg_l + tclabg_q * TLbg1**2
    CTLgb1 = TLgb1 * tclabb_l + tclabb_q * TLgb1**2
    
    # Labor allocation in two sectors
    Lbx1 = Lby1 - TLbg1 + TLgb1
    Lgx1 = Lgy1 - TLgb1 + TLbg1

    Kby1 = (1 - rdshare) * Kb1
    Kgy1 = (1 - rdshare) * Kg1
    
    # Capital transfer between two sectors       
    TKbg1 = ktbg1 * Kby1
    TKgb1 = ktgb1 * Kgy1
    
    CTKbg1 = TKbg1 * tccapg_l + tccapg_q * TKbg1**2
    CTKgb1 = TKgb1 * tccapb_l + tccapb_q * TKgb1**2
    
    # Capital allocation in two sectors
    Kbx1 = Kby1 - TKbg1 + TKgb1
    Kgx1 = Kgy1 - TKgb1 + TKbg1
    
    ################# R&D #####################
    # Labor and Capital input to R&D sector # https://www.oecd-ilibrary.org/docserver/sti_scoreboard-2015-10-en.pdf?expires=1676648756&id=id&accname=guest&checksum=2207CCF67522ECEED1A1C960623BF9E9
    Lrd_b = rdshare * Lb1
    Lrd_g = rdshare * Lg1
    
    Krd_b = rdshare * Kb1
    Krd_g = rdshare * Kg1

    RDin = (Ab1 * Krd_b)**thetab * (Bb1 * Lrd_b)**(1 - thetab) + (Ag1 * Krd_g)**thetag * (Bg1 * Lrd_g)**(1 - thetag)
    RDink = rdk1 * RDin
    rdl1 = 1 - rdk1
    RDinl = rdl1 * RDin

    RDinlb = eb1 * RDinl
    eg1 = 1 - eb1
    RDinlg=eg1 * RDinl

    RDinkb = rb1 * RDink
    rg1 = 1 - rb1
    RDinkg=rg1 * RDink

    # Labor productivity in two sectors
    
    Bb2 = (1 + web0 * RDinlb/RDin)**pbb * Bb1
    Bg2 = (1 + weg0 * RDinlg/RDin)**pbg * Bg1

    # Capital productivity in two sectors
    Ab2 = (1 + wrb0 * RDinkb/RDin)**pab * Ab1
    Ag2 = (1 + wrg0 * RDinkg/RDin)**pag * Ag1

    # Gross output of two sectors
    YGROSSb1 = (Ab1 * Kby1)**thetab * (Bb1 * Lby1/1000)**(1 - thetab)
    YGROSSg1 = (Ag1 * Kgy1)**thetag * (Bg1 * Lgy1/1000)**(1 - thetag)
    
    # Gross emissions from two sectors
    Eb1 = sigb1 * YGROSSb1
    Eg1 = sigg1 * YGROSSg1
    
    # Net emissions from two sectors after abatement 
    EIND1 = Eb1  * (1 - mub1) + Eg1 * (1 - mug1)
    E1 = EIND1 + Etree1

    # Share of damages in two sectors    
    DAMFRACb1 = ab2 * TATM1**ab3
    DAMFRACg1 = ag2 * TATM1**ag3
    
    # Total damages in two sectors 
    DAMAGEb1 = DAMFRACb1 * YGROSSb1
    DAMAGEg1 = DAMFRACg1 * YGROSSg1
    
    # Abatement cost
    gamma1 = gamma0 * ((Ab1 + Ag1 + Bb1 + Bg1)/(ab0 + ag0 + bb0 + bg0))**(-nu) 
    ABTFRACb1 = gamma1 * mub1**gammabar
    ABTFRACg1 = gamma1 * mug1**gammabar
    
    ABTCOSTb1 = ABTFRACb1 * YGROSSb1
    ABTCOSTg1 = ABTFRACg1 * YGROSSg1
    
    # Net output in two sectors             
    Yb1 = max(0, YGROSSb1 - DAMAGEb1 - CTLgb1 - CTKgb1 - ABTCOSTb1)	
    Yg1 = max(0, YGROSSg1 - DAMAGEg1 - CTLbg1 - CTKbg1 - ABTCOSTg1)
    Y1 = Yb1 + Yg1
    
    # Savings in two sectors
    Ib1 = srb1 * Yb1
    Ig1 = srg1 * Yg1
    
    # Consumption in two sectors
    Cb1 = Yb1 - Ib1
    Cg1 = Yg1 - Ig1
    
    # Composite consumption
    C1 = (alphab1 * Cb1**((epsilon - 1)/epsilon) + (1 - alphab1) * Cg1**((epsilon - 1)/epsilon))**(epsilon/(epsilon - 1))

    #C1 = Cb1+Cg1

    # Consumption per capita
    CPC1 = (C1/pop1) * 1000
    
    # Individual utility of consumption
    PERIODU1 = (CPC1**(1 - elasmu) - 1)/ (1 - elasmu) - 1
    
    # Total utility of consumption
    CEMUTOTPER1 = PERIODU1 * pop1 * rr1
    
    # CO2 concentrations
    MAT2 = b11 * MAT1 + b21 * MU1 + E1 * tstep/3.666
    MU2 = b12 * MAT1 + b22 * MU1 + b32 * ML1
    ML2 = b23 * MU1 + b33 * ML1
    
    # Radiative forcing
    FORC2 = fco22x * (np.log(MAT2/mateq))/np.log(2) + forcoth2
    
    # Temperatures
    TATM2 = TATM1 + c1 * (FORC2 - (fco22x/t2xco2) * TATM1 - c3 * (TATM1 - TOCEAN1))
    TOCEAN2 = TOCEAN1 + c4 * (TATM1 - TOCEAN1)
    
    ################# Next-period Capital #####################
    # Capital formation in two sectors
    Kb2 = (1 - dk)**tstep * (Kbx1 + Krd_b) + tstep * Ib1
    Kg2 = (1 - dk)**tstep * (Kgx1 + Krd_g) + tstep * Ig1
    
    ################# Next-period Labor #####################
    # Labor formation in two sectors
    Lb2 = (Lbx1 + Lrd_b) * popr1
    Lg2 = (Lgx1 + Lrd_g) * popr1

    ST2 = [Ab2, Ag2, Bb2, Bg2, Kb2, Kg2, Lb2, Lg2, MAT2, MU2, ML2, TATM2, TOCEAN2]
    LEV2 = [Y1, Yb1, Yg1, YGROSSb1, YGROSSg1, E1, Eb1, Eg1, EIND1, DAMFRACb1, DAMFRACg1, DAMAGEb1, DAMAGEg1, ABTFRACb1, ABTFRACg1, ABTCOSTb1, ABTCOSTg1, Ib1, Ig1, C1, Cb1, Cg1, CPC1, CEMUTOTPER1, PERIODU1, FORC2, gamma1, alphab1, RDin]
    TSF2 = [TLbg1, TLgb1, TKbg1, TKgb1, CTLbg1, CTLgb1, CTKbg1, CTKgb1]
    return ( ST2, LEV2, TSF2 )

# =========================================== Welfare Function (DICE model) =========================================== #
def fDICE(v):  
    W = scale2
    for i in range(T-1):
        if i == 0:
            GET_FORC[i] = fco22x * (np.log(GET_MAT[i]/mateq))/np.log(2) + forcoth[i]
        
        STi = [GET_Abt[i], GET_Agt[i], GET_Bbt[i], GET_Bgt[i], GET_Kbt[i], GET_Kgt[i], GET_Lbt[i], GET_Lgt[i], GET_MAT[i], GET_MU[i], GET_ML[i], GET_TATm[i], GET_TOCEAN[i]]
        EXi = [sigmab[i], sigmag[i], etree[i], l[i], lr[i], rr[i], forcoth[i + 1], i]
        Acti = [v[i], v[T + i], v[2 * T + i], v[3 * T + i], v[4 * T + i], v[5 * T + i], v[6 * T + i], v[7 * T + i], v[8 * T + i], v[9 * T + i], v[10 * T + i]]

        ( STii, LEVii, TSFii ) = state( STi, EXi, Acti )
        [GET_Abt[i + 1], GET_Agt[i + 1], GET_Bbt[i + 1], GET_Bgt[i + 1], GET_Kbt[i + 1], GET_Kgt[i + 1], GET_Lbt[i + 1], GET_Lgt[i + 1], GET_MAT[i + 1], GET_MU[i + 1], GET_ML[i + 1], GET_TATm[i + 1], GET_TOCEAN[i + 1]] = STii
        [GET_Yt[i], GET_Ybt[i], GET_Ygt[i], GET_Gbt[i], GET_Ggt[i], GET_Et[i], GET_Ebt[i], GET_Egt[i], GET_Eind[i], GET_DAMFRACb[i], GET_DAMFRACg[i], GET_DAMAGEb[i], GET_DAMAGEg[i], GET_ABTFRAcb[i], GET_ABTFRAcg[i], GET_ABTCOSb[i], GET_ABTCOSg[i], GET_Ibt[i], GET_Igt[i], GET_Ct[i], GET_Cbt[i], GET_Cgt[i], GET_CPC[i], GET_CEMUTOTPER[i], GET_PERIODU[i], GET_FORC[i + 1], GET_gamma[i], GET_alphab[i], GET_rd[i]] = LEVii
        [GET_TLbgt[i], GET_TLgbt[i], GET_TKbgt[i], GET_TKgbt[i], GET_CTLbgt[i], GET_CTLgbt[i], GET_CTKbgt[i], GET_CTKgbt[i]] = TSFii
        W = W +  tstep * scale1 * GET_CEMUTOTPER[i]
    return -W
    
# ======================================== Optimization Algorithm (DICE model) ======================================== #

ndv = 11 #number of decision variables
if os.path.exists('GETResult.npy'):
    x0=np.load('GETResult.npy')
else:
    x0 = ndv * T * [0]
    x0[(T*0):(T*1)]=[0.04] * T
    x0[(T*1):(T*2)]=[0.03] * T
    x0[(T*2):(T*3)]=[0] * T
    x0[(T*3):(T*4)]=[0] * T
    x0[(T*4):(T*5)]=[0] * T
    x0[(T*5):(T*6)]=[0] * T
    x0[(T*6):(T*7)]=[0.5] * T
    x0[(T*7):(T*8)]=[0.5] * T
    x0[(T*8):(T*9)]=[0.25] * T
    x0[(T*9):(T*10)]=[0.25] * T
    x0[(T*10):(T*11)]=[0.25] * T

# == bounds == #
lb = ndv * T * [0.00]
ub = ndv * T * [0.99]

limmiu = 1.2                                        # Upper limit on abatement control rate after 2150
ub[29:T] = (T - 29) * [limmiu]
ub[T + 29:2 * T] = (T - 29) * [limmiu]

# 2015 initial values

# Abatement rate
lb[0] = 0.03
lb[T] = 0.03
ub[0] = 0.031
ub[T] = 0.031

# Transfer rate
lb[2 * T] = 0.00
lb[3 * T] = 0.00
lb[4 * T] = 0.00
lb[5 * T] = 0.00
ub[2 * T] = 0.001
ub[3 * T] = 0.001
ub[4 * T] = 0.001
ub[5 * T] = 0.001

# Capital R&D rate
lb[6 * T] = 0.50
ub[6 * T] = 0.501

# Brown capital R&D rate
lb[7 * T] = 0.50
ub[7 * T] = 0.501

# Brown labour R&D rate
lb[8 * T] = 0.50
ub[8 * T] = 0.501

# saving rate
lb[9 * T] = 0.26
lb[10 * T] = 0.26
ub[9 * T] = 0.261
ub[10 * T] = 0.261

GET_bnds = Bounds(lb, ub)

# == optimization == # COBYLA
ftol = 1e-15
eps = 1e-7
maxiter = 10000
    
res = minimize(fDICE, x0, method='SLSQP', bounds=GET_bnds, options={'ftol': ftol, 'eps': eps, 'disp': True, 'maxiter': maxiter})
resDICE = res.x

# np.save('GETResult.npy',resDICE)

GET_miub = [resDICE[b] for b in range(0 * T, 1 * T)]
GET_miug = [resDICE[b] for b in range(1 * T, 2 * T)]

GET_kbg = [resDICE[b] for b in range(2 * T, 3 * T)]
GET_kgb = [resDICE[b] for b in range(3 * T, 4 * T)]

GET_lbg = [resDICE[b] for b in range(4 * T, 5 * T)]
GET_lgb = [resDICE[b] for b in range(5 * T, 6 * T)]

GET_rdk = [resDICE[b] for b in range(6 * T, 7 * T)]
GET_rdl = [(1 - GET_rdk[b]) for b in range(T)]

GET_rb = [resDICE[b] for b in range(7 * T, 8 * T)]
GET_rg = [(1 - GET_rb[b]) for b in range(T)]

GET_eb = [resDICE[b] for b in range(8 * T, 9 * T)]
GET_eg = [(1 - GET_eb[b]) for b in range(T)]

GET_srb = [resDICE[b] for b in range(9 * T, 10 * T)]
GET_srg = [resDICE[b] for b in range(10 * T, 11 * T)]

GET_alphag = [(1 - GET_alphab[b]) for b in range(T)]

GET_W = -res.fun

# ======================================== GAMS Results ======================================== #
GET_GAMS_miu = [0.03, 0.187151008, 0.211465058, 0.237698399, 0.265906995, 0.29615028, 0.328490392, 0.362991483, 0.399719103, 0.438739652, 0.480119893, 0.523926564, 0.570226125, 0.619084752, 0.67056881, 0.72474622, 0.781689586, 0.841482679, 0.904248291, 0.970126461, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 0]
GET_GAMS_sav = [0.260591864, 0.257176389, 0.254454369, 0.252145939, 0.250203464, 0.248584521, 0.247251684, 0.246172104, 0.245316979, 0.244660974, 0.244181578, 0.243858369, 0.24367205, 0.243603073, 0.243629419, 0.243722801, 0.243841914, 0.24392016, 0.243868394, 0.243497055, 0.243348473, 0.243504644, 0.243653609, 0.24381924, 0.244028185, 0.244319811, 0.244760999, 0.245470901, 0.246663957, 0.245644746, 0.245817104, 0.245979039, 0.246133054, 0.246280835, 0.246423518, 0.246561869, 0.246696393, 0.246827424, 0.246955175, 0.247079784, 0.24720134, 0.247319903, 0.247435517, 0.247548221, 0.247658052, 0.247765048, 0.247869252, 0.24797071, 0.248069476, 0.248165602, 0.248259149, 0.248350177, 0.24843875, 0.248524932, 0.248608787, 0.248690379, 0.248769773, 0.248847031, 0.248922215, 0.248995383, 0.249066593, 0.249135901, 0.249203359, 0.249269019, 0.249332928, 0.249395133, 0.249455676, 0.249514598, 0.249571935, 0.249627718, 0.249681974, 0.249734719, 0.249785952, 0.249835649, 0.249883738, 0.249930073, 0.249974373, 0.250016129, 0.250054425, 0.250087643, 0.250112945, 0.250125356, 0.25011619, 0.2500703, 0.249961274, 0.249743074, 0.249335355, 0.248597717, 0.24728373, 0.244957948, 0.258278146, 0.258278146, 0.258278146, 0.258278146, 0.258278146, 0.258278146, 0.258278146, 0.258278146, 0.258278146, 0.258278146]
GET_GAMS_TATm = [0.85, 1.016341648, 1.184309383, 1.353597424, 1.523751712, 1.694192862, 1.86424317, 2.033153291, 2.200125802, 2.364333811, 2.524932905, 2.681325332, 2.832682801, 2.978148669, 3.116871822, 3.248021193, 3.37078899, 3.484386887, 3.585082054, 3.672474775, 3.746100953, 3.807907828, 3.860131818, 3.904492253, 3.942329556, 3.974703608, 4.002464023, 4.026300661, 4.046780244, 4.064373184, 4.070147237, 4.066682072, 4.055903972, 4.03926852, 4.017890823, 3.99263807, 3.964195353, 3.933112571, 3.899838003, 3.864742498, 3.828137029, 3.790285558, 3.751414565, 3.711720192, 3.671373665, 3.630525468, 3.589308602, 3.547841162, 3.506228401, 3.464564415, 3.422933515, 3.381411379, 3.340066014, 3.298958567, 3.258144023, 3.217671808, 3.177586299, 3.137927285, 3.098730359, 3.06002727, 3.021846234, 2.984212211, 2.947147155, 2.910670237, 2.874798052, 2.839544803, 2.804922472, 2.770940975, 2.737608312, 2.704930697, 2.672912688, 2.6415573, 2.61086612, 2.58083941, 2.551476202, 2.522774394, 2.494730839, 2.467341424, 2.440601152, 2.414504217, 2.389044072, 2.364213501, 2.34000468, 2.316409238, 2.293418314, 2.271022616, 2.249212463, 2.227977845, 2.20730846, 2.187193759, 2.16762299, 2.148585233, 2.130069436, 2.112064445, 2.094559044, 2.077541972, 2.061001958, 2.044927742, 2.0293081, 2.026026675]
GET_GAMS_Yt = [104.9972283, 125.0094035, 147.2419965, 171.7043232, 198.4314979, 227.4523143, 258.7882001, 292.4498781, 328.4294038, 366.6826898, 407.1501469, 454.3958325, 503.5206083, 554.690929, 608.0341037, 663.6508544, 721.6243151, 782.0269156, 844.9686825, 910.5293712, 980.3476149, 1053.818427, 1130.541323, 1210.485198, 1293.619898, 1379.914571, 1469.336369, 1561.849466, 1657.414325, 1750.924022, 1852.267048, 1957.103686, 2065.325039, 2176.82639, 2291.505827, 2409.262978, 2529.997926, 2653.610317, 2779.998689, 2909.059973, 3040.689177, 3174.779202, 3311.220787, 3449.902538, 3590.711049, 3733.531069, 3878.245724, 4024.736771, 4172.884866, 4322.569861, 4473.6711, 4626.067717, 4779.63894, 4934.264382, 5089.824326, 5246.200001, 5403.273845, 5560.929755, 5719.053319, 5877.532037, 6036.255522, 6195.11569, 6354.006927, 6512.826246, 6671.473421, 6829.851113, 6987.864973, 7145.423735, 7302.439293, 7458.82676, 7614.504521, 7769.394268, 7923.421022, 8076.513148, 8228.602357, 8379.623696, 8529.515535, 8678.219533, 8825.680613, 8971.846916, 9116.669756, 9260.103561, 9402.105824, 9542.637029, 9681.660591, 9819.142781, 9955.052653, 10089.36197, 10222.04511, 10353.07901, 10482.44307, 10610.11906, 10736.09105, 10860.34534, 10982.87033, 11103.65649, 11222.69624, 11339.98388, 11455.82994, 0]
GET_GAMS_Et = [33.38173132, 35.40110156, 37.14478032, 38.56218576, 39.61164224, 40.25983713, 40.48119137, 40.25709617, 39.57518201, 38.42867657, 36.81570289, 34.73872078, 32.20393473, 29.22083043, 25.80163672, 21.96089214, 17.71467689, 13.08004082, 8.072705725, 2.708813465, 1.047980144, 1.024727397, 1.003206411, 0.983060392, 0.96401018, 0.945839289, 0.928386915, 0.911548508, 0.895285822, -16.01341323, -15.86447446, -15.70070902, -15.52262038, -15.33091984, -15.12645612, -14.91016418, -14.68302855, -14.44605658, -14.20025793, -13.94662999, -13.68614793, -13.41975849, -13.14837154, -12.87285974, -12.59405278, -12.31273882, -12.02966317, -11.74552438, -11.46098243, -11.17664854, -10.8930945, -10.61084896, -10.33040021, -10.05220032, -9.776658725, -9.504147532, -9.235011097, -8.9695541, -8.70804864, -8.450740247, -8.197841686, -7.949539776, -7.705997277, -7.46734572, -7.233708486, -7.00517103, -6.78180774, -6.563674275, -6.350805454, -6.143224866, -5.940935607, -5.743934033, -5.552197474, -5.365702793, -5.184399649, -5.008244578, -4.837179116, -4.671140279, -4.510050506, -4.353830404, -4.202398158, -4.05566268, -3.913530776, -3.775898053, -3.642657679, -3.513681178, -3.388829669, -3.267936891, -3.150787237, -3.037088447, -2.926414695, -2.818128611, -2.711189509, -2.603888513, -2.493219728, -2.373613039, -2.233595441, 2.832914097, 4.598556911, 0]

T0 = T - 12
Xaxis = range(2015, 2015 + tstep * T0, tstep)

plt.plot(Xaxis[0:T0], GET_miub[0:T0], '--y', label = "Brown")
plt.plot(Xaxis[0:T0], GET_miug[0:T0], '--g', label = "Green")
plt.plot(Xaxis[0:T0], GET_GAMS_miu[0:T0], 'k', label = "DICE")
plt.xlabel('Time (years)')
plt.ylabel('Abatement Rate')
plt.title('Optimal Abatement Rate')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 1.5, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_kgb[0:T0], '--y', label = "Green -> Brown")
plt.plot(Xaxis[0:T0], GET_kbg[0:T0], '--g', label = "Brown -> Green")
plt.xlabel('Time (years)')
plt.ylabel('Share')
plt.title('Capital Transfer')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_Kbt[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_Kgt[0:T0], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('$trill, 2005')
plt.title('Capital')
plt.legend(loc=1, prop={'size':8})
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_lgb[0:T0], '--y', label = "Green -> Brown")
plt.plot(Xaxis[0:T0], GET_lbg[0:T0], '--g', label = "Brown -> Green")
plt.xlabel('Time (years)')
plt.ylabel('Share')
plt.title('Labor Transfer')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_Lbt[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_Lgt[0:T0], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('million people')
plt.title('Labor')
plt.legend(loc=1, prop={'size':8})
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], [i * j for i, j in zip(GET_rdl[0:T0], (GET_eb[0:T0]))], '--y', label = "R&D Share in Labor Productivity of Brown sector")
plt.plot(Xaxis[0:T0], [i * j for i, j in zip(GET_rdl[0:T0], (GET_eg[0:T0]))], '--g', label = "R&D Share in Labor Productivity of Green sector")
plt.plot(Xaxis[0:T0], [i * j for i, j in zip(GET_rdk[0:T0], (GET_rb[0:T0]))], '-y', label = "R&D Share in Capital Productivity of Brown sector")
plt.plot(Xaxis[0:T0], [i * j for i, j in zip(GET_rdk[0:T0], (GET_rg[0:T0]))], '-g', label = "R&D Share in Capital Productivity of Green sector")
plt.xlabel('Time (years)')
plt.ylabel('R&D Share')
plt.title('R&D Investment')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

# plt.plot(Xaxis[0:T0], GET_rdl[0:T0], '--m', label = "Labor")
# plt.plot(Xaxis[0:T0], GET_rdk[0:T0], '--c', label = "capital")
# plt.xlabel('Time (years)')
# plt.ylabel('Share')
# plt.title('R&D Investment Share')
# plt.legend(loc=1, prop={'size':8})
# plt.ylim(ymax = 1.0, ymin = 0)
# plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
# plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
# plt.show()

# plt.plot(Xaxis[0:T0], GET_rb[0:T0], '--y', label = "Brown sector")
# plt.plot(Xaxis[0:T0], GET_rg[0:T0], '--g', label = "Green sector")
# plt.xlabel('Time (years)')
# plt.ylabel('Share')
# plt.title('R&D Investment Share in Capital Productivit')
# plt.legend(loc=1, prop={'size':8})
# plt.ylim(ymax = 1.0, ymin = 0)
# plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
# plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
# plt.show()

# plt.plot(Xaxis[0:T0], GET_eb[0:T0], '--y', label = "Brown sector")
# plt.plot(Xaxis[0:T0], GET_eg[0:T0], '--g', label = "Green sector")
# plt.xlabel('Time (years)')
# plt.ylabel('Share')
# plt.title('Education Investment Share in Labor Productivity')
# plt.legend(loc=1, prop={'size':8})
# plt.ylim(ymax = 1.0, ymin = 0)
# plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
# plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
# plt.show()

plt.plot(Xaxis[0:T0], GET_srb[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_srg[0:T0], '--g', label = "Green sector")
plt.plot(Xaxis[0:T0], GET_GAMS_sav[0:T0], 'k', label = "DICE")
plt.xlabel('Time (years)')
plt.ylabel('Rate')
plt.title('Saving Rate')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_TATm[0:T0], 'r', label = "GET")
plt.plot(Xaxis[0:T0], GET_GAMS_TATm[0:T0], 'k', label = "DICE")
plt.xlabel('Time (years)')
plt.ylabel('Atmospheric Temperature (Degree C)')
plt.title('Temperature')
plt.legend(loc=1, prop={'size':8})
plt.ylim(ymax = 5, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_DAMAGEb[0:T0], (GET_Gbt[0:T0]))], '--y', label = "CC Damage in Brown sector")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_DAMAGEg[0:T0], (GET_Ggt[0:T0]))], '--g', label = "CC Damage in Green sector")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_ABTCOSb[0:T0], (GET_Gbt[0:T0]))], '--b', label = "Abatement cost in Brown sector")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_ABTCOSg[0:T0], (GET_Ggt[0:T0]))], '--k', label = "Abatement cost in Green sector")
plt.xlabel('Time (years)')
plt.ylabel('Share of output')
plt.title('Climate change costs')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_CTKbgt[0:T0], (GET_Ggt[0:T0]))], '--g', label = "Capital: Brown --> Green")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_CTKgbt[0:T0], (GET_Gbt[0:T0]))], '--y', label = "Capital: Green --> Brown")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_CTLbgt[0:T0], (GET_Ggt[0:T0]))], '-g', label = "Labour: Brown --> Green")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_CTLgbt[0:T0], (GET_Gbt[0:T0]))], '-y', label = "Labour: Green --> Brown")
plt.xlabel('Time (years)')
plt.ylabel('Share of output')
plt.title('Transition costs')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_Ybt[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_Ygt[0:T0], '--g', label = "Green sector")
plt.plot(Xaxis[0:T0], GET_Yt[0:T0], 'r', label = "GET-Y")
plt.plot(Xaxis[0:T0], GET_GAMS_Yt[0:T0], 'k', label = "DICE")
plt.xlabel('Time (years)')
plt.ylabel('GDP ($trill, 2005$)')
plt.title('Economic Output')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 4, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Ybt[0:T0], GET_Yt[0:T0])], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Ygt[0:T0], GET_Yt[0:T0])], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('Share')
plt.title('Share of output')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 4, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

GET_Lt=[i + j for i, j in zip(GET_Lbt[0:T0], (GET_Lgt[0:T0]))]

plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Lbt[0:T0], (GET_Lt[0:T0]))], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Lgt[0:T0], (GET_Lt[0:T0]))], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('')
plt.title('Labor Share')
plt.legend(loc=1, prop={'size':8})
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

GET_Kt=[i + j for i, j in zip(GET_Kbt[0:T0], (GET_Kgt[0:T0]))]

plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Kbt[0:T0], (GET_Kt[0:T0]))], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Kgt[0:T0], (GET_Kt[0:T0]))], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('')
plt.title('Capital Share')
plt.legend(loc=1, prop={'size':8})
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_Abt[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_Agt[0:T0], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('')
plt.title('Capital Productivity')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_Bbt[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_Bgt[0:T0], '--g', label = "Green sector")
plt.xlabel('Time (years)')
plt.ylabel('')
plt.title('Labor Productivity')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 1.0, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

# plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Bbt[0:T0], (GET_Bgt[0:T0]))], '--m', label = "Labor")
# plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Abt[0:T0], (GET_Agt[0:T0]))], '--c', label = "Capital")
# plt.xlabel('Time (years)')
# plt.ylabel('Ratio')
# plt.title('Brown/Green Productivity Ratios')
# plt.legend(loc=1, prop={'size':8})
# plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
# plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
# plt.show()


# plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Lbt[0:T0], (GET_Gbt[0:T0]))], '--y', label = "Brown sector")
# plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Lgt[0:T0], (GET_Ggt[0:T0]))], '--g', label = "Green sector")
# plt.xlabel('Time (years)')
# plt.ylabel('Ratio')
# plt.title('Labour/Gross Output Ratio')
# plt.legend(loc=1, prop={'size':8})
# plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
# plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
# plt.show()

# plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Kbt[0:T0], (GET_Gbt[0:T0]))], '--y', label = "Brown sector")
# plt.plot(Xaxis[0:T0], [i / j for i, j in zip(GET_Kgt[0:T0], (GET_Ggt[0:T0]))], '--g', label = "Green sector")
# plt.xlabel('Time (years)')
# plt.ylabel('Ratio')
# plt.title('Capital/Gross Output Ratio')
# plt.legend(loc=1, prop={'size':8})
# plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
# plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
# plt.show()

plt.plot(Xaxis[0:T0], GET_Ebt[0:T0], '--y', label = "Brown sector")
plt.plot(Xaxis[0:T0], GET_Egt[0:T0], '--g', label = "Green sector")
plt.plot(Xaxis[0:T0], GET_Et[0:T0], 'r', label = "GET-Net emissions")
plt.plot(Xaxis[0:T0], GET_GAMS_Et[0:T0], 'k', label = "DICE")
plt.xlabel('Time (years)')
plt.ylabel('Emissions (GTCO2 per year)')
plt.title('Total CO2 emissions')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 4, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

plt.plot(Xaxis[0:T0], GET_alphab[0:T0], '--y', label = "Brown sector share")
plt.plot(Xaxis[0:T0], GET_alphag[0:T0], '--g', label = "Green sector share")
plt.xlabel('Time (years)')
plt.ylabel('Share')
plt.title('Share in total composite consumption')
plt.legend(loc=1, prop={'size':8})
#plt.ylim(ymax = 4, ymin = 0)
plt.xlim(xmax = max(Xaxis), xmin = min(Xaxis))
plt.xticks(np.arange(min(Xaxis), max(Xaxis), 4 * tstep))
plt.show()

# =============================================== Export into Excel =============================================== #
def output(filename, list1):
    book = xlwt.Workbook(filename)
    varname = ['Abatement rate (Green)', 'Abatement rate (Brown)', 'Atmospheric temperature', 'Total Capital', 'Capital (Green)', 'Capital (Brown)', 'Total Labor', 'Labor (Green)', 'Labor (Brown)', 'Net output', 'Net output (Green)', 'Net output (Brown)', 'net emissions', 'Emissions (Green)', 'Emissions (Brown)', 'Capital transfer (Brown->Green)', 'Capital transfer (Green->Brown)', 'Labor transfer (Brown->Green)', 'Labor transfer (Green->Brown)', 'Capital productivity share (Green)', 'Capital productivity share (Brown)', 'Labor productivity share (Green)', 'Labor productivity share (Brown)', 'Labor R&D share', 'Capital R&D share', 'Consumption (Green)', 'Consumption (Brown)', 'Saving rate (Green)', 'Saving rate (Brown)', 'Consumption share (Green)', 'Consumption share (Brown)', 'Gross output (Green)', 'Gross output (Brown)', 'Capital transfer cost (Brown->Green)', 'Capital transfer cost (Green->Brown)', 'Labor transfer cost (Brown->Green)', 'Labor transfer cost (Green->Brown)', 'Abatement cost (Green)', 'Abatement cost (Brown)', 'Damage cost (Green)', 'Damage cost (Brown)', 'Capital productivity (Green)', 'Capital productivity (Brown)', 'Labor productivity (Green)', 'Labor productivity (Brown)']
    sh = book.add_sheet('Results')
    for ind1 in range(len(list1)):
        sh.write(0, ind1, varname[ind1])
        for t in range(18):
            sh.write(1 + t, ind1, list1[ind1][t])
    book.save(filename)
    
output1 = [GET_miug, GET_miub, GET_TATm, GET_Kt, GET_Kgt, GET_Kbt, GET_Lt, GET_Lgt, GET_Lbt, GET_Yt, GET_Ygt, GET_Ybt, GET_Et, GET_Egt, GET_Ebt, GET_TKbgt, GET_TKgbt, GET_TLbgt, GET_TLgbt, GET_rg, GET_rb, GET_eg, GET_eb, GET_rdl, GET_rdk, GET_Cgt, GET_Cbt, GET_srg, GET_srb, GET_alphag, GET_alphab, GET_Ggt, GET_Gbt, GET_CTKbgt, GET_CTKgbt, GET_CTLbgt, GET_CTLgbt, GET_ABTCOSg, GET_ABTCOSb, GET_DAMAGEg, GET_DAMAGEb, GET_Agt, GET_Abt, GET_Bgt, GET_Bbt]
output('GET_Output_' + GET_case[cs] + '.xls', output1)
