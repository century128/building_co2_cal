import random

import numpy as np
from sympy import *
from matrix_cal import *
from scipy.optimize import minimize
import math
from scipy.optimize import SR1
# 此模块计算投资总费用
# 目标函数 min(LCCT) = lambda_cp*CE_total + CRF*m*C_invest + C_om - C_cer -C_re - zata
# 决策变量
# x = [P_ashp 0,A_stc 1,M_bio 2,P_eb 3,S_pv 4,g_gfb 5,g_cfb 6,g_ofb 7,P_gshp 8,P_wtg 9,M_hy 10,V_wt 11,V_ct 12,p_fc 13,
# E_g 14,m_hy2 15,hy_in 16,hy_out 17,t_wt 18,t_wt2 19,t_ct 20,t_ct2 21,P_ashpQ 22,P_gshpQ 23,P_ru 24,P_el 25,
# P_g_k:26-49,g_gfb_k:50-73,g_cfb_k:74-97,g_ofb_k:98-121,m_hyout_k:122-145,m_hy_k:146-169,p_ashpG_k:170-193,
# p_eb_k:194-217,p_gshp_k:218-241,P_ru_k:242-265,p_hyin_k:266-289,p_hyout_k:290-313,m_bio_k:314-337,g_wt_k:338-361,
# g_wt_in:362-385,g_wt_out:386-409,q_ct_k:410-433,q_ct_in:434-457,q_ct_out:458-481,p_gshpG_k:482-505,p_gshpQ_k:506-529
# g_bhp_k:530-553,g_bhp_in_k:554-577,p_ashpQ_k:578-601,m_hy_purIn:602-625]
#约束条件 电平衡约束，热平衡约束,热水罐约束,冷水罐约束

#----------目标函数-----------------------
def object(x):
    #投资参数按照word模型顺序给出
    #C_invest目标
    c_invest = [3000,1500,360,160,6600,3980,3980,3980,400,3600,3000,500,500,10000,3000]
    invest = 0
    for i in range(14):
        invest += c_invest[i] * x[i]
    invest += x[24] * 1000 #加入制冷机ru
    invest += x[25] * 10000 #加入电解槽el

    #------CE_total目标------------
    CE_total,CE_es,CE_cpesr= 0,0,3

    #CE_an计算
    CE_an = 16068.17

    #CE_es计算 能源系统投产碳排
    CE_ashp = x[0]*2.83/15*2.88#y_e,ashp为15年，GWP为2880
    CE_stc = 25.2*10**(-6)*x[1]
    CE_bb = 1788*10**(-6)*x[2]
    CE_eb = 650*10**(-6)*x[3]
    CE_pv = 34*10**(-6)*x[4]
    CE_cfb = 876*10**(-6)*x[6]
    CE_gfb = 412*10**(-6)*x[5]
    CE_ofb = 579*10**(-6)*x[7]
    CE_gshp = 79*10**(-6)*x[8]
    CE_wtg = 17*10**(-6)*x[9]
    CE_hy = 51*10**(-6)*x[13]
    CE_wt = 500*10**(-6)*x[11]
    CE_ct = 500*10**(-6)*x[12]
    CE_ru = 2000*10**(-6)*x[24]
    CE_r = CE_ashp + CE_ru +CE_gshp
    CE_es = CE_ashp+CE_stc+CE_bb+CE_eb+CE_pv+CE_cfb+CE_gfb+CE_ofb+CE_gshp+CE_wtg+CE_hy+CE_wt+CE_ct+CE_ru+CE_r

    #CE_cpesr计算 运行碳排
    #!!!!此处改成8760的求和
    EF_g,EF_gas,EF_coil,EF_oil = 0.6101,1.535,8.14,2.25
    beta_gfb,beta_cfb,beta_ofb = 0.9,0.84,0.9
    CE_g,CE_coil,CE_gas,CE_oil = 0,0,0,0
    for i in range(24):
        CE_g += x[i+26]
        CE_coil += x[i+74]
        CE_gas += x[i+50]
        CE_oil += x[i+98]
    CE_g = EF_g*CE_g
    CE_coil = EF_coil*CE_coil/beta_cfb
    CE_gas = EF_gas*CE_coil/beta_gfb
    CE_oil = EF_oil*CE_oil/beta_ofb
    CE_cpesr = CE_g + CE_coil + CE_gas + CE_oil

    # EF_coil,EF_gas,EF_oil,EF_bio,EF_g = 8.14,1.535,2.25,1.78,0.5839
    # CE_cpesr = EF_g*x[14] + CE_r + EF_bio*15000*x[2] + EF_coil*x[6] + EF_gas*x[5] + EF_oil*x[7]
    CE_total = CE_an + CE_es + CE_cpesr

    #O&M目标 !!!!!此处改成8760求和
    om_para = [0.58,0.8,1.2,3.5,3.75,40,0.02]#各能源实时价格参数表[lep,bio,coil,gas,oil,hy,sigma]
    p_e,m_bio,m_coal,v_gas,v_oil,m_hy_purIn = 0,0,0,0,0,0
    for i in range(24):
        p_e += x[i+26]
        m_bio += x[i+341]
        m_coal += x[i+74]
        v_gas += x[i+50]
        v_oil += x[i+122]
        m_hy_purIn += x[i+602]
    p_e = om_para[0]*p_e
    m_bio = om_para[1]*m_bio
    m_coal *= (om_para[2]/beta_cfb)
    v_gas *= (om_para[3]/beta_gfb)
    v_oil *= (om_para[4]/beta_ofb)
    m_hy_purIn *= om_para[5]

    C_om = p_e+m_bio+m_coal+v_gas+v_oil+m_hy_purIn+invest*om_para[6]

    #re约束
    r_equ = 0.02 #设备残值系数
    C_re = invest * r_equ

    lamda_cp = 0.06 #碳价
    zeta = 0.1 * x[4] * 300  #新能源补贴金额
    i,m = 0.08,15
    CRF = i*math.pow(1+i,m)/(math.pow(1+i,m)-1)
    return lamda_cp*CE_total + CRF*m*invest + C_om  - C_re - zeta

#---------------------约束条件--------------------------------------
def P_cons(x):
    # p_d = np.ones(24,dtype='float32')
    p_d = [166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,]
    # pv_I_k = np.ones(24,dtype='float32')
    pv_I_k = [0,0,0,0,0,0,0,0,0.066,0.37,0.587,0.711,0.76,0.746,0.667,0.514,0.271,0.019,0,0,0,0,0,0]
    # p_wtg_fai = np.ones(24,dtype='float32')
    p_wtg_fai = [0.042,0.009,0.001,0,0.006,0.021,0.028,0.021,0.026,0.021,0.002,0,0.005,0.012,0.015,0.012,0.012,0.026,0.021,0.025,0.038,0.046,0.045,0.031]
    p_con = np.zeros(24,dtype='float32')
    for i in range(24):
        p_con[i] = x[26 + i] + 0.75 * pv_I_k[i] * x[4] + p_wtg_fai[i] * x[9] + x[i + 290] - p_d[i] - x[170 + i] - x[578 + i] - x[194 + i] - x[482 + i] - x[506 + i] - x[i + 266] - x[i + 242]
    # return x[0] + x[3] + x[6] + x[5] + x[7] + x[9] + (x[15]+x[17]-x[16]-x[10])*P_hy_para[0] + P_para[3] + P_para[1] + P_para[2] + x[2] - x[14] - P_para[0] -P_pv_para[0]*P_pv_para[1]*(1-P_pv_para[2])*x[4] - P_hy_para[1]*x[17]
    return p_con
def G_cons(x):
    # 热守恒约束参数
    # I_k = np.ones(24,dtype='float32')
    I_k = [166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,166.36,]
    # G_d = np.ones(24,dtype='float32')
    G_d = [1136.15, 1247.82, 1320.02, 1406.41, 1463.76, 1507.25, 4433.76, 4486.54, 4267.33, 3051.05, 3450.79, 3178.71, 1228.43, 589.83, 1079.59, 1484.18, 1732.12, 2193.13, 1035.93, 1179.03, 1277.13, 1141.99, 1218.63, 1294.44, ]
    G_w = np.zeros(24,dtype='float32')
    G_con = np.zeros(24,dtype='float32')
    for i in range(24):
        G_con[i] = 4*x[i+170]+0.75*I_k[i]*x[1]+x[314+i]*15000*0.76+x[i+194]*0.97+x[98+i]*16.6+x[i+50]+x[i+74]+x[i+122]+x[482+i]*4+x[i+386]-G_d[i]-G_w[i]-x[i+554]-x[i+362]-0.01*x[i+338]
    return G_con
    #return G_ashp + G_stc + G_bb + G_eb + G_hyOut + G_cfb + G_gfb + G_ofb + G_gshb - G_wt - G_d - G_dhw
def Q_cons(x):
    # 冷守恒约束参数
    Q_con = np.zeros(24,dtype='float32')
    Q_d = np.zeros(24,dtype='float32')
    for i in range(24):
        Q_con[i]=x[i+578]*5+5*x[i+506]+x[i+458]+x[i+242]*2.5-Q_d[i]-x[i+434]-x[i+410]*0.01
    #return Q_ashp + Q_gshp -Q_ct - Q_d
    return Q_con

# ----------国标约束条件---------------
def CE_GB_incons1(x):
    # ------CE_total目标------------

    # CE_an计算
    CE_an = 16068.17

    # CE_es计算 能源系统投产碳排
    CE_ashp = x[0] * 2.83 / 15 * 2.88  # y_e,ashp为15年，GWP为2880
    CE_stc = 25.2 * 10 ** (-6) * x[1]
    CE_bb = 1788 * 10 ** (-6) * x[2]
    CE_eb = 650 * 10 ** (-6) * x[3]
    CE_pv = 34 * 10 ** (-6) * x[4]
    CE_cfb = 876 * 10 ** (-6) * x[6]
    CE_gfb = 412 * 10 ** (-6) * x[5]
    CE_ofb = 579 * 10 ** (-6) * x[7]
    CE_wtg = 17 * 10 ** (-6) * x[9]
    CE_gshp = 79 * 10 ** (-6) * x[8]
    CE_hy = 51 * 10 ** (-6) * x[13]
    CE_wt = 500 * 10 ** (-6) * x[11]
    CE_ct = 500 * 10 ** (-6) * x[12]
    CE_ru = 2000 * 10 ** (-6) * x[24]
    CE_r = CE_ashp + CE_ru + CE_gshp
    CE_es = CE_ashp + CE_stc + CE_bb + CE_eb + CE_pv + CE_cfb + CE_gfb + CE_ofb + CE_gshp + CE_wtg + CE_hy + CE_wt + CE_ct + CE_ru + CE_r

    # CE_cpesr计算 运行碳排
    # !!!!此处改成8760的求和
    EF_g, EF_gas, EF_coil, EF_oil = 0.6101, 1.535, 8.14, 2.25
    beta_gfb, beta_cfb, beta_ofb = 0.9, 0.84, 0.9
    CE_g, CE_coil, CE_gas, CE_oil = 0, 0, 0, 0
    for i in range(24):
        CE_g += x[i + 26]
        CE_coil += x[i + 74]
        CE_gas += x[i + 50]
        CE_oil += x[i + 122]
    CE_g = EF_g * CE_g
    CE_coil = EF_coil * CE_coil / beta_cfb
    CE_gas = EF_gas * CE_coil / beta_gfb
    CE_oil = EF_oil * CE_oil / beta_ofb
    CE_cpesr = CE_g + CE_coil + CE_gas + CE_oil

    building_area = 156027.93  # A
    CE_total_2016 = 35 * building_area
    # EF_coil,EF_gas,EF_oil,EF_bio,EF_g = 8.14,1.535,2.25,1.78,0.5839
    # CE_cpesr = EF_g*x[14] + CE_r + EF_bio*15000*x[2] + EF_coil*x[6] + EF_gas*x[5] + EF_oil*x[7]
    CE_total = CE_an + CE_es + CE_cpesr
    return 0.6 * CE_total_2016 - CE_total
def CE_GB_incons2(x):
    # ------CE_total目标------------
    CE_total, CE_es, CE_cpesr = 0, 0, 3

    # CE_an计算
    CE_an = 16068.17

    # CE_es计算 能源系统投产碳排
    CE_ashp = x[0] * 2.83 / 15 * 2.88  # y_e,ashp为15年，GWP为2880
    CE_stc = 25.2 * 10 ** (-6) * x[1]
    CE_bb = 1788 * 10 ** (-6) * x[2]
    CE_eb = 650 * 10 ** (-6) * x[3]
    CE_pv = 34 * 10 ** (-6) * x[4]
    CE_cfb = 876 * 10 ** (-6) * x[6]
    CE_gfb = 412 * 10 ** (-6) * x[5]
    CE_ofb = 579 * 10 ** (-6) * x[7]
    CE_gshp = 79 * 10 ** (-6) * x[8]
    CE_wtg = 17 * 10 ** (-6) * x[9]
    CE_hy = 51 * 10 ** (-6) * x[13]
    CE_wt = 500 * 10 ** (-6) * x[11]
    CE_ct = 500 * 10 ** (-6) * x[12]
    CE_ru = 2000 * 10 ** (-6) * x[24]
    CE_r = CE_ashp + CE_ru + CE_gshp
    CE_es = CE_ashp + CE_stc + CE_bb + CE_eb + CE_pv + CE_cfb + CE_gfb + CE_ofb + CE_gshp + CE_wtg + CE_hy + CE_wt + CE_ct + CE_ru + CE_r

    # CE_cpesr计算 运行碳排
    # !!!!此处改成8760的求和
    EF_g, EF_gas, EF_coil, EF_oil = 0.6101, 1.535, 8.14, 2.25
    beta_gfb, beta_cfb, beta_ofb = 0.9, 0.84, 0.9
    CE_g, CE_coil, CE_gas, CE_oil = 0, 0, 0, 0
    for i in range(24):
        CE_g += x[i + 26]
        CE_coil += x[i + 74]
        CE_gas += x[i + 50]
        CE_oil += x[i + 122]
    CE_g = EF_g * CE_g
    CE_coil = EF_coil * CE_coil / beta_cfb
    CE_gas = EF_gas * CE_coil / beta_gfb
    CE_oil = EF_oil * CE_oil / beta_ofb
    CE_cpesr = CE_g + CE_coil + CE_gas + CE_oil

    building_area = 156027.93  # A
    CE_total_2016 = 35 * building_area
    # EF_coil,EF_gas,EF_oil,EF_bio,EF_g = 8.14,1.535,2.25,1.78,0.5839
    # CE_cpesr = EF_g*x[14] + CE_r + EF_bio*15000*x[2] + EF_coil*x[6] + EF_gas*x[5] + EF_oil*x[7]
    CE_total = CE_an + CE_es + CE_cpesr
    return 10.5*building_area+CE_total_2016-CE_total

#-----------氢约束----------------
def hy_cons1(x):
    m_hy = np.zeros(23,dtype='float32')
    for i in range(23):
        m_hy[i]=x[i+147]-x[i+146]-x[i+602]-1.399*x[i+266]+1/15*x[i+290]
    return m_hy
def hy_incons1(x):
    p_hy_out = np.zeros(24,dtype='float32')
    for i in range(24):
        p_hy_out[i] = x[13]-x[i+290]
    return p_hy_out
def hy_incons2(x):
    p_hy_in = np.zeros(24,dtype='float32')
    for i in range(24):
        p_hy_in[i] = x[25] - x[i + 266]
    return p_hy_in
def hy_incons3(x):
    m_hy = np.zeros(24,dtype='float32')
    for i in range(24):
        m_hy[i] = x[10] - x[i + 146]
    return m_hy
#-----------地源热泵--------------
def gshp_incons(x):
    P_gshp = np.zeros(24,dtype='float32')
    for i in range(24):
        P_gshp[i]=x[8]-x[i+482]-x[i+506]
    return P_gshp
#-----------土壤储热系统单元-----------------------
def bhp_cons(x):
    g_bhp = np.zeros(23,dtype='float32')
    for i in range(23):
        g_bhp[i] = x[i+531]-x[i+530]-5*x[i+506]*1.2+0.75*4*x[i+482]-x[i+554]
    return g_bhp
#---------------生物质锅炉子系统-----------------
def bio_incons(x):
    m_bio = np.zeros(24,dtype='float32')
    for i in range(24):
        m_bio[i] = x[2]-x[i+314]
    return m_bio
#--------------空气源热泵---------------------
def ashp_incons1(x):
    p_ashp = np.zeros(24,dtype='float32')
    for i in range(24):
        p_ashp[i] = x[0]- x[i+170]
    return p_ashp
def ashp_incons2(x):
    p_ashp = np.zeros(24,dtype='float32')
    for i in range(24):
        p_ashp[i] = x[0]- x[i+578]
    return p_ashp
#--------------电锅炉---------------------
def eb_incons(x):
    p_eb = np.zeros(24,dtype='float32')
    for i in range(24):
        p_eb[i]=x[3]-x[i+194]
    return p_eb
#--------------燃气锅炉-----------------
def gfb_incons(x):
    g_gfb = np.zeros(24,dtype='float32')
    for i in range(24):
        g_gfb[i] = x[5] - x[i + 50]
    return g_gfb
#--------------燃煤锅炉----------------
def cfb_incons(x):
    g_cfb = np.zeros(24,dtype='float32')
    for i in range(24):
        g_cfb[i] = x[6] - x[i + 74]
    return g_cfb
#-------------燃油锅炉-----------------
def ofb_incons(x):
    g_ofb = np.zeros(24,dtype='float32')
    for i in range(24):
        g_ofb[i] = x[7] - x[i + 122]
    return g_ofb
#-------------制冷机组----------------
def ru_incons(x):
    p_ru = np.zeros(24,dtype='float32')
    for i in range(24):
        p_ru[i] = x[24] - x[i + 242]
    return p_ru
#---------热水罐----------------------
def wt_cons(x):
    g_wt = np.zeros(23,dtype='float32')
    for i in range(23):
        g_wt[i]= x[i+339]-x[i+338]-x[i+362]+x[i+386]+0.01*x[i+338]
    return g_wt
def wt_incons(x):
    g_wt = np.zeros(24,dtype='float32')
    for i in range(24):
        g_wt[i] = 1000*4.2*(85-45)*x[11]-x[i+338]
    return g_wt
#------------冷水罐-----------------
def ct_cons(x):
    g_ct = np.zeros(23,dtype='float32')
    for i in range(23):
        g_ct[i]= x[i+411]-x[i+410]-x[i+434]+x[i+458]+0.01*x[i+410]
    return g_ct
def ct_incons(x):
    g_ct = np.zeros(24,dtype='float32')
    for i in range(24):
        g_ct[i] = 1000*4.2*(21-4)*x[12]-x[i+410]
    return g_ct

#设定约束
cons1 = {'type':'eq','fun':P_cons }
cons2 = {'type':'eq','fun':G_cons }
cons3 = {'type':'eq','fun':Q_cons }
cons4 = {'type':'eq','fun':hy_cons1 }
cons5 = {'type':'eq','fun':bhp_cons}
cons6 = {'type':'eq','fun':wt_cons }
cons7 = {'type':'eq','fun':ct_cons }

cons8 = {'type':'ineq','fun':CE_GB_incons1}
cons9 = {'type':'ineq','fun':CE_GB_incons2 }
cons10 = {'type':'ineq','fun':hy_incons1 }
cons11 = {'type':'ineq','fun':hy_incons2 }
cons12 = {'type':'ineq','fun':hy_incons3 }
cons13 = {'type':'ineq','fun':gshp_incons }
cons14 = {'type':'ineq','fun':bio_incons }
cons15 = {'type':'ineq','fun':ashp_incons1 }
cons16 = {'type':'ineq','fun':ashp_incons2 }
cons17 = {'type':'ineq','fun':eb_incons}
cons18 = {'type':'ineq','fun':gfb_incons }
cons19 = {'type':'ineq','fun':ofb_incons }
cons20 = {'type':'ineq','fun':cfb_incons }
cons21 = {'type':'ineq','fun':ru_incons }
cons22 = {'type':'ineq','fun':wt_incons }
cons23 = {'type':'ineq','fun':ct_incons }

cons = [cons1, cons2, cons3, cons4, cons5, cons6, cons7, cons8, cons9,cons10,cons11,cons12,cons13,cons14,cons15,cons16,cons17,cons18,cons19,cons20,cons21,cons22,cons23]
#-----------------------------------变量上下界------------------------------------
# #S_pv bounds
# S_pv_para =[10,1] #Q_hmax,n
# S_pv_ul = S_pv_para[0]/S_pv_para[1] #ul:upper limit
# S_pv = (0,S_pv_ul)
#
# #M_bb bounds
# M_bb_para = [10,1,0.2] #Q_h_max,q_lcv,n_bio
# M_bb_ul = M_bb_para[0]/M_bb_para[1]/(1-M_bb_para[2])
# M_bb = (0,M_bb_ul)
#
# #P_ashp bounds
# P_ashp_para = [10,5] #Q_h_max,COP_ashp
# P_ashp_ll = 0.5*P_ashp_para[0]/P_ashp_para[1] #ll lower limit
# P_ashp_ul = 2*P_ashp_ll
# P_ashp = (P_ashp_ll,P_ashp_ul)
#
# #eb bounds
# P_eb_ul = 100
# P_eb = (0,P_eb_ul)
#
# #gfb bounds
# P_gfb_ul = 100
# P_gfb = (0,P_gfb_ul)
#
# #cfb bounds
# P_cfb_ul = 100
# P_cfb = (0,P_cfb_ul)
#
# #ofb bounds
# P_ofb_ul = 100
# P_ofb = (0,P_ofb_ul)
#
# #t_wt bounds
# t_wt_ll,t_wt_ul = 20,60
# t_wt = (t_wt_ll,t_wt_ul)
#
# #t_cl bounds
# t_ct_ll , t_ct_ul = 0,50
# t_ct = (t_ct_ll,t_ct_ul)
#
# #ashp_q bounds
# ashp_q_para = [0.5,100] #z_q,P_ashp_ul
# ashp_q_ul = ashp_q_para[0]*ashp_q_para[1]
# ashp_q = (0,ashp_q_ul)
#
# #gshp_q_bounds
# gshp_q_para = [0.5,100] #z_q,g_ashp_ul
# gshp_q_ul = gshp_q_para[0]*gshp_q_para[1]
# gshp_q = (0,gshp_q_ul)
#
# b = (0,np.inf)
# bnds = (P_ashp,b,M_bb,P_eb,S_pv,P_gfb,P_cfb,P_ofb,b,b,b,b,b,b,b,b,b,b,t_wt,t_wt,t_ct,t_ct,ashp_q,gshp_q)

# 设定初值
# x0 = np.zeros(626,dtype='float32')
# for i in range(626):
#     x0[i] = 0.001

x0 = np.ones(626)
np.random.seed(0)
a = np.random.rand(626)+x0
#设定边界
b = [(0, np.inf)]
bnds = b * 626
#求解
#res = minimize(object, x0, method='SLSQP', bounds = bnds, constraints=cons, jac="2-point",)
res = minimize(object, x0 = a,bounds=bnds,constraints=cons, jac="3-point",options={ 'disp': True,'maxiter':1000},)
# print(res.message)
# print(res.success)
print(res.fun)