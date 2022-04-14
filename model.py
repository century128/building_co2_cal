import numpy as np
from sympy import *
from matrix_cal import *
from scipy.optimize import minimize
import math
from scipy.optimize import SR1
# 此模块计算投资总费用
#目标函数 min(LCCT) = lambda_cp*CE_total + CRF*m*C_invest + C_om - C_cer -C_re - zata
#决策变量 x = [P_ashp 0,A_stc 1,M_bb 2,P_eb 3,S_pv 4,g_gfb 5,g_cfb 6,g_ofb 7,P_gshp 8,P_wtg 9,M_hy 10,V_wt 11,V_ct 12,p_fc 13,
# E_g 14,m_hy2 15,hy_in 16,hy_out 17,t_wt 18,t_wt2 19,t_ct 20,t_ct2 21,P_ashpQ 22,P_gshpQ 23,P_ru 24,P_el 25]
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
    beta_gfb,beta_cfb,beta_ofb = 9.3,3.513,29.3
    CE_g = x[14]*EF_g
    CE_coil = EF_coil*x[6]/beta_ofb
    CE_gas = EF_gas*x[5]/beta_gfb
    CE_oil = EF_oil*x[7]/beta_ofb
    CE_cpesr = CE_g + CE_coil + CE_gas + CE_oil

    # EF_coil,EF_gas,EF_oil,EF_bio,EF_g = 8.14,1.535,2.25,1.78,0.5839
    # CE_cpesr = EF_g*x[14] + CE_r + EF_bio*15000*x[2] + EF_coil*x[6] + EF_gas*x[5] + EF_oil*x[7]
    CE_total = CE_an + CE_es + CE_cpesr

    #O&M目标 !!!!!此处改成8760求和
    om_para = [0.58,0.8,1.2,3.5,3.75,40,0.02]#各能源实时价格参数表[lep,bio,coil,gas,oil,hy,sigma]
    C_om = om_para[0]*x[15] + om_para[1]*x[2] + om_para[2]*x[6] + om_para[3]*x[5] + om_para[4]*x[7] + om_para[5]*x[2] + om_para[6] * invest

    #re约束
    r_equ = 0.02 #设备残值系数
    C_re = invest * r_equ
    lamda_cp = 0.06 #碳价
    zeta = 0.1 * x[4] * 300  #新能源补贴金额
    i,m = 0.08,15
    CRF = i*math.pow(1+i,m)/(math.pow(1+i,m)-1)
    return lamda_cp*CE_total + CRF*m*invest + C_om  - C_re - zeta

#---------------------约束条件--------------------------------------
def pv_incons(x):
    return 1000000-0.1 * x[4] * 300

def P_cons(x):
    # 电守恒约束参数 [E_wtg,E_l,E_e,E_oe]
    #P_k_d 外界输入
    P_para = [1000,1000,1000,1000]
    P_pv_para = [1,0.21] #I_k,gama,导入当地光照强度I_k
    P_hy_para = [0.715,15] #1/beta_el,beta_fc_g
    p_ashp_para = [4,1]#n_ashp_g,t
    P_gshp_para = []
    return x[0] + x[3] + x[6] + x[5] + x[7] + x[9] + (x[15]+x[17]-x[16]-x[10])*P_hy_para[0] + P_para[3] + P_para[1] + P_para[2] + x[2] - x[14] - P_para[0] -P_pv_para[0]*P_pv_para[1]*(1-P_pv_para[2])*x[4] - P_hy_para[1]*x[17]

def P_incons_hy1(x):#ineq
   E_hy_para = [16.6,1.399,15,1] #beta_fcg,beta_elco,beta_fce,t
   return x[13]-E_hy_para[1]*x[17]

def P_incons_hy2(x):#ineq
   E_hy_para = [0.715,15] #1/beta_el,beta_fc
   return x[14]-(x[15]+x[17]-x[16]-x[10])*E_hy_para[0]

def G_cons(x):
    # 热守恒约束参数
    G_ashp_para = [5]#[n_ashp_q]
    G_ashp = x[0]*G_ashp_para[0]*G_ashp_para[1]

    G_stc_para = [0.75,3.4] #[beta_stc,I_k]
    G_stc = x[4]*G_stc_para[0]*G_stc_para[1]

    G_bb_para = [15000]#q_lcv_bio
    G_bb = x[2]*G_bb_para[0]

    G_eb_para = [0.97]#beta_eb
    G_eb = G_eb_para[0]*x[3]

    G_hy_para = [1] #beta_fc_g
    G_hyOut = G_hy_para[0]*x[17]

    beta_para = [0.88,0.9,0.9,4] #cfb,gfb,ofb,gshb
    G_cfb = beta_para[0]*x[6]
    G_gfb = beta_para[1]*x[5]
    G_ofb = beta_para[2]*x[7]
    G_gshb = beta_para[3]*x[0]

    G_wt_para = [1,1,0.02,1,45] # c,p,u_loss,A_wt,t_envIn
    G_wt = G_wt_para[0]*G_wt_para[1]*x[11]*(x[19]-x[18]) + G_wt_para[2]*G_wt_para[3]*(x[18]-G_wt_para[4])

    G_d,G_dhw = 100,50

    return G_ashp + G_stc + G_bb + G_eb + G_hyOut + G_cfb + G_gfb + G_ofb + G_gshb - G_wt - G_d - G_dhw

def Q_cons(x):
    # 冷守恒约束参数
    Q_para = [1,1,1,1,1,1,1]#n_ashpQ,n_qshp,c_ct,p_ct,u_ctLoss,A_ct,t_env
    Q_ashp = Q_para[0]*x[22]
    Q_gshp = Q_para[1]*x[23]
    Q_ct = Q_para[2]*Q_para[3]*x[12]*(x[20]-x[21])+Q_para[4]*Q_para[5]*(Q_para[6]-x[20])
    Q_d = 1000
    return Q_ashp + Q_gshp -Q_ct - Q_d

#国标约束条件
def CE_GB_incons1(x):
    # CE_total目标
    CE_total, CE_an, CE_es, CE_cpesr = 0, 1, 2, 3
    CE_r = 1  # 冷凝剂碳排
    EF_coil, EF_gas, EF_oil = 1, 1, 1
    CE_cpesr = 0.94 * x[14] + CE_r + EF_coil * x[6] + EF_gas * x[5] + EF_oil * x[7]
    # 此处写CE_cpesr约束
    CE_total_2016 = 100
    CE_total = CE_an + CE_es + CE_cpesr
    return 0.6*CE_total_2016-CE_total

def CE_GB_incons2(x):
    # CE_total目标
    CE_total, CE_an, CE_es, CE_cpesr = 0, 1, 2, 3
    CE_r = 1  # 冷凝剂碳排
    EF_coil, EF_gas, EF_oil = 1, 1, 1
    CE_cpesr = 0.94 * x[14] + CE_r + EF_coil * x[6] + EF_gas * x[5] + EF_oil * x[7]
    # 此处写CE_cpesr约束
    CE_total_2016 = 100
    A = 100
    CE_total = CE_an + CE_es + CE_cpesr
    return 7*A-CE_total+CE_total_2016
'''
def w_cons(x):
    #生活热水约束
    return

def ct_cons(x):
    #热水罐约束
    return
'''
def hy_incons1(x):
    #氢系统约束
    E_hy_para = [1, 15]  # 1/beta_el,beta_fc
    return x[13] - E_hy_para[1]*x[17]

def hy_incons2(x):
    E_hy_para = [1, 15]  # 1/beta_el,beta_fc
    return x[14]-(x[15]+x[17]-x[16]-x[10])*E_hy_para[0]

'''
#供热系统约束
def stc_constraint(x):
    return
'''

#设定约束
cons1 = {'type':'eq','fun':P_cons }
cons2 = {'type':'eq','fun':G_cons }
cons3 = {'type':'eq','fun':Q_cons }

cons4 = {'type':'ineq','fun':P_incons_hy1 }
cons5 = {'type':'ineq','fun':P_incons_hy2 }
cons6 = {'type':'ineq','fun':CE_GB_incons1 }
cons7 = {'type':'ineq','fun':CE_GB_incons2 }
cons8 = {'type':'ineq','fun':hy_incons1 }
cons9 = {'type':'ineq','fun':hy_incons2 }

cons = [cons1, cons2, cons3, cons4, cons5, cons6, cons7, cons8, cons9]
#-----------------------------------变量上下界------------------------------------
#S_pv bounds
S_pv_para =[10,1] #Q_hmax,n
S_pv_ul = S_pv_para[0]/S_pv_para[1] #ul:upper limit
S_pv = (0,S_pv_ul)

#M_bb bounds
M_bb_para = [10,1,0.2] #Q_h_max,q_lcv,n_bio
M_bb_ul = M_bb_para[0]/M_bb_para[1]/(1-M_bb_para[2])
M_bb = (0,M_bb_ul)

#P_ashp bounds
P_ashp_para = [10,5] #Q_h_max,COP_ashp
P_ashp_ll = 0.5*P_ashp_para[0]/P_ashp_para[1] #ll lower limit
P_ashp_ul = 2*P_ashp_ll
P_ashp = (P_ashp_ll,P_ashp_ul)

#eb bounds
P_eb_ul = 100
P_eb = (0,P_eb_ul)

#gfb bounds
P_gfb_ul = 100
P_gfb = (0,P_gfb_ul)

#cfb bounds
P_cfb_ul = 100
P_cfb = (0,P_cfb_ul)

#ofb bounds
P_ofb_ul = 100
P_ofb = (0,P_ofb_ul)

#t_wt bounds
t_wt_ll,t_wt_ul = 20,60
t_wt = (t_wt_ll,t_wt_ul)

#t_cl bounds
t_ct_ll , t_ct_ul = 0,50
t_ct = (t_ct_ll,t_ct_ul)

#ashp_q bounds
ashp_q_para = [0.5,100] #z_q,P_ashp_ul
ashp_q_ul = ashp_q_para[0]*ashp_q_para[1]
ashp_q = (0,ashp_q_ul)

#gshp_q_bounds
gshp_q_para = [0.5,100] #z_q,g_ashp_ul
gshp_q_ul = gshp_q_para[0]*gshp_q_para[1]
gshp_q = (0,gshp_q_ul)

b = (0,np.inf)
bnds = (P_ashp,b,M_bb,P_eb,S_pv,P_gfb,P_cfb,P_ofb,b,b,b,b,b,b,b,b,b,b,t_wt,t_wt,t_ct,t_ct,ashp_q,gshp_q)

# 设定初值
x0 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

#求解
res = minimize(object, x0, method='SLSQP', bounds = bnds, constraints=cons, jac="2-point",)

print(res)
