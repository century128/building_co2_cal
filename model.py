import numpy as np
from sympy import *
from matrix_cal import *
from scipy.optimize import minimize

# 此模块计算投资费用C_invest
#目标函数 CE_toal = sum(c[i]*p[i]),C为单位造价，P为设备容量
#决策变量 x = [P_ashp 0,A_stc 1,M_bb 2,P_eb 3,S_pv 4,P_gfb 5,P_cfb 6,P_ofb 7,P_gshp 8,P_wtg 9,M_hy 10,V_wt 11,V_ct 12,p_fc 13,
# E_g 14,m_hy2 15,hy_in 16,hy_out 17,t_wt 18,t_wt2 19,t_ct 20,t_ct2 21]
#约束条件 电平衡约束，热平衡约束,热水罐约束,冷水罐约束

#目标函数
def object(x):
    #投资参数按照word模型顺序给出
    #invest目标
    c_invest = [2000,1500,360,160,6600,3980,3980,3980,400,5000,12500,500,500,3000]
    invest = 0
    for i in range(13):
        invest += c_invest[i] * x[i]
    invest += c_invest[13]

    #CE_total目标
    CE_total,CE_an,CE_es,CE_cpesr= 0,1,2,3
    CE_r = 1 #冷凝剂碳排
    EF_coil,EF_gas,EF_oil = 1,1,1
    CE_cpesr = 0.94*x[14] + CE_r + EF_coil*x[6] + EF_gas*x[5] + EF_oil*x[7]
    #此处写CE_cpesr约束
    CE_total = CE_an + CE_es + CE_cpesr

    #O&M目标
    om_para = [0.54,0.8,1,4,5,6,7]#各能源实时价格参数表[lep,bio,coil,gas,oil,hy,sigma]
    C_om = om_para[0]*x[15] + om_para[1]*x[2] + om_para[2]*x[6] + om_para[3]*x[5] + om_para[4]*x[7] + om_para[5]*x[2] + om_para[6] * invest

    #CER约束
    lamda_cp,CE_aa, = 1,100
    C_cer = (CE_aa - CE_total)*lamda_cp

    #re约束
    r_equ = 0.2 #设备残值系数
    C_re = invest * r_equ

    zeta = 100 #新能源补贴金额
    CRF,m = 1,1
    return lamda_cp*CE_total + CRF*m*invest + C_om - C_cer - C_re - zeta

def P_cons(x):
    # 电守恒约束参数 [E_wtg,E_l,E_e,E_oe]
    P_para = [1,1,1,1]
    E_pv_para = [1,1,0.5] #I_k,K_e,K_s,
    E_hy_para = [1,15] #1/beta_el,beta_fc
    return x[0] + x[3] + x[6] + x[5] + x[7] + x[9] + (x[15]+x[17]-x[16]-x[10])*E_hy_para[0] + P_para[3] + P_para[1] + P_para[2] + x[2] - x[14] - P_para[0] -E_pv_para[0]*E_pv_para[1]*(1-E_pv_para[2])*x[4] - E_hy_para[1]*x[17]

def P_incons_hy1(x):#ineq
   E_hy_para = [1,15] #1/beta_el,beta_fc
   return x[13]-E_hy_para[1]*x[17]

def P_incons_hy2(x):#ineq
   E_hy_para = [1,15] #1/beta_el,beta_fc
   return x[14]-(x[15]+x[17]-x[16]-x[10])*E_hy_para[0]

def G_cons(x):
    # 热守恒约束参数
    G_ashp_para = [1,1]#[s_ashp,cop_ashp]
    G_ashp = x[0]*G_ashp_para[0]*G_ashp_para[1]

    G_stc_para = [1,1] #[beta_stc,I_k]
    G_stc = x[4]*G_stc_para[0]*G_stc_para[1]

    G_bb_para = [15000]#q_lcv_bio
    G_bb = x[2]*G_bb_para[0]

    G_eb_para = [1]#beta_eb
    G_eb = G_eb_para[0]*x[3]

    G_hy_para = [1] #beta_fc_g
    G_hyOut = G_hy_para[0]*x[17]

    beta_para = [1,1,1,1] #cfb,gfb,ofb,gshb
    G_cfb = beta_para[0]*x[6]
    G_gfb = beta_para[1]*x[5]
    G_ofb = beta_para[2]*x[7]
    G_gshb = beta_para[3]*x[0]

    G_wt_para = [1,1,1,1,1] # c,p,u_loss,A_wt,t_envIn
    G_wt = G_wt_para[0]*G_wt_para[1]*x[11]*(x[19]-x[18]) + G_wt_para[2]*G_wt_para[3]*(x[18]-G_wt_para[4])

    G_d,G_dhw = 100,50

    return G_ashp + G_stc + G_bb + G_eb + G_hyOut + G_cfb + G_gfb + G_ofb + G_gshb - G_wt - G_d - G_dhw

def Q_cons(x):
    # 冷守恒约束参数
    Q_para = [1,1,1,1,1,1,1]#n_ashpQ,n_qshp,c_ct,p_ct,u_ctLoss,A_ct,t_env
    Q_ashp = Q_para[0]*x[0]
    Q_gshp = Q_para[1]*x[8]
    Q_ct = Q_para[2]*Q_para[3]*x[12]*(x[20]-x[21])+Q_para[4]*Q_para[5]*(Q_para[6]-x[20])
    Q_d = 1000
    return Q_ashp + Q_gshp -Q_ct - Q_d

def w_cons(x):
    #生活热水约束
    return

def ct_cons(x):
    #热水罐约束
    return

def hy_incons1(x):
    #氢系统约束
    E_hy_para = [1, 15]  # 1/beta_el,beta_fc
    return x[13] - E_hy_para[1]*x[17]

def hy_incons2(x):
    E_hy_para = [1, 15]  # 1/beta_el,beta_fc
    return x[14]-(x[15]+x[17]-x[16]-x[10])*E_hy_para[0]

#供热系统约束
def stc_constraint(x):
    return

#设定变量上下界约束
b = (1,np.inf)
bnds = (b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b)

# 设定初值
x0 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

'''
#设定约束
cons1 = {'type':'eq','fun':P_constraint }
cons2 = {'type': 'ineq','fun': [function_name]}
...

cons = [cons1 , con2,...]
#求解
res = minimize(object,x0,method='SLSQP',bounds = bnds,constraints=cons )

print(res)
'''