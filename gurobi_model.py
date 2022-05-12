from typing import List

import gurobipy as gp
from gurobipy import GRB
import sys
import numpy as np
import pandas as pd
import math


def plan_carbon_problem(ce_an,p_d,g_d,g_w,q_d,i_pv,fai_wt,period):
    # 国标约束
    building_area = 156027.93  # A
    CE_2016 = 35 * building_area
    A = 156027.93

    # Create a new model
    m = gp.Model("carbon_model")

    # Create variables
    p_ashp = m.addVar(name='p_ashp')
    s_stc = m.addVar(name='s_stc')
    m_bio = m.addVar(name='m_bio')
    p_eb = m.addVar(name='p_eb')
    s_pv = m.addVar(name='s_pv')
    g_gfb = m.addVar(name='g_gfb')
    g_cfb = m.addVar(name='g_cfb')
    g_ofb = m.addVar(name='g_ofb')
    p_gshp = m.addVar(name='p_gshp')
    p_wtg = m.addVar(name='p_wtg')
    m_hy = m.addVar(name='m_hy')
    v_wt = m.addVar(name='v_wt')
    v_ct = m.addVar(name='v_ct')
    p_fc = m.addVar(name='p_fc')
    e_g = m.addVar(name='e_g')
    hy_in = m.addVar(name='hy_in')
    hy_out = m.addVar(name='hy_out')
    t_wt = m.addVar(name='t_wt')
    t_ct = m.addVar(name='t_ct')
    p_ashpQ = m.addVar(name='p_ashpQ')
    p_gshpQ = m.addVar(name='p_gshpQ')
    p_ru = m.addVar(name='p_ru')
    p_el = m.addVar(name='p_el')

    p_g_k = [m.addVar(name=f"p_g_k{i}") for i in range(period)]
    g_gfb_k = [m.addVar(name=f"g_gfb_k{i}") for i in range(period)]
    g_cfb_k = [m.addVar(name=f"g_cfb_k{i}") for i in range(period)]
    g_ofb_k = [m.addVar(name=f"g_ofb_k{i}") for i in range(period)]
    m_hyout_k = [m.addVar(name=f"m_hyout_k{i}") for i in range(period)]
    m_hy_k = [m.addVar(name=f"m_hy_k{i}") for i in range(period)]
    p_ashpG_k = [m.addVar(name=f"p_ashpG_k{i}") for i in range(period)]
    p_eb_k = [m.addVar(name=f"p_eb_k{i}") for i in range(period)]
    p_gshp_k = [m.addVar(name=f"p_gshp_k{i}") for i in range(period)]
    p_ru_k = [m.addVar(name=f"p_ru_k{i}") for i in range(period)]
    p_hyin_k = [m.addVar(name=f"p_hyin_k{i}") for i in range(period)]
    p_hyout_k = [m.addVar(name=f"p_hyout_k{i}") for i in range(period)]
    m_bio_k = [m.addVar(name=f"m_bio_k{i}") for i in range(period)]
    g_wt_k = [m.addVar(name=f"g_wt_k{i}") for i in range(period)]
    g_wt_in_k = [m.addVar(name=f"g_wt_in{i}") for i in range(period)]
    g_wt_out_k = [m.addVar(name=f"g_wt_out{i}") for i in range(period)]
    q_ct_k = [m.addVar(name=f"q_ct_k{i}") for i in range(period)]
    q_ct_in_k = [m.addVar(name=f"q_ct_in{i}") for i in range(period)]
    q_ct_out_k = [m.addVar(name=f"q_ct_out{i}") for i in range(period)]
    p_gshpG_k = [m.addVar(name=f"p_gshpG_k{i}") for i in range(period)]
    p_gshpQ_k = [m.addVar(name=f"p_gshpQ_k{i}") for i in range(period)]
    g_bhp_k = [m.addVar(name=f"g_bhp_k{i}") for i in range(period)]
    g_bhp_in_k = [m.addVar(name=f"g_bhp_in_k{i}") for i in range(period)]
    p_ashpQ_k = [m.addVar(name=f"p_ashpQ_k{i}") for i in range(period)]
    m_hy_purIn_k = [m.addVar(name=f"m_hy_purIn{i}") for i in range(period)]
    p_pv_k = [m.addVar(name=f"p_pv_k{i}") for i in range(period)]
    g_hyout_k = [m.addVar(name=f"g_hyout_k{i}") for i in range(period)]
    p_wtg_k = [m.addVar(name=f"p_wtg_k{i}") for i in range(period)]
    g_stc_k = [m.addVar(name=f"g_stc_k{i}") for i in range(period)]
    g_gshp_k = [m.addVar(name=f"g_gshp_k{i}") for i in range(period)]
    q_gshp_k = [m.addVar(name=f"q_gshp_k{i}") for i in range(period)]
    g_bhp_gshp_k = [m.addVar(name=f"g_bhp_gshp_k{i}") for i in range(period)]
    g_gshp_bhp_k = [m.addVar(name=f"g_gshp_bhp_k{i}") for i in range(period)]
    g_bio_out_k = [m.addVar(name=f"g_bio_k{i}") for i in range(period)]
    g_ashp_k= [m.addVar(name=f"g_ashp_k{i}") for i in range(period)]
    q_ashp_k = [m.addVar(name=f"q_ashp_k{i}") for i in range(period)]
    g_eb_k = [m.addVar(name=f"g_eb_k{i}") for i in range(period)]
    v_gfb_k = [m.addVar(name=f"v_gfb_k{i}") for i in range(period)]
    m_cfb_k = [m.addVar(name=f"m_cfb_k{i}") for i in range(period)]
    v_ofb_k = [m.addVar(name=f"v_ofb_k{i}") for i in range(period)]
    q_ru_k = [m.addVar(name=f"q_ru_k{i}") for i in range(period)]
    # pv
    pv_para = [0.75]  # gama
    p_pv_k = i_pv * s_pv * pv_para[0]
    # hy
    hy_para = [16.6, 1.399, 15, 1]  # b_fcg,b_elco,b_fce,t
    g_hyout_k = hy_para[0] * m_hyout_k
    # wtg
    p_wtg_k = fai_wt * p_wtg
    # stc
    stc_para = [0.75]  # b_stc
    g_stc_k = i_pv * s_stc * stc_para[0]
    # gshp
    gshp_para = [4, 5, 1]  # gshp_g,gshp_q,t
    g_gshp_k = gshp_para[0] * p_gshpG_k * gshp_para[2]
    q_gshp_k = gshp_para[1] * p_gshpQ_k * gshp_para[2]
    # bhp
    bhp_para = [0.75, 1.2]  # 1-1/n_G,1+1/n_Q
    g_bhp_gshp_k = g_gshp_k * bhp_para[0]
    g_gshp_bhp_k = g_gshp_k * bhp_para[1]
    # bio
    bio_para = [15000, 0.76]  # LHV,1-n_bio
    g_bio_k = m_bio_k * bio_para[0]
    g_bio_out_k = g_bio_k * bio_para[1]
    # ashp
    ashp_para = [4, 5, 1]  # n_ashp_g,n_ashp_q,t
    g_ashp_k = p_ashpG_k * ashp_para[0] * ashp_para[2]
    q_ashp_k = p_ashpQ_k * ashp_para[1] * ashp_para[2]
    # eb
    eb_para = [0.97, 1]  # b_eb,t
    g_eb_k = p_eb_k * eb_para[0] *eb_para[1]
    # gfb
    gfb_para = [0.9]  # b_gfb
    v_gfb_k = g_gfb_k / gfb_para[0]
    # cfb
    cfb_para = [0.84]  # b_cfb
    m_cfb_k = g_cfb_k / cfb_para[0]
    # ofb
    ofb_para = [0.9]  # b_ofb
    v_ofb_k = g_ofb_k / ofb_para[0]
    # ru
    ru_para = [2.5]  # n_ru
    q_ru_k = p_ru_k * ru_para[0]
    # wt
    wt_para = [0.01, 1000, 4.2, 85, 45, 21, 4]  # u,p,C,T_WtMax,T_WtMin,T_CtMax,T_CtMin
    # Add constraint
    for i in range(period):
        m.addConstr(p_g_k[i] + p_pv_k[i] + p_wtg_k[i] + p_hyout_k[i] == p_d[i] + p_ashpQ_k[i] + p_ashpG_k[i] + p_eb_k[i] + p_gshp_k[i] + p_hyin_k[i] + p_ru_k[i])
        m.addConstr(g_ashp_k[i] + g_stc_k[i] + g_bio_out_k[i] + g_eb_k[i] + g_hyout_k[i] + g_cfb_k[i] + g_gfb_k[i] + g_ofb_k[i] + g_gshp_k[i] + g_wt_out_k[i] == g_d[i] + g_w[i] +g_bhp_in_k[i] + g_wt_in_k[i] + g_wt_k[i]*wt_para[0])
        m.addConstr(q_ashp_k[i] + q_gshp_k[i] + q_ct_out_k[i] + q_ru_k[i] == q_d[i] + q_ct_in_k[i] + q_ct_k[i]*wt_para[0])
        m.addConstr(p_hyout_k[i]/hy_para[3] <= p_fc)
        m.addConstr(p_hyin_k[i]/hy_para[3] <= p_el)
        m.addConstr(m_hy_k[i] <= m_hy)
        m.addConstr(p_gshpQ_k/gshp_para[2] + p_gshpG_k/gshp_para[2] <= p_gshp)
        m.addConstr(p_gshpQ_k[i] + p_gshpG_k[i] == p_gshp_k[i])
        m.addConstr(m_bio_k[i] <= m_bio)
        m.addConstr(p_ashpG_k[i] <= p_ashp)
        m.addConstr(p_ashpQ_k[i] <= p_ashp)
        m.addConstr(p_eb_k[i] <= p_eb)
        m.addConstr(g_gfb_k[i] <= g_gfb)
        m.addConstr(g_cfb_k[i] <= g_cfb)
        m.addConstr(g_ofb_k[i] <= g_ofb)
        m.addConstr(p_ru_k[i] <= p_ru)
        m.addConstr(g_wt_k[i] <= v_wt*wt_para[1]*wt_para[2]*(wt_para[3] - wt_para[4]))
        m.addConstr(q_ct_k[i] <= v_ct * wt_para[1]*wt_para[2]*(wt_para[5]-wt_para[6]))
    for i in range(period+1):
        m.addConstr(m_hy_k[i+1] == m_hy_k[i] + m_hy_purIn_k[i] + hy_para[1]*p_hyin_k[i] - p_hyout_k[i]/hy_para[2])
        m.addConstr(g_bhp_k[i+1] - g_bhp_k[i] == g_gshp_bhp_k[i] - g_bhp_gshp_k[i] + g_bhp_in_k[i])
        m.addConstr(g_wt_k[i+1] - g_wt_k[i] == g_wt_in_k[i] - g_wt_out_k[i] - g_wt_k[i]*wt_para[0])
        m.addConstr(q_ct_k[i+1] - q_ct_k[i] == q_ct_in_k[i] - q_ct_out_k[i] - wt_para[0]*q_ct_k[i])

    # Set objective
    ## set CE_total
    CE_ashp = p_ashp * 2.83 / 15 * 2.88  # y_e,ashp为15年，GWP为2880
    CE_stc = 25.2 * 10 ** (-6) * s_stc
    CE_bb = 1788 * 10 ** (-6) * m_bio
    CE_eb = 650 * 10 ** (-6) * p_eb
    CE_pv = 34 * 10 ** (-6) * s_pv
    CE_cfb = 876 * 10 ** (-6) * g_cfb
    CE_gfb = 412 * 10 ** (-6) * g_gfb
    CE_ofb = 579 * 10 ** (-6) * g_ofb
    CE_gshp = 79 * 10 ** (-6) * p_gshp
    CE_wtg = 17 * 10 ** (-6) * p_wtg
    CE_hy = 51 * 10 ** (-6) * p_fc
    CE_wt = 500 * 10 ** (-6) * v_wt
    CE_ct = 500 * 10 ** (-6) * v_ct
    CE_ru = 2000 * 10 ** (-6) * p_ru
    CE_r = CE_ashp + CE_ru + CE_gshp
    ce_es = CE_ashp + CE_stc + CE_bb + CE_eb + CE_pv + CE_cfb + CE_gfb + CE_ofb + CE_gshp + CE_wtg + CE_hy + CE_wt + CE_ct + CE_ru + CE_r
    ce_cpesr = gp.quicksum([p_g_k[i]*0.6101 + m_cfb_k[i]*8.14 + v_gfb_k[i]*1.535 + v_ofb_k[i]*2.25 for i in range(period)])
    ce_total = ce_an + ce_es + ce_cpesr
    ## set c_invest
    i_crf, m_crf = 0.08, 15
    CRF = i_crf * math.pow(1 + i_crf, m_crf) / (math.pow(1 + i_crf, m_crf) - 1)
    c_invest = 3000*p_ashp + 1500*s_stc + 360*m_bio + 160*p_eb + 6600*s_pv + 3980*g_gfb + 3980*g_ofb + 3980*g_cfb + 400*p_gshp + 3600*p_wtg + 10000*p_fc + 10000*p_el + 1000*p_ru + 3000*m_hy + 500*v_wt +500*v_ct
    ## set c_om
    # om_para = [0.58, 0.8, 1.2, 3.5, 3.75, 40, 0.02]  各能源实时价格参数表[lep,bio,coil,gas,oil,hy,sigma]
    c_om = gp.quicksum([p_g_k[i]*0.58 + m_bio_k[i]*0.8 + m_cfb_k[i]*1.2 + 3.5*v_gfb_k[i] + 3.75*v_ofb_k[i] + 40*m_hy_purIn_k[i] for i in range(period)])
    c_om += 0.02*c_invest
    ## set c_re
    c_re = c_invest*0.02
    ## set zeta
    zeta = 0.1*300*s_pv
    m.setObjective(0.06*ce_total + CRF**c_invest + c_om - c_re -zeta,GRB.MINIMIZE)

    try:
        m.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")

    for v in m.getVars():
        print('%s %g' % (v.VarName, v.X))
    print('Obj: %g' % m.ObjVal)

if __name__ == '__main__':
    pv = pd.read_csv('pv_38.9869_106.3772.csv')
    i_pv = pv.iloc[:,2]

    wt = pd.read_csv('wind_38.4852_106.2262.csv')
    fai_wt = wt.iloc[:,2]

    info = pd.read_csv('p_d,g_d,q_d.csv',encoding_errors='ignore')
    p_d = info['p_d']
    p_d = p_d.fillna(0)

    g_d = info['g_d']
    g_d = g_d.fillna(0)

    q_d = info['q_d']
    q_d = q_d.fillna(0)

    g_w = info['g_w']
    g_w = g_w.fillna(0)

    plan_carbon_problem(16068.17, p_d, g_d, g_w, q_d, i_pv, fai_wt, 8760)
