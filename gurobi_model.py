import gurobipy as gp
from gurobipy import GRB
import sys
import numpy as np
import pandas as pd

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
    m_hy_purIn = [m.addVar(name=f"m_hy_purIn{i}") for i in range(period)]

    #pv
    pv_para = [0.75]  # gama
    p_pv_k = i_pv * s_pv * pv_para[0]
    #hy
    hy_para = [16.6, 1.399, 15, 1]  # b_fcg,b_elco,b_fce,t
    g_hyout_k = hy_para[0] * m_hyout_k
    #wtg
    p_wtg_k = fai_wt * p_wtg
    # stc
    stc_para = [0.75]  # b_stc
    g_stc_k = i_pv * s_stc * stc_para[0]
    #gshp
    gshp_para = [4, 5, 1]  # gshp_g,gshp_q,t
    g_gshp_k = gshp_para[0] * p_gshpG_k * gshp_para[2]
    q_gshp_k = gshp_para[1] * p_gshpQ_k * gshp_para[2]
    # bhp
    bhp_para = [0.75, 1.2]  # 1-1/n_G,1+1/n_Q
    g_bhp_gshp_k = g_gshp_k * bhp_para[0]
    g_gshp_bhp_k = g_gshp_k * bhp_para[1]
    #bio
    bio_para = [15000, 0.76]  # LHV,1-n_bio
    g_bio_k = m_bio_k * bio_para[0]
    g_bio_out_k = g_bio_k * bio_para[1]
    # ashp
    ashp_para = [4, 5, 1]  # n_ashp_g,n_ashp_q,t
    g_ashp_k = p_ashpG_k * ashp_para[0] * ashp_para[2]
    q_ashp_k = p_ashpQ_k * ashp_para[1] * ashp_para[2]
    #eb
    eb_para = [0.97, 1]  # b_eb,t
    g_eb_k = p_eb_k * eb_para[0] *eb_para[1]
    #gfb
    gfb_para = [0.9]  # b_gfb
    v_gfb_k = g_gfb_k / gfb_para[0]
    # cfb
    cfb_para = [0.84]  # b_cfb
    m_cfb_k = g_cfb_k / cfb_para[0]
    # ofb
    ofb_para = [0.9]  # b_ofb
    v_ofb_k = g_ofb_k / ofb_para[0]
    #ru
    ru_para = [2.5]  # n_ru
    q_ru_k = p_ru_k * ru_para[0]
    # wt
    wt_para = [0.01, 1000, 4.2, 85, 45, 21, 4]  # u,p,C,T_WtMax,T_WtMin,T_CtMax,T_CtMin
    # Add constraint
    for i in range(period):
        m.addConstr(p_g_k[i] + p_pv_k[i] + p_wtg_k[i] + p_hyout_k[i] == p_d[i] + p_ashpQ_k[i] + p_ashpG_k[i] + p_eb_k[i] + p_gshp_k[i] + p_hyin_k[i] + p_ru_k[i])
        m.addConstr(g_ashp_k[i] + g_stc_k[i] + g_bio_out_k[i] + g_eb_k[i] + g_hyout_k[i] + g_cfb_k[i] + g_gfb_k[i] + g_ofb_k[i] + g_gshp_k[i] + g_wt_out_k[i] == g_d[i] + g_w[i] +g_bhp_in_k[i] + g_wt_in_k[i] + g_wt_k[i]*wt_para[0])
