# # # pv
# #     pv_para = [0.75]  # gama
# #     p_pv_k = i_pv * s_pv * pv_para[0]
# #     # hy
# #     hy_para = [16.6, 1.399, 15, 1]  # b_fcg,b_elco,b_fce,t
# #     g_hyout_k = hy_para[0] * m_hyout_k
# #     # wtg
# #     p_wtg_k = fai_wt * p_wtg
# #     # stc
# #     stc_para = [0.75]  # b_stc
# #     g_stc_k = i_pv * s_stc * stc_para[0]
# #     # gshp
# #     gshp_para = [4, 5, 1]  # gshp_g,gshp_q,t
# #     g_gshp_k = gshp_para[0] * p_gshpG_k * gshp_para[2]
# #     q_gshp_k = gshp_para[1] * p_gshpQ_k * gshp_para[2]
# #     # bhp
# #     bhp_para = [0.75, 1.2]  # 1-1/n_G,1+1/n_Q
# #     g_bhp_gshp_k = g_gshp_k * bhp_para[0]
# #     g_gshp_bhp_k = g_gshp_k * bhp_para[1]
# #     # bio
# #     bio_para = [15000, 0.76]  # LHV,1-n_bio
# #     g_bio_k = m_bio_k * bio_para[0]
# #     g_bio_out_k = g_bio_k * bio_para[1]
# #     # ashp
# #     ashp_para = [4, 5, 1]  # n_ashp_g,n_ashp_q,t
# #     g_ashp_k = p_ashpG_k * ashp_para[0] * ashp_para[2]
# #     q_ashp_k = p_ashpQ_k * ashp_para[1] * ashp_para[2]
# #     # eb
# #     eb_para = [0.97, 1]  # b_eb,t
# #     g_eb_k = p_eb_k * eb_para[0] *eb_para[1]
# #     # gfb
# #     gfb_para = [0.9]  # b_gfb
# #     v_gfb_k = g_gfb_k / gfb_para[0]
# #     # cfb
# #     cfb_para = [0.84]  # b_cfb
# #     m_cfb_k = g_cfb_k / cfb_para[0]
# #     # ofb
# #     ofb_para = [0.9]  # b_ofb
# #     v_ofb_k = g_ofb_k / ofb_para[0]
# #     # ru
# #     ru_para = [2.5]  # n_ru
# #     q_ru_k = p_ru_k * ru_para[0]
# #     # wt
# #     wt_para = [0.01, 1000, 4.2, 85, 45, 21, 4]  # u,p,C,T_WtMax,T_WtMin,T_CtMax,T_CtMin
#
# CE_ashp = p_ashp * 2.83 / 15 * 2.88  # y_e,ashp为15年，GWP为2880
#     CE_stc = 25.2 * 10 ** (-6) * s_stc
#     CE_bb = 1788 * 10 ** (-6) * m_bio
#     CE_eb = 650 * 10 ** (-6) * p_eb
#     CE_pv = 34 * 10 ** (-6) * s_pv
#     CE_cfb = 876 * 10 ** (-6) * g_cfb
#     CE_gfb = 412 * 10 ** (-6) * g_gfb
#     CE_ofb = 579 * 10 ** (-6) * g_ofb
#     CE_gshp = 79 * 10 ** (-6) * p_gshp
#     CE_wtg = 17 * 10 ** (-6) * p_wtg
#     CE_hy = 51 * 10 ** (-6) * p_fc
#     CE_wt = 500 * 10 ** (-6) * v_wt
#     CE_ct = 500 * 10 ** (-6) * v_ct
#     CE_ru = 2000 * 10 ** (-6) * p_ru
#     CE_r = CE_ashp + CE_ru + CE_gshp
#     ce_es = CE_ashp + CE_stc + CE_bb + CE_eb + CE_pv + CE_cfb + CE_gfb + CE_ofb + CE_gshp + CE_wtg + CE_hy + CE_wt + CE_ct + CE_ru + CE_r
#     ce_cpesr = gp.quicksum([p_g_k[i]*0.6101 + m_cfb_k[i]*8.14 + v_gfb_k[i]*1.535 + v_ofb_k[i]*2.25 for i in range(period)])
#     ce_total = ce_an + ce_es + ce_cpesr
#     ## set c_invest
#     i_crf, m_crf = 0.08, 15
#     CRF = i_crf * math.pow(1 + i_crf, m_crf) / (math.pow(1 + i_crf, m_crf) - 1)
#     c_invest = 3000*p_ashp + 1500*s_stc + 360*m_bio + 160*p_eb + 6600*s_pv + 3980*g_gfb + 3980*g_ofb + 3980*g_cfb + 400*p_gshp + 3600*p_wtg + 10000*p_fc + 10000*p_el + 1000*p_ru + 3000*m_hy + 500*v_wt +500*v_ct
#     ## set c_om
#     # om_para = [0.58, 0.8, 1.2, 3.5, 3.75, 40, 0.02]  各能源实时价格参数表[lep,bio,coil,gas,oil,hy,sigma]
#     c_om1 = gp.quicksum([p_g_k[i]*0.58 + m_bio_k[i]*0.8 + m_cfb_k[i]*1.2 + 3.5*v_gfb_k[i] + 3.75*v_ofb_k[i] + 40*m_hy_purIn_k[i] for i in range(period)])
#     c_om = c_om1+ 0.02*c_invest
#     ## set c_re
#     c_re = c_invest*0.02
#     ## set zeta
#     zeta = 0.1*300*s_pv