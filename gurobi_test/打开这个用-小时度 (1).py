import pulp
import sys
import time
import xlrd
import xlwt
import numpy as np
MST_HT=20000
P_EL=250000
lambda_ele_in = [0.4096, 0.4096, 0.4096, 0.4096, 0.4096, 0.4096, 0.4096, 0.7704, 1.1313, 1.1313, 1.1313, 1.1313,
                 0.7704, 0.7704, 0.7704, 0.7704, 0.7704, 0.7704, 1.1313, 1.1313, 1.1313, 1.1313, 1.1313, 0.4096] * 7
BETA_g=1.0246#KGCO2/KWh
BETA_h2=0
BETA_r=110#kgco2/GJ
GAMA_g= 0.1229#kgce/kWh
GAMA_h2=4.099#0.3685kgce/m3,当量折标系数+氢气的密度为0.0899 kg/m3，即4.099kgce/kg
lambda_carbon = 0.06
alpha_ele = 1.01
M=100000000
P_steel=464#KWh
P_dq=46.48#KW
G_steel=0.167#GJ
G_blast=0#GJ
G_conver=0#GJ
G_pr=0#GJ
C_cokeovengas=0#kj/m3
C_blastfurnacegas=0#kj/m3
C_convertergas=0#kj/m3
V_cokeovengas=88#m3
V_blastfurnacegas=357#m3
V_convertergas=57#m3
ETA_fc_e=0.3
ETA_fc_g=0.6
ETA_bpt_e=0.3
ETA_bpt_g=0.6
tao=3600
C_h2=142351#kJ/kg
BETA_pv=0.16
S_pv=200#m2
BETA_el=0.022
BETA_co = 0.05
M_INOUT=100000000000000
H_INOUT=100000000000000
#铁矿石，白云石，生石灰，压缩空气，球团矿，废钢，氧气，废钢/贴,焦炭粉、焦炭、无烟煤
M_J=[1192,16.7,71.7,74.5,74.1,59.1,39.8,229.3,34.3,290,124]
ALPHA_J=[0.77,93.2,121.36,0.1,0.94,3.4,92.32,3.4,1.8,1.8,1.4]
LAMDA_J=[0.0012,0.445,0.405,0.127,0.092,0.015,1.267,0.015,3.077,3.077,2.955]
R_J=np.multiply(np.array(M_J),np.array(ALPHA_J))
C_J=np.multiply(np.array(M_J),np.array(LAMDA_J))
R_J=np.sum(R_J)
C_J=np.sum(C_J)
#高炉渣，转炉/电炉渣
M_O=[501,214]
ALPHA_O=[0.45,1.5]
LAMDA_O=[0.55,0.3]
R_O=np.multiply(np.array(M_O),np.array(ALPHA_O))
C_O=np.multiply(np.array(M_O),np.array(LAMDA_O))
R_O=np.sum(R_O)
C_O=np.sum(C_O)
#焦炭粉、焦炭、无烟煤
M_F=[34.3,290,124]
GAMA_F=[0.9714 ,0.9714,0.9428]
E_F=np.multiply(np.array(M_F),np.array(GAMA_F))
E_F=np.sum(E_F)

def operating_problem(dict,period,lambda_h2,lambda_co2, steel_crease_rate):
    y_steel0 = dict['y_steel']
    r_solar = dict['r_solar']
    print(y_steel0)
    print(r_solar)
    p_d=np.zeros(period)
    g_d1=np.zeros(period)
    g_d2_mid=np.zeros(period)
    g_d2=np.zeros(period)
    h_ge=np.zeros(period)
    c_nomald=np.zeros(period)
    print(y_steel0[0])
    print(P_steel - P_dq)
    for i in range(period):
        y_steel[i] = y_steel0[i] * (steel_crease_rate)
        p_d[i] = (P_steel - P_dq) * y_steel[i]
        g_d1[i] = G_steel * y_steel[i]
        g_d2_mid[i] = (G_blast + G_conver - G_pr) * y_steel[i]
        h_ge[i] = (C_cokeovengas * V_cokeovengas + C_blastfurnacegas * V_blastfurnacegas + C_convertergas * V_convertergas) * y_steel[i] * 0.5
        c_nomald[i] = BETA_g * p_d[i] + BETA_r * 0.167 * y_steel[i]
    r_product = (R_J - R_O) * np.sum(y_steel)
    c_product = (C_J - C_O) * np.sum(y_steel)
    e_product = (E_F) * np.sum(y_steel)
    c_nomal=np.sum(c_nomald)
    prob = pulp.LpProblem('Operating Problem', sense=pulp.LpMinimize)
    # p_d = [pulp.LpVariable(f'p_d{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # g_d1 = [pulp.LpVariable(f'g_d1{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # g_d2 = [pulp.LpVariable(f'g_d2{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # g_d2_mid = [pulp.LpVariable(f'g_d2_mid{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # z_gd21 = [pulp.LpVariable(f'z_gd21{i}', cat = pulp.LpBinary) for i in range(period)]
    # z_gd22 = [pulp.LpVariable(f'z_gd22{i}', cat = pulp.LpBinary) for i in range(period)]
    # mid_gd21 = [pulp.LpVariable(f'mid_gd21{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # mid_gd22 = [pulp.LpVariable(f'mid_gd22{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # fin_gd21 = [pulp.LpVariable(f'fin_gd21{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # fin_gd22 = [pulp.LpVariable(f'fin_gd22{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # h_ge = [pulp.LpVariable(f'h_ge{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    # r_product = pulp.LpVariable(f'r_product', lowBound=None, upBound=None, cat=pulp.LpContinuous)
    # c_product = pulp.LpVariable(f'c_product', lowBound=None, upBound=None, cat=pulp.LpContinuous)
    # e_product = pulp.LpVariable(f'e_product', lowBound=None, upBound=None, cat=pulp.LpContinuous)
    p_fc = [pulp.LpVariable(f'p_fc{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    g_fc = [pulp.LpVariable(f'g_fc{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    m_fc = [pulp.LpVariable(f'm_fc{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    h_ge_fc = [pulp.LpVariable(f'h_ge_fc{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    h_in = [pulp.LpVariable(f'h_in{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    h_out = [pulp.LpVariable(f'h_out{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    hs_ge = [pulp.LpVariable(f'hs_ge{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period+1)]
    z_hin = [pulp.LpVariable(f'z_hin{i}', cat=pulp.LpBinary) for i in range(period)]
    z_hout = [pulp.LpVariable(f'z_hout{i}', cat=pulp.LpBinary) for i in range(period)]
    p_bpt = [pulp.LpVariable(f'p_bpt{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    g_bpt = [pulp.LpVariable(f'g_bpt{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    g_fc_bpt = [pulp.LpVariable(f'g_fc_bpt{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    p_pv = [pulp.LpVariable(f'p_pv{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    p_g = [pulp.LpVariable(f'p_g{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    p_el = [pulp.LpVariable(f'p_el{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    p_co = [pulp.LpVariable(f'p_co{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    m_b = [pulp.LpVariable(f'm_b{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    m_out = [pulp.LpVariable(f'm_out{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    m_in = [pulp.LpVariable(f'm_in{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    m_el = [pulp.LpVariable(f'm_el{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period)]
    z_in = [pulp.LpVariable(f'z_in{i}', cat = pulp.LpBinary) for i in range(period)]
    z_out = [pulp.LpVariable(f'z_out{i}', cat = pulp.LpBinary) for i in range(period)]
    m_ht = [pulp.LpVariable(f'm_ht{i}', lowBound = 0, upBound = None, cat = pulp.LpContinuous) for i in range(period+1)]
    r_enegry = pulp.LpVariable(f'r_enegry', lowBound = None, upBound = None, cat = pulp.LpContinuous)
    c_enegry = pulp.LpVariable(f'c_enegry', lowBound = None, upBound = None, cat = pulp.LpContinuous)
    e_enegry = pulp.LpVariable(f'e_enegry', lowBound = None, upBound = None, cat = pulp.LpContinuous)
    r_all = pulp.LpVariable(f'r_all', lowBound = None, upBound = None, cat = pulp.LpContinuous)
    c_all = pulp.LpVariable(f'c_all', lowBound = None, upBound = None, cat = pulp.LpContinuous)
    e_all = pulp.LpVariable(f'e_all', lowBound = None, upBound = None, cat = pulp.LpContinuous)

    prob += (r_all - lambda_co2 * (c_nomal - c_enegry))

    for i in range(period):
        # prob += (p_d[i] == (P_steel - P_dq) * y_steel[i])
        # prob += (g_d1[i] == G_coking * y_steel[i])
        # prob += (g_d2_mid[i] == (G_blast + G_conver - G_pr) * y_steel[i])
        # prob += (z_gd21[i] + z_gd22[i] <= 1)
        # prob += (mid_gd21[i] >= z_gd21[i] * -M)
        # prob += (mid_gd21[i] <= z_gd21[i] * 0)
        # prob += (mid_gd22[i] >= z_gd22[i] * 0)
        # prob += (mid_gd22[i] <= z_gd22[i] * M)
        # prob += (g_d2_mid[i] == mid_gd21[i]  + mid_gd22[i])
        # prob += (g_d2[i] == fin_gd21[i]  + fin_gd22[i])
        # prob += (g_d2_mid[i] - M * (1 - z_gd21[i]) <= fin_gd21[i])
        # prob += (fin_gd21[i] <= g_d2_mid[i] + M * (1 - z_gd21[i]))
        # prob += (fin_gd21[i] <= M * z_gd21[i])
        # prob += (g_d2_mid[i] - M * (1 - z_gd22[i]) <= fin_gd22[i])
        # prob += (fin_gd22[i] <= g_d2_mid[i] + M * (1 - z_gd22[i]))
        # prob += (fin_gd22[i] <= M * z_gd22[i])
        # prob += (h_ge[i] == (C_cokeovengas * V_cokeovengas + C_blastfurnacegas * V_blastfurnacegas + C_convertergas * V_convertergas) * y_steel[i])


        prob += (p_fc[i] * tao == ETA_fc_e * (h_ge[i] + m_fc[i] * C_h2))#kj==kj
        prob += (g_fc[i]  == ETA_fc_g * (h_ge[i] + m_fc[i] * C_h2))#kj==kj
        prob += (p_bpt[i] * tao == ETA_bpt_e * g_fc_bpt[i])#kj==kj
        prob += (g_bpt[i]  == ETA_bpt_g * g_fc_bpt[i])#kj==kj
        prob += (g_fc[i] - g_fc_bpt[i] >= g_d1[i]*1000000)#kj==kj
        prob += (g_bpt[i] >= g_d2[i]*1000000)#kj==kj
        prob += (p_pv[i] == BETA_pv * S_pv * r_solar[i])
        prob += (p_fc[i] + p_bpt[i]+ p_g[i] + p_pv[i] == p_d[i] + p_el[i]+ p_co[i])
        prob += (m_b[i] + m_out[i] + m_el[i] == m_in[i] + m_fc[i])
        prob += (m_el[i] == BETA_el * p_el[i])
        prob += (p_co[i] == BETA_co * m_in[i])
        prob += (z_in[i] + z_out[i] <= 1)
        prob += (m_in[i] >= 0)
        prob += (m_in[i] <= z_in[i] * M_INOUT)
        prob += (m_out[i] >= 0)
        prob += (m_out[i] <= z_out[i] * M_INOUT)
        prob += (m_ht[i+1] == m_ht[i] + m_in[i] - m_out[i])
        prob += (h_ge[i] + h_out[i] == h_ge_fc[i] + h_in[i])
        prob += (z_hin[i] + z_hout[i] <= 1)
        prob += (h_in[i] >= 0)
        prob += (h_in[i] <= z_hin[i] * H_INOUT)
        prob += (h_out[i] >= 0)
        prob += (h_out[i] <= z_hout[i] * H_INOUT)
        prob += (hs_ge[i + 1] == hs_ge[i] + h_in[i] - h_out[i])
        prob += (m_ht[i+1] <= MST_HT)
        prob += (p_el[i] <= P_EL)
    prob += (r_enegry == pulp.lpDot(p_g, lambda_ele_in) + pulp.lpSum(m_b)*lambda_h2)
    prob += (c_enegry == BETA_g * pulp.lpSum(p_g) + BETA_h2 * pulp.lpSum(m_b))
    prob += (e_enegry == GAMA_g * pulp.lpSum(p_g) + GAMA_h2 * pulp.lpSum(m_b))
    prob += (r_all == r_product + r_enegry)
    prob += (c_all == c_product + c_enegry)
    prob += (e_all == e_product + e_enegry)
    prob += (m_ht[period] == m_ht[0])
    prob += (hs_ge[period] == hs_ge[0])


    print("It has reached here")
    prob.solve(pulp.GUROBI_CMD(options=[('MIPgap',0.0001)]))
    print('status:', pulp.LpStatus[prob.status])
    return {'objective':pulp.value(prob.objective),
            '电需求': [p_d[i] for i in range(period)],
            '一级热需求': [g_d1[i] for i in range(period)],
            '二级热需求': [g_d2[i] for i in range(period)],
            '煤气热值': [h_ge[i] for i in range(period)],
            'p_fc': [pulp.value(p_fc[i]) for i in range(period)],
            'p_bpt': [pulp.value(p_bpt[i]) for i in range(period)],
            'p_g': [pulp.value(p_g[i]) for i in range(period)],
            'p_pv': [pulp.value(p_pv[i]) for i in range(period)],
            'p_d': [pulp.value(p_d[i]) for i in range(period)],
            'p_el': [pulp.value(p_el[i]) for i in range(period)],
            'p_co': [pulp.value(p_co[i]) for i in range(period)],
            'g_fc': [pulp.value(g_fc[i]) for i in range(period)],
            'g_fc_bpt': [pulp.value(g_fc_bpt[i]) for i in range(period)],
            'g_d1': [pulp.value(g_d1[i]) for i in range(period)],
            'g_bpt': [pulp.value(g_bpt[i]) for i in range(period)],
            'g_d2': [pulp.value(g_d2[i]) for i in range(period)],
            'm_b': [pulp.value(m_b[i]) for i in range(period)],
            'm_out': [pulp.value(m_out[i]) for i in range(period)],
            'm_el': [pulp.value(m_el[i]) for i in range(period)],
            'm_in': [pulp.value(m_in[i]) for i in range(period)],
            'm_fc': [pulp.value(m_fc[i]) for i in range(period)],
            'h_ge': [pulp.value(h_ge[i]) for i in range(period)],
            'h_out': [pulp.value(h_out[i]) for i in range(period)],
            'h_ge_fc': [pulp.value(h_ge_fc[i]) for i in range(period)],
            'h_in': [pulp.value(h_in[i]) for i in range(period)],
            'm_ht': [pulp.value(m_ht[i]) for i in range(period + 1)],
            'hs_ge': [pulp.value(hs_ge[i]) for i in range(period+1)],
            '生产价格': r_product,
            '生产碳排': c_product,
            '生产能耗': e_product,
            '能源价格': pulp.value(r_enegry),
            '能源碳排': pulp.value(c_enegry),
            '能源能耗': pulp.value(e_enegry),
            '总价格': pulp.value(r_enegry) + r_product,
            '总碳排': pulp.value(c_enegry) + c_product,
            '总能耗': pulp.value(e_enegry) + e_product,
            '日常能源碳排': pulp.value(c_nomal),
            '氢气-买电供电比例':sum((pulp.value(p_fc[i]) + pulp.value(p_bpt[i]))for i in range(period))/sum(pulp.value(p_d[i])for i in range(period)),
            '买氢买电比例': sum((pulp.value(m_b[i]) * lambda_h2)for i in range(period))/pulp.value(r_enegry),
            }
if __name__ == '__main__':
    dataf = xlwt.Workbook()
    dd = dataf.add_sheet('miao')
    dd.write(0, 0, '生产价格')
    dd.write(0, 1, '生产碳排')
    dd.write(0, 2, '生产能耗')
    dd.write(0, 3, '能源价格')
    dd.write(0, 4, '能源碳排')
    dd.write(0, 5, '能源能耗')
    dd.write(0, 6, '总价格')
    dd.write(0, 7, '总碳排')
    dd.write(0, 8, '总能耗')
    dd.write(0, 9, '日常能源碳排')
    dd.write(0, 10, '燃气供电比例')
    dd.write(0, 11, '买氢买电比例')
    dd.write(0, 12, '氢价')
    dd.write(0, 13, '碳价')
    dd.write(0, 14, '产量增量幅度')



    period = 24
    for k in range(10):
        book = xlrd.open_workbook('per.xlsx')
        data = book.sheet_by_index(0)
        pricebook = xlrd.open_workbook('price.xlsx')
        price = pricebook.sheet_by_index(0)
        y_steel = []
        r_solar = []
        dict = {'y_steel': y_steel, 'r_solar': r_solar}
        lambda_h2 = price.cell(k+1, 1).value
        lambda_co2 = price.cell(k+1, 2).value #元/kg
        steel_crease_rate = price.cell(k+1, 3).value / price.cell(1, 3).value
        for l in range(period):
            y_steel.append(data.cell(l, 0).value)
            r_solar.append(data.cell(l, 1).value)
        res = operating_problem(dict, period, lambda_h2, lambda_co2, steel_crease_rate)
        items = list(res.keys())
        print(items[30])
        wb = xlwt.Workbook()
        total = wb.add_sheet('case')
        for i in range(len(items)):
            total.write(0,i,items[i])
            if type(res[items[i]]) == list:
                for j in range(period):
                    total.write(j+2,i,(res[items[i]])[j])
            else:
                total.write(1,i,res[items[i]])
        filename = '氢价' + str(lambda_h2 )+ '碳价' + str(lambda_co2) + '.xls'
        wb.save(filename)
        dd.write(k+1, 12, lambda_h2)
        dd.write(k+1, 13, lambda_co2)
        dd.write(k+1, 14, steel_crease_rate)
        for j in range(12):
            dd.write(k+1,j,res[items[j+28]])
    filename1 = '数据_正常预期碳价' + '.xls'
    dataf.save(filename1)
    print("finished")
