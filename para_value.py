#国标约束
building_area = 156027.93 #A
CE_2016 = 35 * building_area
A = 156027.93

#hy
hy_para = [16.6,1.399,15,1] # b_fcg,b_elco,b_fce,t

#stc
stc_para = [0.75]# b_stc

#gshp
gshp_para = [4,5,1]# gshp_g,gshp_q,t

#bhp
bhp_para = [0.75,1.2] # 1-1/n_G,1+1/n_Q

#bio
bio_para = [15000,0.76] # LHV,1-n_bio

#ashp
ashp_para = [4,5,1] # n_ashp_g,n_ashp_q,t

#eb
eb_para = [0.97,1] # b_eb,t

#gfb
gfb_para = [0.9]# b_gfb

#cfb
cfb_para = [0.84] #b_cfb

#ofb
ofb_para = [0.9] # b_ofb

#ru
ru_para = [2.5] # n_ru

#wt
wt_para = [0.01,1000,4.2,85,45,21,4] # u,p,C,T_WtMax,T_WtMin,T_CtMax,T_CtMin
