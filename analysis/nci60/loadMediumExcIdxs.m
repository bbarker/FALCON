function mediumExcIdxs=loadMediumExcIdxs(model)
mediumExc={'EX_gly(e)','EX_arg_L(e)','EX_asp_L(e)','EX_asn_L(e)','EX_cys_L(e)','EX_glu_L(e)','EX_gln_L(e)','EX_his_L(e)','EX_4hpro(e)','EX_ile_L(e)','EX_leu_L(e)',...
'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)',...
'EX_fol(e)','EX_ncam(e)','EX_bz(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_adpcbl(e)','EX_inost(e)','EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)',...
'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)','EX_co2(e)','EX_h2o(e)','EX_h(e)'};
mediumExcIdxs=[];
for i=1:length(mediumExc)
    mediumExcIdx=find(ismember(model.rxns,mediumExc{i}));
    mediumExcIdxs(end+1)=mediumExcIdx;
end