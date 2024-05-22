function ret = rcSso_NP_Saliva_LoD_variation_fun(Visual_LoD,Sample_Conc,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt)
load('White_Saliva_HRP_Calibration_beta_Cleaned.mat'); %load betas
MW_HRP = Soln_pptys_in(1,3)*1000; %pick a random one to use

Model_output = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only(Sample_Conc,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt);
Cyan_out = sum(Model_output);

Cyan_out = beta_all(1).*(1-exp(-beta_all(2).*Cyan_out)); % Change to cyan intensity

ret = Cyan_out - Visual_LoD; %Objective function

end