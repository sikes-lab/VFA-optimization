function Conc_LoD = rcSso_NP_Saliva_LoD_variation_Solver(Visual_LoD,Sample_Conc_Guess,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt)

% This function takes an input Visual_LoD and operating parameters to
% determine the concentration corresponding to the Visual_LoD
options = optimset('TolFun',1e-16,'TolX',1e-16);

Conc_LoD = fzero(@(Sample_Conc) rcSso_NP_Saliva_LoD_variation_fun(Visual_LoD,Sample_Conc,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt),Sample_Conc_Guess,options);

end