function HRP_LoD = rcSso_NP_Saliva_Neg_Visual_LoD_Solver(Visual_LoD,HRP_Guess,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt)

% This function takes an input Visual_LoD and operating parameters to
% determine the concentration corresponding to the Visual_LoD when nothing
% is in the sample
options = optimset('TolFun',1e-12,'TolX',1e-12);
HRP_LoD = fzero(@(HRP) rcSso_NP_Saliva_LoD_Neg_fun(HRP,Visual_LoD,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt),[0 HRP_Guess*1.5],options);
%Specify bounds to search for LoD because increasing HRP can exceed the visual LoD

end

function ret = rcSso_NP_Saliva_LoD_Neg_fun(HRP,Visual_LoD,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt)
load('White_Saliva_HRP_Calibration_beta_Cleaned'); %load betas
MW_HRP = Soln_pptys_in(1,3)*1000; %pick a random one to use

Soln_pptys_in(1,9) = HRP.*1e-9.*1e6./(Soln_pptys_in(1,3).*1e3); %HRP is already in the desired units

Model_output = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only(0,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,param_pt);
Cyan_out = sum(Model_output);

Cyan_out = beta_all(1).*(1-exp(-beta_all(2).*Cyan_out)); % Change to cyan intensity

ret = Cyan_out - Visual_LoD; %Objective function

end