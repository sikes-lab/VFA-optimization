rcSso_NP_Swab_Sensitivity_variation_k_off(0.2) %Measure Sensitivity @ Cyan Intensity of 0.2

function rcSso_NP_Swab_Sensitivity_variation_k_off(Sensitivity_Location)

% This function plots the variation of LoD with changes in HRP_Conc, given
% a Visual_LoD

% Load relevant files
load('Sorted_NP_Saliva_Data_reorder.mat'); %data containing file
load('White_Saliva_HRP_Calibration_beta_Cleaned.mat'); %load parameters to convert pg on paper to Cyan Intensity
% Fitted graph is y = beta_conc(1)*(1-exp(-beta_conc(1)*m_HRP)) where m_HRP
% is in pg. y = Cyan Intensity @ 1 min

load('White_Saliva_Flow_Rates.mat'); %load flow rate data (key parameters: sample_flow_rate, avg_wash_flow_rate)
load('opt_NP_Saliva_params_all_1layer_Cellulose_NSB.mat')

% Place 'Premix_VFA_parfor_nHRP_only.m' & 'Polson_D_est.m' into the same
% folder as this file

%% Setup fixed parameters
k_offs = linspace(0.01,20,51).*1e-3; %[1/s]

%HRP_Concs = 2;

Soln_pptys = zeros(1,9);

Soln_pptys(1) = 46; %NP MW
Soln_pptys(2) = 53.6; %BA-MBP-Sso Construct
Soln_pptys(3) = 44+44+53; %Molecular weight of NA-HRP (Assume 2 HRPs, NA is 60 kDa, HRP is 44 kDa)
Target_MW = (Soln_pptys(1) + Soln_pptys(2) + Soln_pptys(3))*1e3; %MW of complex (assumed all to have similar diffusivities)

Soln_pptys(4) = (2.14e5).*1e-3; %[m3/(mol.s)] rcSso.BacNP2.1 (E2)
Soln_pptys(5) = (8.47e-4); %[1/s] rcSso.BacNP2.1 (E2) Use the one with His
%Tag (Refer to 20200616_NTU-SMART-MIT meeting)

Soln_pptys(6) = (2.03e5).*1e-3; %[m3/(mol.s)] rcSso.NP-PF2.1 (E1)
Soln_pptys(7) = (6.46e-4); %[1/s] rcSso.NP-PF2.1 (E1) Use the one with His Tag

Soln_pptys(8) = (50e-9)*1e3; %[mol/m3] Used 50 nM Rep in Sample solution
Soln_pptys(9) = (250e-12.*1e3*((100)))./(Soln_pptys(3)); %[mol/m3] Used 250 pM SA-HRP final concentration
% Note that 250 pM is ~25 ng/mL

Soln_pptys(10) = opt_params_glob(2);
Soln_pptys(11) = opt_params_glob(3);

Soln_pptys_in = repmat(Soln_pptys,length(k_offs),1);
Soln_pptys_in(:,7) = (k_offs.'); %convert to SI

% Convert from ugmL to mol/m3

Paper_params = zeros(1,7);
Paper_params(1) = Polson_D_est(Target_MW); %[m2/s] Diffusivity of complex estimated via Polson
Paper_params(2) = 5.5e-6; %[m] Radius of pore, taken at 11um diameter
Paper_params(3) = 180e-6; %[m] Thickness of cellulose filter paper
Paper_params(7) = 0.689; %Void fraction of cellulose paper from Madhu et al
Paper_params(4) = 1.25e-3; %[m] Radius of whole well - 2.5mm for bottom layers, 3.2mm for top layer, but for my sanity let us ignore the top layer first...

Paper_params(5) = (0.5e-6)*(30e-6); %[mol] CBD loaded - 30uM, 0.5uL
Paper_params(6) = opt_params_glob(1);

op_params = cell(1,7);

op_params{1,1} = [sample_flow_rate avg_wash_flow_rate].*1e-6.*1e-3; % [m3/s] 
op_params{1,2} = [20 50].*(1e-6)*1e-3; %[m3] 20 uL sample + 50 (25+25) uL wash

op_params{1,4} = "Premix";
op_params{1,5} = [1*60 25.5]; %Incubated premix for 15 min, followed by 15s with sample

op_params{1,6} = 2; %Premix is 2-fold that used in sample
op_params{1,7} = [1 0]; %Sample, then wash
% M = 26; N = 90; %Choose N such that mesh size is 2 um

% Setup matrices for parfor

Paper_params_in = repmat(Paper_params,length(k_offs),1);

op_params_in = cell(length(k_offs),3);
[op_params_in{:,1}] = deal(op_params{:,1});
[op_params_in{:,2}] = deal(op_params{:,2});

[op_params_in{:,4}] = deal(op_params{:,4});
[op_params_in{:,5}] = deal(op_params{:,5});
[op_params_in{:,6}] = deal(op_params{:,6});
[op_params_in{:,7}] = deal(op_params{:,7});

M = 26; N = 90; %Mesh sizes
InletBC = 'Danckwertz';

Sample_Conc_Guess = 10e-9.*1e3; % [10 nM, convert to nmol/m3]
Sens_All = zeros(size(k_offs));
Sens_Loc = zeros(size(k_offs));

MW_HRP = Soln_pptys_in(1,3)*1000; %pick a random one to use

beta_conc1 = beta_conc(1);
beta_conc2 = beta_conc(2);

parfor i = 1:length(k_offs)    
    %Locate Conc at which sensitivity is at
    Sens_Loc(i) = rcSso_NP_Saliva_LoD_variation_Solver(Sensitivity_Location,Sample_Conc_Guess,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,i);
    
    Sens_Loc_upper = Sens_Loc(i)+0.1e-9.*1e3; %Get upper, shifted by 0.1 nM
    Sens_Loc_lower = Sens_Loc(i)-0.1e-9.*1e3; %Get upper, shifted by 0.1 nM
    
    Sens_upper = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only(Sens_Loc_upper,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,i); 
    Sens_lower = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only(Sens_Loc_lower,Soln_pptys_in,Paper_params_in,op_params_in,M,N,InletBC,i); 
    
    Sens_upper_Cyan = beta_conc1.*(1-exp(-beta_conc2.*sum(Sens_upper)));
    Sens_lower_Cyan = beta_conc1.*(1-exp(-beta_conc2.*sum(Sens_lower)));
    
    Sens_All(i) = (Sens_upper_Cyan-Sens_lower_Cyan)/((Sens_Loc_upper-Sens_Loc_lower).*1e-3.*1e9);
end

figure;

%Sens_All = Sens_All; %Convert to /(ug/mL)
plot(k_offs,Sens_All,'LineWidth',1.5)
xlabel('k_{off} (s^{-1})'); ylabel('Sensitivity (AU\/nM)');
ax = gca;
ax.FontSize = 12;
axis tight;

cWD = pwd; %Current Working Folder

% Navigate to destination folder
cd(strcat(cWD,'\Figures'))
saveas(gcf,'rcSso_NP_Saliva_vary_k_off_Sensitivity_changes.fig')
save('Change_k_off_Sens_Saliva.mat','k_offs','Sens_All');

% Return
cd(cWD)

end