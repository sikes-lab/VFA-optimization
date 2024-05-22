clear; close all; % Clean up
% Load relevant files
load('Sorted_NP_Saliva_Data.mat'); %data containing file
load('White_Saliva_HRP_Calibration_beta_Cleaned.mat'); %load parameters to convert pg on paper to Cyan Intensity
% Fitted graph is y = beta_all(1)*(1-exp(-beta_all(1)*m_HRP)) where m_HRP
% is in pg. y = Cyan Intensity @ 1 min

load('White_Saliva_Flow_Rates.mat'); %load flow rate data (key parameters: sample_flow_rate, avg_wash_flow_rate)

% Place 'Premix_VFA_parfor_nHRP_only.m' & 'Polson_D_est.m' into the same
% folder as this file

%% Extract concentrations for calibration

% Extract for positive
conc_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Calib((Sorted_NP_Saliva_Data.Antigen_Conc_Calib(:,2) == 1),1).*1e-9.*1e3; % was in nM
intensity_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Calib((Sorted_NP_Saliva_Data.Antigen_Conc_Calib(:,2) == 1),3:6);

avg1min_Cyan_pos_Cyan = [repmat(conc_hold,4,1) intensity_hold(:)]; %This flattens intensity so that it fits across all 3 repetitions
pg_hold = log(1-(intensity_hold./beta_all(1)))./(-beta_all(2));
avg1min_Cyan_pos_pg = [repmat(conc_hold,4,1) pg_hold(:)]; %This flattens intensity so that it fits across all 4 repetitions

% Extract for negative
neg_intensity_hold = Sorted_NP_Saliva_Data.HRP_Conc_NSB_Calib((Sorted_NP_Saliva_Data.HRP_Conc_NSB_Calib(:,1) == 25),:);
neg_intensity_hold = mean(neg_intensity_hold((neg_intensity_hold(:,2) == 1),7));
neg_pg = log(1-(neg_intensity_hold./beta_all(1)))./(-beta_all(2));

avg1min_Cyan_pos_pg(:,2) = avg1min_Cyan_pos_pg(:,2)-neg_pg; %get net positive signals

% Extract full set of negative
conc_hold = Sorted_NP_Saliva_Data.HRP_Conc_NSB_Calib((Sorted_NP_Saliva_Data.HRP_Conc_NSB_Calib(:,2) == 1),1); %in ug/mL
conc_hold = conc_hold.*1e-9.*1e6/((44+44+53)*1e3); %Convert concentration to mol/m3
intensity_hold = Sorted_NP_Saliva_Data.HRP_Conc_NSB_Calib((Sorted_NP_Saliva_Data.HRP_Conc_NSB_Calib(:,2) == 1),3:6);

avg1min_Cyan_neg_Cyan = [repmat(conc_hold,4,1) intensity_hold(:)]; %This flattens intensity so that it fits across all 3 repetitions
pg_hold = log(1-(intensity_hold./beta_all(1)))./(-beta_all(2));
avg1min_Cyan_neg_pg = [repmat(conc_hold,4,1) pg_hold(:)]; %This flattens intensity so that it fits across all 3 repetitions

%% Setup fixed parameters

% MW = (53+44+44)*1e3; % Molecular weight of SA-HRP (Assume 2 HRPs, SA is 53 kDa, HRP is 44 kDa)
% Above accounted for in code
Soln_pptys = zeros(1,11);

Soln_pptys(1) = 46; %NP MW
Soln_pptys(2) = 53.6; %BA-MBP-Sso Construct
Soln_pptys(3) = 44+44+60; %Molecular weight of NA-HRP (Assume 2 HRPs, NA is 60 kDa, HRP is 44 kDa)
Target_MW = (Soln_pptys(1) + Soln_pptys(2) + Soln_pptys(3))*1e3; %MW of complex (assumed all to have similar diffusivities)

%Soln_pptys(4) = (3.94e4).*1e-3; %[m3/(mol.s)] rcSso.BacNP2.1 (E2)
%Soln_pptys(5) = (1.73e-4); %[1/s] rcSso.BacNP2.1 (E2)
Soln_pptys(4) = (2.14e5).*1e-3; %[m3/(mol.s)] rcSso.BacNP2.1 (E2)
Soln_pptys(5) = (8.47e-4); %[1/s] rcSso.BacNP2.1 (E2) Use the one with His
%Tag (Refer to 20200616_NTU-SMART-MIT meeting)

%Soln_pptys(6) = (3.59e4).*1e-3; %[m3/(mol.s)] rcSso.NP-PF2.1 (E1)
%Soln_pptys(7) = (1.16e-4); %[1/s] rcSso.NP-PF2.1 (E1) 

Soln_pptys(6) = (2.03e5).*1e-3; %[m3/(mol.s)] rcSso.NP-PF2.1 (E1)
Soln_pptys(7) = (6.46e-4); %[1/s] rcSso.NP-PF2.1 (E1) Use the one with His Tag


Soln_pptys(8) = (50e-9)*1e3; %[mol/m3] Used 50 nM Rep in Sample solution
Soln_pptys(9) = (250e-12.*1e3*((100)))./(Soln_pptys(3)); %[mol/m3] Used 250 pM SA-HRP final concentration
% Convert from ngmL to mol/m3

Paper_params_NSB = zeros(1,7);
Paper_params_NSB(1) = Polson_D_est(Target_MW); %[m2/s] Diffusivity of complex estimated via Polson
Paper_params_NSB(2) = 5.5e-6; %[m] Radius of pore, taken at 11um diameter
Paper_params_NSB(3) = 180e-6; %[m] Thickness of cellulose filter paper
Paper_params_NSB(7) = 0.689; %Void fraction of cellulose paper from Madhu et al
Paper_params_NSB(4) = 1.25e-3; %[m] Radius of whole well - 2.5mm for bottom layers, 3.2mm for top layer, but for my sanity let us ignore the top layer first...

Paper_params = Paper_params_NSB; %Copy over parameters
Paper_params(5) = (0.5e-6)*(30e-6); %[mol] CBD loaded - 30uM, 0.5uL

op_params = cell(1,7);

op_params{1,1} = [sample_flow_rate avg_wash_flow_rate].*1e-6.*1e-3; % [m3/s] 
op_params{1,2} = [20 50].*(1e-6)*1e-3; %[m3] 20 uL sample + 50 (25+25) uL wash

op_params{1,4} = "Premix";
op_params{1,5} = [5*60 20]; %Incubated premix for 5 min, followed by 15s with sample

op_params{1,6} = 2; %Premix is 2-fold that used in sample
op_params{1,7} = [1 0]; %Sample, then wash
% M = 26; N = 90; %Choose N such that mesh size is 2 um

%% Optimize for Accessible CBD
Eff_CBDfrac_guess = 0.5; %Initial guess
A = []; b = []; 
Aeq =[]; beq=[];
lb = 0; ub = 1; %set bounds for accesible CBD frac
options = optimset('TolFun',1e-8,'MaxIter',1e6,'MaxFunEvals',1e6,'Display','iter');
rng default
t0 = tic();
opt_Eff_CBDfrac = fmincon(@(E_CBDfrac) min_Calibrate_CBDfrac(E_CBDfrac,Soln_pptys,Paper_params,op_params,avg1min_Cyan_pos_pg),Eff_CBDfrac_guess,A,b,Aeq,beq,lb,ub,[],options);
toc(t0)

%% Optimize for NSB kinetic constants

%Fill in obtained effective CBD into Paper_params
Paper_params(6) = opt_Eff_CBDfrac; % Freeze opt_Eff_CBDfrac

%Return to actual value
avg1min_Cyan_pos_pg(:,2) = avg1min_Cyan_pos_pg(:,2)+neg_pg;

% Paper_params.Eff_CBDfrac = 0.0210;
K_guess = [2e-3 0.5]; %Initial guess

Param_guess = K_guess;

rng default %reset!
A = []; b = []; 
Aeq =[]; beq=[];
lb = [0 0]; ub = [100 10]; %Focus only on optimization of k_on and k_D NSB

options = optimset('TolFun',1e-8,'MaxIter',1e6,'MaxFunEvals',1e6,'Display','iter');
%opt_params = fmincon(@(E_params) min_Calibrate_NSB(E_params,rate_params,Paper_params,op_params,pg_Data),Param_guess,A,b,Aeq,beq,lb,ub,@(E_params) HRP_constraint(E_params,rate_params,Paper_params,op_params,data),options);
opt_params = fmincon(@(E_params) min_Calibrate_NSB(E_params,Soln_pptys,Paper_params,op_params,avg1min_Cyan_pos_pg,avg1min_Cyan_neg_pg),Param_guess,A,b,Aeq,beq,lb,ub,[],options);
opt_ALL_k_on = opt_params(1);
opt_ALL_k_off = opt_params(2);

%Try Global search
rng default
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',opt_params,... % Use optimized parameters for global search
    'objective',@(E_params) min_Calibrate_NSB(E_params,Soln_pptys,Paper_params,op_params,avg1min_Cyan_pos_pg,avg1min_Cyan_neg_pg),'lb',lb,'ub',ub,'options',options);
glob_min_k_on_K_D = run(gs,problem);

opt_params = [opt_Eff_CBDfrac opt_params]; %just put it in the format i want

opt_params_glob = [opt_Eff_CBDfrac glob_min_k_on_K_D];

save('opt_NP_Saliva_params_NSB.mat','opt_params','opt_params_glob') %save optimized parameters

%% Minimization Objective Functions

%% Fitting Effective CBD fraction
function F = min_Calibrate_CBDfrac(E_CBDfrac,Soln_pptys,Paper_params_SB,op_params,pg_Data)
MW = (Soln_pptys(3))*1e3; % Molecular weight of SA-HRP (Assume 2 HRPs, SA is 53 kDa, HRP is 44 kDa)
Paper_params_SB(6) = E_CBDfrac; %Assign effective CBD fraction

conc = pg_Data(:,1); 
Pure_Adsorbed_Sso = pg_Data(:,2);

Soln_pptys(10) = 0; % placeholders
Soln_pptys(11) = 1; % placeholders

unique_conc = unique(conc); %pull out unique instances to minimize computation time

M = 26; N = 90;
HRP_Adsorbed_model = zeros(size(Pure_Adsorbed_Sso));
HRP_Adsorbed_unique = zeros(size(unique_conc));

% Create rate_params, Paper_params array
Soln_pptys_in = repmat(Soln_pptys,length(unique_conc),1);
Paper_params_in = repmat(Paper_params_SB,length(unique_conc),1);

% Create op_params array
op_params_in = cell(length(unique_conc),7); 
[op_params_in{:,1}] = deal(op_params{1,1});
[op_params_in{:,2}] = deal(op_params{1,2});

[op_params_in{:,4}] = deal(op_params{1,4});
[op_params_in{:,5}] = deal(op_params{1,5});
[op_params_in{:,6}] = deal(op_params{1,6});
[op_params_in{:,7}] = deal(op_params{1,7});

for i = (1:length(unique_conc))
    HRP_Adsorbed = Premix_VFA_parfor_nHRP_only(unique_conc(i),Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);
    HRP_Adsorbed_unique(i) = HRP_Adsorbed(1); %Take specific only
end

for i = 1:length(conc)
    for j = 1:length(unique_conc) %look through and assign (save computational power!)
        if conc(i) == unique_conc(j)
            HRP_Adsorbed_model(i) = HRP_Adsorbed_unique(j);
        end
    end
end

% Convert to pg
HRP_Adsorbed_model = HRP_Adsorbed_model.*1e12.*(MW);
% Using MAE
F = sum((abs(HRP_Adsorbed_model-Pure_Adsorbed_Sso)));

% Using MSE
%F = sum((((HRP_Adsorbed_model-Pure_Adsorbed_Sso))).^2); %compare pg
end

%% Minimization functions
function F = min_Calibrate_NSB(E_params,Soln_pptys,Paper_params,op_params,pg_Data_pos,pg_Data_neg)
% Paper_params.Eff_CBDfrac = E_params(1);
MW = (Soln_pptys(3))*1e3; % Molecular weight of SA-HRP
M = 26; N = 90;
pos_conc = pg_Data_pos(:,1); %Get conc for HRP
neg_conc = pg_Data_neg(:,1); %Get conc for HRP
unique_pos_conc = unique(pos_conc);
unique_neg_conc = unique(neg_conc);

HRP_Adsorbed_model_neg = zeros(size(neg_conc));
HRP_Adsorbed_model_pos = zeros(size(pos_conc));

Soln_pptys(10) = E_params(1); %guessed k_on
Soln_pptys(11) = E_params(2); %guessed k_off

HRP_Adsorbed_neg_unique = zeros(size(unique_neg_conc));
HRP_Adsorbed_pos_unique = zeros(size(unique_pos_conc));

% Create rate_params, Paper_params array
% For positive
Soln_pptys_in = repmat(Soln_pptys,length(unique_pos_conc),1);
Paper_params_in = repmat(Paper_params,length(unique_pos_conc),1);

% Create op_params array
op_params_in = cell(length(unique_pos_conc),7); 
[op_params_in{:,1}] = deal(op_params{1,1});
[op_params_in{:,2}] = deal(op_params{1,2});

[op_params_in{:,4}] = deal(op_params{1,4});
[op_params_in{:,5}] = deal(op_params{1,5});
[op_params_in{:,6}] = deal(op_params{1,6});
[op_params_in{:,7}] = deal(op_params{1,7});

parfor i = (1:length(unique_pos_conc))
    HRP_Adsorbed = Premix_VFA_parfor_nHRP_only(unique_pos_conc(i),Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);
    HRP_Adsorbed_pos_unique(i) = sum(HRP_Adsorbed); %Take both    
end

for i = 1:length(pos_conc)
    for j = 1:length(unique_pos_conc) %look through and assign (save computational power!)
        if pos_conc(i) == unique_pos_conc(j)
            HRP_Adsorbed_model_pos(i) = HRP_Adsorbed_pos_unique(j);
        end
    end
end

% For negative
Soln_pptys_in = repmat(Soln_pptys,length(unique_neg_conc),1);
Soln_pptys_in(:,9) = unique_neg_conc./op_params{1,6}; %Should divide by Dilution factor for Swab!

parfor i = (1:length(unique_neg_conc))
    HRP_Adsorbed = Premix_VFA_parfor_nHRP_only(0,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);
    HRP_Adsorbed_neg_unique(i) = HRP_Adsorbed(2); %Take non-specific only
end

for i = 1:length(neg_conc)
    for j = 1:length(unique_neg_conc) %look through and assign (save computational power!)
        if neg_conc(i) == unique_neg_conc(j)
            HRP_Adsorbed_model_neg(i) = HRP_Adsorbed_neg_unique(j);
        end
    end
end

% Convert to pg
HRP_Adsorbed_model_neg = HRP_Adsorbed_model_neg.*1e12.*(MW); %Converting to pg adsorbed instead of moles
HRP_Adsorbed_model_pos = HRP_Adsorbed_model_pos.*1e12.*(MW); %Converting to pg adsorbed instead of moles

% Using MSE
F = sum((((pg_Data_neg(:,2)-HRP_Adsorbed_model_neg))).^2) + sum((((pg_Data_pos(:,2)-HRP_Adsorbed_model_pos))).^2); %sum up both errors

end

function D = Polson_D_est(MW)

% This function applies Polson's correlation to estimate the Protein
% diffusion coefficient in DILUTE solutions @ 20 C - this is valid for the
% assay which is in the ~uM to nM range

% Input: 
%        MW = Molecular weight of protein in g/mol

% Output: 
%        D = Diffusivity of protein in m^2/s


A = 2.85e-5; %[cm^2 s^(-1) g^(1/3) mol^(-1/3)]
D = A/(MW^(1/3)); %Apply formula, here D = cm^2/s
D = D*1e-4; % Convert to m^2/s

end
