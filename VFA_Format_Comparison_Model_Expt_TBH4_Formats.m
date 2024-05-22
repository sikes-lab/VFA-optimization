%close all; clear; % Clean up
clear;
% This function is to compare the variation of Cyan Intensity with
% concentration of TB using different assay formats

% Load relevant files
load('Sorted_DiffFormat_Data.mat'); %data containing file
load('Black_HRP_Calibration_beta.mat'); %load parameters to convert pg on paper to Cyan Intensity
% Fitted graph is y = beta_conc(1)*(1-exp(-beta_conc(1)*m_HRP)) where m_HRP
% is in pg. y = Cyan Intensity @ 1 min

load('Black_flow_rates.mat'); %load flow rate data (key parameters: sample_flow_4L, wash_flow_4L)
load('opt_TBH4_params_Cellulose_NSB.mat')
load('Simple_fit_efficiency.mat'); %Fitted efficiencies for TopReadout

TBH4_NSB = opt_params_glob;

SA_NSB = opt_params_glob;

% Place 'Premix_VFA_parfor_nHRP_only.m',
% 'Sequential_VFA_parfor_nHRP_only.m',
% 'Top_Readout_Soln_Premix_parfor_nHRP_only.m','Top_Readout_Rxn.m'
% & 'Polson_D_est.m' into the same
% folder as this file.

%% Extract data of interest
% Note that concentrations are in 'pM' assuming 100 kDa. This converts to
% 100 ng/mL per pM. Here, we will convert this to the actual assumed
% concentration

% Extract for Format A
conc_hold = Format.A((Format.A(:,2) == 1),1).*1e-9.*1e3; % was in nM
intensity_hold = Format.A((Format.A(:,2) == 1),6);
std_hold = Format.A((Format.A(:,2) == 1),7);
Format_A_Values = [conc_hold intensity_hold std_hold];

% Extract for Format B
conc_hold = Format.B((Format.B(:,2) == 1),1).*1e-9.*1e3; % was in nM
intensity_hold = Format.B((Format.B(:,2) == 1),6);
std_hold = Format.B((Format.B(:,2) == 1),7);
Format_B_Values = [conc_hold intensity_hold std_hold];

% Extract for Format C
conc_hold = Format.C((Format.C(:,2) == 1),1).*1e-9.*1e3; % was in nM
intensity_hold = Format.C((Format.C(:,2) == 1),6);
std_hold = Format.C((Format.C(:,2) == 1),7);
Format_C_Values = [conc_hold intensity_hold std_hold];


%% Setup fixed parameters
conc_plot = linspace(0,100,51); %Conc range for probing, in nM
tested_concs = Format.A(1:6,1);

Soln_pptys = zeros(1,12);

Soln_pptys(1) = 35; %TB H4 MW
Soln_pptys(2) = 53.6; %BA-MBP-Sso Construct
Soln_pptys(3) = 44+44+53; %Molecular weight of SA-HRP (Assume 2 HRPs, NA is 60 kDa, HRP is 44 kDa)
Target_MW = (Soln_pptys(1) + Soln_pptys(2) + Soln_pptys(3))*1e3; %MW of complex (assumed all to have similar diffusivities)

Soln_pptys(4) = (1.1e5).*1e-3; %[m3/(mol.s)] rcSso.TBH4.E2
Soln_pptys(5) = (4.7e-4); %[1/s] rcSso.TBH4.E2

Soln_pptys(6) = (2.3e5).*1e-3; %[m3/(mol.s)] rcSso.TBH4.E1
Soln_pptys(7) = (7.6e-5); %[1/s] rcSso.TBH4.E1

Soln_pptys(8) = (50e-9)*1e3; %[mol/m3] Used 50 nM Rep in Sample solution
Soln_pptys(9) = (400*1e-12*1e3)*100/(Soln_pptys(3)); %[mol/m3] Used '400pM' SA-HRP in Sample solution
% Convert from 100 kDa to ~141 kDa

Soln_pptys(10) = opt_params_glob(2);
Soln_pptys(11) = opt_params_glob(3);
% Convert from ugmL to mol/m3

Paper_params = zeros(1,8);
Paper_params(1) = Polson_D_est(Target_MW); %[m2/s] Diffusivity of complex estimated via Polson
Paper_params(2) = 5.5e-6; %[m] Radius of pore, taken at 11um diameter
Paper_params(3) = 180e-6; %[m] Thickness of cellulose filter paper
Paper_params(7) = 0.689; %Void fraction of cellulose paper from Madhu et al
Paper_params(4) = 1.25e-3; %[m] Radius of whole well - 2.5mm for bottom layers, 3.2mm for top layer, but for my sanity let us ignore the top layer first...

Paper_params(5) = (0.5e-6)*(30e-6); %[mol] CBD loaded - 30uM, 0.5uL
Paper_params(6) = opt_params_glob(1);

op_params = cell(1,7);

op_params{1,1} = [sample_flow_4L wash_flow_4L]; % [m3/s] 
op_params{1,2} = [10 20].*(1e-6)*1e-3; %[m3] 10 uL sample + 20uL wash

op_params{1,4} = "Premix";
op_params{1,5} = [3600*3 15]; %Incubated premix for 3h, followed by 15s with sample

op_params{1,6} = 10; %Premix is 10-fold that used in sample
op_params{1,7} = [1 0]; %Sample, then wash
% M = 26; N = 90; %Choose N such that mesh size is 2 um

% Setup matrices for parfor

Paper_params_in = repmat(Paper_params,length(conc_plot),1);

Soln_pptys_in = repmat(Soln_pptys,length(conc_plot),1);
op_params_in = cell(length(conc_plot),3);
[op_params_in{:,1}] = deal(op_params{:,1});
[op_params_in{:,2}] = deal(op_params{:,2});

[op_params_in{:,4}] = deal(op_params{:,4});
[op_params_in{:,5}] = deal(op_params{:,5});
[op_params_in{:,6}] = deal(op_params{:,6});
[op_params_in{:,7}] = deal(op_params{:,7});

M = 26; N = 90; %Mesh sizes


%% Format A: Standard Premix (Sample & Reporters) and Application
% Intialize
Cyan_HRP_model_pos_A = zeros(size(conc_plot));
Cyan_HRP_model_neg_A = zeros(size(conc_plot));

MW_HRP = Soln_pptys(3)*1000;

parfor i = 1:length(conc_plot)
    Model_output = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only_3_Stoichio_by10(conc_plot(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);

    Cyan_HRP_model_pos_A(i) = sum(Model_output); %Note that it is the sum for the positive signal.
    Cyan_HRP_model_neg_A(i) = Model_output(2); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_neg_A = beta_conc(1).*(1-exp(-beta_conc(2).*Cyan_HRP_model_neg_A)); % Change to cya
Cyan_HRP_model_pos_A = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_pos_A)));% Change to cyan;

% Get Output for Format A to compare w/ expt
% Setup matrices for parfor

Paper_params_in = repmat(Paper_params,length(conc_plot),1);

Cyan_HRP_model_Expt_pos_A = zeros(size(tested_concs));

parfor i = 1:length(tested_concs)
    Model_output = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only(tested_concs(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);
    Cyan_HRP_model_Expt_pos_A(i) = sum(Model_output); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_Expt_pos_A = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_Expt_pos_A)));% Change to cyan;

%% Format B: Sequential VFA
% Intialize
Cyan_HRP_model_pos_B = zeros(size(conc_plot));
Cyan_HRP_model_neg_B = zeros(size(conc_plot));

% Set Appropriate op_params
op_params = cell(1,7);

op_params{1,1} = [sample_flow_4L sample_flow_4L wash_flow_4L]; % [m3/s] 
op_params{1,2} = [10 10 20].*(1e-6)*1e-3; %[m3] 10 uL sample + 10 uL label + 20uL wash

op_params{1,4} = "Premix";
op_params{1,5} = 3600*3; %Incubated premix for 3h, followed by 15s with sample

op_params{1,6} = 10; %Premix is 10-fold that used in sample
op_params{1,7} = [1 1 0]; %Sample, label then wash

% M = 26; N = 90; %Choose N such that mesh size is 2 um
op_params_in = cell(length(conc_plot),3);
[op_params_in{:,1}] = deal(op_params{:,1});
[op_params_in{:,2}] = deal(op_params{:,2});

[op_params_in{:,4}] = deal(op_params{:,4});
[op_params_in{:,5}] = deal(op_params{:,5});
[op_params_in{:,6}] = deal(op_params{:,6});
[op_params_in{:,7}] = deal(op_params{:,7});

parfor i = 1:length(conc_plot)
    Model_output = MW_HRP.*1e12.*Sequential_VFA_parfor_nHRP_only(conc_plot(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);

    Cyan_HRP_model_pos_B(i) = sum(Model_output); %Note that it is the sum for the positive signal.
    Cyan_HRP_model_neg_B(i) = Model_output(2); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_neg_B = beta_conc(1).*(1-exp(-beta_conc(2).*Cyan_HRP_model_neg_B)); % Change to cya
Cyan_HRP_model_pos_B = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_pos_B)));% Change to cyan;

Cyan_HRP_model_pos_B_Corrected = Cyan_HRP_model_pos_B - 0.04;

% Get Output for Format B to compare w/ expt

Cyan_HRP_model_Expt_pos_B = zeros(size(tested_concs));

parfor i = 1:length(tested_concs)
    Model_output = MW_HRP.*1e12.*Sequential_VFA_parfor_nHRP_only(conc_plot(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);

    Cyan_HRP_model_Expt_pos_B(i) = sum(Model_output); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_Expt_pos_B = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_Expt_pos_B)));% Change to cyan;
Cyan_HRP_model_Expt_pos_B_Corrected = Cyan_HRP_model_Expt_pos_B-0.04;


%% Format C: Premix all
% Intialize
Cyan_HRP_model_pos_C = zeros(size(conc_plot));
Cyan_HRP_model_neg_C = zeros(size(conc_plot));

Soln_pptys_in(:,10:11) = repmat(opt_params_glob(2:3),length(conc_plot),1);

Soln_pptys_in(:,12) = 50e-9*1e3; %50 nM CBD

% Set Appropriate op_params
op_params = cell(1,8);

op_params{1,1} = [sample_flow_4L wash_flow_4L]; % [m3/s] 
op_params{1,2} = [10 20].*(1e-6)*1e-3; %[m3] Shouldn't matter

op_params{1,4} = "Premix";
op_params{1,5} = [3600*3 20]; %Incubated premix for 3h, followed by 15s with sample

op_params{1,6} = 10; %Premix is 10-fold that used in sample
op_params{1,7} = [1 0]; %Shouldn't matter

op_params{1,8} = Avg_Capture_efficiency; %Capture efficiency

% M = 26; N = 90; %Choose N such that mesh size is 2 um
op_params_in = cell(length(conc_plot),3);
[op_params_in{:,1}] = deal(op_params{:,1});
[op_params_in{:,2}] = deal(op_params{:,2});

[op_params_in{:,4}] = deal(op_params{:,4});
[op_params_in{:,5}] = deal(op_params{:,5});
[op_params_in{:,6}] = deal(op_params{:,6});
[op_params_in{:,7}] = deal(op_params{:,7});
[op_params_in{:,8}] = deal(op_params{:,8});

parfor i = 1:length(conc_plot)
    Model_output = MW_HRP.*1e12.*Top_Readout_Soln_Premix_parfor_nHRP_only(conc_plot(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);

    Cyan_HRP_model_pos_C(i) = sum(Model_output); %Note that it is the sum for the positive signal.
    Cyan_HRP_model_neg_C(i) = Model_output(2); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_neg_C = beta_conc(1).*(1-exp(-beta_conc(2).*Cyan_HRP_model_neg_C)); % Change to cya
Cyan_HRP_model_pos_C = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_pos_C)));% Change to cyan;

Cyan_HRP_model_pos_C_Corrected = Cyan_HRP_model_pos_C + 0.065;

% Get Output for Format B to compare w/ expt

Cyan_HRP_model_Expt_pos_C = zeros(size(tested_concs));
Cyan_HRP_model_Expt_neg_C = zeros(size(tested_concs));

parfor i = 1:length(tested_concs)
    Model_output = MW_HRP.*1e12.*Top_Readout_Soln_Premix_parfor_nHRP_only(conc_plot(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);

    Cyan_HRP_model_Expt_pos_C(i) = sum(Model_output); %Note that it is the sum for the positive signal.
    Cyan_HRP_model_Expt_neg_C(i) = Model_output(2); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_Expt_pos_C = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_Expt_pos_C)));% Change to cyan;
Cyan_HRP_model_Expt_pos_C_Corrected = Cyan_HRP_model_Expt_pos_C + 0.065;

%% PLOT
%Now, convert conc to ng/mL

%Plot positive
%conc_plot = conc_plot.*1e9/1000; %Convert to nM
%Comparison plot
figure; 

% Format A
subplot(1,3,1)
plot(conc_plot, Cyan_HRP_model_pos_A,'LineWidth',1.5); 
hold on;
% plot(conc_plot, Cyan_HRP_model_neg_A,'LineWidth',1.5);
errorbar(Format.A(1:6,1),Format.A(1:6,6),Format.A(1:6,7),Format.A(1:6,7),'.','MarkerSize',17,'LineWidth',1.5,'Color','b');
axis tight; box on;
xlabel('[TB H4 Antigen] (nM)'); ylabel('Cyan Intensity(AU)');
xlim([0 100]); ylim([0.08 0.33]);
title('Premixing Reporter and label with flow-based capture');

% Get Deviations

Deviations_FormatA_abs = abs(Cyan_HRP_model_Expt_pos_A-Format.A(1:6,6));
Deviations_FormatA_percent = 100*abs(Cyan_HRP_model_Expt_pos_A-Format.A(1:6,6))./Format.A(1:6,6);

Deviations_FormatA = [Deviations_FormatA_abs Deviations_FormatA_percent];

% Format B
subplot(1,3,2)
plot(conc_plot, Cyan_HRP_model_pos_B,'LineWidth',1.5); 
hold on;
plot(conc_plot, Cyan_HRP_model_neg_B,'LineWidth',1.5); 

plot(conc_plot, Cyan_HRP_model_pos_B_Corrected,'--','LineWidth',1.5); 

errorbar(Format.B(1:6,1),Format.B(1:6,6),Format.B(1:6,7),Format.B(1:6,7),'.','MarkerSize',17,'LineWidth',1.5,'Color','m');
axis tight; box on;
xlabel('[TB H4 Antigen] (nM)'); ylabel('Cyan Intensity(AU)');
xlim([0 100]); ylim([0.08 0.33]);
title('Sequential Sample and Reporter-label Application');

% Get Deviations

Deviations_FormatB_abs = abs(Cyan_HRP_model_Expt_pos_B-Format.B(1:6,6));
Deviations_FormatB_percent = 100*abs(Cyan_HRP_model_Expt_pos_B-Format.B(1:6,6))./Format.B(1:6,6);

Deviations_FormatB_Adj_abs = abs(Cyan_HRP_model_Expt_pos_B_Corrected-Format.B(1:6,6));
Deviations_FormatB_Adj_percent = 100*abs(Cyan_HRP_model_Expt_pos_B_Corrected-Format.B(1:6,6))./Format.B(1:6,6);

Deviations_FormatB = [Deviations_FormatB_abs Deviations_FormatB_percent Deviations_FormatB_Adj_abs Deviations_FormatB_Adj_percent];

% Format C
subplot(1,3,3)
plot(conc_plot(Cyan_HRP_model_pos_C~=0), Cyan_HRP_model_pos_C(Cyan_HRP_model_pos_C~=0),'LineWidth',1.5); 
hold on;
plot(conc_plot, Cyan_HRP_model_neg_C,'LineWidth',1.5); 

plot(conc_plot(Cyan_HRP_model_pos_C_Corrected~=0.05), Cyan_HRP_model_pos_C_Corrected(Cyan_HRP_model_pos_C_Corrected~=0.05),'--','LineWidth',1.5); 

errorbar(Format.C(1:6,1),Format.C(1:6,6),Format.C(1:6,7),Format.C(1:6,7),'.','MarkerSize',17,'LineWidth',1.5,'Color','r');
axis tight; box on;
xlabel('[TB H4 Antigen] (nM)'); ylabel('Cyan Intensity(AU)');
xlim([0 100]); ylim([0.08 0.33]);
title('Full Solution-based mixing');

% Get Deviations

Deviations_FormatC_abs = abs(Cyan_HRP_model_Expt_pos_C-Format.C(1:6,6));
Deviations_FormatC_percent = 100*abs(Cyan_HRP_model_Expt_pos_C-Format.C(1:6,6))./Format.C(1:6,6);

Deviations_FormatC_Adj_abs = abs(Cyan_HRP_model_Expt_pos_C_Corrected-Format.C(1:6,6));
Deviations_FormatC_Adj_percent = 100*abs(Cyan_HRP_model_Expt_pos_C_Corrected-Format.C(1:6,6))./Format.C(1:6,6);

Deviations_FormatC = [Deviations_FormatC_abs Deviations_FormatC_percent Deviations_FormatC_Adj_abs Deviations_FormatC_Adj_percent];

save('rcSso_TBH4_Format_Comparison_Deviations.mat','Deviations_FormatA','Deviations_FormatB','Deviations_FormatC');

%% Save


cWD = pwd; %Current Working Folder

% Navigate to destination folder
cd(strcat(cWD,'\Figures'))
saveas(gcf,'rcSso_TBH4_Format_ABC_vary_Conc_plot_changeNSB.fig')

% Return
cd(cWD)