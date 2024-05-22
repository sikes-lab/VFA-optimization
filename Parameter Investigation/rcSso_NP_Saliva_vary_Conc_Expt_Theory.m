%close all; clear; % Clean up
clear;
% This function is to compare the variation of Cyan Intensity with
% concentration of NP measured by experiment and predicted by model for 
% BA-MBP-rcSso.NP-PF2.1 and rcSso.BacNP2.1-CBD

% Load relevant files
load('Sorted_NP_Saliva_Data.mat'); %data containing file
load('White_Saliva_HRP_Calibration_beta_Cleaned.mat'); %load parameters to convert pg on paper to Cyan Intensity
% Fitted graph is y = beta_conc(1)*(1-exp(-beta_conc(1)*m_HRP)) where m_HRP
% is in pg. y = Cyan Intensity @ 1 min

load('White_Saliva_Flow_Rates.mat'); %load flow rate data (key parameters: sample_flow_rate, avg_wash_flow_rate)
load('opt_NP_Saliva_params_all_1layer_Cellulose_NSB.mat')

% Place 'Premix_VFA_parfor_nHRP_only.m' & 'Polson_D_est.m' into the same
% folder as this file

%% Extract data of interest
% Note that concentrations are in 'pM' assuming 100 kDa. This converts to
% 100 ng/mL per pM. Here, we will convert this to the actual assumed
% concentration

% Calibration data

% Extract for positive
conc_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Calib((Sorted_NP_Saliva_Data.Antigen_Conc_Calib(:,2) == 1),1).*1e-9.*1e3; % was in nM
intensity_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Calib((Sorted_NP_Saliva_Data.Antigen_Conc_Calib(:,2) == 1),7);
std_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Calib((Sorted_NP_Saliva_Data.Antigen_Conc_Calib(:,2) == 1),8);
Cyan_Data_Calib = [conc_hold intensity_hold std_hold];

% Test data

% Extract for positive
conc_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Test((Sorted_NP_Saliva_Data.Antigen_Conc_Test(:,2) == 1),1).*1e-9.*1e3; % was in nM
intensity_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Test((Sorted_NP_Saliva_Data.Antigen_Conc_Test(:,2) == 1),7);
std_hold = Sorted_NP_Saliva_Data.Antigen_Conc_Test((Sorted_NP_Saliva_Data.Antigen_Conc_Test(:,2) == 1),8);
Cyan_Data_Expt = [conc_hold intensity_hold std_hold];

%% Setup fixed parameters
conc_plot = linspace(0,40,41); %Conc range for probing

Soln_pptys = zeros(1,9);

Soln_pptys(1) = 46; %NP MW
Soln_pptys(2) = 53.6; %BA-MBP-Sso Construct
Soln_pptys(3) = 44+44+60; %Molecular weight of NA-HRP (Assume 2 HRPs, NA is 60 kDa, HRP is 44 kDa)
Target_MW = (Soln_pptys(1) + Soln_pptys(2) + Soln_pptys(3))*1e3; %MW of complex (assumed all to have similar diffusivities)

Soln_pptys(4) = (2.14e5).*1e-3; %[m3/(mol.s)] rcSso.BacNP2.1 (E2)
Soln_pptys(5) = (8.47e-4); %[1/s] rcSso.BacNP2.1 (E2) Use the one with His
%Tag (Refer to 20200616_NTU-SMART-MIT meeting)

Soln_pptys(6) = (2.03e5).*1e-3; %[m3/(mol.s)] rcSso.NP-PF2.1 (E1)
Soln_pptys(7) = (6.46e-4); %[1/s] rcSso.NP-PF2.1 (E1) Use the one with His Tag

% Soln_pptys(6) = (1.54e5).*1e-3; %[m3/(mol.s)] rcSso.NP-PF2.1 (E1)
% Soln_pptys(7) = (9.14e-4); %[1/s] rcSso.NP-PF2.1 (E1) Use the one with His Tag

Soln_pptys(8) = (50e-9)*1e3; %[mol/m3] Used 50 nM Rep in Sample solution
Soln_pptys(9) = (250e-12.*1e3*((100)))./(Soln_pptys(3)); %[mol/m3] Used 250 pM SA-HRP final concentration
% Convert from ngmL to mol/m3

Soln_pptys(10) = opt_params_glob(2);
Soln_pptys(11) = opt_params_glob(3);
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
op_params{1,2} = [20 50].*(1e-6)*1e-3; %[m3] 15 uL sample + 50 (25+25) uL wash

op_params{1,4} = "Premix";
op_params{1,5} = [1*60 25.5]; %Incubated premix for 15 min, followed by 15s with sample

op_params{1,6} = 2; %Premix is 15-fold that used in sample
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

% for i = 1:length(conc_plot)
%     op_params_in{i,3} = [conc_plot(i) 0].*1e-9.*1e3; %first convert to [nM] (1e-9), then to [mol/m3] (1e3) inlet concentration of target species    
% end

M = 26; N = 90; %Mesh sizes

% Intialize
Cyan_HRP_model_pos = zeros(size(conc_plot));
Cyan_HRP_model_neg = zeros(size(conc_plot));

MW_HRP = Soln_pptys(3)*1000;

parfor i = 1:length(conc_plot)
    Model_output = MW_HRP.*1e12.*Premix_VFA_parfor_nHRP_only(conc_plot(i).*1e-9.*1e3,Soln_pptys_in,Paper_params_in,op_params_in,M,N,'Danckwertz',i);

    Cyan_HRP_model_pos(i) = sum(Model_output); %Note that it is the sum for the positive signal.
    Cyan_HRP_model_neg(i) = Model_output(2); %Note that it is the sum for the positive signal.
end

Cyan_HRP_model_neg = beta_conc(1).*(1-exp(-beta_conc(2).*Cyan_HRP_model_neg)); % Change to cya
Cyan_HRP_model_pos = beta_conc(1).*(1-exp(-beta_conc(2).*(Cyan_HRP_model_neg+Cyan_HRP_model_pos)));% Change to cyan;

% Cyan_HRP_model_neg = beta_all(1).*(1-exp(-beta_all(2).*Cyan_HRP_model_neg)); % Change to cya
% Cyan_HRP_model_pos = beta_all(1).*(1-exp(-beta_all(2).*(Cyan_HRP_model_neg+Cyan_HRP_model_pos)));% Change to cyan;

%% PLOTTING TIME
%Now, convert conc to ng/mL

%Plot positive
%conc_plot = conc_plot.*1e9/1000; %Convert to nM
%Comparison plot
figure;
semilogx(conc_plot, Cyan_HRP_model_pos,'LineWidth',1.5); 
hold on;
% Plot experimental values NOT used for fitting
errorbar(Cyan_Data_Expt(:,1).*1e9./1000,Cyan_Data_Expt(:,2),Cyan_Data_Expt(:,3),'.','MarkerSize',17,'LineWidth',1.5); %use sample flow - also plot horizontal error bars (need to state both yneg and ypos, xneg, xpos
% Plot values used for fitting
errorbar(Cyan_Data_Calib(:,1).*1e9./1000,Cyan_Data_Calib(:,2),Cyan_Data_Calib(:,3),'.','MarkerSize',17,'LineWidth',1.5); %use sample flow - also plot horizontal error bars (need to state both yneg and ypos, xneg, xpos

semilogx(conc_plot, Cyan_HRP_model_neg,'--','LineWidth',1.5); 
semilogx(conc_plot, Cyan_HRP_model_neg(1).*ones(size(conc_plot)),'LineWidth',1.5); 

title('Influence of concentration on Signal generated'); box on;
xlabel('Concentration (ng/mL)'); ylabel('Cyan Intensity (AU)');

lgd = legend('Theoretical positive signal','Experimental positive signal','Experimental positive signal (for fitting)','Theoretical background signal');
ax = gca;
ax.FontSize = 12;
lgd.FontSize = 12;
axis tight;

cWD = pwd; %Current Working Folder

% Navigate to destination folder
cd(strcat(cWD,'\Figures'))
saveas(gcf,'rcSso_NP_Saliva_vary_Conc_plot.fig')

% Return
cd(cWD)