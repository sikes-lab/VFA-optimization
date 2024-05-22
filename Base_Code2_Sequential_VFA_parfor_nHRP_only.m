function n_HRP = Sequential_VFA_parfor_nHRP_only(c_analyte,Soln_pptys,Paper_params,op_params,M,N,InletBC,param_pt)
% This function codes for the complete assay from premixing in a tube,
% followed by loading onto cellulose paper. The first part, premixing, will
% be modeled via a Batch mixing model. The second part, loading, will apply
% the VFA model, based on a cylindrical pore and surface reactions
%
% NEW in this edition: NSB is calculated for all HRP molecules, alongside
% the specific ones in the during the specific binding step
%
% Inputs:
%  NOTE: ALL VALUES ARE TO BE IN SI UNITS
%
%        c_analyte = Analyte (NP) concentration in solution (input in [mol/m3])
%        Soln_pptys = matrix containing values for analyte and binders (In columns) 
%                       (1) MW_analyte (kDa) 
%                       (2) MW_rep (kDa) - the BA-MBP reporter 
%                       (3) MW_SAHRP 
%                       (4) k_on_capture [m3/(mol.s)]
%                       (5) k_off_capture [1/s]
%                       (6) k_on_rep [m3/(mol.s)]
%                       (7) k_off_rep [1/s]
%                       (8) c_rep [mol/m3]
%                       (9) c_HRP [mol/m3]
%                       (10) k_on_capture_NSB [m3/(mol.s)]
%                       (11) k_off_capture_NSB [1/s]
%                       for each set of parameters (rows)
%
%        Paper_params = matrix containing paper parameters (In columns)
%                       (1) D, diffusivity (SI units)
%                       (2) R, pore radius
%                       (3) L, Thickness of paper
%                       (4) R_well
%                       (5) CBD, CBD loading
%                       (6) eta, effective CBD accessible
%                       (7) epsilon, porosity of paper
%                       for each set of parameters (rows)
%
%        op_params = cell containing VFA operating parameters (in columns)
%                       (1) Q, Flow rate [m3/s]
%                       (2) Vol, Volume of Sample passed [m3]. Sequential - expect 3 steps
%                       (3) sample_conc, sample concentration applied on
%                           each step [mol/m3]
%                       (4) MixingFormat 
%                           "All mix" or "Premix", indicative of whether all
%                           are added together and mixed, or premix, whereby SA-HRP and
%                           BA-MBP-PF2.1 are added together first to give Solution A, followed
%                           by application.
%                       (5) Sample Mixing Times [s] (1x1 vector if using
%                           "Premix"
%                           for each set of parameters (rows)
%                           Note that they are keyed in as a vector
%                           corresponding to the number of sample & wash
%                           application steps.
%                       (6) Dilution of Premix if using
%                       (7) SampleORWash: Check if loading sample or
%                           washing. Sequential - expect 3 steps
%
%        param_pt = set number in question (corresponds to row number in
%                   all parameter inputs)
% Output: 
%        n_HRP = Number of moles of HRP captured on paper

%% Checking important inputs
MixingFormat = op_params{param_pt,4};

% For Steps 0 and 1
if ~strcmp(MixingFormat,'All mix') && ~strcmp(MixingFormat,'Premix') 
    error("Please enter a valid sample mixing format: 'All mix' or 'Premix'.");
end

% For Step 2
No_FlowSteps = length(op_params{param_pt,2}); %counts number of flow steps.

if (length(op_params{param_pt,1}) ~= No_FlowSteps) && (length(op_params{param_pt,7}) ~= No_FlowSteps)
    error("Please ensure that the length of input op_params{?,1} and op_params{?,1} = input op_params{?,2}");
end

%% 0) Premix BA-MBP-Rep and SA-HRP (2nd solution)
% Use only if requested. 
options = odeset('AbsTol',1e-16,'RelTol',1e-10,'NonNegative',true); %Accuracy of ODE solver

if strcmp(MixingFormat,'Premix') 
    if length(op_params{param_pt,5}) ~= 1
        error("Please input a vector for op_params{?,5} of 1 time corresponding to the times used for premixing and sample mixing");
    end
    c0 = [Soln_pptys(param_pt,9) 0 0 0 0 Soln_pptys(param_pt,8)]; % Premixing done @ final concentration
    tspan = [0 op_params{param_pt,5}(1)];
    [~,c] = ode15s(@SA_BA_interaction,tspan,c0,options); %Binding parameters already in code
    c_premix = c(end,:); %to transfer to label application step
end

%% NO LONGER NECESSARY: Mix Solution A with Analyte (NP)

% if strcmp(MixingFormat,'Premix')
%     c_analyte = c_analyte*(1-1/(op_params{param_pt,6}));
%     c_rep = c_premix(end,end);
%     c_HRP = c_premix(end,1);
%     C0 = [c_analyte c_rep c_HRP 0 c_premix(end,2:end-1) zeros(1,10)]; %Initialize for Sample mixing
%     tspan = [0 op_params{param_pt,5}(2)];
%     
%     [~,C] = ode15s(@(t,c) SampleMixFun(t,c,Soln_pptys(param_pt,6),Soln_pptys(param_pt,7)),tspan,C0,options);
%     C_targets = sum(C(end,9:end)); %sums up all signal generating compounds
%         
% elseif strcmp(MixingFormat,'All mix')
%     if length(op_params{param_pt,5}) ~= 1
%         error("Please input only 1 value for op_params.mixtimes corresponding to the time used for sample mixing");
%     end
%     C0 = [c_analyte Soln_pptys.c_rep Soln_pptys.c_HRP zeros(1,15)];
%     tspan = [0 op_params{param_pt,5}];
%     
%     [~,C] = ode15s(@(t,c) SampleMixFun(t,c,Soln_pptys.k_on_rep,Soln_pptys.k_off_rep),tspan,C0,options);
%     C_targets = sum(C(end,9:end)); %sums up all signal generating compounds      
% end
% 
% C_NS_targets = sum(C(end,[3 5:8])); %sums up all signal generating compounds    

% C_targets = C_targets*1e3; % Everything is already in SI units!

    %Here, we first assume that the loss in diffusivity/bulkiness, on
    %average, balances any increase in k_on, such that k_on_capture
    %sufficiently encompasses the capture of all signal generating
    %compounds
    
%% 1) VFA: Sample application, Solution A application & Wash
% load('Optimized_Params_New.mat'); % use values calibrated for rcSso-SA first

% Target_MW = (Soln_pptys(param_pt,1) + Soln_pptys(param_pt,2) + Soln_pptys(param_pt,3))*1e3; %MW of complex (assumed all to have similar diffusivities
% Paper_params(param_pt,1) = Polson_D_est(Target_MW); %[m2/s] Diffusivity of complex estimated via Polson

% Note: Here we assume sample application or wash application

% Add params where required

% Paper_params should be well formatted for this. May be directly inputted
% into function

%M = 26; N = 90; %Spatial discretization. M = r-direction, N = z-direction;

% Storing values to propagate later
Flow_Seq_store = op_params{param_pt,7};
Vol_Flow_store = op_params{param_pt,2};
Flow_speed_store = op_params{param_pt,1}; 

c_analyte = c_analyte*(1-1/(op_params{param_pt,6})); %Dilution factor of sample

op_params{param_pt,3} = c_analyte.*Flow_Seq_store(1); %initialize sample_conc
op_params{param_pt,2} = Vol_Flow_store(1);
op_params{param_pt,1} = Flow_speed_store(1);


hold_D = Paper_params(param_pt,1);

Paper_params(param_pt,1) = Polson_D_est(Soln_pptys(1)*1e3); %Changing D

rate_params = Soln_pptys(:,[4 5 10 11]); % Pull out values relevant for VFA

rate_params(:,(end-1):end) = 0; %Change the last 2 columns to zero

n_temp = VFA_RUN_fun_variableQ_NSB_parfor_nHRP_only(rate_params,Paper_params,op_params,InletBC,M,N,param_pt);

n_Sample = n_temp(1); %Pull out bound sample

C_targets = sum(c_premix(2:5)); %Pull out all that can bind
C_NS_targets = sum(c_premix(1:5)); % NS targets

% Change CBD to n_Sample
Paper_params(param_pt,5) = n_Sample;
hold_eta = Paper_params(param_pt,6);
Paper_params(param_pt,6) = 1; %Because n_Sample is exposed, would expect it to be all effective?

op_params{param_pt,3} = C_targets.*Flow_Seq_store(2:3); %initialize sample_conc
op_params{param_pt,2} = Vol_Flow_store(2:3);
op_params{param_pt,1} = Flow_speed_store(2:3);

rate_params = Soln_pptys(:,[6 7 10 11]); % Pull out values relevant for VFA
rate_params(:,(end-1):end) = 0; %Change the last 2 columns to zero because they will be computed later

Paper_params(param_pt,1) = hold_D; %change back
n_HRP = VFA_RUN_fun_variableQ_NSB_parfor_nHRP_only(rate_params,Paper_params,op_params,InletBC,M,N,param_pt);

%NSB
rate_params = [Soln_pptys(:,[10 11]).*0 Soln_pptys(:,[10 11])]; % Pull out values relevant for VFA, focus on nsb

% Defining parameters
%L = Paper_params(param_pt,3); %[m] Thickness of cellulose filter paper
R_well = Paper_params(param_pt,4); %[m] Radius of whole well
%V_well = L*pi*R_well^2; %[m3] Volume of whole well
Cellulose_load = (87*pi*R_well^2)/162.14; %[mol] LOAD PER SHEET Cellulose loading/Celluose based on 2.5mm diameter well (Mw = 162.14g/mol, density of paper = 87g/m2)

Paper_params(param_pt,6) = hold_eta; %switch back
op_params{param_pt,3} = C_NS_targets.*Flow_Seq_store(2:3); %Calculate for NSB(
%Paper_params(param_pt,6) = Paper_params(param_pt,8); %Switch back for NSB
%Paper_params(param_pt,5) = Cellulose_load;
n_HRP_temp = VFA_RUN_fun_variableQ_NSB_parfor_nHRP_only(rate_params,Paper_params,op_params,InletBC,M,N,param_pt);
n_HRP(2) = n_HRP(2)+n_HRP_temp(2); %sum up non-specifics

end

%% Code for (0) Premix SA-HRP and BA-MBP-PF2.1

function dcdt = SA_BA_interaction(t,c)
%c = conc vector corresponding to 
% (1) S4. (2) (S4+B) (3) (S4+2B) 
% (4)(S4+3B), (5) (S4+4B)
% (6) Binder-BA
k_on_SA = (4.5e7/1e4)*1e-3; %[M^-1 s^-1] --> convert to SI units, assumed from paper
k_off_SA = 5.2e-6; %[s^-1]

dcdt = zeros(length(c),1);

r1 = (1e4)*k_on_SA*c(1)*c(6) - k_off_SA*c(2); %S4 + B -> (S4+B)
r2 = (1e3)*k_on_SA*c(2)*c(6) - 2*k_off_SA*c(3); %(S4+B) + B -> (S4+2B)
r3 = (1e2)*k_on_SA*c(3)*c(6) - 3*k_off_SA*c(4); %(S4+2B) + B -> (S4+3B)
r4 = (1e1)*k_on_SA*c(4)*c(6) - 4*k_off_SA*c(5); %(S4+3B) + B -> (S4+4B)

dcdt(1) = -r1;
dcdt(2) = r1-r2;
dcdt(3) = r2-r3;
dcdt(4) = r3-r4;
dcdt(5) = r4;
dcdt(6) = -r1-r2-r3-r4;

end

%% Code for (1) Mixing with NP
function dCdt = SampleMixFun(t,C,k_on,k_off) %SA_HRP,NP,Binder)

%This function codes for all formats suggested to compare Conjugated

k_on_Binder = k_on; %[M^-1 s^-1] 3.8e5
k_off_Binder = k_off; %[s^-1] 9.36e-4

k_on_SA = (5.1e6/4)*1e-3; %[M^-1 s^-1] For 1 SA to BA
k_off_SA = 5.1e-9; %[s^-1]

dCdt = zeros(length(C),1);

if length(C) ~= 18
    error('Please input C0 with 18 inputs');
end
% In this format, (1) NP, (2) Binder-BA, (3)SA-HRP
%                 (4) NP-Binder-BA, (5) Binder-BA-SA-HRP, 

%                 (6) Binder-BA(2)SA-HRP, (7)Binder-BA(3)SA-HRP
%                 (8) Binder-BA(4)SA-HRP, (9)NP-Binder-BA(2)SA-HRP

%                 (10) NP-Binder-BA(3)SA-HRP, (11)NP-Binder-BA(4)SA-HRP
%                 (12) NP(2)-Binder-BA(2)SA-HRP, (13)NP(2)-Binder-BA(3)SA-HRP
%                 (14) NP(2)-Binder-BA(4)SA-HRP, (15)NP(3)-Binder-BA(3)SA-HRP
%                 (16) NP(3)-Binder-BA(4)SA-HRP, (17)NP(4)-Binder-BA(4)SA-HRP

%                 (18)NP-Binder-BA-SA-HRP

r1 = k_on_Binder*C(1)*C(2) - k_off_Binder*C(4); % (1) NP + (2) Binder-BA -> (4) NP-Binder-BA
r2 = 4*k_on_SA*C(2)*C(3) - k_off_SA*C(5); % (2) Binder-BA + (3) SA-HRP -> (5) Binder-BA-SA-HRP
r3 = 4*k_on_SA*C(4)*C(3) - k_off_SA*C(18); % (4) NP-Binder-BA + (3) SA-HRP -> (18) NP-Binder-BA-SA-HRP
r4 = k_on_Binder*C(1)*C(5) - k_off_Binder*C(18); % (1) NP + (5) Binder-BA-SA-HRP -> (18) NP-Binder-BA-SA-HRP

r5 = 3*k_on_SA*C(2)*C(5) - k_off_SA*C(6); % (5) Binder-BA-SA-HRP + (2) Binder-BA -> (6) Binder-BA(2)-SA-HRP
r6 = 2*k_on_SA*C(2)*C(6) - k_off_SA*C(7); % (6) Binder-BA(2)-SA-HRP + (2) Binder-BA -> (7) Binder-BA(3)-SA-HRP
r7 = 1*k_on_SA*C(2)*C(7) - k_off_SA*C(8); % (7) Binder-BA(3)-SA-HRP + (2) Binder-BA -> (8) Binder-BA(4)-SA-HRP

r8 = 2*k_on_Binder*C(1)*C(6) - k_off_Binder*C(9); % (1) NP + (6) Binder-BA(2)-SA-HRP -> (9) NP-Binder-BA(2)-SA-HRP
r9 = 3*k_on_Binder*C(1)*C(7) - k_off_Binder*C(10); % (1) NP + (7) Binder-BA(3)-SA-HRP -> (10) NP-Binder-BA(3)-SA-HRP
r10 = 4*k_on_Binder*C(1)*C(8) - k_off_Binder*C(11); % (1) NP + (8) Binder-BA(3)-SA-HRP -> (11) NP-Binder-BA(4)-SA-HRP

r11 = k_on_Binder*C(1)*C(9) - k_off_Binder*C(12); % (1) NP + (9) NP-Binder-BA(2)-SA-HRP -> (12) NP(2)-Binder-BA(2)-SA-HRP
r12 = 2*k_on_Binder*C(1)*C(10) - k_off_Binder*C(13); % (1) NP + (10) NP-Binder-BA(3)-SA-HRP -> (13) NP(2)-Binder-BA(3)-SA-HRP
r13 = 3*k_on_Binder*C(1)*C(11) - k_off_Binder*C(14); % (1) NP + (11) NP-Binder-BA(4)-SA-HRP -> (14) NP(2)-Binder-BA(4)-SA-HRP

r14 = k_on_Binder*C(1)*C(13) - k_off_Binder*C(15); % (1) NP + (13) NP(2)-Binder-BA(3)-SA-HRP -> (15) NP(3)-Binder-BA(3)-SA-HRP
r15 = 2*k_on_Binder*C(1)*C(14) - k_off_Binder*C(16); % (1) NP + (14) NP(2)-Binder-BA(4)-SA-HRP -> (16) NP(3)-Binder-BA(4)-SA-HRP

r16 = k_on_Binder*C(1)*C(16) - k_off_Binder*C(17); % (1) NP + (16) NP(3)-Binder-BA(4)-SA-HRP -> (17) NP(4)-Binder-BA(4)-SA-HRP

r17 = 3*k_on_SA*C(4)*C(5) - k_off_SA*C(9); % (4) NP-Binder-BA + (5) Binder-BA-SA-HRP -> (9) NP-Binder-BA(2)-SA-HRP
r18 = 2*k_on_SA*C(4)*C(6) - k_off_SA*C(10); % (4) NP-Binder-BA + (6) Binder-BA(2)-SA-HRP -> (10) NP-Binder-BA(3)-SA-HRP
r19 = 1*k_on_SA*C(4)*C(7) - k_off_SA*C(11); % (4) NP-Binder-BA + (7) Binder-BA(3)-SA-HRP -> (11) NP-Binder-BA(4)-SA-HRP

r20 = 3*k_on_SA*C(4)*C(18) - k_off_SA*C(12); % (4) NP-Binder-BA + (9) NP-Binder-BA-SA-HRP -> (12) NP(2)-Binder-BA(2)-SA-HRP
r21 = 2*k_on_SA*C(4)*C(9) - k_off_SA*C(13); % (4) NP-Binder-BA + (10) NP-Binder-BA(2)-SA-HRP -> (13) NP(2)-Binder-BA(3)-SA-HRP
r22 = 1*k_on_SA*C(4)*C(10) - k_off_SA*C(14); % (4) NP-Binder-BA + (11) NP-Binder-BA(3)-SA-HRP -> (14) NP(2)-Binder-BA(4)-SA-HRP

r23 = 2*k_on_SA*C(4)*C(12) - k_off_SA*C(15); % (4) NP-Binder-BA + (12) NP(2)-Binder-BA(2)-SA-HRP -> (15) NP(3)-Binder-BA(3)-SA-HRP
r24 = 1*k_on_SA*C(4)*C(13) - k_off_SA*C(16); % (4) NP-Binder-BA + (13) NP(2)-Binder-BA(3)-SA-HRP -> (16) NP(3)-Binder-BA(4)-SA-HRP

r25 = 1*k_on_SA*C(4)*C(15) - k_off_SA*C(17); % (4) NP-Binder-BA + (15) NP(3)-Binder-BA(3)-SA-HRP -> (17) NP(4)-Binder-BA(4)-SA-HRP

r26 = 3*k_on_SA*C(2)*C(18) - k_off_SA*C(9); % (2) Binder-BA + (18) NP-Binder-BA-SA-HRP -> (9) NP-Binder-BA(2)-SA-HRP
r27 = 2*k_on_SA*C(2)*C(9) - k_off_SA*C(10); % (2) Binder-BA + (9) NP-Binder-BA(2)-SA-HRP -> (10) NP-Binder-BA(3)-SA-HRP
r28 = 1*k_on_SA*C(2)*C(10) - k_off_SA*C(11); % (2) Binder-BA + (10) NP-Binder-BA(3)-SA-HRP -> (11) NP-Binder-BA(4)-SA-HRP

r29 = 2*k_on_SA*C(2)*C(12) - k_off_SA*C(13); % (2) Binder-BA + (12) NP(2)-Binder-BA(2)-SA-HRP -> (13) NP(2)-Binder-BA(3)-SA-HRP
r30 = 1*k_on_SA*C(2)*C(13) - k_off_SA*C(14); % (2) Binder-BA + (13) NP(2)-Binder-BA(3)-SA-HRP -> (14) NP(2)-Binder-BA(4)-SA-HRP
r31 = 1*k_on_SA*C(2)*C(15) - k_off_SA*C(16); % (2) Binder-BA + (15) NP(3)-Binder-BA(3)-SA-HRP -> (16) NP(3)-Binder-BA(4)-SA-HRP


dCdt(1) = -r1-r4-r8-r9-r10-r11-r12-r13-r14-r15-r16; %NP

% dCdt(2) = sum(C([2 4 5 18]))+ 2*sum(C([6 9 12])) + 3*sum(C([7 10 13 15])) + 4*sum(C([8 11 14 16 17]))-Binder;%-r1-r2-r5-r6-r7-r26-r27-r28-r29-r30-r31; %Binder-BA
dCdt(2) = -r1-r2-r5-r6-r7-r26-r27-r28-r29-r30-r31; %Binder-BA

dCdt(3) = -r2-r3; %SA-HRP

% dCdt(4) = sum(C([1 4 9:11 18]))+ 2*sum(C(12:14)) + 3*sum(C(15:16)) + 4*C(17) - NP;%r1-r3-r17-r18-r19-r20-r21-r22-r23-r24-r25; %NP-Binder-BA
dCdt(4) = r1-r3-r17-r18-r19-r20-r21-r22-r23-r24-r25; %NP-Binder-BA

% dCdt(5) = sum(C(5:end))+C(3)-SA_HRP;%r2-r4-r5-r17; %Binder-BA-SA-HRP
dCdt(5) = r2-r4-r5-r17; %Binder-BA-SA-HRP

dCdt(6) = r5-r6-r8-r18; %Binder-BA(2)-SA-HRP
dCdt(7) = r6-r7-r9-r19; %Binder-BA(3)-SA-HRP
dCdt(8) = r7-r10; %Binder-BA(4)-SA-HRP

dCdt(9) = r8-r11-r21-r27+r17+r26; %NP-Binder-BA(2)-SA-HRP
dCdt(10) = r9-r12-r22-r28+r18+r27; %NP-Binder-BA(3)-SA-HRP
dCdt(11) = r10-r13+r19+r28; %NP-Binder-BA(4)-SA-HRP

dCdt(12) = r11-r23-r29+r20; %NP(2)-Binder-BA(2)-SA-HRP
dCdt(13) = r12-r14-r24-r30+r29+r21; %NP(2)-Binder-BA(3)-SA-HRP
dCdt(14) = r13-r15+r22+r30; %NP(2)-Binder-BA(4)-SA-HRP

dCdt(15) = r14-r25-r31+r23; %NP(3)-Binder-BA(3)-SA-HRP
dCdt(16) = r15-r16+r31+r24; %NP(3)-Binder-BA(4)-SA-HRP

dCdt(17) = r16+r25; %NP(4)-Binder-BA(4)-SA-HRP

dCdt(18) = r3+r4-r20-r26; %NP-Binder-BA-SA-HRP 

end

%% Code for (2) VFA
function n_TOTAL = VFA_RUN_fun_variableQ_NSB_parfor_nHRP_only(rate_params,Paper_params,op_params,InletBC,M,N,param_pt)
% This function codes for running function for VFA - Assumes start from
% CBD loaded VFA
%
% Inputs: 
%   NOTE: ALL VALUES ARE TO BE IN SI UNITS
%        rate_params = matrix containing kinetic parameters (In columns) 
%                       (1) specific k_on
%                       (2) specific K_D 
%
%                       for each set of parameters (rows)
%        Paper_params = matrix containing paper parameters (In columns)
%                       (1) D, diffusivity (SI units)
%                       (2) R, pore radius
%                       (3) L, Thickness of paper
%                       (4) R_well
%                       (5) CBD, CBD loading
%                       (6) eta, effective CBD accessible
%                       (7) epsilon, porosity of paper
%                       for each set of parameters (rows)
%
%        op_params = cell containing operating parameters (in columns)
%                       (1) Q, Flow rate [m3/s]
%                       (2) Vol, Volume of Sample passed [m3]
%                       (3) SB sample_conc, sample concentration applied on
%                           each step [mol/m3]
%                       (4) NSB sample_conc, sample concentration applied on
%                           each step, of stuff i don't want [mol/m3]
%                       for each set of parameters (rows)
%                       Note that they are keyed in as a vector
%                       corresponding to the number of sample & wash
%                       application steps.
%
%        InletBC = 'Danckwertz' or 'Concentration', Inlet flux or
%                   concentration at z=0 BC
%
%        M,N = spatial mesh points in r and z
%
%        param_pt = set number in question (corresponds to row number in
%                   all parameter inputs)

%% Setup fixed parameters
n_TOTAL = zeros(1,2); %prepare 2 slots

k_on = rate_params(param_pt,1); %[m3/(mol.s)] Rate of association.
k_off = rate_params(param_pt,2); %[s^-1] Rate of dissociation. Can take directly from input
k_on_NSB = rate_params(param_pt,3); %[m3/(mol.s)] Rate of association.
k_off_NSB = rate_params(param_pt,4); %[s^-1] Rate of dissociation. Can take directly from input

D = Paper_params(param_pt,1); %[m2/s] Diffusivity of protein
R = Paper_params(param_pt,2); %[m] Radius of pore, taken at 11um diameter
L = Paper_params(param_pt,3); %[m] Thickness of cellulose filter paper
R_well = Paper_params(param_pt,4); %[m] Radius of whole well

V_well = L*pi*R_well^2; %[m3] Volume of whole well
CBD_load = Paper_params(param_pt,5)*Paper_params(param_pt,6); %[mol], based on Seunghyeon's loading
Area_interaction = V_well*Paper_params(param_pt,7)*(2)/(R); %[m2] SA for interaction, assume bundles of cylindrical pores
Cellulose_load = Paper_params(param_pt,6)*(87*pi*(R_well^2))/162.14; %[mol] LOAD PER SHEET Cellulose loading/Celluose based on 2.5mm diameter well (Mw = 162.14g/mol, density of paper = 87g/m2)

c0 = CBD_load/Area_interaction; %[mol/m2] Concentration in terms of coverage per area

c0_NSB = Cellulose_load/Area_interaction; %[mol/m2] Concentration in terms of coverage per area

No_steps = length(op_params{param_pt,2}); %counts number of flow steps.

rmesh = linspace(0,R,M+1);
zmesh = linspace(0,L,N+1);%
lengthR = M+1;
lengthZ = N+1;

% Sparsity pattern
Sbase = kron( eye(lengthZ),spdiags(ones(lengthR,3),-1:1,lengthR,lengthR));
Sright = zeros(lengthR*lengthZ,2*lengthZ);
Sleft = Sright.';
%Add additional 1s - work along diagonal
for ii = 1:(lengthR*lengthZ)
    if (ii-lengthR>=1) && (ii+lengthR<=(lengthR*lengthZ)) %Somewhere in the middle
        Sbase(ii,[ii-lengthR ii+lengthR]) = 1;
    elseif (ii-lengthR<1) && (ii+lengthR<=(lengthR*lengthZ)) %z = 0 edge
        Sbase(ii,ii+lengthR) = 1;
    elseif (ii-lengthR>=1) && (ii+lengthR>(lengthR*lengthZ)) %z=L edge
        Sbase(ii,ii-[1 2].*lengthR) = 1; %Since i use 2nd order backward difference
    end
    if (rem(ii,lengthR) == 1) && (ii>1)
        Sbase(ii,ii-1) = 0; %delete unnecessary values (r = 0 does not depend on r = R)
    end
    if rem(ii,lengthR) == 0
        Sright(ii,[0 1].*lengthZ+ii/lengthR) = 1; %To account for dependence on all surface species (cS and cSNSB)
        
        if ii < (lengthR*lengthZ)
            Sbase(ii,ii+1) = 0; %(r = R does not depend on r = 0)
        end
    end
end

for ii = 1:lengthZ
    %Take it that all surface species depend on all concentrations
    Sleft((lengthZ.*[0 1]+ii),ii*lengthR) = 1; 
end
Scorner = repmat(eye(lengthZ),2,2); %Just take it that all surface species depend upon each other

S = [Sbase Sright; Sleft Scorner];

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'NonNegative',true,'JPattern',S);


%initial conditions
cSoln0 = zeros(lengthR,lengthZ); %nothing inside at first
cSurf0 = zeros(lengthZ,1);
cNSBSurf0 = zeros(lengthZ,1);

c0_in = [cSoln0(:); cSurf0(:); cNSBSurf0(:)];

for k = 1:No_steps
    cIn = op_params{param_pt,3}(k); %[mol/m3] inlet concentration of target species
    t_S = op_params{param_pt,2}(k)/op_params{param_pt,1}(k); % [s] time of flow
    u_superficial = op_params{param_pt,1}(k)/(pi*(R_well)^2); %[m/s] Speed of flow
    u = u_superficial./Paper_params(param_pt,7); %[m/s] True velocity

    if strcmp(InletBC,'Concentration')
        cSoln0 = reshape(c0_in(1:(lengthZ*lengthR)), [lengthR lengthZ]);
        cSoln0(:,1) = cIn;
        c0_in(1:(lengthZ*lengthR)) = cSoln0(:);
    end
    %% Begin ODE15s setup
    tspan = linspace(0,t_S,100);   
    [times,c]  = ode15s(@VFA_fun,tspan,c0_in,options);

    c_Surf = reshape(c(:,((lengthR*lengthZ+1):(lengthR*lengthZ+lengthZ))).', [lengthZ length(times)]);
    c_NSBSurf = reshape(c(:,((lengthR*lengthZ+lengthZ+1):end)).', [lengthZ length(times)]);
    
    AVG_c_Surf = trapz(zmesh,c_Surf(:,end))/(L);
    AVG_c_NSBSurf= trapz(zmesh,c_NSBSurf(:,end))/(L);
    
    n_TOTAL(1)= AVG_c_Surf*Area_interaction; %SB
    n_TOTAL(2)= AVG_c_NSBSurf*Area_interaction; %NSB

    c0_in = c(end,:).'; %Prepare transfer to next step (if any)
end

    %% MOL function for 2D PDE
    function dcdt = VFA_fun(t,c)
    % This function codes for the Time domain ODEs describing the system which
    % is discretized spatially in 2D
    %
    % Inputs: 
    %        c = all grid points: total: (M+1)x(N+1)+(N+1). The first
    %        (M+1)x(N+1) is the soluble species, cSoln, and the next (N+1)
    %        is the surface species, cS

    % Pull out relevant concentrations for ease
    cSoln = reshape(c(1:(lengthR*lengthZ)), [lengthR lengthZ]);
    cS = c((lengthR*lengthZ+1):(lengthR*lengthZ+lengthZ));
    cSNSB = c((lengthR*lengthZ+lengthZ+1):end);
    
    dcSolndt = zeros(size(cSoln));
    dcSdt = zeros(size(cS));
    dcSNSBdt = zeros(size(cSNSB));
    
    delta_r = R/M; %radial mesh size
    delta_z = L/N; %z mesh size

    %Fill interior points
    for i = 2:M %(2:(lengthR-1))
        for j = 2:N %(2:(lengthZ-1))
            dcSolndt(i,j) = (D/(delta_r^2))*(((1/(2*(i-1)))+1)*cSoln(i+1,j) - 2*cSoln(i,j) + (-(1/(2*(i-1)))+1)*cSoln(i-1,j)) - (u/(2*delta_z))*(cSoln(i,j+1)-cSoln(i,j-1));
        end
    end

    %Fill points along z=0 % z=L
    for i = 2:M %(2:(lengthR-1))
        if strcmp(InletBC,'Danckwertz')
            dcSolndt(i,1) = (D/(delta_r^2))*(((1/(2*(i-1)))+1)*cSoln(i+1,1) - 2*cSoln(i,1) + (-(1/(2*(i-1)))+1)*cSoln(i-1,1)) - (((u)^2)/D)*(cSoln(i,1)-cIn);
        end    
        dcSolndt(i,(N+1)) = (D/(delta_r^2))*(((1/(2*(i-1)))+1)*cSoln(i+1,(N+1)) - 2*cSoln(i,(N+1)) + (-(1/(2*(i-1)))+1)*cSoln(i-1,(N+1))) - (u/(delta_z))*(1.5*cSoln(i,N+1)-2*cSoln(i,N)+0.5*cSoln(i,N-1));
    end

    %fill points along r=0 & r=R, except corner points
    for j=2:N %(2:(lengthZ-1))
        dcSolndt(1,j) = (2*D/(delta_r^2))*(cSoln(2,j)-cSoln(1,j)) - (u/(2*delta_z))*(cSoln(1,j+1)-cSoln(1,j-1));
        
        dcSdt(j) = k_on*(c0-cS(j))*cSoln(M+1,j) - k_off*cS(j);
        dcSNSBdt(j) = k_on_NSB*(c0_NSB-cSNSB(j))*cSoln(M+1,j) - k_off_NSB*cSNSB(j);
        Total_Surf_Rate = dcSdt(j) + dcSNSBdt(j);
        dcSolndt(M+1,j) = (D/(delta_r^2))*(2*(cSoln(M,j)-cSoln(M+1,j)) - (2*delta_r/D)*(1+1/(2*M))*Total_Surf_Rate) - (u/(2*delta_z))*(cSoln(M+1,j+1)-cSoln(M+1,j-1));
    end

    %fill corner points at r=0 & r=R

    %At r=0, z=0 & z=L

    if strcmp(InletBC,'Danckwertz')
        dcSolndt(1,1) = (2*D/(delta_r^2))*(cSoln(2,1)-cSoln(1,1)) - (((u)^2)/D)*(cSoln(1,1)-cIn);
    end
    dcSolndt(1,N+1) = (2*D/(delta_r^2))*(cSoln(2,N+1)-cSoln(1,N+1)) - (u/(delta_z))*(1.5*cSoln(1,N+1)-2*cSoln(1,N)+0.5*cSoln(1,N-1)); %Use backward difference here

    %At r=R, z=0 & z=L
    dcSdt(1) = k_on*(c0-cS(1))*cSoln(M+1,1) - k_off*cS(1);
    dcSNSBdt(1) = k_on_NSB*(c0_NSB-cSNSB(1))*cSoln(M+1,1) - k_off_NSB*cSNSB(1);
    Total_Surf_Rate_1 = dcSdt(1) + dcSNSBdt(1);
    
    if strcmp(InletBC,'Danckwertz')
        dcSolndt(M+1,1) = (D/(delta_r^2))*(2*(cSoln(M,1)-cSoln(M+1,1)) - (2*delta_r/D)*(1+1/(2*M))*Total_Surf_Rate_1) - (((u)^2)/D)*(cSoln(M+1,1)-cIn);    
    end
        
    dcSdt(N+1) = k_on*(c0-cS(N+1))*cSoln(M+1,N+1) - k_off*cS(N+1);
    dcSNSBdt(N+1) = k_on_NSB*(c0_NSB-cSNSB(N+1))*cSoln(M+1,N+1) - k_off_NSB*cSNSB(N+1);
    
    Total_Surf_Rate_N1 = dcSdt(N+1) + dcSNSBdt(N+1);
    dcSolndt(M+1,N+1) = (D/(delta_r^2))*(2*(cSoln(M,N+1)-cSoln(M+1,N+1)) - (2*delta_r/D)*(1+1/(2*M))*Total_Surf_Rate_N1)- (u/(delta_z))*(1.5*cSoln(M+1,N+1)-2*cSoln(M+1,N)+0.5*cSoln(M+1,N-1));
    
    %reshape into column - required for ode15s solver
    dcdt = [dcSolndt(:); dcSdt(:); dcSNSBdt(:)];

    end
end
%% Estimation of Diffusivity by Polson's method

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
