function n_HRP = Top_Readout_Soln_Premix_parfor_nHRP_only(c_analyte,Soln_pptys,Paper_params,op_params,M,N,InletBC,param_pt)
% This function codes for the top readout assay from premixing in a tube,
% followed by loading onto cellulose paper. The first part, premixing, will
% be modeled via a Batch mixing model. The second part, loading, will use a
% simple efficiency fitted from the premix model (~53% capture)
%
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
%                       (12) c_CBD [mol/m3]
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
%
%        op_params = cell containing VFA operating parameters (in columns)
%                       (1) Q, Flow rate [m3/s]
%                       (2) Vol, Volume of Sample passed [m3]
%                       (3) sample_conc, sample concentration applied on
%                           each step [mol/m3]
%                       (4) MixingFormat 
%                           "All mix" or "Premix", indicative of whether all
%                           are added together and mixed, or premix, whereby SA-HRP and
%                           BA-MBP-PF2.1 are added together first to give Solution A, followed
%                           by application.
%                       (5) Sample Mixing Times [s] (1x2 vector if using
%                           "Premix"
%                           for each set of parameters (rows)
%                           Note that they are keyed in as a vector
%                           corresponding to the number of sample & wash
%                           application steps.
%                       (6) Dilution of Premix if using
%                       (7) SampleORWash: Check if loading sample or
%                           washing
%                       (8) CaptureEfficiency, Percentage of total target
%                           captured
%
%        param_pt = set number in question (corresponds to row number in
%                   all parameter inputs)
% Output: 
%        n_HRP = Number of moles of HRP captured on paper

%% Checking important inputs
MixingFormat = op_params{param_pt,4};

%% 0) Premix BA-MBP-Rep and SA-HRP
% Use only if requested. 
% Note: There is a dilution depending on the operating parameters premix and sample
options = odeset('AbsTol',1e-14,'RelTol',1e-8,'NonNegative',true); %Accuracy of ODE solver

if strcmp(MixingFormat,'Premix') 
    if length(op_params{param_pt,5}) ~= 2
        error("Please input a vector for op_params{?,5} of 2 times corresponding to the times used for premixing and sample mixing");
    end
    c0 = zeros(1,41);
    c0([2 4 19]) = [Soln_pptys(param_pt,8) Soln_pptys(param_pt,9) Soln_pptys(param_pt,12)].*op_params{param_pt,6}; % Premixing concentration Dilution X of final
    tspan = [0 op_params{param_pt,5}(1)];
    [~,c] = ode15s(@(t,c) Top_Readout_Rxn(t,c,Soln_pptys(param_pt,6),Soln_pptys(param_pt,7),Soln_pptys(param_pt,6),Soln_pptys(param_pt,7)),tspan,c0,options); %Binding parameters already in code
    c_premix = c(end,:)./op_params{param_pt,6}; %to transfer to sample mixing step
end

%% 1) Mix Solution A with Analyte (NP)

if strcmp(MixingFormat,'Premix')
    c_analyte = c_analyte*(1-1/(op_params{param_pt,6}));
    C0 = c_premix; %Initialize for Sample mixing
    C0(1) = c_analyte;
    tspan = [0 op_params{param_pt,5}(2)];
    
    [~,C] = ode15s(@(t,c) Top_Readout_Rxn(t,c,Soln_pptys(param_pt,6),Soln_pptys(param_pt,7),Soln_pptys(param_pt,6),Soln_pptys(param_pt,7)),tspan,C0,options);
    C_targets = sum(C(end,22:41)); %sums up all signal generating compounds
        
elseif strcmp(MixingFormat,'All mix')
    if length(op_params{param_pt,5}) ~= 1
        error("Please input only 1 value for op_params.mixtimes corresponding to the time used for sample mixing");
    end
    C0 = zeros(1,41);
    C0([1 2 4 19]) = [c_analyte Soln_pptys(param_pt,8) Soln_pptys(param_pt,9) Soln_pptys(param_pt,12)].*op_params{param_pt,6}; % Premixing concentration Dilution X of final
    tspan = [0 op_params{param_pt,5}];
    
    [~,C] = ode15s(@(t,c) Top_Readout_Rxn(t,c,Soln_pptys(param_pt,6),Soln_pptys(param_pt,7),Soln_pptys(param_pt,6),Soln_pptys(param_pt,7)),tspan,C0,options);
    C_targets = sum(C(end,22:41)); %sums up all signal generating compounds      
end

C_NS_targets = sum(C(end,4:18)); %sums up all signal generating compounds    

% C_targets = C_targets*1e3; % Everything is already in SI units!

    %Here, we first assume that the loss in diffusivity/bulkiness, on
    %average, balances any increase in k_on, such that k_on_capture
    %sufficiently encompasses the capture of all signal generating
    %compounds
    
%% 2) VFA: Sample application & Wash

%M = 26; N = 90; %Spatial discretization. M = r-direction, N = z-direction;

%op_params{param_pt,3} = C_targets.*op_params{param_pt,7}; %initialize sample_conc

%rate_params = Soln_pptys(:,[4 5 10 11]); % Pull out values relevant for VFA
n_HRP = zeros(1,2);

n_HRP(1) = C_targets.*op_params{param_pt,8}.*op_params{param_pt,2}(1); % Conc*CaptureEfficiency.*SampleVol

%NSB

%C_targets_left =C_targets.*(1-op_params{param_pt,8}).*op_params{param_pt,2};
C_targets_left = C_targets.*(1-op_params{param_pt,8});

op_params{param_pt,3} = (C_targets_left + C_NS_targets).*op_params{param_pt,7}; %initialize sample_conc
rate_params = [Soln_pptys(:,[10 11]).*0 Soln_pptys(:,10) Soln_pptys(:,11)]; % Pull out values relevant for VFA, focus on nsb
%Fudge factor of 1.225 because of more NSB?

% Defining parameters
%L = Paper_params(param_pt,3); %[m] Thickness of cellulose filter paper
R_well = Paper_params(param_pt,4); %[m] Radius of whole well
%V_well = L*pi*R_well^2; %[m3] Volume of whole well
Cellulose_load = (87*pi*R_well^2)/162.14; %[mol] LOAD PER SHEET Cellulose loading/Celluose based on 2.5mm diameter well (Mw = 162.14g/mol, density of paper = 87g/m2)

%Paper_params(param_pt,5) = Cellulose_load;
n_HRP_temp = VFA_RUN_fun_variableQ_NSB_parfor_nHRP_only(rate_params,Paper_params,op_params,InletBC,M,N,param_pt);
if n_HRP(1) == 0 && c_analyte~=0
    n_HRP(2) = 0;
else
    n_HRP(2) = n_HRP(2)+n_HRP_temp(2); %sum up non-specifics
end

end

%% Code for (0) Premix SA-HRP and BA-MBP-PF2.1

function dcdt = SA_BA_interaction(t,c)
%c = conc vector corresponding to 
% (1) S4. (2) (S4+B) (3) (S4+2B) 
% (4)(S4+3B), (5) (S4+4B)
% (6) Binder-BA
k_on_SA = (5.1e6/4)*1e-3; %[M^-1 s^-1] --> convert to SI units, assumed from paper
k_off_SA = 5.1e-9; %[s^-1]

dcdt = zeros(length(c),1);

r1 = 4*k_on_SA*c(1)*c(6) - k_off_SA*c(2); %S4 + B -> (S4+B)
r2 = 3*k_on_SA*c(2)*c(6) - k_off_SA*c(3); %(S4+B) + B -> (S4+2B)
r3 = 2*k_on_SA*c(3)*c(6) - k_off_SA*c(4); %(S4+2B) + B -> (S4+3B)
r4 = k_on_SA*c(4)*c(6) - k_off_SA*c(5); %(S4+3B) + B -> (S4+4B)

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
%% Code for (1) Mixing with NP
function dcdt = Top_Readout_Rxn(t,c,k_on_Rep,k_off_Rep,k_on_CBD,k_off_CBD)
% t = time
% See partner PDF for the numbers of the compounds and reactions
% Rate Constants:
k_on_SA = (5.1e6/4)*1e-3; %[m3 mol^-1 s^-1] 
k_off_SA = 5.1e-9; %[s^-1]

% Rates of Reactions as numbered in PDF:

% Rep Reactions
r1 = k_on_Rep*c(1)*c(2) - k_off_Rep*c(3);

r2 = 4*k_on_SA*c(2)*c(4) - k_off_SA*c(5);

r3 = 4*k_on_SA*c(3)*c(4) - k_off_SA*c(6);
r4 = k_on_Rep*c(1)*c(5) - k_off_Rep*c(6);

r5 = 3*k_on_SA*c(2)*c(5) - 2*k_off_SA*c(7);

r6 = 2*k_on_SA*c(2)*c(7) - 3*k_off_SA*c(8);

r7 = k_on_SA*c(2)*c(8) - 4*k_off_SA*c(9);

r8 = 2*k_on_Rep*c(1)*c(7) - k_off_Rep*c(10);
r9 = 3*k_on_SA*c(6)*c(2) - k_off_SA*c(10);
r10 = 3*k_on_SA*c(3)*c(5) - k_off_SA*c(10);

r11 = 3*k_on_Rep*c(1)*c(8) - k_off_Rep*c(11);
r12 = 2*k_on_SA*c(2)*c(10) - 2*k_off_SA*c(11);
r13 = 2*k_on_SA*c(3)*c(7) - k_off_SA*c(11);

r14 = 4*k_on_Rep*c(1)*c(9) - k_off_Rep*c(12);
r15 = k_on_SA*c(2)*c(11) - 3*k_off_SA*c(12);
r16 = k_on_SA*c(3)*c(8) - k_off_SA*c(12);

r17 = k_on_Rep*c(1)*c(10) - 2*k_off_Rep*c(13);
r18 = 3*k_on_SA*c(3)*c(6) - 2*k_off_SA*c(13);

r19 = 2*k_on_Rep*c(1)*c(11) - 2*k_off_Rep*c(14);
r20 = 2*k_on_SA*c(3)*c(10) - 2*k_off_SA*c(14);
r21 = 2*k_on_SA*c(2)*c(13) - k_off_SA*c(14);

r22 = 3*k_on_Rep*c(1)*c(12)- 2*k_off_Rep*c(15);
r23 = k_on_SA*c(3)*c(11) - 2*k_off_SA*c(15);
r24 = k_on_SA*c(2)*c(14) - k_off_SA*c(15);

r25 = k_on_Rep*c(1)*c(14) - 3*k_off_Rep*c(16);
r26 = 2*k_on_SA*c(3)*c(13) - 3*k_off_SA*c(16);

r27 = k_on_Rep*c(1)*c(18) - 4*k_off_Rep*c(17);
r28 = k_on_SA*c(3)*c(16) - 4*k_off_SA*c(17);

r29 = 2*k_on_Rep*c(1)*c(15) - 3*k_off_Rep*c(18);
r30 = k_on_SA*c(3)*c(14) - 3*k_off_SA*c(18);
r31 = k_on_SA*c(2)*c(16) - k_off_SA*c(18);

% CBD Reactions
r32 = k_on_CBD*c(1)*c(19) - k_off_CBD*c(20);

r33 = k_on_CBD*c(19)*c(3) - k_off_CBD*c(21);
r34 = k_on_Rep*c(20)*c(2) - k_off_Rep*c(21);

r35 = k_on_CBD*c(19)*c(6) - k_off_CBD*c(22);
r36 = k_on_Rep*c(20)*c(5) - k_off_Rep*c(22);
r37 = 4*k_on_SA*c(4)*c(21) - k_off_SA*c(22);

r38 = k_on_CBD*c(19)*c(10) - k_off_CBD*c(23);
r39 = 2*k_on_Rep*c(20)*c(7) - k_off_Rep*c(23);
r40 = 3*k_on_SA*c(21)*c(5) - k_off_SA*c(23);
r41 = 3*k_on_SA*c(2)*c(22) - k_off_SA*c(23);

r42 = k_on_CBD*c(19)*c(11) - k_off_CBD*c(24);
r43 = 3*k_on_Rep*c(20)*c(8) - k_off_Rep*c(24);
r44 = 2*k_on_SA*c(21)*c(7) - k_off_SA*c(24);
r45 = 2*k_on_SA*c(23)*c(2) - 2*k_off_SA*c(24);

r46 = k_on_CBD*c(19)*c(12) - k_off_CBD*c(25);

r47 = 4*k_on_Rep*c(20)*c(9) - k_off_Rep*c(25);

r48 = k_on_SA*c(21)*c(8) - k_off_SA*c(25);
r49 = k_on_SA*c(24)*c(2) - k_off_SA*c(25);

r50 = k_on_CBD*c(19)*c(13) - k_off_CBD*c(26);
r51 = 3*k_on_Rep*c(20)*c(10) - k_off_Rep*c(26);
r52 = 3*k_on_SA*c(21)*c(6) - k_off_SA*c(26);
r53 = 3&k_on_SA*c(22)*c(3) - k_off_SA*c(26);
r54 = k_on_Rep*c(23)*c(1) - k_off_Rep*c(26);

r55 = k_on_CBD*c(19)*c(26) - k_off_CBD*c(27);
r56 = k_on_Rep*c(20)*c(23) - 2*k_off_Rep*c(27);
r57 = 3*k_on_SA*c(21)*c(22) - 2*k_off_SA*c(27);

r58 = k_on_CBD*c(19)*c(14) - k_off_CBD*c(28);
r59 = 2*k_on_Rep*c(20)*c(11) - k_off_Rep*c(28);
r60 = 2*k_on_SA*c(21)*c(10) - k_off_SA*c(28);
r61 = 2*k_on_SA*c(23)*c(3) - k_off_SA*c(28);
r62 = 2*k_on_Rep*c(24)*c(1) - k_off_Rep*c(28);

r63 = k_on_CBD*c(19)*c(28) - 2*k_off_CBD*c(29);
r64 = 2*k_on_Rep*c(20)*c(24) - 2*k_off_Rep*c(29);
r65 = 2*k_on_SA*c(21)*c(23) - 2*k_off_SA*c(29);
r66 = 2*k_on_SA*c(27)*c(2) - k_off_SA*c(29);

r67 = 2*k_on_CBD*c(19)*c(15) - k_off_CBD*c(30);
r68 = 3*k_on_Rep*c(20)*c(12) - k_off_Rep*c(30);
r69 = k_on_SA*c(21)*c(11) - k_off_SA*c(30);
r70 = k_on_SA*c(24)*c(3) - k_off_SA*c(30);
r71 = 3*k_on_Rep*c(25)*c(1) - k_off_Rep*c(30);
r72 = k_on_SA*c(28)*c(2) - k_off_SA*c(30);

r73 = k_on_CBD*c(19)*c(28) - 2*k_off_CBD*c(31);
r74 = 3*k_on_Rep*c(20)*c(25) - 2*k_off_Rep*c(31);
r75 = k_on_SA*c(21)*c(24) - 2*k_off_SA*c(31);
r76 = k_on_SA*c(29)*c(2) - 2*k_off_SA*c(31);

r77 = 3*k_on_CBD*c(19)*c(16) - k_off_CBD*c(32);
r78 = k_on_Rep*c(20)*c(14) - k_off_Rep*c(32);
r79 = 2*k_on_SA*c(21)*c(13) - k_off_SA*c(32);
r80 = 2*k_on_SA*c(26)*c(3) - 2*k_off_SA*c(32);
r81 = k_on_Rep*c(28)*c(1) - 2*k_off_Rep*c(32);

r82 = 2*k_on_Rep*c(19)*c(32) - 2*k_off_Rep*c(33);
r83 = k_on_Rep*c(20)*c(28) - 2*k_off_Rep*c(33);
r84 = 2*k_on_SA*c(21)*c(26) - 2*k_off_SA*c(33);
r85 = 2*k_on_SA*c(27)*c(3) - k_off_SA*c(33);
r86 = k_on_Rep*c(29)*c(1) - k_off_Rep*c(33);

r87 = k_on_CBD*c(19)*c(33) - 3*k_off_CBD*c(34);
r88 = k_on_Rep*c(20)*c(29) - 3*k_off_Rep*c(34);
r89 = 2*k_on_SA*c(21)*c(27) - 3*k_off_SA*c(34);

r90 = k_on_CBD*c(19)*c(17) - k_off_CBD*c(35);
r91 = k_on_Rep*c(20)*c(18) - k_off_Rep*c(35);
r92 = k_on_SA*c(21)*c(16) - k_off_SA*c(35);
r93 = k_on_SA*c(32)*c(3) - 3*k_off_SA*c(35);
r94 = k_on_Rep*c(39)*c(1) - 3*k_off_Rep*c(35);

r95 = k_on_CBD*c(19)*c(35) - k_off_CBD*c(36);
r96 = k_on_Rep*c(20)*c(39) - 2*k_off_Rep*c(36);
r97 = k_on_SA*c(21)*c(32) - 2*k_off_SA*c(36);
r98 = k_on_SA*c(33)*c(3) - 2*k_off_SA*c(36);
r99 = k_on_Rep*c(40)*c(1) - 2*k_off_Rep*c(36);

r100 = 2*k_on_CBD*c(19)*c(36) - 3*k_off_CBD*c(37);
r101 = k_on_Rep*c(20)*c(40) - 3*k_off_Rep*c(37);
r102 = k_on_SA*c(21)*c(33) - 3*k_off_SA*c(37);
r103 = k_on_SA*c(34)*c(3) - k_off_SA*c(37);
r104 = k_on_Rep*c(41)*c(1) - k_off_Rep*c(37);

r105 = k_on_CBD*c(19)*c(37) - 4*k_off_CBD*c(38);
r106 = k_on_Rep*c(20)*c(41) - 4*k_off_Rep*c(38);
r107 = k_on_SA*c(21)*c(34) - 4*k_off_SA*c(38);

r108 = 2*k_on_SA*c(23)*c(2) - k_off_SA*c(28);

r109 = 2*k_on_CBD*c(19)*c(18) - k_off_CBD*c(39);
r110 = 2*k_on_Rep*c(20)*c(15) - k_off_Rep*c(39);
r111 = k_on_SA*c(21)*c(14) - k_off_SA*c(39);
r112 = k_on_SA*c(32)*c(2) -k_off_SA*c(39);
r113 = k_on_SA*c(28)*c(3) - 2*k_off_SA*c(39);
r114 = 2*k_on_Rep*c(30)*c(1) - 2*k_off_Rep*c(39);

r115 = 2*k_on_CBD*c(19)*c(39) - 2*k_off_CBD*c(40);
r116 = 2*k_on_Rep*c(20)*c(30) - 2*k_off_Rep*c(40);
r117 = k_on_SA*c(21)*c(28) - 2*k_off_SA*c(40);
r118 = k_on_SA*c(33)*c(2) - k_off_SA*c(40);
r119 = k_on_SA*c(29)*c(3) - k_off_SA*c(40);
r120 = 2*k_on_Rep*c(31)*c(1) - k_off_Rep*c(40);

r121 = k_on_CBD*c(19)*c(40) - 3*k_off_CBD*c(41);
r122 = 2*k_on_Rep*c(20)*c(31) - 3*k_off_Rep*c(41);
r123 = k_on_SA*c(21)*c(29) - 3*k_off_SA*c(41);
r124 = k_on_SA*c(34)*c(2) - k_off_SA*c(41);

% Concentrations of Molecules (See partner PDF for numbers): 
dcdt = zeros(41,1);
dcdt(1) = -r1 -r4 -r8 -r11 -r14 -r17 -r19 -r22 -r25 -r27 -r29 -r32 -r54 -r62 -r71 -r81 -r86 -r94 -r99 -r104 -r114 -r120; %N
dcdt(2) = -r1 -r2 -r5 -r6 -r7 -r9 -r12 -r15 -r21 -r24 -r31 -r34 -r41 -r45 -r49 -r108 -r66 -r72 -r76 -r112 -r118 -r124; %BA (MBP)
dcdt(3) = r1 -r3 -r10 -r13 -r16 -r18 -r20 -r23 -r26 -r28 -r30 -r33 -r53 -r61 -r70 -r80 -r85 -r93 -r98 -r103 -r113 -r119; % N:BA
dcdt(4) = -r2 -r3 -r37; %HRP
dcdt(5) = r2 -r4 -r5 -r10 -r36 -r40; %BA:HRP
dcdt(6) = r3 +r4 -r9 -r18 -r35 -r52; %N:BA:HRP
dcdt(7) = r5 -r6 -r8 -r13 -r39 -r44; %BA2:HRP
dcdt(8) = r6 -r7 -r11 -r16 -r43 -r48; %BA3:HRP
dcdt(9) = r7 -r14 -r47; %BA4:HRP
dcdt(10) = r8 +r9 +r10 -r12 -r17 -r20 -r38 -r51 -r60; %N:BA2:HRP
dcdt(11) = r11 +r12 +r13 -r15 -r19 -r23 -r42 -r59 -r69; %N:BA3:HRP
dcdt(12) = r14 +r15 +r16 -r22 -r46 -r68; %N:BA4:HRP
dcdt(13) = r17 +r18 -r21 -r26 -r50 -r79; %N2:BA2:HRP
dcdt(14) = r19 +r20 +r21 -r24 -r25 -r30 -r58 -r78 -r111; %N2:BA3:HRP
dcdt(15) = r22 +r23 +r24 -r29 -r67 -r110; %N2:BA4:HRP
dcdt(16) = r25 +r26 -r28 -r31 -r77 -r92; %N3:BA3:HRP
dcdt(17) = r27 +r28 -r90; %N4:BA4:HRP
dcdt(18) = r29 +r30 +r31 -r27 -r91 -r109; %N3:BA4:HRP
dcdt(19) = -r32 -r33 -r35 -r38 -r42 -r46 -r50 -r55 -r58 -r63 -r67 -r73 -r77 -r82 -r87 -r90 -r95 -r100 -r105 -r109 -r115 -r121; %CBD
dcdt(20) =  r32 -r34 -r36 -r39 -r43 -r47 -r51 -r56 -r59 -r64 -r68 -r74 -r78 -r83 -r88 -r91 -r96 -r101 -r106 -r110 -r116 -r122; %CBD:N
dcdt(21) =  r33 +r34 -r37 -r40 -r44 -r48 -r52 -r57 -r60 -r65 -r69 -r75 -r79 -r84 -r89 -r92 -r97 -r102 -r107 -r111 -r117 -r123; %CBD:N:BA
dcdt(22) = r35 +r36 +r37 -r41 -r53 -r57; %CBD:N:BA:HRP
dcdt(23) = r38 +r39 +r40 +r41 -r45 -r54 -r56 -r61 -r65 -r108;  %CBD:N:BA2:HRP
dcdt(24) = r42 +r43 +r44 +r45 -r49 -r62 -r64 -r70 -r75; %CBD:N:BA3:HRP
dcdt(25) = r46 +r47 +r48 +r49 -r71 -r74; %CBD:N:BA4:HRP
dcdt(26) = r50 +r51 +r52 +r53 +r54 -r55 -r80 -r84; %CBD:N2:BA2:HRP
dcdt(27) = r55 +r61 +r57 -r66 -r85 -r89; %CBD2:N2:BA2:HRP
dcdt(28) = r58 +r59 +r60 +r61 +r62 +r108 -r63 -r72 -r73 -r81 -r83 -r113 -r117; %CBD:N2:BA3:HRP
dcdt(29) = r63 +r64 +r65 +r66 -r76 -r86 -r88 -r119 -r123; %CBD2:N2:BA3:HRP
dcdt(30) = r67 +r68 +r69 +r70 +r71 +r72 -r114 -r116; %CBD:N2:BA4:HRP
dcdt(31) = r73 +r74 +r75 +r76 -r120 -r122; %CBD2:N2:BA4:HRP
dcdt(32) = r77 +r78 +r79 +r80 +r81 -r82 -r93 -r97 -r112; %CBD:N3:BA3:HRP
dcdt(33) = r82 +r83 +r84 +r85 +r86 -r87 -r98 -r102 -r118; %CBD2:N3:BA3:HRP
dcdt(34) = r87 +r88 +r89 -r95 -r103 -r107 -r124; %CBD3:N3:BA3:HRP
dcdt(35) = r90 + r91 +r92 +r93 +r94 -r95; %CBD:N4:BA4:HRP
dcdt(36) = r95 + r96 +r97 +r98 +r99 -r100; %CBD2:N4:BA4:HRP
dcdt(37) = r100 +r101 +r102 +r103 +r104 -r105; %CBD3:N4:BA4:HRP
dcdt(38) = r105 +r106 +r107; %CBD4:N4:BA4:HRP
dcdt(39) = r109 +r110 +r111 +r112 +r113 +r114 -r94 -r96 -r115; %CBD:N3:BA4:HRP
dcdt(40) = r115 +r116 +r117 +r118 +r119 +r120 -r99 -r101 -r121; %CBD2:N3:BA4:HRP
dcdt(41) = r121 +r122 +r123 +r124 -r104 -r106; %CBD3:N3:BA4:HRP
 
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

