% This code serves to curve fit all processed data and saves the parameters
% of interest for fitting concentration
close all; clear
load('Processed_Cyan_White_data.mat');

%% Fit with time and concentration

% Note that Column 1 = m_HRP [pg], Column 2 = Rep count,
%           Column 3 = Time [min], Column 4 = Cyan Intensity

beta_guess = [0.65; 0.65*0.01]; % Guess coeff based on previous knowledge: b1*(1-exp(-b2*c*t))

modelfun_all = @(b,x) b(1).*(1-exp(-b(2).*x(:,1).*x(:,2))); %x here contains both time [min] and m_HRP [pg]
%Note .*s which are to ensure everything is multiplied element wise

opt_options = statset('TolFun',1e-12,'TolX',1e-12); %Setting tolerances

X_in = Sorted_Processed_Data(:,[1 3]);
Y_in = Sorted_Processed_Data(:,4);

beta_all = nlinfit(X_in,Y_in,modelfun_all,beta_guess,opt_options);

% Plot fit

conc = linspace(0,200,1000); %From 0 to 200 pg
time = linspace(0,4,400); %From 0 to 4 min
[C,t] = meshgrid(conc,time);
Model_Cyan_Intensity = zeros(size(C));

for i = 1:numel(C)
    Model_Cyan_Intensity(i) = modelfun_all(beta_all,[C(i) t(i)]);
end

f=figure;
s = surf(C,t,Model_Cyan_Intensity); hold on;
s.EdgeColor = 'None';
plot3(Sorted_Processed_Data(:,1),Sorted_Processed_Data(:,3),Sorted_Processed_Data(:,4),'ko','MarkerSize',4,'MarkerFaceColor','k')
xlabel('Amount of label (pg)'); ylabel('Color Development Time (min)'); zlabel('Cyan Intensity')
%lgd = legend('Label Calibration Curve (Saturating exponential)','Experimental Calibration Data');
%title('Cyan Intensity as a function of Label mass and Time')
f.Position(1:2) = f.Position(1:2).*0.75; %115% the original size
f.Position(3:4) = f.Position(3:4).*1.15; %115% the original size
zlim([0 0.5])
ax1 = gca;
ax1.FontSize = 12;
%lgd.FontSize = 12;
x = [0 0 200 200];
y = [1 1 1 1]; %@ 1 min
z = [0 0.8 0.8 0];
patch(x,y,z,'red','FaceAlpha',0.25);



%% Fit @ 1 min

% Note that Column 1 = m_HRP [pg], Column 2 = Rep count,
%           Column 3 = Time [min], Column 4 = Cyan Intensity

beta_guess = beta_all; % Guess coeff based on previous knowledge: b1*(1-exp(-b2*c*t))

modelfun_conc = @(b,x) b(1).*(1-exp(-b(2).*x)); %x here contains only m_HRP [pg]
%Note .*s which are to ensure everything is multiplied element wise

opt_options = statset('TolFun',1e-12,'TolX',1e-12,'MaxIter',1000); %Setting tolerances
options = optimset('TolFun',1e-8,'MaxIter',1e6,'MaxFunEvals',1e6,'Display','iter');

X_in = Sorted_Processed_Data(Sorted_Processed_Data(:,3)==1,1); %only take those at 1 min
Y_in = Sorted_Processed_Data(Sorted_Processed_Data(:,3)==1,4); %only take those at 1 min
A = []; b = []; 
Aeq =[]; beq=[];
lb = [0 0]; ub = [10 0.1];

beta_conc = fmincon(@(b) fit_beta(b,X_in,Y_in),beta_guess,A,b,Aeq,beq,lb,ub,[],options);

%beta_conc = nlinfit(X_in,Y_in,modelfun_conc,beta_guess,opt_options);

% Plot fit

conc = linspace(0,200,1000); %From 0 to 200 pg
Model_Cyan_Intensity_conc = zeros(size(conc));

for i = 1:length(conc)
    Model_Cyan_Intensity_conc(i) = modelfun_conc(beta_conc,conc(i));
end

f=figure; hold on
plot(conc,Model_Cyan_Intensity_conc,'LineWidth',1.5)
plot(X_in,Y_in,'o')
xlabel('Amount of label (pg)'); ylabel('Cyan Intensity')
%title('Cyan Intensity against Label mass at Color Development Time of 1 min')
ylim([0 0.5])
lgd = legend('Label Calibration Curve (Saturating exponential)','Experimental Calibration Data','Location','northwest');
box on
f.Position(1:2) = f.Position(1:2).*0.75; %115% the original size
f.Position(3:4) = f.Position(3:4).*1.15; %115% the original size
ax1 = gca;
ax1.FontSize = 12;
lgd.FontSize = 12;

%Navigate to folder & save
currentFolder = pwd;
Folder_index_start = strfind(currentFolder,'\'); %Find all backslashes
Folder_index_start = Folder_index_start(end); %get the last slash, indicative of where name is

targetFolder = strcat(currentFolder(1:Folder_index_start),'Finalized parameters');
cd(targetFolder)

figure(1)
saveas(gcf,'Saliva_CalibrationCurve_3D.fig')
saveas(gcf,'Saliva_CalibrationCurve_3D.png')
figure(2)
saveas(gcf,'Saliva_CalibrationCurve_1min.fig')
saveas(gcf,'Saliva_CalibrationCurve_1min.png')

cd(currentFolder) %go back

function F = fit_beta(b,X_in,Y_in)
modelfun_conc = @(b,x) b(1).*(1-exp(-b(2).*x)); %x here contains only m_HRP [pg]

for i = 1:length(X_in)
    Y_fit(i) = modelfun_conc(b,X_in(i));
end

%MSE
%F = sum((Y_fit.'-Y_in).^2);

%MAE
F = sum(abs(Y_fit.'-Y_in));

end
