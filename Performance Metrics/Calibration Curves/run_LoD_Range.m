Time0 = tic;
rcSso_NP_Saliva_LoD_SampleVol
Time1 = toc(Time0);
sprintf('Time for LoD_Range_Sample_Vol = %d sec',Time1)

Time2_0 = tic;
rcSso_NP_Saliva_LoD_RepConc
Time2 = toc(Time2_0);
sprintf('Time for LoD_Range_Rep_Conc = %d sec',Time2)

Time3_0 = tic;
rcSso_NP_Saliva_LoD_k_on
Time3 = toc(Time3_0);
sprintf('Time for LoD_Range_k_on = %d sec',Time3)

Time4_0 = tic;
rcSso_NP_Saliva_LoD_k_off
Time4 = toc(Time4_0);
sprintf('Time for LoD_Range_k_off = %d sec',Time4)

Time5_0 = tic;
rcSso_NP_Saliva_LoD_IncubTime
Time5 = toc(Time5_0);
sprintf('Time for LoD_Range_IncubTime = %d sec',Time5)

Time6_0 = tic;
rcSso_NP_Saliva_LoD_HRP_Conc
Time6 = toc(Time6_0);
sprintf('Time for LoD_Range_HRP_Conc = %d sec',Time6)

Times = [Time1 Time2 Time3 Time4 Time5 Time6];

save('Saliva_LoD_Range_RunTimes.mat','Times');
