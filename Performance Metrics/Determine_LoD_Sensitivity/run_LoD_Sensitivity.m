%Run LoD
Time0 = tic;
run_rcSso_NP_Saliva_LoD_variation_HRP_Conc
Time1 = toc(Time0);
sprintf('Time for LoD_HRP_Conc = %d sec',Time1)

Time2_0 = tic;
run_rcSso_NP_Saliva_LoD_variation_Sample_Vol
Time2 = toc(Time2_0);
sprintf('Time for LoD_SampleVol = %d sec',Time2)

Time3_0 = tic;
run_rcSso_NP_Saliva_LoD_variation_k_on
Time3 = toc(Time3_0);
sprintf('Time for LoD_k_on = %d sec',Time3)

Time4_0 = tic;
run_rcSso_NP_Saliva_LoD_variation_k_off
Time4 = toc(Time4_0);
sprintf('Time for LoD_k_off = %d sec',Time4)

Time5_0 = tic;
run_rcSso_NP_Saliva_LoD_variation_IncubTime
Time5 = toc(Time5_0);
sprintf('Time for LoD_IncubTime = %d sec',Time5)

Time6_0 = tic;
run_rcSso_NP_Saliva_LoD_variation_Rep_Conc
Time6 = toc(Time6_0);
sprintf('Time for LoD_Rep_Conc = %d sec',Time6)

%Run Sensitivity
Time7_0 = tic;
run_rcSso_NP_Saliva_Sensitivity_variation_HRP_Conc
Time7 = toc(Time7_0);
sprintf('Time for Sens_HRP_Conc = %d sec',Time7)

Time8_0 = tic;
run_rcSso_NP_Saliva_Sensitivity_variation_Sample_Vol
Time8 = toc(Time8_0);
sprintf('Time for Sens_SampleVol = %d sec',Time8)

Time9_0 = tic;
run_rcSso_NP_Saliva_Sensitivity_variation_k_on
Time9 = toc(Time9_0);
sprintf('Time for Sens_k_on = %d sec',Time9)

Time10_0 = tic;
run_rcSso_NP_Saliva_Sensitivity_variation_k_off
Time10 = toc(Time10_0);
sprintf('Time for Sens_k_off = %d sec',Time10)

Time11_0 = tic;
run_rcSso_NP_Saliva_Sensitivity_variation_IncubTime
Time11 = toc(Time11_0);
sprintf('Time for Sens_IncubTime = %d sec',Time11)

Time12_0 = tic;
run_rcSso_NP_Saliva_Sensitivity_variation_Rep_Conc
Time12 = toc(Time12_0);
sprintf('Time for Sens_Rep_Conc = %d sec',Time12)

Times = [Time1 Time2 Time3 Time4 Time5 Time6 Time7 Time8 Time9 Time10 Time11 Time12];

save('Saliva_LoD_Sens_RunTimes.mat','Times');
