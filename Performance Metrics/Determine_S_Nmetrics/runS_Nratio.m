Time0 = tic;
run_rcSso_NP_Saliva_SN_variation_Sample_Vol
Time1 = toc(Time0);
sprintf('Time for SN_Sample_Vol = %d sec',Time1)

Time2_0 = tic;
run_rcSso_NP_Saliva_SN_variation_Rep_Conc
Time2 = toc(Time2_0);
sprintf('Time for SN_Rep_Conc = %d sec',Time2)

Time3_0 = tic;
run_rcSso_NP_Saliva_SN_variation_k_on
Time3 = toc(Time3_0);
sprintf('Time for SN_k_on = %d sec',Time3)

Time4_0 = tic;
run_rcSso_NP_Saliva_SN_variation_k_off
Time4 = toc(Time4_0);
sprintf('Time for SN_k_off = %d sec',Time4)

Time5_0 = tic;
run_rcSso_NP_Saliva_SN_variation_IncubTime
Time5 = toc(Time5_0);
sprintf('Time for SN_IncubTime = %d sec',Time5)

Time6_0 = tic;
run_rcSso_NP_Saliva_SN_variation_HRP_Conc
Time6 = toc(Time6_0);
sprintf('Time for SN_HRP_Conc = %d sec',Time6)

Times = [Time1 Time2 Time3 Time4 Time5 Time6];

save('Saliva_SN_RunTimes.mat','Times');
