function dcdt = Top_Readout_Rxn(t,c,k_on_Rep,k_off_Rep,k_on_CBD,k_off_CBD)
% t = time
% See partner PDF for the numbers of the compounds and reactions
% Rate Constants:
k_on_SA = (5.1e6/4); %[m3 mol^-1 s^-1] 
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

