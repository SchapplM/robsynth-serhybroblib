%AB=2; BC=3.5; BE=3; BG=4.5; DC=3; EP=3; GP=4; T1A=2; DT2=2; T1T2=2;
%GH=2;  HW=2; OT1=2; 
AB=2; BC=3.5; BE=2; BG=3; DC=4; EP=3; GP=4; T1A=2; DT2=sqrt(3)/2; T1T2=1;
GH=2;  HW=2; OT1=2; 

% konstaner Winkel fuer VGK ABCD , fuer Rechnung der ZW erforderlich
% T1S=((T2T1*T1D)/(T1D+T2A));
% T2S=T2T1-T1S;
% phi1=atan2(T1S,T1D); phi2=atan2(T2S,T2A);
% AS=sqrt(T2A^2+T2S^2);
% DS=sqrt(T1S^2+T1D^2);
% DA=AS+DS;
phi1 = atan2(T1T2,T1A+DT2);       %%hier eigentlich muss phi1 == phi2
phi2 = phi1;
T1S = T1A*tan(phi1);
T2S = T1T2 - T1S;
AS = sqrt(T1A^2+T1S^2);
DS = sqrt(T2S^2+DT2^2);
DA = AS + DS;

% konstaner Winkel fuer VGK BEFG
BT3=2;phi79=acos(BT3/BC); 
phi410=pi/3; phi78=acos(BT3/BE);
