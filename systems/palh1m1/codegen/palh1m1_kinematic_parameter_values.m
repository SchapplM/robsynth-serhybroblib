AB=2; AM=3; BC=3.5; BE=3; BG=4.5; BL=2; DC=3; EP=3; GP=4;  ML=3.5; T1D=2; T2A=2; T2T1=2;
GH=2;  HW=2; OT2=4; 
% konstanter Winkel fuer VGK PINK
phi312=2*pi/3;
% konstaner Winkel fuer VGK Blau , fuer Rechnung der ZW erforderlich
T1S=((T2T1*T1D)/(T1D+T2A));
T2S=T2T1-T1S;
phi1=atan2(T1S,T1D); phi2=atan2(T2S,T2A);
AS=sqrt(T2A^2+T2S^2);
DS=sqrt(T1S^2+T1D^2);
DA=AS+DS;
% konstaner Winkel fuer VGK GRUEN
BT3=2;phi711=acos(BT3/BC); 
phi413=pi/3; phi710=acos(BT3/BE);
