% Testskript für Kinematik von Planarem Pick&Place-Roboter

% Abderahman Bejaoui, abderahman.bejaoui@outlook.de (Hiwi bei M. Schappler)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

clc
clear 

%% ---------------- Parametereingabe -------------------------------------
L1=650; % Lange Seite Parallelogramm 1: LH, BC
L2=800; % Lange Seite Parallelogramm 2
L3=400; % Gemeinsamer Teil von Fünfgelenkkette und Parallelogramm 2
L4=820; % Seite der Fünfgelenkkette verbunden mit Paralellogramm 2: DE
L5=200; % Kurze Seite Parallelogramm 1 (q1): BH,CL
L6=300; % Kurze Seite Parallelogramm 2 (Endeffektor)
e=170; % Abstand der Antriebe
phi05=120*pi/180; % Winkel Parallelogramm 1 gegen Antriebe
phi1=90*pi/180; % Winkel zwischen Parallelogrammen
% q_a=rand(2,1);
q_a=[10;45]*pi/180;
pkin=[L1;L2;L3;L4;L5;L6;e;phi05;phi1];
%% ------------ Koordinaten-Trafo ----------------------------------------
T=TSRTE_fkine_fixb_rotmat_mdh_sym_varpar(q_a, pkin);
%% Schnittgelenk der F�nfgelenkkette
if ((abs(T(:,:,11)-T(:,:,14))) < 1.0e-3)
    fprintf('Trafos im Schnittgelenk von Fuenfgelenkkette A-B-C-E-D sind richtig \n');
else
    error('Trafos im Schnittgelenk von Fuenfgelenkkette A-B-C-E-D sind nicht richtig');
end
%% Schnittgelenk der Viergelenkkette ( stimmen nicht)
if ((abs(T(:,:,12)-T(:,:,15))) < 1.0e-3)
    fprintf('Trafos im Schnittgelenk von Viergelenkkette B-H-L-C sind richtig \n');
else
    error('Trafos im Schnittgelenk von Viergelenkkette B-H-L-C sind nicht richtig');
end
%% Schnittgelenk der Viergelenkkette ( stimmen nicht)
if ((abs(T(:,:,13)-T(:,:,16))) < 1.0e-3)
    fprintf('Trafos im Schnittgelenk von Viergelenkkette C-G-F-P sind richtig\n');
else
    error('Trafos im Schnittgelenk von Viergelenkkette C-G-F-P sind nicht richtig');
end
%% ----------------- Zuweisung der Koordinaten  ---------------------------
 B =transpose(T(1:3,4,2));
 C =transpose(T(1:3,4,3));
 P=transpose(T(1:3,4,4));  
 E=transpose(T(1:3,4,5)); 
 H=transpose(T(1:3,4,6));
 %C=transpose(T6(1:3,4));
 A=transpose(T(1:3,4,8));
 %C=transpose(T8(1:3,4));
 F=transpose(T(1:3,4,10)); 
 D=transpose(T(1:3,4,11));
 L=transpose(T(1:3,4,12));
 G=transpose(T(1:3,4,13));
 %D=transpose(T13(1:3,4));
 %G=transpose(T14(1:3,4));
 %L=transpose(T15(1:3,4));
%% --------------------- Plot ---------------------
figure(2);clf;  daspect([1 1 1]), axis([-1500, 1000, -1500, 1000, -100 , 100])
 hold on , grid on
 %% ---------------  Gelenke und Punkte  plotten --------------------------
 scatter3(B(1),B(2),B(3),'filled','MarkerFaceColor','k');
 scatter3(C(1),C(2),C(3),'filled','MarkerFaceColor','k');
 scatter3(P(1),P(2),P(3),'filled','MarkerFaceColor','k');
 scatter3(E(1),E(2),E(3),'filled','MarkerFaceColor','k');
 scatter3(H(1),H(2),H(3),'filled','MarkerFaceColor','k');
 scatter3(A(1),A(2),A(3),'filled','MarkerFaceColor','k');  
 scatter3(F(1),F(2),F(3),'filled','MarkerFaceColor','k');
 scatter3(D(1),D(2),D(3),'filled','MarkerFaceColor','k');
 scatter3(G(1),G(2),G(3),'filled','MarkerFaceColor','k');
 scatter3(L(1),L(2),L(3),'filled','MarkerFaceColor','k');
 %% ----------------------- Bennenung der Punkte --------------------------
 strA = 'A';
 A_text = text(A(1)+20,A(2)+100,A(3),strA);
 strB = 'B';
 B_text = text(B(1)+20,B(2)+100,B(3)-20,strB);
 strG = 'G';
 G_text = text(G(1)+30,G(2)+10,G(3),strG);
 strD = 'D';
 D_text = text(D(1)-20,D(2)+100,D(3),strD);
 strC = 'C';
 C_text = text(C(1)+20,C(2)+100,C(3),strC);
 strL = 'L';
 L_text = text(L(1),L(2)+100,L(3)+0.2,strL);
 strE = 'E';
 E_text = text(E(1)+10,E(2)+100,E(3)+30,strE);
 strP = 'P';
 P_text = text(P(1)+20,P(2)+20,P(3),strP);
 strH = 'H';
 H_text = text(H(1)+20,H(2)+100,H(3),strH);
 strF = 'F';
 F_text = text(F(1)+100,F(2)+20,F(3),strF);
 %% ------------- Segemente und Verbindungslinien plotten -----------------
 % serielle Hauptstruktur
  BC=[B;C];
  plot3(BC(:,1),BC(:,2),BC(:,3),'Color','b');
  AB=[A;B];
  plot3(AB(:,1),AB(:,2),AB(:,3),'Color','b');
  CE=[C;E];
  plot3(CE(:,1),CE(:,2),CE(:,3),'Color','g');
  ED=[E;D];
  plot3(ED(:,1),ED(:,2),ED(:,3),'Color','b');
  AD=[A;D];
  plot3(AD(:,1),AD(:,2),AD(:,3),'Color','b'); 
  EP=[E;P];
  plot3(EP(:,1),EP(:,2),EP(:,3),'Color','g');
  FP=[P;F];
  plot3(FP(:,1),FP(:,2),FP(:,3),'Color','g'); 
  FG=[F;G];
  plot3(FG(:,1),FG(:,2),FG(:,3),'Color','g');
  GC=[C;G];
  plot3(GC(:,1),GC(:,2),GC(:,3),'Color','g');
  LC=[L;C];
  plot3(LC(:,1),LC(:,2),LC(:,3),'Color','m');
  LH=[L;H];
  plot3(LH(:,1),LH(:,2),LH(:,3),'Color','m');
  HB=[H;B];
  plot3(HB(:,1),HB(:,2),HB(:,3),'Color','m');
%% -------------------- Koordinatensysteme plotten ------------------------
 %trplot(T(:,:,1),'frame','0','rgb','arrow','color','k','lenght',5,'width',0.5); 
 trplot(T(:,:,2),'frame','1','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,3),'frame','2','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,4),'frame','3','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,5),'frame','4','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,6),'frame','5','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,7),'frame','6','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,8),'frame','7','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,9),'frame','8','rgb','arrow','color','k','lenght',5,'width',0.5); 
 trplot(T(:,:,10),'frame','9','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,11),'frame','10','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,12),'frame','11','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,13),'frame','12','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,14),'frame','13','rgb','arrow','color','k','lenght',5,'width',0.5);
 trplot(T(:,:,15),'frame','14','rgb','arrow','color','k','lenght',10,'width',0.5);
 trplot(T(:,:,16),'frame','15','rgb','arrow','color','k','lenght',10,'width',0.5);
 hold off
 
 xlabel('x');ylabel('y');
 view(0,90)