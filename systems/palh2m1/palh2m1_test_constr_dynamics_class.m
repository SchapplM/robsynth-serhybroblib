% Teste Kinematik- und Dynamikfunktionen für palh2m1 mit und ohne Zwangsbedingungen
% Rechnet Kinematikterme von der geschlossenen auf die offene Struktur und
% rechnet Dynamikterme von der offenen auf die geschlossene Struktur um.
% 
% Modell ohne Zwangsbedingungen: S5RRRRR1
% Modell mit Zwangsbedingung: palh2m1

% Quellen:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates
% [UdwadiaKal1992] Udwadia, Firdaus E and Kalaba, Robert E: A new perspective on constrained motion (1992) 
% [Docquier2013] Docquier, Nicolas and Poncelet, Antoine and Fisette, Paul: ROBOTRAN: a powerful symbolic gnerator of multibody models (2013)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Modifikation durch händische Wahl der Parameter und Benutzung der Klassen
% Datum: 2018-12-07 15:09
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

% Moritz Schappler, schappler@irt.uni-hannover.de, 2017-12
% (C) Institut für Regelungstechnik, Universität Hannover

clear
clc

%% Init
% Annahme: Roboter "palh2m1" ist das System mit Zwangsbedingungen, "S5RRRRR1"
% ist das gleiche System ohne Zwangsbedinungen
SName='S5RRRRR1';
RSO = serroblib_create_robot_class(SName);
RSO.descr = 'Palettierroboter 5FG seriell';
serroblib_addtopath({SName});

% Modell des 5FG-Roboters mit einem abhängigen FG als Ersatzmodell
RSC = hybroblib_create_robot_class('palh2m1');
RSC.descr = 'Palettierroboter 5FG seriell, 1FG eliminiert';

TSSo = RSO.gen_testsettings();
TSSc = RSC.gen_testsettings(true); % mit Zwangsbedingungen; zufällige Dynamik-Parameter

% Dynamikparameter angleichen
RSO.DynPar.mges = RSC.DynPar.mges;
RSO.DynPar.rSges = RSC.DynPar.rSges;
RSO.DynPar.mrSges = RSC.DynPar.mrSges;
RSO.DynPar.Icges = RSC.DynPar.Icges;
RSO.DynPar.Ifges = RSC.DynPar.Ifges;
RSO.DynPar.mpv = RSC.DynPar.mpv;

% Kinematik-Parameter
A2Off   = 253e-3;  %[m]
L1      = 1500e-3; %[m]
A3Off   = 0e-3;  %[m]
L2      = 1200e-3; %[m]
A4Off   = 150e-3;  %[m]
pkin = RSO.pkin;
pkin(1) = A2Off; % a2
pkin(2) = sqrt(L1^2+A3Off^2); % a3 % TODO: Stimmt das?
pkin(3) = L2; % a4
pkin(4) = A4Off; % a5
pkin(5) = 0; % d1
pkin(6) = 0; % d5
% Parameter aktualisieren
RSO.update_mdh(pkin);
RSC.update_mdh(pkin);

%% Test der Zwangsbedingungen mit Projektion der Gelenkmomente
% Nach [NakamuraGho1989]

for iq = 1:TSSc.n
  g = TSSc.G(iq,:)';
  
  % Minimalkoordinaten für System mit ZB
  q1 = TSSc.Q(iq,:)';
  qD1 = TSSc.QD(iq,:)';
  qDD1 = TSSc.QDD(iq,:)';

  % Umrechnung von System mit ZB auf offenes System ohne ZB
  W = palh2m1_kinconstr_expl_jacobian_mdh_sym_varpar(q1, RSC.pkin);
  WD = palh2m1_kinconstr_expl_jacobianD_mdh_sym_varpar(q1, qD1, RSC.pkin);
  q = palh2m1_kinconstr_expl_mdh_sym_varpar(q1, RSC.pkin);
  qD = W*qD1;
  qDD = W*qDD1 + WD*qD1;
  
  % Vergleiche Gelenk-Transformationsmatrizen
  Tc = palh2m1_joint_trafo_rotmat_mdh_sym_varpar(q1, RSC.pkin);
  Tu = S5RRRRR1_joint_trafo_rotmat_mdh_sym_varpar(q, RSO.pkin);
  for jj = 1:RSO.NJ 
    test = Tc(:,:,jj)\Tu(:,:,jj) - eye(4);
    if any( abs(test(:)) > 1e-10 )
      [Tc(:,:,jj),Tu(:,:,jj)] %#ok<NOPTS>
      r2eulxyz(Tc(1:3,1:3,jj)' * Tu(1:3,1:3,jj))*180/pi %#ok<NOPTS>
      error('Gelenk-Transformation stimmt nicht');
    end
  end
  % Vergleiche KörperKS-Transformationsmatrizen
  Tcc = palh2m1_fkine_fixb_rotmat_mdh_sym_varpar(q1, RSC.pkin);
  Tcu = S5RRRRR1_fkine_fixb_rotmat_mdh_sym_varpar(q, RSO.pkin);
  for jj = 1:RSC.NL % 1=Basis
    test = Tcc(:,:,jj)\Tcu(:,:,jj) - eye(4);
    if any( abs(test(:)) > 1e-10 )
      [Tcc(:,:,jj),Tcu(:,:,jj)] %#ok<NOPTS>
      r2eulxyz(Tcc(1:3,1:3,jj)' * Tcu(1:3,1:3,jj))*180/pi %#ok<NOPTS>
      error('Transformation stimmt nicht');
    end
  end
  % Vergleiche Geschwindigkeit und Jacobi-Matrizen
  for jj = 1:RSC.NL % 1=Basis
    r_i_i_C = rand(3,1); % Jacobi-Matrix für zufälligen Punkt
    Jgc = palh2m1_jacobig_floatb_twist_sym_varpar(q1, uint8(jj), r_i_i_C, RSC.pkin);
    Jgo = S5RRRRR1_jacobig_floatb_twist_sym_varpar(q, uint8(jj), r_i_i_C, RSO.pkin);
    % Geschwindigkeit des Punktes mit/ohne ZB
    Vc = Jgc*qD1;
    Vo = Jgo*qD;
    if any( abs(Vc-Vo) > 1e-10 )
      error('Geschwindigkeit/Jacobi-Matrix stimmt nicht');
    end
  end

  % Dynamik der offenen Kette
  tau = S5RRRRR1_invdynJ_fixb_slag_vp1(q, qD, qDD, g, ...
    RSO.pkin, RSO.DynPar.mges, RSO.DynPar.rSges, RSO.DynPar.Icges);
  
  % Dynamik der geschlossenen Kette aus Jacobi (Projektion)
  % [NakamuraGho1989], Gl. 5
  tau1 = W' * tau;
  
  % Vergleich mit Lösung mit Elimination
  tau1_test = palh2m1_invdynJ_fixb_slag_vp2(q1, qD1, qDD1, g, ...
    RSC.pkin, RSC.DynPar.mges, RSC.DynPar.mrSges, RSC.DynPar.Ifges);
  
  test = tau1-tau1_test;
  if any( abs(test(:)) > 1e6*eps(max(abs(tau1))) )
    error('Dynamik mit ZB-Jacobi stimmt nicht mit Dynamik aus Eliminations-Ansatz');
  end

end

fprintf('Test der Zwangsbedingungen offen/geschlossen für %d Kombinationen nach [NakamuraGho1989] erfolgreich\n', TSSc.n);

%% Test der Zwangsbedingungen mit Gauß'schem Prinzip des kleinsten Zwangs
% Nach [UdwadiaKal1992]

II = ~logical(RSC.MDH.mu); % Indizes n2 in n
n = RSO.NJ; % Anzahl Gelenke
n1 = RSC.NQJ; % Anzahl unabhängige Gelenke
n2 = sum(II); % Anzahl abhängige Gelenke

for iq = 1:TSSc.n
  g = TSSc.G(iq,:)';
  
  % Minimalkoordinaten für System mit ZB
  q1 = TSSc.Q(iq,:)';
  qD1 = TSSc.QD(iq,:)';

  % Umrechnung von System mit ZB auf offenes System ohne ZB
  W = palh2m1_kinconstr_expl_jacobian_mdh_sym_varpar(q1, RSC.pkin);
  WD = palh2m1_kinconstr_expl_jacobianD_mdh_sym_varpar(q1, qD1, RSC.pkin);
  % Gelenkposition und -geschwindigkeit ohne ZB
  q = palh2m1_kinconstr_expl_mdh_sym_varpar(q1, RSC.pkin);
  qD = W*qD1;

  % Zwangsbedingungen: Bestimme implizite Formulierung nach [Docquier2013]
  % aus expliziter Formulierung (siehe Aufzeichnungen Schappler, 19.12.2017)
  F = W(II,:);
  FD = WD(II,:);
  J = NaN(n2,n); % Aus der Zeitableitung der impliziten ZB nach den Koordinaten des offenen Systems
  JD = NaN(n2,n);
  J(:,~II) = -F;
  J(:,II) = eye(n2,n2);
  JD(:,~II) = -FD;
  JD(:,II) = zeros(n2,n2);
  
  % Zwangsbedingungen in Form aus [UdwadiaKal1992], Gl. (3) konvertieren
  % (mit Koeffizientenvergleich der ZB-Gleichung (siehe [Docquier2013]).
  A = J;
  b = -JD*qD;
  
  % Dynamik ohne ZB ([UdwadiaKal1992], Gl. (1), (2)
  M = S5RRRRR1_inertiaJ_slag_vp1(q, ...
    RSO.pkin, RSO.DynPar.mges, RSO.DynPar.rSges, RSO.DynPar.Icges);
  Q = -S5RRRRR1_invdynJ_fixb_slag_vp1(q, qD, qD*0, g, ...
    RSO.pkin, RSO.DynPar.mges, RSO.DynPar.rSges, RSO.DynPar.Icges);
  a = M \ Q;

  K = sqrtm(M) * pinv( A / sqrtm(M) ); % [UdwadiaKal1992] Text nach Gl. (5a')
  Qc = K*(b - A*(M\Q)); % [UdwadiaKal1992] Gl. (5b)
  
  % Beschleunigung des offenen Systems mit Kräften aus den Zwangsbedingungen
  qDD2 = M\(Q+Qc);
  
  % Referenzlösung: Beschleunigung des geschlossenen Systems (mit ZB)
  M1 = palh2m1_inertiaJ_slag_vp2(q1, ...
    RSC.pkin, RSC.DynPar.mges, RSC.DynPar.mrSges, RSC.DynPar.Ifges);
  tau1 = palh2m1_invdynJ_fixb_slag_vp2(q1, qD1, qD1*0, g, ...
    RSC.pkin, RSC.DynPar.mges, RSC.DynPar.mrSges, RSC.DynPar.Ifges);
  qDD1 = M1 \ (-tau1);
  qDD2_test = W*qDD1 + WD*qD1;
  
  test = qDD2 - qDD2_test;
  if any(abs(test) > 1e-8)
    error('Abgleich nach [UdwadiaKal1992] Fehlgeschlagen');
  end
end

fprintf('Test der Zwangsbedingungen offen/geschlossen für %d Kombinationen nach [UdwadiaKal1992] erfolgreich\n', TSSc.n);
