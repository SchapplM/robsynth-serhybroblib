% Calculate joint inertia based on Composite Rigid Body Algorithm (CRBA) from [Featherstone2008]
%
% Input:
% q [6x1]
%   Joint Angles [rad]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% H [6x6]
%   Inertia Matrix
%
% Sources:
% [Featherstone2008] Roy Featherstone: Rigid Body Dynamics Algorithms (2008)

% TODO: Funktioniert noch nicht für Schubgelenke

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-12
% (C) Institut für Regelungstechnik, Universität Hannover

function H = palh2m2OL_inertiaJ_nCRB_vp1(q, pkin, m_num, rSges_num_mdh, Icges_num_mdh)
%%Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(q) && all(size(q) == [6 1]), ...
  'palh2m2OL_inertiaJ_nCRB_vp1: q has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaJ_nCRB_vp1: Kinematic parameters pkin have to be [5x1] (double)');
assert(isreal(m_num) && all(size(m_num) == [7 1]), ...
  'palh2m2OL_inertiaJ_nCRB_vp1: m_num has to be [7x1] (double)');
assert(isreal(rSges_num_mdh) && all(size(rSges_num_mdh) == [7,3]), ...
  'palh2m2OL_inertiaJ_nCRB_vp1: rSges_num_mdh has to be [7x3] (double)');
assert(isreal(Icges_num_mdh) && all(size(Icges_num_mdh) == [7 6]), ...
  'palh2m2OL_inertiaJ_nCRB_vp1: Icges_num_mdh has to be [7x6] (double)');

%% Init
[v_mdh, sigma] = palh2m2OL_structural_kinematic_parameters();

%% Allgemein
% Vorgänger-Segmente der Gelenke, [Featherstone2008] S. 71
% p = v_mdh; % 0 ist die Basis (pelvis)
% Nachfolger-Segmente der Gelenke
% s = (1:6)';
% parent array (Vorgänger-Segment für alle Segmente)
% [Featherstone2008] Gl. (4.2)
% Hier ist lambda identisch mit v_mdh
% lambda = min([p,s]')'; %#ok<UDIM>
lambda = v_mdh;

%% Kinematik
T_mdh = palh2m2OL_joint_trafo_rotmat_mdh_sym_varpar(q, pkin);

%% Massenmatrix zusammensetzen
H = zeros(6,6);
% compose matrix for all links

I_c_ges = NaN(6,36);
for i = 1:6 % loop through all bodies, [Featherstone2008] T6.2, Line 2
  % Trägheitsmatrix für diesen Körper
  I_c = inertiavector2matrix(Icges_num_mdh(i+1,:)); % Andere Indexnotation beachten
  % [Featherstone2008] Gl. (2.63)
  I_O = [I_c+m_num(i+1)*skew(rSges_num_mdh(i+1,:))*skew(rSges_num_mdh(i+1,:))', m_num(i+1)*skew(rSges_num_mdh(i+1,:)); ...
         m_num(i+1)*skew(rSges_num_mdh(i+1,:))',                                m_num(i+1)*eye(3)];

  % Komposit-Matrix initialisieren
  I_c_ges(i,:) = I_O(:);
end

% loop through all bodies, [Featherstone2008] T6.2, Line 4
% Basis zählt nicht als Körper
for i = 6:-1:1 % Wird die Schleife mit uint8 definiert, funktioniert der generierte mex-Code nicht (Matlab-Bug)
  % Alt. 1. Berechnung mit kumulativen Transformationsmatrizen (intuitiver):
  % R_0_i = T_c_mdh(1:3,1:3,i+1); % Andere Indexnotation beachten. Matlab: 1=Basis, Algorithmus: 0=Basis
  % % Transformation zum Vorgänger-Segment
  % R_0_li = T_c_mdh(1:3,1:3,lambda(i)+1);
  % R_i_li = R_0_i' * R_0_li;
  % r_i_i_li = R_0_i' * r_0_i_li;
  % r_0_i_li = -T_c_mdh(1:3,4,i+1) + T_c_mdh(1:3,4,lambda(i)+1);

  % Alt. 2. Berechnung mit direkten Transformationsmatrizen (schneller):
  % Transformation von diesem zum Vorgänger
  R_i_li = T_mdh(1:3,1:3,i)';
  % Vektor von diesem zum Vorgänger
  r_i_i_li = -R_i_li*T_mdh(1:3,4,i);

  % [Featherstone2008], S. 22
  %   A = i
  %   B = li
  %   r = r_i_i_li
  %   E = R_B_A = R_li_i = R_i_li'
  %  Gl. (2.26)
  X_i_li = [R_i_li, zeros(3,3);skew(r_i_i_li)*R_i_li , R_i_li];

  I_i_c = reshape(I_c_ges(i,:), 6, 6);

  % Komposit-Matrix aufbauen
  if lambda(i) > 0
    I_li_c = reshape(I_c_ges(lambda(i),:), 6, 6);

    % [Featherstone2008] T6.2, Line 7 (s steht für Stern)
    % [Featherstone2008] Gl. (2.25) (S. 22, Formelzeichen A,B,r,E s.o.)
    X_li_i_s = [R_i_li(1:3,1:3)', -R_i_li(1:3,1:3)'*skew(r_i_i_li); zeros(3,3), R_i_li(1:3,1:3)'];
    % Transformiere die Trägheit des aktuellen Körpers i in die
    % Koordinaten des Vorgängers und füge zu dessen Komposit-Trägheit
    % hinzu
    I_li_c_new = I_li_c + X_li_i_s*I_i_c*X_i_li;
    I_c_ges(lambda(i),:) = I_li_c_new(:);
  end

  % Gelenk-Transformationsmatrix
  % Siehe [Featherstone2008] Gl. (3.33), Example 3.2
  % S_i = [0; 0; 1; 0; 0; 0];

  % [Featherstone2008] T6.2, Line 9
  % F = I_i_c*S_i;
  if sigma(i) == 0
    F = I_i_c(:,3); % ausnutzen, dass S_i vorgegeben ist
  else
    F = I_i_c(:,6); % TODO: Ergebnis noch falsch.
  end
  % [Featherstone2008] T6.2, Line 10
  % H(i,i) = S_i' * F; % Diagonalelement der Massenmatrix
  H(i,i) = F(3); % ausnutzen, dass S_i vorgegeben ist

  j = uint8(i); % [Featherstone2008] T6.2, Line 11
  for tmp = 1:6 % Dummy-Schleife für Kompilierbarkeit. Abbruchkriterium s.u.
    if lambda(j) == 0 % [Featherstone2008] T6.2, Line 12, Abbruchbedingung
      break;
    end
    % Alt. 1. Berechnung mit kumulativen Transformationsmatrizen (intuitiver):
    % R_0_j = T_c_mdh(1:3,1:3,j+1);
    % R_0_lj = T_c_mdh(1:3,1:3,lambda(j)+1);
    % R_j_lj = R_0_j' * R_0_lj;
    % r_0_j_lj = -T_c_mdh(1:3,4,j+1) + T_c_mdh(1:3,4,lambda(j)+1);
    % r_j_j_lj_test = R_0_j' * r_0_j_lj;

    % Alt. 2. Berechnung mit direkten Transformationsmatrizen (schneller):
    R_j_lj = T_mdh(1:3,1:3,j)';
    r_j_j_lj = -R_j_lj*T_mdh(1:3,4,j);

    % [Featherstone2008], S. 22
    %   A = j
    %   B = lj
    %   r = r_A_A_B = r_j_j_lj
    %   E = R_B_A = R_lj_j = R_j_lj'
    %  Gl. (2.25): X_B_A_s = X_lj_j_s
    X_lj_j_s = [R_j_lj(1:3,1:3)', -R_j_lj(1:3,1:3)'*skew(r_j_j_lj); zeros(3,3), R_j_lj(1:3,1:3)'];

    % [Featherstone2008] T6.2, Line 13
    F = X_lj_j_s*F; % F war vorher in Koordinaten j, hiernach in lambda(j) (des Vorgängers)
    % [Featherstone2008] T6.2, Line 14
    j = (lambda(j)); % Kette weiter nach hinten durchgehen (bis zur Basis)

    if sigma(i) == 0
      S_j = [0; 0; 1; 0; 0; 0];% Siehe [Featherstone2008] Gl. (3.33), Example 3.2
    else
      S_j = [0; 0; 0; 0; 0; 1]; % TODO: Ergebnis noch falsch.
    end
    % [Featherstone2008] T6.2, Line 15
    H(i,j) = F' * S_j;
    % [Featherstone2008] T6.2, Line 16
    H(j,i) = H(i,j); % Symmetrie ausnutzen
  end
end

