% Explicit kinematic constraints of
% palh1m1TE
% Numeric calculation from joint transformation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% jv [15x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% Berechnungen Schappler; 30.11.2018

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:01
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = mg10hlDE2_kinconstr_expl_mdh_num_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE2_kinconstr_expl_mdh_num_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE2_kinconstr_expl_mdh_num_varpar: pkin has to be [17x1] (double)');

%% Berechnung
% Kinematikparameter und Struktureigenschaften
[beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = mg10hlDE2_pkin2mdhparam(pkin);
[~, sigma_mdh] = mg10hlDE2_structural_kinematic_parameters();

% Kinematik: Einzel-Gelenktransformationen
T = mg10hlDE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin);

jv = zeros(15, 1);
for i = 1:15
  % MDH-Parameter
  sigma = sigma_mdh(i);
  beta = beta_mdh(i);
  b = b_mdh(i);
  alpha = alpha_mdh(i);
  a = a_mdh(i);
  theta = theta_mdh(i);
  d = d_mdh(i);
  q_offset = qoffset_mdh(i);
  % Transformationsschritte außer theta und d
  T_1 = trotz(beta) * transl([0;0;b]) * trotx(alpha) * transl([a;0;0]);
  % MDH-Transformation so umformen, dass nur noch die Transformation
  % in die Richtung des Gelenks (theta und d) übrig bleibt
  if sigma == 0 % Drehgelenk
    % Gelenkvariable theta aus Rotationsmatrix berechnen
    R_ztheta = T_1(1:3,1:3)' * T(1:3,1:3,i);
    theta = atan2(R_ztheta(2,1), R_ztheta(1,1));
    jv(i) = theta - q_offset;
  elseif sigma == 1 % Schubgelenk
    % Gelenkvariable d aus Translations-Teil der Transformationsmatrix
    T_Tzd = invtr(T_1 * trotz(theta)) * T(:,:,i);
    d = T_Tzd(3,4);
    jv(i) = d - q_offset;
  end
end
