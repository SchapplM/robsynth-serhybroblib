% Explicit kinematic constraints of
% palh1m1TE
% Numeric calculation from joint transformation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% jv [15x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% Berechnungen Schappler; 30.11.2018

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = picker2Dm1DE2_kinconstr_expl_mdh_num_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_kinconstr_expl_mdh_num_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_kinconstr_expl_mdh_num_varpar: pkin has to be [9x1] (double)');

%% Berechnung
% Kinematikparameter und Struktureigenschaften
[beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = picker2Dm1DE2_pkin2mdhparam(pkin);
[~, sigma_mdh] = picker2Dm1DE2_structural_kinematic_parameters();

% Kinematik: Einzel-Gelenktransformationen
T = picker2Dm1DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin);

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
