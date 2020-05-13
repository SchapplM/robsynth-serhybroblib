% Calculate homogenous joint transformation matrices for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh2m2DE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:09
% EndTime: 2020-05-03 01:06:09
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (0->0), div. (0->0), fcn. (24->8), ass. (0->9)
t27 = cos(qJ(1));
t26 = cos(qJ(2));
t25 = cos(qJ(3));
t24 = cos(qJ(4));
t23 = sin(qJ(1));
t22 = sin(qJ(2));
t21 = sin(qJ(3));
t20 = sin(qJ(4));
t1 = [t27, -t23, 0, 0; t23, t27, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t22, 0, pkin(1); 0, 0, -1, 0; t22, t26, 0, 0; 0, 0, 0, 1; t26, t22, 0, pkin(4); -t22, t26, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t21, 0, pkin(2); t21, t25, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, t21, 0, pkin(5); -t21, t25, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t20, 0, pkin(3); 0, 0, 1, 0; -t20, -t24, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
