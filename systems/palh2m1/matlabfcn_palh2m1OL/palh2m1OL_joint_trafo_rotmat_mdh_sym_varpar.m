% Calculate homogenous joint transformation matrices for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh2m1OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:58:22
% EndTime: 2020-05-02 23:58:22
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t39 = cos(qJ(1));
t38 = cos(qJ(2));
t37 = cos(qJ(3));
t36 = cos(qJ(4));
t35 = cos(qJ(5));
t34 = sin(qJ(1));
t33 = sin(qJ(2));
t32 = sin(qJ(3));
t31 = sin(qJ(4));
t30 = sin(qJ(5));
t1 = [t39, -t34, 0, 0; t34, t39, 0, 0; 0, 0, 1, pkin(5); 0, 0, 0, 1; t38, -t33, 0, pkin(1); 0, 0, 1, 0; -t33, -t38, 0, 0; 0, 0, 0, 1; t37, -t32, 0, pkin(2); t32, t37, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t36, -t31, 0, pkin(3); t31, t36, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t30, 0, pkin(4); 0, 0, -1, -pkin(6); t30, t35, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
