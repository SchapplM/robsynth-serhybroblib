% Calculate homogenous joint transformation matrices for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = fourbarprisDE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:43:53
% EndTime: 2020-05-07 09:43:53
% DurationCPUTime: 0.05s
% Computational Cost: add. (110->15), mult. (72->25), div. (24->3), fcn. (12->2), ass. (0->22)
t63 = qJ(1) + pkin(3);
t60 = -pkin(2) + t63;
t61 = -pkin(2) - t63;
t43 = sqrt(-(pkin(1) + t61) * (pkin(1) + t60) * (pkin(1) - t60) * (pkin(1) - t61));
t68 = -t43 / 0.2e1;
t67 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t54 = 0.1e1 / pkin(1);
t66 = t54 / 0.2e1;
t48 = 0.1e1 / t63;
t52 = 0.1e1 / pkin(2);
t64 = t48 * t52;
t59 = t54 * t68;
t58 = t64 / 0.2e1;
t57 = t48 * t66;
t56 = t52 * t66;
t51 = (pkin(2) ^ 2);
t55 = t51 + t67;
t53 = pkin(1) ^ 2;
t46 = (t53 + t55) * t56;
t45 = (-t53 + t55) * t57;
t44 = (t51 - t53 - t67) * t58;
t1 = [t45, t43 * t57, 0, pkin(1); t48 * t59, t45, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, qJ(1); -1, 0, 0, 0; 0, -1, 0, 0; 0, 0, 0, 1; t46, t43 * t56, 0, 0; t52 * t59, t46, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t44, t43 * t58, 0, pkin(2); t64 * t68, t44, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, -1, 0, 0; 0, 0, -1, 0; 1, 0, 0, pkin(3); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
