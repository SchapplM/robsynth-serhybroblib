% Calculate homogenous joint transformation matrices for
% fourbarprisTE
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
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = fourbarprisTE_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:05:46
% EndTime: 2020-06-27 17:05:47
% DurationCPUTime: 0.05s
% Computational Cost: add. (110->15), mult. (72->25), div. (24->3), fcn. (12->2), ass. (0->22)
t50 = qJ(1) + pkin(3);
t47 = -pkin(2) + t50;
t48 = -pkin(2) - t50;
t30 = sqrt(-(pkin(1) + t48) * (pkin(1) + t47) * (pkin(1) - t47) * (pkin(1) - t48));
t55 = -t30 / 0.2e1;
t54 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t41 = 0.1e1 / pkin(1);
t53 = t41 / 0.2e1;
t35 = 0.1e1 / t50;
t39 = 0.1e1 / pkin(2);
t51 = t35 * t39;
t46 = t41 * t55;
t45 = t51 / 0.2e1;
t44 = t35 * t53;
t43 = t39 * t53;
t38 = (pkin(2) ^ 2);
t42 = t38 + t54;
t40 = pkin(1) ^ 2;
t33 = (t40 + t42) * t43;
t32 = (-t40 + t42) * t44;
t31 = (t38 - t40 - t54) * t45;
t1 = [t32, t30 * t44, 0, pkin(1); t35 * t46, t32, 0, 0; 0, 0, 1, 0; 0, 0, 1, qJ(1); -1, 0, 0, 0; 0, -1, 0, 0; t33, t30 * t43, 0, 0; t39 * t46, t33, 0, 0; 0, 0, 1, 0; t31, t30 * t45, 0, pkin(2); t51 * t55, t31, 0, 0; 0, 0, 1, 0; 0, -1, 0, 0; 0, 0, -1, 0; 1, 0, 0, pkin(3);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
