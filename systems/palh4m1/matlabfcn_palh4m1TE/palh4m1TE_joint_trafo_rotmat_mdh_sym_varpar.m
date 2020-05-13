% Calculate homogenous joint transformation matrices for
% palh4m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% T_mdh [4x4x9]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh4m1TE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 21:35:24
% EndTime: 2020-04-11 21:35:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (271->53), mult. (274->41), div. (24->3), fcn. (90->10), ass. (0->174)
unknown=NaN(36,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = cos(qJ(5));
t4 = sin(qJ(5));
t5 = pkin(1) * t4;
t7 = 0.2e1 * t5 * pkin(2);
t8 = pkin(1) ^ 2;
t18 = sqrt(-(-t7 + t8 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3))) * (-t7 + t8 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3))));
t21 = t5 - pkin(2);
t22 = pkin(2) ^ 2;
t23 = qJ(2) + pkin(6);
t24 = t23 ^ 2;
t25 = pkin(3) ^ 2;
t26 = -t7 + t8 + t22 + t24 - t25;
t29 = 0.1e1 / t23;
t32 = 0.1e1 / (-t7 + t8 + t22);
t34 = (-t3 * t18 * pkin(1) - t21 * t26) * t29 * t32 / 0.2e1;
t41 = (pkin(1) * t3 * t26 - t21 * t18) * t29 * t32 / 0.2e1;
t42 = 0.1e1 / pkin(3);
t45 = t29 * t42 * t18 / 0.2e1;
t49 = (t7 - t8 - t22 + t24 + t25) * t29 * t42 / 0.2e1;
t50 = cos(qJ(3));
t51 = sin(qJ(3));
t52 = cos(qJ(4));
t53 = sin(qJ(4));
t54 = pkin(2) * t3;
t57 = -pkin(2) * t4 + pkin(1);
t58 = -t7 + t8 + t22 - t24 + t25;
t63 = (-t54 * t18 + t57 * t58) * t42 * t32 / 0.2e1;
t69 = (t57 * t18 + t54 * t58) * t42 * t32 / 0.2e1;
unknown(1,1) = t1;
unknown(1,2) = -t2;
unknown(1,3) = 0.0e0;
unknown(1,4) = 0.0e0;
unknown(2,1) = t2;
unknown(2,2) = t1;
unknown(2,3) = 0.0e0;
unknown(2,4) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = 0.1e1;
unknown(3,4) = pkin(7);
unknown(4,1) = 0.0e0;
unknown(4,2) = 0.0e0;
unknown(4,3) = 0.0e0;
unknown(4,4) = 0.1e1;
unknown(5,1) = t34;
unknown(5,2) = -t41;
unknown(5,3) = 0.0e0;
unknown(5,4) = -pkin(9);
unknown(6,1) = 0.0e0;
unknown(6,2) = 0.0e0;
unknown(6,3) = -0.1e1;
unknown(6,4) = 0.0e0;
unknown(7,1) = t41;
unknown(7,2) = t34;
unknown(7,3) = 0.0e0;
unknown(7,4) = 0.0e0;
unknown(8,1) = 0.0e0;
unknown(8,2) = 0.0e0;
unknown(8,3) = 0.0e0;
unknown(8,4) = 0.1e1;
unknown(9,1) = 0.0e0;
unknown(9,2) = 0.0e0;
unknown(9,3) = 0.1e1;
unknown(9,4) = t23;
unknown(10,1) = -0.1e1;
unknown(10,2) = 0.0e0;
unknown(10,3) = 0.0e0;
unknown(10,4) = 0.0e0;
unknown(11,1) = 0.0e0;
unknown(11,2) = -0.1e1;
unknown(11,3) = 0.0e0;
unknown(11,4) = 0.0e0;
unknown(12,1) = 0.0e0;
unknown(12,2) = 0.0e0;
unknown(12,3) = 0.0e0;
unknown(12,4) = 0.1e1;
unknown(13,1) = -t45;
unknown(13,2) = -t49;
unknown(13,3) = 0.0e0;
unknown(13,4) = 0.0e0;
unknown(14,1) = 0.0e0;
unknown(14,2) = 0.0e0;
unknown(14,3) = -0.1e1;
unknown(14,4) = 0.0e0;
unknown(15,1) = t49;
unknown(15,2) = -t45;
unknown(15,3) = 0.0e0;
unknown(15,4) = 0.0e0;
unknown(16,1) = 0.0e0;
unknown(16,2) = 0.0e0;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.1e1;
unknown(17,1) = t50;
unknown(17,2) = -t51;
unknown(17,3) = 0.0e0;
unknown(17,4) = pkin(4);
unknown(18,1) = t51;
unknown(18,2) = t50;
unknown(18,3) = 0.0e0;
unknown(18,4) = 0.0e0;
unknown(19,1) = 0.0e0;
unknown(19,2) = 0.0e0;
unknown(19,3) = 0.1e1;
unknown(19,4) = 0.0e0;
unknown(20,1) = 0.0e0;
unknown(20,2) = 0.0e0;
unknown(20,3) = 0.0e0;
unknown(20,4) = 0.1e1;
unknown(21,1) = 0.0e0;
unknown(21,2) = 0.0e0;
unknown(21,3) = 0.1e1;
unknown(21,4) = pkin(5);
unknown(22,1) = -t52;
unknown(22,2) = t53;
unknown(22,3) = 0.0e0;
unknown(22,4) = 0.0e0;
unknown(23,1) = -t53;
unknown(23,2) = -t52;
unknown(23,3) = 0.0e0;
unknown(23,4) = 0.0e0;
unknown(24,1) = 0.0e0;
unknown(24,2) = 0.0e0;
unknown(24,3) = 0.0e0;
unknown(24,4) = 0.1e1;
unknown(25,1) = -t4;
unknown(25,2) = -t3;
unknown(25,3) = 0.0e0;
unknown(25,4) = pkin(8);
unknown(26,1) = 0.0e0;
unknown(26,2) = 0.0e0;
unknown(26,3) = -0.1e1;
unknown(26,4) = 0.0e0;
unknown(27,1) = t3;
unknown(27,2) = -t4;
unknown(27,3) = 0.0e0;
unknown(27,4) = 0.0e0;
unknown(28,1) = 0.0e0;
unknown(28,2) = 0.0e0;
unknown(28,3) = 0.0e0;
unknown(28,4) = 0.1e1;
unknown(29,1) = t63;
unknown(29,2) = t69;
unknown(29,3) = 0.0e0;
unknown(29,4) = -pkin(3);
unknown(30,1) = -t69;
unknown(30,2) = t63;
unknown(30,3) = 0.0e0;
unknown(30,4) = 0.0e0;
unknown(31,1) = 0.0e0;
unknown(31,2) = 0.0e0;
unknown(31,3) = 0.1e1;
unknown(31,4) = 0.0e0;
unknown(32,1) = 0.0e0;
unknown(32,2) = 0.0e0;
unknown(32,3) = 0.0e0;
unknown(32,4) = 0.1e1;
unknown(33,1) = 0.1e1;
unknown(33,2) = 0.0e0;
unknown(33,3) = 0.0e0;
unknown(33,4) = pkin(1);
unknown(34,1) = 0.0e0;
unknown(34,2) = 0.1e1;
unknown(34,3) = 0.0e0;
unknown(34,4) = 0.0e0;
unknown(35,1) = 0.0e0;
unknown(35,2) = 0.0e0;
unknown(35,3) = 0.1e1;
unknown(35,4) = 0.0e0;
unknown(36,1) = 0.0e0;
unknown(36,2) = 0.0e0;
unknown(36,3) = 0.0e0;
unknown(36,4) = 0.1e1;
T_ges = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,9);             % numerisch
else,                         T_mdh = sym('xx', [4,4,9]); end % symbolisch

for i = 1:9
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
