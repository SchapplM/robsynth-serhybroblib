% Calculate homogenous joint transformation matrices for
% palh4m1DE2
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
% Datum: 2020-04-11 22:54
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh4m1DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 22:35:49
% EndTime: 2020-04-11 22:35:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (804->56), mult. (874->58), div. (72->9), fcn. (254->16), ass. (0->196)
unknown=NaN(36,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = cos(qJ(5));
t4 = sin(qJ(5));
t5 = t4 * pkin(1);
t7 = 0.2e1 * pkin(2) * t5;
t8 = pkin(1) ^ 2;
t17 = (-t7 + t8 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3))) * (-t7 + t8 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3)));
t18 = sqrt(-t17);
t21 = t5 - pkin(2);
t22 = pkin(2) ^ 2;
t23 = qJ(2) + pkin(6);
t24 = t23 ^ 2;
t25 = pkin(3) ^ 2;
t26 = -t7 + t8 + t22 + t24 - t25;
t28 = -pkin(1) * t18 * t3 - t26 * t21;
t29 = 0.1e1 / t23;
t31 = -t7 + t8 + t22;
t32 = 0.1e1 / t31;
t36 = t26 * t3 * pkin(1) - t18 * t21;
t37 = t36 ^ 2;
t38 = 0.1e1 / t24;
t40 = t31 ^ 2;
t41 = 0.1e1 / t40;
t43 = t28 ^ 2;
t47 = sqrt(t41 * t38 * t37 + t41 * t38 * t43);
t49 = 0.1e1 / t47 * t32;
t50 = t49 * t29 * t28;
t52 = t49 * t29 * t36;
t53 = 0.1e1 / pkin(3);
t55 = 0.1e1 / t25;
t58 = t7 - t8 - t22 + t24 + t25;
t59 = t58 ^ 2;
t63 = sqrt(-t17 * t55 * t38 + t55 * t38 * t59);
t64 = 0.1e1 / t63;
t66 = t64 * t18 * t53 * t29;
t69 = t64 * t53 * t29 * t58;
t70 = cos(qJ(3));
t71 = sin(qJ(3));
t72 = cos(qJ(4));
t73 = sin(qJ(4));
t74 = t3 * pkin(2);
t77 = -t4 * pkin(2) + pkin(1);
t78 = -t7 + t8 + t22 - t24 + t25;
t80 = -t18 * t74 + t78 * t77;
t84 = t18 * t77 + t78 * t74;
t85 = t84 ^ 2;
t88 = t80 ^ 2;
t92 = sqrt(t41 * t55 * t85 + t41 * t55 * t88);
t94 = 0.1e1 / t92 * t32;
t95 = t94 * t53 * t80;
t97 = t94 * t53 * t84;
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
unknown(5,1) = t50;
unknown(5,2) = -t52;
unknown(5,3) = 0.0e0;
unknown(5,4) = -pkin(9);
unknown(6,1) = 0.0e0;
unknown(6,2) = 0.0e0;
unknown(6,3) = -0.1e1;
unknown(6,4) = 0.0e0;
unknown(7,1) = t52;
unknown(7,2) = t50;
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
unknown(13,1) = -t66;
unknown(13,2) = -t69;
unknown(13,3) = 0.0e0;
unknown(13,4) = 0.0e0;
unknown(14,1) = 0.0e0;
unknown(14,2) = 0.0e0;
unknown(14,3) = -0.1e1;
unknown(14,4) = 0.0e0;
unknown(15,1) = t69;
unknown(15,2) = -t66;
unknown(15,3) = 0.0e0;
unknown(15,4) = 0.0e0;
unknown(16,1) = 0.0e0;
unknown(16,2) = 0.0e0;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.1e1;
unknown(17,1) = t70;
unknown(17,2) = -t71;
unknown(17,3) = 0.0e0;
unknown(17,4) = pkin(4);
unknown(18,1) = t71;
unknown(18,2) = t70;
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
unknown(22,1) = -t72;
unknown(22,2) = t73;
unknown(22,3) = 0.0e0;
unknown(22,4) = 0.0e0;
unknown(23,1) = -t73;
unknown(23,2) = -t72;
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
unknown(29,1) = t95;
unknown(29,2) = t97;
unknown(29,3) = 0.0e0;
unknown(29,4) = -pkin(3);
unknown(30,1) = -t97;
unknown(30,2) = t95;
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
