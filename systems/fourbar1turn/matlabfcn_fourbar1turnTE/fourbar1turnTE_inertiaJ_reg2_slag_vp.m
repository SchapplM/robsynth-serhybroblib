% Calculate inertial parameters regressor of joint inertia matrix for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% MM_reg [((2+1)*2/2)x(2*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbar1turnTE_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t75 = pkin(4) ^ 2;
t74 = pkin(3) ^ 2;
t73 = -2 * pkin(2);
t43 = pkin(2) ^ 2;
t44 = pkin(1) ^ 2;
t36 = cos(qJ(2));
t68 = pkin(2) * t36;
t60 = -0.2e1 * pkin(1) * t68 + t44;
t29 = t43 + t60;
t59 = t74 - t75;
t25 = t29 + t59;
t31 = pkin(1) * t36 - pkin(2);
t35 = sin(qJ(2));
t71 = (-pkin(3) - pkin(4));
t23 = ((pkin(2) - t71) * (pkin(2) + t71)) + t60;
t70 = (-pkin(3) + pkin(4));
t24 = ((pkin(2) - t70) * (pkin(2) + t70)) + t60;
t45 = sqrt(-t23 * t24);
t62 = t35 * t45;
t14 = -pkin(1) * t62 - t31 * t25;
t72 = -t14 / 0.2e1;
t69 = pkin(2) * t35;
t67 = t35 * pkin(1);
t57 = pkin(2) * t67;
t66 = 0.1e1 / t45 * (-t23 - t24) * t57;
t27 = 0.1e1 / t29;
t38 = 0.1e1 / pkin(4);
t65 = t27 * t38;
t41 = 0.1e1 / pkin(3);
t64 = t27 * t41;
t28 = 0.1e1 / t29 ^ 2;
t63 = t28 / t75;
t61 = t36 * t45;
t58 = 0.2e1 * t68;
t56 = pkin(2) * t64;
t55 = t28 * t69;
t26 = t29 - t59;
t30 = pkin(1) - t68;
t15 = -pkin(2) * t62 + t30 * t26;
t54 = t15 * t65;
t21 = t26 * t69;
t16 = t30 * t45 + t21;
t53 = t16 * t65;
t51 = t63 / 0.4e1;
t46 = t14 ^ 2;
t10 = 0.1e1 / t46;
t22 = t25 * t67;
t17 = -t31 * t45 + t22;
t13 = t17 ^ 2;
t33 = t35 ^ 2;
t2 = 0.1e1 + (((0.2e1 * t44 * t33 * pkin(2) - t31 * t66) * t27 + ((t36 * t25 + t62) * t27 - 0.2e1 * t17 * t55) * pkin(1)) / t14 - (t22 * t27 + (-t27 * t61 + ((t31 * t73 - t66) * t27 + t14 * t28 * t73) * t35) * pkin(1)) * t17 * t10) * pkin(3) * t29 * t41 / (t13 * t10 + 0.1e1);
t50 = t2 * t56;
t47 = t15 ^ 2;
t34 = t36 ^ 2;
t12 = t16 ^ 2;
t11 = 0.1e1 / t47;
t5 = (t35 * t17 / 0.2e1 + t36 * t72) * t64;
t3 = (t14 * t35 + t17 * t36) * t64 / 0.2e1;
t1 = 0.2e1 * (-((t30 * t66 + (t36 * t26 + t62) * pkin(2)) * t27 / 0.2e1 + (t43 * t33 * t27 - t16 * t55) * pkin(1)) / t15 - (-(t21 + (-t35 * t66 - t61) * pkin(2)) * t27 / 0.2e1 + (t15 * t28 - t27 * t30) * t57) * t16 * t11) * pkin(4) * t29 * t38 / (t12 * t11 + 0.1e1);
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t33, 0.2e1 * t35 * t36, 0, t34, 0, 0, 0, 0, 0, 0, t3 ^ 2, -0.2e1 * t3 * t5, 0, t5 ^ 2, 0, 0, t5 * t58, t3 * t58, 0, t34 * t43, t12 * t51, -t16 * t15 * t63 / 0.2e1, 0, t47 * t51, 0, 0, -pkin(1) * t54, -pkin(1) * t53, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t36, 0, 0, 0, 0, 0, 0, 0, -t3 * t2, 0, t5 * t2, 0, 0, 0, (-t17 * t5 / 0.2e1 + t3 * t72) * t56, 0, 0, 0, t1 * t53 / 0.2e1, 0, -t1 * t54 / 0.2e1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2, -t14 * t50, t17 * t50, 0, (t13 / 0.4e1 + t46 / 0.4e1) * t43 / t74 * t28, 0, 0, 0, 0, 0, t1 ^ 2, 0, 0, 0, 0;];
MM_reg = t4;
