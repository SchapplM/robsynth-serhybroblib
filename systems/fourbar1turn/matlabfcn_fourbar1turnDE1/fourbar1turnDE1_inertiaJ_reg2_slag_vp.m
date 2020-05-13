% Calculate inertial parameters regressor of joint inertia matrix for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbar1turnDE1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_inertiaJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t49 = pkin(2) ^ 2;
t50 = pkin(1) ^ 2;
t42 = cos(qJ(2));
t70 = pkin(2) * t42;
t63 = -0.2e1 * t70;
t65 = pkin(1) * t63 + t50;
t35 = t49 + t65;
t76 = pkin(3) ^ 2;
t77 = pkin(4) ^ 2;
t64 = t76 - t77;
t31 = t35 + t64;
t37 = pkin(1) * t42 - pkin(2);
t41 = sin(qJ(2));
t73 = -pkin(3) - pkin(4);
t29 = (pkin(2) - t73) * (pkin(2) + t73) + t65;
t72 = -pkin(3) + pkin(4);
t30 = (pkin(2) - t72) * (pkin(2) + t72) + t65;
t51 = sqrt(-t29 * t30);
t67 = t41 * t51;
t20 = -pkin(1) * t67 - t31 * t37;
t78 = t20 ^ 2;
t75 = -0.2e1 * pkin(2);
t28 = pkin(1) * t41 * t31;
t23 = -t37 * t51 + t28;
t19 = t23 ^ 2;
t33 = 0.1e1 / t35;
t34 = 0.1e1 / t35 ^ 2;
t47 = 0.1e1 / pkin(3);
t59 = t33 * t47 * ((t19 + t78) / t76 * t34) ^ (-0.1e1 / 0.2e1);
t5 = (-t20 * t42 + t23 * t41) * t59;
t74 = 0.2e1 * t5;
t71 = pkin(2) * t41;
t62 = pkin(1) * t71;
t69 = 0.1e1 / t51 * (-t29 - t30) * t62;
t68 = t34 / t77;
t66 = t42 * t51;
t61 = t34 * t71;
t32 = t35 - t64;
t36 = pkin(1) - t70;
t21 = -pkin(2) * t67 + t32 * t36;
t16 = t21 ^ 2;
t27 = t32 * t71;
t22 = t36 * t51 + t27;
t18 = t22 ^ 2;
t10 = (t16 + t18) * t68;
t44 = 0.1e1 / pkin(4);
t60 = t33 * t44 * t10 ^ (-0.1e1 / 0.2e1);
t58 = 0.1e1 / t10 * t68;
t57 = pkin(2) * t59;
t17 = 0.1e1 / t21 ^ 2;
t39 = t41 ^ 2;
t1 = 0.2e1 * (-((t36 * t69 + (t32 * t42 + t67) * pkin(2)) * t33 / 0.2e1 + (t33 * t39 * t49 - t22 * t61) * pkin(1)) / t21 - (-(t27 + (-t41 * t69 - t66) * pkin(2)) * t33 / 0.2e1 + (t21 * t34 - t33 * t36) * t62) * t22 * t17) * pkin(4) / (t17 * t18 + 0.1e1) * t35 * t44;
t56 = t1 * t60;
t53 = -0.2e1 * pkin(1) * t60;
t15 = 0.1e1 / t78;
t2 = 0.1e1 + (((0.2e1 * t50 * t39 * pkin(2) - t37 * t69) * t33 + ((t31 * t42 + t67) * t33 - 0.2e1 * t23 * t61) * pkin(1)) / t20 - (t28 * t33 + (-t33 * t66 + ((t37 * t75 - t69) * t33 + t20 * t34 * t75) * t41) * pkin(1)) * t23 * t15) * pkin(3) / (t15 * t19 + 0.1e1) * t35 * t47;
t52 = t2 * t57;
t40 = t42 ^ 2;
t3 = (-t20 * t41 - t23 * t42) * t59;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39, 0.2e1 * t41 * t42, 0, t40, 0, 0, 0, 0, 0, 0, t3 ^ 2, t3 * t74, 0, t5 ^ 2, 0, 0, t70 * t74, t3 * t63, 0, t40 * t49, t18 * t58, -0.2e1 * t22 * t21 * t58, 0, t16 * t58, 0, 0, t21 * t53, t22 * t53, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t42, 0, 0, 0, 0, 0, 0, 0, t3 * t2, 0, t5 * t2, 0, 0, 0, (t20 * t3 - t23 * t5) * t57, 0, 0, 0, t22 * t56, 0, -t21 * t56, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2, -0.2e1 * t20 * t52, 0.2e1 * t23 * t52, 0, t49, 0, 0, 0, 0, 0, t1 ^ 2, 0, 0, 0, 0;];
MM_reg = t4;
