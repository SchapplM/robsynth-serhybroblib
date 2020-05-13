% Calculate inertial parameters regressor of gravitation load for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnDE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t83 = pkin(4) ^ 2;
t46 = pkin(2) ^ 2;
t47 = pkin(1) ^ 2;
t39 = cos(qJ(2));
t70 = t39 * pkin(1);
t81 = -2 * pkin(2);
t64 = t70 * t81 + t47;
t32 = t46 + t64;
t29 = 0.1e1 / t32 ^ 2;
t82 = t29 / t83;
t28 = 0.1e1 / t32;
t63 = pkin(3) ^ 2 - t83;
t27 = t32 - t63;
t33 = -pkin(2) * t39 + pkin(1);
t37 = sin(qJ(2));
t78 = -pkin(3) - pkin(4);
t24 = (pkin(2) - t78) * (pkin(2) + t78) + t64;
t77 = -pkin(3) + pkin(4);
t25 = (pkin(2) - t77) * (pkin(2) + t77) + t64;
t49 = sqrt(-t24 * t25);
t67 = t37 * t49;
t16 = -pkin(2) * t67 + t27 * t33;
t72 = pkin(2) * t37;
t22 = t27 * t72;
t17 = t33 * t49 + t22;
t65 = t16 ^ 2 + t17 ^ 2;
t10 = t65 * t82;
t8 = t10 ^ (-0.1e1 / 0.2e1);
t80 = 0.2e1 * t37 ^ 2;
t66 = t49 * t39;
t71 = t37 * pkin(1);
t62 = pkin(2) * t71;
t69 = 0.1e1 / t49 * (-t24 - t25) * t62;
t3 = t22 + (-t66 + (0.2e1 * t33 * pkin(1) - t69) * t37) * pkin(2);
t4 = t33 * t69 + t46 * pkin(1) * t80 + (t27 * t39 + t67) * pkin(2);
t79 = 0.2e1 * (-0.2e1 * t28 * t62 * t65 + t16 * t3 + t17 * t4) * t8 / t10 * t82;
t40 = cos(qJ(1));
t76 = g(1) * t40;
t38 = sin(qJ(1));
t75 = g(2) * t38;
t74 = g(3) * t16;
t73 = g(3) * t17;
t45 = 0.1e1 / pkin(3);
t68 = t28 * t45;
t61 = t29 * t72;
t59 = pkin(1) * t8 * t61;
t31 = t75 + t76;
t58 = g(1) * t38 - g(2) * t40;
t57 = -t76 / 0.2e1 - t75 / 0.2e1;
t56 = t58 * t39;
t55 = 0.2e1 * t31;
t42 = 0.1e1 / pkin(4);
t54 = t58 * t8 * t42 * t28;
t53 = -g(3) * t39 + t31 * t37;
t34 = -pkin(2) + t70;
t26 = t32 + t63;
t23 = t26 * t71;
t18 = -t34 * t49 + t23;
t15 = -pkin(1) * t67 - t26 * t34;
t12 = 0.1e1 / t15 ^ 2;
t7 = qJ(2) + atan2(t18 * t68, t15 * t68);
t6 = cos(t7);
t5 = sin(t7);
t1 = 0.1e1 + (((t47 * pkin(2) * t80 - t34 * t69) * t28 + ((t26 * t39 + t67) * t28 - 0.2e1 * t18 * t61) * pkin(1)) / t15 - (t23 * t28 + (-t28 * t66 + ((t34 * t81 - t69) * t28 + t15 * t29 * t81) * t37) * pkin(1)) * t18 * t12) * pkin(3) / (t12 * t18 ^ 2 + 0.1e1) * t32 * t45;
t2 = [0, 0, 0, 0, 0, 0, t58, t31, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t58 * t37, -t31, 0, 0, 0, 0, 0, 0, 0, -t58 * t6, t58 * t5, -t31, pkin(2) * t56, 0, 0, 0, 0, 0, 0, -t16 * t54, -t17 * t54, -t31, t58 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * t37 + t31 * t39, 0, 0, 0, 0, 0, 0, 0, 0, (g(3) * t6 - t31 * t5) * t1, (-g(3) * t5 - t31 * t6) * t1, 0, t53 * pkin(2), 0, 0, 0, 0, 0, 0, ((-t16 * t55 + 0.2e1 * t73) * t59 + ((-g(3) * t4 + t3 * t31) * t8 + (t73 / 0.2e1 + t57 * t16) * t79) * t28) * t42, ((-t17 * t55 - 0.2e1 * t74) * t59 + ((g(3) * t3 + t31 * t4) * t8 + (-t74 / 0.2e1 + t57 * t17) * t79) * t28) * t42, 0, 0;];
taug_reg = t2;
