% Calculate inertial parameters regressor of gravitation load for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh1m2DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_gravloadJ_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t65 = cos(qJ(3));
t67 = cos(qJ(1));
t83 = sin(qJ(1));
t70 = -g(1) * t67 - g(2) * t83;
t91 = t65 * t70;
t57 = sin(pkin(19));
t58 = cos(pkin(19));
t60 = sin(qJ(3));
t40 = t65 * t57 + t60 * t58;
t41 = t60 * t57 - t65 * t58;
t29 = qJ(2) + atan2(t41, t40);
t27 = sin(t29);
t28 = cos(t29);
t4 = g(3) * t27 - t70 * t28;
t90 = g(3) * t28 + t70 * t27;
t61 = sin(qJ(2));
t66 = cos(qJ(2));
t89 = g(3) * t61 - t70 * t66;
t37 = g(3) * t66 + t70 * t61;
t56 = pkin(18) - pkin(22);
t75 = -pkin(21) + t56;
t49 = -qJ(2) - qJ(3) + t75;
t52 = -qJ(2) + t56;
t76 = pkin(21) - atan2(cos(t52), -sin(t52)) - qJ(2);
t24 = -atan2(-sin(t49), cos(t49)) + t76;
t22 = qJ(1) + t24;
t88 = sin(t22) / 0.2e1;
t23 = -qJ(1) + t24;
t87 = -cos(t23) / 0.2e1;
t86 = g(3) * t65;
t19 = sin(t23);
t84 = g(1) * t87 + g(2) * t19 / 0.2e1;
t45 = -g(1) * t83 + g(2) * t67;
t50 = -pkin(20) + t75;
t47 = sin(t50);
t82 = t45 * t47;
t59 = sin(qJ(4));
t81 = t59 * t45;
t64 = cos(qJ(4));
t80 = t64 * t45;
t79 = t66 * t65;
t78 = -g(1) * t19 / 0.2e1 + g(2) * t87;
t77 = t61 * pkin(1) - pkin(15);
t48 = cos(t50);
t74 = -g(3) * t47 + t48 * t70;
t73 = t61 * t60 - t79;
t20 = cos(t22);
t72 = g(1) * t88 - g(2) * t20 / 0.2e1;
t71 = g(1) * t20 / 0.2e1 + g(2) * t88;
t69 = cos(pkin(17));
t68 = cos(pkin(18));
t63 = sin(pkin(17));
t62 = sin(pkin(18));
t51 = pkin(5) * t60 + pkin(1);
t44 = pkin(15) * t45;
t43 = t63 * t62 + t69 * t68;
t42 = t62 * t69 - t68 * t63;
t38 = t77 * t45;
t36 = -t60 * t70 - t86;
t34 = g(3) * t60 - t91;
t32 = pkin(1) * t89;
t17 = t37 * t68 + t62 * t89;
t16 = t61 * t34 + t36 * t66;
t15 = t34 * t66 - t36 * t61;
t14 = -t62 * t37 + t68 * t89;
t13 = t16 * pkin(5);
t8 = (-pkin(5) * t86 - t70 * t51) * t66 + t61 * (-pkin(5) * t91 + g(3) * t51);
t7 = atan2(t41, -t40) + t29;
t6 = cos(t7);
t5 = sin(t7);
t3 = pkin(2) * t4;
t2 = -cos(t24) * g(3) - t72 + t78;
t1 = sin(t24) * g(3) - t71 + t84;
t9 = [0, 0, 0, 0, 0, 0, -t45, -t70, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t45, t66 * t45, t70, -t44, 0, 0, 0, 0, 0, 0, t73 * t45, t45 * (t66 * t60 + t61 * t65), t70, t38, 0, 0, 0, 0, 0, 0, t45 * t48, -t82, t70, t45 * (t73 * pkin(5) + t77), 0, 0, 0, 0, 0, 0, t48 * t80 + t59 * t70, -t48 * t81 + t64 * t70, t82, t45 * (-pkin(5) * t79 + pkin(9) * t48 + pkin(11) * t47 + t61 * t51 - pkin(15)), 0, 0, 0, 0, 0, 0, (-t42 * t66 + t61 * t43) * t45, t45 * (t61 * t42 + t43 * t66), t70, pkin(14) * t45, 0, 0, 0, 0, 0, 0, t45 * cos(t56), -t45 * sin(t56), t70, t38, 0, 0, 0, 0, 0, 0, t45 * t27, t45 * t28, t70, -t44, 0, 0, 0, 0, 0, 0, -t45 * t5, -t45 * t6, t70, t45 * (pkin(2) * t27 - pkin(15)), 0, 0, 0, 0, 0, 0, t71 + t84, t72 + t78, t70, t45 * (-pkin(4) * sin(t76) + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t37, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, t14 * t69 + t63 * t17, -t14 * t63 + t17 * t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, t4, t90, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t5 + t6 * t70, -g(3) * t6 - t5 * t70, 0, t3, 0, 0, 0, 0, 0, 0, t1, t2, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t59 + t80, t74 * t64 - t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t9;
