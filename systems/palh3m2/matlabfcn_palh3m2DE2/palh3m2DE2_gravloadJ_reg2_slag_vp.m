% Calculate inertial parameters regressor of gravitation load for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh3m2DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_gravloadJ_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t100 = qJ(2) + qJ(3);
t122 = pkin(15) + qJ(2);
t94 = pkin(18) + t122;
t76 = pkin(17) + qJ(3) + t94;
t67 = pkin(16) + t76;
t46 = atan2(-sin(t67), cos(t67));
t120 = t46 + t100;
t42 = qJ(4) + t120;
t38 = qJ(1) + t42;
t43 = -qJ(4) + t120;
t39 = qJ(1) - t43;
t161 = cos(t38) + cos(t39);
t101 = qJ(1) + qJ(4);
t160 = -sin(t101) / 0.2e1 - sin(t39) / 0.4e1 - sin(t38) / 0.4e1;
t103 = sin(qJ(3));
t106 = cos(qJ(3));
t107 = cos(qJ(2));
t105 = sin(qJ(1));
t132 = g(2) * t105;
t108 = cos(qJ(1));
t133 = g(1) * t108;
t65 = t132 + t133;
t159 = (-t106 * g(3) + t103 * t65) * t107;
t102 = qJ(1) - qJ(4);
t40 = qJ(1) + t43;
t41 = qJ(1) - t42;
t117 = cos(t102) / 0.2e1 - cos(t41) / 0.4e1 - cos(t40) / 0.4e1;
t118 = sin(t102) / 0.2e1 - sin(t41) / 0.4e1 - sin(t40) / 0.4e1;
t157 = -t118 * g(1) + t117 * g(2);
t156 = t117 * g(1) + t118 * g(2);
t155 = -g(1) / 0.2e1;
t153 = g(1) / 0.2e1;
t152 = -g(2) / 0.2e1;
t151 = g(2) / 0.2e1;
t123 = qJ(2) + atan2(sin(t94), -cos(t94));
t119 = pkin(17) - t123;
t17 = -atan2(-sin(t76), -cos(t76)) + t119;
t14 = -qJ(1) + t17;
t149 = sin(t14) / 0.2e1;
t13 = qJ(1) + t17;
t148 = cos(t13) / 0.2e1;
t96 = qJ(1) + t100;
t44 = t46 + t96;
t145 = -sin(t44) / 0.2e1;
t97 = qJ(1) - t100;
t45 = -t46 + t97;
t144 = cos(t45) / 0.2e1;
t54 = qJ(1) + t123;
t143 = sin(t54) / 0.2e1;
t55 = qJ(1) - t123;
t142 = -cos(t55) / 0.2e1;
t95 = pkin(14) - t122;
t88 = -qJ(1) + t95;
t141 = -sin(t88) / 0.2e1;
t87 = qJ(1) + t95;
t140 = -cos(t87) / 0.2e1;
t82 = sin(t96);
t139 = -t82 / 0.2e1;
t85 = cos(t97);
t138 = t85 / 0.2e1;
t12 = cos(t14);
t136 = g(2) * t149 + t12 * t155;
t135 = g(1) * t149 + t12 * t151;
t134 = g(3) * cos(t100);
t131 = t103 * g(3);
t129 = t107 * g(3);
t80 = cos(t88);
t128 = g(2) * t141 + t80 * t153;
t127 = g(1) * t141 + t80 * t152;
t84 = cos(t96);
t126 = g(2) * t139 + t84 * t155;
t125 = g(1) * t139 + t84 * t151;
t124 = t107 * pkin(1) + pkin(12);
t104 = sin(qJ(2));
t86 = -t106 * pkin(4) + pkin(1);
t121 = pkin(4) * t103 * t104 + t86 * t107 + pkin(12);
t92 = cos(t101);
t116 = t92 * t151 + t161 * g(2) / 0.4e1 + t160 * g(1);
t115 = t92 * t155 - t161 * g(1) / 0.4e1 + t160 * g(2);
t9 = sin(t13);
t114 = g(1) * t148 + t9 * t151;
t113 = g(2) * t148 + t9 * t155;
t77 = sin(t87);
t112 = g(2) * t140 + t77 * t153;
t111 = g(1) * t140 + t77 * t152;
t83 = sin(t97);
t110 = g(2) * t138 + t83 * t155;
t109 = g(1) * t138 + t83 * t151;
t64 = -g(1) * t105 + g(2) * t108;
t63 = t106 * pkin(8) + t103 * pkin(10);
t62 = -t103 * pkin(8) + t106 * pkin(10);
t56 = t64 * t124;
t53 = t104 * t65 - t129;
t51 = cos(t54);
t50 = sin(t55);
t48 = t53 * pkin(1);
t36 = cos(t44);
t35 = sin(t45);
t16 = -t110 + t125 + t134;
t15 = -g(3) * sin(t100) - t109 + t126;
t4 = (pkin(1) * t132 + t86 * t133) * t104 - pkin(1) * t129 + ((-t106 * t132 - t131) * t104 - t159) * pkin(4);
t3 = (t35 / 0.2e1 + t145) * g(2) + (t144 - t36 / 0.2e1) * g(1);
t2 = sin(t17) * g(3) - t114 + t136;
t1 = cos(t17) * g(3) - t113 + t135;
t5 = [0, 0, 0, 0, 0, 0, -t64, t65, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t64, t104 * t64, -t65, -pkin(12) * t64, 0, 0, 0, 0, 0, 0, t110 + t125, t109 + t126, -t65, -t56, 0, 0, 0, 0, 0, 0, (t144 + t36 / 0.2e1) * g(2) + (-t35 / 0.2e1 + t145) * g(1), t3, -t65, -t64 * t121, 0, 0, 0, 0, 0, 0, t116 - t157, t115 - t156, -t3, -((-t104 * t62 - t63 * t107) * cos(t46) + (t104 * t63 - t62 * t107) * sin(t46) + t121) * t64, 0, 0, 0, 0, 0, 0, t112 + t127, t111 + t128, -t65, pkin(6) * t64, 0, 0, 0, 0, 0, 0, (-t51 / 0.2e1 + t142) * g(2) + (t143 + t50 / 0.2e1) * g(1), (t143 - t50 / 0.2e1) * g(2) + (t51 / 0.2e1 + t142) * g(1), -t65, -t56, 0, 0, 0, 0, 0, 0, t113 + t135, t114 + t136, -t65, -t64 * (pkin(3) * cos(t119) + t124); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * t104 + t107 * t65, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, -g(3) * cos(t95) - t112 + t127, -g(3) * sin(t95) - t111 + t128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, t1, t2, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(t159 + t104 * (t65 * t106 + t131)) * pkin(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, (0.2e1 * t134 + (t84 - t85) * g(2) + (-t82 + t83) * g(1)) * pkin(4) / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (cos(t42) / 0.2e1 - cos(t43) / 0.2e1) * g(3) + t116 + t157, (-sin(t43) / 0.2e1 - sin(t42) / 0.2e1) * g(3) + t115 + t156, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t5;
