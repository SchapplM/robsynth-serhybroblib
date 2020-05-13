% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh3m2DE2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energykin_fixb_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:42
% EndTime: 2020-05-07 04:22:44
% DurationCPUTime: 1.06s
% Computational Cost: add. (548->175), mult. (1369->350), div. (0->0), fcn. (978->20), ass. (0->155)
t83 = sin(pkin(16));
t84 = cos(pkin(16));
t94 = sin(pkin(15));
t99 = cos(pkin(15));
t34 = t83 * t99 + t84 * t94;
t75 = pkin(17) + pkin(18);
t65 = sin(t75);
t186 = t65 * t34;
t103 = qJD(1) ^ 2;
t97 = cos(qJ(3));
t176 = t97 * pkin(4);
t57 = -pkin(1) + t176;
t98 = cos(qJ(2));
t173 = t57 * t98;
t92 = sin(qJ(3));
t93 = sin(qJ(2));
t166 = t93 * t92;
t55 = pkin(4) * t166;
t47 = t55 + pkin(12);
t25 = -t47 + t173;
t154 = t103 * t25;
t80 = t99 ^ 2;
t67 = -t80 / 0.2e1;
t142 = t67 + 0.1e1 / 0.4e1;
t169 = t83 * t84;
t63 = t80 - 0.1e1 / 0.2e1;
t123 = t63 * pkin(10);
t160 = t99 * t94;
t32 = pkin(8) * t160 + t123;
t140 = pkin(10) * t160;
t33 = t63 * pkin(8) - t140;
t73 = t84 ^ 2;
t119 = -t32 * t169 + t33 * t73;
t66 = cos(t75);
t62 = t66 ^ 2;
t152 = t103 * t62;
t44 = t94 * pkin(8) + t99 * pkin(10);
t45 = t99 * pkin(8) - t94 * pkin(10);
t22 = -t83 * t44 + t45 * t84;
t35 = -t83 * t94 + t84 * t99;
t52 = -t160 / 0.2e1;
t131 = pkin(8) * t52 + t142 * pkin(10);
t136 = t33 * t169;
t4 = t32 * t73 + t131 + t136;
t179 = t80 / 0.2e1;
t134 = t179 - 0.1e1 / 0.4e1;
t121 = t134 * pkin(8);
t5 = t140 / 0.2e1 - t121 + t119;
t185 = ((-t25 * t35 + 0.4e1 * t4 * t65) * t103 - qJD(1) * qJD(4) * t22) * t66 - 0.4e1 * t5 * t152 + t154 * t186 + (-t25 * qJD(4) + (-t99 * t45 + 0.2e1 * t119) * qJD(1)) * qJD(1);
t117 = t92 * t34 - t35 * t97;
t158 = pkin(4) * qJD(3);
t17 = t34 * t97 + t92 * t35;
t7 = -pkin(1) * t34 + t17 * pkin(4);
t8 = pkin(1) * t35 + t117 * pkin(4);
t184 = qJD(4) * ((t7 * t98 - t93 * t8) * qJD(2) + (-t117 * t93 + t17 * t98) * t158) * t66;
t91 = sin(qJ(4));
t144 = qJD(1) * t91;
t96 = cos(qJ(4));
t143 = qJD(1) * t96;
t183 = qJD(4) * t65;
t180 = -0.4e1 * t103;
t90 = cos(pkin(17));
t178 = pkin(3) * t90;
t88 = sin(pkin(17));
t177 = t88 * pkin(3);
t175 = t98 * pkin(1);
t76 = qJD(2) ^ 2 / 0.2e1;
t174 = t103 / 0.2e1;
t172 = t65 * t35;
t171 = t65 * t66;
t170 = t66 * t34;
t105 = pkin(10) ^ 2;
t106 = pkin(8) ^ 2;
t72 = -t105 + t106;
t53 = t72 * t80;
t87 = sin(pkin(18));
t168 = t87 * pkin(12);
t89 = cos(pkin(18));
t167 = t87 * t89;
t165 = t93 * t98;
t164 = t94 * t89;
t95 = sin(pkin(14));
t163 = t94 * t95;
t162 = t98 * t92;
t161 = t99 * t87;
t126 = pkin(8) * t140;
t159 = -0.2e1 * t126 + t53;
t157 = pkin(6) * t103;
t12 = t170 + t172;
t155 = t103 * t12 ^ 2;
t56 = pkin(12) + t175;
t153 = t103 * t56;
t82 = qJ(2) + qJ(3);
t69 = cos(t82);
t151 = t103 * t69;
t150 = t103 * t98;
t107 = pkin(4) ^ 2;
t78 = t97 ^ 2;
t149 = t78 * t107;
t100 = cos(pkin(14));
t148 = t95 * t100;
t147 = t99 * t100;
t109 = pkin(1) ^ 2;
t146 = -0.2e1 * pkin(1) * t176 + t109;
t74 = qJD(2) + qJD(3);
t145 = qJD(1) * t74;
t141 = qJD(1) * qJD(2);
t139 = t99 * t177;
t64 = t106 / 0.2e1 - t105 / 0.2e1;
t24 = pkin(8) * t123 + t64 * t160;
t137 = t24 * t169;
t77 = t89 ^ 2;
t135 = t77 * t160;
t133 = pkin(1) * qJD(2) * t74;
t130 = t172 / 0.4e1;
t51 = t160 / 0.4e1;
t129 = t155 / 0.2e1;
t127 = t93 * t141;
t125 = -t106 / 0.4e1 + t105 / 0.4e1;
t39 = t94 * t100 - t99 * t95;
t42 = t147 + t163;
t124 = -t93 * t39 + t42 * t98;
t122 = t97 * t133;
t120 = pkin(1) * t127;
t13 = -t35 * t66 + t186;
t37 = t161 + t164;
t116 = t93 * t99 + t94 * t98;
t115 = t93 * t94 - t98 * t99;
t113 = -(t97 * pkin(1) - pkin(4)) * qJD(2) * t158 + (t107 + t146) * t76 + t107 * qJD(3) ^ 2 / 0.2e1;
t40 = t97 * t93 + t162;
t16 = (pkin(4) * t162 + t57 * t93) * qJD(2) + t40 * t158;
t112 = -t63 * t73 + t160 * t169 - 0.1e1 / 0.4e1;
t38 = -t94 * t87 + t99 * t89;
t111 = -t37 * t177 + t38 * t178 - pkin(12);
t104 = pkin(12) ^ 2;
t79 = t98 ^ 2;
t110 = -0.2e1 * t47 * t173 + t104 + t107 - t149 + (-t107 + t146 + 0.2e1 * t149) * t79 + 0.2e1 * pkin(12) * t55;
t81 = t100 ^ 2;
t71 = t109 * t76;
t68 = sin(t82);
t61 = t74 ^ 2 / 0.2e1;
t46 = t92 * t133;
t41 = -t98 * t97 + t166;
t36 = t56 ^ 2 * t174 + t71;
t23 = -t64 + t159;
t21 = t44 * t84 + t83 * t45;
t20 = t73 * t160 + t63 * t169 + t52;
t19 = t98 * t39 + t93 * t42;
t15 = t116 * t100 + t115 * t95;
t14 = t115 * t100 - t116 * t95;
t10 = -t13 * qJD(1) + qJD(4);
t9 = t111 - t175;
t3 = (t170 / 0.4e1 + t130) * qJD(4) + ((t67 - t112) * t171 + (t62 - 0.1e1 / 0.2e1) * t20) * qJD(1);
t2 = (t93 * t7 + t8 * t98) * qJD(2) - (-t117 * t98 - t93 * t17) * t158;
t1 = [0, 0, 0, 0, 0, t174, 0, 0, 0, 0, t93 ^ 2 * t174, t93 * t150, t127, t79 * t174, t98 * t141, t76, pkin(12) * t150, -t103 * pkin(12) * t93, 0, t104 * t174, t40 ^ 2 * t174, 0.2e1 * ((t78 - 0.1e1 / 0.2e1) * t165 + (t79 - 0.1e1 / 0.2e1) * t97 * t92) * t103, -t40 * t145, t41 ^ 2 * t174, t41 * t145, t61, t41 * t153 - t122, t40 * t153 + t46, -t120, t36, t129, (t20 * t62 - (t179 + t112) * t171 + t73 * t52 - t134 * t169 + t51) * t180, 0, t13 ^ 2 * t174, 0, 0, -t13 * t154, t12 * t154, qJD(1) * t16, t110 * t174 + t113, t96 ^ 2 * t129, -t96 * t91 * t155, 0.4e1 * t3 * t143, t91 ^ 2 * t129, -0.4e1 * t3 * t144, t10 ^ 2 / 0.2e1, t91 * t184 + (-t21 * t143 - t2 * t91) * t183 + t16 * t144 - t185 * t96, t96 * t184 + (t21 * t144 - t2 * t96) * t183 + t16 * t143 + t185 * t91, (t4 * t62 + (t5 * t65 + t25 * t34 / 0.4e1) * t66 + t25 * t130 + t131 * t73 - t136 / 0.2e1 + pkin(8) * t51 + (-0.1e1 / 0.4e1 + t80 / 0.4e1) * pkin(10)) * t180, 0.2e1 * (t23 * t73 - t64 * t80 - t125 + t126 - 0.2e1 * t137) * t152 + ((t24 * t73 + t23 * t169 / 0.2e1 + t125 * t160 - pkin(10) * t121) * t65 - t25 * t22 / 0.4e1) * t66 * t180 - t21 * t65 * t154 + ((0.4e1 * t126 + t72 - 0.2e1 * t53) * t73 + 0.4e1 * t137 + t105 + t110 + t159) * t174 + t113, t15 ^ 2 * t174, 0.4e1 * ((-t63 * t148 + t81 * t160 + t52) * t79 + (t147 * t163 + t63 * t81 + t142) * t165 + t81 * t52 + t134 * t148 + t51) * t103, t19 * t141, t14 ^ 2 * t174, t124 * t141, t76, -t124 * t157, t19 * t157, 0, pkin(6) ^ 2 * t174, t37 ^ 2 * t174, -0.2e1 * t103 * (t63 * t167 + t135 + t52), 0, t38 ^ 2 * t174, 0, 0, -t38 * t153, -t37 * t153, -t120, t36, t68 ^ 2 * t174, t68 * t151, -t68 * t145, t69 ^ 2 * t174, -t69 * t145, t61, t151 * t9 - t122, -t9 * t103 * t68 + t46, -t120, t71 + (t79 * t109 - 0.2e1 * t111 * t175 - 0.4e1 * (t135 * t177 + (t99 * pkin(12) / 0.2e1 + t63 * t87 * t177) * t89 - t94 * (t139 + t168) / 0.2e1) * t178 + 0.2e1 * pkin(3) * (pkin(3) * t161 + t88 * pkin(12)) * t164 + 0.2e1 * t139 * t168 + t104 + (0.4e1 * (-t160 * t167 + t63 * t77 + t142) * t90 ^ 2 + (-0.2e1 * t80 + 0.1e1) * t77 + t80) * pkin(3) ^ 2) * t174;];
T_reg = t1;
