% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh3m2DE1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_energykin_fixb_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:03:37
% EndTime: 2020-05-07 02:03:38
% DurationCPUTime: 1.00s
% Computational Cost: add. (548->173), mult. (1369->346), div. (0->0), fcn. (978->20), ass. (0->153)
t83 = sin(pkin(16));
t84 = cos(pkin(16));
t94 = sin(pkin(15));
t99 = cos(pkin(15));
t34 = t83 * t99 + t84 * t94;
t75 = pkin(17) + pkin(18);
t65 = sin(t75);
t187 = t65 * t34;
t103 = qJD(1) ^ 2;
t97 = cos(qJ(3));
t58 = pkin(4) * t97 - pkin(1);
t98 = cos(qJ(2));
t171 = t58 * t98;
t92 = sin(qJ(3));
t93 = sin(qJ(2));
t162 = t92 * t93;
t55 = pkin(4) * t162;
t47 = t55 + pkin(12);
t25 = -t47 + t171;
t154 = t103 * t25;
t182 = qJD(4) * t65;
t44 = t94 * pkin(8) + t99 * pkin(10);
t45 = t99 * pkin(8) - t94 * pkin(10);
t21 = t44 * t84 + t83 * t45;
t186 = t182 * t21;
t80 = t99 ^ 2;
t67 = -t80 / 0.2e1;
t142 = t67 + 0.1e1 / 0.4e1;
t167 = t83 * t84;
t63 = t80 - 0.1e1 / 0.2e1;
t123 = t63 * pkin(10);
t159 = t99 * t94;
t32 = pkin(8) * t159 + t123;
t140 = pkin(10) * t159;
t33 = t63 * pkin(8) - t140;
t73 = t84 ^ 2;
t119 = -t32 * t167 + t33 * t73;
t66 = cos(t75);
t62 = t66 ^ 2;
t152 = t103 * t62;
t22 = -t83 * t44 + t45 * t84;
t35 = -t83 * t94 + t84 * t99;
t52 = -t159 / 0.2e1;
t131 = pkin(8) * t52 + t142 * pkin(10);
t135 = t33 * t167;
t4 = t32 * t73 + t131 + t135;
t177 = t80 / 0.2e1;
t134 = t177 - 0.1e1 / 0.4e1;
t121 = t134 * pkin(8);
t5 = t140 / 0.2e1 - t121 + t119;
t185 = ((-t25 * t35 + 0.4e1 * t4 * t65) * t103 - qJD(1) * qJD(4) * t22) * t66 - 0.4e1 * t5 * t152 + t154 * t187 + (-t25 * qJD(4) + (-t99 * t45 + 0.2e1 * t119) * qJD(1)) * qJD(1);
t117 = t92 * t34 - t35 * t97;
t157 = pkin(4) * qJD(3);
t17 = t34 * t97 + t92 * t35;
t161 = t92 * t98;
t40 = t97 * t93 + t161;
t181 = ((pkin(4) * t161 + t58 * t93) * qJD(2) + t40 * t157) * qJD(1);
t7 = -pkin(1) * t34 + t17 * pkin(4);
t8 = pkin(1) * t35 + t117 * pkin(4);
t184 = -((t93 * t7 + t8 * t98) * qJD(2) - (-t117 * t98 - t17 * t93) * t157) * t182 + t181 + qJD(4) * ((t7 * t98 - t93 * t8) * qJD(2) + (-t117 * t93 + t17 * t98) * t157) * t66;
t178 = -0.4e1 * t103;
t176 = pkin(1) * t98;
t90 = cos(pkin(17));
t175 = pkin(3) * t90;
t88 = sin(pkin(17));
t174 = t88 * pkin(3);
t173 = t97 * pkin(1);
t76 = qJD(2) ^ 2 / 0.2e1;
t172 = t103 / 0.2e1;
t170 = t65 * t35;
t169 = t65 * t66;
t168 = t66 * t34;
t105 = pkin(10) ^ 2;
t106 = pkin(8) ^ 2;
t72 = -t105 + t106;
t53 = t72 * t80;
t87 = sin(pkin(18));
t166 = t87 * pkin(12);
t89 = cos(pkin(18));
t165 = t87 * t89;
t164 = t87 * t99;
t163 = t89 * t94;
t160 = t93 * t98;
t126 = pkin(8) * t140;
t158 = -0.2e1 * t126 + t53;
t156 = pkin(6) * t103;
t56 = pkin(12) + t176;
t153 = t103 * t56;
t82 = qJ(3) + qJ(2);
t69 = cos(t82);
t151 = t103 * t69;
t150 = t103 * t98;
t12 = t168 + t170;
t149 = t12 ^ 2 * t103;
t107 = pkin(4) ^ 2;
t78 = t97 ^ 2;
t148 = t78 * t107;
t100 = cos(pkin(14));
t95 = sin(pkin(14));
t147 = t95 * t100;
t109 = pkin(1) ^ 2;
t146 = -0.2e1 * pkin(4) * t173 + t109;
t74 = qJD(2) + qJD(3);
t145 = qJD(1) * t74;
t91 = sin(qJ(4));
t144 = qJD(1) * t91;
t96 = cos(qJ(4));
t143 = qJD(1) * t96;
t141 = qJD(1) * qJD(2);
t139 = t99 * t174;
t64 = t106 / 0.2e1 - t105 / 0.2e1;
t24 = pkin(8) * t123 + t64 * t159;
t137 = t24 * t167;
t77 = t89 ^ 2;
t136 = t77 * t159;
t133 = pkin(1) * qJD(2) * t74;
t130 = t170 / 0.4e1;
t51 = t159 / 0.4e1;
t129 = t149 / 0.2e1;
t127 = t93 * t141;
t125 = -t106 / 0.4e1 + t105 / 0.4e1;
t41 = t100 * t94 - t95 * t99;
t42 = t100 * t99 + t95 * t94;
t124 = -t93 * t41 + t42 * t98;
t122 = t97 * t133;
t120 = pkin(1) * t127;
t13 = -t35 * t66 + t187;
t37 = t163 + t164;
t116 = t93 * t99 + t94 * t98;
t115 = t93 * t94 - t98 * t99;
t113 = -(-pkin(4) + t173) * qJD(2) * t157 + (t107 + t146) * t76 + t107 * qJD(3) ^ 2 / 0.2e1;
t112 = -t63 * t73 + t159 * t167 - 0.1e1 / 0.4e1;
t38 = -t87 * t94 + t89 * t99;
t111 = -t37 * t174 + t38 * t175 - pkin(12);
t104 = pkin(12) ^ 2;
t79 = t98 ^ 2;
t110 = -0.2e1 * t47 * t171 + t104 + t107 - t148 + (-t107 + t146 + 0.2e1 * t148) * t79 + 0.2e1 * pkin(12) * t55;
t81 = t100 ^ 2;
t71 = t109 * t76;
t68 = sin(t82);
t61 = t74 ^ 2 / 0.2e1;
t46 = t92 * t133;
t39 = -t97 * t98 + t162;
t36 = t56 ^ 2 * t172 + t71;
t23 = -t64 + t158;
t20 = t73 * t159 + t63 * t167 + t52;
t19 = t98 * t41 + t93 * t42;
t15 = t116 * t100 + t115 * t95;
t14 = t115 * t100 - t116 * t95;
t10 = -t13 * qJD(1) + qJD(4);
t9 = t111 - t176;
t3 = (t168 / 0.4e1 + t130) * qJD(4) + ((t67 - t112) * t169 + (t62 - 0.1e1 / 0.2e1) * t20) * qJD(1);
t1 = [0, 0, 0, 0, 0, t172, 0, 0, 0, 0, t93 ^ 2 * t172, t93 * t150, t127, t79 * t172, t98 * t141, t76, pkin(12) * t150, -t103 * pkin(12) * t93, 0, t104 * t172, t40 ^ 2 * t172, 0.2e1 * ((t78 - 0.1e1 / 0.2e1) * t160 + (t79 - 0.1e1 / 0.2e1) * t97 * t92) * t103, -t40 * t145, t39 ^ 2 * t172, t39 * t145, t61, t39 * t153 - t122, t40 * t153 + t46, -t120, t36, t129, (t20 * t62 - (t177 + t112) * t169 + t73 * t52 - t134 * t167 + t51) * t178, 0, t13 ^ 2 * t172, 0, 0, -t13 * t154, t12 * t154, t181, t110 * t172 + t113, t96 ^ 2 * t129, -t91 * t96 * t149, 0.4e1 * t3 * t143, t91 ^ 2 * t129, -0.4e1 * t3 * t144, t10 ^ 2 / 0.2e1, -t143 * t186 + t184 * t91 - t185 * t96, t144 * t186 + t184 * t96 + t185 * t91, (t4 * t62 + (t5 * t65 + t25 * t34 / 0.4e1) * t66 + t25 * t130 + t131 * t73 - t135 / 0.2e1 + pkin(8) * t51 + (t80 / 0.4e1 - 0.1e1 / 0.4e1) * pkin(10)) * t178, 0.2e1 * (t23 * t73 - t64 * t80 - t125 + t126 - 0.2e1 * t137) * t152 + ((t24 * t73 + t23 * t167 / 0.2e1 + t125 * t159 - pkin(10) * t121) * t65 - t25 * t22 / 0.4e1) * t66 * t178 - t21 * t65 * t154 + ((0.4e1 * t126 + t72 - 0.2e1 * t53) * t73 + 0.4e1 * t137 + t105 + t110 + t158) * t172 + t113, t15 ^ 2 * t172, 0.4e1 * t103 * ((-t63 * t147 + t81 * t159 + t52) * t79 + (t147 * t159 + t63 * t81 + t142) * t160 + t81 * t52 + t134 * t147 + t51), t19 * t141, t14 ^ 2 * t172, t124 * t141, t76, -t124 * t156, t19 * t156, 0, pkin(6) ^ 2 * t172, t37 ^ 2 * t172, -0.2e1 * t103 * (t63 * t165 + t136 + t52), 0, t38 ^ 2 * t172, 0, 0, -t38 * t153, -t37 * t153, -t120, t36, t68 ^ 2 * t172, t68 * t151, -t68 * t145, t69 ^ 2 * t172, -t69 * t145, t61, t9 * t151 - t122, -t9 * t103 * t68 + t46, -t120, t71 + (t79 * t109 - 0.2e1 * t111 * t176 - 0.4e1 * (t136 * t174 + (pkin(12) * t99 / 0.2e1 + t63 * t87 * t174) * t89 - t94 * (t139 + t166) / 0.2e1) * t175 + 0.2e1 * pkin(3) * (pkin(3) * t164 + t88 * pkin(12)) * t163 + 0.2e1 * t139 * t166 + t104 + (0.4e1 * (-t159 * t165 + t63 * t77 + t142) * t90 ^ 2 + (-0.2e1 * t80 + 0.1e1) * t77 + t80) * pkin(3) ^ 2) * t172;];
T_reg = t1;
