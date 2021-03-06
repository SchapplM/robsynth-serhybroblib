% Calculate inertial parameters regressor of coriolis joint torque vector for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% tauc_reg [12x(12*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = picker2Dm1OL_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:45:15
% EndTime: 2020-05-11 05:45:24
% DurationCPUTime: 1.54s
% Computational Cost: add. (1709->202), mult. (3918->308), div. (0->0), fcn. (2176->14), ass. (0->141)
t114 = cos(qJ(2));
t176 = pkin(1) * qJD(2);
t145 = qJD(1) * t176;
t135 = t114 * t145;
t177 = pkin(1) * qJD(1);
t144 = qJD(6) * t177;
t186 = t114 * t144 + t135;
t151 = qJD(2) + qJD(6);
t113 = cos(qJ(3));
t152 = qJD(2) + qJD(3);
t160 = qJD(3) * t113;
t100 = qJD(1) + qJD(2);
t94 = t114 * t177;
t70 = t100 * pkin(2) + t94;
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t146 = t108 * t177;
t84 = t107 * t146;
t28 = -t113 * t135 + t152 * t84 - t70 * t160;
t185 = pkin(1) * t151;
t112 = cos(qJ(4));
t106 = sin(qJ(4));
t136 = t108 * t145;
t139 = t106 * t146;
t179 = qJD(4) * t139 + t106 * t136;
t69 = t100 * pkin(3) + t94;
t26 = (qJD(4) * t69 + t135) * t112 - t179;
t105 = sin(qJ(6));
t111 = cos(qJ(6));
t118 = -t105 * t146 * t151 + t186 * t111;
t101 = sin(qJ(10));
t102 = cos(qJ(10));
t153 = qJD(10) * t106;
t154 = qJD(10) * t102;
t170 = t101 * t106;
t163 = t108 * t112;
t124 = -t106 * t114 - t163;
t60 = t124 * t177;
t167 = t106 * t108;
t123 = t112 * t114 - t167;
t62 = t123 * t177;
t90 = t112 * pkin(3) + pkin(4);
t184 = -t101 * t60 - t102 * t62 + t90 * t154 - (t101 * t153 + (-t102 * t112 + t170) * qJD(4)) * pkin(3);
t155 = qJD(10) * t101;
t169 = t102 * t106;
t183 = -t101 * t62 + t102 * t60 + t90 * t155 + (t102 * t153 + (t101 * t112 + t169) * qJD(4)) * pkin(3);
t103 = sin(qJ(9));
t109 = cos(qJ(9));
t156 = qJD(9) * t109;
t157 = qJD(9) * t103;
t168 = t103 * t107;
t162 = t108 * t113;
t122 = t107 * t114 + t162;
t61 = t122 * t177;
t63 = -t113 * t94 + t84;
t91 = -t113 * pkin(2) + pkin(6);
t182 = -t103 * t61 - t109 * t63 + t91 * t156 - (-t107 * t157 + (t109 * t113 - t168) * qJD(3)) * pkin(2);
t165 = t107 * t109;
t181 = -t103 * t63 + t109 * t61 + t91 * t157 + (-t107 * t156 + (-t103 * t113 - t165) * qJD(3)) * pkin(2);
t175 = pkin(1) * qJD(8);
t54 = t106 * t69 + t112 * t146;
t174 = t101 * t54;
t173 = t102 * t54;
t138 = t113 * t146;
t55 = t107 * t70 + t138;
t172 = t103 * t55;
t171 = t109 * t55;
t166 = t107 * t108;
t164 = t108 * t111;
t161 = qJD(3) * t107;
t159 = qJD(4) * t106;
t158 = qJD(4) * t112;
t116 = (t124 * qJD(2) - t108 * t158) * pkin(1);
t115 = qJD(1) * t116;
t27 = -t69 * t159 + t115;
t52 = t112 * t69 - t139;
t96 = qJD(4) + t100;
t48 = t96 * pkin(4) + t52;
t150 = t101 * t27 + t102 * t26 + t48 * t154;
t29 = qJD(3) * t138 + t107 * t135 + t113 * t136 + t70 * t161;
t53 = -t113 * t70 + t84;
t97 = qJD(1) + t152;
t49 = t97 * pkin(6) + t53;
t149 = t103 * t29 + t109 * t28 + t49 * t156;
t99 = qJD(1) + qJD(8);
t148 = t99 * t177;
t147 = t99 * t175;
t143 = qJD(1) * t175;
t98 = t114 * pkin(1);
t93 = t98 + pkin(2);
t142 = pkin(1) * t166 - t113 * t93;
t141 = t186 * t105 + t111 * t136 + t144 * t164;
t92 = t98 + pkin(3);
t134 = -pkin(1) * t167 + t112 * t92;
t133 = (-qJD(2) + t100) * t177;
t132 = (-qJD(1) - t100) * t176;
t131 = (-pkin(3) * t96 - t69) * qJD(4);
t56 = pkin(4) + t134;
t64 = pkin(1) * t163 + t106 * t92;
t130 = t101 * t64 - t102 * t56;
t129 = t101 * t56 + t102 * t64;
t57 = pkin(6) + t142;
t65 = -pkin(1) * t162 - t107 * t93;
t128 = t103 * t65 - t109 * t57;
t127 = t103 * t57 + t109 * t65;
t2 = t101 * t26 - t102 * t27 + t54 * t154 + t48 * t155;
t4 = t103 * t28 - t109 * t29 - t55 * t156 + t49 * t157;
t126 = t105 * t114 + t164;
t125 = t105 * t108 - t111 * t114;
t121 = t55 * t157 + t149;
t1 = t54 * t155 - t150;
t110 = cos(qJ(8));
t104 = sin(qJ(8));
t95 = qJD(1) + t151;
t89 = qJD(9) + t97;
t88 = qJD(10) + t96;
t86 = t110 * t143;
t85 = t104 * t143;
t59 = t126 * t177;
t58 = t125 * t177;
t47 = t126 * t185;
t46 = t125 * t185;
t41 = t93 * t161 + (t122 * qJD(2) + t108 * t160) * pkin(1);
t40 = -t93 * t160 + (t108 * t161 + (-t113 * t114 + t166) * qJD(2)) * pkin(1);
t39 = -t92 * t159 + t116;
t38 = t92 * t158 + (t123 * qJD(2) - t108 * t159) * pkin(1);
t19 = -t109 * t53 - t172;
t18 = t103 * t53 - t171;
t17 = -t102 * t52 + t174;
t16 = t101 * t52 + t173;
t15 = t103 * t49 - t171;
t14 = -t109 * t49 - t172;
t13 = t101 * t48 + t173;
t12 = -t102 * t48 + t174;
t11 = t58 * t95 + t118;
t10 = -t59 * t95 + t141;
t8 = t127 * qJD(9) + t103 * t40 - t109 * t41;
t7 = t128 * qJD(9) - t103 * t41 - t109 * t40;
t6 = t129 * qJD(10) + t101 * t38 - t102 * t39;
t5 = t130 * qJD(10) - t101 * t39 - t102 * t38;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t132, t114 * t132, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t97 + t29, -t40 * t97 - t28, 0, t29 * t142 + t28 * t65 - t55 * t40 + t53 * t41, 0, 0, 0, 0, 0, 0, t39 * t96 + t27, -t38 * t96 - t26, 0, t27 * t134 + t26 * t64 + t54 * t38 + t52 * t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t95 + t141, -t46 * t95 + t118, 0, -t59 * t46 + t58 * t47 + (t118 * t126 + t141 * t125) * pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * t147 + t85, t110 * t147 + t86, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t89 + t4, -t7 * t89 + t121, 0, t121 * t127 + t128 * t4 + t14 * t8 - t15 * t7, 0, 0, 0, 0, 0, 0, t6 * t88 + t2, -t5 * t88 - t1, 0, -t1 * t129 + t12 * t6 - t13 * t5 + t130 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t133, t114 * t133, 0, 0, 0, 0, 0, 0, 0, 0, (pkin(2) * t161 - t61) * t97 + t29, (pkin(2) * t160 + t63) * t97 - t28, 0, -t53 * t61 + t55 * t63 + (-t107 * t28 - t113 * t29 + (t107 * t53 + t113 * t55) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t106 * t131 - t60 * t96 + t115, t62 * t96 + (t131 - t135) * t112 + t179, 0, -t52 * t60 - t54 * t62 + (t106 * t26 + t112 * t27 + (-t106 * t52 + t112 * t54) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181 * t89 + t4, t182 * t89 + t121, 0, -t121 * (pkin(2) * t165 - t103 * t91) + t4 * (-pkin(2) * t168 - t109 * t91) + t182 * t15 + t181 * t14, 0, 0, 0, 0, 0, 0, t183 * t88 + t2, t184 * t88 - t1, 0, t1 * (-pkin(3) * t169 - t101 * t90) + t2 * (pkin(3) * t170 - t102 * t90) + t184 * t13 + t183 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t97 + t29, t53 * t97 - t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (pkin(6) * t157 - t18) * t89 + t4, t19 * t89 + (pkin(6) * t109 * t89 + t172) * qJD(9) + t149, 0, -t14 * t18 + t15 * t19 + (t121 * t103 - t4 * t109 + (t103 * t14 + t109 * t15) * qJD(9)) * pkin(6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t96 + t27, t52 * t96 - t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (pkin(4) * t155 - t16) * t88 + t2, t17 * t88 + (pkin(4) * t102 * t88 - t174) * qJD(10) + t150, 0, -t12 * t16 + t13 * t17 + (-t1 * t101 - t2 * t102 + (t101 * t12 + t102 * t13) * qJD(10)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * t148 + t85, -t110 * t148 + t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t89 + t4, t14 * t89 + t121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t88 + t2, t12 * t88 - t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg = t3;
