% Calculate inertial parameters regressor of coriolis matrix for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% cmat_reg [(2*2)x(2*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbar1turnTE_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:19:53
% EndTime: 2020-04-12 19:20:01
% DurationCPUTime: 2.34s
% Computational Cost: add. (26783->148), mult. (36879->339), div. (1238->20), fcn. (10269->4), ass. (0->167)
t86 = sin(qJ(2));
t195 = pkin(2) * t86;
t157 = pkin(1) * t195;
t87 = cos(qJ(2));
t194 = pkin(2) * t87;
t156 = pkin(1) * t194;
t96 = pkin(1) ^ 2;
t174 = -0.2e1 * t156 + t96;
t95 = pkin(2) ^ 2;
t79 = t95 + t174;
t77 = 0.1e1 / t79 ^ 2;
t141 = t77 * t157;
t204 = 0.1e1 / t79;
t161 = -t204 / 0.2e1;
t200 = -pkin(3) - pkin(4);
t73 = (pkin(2) - t200) * (pkin(2) + t200) + t174;
t199 = -pkin(3) + pkin(4);
t74 = (pkin(2) - t199) * (pkin(2) + t199) + t174;
t181 = t73 * t74;
t97 = sqrt(-t181);
t177 = t87 * t97;
t124 = pkin(1) * pkin(2) * (-t73 - t74);
t63 = t86 * t124;
t67 = 0.1e1 / t97;
t182 = t67 * t63;
t209 = pkin(3) ^ 2;
t210 = pkin(4) ^ 2;
t173 = t209 - t210;
t76 = t79 - t173;
t69 = t76 * t195;
t82 = pkin(1) - t194;
t29 = t69 + (-t177 + (0.2e1 * t82 * pkin(1) - t182) * t86) * pkin(2);
t178 = t86 * t97;
t65 = pkin(2) * t178;
t58 = t76 * t82 - t65;
t108 = t58 * t141 + t29 * t161;
t90 = 0.1e1 / pkin(4);
t17 = t108 * t90;
t101 = t58 ^ 2;
t53 = 0.1e1 / t101;
t59 = t82 * t97 + t69;
t184 = t53 * t59;
t160 = t204 / 0.2e1;
t197 = pkin(1) * t95;
t85 = t86 ^ 2;
t153 = t85 * t197;
t175 = t76 * t194 + t65;
t31 = t82 * t182 + 0.2e1 * t153 + t175;
t107 = -t59 * t141 + t31 * t160;
t20 = t107 * t90;
t52 = 0.1e1 / t58;
t187 = t20 * t52;
t211 = -0.2e1 * t17 * t184 - 0.2e1 * t187;
t193 = pkin(3) * t79;
t208 = 0.1e1 / t193;
t180 = t85 * t96;
t151 = pkin(2) * t180;
t198 = pkin(1) * t87;
t64 = pkin(1) * t178;
t75 = t79 + t173;
t176 = t75 * t198 + t64;
t83 = -pkin(2) + t198;
t30 = -t83 * t182 + 0.2e1 * t151 + t176;
t71 = pkin(1) * t86 * t75;
t60 = -t83 * t97 + t71;
t183 = t60 * t30;
t158 = -0.2e1 * t83 * pkin(2);
t28 = t71 + (-t177 + (t158 - t182) * t86) * pkin(1);
t57 = -t75 * t83 - t64;
t49 = 0.1e1 / t57;
t99 = t57 ^ 2;
t50 = 0.1e1 / t99;
t186 = t28 * t49 * t50;
t56 = t60 ^ 2;
t40 = t50 * t56 + 0.1e1;
t207 = (t50 * t183 - t56 * t186) / t40 ^ 2;
t205 = 0.1e1 / pkin(3);
t203 = 0.2e1 * t67;
t147 = t95 * t180;
t61 = t87 * t124 - 0.4e1 * t147;
t202 = -t61 / 0.2e1;
t201 = t77 / 0.4e1;
t196 = pkin(2) * t77;
t55 = t59 ^ 2;
t39 = t53 * t55 + 0.1e1;
t35 = 0.1e1 / t39;
t192 = pkin(4) * t35;
t191 = t28 * t205;
t143 = t205 * t160;
t125 = t87 * t143;
t127 = t86 * t143;
t32 = t60 * t125 + t57 * t127;
t120 = t205 * t141;
t110 = t60 * t120;
t118 = t85 * t205 * pkin(1) * t196;
t144 = t205 * t161;
t128 = t60 * t144;
t130 = t30 * t143;
t8 = -t57 * t118 + (t28 * t143 + t128) * t86 + (t57 * t143 - t110 + t130) * t87;
t190 = t32 * t8;
t126 = t87 * t144;
t33 = t57 * t126 + t60 * t127;
t139 = t77 * t156;
t119 = t205 * t139;
t9 = t28 * t126 + (t125 - t118) * t60 + (t130 + (t143 + t119) * t57) * t86;
t189 = t33 * t9;
t188 = t57 * t205;
t185 = t29 * t52 * t53;
t179 = t86 * t63;
t3 = -t32 * t9 - t33 * t8;
t172 = t3 * qJD(1);
t6 = (-t32 * t86 + t8 * t87) * pkin(2);
t171 = t6 * qJD(1);
t7 = (-t33 * t86 + t87 * t9) * pkin(2);
t170 = t7 * qJD(1);
t169 = qJD(2) * t90;
t78 = t204 * t77;
t140 = t78 * t157;
t91 = 0.1e1 / t210;
t10 = (-t59 * t58 * t140 + (t31 * t58 / 0.4e1 + t59 * t29 / 0.4e1) * t77) * t91;
t168 = t10 * qJD(1);
t113 = -t140 / 0.2e1;
t14 = (t58 * t29 * t201 + t101 * t113) * t91;
t167 = t14 * qJD(1);
t15 = (t59 * t31 * t201 + t55 * t113) * t91;
t166 = t15 * qJD(1);
t150 = t96 * t195;
t115 = t77 * t90 * t150;
t132 = pkin(1) * t90 * t160;
t16 = -t58 * t115 + t29 * t132;
t165 = t16 * qJD(1);
t19 = -t59 * t115 + t31 * t132;
t164 = t19 * qJD(1);
t80 = t87 ^ 2 - t85;
t163 = t80 * qJD(1);
t162 = t87 * qJD(2);
t159 = pkin(4) / t160;
t37 = 0.1e1 / t40;
t155 = t37 * t193;
t154 = t79 * t192;
t152 = t86 * t197;
t148 = 0.4e1 * t67 / t181 * t63 ^ 2;
t146 = qJD(1) * t190;
t145 = qJD(1) * t189;
t142 = t86 * t162;
t81 = t86 * t87 * qJD(1);
t138 = t77 * t152;
t137 = t49 * t155;
t136 = t50 * t155;
t133 = t78 * t147;
t131 = t28 * t144;
t129 = t57 * t144;
t123 = t95 * t81;
t122 = -t148 / 0.4e1;
t121 = pkin(3) * t37 * t157;
t116 = t17 * t59 * t159;
t114 = t60 * t136;
t18 = -0.2e1 * t57 * t120 + t204 * t191;
t21 = t30 * t208 - 0.2e1 * t110;
t5 = -t18 * t114 + t21 * t137 + 0.1e1;
t111 = t205 * t5 * t138;
t109 = 0.8e1 * t205 * t133;
t106 = t86 * t122 + (t86 * t202 - t87 * t63) * t203;
t4 = t154 * t211;
t2 = ((t83 * t122 - t71 + (pkin(1) * t97 + 0.6e1 * t150) * t87 + (pkin(1) * t179 + t83 * t202) * t203) * t208 - 0.4e1 * t30 * t120 + (t109 - 0.2e1 * t119) * t60) * t137 - ((0.4e1 * t151 + t176) * t208 + t57 * t109 + ((t87 * t158 + t106) * t208 + (-0.2e1 * t87 * t188 - 0.4e1 * t86 * t191) * t196) * pkin(1)) * t114 + (-t28 * t136 + (0.2e1 * t121 + 0.1e1 / t144 * t207) * t49) * t21 + (-t30 * t136 + (-0.2e1 * t121 * t50 + (t37 * t186 + t50 * t207) / t143) * t60) * t18;
t1 = 0.2e1 * t35 * t116 * t185 + 0.2e1 * t157 * t192 * t211 + 0.2e1 * (t116 * t53 + t159 * t187) / t39 ^ 2 * (t31 * t184 - t55 * t185) + 0.2e1 * ((-t17 * t31 + t20 * t29) * t53 + (-((t82 * t148 / 0.4e1 - t69 + (pkin(2) * t97 + 0.6e1 * t152) * t87 + (pkin(2) * t179 + t82 * t61 / 0.2e1) * t203) * t160 - 0.2e1 * t31 * t141 + (0.4e1 * t133 - t139) * t59) * t52 - ((0.4e1 * t153 + t175) * t161 - 0.4e1 * t58 * t133 + (t106 * t161 + (-t82 * t87 * t204 + (0.2e1 * t29 * t86 + t58 * t87) * t77) * pkin(1)) * pkin(2)) * t184) * t90) * t154;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t80 * qJD(2), 0, -t142, 0, 0, 0, 0, 0, 0, qJD(2) * t190, t3 * qJD(2), 0, qJD(2) * t189, 0, 0, t7 * qJD(2), t6 * qJD(2), 0, -t95 * t142, t15 * qJD(2), -t10 * qJD(2), 0, t14 * qJD(2), 0, 0, -t16 * qJD(2), -t19 * qJD(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t163, t162, -t81, -t86 * qJD(2), 0, 0, 0, 0, 0, t146, t172, (-t2 * t32 - t5 * t8) * qJD(2), t145, (t2 * t33 + t5 * t9) * qJD(2), 0, t170, t171, ((t60 * t205 * t33 + t32 * t188) * t138 + (t144 * t30 * t33 + t128 * t9 + t129 * t8 + t131 * t32) * pkin(2)) * qJD(2), -t123, t166, -t168, (t59 * t1 * t160 + t107 * t4) * t169, t167, (t58 * t1 * t161 + t108 * t4) * t169, 0, -t165, -t164, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t163, 0, t81, 0, 0, 0, 0, 0, 0, -t146, -t172, 0, -t145, 0, 0, -t170, -t171, 0, t123, -t166, t168, 0, -t167, 0, 0, t165, t164, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t2 * qJD(2), (t57 * t111 + (t129 * t2 + t131 * t5) * pkin(2)) * qJD(2), (-t60 * t111 + (t143 * t2 * t60 + t130 * t5) * pkin(2)) * qJD(2), 0, ((t183 / 0.4e1 + t57 * t28 / 0.4e1) * t77 + (-t56 / 0.2e1 - t99 / 0.2e1) * t140) * t95 / t209 * qJD(2), 0, 0, 0, 0, 0, t4 * t1 * qJD(2), 0, 0, 0, 0;];
cmat_reg = t11;
