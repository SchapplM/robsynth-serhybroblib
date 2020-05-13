% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh1m2TE_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_inertiaDJ_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:47:17
% EndTime: 2020-05-01 20:47:33
% DurationCPUTime: 2.43s
% Computational Cost: add. (813->184), mult. (2093->395), div. (0->0), fcn. (1748->22), ass. (0->151)
t103 = cos(qJ(2));
t102 = cos(qJ(3));
t98 = sin(qJ(2));
t152 = t102 * t98;
t97 = sin(qJ(3));
t130 = t97 * t152;
t142 = qJD(2) * t103;
t85 = t102 ^ 2;
t167 = 0.4e1 * t85 - 0.2e1;
t176 = t97 ^ 2 - t85;
t86 = t103 ^ 2;
t181 = -t176 * qJD(3) * (0.4e1 * t86 - 0.2e1) - (0.8e1 * qJD(3) * t130 - t167 * t142) * t103;
t144 = t102 * t103;
t156 = t97 * t98;
t49 = t144 - t156;
t106 = t49 * qJD(3);
t147 = qJD(2) * t98;
t70 = pkin(5) * t97 + pkin(1);
t23 = (t102 * t142 + t106) * pkin(5) - t70 * t147;
t104 = cos(pkin(18));
t91 = sin(pkin(20));
t94 = cos(pkin(20));
t99 = sin(pkin(18));
t42 = -t104 * t91 + t99 * t94;
t45 = t104 * t94 + t99 * t91;
t80 = pkin(22) + pkin(21);
t74 = sin(t80);
t75 = cos(t80);
t21 = t42 * t75 - t74 * t45;
t180 = t21 ^ 2;
t148 = t97 * t103;
t47 = t148 + t152;
t107 = t47 * qJD(3);
t22 = t70 * t142 + (t102 * t147 + t107) * pkin(5);
t172 = 0.2e1 * t22;
t178 = t21 * t172;
t84 = t98 ^ 2;
t155 = t84 - t86;
t100 = sin(pkin(17));
t105 = cos(pkin(17));
t48 = -t104 * t100 + t99 * t105;
t50 = t100 * t99 + t105 * t104;
t112 = t50 * t103 + t98 * t48;
t25 = t112 * qJD(2);
t92 = sin(pkin(19));
t95 = cos(pkin(19));
t113 = t102 * t92 + t97 * t95;
t108 = t113 * qJD(2);
t36 = t113 * qJD(3);
t40 = -t103 * t95 + t92 * t98;
t41 = t92 * t103 + t98 * t95;
t116 = t41 * t102 - t40 * t97;
t34 = -t95 * t142 + t92 * t147;
t35 = t41 * qJD(2);
t174 = t116 * qJD(3) + t35 * t102 - t34 * t97;
t173 = -2 * pkin(15);
t111 = t103 * t48 - t98 * t50;
t24 = t111 * qJD(2);
t171 = 0.2e1 * t24;
t87 = t104 ^ 2;
t169 = -0.2e1 * t87;
t168 = 0.2e1 * t87;
t166 = 0.8e1 * t85 - 0.4e1;
t163 = 0.4e1 * t87 - 0.2e1;
t115 = -t40 * t102 - t41 * t97;
t12 = t115 * qJD(3) - t34 * t102 - t35 * t97;
t162 = pkin(2) * t12;
t161 = pkin(2) * t36;
t114 = t102 * t95 - t97 * t92;
t37 = t114 * qJD(3);
t160 = pkin(2) * t37;
t159 = t74 * t42;
t88 = qJ(3) + qJ(2);
t77 = sin(t88);
t82 = qJD(2) + qJD(3);
t158 = t82 * t77;
t78 = cos(t88);
t56 = t82 * t78;
t157 = t91 * t94;
t68 = pkin(1) * t98 - pkin(15);
t154 = pkin(1) * qJD(3);
t153 = t102 * t97;
t151 = t103 * t98;
t150 = t104 * t99;
t96 = sin(qJ(4));
t145 = qJD(4) * t96;
t73 = t87 - 0.1e1 / 0.2e1;
t101 = cos(qJ(4));
t141 = qJD(4) * t101;
t132 = -t87 / 0.2e1 + 0.1e1 / 0.4e1;
t69 = t75 ^ 2;
t81 = t94 ^ 2;
t140 = -0.8e1 * qJD(4) * ((t73 * t157 + (-t81 + 0.1e1 / 0.2e1) * t150) * t69 + (t150 * t157 + t73 * t81 + t132) * t74 * t75 + t81 * t150 / 0.2e1 + t132 * t157 - t150 / 0.4e1);
t139 = -0.2e1 * t151;
t138 = 0.2e1 * t151;
t135 = pkin(11) * t150;
t137 = pkin(9) * t169 + pkin(9) - 0.2e1 * t135;
t136 = pkin(9) * t150;
t134 = (t73 * pkin(11) - t136) * t157;
t133 = t77 * t56;
t131 = t97 * t154;
t129 = pkin(1) * t142;
t71 = t102 * t154;
t128 = pkin(15) * t147;
t127 = t97 * t144;
t125 = pkin(15) * t142;
t124 = t98 * t142;
t123 = pkin(1) * t86 + pkin(15) * t98 - pkin(1);
t121 = 0.2e1 * t129;
t120 = t98 * t127;
t119 = t150 * t151;
t118 = t180 * t96 * t141;
t32 = -pkin(5) * t144 + t70 * t98;
t89 = sin(pkin(22));
t93 = cos(pkin(22));
t110 = -t104 * t89 + t99 * t93;
t31 = -pkin(15) + t32;
t6 = t22 * (t45 * t75 + t159);
t109 = qJD(2) * t114;
t72 = pkin(1) * t147;
t64 = 0.2e1 * t71;
t63 = -0.2e1 * t131;
t62 = pkin(5) * t71;
t61 = 0.2e1 * t62;
t58 = -0.2e1 * t124;
t57 = 0.2e1 * t124;
t54 = t104 * pkin(9) + t99 * pkin(11);
t53 = -t99 * pkin(9) + t104 * pkin(11);
t52 = t167 * t98;
t51 = t68 * t121;
t46 = 0.2e1 * t155 * qJD(2);
t44 = t93 * t104 + t99 * t89;
t33 = pkin(5) * t152 + t70 * t103;
t29 = t49 * qJD(2) + t106;
t28 = t47 * qJD(2) + t107;
t27 = t49 * t104 + t99 * t47;
t26 = -t47 * t104 + t99 * t49;
t17 = t33 * t104 + t32 * t99;
t16 = t32 * t104 - t33 * t99;
t15 = (cos(pkin(21)) * t44 + t110 * sin(pkin(21))) * pkin(4) + t68;
t14 = -t29 * t104 - t99 * t28;
t13 = -t28 * t104 + t99 * t29;
t10 = t23 * t104 + t22 * t99;
t9 = t22 * t104 - t23 * t99;
t8 = (-t91 * t53 + t54 * t94) * t75 + (-t53 * t94 - t91 * t54) * t74 + t31;
t5 = (t26 * t94 - t27 * t91) * t75 - t74 * (t26 * t91 + t27 * t94);
t4 = (-t91 * t16 + t17 * t94) * t75 - t74 * (t16 * t94 + t91 * t17);
t3 = ((t163 * pkin(9) + 0.4e1 * t135) * t81 - 0.4e1 * t134 + t137) * t69 + (((-t163 * pkin(11) + 0.4e1 * t136) * t81 - 0.4e1 * (t73 * pkin(9) + t135) * t157 + pkin(11) * t168 - 0.2e1 * t136 - pkin(11)) * t74 + t45 * t31) * t75 + t31 * t159 + t137 * t81 + 0.2e1 * t134 + t104 * t54;
t2 = (-t13 * t91 + t14 * t94) * t75 - t74 * (t13 * t94 + t14 * t91);
t1 = (t10 * t94 - t91 * t9) * t75 - t74 * (t91 * t10 + t9 * t94);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t46, 0, t57, 0, 0, -0.2e1 * t125, 0.2e1 * t128, 0, 0, 0.2e1 * t47 * t29, (-t52 - 0.8e1 * t127) * t147 + t181, 0, -0.2e1 * t49 * t28, 0, 0, 0.2e1 * (-t123 * t102 + t68 * t148) * qJD(3) + 0.2e1 * (-pkin(15) * t148 + t68 * t152 + (-t86 * t102 + t97 * t138) * pkin(1)) * qJD(2), 0.2e1 * (t123 * t97 + t68 * t144) * qJD(3) + 0.2e1 * (-pkin(15) * t144 - t68 * t156 + (t102 * t138 + t86 * t97) * pkin(1)) * qJD(2), 0, t51, 0, 0, 0, 0, 0, 0, 0.2e1 * t6, -t178, 0, t31 * t172, -0.2e1 * t118, 0.2e1 * (-t101 ^ 2 + t96 ^ 2) * t180 * qJD(4), t96 * t140, 0.2e1 * t118, t101 * t140, 0, 0.2e1 * t101 * t6 - 0.2e1 * t3 * t145, -0.2e1 * t3 * t141 - 0.2e1 * t6 * t96, t178, t8 * t172, t112 * t171, (0.8e1 * t119 + (-0.8e1 * (t139 * t87 - t155 * t150 + t151) * t100 + 0.4e1 * (t168 * t84 + t169 * t86 - 0.4e1 * t119 - t155) * t105) * t105 - t155 * t163) * qJD(2), 0, -0.2e1 * t111 * t25, 0, 0, 0.2e1 * pkin(14) * t25, pkin(14) * t171, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t121, -0.2e1 * t110 * t129, 0, t51, 0.2e1 * t116 * t12, 0.8e1 * qJD(2) * t120 + t52 * t147 + (0.8e1 * (((t85 - 0.1e1 / 0.2e1) * t139 + t155 * t153) * qJD(2) + (t176 * t151 + (-0.2e1 * t86 + 0.1e1) * t153) * qJD(3)) * t92 + ((-t166 * t84 + (t103 * t166 - 0.16e2 * t130) * t103) * qJD(2) + (-0.16e2 * t120 - t176 * (0.8e1 * t86 - 0.4e1)) * qJD(3)) * t95) * t95 - t181, 0, -0.2e1 * t115 * t174, 0, 0, t174 * t173, t12 * t173, 0, 0, t58, t46, 0, t57, 0, 0, -0.2e1 * t125 + 0.2e1 * (-t86 * t109 + ((t37 + t109) * t98 + (0.2e1 * t108 + t36) * t103) * t98) * pkin(2), 0.2e1 * t128 + 0.2e1 * (-t84 * t108 + ((t36 + t108) * t103 + (0.2e1 * t109 + t37) * t98) * t103) * pkin(2), 0, 0.2e1 * (pkin(15) + (-t47 * t92 + t49 * t95) * pkin(2)) * (-t28 * t95 - t29 * t92) * pkin(2), 0.2e1 * t133, 0.2e1 * (-t77 ^ 2 + t78 ^ 2) * t82, 0, -0.2e1 * t133, 0, 0, -0.2e1 * t129 * t78 + 0.2e1 * t15 * t158, 0.2e1 * t129 * t77 + 0.2e1 * t15 * t56, 0, t15 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, 0, -t142, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t28, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, -t33 * t141 - t23 * t96, -t101 * t23 + t33 * t145, 0, 0, 0, 0, t24, 0, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, t12, 0, -t174, 0, 0, 0, 0, 0, 0, 0, -t147, 0, -t142, 0, 0, 0, -t162, 0, 0, 0, t56, 0, -t158, 0, 0, 0, t72, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t160, 0.2e1 * t161, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(5) * t29, 0, 0, 0, 0, 0, 0, 0, (-t47 * t141 - t29 * t96) * pkin(5), (-t101 * t29 + t47 * t145) * pkin(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t174, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, 0, 0, 0, t56, 0, -t158, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t161, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t131, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t145, 0, t21 * t141, 0, t101 * t22 - t8 * t145, -t8 * t141 - t22 * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t96 - t4 * t141, -t1 * t101 + t145 * t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t5 * t141 + t2 * t96) * pkin(5), (t101 * t2 - t145 * t5) * pkin(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
