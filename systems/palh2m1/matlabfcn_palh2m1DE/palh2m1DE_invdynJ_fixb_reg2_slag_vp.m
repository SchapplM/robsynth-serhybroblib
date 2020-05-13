% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m1DE_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:31
% EndTime: 2020-05-02 23:52:36
% DurationCPUTime: 1.36s
% Computational Cost: add. (680->196), mult. (1544->316), div. (0->0), fcn. (1005->16), ass. (0->152)
t94 = sin(qJ(1));
t98 = cos(qJ(1));
t203 = g(1) * t94 - g(2) * t98;
t96 = cos(qJ(3));
t190 = t96 * pkin(3);
t62 = pkin(2) + t190;
t97 = cos(qJ(2));
t43 = t62 * t97;
t92 = sin(qJ(3));
t93 = sin(qJ(2));
t180 = t93 * t92;
t57 = pkin(3) * t180;
t175 = -t57 + t43;
t210 = pkin(1) + t175;
t19 = pkin(4) + t210;
t135 = 2 * qJD(1);
t189 = t97 * pkin(2);
t64 = pkin(1) + t189;
t213 = t135 * t64;
t161 = qJD(2) * t93;
t134 = pkin(2) * t161;
t152 = qJDD(1) * t64;
t212 = t134 * t135 - 0.2e1 * t152 - t203;
t177 = t97 * t92;
t178 = t96 * t93;
t115 = t177 + t178;
t176 = t97 * t96;
t29 = t176 - t180;
t112 = t29 * qJD(3);
t13 = t29 * qJD(2) + t112;
t141 = pkin(3) * t177;
t23 = t93 * t62 + t141;
t194 = g(2) * t94;
t195 = g(1) * t98;
t42 = t194 + t195;
t160 = qJD(2) * t97;
t9 = t62 * t160 + (-t92 * t161 + t112) * pkin(3);
t211 = t9 * qJD(2) + t23 * qJDD(2) + (qJD(3) * t13 + qJDD(3) * t115) * pkin(3) + t42;
t81 = qJD(2) + qJD(3);
t14 = t81 * t115;
t101 = qJD(2) ^ 2;
t209 = qJDD(2) * t93 + t101 * t97;
t192 = t92 * g(3);
t116 = -t42 * t96 + t192;
t76 = t96 * g(3);
t126 = t92 * t42 + t76;
t102 = qJD(1) ^ 2;
t170 = t102 * t64;
t191 = t96 * pkin(2);
t208 = qJDD(2) * t191 + t115 * t170 - t116 * t93 + t126 * t97;
t89 = qJ(2) + qJ(3);
t75 = cos(t89);
t207 = pkin(3) * t75;
t87 = t97 ^ 2;
t155 = t87 - 0.1e1 / 0.2e1;
t33 = -0.2e1 * t176 * t180;
t86 = t96 ^ 2;
t206 = t33 + (-t92 ^ 2 + t86) * t155;
t85 = t93 ^ 2;
t173 = t85 - t87;
t70 = t86 - 0.1e1 / 0.2e1;
t205 = -t173 * t70 + t33;
t202 = -2 * t102;
t201 = -2 * qJD(1);
t200 = 0.2e1 * qJDD(1);
t199 = pkin(3) * (pkin(3) + t191);
t74 = sin(t89);
t198 = pkin(3) * t74;
t197 = pkin(3) * t92;
t193 = g(3) * t97;
t188 = t102 / 0.2e1;
t172 = pkin(3) * qJD(3);
t24 = t29 * t172;
t187 = -t175 * qJD(2) - t24 + t9;
t95 = cos(qJ(4));
t186 = t211 * t95;
t184 = pkin(3) * t102;
t82 = qJD(1) + qJD(4);
t183 = t82 * t95;
t79 = qJDD(1) + qJDD(4);
t91 = sin(qJ(4));
t182 = t91 * t79;
t181 = t91 * t82;
t179 = t95 * t79;
t169 = t102 * t93;
t104 = pkin(3) ^ 2;
t68 = 0.2e1 * t89;
t168 = t104 * sin(t68);
t103 = 0.2e1 * qJ(2);
t105 = pkin(2) ^ 2;
t167 = t105 * sin(t103);
t8 = -t62 * t161 + (-t115 * qJD(3) - t92 * t160) * pkin(3);
t166 = t8 * qJD(1);
t165 = qJDD(1) / 0.2e1;
t164 = pkin(3) * qJDD(1);
t163 = qJD(1) * t91;
t162 = qJD(1) * t95;
t159 = qJD(4) * t82;
t11 = t23 * qJD(2) + t115 * t172;
t158 = t11 * qJD(4);
t157 = -qJD(1) - t82;
t156 = -qJD(4) + t82;
t154 = qJDD(1) * t115;
t153 = qJDD(1) * t29;
t150 = t19 * qJDD(1);
t149 = t93 * qJDD(1);
t148 = t97 * qJDD(1);
t147 = qJD(1) * qJD(2);
t146 = pkin(2) * t197;
t145 = pkin(3) * t194;
t144 = pkin(2) * t184;
t88 = qJ(3) + t103;
t73 = sin(t88);
t143 = pkin(3) * (qJD(3) + 0.2e1 * qJD(2)) * t73;
t142 = t81 * t198;
t139 = t74 * t184;
t138 = qJD(3) ^ 2 * t197;
t137 = pkin(1) * t200;
t136 = -0.2e1 * t161;
t133 = t102 * t29 * t115;
t132 = t97 * t169;
t130 = t19 * t163;
t129 = 0.2e1 * pkin(2) * t190 + t104 + t105;
t127 = t97 * t147;
t54 = t92 * t145;
t125 = pkin(2) * cos(t88) * t164 - qJD(1) * t81 * t168 + (cos(t68) * t104 + cos(t103) * t105) * t165;
t124 = qJD(2) * (-qJD(3) + t81);
t123 = qJD(3) * (-qJD(2) - t81);
t122 = qJD(4) * t157;
t121 = qJD(1) * t136;
t120 = t8 * t162 + (t150 + t158) * t95 + t211 * t91;
t119 = t93 * t127;
t56 = -0.2e1 * qJD(3) * t146;
t118 = -t116 * t97 - t126 * t93 + t29 * t170;
t114 = t56 / 0.2e1 - qJD(2) * t167;
t113 = t209 * pkin(2) + t42;
t111 = pkin(1) * t102 + t42;
t35 = t168 * t188;
t49 = t73 * t144;
t110 = g(3) * t43 + qJD(2) * t56 + t129 * qJDD(2) + qJDD(3) * t199 + t167 * t188 + t97 * t54 + t35 + t49;
t109 = t203 + t137;
t107 = t35 + t145 * t178 + t49 / 0.2e1 + qJDD(2) * t199 + t92 * t144 / 0.2e1 + t101 * t146 + qJDD(3) * t104 + (pkin(3) * t76 + t54) * t97;
t106 = pkin(1) ^ 2;
t99 = pkin(1) + pkin(4);
t78 = qJDD(2) + qJDD(3);
t71 = pkin(3) * t192;
t50 = pkin(1) * t139;
t38 = t99 * t139;
t26 = t115 * pkin(3);
t12 = (pkin(3) * t176 - t57) * qJD(2) + t24;
t4 = -t154 + (t29 * t81 - t13) * qJD(1);
t1 = [0, 0, 0, 0, 0, qJDD(1), t203, t42, 0, 0, t85 * qJDD(1) + 0.2e1 * t119, -0.2e1 * t173 * t147 + 0.2e1 * t93 * t148, -t209, t87 * qJDD(1) - 0.2e1 * t119, -qJDD(2) * t97 + t101 * t93, 0, pkin(1) * t121 + t109 * t97, -0.2e1 * pkin(1) * t127 - t109 * t93, t42, pkin(1) * t203 + qJDD(1) * t106, -(t13 * t201 - t154) * t115, 0.4e1 * (t155 * t96 * t92 + t93 * t70 * t97) * qJDD(1) + 0.4e1 * t205 * t147 + 0.4e1 * t206 * qJD(1) * qJD(3), -t115 * t78 - t81 * t13, (t14 * t201 + t153) * t29, t81 * t14 - t78 * t29, 0, -t14 * t213 - t212 * t29, t212 * t115 - t13 * t213, t113, (pkin(2) * t121 + t152 + t203) * t64, 0, 0, 0, qJDD(1), 0, 0, t200 * t210 + 0.2e1 * t166 + t203, 0, (t75 * t81 ^ 2 + t74 * t78) * pkin(3) + t113, t137 * t207 + (0.4e1 * pkin(1) * t189 + 0.2e1 * t106 + t129) * t165 + (-pkin(2) * t143 + 0.2e1 * (-t134 - t142) * pkin(1) + t114) * qJD(1) + t125 + t210 * t203, 0, 0, 0, 0, 0, t79, (t8 * t82 + t203) * t95 + (t122 * t91 + t179) * t19 + t120, t95 * t19 * t122 + (-t158 - t203 + t157 * t8 + (-qJDD(1) - t79) * t19) * t91 + t186, 0, t129 * t165 + (0.2e1 * t189 + t99 + 0.2e1 * t207) * qJDD(1) * t99 + (-0.2e1 * t99 * t142 + (t136 * t99 - t143) * pkin(2) + t114) * qJD(1) + t125 + t203 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, t173 * t102, -t149, t132, -t148, qJDD(2), t111 * t93 + t193, -g(3) * t93 + t111 * t97, 0, 0, -t133, t205 * t202, t4, t133, -t153, t78, (t92 * t123 + t29 * t169 + t96 * t78) * pkin(2) + t208, (-t115 * t169 + (-qJDD(2) - t78) * t92 + t96 * t123) * pkin(2) + t118, pkin(2) * t149, qJDD(2) * t105 + (t193 + (t42 + t170) * t93) * pkin(2), 0, 0, 0, 0, 0, 0, t23 * t102, 0, (t93 * pkin(2) + t198) * qJDD(1), t50 - (-t42 * t62 + t71) * t93 + t141 * t195 + (pkin(1) * t169 - t138) * pkin(2) + t110, 0, 0, 0, 0, 0, 0, t23 * t182 + (t23 * t183 + t187 * t91) * t82, t23 * t179 + (-t23 * t181 + t187 * t95) * t82, 0, t38 + t23 * t195 - (-t62 * t194 + t71) * t93 + (t169 * t99 - t138) * pkin(2) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t206 * t202, t4, t133, -t153, t78, t92 * pkin(2) * t124 + t208, (-qJDD(2) * t92 + t96 * t124) * pkin(2) + t118, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t102, 0, t74 * t164, t50 + (-g(3) * t180 + t115 * t195) * pkin(3) + t107, 0, 0, 0, 0, 0, 0, -(t91 * t12 - t162 * t26) * t82 + (t13 * t181 - (-t159 * t95 - t182) * t115) * pkin(3), (-t12 * t95 - t163 * t26) * t82 + (t13 * t183 - (t159 * t91 - t179) * t115) * pkin(3), 0, -g(3) * t57 + t26 * t195 + t107 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -qJD(4) * t130 - (t11 * t95 - t130) * t82 + t203 * t95 + t120, t156 * t19 * t162 + (t11 * t156 - t150 - t166 - t203) * t91 + t186, 0, 0;];
tau_reg = t1;
