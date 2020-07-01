% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tau_reg [5x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m1OL_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:12
% EndTime: 2020-06-30 17:45:21
% DurationCPUTime: 2.97s
% Computational Cost: add. (3605->299), mult. (8429->408), div. (0->0), fcn. (7025->14), ass. (0->176)
t135 = qJ(2) + qJ(3);
t129 = qJ(4) + t135;
t116 = sin(t129);
t117 = cos(t129);
t142 = sin(qJ(1));
t147 = cos(qJ(1));
t174 = g(1) * t147 + g(2) * t142;
t250 = g(3) * t117 + t116 * t174;
t143 = cos(qJ(5));
t208 = qJD(5) * t143;
t139 = sin(qJ(4));
t144 = cos(qJ(4));
t145 = cos(qJ(3));
t140 = sin(qJ(3));
t141 = sin(qJ(2));
t214 = qJD(1) * t141;
t197 = t140 * t214;
t146 = cos(qJ(2));
t213 = qJD(1) * t146;
t87 = -t145 * t213 + t197;
t168 = t140 * t146 + t145 * t141;
t88 = t168 * qJD(1);
t53 = t139 * t88 + t144 * t87;
t259 = t143 * t53;
t267 = t208 - t259;
t229 = pkin(2) * qJD(2);
t202 = t140 * t229;
t180 = t139 * t202;
t240 = t145 * pkin(2);
t122 = qJDD(2) * t240;
t131 = qJDD(2) + qJDD(3);
t74 = t131 * pkin(3) - qJD(3) * t202 + t122;
t191 = -qJD(4) * t180 + t139 * t74;
t241 = g(3) * t116;
t101 = qJD(1) * pkin(1) + pkin(2) * t213;
t67 = -t87 * pkin(3) + t101;
t265 = t174 * t117 - t67 * t53 - t191 - t241;
t212 = qJD(3) * t145;
t195 = qJD(2) * t212;
t206 = qJDD(2) * t140;
t156 = (t195 + t206) * pkin(2);
t132 = qJD(2) + qJD(3);
t93 = t132 * pkin(3) + t145 * t229;
t264 = (qJD(4) * t93 + t156) * t144;
t263 = -t264 + t265;
t171 = t139 * t87 - t144 * t88;
t262 = -t67 * t171 + t250;
t126 = qJD(4) + t132;
t138 = sin(qJ(5));
t46 = t138 * t126 + t143 * t171;
t125 = qJDD(4) + t131;
t210 = qJD(4) * t144;
t211 = qJD(4) * t139;
t205 = t141 * qJDD(1);
t207 = qJD(1) * qJD(2);
t196 = t141 * t207;
t204 = t146 * qJDD(1);
t255 = -t204 + t196;
t42 = qJD(3) * t197 + (-t132 * t213 - t205) * t145 + t255 * t140;
t64 = t132 * t168;
t89 = t140 * t141 - t145 * t146;
t43 = qJD(1) * t64 + qJDD(1) * t89;
t20 = t139 * t43 + t144 * t42 + t87 * t210 + t211 * t88;
t209 = qJD(5) * t138;
t11 = t138 * t125 + t126 * t208 + t143 * t20 - t171 * t209;
t9 = t11 * t138;
t6 = t267 * t46 + t9;
t217 = -qJD(5) + t53;
t21 = qJD(4) * t171 + t139 * t42 - t144 * t43;
t19 = qJDD(5) + t21;
t48 = t217 * t208;
t231 = t138 * t19 - t48;
t5 = -t171 * t46 + t217 * t259 + t231;
t258 = t138 * t217;
t44 = -t143 * t126 + t138 * t171;
t4 = t143 * t19 + t171 * t44 - t217 * t258;
t12 = qJD(5) * t46 - t143 * t125 + t138 * t20;
t1 = t11 * t143 - t138 * t12 + t258 * t46 - t267 * t44;
t198 = t140 * t210;
t261 = pkin(2) * (qJD(2) * (t139 * t212 + t198) + t139 * t206);
t68 = t144 * t93 - t180;
t65 = -t126 * pkin(4) - t68;
t260 = t65 * t53;
t238 = t171 * t53;
t71 = t144 * t74;
t176 = t211 * t93 - t71;
t256 = -t176 + t262;
t15 = t171 ^ 2 - t53 ^ 2;
t30 = pkin(4) * t171 - t53 * pkin(6);
t13 = -t53 * t126 + t20;
t118 = t139 * pkin(3) + pkin(6);
t245 = t88 * pkin(3);
t28 = -t245 + t30;
t252 = (qJD(5) * t118 + t28) * t217;
t200 = pkin(2) * t214;
t120 = pkin(3) + t240;
t223 = t140 * t144;
t216 = pkin(2) * t223 + t139 * t120;
t82 = pkin(6) + t216;
t251 = (qJD(5) * t82 - t200 + t28) * t217;
t239 = t217 * t171;
t24 = -pkin(4) * t53 - pkin(6) * t171 + t67;
t189 = t125 * pkin(6) + qJD(5) * t24 + t191 + t264;
t61 = -t139 * t168 - t144 * t89;
t63 = t132 * t89;
t26 = -qJD(4) * t61 + t139 * t64 + t144 * t63;
t62 = t139 * t89 - t144 * t168;
t121 = t146 * pkin(2) + pkin(1);
t73 = -t89 * pkin(3) + t121;
t29 = t61 * pkin(4) - t62 * pkin(6) + t73;
t33 = -t125 * pkin(4) + t176 + t261;
t248 = qJD(5) * t217 * t29 - t189 * t61 + t65 * t26 + t33 * t62;
t69 = t139 * t93 + t144 * t202;
t66 = t126 * pkin(6) + t69;
t22 = -t138 * t66 + t143 * t24;
t166 = t250 * t143 - t22 * t171 + t65 * t209;
t23 = t138 * t24 + t143 * t66;
t179 = t33 * t138 + t23 * t171 + t65 * t208;
t14 = t126 * t171 - t21;
t237 = t62 * t65;
t232 = t88 * t87;
t136 = qJDD(1) * pkin(1);
t224 = t139 * t140;
t222 = t142 * t138;
t221 = t142 * t143;
t220 = t144 * t120;
t219 = t147 * t138;
t218 = t147 * t143;
t133 = t141 ^ 2;
t215 = -t146 ^ 2 + t133;
t201 = t141 * t229;
t199 = t62 * t209;
t85 = -t255 * pkin(2) + t136;
t36 = -t43 * pkin(3) + t85;
t3 = t21 * pkin(4) - t20 * pkin(6) + t36;
t193 = qJD(5) * t66 - t3;
t184 = -0.2e1 * pkin(1) * t207;
t182 = qJD(2) * (-qJD(3) + t132);
t181 = qJD(3) * (-qJD(2) - t132);
t178 = qJD(2) * t198;
t170 = t139 * t145 + t223;
t83 = t170 * t229;
t177 = pkin(3) * t211 - t83;
t27 = qJD(4) * t62 + t139 * t63 - t144 * t64;
t55 = -t64 * pkin(3) - t201;
t175 = t29 * t19 - (t27 * pkin(4) - t26 * pkin(6) + t55) * t217;
t173 = g(1) * t142 - g(2) * t147;
t172 = t19 * t62 - t217 * t26;
t169 = t144 * t145 - t224;
t167 = -t189 - t241;
t127 = sin(t135);
t128 = cos(t135);
t165 = g(3) * t128 + t101 * t88 + t174 * t127 + t122;
t163 = -g(3) * t127 - t101 * t87 + t174 * t128;
t149 = qJD(1) ^ 2;
t162 = pkin(1) * t149 + t174;
t56 = t120 * t210 + (qJD(3) * t169 - t140 * t211) * pkin(2);
t161 = -t19 * t82 + t217 * t56 - t260;
t160 = t173 + 0.2e1 * t136;
t84 = t169 * t229;
t152 = -t118 * t19 - t260 - (-pkin(3) * t210 + t84) * t217;
t151 = (-pkin(3) * t126 - t93) * qJD(4) - t156;
t148 = qJD(2) ^ 2;
t119 = -t144 * pkin(3) - pkin(4);
t81 = pkin(2) * t224 - pkin(4) - t220;
t80 = t117 * t218 - t222;
t79 = -t117 * t219 - t221;
t78 = -t117 * t221 - t219;
t77 = t117 * t222 - t218;
t70 = -t200 - t245;
t57 = t120 * t211 + (qJD(3) * t170 + t198) * pkin(2);
t47 = -t87 ^ 2 + t88 ^ 2;
t35 = -t88 * t132 + t43;
t34 = -t87 * t132 + t42;
t2 = t143 * t3;
t7 = [qJDD(1), t173, t174, t133 * qJDD(1) + 0.2e1 * t146 * t196, 0.2e1 * t141 * t204 - 0.2e1 * t207 * t215, -qJDD(2) * t141 - t148 * t146, -qJDD(2) * t146 + t148 * t141, 0, t141 * t184 + t146 * t160, -t141 * t160 + t146 * t184, -t168 * t42 - t88 * t63, -t168 * t43 + t42 * t89 + t63 * t87 - t88 * t64, -t131 * t168 + t63 * t132, t89 * t131 + t64 * t132, 0, -t101 * t64 - t121 * t43 + t128 * t173 + t201 * t87 - t85 * t89, t101 * t63 + t121 * t42 - t127 * t173 - t168 * t85 + t201 * t88, t171 * t26 + t20 * t62, -t171 * t27 - t20 * t61 - t62 * t21 + t26 * t53, t62 * t125 + t26 * t126, -t61 * t125 - t27 * t126, 0, t117 * t173 + t73 * t21 + t67 * t27 + t36 * t61 - t53 * t55, -t116 * t173 + t171 * t55 + t73 * t20 + t67 * t26 + t36 * t62, -t46 * t199 + (t11 * t62 + t26 * t46) * t143, (-t138 * t46 - t143 * t44) * t26 + (-t9 - t12 * t143 + (t138 * t44 - t143 * t46) * qJD(5)) * t62, t11 * t61 + t143 * t172 + t199 * t217 + t46 * t27, -t12 * t61 - t138 * t172 - t44 * t27 + t48 * t62, t19 * t61 - t217 * t27, -g(1) * t78 - g(2) * t80 + t2 * t61 + t22 * t27 + ((-t61 * t66 + t237) * qJD(5) + t175) * t143 + t248 * t138, -g(1) * t77 - g(2) * t79 - t23 * t27 + t248 * t143 + (-qJD(5) * t237 + t193 * t61 - t175) * t138; 0, 0, 0, -t141 * t149 * t146, t215 * t149, -t205, -t204, qJDD(2), g(3) * t146 + t141 * t162, -g(3) * t141 + t146 * t162, t232, t47, t34, t35, t131, (t131 * t145 + t140 * t181 - t214 * t87) * pkin(2) + t165, (-t88 * t214 + (-qJDD(2) - t131) * t140 + t145 * t181) * pkin(2) + t163, -t238, t15, t13, t14, t125, t125 * t220 - t57 * t126 + t70 * t53 + (-t178 + (-t195 + (-qJDD(2) - t125) * t140) * t139) * pkin(2) + t256, -t125 * t216 - t56 * t126 - t171 * t70 + t263, t6, t1, t5, t4, t239, t81 * t12 + t57 * t44 + (-t33 + t251) * t143 + t161 * t138 + t166, t81 * t11 + t57 * t46 + t161 * t143 + (-t250 - t251) * t138 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t47, t34, t35, t131, pkin(2) * t140 * t182 + t165, (t145 * t182 - t206) * pkin(2) + t163, -t238, t15, t13, t14, t125, -pkin(2) * t178 + t83 * t126 + t71 + (t125 * t144 - t53 * t88) * pkin(3) + t151 * t139 + t262, t84 * t126 + (-t125 * t139 + t171 * t88) * pkin(3) + t151 * t144 + t265, t6, t1, t5, t4, t239, t119 * t12 + t177 * t44 + (-t33 + t252) * t143 + t152 * t138 + t166, t119 * t11 + t177 * t46 + t152 * t143 + (-t250 - t252) * t138 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t15, t13, t14, t125, t69 * t126 + t256 - t261, t68 * t126 + t263, t6, t1, t5, t4, t239, -pkin(4) * t12 - t33 * t143 + (-t138 * t68 + t143 * t30) * t217 - t69 * t44 - t138 * t260 - t231 * pkin(6) + t166, -pkin(4) * t11 - t69 * t46 + (-pkin(6) * t19 - t217 * t68 - t260) * t143 + (-(pkin(6) * qJD(5) + t30) * t217 - t250) * t138 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t44, -t44 ^ 2 + t46 ^ 2, -t217 * t44 + t11, -t217 * t46 - t12, t19, -g(1) * t79 + g(2) * t77 + t138 * t167 - t208 * t66 - t217 * t23 - t65 * t46 + t2, g(1) * t80 - g(2) * t78 + t138 * t193 + t143 * t167 - t217 * t22 + t65 * t44;];
tau_reg = t7;
