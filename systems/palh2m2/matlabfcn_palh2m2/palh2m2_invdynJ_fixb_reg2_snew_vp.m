% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% palh2m2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = palh2m2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 14:45:30
% EndTime: 2019-06-06 14:45:33
% DurationCPUTime: 3.07s
% Computational Cost: add. (24435->252), mult. (61313->367), div. (0->0), fcn. (35203->8), ass. (0->176)
t196 = sin(qJ(2));
t183 = t196 * qJDD(1);
t200 = cos(qJ(2));
t230 = qJD(1) * qJD(2);
t222 = t200 * t230;
t176 = t183 + t222;
t184 = t200 * qJDD(1);
t223 = t196 * t230;
t177 = t184 - t223;
t149 = t176 * t196 + t177 * t200;
t199 = cos(qJ(3));
t144 = t199 * t149;
t195 = sin(qJ(3));
t190 = t196 ^ 2;
t193 = t200 ^ 2;
t260 = t190 + t193;
t180 = t260 * qJD(1);
t136 = -t180 * qJD(2) + t176 * t200 - t196 * t177;
t231 = qJD(3) * t180;
t214 = t136 + t231;
t124 = -t195 * t214 + t144;
t264 = -t195 * t231 + t124;
t263 = t177 - t223;
t189 = t195 ^ 2;
t192 = t199 ^ 2;
t234 = t189 + t192;
t262 = pkin(3) * t189 + t199 * (pkin(3) * t199 + pkin(5)) + pkin(2);
t203 = qJD(1) ^ 2;
t197 = sin(qJ(1));
t257 = cos(qJ(1));
t205 = t257 * g(1) + t197 * g(2);
t173 = -t203 * pkin(1) - t205;
t158 = g(3) * t200 + t173 * t196;
t181 = t200 * t203 * t196;
t229 = qJDD(2) + t181;
t147 = t229 * pkin(4) - t158;
t159 = -g(3) * t196 + t173 * t200;
t186 = t193 * t203;
t202 = qJD(2) ^ 2;
t151 = (-t202 - t186) * pkin(4) + t159;
t132 = -t196 * t147 + t151 * t200;
t179 = t180 ^ 2;
t127 = -pkin(2) * t179 + t132;
t244 = t147 * t200;
t131 = t151 * t196 + t244;
t116 = t195 * t127 - t199 * t131;
t166 = t199 * t179 * t195;
t160 = qJDD(3) + t166;
t109 = t160 * pkin(5) - t116;
t117 = t199 * t127 + t195 * t131;
t201 = qJD(3) ^ 2;
t239 = t179 * t192;
t165 = -t201 - t239;
t110 = pkin(5) * t165 + t117;
t84 = t199 * t109 + t195 * t110;
t81 = t199 * t84;
t85 = -t195 * t109 + t110 * t199;
t49 = -t195 * t85 + t81;
t261 = t195 * t49;
t156 = t234 * t180;
t152 = qJD(4) + t156;
t258 = t152 ^ 2;
t153 = t156 ^ 2;
t256 = pkin(4) * t153;
t255 = pkin(4) * t179;
t254 = pkin(5) * t153;
t80 = t195 * t84;
t51 = t199 * t85 + t80;
t26 = t195 * t51 + t199 * t49;
t48 = pkin(5) * t49;
t253 = pkin(2) * t26 + t48;
t243 = t149 * t195;
t126 = t199 * t214 + t243;
t105 = t124 * t199 + t126 * t195;
t245 = t105 * t195;
t100 = -t156 * qJD(3) - t195 * t124 + t126 * t199;
t246 = t100 * t199;
t68 = -t245 - t246;
t102 = t199 * t105;
t247 = t100 * t195;
t70 = t102 - t247;
t41 = t195 * t70 + t199 * t68;
t66 = pkin(5) * t68;
t252 = pkin(2) * t41 + t66;
t194 = sin(qJ(4));
t198 = cos(qJ(4));
t75 = -pkin(3) * t153 + t85;
t221 = g(1) * t197 - t257 * g(2);
t172 = qJDD(1) * pkin(1) + t221;
t145 = t263 * pkin(4) + t172;
t133 = pkin(2) * t149 + t145;
t108 = t264 * pkin(5) + t133;
t78 = -pkin(3) * t105 - t108;
t44 = t194 * t75 + t198 * t78;
t45 = -t194 * t78 + t198 * t75;
t22 = t194 * t44 + t198 * t45;
t15 = -t195 * t22 + t81;
t83 = pkin(3) * t84;
t251 = pkin(5) * t15 + t83;
t103 = qJDD(4) + t105;
t216 = -t103 * t194 - t198 * t258;
t249 = t195 * t216;
t248 = t198 * t84;
t242 = t160 * t195;
t161 = qJDD(3) - t166;
t241 = t161 * t199;
t240 = t179 * t189;
t238 = t194 * t195;
t237 = t194 * t199;
t236 = t198 * t100;
t235 = t198 * t199;
t232 = pkin(3) * t100;
t16 = t199 * t22 + t80;
t6 = t15 * t199 + t16 * t195;
t228 = pkin(2) * t6 + t251;
t227 = t100 * t238;
t226 = t100 * t237;
t220 = -t194 * t84 + t198 * t232;
t219 = t194 * t232 + t248;
t218 = t200 * pkin(4) + pkin(1);
t215 = t116 * t195 + t199 * t117;
t213 = t136 + 0.2e1 * t231;
t206 = -t198 * t103 + t258 * t194;
t97 = t100 * t235;
t59 = -t195 * t206 + t97;
t212 = pkin(5) * t59 + t220;
t56 = t226 - t249;
t211 = pkin(5) * t56 + t219;
t138 = t199 * t160 + t165 * t195;
t210 = pkin(2) * t138 - t116;
t20 = t194 * t45 - t198 * t44;
t95 = t195 * t236;
t63 = t199 * t206 + t95;
t36 = t195 * t63 + t199 * t59;
t209 = pkin(2) * t36 + t212;
t90 = t199 * t216;
t60 = t90 + t227;
t33 = t195 * t60 + t199 * t56;
t208 = pkin(2) * t33 + t211;
t207 = t234 * t260;
t74 = t195 * t254 + t84;
t163 = -t201 - t240;
t140 = -t161 * t195 + t163 * t199;
t204 = pkin(2) * t140 - t117;
t191 = t198 ^ 2;
t188 = t194 ^ 2;
t185 = t190 * t203;
t178 = t184 - 0.2e1 * t223;
t175 = t183 + 0.2e1 * t222;
t155 = (t189 - t192) * t179;
t154 = t234 * t179;
t129 = t136 * t199 + t243;
t128 = -t136 * t195 + t144;
t125 = t199 * t213 + t243;
t123 = t195 * t213 - t144;
t114 = -pkin(2) * t123 + t133 * t199;
t113 = -pkin(2) * t125 - t133 * t195;
t112 = t128 * t195 - t129 * t199;
t111 = pkin(2) * t112;
t96 = t194 * t236;
t88 = -t116 * t199 + t117 * t195;
t87 = pkin(2) * t88;
t86 = pkin(2) * t154 + t215;
t79 = (pkin(5) * t199 + pkin(2)) * t108;
t65 = (t188 - t191) * t100;
t64 = (t188 + t191) * t100;
t46 = pkin(2) * t105 + t108 * t189 + t199 * (pkin(5) * t105 + t108 * t199);
t43 = -pkin(2) * t100 - pkin(5) * t246;
t30 = pkin(3) * t216 - t45;
t29 = -pkin(3) * t206 - t44;
t23 = pkin(2) * t153 - t261 + t199 * (t51 + t254);
t8 = pkin(2) * t216 + (pkin(5) * t216 - t195 * t248 + t199 * t30) * t199 + t195 * (t195 * t30 + t84 * t235);
t7 = -pkin(2) * t206 + (-pkin(5) * t206 + t199 * t29 - t84 * t238) * t199 + t195 * (t195 * t29 + t84 * t237);
t2 = t262 * t64;
t1 = t262 * t20;
t3 = [0, 0, 0, 0, 0, qJDD(1), t221, t205, 0, 0, (t176 + t222) * t196, t175 * t200 + t178 * t196, t196 * t229 + t200 * (-t185 + t202), t263 * t200, t196 * (t186 - t202) + t200 * (qJDD(2) - t181), 0, pkin(1) * t178 + t172 * t200, -pkin(1) * t175 - t172 * t196, t196 * t158 + t200 * t159 + pkin(1) * (t185 + t186), pkin(1) * t172, 0, t260 * t136, 0, t260 * t149, 0, 0, t190 * t145 + t200 * (pkin(4) * t149 + t145 * t200) + pkin(1) * t149, -t218 * t136, t190 * t132 + t200 * (t132 * t200 + t255) + pkin(1) * t179, t218 * t145, t260 * (t199 * t231 + t126) * t195, t260 * (-t123 * t195 + t125 * t199), t260 * (t242 + (t201 - t240) * t199), t260 * t264 * t199, t260 * (t241 + (-t201 + t239) * t195), 0, t190 * t114 + t200 * (-pkin(4) * t123 + t114 * t200) - pkin(1) * t123, t190 * t113 + t200 * (-pkin(4) * t125 + t113 * t200) - pkin(1) * t125, t200 * (pkin(4) * t154 + t200 * t86) + pkin(1) * t154 + t86 * t190, (t190 * pkin(2) + t200 * (pkin(2) * t200 + pkin(4)) + pkin(1)) * t133, 0, t260 * (t195 * (t102 + t247) + t199 * (-t245 + t246)), 0, t105 * t207, 0, 0, t200 * (pkin(4) * t105 + t200 * t46) + pkin(1) * t105 + t46 * t190, t200 * (-pkin(4) * t100 + t200 * t43) - pkin(1) * t100 + t43 * t190, t200 * (t200 * t23 + t256) + pkin(1) * t153 + t23 * t190, t190 * t79 + t200 * (pkin(4) * t108 + t200 * t79) + pkin(1) * t108, 0, 0, t260 * (t95 * t195 + t199 * t97), 0, t260 * (t195 * (t90 - t227) + t199 * (-t226 - t249)), t103 * t207, t190 * t7 + t200 * (-pkin(4) * t206 + t200 * t7) - pkin(1) * t206, t190 * t8 + t200 * (pkin(4) * t216 + t200 * t8) + pkin(1) * t216, t200 * (-pkin(4) * t64 - t2 * t200) - pkin(1) * t64 - t2 * t190, t200 * (pkin(4) * t20 + t1 * t200) + pkin(1) * t20 + t1 * t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t186 + t185, t183, t181, t184, qJDD(2), -t158, -t159, 0, 0, 0, -t179, t136, 0, t149, 0, t244 + (t151 + t255) * t196, -t132, pkin(4) * (-t136 * t200 - t149 * t196), pkin(4) * (t131 * t200 - t132 * t196), -t166, t155, t129, t166, t128, qJDD(3), pkin(4) * (-t196 * (t165 * t199 - t242) + t200 * t138) + t210, pkin(4) * (-t196 * (-t163 * t195 - t241) + t200 * t140) + t204, pkin(4) * (-t196 * (t128 * t199 + t129 * t195) + t200 * t112) + t111, pkin(4) * (-t196 * t215 + t200 * t88) + t87, 0, -t153, t100, 0, t105, 0, t196 * t234 * t256 + t74, -t85, pkin(4) * (-t196 * (-t195 * t68 + t199 * t70) + t200 * t41) + t252, pkin(4) * (-t196 * (t199 * t51 - t261) + t200 * t26) + t253, -t96, t65, t216, t96, t206, 0, pkin(4) * (-t196 * (-t195 * t56 + t199 * t60) + t200 * t33) + t208, pkin(4) * (-t196 * (-t195 * t59 + t199 * t63) + t200 * t36) + t209, -t22, pkin(4) * (-t196 * (-t15 * t195 + t16 * t199) + t200 * t6) + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t136, 0, t149, 0, t131, -t132, 0, 0, -t166, t155, t129, t166, t128, qJDD(3), t210, t204, t111, t87, 0, -t153, t100, 0, t105, 0, t74, -t85, t252, t253, -t96, t65, t216, t96, t206, 0, t208, t209, -t22, t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t155, t129, t166, t128, qJDD(3), -t116, -t117, 0, 0, 0, -t153, t100, 0, t105, 0, t74, -t85, t66, t48, -t96, t65, t216, t96, t206, 0, t211, t212, -t22, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t100, 0, t105, 0, t84, -t85, 0, 0, -t96, t65, t216, t96, t206, 0, t219, t220, -t22, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, 0, -t194 * t100, t103, -t44, -t45, 0, 0;];
tauJ_reg  = t3;
