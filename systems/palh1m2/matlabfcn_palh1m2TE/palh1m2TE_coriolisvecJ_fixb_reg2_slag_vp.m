% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh1m2TE_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:47:08
% EndTime: 2020-05-01 20:47:45
% DurationCPUTime: 6.03s
% Computational Cost: add. (10211->276), mult. (20254->522), div. (0->0), fcn. (25215->22), ass. (0->243)
t166 = sin(qJ(3));
t292 = pkin(1) * qJD(2);
t237 = t166 * t292;
t157 = qJD(2) + qJD(3);
t297 = t157 * pkin(5);
t192 = t237 + t297;
t242 = pkin(22) + pkin(21);
t150 = cos(t242);
t218 = sin(t242);
t167 = sin(qJ(2));
t274 = sin(pkin(20));
t277 = cos(pkin(20));
t303 = sin(pkin(18));
t305 = cos(pkin(18));
t127 = t303 * t274 + t277 * t305;
t169 = cos(qJ(3));
t191 = -t305 * t274 + t277 * t303;
t330 = t127 * t169 + t166 * t191;
t170 = cos(qJ(2));
t89 = -t166 * t127 + t169 * t191;
t72 = t89 * t170;
t227 = -t330 * t167 + t72;
t314 = t330 * t170;
t339 = t89 * t167 + t314;
t342 = t339 * t150 + t218 * t227;
t340 = t342 * t192;
t240 = t169 * t292;
t335 = t150 * t227 - t218 * t339;
t341 = t192 * t335;
t18 = t342 * t240 - t341;
t164 = cos(pkin(19));
t275 = sin(pkin(19));
t125 = t169 * t164 - t166 * t275;
t116 = t125 * qJD(3);
t245 = qJD(1) * qJD(2);
t231 = t167 * t245;
t251 = qJD(1) * t170;
t123 = t166 * t164 + t169 * t275;
t115 = t123 * qJD(3);
t271 = t115 * t167;
t324 = pkin(2) * t157;
t333 = t123 * t324;
t287 = t167 * t333;
t299 = pkin(2) * t123;
t97 = t125 * t324;
t338 = ((t271 + (-qJD(2) * t125 - t116) * t170) * pkin(2) - t287) * qJD(1) + t231 * t299 + t97 * t251;
t300 = pkin(2) * t116;
t88 = t123 * t170 + t125 * t167;
t318 = t88 * qJD(1);
t334 = pkin(2) * t318;
t91 = -t123 * t167 + t125 * t170;
t188 = t91 * qJD(1);
t272 = qJD(1) * pkin(15);
t61 = -pkin(2) * t188 - t272;
t66 = t115 * t324;
t111 = qJD(3) * t300;
t112 = qJD(2) * t300;
t67 = t112 + t111;
t336 = -t334 * t61 + t67 * t299 + t333 * t300 + (-t115 * t97 - t125 * t66) * pkin(2);
t320 = t188 * t318;
t332 = t318 * t251;
t331 = t318 * t272;
t65 = -t127 * t150 - t191 * t218;
t171 = cos(pkin(17));
t304 = sin(pkin(17));
t129 = t303 * t171 - t304 * t305;
t131 = t171 * t305 + t304 * t303;
t220 = t167 * t129 + t131 * t170;
t85 = qJD(2) * t220;
t329 = -t188 ^ 2 + t318 ^ 2;
t182 = -qJD(2) * t88 - t116 * t167;
t41 = qJD(1) * t182 - t115 * t251;
t328 = t157 * t318 + t41;
t325 = 0.2e1 * qJD(2);
t252 = qJD(1) * t167;
t120 = -t166 * t252 + t169 * t251;
t130 = -t166 * t167 + t169 * t170;
t103 = t157 * t130;
t84 = qJD(1) * t103;
t50 = -t157 * t120 + t84;
t321 = t188 * t272;
t253 = qJD(1) * t157;
t316 = -qJD(3) + t157;
t183 = qJD(2) * t91 - t271;
t42 = qJD(1) * t183 + t116 * t251;
t315 = -t157 * t188 + t42;
t165 = sin(qJ(4));
t168 = cos(qJ(4));
t147 = pkin(1) * t167 - pkin(15);
t138 = t147 * qJD(1);
t104 = -t120 * pkin(5) + t138;
t64 = t218 * t127 - t191 * t150;
t62 = t64 * qJD(1);
t63 = t65 * qJD(1);
t39 = -t63 * pkin(9) - t62 * pkin(11) + t104;
t11 = -t165 * t18 + t168 * t39;
t249 = qJD(2) * t170;
t250 = qJD(2) * t167;
t76 = t330 * qJD(3);
t281 = t76 * t167;
t82 = t89 * qJD(3);
t280 = t82 * t167;
t37 = -qJD(2) * t339 - t76 * t170 - t280;
t25 = (-qJD(3) * t72 - t89 * t249 + t250 * t330 + t281) * t150 - t218 * t37;
t185 = -t169 * t25 + (-t166 * t342 - t169 * t335) * qJD(3);
t202 = t218 * (qJD(2) * t227 + t82 * t170 - t281);
t24 = (qJD(3) * t314 + t249 * t330 + t250 * t89 + t280) * t150 + t202;
t4 = t24 * t297 + (t166 * t24 + t185) * t292;
t230 = t170 * t245;
t217 = pkin(1) * t230;
t128 = t166 * t170 + t169 * t167;
t102 = t157 * t128;
t83 = qJD(1) * t102;
t57 = t83 * pkin(5) + t217;
t1 = t11 * qJD(4) + t165 * t57 + t168 * t4;
t12 = t165 * t39 + t168 * t18;
t2 = -t12 * qJD(4) - t165 * t4 + t168 * t57;
t204 = t11 * t165 - t12 * t168;
t313 = qJD(4) * t204 - t1 * t165 - t168 * t2;
t312 = 0.2e1 * qJD(1);
t186 = pkin(1) * t335;
t181 = qJD(2) * t186;
t16 = t169 * t181 + t340;
t178 = t166 * t181;
t9 = -t178 + t341;
t309 = t16 * t9;
t308 = t9 * t62;
t149 = pkin(1) * t166 + pkin(5);
t23 = t37 * t150 - t202;
t302 = pkin(1) * t169;
t238 = qJD(3) * t302;
t248 = qJD(3) * t166;
t307 = -t25 * t149 - t186 * t248 + t23 * t302 + t238 * t342 + t18;
t17 = t240 * t335 + t340;
t6 = pkin(1) * t185 + t24 * t149;
t306 = -t17 + t6;
t301 = pkin(1) * t170;
t119 = qJD(1) * t128;
t298 = t119 * pkin(5);
t296 = t25 * t62;
t239 = pkin(1) * t249;
t68 = t102 * pkin(5) + t239;
t295 = t68 * t62;
t162 = qJ(3) + qJ(2);
t153 = sin(t162);
t273 = sin(pkin(22));
t276 = cos(pkin(22));
t124 = -t303 * t273 - t276 * t305;
t126 = t273 * t305 - t303 * t276;
t58 = (-cos(pkin(21)) * t124 - t126 * sin(pkin(21))) * pkin(4) + t147;
t291 = t153 * t58;
t290 = t157 * t58;
t59 = -qJD(4) + t63;
t289 = t165 * t59;
t288 = t165 * t64;
t286 = t168 * t59;
t285 = t168 * t62;
t284 = t168 * t64;
t173 = qJD(1) ^ 2;
t98 = t170 * t129 - t167 * t131;
t283 = t173 * t98;
t282 = t173 * t220;
t86 = qJD(2) * t98;
t270 = t119 * t120;
t154 = cos(t162);
t266 = t154 * t173;
t265 = t157 * t119;
t263 = t165 * t168;
t260 = t170 * t173;
t172 = qJD(2) ^ 2;
t259 = t172 * t167;
t258 = t172 * t170;
t141 = qJD(2) * t238;
t257 = t157 * t238 + t141;
t256 = t153 ^ 2 - t154 ^ 2;
t158 = t165 ^ 2;
t160 = t168 ^ 2;
t255 = t158 - t160;
t254 = t167 ^ 2 - t170 ^ 2;
t247 = qJD(4) * t165;
t246 = qJD(4) * t168;
t244 = t85 * t312;
t243 = t86 * t312;
t241 = pkin(1) * t260;
t236 = qJD(4) * t62 * t64;
t60 = t62 ^ 2;
t235 = t60 * t263;
t234 = t58 * t266;
t233 = t307 * t62;
t232 = t120 * t138;
t229 = qJD(2) * t248;
t226 = (-qJD(4) - t59) * t62;
t43 = -t115 * t170 + t182;
t224 = qJD(1) * t43 + t41;
t70 = -t91 * pkin(2) - pkin(15);
t223 = qJD(1) * t70 + t61;
t219 = t147 * t241;
t216 = t153 * t154 * t253;
t215 = -t130 * pkin(5) + t147;
t214 = t167 * t230;
t3 = -qJD(3) * t178 + t141 * t342 - t192 * t25 + t23 * t240;
t213 = -t16 * t25 + t3 * t342;
t212 = 0.2e1 * t217;
t20 = t149 * t342 + t169 * t186;
t21 = -t149 * t335 + t302 * t342;
t210 = t20 * t62 + t21 * t59;
t209 = t24 * t59 - t296;
t208 = -t335 * t59 + t342 * t62;
t207 = pkin(2) * t224;
t206 = t236 * t263;
t205 = t11 * t168 + t12 * t165;
t198 = t124 * t170 + t126 * t167;
t80 = t198 * qJD(2);
t93 = t167 * t124 - t126 * t170;
t81 = t93 * qJD(2);
t203 = t124 * t81 + t126 * t80;
t196 = -t157 * t240 + t141;
t195 = t316 * t237;
t194 = t307 * t16 + t3 * t20;
t193 = t59 * t6 + t233;
t190 = qJD(4) * (t59 * t64 + t62 * t65);
t177 = -t205 * qJD(4) + t1 * t168 - t165 * t2;
t156 = t157 ^ 2;
t140 = t167 * t260;
t136 = pkin(1) * t157 * t248;
t135 = -0.2e1 * t214;
t134 = 0.2e1 * t214;
t133 = t254 * t173;
t132 = t153 * t266;
t122 = t147 * t212;
t121 = t256 * t173;
t118 = 0.2e1 * t254 * t245;
t105 = pkin(1) * t251 + t298;
t55 = t119 ^ 2 - t120 ^ 2;
t49 = -t128 * t253 + t265;
t44 = t116 * t170 + t183;
t40 = -t65 * pkin(9) - t64 * pkin(11) + t215;
t14 = t165 * t105 + t168 * t17;
t13 = t168 * t105 - t165 * t17;
t10 = -t237 * t342 + t340;
t8 = t168 * t10 + t165 * t298;
t7 = -t165 * t10 + t168 * t298;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t118, -t259, t134, -t258, 0, -0.2e1 * pkin(15) * t230, 0.2e1 * pkin(15) * t231, 0, 0, t119 * t103 + t84 * t128, -t119 * t102 + t103 * t120 - t128 * t83 + t84 * t130, t157 * t103, -t120 * t102 - t83 * t130, -t157 * t102, 0, 0.2e1 * t147 * t83 + (-qJD(1) * t130 - t120) * t239, 0.2e1 * t119 * t239 + 0.2e1 * t147 * t84, (t102 * t169 - t103 * t166 + (-t128 * t169 + t130 * t166) * qJD(3)) * t292, t122, 0, 0, 0, 0, 0, 0, -t57 * t65 - t68 * t63, t57 * t64 + t295, t3 * t64 + t4 * t65, t104 * t68 + t215 * t57, -0.2e1 * t206, 0.2e1 * t255 * t236, t165 * t190, 0.2e1 * t206, t168 * t190, 0, t3 * t288 - t68 * t286 - t2 * t65 + (t16 * t284 + t40 * t289) * qJD(4), t68 * t289 + t3 * t284 + t1 * t65 + (-t16 * t288 + t40 * t286) * qJD(4), (-t158 - t160) * t295 + t313 * t64, t205 * t68 - t313 * t40, t220 * t243, (-t220 * t85 + t86 * t98) * t312, qJD(2) * t86, -t98 * t244, -qJD(2) * t85, 0, pkin(14) * t244, pkin(14) * t243, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t124 * t217, t126 * t212, t203 * t292, t122, t318 * t44 + t42 * t88, t188 * t44 + t318 * t43 + t41 * t88 + t91 * t42, t157 * t44, t188 * t43 + t41 * t91, t43 * t157, 0, t224 * pkin(15), (-qJD(1) * t44 - t42) * pkin(15), 0, 0, t135, t118, -t259, t134, -t258, 0, -t167 * t207 + t223 * t249, -t170 * t207 - t223 * t250, t66 * t167 - t67 * t170 + (-t170 * t97 + t287) * qJD(2), (-t41 * t70 - t43 * t61) * pkin(2), 0.2e1 * t216, -0.2e1 * t256 * t253, t156 * t154, -0.2e1 * t216, -t156 * t153, 0, (t153 * t290 - t154 * t239) * t312, (t153 * t239 + t154 * t290) * t312, t316 * (t153 * t169 - t154 * t166) * t292, t58 * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t133, 0, -t140, 0, 0, pkin(15) * t260, -t173 * pkin(15) * t167, 0, 0, -t270, t55, t50, t270, t49, 0, (-t119 * t147 + t120 * t301) * qJD(1) + t257, -t232 - t136 + (-t119 * t251 - t229) * pkin(1), ((t83 - t265) * t169 - t50 * t166) * pkin(1), -t219, 0, 0, 0, 0, 0, 0, t105 * t63, -t105 * t62, t306 * t63 + t233, -t104 * t105 + t306 * t18 + t4 * t21 + t194, 0, 0, 0, 0, 0, 0, t13 * t59 + t165 * t193 + t210 * t246, -t14 * t59 + t168 * t193 - t210 * t247, (t13 * t168 + t14 * t165) * t62, (t168 * t6 - t14) * t12 + (-t165 * t6 - t13) * t11 + t177 * t21 + t194, -t220 * t283, (t220 ^ 2 - t98 ^ 2) * t173, 0, t98 * t282, 0, 0, -pkin(14) * t282, -pkin(14) * t283, 0, 0, 0, 0, 0, 0, 0, 0, t124 * t241, -t126 * t241, ((-t124 * t93 - t126 * t198) * qJD(2) + t203) * pkin(1) * qJD(1), -t219 + (-t198 * t81 + t80 * t93) * pkin(1) ^ 2 * t325, -t320, t329, t315, t320, t328, 0, t331, t321, 0, 0, t140, -t133, 0, -t140, 0, 0, t111 + 0.2e1 * t112 + (-t167 * t334 - t170 * t61) * qJD(1), t61 * t252 + (-t332 + (t325 + qJD(3)) * t115) * pkin(2), t338, t336, -t132, t121, 0, t132, 0, 0, (t154 * t301 - t291) * t173 + t257, -t234 - t136 + (-t153 * t260 - t229) * pkin(1), 0, -t58 * t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, t55, t50, t270, t49, 0, -t119 * t138 + t196, t195 - t232, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t298, -t62 * t298, -t10 * t63 - t308 + (t24 * t63 - t296) * pkin(5), -t18 * t10 - t309 + (-t104 * t119 + t18 * t24 - t335 * t4 + t213) * pkin(5), 0, 0, 0, 0, 0, 0, -t165 * t308 + t7 * t59 + (t165 * t209 + t208 * t246) * pkin(5), -t9 * t285 - t8 * t59 + (t168 * t209 - t208 * t247) * pkin(5), (t165 * t8 + t168 * t7) * t62, -t11 * t7 - t12 * t8 - t309 + (-t177 * t335 - t204 * t24 + t213) * pkin(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t320, t329, t315, t320, t328, 0, t331, t321, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * qJD(2) - t252 * t334 + t112, -t333 * qJD(2) + (qJD(2) * t115 - t332) * pkin(2), t338, t336, -t132, t121, 0, t132, 0, 0, -t173 * t291 + t196, t195 - t234, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, -t255 * t60, t165 * t226, -t235, t168 * t226, 0, -t12 * t59 - t16 * t285 + t2, -t11 * t59 + (-qJD(4) * t39 - t4) * t168 + (qJD(4) * t18 + t16 * t62 - t57) * t165, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg = t5;
