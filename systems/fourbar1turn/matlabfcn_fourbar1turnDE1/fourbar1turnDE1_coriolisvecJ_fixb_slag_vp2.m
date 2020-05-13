% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbar1turnDE1
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
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:32
% EndTime: 2020-04-12 19:26:05
% DurationCPUTime: 12.10s
% Computational Cost: add. (141585->381), mult. (203602->841), div. (7454->21), fcn. (54637->8), ass. (0->305)
t135 = pkin(2) ^ 2;
t136 = pkin(1) ^ 2;
t125 = cos(qJ(2));
t325 = pkin(2) * t125;
t273 = -0.2e1 * pkin(1) * t325 + t136;
t118 = t135 + t273;
t353 = pkin(3) ^ 2;
t354 = pkin(4) ^ 2;
t272 = t353 - t354;
t113 = t118 + t272;
t120 = pkin(1) * t125 - pkin(2);
t124 = sin(qJ(2));
t338 = -pkin(3) - pkin(4);
t111 = (pkin(2) - t338) * (pkin(2) + t338) + t273;
t337 = pkin(4) - pkin(3);
t112 = (pkin(2) - t337) * (pkin(2) + t337) + t273;
t287 = t111 * t112;
t139 = sqrt(-t287);
t280 = t124 * t139;
t91 = -pkin(1) * t280 - t113 * t120;
t83 = 0.1e1 / t91 ^ 2;
t331 = pkin(1) * t113;
t108 = t124 * t331;
t94 = -t120 * t139 + t108;
t309 = t83 * t94;
t115 = 0.1e1 / t118;
t132 = 0.1e1 / pkin(3);
t116 = 0.1e1 / t118 ^ 2;
t266 = qJD(2) * t124;
t244 = pkin(2) * t266;
t220 = pkin(1) * t244;
t194 = t116 * t220;
t182 = -0.2e1 * t194;
t342 = -0.2e1 * t120;
t261 = pkin(2) * t342;
t224 = t113 + t261;
t103 = 0.1e1 / t139;
t214 = pkin(1) * pkin(2) * (-t111 - t112);
t99 = t124 * t214;
t98 = qJD(2) * t99;
t298 = t103 * t98;
t240 = t124 * t298;
t278 = t125 * t139;
t59 = (-t240 + (t124 * t224 - t278) * qJD(2)) * pkin(1);
t51 = (t115 * t59 + t182 * t91) * t132;
t123 = t124 ^ 2;
t283 = t123 * t136;
t246 = pkin(2) * t283;
t215 = qJD(2) * t246;
t232 = t139 * t266;
t265 = qJD(2) * t125;
t274 = pkin(1) * t232 + t265 * t331;
t289 = t103 * t120;
t62 = -t289 * t98 + 0.2e1 * t215 + t274;
t53 = (t115 * t62 + t182 * t94) * t132;
t82 = 0.1e1 / t91;
t348 = -t51 * t309 + t53 * t82;
t90 = t94 ^ 2;
t80 = t83 * t90 + 0.1e1;
t77 = 0.1e1 / t80;
t161 = t348 * t77;
t114 = t118 - t272;
t119 = pkin(1) - t325;
t92 = -pkin(2) * t280 + t114 * t119;
t87 = 0.1e1 / t92 ^ 2;
t326 = pkin(2) * t124;
t107 = t114 * t326;
t93 = t119 * t139 + t107;
t307 = t87 * t93;
t129 = 0.1e1 / pkin(4);
t335 = -t115 / 0.2e1;
t259 = 0.2e1 * t119 * pkin(1);
t60 = (-t240 + (-t278 + (t114 + t259) * t124) * qJD(2)) * pkin(2);
t50 = (t194 * t92 + t335 * t60) * t129;
t334 = t115 / 0.2e1;
t247 = t135 * t123 * pkin(1);
t216 = qJD(2) * t247;
t275 = (t114 * t265 + t232) * pkin(2);
t290 = t103 * t119;
t61 = t290 * t98 + 0.2e1 * t216 + t275;
t52 = (-t194 * t93 + t334 * t61) * t129;
t89 = t93 ^ 2;
t79 = t87 * t89 + 0.1e1;
t75 = 0.1e1 / t79;
t86 = 0.1e1 / t92;
t163 = t75 * (t50 * t307 + t52 * t86);
t254 = 0.2e1 * t116;
t269 = qJD(1) * t129;
t352 = -0.4e1 * t92;
t281 = t124 * t135;
t248 = pkin(1) * t281;
t133 = 0.1e1 / t353;
t81 = t91 ^ 2;
t304 = t81 + t90;
t74 = t304 * t133 * t116;
t69 = t74 ^ (-0.1e1 / 0.2e1);
t181 = t69 * t248 * t254;
t323 = pkin(4) * t118;
t260 = 0.2e1 * t323;
t263 = qJD(1) * qJD(2);
t228 = t124 * t263;
t351 = pkin(2) * t228;
t350 = t98 * t99;
t252 = t75 * t323;
t221 = t87 * t252;
t349 = t129 * t93 * t221;
t264 = qJD(2) * t132;
t130 = 0.1e1 / t354;
t85 = t92 ^ 2;
t303 = t85 + t89;
t73 = t303 * t130 * t116;
t67 = t73 ^ (-0.1e1 / 0.2e1);
t180 = t67 * t194;
t285 = t115 * t129;
t117 = t115 * t116;
t255 = pkin(1) * t326;
t212 = 0.4e1 * t255;
t185 = t117 * t212;
t157 = t303 * t185;
t38 = ((t60 * t92 + t61 * t93) * t254 - qJD(2) * t157) * t130;
t68 = t67 / t73;
t347 = t38 * t68 * t285 / 0.4e1 + t129 * t180;
t71 = 0.1e1 / t74;
t346 = -0.2e1 * t87;
t345 = 0.4e1 * t93;
t344 = 0.2e1 * t94;
t343 = -0.2e1 * t115;
t127 = qJD(2) ^ 2;
t223 = t116 * t255;
t200 = t92 * t223;
t297 = t103 * t99;
t64 = t107 + (-t278 + (t259 - t297) * t124) * pkin(2);
t54 = (t335 * t64 + t200) * t129;
t199 = t93 * t223;
t66 = t99 * t290 + 0.2e1 * t247 + (t114 * t125 + t280) * pkin(2);
t56 = (t334 * t66 - t199) * t129;
t162 = t75 * (-t307 * t54 - t56 * t86);
t76 = 0.1e1 / t79 ^ 2;
t230 = 0.4e1 * t76 * t307;
t310 = t76 * t86;
t88 = t86 / t85;
t306 = t88 * t89;
t42 = -t306 * t60 + t307 * t61;
t251 = t42 * t310;
t256 = t88 * t345;
t258 = 0.2e1 * t75 * t87;
t288 = 0.4e1 * t103 / t287;
t231 = -t288 / 0.4e1;
t188 = t231 * t350;
t234 = t135 * t283;
t151 = t125 * t214 - 0.4e1 * t234;
t96 = t151 * qJD(2);
t148 = (-0.2e1 * t125 * t298 + (-t103 * t96 + t188) * t124) * t115;
t209 = t117 * t234;
t179 = t209 * t345;
t184 = qJD(2) * t209;
t195 = t129 * t86 * t252;
t198 = 0.6e1 * t125 * t248;
t202 = t119 * t288 / 0.4e1;
t286 = t115 * t125;
t235 = t119 * t286;
t243 = t115 * t298;
t284 = t116 * t124;
t330 = pkin(1) * t116;
t305 = -0.2e1 * ((0.4e1 * t216 + t275) * t335 + t184 * t352 + (-t148 / 0.2e1 + (t60 * t284 + (-t235 + (t124 * t64 + t125 * t92) * t116) * qJD(2)) * pkin(1)) * pkin(2)) * t349 - 0.2e1 * ((qJD(2) * t198 + t202 * t350 + t290 * t96) * t334 + qJD(2) * t179 + ((t243 / 0.2e1 - t61 * t330) * t124 + ((t278 + (-t114 + t297) * t124) * t334 + (-t124 * t66 - t125 * t93) * t330) * qJD(2)) * pkin(2)) * t195;
t3 = (qJD(2) * t162 * t212 + ((t258 * t60 + 0.4e1 * t251) * t56 + (t42 * t230 + (t256 * t60 + t346 * t61) * t75) * t54) * t118) * pkin(4) + t305;
t341 = t3 / 0.2e1;
t70 = t69 * t71;
t340 = t70 / 0.2e1;
t339 = -t93 / 0.2e1;
t213 = 0.2e1 * t255;
t84 = t82 / t81;
t257 = t84 * t344;
t308 = t84 * t90;
t175 = 0.8e1 * t184;
t324 = pkin(3) * t118;
t253 = t82 * t324;
t196 = t132 * t77 * t253;
t282 = t124 * t125;
t197 = 0.6e1 * pkin(2) * t136 * t282;
t229 = pkin(2) * t254;
t327 = pkin(2) * t116;
t262 = -0.2e1 * t327;
t293 = t125 * t94;
t65 = -t99 * t289 + 0.2e1 * t246 + (t113 * t125 + t280) * pkin(1);
t31 = ((qJD(2) * t197 + t120 * t188 - t96 * t289) * t115 + t94 * t175 + ((t62 * t262 + t243) * t124 + ((t278 + (-t113 + t297) * t124) * t115 + (-t124 * t65 - t293) * t229) * qJD(2)) * pkin(1)) * t196;
t294 = t125 * t91;
t329 = pkin(1) * t132;
t63 = t108 + (-t278 + (t261 - t297) * t124) * pkin(1);
t320 = (((0.4e1 * t215 + t274) * t115 + t91 * t175) * t132 + (t148 + (-0.2e1 * t59 * t284 + (t286 * t342 + (-t124 * t63 - t294) * t254) * qJD(2)) * pkin(2)) * t329) * t94;
t78 = 0.1e1 / t80 ^ 2;
t2 = t31 + (t161 * t213 + (-0.2e1 * t348 * t78 * (-t308 * t63 + t309 * t65) + (t51 * t63 * t257 + (-t51 * t65 - t53 * t63 - t320) * t83) * t77) * t118) * pkin(3);
t205 = -0.2e1 * t223;
t55 = (t115 * t63 + t205 * t91) * t132;
t57 = (t115 * t65 + t205 * t94) * t132;
t177 = t309 * t55 - t57 * t82;
t160 = t177 * t77;
t316 = (-t308 * t59 + t309 * t62) * t78;
t4 = t31 + (-qJD(2) * t160 * t213 + (0.2e1 * t177 * t316 + (t55 * t59 * t257 + (-t55 * t62 - t57 * t59 - t320) * t83) * t77) * t118) * pkin(3);
t336 = t2 - t4;
t333 = Ifges(5,5) * t93;
t332 = Ifges(5,6) * t92;
t328 = pkin(2) * t115;
t33 = t161 * t324 + qJD(2);
t296 = t124 * t91;
t186 = t293 + t296;
t241 = t115 * t132 * t69;
t48 = t186 * t241;
t46 = qJD(1) * t48;
t16 = mrSges(4,1) * t33 + mrSges(4,3) * t46;
t322 = t16 * t91;
t211 = qJD(1) * t241;
t295 = t124 * t94;
t58 = t211 * t295;
t47 = -t211 * t294 + t58;
t17 = -mrSges(4,2) * t33 + mrSges(4,3) * t47;
t321 = t17 * t94;
t319 = t33 * Ifges(4,3);
t34 = t163 * t260;
t318 = t34 * Ifges(5,3);
t158 = t304 * t185;
t164 = 0.2e1 * t59 * t91 + 0.2e1 * t62 * t94;
t39 = (-qJD(2) * t158 + t116 * t164) * t133;
t317 = t39 * t70;
t315 = t46 * Ifges(4,5);
t314 = t47 * Ifges(4,6);
t311 = t68 * t92;
t302 = -t91 - t65;
t301 = Ifges(3,4) * t124;
t300 = Ifges(3,4) * t125;
t299 = pkin(2) * qJD(2);
t292 = Ifges(3,5) * qJD(2);
t291 = Ifges(3,6) * qJD(2);
t279 = t125 * t127;
t277 = t127 * t135;
t276 = t139 * t127;
t271 = qJD(1) * t124;
t270 = qJD(1) * t125;
t268 = qJD(1) * t132;
t267 = qJD(2) * t116;
t245 = t115 * t299;
t242 = t67 * t285;
t239 = t68 * t339;
t238 = t91 * t340;
t237 = t94 * t340;
t236 = qJD(2) * t298;
t233 = t125 * t276;
t219 = 0.2e1 * t236;
t41 = ((t63 * t91 + t65 * t94) * t254 - t158) * t133;
t217 = t41 * t237;
t97 = t98 ^ 2;
t208 = t97 * t231;
t207 = -t242 / 0.2e1;
t206 = t242 / 0.2e1;
t192 = t93 * Ifges(5,1) - t92 * Ifges(5,4);
t191 = t93 * Ifges(5,4) - t92 * Ifges(5,2);
t190 = -t332 + t333;
t189 = pkin(1) * t69 * t229;
t187 = -t294 + t295;
t173 = 0.2e1 * t67 * t200;
t172 = -0.2e1 * t67 * t199;
t171 = t60 * t67 - t38 * t311 / 0.2e1;
t170 = t239 * t38 + t61 * t67;
t169 = t238 * t41 - t63 * t69;
t168 = t127 * t181;
t167 = -t114 * t127 + t219;
t165 = (Ifges(3,2) * t125 + t301) * qJD(1);
t159 = 0.2e1 * t180;
t156 = t70 * (-t295 / 0.2e1 + t294 / 0.2e1);
t155 = t192 * t242;
t154 = t191 * t242;
t153 = t190 * t242;
t95 = t151 * t127;
t152 = -t103 * t95 + t208 + t276;
t150 = (t123 * t91 + t282 * t94) * t189;
t149 = (-t123 * t94 + t282 * t91) * t189;
t9 = (qJD(2) * t150 + ((t296 / 0.2e1 + t293 / 0.2e1) * t317 + (qJD(2) * t187 - t124 * t59 - t125 * t62) * t69) * t115) * t132;
t10 = (qJD(2) * t149 + (t39 * t156 + (qJD(2) * t186 + t124 * t62 - t125 * t59) * t69) * t115) * t132;
t122 = Ifges(3,4) * t270;
t110 = Ifges(3,1) * t271 + t122 + t292;
t109 = t165 + t291;
t44 = -t306 * t64 + t307 * t66;
t40 = ((t64 * t92 + t66 * t93) * t254 - t157) * t130;
t37 = -mrSges(4,1) * t47 - mrSges(4,2) * t46;
t28 = (t172 + (t239 * t40 + t66 * t67) * t115) * t269;
t27 = (t94 * t181 + (-t65 * t69 + t217) * t328) * t264;
t26 = (t173 + (-t64 * t67 + t40 * t311 / 0.2e1) * t115) * t269;
t25 = (t169 * t328 + t181 * t91) * t264;
t24 = (qJD(2) * t172 + t115 * t170) * t269;
t23 = (t94 * t168 + (t237 * t39 - t62 * t69) * t245) * t132;
t22 = (qJD(2) * t173 - t115 * t171) * t269;
t21 = (t91 * t168 + (t238 * t39 - t59 * t69) * t245) * t132;
t20 = -t34 * Ifges(5,5) + qJD(1) * t155;
t19 = -t34 * Ifges(5,6) + qJD(1) * t154;
t15 = -t46 * Ifges(4,1) + t47 * Ifges(4,4) + t33 * Ifges(4,5);
t14 = -t46 * Ifges(4,4) + t47 * Ifges(4,2) + t33 * Ifges(4,6);
t13 = t314 - t315 + t319;
t12 = (t149 + (t41 * t156 + ((-t63 + t94) * t125 - t302 * t124) * t69) * t115) * t268;
t11 = t58 + (t150 + (t169 * t124 + (t302 * t69 + t217) * t125) * t115) * t268;
t8 = qJD(1) * t10;
t7 = qJD(1) * t9;
t6 = ((t120 * t208 + t127 * t197 - t95 * t289) * t115 + 0.8e1 * t94 * t127 * t209 + ((t233 + (-t113 * t127 + t219) * t124) * t115 + (-0.4e1 * t266 * t62 - 0.2e1 * t279 * t94) * t327) * pkin(1)) * t196 - 0.2e1 * t53 * t253 * t316 + 0.2e1 * pkin(3) * t220 * t161 + ((t59 * t77 * t84 + t316 * t83) * t51 * t344 + (-t53 * t59 - ((0.8e1 * t117 * t135 * t91 + 0.4e1 * t328) * t132 * t127 * t283 + ((t236 * t343 + (t115 * t224 + t262 * t91) * t127) * t125 + (-0.4e1 * pkin(2) * t267 * t59 + t115 * t152) * t124) * t329) * t94 - t51 * t62) * t77 * t83) * t324;
t5 = -0.2e1 * ((t127 * t198 + t97 * t202 + t95 * t290) * t334 + t127 * t179 + ((t124 * t167 + t233) * t334 + (-0.2e1 * t266 * t61 - t279 * t93) * t330) * pkin(2)) * t195 - 0.2e1 * ((t117 * t136 * t352 + pkin(1) * t343) * t123 * t277 + ((t124 * t152 - t125 * t167) * t335 + (-t127 * t235 + (0.2e1 * t266 * t60 + t279 * t92) * t116) * pkin(1)) * pkin(2)) * t349 + 0.2e1 * (-t50 * t61 + t52 * t60) * t221 - 0.4e1 * pkin(4) * t220 * t163 + 0.2e1 * (t52 * t251 + (t42 * t76 * t87 + t60 * t75 * t88) * t50 * t93) * t260;
t1 = (-t163 * t212 + ((t258 * t64 + 0.4e1 * t310 * t44) * t52 + (t44 * t230 + (t256 * t64 + t346 * t66) * t75) * t50) * t118) * pkin(4) + t305;
t18 = [(Ifges(3,1) * t125 - t301) * t228 + t22 * t154 / 0.2e1 + t5 * t153 / 0.2e1 + t47 * (Ifges(4,4) * t9 + Ifges(4,2) * t10) / 0.2e1 + (-((-t92 * mrSges(5,1) - t93 * mrSges(5,2)) * t159 + (mrSges(5,1) * t171 + mrSges(5,2) * t170) * t115) * t269 + mrSges(5,1) * t22 - mrSges(5,2) * t24) * pkin(1) - t46 * (Ifges(4,1) * t9 + Ifges(4,4) * t10) / 0.2e1 - pkin(2) * (-mrSges(4,1) * t10 + mrSges(4,2) * t9) * t270 - (-mrSges(4,1) * t8 + mrSges(4,2) * t7) * t325 + t127 * (Ifges(3,5) * t125 - Ifges(3,6) * t124) / 0.2e1 + t24 * t155 / 0.2e1 + (qJD(1) * (t124 * Ifges(3,1) + t300) + t110) * t265 / 0.2e1 - (t165 + t109) * t266 / 0.2e1 + (t60 * t207 + t347 * t92) * t19 + (t61 * t206 - t347 * t93) * t20 + t93 * (Ifges(5,1) * t24 + Ifges(5,4) * t22 + Ifges(5,5) * t5) * t206 + t92 * (Ifges(5,4) * t24 + Ifges(5,2) * t22 + Ifges(5,6) * t5) * t207 + t10 * t14 / 0.2e1 + t9 * t15 / 0.2e1 + t37 * t244 - (mrSges(4,2) * t351 - mrSges(4,3) * t21 + Ifges(4,1) * t7 + Ifges(4,4) * t8 + Ifges(4,5) * t6) * t48 + (-0.2e1 * m(4) * t281 - Ifges(3,2) * t124 + t300) * t125 * t263 + t33 * (Ifges(4,5) * t9 + Ifges(4,6) * t10) / 0.2e1 + ((-t10 * t94 + t9 * t91) * mrSges(4,3) * t299 + (-mrSges(4,1) * t351 + mrSges(4,3) * t23 + Ifges(4,4) * t7 + Ifges(4,2) * t8 + Ifges(4,6) * t6) * t187) * t241 + (-t34 * (-t190 * t159 + (Ifges(5,5) * t170 - Ifges(5,6) * t171) * t115) / 0.2e1 + (t93 * (-t192 * t159 + (Ifges(5,1) * t170 - Ifges(5,4) * t171) * t115) / 0.2e1 - t92 * (-t191 * t159 + (Ifges(5,4) * t170 - Ifges(5,2) * t171) * t115) / 0.2e1) * t115 * t67 * t269) * t129; t46 * (Ifges(4,1) * t11 + Ifges(4,4) * t12 + Ifges(4,5) * t2) / 0.2e1 - t47 * (Ifges(4,4) * t11 + Ifges(4,2) * t12 + Ifges(4,6) * t2) / 0.2e1 - t318 * t341 - t2 * t13 / 0.2e1 - t12 * t14 / 0.2e1 - t11 * t15 / 0.2e1 - t25 * t16 - t26 * t19 / 0.2e1 - t27 * t17 - t28 * t20 / 0.2e1 - t33 * (Ifges(4,5) * t11 + Ifges(4,6) * t12 + Ifges(4,3) * t2) / 0.2e1 + t34 * (Ifges(5,5) * t28 + Ifges(5,6) * t26 + Ifges(5,3) * t1) / 0.2e1 + (t341 - t1 / 0.2e1) * (qJD(1) * t153 - t318) + (t319 / 0.2e1 - t315 / 0.2e1 + t314 / 0.2e1 + t13 / 0.2e1) * t4 + m(4) * (-t71 * t158 * t277 + (t71 * t164 - t304 / t74 ^ 2 * t39) * t135 * t267) * t133 / 0.2e1 + ((t321 / 0.2e1 + t322 / 0.2e1) * t317 + (-t59 * t16 - t62 * t17 + (-m(4) * t23 + mrSges(4,2) * t6 - mrSges(4,3) * t8 + (m(4) * t27 - mrSges(4,2) * t336 + mrSges(4,3) * t12) * qJD(2)) * t94 + (-m(4) * t21 - mrSges(4,1) * t6 + mrSges(4,3) * t7 + (m(4) * t25 + mrSges(4,1) * t336 - mrSges(4,3) * t11) * qJD(2)) * t91) * t69) * t328 * t132 + (mrSges(4,1) * t21 - mrSges(4,2) * t23 + Ifges(4,5) * t7 + Ifges(4,6) * t8 + Ifges(4,3) * t6) * (-t160 * t324 + 0.1e1) + (pkin(1) * (-mrSges(5,1) * t26 + mrSges(5,2) * t28) + (t292 / 0.2e1 - t110 / 0.2e1 - t122 / 0.2e1 + pkin(2) * (-mrSges(4,1) * t12 + mrSges(4,2) * t11)) * t125 + (t92 * (Ifges(5,4) * t28 + Ifges(5,2) * t26 + Ifges(5,6) * t1) / 0.2e1 + (Ifges(5,1) * t28 + Ifges(5,4) * t26 + Ifges(5,5) * t1) * t339 + (t333 / 0.2e1 - t332 / 0.2e1) * t3) * t242 + (-t291 / 0.2e1 + t109 / 0.2e1 - pkin(2) * t37 + Ifges(3,4) * t271 / 0.2e1 + (m(4) * t135 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t270) * t124) * qJD(1) + (Ifges(5,5) * t24 + Ifges(5,6) * t22 + Ifges(5,3) * t5) * t162 * t260 + (t321 + t322) * t264 * t181;];
tauc = t18(:);
