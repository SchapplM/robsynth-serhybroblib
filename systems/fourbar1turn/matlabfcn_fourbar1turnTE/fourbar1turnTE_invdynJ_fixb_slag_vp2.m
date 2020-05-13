% Calculate vector of inverse dynamics joint torques for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnTE_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:34
% EndTime: 2020-04-12 19:19:06
% DurationCPUTime: 13.78s
% Computational Cost: add. (94617->428), mult. (134147->913), div. (4076->15), fcn. (35917->6), ass. (0->314)
t131 = pkin(2) ^ 2;
t132 = pkin(1) ^ 2;
t121 = cos(qJ(2));
t333 = pkin(2) * t121;
t278 = -0.2e1 * pkin(1) * t333 + t132;
t113 = t131 + t278;
t391 = pkin(3) ^ 2;
t277 = -pkin(4) ^ 2 + t391;
t106 = t113 + t277;
t115 = pkin(1) * t121 - pkin(2);
t119 = sin(qJ(2));
t363 = -pkin(3) - pkin(4);
t104 = (pkin(2) - t363) * (pkin(2) + t363) + t278;
t362 = pkin(4) - pkin(3);
t105 = (pkin(2) - t362) * (pkin(2) + t362) + t278;
t297 = t104 * t105;
t133 = sqrt(-t297);
t287 = t119 * t133;
t85 = -pkin(1) * t287 - t106 * t115;
t138 = t85 ^ 2;
t78 = 0.1e1 / t138;
t340 = pkin(1) * t106;
t101 = t119 * t340;
t88 = -t115 * t133 + t101;
t322 = t78 * t88;
t110 = 0.1e1 / t113;
t128 = 0.1e1 / pkin(3);
t111 = 0.1e1 / t113 ^ 2;
t293 = t111 * t119;
t246 = pkin(1) * t293;
t212 = qJD(2) * t246;
t195 = pkin(2) * t212;
t183 = -0.2e1 * t195;
t377 = -0.2e1 * t115;
t262 = pkin(2) * t377;
t223 = t106 + t262;
t207 = (-t104 - t105) * pkin(1) * pkin(2);
t92 = t119 * t207;
t91 = qJD(2) * t92;
t96 = 0.1e1 / t133;
t317 = t91 * t96;
t243 = t119 * t317;
t283 = t121 * t133;
t50 = (-t243 + (t119 * t223 - t283) * qJD(2)) * pkin(1);
t39 = (t110 * t50 + t183 * t85) * t128;
t118 = t119 ^ 2;
t290 = t118 * t132;
t234 = qJD(2) * t290;
t209 = pkin(2) * t234;
t305 = t115 * t96;
t271 = qJD(2) * t119;
t233 = t133 * t271;
t270 = qJD(2) * t121;
t310 = pkin(1) * t233 + t270 * t340;
t53 = -t305 * t91 + 0.2e1 * t209 + t310;
t41 = (t110 * t53 + t183 * t88) * t128;
t77 = 0.1e1 / t85;
t392 = -t39 * t322 + t41 * t77;
t84 = t88 ^ 2;
t67 = t78 * t84 + 0.1e1;
t64 = 0.1e1 / t67;
t162 = t392 * t64;
t107 = t113 - t277;
t114 = pkin(1) - t333;
t86 = -pkin(2) * t287 + t107 * t114;
t81 = 0.1e1 / t86 ^ 2;
t334 = pkin(2) * t119;
t100 = t107 * t334;
t87 = t114 * t133 + t100;
t320 = t81 * t87;
t126 = 0.1e1 / pkin(4);
t356 = t110 / 0.2e1;
t260 = 0.2e1 * t114 * pkin(1);
t51 = (-t243 + (-t283 + (t107 + t260) * t119) * qJD(2)) * pkin(2);
t147 = -t195 * t86 + t356 * t51;
t38 = t147 * t126;
t291 = t118 * t131;
t242 = pkin(1) * t291;
t210 = qJD(2) * t242;
t306 = t114 * t96;
t238 = pkin(2) * t270;
t311 = pkin(2) * t233 + t107 * t238;
t52 = t306 * t91 + 0.2e1 * t210 + t311;
t146 = -t195 * t87 + t356 * t52;
t40 = t146 * t126;
t83 = t87 ^ 2;
t66 = t81 * t83 + 0.1e1;
t62 = 0.1e1 / t66;
t80 = 0.1e1 / t86;
t164 = t62 * (t38 * t320 - t40 * t80);
t331 = pkin(4) * t113;
t261 = 0.2e1 * t331;
t294 = t110 * t128;
t230 = -t294 / 0.2e1;
t206 = t85 * t230;
t229 = t294 / 0.2e1;
t301 = t119 * t88;
t61 = t121 * t206 + t229 * t301;
t295 = t110 * t126;
t364 = t87 / 0.2e1;
t365 = t86 / 0.2e1;
t152 = (mrSges(5,1) * t365 + mrSges(5,2) * t364) * t295;
t385 = -m(4) * t333 - mrSges(3,1) * t121 + mrSges(3,2) * t119;
t394 = m(5) * pkin(1) + mrSges(2,1) - t152 - t385;
t390 = -0.4e1 * t86;
t389 = 0.6e1 * t119;
t239 = pkin(2) * t271;
t217 = pkin(1) * t239;
t388 = -0.2e1 * t217;
t387 = 2 * qJDD(2);
t337 = pkin(1) * t128;
t250 = t62 * t331;
t219 = t81 * t250;
t386 = t126 * t87 * t219;
t120 = sin(qJ(1));
t122 = cos(qJ(1));
t384 = g(1) * t122 + g(2) * t120;
t383 = -mrSges(4,3) + mrSges(2,2) - mrSges(3,3) - mrSges(5,3);
t381 = -0.2e1 * t81;
t380 = 0.4e1 * t87;
t379 = 0.2e1 * t88;
t378 = t126 ^ 2;
t376 = qJD(1) ^ 2;
t124 = qJD(2) ^ 2;
t222 = pkin(2) * t246;
t357 = -t110 / 0.2e1;
t316 = t96 * t92;
t55 = t100 + (-t283 + (t260 - t316) * t119) * pkin(2);
t157 = t222 * t86 + t357 * t55;
t45 = t157 * t126;
t59 = t92 * t306 + 0.2e1 * t242 + (t107 * t121 + t287) * pkin(2);
t156 = -t222 * t87 + t356 * t59;
t48 = t156 * t126;
t163 = t62 * (-t320 * t45 - t48 * t80);
t63 = 0.1e1 / t66 ^ 2;
t228 = 0.4e1 * t63 * t320;
t82 = t80 * t81;
t319 = t82 * t83;
t29 = -t319 * t51 + t320 * t52;
t323 = t63 * t80;
t249 = t29 * t323;
t257 = t82 * t380;
t259 = 0.2e1 * t62 * t81;
t298 = 0.4e1 * t96 / t297;
t216 = t91 * t92 * t298;
t193 = -t216 / 0.4e1;
t315 = -qJD(2) / 0.2e1;
t355 = -t119 / 0.2e1;
t215 = -0.4e1 * t131 * t290;
t89 = (t121 * t207 + t215) * qJD(2);
t143 = (t119 * t193 + 0.2e1 * (t89 * t355 + (-t91 / 0.2e1 + t92 * t315) * t121) * t96) * t110;
t112 = t110 * t111;
t292 = t112 * t131;
t185 = t234 * t292;
t196 = t126 * t80 * t250;
t284 = t121 * t131;
t235 = t119 * t284;
t197 = 0.6e1 * pkin(1) * t235;
t241 = t110 * t317;
t286 = t121 * t110;
t338 = pkin(1) * t111;
t314 = -0.2e1 * ((0.4e1 * t210 + t311) * t357 + t185 * t390 + (-t143 / 0.2e1 + (t51 * t293 + (-t114 * t286 + (t119 * t55 + t121 * t86) * t111) * qJD(2)) * pkin(1)) * pkin(2)) * t386 - 0.2e1 * ((t114 * t216 / 0.4e1 + t89 * t306 + qJD(2) * t197) * t356 + t185 * t380 + ((t241 / 0.2e1 - t52 * t338) * t119 + ((t283 + (-t107 + t316) * t119) * t356 + (-t119 * t59 - t121 * t87) * t338) * qJD(2)) * pkin(2)) * t196;
t3 = (0.4e1 * t163 * t217 + ((t259 * t51 + 0.4e1 * t249) * t48 + (t29 * t228 + (t257 * t51 + t381 * t52) * t62) * t45) * t113) * pkin(4) + t314;
t375 = t3 / 0.2e1;
t366 = -t86 / 0.2e1;
t155 = (Ifges(5,1) * t364 + Ifges(5,4) * t366) * t295;
t20 = t164 * t261;
t17 = Ifges(5,5) * t20 + qJD(1) * t155;
t371 = -t17 / 0.2e1;
t368 = -t85 / 0.2e1;
t367 = t85 / 0.2e1;
t267 = qJD(1) * qJD(2);
t108 = qJDD(1) * t121 - t119 * t267;
t359 = -t108 / 0.2e1;
t227 = t121 * t267;
t109 = qJDD(1) * t119 + t227;
t358 = t109 / 0.2e1;
t354 = t119 / 0.2e1;
t353 = -t121 / 0.2e1;
t352 = t121 / 0.2e1;
t351 = mrSges(4,1) * t61;
t299 = t121 * t88;
t302 = t119 * t85;
t168 = t299 / 0.2e1 + t302 / 0.2e1;
t60 = t168 * t294;
t350 = mrSges(4,2) * t60;
t347 = Ifges(5,5) * t87;
t344 = Ifges(5,6) * t86;
t343 = Ifges(5,3) * t20;
t339 = pkin(1) * t110;
t336 = pkin(2) * t110;
t335 = pkin(2) * t111;
t332 = pkin(3) * t113;
t179 = 0.8e1 * t185;
t252 = 0.2e1 * t111;
t300 = t121 * t85;
t54 = t101 + (-t283 + (t262 - t316) * t119) * pkin(1);
t304 = t119 * t54;
t330 = (((0.4e1 * t209 + t310) * t110 + t85 * t179) * t128 + (t143 + (-0.2e1 * t50 * t293 + (t286 * t377 + (-t300 - t304) * t252) * qJD(2)) * pkin(2)) * t337) * t88;
t318 = t88 * t53;
t79 = t77 * t78;
t321 = t79 * t84;
t65 = 0.1e1 / t67 ^ 2;
t328 = (t318 * t78 - t321 * t50) * t65;
t205 = t88 * t230;
t313 = (t119 * t205 + t229 * t300) * t120;
t312 = (t299 + t302) * t122 * t229;
t57 = t61 * qJD(1);
t309 = Ifges(3,4) * t119;
t308 = Ifges(3,4) * t121;
t307 = pkin(1) * qJD(1);
t58 = -t92 * t305 + 0.2e1 * pkin(2) * t290 + (t106 * t121 + t287) * pkin(1);
t303 = t119 * t58;
t296 = t110 * t115;
t289 = t119 * t121;
t288 = t119 * t131;
t285 = t121 * t124;
t282 = t124 * t132;
t281 = t124 * t133;
t280 = t53 * qJD(2);
t276 = qJD(1) * t110;
t275 = qJD(1) * t119;
t274 = qJD(1) * t121;
t273 = qJD(1) * t128;
t272 = qJD(2) * t111;
t269 = qJD(2) * t128;
t268 = t88 * qJDD(2);
t266 = qJDD(1) * t110;
t265 = qJDD(2) * t106;
t264 = qJDD(2) * t107;
t263 = -0.2e1 * t335;
t258 = t79 * t379;
t255 = m(4) * t334;
t254 = pkin(1) * t335;
t253 = pkin(1) * t334;
t251 = t77 * t332;
t245 = pkin(1) * t288;
t244 = pkin(2) * t294;
t240 = pkin(2) * t274;
t237 = -t336 / 0.2e1;
t236 = t124 * t291;
t232 = -t295 / 0.4e1;
t231 = t295 / 0.4e1;
t226 = -t276 / 0.4e1;
t225 = t276 / 0.4e1;
t221 = t64 * t251;
t214 = t111 * t245;
t213 = t124 * t245;
t211 = qJD(2) * t244;
t208 = pkin(2) * t230;
t204 = -0.2e1 * t222;
t194 = t111 * t213;
t190 = Ifges(3,2) * t121 + t309;
t189 = Ifges(3,5) * t121 - Ifges(3,6) * t119;
t188 = qJD(2) * t208;
t187 = t211 / 0.2e1;
t186 = -t110 * t114 + t111 * t86;
t47 = (t54 * t110 + t204 * t85) * t128;
t49 = (t58 * t110 + t204 * t88) * t128;
t181 = t322 * t47 - t49 * t77;
t180 = t126 * t195;
t176 = t128 * t131 * t212;
t175 = t85 * t187;
t173 = t118 * t85 + t289 * t88;
t172 = -t118 * t88 + t289 * t85;
t171 = t54 * t353 + t303 / 0.2e1;
t170 = t58 * t352 + t304 / 0.2e1;
t169 = -t300 / 0.2e1 + t301 / 0.2e1;
t167 = t91 ^ 2 * t298 / 0.4e1 + t96 * (t124 * t215 + (qJDD(2) * t119 + t285) * t207);
t166 = 0.2e1 * qJD(2) * t317 + qJDD(2) * t133;
t165 = t119 * (Ifges(3,1) * t121 - t309);
t161 = t181 * t64;
t160 = (t353 * t53 + t355 * t50) * t110;
t159 = (t353 * t50 + t354 * t53) * t110;
t158 = t169 * t110;
t154 = (Ifges(5,4) * t364 + Ifges(5,2) * t366) * t295;
t153 = (t347 / 0.2e1 - t344 / 0.2e1) * t295;
t151 = -t167 + t281;
t150 = t173 * t254;
t149 = t172 * t254;
t148 = -t107 * t124 + t166;
t145 = t168 + t171;
t144 = (mrSges(4,1) * t172 - mrSges(4,2) * t173) * t254;
t117 = Ifges(3,4) * t274;
t103 = Ifges(3,1) * t275 + Ifges(3,5) * qJD(2) + t117;
t102 = Ifges(3,6) * qJD(2) + qJD(1) * t190;
t56 = qJD(1) * t60;
t46 = qJD(1) * t48;
t44 = (t214 * t88 + t237 * t58) * t269;
t43 = qJD(1) * t45;
t42 = (t214 * t85 + t237 * t54) * t269;
t37 = -mrSges(4,1) * t57 - mrSges(4,2) * t56;
t36 = (qJD(1) * t146 + t266 * t364) * t126;
t35 = (t88 * t194 + (-t280 / 0.2e1 - t268 / 0.2e1) * t336) * t128;
t34 = (-qJD(1) * t147 + t266 * t366) * t126;
t33 = (t85 * t194 + (qJDD(2) * t368 + t315 * t50) * t336) * t128;
t31 = -t319 * t55 + t320 * t59;
t28 = (t149 + ((t88 / 0.2e1 - t54 / 0.2e1) * t121 + (t58 / 0.2e1 + t367) * t119) * t110) * t273;
t27 = (-t170 * t110 + t150) * t273 + t57;
t26 = (t159 + (t110 * t168 + t149) * qJD(2)) * t128;
t25 = (t160 + (t158 + t150) * qJD(2)) * t128;
t24 = ((t358 * t88 + t359 * t85) * t110 + (qJD(2) * t149 + t159) * qJD(1)) * t128;
t23 = ((t109 * t368 + t359 * t88) * t110 + (qJD(2) * t150 + t160) * qJD(1)) * t128;
t19 = t162 * t332 + qJD(2);
t16 = Ifges(5,6) * t20 + qJD(1) * t154;
t12 = ((t132 * t238 * t389 + t115 * t193 - t305 * t89) * t110 + t88 * t179 + ((t263 * t53 + t241) * t119 + ((t283 + (-t106 + t316) * t119) * t110 + (-t299 - t303) * pkin(2) * t252) * qJD(2)) * pkin(1)) * t128 * t221;
t9 = -t56 * Ifges(4,1) + t57 * Ifges(4,4) + t19 * Ifges(4,5);
t8 = -t56 * Ifges(4,4) + t57 * Ifges(4,2) + t19 * Ifges(4,6);
t6 = qJDD(2) + ((-t167 * t296 + (0.8e1 * t88 * t112 * t236 + (t118 * t387 + t285 * t389) * t336) * t132) * t128 + (((t265 + t281) * t110 + t88 * t124 * t263) * t121 + ((-t106 * t124 + t166) * t110 + (-0.2e1 * t268 - 0.4e1 * t280) * t335) * t119) * t337) * t221 - 0.2e1 * t41 * t251 * t328 + 0.2e1 * pkin(3) * t217 * t162 + ((t50 * t64 * t79 + t328 * t78) * t39 * t379 + (-t41 * t50 - ((0.8e1 * t292 * t85 + 0.4e1 * t336) * t128 * t118 * t282 + ((-t166 * t110 + (t110 * t223 + t263 * t85) * t124) * t121 + ((t151 + t265) * t110 + (-0.4e1 * t50 * t272 + (-t111 * t85 - t296) * t387) * pkin(2)) * t119) * t337) * t88 - t39 * t53) * t64 * t78) * t332;
t5 = -0.2e1 * ((t167 * t114 + t124 * t197) * t356 + (t112 * t282 * t380 + qJDD(2) * t339) * t291 + (((t264 + t281) * t121 + t148 * t119) * t356 + (-t87 * t285 + (-0.2e1 * qJD(2) * t52 - qJDD(2) * t87) * t119) * t338) * pkin(2)) * t196 - 0.2e1 * ((t112 * t132 * t390 - 0.2e1 * t339) * t236 + ((pkin(1) * t124 * t186 + t148 * t356) * t121 + ((t151 + t264) * t357 + (qJDD(2) * t186 + 0.2e1 * t272 * t51) * pkin(1)) * t119) * pkin(2)) * t386 + 0.2e1 * (t38 * t52 + t40 * t51) * t219 - 0.2e1 * pkin(4) * t388 * t164 + 0.2e1 * (t40 * t249 + (-t29 * t63 * t81 - t51 * t62 * t82) * t38 * t87) * t261;
t4 = t12 + (t161 * t388 + (0.2e1 * t181 * t328 + (t47 * t50 * t258 + (-t47 * t53 - t49 * t50 - t330) * t78) * t64) * t113) * pkin(3);
t2 = t12 + (0.2e1 * t162 * t253 + (-0.2e1 * t392 * t65 * (-t321 * t54 + t322 * t58) + (t39 * t54 * t258 + (-t39 * t58 - t41 * t54 - t330) * t78) * t64) * t113) * pkin(3);
t1 = (0.4e1 * t164 * t253 + ((t259 * t55 + 0.4e1 * t31 * t323) * t40 - (t31 * t228 + (t257 * t55 + t381 * t59) * t62) * t38) * t113) * pkin(4) + t314;
t7 = [(t378 * (Ifges(5,1) * t146 - Ifges(5,4) * t147) * t225 + t180 * t371 + (Ifges(5,1) * t36 + Ifges(5,4) * t34 + Ifges(5,5) * t5) * t231) * t87 + (t378 * (Ifges(5,4) * t146 - Ifges(5,2) * t147) * t226 + (Ifges(5,4) * t36 + Ifges(5,2) * t34 + Ifges(5,6) * t5) * t232) * t86 + (mrSges(4,3) * t35 + Ifges(4,4) * t23 + Ifges(4,2) * t24 + Ifges(4,6) * t6) * t61 - (-mrSges(4,3) * t33 + Ifges(4,1) * t23 + Ifges(4,4) * t24 + Ifges(4,5) * t6) * t60 + (t165 + t121 * (-Ifges(3,2) * t119 + t308)) * t267 / 0.2e1 + (m(4) * t284 + Ifges(3,2) * t352 + Ifges(3,4) * t354 + t190 / 0.2e1 - pkin(2) * (-t350 - t351)) * t108 + (t180 * t365 + t232 * t51) * t16 + (t20 * (Ifges(5,5) * t146 - Ifges(5,6) * t147) / 0.2e1 - (mrSges(5,1) * t147 + mrSges(5,2) * t146) * t307) * t126 + t103 * t270 / 0.2e1 - t102 * t271 / 0.2e1 + t19 * (Ifges(4,5) * t25 + Ifges(4,6) * t26) / 0.2e1 + t26 * t8 / 0.2e1 + qJDD(2) * (Ifges(3,5) * t119 + Ifges(3,6) * t121) / 0.2e1 - (-mrSges(4,1) * t26 + mrSges(4,2) * t25) * t240 + (Ifges(3,4) * t109 + Ifges(3,6) * qJDD(2)) * t352 + t124 * t189 / 0.2e1 - (-mrSges(4,1) * t24 + mrSges(4,2) * t23) * t333 + (Ifges(3,1) * t109 + Ifges(3,5) * qJDD(2)) * t354 + t5 * t153 / 0.2e1 + t57 * (Ifges(4,4) * t25 + Ifges(4,2) * t26) / 0.2e1 - t56 * (Ifges(4,1) * t25 + Ifges(4,4) * t26) / 0.2e1 + t25 * t9 / 0.2e1 - pkin(1) * (-mrSges(5,1) * t34 + mrSges(5,2) * t36) + (-t312 * mrSges(4,2) + t383 * t120 + (-mrSges(4,1) * t128 * t158 - t394) * t122) * g(2) + (-t313 * mrSges(4,1) + t383 * t122 + (t350 + t394) * t120) * g(1) - m(4) * t227 * t288 + (m(5) * t132 - pkin(1) * t152 + Ifges(2,3)) * qJDD(1) + (t119 * Ifges(3,1) + t308) * t358 + t52 * t17 * t231 + t37 * t239 + (t188 * t26 * t88 + t175 * t25) * mrSges(4,3) + t34 * t154 / 0.2e1 + t36 * t155 / 0.2e1; (-t351 - ((mrSges(4,1) * t173 + mrSges(4,2) * t172) * t254 + (-mrSges(4,1) * t170 + mrSges(4,2) * t145) * t110) * t128 + t385) * g(3) + ((-mrSges(5,1) * t156 - mrSges(5,2) * t157) * g(3) - t384 * (mrSges(5,1) * t157 - mrSges(5,2) * t156) + (t3 * t347 + t86 * (Ifges(5,4) * t46 + Ifges(5,2) * t43 + Ifges(5,6) * t1)) * t225 + (t3 * t344 + t87 * (Ifges(5,1) * t46 + Ifges(5,4) * t43 + Ifges(5,5) * t1)) * t226) * t126 - (-Ifges(3,2) * t275 + t103 + t117) * t274 / 0.2e1 + (t176 * t88 + t208 * t53 - t44) * (-mrSges(4,2) * t19 + mrSges(4,3) * t57) + (t375 - t1 / 0.2e1) * (qJD(1) * t153 + t343) + t384 * (mrSges(3,1) * t119 + mrSges(3,2) * t121) + (mrSges(4,3) * t28 + (-t2 + t4) * mrSges(4,2)) * t88 * t187 + (t176 * t85 + t208 * t50 - t42) * (mrSges(4,1) * t19 + mrSges(4,3) * t56) + (t4 - t2 / 0.2e1) * (-t56 * Ifges(4,5) + t57 * Ifges(4,6) + t19 * Ifges(4,3)) - g(1) * (t312 * mrSges(4,1) + (-t255 + (t144 + (t171 * mrSges(4,1) + (-t169 + t170) * mrSges(4,2)) * t110) * t128) * t122) - g(2) * (t313 * mrSges(4,2) + (-t255 + (t144 + (mrSges(4,1) * t145 + mrSges(4,2) * t170) * t110) * t128) * t120) + (-t37 * t275 + (-mrSges(4,2) * t6 + mrSges(4,3) * t24) * t205 + (mrSges(4,1) * t6 - mrSges(4,3) * t23) * t206) * pkin(2) + (-mrSges(5,1) * t43 + mrSges(5,2) * t46) * t307 + (t33 * mrSges(4,1) - t35 * mrSges(4,2) + Ifges(4,5) * t23 + Ifges(4,6) * t24 + Ifges(4,3) * t6) * (-t161 * t332 + 0.1e1) + t102 * t275 / 0.2e1 - t27 * t9 / 0.2e1 - t19 * (Ifges(4,5) * t27 + Ifges(4,6) * t28 + Ifges(4,3) * t2) / 0.2e1 - t28 * t8 / 0.2e1 - t376 * t165 / 0.2e1 + Ifges(3,6) * t108 + Ifges(3,5) * t109 - t57 * (Ifges(4,4) * t27 + Ifges(4,2) * t28 + Ifges(4,6) * t2) / 0.2e1 + t56 * (Ifges(4,1) * t27 + Ifges(4,4) * t28 + Ifges(4,5) * t2) / 0.2e1 - t189 * t267 / 0.2e1 - t43 * t16 / 0.2e1 - t20 * (Ifges(5,5) * t46 + Ifges(5,6) * t43 + Ifges(5,3) * t1) / 0.2e1 + (Ifges(5,5) * t36 + Ifges(5,6) * t34 + Ifges(5,3) * t5) * t163 * t261 + (Ifges(3,3) * qJDD(2)) + ((-t33 * t85 - t35 * t88) * t244 / 0.2e1 + ((t318 / 0.2e1 + t50 * t367) * t131 * t272 + (-t138 - t84) * pkin(2) * t112 * t213) / t391 / 0.2e1 + t376 * t235 - (-t42 * t85 - t44 * t88) * t211 / 0.2e1) * m(4) + t46 * t371 + t343 * t375 + (-mrSges(4,1) * t28 + mrSges(4,2) * t27) * t240 + (mrSges(4,1) * t2 - mrSges(4,3) * t27) * t175 + t85 * t4 * mrSges(4,1) * t188;];
tau = t7;
