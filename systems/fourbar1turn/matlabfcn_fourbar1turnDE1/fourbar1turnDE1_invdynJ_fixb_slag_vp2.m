% Calculate vector of inverse dynamics joint torques for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnDE1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:33
% EndTime: 2020-04-12 19:26:22
% DurationCPUTime: 22.23s
% Computational Cost: add. (153153->445), mult. (220445->960), div. (8234->21), fcn. (59337->10), ass. (0->334)
t433 = 2 * pkin(2);
t143 = pkin(2) ^ 2;
t144 = pkin(1) ^ 2;
t132 = cos(qJ(2));
t375 = pkin(2) * t132;
t314 = -0.2e1 * pkin(1) * t375 + t144;
t124 = t143 + t314;
t121 = 0.1e1 / t124;
t390 = -t121 / 0.2e1;
t125 = pkin(1) - t375;
t432 = t121 * t125;
t425 = pkin(3) ^ 2;
t426 = pkin(4) ^ 2;
t313 = t425 - t426;
t117 = t124 + t313;
t126 = pkin(1) * t132 - pkin(2);
t130 = sin(qJ(2));
t396 = (-pkin(3) - pkin(4));
t115 = ((pkin(2) - t396) * (pkin(2) + t396)) + t314;
t395 = (pkin(4) - pkin(3));
t116 = ((pkin(2) - t395) * (pkin(2) + t395)) + t314;
t334 = t115 * t116;
t147 = sqrt(-t334);
t324 = t130 * t147;
t96 = -pkin(1) * t324 - t117 * t126;
t86 = t96 ^ 2;
t87 = 0.1e1 / t96;
t89 = t87 / t86;
t382 = pkin(1) * t117;
t112 = t130 * t382;
t99 = -t126 * t147 + t112;
t95 = t99 ^ 2;
t358 = t89 * t95;
t88 = 0.1e1 / t96 ^ 2;
t359 = t88 * t99;
t299 = -0.2e1 * pkin(2) * t126;
t262 = t117 + t299;
t246 = (-t115 - t116) * pkin(1) * pkin(2);
t103 = t130 * t246;
t102 = qJD(2) * t103;
t107 = 0.1e1 / t147;
t339 = t102 * t107;
t273 = t130 * t339;
t321 = t132 * t147;
t63 = (-t273 + (t130 * t262 - t321) * qJD(2)) * pkin(1);
t307 = qJD(2) * t130;
t268 = t147 * t307;
t129 = t130 ^ 2;
t326 = t129 * t144;
t269 = qJD(2) * t326;
t306 = qJD(2) * t132;
t336 = t107 * t126;
t66 = pkin(1) * t268 - t102 * t336 + t269 * t433 + t306 * t382;
t84 = t88 * t95 + 0.1e1;
t82 = 0.1e1 / t84 ^ 2;
t367 = (-t358 * t63 + t359 * t66) * t82;
t431 = -0.2e1 * t367;
t389 = t121 / 0.2e1;
t140 = 0.1e1 / pkin(3);
t122 = 0.1e1 / t124 ^ 2;
t329 = t122 * t130;
t284 = pkin(1) * t329;
t261 = pkin(2) * t284;
t240 = -0.2e1 * t261;
t338 = t107 * t103;
t67 = t112 + (-t321 + (t299 - t338) * t130) * pkin(1);
t55 = (t121 * t67 + t240 * t96) * t140;
t69 = -t103 * t336 + t326 * t433 + (t117 * t132 + t324) * pkin(1);
t57 = (t121 * t69 + t240 * t99) * t140;
t206 = t359 * t55 - t57 * t87;
t81 = 0.1e1 / t84;
t177 = t206 * t81;
t281 = pkin(2) * t307;
t256 = pkin(1) * t281;
t230 = t122 * t256;
t209 = -0.2e1 * t230;
t51 = (t121 * t63 + t209 * t96) * t140;
t53 = (t121 * t66 + t209 * t99) * t140;
t428 = -t51 * t359 + t53 * t87;
t178 = t428 * t81;
t376 = pkin(2) * t130;
t291 = pkin(1) * t376;
t245 = 0.2e1 * t291;
t406 = 0.2e1 * t99;
t421 = ((qJD(2) * t177 + t178) * t245 + (-0.2e1 * t428 * t82 * (-t358 * t67 + t359 * t69) + t206 * t431 + ((-t51 * t69 - t53 * t67 + t55 * t66 + t57 * t63) * t88 + (t51 * t67 - t55 * t63) * t89 * t406) * t81) * t124) * pkin(3);
t118 = t124 - t313;
t97 = -pkin(2) * t324 + t118 * t125;
t92 = 0.1e1 / t97 ^ 2;
t111 = t118 * t376;
t98 = t125 * t147 + t111;
t357 = t92 * t98;
t137 = 0.1e1 / pkin(4);
t297 = 0.2e1 * t125 * pkin(1);
t64 = (-t273 + (-t321 + (t118 + t297) * t130) * qJD(2)) * pkin(2);
t50 = (t230 * t97 + t390 * t64) * t137;
t327 = t129 * t143;
t283 = pkin(1) * t327;
t248 = qJD(2) * t283;
t317 = (t118 * t306 + t268) * pkin(2);
t337 = t107 * t125;
t65 = t102 * t337 + 0.2e1 * t248 + t317;
t52 = (-t230 * t98 + t389 * t65) * t137;
t94 = t98 ^ 2;
t83 = t92 * t94 + 0.1e1;
t79 = 0.1e1 / t83;
t91 = 0.1e1 / t97;
t180 = t79 * (t50 * t357 + t52 * t91);
t373 = pkin(4) * t124;
t298 = 0.2e1 * t373;
t340 = t132 * t99;
t343 = t130 * t96;
t218 = t340 + t343;
t330 = t121 * t140;
t141 = 0.1e1 / t425;
t351 = t86 + t95;
t78 = t351 * t141 * t122;
t73 = t78 ^ (-0.1e1 / 0.2e1);
t278 = t73 * t330;
t48 = t218 * t278;
t424 = -0.4e1 * t97;
t422 = 2 * qJDD(2);
t303 = qJD(1) * qJD(2);
t119 = qJDD(1) * t132 - t130 * t303;
t420 = t119 / 0.2e1;
t379 = pkin(1) * t140;
t419 = pkin(2) * t119;
t418 = Ifges(3,2) * t132;
t417 = t122 * t245;
t288 = t79 * t373;
t258 = t92 * t288;
t416 = t137 * t98 * t258;
t208 = 0.2e1 * t143 * t73 * t284;
t305 = qJD(2) * t140;
t123 = t121 * t122;
t244 = 0.4e1 * t291;
t217 = t123 * t244;
t175 = t351 * t217;
t181 = 0.2e1 * t63 * t96 + 0.2e1 * t66 * t99;
t39 = (-qJD(2) * t175 + t122 * t181) * t141;
t75 = 0.1e1 / t78;
t74 = t73 * t75;
t369 = t39 * t74;
t415 = pkin(2) * t330 * t369 / 0.2e1 + t208 * t305;
t414 = -m(4) * t375 - mrSges(3,1) * t132 + mrSges(3,2) * t130;
t131 = sin(qJ(1));
t133 = cos(qJ(1));
t413 = g(1) * t133 + g(2) * t131;
t331 = t121 * t137;
t138 = 0.1e1 / t426;
t90 = t97 ^ 2;
t350 = t90 + t94;
t174 = t350 * t217;
t290 = 0.2e1 * t122;
t38 = ((t64 * t97 + t65 * t98) * t290 - qJD(2) * t174) * t138;
t77 = t350 * t138 * t122;
t71 = t77 ^ (-0.1e1 / 0.2e1);
t72 = t71 / t77;
t412 = t38 * t72 * t331 / 0.4e1 + t137 * t71 * t230;
t411 = -mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + mrSges(2,2);
t410 = m(5) * pkin(1) + mrSges(2,1) - t414;
t352 = t69 + t96;
t353 = t67 - t99;
t409 = mrSges(4,1) * t352 + mrSges(4,2) * t353;
t408 = -0.2e1 * t92;
t407 = 0.4e1 * t98;
t405 = t137 ^ 2;
t403 = qJD(1) ^ 2;
t135 = qJD(2) ^ 2;
t234 = t97 * t261;
t68 = t111 + (-t321 + (t297 - t338) * t130) * pkin(2);
t54 = (t390 * t68 + t234) * t137;
t233 = t98 * t261;
t70 = t103 * t337 + 0.2e1 * t283 + (t118 * t132 + t324) * pkin(2);
t56 = (t389 * t70 - t233) * t137;
t179 = t79 * (-t357 * t54 - t56 * t91);
t80 = 0.1e1 / t83 ^ 2;
t266 = 0.4e1 * t80 * t357;
t360 = t80 * t91;
t93 = t91 / t90;
t356 = t93 * t94;
t42 = -t356 * t64 + t357 * t65;
t287 = t42 * t360;
t294 = t93 * t407;
t296 = 0.2e1 * t79 * t92;
t254 = -0.4e1 * t143 * t326;
t100 = (t132 * t246 + t254) * qJD(2);
t328 = t123 * t143;
t214 = t269 * t328;
t231 = t137 * t91 * t288;
t322 = t132 * t143;
t272 = t130 * t322;
t232 = 0.6e1 * pkin(1) * t272;
t335 = 0.4e1 * t107 / t334;
t243 = t102 * t103 * t335;
t380 = pkin(1) * t122;
t355 = -0.2e1 * ((0.4e1 * t248 + t317) * t390 + t214 * t424 + ((-0.2e1 * t132 * t339 + (-t100 * t107 - t243 / 0.4e1) * t130) * t390 + (t64 * t329 + (-t132 * t432 + (t130 * t68 + t132 * t97) * t122) * qJD(2)) * pkin(1)) * pkin(2)) * t416 - 0.2e1 * ((t125 * t243 / 0.4e1 + t100 * t337 + qJD(2) * t232) * t389 + t214 * t407 + ((t339 * t389 - t65 * t380) * t130 + ((t321 + (-t118 + t338) * t130) * t389 + (-t130 * t70 - t132 * t98) * t380) * qJD(2)) * pkin(2)) * t231;
t3 = (qJD(2) * t179 * t244 + ((t296 * t64 + 0.4e1 * t287) * t56 + (t42 * t266 + (t294 * t64 + t408 * t65) * t79) * t54) * t124) * pkin(4) + t355;
t402 = t3 / 0.2e1;
t387 = Ifges(5,5) * t98;
t385 = Ifges(5,6) * t97;
t381 = pkin(1) * t121;
t378 = pkin(2) * t121;
t377 = pkin(2) * t122;
t374 = pkin(3) * t124;
t34 = t180 * t298;
t370 = t34 * Ifges(5,3);
t41 = ((t67 * t96 + t69 * t99) * t290 - t175) * t141;
t368 = t41 * t74;
t362 = t72 * t97;
t361 = t74 * t99;
t354 = t133 * t48;
t349 = Ifges(3,4) * t130;
t348 = Ifges(3,4) * t132;
t347 = t121 * t71;
t346 = t121 * t73;
t342 = t130 * t99;
t341 = t132 * t96;
t333 = t121 * t126;
t325 = t130 * t132;
t323 = t132 * t135;
t320 = t135 * t143;
t319 = t135 * t144;
t318 = t147 * t135;
t312 = qJD(1) * t130;
t311 = qJD(1) * t132;
t310 = qJD(1) * t137;
t309 = qJD(1) * t140;
t308 = qJD(2) * t122;
t304 = qJDD(2) * t99;
t302 = qJDD(2) * t117;
t301 = qJDD(2) * t118;
t300 = -0.2e1 * t377;
t292 = m(4) * t376;
t282 = pkin(2) * t311;
t279 = t71 * t331;
t277 = -t72 * t98 / 0.2e1;
t276 = t74 * t96 / 0.2e1;
t274 = t129 * t320;
t271 = -t347 / 0.2e1;
t270 = t347 / 0.2e1;
t267 = qJDD(1) * t347;
t265 = t132 * t303;
t264 = pkin(1) * t73 * t377;
t257 = t303 / 0.2e1;
t255 = pkin(2) * t278;
t252 = t96 * t278;
t251 = t99 * t278;
t250 = t41 * t361 / 0.2e1;
t242 = t137 * t271;
t241 = t137 * t270;
t237 = qJD(1) * t271;
t236 = qJD(1) * t270;
t229 = t97 * mrSges(5,1) + t98 * mrSges(5,2);
t228 = t98 * Ifges(5,1) - t97 * Ifges(5,4);
t227 = t98 * Ifges(5,4) - t97 * Ifges(5,2);
t226 = -t385 + t387;
t225 = qJD(2) * t255;
t62 = t130 * t251;
t224 = t132 * t252;
t221 = t349 + t418;
t220 = Ifges(3,5) * t132 - Ifges(3,6) * t130;
t219 = -t341 + t342;
t216 = t122 * t97 - t432;
t215 = 0.2e1 * t255;
t205 = t96 * t225;
t204 = t99 * t225;
t200 = 0.2e1 * t71 * t234;
t199 = -0.2e1 * t71 * t233;
t198 = t71 * t417;
t197 = t129 * t96 + t325 * t99;
t196 = t229 * t71;
t195 = -t342 / 0.2e1 + t341 / 0.2e1;
t194 = -t340 / 0.2e1 - t343 / 0.2e1;
t193 = t276 * t39 - t63 * t73;
t192 = t64 * t71 - t38 * t362 / 0.2e1;
t191 = t277 * t38 + t65 * t71;
t190 = t66 * t73 - t39 * t361 / 0.2e1;
t189 = t276 * t41 - t67 * t73;
t40 = ((t68 * t97 + t70 * t98) * t290 - t174) * t138;
t188 = -t68 * t71 + t40 * t362 / 0.2e1;
t187 = t277 * t40 + t70 * t71;
t186 = t135 * t208;
t185 = t130 * (Ifges(3,1) * t132 - t349);
t183 = t102 ^ 2 * t335 / 0.4e1 + t107 * (t135 * t254 + (qJDD(2) * t130 + t323) * t246);
t182 = 0.2e1 * qJD(2) * t339 + qJDD(2) * t147;
t176 = qJD(2) * t198;
t172 = t74 * t195;
t171 = 0.2e1 * t197;
t170 = -0.2e1 * t129 * t99 + 0.2e1 * t325 * t96;
t169 = t228 * t279;
t168 = t227 * t279;
t167 = t226 * t279;
t166 = t196 * t331;
t165 = -t183 + t318;
t164 = -t118 * t135 + t182;
t163 = t171 * t264;
t162 = t170 * t264;
t161 = qJD(2) * t163;
t160 = qJD(2) * t162;
t159 = (mrSges(4,1) * t195 + mrSges(4,2) * t194) * t368;
t157 = (mrSges(4,1) * t170 - 0.2e1 * mrSges(4,2) * t197) * t264;
t128 = Ifges(3,4) * t311;
t120 = qJDD(1) * t130 + t265;
t114 = Ifges(3,1) * t312 + Ifges(3,5) * qJD(2) + t128;
t113 = Ifges(3,6) * qJD(2) + qJD(1) * t221;
t61 = qJD(1) * t62;
t58 = t131 * t224;
t47 = -qJD(1) * t224 + t61;
t46 = qJD(1) * t48;
t44 = -t356 * t68 + t357 * t70;
t37 = -mrSges(4,1) * t47 - mrSges(4,2) * t46;
t33 = t178 * t374 + qJD(2);
t28 = (t187 * t121 + t199) * t310;
t27 = (t99 * t208 + (-t69 * t73 + t250) * t378) * t305;
t26 = (t188 * t121 + t200) * t310;
t25 = (t189 * t378 + t208 * t96) * t305;
t24 = -t34 * Ifges(5,5) + qJD(1) * t169;
t23 = -t34 * Ifges(5,6) + qJD(1) * t168;
t19 = (t98 * t267 + (qJD(2) * t199 + t121 * t191) * qJD(1)) * t137;
t18 = (t99 * t186 + (-qJD(2) * t190 - t304 * t73) * t378) * t140;
t17 = (-t97 * t267 + (qJD(2) * t200 - t121 * t192) * qJD(1)) * t137;
t16 = (t96 * t186 + (-t96 * t73 * qJDD(2) + qJD(2) * t193) * t378) * t140;
t15 = -t46 * Ifges(4,1) + t47 * Ifges(4,4) + t33 * Ifges(4,5);
t14 = -t46 * Ifges(4,4) + t47 * Ifges(4,2) + t33 * Ifges(4,6);
t12 = (t162 + (t41 * t172 + (t352 * t130 - t353 * t132) * t73) * t121) * t309;
t11 = t61 + (t163 + (t189 * t130 + (-t352 * t73 + t250) * t132) * t121) * t309;
t10 = (t160 + (t39 * t172 + (qJD(2) * t218 + t130 * t66 - t132 * t63) * t73) * t121) * t140;
t9 = (t161 + (-t194 * t369 + (qJD(2) * t219 - t130 * t63 - t132 * t66) * t73) * t121) * t140;
t8 = ((-t119 * t96 + t120 * t99) * t346 + (t160 + (t130 * t190 + t132 * t193) * t121) * qJD(1)) * t140;
t7 = ((-t119 * t99 - t120 * t96) * t346 + (t161 + (t130 * t193 - t132 * t190) * t121) * qJD(1)) * t140;
t6 = 0.2e1 * pkin(3) * t256 * t178 + qJDD(2) + ((t89 * t81 * t63 + t88 * t367) * t51 * t406 + (-t53 * t63 - ((0.8e1 * t328 * t96 + 0.4e1 * t378) * t140 * t129 * t319 + ((-t182 * t121 + (t121 * t262 + t300 * t96) * t135) * t132 + ((t165 + t302) * t121 + (-0.4e1 * t63 * t308 + (-t122 * t96 - t333) * t422) * pkin(2)) * t130) * t379) * t99 - t51 * t66) * t81 * t88 + (((-t183 * t333 + (0.8e1 * t99 * t123 * t274 + (t129 * t422 + 0.6e1 * t130 * t323) * t378) * t144) * t140 + (((t302 + t318) * t121 + t99 * t135 * t300) * t132 + ((-t117 * t135 + t182) * t121 + (-0.4e1 * qJD(2) * t66 - 0.2e1 * t304) * t377) * t130) * t379) * t81 + t53 * t431) * t87) * t374;
t5 = -0.2e1 * ((t125 * t183 + t135 * t232) * t389 + (t123 * t319 * t407 + qJDD(2) * t381) * t327 + (((t301 + t318) * t132 + t164 * t130) * t389 + (-t98 * t323 + (-0.2e1 * qJD(2) * t65 - qJDD(2) * t98) * t130) * t380) * pkin(2)) * t231 - 0.2e1 * ((t123 * t144 * t424 - 0.2e1 * t381) * t274 + ((pkin(1) * t135 * t216 + t164 * t389) * t132 + ((t165 + t301) * t390 + (qJDD(2) * t216 + 0.2e1 * t308 * t64) * pkin(1)) * t130) * pkin(2)) * t416 + 0.2e1 * (-t50 * t65 + t52 * t64) * t258 - 0.4e1 * pkin(4) * t256 * t180 + 0.2e1 * (t52 * t287 + (t92 * t80 * t42 + t93 * t79 * t64) * t50 * t98) * t298;
t1 = (-t180 * t244 + ((t296 * t68 + 0.4e1 * t360 * t44) * t52 + (t44 * t266 + (t294 * t68 + t408 * t70) * t79) * t50) * t124) * pkin(4) + t355;
t2 = [(-t58 * mrSges(4,1) + t411 * t133 + (-t166 - (-mrSges(4,1) * t342 - mrSges(4,2) * t218) * t278 + t410) * t131) * g(1) + (-t354 * mrSges(4,2) + (-mrSges(4,1) * t219 * t278 + t229 * t279 - t410) * t133 + t411 * t131) * g(2) + (mrSges(5,1) * t17 - mrSges(5,2) * t19 - qJDD(1) * t166 - (-t229 * t176 + (mrSges(5,1) * t192 + mrSges(5,2) * t191) * t121) * t310) * pkin(1) + (m(5) * t144 + Ifges(2,3)) * qJDD(1) + (-t10 * t204 + t205 * t9) * mrSges(4,3) + (t132 * t348 + t185) * t257 + t5 * t167 / 0.2e1 + t33 * (Ifges(4,5) * t9 + Ifges(4,6) * t10) / 0.2e1 + t9 * t15 / 0.2e1 + t19 * t169 / 0.2e1 + t10 * t14 / 0.2e1 - (-mrSges(4,1) * t8 + mrSges(4,2) * t7) * t375 + t37 * t281 + t120 * t348 / 0.2e1 + t132 * (Ifges(3,4) * t120 + Ifges(3,2) * t119) / 0.2e1 + (-m(4) * t143 * t265 + Ifges(3,1) * t120 + Ifges(3,4) * t420 - t257 * t418) * t130 + t135 * t220 / 0.2e1 + t65 * t24 * t241 + t64 * t23 * t242 + ((Ifges(5,4) * t19 + Ifges(5,2) * t17 + Ifges(5,6) * t5) * t242 + t405 * (-t227 * t176 + (Ifges(5,4) * t191 - Ifges(5,2) * t192) * t121) * t237 + t412 * t23) * t97 + (t405 * (-t228 * t176 + (Ifges(5,1) * t191 - Ifges(5,4) * t192) * t121) * t236 + (Ifges(5,1) * t19 + Ifges(5,4) * t17 + Ifges(5,5) * t5) * t241 - t412 * t24) * t98 - (-mrSges(4,1) * t10 + mrSges(4,2) * t9) * t282 + qJDD(2) * (Ifges(3,5) * t130 + Ifges(3,6) * t132) - t46 * (Ifges(4,1) * t9 + Ifges(4,4) * t10) / 0.2e1 + t47 * (Ifges(4,4) * t9 + Ifges(4,2) * t10) / 0.2e1 - t34 * (-t226 * t176 + (Ifges(5,5) * t191 - Ifges(5,6) * t192) * t121) * t137 / 0.2e1 + t17 * t168 / 0.2e1 + t114 * t306 / 0.2e1 - t113 * t307 / 0.2e1 - (-mrSges(4,2) * t419 - mrSges(4,3) * t16 + Ifges(4,1) * t7 + Ifges(4,4) * t8 + Ifges(4,5) * t6) * t48 + (mrSges(4,1) * t419 + mrSges(4,3) * t18 + Ifges(4,4) * t7 + Ifges(4,2) * t8 + Ifges(4,6) * t6) * (t62 - t224) + m(4) * t119 * t322 + t221 * t420; (-t413 * (t196 * t417 + (mrSges(5,1) * t188 - mrSges(5,2) * t187) * t121) + (t97 * (Ifges(5,4) * t28 + Ifges(5,2) * t26 + Ifges(5,6) * t1) + t3 * t387) * t236 + (t98 * (Ifges(5,1) * t28 + Ifges(5,4) * t26 + Ifges(5,5) * t1) + t3 * t385) * t237 + (-(-mrSges(5,1) * t98 + mrSges(5,2) * t97) * t198 - (mrSges(5,1) * t187 + mrSges(5,2) * t188) * t121) * g(3)) * t137 + (-t62 * mrSges(4,1) - ((mrSges(4,1) * t171 + mrSges(4,2) * t170) * t264 + ((-mrSges(4,1) * t194 + mrSges(4,2) * t195) * t368 + ((-t67 * mrSges(4,1) + mrSges(4,2) * t352) * t130 - t409 * t132) * t73) * t121) * t140 + t414) * g(3) + (-t1 / 0.2e1 + t402) * (qJD(1) * t167 - t370) + (t16 * mrSges(4,1) - t18 * mrSges(4,2) + Ifges(4,5) * t7 + Ifges(4,6) * t8 + Ifges(4,3) * t6) * (-t177 * t374 + 0.1e1) + (-t255 * t63 + t415 * t96 - t25) * (mrSges(4,1) * t33 + mrSges(4,3) * t46) - g(2) * (t58 * mrSges(4,2) + (-t292 + (t157 + (t159 + ((-mrSges(4,1) * t353 + t69 * mrSges(4,2)) * t132 + t409 * t130) * t73) * t121) * t140) * t131) - t47 * (Ifges(4,4) * t11 + Ifges(4,2) * t12) / 0.2e1 + (-(-mrSges(4,2) * t6 + mrSges(4,3) * t8) * t251 - (mrSges(4,1) * t6 - mrSges(4,3) * t7) * t252 - t37 * t312) * pkin(2) + t421 * (t46 * Ifges(4,5) - t47 * Ifges(4,6) - t33 * Ifges(4,3)) - (-Ifges(3,2) * t312 + t114 + t128) * t311 / 0.2e1 + (t403 * t272 - (-t25 * t96 - t27 * t99) * qJD(2) * t215 / 0.2e1 + (-t16 * t96 - t18 * t99) * t215 / 0.2e1 + (-t75 * t175 * t320 + (t75 * t181 - t351 / t78 ^ 2 * t39) * t143 * t308) * t141 / 0.2e1) * m(4) - t220 * t303 / 0.2e1 + Ifges(3,6) * t119 + Ifges(3,5) * t120 - t33 * (Ifges(4,5) * t11 + Ifges(4,6) * t12) / 0.2e1 - t28 * t24 / 0.2e1 - t11 * t15 / 0.2e1 - t12 * t14 / 0.2e1 - t26 * t23 / 0.2e1 + t113 * t312 / 0.2e1 + (-mrSges(4,1) * t12 + mrSges(4,2) * t11) * t282 + (t421 * mrSges(4,1) - mrSges(4,3) * t11) * t205 + (-t421 * mrSges(4,2) + mrSges(4,3) * t12) * t204 - t403 * t185 / 0.2e1 + (-t255 * t66 + t415 * t99 - t27) * (-mrSges(4,2) * t33 + mrSges(4,3) * t47) + (Ifges(5,5) * t19 + Ifges(5,6) * t17 + Ifges(5,3) * t5) * t179 * t298 + t413 * (mrSges(3,1) * t130 + mrSges(3,2) * t132) - t370 * t402 + t46 * (Ifges(4,1) * t11 + Ifges(4,4) * t12) / 0.2e1 - g(1) * (t354 * mrSges(4,1) + (-t292 + (t157 + (t159 + ((t130 * t69 - t132 * t67) * mrSges(4,1) + (t130 * t67 + t132 * t69 - t219) * mrSges(4,2)) * t73) * t121) * t140) * t133) + t34 * (Ifges(5,5) * t28 + Ifges(5,6) * t26 + Ifges(5,3) * t1) / 0.2e1 + qJD(1) * pkin(1) * (-mrSges(5,1) * t26 + mrSges(5,2) * t28) + (Ifges(3,3) * qJDD(2));];
tau = t2;
