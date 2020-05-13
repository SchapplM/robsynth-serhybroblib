% Calculate vector of inverse dynamics joint torques for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnDE2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:30
% EndTime: 2020-04-12 19:34:18
% DurationCPUTime: 20.68s
% Computational Cost: add. (149143->419), mult. (214385->909), div. (7926->21), fcn. (57649->13), ass. (0->323)
t420 = 2 * pkin(2);
t142 = pkin(2) ^ 2;
t143 = pkin(1) ^ 2;
t131 = cos(qJ(2));
t359 = pkin(2) * t131;
t306 = -0.2e1 * pkin(1) * t359 + t143;
t123 = t142 + t306;
t120 = 0.1e1 / t123;
t377 = -t120 / 0.2e1;
t124 = pkin(1) - t359;
t419 = t120 * t124;
t410 = pkin(3) ^ 2;
t411 = pkin(4) ^ 2;
t305 = t410 - t411;
t116 = t123 + t305;
t125 = pkin(1) * t131 - pkin(2);
t129 = sin(qJ(2));
t383 = (-pkin(3) - pkin(4));
t114 = ((pkin(2) - t383) * (pkin(2) + t383)) + t306;
t382 = (pkin(4) - pkin(3));
t115 = ((pkin(2) - t382) * (pkin(2) + t382)) + t306;
t326 = t114 * t115;
t146 = sqrt(-t326);
t316 = t129 * t146;
t95 = -pkin(1) * t316 - t116 * t125;
t85 = t95 ^ 2;
t86 = 0.1e1 / t95;
t88 = t86 / t85;
t366 = pkin(1) * t116;
t111 = t129 * t366;
t98 = -t125 * t146 + t111;
t94 = t98 ^ 2;
t346 = t88 * t94;
t87 = 0.1e1 / t95 ^ 2;
t347 = t87 * t98;
t291 = -0.2e1 * pkin(2) * t125;
t255 = t116 + t291;
t416 = pkin(1) * pkin(2);
t242 = (-t114 - t115) * t416;
t102 = t129 * t242;
t101 = qJD(2) * t102;
t106 = 0.1e1 / t146;
t331 = t101 * t106;
t266 = t129 * t331;
t313 = t131 * t146;
t59 = (-t266 + (t129 * t255 - t313) * qJD(2)) * pkin(1);
t299 = qJD(2) * t129;
t261 = t146 * t299;
t128 = t129 ^ 2;
t318 = t128 * t143;
t262 = qJD(2) * t318;
t298 = qJD(2) * t131;
t328 = t106 * t125;
t62 = pkin(1) * t261 - t101 * t328 + t262 * t420 + t298 * t366;
t83 = t87 * t94 + 0.1e1;
t81 = 0.1e1 / t83 ^ 2;
t353 = (-t346 * t59 + t347 * t62) * t81;
t418 = -0.2e1 * t353;
t376 = t120 / 0.2e1;
t139 = 0.1e1 / pkin(3);
t121 = 0.1e1 / t123 ^ 2;
t274 = pkin(2) * t299;
t249 = pkin(1) * t274;
t225 = t121 * t249;
t203 = -0.2e1 * t225;
t51 = (t120 * t59 + t203 * t95) * t139;
t53 = (t120 * t62 + t203 * t98) * t139;
t412 = -t51 * t347 + t53 * t86;
t80 = 0.1e1 / t83;
t174 = t412 * t80;
t321 = t121 * t129;
t277 = pkin(1) * t321;
t254 = pkin(2) * t277;
t235 = -0.2e1 * t254;
t330 = t106 * t102;
t63 = t111 + (-t313 + (t291 - t330) * t129) * pkin(1);
t55 = (t120 * t63 + t235 * t95) * t139;
t65 = -t102 * t328 + t318 * t420 + (t116 * t131 + t316) * pkin(1);
t57 = (t120 * t65 + t235 * t98) * t139;
t200 = t347 * t55 - t57 * t86;
t175 = t80 * t200;
t360 = pkin(2) * t129;
t284 = pkin(1) * t360;
t241 = 0.2e1 * t284;
t393 = 0.2e1 * t98;
t406 = ((qJD(2) * t175 + t174) * t241 + (-0.2e1 * t412 * t81 * (-t346 * t63 + t347 * t65) + t200 * t418 + ((-t51 * t65 - t53 * t63 + t55 * t62 + t57 * t59) * t87 + (t51 * t63 - t55 * t59) * t88 * t393) * t80) * t123) * pkin(3);
t117 = t123 - t305;
t96 = -pkin(2) * t316 + t117 * t124;
t91 = 0.1e1 / t96 ^ 2;
t110 = t117 * t360;
t97 = t124 * t146 + t110;
t345 = t91 * t97;
t136 = 0.1e1 / pkin(4);
t289 = 0.2e1 * t124 * pkin(1);
t60 = (-t266 + (-t313 + (t117 + t289) * t129) * qJD(2)) * pkin(2);
t50 = (t225 * t96 + t377 * t60) * t136;
t319 = t128 * t142;
t276 = pkin(1) * t319;
t244 = qJD(2) * t276;
t309 = (t117 * t298 + t261) * pkin(2);
t329 = t106 * t124;
t61 = t101 * t329 + 0.2e1 * t244 + t309;
t52 = (-t225 * t97 + t376 * t61) * t136;
t93 = t97 ^ 2;
t82 = t91 * t93 + 0.1e1;
t78 = 0.1e1 / t82;
t90 = 0.1e1 / t96;
t177 = t78 * (t50 * t345 + t52 * t90);
t357 = pkin(4) * t123;
t290 = 0.2e1 * t357;
t414 = -m(4) * t359 - mrSges(3,1) * t131 + mrSges(3,2) * t129;
t409 = -0.4e1 * t96;
t407 = 2 * qJDD(2);
t295 = qJD(1) * qJD(2);
t118 = qJDD(1) * t131 - t129 * t295;
t405 = t118 / 0.2e1;
t363 = pkin(1) * t139;
t404 = pkin(2) * t118;
t403 = t131 * Ifges(3,2);
t402 = t121 * t241;
t281 = t78 * t357;
t251 = t91 * t281;
t401 = t136 * t97 * t251;
t140 = 0.1e1 / t410;
t342 = t85 + t94;
t77 = t342 * t140 * t121;
t72 = t77 ^ (-0.1e1 / 0.2e1);
t202 = 0.2e1 * t142 * t72 * t277;
t297 = qJD(2) * t139;
t322 = t120 * t139;
t122 = t120 * t121;
t240 = 0.4e1 * t284;
t211 = t122 * t240;
t172 = t342 * t211;
t178 = 0.2e1 * t59 * t95 + 0.2e1 * t62 * t98;
t39 = (-qJD(2) * t172 + t121 * t178) * t140;
t74 = 0.1e1 / t77;
t73 = t72 * t74;
t354 = t39 * t73;
t400 = pkin(2) * t322 * t354 / 0.2e1 + t202 * t297;
t323 = t120 * t136;
t137 = 0.1e1 / t411;
t89 = t96 ^ 2;
t341 = t89 + t93;
t171 = t341 * t211;
t283 = 0.2e1 * t121;
t38 = ((t60 * t96 + t61 * t97) * t283 - qJD(2) * t171) * t137;
t76 = t341 * t137 * t121;
t70 = t76 ^ (-0.1e1 / 0.2e1);
t71 = t70 / t76;
t399 = t136 * t70 * t225 + t38 * t71 * t323 / 0.4e1;
t397 = -mrSges(4,3) - mrSges(5,3) - mrSges(3,3) + mrSges(2,2);
t69 = qJ(2) + atan2(t98 * t322, t95 * t322);
t67 = sin(t69);
t68 = cos(t69);
t224 = -mrSges(4,1) * t68 + mrSges(4,2) * t67;
t396 = m(5) * pkin(1) + mrSges(2,1) + t224 - t414;
t395 = -0.2e1 * t91;
t394 = 0.4e1 * t97;
t392 = t136 ^ 2;
t390 = qJD(1) ^ 2;
t134 = qJD(2) ^ 2;
t229 = t96 * t254;
t64 = t110 + (-t313 + (t289 - t330) * t129) * pkin(2);
t54 = (t377 * t64 + t229) * t136;
t228 = t97 * t254;
t66 = t102 * t329 + 0.2e1 * t276 + (t117 * t131 + t316) * pkin(2);
t56 = (t376 * t66 - t228) * t136;
t176 = t78 * (-t345 * t54 - t56 * t90);
t79 = 0.1e1 / t82 ^ 2;
t259 = 0.4e1 * t79 * t345;
t348 = t79 * t90;
t92 = t90 / t89;
t344 = t92 * t93;
t42 = -t344 * t60 + t345 * t61;
t280 = t42 * t348;
t286 = t92 * t394;
t288 = 0.2e1 * t78 * t91;
t320 = t122 * t142;
t208 = t262 * t320;
t226 = t136 * t90 * t281;
t314 = t131 * t142;
t265 = t129 * t314;
t227 = 0.6e1 * pkin(1) * t265;
t327 = 0.4e1 * t106 / t326;
t238 = t101 * t102 * t327;
t364 = pkin(1) * t121;
t247 = -0.4e1 * t142 * t318;
t99 = (t131 * t242 + t247) * qJD(2);
t343 = -0.2e1 * ((0.4e1 * t244 + t309) * t377 + t208 * t409 + ((-0.2e1 * t131 * t331 + (-t106 * t99 - t238 / 0.4e1) * t129) * t377 + (t60 * t321 + (-t131 * t419 + (t129 * t64 + t131 * t96) * t121) * qJD(2)) * pkin(1)) * pkin(2)) * t401 - 0.2e1 * ((t124 * t238 / 0.4e1 + t99 * t329 + qJD(2) * t227) * t376 + t208 * t394 + ((t331 * t376 - t61 * t364) * t129 + ((t313 + (-t117 + t330) * t129) * t376 + (-t129 * t66 - t131 * t97) * t364) * qJD(2)) * pkin(2)) * t226;
t3 = (qJD(2) * t176 * t240 + ((t288 * t60 + 0.4e1 * t280) * t56 + (t42 * t259 + (t286 * t60 + t395 * t61) * t78) * t54) * t123) * pkin(4) + t343;
t389 = t3 / 0.2e1;
t373 = Ifges(5,5) * t97;
t370 = Ifges(5,6) * t96;
t365 = pkin(1) * t120;
t362 = pkin(2) * t120;
t361 = pkin(2) * t121;
t358 = pkin(3) * t123;
t34 = t177 * t290;
t355 = t34 * Ifges(5,3);
t350 = t71 * t96;
t349 = t73 * t98;
t340 = -t95 - t65;
t339 = Ifges(3,4) * t129;
t338 = Ifges(3,4) * t131;
t337 = t120 * t70;
t336 = t120 * t72;
t335 = t129 * t98;
t334 = t131 * t95;
t333 = t95 * t129;
t332 = t98 * t131;
t325 = t120 * t125;
t317 = t129 * t131;
t315 = t131 * t134;
t312 = t134 * t142;
t311 = t134 * t143;
t310 = t146 * t134;
t304 = qJD(1) * t129;
t303 = qJD(1) * t131;
t302 = qJD(1) * t136;
t301 = qJD(1) * t139;
t300 = qJD(2) * t121;
t296 = qJDD(2) * t98;
t294 = qJDD(2) * t116;
t293 = qJDD(2) * t117;
t292 = -0.2e1 * t361;
t275 = pkin(2) * t303;
t272 = t70 * t323;
t271 = t72 * t322;
t270 = -t71 * t97 / 0.2e1;
t269 = t73 * t95 / 0.2e1;
t267 = t128 * t312;
t264 = -t337 / 0.2e1;
t263 = t337 / 0.2e1;
t260 = qJDD(1) * t337;
t257 = t131 * t295;
t250 = t295 / 0.2e1;
t248 = pkin(2) * t271;
t41 = ((t63 * t95 + t65 * t98) * t283 - t172) * t140;
t245 = t41 * t349 / 0.2e1;
t239 = qJD(1) * t271;
t237 = t136 * t264;
t236 = t136 * t263;
t232 = qJD(1) * t264;
t231 = qJD(1) * t263;
t223 = mrSges(5,1) * t96 + mrSges(5,2) * t97;
t222 = Ifges(5,1) * t97 - Ifges(5,4) * t96;
t221 = t97 * Ifges(5,4) - t96 * Ifges(5,2);
t220 = -t370 + t373;
t219 = qJD(2) * t248;
t218 = t72 * t283 * t416;
t215 = t339 + t403;
t214 = Ifges(3,5) * t131 - Ifges(3,6) * t129;
t213 = -t334 + t335;
t212 = t332 + t333;
t210 = t121 * t96 - t419;
t209 = 0.2e1 * t248;
t199 = t95 * t219;
t198 = t98 * t219;
t194 = 0.2e1 * t70 * t229;
t193 = -0.2e1 * t70 * t228;
t192 = t70 * t402;
t191 = t223 * t70;
t190 = t269 * t39 - t59 * t72;
t189 = t60 * t70 - t38 * t350 / 0.2e1;
t188 = t270 * t38 + t61 * t70;
t187 = t62 * t72 - t39 * t349 / 0.2e1;
t186 = t269 * t41 - t63 * t72;
t40 = ((t64 * t96 + t66 * t97) * t283 - t171) * t137;
t185 = -t64 * t70 + t40 * t350 / 0.2e1;
t184 = t270 * t40 + t66 * t70;
t183 = t134 * t202;
t182 = t129 * (Ifges(3,1) * t131 - t339);
t180 = t101 ^ 2 * t327 / 0.4e1 + t106 * (t134 * t247 + (qJDD(2) * t129 + t315) * t242);
t179 = 0.2e1 * qJD(2) * t331 + qJDD(2) * t146;
t173 = qJD(2) * t192;
t169 = t73 * (-t335 / 0.2e1 + t334 / 0.2e1);
t167 = t222 * t272;
t166 = t221 * t272;
t165 = t220 * t272;
t164 = t191 * t323;
t48 = t212 * t271;
t163 = -t180 + t310;
t162 = -t117 * t134 + t179;
t160 = (t128 * t95 + t317 * t98) * t218;
t159 = (-t128 * t98 + t317 * t95) * t218;
t158 = qJD(2) * t160;
t157 = qJD(2) * t159;
t132 = cos(qJ(1));
t130 = sin(qJ(1));
t127 = Ifges(3,4) * t303;
t119 = qJDD(1) * t129 + t257;
t113 = Ifges(3,1) * t304 + Ifges(3,5) * qJD(2) + t127;
t112 = Ifges(3,6) * qJD(2) + qJD(1) * t215;
t58 = t239 * t335;
t47 = -t239 * t334 + t58;
t46 = qJD(1) * t48;
t44 = -t344 * t64 + t345 * t66;
t37 = -mrSges(4,1) * t47 - mrSges(4,2) * t46;
t36 = -t175 * t358 + 0.1e1;
t33 = t174 * t358 + qJD(2);
t28 = (t184 * t120 + t193) * t302;
t27 = (t98 * t202 + (-t65 * t72 + t245) * t362) * t297;
t26 = (t185 * t120 + t194) * t302;
t25 = (t186 * t362 + t202 * t95) * t297;
t24 = -t34 * Ifges(5,5) + qJD(1) * t167;
t23 = -t34 * Ifges(5,6) + qJD(1) * t166;
t19 = (t97 * t260 + (qJD(2) * t193 + t120 * t188) * qJD(1)) * t136;
t18 = (t98 * t183 + (-qJD(2) * t187 - t296 * t72) * t362) * t139;
t17 = (-t96 * t260 + (qJD(2) * t194 - t120 * t189) * qJD(1)) * t136;
t16 = (t95 * t183 + (-t95 * t72 * qJDD(2) + qJD(2) * t190) * t362) * t139;
t15 = -t46 * Ifges(4,1) + t47 * Ifges(4,4) + t33 * Ifges(4,5);
t14 = -Ifges(4,4) * t46 + Ifges(4,2) * t47 + Ifges(4,6) * t33;
t12 = (t159 + (t41 * t169 + ((-t63 + t98) * t131 - t340 * t129) * t72) * t120) * t301;
t11 = t58 + (t160 + (t186 * t129 + (t340 * t72 + t245) * t131) * t120) * t301;
t10 = (t157 + (t39 * t169 + (qJD(2) * t212 + t62 * t129 - t59 * t131) * t72) * t120) * t139;
t9 = (t158 + ((t333 / 0.2e1 + t332 / 0.2e1) * t354 + (qJD(2) * t213 - t59 * t129 - t62 * t131) * t72) * t120) * t139;
t8 = ((-t118 * t95 + t119 * t98) * t336 + (t157 + (t129 * t187 + t131 * t190) * t120) * qJD(1)) * t139;
t7 = ((-t118 * t98 - t119 * t95) * t336 + (t158 + (t129 * t190 - t131 * t187) * t120) * qJD(1)) * t139;
t6 = 0.2e1 * pkin(3) * t249 * t174 + qJDD(2) + ((t59 * t80 * t88 + t353 * t87) * t51 * t393 + (-t53 * t59 - ((0.8e1 * t320 * t95 + 0.4e1 * t362) * t139 * t128 * t311 + ((-t179 * t120 + (t120 * t255 + t292 * t95) * t134) * t131 + ((t163 + t294) * t120 + (-0.4e1 * t59 * t300 + (-t121 * t95 - t325) * t407) * pkin(2)) * t129) * t363) * t98 - t51 * t62) * t80 * t87 + (((-t180 * t325 + (0.8e1 * t98 * t122 * t267 + (t128 * t407 + 0.6e1 * t129 * t315) * t362) * t143) * t139 + (((t294 + t310) * t120 + t98 * t134 * t292) * t131 + ((-t116 * t134 + t179) * t120 + (-0.4e1 * qJD(2) * t62 - 0.2e1 * t296) * t361) * t129) * t363) * t80 + t53 * t418) * t86) * t358;
t5 = -0.2e1 * ((t180 * t124 + t134 * t227) * t376 + (t122 * t311 * t394 + qJDD(2) * t365) * t319 + (((t293 + t310) * t131 + t162 * t129) * t376 + (-t97 * t315 + (-0.2e1 * qJD(2) * t61 - qJDD(2) * t97) * t129) * t364) * pkin(2)) * t226 - 0.2e1 * ((t122 * t143 * t409 - 0.2e1 * t365) * t267 + ((pkin(1) * t134 * t210 + t162 * t376) * t131 + ((t163 + t293) * t377 + (qJDD(2) * t210 + 0.2e1 * t300 * t60) * pkin(1)) * t129) * pkin(2)) * t401 + 0.2e1 * (-t50 * t61 + t52 * t60) * t251 - 0.4e1 * pkin(4) * t249 * t177 + 0.2e1 * (t52 * t280 + (t42 * t79 * t91 + t60 * t78 * t92) * t50 * t97) * t290;
t1 = (-t177 * t240 + ((t288 * t64 + 0.4e1 * t348 * t44) * t52 + (t44 * t259 + (t286 * t64 + t395 * t66) * t78) * t50) * t123) * pkin(4) + t343;
t2 = [(-t10 * t198 + t199 * t9) * mrSges(4,3) + (-m(4) * t142 * t257 + Ifges(3,1) * t119 + Ifges(3,4) * t405 - t250 * t403) * t129 - (-mrSges(4,2) * t404 - mrSges(4,3) * t16 + Ifges(4,1) * t7 + Ifges(4,4) * t8 + Ifges(4,5) * t6) * t48 + (t392 * (-t221 * t173 + (Ifges(5,4) * t188 - Ifges(5,2) * t189) * t120) * t232 + (Ifges(5,4) * t19 + Ifges(5,2) * t17 + Ifges(5,6) * t5) * t237 + t399 * t23) * t96 + (t392 * (-t222 * t173 + (Ifges(5,1) * t188 - Ifges(5,4) * t189) * t120) * t231 + (Ifges(5,1) * t19 + Ifges(5,4) * t17 + Ifges(5,5) * t5) * t236 - t399 * t24) * t97 + t131 * (Ifges(3,4) * t119 + Ifges(3,2) * t118) / 0.2e1 - t34 * (-t220 * t173 + (Ifges(5,5) * t188 - Ifges(5,6) * t189) * t120) * t136 / 0.2e1 - (-mrSges(4,1) * t8 + mrSges(4,2) * t7) * t359 + ((t223 * t272 - t396) * t132 + t397 * t130) * g(2) + (t397 * t132 + (-t164 + t396) * t130) * g(1) + t17 * t166 / 0.2e1 + t5 * t165 / 0.2e1 + t134 * t214 / 0.2e1 + t119 * t338 / 0.2e1 + t113 * t298 / 0.2e1 + m(4) * t118 * t314 - t112 * t299 / 0.2e1 + (m(5) * t143 + Ifges(2,3)) * qJDD(1) + (-(-t223 * t173 + (mrSges(5,1) * t189 + mrSges(5,2) * t188) * t120) * t302 - qJDD(1) * t164 + mrSges(5,1) * t17 - mrSges(5,2) * t19) * pkin(1) + t19 * t167 / 0.2e1 + t215 * t405 + (mrSges(4,1) * t404 + mrSges(4,3) * t18 + Ifges(4,4) * t7 + Ifges(4,2) * t8 + Ifges(4,6) * t6) * t213 * t271 + t33 * (Ifges(4,5) * t9 + Ifges(4,6) * t10) / 0.2e1 + t9 * t15 / 0.2e1 + t10 * t14 / 0.2e1 - t46 * (Ifges(4,1) * t9 + Ifges(4,4) * t10) / 0.2e1 + t47 * (Ifges(4,4) * t9 + Ifges(4,2) * t10) / 0.2e1 + qJDD(2) * (Ifges(3,5) * t129 + Ifges(3,6) * t131) - (-mrSges(4,1) * t10 + mrSges(4,2) * t9) * t275 + t37 * t274 + (t131 * t338 + t182) * t250 + t61 * t24 * t236 + t60 * t23 * t237; (g(1) * t132 + g(2) * t130) * (-t136 * (t191 * t402 + (mrSges(5,1) * t185 - mrSges(5,2) * t184) * t120) + mrSges(3,1) * t129 + mrSges(3,2) * t131 + m(4) * t360 - (mrSges(4,1) * t67 + mrSges(4,2) * t68) * t36) + (-t248 * t59 + t400 * t95 - t25) * (mrSges(4,1) * t33 + mrSges(4,3) * t46) + t406 * (Ifges(4,5) * t46 - Ifges(4,6) * t47 - Ifges(4,3) * t33) - (-Ifges(3,2) * t304 + t113 + t127) * t303 / 0.2e1 + ((-t16 * t95 - t18 * t98) * t209 / 0.2e1 + (-t74 * t172 * t312 + (t74 * t178 - t342 / t77 ^ 2 * t39) * t142 * t300) * t140 / 0.2e1 + t390 * t265 - (-t25 * t95 - t27 * t98) * qJD(2) * t209 / 0.2e1) * m(4) + (mrSges(4,1) * t406 - mrSges(4,3) * t11) * t199 + (-mrSges(4,2) * t406 + mrSges(4,3) * t12) * t198 + (t16 * mrSges(4,1) - t18 * mrSges(4,2)) * t36 + (-t248 * t62 + t400 * t98 - t27) * (-mrSges(4,2) * t33 + mrSges(4,3) * t47) + (-t98 * (-mrSges(4,2) * t6 + mrSges(4,3) * t8) - t95 * (mrSges(4,1) * t6 - mrSges(4,3) * t7)) * t248 + (t389 - t1 / 0.2e1) * (qJD(1) * t165 - t355) + t46 * (Ifges(4,1) * t11 + Ifges(4,4) * t12) / 0.2e1 + (-pkin(2) * t37 + t112 / 0.2e1) * t304 - t47 * (Ifges(4,4) * t11 + Ifges(4,2) * t12) / 0.2e1 + ((t3 * t373 + t96 * (Ifges(5,4) * t28 + Ifges(5,2) * t26 + Ifges(5,6) * t1)) * t231 + (t3 * t370 + t97 * (Ifges(5,1) * t28 + Ifges(5,4) * t26 + Ifges(5,5) * t1)) * t232 + (-(-mrSges(5,1) * t97 + mrSges(5,2) * t96) * t192 - (mrSges(5,1) * t184 + mrSges(5,2) * t185) * t120) * g(3)) * t136 - t390 * t182 / 0.2e1 + (Ifges(5,5) * t19 + Ifges(5,6) * t17 + Ifges(5,3) * t5) * t176 * t290 - t214 * t295 / 0.2e1 + t34 * (Ifges(5,5) * t28 + Ifges(5,6) * t26 + Ifges(5,3) * t1) / 0.2e1 - t355 * t389 - t28 * t24 / 0.2e1 - t26 * t23 / 0.2e1 - t11 * t15 / 0.2e1 - t12 * t14 / 0.2e1 + (Ifges(3,3) * qJDD(2)) + Ifges(3,6) * t118 + Ifges(3,5) * t119 + (-t224 * t36 + t414) * g(3) - t33 * (Ifges(4,5) * t11 + Ifges(4,6) * t12) / 0.2e1 + t36 * (Ifges(4,5) * t7 + Ifges(4,6) * t8 + Ifges(4,3) * t6) + (-mrSges(4,1) * t12 + mrSges(4,2) * t11) * t275 + qJD(1) * pkin(1) * (-mrSges(5,1) * t26 + mrSges(5,2) * t28);];
tau = t2;
