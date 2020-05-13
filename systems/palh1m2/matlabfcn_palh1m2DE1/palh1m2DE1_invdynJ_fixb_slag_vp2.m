% Calculate vector of inverse dynamics joint torques for
% palh1m2DE1
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
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m2DE1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:38
% EndTime: 2020-05-01 20:58:39
% DurationCPUTime: 28.63s
% Computational Cost: add. (11598->626), mult. (21976->941), div. (0->0), fcn. (25894->24), ass. (0->313)
t269 = cos(qJ(3));
t373 = qJD(2) * t269;
t358 = pkin(1) * t373;
t251 = pkin(22) + pkin(21);
t243 = sin(t251);
t244 = cos(t251);
t256 = sin(pkin(20));
t260 = cos(pkin(20));
t266 = sin(pkin(18));
t272 = cos(pkin(18));
t206 = -t272 * t256 + t260 * t266;
t209 = t256 * t266 + t260 * t272;
t263 = sin(qJ(3));
t135 = t206 * t269 - t209 * t263;
t270 = cos(qJ(2));
t117 = t135 * t270;
t264 = sin(qJ(2));
t306 = t263 * t206 + t209 * t269;
t345 = -t306 * t264 + t117;
t518 = t135 * t264 + t270 * t306;
t515 = t345 * t243 + t244 * t518;
t252 = qJD(3) + qJD(2);
t421 = pkin(1) * qJD(2);
t359 = t263 * t421;
t304 = pkin(5) * t252 + t359;
t318 = -t518 * t243 + t244 * t345;
t523 = t304 * t318;
t18 = t515 * t358 - t523;
t127 = t135 * qJD(3);
t122 = t306 * qJD(3);
t519 = qJD(2) * t518 + t122 * t270 + t127 * t264;
t521 = -qJD(2) * t345 + t122 * t264;
t24 = t244 * t519 + (t127 * t270 - t521) * t243;
t485 = Ifges(3,4) + Ifges(10,4);
t520 = t515 * t304;
t514 = Ifges(3,1) + Ifges(10,1);
t513 = Ifges(3,2) + Ifges(10,2);
t239 = pkin(1) * t263 + pkin(5);
t25 = (-qJD(3) * t117 + t521) * t244 + t519 * t243;
t302 = pkin(1) * t318;
t433 = pkin(1) * t269;
t362 = t515 * t433;
t434 = -t24 * t433 - t25 * t239 + (-t263 * t302 + t362) * qJD(3) + t18;
t517 = t485 * t270;
t516 = t485 * t264;
t484 = Ifges(10,5) + Ifges(3,5);
t483 = Ifges(10,6) + Ifges(3,6);
t257 = sin(pkin(19));
t261 = cos(pkin(19));
t204 = t257 * t269 + t261 * t263;
t207 = -t263 * t257 + t261 * t269;
t137 = -t204 * t264 + t207 * t270;
t288 = t137 * qJD(1);
t113 = Ifges(9,4) * t288;
t512 = Ifges(9,2) * t288;
t511 = -t513 * t264 + t517;
t510 = t270 * t514 - t516;
t494 = m(5) + m(6);
t508 = m(8) + m(4);
t267 = sin(pkin(17));
t273 = cos(pkin(17));
t213 = t266 * t273 - t272 * t267;
t215 = t266 * t267 + t272 * t273;
t145 = t213 * t270 - t215 * t264;
t507 = -t145 / 0.2e1;
t446 = t145 / 0.2e1;
t435 = t270 / 0.2e1;
t134 = t204 * t270 + t207 * t264;
t463 = t134 * qJD(1);
t505 = t463 / 0.2e1;
t17 = t318 * t358 + t520;
t6 = t239 * t24 + (-t25 * t269 + (-t263 * t515 - t269 * t318) * qJD(3)) * pkin(1);
t504 = -t17 + t6;
t487 = qJD(2) / 0.2e1;
t503 = Ifges(7,4) * t145;
t502 = Ifges(9,4) * t463;
t398 = mrSges(11,1) + mrSges(4,1);
t397 = -mrSges(4,2) - mrSges(11,2);
t285 = qJD(2) * t302;
t16 = t269 * t285 + t520;
t20 = t239 * t515 + t269 * t302;
t366 = qJD(2) * qJD(3);
t202 = (-qJDD(2) * t269 + t263 * t366) * pkin(1);
t364 = qJDD(2) * t263;
t203 = (t269 * t366 + t364) * pkin(1);
t250 = qJDD(2) + qJDD(3);
t3 = -t24 * t358 - t25 * t304 + t515 * (pkin(5) * t250 + t203) - t318 * t202;
t500 = t434 * t16 + t20 * t3;
t265 = sin(qJ(1));
t271 = cos(qJ(1));
t338 = g(1) * t271 + g(2) * t265;
t180 = t266 * g(3) + t272 * t338;
t181 = g(3) * t272 - t266 * t338;
t254 = sin(pkin(22));
t258 = cos(pkin(22));
t102 = t180 * t254 + t181 * t258;
t255 = sin(pkin(21));
t259 = cos(pkin(21));
t311 = t180 * t258 - t181 * t254;
t498 = t256 * (t102 * t259 + t255 * t311) - t260 * (-t102 * t255 + t259 * t311);
t262 = sin(qJ(4));
t268 = cos(qJ(4));
t336 = mrSges(6,1) * t268 - mrSges(6,2) * t262;
t210 = pkin(9) * m(6) + mrSges(5,1) + t336;
t245 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t155 = m(11) * pkin(4) + t210 * t260 - t245 * t256;
t305 = t210 * t256 + t245 * t260;
t497 = t155 * t255 + t259 * t305 - mrSges(8,2);
t212 = -t263 * t264 + t269 * t270;
t197 = t212 * qJD(1);
t214 = t263 * t270 + t264 * t269;
t196 = t214 * qJD(1);
t416 = Ifges(4,4) * t196;
t109 = t197 * Ifges(4,2) + Ifges(4,6) * t252 + t416;
t182 = Ifges(4,4) * t197;
t110 = t196 * Ifges(4,1) + Ifges(4,5) * t252 + t182;
t346 = pkin(2) * t252;
t142 = t204 * t346;
t144 = t207 * t346;
t253 = qJ(3) + qJ(2);
t247 = cos(t253);
t246 = sin(t253);
t394 = Ifges(11,4) * t246;
t313 = Ifges(11,2) * t247 + t394;
t168 = Ifges(11,6) * t252 + qJD(1) * t313;
t377 = qJD(1) * t247;
t228 = Ifges(11,4) * t377;
t378 = qJD(1) * t246;
t169 = Ifges(11,1) * t378 + Ifges(11,5) * t252 + t228;
t376 = qJD(1) * t252;
t178 = qJDD(1) * t246 + t247 * t376;
t179 = qJDD(1) * t247 - t246 * t376;
t187 = t204 * qJD(3);
t188 = t207 * qJD(3);
t367 = qJD(1) * qJD(2);
t218 = qJDD(1) * t270 - t264 * t367;
t353 = t270 * t367;
t219 = -qJDD(1) * t264 - t353;
t419 = mrSges(10,3) * t264;
t221 = -qJD(2) * mrSges(10,2) - qJD(1) * t419;
t375 = qJD(1) * t270;
t222 = qJD(2) * mrSges(10,1) - mrSges(10,3) * t375;
t291 = t246 * (Ifges(11,1) * t247 - t394);
t312 = Ifges(11,5) * t247 - Ifges(11,6) * t246;
t395 = mrSges(11,3) * t247;
t337 = t359 * t395;
t342 = mrSges(4,3) * t358;
t343 = mrSges(4,3) * t359;
t432 = pkin(2) * t204;
t439 = -t252 / 0.2e1;
t444 = t196 / 0.2e1;
t453 = pkin(2) * m(10);
t454 = qJD(1) ^ 2;
t310 = -t187 * t270 - t188 * t264;
t56 = qJD(1) * t310 + t204 * t219 + t207 * t218;
t309 = -t187 * t264 + t188 * t270;
t57 = qJD(1) * t309 + t204 * t218 - t207 * t219;
t65 = Ifges(9,6) * t252 + t502 + t512;
t66 = Ifges(9,1) * t463 + Ifges(9,5) * t252 + t113;
t80 = (-t187 * t252 + t207 * t250) * pkin(2);
t81 = (t188 * t252 + t204 * t250) * pkin(2);
t287 = t214 * qJD(3);
t95 = -qJD(1) * t287 + t218 * t269 + t219 * t263;
t286 = t212 * qJD(3);
t96 = qJD(1) * t286 + t218 * t263 - t219 * t269;
t495 = -(-Ifges(11,2) * t378 + t169 + t228) * t377 / 0.2e1 - (-Ifges(4,2) * t196 + t110 + t182) * t197 / 0.2e1 - (-Ifges(9,2) * t463 + t113 + t66) * t288 / 0.2e1 + (Ifges(4,3) + Ifges(9,3) + Ifges(11,3)) * t250 + (t207 * (-qJDD(2) * mrSges(10,2) + mrSges(10,3) * t219) - t187 * t221 + t188 * t222) * pkin(2) - t312 * t376 / 0.2e1 + t65 * t505 + t202 * t397 - t196 * (Ifges(4,1) * t197 - t416) / 0.2e1 + (qJDD(2) * mrSges(10,1) - mrSges(10,3) * t218) * t432 + t109 * t444 + (Ifges(4,5) * t197 + Ifges(9,5) * t288 - Ifges(4,6) * t196 - Ifges(9,6) * t463) * t439 - t454 * t291 / 0.2e1 + (t142 * t188 - t144 * t187 + t204 * t81 + t207 * t80) * t453 + Ifges(11,5) * t178 + Ifges(11,6) * t179 + Ifges(4,6) * t95 + Ifges(4,5) * t96 + Ifges(9,6) * t56 + Ifges(9,5) * t57 - t196 * t342 + t197 * t343 + qJD(1) * t337 + t168 * t378 / 0.2e1 + t398 * t203 - (Ifges(9,1) * t288 - t502) * t463 / 0.2e1;
t493 = m(9) + m(3);
t492 = pkin(14) * m(7);
t491 = t179 / 0.2e1;
t440 = t250 / 0.2e1;
t490 = t252 / 0.2e1;
t344 = t213 * t264 + t215 * t270;
t482 = Ifges(7,4) * t344;
t489 = Ifges(7,2) * t446 + t482 / 0.2e1;
t488 = pkin(2) * t463;
t486 = m(11) + t508 + t494;
t111 = -t206 * t244 + t209 * t243;
t106 = t111 * qJD(1);
t307 = t243 * t206 + t209 * t244;
t107 = t307 * qJD(1);
t237 = pkin(1) * t264 - pkin(15);
t227 = t237 * qJD(1);
t159 = -pkin(5) * t197 + t227;
t49 = pkin(9) * t107 - pkin(11) * t106 + t159;
t11 = -t18 * t262 + t268 * t49;
t12 = t18 * t268 + t262 * t49;
t319 = -t11 * t262 + t12 * t268;
t478 = t319 * mrSges(6,3);
t64 = mrSges(5,1) * t107 + mrSges(5,2) * t106;
t476 = -m(5) * t159 - t64;
t473 = (mrSges(11,3) * t377 + mrSges(4,3) * t197 + t252 * t397) * t263;
t469 = qJD(1) * t511 + qJD(2) * t483;
t468 = qJD(1) * t510 + qJD(2) * t484;
t467 = -t264 * t484 - t270 * t483;
t101 = qJD(4) + t107;
t420 = mrSges(6,3) * t106;
t62 = -mrSges(6,2) * t101 - t262 * t420;
t63 = mrSges(6,1) * t101 - t268 * t420;
t465 = -t262 * t62 - t268 * t63;
t4 = -t515 * t202 + (t24 * t252 - t250 * t318) * pkin(5) + (-t318 * t364 + (t24 * t263 + (-qJD(3) * t318 - t25) * t269) * qJD(2)) * pkin(1);
t103 = t111 * qJDD(1);
t104 = t307 * qJDD(1);
t226 = t237 * qJDD(1);
t231 = pkin(1) * t353;
t189 = t231 + t226;
t83 = -pkin(5) * t95 + t189;
t44 = pkin(9) * t104 - pkin(11) * t103 + t83;
t1 = qJD(4) * t11 + t262 * t44 + t268 * t4;
t2 = -qJD(4) * t12 - t262 * t4 + t268 * t44;
t464 = t1 * t268 - t2 * t262;
t462 = 0.2e1 * t440;
t314 = -mrSges(11,1) * t247 + mrSges(11,2) * t246;
t205 = -t266 * t254 - t258 * t272;
t208 = t254 * t272 - t266 * t258;
t332 = -mrSges(8,1) * t205 + mrSges(8,2) * t208;
t461 = -mrSges(4,1) * t197 + mrSges(4,2) * t196 + (t314 + t332) * qJD(1);
t248 = mrSges(9,1) + t453;
t191 = mrSges(9,2) * t261 + t248 * t257 - t397;
t282 = pkin(5) * t494 - t257 * mrSges(9,2) + t248 * t261 + t398;
t459 = pkin(1) * t486 + t191 * t269 + t263 * t282 + mrSges(3,1) + mrSges(10,1);
t458 = -t191 * t263 + t269 * t282 - mrSges(3,2) - mrSges(10,2);
t457 = -pkin(15) * (mrSges(3,1) * t270 - mrSges(3,2) * t264) - (-t513 * t270 - t516) * t264 / 0.2e1 + (-t264 * t514 - t517) * t435;
t99 = (-t205 * t259 - t208 * t255) * pkin(4) + t237;
t400 = t99 * qJD(1) * (mrSges(11,1) * t246 + mrSges(11,2) * t247);
t456 = pkin(15) * (mrSges(9,1) * t463 + mrSges(9,2) * t288) - t237 * (mrSges(4,1) * t196 + mrSges(4,2) * t197) - t400;
t451 = Ifges(7,6) * t487 + qJD(1) * t489;
t9 = -t263 * t285 + t523;
t449 = t16 * t9;
t445 = t344 / 0.2e1;
t331 = mrSges(10,1) * t264 + mrSges(10,2) * t270;
t217 = t331 * qJD(1);
t431 = pkin(2) * t217;
t430 = pkin(5) * t196;
t429 = m(11) * t99;
t392 = qJD(1) * pkin(15);
t105 = -pkin(2) * t288 - t392;
t427 = m(10) * t105;
t425 = t111 * t3;
t415 = Ifges(6,4) * t262;
t414 = Ifges(6,4) * t268;
t329 = Ifges(6,1) * t268 - t415;
t408 = t262 * (t101 * Ifges(6,5) + t106 * t329);
t325 = -Ifges(6,2) * t262 + t414;
t406 = t268 * (t101 * Ifges(6,6) + t106 * t325);
t396 = mrSges(11,3) * t246;
t393 = Ifges(11,4) * t247;
t390 = t111 * t262;
t389 = t111 * t268;
t354 = mrSges(11,3) * t378;
t380 = -mrSges(4,3) * t196 + t252 * t398 - t354;
t374 = qJD(2) * t264;
t372 = qJD(2) * t270;
t371 = qJD(4) * t111;
t370 = qJD(4) * t262;
t369 = qJD(4) * t268;
t368 = qJDD(1) * pkin(15);
t100 = qJDD(4) + t104;
t60 = t103 * t268 - t106 * t370;
t61 = -t103 * t262 - t106 * t369;
t361 = Ifges(6,5) * t60 + Ifges(6,6) * t61 + Ifges(6,3) * t100;
t357 = pkin(1) * t372;
t341 = t373 * t396;
t171 = -pkin(5) * t212 + t237;
t335 = mrSges(6,1) * t262 + mrSges(6,2) * t268;
t334 = mrSges(7,1) * t272 + mrSges(7,2) * t266;
t333 = mrSges(7,1) * t266 - mrSges(7,2) * t272;
t328 = Ifges(7,1) * t344 + t503;
t320 = t11 * t268 + t12 * t262;
t317 = -t262 * t63 + t268 * t62;
t308 = t205 * t270 + t208 * t264;
t139 = t205 * t264 - t208 * t270;
t303 = -mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3) + mrSges(11,3) + t335;
t300 = t16 * t336;
t299 = t101 * (-Ifges(6,5) * t262 - Ifges(6,6) * t268);
t298 = t105 * (mrSges(10,1) * t270 - mrSges(10,2) * t264);
t297 = t262 * (-Ifges(6,2) * t268 - t415);
t294 = t268 * (-Ifges(6,1) * t262 - t414);
t45 = mrSges(6,1) * t100 - mrSges(6,3) * t60;
t46 = -mrSges(6,2) * t100 + mrSges(6,3) * t61;
t281 = qJD(4) * t465 - t262 * t45 + t268 * t46;
t223 = mrSges(7,1) * t267 + mrSges(7,2) * t273;
t224 = mrSges(7,1) * t273 - mrSges(7,2) * t267;
t85 = t155 * t259 - t255 * t305 + mrSges(8,1);
t279 = t492 + t459 * t264 + (t223 * t264 - t270 * t224 + t85 * t254 + t258 * t497) * t266 - t458 * t270 + (t270 * t223 + t224 * t264 - t254 * t497 + t85 * t258) * t272 - mrSges(2,1) + (-m(10) - t486 - t493) * pkin(15);
t225 = g(1) * t265 - g(2) * t271;
t162 = pkin(1) * t375 + t430;
t152 = -qJD(2) * t214 - t287;
t151 = qJD(2) * t212 + t286;
t143 = t252 * t432;
t130 = t145 * qJD(2);
t129 = t344 * qJD(2);
t126 = t139 * qJD(2);
t125 = t308 * qJD(2);
t90 = qJDD(1) * t99 + t231;
t86 = -g(3) * t282 + t191 * t338;
t84 = g(3) * t191 + t282 * t338;
t82 = -t267 * t334 + t273 * t333 + t458;
t79 = -qJD(1) * t129 + qJDD(1) * t145;
t78 = qJD(1) * t130 + qJDD(1) * t344;
t73 = Ifges(7,5) * qJD(2) + qJD(1) * t328;
t69 = (-qJD(2) * t125 - qJDD(2) * t139) * pkin(1);
t68 = (qJD(2) * t126 - qJDD(2) * t308) * pkin(1);
t67 = t267 * t333 + t273 * t334 + t459;
t58 = t335 * t106;
t51 = -pkin(2) * t56 - t368;
t27 = -mrSges(6,1) * t61 + mrSges(6,2) * t60;
t21 = -t239 * t318 + t362;
t14 = t162 * t262 + t17 * t268;
t13 = t162 * t268 - t17 * t262;
t10 = -t359 * t515 + t520;
t8 = t10 * t268 + t262 * t430;
t7 = -t10 * t262 + t268 * t430;
t5 = [(-mrSges(4,1) * t95 + mrSges(4,2) * t96 + t508 * (t189 + t231)) * t237 + t178 * t393 / 0.2e1 + t103 * t111 * Ifges(5,1) + (t142 * t374 - t144 * t372 - t270 * t81) * mrSges(10,3) + (-mrSges(9,2) * pkin(15) + Ifges(9,1) * t134 + Ifges(9,4) * t137) * t57 + (mrSges(9,1) * pkin(15) + Ifges(9,4) * t134 + Ifges(9,2) * t137) * t56 + (0.2e1 * Ifges(9,5) * t134 + Ifges(11,6) * t247 + 0.2e1 * Ifges(9,6) * t137) * t440 + t307 * t361 / 0.2e1 + (-t307 * t4 + t425) * mrSges(5,3) + t100 * (Ifges(6,3) * t307 + (Ifges(6,5) * t268 - Ifges(6,6) * t262) * t111) / 0.2e1 + (m(6) * (qJD(4) * t319 + t1 * t262 + t2 * t268) - t63 * t370 + t268 * t45 + t262 * t46 + t62 * t369) * (pkin(9) * t307 - pkin(11) * t111 + t171) + t2 * (mrSges(6,1) * t307 - mrSges(6,3) * t389) + t1 * (-mrSges(6,2) * t307 - mrSges(6,3) * t390) + t60 * (Ifges(6,5) * t307 + t111 * t329) / 0.2e1 + t61 * (Ifges(6,6) * t307 + t111 * t325) / 0.2e1 + t83 * (mrSges(5,1) * t307 + mrSges(5,2) * t111) + (-t103 * t307 - t104 * t111) * Ifges(5,4) + (Ifges(7,5) * t344 + Ifges(7,6) * t145 + t484 * t270) * qJDD(2) / 0.2e1 + (Ifges(8,1) * t208 ^ 2 + Ifges(2,3) + (0.2e1 * Ifges(8,4) * t208 + Ifges(8,2) * t205) * t205 + (-mrSges(7,1) * t145 + mrSges(7,2) * t344 + t492) * pkin(14) + t493 * pkin(15) ^ 2) * qJDD(1) + (pkin(1) * t341 + t312 * t490 - t337 + t400) * t252 + (Ifges(11,1) * t178 + Ifges(11,4) * t491 + t462 * Ifges(11,5) + t168 * t439) * t246 + t461 * t357 + (t231 + t90) * t429 - t468 * t374 / 0.2e1 - t469 * t372 / 0.2e1 + (-t189 * mrSges(4,1) + t202 * mrSges(4,3) + Ifges(4,4) * t96 + Ifges(4,2) * t95 + Ifges(4,6) * t462) * t212 + (t189 * mrSges(4,2) - t203 * mrSges(4,3) + Ifges(4,1) * t96 + Ifges(4,4) * t95 + Ifges(4,5) * t462) * t214 + (Ifges(4,5) * t151 + Ifges(4,6) * t152 + t247 * t169) * t490 - t80 * t419 + t335 * t425 + (m(10) * t51 - mrSges(10,1) * t219 + mrSges(10,2) * t218) * (-pkin(2) * t137 - pkin(15)) + (t106 * t294 + t299) * t371 / 0.2e1 - (t106 * t297 + t406 + t408) * t371 / 0.2e1 - t152 * t342 + (Ifges(6,1) * t60 + Ifges(6,4) * t61 + Ifges(6,5) * t100) * t389 / 0.2e1 - (Ifges(6,4) * t60 + Ifges(6,2) * t61 + Ifges(6,6) * t100) * t390 / 0.2e1 - t203 * t396 + t457 * t367 + (Ifges(4,1) * t151 + Ifges(4,4) * t152) * t444 + (Ifges(7,1) * t78 + Ifges(7,4) * t79 + Ifges(7,5) * qJDD(2)) * t445 + (Ifges(7,4) * t78 + Ifges(7,2) * t79 + Ifges(7,6) * qJDD(2)) * t446 + (t226 + t189) * t332 + (mrSges(9,1) * t137 - mrSges(3,2) * t270 - mrSges(9,2) * t134) * t368 + (-mrSges(4,1) * t152 + mrSges(4,2) * t151) * t227 + t300 * t371 + t202 * t395 + (t247 * (-Ifges(11,2) * t246 + t393) + t291) * t376 / 0.2e1 + (m(6) * t320 - t465 - t476) * (-pkin(5) * t152 + t357) - t371 * t478 - t151 * t343 - t129 * t451 + (t205 * t68 - t208 * t69) * mrSges(8,3) + (pkin(14) * (mrSges(7,1) * t129 + mrSges(7,2) * t130) + (Ifges(7,1) * t130 - Ifges(7,4) * t129) * t445 + (Ifges(7,4) * t130 - Ifges(7,2) * t129) * t446) * qJD(1) + (Ifges(7,5) * t130 - Ifges(7,6) * t129) * t487 + t51 * t331 + t78 * t328 / 0.2e1 + t247 * (Ifges(11,4) * t178 + Ifges(11,2) * t179 + Ifges(11,6) * t250) / 0.2e1 + t90 * t314 - pkin(15) * (-mrSges(3,1) * t219 + mrSges(3,2) * t218) - t265 * (g(1) * t279 + g(2) * t303) + (-g(1) * t303 + g(2) * t279) * t271 + t197 * (Ifges(4,4) * t151 + Ifges(4,2) * t152) / 0.2e1 + t99 * (-mrSges(11,1) * t179 + mrSges(11,2) * t178) + t152 * t109 / 0.2e1 + t151 * t110 / 0.2e1 + t130 * t73 / 0.2e1 + pkin(14) * (-mrSges(7,1) * t79 + mrSges(7,2) * t78) + (m(5) * t83 + mrSges(5,1) * t104 + mrSges(5,2) * t103) * t171 + (-mrSges(3,1) * t368 - t483 * qJDD(2) - t513 * t219 / 0.2e1 - t485 * t218 / 0.2e1) * t264 + (t484 * qJDD(2) + t218 * t514 + t485 * t219) * t435 + (Ifges(9,1) * t505 - mrSges(9,2) * t392 + t113 / 0.2e1 + t66 / 0.2e1 + Ifges(9,5) * t490) * (qJD(2) * t137 + t309) + (Ifges(9,4) * t505 + mrSges(9,1) * t392 - t431 + t512 / 0.2e1 - pkin(2) * t427 + t65 / 0.2e1 + Ifges(9,6) * t490) * (-qJD(2) * t134 + t310) + (t467 * t487 + t298) * qJD(2) + t510 * t218 / 0.2e1 + t511 * t219 / 0.2e1 + t79 * t489 + t313 * t491 + t104 * Ifges(5,2) * t307; t483 * t219 + t484 * t218 + (0.2e1 * (m(11) / 0.2e1 + m(4) / 0.2e1) * (-t202 * t269 + t203 * t263) + (-t341 - t461 * t270 + (t125 * t208 + t126 * t205 + (-t139 * t205 - t208 * t308) * qJD(2)) * mrSges(8,3)) * qJD(1) + (-mrSges(11,3) * t179 - mrSges(4,3) * t95 - t250 * t397) * t269 + (-mrSges(11,3) * t178 - mrSges(4,3) * t96 + t250 * t398) * t263 + (t139 * t208 - t205 * t308) * qJDD(1) * mrSges(8,3) + (-t508 * t237 - t429) * t454 * t270 + (t269 * t380 + t473) * qJD(3) + (-t308 * t68 - t139 * t69 + (t125 * t139 - t126 * t308) * t421) * m(8)) * pkin(1) + (-t298 + t73 * t507 + t344 * t451 + (-pkin(14) * (mrSges(7,1) * t344 + mrSges(7,2) * t145) + (-Ifges(7,2) * t344 + t503) * t507 - t344 * (Ifges(7,1) * t145 - t482) / 0.2e1 - t457) * qJD(1) + (-t142 * t264 + t144 * t270) * mrSges(10,3) + t468 * t264 / 0.2e1 + t469 * t435 - (Ifges(7,5) * t145 - Ifges(7,6) * t344 + t467) * qJD(2) / 0.2e1 + t456) * qJD(1) - (t217 + t427) * t488 + (t319 * t6 + (-qJD(4) * t320 + t464) * t21 - t11 * t13 - t12 * t14 + t500) * m(6) + t281 * t21 + t434 * t58 + (Ifges(7,3) + Ifges(10,3) + Ifges(3,3)) * qJDD(2) + t495 + (-t159 * t162 + t18 * t504 + t21 * t4 + t500) * m(5) + (g(3) * t67 + t338 * t82) * t264 + (-g(3) * t82 + t338 * t67) * t270 + t317 * t6 - t162 * t64 - t80 * mrSges(10,2) + t81 * mrSges(10,1) + Ifges(7,5) * t78 + Ifges(7,6) * t79 - t14 * t62 - t13 * t63 + t20 * t27 + (t20 * t103 - t21 * t104 + t434 * t106 - t107 * t504) * mrSges(5,3); t456 * qJD(1) - m(10) * (t105 * t488 + (t142 - t143) * t144) + (-t473 + (-t354 - t380) * t269) * t421 - t463 * t431 - m(5) * (t10 * t18 + t449) - m(6) * (t11 * t7 + t12 * t8 + t449) + (t263 * t84 + t269 * t86) * t270 + (-t263 * t86 + t269 * t84) * t264 + t143 * t221 - t144 * t222 - t8 * t62 - t7 * t63 - t9 * t58 + (t10 * t107 - t106 * t9) * mrSges(5,3) + (-t25 * t58 + t515 * t27 + (t103 * t515 - t106 * t25) * mrSges(5,3) + t494 * (-t16 * t25 + t3 * t515) + (m(5) * t18 + m(6) * t319 - t107 * mrSges(5,3) + t317) * t24 + t476 * t196 - (t281 + m(6) * (-t11 * t369 - t12 * t370 + t464) + m(5) * t4 - t104 * mrSges(5,3)) * t318) * pkin(5) + t495; -t11 * t62 + t12 * t63 + (t262 * t225 + t268 * t498 - t1) * mrSges(6,2) + (-t225 * t268 + t262 * t498 + t2) * mrSges(6,1) + (t406 / 0.2e1 + t408 / 0.2e1 - t299 / 0.2e1 - t300 + (t297 / 0.2e1 - t294 / 0.2e1) * t106 + t478) * t106 + t361;];
tau = t5;
