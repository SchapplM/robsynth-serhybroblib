% Calculate joint inertia matrix for
% palh1m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_inertiaJ_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_inertiaJ_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1IC_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1IC_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:03:03
% EndTime: 2020-04-15 20:03:05
% DurationCPUTime: 1.89s
% Computational Cost: add. (4186->293), mult. (5966->443), div. (156->8), fcn. (5676->32), ass. (0->172)
t345 = pkin(20) + qJ(7) + qJ(2);
t326 = sin(t345);
t328 = cos(t345);
t362 = sin(qJ(6));
t370 = cos(qJ(6));
t272 = 0.1e1 / (-t326 * t362 - t328 * t370) / pkin(7) / pkin(3);
t366 = sin(qJ(2));
t350 = t366 * pkin(1);
t307 = pkin(3) * t326 + t350;
t374 = cos(qJ(2));
t308 = pkin(1) * t374 + pkin(3) * t328;
t233 = (t307 * t362 + t308 * t370) * t272 * pkin(7);
t356 = sin(pkin(19));
t357 = cos(pkin(19));
t358 = cos(qJ(10));
t417 = sin(qJ(10));
t295 = t356 * t358 - t357 * t417;
t296 = -t356 * t417 - t357 * t358;
t361 = sin(qJ(7));
t426 = pkin(1) * t361;
t314 = pkin(4) * t356 + t426;
t369 = cos(qJ(7));
t424 = pkin(1) * t369;
t316 = pkin(4) * t357 + t424;
t245 = -t295 * t314 + t296 * t316;
t244 = t245 * mrSges(11,1);
t256 = (-t295 * t356 + t296 * t357) * pkin(4);
t254 = t256 * mrSges(11,1);
t405 = Ifges(11,3) + t254;
t246 = t295 * t316 + t296 * t314;
t408 = t246 * mrSges(11,2);
t257 = (t295 * t357 + t296 * t356) * pkin(4);
t434 = t257 * mrSges(11,2);
t438 = Ifges(11,3) + t244 - t408 + t233 * (t405 - t434);
t351 = -qJ(7) + pkin(19);
t343 = -qJ(10) + t351;
t323 = sin(t343);
t286 = t323 * pkin(8) - pkin(4) * sin(t351);
t324 = cos(t343);
t287 = t324 * pkin(8) - pkin(4) * cos(t351);
t344 = qJ(3) + qJ(4) + pkin(18);
t325 = sin(t344);
t327 = cos(t344);
t270 = 0.1e1 / (-t323 * t325 + t324 * t327) / pkin(10) / pkin(8);
t384 = t270 * t233;
t212 = (-t286 * t325 + t287 * t327) * pkin(10) * t384;
t207 = -Ifges(11,3) * t212 + t438;
t365 = sin(qJ(3));
t305 = pkin(5) * t365 + pkin(10) * t325;
t373 = cos(qJ(3));
t306 = pkin(5) * t373 + pkin(10) * t327;
t227 = (-t305 * t327 + t306 * t325) * t270 * pkin(10);
t425 = pkin(1) * t365;
t333 = pkin(5) + t425;
t364 = sin(qJ(4));
t372 = cos(qJ(4));
t423 = pkin(1) * t373;
t281 = t333 * t372 + t364 * t423;
t279 = -pkin(9) - t281;
t282 = t364 * t333 - t372 * t423;
t280 = pkin(11) + t282;
t421 = pkin(5) * t372;
t334 = -pkin(9) - t421;
t336 = mrSges(4,1) * t425;
t340 = mrSges(4,2) * t423;
t371 = cos(qJ(5));
t363 = sin(qJ(5));
t412 = Ifges(6,4) * t363;
t317 = Ifges(6,2) * t371 + t412;
t411 = Ifges(6,4) * t371;
t318 = Ifges(6,1) * t363 + t411;
t396 = t371 * t317 + t363 * t318 + Ifges(5,3);
t394 = Ifges(4,3) + t396;
t422 = pkin(5) * t364;
t332 = pkin(11) + t422;
t353 = t363 ^ 2;
t354 = t371 ^ 2;
t399 = t353 + t354;
t395 = t399 * t332;
t355 = qJ(8) + qJ(9);
t346 = sin(t355);
t360 = sin(qJ(8));
t310 = pkin(2) * t360 - pkin(12) * t346;
t347 = cos(t355);
t368 = cos(qJ(8));
t311 = -pkin(2) * t368 + pkin(12) * t347;
t352 = qJ(3) + pkin(17);
t341 = sin(t352);
t342 = cos(t352);
t420 = pkin(6) / (t310 * t347 + t311 * t346) / pkin(12);
t224 = (-t310 * t341 + t311 * t342) * t420;
t228 = (-t341 * t346 - t342 * t347) * pkin(12) * t420;
t359 = sin(qJ(9));
t415 = mrSges(10,2) * t359;
t367 = cos(qJ(9));
t416 = mrSges(10,1) * t367;
t309 = Ifges(10,3) + (t415 - t416) * pkin(2);
t392 = Ifges(10,3) + Ifges(9,3) + (0.2e1 * t415 - 0.2e1 * t416 + m(10) * (t359 ^ 2 + t367 ^ 2) * pkin(2)) * pkin(2);
t404 = t224 * t309 + t228 * t392;
t391 = mrSges(6,1) * t371 - mrSges(6,2) * t363;
t285 = t334 * t391;
t414 = mrSges(6,3) * t353;
t312 = t332 * t414;
t413 = mrSges(6,3) * t354;
t313 = t332 * t413;
t339 = mrSges(5,1) * t421;
t431 = -t285 + t312 + t313 + t339;
t271 = t279 * t391;
t273 = t280 * t414;
t274 = t280 * t413;
t277 = t281 * mrSges(5,1);
t432 = -t271 + t273 + t274 + t277;
t437 = m(6) * (t279 * t334 + t280 * t395) + m(5) * (t281 * t372 + t282 * t364) * pkin(5) + t336 + t340 + (-t282 - t422) * mrSges(5,2) + t394 + t404 + t431 + t432 + t207 * t227;
t433 = 0.2e1 * mrSges(8,1) * t424 - 0.2e1 * mrSges(8,2) * t426 + 0.2e1 * t244;
t304 = -t365 * t366 + t373 * t374;
t329 = -pkin(15) + t350;
t275 = -pkin(5) * t304 + t329;
t429 = 0.2e1 * t275;
t418 = pkin(9) * t391;
t410 = t282 * mrSges(5,2);
t409 = Ifges(11,3) + Ifges(8,3);
t301 = t365 * t374 + t366 * t373;
t269 = t301 * t372 + t304 * t364;
t407 = t269 * t363;
t406 = t269 * t371;
t300 = -t361 * t374 - t366 * t369;
t303 = -t361 * t366 + t369 * t374;
t242 = -t295 * t303 + t296 * t300;
t243 = t295 * t300 + t296 * t303;
t219 = Ifges(11,5) * t243 + Ifges(11,6) * t242;
t267 = t301 * t364 - t372 * t304;
t402 = Ifges(6,5) * t406 + Ifges(6,3) * t267;
t299 = -t360 * t374 - t366 * t368;
t302 = -t360 * t366 + t368 * t374;
t266 = -t299 * t367 + t302 * t359;
t268 = -t299 * t359 - t302 * t367;
t401 = Ifges(10,5) * t268 + Ifges(10,6) * t266;
t400 = Ifges(6,5) * t363 + Ifges(6,6) * t371;
t397 = mrSges(5,2) * t422;
t393 = Ifges(8,5) * t303 + Ifges(8,6) * t300 + t219;
t390 = mrSges(6,1) * t363 + mrSges(6,2) * t371;
t236 = -mrSges(6,2) * t267 - mrSges(6,3) * t407;
t237 = mrSges(6,1) * t267 - mrSges(6,3) * t406;
t388 = t236 * t371 - t237 * t363;
t222 = Ifges(6,6) * t267 + (-Ifges(6,2) * t363 + t411) * t269;
t223 = Ifges(6,5) * t267 + (Ifges(6,1) * t371 - t412) * t269;
t386 = t363 * t223 / 0.2e1 + t371 * t222 / 0.2e1 - t317 * t407 / 0.2e1 + t318 * t406 / 0.2e1 + Ifges(5,5) * t269 + (t400 / 0.2e1 - Ifges(5,6)) * t267;
t335 = pkin(11) * t414;
t338 = pkin(11) * t413;
t385 = t335 + t338 + t396 + t418;
t383 = Ifges(4,5) * t301 + Ifges(4,6) * t304 + t386;
t381 = Ifges(9,5) * t302 + Ifges(9,6) * t299 + t401 + (-t266 * t359 + t268 * t367) * pkin(2) * mrSges(10,3);
t213 = m(6) * (pkin(11) * t280 * t399 - pkin(9) * t279) + t385 - t410 + t432;
t216 = m(6) * (-pkin(9) * t334 + pkin(11) * t395) + t385 - t397 + t431;
t283 = -pkin(2) * t299 - pkin(15);
t276 = -pkin(11) * t390 + t400;
t253 = (-t300 * t357 - t303 * t356) * pkin(4) + t329;
t239 = 0.2e1 * t335 + 0.2e1 * t338 + 0.2e1 * t418 + m(6) * (pkin(11) ^ 2 * t399 + pkin(9) ^ 2) + t396;
t235 = t390 * t269;
t231 = (t307 * t328 - t308 * t326) * t272 * pkin(3);
t230 = pkin(9) * t267 - pkin(11) * t269 + t275;
t225 = (t305 * t323 - t306 * t324) * t270 * pkin(8);
t215 = t225 * t276 - t332 * t390 + t400;
t211 = (t286 * t324 - t287 * t323) * pkin(8) * t384;
t210 = -Ifges(6,6) * t407 + t230 * t391 + t402;
t209 = t225 * t239 + t216;
t208 = -t211 * t276 - t280 * t390 + t400;
t206 = -t211 * t239 + t213;
t205 = -pkin(9) * t235 + pkin(11) * t388 + t386;
t204 = t228 * t381 + t383 + (-t267 * t364 - t269 * t372) * mrSges(5,3) * pkin(5) + t388 * t332 + t225 * t205 + t227 * t219 + t224 * t401 + t334 * t235;
t203 = t381 + t383 + t388 * t280 + ((t300 * t361 - t303 * t369) * mrSges(8,3) + (-t301 * t365 - t304 * t373) * mrSges(4,3)) * pkin(1) + (-t267 * t282 - t269 * t281) * mrSges(5,3) + t233 * t393 + (-t245 * t243 + t246 * t242 + t233 * (t242 * t257 - t243 * t256)) * mrSges(11,3) - t211 * t205 - t212 * t219 + t279 * t235 - Ifges(3,6) * t366 + t231 * (Ifges(7,5) * t362 + Ifges(7,6) * t370) + Ifges(3,5) * t374 + t393;
t1 = [(mrSges(5,2) * t429 + Ifges(5,1) * t269 - t222 * t363 + t223 * t371 + (-Ifges(6,6) * t363 - (2 * Ifges(5,4))) * t267) * t269 + (mrSges(5,1) * t429 + Ifges(5,2) * t267 + t402) * t267 + (m(6) * t399 * t230 + 0.2e1 * t363 * t236 + 0.2e1 * t371 * t237) * t230 + 0.2e1 * (-mrSges(4,1) * t304 - mrSges(8,1) * t300 + mrSges(4,2) * t301 + mrSges(8,2) * t303) * t329 - 0.2e1 * (mrSges(3,1) * t366 - mrSges(9,1) * t299 + mrSges(3,2) * t374 + mrSges(9,2) * t302) * pkin(15) + (m(9) + m(3)) * pkin(15) ^ 2 + (m(8) + m(4)) * t329 ^ 2 + (0.2e1 * Ifges(8,4) * t303 + Ifges(8,2) * t300) * t300 + (-0.2e1 * Ifges(3,4) * t374 + Ifges(3,2) * t366) * t366 + Ifges(9,1) * t302 ^ 2 + Ifges(8,1) * t303 ^ 2 + Ifges(4,2) * t304 ^ 2 + Ifges(3,1) * t374 ^ 2 + (Ifges(4,1) * t301 + 0.2e1 * Ifges(4,4) * t304) * t301 + (0.2e1 * Ifges(9,4) * t302 + Ifges(9,2) * t299) * t299 + Ifges(2,3) + t242 * (Ifges(11,4) * t243 + Ifges(11,2) * t242) + t243 * (Ifges(11,1) * t243 + Ifges(11,4) * t242) + m(11) * t253 ^ 2 + 0.2e1 * t253 * (-mrSges(11,1) * t242 + mrSges(11,2) * t243) + t268 * (Ifges(10,1) * t268 + Ifges(10,4) * t266) + t266 * (Ifges(10,4) * t268 + Ifges(10,2) * t266) + m(5) * t275 ^ 2 + m(10) * t283 ^ 2 + 0.2e1 * t283 * (-mrSges(10,1) * t266 + mrSges(10,2) * t268) + m(7) * pkin(14) ^ 2 + t370 * (Ifges(7,4) * t362 + Ifges(7,2) * t370) + t362 * (Ifges(7,1) * t362 + Ifges(7,4) * t370) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t370 + mrSges(7,2) * t362), t203, t204, t210; t203, -0.2e1 * t410 + t409 - 0.2e1 * t408 + t394 + t392 + m(6) * (t280 ^ 2 * t399 + t279 ^ 2) + (m(8) * (t361 ^ 2 + t369 ^ 2) + m(4) * (t365 ^ 2 + t373 ^ 2)) * pkin(1) ^ 2 + (0.2e1 * (-t246 - t257) * mrSges(11,2) + 0.2e1 * Ifges(8,3) + 0.2e1 * m(11) * (t245 * t256 + t246 * t257) + 0.2e1 * t405 + (m(11) * (t256 ^ 2 + t257 ^ 2) + 0.2e1 * t254 + t409 - 0.2e1 * t434) * t233 + t433) * t233 + t433 - 0.2e1 * t271 + 0.2e1 * t336 + 0.2e1 * t340 + m(5) * (t281 ^ 2 + t282 ^ 2) + m(11) * (t245 ^ 2 + t246 ^ 2) - (t207 + t438) * t212 - (t206 + t213) * t211 + 0.2e1 * t273 + 0.2e1 * t274 + t231 ^ 2 * Ifges(7,3) + Ifges(3,3) + 0.2e1 * t277, t206 * t225 - t211 * t216 + t437, t208; t204, -t209 * t211 + t225 * t213 + t437, 0.2e1 * t312 + 0.2e1 * t313 - 0.2e1 * t397 + 0.2e1 * t339 - 0.2e1 * t285 + t404 * t228 + (Ifges(10,3) * t224 + t228 * t309) * t224 + t227 ^ 2 * Ifges(11,3) + (t216 + t209) * t225 + m(6) * (t332 ^ 2 * t399 + t334 ^ 2) + m(5) * (t364 ^ 2 + t372 ^ 2) * pkin(5) ^ 2 + t394, t215; t210, t208, t215, Ifges(6,3);];
Mq = t1;
