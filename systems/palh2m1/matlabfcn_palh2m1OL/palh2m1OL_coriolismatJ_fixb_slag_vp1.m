% Calculate matrix of centrifugal and coriolis load on the joints for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh2m1OL_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1OL_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:21:55
% EndTime: 2020-05-03 00:22:07
% DurationCPUTime: 7.60s
% Computational Cost: add. (3248->509), mult. (4579->721), div. (0->0), fcn. (690->68), ass. (0->392)
t326 = 2 * qJ(4);
t328 = 2 * qJ(2);
t471 = 2 * qJ(3) + t328;
t243 = t326 + qJ(5) + t471;
t206 = cos(t243);
t315 = rSges(6,3) + pkin(6);
t278 = t315 * rSges(6,1);
t318 = rSges(6,2) * pkin(4);
t150 = (t278 - t318) * m(6) - Icges(6,5);
t593 = t150 / 0.2e1;
t424 = t206 * t593;
t604 = t150 * qJD(1) * t206;
t426 = -qJ(5) + t471;
t247 = qJ(4) + t426;
t200 = sin(t247);
t297 = qJ(4) - qJ(5);
t259 = sin(t297);
t603 = t200 + t259;
t210 = cos(t247);
t266 = cos(t297);
t602 = t210 + t266;
t527 = m(6) * t315;
t540 = rSges(5,2) * m(5);
t190 = -t527 + t540;
t299 = t328 + qJ(3);
t273 = qJ(4) + t299;
t226 = sin(t273);
t233 = cos(t273);
t238 = m(5) * rSges(5,1) + pkin(4) * m(6);
t298 = qJ(3) + qJ(4);
t260 = sin(t298);
t267 = cos(t298);
t601 = (t190 * (t233 + t267) + t238 * (t226 + t260)) * pkin(2) * qJD(1);
t296 = qJ(4) + qJ(5);
t271 = qJ(3) + t296;
t250 = t328 + t271;
t203 = sin(t250);
t224 = sin(t271);
t600 = t203 + t224;
t272 = qJ(3) + t297;
t251 = t328 + t272;
t204 = sin(t251);
t225 = sin(t272);
t599 = t204 + t225;
t213 = cos(t250);
t231 = cos(t271);
t598 = t213 + t231;
t214 = cos(t251);
t232 = cos(t272);
t597 = t214 + t232;
t248 = qJ(2) + t271;
t184 = 2 * t248;
t171 = sin(t184);
t249 = qJ(2) + t272;
t185 = 2 * t249;
t172 = sin(t185);
t329 = rSges(6,2) ^ 2;
t330 = rSges(6,1) ^ 2;
t177 = m(6) * (-t329 + t330) - Icges(6,1) + Icges(6,2);
t277 = qJ(2) + t298;
t218 = 2 * t277;
t179 = sin(t218);
t316 = Icges(6,2) / 0.2e1;
t528 = t330 / 0.2e1;
t346 = -(rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - (t528 + t329 / 0.2e1 - (pkin(4) + t315) * (-pkin(4) + t315)) * m(6) + Icges(5,1) - Icges(5,2) - Icges(6,3);
t308 = sin(qJ(4));
t485 = t238 * t308;
t312 = cos(qJ(4));
t488 = t190 * t312;
t48 = t485 + t488;
t265 = cos(t296);
t513 = rSges(6,2) * t265;
t258 = sin(t296);
t517 = rSges(6,1) * t258;
t279 = rSges(6,2) * t315;
t542 = rSges(6,1) * pkin(4);
t153 = m(6) * (-t279 + t542) + Icges(6,6);
t244 = t326 + t426;
t197 = sin(t244);
t565 = -t153 * t197 / 0.2e1;
t274 = qJ(4) + t471;
t234 = cos(t274);
t462 = qJD(1) * t190;
t72 = pkin(3) * t234 * t462;
t227 = sin(t274);
t460 = qJD(1) * t238;
t93 = pkin(3) * t227 * t460;
t594 = (pkin(3) * (t48 / 0.2e1 + m(6) * (t513 / 0.4e1 + t517 / 0.4e1)) + (t171 / 0.8e1 + t172 / 0.8e1) * t177 - t424 - t565 - (Icges(6,1) / 0.2e1 + t316 + t346) * t179 / 0.2e1) * qJD(1) + t72 / 0.2e1 + t93 / 0.2e1;
t230 = sin(t277);
t592 = t230 / 0.2e1;
t301 = qJ(2) + qJ(3);
t239 = 2 * t301;
t549 = m(5) + m(6);
t496 = (pkin(3) ^ 2 * t549 + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(4,1) + Icges(4,2)) * sin(t239);
t427 = -t496 / 0.2e1;
t246 = qJ(5) + t274;
t199 = sin(t246);
t209 = cos(t246);
t519 = m(6) * qJD(1);
t434 = t519 / 0.4e1;
t381 = rSges(6,2) * t434;
t356 = pkin(3) * t381;
t386 = rSges(6,1) * t434;
t364 = pkin(3) * t386;
t591 = t199 * t364 + t209 * t356;
t325 = 2 * qJ(5);
t245 = t325 + t277;
t198 = sin(t245);
t252 = -2 * qJ(5) + t277;
t205 = sin(t252);
t590 = t198 + t205;
t589 = t200 + t199;
t295 = -qJ(5) + qJ(2);
t275 = qJ(3) + t295;
t228 = sin(t275);
t300 = qJ(2) + qJ(5);
t276 = qJ(3) + t300;
t229 = sin(t276);
t588 = t228 + t229;
t208 = cos(t245);
t546 = m(6) * rSges(6,2);
t240 = rSges(6,1) * t546 - Icges(6,4);
t420 = t240 * t208 / 0.2e1;
t215 = cos(t252);
t484 = t240 * t215;
t587 = t420 - t484 / 0.2e1;
t586 = t513 + t517;
t268 = cos(t299);
t548 = m(4) * rSges(4,2);
t453 = t268 * t548;
t319 = rSges(4,1) * m(4);
t323 = pkin(3) * m(6);
t183 = m(5) * pkin(3) + t319 + t323;
t261 = sin(t299);
t490 = t183 * t261;
t309 = sin(qJ(3));
t313 = cos(qJ(3));
t65 = t183 * t309 + t313 * t548;
t585 = qJD(1) * (t427 - (t490 / 0.2e1 + t65 / 0.2e1 + t453 / 0.2e1) * pkin(2));
t535 = t177 / 0.4e1;
t583 = t590 * t535;
t435 = -t519 / 0.4e1;
t383 = rSges(6,2) * t435;
t358 = pkin(3) * t383;
t582 = t602 * t358 + t603 * t364;
t581 = 0.2e1 * m(5);
t303 = qJD(5) / 0.2e1;
t580 = -t303 / 0.2e1;
t510 = pkin(1) * qJD(1);
t289 = qJD(3) + qJD(4);
t256 = qJD(2) + t289;
t524 = pkin(4) * t256;
t167 = t510 + t524;
t174 = cos(t185);
t495 = Icges(6,5) * t256;
t193 = -0.8e1 * t495;
t459 = qJD(1) * t240;
t398 = t459 / 0.4e1;
t212 = cos(t249);
t412 = t212 / 0.16e2;
t476 = t256 * t315;
t441 = rSges(6,1) * t476;
t552 = 0.8e1 * m(6);
t579 = ((rSges(6,2) * t167 + t441) * t552 + t193) * t412 + t174 * t398;
t257 = sin(t295);
t264 = cos(t295);
t447 = -t546 / 0.2e1;
t380 = qJD(5) * t447;
t547 = m(6) * rSges(6,1);
t448 = t547 / 0.2e1;
t385 = qJD(5) * t448;
t578 = (t257 * t385 + t264 * t380) * pkin(2);
t520 = m(4) + t549;
t541 = rSges(3,2) * m(3);
t353 = -((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + t520 * pkin(2) ^ 2) * sin(t328) / 0.2e1 - (rSges(3,1) * t541 - Icges(3,4)) * cos(t328);
t168 = t510 - t524;
t576 = rSges(6,2) * t168 + t441;
t470 = t329 + t330;
t506 = pkin(4) * qJD(1);
t575 = -Icges(6,3) * qJD(5) - (0.2e1 * pkin(1) * t506 + t470 * qJD(5)) * m(6);
t269 = cos(t300);
t574 = (-t269 / 0.2e1 + t264 / 0.2e1) * qJD(1);
t262 = sin(t300);
t573 = (-t262 / 0.2e1 - t257 / 0.2e1) * qJD(1);
t211 = cos(t248);
t293 = cos(t325);
t483 = t240 * t293;
t572 = qJD(1) * t211 * t593 + t256 * t483;
t504 = qJD(5) * pkin(4);
t219 = t504 + t510;
t501 = t219 * rSges(6,1);
t166 = 0.8e1 * t501;
t455 = qJD(5) * t315;
t437 = rSges(6,2) * t455;
t220 = -0.8e1 * t437;
t458 = qJD(5) * Icges(6,5);
t284 = 0.8e1 * t458;
t469 = Icges(6,6) * qJD(5);
t285 = 0.8e1 * t469;
t400 = -t459 / 0.4e1;
t202 = sin(t249);
t414 = t202 / 0.16e2;
t439 = rSges(6,1) * t455;
t514 = rSges(6,2) * t219;
t571 = ((-t439 - t514) * t552 + t284) * t412 + ((t166 + t220) * m(6) + t285) * t414 + t174 * t400;
t235 = cos(t275);
t290 = qJD(2) + qJD(3);
t526 = pkin(3) * t290;
t551 = m(6) / 0.2e1;
t391 = t526 * t551;
t373 = rSges(6,2) * t391;
t452 = -t323 / 0.2e1;
t392 = t290 * t452;
t374 = rSges(6,1) * t392;
t464 = qJD(1) * t177;
t404 = -t464 / 0.8e1;
t570 = t172 * t404 + t228 * t374 + t235 * t373;
t482 = (rSges(4,2) * t319 - Icges(4,4)) * cos(t239);
t569 = qJD(1) * t482 + t72 + t93;
t445 = rSges(6,2) * t476;
t169 = 0.8e1 * t445;
t494 = Icges(6,6) * t256;
t191 = -0.8e1 * t494;
t508 = pkin(2) * qJD(2);
t390 = t508 * t551;
t360 = rSges(6,2) * t390;
t449 = -t547 / 0.2e1;
t568 = t264 * t360 + t257 * t449 * t508 + ((-0.8e1 * rSges(6,1) * t167 + t169) * m(6) + t191) * t414;
t355 = pkin(3) * t380;
t363 = pkin(3) * t385;
t402 = t464 / 0.8e1;
t566 = t172 * t402 + t228 * t363 + t235 * t355;
t236 = cos(t276);
t201 = sin(t248);
t399 = t459 / 0.2e1;
t401 = -t459 / 0.2e1;
t465 = qJD(1) * t153;
t407 = -t465 / 0.2e1;
t152 = (t279 + t542) * m(6) - Icges(6,6);
t466 = qJD(1) * t152;
t409 = -t466 / 0.2e1;
t149 = (t278 + t318) * m(6) - Icges(6,5);
t467 = qJD(1) * t149;
t410 = t467 / 0.2e1;
t489 = (t470 * m(6) + Icges(6,3)) * t230;
t421 = -t489 / 0.2e1;
t307 = sin(qJ(5));
t311 = cos(qJ(5));
t175 = rSges(6,1) * t307 + rSges(6,2) * t311;
t525 = pkin(4) * t175;
t454 = m(6) * t525;
t291 = sin(t325);
t491 = t177 * t291;
t338 = qJD(1) * t421 + t201 * t409 + t202 * t407 + t208 * t401 + t212 * t410 + t215 * t399 + t572 - t590 * t464 / 0.4e1 + (t454 + t491 / 0.2e1) * t256;
t446 = t546 / 0.2e1;
t382 = qJD(1) * t446;
t357 = pkin(3) * t382;
t436 = -t519 / 0.2e1;
t384 = rSges(6,2) * t436;
t359 = pkin(3) * t384;
t388 = rSges(6,1) * t436;
t367 = pkin(3) * t388;
t562 = t235 * t357 + t236 * t359 + t259 * t374 + t266 * t373 + t588 * t367 + t586 * t391 + t338;
t361 = pkin(2) * t381;
t368 = pkin(2) * t386;
t559 = t598 * t361 + t600 * t368;
t362 = pkin(2) * t383;
t387 = rSges(6,1) * t435;
t369 = pkin(2) * t387;
t558 = t598 * t362 + t600 * t369;
t557 = t597 * t361 + t599 * t369;
t556 = t597 * t362 + t599 * t368;
t192 = 0.8e1 * t495;
t286 = -0.8e1 * t510;
t413 = t211 / 0.16e2;
t415 = t201 / 0.16e2;
t555 = (-t552 * t576 + t192) * t413 + (((t286 + 0.8e1 * t524) * rSges(6,1) + t169) * m(6) + t191) * t415 + (((t286 - 0.8e1 * t524) * rSges(6,1) + t169) * m(6) + t191) * t414 - t601 / 0.2e1;
t196 = sin(t243);
t207 = cos(t244);
t287 = qJD(4) + qJD(5);
t254 = qJD(3) + t287;
t221 = qJD(2) + t254;
t237 = cos(t277);
t492 = (rSges(5,1) * t540 - pkin(4) * t527 - Icges(5,4)) * cos(t218);
t45 = qJD(1) * t492;
t450 = rSges(5,1) * t510;
t511 = rSges(5,3) * t256;
t288 = qJD(4) - qJD(5);
t255 = qJD(3) + t288;
t222 = qJD(2) + t255;
t530 = t222 / 0.2e1;
t554 = ((-rSges(5,2) * t511 - t450) * t581 + 0.2e1 * t256 * Icges(5,6) + t575) * t592 + t207 * t410 + t196 * t409 + t221 * t420 + t484 * t530 + (t510 * t527 + (rSges(5,1) * t511 - rSges(5,2) * t510) * m(5) - Icges(5,5) * t256) * t237 - t45;
t553 = t556 + t559 + t601 / 0.2e1;
t550 = m(6) * pkin(2);
t545 = m(6) * t576;
t544 = m(6) * (t439 - t514);
t543 = rSges(3,1) * m(3);
t539 = rSges(3,3) * m(3);
t538 = rSges(4,3) * m(4);
t537 = rSges(5,3) * m(5);
t347 = t175 * qJD(5) * m(6);
t23 = qJD(4) * t190 + t347;
t536 = pkin(3) * t23;
t283 = qJD(2) + qJD(3) / 0.2e1;
t305 = qJD(4) / 0.2e1;
t223 = t305 + t283;
t534 = t580 - t223 / 0.2e1;
t396 = t305 + t290;
t533 = t580 - t396 / 0.2e1;
t532 = -t221 / 0.2e1;
t531 = -t222 / 0.2e1;
t529 = -t289 / 0.2e1;
t522 = qJD(1) / 0.2e1;
t521 = -qJD(3) / 0.2e1;
t518 = rSges(6,1) * t168;
t516 = rSges(6,1) * t262;
t512 = rSges(6,2) * t269;
t507 = pkin(2) * qJD(5);
t505 = qJD(4) * pkin(3);
t317 = -Icges(6,1) / 0.2e1;
t503 = (t317 - Icges(6,2) / 0.2e1 - t346) * t179;
t19 = t190 * t256 + t347;
t502 = t19 * t313;
t70 = (rSges(6,1) * t311 - rSges(6,2) * t307) * m(6) + t238;
t500 = t256 * t70;
t499 = t308 * t70;
t498 = t313 * t70;
t497 = t70 * t312;
t487 = t221 * t198;
t486 = t222 * t205;
t480 = t254 * t224;
t479 = t254 * t231;
t478 = t255 * t225;
t477 = t255 * t232;
t475 = t287 * t258;
t474 = t287 * t265;
t473 = t288 * t259;
t472 = t288 * t266;
t463 = qJD(1) * t183;
t457 = qJD(5) * t291;
t456 = qJD(5) * t293;
t451 = qJD(1) * t548;
t444 = t309 * t508;
t443 = t313 * t508;
t442 = t175 * t504;
t440 = t256 * t499;
t438 = qJD(4) * t499;
t431 = t507 / 0.2e1;
t430 = pkin(3) * qJD(1) / 0.4e1;
t429 = -t503 / 0.2e1;
t155 = t240 * t456;
t422 = -t491 / 0.2e1;
t411 = -t467 / 0.2e1;
t408 = t466 / 0.2e1;
t406 = t465 / 0.2e1;
t403 = t464 / 0.4e1;
t395 = pkin(1) * t451;
t394 = pkin(3) * t448;
t372 = pkin(3) * t229 * t449;
t371 = pkin(3) * t236 * t446;
t366 = pkin(3) * t387;
t365 = qJD(1) * t394;
t310 = sin(qJ(2));
t314 = cos(qJ(2));
t352 = (t520 * pkin(2) + t543) * t310 + t314 * t541;
t86 = -t177 * t457 / 0.2e1;
t350 = -m(6) * t442 - t155 + t86;
t173 = cos(t184);
t60 = t173 * t400;
t349 = t199 * t366 + t209 * t358 + t558 + t60;
t348 = (-t454 - t483) * qJD(5);
t85 = qJD(5) * t422;
t344 = t85 + t348;
t343 = t86 + t348;
t341 = t602 * t356 + t603 * t366 + t557 + t579;
t59 = t173 * t398;
t340 = (t284 - 0.8e1 * t544) * t413 + ((t166 + 0.8e1 * t437) * m(6) - 0.8e1 * t469) * t415 + (t450 * t581 - t575) * t592 + t207 * t411 + t196 * t408 + t45 + t59 + pkin(1) * t237 * t462 + t571 + t587 * qJD(5);
t339 = -t307 * ((rSges(6,1) * t524 - t437) * m(6) + t469) + ((-rSges(6,2) * t524 - t439) * t552 + t284) * t311 / 0.8e1 + t201 * t408 + t212 * t411 + t202 * t406 + t256 * t422 + t489 * t522 + t208 * t399 + t215 * t401 - t572 + t590 * t403;
t336 = t235 * t359 + t236 * t357 + t339 + t586 * t221 * t452 + t588 * t365 + (pkin(3) * t266 * t447 + t259 * t394) * t222;
t25 = t171 * t404;
t263 = sin(t301);
t270 = cos(t301);
t335 = (-t290 * (rSges(4,2) * t538 - Icges(4,6)) + (-t549 * pkin(3) - t319) * t510) * t263 + (-t395 + (rSges(4,1) * t538 + pkin(3) * t537 - Icges(4,5)) * t290) * t270 - t177 * t486 / 0.4e1 + t554 + rSges(6,1) * t229 * t391 + t290 * t371 + t210 * t357 + t209 * t359 - t569 + t570 + t60 + t25 + t197 * t407 + t487 * t535 + t579 + t589 * t367 + (t424 + t429) * qJD(1);
t24 = t171 * t402;
t334 = -t604 / 0.2e1 + t270 * t395 + t209 * t357 + t210 * t359 + t229 * t363 + t566 + t24 + t340 + pkin(1) * t263 * t463 + t503 * t522 + t197 * t406 + t569 + t589 * t365 + (t371 + t583) * qJD(5);
t333 = t604 / 0.4e1 + t197 * t465 / 0.4e1 + t291 * t403 + t293 * t399 - t207 * t467 / 0.4e1 - t196 * t466 / 0.4e1 + (t583 + t587) * t256;
t304 = -qJD(5) / 0.2e1;
t187 = t304 + t396;
t182 = t304 + t223;
t126 = t190 * t443;
t49 = t256 * t421;
t1 = [(t427 - t482) * t290 + (t429 - t492) * t256 + (t149 * t207 / 0.2e1 + t565) * (t304 + t256) + (t424 - t152 * t196 / 0.2e1) * (t303 + t256) + t353 * qJD(2) + (-t221 * t173 / 0.4e1 + t222 * t174 / 0.4e1 + t456 / 0.2e1) * t240 + (-t221 * t171 / 0.8e1 - t222 * t172 / 0.8e1 + t457 / 0.4e1) * t177 + ((-t190 * t234 - t238 * t227) * t396 + (-t485 / 0.2e1 - t488 / 0.2e1) * qJD(4)) * pkin(3) + ((-t183 * t263 - t270 * t548) * t290 + (-t190 * t237 - t230 * t238) * t256 - t352 * qJD(2)) * pkin(1) + ((-t223 * t226 + t260 * t529) * t238 + (-t223 * t233 + t267 * t529) * t190 + (-t283 * t261 + t309 * t521) * t183 + (-t283 * t268 + t313 * t521) * t548) * pkin(2) + (-t442 / 0.2e1 + ((t211 * t532 + t212 * t530) * pkin(1) + (t209 * t533 + t187 * t210 / 0.2e1 - t474 / 0.4e1 + t472 / 0.4e1) * pkin(3) + (t213 * t534 + t182 * t214 / 0.2e1 - t479 / 0.4e1 + t477 / 0.4e1) * pkin(2)) * rSges(6,2) + ((t201 * t532 + t202 * t531) * pkin(1) + (t199 * t533 - t187 * t200 / 0.2e1 - t475 / 0.4e1 - t473 / 0.4e1) * pkin(3) + (t203 * t534 - t182 * t204 / 0.2e1 - t480 / 0.4e1 - t478 / 0.4e1) * pkin(2)) * rSges(6,1)) * m(6), (-t510 * t541 + (rSges(3,1) * t539 - Icges(3,5)) * qJD(2)) * t314 + t335 + (-qJD(2) * (rSges(3,2) * t539 - Icges(3,6)) - t543 * t510) * t310 + t390 * t516 + ((t169 - 0.8e1 * t518) * m(6) + t191) * t415 + (t192 - 0.8e1 * t545) * t413 + t269 * t360 + (t427 + t353) * qJD(1) + t568 + ((t203 + t204) * t388 + (t537 + t538) * qJD(2) * t314 - t233 * t462 + t213 * t384 + t214 * t382 - t520 * t510 * t310 - t226 * t460 - t261 * t463 - t268 * t451) * pkin(2), t335 + t555 + t557 + t558 + t585, -(-t487 / 0.4e1 + t486 / 0.4e1) * t177 + t554 + t341 + t349 + t555 - t594, -m(6) * (t175 * (pkin(1) * qJD(5) + t506 / 0.2e1) + (t265 * t430 + t269 * t431) * rSges(6,2) + (t258 * t430 + t262 * t431) * rSges(6,1)) + ((t220 - 0.8e1 * t501) * m(6) + t285) * t415 + t556 + qJD(5) * t372 + t236 * t355 + t566 + t25 + t333 + t49 + t349 + (-t458 + t544) * t211 / 0.2e1 + t571 + t578 + t582; (pkin(1) * t352 + t496 / 0.2e1 + (t453 + t490 + t190 * t233 + t238 * t226 + ((t213 / 0.2e1 - t214 / 0.2e1) * rSges(6,2) + (t203 / 0.2e1 + t204 / 0.2e1) * rSges(6,1)) * m(6)) * pkin(2) - t353) * qJD(1) + (t512 / 0.2e1 + t516 / 0.2e1) * m(6) * t507 + t334 + t578, -t155 + t85 - t48 * t505 + ((-t190 * t267 - t238 * t260) * t289 - t65 * qJD(3)) * pkin(2) + (-t442 + ((t472 / 0.2e1 - t474 / 0.2e1) * pkin(3) + (t477 / 0.2e1 - t479 / 0.2e1) * pkin(2)) * rSges(6,2) + ((-t473 / 0.2e1 - t475 / 0.2e1) * pkin(3) + (-t478 / 0.2e1 - t480 / 0.2e1) * pkin(2)) * rSges(6,1)) * m(6), -pkin(3) * t438 - t536 * t312 + t344 + (-(t183 * t290 - t19 * t308 + t256 * t497) * t309 - t502 * t312 - (t290 * t548 + t440) * t313) * pkin(2), (-pkin(3) * t500 + (t19 * t309 - t256 * t498) * pkin(2)) * t308 + (-pkin(3) * t19 + (-t309 * t500 - t502) * pkin(2)) * t312 + t343, t336 + ((t231 * t532 + t232 * t531 - t574) * rSges(6,2) + (t224 * t532 + t225 * t530 - t573) * rSges(6,1)) * t550; t334 + t553 - t585, (-t190 * t308 + t183 + t497) * t444 + (t126 - t536) * t312 + (t443 - t505) * t499 + t443 * t548 + t344, (-t23 * t312 - t438) * pkin(3) + t350, (-t19 * t312 - t440) * pkin(3) + t350, t336; t553 - (-t205 / 0.4e1 - t198 / 0.4e1) * t177 * qJD(5) + t340 + t582 + t591 + t594, (t70 * t526 + (-t190 * t309 + t498) * t508) * t308 + (t190 * t526 + t70 * t444 + t126) * t312 + t343, (t488 + t499) * t526 + t350, -0.2e1 * qJD(5) * (Icges(6,4) / 0.2e1 + ((t317 + t316) * t307 + t240 * t311) * t311 + (((t528 - t329 / 0.2e1) * t307 + t318 / 0.2e1) * t311 + (pkin(4) * t307 / 0.2e1 - rSges(6,2) / 0.2e1) * rSges(6,1)) * m(6)), t339; (t193 + 0.8e1 * t545) * t413 + t559 + rSges(6,2) * t236 * t392 + t290 * t372 + t24 + t59 - t333 + t341 + t49 + 0.8e1 * ((-t445 + t518) * m(6) + t494) * t415 + ((-t512 - t516) * t508 + (t525 + (t513 / 0.2e1 + t517 / 0.2e1) * pkin(3)) * qJD(1)) * t551 + t568 + t570 + t591, (((t232 / 0.2e1 + t231 / 0.2e1) * qJD(2) + t574) * rSges(6,2) + ((t224 / 0.2e1 - t225 / 0.2e1) * qJD(2) + t573) * rSges(6,1)) * t550 + t562, t562, t338, 0;];
Cq = t1;
