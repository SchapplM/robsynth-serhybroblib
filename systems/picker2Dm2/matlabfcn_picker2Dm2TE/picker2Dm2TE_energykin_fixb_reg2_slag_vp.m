% Calculate inertial parameters regressor of fixed base kinetic energy for
% picker2Dm2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% T_reg [1x(2*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = picker2Dm2TE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_energykin_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'picker2Dm2TE_energykin_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 11:52:19
% EndTime: 2020-05-09 11:52:40
% DurationCPUTime: 20.41s
% Computational Cost: add. (279787->647), mult. (868886->1098), div. (3900->23), fcn. (151278->10), ass. (0->481)
t472 = 4 * pkin(1);
t317 = (pkin(3) ^ 2);
t590 = -2 * t317;
t324 = (pkin(7) ^ 2);
t308 = -2 * t324;
t275 = cos(pkin(8));
t280 = cos(qJ(1));
t540 = sin(pkin(8));
t412 = t540 * t280;
t278 = sin(qJ(1));
t482 = qJD(1) * t278;
t153 = -qJD(1) * t412 + t275 * t482;
t166 = t275 * t280 + t278 * t540;
t154 = t166 * qJD(1);
t164 = t275 * t278 - t412;
t159 = t164 ^ 2;
t337 = t166 ^ 2;
t571 = 0.1e1 / t337;
t82 = (t153 * t164 * t571 + t154 / t166) / (t159 * t571 + 0.1e1);
t281 = cos(pkin(9));
t279 = cos(qJ(2));
t520 = t278 * t279;
t455 = pkin(3) * t520;
t405 = pkin(1) * t455;
t195 = -0.2e1 * t405;
t277 = sin(qJ(2));
t569 = 0.2e1 * t277;
t476 = pkin(7) * t569;
t215 = pkin(3) * t476;
t322 = pkin(1) ^ 2;
t237 = t322 + t324;
t421 = t317 + t237;
t231 = pkin(3) * t277;
t213 = t231 + pkin(7);
t531 = t213 * t280;
t457 = pkin(1) * t531;
t120 = t195 + t215 + t421 + 0.2e1 * t457;
t118 = 0.1e1 / t120;
t318 = 0.1e1 / pkin(3);
t554 = t166 / 0.2e1;
t555 = t164 / 0.2e1;
t294 = 2 * t317;
t300 = 3 * t322;
t312 = (pkin(4) ^ 2);
t484 = t324 - t312;
t208 = t294 + t300 + t484;
t393 = -0.4e1 * t405;
t139 = t208 + t215 + t393;
t518 = t279 * t280;
t149 = pkin(3) * t518 + t213 * t278;
t194 = -pkin(1) + t455;
t220 = t322 + t484;
t242 = t277 ^ 2;
t177 = t215 + t220;
t143 = t195 + t177;
t582 = 4 * t322;
t247 = t317 * t582;
t511 = t317 * t324;
t218 = t247 - 4 * t511;
t320 = t322 ^ 2;
t439 = -0.4e1 * pkin(3) * pkin(7) * t220;
t488 = -t317 + t324;
t525 = t242 * t317;
t168 = t215 + t488 + 0.2e1 * t525;
t246 = t280 ^ 2;
t532 = t168 * t246;
t585 = 0.2e1 * pkin(3);
t86 = t218 * t242 + t277 * t439 - t320 - (t324 - (t585 + pkin(4)) * pkin(4)) * (t324 + (t585 - pkin(4)) * pkin(4)) + (t308 + (2 * t312) - (4 * t317) - 0.4e1 * t532) * t322 + (-t143 * t531 + t177 * t455) * t472;
t325 = sqrt(t86);
t73 = -t139 * t531 + t149 * t325 + (t194 * t476 + t208 * t520) * pkin(3) + (-0.2e1 * t532 + (0.2e1 * t242 - 0.4e1) * t317 - t220) * pkin(1);
t176 = t294 + t177;
t233 = pkin(1) * t280;
t468 = 0.2e1 * t233;
t105 = t168 * t468 + t176 * t213;
t109 = t176 * t280 + (0.4e1 * t246 - 0.2e1) * t213 * pkin(1);
t147 = -t194 + t531;
t232 = pkin(3) * t279;
t76 = t105 * t278 + t109 * t232 + t147 * t325;
t357 = t554 * t76 + t555 * t73;
t119 = 0.1e1 / t120 ^ 2;
t466 = 0.2e1 * t232;
t438 = pkin(7) * t466;
t206 = qJD(2) * t438;
t521 = t277 * t278;
t564 = pkin(1) * pkin(3);
t406 = t521 * t564;
t381 = qJD(2) * t406;
t480 = qJD(1) * t280;
t414 = t279 * t480;
t382 = t414 * t564;
t587 = -t382 + t381;
t501 = 0.2e1 * t587;
t430 = t206 + t501;
t478 = qJD(2) * t279;
t445 = pkin(3) * t478;
t473 = 2 * pkin(1);
t538 = t119 * ((-t213 * t482 + t280 * t445) * t473 + t430);
t557 = -t153 / 0.2e1;
t562 = t73 / 0.2e1;
t477 = qJD(2) * t317;
t386 = t277 * t279 * t477;
t378 = 0.4e1 * t386;
t155 = t206 + t378;
t418 = t278 * t480;
t390 = t168 * t418;
t533 = t155 * t246;
t561 = pkin(7) * t278;
t541 = ((0.8e1 * t390 - 0.4e1 * t533) * t322 + (t218 * t569 + t439) * t478 + (0.2e1 * t279 ^ 2 * t477 * t561 + (t143 * t482 - t280 * t430) * t213 + (t177 * t414 + (-t143 * t518 - t177 * t521) * qJD(2)) * pkin(3)) * t472) / t325;
t440 = t541 / 0.2e1;
t481 = qJD(1) * t279;
t483 = 0.2e1 * pkin(7);
t517 = t280 * t325;
t519 = t278 * t325;
t550 = pkin(1) * t246;
t560 = pkin(7) * t280;
t59 = t147 * t440 + (t105 * t280 - t213 * t519) * qJD(1) + (t155 * t280 - t168 * t482) * t278 * t473 + ((-t517 + (-t176 - 0.8e1 * t457) * t278) * t481 + ((-t109 + t519) * t277 + (t517 + (t213 * t483 + t176) * t278 + (-pkin(1) + t560 + 0.2e1 * t550) * t279 * t585) * t279) * qJD(2)) * pkin(3);
t479 = qJD(2) * t277;
t416 = t278 * t479;
t352 = t414 - t416;
t347 = 0.4e1 * t352;
t345 = t347 * t233;
t367 = t277 * t280 - t520;
t60 = t149 * t440 + (-t206 * t280 + (t139 * t278 + t517) * qJD(1)) * t213 + (t378 + 0.4e1 * t390 - 0.2e1 * t533) * pkin(1) + ((t208 * t280 - t519) * t481 + (-t139 * t518 - t208 * t521 - t325 * t367) * qJD(2) + (t194 * t478 + t231 * t352) * t483 + t213 * t345) * pkin(3);
t37 = (-t357 * t538 + (t154 * t562 + t554 * t59 + t555 * t60 + t557 * t76) * t118) * t318;
t556 = -t164 / 0.2e1;
t356 = t554 * t73 + t556 * t76;
t38 = (-t356 * t538 + (t60 * t554 + t73 * t557 + t59 * t556 - t76 * t154 / 0.2e1) * t118) * t318;
t551 = sin(pkin(9));
t537 = t118 * t318;
t350 = t356 * t537;
t67 = t357 * t537;
t56 = t281 * t350 + t551 * t67;
t54 = t281 * t67 - t350 * t551;
t570 = 0.1e1 / t56 ^ 2;
t576 = 0.1e1 / (t54 ^ 2 * t570 + 0.1e1);
t589 = t576 * (-t281 * t37 + t38 * t551) / t56;
t588 = -0.6e1 * t382 + 0.6e1 * t381;
t214 = t233 + pkin(7);
t180 = t231 + t214;
t265 = 0.4e1 / 0.3e1 * t317;
t259 = -t312 / 0.3e1;
t425 = t259 + t237;
t197 = t265 + t425;
t299 = 6 * t322;
t261 = -0.2e1 / 0.3e1 * t312;
t266 = 0.2e1 / 0.3e1 * t317;
t307 = 2 * t324;
t424 = t261 + t266 + t307;
t330 = t317 ^ 2;
t423 = t261 + t237;
t502 = (t266 + t423) * t237 + t330;
t101 = t197 * t393 + (t299 + t424) * t317 + t502;
t269 = -t317 / 0.3e1;
t227 = t269 + t324;
t170 = t227 * t195;
t516 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t122 = t197 * t516 + t170;
t248 = 0.10e2 / 0.3e1 * t322;
t123 = (t248 + t424) * t317 + t502;
t216 = pkin(7) * t468;
t239 = -3 * t322 + t324;
t509 = t322 * t246;
t459 = 0.4e1 * t509;
t169 = t216 + t459 + t239;
t236 = -3 * t317 + t324;
t245 = t280 * t246;
t326 = pkin(1) * t322;
t523 = t245 * t326;
t464 = pkin(7) * t523;
t413 = 0.8e1 * t464;
t185 = t236 * t413;
t209 = -t312 + t421;
t230 = t237 ^ 2;
t282 = 15 * t320;
t289 = 18 * t324;
t290 = -2 * t312;
t292 = -6 * t312;
t323 = t324 ^ 2;
t302 = 3 * t323;
t329 = pkin(3) * t317;
t314 = t329 ^ 2;
t450 = 0.12e2 * t509;
t451 = 0.12e2 * t525;
t471 = 6 * pkin(1);
t306 = 3 * t324;
t496 = 15 * t322 + t306;
t241 = t277 * t242;
t526 = t241 * t329;
t499 = t322 / 0.3e1 + t324;
t146 = -0.4e1 / 0.9e1 * t405 + 0.4e1 / 0.9e1 * t317 - t312 / 0.9e1 + t499;
t257 = -t312 / 0.6e1;
t371 = t324 - t405;
t157 = t257 + t266 + t371;
t229 = -t322 / 0.3e1 + t324;
t410 = -0.2e1 * t455;
t475 = 0.4e1 * t560;
t515 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t88 = 0.4e1 * t464 + 0.6e1 * t146 * t509 + t197 * t515 + (t157 * t475 + t229 * t410) * pkin(1);
t151 = t195 + t197;
t217 = t488 * t582;
t467 = 0.4e1 * t233;
t90 = pkin(7) * t151 * t467 + t217 * t246 + t101;
t70 = t88 * t451 + t185 + t122 * t450 + t314 + ((-t312 + t317 + t496) * t330) + ((t282 + (t289 + t292 + 6 * t317) * t322 + t302 + (t290 + t294) * t324) * t317) + (t230 * t209) + (0.8e1 * t169 * t526 + 0.6e1 * t231 * t90) * t214 + (t101 * t560 - t123 * t455) * t471;
t492 = t300 + t324;
t221 = t317 + t492;
t184 = t221 * t232;
t547 = t278 * pkin(1);
t188 = t232 - t547;
t454 = t322 * t232;
t222 = 3 * t317 + t237;
t529 = t222 * t278;
t100 = -0.2e1 * t246 * t454 + t184 + (0.2e1 * t188 * t560 - t529) * pkin(1);
t358 = -t278 * t326 + t454;
t178 = 0.2e1 * t358;
t235 = t317 + 2 * t322;
t238 = -t322 + t324;
t103 = t238 * t232 + t178 * t246 + (t235 * t278 + t280 * t438) * pkin(1);
t138 = -pkin(1) * t529 + t184;
t566 = 4 * t320;
t581 = 8 * t322;
t212 = pkin(3) * t566 + t329 * t581;
t433 = t278 * t516;
t565 = 4 * t326;
t140 = t212 * t279 + t433 * t565;
t298 = 5 * t320;
t485 = t323 + t330;
t284 = 10 * t322;
t495 = t284 + t307;
t507 = t324 * t322;
t163 = t317 * t495 + t298 + t485 + 6 * t507;
t293 = 5 * t330;
t304 = 6 * t324;
t174 = t293 + (t284 + t304) * t317 + t230;
t435 = t214 * t526;
t402 = -0.8e1 * t435;
t463 = -0.4e1 * t525;
t530 = t214 * t277;
t77 = t103 * t463 + t140 * t246 + (-0.4e1 * t100 * t530 + (-t163 + t413) * t279) * pkin(3) + (-0.4e1 * t138 * t560 + (t174 + t402) * t278) * pkin(1);
t64 = t180 * t70 + t325 * t77;
t586 = t64 ^ 2;
t311 = t312 ^ 2;
t487 = t320 + t323;
t491 = t307 - t312;
t508 = t324 * t312;
t355 = (t491 * t322) + t311 / 0.6e1 + t487 - t508;
t351 = 0.5e1 / 0.6e1 * t330 + t355;
t260 = -t312 / 0.2e1;
t202 = t260 + t421;
t411 = -0.4e1 * t455;
t384 = t202 * t411;
t372 = pkin(1) * t384;
t104 = t372 + ((t299 + t491) * t317) + t351;
t124 = t202 * t516 + t170;
t469 = pkin(7) * t233;
t437 = 0.6e1 * t469;
t131 = (t248 + t491) * t317 + t351;
t262 = -0.3e1 / 0.2e1 * t312;
t500 = t311 / 0.2e1 - t330 / 0.2e1;
t391 = -(3 * t508) + t302 + t500;
t504 = t237 * ((t262 + t307) * t322 - 0.3e1 / 0.2e1 * t508 + t487 + t500) + t314;
t87 = -0.6e1 * t131 * t405 + (t282 + ((t289 - 9 * t312) * t322) + t391) * t317 + (t262 + t496) * t330 + t504;
t79 = t104 * t437 + t124 * t450 + t185 + t87;
t584 = 0.8e1 * t79;
t267 = t317 / 0.3e1;
t427 = t312 / 0.3e1 + t267 + t307;
t583 = 0.10e2 / 0.3e1 * t330 + (-t322 + t427) * t590 + (t269 + t425) * t308;
t258 = -t312 / 0.4e1;
t580 = t258 + t317 / 0.2e1;
t579 = t306 - t312 - t317;
t383 = qJD(1) * t406;
t417 = t214 * t478;
t396 = pkin(3) * t417;
t578 = t383 - t396;
t577 = qJD(1) - qJD(2);
t573 = -0.8e1 * pkin(7);
t572 = -4 * pkin(1);
t568 = -0.8e1 * t280;
t567 = -0.2e1 * t280;
t62 = 0.1e1 / t64;
t563 = t62 / 0.4e1;
t559 = t118 / 0.2e1;
t249 = -0.20e2 / 0.3e1 * t322;
t305 = 4 * t324;
t428 = 0.2e1 / 0.3e1 * t312 + t266 + t305;
t429 = 0.4e1 / 0.3e1 * t312 + t265 + t308;
t558 = t330 / 0.2e1 - (t249 + t428) * t317 / 0.2e1 + 0.3e1 / 0.2e1 * t320 - t429 * t322 / 0.2e1 - t323 / 0.2e1;
t553 = -t325 / 0.4e1;
t552 = t325 / 0.4e1;
t549 = pkin(3) * t318;
t349 = t352 * pkin(3);
t128 = t322 * t349 * t475;
t448 = pkin(1) * t482;
t175 = t445 - t448;
t404 = pkin(7) * t448;
t207 = -0.2e1 * t404;
t243 = t278 ^ 2;
t370 = t242 * t329 * t417;
t359 = -0.24e2 * t370;
t362 = t229 * t381;
t363 = t229 * t382;
t522 = t246 * t326;
t465 = pkin(7) * t522;
t379 = t465 * t482;
t364 = -0.24e2 * t379;
t385 = t322 * t418;
t374 = -0.24e2 * t385;
t446 = pkin(3) * t479;
t375 = t446 * t523;
t420 = qJD(1) * t526;
t397 = pkin(1) * t420;
t377 = t278 * t397;
t388 = t480 * t516;
t392 = -0.6e1 * t404;
t447 = pkin(1) * t480;
t394 = t222 * t447;
t401 = -0.2e1 * t418;
t407 = t509 * t564;
t415 = t280 * t479;
t456 = pkin(3) * t530;
t503 = 0.4e1 * t587 * t197;
t344 = pkin(1) * t227 * t349 * t509;
t505 = t236 * t364 - 0.24e2 * t344;
t545 = ((t214 * t397 * t568 + t243 * t420 * t581 + t359 * t547 + (0.2e1 * (-t322 * t446 - t326 * t480) * t246 + t178 * t401 - t238 * t446 + (t235 * t480 + (-t278 * t481 - t415) * pkin(7) * t585) * pkin(1)) * t463 - 0.8e1 * t103 * t386 - 0.4e1 * (-t394 + (0.4e1 * t279 * t385 + (-t221 + 0.2e1 * t509) * t479) * pkin(3) + ((-t446 - t447) * t233 - t188 * t448) * t483) * t456 + t375 * t573 + t364 * t232 + (-t212 * t479 + t388 * t565) * t246 + t140 * t401 - 0.4e1 * (-t221 * t446 - t394) * t469 + 0.4e1 * t138 * t404 + t163 * t446 + t174 * t447 + 0.4e1 * t578 * t100) * t325 + t77 * t440 + (0.8e1 * (t207 - 0.8e1 * t385) * t435 + (-0.12e2 * t379 + 0.6e1 * (-0.4e1 / 0.9e1 * t414 + 0.4e1 / 0.9e1 * t416) * t407 - 0.12e2 * t146 * t385 - t128 - 0.4e1 * t157 * t404 - 0.2e1 * t363 + 0.2e1 * t362) * t451 + 0.24e2 * t88 * t386 + 0.6e1 * (t217 * t401 + (-t151 * t482 + t280 * t501) * pkin(7) * t472 + t503) * t456 + t122 * t374 + t503 * t437 + t101 * t392 + t505 - 0.6e1 * t578 * t90 + (-0.8e1 * t377 + 0.24e2 * t370) * t169 + t588 * t123) * t180 + t70 * t175) / t586;
t334 = pkin(7) * t324;
t219 = -0.12e2 * pkin(7) * t326 + t334 * t472;
t226 = -8 * t320 + 12 * t507;
t409 = 0.16e2 * t464;
t244 = t246 ^ 2;
t510 = t320 * t244;
t462 = 0.8e1 * t510;
t111 = t219 * t280 + t226 * t246 + t409 + t462 + t487 - (6 * t507);
t126 = t195 * t516 + t202 * t236;
t186 = 16 * (t485 - 6 * t511) * t320;
t234 = -30 * t312 + 60 * t324;
t240 = t242 ^ 2;
t490 = t311 - t330;
t361 = 6 * t323 + t490 - 6 * t508;
t408 = 0.32e2 * t464;
t528 = t230 * (-t317 + t220);
t152 = -t330 / 0.6e1 + t355;
t273 = t322 / 0.2e1;
t498 = t273 + t324;
t156 = -0.2e1 / 0.3e1 * t405 + t258 + t498;
t224 = (t305 + t312) * t322;
t270 = -0.2e1 / 0.3e1 * t317;
t228 = t270 + t324;
t474 = 0.6e1 * t560;
t191 = t324 + t317 / 0.4e1 + t322 / 0.4e1 - t312 / 0.8e1;
t497 = 0.4e1 / 0.7e1 * t324 - t312 / 0.7e1;
t97 = -0.32e2 / 0.21e2 * t191 * t405 + 0.5e1 / 0.42e2 * t330 + (0.16e2 / 0.21e2 * t322 + t497) * t317 + t320 / 0.7e1 + t497 * t322 + t323 - 0.3e1 / 0.7e1 * t508 + t311 / 0.42e2;
t193 = t499 + t580;
t271 = 0.4e1 / 0.3e1 * t322;
t99 = -0.8e1 / 0.3e1 * t193 * t405 + 0.5e1 / 0.18e2 * t330 + (t271 + t259) * t317 + t323 - t320 / 0.3e1 + t311 / 0.18e2 + (t265 + 0.2e1 / 0.3e1 * t322 + t261) * t324;
t78 = t228 * t462 + t156 * t409 + 0.14e2 * t97 * t509 + (t238 * t330) + (t224 - 0.10e2 / 0.3e1 * t320 + (2 * t323) - t508) * t317 + t152 * t515 + (t229 * t384 + t474 * t99) * pkin(1);
t303 = 8 * t324;
t144 = t326 * t411 + t247 + t566 + ((t290 + t303) * t322);
t150 = -t322 + t371 + t580;
t89 = t413 + t144 * t246 + t202 * t239 + (t150 * t475 + t410 * t515) * pkin(1);
t486 = t323 - t320;
t91 = t227 * t372 - t314 + (-t248 - t484) * t330 + (t224 + t330 / 0.6e1 - t311 / 0.6e1 + t486) * t317 + t152 * t324;
t291 = -5 * t312;
t297 = 7 * t320;
t95 = (t262 + t306 + (7 * t322)) * t330 + (t297 + ((t291 + 10 * t324) * t322) + t391) * t317 + t504;
t58 = t186 * t244 + t126 * t408 + 0.24e2 * t91 * t509 + (t290 + t305 + 28 * t322) * t314 + (t209 * t528) + ((t234 * t320) + 0.24e2 * t78 * t242 + (t292 * t323) + (t299 * t361) + (t307 * t490) + (28 * t326 ^ 2) + 0.4e1 * t334 ^ 2) * t317 + (t231 * t584 + 0.32e2 * t526 * t89) * t214 + 0.8e1 * (-t455 * t95 + t560 * t87) * pkin(1) + (0.16e2 * t111 * t240 + (t234 * t322) + (70 * t320) + t330 + t361) * t330;
t205 = t579 * t284;
t493 = t291 - 5 * t317;
t102 = 0.7e1 * t314 + ((35 * t322 + 15 * t324 + t493) * t330) + ((21 * t320 + t205 + 9 * t323 + (t292 - 6 * t317) * t324) * t317) + t528;
t190 = t324 + 0.5e1 / 0.2e1 * t317 + 0.3e1 / 0.2e1 * t322 + t260;
t444 = (pkin(1) * t236) / 0.2e1;
t135 = t190 * t232 + t278 * t444;
t158 = 0.4e1 / 0.3e1 * t509 + t216 + t229;
t527 = t240 * t330;
t400 = -0.24e2 * t158 * t527;
t452 = -0.12e2 * t525;
t461 = -0.6e1 * t509;
t189 = 0.7e1 / 0.6e1 * t317 + t257 + t498;
t426 = t257 + t267 + t324;
t192 = t271 + t426;
t121 = -t189 * t547 + t192 * t232;
t127 = (t317 * t238) - 0.5e1 / 0.3e1 * t320 + t427 * t322 + t324 * (t259 + t227);
t199 = t322 + t426;
t223 = t294 + t238;
t133 = t199 * t232 - t223 * t547 / 0.2e1;
t524 = t245 * t320;
t80 = -0.4e1 * t524 * t561 + t121 * t459 + (-0.8e1 / 0.3e1 * t510 + t127) * t232 + (t133 * t475 + t278 * t558) * pkin(1);
t494 = t290 + t590;
t422 = t304 + t494;
t115 = t330 + (t261 + t270 + t495) * t317 + t298 + (t422 * t322) + t324 * (t261 + t228);
t108 = t115 * t232;
t196 = 0.8e1 / 0.3e1 * t317 + t425;
t198 = t259 + t266 + t492;
t125 = -t196 * t547 + t198 * t232;
t203 = 0.5e1 / 0.6e1 * t317 + t273 + t257;
t137 = pkin(1) * t433 + t203 * t466;
t453 = t326 * t232;
t460 = -0.4e1 * t509;
t132 = t293 + ((t284 + t422) * t317) + (t270 + t423) * t237;
t536 = t132 * t278;
t83 = t245 * t453 * t573 + t137 * t460 + t108 + (t125 * t475 - t536) * pkin(1);
t92 = -pkin(1) * t536 + t108;
t130 = -(3 * t330) + (t249 + t429) * t317 + t428 * t322 + t486;
t93 = t130 * t232 + t547 * t583;
t94 = t314 + ((21 * t322 + t579) * t330) + ((t324 * t494 + t205 + t302 + 35 * t320) * t317) + ((t297 + (t303 + t493) * t322 + t324 * (-t317 + t484)) * t237);
t179 = 0.4e1 * t358;
t187 = t232 + 0.2e1 * t547;
t201 = t260 + t221;
t96 = t239 * t232 + t179 * t246 + (t187 * t560 + t201 * t278) * t473;
t65 = t135 * t409 + t96 * t402 + t80 * t452 + t93 * t461 + (-0.6e1 * t83 * t530 + (0.24e2 * t227 * t510 - t94) * t279) * pkin(3) + (-0.6e1 * t92 * t560 + (t102 + t400) * t278) * pkin(1);
t48 = t180 * t58 + t325 * t65;
t543 = 0.1e1 / t48 ^ 2 * t586;
t542 = t48 * t62;
t432 = -t538 / 0.2e1;
t343 = t73 ^ 2;
t72 = 0.1e1 / t343;
t75 = t76 ^ 2;
t33 = qJD(1) + ((t432 * t76 + t559 * t59) / t562 - 0.2e1 * (t432 * t73 + t559 * t60) * t76 * t72) * t120 / (t72 * t75 + 0.1e1) * t549;
t539 = pkin(1) * qJD(1);
t449 = t118 * t539;
t398 = t318 * t449;
t71 = t398 * t562;
t32 = pkin(2) * t33 + t71;
t373 = -t398 / 0.2e1;
t360 = t76 * t373;
t49 = t54 * t360;
t19 = -t32 * t56 + t49;
t313 = 0.1e1 / pkin(4);
t319 = 0.1e1 / pkin(3) ^ 2;
t380 = t319 * t76 * t449;
t431 = (pkin(3) * t33 + t71) * t318 / 0.2e1;
t512 = t313 * t325;
t18 = t313 * t431 * t542 - t380 * t512 / 0.4e1;
t309 = qJD(1) ^ 2;
t514 = t309 * t322;
t513 = t313 * t319;
t506 = t588 * t131;
t81 = qJD(1) - t82;
t458 = t81 * t539;
t443 = -t545 / 0.4e1;
t442 = t73 * t563;
t441 = t76 * t563;
t436 = t118 * t513;
t434 = t278 * t527;
t419 = t245 * t482;
t399 = t33 + t589;
t395 = pkin(3) * t415;
t389 = t320 * t419;
t387 = t241 * t330 * t478;
t346 = t352 * t473;
t348 = -t115 * t446 - t132 * t447;
t365 = -0.48e2 * t379;
t376 = t446 * t510;
t26 = (t400 * t447 - 0.24e2 * pkin(1) * (-0.8e1 / 0.3e1 * t385 + t207) * t434 - 0.96e2 * t158 * t387 * t547 + ((-t239 + t460) * t446 + 0.2e1 * (pkin(1) * t201 - t179 * t278 - 0.2e1 * t522) * t480 + (-0.2e1 * t395 + (-0.2e1 * t187 * t278 + 0.4e1 * t550) * qJD(1)) * pkin(1) * pkin(7)) * t402 + (0.8e1 / 0.3e1 * t376 + (-t189 * t447 - t192 * t446) * t459 - t127 * t446 + t447 * t558 + (0.32e2 / 0.3e1 * t524 * t232 + t322 * t121 * t568) * t482 + (t199 * t395 * t572 + ((0.12e2 * t243 * t246 - 0.4e1 * t244) * t320 + (-0.4e1 * t133 * t278 - 0.2e1 * t223 * t550) * pkin(1)) * qJD(1)) * pkin(7)) * t452 - 0.24e2 * t80 * t386 - 0.6e1 * ((-0.4e1 * (pkin(1) * t388 - 0.2e1 * t203 * t446) * t246 + 0.8e1 * t137 * t418) * t322 + (0.8e1 * t375 + (-t196 * t447 - t198 * t446) * t467 + (t125 * t572 + 0.24e2 * t246 * t453) * t482) * pkin(7) + t348) * t456 + (-t190 * t446 + t444 * t480) * t409 + t135 * t365 + (-t130 * t446 + t447 * t583) * t461 + 0.12e2 * t93 * t385 - 0.6e1 * t348 * t469 + 0.6e1 * t92 * t404 + t94 * t446 + t102 * t447 + (0.8e1 * t377 + t359) * t96 + 0.6e1 * t578 * t83 + (-0.96e2 * t232 * t389 - 0.24e2 * t376) * t227) * t325 + t65 * t440 + (0.16e2 * (t226 * t567 - t219 - 0.48e2 * t465 - 0.32e2 * t524) * qJD(1) * t434 + 0.64e2 * t111 * t387 + 0.32e2 * (-t128 + (t144 * t567 + (t150 * t572 - 0.24e2 * t522) * pkin(7)) * t482 + (-t346 * t515 - t347 * t522) * pkin(3)) * t435 + 0.24e2 * (-0.32e2 * t228 * t389 + (-0.2e1 / 0.3e1 * t414 + 0.2e1 / 0.3e1 * t416) * t409 * t564 + t156 * t365 - 0.28e2 * t97 * t385 + (-0.8e1 / 0.3e1 * t414 + 0.8e1 / 0.3e1 * t416) * t193 * pkin(3) * t322 * t474 + t99 * t392 + 0.4e1 * (-t363 + t362) * t202 - 0.64e2 / 0.3e1 * t352 * t191 * t407) * t525 + 0.48e2 * t78 * t386 - 0.8e1 * t79 * t383 + 0.8e1 * (t124 * t374 + (-pkin(3) * t202 * t345 - t104 * t482) * pkin(7) * t471 + t505 + t506) * t456 + t396 * t584 - 0.4e1 * t186 * t419 - pkin(3) * t346 * t408 * t516 - 0.96e2 * t126 * t379 - 0.96e2 * t202 * t344 - 0.48e2 * t91 * t385 + 0.8e1 * t506 * t469 - 0.8e1 * t87 * t404 + 0.8e1 * t587 * t95 + (-0.32e2 * t377 + 0.96e2 * t370) * t89) * t180 + t58 * t175;
t9 = t33 + (0.1e1 / t48 * t64 * t440 - 0.2e1 * (t26 * t62 / 0.2e1 - t48 * t545 / 0.2e1) * pkin(4) * t512 * t543 * t549) / (t543 * t86 + 0.1e1);
t366 = t518 + t521;
t354 = t441 * t48 + t552 * t73;
t353 = t442 * t48 + t553 * t76;
t114 = t577 * t366;
t113 = t577 * t367;
t50 = t56 * t360;
t40 = t54 * t71 + t50;
t39 = t373 * t56 * t73 + t49;
t36 = t354 * t436;
t35 = t353 * t436;
t24 = -t35 * t367 + t36 * t366;
t23 = -t35 * t366 - t36 * t367;
t22 = 0.1e1 / t24 ^ 2;
t20 = t32 * t54 + t50;
t17 = (t380 * t542 / 0.4e1 + t325 * t431) * t313;
t13 = (-t354 * t538 + (t60 * t552 + t73 * t541 / 0.8e1 + t26 * t441 + (t443 * t76 + t563 * t59) * t48) * t118) * t513;
t12 = (-t353 * t538 + (t26 * t442 + t59 * t553 - t76 * t541 / 0.8e1 + (t443 * t73 + t563 * t60) * t48) * t118) * t513;
t11 = (t281 * t38 + t37 * t551) * t54 * t570 * t576 + t399;
t10 = t11 ^ 2 / 0.2e1;
t8 = pkin(6) * t11 + t19;
t7 = pkin(4) * t9 + t18;
t6 = t20 * t56 + t54 * t8;
t5 = -t20 * t54 + t56 * t8;
t4 = t399 - t589;
t3 = t17 * t24 - t23 * t7;
t2 = t17 * t23 + t24 * t7;
t1 = (-(-t113 * t35 + t114 * t36 - t12 * t366 - t13 * t367) / t24 - (-t113 * t36 - t114 * t35 + t12 * t367 - t13 * t366) * t23 * t22) / (t22 * t23 ^ 2 + 0.1e1) + t9;
t14 = [0, 0, 0, 0, 0, t309 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 / 0.2e1, t33 * t71, t33 * t360, 0, (t75 / 0.8e1 + t343 / 0.8e1) * t319 * t119 * t514, 0, 0, 0, 0, 0, t10, t19 * t11, -t20 * t11, 0, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t9 ^ 2 / 0.2e1, t18 * t9, -t17 * t9, 0, t17 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t82 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t39 * t11, -t40 * t11, 0, t40 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 ^ 2 / 0.2e1, -t166 * t458, -t164 * t458, 0, (t159 / 0.2e1 + t337 / 0.2e1) * t514, 0, 0, 0, 0, 0, t4 ^ 2 / 0.2e1, t5 * t4, -t6 * t4, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t1 ^ 2 / 0.2e1, t2 * t1, -t3 * t1, 0, t3 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t14;
