% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1turnOL_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:48
% EndTime: 2020-04-12 19:41:04
% DurationCPUTime: 12.02s
% Computational Cost: add. (5978->575), mult. (10706->863), div. (0->0), fcn. (8532->8), ass. (0->358)
t264 = qJ(2) + qJ(3);
t248 = sin(t264);
t249 = cos(t264);
t342 = rSges(4,1) * t248 + rSges(4,2) * t249;
t267 = sin(qJ(1));
t328 = Icges(4,5) * t248 + Icges(4,6) * t249;
t163 = t328 * t267;
t270 = cos(qJ(1));
t164 = t328 * t270;
t438 = t249 * t267;
t440 = t248 * t267;
t454 = Icges(4,6) * t270;
t131 = Icges(4,4) * t438 - Icges(4,2) * t440 + t454;
t463 = Icges(4,4) * t249;
t335 = Icges(4,1) * t248 + t463;
t514 = -t335 * t267 - t131;
t332 = -Icges(4,2) * t248 + t463;
t132 = Icges(4,6) * t267 - t270 * t332;
t513 = -t335 * t270 + t132;
t464 = Icges(4,4) * t248;
t331 = Icges(4,2) * t249 + t464;
t336 = Icges(4,1) * t249 - t464;
t512 = -t331 + t336;
t511 = -t332 - t335;
t261 = qJD(2) + qJD(3);
t510 = t342 * t261;
t509 = 2 * qJD(2);
t508 = 2 * qJD(4);
t233 = rSges(4,2) * t440;
t387 = rSges(4,1) * t438;
t135 = rSges(4,3) * t270 - t233 + t387;
t209 = t261 * t267;
t210 = t261 * t270;
t478 = rSges(4,2) * t248;
t480 = rSges(4,1) * t249;
t343 = -t478 + t480;
t297 = rSges(4,3) * t267 - t270 * t343;
t474 = pkin(2) * qJD(2);
t269 = cos(qJ(2));
t402 = t267 ^ 2 + t270 ^ 2;
t506 = t269 * t402;
t50 = -t209 * t135 + t210 * t297 + t474 * t506;
t507 = qJD(2) * t50;
t268 = cos(qJ(4));
t256 = Icges(5,4) * t268;
t265 = sin(qJ(4));
t330 = -Icges(5,2) * t265 + t256;
t219 = Icges(5,1) * t265 + t256;
t257 = Icges(3,4) * t269;
t266 = sin(qJ(2));
t333 = -Icges(3,2) * t266 + t257;
t221 = Icges(3,1) * t266 + t257;
t169 = t342 * t267;
t170 = t342 * t270;
t351 = qJD(1) * t261;
t195 = t267 * t351;
t196 = t270 * t351;
t271 = qJD(2) ^ 2;
t362 = t266 * t402;
t395 = qJD(1) * t270;
t350 = rSges(4,3) * t395 + qJD(1) * t387 + t270 * t510;
t78 = -qJD(1) * t233 + t350;
t294 = qJD(1) * t297;
t79 = t169 * t261 + t294;
t25 = -pkin(2) * t271 * t362 - t196 * t135 - t195 * t297 + t209 * t79 + t210 * t78;
t146 = t343 * t261;
t272 = qJD(1) ^ 2;
t361 = (-t271 - t272) * t269;
t389 = qJD(1) * qJD(2);
t379 = t266 * t389;
t36 = -qJD(1) * t79 + t146 * t210 - t342 * t195 + (0.2e1 * t267 * t379 + t270 * t361) * pkin(2);
t37 = qJD(1) * t78 + t146 * t209 + t342 * t196 + (t267 * t361 - 0.2e1 * t270 * t379) * pkin(2);
t468 = t37 * t267;
t484 = pkin(2) * t269;
t236 = t395 * t484;
t385 = t266 * t474;
t63 = -t209 * t342 + t267 * t385 - t236 - t294;
t393 = qJD(2) * t270;
t396 = qJD(1) * t267;
t64 = -pkin(2) * (t266 * t393 + t269 * t396) + qJD(1) * t135 + t210 * t342;
t505 = t25 * (-t267 * t135 + t270 * t297) - (-t468 - t36 * t270 + (t267 * t64 + t270 * t63) * qJD(1)) * t342 + t63 * (qJD(1) * t170 + t209 * t343) + (-t135 * t395 + t270 * t78 + (-t294 + t79) * t267 - t209 * t169 - t170 * t210) * t50;
t214 = Icges(3,5) * t269 - Icges(3,6) * t266;
t213 = Icges(3,5) * t266 + Icges(3,6) * t269;
t306 = qJD(2) * t213;
t465 = Icges(3,4) * t266;
t222 = Icges(3,1) * t269 - t465;
t158 = Icges(3,5) * t267 + t222 * t270;
t154 = Icges(3,6) * t267 + t270 * t333;
t445 = t154 * t266;
t321 = -t158 * t269 + t445;
t452 = Icges(3,3) * t270;
t504 = -t270 * t306 + (-t214 * t267 + t321 + t452) * qJD(1);
t434 = t266 * t267;
t242 = Icges(3,4) * t434;
t431 = t267 * t269;
t461 = Icges(3,5) * t270;
t157 = Icges(3,1) * t431 - t242 - t461;
t455 = Icges(3,6) * t270;
t153 = Icges(3,4) * t431 - Icges(3,2) * t434 - t455;
t446 = t153 * t266;
t322 = -t157 * t269 + t446;
t150 = Icges(3,3) * t267 + t214 * t270;
t399 = qJD(1) * t150;
t503 = qJD(1) * t322 - t267 * t306 + t399;
t212 = Icges(5,5) * t268 - Icges(5,6) * t265;
t211 = Icges(5,5) * t265 + Icges(5,6) * t268;
t303 = qJD(4) * t211;
t462 = Icges(5,4) * t265;
t220 = Icges(5,1) * t268 - t462;
t156 = Icges(5,5) * t267 + t220 * t270;
t152 = Icges(5,6) * t267 + t270 * t330;
t447 = t152 * t265;
t323 = -t156 * t268 + t447;
t450 = Icges(5,3) * t270;
t502 = -t270 * t303 + (-t212 * t267 + t323 + t450) * qJD(1);
t436 = t265 * t267;
t241 = Icges(5,4) * t436;
t432 = t267 * t268;
t458 = Icges(5,5) * t270;
t155 = Icges(5,1) * t432 - t241 - t458;
t453 = Icges(5,6) * t270;
t151 = Icges(5,4) * t432 - Icges(5,2) * t436 - t453;
t448 = t151 * t265;
t324 = -t155 * t268 + t448;
t148 = Icges(5,3) * t267 + t212 * t270;
t400 = qJD(1) * t148;
t501 = qJD(1) * t324 - t267 * t303 + t400;
t329 = Icges(4,5) * t249 - Icges(4,6) * t248;
t439 = t248 * t270;
t232 = Icges(4,4) * t439;
t437 = t249 * t270;
t460 = Icges(4,5) * t267;
t134 = -Icges(4,1) * t437 + t232 + t460;
t449 = t134 * t249;
t451 = Icges(4,3) * t270;
t500 = -t261 * t164 + (t132 * t248 - t267 * t329 - t449 - t451) * qJD(1);
t231 = Icges(4,4) * t440;
t459 = Icges(4,5) * t270;
t133 = Icges(4,1) * t438 - t231 + t459;
t326 = t131 * t248 - t133 * t249;
t130 = Icges(4,3) * t267 - t270 * t329;
t401 = qJD(1) * t130;
t499 = qJD(1) * t326 + t163 * t261 + t401;
t149 = Icges(3,5) * t431 - Icges(3,6) * t434 - t452;
t57 = -t270 * t149 - t267 * t322;
t147 = Icges(5,5) * t432 - Icges(5,6) * t436 - t450;
t55 = -t270 * t147 - t267 * t324;
t320 = -t248 * t331 + t249 * t335;
t498 = qJD(1) * t320 + t329 * t261;
t215 = Icges(5,2) * t268 + t462;
t319 = t215 * t265 - t219 * t268;
t497 = t319 * qJD(1) + t212 * qJD(4);
t217 = Icges(3,2) * t269 + t465;
t318 = t217 * t266 - t269 * t221;
t496 = t318 * qJD(1) + t214 * qJD(2);
t495 = t267 * (-t217 * t270 + t158) - t270 * (-Icges(3,2) * t431 + t157 - t242);
t494 = t267 * (-t215 * t270 + t156) - t270 * (-Icges(5,2) * t432 + t155 - t241);
t493 = qJD(1) * t511 + t209 * (Icges(4,2) * t437 + t134 + t232) + t210 * (-Icges(4,2) * t438 + t133 - t231);
t492 = t195 / 0.2e1;
t491 = t196 / 0.2e1;
t490 = -t209 / 0.2e1;
t489 = t209 / 0.2e1;
t488 = -t210 / 0.2e1;
t487 = t210 / 0.2e1;
t486 = t267 / 0.2e1;
t485 = -t270 / 0.2e1;
t483 = -qJD(1) / 0.2e1;
t482 = qJD(1) / 0.2e1;
t479 = rSges(5,1) * t268;
t476 = rSges(5,2) * t265;
t475 = rSges(5,2) * t268;
t223 = rSges(5,1) * t265 + t475;
t185 = t223 * t270;
t258 = t267 * rSges(5,3);
t430 = t268 * t270;
t403 = rSges(5,1) * t430 + t258;
t435 = t265 * t270;
t161 = -rSges(5,2) * t435 + t403;
t392 = qJD(4) * t267;
t382 = t223 * t392;
t93 = -t382 + (pkin(1) * t270 + t161) * qJD(1);
t473 = t185 * t93;
t472 = t267 * t63;
t243 = rSges(5,2) * t436;
t404 = t270 * rSges(5,3) + t243;
t159 = rSges(5,1) * t432 - t404;
t391 = qJD(4) * t270;
t381 = t223 * t391;
t92 = -t381 + (-pkin(1) * t267 - t159) * qJD(1);
t471 = t267 * t92;
t470 = t270 * rSges(3,3);
t469 = t270 * t92;
t444 = t211 * t267;
t443 = t211 * t270;
t442 = t213 * t267;
t441 = t213 * t270;
t433 = t266 * t270;
t429 = t269 * t270;
t129 = Icges(4,5) * t438 - Icges(4,6) * t440 + t451;
t428 = t270 * t129;
t66 = t331 * t440 - t335 * t438 - t164;
t425 = t66 * qJD(1);
t88 = -t267 * t319 - t443;
t424 = t88 * qJD(1);
t89 = -t267 * t318 - t441;
t423 = t89 * qJD(1);
t422 = t267 * t129 + t131 * t439;
t421 = t267 * t130 + t132 * t439;
t420 = -t267 * t147 - t155 * t430;
t419 = t267 * t148 + t156 * t430;
t418 = -t267 * t149 - t157 * t429;
t417 = t267 * t150 + t158 * t429;
t409 = -t215 + t220;
t408 = t219 + t330;
t407 = -t217 + t222;
t406 = t221 + t333;
t405 = rSges(5,3) * t395 + qJD(1) * t243;
t344 = rSges(3,1) * t269 - rSges(3,2) * t266;
t162 = rSges(3,3) * t267 + t270 * t344;
t398 = qJD(1) * t162;
t397 = qJD(1) * t212;
t394 = qJD(2) * t267;
t390 = t214 * qJD(1);
t224 = rSges(3,1) * t266 + rSges(3,2) * t269;
t384 = t224 * t394;
t383 = t224 * t393;
t380 = -pkin(1) - t479;
t378 = t396 / 0.2e1;
t377 = t395 / 0.2e1;
t376 = -t394 / 0.2e1;
t373 = t393 / 0.2e1;
t372 = -t392 / 0.2e1;
t369 = t391 / 0.2e1;
t366 = -(-t270 * t336 + t460) * qJD(1) + t514 * t261;
t365 = -(t267 * t336 + t459) * qJD(1) + t513 * t261;
t310 = t261 * t331;
t364 = -qJD(1) * t132 + t133 * t261 - t267 * t310;
t363 = t134 * t261 + t270 * t310 + (t267 * t332 + t454) * qJD(1);
t112 = t132 * t440;
t360 = t270 * t130 - t112;
t121 = t156 * t432;
t359 = t270 * t148 - t121;
t122 = t158 * t431;
t358 = t270 * t150 - t122;
t357 = t129 + t449;
t356 = t511 * t261;
t355 = t512 * t261;
t354 = -t147 + t447;
t353 = -t149 + t445;
t352 = -qJD(1) * t169 + t210 * t343;
t341 = -t476 + t479;
t339 = -t267 * t93 - t469;
t205 = t341 * qJD(4);
t338 = -pkin(1) * t272 - qJD(4) * t205;
t110 = t384 - t398;
t160 = rSges(3,1) * t431 - rSges(3,2) * t434 - t470;
t111 = -qJD(1) * t160 - t383;
t327 = t110 * t267 - t111 * t270;
t68 = t131 * t249 + t133 * t248;
t84 = t151 * t268 + t155 * t265;
t85 = t152 * t268 + t156 * t265;
t86 = t153 * t269 + t157 * t266;
t87 = t154 * t269 + t158 * t266;
t183 = t223 * t267;
t316 = qJD(2) * t224;
t56 = -t152 * t436 - t359;
t314 = (t267 * t56 - t270 * t55) * qJD(4);
t58 = -t154 * t434 - t358;
t313 = (t267 * t58 - t270 * t57) * qJD(2);
t59 = -t151 * t435 - t420;
t60 = -t152 * t435 + t419;
t312 = (t267 * t60 - t270 * t59) * qJD(4);
t61 = -t153 * t433 - t418;
t62 = -t154 * t433 + t417;
t311 = (t267 * t62 - t270 * t61) * qJD(2);
t309 = t326 * t267;
t308 = qJD(2) * t221;
t307 = qJD(2) * t217;
t305 = qJD(4) * t219;
t304 = qJD(4) * t215;
t80 = (t159 * t267 + t161 * t270) * qJD(4);
t81 = (t160 * t267 + t162 * t270) * qJD(2);
t301 = -qJD(1) * t329 - t163 * t210 + t164 * t209;
t300 = t151 * t270 - t152 * t267;
t299 = t153 * t270 - t154 * t267;
t285 = -qJD(1) * t129 - t248 * t364 + t249 * t366;
t10 = t285 * t267 - t270 * t499;
t284 = t248 * t363 + t249 * t365 + t401;
t11 = t284 * t267 + t270 * t500;
t46 = -t309 + t428;
t47 = -t134 * t438 - t360;
t23 = t209 * t47 - t210 * t46 - t425;
t48 = t133 * t437 - t422;
t49 = -t134 * t437 + t421;
t67 = t270 * t320 - t163;
t65 = t67 * qJD(1);
t24 = t209 * t49 - t210 * t48 + t65;
t290 = qJD(1) * t512 + t209 * t513 - t210 * t514;
t274 = t248 * t493 + t290 * t249;
t283 = -qJD(1) * t328 + t248 * t356 + t249 * t355;
t32 = -t267 * t498 + t283 * t270;
t33 = t283 * t267 + t270 * t498;
t34 = t248 * t366 + t249 * t364;
t35 = t248 * t365 - t249 * t363;
t69 = -t132 * t249 - t134 * t248;
t8 = t267 * t499 + t285 * t270;
t9 = -t267 * t500 + t284 * t270;
t298 = (qJD(1) * t32 + t195 * t48 + t196 * t49 + t209 * t9 - t210 * t8) * t486 + (t290 * t248 - t249 * t493) * t483 + t23 * t378 + t24 * t377 + (qJD(1) * t33 - t10 * t210 + t11 * t209 + t195 * t46 + t196 * t47) * t485 + (t267 * t9 - t270 * t8 + (t267 * t48 + t270 * t49) * qJD(1)) * t489 + (t267 * t47 - t270 * t46) * t492 + (t267 * t49 - t270 * t48) * t491 + (-t10 * t270 + t11 * t267 + (t267 * t46 + t270 * t47) * qJD(1)) * t488 + (t267 * t35 - t270 * t34 + (t267 * t68 + t270 * t69) * qJD(1)) * t482 + (t267 * t301 + t270 * t274) * t490 + (t267 * t274 - t270 * t301) * t487;
t296 = (-t265 * t408 + t268 * t409) * qJD(1);
t295 = (-t266 * t406 + t269 * t407) * qJD(1);
t103 = qJD(1) * t156 - t267 * t305;
t99 = qJD(1) * t152 - t267 * t304;
t282 = qJD(1) * t147 - qJD(4) * t84 + t103 * t268 - t265 * t99;
t102 = -t270 * t305 + (-t220 * t267 + t458) * qJD(1);
t98 = -t270 * t304 + (-t267 * t330 + t453) * qJD(1);
t281 = -qJD(4) * t85 + t102 * t268 - t265 * t98 + t400;
t101 = qJD(1) * t154 - t267 * t307;
t105 = qJD(1) * t158 - t267 * t308;
t280 = qJD(1) * t149 - qJD(2) * t86 - t101 * t266 + t105 * t269;
t100 = -t270 * t307 + (-t267 * t333 + t455) * qJD(1);
t104 = -t270 * t308 + (-t222 * t267 + t461) * qJD(1);
t279 = -qJD(2) * t87 - t100 * t266 + t104 * t269 + t399;
t199 = t330 * qJD(4);
t201 = t220 * qJD(4);
t278 = qJD(1) * t211 - t199 * t265 + t201 * t268 + (-t215 * t268 - t219 * t265) * qJD(4);
t200 = t333 * qJD(2);
t202 = t222 * qJD(2);
t277 = qJD(1) * t213 - t200 * t266 + t202 * t269 + (-t217 * t269 - t221 * t266) * qJD(2);
t276 = -t265 * t494 + t300 * t268;
t275 = -t266 * t495 + t299 * t269;
t206 = t344 * qJD(2);
t186 = t224 * t270;
t184 = t224 * t267;
t109 = -t267 * t316 + t398;
t108 = -qJD(4) * t183 + (t270 * t341 + t258) * qJD(1);
t107 = -t270 * t316 + (-t267 * t344 + t470) * qJD(1);
t106 = -t391 * t475 + (-t265 * t391 - t268 * t396) * rSges(5,1) + t405;
t91 = -t270 * t318 + t442;
t90 = -t270 * t319 + t444;
t83 = t91 * qJD(1);
t82 = t90 * qJD(1);
t54 = -t206 * t394 + (t107 - t383) * qJD(1);
t53 = -t206 * t393 + (-t109 + t384) * qJD(1);
t52 = t338 * t270 + (-t108 + t382) * qJD(1);
t51 = t338 * t267 + (t106 - t381) * qJD(1);
t45 = t277 * t267 - t270 * t496;
t44 = t278 * t267 - t270 * t497;
t43 = t267 * t496 + t277 * t270;
t42 = t267 * t497 + t278 * t270;
t41 = -qJD(2) * t321 + t100 * t269 + t104 * t266;
t40 = -t322 * qJD(2) + t101 * t269 + t105 * t266;
t39 = -qJD(4) * t323 + t102 * t265 + t268 * t98;
t38 = -t324 * qJD(4) + t103 * t265 + t268 * t99;
t29 = t83 + t311;
t28 = t82 + t312;
t27 = t313 + t423;
t26 = t314 + t424;
t1 = [m(3) * (-t107 * t110 - t109 * t111 - t160 * t53 + t162 * t54) + t68 * t492 + (t65 + (t47 + (-t133 * t270 + t134 * t267) * t249 + t360 + t422) * t210 + (t131 * t440 - t428 + t46 + (-t133 * t267 - t134 * t270) * t249 + t421) * t209) * t487 + m(5) * (t52 * (t267 * t380 + t404) + t51 * ((pkin(1) - t476) * t270 + t403) + t93 * t405 + (t223 * t471 - t473) * qJD(4)) + (t83 + ((t58 - t122 + (t150 + t446) * t270 + t418) * t270 + t417 * t267) * qJD(2)) * t373 + (t82 + ((t56 - t121 + (t148 + t448) * t270 + t420) * t270 + t419 * t267) * qJD(4)) * t369 - t195 * t66 / 0.2e1 + (t69 + t67) * t491 + (t425 + (t49 - t309 - t421) * t210 + (t357 * t267 - t112 + t48) * t209 + ((t130 + t326) * t209 + t357 * t210) * t270 + t23) * t490 + (t35 + t32) * t489 + (t27 - t423 + ((t270 * t353 - t417 + t62) * t270 + (t267 * t353 + t358 + t61) * t267) * qJD(2)) * t376 + (t41 + t43) * t394 / 0.2e1 + (t26 - t424 + ((t270 * t354 - t419 + t60) * t270 + (t267 * t354 + t359 + t59) * t267) * qJD(4)) * t372 + (t39 + t42) * t392 / 0.2e1 + (-qJD(2) * t318 + t200 * t269 + t202 * t266 - qJD(4) * t319 + t199 * t268 + t201 * t265 + t248 * t355 - t249 * t356 + m(5) * ((-pkin(1) - t341) * t469 + (-t92 * rSges(5,3) + t380 * t93) * t267)) * qJD(1) + m(4) * (-t36 * t233 - t64 * t236 - t63 * t350 + (t36 * rSges(4,3) + t37 * (-t343 + t484) + t63 * t385 + t64 * t343 * qJD(1)) * t270 + (t36 * (t480 - t484) + t64 * (t385 - t510) + t37 * rSges(4,3) + (-t64 * rSges(4,3) - t63 * (-t478 - t484)) * qJD(1)) * t267) + (t34 + t33 + t24) * t488 - (t40 + t45 + t29) * t393 / 0.2e1 - (t38 + t44 + t28) * t391 / 0.2e1 + ((t86 + t89) * t267 + (t87 + t91) * t270) * t389 / 0.2e1 + ((t84 + t88) * t267 + (t85 + t90) * t270) * qJD(4) * t482; ((t266 * t407 + t269 * t406) * qJD(1) + (t299 * t266 + t269 * t495) * qJD(2)) * t483 + ((-t393 * t442 - t390) * t270 + (t295 + (t270 * t441 + t275) * qJD(2)) * t267) * t373 + ((-t394 * t441 + t390) * t267 + (t295 + (t267 * t442 + t275) * qJD(2)) * t270) * t376 + t298 + (t267 * t41 - t270 * t40 + (t86 * t267 + t270 * t87) * qJD(1)) * t482 + (qJD(1) * t43 + (-(t267 * t503 + t280 * t270) * t270 + (t267 * t504 + t279 * t270) * t267 + (t61 * t267 + t62 * t270) * qJD(1)) * t509) * t486 + (qJD(1) * t45 + (-(t280 * t267 - t270 * t503) * t270 + (t279 * t267 - t270 * t504) * t267 + (t57 * t267 + t58 * t270) * qJD(1)) * t509) * t485 + (t313 + t27) * t378 + (t311 + t29) * t377 + (-t146 * t472 + (t362 * t507 + t25 * t506 + (-t63 * t395 - t468 + (qJD(1) * t63 - t36) * t270 - t402 * t507) * t266) * pkin(2) + (t146 * t270 - t352) * t64 + t505) * m(4) + (-(t110 * t186 + t111 * t184) * qJD(1) - (t81 * (-t184 * t267 - t186 * t270) + t327 * t344) * qJD(2) + 0.2e1 * t81 * (t107 * t270 + t109 * t267 + (t160 * t270 - t162 * t267) * qJD(1)) + t327 * t206 + (-t54 * t267 - t53 * t270 + (t110 * t270 + t111 * t267) * qJD(1)) * t224) * m(3); t298 + (-(-t270 * t64 + t472) * t146 - t64 * t352 + t505) * m(4); (t267 * t39 - t270 * t38 + (t84 * t267 + t270 * t85) * qJD(1)) * t482 + ((-t392 * t443 + t397) * t267 + (t296 + (t267 * t444 + t276) * qJD(4)) * t270) * t372 + ((-t391 * t444 - t397) * t270 + (t296 + (t270 * t443 + t276) * qJD(4)) * t267) * t369 + ((t265 * t409 + t268 * t408) * qJD(1) + (t300 * t265 + t268 * t494) * qJD(4)) * t483 + (qJD(1) * t42 + (-(t267 * t501 + t282 * t270) * t270 + (t267 * t502 + t281 * t270) * t267 + (t59 * t267 + t60 * t270) * qJD(1)) * t508) * t486 + (qJD(1) * t44 + (-(t282 * t267 - t270 * t501) * t270 + (t281 * t267 - t270 * t502) * t267 + (t55 * t267 + t56 * t270) * qJD(1)) * t508) * t485 + (t26 + t314) * t378 + (t28 + t312) * t377 + (0.2e1 * t80 * (t106 * t270 + t108 * t267 + (t159 * t270 - t161 * t267) * qJD(1)) + t339 * t205 + (-t51 * t267 - t52 * t270 + (-t270 * t93 + t471) * qJD(1)) * t223 - (t183 * t92 - t473) * qJD(1) - (t80 * (-t183 * t267 - t185 * t270) + t339 * t341) * qJD(4)) * m(5); 0;];
tauc = t1(:);
