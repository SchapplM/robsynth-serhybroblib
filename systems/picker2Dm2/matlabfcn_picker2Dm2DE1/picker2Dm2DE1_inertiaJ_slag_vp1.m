% Calculate joint inertia matrix for
% picker2Dm2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 18:54
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm2DE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE1_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE1_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2DE1_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2DE1_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 14:31:00
% EndTime: 2020-05-09 14:32:10
% DurationCPUTime: 58.65s
% Computational Cost: add. (1078703->795), mult. (3282474->1302), div. (17297->22), fcn. (553974->31), ass. (0->585)
t363 = sin(pkin(9));
t366 = cos(pkin(9));
t401 = (pkin(3) ^ 2);
t379 = 2 * t401;
t406 = (pkin(1) ^ 2);
t385 = 3 * t406;
t396 = (pkin(4) ^ 2);
t408 = (pkin(7) ^ 2);
t582 = (t408 - t396);
t287 = t379 + t385 + t582;
t361 = sin(qJ(2));
t684 = 0.2e1 * t361;
t580 = pkin(7) * t684;
t296 = pkin(3) * t580;
t362 = sin(qJ(1));
t364 = cos(qJ(2));
t619 = t362 * t364;
t290 = pkin(3) * t619;
t478 = pkin(1) * t290;
t462 = -0.4e1 * t478;
t216 = t287 + t296 + t462;
t314 = pkin(3) * t361;
t294 = t314 + pkin(7);
t365 = cos(qJ(1));
t274 = t294 * t365;
t275 = t290 - pkin(1);
t303 = t406 + t582;
t326 = t361 ^ 2;
t586 = -t401 + t408;
t625 = t326 * t401;
t242 = t296 + t586 + 0.2e1 * t625;
t330 = t365 ^ 2;
t635 = t242 * t330;
t617 = t364 * t365;
t539 = pkin(3) * t617;
t633 = t294 * t362;
t225 = t539 + t633;
t253 = t296 + t303;
t277 = -0.2e1 * t478;
t219 = t277 + t253;
t705 = 4 * t406;
t331 = t401 * t705;
t610 = t401 * t408;
t301 = t331 - 4 * t610;
t393 = -2 * t408;
t394 = 0.2e1 * pkin(3);
t404 = t406 ^ 2;
t527 = -0.4e1 * pkin(3) * pkin(7) * t303;
t708 = 0.4e1 * pkin(1);
t156 = t301 * t326 + t361 * t527 - t404 - (t408 - (t394 + pkin(4)) * pkin(4)) * (t408 + (t394 - pkin(4)) * pkin(4)) + (t393 + (2 * t396) - (4 * t401) - 0.4e1 * t635) * t406 + (-t219 * t274 + t253 * t290) * t708;
t409 = sqrt(t156);
t637 = t225 * t409;
t143 = -t216 * t274 + t637 + (t275 * t580 + t287 * t619) * pkin(3) + (-0.2e1 * t635 + (0.2e1 * t326 - 0.4e1) * t401 - t303) * pkin(1);
t251 = t379 + t253;
t317 = pkin(1) * t365;
t546 = t242 * t317;
t181 = t251 * t294 + 0.2e1 * t546;
t183 = t251 * t365 + (0.4e1 * t330 - 0.2e1) * t294 * pkin(1);
t599 = -t290 + t274;
t223 = pkin(1) + t599;
t315 = pkin(3) * t364;
t144 = t181 * t362 + t183 * t315 + t223 * t409;
t359 = cos(pkin(8));
t653 = sin(pkin(8));
t239 = t359 * t365 + t362 * t653;
t664 = t239 / 0.2e1;
t237 = t359 * t362 - t365 * t653;
t665 = t237 / 0.2e1;
t441 = t143 * t665 + t144 * t664;
t321 = t406 + t408;
t507 = t401 + t321;
t544 = pkin(1) * t274;
t197 = t277 + t296 + t507 + 0.2e1 * t544;
t193 = 0.1e1 / t197;
t402 = 0.1e1 / pkin(3);
t639 = t193 * t402;
t128 = t441 * t639;
t640 = t144 * t237;
t140 = -t639 * t640 / 0.2e1;
t495 = t639 / 0.2e1;
t460 = t143 * t495;
t129 = t239 * t460 + t140;
t118 = t128 * t363 + t129 * t366;
t117 = -t128 * t366 + t129 * t363;
t689 = 0.1e1 / t118 ^ 2;
t694 = 0.1e1 / (t117 ^ 2 * t689 + 0.1e1);
t714 = 0.1e1 / t118 * t694;
t618 = t362 * t365;
t634 = t294 * t330;
t704 = 8 * t406;
t160 = (t242 * t618 + t315 * t634) * t704 + (t219 * t633 + t253 * t539) * t708;
t155 = 0.1e1 / t409;
t669 = t155 / 0.2e1;
t496 = t225 * t669;
t674 = pkin(7) * t365;
t579 = 0.2e1 * t674;
t620 = t361 * t401;
t136 = t599 * t409 + t160 * t496 + (t216 * t294 + 0.4e1 * t546) * t362 + (t579 * t620 + (t287 * t365 + t634 * t708) * pkin(3)) * t364;
t576 = 0.2e1 * pkin(1);
t218 = t225 * t576;
t194 = 0.1e1 / t197 ^ 2;
t434 = t194 * t441;
t327 = t362 ^ 2;
t497 = t223 * t669;
t692 = -0.2e1 * pkin(1);
t134 = -t637 + t160 * t497 + t242 * t327 * t692 + t181 * t365 + (-t251 - 0.8e1 * t544) * t290;
t703 = t143 / 0.2e1;
t489 = t703 + t134 / 0.2e1;
t97 = t140 + (t218 * t434 + (t136 * t665 + t239 * t489) * t193) * t402;
t433 = t194 * (-t143 * t239 / 0.2e1 + t640 / 0.2e1);
t98 = (-t218 * t433 + ((t136 / 0.2e1 - t144 / 0.2e1) * t239 - t489 * t237) * t193) * t402;
t716 = (t363 * t98 - t366 * t97) * t714;
t566 = 0.2e1 * t315;
t297 = pkin(7) * t566;
t616 = t364 * t401;
t518 = t361 * t616;
t248 = t297 + 0.4e1 * t518;
t608 = t406 * t330;
t551 = -0.4e1 * t608;
t675 = pkin(7) * t362;
t561 = t401 * t675;
t622 = t361 * t362;
t541 = pkin(3) * t622;
t479 = pkin(1) * t541;
t600 = 0.2e1 * t479 + t297;
t152 = t248 * t551 + (t301 * t684 + t527) * t364 + (-t600 * t274 + 0.2e1 * t364 ^ 2 * t561 + (-t219 * t617 - t253 * t622) * pkin(3)) * t708;
t572 = 0.2e1 * t317;
t660 = pkin(1) * t330;
t693 = 0.2e1 * pkin(7);
t132 = t152 * t497 + t248 * t362 * t572 + ((t362 * t409 - t183) * t361 + (t365 * t409 + (t294 * t693 + t251) * t362 + (-pkin(1) + 0.2e1 * t660 + t674) * t566) * t364) * pkin(3);
t621 = t361 * t365;
t540 = pkin(3) * t621;
t133 = (t290 - t540) * t409 + t152 * t496 - 0.2e1 * t248 * t660 - (t297 + 0.4e1 * t479) * t274 - 0.2e1 * t326 * t561 - t287 * t541 + (t620 * t708 + (-t216 * t365 + t275 * t693) * pkin(3)) * t364;
t480 = pkin(1) * t539;
t214 = 0.2e1 * t480 + t600;
t102 = (-t214 * t434 + (t132 * t664 + t133 * t665) * t193) * t402;
t103 = (t214 * t433 + (t133 * t664 - t132 * t237 / 0.2e1) * t193) * t402;
t715 = t714 * (-t102 * t366 + t103 * t363);
t567 = 0.8e1 * t314;
t325 = t361 * t326;
t413 = pkin(3) * t401;
t627 = t325 * t413;
t535 = 0.32e2 * t627;
t245 = t617 + t622;
t246 = t619 - t621;
t271 = t408 + t401 / 0.4e1 + t406 / 0.4e1 - t396 / 0.8e1;
t395 = t396 ^ 2;
t407 = t408 ^ 2;
t414 = t401 ^ 2;
t595 = 0.4e1 / 0.7e1 * t408 - t396 / 0.7e1;
t607 = t408 * t396;
t174 = -0.32e2 / 0.21e2 * t271 * t478 + 0.5e1 / 0.42e2 * t414 + (0.16e2 / 0.21e2 * t406 + t595) * t401 + t404 / 0.7e1 + t595 * t406 + t407 - 0.3e1 / 0.7e1 * t607 + t395 / 0.42e2;
t597 = t406 / 0.3e1 + t408;
t342 = -t396 / 0.4e1;
t697 = t342 + t401 / 0.2e1;
t273 = t597 + t697;
t343 = -t396 / 0.3e1;
t345 = -0.2e1 / 0.3e1 * t396;
t349 = 0.4e1 / 0.3e1 * t401;
t355 = 0.4e1 / 0.3e1 * t406;
t175 = -0.8e1 / 0.3e1 * t273 * t478 + 0.5e1 / 0.18e2 * t414 + (t355 + t343) * t401 + t407 - t404 / 0.3e1 + t395 / 0.18e2 + (t349 + 0.2e1 / 0.3e1 * t406 + t345) * t408;
t585 = t404 + t407;
t392 = 2 * t408;
t589 = t392 - t396;
t439 = (t589 * t406) + t395 / 0.6e1 + t585 - t607;
t228 = -t414 / 0.6e1 + t439;
t357 = t406 / 0.2e1;
t596 = t357 + t408;
t230 = -0.2e1 / 0.3e1 * t478 + t342 + t596;
t390 = 4 * t408;
t307 = (t390 + t396) * t406;
t354 = -0.2e1 / 0.3e1 * t401;
t311 = t354 + t408;
t312 = -t406 / 0.3e1 + t408;
t322 = -t406 + t408;
t344 = -t396 / 0.2e1;
t284 = t344 + t507;
t483 = -0.4e1 * t290;
t452 = t284 * t483;
t329 = t365 * t330;
t410 = pkin(1) * t406;
t623 = t329 * t410;
t559 = pkin(7) * t623;
t484 = 0.16e2 * t559;
t328 = t330 ^ 2;
t609 = t404 * t328;
t552 = 0.8e1 * t609;
t613 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t146 = t311 * t552 + t230 * t484 + 0.14e2 * t174 * t608 + (t322 * t414) + (t307 - 0.10e2 / 0.3e1 * t404 + (2 * t407) - t607) * t401 + t228 * t613 + (0.6e1 * t175 * t674 + t312 * t452) * pkin(1);
t346 = -0.3e1 / 0.2e1 * t396;
t367 = 15 * t404;
t374 = 18 * t408;
t387 = 3 * t407;
t598 = t395 / 0.2e1 - t414 / 0.2e1;
t461 = -(3 * t607) + t387 + t598;
t391 = 3 * t408;
t594 = 15 * t406 + t391;
t398 = t413 ^ 2;
t603 = t321 * ((t346 + t392) * t406 - 0.3e1 / 0.2e1 * t607 + t585 + t598) + t398;
t332 = 0.10e2 / 0.3e1 * t406;
t435 = 0.5e1 / 0.6e1 * t414 + t439;
t207 = (t332 + t589) * t401 + t435;
t706 = -0.6e1 * t207;
t159 = t478 * t706 + (t367 + ((t374 - 9 * t396) * t406) + t461) * t401 + (t346 + t594) * t414 + t603;
t353 = -t401 / 0.3e1;
t310 = t353 + t408;
t448 = pkin(1) * t452;
t584 = t407 - t404;
t168 = t310 * t448 - t398 + (-t332 - t582) * t414 + (t307 + t414 / 0.6e1 - t395 / 0.6e1 + t584) * t401 + t228 * t408;
t376 = -5 * t396;
t382 = 7 * t404;
t172 = (t346 + t391 + (7 * t406)) * t414 + (t382 + ((t376 + 10 * t408) * t406) + t461) * t401 + t603;
t418 = pkin(7) * t408;
t302 = -0.12e2 * pkin(7) * t410 + t418 * t708;
t606 = t408 * t406;
t309 = -8 * t404 + 12 * t606;
t184 = t302 * t365 + t309 * t330 + t484 + t552 + t585 - (6 * t606);
t320 = -3 * t401 + t408;
t614 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t203 = t277 * t614 + t284 * t320;
t583 = t407 + t414;
t264 = 16 * (t583 - 6 * t610) * t404;
t288 = -t396 + t507;
t295 = t317 + pkin(7);
t318 = -30 * t396 + 60 * t408;
t324 = t326 ^ 2;
t375 = -2 * t396;
t377 = -6 * t396;
t384 = 6 * t406;
t588 = t395 - t414;
t444 = 6 * t407 + t588 - 6 * t607;
t313 = t321 ^ 2;
t630 = t313 * (-t401 + t303);
t180 = t448 + ((t384 + t589) * t401) + t435;
t244 = t310 * t277;
t201 = t284 * t614 + t244;
t490 = 0.8e1 * t559;
t263 = t320 * t490;
t574 = pkin(7) * t317;
t525 = 0.6e1 * t574;
t529 = 0.12e2 * t608;
t147 = t180 * t525 + t201 * t529 + t159 + t263;
t388 = 8 * t408;
t220 = t410 * t483 + t331 + (4 * t404) + ((t375 + t388) * t406);
t447 = t408 - t478;
t226 = -t406 + t447 + t697;
t323 = -3 * t406 + t408;
t482 = -0.2e1 * t290;
t578 = 0.4e1 * t674;
t162 = t490 + t220 * t330 + t284 * t323 + (t226 * t578 + t482 * t613) * pkin(1);
t695 = t147 * t567 + t162 * t535;
t122 = t264 * t328 + 0.32e2 * t203 * t559 + 0.24e2 * t168 * t608 + (t375 + t390 + 28 * t406) * t398 + (t288 * t630) + (0.24e2 * t146 * t326 + (t318 * t404) + (t377 * t407) + (t384 * t444) + (t392 * t588) + 0.28e2 * t410 ^ 2 + 0.4e1 * t418 ^ 2) * t401 + t695 * t295 + 0.8e1 * (t159 * t674 - t172 * t290) * pkin(1) + (0.16e2 * t184 * t324 + (t318 * t406) + (70 * t404) + t414 + t444) * t414;
t341 = -t396 / 0.6e1;
t269 = 0.7e1 / 0.6e1 * t401 + t341 + t596;
t351 = t401 / 0.3e1;
t512 = t341 + t351 + t408;
t272 = t355 + t512;
t316 = pkin(1) * t362;
t198 = -t269 * t316 + t272 * t315;
t281 = t406 + t512;
t306 = t379 + t322;
t209 = t281 * t315 - t306 * t316 / 0.2e1;
t513 = t396 / 0.3e1 + t351 + t392;
t459 = -0.8e1 / 0.3e1 * t609 + (t401 * t322) - 0.5e1 / 0.3e1 * t404 + t513 * t406 + t408 * (t343 + t310);
t549 = 0.4e1 * t608;
t673 = pkin(7) * t404;
t577 = -0.4e1 * t673;
t333 = -0.20e2 / 0.3e1 * t406;
t350 = 0.2e1 / 0.3e1 * t401;
t514 = 0.2e1 / 0.3e1 * t396 + t350 + t390;
t515 = 0.4e1 / 0.3e1 * t396 + t349 + t393;
t666 = t414 / 0.2e1 - (t333 + t514) * t401 / 0.2e1 + 0.3e1 / 0.2e1 * t404 - t515 * t406 / 0.2e1 - t407 / 0.2e1;
t150 = t362 * t329 * t577 + t198 * t549 + t459 * t315 + (t209 * t578 + t362 * t666) * pkin(1);
t383 = 5 * t404;
t389 = 6 * t408;
t592 = t375 - 2 * t401;
t508 = t389 + t592;
t369 = 10 * t406;
t593 = t369 + t392;
t190 = t414 + (t345 + t354 + t593) * t401 + t383 + (t508 * t406) + t408 * (t345 + t311);
t182 = t190 * t315;
t511 = t343 + t321;
t278 = 0.8e1 / 0.3e1 * t401 + t511;
t590 = t385 + t408;
t280 = t343 + t350 + t590;
t202 = -t278 * t316 + t280 * t315;
t285 = 0.5e1 / 0.6e1 * t401 + t357 + t341;
t213 = t285 * t566 + t316 * t614;
t537 = t410 * t315;
t466 = t329 * t537;
t378 = 5 * t414;
t509 = t345 + t321;
t208 = t378 + ((t369 + t508) * t401) + (t354 + t509) * t321;
t638 = t208 * t362;
t151 = -0.8e1 * pkin(7) * t466 + t213 * t551 + t182 + (t202 * t578 - t638) * pkin(1);
t169 = -pkin(1) * t638 + t182;
t206 = -(3 * t414) + (t333 + t515) * t401 + t514 * t406 + t584;
t210 = -0.5e1 / 0.3e1 * t414 + (-t406 + t513) * t401 + t408 * (t353 + t511);
t573 = -0.2e1 * t316;
t170 = t206 * t315 + t210 * t573;
t538 = t406 * t315;
t443 = -t362 * t410 + t538;
t255 = 0.4e1 * t443;
t265 = 0.2e1 * t316 + t315;
t304 = t401 + t590;
t283 = t344 + t304;
t173 = t323 * t315 + t255 * t330 + (t265 * t674 + t283 * t362) * t576;
t270 = t408 + 0.5e1 / 0.2e1 * t401 + 0.3e1 / 0.2e1 * t406 + t344;
t629 = t320 * t362;
t211 = t270 * t315 + pkin(1) * t629 / 0.2e1;
t696 = t391 - t396 - t401;
t286 = t696 * t369;
t591 = t376 - 5 * t401;
t451 = 0.24e2 * t310 * t609 - t398 - ((21 * t406 + t696) * t414) - ((t408 * t592 + t286 + t387 + 35 * t404) * t401) - ((t382 + (t388 + t591) * t406 + t408 * (-t401 + t582)) * t321);
t521 = t295 * t627;
t469 = -0.8e1 * t521;
t534 = -0.12e2 * t625;
t632 = t295 * t361;
t299 = pkin(7) * t572;
t232 = 0.4e1 / 0.3e1 * t608 + t299 + t312;
t628 = t324 * t414;
t698 = 0.7e1 * t398 + ((35 * t406 + 15 * t408 + t591) * t414) + ((21 * t404 + t286 + 9 * t407 + (t377 - 6 * t401) * t408) * t401) + t630 - 0.24e2 * t232 * t628;
t127 = t211 * t484 + t173 * t469 + t150 * t534 - 0.6e1 * t170 * t608 + (-0.6e1 * t151 * t632 + t364 * t451) * pkin(3) + (-0.6e1 * t169 * t674 + t362 * t698) * pkin(1);
t256 = t314 + t295;
t107 = t122 * t256 + t127 * t409;
t222 = -0.4e1 / 0.9e1 * t478 + 0.4e1 / 0.9e1 * t401 - t396 / 0.9e1 + t597;
t231 = t341 + t350 + t447;
t279 = t349 + t511;
t548 = 0.6e1 * t608;
t161 = 0.4e1 * t559 + t222 * t548 + t279 * t613 + (t231 * t578 + t312 * t482) * pkin(1);
t510 = t345 + t350 + t392;
t601 = (t350 + t509) * t321 + t414;
t177 = t279 * t462 + (t384 + t510) * t401 + t601;
t199 = t279 * t614 + t244;
t200 = (t332 + t510) * t401 + t601;
t227 = t277 + t279;
t300 = t586 * t705;
t571 = 0.4e1 * t317;
t163 = pkin(7) * t227 * t571 + t300 * t330 + t177;
t243 = t299 + t549 + t323;
t557 = 0.8e1 * t627;
t568 = 0.6e1 * t314;
t438 = t163 * t568 + t243 * t557;
t533 = 0.12e2 * t625;
t707 = 0.6e1 * pkin(1);
t137 = t161 * t533 + t263 + t199 * t529 + t398 + ((-t396 + t401 + t594) * t414) + ((t367 + (t374 + t377 + 6 * t401) * t406 + t387 + (t375 + t379) * t408) * t401) + (t313 * t288) + t438 * t295 + (t177 * t674 - t200 * t290) * t707;
t254 = 0.2e1 * t443;
t319 = t401 + 2 * t406;
t179 = t322 * t315 + t254 * t330 + (t297 * t365 + t319 * t362) * pkin(1);
t262 = t304 * t315;
t305 = 3 * t401 + t321;
t631 = t305 * t362;
t215 = -pkin(1) * t631 + t262;
t547 = t413 * t704;
t659 = pkin(3) * t404;
t293 = t547 + 0.4e1 * t659;
t517 = t410 * t614;
t217 = t293 * t364 + 0.4e1 * t362 * t517;
t236 = t401 * t593 + t383 + t583 + 6 * t606;
t249 = t378 + (t369 + t389) * t401 + t313;
t556 = -0.4e1 * t625;
t266 = -t316 + t315;
t467 = t330 * t538;
t176 = -0.2e1 * t467 + t262 + (t266 * t579 - t631) * pkin(1);
t687 = -0.4e1 * t176;
t145 = t179 * t556 + t217 * t330 + (t632 * t687 + (-t236 + t490) * t364) * pkin(3) + (-0.4e1 * t215 * t674 + (t249 + t469) * t362) * pkin(1);
t126 = t137 * t256 + t145 * t409;
t124 = 0.1e1 / t126;
t671 = t124 / 0.4e1;
t504 = t144 * t671;
t662 = t409 / 0.4e1;
t437 = t107 * t504 + t143 * t662;
t397 = 0.1e1 / pkin(4);
t611 = t397 / pkin(3) ^ 2;
t522 = t193 * t611;
t84 = t437 * t522;
t82 = t246 * t84;
t505 = t143 * t671;
t663 = -t409 / 0.4e1;
t436 = t107 * t505 + t144 * t663;
t83 = t436 * t522;
t61 = -t245 * t83 + t82;
t81 = t245 * t84;
t62 = -t246 * t83 - t81;
t48 = atan2(t61, t62);
t45 = sin(t48);
t46 = cos(t48);
t429 = atan2(t144 * t495, t460);
t427 = sin(t429);
t428 = cos(t429);
t120 = -t362 * t428 - t365 * t427;
t121 = t362 * t427 - t365 * t428;
t612 = t397 * t402;
t492 = t612 / 0.2e1;
t87 = atan2(t409 * t492, t107 * t124 * t492);
t85 = sin(t87);
t86 = cos(t87);
t481 = -t120 * t85 + t121 * t86;
t51 = t120 * t86 + t121 * t85;
t713 = t45 * t51 - t46 * t481;
t712 = t45 * t481 + t46 * t51;
t95 = atan2(t117, t118);
t89 = sin(t95);
t91 = cos(t95);
t450 = t120 * t89 - t121 * t91;
t57 = t120 * t91 + t121 * t89;
t96 = atan2(t117, -t118);
t90 = sin(t96);
t92 = cos(t96);
t711 = -t450 * t92 - t57 * t90;
t710 = -t450 * t90 + t57 * t92;
t709 = t126 ^ 2;
t686 = 0.8e1 * t256;
t683 = -0.2e1 * t365;
t682 = -0.6e1 * t409;
t681 = pkin(1) * pkin(7);
t680 = pkin(4) * t51;
t679 = pkin(4) * t481;
t678 = pkin(6) * t57;
t677 = pkin(6) * t450;
t676 = pkin(7) * t330;
t672 = t124 / 0.2e1;
t125 = 0.1e1 / t709;
t670 = -t125 / 0.4e1;
t668 = t193 / 0.2e1;
t667 = -t194 / 0.2e1;
t661 = pkin(1) * t312;
t658 = pkin(3) * t409;
t657 = t172 * pkin(1);
t60 = 0.1e1 / t62 ^ 2;
t656 = t60 * t61;
t453 = pkin(7) * t467;
t268 = -0.4e1 * t453;
t298 = pkin(7) * t573;
t454 = t537 * t676;
t463 = t329 * t517;
t491 = pkin(7) * t549;
t498 = t160 * t669;
t615 = t365 * t406;
t516 = t362 * t615;
t531 = -0.24e2 * t615;
t555 = 0.8e1 * t625;
t558 = -0.8e1 * t627;
t605 = t410 * t330;
t562 = pkin(7) * t605;
t569 = -0.4e1 * t314;
t485 = -0.24e2 * t562;
t602 = -0.24e2 * t310 * t466 + t485 * t629;
t108 = (t327 * t325 * t547 + (t317 * t319 - 0.2e1 * t623) * t556 + 0.4e1 * t463 + t305 * t491 + t249 * t317) * t409 + t145 * t498 + ((-0.8e1 / 0.3e1 * t466 + t268 - 0.2e1 * t312 * t480) * t533 - 0.24e2 * t279 * t453 - 0.6e1 * t200 * t480 + t602) * t256 + ((t254 * t365 * t555 + t217 * t683 - 0.24e2 * t454) * t409 + (0.12e2 * (-t222 * t615 - t562) * t533 + t199 * t531) * t256 + (0.4e1 * t176 * t361 * t658 - t438 * t256 - t137 + ((t315 * t555 + 0.4e1 * t215) * t409 + (-0.48e2 * t231 * t625 - 0.6e1 * t177) * t256) * pkin(7)) * pkin(1)) * t362 + ((t558 * t317 + ((0.4e1 * t290 * t365 - 0.2e1 * t676) * t406 + (-0.2e1 * t266 * t675 - t305 * t365) * pkin(1)) * t569) * t409 + ((t298 - 0.8e1 * t516) * t557 + (-0.8e1 * t453 - 0.2e1 * t300 * t618 + (-t227 * t675 - t279 * t539) * t708) * t568) * t256) * t295;
t431 = t194 * t437;
t501 = t143 * t155 / 0.8e1;
t502 = t144 * t670;
t488 = 0.32e2 / 0.3e1 * t404;
t455 = t329 * t488;
t456 = 0.64e2 / 0.3e1 * t271 * t410;
t464 = 0.32e2 * t521;
t471 = t614 * t673;
t486 = -0.48e2 * t562;
t530 = -0.32e2 * t404 * t329;
t532 = 0.24e2 * t625;
t536 = -0.96e2 * t310 * t329;
t543 = pkin(1) * t613;
t545 = t284 * t661;
t550 = -0.2e1 * t608;
t553 = 0.8e1 * t615;
t554 = -0.4e1 * t623;
t563 = pkin(7) * t608;
t570 = -0.6e1 * t314;
t644 = -0.4e1 * t681;
t645 = -0.6e1 * t681;
t73 = ((t283 * t572 + t491 + t554) * t469 + (0.12e2 * t327 * t330 * t673 + t269 * t554 + t328 * t577) * t534 + 0.12e2 * t210 * t623 + (t534 * t666 + t698) * t317 + (t306 * t534 * t550 + t208 * t548 + t320 * t552) * pkin(7)) * t409 + t127 * t498 + t268 * t256 * t464 + (((pkin(7) * t278 * t551 - t208 * t317 - 0.4e1 * t463) * t682 + t602 * t686) * t632 + ((-pkin(7) * t328 * t488 - 0.16e2 * t273 * t563 - t329 * t456 - 0.4e1 * t365 * t545) * t532 - 0.64e2 * t328 * t471 + t410 * t284 * t536 - 0.48e2 * t207 * t563 - 0.8e1 * t365 * t657 + ((t543 * t683 + t554) * t535 + (-0.24e2 * t284 * t563 + t317 * t706) * t567) * t295) * t256 * t364) * pkin(3) + ((0.16e2 * (t309 * t683 - t302 + t486 + t530) * t628 + (t220 * t683 + t226 * t644 + t485) * t464 + (-0.28e2 * t174 * t615 + t175 * t645 + t230 * t486 + t311 * t530) * t532 + t295 * (t180 * t645 + t201 * t531) * t567 - 0.4e1 * t264 * t329 - 0.96e2 * t203 * t562 - 0.48e2 * t168 * t615 - 0.8e1 * t159 * t681) * t256 + ((-0.8e1 * t198 * t615 + t315 * t455) * t534 + t364 * t536 * t659 + t211 * t486 + 0.12e2 * t170 * t615 + (0.2e1 * (-t255 * t365 - t265 * t681) * t558 + (t202 * t644 + t213 * t553 + 0.24e2 * t454) * t570) * t295) * t409 + (-t695 * t256 - t122 + (t151 * t568 + t173 * t557 + (-0.24e2 * t298 + 0.64e2 * t516) * t628 + (0.48e2 * t209 * t625 + 0.6e1 * t169) * pkin(7)) * t409) * pkin(1)) * t362;
t655 = (t218 * t431 + (t136 * t662 + t160 * t501 + t73 * t504 + (t108 * t502 + t134 * t671) * t107) * t193) * t611 + t83;
t519 = t310 * t605;
t229 = 0.24e2 * t519 * t541;
t487 = pkin(7) * t540;
t446 = t406 * t362 * t487;
t252 = 0.4e1 * t446;
t267 = t487 * t692;
t499 = t152 * t669;
t624 = t326 * t413;
t520 = t295 * t624;
t109 = (t267 * t556 + (-0.8e1 * t179 * t616 - t293 * t330 + ((-t322 + t550) * t556 + t236 + (t304 * t571 - 0.8e1 * t623) * pkin(7)) * pkin(3)) * t361 + ((t267 + (-t304 + 0.2e1 * t608) * t314) * t569 + (pkin(3) * t687 - 0.24e2 * t316 * t624) * t364) * t295) * t409 + t145 * t499 + (0.24e2 * t243 * t364 * t520 + (t252 + (0.8e1 / 0.3e1 * t605 + 0.2e1 * t661) * t541) * t533 + 0.24e2 * t161 * t518 + t229 + 0.24e2 * t279 * t446 + 0.6e1 * t200 * t479 + 0.6e1 * ((pkin(7) * t553 + t279 * t708) * t362 * t625 + t163 * t315) * t295) * t256 + t137 * t315;
t526 = -0.4e1 * t574;
t560 = pkin(7) * t615;
t626 = t325 * t414;
t74 = t127 * t499 + (0.32e2 * t252 * t256 - 0.8e1 * t267 * t409) * t521 + (0.24e2 * (-0.4e1 * t232 * t316 * t626 - t150 * t620 - t173 * t520) * t409 + (0.48e2 * t146 * t620 + 0.96e2 * t162 * t520 + 0.64e2 * t184 * t626) * t256) * t364 + ((t122 + (t147 * t686 + t151 * t682) * t295) * t364 + (t295 * t229 * t686 + ((t272 * t551 + t281 * t526 - t459) * t534 - 0.16e2 * t270 * t559 + t206 * t548 + t190 * t525 + ((-t323 + t551) * t558 + (t280 * t526 + 0.8e1 * t285 * t608 - t190 + t490) * t570) * t295 - t451) * t409 + ((pkin(7) * t455 + 0.16e2 * t273 * t560 + t330 * t456 + 0.4e1 * t545) * t532 + 0.64e2 * t329 * t471 + 0.96e2 * t284 * t519 + 0.48e2 * t207 * t560 + 0.8e1 * t657 + ((0.2e1 * t543 + 0.4e1 * t605) * t535 + (t207 * t707 + 0.24e2 * t284 * t560) * t567) * t295) * t256 * t362) * t361) * pkin(3);
t654 = (-t214 * t431 + (t133 * t662 + t152 * t501 + t74 * t504 + (t109 * t502 + t132 * t671) * t107) * t193) * t611 - t83;
t142 = 0.1e1 / t143 ^ 2;
t468 = pkin(3) / (t142 * t144 ^ 2 + 0.1e1) * t197 * t402;
t442 = -0.2e1 * t142 * t144 * t468;
t445 = t468 / t703;
t494 = t214 * t667;
t79 = (t132 * t668 + t144 * t494) * t445 + (t133 * t668 + t143 * t494) * t442;
t649 = t120 * t79;
t493 = t218 * t667;
t80 = 0.1e1 + (t134 * t668 - t144 * t493) * t445 + (t136 * t668 - t143 * t493) * t442;
t648 = t120 * t80;
t647 = t121 * t79;
t646 = t121 * t80;
t643 = pkin(3) * t646 - t317;
t642 = pkin(2) * t646 - t317;
t641 = 0.1e1 / t107 ^ 2 * t709;
t636 = t237 ^ 2 / t239 ^ 2;
t604 = Icges(7,3) + Icges(4,3);
t565 = pkin(2) * t649;
t564 = pkin(3) * t649;
t528 = t117 * t689 * t694;
t524 = t79 + t715;
t506 = -t107 * t125 / 0.2e1;
t503 = t143 * t670;
t500 = -t144 * t155 / 0.8e1;
t472 = t80 + t716;
t88 = 0.1e1 / (t156 * t641 + 0.1e1);
t432 = -0.2e1 * pkin(4) * t88 * t612 * t641 * t658;
t449 = 0.1e1 / t107 * t126 * t88 * t669;
t28 = t152 * t449 + (t109 * t506 + t672 * t74) * t432 + t79;
t458 = -pkin(2) * t648 + t316;
t457 = -pkin(3) * t648 + t316;
t27 = t160 * t449 + (t108 * t506 + t672 * t73) * t432 + t80;
t430 = t194 * t436;
t261 = -rSges(2,1) * t365 + rSges(2,2) * t362;
t260 = rSges(8,1) * t364 - rSges(8,2) * t361;
t259 = rSges(2,1) * t362 + rSges(2,2) * t365;
t258 = rSges(8,1) * t361 + rSges(8,2) * t364;
t196 = atan2(-t237, -t239);
t195 = atan2(-t237, t239);
t191 = 0.1e1 / (0.1e1 + t636);
t189 = cos(t196);
t188 = cos(t195);
t187 = sin(t196);
t186 = sin(t195);
t167 = -t186 * t362 + t188 * t365;
t166 = t186 * t365 + t188 * t362;
t165 = -t187 * t653 + t359 * t189;
t164 = t359 * t187 + t189 * t653;
t153 = -t191 * t636 - t191 + 0.1e1;
t149 = rSges(6,1) * t165 - rSges(6,2) * t164;
t148 = rSges(6,1) * t164 + rSges(6,2) * t165;
t139 = -t317 + t153 * (rSges(9,1) * t167 - rSges(9,2) * t166);
t138 = t316 - t153 * (rSges(9,1) * t166 + rSges(9,2) * t167);
t100 = rSges(3,1) * t121 - rSges(3,2) * t120;
t99 = rSges(3,1) * t120 + rSges(3,2) * t121;
t68 = pkin(2) * t647;
t67 = pkin(3) * t647;
t64 = t100 * t80 - t317;
t63 = -t80 * t99 + t316;
t59 = 0.1e1 / t62;
t47 = 0.1e1 / (t60 * t61 ^ 2 + 0.1e1);
t44 = rSges(4,1) * t450 + rSges(4,2) * t57;
t43 = rSges(7,1) * t450 + rSges(7,2) * t57;
t42 = -rSges(4,1) * t57 + rSges(4,2) * t450;
t41 = -rSges(7,1) * t57 + rSges(7,2) * t450;
t38 = rSges(5,1) * t481 - rSges(5,2) * t51;
t37 = rSges(5,1) * t51 + rSges(5,2) * t481;
t35 = (-t214 * t430 + (t74 * t505 + t132 * t663 + t152 * t500 + (t109 * t503 + t133 * t671) * t107) * t193) * t611;
t33 = (t218 * t430 + (t73 * t505 + t134 * t663 + t160 * t500 + (t108 * t503 + t136 * t671) * t107) * t193) * t611;
t32 = -t528 * (t102 * t363 + t103 * t366) + t524;
t30 = -t528 * (t363 * t97 + t366 * t98) + t472;
t26 = t524 - t715;
t25 = t30 * t43 - t317;
t24 = -t30 * t41 + t316;
t23 = rSges(10,1) * t711 - rSges(10,2) * t710;
t22 = rSges(10,1) * t710 + rSges(10,2) * t711;
t21 = t32 * t44 + t68;
t20 = -t32 * t42 - t565;
t19 = t472 - t716;
t18 = t30 * t44 + t642;
t17 = -t30 * t42 + t458;
t16 = t28 * t38 + t67;
t15 = -t28 * t37 - t564;
t14 = t27 * t38 + t643;
t13 = -t27 * t37 + t457;
t12 = rSges(11,1) * t713 + rSges(11,2) * t712;
t11 = -rSges(11,1) * t712 + rSges(11,2) * t713;
t10 = t23 * t26 + t32 * t677 + t68;
t9 = -t22 * t26 + t32 * t678 - t565;
t8 = t19 * t23 + t30 * t677 + t642;
t7 = -t19 * t22 + t30 * t678 + t458;
t6 = ((-t245 * t35 + t246 * t654 - t81) * t59 - ((-t35 - t84) * t246 - t654 * t245) * t656) * t47 + t28;
t5 = ((t655 * t246 + (-t33 + t84) * t245) * t59 - (-t245 * t655 - t246 * t33 + t82) * t656) * t47 + t27;
t4 = t12 * t6 + t28 * t679 + t67;
t3 = -t11 * t6 - t28 * t680 - t564;
t2 = t12 * t5 + t27 * t679 + t643;
t1 = -t11 * t5 - t27 * t680 + t457;
t29 = [t80 ^ 2 * Icges(3,3) + t153 ^ 2 * Icges(9,3) + t5 ^ 2 * Icges(11,3) + t19 ^ 2 * Icges(10,3) + t27 ^ 2 * Icges(5,3) + (t7 ^ 2 + t8 ^ 2) * m(10) + (t13 ^ 2 + t14 ^ 2) * m(5) + (t63 ^ 2 + t64 ^ 2) * m(3) + (t17 ^ 2 + t18 ^ 2) * m(4) + (t1 ^ 2 + t2 ^ 2) * m(11) + (t24 ^ 2 + t25 ^ 2) * m(7) + Icges(2,3) + t604 * t30 ^ 2 + Icges(6,3) + m(9) * (t138 ^ 2 + t139 ^ 2) + m(2) * (t259 ^ 2 + t261 ^ 2) + m(6) * (t148 ^ 2 + t149 ^ 2); (t10 * t8 + t7 * t9) * m(10) + (t13 * t15 + t14 * t16) * m(5) + (t17 * t20 + t18 * t21) * m(4) + t5 * Icges(11,3) * t6 + t19 * Icges(10,3) * t26 + t27 * Icges(5,3) * t28 + (t1 * t3 + t2 * t4) * m(11) + (t80 * Icges(3,3) + (t100 * t64 - t63 * t99) * m(3)) * t79 + ((-t24 * t41 + t25 * t43) * m(7) + t604 * t30) * t32; (t3 ^ 2 + t4 ^ 2) * m(11) + (t10 ^ 2 + t9 ^ 2) * m(10) + m(8) * (t258 ^ 2 + t260 ^ 2) + (t15 ^ 2 + t16 ^ 2) * m(5) + t6 ^ 2 * Icges(11,3) + t28 ^ 2 * Icges(5,3) + t26 ^ 2 * Icges(10,3) + (t20 ^ 2 + t21 ^ 2) * m(4) + Icges(8,3) + (Icges(3,3) + (t100 ^ 2 + t99 ^ 2) * m(3)) * t79 ^ 2 + ((t41 ^ 2 + t43 ^ 2) * m(7) + t604) * t32 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t29(1), t29(2); t29(2), t29(3);];
Mq = res;
