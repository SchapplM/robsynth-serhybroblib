% Calculate Gravitation load on the joints for
% picker2Dm1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm1DE2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE2_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1DE2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 20:28:31
% EndTime: 2020-05-10 20:34:34
% DurationCPUTime: 68.18s
% Computational Cost: add. (1294550->820), mult. (3572958->1396), div. (46468->42), fcn. (939892->38), ass. (0->619)
t734 = 4 * pkin(1);
t344 = sin(qJ(2));
t297 = pkin(3) * t344;
t578 = 0.8e1 * t297;
t309 = t344 ^ 2;
t308 = t344 * t309;
t391 = pkin(3) ^ 2;
t408 = pkin(3) * t391;
t638 = t308 * t408;
t541 = 0.32e2 * t638;
t386 = pkin(4) ^ 2;
t345 = sin(qJ(1));
t347 = cos(qJ(2));
t630 = t345 * t347;
t273 = pkin(3) * t630;
t475 = pkin(1) * t273;
t397 = pkin(1) ^ 2;
t400 = pkin(7) ^ 2;
t607 = t397 / 0.3e1 + t400;
t211 = -0.4e1 / 0.9e1 * t475 + 0.4e1 / 0.9e1 * t391 - t386 / 0.9e1 + t607;
t324 = -t386 / 0.6e1;
t333 = 0.2e1 / 0.3e1 * t391;
t445 = t400 - t475;
t223 = t324 + t333 + t445;
t332 = 0.4e1 / 0.3e1 * t391;
t304 = t397 + t400;
t326 = -t386 / 0.3e1;
t509 = t326 + t304;
t262 = t332 + t509;
t295 = -t397 / 0.3e1 + t400;
t479 = -0.2e1 * t273;
t348 = cos(qJ(1));
t313 = t348 ^ 2;
t618 = t397 * t313;
t557 = 0.6e1 * t618;
t312 = t348 * t313;
t404 = pkin(1) * t397;
t634 = t312 * t404;
t568 = pkin(7) * t634;
t701 = pkin(7) * t348;
t589 = 0.4e1 * t701;
t624 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t153 = 0.4e1 * t568 + t211 * t557 + t262 * t624 + (t223 * t589 + t295 * t479) * pkin(1);
t367 = 6 * t397;
t460 = -0.4e1 * t475;
t328 = -0.2e1 / 0.3e1 * t386;
t376 = 0.2e1 * t400;
t508 = t328 + t333 + t376;
t409 = t391 ^ 2;
t507 = t328 + t304;
t611 = (t333 + t507) * t304 + t409;
t167 = t262 * t460 + (t367 + t508) * t391 + t611;
t260 = -0.2e1 * t475;
t336 = -t391 / 0.3e1;
t293 = t336 + t400;
t231 = t293 * t260;
t625 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t181 = t262 * t625 + t231;
t315 = 0.10e2 / 0.3e1 * t397;
t182 = (t315 + t508) * t391 + t611;
t303 = -0.3e1 * t391 + t400;
t487 = 0.8e1 * t568;
t246 = t303 * t487;
t505 = t391 + t304;
t271 = -t386 + t505;
t300 = pkin(1) * t348;
t278 = t300 + pkin(7);
t296 = t304 ^ 2;
t395 = t397 ^ 2;
t350 = 15 * t395;
t357 = 0.18e2 * t400;
t358 = -0.2e1 * t386;
t360 = -0.6e1 * t386;
t362 = 0.2e1 * t391;
t399 = t400 ^ 2;
t371 = 0.3e1 * t399;
t388 = t408 ^ 2;
t219 = t260 + t262;
t596 = -t391 + t400;
t731 = 4 * t397;
t283 = t596 * t731;
t582 = 0.4e1 * t300;
t155 = pkin(7) * t219 * t582 + t283 * t313 + t167;
t583 = 0.2e1 * t300;
t282 = pkin(7) * t583;
t306 = -(3 * t397) + t400;
t558 = 0.4e1 * t618;
t230 = t282 + t558 + t306;
t566 = 0.8e1 * t638;
t579 = 0.6e1 * t297;
t435 = t155 * t579 + t230 * t566;
t535 = 0.12e2 * t618;
t636 = t309 * t391;
t539 = 0.12e2 * t636;
t375 = 0.3e1 * t400;
t604 = (15 * t397) + t375;
t733 = 6 * pkin(1);
t120 = t153 * t539 + t246 + t181 * t535 + t388 + (-t386 + t391 + t604) * t409 + (t350 + (t357 + t360 + 0.6e1 * t391) * t397 + t371 + (t358 + t362) * t400) * t391 + t296 * t271 + t435 * t278 + (t167 * t701 - t182 * t273) * t733;
t298 = pkin(3) * t347;
t545 = t397 * t298;
t440 = -t345 * t404 + t545;
t241 = 0.2e1 * t440;
t577 = 0.2e1 * t298;
t280 = pkin(7) * t577;
t369 = 2 * t397;
t302 = t369 + t391;
t305 = -t397 + t400;
t169 = t305 * t298 + t241 * t313 + (t280 * t348 + t302 * t345) * pkin(1);
t368 = 3 * t397;
t600 = t368 + t400;
t287 = t391 + t600;
t245 = t287 * t298;
t288 = 0.3e1 * t391 + t304;
t642 = t288 * t345;
t201 = -pkin(1) * t642 + t245;
t730 = 8 * t397;
t556 = t408 * t730;
t685 = pkin(3) * t395;
t276 = t556 + 0.4e1 * t685;
t521 = t404 * t625;
t204 = t276 * t347 + 0.4e1 * t345 * t521;
t366 = 5 * t395;
t593 = t399 + t409;
t352 = 10 * t397;
t603 = t352 + t376;
t616 = t400 * t397;
t225 = t391 * t603 + t366 + t593 + 0.6e1 * t616;
t361 = 0.5e1 * t409;
t373 = 0.6e1 * t400;
t236 = t361 + (t352 + t373) * t391 + t296;
t525 = t278 * t638;
t473 = -0.8e1 * t525;
t565 = -0.4e1 * t636;
t643 = t278 * t344;
t299 = pkin(1) * t345;
t249 = -t299 + t298;
t471 = t313 * t545;
t590 = 0.2e1 * t701;
t166 = -0.2e1 * t471 + t245 + (t249 * t590 - t642) * pkin(1);
t719 = -0.4e1 * t166;
t133 = t169 * t565 + t204 * t313 + (t643 * t719 + (-t225 + t487) * t347) * pkin(3) + (-0.4e1 * t201 * t701 + (t236 + t473) * t345) * pkin(1);
t243 = t297 + t278;
t716 = 0.2e1 * t344;
t591 = pkin(7) * t716;
t279 = pkin(3) * t591;
t592 = t400 - t386;
t286 = t397 + t592;
t240 = t279 + t286;
t206 = t260 + t240;
t277 = t297 + pkin(7);
t257 = t277 * t348;
t314 = t391 * t731;
t620 = t391 * t400;
t284 = t314 - 0.4e1 * t620;
t377 = -0.2e1 * t400;
t378 = 0.2e1 * pkin(3);
t531 = -0.4e1 * pkin(3) * pkin(7) * t286;
t229 = t279 + t596 + 0.2e1 * t636;
t646 = t229 * t313;
t143 = t284 * t309 + t344 * t531 - t395 - (t400 - (t378 + pkin(4)) * pkin(4)) * (t400 + (t378 - pkin(4)) * pkin(4)) + (t377 + 0.2e1 * t386 - 0.4e1 * t391 - 0.4e1 * t646) * t397 + (-t206 * t257 + t240 * t273) * t734;
t401 = sqrt(t143);
t109 = t120 * t243 + t133 * t401;
t736 = t109 ^ 2;
t724 = -2 * pkin(1);
t735 = 2 * pkin(1);
t379 = 2 * pkin(2);
t385 = t386 ^ 2;
t595 = t395 + t399;
t599 = t376 - t386;
t617 = t400 * t386;
t436 = t599 * t397 + t385 / 0.6e1 + t595 - t617;
t430 = 0.5e1 / 0.6e1 * t409 + t436;
t189 = (t315 + t599) * t391 + t430;
t732 = -0.6e1 * t189;
t224 = 0.4e1 / 0.3e1 * t618 + t282 + t295;
t727 = t375 - t386 - t391;
t269 = t727 * t352;
t359 = -0.5e1 * t386;
t601 = t359 - 0.5e1 * t391;
t307 = t309 ^ 2;
t639 = t307 * t409;
t641 = t296 * (-t391 + t286);
t729 = 0.7e1 * t388 + ((35 * t397) + 0.15e2 * t400 + t601) * t409 + ((21 * t395) + t269 + 0.9e1 * t399 + (t360 - 0.6e1 * t391) * t400) * t391 + t641 - 0.24e2 * t224 * t639;
t325 = -t386 / 0.4e1;
t728 = t325 + t391 / 0.2e1;
t329 = -0.3e1 / 0.2e1 * t386;
t608 = t385 / 0.2e1 - t409 / 0.2e1;
t459 = -0.3e1 * t617 + t371 + t608;
t614 = t304 * ((t329 + t376) * t397 - 0.3e1 / 0.2e1 * t617 + t595 + t608) + t388;
t150 = t475 * t732 + (t350 + (t357 - 0.9e1 * t386) * t397 + t459) * t391 + (t329 + t604) * t409 + t614;
t327 = -t386 / 0.2e1;
t267 = t327 + t505;
t480 = -0.4e1 * t273;
t451 = t267 * t480;
t446 = pkin(1) * t451;
t170 = t446 + (t367 + t599) * t391 + t430;
t183 = t267 * t625 + t231;
t585 = pkin(7) * t300;
t529 = 0.6e1 * t585;
t137 = t170 * t529 + t183 * t535 + t150 + t246;
t372 = 0.8e1 * t400;
t207 = t404 * t480 + t314 + (4 * t395) + (t358 + t372) * t397;
t218 = -t397 + t445 + t728;
t154 = t487 + t207 * t313 + t267 * t306 + (t218 * t589 + t479 * t624) * pkin(1);
t726 = t137 * t578 + t154 * t541;
t725 = 0.2e1 * pkin(7);
t238 = t362 + t240;
t553 = t229 * t300;
t171 = t238 * t277 + 0.2e1 * t553;
t173 = t238 * t348 + (0.4e1 * t313 - 0.2e1) * t277 * pkin(1);
t609 = -t273 + t257;
t212 = pkin(1) + t609;
t132 = t171 * t345 + t173 * t298 + t212 * t401;
t342 = cos(pkin(8));
t676 = sin(pkin(8));
t227 = -t342 * t348 - t345 * t676;
t683 = pkin(5) * t227;
t217 = t683 * t724;
t383 = pkin(5) ^ 2;
t208 = t217 + 0.2e1 * t383;
t215 = -pkin(1) + t683;
t226 = t342 * t345 - t348 * t676;
t613 = t217 + t383;
t194 = -(t735 + pkin(5)) * pkin(5) + t613;
t195 = pkin(5) * (t735 - pkin(5)) + t613;
t402 = sqrt(-t194 * t195);
t148 = -pkin(5) * t226 * t208 - t215 * t402;
t658 = t132 * t148;
t270 = t362 + t368 + t592;
t202 = t270 + t279 + t460;
t258 = t273 - pkin(1);
t628 = t347 * t348;
t546 = pkin(3) * t628;
t644 = t277 * t345;
t214 = t546 + t644;
t649 = t214 * t401;
t131 = -t202 * t257 + t649 + (t258 * t591 + t270 * t630) * pkin(3) + (-0.2e1 * t646 + (0.2e1 * t309 - 0.4e1) * t391 - t286) * pkin(1);
t648 = t226 * t402;
t146 = pkin(5) * t648 - t208 * t215;
t660 = t131 * t146;
t438 = t660 / 0.4e1 + t658 / 0.4e1;
t551 = pkin(1) * t257;
t179 = t260 + t279 + t505 + 0.2e1 * t551;
t177 = 0.1e1 / t179;
t392 = 0.1e1 / pkin(3);
t203 = t397 + t613;
t199 = 0.1e1 / t203;
t652 = t199 / pkin(5);
t526 = t392 * t652;
t465 = t177 * t526;
t104 = t438 * t465;
t657 = t146 * t132;
t659 = t131 * t148;
t437 = t659 / 0.4e1 - t657 / 0.4e1;
t105 = t437 * t465;
t346 = sin(pkin(9));
t349 = cos(pkin(9));
t98 = t104 * t349 - t105 * t346;
t722 = 0.1e1 / t98;
t721 = 0.1e1 / t146;
t720 = 0.1e1 / t98 ^ 2;
t718 = 0.8e1 * t243;
t715 = -0.2e1 * t348;
t381 = (pkin(6) ^ 2);
t714 = 2 * t381;
t713 = -0.6e1 * t401;
t712 = pkin(1) * pkin(7);
t704 = pkin(6) * t98;
t93 = t704 * t379;
t663 = t381 + t93;
t86 = -(t379 + pkin(6)) * pkin(6) + t663;
t87 = pkin(6) * (t379 - pkin(6)) + t663;
t679 = t86 * t87;
t403 = sqrt(-t679);
t64 = 0.1e1 / t403;
t711 = -t64 / 0.2e1;
t254 = t400 + t391 / 0.4e1 + t397 / 0.4e1 - t386 / 0.8e1;
t605 = 0.4e1 / 0.7e1 * t400 - t386 / 0.7e1;
t164 = -0.32e2 / 0.21e2 * t254 * t475 + 0.5e1 / 0.42e2 * t409 + (0.16e2 / 0.21e2 * t397 + t605) * t391 + t395 / 0.7e1 + t605 * t397 + t399 - 0.3e1 / 0.7e1 * t617 + t385 / 0.42e2;
t256 = t607 + t728;
t338 = 0.4e1 / 0.3e1 * t397;
t165 = -0.8e1 / 0.3e1 * t256 * t475 + 0.5e1 / 0.18e2 * t409 + (t338 + t326) * t391 + t399 - t395 / 0.3e1 + t385 / 0.18e2 + (t332 + 0.2e1 / 0.3e1 * t397 + t328) * t400;
t220 = -t409 / 0.6e1 + t436;
t340 = t397 / 0.2e1;
t606 = t340 + t400;
t222 = -0.2e1 / 0.3e1 * t475 + t325 + t606;
t374 = 0.4e1 * t400;
t290 = (t374 + t386) * t397;
t337 = -0.2e1 / 0.3e1 * t391;
t294 = t337 + t400;
t478 = 0.16e2 * t568;
t311 = t313 ^ 2;
t619 = t395 * t311;
t561 = 0.8e1 * t619;
t136 = t294 * t561 + t222 * t478 + 0.14e2 * t164 * t618 + t305 * t409 + (t290 - 0.10e2 / 0.3e1 * t395 + 0.2e1 * t399 - t617) * t391 + t220 * t624 + (0.6e1 * t165 * t701 + t295 * t451) * pkin(1);
t594 = t399 - t395;
t156 = t293 * t446 - t388 + (-t315 - t592) * t409 + (t290 + t409 / 0.6e1 - t385 / 0.6e1 + t594) * t391 + t220 * t400;
t365 = 7 * t395;
t160 = (t329 + t375 + (7 * t397)) * t409 + (t365 + (t359 + 0.10e2 * t400) * t397 + t459) * t391 + t614;
t413 = pkin(7) * t400;
t285 = -0.12e2 * pkin(7) * t404 + t413 * t734;
t292 = -(8 * t395) + 0.12e2 * t616;
t174 = t285 * t348 + t292 * t313 + t478 + t561 + t595 - 0.6e1 * t616;
t185 = t260 * t625 + t267 * t303;
t247 = 0.16e2 * (t593 - 0.6e1 * t620) * t395;
t301 = -0.30e2 * t386 + 0.60e2 * t400;
t598 = t385 - t409;
t441 = 0.6e1 * t399 + t598 - 0.6e1 * t617;
t103 = t247 * t311 + 0.32e2 * t185 * t568 + 0.24e2 * t156 * t618 + (t358 + t374 + (28 * t397)) * t388 + t271 * t641 + (0.24e2 * t136 * t309 + t301 * t395 + t399 * t360 + t441 * t367 + t598 * t376 + (28 * t404 ^ 2) + 0.4e1 * t413 ^ 2) * t391 + t726 * t278 + 0.8e1 * (t150 * t701 - t160 * t273) * pkin(1) + (0.16e2 * t174 * t307 + t301 * t397 + (70 * t395) + t409 + t441) * t409;
t252 = 0.7e1 / 0.6e1 * t391 + t324 + t606;
t334 = t391 / 0.3e1;
t510 = t324 + t334 + t400;
t255 = t338 + t510;
t180 = -t252 * t299 + t255 * t298;
t264 = t397 + t510;
t289 = t362 + t305;
t191 = t264 * t298 - t289 * t299 / 0.2e1;
t511 = t386 / 0.3e1 + t334 + t376;
t458 = -0.8e1 / 0.3e1 * t619 + t391 * t305 - 0.5e1 / 0.3e1 * t395 + t511 * t397 + t400 * (t326 + t293);
t700 = pkin(7) * t395;
t588 = -0.4e1 * t700;
t316 = -0.20e2 / 0.3e1 * t397;
t512 = 0.2e1 / 0.3e1 * t386 + t333 + t374;
t513 = 0.4e1 / 0.3e1 * t386 + t332 + t377;
t692 = t409 / 0.2e1 - (t316 + t512) * t391 / 0.2e1 + 0.3e1 / 0.2e1 * t395 - t513 * t397 / 0.2e1 - t399 / 0.2e1;
t138 = t345 * t312 * t588 + t180 * t558 + t458 * t298 + (t191 * t589 + t345 * t692) * pkin(1);
t602 = t358 - 0.2e1 * t391;
t506 = t373 + t602;
t176 = t409 + (t328 + t337 + t603) * t391 + t366 + t506 * t397 + t400 * (t328 + t294);
t172 = t176 * t298;
t261 = 0.8e1 / 0.3e1 * t391 + t509;
t263 = t326 + t333 + t600;
t184 = -t261 * t299 + t263 * t298;
t268 = 0.5e1 / 0.6e1 * t391 + t340 + t324;
t197 = t268 * t577 + t299 * t625;
t544 = t404 * t298;
t470 = t312 * t544;
t560 = -0.4e1 * t618;
t190 = t361 + (t352 + t506) * t391 + (t337 + t507) * t304;
t653 = t190 * t345;
t139 = -0.8e1 * pkin(7) * t470 + t197 * t560 + t172 + (t184 * t589 - t653) * pkin(1);
t157 = -pkin(1) * t653 + t172;
t188 = -0.3e1 * t409 + (t316 + t513) * t391 + t512 * t397 + t594;
t192 = -0.5e1 / 0.3e1 * t409 + (-t397 + t511) * t391 + t400 * (t336 + t509);
t584 = -0.2e1 * t299;
t158 = t188 * t298 + t192 * t584;
t242 = 0.4e1 * t440;
t248 = 0.2e1 * t299 + t298;
t266 = t327 + t287;
t161 = t306 * t298 + t242 * t313 + (t248 * t701 + t266 * t345) * t735;
t253 = t400 + 0.5e1 / 0.2e1 * t391 + 0.3e1 / 0.2e1 * t397 + t327;
t640 = t303 * t345;
t193 = t253 * t298 + pkin(1) * t640 / 0.2e1;
t448 = 0.24e2 * t293 * t619 - t388 - ((21 * t397) + t727) * t409 - (t400 * t602 + t269 + t371 + (35 * t395)) * t391 - (t365 + (t372 + t601) * t397 + t400 * (-t391 + t592)) * t304;
t540 = -0.12e2 * t636;
t110 = t193 * t478 + t161 * t473 + t138 * t540 - 0.6e1 * t158 * t618 + (-0.6e1 * t139 * t643 + t347 * t448) * pkin(3) + (-0.6e1 * t157 * t701 + t729 * t345) * pkin(1);
t83 = t103 * t243 + t110 * t401;
t710 = t83 / 0.4e1;
t407 = pkin(2) ^ 2;
t90 = t407 + t663;
t709 = pkin(2) / t90 ^ 2;
t693 = t177 / 0.2e1;
t496 = t392 * t693;
t113 = qJ(1) + atan2(t132 * t496, t131 * t496);
t387 = 0.1e1 / pkin(4);
t622 = t387 * t392;
t489 = t622 / 0.2e1;
t107 = 0.1e1 / t109;
t673 = t107 * t83;
t70 = atan2(t401 * t489, t489 * t673) + t113;
t66 = sin(t70);
t708 = pkin(4) * t66;
t67 = cos(t70);
t707 = pkin(4) * t67;
t100 = t104 * t346 + t105 * t349;
t682 = pkin(6) * t100;
t91 = t93 + t714;
t92 = -pkin(2) - t704;
t39 = -t403 * t682 - t91 * t92;
t40 = -t403 * t92 + t682 * t91;
t382 = 0.1e1 / pkin(6);
t515 = t382 / t90 / 0.2e1;
t20 = atan2(t40 * t515, t39 * t515) + t113;
t18 = sin(t20);
t706 = pkin(6) * t18;
t19 = cos(t20);
t705 = pkin(6) * t19;
t703 = pkin(7) * t313;
t702 = pkin(7) * t345;
t699 = t107 / 0.2e1;
t629 = t345 * t348;
t645 = t277 * t313;
t151 = (t229 * t629 + t298 * t645) * t730 + (t206 * t644 + t240 * t546) * t734;
t310 = t345 ^ 2;
t142 = 0.1e1 / t401;
t695 = t142 / 0.2e1;
t499 = t212 * t695;
t117 = -t649 + t151 * t499 + t229 * t310 * t724 + t171 * t348 + (-t238 - 0.8e1 * t551) * t273;
t698 = t117 / 0.4e1;
t697 = t131 / 0.4e1;
t696 = t132 / 0.4e1;
t694 = t146 / 0.4e1;
t691 = -t401 / 0.4e1;
t690 = t401 / 0.4e1;
t689 = t403 / 0.2e1;
t688 = 0.2e1 * t226 ^ 2;
t687 = pkin(1) * t295;
t686 = pkin(1) * t313;
t684 = pkin(3) * t401;
t681 = t160 * pkin(1);
t632 = t344 * t348;
t233 = t630 - t632;
t633 = t344 * t345;
t232 = t628 + t633;
t518 = t673 / 0.4e1;
t434 = t131 * t690 + t132 * t518;
t621 = t387 / pkin(3) ^ 2;
t527 = t177 * t621;
t72 = t434 * t527;
t68 = t232 * t72;
t433 = t131 * t518 + t132 * t691;
t71 = t433 * t527;
t44 = -t233 * t71 - t68;
t42 = 0.1e1 / t44 ^ 2;
t69 = t233 * t72;
t43 = -t232 * t71 + t69;
t680 = t42 * t43;
t498 = t214 * t695;
t631 = t344 * t391;
t119 = t609 * t401 + t151 * t498 + (t202 * t277 + 0.4e1 * t553) * t345 + (t590 * t631 + (t270 * t348 + t645 * t734) * pkin(3)) * t347;
t205 = t214 * t735;
t178 = 0.1e1 / t179 ^ 2;
t428 = t178 * t434;
t672 = 0.1e1 / t736 * t83;
t516 = -t672 / 0.4e1;
t461 = t132 * t516;
t452 = pkin(7) * t471;
t251 = -0.4e1 * t452;
t281 = pkin(7) * t584;
t453 = t544 * t703;
t486 = 0.32e2 / 0.3e1 * t395;
t454 = t312 * t486;
t455 = 0.64e2 / 0.3e1 * t254 * t404;
t463 = t312 * t521;
t467 = 0.32e2 * t525;
t474 = t625 * t700;
t615 = t404 * t313;
t572 = pkin(7) * t615;
t481 = -0.24e2 * t572;
t482 = -0.48e2 * t572;
t488 = pkin(7) * t558;
t500 = t151 * t695;
t626 = t348 * t397;
t520 = t345 * t626;
t536 = -0.24e2 * t626;
t537 = -0.32e2 * t312 * t395;
t538 = 0.24e2 * t636;
t542 = -0.96e2 * t293 * t312;
t550 = pkin(1) * t624;
t552 = t267 * t687;
t559 = -0.2e1 * t618;
t562 = 0.8e1 * t626;
t563 = -0.4e1 * t634;
t567 = -0.8e1 * t638;
t573 = pkin(7) * t618;
t581 = -0.6e1 * t297;
t612 = -0.24e2 * t293 * t470 + t481 * t640;
t666 = -0.4e1 * t712;
t667 = -0.6e1 * t712;
t50 = ((t266 * t583 + t488 + t563) * t473 + (0.12e2 * t310 * t313 * t700 + t252 * t563 + t311 * t588) * t540 + 0.12e2 * t192 * t634 + (t692 * t540 + t729) * t300 + (t289 * t559 * t540 + t190 * t557 + t303 * t561) * pkin(7)) * t401 + t110 * t500 + t251 * t243 * t467 + (((pkin(7) * t261 * t560 - t190 * t300 - 0.4e1 * t463) * t713 + t612 * t718) * t643 + ((-pkin(7) * t311 * t486 - 0.16e2 * t256 * t573 - t312 * t455 - 0.4e1 * t348 * t552) * t538 - 0.64e2 * t311 * t474 + t404 * t267 * t542 - 0.48e2 * t189 * t573 - 0.8e1 * t348 * t681 + ((t550 * t715 + t563) * t541 + (-0.24e2 * t267 * t573 + t300 * t732) * t578) * t278) * t243 * t347) * pkin(3) + ((0.16e2 * (t292 * t715 - t285 + t482 + t537) * t639 + (t207 * t715 + t218 * t666 + t481) * t467 + (-0.28e2 * t164 * t626 + t165 * t667 + t222 * t482 + t294 * t537) * t538 + t278 * (t170 * t667 + t183 * t536) * t578 - 0.4e1 * t247 * t312 - 0.96e2 * t185 * t572 - 0.48e2 * t156 * t626 - 0.8e1 * t150 * t712) * t243 + ((-0.8e1 * t180 * t626 + t298 * t454) * t540 + t347 * t542 * t685 + t193 * t482 + 0.12e2 * t158 * t626 + (0.2e1 * (-t242 * t348 - t248 * t712) * t567 + (t184 * t666 + t197 * t562 + 0.24e2 * t453) * t581) * t278) * t401 + (-t726 * t243 - t103 + (t139 * t579 + t161 * t566 + (-0.24e2 * t281 + 0.64e2 * t520) * t639 + (0.48e2 * t191 * t636 + 0.6e1 * t157) * pkin(7)) * t401) * pkin(1)) * t345;
t503 = t131 * t142 / 0.8e1;
t477 = pkin(1) * t546;
t564 = 0.8e1 * t636;
t580 = -0.4e1 * t297;
t84 = (t310 * t308 * t556 + (t300 * t302 - 0.2e1 * t634) * t565 + 0.4e1 * t463 + t288 * t488 + t236 * t300) * t401 + t133 * t500 + ((-0.8e1 / 0.3e1 * t470 + t251 - 0.2e1 * t295 * t477) * t539 - 0.24e2 * t262 * t452 - 0.6e1 * t182 * t477 + t612) * t243 + ((t241 * t348 * t564 + t204 * t715 - 0.24e2 * t453) * t401 + (0.12e2 * (-t211 * t626 - t572) * t539 + t181 * t536) * t243 + (0.4e1 * t166 * t344 * t684 - t435 * t243 - t120 + ((t298 * t564 + 0.4e1 * t201) * t401 + (-0.48e2 * t223 * t636 - 0.6e1 * t167) * t243) * pkin(7)) * pkin(1)) * t345 + ((t567 * t300 + ((0.4e1 * t273 * t348 - 0.2e1 * t703) * t397 + (-0.2e1 * t249 * t702 - t288 * t348) * pkin(1)) * t580) * t401 + ((t281 - 0.8e1 * t520) * t566 + (-0.8e1 * t452 - 0.2e1 * t283 * t629 + (-t219 * t702 - t262 * t546) * t734) * t579) * t243) * t278;
t678 = (t205 * t428 + (t119 * t690 + t151 * t503 + t84 * t461 + (t50 * t696 + t83 * t698) * t107) * t177) * t621 + t71;
t627 = t347 * t391;
t522 = t344 * t627;
t235 = t280 + 0.4e1 * t522;
t571 = t391 * t702;
t548 = pkin(3) * t633;
t476 = pkin(1) * t548;
t610 = 0.2e1 * t476 + t280;
t140 = t235 * t560 + (t284 * t716 + t531) * t347 + (-t610 * t257 + 0.2e1 * t347 ^ 2 * t571 + (-t206 * t628 - t240 * t633) * pkin(3)) * t734;
t115 = t140 * t499 + t235 * t345 * t583 + ((t345 * t401 - t173) * t344 + (t348 * t401 + (t277 * t725 + t238) * t345 + (-pkin(1) + 0.2e1 * t686 + t701) * t577) * t347) * pkin(3);
t547 = pkin(3) * t632;
t116 = (t273 - t547) * t401 + t140 * t498 - 0.2e1 * t235 * t686 - (t280 + 0.4e1 * t476) * t257 - 0.2e1 * t309 * t571 - t270 * t548 + (t631 * t734 + (-t202 * t348 + t258 * t725) * pkin(3)) * t347;
t198 = 0.2e1 * t477 + t610;
t523 = t293 * t615;
t221 = 0.24e2 * t523 * t548;
t484 = pkin(7) * t547;
t443 = t397 * t345 * t484;
t239 = 0.4e1 * t443;
t250 = t484 * t724;
t501 = t140 * t695;
t635 = t309 * t408;
t524 = t278 * t635;
t530 = -0.4e1 * t585;
t570 = pkin(7) * t626;
t637 = t308 * t409;
t60 = t110 * t501 + (0.32e2 * t239 * t243 - 0.8e1 * t250 * t401) * t525 + (0.24e2 * (-0.4e1 * t224 * t299 * t637 - t138 * t631 - t161 * t524) * t401 + (0.48e2 * t136 * t631 + 0.96e2 * t154 * t524 + 0.64e2 * t174 * t637) * t243) * t347 + ((t103 + (t137 * t718 + t139 * t713) * t278) * t347 + (t278 * t221 * t718 + ((t255 * t560 + t264 * t530 - t458) * t540 - 0.16e2 * t253 * t568 + t188 * t557 + t176 * t529 + ((-t306 + t560) * t567 + (t263 * t530 + 0.8e1 * t268 * t618 - t176 + t487) * t581) * t278 - t448) * t401 + ((pkin(7) * t454 + 0.16e2 * t256 * t570 + t313 * t455 + 0.4e1 * t552) * t538 + 0.64e2 * t312 * t474 + 0.96e2 * t267 * t523 + 0.48e2 * t189 * t570 + 0.8e1 * t681 + ((0.2e1 * t550 + 0.4e1 * t615) * t541 + (t189 * t733 + 0.24e2 * t267 * t570) * t578) * t278) * t243 * t345) * t344) * pkin(3);
t85 = (t250 * t565 + (-0.8e1 * t169 * t627 - t276 * t313 + ((-t305 + t559) * t565 + t225 + (t287 * t582 - 0.8e1 * t634) * pkin(7)) * pkin(3)) * t344 + ((t250 + (-t287 + 0.2e1 * t618) * t297) * t580 + (pkin(3) * t719 - 0.24e2 * t299 * t635) * t347) * t278) * t401 + t133 * t501 + (0.24e2 * t230 * t347 * t524 + (t239 + (0.8e1 / 0.3e1 * t615 + 0.2e1 * t687) * t548) * t539 + 0.24e2 * t153 * t522 + t221 + 0.24e2 * t262 * t443 + 0.6e1 * t182 * t476 + 0.6e1 * ((pkin(7) * t562 + t262 * t734) * t345 * t636 + t155 * t298) * t278) * t243 + t120 * t298;
t677 = (-t198 * t428 + (t116 * t690 + t140 * t503 + t85 * t461 + (t115 * t710 + t60 * t696) * t107) * t177) * t621 - t71;
t675 = t100 * t720;
t674 = t736 / t83 ^ 2;
t111 = sin(t113);
t130 = 0.1e1 / t131 ^ 2;
t472 = pkin(3) / (t130 * t132 ^ 2 + 0.1e1) * t179 * t392;
t439 = -0.2e1 * t130 * t132 * t472;
t442 = 0.2e1 / t131 * t472;
t655 = t178 * t198;
t495 = -t655 / 0.2e1;
t661 = t116 * t177;
t662 = t115 * t177;
t61 = (t662 / 0.2e1 + t132 * t495) * t442 + (t661 / 0.2e1 + t131 * t495) * t439;
t671 = t111 * t61;
t654 = t178 * t205;
t493 = t654 / 0.2e1;
t62 = 0.1e1 + (t117 * t693 + t132 * t493) * t442 + (t119 * t693 + t131 * t493) * t439;
t670 = t111 * t62;
t112 = cos(t113);
t669 = t112 * t61;
t668 = t112 * t62;
t665 = pkin(3) * t670 + t299;
t664 = pkin(2) * t670 + t299;
t485 = t226 * t735;
t656 = 0.1e1 / t402 * (t194 + t195) * pkin(5) * t485;
t651 = t199 / pkin(1);
t200 = 0.1e1 / t203 ^ 2;
t650 = t200 * t226;
t647 = t227 * t402;
t623 = t382 / pkin(2);
t497 = t656 / 0.2e1;
t121 = (-t647 + (t215 * t735 - t208 + t497) * t226) * pkin(5);
t122 = -t215 * t656 / 0.2e1 + t383 * pkin(1) * t688 + (t208 * t227 - t648) * pkin(5);
t464 = t652 * t654;
t554 = pkin(1) * t650;
t51 = (t438 * t464 + ((t660 / 0.2e1 + t658 / 0.2e1) * t554 + (t119 * t694 + t121 * t697 + t122 * t696 + t148 * t698) * t652) * t177) * t392;
t52 = (-t437 * t464 + ((-t659 / 0.2e1 + t657 / 0.2e1) * t554 + (-t119 * t148 / 0.4e1 - t131 * t122 / 0.4e1 + t121 * t696 + t117 * t694) * t652) * t177) * t392;
t35 = t346 * t52 + t349 * t51;
t576 = t35 * t709;
t494 = -t655 / 0.4e1;
t431 = t661 / 0.4e1 + t131 * t494;
t432 = t662 / 0.4e1 + t132 * t494;
t75 = (t146 * t431 + t148 * t432) * t526;
t76 = (t146 * t432 - t148 * t431) * t526;
t48 = t346 * t76 + t349 * t75;
t575 = t48 * t709;
t38 = 0.1e1 / t39 ^ 2;
t574 = pkin(6) / (t38 * t40 ^ 2 + 0.1e1) * t90;
t569 = g(2) * t669;
t543 = pkin(5) * t650;
t534 = 0.1e1 / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t407 * t720 * t679 / 0.4e1) * t623;
t533 = t92 * t711;
t532 = -t64 * t722 / 0.4e1;
t519 = t100 * t711;
t517 = -t672 / 0.2e1;
t514 = t720 * t689;
t504 = -0.2e1 * pkin(2) * t92 + t91;
t502 = -t132 * t142 / 0.8e1;
t492 = t652 / 0.2e1;
t491 = -t651 / 0.2e1;
t490 = t651 / 0.2e1;
t483 = pkin(2) * t100 * t714;
t444 = pkin(6) * (-t86 - t87) * t379;
t34 = t48 * t444;
t449 = -0.2e1 * t38 * t40 * t574;
t466 = 0.2e1 / t39 * t574;
t469 = pkin(6) * t515;
t49 = t346 * t75 - t349 * t76;
t4 = ((t34 * t533 + t48 * t483 + (t403 * t48 + t49 * t91) * pkin(6)) * t515 - t40 * t575) * t466 + (-t39 * t575 + (t34 * t519 - t403 * t49 + t48 * t504) * t469) * t449 + t61;
t74 = 0.1e1 / (t143 * t674 + 0.1e1);
t429 = -0.2e1 * pkin(4) * t74 * t622 * t674 * t684;
t447 = t109 * t74 / t83 * t695;
t17 = t140 * t447 + (t85 * t517 + t60 * t699) * t429 + t61;
t21 = t35 * t444;
t36 = t346 * t51 - t349 * t52;
t2 = ((t21 * t533 + t35 * t483 + (t35 * t403 + t36 * t91) * pkin(6)) * t515 - t40 * t576) * t466 + (-t39 * t576 + (t21 * t519 + t35 * t504 - t403 * t36) * t469) * t449 + t62;
t462 = t131 * t516;
t16 = t151 * t447 + (t50 * t699 + t84 * t517) * t429 + t62;
t457 = -pkin(2) * t668 - t300;
t456 = -pkin(3) * t668 - t300;
t450 = g(1) * t299 - g(2) * t300;
t427 = t178 * t433;
t426 = g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (rSges(4,1) * t19 - rSges(4,2) * t18);
t425 = g(1) * (rSges(5,1) * t66 + rSges(5,2) * t67) + g(2) * (-rSges(5,1) * t67 + rSges(5,2) * t66);
t59 = atan2(t100, t98) + t113;
t57 = sin(t59);
t58 = cos(t59);
t424 = g(1) * (-rSges(7,1) * t57 - rSges(7,2) * t58) + g(2) * (rSges(7,1) * t58 - rSges(7,2) * t57);
t14 = atan2(t623 * t689, -t98) + t20;
t12 = sin(t14);
t13 = cos(t14);
t423 = g(1) * (rSges(10,1) * t12 + rSges(10,2) * t13) + g(2) * (-rSges(10,1) * t13 + rSges(10,2) * t12);
t30 = atan2(t43, t44) + t70;
t28 = sin(t30);
t29 = cos(t30);
t422 = g(1) * (-rSges(11,1) * t28 - rSges(11,2) * t29) + g(2) * (rSges(11,1) * t29 - rSges(11,2) * t28);
t421 = g(1) * (rSges(3,1) * t111 + rSges(3,2) * t112) + g(2) * (-rSges(3,1) * t112 + rSges(3,2) * t111);
t216 = -pkin(1) * t227 + pkin(5);
t209 = t217 + t369;
t149 = -pkin(1) * t226 * t209 + t216 * t402;
t147 = pkin(1) * t648 + t209 * t216;
t145 = 0.1e1 / t147 ^ 2;
t144 = 0.1e1 / t146 ^ 2;
t128 = qJ(1) + atan2(t148 * t492, t146 * t492);
t127 = pkin(8) + atan2(t149 * t490, t147 * t491);
t126 = cos(t128);
t125 = sin(t128);
t124 = cos(t127);
t123 = sin(t127);
t65 = 0.1e1 / (t100 ^ 2 * t720 + 0.1e1);
t54 = pkin(2) * t671;
t53 = pkin(3) * t671;
t41 = 0.1e1 / t44;
t33 = 0.1e1 / (t42 * t43 ^ 2 + 0.1e1);
t25 = (-t198 * t427 + (t85 * t462 + t115 * t691 + t140 * t502 + (t116 * t710 + t60 * t697) * t107) * t177) * t621;
t23 = (t205 * t427 + (t84 * t462 + t117 * t691 + t151 * t502 + (t119 * t710 + t50 * t697) * t107) * t177) * t621;
t1 = [-m(2) * (g(1) * (rSges(2,1) * t345 + rSges(2,2) * t348) + g(2) * (-rSges(2,1) * t348 + rSges(2,2) * t345)) - m(3) * (t421 * t62 + t450) - m(4) * (g(1) * t664 + g(2) * t457 + t2 * t426) - m(5) * (g(1) * t665 + g(2) * t456 + t16 * t425) + m(6) * (g(1) * (-rSges(6,1) * t123 - rSges(6,2) * t124) + g(2) * (rSges(6,1) * t124 - rSges(6,2) * t123)) * (-((t216 * t497 + t397 * pkin(5) * t688 + (t209 * t227 - t648) * pkin(1)) * t490 + t149 * t543) / t147 - (-t147 * t543 + (-t647 + (-0.2e1 * t216 * pkin(5) - t209 + t497) * t226) * pkin(1) * t491) * t149 * t145) / (t145 * t149 ^ 2 + 0.1e1) * t203 * t724 - m(7) * (t424 * ((-t35 * t675 + t36 * t722) * t65 + t62) + t450) - m(9) * ((g(1) * (-rSges(9,1) * t125 - rSges(9,2) * t126) + g(2) * (rSges(9,1) * t126 - rSges(9,2) * t125)) * (0.1e1 + (t122 * t721 * t652 + (-t121 * t144 * t652 + (-t146 * t144 + t721) * t200 * t485) * t148) * t203 / (t144 * t148 ^ 2 + 0.1e1) * pkin(5)) + t450) - m(10) * (g(1) * (-t2 * t706 + t664) + g(2) * (t2 * t705 + t457) + t423 * ((t21 * t532 + t35 * t514) * t534 + t2)) - m(11) * (g(1) * (t16 * t708 + t665) + g(2) * (-t16 * t707 + t456) + t422 * (((t678 * t233 + (-t23 + t72) * t232) * t41 - (-t23 * t233 - t232 * t678 + t69) * t680) * t33 + t16)), -m(3) * t421 * t61 - m(4) * (-pkin(2) * t569 + g(1) * t54 + t4 * t426) - m(5) * (-pkin(3) * t569 + g(1) * t53 + t17 * t425) - m(7) * t424 * ((-t48 * t675 + t49 * t722) * t65 + t61) - m(8) * (g(1) * (rSges(8,1) * t347 - rSges(8,2) * t344) + g(2) * (rSges(8,1) * t344 + rSges(8,2) * t347)) - m(10) * (g(1) * (-t4 * t706 + t54) + g(2) * (-pkin(2) * t669 + t4 * t705) + t423 * ((t34 * t532 + t48 * t514) * t534 + t4)) - m(11) * (g(1) * (t17 * t708 + t53) + g(2) * (-pkin(3) * t669 - t17 * t707) + t422 * (((-t232 * t25 + t233 * t677 - t68) * t41 - ((-t25 - t72) * t233 - t677 * t232) * t680) * t33 + t17))];
taug = t1(:);