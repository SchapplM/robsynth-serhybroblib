% Calculate Gravitation load on the joints for
% picker2Dm2DE1
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
% Datum: 2020-05-09 18:54
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 14:32:11
% EndTime: 2020-05-09 14:39:17
% DurationCPUTime: 78.78s
% Computational Cost: add. (1596568->800), mult. (4827102->1318), div. (26008->22), fcn. (841754->28), ass. (0->573)
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t344 = pkin(1) ^ 2;
t355 = pkin(4) ^ 2;
t357 = pkin(7) ^ 2;
t594 = t357 - t355;
t535 = t344 + t594;
t309 = cos(qJ(2));
t625 = t307 * t309;
t568 = pkin(1) * t625;
t658 = sin(qJ(2));
t569 = t658 * pkin(7);
t419 = pkin(1) * (0.2e1 * (-t568 + t569) * pkin(3) + t535);
t554 = pkin(3) * t658;
t494 = t554 + pkin(7);
t417 = t494 * t419;
t516 = pkin(7) * t554;
t466 = t357 + 0.2e1 * t516;
t450 = t344 + t466;
t442 = -t355 + t450;
t665 = pkin(1) * pkin(3);
t433 = t442 * t665;
t424 = 0.4e1 * t309 * t433;
t350 = pkin(3) ^ 2;
t365 = t658 ^ 2;
t617 = t350 * t365;
t571 = 0.2e1 * t617;
t434 = -t350 + t466 + t571;
t427 = t344 * t434;
t263 = pkin(3) * t309;
t277 = t310 ^ 2;
t455 = t277 * t494;
t441 = t455 * t263;
t670 = 0.8e1 * t344;
t675 = 0.4e1 * t307;
t136 = t441 * t670 + t417 * t675 + (0.8e1 * t307 * t427 + t424) * t310;
t337 = 0.1e1 / pkin(4);
t342 = t344 ^ 2;
t666 = pkin(3) * pkin(7);
t456 = -0.4e1 * t535 * t666;
t599 = t344 - t357;
t667 = 0.4e1 * t350;
t464 = t599 * t667;
t668 = 0.2e1 * t350;
t533 = t668 + t594;
t683 = 0.1e1 / pkin(3);
t643 = t683 / 0.2e1;
t713 = 0.1e1 / t643;
t380 = -0.4e1 * t277 * t427 - 0.4e1 * t310 * t417 + t307 * t424 + t365 * t464 + t658 * t456 - t342 - 0.2e1 * t533 * t344 - (t357 - (t713 + pkin(4)) * pkin(4)) * (t357 + (t713 - pkin(4)) * pkin(4));
t346 = sqrt(t380);
t644 = -t683 / 0.2e1;
t564 = pkin(3) * t625;
t515 = pkin(1) * t564;
t609 = t344 / 0.3e1 + t357;
t189 = -0.4e1 / 0.9e1 * t515 + 0.4e1 / 0.9e1 * t350 - t355 / 0.9e1 + t609;
t287 = -t355 / 0.6e1;
t296 = 0.2e1 / 0.3e1 * t350;
t465 = t357 - t515;
t196 = t287 + t296 + t465;
t295 = 0.4e1 / 0.3e1 * t350;
t269 = t344 + t357;
t289 = -t355 / 0.3e1;
t539 = t289 + t269;
t234 = t295 + t539;
t261 = -t344 / 0.3e1 + t357;
t518 = -0.2e1 * t564;
t618 = t344 * t277;
t574 = 0.6e1 * t618;
t276 = t310 * t277;
t347 = pkin(1) * t344;
t627 = t276 * t347;
t584 = pkin(7) * t627;
t662 = pkin(7) * t310;
t592 = 0.4e1 * t662;
t620 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t137 = 0.4e1 * t584 + t189 * t574 + t234 * t620 + (t196 * t592 + t261 * t518) * pkin(1);
t328 = 0.6e1 * t344;
t291 = -0.2e1 / 0.3e1 * t355;
t335 = 0.2e1 * t357;
t538 = t291 + t296 + t335;
t352 = t350 ^ 2;
t537 = t291 + t269;
t612 = (t296 + t537) * t269 + t352;
t149 = -0.4e1 * t234 * t515 + (t328 + t538) * t350 + t612;
t299 = -t350 / 0.3e1;
t259 = t299 + t357;
t495 = -0.2e1 * t515;
t208 = t259 * t495;
t621 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t168 = t234 * t621 + t208;
t278 = 0.10e2 / 0.3e1 * t344;
t169 = (t278 + t538) * t350 + t612;
t268 = -0.3e1 * t350 + t357;
t525 = 0.8e1 * t584;
t221 = t268 * t525;
t242 = t350 + t535;
t265 = pkin(1) * t310;
t247 = t265 + pkin(7);
t262 = t269 ^ 2;
t312 = 0.15e2 * t342;
t319 = 0.18e2 * t357;
t320 = -0.2e1 * t355;
t322 = -0.6e1 * t355;
t345 = t357 ^ 2;
t330 = 0.3e1 * t345;
t351 = pkin(3) * t350;
t338 = t351 ^ 2;
t192 = t234 + t495;
t251 = 0.4e1 * (-t350 + t357) * t344;
t589 = pkin(7) * t265;
t139 = 0.4e1 * t192 * t589 + t251 * t277 + t149;
t587 = 0.2e1 * t265;
t250 = pkin(7) * t587;
t271 = -0.3e1 * t344 + t357;
t575 = 0.4e1 * t618;
t207 = t250 + t575 + t271;
t512 = 0.6e1 * t554;
t273 = t658 * t365;
t630 = t273 * t351;
t580 = 0.8e1 * t630;
t439 = t139 * t512 + t207 * t580;
t557 = 0.12e2 * t617;
t559 = 0.12e2 * t618;
t590 = 0.6e1 * pkin(1);
t334 = 0.3e1 * t357;
t606 = 0.15e2 * t344 + t334;
t122 = t137 * t557 + t221 + t168 * t559 + t338 + (t350 - t355 + t606) * t352 + (t312 + (t319 + t322 + 0.6e1 * t350) * t344 + t330 + (t320 + t668) * t357) * t350 + t262 * t242 + t439 * t247 + (t149 * t662 - t169 * t564) * t590;
t565 = t344 * t263;
t452 = -t307 * t347 + t565;
t216 = 0.2e1 * t452;
t586 = t309 * t713;
t248 = pkin(7) * t586;
t267 = 0.2e1 * t344 + t350;
t151 = -t599 * t263 + t216 * t277 + (t310 * t248 + t267 * t307) * pkin(1);
t671 = 0.3e1 * t344;
t595 = t357 + t671;
t253 = t350 + t595;
t220 = t253 * t263;
t254 = 0.3e1 * t350 + t269;
t634 = t254 * t307;
t184 = -pkin(1) * t634 + t220;
t573 = t351 * t670;
t672 = 0.4e1 * t342;
t246 = pkin(3) * t672 + t573;
t546 = t347 * t621;
t185 = t246 * t309 + t546 * t675;
t327 = 0.5e1 * t342;
t597 = t345 + t352;
t314 = 0.10e2 * t344;
t605 = t314 + t335;
t616 = t357 * t344;
t201 = t605 * t350 + t327 + t597 + 0.6e1 * t616;
t323 = 0.5e1 * t352;
t332 = 0.6e1 * t357;
t213 = t323 + (t314 + t332) * t350 + t262;
t550 = t247 * t630;
t510 = -0.8e1 * t550;
t532 = t247 * t658;
t572 = -0.4e1 * t617;
t264 = pkin(1) * t307;
t224 = -t264 + t263;
t509 = t277 * t565;
t593 = 0.2e1 * t662;
t148 = -0.2e1 * t509 + t220 + (t224 * t593 - t634) * pkin(1);
t711 = -0.4e1 * t148;
t126 = t151 * t572 + t185 * t277 + (t532 * t711 + (-t201 + t525) * t309) * pkin(3) + (-0.4e1 * t184 * t662 + (t213 + t510) * t307) * pkin(1);
t218 = t265 + t494;
t109 = t122 * t218 + t126 * t346;
t714 = t109 ^ 2;
t229 = t357 + t350 / 0.4e1 + t344 / 0.4e1 - t355 / 0.8e1;
t336 = t355 ^ 2;
t607 = 0.4e1 / 0.7e1 * t357 - t355 / 0.7e1;
t615 = t357 * t355;
t146 = -0.32e2 / 0.21e2 * t229 * t515 + 0.5e1 / 0.42e2 * t352 + (0.16e2 / 0.21e2 * t344 + t607) * t350 + t342 / 0.7e1 + t607 * t344 + t345 - 0.3e1 / 0.7e1 * t615 + t336 / 0.42e2;
t288 = -t355 / 0.4e1;
t687 = t288 + t350 / 0.2e1;
t231 = t609 + t687;
t301 = 0.4e1 / 0.3e1 * t344;
t147 = -0.8e1 / 0.3e1 * t231 * t515 + 0.5e1 / 0.18e2 * t352 + (t301 + t289) * t350 + t345 - t342 / 0.3e1 + t336 / 0.18e2 + (t295 + 0.2e1 / 0.3e1 * t344 + t291) * t357;
t600 = t342 + t345;
t602 = t335 - t355;
t445 = t602 * t344 + t336 / 0.6e1 + t600 - t615;
t193 = -t352 / 0.6e1 + t445;
t303 = t344 / 0.2e1;
t608 = t303 + t357;
t195 = -0.2e1 / 0.3e1 * t515 + t288 + t608;
t333 = 0.4e1 * t357;
t256 = (t333 + t355) * t344;
t300 = -0.2e1 / 0.3e1 * t350;
t260 = t300 + t357;
t290 = -t355 / 0.2e1;
t239 = t290 + t350 + t269;
t709 = -0.4e1 * t239;
t482 = t564 * t709;
t521 = 0.16e2 * t584;
t275 = t277 ^ 2;
t629 = t275 * t342;
t579 = 0.8e1 * t629;
t127 = t260 * t579 + t195 * t521 + 0.14e2 * t146 * t618 - t599 * t352 + (t256 - 0.10e2 / 0.3e1 * t342 + 0.2e1 * t345 - t615) * t350 + t193 * t620 + (0.6e1 * t147 * t662 + t261 * t482) * pkin(1);
t440 = 0.5e1 / 0.6e1 * t352 + t445;
t176 = (t278 + t602) * t350 + t440;
t292 = -0.3e1 / 0.2e1 * t355;
t610 = t336 / 0.2e1 - t352 / 0.2e1;
t493 = -0.3e1 * t615 + t330 + t610;
t614 = t269 * ((t292 + t335) * t344 - 0.3e1 / 0.2e1 * t615 + t600 + t610) + t338;
t135 = -0.6e1 * t176 * t515 + (t312 + (t319 - 0.9e1 * t355) * t344 + t493) * t350 + (t292 + t606) * t352 + t614;
t467 = pkin(1) * t482;
t598 = t345 - t342;
t140 = t259 * t467 - t338 + (-t278 - t594) * t352 + (t256 + t352 / 0.6e1 - t336 / 0.6e1 + t598) * t350 + t193 * t357;
t321 = -0.5e1 * t355;
t326 = 0.7e1 * t342;
t144 = (t292 + t334 + 0.7e1 * t344) * t352 + (t326 + (t321 + 0.10e2 * t357) * t344 + t493) * t350 + t614;
t358 = pkin(7) * t357;
t681 = 0.4e1 * pkin(1);
t252 = -0.12e2 * pkin(7) * t347 + t358 * t681;
t673 = -0.8e1 * t342;
t258 = t673 + 0.12e2 * t616;
t154 = t252 * t310 + t258 * t277 + t521 + t579 + t600 - 0.6e1 * t616;
t172 = t239 * t268 + t495 * t621;
t669 = -0.6e1 * t350;
t222 = 0.16e2 * (t357 * t669 + t597) * t342;
t266 = -0.30e2 * t355 + 0.60e2 * t357;
t272 = t365 ^ 2;
t152 = t467 + (t328 + t602) * t350 + t440;
t170 = t239 * t621 + t208;
t128 = 0.6e1 * t152 * t589 + t170 * t559 + t135 + t221;
t331 = 0.8e1 * t357;
t507 = t347 * t564;
t187 = -0.4e1 * t507 + t672 + (t667 + t320 + t331) * t344;
t191 = -t344 + t465 + t687;
t138 = t525 + t187 * t277 + t239 * t271 + (t191 * t592 + t518 * t620) * pkin(1);
t511 = 0.8e1 * t554;
t562 = 0.32e2 * t630;
t437 = t128 * t511 + t138 * t562;
t601 = t336 - t352;
t453 = 0.6e1 * t345 + t601 - 0.6e1 * t615;
t534 = -t350 + t594;
t633 = t262 * (t344 + t534);
t105 = t222 * t275 + 0.32e2 * t172 * t584 + 0.24e2 * t140 * t618 + (t320 + t333 + 0.28e2 * t344) * t338 + t242 * t633 + (0.24e2 * t127 * t365 + t266 * t342 + t345 * t322 + t453 * t328 + t601 * t335 + 0.28e2 * t347 ^ 2 + 0.4e1 * t358 ^ 2) * t350 + t437 * t247 + 0.8e1 * (t135 * t662 - t144 * t564) * pkin(1) + (0.16e2 * t154 * t272 + t266 * t344 + 0.70e2 * t342 + t352 + t453) * t352;
t227 = 0.7e1 / 0.6e1 * t350 + t287 + t608;
t297 = t350 / 0.3e1;
t540 = t287 + t297 + t357;
t230 = t301 + t540;
t167 = -t227 * t264 + t230 * t263;
t541 = t355 / 0.3e1 + t297 + t335;
t173 = -t350 * t599 - 0.5e1 / 0.3e1 * t342 + t541 * t344 + t357 * (t289 + t259);
t236 = t344 + t540;
t255 = t668 - t599;
t178 = t236 * t263 - t255 * t264 / 0.2e1;
t628 = t276 * t342;
t279 = -0.20e2 / 0.3e1 * t344;
t542 = 0.2e1 / 0.3e1 * t355 + t296 + t333;
t543 = 0.4e1 / 0.3e1 * t355 + t295 - 0.2e1 * t357;
t660 = t352 / 0.2e1 - (t279 + t542) * t350 / 0.2e1 + 0.3e1 / 0.2e1 * t342 - t543 * t344 / 0.2e1 - t345 / 0.2e1;
t663 = pkin(7) * t307;
t129 = -0.4e1 * t628 * t663 + t167 * t575 + (-0.8e1 / 0.3e1 * t629 + t173) * t263 + (t178 * t592 + t307 * t660) * pkin(1);
t604 = t320 - 0.2e1 * t350;
t536 = t332 + t604;
t160 = t352 + (t291 + t300 + t605) * t350 + t327 + t536 * t344 + t357 * (t291 + t260);
t153 = t160 * t263;
t177 = t323 + (t314 + t536) * t350 + (t300 + t537) * t269;
t636 = t177 * t307;
t141 = -pkin(1) * t636 + t153;
t175 = -0.3e1 * t352 + (t279 + t543) * t350 + t542 * t344 + t598;
t179 = -0.5e1 / 0.3e1 * t352 + (-t344 + t541) * t350 + t357 * (t299 + t539);
t588 = -0.2e1 * t264;
t142 = t175 * t263 + t179 * t588;
t686 = t334 - t350 - t355;
t241 = t686 * t314;
t603 = t321 - 0.5e1 * t350;
t143 = t338 + (0.21e2 * t344 + t686) * t352 + (t604 * t357 + t241 + t330 + 0.35e2 * t342) * t350 + (t326 + (t331 + t603) * t344 + t357 * t534) * t269;
t217 = 0.4e1 * t452;
t223 = 0.2e1 * t264 + t263;
t238 = t290 + t253;
t591 = 0.2e1 * pkin(1);
t145 = t271 * t263 + t217 * t277 + (t223 * t662 + t238 * t307) * t591;
t228 = t357 + 0.5e1 / 0.2e1 * t350 + 0.3e1 / 0.2e1 * t344 + t290;
t632 = t268 * t307;
t180 = t228 * t263 + pkin(1) * t632 / 0.2e1;
t558 = -0.12e2 * t617;
t197 = 0.4e1 / 0.3e1 * t618 + t250 + t261;
t631 = t272 * t352;
t688 = 0.7e1 * t338 + (0.35e2 * t344 + 0.15e2 * t357 + t603) * t352 + (0.21e2 * t342 + t241 + 0.9e1 * t345 + (t322 + t669) * t357) * t350 + t633 - 0.24e2 * t197 * t631;
t233 = 0.8e1 / 0.3e1 * t350 + t539;
t235 = t289 + t296 + t595;
t171 = -t233 * t264 + t235 * t263;
t240 = 0.5e1 / 0.6e1 * t350 + t303 + t287;
t548 = t307 * t621;
t182 = pkin(1) * t548 + t240 * t586;
t566 = t276 * t263;
t508 = t347 * t566;
t576 = -0.4e1 * t618;
t130 = -0.8e1 * pkin(7) * t508 + t182 * t576 + t153 + (t171 * t592 - t636) * pkin(1);
t712 = -0.6e1 * t130;
t110 = t180 * t521 + t145 * t510 + t129 * t558 - 0.6e1 * t142 * t618 + (t532 * t712 + (0.24e2 * t259 * t629 - t143) * t309) * pkin(3) + (-0.6e1 * t141 * t662 + t688 * t307) * pkin(1);
t92 = t105 * t218 + t110 * t346;
t647 = t714 / t92 ^ 2;
t77 = 0.1e1 / (t380 * t647 + 0.1e1);
t444 = pkin(4) * t337 * t346 * t77 / t644 * t647;
t133 = 0.1e1 / t346;
t661 = t133 / 0.2e1;
t470 = t109 * t77 / t92 * t661;
t108 = 0.1e1 / t714;
t645 = t108 * t92;
t506 = t644 * t645;
t107 = 0.1e1 / t109;
t552 = t107 * t643;
t483 = pkin(7) * t509;
t226 = -0.4e1 * t483;
t249 = pkin(7) * t588;
t274 = t307 ^ 2;
t522 = 0.32e2 / 0.3e1 * t342;
t484 = t276 * t522;
t485 = 0.64e2 / 0.3e1 * t229 * t347;
t496 = t276 * t546;
t626 = t277 * t347;
t582 = pkin(7) * t626;
t520 = -0.48e2 * t582;
t524 = pkin(7) * t575;
t526 = t136 * t661;
t622 = t310 * t344;
t544 = t307 * t622;
t556 = 0.24e2 * t617;
t560 = -0.24e2 * t626;
t561 = -0.32e2 * t628;
t623 = t309 * t310;
t563 = pkin(3) * t623;
t624 = t307 * t310;
t577 = -0.2e1 * t624;
t578 = -0.4e1 * t627;
t581 = -0.8e1 * t630;
t583 = pkin(7) * t618;
t585 = pkin(7) * t629;
t519 = pkin(7) * t560;
t708 = -0.24e2 * t259;
t613 = t508 * t708 + t519 * t632;
t664 = pkin(7) * t277;
t674 = -0.2e1 * t310;
t63 = ((0.12e2 * t274 * t342 * t664 + t227 * t578 - 0.2e1 * t255 * t583 - 0.4e1 * t585) * t558 + 0.12e2 * t179 * t627 + (t660 * t558 + t688) * t265 + (t177 * t574 + t268 * t579) * pkin(7)) * t346 + t110 * t526 + ((-pkin(7) * t275 * t522 + t261 * t265 * t709 - 0.16e2 * t231 * t583 - t276 * t485) * t556 - 0.64e2 * t585 * t621 - 0.96e2 * t259 * t239 * t627 - 0.48e2 * t176 * t583 - 0.8e1 * t144 * t265) * t218 * t263 + (0.12e2 * (-(-0.8e1 * t167 * t622 + t484 * t263) * t617 + t259 * t566 * t673 - 0.4e1 * t180 * t582 + t142 * t622) * t346 + (0.16e2 * (t258 * t674 - t252 + t520 + t561) * t631 + (-0.28e2 * t146 * t622 + t195 * t520 + t260 * t561) * t556 - 0.4e1 * t222 * t276 - 0.96e2 * t172 * t582 - 0.48e2 * t140 * t622) * t218 + (-t437 * t218 - t105 + (t130 * t512 + t145 * t580 + (-0.24e2 * t249 + 0.64e2 * t544) * t631) * t346 + ((0.48e2 * t178 * t617 + 0.6e1 * t141) * t346 + (-0.144e3 * t147 * t617 - 0.8e1 * t135) * t218) * pkin(7)) * pkin(1)) * t307 + (((t524 + t217 * t577 + t578 + (-t223 * t663 + t238 * t310) * t591) * t581 - 0.6e1 * (-0.4e1 * t496 + (t182 * t307 * t670 - pkin(1) * t177) * t310 + (-0.4e1 * t171 * t264 + (-0.4e1 * t233 * t344 + 0.24e2 * t507) * t277) * pkin(7)) * t554) * t346 + ((t226 + (-0.2e1 * t620 * t265 + t578) * t263 + (t187 * t674 + (-0.4e1 * pkin(1) * t191 + t560) * pkin(7)) * t307) * t562 + (0.24e2 * (-t239 * t263 * t664 - t170 * t624) * t344 + (-t152 * t663 - t176 * t563) * t590 + t613) * t511) * t218) * t247;
t513 = 0.4e1 * t554;
t514 = pkin(1) * t563;
t570 = 0.8e1 * t617;
t93 = (t274 * t273 * t573 + (t267 * t265 - 0.2e1 * t627) * t572 + 0.4e1 * t496 + t254 * t524 + t213 * t265) * t346 + t126 * t526 + ((-0.8e1 / 0.3e1 * t508 + t226 - 0.2e1 * t261 * t514) * t557 - 0.24e2 * t234 * t483 - 0.6e1 * t169 * t514 + t613) * t218 + ((t216 * t310 * t570 + t185 * t674 + t519 * t263) * t346 + (0.12e2 * (-t189 * t622 - t582) * t557 - 0.24e2 * t168 * t622) * t218 + (t148 * t346 * t513 - t439 * t218 - t122 + ((t570 * t263 + 0.4e1 * t184) * t346 + (-0.48e2 * t196 * t617 - 0.6e1 * t149) * t218) * pkin(7)) * pkin(1)) * t307 + ((t581 * t265 - 0.4e1 * ((t563 * t675 - 0.2e1 * t664) * t344 + (-0.2e1 * t224 * t663 - t254 * t310) * pkin(1)) * t554) * t346 + ((t249 - 0.8e1 * t544) * t580 + (-0.8e1 * t483 + t251 * t577 + (-t192 * t663 - t234 * t563) * t681) * t512) * t218) * t247;
t20 = t136 * t470 + (t93 * t506 + t63 * t552) * t444;
t492 = t671 + t533;
t421 = (-0.4e1 * t568 + 0.2e1 * t569) * pkin(3) + t492;
t418 = t494 * t421;
t430 = pkin(1) * t434;
t425 = -0.2e1 * t430;
t446 = pkin(7) * (-pkin(1) + t564) * t713;
t458 = pkin(3) * t492;
t448 = t309 * t458;
t431 = t307 * t494 + t563;
t689 = t431 * t346;
t407 = t277 * t425 + t307 * t448 - t310 * t418 + t689 + t658 * t446 + (t571 - t667 - t535) * pkin(1);
t454 = t494 * t310;
t432 = t454 - t564;
t166 = t432 * t591 + t350 + t450;
t679 = 0.1e1 / t166;
t400 = t679 * t407;
t392 = t400 * t643;
t423 = t310 * t430;
t438 = t668 + t442;
t415 = t438 * t494 + 0.2e1 * t423;
t682 = -0.2e1 * pkin(1);
t416 = pkin(3) * (t310 * t438 + t455 * t681 + t494 * t682);
t426 = pkin(1) + t432;
t412 = t307 * t415 + t309 * t416 + t346 * t426;
t411 = t683 * t412;
t408 = t411 / 0.2e1;
t403 = t679 * t408;
t378 = atan2(t403, t392);
t376 = sin(t378);
t377 = cos(t378);
t104 = t307 * t376 - t310 * t377;
t420 = t426 * t661;
t449 = pkin(1) * t454;
t118 = -t689 + t136 * t420 + t274 * t425 + t415 * t310 + (-0.8e1 * t449 - t438) * t564;
t422 = t431 * t661;
t490 = t350 * t658 * t309;
t119 = t136 * t422 + t307 * t418 + t310 * t448 + t346 * t432 + t423 * t675 + t441 * t681 + t490 * t593;
t405 = 0.1e1 / t407 ^ 2;
t379 = 0.1e1 / (t405 * t412 ^ 2 + 0.1e1);
t374 = t379 * t405 * t412;
t375 = t379 / t407;
t372 = t118 * t375 - t119 * t374;
t368 = t372 * t376;
t369 = t372 * t377;
t41 = -t307 * t368 + t310 * t369 - t104;
t103 = -t307 * t377 - t310 * t376;
t42 = t307 * t369 + t310 * t368 - t103;
t551 = t337 * t643;
t646 = t107 * t92;
t76 = atan2(t346 * t551, t551 * t646);
t74 = sin(t76);
t75 = cos(t76);
t517 = -t103 * t74 + t104 * t75;
t3 = -t20 * t517 + t41 * t75 - t42 * t74;
t529 = t658 * t307;
t209 = t529 + t623;
t528 = t658 * t310;
t210 = -t528 + t625;
t396 = t400 / 0.4e1;
t410 = t679 * t412;
t619 = t337 / t350;
t73 = (t346 * t396 + t410 * t646 / 0.4e1) * t619;
t71 = t210 * t73;
t391 = t107 * t396;
t659 = -t346 / 0.4e1;
t72 = (t92 * t391 + t410 * t659) * t619;
t59 = -t209 * t72 + t71;
t70 = t209 * t73;
t60 = -t210 * t72 - t70;
t34 = atan2(t59, t60);
t31 = sin(t34);
t32 = cos(t34);
t50 = t103 * t75 + t104 * t74;
t429 = -t20 * t50 + t41 * t74 + t42 * t75;
t728 = -t3 * t31 - t32 * t429;
t727 = t3 * t32 - t31 * t429;
t473 = 0.4e1 * t490;
t212 = t473 + t248;
t479 = t529 * t665;
t611 = 0.2e1 * t479 + t248;
t690 = t350 * t309 ^ 2;
t131 = 0.8e1 * pkin(7) * t264 * t690 + t212 * t576 - 0.4e1 * t419 * t563 - 0.4e1 * t433 * t529 - 0.4e1 * t611 * t449 + (0.2e1 * t464 * t658 + t456) * t309;
t530 = t347 * t658;
t488 = t277 * t530;
t657 = pkin(3) * t307;
t443 = t259 * t488 * t657;
t194 = 0.24e2 * t443;
t491 = t344 * t528;
t474 = pkin(7) * t491;
t447 = t474 * t657;
t215 = 0.4e1 * t447;
t503 = pkin(1) * t528;
t463 = t503 * t666;
t225 = -0.2e1 * t463;
t481 = t276 * t516;
t457 = t347 * t481;
t635 = t247 * t309;
t497 = t351 * t365 * t635;
t471 = -0.24e2 * t497;
t531 = t277 * t658;
t489 = t344 * t531;
t472 = -0.4e1 * t489;
t480 = -0.4e1 * t503;
t487 = t658 * t629;
t501 = t350 * t532;
t502 = pkin(3) * t532;
t555 = pkin(1) * t658;
t504 = t261 * t555;
t527 = t131 * t661;
t547 = t307 * t617;
t549 = t273 * t309 * t352;
t567 = pkin(3) * t635;
t710 = 0.6e1 * t176;
t66 = (-0.96e2 * t197 * t549 * t264 + t225 * t510 + t145 * t471 - 0.24e2 * t129 * t490 - 0.6e1 * (0.8e1 * t240 * t489 - t658 * t160 + (t235 * t480 + 0.8e1 * t276 * t530) * pkin(7)) * t501 + t567 * t712 - 0.16e2 * t228 * t457 + 0.6e1 * t160 * t463 + t143 * t554 + ((-t658 * t271 + t472) * t510 + (0.8e1 / 0.3e1 * t487 + t230 * t472 + t236 * pkin(7) * t480 - t173 * t658) * t558 + t487 * t708 + 0.6e1 * t175 * t489) * pkin(3)) * t346 + t110 * t527 + 0.8e1 * (0.8e1 * t154 * t549 + 0.4e1 * (t215 + (0.2e1 * t555 * t620 + 0.4e1 * t488) * t657) * t550 + 0.12e2 * t138 * t497 + 0.3e1 * (t485 * t531 + 0.4e1 * t239 * t504 + (0.16e2 * t231 * t491 + t658 * t484) * pkin(7)) * pkin(3) * t547 + 0.6e1 * t127 * t490 + (t194 + (0.24e2 * t239 * t474 + t555 * t710) * t657) * t502 + t128 * t567 + 0.8e1 * t342 * t481 * t548 + 0.12e2 * t239 * t443 + t447 * t710 + t144 * t479) * t218 + t105 * t263;
t94 = (t471 * t264 + (t225 + (t599 * t658 - 0.2e1 * t489) * pkin(3)) * t572 - 0.8e1 * t151 * t490 - 0.4e1 * (t225 + (-t658 * t253 + 0.2e1 * t489) * pkin(3)) * t502 + t567 * t711 - 0.8e1 * t457 - t246 * t531 + 0.4e1 * t253 * t463 + t201 * t554) * t346 + t126 * t527 + (0.24e2 * t207 * t497 + (t215 + (0.8e1 / 0.3e1 * t488 + 0.2e1 * t504) * t657) * t557 + 0.24e2 * t137 * t490 + 0.6e1 * (0.4e1 * t234 * t555 + 0.8e1 * t474) * t307 * t501 + 0.6e1 * t139 * t567 + t194 + 0.24e2 * t234 * t447 + 0.6e1 * t169 * t479) * t218 + t122 * t263;
t21 = t131 * t470 + (t94 * t506 + t66 * t552) * t444;
t656 = pkin(3) * t346;
t116 = t209 * t656 + t131 * t420 + (t212 * t587 + (t344 + (t513 + 0.2e1 * pkin(7)) * pkin(7) + t533) * t263) * t307 - t658 * t416 + (t593 + (0.4e1 * t277 - 0.2e1) * pkin(1)) * t690;
t117 = t210 * t656 + t131 * t422 + t212 * t277 * t682 - t421 * t563 - (0.4e1 * t479 + t248) * t454 + pkin(1) * t473 - 0.2e1 * pkin(7) * t547 + t309 * t446 - t458 * t529;
t373 = t116 * t375 - t117 * t374;
t370 = t373 * t376;
t371 = t373 * t377;
t47 = t307 * t370 - t310 * t371;
t48 = t307 * t371 + t310 * t370;
t428 = -t21 * t50 - t47 * t74 + t48 * t75;
t9 = -t21 * t517 - t47 * t75 - t48 * t74;
t726 = -t31 * t9 - t32 * t428;
t725 = -t31 * t428 + t32 * t9;
t308 = sin(pkin(9));
t311 = cos(pkin(9));
t305 = cos(pkin(8));
t652 = sin(pkin(8));
t202 = t305 * t307 - t652 * t310;
t204 = t305 * t310 + t652 * t307;
t111 = t202 * t392 + t204 * t403;
t409 = -t411 / 0.2e1;
t691 = t202 * t679;
t112 = t204 * t392 + t409 * t691;
t101 = t111 * t308 + t112 * t311;
t680 = 0.1e1 / t101;
t100 = -t111 * t311 + t112 * t308;
t677 = 0.1e1 / t101 ^ 2;
t693 = t100 * t677;
t398 = t407 * t644;
t163 = 0.1e1 / t166 ^ 2;
t186 = t431 * t591;
t637 = t163 * t186;
t387 = t398 * t637;
t404 = t204 * t409;
t553 = t679 * t643;
t498 = t204 * t553;
t499 = t202 * t553;
t86 = t118 * t498 + t119 * t499 - t202 * t387 - t404 * t637 + t112;
t402 = t202 * t408;
t500 = t644 * t691;
t87 = t118 * t500 + t119 * t498 - t204 * t387 - t402 * t637 - t111;
t716 = (t308 * t86 + t311 * t87) * t693 - (t308 * t87 - t311 * t86) * t680;
t82 = 0.1e1 / (t100 ^ 2 * t677 + 0.1e1);
t25 = t716 * t82;
t84 = atan2(t100, t101);
t78 = sin(t84);
t80 = cos(t84);
t468 = t103 * t78 - t104 * t80;
t14 = -t25 * t468 + t41 * t80 - t42 * t78;
t648 = t104 * t78;
t436 = -t41 * t78 - t25 * t648 + (-t103 * t25 - t42) * t80;
t85 = atan2(t100, -t101);
t79 = sin(t85);
t81 = cos(t85);
t724 = t14 * t79 - t436 * t81;
t723 = t14 * t81 + t436 * t79;
t183 = 0.2e1 * t514 + t611;
t638 = t163 * t183;
t88 = t117 * t499 + t116 * t498 + (t202 * t398 + t404) * t638;
t89 = t117 * t498 + t116 * t500 + (t204 * t398 + t402) * t638;
t715 = (t308 * t88 + t311 * t89) * t693 - (t308 * t89 - t311 * t88) * t680;
t29 = t715 * t82;
t16 = -t29 * t468 - t47 * t80 - t48 * t78;
t435 = t47 * t78 - t29 * t648 + (-t103 * t29 - t48) * t80;
t722 = t16 * t79 - t435 * t81;
t721 = t16 * t81 + t435 * t79;
t720 = t31 * t50 - t32 * t517;
t719 = t31 * t517 + t32 * t50;
t55 = t103 * t80 + t648;
t718 = t468 * t81 + t55 * t79;
t717 = t468 * t79 - t55 * t81;
t692 = t133 * t679;
t397 = t407 * t619;
t390 = -t397 / 0.4e1;
t386 = t163 * t390;
t406 = t412 * t619;
t401 = t406 / 0.4e1;
t685 = t163 * t346 * t401 + t386 * t646;
t399 = -t92 * t406 / 0.4e1;
t684 = t107 * t163 * t399 + t346 * t386;
t58 = 0.1e1 / t60 ^ 2;
t655 = t58 * t59;
t384 = t397 * t692 / 0.8e1;
t388 = t108 * t679 * t399;
t395 = t107 * t679 * t401;
t545 = t679 * t619;
t486 = t545 / 0.4e1;
t451 = t486 * t646;
t461 = t346 * t486;
t654 = t118 * t451 + t119 * t461 + t136 * t384 - t186 * t684 + t388 * t93 + t395 * t63 + t72;
t653 = t116 * t451 + t117 * t461 + t131 * t384 + t183 * t684 + t388 * t94 + t395 * t66 - t72;
t642 = t42 * pkin(3) + t264;
t641 = t42 * pkin(2) + t264;
t640 = -t41 * pkin(3) - t265;
t639 = -t41 * pkin(2) - t265;
t462 = t545 * t659;
t414 = t82 * (g(1) * (rSges(10,1) * t717 + rSges(10,2) * t718) + g(2) * (-rSges(10,1) * t718 + rSges(10,2) * t717));
t413 = 0.1e1 / (t58 * t59 ^ 2 + 0.1e1) * (g(1) * (t719 * rSges(11,1) - rSges(11,2) * t720) + g(2) * (rSges(11,1) * t720 + t719 * rSges(11,2)));
t394 = -t406 * t692 / 0.8e1;
t385 = t391 * t619;
t381 = t679 * t390 * t645;
t165 = atan2(-t202, -t204);
t159 = cos(t165);
t157 = sin(t165);
t124 = -t157 * t305 - t652 * t159;
t123 = -t652 * t157 + t159 * t305;
t57 = 0.1e1 / t60;
t46 = t48 * pkin(2);
t45 = t48 * pkin(3);
t44 = t47 * pkin(2);
t43 = t47 * pkin(3);
t27 = t116 * t462 + t117 * t451 + t131 * t394 + t183 * t685 + t381 * t94 + t385 * t66;
t23 = t118 * t462 + t119 * t451 + t136 * t394 - t186 * t685 + t381 * t93 + t385 * t63;
t1 = [-m(2) * (g(1) * (rSges(2,1) * t307 + rSges(2,2) * t310) + g(2) * (-rSges(2,1) * t310 + rSges(2,2) * t307)) - m(3) * (g(1) * (rSges(3,1) * t42 + rSges(3,2) * t41 + t264) + g(2) * (-rSges(3,1) * t41 + rSges(3,2) * t42 - t265)) - m(4) * (g(1) * (rSges(4,1) * t436 - rSges(4,2) * t14 + t641) + g(2) * (rSges(4,1) * t14 + rSges(4,2) * t436 + t639)) - m(5) * (g(1) * (rSges(5,1) * t429 + rSges(5,2) * t3 + t642) + g(2) * (-rSges(5,1) * t3 + rSges(5,2) * t429 + t640)) - m(6) * (g(1) * (rSges(6,1) * t124 - rSges(6,2) * t123) + g(2) * (rSges(6,1) * t123 + rSges(6,2) * t124)) - m(7) * (g(1) * (rSges(7,1) * t436 - rSges(7,2) * t14 + t264) + g(2) * (rSges(7,1) * t14 + rSges(7,2) * t436 - t265)) - m(9) * (g(1) * t264 - g(2) * t265) - m(10) * (g(1) * (t436 * pkin(6) + t724 * rSges(10,1) + t723 * rSges(10,2) + t641) + g(2) * (t14 * pkin(6) - t723 * rSges(10,1) + t724 * rSges(10,2) + t639) + t716 * t414) - m(11) * (g(1) * (t429 * pkin(4) + t728 * rSges(11,1) - t727 * rSges(11,2) + t642) + g(2) * (-t3 * pkin(4) + t727 * rSges(11,1) + t728 * rSges(11,2) + t640) + ((t654 * t210 + (-t23 + t73) * t209) * t57 - (-t654 * t209 - t210 * t23 + t71) * t655) * t413), -m(3) * (g(1) * (rSges(3,1) * t48 - rSges(3,2) * t47) + g(2) * (rSges(3,1) * t47 + rSges(3,2) * t48)) - m(4) * (g(1) * (rSges(4,1) * t435 - rSges(4,2) * t16 + t46) + g(2) * (rSges(4,1) * t16 + rSges(4,2) * t435 + t44)) - m(5) * (g(1) * (rSges(5,1) * t428 + rSges(5,2) * t9 + t45) + g(2) * (-rSges(5,1) * t9 + rSges(5,2) * t428 + t43)) - m(7) * (g(1) * (rSges(7,1) * t435 - rSges(7,2) * t16) + g(2) * (rSges(7,1) * t16 + rSges(7,2) * t435)) - m(8) * (g(1) * (t309 * rSges(8,1) - t658 * rSges(8,2)) + g(2) * (t658 * rSges(8,1) + t309 * rSges(8,2))) - m(10) * (g(1) * (t435 * pkin(6) + t722 * rSges(10,1) + t721 * rSges(10,2) + t46) + g(2) * (t16 * pkin(6) - t721 * rSges(10,1) + t722 * rSges(10,2) + t44) + t715 * t414) - m(11) * (g(1) * (t428 * pkin(4) + t726 * rSges(11,1) - t725 * rSges(11,2) + t45) + g(2) * (-t9 * pkin(4) + t725 * rSges(11,1) + t726 * rSges(11,2) + t43) + ((-t209 * t27 + t653 * t210 - t70) * t57 - ((-t27 - t73) * t210 - t653 * t209) * t655) * t413)];
taug = t1(:);
