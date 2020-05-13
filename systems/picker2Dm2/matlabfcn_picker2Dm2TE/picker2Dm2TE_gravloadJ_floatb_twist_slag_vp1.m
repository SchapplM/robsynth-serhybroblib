% Calculate Gravitation load on the joints for
% picker2Dm2TE
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
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm2TE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2TE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2TE_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2TE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 10:37:22
% EndTime: 2020-05-09 10:40:22
% DurationCPUTime: 33.46s
% Computational Cost: add. (438790->796), mult. (1346430->1292), div. (6812->10), fcn. (238044->10), ass. (0->541)
t637 = 0.1e1 / pkin(3);
t608 = -t637 / 0.2e1;
t606 = t637 / 0.2e1;
t615 = sin(qJ(2));
t616 = sin(qJ(1));
t461 = t615 * t616;
t257 = cos(qJ(2));
t258 = cos(qJ(1));
t596 = t257 * t258;
t162 = t461 + t596;
t507 = t615 * t258;
t513 = t257 * t616;
t163 = -t507 + t513;
t227 = t258 ^ 2;
t303 = (pkin(4) ^ 2);
t305 = (pkin(7) ^ 2);
t298 = (pkin(3) ^ 2);
t627 = 2 * t298;
t567 = t305 + t627;
t517 = -t303 + t567;
t294 = (pkin(1) ^ 2);
t630 = 3 * t294;
t471 = t630 + t517;
t482 = pkin(1) * t513;
t546 = t615 * pkin(7);
t498 = 0.2e1 * t546;
t381 = (-0.4e1 * t482 + t498) * pkin(3) + t471;
t533 = pkin(3) * t615;
t473 = t533 + pkin(7);
t378 = t473 * t381;
t441 = pkin(3) * t498 + t305;
t311 = t615 ^ 2;
t590 = t298 * t311;
t549 = 0.2e1 * t590;
t397 = -t298 + t441 + t549;
t394 = pkin(1) * t397;
t393 = -0.2e1 * t394;
t479 = pkin(3) * t513;
t670 = 0.1e1 / t606;
t406 = pkin(7) * (t479 - pkin(1)) * t670;
t437 = pkin(3) * t471;
t408 = t616 * t437;
t569 = t303 - t305;
t519 = t294 - t569;
t626 = 4 * t298;
t290 = t294 ^ 2;
t379 = pkin(1) * (0.2e1 * (-t482 + t546) * pkin(3) + t519);
t377 = t473 * t379;
t424 = t294 + t441;
t411 = -t303 + t424;
t395 = pkin(1) * pkin(3) * t411;
t390 = t616 * t395;
t607 = -t637 / 0.4e1;
t435 = pkin(7) * t519 / t607;
t571 = t294 - t305;
t440 = t571 * t626;
t591 = t294 * t227;
t553 = -0.4e1 * t591;
t633 = 0.4e1 * t257;
t664 = -2 * t294;
t293 = sqrt(t397 * t553 - 0.4e1 * t258 * t377 + t390 * t633 + t311 * t440 + t615 * t435 - t290 + (t517 * t664) - (t305 - (t670 + pkin(4)) * pkin(4)) * (t305 + (t670 - pkin(4)) * pkin(4)));
t414 = t616 * t473;
t543 = pkin(3) * t596;
t392 = t414 + t543;
t661 = t392 * t293;
t370 = t227 * t393 + t257 * t408 - t258 * t378 + t661 + t615 * t406 + (t549 - t626 - t519) * pkin(1);
t430 = t473 * t258;
t391 = -t479 + t430;
t562 = 0.2e1 * pkin(1);
t376 = t391 * t562 + t298 + t424;
t374 = 0.1e1 / t376;
t360 = t374 * t370;
t357 = t360 / 0.4e1;
t216 = pkin(1) * t258;
t170 = t216 + t473;
t455 = pkin(1) * t479;
t583 = t294 / 0.3e1 + t305;
t149 = -0.4e1 / 0.9e1 * t455 + 0.4e1 / 0.9e1 * t298 - t303 / 0.9e1 + t583;
t237 = -t303 / 0.6e1;
t246 = 0.2e1 / 0.3e1 * t298;
t422 = t305 - t455;
t156 = t246 + t237 + t422;
t245 = 0.4e1 / 0.3e1 * t298;
t220 = t294 + t305;
t239 = -t303 / 0.3e1;
t523 = t239 + t220;
t185 = t245 + t523;
t212 = -t294 / 0.3e1 + t305;
t456 = -0.2e1 * t479;
t226 = t258 * t227;
t295 = pkin(1) * t294;
t597 = t226 * t295;
t556 = pkin(7) * t597;
t623 = pkin(7) * t258;
t564 = 0.4e1 * t623;
t594 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t107 = 0.4e1 * t556 + 0.6e1 * t149 * t591 + t185 * t594 + (t156 * t564 + t212 * t456) * pkin(1);
t276 = 6 * t294;
t241 = -0.2e1 / 0.3e1 * t303;
t283 = 2 * t305;
t522 = t241 + t246 + t283;
t300 = t298 ^ 2;
t521 = t241 + t220;
t585 = (t246 + t521) * t220 + t300;
t119 = -0.4e1 * t185 * t455 + (t276 + t522) * t298 + t585;
t249 = -t298 / 0.3e1;
t210 = t249 + t305;
t438 = -0.2e1 * t455;
t161 = t210 * t438;
t595 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t130 = t185 * t595 + t161;
t228 = 0.10e2 / 0.3e1 * t294;
t131 = (t228 + t522) * t298 + t585;
t560 = 0.2e1 * t216;
t201 = pkin(7) * t560;
t222 = -3 * t294 + t305;
t552 = 0.4e1 * t591;
t160 = t201 + t552 + t222;
t219 = -3 * t298 + t305;
t501 = 0.8e1 * t556;
t173 = t219 * t501;
t193 = t298 + t519;
t198 = t216 + pkin(7);
t213 = t220 ^ 2;
t260 = 15 * t290;
t267 = 18 * t305;
t268 = -2 * t303;
t270 = -6 * t303;
t292 = t305 ^ 2;
t278 = 3 * t292;
t299 = pkin(3) * t298;
t286 = t299 ^ 2;
t538 = 0.12e2 * t590;
t540 = 0.12e2 * t591;
t561 = 0.6e1 * pkin(1);
t282 = 3 * t305;
t580 = 15 * t294 + t282;
t224 = t615 * t311;
t598 = t224 * t299;
t152 = t438 + t185;
t202 = 4 * (-t298 + t305) * t294;
t559 = 0.4e1 * t216;
t109 = pkin(7) * t152 * t559 + t202 * t227 + t119;
t668 = 0.6e1 * t109;
t89 = t107 * t538 + t173 + t130 * t540 + t286 + ((t298 - t303 + t580) * t300) + ((t260 + (t267 + t270 + 6 * t298) * t294 + t278 + (t268 + t627) * t305) * t298) + (t213 * t193) + (0.8e1 * t160 * t598 + t533 * t668) * t198 + (t119 * t623 - t131 * t479) * t561;
t508 = t616 * t295;
t214 = pkin(3) * t257;
t542 = t294 * t214;
t412 = -t508 + t542;
t168 = 0.2e1 * t412;
t558 = t257 * t670;
t199 = pkin(7) * t558;
t218 = 2 * t294 + t298;
t121 = -t571 * t214 + t168 * t227 + (t199 * t258 + t218 * t616) * pkin(1);
t568 = t305 + t630;
t204 = t298 + t568;
t172 = t204 * t214;
t205 = 3 * t298 + t220;
t509 = t616 * t205;
t145 = -pkin(1) * t509 + t172;
t629 = 8 * t294;
t551 = t299 * t629;
t632 = 4 * t290;
t197 = pkin(3) * t632 + t551;
t468 = t616 * t595;
t146 = t197 * t257 + 0.4e1 * t295 * t468;
t275 = 5 * t290;
t572 = t292 + t300;
t262 = 10 * t294;
t579 = t262 + t283;
t589 = t305 * t294;
t158 = t298 * t579 + t275 + t572 + 6 * t589;
t271 = 5 * t300;
t280 = 6 * t305;
t165 = t271 + (t262 + t280) * t298 + t213;
t469 = t616 * t598;
t451 = -0.8e1 * t469;
t516 = t198 * t615;
t550 = -0.4e1 * t590;
t215 = pkin(1) * t616;
t176 = -t215 + t214;
t487 = t227 * t542;
t565 = 0.2e1 * t623;
t118 = -0.2e1 * t487 + t172 + (t176 * t565 - t509) * pkin(1);
t667 = -0.4e1 * t118;
t98 = t121 * t550 + t146 * t227 + (t516 * t667 + (-t158 + t501) * t257) * pkin(3) + (-0.4e1 * t145 * t623 + t165 * t616 + t198 * t451) * pkin(1);
t79 = t89 * t170 + t293 * t98;
t77 = 0.1e1 / t79;
t353 = t77 * t357;
t399 = t627 + t411;
t372 = 0.2e1 * t258 * t394 + t399 * t473;
t431 = t227 * t473;
t635 = 0.4e1 * pkin(1);
t636 = -0.2e1 * pkin(1);
t373 = pkin(3) * (t258 * t399 + t431 * t635 + t473 * t636);
t388 = pkin(1) + t391;
t371 = t257 * t373 + t293 * t388 + t372 * t616;
t368 = t374 * t371;
t285 = 0.1e1 / pkin(4);
t593 = t285 / t298;
t620 = -t293 / 0.4e1;
t284 = t303 ^ 2;
t574 = t290 + t292;
t576 = t283 - t303;
t588 = t305 * t303;
t415 = (t576 * t294) + t284 / 0.6e1 + t574 - t588;
t403 = 0.5e1 / 0.6e1 * t300 + t415;
t138 = (t228 + t576) * t298 + t403;
t242 = -0.3e1 / 0.2e1 * t303;
t584 = t284 / 0.2e1 - t300 / 0.2e1;
t472 = -(3 * t588) + t278 + t584;
t587 = t220 * ((t242 + t283) * t294 - 0.3e1 / 0.2e1 * t588 + t574 + t584) + t286;
t105 = -0.6e1 * t138 * t455 + (t260 + ((t267 - 9 * t303) * t294) + t472) * t298 + (t242 + t580) * t300 + t587;
t240 = -t303 / 0.2e1;
t190 = t240 + t298 + t220;
t665 = -0.4e1 * t190;
t434 = t479 * t665;
t423 = pkin(1) * t434;
t122 = t423 + ((t276 + t576) * t298) + t403;
t132 = t190 * t595 + t161;
t563 = 0.6e1 * t623;
t100 = pkin(1) * t122 * t563 + t132 * t540 + t105 + t173;
t279 = 8 * t305;
t446 = t295 * t479;
t147 = -0.4e1 * t446 + t632 + ((t626 + t268 + t279) * t294);
t238 = -t303 / 0.4e1;
t658 = t238 + t298 / 0.2e1;
t151 = -t294 + t422 + t658;
t108 = t501 + t147 * t227 + t222 * t190 + (t151 * t564 + t456 * t594) * pkin(1);
t153 = -t300 / 0.6e1 + t415;
t281 = 4 * t305;
t207 = (t281 + t303) * t294;
t573 = t292 - t290;
t110 = t210 * t423 - t286 + (-t228 + t569) * t300 + (t207 + t300 / 0.6e1 - t284 / 0.6e1 + t573) * t298 + t153 * t305;
t269 = -5 * t303;
t274 = 7 * t290;
t114 = (t242 + t282 + (7 * t294)) * t300 + (t274 + ((t269 + 10 * t305) * t294) + t472) * t298 + t587;
t306 = pkin(7) * t305;
t203 = -0.12e2 * pkin(7) * t295 + t306 * t635;
t209 = -8 * t290 + 12 * t589;
t497 = 0.16e2 * t556;
t225 = t227 ^ 2;
t592 = t290 * t225;
t554 = 0.8e1 * t592;
t126 = t203 * t258 + t209 * t227 + t497 + t554 + t574 - (6 * t589);
t134 = t219 * t190 + t438 * t595;
t628 = -6 * t298;
t174 = 16 * (t305 * t628 + t572) * t290;
t217 = -30 * t303 + 60 * t305;
t223 = t311 ^ 2;
t575 = t284 - t300;
t428 = 6 * t292 + t575 - 6 * t588;
t494 = 0.8e1 * t533;
t541 = 0.32e2 * t598;
t518 = -t298 - t569;
t600 = t213 * (t294 + t518);
t181 = t305 + t298 / 0.4e1 + t294 / 0.4e1 - t303 / 0.8e1;
t581 = 0.4e1 / 0.7e1 * t305 - t303 / 0.7e1;
t116 = -0.32e2 / 0.21e2 * t181 * t455 + 0.5e1 / 0.42e2 * t300 + (0.16e2 / 0.21e2 * t294 + t581) * t298 + t290 / 0.7e1 + t581 * t294 + t292 - 0.3e1 / 0.7e1 * t588 + t284 / 0.42e2;
t183 = t583 + t658;
t251 = 0.4e1 / 0.3e1 * t294;
t117 = -0.8e1 / 0.3e1 * t183 * t455 + 0.5e1 / 0.18e2 * t300 + (t251 + t239) * t298 + t292 - t290 / 0.3e1 + t284 / 0.18e2 + (t245 + 0.2e1 / 0.3e1 * t294 + t241) * t305;
t253 = t294 / 0.2e1;
t582 = t253 + t305;
t155 = -0.2e1 / 0.3e1 * t455 + t238 + t582;
t250 = -0.2e1 / 0.3e1 * t298;
t211 = t250 + t305;
t99 = t211 * t554 + t155 * t497 + 0.14e2 * t116 * t591 - (t571 * t300) + (t207 - 0.10e2 / 0.3e1 * t290 + (2 * t292) - t588) * t298 + t153 * t594 + (t117 * t563 + t212 * t434) * pkin(1);
t76 = t174 * t225 + 0.32e2 * t134 * t556 + 0.24e2 * t110 * t591 + (t268 + t281 + 28 * t294) * t286 + (t193 * t600) + ((t217 * t290) + (t270 * t292) + (t428 * t276) + (t575 * t283) + 0.28e2 * t295 ^ 2 + 0.4e1 * t306 ^ 2 + 0.24e2 * t99 * t311) * t298 + (t100 * t494 + t108 * t541) * t198 + 0.8e1 * (t105 * t623 - t114 * t479) * pkin(1) + (0.16e2 * t126 * t223 + (t217 * t294) + (70 * t290) + t300 + t428) * t300;
t179 = 0.7e1 / 0.6e1 * t298 + t237 + t582;
t247 = t298 / 0.3e1;
t524 = t237 + t247 + t305;
t182 = t251 + t524;
t129 = -t179 * t215 + t182 * t214;
t525 = t303 / 0.3e1 + t247 + t283;
t135 = -(t298 * t571) - 0.5e1 / 0.3e1 * t290 + t525 * t294 + t305 * (t239 + t210);
t229 = -0.20e2 / 0.3e1 * t294;
t526 = 0.2e1 / 0.3e1 * t303 + t246 + t281;
t527 = 0.4e1 / 0.3e1 * t303 + t245 - (2 * t305);
t136 = -t300 + (t229 + t526) * t298 - (3 * t290) + t527 * t294 + t292;
t187 = t294 + t524;
t206 = -t294 + t567;
t140 = t187 * t214 - t206 * t215 / 0.2e1;
t515 = t226 * t616;
t467 = t290 * t515;
t101 = -0.4e1 * pkin(7) * t467 + t129 * t552 + (-0.8e1 / 0.3e1 * t592 + t135) * t214 + (t140 * t564 - t616 * t136 / 0.2e1) * pkin(1);
t578 = t268 - 2 * t298;
t520 = t280 + t578;
t128 = t300 + (t241 + t250 + t579) * t298 + t275 + (t520 * t294) + t305 * (t241 + t211);
t125 = t128 * t214;
t139 = t271 + ((t262 + t520) * t298) + (t250 + t521) * t220;
t510 = t616 * t139;
t111 = -pkin(1) * t510 + t125;
t137 = -(3 * t300) + (t229 + t527) * t298 + t526 * t294 + t573;
t141 = -0.5e1 / 0.3e1 * t300 + (-t294 + t525) * t298 + t305 * (t249 + t523);
t495 = -0.2e1 * t215;
t112 = t137 * t214 + t141 * t495;
t657 = t282 - t298 - t303;
t192 = t657 * t262;
t577 = t269 - 5 * t298;
t113 = t286 + ((21 * t294 + t657) * t300) + ((t305 * t578 + t192 + t278 + 35 * t290) * t298) + ((t274 + (t279 + t577) * t294 + t305 * t518) * t220);
t169 = 0.4e1 * t412;
t175 = 0.2e1 * t215 + t214;
t189 = t240 + t204;
t115 = t222 * t214 + t169 * t227 + (t175 * t623 + t189 * t616) * t562;
t120 = 0.7e1 * t286 + ((35 * t294 + 15 * t305 + t577) * t300) + ((21 * t290 + t192 + 9 * t292 + (t270 + t628) * t305) * t298) + t600;
t180 = t305 + 0.5e1 / 0.2e1 * t298 + 0.3e1 / 0.2e1 * t294 + t240;
t142 = t180 * t214 + t219 * t215 / 0.2e1;
t157 = 0.4e1 / 0.3e1 * t591 + t201 + t212;
t599 = t223 * t300;
t470 = t616 * t599;
t445 = -0.24e2 * t470;
t529 = t198 * t598;
t493 = -0.8e1 * t529;
t539 = -0.12e2 * t590;
t184 = 0.8e1 / 0.3e1 * t298 + t523;
t186 = t239 + t246 + t568;
t133 = -t184 * t215 + t186 * t214;
t191 = 0.5e1 / 0.6e1 * t298 + t253 + t237;
t144 = pkin(1) * t468 + t191 * t558;
t544 = t226 * t214;
t488 = t295 * t544;
t102 = -0.8e1 * pkin(7) * t488 + t144 * t553 + t125 + (t133 * t564 - t510) * pkin(1);
t669 = -0.6e1 * t102;
t80 = t142 * t497 + t115 * t493 + t101 * t539 - 0.6e1 * t112 * t591 + (t516 * t669 + (0.24e2 * t210 * t592 - t113) * t257) * pkin(3) + (-0.6e1 * t111 * t623 + t120 * t616 + t157 * t445) * pkin(1);
t64 = t76 * t170 + t293 * t80;
t50 = (t353 * t64 + t368 * t620) * t593;
t364 = t77 * t368 / 0.4e1;
t51 = (t293 * t357 + t364 * t64) * t593;
t25 = -t162 * t50 + t163 * t51;
t26 = t162 * t51 + t163 * t50;
t465 = t298 * t615 * t257;
t449 = 0.4e1 * t465;
t426 = t449 + t199;
t410 = t227 * t426;
t442 = pkin(3) * t461;
t427 = pkin(1) * t442;
t547 = pkin(7) * t616;
t489 = t298 * t547;
t631 = -4 * t294;
t634 = t257 ^ 2;
t103 = -0.4e1 * t379 * t543 - 0.4e1 * t615 * t390 + t410 * t631 + (0.2e1 * t615 * t440 + t435) * t257 + (-0.4e1 * (0.2e1 * t427 + t199) * t430 + 0.8e1 * t634 * t489) * pkin(1);
t104 = 0.1e1 / t293;
t622 = t104 / 0.2e1;
t380 = t388 * t622;
t396 = pkin(3) * t162;
t340 = t374 * (t293 * t396 + t103 * t380 + (t426 * t560 + (t294 + (0.4e1 * t533 + 0.2e1 * pkin(7)) * pkin(7) + t517) * t214) * t616 - t615 * t373 + t298 * t634 * (t565 + (0.4e1 * t227 - 0.2e1) * pkin(1)));
t330 = t340 * t593;
t375 = 0.1e1 / t376 ^ 2;
t369 = t375 * t371;
t387 = t396 * t562 + t199;
t348 = t387 * t369;
t344 = t348 * t593;
t354 = t360 * t593;
t351 = t104 * t354 / 0.8e1;
t365 = t368 * t593;
t78 = 0.1e1 / t79 ^ 2;
t663 = t64 * t78;
t649 = -t663 / 0.4e1;
t359 = t365 * t649;
t363 = t364 * t593;
t420 = t227 * t295 * t461;
t614 = pkin(3) * t210;
t402 = t420 * t614;
t154 = 0.24e2 * t402;
t421 = t258 * t294 * t461;
t413 = pkin(7) * t421;
t407 = pkin(3) * t413;
t167 = 0.4e1 * t407;
t483 = pkin(1) * t507;
t625 = pkin(3) * pkin(7);
t439 = t483 * t625;
t177 = -0.2e1 * t439;
t405 = 0.24e2 * t413;
t416 = 0.32e2 / 0.3e1 * t467;
t443 = pkin(1) * t461;
t425 = t212 * t443;
t464 = t615 * t597;
t436 = t464 * t625;
t601 = t198 * t257;
t475 = t299 * t311 * t601;
t447 = -0.24e2 * t475;
t514 = t227 * t615;
t463 = t294 * t514;
t448 = -0.4e1 * t463;
t459 = -0.4e1 * t483;
t462 = t615 * t592;
t480 = t298 * t516;
t481 = pkin(3) * t516;
t506 = t103 * t622;
t528 = t226 * t595;
t545 = pkin(3) * t601;
t548 = 0.3e1 * t590;
t557 = 0.8e1 * t170;
t648 = 0.64e2 / 0.3e1 * t181;
t662 = t224 * t300;
t666 = 0.6e1 * t138;
t38 = (-0.96e2 * t157 * t482 * t662 + t177 * t493 + t115 * t447 - 0.24e2 * t101 * t465 - 0.6e1 * (0.8e1 * t191 * t463 - t615 * t128 + (t186 * t459 + 0.8e1 * t464) * pkin(7)) * t480 + t545 * t669 - 0.24e2 * t462 * t614 - 0.16e2 * t180 * t436 + 0.6e1 * t128 * t439 + t113 * t533 + ((-t222 * t615 + t448) * t493 + (0.8e1 / 0.3e1 * t462 + t182 * t448 + t187 * pkin(7) * t459 - t135 * t615) * t539 + 0.6e1 * t137 * t463) * pkin(3)) * t293 + t80 * t506 + (0.8e1 * t126 * t257 * t662 + 0.4e1 * (t167 + (0.2e1 * t443 * t594 + 0.4e1 * t420) * pkin(3)) * t529 + 0.12e2 * t108 * t475 + (t420 * t648 + 0.4e1 * t190 * t425 + (0.16e2 * t183 * t421 + t416 * t615) * pkin(7)) * pkin(3) * t548 + 0.6e1 * t99 * t465 + (t154 + (t190 * t405 + t443 * t666) * pkin(3)) * t481 + t100 * t545 + 0.8e1 * t290 * pkin(7) * t442 * t528 + 0.12e2 * t190 * t402 + t407 * t666 + t114 * t427) * t557 + t76 * t214;
t613 = t64 * t77;
t382 = t392 * t622;
t419 = 0.4e1 * t427;
t320 = t374 * (t163 * pkin(3) * t293 + t103 * t382 + t410 * t636 - t381 * t543 - (t419 + t199) * t430 + pkin(1) * t449 - 0.2e1 * t311 * t489 + t257 * t406 - t615 * t408);
t361 = t375 * t370;
t338 = t387 * t361;
t503 = t593 / 0.4e1;
t504 = -t593 / 0.4e1;
t656 = t320 * t503 + t338 * t504;
t458 = -0.4e1 * t481;
t73 = (t447 * t215 + (t177 + (t571 * t615 - 0.2e1 * t463) * pkin(3)) * t550 - 0.8e1 * t121 * t465 + (t177 + (-t204 * t615 + 0.2e1 * t463) * pkin(3)) * t458 + t545 * t667 - 0.8e1 * t436 - t197 * t514 + 0.4e1 * t204 * t439 + t158 * t533) * t293 + t98 * t506 + (0.24e2 * t160 * t475 + (t167 + (0.8e1 / 0.3e1 * t420 + 0.2e1 * t425) * pkin(3)) * t538 + 0.24e2 * t107 * t465 + 0.48e2 * t413 * t480 + t545 * t668 + t154 + 0.6e1 * t131 * t427 + (pkin(3) * t405 + 0.24e2 * t443 * t480) * t185) * t170 + t89 * t214;
t610 = t103 * t351 + t73 * t359 + t38 * t363 + (t330 / 0.4e1 - t344 / 0.4e1) * t613 + t656 * t293 - t50;
t349 = t354 * t649;
t352 = t353 * t593;
t362 = -t104 * t365 / 0.8e1;
t619 = t293 / 0.4e1;
t611 = t103 * t362 + t330 * t620 + t344 * t619 + t349 * t73 + t352 * t38 + t656 * t613 + t51;
t3 = -t162 * t611 + t163 * t610;
t358 = t637 * t360;
t355 = -t358 / 0.2e1;
t367 = t637 * t368;
t366 = t367 / 0.2e1;
t86 = t258 * t355 + t366 * t616;
t432 = t285 * t86 * t608 * t663;
t536 = t77 * t606;
t484 = t285 * t536;
t453 = t86 * t484;
t454 = t64 * t484;
t532 = t104 * t607;
t477 = t103 * t532;
t531 = t293 * t608;
t337 = t637 * t340;
t346 = t637 * t348;
t621 = -t258 / 0.2e1;
t679 = t320 * t606 + t338 * t608;
t70 = t337 * t621 + t258 * t346 / 0.2e1 - t679 * t616;
t650 = t337 / 0.2e1 - t346 / 0.2e1;
t71 = -t258 * t679 + t650 * t616;
t85 = t355 * t616 + t367 * t621;
t383 = t71 * t454 + t38 * t453 + t73 * t432 + (t477 * t85 + t531 * t70) * t285;
t4 = t162 * t610 + t163 * t611;
t429 = t285 * t531 * t85 + t64 * t453;
t530 = t293 * t606;
t534 = t85 * t606;
t659 = (t530 * t86 + t534 * t613) * t285;
t535 = t85 * t608;
t492 = t77 * t535;
t680 = (-t38 * t492 + (t535 * t73 * t78 + t536 * t70) * t64 - t477 * t86 - t531 * t71) * t285;
t685 = t25 * t680 + t26 * t383 + t3 * t659 + t4 * t429;
t389 = t397 * t616;
t409 = t431 * t214;
t106 = t409 * t629 + 0.4e1 * t616 * t377 + (t389 * t629 + t395 * t633) * t258;
t566 = t616 ^ 2;
t323 = t374 * (-t661 + t106 * t380 + t566 * t393 + t372 * t258 + (-0.8e1 * t216 * t414 - t399 * t616) * t214);
t321 = t323 * t593;
t660 = t392 * t562;
t347 = t369 * t660;
t342 = t347 * t593;
t460 = pkin(7) * t487;
t178 = -0.4e1 * t460;
t200 = pkin(7) * t495;
t444 = -0.32e2 * t467;
t511 = t294 * t616;
t466 = t258 * t511;
t450 = -0.8e1 * t466;
t452 = t210 * t488;
t490 = t227 * t547;
t457 = -0.24e2 * t490;
t474 = t295 * t528;
t478 = -0.48e2 * t508;
t512 = t258 * t616;
t485 = -0.2e1 * t512;
t486 = t212 * t543;
t499 = -0.4e1 * t547;
t502 = pkin(7) * t552;
t505 = t106 * t622;
t537 = 0.24e2 * t590;
t555 = -0.4e1 * t597;
t433 = t295 * t457;
t586 = t219 * t433 - 0.24e2 * t452;
t624 = pkin(7) * t227;
t36 = 0.12e2 * (-(t129 * t450 + t179 * t555 + t416 * t214) * t590 - 0.8e1 * t210 * t467 * t214 + t141 * t597 + t112 * t466) * t293 + t80 * t505 + (0.16e2 * (-t203 * t616 + t209 * t485 + t444) * t599 + (-0.28e2 * t116 * t466 + t211 * t444 - t488 * t648) * t537 - 0.4e1 * t174 * t515 - 0.96e2 * t190 * t452 - 0.48e2 * t110 * t466) * t170 + ((((8 * t219) + 0.48e2 * t590) * t293 + (-0.256e3 * t590 - 0.64e2 * t595) * t170 * t214) * t592 + ((((t206 * t664) + 0.12e2 * t290 * t566) * t539 + t142 * t478 + t139 * t276) * t293 + (-0.768e3 * t295 * t470 + (t155 * t478 - 0.16e2 * t183 * t542) * t537 - 0.96e2 * t134 * t508 - 0.48e2 * t138 * t542) * t170) * t227) * pkin(7) + (((-0.8e1 / 0.3e1 * t466 + t200) * t445 + 0.8e1 * t115 * t469 + (t136 * t621 + t140 * t499) * t539 + 0.6e1 * t102 * t442 + 0.6e1 * t111 * t547 + (-0.24e2 * t157 * t599 + t120) * t258) * t293 + (-0.4e1 * t108 * t469 + (-0.6e1 * t117 * t547 + t486 * t665) * t548 - t100 * t442 - t105 * t547 - t114 * t543) * t557 - t76 * t616) * pkin(1) + ((-0.8e1 * (t169 * t485 + t502 + t555 + (-t175 * t547 + t189 * t258) * t562) * t598 - 0.6e1 * (-0.4e1 * t474 + (-pkin(1) * t139 + 0.8e1 * t144 * t511) * t258 + (-0.4e1 * t133 * t215 + (t184 * t631 + 0.24e2 * t446) * t227) * pkin(7)) * t533) * t293 + ((t147 * t485 + t178 + (t457 - 0.4e1 * t544) * t295 + (t151 * t499 - 0.2e1 * t543 * t594) * pkin(1)) * t541 + (0.24e2 * (-t190 * t214 * t624 - t132 * t512) * t294 + (-t122 * t547 - t138 * t543) * t561 + t586) * t494) * t170) * t198;
t496 = pkin(7) * t215;
t65 = (t566 * t224 * t551 + (t168 * t485 - 0.2e1 * t597 + (pkin(7) * t456 + t218 * t258) * pkin(1)) * t550 + t118 * t419 + ((0.4e1 * t258 * t479 - 0.2e1 * t624) * t294 + (-0.2e1 * t176 * t547 - t205 * t258) * pkin(1)) * t458 + t433 * t214 + 0.4e1 * t474 + t146 * t485 + t205 * t502 + 0.4e1 * t145 * t496 + (t493 + t165) * t216) * t293 + t98 * t505 + (0.8e1 * (t450 + t200) * t529 + (-0.12e2 * t149 * t466 + t178 + (-0.12e2 * t490 - 0.8e1 / 0.3e1 * t544) * t295) * t538 - 0.6e1 * t109 * t427 + 0.6e1 * (-0.8e1 * t460 + t202 * t485 + (-t152 * t547 - t185 * t543) * t635) * t481 - 0.24e2 * t130 * t466 - 0.24e2 * t185 * t460 - 0.6e1 * t119 * t496 + t586 + (t160 * t451 + (t156 * t499 - 0.2e1 * t486) * t538 - 0.6e1 * t131 * t543) * pkin(1)) * t170 - t89 * t215;
t336 = t361 * t660;
t341 = t374 * (t106 * t382 + t293 * t391 + t378 * t616 + t389 * t559 + t409 * t635 + t437 * t596 + t465 * t565);
t653 = -t336 * t504 + t341 * t503;
t612 = t106 * t351 + t65 * t359 + t36 * t363 + (t321 / 0.4e1 + t342 / 0.4e1) * t613 + t653 * t293 + t50;
t618 = -t106 * t362 - t321 * t620 + t342 * t619 - t349 * t65 - t352 * t36 - t653 * t613 + t51;
t1 = t162 * t618 + t163 * t612;
t2 = t162 * t612 - t163 * t618;
t476 = t106 * t532;
t651 = (t336 + t341) * t606;
t654 = t323 * t606 - t347 * t608;
t60 = -t258 * t651 + t654 * t616 - t85;
t61 = t654 * t258 + t651 * t616 - t86;
t384 = t60 * t454 + t36 * t453 + t65 * t432 + (t476 * t85 + t530 * t61) * t285;
t681 = (t36 * t492 + (t534 * t65 * t78 + t536 * t61) * t64 + t476 * t86 + t531 * t60) * t285;
t684 = t1 * t659 + t2 * t429 - t25 * t681 + t26 * t384;
t683 = -t25 * t383 + t26 * t680 - t3 * t429 + t4 * t659;
t682 = t1 * t429 - t2 * t659 + t25 * t384 + t26 * t681;
t259 = cos(pkin(9));
t609 = cos(pkin(8));
t195 = t609 * t616;
t255 = sin(pkin(8));
t491 = t255 * t258 - t195;
t350 = t491 * t358;
t404 = t255 * t616 + t258 * t609;
t314 = t350 / 0.2e1 + t654 * t491 + (-t366 + t651) * t404;
t356 = t491 * t366;
t55 = t356 - t651 * t491 + (-t355 + t654) * t404;
t617 = sin(pkin(9));
t33 = t259 * t55 - t314 * t617;
t34 = t259 * t314 + t55 * t617;
t316 = t404 * t358 / 0.2e1 + t356;
t84 = -t350 / 0.2e1 + t404 * t366;
t74 = t259 * t84 - t316 * t617;
t75 = t259 * t316 + t617 * t84;
t20 = t33 * t86 - t34 * t85 + t60 * t74 + t61 * t75;
t418 = -t33 * t85 - t86 * t34 - t60 * t75 + t61 * t74;
t46 = t74 * t86 - t75 * t85;
t500 = -t74 * t85 - t86 * t75;
t674 = -t20 * t74 - t33 * t46 + t34 * t500 + t418 * t75;
t673 = t20 * t75 + t33 * t500 + t34 * t46 + t418 * t74;
t315 = t679 * t404 + t650 * t491;
t63 = t650 * t404 - t491 * t679;
t43 = t259 * t63 - t315 * t617;
t44 = t259 * t315 + t617 * t63;
t22 = t43 * t86 - t44 * t85 - t70 * t75 + t71 * t74;
t417 = -t43 * t85 - t86 * t44 - t70 * t74 - t71 * t75;
t672 = -t22 * t74 + t417 * t75 - t43 * t46 + t44 * t500;
t671 = t22 * t75 + t417 * t74 + t43 * t500 + t44 * t46;
t605 = t60 * pkin(3) + t215;
t604 = t60 * pkin(2) + t215;
t603 = -t61 * pkin(3) - t216;
t602 = -t61 * pkin(2) - t216;
t124 = t255 ^ 2 * t616 + t195 * t609;
t123 = -t255 * t491 - t404 * t609;
t69 = t71 * pkin(2);
t68 = t71 * pkin(3);
t67 = t70 * pkin(2);
t66 = t70 * pkin(3);
t5 = [-m(2) * (g(1) * (rSges(2,1) * t616 + t258 * rSges(2,2)) + g(2) * (-t258 * rSges(2,1) + rSges(2,2) * t616)) - m(3) * (g(1) * (rSges(3,1) * t60 + rSges(3,2) * t61 + t215) + g(2) * (-rSges(3,1) * t61 + rSges(3,2) * t60 - t216)) - m(4) * (g(1) * (rSges(4,1) * t418 - rSges(4,2) * t20 + t604) + g(2) * (rSges(4,1) * t20 + rSges(4,2) * t418 + t602)) - m(5) * (g(1) * (rSges(5,1) * t384 + rSges(5,2) * t681 + t605) + g(2) * (-rSges(5,1) * t681 + rSges(5,2) * t384 + t603)) - m(6) * (g(1) * (rSges(6,1) * t124 - rSges(6,2) * t123) + g(2) * (rSges(6,1) * t123 + rSges(6,2) * t124)) - m(7) * (g(1) * (rSges(7,1) * t418 - rSges(7,2) * t20 + t215) + g(2) * (rSges(7,1) * t20 + rSges(7,2) * t418 - t216)) - m(9) * (g(1) * t215 - g(2) * t216) - m(10) * (g(1) * (t418 * pkin(6) + rSges(10,1) * t674 - rSges(10,2) * t673 + t604) + g(2) * (t20 * pkin(6) + rSges(10,1) * t673 + rSges(10,2) * t674 + t602)) - m(11) * (g(1) * (t384 * pkin(4) + t684 * rSges(11,1) + t682 * rSges(11,2) + t605) + g(2) * (-t681 * pkin(4) - t682 * rSges(11,1) + t684 * rSges(11,2) + t603)), -m(3) * (g(1) * (rSges(3,1) * t71 - rSges(3,2) * t70) + g(2) * (rSges(3,1) * t70 + rSges(3,2) * t71)) - m(4) * (g(1) * (rSges(4,1) * t417 - rSges(4,2) * t22 + t69) + g(2) * (rSges(4,1) * t22 + rSges(4,2) * t417 + t67)) - m(5) * (g(1) * (rSges(5,1) * t383 - rSges(5,2) * t680 + t68) + g(2) * (rSges(5,1) * t680 + rSges(5,2) * t383 + t66)) - m(7) * (g(1) * (rSges(7,1) * t417 - rSges(7,2) * t22) + g(2) * (rSges(7,1) * t22 + rSges(7,2) * t417)) - m(8) * (g(1) * (t257 * rSges(8,1) - rSges(8,2) * t615) + g(2) * (rSges(8,1) * t615 + t257 * rSges(8,2))) - m(10) * (g(1) * (t417 * pkin(6) + rSges(10,1) * t672 - rSges(10,2) * t671 + t69) + g(2) * (t22 * pkin(6) + rSges(10,1) * t671 + rSges(10,2) * t672 + t67)) - m(11) * (g(1) * (t383 * pkin(4) + t685 * rSges(11,1) - t683 * rSges(11,2) + t68) + g(2) * (t680 * pkin(4) + t683 * rSges(11,1) + t685 * rSges(11,2) + t66))];
taug = t5(:);
