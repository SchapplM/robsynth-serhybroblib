% Calculate joint inertia matrix for
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1TE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1TE_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1TE_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1TE_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:14:20
% EndTime: 2020-04-27 09:14:37
% DurationCPUTime: 17.38s
% Computational Cost: add. (279048->649), mult. (833036->1035), div. (2064->14), fcn. (102198->6), ass. (0->459)
t261 = pkin(5) ^ 2;
t260 = t261 ^ 2;
t264 = pkin(4) ^ 2;
t263 = t264 ^ 2;
t437 = t260 - t263;
t529 = 4 * pkin(3);
t265 = 0.1e1 / pkin(4);
t231 = cos(qJ(1));
t490 = pkin(2) * t231;
t163 = pkin(1) - t490;
t230 = cos(qJ(2));
t147 = t163 * t230;
t269 = pkin(3) ^ 2;
t249 = 3 * t269;
t273 = pkin(2) ^ 2;
t275 = pkin(1) ^ 2;
t430 = t275 - t261;
t343 = t273 + t430;
t148 = t249 + t264 + t343;
t228 = sin(qJ(2));
t229 = sin(qJ(1));
t468 = t228 * t229;
t158 = pkin(2) * t468;
t150 = t158 - pkin(3);
t344 = t269 + t430;
t156 = t264 + t344;
t252 = 0.3e1 * t273;
t196 = t231 ^ 2;
t470 = t196 * t273;
t405 = 0.2e1 * t470;
t512 = -0.2e1 * t231;
t424 = pkin(1) * t512;
t413 = pkin(1) * t490;
t168 = -0.2e1 * t413;
t434 = -t273 + t275;
t115 = t168 + t405 + t434;
t193 = t230 ^ 2;
t479 = t115 * t193;
t466 = t230 * t229;
t391 = pkin(2) * t466;
t105 = t163 * t228 + t391;
t313 = t269 + t343;
t149 = -t264 + t313;
t112 = t168 + t149;
t510 = 4 * t269;
t197 = t273 * t510;
t456 = t273 * t275;
t169 = t197 - 0.4e1 * t456;
t243 = 0.2e1 * t264;
t259 = -0.2e1 * t275;
t267 = t269 ^ 2;
t502 = -pkin(4) + pkin(5);
t503 = -pkin(4) - pkin(5);
t335 = pkin(3) * t158;
t151 = -0.2e1 * t335;
t92 = t151 + t112;
t276 = sqrt(t169 * t196 + 0.4e1 * t149 * t413 - t267 - (t275 + (pkin(2) - t503) * (pkin(2) + t503)) * (t275 + (pkin(2) - t502) * (pkin(2) + t502)) + (t243 + t259 + 0.2e1 * t261 - 0.6e1 * t273 - 0.4e1 * t479) * t269 + (t112 * t158 - t147 * t92) * t529);
t480 = t105 * t276;
t318 = -0.4e1 * t335;
t93 = t148 + t168 + t318;
t37 = -t93 * t147 + t480 + (t148 * t468 + t150 * t424) * pkin(2) + (-t156 - t252 + t405 - 0.2e1 * t479) * pkin(3);
t294 = t37 ^ 2;
t36 = 0.1e1 / t294;
t448 = -t158 + t147;
t103 = pkin(3) + t448;
t184 = pkin(2) * t229;
t111 = t168 + t264 + t313;
t183 = pkin(3) * t230;
t387 = t115 * t183;
t67 = t111 * t163 + 0.2e1 * t387;
t489 = pkin(3) * t163;
t68 = t111 * t230 + (0.4e1 * t193 - 0.2e1) * t489;
t39 = t103 * t276 + t184 * t68 + t228 * t67;
t38 = t39 ^ 2;
t433 = t273 + t275;
t345 = t269 + t433;
t386 = pkin(3) * t147;
t78 = t151 + t168 + t345 + 0.2e1 * t386;
t331 = pkin(4) * t265 / (t36 * t38 + 0.1e1) * t78;
t304 = t36 * t39 * t331;
t190 = t228 ^ 2;
t47 = 0.1e1 / t276;
t506 = t47 / 0.2e1;
t358 = t103 * t506;
t464 = t230 * t269;
t402 = 0.8e1 * t464;
t414 = 0.4e1 * t183;
t417 = 0.4e1 * t489;
t457 = t269 * t193;
t49 = (t115 * t402 + t417 * t92) * t228 + (t112 * t414 + 0.8e1 * t163 * t457) * t184;
t516 = -2 * pkin(3);
t31 = -t480 + t49 * t358 + t115 * t190 * t516 + t67 * t230 + (-t111 - 0.8e1 * t386) * t158;
t310 = 0.1e1 / t37 * t331;
t357 = t105 * t506;
t463 = t230 * t273;
t32 = t448 * t276 + t49 * t357 + (t163 * t93 + 0.4e1 * t387) * t228 + (t424 * t463 + (t148 * t230 + t193 * t417) * pkin(2)) * t229;
t76 = 0.1e1 / t78 ^ 2;
t504 = -t76 / 0.2e1;
t425 = 2 * pkin(3);
t98 = t105 * t425;
t364 = t98 * t504;
t75 = 0.1e1 / t78;
t505 = t75 / 0.2e1;
t17 = 0.1e1 + 0.2e1 * (t31 * t505 - t364 * t39) * t310 - 0.2e1 * (t32 * t505 - t364 * t37) * t304;
t530 = t17 * Ifges(5,3);
t519 = 0.2e1 * pkin(1);
t199 = 0.10e2 / 0.3e1 * t269;
t211 = -0.2e1 / 0.3e1 * t261;
t216 = 0.2e1 / 0.3e1 * t264;
t258 = 0.2e1 * t275;
t348 = t211 + t216 + t258;
t186 = t269 + t275;
t281 = t273 ^ 2;
t347 = t211 + t186;
t450 = (t216 + t347) * t186 + t281;
t80 = (t199 + t348) * t273 + t450;
t528 = -0.6e1 * t80;
t274 = t275 ^ 2;
t435 = t267 + t274;
t438 = t258 - t261;
t455 = t275 * t261;
t520 = -t260 / 0.6e1 + t263 / 0.6e1;
t106 = t438 * t269 + t435 - t455 - t520;
t295 = t106 + t281;
t82 = (t199 + t438) * t273 + t295;
t527 = -0.6e1 * t82;
t423 = 0.2e1 * t184;
t248 = 6 * t269;
t209 = -t261 / 0.3e1;
t217 = t264 / 0.3e1;
t132 = t209 + t217 + t345;
t526 = -0.24e2 * t132;
t238 = 10 * t269;
t412 = pkin(1) * t183;
t166 = 0.2e1 * t412;
t178 = -t269 / 0.3e1 + t275;
t110 = 0.4e1 / 0.3e1 * t457 + t166 + t178;
t257 = 0.3e1 * t275;
t521 = t257 - t261 - t264;
t155 = t521 * t238;
t242 = -0.6e1 * t261;
t280 = pkin(2) * t273;
t270 = t280 ^ 2;
t241 = -0.5e1 * t261;
t439 = t241 - 0.5e1 * t264;
t194 = t196 ^ 2;
t473 = t194 * t281;
t181 = t186 ^ 2;
t476 = t181 * (-t264 + t344);
t525 = 0.7e1 * t270 + ((35 * t269) + 0.15e2 * t275 + t439) * t281 + ((21 * t267) + t155 + 0.9e1 * t274 + (t242 - 0.6e1 * t264) * t275) * t273 + t476 - 0.24e2 * t110 * t473;
t524 = t228 * t273;
t188 = -(3 * t269) + t275;
t398 = 0.4e1 * t457;
t523 = t188 + t398;
t522 = t209 - t264 / 0.3e1;
t518 = -0.4e1 * pkin(2);
t517 = -4 * pkin(3);
t189 = -0.3e1 * t273 + t275;
t192 = t230 * t193;
t284 = pkin(3) * t269;
t474 = t192 * t284;
t396 = pkin(1) * t474;
t341 = 0.8e1 * t396;
t137 = t189 * t341;
t355 = 0.6e1 * t412;
t377 = 0.12e2 * t457;
t212 = -0.3e1 / 0.2e1 * t261;
t236 = 15 * t267;
t239 = 0.18e2 * t275;
t253 = 0.3e1 * t274;
t447 = t260 / 0.2e1 - t263 / 0.2e1;
t316 = -0.3e1 * t455 + t253 + t447;
t442 = (15 * t269) + t257;
t481 = t270 + t186 * ((t212 + t258) * t269 - 0.3e1 / 0.2e1 * t455 + t435 + t447);
t50 = t335 * t527 + (t236 + (t239 - 0.9e1 * t261) * t269 + t316) * t273 + (t212 + t442) * t281 + t481;
t210 = -t261 / 0.2e1;
t154 = t210 + t345;
t305 = t154 * t318;
t65 = t305 + (t248 + t438) * t273 + t295;
t179 = t275 - t273 / 0.3e1;
t117 = t179 * t151;
t460 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t79 = t154 * t460 + t117;
t42 = t355 * t65 + t377 * t79 + t137 + t50;
t515 = 0.8e1 * t42;
t133 = t163 + t183;
t514 = -0.8e1 * t133;
t513 = -0.2e1 * t230;
t511 = 4 * t267;
t509 = pkin(1) * pkin(3);
t114 = t166 + t523;
t240 = -0.2e1 * t261;
t164 = pkin(1) + t183;
t195 = t231 * t196;
t472 = t195 * t280;
t363 = t164 * t472;
t385 = pkin(3) * t468;
t478 = t164 * t231;
t207 = -t261 / 0.6e1;
t226 = t273 / 0.2e1;
t444 = t226 + t275;
t300 = -t335 + t444;
t104 = t264 / 0.6e1 + t207 + t300;
t426 = 0.4e1 * pkin(1);
t461 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t223 = t269 / 0.3e1;
t96 = -0.4e1 / 0.9e1 * t335 + t275 + t273 / 0.3e1 + t223 + t264 / 0.9e1 - t261 / 0.9e1;
t51 = t178 * t151 + t132 * t461 + 0.6e1 * t96 * t457 + (t104 * t183 + t474) * t426;
t102 = t151 + t132;
t170 = t434 * t510;
t356 = 0.4e1 * t412;
t62 = t132 * t318 + (t248 + t348) * t273 + t450;
t52 = t102 * t356 + t170 * t193 + t62;
t77 = t132 * t460 + t117;
t34 = -0.8e1 * t114 * t363 + t137 + t77 * t377 + t62 * t355 + t270 + (-t261 + t264 + t442) * t281 + t181 * t156 + (0.12e2 * t51 * t196 + t236 + (t239 + t242 + 0.6e1 * t264) * t269 + t253 + (t240 + t243) * t275) * t273 + 0.6e1 * (-t385 * t80 - t478 * t52) * pkin(2);
t247 = 5 * t267;
t431 = t274 + t281;
t441 = t238 + t258;
t454 = t275 * t269;
t113 = t441 * t273 + t247 + t431 + 0.6e1 * t454;
t251 = 0.5e1 * t281;
t255 = 0.6e1 * t275;
t122 = t251 + (t238 + t255) * t273 + t181;
t326 = 0.8e1 * t363;
t406 = -0.4e1 * t470;
t493 = pkin(1) * t230;
t174 = t249 + t433;
t138 = t174 * t184;
t486 = pkin(3) * t228;
t145 = t184 - t486;
t390 = t269 * t184;
t325 = t193 * t390;
t172 = t252 + t186;
t477 = t172 * t228;
t63 = -0.2e1 * t325 + t138 + (0.2e1 * t145 * t493 - t477) * pkin(3);
t301 = -t228 * t284 + t390;
t135 = 0.2e1 * t301;
t185 = (2 * t269) + t273;
t187 = -t269 + t275;
t66 = t185 * t486 + t135 * t193 + (t187 + t166) * t184;
t95 = -pkin(3) * t477 + t138;
t171 = pkin(2) * t511 + 0.8e1 * t269 * t280;
t360 = t284 * t460;
t97 = t171 * t229 + 0.4e1 * t228 * t360;
t40 = t66 * t406 + t193 * t97 + (-0.4e1 * t95 * t493 + (t122 + t326) * t228) * pkin(3) + (0.4e1 * t63 * t478 + (-t113 + t341) * t229) * pkin(2);
t27 = t133 * t34 + t276 * t40;
t508 = 0.1e1 / t27 / 0.4e1;
t507 = -0.1e1 / t27 ^ 2 / 0.4e1;
t465 = t230 * t231;
t120 = -t465 - t468;
t467 = t228 * t231;
t121 = -t466 + t467;
t139 = 0.16e2 * (t431 - 0.6e1 * t456) * t267;
t182 = -0.30e2 * t261 + 0.60e2 * t275;
t191 = t193 ^ 2;
t256 = 0.4e1 * t275;
t277 = pkin(1) * t275;
t302 = 0.6e1 * t274 + t437 - 0.6e1 * t455;
t321 = -0.32e2 * t363;
t208 = -t261 / 0.4e1;
t224 = t269 / 0.2e1;
t109 = -0.2e1 / 0.3e1 * t335 + t275 + t224 + t208;
t175 = (t256 + t261) * t269;
t180 = t275 - 0.2e1 / 0.3e1 * t273;
t383 = 0.16e2 * t474;
t458 = t267 * t191;
t401 = 0.8e1 * t458;
t143 = t275 + t273 / 0.4e1 + t269 / 0.4e1 - t261 / 0.8e1;
t443 = 0.4e1 / 0.7e1 * t275 - t261 / 0.7e1;
t60 = -0.32e2 / 0.21e2 * t143 * t335 + t281 / 0.7e1 + (0.16e2 / 0.21e2 * t269 + t443) * t273 + t267 / 0.7e1 + t443 * t269 + t274 - 0.3e1 / 0.7e1 * t455 + t437 / 0.42e2;
t144 = t208 + t223 + t444;
t222 = 0.4e1 / 0.3e1 * t269;
t499 = 0.4e1 / 0.3e1 * t273;
t61 = -0.8e1 / 0.3e1 * t144 * t335 + t281 / 0.3e1 + (t222 + t209) * t273 + t274 - t267 / 0.3e1 + (t499 + 0.2e1 / 0.3e1 * t269 + t211) * t275 + t437 / 0.18e2;
t41 = t180 * t401 + 0.14e2 * t60 * t457 + t178 * t305 + t187 * t281 + (t175 - 0.10e2 / 0.3e1 * t267 + 0.2e1 * t274 - t455) * t273 + t106 * t461 + (t109 * t383 + 0.6e1 * t183 * t61) * pkin(1);
t254 = 0.8e1 * t275;
t101 = -0.4e1 * t284 * t158 + t197 + t511 + (t240 + t254) * t269;
t107 = t208 - t269 + t300;
t409 = 0.8e1 * t474;
t53 = t151 * t461 + t101 * t193 + t154 * t188 + (t107 * t414 + t409) * pkin(1);
t432 = t274 - t267;
t54 = t179 * t305 - t270 + (-t199 - t430) * t281 + (t175 + t432 + t520) * t273 + t106 * t275;
t246 = 7 * t267;
t59 = (t212 + t257 + (7 * t269)) * t281 + (t246 + (t241 + 0.10e2 * t275) * t269 + t316) * t273 + t481;
t162 = -0.12e2 * pkin(1) * t284 + t277 * t529;
t177 = -(8 * t267) + 0.12e2 * t454;
t328 = pkin(1) * t383;
t71 = t162 * t230 + t177 * t193 + t328 + t401 + t435 - 0.6e1 * t454;
t83 = t151 * t460 + t154 * t189;
t24 = t53 * t321 + t139 * t191 + 0.24e2 * t54 * t457 + (t240 + t256 + (28 * t269)) * t270 + t156 * t476 + (t182 * t267 + 0.24e2 * t41 * t196 + t242 * t274 + t302 * t248 + t437 * t258 + 0.4e1 * t277 ^ 2 + (28 * t284 ^ 2)) * t273 + 0.8e1 * (-t385 * t59 - t42 * t478) * pkin(2) + (0.8e1 * t183 * t50 + 0.32e2 * t474 * t83) * pkin(1) + (t182 * t269 + 0.16e2 * t194 * t71 + (70 * t267) + t281 + t302) * t281;
t440 = t240 - 0.2e1 * t264;
t306 = 0.24e2 * t179 * t458 - t270 - ((21 * t269) + t521) * t281 - (t275 * t440 + t155 + t253 + (35 * t267)) * t273 - (t246 + (t254 + t439) * t269 + t275 * (-t264 + t430)) * t186;
t381 = -0.12e2 * t470;
t400 = -0.6e1 * t457;
t350 = t275 + t522;
t352 = t261 / 0.3e1 + t217 + t258;
t317 = -0.8e1 / 0.3e1 * t458 + t273 * t187 - 0.5e1 / 0.3e1 * t267 + t352 * t269 + t275 * t350;
t374 = -t486 / 0.2e1;
t475 = t192 * t267;
t446 = t207 - t264 / 0.6e1;
t351 = t275 + t446;
t125 = t499 + t224 + t351;
t315 = t226 + t351;
t126 = t222 + t315;
t72 = -t125 * t486 + t126 * t184;
t128 = t269 + t315;
t173 = 0.2e1 * t273 + t187;
t81 = t128 * t184 + t173 * t374;
t200 = -0.20e2 / 0.3e1 * t269;
t353 = 0.2e1 / 0.3e1 * t261 + t216 + t256;
t354 = 0.4e1 / 0.3e1 * t261 + 0.4e1 / 0.3e1 * t264 + t259;
t85 = -t281 + (t200 + t353) * t273 - (3 * t267) + t354 * t269 + t274;
t43 = t72 * t398 + t85 * t374 + t317 * t184 + (t183 * t81 - t228 * t475) * t426;
t389 = t284 * t184;
t324 = t192 * t389;
t399 = -0.4e1 * t457;
t346 = t255 + t440;
t221 = -0.2e1 / 0.3e1 * t264;
t445 = t211 + t221;
t73 = t281 + (t441 + t445) * t273 + t247 + t346 * t269 + t275 * (t275 + t445);
t88 = t251 + (t238 + t346) * t273 + (t221 + t347) * t186;
t55 = t73 * t184 - t486 * t88;
t314 = t269 + t350;
t130 = t252 + t314;
t131 = t174 + t522;
t74 = -t130 * t486 + t131 * t184;
t142 = t224 + t273 + t446;
t89 = t142 * t423 + t460 * t486;
t44 = t89 * t399 + (t414 * t74 - 0.8e1 * t324) * pkin(1) + t55;
t416 = -0.2e1 * t486;
t86 = -0.3e1 * t281 + (t200 + t354) * t273 + t353 * t269 + t432;
t90 = -0.5e1 / 0.3e1 * t281 + (-t269 + t352) * t273 + t275 * t314;
t56 = t184 * t86 + t416 * t90;
t349 = t210 - t264 / 0.2e1 + t275;
t127 = 0.3e1 / 0.2e1 * t273 + t249 + t349;
t134 = 0.4e1 * t301;
t146 = t184 + 0.2e1 * t486;
t57 = t188 * t184 + t134 * t193 + (t127 * t228 + t146 * t493) * t425;
t129 = t252 + 0.3e1 / 0.2e1 * t269 + t349;
t87 = t129 * t184 + t189 * t486 / 0.2e1;
t28 = t87 * t328 + t57 * t326 + t56 * t400 + t43 * t381 + (t525 * t228 - 0.6e1 * t55 * t493) * pkin(3) + (t229 * t306 + 0.6e1 * t44 * t478) * pkin(2);
t22 = t133 * t24 + t276 * t28;
t372 = t39 * t508;
t497 = t276 / 0.4e1;
t299 = t22 * t372 + t37 * t497;
t266 = 0.1e1 / pkin(4) ^ 2;
t459 = 0.1e1 / pkin(5) * t266;
t375 = t75 * t459;
t20 = t299 * t375;
t18 = t20 * t121;
t453 = t284 * t193;
t376 = -0.24e2 * t453;
t320 = t228 * t376;
t108 = t179 * t320 * t490;
t388 = pkin(2) * t465;
t336 = pkin(1) * t388;
t303 = t269 * t228 * t336;
t124 = -0.4e1 * t303;
t141 = t336 * t425;
t492 = pkin(1) * t267;
t307 = -0.64e2 * t460 * t492;
t342 = 0.32e2 / 0.3e1 * t267;
t311 = t192 * t342;
t312 = 0.64e2 / 0.3e1 * t143 * t284;
t322 = -0.96e2 * t154 * t179 * t284;
t491 = pkin(1) * t269;
t329 = -0.24e2 * t154 * t491;
t330 = -0.16e2 * t144 * t491;
t337 = t461 * t516;
t488 = pkin(3) * t178;
t338 = -0.4e1 * t154 * t488;
t339 = -0.48e2 * t82 * t491;
t469 = t196 * t280;
t362 = t164 * t469;
t167 = pkin(1) * t423;
t462 = t231 * t273;
t361 = t229 * t462;
t119 = t167 - 0.4e1 * t361;
t393 = pkin(1) * t524;
t392 = pkin(2) * t467;
t334 = pkin(3) * t392;
t449 = -0.2e1 * t334 + t167;
t45 = t119 * t399 + (pkin(1) * t149 * t518 + t169 * t512) * t229 + (-t449 * t147 + 0.2e1 * t229 ^ 2 * t393 + (t112 * t467 - t466 * t92) * pkin(2)) * t529;
t367 = t45 * t506;
t379 = 0.24e2 * t470;
t382 = -0.32e2 * t472;
t407 = 0.8e1 * t472;
t419 = 0.6e1 * t490;
t422 = -0.8e1 * t490;
t427 = -0.8e1 * t59 * pkin(3);
t428 = pkin(3) * t527;
t471 = t195 * t281;
t15 = t28 * t367 + (-0.32e2 * t124 * t133 + 0.8e1 * t141 * t276) * t363 + (0.24e2 * (0.4e1 * t110 * t471 * t486 - t362 * t57 + t43 * t462) * t276 + (0.96e2 * t362 * t53 - 0.48e2 * t41 * t462 - 0.64e2 * t471 * t71) * t133) * t229 + ((t24 + (t133 * t515 - 0.6e1 * t276 * t44) * t164) * t229 + (((t126 * t398 + t128 * t356 + t317) * t381 + t129 * t328 + t86 * t400 - 0.6e1 * t73 * t412 + t306) * t276 + (t108 * t514 + (t523 * t407 + (t131 * t356 - 0.8e1 * t142 * t457 - 0.8e1 * t396 + t73) * t419) * t276) * t164 + ((-pkin(1) * t311 - t193 * t312 + t230 * t330 + t338) * t379 + t192 * t307 + t193 * t322 + t230 * t339 + t427 + ((t337 - 0.4e1 * t453) * t382 + (t230 * t329 + t428) * t422) * t164) * t133 * t228) * t231) * pkin(2);
t373 = t37 * t508;
t498 = -t276 / 0.4e1;
t298 = t22 * t373 + t39 * t498;
t19 = t298 * t375;
t380 = 0.12e2 * t470;
t403 = -0.8e1 * t464;
t404 = 0.8e1 * t229 * t273;
t418 = t132 * t517;
t420 = 0.4e1 * t490;
t429 = t63 * t518;
t23 = (t141 * t406 + (t66 * t404 + t171 * t193 + ((t187 + 0.2e1 * t457) * t406 - t113 + (-0.4e1 * t174 * t183 + t409) * pkin(1)) * pkin(2)) * t231 + ((t141 + (t174 - 0.2e1 * t457) * t490) * t420 + (-0.24e2 * t469 * t486 + t429) * t229) * t164) * t276 + t40 * t367 + (0.24e2 * t114 * t229 * t362 + (t124 + (-0.8e1 / 0.3e1 * t453 - 0.2e1 * t488) * t392) * t380 - 0.24e2 * t51 * t361 + t108 + t303 * t526 + t334 * t528 + 0.6e1 * (-(pkin(1) * t403 + t418) * t196 * t524 + t52 * t184) * t164) * t133 + t34 * t184;
t415 = 0.2e1 * t183;
t487 = pkin(3) * t193;
t29 = t45 * t358 + t119 * t228 * t415 + ((-t228 * t276 + t68) * t231 + (t230 * t276 + (t163 * t519 + t111) * t228 + (-pkin(3) + 0.2e1 * t487 + t493) * t423) * t229) * pkin(2);
t297 = t76 * t299;
t30 = (t158 + t388) * t276 + t45 * t357 - 0.2e1 * t119 * t487 - (t167 - 0.4e1 * t334) * t147 - 0.2e1 * t196 * t393 + t148 * t392 + (t462 * t517 + (t150 * t519 - t230 * t93) * pkin(2)) * t229;
t369 = t37 * t47 / 0.8e1;
t370 = t39 * t507;
t333 = pkin(3) * t391;
t94 = 0.2e1 * t333 + t449;
t495 = t19 - (-t94 * t297 + (t30 * t497 + t45 * t369 + t15 * t372 + (t23 * t370 + t29 * t508) * t22) * t75) * t459;
t296 = t76 * t298;
t368 = -t39 * t47 / 0.8e1;
t371 = t37 * t507;
t5 = (-t94 * t296 + (t15 * t373 + t29 * t498 + t45 * t368 + (t23 * t371 + t30 * t508) * t22) * t75) * t459;
t10 = -t120 * t20 + t121 * t19;
t12 = -t120 * t19 - t18;
t289 = t12 ^ 2;
t9 = 0.1e1 / t289;
t500 = t10 * t9;
t7 = 0.1e1 / (t10 ^ 2 * t9 + 0.1e1);
t8 = 0.1e1 / t12;
t2 = 0.1e1 + ((t120 * t495 + t121 * t5 + t18) * t8 - (t495 * t121 + (-t20 - t5) * t120) * t500) * t7;
t501 = Ifges(3,3) * t2;
t308 = pkin(1) * t325;
t140 = -0.4e1 * t308;
t165 = pkin(1) * t416;
t309 = pkin(1) * t193 * t389;
t319 = t192 * t360;
t394 = pkin(1) * t453;
t327 = -0.48e2 * t394;
t395 = pkin(1) * t457;
t340 = 0.4e1 * t395;
t359 = t228 * t464;
t366 = t49 * t506;
t378 = -0.24e2 * t464;
t384 = -0.32e2 * t475;
t397 = pkin(1) * t458;
t408 = -0.8e1 * t472;
t410 = -0.4e1 * t474;
t421 = -0.6e1 * t490;
t451 = pkin(1) * t189 * t320 - 0.24e2 * t179 * t324;
t452 = -0.4e1 * t509;
t14 = ((t127 * t415 + t340 + t410) * t326 + (0.12e2 * t190 * t193 * t492 + t125 * t410 - 0.2e1 * t173 * t395 - 0.4e1 * t397) * t381 + 0.8e1 * t189 * t397 + 0.12e2 * t90 * t474 + 0.6e1 * t88 * t395 + (-t85 * t381 / 0.2e1 + t525) * t183) * t276 + t28 * t366 + t140 * t133 * t321 + ((0.6e1 * (-0.4e1 * t130 * t395 - t183 * t88 - 0.4e1 * t319) * t276 + t451 * t514) * t478 + ((-pkin(1) * t191 * t342 - t192 * t312 + t193 * t330 + t230 * t338) * t379 + t191 * t307 + t192 * t322 + t193 * t339 + t230 * t427 + ((t230 * t337 + t410) * t382 + (t193 * t329 + t230 * t428) * t422) * t164) * t133 * t229) * pkin(2) + (0.12e2 * (-(t184 * t311 + t403 * t72) * t470 - 0.8e1 * t179 * t475 * t184 - 0.4e1 * t87 * t394 + t56 * t464) * t276 + (0.16e2 * (t177 * t513 - t162 + t327 + t384) * t473 + (t109 * t327 + t180 * t384 - 0.28e2 * t464 * t60) * t379 - 0.4e1 * t139 * t192 - 0.96e2 * t83 * t394 - 0.48e2 * t54 * t464) * t133 + ((0.32e2 * t53 * t472 + t490 * t515) * t133 - t24 + (t44 * t421 + t57 * t408 + (-0.24e2 * t165 + 0.64e2 * t359) * t473) * t276 + ((0.48e2 * t470 * t81 + 0.6e1 * t55) * t276 + (-0.144e3 * t470 * t61 - 0.8e1 * t50) * t133) * pkin(1)) * pkin(3) + ((0.2e1 * (-t134 * t230 - t146 * t509) * t407 + (t402 * t89 + t452 * t74 + 0.24e2 * t309) * t419) * t276 + ((pkin(1) * t376 + t101 * t513 + t107 * t452) * t382 + (t378 * t79 - 0.6e1 * t509 * t65) * t422) * t133) * t164) * t228;
t21 = (t269 * t190 * t408 + (t183 * t185 - 0.2e1 * t474) * t406 + 0.4e1 * t319 + t172 * t340 + t122 * t183) * t276 + t40 * t366 + ((-0.8e1 / 0.3e1 * t324 + t140 - 0.2e1 * t178 * t333) * t380 + t308 * t526 + t333 * t528 + t451) * t133 + ((0.8e1 * t135 * t196 * t463 + t97 * t513 - 0.24e2 * t309) * t276 + (0.12e2 * (-t464 * t96 - t394) * t380 + t77 * t378) * t133 + (t231 * t276 * t429 + (t114 * t407 + t52 * t419) * t133 - t34 + ((pkin(2) * t196 * t404 + 0.4e1 * t95) * t276 + (-0.48e2 * t104 * t470 - 0.6e1 * t62) * t133) * pkin(1)) * pkin(3)) * t228 + ((t407 * t183 + ((-pkin(3) * t172 + t158 * t510) * t230 + (-t145 * t486 - t457) * t519) * t420) * t276 + ((t165 - 0.8e1 * t359) * t408 + ((-0.2e1 * t170 * t228 + t184 * t418) * t230 + (-0.4e1 * t102 * t486 - 0.8e1 * t325) * pkin(1)) * t421) * t133) * t164;
t496 = -t19 - (t98 * t297 + (t32 * t497 + t49 * t369 + t14 * t372 + (t21 * t370 + t31 * t508) * t22) * t75) * t459;
t494 = t20 - (t98 * t296 + (t14 * t373 + t31 * t498 + t49 * t368 + (t21 * t371 + t32 * t508) * t22) * t75) * t459;
t485 = m(3) * t273;
t484 = mrSges(3,1) * t12;
t483 = t37 * mrSges(5,1);
t482 = t39 * mrSges(5,2);
t411 = pkin(3) * t265 * t75;
t365 = t94 * t504;
t16 = 0.2e1 * (t29 * t505 + t365 * t39) * t310 - 0.2e1 * (t30 * t505 + t365 * t37) * t304;
t1 = ((t120 * t496 - t121 * t494) * t8 - (t120 * t494 + t121 * t496) * t500) * t7;
t3 = [t289 * t485 + t16 ^ 2 * Ifges(5,3) + Ifges(2,3) + (-0.2e1 * pkin(2) * t484 + t501) * t2 + (0.2e1 * pkin(2) * mrSges(3,2) * t2 + t10 * t485) * t10; (t530 + (-t482 / 0.2e1 + t483 / 0.2e1) * t411) * t16 + (t501 + (mrSges(3,2) * t10 - t484) * pkin(2)) * t1; t1 ^ 2 * Ifges(3,3) + Ifges(4,3) + m(5) * (t38 / 0.2e1 + t294 / 0.2e1) * t76 * t266 * t224 + ((-t482 + t483) * t411 + t530) * t17;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t3(1), t3(2); t3(2), t3(3);];
Mq = res;
