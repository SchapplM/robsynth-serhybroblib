% Calculate joint inertia matrix for
% palh1m2OL
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
% Mq [13x13]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_inertiaJ_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_inertiaJ_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2OL_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2OL_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:37
% EndTime: 2020-05-02 21:17:42
% DurationCPUTime: 5.16s
% Computational Cost: add. (10301->592), mult. (9770->894), div. (0->0), fcn. (9739->27), ass. (0->336)
t306 = sin(qJ(1));
t311 = cos(qJ(1));
t301 = qJ(2) + qJ(3);
t286 = sin(t301);
t289 = cos(t301);
t353 = Icges(4,5) * t289 - Icges(4,6) * t286;
t167 = -Icges(4,3) * t311 + t353 * t306;
t168 = Icges(4,3) * t306 + t353 * t311;
t298 = t311 ^ 2;
t414 = Icges(4,4) * t289;
t360 = -Icges(4,2) * t286 + t414;
t174 = Icges(4,6) * t306 + t360 * t311;
t415 = Icges(4,4) * t286;
t367 = Icges(4,1) * t289 - t415;
t180 = Icges(4,5) * t306 + t367 * t311;
t336 = -t174 * t286 + t180 * t289;
t173 = -Icges(4,6) * t311 + t360 * t306;
t179 = -Icges(4,5) * t311 + t367 * t306;
t337 = t173 * t286 - t179 * t289;
t291 = qJ(4) + t301;
t279 = cos(t291);
t308 = cos(qJ(5));
t394 = t308 * t311;
t303 = sin(qJ(5));
t396 = t306 * t303;
t201 = -t279 * t396 - t394;
t395 = t306 * t308;
t398 = t303 * t311;
t202 = t279 * t395 - t398;
t277 = sin(t291);
t402 = t277 * t306;
t91 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t402;
t93 = Icges(6,4) * t202 + Icges(6,2) * t201 + Icges(6,6) * t402;
t95 = Icges(6,1) * t202 + Icges(6,4) * t201 + Icges(6,5) * t402;
t19 = t201 * t93 + t202 * t95 + t91 * t402;
t203 = -t279 * t398 + t395;
t204 = t279 * t394 + t396;
t401 = t277 * t311;
t92 = Icges(6,5) * t204 + Icges(6,6) * t203 + Icges(6,3) * t401;
t94 = Icges(6,4) * t204 + Icges(6,2) * t203 + Icges(6,6) * t401;
t96 = Icges(6,1) * t204 + Icges(6,4) * t203 + Icges(6,5) * t401;
t20 = t201 * t94 + t202 * t96 + t92 * t402;
t10 = -t19 * t311 + t20 * t306;
t352 = Icges(5,5) * t279 - Icges(5,6) * t277;
t144 = -Icges(5,3) * t311 + t352 * t306;
t145 = Icges(5,3) * t306 + t352 * t311;
t412 = Icges(5,4) * t279;
t359 = -Icges(5,2) * t277 + t412;
t149 = Icges(5,6) * t306 + t359 * t311;
t413 = Icges(5,4) * t277;
t366 = Icges(5,1) * t279 - t413;
t153 = Icges(5,5) * t306 + t366 * t311;
t342 = -t149 * t277 + t153 * t279;
t148 = -Icges(5,6) * t311 + t359 * t306;
t152 = -Icges(5,5) * t311 + t366 * t306;
t343 = t148 * t277 - t152 * t279;
t426 = -t298 * t144 - (t342 * t306 + (-t145 + t343) * t311) * t306 - t10;
t440 = -t298 * t167 - (t336 * t306 + (-t168 + t337) * t311) * t306 + t426;
t257 = rSges(6,1) * t303 + rSges(6,2) * t308;
t258 = rSges(6,1) * t308 - rSges(6,2) * t303;
t400 = t279 * t311;
t104 = rSges(6,3) * t401 + t306 * t257 + t258 * t400;
t439 = pkin(9) * t400 + pkin(11) * t401 + t104;
t438 = t306 * t311;
t299 = qJ(2) + qJ(8);
t290 = qJ(9) + t299;
t276 = sin(t290);
t278 = cos(t290);
t219 = rSges(10,1) * t276 + rSges(10,2) * t278;
t284 = sin(t299);
t437 = -pkin(2) * t284 + t219;
t300 = qJ(2) + qJ(7);
t285 = sin(t300);
t288 = cos(t300);
t369 = rSges(8,1) * t285 + rSges(8,2) * t288;
t436 = -rSges(8,3) * t311 - t369 * t306;
t424 = rSges(4,2) * t286;
t425 = rSges(4,1) * t289;
t435 = (-t424 + t425) * t306 - rSges(4,3) * t311;
t297 = t306 ^ 2;
t434 = t306 / 0.2e1;
t433 = -t311 / 0.2e1;
t21 = t203 * t93 + t204 * t95 + t91 * t401;
t22 = t203 * t94 + t204 * t96 + t92 * t401;
t11 = -t21 * t311 + t22 * t306;
t432 = (t297 * t145 + t11 + (t343 * t311 + (-t144 + t342) * t306) * t311) * t306;
t305 = sin(qJ(2));
t431 = pkin(1) * t305;
t310 = cos(qJ(2));
t430 = pkin(1) * t310;
t296 = -qJ(7) + pkin(19);
t281 = sin(t296);
t282 = cos(t296);
t428 = pkin(4) * (t281 * t305 + t282 * t310);
t304 = sin(qJ(3));
t309 = cos(qJ(3));
t397 = t305 * t309;
t427 = pkin(5) * (t304 * t310 + t397);
t422 = rSges(5,3) * t311;
t31 = -t279 * t91 + (-t303 * t93 + t308 * t95) * t277;
t420 = t31 * t311;
t348 = Icges(10,5) * t276 + Icges(10,6) * t278;
t142 = -Icges(10,3) * t311 + t348 * t306;
t143 = Icges(10,3) * t306 + t348 * t311;
t405 = Icges(10,4) * t276;
t355 = Icges(10,2) * t278 + t405;
t147 = Icges(10,6) * t306 + t355 * t311;
t404 = Icges(10,4) * t278;
t362 = Icges(10,1) * t276 + t404;
t151 = Icges(10,5) * t306 + t362 * t311;
t344 = t147 * t278 + t151 * t276;
t146 = -Icges(10,6) * t311 + t355 * t306;
t150 = -Icges(10,5) * t311 + t362 * t306;
t345 = -t146 * t278 - t150 * t276;
t17 = t298 * t142 + (t344 * t306 + (-t143 + t345) * t311) * t306;
t419 = t311 * t17;
t32 = -t279 * t92 + (-t303 * t94 + t308 * t96) * t277;
t418 = t32 * t306;
t417 = Icges(3,4) * t305;
t416 = Icges(3,4) * t310;
t302 = sin(qJ(6));
t411 = Icges(7,4) * t302;
t307 = cos(qJ(6));
t410 = Icges(7,4) * t307;
t409 = Icges(8,4) * t285;
t408 = Icges(8,4) * t288;
t407 = Icges(9,4) * t284;
t287 = cos(t299);
t406 = Icges(9,4) * t287;
t403 = t257 * t311;
t128 = -Icges(6,6) * t279 + (Icges(6,4) * t308 - Icges(6,2) * t303) * t277;
t399 = t303 * t128;
t321 = rSges(5,1) * t400 - rSges(5,2) * t401 + t306 * rSges(5,3);
t370 = rSges(5,1) * t279 - rSges(5,2) * t277;
t85 = t306 * (t370 * t306 - t422) + t311 * t321;
t132 = -rSges(6,3) * t279 + t258 * t277;
t393 = -pkin(9) * t277 + pkin(11) * t279 - t132;
t292 = t306 * rSges(8,3);
t88 = t306 * t436 + t311 * (-t369 * t311 + t292);
t390 = t306 * rSges(4,3) + t311 * t425;
t89 = t306 * t435 + t311 * (-t311 * t424 + t390);
t280 = pkin(5) * t304 + pkin(1);
t208 = pkin(5) * t397 + t280 * t310;
t216 = rSges(5,1) * t277 + rSges(5,2) * t279;
t392 = -t208 - t216;
t389 = t297 + t298;
t283 = -qJ(10) + t296;
t273 = -qJ(2) + t283;
t266 = sin(t273);
t388 = Icges(11,4) * t266;
t267 = cos(t273);
t387 = Icges(11,4) * t267;
t386 = t306 * (t297 * t168 + (t337 * t311 + (-t167 + t336) * t306) * t311) + t432;
t385 = -t208 + t393;
t274 = pkin(5) * t289;
t76 = t274 + t85;
t233 = rSges(4,1) * t286 + rSges(4,2) * t289;
t384 = -t233 - t430;
t235 = rSges(8,1) * t288 - rSges(8,2) * t285;
t383 = -t235 - t430;
t217 = rSges(10,1) * t278 - t276 * rSges(10,2);
t382 = -pkin(2) * t287 + t217;
t381 = -t216 - t427;
t322 = -Icges(11,5) * t266 + Icges(11,6) * t267;
t121 = -Icges(11,3) * t311 + t322 * t306;
t122 = Icges(11,3) * t306 + t322 * t311;
t323 = Icges(11,2) * t267 - t388;
t124 = Icges(11,6) * t306 + t323 * t311;
t324 = -Icges(11,1) * t266 + t387;
t126 = Icges(11,5) * t306 + t324 * t311;
t346 = t124 * t267 - t126 * t266;
t123 = -Icges(11,6) * t311 + t323 * t306;
t125 = -Icges(11,5) * t311 + t324 * t306;
t347 = -t123 * t267 + t125 * t266;
t13 = t306 * (t297 * t122 + (t347 * t311 + (-t121 + t346) * t306) * t311);
t14 = t298 * t121 + (t346 * t306 + (-t122 + t347) * t311) * t306;
t380 = -t311 * t14 + t13;
t15 = t306 * (t297 * t143 + (t345 * t311 + (-t142 + t344) * t306) * t311);
t379 = t15 - t419;
t133 = t382 * t306;
t134 = t382 * t311;
t349 = -Icges(9,5) * t284 - Icges(9,6) * t287;
t163 = -Icges(9,3) * t311 + t349 * t306;
t164 = Icges(9,3) * t306 + t349 * t311;
t234 = rSges(9,1) * t287 - rSges(9,2) * t284;
t236 = -rSges(9,1) * t284 - rSges(9,2) * t287;
t356 = -Icges(9,2) * t287 - t407;
t170 = Icges(9,6) * t306 + t356 * t311;
t363 = -Icges(9,1) * t284 - t406;
t176 = Icges(9,5) * t306 + t363 * t311;
t340 = -t170 * t287 - t176 * t284;
t169 = -Icges(9,6) * t311 + t356 * t306;
t175 = -Icges(9,5) * t311 + t363 * t306;
t341 = t169 * t287 + t175 * t284;
t378 = t15 + t306 * (t297 * t164 + (t341 * t311 + (-t163 + t340) * t306) * t311) + m(10) * (t133 ^ 2 + t134 ^ 2 + t437 ^ 2) + m(9) * (t389 * t234 ^ 2 + t236 ^ 2);
t184 = -Icges(11,5) * t267 - Icges(11,6) * t266;
t185 = -Icges(11,2) * t266 - t387;
t186 = -Icges(11,1) * t267 - t388;
t335 = t185 * t267 - t186 * t266;
t377 = (-t124 * t266 - t126 * t267 + t306 * t184 + t335 * t311) * t434 + (-t123 * t266 - t125 * t267 - t184 * t311 + t335 * t306) * t433;
t210 = -Icges(10,5) * t278 + Icges(10,6) * t276;
t212 = Icges(10,2) * t276 - t404;
t214 = -Icges(10,1) * t278 + t405;
t330 = t212 * t278 + t214 * t276;
t376 = (t147 * t276 - t151 * t278 + t306 * t210 + t330 * t311) * t434 + (t146 * t276 - t150 * t278 - t210 * t311 + t330 * t306) * t433;
t199 = -rSges(11,1) * t266 + rSges(11,2) * t267;
t127 = -Icges(6,3) * t279 + (Icges(6,5) * t308 - Icges(6,6) * t303) * t277;
t129 = -Icges(6,5) * t279 + (Icges(6,1) * t308 - Icges(6,4) * t303) * t277;
t39 = t127 * t402 + t128 * t201 + t129 * t202;
t4 = -t39 * t279 + (t19 * t306 + t20 * t311) * t277;
t40 = t127 * t401 + t203 * t128 + t204 * t129;
t5 = -t40 * t279 + (t21 * t306 + t22 * t311) * t277;
t375 = t10 * t402 / 0.2e1 + t4 * t433 + t5 * t434 - t279 * (t418 - t420) / 0.2e1 + t11 * t401 / 0.2e1;
t103 = -t403 + (rSges(6,3) * t277 + t258 * t279) * t306;
t38 = t306 * t103 + t297 * (pkin(9) * t279 + pkin(11) * t277) + t439 * t311;
t374 = t297 / 0.2e1 + t298 / 0.2e1;
t373 = t393 - t427;
t372 = -t428 - t430;
t137 = t199 - pkin(4) * sin(qJ(2) - t296);
t254 = rSges(3,1) * t305 + rSges(3,2) * t310;
t256 = rSges(7,1) * t307 - rSges(7,2) * t302;
t34 = t274 + t38;
t368 = -Icges(3,1) * t305 - t416;
t365 = Icges(7,1) * t307 - t411;
t364 = -Icges(8,1) * t285 - t408;
t361 = -Icges(3,2) * t310 - t417;
t358 = -Icges(7,2) * t302 + t410;
t357 = -Icges(8,2) * t288 - t409;
t354 = -Icges(3,5) * t305 - Icges(3,6) * t310;
t351 = Icges(7,5) * t307 - Icges(7,6) * t302;
t350 = -Icges(8,5) * t285 - Icges(8,6) * t288;
t171 = -Icges(8,6) * t311 + t357 * t306;
t177 = -Icges(8,5) * t311 + t364 * t306;
t339 = t171 * t288 + t177 * t285;
t172 = Icges(8,6) * t306 + t357 * t311;
t178 = Icges(8,5) * t306 + t364 * t311;
t338 = -t172 * t288 - t178 * t285;
t213 = Icges(5,2) * t279 + t413;
t215 = Icges(5,1) * t277 + t412;
t329 = -t213 * t277 + t215 * t279;
t227 = -Icges(9,2) * t284 + t406;
t230 = Icges(9,1) * t287 - t407;
t328 = -t227 * t287 - t230 * t284;
t228 = -Icges(8,2) * t285 + t408;
t231 = Icges(8,1) * t288 - t409;
t327 = -t228 * t288 - t231 * t285;
t229 = Icges(4,2) * t289 + t415;
t232 = Icges(4,1) * t286 + t414;
t326 = -t229 * t286 + t232 * t289;
t325 = t426 * t311 + t432;
t165 = -Icges(8,3) * t311 + t350 * t306;
t166 = Icges(8,3) * t306 + t350 * t311;
t24 = t306 * (t297 * t166 + (t339 * t311 + (-t165 + t338) * t306) * t311);
t27 = t298 * t165 + (t338 * t306 + (-t166 + t339) * t311) * t306;
t320 = t13 + t24 + (-t27 - t14) * t311;
t319 = -pkin(14) + t256;
t211 = Icges(5,5) * t277 + Icges(5,6) * t279;
t318 = t418 / 0.2e1 - t420 / 0.2e1 + (t149 * t279 + t153 * t277 + t306 * t211 + t329 * t311 + t40) * t434 + (t148 * t279 + t152 * t277 - t211 * t311 + t329 * t306 + t39) * t433;
t225 = Icges(8,5) * t288 - Icges(8,6) * t285;
t317 = t377 + (-t172 * t285 + t178 * t288 + t306 * t225 + t327 * t311) * t434 + (-t171 * t285 + t177 * t288 - t225 * t311 + t327 * t306) * t433;
t316 = t440 * t311 + t386;
t245 = rSges(11,1) * t305 + rSges(11,2) * t310;
t246 = rSges(11,1) * t310 - rSges(11,2) * t305;
t271 = sin(t283);
t272 = cos(t283);
t314 = -pkin(4) * (t281 * t310 - t282 * t305) - t245 * t272 + t246 * t271;
t107 = rSges(10,3) * t311 + (-pkin(15) - t437) * t306;
t295 = t311 * pkin(15);
t108 = t306 * rSges(10,3) + t311 * t437 + t295;
t135 = rSges(9,3) * t311 + (-pkin(15) - t236) * t306;
t136 = t306 * rSges(9,3) + t236 * t311 + t295;
t224 = Icges(9,5) * t287 - Icges(9,6) * t284;
t313 = m(10) * (t107 * t134 + t108 * t133) + m(9) * (-t135 * t311 - t136 * t306) * t234 + t376 + (-t170 * t284 + t176 * t287 + t306 * t224 + t328 * t311) * t434 + (-t169 * t284 + t175 * t287 - t224 * t311 + t328 * t306) * t433;
t226 = Icges(4,5) * t286 + Icges(4,6) * t289;
t312 = t318 + (t174 * t289 + t180 * t286 + t306 * t226 + t326 * t311) * t434 + (t173 * t289 + t179 * t286 - t226 * t311 + t326 * t306) * t433;
t275 = -pkin(15) + t431;
t263 = t275 * t306;
t261 = rSges(2,1) * t311 - t306 * rSges(2,2);
t260 = rSges(3,1) * t310 - rSges(3,2) * t305;
t255 = -t306 * rSges(2,1) - rSges(2,2) * t311;
t253 = rSges(7,1) * t302 + rSges(7,2) * t307;
t205 = pkin(5) * t309 * t310 - t280 * t305 + pkin(15);
t190 = Icges(3,3) * t306 + t354 * t311;
t189 = -Icges(3,3) * t311 + t354 * t306;
t188 = Icges(7,3) * t306 + t351 * t311;
t187 = -Icges(7,3) * t311 + t351 * t306;
t183 = t205 * t311;
t181 = rSges(3,3) * t311 + (-pkin(15) + t254) * t306;
t162 = t306 * rSges(7,3) + t319 * t311;
t161 = t306 * rSges(3,3) - t254 * t311 + t295;
t160 = rSges(7,3) * t311 - t319 * t306;
t157 = t383 * t311;
t156 = t384 * t311;
t155 = t383 * t306;
t154 = t384 * t306;
t120 = t137 - t431;
t119 = (-t275 - t424) * t311 + t390;
t118 = t292 + (-t275 - t369) * t311;
t117 = t263 - t435;
t116 = t263 - t436;
t115 = t245 * t271 + t246 * t272;
t113 = t381 * t311;
t112 = t381 * t306;
t111 = t115 * t311;
t110 = t115 * t306;
t109 = t277 * t308 * t129;
t106 = t392 * t311;
t105 = t392 * t306;
t100 = t183 + t321;
t99 = t422 + (-t205 - t370) * t306;
t98 = t393 * t311;
t97 = t393 * t306;
t87 = t89 - t431;
t86 = t88 - t431;
t84 = -t311 * t428 + t111;
t83 = -t306 * t428 + t110;
t69 = t372 * t311 + t111;
t68 = t372 * t306 + t110;
t67 = t373 * t311;
t66 = t373 * t306;
t65 = t76 - t431;
t60 = t385 * t311;
t59 = t385 * t306;
t58 = t306 * rSges(11,3) + (-t275 - t314) * t311;
t57 = rSges(11,3) * t311 + t314 * t306 + t263;
t50 = t183 + t439;
t49 = t403 + (-t205 + (-pkin(9) - t258) * t279 + (-pkin(11) - rSges(6,3)) * t277) * t306;
t47 = -t279 * t104 - t132 * t401;
t46 = t103 * t279 + t132 * t402;
t42 = (t103 * t311 - t104 * t306) * t277;
t41 = -t279 * t127 - t277 * t399 + t109;
t33 = t34 - t431;
t26 = t298 * t163 + (t340 * t306 + (-t164 + t341) * t311) * t306;
t6 = m(10) * (t219 * t437 + (t133 * t306 + t134 * t311) * t217) + t379;
t1 = (-t17 - t26) * t311 + t378;
t2 = [t302 * (Icges(7,1) * t302 + t410) + t307 * (Icges(7,2) * t307 + t411) - t305 * (-Icges(3,2) * t305 + t416) + t310 * (Icges(3,1) * t310 - t417) + (t215 - t399) * t277 + t109 + (-t127 + t213) * t279 + m(6) * (t49 ^ 2 + t50 ^ 2) + m(11) * (t57 ^ 2 + t58 ^ 2) + m(5) * (t100 ^ 2 + t99 ^ 2) + m(10) * (t107 ^ 2 + t108 ^ 2) + m(8) * (t116 ^ 2 + t118 ^ 2) + m(4) * (t117 ^ 2 + t119 ^ 2) + m(9) * (t135 ^ 2 + t136 ^ 2) + m(7) * (t160 ^ 2 + t162 ^ 2) + m(3) * (t161 ^ 2 + t181 ^ 2) + m(2) * (t255 ^ 2 + t261 ^ 2) + t289 * t229 + t286 * t232 + t287 * t230 + t288 * t231 - t284 * t227 - t285 * t228 - t278 * t214 + t276 * t212 - t266 * t185 - t267 * t186 + Icges(2,3); t312 + m(3) * (-t161 * t306 - t181 * t311) * t260 + t374 * (Icges(3,5) * t310 - Icges(3,6) * t305) + t317 + m(6) * (t49 * t60 + t50 * t59) + m(11) * (t57 * t69 + t58 * t68) + m(5) * (t100 * t105 + t106 * t99) + m(4) * (t117 * t156 + t119 * t154) + m(8) * (t116 * t157 + t118 * t155) + t313 + (-(-Icges(3,6) * t311 + t361 * t306) * t305 + (-Icges(3,5) * t311 + t368 * t306) * t310) * t433 + (-(Icges(3,6) * t306 + t361 * t311) * t305 + (Icges(3,5) * t306 + t368 * t311) * t310) * t434; m(6) * (t33 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t105 ^ 2 + t106 ^ 2 + t65 ^ 2) + m(4) * (t154 ^ 2 + t156 ^ 2 + t87 ^ 2) + m(11) * (t120 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(8) * (t155 ^ 2 + t157 ^ 2 + t86 ^ 2) + t378 + t306 * t297 * t190 + t24 + m(3) * (t389 * t260 ^ 2 + t254 ^ 2) - t419 + t380 + t386 + (-t298 * t189 - t26 - t27 + (-t306 * t189 + t311 * t190) * t306 + t440) * t311; t312 + m(6) * (t49 * t67 + t50 * t66) + m(5) * (t100 * t112 + t113 * t99) + m(4) * (-t117 * t311 - t119 * t306) * t233; m(6) * (t33 * t34 + t59 * t66 + t60 * t67) + m(5) * (t105 * t112 + t106 * t113 + t65 * t76) + m(4) * (t89 * t87 + (-t154 * t306 - t156 * t311) * t233) + t316; m(6) * (t34 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2 + t76 ^ 2) + m(4) * (t389 * t233 ^ 2 + t89 ^ 2) + t316; m(6) * (t49 * t98 + t50 * t97) + m(5) * (-t100 * t306 - t311 * t99) * t216 + t318; m(6) * (t33 * t38 + t59 * t97 + t60 * t98) + m(5) * (t85 * t65 + (-t105 * t306 - t106 * t311) * t216) + t325; m(6) * (t34 * t38 + t66 * t97 + t67 * t98) + m(5) * (t85 * t76 + (-t112 * t306 - t113 * t311) * t216) + t325; m(5) * (t389 * t216 ^ 2 + t85 ^ 2) + m(6) * (t38 ^ 2 + t97 ^ 2 + t98 ^ 2) + t325; -t41 * t279 + m(6) * (t46 * t49 + t47 * t50) + ((t32 / 0.2e1 + t40 / 0.2e1) * t311 + (t31 / 0.2e1 + t39 / 0.2e1) * t306) * t277; m(6) * (t33 * t42 + t46 * t60 + t47 * t59) + t375; m(6) * (t34 * t42 + t46 * t67 + t47 * t66) + t375; m(6) * (t38 * t42 + t46 * t98 + t47 * t97) + t375; t279 ^ 2 * t41 + m(6) * (t42 ^ 2 + t46 ^ 2 + t47 ^ 2) + (t311 * t5 + t306 * t4 - t279 * (t306 * t31 + t311 * t32)) * t277; ((Icges(7,6) * t306 + t358 * t311) * t307 + (Icges(7,5) * t306 + t365 * t311) * t302) * t434 + ((-Icges(7,6) * t311 + t358 * t306) * t307 + (-Icges(7,5) * t311 + t365 * t306) * t302) * t433 + m(7) * (-t160 * t311 - t162 * t306) * t253 + t374 * (Icges(7,5) * t302 + Icges(7,6) * t307); 0; 0; 0; 0; t306 * (-t187 * t438 + t297 * t188) - t311 * (t298 * t187 - t188 * t438) + m(7) * (t253 ^ 2 * t389 + t256 ^ 2); m(11) * (t57 * t84 + t58 * t83) + m(8) * (-t116 * t311 - t118 * t306) * t235 + t317; m(11) * (t120 * t137 + t68 * t83 + t69 * t84) + m(8) * (t88 * t86 + (-t155 * t306 - t157 * t311) * t235) + t320; 0; 0; 0; 0; m(8) * (t235 ^ 2 * t389 + t88 ^ 2) + m(11) * (t137 ^ 2 + t83 ^ 2 + t84 ^ 2) + t320; t313; t1; 0; 0; 0; 0; 0; t1; m(10) * (t107 * t311 + t108 * t306) * t217 + t376; t6; 0; 0; 0; 0; 0; t6; m(10) * (t217 ^ 2 * t389 + t219 ^ 2) + t379; m(11) * (t306 * t58 + t311 * t57) * t115 + t377; m(11) * (t199 * t120 + (t306 * t68 + t311 * t69) * t115) + t380; 0; 0; 0; 0; m(11) * (t199 * t137 + (t306 * t83 + t311 * t84) * t115) + t380; 0; 0; m(11) * (t115 ^ 2 * t389 + t199 ^ 2) + t380; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_13_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11), t2(16), t2(22), t2(29), t2(37), t2(46), t2(56), t2(67), t2(79); t2(2), t2(3), t2(5), t2(8), t2(12), t2(17), t2(23), t2(30), t2(38), t2(47), t2(57), t2(68), t2(80); t2(4), t2(5), t2(6), t2(9), t2(13), t2(18), t2(24), t2(31), t2(39), t2(48), t2(58), t2(69), t2(81); t2(7), t2(8), t2(9), t2(10), t2(14), t2(19), t2(25), t2(32), t2(40), t2(49), t2(59), t2(70), t2(82); t2(11), t2(12), t2(13), t2(14), t2(15), t2(20), t2(26), t2(33), t2(41), t2(50), t2(60), t2(71), t2(83); t2(16), t2(17), t2(18), t2(19), t2(20), t2(21), t2(27), t2(34), t2(42), t2(51), t2(61), t2(72), t2(84); t2(22), t2(23), t2(24), t2(25), t2(26), t2(27), t2(28), t2(35), t2(43), t2(52), t2(62), t2(73), t2(85); t2(29), t2(30), t2(31), t2(32), t2(33), t2(34), t2(35), t2(36), t2(44), t2(53), t2(63), t2(74), t2(86); t2(37), t2(38), t2(39), t2(40), t2(41), t2(42), t2(43), t2(44), t2(45), t2(54), t2(64), t2(75), t2(87); t2(46), t2(47), t2(48), t2(49), t2(50), t2(51), t2(52), t2(53), t2(54), t2(55), t2(65), t2(76), t2(88); t2(56), t2(57), t2(58), t2(59), t2(60), t2(61), t2(62), t2(63), t2(64), t2(65), t2(66), t2(77), t2(89); t2(67), t2(68), t2(69), t2(70), t2(71), t2(72), t2(73), t2(74), t2(75), t2(76), t2(77), t2(78), t2(90); t2(79), t2(80), t2(81), t2(82), t2(83), t2(84), t2(85), t2(86), t2(87), t2(88), t2(89), t2(90), t2(91);];
Mq = res;
