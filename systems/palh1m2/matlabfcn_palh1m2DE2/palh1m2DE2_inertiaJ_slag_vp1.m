% Calculate joint inertia matrix for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_inertiaJ_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE2_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:21
% EndTime: 2020-05-02 20:58:28
% DurationCPUTime: 6.77s
% Computational Cost: add. (17033->503), mult. (24824->721), div. (0->0), fcn. (37541->74), ass. (0->299)
t296 = sin(pkin(19));
t299 = cos(pkin(19));
t301 = sin(qJ(3));
t307 = cos(qJ(3));
t191 = t296 * t307 + t299 * t301;
t193 = t296 * t301 - t299 * t307;
t146 = qJ(2) + atan2(t193, t191);
t142 = sin(t146);
t143 = cos(t146);
t289 = qJ(2) + qJ(3);
t268 = sin(t289);
t272 = cos(t289);
t304 = sin(pkin(18));
t310 = cos(pkin(18));
t200 = -t301 * t310 + t304 * t307;
t201 = t301 * t304 + t307 * t310;
t302 = sin(qJ(2));
t308 = cos(qJ(2));
t135 = -t200 * t308 + t201 * t302;
t136 = t200 * t302 + t201 * t308;
t294 = sin(pkin(22));
t297 = cos(pkin(22));
t194 = t294 * t304 + t297 * t310;
t195 = t294 * t310 - t297 * t304;
t284 = pkin(22) + pkin(21);
t263 = sin(t284);
t264 = cos(t284);
t46 = -qJ(2) - atan2(t194 * t308 - t195 * t302, t194 * t302 + t195 * t308) + pkin(21) - atan2(t135 * t264 + t136 * t263, -t135 * t263 + t136 * t264);
t44 = sin(t46);
t45 = cos(t46);
t453 = -Icges(11,5) * t44 + Icges(4,5) * t272 - Icges(9,5) * t142 + Icges(11,6) * t45 - Icges(4,6) * t268 - Icges(9,6) * t143;
t452 = Icges(11,3) + Icges(4,3) + Icges(9,3);
t303 = sin(qJ(1));
t309 = cos(qJ(1));
t451 = t452 * t303 + t453 * t309;
t450 = -t453 * t303 + t452 * t309;
t449 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t392 = Icges(4,4) * t272;
t344 = -Icges(4,2) * t268 + t392;
t156 = -Icges(4,6) * t309 + t344 * t303;
t393 = Icges(4,4) * t268;
t350 = Icges(4,1) * t272 - t393;
t158 = -Icges(4,5) * t309 + t350 * t303;
t397 = Icges(11,4) * t44;
t346 = Icges(11,2) * t45 - t397;
t26 = -Icges(11,6) * t309 + t346 * t303;
t396 = Icges(11,4) * t45;
t352 = -Icges(11,1) * t44 + t396;
t28 = -Icges(11,5) * t309 + t352 * t303;
t390 = Icges(9,4) * t142;
t342 = -Icges(9,2) * t143 - t390;
t89 = -Icges(9,6) * t309 + t342 * t303;
t389 = Icges(9,4) * t143;
t348 = -Icges(9,1) * t142 - t389;
t91 = -Icges(9,5) * t309 + t348 * t303;
t448 = t142 * t91 + t143 * t89 + t156 * t268 - t158 * t272 - t26 * t45 + t28 * t44;
t114 = atan2(t193, -t191) + t146;
t112 = sin(t114);
t113 = cos(t114);
t305 = sin(pkin(17));
t311 = cos(pkin(17));
t199 = t304 * t311 - t305 * t310;
t203 = t304 * t305 + t310 * t311;
t137 = t199 * t302 + t203 * t308;
t368 = t199 * t308 - t203 * t302;
t447 = -Icges(3,5) * t302 + Icges(7,5) * t368 + Icges(10,5) * t112 - Icges(3,6) * t308 - Icges(7,6) * t137 + Icges(10,6) * t113;
t157 = Icges(4,6) * t303 + t344 * t309;
t159 = Icges(4,5) * t303 + t350 * t309;
t27 = Icges(11,6) * t303 + t346 * t309;
t29 = Icges(11,5) * t303 + t352 * t309;
t90 = Icges(9,6) * t303 + t342 * t309;
t92 = Icges(9,5) * t303 + t348 * t309;
t446 = t142 * t92 + t143 * t90 + t157 * t268 - t159 * t272 - t27 * t45 + t29 * t44;
t290 = qJ(1) + qJ(4);
t273 = cos(t290);
t369 = pkin(18) - t284;
t256 = -pkin(20) + t369;
t324 = qJ(4) + t256;
t234 = qJ(1) + t324;
t325 = -qJ(4) + t256;
t237 = -qJ(1) + t325;
t443 = cos(t237) + cos(t234);
t445 = t273 / 0.2e1 + t443 / 0.4e1;
t260 = pkin(5) * t301 + pkin(1);
t378 = t307 * t308;
t437 = pkin(5) * t378 - t260 * t302;
t444 = t437 ^ 2;
t252 = t369 - t289;
t285 = pkin(18) - pkin(22);
t267 = -qJ(2) + t285;
t374 = pkin(21) - atan2(cos(t267), -sin(t267));
t129 = -atan2(-sin(t252), cos(t252)) + t374;
t442 = cos(t129) * (rSges(11,1) * t308 - rSges(11,2) * t302) + sin(t129) * (rSges(11,1) * t302 + rSges(11,2) * t308);
t441 = -Icges(11,5) * t45 + Icges(4,5) * t268 + Icges(9,5) * t143 - Icges(11,6) * t44 + Icges(4,6) * t272 - Icges(9,6) * t142;
t101 = -Icges(9,2) * t142 + t389;
t102 = Icges(9,1) * t143 - t390;
t180 = Icges(4,2) * t272 + t393;
t181 = Icges(4,1) * t268 + t392;
t31 = -Icges(11,2) * t44 - t396;
t32 = -Icges(11,1) * t45 - t397;
t440 = -t101 * t143 - t102 * t142 - t180 * t268 + t181 * t272 + t31 * t45 - t32 * t44;
t208 = sin(t234);
t211 = sin(t237);
t269 = sin(t290);
t439 = t269 / 0.2e1 - t211 / 0.4e1 + t208 / 0.4e1;
t438 = Icges(7,4) * t368;
t128 = -qJ(2) + t129;
t78 = rSges(11,2) * cos(t128) - rSges(11,1) * sin(t128);
t436 = rSges(5,1) * cos(t256) - rSges(5,2) * sin(t256) - pkin(15) - t437;
t288 = t309 ^ 2;
t435 = t450 * t288 + (t446 * t303 + (-t448 + t451) * t309) * t303;
t434 = t449 * t303 + t447 * t309;
t433 = -t447 * t303 + t449 * t309;
t406 = rSges(4,2) * t268;
t407 = rSges(4,1) * t272;
t432 = (-t406 + t407) * t303 - rSges(4,3) * t309;
t235 = qJ(1) + t325;
t209 = sin(t235);
t236 = -qJ(1) + t324;
t210 = sin(t236);
t291 = qJ(1) - qJ(4);
t270 = sin(t291);
t362 = t270 / 0.2e1 + t210 / 0.4e1 - t209 / 0.4e1;
t213 = cos(t235);
t214 = cos(t236);
t274 = cos(t291);
t363 = t214 / 0.4e1 - t274 / 0.2e1 + t213 / 0.4e1;
t431 = -t363 * rSges(6,1) + t362 * rSges(6,2);
t430 = -t362 * rSges(6,1) - t363 * rSges(6,2);
t323 = -t445 * rSges(6,1) + t439 * rSges(6,2);
t66 = t323 - t431;
t429 = 0.2e1 * t66;
t322 = t439 * rSges(6,1) + t445 * rSges(6,2);
t67 = t322 - t430;
t428 = 0.2e1 * t67;
t286 = t303 ^ 2;
t427 = -pkin(1) / 0.2e1;
t426 = m(6) / 0.2e1;
t247 = -qJ(1) + t256;
t230 = sin(t247);
t421 = -t230 / 0.2e1;
t246 = qJ(1) + t256;
t231 = cos(t246);
t420 = -t231 / 0.2e1;
t418 = t303 / 0.2e1;
t417 = -t309 / 0.2e1;
t293 = qJ(1) - qJ(2);
t416 = pkin(1) * sin(t293);
t415 = pkin(1) * cos(t293);
t414 = pkin(1) * t302;
t413 = pkin(1) * t308;
t412 = pkin(2) * t143;
t202 = -t301 * t302 + t378;
t411 = pkin(5) * t202;
t277 = qJ(1) - t289;
t410 = pkin(5) * sin(t277);
t409 = pkin(5) * cos(t277);
t295 = sin(pkin(20));
t298 = cos(pkin(20));
t192 = -t295 * t310 + t298 * t304;
t196 = t295 * t304 + t298 * t310;
t130 = t192 * t307 - t196 * t301;
t131 = t192 * t301 + t196 * t307;
t332 = t302 * t130 + t131 * t308;
t333 = t130 * t308 - t302 * t131;
t39 = atan2(t263 * t332 - t333 * t264, -t333 * t263 - t332 * t264) + t289;
t37 = sin(t39);
t229 = sin(t246);
t232 = cos(t247);
t402 = t303 * rSges(6,1);
t261 = -t402 / 0.2e1;
t280 = t303 * rSges(6,2);
t282 = t309 * rSges(6,2);
t400 = t309 * rSges(6,1);
t373 = t400 / 0.2e1;
t41 = (t282 + t402) * t274 + (t280 - t400) * t270 + (t282 - t402) * t273 + t269 * (t280 + t400) + (-t210 + t209) * (t373 - t280 / 0.2e1) + (-t211 + t208) * (t373 + t280 / 0.2e1) + (t214 + t213) * (-t282 / 0.2e1 + t261) + t443 * (t282 / 0.2e1 + t261) + ((-t231 + t232) * t309 + (-t229 - t230) * t303) * rSges(6,3);
t408 = t37 * t41;
t404 = rSges(7,3) * t309;
t300 = sin(qJ(4));
t306 = cos(qJ(4));
t38 = cos(t39);
t20 = -Icges(6,6) * t38 + (Icges(6,4) * t306 - Icges(6,2) * t300) * t37;
t403 = t20 * t300;
t401 = t303 * t37;
t399 = t309 * t37;
t395 = Icges(3,4) * t302;
t394 = Icges(3,4) * t308;
t391 = Icges(7,4) * t137;
t388 = Icges(10,4) * t112;
t387 = Icges(10,4) * t113;
t384 = t300 * t303;
t383 = t300 * t309;
t382 = t302 * t307;
t381 = t303 * t306;
t380 = t303 * t309;
t379 = t306 * t309;
t227 = -rSges(7,1) * t304 + rSges(7,2) * t310;
t228 = rSges(7,1) * t310 + rSges(7,2) * t304;
t148 = -t227 * t311 - t228 * t305;
t326 = t227 * t305 - t228 * t311;
t110 = t148 * t308 + t326 * t302;
t377 = t303 * rSges(7,3) + t110 * t309;
t75 = t442 * t303;
t76 = t442 * t309;
t376 = t303 * rSges(4,3) + t309 * t407;
t116 = t303 * t432 + t309 * (-t309 * t406 + t376);
t276 = qJ(1) + t289;
t240 = -pkin(5) * sin(t276) / 0.2e1;
t171 = t410 / 0.2e1 + t240;
t242 = pkin(5) * cos(t276) / 0.2e1;
t172 = -t409 / 0.2e1 + t242;
t375 = t286 + t288;
t103 = -rSges(9,1) * t142 - rSges(9,2) * t143;
t359 = -rSges(9,1) * t143 + rSges(9,2) * t142;
t95 = t359 * t303;
t96 = t359 * t309;
t372 = m(9) * (t103 ^ 2 + t95 ^ 2 + t96 ^ 2) + (t451 * t286 + ((-t446 + t450) * t303 + t448 * t309) * t309) * t303;
t190 = rSges(4,1) * t268 + rSges(4,2) * t272;
t371 = -t190 - t413;
t370 = 0.2e1 * t426;
t222 = rSges(3,1) * t302 + rSges(3,2) * t308;
t265 = sin(t285);
t266 = cos(t285);
t360 = rSges(8,1) * t266 - rSges(8,2) * t265;
t351 = -Icges(3,1) * t302 - t394;
t349 = Icges(7,1) * t368 - t391;
t347 = Icges(10,1) * t112 + t387;
t345 = -Icges(3,2) * t308 - t395;
t343 = -Icges(7,2) * t137 + t438;
t341 = Icges(10,2) * t113 + t388;
t121 = t309 * rSges(5,3) + t436 * t303;
t122 = t303 * rSges(5,3) - t436 * t309;
t321 = m(5) * (-t121 * t309 - t122 * t303);
t58 = -pkin(2) * t142 + rSges(10,1) * t112 + rSges(10,2) * t113;
t319 = rSges(10,1) * t113 - rSges(10,2) * t112 - t412;
t318 = pkin(4) * sin(qJ(2) - t374) - t78;
t317 = t435 * t309 + t372;
t79 = rSges(9,3) * t309 + (-pkin(15) - t103) * t303;
t281 = t309 * pkin(15);
t80 = rSges(9,3) * t303 + t103 * t309 + t281;
t316 = m(9) * (t79 * t96 + t80 * t95) + (-t142 * t90 + t143 * t92 + t157 * t272 + t159 * t268 - t27 * t44 - t29 * t45 + t441 * t303 + t440 * t309) * t418 + (-t142 * t89 + t143 * t91 + t156 * t272 + t158 * t268 - t26 * t44 - t28 * t45 + t440 * t303 - t441 * t309) * t417;
t313 = pkin(5) ^ 2;
t312 = pkin(11) + rSges(6,3);
t292 = qJ(1) + qJ(2);
t257 = -pkin(15) + t414;
t254 = cos(t292) * t427;
t253 = sin(t292) * t427;
t233 = t257 * t303;
t226 = rSges(2,1) * t309 - rSges(2,2) * t303;
t225 = rSges(3,1) * t308 - rSges(3,2) * t302;
t223 = -rSges(2,1) * t303 - rSges(2,2) * t309;
t198 = t301 * t308 + t382;
t177 = t313 * t202 ^ 2;
t173 = pkin(5) * t382 + t260 * t308;
t160 = rSges(3,3) * t309 + (-pkin(15) + t222) * t303;
t153 = rSges(3,3) * t303 - t222 * t309 + t281;
t152 = t371 * t309;
t151 = t371 * t303;
t145 = -t416 / 0.2e1 + t253 + t172;
t144 = t254 - t415 / 0.2e1 + t171;
t140 = t233 - t432;
t139 = (-t257 - t406) * t309 + t376;
t133 = rSges(8,3) * t309 + t360 * t303 + t233;
t132 = rSges(8,3) * t303 + (-t257 - t360) * t309;
t115 = t116 - t414;
t111 = t148 * t302 - t326 * t308;
t94 = -pkin(14) * t309 + t377;
t93 = t404 + (pkin(14) - t110) * t303;
t77 = t78 - t414;
t74 = -t309 * t413 + t76;
t73 = -t303 * t413 + t75;
t69 = rSges(11,3) * t303 + (-t257 - t318) * t309;
t68 = rSges(11,3) * t309 + t318 * t303 + t233;
t59 = t303 * (t110 * t303 - t404) + t309 * t377;
t57 = t319 * t309;
t56 = t319 * t303;
t55 = rSges(10,3) * t303 + t58 * t309 + t281;
t54 = rSges(10,3) * t309 + (-pkin(15) - t58) * t303;
t51 = -t410 / 0.2e1 + t240 + t254 + t415 / 0.2e1 - t303 * pkin(15) + (t420 + t232 / 0.2e1) * t312 + (t229 / 0.2e1 + t421) * pkin(9) + t322 + t430;
t50 = t409 / 0.2e1 + t242 + t416 / 0.2e1 + t253 + t281 + (-t229 / 0.2e1 + t421) * t312 + (t420 - t232 / 0.2e1) * pkin(9) + t323 + t431;
t36 = t38 * t379 + t384;
t35 = -t38 * t383 + t381;
t34 = t38 * t381 - t383;
t33 = -t38 * t384 - t379;
t21 = -Icges(6,5) * t38 + (Icges(6,1) * t306 - Icges(6,4) * t300) * t37;
t19 = -Icges(6,3) * t38 + (Icges(6,5) * t306 - Icges(6,6) * t300) * t37;
t18 = Icges(6,1) * t36 + Icges(6,4) * t35 + Icges(6,5) * t399;
t17 = Icges(6,1) * t34 + Icges(6,4) * t33 + Icges(6,5) * t401;
t16 = Icges(6,4) * t36 + Icges(6,2) * t35 + Icges(6,6) * t399;
t15 = Icges(6,4) * t34 + Icges(6,2) * t33 + Icges(6,6) * t401;
t14 = Icges(6,5) * t36 + Icges(6,6) * t35 + Icges(6,3) * t399;
t13 = Icges(6,5) * t34 + Icges(6,6) * t33 + Icges(6,3) * t401;
t12 = t37 * t306 * t21;
t5 = -t38 * t19 - t37 * t403 + t12;
t4 = t19 * t399 + t20 * t35 + t21 * t36;
t3 = t19 * t401 + t20 * t33 + t21 * t34;
t2 = -t14 * t38 + (-t16 * t300 + t18 * t306) * t37;
t1 = -t13 * t38 + (-t15 * t300 + t17 * t306) * t37;
t6 = [m(2) * (t223 ^ 2 + t226 ^ 2) + m(3) * (t153 ^ 2 + t160 ^ 2) + m(4) * (t139 ^ 2 + t140 ^ 2) + m(8) * (t132 ^ 2 + t133 ^ 2) + m(5) * (t121 ^ 2 + t122 ^ 2) + m(7) * (t93 ^ 2 + t94 ^ 2) + m(9) * (t79 ^ 2 + t80 ^ 2) + m(11) * (t68 ^ 2 + t69 ^ 2) + m(10) * (t54 ^ 2 + t55 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2) + t368 * (Icges(7,2) * t368 + t391) + (Icges(5,2) * t38 - t19) * t38 + (Icges(8,1) * t265 + 0.2e1 * Icges(8,4) * t266) * t265 + t266 ^ 2 * Icges(8,2) + t137 * (Icges(7,1) * t137 + t438) + t12 + t268 * t181 + t272 * t180 + t143 * t102 - t142 * t101 - t44 * t31 - t45 * t32 + (Icges(5,1) * t37 + 0.2e1 * Icges(5,4) * t38 - t403) * t37 - t302 * (-Icges(3,2) * t302 + t394) + t308 * (Icges(3,1) * t308 - t395) - t113 * (-Icges(10,1) * t113 + t388) + t112 * (Icges(10,2) * t112 - t387) + Icges(2,3); t316 + m(8) * (-t132 * t303 - t133 * t309) * t413 + m(7) * (-t303 * t94 - t309 * t93) * t111 + t173 * t321 + m(6) * (t144 * t51 + t145 * t50) + m(4) * (t139 * t151 + t140 * t152) + m(11) * (t68 * t74 + t69 * t73) + m(10) * (t54 * t57 + t55 * t56) + m(3) * (-t153 * t303 - t160 * t309) * t225 + (t112 * (Icges(10,6) * t303 + t341 * t309) - t113 * (Icges(10,5) * t303 + t347 * t309) + t137 * (Icges(7,5) * t303 + t349 * t309) - (Icges(3,6) * t303 + t345 * t309) * t302 + (Icges(3,5) * t303 + t351 * t309) * t308 + t368 * (Icges(7,6) * t303 + t343 * t309)) * t418 + (t112 * (-Icges(10,6) * t309 + t341 * t303) - t113 * (-Icges(10,5) * t309 + t347 * t303) + t137 * (-Icges(7,5) * t309 + t349 * t303) - (-Icges(3,6) * t309 + t345 * t303) * t302 + (-Icges(3,5) * t309 + t351 * t303) * t308 + t368 * (-Icges(7,6) * t309 + t343 * t303)) * t417 + (Icges(3,5) * t308 + Icges(7,5) * t137 - Icges(10,5) * t113 - Icges(3,6) * t302 + Icges(7,6) * t368 + Icges(10,6) * t112) * (t286 / 0.2e1 + t288 / 0.2e1); (t144 ^ 2 + t145 ^ 2 + t444) * t370 + m(10) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(4) * (t115 ^ 2 + t151 ^ 2 + t152 ^ 2) + m(11) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) + t372 + m(8) * (t375 * t308 ^ 2 + t302 ^ 2) * pkin(1) ^ 2 + m(3) * (t375 * t225 ^ 2 + t222 ^ 2) + m(7) * (t375 * t111 ^ 2 + t59 ^ 2) + m(5) * (t375 * t173 ^ 2 + t444) + t434 * t303 * t286 + (t433 * t288 + (t433 * t303 + t434 * t309) * t303 + t435) * t309; t316 + m(11) * (t68 * t76 + t69 * t75) + m(6) * (t171 * t51 + t172 * t50) + t198 * pkin(5) * t321 + m(10) * (-t303 * t55 - t309 * t54) * t412 + m(4) * (-t139 * t303 - t140 * t309) * t190; (t144 * t171 + t145 * t172 + t411 * t437) * t370 + m(11) * (t73 * t75 + t74 * t76 + t77 * t78) + m(4) * (t115 * t116 + (-t151 * t303 - t152 * t309) * t190) + m(5) * (t375 * t198 * t173 + t202 * t437) * pkin(5) + m(10) * (-t142 * t58 + (-t303 * t56 - t309 * t57) * t143) * pkin(2) + t317; m(4) * (t375 * t190 ^ 2 + t116 ^ 2) + m(6) * (t171 ^ 2 + t172 ^ 2 + t177) + m(5) * (t375 * t313 * t198 ^ 2 + t177) + m(11) * (t75 ^ 2 + t76 ^ 2 + t78 ^ 2) + m(10) * (t375 * t143 ^ 2 + t142 ^ 2) * pkin(2) ^ 2 + t317; m(6) * (t50 * t66 + t51 * t67) - t5 * t38 + ((t2 / 0.2e1 + t4 / 0.2e1) * t309 + (t3 / 0.2e1 + t1 / 0.2e1) * t303) * t37; (t144 * t428 + t145 * t429 - t408 * t437) * t426; (t171 * t428 + t172 * t429 - t408 * t411) * t426; t38 ^ 2 * t5 + (t66 ^ 2 + t67 ^ 2) * t370 + ((t303 * ((t16 * t33 + t18 * t34) * t309 + (t15 * t33 + t17 * t34) * t303) + m(6) * t41 ^ 2 / 0.4e1 + t309 * ((t16 * t35 + t18 * t36) * t309 + (t15 * t35 + t17 * t36) * t303) + (t303 * (t13 * t286 + t14 * t380) + t309 * (t13 * t380 + t14 * t288)) * t37) * t37 + ((-t2 - t4) * t309 + (-t1 - t3) * t303) * t38) * t37;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
