% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m1DE_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:26
% DurationCPUTime: 13.35s
% Computational Cost: add. (6134->515), mult. (10314->749), div. (0->0), fcn. (8068->10), ass. (0->321)
t232 = sin(qJ(1));
t231 = sin(qJ(2));
t234 = cos(qJ(3));
t230 = sin(qJ(3));
t235 = cos(qJ(2));
t373 = t230 * t235;
t176 = t231 * t234 + t373;
t405 = pkin(3) * qJD(3);
t340 = t176 * t405;
t146 = t232 * t340;
t214 = pkin(3) * t234 + pkin(2);
t148 = pkin(3) * t373 + t214 * t231;
t349 = qJD(2) * t232;
t361 = t148 * t349 + t146;
t374 = t230 * t231;
t147 = -pkin(3) * t374 + t214 * t235;
t291 = pkin(1) + t147;
t141 = rSges(5,1) + t291;
t236 = cos(qJ(1));
t100 = -rSges(5,3) * t232 + t141 * t236;
t453 = qJD(1) * t100;
t53 = -t361 + t453;
t227 = qJD(1) + qJD(4);
t196 = rSges(6,1) * t236 - rSges(6,2) * t232;
t198 = rSges(6,1) * t232 + rSges(6,2) * t236;
t229 = sin(qJ(4));
t233 = cos(qJ(4));
t307 = -t196 * t233 + t198 * t229;
t454 = t227 * t307;
t226 = qJD(2) + qJD(3);
t228 = qJ(2) + qJ(3);
t219 = sin(t228);
t220 = cos(t228);
t296 = rSges(4,1) * t219 + rSges(4,2) * t220;
t439 = t296 * t226;
t101 = -rSges(5,3) * t236 - t141 * t232;
t452 = qJD(1) * t101;
t348 = qJD(2) * t235;
t448 = pkin(3) * (qJD(3) * t176 + t230 * t348);
t377 = t220 * t232;
t380 = t219 * t232;
t392 = Icges(4,6) * t236;
t112 = Icges(4,4) * t377 - Icges(4,2) * t380 + t392;
t396 = Icges(4,4) * t220;
t165 = -Icges(4,1) * t219 - t396;
t447 = -t165 * t232 + t112;
t287 = -Icges(4,2) * t219 + t396;
t113 = -Icges(4,6) * t232 + t236 * t287;
t446 = -t165 * t236 + t113;
t397 = Icges(4,4) * t219;
t289 = Icges(4,1) * t220 - t397;
t115 = -Icges(4,5) * t232 + t236 * t289;
t163 = -Icges(4,2) * t220 - t397;
t445 = -t163 * t236 - t115;
t444 = t163 + t289;
t443 = t287 - t165;
t192 = -t229 * rSges(6,1) - t233 * rSges(6,2);
t199 = rSges(6,1) * t233 - rSges(6,2) * t229;
t278 = t192 * t236 - t199 * t232;
t442 = 0.2e1 * qJD(2);
t404 = t231 * rSges(3,2);
t237 = qJD(2) ^ 2;
t417 = pkin(2) * t231;
t217 = t237 * t417;
t149 = pkin(3) * t219 * t226 ^ 2 + t217;
t416 = pkin(2) * t235;
t441 = t149 * (-pkin(3) * t220 - t416);
t177 = t234 * t235 - t374;
t339 = t177 * t405;
t273 = -qJD(2) * t147 - t339;
t352 = qJD(1) * t232;
t437 = t148 * t352 + t273 * t236;
t436 = qJD(2) * t177;
t186 = -Icges(3,5) * t231 - Icges(3,6) * t235;
t261 = qJD(2) * t186;
t399 = Icges(3,4) * t231;
t290 = Icges(3,1) * t235 - t399;
t131 = -Icges(3,5) * t232 + t236 * t290;
t398 = Icges(3,4) * t235;
t288 = -Icges(3,2) * t231 + t398;
t129 = -Icges(3,6) * t232 + t236 * t288;
t388 = t129 * t231;
t281 = -t131 * t235 + t388;
t286 = Icges(3,5) * t235 - Icges(3,6) * t231;
t391 = Icges(3,3) * t236;
t435 = -t236 * t261 + (t232 * t286 + t281 + t391) * qJD(1);
t370 = t232 * t235;
t372 = t231 * t232;
t393 = Icges(3,6) * t236;
t128 = Icges(3,4) * t370 - Icges(3,2) * t372 + t393;
t211 = Icges(3,4) * t372;
t395 = Icges(3,5) * t236;
t130 = Icges(3,1) * t370 - t211 + t395;
t282 = t128 * t231 - t130 * t235;
t127 = -Icges(3,3) * t232 + t236 * t286;
t354 = qJD(1) * t127;
t434 = qJD(1) * t282 - t232 * t261 - t354;
t161 = -Icges(4,5) * t219 - Icges(4,6) * t220;
t136 = t161 * t236;
t285 = Icges(4,5) * t220 - Icges(4,6) * t219;
t389 = t113 * t219;
t390 = Icges(4,3) * t236;
t433 = -t226 * t136 + (-t115 * t220 + t232 * t285 + t389 + t390) * qJD(1);
t135 = t161 * t232;
t202 = Icges(4,4) * t380;
t394 = Icges(4,5) * t236;
t114 = Icges(4,1) * t377 - t202 + t394;
t284 = t112 * t219 - t114 * t220;
t111 = -Icges(4,3) * t232 + t236 * t285;
t355 = qJD(1) * t111;
t432 = qJD(1) * t284 - t135 * t226 - t355;
t280 = t163 * t219 - t165 * t220;
t431 = qJD(1) * t280 + t285 * t226;
t188 = -Icges(3,2) * t235 - t399;
t190 = -Icges(3,1) * t231 - t398;
t279 = t188 * t231 - t235 * t190;
t430 = t279 * qJD(1) + t286 * qJD(2);
t428 = -t128 * t236 + t129 * t232;
t185 = t226 * t236;
t375 = t226 * t232;
t427 = qJD(1) * t444 + t185 * t447 - t375 * t446;
t204 = rSges(4,2) * t380;
t342 = rSges(4,1) * t377;
t116 = rSges(4,3) * t236 - t204 + t342;
t376 = t220 * t236;
t205 = rSges(4,1) * t376;
t379 = t219 * t236;
t408 = rSges(4,3) * t232;
t117 = -rSges(4,2) * t379 + t205 - t408;
t143 = t296 * t232;
t144 = t296 * t236;
t341 = pkin(2) * t348;
t183 = qJD(1) * t204;
t74 = -qJD(1) * t342 + t183 + (-rSges(4,3) * qJD(1) - t439) * t236;
t410 = rSges(4,2) * t219;
t411 = rSges(4,1) * t220;
t297 = -t410 + t411;
t75 = -t226 * t143 + (t236 * t297 - t408) * qJD(1);
t426 = (-t232 * t75 + t117 * t352 + (-qJD(1) * t116 - t74) * t236 - t375 * t143 - t144 * t185) * (-t116 * t375 - t117 * t185 - t341);
t305 = qJD(1) * t226;
t169 = t232 * t305;
t425 = -t169 / 0.2e1;
t170 = t236 * t305;
t424 = -t170 / 0.2e1;
t423 = t375 / 0.2e1;
t422 = -t375 / 0.2e1;
t421 = -t185 / 0.2e1;
t420 = t185 / 0.2e1;
t419 = -t232 / 0.2e1;
t418 = t236 / 0.2e1;
t415 = -qJD(1) / 0.2e1;
t414 = qJD(1) / 0.2e1;
t412 = rSges(3,1) * t235;
t409 = rSges(3,3) * t236;
t403 = t232 * rSges(3,3);
t400 = t236 * t111 + t115 * t377;
t145 = pkin(4) + t291;
t387 = t145 * t232;
t386 = t145 * t236;
t385 = t149 * t220;
t152 = t186 * t232;
t153 = t186 * t236;
t215 = pkin(1) + t416;
t238 = qJD(1) ^ 2;
t382 = t215 * t238;
t371 = t231 * t236;
t369 = t235 * t236;
t368 = t235 * t237;
t118 = t236 * t127;
t61 = -t232 * t280 + t136;
t367 = t61 * qJD(1);
t79 = -t232 * t279 + t153;
t366 = t79 * qJD(1);
t365 = t131 * t370 + t118;
t345 = qJD(1) * qJD(2);
t332 = t232 * t345;
t362 = qJD(1) * t146 + t148 * t332;
t357 = t188 + t290;
t356 = t288 - t190;
t353 = qJD(1) * t286;
t351 = qJD(1) * t236;
t350 = qJD(2) * t231;
t347 = qJD(2) * t236;
t335 = t231 * t349;
t334 = t231 * t347;
t333 = -pkin(1) - t412;
t331 = -t352 / 0.2e1;
t330 = -t351 / 0.2e1;
t328 = t349 / 0.2e1;
t327 = -t347 / 0.2e1;
t212 = rSges(3,2) * t372;
t276 = -rSges(3,1) * t370 - t409;
t325 = -pkin(1) * t232 + t212 + t276;
t320 = rSges(3,1) * t369 - t403;
t324 = -rSges(3,2) * t371 + pkin(1) * t236 + t320;
t323 = -t296 - t417;
t321 = -t215 - t411;
t319 = -qJD(1) * t115 + t226 * t447;
t318 = -(-t232 * t289 - t394) * qJD(1) + t446 * t226;
t317 = -qJD(1) * t113 - t114 * t226 - t163 * t375;
t316 = -(-t232 * t287 - t392) * qJD(1) + t445 * t226;
t110 = Icges(4,5) * t377 - Icges(4,6) * t380 + t390;
t314 = t110 + t389;
t313 = t443 * t226;
t312 = t444 * t226;
t311 = -qJD(1) * t144 - t297 * t375;
t126 = Icges(3,5) * t370 - Icges(3,6) * t372 + t391;
t310 = -t126 - t388;
t309 = t278 * qJD(4);
t104 = t192 * t232 + t199 * t236;
t306 = qJD(1) * t143 - t185 * t297;
t125 = t297 * t226;
t302 = -t125 - t341;
t301 = -t112 * t379 + t114 * t376;
t300 = -t128 * t371 + t130 * t369;
t193 = -t404 + t412;
t298 = t231 * rSges(3,1) + t235 * rSges(3,2);
t103 = t196 * t229 + t198 * t233;
t272 = -qJD(2) * t148 - t340;
t258 = t272 * t236;
t40 = -t103 * t227 - t145 * t352 + t258;
t277 = -t145 * t351 + t361;
t41 = t104 * t227 - t277;
t295 = -t232 * t41 - t236 * t40;
t54 = t258 + t452;
t294 = -t232 * t53 - t236 * t54;
t179 = rSges(3,1) * t348 - rSges(3,2) * t350;
t292 = -pkin(1) * t238 - qJD(2) * t179;
t65 = -t112 * t220 - t114 * t219;
t77 = -t128 * t235 - t130 * t231;
t78 = -t129 * t235 - t131 * t231;
t260 = t177 * qJD(3);
t91 = t214 * t348 + (-t230 * t350 + t260) * pkin(3);
t96 = t260 + t436;
t275 = -qJD(2) * t91 - t405 * t96;
t274 = 0.2e1 * t298;
t271 = -pkin(3) * t436 - t339;
t270 = -t111 + t284;
t48 = t126 * t236 - t232 * t282;
t49 = -t129 * t372 + t365;
t267 = (-t232 * t49 + t236 * t48) * qJD(2);
t50 = -t126 * t232 + t300;
t108 = t131 * t369;
t51 = -t127 * t232 - t129 * t371 + t108;
t266 = (-t232 * t51 + t236 * t50) * qJD(2);
t95 = t115 * t376;
t265 = -t236 * t314 + t95;
t263 = qJD(2) * t190;
t262 = qJD(2) * t188;
t257 = -qJD(1) * t285 + t135 * t185 - t136 * t375;
t256 = -(-t188 * t236 - t131) * t232 + (Icges(3,2) * t370 - t130 + t211) * t236;
t255 = -t145 * t238 + t275;
t246 = -qJD(1) * t110 + t219 * t317 - t220 * t319;
t10 = t246 * t232 - t236 * t432;
t245 = t219 * t316 - t220 * t318 - t355;
t11 = t245 * t232 - t236 * t433;
t42 = t110 * t236 - t232 * t284;
t43 = -t113 * t380 + t400;
t19 = t185 * t42 - t375 * t43 + t367;
t44 = -t110 * t232 + t301;
t45 = -t111 * t232 - t113 * t379 + t95;
t62 = -t236 * t280 - t135;
t59 = t62 * qJD(1);
t20 = t185 * t44 - t375 * t45 + t59;
t249 = -t445 * t375 + (Icges(4,2) * t377 - t114 + t202) * t185 + t443 * qJD(1);
t239 = t249 * t219 - t220 * t427;
t244 = -qJD(1) * t161 + t219 * t313 - t220 * t312;
t28 = t232 * t431 + t244 * t236;
t29 = t244 * t232 - t236 * t431;
t32 = t219 * t319 + t220 * t317;
t33 = t219 * t318 + t220 * t316;
t66 = -t113 * t220 - t115 * t219;
t8 = t232 * t432 + t246 * t236;
t9 = t232 * t433 + t245 * t236;
t254 = (qJD(1) * t28 - t169 * t44 - t170 * t45 + t185 * t8 - t375 * t9) * t419 + (t219 * t427 + t249 * t220) * t415 + t19 * t331 + t20 * t330 + (qJD(1) * t29 + t10 * t185 - t11 * t375 - t169 * t42 - t170 * t43) * t418 + (-t232 * t43 + t236 * t42) * t425 + (-t232 * t45 + t236 * t44) * t424 + (-t232 * t9 + t236 * t8 + (-t232 * t44 - t236 * t45) * qJD(1)) * t422 + (t10 * t236 - t11 * t232 + (-t232 * t42 - t236 * t43) * qJD(1)) * t420 + (-t232 * t33 + t236 * t32 + (-t232 * t65 - t236 * t66) * qJD(1)) * t414 + (-t232 * t257 + t236 * t239) * t423 + (t232 * t239 + t236 * t257) * t421;
t253 = (t231 * t356 - t235 * t357) * qJD(1);
t158 = t176 * pkin(3);
t252 = t158 * t352 + t236 * t271;
t63 = qJD(1) * t278 + t309;
t90 = -t214 * t350 - t448;
t21 = t227 * t63 + t255 * t232 + (t272 + t90) * t351;
t22 = t227 * t454 + t236 * t255 - t352 * t90 + t362;
t251 = -t21 * t232 - t22 * t236 + (t232 * t40 - t236 * t41) * qJD(1);
t55 = t214 * t335 + t232 * t448 - t453;
t30 = qJD(1) * t55 + t236 * t275 + t362;
t89 = t90 * t236;
t56 = t89 + t452;
t31 = t275 * t232 + (t56 + t258) * qJD(1);
t250 = -t232 * t31 - t236 * t30 + (t232 * t54 - t236 * t53) * qJD(1);
t57 = -pkin(2) * t334 - t185 * t296 + (-t215 * t232 - t116) * qJD(1);
t86 = qJD(1) * t129 + t232 * t262;
t88 = qJD(1) * t131 + t232 * t263;
t243 = -qJD(1) * t126 + qJD(2) * t77 - t231 * t86 + t235 * t88;
t85 = t236 * t262 + (-t232 * t288 - t393) * qJD(1);
t87 = t236 * t263 + (-t232 * t290 - t395) * qJD(1);
t242 = qJD(2) * t78 - t231 * t85 + t235 * t87 - t354;
t174 = t288 * qJD(2);
t175 = t290 * qJD(2);
t241 = -qJD(1) * t186 + t174 * t231 - t175 * t235 + (-t188 * t235 - t190 * t231) * qJD(2);
t240 = t256 * t231 + t235 * t428;
t207 = pkin(2) * t335;
t206 = t352 * t404;
t82 = t324 * qJD(1) - t298 * t349;
t81 = t325 * qJD(1) - t298 * t347;
t80 = -t236 * t279 - t152;
t76 = t80 * qJD(1);
t67 = -t158 * t351 + t232 * t271;
t60 = -t148 * t351 + t232 * t273;
t58 = -t296 * t375 - t207 + (t215 * t236 + t117) * qJD(1);
t47 = t292 * t236 + ((-t193 * t236 + t403) * qJD(1) + t274 * t349) * qJD(1);
t46 = t292 * t232 + (qJD(1) * t276 - t274 * t347 + t206) * qJD(1);
t39 = t241 * t232 - t236 * t430;
t38 = t232 * t430 + t241 * t236;
t37 = qJD(2) * t281 - t231 * t87 - t235 * t85;
t36 = qJD(2) * t282 - t231 * t88 - t235 * t86;
t35 = -t236 * t382 - qJD(1) * t75 - t125 * t185 + t296 * t169 + (0.2e1 * t231 * t332 - t236 * t368) * pkin(2);
t34 = -t232 * t382 + qJD(1) * t74 - t125 * t375 - t296 * t170 + (-0.2e1 * qJD(1) * t334 - t232 * t368) * pkin(2);
t25 = -t116 * t170 + t117 * t169 - t185 * t74 - t375 * t75 + t217;
t24 = t76 + t266;
t23 = t267 + t366;
t1 = [(qJD(2) * t279 + t174 * t235 + t175 * t231) * qJD(1) + (t59 + (t301 - t43 + t400) * t185 - (t42 + t265) * t375 + (-t185 * t314 - t270 * t375) * t232) * t421 + m(4) * (t58 * t183 + t35 * t204 + t34 * t205 + t57 * t207 + (t35 * t321 - t34 * rSges(4,3) + t57 * t439 + (t57 * rSges(4,3) + t321 * t58) * qJD(1)) * t232 + (-t35 * rSges(4,3) + t34 * (t215 - t410) + t58 * (-pkin(2) * t350 - t439) + (t57 * (-t215 - t297) - t58 * rSges(4,3)) * qJD(1)) * t236) + (t76 + ((t300 - t49 + t365) * t236 + (-t108 - t48 + (t127 - t282) * t232) * t232) * qJD(2)) * t327 + (t219 * t312 + t220 * t313) * qJD(1) + (t65 + t61) * t425 + (t66 + t62) * t424 + (-t367 + (-t45 + t265) * t185 - (t236 * t284 - t400 + t44) * t375 + (t185 * t270 - t314 * t375) * t232 + t19) * t423 + (t33 + t28) * t422 - (t37 + t38) * t349 / 0.2e1 + (-t366 + ((t236 * t310 + t108 - t51) * t236 + (t232 * t310 - t118 + t365 - t50) * t232) * qJD(2) + t23) * t328 + (t22 * (-t103 - t387) + t21 * (t104 + t386) + (t309 + t89 + (t278 - t387) * qJD(1) - (-qJD(1) * t145 - t199 * t227) * t232 - (t192 * t227 + t272) * t236) * t41 + (-t386 * qJD(1) - t232 * t90 - t277) * t40) * m(6) + (t100 * t31 + t101 * t30 + (t55 + t53) * t54 + (t56 + t141 * t352 - (-rSges(5,3) * qJD(1) + t272) * t236) * t53) * m(5) + (t47 * (t232 * t333 + t212 - t409) + t46 * ((pkin(1) - t404) * t236 + t320) + t82 * t206 + (t324 * t81 - t325 * t82 + (t81 * rSges(3,3) + t333 * t82) * t232 + (t81 * (-pkin(1) - t193) - t82 * rSges(3,3)) * t236) * qJD(1)) * m(3) + (t32 + t29 + t20) * t420 + (t36 + t39 + t24) * t347 / 0.2e1 - ((t77 + t79) * t232 + (t78 + t80) * t236) * t345 / 0.2e1; ((t231 * t357 + t235 * t356) * qJD(1) + (-t231 * t428 + t256 * t235) * qJD(2)) * t415 + (-t232 * t37 + t236 * t36 + (-t77 * t232 - t236 * t78) * qJD(1)) * t414 + t254 + ((t152 * t347 - t353) * t236 + (t253 + (-t236 * t153 + t240) * qJD(2)) * t232) * t327 + ((t153 * t349 + t353) * t232 + (t253 + (-t232 * t152 + t240) * qJD(2)) * t236) * t328 + (qJD(1) * t38 + ((t232 * t434 + t243 * t236) * t236 - (t232 * t435 + t242 * t236) * t232 + (-t50 * t232 - t51 * t236) * qJD(1)) * t442) * t419 + (qJD(1) * t39 + ((t243 * t232 - t236 * t434) * t236 - (t242 * t232 - t236 * t435) * t232 + (-t48 * t232 - t49 * t236) * qJD(1)) * t442) * t418 + (t267 + t23) * t331 + (t24 + t266) * t330 + (t251 * t148 + t295 * t91 - t40 * t437 - t41 * t60 + t441) * m(6) + (t250 * t148 + t294 * t91 - t437 * t54 - t53 * t60 + t441) * m(5) + (-t25 * t416 + (qJD(1) * t296 * t57 - t25 * t116 + t302 * t58 + t323 * t34) * t232 + (-t25 * t117 + t57 * t302 + (t58 * qJD(1) + t35) * t323) * t236 - t57 * (-t236 * t341 + t306) - t58 * ((-t231 * t351 - t232 * t348) * pkin(2) + t311) + t426) * m(4) + ((-t232 * t82 - t236 * t81) * t179 - (t46 * t232 + t47 * t236) * t298 + (-t298 * t237 + t81 * t347 + t82 * t349) * t193) * m(3); t254 + (t25 * (-t116 * t232 - t117 * t236) - (t232 * t58 + t236 * t57) * t125 - (t34 * t232 + t35 * t236 + (-t232 * t57 + t236 * t58) * qJD(1)) * t296 - t57 * t306 - t58 * t311 + t426) * m(4) + ((t176 * t251 + t295 * t96 - t385) * pkin(3) - t40 * t252 - t41 * t67) * m(6) + ((t176 * t250 + t294 * t96 - t385) * pkin(3) - t54 * t252 - t53 * t67) * m(5); (-t103 * t22 + t104 * t21 + t40 * t454 + t41 * t63 - (t41 * t278 + t40 * t307) * t227) * m(6);];
tauc = t1(:);
