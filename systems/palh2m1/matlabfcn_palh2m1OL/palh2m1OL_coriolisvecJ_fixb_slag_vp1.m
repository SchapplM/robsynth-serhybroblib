% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m1OL_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:21:04
% EndTime: 2020-05-03 00:21:19
% DurationCPUTime: 7.37s
% Computational Cost: add. (2333->422), mult. (3150->642), div. (0->0), fcn. (1050->70), ass. (0->325)
t237 = qJD(2) + qJD(3);
t252 = qJD(4) / 0.2e1;
t188 = t252 + t237;
t251 = -qJD(5) / 0.2e1;
t113 = t251 + t188;
t272 = 2 * qJ(2);
t393 = 2 * qJ(3) + t272;
t354 = qJ(5) + t393;
t181 = qJ(4) + t354;
t126 = sin(t181);
t353 = -qJ(5) + t393;
t182 = qJ(4) + t353;
t127 = sin(t182);
t268 = qJD(1) ^ 2;
t450 = t268 / 0.2e1;
t245 = qJ(3) + qJ(4);
t220 = qJ(2) + t245;
t149 = 2 * t220;
t264 = rSges(6,3) + pkin(6);
t273 = rSges(6,2) ^ 2;
t274 = rSges(6,1) ^ 2;
t350 = t274 / 0.2e1 + t273 / 0.2e1;
t515 = sin(t149) * ((-(pkin(4) + t264) * (-pkin(4) + t264) + t350) * m(6) + (rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - Icges(5,1) - Icges(6,1) / 0.2e1 + Icges(5,2) - Icges(6,2) / 0.2e1 + Icges(6,3));
t13 = t450 * t515;
t136 = cos(t181);
t137 = cos(t182);
t242 = -qJ(5) + qJ(2);
t217 = qJ(3) + t242;
t160 = sin(t217);
t247 = qJ(2) + qJ(5);
t218 = qJ(3) + t247;
t161 = sin(t218);
t168 = cos(t218);
t248 = qJ(2) + qJ(3);
t213 = cos(t248);
t378 = qJD(5) * qJD(1);
t358 = m(6) * t378;
t329 = rSges(6,2) * t358;
t302 = pkin(3) * t329;
t244 = qJ(4) - qJ(5);
t209 = cos(t244);
t233 = qJD(4) - qJD(5);
t467 = m(6) * pkin(3);
t374 = t233 * t467;
t314 = rSges(6,2) * t209 * t374;
t464 = m(5) * rSges(5,1);
t171 = pkin(4) * m(6) + t464;
t409 = t171 * sin(qJ(4));
t327 = qJD(4) * pkin(3) * t409;
t330 = rSges(6,1) * t358;
t446 = m(6) * t268;
t331 = rSges(6,1) * t446 / 0.2e1;
t371 = rSges(6,2) * t446;
t466 = m(4) * rSges(4,2);
t372 = t268 * t466;
t174 = 2 * t248;
t265 = rSges(4,1) * m(4);
t405 = (rSges(4,2) * t265 - Icges(4,4)) * cos(t174);
t426 = ((rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(4,1) + Icges(4,2) + (m(6) + m(5)) * pkin(3) ^ 2) * sin(t174);
t471 = -0.2e1 * t188;
t109 = m(5) * pkin(3) + t265 + t467;
t206 = sin(t248);
t508 = t206 * pkin(1) * t109;
t216 = qJ(4) + t393;
t159 = sin(t216);
t408 = t171 * t268;
t74 = pkin(3) * t159 * t408;
t516 = t13 + pkin(1) * t213 * t372 + t113 * t314 + t168 * t302 + t327 * t471 + t426 * t450 + t74 + (t508 + t405) * t268 + ((t136 / 0.2e1 - t137 / 0.2e1) * t371 + (t126 + t127) * t331 + (t160 + t161) * t330) * pkin(3);
t270 = 2 * qJ(4);
t179 = t270 + t353;
t225 = rSges(6,2) * t264;
t456 = pkin(4) * rSges(6,1);
t479 = t225 - t456;
t512 = t479 * m(6) - Icges(6,6);
t514 = sin(t179) * t512;
t257 = sin(qJ(2));
t509 = qJD(2) ^ 2;
t513 = t257 * t509;
t511 = t237 ^ 2;
t510 = -2 * Icges(6,4);
t266 = m(4) + m(5);
t462 = rSges(3,1) * m(3);
t468 = m(6) * pkin(2);
t507 = t257 * pkin(1) * (pkin(2) * t266 + t462 + t468);
t506 = Icges(6,1) - Icges(6,2);
t224 = t264 * rSges(6,1);
t457 = rSges(6,2) * pkin(4);
t115 = t224 - t457;
t505 = m(6) * t115 - Icges(6,5);
t114 = t224 + t457;
t504 = m(6) * t114 - Icges(6,5);
t116 = t225 + t456;
t503 = -t116 * m(6) + Icges(6,6);
t246 = t272 + qJ(3);
t219 = qJ(4) + t246;
t183 = qJ(5) + t219;
t128 = sin(t183);
t184 = -qJ(5) + t219;
t129 = sin(t184);
t502 = t129 + t128;
t243 = qJ(4) + qJ(5);
t201 = sin(t243);
t232 = qJD(4) + qJD(5);
t375 = t232 * t467;
t312 = rSges(6,1) * t201 * t375;
t208 = cos(t243);
t313 = rSges(6,2) * t208 * t375;
t501 = -t312 / 0.2e1 - t313 / 0.2e1;
t190 = qJD(3) + t232;
t376 = t190 * t468;
t214 = qJ(3) + t243;
t483 = rSges(6,2) * cos(t214);
t316 = t376 * t483;
t157 = sin(t214);
t317 = rSges(6,1) * t157 * t376;
t500 = -t316 / 0.2e1 - t317 / 0.2e1;
t256 = sin(qJ(3));
t261 = cos(qJ(3));
t499 = (-t256 * t109 - t261 * t466) * pkin(2) * qJD(3);
t447 = m(6) * t264;
t459 = rSges(5,2) * m(5);
t119 = -t447 + t459;
t210 = cos(t245);
t234 = qJD(3) + qJD(4);
t480 = t171 * sin(t245);
t498 = (-t210 * t119 - t480) * pkin(2) * t234;
t254 = sin(qJ(5));
t259 = cos(qJ(5));
t497 = -t254 * rSges(6,1) - t259 * rSges(6,2);
t194 = qJD(2) + t234;
t235 = qJD(2) + qJD(5);
t192 = qJD(3) + t235;
t496 = t192 * (rSges(6,1) * t161 + rSges(6,2) * t168);
t205 = sin(t247);
t212 = cos(t247);
t495 = t235 * (rSges(6,1) * t205 + rSges(6,2) * t212);
t104 = (-t273 + t274) * m(6) - t506;
t269 = 2 * qJ(5);
t180 = t269 + t220;
t125 = sin(t180);
t187 = -2 * qJ(5) + t220;
t132 = sin(t187);
t135 = cos(t180);
t185 = qJ(2) + t214;
t140 = cos(t185);
t215 = qJ(3) + t244;
t186 = qJ(2) + t215;
t141 = cos(t186);
t163 = sin(t220);
t170 = cos(t220);
t458 = rSges(6,2) * m(6);
t175 = rSges(6,1) * t458 - Icges(6,4);
t249 = pkin(1) * qJD(1);
t195 = t249 * t464;
t196 = Icges(6,5) * t378;
t226 = t273 + t274;
t433 = pkin(4) * qJD(1);
t230 = pkin(1) * t433;
t240 = cos(t269);
t310 = t104 * t378 / 0.2e1;
t349 = t175 * t378;
t366 = rSges(6,2) * t249;
t391 = Icges(6,3) * qJD(5);
t407 = t175 * t194;
t422 = (-rSges(5,1) * t459 + pkin(4) * t447 + Icges(5,4)) * cos(t149);
t178 = t270 + t354;
t423 = t503 * sin(t178);
t424 = t505 * cos(t178);
t425 = t504 * cos(t179);
t432 = qJD(1) * m(6);
t445 = pkin(1) * t268;
t469 = -0.2e1 * qJD(5);
t110 = 2 * t185;
t96 = sin(t110);
t111 = 2 * t186;
t97 = sin(t111);
t98 = cos(t110);
t99 = cos(t111);
t494 = ((qJD(5) * t226 + t230) * m(6) + t391 + t195) * qJD(1) * t163 + (-(t366 / 0.2e1 + qJD(5) * t114) * t432 + t196) * t141 + (-(-t366 / 0.2e1 + t115 * qJD(5)) * t432 + t196) * t140 + t119 * t170 * t445 + t240 * t407 * t469 + t135 * t349 + (t125 + t132) * t310 + (-t514 - t423) * t450 + (-t422 - t424 / 0.2e1 - t425 / 0.2e1 + (t98 / 0.4e1 - t99 / 0.4e1) * t175 + (t96 + t97) * t104 / 0.8e1) * t268;
t493 = 2 * Icges(6,4);
t253 = qJD(3) / 0.2e1;
t439 = qJD(2) / 0.2e1;
t488 = pkin(4) * (t439 + t253 + t252);
t486 = qJD(1) / 0.2e1;
t238 = sin(t269);
t381 = qJD(5) * t264;
t418 = t104 * t194;
t442 = pkin(4) * t194;
t474 = 0.2e1 * rSges(6,2);
t478 = -t238 * t418 - t259 * ((rSges(6,1) * t381 + t442 * t474) * m(6) - Icges(6,5) * qJD(5));
t263 = cos(qJ(1));
t382 = qJD(5) * t163;
t258 = sin(qJ(1));
t71 = t194 * t258;
t25 = t263 * t382 - t71;
t72 = t194 * t263;
t26 = t258 * t382 + t72;
t410 = t163 * t263;
t412 = t163 * t258;
t397 = t259 * t263;
t399 = t258 * t254;
t53 = -t170 * t399 + t397;
t398 = t258 * t259;
t401 = t254 * t263;
t54 = t170 * t398 + t401;
t55 = -t170 * t401 - t398;
t56 = t170 * t397 - t399;
t95 = qJD(5) * t170 + qJD(1);
t278 = (-Icges(6,6) * t410 + t506 * t55 + t56 * t510) * t25 - (Icges(6,6) * t412 + t54 * t493 - t506 * t53) * t26 - (Icges(6,6) * t170 + (-t506 * t254 + t259 * t510) * t163) * t95;
t475 = 0.16e2 * m(6);
t473 = -0.2e1 * t119;
t472 = -0.2e1 * t171;
t470 = -0.2e1 * t249;
t465 = m(4) * rSges(4,3);
t463 = m(5) * rSges(5,3);
t461 = rSges(6,1) * m(6);
t460 = rSges(3,2) * m(3);
t148 = t268 + 0.2e1 * t511;
t455 = t148 / 0.2e1;
t152 = qJD(2) + t190;
t454 = -t152 / 0.2e1;
t189 = t268 + 0.2e1 * t509;
t452 = t189 / 0.2e1;
t443 = pkin(3) * t119;
t250 = qJD(5) / 0.2e1;
t431 = t470 * rSges(6,1);
t427 = ((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + (m(6) + t266) * pkin(2) ^ 2) * sin(t272);
t420 = t497 * t258;
t419 = t497 * t263;
t112 = t250 + t188;
t416 = t112 * t232;
t414 = t152 * t264;
t191 = qJD(3) + t233;
t153 = qJD(2) + t191;
t413 = t153 * t264;
t404 = (rSges(3,1) * t460 - Icges(3,4)) * cos(t272);
t403 = t194 * t264;
t402 = t213 * t237;
t262 = cos(qJ(2));
t400 = t256 * t262;
t166 = cos(t216);
t396 = t268 * t166;
t204 = sin(t246);
t395 = t268 * t204;
t394 = pkin(4) * t250 + t249;
t392 = Icges(6,6) * qJD(5);
t388 = qJD(1) * t175;
t387 = qJD(1) * t194;
t385 = qJD(2) * t262;
t260 = cos(qJ(4));
t384 = qJD(4) * t260;
t383 = qJD(5) * t479;
t380 = t116 * qJD(5);
t228 = qJD(2) + t253;
t377 = m(6) * t456;
t236 = qJD(2) - qJD(5);
t370 = pkin(3) * t432;
t369 = pkin(2) * t432;
t368 = rSges(6,1) * t442;
t367 = t468 / 0.4e1;
t229 = rSges(6,1) * t249;
t169 = cos(t219);
t365 = t268 * pkin(2) * t169;
t364 = t463 / 0.2e1;
t363 = -t461 / 0.2e1;
t359 = pkin(2) * t452;
t357 = rSges(6,2) * t381;
t351 = t407 / 0.2e1;
t155 = t252 + t228;
t348 = pkin(2) * t371;
t346 = pkin(2) * cos(t215) * t458;
t345 = -qJD(1) * t104 / 0.4e1;
t107 = t250 + t155;
t341 = t107 * t369;
t108 = t251 + t155;
t340 = t108 * t369;
t339 = t112 * t370;
t338 = t113 * t370;
t337 = t189 * t367;
t336 = t236 * t468 / 0.2e1;
t193 = qJD(3) + t236;
t335 = t193 * t467 / 0.2e1;
t162 = sin(t219);
t332 = pkin(2) * t162 * t408;
t328 = t384 * t443;
t138 = cos(t183);
t321 = t138 * t348;
t139 = cos(t184);
t320 = t139 * t348;
t158 = sin(t215);
t319 = pkin(2) * t191 * t158 * t461;
t202 = sin(t244);
t318 = rSges(6,1) * t202 * t374;
t315 = t191 * t346;
t299 = pkin(4) * t170 + pkin(6) * t163;
t101 = rSges(6,1) * t259 - rSges(6,2) * t254;
t293 = rSges(6,3) * t163 + t101 * t170;
t207 = cos(t242);
t291 = rSges(6,2) * t207 * t336;
t167 = cos(t217);
t290 = rSges(6,2) * t167 * t335;
t200 = sin(t242);
t289 = t236 * pkin(2) * t200 * t363;
t285 = t193 * pkin(3) * t160 * t363;
t173 = pkin(3) * t261 + pkin(2);
t283 = -pkin(3) * qJD(3) * (t257 * t261 + t400) - qJD(2) * (pkin(3) * t400 + t173 * t257);
t282 = (Icges(6,5) * t55 - Icges(6,6) * t56) * t25 + (Icges(6,5) * t53 - Icges(6,6) * t54) * t26 + (Icges(6,5) * t254 + Icges(6,6) * t259) * t163 * t95;
t281 = t163 * t282;
t279 = (Icges(6,5) * t410 + t55 * t493 + t506 * t56) * t25 + (Icges(6,5) * t412 + t53 * t493 + t506 * t54) * t26 + (Icges(6,5) * t170 + (t254 * t493 - t259 * t506) * t163) * t95;
t130 = sin(t185);
t131 = sin(t186);
t142 = cos(t187);
t277 = -t320 / 0.4e1 + t321 / 0.4e1 + (-t142 * t388 - 0.2e1 * t254 * ((t368 - t357 / 0.2e1) * m(6) + t392 / 0.2e1) + t478) * qJD(5) + t337 * t483 + t359 * t480 + t332 / 0.2e1 - t189 * t346 / 0.4e1 + (((t229 - 0.2e1 * t383) * t475 + 0.32e2 * t392) * t131 + ((t229 + 0.2e1 * t380) * t475 - 0.32e2 * t392) * t130) * qJD(1) / 0.32e2 + (t365 / 0.2e1 + t210 * t359) * t119 + t494 + ((t158 + t157) * t337 + t502 * t268 * t367) * rSges(6,1);
t211 = cos(t246);
t70 = -pkin(4) * t163 + pkin(6) * t170;
t60 = -pkin(3) * t256 * t257 + t173 * t262 + pkin(1);
t45 = t299 * t263;
t44 = t299 * t258;
t40 = t497 * t163;
t24 = rSges(6,3) * t170 - t101 * t163;
t20 = -t101 * t258 + t170 * t419;
t19 = t101 * t263 + t170 * t420;
t6 = t263 * t293 + t420;
t5 = t258 * t293 - t419;
t1 = [(t315 + t314) * t486 + (((rSges(6,1) * t414 + (t394 + t488) * t474) * m(6) - Icges(6,5) * t152) * t141 + ((rSges(6,2) * t414 + t431) * m(6) + (-Icges(6,6) - t377) * t152) * t131 + t99 * t388) * t153 / 0.2e1 - 0.2e1 * t511 * (rSges(4,2) * t465 / 0.2e1 - Icges(4,6) / 0.2e1) * t206 + (t96 * t345 + ((rSges(6,2) * t413 + t431) * m(6) + (-Icges(6,6) + t377) * t153) * t130 / 0.2e1) * t152 + Icges(3,6) * t513 + t238 * t310 + ((t463 + t465) * pkin(2) + rSges(3,3) * t462 - Icges(3,5)) * qJD(2) * t385 + (((t394 - t488) * t474 * m(6) - Icges(6,5) * t153) * t140 + t98 * t388) * t454 + (-t136 * t339 + t137 * t338 - t138 * t341 + t139 * t340) * rSges(6,2) + (t470 * t466 + t237 * (rSges(4,3) * t265 + pkin(3) * t463 - Icges(4,5))) * t402 + (t290 + t285) * t192 + t497 * m(6) * qJD(5) * (pkin(1) * qJD(5) + t433) + t153 * t97 * t345 + t240 * t349 + (t413 * m(6) * t140 * t454 - t126 * t339 - t127 * t338 - t128 * t341 - t129 * t340) * rSges(6,1) + (((t162 * t472 + t169 * t473) * t155 - 0.2e1 * (t109 * t204 + t211 * t466) * t228) * pkin(2) - t319 / 0.2e1 - t318 / 0.2e1 - t327 - t328 + t498 + t499 + t500 + t501 + (t423 + t424) * (t250 + t194) + (t159 * t472 + t166 * t473) * pkin(3) * t188 + (t425 + t514) * (t251 + t194) + (-t426 - 0.2e1 * t405 - 0.2e1 * t508) * t237 + (-t427 - 0.2e1 * t404 - 0.2e1 * t507) * qJD(2)) * qJD(1) - 0.2e1 * ((t350 * qJD(5) + t230) * m(6) + t195 + t391 / 0.2e1 + t194 * (rSges(5,2) * t364 - Icges(5,6) / 0.2e1)) * t194 * t163 + (t135 * t351 + t125 * t418 / 0.4e1) * (0.2e1 * qJD(5) + t194) + (t142 * t351 - t132 * t418 / 0.4e1) * (t469 + t194) + t336 * t495 + t335 * t496 + (t289 + t291) * t235 + 0.2e1 * (-t119 * t249 + t194 * (rSges(5,1) * t364 - Icges(5,5) / 0.2e1)) * t194 * t170 + (-rSges(3,3) * t513 + t470 * t385) * t460 + (0.2e1 * t422 - t515) * t387; (-t316 - t317) * t107 + (((t229 / 0.2e1 + t380) * m(6) - t392) * t130 + (-(-t229 / 0.2e1 + t383) * m(6) + t392) * t131) * qJD(1) + t427 * t450 - t320 / 0.2e1 + t321 / 0.2e1 + (-t312 - t313) * t112 + (((t357 - 0.2e1 * t368) * m(6) - t392) * t254 + t478) * qJD(5) + t396 * t443 + t494 - t142 * t349 + (t315 - t319) * t108 + (t404 + t507) * t268 + 0.2e1 * t499 * t228 + ((-t207 + t212) * t329 + t502 * t331 + (t200 + t205) * t330 + t395 * t109 + t211 * t372) * pkin(2) - t167 * t302 - t113 * t318 + t119 * t365 + 0.2e1 * t498 * t155 + t516 + t262 * t445 * t460 + t332 + t328 * t471; ((t384 * t471 + t396) * t119 + ((-t167 * t378 - t208 * t416) * rSges(6,2) + (-t113 * t202 * t233 - t201 * t416) * rSges(6,1)) * m(6)) * pkin(3) + ((t256 * t452 + t395 / 0.2e1) * t109 + (t211 * t450 + t261 * t452) * t466) * pkin(2) + t277 + t516; t277 + (t409 * t455 + (t260 * t455 + t396 / 0.2e1) * t119 + (((t136 / 0.4e1 - t137 / 0.4e1) * t268 + (t208 / 0.4e1 - t209 / 0.4e1) * t148) * rSges(6,2) + ((t127 / 0.4e1 + t126 / 0.4e1) * t268 + (t201 / 0.4e1 + t202 / 0.4e1) * t148) * rSges(6,1)) * m(6)) * pkin(3) + t13 + t74 / 0.2e1; -t95 * (t282 * t170 + (t279 * t254 - t278 * t259) * t163) / 0.2e1 + t319 * t439 - t25 * (t263 * t281 + t278 * t56 + t279 * t55) / 0.2e1 - t26 * (t258 * t281 + t278 * t54 + t279 * t53) / 0.2e1 - m(6) * ((t24 * t26 - t5 * t95 + t70 * t72 + t283 * t263 + (-t258 * t60 - t44) * qJD(1)) * (-t19 * t95 - t26 * t40) + (-t24 * t25 + t6 * t95 + t70 * t71 + t283 * t258 + (t263 * t60 + t45) * qJD(1)) * (t20 * t95 + t25 * t40) + (-pkin(2) * t385 - pkin(3) * t402 + t25 * t5 - t26 * t6 - t44 * t71 - t45 * t72) * (t19 * t25 - t20 * t26)) + qJD(1) * t289 + qJD(1) * t290 + qJD(1) * t291 + qJD(1) * t285 - (m(6) * t226 + Icges(6,3)) * t163 * t387 + (((-rSges(6,1) * t403 - t366) * m(6) + Icges(6,5) * t194) * t259 + ((rSges(6,2) * t403 - t229) * m(6) - Icges(6,6) * t194) * t254) * qJD(5) + (-t315 / 0.2e1 + t500) * qJD(2) + (-t314 / 0.2e1 + t318 / 0.2e1 + t501) * t237 - (pkin(2) * t495 + pkin(3) * t496) * t432 / 0.2e1 + ((t503 * t130 + t505 * t140) * t152 + (t512 * t131 + t504 * t141) * t153) * t486;];
tauc = t1(:);
