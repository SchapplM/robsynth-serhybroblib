% Calculate kinetic energy for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:57:41
% EndTime: 2020-05-07 01:57:50
% DurationCPUTime: 7.12s
% Computational Cost: add. (2640->494), mult. (4307->721), div. (0->0), fcn. (4945->22), ass. (0->272)
t469 = Icges(3,3) + Icges(7,3);
t332 = sin(pkin(18));
t334 = cos(pkin(18));
t340 = sin(pkin(15));
t346 = cos(pkin(15));
t258 = t332 * t346 + t334 * t340;
t371 = t332 * t340 - t334 * t346;
t330 = sin(pkin(16));
t331 = cos(pkin(16));
t256 = t330 * t346 + t331 * t340;
t257 = -t330 * t340 + t331 * t346;
t328 = pkin(17) + pkin(18);
t316 = sin(t328);
t317 = cos(t328);
t401 = t256 * t317 + t257 * t316;
t461 = t316 * t256 - t257 * t317;
t468 = -Icges(5,5) * t461 - Icges(8,5) * t371 + Icges(5,6) * t401 + Icges(8,6) * t258;
t341 = sin(pkin(14));
t347 = cos(pkin(14));
t266 = t340 * t347 - t341 * t346;
t267 = t341 * t340 + t346 * t347;
t338 = sin(qJ(2));
t344 = cos(qJ(2));
t205 = t266 * t344 + t267 * t338;
t206 = -t266 * t338 + t267 * t344;
t467 = Icges(3,5) * t344 + Icges(7,5) * t206 - Icges(3,6) * t338 - Icges(7,6) * t205;
t466 = Icges(2,2) + Icges(5,3) + Icges(8,3);
t343 = cos(qJ(3));
t411 = pkin(4) * t343 - pkin(1);
t337 = sin(qJ(3));
t427 = t337 * t338;
t424 = pkin(4) * t427 - t411 * t344;
t246 = pkin(12) + t424;
t464 = Icges(7,4) * t205;
t250 = t258 * pkin(3);
t251 = t371 * pkin(3);
t333 = sin(pkin(17));
t335 = cos(pkin(17));
t463 = -t250 * t335 + t251 * t333 - pkin(11);
t451 = pkin(1) * t344;
t312 = pkin(12) + t451;
t462 = Icges(5,4) * t461;
t298 = -rSges(8,1) * t340 + rSges(8,2) * t346;
t299 = rSges(8,1) * t346 + rSges(8,2) * t340;
t425 = -t298 * t332 - t299 * t334;
t339 = sin(qJ(1));
t345 = cos(qJ(1));
t440 = Icges(7,4) * t206;
t384 = -Icges(7,2) * t205 + t440;
t157 = -Icges(7,6) * t345 + t339 * t384;
t158 = Icges(7,6) * t339 + t345 * t384;
t390 = Icges(7,1) * t206 - t464;
t159 = -Icges(7,5) * t345 + t339 * t390;
t160 = Icges(7,5) * t339 + t345 * t390;
t162 = Icges(7,2) * t206 + t464;
t163 = Icges(7,1) * t205 + t440;
t444 = Icges(3,4) * t344;
t387 = -Icges(3,2) * t338 + t444;
t237 = -Icges(3,6) * t345 + t339 * t387;
t238 = Icges(3,6) * t339 + t345 * t387;
t445 = Icges(3,4) * t338;
t393 = Icges(3,1) * t344 - t445;
t239 = -Icges(3,5) * t345 + t339 * t393;
t240 = Icges(3,5) * t339 + t345 * t393;
t282 = Icges(3,2) * t344 + t445;
t285 = Icges(3,1) * t338 + t444;
t423 = qJD(2) * t345;
t305 = V_base(5) - t423;
t319 = qJD(2) * t339;
t306 = V_base(4) + t319;
t318 = V_base(6) + qJD(1);
t460 = (-t162 * t205 + t163 * t206 - t282 * t338 + t285 * t344) * t318 + (-t158 * t205 + t160 * t206 - t238 * t338 + t240 * t344) * t306 + (-t157 * t205 + t159 * t206 - t237 * t338 + t239 * t344) * t305;
t459 = (Icges(3,5) * t338 + Icges(7,5) * t205 + Icges(3,6) * t344 + Icges(7,6) * t206) * t318 + (t469 * t339 + t467 * t345) * t306 + (t467 * t339 - t469 * t345) * t305;
t327 = qJD(2) + qJD(3);
t458 = -t327 * t345 + V_base(5);
t397 = rSges(3,1) * t344 - rSges(3,2) * t338;
t370 = pkin(12) + t397;
t457 = t339 * rSges(3,3) + t345 * t370;
t322 = Icges(2,4) * t345;
t446 = Icges(2,4) * t339;
t456 = (-Icges(5,5) * t401 - Icges(8,5) * t258 - Icges(5,6) * t461 - Icges(8,6) * t371) * t318 + (t468 * t339 + t466 * t345 + t446) * V_base(5) + (-t466 * t339 + t468 * t345 + t322) * V_base(4);
t385 = -Icges(5,2) * t401 + t462;
t148 = -Icges(5,6) * t345 + t339 * t385;
t149 = Icges(5,6) * t339 + t345 * t385;
t441 = Icges(5,4) * t401;
t391 = Icges(5,1) * t461 - t441;
t150 = -Icges(5,5) * t345 + t339 * t391;
t151 = Icges(5,5) * t339 + t345 * t391;
t153 = Icges(5,2) * t461 + t441;
t154 = Icges(5,1) * t401 + t462;
t438 = Icges(8,4) * t371;
t383 = -Icges(8,2) * t258 + t438;
t186 = -Icges(8,6) * t345 + t339 * t383;
t187 = Icges(8,6) * t339 + t345 * t383;
t439 = Icges(8,4) * t258;
t389 = Icges(8,1) * t371 - t439;
t188 = -Icges(8,5) * t345 + t339 * t389;
t189 = Icges(8,5) * t339 + t345 * t389;
t209 = Icges(8,2) * t371 + t439;
t210 = Icges(8,1) * t258 + t438;
t455 = (-t153 * t401 + t154 * t461 - t209 * t258 + t210 * t371) * t318 + (Icges(2,1) * t339 - t148 * t401 + t150 * t461 - t186 * t258 + t188 * t371 + t322) * V_base(5) + (Icges(2,1) * t345 - t149 * t401 + t151 * t461 - t187 * t258 + t189 * t371 - t446) * V_base(4);
t450 = pkin(1) * qJD(2);
t449 = pkin(4) * qJD(3);
t447 = t345 * rSges(3,3);
t264 = -t343 * t344 + t427;
t443 = Icges(4,4) * t264;
t426 = t337 * t344;
t265 = t338 * t343 + t426;
t442 = Icges(4,4) * t265;
t329 = qJ(3) + qJ(2);
t320 = sin(t329);
t437 = Icges(9,4) * t320;
t321 = cos(t329);
t436 = Icges(9,4) * t321;
t288 = pkin(8) * t330 + pkin(10) * t331;
t289 = pkin(8) * t331 - pkin(10) * t330;
t217 = t288 * t346 + t289 * t340;
t218 = -t288 * t340 + t289 * t346;
t435 = (-t217 * t316 + t218 * t317) * t318;
t434 = t401 * t339;
t342 = cos(qJ(4));
t433 = t401 * t342;
t432 = t401 * t345;
t336 = sin(qJ(4));
t431 = t401 * t336;
t429 = t411 * t338;
t422 = qJD(4) * t401;
t365 = t258 * t333;
t421 = ((-t335 * t371 - t365) * pkin(3) - t312) * qJD(1);
t420 = t312 * qJD(1);
t419 = t344 * t450 + V_base(3);
t324 = V_base(5) * pkin(11);
t418 = t324 + V_base(1);
t417 = pkin(4) * t426;
t414 = t338 * t450;
t413 = t265 * t449;
t325 = V_base(4) * pkin(12);
t400 = pkin(4) * t337 * V_base(4);
t409 = t344 * V_base(4);
t412 = t338 * t400 - t411 * t409 + t325;
t410 = t338 * V_base(5);
t408 = t345 * V_base(5);
t407 = t345 * V_base(6);
t403 = t251 * t335 + pkin(12);
t406 = -pkin(3) * t365 - t403 - t451;
t404 = -pkin(11) + t429;
t399 = t345 * t413 - t410 * t411 + t418;
t398 = t318 * t246;
t396 = rSges(6,1) * t342 - rSges(6,2) * t336;
t395 = -rSges(9,1) * t321 + rSges(9,2) * t320;
t394 = -V_base(4) * pkin(11) + V_base(2);
t392 = Icges(4,1) * t264 + t442;
t388 = -Icges(9,1) * t321 + t437;
t386 = Icges(4,2) * t265 + t443;
t382 = Icges(9,2) * t320 - t436;
t380 = Icges(4,5) * t264 + Icges(4,6) * t265;
t376 = -Icges(9,5) * t321 + Icges(9,6) * t320;
t178 = rSges(6,3) * t257 + t256 * t396;
t179 = -rSges(6,3) * t256 + t257 * t396;
t375 = t178 * t316 - t179 * t317;
t275 = rSges(5,1) * t330 - rSges(5,2) * t331;
t276 = rSges(5,1) * t331 + rSges(5,2) * t330;
t215 = t275 * t346 + t276 * t340;
t216 = -t275 * t340 + t276 * t346;
t374 = -t215 * t316 + t216 * t317;
t292 = -t346 * rSges(7,1) + rSges(7,2) * t340;
t295 = rSges(7,1) * t340 + rSges(7,2) * t346;
t372 = t341 * t292 + t295 * t347;
t369 = t424 * qJD(2) + t264 * t449 + V_base(3);
t368 = t306 * t338;
t367 = t336 * t461;
t366 = t461 * t342;
t364 = pkin(1) * t410 - t345 * t414 + V_base(1);
t363 = t339 * V_base(4) - t408;
t263 = qJD(3) * t339 + t306;
t361 = (-Icges(4,3) * t345 + t339 * t380) * t458 + (Icges(4,3) * t339 + t345 * t380) * t263 + (-Icges(4,5) * t265 + Icges(4,6) * t264) * t318;
t278 = t327 * t339 + V_base(4);
t360 = (-Icges(9,3) * t345 + t339 * t376) * t458 + (Icges(9,3) * t339 + t345 * t376) * t278 + (-Icges(9,5) * t320 - Icges(9,6) * t321) * t318;
t358 = t324 + t364;
t357 = t246 * V_base(5);
t247 = -t417 - t429;
t356 = -t247 * t319 + t339 * t413 + t344 * t400 + V_base(2) + (qJD(1) * t345 + t407) * t246;
t194 = -Icges(4,6) * t345 + t339 * t386;
t195 = Icges(4,6) * t339 + t345 * t386;
t196 = -Icges(4,5) * t345 + t339 * t392;
t197 = Icges(4,5) * t339 + t345 * t392;
t212 = Icges(4,2) * t264 - t442;
t213 = -Icges(4,1) * t265 + t443;
t352 = (t195 * t265 + t197 * t264) * t263 + (t194 * t265 + t196 * t264) * t458 + (t212 * t265 + t213 * t264) * t318;
t228 = -Icges(9,6) * t345 + t339 * t382;
t229 = Icges(9,6) * t339 + t345 * t382;
t230 = -Icges(9,5) * t345 + t339 * t388;
t231 = Icges(9,5) * t339 + t345 * t388;
t253 = -Icges(9,2) * t321 - t437;
t254 = -Icges(9,1) * t320 - t436;
t351 = (t229 * t320 - t231 * t321) * t278 + (t228 * t320 - t230 * t321) * t458 + (t253 * t320 - t254 * t321) * t318;
t309 = pkin(1) * t409;
t303 = pkin(8) * t346 - pkin(10) * t340;
t302 = pkin(8) * t340 + pkin(10) * t346;
t297 = rSges(5,1) * t346 + rSges(5,2) * t340;
t296 = -rSges(5,1) * t340 + rSges(5,2) * t346;
t294 = rSges(2,1) * t345 - rSges(2,2) * t339;
t293 = rSges(6,1) * t336 + rSges(6,2) * t342;
t291 = rSges(2,1) * t339 + rSges(2,2) * t345;
t290 = rSges(3,1) * t338 + rSges(3,2) * t344;
t281 = Icges(2,5) * t345 - Icges(2,6) * t339;
t280 = Icges(2,5) * t339 + Icges(2,6) * t345;
t274 = t345 * t420;
t273 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t272 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t271 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t269 = t312 * V_base(6);
t255 = -rSges(9,1) * t320 - rSges(9,2) * t321;
t233 = rSges(9,3) * t339 + t345 * t395;
t232 = -rSges(9,3) * t345 + t339 * t395;
t224 = -t292 * t347 + t341 * t295;
t223 = (-rSges(4,1) * t344 + rSges(4,2) * t338) * t337 - t343 * (rSges(4,1) * t338 + rSges(4,2) * t344);
t222 = (-rSges(4,1) * t343 + rSges(4,2) * t337) * t344 + t338 * (rSges(4,1) * t337 + rSges(4,2) * t343);
t221 = V_base(5) * rSges(2,3) - t291 * t318 + t418;
t220 = t294 * t318 + V_base(2) + (-pkin(11) - rSges(2,3)) * V_base(4);
t219 = -t298 * t334 + t299 * t332;
t214 = t291 * V_base(4) - t294 * V_base(5) + V_base(3);
t204 = rSges(4,3) * t339 + t222 * t345;
t203 = -rSges(4,3) * t345 + t222 * t339;
t191 = t250 * t333 + t403;
t182 = -qJD(4) * t461 + t318;
t181 = t345 * t422 + V_base(4);
t180 = t339 * t422 + V_base(5);
t177 = t336 * t339 + t345 * t366;
t176 = t339 * t342 - t345 * t367;
t175 = -t336 * t345 + t339 * t366;
t174 = -t339 * t367 - t342 * t345;
t173 = t305 * t290 + (-t339 * t370 + t447) * t318 + t418;
t172 = -t306 * t290 + t318 * t457 + t394;
t171 = t224 * t344 - t338 * t372;
t170 = t338 * t224 + t344 * t372;
t168 = t217 * t317 + t218 * t316;
t167 = t215 * t317 + t216 * t316;
t166 = t397 * qJD(2) + V_base(3) + t325 * t339 + (t339 * t397 - t447) * V_base(4) - t457 * V_base(5);
t165 = rSges(7,3) * t339 + t171 * t345;
t164 = -rSges(7,3) * t345 + t171 * t339;
t145 = (-rSges(8,3) * V_base(4) + (-t312 - t425) * V_base(5)) * t345 + (t309 + (pkin(12) + t425) * V_base(4) - V_base(5) * rSges(8,3)) * t339 + t419;
t144 = t178 * t317 + t179 * t316;
t143 = t274 + V_base(2) - t339 * t414 + t312 * t407 + t318 * (t339 * rSges(8,3) + t345 * t425) + (-pkin(1) * t338 - pkin(11) - t219) * V_base(4);
t142 = V_base(5) * t219 + t358 + (t345 * rSges(8,3) + ((rSges(8,1) * t334 + rSges(8,2) * t332) * t346 - t340 * (rSges(8,1) * t332 - rSges(8,2) * t334) - t312) * t339) * t318;
t141 = Icges(6,1) * t433 - Icges(6,4) * t431 - Icges(6,5) * t461;
t140 = Icges(6,4) * t433 - Icges(6,2) * t431 - Icges(6,6) * t461;
t139 = Icges(6,5) * t433 - Icges(6,6) * t431 - Icges(6,3) * t461;
t138 = -pkin(1) * t368 + t204 * t318 - t223 * t263 + t269 * t345 + t274 + t394;
t137 = -t203 * t318 + t223 * t458 + (-t269 - t420) * t339 + t358;
t136 = t293 * t339 + t345 * t375;
t135 = -t293 * t345 + t339 * t375;
t134 = t203 * t263 - t204 * t458 + t312 * t363 + t419;
t133 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t432;
t132 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t434;
t131 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t432;
t130 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t434;
t129 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t432;
t128 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t434;
t127 = V_base(5) * pkin(13) + t170 * t305 + (pkin(6) * t339 - t164) * t318 + t418;
t126 = -t170 * t306 + V_base(2) + (-pkin(11) - pkin(13)) * V_base(4) + (-pkin(6) * t345 + t165) * t318;
t125 = -pkin(6) * t363 + t164 * t306 - t165 * t305 + V_base(3);
t124 = (t191 * V_base(4) + t309) * t339 + t278 * t232 - t458 * t233 + t406 * t408 + t419;
t123 = (-rSges(5,3) * V_base(4) - t357) * t345 + (-V_base(5) * rSges(5,3) + t412) * t339 - ((t296 * t330 + t297 * t331) * t317 + (t296 * t331 - t297 * t330) * t316) * t363 + t369;
t122 = V_base(2) + t463 * V_base(4) - t278 * t255 + t318 * t233 + (t191 * V_base(6) - t421) * t345 + (t344 * t407 - t368) * pkin(1);
t121 = -t463 * V_base(5) + t458 * t255 - t318 * t232 + (t406 * V_base(6) + t421) * t339 + t364;
t120 = -(-t339 * rSges(5,3) + t345 * t374) * t318 + (-t167 + t404) * V_base(4) + t356;
t119 = (t167 - t417) * V_base(5) + (rSges(5,3) * t318 - qJD(2) * t247) * t345 + (t318 * t374 - t398) * t339 + t399;
t118 = -t357 * t345 + t412 * t339 + t181 * t135 - t180 * t136 - ((-t302 * t330 + t303 * t331) * t317 - (t302 * t331 + t303 * t330) * t316) * t363 + t369;
t117 = -t345 * t435 + t136 * t182 - t144 * t181 + (-t168 + t404) * V_base(4) + t356;
t116 = -t247 * t423 - t135 * t182 + t144 * t180 + (t168 - t417) * V_base(5) + (-t398 + t435) * t339 + t399;
t1 = t278 * (t360 * t339 + t351 * t345) / 0.2e1 + t263 * (t361 * t339 + t352 * t345) / 0.2e1 + t181 * ((t129 * t432 + t176 * t131 + t177 * t133) * t181 + (t128 * t432 + t130 * t176 + t132 * t177) * t180 + (t139 * t432 + t140 * t176 + t141 * t177) * t182) / 0.2e1 + t182 * ((-t129 * t461 - t131 * t431 + t133 * t433) * t181 + (-t128 * t461 - t130 * t431 + t132 * t433) * t180 + (-t139 * t461 - t140 * t431 + t141 * t433) * t182) / 0.2e1 + t180 * ((t129 * t434 + t131 * t174 + t133 * t175) * t181 + (t128 * t434 + t174 * t130 + t175 * t132) * t180 + (t139 * t434 + t140 * t174 + t141 * t175) * t182) / 0.2e1 + m(1) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + m(3) * (t166 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(4) * (t134 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(8) * (t142 ^ 2 + t143 ^ 2 + t145 ^ 2) / 0.2e1 + m(2) * (t214 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t123 ^ 2) / 0.2e1 + m(9) * (t121 ^ 2 + t122 ^ 2 + t124 ^ 2) / 0.2e1 + m(7) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(6) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + (t460 * t339 - t459 * t345) * t305 / 0.2e1 + (t459 * t339 + t460 * t345) * t306 / 0.2e1 + (Icges(1,1) * V_base(4) + t281 * t318 - t339 * t456 + t345 * t455) * V_base(4) / 0.2e1 + (Icges(1,2) * V_base(5) + t280 * t318 + t339 * t455 + t345 * t456) * V_base(5) / 0.2e1 + ((-t229 * t321 - t231 * t320) * t278 + (t195 * t264 - t197 * t265) * t263 + (t158 * t206 + t160 * t205 + t238 * t344 + t240 * t338) * t306 + (t157 * t206 + t159 * t205 + t237 * t344 + t239 * t338) * t305 + (t148 * t461 + t150 * t401 + t186 * t371 + t188 * t258 + t280) * V_base(5) + (t149 * t461 + t151 * t401 + t187 * t371 + t189 * t258 + t281) * V_base(4) + (t153 * t461 + t154 * t401 + t206 * t162 + t205 * t163 + t209 * t371 + t258 * t210 + t264 * t212 - t265 * t213 - t321 * t253 - t320 * t254 + t344 * t282 + t338 * t285 + Icges(2,3)) * t318 + (t194 * t264 - t196 * t265 - t228 * t321 - t230 * t320) * t458) * t318 / 0.2e1 + V_base(5) * V_base(4) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((-t360 - t361) * t345 + (t351 + t352) * t339) * t458 / 0.2e1;
T = t1;
