% Calculate kinetic energy for
% palh3m2TE
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2TE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2TE_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:48
% EndTime: 2020-05-07 01:41:57
% DurationCPUTime: 6.91s
% Computational Cost: add. (2640->494), mult. (4307->724), div. (0->0), fcn. (4945->22), ass. (0->270)
t463 = Icges(3,3) + Icges(7,3);
t327 = sin(pkin(16));
t328 = cos(pkin(16));
t337 = sin(pkin(15));
t343 = cos(pkin(15));
t254 = t327 * t343 + t328 * t337;
t255 = -t327 * t337 + t328 * t343;
t325 = pkin(17) + pkin(18);
t313 = sin(t325);
t314 = cos(t325);
t200 = t254 * t314 + t255 * t313;
t329 = sin(pkin(18));
t331 = cos(pkin(18));
t256 = t329 * t343 + t331 * t337;
t368 = t329 * t337 - t331 * t343;
t455 = t313 * t254 - t255 * t314;
t462 = -Icges(5,5) * t455 - Icges(8,5) * t368 + Icges(5,6) * t200 + Icges(8,6) * t256;
t338 = sin(pkin(14));
t344 = cos(pkin(14));
t264 = t337 * t344 - t343 * t338;
t265 = t338 * t337 + t343 * t344;
t335 = sin(qJ(2));
t341 = cos(qJ(2));
t204 = t264 * t341 + t265 * t335;
t205 = -t264 * t335 + t265 * t341;
t461 = Icges(3,5) * t341 + Icges(7,5) * t205 - Icges(3,6) * t335 - Icges(7,6) * t204;
t460 = Icges(2,2) + Icges(5,3) + Icges(8,3);
t340 = cos(qJ(3));
t407 = pkin(4) * t340 - pkin(1);
t334 = sin(qJ(3));
t423 = t334 * t335;
t420 = pkin(4) * t423 - t407 * t341;
t244 = pkin(12) + t420;
t458 = Icges(7,4) * t204;
t248 = t256 * pkin(3);
t249 = t368 * pkin(3);
t330 = sin(pkin(17));
t332 = cos(pkin(17));
t457 = -t248 * t332 + t249 * t330 - pkin(11);
t445 = pkin(1) * t341;
t309 = pkin(12) + t445;
t456 = Icges(5,4) * t455;
t296 = -rSges(8,1) * t337 + rSges(8,2) * t343;
t297 = rSges(8,1) * t343 + rSges(8,2) * t337;
t421 = -t296 * t329 - t297 * t331;
t336 = sin(qJ(1));
t342 = cos(qJ(1));
t434 = Icges(7,4) * t205;
t381 = -Icges(7,2) * t204 + t434;
t157 = -Icges(7,6) * t342 + t336 * t381;
t158 = Icges(7,6) * t336 + t342 * t381;
t387 = Icges(7,1) * t205 - t458;
t159 = -Icges(7,5) * t342 + t336 * t387;
t160 = Icges(7,5) * t336 + t342 * t387;
t162 = Icges(7,2) * t205 + t458;
t163 = Icges(7,1) * t204 + t434;
t438 = Icges(3,4) * t341;
t384 = -Icges(3,2) * t335 + t438;
t236 = -Icges(3,6) * t342 + t336 * t384;
t237 = Icges(3,6) * t336 + t342 * t384;
t439 = Icges(3,4) * t335;
t390 = Icges(3,1) * t341 - t439;
t238 = -Icges(3,5) * t342 + t336 * t390;
t239 = Icges(3,5) * t336 + t342 * t390;
t280 = Icges(3,2) * t341 + t439;
t283 = Icges(3,1) * t335 + t438;
t419 = qJD(2) * t342;
t303 = V_base(5) - t419;
t316 = qJD(2) * t336;
t304 = V_base(4) + t316;
t315 = V_base(6) + qJD(1);
t454 = (-t162 * t204 + t163 * t205 - t280 * t335 + t283 * t341) * t315 + (-t158 * t204 + t160 * t205 - t237 * t335 + t239 * t341) * t304 + (-t157 * t204 + t159 * t205 - t236 * t335 + t238 * t341) * t303;
t453 = (Icges(3,5) * t335 + Icges(7,5) * t204 + Icges(3,6) * t341 + Icges(7,6) * t205) * t315 + (t463 * t336 + t461 * t342) * t304 + (t461 * t336 - t463 * t342) * t303;
t324 = qJD(2) + qJD(3);
t452 = -t324 * t342 + V_base(5);
t394 = rSges(3,1) * t341 - rSges(3,2) * t335;
t367 = pkin(12) + t394;
t451 = t336 * rSges(3,3) + t342 * t367;
t319 = Icges(2,4) * t342;
t440 = Icges(2,4) * t336;
t450 = (-Icges(5,5) * t200 - Icges(8,5) * t256 - Icges(5,6) * t455 - Icges(8,6) * t368) * t315 + (t462 * t336 + t460 * t342 + t440) * V_base(5) + (-t460 * t336 + t462 * t342 + t319) * V_base(4);
t382 = -Icges(5,2) * t200 + t456;
t148 = -Icges(5,6) * t342 + t336 * t382;
t149 = Icges(5,6) * t336 + t342 * t382;
t435 = Icges(5,4) * t200;
t388 = Icges(5,1) * t455 - t435;
t150 = -Icges(5,5) * t342 + t336 * t388;
t151 = Icges(5,5) * t336 + t342 * t388;
t153 = Icges(5,2) * t455 + t435;
t154 = Icges(5,1) * t200 + t456;
t432 = Icges(8,4) * t368;
t380 = -Icges(8,2) * t256 + t432;
t186 = -Icges(8,6) * t342 + t336 * t380;
t187 = Icges(8,6) * t336 + t342 * t380;
t433 = Icges(8,4) * t256;
t386 = Icges(8,1) * t368 - t433;
t188 = -Icges(8,5) * t342 + t336 * t386;
t189 = Icges(8,5) * t336 + t342 * t386;
t208 = Icges(8,2) * t368 + t433;
t209 = Icges(8,1) * t256 + t432;
t449 = (-t153 * t200 + t154 * t455 - t208 * t256 + t209 * t368) * t315 + (Icges(2,1) * t336 - t148 * t200 + t150 * t455 - t186 * t256 + t188 * t368 + t319) * V_base(5) + (Icges(2,1) * t342 - t149 * t200 + t151 * t455 - t187 * t256 + t189 * t368 - t440) * V_base(4);
t444 = pkin(1) * qJD(2);
t443 = pkin(4) * qJD(3);
t441 = t342 * rSges(3,3);
t262 = -t340 * t341 + t423;
t437 = Icges(4,4) * t262;
t422 = t334 * t341;
t263 = t335 * t340 + t422;
t436 = Icges(4,4) * t263;
t326 = qJ(3) + qJ(2);
t317 = sin(t326);
t431 = Icges(9,4) * t317;
t318 = cos(t326);
t430 = Icges(9,4) * t318;
t286 = pkin(8) * t327 + pkin(10) * t328;
t287 = pkin(8) * t328 - pkin(10) * t327;
t216 = t286 * t343 + t287 * t337;
t217 = -t286 * t337 + t287 * t343;
t429 = (-t216 * t313 + t217 * t314) * t315;
t428 = t200 * t336;
t427 = t200 * t342;
t425 = t407 * t335;
t418 = qJD(4) * t200;
t362 = t256 * t330;
t417 = ((-t332 * t368 - t362) * pkin(3) - t309) * qJD(1);
t416 = t309 * qJD(1);
t415 = t341 * t444 + V_base(3);
t321 = V_base(5) * pkin(11);
t414 = t321 + V_base(1);
t413 = pkin(4) * t422;
t410 = t335 * t444;
t409 = t263 * t443;
t322 = V_base(4) * pkin(12);
t397 = pkin(4) * t334 * V_base(4);
t405 = t341 * V_base(4);
t408 = t335 * t397 - t407 * t405 + t322;
t406 = t335 * V_base(5);
t404 = t342 * V_base(5);
t403 = t342 * V_base(6);
t399 = t249 * t332 + pkin(12);
t402 = -pkin(3) * t362 - t399 - t445;
t400 = -pkin(11) + t425;
t396 = t342 * t409 - t406 * t407 + t414;
t395 = t315 * t244;
t333 = sin(qJ(4));
t339 = cos(qJ(4));
t393 = rSges(6,1) * t339 - rSges(6,2) * t333;
t392 = -rSges(9,1) * t318 + rSges(9,2) * t317;
t391 = -V_base(4) * pkin(11) + V_base(2);
t389 = Icges(4,1) * t262 + t436;
t385 = -Icges(9,1) * t318 + t431;
t383 = Icges(4,2) * t263 + t437;
t379 = Icges(9,2) * t317 - t430;
t377 = Icges(4,5) * t262 + Icges(4,6) * t263;
t373 = -Icges(9,5) * t318 + Icges(9,6) * t317;
t178 = rSges(6,3) * t255 + t254 * t393;
t179 = -rSges(6,3) * t254 + t255 * t393;
t372 = t178 * t313 - t179 * t314;
t273 = rSges(5,1) * t327 - rSges(5,2) * t328;
t274 = rSges(5,1) * t328 + rSges(5,2) * t327;
t214 = t273 * t343 + t274 * t337;
t215 = -t273 * t337 + t274 * t343;
t371 = -t214 * t313 + t215 * t314;
t290 = -t343 * rSges(7,1) + rSges(7,2) * t337;
t293 = rSges(7,1) * t337 + rSges(7,2) * t343;
t369 = t338 * t290 + t293 * t344;
t366 = t420 * qJD(2) + t262 * t443 + V_base(3);
t365 = t304 * t335;
t364 = t333 * t455;
t363 = t455 * t339;
t361 = pkin(1) * t406 - t342 * t410 + V_base(1);
t360 = t336 * V_base(4) - t404;
t261 = qJD(3) * t336 + t304;
t358 = (-Icges(4,3) * t342 + t336 * t377) * t452 + (Icges(4,3) * t336 + t342 * t377) * t261 + (-Icges(4,5) * t263 + Icges(4,6) * t262) * t315;
t276 = t324 * t336 + V_base(4);
t357 = (-Icges(9,3) * t342 + t336 * t373) * t452 + (Icges(9,3) * t336 + t342 * t373) * t276 + (-Icges(9,5) * t317 - Icges(9,6) * t318) * t315;
t355 = t321 + t361;
t354 = t244 * V_base(5);
t245 = -t413 - t425;
t353 = -t245 * t316 + t336 * t409 + t341 * t397 + V_base(2) + (qJD(1) * t342 + t403) * t244;
t194 = -Icges(4,6) * t342 + t336 * t383;
t195 = Icges(4,6) * t336 + t342 * t383;
t196 = -Icges(4,5) * t342 + t336 * t389;
t197 = Icges(4,5) * t336 + t342 * t389;
t211 = Icges(4,2) * t262 - t436;
t212 = -Icges(4,1) * t263 + t437;
t349 = (t195 * t263 + t197 * t262) * t261 + (t194 * t263 + t196 * t262) * t452 + (t211 * t263 + t212 * t262) * t315;
t227 = -Icges(9,6) * t342 + t336 * t379;
t228 = Icges(9,6) * t336 + t342 * t379;
t229 = -Icges(9,5) * t342 + t336 * t385;
t230 = Icges(9,5) * t336 + t342 * t385;
t251 = -Icges(9,2) * t318 - t431;
t252 = -Icges(9,1) * t317 - t430;
t348 = (t228 * t317 - t230 * t318) * t276 + (t227 * t317 - t229 * t318) * t452 + (t251 * t317 - t252 * t318) * t315;
t307 = pkin(1) * t405;
t301 = pkin(8) * t343 - pkin(10) * t337;
t300 = pkin(8) * t337 + pkin(10) * t343;
t295 = rSges(5,1) * t343 + rSges(5,2) * t337;
t294 = -rSges(5,1) * t337 + rSges(5,2) * t343;
t292 = rSges(2,1) * t342 - rSges(2,2) * t336;
t291 = rSges(6,1) * t333 + rSges(6,2) * t339;
t289 = rSges(2,1) * t336 + rSges(2,2) * t342;
t288 = rSges(3,1) * t335 + rSges(3,2) * t341;
t279 = Icges(2,5) * t342 - Icges(2,6) * t336;
t278 = Icges(2,5) * t336 + Icges(2,6) * t342;
t272 = t342 * t416;
t271 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t270 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t269 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t267 = t309 * V_base(6);
t253 = -rSges(9,1) * t317 - rSges(9,2) * t318;
t232 = rSges(9,3) * t336 + t342 * t392;
t231 = -rSges(9,3) * t342 + t336 * t392;
t223 = -t290 * t344 + t338 * t293;
t222 = (-rSges(4,1) * t341 + rSges(4,2) * t335) * t334 - t340 * (rSges(4,1) * t335 + rSges(4,2) * t341);
t221 = (-rSges(4,1) * t340 + rSges(4,2) * t334) * t341 + t335 * (rSges(4,1) * t334 + rSges(4,2) * t340);
t220 = V_base(5) * rSges(2,3) - t289 * t315 + t414;
t219 = t292 * t315 + V_base(2) + (-pkin(11) - rSges(2,3)) * V_base(4);
t218 = -t296 * t331 + t297 * t329;
t213 = t289 * V_base(4) - t292 * V_base(5) + V_base(3);
t203 = rSges(4,3) * t336 + t221 * t342;
t202 = -rSges(4,3) * t342 + t221 * t336;
t191 = t248 * t330 + t399;
t182 = -qJD(4) * t455 + t315;
t181 = t342 * t418 + V_base(4);
t180 = t336 * t418 + V_base(5);
t177 = t333 * t336 + t342 * t363;
t176 = t336 * t339 - t342 * t364;
t175 = -t333 * t342 + t336 * t363;
t174 = -t336 * t364 - t339 * t342;
t173 = t303 * t288 + (-t336 * t367 + t441) * t315 + t414;
t172 = -t304 * t288 + t451 * t315 + t391;
t171 = t223 * t341 - t335 * t369;
t170 = t223 * t335 + t341 * t369;
t168 = t216 * t314 + t217 * t313;
t167 = t214 * t314 + t215 * t313;
t166 = t394 * qJD(2) + V_base(3) + t322 * t336 + (t336 * t394 - t441) * V_base(4) - t451 * V_base(5);
t165 = rSges(7,3) * t336 + t171 * t342;
t164 = -rSges(7,3) * t342 + t171 * t336;
t145 = (-rSges(8,3) * V_base(4) + (-t309 - t421) * V_base(5)) * t342 + (t307 + (pkin(12) + t421) * V_base(4) - V_base(5) * rSges(8,3)) * t336 + t415;
t144 = t178 * t314 + t179 * t313;
t143 = t272 + V_base(2) - t336 * t410 + t309 * t403 + t315 * (t336 * rSges(8,3) + t342 * t421) + (-pkin(1) * t335 - pkin(11) - t218) * V_base(4);
t142 = V_base(5) * t218 + t355 + (t342 * rSges(8,3) + ((rSges(8,1) * t331 + rSges(8,2) * t329) * t343 - t337 * (rSges(8,1) * t329 - rSges(8,2) * t331) - t309) * t336) * t315;
t141 = -Icges(6,5) * t455 + (Icges(6,1) * t339 - Icges(6,4) * t333) * t200;
t140 = -Icges(6,6) * t455 + (Icges(6,4) * t339 - Icges(6,2) * t333) * t200;
t139 = -Icges(6,3) * t455 + (Icges(6,5) * t339 - Icges(6,6) * t333) * t200;
t138 = -pkin(1) * t365 + t203 * t315 - t222 * t261 + t267 * t342 + t272 + t391;
t137 = -t202 * t315 + t222 * t452 + (-t267 - t416) * t336 + t355;
t136 = -t291 * t342 + t336 * t372;
t135 = t291 * t336 + t342 * t372;
t134 = t202 * t261 - t203 * t452 + t309 * t360 + t415;
t133 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t427;
t132 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t428;
t131 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t427;
t130 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t428;
t129 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t427;
t128 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t428;
t127 = V_base(5) * pkin(13) + t170 * t303 + (pkin(6) * t336 - t164) * t315 + t414;
t126 = -t170 * t304 + V_base(2) + (-pkin(11) - pkin(13)) * V_base(4) + (-pkin(6) * t342 + t165) * t315;
t125 = -pkin(6) * t360 + t164 * t304 - t165 * t303 + V_base(3);
t124 = (t191 * V_base(4) + t307) * t336 + t276 * t231 - t452 * t232 + t402 * t404 + t415;
t123 = (-rSges(5,3) * V_base(4) - t354) * t342 + (-V_base(5) * rSges(5,3) + t408) * t336 - ((t294 * t327 + t295 * t328) * t314 + (t294 * t328 - t295 * t327) * t313) * t360 + t366;
t122 = V_base(2) + t457 * V_base(4) - t276 * t253 + t315 * t232 + (t191 * V_base(6) - t417) * t342 + (t341 * t403 - t365) * pkin(1);
t121 = -t457 * V_base(5) + t452 * t253 - t315 * t231 + (t402 * V_base(6) + t417) * t336 + t361;
t120 = -t315 * (-rSges(5,3) * t336 + t342 * t371) + (-t167 + t400) * V_base(4) + t353;
t119 = (t167 - t413) * V_base(5) + (rSges(5,3) * t315 - qJD(2) * t245) * t342 + (t315 * t371 - t395) * t336 + t396;
t118 = -t354 * t342 + t408 * t336 + t181 * t136 - t180 * t135 - ((-t300 * t327 + t301 * t328) * t314 - (t300 * t328 + t301 * t327) * t313) * t360 + t366;
t117 = -t342 * t429 + t135 * t182 - t144 * t181 + (-t168 + t400) * V_base(4) + t353;
t116 = -t245 * t419 - t136 * t182 + t144 * t180 + (t168 - t413) * V_base(5) + (-t395 + t429) * t336 + t396;
t1 = t276 * (t357 * t336 + t348 * t342) / 0.2e1 + t261 * (t358 * t336 + t349 * t342) / 0.2e1 + t182 * (-(t128 * t180 + t129 * t181 + t139 * t182) * t455 + ((-t131 * t333 + t133 * t339) * t181 + (-t130 * t333 + t132 * t339) * t180 + (-t140 * t333 + t141 * t339) * t182) * t200) / 0.2e1 + t181 * ((t129 * t427 + t176 * t131 + t177 * t133) * t181 + (t128 * t427 + t130 * t176 + t132 * t177) * t180 + (t139 * t427 + t140 * t176 + t141 * t177) * t182) / 0.2e1 + t180 * ((t129 * t428 + t131 * t174 + t133 * t175) * t181 + (t128 * t428 + t174 * t130 + t175 * t132) * t180 + (t139 * t428 + t140 * t174 + t141 * t175) * t182) / 0.2e1 + m(2) * (t213 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(1) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(3) * (t166 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(8) * (t142 ^ 2 + t143 ^ 2 + t145 ^ 2) / 0.2e1 + m(4) * (t134 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(7) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t123 ^ 2) / 0.2e1 + m(9) * (t121 ^ 2 + t122 ^ 2 + t124 ^ 2) / 0.2e1 + m(6) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + (t454 * t336 - t453 * t342) * t303 / 0.2e1 + (t453 * t336 + t454 * t342) * t304 / 0.2e1 + (Icges(1,1) * V_base(4) + t279 * t315 - t450 * t336 + t449 * t342) * V_base(4) / 0.2e1 + (Icges(1,2) * V_base(5) + t278 * t315 + t449 * t336 + t450 * t342) * V_base(5) / 0.2e1 + ((-t228 * t318 - t230 * t317) * t276 + (t195 * t262 - t197 * t263) * t261 + (t158 * t205 + t160 * t204 + t237 * t341 + t239 * t335) * t304 + (t157 * t205 + t159 * t204 + t236 * t341 + t238 * t335) * t303 + (t148 * t455 + t150 * t200 + t186 * t368 + t188 * t256 + t278) * V_base(5) + (t149 * t455 + t151 * t200 + t187 * t368 + t189 * t256 + t279) * V_base(4) + (t153 * t455 + t154 * t200 + t162 * t205 + t163 * t204 + t208 * t368 + t209 * t256 + t211 * t262 - t212 * t263 - t251 * t318 - t252 * t317 + t280 * t341 + t283 * t335 + Icges(2,3)) * t315 + (t194 * t262 - t196 * t263 - t227 * t318 - t229 * t317) * t452) * t315 / 0.2e1 + V_base(5) * V_base(4) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((-t357 - t358) * t342 + (t348 + t349) * t336) * t452 / 0.2e1;
T = t1;
