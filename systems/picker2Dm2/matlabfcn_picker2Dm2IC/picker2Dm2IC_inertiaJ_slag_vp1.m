% Calculate joint inertia matrix for
% picker2Dm2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm2IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_inertiaJ_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2IC_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2IC_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:53
% EndTime: 2020-05-11 09:20:55
% DurationCPUTime: 1.23s
% Computational Cost: add. (9033->294), mult. (4177->344), div. (313->10), fcn. (3150->59), ass. (0->189)
t319 = qJ(1) + qJ(2);
t307 = sin(t319);
t309 = cos(t319);
t321 = sin(qJ(7));
t325 = cos(qJ(7));
t374 = (t307 * t321 + t309 * t325) * pkin(3);
t326 = cos(qJ(1));
t373 = pkin(1) * t326;
t328 = 0.1e1 / pkin(4);
t372 = pkin(1) * t328;
t357 = qJ(3) - qJ(6);
t345 = -qJ(9) - t357;
t203 = 0.1e1 / ((-cos(qJ(4) - t357) + cos(qJ(4) + t357)) * pkin(6) + (cos(qJ(4) + t345) - cos(qJ(4) - t345)) * pkin(2));
t371 = pkin(2) * t203;
t370 = pkin(2) * t309;
t369 = pkin(3) * t309;
t358 = qJ(2) + qJ(4);
t313 = qJ(1) + t358;
t297 = cos(t313);
t368 = pkin(4) * t297;
t327 = 0.1e1 / pkin(5);
t318 = qJ(1) + qJ(8);
t316 = pkin(8) + qJ(5);
t338 = t316 + t357;
t330 = t338 - t318;
t339 = t316 - t357;
t331 = t339 - t318;
t367 = 0.1e1 / (pkin(6) * (cos(t331) - cos(t330)) + (-cos(qJ(9) - t331) + cos(qJ(9) + t330)) * pkin(2)) * t327;
t329 = 0.1e1 / pkin(3);
t366 = t203 * t329;
t294 = sin(t313);
t335 = t294 * t309 - t297 * t307;
t217 = 0.1e1 / t335;
t334 = t294 * t321 + t297 * t325;
t365 = t217 * t334;
t305 = sin(t358);
t322 = sin(qJ(4));
t249 = -pkin(1) * t305 - pkin(3) * t322;
t317 = 0.1e1 / t322;
t364 = t249 * t317;
t363 = t249 * t329;
t323 = sin(qJ(2));
t362 = (pkin(3) * t323 + pkin(4) * t305) * t329;
t361 = t317 * t328;
t359 = Icges(4,3) + Icges(10,3);
t312 = qJ(6) + t319;
t293 = sin(t312);
t296 = cos(t312);
t237 = -rSges(7,1) * t293 - rSges(7,2) * t296;
t324 = sin(qJ(1));
t315 = t324 * pkin(1);
t222 = t237 + t315;
t240 = t296 * rSges(7,1) - rSges(7,2) * t293;
t223 = t240 - t373;
t186 = m(7) * (t222 * t237 + t223 * t240);
t185 = Icges(7,3) + t186;
t346 = qJ(9) + t319;
t302 = qJ(3) + t346;
t281 = sin(t302);
t284 = cos(t302);
t229 = t281 * rSges(10,1) + t284 * rSges(10,2);
t314 = qJ(3) + t319;
t295 = sin(t314);
t215 = -pkin(6) * t295 + t229;
t292 = pkin(2) * t307;
t207 = t215 + t292;
t201 = t207 + t315;
t230 = -rSges(10,1) * t284 + t281 * rSges(10,2);
t298 = cos(t314);
t216 = pkin(6) * t298 + t230;
t208 = t216 - t370;
t202 = t208 - t373;
t172 = Icges(10,3) + m(10) * (t201 * t229 + t202 * t230);
t301 = -qJ(9) + t313;
t265 = t301 - t357;
t300 = qJ(4) + t346;
t266 = t300 + t357;
t356 = sin(t266) - sin(t265);
t355 = cos(t266) - cos(t265);
t238 = t294 * rSges(5,1) + t297 * rSges(5,2);
t354 = -cos(0.2e1 * t318) + cos(0.2e1 * t316);
t353 = sin(t301) - sin(t300);
t352 = cos(t301) - cos(t300);
t246 = t307 * rSges(3,1) + t309 * rSges(3,2);
t299 = qJ(10) + t313;
t275 = sin(t299);
t277 = cos(t299);
t224 = -rSges(11,1) * t275 - rSges(11,2) * t277;
t274 = pkin(4) * t294;
t209 = t224 + t274;
t291 = pkin(3) * t307;
t204 = t209 + t291;
t199 = t204 + t315;
t226 = t277 * rSges(11,1) - rSges(11,2) * t275;
t210 = t226 - t368;
t205 = t210 - t369;
t200 = t205 - t373;
t170 = Icges(11,3) + m(11) * (t199 * t224 + t200 * t226);
t351 = Icges(5,3) + Icges(11,3);
t189 = -t355 * t321 + t356 * t325;
t350 = t189 * t371;
t344 = -qJ(1) + t316;
t340 = qJ(4) + t344;
t336 = -qJ(8) + t340;
t341 = -qJ(4) + t344;
t337 = -qJ(8) + t341;
t349 = (pkin(3) * (cos(t341) - cos(t340)) + (cos(qJ(2) - t337) - cos(qJ(2) + t336)) * pkin(5)) / (cos(t337) - cos(t336)) * t327;
t206 = -t354 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1);
t236 = 0.1e1 / t354;
t348 = t206 * t327 * t236;
t347 = t317 * t363;
t218 = t291 + t238;
t248 = -rSges(3,1) * t309 + t307 * rSges(3,2);
t241 = -rSges(5,1) * t297 + t294 * rSges(5,2);
t242 = t298 * rSges(4,1) - rSges(4,2) * t295;
t306 = sin(t318);
t308 = cos(t318);
t247 = t308 * rSges(9,1) - rSges(9,2) * t306;
t343 = pkin(1) * t349;
t342 = t329 * t349;
t239 = -rSges(4,1) * t295 - rSges(4,2) * t298;
t220 = t239 + t292;
t212 = t220 + t315;
t221 = t242 - t370;
t214 = t221 - t373;
t161 = m(4) * (t212 * t239 + t214 * t242) + m(10) * (t201 * t215 + t202 * t216) + t359;
t211 = t315 + t218;
t219 = t241 - t369;
t213 = t219 - t373;
t160 = m(5) * (t211 * t238 + t213 * t241) + m(11) * (t199 * t209 + t200 * t210) + t351;
t245 = -rSges(9,1) * t306 - rSges(9,2) * t308;
t332 = Icges(3,3) + Icges(7,3) + t351 + t359;
t233 = t315 + t246;
t235 = t248 - t373;
t149 = m(5) * (t211 * t218 + t213 * t219) + m(4) * (t212 * t220 + t214 * t221) + m(3) * (t233 * t246 + t235 * t248) + m(11) * (t199 * t204 + t200 * t205) + m(10) * (t201 * t207 + t202 * t208) + t186 + t332;
t320 = sin(qJ(8));
t311 = -qJ(9) + t316;
t310 = qJ(9) + t316;
t304 = cos(t316);
t303 = sin(t316);
t286 = qJ(9) + t338;
t285 = -qJ(9) + t339;
t276 = -sin(qJ(8) - t344);
t260 = -rSges(2,1) * t326 + rSges(2,2) * t324;
t259 = rSges(8,1) * t325 - rSges(8,2) * t321;
t258 = rSges(2,1) * t324 + rSges(2,2) * t326;
t257 = rSges(8,1) * t321 + rSges(8,2) * t325;
t252 = -pkin(5) * t306 + t315;
t251 = -pkin(5) * t308 + t373;
t244 = rSges(6,1) * t304 - rSges(6,2) * t303;
t243 = -rSges(6,1) * t303 - rSges(6,2) * t304;
t234 = t247 - t373;
t232 = t245 + t315;
t228 = t368 + t369 + t373;
t227 = t315 + t291 + t274;
t198 = pkin(4) * t334 + t374;
t197 = m(7) * (t237 ^ 2 + t240 ^ 2);
t196 = Icges(7,3) + t197;
t194 = t352 * t321 - t353 * t325;
t193 = pkin(4) * t335 - t374;
t192 = Icges(10,3) + m(10) * (t229 ^ 2 + t230 ^ 2);
t191 = Icges(11,3) + m(11) * (t224 ^ 2 + t226 ^ 2);
t182 = Icges(10,3) + m(10) * (t215 * t229 + t216 * t230);
t181 = m(10) * (t215 ^ 2 + t216 ^ 2) + m(4) * (t239 ^ 2 + t242 ^ 2) + t359;
t180 = Icges(11,3) + m(11) * (t209 * t224 + t210 * t226);
t179 = m(11) * (t209 ^ 2 + t210 ^ 2) + m(5) * (t238 ^ 2 + t241 ^ 2) + t351;
t176 = m(10) * (t207 * t229 + t208 * t230) + Icges(10,3);
t175 = m(11) * (t204 * t224 + t205 * t226) + Icges(11,3);
t173 = t194 * t196 * t371;
t165 = t176 * t365;
t164 = -t196 * t365 + t173;
t163 = m(10) * (t207 * t215 + t208 * t216) + m(4) * (t220 * t239 + t221 * t242) + t359;
t162 = m(11) * (t204 * t209 + t205 * t210) + m(5) * (t218 * t238 + t219 * t241) + t351;
t159 = ((t353 * t227 + t352 * t228) * t366 + (-(sin(t311) - sin(t310)) * t252 - (cos(t311) - cos(t310)) * t251) * t367) * pkin(2);
t158 = ((-t356 * t227 - t355 * t228) * t366 + ((sin(t286) - sin(t285)) * t252 + (cos(t286) - cos(t285)) * t251) * t367) * pkin(2);
t157 = t197 + m(3) * (t246 ^ 2 + t248 ^ 2) + m(11) * (t204 ^ 2 + t205 ^ 2) + m(10) * (t207 ^ 2 + t208 ^ 2) + m(5) * (t218 ^ 2 + t219 ^ 2) + m(4) * (t220 ^ 2 + t221 ^ 2) + t332;
t156 = t159 * t196;
t155 = (t175 * t363 + (t180 * t362 - t191 * t323) * t372) * t317 + t170;
t154 = t182 * t350 - t192 * t365 - t165;
t153 = (-t175 * t334 + (t180 * t198 + t191 * t193) * t328) * t217;
t152 = t196 * t347 + t156 + t185;
t151 = t181 * t350 - (t163 + t182) * t365;
t150 = (-t162 * t334 + (t179 * t198 + t180 * t193) * t328) * t217;
t148 = (t162 * t363 + (t179 * t362 - t180 * t323) * t372) * t317 + t160;
t147 = t158 * t182 + (t176 * t364 - t192 * t343) * t329 + t172;
t146 = t158 * t181 + (t163 * t364 - t182 * t343) * t329 + t161;
t145 = t163 * t350 - t165 + t173 + (-t157 * t334 + (t162 * t198 + t175 * t193) * t328) * t217;
t144 = t156 + (-t176 * t342 + (t162 * t362 - t175 * t323) * t361) * pkin(1) + t158 * t163 + t149 + t157 * t347;
t1 = [t332 + t206 ^ 2 * t327 ^ 2 * t236 ^ 2 * m(9) * (t245 ^ 2 + t247 ^ 2) + ((-t147 - t172) * t342 + ((-t155 - t170) * t323 + (t148 + t160) * t362) * t361) * pkin(1) + m(2) * (t258 ^ 2 + t260 ^ 2) + m(9) * (t232 ^ 2 + t234 ^ 2) + m(3) * (t233 ^ 2 + t235 ^ 2) + m(7) * (t222 ^ 2 + t223 ^ 2) + m(11) * (t199 ^ 2 + t200 ^ 2) + m(10) * (t201 ^ 2 + t202 ^ 2) + m(5) * (t211 ^ 2 + t213 ^ 2) + m(4) * (t212 ^ 2 + t214 ^ 2) + (t185 + t152) * t159 + (0.2e1 * m(9) * (t232 * t245 + t234 * t247) + (0.2e1 + t348) * Icges(9,3)) * t348 + (t146 + t161) * t158 + (t144 + t149) * t347 + Icges(2,3) + Icges(9,3) + (Icges(6,3) + m(6) * (t243 ^ 2 + t244 ^ 2)) * t320 ^ 2 / t276 ^ 2, (t146 * t189 + t152 * t194) * t371 + ((t148 * t198 + t155 * t193) * t328 - (t144 + t147) * t334) * t217; t145 * t347 + t151 * t158 + t164 * t159 + (t161 * t189 + t185 * t194) * t371 + (-t154 * t342 + (t150 * t362 - t153 * t323) * t361) * pkin(1) + ((t160 * t198 + t170 * t193) * t328 - (t149 + t172) * t334) * t217, Icges(8,3) + m(8) * (t257 ^ 2 + t259 ^ 2) + (t151 * t189 + t164 * t194) * t371 + ((t150 * t198 + t153 * t193) * t328 - (t145 + t154) * t334) * t217;];
Mq = t1;
