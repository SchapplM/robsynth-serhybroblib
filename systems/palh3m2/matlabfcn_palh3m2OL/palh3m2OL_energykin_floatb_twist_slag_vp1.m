% Calculate kinetic energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(6,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: qJD has to be [10x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:24
% EndTime: 2020-05-07 04:33:30
% DurationCPUTime: 4.94s
% Computational Cost: add. (2121->428), mult. (2110->650), div. (0->0), fcn. (1832->18), ass. (0->229)
t280 = sin(qJ(1));
t284 = cos(qJ(1));
t367 = -t280 * V_base(4) + V_base(5) * t284;
t279 = sin(qJ(2));
t365 = pkin(1) * t279;
t275 = qJ(2) + qJ(7);
t270 = pkin(15) - t275;
t364 = pkin(3) * sin(t270);
t276 = qJ(2) + qJ(3);
t267 = sin(t276);
t363 = pkin(4) * t267;
t362 = pkin(12) * t280;
t283 = cos(qJ(2));
t361 = t283 * pkin(1);
t360 = Icges(2,4) * t280;
t359 = Icges(3,4) * t279;
t358 = Icges(3,4) * t283;
t357 = Icges(4,4) * t267;
t269 = cos(t276);
t356 = Icges(4,4) * t269;
t272 = qJ(4) + t276;
t258 = sin(t272);
t355 = Icges(5,4) * t258;
t259 = cos(t272);
t354 = Icges(5,4) * t259;
t277 = sin(qJ(6));
t353 = Icges(7,4) * t277;
t281 = cos(qJ(6));
t352 = Icges(7,4) * t281;
t266 = sin(t275);
t351 = Icges(8,4) * t266;
t268 = cos(t275);
t350 = Icges(8,4) * t268;
t260 = -qJ(8) + t270;
t253 = sin(t260);
t349 = Icges(9,4) * t253;
t254 = cos(t260);
t348 = Icges(9,4) * t254;
t347 = t258 * t280;
t346 = t258 * t284;
t278 = sin(qJ(5));
t345 = t278 * t284;
t344 = t280 * t278;
t282 = cos(qJ(5));
t343 = t280 * t282;
t342 = t282 * t284;
t341 = pkin(3) * cos(t270);
t340 = pkin(4) * t269;
t339 = qJD(5) * t258;
t338 = -qJD(2) - qJD(3);
t337 = -qJD(2) - qJD(7);
t336 = V_base(5) * pkin(11) + V_base(1);
t252 = qJD(2) * t280 + V_base(4);
t261 = V_base(6) + qJD(1);
t211 = t361 * t280;
t331 = -t211 - t362;
t250 = -qJD(2) * t284 + V_base(5);
t330 = t250 * t365 + t336;
t224 = qJD(7) * t280 + t252;
t225 = qJD(3) * t280 + t252;
t155 = t340 * t280;
t329 = t155 + t331;
t328 = -pkin(8) * t259 - pkin(10) * t258;
t327 = rSges(3,1) * t283 - rSges(3,2) * t279;
t326 = -rSges(4,1) * t269 + rSges(4,2) * t267;
t325 = -rSges(5,1) * t259 + rSges(5,2) * t258;
t324 = rSges(7,1) * t281 - rSges(7,2) * t277;
t323 = rSges(8,1) * t268 - rSges(8,2) * t266;
t322 = -rSges(9,1) * t254 - rSges(9,2) * t253;
t205 = qJD(4) * t280 + t225;
t321 = Icges(3,1) * t283 - t359;
t320 = -Icges(4,1) * t269 + t357;
t319 = -Icges(5,1) * t259 + t355;
t318 = Icges(7,1) * t281 - t353;
t317 = Icges(8,1) * t268 - t351;
t316 = -Icges(9,1) * t254 - t349;
t315 = -Icges(3,2) * t279 + t358;
t314 = Icges(4,2) * t267 - t356;
t313 = Icges(5,2) * t258 - t354;
t312 = -Icges(7,2) * t277 + t352;
t311 = -Icges(8,2) * t266 + t350;
t310 = -Icges(9,2) * t253 - t348;
t309 = Icges(3,5) * t283 - Icges(3,6) * t279;
t308 = -Icges(4,5) * t269 + Icges(4,6) * t267;
t307 = -Icges(5,5) * t259 + Icges(5,6) * t258;
t306 = Icges(7,5) * t281 - Icges(7,6) * t277;
t305 = Icges(8,5) * t268 - Icges(8,6) * t266;
t304 = -Icges(9,5) * t254 - Icges(9,6) * t253;
t303 = t261 * t284 * pkin(12) - V_base(4) * pkin(11) + V_base(2);
t223 = t338 * t284 + V_base(5);
t302 = -t223 * t363 + t330;
t301 = -t367 * pkin(12) + V_base(3);
t203 = V_base(5) + (-qJD(4) + t338) * t284;
t202 = V_base(5) + (-qJD(8) + t337) * t284;
t204 = qJD(8) * t280 + t224;
t300 = (-Icges(9,3) * t284 + t280 * t304) * t202 + (Icges(9,3) * t280 + t284 * t304) * t204 + (Icges(9,5) * t253 - Icges(9,6) * t254) * t261;
t299 = (-Icges(5,3) * t284 + t280 * t307) * t203 + (Icges(5,3) * t280 + t284 * t307) * t205 + (-Icges(5,5) * t258 - Icges(5,6) * t259) * t261;
t222 = t337 * t284 + V_base(5);
t298 = (-Icges(8,3) * t284 + t280 * t305) * t222 + (Icges(8,3) * t280 + t284 * t305) * t224 + (Icges(8,5) * t266 + Icges(8,6) * t268) * t261;
t297 = (-Icges(4,3) * t284 + t280 * t308) * t223 + (Icges(4,3) * t280 + t284 * t308) * t225 + (-Icges(4,5) * t267 - Icges(4,6) * t269) * t261;
t249 = -qJD(6) * t284 + V_base(5);
t251 = qJD(6) * t280 + V_base(4);
t296 = (-Icges(7,3) * t284 + t280 * t306) * t249 + (Icges(7,3) * t280 + t284 * t306) * t251 + (Icges(7,5) * t277 + Icges(7,6) * t281) * t261;
t295 = (-Icges(3,3) * t284 + t280 * t309) * t250 + (Icges(3,3) * t280 + t284 * t309) * t252 + (Icges(3,5) * t279 + Icges(3,6) * t283) * t261;
t212 = t361 * t284;
t294 = t261 * t212 - t252 * t365 + t303;
t293 = t252 * t211 - t250 * t212 + t301;
t156 = t340 * t284;
t292 = -t261 * t156 + t225 * t363 + t294;
t291 = -t225 * t155 + t223 * t156 + t293;
t133 = -Icges(9,6) * t284 + t280 * t310;
t134 = Icges(9,6) * t280 + t284 * t310;
t135 = -Icges(9,5) * t284 + t280 * t316;
t136 = Icges(9,5) * t280 + t284 * t316;
t194 = -Icges(9,2) * t254 + t349;
t195 = Icges(9,1) * t253 - t348;
t290 = (-t134 * t253 - t136 * t254) * t204 + (-t133 * t253 - t135 * t254) * t202 + (-t194 * t253 - t195 * t254) * t261;
t149 = -Icges(5,6) * t284 + t280 * t313;
t150 = Icges(5,6) * t280 + t284 * t313;
t151 = -Icges(5,5) * t284 + t280 * t319;
t152 = Icges(5,5) * t280 + t284 * t319;
t207 = -Icges(5,2) * t259 - t355;
t208 = -Icges(5,1) * t258 - t354;
t289 = (t150 * t258 - t152 * t259) * t205 + (t149 * t258 - t151 * t259) * t203 + (t207 * t258 - t208 * t259) * t261;
t161 = -Icges(8,6) * t284 + t280 * t311;
t162 = Icges(8,6) * t280 + t284 * t311;
t165 = -Icges(8,5) * t284 + t280 * t317;
t166 = Icges(8,5) * t280 + t284 * t317;
t215 = Icges(8,2) * t268 + t351;
t217 = Icges(8,1) * t266 + t350;
t288 = (-t162 * t266 + t166 * t268) * t224 + (-t161 * t266 + t165 * t268) * t222 + (-t215 * t266 + t217 * t268) * t261;
t163 = -Icges(4,6) * t284 + t280 * t314;
t164 = Icges(4,6) * t280 + t284 * t314;
t167 = -Icges(4,5) * t284 + t280 * t320;
t168 = Icges(4,5) * t280 + t284 * t320;
t216 = -Icges(4,2) * t269 - t357;
t218 = -Icges(4,1) * t267 - t356;
t287 = (t164 * t267 - t168 * t269) * t225 + (t163 * t267 - t167 * t269) * t223 + (t216 * t267 - t218 * t269) * t261;
t181 = -Icges(7,6) * t284 + t280 * t312;
t182 = Icges(7,6) * t280 + t284 * t312;
t185 = -Icges(7,5) * t284 + t280 * t318;
t186 = Icges(7,5) * t280 + t284 * t318;
t236 = Icges(7,2) * t281 + t353;
t240 = Icges(7,1) * t277 + t352;
t286 = (-t182 * t277 + t186 * t281) * t251 + (-t181 * t277 + t185 * t281) * t249 + (-t236 * t277 + t240 * t281) * t261;
t183 = -Icges(3,6) * t284 + t280 * t315;
t184 = Icges(3,6) * t280 + t284 * t315;
t187 = -Icges(3,5) * t284 + t280 * t321;
t188 = Icges(3,5) * t280 + t284 * t321;
t237 = Icges(3,2) * t283 + t359;
t241 = Icges(3,1) * t279 + t358;
t285 = (-t184 * t279 + t188 * t283) * t252 + (-t183 * t279 + t187 * t283) * t250 + (-t237 * t279 + t241 * t283) * t261;
t271 = Icges(2,4) * t284;
t247 = rSges(2,1) * t284 - t280 * rSges(2,2);
t246 = t280 * rSges(2,1) + rSges(2,2) * t284;
t245 = rSges(3,1) * t279 + rSges(3,2) * t283;
t244 = rSges(7,1) * t277 + rSges(7,2) * t281;
t243 = Icges(2,1) * t284 - t360;
t242 = Icges(2,1) * t280 + t271;
t239 = -Icges(2,2) * t280 + t271;
t238 = Icges(2,2) * t284 + t360;
t231 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t230 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t229 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t226 = qJD(5) * t259 + t261;
t221 = -rSges(4,1) * t267 - rSges(4,2) * t269;
t220 = rSges(8,1) * t266 + rSges(8,2) * t268;
t210 = -pkin(8) * t258 + pkin(10) * t259;
t209 = -rSges(5,1) * t258 - rSges(5,2) * t259;
t201 = -t259 * t342 + t344;
t200 = t259 * t345 + t343;
t199 = -t259 * t343 - t345;
t198 = t259 * t344 - t342;
t197 = rSges(9,1) * t253 - rSges(9,2) * t254;
t192 = t280 * rSges(3,3) + t284 * t327;
t191 = t280 * rSges(7,3) + t284 * t324;
t190 = -rSges(3,3) * t284 + t280 * t327;
t189 = -rSges(7,3) * t284 + t280 * t324;
t175 = t328 * t284;
t174 = t328 * t280;
t173 = t280 * rSges(4,3) + t284 * t326;
t172 = t280 * rSges(8,3) + t284 * t323;
t171 = -rSges(4,3) * t284 + t280 * t326;
t170 = -rSges(8,3) * t284 + t280 * t323;
t154 = t280 * rSges(5,3) + t284 * t325;
t153 = -rSges(5,3) * t284 + t280 * t325;
t146 = V_base(5) * rSges(2,3) - t246 * t261 + t336;
t145 = t247 * t261 + V_base(2) + (-pkin(11) - rSges(2,3)) * V_base(4);
t144 = t341 * t284;
t143 = t341 * t280;
t142 = -t284 * t339 + t205;
t141 = -t280 * t339 + t203;
t140 = t246 * V_base(4) - t247 * V_base(5) + V_base(3);
t139 = t280 * rSges(9,3) + t284 * t322;
t138 = -rSges(9,3) * t284 + t280 * t322;
t137 = rSges(6,3) * t259 + (-rSges(6,1) * t282 + rSges(6,2) * t278) * t258;
t130 = Icges(6,5) * t259 + (-Icges(6,1) * t282 + Icges(6,4) * t278) * t258;
t129 = Icges(6,6) * t259 + (-Icges(6,4) * t282 + Icges(6,2) * t278) * t258;
t128 = Icges(6,3) * t259 + (-Icges(6,5) * t282 + Icges(6,6) * t278) * t258;
t125 = t201 * rSges(6,1) + t200 * rSges(6,2) - rSges(6,3) * t346;
t124 = rSges(6,1) * t199 + rSges(6,2) * t198 - rSges(6,3) * t347;
t123 = Icges(6,1) * t201 + Icges(6,4) * t200 - Icges(6,5) * t346;
t122 = Icges(6,1) * t199 + Icges(6,4) * t198 - Icges(6,5) * t347;
t121 = Icges(6,4) * t201 + Icges(6,2) * t200 - Icges(6,6) * t346;
t120 = Icges(6,4) * t199 + Icges(6,2) * t198 - Icges(6,6) * t347;
t119 = Icges(6,5) * t201 + Icges(6,6) * t200 - Icges(6,3) * t346;
t118 = Icges(6,5) * t199 + Icges(6,6) * t198 - Icges(6,3) * t347;
t117 = t245 * t250 + (-t190 - t362) * t261 + t336;
t116 = t192 * t261 - t245 * t252 + t303;
t115 = V_base(5) * pkin(13) + t244 * t249 + (pkin(6) * t280 - t189) * t261 + t336;
t114 = -t251 * t244 + V_base(2) + (-pkin(11) - pkin(13)) * V_base(4) + (-pkin(6) * t284 + t191) * t261;
t113 = t367 * pkin(6) + t251 * t189 - t249 * t191 + V_base(3);
t112 = t252 * t190 - t250 * t192 + t301;
t111 = t221 * t223 + (-t171 + t331) * t261 + t330;
t110 = t220 * t222 + (-t170 + t331) * t261 + t330;
t109 = t173 * t261 - t221 * t225 + t294;
t108 = t172 * t261 - t220 * t224 + t294;
t107 = t225 * t171 - t223 * t173 + t293;
t106 = t224 * t170 - t222 * t172 + t293;
t105 = t203 * t209 + (-t153 + t329) * t261 + t302;
t104 = t154 * t261 - t205 * t209 + t292;
t103 = -t222 * t364 + t197 * t202 + (-t138 - t143 + t331) * t261 + t330;
t102 = t224 * t364 - t197 * t204 + (t139 + t144) * t261 + t294;
t101 = t205 * t153 - t203 * t154 + t291;
t100 = t204 * t138 - t202 * t139 + t224 * t143 - t222 * t144 + t293;
t99 = -t124 * t226 + t137 * t141 + t203 * t210 + (-t174 + t329) * t261 + t302;
t98 = t125 * t226 - t137 * t142 + t175 * t261 - t205 * t210 + t292;
t97 = t142 * t124 - t141 * t125 + t205 * t174 - t203 * t175 + t291;
t1 = ((-t280 * t238 + t242 * t284 + Icges(1,4)) * V_base(5) + (-t280 * t239 + t243 * t284 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t184 * t283 + t188 * t279) * t252 + (t183 * t283 + t187 * t279) * t250 + (t182 * t281 + t186 * t277) * t251 + (t181 * t281 + t185 * t277) * t249 + (-t164 * t269 - t168 * t267) * t225 + (-t163 * t269 - t167 * t267) * t223 + (t162 * t268 + t166 * t266) * t224 + (t161 * t268 + t165 * t266) * t222 + (-t150 * t259 - t152 * t258) * t205 + (-t149 * t259 - t151 * t258) * t203 + (-t134 * t254 + t136 * t253) * t204 + (-t133 * t254 + t135 * t253) * t202 + (-t194 * t254 + t195 * t253 - t207 * t259 - t208 * t258 + t215 * t268 - t216 * t269 + t217 * t266 - t218 * t267 + t236 * t281 + t237 * t283 + t240 * t277 + t241 * t279 + Icges(2,3)) * t261) * t261 / 0.2e1 + V_base(5) * t261 * (Icges(2,5) * t280 + Icges(2,6) * t284) + t261 * V_base(4) * (Icges(2,5) * t284 - Icges(2,6) * t280) + ((t238 * t284 + t280 * t242 + Icges(1,2)) * V_base(5) + (t239 * t284 + t280 * t243 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t251 * (t280 * t296 + t284 * t286) / 0.2e1 + t249 * (t280 * t286 - t284 * t296) / 0.2e1 + t225 * (t280 * t297 + t284 * t287) / 0.2e1 + t223 * (t280 * t287 - t284 * t297) / 0.2e1 + t224 * (t280 * t298 + t284 * t288) / 0.2e1 + t222 * (t280 * t288 - t284 * t298) / 0.2e1 + t205 * (t280 * t299 + t284 * t289) / 0.2e1 + t203 * (t280 * t289 - t284 * t299) / 0.2e1 + t204 * (t280 * t300 + t284 * t290) / 0.2e1 + t202 * (t280 * t290 - t284 * t300) / 0.2e1 + t252 * (t280 * t295 + t284 * t285) / 0.2e1 + t250 * (t280 * t285 - t284 * t295) / 0.2e1 + t141 * ((-t119 * t347 + t121 * t198 + t123 * t199) * t142 + (-t118 * t347 + t198 * t120 + t199 * t122) * t141 + (-t128 * t347 + t129 * t198 + t130 * t199) * t226) / 0.2e1 + t142 * ((-t119 * t346 + t200 * t121 + t201 * t123) * t142 + (-t118 * t346 + t200 * t120 + t201 * t122) * t141 + (-t128 * t346 + t200 * t129 + t201 * t130) * t226) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + t226 * ((t118 * t141 + t119 * t142 + t128 * t226) * t259 + ((t121 * t278 - t123 * t282) * t142 + (t120 * t278 - t122 * t282) * t141 + (t129 * t278 - t130 * t282) * t226) * t258) / 0.2e1 + m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + m(2) * (t140 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(8) * (t106 ^ 2 + t108 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t109 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(3) * (t112 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(9) * (t100 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1;
T = t1;
