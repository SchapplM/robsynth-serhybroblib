% Calculate joint inertia matrix for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [10x10]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_inertiaJ_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_inertiaJ_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1OL_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1OL_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:04:00
% EndTime: 2020-04-20 17:04:34
% DurationCPUTime: 4.83s
% Computational Cost: add. (9218->498), mult. (8549->735), div. (0->0), fcn. (8722->18), ass. (0->274)
t246 = sin(qJ(1));
t250 = cos(qJ(1));
t242 = qJ(2) + qJ(3);
t227 = sin(t242);
t229 = cos(t242);
t281 = -Icges(4,5) * t229 + Icges(4,6) * t227;
t142 = -Icges(4,3) * t250 + t246 * t281;
t143 = Icges(4,3) * t246 + t250 * t281;
t240 = t250 ^ 2;
t336 = Icges(4,4) * t229;
t287 = Icges(4,2) * t227 - t336;
t147 = Icges(4,6) * t246 + t250 * t287;
t337 = Icges(4,4) * t227;
t293 = -Icges(4,1) * t229 + t337;
t151 = Icges(4,5) * t246 + t250 * t293;
t269 = t147 * t227 - t151 * t229;
t146 = -Icges(4,6) * t250 + t246 * t287;
t150 = -Icges(4,5) * t250 + t246 * t293;
t270 = -t146 * t227 + t150 * t229;
t231 = qJ(4) + t242;
t223 = sin(t231);
t224 = cos(t231);
t280 = -Icges(5,5) * t224 + Icges(5,6) * t223;
t124 = -Icges(5,3) * t250 + t246 * t280;
t125 = Icges(5,3) * t246 + t250 * t280;
t334 = Icges(5,4) * t224;
t286 = Icges(5,2) * t223 - t334;
t127 = Icges(5,6) * t246 + t250 * t286;
t335 = Icges(5,4) * t223;
t292 = -Icges(5,1) * t224 + t335;
t129 = Icges(5,5) * t246 + t250 * t292;
t273 = t127 * t223 - t129 * t224;
t126 = -Icges(5,6) * t250 + t246 * t286;
t128 = -Icges(5,5) * t250 + t246 * t292;
t274 = -t126 * t223 + t128 * t224;
t248 = cos(qJ(5));
t320 = t248 * t250;
t244 = sin(qJ(5));
t322 = t246 * t244;
t171 = t224 * t322 - t320;
t321 = t246 * t248;
t323 = t244 * t250;
t172 = -t224 * t321 - t323;
t326 = t223 * t246;
t72 = Icges(6,5) * t172 + Icges(6,6) * t171 - Icges(6,3) * t326;
t74 = Icges(6,4) * t172 + Icges(6,2) * t171 - Icges(6,6) * t326;
t76 = Icges(6,1) * t172 + Icges(6,4) * t171 - Icges(6,5) * t326;
t17 = t171 * t74 + t172 * t76 - t72 * t326;
t173 = t224 * t323 + t321;
t174 = -t224 * t320 + t322;
t325 = t223 * t250;
t73 = Icges(6,5) * t174 + Icges(6,6) * t173 - Icges(6,3) * t325;
t75 = Icges(6,4) * t174 + Icges(6,2) * t173 - Icges(6,6) * t325;
t77 = Icges(6,1) * t174 + Icges(6,4) * t173 - Icges(6,5) * t325;
t18 = t171 * t75 + t172 * t77 - t73 * t326;
t8 = -t17 * t250 + t18 * t246;
t360 = -t240 * t124 - (t273 * t246 + (-t125 + t274) * t250) * t246 - t8;
t367 = -t240 * t142 - (t269 * t246 + (-t143 + t270) * t250) * t246 + t360;
t239 = t246 ^ 2;
t313 = t239 + t240;
t241 = qJ(2) + qJ(7);
t226 = sin(t241);
t228 = cos(t241);
t296 = rSges(8,1) * t228 - rSges(8,2) * t226;
t366 = t246 * t250;
t357 = pkin(8) * t224;
t365 = -(pkin(10) + rSges(6,3)) * t223 - t357;
t364 = t246 / 0.2e1;
t363 = -t250 / 0.2e1;
t19 = t173 * t74 + t174 * t76 - t72 * t325;
t20 = t173 * t75 + t174 * t77 - t73 * t325;
t9 = -t19 * t250 + t20 * t246;
t361 = (t239 * t125 + t9 + (t274 * t250 + (-t124 + t273) * t246) * t250) * t246;
t245 = sin(qJ(2));
t359 = pkin(1) * t245;
t230 = pkin(15) - t241;
t358 = pkin(3) * sin(t230);
t104 = Icges(6,3) * t224 + (-Icges(6,5) * t248 + Icges(6,6) * t244) * t223;
t105 = Icges(6,6) * t224 + (-Icges(6,4) * t248 + Icges(6,2) * t244) * t223;
t356 = t223 * t244 * t105 + t224 * t104;
t249 = cos(qJ(2));
t355 = rSges(3,1) * t249;
t354 = rSges(4,1) * t229;
t353 = rSges(5,1) * t224;
t247 = cos(qJ(6));
t352 = rSges(7,1) * t247;
t350 = rSges(3,2) * t245;
t243 = sin(qJ(6));
t349 = rSges(7,2) * t243;
t347 = rSges(4,3) * t250;
t346 = rSges(5,3) * t250;
t345 = rSges(8,3) * t250;
t344 = rSges(9,3) * t250;
t343 = t250 * rSges(3,3);
t342 = t250 * rSges(7,3);
t27 = t224 * t72 + (t244 * t74 - t248 * t76) * t223;
t341 = t27 * t250;
t28 = t224 * t73 + (t244 * t75 - t248 * t77) * t223;
t340 = t28 * t246;
t222 = t249 * pkin(1) + pkin(12);
t339 = Icges(3,4) * t245;
t338 = Icges(3,4) * t249;
t333 = Icges(7,4) * t243;
t332 = Icges(7,4) * t247;
t331 = Icges(8,4) * t226;
t330 = Icges(8,4) * t228;
t225 = -qJ(8) + t230;
t217 = sin(t225);
t329 = Icges(9,4) * t217;
t218 = cos(t225);
t328 = Icges(9,4) * t218;
t106 = Icges(6,5) * t224 + (-Icges(6,1) * t248 + Icges(6,4) * t244) * t223;
t327 = t106 * t248;
t324 = t227 * t250;
t295 = -rSges(9,1) * t218 - rSges(9,2) * t217;
t252 = t246 * rSges(9,3) + t250 * t295;
t60 = t246 * (t246 * t295 - t344) + t250 * t252;
t258 = rSges(5,2) * t325 + t246 * rSges(5,3) - t250 * t353;
t299 = rSges(5,2) * t223 - t353;
t67 = t246 * (t246 * t299 - t346) + t250 * t258;
t196 = -pkin(4) * t229 + t222;
t184 = t250 * t196;
t212 = t250 * t222;
t319 = t239 * (t196 - t222) + t250 * (-t212 + t184);
t117 = rSges(6,3) * t224 + (-rSges(6,1) * t248 + rSges(6,2) * t244) * t223;
t318 = pkin(8) * t223 - pkin(10) * t224 - t117;
t257 = t246 * rSges(8,3) + t296 * t250;
t70 = t246 * (t246 * t296 - t345) + t250 * t257;
t259 = rSges(4,2) * t324 + t246 * rSges(4,3) - t250 * t354;
t300 = rSges(4,2) * t227 - t354;
t71 = t246 * (t246 * t300 - t347) + t250 * t259;
t317 = t174 * rSges(6,1) + t173 * rSges(6,2);
t316 = t239 * (-pkin(12) + t222) + t250 * (-pkin(12) * t250 + t212);
t315 = t246 * rSges(7,3) + t250 * t352;
t314 = t246 * rSges(3,3) + t250 * t355;
t312 = t246 * (t239 * t143 + (t270 * t250 + (-t142 + t269) * t246) * t250) + t361;
t170 = rSges(9,1) * t217 - rSges(9,2) * t218;
t311 = -t170 - t359;
t182 = -rSges(5,1) * t223 - rSges(5,2) * t224;
t310 = -t182 - t359;
t194 = rSges(8,1) * t226 + rSges(8,2) * t228;
t309 = -t194 - t359;
t195 = -rSges(4,1) * t227 - rSges(4,2) * t229;
t308 = -t195 - t359;
t277 = -Icges(9,5) * t218 - Icges(9,6) * t217;
t107 = -Icges(9,3) * t250 + t246 * t277;
t108 = Icges(9,3) * t246 + t250 * t277;
t283 = -Icges(9,2) * t217 - t328;
t110 = Icges(9,6) * t246 + t250 * t283;
t289 = -Icges(9,1) * t218 - t329;
t112 = Icges(9,5) * t246 + t250 * t289;
t275 = -t110 * t217 - t112 * t218;
t109 = -Icges(9,6) * t250 + t246 * t283;
t111 = -Icges(9,5) * t250 + t246 * t289;
t276 = t109 * t217 + t111 * t218;
t12 = t246 * (t239 * t108 + (t276 * t250 + (-t107 + t275) * t246) * t250);
t13 = t240 * t107 + (t275 * t246 + (-t108 + t276) * t250) * t246;
t307 = -t250 * t13 + t12;
t167 = Icges(9,5) * t217 - Icges(9,6) * t218;
t168 = -Icges(9,2) * t218 + t329;
t169 = Icges(9,1) * t217 - t328;
t264 = -t168 * t217 - t169 * t218;
t306 = (-t110 * t218 + t112 * t217 + t246 * t167 + t250 * t264) * t364 + (-t109 * t218 + t111 * t217 - t167 * t250 + t246 * t264) * t363;
t78 = t318 * t246;
t79 = t318 * t250;
t298 = -rSges(6,1) * t172 - rSges(6,2) * t171;
t80 = -rSges(6,3) * t326 - t298;
t81 = -rSges(6,3) * t325 + t317;
t30 = t246 * t80 + t250 * t81 + t313 * (-pkin(10) * t223 - t357);
t34 = -t104 * t326 + t105 * t171 + t106 * t172;
t3 = t34 * t224 + (-t17 * t246 - t18 * t250) * t223;
t35 = -t104 * t325 + t173 * t105 + t174 * t106;
t4 = t35 * t224 + (-t19 * t246 - t20 * t250) * t223;
t305 = t3 * t363 - t8 * t326 / 0.2e1 + t4 * t364 + t224 * (t340 - t341) / 0.2e1 - t9 * t325 / 0.2e1;
t193 = pkin(3) * cos(t230) + t222;
t177 = t250 * t193;
t36 = t239 * (t193 - t222) + t250 * (-t212 + t177) + t60;
t39 = t319 + t67;
t304 = t240 / 0.2e1 + t239 / 0.2e1;
t303 = t318 - t359;
t301 = -t350 + t355;
t297 = -t349 + t352;
t294 = Icges(3,1) * t249 - t339;
t291 = Icges(7,1) * t247 - t333;
t290 = Icges(8,1) * t228 - t331;
t288 = -Icges(3,2) * t245 + t338;
t285 = -Icges(7,2) * t243 + t332;
t284 = -Icges(8,2) * t226 + t330;
t282 = Icges(3,5) * t249 - Icges(3,6) * t245;
t279 = Icges(7,5) * t247 - Icges(7,6) * t243;
t278 = Icges(8,5) * t228 - Icges(8,6) * t226;
t144 = -Icges(8,6) * t250 + t246 * t284;
t148 = -Icges(8,5) * t250 + t246 * t290;
t272 = t144 * t226 - t148 * t228;
t145 = Icges(8,6) * t246 + t250 * t284;
t149 = Icges(8,5) * t246 + t250 * t290;
t271 = -t145 * t226 + t149 * t228;
t180 = -Icges(5,2) * t224 - t335;
t181 = -Icges(5,1) * t223 - t334;
t263 = t180 * t223 - t181 * t224;
t189 = Icges(8,2) * t228 + t331;
t191 = Icges(8,1) * t226 + t330;
t262 = -t189 * t226 + t191 * t228;
t190 = -Icges(4,2) * t229 - t337;
t192 = -Icges(4,1) * t227 - t336;
t261 = t190 * t227 - t192 * t229;
t260 = t360 * t250 + t361;
t140 = -Icges(8,3) * t250 + t246 * t278;
t141 = Icges(8,3) * t246 + t250 * t278;
t21 = t246 * (t239 * t141 + (t272 * t250 + (-t140 + t271) * t246) * t250);
t23 = t240 * t140 + (t271 * t246 + (-t141 + t272) * t250) * t246;
t256 = t12 + t21 + (-t13 - t23) * t250;
t14 = t30 + t319;
t179 = -Icges(5,5) * t223 - Icges(5,6) * t224;
t255 = -t341 / 0.2e1 + t340 / 0.2e1 + (-t127 * t224 - t129 * t223 + t246 * t179 + t250 * t263 + t35) * t364 + (-t126 * t224 - t128 * t223 - t179 * t250 + t246 * t263 + t34) * t363;
t187 = Icges(8,5) * t226 + Icges(8,6) * t228;
t254 = t306 + (t145 * t228 + t149 * t226 + t246 * t187 + t250 * t262) * t364 + (t144 * t228 + t148 * t226 - t187 * t250 + t246 * t262) * t363;
t253 = t367 * t250 + t312;
t188 = -Icges(4,5) * t227 - Icges(4,6) * t229;
t251 = t255 + (-t147 * t229 - t151 * t227 + t246 * t188 + t250 * t261) * t364 + (-t146 * t229 - t150 * t227 - t188 * t250 + t246 * t261) * t363;
t216 = pkin(4) * t324;
t215 = t246 * pkin(4) * t227;
t211 = rSges(2,1) * t250 - t246 * rSges(2,2);
t210 = -t246 * rSges(2,1) - rSges(2,2) * t250;
t209 = rSges(3,1) * t245 + rSges(3,2) * t249;
t208 = rSges(7,1) * t243 + rSges(7,2) * t247;
t199 = t250 * t358;
t197 = t246 * t358;
t157 = Icges(3,3) * t246 + t250 * t282;
t156 = -Icges(3,3) * t250 + t246 * t282;
t155 = Icges(7,3) * t246 + t250 * t279;
t154 = -Icges(7,3) * t250 + t246 * t279;
t139 = (-pkin(6) - t349) * t250 + t315;
t138 = (pkin(12) - t350) * t250 + t314;
t137 = t342 + (pkin(6) - t297) * t246;
t136 = t343 + (-pkin(12) - t301) * t246;
t133 = t308 * t250;
t132 = t309 * t250;
t131 = t308 * t246;
t130 = t309 * t246;
t119 = -t182 * t250 + t216;
t118 = -t182 * t246 + t215;
t99 = t212 + t259;
t98 = t212 + t257;
t97 = t347 + (-t222 - t300) * t246;
t96 = t345 + (-t222 - t296) * t246;
t95 = -t170 * t250 + t199;
t94 = -t170 * t246 + t197;
t93 = t250 * t310 + t216;
t92 = t246 * t310 + t215;
t89 = t250 * t311 + t199;
t88 = t246 * t311 + t197;
t87 = t184 + t258;
t86 = t346 + (-t196 - t299) * t246;
t85 = t177 + t252;
t84 = t344 + (-t193 - t295) * t246;
t83 = t250 * (-t250 * t350 + t314) + (t246 * t301 - t343) * t246;
t82 = t250 * (-t250 * t349 + t315) + (t246 * t297 - t342) * t246;
t66 = t216 + t79;
t65 = t215 + t78;
t55 = t250 * t303 + t216;
t54 = t246 * t303 + t215;
t45 = t365 * t250 + t184 + t317;
t44 = (-t196 - t365) * t246 + t298;
t43 = t316 + t71;
t42 = t316 + t70;
t41 = t117 * t325 + t224 * t81;
t40 = -t117 * t326 - t224 * t80;
t38 = (-t223 * t327 + t356) * t224;
t37 = (t246 * t81 - t250 * t80) * t223;
t31 = t39 + t316;
t29 = t36 + t316;
t11 = t14 + t316;
t1 = [t243 * (Icges(7,1) * t243 + t332) + t247 * (Icges(7,2) * t247 + t333) + t245 * (Icges(3,1) * t245 + t338) + t249 * (Icges(3,2) * t249 + t339) + (-t181 - t327) * t223 + m(4) * (t97 ^ 2 + t99 ^ 2) + m(3) * (t136 ^ 2 + t138 ^ 2) + m(2) * (t210 ^ 2 + t211 ^ 2) + m(9) * (t84 ^ 2 + t85 ^ 2) + m(8) * (t96 ^ 2 + t98 ^ 2) + m(7) * (t137 ^ 2 + t139 ^ 2) + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t86 ^ 2 + t87 ^ 2) + t356 - t227 * t192 + t228 * t189 - t229 * t190 - t224 * t180 + t226 * t191 + t217 * t169 - t218 * t168 + Icges(2,3); m(3) * (-t136 * t250 - t138 * t246) * t209 + ((-Icges(3,6) * t250 + t246 * t288) * t249 + (-Icges(3,5) * t250 + t246 * t294) * t245) * t363 + t251 + t304 * (Icges(3,5) * t245 + Icges(3,6) * t249) + ((Icges(3,6) * t246 + t250 * t288) * t249 + (Icges(3,5) * t246 + t250 * t294) * t245) * t364 + t254 + m(8) * (t130 * t98 + t132 * t96) + m(4) * (t131 * t99 + t133 * t97) + m(5) * (t86 * t93 + t87 * t92) + m(9) * (t84 * t89 + t85 * t88) + m(6) * (t44 * t55 + t45 * t54); t21 + m(8) * (t130 ^ 2 + t132 ^ 2 + t42 ^ 2) + m(6) * (t11 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(9) * (t29 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(3) * (t209 ^ 2 * t313 + t83 ^ 2) + t246 * t239 * t157 + m(5) * (t31 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(4) * (t131 ^ 2 + t133 ^ 2 + t43 ^ 2) + t307 + t312 + (-t240 * t156 - t23 + (-t246 * t156 + t250 * t157) * t246 + t367) * t250; m(6) * (t44 * t66 + t45 * t65) + m(5) * (t118 * t87 + t119 * t86) + m(4) * (-t246 * t99 - t250 * t97) * t195 + t251; m(5) * (t118 * t92 + t119 * t93 + t31 * t39) + m(4) * (t71 * t43 + (-t131 * t246 - t133 * t250) * t195) + m(6) * (t11 * t14 + t54 * t65 + t55 * t66) + t253; m(6) * (t14 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t118 ^ 2 + t119 ^ 2 + t39 ^ 2) + m(4) * (t195 ^ 2 * t313 + t71 ^ 2) + t253; m(6) * (t44 * t79 + t45 * t78) + m(5) * (-t246 * t87 - t250 * t86) * t182 + t255; m(5) * (t67 * t31 + (-t246 * t92 - t250 * t93) * t182) + m(6) * (t11 * t30 + t54 * t78 + t55 * t79) + t260; m(6) * (t14 * t30 + t65 * t78 + t66 * t79) + m(5) * (t67 * t39 + (-t118 * t246 - t119 * t250) * t182) + t260; m(6) * (t30 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t182 ^ 2 * t313 + t67 ^ 2) + t260; t38 + m(6) * (t40 * t44 + t41 * t45) + ((-t28 / 0.2e1 - t35 / 0.2e1) * t250 + (-t27 / 0.2e1 - t34 / 0.2e1) * t246) * t223; m(6) * (t11 * t37 + t40 * t55 + t41 * t54) + t305; m(6) * (t14 * t37 + t40 * t66 + t41 * t65) + t305; m(6) * (t30 * t37 + t40 * t79 + t41 * t78) + t305; t224 * t38 + m(6) * (t37 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t224 * (-t246 * t27 - t250 * t28) - t250 * t4 - t246 * t3) * t223; (t247 * (-Icges(7,6) * t250 + t246 * t285) + t243 * (-Icges(7,5) * t250 + t246 * t291)) * t363 + (t247 * (Icges(7,6) * t246 + t250 * t285) + t243 * (Icges(7,5) * t246 + t250 * t291)) * t364 + m(7) * (-t137 * t250 - t139 * t246) * t208 + t304 * (Icges(7,5) * t243 + Icges(7,6) * t247); 0; 0; 0; 0; t246 * (-t154 * t366 + t239 * t155) - t250 * (t240 * t154 - t155 * t366) + m(7) * (t208 ^ 2 * t313 + t82 ^ 2); m(9) * (t84 * t95 + t85 * t94) + m(8) * (-t246 * t98 - t250 * t96) * t194 + t254; m(8) * (t70 * t42 + (-t130 * t246 - t132 * t250) * t194) + m(9) * (t29 * t36 + t88 * t94 + t89 * t95) + t256; 0; 0; 0; 0; m(9) * (t36 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(8) * (t194 ^ 2 * t313 + t70 ^ 2) + t256; m(9) * (-t246 * t85 - t250 * t84) * t170 + t306; m(9) * (t60 * t29 + (-t246 * t88 - t250 * t89) * t170) + t307; 0; 0; 0; 0; m(9) * (t60 * t36 + (-t246 * t94 - t250 * t95) * t170) + t307; m(9) * (t170 ^ 2 * t313 + t60 ^ 2) + t307; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_10_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55);];
Mq = res;
