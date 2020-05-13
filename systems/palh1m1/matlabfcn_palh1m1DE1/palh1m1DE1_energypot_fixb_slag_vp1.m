% Calculate potential energy for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE1_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:47:45
% EndTime: 2020-04-13 14:47:53
% DurationCPUTime: 7.70s
% Computational Cost: add. (159809->229), mult. (240221->333), div. (11448->9), fcn. (152251->46), ass. (0->143)
t336 = pkin(12) + rSges(6,3);
t349 = -pkin(2) - pkin(13);
t348 = -pkin(2) + pkin(13);
t347 = -pkin(3) - pkin(8);
t346 = -pkin(8) + pkin(3);
t345 = -pkin(9) - pkin(11);
t344 = -pkin(9) + pkin(11);
t300 = pkin(1) ^ 2;
t274 = sin(qJ(2));
t276 = sin(pkin(19));
t280 = cos(qJ(2));
t282 = cos(pkin(19));
t254 = t274 * t282 - t276 * t280;
t331 = pkin(7) * t254;
t322 = -0.2e1 * pkin(1) * t331 + t300;
t315 = pkin(7) ^ 2 + t322;
t319 = pkin(3) ^ 2 - pkin(8) ^ 2;
t237 = t315 + t319;
t248 = pkin(1) - t331;
t234 = sqrt(-((pkin(7) - t346) * (pkin(7) + t346) + t322) * ((pkin(7) - t347) * (pkin(7) + t347) + t322));
t257 = t274 * t276 + t280 * t282;
t327 = t234 * t257;
t227 = -pkin(7) * t327 + t237 * t248;
t343 = -t227 / 0.2e1;
t228 = pkin(7) * t237 * t257 + t234 * t248;
t342 = t228 / 0.2e1;
t341 = sin(pkin(23)) / 0.2e1;
t340 = sin(pkin(21)) / 0.2e1;
t279 = cos(qJ(3));
t339 = -t279 / 0.2e1;
t338 = cos(pkin(18)) / 0.2e1;
t337 = 0.1e1 / pkin(2) / 0.2e1;
t335 = pkin(1) * t274;
t273 = sin(qJ(3));
t239 = 0.1e1 / t315;
t325 = t239 / pkin(3);
t224 = (t228 * t339 + t273 * t343) * t325;
t225 = (t227 * t339 + t273 * t342) * t325;
t263 = pkin(23) + pkin(22);
t258 = sin(t263);
t259 = cos(t263);
t213 = t224 * t259 + t225 * t258;
t334 = pkin(5) * t213;
t267 = sin(pkin(20));
t271 = cos(pkin(20));
t333 = pkin(6) * (-t267 * t279 - t271 * t273);
t332 = pkin(6) * (t267 * t273 - t271 * t279);
t330 = t280 * pkin(1) + pkin(14);
t324 = -0.2e1 * pkin(4) * t334 + pkin(5) ^ 2;
t194 = sqrt(-((pkin(4) - t344) * (pkin(4) + t344) + t324) * ((pkin(4) - t345) * (pkin(4) + t345) + t324));
t214 = -t224 * t258 + t225 * t259;
t329 = t194 * t214;
t316 = pkin(4) ^ 2 + t324;
t209 = 0.1e1 / t316;
t328 = t209 / pkin(11);
t326 = t239 / pkin(8);
t318 = pkin(1) * t333;
t247 = -0.2e1 * t318;
t293 = pkin(6) ^ 2;
t323 = t247 + t293;
t321 = pkin(9) ^ 2 - pkin(11) ^ 2;
t320 = t293 + t300;
t255 = t273 * t280 + t274 * t279;
t317 = t255 * pkin(5) + t330;
t314 = -pkin(13) ^ 2 + t320;
t313 = t209 / pkin(9) / 0.2e1;
t312 = 0.1e1 / (t247 + t320) * t337;
t311 = 0.1e1 / pkin(13) * t337;
t275 = sin(qJ(1));
t260 = t275 * pkin(16);
t310 = -t275 * t335 + t260;
t281 = cos(qJ(1));
t261 = t281 * pkin(16);
t309 = -t281 * t335 + t261;
t308 = -rSges(3,1) * t274 - rSges(3,2) * t280;
t268 = cos(pkin(23));
t219 = atan2((t227 * t341 + t268 * t342) * t325, (t228 * t341 + t268 * t343) * t325);
t215 = sin(t219);
t216 = cos(t219);
t199 = -t215 * t280 - t216 * t274;
t307 = t215 * t274 - t216 * t280;
t233 = sqrt(-((pkin(1) - t348) * (pkin(1) + t348) + t323) * ((pkin(1) - t349) * (pkin(1) + t349) + t323));
t298 = pkin(2) ^ 2;
t235 = t247 + t298 + t314;
t246 = -pkin(1) + t333;
t223 = atan2((-t233 * t246 + t235 * t332) * t312, (-t233 * t332 - t235 * t246) * t312);
t221 = sin(t223);
t222 = cos(t223);
t207 = -t221 * t280 - t222 * t274;
t306 = t221 * t274 - t222 * t280;
t256 = -t273 * t274 + t279 * t280;
t242 = t256 * t275;
t305 = t242 * pkin(5) + t310;
t244 = t256 * t281;
t304 = t244 * pkin(5) + t309;
t236 = t315 - t319;
t249 = pkin(1) * t254 - pkin(7);
t226 = -pkin(1) * t327 - t236 * t249;
t229 = pkin(1) * t236 * t257 - t234 * t249;
t277 = sin(pkin(18));
t220 = atan2((t229 * t338 + t226 * t277 / 0.2e1) * t326, (t226 * t338 - t277 * t229 / 0.2e1) * t326);
t217 = sin(t220);
t218 = cos(t220);
t303 = rSges(7,1) * t218 - rSges(7,2) * t217 - pkin(15);
t205 = t316 + t321;
t210 = -pkin(4) + t334;
t302 = atan2((pkin(5) * t205 * t214 - t194 * t210) * t313, (-pkin(5) * t329 - t205 * t210) * t313);
t301 = sin(t302);
t278 = cos(qJ(4));
t272 = sin(qJ(4));
t270 = cos(pkin(21));
t269 = cos(pkin(22));
t265 = sin(pkin(22));
t245 = t255 * t281;
t243 = t255 * t275;
t232 = atan2(t233 * t311, (t298 - t314 + 0.2e1 * t318) * t311);
t231 = cos(t232);
t230 = sin(t232);
t211 = -pkin(4) * t213 + pkin(5);
t206 = t316 - t321;
t204 = t207 * t281;
t203 = t306 * t281;
t202 = t207 * t275;
t201 = t306 * t275;
t198 = t199 * t281;
t197 = t307 * t281;
t196 = t199 * t275;
t195 = t307 * t275;
t193 = pkin(4) * t206 * t214 + t194 * t211;
t192 = -pkin(4) * t329 + t206 * t211;
t191 = cos(t302);
t189 = atan2((t192 * t340 + t193 * t270 / 0.2e1) * t328, (-t192 * t270 / 0.2e1 + t193 * t340) * t328);
t188 = cos(t189);
t187 = sin(t189);
t186 = -t269 * t191 - t265 * t301;
t185 = t191 * t265 - t269 * t301;
t181 = t187 * t256 + t188 * t255;
t180 = t187 * t255 - t256 * t188;
t179 = -t187 * t245 + t188 * t244;
t178 = t187 * t244 + t245 * t188;
t177 = -t187 * t243 + t188 * t242;
t176 = t187 * t242 + t243 * t188;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t281 - rSges(2,2) * t275) + g(2) * (rSges(2,1) * t275 + rSges(2,2) * t281) + g(3) * (pkin(14) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t275 + t281 * t308 + t261) + g(2) * (-rSges(3,3) * t281 + t275 * t308 + t260) + g(3) * (rSges(3,1) * t280 - rSges(3,2) * t274 + pkin(14))) - m(4) * (g(1) * (rSges(4,1) * t244 - rSges(4,2) * t245 + rSges(4,3) * t275 + t309) + g(2) * (rSges(4,1) * t242 - rSges(4,2) * t243 - rSges(4,3) * t281 + t310) + g(3) * (rSges(4,1) * t255 + rSges(4,2) * t256 + t330)) - m(5) * (g(1) * (rSges(5,1) * t179 - rSges(5,2) * t178 + rSges(5,3) * t275 + t304) + g(2) * (rSges(5,1) * t177 - rSges(5,2) * t176 - rSges(5,3) * t281 + t305) + g(3) * (rSges(5,1) * t181 - rSges(5,2) * t180 + t317)) - m(6) * (g(1) * (t179 * pkin(10) + (t179 * t278 + t272 * t275) * rSges(6,1) + (-t179 * t272 + t275 * t278) * rSges(6,2) + t336 * t178 + t304) + g(2) * (t177 * pkin(10) + (t177 * t278 - t272 * t281) * rSges(6,1) + (-t177 * t272 - t278 * t281) * rSges(6,2) + t336 * t176 + t305) + (t317 + (t278 * rSges(6,1) - t272 * rSges(6,2) + pkin(10)) * t181 + t336 * t180) * g(3)) - m(7) * (g(3) * (rSges(7,1) * t217 + rSges(7,2) * t218 + pkin(14) - pkin(17)) + (-g(2) * rSges(7,3) + g(1) * t303) * t281 + (g(1) * rSges(7,3) + g(2) * t303) * t275) - m(8) * (g(1) * (rSges(8,1) * t198 + rSges(8,2) * t197 + rSges(8,3) * t275 + t309) + g(2) * (rSges(8,1) * t196 + rSges(8,2) * t195 - rSges(8,3) * t281 + t310) + g(3) * (-rSges(8,1) * t307 + rSges(8,2) * t199 + t330)) - m(9) * (g(1) * (rSges(9,1) * t204 + rSges(9,2) * t203 + rSges(9,3) * t275 + t261) + g(2) * (rSges(9,1) * t202 + rSges(9,2) * t201 - rSges(9,3) * t281 + t260) + g(3) * (-rSges(9,1) * t306 + rSges(9,2) * t207 + pkin(14))) - m(10) * (g(1) * (t204 * pkin(2) + t261 + (-t203 * t230 - t204 * t231) * rSges(10,1) + (-t203 * t231 + t204 * t230) * rSges(10,2) + t275 * rSges(10,3)) + g(2) * (t202 * pkin(2) + t260 + (-t201 * t230 - t202 * t231) * rSges(10,1) + (-t201 * t231 + t202 * t230) * rSges(10,2) - t281 * rSges(10,3)) + g(3) * (-t306 * pkin(2) + pkin(14) + (-t207 * t230 + t231 * t306) * rSges(10,1) + (-t207 * t231 - t230 * t306) * rSges(10,2))) - m(11) * (g(1) * ((t185 * t197 + t186 * t198) * rSges(11,1) + (-t185 * t198 + t186 * t197) * rSges(11,2) + t275 * rSges(11,3) + t309) + g(2) * ((t185 * t195 + t186 * t196) * rSges(11,1) + (-t185 * t196 + t186 * t195) * rSges(11,2) - t281 * rSges(11,3) + t310) + g(3) * ((t185 * t199 - t186 * t307) * rSges(11,1) + (t185 * t307 + t186 * t199) * rSges(11,2) + t330) + (g(1) * (-t197 * t265 + t198 * t269) + g(2) * (-t195 * t265 + t196 * t269) + g(3) * (-t199 * t265 - t269 * t307)) * pkin(4));
U = t1;
