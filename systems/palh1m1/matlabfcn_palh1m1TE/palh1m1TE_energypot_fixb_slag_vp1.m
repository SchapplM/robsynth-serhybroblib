% Calculate potential energy for
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1TE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1TE_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:10:28
% EndTime: 2020-04-12 20:10:34
% DurationCPUTime: 6.28s
% Computational Cost: add. (80051->229), mult. (120425->343), div. (5724->9), fcn. (76231->28), ass. (0->140)
t325 = -pkin(12) - rSges(6,3);
t235 = sin(pkin(20));
t237 = cos(pkin(20));
t239 = sin(qJ(3));
t244 = cos(qJ(3));
t322 = pkin(6) * (-t235 * t244 - t237 * t239);
t305 = pkin(1) * t322;
t221 = -0.2e1 * t305;
t252 = pkin(6) ^ 2;
t310 = t221 + t252;
t330 = pkin(13) - pkin(2);
t335 = -pkin(2) - pkin(13);
t207 = sqrt(-((pkin(1) - t330) * (pkin(1) + t330) + t310) * ((pkin(1) - t335) * (pkin(1) + t335) + t310));
t253 = pkin(2) ^ 2;
t255 = pkin(1) ^ 2;
t308 = t252 + t255;
t298 = -pkin(13) ^ 2 + t308;
t210 = t221 + t253 + t298;
t220 = -pkin(1) + t322;
t321 = pkin(6) * (t235 * t239 - t237 * t244);
t205 = -t207 * t321 - t210 * t220;
t206 = -t207 * t220 + t210 * t321;
t245 = cos(qJ(2));
t254 = 0.1e1 / pkin(2);
t314 = 0.1e1 / (t221 + t308) * t254;
t240 = sin(qJ(2));
t326 = -t240 / 0.2e1;
t337 = (t245 * t205 / 0.2e1 + t206 * t326) * t314;
t336 = 0.1e1 / pkin(3);
t334 = pkin(3) - pkin(8);
t333 = -pkin(8) - pkin(3);
t332 = -pkin(9) - pkin(11);
t331 = -pkin(9) + pkin(11);
t329 = -t207 / 0.2e1;
t328 = t207 / 0.2e1;
t327 = -t253 / 0.2e1 + t298 / 0.2e1 - t305;
t324 = sin(pkin(18));
t323 = pkin(1) * t240;
t242 = sin(pkin(19));
t247 = cos(pkin(19));
t227 = t240 * t247 - t242 * t245;
t320 = pkin(7) * t227;
t319 = t245 * pkin(1) + pkin(14);
t318 = cos(pkin(21));
t317 = cos(pkin(23));
t316 = sin(pkin(21));
t315 = sin(pkin(23));
t309 = -0.2e1 * pkin(1) * t320 + t255;
t299 = pkin(7) ^ 2 + t309;
t213 = 0.1e1 / t299;
t313 = t213 / pkin(8);
t208 = sqrt(-((pkin(7) - t334) * (pkin(7) + t334) + t309) * ((pkin(7) - t333) * (pkin(7) + t333) + t309));
t230 = t240 * t242 + t245 * t247;
t312 = t230 * t208;
t311 = 0.1e1 / pkin(13) * t254;
t307 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t306 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t304 = pkin(23) + pkin(22);
t303 = cos(pkin(18)) / 0.2e1;
t302 = rSges(10,1) * t311;
t301 = rSges(10,2) * t311;
t228 = t239 * t245 + t240 * t244;
t300 = t228 * pkin(5) + t319;
t297 = pkin(1) - t320;
t296 = cos(t304);
t295 = sin(t304);
t241 = sin(qJ(1));
t231 = t241 * pkin(16);
t294 = -t241 * t323 + t231;
t246 = cos(qJ(1));
t232 = t246 * pkin(16);
t293 = -t246 * t323 + t232;
t292 = -rSges(3,1) * t240 - rSges(3,2) * t245;
t291 = t299 - t307;
t282 = t336 * (-pkin(7) * t312 + t291 * t297);
t280 = -t282 / 0.2e1;
t283 = t336 * (pkin(7) * t230 * t291 + t208 * t297);
t281 = t283 / 0.2e1;
t197 = (t280 * t317 + t281 * t315) * t213;
t198 = (t317 * t281 + t315 * t282 / 0.2e1) * t213;
t189 = t197 * t245 - t198 * t240;
t190 = -t197 * t240 - t198 * t245;
t229 = -t239 * t240 + t244 * t245;
t216 = t229 * t241;
t290 = t216 * pkin(5) + t294;
t218 = t229 * t246;
t289 = t218 * pkin(5) + t293;
t211 = t299 + t307;
t222 = pkin(1) * t227 - pkin(7);
t286 = -pkin(1) * t312 - t211 * t222;
t287 = pkin(1) * t211 * t230 - t208 * t222;
t199 = (t286 * t303 - t324 * t287 / 0.2e1) * t313;
t200 = (t287 * t303 + t286 * t324 / 0.2e1) * t313;
t288 = rSges(7,1) * t199 - rSges(7,2) * t200 - pkin(15);
t196 = (-t245 * t206 / 0.2e1 + t205 * t326) * t314;
t279 = t239 * t280 - t244 * t283 / 0.2e1;
t278 = t239 * t281 + t244 * t280;
t277 = t213 * (t278 * t295 + t279 * t296);
t276 = t213 * (t278 * t296 - t279 * t295);
t275 = pkin(5) * t277;
t274 = pkin(4) - t275;
t273 = -pkin(4) * t277 + pkin(5);
t272 = -0.2e1 * pkin(4) * t275 + pkin(5) ^ 2;
t271 = pkin(4) ^ 2 + t272;
t270 = 0.1e1 / t271;
t269 = 0.1e1 / pkin(9) * t270;
t268 = 0.1e1 / pkin(11) * t270;
t267 = t271 + t306;
t266 = t271 - t306;
t265 = sqrt(-((pkin(4) - t331) * (pkin(4) + t331) + t272) * ((pkin(4) - t332) * (pkin(4) + t332) + t272));
t264 = t265 * t276;
t263 = (-pkin(5) * t264 + t266 * t274) * t269;
t262 = (-pkin(4) * t264 + t267 * t273) * t268;
t261 = -(pkin(5) * t266 * t276 + t265 * t274) * t269 / 0.2e1;
t260 = (pkin(4) * t267 * t276 + t265 * t273) * t268 / 0.2e1;
t243 = cos(qJ(4));
t238 = sin(qJ(4));
t236 = cos(pkin(22));
t234 = sin(pkin(22));
t219 = t228 * t246;
t217 = t228 * t241;
t194 = t246 * t337;
t193 = t246 * t196;
t192 = t241 * t337;
t191 = t241 * t196;
t188 = t189 * t246;
t187 = t190 * t246;
t186 = t189 * t241;
t185 = t190 * t241;
t180 = -t318 * t262 / 0.2e1 + t316 * t260;
t179 = t316 * t262 / 0.2e1 + t318 * t260;
t178 = -t236 * t263 / 0.2e1 + t234 * t261;
t177 = t234 * t263 / 0.2e1 + t236 * t261;
t176 = t179 * t229 + t180 * t228;
t175 = -t179 * t228 + t180 * t229;
t174 = -t179 * t219 + t180 * t218;
t173 = -t179 * t218 - t180 * t219;
t172 = -t179 * t217 + t180 * t216;
t171 = -t179 * t216 - t180 * t217;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t246 - rSges(2,2) * t241) + g(2) * (rSges(2,1) * t241 + rSges(2,2) * t246) + g(3) * (pkin(14) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t241 + t246 * t292 + t232) + g(2) * (-rSges(3,3) * t246 + t241 * t292 + t231) + g(3) * (rSges(3,1) * t245 - rSges(3,2) * t240 + pkin(14))) - m(4) * (g(1) * (rSges(4,1) * t218 - rSges(4,2) * t219 + rSges(4,3) * t241 + t293) + g(2) * (rSges(4,1) * t216 - rSges(4,2) * t217 - rSges(4,3) * t246 + t294) + g(3) * (rSges(4,1) * t228 + rSges(4,2) * t229 + t319)) - m(5) * (g(1) * (rSges(5,1) * t174 + rSges(5,2) * t173 + rSges(5,3) * t241 + t289) + g(2) * (rSges(5,1) * t172 + rSges(5,2) * t171 - rSges(5,3) * t246 + t290) + g(3) * (rSges(5,1) * t176 + rSges(5,2) * t175 + t300)) - m(6) * (g(1) * (t174 * pkin(10) + (t174 * t243 + t238 * t241) * rSges(6,1) + (-t174 * t238 + t241 * t243) * rSges(6,2) + t325 * t173 + t289) + g(2) * (t172 * pkin(10) + (t172 * t243 - t238 * t246) * rSges(6,1) + (-t172 * t238 - t243 * t246) * rSges(6,2) + t325 * t171 + t290) + (t300 + (t243 * rSges(6,1) - t238 * rSges(6,2) + pkin(10)) * t176 + t325 * t175) * g(3)) - m(7) * (g(3) * (rSges(7,1) * t200 + rSges(7,2) * t199 + pkin(14) - pkin(17)) + (-g(2) * rSges(7,3) + g(1) * t288) * t246 + (g(1) * rSges(7,3) + g(2) * t288) * t241) - m(8) * (g(1) * (rSges(8,1) * t187 - rSges(8,2) * t188 + rSges(8,3) * t241 + t293) + g(2) * (rSges(8,1) * t185 - rSges(8,2) * t186 - rSges(8,3) * t246 + t294) + g(3) * (rSges(8,1) * t189 + rSges(8,2) * t190 + t319)) - m(9) * (g(1) * (rSges(9,1) * t193 - rSges(9,2) * t194 + rSges(9,3) * t241 + t232) + g(2) * (rSges(9,1) * t191 - rSges(9,2) * t192 - rSges(9,3) * t246 + t231) + g(3) * (rSges(9,1) * t337 + rSges(9,2) * t196 + pkin(14))) - m(10) * (g(1) * (t193 * pkin(2) + t232 + (t193 * t327 - t194 * t329) * t302 + (t193 * t328 - t194 * t327) * t301 + t241 * rSges(10,3)) + g(2) * (t191 * pkin(2) + t231 + (t191 * t327 - t192 * t329) * t302 + (t191 * t328 - t192 * t327) * t301 - t246 * rSges(10,3)) + g(3) * (t337 * pkin(2) + pkin(14) + (t196 * t329 + t327 * t337) * t302 + (t196 * t327 + t328 * t337) * t301)) - m(11) * (g(1) * ((-t177 * t188 + t178 * t187) * rSges(11,1) + (-t177 * t187 - t178 * t188) * rSges(11,2) + t241 * rSges(11,3) + t293) + g(2) * ((-t177 * t186 + t178 * t185) * rSges(11,1) + (-t177 * t185 - t178 * t186) * rSges(11,2) - t246 * rSges(11,3) + t294) + g(3) * ((t177 * t190 + t178 * t189) * rSges(11,1) + (-t177 * t189 + t178 * t190) * rSges(11,2) + t319) + (g(1) * (t187 * t236 + t188 * t234) + g(2) * (t185 * t236 + t186 * t234) + g(3) * (t189 * t236 - t190 * t234)) * pkin(4));
U = t1;
