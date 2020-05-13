% Calculate inertial parameters regressor of potential energy for
% palh3m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m1DE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energypot_fixb_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 05:12:55
% EndTime: 2020-04-19 05:13:03
% DurationCPUTime: 6.87s
% Computational Cost: add. (157217->132), mult. (236326->220), div. (11232->6), fcn. (149713->36), ass. (0->112)
t319 = cos(qJ(2));
t373 = pkin(1) * t319 + pkin(13);
t372 = -pkin(6) - pkin(2);
t371 = -pkin(6) + pkin(2);
t370 = -pkin(8) - pkin(10);
t369 = -pkin(8) + pkin(10);
t368 = g(3) * pkin(12);
t313 = sin(qJ(2));
t315 = sin(pkin(16));
t321 = cos(pkin(16));
t293 = t313 * t315 - t319 * t321;
t359 = pkin(5) * t293;
t350 = (-0.2e1 * t359 + pkin(1)) * pkin(1);
t341 = pkin(5) ^ 2 + t350;
t346 = pkin(2) ^ 2 - pkin(6) ^ 2;
t280 = t341 + t346;
t288 = pkin(1) - t359;
t277 = sqrt(-((pkin(5) - t371) * (pkin(5) + t371) + t350) * ((pkin(5) - t372) * (pkin(5) + t372) + t350));
t295 = t313 * t321 + t319 * t315;
t354 = t277 * t295;
t274 = -pkin(5) * t354 + t288 * t280;
t367 = -t274 / 0.2e1;
t275 = pkin(5) * t295 * t280 + t288 * t277;
t366 = t275 / 0.2e1;
t365 = sin(pkin(17)) / 0.2e1;
t364 = sin(pkin(19)) / 0.2e1;
t312 = sin(qJ(3));
t363 = t312 / 0.2e1;
t362 = cos(pkin(15)) / 0.2e1;
t318 = cos(qJ(3));
t281 = 0.1e1 / t341;
t352 = t281 / pkin(2);
t271 = (t275 * t363 + t318 * t367) * t352;
t272 = (t274 * t363 + t318 * t366) * t352;
t304 = pkin(18) + pkin(19);
t299 = sin(t304);
t300 = cos(t304);
t264 = -t300 * t271 - t299 * t272;
t360 = pkin(4) * t264;
t351 = -0.2e1 * pkin(3) * t360 + pkin(4) ^ 2;
t344 = pkin(3) ^ 2 + t351;
t347 = pkin(8) ^ 2 - pkin(10) ^ 2;
t258 = t344 - t347;
t261 = -pkin(3) * t264 + pkin(4);
t250 = sqrt(-((pkin(3) - t369) * (pkin(3) + t369) + t351) * ((pkin(3) - t370) * (pkin(3) + t370) + t351));
t263 = t299 * t271 - t300 * t272;
t356 = t250 * t263;
t248 = -pkin(3) * t356 + t261 * t258;
t249 = pkin(3) * t263 * t258 + t261 * t250;
t306 = cos(pkin(17));
t259 = 0.1e1 / t344;
t355 = t259 / pkin(10);
t245 = atan2((t249 * t306 / 0.2e1 + t248 * t365) * t355, (-t248 * t306 / 0.2e1 + t249 * t365) * t355);
t243 = sin(t245);
t244 = cos(t245);
t292 = t313 * t312 - t319 * t318;
t338 = t319 * t312 + t313 * t318;
t237 = t292 * t243 - t244 * t338;
t358 = g(3) * t237;
t357 = t313 * pkin(1) + pkin(12);
t353 = t281 / pkin(6);
t314 = sin(qJ(1));
t349 = t373 * t314;
t320 = cos(qJ(1));
t348 = t373 * t320;
t345 = -pkin(4) * t338 + t357;
t285 = t292 * t314;
t343 = t285 * pkin(4) + t349;
t287 = t292 * t320;
t342 = t287 * pkin(4) + t348;
t340 = t259 / pkin(8) / 0.2e1;
t339 = g(1) * t320 + g(2) * t314;
t309 = cos(pkin(19));
t269 = atan2((t274 * t364 + t309 * t366) * t352, (t275 * t364 + t309 * t367) * t352);
t265 = sin(t269);
t266 = cos(t269);
t256 = t319 * t265 + t313 * t266;
t255 = -t313 * t265 + t319 * t266;
t284 = t338 * t314;
t232 = t285 * t243 - t284 * t244;
t286 = t338 * t320;
t234 = t287 * t243 - t286 * t244;
t236 = -t243 * t338 - t292 * t244;
t337 = g(1) * t234 + g(2) * t232 + g(3) * t236;
t278 = -g(1) * t348 - g(2) * t349 - g(3) * t357;
t257 = t344 + t347;
t260 = -pkin(3) + t360;
t336 = atan2((pkin(4) * t263 * t257 - t260 * t250) * t340, (-pkin(4) * t356 - t260 * t257) * t340);
t335 = sin(t336);
t317 = cos(qJ(4));
t316 = sin(pkin(15));
t311 = sin(qJ(4));
t310 = cos(pkin(18));
t308 = sin(pkin(18));
t296 = -g(1) * t314 + g(2) * t320;
t289 = pkin(1) * t293 - pkin(5);
t279 = t341 - t346;
t276 = pkin(1) * t295 * t279 - t289 * t277;
t273 = -pkin(1) * t354 - t289 * t279;
t270 = atan2((t276 * t362 - t273 * t316 / 0.2e1) * t353, (t273 * t362 + t276 * t316 / 0.2e1) * t353);
t268 = cos(t270);
t267 = sin(t270);
t254 = t255 * t320;
t253 = t256 * t320;
t252 = t255 * t314;
t251 = t256 * t314;
t247 = cos(t336);
t242 = -t310 * t247 - t308 * t335;
t241 = t308 * t247 - t310 * t335;
t235 = t286 * t243 + t287 * t244;
t233 = t284 * t243 + t285 * t244;
t1 = [0, 0, 0, 0, 0, 0, -t339, -t296, -g(3), -t368, 0, 0, 0, 0, 0, 0, -g(3) * t313 - t319 * t339, -g(3) * t319 + t313 * t339, t296, -pkin(13) * t339 - t368, 0, 0, 0, 0, 0, 0, -g(1) * t287 - g(2) * t285 + g(3) * t338, -g(1) * t286 - g(2) * t284 - g(3) * t292, t296, t278, 0, 0, 0, 0, 0, 0, -g(1) * t235 - g(2) * t233 - t358, t337, t296, -g(1) * t342 - g(2) * t343 - g(3) * t345, 0, 0, 0, 0, 0, 0, -t317 * t358 - g(2) * (t233 * t317 - t320 * t311) - g(1) * (t235 * t317 + t314 * t311), t311 * t358 - g(2) * (-t233 * t311 - t320 * t317) - g(1) * (-t235 * t311 + t314 * t317), -t337, -g(3) * (t237 * pkin(9) + t236 * pkin(11) + t345) - g(2) * (t233 * pkin(9) + t232 * pkin(11) + t343) - g(1) * (t235 * pkin(9) + t234 * pkin(11) + t342), 0, 0, 0, 0, 0, 0, -g(3) * t267 - t268 * t339, -g(3) * t268 + t267 * t339, t296, -g(3) * (pkin(14) + pkin(12)) + t339 * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t254 - g(2) * t252 - g(3) * t256, g(1) * t253 + g(2) * t251 - g(3) * t255, t296, t278, 0, 0, 0, 0, 0, 0, -g(3) * (t255 * t241 + t256 * t242) - g(2) * (-t251 * t241 + t252 * t242) - g(1) * (-t253 * t241 + t254 * t242), -g(3) * (-t256 * t241 + t255 * t242) - g(2) * (-t252 * t241 - t251 * t242) - g(1) * (-t254 * t241 - t253 * t242), t296, (-g(3) * (-t255 * t308 + t256 * t310) - g(2) * (t251 * t308 + t252 * t310) - g(1) * (t253 * t308 + t254 * t310)) * pkin(3) + t278;];
U_reg = t1;
