% Calculate inertial parameters regressor of potential energy for
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m1TE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energypot_fixb_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 00:30:20
% EndTime: 2020-04-18 00:30:26
% DurationCPUTime: 5.13s
% Computational Cost: add. (78713->132), mult. (118414->223), div. (5616->6), fcn. (74977->24), ass. (0->107)
t286 = cos(qJ(2));
t357 = pkin(1) * t286 + pkin(13);
t356 = 0.1e1 / pkin(2);
t355 = -pkin(6) - pkin(2);
t354 = -pkin(6) + pkin(2);
t353 = -pkin(8) - pkin(10);
t352 = -pkin(8) + pkin(10);
t351 = g(3) * pkin(12);
t350 = sin(pkin(15));
t281 = sin(qJ(2));
t283 = sin(pkin(16));
t288 = cos(pkin(16));
t268 = t281 * t283 - t286 * t288;
t348 = pkin(5) * t268;
t339 = (-0.2e1 * t348 + pkin(1)) * pkin(1);
t329 = pkin(5) ^ 2 + t339;
t257 = 0.1e1 / t329;
t280 = sin(qJ(3));
t285 = cos(qJ(3));
t336 = -pkin(2) ^ 2 + pkin(6) ^ 2;
t324 = t329 - t336;
t328 = pkin(1) - t348;
t254 = sqrt(-((pkin(5) - t354) * (pkin(5) + t354) + t339) * ((pkin(5) - t355) * (pkin(5) + t355) + t339));
t270 = t281 * t288 + t283 * t286;
t340 = t270 * t254;
t319 = t356 * (-pkin(5) * t340 + t328 * t324);
t317 = -t319 / 0.2e1;
t318 = t356 * (pkin(5) * t270 * t324 + t328 * t254) / 0.2e1;
t314 = t280 * t318 + t285 * t317;
t316 = t319 / 0.2e1;
t315 = t280 * t316 + t285 * t318;
t334 = pkin(18) + pkin(19);
t326 = sin(t334);
t327 = cos(t334);
t313 = t257 * (-t327 * t314 - t326 * t315);
t311 = pkin(4) * t313;
t308 = -0.2e1 * pkin(3) * t311 + pkin(4) ^ 2;
t301 = sqrt(-((pkin(3) - t352) * (pkin(3) + t352) + t308) * ((pkin(3) - t353) * (pkin(3) + t353) + t308));
t307 = pkin(3) ^ 2 + t308;
t335 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t303 = t307 + t335;
t306 = 0.1e1 / t307;
t304 = 0.1e1 / pkin(10) * t306;
t309 = -pkin(3) * t313 + pkin(4);
t312 = (t326 * t314 - t327 * t315) * t257;
t296 = (pkin(3) * t303 * t312 + t309 * t301) * t304 / 0.2e1;
t300 = t301 * t312;
t298 = (-pkin(3) * t300 + t309 * t303) * t304;
t342 = sin(pkin(17));
t343 = cos(pkin(17));
t232 = -t343 * t298 / 0.2e1 + t342 * t296;
t233 = t343 * t296 + t342 * t298 / 0.2e1;
t267 = t280 * t281 - t285 * t286;
t323 = t280 * t286 + t281 * t285;
t230 = -t232 * t323 + t233 * t267;
t347 = g(3) * t230;
t346 = cos(pkin(19));
t345 = sin(pkin(19));
t344 = t281 * pkin(1) + pkin(12);
t341 = t257 / pkin(6);
t282 = sin(qJ(1));
t338 = t357 * t282;
t287 = cos(qJ(1));
t337 = t357 * t287;
t333 = cos(pkin(15)) / 0.2e1;
t332 = -pkin(4) * t323 + t344;
t261 = t267 * t282;
t331 = t261 * pkin(4) + t338;
t263 = t267 * t287;
t330 = t263 * pkin(4) + t337;
t325 = g(1) * t287 + g(2) * t282;
t246 = (t346 * t317 + t345 * t318) * t257;
t247 = (t345 * t316 + t346 * t318) * t257;
t245 = t246 * t286 - t247 * t281;
t244 = t246 * t281 + t247 * t286;
t256 = t329 + t336;
t264 = pkin(1) * t268 - pkin(5);
t322 = pkin(1) * t256 * t270 - t254 * t264;
t321 = -pkin(1) * t340 - t256 * t264;
t260 = t323 * t282;
t227 = t232 * t260 - t233 * t261;
t262 = t323 * t287;
t229 = t232 * t262 - t233 * t263;
t231 = t232 * t267 + t233 * t323;
t320 = g(1) * t229 + g(2) * t227 + g(3) * t231;
t255 = -g(1) * t337 - g(2) * t338 - g(3) * t344;
t310 = pkin(3) - t311;
t305 = 0.1e1 / pkin(8) * t306;
t302 = t307 - t335;
t299 = (-pkin(4) * t300 + t310 * t302) * t305;
t297 = -(pkin(4) * t302 * t312 + t310 * t301) * t305 / 0.2e1;
t284 = cos(qJ(4));
t279 = sin(qJ(4));
t278 = cos(pkin(18));
t277 = sin(pkin(18));
t271 = -g(1) * t282 + g(2) * t287;
t249 = (t322 * t333 - t321 * t350 / 0.2e1) * t341;
t248 = (t321 * t333 + t322 * t350 / 0.2e1) * t341;
t243 = t244 * t287;
t242 = t245 * t287;
t241 = t244 * t282;
t240 = t245 * t282;
t235 = -t278 * t299 / 0.2e1 + t277 * t297;
t234 = t277 * t299 / 0.2e1 + t278 * t297;
t228 = t232 * t263 + t233 * t262;
t226 = t232 * t261 + t233 * t260;
t1 = [0, 0, 0, 0, 0, 0, -t325, -t271, -g(3), -t351, 0, 0, 0, 0, 0, 0, -g(3) * t281 - t325 * t286, -g(3) * t286 + t325 * t281, t271, -t325 * pkin(13) - t351, 0, 0, 0, 0, 0, 0, -g(1) * t263 - g(2) * t261 + g(3) * t323, -g(1) * t262 - g(2) * t260 - g(3) * t267, t271, t255, 0, 0, 0, 0, 0, 0, -g(1) * t228 - g(2) * t226 - t347, -t320, t271, -g(1) * t330 - g(2) * t331 - g(3) * t332, 0, 0, 0, 0, 0, 0, -t284 * t347 - g(2) * (t226 * t284 - t279 * t287) - g(1) * (t228 * t284 + t279 * t282), t279 * t347 - g(2) * (-t226 * t279 - t284 * t287) - g(1) * (-t228 * t279 + t282 * t284), t320, -g(3) * (pkin(9) * t230 - pkin(11) * t231 + t332) - g(2) * (pkin(9) * t226 - pkin(11) * t227 + t331) - g(1) * (pkin(9) * t228 - pkin(11) * t229 + t330), 0, 0, 0, 0, 0, 0, -g(3) * t249 - t325 * t248, -g(3) * t248 + t325 * t249, t271, -g(3) * (pkin(14) + pkin(12)) + t325 * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t242 - g(2) * t240 - g(3) * t244, g(1) * t243 + g(2) * t241 - g(3) * t245, t271, t255, 0, 0, 0, 0, 0, 0, -g(3) * (t234 * t245 + t235 * t244) - g(2) * (-t234 * t241 + t235 * t240) - g(1) * (-t234 * t243 + t235 * t242), -g(3) * (-t234 * t244 + t235 * t245) - g(2) * (-t234 * t240 - t235 * t241) - g(1) * (-t234 * t242 - t235 * t243), t271, (-g(3) * (t244 * t278 - t245 * t277) - g(2) * (t240 * t278 + t241 * t277) - g(1) * (t242 * t278 + t243 * t277)) * pkin(3) + t255;];
U_reg = t1;
