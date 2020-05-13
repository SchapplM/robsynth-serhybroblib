% Calculate inertial parameters regressor of potential energy for
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m1DE2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energypot_fixb_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 05:40:08
% EndTime: 2020-04-20 05:40:12
% DurationCPUTime: 3.86s
% Computational Cost: add. (69629->102), mult. (104487->154), div. (4968->6), fcn. (66105->38), ass. (0->94)
t293 = sin(qJ(1));
t299 = cos(qJ(1));
t315 = g(1) * t299 + g(2) * t293;
t350 = -pkin(6) - pkin(2);
t349 = -pkin(6) + pkin(2);
t348 = -pkin(8) - pkin(10);
t347 = -pkin(8) + pkin(10);
t346 = g(3) * pkin(12);
t292 = sin(qJ(2));
t294 = sin(pkin(16));
t298 = cos(qJ(2));
t300 = cos(pkin(16));
t273 = t292 * t294 - t298 * t300;
t338 = pkin(5) * t273;
t322 = (-0.2e1 * t338 + pkin(1)) * pkin(1);
t318 = pkin(5) ^ 2 + t322;
t320 = pkin(2) ^ 2 - pkin(6) ^ 2;
t268 = t318 + t320;
t270 = pkin(1) - t338;
t265 = sqrt(-((pkin(5) - t349) * (pkin(5) + t349) + t322) * ((pkin(5) - t350) * (pkin(5) + t350) + t322));
t274 = t292 * t300 + t294 * t298;
t330 = t265 * t274;
t262 = -pkin(5) * t330 + t268 * t270;
t345 = -t262 / 0.2e1;
t263 = pkin(5) * t268 * t274 + t265 * t270;
t344 = t263 / 0.2e1;
t343 = sin(pkin(17)) / 0.2e1;
t342 = sin(pkin(19)) / 0.2e1;
t341 = sin(qJ(3)) / 0.2e1;
t340 = cos(pkin(15)) / 0.2e1;
t297 = cos(qJ(3));
t269 = 0.1e1 / t318;
t328 = t269 / pkin(2);
t259 = (t263 * t341 + t297 * t345) * t328;
t260 = (t262 * t341 + t297 * t344) * t328;
t284 = pkin(18) + pkin(19);
t278 = sin(t284);
t279 = cos(t284);
t250 = -t259 * t279 - t260 * t278;
t339 = pkin(4) * t250;
t323 = -0.2e1 * pkin(3) * t339 + pkin(4) ^ 2;
t319 = pkin(3) ^ 2 + t323;
t321 = pkin(8) ^ 2 - pkin(10) ^ 2;
t243 = t319 - t321;
t246 = -pkin(3) * t250 + pkin(4);
t241 = sqrt(-((pkin(3) - t347) * (pkin(3) + t347) + t323) * ((pkin(3) - t348) * (pkin(3) + t348) + t323));
t249 = t259 * t278 - t260 * t279;
t332 = t241 * t249;
t239 = -pkin(3) * t332 + t243 * t246;
t240 = pkin(3) * t243 * t249 + t241 * t246;
t285 = qJ(2) + qJ(3);
t287 = cos(pkin(17));
t244 = 0.1e1 / t319;
t331 = t244 / pkin(10);
t235 = atan2((t240 * t287 / 0.2e1 + t239 * t343) * t331, (-t239 * t287 / 0.2e1 + t240 * t343) * t331) + t285;
t233 = sin(t235);
t335 = g(3) * t233;
t334 = t292 * pkin(1) + pkin(12);
t333 = t298 * pkin(1) + pkin(13);
t329 = t269 / pkin(6);
t290 = sin(qJ(4));
t327 = t293 * t290;
t296 = cos(qJ(4));
t326 = t293 * t296;
t325 = t299 * t290;
t324 = t299 * t296;
t289 = cos(pkin(19));
t254 = qJ(2) + atan2((t262 * t342 + t289 * t344) * t328, (t263 * t342 + t289 * t345) * t328);
t317 = t244 / pkin(8) / 0.2e1;
t253 = pkin(18) - t254;
t280 = sin(t285);
t316 = -pkin(4) * t280 + t334;
t295 = sin(pkin(15));
t281 = cos(t285);
t276 = -g(1) * t293 + g(2) * t299;
t275 = -pkin(4) * t281 + t333;
t271 = pkin(1) * t273 - pkin(5);
t267 = t318 - t320;
t266 = -g(3) * t334 - t315 * t333;
t264 = pkin(1) * t267 * t274 - t265 * t271;
t261 = -pkin(1) * t330 - t267 * t271;
t258 = atan2((t264 * t340 - t261 * t295 / 0.2e1) * t329, (t261 * t340 + t264 * t295 / 0.2e1) * t329);
t256 = cos(t258);
t255 = sin(t258);
t252 = cos(t254);
t251 = sin(t254);
t245 = -pkin(3) + t339;
t242 = t319 + t321;
t238 = -atan2((pkin(4) * t242 * t249 - t241 * t245) * t317, (-pkin(4) * t332 - t242 * t245) * t317) + t253;
t237 = cos(t238);
t236 = sin(t238);
t234 = cos(t235);
t232 = -g(3) * t234 + t233 * t315;
t1 = [0, 0, 0, 0, 0, 0, -t315, -t276, -g(3), -t346, 0, 0, 0, 0, 0, 0, -g(3) * t292 - t298 * t315, -g(3) * t298 + t292 * t315, t276, -pkin(13) * t315 - t346, 0, 0, 0, 0, 0, 0, g(3) * t280 + t281 * t315, g(3) * t281 - t280 * t315, t276, t266, 0, 0, 0, 0, 0, 0, t234 * t315 + t335, -t232, t276, -g(3) * t316 - t275 * t315, 0, 0, 0, 0, 0, 0, -g(2) * (-t234 * t326 - t325) - g(1) * (-t234 * t324 + t327) + t296 * t335, -g(2) * (t234 * t327 - t324) - g(1) * (t234 * t325 + t326) - t290 * t335, t232, -g(3) * (-pkin(9) * t233 + pkin(11) * t234 + t316) + t315 * (pkin(9) * t234 + pkin(11) * t233 - t275), 0, 0, 0, 0, 0, 0, -g(3) * t255 - t256 * t315, -g(3) * t256 + t255 * t315, t276, -g(3) * (pkin(14) + pkin(12)) + t315 * pkin(7), 0, 0, 0, 0, 0, 0, -g(3) * t251 - t252 * t315, -g(3) * t252 + t251 * t315, t276, t266, 0, 0, 0, 0, 0, 0, -g(3) * t236 + t237 * t315, g(3) * t237 + t236 * t315, t276, -g(3) * (-pkin(3) * sin(t253) + t334) - t315 * (pkin(3) * cos(t253) + t333);];
U_reg = t1;
