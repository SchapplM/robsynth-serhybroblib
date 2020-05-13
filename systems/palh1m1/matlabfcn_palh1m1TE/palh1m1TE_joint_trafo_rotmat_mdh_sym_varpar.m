% Calculate homogenous joint transformation matrices for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_mdh [4x4x16]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh1m1TE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:08:59
% EndTime: 2020-04-12 20:09:00
% DurationCPUTime: 1.22s
% Computational Cost: add. (13522->109), mult. (20282->177), div. (992->9), fcn. (12812->28), ass. (0->145)
t374 = -pkin(2) - pkin(13);
t373 = -pkin(2) + pkin(13);
t372 = -pkin(8) - pkin(3);
t371 = -pkin(8) + pkin(3);
t370 = -pkin(9) - pkin(11);
t369 = -pkin(9) + pkin(11);
t289 = sin(pkin(23));
t368 = t289 / 0.2e1;
t291 = sin(pkin(21));
t367 = -t291 / 0.2e1;
t366 = t291 / 0.2e1;
t292 = sin(pkin(20));
t365 = t292 / 0.2e1;
t293 = cos(pkin(23));
t364 = -t293 / 0.2e1;
t363 = t293 / 0.2e1;
t294 = cos(pkin(22));
t362 = -t294 / 0.2e1;
t295 = cos(pkin(21));
t361 = -t295 / 0.2e1;
t360 = t295 / 0.2e1;
t304 = cos(qJ(3));
t359 = -t304 / 0.2e1;
t358 = cos(pkin(18)) / 0.2e1;
t315 = pkin(8) ^ 2;
t321 = pkin(3) ^ 2;
t299 = sin(qJ(2));
t301 = sin(pkin(19));
t305 = cos(qJ(2));
t307 = cos(pkin(19));
t284 = t299 * t307 - t305 * t301;
t355 = pkin(7) * t284;
t334 = pkin(1) * t355;
t281 = -0.2e1 * t334;
t325 = pkin(1) ^ 2;
t339 = pkin(7) ^ 2 + t325;
t331 = t281 + t339;
t273 = -t315 + t321 + t331;
t279 = pkin(1) - t355;
t340 = t281 + t325;
t267 = sqrt(-((pkin(7) - t371) * (pkin(7) + t371) + t340) * ((pkin(7) - t372) * (pkin(7) + t372) + t340));
t285 = t299 * t301 + t305 * t307;
t350 = t267 * t285;
t261 = -pkin(7) * t350 + t279 * t273;
t262 = pkin(7) * t285 * t273 + t279 * t267;
t298 = sin(qJ(3));
t275 = 0.1e1 / t331;
t322 = 0.1e1 / pkin(3);
t346 = t275 * t322;
t252 = (-t298 * t261 / 0.2e1 + t262 * t359) * t346;
t253 = (t261 * t359 + t298 * t262 / 0.2e1) * t346;
t288 = pkin(23) + pkin(22);
t286 = sin(t288);
t287 = cos(t288);
t246 = t287 * t252 + t286 * t253;
t357 = pkin(5) * t246;
t296 = cos(pkin(20));
t282 = -t304 * t292 - t298 * t296;
t356 = pkin(6) * t282;
t336 = pkin(4) * t357;
t245 = -0.2e1 * t336;
t319 = pkin(5) ^ 2;
t342 = t245 + t319;
t238 = sqrt(-((pkin(4) - t369) * (pkin(4) + t369) + t342) * ((pkin(4) - t370) * (pkin(4) + t370) + t342));
t247 = -t286 * t252 + t287 * t253;
t354 = t238 * t247;
t337 = pkin(4) ^ 2 + t319;
t333 = t245 + t337;
t242 = 0.1e1 / t333;
t312 = 0.1e1 / pkin(11);
t353 = t242 * t312;
t314 = 0.1e1 / pkin(9);
t352 = t242 * t314;
t335 = pkin(1) * t356;
t278 = -0.2e1 * t335;
t318 = pkin(6) ^ 2;
t341 = t278 + t318;
t266 = sqrt(-((pkin(1) - t373) * (pkin(1) + t373) + t341) * ((pkin(1) - t374) * (pkin(1) + t374) + t341));
t283 = t298 * t292 - t304 * t296;
t351 = t266 * t283;
t338 = t318 + t325;
t332 = t278 + t338;
t274 = 0.1e1 / t332;
t310 = 0.1e1 / pkin(13);
t349 = t274 * t310;
t324 = 0.1e1 / pkin(2);
t348 = t274 * t324;
t316 = 0.1e1 / pkin(8);
t347 = t275 * t316;
t345 = t310 * t324;
t344 = t312 * t314;
t343 = t316 * t322;
t330 = -t321 + t339;
t323 = pkin(2) ^ 2;
t329 = -t323 + t338;
t313 = pkin(9) ^ 2;
t328 = -t313 + t337;
t327 = t348 / 0.2e1;
t326 = -t345 / 0.2e1;
t311 = pkin(11) ^ 2;
t309 = pkin(13) ^ 2;
t306 = cos(qJ(1));
t303 = cos(qJ(4));
t302 = sin(pkin(18));
t300 = sin(qJ(1));
t297 = sin(qJ(4));
t290 = sin(pkin(22));
t280 = pkin(1) * t284 - pkin(7);
t277 = -pkin(1) * t282 + pkin(6);
t276 = -pkin(1) + t356;
t272 = t281 + t315 + t330;
t271 = t315 - t330 + 0.2e1 * t334;
t270 = t278 + t309 + t329;
t269 = -t309 + t323 + t332;
t268 = (t309 - t329 + 0.2e1 * t335) * t326;
t265 = (t267 * t368 + t271 * t363) * t343;
t264 = (t267 * t364 + t271 * t368) * t343;
t263 = pkin(1) * t285 * t272 - t280 * t267;
t260 = -pkin(1) * t350 - t280 * t272;
t259 = pkin(1) * t283 * t270 + t277 * t266;
t258 = pkin(6) * t283 * t269 - t276 * t266;
t257 = -pkin(1) * t351 + t277 * t270;
t256 = (-pkin(6) * t351 - t276 * t269) * t327;
t255 = (t263 * t358 + t260 * t302 / 0.2e1) * t347;
t254 = (t260 * t358 - t302 * t263 / 0.2e1) * t347;
t251 = (t261 * t368 + t262 * t363) * t346;
t250 = (t261 * t364 + t262 * t368) * t346;
t249 = (t296 * t257 / 0.2e1 + t259 * t365) * t349;
t248 = (t257 * t365 - t296 * t259 / 0.2e1) * t349;
t244 = -pkin(4) * t246 + pkin(5);
t243 = -pkin(4) + t357;
t241 = t245 + t311 + t328;
t240 = -t311 + t313 + t333;
t239 = t311 - t328 + 0.2e1 * t336;
t237 = (t238 * t367 + t239 * t361) * t344;
t236 = (t238 * t360 + t239 * t367) * t344;
t235 = pkin(4) * t247 * t241 + t244 * t238;
t234 = pkin(5) * t247 * t240 - t243 * t238;
t233 = -pkin(4) * t354 + t244 * t241;
t232 = -pkin(5) * t354 - t243 * t240;
t231 = (t233 * t361 + t235 * t366) * t353;
t230 = (t233 * t366 + t235 * t360) * t353;
t229 = (t232 * t362 - t290 * t234 / 0.2e1) * t352;
t228 = (t290 * t232 / 0.2e1 + t234 * t362) * t352;
t1 = [t306, -t300, 0, 0; t300, t306, 0, 0; 0, 0, 1, pkin(14); 0, 0, 0, 1; -t299, -t305, 0, pkin(16); 0, 0, -1, 0; t305, -t299, 0, 0; 0, 0, 0, 1; t298, t304, 0, pkin(1); -t304, t298, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t231, -t230, 0, pkin(5); t230, t231, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t303, -t297, 0, pkin(10); 0, 0, -1, -pkin(12); t297, t303, 0, 0; 0, 0, 0, 1; t254, -t255, 0, -pkin(15); 0, 0, -1, 0; t255, t254, 0, -pkin(17); 0, 0, 0, 1; t250, -t251, 0, pkin(1); t251, t250, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t256, -t258 * t348 / 0.2e1, 0, 0; t258 * t327, t256, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t268, t266 * t345 / 0.2e1, 0, pkin(2); t266 * t326, t268, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t229, -t228, 0, t294 * pkin(4); t228, t229, 0, -t290 * pkin(4); 0, 0, 1, 0; 0, 0, 0, 1; t265, -t264, 0, t293 * pkin(3); t264, t265, 0, t289 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; t249, -t248, 0, t296 * pkin(6); t248, t249, 0, t292 * pkin(6); 0, 0, 1, 0; 0, 0, 0, 1; t237, -t236, 0, t295 * pkin(11); t236, t237, 0, t291 * pkin(11); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(8); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(13); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(9); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,16);             % numerisch
else,                         T_mdh = sym('xx', [4,4,16]); end % symbolisch

for i = 1:16
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
