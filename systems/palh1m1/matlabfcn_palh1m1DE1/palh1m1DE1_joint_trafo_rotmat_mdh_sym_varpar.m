% Calculate homogenous joint transformation matrices for
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh1m1DE1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:45:58
% EndTime: 2020-04-13 14:46:00
% DurationCPUTime: 1.96s
% Computational Cost: add. (27025->113), mult. (40524->168), div. (1984->9), fcn. (25672->55), ass. (0->150)
t380 = -pkin(2) - pkin(13);
t379 = -pkin(2) + pkin(13);
t378 = -pkin(8) - pkin(3);
t377 = -pkin(8) + pkin(3);
t376 = -pkin(9) - pkin(11);
t375 = -pkin(9) + pkin(11);
t327 = pkin(8) ^ 2;
t333 = pkin(3) ^ 2;
t311 = sin(qJ(2));
t313 = sin(pkin(19));
t317 = cos(qJ(2));
t319 = cos(pkin(19));
t296 = t311 * t319 - t317 * t313;
t365 = pkin(7) * t296;
t347 = pkin(1) * t365;
t293 = -0.2e1 * t347;
t337 = pkin(1) ^ 2;
t352 = pkin(7) ^ 2 + t337;
t344 = t293 + t352;
t285 = -t327 + t333 + t344;
t291 = pkin(1) - t365;
t353 = t293 + t337;
t281 = sqrt(-((pkin(7) - t377) * (pkin(7) + t377) + t353) * ((pkin(7) - t378) * (pkin(7) + t378) + t353));
t297 = t311 * t313 + t317 * t319;
t361 = t281 * t297;
t271 = -pkin(7) * t361 + t291 * t285;
t374 = -t271 / 0.2e1;
t272 = pkin(7) * t297 * t285 + t291 * t281;
t373 = t272 / 0.2e1;
t301 = sin(pkin(23));
t372 = t301 / 0.2e1;
t303 = sin(pkin(21));
t371 = t303 / 0.2e1;
t316 = cos(qJ(3));
t370 = -t316 / 0.2e1;
t369 = cos(pkin(18)) / 0.2e1;
t368 = 0.1e1 / pkin(2) / 0.2e1;
t310 = sin(qJ(3));
t287 = 0.1e1 / t344;
t334 = 0.1e1 / pkin(3);
t358 = t287 * t334;
t266 = (t272 * t370 + t310 * t374) * t358;
t267 = (t271 * t370 + t310 * t373) * t358;
t300 = pkin(23) + pkin(22);
t298 = sin(t300);
t299 = cos(t300);
t252 = t299 * t266 + t298 * t267;
t367 = pkin(5) * t252;
t304 = sin(pkin(20));
t308 = cos(pkin(20));
t294 = -t316 * t304 - t310 * t308;
t366 = pkin(6) * t294;
t349 = pkin(4) * t367;
t251 = -0.2e1 * t349;
t331 = pkin(5) ^ 2;
t355 = t251 + t331;
t243 = sqrt(-((pkin(4) - t375) * (pkin(4) + t375) + t355) * ((pkin(4) - t376) * (pkin(4) + t376) + t355));
t253 = -t298 * t266 + t299 * t267;
t364 = t243 * t253;
t350 = pkin(4) ^ 2 + t331;
t346 = t251 + t350;
t248 = 0.1e1 / t346;
t324 = 0.1e1 / pkin(11);
t363 = t248 * t324;
t348 = pkin(1) * t366;
t290 = -0.2e1 * t348;
t330 = pkin(6) ^ 2;
t354 = t290 + t330;
t280 = sqrt(-((pkin(1) - t379) * (pkin(1) + t379) + t354) * ((pkin(1) - t380) * (pkin(1) + t380) + t354));
t295 = t310 * t304 - t316 * t308;
t362 = t280 * t295;
t351 = t330 + t337;
t345 = t290 + t351;
t286 = 0.1e1 / t345;
t322 = 0.1e1 / pkin(13);
t360 = t286 * t322;
t328 = 0.1e1 / pkin(8);
t359 = t287 * t328;
t326 = 0.1e1 / pkin(9);
t357 = t324 * t326;
t356 = t328 * t334;
t343 = -t333 + t352;
t335 = pkin(2) ^ 2;
t342 = -t335 + t351;
t325 = pkin(9) ^ 2;
t341 = -t325 + t350;
t340 = t248 * t326 / 0.2e1;
t339 = t286 * t368;
t338 = t322 * t368;
t323 = pkin(11) ^ 2;
t321 = pkin(13) ^ 2;
t318 = cos(qJ(1));
t315 = cos(qJ(4));
t314 = sin(pkin(18));
t312 = sin(qJ(1));
t309 = sin(qJ(4));
t307 = cos(pkin(21));
t306 = cos(pkin(22));
t305 = cos(pkin(23));
t302 = sin(pkin(22));
t292 = pkin(1) * t296 - pkin(7);
t289 = -pkin(1) * t294 + pkin(6);
t288 = -pkin(1) + t366;
t284 = t293 + t327 + t343;
t283 = t290 + t321 + t342;
t282 = -t321 + t335 + t345;
t279 = atan2(t281 * t356 / 0.2e1, -(t327 - t343 + 0.2e1 * t347) * t356 / 0.2e1);
t278 = cos(t279);
t277 = sin(t279);
t276 = atan2(t280 * t338, (t321 - t342 + 0.2e1 * t348) * t338);
t275 = cos(t276);
t274 = sin(t276);
t273 = pkin(1) * t297 * t284 - t292 * t281;
t270 = -pkin(1) * t361 - t292 * t284;
t269 = t301 * t277 - t305 * t278;
t268 = -t305 * t277 - t301 * t278;
t265 = atan2((pkin(1) * t295 * t283 + t289 * t280) * t360 / 0.2e1, -(-pkin(1) * t362 + t289 * t283) * t360 / 0.2e1);
t264 = atan2((pkin(6) * t295 * t282 - t288 * t280) * t339, (-pkin(6) * t362 - t288 * t282) * t339);
t263 = cos(t265);
t262 = cos(t264);
t261 = sin(t265);
t260 = sin(t264);
t259 = atan2((t273 * t369 + t270 * t314 / 0.2e1) * t359, (t270 * t369 - t314 * t273 / 0.2e1) * t359);
t258 = atan2((t271 * t372 + t305 * t373) * t358, (t272 * t372 + t305 * t374) * t358);
t257 = cos(t259);
t256 = sin(t259);
t255 = cos(t258);
t254 = sin(t258);
t250 = -pkin(4) * t252 + pkin(5);
t249 = -pkin(4) + t367;
t247 = t304 * t261 - t308 * t263;
t246 = -t308 * t261 - t304 * t263;
t245 = t251 + t323 + t341;
t244 = -t323 + t325 + t346;
t242 = atan2(t243 * t357 / 0.2e1, -(t323 - t341 + 0.2e1 * t349) * t357 / 0.2e1);
t241 = cos(t242);
t240 = sin(t242);
t239 = pkin(4) * t253 * t245 + t250 * t243;
t238 = -pkin(4) * t364 + t250 * t245;
t237 = -t303 * t240 + t307 * t241;
t236 = t307 * t240 + t303 * t241;
t235 = atan2((pkin(5) * t253 * t244 - t249 * t243) * t340, (-pkin(5) * t364 - t249 * t244) * t340);
t234 = cos(t235);
t233 = sin(t235);
t232 = atan2((t238 * t371 + t239 * t307 / 0.2e1) * t363, (-t238 * t307 / 0.2e1 + t239 * t371) * t363);
t231 = cos(t232);
t230 = sin(t232);
t229 = -t302 * t233 - t306 * t234;
t228 = -t306 * t233 + t302 * t234;
t1 = [t318, -t312, 0, 0; t312, t318, 0, 0; 0, 0, 1, pkin(14); 0, 0, 0, 1; -t311, -t317, 0, pkin(16); 0, 0, -1, 0; t317, -t311, 0, 0; 0, 0, 0, 1; t310, t316, 0, pkin(1); -t316, t310, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t231, -t230, 0, pkin(5); t230, t231, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t315, -t309, 0, pkin(10); 0, 0, -1, -pkin(12); t309, t315, 0, 0; 0, 0, 0, 1; t257, -t256, 0, -pkin(15); 0, 0, -1, 0; t256, t257, 0, -pkin(17); 0, 0, 0, 1; t255, -t254, 0, pkin(1); t254, t255, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t262, -t260, 0, 0; t260, t262, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t275, t274, 0, pkin(2); -t274, -t275, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t229, -t228, 0, t306 * pkin(4); t228, t229, 0, -t302 * pkin(4); 0, 0, 1, 0; 0, 0, 0, 1; t269, -t268, 0, t305 * pkin(3); t268, t269, 0, t301 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; t247, -t246, 0, t308 * pkin(6); t246, t247, 0, t304 * pkin(6); 0, 0, 1, 0; 0, 0, 0, 1; t237, -t236, 0, t307 * pkin(11); t236, t237, 0, t303 * pkin(11); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(8); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(13); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(9); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,16);             % numerisch
else,                         T_mdh = sym('xx', [4,4,16]); end % symbolisch

for i = 1:16
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
