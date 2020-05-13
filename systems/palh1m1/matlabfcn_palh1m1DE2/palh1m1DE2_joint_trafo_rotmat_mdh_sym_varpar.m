% Calculate homogenous joint transformation matrices for
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh1m1DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:00:34
% EndTime: 2020-04-14 20:00:36
% DurationCPUTime: 2.08s
% Computational Cost: add. (27025->113), mult. (40524->168), div. (1984->9), fcn. (25672->55), ass. (0->150)
t355 = -pkin(2) - pkin(13);
t354 = -pkin(2) + pkin(13);
t353 = -pkin(8) - pkin(3);
t352 = -pkin(8) + pkin(3);
t351 = -pkin(9) - pkin(11);
t350 = -pkin(9) + pkin(11);
t302 = pkin(8) ^ 2;
t308 = pkin(3) ^ 2;
t286 = sin(qJ(2));
t288 = sin(pkin(19));
t292 = cos(qJ(2));
t294 = cos(pkin(19));
t271 = t286 * t294 - t292 * t288;
t340 = pkin(7) * t271;
t322 = pkin(1) * t340;
t268 = -0.2e1 * t322;
t312 = pkin(1) ^ 2;
t327 = pkin(7) ^ 2 + t312;
t319 = t268 + t327;
t260 = -t302 + t308 + t319;
t266 = pkin(1) - t340;
t328 = t268 + t312;
t256 = sqrt(-((pkin(7) - t352) * (pkin(7) + t352) + t328) * ((pkin(7) - t353) * (pkin(7) + t353) + t328));
t272 = t286 * t288 + t292 * t294;
t336 = t256 * t272;
t246 = -pkin(7) * t336 + t266 * t260;
t349 = -t246 / 0.2e1;
t247 = pkin(7) * t272 * t260 + t266 * t256;
t348 = t247 / 0.2e1;
t276 = sin(pkin(23));
t347 = t276 / 0.2e1;
t278 = sin(pkin(21));
t346 = t278 / 0.2e1;
t291 = cos(qJ(3));
t345 = -t291 / 0.2e1;
t344 = cos(pkin(18)) / 0.2e1;
t343 = 0.1e1 / pkin(2) / 0.2e1;
t285 = sin(qJ(3));
t262 = 0.1e1 / t319;
t309 = 0.1e1 / pkin(3);
t333 = t262 * t309;
t241 = (t247 * t345 + t285 * t349) * t333;
t242 = (t246 * t345 + t285 * t348) * t333;
t275 = pkin(23) + pkin(22);
t273 = sin(t275);
t274 = cos(t275);
t227 = t274 * t241 + t273 * t242;
t342 = pkin(5) * t227;
t279 = sin(pkin(20));
t283 = cos(pkin(20));
t269 = -t291 * t279 - t285 * t283;
t341 = pkin(6) * t269;
t324 = pkin(4) * t342;
t226 = -0.2e1 * t324;
t306 = pkin(5) ^ 2;
t330 = t226 + t306;
t218 = sqrt(-((pkin(4) - t350) * (pkin(4) + t350) + t330) * ((pkin(4) - t351) * (pkin(4) + t351) + t330));
t228 = -t273 * t241 + t274 * t242;
t339 = t218 * t228;
t325 = pkin(4) ^ 2 + t306;
t321 = t226 + t325;
t223 = 0.1e1 / t321;
t299 = 0.1e1 / pkin(11);
t338 = t223 * t299;
t323 = pkin(1) * t341;
t265 = -0.2e1 * t323;
t305 = pkin(6) ^ 2;
t329 = t265 + t305;
t255 = sqrt(-((pkin(1) - t354) * (pkin(1) + t354) + t329) * ((pkin(1) - t355) * (pkin(1) + t355) + t329));
t270 = t285 * t279 - t291 * t283;
t337 = t255 * t270;
t326 = t305 + t312;
t320 = t265 + t326;
t261 = 0.1e1 / t320;
t297 = 0.1e1 / pkin(13);
t335 = t261 * t297;
t303 = 0.1e1 / pkin(8);
t334 = t262 * t303;
t301 = 0.1e1 / pkin(9);
t332 = t299 * t301;
t331 = t303 * t309;
t318 = -t308 + t327;
t310 = pkin(2) ^ 2;
t317 = -t310 + t326;
t300 = pkin(9) ^ 2;
t316 = -t300 + t325;
t315 = t223 * t301 / 0.2e1;
t314 = t261 * t343;
t313 = t297 * t343;
t298 = pkin(11) ^ 2;
t296 = pkin(13) ^ 2;
t293 = cos(qJ(1));
t290 = cos(qJ(4));
t289 = sin(pkin(18));
t287 = sin(qJ(1));
t284 = sin(qJ(4));
t282 = cos(pkin(21));
t281 = cos(pkin(22));
t280 = cos(pkin(23));
t277 = sin(pkin(22));
t267 = pkin(1) * t271 - pkin(7);
t264 = -pkin(1) * t269 + pkin(6);
t263 = -pkin(1) + t341;
t259 = t268 + t302 + t318;
t258 = t265 + t296 + t317;
t257 = -t296 + t310 + t320;
t254 = atan2(t256 * t331 / 0.2e1, -(t302 - t318 + 0.2e1 * t322) * t331 / 0.2e1);
t253 = cos(t254);
t252 = sin(t254);
t251 = atan2(t255 * t313, (t296 - t317 + 0.2e1 * t323) * t313);
t250 = cos(t251);
t249 = sin(t251);
t248 = pkin(1) * t272 * t259 - t267 * t256;
t245 = -pkin(1) * t336 - t267 * t259;
t244 = t276 * t252 - t280 * t253;
t243 = -t280 * t252 - t276 * t253;
t240 = atan2((pkin(1) * t270 * t258 + t264 * t255) * t335 / 0.2e1, -(-pkin(1) * t337 + t264 * t258) * t335 / 0.2e1);
t239 = atan2((pkin(6) * t270 * t257 - t263 * t255) * t314, (-pkin(6) * t337 - t263 * t257) * t314);
t238 = cos(t240);
t237 = cos(t239);
t236 = sin(t240);
t235 = sin(t239);
t234 = atan2((t248 * t344 + t245 * t289 / 0.2e1) * t334, (t245 * t344 - t289 * t248 / 0.2e1) * t334);
t233 = atan2((t246 * t347 + t280 * t348) * t333, (t247 * t347 + t280 * t349) * t333);
t232 = cos(t234);
t231 = sin(t234);
t230 = cos(t233);
t229 = sin(t233);
t225 = -pkin(4) * t227 + pkin(5);
t224 = -pkin(4) + t342;
t222 = t279 * t236 - t283 * t238;
t221 = -t283 * t236 - t279 * t238;
t220 = t226 + t298 + t316;
t219 = -t298 + t300 + t321;
t217 = atan2(t218 * t332 / 0.2e1, -(t298 - t316 + 0.2e1 * t324) * t332 / 0.2e1);
t216 = cos(t217);
t215 = sin(t217);
t214 = pkin(4) * t228 * t220 + t225 * t218;
t213 = -pkin(4) * t339 + t225 * t220;
t212 = -t278 * t215 + t282 * t216;
t211 = t282 * t215 + t278 * t216;
t210 = atan2((pkin(5) * t228 * t219 - t224 * t218) * t315, (-pkin(5) * t339 - t224 * t219) * t315);
t209 = cos(t210);
t208 = sin(t210);
t207 = atan2((t213 * t346 + t214 * t282 / 0.2e1) * t338, (-t213 * t282 / 0.2e1 + t214 * t346) * t338);
t206 = cos(t207);
t205 = sin(t207);
t204 = -t277 * t208 - t281 * t209;
t203 = -t281 * t208 + t277 * t209;
t1 = [t293, -t287, 0, 0; t287, t293, 0, 0; 0, 0, 1, pkin(14); 0, 0, 0, 1; -t286, -t292, 0, pkin(16); 0, 0, -1, 0; t292, -t286, 0, 0; 0, 0, 0, 1; t285, t291, 0, pkin(1); -t291, t285, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t206, -t205, 0, pkin(5); t205, t206, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t290, -t284, 0, pkin(10); 0, 0, -1, -pkin(12); t284, t290, 0, 0; 0, 0, 0, 1; t232, -t231, 0, -pkin(15); 0, 0, -1, 0; t231, t232, 0, -pkin(17); 0, 0, 0, 1; t230, -t229, 0, pkin(1); t229, t230, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t237, -t235, 0, 0; t235, t237, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t250, t249, 0, pkin(2); -t249, -t250, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t204, -t203, 0, t281 * pkin(4); t203, t204, 0, -t277 * pkin(4); 0, 0, 1, 0; 0, 0, 0, 1; t244, -t243, 0, t280 * pkin(3); t243, t244, 0, t276 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; t222, -t221, 0, t283 * pkin(6); t221, t222, 0, t279 * pkin(6); 0, 0, 1, 0; 0, 0, 0, 1; t212, -t211, 0, t282 * pkin(11); t211, t212, 0, t278 * pkin(11); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(8); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(13); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(9); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,16);             % numerisch
else,                         T_mdh = sym('xx', [4,4,16]); end % symbolisch

for i = 1:16
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
