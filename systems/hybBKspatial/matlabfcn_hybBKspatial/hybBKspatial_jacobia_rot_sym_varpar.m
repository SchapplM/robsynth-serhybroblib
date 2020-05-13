% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% hybBKspatial
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in hybBKspatial_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED,L1,L2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 19:31
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = hybBKspatial_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'hybBKspatial_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'hybBKspatial_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'hybBKspatial_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:49
	% EndTime: 2020-04-11 19:22:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:49
	% EndTime: 2020-04-11 19:22:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:49
	% EndTime: 2020-04-11 19:22:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:23:20
	% EndTime: 2020-04-11 19:23:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:49
	% EndTime: 2020-04-11 19:22:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:51
	% EndTime: 2020-04-11 19:22:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:52
	% EndTime: 2020-04-11 19:22:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:53
	% EndTime: 2020-04-11 19:23:04
	% DurationCPUTime: 5.59s
	% Computational Cost: add. (99532->167), mult. (213928->331), div. (2707->13), fcn. (55037->15), ass. (0->161)
	t257 = sin(qJ(3));
	t258 = sin(qJ(2));
	t317 = t257 * t258;
	t247 = pkin(2) * t317;
	t292 = pkin(3) * t247;
	t246 = -0.2e1 * t292;
	t269 = pkin(3) ^ 2;
	t270 = pkin(2) ^ 2;
	t263 = cos(qJ(2));
	t334 = pkin(2) * t263;
	t248 = pkin(1) - t334;
	t262 = cos(qJ(3));
	t242 = t248 * t262;
	t302 = pkin(3) * t242;
	t306 = pkin(1) * t334;
	t250 = -0.2e1 * t306;
	t271 = pkin(1) ^ 2;
	t309 = t250 + t271;
	t231 = t246 + t269 + t270 + 0.2e1 * t302 + t309;
	t230 = 0.1e1 / t231 ^ 2;
	t349 = -0.2e1 * t230;
	t267 = pkin(4) ^ 2;
	t308 = -pkin(5) ^ 2 + t271;
	t299 = t270 + t308;
	t290 = t269 + t299;
	t238 = t250 + t267 + t290;
	t254 = t263 ^ 2;
	t343 = 0.2e1 * t254;
	t240 = t270 * t343 - t270 + t309;
	t320 = t240 * t262;
	t303 = pkin(3) * t320;
	t227 = t248 * t238 + 0.2e1 * t303;
	t253 = t262 ^ 2;
	t333 = pkin(3) * t248;
	t228 = t238 * t262 + (0.4e1 * t253 - 0.2e1) * t333;
	t310 = -t247 + t242;
	t236 = pkin(3) + t310;
	t244 = -t267 + t290;
	t239 = t250 + t244;
	t232 = t246 + t239;
	t251 = (t269 - t271) * t270;
	t345 = -0.2e1 * t240;
	t288 = t253 * t345 - t308;
	t341 = -pkin(4) + pkin(5);
	t342 = -pkin(4) - pkin(5);
	t272 = sqrt(0.4e1 * t251 * t254 + 0.4e1 * t244 * t306 - (t271 + (pkin(2) - t342) * (pkin(2) + t342)) * (t271 + (pkin(2) - t341) * (pkin(2) + t341)) + 0.4e1 * (-t232 * t242 + t239 * t247) * pkin(3) + (0.2e1 * t267 - 0.6e1 * t270 + 0.2e1 * t288 - t269) * t269);
	t335 = pkin(2) * t258;
	t221 = t227 * t257 + t228 * t335 + t236 * t272;
	t243 = 0.3e1 * t269 + t267 + t299;
	t233 = t243 + t250 - 0.4e1 * t292;
	t315 = t258 * t262;
	t237 = pkin(2) * t315 + t257 * t248;
	t321 = t237 * t272;
	t336 = pkin(1) * (t247 - pkin(3));
	t220 = -t233 * t242 + t321 + (t243 * t317 - 0.2e1 * t263 * t336) * pkin(2) + (-t267 - t269 + (t343 - 0.3e1) * t270 + t288) * pkin(3);
	t318 = t257 * t220;
	t337 = t262 / 0.2e1;
	t282 = t221 * t337 + t318 / 0.2e1;
	t229 = 0.1e1 / t231;
	t268 = 0.1e1 / pkin(4);
	t322 = t229 * t268;
	t348 = t282 * t322;
	t339 = t257 / 0.2e1;
	t297 = t221 * t339;
	t338 = -t262 / 0.2e1;
	t283 = t220 * t338 + t297;
	t264 = cos(qJ(1));
	t312 = t264 * t268;
	t207 = t283 * t229 * t312;
	t256 = sin(qJ(4));
	t261 = cos(qJ(4));
	t276 = t264 * t348;
	t199 = t207 * t261 + t256 * t276;
	t255 = sin(qJ(5));
	t259 = sin(qJ(1));
	t260 = cos(qJ(5));
	t193 = t199 * t260 - t259 * t255;
	t191 = 0.1e1 / t193 ^ 2;
	t192 = t199 * t255 + t259 * t260;
	t294 = t191 * t192 ^ 2 + 0.1e1;
	t181 = 0.1e1 / t294;
	t190 = 0.1e1 / t193;
	t327 = t191 * t192;
	t281 = (-t255 * t190 + t260 * t327) * t181;
	t347 = -0.4e1 * pkin(2);
	t346 = 0.2e1 * t230;
	t344 = 0.2e1 * t248;
	t224 = 0.1e1 / t272;
	t319 = t253 * t269;
	t331 = pkin(3) * t262;
	t305 = 0.2e1 * t224 * ((t232 * t333 + 0.2e1 * t269 * t320) * t257 + (t239 * t331 + t319 * t344) * t335);
	t307 = -0.2e1 * pkin(1) * t270;
	t313 = t262 * t263;
	t332 = pkin(3) * t253;
	t213 = t310 * t272 + t237 * t305 + (t248 * t233 + 0.4e1 * t303) * t257 + (t307 * t313 + (t243 * t262 + 0.4e1 * t248 * t332) * pkin(2)) * t258;
	t340 = -t213 / 0.2e1;
	t277 = t259 * t348;
	t314 = t259 * t268;
	t300 = t229 * t314;
	t206 = t283 * t300;
	t325 = t206 * t256;
	t195 = -t261 * t277 + t325;
	t289 = t322 * t338;
	t210 = t220 * t289 + t297 * t322;
	t202 = t210 * t261 + t256 * t348;
	t185 = atan2(t195, t202);
	t182 = sin(t185);
	t183 = cos(t185);
	t174 = t182 * t195 + t183 * t202;
	t173 = 0.1e1 / t174 ^ 2;
	t197 = t207 * t256 - t261 * t276;
	t330 = t173 * t197;
	t329 = t173 * t197 ^ 2;
	t328 = t183 * t195;
	t201 = 0.1e1 / t202 ^ 2;
	t326 = t195 * t201;
	t324 = t206 * t261;
	t323 = t210 * t256;
	t316 = t257 * t263;
	t311 = t259 * t221 * t289 - t300 * t318 / 0.2e1;
	t304 = -0.4e1 * t263 * t270;
	t301 = pkin(2) * t316;
	t252 = pkin(1) * t335;
	t249 = 0.2e1 * t252;
	t241 = t258 * t304 + t249;
	t291 = pkin(3) * t301;
	t293 = -0.8e1 * t302;
	t298 = (-0.4e1 * t241 * t319 + (t252 - t291) * t293 + 0.4e1 * t239 * t291 + (t232 * t331 * t347 - 0.8e1 * t251 * t263 + (0.8e1 * pkin(3) * t270 * t317 + t244 * t347) * pkin(1)) * t258) * t224 / 0.2e1;
	t212 = -t321 + t236 * t305 + t227 * t262 + ((-t238 + t293) * t335 + pkin(3) * t257 * t345) * t257;
	t296 = t212 / 0.2e1 + t220 / 0.2e1;
	t295 = t221 / 0.2e1 + t340;
	t287 = -t182 * t202 + t328;
	t208 = t236 * t298 + 0.2e1 * t241 * t257 * t331 + ((-t257 * t272 + t228) * t263 + (t262 * t272 + (pkin(1) * t344 + t238) * t257 + 0.2e1 * (pkin(1) * t262 - pkin(3) + 0.2e1 * t332) * t335) * t258) * pkin(2);
	t211 = (pkin(2) * t313 + t247) * t272 + t237 * t298 - 0.2e1 * t241 * t332 - (t249 - 0.4e1 * t291) * t242 + t254 * t257 * t307 + t243 * t301 + (pkin(3) * t304 + (-t233 * t262 + 0.2e1 * t336) * pkin(2)) * t258;
	t285 = t208 * t339 + t211 * t338;
	t284 = t208 * t337 + t211 * t339;
	t280 = t283 * t349;
	t279 = t282 * t349;
	t275 = (-t283 * t256 + t282 * t261) * t346;
	t274 = (-t282 * t256 - t283 * t261) * t346;
	t235 = t237 * pkin(3);
	t234 = t252 + (t315 - t316) * pkin(3) * pkin(2);
	t203 = t261 * t348 - t323;
	t200 = 0.1e1 / t202;
	t196 = t256 * t277 + t324;
	t194 = t311 * t256 - t324;
	t189 = (t284 * t229 + t234 * t279) * t312;
	t188 = (t285 * t229 + t234 * t280) * t312;
	t187 = (-t235 * t279 + (-t295 * t257 + t296 * t262) * t229) * t312;
	t186 = (-t235 * t280 + (t296 * t257 + t295 * t262) * t229) * t312;
	t184 = 0.1e1 / (t195 ^ 2 * t201 + 0.1e1);
	t180 = ((t284 * t256 + t285 * t261) * t229 + t234 * t274) * t268;
	t178 = ((t285 * t256 - t284 * t261) * t229 + t234 * t275) * t314;
	t177 = -t323 + (-t235 * t274 + (-(t212 * t338 + t257 * t340) * t256 + (t212 * t339 + t213 * t338 + t282) * t261) * t229) * t268;
	t175 = (-t235 * t275 + ((t295 * t256 - t296 * t261) * t262 + (t296 * t256 + t295 * t261) * t257) * t229) * t314;
	t172 = 0.1e1 / t174;
	t171 = 0.1e1 / (0.1e1 + t329);
	t170 = (t196 * t200 - t203 * t326) * t184;
	t169 = (t178 * t200 - t180 * t326) * t184;
	t168 = (t175 * t200 - t177 * t326) * t184;
	t1 = [t197 * t200 * t184, t169, t168, t170, 0, 0; ((t311 * t261 + t325) * t172 + (t182 + (t200 * t328 - t182) * t184) * t329) * t171, ((-t188 * t256 + t189 * t261) * t172 + (t287 * t169 + t178 * t182 + t180 * t183) * t330) * t171, ((-t186 * t256 + t187 * t261) * t172 + (t287 * t168 + t175 * t182 + t177 * t183) * t330) * t171, (-t199 * t172 + (t287 * t170 + t182 * t196 + t183 * t203) * t330) * t171, 0, 0; ((t194 * t255 + t264 * t260) * t190 - (t194 * t260 - t264 * t255) * t327) * t181, -(t188 * t261 + t189 * t256) * t281, -(t186 * t261 + t187 * t256) * t281, t197 * t281, t294 * t181, 0;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-11 19:22:55
	% EndTime: 2020-04-11 19:23:09
	% DurationCPUTime: 7.63s
	% Computational Cost: add. (138152->185), mult. (296370->380), div. (3813->15), fcn. (76792->17), ass. (0->174)
	t300 = cos(qJ(2));
	t372 = pkin(2) * t300;
	t341 = pkin(1) * t372;
	t283 = -0.2e1 * t341;
	t304 = pkin(4) ^ 2;
	t306 = pkin(3) ^ 2;
	t307 = pkin(2) ^ 2;
	t308 = pkin(1) ^ 2;
	t343 = -pkin(5) ^ 2 + t308;
	t332 = t307 + t343;
	t322 = t306 + t332;
	t270 = t283 + t304 + t322;
	t281 = pkin(1) - t372;
	t344 = t283 + t308;
	t289 = t300 ^ 2;
	t380 = 0.2e1 * t289;
	t272 = t307 * t380 - t307 + t344;
	t299 = cos(qJ(3));
	t358 = t272 * t299;
	t337 = pkin(3) * t358;
	t259 = t270 * t281 + 0.2e1 * t337;
	t288 = t299 ^ 2;
	t371 = pkin(3) * t281;
	t260 = t270 * t299 + (0.4e1 * t288 - 0.2e1) * t371;
	t274 = t281 * t299;
	t293 = sin(qJ(3));
	t294 = sin(qJ(2));
	t353 = t293 * t294;
	t279 = pkin(2) * t353;
	t345 = -t279 + t274;
	t268 = pkin(3) + t345;
	t276 = -t304 + t322;
	t271 = t283 + t276;
	t325 = pkin(3) * t279;
	t278 = -0.2e1 * t325;
	t264 = t278 + t271;
	t284 = (t306 - t308) * t307;
	t382 = -0.2e1 * t272;
	t320 = t288 * t382 - t343;
	t378 = -pkin(4) + pkin(5);
	t379 = -pkin(4) - pkin(5);
	t309 = sqrt(0.4e1 * t284 * t289 + 0.4e1 * t276 * t341 - (t308 + (pkin(2) - t379) * (pkin(2) + t379)) * (t308 + (pkin(2) - t378) * (pkin(2) + t378)) + 0.4e1 * (-t264 * t274 + t271 * t279) * pkin(3) + (0.2e1 * t304 - 0.6e1 * t307 + 0.2e1 * t320 - t306) * t306);
	t373 = pkin(2) * t294;
	t253 = t259 * t293 + t260 * t373 + t268 * t309;
	t275 = 0.3e1 * t306 + t304 + t332;
	t265 = t275 + t283 - 0.4e1 * t325;
	t351 = t294 * t299;
	t269 = pkin(2) * t351 + t281 * t293;
	t359 = t269 * t309;
	t374 = pkin(1) * (t279 - pkin(3));
	t252 = -t265 * t274 + t359 + (t275 * t353 - 0.2e1 * t300 * t374) * pkin(2) + (-t304 - t306 + (t380 - 0.3e1) * t307 + t320) * pkin(3);
	t375 = -t299 / 0.2e1;
	t330 = t252 * t375;
	t376 = t293 / 0.2e1;
	t316 = t253 * t376 + t330;
	t336 = pkin(3) * t274;
	t263 = t278 + t306 + t307 + 0.2e1 * t336 + t344;
	t261 = 0.1e1 / t263;
	t295 = sin(qJ(1));
	t305 = 0.1e1 / pkin(4);
	t350 = t295 * t305;
	t334 = t261 * t350;
	t234 = t316 * t334;
	t348 = t299 * t253;
	t377 = -t293 / 0.2e1;
	t315 = -t348 / 0.2e1 + t252 * t377;
	t235 = t315 * t334;
	t292 = sin(qJ(4));
	t298 = cos(qJ(4));
	t225 = -t234 * t298 + t235 * t292;
	t291 = sin(qJ(5));
	t297 = cos(qJ(5));
	t301 = cos(qJ(1));
	t385 = -t225 * t291 - t297 * t301;
	t220 = -t225 * t297 + t291 * t301;
	t384 = -0.4e1 * pkin(2);
	t262 = 0.1e1 / t263 ^ 2;
	t383 = 0.2e1 * t262;
	t381 = 0.2e1 * t281;
	t370 = pkin(3) * t288;
	t369 = pkin(3) * t299;
	t360 = t261 * t305;
	t239 = t315 * t360;
	t321 = t360 * t376;
	t240 = t253 * t321 + t330 * t360;
	t233 = t239 * t298 + t240 * t292;
	t362 = t233 * t297;
	t212 = atan2(t220, -t362);
	t210 = sin(t212);
	t211 = cos(t212);
	t195 = t210 * t220 - t211 * t362;
	t194 = 0.1e1 / t195 ^ 2;
	t346 = t301 * t305;
	t333 = t261 * t346;
	t236 = t316 * t333;
	t237 = t301 * t252 * t321 + t333 * t348 / 0.2e1;
	t229 = t236 * t298 + t237 * t292;
	t222 = t229 * t297 - t295 * t291;
	t368 = t194 * t222;
	t367 = t194 * t222 ^ 2;
	t221 = t229 * t291 + t295 * t297;
	t228 = -t236 * t292 + t237 * t298;
	t290 = sin(qJ(6));
	t296 = cos(qJ(6));
	t209 = t221 * t296 + t228 * t290;
	t206 = 0.1e1 / t209 ^ 2;
	t207 = -t221 * t290 + t228 * t296;
	t366 = t206 * t207;
	t365 = t210 * t233;
	t364 = t211 * t220;
	t230 = 0.1e1 / t233;
	t286 = 0.1e1 / t297;
	t363 = t230 * t286;
	t361 = t240 * t298;
	t357 = t288 * t306;
	t356 = t290 * t291;
	t355 = t291 * t296;
	t352 = t293 * t300;
	t347 = t299 * t300;
	t342 = -0.2e1 * pkin(1) * t307;
	t256 = 0.1e1 / t309;
	t340 = 0.2e1 * t256 * ((t264 * t371 + 0.2e1 * t306 * t358) * t293 + (t271 * t369 + t357 * t381) * t373);
	t339 = -0.4e1 * t300 * t307;
	t338 = pkin(2) * t352;
	t231 = 0.1e1 / t233 ^ 2;
	t335 = t220 * t231 * t286;
	t285 = pkin(1) * t373;
	t282 = 0.2e1 * t285;
	t273 = t294 * t339 + t282;
	t324 = pkin(3) * t338;
	t326 = -0.8e1 * t336;
	t331 = (-0.4e1 * t273 * t357 + (t285 - t324) * t326 + 0.4e1 * t271 * t324 + (t264 * t369 * t384 - 0.8e1 * t284 * t300 + (0.8e1 * pkin(3) * t307 * t353 + t276 * t384) * pkin(1)) * t294) * t256 / 0.2e1;
	t242 = -t359 + t268 * t340 + t259 * t299 + ((-t270 + t326) * t373 + pkin(3) * t293 * t382) * t293;
	t329 = t242 / 0.2e1 + t252 / 0.2e1;
	t243 = t345 * t309 + t269 * t340 + (t281 * t265 + 0.4e1 * t337) * t293 + (t342 * t347 + (t275 * t299 + 0.4e1 * t281 * t370) * pkin(2)) * t294;
	t328 = t253 / 0.2e1 - t243 / 0.2e1;
	t327 = t206 * t207 ^ 2 + 0.1e1;
	t323 = t364 * t368;
	t226 = -t234 * t292 - t235 * t298;
	t238 = t268 * t331 + 0.2e1 * t273 * t293 * t369 + ((-t293 * t309 + t260) * t300 + (t299 * t309 + (pkin(1) * t381 + t270) * t293 + 0.2e1 * (pkin(1) * t299 - pkin(3) + 0.2e1 * t370) * t373) * t294) * pkin(2);
	t241 = (pkin(2) * t347 + t279) * t309 + t269 * t331 - 0.2e1 * t273 * t370 - (t282 - 0.4e1 * t324) * t274 + t289 * t293 * t342 + t275 * t338 + (pkin(3) * t339 + (-t265 * t299 + 0.2e1 * t374) * pkin(2)) * t294;
	t319 = t238 * t376 + t241 * t375;
	t318 = t299 * t238 / 0.2e1 + t241 * t376;
	t317 = t242 * t376 + t243 * t375;
	t314 = -0.2e1 * t262 * t316;
	t313 = t315 * t383;
	t312 = (-t316 * t292 - t315 * t298) * t383;
	t311 = (t315 * t292 - t316 * t298) * t383;
	t287 = 0.1e1 / t297 ^ 2;
	t267 = t269 * pkin(3);
	t266 = t285 + (t351 - t352) * pkin(3) * pkin(2);
	t232 = -t239 * t292 + t361;
	t217 = (t318 * t261 + t266 * t313) * t346;
	t216 = (t319 * t261 + t266 * t314) * t346;
	t215 = (-t267 * t313 + (-t328 * t293 + t329 * t299) * t261) * t346;
	t214 = (t317 * t261 - t267 * t314) * t346 + t237;
	t213 = 0.1e1 / (t220 ^ 2 * t231 * t287 + 0.1e1);
	t205 = 0.1e1 / t209;
	t204 = ((t319 * t292 - t318 * t298) * t261 + t266 * t312) * t305;
	t203 = t216 * t298 + t217 * t292;
	t202 = -t216 * t292 + t217 * t298;
	t201 = ((t318 * t292 + t319 * t298) * t261 + t266 * t311) * t350;
	t200 = t361 + (-t267 * t312 + ((t242 * t375 + t243 * t377) * t298 + (-t315 + t317) * t292) * t261) * t305;
	t199 = t214 * t298 + t215 * t292;
	t198 = -t214 * t292 + t215 * t298;
	t197 = (-t267 * t311 + ((t329 * t292 + t328 * t298) * t299 + (-t328 * t292 + t329 * t298) * t293) * t261) * t350;
	t196 = 0.1e1 / t327;
	t193 = 0.1e1 / t195;
	t192 = 0.1e1 / (0.1e1 + t367);
	t191 = (-t220 * t287 * t291 + t286 * t385) * t230 * t213;
	t190 = (-t226 * t230 + t232 * t335) * t213;
	t189 = (-t201 * t230 + t204 * t335) * t213;
	t188 = (-t197 * t230 + t200 * t335) * t213;
	t1 = [-t222 * t213 * t363, t189, t188, t190, t191, 0; (t220 * t193 + (t210 + (-t363 * t364 - t210) * t213) * t367) * t192, (t189 * t323 + (-t203 * t193 + (t189 * t365 + t201 * t210 - t204 * t211) * t368) * t297) * t192, (t188 * t323 + (-t199 * t193 + (t188 * t365 + t197 * t210 - t200 * t211) * t368) * t297) * t192, (t190 * t323 + (-t228 * t193 + (t190 * t365 + t210 * t226 - t211 * t232) * t368) * t297) * t192, (t221 * t193 + (t211 * t233 * t291 - t210 * t385 + (t210 * t362 + t364) * t191) * t368) * t192, 0; ((t226 * t296 - t290 * t385) * t205 + (-t226 * t290 - t296 * t385) * t366) * t196, ((-t202 * t296 + t203 * t356) * t205 + (t202 * t290 + t203 * t355) * t366) * t196, ((-t198 * t296 + t199 * t356) * t205 + (t198 * t290 + t199 * t355) * t366) * t196, ((t228 * t356 + t229 * t296) * t205 + (t228 * t355 - t229 * t290) * t366) * t196, (t290 * t205 + t296 * t366) * t222 * t196, t327 * t196;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end