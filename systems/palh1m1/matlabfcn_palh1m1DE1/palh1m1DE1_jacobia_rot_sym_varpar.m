% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh1m1DE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh1m1DE1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_jacobia_rot_sym_varpar: pkin has to be [23x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:44:06
	% EndTime: 2020-04-14 18:44:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:45:13
	% EndTime: 2020-04-14 18:47:43
	% DurationCPUTime: 87.56s
	% Computational Cost: add. (2346557->115), mult. (3533818->242), div. (161853->20), fcn. (2237868->24), ass. (0->151)
	t324 = pkin(7) ^ 2;
	t316 = sin(qJ(2));
	t320 = cos(qJ(2));
	t322 = cos(pkin(19));
	t387 = sin(pkin(19));
	t301 = t316 * t322 - t320 * t387;
	t384 = pkin(7) * t301;
	t396 = -2 * pkin(1);
	t359 = (pkin(1) ^ 2) + t384 * t396;
	t293 = t324 + t359;
	t290 = pkin(3) ^ 2 - pkin(8) ^ 2 + t293;
	t298 = pkin(1) - t384;
	t303 = t316 * t387 + t320 * t322;
	t393 = -pkin(8) - pkin(3);
	t288 = (pkin(7) - t393) * (pkin(7) + t393) + t359;
	t392 = -pkin(8) + pkin(3);
	t289 = (pkin(7) - t392) * (pkin(7) + t392) + t359;
	t330 = sqrt(-t289 * t288);
	t383 = pkin(7) * t303;
	t358 = pkin(1) * t383;
	t370 = 0.2e1 / t330 * (t288 + t289) * t358;
	t276 = (t301 * t330 + (-t370 / 0.2e1 - t290 + t298 * t396) * t303) * pkin(7);
	t368 = t303 * t330;
	t277 = t298 * t370 / 0.2e1 + t324 * t303 ^ 2 * t396 + (-t301 * t290 - t368) * pkin(7);
	t291 = 0.1e1 / t293;
	t315 = sin(qJ(3));
	t327 = 0.1e1 / pkin(3);
	t348 = 0.1e1 / t293 ^ 2 * t358;
	t284 = t290 * t383 + t298 * t330;
	t319 = cos(qJ(3));
	t363 = t319 * t284;
	t283 = -pkin(7) * t368 + t290 * t298;
	t367 = t315 * t283;
	t388 = -t319 / 0.2e1;
	t262 = ((-t315 * t276 / 0.2e1 + t277 * t388) * t291 + (-t363 - t367) * t348) * t327;
	t364 = t319 * t283;
	t366 = t315 * t284;
	t263 = ((t276 * t388 + t315 * t277 / 0.2e1) * t291 + (-t364 + t366) * t348) * t327;
	t311 = pkin(23) + pkin(22);
	t309 = sin(t311);
	t310 = cos(t311);
	t258 = t262 * t310 + t263 * t309;
	t369 = t291 * t327;
	t279 = (t367 / 0.2e1 + t363 / 0.2e1) * t369;
	t280 = (-t364 / 0.2e1 + t366 / 0.2e1) * t369;
	t272 = t279 * t310 - t280 * t309;
	t386 = pkin(4) * t272;
	t395 = -2 * pkin(5);
	t360 = (pkin(5) ^ 2) - t386 * t395;
	t391 = -pkin(9) - pkin(11);
	t264 = (pkin(4) - t391) * (pkin(4) + t391) + t360;
	t390 = -pkin(9) + pkin(11);
	t265 = (pkin(4) - t390) * (pkin(4) + t390) + t360;
	t394 = pkin(4) * pkin(5);
	t332 = 0.2e1 * (t264 + t265) * t394;
	t254 = t258 * t332;
	t259 = -t262 * t309 + t263 * t310;
	t329 = sqrt(-t265 * t264);
	t326 = pkin(4) ^ 2;
	t269 = t326 + t360;
	t266 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t269;
	t270 = pkin(5) + t386;
	t349 = t270 * t395 - t266;
	t261 = 0.1e1 / t329;
	t344 = t279 * t309 + t310 * t280;
	t355 = -t261 * t344 / 0.2e1;
	t236 = (t254 * t355 + t349 * t258 - t259 * t329) * pkin(4);
	t350 = t326 * t344 * t395;
	t356 = t261 * t270 / 0.2e1;
	t237 = t254 * t356 + t258 * t350 + (-t258 * t329 + t259 * t266) * pkin(4);
	t385 = pkin(4) * t344;
	t255 = t266 * t270 - t329 * t385;
	t256 = t266 * t385 + t270 * t329;
	t313 = cos(pkin(21));
	t267 = 0.1e1 / t269;
	t323 = 0.1e1 / pkin(11);
	t371 = t267 * t323;
	t312 = sin(pkin(21));
	t389 = t312 / 0.2e1;
	t253 = (-t255 * t313 / 0.2e1 + t256 * t389) * t371;
	t250 = 0.1e1 / t253;
	t357 = 0.1e1 / t269 ^ 2 * t394;
	t347 = t258 * t357;
	t342 = t256 * t347;
	t343 = t255 * t347;
	t372 = t267 * t313;
	t352 = t372 / 0.2e1;
	t353 = -t372 / 0.2e1;
	t354 = t267 * t389;
	t251 = 0.1e1 / t253 ^ 2;
	t252 = (t255 * t389 + t256 * t313 / 0.2e1) * t371;
	t373 = t251 * t252;
	t374 = 0.1e1 / (t251 * t252 ^ 2 + 0.1e1) * t323;
	t403 = ((t236 * t354 + t237 * t352 + t312 * t343 + t313 * t342) * t250 - (t236 * t353 + t237 * t354 + t312 * t342 - t313 * t343) * t373) * t374 + 0.1e1;
	t257 = t344 * t332;
	t248 = (t257 * t355 - t272 * t329 + t344 * t349) * pkin(4);
	t249 = t257 * t356 + t344 * t350 + (t272 * t266 - t329 * t344) * pkin(4);
	t346 = t344 * t357;
	t340 = t313 * t346;
	t341 = t312 * t346;
	t402 = ((t248 * t354 + t249 * t352 + t255 * t341 + t256 * t340) * t250 - (t248 * t353 + t249 * t354 - t255 * t340 + t256 * t341) * t373) * t374 + 0.1e1;
	t247 = atan2(t252, t253);
	t244 = sin(t247);
	t245 = cos(t247);
	t302 = -t315 * t316 + t319 * t320;
	t333 = t320 * t315 + t316 * t319;
	t235 = t244 * t333 - t245 * t302;
	t317 = sin(qJ(1));
	t294 = t302 * t317;
	t295 = t333 * t317;
	t337 = -t244 * t294 - t295 * t245;
	t223 = atan2(t337, t235);
	t220 = sin(t223);
	t221 = cos(t223);
	t216 = t220 * t337 + t221 * t235;
	t214 = 0.1e1 / t216;
	t321 = cos(qJ(1));
	t296 = t302 * t321;
	t297 = t333 * t321;
	t399 = -t297 * t244 + t245 * t296;
	t401 = t214 * t399;
	t314 = sin(qJ(4));
	t318 = cos(qJ(4));
	t227 = t314 * t317 + t318 * t399;
	t225 = 0.1e1 / t227 ^ 2;
	t226 = t314 * t399 - t317 * t318;
	t351 = t225 * t226 ^ 2 + 0.1e1;
	t219 = 0.1e1 / t351;
	t224 = 0.1e1 / t227;
	t379 = t225 * t226;
	t398 = -t296 * t244 - t297 * t245;
	t400 = (t314 * t224 - t318 * t379) * t219 * t398;
	t228 = t244 * t295 - t294 * t245;
	t397 = t302 * t244 + t245 * t333;
	t215 = 0.1e1 / t216 ^ 2;
	t382 = t215 * t398;
	t381 = t215 * t398 ^ 2;
	t380 = t221 * t337;
	t234 = 0.1e1 / t235 ^ 2;
	t378 = t337 * t234;
	t339 = -t220 * t235 + t380;
	t233 = 0.1e1 / t235;
	t222 = 0.1e1 / (t234 * t337 ^ 2 + 0.1e1);
	t213 = 0.1e1 / (0.1e1 + t381);
	t212 = t402 * t397;
	t210 = t402 * t228;
	t209 = t403 * t397;
	t207 = t403 * t228;
	t206 = (t210 * t233 - t212 * t378) * t222;
	t205 = (t207 * t233 - t209 * t378) * t222;
	t1 = [t398 * t233 * t222, t205, t206, 0; (t337 * t214 - (-t220 + (-t233 * t380 + t220) * t222) * t381) * t213, ((t339 * t205 + t207 * t220 + t209 * t221) * t382 + t403 * t401) * t213, ((t339 * t206 + t210 * t220 + t212 * t221) * t382 + t402 * t401) * t213, 0; ((t228 * t314 - t318 * t321) * t224 - (t228 * t318 + t314 * t321) * t379) * t219, t403 * t400, t402 * t400, t351 * t219;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:59
	% EndTime: 2020-04-14 18:42:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:01
	% EndTime: 2020-04-14 18:43:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:59
	% EndTime: 2020-04-14 18:42:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:02
	% EndTime: 2020-04-14 18:43:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:44:14
	% EndTime: 2020-04-14 18:44:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end