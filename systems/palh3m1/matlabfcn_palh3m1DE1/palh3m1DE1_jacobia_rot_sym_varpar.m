% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m1DE1
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
%   Wie in palh3m1DE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh3m1DE1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_jacobia_rot_sym_varpar: pkin has to be [19x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:19:02
	% EndTime: 2020-04-19 18:19:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:20:10
	% EndTime: 2020-04-19 18:22:34
	% DurationCPUTime: 88.81s
	% Computational Cost: add. (2346557->117), mult. (3533818->246), div. (161853->20), fcn. (2237868->24), ass. (0->153)
	t321 = pkin(5) ^ 2;
	t316 = sin(pkin(16));
	t319 = cos(qJ(2));
	t387 = sin(qJ(2));
	t389 = cos(pkin(16));
	t300 = t387 * t316 - t319 * t389;
	t384 = pkin(5) * t300;
	t398 = -2 * pkin(1);
	t360 = (pkin(1) ^ 2) + t384 * t398;
	t292 = t321 + t360;
	t289 = pkin(2) ^ 2 - pkin(6) ^ 2 + t292;
	t297 = pkin(1) - t384;
	t302 = t319 * t316 + t387 * t389;
	t395 = -pkin(6) - pkin(2);
	t287 = (pkin(5) - t395) * (pkin(5) + t395) + t360;
	t394 = -pkin(6) + pkin(2);
	t288 = (pkin(5) - t394) * (pkin(5) + t394) + t360;
	t327 = sqrt(-t288 * t287);
	t383 = pkin(5) * t302;
	t359 = pkin(1) * t383;
	t369 = 0.2e1 / t327 * (t287 + t288) * t359;
	t276 = (t300 * t327 + (-t369 / 0.2e1 - t289 + t297 * t398) * t302) * pkin(5);
	t364 = t302 * t327;
	t277 = t297 * t369 / 0.2e1 + t321 * t302 ^ 2 * t398 + (-t300 * t289 - t364) * pkin(5);
	t290 = 0.1e1 / t292;
	t318 = cos(qJ(3));
	t324 = 0.1e1 / pkin(2);
	t347 = 0.1e1 / t292 ^ 2 * t359;
	t283 = t289 * t383 + t297 * t327;
	t370 = t283 * t318;
	t282 = -pkin(5) * t364 + t297 * t289;
	t314 = sin(qJ(3));
	t373 = t282 * t314;
	t390 = t314 / 0.2e1;
	t262 = ((t277 * t318 / 0.2e1 + t276 * t390) * t290 + (t370 + t373) * t347) * t324;
	t371 = t283 * t314;
	t372 = t282 * t318;
	t263 = ((-t276 * t318 / 0.2e1 + t277 * t390) * t290 + (t371 - t372) * t347) * t324;
	t310 = pkin(18) + pkin(19);
	t308 = sin(t310);
	t309 = cos(t310);
	t259 = -t308 * t262 - t309 * t263;
	t368 = t290 * t324;
	t278 = (-t372 / 0.2e1 + t371 / 0.2e1) * t368;
	t279 = (t370 / 0.2e1 + t373 / 0.2e1) * t368;
	t273 = t309 * t278 + t308 * t279;
	t385 = pkin(3) * t273;
	t397 = -2 * pkin(4);
	t361 = (pkin(4) ^ 2) - t385 * t397;
	t393 = -pkin(8) - pkin(10);
	t264 = (pkin(3) - t393) * (pkin(3) + t393) + t361;
	t392 = -pkin(8) + pkin(10);
	t265 = (pkin(3) - t392) * (pkin(3) + t392) + t361;
	t396 = pkin(3) * pkin(4);
	t330 = 0.2e1 * (t264 + t265) * t396;
	t254 = t259 * t330;
	t258 = -t309 * t262 + t308 * t263;
	t326 = sqrt(-t265 * t264);
	t323 = pkin(3) ^ 2;
	t269 = t323 + t361;
	t266 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t269;
	t270 = pkin(4) + t385;
	t348 = t270 * t397 - t266;
	t261 = 0.1e1 / t326;
	t331 = t308 * t278 - t309 * t279;
	t354 = -t261 * t331 / 0.2e1;
	t236 = (t254 * t354 - t258 * t326 + t348 * t259) * pkin(3);
	t349 = t323 * t331 * t397;
	t355 = t261 * t270 / 0.2e1;
	t237 = t254 * t355 + t259 * t349 + (t258 * t266 - t259 * t326) * pkin(3);
	t386 = pkin(3) * t331;
	t255 = t270 * t266 - t326 * t386;
	t256 = t266 * t386 + t270 * t326;
	t312 = cos(pkin(17));
	t267 = 0.1e1 / t269;
	t320 = 0.1e1 / pkin(10);
	t374 = t267 * t320;
	t311 = sin(pkin(17));
	t391 = t311 / 0.2e1;
	t252 = (-t255 * t312 / 0.2e1 + t256 * t391) * t374;
	t250 = 0.1e1 / t252;
	t358 = 0.1e1 / t269 ^ 2 * t396;
	t345 = t312 * t358;
	t339 = t259 * t345;
	t346 = t311 * t358;
	t341 = t259 * t346;
	t375 = t267 * t312;
	t351 = t375 / 0.2e1;
	t352 = -t375 / 0.2e1;
	t353 = t267 * t391;
	t251 = 0.1e1 / t252 ^ 2;
	t253 = (t256 * t312 / 0.2e1 + t255 * t391) * t374;
	t376 = t251 * t253;
	t377 = 0.1e1 / (t253 ^ 2 * t251 + 0.1e1) * t320;
	t406 = ((t236 * t353 + t237 * t351 + t255 * t341 + t256 * t339) * t250 - (t236 * t352 + t237 * t353 - t255 * t339 + t256 * t341) * t376) * t377 + 0.1e1;
	t257 = t331 * t330;
	t248 = (t257 * t354 - t273 * t326 + t331 * t348) * pkin(3);
	t249 = t257 * t355 + t331 * t349 + (t273 * t266 - t326 * t331) * pkin(3);
	t338 = t331 * t345;
	t340 = t331 * t346;
	t405 = ((t248 * t353 + t249 * t351 + t255 * t340 + t256 * t338) * t250 - (t248 * t352 + t249 * t353 - t255 * t338 + t256 * t340) * t376) * t377 + 0.1e1;
	t247 = atan2(t253, t252);
	t244 = sin(t247);
	t245 = cos(t247);
	t299 = t387 * t314 - t319 * t318;
	t329 = t319 * t314 + t387 * t318;
	t235 = -t244 * t329 - t299 * t245;
	t315 = sin(qJ(1));
	t293 = t329 * t315;
	t294 = t299 * t315;
	t336 = -t294 * t244 + t293 * t245;
	t223 = atan2(t336, t235);
	t220 = sin(t223);
	t221 = cos(t223);
	t216 = t220 * t336 + t221 * t235;
	t214 = 0.1e1 / t216;
	t388 = cos(qJ(1));
	t342 = t388 * t387;
	t357 = t319 * t388;
	t295 = -t314 * t342 + t318 * t357;
	t296 = t314 * t357 + t318 * t342;
	t400 = t296 * t244 - t295 * t245;
	t404 = t214 * t400;
	t313 = sin(qJ(4));
	t317 = cos(qJ(4));
	t227 = t315 * t313 + t317 * t400;
	t225 = 0.1e1 / t227 ^ 2;
	t226 = t313 * t400 - t315 * t317;
	t350 = t225 * t226 ^ 2 + 0.1e1;
	t219 = 0.1e1 / t350;
	t224 = 0.1e1 / t227;
	t379 = t225 * t226;
	t401 = t295 * t244 + t296 * t245;
	t403 = (t313 * t224 - t317 * t379) * t219 * t401;
	t402 = t293 * t244 + t294 * t245;
	t399 = t299 * t244 - t245 * t329;
	t215 = 0.1e1 / t216 ^ 2;
	t382 = t215 * t401;
	t381 = t215 * t401 ^ 2;
	t380 = t221 * t336;
	t234 = 0.1e1 / t235 ^ 2;
	t378 = t336 * t234;
	t337 = -t220 * t235 + t380;
	t233 = 0.1e1 / t235;
	t222 = 0.1e1 / (t234 * t336 ^ 2 + 0.1e1);
	t213 = 0.1e1 / (0.1e1 + t381);
	t212 = t405 * t399;
	t210 = t405 * t402;
	t209 = t406 * t399;
	t207 = t406 * t402;
	t206 = (-t210 * t233 - t212 * t378) * t222;
	t205 = (-t207 * t233 - t209 * t378) * t222;
	t1 = [t401 * t233 * t222, t205, t206, 0; (t336 * t214 - (-t220 + (-t233 * t380 + t220) * t222) * t381) * t213, ((t337 * t205 - t207 * t220 + t209 * t221) * t382 + t406 * t404) * t213, ((t337 * t206 - t210 * t220 + t212 * t221) * t382 + t405 * t404) * t213, 0; ((-t313 * t402 - t388 * t317) * t224 - (t388 * t313 - t317 * t402) * t379) * t219, t406 * t403, t405 * t403, t350 * t219;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:15
	% EndTime: 2020-04-19 18:18:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:16
	% EndTime: 2020-04-19 18:18:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:19:11
	% EndTime: 2020-04-19 18:19:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end