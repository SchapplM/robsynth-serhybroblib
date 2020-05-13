% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m1TE
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
%   Wie in palh1m1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh1m1TE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1TE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_jacobia_rot_sym_varpar: pkin has to be [23x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:11
	% EndTime: 2020-04-13 14:18:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:29
	% EndTime: 2020-04-13 14:19:10
	% DurationCPUTime: 40.82s
	% Computational Cost: add. (896623->129), mult. (1352046->244), div. (61037->17), fcn. (857228->21), ass. (0->145)
	t335 = sin(qJ(2));
	t336 = sin(pkin(19));
	t338 = cos(qJ(2));
	t339 = cos(pkin(19));
	t243 = t335 * t336 + t338 * t339;
	t304 = t335 * t339 - t338 * t336;
	t303 = pkin(7) * t304;
	t300 = (-0.2e1 * t303 + pkin(1)) * pkin(1);
	t342 = -pkin(8) + pkin(3);
	t343 = -pkin(8) - pkin(3);
	t254 = sqrt(-((pkin(7) - t342) * (pkin(7) + t342) + t300) * ((pkin(7) - t343) * (pkin(7) + t343) + t300));
	t351 = pkin(7) ^ 2;
	t236 = t300 + t351;
	t255 = pkin(3) ^ 2;
	t359 = -pkin(8) ^ 2 + t255;
	t299 = t236 + t359;
	t298 = pkin(7) * t299;
	t301 = -t303 + pkin(1);
	t352 = 0.1e1 / pkin(3);
	t297 = t352 * (t243 * t298 + t301 * t254);
	t345 = 0.1e1 / t236;
	t294 = t345 * t297 / 0.2e1;
	t331 = pkin(7) * t243;
	t296 = t352 * (-t254 * t331 + t301 * t299);
	t295 = t345 * t296;
	t334 = sin(qJ(3));
	t337 = cos(qJ(3));
	t290 = t334 * t295 / 0.2e1 + t337 * t294;
	t361 = -t337 / 0.2e1;
	t291 = t334 * t294 + t295 * t361;
	t314 = pkin(23) + pkin(22);
	t310 = sin(t314);
	t311 = cos(t314);
	t283 = t311 * t290 - t310 * t291;
	t281 = pkin(4) * t283;
	t279 = t281 + pkin(5);
	t365 = t279 / 0.2e1;
	t353 = t310 * t290 + t311 * t291;
	t282 = pkin(4) * t353;
	t364 = -t282 / 0.2e1;
	t328 = sin(pkin(21));
	t363 = t328 / 0.2e1;
	t329 = cos(pkin(21));
	t362 = -t329 / 0.2e1;
	t347 = -0.2e1 * pkin(5);
	t229 = pkin(5) ^ 2 + (-t283 * t347 + pkin(4)) * pkin(4);
	t344 = pkin(4) * pkin(5);
	t360 = 0.1e1 / t229 ^ 2 * t344;
	t278 = (0.2e1 * t281 + pkin(5)) * pkin(5);
	t340 = -pkin(9) + pkin(11);
	t341 = -pkin(9) - pkin(11);
	t253 = sqrt(-((pkin(4) - t340) * (pkin(4) + t340) + t278) * ((pkin(4) - t341) * (pkin(4) + t341) + t278));
	t349 = pkin(9) ^ 2;
	t276 = pkin(11) ^ 2 + t229 - t349;
	t348 = 0.1e1 / pkin(11);
	t269 = t348 * (-t253 * t282 + t279 * t276);
	t265 = t269 * t360;
	t275 = pkin(4) * t276;
	t270 = t348 * (t279 * t253 + t275 * t353);
	t266 = t270 * t360;
	t358 = t328 * t265 + t329 * t266;
	t357 = -t329 * t265 + t328 * t266;
	t356 = -0.2e1 * t279 * t344 - t275;
	t242 = -t335 * t334 + t338 * t337;
	t332 = pkin(4) * t253;
	t355 = pkin(4) ^ 2 * t347 * t353 - t332;
	t354 = 0.4e1 / t253 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t278 - t349) * t344;
	t346 = 0.1e1 / t229;
	t333 = t352 * t345;
	t330 = t348 * t346;
	t267 = t346 * t270 / 0.2e1;
	t268 = t346 * t269;
	t223 = t329 * t267 + t268 * t363;
	t224 = t328 * t267 + t268 * t362;
	t250 = sin(qJ(1));
	t237 = t242 * t250;
	t241 = -t338 * t334 - t335 * t337;
	t238 = t241 * t250;
	t213 = -t223 * t237 + t224 * t238;
	t219 = -t223 * t241 - t224 * t242;
	t203 = atan2(t213, t219);
	t200 = sin(t203);
	t201 = cos(t203);
	t193 = t200 * t213 + t201 * t219;
	t192 = 0.1e1 / t193 ^ 2;
	t252 = cos(qJ(1));
	t239 = t242 * t252;
	t222 = t239 * t223;
	t240 = t241 * t252;
	t214 = t224 * t240 - t222;
	t327 = t192 * t214;
	t326 = t192 * t214 ^ 2;
	t325 = t201 * t213;
	t216 = t223 * t240 + t224 * t239;
	t249 = sin(qJ(4));
	t251 = cos(qJ(4));
	t211 = t216 * t251 + t249 * t250;
	t209 = 0.1e1 / t211 ^ 2;
	t210 = t216 * t249 - t250 * t251;
	t324 = t209 * t210;
	t218 = 0.1e1 / t219 ^ 2;
	t323 = t213 * t218;
	t316 = pkin(1) * t331;
	t322 = 0.4e1 / t254 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t300 - t255) * t316;
	t231 = (t304 * t254 + (-t322 / 0.2e1 - t351 - t359) * t243) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t303) * t316;
	t232 = t301 * t322 / 0.2e1 - t304 * t298 + (-t254 - 0.2e1 * t316) * t331;
	t235 = 0.1e1 / t236 ^ 2;
	t292 = t296 * t316;
	t293 = t297 * t316;
	t306 = t333 * t361;
	t309 = t334 * t333;
	t286 = -t231 * t309 / 0.2e1 + t232 * t306 + (-t334 * t292 - t337 * t293) * t235;
	t287 = t231 * t306 + t232 * t309 / 0.2e1 + (-t337 * t292 + t334 * t293) * t235;
	t227 = t311 * t286 + t310 * t287;
	t271 = pkin(4) * (-t310 * t286 + t311 * t287);
	t273 = t227 * t354;
	t313 = t330 / 0.2e1;
	t257 = (t355 * t227 + t276 * t271 + t273 * t365) * t313;
	t258 = (t356 * t227 - t253 * t271 + t273 * t364) * t330;
	t321 = t358 * t227 + t329 * t257 + t258 * t363 + t224;
	t198 = t357 * t227 + t328 * t257 + t258 * t362;
	t320 = -t198 + t223;
	t272 = t353 * t354;
	t259 = (t272 * t365 + t283 * t275 + t353 * t355) * t313;
	t260 = (t272 * t364 - t283 * t332 + t353 * t356) * t330;
	t319 = t329 * t259 + t260 * t363 + t353 * t358 + t224;
	t205 = t328 * t259 + t260 * t362 + t353 * t357;
	t318 = -t205 + t223;
	t312 = t209 * t210 ^ 2 + 0.1e1;
	t305 = -t200 * t219 + t325;
	t199 = 0.1e1 / t312;
	t208 = 0.1e1 / t211;
	t302 = (t208 * t249 - t251 * t324) * t199;
	t217 = 0.1e1 / t219;
	t212 = -t223 * t238 - t224 * t237;
	t202 = 0.1e1 / (t213 ^ 2 * t218 + 0.1e1);
	t196 = -t319 * t241 + t318 * t242;
	t194 = -t319 * t237 - t318 * t238;
	t191 = 0.1e1 / t193;
	t190 = -t321 * t241 + t320 * t242;
	t188 = -t321 * t237 - t320 * t238;
	t187 = 0.1e1 / (0.1e1 + t326);
	t186 = (t194 * t217 - t196 * t323) * t202;
	t185 = (t188 * t217 - t190 * t323) * t202;
	t1 = [t214 * t217 * t202, t185, t186, 0; (t213 * t191 + (t200 + (t217 * t325 - t200) * t202) * t326) * t187, ((t321 * t239 + t320 * t240) * t191 + (t305 * t185 + t188 * t200 + t190 * t201) * t327) * t187, ((t319 * t239 + t318 * t240) * t191 + (t305 * t186 + t194 * t200 + t196 * t201) * t327) * t187, 0; ((t212 * t249 - t251 * t252) * t208 - (t212 * t251 + t249 * t252) * t324) * t199, (t198 * t239 + t321 * t240 - t222) * t302, (t205 * t239 + t319 * t240 - t222) * t302, t312 * t199;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:16
	% EndTime: 2020-04-13 14:18:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end