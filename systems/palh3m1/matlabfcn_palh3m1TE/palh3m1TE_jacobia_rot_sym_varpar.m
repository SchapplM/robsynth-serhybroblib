% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m1TE
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
%   Wie in palh3m1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh3m1TE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1TE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_jacobia_rot_sym_varpar: pkin has to be [19x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:20
	% EndTime: 2020-04-18 09:52:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:49
	% EndTime: 2020-04-18 09:53:57
	% DurationCPUTime: 40.83s
	% Computational Cost: add. (896623->129), mult. (1352046->244), div. (61037->17), fcn. (857228->21), ass. (0->145)
	t336 = sin(qJ(2));
	t337 = sin(pkin(16));
	t339 = cos(qJ(2));
	t341 = cos(pkin(16));
	t243 = t336 * t341 + t339 * t337;
	t303 = t336 * t337 - t339 * t341;
	t302 = pkin(5) * t303;
	t299 = (-0.2e1 * t302 + pkin(1)) * pkin(1);
	t344 = -pkin(6) + pkin(2);
	t345 = -pkin(6) - pkin(2);
	t253 = sqrt(-((pkin(5) - t344) * (pkin(5) + t344) + t299) * ((pkin(5) - t345) * (pkin(5) + t345) + t299));
	t353 = pkin(5) ^ 2;
	t236 = t299 + t353;
	t254 = pkin(2) ^ 2;
	t362 = -pkin(6) ^ 2 + t254;
	t298 = t236 + t362;
	t297 = pkin(5) * t298;
	t300 = -t302 + pkin(1);
	t354 = 0.1e1 / pkin(2);
	t296 = t354 * (t243 * t297 + t300 * t253);
	t347 = 0.1e1 / t236;
	t293 = t347 * t296 / 0.2e1;
	t332 = pkin(5) * t243;
	t295 = t354 * (-t253 * t332 + t300 * t298);
	t294 = t347 * t295;
	t335 = sin(qJ(3));
	t338 = cos(qJ(3));
	t289 = -t338 * t294 / 0.2e1 + t335 * t293;
	t364 = t335 / 0.2e1;
	t290 = t338 * t293 + t294 * t364;
	t317 = pkin(18) + pkin(19);
	t311 = sin(t317);
	t312 = cos(t317);
	t283 = t312 * t289 + t311 * t290;
	t281 = pkin(3) * t283;
	t278 = t281 + pkin(4);
	t368 = t278 / 0.2e1;
	t355 = t311 * t289 - t312 * t290;
	t280 = pkin(3) * t355;
	t367 = -t280 / 0.2e1;
	t329 = sin(pkin(17));
	t366 = t329 / 0.2e1;
	t330 = cos(pkin(17));
	t365 = -t330 / 0.2e1;
	t349 = -0.2e1 * pkin(4);
	t229 = pkin(4) ^ 2 + (-t283 * t349 + pkin(3)) * pkin(3);
	t346 = pkin(3) * pkin(4);
	t363 = 0.1e1 / t229 ^ 2 * t346;
	t277 = (0.2e1 * t281 + pkin(4)) * pkin(4);
	t342 = -pkin(8) + pkin(10);
	t343 = -pkin(8) - pkin(10);
	t252 = sqrt(-((pkin(3) - t342) * (pkin(3) + t342) + t277) * ((pkin(3) - t343) * (pkin(3) + t343) + t277));
	t351 = pkin(8) ^ 2;
	t275 = pkin(10) ^ 2 + t229 - t351;
	t350 = 0.1e1 / pkin(10);
	t268 = t350 * (-t252 * t280 + t278 * t275);
	t264 = t268 * t363;
	t274 = pkin(3) * t275;
	t269 = t350 * (t278 * t252 + t274 * t355);
	t265 = t269 * t363;
	t361 = -t330 * t264 + t329 * t265;
	t360 = t329 * t264 + t330 * t265;
	t359 = -0.2e1 * t278 * t346 - t274;
	t241 = t336 * t335 - t339 * t338;
	t358 = t339 * t335 + t336 * t338;
	t333 = pkin(3) * t252;
	t357 = pkin(3) ^ 2 * t349 * t355 - t333;
	t356 = 0.4e1 / t252 * ((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t277 - t351) * t346;
	t348 = 0.1e1 / t229;
	t340 = cos(qJ(1));
	t334 = t354 * t347;
	t331 = t350 * t348;
	t266 = t348 * t269 / 0.2e1;
	t267 = t348 * t268;
	t223 = t329 * t266 + t267 * t365;
	t224 = t330 * t266 + t267 * t366;
	t250 = sin(qJ(1));
	t237 = t358 * t250;
	t238 = t241 * t250;
	t213 = t223 * t237 - t224 * t238;
	t219 = -t223 * t241 - t224 * t358;
	t203 = atan2(t213, t219);
	t200 = sin(t203);
	t201 = cos(t203);
	t193 = t200 * t213 + t201 * t219;
	t192 = 0.1e1 / t193 ^ 2;
	t239 = t241 * t340;
	t222 = t239 * t224;
	t240 = t358 * t340;
	t215 = t223 * t240 - t222;
	t328 = t192 * t215;
	t327 = t192 * t215 ^ 2;
	t326 = t201 * t213;
	t214 = t223 * t239 + t224 * t240;
	t249 = sin(qJ(4));
	t251 = cos(qJ(4));
	t211 = t214 * t251 + t250 * t249;
	t209 = 0.1e1 / t211 ^ 2;
	t210 = t214 * t249 - t250 * t251;
	t325 = t209 * t210;
	t218 = 0.1e1 / t219 ^ 2;
	t324 = t213 * t218;
	t316 = pkin(1) * t332;
	t323 = 0.4e1 / t253 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t299 - t254) * t316;
	t231 = (t303 * t253 + (-t323 / 0.2e1 - t353 - t362) * t243) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t302) * t316;
	t232 = t300 * t323 / 0.2e1 - t303 * t297 + (-t253 - 0.2e1 * t316) * t332;
	t235 = 0.1e1 / t236 ^ 2;
	t291 = t295 * t316;
	t292 = t296 * t316;
	t305 = t334 * t364;
	t310 = t338 * t334;
	t285 = t232 * t310 / 0.2e1 + t231 * t305 + (t335 * t291 + t338 * t292) * t235;
	t286 = -t231 * t310 / 0.2e1 + t232 * t305 + (-t338 * t291 + t335 * t292) * t235;
	t227 = -t311 * t285 - t312 * t286;
	t270 = pkin(3) * (-t312 * t285 + t311 * t286);
	t272 = t227 * t356;
	t314 = t331 / 0.2e1;
	t256 = (t357 * t227 + t275 * t270 + t272 * t368) * t314;
	t257 = (t359 * t227 - t252 * t270 + t272 * t367) * t331;
	t197 = t361 * t227 + t329 * t256 + t257 * t365;
	t322 = t197 - t224;
	t321 = t360 * t227 + t330 * t256 + t257 * t366 + t223;
	t271 = t355 * t356;
	t258 = (t271 * t368 + t283 * t274 + t355 * t357) * t314;
	t259 = (t271 * t367 - t283 * t333 + t355 * t359) * t331;
	t204 = t329 * t258 + t259 * t365 + t355 * t361;
	t320 = t204 - t224;
	t319 = t330 * t258 + t259 * t366 + t355 * t360 + t223;
	t313 = t209 * t210 ^ 2 + 0.1e1;
	t304 = -t200 * t219 + t326;
	t199 = 0.1e1 / t313;
	t208 = 0.1e1 / t211;
	t301 = (t208 * t249 - t251 * t325) * t199;
	t217 = 0.1e1 / t219;
	t212 = -t223 * t238 - t224 * t237;
	t202 = 0.1e1 / (t213 ^ 2 * t218 + 0.1e1);
	t196 = -t320 * t241 - t319 * t358;
	t194 = t320 * t237 - t319 * t238;
	t191 = 0.1e1 / t193;
	t190 = -t322 * t241 - t321 * t358;
	t188 = t322 * t237 - t321 * t238;
	t187 = 0.1e1 / (0.1e1 + t327);
	t186 = (t194 * t217 - t196 * t324) * t202;
	t185 = (t188 * t217 - t190 * t324) * t202;
	t1 = [t215 * t217 * t202, t185, t186, 0; (t213 * t191 + (t200 + (t217 * t326 - t200) * t202) * t327) * t187, ((t321 * t239 - t322 * t240) * t191 + (t304 * t185 + t188 * t200 + t190 * t201) * t328) * t187, ((t319 * t239 - t320 * t240) * t191 + (t304 * t186 + t194 * t200 + t196 * t201) * t328) * t187, 0; ((t212 * t249 - t340 * t251) * t208 - (t212 * t251 + t340 * t249) * t325) * t199, (t197 * t239 + t321 * t240 - t222) * t301, (t204 * t239 + t319 * t240 - t222) * t301, t313 * t199;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:03
	% EndTime: 2020-04-18 09:52:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:04
	% EndTime: 2020-04-18 09:52:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:26
	% EndTime: 2020-04-18 09:52:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end