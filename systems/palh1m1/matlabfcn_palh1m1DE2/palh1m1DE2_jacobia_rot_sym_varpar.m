% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m1DE2
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
%   Wie in palh1m1DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh1m1DE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_jacobia_rot_sym_varpar: pkin has to be [23x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:50:09
	% EndTime: 2020-04-15 18:50:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:50:32
	% EndTime: 2020-04-15 18:51:36
	% DurationCPUTime: 37.75s
	% Computational Cost: add. (1058724->109), mult. (1593950->224), div. (72763->20), fcn. (1009097->24), ass. (0->145)
	t313 = -2 * pkin(1);
	t312 = -2 * pkin(5);
	t311 = pkin(4) * pkin(5);
	t310 = -pkin(8) - pkin(3);
	t309 = -pkin(8) + pkin(3);
	t308 = (-pkin(9) - pkin(11));
	t307 = (-pkin(9) + pkin(11));
	t231 = sin(pkin(21));
	t306 = t231 / 0.2e1;
	t238 = cos(qJ(3));
	t305 = -t238 / 0.2e1;
	t304 = cos(qJ(2));
	t303 = sin(pkin(19));
	t242 = pkin(7) ^ 2;
	t235 = sin(qJ(2));
	t240 = cos(pkin(19));
	t225 = t235 * t240 - t304 * t303;
	t300 = pkin(7) * t225;
	t272 = (pkin(1) ^ 2) + t300 * t313;
	t222 = t242 + t272;
	t219 = pkin(3) ^ 2 - pkin(8) ^ 2 + t222;
	t223 = pkin(1) - t300;
	t217 = (pkin(7) - t310) * (pkin(7) + t310) + t272;
	t218 = (pkin(7) - t309) * (pkin(7) + t309) + t272;
	t248 = sqrt(-t218 * t217);
	t226 = t235 * t303 + t304 * t240;
	t299 = pkin(7) * t226;
	t213 = t219 * t299 + t223 * t248;
	t276 = t238 * t213;
	t283 = t226 * t248;
	t212 = -pkin(7) * t283 + t223 * t219;
	t234 = sin(qJ(3));
	t282 = t234 * t212;
	t220 = 0.1e1 / t222;
	t245 = 0.1e1 / pkin(3);
	t284 = t220 * t245;
	t208 = (t282 / 0.2e1 + t276 / 0.2e1) * t284;
	t277 = t238 * t212;
	t281 = t234 * t213;
	t209 = (-t277 / 0.2e1 + t281 / 0.2e1) * t284;
	t230 = pkin(23) + pkin(22);
	t228 = sin(t230);
	t229 = cos(t230);
	t201 = t229 * t208 - t228 * t209;
	t302 = pkin(4) * t201;
	t258 = t228 * t208 + t229 * t209;
	t301 = pkin(4) * t258;
	t244 = pkin(4) ^ 2;
	t273 = (pkin(5) ^ 2) - t302 * t312;
	t198 = t244 + t273;
	t195 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t198;
	t199 = pkin(5) + t302;
	t193 = ((pkin(4) - t308) * (pkin(4) + t308)) + t273;
	t194 = ((pkin(4) - t307) * (pkin(4) + t307)) + t273;
	t247 = sqrt(-t194 * t193);
	t184 = t199 * t195 - t247 * t301;
	t185 = t195 * t301 + t199 * t247;
	t232 = cos(pkin(21));
	t196 = 0.1e1 / t198;
	t241 = 0.1e1 / pkin(11);
	t286 = t196 * t241;
	t181 = (t184 * t306 + t185 * t232 / 0.2e1) * t286;
	t182 = (-t184 * t232 / 0.2e1 + t185 * t306) * t286;
	t175 = qJ(2) + qJ(3) + atan2(t181, t182);
	t173 = sin(t175);
	t170 = t173 ^ 2;
	t174 = cos(t175);
	t236 = sin(qJ(1));
	t280 = t236 * t173;
	t162 = atan2(-t280, -t174);
	t160 = sin(t162);
	t161 = cos(t162);
	t158 = -t160 * t280 - t161 * t174;
	t157 = 0.1e1 / t158 ^ 2;
	t239 = cos(qJ(1));
	t296 = t157 * t239 ^ 2;
	t155 = 0.1e1 / (t170 * t296 + 0.1e1);
	t298 = t155 * t239;
	t297 = t157 * t173;
	t295 = t160 * t174;
	t292 = t170 / t174 ^ 2;
	t163 = 0.1e1 / (t236 ^ 2 * t292 + 0.1e1);
	t294 = t163 * t236;
	t237 = cos(qJ(4));
	t274 = t239 * t237;
	t233 = sin(qJ(4));
	t279 = t236 * t233;
	t169 = t174 * t274 + t279;
	t167 = 0.1e1 / t169 ^ 2;
	t275 = t239 * t233;
	t278 = t236 * t237;
	t168 = t174 * t275 - t278;
	t293 = t167 * t168;
	t291 = t173 * t239;
	t156 = 0.1e1 / t158;
	t290 = t174 * t156;
	t180 = 0.1e1 / t182 ^ 2;
	t289 = 0.1e1 / (t181 ^ 2 * t180 + 0.1e1) * t241;
	t288 = t180 * t181;
	t287 = t196 * t232;
	t271 = pkin(1) * t299;
	t285 = 0.2e1 / t248 * (t217 + t218) * t271;
	t270 = 0.1e1 / t198 ^ 2 * t311;
	t190 = 0.1e1 / t247;
	t269 = t190 * t199 / 0.2e1;
	t268 = -t190 * t258 / 0.2e1;
	t267 = t196 * t306;
	t266 = -t287 / 0.2e1;
	t265 = t287 / 0.2e1;
	t264 = t167 * t168 ^ 2 + 0.1e1;
	t263 = t258 * t244 * t312;
	t262 = t199 * t312 - t195;
	t261 = 0.1e1 / t222 ^ 2 * t271;
	t260 = t231 * t270;
	t259 = t232 * t270;
	t257 = t184 * t259;
	t256 = t185 * t259;
	t205 = (t225 * t248 + (-t285 / 0.2e1 - t219 + t223 * t313) * t226) * pkin(7);
	t206 = t223 * t285 / 0.2e1 + t242 * t226 ^ 2 * t313 + (-t225 * t219 - t283) * pkin(7);
	t191 = ((-t234 * t205 / 0.2e1 + t206 * t305) * t220 + (-t276 - t282) * t261) * t245;
	t192 = ((t205 * t305 + t234 * t206 / 0.2e1) * t220 + (-t277 + t281) * t261) * t245;
	t187 = t229 * t191 + t228 * t192;
	t255 = t187 * t260;
	t254 = t258 * t260;
	t253 = 0.2e1 * (t193 + t194) * t311;
	t252 = (0.1e1 + t292) * t294;
	t251 = -t161 * t280 + t295;
	t250 = t161 * t173 - t236 * t295;
	t159 = 0.1e1 / t264;
	t166 = 0.1e1 / t169;
	t249 = (-t166 * t233 + t237 * t293) * t159 * t291;
	t188 = -t228 * t191 + t229 * t192;
	t186 = t258 * t253;
	t183 = t187 * t253;
	t179 = 0.1e1 / t182;
	t178 = t186 * t269 + t258 * t263 + (t201 * t195 - t247 * t258) * pkin(4);
	t177 = (t186 * t268 - t201 * t247 + t258 * t262) * pkin(4);
	t171 = 0.1e1 / t174;
	t165 = t183 * t269 + t187 * t263 + (-t187 * t247 + t188 * t195) * pkin(4);
	t164 = (t183 * t268 + t262 * t187 - t188 * t247) * pkin(4);
	t154 = 0.1e1 + ((t177 * t267 + t178 * t265 + t184 * t254 + t256 * t258) * t179 - (t177 * t266 + t178 * t267 + t185 * t254 - t257 * t258) * t288) * t289;
	t153 = 0.1e1 + ((t164 * t267 + t165 * t265 + t184 * t255 + t187 * t256) * t179 - (t164 * t266 + t165 * t267 + t185 * t255 - t187 * t257) * t288) * t289;
	t152 = t154 * t252;
	t151 = t153 * t252;
	t1 = [t171 * t163 * t291, t151, t152, 0; (-t156 * t280 - (-t161 * t170 * t171 * t294 + (t163 - 0.1e1) * t173 * t160) * t173 * t296) * t155, (t153 * t290 - (t251 * t151 + t250 * t153) * t297) * t298, (t154 * t290 - (t251 * t152 + t250 * t154) * t297) * t298, 0; ((-t174 * t279 - t274) * t166 - (-t174 * t278 + t275) * t293) * t159, t153 * t249, t154 * t249, t264 * t159;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:56
	% EndTime: 2020-04-15 18:49:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end