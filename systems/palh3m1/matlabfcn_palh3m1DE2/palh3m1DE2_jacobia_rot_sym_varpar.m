% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m1DE2
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
%   Wie in palh3m1DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh3m1DE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_jacobia_rot_sym_varpar: pkin has to be [19x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:19
	% EndTime: 2020-04-20 16:21:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:54
	% EndTime: 2020-04-20 16:23:20
	% DurationCPUTime: 38.86s
	% Computational Cost: add. (1058654->107), mult. (1593950->224), div. (72763->20), fcn. (1009097->24), ass. (0->145)
	t314 = -2 * pkin(1);
	t313 = -2 * pkin(4);
	t312 = pkin(3) * pkin(4);
	t311 = -pkin(6) - pkin(2);
	t310 = -pkin(6) + pkin(2);
	t309 = (-pkin(8) - pkin(10));
	t308 = (-pkin(8) + pkin(10));
	t232 = sin(pkin(17));
	t307 = t232 / 0.2e1;
	t235 = sin(qJ(3));
	t306 = t235 / 0.2e1;
	t305 = cos(pkin(16));
	t304 = cos(qJ(2));
	t243 = pkin(5) ^ 2;
	t236 = sin(qJ(2));
	t238 = sin(pkin(16));
	t226 = t236 * t238 - t304 * t305;
	t301 = pkin(5) * t226;
	t273 = (pkin(1) ^ 2) + t301 * t314;
	t223 = t243 + t273;
	t221 = 0.1e1 / t223;
	t246 = 0.1e1 / pkin(2);
	t281 = t221 * t246;
	t220 = pkin(2) ^ 2 - pkin(6) ^ 2 + t223;
	t224 = pkin(1) - t301;
	t218 = (pkin(5) - t311) * (pkin(5) + t311) + t273;
	t219 = (pkin(5) - t310) * (pkin(5) + t310) + t273;
	t249 = sqrt(-t219 * t218);
	t227 = t236 * t305 + t238 * t304;
	t300 = pkin(5) * t227;
	t214 = t220 * t300 + t224 * t249;
	t284 = t214 * t235;
	t280 = t227 * t249;
	t213 = -pkin(5) * t280 + t220 * t224;
	t240 = cos(qJ(3));
	t285 = t213 * t240;
	t209 = (-t285 / 0.2e1 + t284 / 0.2e1) * t281;
	t283 = t214 * t240;
	t286 = t213 * t235;
	t210 = (t283 / 0.2e1 + t286 / 0.2e1) * t281;
	t231 = pkin(18) + pkin(19);
	t229 = sin(t231);
	t230 = cos(t231);
	t255 = t209 * t229 - t210 * t230;
	t303 = pkin(3) * t255;
	t204 = t209 * t230 + t210 * t229;
	t302 = pkin(3) * t204;
	t245 = pkin(3) ^ 2;
	t274 = (pkin(4) ^ 2) - t302 * t313;
	t200 = t245 + t274;
	t197 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t200;
	t201 = pkin(4) + t302;
	t195 = ((pkin(3) - t309) * (pkin(3) + t309)) + t274;
	t196 = ((pkin(3) - t308) * (pkin(3) + t308)) + t274;
	t248 = sqrt(-t196 * t195);
	t186 = t197 * t201 - t248 * t303;
	t187 = t197 * t303 + t201 * t248;
	t233 = cos(pkin(17));
	t198 = 0.1e1 / t200;
	t242 = 0.1e1 / pkin(10);
	t287 = t198 * t242;
	t183 = (-t186 * t233 / 0.2e1 + t187 * t307) * t287;
	t184 = (t187 * t233 / 0.2e1 + t186 * t307) * t287;
	t177 = qJ(2) + qJ(3) + atan2(t184, t183);
	t175 = sin(t177);
	t172 = t175 ^ 2;
	t176 = cos(t177);
	t237 = sin(qJ(1));
	t277 = t237 * t175;
	t164 = atan2(t277, t176);
	t161 = sin(t164);
	t162 = cos(t164);
	t159 = t161 * t277 + t162 * t176;
	t158 = 0.1e1 / t159 ^ 2;
	t241 = cos(qJ(1));
	t297 = t158 * t241 ^ 2;
	t156 = 0.1e1 / (t172 * t297 + 0.1e1);
	t299 = t156 * t241;
	t298 = t158 * t175;
	t296 = t161 * t176;
	t293 = t172 / t176 ^ 2;
	t163 = 0.1e1 / (t237 ^ 2 * t293 + 0.1e1);
	t295 = t163 * t237;
	t239 = cos(qJ(4));
	t275 = t239 * t241;
	t234 = sin(qJ(4));
	t279 = t234 * t237;
	t171 = -t176 * t275 + t279;
	t168 = 0.1e1 / t171 ^ 2;
	t276 = t237 * t239;
	t278 = t234 * t241;
	t169 = t176 * t278 + t276;
	t294 = t168 * t169;
	t292 = t175 * t241;
	t157 = 0.1e1 / t159;
	t291 = t176 * t157;
	t182 = 0.1e1 / t183 ^ 2;
	t290 = 0.1e1 / (t182 * t184 ^ 2 + 0.1e1) * t242;
	t289 = t182 * t184;
	t288 = t198 * t233;
	t272 = pkin(1) * t300;
	t282 = 0.2e1 / t249 * (t218 + t219) * t272;
	t271 = 0.1e1 / t200 ^ 2 * t312;
	t192 = 0.1e1 / t248;
	t270 = t192 * t201 / 0.2e1;
	t269 = -t192 * t255 / 0.2e1;
	t268 = t198 * t307;
	t267 = -t288 / 0.2e1;
	t266 = t288 / 0.2e1;
	t265 = t168 * t169 ^ 2 + 0.1e1;
	t264 = t245 * t255 * t313;
	t263 = t201 * t313 - t197;
	t262 = 0.1e1 / t223 ^ 2 * t272;
	t207 = (t226 * t249 + (-t282 / 0.2e1 - t220 + t224 * t314) * t227) * pkin(5);
	t208 = t224 * t282 / 0.2e1 + t243 * t227 ^ 2 * t314 + (-t226 * t220 - t280) * pkin(5);
	t193 = ((t208 * t240 / 0.2e1 + t207 * t306) * t221 + (t283 + t286) * t262) * t246;
	t194 = ((-t207 * t240 / 0.2e1 + t208 * t306) * t221 + (t284 - t285) * t262) * t246;
	t190 = -t193 * t229 - t194 * t230;
	t261 = t190 * t271;
	t260 = t255 * t271;
	t259 = t186 * t261;
	t258 = t187 * t261;
	t257 = t232 * t260;
	t256 = t233 * t260;
	t254 = 0.2e1 * (t195 + t196) * t312;
	t253 = (0.1e1 + t293) * t295;
	t252 = t162 * t277 - t296;
	t251 = -t162 * t175 + t237 * t296;
	t160 = 0.1e1 / t265;
	t167 = 0.1e1 / t171;
	t250 = (t167 * t234 + t239 * t294) * t160 * t292;
	t189 = -t193 * t230 + t194 * t229;
	t188 = t255 * t254;
	t185 = t190 * t254;
	t181 = 0.1e1 / t183;
	t180 = t188 * t270 + t255 * t264 + (t204 * t197 - t248 * t255) * pkin(3);
	t179 = (t188 * t269 - t204 * t248 + t255 * t263) * pkin(3);
	t173 = 0.1e1 / t176;
	t166 = t185 * t270 + t190 * t264 + (t189 * t197 - t190 * t248) * pkin(3);
	t165 = (t185 * t269 - t189 * t248 + t190 * t263) * pkin(3);
	t155 = 0.1e1 + ((t179 * t268 + t180 * t266 + t186 * t257 + t187 * t256) * t181 - (t179 * t267 + t180 * t268 - t186 * t256 + t187 * t257) * t289) * t290;
	t154 = 0.1e1 + ((t165 * t268 + t166 * t266 + t232 * t259 + t233 * t258) * t181 - (t165 * t267 + t166 * t268 + t232 * t258 - t233 * t259) * t289) * t290;
	t153 = t155 * t253;
	t152 = t154 * t253;
	t1 = [t173 * t163 * t292, t152, t153, 0; (t157 * t277 + (t162 * t172 * t173 * t295 + (-t163 + 0.1e1) * t175 * t161) * t175 * t297) * t156, (-t154 * t291 + (t152 * t252 + t154 * t251) * t298) * t299, (-t155 * t291 + (t153 * t252 + t155 * t251) * t298) * t299, 0; ((t176 * t279 - t275) * t167 + (t176 * t276 + t278) * t294) * t160, t154 * t250, t155 * t250, t265 * t160;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:43
	% EndTime: 2020-04-20 16:20:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:43
	% EndTime: 2020-04-20 16:20:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:03
	% EndTime: 2020-04-20 16:21:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end