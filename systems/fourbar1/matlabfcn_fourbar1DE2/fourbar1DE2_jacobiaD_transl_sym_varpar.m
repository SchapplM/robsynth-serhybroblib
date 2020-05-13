% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% fourbar1DE2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% JaD_transl [3x1]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = fourbar1DE2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_jacobiaD_transl_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE2_jacobiaD_transl_sym_varpar: qJD has to be [1x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fourbar1DE2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1DE2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_jacobiaD_transl_sym_varpar: pkin has to be [4x1] (double)');
JaD_transl=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 20:05:27
	% EndTime: 2020-04-24 20:05:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 20:05:28
	% EndTime: 2020-04-24 20:05:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1); (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1); 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 20:05:29
	% EndTime: 2020-04-24 20:05:29
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (1160->54), mult. (1736->97), div. (40->5), fcn. (436->4), ass. (0->58)
	t252 = pkin(2) ^ 2;
	t253 = pkin(1) ^ 2;
	t249 = cos(qJ(1));
	t292 = pkin(2) * t249;
	t298 = -2 * pkin(1);
	t284 = t292 * t298 + t253;
	t244 = t252 + t284;
	t241 = 0.1e1 / t244;
	t296 = pkin(1) * pkin(2);
	t251 = 0.1e1 / pkin(3);
	t302 = 0.1e1 / t244 ^ 2 * t251;
	t276 = t296 * t302;
	t271 = qJD(1) * t276;
	t248 = sin(qJ(1));
	t247 = t248 ^ 2;
	t283 = qJD(1) * t252;
	t275 = 0.4e1 * t247 * t283;
	t305 = t241 * t253 * t275 * t302 - t249 * t271;
	t295 = -pkin(3) - pkin(4);
	t238 = (pkin(2) - t295) * (pkin(2) + t295) + t284;
	t294 = -pkin(3) + pkin(4);
	t239 = (pkin(2) - t294) * (pkin(2) + t294) + t284;
	t274 = (-t238 - t239) * t296;
	t229 = t248 * t274;
	t228 = qJD(1) * t229;
	t240 = pkin(3) ^ 2 - pkin(4) ^ 2 + t244;
	t291 = qJD(1) * pkin(2);
	t280 = t249 * t291;
	t236 = t240 * t280;
	t293 = pkin(2) * t248;
	t237 = t240 * t293;
	t245 = -pkin(1) + t292;
	t289 = t238 * t239;
	t254 = sqrt(-t289);
	t287 = t248 * t254;
	t231 = pkin(2) * t287;
	t286 = t252 * t247;
	t262 = 0.2e1 * pkin(1) * t286 - t231;
	t266 = t248 * t271;
	t272 = t248 * t276;
	t233 = 0.1e1 / t254;
	t290 = t233 * t245;
	t304 = -(t262 * qJD(1) + t228 * t290 + t236) * t272 + t305 * (t245 * t254 + t237) - (t229 * t290 + t240 * t292 + t262) * t266;
	t303 = t233 * t229;
	t223 = t233 * t228 * t293;
	t227 = (t249 * t274 - 0.4e1 * t253 * t286) * qJD(1);
	t268 = 0.1e1 / t289 * t228 * t303;
	t278 = t228 / 0.2e1;
	t269 = 0.2e1 * t278;
	t279 = t241 * t251 / 0.2e1;
	t281 = t248 * t291;
	t282 = t245 * t298;
	t285 = -t240 * t281 - t254 * t280;
	t297 = 0.2e1 * t233;
	t300 = -(t237 + (t249 * t254 + (t282 + t303) * t248) * pkin(2)) * t266 - (t281 * t282 + t223 - t285) * t272 + (pkin(1) * t275 + t236 + (t248 * t268 + (t249 * t282 - t287) * qJD(1) + (t248 * t227 / 0.2e1 + t269 * t249) * t297) * pkin(2)) * t279 + t305 * (-t245 * t240 + t231);
	t288 = t245 * t227;
	t258 = 0.6e1 * t249 * t248 * pkin(1) * t283 + t245 * t268 + t285;
	t1 = [-t280 + (((t288 / 0.2e1 - t269 * t293) * t297 + t258) * t279 + t304) * r_i_i_C(2) + t300 * r_i_i_C(1); -t281 + t300 * r_i_i_C(2) + ((t223 + (t278 * t293 - t288 / 0.2e1) * t297 - t258) * t279 - t304) * r_i_i_C(1); 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 20:05:29
	% EndTime: 2020-04-24 20:05:29
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (1158->56), mult. (1732->99), div. (40->5), fcn. (434->4), ass. (0->56)
	t247 = pkin(2) ^ 2;
	t248 = pkin(1) ^ 2;
	t244 = cos(qJ(1));
	t287 = pkin(2) * t244;
	t279 = -0.2e1 * pkin(1) * t287 + t248;
	t239 = t247 + t279;
	t236 = 0.1e1 / t239;
	t291 = pkin(1) * pkin(2);
	t246 = 0.1e1 / pkin(4);
	t296 = 0.1e1 / t239 ^ 2 * t246;
	t270 = t291 * t296;
	t265 = qJD(1) * t270;
	t243 = sin(qJ(1));
	t282 = t247 * t243 ^ 2;
	t273 = qJD(1) * t282;
	t299 = 0.4e1 * t236 * t248 * t273 * t296 - t244 * t265;
	t290 = -pkin(3) - pkin(4);
	t233 = (pkin(2) - t290) * (pkin(2) + t290) + t279;
	t289 = -pkin(3) + pkin(4);
	t234 = (pkin(2) - t289) * (pkin(2) + t289) + t279;
	t285 = t233 * t234;
	t249 = sqrt(-t285);
	t281 = t249 * t243;
	t235 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t239;
	t284 = t235 * t244;
	t298 = (t281 + t284) * pkin(2);
	t268 = (-t233 - t234) * t291;
	t227 = t243 * t268;
	t231 = 0.1e1 / t249;
	t297 = t231 * t227;
	t225 = (t244 * t268 - 0.4e1 * t248 * t282) * qJD(1);
	t240 = -pkin(1) + t287;
	t226 = qJD(1) * t227;
	t264 = 0.1e1 / t285 * t226 * t297;
	t295 = (t231 * t225 + t264) * t240;
	t229 = pkin(2) * t281;
	t269 = -t235 + t297;
	t277 = 0.2e1 * t240 * pkin(1);
	t259 = t269 + t277;
	t261 = t243 * t265;
	t266 = t243 * t270;
	t272 = t236 * t246 / 0.2e1;
	t278 = qJD(1) * t243;
	t283 = t244 * t249;
	t228 = pkin(2) * qJD(1) * t283;
	t288 = pkin(2) * t243;
	t280 = t231 * t226 * t288 + t228;
	t292 = 0.2e1 * t231;
	t294 = -(t259 * t243 + t283) * pkin(2) * t261 - ((-t235 + t277) * pkin(2) * t278 + t280) * t266 + (-0.4e1 * pkin(1) * t273 + (t243 * t264 + (t244 * t226 / 0.2e1 + t243 * t225 / 0.2e1) * t292 + (t259 * t244 - t281) * qJD(1)) * pkin(2)) * t272 + t299 * (t235 * t240 + t229);
	t276 = pkin(1) * t282;
	t271 = 0.2e1 * t276;
	t286 = t231 * t240;
	t293 = (pkin(2) * t284 - t227 * t286 + t229 + t271) * t261 - t299 * (t235 * t288 - t240 * t249);
	t275 = t247 * t244 * pkin(1);
	t274 = t226 * t286;
	t1 = [((-t228 + t295 + (-0.6e1 * qJD(1) * t275 + (qJD(1) * t235 - t226 * t292) * pkin(2)) * t243) * t272 - (t274 + (-0.2e1 * t276 - t298) * qJD(1)) * t266 + t293) * r_i_i_C(2) + t294 * r_i_i_C(1); t294 * r_i_i_C(2) + ((-t295 + (t269 * pkin(2) + 0.6e1 * t275) * t278 + t280) * t272 - (-t274 + (t271 + t298) * qJD(1)) * t266 - t293) * r_i_i_C(1); 0;];
	JaD_transl = t1;
end