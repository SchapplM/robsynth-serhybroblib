% Calculate inertial parameters regressor of joint inertia matrix for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh3m2TE_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_inertiaJ_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t237 = cos(pkin(16));
t229 = t237 ^ 2;
t254 = pkin(10) ^ 2;
t255 = pkin(8) ^ 2;
t228 = -t254 + t255;
t250 = cos(pkin(15));
t234 = t250 ^ 2;
t203 = t228 * t234;
t245 = sin(pkin(15));
t284 = t250 * t245;
t279 = pkin(10) * t284;
t269 = pkin(8) * t279;
t316 = -0.4e1 * t269 + 0.2e1 * t203 - t228;
t318 = t316 * t229;
t276 = t229 * t284;
t317 = -0.4e1 * t276 + 0.2e1 * t284;
t314 = 2 * pkin(6);
t313 = 2 * pkin(12);
t248 = cos(qJ(3));
t213 = pkin(4) * t248 - pkin(1);
t249 = cos(qJ(2));
t195 = t213 * t249;
t243 = sin(qJ(3));
t244 = sin(qJ(2));
t287 = t243 * t244;
t204 = pkin(4) * t287;
t178 = t204 - t195;
t176 = pkin(12) + t178;
t312 = -0.2e1 * t176;
t311 = 0.2e1 * t176;
t210 = pkin(1) * t249 + pkin(12);
t310 = 0.2e1 * t210;
t235 = qJ(3) + qJ(2);
t225 = sin(t235);
t309 = 0.2e1 * t225;
t247 = cos(qJ(4));
t308 = -0.2e1 * t247;
t307 = -t234 / 0.2e1;
t233 = t249 ^ 2;
t306 = 0.4e1 * t233 - 0.2e1;
t305 = -0.2e1 * t234 + 0.1e1;
t304 = 0.4e1 * t234 - 0.2e1;
t303 = 0.8e1 * t234 - 0.4e1;
t302 = pkin(1) * t244;
t296 = t248 * pkin(1);
t301 = pkin(4) * (-pkin(4) + t296);
t242 = sin(qJ(4));
t300 = pkin(4) * t242;
t299 = pkin(4) * t247;
t298 = pkin(8) * t234;
t297 = pkin(10) * t234;
t236 = sin(pkin(16));
t183 = -t236 * t250 - t237 * t245;
t230 = pkin(17) + pkin(18);
t217 = sin(t230);
t295 = t183 * t217;
t218 = cos(t230);
t294 = t217 * t218;
t232 = t248 ^ 2;
t256 = pkin(4) ^ 2;
t293 = t232 * t256;
t292 = t236 * t237;
t238 = sin(pkin(18));
t240 = cos(pkin(18));
t291 = t238 * t240;
t290 = t238 * t250;
t239 = sin(pkin(17));
t289 = t239 * t210;
t288 = t240 * t245;
t286 = t243 * t249;
t285 = t244 * t249;
t283 = -0.2e1 * t269 + t203;
t220 = -0.2e1 * t296;
t258 = pkin(1) ^ 2;
t282 = pkin(4) * t220 + t258;
t216 = t234 - 0.1e1 / 0.2e1;
t281 = t217 * t311;
t280 = pkin(8) * t284;
t278 = pkin(3) * t290;
t277 = 0.2e1 * t298 - 0.2e1 * t279 - pkin(8);
t270 = t216 * pkin(10);
t271 = t255 / 0.2e1 - t254 / 0.2e1;
t275 = (pkin(8) * t270 + t271 * t284) * t292;
t274 = (t216 * pkin(8) - t279) * t292;
t273 = (t270 + t280) * t292;
t272 = 0.4e1 * t297 + 0.4e1 * t280 - 0.2e1 * pkin(10);
t184 = -t236 * t245 + t237 * t250;
t268 = t184 * t218 + t295;
t267 = t272 * t229 + 0.4e1 * t274 - 0.2e1 * t280 - 0.2e1 * t297;
t188 = -t248 * t249 + t287;
t189 = t248 * t244 + t286;
t265 = t188 * t250 + t245 * t189;
t246 = sin(pkin(14));
t251 = cos(pkin(14));
t190 = t251 * t245 - t250 * t246;
t191 = t246 * t245 + t251 * t250;
t264 = t244 * t190 - t191 * t249;
t192 = t245 * pkin(8) + t250 * pkin(10);
t193 = t250 * pkin(8) - t245 * pkin(10);
t263 = t236 * t192 - t193 * t237;
t186 = t288 + t290;
t185 = t238 * t245 - t240 * t250;
t262 = pkin(4) * t286 + t213 * t244;
t261 = (-t303 * t229 + t304) * t217 * pkin(8);
t252 = pkin(12) ^ 2;
t260 = -0.2e1 * (t204 + pkin(12)) * t195 + (-t256 + t282 + 0.2e1 * t293) * t233 + t204 * t313 + t252 + t256 - t293;
t259 = -0.1e1 / 0.4e1 + t284 * t292 - t216 * t229;
t241 = cos(pkin(17));
t231 = t240 ^ 2;
t227 = t243 * pkin(1);
t226 = cos(t235);
t224 = t234 / 0.2e1;
t219 = 0.2e1 * t227;
t215 = t218 ^ 2;
t206 = t210 ^ 2;
t205 = 0.2e1 * t285;
t196 = t256 + t282;
t171 = t192 * t237 + t236 * t193;
t169 = -t245 * t188 + t189 * t250;
t168 = t249 * t190 + t244 * t191;
t167 = -t218 * t183 + t217 * t184;
t165 = t167 ^ 2;
t164 = t268 ^ 2;
t163 = t178 * t250 + t245 * t262;
t162 = t178 * t245 - t250 * t262;
t161 = (t185 * t241 + t186 * t239) * pkin(3) + t210;
t160 = t171 * t217 + t263 * t218 + t176;
t159 = (t216 * t292 + (t229 - 0.1e1 / 0.2e1) * t284) * t215 - (t224 + t259) * t294 - t276 / 0.2e1 + (t307 + 0.1e1 / 0.4e1) * t292 + t284 / 0.4e1;
t158 = (t169 * t237 - t265 * t236) * t218 - (t169 * t236 + t265 * t237) * t217;
t157 = (t162 * t237 + t163 * t236) * t218 + (-t236 * t162 + t163 * t237) * t217;
t156 = ((-t304 * pkin(8) + 0.4e1 * t279) * t229 + 0.4e1 * t273 + t277) * t215 + ((pkin(10) + t267) * t217 + t184 * t176) * t218 + t176 * t295 + t277 * t229 - 0.2e1 * t273 - t298 + t279;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t244 ^ 2, t205, 0, t233, 0, 0, t249 * t313, -0.2e1 * pkin(12) * t244, 0, t252, t189 ^ 2, (0.4e1 * t232 - 0.2e1) * t285 + t306 * t248 * t243, 0, t188 ^ 2, 0, 0, t188 * t310, t189 * t310, 0, t206, t165, ((-0.8e1 * t229 + 0.4e1) * t284 - t303 * t292) * t215 - 0.8e1 * (t307 - t259) * t294 + t304 * t292 - t317, 0, t164, 0, 0, -t268 * t311, t167 * t312, 0, t260, t247 ^ 2 * t165, t165 * t242 * t308, 0.8e1 * t247 * t159, t242 ^ 2 * t165, -0.8e1 * t242 * t159, t164, t156 * t308, 0.2e1 * t156 * t242, ((-t303 * pkin(10) - 0.8e1 * t280) * t229 - 0.8e1 * t274 + t272) * t215 + ((0.8e1 * pkin(10) * t276 + 0.8e1 * t273 - 0.4e1 * t279) * t217 + t183 * t312 + t261) * t218 + t184 * t281 + 0.2e1 * pkin(10) + t267, (-0.8e1 * t275 - t316 + 0.2e1 * t318) * t215 + ((-0.4e1 * (-t271 + t283) * t292 + t317 * t228) * t217 + t263 * t311 + pkin(10) * t261) * t218 + t171 * t281 - t318 + 0.4e1 * t275 + t254 + t260 + t283, t168 ^ 2, -0.4e1 * t234 * t285 - t306 * t284 + t205 + (-0.8e1 * ((t233 - 0.1e1 / 0.2e1) * t234 - t284 * t285 - t233 / 0.2e1 + 0.1e1 / 0.4e1) * t246 + ((0.8e1 * t233 - 0.4e1) * t284 + t303 * t285) * t251) * t251, 0, t264 ^ 2, 0, 0, t264 * t314, t168 * t314, 0, pkin(6) ^ 2, t186 ^ 2, (-0.4e1 * t231 + 0.2e1) * t284 - t304 * t291, 0, t185 ^ 2, 0, 0, t185 * t310, -0.2e1 * t210 * t186, 0, t206, t225 ^ 2, t226 * t309, 0, t226 ^ 2, 0, 0, -0.2e1 * t161 * t226, t161 * t309, 0, 0.2e1 * t278 * t289 + t206 + 0.2e1 * ((t210 * t185 + (t305 * t291 + (-0.2e1 * t231 + 0.1e1) * t284) * t239 * pkin(3)) * t241 + (t278 + t289) * t288) * pkin(3) + (-0.4e1 * (-t216 * t231 + t284 * t291 + t224 - 0.1e1 / 0.4e1) * t241 ^ 2 + t305 * t231 + t234) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, 0, t249, 0, 0, 0, 0, 0, 0, 0, -t189, 0, t188, 0, 0, 0, -t302, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, 0, 0, 0, 0, 0, 0, 0, t242 * t262, t262 * t247, 0, 0, 0, 0, t168, 0, -t264, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t302, 0, 0, 0, -t225, 0, -t226, 0, 0, 0, -t302, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t220, t219, 0, t258, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, 0, 0, 0, 0, 0, 1, t220, t219, 0, t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, t188, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pkin(4) * t189, 0, 0, 0, 0, 0, 0, 0, t189 * t300, t189 * t299, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, 0, -t226, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t296, t227, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t296, t227, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167 * t247, 0, -t167 * t242, t268, -t247 * t160, t160 * t242, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157 * t242, -t247 * t157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158 * t300, t158 * t299, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MM_reg = t1;
