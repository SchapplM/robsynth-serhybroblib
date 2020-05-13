% Calculate potential energy for
% picker2Dm1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 20:21:18
% EndTime: 2020-05-10 20:21:33
% DurationCPUTime: 7.87s
% Computational Cost: add. (65214->507), mult. (176706->627), div. (2626->12), fcn. (45644->38), ass. (0->267)
t329 = 4 * pkin(1);
t306 = 0.1e1 / pkin(2) / 0.2e1;
t168 = 0.1e1 / t306;
t175 = pkin(4) ^ 2;
t117 = -t175 / 0.4e1;
t180 = pkin(3) ^ 2;
t328 = t117 + t180 / 0.2e1;
t327 = 2 * pkin(7);
t134 = sin(pkin(8));
t135 = cos(pkin(8));
t138 = sin(qJ(1));
t141 = cos(qJ(1));
t61 = t134 * t141 - t135 * t138;
t326 = 0.2e1 * t61;
t107 = t141 ^ 2;
t325 = -0.2e1 * t107;
t153 = 0.2e1 * t180;
t324 = 0.4e1 * t180;
t185 = pkin(1) ^ 2;
t183 = t185 ^ 2;
t323 = 4 * t183;
t322 = 2 * t185;
t157 = 6 * t185;
t188 = pkin(7) ^ 2;
t165 = 2 * t188;
t193 = t180 ^ 2;
t152 = 0.5e1 * t193;
t169 = 2 * pkin(1);
t172 = pkin(5) ^ 2;
t62 = -t134 * t138 - t135 * t141;
t320 = pkin(5) * t62;
t270 = pkin(1) * t320;
t304 = t172 - 0.2e1 * t270;
t46 = sqrt(-(-(t169 + pkin(5)) * pkin(5) + t304) * (pkin(5) * (t169 - pkin(5)) + t304));
t315 = t46 * t61;
t53 = t172 - t270;
t57 = -pkin(1) + t320;
t44 = -pkin(5) * t315 - 0.2e1 * t53 * t57;
t321 = t44 / 0.4e1;
t95 = pkin(1) * t141;
t82 = t95 + pkin(7);
t139 = sin(pkin(9));
t142 = cos(pkin(9));
t181 = 0.1e1 / pkin(3);
t99 = t185 + t188;
t232 = t180 + t99;
t268 = 0.2e1 * t95;
t140 = cos(qJ(2));
t298 = t138 * t140;
t257 = pkin(3) * t298;
t227 = pkin(1) * t257;
t75 = -0.2e1 * t227;
t137 = sin(qJ(2));
t93 = pkin(3) * t137;
t81 = t93 + pkin(7);
t83 = t93 * t327;
t49 = 0.1e1 / (t268 * t81 + t232 + t75 + t83);
t311 = t181 * t49;
t52 = 0.1e1 / (t185 + t304);
t312 = 0.1e1 / pkin(5) * t52;
t223 = t311 * t312;
t215 = -pkin(1) + t257;
t158 = 3 * t185;
t278 = t175 - t188;
t216 = t153 + t158 - t278;
t219 = -0.4e1 * t227;
t233 = -t175 + t99;
t297 = t140 * t141;
t256 = pkin(3) * t297;
t104 = t137 ^ 2;
t299 = t104 * t180;
t261 = 0.2e1 * t299;
t313 = t141 * t81;
t167 = 0.2e1 * pkin(3);
t222 = t83 + t233;
t293 = t185 * t107;
t259 = -0.4e1 * t293;
t265 = -0.4e1 * t93;
t275 = t185 - t188;
t277 = -t180 + t188;
t63 = t83 + t261 + t277;
t43 = sqrt(t63 * t259 + 0.4e1 * t275 * t299 + pkin(7) * t233 * t265 - t183 + (-0.2e1 * t180 + t278) * t322 - (t188 - (t167 + pkin(4)) * pkin(4)) * (t188 + (t167 - pkin(4)) * pkin(4)) + (-(t75 + t222) * t313 + t222 * t257) * t329);
t41 = (t138 * t81 + t256) * t43 - (t216 + t83 + t219) * t313 + t215 * t83 + t216 * t257 + (t325 * t63 - t233 + t261 - t324) * pkin(1);
t68 = t153 + t222;
t94 = pkin(3) * t140;
t42 = (-t215 + t313) * t43 + (t268 * t63 + t68 * t81) * t138 + (t141 * t68 + (0.4e1 * t107 - 0.2e1) * t81 * pkin(1)) * t94;
t45 = pkin(5) * t326 * t53 - t46 * t57;
t28 = (t41 * t321 + t42 * t45 / 0.4e1) * t223;
t29 = (-t41 * t45 / 0.4e1 + t42 * t321) * t223;
t26 = t139 * t29 + t142 * t28;
t319 = pkin(6) * t26;
t27 = t139 * t28 - t142 * t29;
t318 = pkin(6) * t27;
t317 = pkin(1) * t138;
t101 = -3 * t185 + t188;
t106 = t141 * t107;
t116 = -t175 / 0.6e1;
t118 = -t175 / 0.3e1;
t119 = -t175 / 0.2e1;
t125 = 0.2e1 / 0.3e1 * t180;
t128 = -t180 / 0.3e1;
t130 = 0.4e1 / 0.3e1 * t185;
t132 = t185 / 0.2e1;
t143 = 15 * t183;
t144 = 15 * t185;
t145 = 10 * t185;
t150 = -0.2e1 * t175;
t151 = -0.5e1 * t175;
t154 = 7 * t183;
t155 = 5 * t183;
t156 = 7 * t185;
t187 = t188 ^ 2;
t161 = 3 * t187;
t162 = 8 * t188;
t163 = 4 * t188;
t174 = t175 ^ 2;
t192 = pkin(3) * t180;
t177 = t192 ^ 2;
t189 = pkin(1) * t185;
t197 = pkin(7) * t188;
t276 = t183 + t187;
t281 = t165 - t175;
t292 = t188 * t175;
t205 = t281 * t185 + t174 / 0.6e1 + t276 - t292;
t204 = 0.5e1 / 0.6e1 * t193 + t205;
t208 = t188 - t227;
t78 = t119 + t232;
t209 = t78 * t219;
t289 = t174 / 0.2e1 - t193 / 0.2e1;
t218 = -0.3e1 * t292 + t161 + t289;
t220 = -0.6e1 * t227;
t121 = -0.3e1 / 0.2e1 * t175;
t164 = 3 * t188;
t288 = t121 + t164;
t303 = t177 + t99 * ((t121 + t165) * t185 - 0.3e1 / 0.2e1 * t292 + t276 + t289);
t91 = 0.10e2 / 0.3e1 * t185;
t212 = ((t91 + t281) * t180 + t204) * t220 + (t143 + (-0.9e1 * t175 + (18 * t188)) * t185 + t218) * t180 + (t144 + t288) * t193 + t303;
t282 = t158 + t188;
t235 = t180 + t282;
t213 = t235 * t94 - t317 * (0.3e1 * t180 + t99);
t120 = -0.2e1 / 0.3e1 * t175;
t129 = -0.2e1 / 0.3e1 * t180;
t236 = t120 + t99;
t283 = t145 + t165;
t287 = t129 + t188;
t97 = -t175 - t180;
t85 = t164 + t97;
t307 = t85 * t185;
t214 = -t317 * (t152 + ((5 * t185) + t85) * t153 + (t129 + t236) * t99) + (t193 + (t120 + t129 + t283) * t180 + t155 + 0.2e1 * t307 + t188 * (t120 + t287)) * t94;
t280 = t174 - t193;
t217 = -0.6e1 * t292 + (6 * t187) + t280;
t254 = t189 * t94;
t224 = t106 * t254;
t294 = t183 * t107 ^ 2;
t225 = t294 * t94;
t237 = t120 + t125 + t165;
t302 = t193 + (t125 + t236) * t99;
t124 = 0.4e1 / 0.3e1 * t180;
t238 = t118 + t99;
t76 = t124 + t238;
t226 = t76 * t219 + t302 + (t157 + t237) * t180;
t290 = t189 * t106;
t252 = 0.16e2 * t290;
t228 = pkin(7) * t252;
t229 = 0.20e2 / 0.3e1 * t185;
t263 = pkin(7) * t290;
t230 = 0.8e1 * t263;
t279 = -t175 + t180;
t234 = t164 + t279;
t126 = t180 / 0.3e1;
t239 = t116 + t126 + t188;
t240 = t175 / 0.3e1 + t126 + t165;
t241 = 0.2e1 / 0.3e1 * t175 + t125 + t163;
t242 = 0.4e1 / 0.3e1 * t175 + t124 - (2 * t188);
t296 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t246 = t138 * t296;
t266 = 0.6e1 * t95;
t247 = pkin(7) * t266;
t267 = 0.4e1 * t95;
t248 = pkin(7) * t267;
t250 = -t317 / 0.2e1;
t253 = 0.12e2 * t293;
t255 = t185 * t94;
t258 = 0.4e1 * t293;
t260 = 0.8e1 * t294;
t300 = t137 * t104 * t192;
t262 = -0.8e1 * t300;
t269 = 0.2e1 * t317;
t271 = pkin(7) * t95;
t272 = 4 * pkin(7);
t273 = t187 + t193;
t274 = t187 - t183;
t284 = 0.4e1 / 0.7e1 * t188 - t175 / 0.7e1;
t285 = t132 + t188;
t286 = t185 / 0.3e1 + t188;
t291 = t188 * t185;
t295 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t301 = t104 ^ 2 * t193;
t308 = (-t138 * t189 + t255) * t107;
t309 = t188 * t97;
t92 = t99 ^ 2;
t314 = t92 * (-t180 + t233);
t60 = -t193 / 0.6e1 + t205;
t89 = t128 + t188;
t64 = t89 * t75;
t70 = t93 + t82;
t98 = -0.3e1 * t180 + t188;
t73 = t98 * t230;
t74 = 0.10e2 * t307;
t79 = -t175 + t232;
t84 = pkin(7) * t268;
t87 = (t163 + t175) * t185;
t90 = -t185 / 0.3e1 + t188;
t96 = -0.30e2 * t175 + (60 * t188);
t316 = ((-0.24e2 * (0.4e1 / 0.3e1 * t293 + t84 + t90) * t301 * t317 - 0.12e2 * (-0.8e1 / 0.3e1 * t225 + ((t130 + t239) * t94 - (0.7e1 / 0.6e1 * t180 + t116 + t285) * t317) * t258 + (-t180 * t275 - 0.5e1 / 0.3e1 * t183 + t240 * t185 + t188 * (t118 + t89)) * t94 + (-t193 + (-t229 + t241) * t180 - (3 * t183) + t242 * t185 + t187) * t250 + (-t138 * t183 * t106 + ((t185 + t239) * t94 + (t153 - t275) * t250) * t95) * t272) * t299 + 0.24e2 * t89 * t225 + ((t188 + 0.5e1 / 0.2e1 * t180 + 0.3e1 / 0.2e1 * t185 + t119) * t94 + t98 * t317 / 0.2e1) * t228 - 0.6e1 * ((-0.3e1 * t193 + (-t229 + t242) * t180 + t241 * t185 + t274) * t94 - 0.2e1 * (-0.5e1 / 0.3e1 * t193 + (-t185 + t240) * t180 + t188 * (t128 + t238)) * t317) * t293 - 0.6e1 * t214 * t271 - (t177 + ((21 * t185) + t85) * t193 + (t161 + (35 * t183) + t74 + 0.2e1 * t309) * t180 + (t154 + (t151 + t162 - 0.5e1 * t180) * t185 + t188 * (-t175 + t277)) * t99) * t94 + (0.7e1 * t177 + (t156 + t85) * t152 + ((21 * t183) + (9 * t187) + t74 + 0.6e1 * t309) * t180 + t314) * t317) * t43 + (0.16e2 * (t260 + t228 + (-8 * t183 + 12 * t291) * t107 + (-12 * pkin(7) * t189 + t197 * t329) * t141 - (6 * t291) + t276) * t301 + 0.24e2 * (t287 * t260 + 0.14e2 * (-0.32e2 / 0.21e2 * (t188 + t180 / 0.4e1 + t185 / 0.4e1 - t175 / 0.8e1) * t227 + 0.5e1 / 0.42e2 * t193 + (0.16e2 / 0.21e2 * t185 + t284) * t180 + t183 / 0.7e1 + t284 * t185 + t187 - 0.3e1 / 0.7e1 * t292 + t174 / 0.42e2) * t293 + t90 * t209 - t275 * t193 + (-0.10e2 / 0.3e1 * t183 + (2 * t187) - t292 + t87) * t180 + t60 * t295 + ((-0.2e1 / 0.3e1 * t227 + t117 + t285) * t252 + (-0.8e1 / 0.3e1 * (t286 + t328) * t227 + 0.5e1 / 0.18e2 * t193 + (0.4e1 / 0.3e1 * t188 + t130 + t118) * t180 + t187 + 0.2e1 / 0.3e1 * t291 - 0.2e1 / 0.3e1 * t292 - t183 / 0.3e1 + t174 / 0.18e2) * t266) * pkin(7)) * t299 + 0.16e2 * (-0.6e1 * t188 * t180 + t273) * t294 + 0.32e2 * (t296 * t75 + t78 * t98) * t263 + 0.24e2 * (t89 * t209 - t177 + (-t91 + t278) * t193 + (t87 + t193 / 0.6e1 - t174 / 0.6e1 + t274) * t180 + t60 * t188) * t293 + 0.8e1 * t212 * t271 - 0.8e1 * ((t156 + t288) * t193 + (t154 + (t151 + (10 * t188)) * t185 + t218) * t180 + t303) * t227 + t193 ^ 2 + (t150 + t163 + (28 * t185)) * t177 + (t185 * t96 + (70 * t183) + t217) * t193 + (t217 * t157 + t280 * t165 - 0.6e1 * t187 * t175 + t96 * t183 + (28 * t189 ^ 2) + (4 * t197 ^ 2)) * t180 + t79 * t314) * t70 + (((0.4e1 * t308 + (t94 + t269) * t84 + t101 * t94 + (t119 + t235) * t269) * t262 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t180 + t132 + t116) * t94 + pkin(1) * t246) * t259 + (-0.8e1 * t224 + ((t118 + t125 + t282) * t94 - (0.8e1 / 0.3e1 * t180 + t238) * t317) * t267) * pkin(7) + t214) * t93) * t43 + (0.32e2 * (t230 + (-0.4e1 * t138 * t254 + t323 + (t324 + t150 + t162) * t185) * t107 + (-t185 + t208 + t328) * t248 + t75 * t295 + t101 * t78) * t300 + 0.8e1 * (t73 + (t296 * t78 + t64) * t253 + (t209 + (t157 + t281) * t180 + t204) * t247 + t212) * t93) * t70) * t82) / ((-0.4e1 * (-t275 * t94 + 0.2e1 * t308 + (t256 * t327 + t138 * (t180 + t322)) * pkin(1)) * t299 + 0.8e1 * pkin(7) * t224 + ((pkin(3) * t323 + 0.8e1 * t185 * t192) * t140 + 0.4e1 * t189 * t246) * t107 - 0.4e1 * t213 * t271 - (t180 * t283 + t155 + t273 + (6 * t291)) * t94 + (t152 + (t145 + 6 * t188) * t180 + t92) * t317) * t43 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t227 + 0.4e1 / 0.9e1 * t180 - t175 / 0.9e1 + t286) * t293 + t90 * t75 + t76 * t295 + (t290 + (t116 + t125 + t208) * t95) * t272) * t299 + t73 + (t296 * t76 + t64) * t253 + t226 * t247 + ((t91 + t237) * t180 + t302) * t220 + t177 + (t144 + t234) * t193 + (t234 * t157 + t279 * t165 + t143 + t161) * t180 + t92 * t79) * t70 + ((t262 * t317 + (t255 * t325 + (t94 - t317) * t84 + t213) * t265) * t43 + (0.8e1 * (t84 + t258 + t101) * t300 + 0.6e1 * (t277 * t258 + (t75 + t76) * t248 + t226) * t93) * t70) * t82);
t310 = 0.1e1 / pkin(1) * t52;
t243 = t311 / 0.2e1;
t33 = qJ(1) + atan2(t42 * t243, t41 * t243);
t170 = pkin(6) ^ 2;
t25 = t168 * t319;
t305 = t170 + t25;
t264 = r_base(1) - t95;
t14 = sqrt(-(-(t168 + pkin(6)) * pkin(6) + t305) * (pkin(6) * (t168 - pkin(6)) + t305));
t23 = t25 + 0.2e1 * t170;
t24 = -pkin(2) - t319;
t171 = 0.1e1 / pkin(6);
t245 = t171 / (pkin(2) ^ 2 + t305) / 0.2e1;
t6 = atan2((-t14 * t24 + t23 * t318) * t245, (-t14 * t318 - t23 * t24) * t245) + t33;
t176 = 0.1e1 / pkin(4);
t231 = t176 * t181 / 0.2e1;
t17 = atan2(t43 * t231, t231 * t316) + t33;
t251 = t176 / pkin(3) ^ 2 * t49;
t249 = t316 / 0.4e1;
t244 = t312 / 0.2e1;
t221 = r_base(2) - t317;
t32 = cos(t33);
t211 = -pkin(2) * t32 + t264;
t210 = -pkin(3) * t32 + t264;
t31 = sin(t33);
t207 = -pkin(2) * t31 + t221;
t206 = -pkin(3) * t31 + t221;
t66 = -t137 * t141 + t298;
t65 = -t137 * t138 - t297;
t58 = -pkin(1) * t62 + pkin(5);
t54 = t185 - t270;
t40 = qJ(1) + atan2(t45 * t244, t44 * t244);
t39 = pkin(8) + atan2((pkin(1) * t326 * t54 + t46 * t58) * t310 / 0.2e1, -(-pkin(1) * t315 + 0.2e1 * t54 * t58) * t310 / 0.2e1);
t38 = cos(t40);
t37 = sin(t40);
t36 = cos(t39);
t35 = sin(t39);
t19 = (t41 * t43 / 0.4e1 + t42 * t249) * t251;
t18 = (t41 * t249 - t42 * t43 / 0.4e1) * t251;
t16 = cos(t17);
t15 = sin(t17);
t13 = atan2(t27, t26) + t33;
t12 = cos(t13);
t11 = sin(t13);
t10 = atan2(t18 * t65 + t19 * t66, -t18 * t66 + t19 * t65) + t17;
t9 = cos(t10);
t8 = sin(t10);
t5 = cos(t6);
t4 = sin(t6);
t3 = atan2(t171 * t14 * t306, -t26) + t6;
t2 = cos(t3);
t1 = sin(t3);
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-rSges(2,1) * t141 + rSges(2,2) * t138 + r_base(1)) + g(2) * (-rSges(2,1) * t138 - rSges(2,2) * t141 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,1) * t32 + rSges(3,2) * t31 + t264) + g(2) * (-rSges(3,1) * t31 - rSges(3,2) * t32 + t221) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t5 - rSges(4,2) * t4 + t211) + g(2) * (rSges(4,1) * t4 + rSges(4,2) * t5 + t207) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (-rSges(5,1) * t16 + rSges(5,2) * t15 + t210) + g(2) * (-rSges(5,1) * t15 - rSges(5,2) * t16 + t206) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (pkin(5) * t135 + rSges(6,1) * t36 - rSges(6,2) * t35 + r_base(1)) + g(2) * (pkin(5) * t134 + rSges(6,1) * t35 + rSges(6,2) * t36 + r_base(2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (rSges(7,1) * t12 - rSges(7,2) * t11 + t264) + g(2) * (rSges(7,1) * t11 + rSges(7,2) * t12 + t221) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (rSges(8,1) * t137 + rSges(8,2) * t140 + pkin(7) + r_base(1)) + g(2) * (-rSges(8,1) * t140 + rSges(8,2) * t137 + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (rSges(9,1) * t38 - rSges(9,2) * t37 + t264) + g(2) * (rSges(9,1) * t37 + rSges(9,2) * t38 + t221) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (pkin(6) * t5 - rSges(10,1) * t2 + rSges(10,2) * t1 + t211) + g(2) * (pkin(6) * t4 - rSges(10,1) * t1 - rSges(10,2) * t2 + t207) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (-pkin(4) * t16 + rSges(11,1) * t9 - rSges(11,2) * t8 + t210) + g(2) * (-pkin(4) * t15 + rSges(11,1) * t8 + rSges(11,2) * t9 + t206) + g(3) * (r_base(3) + rSges(11,3)));
U = t7;
