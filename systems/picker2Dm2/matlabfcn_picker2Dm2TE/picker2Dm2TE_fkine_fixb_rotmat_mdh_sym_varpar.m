% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% picker2Dm2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% T_c_mdh [4x4x(15+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   11:  mdh base (link 0) -> mdh frame (11-1), link (11-1)
%   ...
%   15+1:  mdh base (link 0) -> mdh frame (15)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = picker2Dm2TE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 10:31:38
% EndTime: 2020-05-09 10:31:47
% DurationCPUTime: 7.51s
% Computational Cost: add. (61477->430), mult. (172690->535), div. (1582->5), fcn. (36146->10), ass. (0->238)
t130 = cos(pkin(9));
t162 = 0.1e1 / pkin(3);
t127 = sin(qJ(1));
t161 = pkin(3) ^ 2;
t166 = (pkin(1) ^ 2);
t168 = pkin(7) ^ 2;
t89 = t166 + t168;
t226 = t161 + t89;
t129 = cos(qJ(1));
t126 = sin(qJ(2));
t83 = pkin(3) * t126;
t307 = t83 + pkin(7);
t240 = t307 * t129;
t128 = cos(qJ(2));
t290 = t127 * t128;
t247 = pkin(3) * t290;
t220 = pkin(1) * t247;
t61 = -0.2e1 * t220;
t312 = 0.2e1 * pkin(7);
t71 = t83 * t312;
t192 = 0.1e1 / (0.2e1 * pkin(1) * t240 + t226 + t61 + t71);
t94 = t126 ^ 2;
t300 = t161 * t94;
t256 = 0.2e1 * t300;
t270 = -t161 + t168;
t199 = t71 + t256 + t270;
t97 = t129 ^ 2;
t194 = t199 * t97;
t210 = -pkin(1) + t247;
t145 = 3 * t166;
t156 = pkin(4) ^ 2;
t271 = t156 - t168;
t227 = -0.2e1 * t161 + t271;
t211 = t145 - t227;
t212 = -0.4e1 * t220;
t228 = -t156 + t89;
t289 = t128 * t129;
t246 = pkin(3) * t289;
t310 = 0.4e1 * t161;
t154 = 0.2e1 * pkin(3);
t164 = t166 ^ 2;
t214 = t71 + t228;
t258 = -0.4e1 * t83;
t268 = t166 - t168;
t308 = 2 * t166;
t321 = 0.4e1 * pkin(1);
t35 = sqrt(-0.4e1 * t166 * t194 + 0.4e1 * t268 * t300 + pkin(7) * t228 * t258 - t164 + t227 * t308 - (t168 - (t154 + pkin(4)) * pkin(4)) * (t168 + (t154 - pkin(4)) * pkin(4)) + (-(t61 + t214) * t240 + t214 * t247) * t321);
t188 = t192 * ((t127 * t307 + t246) * t35 - (t211 + t71 + t212) * t240 + t210 * t71 + t211 * t247 + (-0.2e1 * t194 + t256 - t310 - t228) * pkin(1));
t186 = t188 / 0.2e1;
t311 = 0.2e1 * t161;
t202 = t311 + t214;
t85 = pkin(1) * t129;
t261 = 0.2e1 * t85;
t84 = pkin(3) * t128;
t191 = t192 * ((t240 - t210) * t35 + (t199 * t261 + t202 * t307) * t127 + (t129 * t202 + (0.4e1 * t97 - 0.2e1) * pkin(1) * t307) * t84);
t190 = -t191 / 0.2e1;
t124 = cos(pkin(8));
t301 = sin(pkin(8));
t49 = t124 * t127 - t129 * t301;
t50 = t124 * t129 + t127 * t301;
t184 = t162 * (t186 * t50 + t190 * t49);
t189 = t191 / 0.2e1;
t26 = (t186 * t49 + t189 * t50) * t162;
t306 = sin(pkin(9));
t22 = t130 * t26 - t184 * t306;
t23 = t130 * t184 + t26 * t306;
t187 = -t188 / 0.2e1;
t31 = (t127 * t187 + t129 * t190) * t162;
t32 = (t127 * t189 + t129 * t187) * t162;
t14 = t22 * t32 - t23 * t31;
t224 = -t22 * t31 - t32 * t23;
t326 = -t14 * t22 + t224 * t23;
t2 = t14 * t23 + t22 * t224;
t106 = -t156 / 0.6e1;
t107 = -t156 / 0.4e1;
t108 = -t156 / 0.3e1;
t109 = -t156 / 0.2e1;
t115 = 0.2e1 / 0.3e1 * t161;
t118 = -t161 / 0.3e1;
t120 = 0.4e1 / 0.3e1 * t166;
t122 = t166 / 0.2e1;
t131 = 15 * t164;
t132 = 15 * t166;
t133 = 10 * t166;
t138 = -0.2e1 * t156;
t139 = -0.5e1 * t156;
t173 = t161 ^ 2;
t140 = 0.5e1 * t173;
t141 = 7 * t164;
t142 = 5 * t164;
t143 = 7 * t166;
t144 = 6 * t166;
t167 = t168 ^ 2;
t148 = 0.3e1 * t167;
t149 = 0.8e1 * t168;
t150 = 0.4e1 * t168;
t152 = 0.2e1 * t168;
t155 = t156 ^ 2;
t172 = pkin(3) * t161;
t158 = t172 ^ 2;
t169 = pkin(1) * t166;
t177 = pkin(7) * t168;
t269 = t164 + t167;
t274 = t152 - t156;
t284 = t168 * t156;
t196 = t274 * t166 + t155 / 0.6e1 + t269 - t284;
t193 = 0.5e1 / 0.6e1 * t173 + t196;
t200 = t168 - t220;
t64 = t109 + t226;
t203 = t64 * t212;
t213 = -0.6e1 * t220;
t282 = t155 / 0.2e1 - t173 / 0.2e1;
t216 = t148 - 0.3e1 * t284 + t282;
t111 = -0.3e1 / 0.2e1 * t156;
t151 = 0.3e1 * t168;
t281 = t111 + t151;
t292 = t158 + t89 * ((t111 + t152) * t166 - 0.3e1 / 0.2e1 * t284 + t269 + t282);
t81 = 0.10e2 / 0.3e1 * t166;
t206 = ((t81 + t274) * t161 + t193) * t213 + (t131 + (-0.9e1 * t156 + 0.18e2 * t168) * t166 + t216) * t161 + (t132 + t281) * t173 + t292;
t275 = t145 + t168;
t230 = t161 + t275;
t305 = pkin(1) * t127;
t208 = t230 * t84 - t305 * (0.3e1 * t161 + t89);
t110 = -0.2e1 / 0.3e1 * t156;
t119 = -0.2e1 / 0.3e1 * t161;
t231 = t110 + t89;
t276 = t133 + t152;
t280 = t119 + t168;
t87 = -t156 - t161;
t73 = t151 + t87;
t293 = t73 * t166;
t209 = -t305 * (t140 + ((5 * t166) + t73) * t311 + (t119 + t231) * t89) + (t173 + (t110 + t119 + t276) * t161 + t142 + 0.2e1 * t293 + t168 * (t110 + t280)) * t84;
t273 = t155 - t173;
t215 = 0.6e1 * t167 - 0.6e1 * t284 + t273;
t232 = t110 + t115 + t152;
t291 = t173 + (t115 + t231) * t89;
t114 = 0.4e1 / 0.3e1 * t161;
t233 = t108 + t89;
t62 = t114 + t233;
t217 = t62 * t212 + t291 + (t144 + t232) * t161;
t96 = t129 * t97;
t296 = t169 * t96;
t218 = t296 * t84;
t299 = t164 * t97 ^ 2;
t219 = t299 * t84;
t223 = 0.20e2 / 0.3e1 * t166;
t249 = 0.16e2 * t296;
t225 = pkin(7) * t249;
t272 = -t156 + t161;
t229 = t151 + t272;
t116 = t161 / 0.3e1;
t234 = t106 + t116 + t168;
t235 = t156 / 0.3e1 + t116 + t152;
t236 = 0.2e1 / 0.3e1 * t156 + t115 + t150;
t237 = 0.4e1 / 0.3e1 * t156 + t114 - 0.2e1 * t168;
t257 = pkin(7) * t296;
t238 = 0.8e1 * t257;
t288 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t239 = t127 * t288;
t259 = 0.6e1 * t85;
t241 = pkin(7) * t259;
t260 = 0.4e1 * t85;
t242 = pkin(7) * t260;
t244 = -t305 / 0.2e1;
t245 = t166 * t84;
t298 = t166 * t97;
t250 = 0.12e2 * t298;
t295 = t172 * t126 * t94;
t253 = -0.8e1 * t295;
t254 = 0.4e1 * t298;
t255 = 0.8e1 * t299;
t262 = 0.2e1 * t305;
t263 = pkin(7) * t85;
t265 = 0.4e1 * pkin(7);
t266 = t167 + t173;
t267 = t167 - t164;
t277 = 0.4e1 / 0.7e1 * t168 - t156 / 0.7e1;
t278 = t122 + t168;
t279 = t166 / 0.3e1 + t168;
t283 = t168 * t166;
t287 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t294 = t173 * t94 ^ 2;
t297 = t168 * t87;
t82 = t89 ^ 2;
t302 = t82 * (-t161 + t228);
t303 = (-t127 * t169 + t245) * t97;
t309 = 4 * t164;
t313 = t107 + t161 / 0.2e1;
t48 = -t173 / 0.6e1 + t196;
t77 = t118 + t168;
t51 = t77 * t61;
t70 = t85 + pkin(7);
t56 = t83 + t70;
t88 = -0.3e1 * t161 + t168;
t59 = t88 * t238;
t60 = 0.10e2 * t293;
t65 = -t156 + t226;
t72 = pkin(7) * t261;
t75 = (t150 + t156) * t166;
t78 = -t166 / 0.3e1 + t168;
t86 = -0.30e2 * t156 + 0.60e2 * t168;
t91 = -(3 * t166) + t168;
t304 = ((-0.24e2 * (0.4e1 / 0.3e1 * t298 + t72 + t78) * t294 * t305 - 0.12e2 * (-0.8e1 / 0.3e1 * t219 + ((t120 + t234) * t84 - (0.7e1 / 0.6e1 * t161 + t106 + t278) * t305) * t254 + (-t161 * t268 - 0.5e1 / 0.3e1 * t164 + t235 * t166 + t168 * (t108 + t77)) * t84 + (-t173 + (-t223 + t236) * t161 - (3 * t164) + t237 * t166 + t167) * t244 + (-t127 * t164 * t96 + ((t166 + t234) * t84 + (t311 - t268) * t244) * t85) * t265) * t300 + 0.24e2 * t77 * t219 + ((t168 + 0.5e1 / 0.2e1 * t161 + 0.3e1 / 0.2e1 * t166 + t109) * t84 + t88 * t305 / 0.2e1) * t225 - 0.6e1 * ((-0.3e1 * t173 + (-t223 + t237) * t161 + t236 * t166 + t267) * t84 - 0.2e1 * (-0.5e1 / 0.3e1 * t173 + (-t166 + t235) * t161 + t168 * (t118 + t233)) * t305) * t298 - 0.6e1 * t209 * t263 - (t158 + ((21 * t166) + t73) * t173 + (t148 + (35 * t164) + t60 + 0.2e1 * t297) * t161 + (t141 + (t139 + t149 - 0.5e1 * t161) * t166 + t168 * (-t156 + t270)) * t89) * t84 + (0.7e1 * t158 + (t143 + t73) * t140 + ((21 * t164) + 0.9e1 * t167 + t60 + 0.6e1 * t297) * t161 + t302) * t305) * t35 + (0.16e2 * (t255 + t225 + (-(8 * t164) + 0.12e2 * t283) * t97 + (-0.12e2 * pkin(7) * t169 + t177 * t321) * t129 - 0.6e1 * t283 + t269) * t294 + 0.24e2 * (t280 * t255 + 0.14e2 * (-0.32e2 / 0.21e2 * (t168 + t161 / 0.4e1 + t166 / 0.4e1 - t156 / 0.8e1) * t220 + 0.5e1 / 0.42e2 * t173 + (0.16e2 / 0.21e2 * t166 + t277) * t161 + t164 / 0.7e1 + t277 * t166 + t167 - 0.3e1 / 0.7e1 * t284 + t155 / 0.42e2) * t298 + t78 * t203 - t268 * t173 + (-0.10e2 / 0.3e1 * t164 + 0.2e1 * t167 - t284 + t75) * t161 + t48 * t287 + ((-0.2e1 / 0.3e1 * t220 + t107 + t278) * t249 + (-0.8e1 / 0.3e1 * (t279 + t313) * t220 + 0.5e1 / 0.18e2 * t173 + (0.4e1 / 0.3e1 * t168 + t120 + t108) * t161 + t167 + 0.2e1 / 0.3e1 * t283 - 0.2e1 / 0.3e1 * t284 - t164 / 0.3e1 + t155 / 0.18e2) * t259) * pkin(7)) * t300 + 0.16e2 * (-0.6e1 * t168 * t161 + t266) * t299 + 0.32e2 * (t288 * t61 + t64 * t88) * t257 + 0.24e2 * (t77 * t203 - t158 + (-t81 + t271) * t173 + (t75 + t173 / 0.6e1 - t155 / 0.6e1 + t267) * t161 + t48 * t168) * t298 + 0.8e1 * t206 * t263 - 0.8e1 * ((t143 + t281) * t173 + (t141 + (t139 + 0.10e2 * t168) * t166 + t216) * t161 + t292) * t220 + t173 ^ 2 + (t138 + t150 + (28 * t166)) * t158 + (t166 * t86 + (70 * t164) + t215) * t173 + (t215 * t144 + t273 * t152 - 0.6e1 * t167 * t156 + t86 * t164 + 0.28e2 * t169 ^ 2 + 0.4e1 * t177 ^ 2) * t161 + t65 * t302) * t56 + (((0.4e1 * t303 + (t84 + t262) * t72 + t91 * t84 + (t109 + t230) * t262) * t253 - 0.6e1 * (-0.4e1 * (0.2e1 * (0.5e1 / 0.6e1 * t161 + t122 + t106) * t84 + pkin(1) * t239) * t298 + (-0.8e1 * t218 + ((t108 + t115 + t275) * t84 - (0.8e1 / 0.3e1 * t161 + t233) * t305) * t260) * pkin(7) + t209) * t83) * t35 + (0.32e2 * (t238 + (-0.4e1 * t169 * t247 + t309 + (t310 + t138 + t149) * t166) * t97 + (-t166 + t200 + t313) * t242 + t61 * t287 + t91 * t64) * t295 + 0.8e1 * (t59 + (t288 * t64 + t51) * t250 + (t203 + (t144 + t274) * t161 + t193) * t241 + t206) * t83) * t56) * t70) / ((-0.4e1 * (-t268 * t84 + 0.2e1 * t303 + (t246 * t312 + t127 * (t161 + t308)) * pkin(1)) * t300 + 0.8e1 * pkin(7) * t218 + ((pkin(3) * t309 + 0.8e1 * t166 * t172) * t128 + 0.4e1 * t169 * t239) * t97 - 0.4e1 * t208 * t263 - (t161 * t276 + t142 + t266 + 0.6e1 * t283) * t84 + (t140 + (t133 + 0.6e1 * t168) * t161 + t82) * t305) * t35 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t220 + 0.4e1 / 0.9e1 * t161 - t156 / 0.9e1 + t279) * t298 + t78 * t61 + t62 * t287 + (t296 + (t106 + t115 + t200) * t85) * t265) * t300 + t59 + (t288 * t62 + t51) * t250 + t217 * t241 + ((t81 + t232) * t161 + t291) * t213 + t158 + (t132 + t229) * t173 + (t229 * t144 + t272 * t152 + t131 + t148) * t161 + t82 * t65) * t56 + ((t253 * t305 + (-0.2e1 * t97 * t245 + (t84 - t305) * t72 + t208) * t258) * t35 + (0.8e1 * (t72 + t254 + t91) * t295 + 0.6e1 * (t270 * t254 + (t61 + t62) * t242 + t217) * t83) * t56) * t70);
t243 = t304 / 0.2e1;
t157 = 0.1e1 / pkin(4);
t286 = t157 * t162;
t197 = (-t31 * t35 / 0.2e1 + t32 * t243) * t286;
t314 = (t31 * t243 + t32 * t35 / 0.2e1) * t286;
t185 = t188 / 0.4e1;
t285 = t157 / pkin(3) ^ 2;
t17 = (t185 * t304 - t35 * t191 / 0.4e1) * t285;
t18 = (t35 * t185 + t191 * t304 / 0.4e1) * t285;
t52 = -t126 * t127 - t289;
t53 = -t126 * t129 + t290;
t5 = t17 * t52 + t18 * t53;
t6 = t17 * t53 - t18 * t52;
t325 = t197 * t6 + t314 * t5;
t324 = -t197 * t5 + t314 * t6;
t201 = t127 * t49 + t129 * t50;
t41 = t127 * t50 - t129 * t49;
t323 = -t201 * t50 + t41 * t49;
t322 = t201 * t49 + t41 * t50;
t264 = pkin(7) + 0;
t252 = t301 * pkin(5) + 0;
t251 = t124 * pkin(5) + 0;
t68 = -t85 + 0;
t222 = t32 * pkin(3) + t68;
t221 = t32 * pkin(2) + t68;
t67 = 0 - t305;
t12 = t224 * pkin(6);
t207 = t12 + t221;
t205 = t31 * pkin(3) + t67;
t204 = t31 * pkin(2) + t67;
t11 = t14 * pkin(6);
t198 = t11 + t204;
t39 = -t124 * t50 + t301 * t49;
t38 = -t124 * t49 - t301 * t50;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t129, t127, 0, 0; -t127, -t129, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t31, 0, t68; t31, t32, 0, t67; 0, 0, 1, 0; 0, 0, 0, 1; t224, -t14, 0, t221; t14, t224, 0, t204; 0, 0, 1, 0; 0, 0, 0, 1; t197, -t314, 0, t222; t314, t197, 0, t205; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t38, 0, t251; t38, t39, 0, t252; 0, 0, 1, 0; 0, 0, 0, 1; t224, -t14, 0, t68; t14, t224, 0, t67; 0, 0, 1, 0; 0, 0, 0, 1; t126, t128, 0, t264; -t128, t126, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t201, -t41, 0, t68; t41, t201, 0, t67; 0, 0, 1, 0; 0, 0, 0, 1; t326, -t2, 0, t207; t2, t326, 0, t198; 0, 0, 1, 0; 0, 0, 0, 1; t325, -t324, 0, pkin(4) * t197 + t222; t324, t325, 0, pkin(4) * t314 + t205; 0, 0, 1, 0; 0, 0, 0, 1; t323, t322, 0, pkin(5) * t201 + t68; -t322, t323, 0, pkin(5) * t41 + t67; 0, 0, 1, 0; 0, 0, 0, 1; t326, -t2, 0, t12 + t68; t2, t326, 0, t11 + t67; 0, 0, 1, 0; 0, 0, 0, 1; t126, t128, 0, t83 + t264; -t128, t126, 0, -t84 + 0; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t38, 0, pkin(1) * t39 + t251; t38, t39, 0, pkin(1) * t38 + t252; 0, 0, 1, 0; 0, 0, 0, 1; t326, -t2, 0, pkin(2) * t326 + t207; t2, t326, 0, pkin(2) * t2 + t198; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,15+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,15+1]); end % symbolisch
for i = 1:15+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end