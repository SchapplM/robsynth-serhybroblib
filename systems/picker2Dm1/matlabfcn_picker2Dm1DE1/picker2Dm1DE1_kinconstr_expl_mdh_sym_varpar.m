% Explicit kinematic constraints of
% picker2Dm1DE1
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
% jv [15x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 19:54
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = picker2Dm1DE1_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE1_kinconstr_expl_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE1_kinconstr_expl_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 09:01:41
% EndTime: 2020-05-10 09:01:45
% DurationCPUTime: 3.95s
% Computational Cost: add. (12598->405), mult. (34132->552), div. (521->11), fcn. (8960->24), ass. (0->244)
t301 = 4 * pkin(1);
t155 = pkin(3) ^ 2;
t150 = pkin(4) ^ 2;
t92 = -t150 / 0.4e1;
t300 = t155 / 0.2e1 + t92;
t299 = 2 * pkin(7);
t143 = 2 * pkin(2);
t109 = sin(pkin(8));
t110 = cos(pkin(8));
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t36 = t109 * t116 - t110 * t113;
t298 = 0.2e1 * t36;
t82 = t116 ^ 2;
t297 = -0.2e1 * t82;
t128 = 0.2e1 * t155;
t296 = 0.4e1 * t155;
t162 = pkin(1) ^ 2;
t160 = t162 ^ 2;
t295 = 4 * t160;
t294 = 2 * t162;
t132 = 6 * t162;
t165 = pkin(7) ^ 2;
t140 = 2 * t165;
t170 = t155 ^ 2;
t127 = 0.5e1 * t170;
t147 = pkin(5) ^ 2;
t37 = -t109 * t113 - t110 * t116;
t292 = pkin(5) * t37;
t241 = pkin(1) * t292;
t28 = t147 - t241;
t144 = 2 * pkin(1);
t267 = t147 - 0.2e1 * t241;
t21 = sqrt(-(-(t144 + pkin(5)) * pkin(5) + t267) * (pkin(5) * (t144 - pkin(5)) + t267));
t285 = t21 * t36;
t32 = -pkin(1) + t292;
t19 = -pkin(5) * t285 - 0.2e1 * t28 * t32;
t293 = t19 / 0.4e1;
t70 = pkin(1) * t116;
t57 = t70 + pkin(7);
t114 = sin(pkin(9));
t117 = cos(pkin(9));
t142 = 2 * pkin(3);
t74 = t162 + t165;
t204 = -t150 + t74;
t112 = sin(qJ(2));
t68 = pkin(3) * t112;
t58 = t68 * t299;
t192 = t58 + t204;
t115 = cos(qJ(2));
t263 = t113 * t115;
t227 = pkin(3) * t263;
t276 = t162 * t82;
t232 = -0.4e1 * t276;
t236 = -0.4e1 * t68;
t246 = t162 - t165;
t249 = t150 - t165;
t79 = t112 ^ 2;
t279 = t155 * t79;
t56 = t68 + pkin(7);
t281 = t116 * t56;
t234 = 0.2e1 * t279;
t248 = -t155 + t165;
t38 = t58 + t234 + t248;
t199 = pkin(1) * t227;
t50 = -0.2e1 * t199;
t18 = sqrt(t38 * t232 + 0.4e1 * t246 * t279 + pkin(7) * t204 * t236 - t160 + (-0.2e1 * t155 + t249) * t294 - ((t165 - (t142 + pkin(4)) * pkin(4)) * (t165 + (t142 - pkin(4)) * pkin(4))) + (-(t50 + t192) * t281 + t192 * t227) * t301);
t188 = -pkin(1) + t227;
t133 = 3 * t162;
t189 = t128 + t133 - t249;
t190 = -0.4e1 * t199;
t262 = t115 * t116;
t226 = pkin(3) * t262;
t16 = (t113 * t56 + t226) * t18 - (t189 + t58 + t190) * t281 + t188 * t58 + t189 * t227 + (t38 * t297 - t204 + t234 - t296) * pkin(1);
t239 = 0.2e1 * t70;
t43 = t128 + t192;
t69 = pkin(3) * t115;
t17 = (-t188 + t281) * t18 + (t239 * t38 + t43 * t56) * t113 + (t116 * t43 + (0.4e1 * t82 - 0.2e1) * t56 * pkin(1)) * t69;
t156 = 1 / pkin(3);
t203 = t155 + t74;
t24 = 0.1e1 / (t239 * t56 + t203 + t50 + t58);
t278 = t156 * t24;
t148 = 0.1e1 / pkin(5);
t27 = 0.1e1 / (t162 + t267);
t280 = t148 * t27;
t194 = t278 * t280;
t20 = pkin(5) * t28 * t298 - t21 * t32;
t13 = (t16 * t293 + t17 * t20 / 0.4e1) * t194;
t14 = (-t16 * t20 / 0.4e1 + t17 * t293) * t194;
t11 = t114 * t14 + t117 * t13;
t291 = pkin(6) * t11;
t10 = t291 * t143;
t145 = pkin(6) ^ 2;
t269 = t10 + t145;
t1 = sqrt(-(-(t143 + pkin(6)) * pkin(6) + t269) * (pkin(6) * (t143 - pkin(6)) + t269));
t12 = t114 * t13 - t117 * t14;
t290 = t1 * t12;
t100 = 0.2e1 / 0.3e1 * t155;
t103 = -t155 / 0.3e1;
t105 = 0.4e1 / 0.3e1 * t162;
t107 = t162 / 0.2e1;
t118 = 15 * t160;
t119 = 15 * t162;
t120 = 10 * t162;
t125 = -0.2e1 * t150;
t126 = -0.5e1 * t150;
t129 = 7 * t160;
t130 = 5 * t160;
t131 = 7 * t162;
t164 = t165 ^ 2;
t136 = 3 * t164;
t137 = 8 * t165;
t138 = 4 * t165;
t149 = t150 ^ 2;
t169 = pkin(3) * t155;
t152 = t169 ^ 2;
t166 = pkin(1) * t162;
t174 = pkin(7) * t165;
t247 = t160 + t164;
t252 = t140 - t150;
t259 = t165 * t150;
t182 = t252 * t162 + t149 / 0.6e1 + t247 - t259;
t181 = t182 + 0.5e1 / 0.6e1 * t170;
t183 = t165 - t199;
t94 = -t150 / 0.2e1;
t53 = t94 + t203;
t184 = t53 * t190;
t191 = -0.6e1 * t199;
t282 = t149 / 0.2e1 - t170 / 0.2e1;
t196 = t136 - 0.3e1 * t259 + t282;
t139 = 3 * t165;
t96 = -0.3e1 / 0.2e1 * t150;
t264 = t96 + t139;
t266 = t152 + t74 * ((t96 + t140) * t162 - 0.3e1 / 0.2e1 * t259 + t247 + t282);
t66 = 0.10e2 / 0.3e1 * t162;
t185 = ((t66 + t252) * t155 + t181) * t191 + (t118 + (-0.9e1 * t150 + (18 * t165)) * t162 + t196) * t155 + (t119 + t264) * t170 + t266;
t253 = t133 + t165;
t206 = t155 + t253;
t287 = pkin(1) * t113;
t186 = t206 * t69 - t287 * (0.3e1 * t155 + t74);
t104 = -0.2e1 / 0.3e1 * t155;
t95 = -0.2e1 / 0.3e1 * t150;
t212 = t95 + t74;
t254 = t120 + t140;
t257 = t104 + t165;
t72 = -t150 - t155;
t60 = t139 + t72;
t270 = t60 * t162;
t187 = -(t127 + ((5 * t162) + t60) * t128 + (t104 + t212) * t74) * t287 + (t170 + (t104 + t95 + t254) * t155 + t130 + 0.2e1 * t270 + t165 * (t95 + t257)) * t69;
t251 = t149 - t170;
t193 = (6 * t164) - 0.6e1 * t259 + t251;
t216 = t100 + t140 + t95;
t265 = t170 + (t100 + t212) * t74;
t93 = -t150 / 0.3e1;
t211 = t93 + t74;
t99 = 0.4e1 / 0.3e1 * t155;
t51 = t99 + t211;
t195 = t51 * t190 + t265 + (t132 + t216) * t155;
t224 = t166 * t69;
t81 = t116 * t82;
t197 = t81 * t224;
t277 = t160 * t82 ^ 2;
t198 = t277 * t69;
t200 = 0.20e2 / 0.3e1 * t162;
t273 = t166 * t81;
t228 = 0.16e2 * t273;
t201 = pkin(7) * t228;
t250 = -t150 + t155;
t205 = t139 + t250;
t235 = pkin(7) * t273;
t207 = 0.8e1 * t235;
t261 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t210 = t113 * t261;
t101 = t155 / 0.3e1;
t91 = -t150 / 0.6e1;
t213 = t101 + t91 + t165;
t214 = t101 + t150 / 0.3e1 + t140;
t215 = t100 + 0.2e1 / 0.3e1 * t150 + t138;
t237 = 0.6e1 * t70;
t217 = pkin(7) * t237;
t238 = 0.4e1 * t70;
t218 = pkin(7) * t238;
t220 = -t287 / 0.2e1;
t222 = 0.4e1 / 0.3e1 * t150 + t99 - (2 * t165);
t225 = t162 * t69;
t229 = 0.12e2 * t276;
t272 = t169 * t112 * t79;
t230 = -0.8e1 * t272;
t231 = 0.4e1 * t276;
t233 = 0.8e1 * t277;
t240 = 0.2e1 * t287;
t242 = pkin(7) * t70;
t243 = 4 * pkin(7);
t244 = t164 + t170;
t245 = t164 - t160;
t255 = t107 + t165;
t256 = t162 / 0.3e1 + t165;
t258 = t165 * t162;
t260 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t268 = 0.4e1 / 0.7e1 * t165 - t150 / 0.7e1;
t271 = t170 * t79 ^ 2;
t274 = t165 * t72;
t67 = t74 ^ 2;
t283 = t67 * (-t155 + t204);
t284 = (-t113 * t166 + t225) * t82;
t35 = -t170 / 0.6e1 + t182;
t64 = t103 + t165;
t39 = t64 * t50;
t45 = t68 + t57;
t73 = -0.3e1 * t155 + t165;
t48 = t73 * t207;
t49 = 0.10e2 * t270;
t54 = -t150 + t203;
t59 = pkin(7) * t239;
t62 = (t138 + t150) * t162;
t65 = -t162 / 0.3e1 + t165;
t71 = -0.30e2 * t150 + (60 * t165);
t76 = -3 * t162 + t165;
t289 = 0.1e1 / ((-0.4e1 * (-t246 * t69 + 0.2e1 * t284 + (t226 * t299 + t113 * (t155 + t294)) * pkin(1)) * t279 + 0.8e1 * pkin(7) * t197 + (((pkin(3) * t295) + 0.8e1 * t162 * t169) * t115 + 0.4e1 * t166 * t210) * t82 - 0.4e1 * t186 * t242 - (t155 * t254 + t130 + t244 + (6 * t258)) * t69 + (t127 + (t120 + 6 * t165) * t155 + t67) * t287) * t18 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t199 + 0.4e1 / 0.9e1 * t155 - t150 / 0.9e1 + t256) * t276 + t65 * t50 + t51 * t260 + (t273 + (t100 + t183 + t91) * t70) * t243) * t279 + t48 + (t261 * t51 + t39) * t229 + t195 * t217 + ((t66 + t216) * t155 + t265) * t191 + t152 + (t119 + t205) * t170 + (t132 * t205 + t140 * t250 + t118 + t136) * t155 + t67 * t54) * t45 + ((t230 * t287 + (t225 * t297 + (t69 - t287) * t59 + t186) * t236) * t18 + (0.8e1 * (t59 + t231 + t76) * t272 + 0.6e1 * (t248 * t231 + (t50 + t51) * t218 + t195) * t68) * t45) * t57) * ((-0.24e2 * (0.4e1 / 0.3e1 * t276 + t59 + t65) * t271 * t287 - 0.12e2 * (-0.8e1 / 0.3e1 * t198 + ((t105 + t213) * t69 - (0.7e1 / 0.6e1 * t155 + t91 + t255) * t287) * t231 + (-t155 * t246 - 0.5e1 / 0.3e1 * t160 + t214 * t162 + t165 * (t93 + t64)) * t69 + (-t170 + (-t200 + t215) * t155 - (3 * t160) + t222 * t162 + t164) * t220 + (-t113 * t160 * t81 + ((t162 + t213) * t69 + (t128 - t246) * t220) * t70) * t243) * t279 + 0.24e2 * t64 * t198 + ((t165 + 0.5e1 / 0.2e1 * t155 + 0.3e1 / 0.2e1 * t162 + t94) * t69 + t73 * t287 / 0.2e1) * t201 - 0.6e1 * ((-0.3e1 * t170 + (-t200 + t222) * t155 + t215 * t162 + t245) * t69 - 0.2e1 * (-0.5e1 / 0.3e1 * t170 + (-t162 + t214) * t155 + t165 * (t103 + t211)) * t287) * t276 - 0.6e1 * t187 * t242 - (t152 + ((21 * t162) + t60) * t170 + (t136 + (35 * t160) + t49 + 0.2e1 * t274) * t155 + (t129 + (t126 + t137 - 0.5e1 * t155) * t162 + t165 * (-t150 + t248)) * t74) * t69 + (0.7e1 * t152 + (t131 + t60) * t127 + ((21 * t160) + (9 * t164) + t49 + 0.6e1 * t274) * t155 + t283) * t287) * t18 + (0.16e2 * (t233 + t201 + (-8 * t160 + 12 * t258) * t82 + (-12 * pkin(7) * t166 + t174 * t301) * t116 - (6 * t258) + t247) * t271 + 0.24e2 * (t257 * t233 + 0.14e2 * (-0.32e2 / 0.21e2 * (t165 + t155 / 0.4e1 + t162 / 0.4e1 - t150 / 0.8e1) * t199 + 0.5e1 / 0.42e2 * t170 + (0.16e2 / 0.21e2 * t162 + t268) * t155 + t160 / 0.7e1 + t268 * t162 + t164 - 0.3e1 / 0.7e1 * t259 + t149 / 0.42e2) * t276 + t65 * t184 - t246 * t170 + (-0.10e2 / 0.3e1 * t160 + (2 * t164) - t259 + t62) * t155 + t35 * t260 + ((-0.2e1 / 0.3e1 * t199 + t92 + t255) * t228 + (-0.8e1 / 0.3e1 * (t256 + t300) * t199 + 0.5e1 / 0.18e2 * t170 + (0.4e1 / 0.3e1 * t165 + t105 + t93) * t155 + t164 + 0.2e1 / 0.3e1 * t258 - 0.2e1 / 0.3e1 * t259 - t160 / 0.3e1 + t149 / 0.18e2) * t237) * pkin(7)) * t279 + 0.16e2 * (-0.6e1 * t165 * t155 + t244) * t277 + 0.32e2 * (t261 * t50 + t53 * t73) * t235 + 0.24e2 * (t64 * t184 - t152 + (-t66 + t249) * t170 + (t62 + t170 / 0.6e1 - t149 / 0.6e1 + t245) * t155 + t35 * t165) * t276 + 0.8e1 * t185 * t242 - 0.8e1 * ((t131 + t264) * t170 + (t129 + (t126 + (10 * t165)) * t162 + t196) * t155 + t266) * t199 + t170 ^ 2 + (t125 + t138 + (28 * t162)) * t152 + (t162 * t71 + (70 * t160) + t193) * t170 + (t132 * t193 + t140 * t251 - 0.6e1 * t164 * t150 + t71 * t160 + (28 * t166 ^ 2) + (4 * t174 ^ 2)) * t155 + t54 * t283) * t45 + (((0.4e1 * t284 + (t69 + t240) * t59 + t76 * t69 + (t94 + t206) * t240) * t230 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t155 + t107 + t91) * t69 + pkin(1) * t210) * t232 + (-0.8e1 * t197 + ((t100 + t93 + t253) * t69 - (0.8e1 / 0.3e1 * t155 + t211) * t287) * t238) * pkin(7) + t187) * t68) * t18 + (0.32e2 * (t207 + (-0.4e1 * t113 * t224 + t295 + (t296 + t125 + t137) * t162) * t82 + (-t162 + t183 + t300) * t218 + t50 * t260 + t76 * t53) * t272 + 0.8e1 * (t48 + (t261 * t53 + t39) * t229 + (t184 + (t132 + t252) * t155 + t181) * t217 + t185) * t68) * t45) * t57);
t288 = 0.1e1 / pkin(6) / 0.2e1;
t159 = 1 / pkin(2);
t158 = pkin(2) ^ 2;
t5 = 0.1e1 / (t158 + t269);
t286 = t159 * t5;
t163 = 1 / pkin(1);
t275 = t163 * t27;
t223 = t289 / 0.4e1;
t151 = 1 / pkin(4);
t221 = t151 / (pkin(3) ^ 2) * t24;
t219 = t5 * t288;
t209 = t280 / 0.2e1;
t208 = t278 / 0.2e1;
t202 = (t151 * t156) / 0.2e1;
t41 = -t112 * t116 + t263;
t40 = -t112 * t113 - t262;
t33 = -pkin(1) * t37 + pkin(5);
t29 = t162 - t241;
t9 = pkin(2) * t11 + pkin(6);
t8 = -pkin(2) - t291;
t7 = t10 + (2 * t158);
t6 = t10 + 0.2e1 * t145;
t3 = (t16 * t18 / 0.4e1 + t17 * t223) * t221;
t2 = (t16 * t223 - t17 * t18 / 0.4e1) * t221;
t4 = [qJ(1); atan2(t17 * t208, t16 * t208); atan2((pkin(6) * t12 * t6 - t1 * t8) * t219, (-pkin(6) * t290 - t6 * t8) * t219); atan2(t18 * t202, t202 * t289); atan2((pkin(1) * t29 * t298 + t21 * t33) * t275 / 0.2e1, -(-pkin(1) * t285 + 0.2e1 * t29 * t33) * t275 / 0.2e1); atan2(t12, t11); qJ(2); atan2(t20 * t209, t19 * t209); atan2(t159 * t1 * t288, -t11); atan2(t2 * t40 + t3 * t41, -t2 * t41 + t3 * t40); atan2(t148 * t163 * t21 / 0.2e1, t37); atan2((pkin(2) * t12 * t7 + t1 * t9) * t286 / 0.2e1, -(-pkin(2) * t290 + t7 * t9) * t286 / 0.2e1); 0; 0; 0;];
jv = t4(:);