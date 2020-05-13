% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh1m1TE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1TE_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energykin_fixb_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 04:17:03
% EndTime: 2020-04-13 04:17:41
% DurationCPUTime: 38.48s
% Computational Cost: add. (989643->233), mult. (1534373->561), div. (63348->34), fcn. (957108->26), ass. (0->254)
t302 = pkin(13) ^ 2;
t159 = sin(pkin(20));
t163 = cos(pkin(20));
t165 = sin(qJ(3));
t169 = cos(qJ(3));
t207 = t169 * t159 + t165 * t163;
t270 = pkin(6) * t207;
t227 = pkin(1) * t270;
t133 = 0.2e1 * t227;
t185 = pkin(6) ^ 2;
t236 = t133 + t185;
t294 = -pkin(2) - pkin(13);
t117 = (pkin(1) - t294) * (pkin(1) + t294) + t236;
t293 = -pkin(2) + pkin(13);
t118 = (pkin(1) - t293) * (pkin(1) + t293) + t236;
t252 = t118 * t117;
t193 = sqrt(-t252);
t273 = t193 / 0.2e1;
t190 = pkin(2) ^ 2;
t192 = pkin(1) ^ 2;
t233 = t185 + t192;
t218 = t233 - t302;
t123 = t133 + t190 + t218;
t132 = -pkin(1) - t270;
t144 = t165 * t159 - t169 * t163;
t269 = pkin(6) * t144;
t104 = t123 * t269 - t132 * t193;
t301 = pkin(6) * t104;
t299 = -0.2e1 * pkin(1);
t186 = pkin(5) ^ 2;
t154 = pkin(23) + pkin(22);
t151 = sin(t154);
t152 = cos(t154);
t184 = pkin(7) ^ 2;
t166 = sin(qJ(2));
t167 = sin(pkin(19));
t170 = cos(qJ(2));
t171 = cos(pkin(19));
t145 = t166 * t171 - t170 * t167;
t268 = pkin(7) * t145;
t235 = t268 * t299 + t192;
t131 = t184 + t235;
t232 = pkin(3) ^ 2 - pkin(8) ^ 2;
t125 = t131 + t232;
t136 = pkin(1) - t268;
t146 = t166 * t167 + t170 * t171;
t292 = -pkin(8) - pkin(3);
t119 = (pkin(7) - t292) * (pkin(7) + t292) + t235;
t291 = -pkin(8) + pkin(3);
t120 = (pkin(7) - t291) * (pkin(7) + t291) + t235;
t194 = sqrt(-t120 * t119);
t107 = pkin(7) * t146 * t125 + t136 * t194;
t240 = t169 * t107;
t247 = t146 * t194;
t106 = -pkin(7) * t247 + t136 * t125;
t243 = t165 * t106;
t203 = t243 / 0.2e1 + t240 / 0.2e1;
t129 = 0.1e1 / t131;
t189 = 0.1e1 / pkin(3);
t249 = t129 * t189;
t94 = t203 * t249;
t241 = t169 * t106;
t242 = t165 * t107;
t202 = -t241 / 0.2e1 + t242 / 0.2e1;
t95 = t202 * t249;
t69 = t151 * t95 - t152 * t94;
t285 = pkin(5) * t69;
t258 = -0.2e1 * pkin(4) * t285 + t186;
t290 = -pkin(9) - pkin(11);
t59 = (pkin(4) - t290) * (pkin(4) + t290) + t258;
t289 = -pkin(9) + pkin(11);
t60 = (pkin(4) - t289) * (pkin(4) + t289) + t258;
t195 = sqrt(-t60 * t59);
t70 = t151 * t94 + t152 * t95;
t259 = t195 * t70;
t234 = pkin(9) ^ 2 - pkin(11) ^ 2;
t187 = pkin(4) ^ 2;
t65 = t187 + t258;
t61 = t65 + t234;
t66 = -pkin(4) + t285;
t42 = -pkin(5) * t259 - t66 * t61;
t298 = t42 / 0.2e1;
t63 = 0.1e1 / t65;
t297 = t63 / 0.2e1;
t142 = t146 * qJD(2);
t141 = t145 * qJD(2);
t228 = pkin(1) * pkin(7) * t142;
t253 = 0.2e1 * (t119 + t120) * t228 / t194;
t214 = -t253 / 0.2e1;
t201 = t141 * t194 + t146 * t214;
t80 = ((t136 * t299 - t125) * t142 + t201) * pkin(7);
t296 = -t80 / 0.2e1;
t226 = -0.2e1 * t142 * t146;
t248 = t142 * t194;
t81 = t136 * t253 / 0.2e1 + t184 * pkin(1) * t226 + (-t141 * t125 - t248) * pkin(7);
t295 = t81 / 0.2e1;
t160 = cos(pkin(23));
t212 = 0.1e1 / t131 ^ 2 * t228;
t244 = t160 * t107;
t245 = t160 * t106;
t156 = sin(pkin(23));
t246 = t156 * t107;
t255 = t106 * t156;
t280 = t156 / 0.2e1;
t90 = (-t245 / 0.2e1 + t246 / 0.2e1) * t249;
t197 = t90 ^ 2;
t88 = 0.1e1 / t197;
t91 = (t244 / 0.2e1 + t255 / 0.2e1) * t249;
t89 = t91 ^ 2;
t35 = qJD(2) + (((t160 * t295 + t80 * t280) * t129 + (t244 + t255) * t212) / t90 - ((t160 * t296 + t81 * t280) * t129 + (-t245 + t246) * t212) * t91 * t88) * t189 / (t89 * t88 + 0.1e1);
t288 = pkin(4) * t35;
t277 = -t169 / 0.2e1;
t50 = ((-t241 + t242) * t212 + (t203 * qJD(3) + t165 * t295 + t80 * t277) * t129) * t189;
t51 = ((-t240 - t243) * t212 + (t202 * qJD(3) + t165 * t296 + t81 * t277) * t129) * t189;
t45 = t151 * t50 + t152 * t51;
t287 = pkin(4) * t45;
t284 = pkin(5) * t70;
t43 = -t195 * t66 + t61 * t284;
t286 = pkin(5) * t43;
t122 = t190 - t218 - 0.2e1 * t227;
t283 = -t122 / 0.2e1;
t282 = t122 / 0.2e1;
t128 = t133 + t233;
t126 = 0.1e1 / t128;
t281 = t126 / 0.2e1;
t158 = sin(pkin(21));
t279 = t158 / 0.2e1;
t161 = cos(pkin(22));
t278 = -t161 / 0.2e1;
t172 = cos(pkin(18));
t276 = t172 / 0.2e1;
t174 = qJD(1) ^ 2;
t275 = t174 / 0.2e1;
t274 = -t193 / 0.2e1;
t272 = cos(qJ(4));
t135 = t144 * qJD(3);
t271 = pkin(1) * t135;
t230 = pkin(5) * t287;
t267 = 0.2e1 * (t59 + t60) * t230 / t195;
t266 = pkin(1) * qJD(2);
t162 = cos(pkin(21));
t265 = t162 * t63;
t124 = t131 - t232;
t137 = pkin(1) * t145 - pkin(7);
t105 = -pkin(1) * t247 - t137 * t124;
t108 = pkin(1) * t146 * t124 - t137 * t194;
t168 = sin(pkin(18));
t183 = 0.1e1 / pkin(8);
t250 = t129 * t183;
t97 = (t108 * t276 + t105 * t168 / 0.2e1) * t250;
t264 = t174 * t97;
t176 = 0.1e1 / pkin(13);
t103 = -t132 * t123 - t193 * t269;
t102 = 0.1e1 / t103 ^ 2;
t127 = 0.1e1 / t128 ^ 2;
t134 = t207 * qJD(3);
t191 = 0.1e1 / pkin(2);
t229 = pkin(6) * t271;
t254 = (t117 + t118) * t229 / t273;
t215 = -t254 / 0.2e1;
t49 = qJD(2) + 0.2e1 * (-t102 * ((-t135 * t123 - t134 * t193 + t144 * t215) * t281 + (t103 * t127 + t126 * t132) * t271) * t301 + 0.1e1 / t103 * ((t132 * t215 + (t134 * t123 - t135 * t193) * pkin(6)) * t281 + (-t126 * t144 * t185 + t127 * t301) * t271)) / (t104 ^ 2 * t102 + 0.1e1) * t128 * pkin(2) * t191;
t263 = t176 * t49;
t179 = 0.1e1 / pkin(11);
t262 = t179 * t63;
t181 = 0.1e1 / pkin(9);
t261 = t181 * t63;
t260 = t195 * t45;
t257 = qJD(1) * pkin(16);
t208 = t108 * t212;
t209 = t105 * t212;
t213 = t129 * t276;
t251 = t129 * t168;
t79 = ((0.2e1 * t137 * pkin(7) - t124) * t142 + t201) * pkin(1);
t82 = t137 * t214 + t192 * pkin(7) * t226 + (-t141 * t124 - t248) * pkin(1);
t96 = (t105 * t276 - t168 * t108 / 0.2e1) * t250;
t198 = t96 ^ 2;
t92 = 0.1e1 / t198;
t93 = t97 ^ 2;
t36 = ((t82 * t213 + t172 * t208 + t79 * t251 / 0.2e1 + t168 * t209) / t96 - (t79 * t213 + t172 * t209 - t82 * t251 / 0.2e1 - t168 * t208) * t97 * t92) / (t93 * t92 + 0.1e1) * t183;
t256 = qJD(1) * t36;
t239 = t170 * t174;
t173 = qJD(2) ^ 2;
t238 = t173 * t192;
t237 = t176 * t191;
t231 = qJD(1) * qJD(2);
t155 = qJD(2) + qJD(3);
t196 = t122 ^ 2;
t121 = 0.1e1 / t196;
t46 = (0.1e1 / t122 * t254 / 0.2e1 - 0.2e1 * t193 * t121 * t229) / (-t121 * t252 + 0.1e1) + t49;
t225 = t46 * t263;
t224 = t90 * t266;
t223 = t91 * t266;
t222 = t169 * t266;
t221 = t165 * t266;
t220 = -t267 / 0.2e1;
t219 = t63 * t279;
t217 = qJD(1) * t126 * t191;
t64 = 0.1e1 / t65 ^ 2;
t216 = t64 * t230;
t149 = t166 * qJD(1) * pkin(1) - t257;
t139 = (t165 * t170 + t166 * t169) * qJD(1);
t140 = (-t165 * t166 + t169 * t170) * qJD(1);
t62 = t65 - t234;
t67 = -pkin(4) * t69 + pkin(5);
t205 = -pkin(4) * t259 + t67 * t62;
t206 = pkin(4) * t70 * t62 + t195 * t67;
t29 = (t205 * t279 + t206 * t162 / 0.2e1) * t262;
t30 = (-t205 * t162 / 0.2e1 + t206 * t279) * t262;
t20 = -t29 * t139 + t30 * t140;
t148 = t155 * pkin(5) + t221;
t23 = t29 * t148 - t30 * t222;
t116 = -t140 * pkin(5) + t149;
t44 = -t151 * t51 + t152 * t50;
t204 = -t44 * t195 + t70 * t220;
t22 = t30 * t148 + t29 * t222;
t200 = t206 * t216;
t199 = t205 * t216;
t164 = sin(qJ(4));
t157 = sin(pkin(22));
t153 = pkin(16) ^ 2 * t275;
t147 = t149 ^ 2 / 0.2e1;
t85 = (t103 * t166 + t104 * t170) * t217 / 0.2e1;
t84 = (t103 * t170 / 0.2e1 - t104 * t166 / 0.2e1) * t217;
t83 = t85 * pkin(2) - t257;
t72 = (-t166 * t90 - t170 * t91) * qJD(1);
t71 = (-t166 * t91 + t170 * t90) * qJD(1);
t58 = (t84 * t273 + t85 * t282) * t237;
t56 = (t274 * t85 + t282 * t84) * t237;
t54 = (-t157 * t71 - t161 * t72) * pkin(4) + t149;
t48 = t49 ^ 2;
t41 = 0.1e1 / t42 ^ 2;
t34 = t157 * t288 + t223;
t33 = t161 * t288 + t224;
t28 = (t42 * t278 - t157 * t43 / 0.2e1) * t261;
t27 = (t157 * t298 + t43 * t278) * t261;
t26 = 0.1e1 / t30 ^ 2;
t21 = t30 * t139 + t29 * t140;
t18 = -qJD(4) + t20;
t17 = -t27 * t71 + t28 * t72;
t16 = t27 * t72 + t28 * t71;
t15 = t67 * t267 / 0.2e1 - 0.2e1 * t187 * t45 * t284 + (t44 * t62 - t260) * pkin(4);
t14 = ((-0.2e1 * t67 * pkin(5) - t62) * t45 + t204) * pkin(4);
t13 = -t27 * t34 + t28 * t33;
t12 = t27 * t33 + t28 * t34;
t11 = -t20 * pkin(10) - t21 * pkin(12) + t116;
t10 = t35 + (((t66 * t220 + (t44 * t61 - t260) * pkin(5)) * t297 + (-t186 * t63 * t70 + t64 * t286) * t287) / t298 - 0.2e1 * ((-t45 * t61 + t204) * t297 + (t42 * t64 + t63 * t66) * t287) * t41 * t286) * pkin(9) * t181 / (t41 * t43 ^ 2 + 0.1e1) * t65;
t9 = ((t14 * t219 + t158 * t199 + t15 * t265 / 0.2e1 + t162 * t200) / t30 - (-t14 * t265 / 0.2e1 - t162 * t199 + t15 * t219 + t158 * t200) * t29 * t26) / (t26 * t29 ^ 2 + 0.1e1) * t179 + t155;
t7 = t9 * pkin(12) + t23;
t6 = -t9 * pkin(10) - t22;
t5 = t164 * t9 + t272 * t21;
t3 = t164 * t21 - t272 * t9;
t2 = t164 * t11 + t272 * t7;
t1 = t272 * t11 - t164 * t7;
t4 = [0, 0, 0, 0, 0, t275, 0, 0, 0, 0, t170 ^ 2 * t275, -t166 * t239, t170 * t231, t166 ^ 2 * t275, -t166 * t231, t173 / 0.2e1, -t174 * pkin(16) * t166, -pkin(16) * t239, 0, t153, t139 ^ 2 / 0.2e1, t139 * t140, t139 * t155, t140 ^ 2 / 0.2e1, t140 * t155, t155 ^ 2 / 0.2e1, -t149 * t140 + t155 * t221, t149 * t139 + t155 * t222, (-t139 * t165 - t140 * t169) * t266, t147 + (t169 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1) * t238, t21 ^ 2 / 0.2e1, t20 * t21, t9 * t21, t20 ^ 2 / 0.2e1, t9 * t20, t9 ^ 2 / 0.2e1, -t116 * t20 + t22 * t9, t116 * t21 - t23 * t9, t20 * t23 - t21 * t22, t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t116 ^ 2 / 0.2e1, t5 ^ 2 / 0.2e1, -t3 * t5, -t5 * t18, t3 ^ 2 / 0.2e1, t3 * t18, t18 ^ 2 / 0.2e1, -t1 * t18 + t3 * t6, t18 * t2 + t5 * t6, -t1 * t5 - t2 * t3, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t93 * t275, t96 * t264, t97 * t256, t198 * t275, t96 * t256, t36 ^ 2 / 0.2e1, -t174 * pkin(15) * t96, pkin(15) * t264, 0, pkin(15) ^ 2 * t275, t71 ^ 2 / 0.2e1, t71 * t72, t71 * t35, t72 ^ 2 / 0.2e1, t35 * t72, t35 ^ 2 / 0.2e1, -t149 * t72 + t224 * t35, t149 * t71 - t223 * t35, (-t71 * t90 + t72 * t91) * t266, t147 + (t89 / 0.2e1 + t197 / 0.2e1) * t238, t84 ^ 2 / 0.2e1, -t85 * t84, t49 * t84, t85 ^ 2 / 0.2e1, -t49 * t85, t48 / 0.2e1, -t85 * t257, -t84 * t257, 0, t153, t56 ^ 2 / 0.2e1, -t58 * t56, -t46 * t56, t58 ^ 2 / 0.2e1, t46 * t58, t46 ^ 2 / 0.2e1, t225 * t283 - t83 * t58, t225 * t273 - t83 * t56, (t58 * t274 + t56 * t283) * t263, t83 ^ 2 / 0.2e1 + (-t252 / 0.8e1 + t196 / 0.8e1) * t48 / t302, t16 ^ 2 / 0.2e1, t17 * t16, t10 * t16, t17 ^ 2 / 0.2e1, t10 * t17, t10 ^ 2 / 0.2e1, t10 * t13 - t17 * t54, -t10 * t12 + t16 * t54, t12 * t17 - t13 * t16, t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1;];
T_reg = t4;
