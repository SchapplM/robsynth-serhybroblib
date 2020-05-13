% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh1m1DE2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energykin_fixb_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 06:37:03
% EndTime: 2020-04-15 06:37:46
% DurationCPUTime: 42.85s
% Computational Cost: add. (1187200->232), mult. (1831044->561), div. (77522->33), fcn. (1145176->44), ass. (0->260)
t171 = sin(pkin(20));
t175 = cos(pkin(20));
t177 = sin(qJ(3));
t181 = cos(qJ(3));
t211 = t171 * t181 + t175 * t177;
t274 = pkin(6) * t211;
t233 = pkin(1) * t274;
t145 = 0.2e1 * t233;
t201 = pkin(2) ^ 2;
t196 = pkin(6) ^ 2;
t203 = pkin(1) ^ 2;
t240 = t196 + t203;
t224 = -pkin(13) ^ 2 + t240;
t135 = t145 + t201 + t224;
t144 = -pkin(1) - t274;
t243 = t145 + t196;
t296 = -pkin(2) - pkin(13);
t129 = (pkin(1) - t296) * (pkin(1) + t296) + t243;
t295 = -pkin(2) + pkin(13);
t130 = (pkin(1) - t295) * (pkin(1) + t295) + t243;
t258 = t130 * t129;
t204 = sqrt(-t258);
t156 = t171 * t177 - t175 * t181;
t273 = pkin(6) * t156;
t113 = t135 * t273 - t144 * t204;
t302 = pkin(6) * t113;
t300 = -0.2e1 * pkin(1);
t198 = pkin(4) ^ 2;
t197 = pkin(5) ^ 2;
t195 = pkin(7) ^ 2;
t178 = sin(qJ(2));
t179 = sin(pkin(19));
t182 = cos(qJ(2));
t183 = cos(pkin(19));
t157 = t178 * t183 - t179 * t182;
t272 = pkin(7) * t157;
t242 = t272 * t300 + t203;
t143 = t195 + t242;
t239 = pkin(3) ^ 2 - pkin(8) ^ 2;
t137 = t143 + t239;
t148 = pkin(1) - t272;
t158 = t178 * t179 + t182 * t183;
t294 = -pkin(8) - pkin(3);
t131 = (pkin(7) - t294) * (pkin(7) + t294) + t242;
t293 = -pkin(8) + pkin(3);
t132 = (pkin(7) - t293) * (pkin(7) + t293) + t242;
t205 = sqrt(-t132 * t131);
t116 = pkin(7) * t137 * t158 + t148 * t205;
t246 = t181 * t116;
t253 = t158 * t205;
t115 = -pkin(7) * t253 + t137 * t148;
t249 = t177 * t115;
t209 = t249 / 0.2e1 + t246 / 0.2e1;
t141 = 0.1e1 / t143;
t200 = 0.1e1 / pkin(3);
t255 = t141 * t200;
t105 = t209 * t255;
t247 = t181 * t115;
t248 = t177 * t116;
t208 = -t247 / 0.2e1 + t248 / 0.2e1;
t106 = t208 * t255;
t166 = pkin(23) + pkin(22);
t163 = sin(t166);
t164 = cos(t166);
t81 = -t105 * t164 + t106 * t163;
t285 = pkin(5) * t81;
t264 = -0.2e1 * pkin(4) * t285 + t197;
t77 = t198 + t264;
t75 = 0.1e1 / t77;
t299 = t75 / 0.2e1;
t154 = t158 * qJD(2);
t153 = t157 * qJD(2);
t235 = pkin(1) * pkin(7) * t154;
t259 = 0.2e1 * (t131 + t132) * t235 / t205;
t221 = -t259 / 0.2e1;
t207 = t153 * t205 + t158 * t221;
t94 = ((t148 * t300 - t137) * t154 + t207) * pkin(7);
t298 = -t94 / 0.2e1;
t232 = -0.2e1 * t154 * t158;
t254 = t154 * t205;
t95 = t148 * t259 / 0.2e1 + t195 * pkin(1) * t232 + (-t137 * t153 - t254) * pkin(7);
t297 = t95 / 0.2e1;
t292 = -pkin(9) - pkin(11);
t291 = -pkin(9) + pkin(11);
t192 = 0.1e1 / pkin(9);
t225 = t192 * t299;
t68 = (pkin(4) - t292) * (pkin(4) + t292) + t264;
t69 = (pkin(4) - t291) * (pkin(4) + t291) + t264;
t206 = sqrt(-t69 * t68);
t82 = t105 * t163 + t106 * t164;
t265 = t206 * t82;
t241 = pkin(9) ^ 2 - pkin(11) ^ 2;
t73 = t77 + t241;
t78 = -pkin(4) + t285;
t45 = -pkin(5) * t265 - t73 * t78;
t284 = pkin(5) * t82;
t47 = -t206 * t78 + t73 * t284;
t37 = atan2(t47 * t225, t45 * t225);
t290 = sin(t37);
t112 = -t135 * t144 - t204 * t273;
t111 = 0.1e1 / t112 ^ 2;
t140 = t145 + t240;
t138 = 0.1e1 / t140;
t139 = 0.1e1 / t140 ^ 2;
t146 = t211 * qJD(3);
t147 = t156 * qJD(3);
t202 = 0.1e1 / pkin(2);
t275 = pkin(1) * t147;
t236 = pkin(6) * t275;
t260 = 0.2e1 * (t129 + t130) * t236 / t204;
t222 = -t260 / 0.2e1;
t283 = t138 / 0.2e1;
t55 = qJD(2) + 0.2e1 * (-t111 * ((-t147 * t135 - t146 * t204 + t156 * t222) * t283 + (t112 * t139 + t138 * t144) * t275) * t302 + 0.1e1 / t112 * ((t144 * t222 + (t135 * t146 - t147 * t204) * pkin(6)) * t283 + (-t138 * t156 * t196 + t139 * t302) * t275)) * t140 / (t111 * t113 ^ 2 + 0.1e1) * pkin(2) * t202;
t289 = pkin(2) * t55;
t172 = cos(pkin(23));
t251 = t172 * t115;
t168 = sin(pkin(23));
t252 = t168 * t116;
t102 = (-t251 / 0.2e1 + t252 / 0.2e1) * t255;
t101 = 0.1e1 / t102 ^ 2;
t250 = t172 * t116;
t261 = t115 * t168;
t103 = (t250 / 0.2e1 + t261 / 0.2e1) * t255;
t217 = 0.1e1 / t143 ^ 2 * t235;
t282 = t168 / 0.2e1;
t42 = qJD(2) + (((t172 * t297 + t94 * t282) * t141 + (t250 + t261) * t217) / t102 - ((t172 * t298 + t95 * t282) * t141 + (-t251 + t252) * t217) * t103 * t101) * t200 / (t101 * t103 ^ 2 + 0.1e1);
t288 = pkin(4) * t42;
t280 = -t181 / 0.2e1;
t59 = ((-t247 + t248) * t217 + (t209 * qJD(3) + t177 * t297 + t94 * t280) * t141) * t200;
t60 = ((-t246 - t249) * t217 + (t208 * qJD(3) + t177 * t298 + t95 * t280) * t141) * t200;
t50 = t163 * t59 + t164 * t60;
t287 = pkin(4) * t50;
t286 = pkin(5) * t47;
t170 = sin(pkin(21));
t281 = t170 / 0.2e1;
t184 = cos(pkin(18));
t279 = t184 / 0.2e1;
t186 = qJD(1) ^ 2;
t278 = t186 / 0.2e1;
t277 = t202 / 0.2e1;
t276 = cos(qJ(4));
t237 = pkin(5) * t287;
t271 = 0.2e1 * (t68 + t69) * t237 / t206;
t270 = pkin(1) * qJD(2);
t174 = cos(pkin(21));
t269 = t174 * t75;
t136 = t143 - t239;
t149 = pkin(1) * t157 - pkin(7);
t114 = -pkin(1) * t253 - t136 * t149;
t117 = pkin(1) * t136 * t158 - t149 * t205;
t180 = sin(pkin(18));
t194 = 0.1e1 / pkin(8);
t256 = t141 * t194;
t107 = (t114 * t279 - t180 * t117 / 0.2e1) * t256;
t108 = (t117 * t279 + t114 * t180 / 0.2e1) * t256;
t90 = atan2(t108, t107);
t87 = cos(t90);
t268 = t186 * t87;
t190 = 0.1e1 / pkin(11);
t267 = t190 * t75;
t266 = t206 * t50;
t263 = qJD(1) * pkin(16);
t104 = 0.1e1 / t107 ^ 2;
t212 = t117 * t217;
t213 = t114 * t217;
t219 = t141 * t279;
t257 = t141 * t180;
t93 = ((0.2e1 * pkin(7) * t149 - t136) * t154 + t207) * pkin(1);
t96 = t149 * t221 + t203 * pkin(7) * t232 + (-t136 * t153 - t254) * pkin(1);
t43 = ((t96 * t219 + t184 * t212 + t93 * t257 / 0.2e1 + t180 * t213) / t107 - (t93 * t219 + t184 * t213 - t96 * t257 / 0.2e1 - t180 * t212) * t108 * t104) / (t104 * t108 ^ 2 + 0.1e1) * t194;
t262 = qJD(1) * t43;
t245 = t182 * t186;
t185 = qJD(2) ^ 2;
t244 = t185 * t203;
t238 = qJD(1) * qJD(2);
t167 = qJD(2) + qJD(3);
t134 = t201 - t224 - 0.2e1 * t233;
t133 = 0.1e1 / t134 ^ 2;
t51 = (0.1e1 / t134 * t260 / 0.2e1 - 0.2e1 * t204 * t133 * t236) / (-t133 * t258 + 0.1e1) + t55;
t234 = t51 * t289;
t88 = atan2(t103, t102);
t83 = sin(t88);
t231 = t83 * t270;
t84 = cos(t88);
t230 = t84 * t270;
t229 = t177 * t270;
t228 = t181 * t270;
t227 = -t271 / 0.2e1;
t226 = t75 * t281;
t76 = 0.1e1 / t77 ^ 2;
t223 = t76 * t237;
t220 = t138 * t277;
t218 = 0.1e1 / pkin(13) * t277;
t151 = (t177 * t182 + t178 * t181) * qJD(1);
t152 = (-t177 * t178 + t181 * t182) * qJD(1);
t74 = t77 - t241;
t79 = -pkin(4) * t81 + pkin(5);
t46 = -pkin(4) * t265 + t74 * t79;
t48 = pkin(4) * t74 * t82 + t206 * t79;
t34 = (t46 * t281 + t48 * t174 / 0.2e1) * t267;
t35 = (-t46 * t174 / 0.2e1 + t48 * t281) * t267;
t31 = atan2(t34, t35);
t28 = sin(t31);
t29 = cos(t31);
t17 = t151 * t28 - t29 * t152;
t161 = t178 * qJD(1) * pkin(1) - t263;
t215 = t170 * t223;
t214 = t174 * t223;
t160 = pkin(5) * t167 + t229;
t21 = t28 * t160 - t29 * t228;
t128 = -pkin(5) * t152 + t161;
t49 = -t163 * t60 + t164 * t59;
t210 = -t49 * t206 + t82 * t227;
t20 = t160 * t29 + t28 * t228;
t176 = sin(qJ(4));
t173 = cos(pkin(22));
t169 = sin(pkin(22));
t165 = pkin(16) ^ 2 * t278;
t159 = t161 ^ 2 / 0.2e1;
t122 = atan2(t204 * t218, t134 * t218);
t120 = cos(t122);
t119 = sin(t122);
t99 = atan2(t113 * t220, t112 * t220);
t98 = cos(t99);
t97 = sin(t99);
t86 = sin(t90);
t72 = (-t178 * t97 + t182 * t98) * qJD(1);
t70 = (-t178 * t98 - t182 * t97) * qJD(1);
t67 = -pkin(2) * t70 - t263;
t65 = (-t178 * t83 + t182 * t84) * qJD(1);
t63 = (-t178 * t84 - t182 * t83) * qJD(1);
t57 = t119 * t70 + t120 * t72;
t56 = t119 * t72 - t120 * t70;
t54 = t55 ^ 2;
t52 = (-t169 * t65 - t173 * t63) * pkin(4) + t161;
t44 = 0.1e1 / t45 ^ 2;
t41 = t173 * t288 + t230;
t40 = t169 * t288 + t231;
t36 = cos(t37);
t33 = 0.1e1 / t35 ^ 2;
t27 = -t169 * t290 - t173 * t36;
t26 = t169 * t36 - t173 * t290;
t23 = t79 * t271 / 0.2e1 - 0.2e1 * t198 * t50 * t284 + (t49 * t74 - t266) * pkin(4);
t22 = ((-0.2e1 * pkin(5) * t79 - t74) * t50 + t210) * pkin(4);
t19 = t151 * t29 + t152 * t28;
t16 = qJD(4) + t17;
t15 = t26 * t63 + t27 * t65;
t14 = -t26 * t65 + t27 * t63;
t13 = t26 * t41 + t27 * t40;
t12 = -t26 * t40 + t27 * t41;
t11 = pkin(10) * t17 - pkin(12) * t19 + t128;
t10 = t42 + 0.2e1 * (((t78 * t227 + (t49 * t73 - t266) * pkin(5)) * t299 + (-t197 * t75 * t82 + t76 * t286) * t287) / t45 - ((-t50 * t73 + t210) * t299 + (t45 * t76 + t75 * t78) * t287) * t44 * t286) * pkin(9) * t192 / (t44 * t47 ^ 2 + 0.1e1) * t77;
t9 = ((t22 * t226 + t46 * t215 + t23 * t269 / 0.2e1 + t48 * t214) / t35 - (-t22 * t269 / 0.2e1 - t46 * t214 + t23 * t226 + t48 * t215) * t34 * t33) / (t33 * t34 ^ 2 + 0.1e1) * t190 + t167;
t7 = pkin(12) * t9 + t21;
t6 = -pkin(10) * t9 - t20;
t5 = t176 * t9 + t276 * t19;
t3 = t176 * t19 - t276 * t9;
t2 = t176 * t11 + t276 * t7;
t1 = t276 * t11 - t176 * t7;
t4 = [0, 0, 0, 0, 0, t278, 0, 0, 0, 0, t182 ^ 2 * t278, -t178 * t245, t182 * t238, t178 ^ 2 * t278, -t178 * t238, t185 / 0.2e1, -t186 * pkin(16) * t178, -pkin(16) * t245, 0, t165, t151 ^ 2 / 0.2e1, t151 * t152, t151 * t167, t152 ^ 2 / 0.2e1, t152 * t167, t167 ^ 2 / 0.2e1, -t152 * t161 + t167 * t229, t151 * t161 + t167 * t228, (-t151 * t177 - t152 * t181) * t270, t159 + (t181 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1) * t244, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t9, t17 ^ 2 / 0.2e1, -t17 * t9, t9 ^ 2 / 0.2e1, t128 * t17 + t20 * t9, t128 * t19 - t21 * t9, -t17 * t21 - t19 * t20, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1, t5 ^ 2 / 0.2e1, -t5 * t3, t5 * t16, t3 ^ 2 / 0.2e1, -t16 * t3, t16 ^ 2 / 0.2e1, t1 * t16 + t3 * t6, -t16 * t2 + t5 * t6, -t1 * t5 - t2 * t3, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t86 ^ 2 * t278, t86 * t268, t86 * t262, t87 ^ 2 * t278, t87 * t262, t43 ^ 2 / 0.2e1, -pkin(15) * t268, t186 * pkin(15) * t86, 0, pkin(15) ^ 2 * t278, t65 ^ 2 / 0.2e1, t63 * t65, t42 * t65, t63 ^ 2 / 0.2e1, t42 * t63, t42 ^ 2 / 0.2e1, -t161 * t63 + t230 * t42, t161 * t65 - t231 * t42, (t63 * t83 - t65 * t84) * t270, t159 + (t83 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1) * t244, t72 ^ 2 / 0.2e1, t72 * t70, t72 * t55, t70 ^ 2 / 0.2e1, t55 * t70, t54 / 0.2e1, t70 * t263, -t72 * t263, 0, t165, t57 ^ 2 / 0.2e1, -t56 * t57, -t51 * t57, t56 ^ 2 / 0.2e1, t51 * t56, t51 ^ 2 / 0.2e1, -t120 * t234 - t56 * t67, t119 * t234 - t57 * t67, (-t119 * t56 - t120 * t57) * t289, t67 ^ 2 / 0.2e1 + (t119 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1) * t54 * t201, t15 ^ 2 / 0.2e1, t15 * t14, t10 * t15, t14 ^ 2 / 0.2e1, t10 * t14, t10 ^ 2 / 0.2e1, t10 * t12 - t14 * t52, -t10 * t13 + t15 * t52, -t12 * t15 + t13 * t14, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1;];
T_reg = t4;
