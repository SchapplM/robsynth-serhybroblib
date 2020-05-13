% Calculate kinetic energy for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1TE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:10:21
% EndTime: 2020-04-12 20:11:04
% DurationCPUTime: 41.90s
% Computational Cost: add. (999549->389), mult. (1549325->684), div. (63732->33), fcn. (966768->28), ass. (0->269)
t192 = sin(pkin(20));
t196 = cos(pkin(20));
t198 = sin(qJ(3));
t204 = cos(qJ(3));
t236 = t192 * t204 + t196 * t198;
t287 = pkin(6) * t236;
t251 = pkin(1) * t287;
t163 = 0.2e1 * t251;
t219 = pkin(6) ^ 2;
t259 = t163 + t219;
t317 = -pkin(2) - pkin(13);
t141 = (pkin(1) - t317) * (pkin(1) + t317) + t259;
t312 = pkin(13) - pkin(2);
t142 = (pkin(1) - t312) * (pkin(1) + t312) + t259;
t274 = t142 * t141;
t227 = sqrt(-t274);
t290 = t227 / 0.2e1;
t322 = -0.2e1 * pkin(1);
t220 = pkin(5) ^ 2;
t218 = pkin(7) ^ 2;
t226 = pkin(1) ^ 2;
t199 = sin(qJ(2));
t201 = sin(pkin(19));
t205 = cos(qJ(2));
t207 = cos(pkin(19));
t175 = t199 * t207 - t201 * t205;
t285 = pkin(7) * t175;
t258 = t285 * t322 + t226;
t159 = t218 + t258;
t255 = pkin(3) ^ 2 - pkin(8) ^ 2;
t151 = t159 + t255;
t166 = pkin(1) - t285;
t177 = t199 * t201 + t205 * t207;
t315 = pkin(7) + pkin(8);
t316 = pkin(7) - pkin(8);
t143 = (pkin(3) + t315) * (-pkin(3) + t316) + t258;
t144 = (-pkin(3) + t315) * (pkin(3) + t316) + t258;
t228 = sqrt(-t144 * t143);
t268 = t177 * t228;
t118 = -pkin(7) * t268 + t151 * t166;
t262 = t204 * t118;
t119 = pkin(7) * t151 * t177 + t166 * t228;
t263 = t198 * t119;
t231 = -t262 / 0.2e1 + t263 / 0.2e1;
t157 = 0.1e1 / t159;
t223 = 0.1e1 / pkin(3);
t270 = t157 * t223;
t100 = t231 * t270;
t188 = pkin(23) + pkin(22);
t185 = sin(t188);
t186 = cos(t188);
t261 = t204 * t119;
t264 = t198 * t118;
t232 = t264 / 0.2e1 + t261 / 0.2e1;
t99 = t232 * t270;
t78 = t100 * t185 - t186 * t99;
t307 = pkin(5) * t78;
t278 = -0.2e1 * pkin(4) * t307 + t220;
t314 = -pkin(9) - pkin(11);
t66 = (pkin(4) - t314) * (pkin(4) + t314) + t278;
t313 = pkin(11) - pkin(9);
t67 = (pkin(4) - t313) * (pkin(4) + t313) + t278;
t229 = sqrt(-t67 * t66);
t79 = t100 * t186 + t185 * t99;
t279 = t229 * t79;
t257 = pkin(9) ^ 2 - pkin(11) ^ 2;
t221 = pkin(4) ^ 2;
t74 = t221 + t278;
t68 = t74 + t257;
t75 = -pkin(4) + t307;
t40 = -pkin(5) * t279 - t68 * t75;
t321 = t40 / 0.2e1;
t70 = 0.1e1 / t74;
t320 = t70 / 0.2e1;
t171 = t177 * qJD(2);
t170 = t175 * qJD(2);
t252 = pkin(1) * pkin(7) * t171;
t275 = 0.2e1 * (t143 + t144) * t252 / t228;
t245 = -t275 / 0.2e1;
t230 = t170 * t228 + t177 * t245;
t92 = ((t166 * t322 - t151) * t171 + t230) * pkin(7);
t319 = -t92 / 0.2e1;
t250 = -0.2e1 * t171 * t177;
t269 = t171 * t228;
t93 = t166 * t275 / 0.2e1 + t218 * pkin(1) * t250 + (-t151 * t170 - t269) * pkin(7);
t318 = t93 / 0.2e1;
t200 = sin(qJ(1));
t206 = cos(qJ(1));
t176 = -t200 * V_base(4) + t206 * V_base(5);
t172 = qJD(2) - t176;
t193 = cos(pkin(23));
t243 = 0.1e1 / t159 ^ 2 * t252;
t265 = t193 * t119;
t266 = t193 * t118;
t189 = sin(pkin(23));
t267 = t189 * t119;
t277 = t118 * t189;
t298 = t189 / 0.2e1;
t96 = (-t266 / 0.2e1 + t267 / 0.2e1) * t270;
t95 = 0.1e1 / t96 ^ 2;
t97 = (t265 / 0.2e1 + t277 / 0.2e1) * t270;
t33 = t172 + (-((t193 * t319 + t298 * t93) * t157 + (-t266 + t267) * t243) * t97 * t95 + ((t193 * t318 + t298 * t92) * t157 + (t265 + t277) * t243) / t96) * t223 / (t95 * t97 ^ 2 + 0.1e1);
t311 = pkin(4) * t33;
t293 = -t204 / 0.2e1;
t51 = ((-t262 + t263) * t243 + (qJD(3) * t232 + t198 * t318 + t293 * t92) * t157) * t223;
t52 = ((-t261 - t264) * t243 + (qJD(3) * t231 + t198 * t319 + t293 * t93) * t157) * t223;
t45 = t185 * t51 + t186 * t52;
t310 = pkin(4) * t45;
t309 = pkin(4) * t79;
t306 = pkin(5) * t79;
t41 = -t229 * t75 + t306 * t68;
t308 = pkin(5) * t41;
t224 = pkin(2) ^ 2;
t256 = t219 + t226;
t248 = -pkin(13) ^ 2 + t256;
t149 = t163 + t224 + t248;
t162 = -pkin(1) - t287;
t174 = t192 * t198 - t196 * t204;
t286 = pkin(6) * t174;
t115 = -t149 * t162 - t227 * t286;
t305 = t115 / 0.2e1;
t116 = t149 * t286 - t162 * t227;
t304 = -t116 / 0.2e1;
t182 = pkin(14) * V_base(5) + V_base(1);
t183 = -pkin(14) * V_base(4) + V_base(2);
t153 = t206 * t182 + t200 * t183;
t160 = -pkin(16) * t176 + V_base(3);
t136 = -t153 * t205 - t160 * t199;
t303 = t136 / 0.2e1;
t178 = t200 * V_base(5) + t206 * V_base(4);
t187 = V_base(6) + qJD(1);
t145 = -t178 * t199 + t187 * t205;
t302 = t145 / 0.2e1;
t146 = -t178 * t205 - t187 * t199;
t301 = t146 / 0.2e1;
t148 = t224 - t248 - 0.2e1 * t251;
t300 = -t148 / 0.2e1;
t156 = t163 + t256;
t154 = 0.1e1 / t156;
t299 = t154 / 0.2e1;
t191 = sin(pkin(21));
t297 = t191 / 0.2e1;
t194 = cos(pkin(22));
t296 = -t194 / 0.2e1;
t195 = cos(pkin(21));
t295 = -t195 / 0.2e1;
t294 = t195 / 0.2e1;
t208 = cos(pkin(18));
t292 = t208 / 0.2e1;
t291 = -t227 / 0.2e1;
t165 = t174 * qJD(3);
t289 = pkin(1) * t165;
t288 = pkin(6) * t116;
t254 = pkin(5) * t310;
t284 = 0.2e1 * (t66 + t67) * t254 / t229;
t125 = pkin(1) * t172 + t136;
t135 = -t153 * t199 + t160 * t205;
t112 = t198 * t125 + t204 * t135;
t169 = qJD(3) + t172;
t106 = pkin(5) * t169 + t112;
t111 = -t125 * t204 + t135 * t198;
t69 = t74 - t257;
t76 = -pkin(4) * t78 + pkin(5);
t234 = -pkin(4) * t279 + t69 * t76;
t235 = t229 * t76 + t309 * t69;
t213 = 0.1e1 / pkin(11);
t282 = t213 * t70;
t27 = (t234 * t297 + t235 * t294) * t282;
t28 = (t234 * t295 + t235 * t297) * t282;
t17 = t27 * t106 + t28 * t111;
t60 = t97 * t125 + t96 * t135;
t61 = t96 * t125 - t97 * t135;
t215 = 0.1e1 / pkin(9);
t281 = t215 * t70;
t280 = t229 * t45;
t253 = pkin(6) * t289;
t276 = (t141 + t142) * t253 / t290;
t225 = 0.1e1 / pkin(2);
t273 = t154 * t225;
t202 = sin(pkin(18));
t272 = t157 * t202;
t217 = 0.1e1 / pkin(8);
t271 = t157 * t217;
t260 = 0.1e1 / pkin(13) * t225;
t83 = (t115 * t303 + t135 * t304) * t273;
t249 = -t284 / 0.2e1;
t71 = 0.1e1 / t74 ^ 2;
t247 = t71 * t254;
t246 = -t276 / 0.2e1;
t244 = t157 * t292;
t152 = -t200 * t182 + t183 * t206;
t240 = t202 * t243;
t239 = t208 * t243;
t16 = t106 * t28 - t111 * t27;
t126 = t145 * t198 - t146 * t204;
t127 = t145 * t204 + t146 * t198;
t19 = -t126 * t27 + t127 * t28;
t238 = t191 * t69 + t195 * t229;
t237 = t191 * t229 - t195 * t69;
t140 = -pkin(16) * t187 - t152;
t44 = -t185 * t52 + t186 * t51;
t233 = -t44 * t229 + t249 * t79;
t114 = 0.1e1 / t115 ^ 2;
t155 = 0.1e1 / t156 ^ 2;
t164 = t236 * qJD(3);
t48 = t172 + (-0.2e1 * ((-t165 * t149 - t164 * t227 + t174 * t246) * t299 + (t115 * t155 + t154 * t162) * t289) * t114 * t288 + ((t162 * t246 + (t149 * t164 - t165 * t227) * pkin(6)) * t299 + (-t154 * t174 * t219 + t155 * t288) * t289) / t305) * pkin(2) / (t114 * t116 ^ 2 + 0.1e1) * t156 * t225;
t129 = -pkin(1) * t146 + t140;
t113 = -pkin(5) * t127 + t129;
t209 = V_base(3) ^ 2;
t203 = cos(qJ(4));
t197 = sin(qJ(4));
t190 = sin(pkin(22));
t167 = pkin(1) * t175 - pkin(7);
t161 = pkin(15) * t176 + V_base(3);
t150 = t159 - t255;
t147 = 0.1e1 / t148 ^ 2;
t139 = t140 ^ 2;
t138 = -pkin(17) * t176 + t153;
t137 = pkin(15) * t187 - pkin(17) * t178 - t152;
t128 = t129 ^ 2;
t120 = pkin(1) * t150 * t177 - t167 * t228;
t117 = -pkin(1) * t268 - t150 * t167;
t102 = (t120 * t292 + t117 * t202 / 0.2e1) * t271;
t101 = (t117 * t292 - t202 * t120 / 0.2e1) * t271;
t98 = 0.1e1 / t101 ^ 2;
t94 = t167 * t245 + t226 * pkin(7) * t250 + (-t150 * t170 - t269) * pkin(1);
t91 = ((0.2e1 * pkin(7) * t167 - t150) * t171 + t230) * pkin(1);
t88 = (t115 * t301 + t145 * t304) * t273;
t87 = (t115 * t302 + t116 * t301) * t273;
t84 = -pkin(2) * t88 + t140;
t82 = (t116 * t303 + t135 * t305) * t273;
t73 = t101 * t187 - t102 * t178;
t72 = t101 * t178 + t102 * t187;
t65 = -t145 * t97 + t146 * t96;
t64 = t145 * t96 + t146 * t97;
t63 = t101 * t161 - t102 * t138;
t62 = t101 * t138 + t102 * t161;
t59 = (t290 * t87 + t300 * t88) * t260;
t58 = (t291 * t88 + t300 * t87) * t260;
t50 = (-t190 * t64 - t194 * t65) * pkin(4) + t129;
t47 = pkin(2) * t48 + t83;
t46 = (0.1e1 / t148 * t276 / 0.2e1 - 0.2e1 * t227 * t147 * t253) / (-t147 * t274 + 0.1e1) + t48;
t43 = (t290 * t82 + t300 * t47) * t260;
t42 = (t291 * t47 + t300 * t82) * t260;
t39 = 0.1e1 / t40 ^ 2;
t34 = ((t94 * t244 + t120 * t239 + t91 * t272 / 0.2e1 + t117 * t240) / t101 - (t91 * t244 + t117 * t239 - t94 * t272 / 0.2e1 - t120 * t240) * t102 * t98) / (t102 ^ 2 * t98 + 0.1e1) * t217 - t176;
t32 = t194 * t311 + t61;
t31 = t190 * t311 + t60;
t26 = (t40 * t296 - t190 * t41 / 0.2e1) * t281;
t25 = (t190 * t321 + t296 * t41) * t281;
t24 = 0.1e1 / t28 ^ 2;
t20 = t126 * t28 + t127 * t27;
t18 = qJD(4) - t19;
t15 = -t25 * t64 + t26 * t65;
t14 = t25 * t65 + t26 * t64;
t13 = t76 * t284 / 0.2e1 - 0.2e1 * t221 * t45 * t306 + (t44 * t69 - t280) * pkin(4);
t12 = ((-0.2e1 * pkin(5) * t76 - t69) * t45 + t233) * pkin(4);
t11 = -t25 * t31 + t26 * t32;
t10 = t25 * t32 + t26 * t31;
t9 = -pkin(10) * t19 - pkin(12) * t20 + t113;
t8 = t33 + (((t75 * t249 + (t44 * t68 - t280) * pkin(5)) * t320 + (-t220 * t70 * t79 + t308 * t71) * t310) / t321 - 0.2e1 * ((-t45 * t68 + t233) * t320 + (t40 * t71 + t70 * t75) * t310) * t39 * t308) * pkin(9) * t215 / (t39 * t41 ^ 2 + 0.1e1) * t74;
t7 = t169 + (((t12 * t297 + t13 * t294) * t70 + (-t237 * t309 + t238 * t76) * t247) / t28 - ((t12 * t295 + t13 * t297) * t70 + (t237 * t76 + t238 * t309) * t247) * t27 * t24) / (t24 * t27 ^ 2 + 0.1e1) * t213;
t6 = t197 * t7 + t20 * t203;
t5 = -t197 * t20 + t203 * t7;
t4 = pkin(12) * t7 + t17;
t3 = -pkin(10) * t7 - t16;
t2 = t197 * t9 + t203 * t4;
t1 = -t197 * t4 + t203 * t9;
t21 = (t129 * mrSges(8,2) - t61 * mrSges(8,3) + Ifges(8,4) * t65 + Ifges(8,1) * t64 / 0.2e1) * t64 + (t129 * mrSges(4,2) - t112 * mrSges(4,3) + Ifges(4,4) * t127 + Ifges(4,5) * t169 + Ifges(4,1) * t126 / 0.2e1) * t126 + (t152 * mrSges(2,1) - t153 * mrSges(2,2) + Ifges(2,3) * t187 / 0.2e1) * t187 + (V_base(3) * mrSges(2,2) - t152 * mrSges(2,3) + Ifges(2,5) * t187 + Ifges(2,1) * t178 / 0.2e1) * t178 + (-t113 * mrSges(5,1) + t17 * mrSges(5,3) + Ifges(5,4) * t20 + Ifges(5,6) * t7 + Ifges(5,2) * t19 / 0.2e1) * t19 + (t11 * mrSges(11,1) - t10 * mrSges(11,2) + Ifges(11,3) * t8 / 0.2e1) * t8 + (-t137 * mrSges(7,1) + t62 * mrSges(7,3) + Ifges(7,2) * t73 / 0.2e1) * t73 + (t113 * mrSges(5,2) - t16 * mrSges(5,3) + Ifges(5,5) * t7 + Ifges(5,1) * t20 / 0.2e1) * t20 + (t137 * mrSges(7,2) - t63 * mrSges(7,3) + Ifges(7,4) * t73 + Ifges(7,1) * t72 / 0.2e1) * t72 + (t43 * mrSges(10,1) - t42 * mrSges(10,2) + Ifges(10,5) * t58 + Ifges(10,6) * t59 + Ifges(10,3) * t46 / 0.2e1) * t46 + (t136 * mrSges(3,1) - t135 * mrSges(3,2) + Ifges(3,3) * t172 / 0.2e1) * t172 + (-t140 * mrSges(9,1) + t82 * mrSges(9,3) + Ifges(9,2) * t88 / 0.2e1) * t88 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t140 * mrSges(3,1) + t135 * mrSges(3,3) + Ifges(3,2) * t301 + Ifges(3,6) * t172) * t146 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t6 / 0.2e1) * t6 + (-t50 * mrSges(11,1) + t10 * mrSges(11,3) + Ifges(11,6) * t8 + Ifges(11,2) * t15 / 0.2e1) * t15 + (t83 * mrSges(9,1) - t82 * mrSges(9,2) + Ifges(9,5) * t87 + Ifges(9,6) * t88 + Ifges(9,3) * t48 / 0.2e1) * t48 + (-t84 * mrSges(10,1) + t42 * mrSges(10,3) + Ifges(10,2) * t59 / 0.2e1) * t59 + (t63 * mrSges(7,1) - t62 * mrSges(7,2) + Ifges(7,5) * t72 + Ifges(7,6) * t73 + Ifges(7,3) * t34 / 0.2e1) * t34 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t129 * mrSges(8,1) + t60 * mrSges(8,3) + Ifges(8,2) * t65 / 0.2e1) * t65 + (t61 * mrSges(8,1) - t60 * mrSges(8,2) + Ifges(8,5) * t64 + Ifges(8,6) * t65 + Ifges(8,3) * t33 / 0.2e1) * t33 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t6 + Ifges(6,6) * t5 + Ifges(6,3) * t18 / 0.2e1) * t18 + (t112 * mrSges(4,1) - t111 * mrSges(4,2) + Ifges(4,3) * t169 / 0.2e1) * t169 + (t140 * mrSges(3,2) - t136 * mrSges(3,3) + Ifges(3,1) * t302 + Ifges(3,4) * t146 + Ifges(3,5) * t172) * t145 + m(11) * (t10 ^ 2 + t11 ^ 2 + t50 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(2,1) + t153 * mrSges(2,3) + Ifges(2,4) * t178 + Ifges(2,6) * t187 + Ifges(2,2) * t176 / 0.2e1) * t176 + (-t129 * mrSges(4,1) + t111 * mrSges(4,3) + Ifges(4,6) * t169 + Ifges(4,2) * t127 / 0.2e1) * t127 + (t16 * mrSges(5,1) - t17 * mrSges(5,2) + Ifges(5,3) * t7 / 0.2e1) * t7 + (t84 * mrSges(10,2) - t43 * mrSges(10,3) + Ifges(10,4) * t59 + Ifges(10,1) * t58 / 0.2e1) * t58 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t209) / 0.2e1 + m(2) * (t152 ^ 2 + t153 ^ 2 + t209) / 0.2e1 + m(3) * (t135 ^ 2 + t136 ^ 2 + t139) / 0.2e1 + m(9) * (t82 ^ 2 + t83 ^ 2 + t139) / 0.2e1 + (t140 * mrSges(9,2) - t83 * mrSges(9,3) + Ifges(9,4) * t88 + Ifges(9,1) * t87 / 0.2e1) * t87 + m(4) * (t111 ^ 2 + t112 ^ 2 + t128) / 0.2e1 + m(8) * (t60 ^ 2 + t61 ^ 2 + t128) / 0.2e1 + m(7) * (t137 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t16 ^ 2 + t17 ^ 2) / 0.2e1 + m(10) * (t42 ^ 2 + t43 ^ 2 + t84 ^ 2) / 0.2e1 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t6 + Ifges(6,2) * t5 / 0.2e1) * t5 + (t50 * mrSges(11,2) - t11 * mrSges(11,3) + Ifges(11,4) * t15 + Ifges(11,5) * t8 + Ifges(11,1) * t14 / 0.2e1) * t14;
T = t21;
