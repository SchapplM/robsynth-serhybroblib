% Calculate joint inertia matrix for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1TE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_inertiaJ_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1TE_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1TE_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:38:58
% EndTime: 2020-04-17 15:48:59
% DurationCPUTime: 112.36s
% Computational Cost: add. (3195589->322), mult. (4826472->585), div. (212122->22), fcn. (3059918->22), ass. (0->250)
t152 = sin(qJ(4));
t157 = cos(qJ(4));
t309 = Ifges(6,5) * t152 + Ifges(6,6) * t157;
t308 = 2 * pkin(3);
t168 = pkin(4) ^ 2;
t167 = pkin(5) ^ 2;
t172 = pkin(1) ^ 2;
t154 = sin(qJ(2));
t155 = sin(pkin(16));
t159 = cos(qJ(2));
t278 = cos(pkin(16));
t134 = t154 * t155 - t159 * t278;
t276 = pkin(5) * t134;
t304 = -2 * pkin(1);
t225 = t276 * t304 + t172;
t126 = t167 + t225;
t124 = 0.1e1 / t126;
t171 = 0.1e1 / pkin(2);
t228 = t124 * t171;
t221 = pkin(2) ^ 2 - pkin(6) ^ 2;
t122 = t126 + t221;
t127 = pkin(1) - t276;
t295 = pkin(5) + pkin(6);
t296 = pkin(5) - pkin(6);
t119 = (pkin(2) + t295) * (-pkin(2) + t296) + t225;
t120 = (-pkin(2) + t295) * (pkin(2) + t296) + t225;
t173 = sqrt(-t120 * t119);
t136 = t154 * t278 + t159 * t155;
t275 = pkin(5) * t136;
t114 = t122 * t275 + t127 * t173;
t153 = sin(qJ(3));
t233 = t114 * t153;
t226 = t136 * t173;
t113 = -pkin(5) * t226 + t122 * t127;
t158 = cos(qJ(3));
t236 = t113 * t158;
t107 = (-t236 / 0.2e1 + t233 / 0.2e1) * t228;
t232 = t114 * t158;
t237 = t113 * t153;
t108 = (t232 / 0.2e1 + t237 / 0.2e1) * t228;
t143 = pkin(18) + pkin(19);
t141 = sin(t143);
t142 = cos(t143);
t90 = t107 * t142 + t108 * t141;
t287 = pkin(4) * t90;
t240 = t287 * t308 + t168;
t294 = (-pkin(8) - pkin(10));
t79 = ((pkin(3) - t294) * (pkin(3) + t294)) + t240;
t293 = (pkin(10) - pkin(8));
t80 = ((pkin(3) - t293) * (pkin(3) + t293)) + t240;
t174 = sqrt(-t80 * t79);
t187 = t107 * t141 - t108 * t142;
t306 = t174 * t187;
t305 = -t158 * mrSges(4,1) + t153 * mrSges(4,2);
t303 = -0.2e1 * t136 ^ 2;
t147 = cos(pkin(17));
t169 = pkin(3) ^ 2;
t85 = t169 + t240;
t84 = 0.1e1 / t85 ^ 2;
t289 = pkin(4) * t84;
t219 = pkin(3) * t289;
t206 = t147 * t219;
t218 = pkin(1) * t275;
t204 = 0.1e1 / t126 ^ 2 * t218;
t281 = t153 / 0.2e1;
t231 = 0.2e1 / t173 * (t119 + t120) * t218;
t100 = t127 * t231 / 0.2e1 + t167 * pkin(1) * t303 + (-t122 * t134 - t226) * pkin(5);
t286 = t100 / 0.2e1;
t210 = -t231 / 0.2e1;
t227 = t134 * t173;
t98 = (t227 + (t127 * t304 - t122 + t210) * t136) * pkin(5);
t77 = ((t158 * t286 + t281 * t98) * t124 + (t232 + t237) * t204) * t171;
t297 = -t98 / 0.2e1;
t78 = ((t100 * t281 + t158 * t297) * t124 + (t233 - t236) * t204) * t171;
t71 = -t141 * t77 - t142 * t78;
t198 = t71 * t206;
t146 = sin(pkin(17));
t207 = t146 * t219;
t200 = t71 * t207;
t83 = 0.1e1 / t85;
t250 = t147 * t83;
t212 = t250 / 0.2e1;
t213 = -t250 / 0.2e1;
t285 = t146 / 0.2e1;
t214 = t83 * t285;
t162 = 0.1e1 / pkin(10);
t245 = t162 * t83;
t222 = pkin(8) ^ 2 - pkin(10) ^ 2;
t82 = t85 - t222;
t87 = pkin(3) * t90 + pkin(4);
t66 = -pkin(3) * t306 + t82 * t87;
t68 = pkin(3) * t187 * t82 + t174 * t87;
t49 = (-t66 * t147 / 0.2e1 + t68 * t285) * t245;
t176 = t49 ^ 2;
t47 = 0.1e1 / t176;
t50 = (t68 * t147 / 0.2e1 + t66 * t285) * t245;
t48 = t50 ^ 2;
t246 = t162 / (t47 * t48 + 0.1e1);
t259 = t47 * t50;
t74 = 0.1e1 / t174;
t299 = -t74 / 0.2e1;
t215 = t187 * t299;
t186 = pkin(4) * (t79 + t80) * t308;
t54 = t71 * t186;
t70 = t141 * t78 - t142 * t77;
t181 = -t70 * t174 + t215 * t54;
t211 = -0.2e1 * pkin(4) * t87 - t82;
t34 = (t211 * t71 + t181) * pkin(3);
t288 = pkin(4) * t187;
t208 = -0.2e1 * t169 * t288;
t216 = t74 * t87 / 0.2e1;
t243 = t174 * t71;
t35 = t54 * t216 + t71 * t208 + (t70 * t82 - t243) * pkin(3);
t46 = 0.1e1 / t49;
t12 = 0.1e1 + ((t198 * t68 + t200 * t66 + t212 * t35 + t214 * t34) * t46 - (-t198 * t66 + t200 * t68 + t213 * t34 + t214 * t35) * t259) * t246;
t140 = -pkin(1) * t158 + pkin(4);
t277 = pkin(1) * t153;
t41 = t140 * t49 + t277 * t50;
t9 = -pkin(9) * t12 - t41;
t301 = m(6) * t9;
t133 = t153 * t154 - t158 * t159;
t135 = -t153 * t159 - t154 * t158;
t39 = t133 * t49 - t135 * t50;
t300 = -t39 / 0.2e1;
t298 = t83 / 0.2e1;
t148 = sin(pkin(19));
t235 = t114 * t148;
t150 = cos(pkin(19));
t238 = t113 * t150;
t104 = (-t238 / 0.2e1 + t235 / 0.2e1) * t228;
t175 = t104 ^ 2;
t101 = 0.1e1 / t175;
t234 = t114 * t150;
t239 = t113 * t148;
t105 = (t234 / 0.2e1 + t239 / 0.2e1) * t228;
t102 = t105 ^ 2;
t284 = t148 / 0.2e1;
t57 = 0.1e1 + (-((t100 * t284 + t150 * t297) * t124 + (t235 - t238) * t204) * t105 * t101 + ((t150 * t286 + t284 * t98) * t124 + (t234 + t239) * t204) / t104) * t171 / (t101 * t102 + 0.1e1);
t292 = pkin(3) * t57;
t291 = pkin(4) * t49;
t290 = pkin(4) * t50;
t151 = cos(pkin(18));
t283 = -t151 / 0.2e1;
t282 = -t152 / 0.2e1;
t280 = t157 / 0.2e1;
t160 = cos(pkin(15));
t279 = t160 / 0.2e1;
t149 = sin(pkin(18));
t164 = 0.1e1 / pkin(8);
t244 = t164 * t83;
t81 = t85 + t222;
t86 = -pkin(3) - t287;
t65 = -pkin(4) * t306 - t81 * t86;
t67 = -t174 * t86 + t288 * t81;
t51 = (t149 * t65 / 0.2e1 + t67 * t283) * t244;
t52 = (t65 * t283 - t149 * t67 / 0.2e1) * t244;
t55 = pkin(1) * t104 + t151 * t292;
t56 = pkin(1) * t105 + t149 * t292;
t31 = -t51 * t56 + t52 * t55;
t274 = mrSges(9,1) * t31;
t273 = mrSges(5,3) * t39;
t272 = Ifges(6,6) * t39;
t197 = t187 * t206;
t199 = t187 * t207;
t69 = t187 * t186;
t180 = -t90 * t174 + t215 * t69;
t43 = (t187 * t211 + t180) * pkin(3);
t44 = t69 * t216 + t187 * t208 + (t82 * t90 - t306) * pkin(3);
t18 = 0.1e1 + ((t197 * t68 + t199 * t66 + t212 * t44 + t214 * t43) * t46 - (-t197 * t66 + t199 * t68 + t213 * t43 + t214 * t44) * t259) * t246;
t271 = Ifges(5,3) * t18;
t178 = pkin(3) * (-t168 * t187 * t83 + t289 * t67);
t64 = 0.1e1 / t65 ^ 2;
t203 = pkin(8) * t164 / (t64 * t67 ^ 2 + 0.1e1) * t85;
t183 = pkin(4) * t64 * t67 * t203;
t185 = pkin(3) * (t65 * t84 + t83 * t86);
t196 = 0.1e1 / t65 * t203;
t217 = t86 * t299;
t19 = 0.2e1 * ((t54 * t217 + (t70 * t81 - t243) * pkin(4)) * t298 + t71 * t178) * t196 - 0.2e1 * ((-t71 * t81 + t181) * t298 + t71 * t185) * t183 + t57;
t270 = Ifges(9,3) * t19;
t253 = Ifges(6,4) * t157;
t191 = Ifges(6,1) * t152 + t253;
t3 = t191 * t12;
t269 = t152 * t3;
t6 = t191 * t18;
t268 = t152 * t6;
t254 = Ifges(6,4) * t152;
t190 = Ifges(6,2) * t157 + t254;
t2 = t190 * t12;
t267 = t157 * t2;
t5 = t190 * t18;
t266 = t157 * t5;
t40 = t50 * t140 - t277 * t49;
t261 = t40 * mrSges(5,2);
t260 = t41 * mrSges(5,1);
t258 = t309 * t12;
t257 = t309 * t18;
t38 = t133 * t50 + t135 * t49;
t248 = t157 * t38;
t256 = Ifges(6,5) * t248 - Ifges(6,3) * t39;
t255 = mrSges(6,1) * t157;
t156 = sin(pkin(15));
t230 = t124 * t156;
t166 = 0.1e1 / pkin(6);
t229 = t124 * t166;
t224 = Ifges(4,5) * t135 + Ifges(4,6) * t133;
t223 = t152 ^ 2 + t157 ^ 2;
t209 = t124 * t279;
t138 = -pkin(1) * t159 - pkin(13);
t205 = mrSges(6,3) * t223;
t201 = 0.2e1 * t223;
t121 = t126 - t221;
t128 = pkin(1) * t134 - pkin(5);
t112 = -pkin(1) * t226 - t121 * t128;
t195 = t112 * t204;
t115 = pkin(1) * t136 * t121 - t128 * t173;
t194 = t115 * t204;
t193 = mrSges(6,2) * t152 - t255;
t192 = mrSges(6,1) * t152 + mrSges(6,2) * t157;
t27 = -mrSges(6,3) * t152 * t38 + mrSges(6,2) * t39;
t28 = -mrSges(6,1) * t39 - mrSges(6,3) * t248;
t189 = -t152 * t28 + t157 * t27;
t188 = mrSges(6,3) * t201;
t123 = -pkin(4) * t133 + t138;
t184 = (mrSges(5,1) * t49 - mrSges(5,2) * t50) * pkin(4);
t182 = m(6) * t201 / 0.2e1;
t93 = t104 * t154 + t105 * t159;
t94 = t104 * t159 - t105 * t154;
t32 = t51 * t94 + t52 * t93;
t33 = -t51 * t93 + t52 * t94;
t179 = t32 * Ifges(9,5) + t33 * Ifges(9,6);
t21 = -t272 + (-Ifges(6,2) * t152 + t253) * t38;
t22 = -Ifges(6,5) * t39 + (Ifges(6,1) * t157 - t254) * t38;
t177 = t21 * t280 + t152 * t22 / 0.2e1 + t38 * Ifges(5,5) + t39 * Ifges(5,6);
t109 = (t115 * t279 - t112 * t156 / 0.2e1) * t229;
t106 = (t112 * t279 + t115 * t156 / 0.2e1) * t229;
t103 = 0.1e1 / t106 ^ 2;
t99 = t128 * t210 + t172 * pkin(5) * t303 + (-t121 * t134 - t226) * pkin(1);
t97 = (t227 + (0.2e1 * t128 * pkin(5) - t121 + t210) * t136) * pkin(1);
t75 = (-t149 * t93 - t151 * t94) * pkin(3) + t138;
t58 = ((t99 * t209 + t160 * t194 - t97 * t230 / 0.2e1 - t156 * t195) / t106 - (t97 * t209 + t160 * t195 + t99 * t230 / 0.2e1 + t156 * t194) * t109 * t103) / (t103 * t109 ^ 2 + 0.1e1) * t166;
t30 = t51 * t55 + t52 * t56;
t26 = t192 * t38;
t25 = -pkin(9) * t39 - pkin(11) * t38 + t123;
t20 = 0.2e1 * ((t69 * t217 + (t90 * t81 - t306) * pkin(4)) * t298 + t187 * t178) * t196 - 0.2e1 * ((-t187 * t81 + t180) * t298 + t187 * t185) * t183;
t15 = pkin(11) * t18 + t290;
t14 = -pkin(9) * t18 - t291;
t8 = pkin(11) * t12 + t40;
t4 = t193 * t18;
t1 = t193 * t12;
t7 = [(m(4) + m(8)) * t138 ^ 2 + m(7) * pkin(7) ^ 2 + m(3) * pkin(13) ^ 2 + t154 * (Ifges(3,1) * t154 + Ifges(3,4) * t159) - 0.2e1 * pkin(13) * (-mrSges(3,1) * t159 + mrSges(3,2) * t154) + t159 * (Ifges(3,4) * t154 + Ifges(3,2) * t159) + m(5) * t123 ^ 2 + Ifges(4,1) * t135 ^ 2 + 0.2e1 * pkin(7) * (-mrSges(7,1) * t106 + mrSges(7,2) * t109) + t106 * (Ifges(7,4) * t109 + Ifges(7,2) * t106) + t109 * (Ifges(7,1) * t109 + Ifges(7,4) * t106) + Ifges(8,2) * t94 ^ 2 + m(9) * t75 ^ 2 + 0.2e1 * t75 * (-mrSges(9,1) * t33 + mrSges(9,2) * t32) + t33 * (Ifges(9,4) * t32 + Ifges(9,2) * t33) + t32 * (Ifges(9,1) * t32 + Ifges(9,4) * t33) + (-0.2e1 * mrSges(5,1) * t123 + Ifges(5,2) * t39 - t256) * t39 + (0.2e1 * mrSges(5,2) * t123 + Ifges(5,1) * t38 + 0.2e1 * Ifges(5,4) * t39 + t157 * t22 + (-t21 + t272) * t152) * t38 + Ifges(2,3) + (Ifges(8,1) * t93 + 0.2e1 * Ifges(8,4) * t94) * t93 + (0.2e1 * Ifges(4,4) * t135 + Ifges(4,2) * t133) * t133 + 0.2e1 * (-mrSges(4,1) * t133 - mrSges(8,1) * t94 + mrSges(4,2) * t135 + mrSges(8,2) * t93) * t138 + (0.2e1 * t152 * t27 + 0.2e1 * t157 * t28 + t182 * t25) * t25; t9 * t26 + t258 * t300 + Ifges(3,5) * t154 + Ifges(3,6) * t159 + t40 * t273 + t189 * t8 + (t30 * t33 - t31 * t32) * mrSges(9,3) + (-t41 * mrSges(5,3) + t2 * t282 + t280 * t3) * t38 + t179 * t19 + ((-t104 * t93 + t105 * t94) * mrSges(8,3) + (-t133 * t153 + t135 * t158) * mrSges(4,3)) * pkin(1) + t177 * t12 + t58 * (Ifges(7,5) * t109 + Ifges(7,6) * t106) + t57 * (Ifges(8,5) * t93 + Ifges(8,6) * t94) + t224; t57 ^ 2 * Ifges(8,3) + t58 ^ 2 * Ifges(7,3) + (t40 ^ 2 + t41 ^ 2) * m(5) + m(9) * t31 ^ 2 + Ifges(3,3) + Ifges(4,3) + (0.2e1 * t1 + t301) * t9 + t8 ^ 2 * t182 + (t270 + 0.2e1 * t274) * t19 + (m(9) * t30 - 0.2e1 * mrSges(9,2) * t19) * t30 + (m(8) * (t102 + t175) + m(4) * (t153 ^ 2 + t158 ^ 2)) * t172 + 0.2e1 * ((mrSges(8,1) * t104 - mrSges(8,2) * t105) * t57 + t305) * pkin(1) + (t12 * Ifges(5,3) + t188 * t8 + 0.2e1 * t260 - 0.2e1 * t261 + t267 + t269) * t12; t14 * t26 + t257 * t300 + t273 * t290 + t189 * t15 + (-mrSges(5,3) * t291 + t280 * t6 + t282 * t5) * t38 + t179 * t20 + t177 * t18 + t224; t9 * t4 + Ifges(4,3) + t8 * t15 * t182 + (t1 + t301) * t14 + (t40 * t50 + t41 * t49) * pkin(4) * m(5) + t305 * pkin(1) + (-mrSges(9,2) * t30 + t270 + t274) * t20 + (t267 / 0.2e1 + t269 / 0.2e1 - t261 + t260 + t8 * t205) * t18 + (t271 + t266 / 0.2e1 + t268 / 0.2e1 + t15 * t205 + t184) * t12; t20 ^ 2 * Ifges(9,3) + Ifges(4,3) + (t176 + t48) * t168 * m(5) + (m(6) * t14 + 0.2e1 * t4) * t14 + t15 ^ 2 * t182 + (t15 * t188 + 0.2e1 * t184 + t266 + t268 + t271) * t18; t25 * t255 + (-mrSges(6,2) * t25 - Ifges(6,6) * t38) * t152 + t256; -t192 * t8 + t258; -t15 * t192 + t257; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
