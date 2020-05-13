% Calculate joint inertia matrix for
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1DE2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_inertiaJ_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE2_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1DE2_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:59:41
% EndTime: 2020-04-19 20:10:15
% DurationCPUTime: 127.18s
% Computational Cost: add. (3616052->322), mult. (5457796->584), div. (242190->22), fcn. (3459630->34), ass. (0->259)
t158 = sin(qJ(4));
t163 = cos(qJ(4));
t315 = Ifges(6,5) * t158 + Ifges(6,6) * t163;
t314 = 2 * pkin(3);
t174 = pkin(4) ^ 2;
t173 = pkin(5) ^ 2;
t178 = pkin(1) ^ 2;
t160 = sin(qJ(2));
t161 = sin(pkin(16));
t165 = cos(qJ(2));
t283 = cos(pkin(16));
t140 = t160 * t161 - t165 * t283;
t282 = pkin(5) * t140;
t310 = -2 * pkin(1);
t231 = t282 * t310 + t178;
t133 = t173 + t231;
t131 = 0.1e1 / t133;
t177 = 0.1e1 / pkin(2);
t234 = t131 * t177;
t227 = pkin(2) ^ 2 - pkin(6) ^ 2;
t129 = t133 + t227;
t134 = pkin(1) - t282;
t302 = -pkin(2) - pkin(6);
t126 = (pkin(5) - t302) * (pkin(5) + t302) + t231;
t301 = -pkin(6) + pkin(2);
t127 = (pkin(5) - t301) * (pkin(5) + t301) + t231;
t179 = sqrt(-t127 * t126);
t141 = t160 * t283 + t161 * t165;
t281 = pkin(5) * t141;
t121 = t129 * t281 + t134 * t179;
t159 = sin(qJ(3));
t239 = t121 * t159;
t232 = t141 * t179;
t120 = -pkin(5) * t232 + t129 * t134;
t164 = cos(qJ(3));
t242 = t120 * t164;
t114 = (-t242 / 0.2e1 + t239 / 0.2e1) * t234;
t238 = t121 * t164;
t243 = t120 * t159;
t115 = (t238 / 0.2e1 + t243 / 0.2e1) * t234;
t149 = pkin(18) + pkin(19);
t147 = sin(t149);
t148 = cos(t149);
t94 = t114 * t148 + t115 * t147;
t292 = t94 * pkin(4);
t246 = t292 * t314 + t174;
t300 = (-pkin(8) - pkin(10));
t83 = ((pkin(3) - t300) * (pkin(3) + t300)) + t246;
t299 = (-pkin(8) + pkin(10));
t84 = ((pkin(3) - t299) * (pkin(3) + t299)) + t246;
t180 = sqrt(-t84 * t83);
t192 = t114 * t147 - t115 * t148;
t312 = t180 * t192;
t311 = -t164 * mrSges(4,1) + t159 * mrSges(4,2);
t144 = -t165 * pkin(1) - pkin(13);
t191 = t159 * t160 - t164 * t165;
t130 = -pkin(4) * t191 + t144;
t309 = 0.2e1 * t130;
t308 = -0.2e1 * t141 ^ 2;
t307 = 0.2e1 * t144;
t153 = cos(pkin(17));
t175 = pkin(3) ^ 2;
t89 = t175 + t246;
t88 = 0.1e1 / t89 ^ 2;
t294 = pkin(4) * t88;
t225 = pkin(3) * t294;
t211 = t153 * t225;
t224 = pkin(1) * t281;
t237 = 0.2e1 / t179 * (t126 + t127) * t224;
t215 = -t237 / 0.2e1;
t233 = t140 * t179;
t106 = (t233 + (t134 * t310 - t129 + t215) * t141) * pkin(5);
t209 = 0.1e1 / t133 ^ 2 * t224;
t286 = t159 / 0.2e1;
t108 = t134 * t237 / 0.2e1 + t173 * pkin(1) * t308 + (-t129 * t140 - t232) * pkin(5);
t290 = t108 / 0.2e1;
t81 = ((t106 * t286 + t164 * t290) * t131 + (t238 + t243) * t209) * t177;
t291 = -t106 / 0.2e1;
t82 = ((t108 * t286 + t164 * t291) * t131 + (t239 - t242) * t209) * t177;
t74 = -t147 * t81 - t148 * t82;
t203 = t74 * t211;
t152 = sin(pkin(17));
t212 = t152 * t225;
t205 = t74 * t212;
t87 = 0.1e1 / t89;
t255 = t153 * t87;
t218 = t255 / 0.2e1;
t219 = -t255 / 0.2e1;
t289 = t152 / 0.2e1;
t220 = t87 * t289;
t168 = 0.1e1 / pkin(10);
t250 = t168 * t87;
t228 = pkin(8) ^ 2 - pkin(10) ^ 2;
t86 = t89 - t228;
t91 = pkin(3) * t94 + pkin(4);
t68 = -pkin(3) * t312 + t86 * t91;
t70 = pkin(3) * t192 * t86 + t180 * t91;
t55 = (-t68 * t153 / 0.2e1 + t70 * t289) * t250;
t54 = 0.1e1 / t55 ^ 2;
t56 = (t70 * t153 / 0.2e1 + t68 * t289) * t250;
t251 = t168 / (t54 * t56 ^ 2 + 0.1e1);
t264 = t54 * t56;
t77 = 0.1e1 / t180;
t304 = -t77 / 0.2e1;
t221 = t192 * t304;
t190 = pkin(4) * (t83 + t84) * t314;
t60 = t74 * t190;
t73 = t147 * t82 - t148 * t81;
t185 = -t180 * t73 + t221 * t60;
t216 = -0.2e1 * pkin(4) * t91 - t86;
t40 = (t216 * t74 + t185) * pkin(3);
t293 = pkin(4) * t192;
t213 = -0.2e1 * t175 * t293;
t222 = t77 * t91 / 0.2e1;
t249 = t180 * t74;
t41 = t60 * t222 + t74 * t213 + (t73 * t86 - t249) * pkin(3);
t53 = 0.1e1 / t55;
t12 = 0.1e1 + ((t203 * t70 + t205 * t68 + t218 * t41 + t220 * t40) * t53 - (-t203 * t68 + t205 * t70 + t219 * t40 + t220 * t41) * t264) * t251;
t146 = -pkin(1) * t164 + pkin(4);
t273 = t159 * pkin(1);
t49 = atan2(t56, t55);
t46 = sin(t49);
t47 = cos(t49);
t38 = t146 * t47 + t273 * t46;
t8 = -pkin(9) * t12 - t38;
t306 = m(6) * t8;
t142 = -t159 * t165 - t160 * t164;
t36 = t142 * t46 - t191 * t47;
t305 = t36 / 0.2e1;
t303 = t87 / 0.2e1;
t170 = 0.1e1 / pkin(8);
t217 = t170 * t303;
t85 = t89 + t228;
t90 = -pkin(3) - t292;
t67 = -pkin(4) * t312 - t85 * t90;
t69 = -t180 * t90 + t293 * t85;
t58 = atan2(t69 * t217, t67 * t217);
t298 = sin(t58);
t154 = sin(pkin(19));
t241 = t121 * t154;
t156 = cos(pkin(19));
t244 = t120 * t156;
t111 = (-t244 / 0.2e1 + t241 / 0.2e1) * t234;
t109 = 0.1e1 / t111 ^ 2;
t240 = t121 * t156;
t245 = t120 * t154;
t112 = (t240 / 0.2e1 + t245 / 0.2e1) * t234;
t288 = t154 / 0.2e1;
t63 = 0.1e1 + (-((t108 * t288 + t156 * t291) * t131 + (t241 - t244) * t209) * t112 * t109 + ((t106 * t288 + t156 * t290) * t131 + (t240 + t245) * t209) / t111) * t177 / (t109 * t112 ^ 2 + 0.1e1);
t297 = pkin(3) * t63;
t296 = pkin(4) * t46;
t295 = pkin(4) * t47;
t287 = -t158 / 0.2e1;
t285 = t163 / 0.2e1;
t166 = cos(pkin(15));
t284 = t166 / 0.2e1;
t155 = sin(pkin(18));
t157 = cos(pkin(18));
t57 = cos(t58);
t44 = t155 * t57 - t157 * t298;
t45 = -t155 * t298 - t157 * t57;
t103 = atan2(t112, t111);
t97 = sin(t103);
t61 = pkin(1) * t97 + t155 * t297;
t98 = cos(t103);
t62 = pkin(1) * t98 + t157 * t297;
t31 = t44 * t62 + t45 * t61;
t280 = mrSges(9,2) * t31;
t279 = mrSges(5,3) * t36;
t278 = Ifges(6,6) * t36;
t202 = t192 * t211;
t204 = t192 * t212;
t71 = t192 * t190;
t184 = -t180 * t94 + t221 * t71;
t50 = (t192 * t216 + t184) * pkin(3);
t51 = t71 * t222 + t192 * t213 + (t86 * t94 - t312) * pkin(3);
t18 = 0.1e1 + ((t202 * t70 + t204 * t68 + t218 * t51 + t220 * t50) * t53 - (-t202 * t68 + t204 * t70 + t219 * t50 + t220 * t51) * t264) * t251;
t277 = Ifges(5,3) * t18;
t182 = pkin(3) * (-t174 * t192 * t87 + t294 * t69);
t66 = 0.1e1 / t67 ^ 2;
t208 = pkin(8) * t170 / (t66 * t69 ^ 2 + 0.1e1) * t89;
t187 = pkin(4) * t66 * t69 * t208;
t189 = pkin(3) * (t67 * t88 + t87 * t90);
t201 = 0.1e1 / t67 * t208;
t223 = t90 * t304;
t21 = 0.2e1 * ((t60 * t223 + (t73 * t85 - t249) * pkin(4)) * t303 + t74 * t182) * t201 - 0.2e1 * ((-t74 * t85 + t185) * t303 + t74 * t189) * t187 + t63;
t276 = Ifges(9,3) * t21;
t258 = Ifges(6,4) * t163;
t196 = Ifges(6,1) * t158 + t258;
t3 = t196 * t12;
t275 = t158 * t3;
t6 = t196 * t18;
t274 = t158 * t6;
t259 = Ifges(6,4) * t158;
t195 = Ifges(6,2) * t163 + t259;
t2 = t195 * t12;
t272 = t163 * t2;
t5 = t195 * t18;
t271 = t163 * t5;
t266 = t38 * mrSges(5,1);
t39 = t146 * t46 - t273 * t47;
t265 = t39 * mrSges(5,2);
t263 = t315 * t12;
t262 = t315 * t18;
t37 = t142 * t47 + t191 * t46;
t253 = t163 * t37;
t261 = Ifges(6,5) * t253 + Ifges(6,3) * t36;
t260 = mrSges(6,1) * t163;
t162 = sin(pkin(15));
t236 = t131 * t162;
t172 = 0.1e1 / pkin(6);
t235 = t131 * t172;
t230 = Ifges(4,5) * t142 + Ifges(4,6) * t191;
t229 = t158 ^ 2 + t163 ^ 2;
t214 = t131 * t284;
t210 = mrSges(6,3) * t229;
t206 = 0.2e1 * t229;
t128 = t133 - t227;
t135 = pkin(1) * t140 - pkin(5);
t119 = -pkin(1) * t232 - t128 * t135;
t200 = t119 * t209;
t122 = pkin(1) * t141 * t128 - t135 * t179;
t199 = t122 * t209;
t198 = mrSges(6,2) * t158 - t260;
t197 = mrSges(6,1) * t158 + mrSges(6,2) * t163;
t26 = -mrSges(6,3) * t158 * t37 - mrSges(6,2) * t36;
t27 = mrSges(6,1) * t36 - mrSges(6,3) * t253;
t194 = -t158 * t27 + t163 * t26;
t193 = mrSges(6,3) * t206;
t188 = (mrSges(5,1) * t47 - mrSges(5,2) * t46) * pkin(4);
t186 = m(6) * t206 / 0.2e1;
t78 = -t160 * t97 + t165 * t98;
t79 = t160 * t98 + t165 * t97;
t32 = -t44 * t79 + t45 * t78;
t33 = t44 * t78 + t45 * t79;
t183 = t33 * Ifges(9,5) + t32 * Ifges(9,6);
t19 = t278 + (-Ifges(6,2) * t158 + t258) * t37;
t20 = Ifges(6,5) * t36 + (Ifges(6,1) * t163 - t259) * t37;
t181 = t37 * Ifges(5,5) - t36 * Ifges(5,6) + t19 * t285 + t158 * t20 / 0.2e1;
t116 = (t122 * t284 - t119 * t162 / 0.2e1) * t235;
t113 = (t119 * t284 + t122 * t162 / 0.2e1) * t235;
t110 = 0.1e1 / t113 ^ 2;
t107 = t135 * t215 + t178 * pkin(5) * t308 + (-t128 * t140 - t232) * pkin(1);
t105 = (t233 + (0.2e1 * pkin(5) * t135 - t128 + t215) * t141) * pkin(1);
t104 = atan2(t116, t113);
t101 = cos(t104);
t100 = sin(t104);
t72 = (-t155 * t79 - t157 * t78) * pkin(3) + t144;
t64 = ((t107 * t214 + t166 * t199 - t105 * t236 / 0.2e1 - t162 * t200) / t113 - (t105 * t214 + t166 * t200 + t107 * t236 / 0.2e1 + t162 * t199) * t116 * t110) / (t110 * t116 ^ 2 + 0.1e1) * t172;
t30 = -t44 * t61 + t45 * t62;
t29 = 0.2e1 * ((t71 * t223 + (t85 * t94 - t312) * pkin(4)) * t303 + t192 * t182) * t201 - 0.2e1 * ((-t192 * t85 + t184) * t303 + t192 * t189) * t187;
t25 = t197 * t37;
t24 = t36 * pkin(9) - t37 * pkin(11) + t130;
t15 = -pkin(9) * t18 - t295;
t14 = pkin(11) * t18 + t296;
t9 = pkin(11) * t12 + t39;
t4 = t198 * t18;
t1 = t198 * t12;
t7 = [m(3) * pkin(13) ^ 2 + m(7) * pkin(7) ^ 2 - 0.2e1 * pkin(13) * (-mrSges(3,1) * t165 + mrSges(3,2) * t160) + t165 * (Ifges(3,4) * t160 + Ifges(3,2) * t165) + t160 * (Ifges(3,1) * t160 + Ifges(3,4) * t165) + m(5) * t130 ^ 2 + 0.2e1 * pkin(7) * (-mrSges(7,1) * t101 + mrSges(7,2) * t100) + t100 * (Ifges(7,1) * t100 + Ifges(7,4) * t101) + t101 * (Ifges(7,4) * t100 + Ifges(7,2) * t101) + (mrSges(5,2) * t309 + Ifges(5,1) * t37 - 0.2e1 * Ifges(5,4) * t36 + t163 * t20 + (-t19 - t278) * t158) * t37 + (mrSges(5,1) * t309 + Ifges(5,2) * t36 + t261) * t36 + (m(8) + m(4)) * t144 ^ 2 + t33 * (Ifges(9,1) * t33 + Ifges(9,4) * t32) + t32 * (Ifges(9,4) * t33 + Ifges(9,2) * t32) + m(9) * t72 ^ 2 + 0.2e1 * t72 * (-mrSges(9,1) * t32 + mrSges(9,2) * t33) + Ifges(2,3) + (mrSges(8,2) * t307 + Ifges(8,1) * t79) * t79 + (-mrSges(8,1) * t307 + 0.2e1 * Ifges(8,4) * t79 + Ifges(8,2) * t78) * t78 + (0.2e1 * t158 * t26 + 0.2e1 * t163 * t27 + t186 * t24) * t24 + (-mrSges(4,1) * t307 + Ifges(4,2) * t191) * t191 + (mrSges(4,2) * t307 + Ifges(4,1) * t142 + 0.2e1 * Ifges(4,4) * t191) * t142; -t39 * t279 + t263 * t305 + Ifges(3,6) * t165 + Ifges(3,5) * t160 + t8 * t25 + t194 * t9 + (-t30 * t33 + t31 * t32) * mrSges(9,3) + (-t38 * mrSges(5,3) + t2 * t287 + t285 * t3) * t37 + t183 * t21 + ((t78 * t97 - t79 * t98) * mrSges(8,3) + (t164 * t142 - t159 * t191) * mrSges(4,3)) * pkin(1) + t181 * t12 + t64 * (Ifges(7,5) * t100 + Ifges(7,6) * t101) + t63 * (Ifges(8,5) * t79 + Ifges(8,6) * t78) + t230; (t38 ^ 2 + t39 ^ 2) * m(5) + m(9) * t31 ^ 2 + t63 ^ 2 * Ifges(8,3) + t64 ^ 2 * Ifges(7,3) + Ifges(3,3) + Ifges(4,3) + (0.2e1 * t1 + t306) * t8 + t9 ^ 2 * t186 + (t276 - 0.2e1 * t280) * t21 + (m(9) * t30 + 0.2e1 * mrSges(9,1) * t21) * t30 + (m(4) * (t159 ^ 2 + t164 ^ 2) + m(8) * (t97 ^ 2 + t98 ^ 2)) * t178 + 0.2e1 * ((mrSges(8,1) * t98 - mrSges(8,2) * t97) * t63 + t311) * pkin(1) + (Ifges(5,3) * t12 + t193 * t9 - 0.2e1 * t265 + 0.2e1 * t266 + t272 + t275) * t12; t15 * t25 + t262 * t305 - t279 * t296 + t194 * t14 + (-mrSges(5,3) * t295 + t285 * t6 + t287 * t5) * t37 + t183 * t29 + t181 * t18 + t230; t8 * t4 + Ifges(4,3) + (t1 + t306) * t15 + t9 * t14 * t186 + (t38 * t47 + t39 * t46) * pkin(4) * m(5) + t311 * pkin(1) + (mrSges(9,1) * t30 + t276 - t280) * t29 + (-t265 + t266 + t275 / 0.2e1 + t272 / 0.2e1 + t9 * t210) * t18 + (t274 / 0.2e1 + t271 / 0.2e1 + t277 + t14 * t210 + t188) * t12; t29 ^ 2 * Ifges(9,3) + Ifges(4,3) + (t46 ^ 2 + t47 ^ 2) * t174 * m(5) + (m(6) * t15 + 0.2e1 * t4) * t15 + t14 ^ 2 * t186 + (t14 * t193 + 0.2e1 * t188 + t271 + t274 + t277) * t18; t24 * t260 + (-mrSges(6,2) * t24 - Ifges(6,6) * t37) * t158 + t261; -t197 * t9 + t263; -t14 * t197 + t262; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
