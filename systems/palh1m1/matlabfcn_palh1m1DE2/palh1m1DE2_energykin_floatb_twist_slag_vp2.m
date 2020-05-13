% Calculate kinetic energy for
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1DE2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:02:08
% EndTime: 2020-04-14 20:03:01
% DurationCPUTime: 46.29s
% Computational Cost: add. (1189089->387), mult. (1834017->669), div. (77340->33), fcn. (1147456->46), ass. (0->277)
t321 = -2 * pkin(1);
t231 = pkin(4) ^ 2;
t230 = pkin(5) ^ 2;
t228 = pkin(7) ^ 2;
t236 = pkin(1) ^ 2;
t209 = sin(qJ(2));
t211 = sin(pkin(19));
t215 = cos(qJ(2));
t217 = cos(pkin(19));
t185 = t209 * t217 - t211 * t215;
t295 = pkin(7) * t185;
t267 = t295 * t321 + t236;
t169 = t228 + t267;
t264 = pkin(3) ^ 2 - pkin(8) ^ 2;
t161 = t169 + t264;
t176 = pkin(1) - t295;
t187 = t209 * t211 + t215 * t217;
t317 = pkin(7) + pkin(8);
t318 = pkin(7) - pkin(8);
t153 = (pkin(3) + t317) * (-pkin(3) + t318) + t267;
t154 = (-pkin(3) + t317) * (pkin(3) + t318) + t267;
t238 = sqrt(-t154 * t153);
t126 = pkin(7) * t161 * t187 + t176 * t238;
t214 = cos(qJ(3));
t269 = t214 * t126;
t276 = t187 * t238;
t125 = -pkin(7) * t276 + t161 * t176;
t208 = sin(qJ(3));
t272 = t208 * t125;
t242 = t272 / 0.2e1 + t269 / 0.2e1;
t167 = 0.1e1 / t169;
t233 = 0.1e1 / pkin(3);
t278 = t167 * t233;
t112 = t242 * t278;
t270 = t214 * t125;
t271 = t208 * t126;
t241 = -t270 / 0.2e1 + t271 / 0.2e1;
t113 = t241 * t278;
t198 = pkin(23) + pkin(22);
t195 = sin(t198);
t196 = cos(t198);
t90 = -t112 * t196 + t113 * t195;
t309 = pkin(5) * t90;
t285 = -0.2e1 * pkin(4) * t309 + t230;
t86 = t231 + t285;
t84 = 0.1e1 / t86;
t320 = t84 / 0.2e1;
t319 = (-pkin(2) - pkin(13));
t316 = -pkin(9) - pkin(11);
t315 = pkin(11) - pkin(9);
t314 = (pkin(13) - pkin(2));
t225 = 0.1e1 / pkin(9);
t257 = t225 * t320;
t80 = (pkin(4) - t316) * (pkin(4) + t316) + t285;
t81 = (pkin(4) - t315) * (pkin(4) + t315) + t285;
t239 = sqrt(-t81 * t80);
t91 = t112 * t195 + t113 * t196;
t290 = t239 * t91;
t266 = pkin(9) ^ 2 - pkin(11) ^ 2;
t82 = t86 + t266;
t87 = -pkin(4) + t309;
t44 = -pkin(5) * t290 - t82 * t87;
t308 = pkin(5) * t91;
t46 = -t239 * t87 + t308 * t82;
t34 = atan2(t46 * t257, t44 * t257);
t313 = sin(t34);
t181 = t187 * qJD(2);
t180 = t185 * qJD(2);
t261 = pkin(1) * pkin(7) * t181;
t282 = 0.2e1 * (t153 + t154) * t261 / t238;
t253 = -t282 / 0.2e1;
t240 = t180 * t238 + t187 * t253;
t102 = ((t176 * t321 - t161) * t181 + t240) * pkin(7);
t259 = -0.2e1 * t181 * t187;
t277 = t181 * t238;
t103 = t176 * t282 / 0.2e1 + t228 * pkin(1) * t259 + (-t161 * t180 - t277) * pkin(7);
t203 = cos(pkin(23));
t274 = t203 * t125;
t199 = sin(pkin(23));
t275 = t199 * t126;
t109 = (-t274 / 0.2e1 + t275 / 0.2e1) * t278;
t108 = 0.1e1 / t109 ^ 2;
t273 = t203 * t126;
t284 = t125 * t199;
t110 = (t273 / 0.2e1 + t284 / 0.2e1) * t278;
t210 = sin(qJ(1));
t216 = cos(qJ(1));
t186 = -t210 * V_base(4) + t216 * V_base(5);
t182 = qJD(2) - t186;
t249 = 0.1e1 / t169 ^ 2 * t261;
t304 = t199 / 0.2e1;
t306 = t103 / 0.2e1;
t307 = -t102 / 0.2e1;
t39 = t182 + (-((t103 * t304 + t203 * t307) * t167 + (-t274 + t275) * t249) * t110 * t108 + ((t102 * t304 + t203 * t306) * t167 + (t273 + t284) * t249) / t109) * t233 / (t108 * t110 ^ 2 + 0.1e1);
t312 = pkin(4) * t39;
t302 = -t214 / 0.2e1;
t57 = ((-t270 + t271) * t249 + (qJD(3) * t242 + t102 * t302 + t208 * t306) * t167) * t233;
t58 = ((-t269 - t272) * t249 + (qJD(3) * t241 + t103 * t302 + t208 * t307) * t167) * t233;
t50 = t195 * t57 + t196 * t58;
t311 = pkin(4) * t50;
t310 = pkin(5) * t46;
t202 = sin(pkin(20));
t206 = cos(pkin(20));
t244 = t202 * t214 + t206 * t208;
t297 = pkin(6) * t244;
t260 = pkin(1) * t297;
t173 = 0.2e1 * t260;
t229 = pkin(6) ^ 2;
t265 = t229 + t236;
t166 = t173 + t265;
t164 = 0.1e1 / t166;
t305 = t164 / 0.2e1;
t201 = sin(pkin(21));
t303 = t201 / 0.2e1;
t218 = cos(pkin(18));
t301 = t218 / 0.2e1;
t235 = 0.1e1 / pkin(2);
t300 = t235 / 0.2e1;
t184 = t202 * t208 - t206 * t214;
t175 = t184 * qJD(3);
t299 = pkin(1) * t175;
t234 = pkin(2) ^ 2;
t256 = -pkin(13) ^ 2 + t265;
t159 = t173 + t234 + t256;
t172 = -pkin(1) - t297;
t268 = t173 + t229;
t151 = ((pkin(1) - t319) * (pkin(1) + t319)) + t268;
t152 = ((pkin(1) - t314) * (pkin(1) + t314)) + t268;
t281 = t152 * t151;
t237 = sqrt(-t281);
t296 = pkin(6) * t184;
t123 = t159 * t296 - t172 * t237;
t298 = pkin(6) * t123;
t263 = pkin(5) * t311;
t294 = 0.2e1 * (t80 + t81) * t263 / t239;
t192 = pkin(14) * V_base(5) + V_base(1);
t193 = -pkin(14) * V_base(4) + V_base(2);
t163 = t216 * t192 + t210 * t193;
t170 = -pkin(16) * t186 + V_base(3);
t146 = -t163 * t215 - t170 * t209;
t135 = pkin(1) * t182 + t146;
t145 = -t163 * t209 + t170 * t215;
t119 = t208 * t135 + t214 * t145;
t179 = qJD(3) + t182;
t117 = pkin(5) * t179 + t119;
t118 = -t135 * t214 + t145 * t208;
t205 = cos(pkin(21));
t83 = t86 - t266;
t88 = -pkin(4) * t90 + pkin(5);
t47 = pkin(4) * t83 * t91 + t239 * t88;
t286 = t47 * t205;
t45 = -pkin(4) * t290 + t83 * t88;
t289 = t45 * t201;
t223 = 0.1e1 / pkin(11);
t292 = t223 * t84;
t31 = (t289 / 0.2e1 + t286 / 0.2e1) * t292;
t287 = t47 * t201;
t288 = t45 * t205;
t32 = (-t288 / 0.2e1 + t287 / 0.2e1) * t292;
t28 = atan2(t31, t32);
t25 = sin(t28);
t26 = cos(t28);
t15 = t25 * t117 + t26 * t118;
t97 = atan2(t110, t109);
t92 = sin(t97);
t93 = cos(t97);
t60 = t92 * t135 + t93 * t145;
t291 = t239 * t50;
t262 = pkin(6) * t299;
t283 = 0.2e1 * (t151 + t152) * t262 / t237;
t212 = sin(pkin(18));
t280 = t167 * t212;
t227 = 0.1e1 / pkin(8);
t279 = t167 * t227;
t258 = -t294 / 0.2e1;
t85 = 0.1e1 / t86 ^ 2;
t255 = t85 * t263;
t254 = -t283 / 0.2e1;
t252 = t164 * t300;
t251 = t167 * t301;
t250 = 0.1e1 / pkin(13) * t300;
t59 = t93 * t135 - t145 * t92;
t122 = -t159 * t172 - t237 * t296;
t107 = atan2(t123 * t252, t122 * t252);
t105 = sin(t107);
t106 = cos(t107);
t72 = -t105 * t145 + t106 * t146;
t162 = -t210 * t192 + t193 * t216;
t246 = t212 * t249;
t245 = t218 * t249;
t14 = t117 * t26 - t118 * t25;
t188 = t210 * V_base(5) + t216 * V_base(4);
t197 = V_base(6) + qJD(1);
t155 = -t188 * t209 + t197 * t215;
t156 = -t188 * t215 - t197 * t209;
t136 = t155 * t208 - t156 * t214;
t137 = t155 * t214 + t156 * t208;
t17 = -t136 * t25 + t137 * t26;
t150 = -pkin(16) * t197 - t162;
t49 = -t195 * t58 + t196 * t57;
t243 = -t49 * t239 + t258 * t91;
t121 = 0.1e1 / t122 ^ 2;
t165 = 0.1e1 / t166 ^ 2;
t174 = t244 * qJD(3);
t53 = t182 + 0.2e1 * (-((-t175 * t159 - t174 * t237 + t184 * t254) * t305 + (t122 * t165 + t164 * t172) * t299) * t121 * t298 + ((t172 * t254 + (t159 * t174 - t175 * t237) * pkin(6)) * t305 + (-t164 * t184 * t229 + t165 * t298) * t299) / t122) * pkin(2) / (t121 * t123 ^ 2 + 0.1e1) * t166 * t235;
t139 = -pkin(1) * t156 + t150;
t120 = -pkin(5) * t137 + t139;
t219 = V_base(3) ^ 2;
t213 = cos(qJ(4));
t207 = sin(qJ(4));
t204 = cos(pkin(22));
t200 = sin(pkin(22));
t177 = pkin(1) * t185 - pkin(7);
t171 = pkin(15) * t186 + V_base(3);
t160 = t169 - t264;
t158 = t234 - t256 - 0.2e1 * t260;
t157 = 0.1e1 / t158 ^ 2;
t149 = t150 ^ 2;
t148 = -pkin(17) * t186 + t163;
t147 = pkin(15) * t197 - pkin(17) * t188 - t162;
t138 = t139 ^ 2;
t132 = atan2(t237 * t250, t158 * t250);
t130 = cos(t132);
t129 = sin(t132);
t127 = pkin(1) * t160 * t187 - t177 * t238;
t124 = -pkin(1) * t276 - t160 * t177;
t115 = (t127 * t301 + t124 * t212 / 0.2e1) * t279;
t114 = (t124 * t301 - t212 * t127 / 0.2e1) * t279;
t111 = 0.1e1 / t114 ^ 2;
t104 = t177 * t253 + t236 * pkin(7) * t259 + (-t160 * t180 - t277) * pkin(1);
t101 = ((0.2e1 * pkin(7) * t177 - t160) * t181 + t240) * pkin(1);
t99 = atan2(t115, t114);
t96 = cos(t99);
t95 = sin(t99);
t78 = t105 * t156 + t106 * t155;
t77 = -t105 * t155 + t106 * t156;
t74 = -pkin(2) * t77 + t150;
t73 = t105 * t146 + t106 * t145;
t68 = t188 * t96 + t197 * t95;
t67 = -t188 * t95 + t197 * t96;
t64 = t155 * t93 + t156 * t92;
t63 = -t155 * t92 + t156 * t93;
t62 = t148 * t96 + t171 * t95;
t61 = -t148 * t95 + t171 * t96;
t56 = -t129 * t77 - t130 * t78;
t55 = t129 * t78 - t130 * t77;
t52 = (0.1e1 / t158 * t283 / 0.2e1 - 0.2e1 * t237 * t157 * t262) / (-t157 * t281 + 0.1e1) + t53;
t51 = (-t200 * t64 - t204 * t63) * pkin(4) + t139;
t48 = pkin(2) * t53 + t72;
t43 = 0.1e1 / t44 ^ 2;
t42 = -t129 * t48 - t130 * t73;
t41 = t129 * t73 - t130 * t48;
t40 = ((t104 * t251 + t127 * t245 + t101 * t280 / 0.2e1 + t124 * t246) / t114 - (t101 * t251 + t124 * t245 - t104 * t280 / 0.2e1 - t127 * t246) * t115 * t111) / (t111 * t115 ^ 2 + 0.1e1) * t227 - t186;
t38 = t200 * t312 + t60;
t37 = t204 * t312 + t59;
t33 = cos(t34);
t30 = 0.1e1 / t32 ^ 2;
t24 = -t200 * t313 - t204 * t33;
t23 = t200 * t33 - t204 * t313;
t20 = t88 * t294 / 0.2e1 - 0.2e1 * t231 * t50 * t308 + (t49 * t83 - t291) * pkin(4);
t19 = ((-0.2e1 * pkin(5) * t88 - t83) * t50 + t243) * pkin(4);
t18 = t136 * t26 + t137 * t25;
t16 = qJD(4) - t17;
t13 = t23 * t63 + t24 * t64;
t12 = -t23 * t64 + t24 * t63;
t11 = t23 * t37 + t24 * t38;
t10 = -t23 * t38 + t24 * t37;
t9 = -pkin(10) * t17 - pkin(12) * t18 + t120;
t8 = t39 + 0.2e1 * (((t87 * t258 + (t49 * t82 - t291) * pkin(5)) * t320 + (-t230 * t84 * t91 + t310 * t85) * t311) / t44 - ((-t50 * t82 + t243) * t320 + (t44 * t85 + t84 * t87) * t311) * t43 * t310) * pkin(9) * t225 / (t43 * t46 ^ 2 + 0.1e1) * t86;
t7 = t179 + (((t19 * t303 + t20 * t205 / 0.2e1) * t84 + (t286 + t289) * t255) / t32 - ((-t19 * t205 / 0.2e1 + t20 * t303) * t84 + (t287 - t288) * t255) * t31 * t30) * t223 / (t30 * t31 ^ 2 + 0.1e1);
t6 = t18 * t213 + t207 * t7;
t5 = -t18 * t207 + t213 * t7;
t4 = pkin(12) * t7 + t15;
t3 = -pkin(10) * t7 - t14;
t2 = t207 * t9 + t213 * t4;
t1 = -t207 * t4 + t213 * t9;
t21 = (V_base(3) * mrSges(2,2) - t162 * mrSges(2,3) + Ifges(2,5) * t197 + Ifges(2,1) * t188 / 0.2e1) * t188 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t6 + Ifges(6,6) * t5 + Ifges(6,3) * t16 / 0.2e1) * t16 + (-t139 * mrSges(4,1) + t118 * mrSges(4,3) + Ifges(4,6) * t179 + Ifges(4,2) * t137 / 0.2e1) * t137 + (-t139 * mrSges(8,1) + t60 * mrSges(8,3) + Ifges(8,4) * t64 + Ifges(8,2) * t63 / 0.2e1) * t63 + (-t120 * mrSges(5,1) + t15 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t7 + Ifges(5,2) * t17 / 0.2e1) * t17 + (t51 * mrSges(11,2) - t10 * mrSges(11,3) + Ifges(11,5) * t8 + Ifges(11,1) * t13 / 0.2e1) * t13 + (t146 * mrSges(3,1) - t145 * mrSges(3,2) + Ifges(3,3) * t182 / 0.2e1) * t182 + (t147 * mrSges(7,2) - t61 * mrSges(7,3) + Ifges(7,1) * t68 / 0.2e1) * t68 + (t150 * mrSges(9,2) - t72 * mrSges(9,3) + Ifges(9,1) * t78 / 0.2e1) * t78 + m(7) * (t147 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(3) * (t145 ^ 2 + t146 ^ 2 + t149) / 0.2e1 + m(9) * (t72 ^ 2 + t73 ^ 2 + t149) / 0.2e1 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t6 / 0.2e1) * t6 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t150 * mrSges(9,1) + t73 * mrSges(9,3) + Ifges(9,4) * t78 + Ifges(9,2) * t77 / 0.2e1) * t77 + m(5) * (t120 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + (t14 * mrSges(5,1) - t15 * mrSges(5,2) + Ifges(5,3) * t7 / 0.2e1) * t7 + (t119 * mrSges(4,1) - t118 * mrSges(4,2) + Ifges(4,3) * t179 / 0.2e1) * t179 + m(4) * (t118 ^ 2 + t119 ^ 2 + t138) / 0.2e1 + m(8) * (t59 ^ 2 + t60 ^ 2 + t138) / 0.2e1 + (t150 * mrSges(3,2) - t146 * mrSges(3,3) + Ifges(3,4) * t156 + Ifges(3,5) * t182 + Ifges(3,1) * t155 / 0.2e1) * t155 + (t10 * mrSges(11,1) - t11 * mrSges(11,2) + Ifges(11,3) * t8 / 0.2e1) * t8 + (t139 * mrSges(8,2) - t59 * mrSges(8,3) + Ifges(8,1) * t64 / 0.2e1) * t64 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t6 + Ifges(6,2) * t5 / 0.2e1) * t5 + (t74 * mrSges(10,2) - t41 * mrSges(10,3) + Ifges(10,1) * t56 / 0.2e1) * t56 + (-t51 * mrSges(11,1) + t11 * mrSges(11,3) + Ifges(11,4) * t13 + Ifges(11,6) * t8 + Ifges(11,2) * t12 / 0.2e1) * t12 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t139 * mrSges(4,2) - t119 * mrSges(4,3) + Ifges(4,4) * t137 + Ifges(4,5) * t179 + Ifges(4,1) * t136 / 0.2e1) * t136 + (t61 * mrSges(7,1) - t62 * mrSges(7,2) + Ifges(7,5) * t68 + Ifges(7,6) * t67 + Ifges(7,3) * t40 / 0.2e1) * t40 + (t120 * mrSges(5,2) - t14 * mrSges(5,3) + Ifges(5,5) * t7 + Ifges(5,1) * t18 / 0.2e1) * t18 + (-V_base(3) * mrSges(2,1) + t163 * mrSges(2,3) + Ifges(2,4) * t188 + Ifges(2,6) * t197 + Ifges(2,2) * t186 / 0.2e1) * t186 + (t59 * mrSges(8,1) - t60 * mrSges(8,2) + Ifges(8,5) * t64 + Ifges(8,6) * t63 + Ifges(8,3) * t39 / 0.2e1) * t39 + (-t147 * mrSges(7,1) + t62 * mrSges(7,3) + Ifges(7,4) * t68 + Ifges(7,2) * t67 / 0.2e1) * t67 + (-t74 * mrSges(10,1) + t42 * mrSges(10,3) + Ifges(10,4) * t56 + Ifges(10,2) * t55 / 0.2e1) * t55 + (t72 * mrSges(9,1) - t73 * mrSges(9,2) + Ifges(9,5) * t78 + Ifges(9,6) * t77 + Ifges(9,3) * t53 / 0.2e1) * t53 + m(10) * (t41 ^ 2 + t42 ^ 2 + t74 ^ 2) / 0.2e1 + (-t150 * mrSges(3,1) + t145 * mrSges(3,3) + Ifges(3,6) * t182 + Ifges(3,2) * t156 / 0.2e1) * t156 + (t41 * mrSges(10,1) - t42 * mrSges(10,2) + Ifges(10,5) * t56 + Ifges(10,6) * t55 + Ifges(10,3) * t52 / 0.2e1) * t52 + m(11) * (t10 ^ 2 + t11 ^ 2 + t51 ^ 2) / 0.2e1 + (t162 * mrSges(2,1) - t163 * mrSges(2,2) + Ifges(2,3) * t197 / 0.2e1) * t197 + m(2) * (t162 ^ 2 + t163 ^ 2 + t219) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t219) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1;
T = t21;
