% Calculate Gravitation load on the joints for
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1DE2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:05:29
% EndTime: 2020-04-14 20:10:18
% DurationCPUTime: 45.82s
% Computational Cost: add. (1330643->276), mult. (2009044->458), div. (88292->33), fcn. (1273741->48), ass. (0->268)
t353 = pkin(2) * m(10);
t352 = pkin(4) * m(11);
t150 = sin(qJ(2));
t155 = cos(qJ(2));
t173 = pkin(2) ^ 2;
t168 = pkin(6) ^ 2;
t175 = pkin(1) ^ 2;
t274 = t168 + t175;
t257 = -pkin(13) ^ 2 + t274;
t145 = sin(pkin(20));
t147 = cos(pkin(20));
t149 = sin(qJ(3));
t154 = cos(qJ(3));
t128 = t145 * t154 + t147 * t149;
t309 = pkin(6) * t128;
t268 = pkin(1) * t309;
t113 = t173 - t257 - 0.2e1 * t268;
t124 = 0.2e1 * t268;
t277 = t124 + t168;
t332 = -pkin(2) - pkin(13);
t108 = (pkin(1) - t332) * (pkin(1) + t332) + t277;
t327 = pkin(13) - pkin(2);
t109 = (pkin(1) - t327) * (pkin(1) + t327) + t277;
t288 = t109 * t108;
t176 = sqrt(-t288);
t174 = 0.1e1 / pkin(2);
t315 = t174 / 0.2e1;
t251 = 0.1e1 / pkin(13) * t315;
t119 = t124 + t274;
t117 = 0.1e1 / t119;
t253 = t117 * t315;
t114 = t124 + t173 + t257;
t123 = -pkin(1) - t309;
t129 = t145 * t149 - t147 * t154;
t284 = t129 * t176;
t267 = pkin(6) * t284;
t95 = -t114 * t123 - t267;
t308 = pkin(6) * t129;
t96 = t114 * t308 - t123 * t176;
t81 = qJ(2) + atan2(t96 * t253, t95 * t253);
t75 = atan2(t176 * t251, t113 * t251) + t81;
t73 = sin(t75);
t74 = cos(t75);
t225 = mrSges(10,1) * t73 + mrSges(10,2) * t74;
t143 = qJ(2) + qJ(3);
t140 = sin(t143);
t141 = cos(t143);
t243 = t141 * mrSges(4,1) - mrSges(4,2) * t140;
t78 = sin(t81);
t272 = t78 * t353;
t351 = -mrSges(3,1) * t150 - mrSges(3,2) * t155 + t225 + t243 - t272;
t345 = -m(6) - m(5);
t146 = cos(pkin(21));
t162 = 0.1e1 / pkin(11);
t170 = pkin(4) ^ 2;
t169 = pkin(5) ^ 2;
t142 = pkin(23) + pkin(22);
t138 = sin(t142);
t139 = cos(t142);
t167 = pkin(7) ^ 2;
t157 = cos(pkin(19));
t313 = sin(pkin(19));
t130 = t150 * t157 - t155 * t313;
t307 = pkin(7) * t130;
t337 = -0.2e1 * pkin(1);
t276 = t307 * t337 + t175;
t122 = t167 + t276;
t120 = 0.1e1 / t122;
t172 = 0.1e1 / pkin(3);
t285 = t120 * t172;
t318 = t149 / 0.2e1;
t273 = pkin(3) ^ 2 - pkin(8) ^ 2;
t116 = t122 + t273;
t125 = pkin(1) - t307;
t331 = -pkin(3) - pkin(8);
t110 = (pkin(7) - t331) * (pkin(7) + t331) + t276;
t330 = -pkin(8) + pkin(3);
t111 = (pkin(7) - t330) * (pkin(7) + t330) + t276;
t177 = sqrt(-t111 * t110);
t131 = t150 * t313 + t155 * t157;
t306 = pkin(7) * t131;
t99 = t116 * t306 + t125 * t177;
t333 = t99 / 0.2e1;
t282 = t131 * t177;
t98 = -pkin(7) * t282 + t116 * t125;
t87 = (t154 * t333 + t98 * t318) * t285;
t317 = -t154 / 0.2e1;
t88 = (t98 * t317 + t99 * t318) * t285;
t60 = t138 * t88 - t139 * t87;
t323 = pkin(5) * t60;
t293 = -0.2e1 * pkin(4) * t323 + t169;
t54 = t170 + t293;
t52 = 0.1e1 / t54;
t297 = t162 * t52;
t144 = sin(pkin(21));
t319 = t144 / 0.2e1;
t329 = -pkin(9) - pkin(11);
t48 = (pkin(4) - t329) * (pkin(4) + t329) + t293;
t328 = pkin(11) - pkin(9);
t49 = (pkin(4) - t328) * (pkin(4) + t328) + t293;
t178 = sqrt(-t49 * t48);
t246 = t138 * t87 + t139 * t88;
t344 = t178 * t246;
t275 = pkin(9) ^ 2 - pkin(11) ^ 2;
t51 = t54 - t275;
t56 = -pkin(4) * t60 + pkin(5);
t32 = -pkin(4) * t344 + t51 * t56;
t34 = pkin(4) * t246 * t51 + t178 * t56;
t22 = (t32 * t319 + t34 * t146 / 0.2e1) * t297;
t23 = (-t32 * t146 / 0.2e1 + t34 * t319) * t297;
t13 = atan2(t22, t23) + t143;
t11 = sin(t13);
t12 = cos(t13);
t148 = sin(qJ(4));
t153 = cos(qJ(4));
t201 = m(6) * pkin(10) + t153 * mrSges(6,1) - t148 * mrSges(6,2);
t263 = m(6) * pkin(12) + mrSges(6,3);
t350 = (mrSges(5,2) - t263) * t12 + (mrSges(5,1) + t201) * t11;
t151 = sin(qJ(1));
t156 = cos(qJ(1));
t347 = g(1) * t156 + g(2) * t151;
t349 = m(8) + m(4) + m(11);
t291 = sin(pkin(23));
t258 = t291 / 0.2e1;
t292 = cos(pkin(23));
t84 = (-t292 * t98 / 0.2e1 + t99 * t258) * t285;
t85 = (t98 * t258 + t292 * t333) * t285;
t67 = qJ(2) + atan2(t85, t84);
t64 = pkin(22) - t67;
t65 = sin(t67);
t66 = cos(t67);
t348 = sin(t64) * t352 - mrSges(8,1) * t65 - mrSges(8,2) * t66;
t232 = mrSges(5,1) * t12 - mrSges(5,2) * t11;
t343 = -t263 * t11 - t201 * t12 - t232;
t79 = cos(t81);
t227 = -mrSges(9,1) * t78 - mrSges(9,2) * t79;
t342 = -g(3) * t227 + t347 * (mrSges(9,1) * t79 - mrSges(9,2) * t78);
t340 = -mrSges(5,3) - mrSges(10,3) - mrSges(11,3) - mrSges(4,3) - mrSges(9,3) - mrSges(7,3) - mrSges(8,3) - mrSges(3,3) + mrSges(2,2);
t135 = pkin(5) * t141;
t311 = pkin(1) * t150;
t247 = t135 - t311;
t132 = pkin(16) + t247;
t164 = 0.1e1 / pkin(9);
t334 = t52 / 0.2e1;
t259 = t164 * t334;
t50 = t54 + t275;
t55 = -pkin(4) + t323;
t31 = -pkin(5) * t344 - t50 * t55;
t322 = pkin(5) * t246;
t33 = -t178 * t55 + t50 * t322;
t19 = -atan2(t33 * t259, t31 * t259) + t64;
t17 = sin(t19);
t18 = cos(t19);
t223 = mrSges(11,1) * t17 - mrSges(11,2) * t18;
t115 = t122 - t273;
t126 = pkin(1) * t130 - pkin(7);
t100 = pkin(1) * t131 * t115 - t126 * t177;
t152 = sin(pkin(18));
t166 = 0.1e1 / pkin(8);
t286 = t120 * t166;
t158 = cos(pkin(18));
t316 = t158 / 0.2e1;
t97 = -pkin(1) * t282 - t115 * t126;
t89 = (t97 * t316 - t152 * t100 / 0.2e1) * t286;
t90 = (t100 * t316 + t97 * t152 / 0.2e1) * t286;
t72 = atan2(t90, t89);
t68 = sin(t72);
t69 = cos(t72);
t230 = mrSges(7,1) * t69 - mrSges(7,2) * t68;
t338 = m(6) * (pkin(10) * t12 + pkin(12) * t11 + t132) + mrSges(2,1) - t223 + m(5) * t132 + t232 + t227 - m(7) * pkin(15) + t230 + t11 * mrSges(6,3) + t349 * (pkin(16) - t311) + (m(3) + m(10) + m(9)) * pkin(16) + t348 + t351;
t336 = 0.1e1 / t84 ^ 2;
t42 = 0.1e1 / t178;
t335 = -t42 / 0.2e1;
t53 = 0.1e1 / t54 ^ 2;
t324 = pkin(5) * t53;
t321 = pkin(6) * t96;
t320 = t117 / 0.2e1;
t314 = -0.2e1 * t131 ^ 2;
t312 = pkin(1) * t129;
t310 = pkin(5) * t140;
t302 = t155 * pkin(1);
t21 = 0.1e1 / t23 ^ 2;
t301 = t21 * t22;
t118 = 0.1e1 / t119 ^ 2;
t290 = 0.2e1 / t176 * (t108 + t109) * pkin(1) * t308;
t255 = -t290 / 0.2e1;
t94 = 0.1e1 / t95 ^ 2;
t39 = 0.2e1 * (((t123 * t255 + (t114 * t128 - t284) * pkin(6)) * t320 + (-t117 * t129 * t168 + t118 * t321) * t312) / t95 - ((-t128 * t176 + (t255 - t114) * t129) * t320 + (t117 * t123 + t118 * t95) * t312) * t94 * t321) * pkin(2) * t119 * t174 / (t94 * t96 ^ 2 + 0.1e1);
t299 = 0.1e1 / (t21 * t22 ^ 2 + 0.1e1) * t162;
t298 = t146 * t52;
t269 = pkin(1) * t306;
t289 = 0.2e1 / t177 * (t110 + t111) * t269;
t183 = (t125 * t289 / 0.2e1 + t167 * pkin(1) * t314 + (-t116 * t130 - t282) * pkin(7)) * t285;
t254 = -t289 / 0.2e1;
t283 = t130 * t177;
t185 = pkin(7) * (t283 + (t125 * t337 - t116 + t254) * t131) * t285;
t184 = -t185 / 0.2e1;
t245 = 0.1e1 / t122 ^ 2 * t269;
t234 = t172 * t245;
t216 = t154 * t234;
t217 = t149 * t234;
t45 = t149 * t184 + t183 * t317 - t99 * t216 - t98 * t217;
t181 = t183 / 0.2e1;
t46 = t149 * t181 + t154 * t184 - t98 * t216 + t99 * t217;
t37 = t138 * t46 + t139 * t45;
t296 = t178 * t37;
t287 = t120 * t152;
t281 = t148 * t151;
t280 = t148 * t156;
t279 = t151 * t153;
t278 = t153 * t156;
t271 = t79 * t353;
t270 = pkin(4) * t324;
t266 = t55 * t335;
t265 = t42 * t56 / 0.2e1;
t264 = t246 * t335;
t262 = t52 * t319;
t261 = -t298 / 0.2e1;
t260 = t298 / 0.2e1;
t256 = -0.2e1 * pkin(5) * t56 - t51;
t252 = t120 * t316;
t250 = -0.2e1 * t170 * t322;
t249 = t144 * t270;
t248 = t146 * t270;
t30 = 0.1e1 / t31 ^ 2;
t244 = pkin(9) * t164 / (t30 * t33 ^ 2 + 0.1e1) * t54;
t240 = t37 * t249;
t239 = t246 * t249;
t238 = t37 * t248;
t237 = t246 * t248;
t236 = t97 * t245;
t235 = 0.1e1 / t31 * t244;
t233 = t100 * t245;
t226 = mrSges(10,1) * t74 - mrSges(10,2) * t73;
t224 = -mrSges(11,1) * t18 - mrSges(11,2) * t17;
t220 = mrSges(4,1) * t140 + mrSges(4,2) * t141;
t219 = 0.2e1 * pkin(4) * pkin(5) * (t48 + t49);
t215 = pkin(4) * (t31 * t53 + t52 * t55);
t213 = pkin(5) * t30 * t33 * t244;
t25 = t37 * t219;
t38 = -t138 * t45 + t139 * t46;
t210 = -t38 * t178 + t25 * t264;
t35 = t246 * t219;
t209 = t178 * t60 + t35 * t264;
t203 = t292 * t234;
t202 = t291 * t234;
t194 = pkin(4) * (-t169 * t246 * t52 + t33 * t324);
t27 = -0.1e1 + ((t291 * t181 + t292 * t184 + t99 * t202 - t98 * t203) * t85 * t336 - (t292 * t181 + t185 * t258 + t98 * t202 + t99 * t203) / t84) / (t85 ^ 2 * t336 + 0.1e1);
t112 = 0.1e1 / t113 ^ 2;
t86 = 0.1e1 / t89 ^ 2;
t80 = t126 * t254 + t175 * pkin(7) * t314 + (-t115 * t130 - t282) * pkin(1);
t76 = (t283 + (0.2e1 * t126 * pkin(7) - t115 + t254) * t131) * pkin(1);
t36 = (0.1e1 / t113 * t290 / 0.2e1 + t112 * t267 * t337) / (-t112 * t288 + 0.1e1) + t39;
t28 = ((t80 * t252 + t158 * t233 + t76 * t287 / 0.2e1 + t152 * t236) / t89 - (t76 * t252 + t158 * t236 - t80 * t287 / 0.2e1 - t152 * t233) * t90 * t86) / (t86 * t90 ^ 2 + 0.1e1) * t166;
t20 = 0.1e1 / t23;
t16 = t35 * t265 + t246 * t250 + (-t51 * t60 - t344) * pkin(4);
t15 = (t246 * t256 + t209) * pkin(4);
t10 = t12 * t278 + t281;
t9 = -t12 * t280 + t279;
t8 = -t12 * t279 + t280;
t7 = t12 * t281 + t278;
t6 = t25 * t265 + t37 * t250 + (t38 * t51 - t296) * pkin(4);
t5 = (t256 * t37 + t210) * pkin(4);
t4 = -0.2e1 * ((t35 * t266 + (-t50 * t60 - t344) * pkin(5)) * t334 + t246 * t194) * t235 + 0.2e1 * ((-t246 * t50 + t209) * t334 + t246 * t215) * t213;
t3 = -0.2e1 * ((t25 * t266 + (t38 * t50 - t296) * pkin(5)) * t334 + t37 * t194) * t235 + 0.2e1 * ((-t37 * t50 + t210) * t334 + t37 * t215) * t213 + t27;
t2 = 0.1e1 + ((t15 * t262 + t16 * t260 + t34 * t237 + t32 * t239) * t20 - (t15 * t261 + t16 * t262 - t32 * t237 + t34 * t239) * t301) * t299;
t1 = 0.1e1 + ((t34 * t238 + t32 * t240 + t6 * t260 + t5 * t262) * t20 - (-t32 * t238 + t34 * t240 + t5 * t261 + t6 * t262) * t301) * t299;
t14 = [(-t10 * mrSges(6,1) - t9 * mrSges(6,2) + t340 * t151 - t338 * t156) * g(2) + (-t8 * mrSges(6,1) - t7 * mrSges(6,2) + t338 * t151 + t340 * t156) * g(1), (t343 * t1 - t223 * t3 - t230 * t28 + t345 * t247 + t348 * t27 + t349 * t311 - t351) * g(3) + t342 + t347 * (-t224 * t3 + t220 - t226 + t271 - (-mrSges(7,1) * t68 - mrSges(7,2) * t69) * t28 + mrSges(3,1) * t155 - mrSges(3,2) * t150 + (-cos(t64) * t352 - mrSges(8,1) * t66 + mrSges(8,2) * t65) * t27 + t349 * t302 + t345 * (-t302 - t310) + t350 * t1), t342 * t39 + (t345 * t135 + t343 * t2 - t223 * t4 - t225 * t36 + t39 * t272 - t243) * g(3) + t347 * (t350 * t2 - t224 * t4 - t226 * t36 + t39 * t271 - t345 * t310 + t220), -g(1) * (mrSges(6,1) * t9 - mrSges(6,2) * t10) - g(2) * (-mrSges(6,1) * t7 + mrSges(6,2) * t8) - g(3) * (-mrSges(6,1) * t148 - mrSges(6,2) * t153) * t11];
taug = t14(:);
