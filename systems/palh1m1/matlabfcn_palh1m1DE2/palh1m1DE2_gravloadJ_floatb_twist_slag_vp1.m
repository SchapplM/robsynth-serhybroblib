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
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1DE2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:05:29
% EndTime: 2020-04-14 20:10:15
% DurationCPUTime: 45.70s
% Computational Cost: add. (1330644->316), mult. (2009025->550), div. (88292->33), fcn. (1273741->48), ass. (0->272)
t156 = sin(qJ(1));
t161 = cos(qJ(1));
t334 = g(1) * t161 + g(2) * t156;
t338 = g(2) * t161;
t174 = pkin(5) ^ 2;
t147 = pkin(23) + pkin(22);
t142 = sin(t147);
t143 = cos(t147);
t159 = cos(qJ(3));
t172 = pkin(7) ^ 2;
t180 = pkin(1) ^ 2;
t155 = sin(qJ(2));
t160 = cos(qJ(2));
t162 = cos(pkin(19));
t305 = sin(pkin(19));
t133 = t155 * t162 - t160 * t305;
t299 = pkin(7) * t133;
t332 = -2 * pkin(1);
t267 = t299 * t332 + t180;
t122 = t172 + t267;
t120 = 0.1e1 / t122;
t177 = 0.1e1 / pkin(3);
t276 = t120 * t177;
t154 = sin(qJ(3));
t311 = t154 / 0.2e1;
t264 = pkin(3) ^ 2 - pkin(8) ^ 2;
t116 = t122 + t264;
t125 = pkin(1) - t299;
t325 = pkin(7) + pkin(8);
t326 = pkin(7) - pkin(8);
t110 = (pkin(3) + t325) * (-pkin(3) + t326) + t267;
t111 = (-pkin(3) + t325) * (pkin(3) + t326) + t267;
t182 = sqrt(-t111 * t110);
t134 = t155 * t305 + t160 * t162;
t298 = pkin(7) * t134;
t99 = t116 * t298 + t125 * t182;
t328 = t99 / 0.2e1;
t273 = t134 * t182;
t98 = -pkin(7) * t273 + t116 * t125;
t87 = (t159 * t328 + t311 * t98) * t276;
t310 = -t159 / 0.2e1;
t88 = (t310 * t98 + t311 * t99) * t276;
t60 = t142 * t88 - t143 * t87;
t317 = pkin(5) * t60;
t284 = -0.2e1 * pkin(4) * t317 + t174;
t324 = -pkin(9) - pkin(11);
t48 = (pkin(4) - t324) * (pkin(4) + t324) + t284;
t323 = pkin(11) - pkin(9);
t49 = (pkin(4) - t323) * (pkin(4) + t323) + t284;
t183 = sqrt(-t49 * t48);
t240 = t142 * t87 + t143 * t88;
t337 = t183 * t240;
t148 = qJ(2) + qJ(3);
t144 = sin(t148);
t145 = cos(t148);
t336 = t145 * rSges(4,1) - rSges(4,2) * t144;
t178 = pkin(2) ^ 2;
t173 = pkin(6) ^ 2;
t265 = t173 + t180;
t251 = -pkin(13) ^ 2 + t265;
t150 = sin(pkin(20));
t152 = cos(pkin(20));
t131 = t150 * t159 + t152 * t154;
t301 = pkin(6) * t131;
t261 = pkin(1) * t301;
t113 = t178 - t251 - 0.2e1 * t261;
t124 = 0.2e1 * t261;
t268 = t124 + t173;
t327 = -pkin(2) - pkin(13);
t108 = (pkin(1) - t327) * (pkin(1) + t327) + t268;
t322 = pkin(13) - pkin(2);
t109 = (pkin(1) - t322) * (pkin(1) + t322) + t268;
t279 = t109 * t108;
t181 = sqrt(-t279);
t179 = 0.1e1 / pkin(2);
t308 = t179 / 0.2e1;
t245 = 0.1e1 / pkin(13) * t308;
t119 = t124 + t265;
t117 = 0.1e1 / t119;
t247 = t117 * t308;
t114 = t124 + t178 + t251;
t123 = -pkin(1) - t301;
t132 = t150 * t154 - t152 * t159;
t275 = t132 * t181;
t260 = pkin(6) * t275;
t95 = -t114 * t123 - t260;
t300 = pkin(6) * t132;
t96 = t114 * t300 - t123 * t181;
t81 = qJ(2) + atan2(t96 * t247, t95 * t247);
t75 = atan2(t181 * t245, t113 * t245) + t81;
t73 = sin(t75);
t74 = cos(t75);
t218 = rSges(10,1) * t73 + rSges(10,2) * t74;
t78 = sin(t81);
t321 = pkin(2) * t78;
t335 = t218 - t321;
t151 = cos(pkin(21));
t167 = 0.1e1 / pkin(11);
t175 = pkin(4) ^ 2;
t54 = t175 + t284;
t52 = 0.1e1 / t54;
t288 = t167 * t52;
t149 = sin(pkin(21));
t312 = t149 / 0.2e1;
t266 = pkin(9) ^ 2 - pkin(11) ^ 2;
t51 = t54 - t266;
t56 = -pkin(4) * t60 + pkin(5);
t32 = -pkin(4) * t337 + t51 * t56;
t34 = pkin(4) * t240 * t51 + t183 * t56;
t22 = (t32 * t312 + t34 * t151 / 0.2e1) * t288;
t23 = (-t32 * t151 / 0.2e1 + t34 * t312) * t288;
t13 = atan2(t22, t23) + t148;
t11 = sin(t13);
t12 = cos(t13);
t306 = pkin(12) + rSges(6,3);
t333 = pkin(10) * t12 + t11 * t306;
t282 = sin(pkin(23));
t252 = t282 / 0.2e1;
t283 = cos(pkin(23));
t84 = (-t283 * t98 / 0.2e1 + t99 * t252) * t276;
t331 = 0.1e1 / t84 ^ 2;
t42 = 0.1e1 / t183;
t330 = -t42 / 0.2e1;
t329 = t52 / 0.2e1;
t79 = cos(t81);
t320 = pkin(2) * t79;
t85 = (t252 * t98 + t283 * t328) * t276;
t67 = qJ(2) + atan2(t85, t84);
t64 = pkin(22) - t67;
t319 = pkin(4) * sin(t64);
t53 = 0.1e1 / t54 ^ 2;
t318 = pkin(5) * t53;
t316 = pkin(5) * t240;
t315 = pkin(6) * t96;
t313 = t117 / 0.2e1;
t163 = cos(pkin(18));
t309 = t163 / 0.2e1;
t307 = -0.2e1 * t134 ^ 2;
t304 = pkin(1) * t132;
t303 = pkin(1) * t155;
t302 = pkin(5) * t144;
t140 = pkin(16) - t303;
t296 = t140 * t338;
t295 = pkin(16) * t338;
t293 = t160 * pkin(1);
t21 = 0.1e1 / t23 ^ 2;
t292 = t21 * t22;
t118 = 0.1e1 / t119 ^ 2;
t281 = 0.2e1 / t181 * (t108 + t109) * pkin(1) * t300;
t249 = -t281 / 0.2e1;
t94 = 0.1e1 / t95 ^ 2;
t39 = 0.2e1 * (((t123 * t249 + (t114 * t131 - t275) * pkin(6)) * t313 + (-t117 * t132 * t173 + t118 * t315) * t304) / t95 - ((-t131 * t181 + (t249 - t114) * t132) * t313 + (t117 * t123 + t118 * t95) * t304) * t94 * t315) * pkin(2) * t119 * t179 / (t94 * t96 ^ 2 + 0.1e1);
t290 = 0.1e1 / (t21 * t22 ^ 2 + 0.1e1) * t167;
t289 = t151 * t52;
t262 = pkin(1) * t298;
t280 = 0.2e1 / t182 * (t110 + t111) * t262;
t188 = (t125 * t280 / 0.2e1 + t172 * pkin(1) * t307 + (-t116 * t133 - t273) * pkin(7)) * t276;
t248 = -t280 / 0.2e1;
t274 = t133 * t182;
t190 = pkin(7) * (t274 + (t125 * t332 - t116 + t248) * t134) * t276;
t189 = -t190 / 0.2e1;
t239 = 0.1e1 / t122 ^ 2 * t262;
t228 = t177 * t239;
t209 = t159 * t228;
t210 = t154 * t228;
t45 = t154 * t189 + t188 * t310 - t209 * t99 - t210 * t98;
t187 = t188 / 0.2e1;
t46 = t154 * t187 + t159 * t189 - t209 * t98 + t210 * t99;
t37 = t142 * t46 + t143 * t45;
t287 = t183 * t37;
t157 = sin(pkin(18));
t278 = t120 * t157;
t171 = 0.1e1 / pkin(8);
t277 = t120 * t171;
t153 = sin(qJ(4));
t272 = t153 * t156;
t271 = t153 * t161;
t158 = cos(qJ(4));
t270 = t156 * t158;
t269 = t158 * t161;
t263 = pkin(4) * t318;
t55 = -pkin(4) + t317;
t259 = t55 * t330;
t258 = t42 * t56 / 0.2e1;
t257 = t240 * t330;
t256 = t52 * t312;
t255 = -t289 / 0.2e1;
t254 = t289 / 0.2e1;
t169 = 0.1e1 / pkin(9);
t253 = t169 * t329;
t250 = -0.2e1 * pkin(5) * t56 - t51;
t246 = t120 * t309;
t244 = -0.2e1 * t175 * t316;
t243 = t149 * t263;
t242 = t151 * t263;
t139 = pkin(5) * t145;
t241 = t139 - t303;
t50 = t54 + t266;
t31 = -pkin(5) * t337 - t50 * t55;
t30 = 0.1e1 / t31 ^ 2;
t33 = -t183 * t55 + t316 * t50;
t238 = pkin(9) * t169 / (t30 * t33 ^ 2 + 0.1e1) * t54;
t234 = t37 * t243;
t233 = t240 * t243;
t232 = t37 * t242;
t231 = t240 * t242;
t115 = t122 - t264;
t126 = pkin(1) * t133 - pkin(7);
t97 = -pkin(1) * t273 - t115 * t126;
t230 = t97 * t239;
t229 = 0.1e1 / t31 * t238;
t100 = pkin(1) * t134 * t115 - t126 * t182;
t227 = t100 * t239;
t226 = rSges(5,1) * t12 - rSges(5,2) * t11;
t225 = -rSges(5,1) * t11 - rSges(5,2) * t12;
t89 = (t97 * t309 - t157 * t100 / 0.2e1) * t277;
t90 = (t100 * t309 + t97 * t157 / 0.2e1) * t277;
t72 = atan2(t90, t89);
t68 = sin(t72);
t69 = cos(t72);
t224 = rSges(7,1) * t69 - rSges(7,2) * t68;
t65 = sin(t67);
t66 = cos(t67);
t222 = rSges(8,1) * t65 + rSges(8,2) * t66;
t220 = -rSges(9,1) * t78 - rSges(9,2) * t79;
t219 = rSges(10,1) * t74 - rSges(10,2) * t73;
t19 = -atan2(t33 * t253, t31 * t253) + t64;
t17 = sin(t19);
t18 = cos(t19);
t217 = -rSges(11,1) * t18 - rSges(11,2) * t17;
t216 = rSges(11,1) * t17 - rSges(11,2) * t18;
t214 = -rSges(3,1) * t155 - rSges(3,2) * t160;
t212 = -rSges(4,1) * t144 - rSges(4,2) * t145;
t211 = 0.2e1 * pkin(4) * pkin(5) * (t48 + t49);
t208 = -pkin(15) + t224;
t207 = pkin(16) + t335;
t206 = rSges(6,1) * t158 - rSges(6,2) * t153 + pkin(10);
t205 = t216 - t140 - t319;
t204 = pkin(4) * (t31 * t53 + t52 * t55);
t203 = pkin(5) * t30 * t33 * t238;
t25 = t37 * t211;
t38 = -t142 * t45 + t143 * t46;
t202 = -t183 * t38 + t25 * t257;
t35 = t240 * t211;
t201 = t183 * t60 + t257 * t35;
t199 = t283 * t228;
t198 = t282 * t228;
t196 = pkin(4) * (-t174 * t240 * t52 + t318 * t33);
t27 = -0.1e1 + (t331 * (t187 * t282 + t189 * t283 + t198 * t99 - t199 * t98) * t85 - (t187 * t283 + t190 * t252 + t198 * t98 + t199 * t99) / t84) / (t331 * t85 ^ 2 + 0.1e1);
t191 = g(3) * t241 + t334 * (-t293 - t302);
t186 = m(9) * (g(3) * t220 + t334 * (-rSges(9,1) * t79 + rSges(9,2) * t78));
t185 = (g(3) * t206 + t306 * t334) * t12 + (g(3) * t306 - t206 * t334) * t11;
t135 = pkin(16) + t241;
t128 = t161 * t135;
t112 = 0.1e1 / t113 ^ 2;
t86 = 0.1e1 / t89 ^ 2;
t80 = t126 * t248 + t180 * pkin(7) * t307 + (-t115 * t133 - t273) * pkin(1);
t76 = (t274 + (0.2e1 * pkin(7) * t126 - t115 + t248) * t134) * pkin(1);
t36 = (0.1e1 / t113 * t281 / 0.2e1 + t112 * t260 * t332) / (-t112 * t279 + 0.1e1) + t39;
t20 = 0.1e1 / t23;
t16 = t35 * t258 + t240 * t244 + (-t51 * t60 - t337) * pkin(4);
t15 = (t240 * t250 + t201) * pkin(4);
t10 = t12 * t269 + t272;
t9 = -t12 * t271 + t270;
t8 = -t12 * t270 + t271;
t7 = t12 * t272 + t269;
t6 = t25 * t258 + t37 * t244 + (t38 * t51 - t287) * pkin(4);
t5 = (t250 * t37 + t202) * pkin(4);
t3 = -0.2e1 * ((t25 * t259 + (t38 * t50 - t287) * pkin(5)) * t329 + t37 * t196) * t229 + 0.2e1 * ((-t37 * t50 + t202) * t329 + t37 * t204) * t203 + t27;
t2 = 0.1e1 + ((t15 * t256 + t16 * t254 + t231 * t34 + t233 * t32) * t20 - (t15 * t255 + t16 * t256 - t231 * t32 + t233 * t34) * t292) * t290;
t1 = 0.1e1 + ((t232 * t34 + t234 * t32 + t254 * t6 + t256 * t5) * t20 - (-t232 * t32 + t234 * t34 + t255 * t5 + t256 * t6) * t292) * t290;
t4 = [-m(2) * (g(1) * (-rSges(2,1) * t156 - rSges(2,2) * t161) + g(2) * (rSges(2,1) * t161 - rSges(2,2) * t156)) - m(3) * (t295 + (g(1) * rSges(3,3) + g(2) * t214) * t161 + (g(1) * (-pkin(16) - t214) + g(2) * rSges(3,3)) * t156) - m(4) * (t296 + (g(1) * rSges(4,3) + g(2) * t336) * t161 + (g(1) * (-t140 - t336) + g(2) * rSges(4,3)) * t156) - m(5) * (g(2) * t128 + (g(1) * rSges(5,3) + g(2) * t226) * t161 + (g(1) * (-t135 - t226) + g(2) * rSges(5,3)) * t156) - m(6) * ((rSges(6,1) * t10 + rSges(6,2) * t9 + t161 * t333 + t128) * g(2) + (rSges(6,1) * t8 + rSges(6,2) * t7 + (-t135 - t333) * t156) * g(1)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t208) * t161 + (g(2) * rSges(7,3) - g(1) * t208) * t156) - m(8) * (t296 + (g(1) * rSges(8,3) - g(2) * t222) * t161 + (g(1) * (-t140 + t222) + g(2) * rSges(8,3)) * t156) - m(9) * (t295 + (g(1) * rSges(9,3) + g(2) * t220) * t161 + (g(1) * (-pkin(16) - t220) + g(2) * rSges(9,3)) * t156) - m(10) * ((g(1) * rSges(10,3) + g(2) * t207) * t161 + (g(2) * rSges(10,3) - g(1) * t207) * t156) - m(11) * ((g(1) * rSges(11,3) - g(2) * t205) * t161 + (g(2) * rSges(11,3) + g(1) * t205) * t156), -m(3) * (g(3) * t214 + t334 * (-rSges(3,1) * t160 + rSges(3,2) * t155)) - m(4) * (g(3) * (t336 - t303) + t334 * (t212 - t293)) - m(5) * ((g(3) * t226 + t225 * t334) * t1 + t191) - m(6) * (t1 * t185 + t191) - m(7) * (g(3) * t224 + t334 * (-rSges(7,1) * t68 - rSges(7,2) * t69)) * ((t80 * t246 + t163 * t227 + t76 * t278 / 0.2e1 + t157 * t230) / t89 - (t76 * t246 + t163 * t230 - t80 * t278 / 0.2e1 - t157 * t227) * t90 * t86) / (t86 * t90 ^ 2 + 0.1e1) * t171 - m(8) * (g(3) * (t222 * t27 - t303) + t334 * (-t293 + (rSges(8,1) * t66 - rSges(8,2) * t65) * t27)) - t186 - m(10) * (g(3) * t335 + t334 * (t219 - t320)) - m(11) * (g(3) * (t216 * t3 - t27 * t319 - t303) + t334 * (t217 * t3 + pkin(4) * t27 * cos(t64) - t293)), -m(4) * (g(3) * t336 + t334 * t212) - m(5) * (g(3) * (t2 * t226 + t139) + t334 * (t2 * t225 - t302)) - m(6) * (g(3) * t139 + t185 * t2 - t334 * t302) - t39 * t186 - m(10) * (g(3) * (t218 * t36 - t321 * t39) + t334 * (t219 * t36 - t320 * t39)) - 0.2e1 * m(11) * (g(3) * t216 + t217 * t334) * (-((t35 * t259 + (-t50 * t60 - t337) * pkin(5)) * t329 + t240 * t196) * t229 + ((-t240 * t50 + t201) * t329 + t240 * t204) * t203), -m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t10) + g(2) * (-rSges(6,1) * t7 + rSges(6,2) * t8) + g(3) * (-rSges(6,1) * t153 - rSges(6,2) * t158) * t11)];
taug = t4(:);
