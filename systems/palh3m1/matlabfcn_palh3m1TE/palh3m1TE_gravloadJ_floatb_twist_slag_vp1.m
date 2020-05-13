% Calculate Gravitation load on the joints for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1TE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:22:20
% EndTime: 2020-04-17 15:27:21
% DurationCPUTime: 48.49s
% Computational Cost: add. (941692->370), mult. (1427729->599), div. (59424->10), fcn. (908093->24), ass. (0->248)
t293 = sin(qJ(2));
t294 = sin(pkin(16));
t296 = cos(qJ(2));
t297 = cos(pkin(16));
t222 = t293 * t294 - t296 * t297;
t218 = pkin(5) * t222;
t212 = (-0.2e1 * t218 + pkin(1)) * pkin(1);
t311 = pkin(5) ^ 2;
t104 = t212 + t311;
t137 = pkin(2) ^ 2;
t265 = pkin(6) ^ 2 - t137;
t209 = t104 - t265;
t213 = -t218 + pkin(1);
t116 = t293 * t297 + t294 * t296;
t304 = pkin(5) + pkin(6);
t305 = pkin(5) - pkin(6);
t136 = sqrt(-((-pkin(2) + t304) * (pkin(2) + t305) + t212) * ((pkin(2) + t304) * (-pkin(2) + t305) + t212));
t272 = t116 * t136;
t256 = pkin(5) * t272;
t313 = 0.1e1 / pkin(2);
t203 = t313 * (t209 * t213 - t256);
t307 = 0.1e1 / t104;
t202 = t307 * t203;
t200 = -t202 / 0.2e1;
t207 = pkin(5) * t209;
t204 = t313 * (t116 * t207 + t136 * t213);
t262 = t307 / 0.2e1;
t201 = t204 * t262;
t292 = sin(qJ(3));
t295 = cos(qJ(3));
t195 = t200 * t295 + t201 * t292;
t199 = t202 / 0.2e1;
t196 = t199 * t292 + t201 * t295;
t259 = pkin(18) + pkin(19);
t243 = sin(t259);
t244 = cos(t259);
t315 = t195 * t244 + t196 * t243;
t187 = pkin(3) * t315;
t179 = (0.2e1 * t187 + pkin(4)) * pkin(4);
t302 = -pkin(8) + pkin(10);
t303 = -pkin(8) - pkin(10);
t135 = sqrt(-((pkin(3) - t302) * (pkin(3) + t302) + t179) * ((pkin(3) - t303) * (pkin(3) + t303) + t179));
t314 = t195 * t243 - t196 * t244;
t336 = t314 * t135;
t339 = pkin(3) * t336;
t338 = pkin(4) * t336;
t139 = pkin(8) ^ 2;
t264 = pkin(10) ^ 2 - t139;
t186 = pkin(4) * t315;
t312 = pkin(4) ^ 2;
t66 = t312 + (0.2e1 * t186 + pkin(3)) * pkin(3);
t176 = t66 + t264;
t181 = t187 + pkin(4);
t309 = 0.1e1 / pkin(10);
t161 = t309 * (t176 * t181 - t339);
t174 = pkin(3) * t176;
t162 = t309 * (t135 * t181 + t174 * t314);
t273 = sin(pkin(17));
t274 = cos(pkin(17));
t337 = t273 * t161 + t274 * t162;
t277 = sin(pkin(19));
t278 = cos(pkin(19));
t83 = t200 * t278 + t201 * t277;
t84 = t199 * t277 + t201 * t278;
t73 = -t293 * t84 + t296 * t83;
t128 = sin(qJ(4));
t131 = cos(qJ(4));
t132 = cos(qJ(1));
t115 = -t292 * t296 - t293 * t295;
t129 = sin(qJ(1));
t105 = t115 * t129;
t114 = t292 * t293 - t295 * t296;
t106 = t114 * t129;
t156 = t274 * t161;
t157 = t273 * t162;
t308 = 0.1e1 / t66;
t269 = t308 / 0.2e1;
t270 = -t308 / 0.2e1;
t43 = t156 * t270 + t157 * t269;
t44 = t337 * t269;
t29 = -t105 * t44 + t106 * t43;
t335 = t132 * t128 - t131 * t29;
t334 = t128 * t29 + t131 * t132;
t333 = -0.2e1 * t314;
t215 = t222 * t136;
t291 = pkin(1) * t116;
t258 = pkin(5) * t291;
t287 = 0.4e1 / t136 * (t304 * t305 - t137 + t212) * t258;
t252 = -t287 / 0.2e1;
t78 = (t215 + (t252 + t265 - t311) * t116) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t218) * t258;
t332 = -t78 / 0.2e1;
t260 = -0.2e1 * pkin(1) * t116 ^ 2;
t79 = -t256 + t213 * t287 / 0.2e1 - t222 * t207 + t311 * t260;
t331 = t79 / 0.2e1;
t126 = sin(pkin(18));
t330 = t126 / 0.2e1;
t127 = cos(pkin(18));
t329 = -t127 / 0.2e1;
t328 = -t314 / 0.2e1;
t327 = t273 / 0.2e1;
t326 = -t274 / 0.2e1;
t298 = -pkin(11) - rSges(6,3);
t197 = t203 * t258;
t198 = t204 * t258;
t253 = t313 * t262;
t230 = t277 * t253;
t290 = t313 * t307;
t238 = t278 * t290;
t99 = 0.1e1 / t104 ^ 2;
t61 = t238 * t331 + t78 * t230 + (t197 * t277 + t198 * t278) * t99;
t62 = t238 * t332 + t79 * t230 + (-t197 * t278 + t198 * t277) * t99;
t58 = t293 * t62 + t296 * t61 + t73;
t306 = pkin(3) * pkin(4);
t261 = 0.1e1 / t66 ^ 2 * t306;
t325 = (-t156 + t157) * t261;
t324 = t337 * t261;
t316 = 0.4e1 / t135 * ((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t179 - t139) * t306;
t170 = t314 * t316;
t323 = -t135 * t315 + t170 * t328;
t231 = t292 * t253;
t242 = t295 * t290;
t191 = t242 * t331 + t78 * t231 + (t197 * t292 + t198 * t295) * t99;
t192 = t242 * t332 + t79 * t231 + (-t197 * t295 + t198 * t292) * t99;
t167 = -t191 * t244 + t192 * t243;
t60 = -t191 * t243 - t192 * t244;
t171 = t60 * t316;
t322 = -t135 * t167 + t171 * t328;
t175 = t66 - t264;
t173 = pkin(4) * t175;
t180 = -t186 - pkin(3);
t321 = 0.2e1 * t180 * t306 - t173;
t320 = -0.2e1 * t181 * t306 - t174;
t310 = 0.1e1 / pkin(8);
t163 = t310 * (-t135 * t180 + t173 * t314);
t159 = t163 * t261;
t53 = -t175 * t180 - t338;
t299 = t53 * t310;
t239 = t261 * t299;
t319 = -t126 * t159 - t127 * t239;
t318 = t126 * t239 - t127 * t159;
t317 = g(1) * t132 + g(2) * t129;
t301 = t310 * t308;
t300 = t309 * t308;
t169 = t171 / 0.2e1;
t182 = pkin(4) * pkin(3) ^ 2 * t333;
t254 = t309 * t269;
t276 = t60 * t135;
t141 = (-pkin(3) * t276 + t167 * t174 + t169 * t181 + t182 * t60) * t254;
t143 = (pkin(3) * t322 + t320 * t60) * t300;
t13 = t273 * t141 + t143 * t326 + t325 * t60;
t286 = t13 - t44;
t14 = t274 * t141 + t143 * t327 + t324 * t60;
t285 = t14 + t43;
t168 = t170 / 0.2e1;
t145 = (t168 * t181 + t174 * t315 + t182 * t314 - t339) * t254;
t147 = (pkin(3) * t323 + t314 * t320) * t300;
t17 = t273 * t145 + t147 * t326 + t314 * t325;
t284 = t17 - t44;
t18 = t274 * t145 + t147 * t327 + t314 * t324;
t283 = t18 + t43;
t282 = -t105 * t43 - t106 * t44;
t107 = t114 * t132;
t108 = t115 * t132;
t31 = -t107 * t44 - t108 * t43;
t281 = t114 * t43 - t115 * t44;
t280 = -rSges(4,1) * t105 - rSges(4,2) * t106;
t279 = -rSges(4,1) * t108 - rSges(4,2) * t107;
t268 = rSges(4,1) * t114 - rSges(4,2) * t115;
t112 = t114 * pkin(4);
t125 = t296 * pkin(1);
t267 = t112 + t125;
t124 = t132 * pkin(13);
t266 = t125 * t132 + t124;
t255 = t310 * t270;
t251 = t296 * t84;
t250 = t293 * t83;
t248 = pkin(4) * t107 + t266;
t247 = t129 * t296;
t246 = t129 * t293;
t130 = sin(pkin(15));
t245 = t130 * t262;
t241 = pkin(1) * t246;
t240 = t132 * t293 * pkin(1);
t134 = 0.1e1 / pkin(6);
t235 = t134 * t99 * t258;
t229 = t130 * t235;
t133 = cos(pkin(15));
t228 = t133 * t235;
t100 = t105 * pkin(4);
t227 = -t100 - t241;
t102 = t108 * pkin(4);
t226 = -t102 - t240;
t109 = pkin(1) * t222 - pkin(5);
t94 = t104 + t265;
t90 = -pkin(1) * t272 - t109 * t94;
t91 = -t109 * t136 + t291 * t94;
t225 = -pkin(7) + (rSges(7,1) * (t133 * t262 * t90 + t245 * t91) + rSges(7,2) * (-t91 * t307 * t133 / 0.2e1 + t90 * t245)) * t134;
t224 = rSges(6,1) * t131 - rSges(6,2) * t128 + pkin(9);
t223 = (-t125 - pkin(13)) * t129;
t221 = g(1) * t224;
t220 = g(2) * t224;
t219 = g(3) * t224;
t217 = rSges(3,1) * t296 - rSges(3,2) * t293;
t72 = t250 + t251;
t214 = -pkin(4) * t106 + t223;
t210 = -t293 * t61 + t296 * t62 - t251;
t59 = -t250 + t210;
t206 = pkin(1) * t134 * t307 * (t215 + (0.2e1 * t109 * pkin(5) + t252 - t94) * t116);
t205 = t134 * (t109 * t252 + (pkin(5) * t260 - t222 * t94 - t272) * pkin(1)) * t262;
t183 = pkin(3) * t312 * t333;
t160 = t163 * t270;
t148 = (pkin(4) * t323 + t314 * t321) * t301;
t146 = (-t168 * t180 + t173 * t315 + t183 * t314 - t338) * t255;
t144 = (pkin(4) * t322 + t321 * t60) * t301;
t142 = (-pkin(4) * t276 + t167 * t173 - t169 * t180 + t183 * t60) * t255;
t80 = t83 * t246;
t70 = t72 * t132;
t69 = t73 * t132;
t68 = -t247 * t84 - t80;
t67 = -t246 * t84 + t247 * t83;
t64 = t133 * t206 / 0.2e1 + t90 * t228 + t130 * t205 + t91 * t229;
t63 = t133 * t205 + t91 * t228 - t130 * t206 / 0.2e1 - t90 * t229;
t57 = t58 * t132;
t56 = t59 * t132;
t55 = t58 * t129;
t54 = t129 * t210 - t80;
t46 = t127 * t255 * t53 + t126 * t160;
t45 = t126 * t269 * t299 + t127 * t160;
t30 = t107 * t43 - t108 * t44;
t26 = t128 * t129 + t131 * t30;
t25 = -t128 * t30 + t129 * t131;
t20 = t126 * t146 + t148 * t329 + t314 * t319;
t19 = t127 * t146 + t148 * t330 + t314 * t318;
t16 = t126 * t142 + t144 * t329 + t319 * t60;
t15 = t127 * t142 + t144 * t330 + t318 * t60;
t12 = t114 * t284 - t115 * t283;
t11 = t114 * t18 + t115 * t17 + t281;
t10 = -t107 * t283 - t108 * t284;
t9 = t107 * t17 - t108 * t18 + t31;
t8 = -t105 * t284 - t106 * t283;
t7 = -t105 * t18 + t106 * t17 + t282;
t6 = t114 * t286 - t115 * t285;
t5 = t114 * t14 + t115 * t13 + t281;
t4 = -t107 * t285 - t108 * t286;
t3 = t107 * t13 - t108 * t14 + t31;
t2 = -t105 * t286 - t106 * t285;
t1 = -t105 * t14 + t106 * t13 + t282;
t21 = [-m(2) * (g(1) * (-rSges(2,1) * t129 - rSges(2,2) * t132) + g(2) * (rSges(2,1) * t132 - rSges(2,2) * t129)) - m(3) * (g(2) * t124 + (rSges(3,3) * g(1) + g(2) * t217) * t132 + (g(1) * (-pkin(13) - t217) + g(2) * rSges(3,3)) * t129) - m(4) * (g(1) * (-rSges(4,1) * t106 + rSges(4,2) * t105 + rSges(4,3) * t132 + t223) + g(2) * (rSges(4,1) * t107 - rSges(4,2) * t108 + rSges(4,3) * t129 + t266)) - m(5) * (g(1) * (-rSges(5,1) * t29 - rSges(5,2) * t282 + t132 * rSges(5,3) + t214) + g(2) * (rSges(5,1) * t30 + rSges(5,2) * t31 + rSges(5,3) * t129 + t248)) - m(6) * (g(1) * (-t29 * pkin(9) + rSges(6,1) * t335 + rSges(6,2) * t334 - t298 * t282 + t214) + g(2) * (pkin(9) * t30 + rSges(6,1) * t26 + rSges(6,2) * t25 + t298 * t31 + t248)) - m(7) * ((rSges(7,3) * g(1) + g(2) * t225) * t132 + (rSges(7,3) * g(2) - g(1) * t225) * t129) - m(8) * (g(1) * (-rSges(8,1) * t67 - rSges(8,2) * t68 + rSges(8,3) * t132 + t223) + g(2) * (rSges(8,1) * t69 - rSges(8,2) * t70 + rSges(8,3) * t129 + t266)) - m(9) * (g(1) * (-pkin(1) * t247 - t129 * pkin(13) + (-t45 * t68 - t46 * t67) * rSges(9,1) + (t45 * t67 - t46 * t68) * rSges(9,2) + t132 * rSges(9,3)) + g(2) * ((-t45 * t70 + t46 * t69) * rSges(9,1) + (-t45 * t69 - t46 * t70) * rSges(9,2) + t129 * rSges(9,3) + t266) + (g(1) * (t126 * t68 - t127 * t67) + g(2) * (t126 * t70 + t127 * t69)) * pkin(3)), -m(3) * (g(3) * t217 + t317 * (-rSges(3,1) * t293 - rSges(3,2) * t296)) - m(4) * (g(1) * (-t240 + t279) + g(2) * (-t241 + t280) + g(3) * (t125 + t268)) - m(5) * (g(1) * (rSges(5,1) * t3 + rSges(5,2) * t4 + t226) + g(2) * (rSges(5,1) * t1 + rSges(5,2) * t2 + t227) + g(3) * (rSges(5,1) * t5 + rSges(5,2) * t6 + t267)) - m(6) * (g(1) * (t298 * t4 + t226) + g(2) * (t2 * t298 + t227) + g(3) * (t298 * t6 + t267) + t5 * t219 + t3 * t221 + t1 * t220) - m(7) * (g(3) * (rSges(7,1) * t63 + rSges(7,2) * t64) + t317 * (rSges(7,1) * t64 - rSges(7,2) * t63)) - m(8) * (g(1) * (rSges(8,1) * t56 - rSges(8,2) * t57 - t240) + g(2) * (rSges(8,1) * t54 - rSges(8,2) * t55 - t241) + g(3) * (rSges(8,1) * t58 + rSges(8,2) * t59 + t125)) - m(9) * (g(1) * (-t240 + (-t15 * t70 + t16 * t69 - t45 * t57 + t46 * t56) * rSges(9,1) + (-t15 * t69 - t16 * t70 - t45 * t56 - t46 * t57) * rSges(9,2)) + g(2) * (-t241 + (t15 * t68 + t16 * t67 - t45 * t55 + t46 * t54) * rSges(9,1) + (-t15 * t67 + t16 * t68 - t45 * t54 - t46 * t55) * rSges(9,2)) + g(3) * (t125 + (t15 * t73 + t16 * t72 + t45 * t59 + t46 * t58) * rSges(9,1) + (-t15 * t72 + t16 * t73 - t45 * t58 + t46 * t59) * rSges(9,2)) + (g(1) * (t126 * t57 + t127 * t56) + g(2) * (t126 * t55 + t127 * t54) + g(3) * (-t126 * t59 + t127 * t58)) * pkin(3)), -m(4) * (g(1) * t279 + g(2) * t280 + g(3) * t268) - m(5) * (g(1) * (rSges(5,1) * t9 + rSges(5,2) * t10 - t102) + g(2) * (rSges(5,1) * t7 + rSges(5,2) * t8 - t100) + g(3) * (rSges(5,1) * t11 + rSges(5,2) * t12 + t112)) - m(6) * (g(1) * (t10 * t298 - t102) + g(2) * (t298 * t8 - t100) + g(3) * (t12 * t298 + t112) + t9 * t221 + t7 * t220 + t11 * t219) - m(9) * (g(1) * ((-t19 * t70 + t20 * t69) * rSges(9,1) + (-t19 * t69 - t20 * t70) * rSges(9,2)) + g(2) * ((t19 * t68 + t20 * t67) * rSges(9,1) + (-t19 * t67 + t20 * t68) * rSges(9,2)) + g(3) * ((t19 * t73 + t20 * t72) * rSges(9,1) + (-t19 * t72 + t20 * t73) * rSges(9,2))), -m(6) * (g(1) * (rSges(6,1) * t25 - rSges(6,2) * t26) + g(2) * (-rSges(6,1) * t334 + rSges(6,2) * t335) + g(3) * (-rSges(6,1) * t128 - rSges(6,2) * t131) * (t114 * t44 + t115 * t43))];
taug = t21(:);
