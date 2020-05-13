% Calculate Gravitation load on the joints for
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1DE1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:30:11
% EndTime: 2020-04-18 10:42:05
% DurationCPUTime: 119.08s
% Computational Cost: add. (3241904->361), mult. (4894387->612), div. (216516->22), fcn. (3103493->36), ass. (0->271)
t168 = sin(qJ(4));
t173 = cos(qJ(4));
t176 = cos(qJ(1));
t170 = sin(qJ(1));
t169 = sin(qJ(2));
t174 = cos(qJ(3));
t175 = cos(qJ(2));
t309 = sin(qJ(3));
t216 = t169 * t174 + t175 * t309;
t137 = t216 * t170;
t252 = t169 * t309;
t148 = -t174 * t175 + t252;
t138 = t148 * t170;
t163 = cos(pkin(17));
t178 = 0.1e1 / pkin(10);
t183 = pkin(3) ^ 2;
t182 = pkin(4) ^ 2;
t181 = pkin(5) ^ 2;
t186 = pkin(1) ^ 2;
t171 = sin(pkin(16));
t310 = cos(pkin(16));
t149 = t169 * t171 - t175 * t310;
t308 = pkin(5) * t149;
t326 = -2 * pkin(1);
t277 = t308 * t326 + t186;
t136 = t181 + t277;
t273 = pkin(2) ^ 2 - pkin(6) ^ 2;
t125 = t136 + t273;
t141 = pkin(1) - t308;
t321 = -pkin(6) - pkin(2);
t122 = (pkin(5) - t321) * (pkin(5) + t321) + t277;
t322 = pkin(2) - pkin(6);
t123 = (pkin(5) - t322) * (pkin(5) + t322) + t277;
t187 = sqrt(-t123 * t122);
t151 = t169 * t310 + t171 * t175;
t307 = pkin(5) * t151;
t117 = t125 * t307 + t141 * t187;
t253 = t117 * t309;
t130 = 0.1e1 / t136;
t185 = 0.1e1 / pkin(2);
t287 = t130 * t185;
t285 = t151 * t187;
t116 = -pkin(5) * t285 + t125 * t141;
t293 = t116 * t174;
t110 = (-t293 / 0.2e1 + t253 / 0.2e1) * t287;
t254 = t116 * t309;
t292 = t117 * t174;
t111 = (t292 / 0.2e1 + t254 / 0.2e1) * t287;
t161 = pkin(18) + pkin(19);
t157 = sin(t161);
t158 = cos(t161);
t89 = t110 * t158 + t111 * t157;
t215 = pkin(4) * t89;
t294 = 0.2e1 * pkin(3) * t215 + t182;
t85 = t183 + t294;
t83 = 0.1e1 / t85;
t296 = t178 * t83;
t162 = sin(pkin(17));
t314 = t162 / 0.2e1;
t320 = -pkin(8) - pkin(10);
t80 = (pkin(3) - t320) * (pkin(3) + t320) + t294;
t319 = pkin(10) - pkin(8);
t81 = (pkin(3) - t319) * (pkin(3) + t319) + t294;
t188 = sqrt(-t81 * t80);
t225 = t110 * t157 - t111 * t158;
t295 = t188 * t225;
t272 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t82 = t85 + t272;
t86 = pkin(3) * t89 + pkin(4);
t65 = -pkin(3) * t295 + t82 * t86;
t317 = pkin(3) * t225;
t66 = t188 * t86 + t317 * t82;
t59 = (-t65 * t163 / 0.2e1 + t66 * t314) * t296;
t60 = (t66 * t163 / 0.2e1 + t65 * t314) * t296;
t47 = atan2(t60, t59);
t45 = cos(t47);
t35 = t138 * t45;
t44 = sin(t47);
t26 = -t137 * t44 - t35;
t338 = t176 * t168 + t173 * t26;
t337 = t168 * t26 - t173 * t176;
t323 = pkin(3) * pkin(4);
t336 = 0.2e1 * t323;
t71 = 0.1e1 / t188;
t335 = -t71 / 0.2e1;
t271 = 0.1e1 / t85 ^ 2 * t323;
t246 = t163 * t271;
t270 = pkin(1) * t307;
t291 = 0.2e1 / t187 * (t122 + t123) * t270;
t251 = -t291 / 0.2e1;
t286 = t149 * t187;
t102 = (t286 + (t141 * t326 - t125 + t251) * t151) * pkin(5);
t245 = 0.1e1 / t136 ^ 2 * t270;
t266 = t309 / 0.2e1;
t324 = -0.2e1 * t151 ^ 2;
t104 = t141 * t291 / 0.2e1 + t181 * pkin(1) * t324 + (-t125 * t149 - t285) * pkin(5);
t315 = t104 / 0.2e1;
t78 = ((t102 * t266 + t174 * t315) * t130 + (t254 + t292) * t245) * t185;
t316 = -t102 / 0.2e1;
t79 = ((t104 * t266 + t174 * t316) * t130 + (t253 - t293) * t245) * t185;
t69 = -t157 * t78 - t158 * t79;
t241 = t69 * t246;
t247 = t162 * t271;
t243 = t69 * t247;
t301 = t163 * t83;
t257 = t301 / 0.2e1;
t258 = -t301 / 0.2e1;
t259 = t83 * t314;
t58 = 0.1e1 / t59 ^ 2;
t297 = t178 / (t58 * t60 ^ 2 + 0.1e1);
t304 = t58 * t60;
t262 = t225 * t335;
t224 = (t80 + t81) * t336;
t62 = t69 * t224;
t68 = t157 * t79 - t158 * t78;
t219 = -t188 * t68 + t262 * t62;
t325 = -0.2e1 * pkin(4);
t255 = t325 * t86 - t82;
t31 = (t255 * t69 + t219) * pkin(3);
t248 = t183 * t225 * t325;
t263 = t71 * t86 / 0.2e1;
t32 = t62 * t263 + t69 * t248 + (-t188 * t69 + t68 * t82) * pkin(3);
t57 = 0.1e1 / t59;
t15 = ((t241 * t66 + t243 * t65 + t257 * t32 + t259 * t31) * t57 - (-t241 * t65 + t243 * t66 + t258 * t31 + t259 * t32) * t304) * t297;
t334 = t15 + 0.1e1;
t240 = t225 * t246;
t242 = t225 * t247;
t67 = t225 * t224;
t218 = -t188 * t89 + t262 * t67;
t54 = (t225 * t255 + t218) * pkin(3);
t55 = t67 * t263 + t225 * t248 + (t82 * t89 - t295) * pkin(3);
t18 = ((t240 * t66 + t242 * t65 + t257 * t55 + t259 * t54) * t57 - (-t240 * t65 + t242 * t66 + t258 * t54 + t259 * t55) * t304) * t297;
t333 = t18 + 0.1e1;
t311 = pkin(11) + rSges(6,3);
t303 = t83 * t85;
t166 = cos(pkin(19));
t164 = sin(pkin(19));
t313 = t164 / 0.2e1;
t107 = (-t116 * t166 / 0.2e1 + t117 * t313) * t287;
t108 = (t117 * t166 / 0.2e1 + t116 * t313) * t287;
t99 = atan2(t108, t107);
t93 = sin(t99);
t298 = t169 * t93;
t94 = cos(t99);
t92 = t175 * t94;
t76 = t92 - t298;
t211 = -t215 - pkin(3);
t230 = t85 - t272;
t207 = -pkin(4) * t295 - t211 * t230;
t206 = 0.1e1 / t207 ^ 2;
t220 = pkin(4) * t230;
t208 = -t188 * t211 + t220 * t225;
t205 = 0.1e1 / (t206 * t208 ^ 2 + 0.1e1);
t199 = t205 * t206 * t208;
t200 = t205 / t207;
t209 = t211 * t336 - t220;
t210 = t211 * t335;
t223 = -pkin(4) * t188 - 0.2e1 * t182 * t317;
t332 = (t210 * t67 + t220 * t89 + t223 * t225) * t200 - (pkin(4) * t218 + t209 * t225) * t199;
t331 = (t210 * t62 + t220 * t68 + t223 * t69) * t200 - (pkin(4) * t219 + t209 * t69) * t199;
t165 = sin(pkin(18));
t167 = cos(pkin(18));
t267 = t83 / pkin(8) / 0.2e1;
t204 = atan2(t208 * t267, t207 * t267);
t203 = sin(t204);
t61 = cos(t204);
t330 = -t165 * t61 + t167 * t203;
t43 = -t165 * t203 - t167 * t61;
t329 = g(1) * t176 + g(2) * t170;
t328 = t303 * t43;
t327 = t330 * t303;
t105 = 0.1e1 / t107 ^ 2;
t235 = t117 * t245;
t236 = t116 * t245;
t250 = t130 * t313;
t290 = t130 * t166;
t63 = ((t102 * t250 + t164 * t236 + t166 * t235 + t290 * t315) / t107 - (t104 * t250 + t164 * t235 - t166 * t236 + t290 * t316) * t108 * t105) / (t105 * t108 ^ 2 + 0.1e1) * t185;
t318 = t63 + 0.1e1;
t177 = cos(pkin(15));
t312 = t177 / 0.2e1;
t160 = t175 * pkin(1);
t280 = t175 * t176;
t139 = t174 * t280 - t176 * t252;
t140 = t216 * t176;
t227 = t139 * t44 + t140 * t45;
t226 = t148 * t45 + t216 * t44;
t33 = t137 * t45;
t37 = t139 * t45;
t41 = t216 * t45;
t172 = sin(pkin(15));
t289 = t130 * t172;
t180 = 0.1e1 / pkin(6);
t288 = t130 * t180;
t284 = t169 * t170;
t283 = t170 * t175;
t279 = rSges(4,1) * t137 - rSges(4,2) * t138;
t278 = rSges(4,1) * t140 + rSges(4,2) * t139;
t276 = rSges(4,1) * t148 + rSges(4,2) * t216;
t146 = t148 * pkin(4);
t275 = t146 + t160;
t159 = t176 * pkin(13);
t274 = pkin(1) * t280 + t159;
t269 = pkin(1) * t284;
t268 = t176 * t169 * pkin(1);
t265 = t334 * t44;
t264 = t333 * t44;
t261 = t169 * t318;
t260 = t175 * t318;
t256 = -pkin(4) * t139 + t274;
t249 = t130 * t312;
t244 = t94 * t260;
t132 = t137 * pkin(4);
t239 = t132 - t269;
t134 = t140 * pkin(4);
t238 = t134 - t268;
t124 = t136 - t273;
t142 = pkin(1) * t149 - pkin(5);
t115 = -pkin(1) * t285 - t124 * t142;
t237 = t115 * t245;
t118 = pkin(1) * t151 * t124 - t142 * t187;
t234 = t118 * t245;
t109 = (t115 * t312 + t118 * t172 / 0.2e1) * t288;
t112 = (t118 * t312 - t115 * t172 / 0.2e1) * t288;
t100 = atan2(t112, t109);
t96 = sin(t100);
t97 = cos(t100);
t233 = rSges(7,1) * t97 - rSges(7,2) * t96;
t231 = (-pkin(13) - t160) * t170;
t229 = rSges(3,1) * t175 - rSges(3,2) * t169;
t77 = t169 * t94 + t175 * t93;
t222 = -pkin(7) + t233;
t221 = rSges(6,1) * t173 - rSges(6,2) * t168 + pkin(9);
t217 = -pkin(4) * t138 + t231;
t214 = g(1) * t221;
t213 = g(2) * t221;
t212 = g(3) * t221;
t52 = -t260 * t93 - t261 * t94;
t106 = 0.1e1 / t109 ^ 2;
t103 = t142 * t251 + t186 * pkin(5) * t324 + (-t124 * t149 - t285) * pkin(1);
t101 = (t286 + (0.2e1 * pkin(5) * t142 - t124 + t251) * t151) * pkin(1);
t91 = t93 * t284;
t75 = t76 * t176;
t74 = t77 * t176;
t73 = t283 * t94 - t91;
t72 = t77 * t170;
t53 = -t261 * t93 + t63 * t92 + t92;
t51 = t52 * t176;
t50 = (t298 * t318 - t244) * t176;
t49 = t52 * t170;
t48 = t91 + (t298 * t63 - t244) * t170;
t29 = t140 * t44 - t37;
t25 = -t138 * t44 + t33;
t24 = t168 * t170 + t173 * t29;
t23 = -t168 * t29 + t170 * t173;
t17 = t332 * t327;
t16 = t332 * t328;
t14 = t331 * t327;
t13 = t331 * t328;
t12 = t333 * t226;
t11 = t148 * t264 - t18 * t41 - t41;
t10 = t333 * t227;
t9 = t140 * t264 - t18 * t37 - t37;
t8 = -t138 * t264 + t18 * t33 + t33;
t7 = t137 * t264 + t18 * t35 + t35;
t6 = t334 * t226;
t5 = t148 * t265 - t15 * t41 - t41;
t4 = t334 * t227;
t3 = t140 * t265 - t15 * t37 - t37;
t2 = -t138 * t265 + t15 * t33 + t33;
t1 = t137 * t265 + t15 * t35 + t35;
t19 = [-m(2) * (g(1) * (-rSges(2,1) * t170 - rSges(2,2) * t176) + g(2) * (rSges(2,1) * t176 - rSges(2,2) * t170)) - m(3) * (g(2) * t159 + (rSges(3,3) * g(1) + g(2) * t229) * t176 + (g(1) * (-pkin(13) - t229) + g(2) * rSges(3,3)) * t170) - m(4) * (g(1) * (-rSges(4,1) * t138 - rSges(4,2) * t137 + rSges(4,3) * t176 + t231) + g(2) * (-rSges(4,1) * t139 + rSges(4,2) * t140 + rSges(4,3) * t170 + t274)) - m(5) * (g(1) * (rSges(5,1) * t26 - rSges(5,2) * t25 + rSges(5,3) * t176 + t217) + g(2) * (rSges(5,1) * t29 + rSges(5,2) * t227 + rSges(5,3) * t170 + t256)) - m(6) * (g(1) * (t26 * pkin(9) + rSges(6,1) * t338 - rSges(6,2) * t337 + t311 * t25 + t217) + g(2) * (pkin(9) * t29 + rSges(6,1) * t24 + rSges(6,2) * t23 - t227 * t311 + t256)) - m(7) * ((rSges(7,3) * g(1) + g(2) * t222) * t176 + (rSges(7,3) * g(2) - g(1) * t222) * t170) - m(8) * (g(1) * (-rSges(8,1) * t73 + rSges(8,2) * t72 + rSges(8,3) * t176 + t231) + g(2) * (rSges(8,1) * t75 - rSges(8,2) * t74 + rSges(8,3) * t170 + t274)) - m(9) * (g(1) * (-pkin(1) * t283 - t170 * pkin(13) + (-t330 * t72 - t43 * t73) * rSges(9,1) + (-t330 * t73 + t43 * t72) * rSges(9,2) + t176 * rSges(9,3)) + g(2) * ((t330 * t74 + t43 * t75) * rSges(9,1) + (t330 * t75 - t43 * t74) * rSges(9,2) + t170 * rSges(9,3) + t274) + (g(1) * (-t165 * t72 - t167 * t73) + g(2) * (t165 * t74 + t167 * t75)) * pkin(3)), -m(3) * (g(3) * t229 + t329 * (-rSges(3,1) * t169 - rSges(3,2) * t175)) - m(4) * (g(1) * (-t268 + t278) + g(2) * (-t269 + t279) + g(3) * (t160 + t276)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + t238) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t1 + t239) + g(3) * (rSges(5,1) * t6 - rSges(5,2) * t5 + t275)) - m(6) * (g(1) * (t3 * t311 + t238) + g(2) * (t1 * t311 + t239) + g(3) * (t311 * t5 + t275) + t6 * t212 + t4 * t214 + t2 * t213) - m(7) * (g(3) * t233 + t329 * (-rSges(7,1) * t96 - rSges(7,2) * t97)) * ((t103 * t249 + t177 * t234 - t101 * t289 / 0.2e1 - t172 * t237) / t109 - (t101 * t249 + t177 * t237 + t103 * t289 / 0.2e1 + t172 * t234) * t112 * t106) / (t106 * t112 ^ 2 + 0.1e1) * t180 - m(8) * (g(1) * (rSges(8,1) * t51 + rSges(8,2) * t50 - t268) + g(2) * (rSges(8,1) * t49 + rSges(8,2) * t48 - t269) + g(3) * (rSges(8,1) * t53 + rSges(8,2) * t52 + t160)) - m(9) * (g(1) * (-t268 + (-t13 * t74 + t14 * t75 - t330 * t50 + t43 * t51) * rSges(9,1) + (-t13 * t75 - t14 * t74 + t330 * t51 + t43 * t50) * rSges(9,2)) + g(2) * (-t269 + (-t13 * t72 + t14 * t73 - t330 * t48 + t43 * t49) * rSges(9,1) + (-t13 * t73 - t14 * t72 + t330 * t49 + t43 * t48) * rSges(9,2)) + g(3) * (t160 + (t13 * t76 + t14 * t77 - t330 * t52 + t43 * t53) * rSges(9,1) + (-t13 * t77 + t14 * t76 + t330 * t53 + t43 * t52) * rSges(9,2)) + (g(1) * (-t165 * t50 + t167 * t51) + g(2) * (-t165 * t48 + t167 * t49) + g(3) * (-t165 * t52 + t167 * t53)) * pkin(3)), -m(4) * (g(1) * t278 + g(2) * t279 + g(3) * t276) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t134) + g(2) * (rSges(5,1) * t8 - rSges(5,2) * t7 + t132) + g(3) * (rSges(5,1) * t12 - rSges(5,2) * t11 + t146)) - m(6) * (g(1) * (t311 * t9 + t134) + g(2) * (t311 * t7 + t132) + g(3) * (t11 * t311 + t146) + t8 * t213 + t12 * t212 + t10 * t214) - m(9) * (g(1) * ((-t16 * t74 + t17 * t75) * rSges(9,1) + (-t16 * t75 - t17 * t74) * rSges(9,2)) + g(2) * ((-t16 * t72 + t17 * t73) * rSges(9,1) + (-t16 * t73 - t17 * t72) * rSges(9,2)) + g(3) * ((t16 * t76 + t17 * t77) * rSges(9,1) + (-t16 * t77 + t17 * t76) * rSges(9,2))), -m(6) * (g(1) * (rSges(6,1) * t23 - rSges(6,2) * t24) + g(2) * (rSges(6,1) * t337 + rSges(6,2) * t338) + g(3) * (-rSges(6,1) * t168 - rSges(6,2) * t173) * (t148 * t44 - t41))];
taug = t19(:);
