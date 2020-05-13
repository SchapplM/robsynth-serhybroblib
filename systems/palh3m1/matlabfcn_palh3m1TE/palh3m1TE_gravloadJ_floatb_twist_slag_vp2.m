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
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1TE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1TE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:22:20
% EndTime: 2020-04-17 15:27:29
% DurationCPUTime: 48.11s
% Computational Cost: add. (941691->327), mult. (1427744->522), div. (59424->10), fcn. (908093->24), ass. (0->233)
t291 = sin(qJ(2));
t292 = sin(pkin(16));
t294 = cos(qJ(2));
t295 = cos(pkin(16));
t224 = t291 * t292 - t294 * t295;
t221 = pkin(5) * t224;
t214 = (-0.2e1 * t221 + pkin(1)) * pkin(1);
t308 = pkin(5) ^ 2;
t104 = t214 + t308;
t139 = pkin(2) ^ 2;
t264 = pkin(6) ^ 2 - t139;
t211 = t104 - t264;
t216 = -t221 + pkin(1);
t116 = t291 * t295 + t292 * t294;
t301 = pkin(5) + pkin(6);
t302 = pkin(5) - pkin(6);
t138 = sqrt(-((-pkin(2) + t301) * (pkin(2) + t302) + t214) * ((pkin(2) + t301) * (-pkin(2) + t302) + t214));
t270 = t116 * t138;
t255 = pkin(5) * t270;
t310 = 0.1e1 / pkin(2);
t205 = t310 * (t211 * t216 - t255);
t304 = 0.1e1 / t104;
t204 = t304 * t205;
t202 = -t204 / 0.2e1;
t209 = pkin(5) * t211;
t206 = t310 * (t116 * t209 + t138 * t216);
t261 = t304 / 0.2e1;
t203 = t206 * t261;
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t197 = t202 * t293 + t203 * t290;
t201 = t204 / 0.2e1;
t198 = t201 * t290 + t203 * t293;
t258 = pkin(18) + pkin(19);
t242 = sin(t258);
t243 = cos(t258);
t314 = t197 * t243 + t242 * t198;
t189 = pkin(3) * t314;
t181 = (0.2e1 * t189 + pkin(4)) * pkin(4);
t299 = pkin(10) - pkin(8);
t300 = -pkin(8) - pkin(10);
t137 = sqrt(-((pkin(3) - t299) * (pkin(3) + t299) + t181) * ((pkin(3) - t300) * (pkin(3) + t300) + t181));
t313 = t197 * t242 - t243 * t198;
t340 = t313 * t137;
t344 = pkin(3) * t340;
t343 = pkin(4) * t340;
t141 = pkin(8) ^ 2;
t263 = pkin(10) ^ 2 - t141;
t188 = pkin(4) * t314;
t309 = pkin(4) ^ 2;
t66 = t309 + (0.2e1 * t188 + pkin(3)) * pkin(3);
t178 = t66 + t263;
t183 = t189 + pkin(4);
t306 = 0.1e1 / pkin(10);
t163 = t306 * (t178 * t183 - t344);
t176 = pkin(3) * t178;
t164 = t306 * (t137 * t183 + t176 * t313);
t271 = sin(pkin(17));
t272 = cos(pkin(17));
t342 = t271 * t163 + t272 * t164;
t127 = t294 * pkin(1);
t341 = pkin(13) + t127;
t275 = sin(pkin(19));
t276 = cos(pkin(19));
t83 = t202 * t276 + t203 * t275;
t84 = t201 * t275 + t203 * t276;
t73 = -t291 * t84 + t294 * t83;
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t134 = cos(qJ(1));
t115 = -t290 * t294 - t291 * t293;
t131 = sin(qJ(1));
t105 = t115 * t131;
t114 = t291 * t290 - t294 * t293;
t106 = t114 * t131;
t158 = t272 * t163;
t159 = t271 * t164;
t305 = 0.1e1 / t66;
t267 = t305 / 0.2e1;
t268 = -t305 / 0.2e1;
t43 = t158 * t268 + t159 * t267;
t44 = t342 * t267;
t29 = -t105 * t44 + t106 * t43;
t339 = t134 * t130 - t133 * t29;
t338 = t130 * t29 + t133 * t134;
t337 = -0.2e1 * t313;
t336 = pkin(3) * m(9);
t218 = t224 * t138;
t289 = pkin(1) * t116;
t257 = pkin(5) * t289;
t285 = 0.4e1 / t138 * (t301 * t302 - t139 + t214) * t257;
t251 = -t285 / 0.2e1;
t78 = (t218 + (t251 + t264 - t308) * t116) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t221) * t257;
t335 = -t78 / 0.2e1;
t259 = -0.2e1 * pkin(1) * t116 ^ 2;
t79 = -t255 + t216 * t285 / 0.2e1 - t224 * t209 + t308 * t259;
t334 = t79 / 0.2e1;
t333 = -m(6) - m(5);
t128 = sin(pkin(18));
t332 = t128 / 0.2e1;
t129 = cos(pkin(18));
t331 = -t129 / 0.2e1;
t330 = -t313 / 0.2e1;
t329 = t271 / 0.2e1;
t328 = -t272 / 0.2e1;
t327 = -pkin(9) * m(6) - mrSges(6,1) * t133 + mrSges(6,2) * t130 - mrSges(5,1);
t199 = t205 * t257;
t200 = t206 * t257;
t252 = t310 * t261;
t229 = t275 * t252;
t288 = t310 * t304;
t237 = t276 * t288;
t99 = 0.1e1 / t104 ^ 2;
t61 = t237 * t334 + t78 * t229 + (t199 * t275 + t200 * t276) * t99;
t62 = t237 * t335 + t79 * t229 + (-t199 * t276 + t200 * t275) * t99;
t58 = t291 * t62 + t294 * t61 + t73;
t303 = pkin(3) * pkin(4);
t260 = 0.1e1 / t66 ^ 2 * t303;
t326 = (-t158 + t159) * t260;
t325 = t342 * t260;
t315 = 0.4e1 / t137 * ((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t181 - t141) * t303;
t172 = t313 * t315;
t324 = -t137 * t314 + t172 * t330;
t230 = t290 * t252;
t241 = t293 * t288;
t193 = t241 * t334 + t78 * t230 + (t199 * t290 + t200 * t293) * t99;
t194 = t241 * t335 + t79 * t230 + (-t199 * t293 + t200 * t290) * t99;
t169 = -t193 * t243 + t194 * t242;
t60 = -t193 * t242 - t194 * t243;
t173 = t60 * t315;
t323 = -t137 * t169 + t173 * t330;
t177 = t66 - t263;
t175 = pkin(4) * t177;
t182 = -t188 - pkin(3);
t322 = 0.2e1 * t182 * t303 - t175;
t321 = -0.2e1 * t183 * t303 - t176;
t307 = 0.1e1 / pkin(8);
t165 = t307 * (-t137 * t182 + t175 * t313);
t161 = t165 * t260;
t53 = -t177 * t182 - t343;
t296 = t53 * t307;
t238 = t260 * t296;
t320 = -t128 * t161 - t129 * t238;
t319 = t128 * t238 - t129 * t161;
t317 = m(9) + m(4) + m(8);
t135 = cos(pkin(15));
t136 = 0.1e1 / pkin(6);
t219 = mrSges(3,1) * t294 - mrSges(3,2) * t291;
t132 = sin(pkin(15));
t244 = t132 * t261;
t109 = pkin(1) * t224 - pkin(5);
t94 = t104 + t264;
t90 = -pkin(1) * t270 - t109 * t94;
t91 = -t109 * t138 + t289 * t94;
t316 = m(3) * pkin(13) + mrSges(2,1) + t219 - pkin(7) * m(7) - (-mrSges(7,1) * (t135 * t261 * t90 + t244 * t91) - mrSges(7,2) * (-t91 * t304 * t135 / 0.2e1 + t90 * t244)) * t136;
t312 = m(6) * pkin(11) - mrSges(5,2) + mrSges(6,3);
t311 = -mrSges(7,3) - mrSges(4,3) - mrSges(8,3) + mrSges(2,2) - mrSges(5,3) - mrSges(9,3) - mrSges(3,3);
t298 = t307 * t305;
t297 = t306 * t305;
t100 = t105 * pkin(4);
t108 = t115 * t134;
t102 = t108 * pkin(4);
t112 = t114 * pkin(4);
t171 = t173 / 0.2e1;
t184 = pkin(4) * pkin(3) ^ 2 * t337;
t253 = t306 * t267;
t274 = t137 * t60;
t143 = (-pkin(3) * t274 + t169 * t176 + t171 * t183 + t184 * t60) * t253;
t145 = (t323 * pkin(3) + t321 * t60) * t297;
t13 = t271 * t143 + t145 * t328 + t326 * t60;
t284 = t13 - t44;
t14 = t272 * t143 + t145 * t329 + t325 * t60;
t283 = t14 + t43;
t170 = t172 / 0.2e1;
t147 = (t170 * t183 + t176 * t314 + t184 * t313 - t344) * t253;
t149 = (t324 * pkin(3) + t313 * t321) * t297;
t17 = t271 * t147 + t149 * t328 + t313 * t326;
t282 = t17 - t44;
t18 = t272 * t147 + t149 * t329 + t313 * t325;
t281 = t18 + t43;
t280 = -t105 * t43 - t106 * t44;
t107 = t114 * t134;
t31 = -t107 * t44 - t108 * t43;
t279 = t114 * t43 - t115 * t44;
t278 = -t105 * mrSges(4,1) - t106 * mrSges(4,2);
t277 = -t108 * mrSges(4,1) - t107 * mrSges(4,2);
t266 = t114 * mrSges(4,1) - t115 * mrSges(4,2);
t118 = t341 * t134;
t254 = t307 * t268;
t250 = t294 * t84;
t249 = t291 * t83;
t247 = t107 * pkin(4) + t118;
t246 = t131 * t294;
t245 = t131 * t291;
t240 = pkin(1) * t245;
t239 = t291 * t134 * pkin(1);
t234 = t136 * t99 * t257;
t228 = t132 * t234;
t227 = t135 * t234;
t117 = t341 * t131;
t72 = t249 + t250;
t217 = -t106 * pkin(4) - t117;
t212 = -t291 * t61 + t294 * t62 - t250;
t59 = -t249 + t212;
t208 = pkin(1) * t136 * t304 * (t218 + (0.2e1 * t109 * pkin(5) + t251 - t94) * t116);
t207 = t136 * (t109 * t251 + (pkin(5) * t259 - t224 * t94 - t270) * pkin(1)) * t261;
t185 = pkin(3) * t309 * t337;
t162 = t165 * t268;
t150 = (t324 * pkin(4) + t313 * t322) * t298;
t148 = (-t170 * t182 + t175 * t314 + t185 * t313 - t343) * t254;
t146 = (t323 * pkin(4) + t322 * t60) * t298;
t144 = (-pkin(4) * t274 + t169 * t175 - t171 * t182 + t185 * t60) * t254;
t80 = t83 * t245;
t70 = t72 * t134;
t69 = t73 * t134;
t68 = -t246 * t84 - t80;
t67 = -t245 * t84 + t246 * t83;
t64 = t135 * t208 / 0.2e1 + t90 * t227 + t132 * t207 + t91 * t228;
t63 = t135 * t207 + t91 * t227 - t132 * t208 / 0.2e1 - t90 * t228;
t57 = t58 * t134;
t56 = t59 * t134;
t55 = t58 * t131;
t54 = t131 * t212 - t80;
t46 = t129 * t254 * t53 + t128 * t162;
t45 = t128 * t267 * t296 + t129 * t162;
t30 = t107 * t43 - t108 * t44;
t26 = t130 * t131 + t133 * t30;
t25 = -t130 * t30 + t131 * t133;
t20 = t128 * t148 + t150 * t331 + t313 * t320;
t19 = t129 * t148 + t150 * t332 + t313 * t319;
t16 = t128 * t144 + t146 * t331 + t320 * t60;
t15 = t129 * t144 + t146 * t332 + t319 * t60;
t1 = [(-t25 * mrSges(6,2) - mrSges(4,1) * t107 + mrSges(4,2) * t108 - mrSges(8,1) * t69 + mrSges(8,2) * t70 - m(6) * (pkin(9) * t30 + t247) - (t128 * t70 + t129 * t69) * t336 - t26 * mrSges(6,1) - m(5) * t247 - (-t45 * t70 + t46 * t69) * mrSges(9,1) - (-t45 * t69 - t46 * t70) * mrSges(9,2) - t30 * mrSges(5,1) + t312 * t31 - t317 * t118 - t316 * t134 + t311 * t131) * g(2) + (mrSges(4,1) * t106 - mrSges(4,2) * t105 + mrSges(8,1) * t67 + mrSges(8,2) * t68 - m(5) * t217 - t338 * mrSges(6,2) - (-t45 * t68 - t46 * t67) * mrSges(9,1) - (t45 * t67 - t46 * t68) * mrSges(9,2) - t339 * mrSges(6,1) - m(6) * (-pkin(9) * t29 + t217) + t29 * mrSges(5,1) - (t128 * t68 - t129 * t67) * t336 - t312 * t280 + t317 * t117 + t316 * t131 + t311 * t134) * g(1), (-(t128 * t55 + t129 * t54) * t336 - (t15 * t68 + t16 * t67 - t45 * t55 + t46 * t54) * mrSges(9,1) - (-t15 * t67 + t16 * t68 - t45 * t54 - t46 * t55) * mrSges(9,2) - t278 - t54 * mrSges(8,1) + t55 * mrSges(8,2) + t317 * t240 + t333 * (-t100 - t240) + t312 * (-t105 * t284 - t106 * t283) + t327 * (-t105 * t14 + t106 * t13 + t280)) * g(2) + (-(t128 * t57 + t129 * t56) * t336 - (-t15 * t70 + t16 * t69 - t45 * t57 + t46 * t56) * mrSges(9,1) - (-t15 * t69 - t16 * t70 - t45 * t56 - t46 * t57) * mrSges(9,2) - t277 - t56 * mrSges(8,1) + t57 * mrSges(8,2) + t333 * (-t102 - t239) + t312 * (-t107 * t283 - t108 * t284) + t327 * (t107 * t13 - t108 * t14 + t31) + t317 * t239) * g(1) + (-(-t128 * t59 + t129 * t58) * t336 - (t15 * t73 + t16 * t72 + t45 * t59 + t46 * t58) * mrSges(9,1) - (-t15 * t72 + t16 * t73 - t45 * t58 + t46 * t59) * mrSges(9,2) - t266 - t219 - mrSges(7,1) * t63 - mrSges(7,2) * t64 - t58 * mrSges(8,1) - t59 * mrSges(8,2) + t333 * (t112 + t127) + t312 * (t114 * t284 - t115 * t283) + t327 * (t114 * t14 + t115 * t13 + t279) - t317 * t127) * g(3) + (-mrSges(3,1) * t291 + mrSges(7,1) * t64 - mrSges(3,2) * t294 - mrSges(7,2) * t63) * (-g(1) * t134 - g(2) * t131), (-(t19 * t73 + t20 * t72) * mrSges(9,1) - (-t19 * t72 + t20 * t73) * mrSges(9,2) - t266 + t333 * t112 + t312 * (t114 * t282 - t115 * t281) + t327 * (t114 * t18 + t115 * t17 + t279)) * g(3) + (-(t19 * t68 + t20 * t67) * mrSges(9,1) - (-t19 * t67 + t20 * t68) * mrSges(9,2) - t278 - t333 * t100 + t312 * (-t105 * t282 - t106 * t281) + t327 * (-t105 * t18 + t106 * t17 + t280)) * g(2) + (-(-t19 * t70 + t20 * t69) * mrSges(9,1) - (-t19 * t69 - t20 * t70) * mrSges(9,2) - t277 + t327 * (t107 * t17 - t108 * t18 + t31) - t333 * t102 + t312 * (-t107 * t281 - t108 * t282)) * g(1), -g(2) * (-t338 * mrSges(6,1) + t339 * mrSges(6,2)) - g(1) * (mrSges(6,1) * t25 - mrSges(6,2) * t26) - g(3) * (-mrSges(6,1) * t130 - mrSges(6,2) * t133) * (t114 * t44 + t115 * t43)];
taug = t1(:);
