% Calculate joint inertia matrix for
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1DE2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_inertiaJ_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE2_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1DE2_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:29:00
% EndTime: 2020-04-14 20:40:16
% DurationCPUTime: 129.40s
% Computational Cost: add. (3632458->395), mult. (5482303->700), div. (243212->33), fcn. (3475166->44), ass. (0->318)
t396 = 2 * pkin(1);
t204 = sin(qJ(4));
t208 = cos(qJ(4));
t395 = Ifges(6,5) * t204 + Ifges(6,6) * t208;
t394 = 2 * pkin(4);
t199 = sin(pkin(20));
t203 = cos(pkin(20));
t205 = sin(qJ(3));
t209 = cos(qJ(3));
t180 = t199 * t209 + t203 * t205;
t352 = pkin(6) * t180;
t287 = pkin(1) * t352;
t174 = 0.2e1 * t287;
t228 = pkin(2) ^ 2;
t223 = pkin(6) ^ 2;
t230 = pkin(1) ^ 2;
t292 = t223 + t230;
t282 = -pkin(13) ^ 2 + t292;
t163 = t174 + t228 + t282;
t173 = -pkin(1) - t352;
t181 = t199 * t205 - t203 * t209;
t297 = t174 + t223;
t380 = -pkin(2) - pkin(13);
t157 = (pkin(1) - t380) * (pkin(1) + t380) + t297;
t375 = pkin(13) - pkin(2);
t158 = (pkin(1) - t375) * (pkin(1) + t375) + t297;
t312 = t158 * t157;
t231 = sqrt(-t312);
t308 = t181 * t231;
t286 = pkin(6) * t308;
t139 = -t163 * t173 - t286;
t351 = pkin(6) * t181;
t140 = t163 * t351 - t173 * t231;
t169 = t174 + t292;
t166 = 0.1e1 / t169;
t229 = 0.1e1 / pkin(2);
t360 = t229 / 0.2e1;
t275 = t166 * t360;
t126 = atan2(t140 * t275, t139 * t275);
t123 = sin(t126);
t124 = cos(t126);
t206 = sin(qJ(2));
t210 = cos(qJ(2));
t100 = -t123 * t206 + t124 * t210;
t162 = t228 - t282 - 0.2e1 * t287;
t273 = 0.1e1 / pkin(13) * t360;
t149 = atan2(t231 * t273, t162 * t273);
t147 = sin(t149);
t148 = cos(t149);
t99 = -t123 * t210 - t124 * t206;
t81 = t100 * t147 - t148 * t99;
t82 = -t100 * t148 - t147 * t99;
t393 = Ifges(10,5) * t82 + Ifges(10,6) * t81;
t224 = pkin(5) ^ 2;
t222 = pkin(7) ^ 2;
t211 = cos(pkin(19));
t359 = sin(pkin(19));
t184 = t206 * t211 - t210 * t359;
t350 = pkin(7) * t184;
t388 = -2 * pkin(1);
t296 = t350 * t388 + t230;
t172 = t222 + t296;
t291 = pkin(3) ^ 2 - pkin(8) ^ 2;
t165 = t172 + t291;
t175 = pkin(1) - t350;
t378 = pkin(7) + pkin(8);
t379 = pkin(7) - pkin(8);
t159 = (pkin(3) + t378) * (-pkin(3) + t379) + t296;
t160 = (-pkin(3) + t378) * (pkin(3) + t379) + t296;
t232 = sqrt(-t160 * t159);
t185 = t206 * t359 + t210 * t211;
t349 = pkin(7) * t185;
t143 = t165 * t349 + t175 * t232;
t299 = t209 * t143;
t306 = t185 * t232;
t142 = -pkin(7) * t306 + t165 * t175;
t302 = t205 * t142;
t170 = 0.1e1 / t172;
t227 = 0.1e1 / pkin(3);
t309 = t170 * t227;
t131 = (t302 / 0.2e1 + t299 / 0.2e1) * t309;
t300 = t209 * t142;
t301 = t205 * t143;
t132 = (-t300 / 0.2e1 + t301 / 0.2e1) * t309;
t193 = pkin(23) + pkin(22);
t190 = sin(t193);
t191 = cos(t193);
t107 = t131 * t191 - t132 * t190;
t355 = t107 * pkin(5);
t298 = t355 * t394 + t224;
t377 = (-pkin(9) - pkin(11));
t92 = ((pkin(4) - t377) * (pkin(4) + t377)) + t298;
t376 = (pkin(11) - pkin(9));
t93 = ((pkin(4) - t376) * (pkin(4) + t376)) + t298;
t233 = sqrt(-t93 * t92);
t264 = t131 * t190 + t191 * t132;
t391 = t264 * t233;
t390 = t205 * mrSges(4,1) + t209 * mrSges(4,2);
t389 = Ifges(9,5) * t100 + Ifges(9,6) * t99 + (-t147 * t81 + t148 * t82) * mrSges(10,3) * pkin(2);
t187 = t206 * pkin(1) - pkin(16);
t248 = t205 * t206 - t209 * t210;
t168 = pkin(5) * t248 + t187;
t386 = 0.2e1 * t168;
t385 = -0.2e1 * t185 ^ 2;
t384 = 0.2e1 * t187;
t202 = cos(pkin(21));
t225 = pkin(4) ^ 2;
t103 = t225 + t298;
t102 = 0.1e1 / t103 ^ 2;
t356 = pkin(5) * t102;
t288 = pkin(4) * t356;
t267 = t202 * t288;
t289 = pkin(1) * t349;
t313 = 0.2e1 / t232 * (t159 + t160) * t289;
t122 = t175 * t313 / 0.2e1 + t222 * pkin(1) * t385 + (-t165 * t184 - t306) * pkin(7);
t269 = 0.1e1 / t172 ^ 2 * t289;
t362 = -t209 / 0.2e1;
t276 = -t313 / 0.2e1;
t307 = t184 * t232;
t120 = (t307 + (t175 * t388 - t165 + t276) * t185) * pkin(7);
t369 = -t120 / 0.2e1;
t89 = ((t122 * t362 + t205 * t369) * t170 + (-t299 - t302) * t269) * t227;
t368 = t122 / 0.2e1;
t90 = ((t120 * t362 + t205 * t368) * t170 + (-t300 + t301) * t269) * t227;
t74 = t190 * t90 + t191 * t89;
t260 = t74 * t267;
t198 = sin(pkin(21));
t268 = t198 * t288;
t261 = t74 * t268;
t101 = 0.1e1 / t103;
t319 = t101 * t202;
t279 = t319 / 0.2e1;
t280 = -t319 / 0.2e1;
t365 = t198 / 0.2e1;
t281 = t101 * t365;
t217 = 0.1e1 / pkin(11);
t318 = t101 * t217;
t105 = pkin(4) * t107 + pkin(5);
t293 = pkin(9) ^ 2 - pkin(11) ^ 2;
t98 = t103 - t293;
t68 = -pkin(4) * t391 + t105 * t98;
t70 = pkin(4) * t264 * t98 + t105 * t233;
t56 = (-t68 * t202 / 0.2e1 + t70 * t365) * t318;
t54 = 0.1e1 / t56 ^ 2;
t55 = (t68 * t365 + t70 * t202 / 0.2e1) * t318;
t321 = t217 / (t54 * t55 ^ 2 + 0.1e1);
t337 = t54 * t55;
t84 = 0.1e1 / t233;
t381 = -t84 / 0.2e1;
t283 = t264 * t381;
t247 = pkin(5) * (t92 + t93) * t394;
t60 = t74 * t247;
t75 = -t190 * t89 + t191 * t90;
t241 = -t75 * t233 + t283 * t60;
t272 = -0.2e1 * pkin(5) * t105 - t98;
t40 = (t272 * t74 + t241) * pkin(4);
t354 = pkin(5) * t264;
t271 = -0.2e1 * t225 * t354;
t284 = t105 * t84 / 0.2e1;
t320 = t233 * t74;
t41 = t60 * t284 + t74 * t271 + (t75 * t98 - t320) * pkin(4);
t53 = 0.1e1 / t56;
t12 = 0.1e1 + ((t260 * t70 + t261 * t68 + t279 * t41 + t281 * t40) * t53 - (-t260 * t68 + t261 * t70 + t280 * t40 + t281 * t41) * t337) * t321;
t189 = pkin(1) * t205 + pkin(5);
t357 = pkin(1) * t209;
t49 = atan2(t55, t56);
t46 = sin(t49);
t47 = cos(t49);
t38 = t189 * t47 + t357 * t46;
t8 = -pkin(10) * t12 - t38;
t383 = m(6) * t8;
t183 = t205 * t210 + t206 * t209;
t36 = t183 * t46 + t248 * t47;
t382 = t36 / 0.2e1;
t219 = 0.1e1 / pkin(9);
t370 = t101 / 0.2e1;
t278 = t219 * t370;
t104 = -pkin(4) - t355;
t97 = t103 + t293;
t67 = -pkin(5) * t391 - t104 * t97;
t69 = -t104 * t233 + t354 * t97;
t58 = atan2(t69 * t278, t67 * t278);
t374 = sin(t58);
t200 = cos(pkin(23));
t304 = t200 * t142;
t196 = sin(pkin(23));
t305 = t196 * t143;
t128 = (-t304 / 0.2e1 + t305 / 0.2e1) * t309;
t127 = 0.1e1 / t128 ^ 2;
t303 = t200 * t143;
t315 = t142 * t196;
t129 = (t303 / 0.2e1 + t315 / 0.2e1) * t309;
t366 = t196 / 0.2e1;
t63 = 0.1e1 + (-((t122 * t366 + t200 * t369) * t170 + (-t304 + t305) * t269) * t129 * t127 + ((t120 * t366 + t200 * t368) * t170 + (t303 + t315) * t269) / t128) / (t127 * t129 ^ 2 + 0.1e1) * t227;
t373 = pkin(4) * t63;
t372 = pkin(5) * t46;
t371 = pkin(5) * t47;
t367 = t166 / 0.2e1;
t364 = -t204 / 0.2e1;
t363 = t208 / 0.2e1;
t212 = cos(pkin(18));
t361 = t212 / 0.2e1;
t358 = pkin(1) * t181;
t353 = pkin(6) * t140;
t348 = mrSges(5,1) * t38;
t39 = t46 * t189 - t357 * t47;
t347 = mrSges(5,2) * t39;
t346 = mrSges(5,3) * t36;
t343 = Ifges(6,6) * t36;
t256 = t264 * t267;
t257 = t264 * t268;
t71 = t264 * t247;
t240 = -t107 * t233 + t283 * t71;
t50 = (t264 * t272 + t240) * pkin(4);
t51 = t71 * t284 + t264 * t271 + (t107 * t98 - t391) * pkin(4);
t18 = 0.1e1 + ((t256 * t70 + t257 * t68 + t279 * t51 + t281 * t50) * t53 - (-t256 * t68 + t257 * t70 + t280 * t50 + t281 * t51) * t337) * t321;
t342 = Ifges(5,3) * t18;
t331 = Ifges(6,4) * t204;
t251 = Ifges(6,2) * t208 + t331;
t2 = t251 * t12;
t341 = t2 * t208;
t330 = Ifges(6,4) * t208;
t252 = Ifges(6,1) * t204 + t330;
t3 = t252 * t12;
t340 = t204 * t3;
t6 = t252 * t18;
t339 = t204 * t6;
t5 = t251 * t18;
t338 = t208 * t5;
t336 = t395 * t12;
t335 = t395 * t18;
t37 = t47 * t183 - t248 * t46;
t323 = t208 * t37;
t334 = Ifges(6,5) * t323 + Ifges(6,3) * t36;
t138 = 0.1e1 / t139 ^ 2;
t167 = 0.1e1 / t169 ^ 2;
t314 = 0.1e1 / t231 * (t157 + t158) * t351 * t396;
t277 = -t314 / 0.2e1;
t77 = 0.2e1 * (((t173 * t277 + (t163 * t180 - t308) * pkin(6)) * t367 + (-t166 * t181 * t223 + t167 * t353) * t358) / t139 - ((-t180 * t231 + (t277 - t163) * t181) * t367 + (t139 * t167 + t166 * t173) * t358) * t138 * t353) * pkin(2) / (t138 * t140 ^ 2 + 0.1e1) * t169 * t229;
t333 = mrSges(6,1) * t208;
t197 = sin(pkin(22));
t201 = cos(pkin(22));
t57 = cos(t58);
t44 = t197 * t57 - t201 * t374;
t45 = -t197 * t374 - t201 * t57;
t116 = atan2(t129, t128);
t111 = sin(t116);
t61 = pkin(1) * t111 + t197 * t373;
t112 = cos(t116);
t62 = pkin(1) * t112 + t201 * t373;
t31 = t44 * t62 + t45 * t61;
t332 = mrSges(11,2) * t31;
t235 = pkin(4) * (-t101 * t224 * t264 + t356 * t69);
t242 = pkin(4) * (t101 * t104 + t102 * t67);
t66 = 0.1e1 / t67 ^ 2;
t265 = pkin(9) * t103 * t219 / (t66 * t69 ^ 2 + 0.1e1);
t245 = pkin(5) * t66 * t69 * t265;
t255 = 0.1e1 / t67 * t265;
t285 = t104 * t381;
t21 = 0.2e1 * ((t60 * t285 + (t75 * t97 - t320) * pkin(5)) * t370 + t74 * t235) * t255 - 0.2e1 * ((-t74 * t97 + t241) * t370 + t74 * t242) * t245 + t63;
t325 = Ifges(11,3) * t21;
t207 = sin(pkin(18));
t311 = t170 * t207;
t221 = 0.1e1 / pkin(8);
t310 = t170 * t221;
t295 = Ifges(4,5) * t183 - Ifges(4,6) * t248;
t294 = t204 ^ 2 + t208 ^ 2;
t274 = t170 * t361;
t270 = mrSges(6,3) * t294;
t262 = 0.2e1 * t294;
t259 = t207 * t269;
t258 = t212 * t269;
t254 = mrSges(6,2) * t204 - t333;
t253 = mrSges(6,1) * t204 + mrSges(6,2) * t208;
t26 = -mrSges(6,3) * t204 * t37 - mrSges(6,2) * t36;
t27 = mrSges(6,1) * t36 - mrSges(6,3) * t323;
t250 = -t204 * t27 + t208 * t26;
t249 = mrSges(6,3) * t262;
t246 = (mrSges(5,1) * t47 - mrSges(5,2) * t46) * pkin(5);
t244 = m(6) * t262 / 0.2e1;
t243 = (-mrSges(10,1) * t148 + mrSges(10,2) * t147) * pkin(2);
t85 = -t111 * t210 - t112 * t206;
t86 = -t111 * t206 + t112 * t210;
t32 = -t44 * t86 + t45 * t85;
t33 = t44 * t85 + t45 * t86;
t238 = Ifges(11,5) * t33 + Ifges(11,6) * t32;
t237 = 0.2e1 * t243;
t236 = Ifges(9,3) + m(10) * (t147 ^ 2 + t148 ^ 2) * t228;
t19 = t343 + (-Ifges(6,2) * t204 + t330) * t37;
t20 = Ifges(6,5) * t36 + (Ifges(6,1) * t208 - t331) * t37;
t234 = -Ifges(5,6) * t36 + Ifges(5,5) * t37 + t204 * t20 / 0.2e1 + t19 * t363;
t176 = pkin(1) * t184 - pkin(7);
t164 = t172 - t291;
t161 = 0.1e1 / t162 ^ 2;
t144 = pkin(1) * t185 * t164 - t176 * t232;
t141 = -pkin(1) * t306 - t164 * t176;
t134 = (t144 * t361 + t141 * t207 / 0.2e1) * t310;
t133 = (t141 * t361 - t207 * t144 / 0.2e1) * t310;
t130 = 0.1e1 / t133 ^ 2;
t121 = t176 * t276 + t230 * pkin(7) * t385 + (-t164 * t184 - t306) * pkin(1);
t119 = (t307 + (0.2e1 * t176 * pkin(7) - t164 + t276) * t185) * pkin(1);
t118 = atan2(t134, t133);
t115 = cos(t118);
t114 = sin(t118);
t94 = -pkin(2) * t99 - pkin(16);
t73 = (0.1e1 / t162 * t314 / 0.2e1 + t161 * t286 * t388) / (-t161 * t312 + 0.1e1) + t77;
t72 = (-t197 * t86 - t201 * t85) * pkin(4) + t187;
t64 = ((t121 * t274 + t144 * t258 + t119 * t311 / 0.2e1 + t141 * t259) / t133 - (t119 * t274 + t141 * t258 - t121 * t311 / 0.2e1 - t144 * t259) * t134 * t130) / (t130 * t134 ^ 2 + 0.1e1) * t221;
t30 = -t44 * t61 + t45 * t62;
t29 = 0.2e1 * ((t71 * t285 + (t107 * t97 - t391) * pkin(5)) * t370 + t264 * t235) * t255 - 0.2e1 * ((-t264 * t97 + t240) * t370 + t264 * t242) * t245;
t25 = t253 * t37;
t24 = t36 * pkin(10) - t37 * pkin(12) + t168;
t15 = -pkin(10) * t18 - t371;
t14 = pkin(12) * t18 + t372;
t9 = pkin(12) * t12 + t39;
t4 = t254 * t18;
t1 = t254 * t12;
t7 = [(mrSges(8,2) * t384 + Ifges(8,1) * t86) * t86 + (0.2e1 * t204 * t26 + 0.2e1 * t208 * t27 + t244 * t24) * t24 + Ifges(3,1) * t210 ^ 2 - 0.2e1 * (mrSges(3,1) * t206 - mrSges(9,1) * t99 + mrSges(3,2) * t210 + mrSges(9,2) * t100) * pkin(16) + m(5) * t168 ^ 2 + (mrSges(5,2) * t386 + Ifges(5,1) * t37 - 0.2e1 * Ifges(5,4) * t36 + t208 * t20 + (-t19 - t343) * t204) * t37 + (mrSges(5,1) * t386 + Ifges(5,2) * t36 + t334) * t36 + (-0.2e1 * Ifges(3,4) * t210 + Ifges(3,2) * t206) * t206 + Ifges(9,2) * t99 ^ 2 + (Ifges(9,1) * t100 + 0.2e1 * Ifges(9,4) * t99) * t100 + (mrSges(4,1) * t384 + Ifges(4,2) * t248) * t248 + m(7) * pkin(15) ^ 2 + 0.2e1 * pkin(15) * (-mrSges(7,1) * t115 + mrSges(7,2) * t114) + t115 * (Ifges(7,4) * t114 + Ifges(7,2) * t115) + t114 * (Ifges(7,1) * t114 + Ifges(7,4) * t115) + m(10) * t94 ^ 2 + 0.2e1 * t94 * (-mrSges(10,1) * t81 + mrSges(10,2) * t82) + (-mrSges(8,1) * t384 + 0.2e1 * Ifges(8,4) * t86 + Ifges(8,2) * t85) * t85 + t81 * (Ifges(10,4) * t82 + Ifges(10,2) * t81) + t82 * (Ifges(10,1) * t82 + Ifges(10,4) * t81) + m(11) * t72 ^ 2 + 0.2e1 * t72 * (-mrSges(11,1) * t32 + mrSges(11,2) * t33) + (m(8) + m(4)) * t187 ^ 2 + (m(9) + m(3)) * pkin(16) ^ 2 + t32 * (Ifges(11,4) * t33 + Ifges(11,2) * t32) + t33 * (Ifges(11,1) * t33 + Ifges(11,4) * t32) + (mrSges(4,2) * t384 + Ifges(4,1) * t183 - 0.2e1 * Ifges(4,4) * t248) * t183 + Ifges(2,3); -t39 * t346 + t8 * t25 + t336 * t382 - Ifges(3,6) * t206 + Ifges(3,5) * t210 + t250 * t9 + (-t30 * t33 + t31 * t32) * mrSges(11,3) + (-t38 * mrSges(5,3) + t2 * t364 + t3 * t363) * t37 + t238 * t21 + ((t111 * t85 - t112 * t86) * mrSges(8,3) + (-t205 * t183 + t209 * t248) * mrSges(4,3)) * pkin(1) + t234 * t12 + t64 * (Ifges(7,5) * t114 + Ifges(7,6) * t115) + t63 * (Ifges(8,5) * t86 + Ifges(8,6) * t85) + t295 + t389 + t393; m(11) * t31 ^ 2 + (t38 ^ 2 + t39 ^ 2) * m(5) + t64 ^ 2 * Ifges(7,3) + t63 ^ 2 * Ifges(8,3) + Ifges(3,3) + Ifges(4,3) + Ifges(10,3) + (0.2e1 * t1 + t383) * t8 + t9 ^ 2 * t244 + (t325 - 0.2e1 * t332) * t21 + t237 + (m(11) * t30 + 0.2e1 * mrSges(11,1) * t21) * t30 + (m(8) * (t111 ^ 2 + t112 ^ 2) + m(4) * (t205 ^ 2 + t209 ^ 2)) * t230 + ((mrSges(8,1) * t112 - mrSges(8,2) * t111) * t63 + t390) * t396 + (Ifges(5,3) * t12 + t249 * t9 + t340 + t341 - 0.2e1 * t347 + 0.2e1 * t348) * t12 + t236; t15 * t25 + t335 * t382 - t346 * t372 + t73 * t393 + t250 * t14 + (-mrSges(5,3) * t371 + t363 * t6 + t364 * t5) * t37 + t238 * t29 + t389 * t77 + t234 * t18 + t295; t8 * t4 + Ifges(4,3) + (t1 + t383) * t15 + t9 * t14 * t244 + (t38 * t47 + t39 * t46) * pkin(5) * m(5) + t390 * pkin(1) + (Ifges(10,3) + t243) * t73 + (mrSges(11,1) * t30 + t325 - t332) * t29 + (t243 + t236) * t77 + (t348 - t347 + t340 / 0.2e1 + t341 / 0.2e1 + t9 * t270) * t18 + (t342 + t339 / 0.2e1 + t338 / 0.2e1 + t14 * t270 + t246) * t12; t29 ^ 2 * Ifges(11,3) + Ifges(4,3) + (t46 ^ 2 + t47 ^ 2) * t224 * m(5) + (m(6) * t15 + 0.2e1 * t4) * t15 + t14 ^ 2 * t244 + t236 * t77 ^ 2 + (Ifges(10,3) * t73 + t237 * t77) * t73 + (t14 * t249 + 0.2e1 * t246 + t338 + t339 + t342) * t18; t24 * t333 + (-mrSges(6,2) * t24 - Ifges(6,6) * t37) * t204 + t334; -t253 * t9 + t336; -t14 * t253 + t335; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
