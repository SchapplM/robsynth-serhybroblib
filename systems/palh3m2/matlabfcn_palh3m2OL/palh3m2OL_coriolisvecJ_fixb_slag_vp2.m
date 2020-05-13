% Calculate vector of centrifugal and Coriolis load on the joints for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% tauc [10x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh3m2OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_coriolisvecJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:43
% EndTime: 2020-05-07 04:34:46
% DurationCPUTime: 9.84s
% Computational Cost: add. (7859->581), mult. (20550->862), div. (0->0), fcn. (16438->16), ass. (0->274)
t243 = sin(qJ(4));
t244 = sin(qJ(3));
t353 = pkin(1) * qJD(2);
t306 = t244 * t353;
t221 = t243 * t306;
t249 = cos(qJ(4));
t307 = qJD(2) + qJD(3);
t297 = t307 * pkin(4);
t250 = cos(qJ(3));
t303 = t250 * t353;
t263 = -t297 + t303;
t175 = -t249 * t263 + t221;
t237 = qJD(4) + t307;
t170 = -t237 * pkin(8) - t175;
t242 = sin(qJ(5));
t248 = cos(qJ(5));
t289 = mrSges(6,1) * t242 + mrSges(6,2) * t248;
t267 = t170 * t289;
t245 = sin(qJ(2));
t251 = cos(qJ(2));
t213 = t244 * t245 - t250 * t251;
t204 = t213 * qJD(1);
t276 = t244 * t251 + t250 * t245;
t206 = t276 * qJD(1);
t294 = t249 * t204 + t206 * t243;
t280 = t204 * t243 - t249 * t206;
t349 = Ifges(5,1) * t280;
t365 = -t248 / 0.2e1;
t368 = t242 / 0.2e1;
t114 = t237 * t248 - t242 * t280;
t138 = qJD(5) - t294;
t115 = t237 * t242 + t248 * t280;
t333 = t115 * Ifges(6,4);
t38 = t114 * Ifges(6,2) + t138 * Ifges(6,6) + t333;
t113 = Ifges(6,4) * t114;
t39 = t115 * Ifges(6,1) + t138 * Ifges(6,5) + t113;
t416 = t175 * mrSges(5,3) + t39 * t365 + t38 * t368 - t267 - t349 / 0.2e1 - Ifges(5,5) * t237 - Ifges(5,4) * t294;
t332 = t294 * Ifges(5,2);
t379 = -t138 / 0.2e1;
t380 = -t115 / 0.2e1;
t381 = -t114 / 0.2e1;
t207 = t243 * t263;
t176 = t249 * t306 + t207;
t171 = pkin(10) * t237 - t176;
t228 = -t251 * pkin(1) - pkin(12);
t217 = t228 * qJD(1);
t174 = -t204 * pkin(4) + t217;
t62 = -pkin(8) * t294 - pkin(10) * t280 + t174;
t40 = -t171 * t242 + t248 * t62;
t41 = t171 * t248 + t242 * t62;
t415 = 0.2e1 * Ifges(6,5) * t380 + 0.2e1 * Ifges(6,6) * t381 + 0.2e1 * Ifges(6,3) * t379 + Ifges(5,6) * t237 + Ifges(5,4) * t280 - t40 * mrSges(6,1) + t41 * mrSges(6,2) + t332 / 0.2e1;
t240 = sin(qJ(7));
t246 = cos(qJ(7));
t212 = -t240 * t245 + t246 * t251;
t203 = t212 * qJD(1);
t214 = t240 * t251 + t246 * t245;
t205 = t214 * qJD(1);
t239 = sin(pkin(15));
t354 = cos(pkin(15));
t362 = sin(qJ(8));
t363 = cos(qJ(8));
t211 = t239 * t362 + t354 * t363;
t261 = -t239 * t363 + t354 * t362;
t97 = -t203 * t261 - t205 * t211;
t413 = Ifges(9,4) * t97;
t316 = t244 * t249;
t278 = t243 * t250 + t316;
t310 = qJD(4) * t249;
t258 = (qJD(3) * t278 + t244 * t310) * pkin(1);
t124 = qJD(2) * t258 + qJD(4) * t207;
t168 = t307 * t213;
t155 = t168 * qJD(1);
t169 = t307 * t276;
t156 = t169 * qJD(1);
t57 = qJD(4) * t294 + t155 * t249 + t156 * t243;
t26 = qJD(5) * t114 + t248 * t57;
t27 = -qJD(5) * t115 - t242 * t57;
t8 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t411 = m(6) * t124 - t8;
t410 = t242 * t41;
t408 = (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t280;
t295 = -t211 * t203 + t205 * t261;
t238 = qJD(2) + qJD(7);
t166 = t238 * t212;
t153 = t166 * qJD(1);
t167 = t238 * t214;
t154 = t167 * qJD(1);
t197 = t211 * qJD(8);
t198 = t261 * qJD(8);
t31 = -t153 * t211 + t154 * t261 - t197 * t203 + t198 * t205;
t32 = t153 * t261 + t154 * t211 + t197 * t205 + t198 * t203;
t236 = qJD(8) + t238;
t92 = Ifges(9,4) * t295;
t45 = Ifges(9,1) * t97 + Ifges(9,5) * t236 + t92;
t361 = pkin(3) * t239;
t208 = t238 * t361 + t240 * t353;
t302 = t354 * pkin(3);
t305 = t246 * t353;
t209 = t238 * t302 + t305;
t271 = pkin(1) * (-t211 * t246 + t240 * t261);
t264 = qJD(7) * t271;
t68 = qJD(2) * t264 - t197 * t209 + t198 * t208;
t272 = pkin(1) * (t211 * t240 + t246 * t261);
t265 = qJD(7) * t272;
t69 = qJD(2) * t265 + t197 * t208 + t198 * t209;
t405 = t69 * mrSges(9,1) - t68 * mrSges(9,2) + Ifges(9,5) * t31 + Ifges(9,6) * t32 - (-Ifges(9,2) * t97 + t45 + t92) * t295 / 0.2e1;
t356 = t295 * mrSges(9,3);
t401 = (-t197 * t354 + t198 * t239) * pkin(3) - qJD(2) * t271;
t400 = (t197 * t239 + t198 * t354) * pkin(3) - qJD(2) * t272;
t399 = -t242 * t40 + t248 * t41;
t313 = qJD(1) * t245;
t234 = pkin(1) * t313;
t226 = qJD(2) * t234;
t125 = -pkin(4) * t156 + t226;
t58 = qJD(4) * t280 + t155 * t243 - t249 * t156;
t11 = pkin(8) * t58 - pkin(10) * t57 + t125;
t317 = t243 * t244;
t277 = -t249 * t250 + t317;
t266 = t277 * qJD(3);
t123 = t297 * t310 + (qJD(4) * t277 + t266) * t353;
t2 = qJD(5) * t40 + t11 * t242 + t123 * t248;
t323 = qJD(5) * t41;
t3 = t11 * t248 - t123 * t242 - t323;
t398 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t26 + Ifges(6,6) * t27;
t382 = m(4) + m(8);
t396 = -mrSges(4,1) * t204 - mrSges(8,1) * t203 - mrSges(4,2) * t206 + mrSges(8,2) * t205 + t217 * t382;
t85 = pkin(8) * t280 - pkin(10) * t294;
t284 = Ifges(6,5) * t248 - Ifges(6,6) * t242;
t268 = t138 * t284;
t344 = Ifges(6,4) * t242;
t288 = Ifges(6,1) * t248 - t344;
t269 = t115 * t288;
t343 = Ifges(6,4) * t248;
t286 = -Ifges(6,2) * t242 + t343;
t270 = t114 * t286;
t282 = t248 * t40 + t410;
t395 = t282 * mrSges(6,3) - t270 / 0.2e1 - t269 / 0.2e1 - t268 / 0.2e1 - t174 * mrSges(5,2) + t416;
t394 = -0.2e1 * pkin(12);
t393 = t26 / 0.2e1;
t392 = t27 / 0.2e1;
t44 = Ifges(9,2) * t295 + Ifges(9,6) * t236 + t413;
t390 = -t44 / 0.2e1;
t389 = t58 / 0.2e1;
t385 = t295 / 0.2e1;
t384 = -t97 / 0.2e1;
t383 = t97 / 0.2e1;
t377 = t203 / 0.2e1;
t375 = t204 / 0.2e1;
t374 = t205 / 0.2e1;
t373 = -t206 / 0.2e1;
t372 = t206 / 0.2e1;
t371 = -t236 / 0.2e1;
t241 = sin(qJ(6));
t370 = -t241 / 0.2e1;
t369 = -t242 / 0.2e1;
t367 = t245 / 0.2e1;
t247 = cos(qJ(6));
t366 = t247 / 0.2e1;
t364 = t248 / 0.2e1;
t360 = pkin(4) * t206;
t358 = t2 * t248;
t357 = t242 * t3;
t352 = mrSges(6,3) * t248;
t351 = mrSges(8,3) * t203;
t350 = mrSges(8,3) * t205;
t348 = Ifges(3,4) * t245;
t347 = Ifges(3,4) * t251;
t346 = Ifges(4,4) * t206;
t342 = Ifges(7,4) * t241;
t341 = Ifges(3,5) * t251;
t338 = Ifges(3,6) * t245;
t103 = t208 * t211 + t209 * t261;
t334 = t103 * t97;
t330 = t204 * mrSges(4,3);
t329 = t205 * Ifges(8,4);
t328 = t206 * mrSges(4,3);
t325 = t251 * Ifges(3,2);
t324 = mrSges(5,1) * t237 + mrSges(6,1) * t114 - mrSges(6,2) * t115 - mrSges(5,3) * t280;
t322 = t280 * t176;
t321 = (mrSges(4,1) * t307 + t328) * t244;
t318 = t241 * t247;
t233 = qJD(1) * t347;
t315 = t251 * (Ifges(3,1) * t313 + Ifges(3,5) * qJD(2) + t233);
t312 = qJD(2) * t245;
t311 = qJD(4) * t243;
t309 = qJD(5) * t242;
t308 = qJD(5) * t248;
t304 = Ifges(7,2) * t318;
t235 = pkin(1) * t312;
t144 = -pkin(4) * t169 + t235;
t231 = -pkin(1) * t250 + pkin(4);
t194 = -pkin(1) * t316 + t243 * t231;
t291 = -t2 * t242 - t248 * t3;
t290 = mrSges(6,1) * t248 - mrSges(6,2) * t242;
t287 = Ifges(6,1) * t242 + t343;
t285 = Ifges(6,2) * t248 + t344;
t283 = Ifges(6,5) * t242 + Ifges(6,6) * t248;
t279 = t249 * t213 + t243 * t276;
t165 = t213 * t243 - t249 * t276;
t184 = -t213 * pkin(4) + t228;
t193 = pkin(1) * t317 + t231 * t249;
t275 = m(6) * t282;
t274 = pkin(6) * (mrSges(7,1) * t241 + mrSges(7,2) * t247);
t232 = Ifges(7,4) * t247 * qJD(1);
t273 = (Ifges(7,6) * qJD(6) + (Ifges(7,2) * t247 + t342) * qJD(1)) * t370 + (Ifges(7,1) * qJD(1) * t241 + Ifges(7,5) * qJD(6) + t232) * t366;
t70 = -t360 + t85;
t262 = (Ifges(7,5) * t366 + Ifges(7,6) * t370) * qJD(6);
t135 = (-t203 * t239 + t205 * t354) * pkin(3);
t260 = -qJD(5) * t282 - t357;
t259 = t260 * mrSges(6,3);
t6 = t26 * Ifges(6,4) + t27 * Ifges(6,2) + t58 * Ifges(6,6);
t7 = t26 * Ifges(6,1) + t27 * Ifges(6,4) + t58 * Ifges(6,5);
t257 = -t123 * mrSges(5,2) + t2 * t352 + t285 * t392 + t287 * t393 - t38 * t309 / 0.2e1 + t6 * t364 + t7 * t368 + qJD(5) * t267 + t39 * t308 / 0.2e1 + t283 * t389 - Ifges(5,6) * t58 + Ifges(5,5) * t57 + (t290 + mrSges(5,1)) * t124 + (t268 + t270 + t269) * qJD(5) / 0.2e1;
t102 = t208 * t261 - t209 * t211;
t131 = t203 * Ifges(8,2) + t238 * Ifges(8,6) + t329;
t191 = Ifges(8,4) * t203;
t133 = t205 * Ifges(8,1) + t238 * Ifges(8,5) + t191;
t256 = -t97 * t390 + t131 * t374 - t205 * (Ifges(8,1) * t203 - t329) / 0.2e1 - t238 * (Ifges(8,5) * t203 - Ifges(8,6) * t205) / 0.2e1 - Ifges(8,6) * t154 + Ifges(8,5) * t153 + t305 * t351 + (Ifges(9,1) * t295 - t413) * t384 + (Ifges(9,5) * t295 - Ifges(9,6) * t97) * t371 + t102 * t356 - (-Ifges(8,2) * t205 + t133 + t191) * t203 / 0.2e1 + t405;
t10 = -mrSges(6,2) * t58 + mrSges(6,3) * t27;
t71 = -mrSges(6,2) * t138 + mrSges(6,3) * t114;
t72 = mrSges(6,1) * t138 - mrSges(6,3) * t115;
t9 = mrSges(6,1) * t58 - mrSges(6,3) * t26;
t255 = m(6) * (-t308 * t40 - t309 * t41 - t357 + t358) + t248 * t10 - t242 * t9 - t71 * t309 - t72 * t308;
t254 = -t174 * mrSges(5,1) - t176 * mrSges(5,3) + t415;
t132 = Ifges(4,2) * t204 + t307 * Ifges(4,6) - t346;
t192 = Ifges(4,4) * t204;
t134 = -Ifges(4,1) * t206 + t307 * Ifges(4,5) + t192;
t253 = (Ifges(4,1) * t204 + t346) * t372 + t132 * t373 + t306 * t328 + t257 + Ifges(4,5) * t155 + Ifges(4,6) * t156 - t307 * (Ifges(4,5) * t204 + Ifges(4,6) * t206) / 0.2e1 - (Ifges(4,2) * t206 + t134 + t192) * t204 / 0.2e1 + (t244 * mrSges(4,1) + t250 * mrSges(4,2)) * qJD(3) * t353 + t415 * t280 + (mrSges(6,3) * t410 + t284 * t379 + t286 * t381 + t288 * t380 + t40 * t352 + t408 + t416) * t294;
t219 = t246 * pkin(1) + t302;
t218 = pkin(1) * t240 + t361;
t200 = Ifges(3,6) * qJD(2) + (t325 + t348) * qJD(1);
t190 = -t249 * t303 + t221;
t189 = t278 * t353;
t188 = pkin(10) + t194;
t187 = -pkin(8) - t193;
t180 = mrSges(8,1) * t238 - t350;
t179 = -mrSges(4,2) * t307 + t330;
t178 = -mrSges(8,2) * t238 + t351;
t177 = t234 - t360;
t160 = -mrSges(4,1) * t206 + mrSges(4,2) * t204;
t159 = mrSges(8,1) * t205 + mrSges(8,2) * t203;
t158 = (-t211 * t239 - t261 * t354) * pkin(3);
t157 = (-t211 * t354 + t239 * t261) * pkin(3);
t152 = -t231 * t311 + t258;
t151 = t231 * t310 + (t244 * t311 + t266) * pkin(1);
t136 = (-t212 * t354 - t214 * t239) * pkin(3) + t228;
t122 = t234 + t135;
t121 = -t211 * t218 - t219 * t261;
t120 = -t211 * t219 + t218 * t261;
t118 = -mrSges(5,2) * t237 + mrSges(5,3) * t294;
t116 = t217 + (-t203 * t354 - t205 * t239) * pkin(3);
t108 = -t211 * t214 - t212 * t261;
t107 = -t211 * t212 + t214 * t261;
t88 = mrSges(9,1) * t236 - mrSges(9,3) * t97;
t87 = -mrSges(9,2) * t236 + t356;
t86 = t235 + (-t166 * t239 + t167 * t354) * pkin(3);
t84 = -mrSges(5,1) * t294 + mrSges(5,2) * t280;
t83 = mrSges(5,1) * t280 + mrSges(5,2) * t294;
t76 = t226 + (-t153 * t239 + t154 * t354) * pkin(3);
t74 = t197 * t218 + t198 * t219 + t265;
t73 = -t197 * t219 + t198 * t218 + t264;
t67 = qJD(4) * t165 + t168 * t243 - t249 * t169;
t66 = qJD(4) * t279 + t168 * t249 + t169 * t243;
t65 = t234 + t70;
t61 = t175 * t248 + t242 * t85;
t60 = -t175 * t242 + t248 * t85;
t54 = Ifges(6,3) * t58;
t51 = t190 * t248 + t242 * t70;
t50 = -t190 * t242 + t248 * t70;
t49 = -mrSges(9,1) * t295 + mrSges(9,2) * t97;
t48 = mrSges(9,1) * t97 + mrSges(9,2) * t295;
t36 = t166 * t261 + t167 * t211 + t197 * t214 + t198 * t212;
t35 = -t166 * t211 + t167 * t261 - t197 * t212 + t198 * t214;
t1 = [(t396 * t245 + (-t166 * t246 - t167 * t240 + (t212 * t246 + t214 * t240) * qJD(7)) * mrSges(8,3) + (t168 * t250 - t169 * t244 + (-t213 * t250 + t244 * t276) * qJD(3)) * mrSges(4,3)) * t353 + t238 * (Ifges(8,5) * t166 - Ifges(8,6) * t167) / 0.2e1 + (t116 * t86 + t136 * t76) * m(9) + (t349 / 0.2e1 - t395) * t66 + (t125 * t184 + t144 * t174) * m(5) + (t7 * t364 + t6 * t369 + Ifges(5,1) * t57 - Ifges(5,4) * t58 + t125 * mrSges(5,2) + t288 * t393 + t286 * t392 + t284 * t389 + (-mrSges(5,3) - t289) * t124 + t291 * mrSges(6,3) + (-mrSges(6,3) * t399 + t170 * t290 + t283 * t379 + t285 * t381 + t287 * t380 + t38 * t365 + t39 * t369) * qJD(5)) * t165 + (t262 + t273) * qJD(6) + (-t212 * t154 - t167 * t377) * Ifges(8,2) + (t212 * t153 - t214 * t154 + t166 * t377 - t167 * t374) * Ifges(8,4) + (-mrSges(4,1) * t156 + mrSges(8,1) * t154 + mrSges(4,2) * t155 + mrSges(8,2) * t153) * t228 + (-t102 * t35 - t103 * t36 + t107 * t68 - t108 * t69) * mrSges(9,3) + (m(6) * (t308 * t41 - t309 * t40 - t291) + t248 * t9 + t242 * t10 - t72 * t309 + t71 * t308) * (-pkin(8) * t279 - t165 * pkin(10) + t184) - (-t123 * mrSges(5,3) + t54 / 0.2e1 - Ifges(5,4) * t57 + t125 * mrSges(5,1) + (Ifges(5,2) + Ifges(6,3) / 0.2e1) * t58 + t398) * t279 + (-t155 * t276 + t168 * t373) * Ifges(4,1) + (t213 * t155 - t156 * t276 + t168 * t375 + t169 * t373) * Ifges(4,4) + ((0.2e1 * t274 - 0.3e1 / 0.2e1 * t304 + 0.3e1 / 0.2e1 * Ifges(7,1) * t318 + (0.3e1 / 0.2e1 * t247 ^ 2 - 0.3e1 / 0.2e1 * t241 ^ 2) * Ifges(7,4)) * qJD(6) + ((mrSges(3,2) * t394 + 0.3e1 / 0.2e1 * t347) * t251 + (mrSges(3,1) * t394 - 0.3e1 / 0.2e1 * t348 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t251 + (-mrSges(4,1) * t213 - mrSges(8,1) * t212 - mrSges(4,2) * t276 + mrSges(8,2) * t214 + t228 * t382) * pkin(1)) * t245) * qJD(2)) * qJD(1) + (-mrSges(4,1) * t169 + mrSges(8,1) * t167 + mrSges(4,2) * t168 + mrSges(8,2) * t166) * t217 + t36 * t44 / 0.2e1 + t35 * t45 / 0.2e1 + t86 * t49 + t76 * (-mrSges(9,1) * t107 + mrSges(9,2) * t108) + t116 * (-mrSges(9,1) * t36 + mrSges(9,2) * t35) + t136 * (-mrSges(9,1) * t32 + mrSges(9,2) * t31) + (t242 * t71 + t248 * t72 + t275) * (pkin(8) * t67 - pkin(10) * t66 + t144) + t144 * t84 + t166 * t133 / 0.2e1 - t167 * t131 / 0.2e1 + t168 * t134 / 0.2e1 + t169 * t132 / 0.2e1 + t184 * (mrSges(5,1) * t58 + mrSges(5,2) * t57) + t307 * (Ifges(4,5) * t168 + Ifges(4,6) * t169) / 0.2e1 - t200 * t312 / 0.2e1 + qJD(2) * t315 / 0.2e1 + t236 * (Ifges(9,5) * t35 + Ifges(9,6) * t36) / 0.2e1 + (-t332 / 0.2e1 - t254) * t67 + qJD(2) ^ 2 * (-t338 + t341) / 0.2e1 + (t214 * t153 + t166 * t374) * Ifges(8,1) + (t213 * t156 + t169 * t375) * Ifges(4,2) + (t108 * t31 + t35 * t383) * Ifges(9,1) + (t107 * t32 + t36 * t385) * Ifges(9,2) + (t107 * t31 + t108 * t32 + t35 * t385 + t36 * t383) * Ifges(9,4); ((t178 * t246 - t180 * t240) * qJD(7) + (-t179 * t250 + t321) * qJD(3) + (-t153 * t246 - t154 * t240) * mrSges(8,3) + (t155 * t250 - t156 * t244) * mrSges(4,3) + (-t250 * t330 + t240 * t350 + (-mrSges(8,1) * t240 - mrSges(8,2) * t246) * qJD(7)) * qJD(2) - t396 * t313) * pkin(1) + (-t251 * t233 / 0.2e1 - t315 / 0.2e1 + t200 * t367 + (t325 * t367 + pkin(12) * (mrSges(3,1) * t245 + mrSges(3,2) * t251) - t245 * (Ifges(3,1) * t251 - t348) / 0.2e1) * qJD(1) + (t341 / 0.2e1 - t338 / 0.2e1) * qJD(2)) * qJD(1) + t324 * t152 + (-t159 - t160) * t217 + (-t120 * t31 + t121 * t32 - t334) * mrSges(9,3) + (-t193 * t57 - t194 * t58 - t322) * mrSges(5,3) + (-t151 * t72 - t65 * t71 + (-qJD(5) * t71 - t9) * t188 + (-t3 - t323) * mrSges(6,3)) * t242 + (t188 * t10 + t151 * t71 - t65 * t72 + (-t40 * mrSges(6,3) - t188 * t72) * qJD(5)) * t248 + m(5) * (t123 * t194 + t124 * t193 - t151 * t176 + t152 * t175) + m(9) * (t102 * t74 - t103 * t73 + t120 * t69 + t121 * t68) + (-m(9) * t122 - t48) * t116 + t253 + t256 + (-m(5) * t177 - t83) * t174 - t65 * t275 + t73 * t87 + t74 * t88 - t122 * t49 + t151 * t118 - t177 * t84 + t187 * t8 + (-t124 * t187 - t152 * t170 + (t260 + t358) * t188 + t399 * t151) * m(6); (t206 * t84 + (-t243 * t58 - t249 * t57) * mrSges(5,3) + (-t324 * t243 + (-t242 * t72 + t248 * t71 + t118) * t249 + m(6) * (t170 * t243 + t249 * t399)) * qJD(4) + (t123 * t243 + t124 * t249 + 0.2e1 * t174 * t372 + (-t175 * t243 - t176 * t249) * qJD(4)) * m(5)) * pkin(4) + t255 * (pkin(4) * t243 + pkin(10)) - t324 * t189 + t259 + t253 - m(5) * (t175 * t189 - t176 * t190) - mrSges(5,3) * t322 - m(6) * (-t170 * t189 + t40 * t50 + t41 * t51) - t51 * t71 - t50 * t72 - t174 * t83 - t190 * t118 - t217 * t160 + (-t321 + (t179 - t330) * t250) * t353 - t411 * (-pkin(4) * t249 - pkin(8)); t259 - t324 * t176 - m(6) * (-t170 * t176 + t40 * t60 + t41 * t61) + t254 * t280 + t255 * pkin(10) + (t408 + t395) * t294 + t257 - t61 * t71 - t60 * t72 - t175 * t118 + t411 * pkin(8); t54 - t170 * (mrSges(6,1) * t115 + mrSges(6,2) * t114) + (Ifges(6,1) * t114 - t333) * t380 + t115 * t38 / 0.2e1 + (Ifges(6,5) * t114 - Ifges(6,6) * t115) * t379 - t40 * t71 + t41 * t72 + (t114 * t40 + t115 * t41) * mrSges(6,3) + (-Ifges(6,2) * t115 + t113 + t39) * t381 + t398; (-t247 * t232 / 0.2e1 + (-t274 + (Ifges(7,1) * t247 - t342) * t370 + t304 / 0.2e1) * qJD(1) + t262 - t273) * qJD(1); t256 + t400 * t88 + t401 * t87 + (-t157 * t31 + t158 * t32 - t334) * mrSges(9,3) + ((-mrSges(8,2) * qJD(7) - t178) * t246 + (-mrSges(8,1) * qJD(7) + t180 + t350) * t240) * t353 - t116 * t48 - t135 * t49 - t217 * t159 + (t102 * t400 - t103 * t401 - t116 * t135 + t157 * t69 + t158 * t68) * m(9); -t103 * t88 + (-t87 + t356) * t102 + (-t116 * mrSges(9,2) + Ifges(9,1) * t384 + Ifges(9,5) * t371) * t295 - (t116 * mrSges(9,1) + mrSges(9,3) * t103 + Ifges(9,4) * t384 + Ifges(9,6) * t371 + t390) * t97 + t405; 0; 0;];
tauc = t1(:);
