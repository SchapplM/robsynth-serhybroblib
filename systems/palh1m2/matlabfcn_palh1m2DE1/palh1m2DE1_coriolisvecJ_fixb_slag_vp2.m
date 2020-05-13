% Calculate vector of centrifugal and Coriolis load on the joints for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh1m2DE1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_coriolisvecJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_coriolisvecJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:38
% EndTime: 2020-05-01 20:58:05
% DurationCPUTime: 9.64s
% Computational Cost: add. (9991->436), mult. (19784->688), div. (0->0), fcn. (23751->22), ass. (0->235)
t198 = cos(qJ(3));
t302 = pkin(1) * qJD(2);
t252 = t198 * t302;
t186 = pkin(22) + pkin(21);
t181 = sin(t186);
t182 = cos(t186);
t196 = sin(qJ(2));
t192 = cos(pkin(20));
t275 = sin(pkin(20));
t316 = sin(pkin(18));
t318 = cos(pkin(18));
t157 = t192 * t316 - t275 * t318;
t160 = t192 * t318 + t275 * t316;
t195 = sin(qJ(3));
t221 = t195 * t157 + t160 * t198;
t105 = t157 * t198 - t160 * t195;
t199 = cos(qJ(2));
t86 = t105 * t199;
t353 = t196 * t221 - t86;
t375 = t105 * t196 + t199 * t221;
t386 = -t353 * t181 + t182 * t375;
t211 = t375 * t181 + t182 * t353;
t253 = t195 * t302;
t187 = qJD(3) + qJD(2);
t312 = pkin(5) * t187;
t219 = t253 + t312;
t387 = t211 * t219;
t18 = t386 * t252 + t387;
t389 = t219 * t386;
t91 = t221 * qJD(3);
t96 = t105 * qJD(3);
t380 = qJD(2) * t375 + t96 * t196 + t199 * t91;
t384 = qJD(2) * t353 + t91 * t196;
t24 = t380 * t182 + (t199 * t96 - t384) * t181;
t180 = pkin(1) * t195 + pkin(5);
t218 = pkin(1) * t211;
t25 = (-qJD(3) * t86 + t384) * t182 + t380 * t181;
t315 = pkin(1) * t198;
t250 = qJD(3) * t315;
t319 = qJD(3) * t195 * t218 - t25 * t180 - t24 * t315 + t250 * t386 + t18;
t193 = cos(pkin(19));
t276 = sin(pkin(19));
t158 = t198 * t193 - t195 * t276;
t358 = pkin(2) * t187;
t114 = t158 * t358;
t190 = qJ(3) + qJ(2);
t184 = cos(t190);
t270 = Ifges(11,6) * t187;
t183 = sin(t190);
t273 = Ifges(11,4) * t183;
t134 = t270 + (Ifges(11,2) * t184 + t273) * qJD(1);
t258 = qJD(1) * t184;
t172 = Ifges(11,4) * t258;
t259 = qJD(1) * t183;
t135 = Ifges(11,1) * t259 + Ifges(11,5) * t187 + t172;
t162 = -t195 * t196 + t198 * t199;
t152 = t162 * qJD(1);
t139 = Ifges(4,4) * t152;
t155 = t195 * t193 + t198 * t276;
t144 = t155 * qJD(3);
t145 = t158 * qJD(3);
t164 = t195 * t199 + t196 * t198;
t151 = t164 * qJD(1);
t228 = mrSges(10,1) * t196 + mrSges(10,2) * t199;
t165 = t228 * qJD(1);
t245 = mrSges(10,3) * qJD(1) * t196;
t167 = -qJD(2) * mrSges(10,2) - t245;
t256 = qJD(1) * t199;
t244 = mrSges(10,3) * t256;
t168 = qJD(2) * mrSges(10,1) - t244;
t107 = -t155 * t196 + t158 * t199;
t213 = t107 * qJD(1);
t235 = qJD(2) * t250;
t271 = Ifges(11,6) * t183;
t286 = t151 * Ifges(4,4);
t301 = mrSges(4,3) * t151;
t325 = -t183 / 0.2e1;
t328 = t151 / 0.2e1;
t104 = t155 * t199 + t158 * t196;
t346 = t104 * qJD(1);
t338 = t346 / 0.2e1;
t343 = qJD(1) ^ 2;
t363 = mrSges(11,3) * t258 + mrSges(4,3) * t152;
t366 = t155 * t358;
t369 = mrSges(11,1) + mrSges(4,1);
t370 = Ifges(9,4) * t346;
t378 = t184 / 0.2e1;
t372 = Ifges(11,5) * t378;
t49 = Ifges(9,2) * t213 + Ifges(9,6) * t187 + t370;
t80 = Ifges(9,4) * t213;
t50 = Ifges(9,1) * t346 + Ifges(9,5) * t187 + t80;
t55 = -qJD(2) * t104 - t144 * t199 - t145 * t196;
t53 = t55 * qJD(1);
t56 = qJD(2) * t107 - t144 * t196 + t145 * t199;
t54 = t56 * qJD(1);
t70 = -pkin(2) * t213 - qJD(1) * pkin(15);
t74 = t152 * Ifges(4,2) + t187 * Ifges(4,6) + t286;
t75 = t151 * Ifges(4,1) + t187 * Ifges(4,5) + t139;
t78 = t144 * t358;
t79 = t145 * t358;
t120 = t187 * t164;
t97 = t120 * qJD(1);
t119 = t187 * t162;
t98 = t119 * qJD(1);
t385 = t134 * t259 / 0.2e1 + t343 * (Ifges(11,1) * t184 - t273) * t325 + Ifges(9,6) * t53 + Ifges(9,5) * t54 + t74 * t328 - t151 * (Ifges(4,1) * t152 - t286) / 0.2e1 - Ifges(4,6) * t97 + Ifges(4,5) * t98 - t252 * t301 + (t271 / 0.2e1 + t372) * qJD(1) * t187 - (-Ifges(4,2) * t151 + t139 + t75) * t152 / 0.2e1 + t363 * t253 - (-Ifges(11,2) * t259 + t135 + t172) * t258 / 0.2e1 + t369 * t235 - (Ifges(9,1) * t213 - t370) * t346 / 0.2e1 + t49 * t338 - (Ifges(4,5) * t152 + Ifges(9,5) * t213 - Ifges(4,6) * t151 - Ifges(9,6) * t346) * t187 / 0.2e1 - (-Ifges(9,2) * t346 + t50 + t80) * t213 / 0.2e1 + (t346 * (-m(10) * t70 - t165) + t145 * t168 + t155 * qJD(2) * t245 + m(10) * (-t114 * t144 + t145 * t366 + t155 * t79 - t158 * t78) - t144 * t167) * pkin(2);
t17 = -t211 * t252 + t389;
t207 = -t198 * t25 + (-t195 * t386 + t198 * t211) * qJD(3);
t6 = pkin(1) * t207 + t180 * t24;
t377 = -t17 + t6;
t212 = qJD(2) * t218;
t16 = -t198 * t212 + t389;
t20 = t180 * t386 - t198 * t218;
t210 = t195 * t212;
t3 = qJD(3) * t210 - t219 * t25 + t235 * t386 - t24 * t252;
t376 = t16 * t319 + t20 * t3;
t200 = cos(pkin(17));
t317 = sin(pkin(17));
t163 = t200 * t316 - t317 * t318;
t209 = t200 * t318 + t316 * t317;
t115 = t163 * t199 - t196 * t209;
t332 = -t115 / 0.2e1;
t371 = Ifges(7,4) * t115;
t368 = mrSges(11,2) + mrSges(4,2);
t274 = sin(pkin(22));
t277 = cos(pkin(22));
t156 = -t274 * t316 - t277 * t318;
t159 = t274 * t318 - t277 * t316;
t365 = -mrSges(11,1) * t184 - mrSges(8,1) * t156 + mrSges(11,2) * t183 + mrSges(8,2) * t159;
t194 = sin(qJ(4));
t197 = cos(qJ(4));
t178 = pkin(1) * t196 - pkin(15);
t125 = -pkin(5) * t152 + t178 * qJD(1);
t76 = -t157 * t182 + t160 * t181;
t71 = t76 * qJD(1);
t222 = t181 * t157 + t160 * t182;
t72 = t222 * qJD(1);
t42 = pkin(9) * t72 - pkin(11) * t71 + t125;
t11 = -t18 * t194 + t197 * t42;
t12 = t18 * t197 + t194 * t42;
t361 = (-t11 * t197 - t12 * t194) * qJD(4);
t146 = t209 * t199;
t255 = qJD(2) * t196;
t99 = -qJD(2) * t146 - t163 * t255;
t360 = 0.3e1 / 0.2e1 * t99;
t100 = t115 * qJD(2);
t359 = 0.3e1 / 0.2e1 * t100;
t236 = t163 * t196 + t146;
t357 = Ifges(7,4) * t236;
t68 = (-cos(pkin(21)) * t156 - t159 * sin(pkin(21))) * pkin(4) + t178;
t354 = t68 * (mrSges(11,1) * t183 + mrSges(11,2) * t184);
t351 = t178 * (m(4) + m(8));
t348 = (-mrSges(4,1) * t152 + mrSges(4,2) * t151 + qJD(1) * t365) * t199;
t226 = -t11 * t194 + t12 * t197;
t4 = t24 * t312 + (t195 * t24 + t207) * t302;
t254 = qJD(2) * t199;
t251 = pkin(1) * t254;
t67 = pkin(5) * t97 + qJD(1) * t251;
t1 = qJD(4) * t11 + t194 * t67 + t197 * t4;
t2 = -qJD(4) * t12 - t194 * t4 + t197 * t67;
t347 = t1 * t197 - t194 * t2;
t345 = pkin(15) * (mrSges(9,1) * t346 + mrSges(9,2) * t213) - t178 * (mrSges(4,1) * t151 + mrSges(4,2) * t152) - qJD(1) * t354 - t183 * t270;
t341 = t55 / 0.2e1;
t340 = t56 / 0.2e1;
t336 = t213 / 0.2e1;
t9 = t210 - t387;
t333 = t16 * t9;
t331 = t236 / 0.2e1;
t330 = t119 / 0.2e1;
t329 = -t120 / 0.2e1;
t326 = t152 / 0.2e1;
t323 = t194 / 0.2e1;
t314 = pkin(2) * t158;
t313 = pkin(5) * t151;
t311 = m(11) * t68;
t310 = mrSges(6,3) * t71;
t69 = qJD(4) + t72;
t309 = Ifges(6,5) * t69;
t308 = Ifges(7,6) * t99;
t304 = t69 * Ifges(6,6);
t299 = Ifges(3,4) * t196;
t298 = Ifges(3,4) * t199;
t297 = Ifges(6,4) * t194;
t296 = Ifges(6,4) * t197;
t295 = Ifges(10,4) * t196;
t294 = Ifges(10,4) * t199;
t293 = Ifges(7,5) * t100;
t290 = t115 * Ifges(7,2);
t289 = t236 * Ifges(7,1);
t285 = t197 * t71;
t234 = -pkin(5) * t162 + t178;
t43 = pkin(9) * t222 - pkin(11) * t76 + t234;
t268 = qJD(4) * t43;
t267 = qJD(4) * t71;
t262 = t183 * t198;
t242 = mrSges(11,3) * t259;
t261 = t187 * t369 - t242 - t301;
t260 = t187 * t368 - t363;
t249 = t222 * t267;
t248 = -0.3e1 / 0.2e1 * Ifges(10,4) - 0.3e1 / 0.2e1 * Ifges(3,4);
t247 = -Ifges(10,5) / 0.2e1 - Ifges(3,5) / 0.2e1;
t246 = -Ifges(3,6) / 0.2e1 - Ifges(10,6) / 0.2e1;
t233 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t231 = -t16 * t25 + t3 * t386;
t230 = mrSges(6,1) * t197 - mrSges(6,2) * t194;
t223 = t156 * t199 + t159 * t196;
t109 = t156 * t196 - t159 * t199;
t220 = -t184 * t195 + t262;
t94 = t223 * qJD(2);
t95 = t109 * qJD(2);
t217 = (t156 * t95 + t159 * t94) * mrSges(8,3);
t215 = t220 * t187;
t46 = -mrSges(6,2) * t69 - t194 * t310;
t47 = mrSges(6,1) * t69 - mrSges(6,3) * t285;
t206 = qJD(4) * (-t194 * t46 - t197 * t47 + (-t194 ^ 2 - t197 ^ 2) * t310);
t150 = Ifges(3,5) * qJD(2) + (t199 * Ifges(3,1) - t299) * qJD(1);
t149 = Ifges(10,5) * qJD(2) + (t199 * Ifges(10,1) - t295) * qJD(1);
t148 = Ifges(3,6) * qJD(2) + (-Ifges(3,2) * t196 + t298) * qJD(1);
t147 = Ifges(10,6) * qJD(2) + (-Ifges(10,2) * t196 + t294) * qJD(1);
t128 = pkin(1) * t256 + t313;
t84 = -pkin(2) * t107 - pkin(15);
t81 = pkin(5) * t120 + t251;
t58 = Ifges(7,5) * qJD(2) + (t289 + t371) * qJD(1);
t57 = Ifges(7,6) * qJD(2) + qJD(1) * (t290 + t357);
t48 = mrSges(5,1) * t72 + mrSges(5,2) * t71;
t45 = (mrSges(6,1) * t194 + mrSges(6,2) * t197) * t71;
t44 = t230 * t267;
t41 = t309 + (t197 * Ifges(6,1) - t297) * t71;
t40 = t304 + (-Ifges(6,2) * t194 + t296) * t71;
t21 = t180 * t211 + t315 * t386;
t14 = t128 * t194 + t17 * t197;
t13 = t128 * t197 - t17 * t194;
t10 = -t253 * t386 + t389;
t8 = t10 * t197 + t194 * t313;
t7 = -t10 * t194 + t197 * t313;
t5 = [(-t53 * t228 - t55 * t165 + (-t53 * t84 - t55 * t70) * m(10)) * pkin(2) - (-mrSges(5,1) * t67 + mrSges(5,3) * t4 - t233) * t222 + (t196 * t78 - t199 * t79) * mrSges(10,3) + (-t47 * t268 + m(6) * (t1 * t43 - t11 * t268 + t12 * t81) + t81 * t46 - Ifges(6,5) * t249 + (t3 * mrSges(6,1) - t1 * mrSges(6,3) + (-t41 / 0.2e1 + t11 * mrSges(6,3) - t309 / 0.2e1 - t16 * mrSges(6,2) + (0.3e1 / 0.2e1 * t297 + (0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(6,1)) * t197) * t71) * qJD(4)) * t76) * t194 + (-Ifges(6,6) * t249 + m(6) * (t11 * t81 + t12 * t268 + t2 * t43) + t81 * t47 + t46 * t268 + (t3 * mrSges(6,2) - t2 * mrSges(6,3) + (-t304 / 0.2e1 + t16 * mrSges(6,1) - t12 * mrSges(6,3) - t40 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,4) * t285) * qJD(4)) * t76) * t197 + t74 * t329 + t75 * t330 + (t293 / 0.2e1 + t308 / 0.2e1 + (-t147 / 0.2e1 - t148 / 0.2e1 + t70 * mrSges(10,1) - t114 * mrSges(10,3) + t246 * qJD(2)) * t199 + (-t149 / 0.2e1 - t150 / 0.2e1 - t70 * mrSges(10,2) + t366 * mrSges(10,3) + t247 * qJD(2)) * t196 + (t217 + t348 + (-t119 * t195 + t120 * t198 + (t162 * t195 - t164 * t198) * qJD(3)) * mrSges(4,3) + (-qJD(3) * t220 + t215) * mrSges(11,3)) * pkin(1)) * qJD(2) + (-t120 * t326 - t97 * t162) * Ifges(4,2) + (t119 * t326 - t120 * t328 + t98 * t162 - t97 * t164) * Ifges(4,4) + t178 * (mrSges(4,1) * t97 + mrSges(4,2) * t98) + (t125 * t81 + t234 * t67) * m(5) + (t67 * mrSges(5,2) + t3 * mrSges(5,3)) * t76 + t100 * t58 / 0.2e1 + t99 * t57 / 0.2e1 + t81 * t48 - pkin(15) * (-mrSges(9,1) * t53 + mrSges(9,2) * t54) + (t135 * t378 + t134 * t325 + Ifges(4,5) * t330 + Ifges(4,6) * t329 + Ifges(9,5) * t340 + Ifges(9,6) * t341 + (t372 - t271 / 0.2e1) * t187) * t187 + t50 * t340 + t49 * t341 + (t178 * (mrSges(4,1) * t120 + mrSges(4,2) * t119) + t289 * t359 + t290 * t360 - pkin(15) * (-mrSges(9,1) * t55 + mrSges(9,2) * t56) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t99 + mrSges(7,2) * t100) + (0.2e1 * pkin(15) * mrSges(3,2) - t84 * mrSges(10,2) - t196 * t248) * t255 + (t115 * t359 + t236 * t360) * Ifges(7,4) + (t84 * mrSges(10,1) - 0.2e1 * pkin(15) * mrSges(3,1) + t248 * t199 + (0.3e1 / 0.2e1 * Ifges(10,2) + 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(10,1)) * t196 + (-mrSges(4,1) * t162 + mrSges(4,2) * t164 + 0.2e1 * t311 + 0.2e1 * t351 + t365) * pkin(1)) * t254 + (0.2e1 * t354 + (0.3e1 / 0.2e1 * Ifges(11,1) - 0.3e1 / 0.2e1 * Ifges(11,2)) * t183 * t184 + (0.3e1 / 0.2e1 * t184 ^ 2 - 0.3e1 / 0.2e1 * t183 ^ 2) * Ifges(11,4)) * t187) * qJD(1) + (t119 * t328 + t98 * t164) * Ifges(4,1) + (t53 * t107 + t336 * t55) * Ifges(9,2) + (t54 * t104 + t338 * t56) * Ifges(9,1) + (t53 * t104 + t54 * t107 + t336 * t56 + t338 * t55) * Ifges(9,4); (t58 * t332 + t57 * t331 + t293 + t308 - t70 * (mrSges(10,1) * t199 - mrSges(10,2) * t196) + (pkin(15) * (mrSges(3,1) * t199 - mrSges(3,2) * t196) - t236 * (Ifges(7,1) * t115 - t357) / 0.2e1 - pkin(14) * (mrSges(7,1) * t236 + mrSges(7,2) * t115) + (-Ifges(7,2) * t236 + t371) * t332 - (-t294 - t298 + (-Ifges(3,1) - Ifges(10,1)) * t196) * t199 / 0.2e1) * qJD(1) + (t114 * t199 - t196 * t366) * mrSges(10,3) + (Ifges(7,5) * t332 + Ifges(7,6) * t331 + t247 * t196 + (-mrSges(10,3) * t314 + t246) * t199) * qJD(2) + t345 + (t149 + t150 + (-t295 - t299 + (-Ifges(3,2) - Ifges(10,2)) * t199) * qJD(1)) * t196 / 0.2e1 + (t147 + t148) * t199 / 0.2e1) * qJD(1) + (t319 * t71 - t377 * t72) * mrSges(5,3) + t21 * t206 + (-t194 * t6 - t13) * t47 - t128 * t48 + (t197 * t6 - t14) * t46 + t78 * mrSges(10,2) + t79 * mrSges(10,1) + t20 * t44 + t319 * t45 + (t226 * t6 + (t361 + t347) * t21 - t11 * t13 - t12 * t14 + t376) * m(6) + (-t125 * t128 + t18 * t377 + t21 * t4 + t376) * m(5) + ((mrSges(11,3) * t215 + t217 - t348 + (-mrSges(11,3) * t262 + (-t109 * t156 - t159 * t223) * mrSges(8,3)) * qJD(2)) * qJD(1) + (-t195 * t98 + t198 * t97) * mrSges(4,3) + (-t311 - t351) * t343 * t199 + (t261 * t198 + (-qJD(2) * t368 - t260) * t195) * qJD(3) + 0.2e1 * m(8) * (t109 * t94 - t223 * t95) * t302) * pkin(1) + t385; -m(5) * (t10 * t18 + t333) + (t10 * t72 - t71 * t9) * mrSges(5,3) - m(6) * (t11 * t7 + t12 * t8 + t333) + (-t151 * t48 + t386 * t44 + (-mrSges(5,3) * t71 - t45) * t25 + m(6) * (t211 * t361 + t231) + (m(6) * t347 + t206) * t211 + (m(6) * t226 - t72 * mrSges(5,3) - t194 * t47 + t197 * t46) * t24 + (-t125 * t151 + t18 * t24 + t211 * t4 + t231) * m(5)) * pkin(5) + t345 * qJD(1) + (-t244 * t314 + ((-t242 - t261) * t198 + (-qJD(3) * t368 + t260) * t195) * pkin(1)) * qJD(2) + t366 * t167 - t114 * t168 - t8 * t46 - t7 * t47 - t9 * t45 + t385; -t11 * t46 + t12 * t47 + (t197 * t40 / 0.2e1 + t41 * t323 - t16 * t230 + ((-Ifges(6,2) * t197 - t297) * t323 - t197 * (-Ifges(6,1) * t194 - t296) / 0.2e1) * t71 + t226 * mrSges(6,3) + (-t69 / 0.2e1 + qJD(4)) * (-Ifges(6,5) * t194 - Ifges(6,6) * t197)) * t71 + t233;];
tauc = t5(:);
