% Calculate vector of centrifugal and Coriolis load on the joints for
% palh1m2DE2
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
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh1m2DE2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:29
% EndTime: 2020-05-02 20:59:04
% DurationCPUTime: 9.67s
% Computational Cost: add. (9991->436), mult. (19784->688), div. (0->0), fcn. (23751->22), ass. (0->235)
t198 = cos(qJ(3));
t300 = pkin(1) * qJD(2);
t251 = t198 * t300;
t184 = pkin(22) + pkin(21);
t179 = sin(t184);
t180 = cos(t184);
t195 = sin(qJ(2));
t191 = cos(pkin(20));
t196 = sin(pkin(18));
t273 = sin(pkin(20));
t315 = cos(pkin(18));
t156 = t191 * t196 - t273 * t315;
t160 = t191 * t315 + t196 * t273;
t194 = sin(qJ(3));
t221 = t194 * t156 + t160 * t198;
t105 = t156 * t198 - t160 * t194;
t199 = cos(qJ(2));
t86 = t105 * t199;
t350 = t195 * t221 - t86;
t372 = t195 * t105 + t199 * t221;
t383 = -t179 * t350 + t372 * t180;
t210 = t372 * t179 + t180 * t350;
t252 = t194 * t300;
t185 = qJD(2) + qJD(3);
t310 = pkin(5) * t185;
t219 = t252 + t310;
t384 = t210 * t219;
t18 = t383 * t251 + t384;
t386 = t219 * t383;
t91 = t221 * qJD(3);
t96 = t105 * qJD(3);
t377 = qJD(2) * t372 + t195 * t96 + t199 * t91;
t381 = qJD(2) * t350 + t195 * t91;
t24 = t377 * t180 + t179 * (t199 * t96 - t381);
t178 = pkin(1) * t194 + pkin(5);
t218 = pkin(1) * t210;
t25 = (-qJD(3) * t86 + t381) * t180 + t377 * t179;
t313 = pkin(1) * t198;
t250 = qJD(3) * t313;
t316 = qJD(3) * t194 * t218 - t25 * t178 - t24 * t313 + t250 * t383 + t18;
t192 = cos(pkin(19));
t274 = sin(pkin(19));
t157 = t198 * t192 - t194 * t274;
t355 = pkin(2) * t185;
t114 = t157 * t355;
t188 = qJ(2) + qJ(3);
t182 = cos(t188);
t269 = Ifges(11,6) * t185;
t181 = sin(t188);
t272 = Ifges(11,4) * t181;
t134 = t269 + (Ifges(11,2) * t182 + t272) * qJD(1);
t257 = qJD(1) * t182;
t172 = Ifges(11,4) * t257;
t258 = qJD(1) * t181;
t135 = Ifges(11,1) * t258 + Ifges(11,5) * t185 + t172;
t164 = -t194 * t195 + t198 * t199;
t152 = t164 * qJD(1);
t139 = Ifges(4,4) * t152;
t155 = t194 * t192 + t198 * t274;
t144 = t155 * qJD(3);
t145 = t157 * qJD(3);
t163 = t194 * t199 + t195 * t198;
t151 = t163 * qJD(1);
t227 = mrSges(10,1) * t195 + mrSges(10,2) * t199;
t165 = t227 * qJD(1);
t244 = mrSges(10,3) * qJD(1) * t195;
t167 = -qJD(2) * mrSges(10,2) - t244;
t255 = qJD(1) * t199;
t243 = mrSges(10,3) * t255;
t168 = qJD(2) * mrSges(10,1) - t243;
t107 = -t155 * t195 + t157 * t199;
t213 = t107 * qJD(1);
t234 = qJD(2) * t250;
t270 = Ifges(11,6) * t181;
t284 = t151 * Ifges(4,4);
t299 = mrSges(4,3) * t151;
t322 = -t181 / 0.2e1;
t325 = t151 / 0.2e1;
t104 = t155 * t199 + t157 * t195;
t343 = t104 * qJD(1);
t335 = t343 / 0.2e1;
t340 = qJD(1) ^ 2;
t360 = mrSges(11,3) * t257 + mrSges(4,3) * t152;
t363 = t155 * t355;
t366 = mrSges(11,1) + mrSges(4,1);
t367 = Ifges(9,4) * t343;
t375 = t182 / 0.2e1;
t369 = Ifges(11,5) * t375;
t49 = Ifges(9,2) * t213 + Ifges(9,6) * t185 + t367;
t80 = Ifges(9,4) * t213;
t50 = Ifges(9,1) * t343 + Ifges(9,5) * t185 + t80;
t55 = -qJD(2) * t104 - t144 * t199 - t145 * t195;
t53 = t55 * qJD(1);
t56 = qJD(2) * t107 - t144 * t195 + t145 * t199;
t54 = t56 * qJD(1);
t70 = -pkin(2) * t213 - qJD(1) * pkin(15);
t74 = t152 * Ifges(4,2) + t185 * Ifges(4,6) + t284;
t75 = t151 * Ifges(4,1) + t185 * Ifges(4,5) + t139;
t78 = t144 * t355;
t79 = t145 * t355;
t120 = t185 * t163;
t97 = t120 * qJD(1);
t119 = t185 * t164;
t98 = t119 * qJD(1);
t382 = t134 * t258 / 0.2e1 + t340 * (Ifges(11,1) * t182 - t272) * t322 + Ifges(9,6) * t53 + Ifges(9,5) * t54 + t74 * t325 - t151 * (Ifges(4,1) * t152 - t284) / 0.2e1 - Ifges(4,6) * t97 + Ifges(4,5) * t98 - t251 * t299 + (t270 / 0.2e1 + t369) * qJD(1) * t185 - (-Ifges(4,2) * t151 + t139 + t75) * t152 / 0.2e1 + t360 * t252 - (-Ifges(11,2) * t258 + t135 + t172) * t257 / 0.2e1 + t366 * t234 - (Ifges(9,1) * t213 - t367) * t343 / 0.2e1 + t49 * t335 - (Ifges(4,5) * t152 + Ifges(9,5) * t213 - Ifges(4,6) * t151 - Ifges(9,6) * t343) * t185 / 0.2e1 - (-Ifges(9,2) * t343 + t50 + t80) * t213 / 0.2e1 + (t343 * (-m(10) * t70 - t165) + t145 * t168 + t155 * qJD(2) * t244 + m(10) * (-t114 * t144 + t145 * t363 + t155 * t79 - t157 * t78) - t144 * t167) * pkin(2);
t17 = -t210 * t251 + t386;
t207 = -t198 * t25 + (-t194 * t383 + t198 * t210) * qJD(3);
t6 = pkin(1) * t207 + t178 * t24;
t374 = -t17 + t6;
t212 = qJD(2) * t218;
t16 = -t198 * t212 + t386;
t20 = t178 * t383 - t198 * t218;
t209 = t194 * t212;
t3 = qJD(3) * t209 - t219 * t25 + t234 * t383 - t24 * t251;
t373 = t316 * t16 + t20 * t3;
t200 = cos(pkin(17));
t314 = sin(pkin(17));
t162 = t196 * t200 - t314 * t315;
t211 = t196 * t314 + t200 * t315;
t115 = t162 * t199 - t195 * t211;
t329 = -t115 / 0.2e1;
t368 = Ifges(7,4) * t115;
t365 = mrSges(11,2) + mrSges(4,2);
t189 = sin(pkin(22));
t275 = cos(pkin(22));
t158 = t196 * t189 + t275 * t315;
t159 = t189 * t315 - t196 * t275;
t362 = -mrSges(11,1) * t182 + mrSges(8,1) * t158 + mrSges(11,2) * t181 + mrSges(8,2) * t159;
t193 = sin(qJ(4));
t197 = cos(qJ(4));
t176 = t195 * pkin(1) - pkin(15);
t125 = -pkin(5) * t152 + t176 * qJD(1);
t76 = -t156 * t180 + t160 * t179;
t71 = t76 * qJD(1);
t222 = t156 * t179 + t160 * t180;
t72 = t222 * qJD(1);
t42 = pkin(9) * t72 - pkin(11) * t71 + t125;
t11 = -t18 * t193 + t197 * t42;
t12 = t18 * t197 + t193 * t42;
t358 = (-t11 * t197 - t12 * t193) * qJD(4);
t146 = t211 * t199;
t254 = qJD(2) * t195;
t99 = -qJD(2) * t146 - t162 * t254;
t357 = 0.3e1 / 0.2e1 * t99;
t100 = t115 * qJD(2);
t356 = 0.3e1 / 0.2e1 * t100;
t235 = t162 * t195 + t146;
t354 = Ifges(7,4) * t235;
t68 = (cos(pkin(21)) * t158 - t159 * sin(pkin(21))) * pkin(4) + t176;
t351 = t68 * (mrSges(11,1) * t181 + mrSges(11,2) * t182);
t348 = t176 * (m(8) + m(4));
t345 = (-mrSges(4,1) * t152 + mrSges(4,2) * t151 + qJD(1) * t362) * t199;
t225 = -t11 * t193 + t12 * t197;
t4 = t24 * t310 + (t194 * t24 + t207) * t300;
t253 = qJD(2) * t199;
t249 = pkin(1) * t253;
t67 = pkin(5) * t97 + qJD(1) * t249;
t1 = qJD(4) * t11 + t193 * t67 + t197 * t4;
t2 = -qJD(4) * t12 - t193 * t4 + t197 * t67;
t344 = t1 * t197 - t193 * t2;
t342 = pkin(15) * (mrSges(9,1) * t343 + mrSges(9,2) * t213) - t176 * (mrSges(4,1) * t151 + mrSges(4,2) * t152) - qJD(1) * t351 - t181 * t269;
t338 = t55 / 0.2e1;
t337 = t56 / 0.2e1;
t333 = t213 / 0.2e1;
t9 = t209 - t384;
t330 = t16 * t9;
t328 = t235 / 0.2e1;
t327 = t119 / 0.2e1;
t326 = -t120 / 0.2e1;
t323 = t152 / 0.2e1;
t320 = t193 / 0.2e1;
t312 = pkin(2) * t157;
t311 = pkin(5) * t151;
t309 = m(11) * t68;
t308 = mrSges(6,3) * t71;
t69 = qJD(4) + t72;
t307 = Ifges(6,5) * t69;
t306 = Ifges(6,6) * t69;
t305 = Ifges(7,6) * t99;
t297 = Ifges(3,4) * t195;
t296 = Ifges(3,4) * t199;
t295 = Ifges(6,4) * t193;
t294 = Ifges(6,4) * t197;
t293 = Ifges(10,4) * t195;
t292 = Ifges(10,4) * t199;
t291 = Ifges(7,5) * t100;
t288 = t115 * Ifges(7,2);
t287 = t235 * Ifges(7,1);
t279 = t197 * t71;
t233 = -pkin(5) * t164 + t176;
t43 = pkin(9) * t222 - pkin(11) * t76 + t233;
t267 = qJD(4) * t43;
t266 = qJD(4) * t71;
t262 = t181 * t198;
t241 = mrSges(11,3) * t258;
t260 = t185 * t366 - t241 - t299;
t259 = t185 * t365 - t360;
t248 = t222 * t266;
t247 = 0.3e1 / 0.2e1 * Ifges(10,4) + 0.3e1 / 0.2e1 * Ifges(3,4);
t246 = -Ifges(10,5) / 0.2e1 - Ifges(3,5) / 0.2e1;
t245 = -Ifges(10,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t232 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t230 = -t16 * t25 + t3 * t383;
t229 = mrSges(6,1) * t197 - mrSges(6,2) * t193;
t108 = t158 * t199 - t159 * t195;
t109 = -t158 * t195 - t159 * t199;
t220 = -t182 * t194 + t262;
t94 = t108 * qJD(2);
t95 = t109 * qJD(2);
t217 = (-t158 * t95 - t159 * t94) * mrSges(8,3);
t215 = t220 * t185;
t46 = -mrSges(6,2) * t69 - t193 * t308;
t47 = mrSges(6,1) * t69 - mrSges(6,3) * t279;
t206 = qJD(4) * (-t193 * t46 - t197 * t47 + (-t193 ^ 2 - t197 ^ 2) * t308);
t150 = Ifges(3,5) * qJD(2) + (t199 * Ifges(3,1) - t297) * qJD(1);
t149 = Ifges(10,5) * qJD(2) + (t199 * Ifges(10,1) - t293) * qJD(1);
t148 = Ifges(3,6) * qJD(2) + (-t195 * Ifges(3,2) + t296) * qJD(1);
t147 = Ifges(10,6) * qJD(2) + (-t195 * Ifges(10,2) + t292) * qJD(1);
t128 = pkin(1) * t255 + t311;
t84 = -pkin(2) * t107 - pkin(15);
t81 = pkin(5) * t120 + t249;
t58 = Ifges(7,5) * qJD(2) + (t287 + t368) * qJD(1);
t57 = Ifges(7,6) * qJD(2) + qJD(1) * (t288 + t354);
t48 = mrSges(5,1) * t72 + mrSges(5,2) * t71;
t45 = (mrSges(6,1) * t193 + mrSges(6,2) * t197) * t71;
t44 = t229 * t266;
t41 = t307 + (t197 * Ifges(6,1) - t295) * t71;
t40 = t306 + (-t193 * Ifges(6,2) + t294) * t71;
t21 = t178 * t210 + t313 * t383;
t14 = t128 * t193 + t17 * t197;
t13 = t128 * t197 - t17 * t193;
t10 = -t252 * t383 + t386;
t8 = t10 * t197 + t193 * t311;
t7 = -t10 * t193 + t197 * t311;
t5 = [(t125 * t81 + t233 * t67) * m(5) - (-mrSges(5,1) * t67 + mrSges(5,3) * t4 - t232) * t222 + t50 * t337 + t49 * t338 + t74 * t326 + t75 * t327 + (t195 * t78 - t199 * t79) * mrSges(10,3) + (-Ifges(6,6) * t248 + m(6) * (t11 * t81 + t12 * t267 + t2 * t43) + t81 * t47 + t46 * t267 + (t3 * mrSges(6,2) - t2 * mrSges(6,3) + (-t12 * mrSges(6,3) + t16 * mrSges(6,1) - t306 / 0.2e1 - t40 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,4) * t279) * qJD(4)) * t76) * t197 + (m(6) * (t1 * t43 - t11 * t267 + t12 * t81) + t81 * t46 - Ifges(6,5) * t248 - t47 * t267 + (t3 * mrSges(6,1) - t1 * mrSges(6,3) + (t11 * mrSges(6,3) - t307 / 0.2e1 - t16 * mrSges(6,2) - t41 / 0.2e1 + (0.3e1 / 0.2e1 * t295 + (0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(6,1)) * t197) * t71) * qJD(4)) * t76) * t193 + (t291 / 0.2e1 + t305 / 0.2e1 + (-t147 / 0.2e1 - t148 / 0.2e1 + t70 * mrSges(10,1) - t114 * mrSges(10,3) + t245 * qJD(2)) * t199 + (-t149 / 0.2e1 - t150 / 0.2e1 - t70 * mrSges(10,2) + t363 * mrSges(10,3) + t246 * qJD(2)) * t195 + (t217 + t345 + (-t119 * t194 + t120 * t198 + (-t163 * t198 + t164 * t194) * qJD(3)) * mrSges(4,3) + (-qJD(3) * t220 + t215) * mrSges(11,3)) * pkin(1)) * qJD(2) + (-t120 * t323 - t97 * t164) * Ifges(4,2) + (t119 * t323 - t120 * t325 - t97 * t163 + t98 * t164) * Ifges(4,4) + t176 * (mrSges(4,1) * t97 + mrSges(4,2) * t98) + (-t53 * t227 - t55 * t165 + (-t53 * t84 - t55 * t70) * m(10)) * pkin(2) + (t135 * t375 + t134 * t322 + Ifges(4,5) * t327 + Ifges(4,6) * t326 + Ifges(9,5) * t337 + Ifges(9,6) * t338 + (t369 - t270 / 0.2e1) * t185) * t185 + (t176 * (mrSges(4,1) * t120 + mrSges(4,2) * t119) + t288 * t357 + t287 * t356 - pkin(15) * (-mrSges(9,1) * t55 + mrSges(9,2) * t56) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t99 + mrSges(7,2) * t100) + (0.2e1 * pkin(15) * mrSges(3,2) - t84 * mrSges(10,2) + t195 * t247) * t254 + (t115 * t356 + t235 * t357) * Ifges(7,4) + (t84 * mrSges(10,1) - 0.2e1 * pkin(15) * mrSges(3,1) - t247 * t199 + (-0.3e1 / 0.2e1 * Ifges(10,1) + 0.3e1 / 0.2e1 * Ifges(10,2) + 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(3,1)) * t195 + (-mrSges(4,1) * t164 + mrSges(4,2) * t163 + 0.2e1 * t309 + 0.2e1 * t348 + t362) * pkin(1)) * t253 + (0.2e1 * t351 + (-0.3e1 / 0.2e1 * Ifges(11,2) + 0.3e1 / 0.2e1 * Ifges(11,1)) * t181 * t182 + (-0.3e1 / 0.2e1 * t181 ^ 2 + 0.3e1 / 0.2e1 * t182 ^ 2) * Ifges(11,4)) * t185) * qJD(1) + t99 * t57 / 0.2e1 + t100 * t58 / 0.2e1 + t81 * t48 - pkin(15) * (-mrSges(9,1) * t53 + mrSges(9,2) * t54) + (t107 * t53 + t333 * t55) * Ifges(9,2) + (t54 * t104 + t335 * t56) * Ifges(9,1) + (t104 * t53 + t54 * t107 + t333 * t56 + t335 * t55) * Ifges(9,4) + (t119 * t325 + t98 * t163) * Ifges(4,1) + (t67 * mrSges(5,2) + t3 * mrSges(5,3)) * t76; (t316 * t71 - t374 * t72) * mrSges(5,3) + (t58 * t329 + t57 * t328 + t305 + t291 - t70 * (mrSges(10,1) * t199 - mrSges(10,2) * t195) + (pkin(15) * (mrSges(3,1) * t199 - mrSges(3,2) * t195) - t235 * (Ifges(7,1) * t115 - t354) / 0.2e1 - pkin(14) * (mrSges(7,1) * t235 + mrSges(7,2) * t115) + (-Ifges(7,2) * t235 + t368) * t329 - (-t292 - t296 + (-Ifges(3,1) - Ifges(10,1)) * t195) * t199 / 0.2e1) * qJD(1) + (t114 * t199 - t195 * t363) * mrSges(10,3) + (Ifges(7,5) * t329 + Ifges(7,6) * t328 + t246 * t195 + (-mrSges(10,3) * t312 + t245) * t199) * qJD(2) + t342 + (t149 + t150 + (-t293 - t297 + (-Ifges(3,2) - Ifges(10,2)) * t199) * qJD(1)) * t195 / 0.2e1 + (t147 + t148) * t199 / 0.2e1) * qJD(1) + (-t193 * t6 - t13) * t47 + (t197 * t6 - t14) * t46 + t21 * t206 - t128 * t48 + t78 * mrSges(10,2) + t79 * mrSges(10,1) + t20 * t44 + t316 * t45 + (-t11 * t13 - t12 * t14 + t225 * t6 + (t358 + t344) * t21 + t373) * m(6) + (-t125 * t128 + t374 * t18 + t21 * t4 + t373) * m(5) + ((mrSges(11,3) * t215 + t217 - t345 + (-mrSges(11,3) * t262 + (t108 * t159 + t109 * t158) * mrSges(8,3)) * qJD(2)) * qJD(1) + (-t194 * t98 + t198 * t97) * mrSges(4,3) + (-t309 - t348) * t340 * t199 + (t260 * t198 + (-qJD(2) * t365 - t259) * t194) * qJD(3) + 0.2e1 * m(8) * (t108 * t95 - t109 * t94) * t300) * pkin(1) + t382; (t10 * t72 - t71 * t9) * mrSges(5,3) + t342 * qJD(1) - m(6) * (t11 * t7 + t12 * t8 + t330) - m(5) * (t10 * t18 + t330) + (-t151 * t48 + t383 * t44 + (-t71 * mrSges(5,3) - t45) * t25 + m(6) * (t210 * t358 + t230) + (m(6) * t344 + t206) * t210 + (m(6) * t225 - mrSges(5,3) * t72 - t193 * t47 + t197 * t46) * t24 + (-t125 * t151 + t18 * t24 + t210 * t4 + t230) * m(5)) * pkin(5) + t363 * t167 - t114 * t168 - t9 * t45 - t8 * t46 - t7 * t47 + (-t243 * t312 + ((-t241 - t260) * t198 + (-qJD(3) * t365 + t259) * t194) * pkin(1)) * qJD(2) + t382; -t11 * t46 + t12 * t47 + (t197 * t40 / 0.2e1 + t41 * t320 - t16 * t229 + ((-Ifges(6,2) * t197 - t295) * t320 - t197 * (-Ifges(6,1) * t193 - t294) / 0.2e1) * t71 + t225 * mrSges(6,3) + (-t69 / 0.2e1 + qJD(4)) * (-Ifges(6,5) * t193 - Ifges(6,6) * t197)) * t71 + t232;];
tauc = t5(:);
