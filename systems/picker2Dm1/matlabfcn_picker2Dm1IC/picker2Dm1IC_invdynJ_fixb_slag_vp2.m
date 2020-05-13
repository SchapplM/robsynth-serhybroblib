% Calculate vector of inverse dynamics joint torques with ic for
% picker2Dm1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:55
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = picker2Dm1IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1IC_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:54:53
% EndTime: 2020-05-11 05:55:00
% DurationCPUTime: 6.94s
% Computational Cost: add. (4956->491), mult. (8356->609), div. (37->10), fcn. (4350->70), ass. (0->266)
t201 = sin(qJ(8));
t209 = cos(qJ(8));
t322 = (qJD(1) + qJD(8)) * (mrSges(9,1) * t201 + mrSges(9,2) * t209);
t191 = qJDD(1) + qJDD(2);
t172 = qJDD(3) + t191;
t214 = cos(qJ(2));
t308 = pkin(1) * qJD(1);
t167 = t214 * t308;
t193 = qJD(1) + qJD(2);
t101 = pkin(2) * t193 + t167;
t205 = sin(qJ(3));
t213 = cos(qJ(3));
t206 = sin(qJ(2));
t261 = t206 * t308;
t248 = t213 * t261;
t273 = qJD(3) * t205;
t189 = t214 * pkin(1);
t314 = pkin(1) * t206;
t260 = qJD(2) * t314;
t98 = -qJD(1) * t260 + qJDD(1) * t189;
t67 = pkin(2) * t191 + t98;
t274 = qJD(2) * t214;
t99 = (qJD(1) * t274 + qJDD(1) * t206) * pkin(1);
t25 = qJD(3) * t248 + t101 * t273 + t205 * t99 - t213 * t67;
t321 = t25 * mrSges(4,1) + Ifges(4,3) * t172;
t173 = qJD(6) + t193;
t302 = mrSges(7,2) * t173;
t305 = mrSges(7,1) * t173;
t203 = sin(qJ(6));
t211 = cos(qJ(6));
t227 = t203 * t206 - t211 * t214;
t94 = t227 * pkin(1);
t73 = qJD(1) * t94;
t228 = t203 * t214 + t206 * t211;
t76 = t228 * t308;
t320 = t73 * t302 - t76 * t305;
t100 = pkin(3) * t193 + t167;
t204 = sin(qJ(4));
t212 = cos(qJ(4));
t60 = t212 * t100 - t204 * t261;
t318 = pkin(1) * (qJD(2) + qJD(6));
t317 = pkin(2) * m(4);
t316 = pkin(3) * m(5);
t215 = cos(qJ(1));
t313 = pkin(1) * t215;
t285 = qJ(2) + qJ(4);
t186 = qJ(1) + t285;
t158 = qJ(10) + t186;
t127 = sin(t158);
t312 = g(1) * t127;
t197 = qJ(1) + qJ(2);
t187 = qJ(3) + t197;
t154 = sin(t187);
t311 = g(1) * t154;
t207 = sin(qJ(1));
t188 = t207 * pkin(1);
t171 = qJDD(4) + t191;
t62 = t100 * t204 + t212 * t261;
t66 = pkin(3) * t191 + t98;
t23 = -qJD(4) * t62 - t204 * t99 + t212 * t66;
t18 = pkin(4) * t171 + t23;
t198 = sin(qJ(10));
t199 = cos(qJ(10));
t22 = qJD(4) * t60 + t204 * t66 + t212 * t99;
t298 = t199 * t62;
t174 = qJD(4) + t193;
t57 = pkin(4) * t174 + t60;
t27 = t198 * t57 + t298;
t10 = qJD(10) * t27 - t18 * t199 + t198 * t22;
t130 = qJDD(10) + t171;
t310 = t10 * mrSges(11,1) + Ifges(11,3) * t130;
t19 = pkin(6) * t172 + t25;
t200 = sin(qJ(9));
t208 = cos(qJ(9));
t116 = t205 * t261;
t61 = -t101 * t213 + t116;
t24 = qJD(3) * t61 - t205 * t67 - t213 * t99;
t63 = t101 * t205 + t248;
t296 = t208 * t63;
t175 = qJD(3) + t193;
t58 = pkin(6) * t175 + t61;
t29 = t200 * t58 - t296;
t12 = qJD(9) * t29 - t19 * t208 + t200 * t24;
t142 = qJDD(9) + t172;
t309 = t12 * mrSges(10,1) + Ifges(10,3) * t142;
t307 = mrSges(4,1) * t175;
t306 = mrSges(5,1) * t174;
t145 = qJD(9) + t175;
t303 = mrSges(10,1) * t145;
t300 = t154 * mrSges(4,2);
t299 = t198 * t62;
t297 = t200 * t63;
t216 = 0.1e1 / pkin(5);
t196 = qJ(1) + qJ(8);
t194 = pkin(8) + qJ(5);
t284 = qJ(3) - qJ(6);
t244 = t194 + t284;
t220 = t244 - t196;
t245 = t194 - t284;
t221 = t245 - t196;
t295 = t216 / (pkin(6) * (cos(t221) - cos(t220)) + (-cos(qJ(9) - t221) + cos(qJ(9) + t220)) * pkin(2));
t218 = 0.1e1 / pkin(3);
t257 = -qJ(9) - t284;
t56 = 0.1e1 / ((-cos(qJ(4) - t284) + cos(qJ(4) + t284)) * pkin(6) + (cos(qJ(4) + t257) - cos(qJ(4) - t257)) * pkin(2));
t294 = t218 * t56;
t293 = t198 * t204;
t292 = t199 * t204;
t291 = t200 * t205;
t290 = t204 * t206;
t289 = t205 * t206;
t288 = t205 * t208;
t287 = t206 * t212;
t286 = t206 * t213;
t165 = -qJ(9) + t186;
t117 = t165 - t284;
t258 = qJ(9) + t197;
t164 = qJ(4) + t258;
t118 = t164 + t284;
t283 = sin(t118) - sin(t117);
t282 = cos(t118) - cos(t117);
t166 = qJ(3) + t258;
t136 = sin(t166);
t139 = cos(t166);
t281 = t136 * mrSges(10,1) + t139 * mrSges(10,2);
t153 = sin(t186);
t156 = cos(t186);
t280 = t153 * mrSges(5,1) + t156 * mrSges(5,2);
t279 = -cos(0.2e1 * t196) + cos(0.2e1 * t194);
t278 = sin(t165) - sin(t164);
t277 = cos(t165) - cos(t164);
t180 = sin(t197);
t182 = cos(t197);
t276 = -t180 * mrSges(3,1) - t182 * mrSges(3,2);
t150 = pkin(3) * t180;
t275 = t150 + t188;
t272 = qJD(3) * t213;
t271 = qJD(4) * t204;
t270 = qJD(4) * t212;
t269 = qJD(9) * t200;
t268 = qJD(9) * t208;
t267 = qJD(10) * t198;
t266 = qJD(10) * t199;
t265 = qJD(10) * t204;
t264 = qJD(1) * qJD(8);
t195 = 0.1e1 / t204;
t217 = 0.1e1 / pkin(4);
t263 = pkin(1) * t195 * t217;
t256 = -qJ(1) + t194;
t151 = pkin(2) * t180;
t255 = -pkin(6) * t154 + t151;
t119 = t153 * mrSges(5,2);
t254 = -t156 * mrSges(5,1) + t119;
t253 = -t139 * mrSges(10,1) + t136 * mrSges(10,2);
t185 = qJ(6) + t197;
t152 = sin(t185);
t155 = cos(t185);
t90 = t155 * mrSges(7,1) - mrSges(7,2) * t152;
t179 = sin(t196);
t181 = cos(t196);
t252 = t181 * mrSges(9,1) - mrSges(9,2) * t179;
t251 = t23 * mrSges(5,1) + Ifges(5,3) * t171 + t310;
t128 = cos(t158);
t109 = t128 * mrSges(11,1);
t250 = t127 * mrSges(11,2) - t109;
t162 = t189 + pkin(2);
t84 = pkin(1) * t289 - t162 * t213;
t247 = -qJ(4) + t256;
t246 = qJ(4) + t256;
t170 = qJDD(6) + t191;
t35 = qJD(6) * t73 - t203 * t98 - t211 * t99;
t36 = qJD(6) * t76 + t203 * t99 - t211 * t98;
t243 = t36 * mrSges(7,1) - t35 * mrSges(7,2) + Ifges(7,3) * t170;
t190 = qJDD(1) + qJDD(8);
t96 = (-qJDD(1) * t209 + t201 * t264) * pkin(1);
t97 = (-qJDD(1) * t201 - t209 * t264) * pkin(1);
t242 = t96 * mrSges(9,1) - t97 * mrSges(9,2) + Ifges(9,3) * t190;
t161 = t189 + pkin(3);
t83 = -pkin(1) * t290 + t212 * t161;
t157 = cos(t187);
t123 = t157 * mrSges(4,1);
t241 = -t123 - t253;
t240 = -pkin(2) * t182 - t313;
t239 = -pkin(3) * t182 - t313;
t238 = -t154 * mrSges(4,1) - t157 * mrSges(4,2);
t89 = -mrSges(7,1) * t152 - mrSges(7,2) * t155;
t237 = -mrSges(9,1) * t179 - mrSges(9,2) * t181;
t26 = -t199 * t57 + t299;
t71 = pkin(4) + t83;
t86 = pkin(1) * t287 + t161 * t204;
t40 = t198 * t86 - t199 * t71;
t236 = t198 * t71 + t199 * t86;
t28 = -t208 * t58 - t297;
t72 = pkin(6) + t84;
t87 = -pkin(1) * t286 - t162 * t205;
t44 = t200 * t87 - t208 * t72;
t235 = t200 * t72 + t208 * t87;
t234 = -qJ(8) + t247;
t233 = -qJ(8) + t246;
t232 = t127 * mrSges(11,1) + t128 * mrSges(11,2);
t231 = t153 * t182 - t156 * t180;
t202 = sin(qJ(7));
t210 = cos(qJ(7));
t230 = t153 * t202 + t156 * t210;
t226 = -t204 * t214 - t287;
t225 = t212 * t214 - t290;
t224 = t205 * t214 + t286;
t9 = qJD(10) * t26 - t18 * t198 - t199 * t22;
t223 = g(1) * t128 + g(2) * t127 - t9;
t222 = -g(1) * t281 + t309;
t219 = t98 * mrSges(3,1) - t99 * mrSges(3,2) + Ifges(3,3) * t191 + t243 + t251 + t309 + t321;
t184 = -qJ(9) + t194;
t183 = qJ(9) + t194;
t178 = sin(t285);
t177 = cos(t194);
t176 = sin(t194);
t160 = -pkin(2) * t213 + pkin(6);
t159 = pkin(3) * t212 + pkin(4);
t146 = t180 * mrSges(3,2);
t144 = qJD(10) + t174;
t141 = qJ(9) + t244;
t140 = -qJ(9) + t245;
t126 = pkin(4) * t153;
t125 = pkin(6) * t157;
t103 = -pkin(5) * t179 + t188;
t102 = -pkin(5) * t181 + t313;
t95 = t228 * pkin(1);
t85 = pkin(2) * t288 - t160 * t200;
t82 = -pkin(2) * t291 - t160 * t208;
t80 = -pkin(3) * t292 - t159 * t198;
t79 = pkin(3) * t293 - t159 * t199;
t78 = -t167 * t213 + t116;
t77 = t225 * t308;
t75 = t224 * t308;
t74 = t226 * t308;
t70 = t126 + t275;
t69 = pkin(4) * t156 - t239;
t55 = t228 * t318;
t54 = t227 * t318;
t53 = t162 * t273 + (qJD(2) * t224 + t206 * t272) * pkin(1);
t52 = -t162 * t272 + (t206 * t273 + (-t213 * t214 + t289) * qJD(2)) * pkin(1);
t51 = -t161 * t271 + (qJD(2) * t226 - t206 * t270) * pkin(1);
t50 = t161 * t270 + (qJD(2) * t225 - t206 * t271) * pkin(1);
t49 = t160 * t269 + (-t205 * t268 + (-t200 * t213 - t288) * qJD(3)) * pkin(2);
t48 = -t160 * t268 + (-t205 * t269 + (t208 * t213 - t291) * qJD(3)) * pkin(2);
t47 = t159 * t267 + (t199 * t265 + (t198 * t212 + t292) * qJD(4)) * pkin(3);
t46 = -t159 * t266 + (t198 * t265 + (-t199 * t212 + t293) * qJD(4)) * pkin(3);
t43 = -t200 * t75 - t208 * t78;
t42 = t200 * t78 - t208 * t75;
t39 = -t198 * t74 - t199 * t77;
t38 = t198 * t77 - t199 * t74;
t34 = -t208 * t61 - t297;
t33 = t200 * t61 - t296;
t32 = -t199 * t60 + t299;
t31 = t198 * t60 + t298;
t16 = qJD(9) * t235 + t200 * t52 - t208 * t53;
t15 = qJD(9) * t44 - t200 * t53 - t208 * t52;
t14 = qJD(10) * t236 + t198 * t50 - t199 * t51;
t13 = qJD(10) * t40 - t198 * t51 - t199 * t50;
t11 = qJD(9) * t28 - t19 * t200 - t208 * t24;
t6 = -g(1) * t89 - g(2) * t90 + t243 + t320;
t5 = -t29 * t303 - g(2) * t253 + (t145 * t28 - t11) * mrSges(10,2) + t222;
t4 = -g(2) * t109 + (-t144 * t27 + t312) * mrSges(11,1) + (t144 * t26 + t223) * mrSges(11,2) + t310;
t3 = -m(10) * (t28 * t33 - t29 * t34) - t33 * t303 + t241 * g(2) + (t145 * t34 - t11) * mrSges(10,2) + (-t175 * t63 + t311) * mrSges(4,1) + (g(1) * t157 + g(2) * t154 + t175 * t61 - t24) * mrSges(4,2) + ((t142 * t200 + t145 * t268) * mrSges(10,2) + (-t142 * t208 + t145 * t269) * mrSges(10,1) + (-g(2) * t157 - t11 * t200 - t12 * t208 + t268 * t29 + t269 * t28 + t311) * m(10)) * pkin(6) + t222 + t321;
t2 = -m(11) * (t26 * t31 - t27 * t32) + t62 * t306 - g(1) * t280 + (-t109 - t254) * g(2) + (t174 * t60 - t22) * mrSges(5,2) + (-t144 * t31 + t312) * mrSges(11,1) + (t144 * t32 + t223) * mrSges(11,2) + ((t130 * t198 + t144 * t266) * mrSges(11,2) + (-t130 * t199 + t144 * t267) * mrSges(11,1) + (-g(1) * t153 + g(2) * t156 - t10 * t199 - t198 * t9 + t26 * t267 + t266 * t27) * m(11)) * pkin(4) + t251;
t1 = (t204 * t22 + t212 * t23 + (-t204 * t60 + t212 * t62) * qJD(4)) * t316 + (-t205 * t24 - t213 * t25 + (t205 * t61 + t213 * t63) * qJD(3)) * t317 + (t174 * t77 - t22) * mrSges(5,2) + (-t89 - m(11) * (t150 + t126) - m(10) * t255 + (-t316 - t317) * t180 + t232 - t238 + t276 - t280 - t281) * g(1) + (t175 * t78 - t24) * mrSges(4,2) + (-t80 * t130 - t9 + (t39 - t46) * t144) * mrSges(11,2) + t320 + t219 + (mrSges(3,1) * t206 + mrSges(3,2) * t214) * t193 * t308 + ((-t171 * t204 - t174 * t270) * mrSges(5,2) + (t171 * t212 - t174 * t271) * mrSges(5,1)) * pkin(3) + ((t172 * t205 + t175 * t272) * mrSges(4,2) + (-t172 * t213 + t175 * t273) * mrSges(4,1)) * pkin(2) + (t82 * t142 + (-t42 + t49) * t145) * mrSges(10,1) + (t79 * t130 + (-t38 + t47) * t144) * mrSges(11,1) - m(5) * (t60 * t74 + t62 * t77) - m(4) * (t61 * t75 - t63 * t78) + m(11) * (t10 * t79 + t26 * t47 - t27 * t46 + t80 * t9) + m(10) * (t11 * t85 + t12 * t82 + t28 * t49 - t29 * t48) - m(11) * (t26 * t38 - t27 * t39) - m(10) * (t28 * t42 - t29 * t43) + (-m(10) * t125 + t300 - t119 - t146 - t90 + (pkin(4) * m(11) + mrSges(5,1)) * t156 + (mrSges(3,1) + (m(11) + m(5)) * pkin(3) + (m(4) + m(10)) * pkin(2)) * t182 + t241 + t250) * g(2) - t74 * t306 - t75 * t307 + (-t85 * t142 - t11 + (t43 - t48) * t145) * mrSges(10,2);
t7 = [-g(1) * (m(11) * t70 - t232) + t83 * t171 * mrSges(5,1) + t84 * t172 * mrSges(4,1) + t44 * t142 * mrSges(10,1) + m(9) * (-t201 * t97 - t209 * t96) * pkin(1) + m(3) * (t206 * t99 + t214 * t98) * pkin(1) + (-mrSges(9,1) * t209 + mrSges(9,2) * t201) * pkin(1) * t190 + t16 * t303 + t55 * t305 + t51 * t306 + t53 * t307 + ((-pkin(1) * t178 - pkin(3) * t204) * t195 * t1 + (pkin(3) * t206 + pkin(4) * t178) * t2 * t263 - pkin(1) * (pkin(3) * (cos(t247) - cos(t246)) + (cos(-qJ(2) + t234) - cos(qJ(2) + t233)) * pkin(5)) * t216 / (cos(t234) - cos(t233)) * t5) * t218 + (((t277 * t69 + t278 * t70) * t294 + (-(sin(t184) - sin(t183)) * t103 - (cos(t184) - cos(t183)) * t102) * t295) * t6 + ((-t282 * t69 - t283 * t70) * t294 + ((sin(t141) - sin(t140)) * t103 + (cos(t141) - cos(t140)) * t102) * t295) * t3) * pkin(2) + m(11) * (t10 * t40 - t13 * t27 + t14 * t26 - t236 * t9) + (-t13 * t144 + t130 * t236 - t9) * mrSges(11,2) + (t142 * t235 - t145 * t15 - t11) * mrSges(10,2) + m(10) * (-t11 * t235 + t12 * t44 - t15 * t29 + t16 * t28) - g(2) * (-m(3) * t313 + t146) + (g(2) * t182 + t189 * t191 - t193 * t260) * mrSges(3,1) + (-pkin(1) * t193 * t274 - t191 * t314) * mrSges(3,2) + (-t171 * t86 - t174 * t50 - t22) * mrSges(5,2) + (mrSges(7,1) * t94 + mrSges(7,2) * t95) * t170 + m(7) * (-t35 * t95 + t36 * t94 - t54 * t76 + t55 * t73) + (t130 * t40 + t14 * t144) * mrSges(11,1) + pkin(1) * qJD(8) * t322 + m(4) * (t24 * t87 + t25 * t84 - t52 * t63 + t53 * t61) + m(5) * (t22 * t86 + t23 * t83 + t50 * t62 + t51 * t60) - g(2) * (-m(9) * t313 + t252) - g(2) * (-m(7) * t313 + t90) + (-t172 * t87 - t175 * t52 - t24) * mrSges(4,2) + t219 + t201 / sin(-qJ(8) + t256) * (Ifges(6,3) * qJDD(5) - g(1) * (-mrSges(6,1) * t176 - mrSges(6,2) * t177) - g(2) * (mrSges(6,1) * t177 - mrSges(6,2) * t176)) - t206 * t4 * t263 - g(1) * (m(4) * (t151 + t188) + t238) - g(1) * (mrSges(2,1) * t207 + mrSges(2,2) * t215) - g(2) * (-mrSges(2,1) * t215 + mrSges(2,2) * t207) + Ifges(2,3) * qJDD(1) + t242 + (-t279 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) * t216 / t279 * (-g(1) * t237 - g(2) * t252 - t308 * t322 + t242) - g(2) * (-m(11) * t69 - t250) - g(2) * (m(10) * (t125 + t240) + t253) - g(2) * (m(5) * t239 + t254) - g(1) * (m(5) * t275 + t280) - g(1) * (m(10) * (t188 + t255) + t281) - g(2) * (m(4) * t240 + t123 - t300) - t54 * t302 - g(1) * (m(9) * t188 + t237) - g(1) * (m(7) * t188 + t89) - g(1) * (m(3) * t188 - t276); Ifges(8,3) * qJDD(7) - g(1) * (mrSges(8,1) * t210 - mrSges(8,2) * t202) - g(2) * (mrSges(8,1) * t202 + mrSges(8,2) * t210) + ((-t278 * t6 + t283 * t3) * t210 + (t277 * t6 - t282 * t3) * t202) * t56 * pkin(2) + (-(t1 + t5) * t230 + ((t2 * t230 + t231 * t4) * pkin(4) + (t2 - t4) * pkin(3) * (t180 * t202 + t182 * t210)) * t217) / t231;];
tau = t7(:);
