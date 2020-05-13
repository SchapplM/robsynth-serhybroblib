% Calculate vector of cutting forces with Newton-Euler
% palh2m1OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = palh2m1OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:17:25
% EndTime: 2020-05-03 00:17:43
% DurationCPUTime: 8.55s
% Computational Cost: add. (4537->417), mult. (5801->550), div. (0->0), fcn. (3750->10), ass. (0->230)
t190 = m(5) + m(6);
t176 = m(4) + t190;
t173 = m(3) + t176;
t183 = sin(qJ(2));
t188 = cos(qJ(2));
t312 = pkin(2) * t176;
t187 = cos(qJ(3));
t180 = sin(qJ(5));
t172 = mrSges(6,2) * t180;
t185 = cos(qJ(5));
t294 = t185 * mrSges(6,1);
t152 = -t172 + t294;
t315 = pkin(4) * m(6);
t135 = mrSges(5,1) + t152 + t315;
t181 = sin(qJ(4));
t275 = t181 * t135;
t314 = pkin(6) * m(6);
t169 = -mrSges(5,2) + mrSges(6,3) + t314;
t186 = cos(qJ(4));
t277 = t169 * t186;
t112 = t277 - t275;
t111 = mrSges(4,2) - t112;
t182 = sin(qJ(3));
t353 = t111 * t182;
t192 = pkin(3) * m(5);
t195 = m(6) * pkin(3);
t278 = t169 * t181;
t280 = t135 * t186;
t326 = t280 + t278;
t89 = mrSges(4,1) + t192 + t195 + t326;
t355 = t89 * t187 - t353;
t42 = t355 + t312;
t39 = mrSges(3,1) + t42;
t77 = t182 * t89;
t97 = t111 * t187;
t349 = -t97 - t77;
t46 = mrSges(3,2) - t349;
t360 = -t183 * t46 + t39 * t188;
t357 = pkin(1) * t173 + t360;
t7 = mrSges(2,1) + t357;
t291 = t39 * t183;
t13 = t46 * t188 + t291;
t296 = t183 * t349;
t19 = t188 * t355 + t296;
t197 = qJD(4) ^ 2;
t198 = qJD(3) ^ 2;
t199 = qJD(2) ^ 2;
t200 = qJD(1) ^ 2;
t105 = t112 * t187;
t273 = t182 * t326;
t56 = t105 - t273;
t295 = t183 * t56;
t341 = t112 * t182;
t344 = t326 * t187 + t341;
t27 = t344 * t188 + t295;
t290 = t344 * t183;
t354 = t56 * t188 - t290;
t297 = t183 * t355;
t358 = t188 * t349 - t297;
t365 = -t13 * qJDD(2) + t358 * qJDD(3) + t354 * qJDD(4) - t19 * t198 - t27 * t197 - t360 * t199 - t7 * t200;
t242 = -0.2e1 * t278;
t109 = -0.2e1 * t280 + t242;
t76 = -0.2e1 * mrSges(4,1) - 0.2e1 * t192 - 0.2e1 * t195 + t109;
t356 = t76 * t187 + 0.2e1 * t353;
t184 = sin(qJ(1));
t189 = cos(qJ(1));
t352 = g(1) * t189 + t184 * g(2);
t123 = 0.2e1 * t275;
t241 = 0.2e1 * t277;
t94 = t241 - t123 - 0.2e1 * mrSges(4,2);
t82 = t94 * t187;
t41 = (t82 - 0.2e1 * t77) * t188;
t351 = (t41 - 0.2e1 * t297) * qJD(3);
t319 = -0.2e1 * t183;
t108 = t241 - 0.2e1 * t275;
t91 = t108 * t187;
t53 = (t91 - 0.2e1 * t273) * t188;
t151 = t180 * mrSges(6,1) + mrSges(6,2) * t185;
t255 = qJD(5) * t151;
t272 = t182 * t181;
t137 = -t187 * t186 + t272;
t271 = t182 * t186;
t138 = t181 * t187 + t271;
t84 = t137 * t188 + t183 * t138;
t224 = t84 * t255;
t67 = 0.2e1 * t224;
t350 = t67 + (t319 * t344 + t53) * qJD(4);
t178 = t187 ^ 2;
t244 = 0.2e1 * t178;
t174 = pkin(3) * t190;
t100 = t174 + t326;
t346 = t100 * t187 + t341;
t345 = t109 * t187 - 0.2e1 * t341;
t339 = qJD(1) * t152;
t338 = (t53 - 0.2e1 * t290) * qJD(4);
t147 = mrSges(5,3) + t151;
t337 = t147 * qJDD(1) + t352 * t190;
t223 = qJD(5) * t339;
t218 = -0.2e1 * t223;
t269 = t183 * t182;
t332 = 0.2e1 * t187;
t335 = (pkin(3) * t269 - pkin(1) / 0.2e1) * t332;
t261 = t190 * t186;
t276 = t176 * t187;
t334 = pkin(2) * t276 + t89;
t196 = qJD(5) ^ 2;
t327 = -t151 * qJDD(5) - t152 * t196;
t177 = t186 ^ 2;
t281 = t135 * t177;
t333 = (t278 + t174 / 0.2e1) * t186 + t281;
t270 = t182 * t187;
t307 = t182 * pkin(2);
t81 = t176 * t307 - t111;
t330 = 0.2e1 * t89 * t270 + t81;
t160 = t169 * t177;
t148 = 0.2e1 * t160;
t328 = t148 - t169;
t325 = qJD(1) * (-qJD(3) * t358 - qJD(4) * t354 - t224);
t145 = mrSges(4,3) + t147;
t323 = -t145 * qJDD(1) - t176 * t352;
t230 = 0.2e1 * t255;
t219 = qJD(4) * t230;
t285 = qJD(4) * t27;
t316 = 0.2e1 * qJD(3);
t317 = 2 * qJD(2);
t86 = t183 * t137 - t138 * t188;
t322 = t218 + (-qJD(3) * t19 - t285) * t317 - t285 * t316 - t86 * t219 + t365;
t120 = t186 * t275;
t321 = -0.2e1 * t120;
t126 = t182 * t135;
t320 = -0.2e1 * t126;
t318 = -0.2e1 * t186;
t313 = pkin(1) * t176;
t311 = pkin(2) * t190;
t309 = g(3) * t176;
t308 = g(3) * t190;
t306 = t183 * g(3);
t304 = t188 * g(3);
t303 = -0.2e1 * t224 - t338;
t226 = t86 * t255;
t302 = (t345 * t188 - 0.2e1 * t295) * qJD(4) - 0.2e1 * t226;
t288 = (t109 * t182 + 0.2e1 * t105) * qJD(4) + t137 * t230;
t287 = t345 * qJD(4) + t138 * t230;
t286 = qJD(1) * t13;
t284 = t100 * t182;
t136 = mrSges(3,3) + t145;
t134 = mrSges(2,2) + t136;
t282 = t134 * t189;
t274 = t181 * t190;
t268 = t184 * t134;
t267 = t184 * t151;
t266 = t184 * t152;
t265 = t188 * t182;
t264 = t188 * t200;
t263 = t189 * t151;
t262 = t189 * t152;
t231 = -0.2e1 * t255;
t133 = t186 * t231;
t260 = t108 * qJD(4) + t133;
t259 = t109 * qJD(4) + t181 * t230;
t258 = t120 - t160;
t257 = t152 * qJDD(5) - t151 * t196;
t175 = qJD(2) + qJD(3);
t171 = qJD(4) + t175;
t256 = qJD(1) * t171;
t250 = t84 * qJDD(1);
t249 = qJD(1) * qJD(2);
t248 = qJD(4) * qJD(3);
t243 = -0.2e1 * t281;
t240 = 0.4e1 * t178 * t183;
t239 = qJDD(2) + qJDD(3);
t238 = pkin(3) * t274;
t233 = t190 * t307;
t232 = pkin(1) * t269;
t225 = t177 * t126;
t222 = t183 * t249;
t220 = t257 + ((t319 * t355 + t41) * qJD(3) + t350) * qJD(1);
t170 = m(2) + t173;
t217 = 0.2e1 * t223;
t215 = g(1) * t184 - t189 * g(2);
t212 = -t188 * t187 + t269;
t130 = t212 * qJD(1);
t139 = t187 * t183 + t265;
t131 = t139 * qJD(1);
t74 = t186 * t130 + t181 * t131;
t208 = t136 * qJDD(1) + t217;
t207 = t200 * pkin(1) + t352;
t142 = t169 + t238;
t206 = qJDD(1) * pkin(1) + t215;
t153 = t181 * t233;
t205 = t135 - t153;
t204 = -t142 + t148;
t203 = (t174 + 0.2e1 * t278) * t186 + 0.2e1 * t281 - t135;
t201 = t13 * t199 - t19 * qJDD(3) - t358 * t198 - t27 * qJDD(4) - t354 * t197 + (t303 - t351) * qJD(2) + t303 * qJD(3) - t84 * t219 - t360 * qJDD(2) + t327 * t86;
t179 = t188 ^ 2;
t168 = qJDD(4) + t239;
t167 = t171 ^ 2;
t164 = m(1) + t170;
t162 = -pkin(3) + t232;
t156 = t183 * pkin(1) - t182 * pkin(3);
t155 = -t182 * pkin(1) + t183 * pkin(3);
t143 = -t188 * qJDD(1) + t222;
t119 = pkin(3) * t261 + t135;
t103 = t306 - t188 * t207 + (-t179 * t200 - t199) * pkin(2);
t93 = (t138 * pkin(2) + pkin(3) * t181) * t190 + t169;
t88 = t304 + t207 * t183 + (t183 * t264 + qJDD(2)) * pkin(2);
t75 = t181 * t130 - t186 * t131;
t73 = (pkin(2) * t187 + pkin(3)) * t261 + t205;
t71 = t76 * t182;
t70 = qJD(5) - t74;
t63 = t184 * t226;
t62 = t189 * t226;
t61 = t180 * t171 + t185 * t75;
t60 = t185 * t171 - t180 * t75;
t48 = t311 + t346;
t45 = (t91 - 0.2e1 * t284) * t188;
t38 = t86 * qJDD(1) + t84 * t256;
t37 = t86 * t256 + qJDD(5) - t250;
t34 = -t74 * pkin(4) - t75 * pkin(6);
t33 = t187 * t103 + t182 * t88 + (-t130 ^ 2 - t175 ^ 2) * pkin(3);
t32 = t70 * mrSges(6,1) - t61 * mrSges(6,3);
t31 = -t70 * mrSges(6,2) + t60 * mrSges(6,3);
t30 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t29 = -t182 * t103 + t187 * t88 + (-t131 * t130 + t239) * pkin(3);
t17 = t60 * qJD(5) + t180 * t168 + t185 * t38;
t16 = -t61 * qJD(5) + t185 * t168 - t180 * t38;
t6 = ((t82 + t71 - 0.2e1 * mrSges(3,2)) * t188 - 0.2e1 * t291) * t249;
t2 = -t167 * pkin(4) + t168 * pkin(6) + t181 * t29 + t186 * t33 + t74 * t34;
t1 = (-t74 * t171 - t38) * pkin(6) + (-t250 + (qJD(1) * t86 + t75) * t171) * pkin(4) + (t182 * (-t183 * qJDD(1) - t188 * t249) - t187 * t143 - (qJD(3) + t175) * t131) * pkin(3) + (-t143 - t222) * pkin(2) + t206;
t3 = [(t262 * t84 + t267) * t196 - t164 * g(1) + 0.2e1 * t184 * t325 - t62 * t316 + (t184 * t286 - t62) * t317 + (t263 * t84 - t266) * qJDD(5) + t268 * t200 + (-t7 * t184 - t282) * qJDD(1) + t322 * t189, -t134 * qJDD(1) + t302 * qJD(3) + ((t356 * t188 - 0.2e1 * t296) * qJD(3) + t302) * qJD(2) - t327 * t84 - t352 * t170 + 0.2e1 * (-qJD(4) * t151 * t86 - t339) * qJD(5) + t365, -t355 * t198 - t39 * t199 + t56 * qJDD(4) - t344 * t197 - t46 * qJDD(2) + t349 * qJDD(3) + t173 * t306 + (t356 * qJD(3) + t287) * qJD(2) + t287 * qJD(3) + t138 * t219 - t327 * t137 + (-t173 * t352 - t357 * t200 - t208) * t188, -t89 * t198 + t112 * qJDD(4) - t326 * t197 - t334 * t199 + ((t76 * t178 + (-t182 * t94 - t312) * t187 + t89) * t179 + (-pkin(1) * t276 + (t111 * t244 + t330) * t183) * t188 + t89 * t178 - t111 * t270 + t176 * t232 - t89) * t200 - t111 * qJDD(3) + t81 * qJDD(2) + t181 * t219 + t259 * qJD(3) + (t76 * qJD(3) + t259) * qJD(2) + t139 * t309 + t327 * t186 - (t218 + t323) * t212, t93 * qJDD(2) - t86 * t308 + t142 * qJDD(3) - t135 * t197 + (-qJD(3) * t119 - qJD(4) * t135) * t317 + (((t135 + t243 + (-t174 + t242) * t186) * t244 + (-pkin(2) * t261 + (0.4e1 * t120 + 0.2e1 * t142 - 0.4e1 * t160) * t182) * t187 + t153 + t203) * t179 + t203 * t178 + (t156 * t274 + (t321 + t328) * t182) * t187 - t186 * (-t162 * t190 + t326) + ((t238 / 0.2e1 + t314 / 0.2e1 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + t258) * t240 + t261 * t335 - t155 * t274 + ((0.4e1 * t225 + t320 + (0.4e1 * t169 * t271 + t311) * t181) * t187 + (t275 - t233 / 0.2e1) * t318 + t328) * t183) * t188) * t200 - t119 * t198 + t169 * qJDD(4) - t73 * t199 - 0.2e1 * t135 * t248 + t327 + (-t218 + t337) * t84, m(6) * (t180 * t1 + t185 * t2) + t16 * mrSges(6,3) - t37 * mrSges(6,2) + t60 * t30 - t70 * t32; (t266 * t84 - t263) * t196 - 0.2e1 * t189 * t325 + (t267 * t84 + t262) * qJDD(5) - t282 * t200 + (t7 * t189 - t268) * qJDD(1) + (-t189 * t286 - t63) * t317 - t63 * t316 - t164 * g(2) + t322 * t184, t7 * qJDD(1) + t6 - t134 * t200 + t215 * t170 + (t338 + t67 + t351) * qJD(1) + t257, t355 * qJDD(3) + t349 * t198 + t344 * qJDD(4) + t56 * t197 + t39 * qJDD(2) + t137 * t219 + t173 * t304 + ((t71 - 0.2e1 * t97) * qJD(3) + t288) * qJD(2) + t288 * qJD(3) + t327 * t138 + (t207 * t173 + t39 * t264 + t208) * t183 + (-t199 + (t179 - 0.1e1) * t200) * t46, t334 * qJDD(2) + t81 * t199 + qJD(4) * t133 + (t94 * qJD(3) + t260) * qJD(2) + t260 * qJD(3) - t111 * t198 + t326 * qJDD(4) + t112 * t197 - t212 * t309 + ((-t178 * t94 + t330) * t179 + t265 * t313 + t349 * t187 + ((t244 - 0.1e1) * t188 * t89 + ((-0.2e1 * t353 + t312) * t188 + t313) * t187) * t183) * t200 + t89 * qJDD(3) + t327 * t181 + (t217 - t323) * t139, (((0.2e1 * t120 + t142 - 0.2e1 * t160) * t244 + (pkin(2) * t274 + 0.4e1 * t333 * t182 + t320) * t187 + (-t123 + t233) * t186 + t204) * t179 + ((-t294 / 0.2e1 + t172 / 0.2e1 - t315 / 0.2e1 - mrSges(5,1) / 0.2e1 + t333) * t240 - t274 * t335 - t155 * t261 + (((t135 * t272 - t311 / 0.4e1) * t318 + (0.2e1 * t177 - 0.1e1) * t182 * t169) * t332 + t243 + t186 * t242 + t205) * t183) * t188 + (t321 + t204) * t178 + (-0.2e1 * t225 + (t156 * t190 + t182 * t242) * t186 + t126) * t187 - t162 * t274 + t169 + t258) * t200 + t171 * t231 + (qJD(3) * t142 + qJD(4) * t169) * t317 + t93 * t199 + t142 * t198 - t84 * t308 + 0.2e1 * t169 * t248 + t73 * qJDD(2) + t119 * qJDD(3) + t135 * qJDD(4) + t169 * t197 - (t217 + t337) * t86, m(6) * (t185 * t1 - t180 * t2) - t17 * mrSges(6,3) + t37 * mrSges(6,1) - t61 * t30 + t70 * t31; -t164 * g(3) + t201, -t170 * g(3) + t201, t357 * qJDD(1) - t136 * t200 + t215 * t173 + t220 + t6, (t42 * t188 + t296) * qJDD(1) + (t42 * t319 + t41) * t249 - t145 * t200 + t206 * t176 + t220, (t48 * t188 + (t105 - t284) * t183) * qJDD(1) - t147 * t200 + t206 * t190 + ((t48 * t319 + t45) * qJD(2) + (t346 * t319 + t45) * qJD(3) + t350) * qJD(1) + t257, m(6) * (-t168 * pkin(4) - t167 * pkin(6) + t181 * t33 - t186 * t29 + t75 * t34) + t17 * mrSges(6,2) - t16 * mrSges(6,1) + t61 * t32 - t60 * t31;];
f_new = t3;
