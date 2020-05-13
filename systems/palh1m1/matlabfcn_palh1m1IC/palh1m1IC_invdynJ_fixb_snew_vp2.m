% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% palh1m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% qJDD [13x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m1IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1IC_invdynJ_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:03:05
% EndTime: 2020-04-15 20:03:10
% DurationCPUTime: 4.78s
% Computational Cost: add. (23592->555), mult. (51578->734), div. (20->8), fcn. (39676->34), ass. (0->228)
t264 = sin(qJ(7));
t269 = sin(qJ(2));
t273 = cos(qJ(7));
t278 = cos(qJ(2));
t199 = (-t264 * t269 + t273 * t278) * qJD(1);
t298 = qJD(1) * qJD(2);
t296 = t278 * t298;
t214 = -qJDD(1) * t269 - t296;
t216 = qJDD(1) * t278 - t269 * t298;
t135 = -qJD(7) * t199 + t214 * t273 - t216 * t264;
t268 = sin(qJ(3));
t277 = cos(qJ(3));
t200 = (-t268 * t269 + t277 * t278) * qJD(1);
t136 = qJD(3) * t200 - t214 * t277 + t216 * t268;
t196 = (-t264 * t278 - t269 * t273) * qJD(1);
t138 = qJD(7) * t196 + t214 * t264 + t216 * t273;
t197 = (t268 * t278 + t269 * t277) * qJD(1);
t139 = -qJD(3) * t197 + t214 * t268 + t216 * t277;
t270 = sin(qJ(1));
t279 = cos(qJ(1));
t295 = g(1) * t270 - t279 * g(2);
t209 = -qJDD(1) * pkin(15) - t295;
t169 = t209 + (-t214 + t296) * pkin(1);
t254 = qJD(2) + qJD(7);
t180 = -mrSges(8,2) * t254 + mrSges(8,3) * t196;
t255 = qJD(2) + qJD(3);
t181 = mrSges(4,1) * t255 - mrSges(4,3) * t197;
t182 = mrSges(8,1) * t254 - mrSges(8,3) * t199;
t183 = -mrSges(4,2) * t255 + mrSges(4,3) * t200;
t267 = sin(qJ(4));
t276 = cos(qJ(4));
t162 = t197 * t276 + t200 * t267;
t241 = qJD(4) + t255;
t266 = sin(qJ(5));
t275 = cos(qJ(5));
t126 = t162 * t275 + t241 * t266;
t160 = -t197 * t267 + t200 * t276;
t158 = qJD(5) - t160;
t78 = -qJD(4) * t162 - t136 * t267 + t139 * t276;
t80 = qJD(4) * t160 + t136 * t276 + t139 * t267;
t94 = t169 + (t197 * t255 - t139) * pkin(5);
t44 = (-t160 * t241 - t80) * pkin(11) + (t162 * t241 - t78) * pkin(9) + t94;
t110 = -pkin(9) * t160 - pkin(11) * t162;
t234 = t241 ^ 2;
t252 = qJDD(2) + qJDD(3);
t236 = qJDD(4) + t252;
t280 = qJD(1) ^ 2;
t293 = -g(1) * t279 - g(2) * t270;
t211 = -pkin(15) * t280 + t293;
t187 = t269 * g(3) - t211 * t278;
t173 = (-t269 * t278 * t280 + qJDD(2)) * pkin(1) + t187;
t185 = -g(3) * t278 - t211 * t269;
t307 = t269 ^ 2;
t175 = (-qJD(2) ^ 2 - t280 * t307) * pkin(1) + t185;
t119 = t268 * t173 + t277 * t175;
t91 = (t197 * t200 + t252) * pkin(5) + t119;
t117 = -t173 * t277 + t268 * t175;
t98 = (-t200 ^ 2 - t255 ^ 2) * pkin(5) + t117;
t55 = t267 * t91 + t276 * t98;
t46 = -pkin(9) * t234 + pkin(11) * t236 + t110 * t160 + t55;
t34 = -t266 * t46 + t275 * t44;
t125 = -t162 * t266 + t241 * t275;
t52 = qJD(5) * t125 + t236 * t266 + t275 * t80;
t75 = qJDD(5) - t78;
t92 = -mrSges(6,1) * t125 + mrSges(6,2) * t126;
t95 = -mrSges(6,2) * t158 + mrSges(6,3) * t125;
t26 = m(6) * t34 + mrSges(6,1) * t75 - mrSges(6,3) * t52 - t126 * t92 + t158 * t95;
t35 = t266 * t44 + t275 * t46;
t51 = -qJD(5) * t126 + t236 * t275 - t266 * t80;
t96 = mrSges(6,1) * t158 - mrSges(6,3) * t126;
t27 = m(6) * t35 - mrSges(6,2) * t75 + mrSges(6,3) * t51 + t125 * t92 - t158 * t96;
t12 = t275 * t26 + t266 * t27;
t143 = -mrSges(5,2) * t241 + mrSges(5,3) * t160;
t145 = mrSges(5,1) * t241 - mrSges(5,3) * t162;
t289 = -m(5) * t94 + t78 * mrSges(5,1) - t80 * mrSges(5,2) + t160 * t143 - t162 * t145 - t12;
t259 = sin(pkin(19));
t260 = cos(pkin(19));
t261 = cos(qJ(10));
t299 = sin(qJ(10));
t203 = t259 * t261 - t260 * t299;
t204 = -t259 * t299 - t260 * t261;
t120 = t196 * t204 - t199 * t203;
t237 = qJD(10) + t254;
t114 = -mrSges(11,2) * t237 + mrSges(11,3) * t120;
t121 = t196 * t203 + t199 * t204;
t115 = mrSges(11,1) * t237 - mrSges(11,3) * t121;
t60 = -qJD(10) * t121 + t135 * t204 - t138 * t203;
t61 = qJD(10) * t120 + t135 * t203 + t138 * t204;
t63 = (-t135 * t260 - t138 * t259 + (-t196 * t259 + t199 * t260) * t254) * pkin(4) + t169;
t291 = -m(11) * t63 + t60 * mrSges(11,1) - t61 * mrSges(11,2) + t120 * t114 - t121 * t115;
t308 = t139 * mrSges(4,1) + t135 * mrSges(8,1) - t136 * mrSges(4,2) - t138 * mrSges(8,2) + t196 * t180 - t197 * t181 - t199 * t182 + t200 * t183 + t289 + t291 + (-m(4) - m(8)) * t169;
t102 = Ifges(5,4) * t162 + Ifges(5,2) * t160 + Ifges(5,6) * t241;
t104 = Ifges(5,1) * t162 + Ifges(5,4) * t160 + Ifges(5,5) * t241;
t54 = -t267 * t98 + t276 * t91;
t45 = -pkin(9) * t236 - pkin(11) * t234 + t110 * t162 - t54;
t67 = Ifges(6,5) * t126 + Ifges(6,6) * t125 + Ifges(6,3) * t158;
t69 = Ifges(6,1) * t126 + Ifges(6,4) * t125 + Ifges(6,5) * t158;
t16 = -mrSges(6,1) * t45 + mrSges(6,3) * t35 + Ifges(6,4) * t52 + Ifges(6,2) * t51 + Ifges(6,6) * t75 - t126 * t67 + t158 * t69;
t68 = Ifges(6,4) * t126 + Ifges(6,2) * t125 + Ifges(6,6) * t158;
t17 = mrSges(6,2) * t45 - mrSges(6,3) * t34 + Ifges(6,1) * t52 + Ifges(6,4) * t51 + Ifges(6,5) * t75 + t125 * t67 - t158 * t68;
t290 = -m(6) * t45 + t51 * mrSges(6,1) - mrSges(6,2) * t52 + t125 * t95 - t126 * t96;
t294 = -t26 * t266 + t275 * t27;
t3 = pkin(9) * t290 + pkin(11) * t294 + mrSges(5,1) * t54 - mrSges(5,2) * t55 + Ifges(5,5) * t80 + Ifges(5,6) * t78 + Ifges(5,3) * t236 + t162 * t102 - t160 * t104 + t275 * t16 + t266 * t17;
t306 = pkin(8) * t3;
t251 = qJDD(2) + qJDD(7);
t230 = qJDD(10) + t251;
t116 = t273 * t173 - t175 * t264;
t148 = (-t196 * t260 - t199 * t259) * pkin(4);
t249 = t254 ^ 2;
t81 = -t148 * t199 + (t249 * t259 + t251 * t260) * pkin(4) + t116;
t118 = t264 * t173 + t273 * t175;
t82 = t148 * t196 + (-t249 * t260 + t251 * t259) * pkin(4) + t118;
t48 = -t203 * t82 + t204 * t81;
t49 = t203 * t81 + t204 * t82;
t84 = Ifges(11,4) * t121 + Ifges(11,2) * t120 + Ifges(11,6) * t237;
t85 = Ifges(11,1) * t121 + Ifges(11,4) * t120 + Ifges(11,5) * t237;
t28 = mrSges(11,1) * t48 - mrSges(11,2) * t49 + Ifges(11,5) * t61 + Ifges(11,6) * t60 + Ifges(11,3) * t230 - t120 * t85 + t121 * t84;
t305 = pkin(10) * t28;
t109 = -mrSges(5,1) * t160 + mrSges(5,2) * t162;
t10 = m(5) * t55 - mrSges(5,2) * t236 + mrSges(5,3) * t78 + t109 * t160 - t145 * t241 + t294;
t20 = m(5) * t54 + mrSges(5,1) * t236 - mrSges(5,3) * t80 - t109 * t162 + t143 * t241 + t290;
t304 = t267 * t10 + t276 * t20;
t303 = pkin(4) * t259;
t302 = pkin(4) * t260;
t89 = -mrSges(11,1) * t120 + mrSges(11,2) * t121;
t40 = m(11) * t48 + mrSges(11,1) * t230 - mrSges(11,3) * t61 + t114 * t237 - t121 * t89;
t41 = m(11) * t49 - mrSges(11,2) * t230 + mrSges(11,3) * t60 - t115 * t237 + t120 * t89;
t301 = t203 * t41 + t204 * t40;
t300 = -t203 * t40 + t204 * t41;
t263 = sin(qJ(8));
t272 = cos(qJ(8));
t141 = t272 * t185 + t263 * t187;
t297 = qJD(1) * qJD(6);
t256 = -qJ(7) + pkin(19);
t253 = qJD(2) + qJD(8);
t250 = qJDD(2) + qJDD(8);
t140 = -t185 * t263 + t272 * t187;
t195 = (-t263 * t278 - t269 * t272) * qJD(1);
t198 = (-t263 * t269 + t272 * t278) * qJD(1);
t262 = sin(qJ(9));
t271 = cos(qJ(9));
t159 = -t195 * t271 + t198 * t262;
t161 = -t195 * t262 - t198 * t271;
t240 = qJD(9) + t253;
t101 = Ifges(10,4) * t161 + Ifges(10,2) * t159 + Ifges(10,6) * t240;
t103 = Ifges(10,1) * t161 + Ifges(10,4) * t159 + Ifges(10,5) * t240;
t235 = qJDD(9) + t250;
t105 = (t195 * t198 + t250) * pkin(2) + t140;
t111 = (-t195 ^ 2 - t253 ^ 2) * pkin(2) + t141;
t65 = -t105 * t271 + t111 * t262;
t66 = -t105 * t262 - t111 * t271;
t134 = -qJD(8) * t198 + t214 * t272 - t216 * t263;
t137 = qJD(8) * t195 + t214 * t263 + t216 * t272;
t77 = -qJD(9) * t161 - t134 * t271 + t137 * t262;
t79 = qJD(9) * t159 - t134 * t262 - t137 * t271;
t288 = mrSges(10,1) * t65 - mrSges(10,2) * t66 + Ifges(10,5) * t79 + Ifges(10,6) * t77 + Ifges(10,3) * t235 + t161 * t101 - t159 * t103;
t107 = (t198 * t253 - t134) * pkin(2) + t209;
t142 = -mrSges(10,2) * t240 + mrSges(10,3) * t159;
t144 = mrSges(10,1) * t240 - mrSges(10,3) * t161;
t287 = m(10) * t107 - t77 * mrSges(10,1) + t79 * mrSges(10,2) - t159 * t142 + t161 * t144;
t286 = mrSges(6,1) * t34 - mrSges(6,2) * t35 + Ifges(6,5) * t52 + Ifges(6,6) * t51 + Ifges(6,3) * t75 - t125 * t69 + t126 * t68;
t108 = -mrSges(10,1) * t159 + mrSges(10,2) * t161;
t152 = Ifges(9,4) * t198 + Ifges(9,2) * t195 + Ifges(9,6) * t253;
t155 = Ifges(9,1) * t198 + Ifges(9,4) * t195 + Ifges(9,5) * t253;
t284 = -mrSges(9,2) * t141 - t195 * t155 + t198 * t152 + Ifges(9,6) * t134 + Ifges(9,5) * t137 + mrSges(9,1) * t140 + Ifges(9,3) * t250 + t288 + pkin(2) * (-t262 * (m(10) * t66 - mrSges(10,2) * t235 + mrSges(10,3) * t77 + t108 * t159 - t144 * t240) - t271 * (m(10) * t65 + mrSges(10,1) * t235 - mrSges(10,3) * t79 - t108 * t161 + t142 * t240));
t153 = Ifges(8,4) * t199 + Ifges(8,2) * t196 + Ifges(8,6) * t254;
t156 = Ifges(8,1) * t199 + Ifges(8,4) * t196 + Ifges(8,5) * t254;
t283 = mrSges(8,1) * t116 - mrSges(8,2) * t118 + Ifges(8,5) * t138 + Ifges(8,6) * t135 + Ifges(8,3) * t251 + t199 * t153 - t196 * t156 + t300 * t303 + t301 * t302 + t28;
t154 = Ifges(4,4) * t197 + Ifges(4,2) * t200 + Ifges(4,6) * t255;
t157 = Ifges(4,1) * t197 + Ifges(4,4) * t200 + Ifges(4,5) * t255;
t281 = pkin(5) * t304 + mrSges(4,1) * t119 - mrSges(4,2) * t117 + Ifges(4,5) * t136 + Ifges(4,6) * t139 + Ifges(4,3) * t252 + t197 * t154 - t200 * t157 + t3;
t274 = cos(qJ(6));
t265 = sin(qJ(6));
t258 = qJ(8) + qJ(9);
t257 = qJ(3) + pkin(17);
t246 = cos(t258);
t245 = sin(t258);
t244 = pkin(20) + qJ(7) + qJ(2);
t243 = qJ(3) + qJ(4) + pkin(18);
t242 = -qJ(10) + t256;
t239 = cos(t257);
t238 = sin(t257);
t229 = cos(t244);
t228 = cos(t243);
t227 = sin(t244);
t226 = sin(t243);
t225 = cos(t242);
t224 = sin(t242);
t218 = -pkin(2) * t272 + pkin(12) * t246;
t217 = pkin(2) * t263 - pkin(12) * t245;
t215 = qJDD(1) * t274 - t265 * t297;
t213 = qJDD(1) * t265 + t274 * t297;
t212 = pkin(14) * t280 + t293;
t210 = qJDD(1) * pkin(14) - t295;
t208 = pkin(1) * t278 + pkin(3) * t229;
t207 = pkin(1) * t269 + pkin(3) * t227;
t206 = pkin(5) * t277 + pkin(10) * t228;
t205 = pkin(5) * t268 + pkin(10) * t226;
t202 = t225 * pkin(8) - pkin(4) * cos(t256);
t201 = t224 * pkin(8) - pkin(4) * sin(t256);
t194 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t278 - Ifges(3,4) * t269) * qJD(1);
t193 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t265 + Ifges(7,4) * t274) * qJD(1);
t192 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t278 - Ifges(3,2) * t269) * qJD(1);
t191 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t265 + Ifges(7,2) * t274) * qJD(1);
t186 = -g(3) * t265 + t212 * t274;
t184 = -g(3) * t274 - t212 * t265;
t174 = 0.1e1 / (-t224 * t226 + t225 * t228) / pkin(10) / pkin(8);
t168 = -mrSges(8,1) * t196 + mrSges(8,2) * t199;
t167 = -mrSges(4,1) * t200 + mrSges(4,2) * t197;
t151 = Ifges(4,5) * t197 + Ifges(4,6) * t200 + Ifges(4,3) * t255;
t150 = Ifges(8,5) * t199 + Ifges(8,6) * t196 + Ifges(8,3) * t254;
t149 = Ifges(9,5) * t198 + Ifges(9,6) * t195 + Ifges(9,3) * t253;
t100 = Ifges(5,5) * t162 + Ifges(5,6) * t160 + Ifges(5,3) * t241;
t99 = Ifges(10,5) * t161 + Ifges(10,6) * t159 + Ifges(10,3) * t240;
t83 = Ifges(11,5) * t121 + Ifges(11,6) * t120 + Ifges(11,3) * t237;
t43 = mrSges(10,2) * t107 - mrSges(10,3) * t65 + Ifges(10,1) * t79 + Ifges(10,4) * t77 + Ifges(10,5) * t235 - t101 * t240 + t159 * t99;
t42 = -mrSges(10,1) * t107 + mrSges(10,3) * t66 + Ifges(10,4) * t79 + Ifges(10,2) * t77 + Ifges(10,6) * t235 + t103 * t240 - t161 * t99;
t30 = mrSges(11,2) * t63 - mrSges(11,3) * t48 + Ifges(11,1) * t61 + Ifges(11,4) * t60 + Ifges(11,5) * t230 + t120 * t83 - t237 * t84;
t29 = -mrSges(11,1) * t63 + mrSges(11,3) * t49 + Ifges(11,4) * t61 + Ifges(11,2) * t60 + Ifges(11,6) * t230 - t121 * t83 + t237 * t85;
t18 = mrSges(9,2) * t209 - mrSges(9,3) * t140 + Ifges(9,1) * t137 + Ifges(9,4) * t134 + Ifges(9,5) * t250 + t149 * t195 - t152 * t253 + t262 * t42 - t271 * t43;
t13 = -pkin(2) * t287 - mrSges(9,1) * t209 + mrSges(9,3) * t141 + Ifges(9,4) * t137 + Ifges(9,2) * t134 + Ifges(9,6) * t250 - t198 * t149 + t253 * t155 - t262 * t43 - t271 * t42;
t8 = mrSges(8,2) * t169 - mrSges(8,3) * t116 + Ifges(8,1) * t138 + Ifges(8,4) * t135 + Ifges(8,5) * t251 + t150 * t196 - t153 * t254 - t203 * t29 + t204 * t30 + t291 * t303;
t7 = -mrSges(8,1) * t169 + mrSges(8,3) * t118 + Ifges(8,4) * t138 + Ifges(8,2) * t135 + Ifges(8,6) * t251 - t150 * t199 + t156 * t254 + t203 * t30 + t204 * t29 + t291 * t302;
t5 = -pkin(9) * t12 - mrSges(5,1) * t94 + mrSges(5,3) * t55 + Ifges(5,4) * t80 + Ifges(5,2) * t78 + Ifges(5,6) * t236 - t100 * t162 + t104 * t241 - t286;
t4 = -pkin(11) * t12 + mrSges(5,2) * t94 - mrSges(5,3) * t54 + Ifges(5,1) * t80 + Ifges(5,4) * t78 + Ifges(5,5) * t236 + t100 * t160 - t102 * t241 - t16 * t266 + t17 * t275;
t2 = mrSges(4,2) * t169 - mrSges(4,3) * t119 + Ifges(4,1) * t136 + Ifges(4,4) * t139 + Ifges(4,5) * t252 + t151 * t200 - t154 * t255 - t267 * t5 + t276 * t4;
t1 = pkin(5) * t289 - mrSges(4,1) * t169 + mrSges(4,3) * t117 + Ifges(4,4) * t136 + Ifges(4,2) * t139 + Ifges(4,6) * t252 - t197 * t151 + t255 * t157 + t267 * t4 + t276 * t5;
t6 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t295 - mrSges(2,2) * t293 + t278 * (mrSges(3,2) * t209 - mrSges(3,3) * t187 + Ifges(3,1) * t216 + Ifges(3,4) * t214 + Ifges(3,5) * qJDD(2) - qJD(2) * t192 + t1 * t277 - t13 * t263 + t18 * t272 + t2 * t268 - t264 * t7 + t273 * t8) - t269 * (pkin(1) * t308 - mrSges(3,1) * t209 + mrSges(3,3) * t185 + Ifges(3,4) * t216 + Ifges(3,2) * t214 + Ifges(3,6) * qJDD(2) + qJD(2) * t194 + t268 * t1 + t272 * t13 + t263 * t18 - t277 * t2 + t264 * t8 + t273 * t7) + pkin(15) * (-t198 * (mrSges(9,1) * t253 - mrSges(9,3) * t198) + t195 * (-mrSges(9,2) * t253 + mrSges(9,3) * t195) - t216 * mrSges(3,2) + t214 * mrSges(3,1) - t137 * mrSges(9,2) + t134 * mrSges(9,1) - t287 + (-m(9) - m(3)) * t209 + t308) + t265 * (mrSges(7,2) * t210 - mrSges(7,3) * t184 + Ifges(7,1) * t213 + Ifges(7,4) * t215 + Ifges(7,5) * qJDD(6) - qJD(6) * t191) + t274 * (-mrSges(7,1) * t210 + mrSges(7,3) * t186 + Ifges(7,4) * t213 + Ifges(7,2) * t215 + Ifges(7,6) * qJDD(6) + qJD(6) * t193) - pkin(14) * (-m(7) * t210 + t215 * mrSges(7,1) - t213 * mrSges(7,2)) + ((-pkin(14) * (t265 ^ 2 + t274 ^ 2) * mrSges(7,3) + pkin(15) * (t278 ^ 2 + t307) * mrSges(3,3)) * qJD(1) - pkin(14) * (-mrSges(7,1) * t265 - mrSges(7,2) * t274) * qJD(6) + pkin(15) * (-mrSges(3,1) * t278 + mrSges(3,2) * t269) * qJD(2)) * qJD(1); (-t277 * (m(4) * t117 - mrSges(4,2) * t252 + mrSges(4,3) * t139 + t10 * t276 + t167 * t200 - t181 * t255 - t20 * t267) + t268 * (m(4) * t119 + mrSges(4,1) * t252 - mrSges(4,3) * t136 - t167 * t197 + t183 * t255 + t304) + t264 * (m(8) * t118 - mrSges(8,2) * t251 + mrSges(8,3) * t135 + t168 * t196 - t182 * t254 + t300) + t273 * (m(8) * t116 + mrSges(8,1) * t251 - mrSges(8,3) * t138 - t168 * t199 + t180 * t254 + t301)) * pkin(1) + (t278 * t192 + t269 * t194) * qJD(1) + ((t207 * t229 - t208 * t227) * pkin(3) * (mrSges(7,1) * t184 - mrSges(7,2) * t186 + Ifges(7,5) * t213 + Ifges(7,6) * t215 + Ifges(7,3) * qJDD(6) + (t265 * t191 - t274 * t193) * qJD(1)) + (t283 + (-(-t201 * t226 + t202 * t228) * t305 - (t201 * t225 - t202 * t224) * t306) * t174) * (t207 * t265 + t208 * t274) * pkin(7)) / (-t227 * t265 - t229 * t274) / pkin(7) / pkin(3) + Ifges(3,5) * t216 + Ifges(3,6) * t214 - mrSges(3,2) * t185 + mrSges(3,1) * t187 + t281 + t283 + t284 + Ifges(3,3) * qJDD(2); ((-t205 * t228 + t206 * t226) * t305 + (t205 * t224 - t206 * t225) * t306) * t174 + ((-t217 * t238 + t218 * t239) * t288 + (-t238 * t245 - t239 * t246) * t284 * pkin(12)) / (t217 * t246 + t218 * t245) / pkin(12) * pkin(6) + t281; t286;];
tau = t6(:);
