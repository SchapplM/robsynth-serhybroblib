% Calculate vector of cutting torques with Newton-Euler for
% picker2Dm2OL
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
% m [3x11]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = picker2Dm2OL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2OL_invdynm_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:19:03
% EndTime: 2020-05-09 23:19:08
% DurationCPUTime: 1.42s
% Computational Cost: add. (13745->263), mult. (17322->322), div. (0->0), fcn. (9414->22), ass. (0->128)
t265 = pkin(5) * m(6);
t264 = -m(4) - m(10);
t263 = -m(11) - m(5);
t236 = sin(qJ(1));
t245 = cos(qJ(1));
t187 = -t236 * g(1) + t245 * g(2);
t182 = qJDD(1) * pkin(1) + t187;
t189 = t245 * g(1) + t236 * g(2);
t248 = qJD(1) ^ 2;
t183 = -t248 * pkin(1) + t189;
t235 = sin(qJ(2));
t244 = cos(qJ(2));
t169 = t244 * t182 - t235 * t183;
t218 = qJDD(1) + qJDD(2);
t162 = t218 * pkin(3) + t169;
t171 = t235 * t182 + t244 * t183;
t220 = qJD(1) + qJD(2);
t216 = t220 ^ 2;
t164 = -t216 * pkin(3) + t171;
t233 = sin(qJ(4));
t242 = cos(qJ(4));
t145 = t242 * t162 - t233 * t164;
t210 = qJD(4) + t220;
t204 = t210 ^ 2;
t207 = qJDD(4) + t218;
t137 = t207 * pkin(4) + t145;
t147 = t233 * t162 + t242 * t164;
t139 = -t204 * pkin(4) + t147;
t224 = sin(qJ(10));
t226 = cos(qJ(10));
t130 = -t226 * t137 + t224 * t139;
t192 = qJDD(10) + t207;
t199 = qJD(10) + t210;
t196 = t199 ^ 2;
t124 = m(11) * t130 + t192 * mrSges(11,1) - t196 * mrSges(11,2);
t131 = -t224 * t137 - t226 * t139;
t125 = m(11) * t131 - t196 * mrSges(11,1) - t192 * mrSges(11,2);
t254 = -t226 * t124 - t224 * t125;
t114 = m(5) * t145 + t207 * mrSges(5,1) - t204 * mrSges(5,2) + t254;
t115 = m(5) * t147 - t204 * mrSges(5,1) - t207 * mrSges(5,2) + t224 * t124 - t226 * t125;
t262 = t242 * t114 + t233 * t115;
t211 = qJD(3) + t220;
t163 = t218 * pkin(2) + t169;
t165 = -t216 * pkin(2) + t171;
t234 = sin(qJ(3));
t243 = cos(qJ(3));
t146 = -t243 * t163 + t234 * t165;
t208 = qJDD(3) + t218;
t225 = sin(pkin(8));
t227 = cos(pkin(8));
t232 = sin(qJ(5));
t241 = cos(qJ(5));
t180 = t225 * t241 + t227 * t232;
t181 = -t225 * t232 + t227 * t241;
t173 = t180 * g(1) - t181 * g(2);
t174 = -t181 * g(1) - t180 * g(2);
t261 = mrSges(6,1) * t173 - mrSges(6,2) * t174 + Ifges(6,3) * qJDD(5);
t231 = sin(qJ(6));
t240 = cos(qJ(6));
t150 = -t240 * t169 + t231 * t171;
t151 = -t231 * t169 - t240 * t171;
t206 = qJDD(6) + t218;
t260 = mrSges(7,1) * t150 - mrSges(7,2) * t151 + Ifges(7,3) * t206;
t230 = sin(qJ(7));
t239 = cos(qJ(7));
t186 = -t230 * g(1) + t239 * g(2);
t188 = -t239 * g(1) - t230 * g(2);
t259 = mrSges(8,1) * t188 - mrSges(8,2) * t186 + Ifges(8,3) * qJDD(7);
t229 = sin(qJ(8));
t238 = cos(qJ(8));
t168 = -t238 * t182 + t229 * t183;
t170 = -t229 * t182 - t238 * t183;
t217 = qJDD(1) + qJDD(8);
t258 = mrSges(9,1) * t168 - mrSges(9,2) * t170 + Ifges(9,3) * t217;
t138 = t208 * pkin(6) + t146;
t148 = -t234 * t163 - t243 * t165;
t205 = t211 ^ 2;
t140 = -t205 * pkin(6) + t148;
t228 = sin(qJ(9));
t237 = cos(qJ(9));
t132 = -t237 * t138 + t228 * t140;
t133 = -t228 * t138 - t237 * t140;
t197 = qJDD(9) + t208;
t257 = mrSges(10,1) * t132 - mrSges(10,2) * t133 + Ifges(10,3) * t197;
t256 = mrSges(11,1) * t130 - mrSges(11,2) * t131 + Ifges(11,3) * t192;
t200 = qJD(9) + t211;
t198 = t200 ^ 2;
t126 = m(10) * t132 + t197 * mrSges(10,1) - t198 * mrSges(10,2);
t127 = m(10) * t133 - t198 * mrSges(10,1) - t197 * mrSges(10,2);
t253 = -t237 * t126 - t228 * t127;
t116 = m(4) * t146 + t208 * mrSges(4,1) - t205 * mrSges(4,2) + t253;
t117 = m(4) * t148 - t205 * mrSges(4,1) - t208 * mrSges(4,2) + t228 * t126 - t237 * t127;
t255 = -t243 * t116 - t234 * t117;
t252 = pkin(6) * t253 + mrSges(4,1) * t146 - mrSges(4,2) * t148 + Ifges(4,3) * t208 + t257;
t251 = pkin(4) * t254 + mrSges(5,1) * t145 - mrSges(5,2) * t147 + Ifges(5,3) * t207 + t256;
t250 = pkin(2) * t255 + pkin(3) * t262 + mrSges(3,1) * t169 - mrSges(3,2) * t171 + Ifges(3,3) * t218 + t251 + t252 + t260;
t209 = qJD(6) + t220;
t203 = t209 ^ 2;
t141 = m(7) * t150 + t206 * mrSges(7,1) - t203 * mrSges(7,2);
t142 = m(7) * t151 - t203 * mrSges(7,1) - t206 * mrSges(7,2);
t219 = qJD(1) + qJD(8);
t215 = t219 ^ 2;
t249 = mrSges(2,1) * t187 - mrSges(2,2) * t189 + Ifges(2,3) * qJDD(1) + t250 + t258 + (t235 * (m(3) * t171 - t216 * mrSges(3,1) - t218 * mrSges(3,2) - t233 * t114 + t242 * t115 + t234 * t116 - t243 * t117 + t231 * t141 - t240 * t142) + t244 * (m(3) * t169 + t218 * mrSges(3,1) - t216 * mrSges(3,2) - t240 * t141 - t231 * t142 + t255 + t262) - t229 * (m(9) * t170 - t215 * mrSges(9,1) - t217 * mrSges(9,2)) - t238 * (m(9) * t168 + t217 * mrSges(9,1) - t215 * mrSges(9,2))) * pkin(1);
t247 = qJD(5) ^ 2;
t246 = qJD(7) ^ 2;
t176 = -mrSges(8,2) * g(3) - mrSges(8,3) * t188 + Ifges(8,5) * qJDD(7) - t246 * Ifges(8,6);
t175 = mrSges(8,1) * g(3) + mrSges(8,3) * t186 + t246 * Ifges(8,5) + Ifges(8,6) * qJDD(7);
t161 = m(6) * t174 - t247 * mrSges(6,1) - qJDD(5) * mrSges(6,2);
t160 = m(6) * t173 + qJDD(5) * mrSges(6,1) - t247 * mrSges(6,2);
t155 = -mrSges(6,2) * g(3) - mrSges(6,3) * t173 + Ifges(6,5) * qJDD(5) - t247 * Ifges(6,6);
t154 = mrSges(6,1) * g(3) + mrSges(6,3) * t174 + t247 * Ifges(6,5) + Ifges(6,6) * qJDD(5);
t153 = -mrSges(9,2) * g(3) - mrSges(9,3) * t168 + Ifges(9,5) * t217 - t215 * Ifges(9,6);
t152 = mrSges(9,1) * g(3) + mrSges(9,3) * t170 + t215 * Ifges(9,5) + Ifges(9,6) * t217;
t136 = -mrSges(7,2) * g(3) - mrSges(7,3) * t150 + Ifges(7,5) * t206 - t203 * Ifges(7,6);
t135 = mrSges(7,1) * g(3) + mrSges(7,3) * t151 + t203 * Ifges(7,5) + Ifges(7,6) * t206;
t123 = -mrSges(10,2) * g(3) - mrSges(10,3) * t132 + Ifges(10,5) * t197 - t198 * Ifges(10,6);
t122 = mrSges(10,1) * g(3) + mrSges(10,3) * t133 + t198 * Ifges(10,5) + Ifges(10,6) * t197;
t121 = -mrSges(11,2) * g(3) - mrSges(11,3) * t130 + Ifges(11,5) * t192 - t196 * Ifges(11,6);
t120 = mrSges(11,1) * g(3) + mrSges(11,3) * t131 + t196 * Ifges(11,5) + Ifges(11,6) * t192;
t111 = -mrSges(4,2) * g(3) - mrSges(4,3) * t146 + Ifges(4,5) * t208 - t205 * Ifges(4,6) + t228 * t122 - t237 * t123;
t110 = -mrSges(5,2) * g(3) - mrSges(5,3) * t145 + Ifges(5,5) * t207 - t204 * Ifges(5,6) + t224 * t120 - t226 * t121;
t109 = mrSges(4,3) * t148 + t205 * Ifges(4,5) + Ifges(4,6) * t208 - t237 * t122 - t228 * t123 + (pkin(6) * m(10) + mrSges(4,1)) * g(3);
t108 = mrSges(5,3) * t147 + t204 * Ifges(5,5) + Ifges(5,6) * t207 - t226 * t120 - t224 * t121 + (pkin(4) * m(11) + mrSges(5,1)) * g(3);
t105 = -mrSges(3,2) * g(3) - mrSges(3,3) * t169 + Ifges(3,5) * t218 - t216 * Ifges(3,6) - t233 * t108 + t234 * t109 + t242 * t110 - t243 * t111 + t231 * t135 - t240 * t136;
t104 = mrSges(3,3) * t171 + t216 * Ifges(3,5) + Ifges(3,6) * t218 + t242 * t108 - t243 * t109 + t233 * t110 - t234 * t111 - t240 * t135 - t231 * t136 + (-pkin(2) * t264 - pkin(3) * t263 + mrSges(3,1)) * g(3);
t102 = -mrSges(2,2) * g(3) - mrSges(2,3) * t187 + Ifges(2,5) * qJDD(1) - t248 * Ifges(2,6) - t235 * t104 + t244 * t105 + t229 * t152 - t238 * t153;
t101 = mrSges(2,3) * t189 + t248 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t244 * t104 + t235 * t105 - t238 * t152 - t229 * t153 + (mrSges(2,1) + (m(3) + m(7) + m(9) - t263 - t264) * pkin(1)) * g(3);
t1 = [mrSges(1,3) * g(2) + t236 * t101 - t245 * t102 - t180 * t154 + t181 * t155 + t239 * t175 + t230 * t176 + (-t225 * t265 - mrSges(1,2)) * g(3), t102, t105, t111, t110, t155, t136, t176, t153, t123, t121, 0, 0, 0, 0, 0; -mrSges(1,3) * g(1) - t245 * t101 - t236 * t102 + t181 * t154 + t180 * t155 + t230 * t175 - t239 * t176 + (m(8) * pkin(7) + t227 * t265 + mrSges(1,1)) * g(3), t101, t104, t109, t108, t154, t135, t175, t152, t122, t120, 0, 0, 0, 0, 0; t249 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-t225 * (-t180 * t160 + t181 * t161) + t227 * (t181 * t160 + t180 * t161)) * pkin(5) + pkin(7) * (-t239 * (m(8) * t186 - t246 * mrSges(8,1) - qJDD(7) * mrSges(8,2)) + t230 * (m(8) * t188 + qJDD(7) * mrSges(8,1) - t246 * mrSges(8,2))) + t259 + t261, t249, t250, t252, t251, t261, t260, t259, t258, t257, t256, 0, 0, 0, 0, 0;];
m_new = t1;
