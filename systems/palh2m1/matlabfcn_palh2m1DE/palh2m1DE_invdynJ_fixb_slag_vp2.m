% Calculate vector of inverse dynamics joint torques for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m1DE_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:14
% DurationCPUTime: 2.31s
% Computational Cost: add. (1793->293), mult. (2328->411), div. (0->0), fcn. (1097->18), ass. (0->177)
t147 = qJD(1) ^ 2;
t138 = cos(qJ(3));
t104 = pkin(3) * t138 + pkin(2);
t139 = cos(qJ(2));
t134 = sin(qJ(3));
t135 = sin(qJ(2));
t197 = t134 * t135;
t50 = -pkin(3) * t197 + t104 * t139 + pkin(1);
t237 = -m(5) * t50 - mrSges(5,1);
t246 = t147 * t237;
t245 = m(5) + m(6);
t131 = qJ(4) + qJ(2);
t120 = qJ(3) + t131;
t101 = sin(t120);
t195 = -qJ(4) + qJ(2);
t180 = qJ(3) + t195;
t172 = sin(t180);
t129 = qJD(2) + qJD(3);
t175 = qJD(1) + t129;
t151 = t175 * t172;
t150 = -t151 / 0.2e1;
t176 = qJD(1) - t129;
t157 = -t176 * t101 / 0.2e1;
t127 = qJD(4) + qJD(2);
t113 = qJD(3) + t127;
t219 = pkin(3) * t113;
t181 = t219 / 0.2e1;
t244 = t101 * t181 - (t150 + t157) * pkin(3);
t102 = cos(t120);
t103 = cos(t180);
t168 = t103 * t175;
t158 = -t168 / 0.2e1;
t182 = -t219 / 0.2e1;
t243 = pkin(3) * t158 - t102 * t182;
t115 = sin(t131);
t178 = sin(t195);
t186 = qJD(1) + qJD(2);
t155 = t186 * t178;
t152 = -t155 / 0.2e1;
t187 = qJD(1) - qJD(2);
t165 = -t187 * t115 / 0.2e1;
t225 = pkin(2) * t127;
t183 = t225 / 0.2e1;
t128 = -qJD(4) + qJD(2);
t114 = qJD(3) + t128;
t161 = -t172 / 0.2e1;
t92 = pkin(3) * t161;
t68 = t114 * t92;
t171 = -t178 / 0.2e1;
t95 = pkin(2) * t171;
t213 = t128 * t95 + t68;
t242 = t115 * t183 - (t152 + t165) * pkin(2) + t213 + t244;
t117 = cos(t131);
t118 = cos(t195);
t174 = t118 * t186;
t166 = -t174 / 0.2e1;
t184 = -t225 / 0.2e1;
t167 = t117 * t187 / 0.2e1;
t159 = t102 * t176 / 0.2e1;
t73 = pkin(3) * t159;
t211 = pkin(2) * t167 + t73;
t226 = pkin(2) * t118;
t220 = pkin(3) * t103;
t69 = -t114 * t220 / 0.2e1;
t212 = t69 - t128 * t226 / 0.2e1;
t241 = pkin(2) * t166 - t117 * t184 + t211 - t212 + t243;
t240 = t68 + t244;
t239 = -t69 + t73 + t243;
t238 = qJDD(1) / 0.2e1;
t196 = t134 * t139;
t163 = t135 * t138 + t196;
t154 = t163 * qJD(3);
t173 = t245 * pkin(3) + mrSges(4,1);
t63 = mrSges(4,2) * t138 + t134 * t173 + mrSges(3,2);
t47 = t138 * t173 - mrSges(4,2) * t134 + mrSges(3,1) + (m(4) + t245) * pkin(2);
t62 = t163 * qJD(1);
t235 = -t62 / 0.2e1;
t123 = qJDD(2) + qJDD(3);
t233 = t238 - t123 / 0.2e1;
t232 = t238 - qJDD(2) / 0.2e1;
t231 = t135 / 0.2e1;
t162 = -t138 * t139 + t197;
t61 = t162 * qJD(1);
t230 = mrSges(4,3) * t61;
t229 = mrSges(4,3) * t62;
t228 = pkin(2) * t115;
t227 = pkin(2) * t117;
t224 = pkin(2) * t135;
t223 = pkin(2) * t139;
t222 = pkin(3) * t101;
t221 = pkin(3) * t102;
t132 = qJ(2) + qJ(3);
t116 = sin(t132);
t218 = pkin(3) * t116;
t119 = cos(t132);
t217 = pkin(3) * t119;
t122 = t129 ^ 2;
t146 = qJD(2) ^ 2;
t33 = -qJDD(2) * t223 + t122 * t218 - t123 * t217 + t146 * t224;
t216 = t33 * (t217 + t223);
t215 = t62 * Ifges(4,4);
t124 = qJDD(1) + qJDD(4);
t110 = qJDD(1) + t123;
t125 = qJDD(1) + qJDD(2);
t133 = sin(qJ(4));
t137 = cos(qJ(4));
t142 = pkin(1) + pkin(4);
t189 = qJD(1) * qJD(4);
t3 = (qJDD(1) * t137 - t133 * t189) * t142 + (t114 * t150 + t113 * t157 + t110 * t103 / 0.2e1 + t102 * t233) * pkin(3) + (t128 * t152 + t127 * t165 + t125 * t118 / 0.2e1 + t117 * t232) * pkin(2);
t214 = t3 * mrSges(6,1) + Ifges(6,3) * t124;
t100 = mrSges(6,1) * g(2) - mrSges(6,2) * g(1);
t99 = mrSges(6,1) * g(1) + mrSges(6,2) * g(2);
t210 = -t100 * t137 + t99 * t133;
t209 = t92 + t95;
t94 = t220 / 0.2e1;
t208 = t94 + t226 / 0.2e1;
t207 = Ifges(3,4) * t135;
t206 = Ifges(3,4) * t139;
t205 = Ifges(3,5) * t139;
t204 = Ifges(3,6) * t135;
t203 = pkin(2) * qJD(2);
t199 = pkin(1) * qJDD(1);
t198 = t104 * t135;
t105 = pkin(1) + t223;
t194 = qJD(1) * t105;
t193 = qJD(1) * t142;
t192 = qJD(4) * t142;
t191 = qJDD(1) * mrSges(5,3);
t190 = qJD(1) * qJD(2);
t188 = qJD(2) * qJD(3);
t185 = t33 * t217;
t177 = t135 * t190;
t136 = sin(qJ(1));
t140 = cos(qJ(1));
t170 = g(1) * t140 + g(2) * t136;
t169 = g(1) * t136 - g(2) * t140;
t164 = -t100 * t133 - t137 * t99;
t160 = (Ifges(3,6) * qJD(2) + (-t139 * Ifges(3,2) - t207) * qJD(1)) * t231 - t139 * (Ifges(3,5) * qJD(2) + (-t135 * Ifges(3,1) - t206) * qJD(1)) / 0.2e1;
t153 = t162 * qJD(3);
t79 = -qJDD(1) * t139 + t177;
t80 = -qJDD(1) * t135 - t139 * t190;
t15 = qJD(1) * t153 + t134 * t79 + t138 * t80;
t16 = qJD(1) * t154 - t134 * t80 + t138 * t79;
t22 = t61 * Ifges(4,2) + t129 * Ifges(4,6) - t215;
t51 = Ifges(4,4) * t61;
t23 = -t62 * Ifges(4,1) + t129 * Ifges(4,5) + t51;
t74 = (qJDD(2) * t138 - t134 * t188) * pkin(2);
t75 = (qJDD(2) * t134 + t138 * t188) * pkin(2);
t149 = -t75 * mrSges(4,2) + t22 * t235 + Ifges(4,3) * t123 + t62 * (Ifges(4,1) * t61 + t215) / 0.2e1 + Ifges(4,6) * t16 + Ifges(4,5) * t15 - t129 * (Ifges(4,5) * t61 + Ifges(4,6) * t62) / 0.2e1 + t138 * t203 * t230 + t74 * mrSges(4,1) - (Ifges(4,2) * t62 + t23 + t51) * t61 / 0.2e1;
t148 = (-Ifges(3,1) * t139 + t207) * t231 + t139 * (Ifges(3,2) * t135 - t206) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t135 - mrSges(3,2) * t139);
t130 = qJD(1) + qJD(4);
t112 = mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t93 = -t221 / 0.2e1;
t91 = -t222 / 0.2e1;
t76 = t142 * m(6) + mrSges(2,1) + mrSges(5,1) + (m(5) + m(3) + m(4)) * pkin(1);
t58 = -pkin(2) * t177 + qJDD(1) * t105;
t53 = t94 + t93;
t52 = t92 + t91;
t49 = mrSges(4,2) * t170 + g(3) * t173;
t46 = mrSges(4,1) * t129 + t229;
t45 = -mrSges(4,2) * t129 + t230;
t44 = -g(3) * mrSges(4,2) + t170 * t173;
t37 = qJD(2) * t163 + t154;
t36 = qJD(2) * t162 + t153;
t35 = t93 - t227 / 0.2e1 + t208;
t34 = t91 - t228 / 0.2e1 + t209;
t32 = (-t116 * t123 - t119 * t122) * pkin(3) + (-qJDD(2) * t135 - t139 * t146) * pkin(2);
t28 = -mrSges(4,1) * t61 - mrSges(4,2) * t62;
t27 = -mrSges(4,1) * t62 + mrSges(4,2) * t61;
t25 = t221 / 0.2e1 + t227 / 0.2e1 + t137 * t142 + t208;
t24 = t222 / 0.2e1 + t228 / 0.2e1 + t133 * t142 + t209;
t10 = pkin(3) * t168 / 0.2e1 + pkin(2) * t174 / 0.2e1 + t137 * t193 + t211;
t9 = -t133 * t193 + (t151 / 0.2e1 + t157) * pkin(3) + (t155 / 0.2e1 + t165) * pkin(2);
t8 = t102 * t181 + t117 * t183 + t137 * t192 + t212;
t7 = t101 * t182 + t115 * t184 - t133 * t192 + t213;
t2 = (qJDD(1) * t133 + t137 * t189) * t142 + (t101 * t233 + t110 * t161 + t113 * t159 + t114 * t158) * pkin(3) + (t115 * t232 + t125 * t171 + t127 * t167 + t128 * t166) * pkin(2);
t1 = [(m(4) * t58 - mrSges(4,1) * t16 + mrSges(4,2) * t15) * t105 + (-mrSges(4,1) * t58 + mrSges(4,3) * t75 + Ifges(4,4) * t15 + Ifges(4,2) * t16 + Ifges(4,6) * t123) * t162 - (mrSges(4,2) * t58 - mrSges(4,3) * t74 + Ifges(4,1) * t15 + Ifges(4,4) * t16 + Ifges(4,5) * t123) * t163 + (g(1) * t76 + g(2) * t112 - t164) * t136 + m(6) * (t10 * t7 + t2 * t24 + t25 * t3 - t8 * t9) + t129 * (Ifges(4,5) * t36 + Ifges(4,6) * t37) / 0.2e1 + (-t124 * t24 - t130 * t8 - t2) * mrSges(6,2) + pkin(1) * (-mrSges(3,1) * t79 + mrSges(3,2) * t80) + (Ifges(4,1) * t36 + Ifges(4,4) * t37) * t235 + t61 * (Ifges(4,4) * t36 + Ifges(4,2) * t37) / 0.2e1 - t32 * mrSges(5,3) + t36 * t23 / 0.2e1 + t37 * t22 / 0.2e1 + (t124 * t25 + t130 * t7) * mrSges(6,1) + (mrSges(3,1) * t199 - Ifges(3,4) * t80 - Ifges(3,2) * t79 - Ifges(3,6) * qJDD(2) + t169 * t47) * t139 + (-mrSges(3,2) * t199 - Ifges(3,1) * t80 - Ifges(3,4) * t79 - Ifges(3,5) * qJDD(2) - t169 * t63) * t135 + ((-t205 / 0.2e1 + t204 / 0.2e1) * qJD(2) + ((-m(4) * t194 - t28) * t135 + (t134 * t37 - t138 * t36) * mrSges(4,3)) * pkin(2) + t160) * qJD(2) + (g(1) * t112 - g(2) * t76 + t210) * t140 + t214 + (m(3) * pkin(1) ^ 2 + Ifges(5,2) + Ifges(2,3) + (mrSges(5,1) - t237) * t50) * qJDD(1) + (-0.2e1 * t237 * (-qJD(2) * t198 + (-qJD(2) * t196 - t154) * pkin(3)) + t105 * (-mrSges(4,1) * t37 + mrSges(4,2) * t36) - t148 * qJD(2)) * qJD(1); (mrSges(6,1) * t35 - mrSges(6,2) * t34) * t124 + ((mrSges(4,1) * t123 - mrSges(4,3) * t15 + qJD(3) * t45) * t138 + (-mrSges(4,2) * t123 - qJD(3) * t46 + (-qJD(2) * t62 + t16) * mrSges(4,3)) * t134 + (t105 * t147 * t135 + t134 * t75 + t138 * t74) * m(4)) * pkin(2) - t216 * m(5) + (mrSges(6,1) * t242 + mrSges(6,2) * t241) * t130 + t149 + t139 * (t47 * g(3) + t170 * t63) + (-t63 * g(3) + t170 * t47) * t135 + Ifges(3,6) * t79 + Ifges(3,5) * t80 + (-qJD(2) * (t204 - t205) / 0.2e1 + t28 * t224 - t105 * t27 + t148 * qJD(1) - t160) * qJD(1) + Ifges(3,3) * qJDD(2) + (m(5) * t32 - t191) * (-t218 - t224) + (-pkin(3) * t196 - t198) * t246 + (t10 * t242 + t2 * t34 + t241 * t9 + t3 * t35 - t216) * m(6); (mrSges(6,1) * t240 + mrSges(6,2) * t239) * t130 + (-t138 * t45 + (t46 - t229) * t134) * t203 + m(5) * (-t218 * t32 - t185) + (mrSges(6,1) * t53 - mrSges(6,2) * t52) * t124 + t149 + (t134 * t44 + t138 * t49) * t139 + t135 * (-t134 * t49 + t138 * t44) + t191 * t218 - t27 * t194 - t163 * pkin(3) * t246 + (t10 * t240 + t2 * t52 + t239 * t9 + t3 * t53 - t185) * m(6); -t9 * t130 * mrSges(6,1) + t210 * t140 - t164 * t136 + (t10 * t130 - t2) * mrSges(6,2) + t214;];
tau = t1;
