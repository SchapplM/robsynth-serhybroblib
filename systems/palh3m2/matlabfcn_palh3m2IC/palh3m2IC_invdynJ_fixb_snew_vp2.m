% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% qJDD [10x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m2IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2IC_invdynJ_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 05:00:13
% EndTime: 2020-05-07 05:00:17
% DurationCPUTime: 3.85s
% Computational Cost: add. (21630->460), mult. (47387->573), div. (14->8), fcn. (36485->39), ass. (0->180)
t232 = qJ(4) + pkin(14);
t230 = (qJ(3) + t232);
t226 = (pkin(15) + t230);
t238 = (qJ(2) - qJ(6));
t223 = (t226 - t238);
t218 = -2 * qJ(7) - pkin(16) + t223;
t224 = (t226 + t238);
t219 = pkin(16) + t224;
t252 = -cos((qJ(8) - t218)) + cos((qJ(8) - t219));
t191 = sin(qJ(7));
t196 = sin(qJ(2));
t199 = cos(qJ(7));
t204 = cos(qJ(2));
t154 = (t191 * t204 + t196 * t199) * qJD(1);
t234 = qJD(1) * qJD(2);
t165 = qJDD(1) * t196 + t204 * t234;
t231 = t196 * t234;
t167 = qJDD(1) * t204 - t231;
t110 = -qJD(7) * t154 - t165 * t191 + t167 * t199;
t195 = sin(qJ(3));
t203 = cos(qJ(3));
t155 = (-t195 * t204 - t196 * t203) * qJD(1);
t111 = -qJD(3) * t155 + t165 * t195 - t167 * t203;
t152 = (-t191 * t196 + t199 * t204) * qJD(1);
t112 = qJD(7) * t152 + t165 * t199 + t167 * t191;
t153 = (t195 * t196 - t203 * t204) * qJD(1);
t113 = qJD(3) * t153 - t165 * t203 - t167 * t195;
t197 = sin(qJ(1));
t205 = cos(qJ(1));
t229 = g(1) * t197 - t205 * g(2);
t160 = -qJDD(1) * pkin(12) - t229;
t130 = t160 + (-t167 + t231) * pkin(1);
t187 = qJD(2) + qJD(7);
t136 = -mrSges(8,2) * t187 + mrSges(8,3) * t152;
t188 = qJD(2) + qJD(3);
t137 = -mrSges(4,2) * t188 + mrSges(4,3) * t153;
t138 = mrSges(8,1) * t187 - mrSges(8,3) * t154;
t139 = mrSges(4,1) * t188 - mrSges(4,3) * t155;
t194 = sin(qJ(4));
t202 = cos(qJ(4));
t124 = t153 * t202 - t155 * t194;
t182 = qJD(4) + t188;
t114 = -mrSges(5,2) * t182 + mrSges(5,3) * t124;
t125 = t153 * t194 + t155 * t202;
t115 = mrSges(5,1) * t182 - mrSges(5,3) * t125;
t193 = sin(qJ(5));
t201 = cos(qJ(5));
t105 = t125 * t201 + t182 * t193;
t123 = qJD(5) - t124;
t67 = -qJD(4) * t125 + t111 * t202 - t113 * t194;
t68 = qJD(4) * t124 + t111 * t194 + t113 * t202;
t81 = t130 + (t155 * t188 - t111) * pkin(4);
t39 = (-t124 * t182 - t68) * pkin(10) + (t125 * t182 - t67) * pkin(8) + t81;
t178 = t182 ^ 2;
t186 = qJDD(2) + qJDD(3);
t180 = qJDD(4) + t186;
t206 = qJD(1) ^ 2;
t225 = -g(1) * t205 - g(2) * t197;
t162 = -pkin(12) * t206 + t225;
t141 = -g(3) * t204 - t162 * t196;
t133 = (t196 * t204 * t206 + qJDD(2)) * pkin(1) + t141;
t143 = -g(3) * t196 + t204 * t162;
t250 = t204 ^ 2;
t134 = (-qJD(2) ^ 2 - t206 * t250) * pkin(1) + t143;
t95 = -t133 * t203 + t195 * t134;
t78 = (t153 * t155 + t186) * pkin(4) + t95;
t97 = -t133 * t195 - t134 * t203;
t85 = (-t153 ^ 2 - t188 ^ 2) * pkin(4) + t97;
t50 = t194 * t78 + t202 * t85;
t91 = -pkin(8) * t124 - pkin(10) * t125;
t41 = -pkin(8) * t178 + pkin(10) * t180 + t124 * t91 + t50;
t30 = -t193 * t41 + t201 * t39;
t104 = -t125 * t193 + t182 * t201;
t47 = qJD(5) * t104 + t180 * t193 + t201 * t68;
t65 = qJDD(5) - t67;
t79 = -mrSges(6,1) * t104 + mrSges(6,2) * t105;
t82 = -mrSges(6,2) * t123 + mrSges(6,3) * t104;
t24 = m(6) * t30 + mrSges(6,1) * t65 - mrSges(6,3) * t47 - t105 * t79 + t123 * t82;
t31 = t193 * t39 + t201 * t41;
t46 = -qJD(5) * t105 + t180 * t201 - t193 * t68;
t83 = mrSges(6,1) * t123 - mrSges(6,3) * t105;
t25 = m(6) * t31 - mrSges(6,2) * t65 + mrSges(6,3) * t46 + t104 * t79 - t123 * t83;
t12 = t193 * t25 + t201 * t24;
t216 = -m(5) * t81 + t67 * mrSges(5,1) - t68 * mrSges(5,2) + t124 * t114 - t125 * t115 - t12;
t189 = sin(pkin(15));
t190 = cos(pkin(15));
t198 = cos(qJ(8));
t246 = sin(qJ(8));
t157 = t189 * t198 - t190 * t246;
t158 = -t189 * t246 - t190 * t198;
t100 = t152 * t158 - t154 * t157;
t101 = t152 * t157 + t154 * t158;
t56 = -qJD(8) * t101 + t110 * t158 - t112 * t157;
t57 = qJD(8) * t100 + t110 * t157 + t112 * t158;
t58 = (-t110 * t190 - t112 * t189 + (-t152 * t189 + t154 * t190) * t187) * pkin(3) + t130;
t181 = qJD(8) + t187;
t98 = -mrSges(9,2) * t181 + mrSges(9,3) * t100;
t99 = mrSges(9,1) * t181 - mrSges(9,3) * t101;
t222 = -m(9) * t58 + t56 * mrSges(9,1) - t57 * mrSges(9,2) + t100 * t98 - t101 * t99;
t251 = t111 * mrSges(4,1) + t110 * mrSges(8,1) - t113 * mrSges(4,2) - t112 * mrSges(8,2) + t152 * t136 + t153 * t137 - t154 * t138 - t155 * t139 + t216 + t222 + (-m(4) - m(8)) * t130;
t228 = -t193 * t24 + t201 * t25;
t90 = -mrSges(5,1) * t124 + mrSges(5,2) * t125;
t10 = m(5) * t50 - mrSges(5,2) * t180 + mrSges(5,3) * t67 - t115 * t182 + t124 * t90 + t228;
t49 = -t194 * t85 + t202 * t78;
t40 = -pkin(8) * t180 - pkin(10) * t178 + t125 * t91 - t49;
t217 = -m(6) * t40 + t46 * mrSges(6,1) - mrSges(6,2) * t47 + t104 * t82 - t105 * t83;
t18 = m(5) * t49 + mrSges(5,1) * t180 - mrSges(5,3) * t68 + t114 * t182 - t125 * t90 + t217;
t247 = t194 * t10 + t202 * t18;
t245 = pkin(3) * t189;
t244 = pkin(3) * t190;
t59 = Ifges(6,5) * t105 + Ifges(6,6) * t104 + Ifges(6,3) * t123;
t61 = Ifges(6,1) * t105 + Ifges(6,4) * t104 + Ifges(6,5) * t123;
t15 = -mrSges(6,1) * t40 + mrSges(6,3) * t31 + Ifges(6,4) * t47 + Ifges(6,2) * t46 + Ifges(6,6) * t65 - t105 * t59 + t123 * t61;
t60 = Ifges(6,4) * t105 + Ifges(6,2) * t104 + Ifges(6,6) * t123;
t16 = mrSges(6,2) * t40 - mrSges(6,3) * t30 + Ifges(6,1) * t47 + Ifges(6,4) * t46 + Ifges(6,5) * t65 + t104 * t59 - t123 * t60;
t87 = Ifges(5,4) * t125 + Ifges(5,2) * t124 + Ifges(5,6) * t182;
t88 = Ifges(5,1) * t125 + Ifges(5,4) * t124 + Ifges(5,5) * t182;
t3 = pkin(8) * t217 + pkin(10) * t228 + mrSges(5,1) * t49 - mrSges(5,2) * t50 + Ifges(5,5) * t68 + Ifges(5,6) * t67 + Ifges(5,3) * t180 - t124 * t88 + t125 * t87 + t201 * t15 + t193 * t16;
t243 = 0.1e1 / pkin(9) * t3;
t185 = qJDD(2) + qJDD(7);
t179 = qJDD(8) + t185;
t122 = (-t152 * t190 - t154 * t189) * pkin(3);
t184 = t187 ^ 2;
t94 = t199 * t133 - t134 * t191;
t69 = -t122 * t154 + (t184 * t189 + t185 * t190) * pkin(3) + t94;
t96 = t191 * t133 + t199 * t134;
t70 = t122 * t152 + (-t184 * t190 + t185 * t189) * pkin(3) + t96;
t43 = -t157 * t70 + t158 * t69;
t76 = -mrSges(9,1) * t100 + mrSges(9,2) * t101;
t37 = m(9) * t43 + mrSges(9,1) * t179 - mrSges(9,3) * t57 - t101 * t76 + t181 * t98;
t44 = t157 * t69 + t158 * t70;
t38 = m(9) * t44 - mrSges(9,2) * t179 + mrSges(9,3) * t56 + t100 * t76 - t181 * t99;
t242 = t157 * t38 + t158 * t37;
t241 = -t157 * t37 + t158 * t38;
t72 = Ifges(9,4) * t101 + Ifges(9,2) * t100 + Ifges(9,6) * t181;
t73 = Ifges(9,1) * t101 + Ifges(9,4) * t100 + Ifges(9,5) * t181;
t26 = mrSges(9,1) * t43 - mrSges(9,2) * t44 + Ifges(9,5) * t57 + Ifges(9,6) * t56 + Ifges(9,3) * t179 - t100 * t73 + t101 * t72;
t240 = 0.1e1 / pkin(7) * t26;
t237 = qJ(7) + qJ(8);
t235 = qJ(7) + pkin(16);
t233 = qJD(1) * qJD(6);
t227 = t235 + t238;
t221 = -qJ(7) + t223;
t220 = -qJ(7) + t224;
t215 = mrSges(6,1) * t30 - mrSges(6,2) * t31 + Ifges(6,5) * t47 + Ifges(6,6) * t46 + Ifges(6,3) * t65 - t104 * t61 + t105 * t60;
t118 = Ifges(8,4) * t154 + Ifges(8,2) * t152 + Ifges(8,6) * t187;
t120 = Ifges(8,1) * t154 + Ifges(8,4) * t152 + Ifges(8,5) * t187;
t213 = mrSges(8,1) * t94 - mrSges(8,2) * t96 + Ifges(8,5) * t112 + Ifges(8,6) * t110 + Ifges(8,3) * t185 + t154 * t118 - t152 * t120 + t241 * t245 + t242 * t244 + t26;
t119 = Ifges(4,4) * t155 + Ifges(4,2) * t153 + Ifges(4,6) * t188;
t121 = Ifges(4,1) * t155 + Ifges(4,4) * t153 + Ifges(4,5) * t188;
t211 = pkin(4) * t247 + mrSges(4,1) * t95 - mrSges(4,2) * t97 + Ifges(4,5) * t113 + Ifges(4,6) * t111 + Ifges(4,3) * t186 + t155 * t119 - t153 * t121 + t3;
t210 = 0.1e1 / pkin(2);
t200 = cos(qJ(6));
t192 = sin(qJ(6));
t171 = sin(t227);
t166 = qJDD(1) * t200 - t192 * t233;
t164 = qJDD(1) * t192 + t200 * t233;
t163 = pkin(6) * t206 + t225;
t161 = qJDD(1) * pkin(6) - t229;
t151 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t196 + Ifges(3,4) * t204) * qJD(1);
t150 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t192 + Ifges(7,4) * t200) * qJD(1);
t149 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t196 + Ifges(3,2) * t204) * qJD(1);
t148 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t192 + Ifges(7,2) * t200) * qJD(1);
t142 = -g(3) * t192 + t163 * t200;
t140 = -g(3) * t200 - t163 * t192;
t129 = -mrSges(4,1) * t153 + mrSges(4,2) * t155;
t128 = -mrSges(8,1) * t152 + mrSges(8,2) * t154;
t117 = Ifges(4,5) * t155 + Ifges(4,6) * t153 + Ifges(4,3) * t188;
t116 = Ifges(8,5) * t154 + Ifges(8,6) * t152 + Ifges(8,3) * t187;
t86 = Ifges(5,5) * t125 + Ifges(5,6) * t124 + Ifges(5,3) * t182;
t71 = Ifges(9,5) * t101 + Ifges(9,6) * t100 + Ifges(9,3) * t181;
t28 = mrSges(9,2) * t58 - mrSges(9,3) * t43 + Ifges(9,1) * t57 + Ifges(9,4) * t56 + Ifges(9,5) * t179 + t100 * t71 - t181 * t72;
t27 = -mrSges(9,1) * t58 + mrSges(9,3) * t44 + Ifges(9,4) * t57 + Ifges(9,2) * t56 + Ifges(9,6) * t179 - t101 * t71 + t181 * t73;
t8 = mrSges(8,2) * t130 - mrSges(8,3) * t94 + Ifges(8,1) * t112 + Ifges(8,4) * t110 + Ifges(8,5) * t185 + t116 * t152 - t118 * t187 - t157 * t27 + t158 * t28 + t222 * t245;
t7 = -mrSges(8,1) * t130 + mrSges(8,3) * t96 + Ifges(8,4) * t112 + Ifges(8,2) * t110 + Ifges(8,6) * t185 - t116 * t154 + t120 * t187 + t157 * t28 + t158 * t27 + t222 * t244;
t5 = -pkin(8) * t12 - mrSges(5,1) * t81 + mrSges(5,3) * t50 + Ifges(5,4) * t68 + Ifges(5,2) * t67 + Ifges(5,6) * t180 - t125 * t86 + t182 * t88 - t215;
t4 = -pkin(10) * t12 + mrSges(5,2) * t81 - mrSges(5,3) * t49 + Ifges(5,1) * t68 + Ifges(5,4) * t67 + Ifges(5,5) * t180 + t124 * t86 - t15 * t193 + t16 * t201 - t182 * t87;
t2 = mrSges(4,2) * t130 - mrSges(4,3) * t95 + Ifges(4,1) * t113 + Ifges(4,4) * t111 + Ifges(4,5) * t186 + t117 * t153 - t119 * t188 - t194 * t5 + t202 * t4;
t1 = pkin(4) * t216 - mrSges(4,1) * t130 + mrSges(4,3) * t97 + Ifges(4,4) * t113 + Ifges(4,2) * t111 + Ifges(4,6) * t186 - t155 * t117 + t188 * t121 + t194 * t4 + t202 * t5;
t6 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t229 - mrSges(2,2) * t225 + t196 * (mrSges(3,2) * t160 - mrSges(3,3) * t141 + Ifges(3,1) * t165 + Ifges(3,4) * t167 + Ifges(3,5) * qJDD(2) - qJD(2) * t149 + t1 * t195 - t191 * t7 + t199 * t8 - t2 * t203) + t204 * (t251 * pkin(1) - mrSges(3,1) * t160 + mrSges(3,3) * t143 + Ifges(3,4) * t165 + Ifges(3,2) * t167 + Ifges(3,6) * qJDD(2) + qJD(2) * t151 - t203 * t1 + t191 * t8 - t195 * t2 + t199 * t7) + pkin(12) * (-m(3) * t160 + t167 * mrSges(3,1) - t165 * mrSges(3,2) + t251) + t192 * (mrSges(7,2) * t161 - mrSges(7,3) * t140 + Ifges(7,1) * t164 + Ifges(7,4) * t166 + Ifges(7,5) * qJDD(6) - qJD(6) * t148) + t200 * (-mrSges(7,1) * t161 + mrSges(7,3) * t142 + Ifges(7,4) * t164 + Ifges(7,2) * t166 + Ifges(7,6) * qJDD(6) + qJD(6) * t150) - pkin(6) * (-m(7) * t161 + t166 * mrSges(7,1) - t164 * mrSges(7,2)) + ((-pkin(6) * (t192 ^ 2 + t200 ^ 2) * mrSges(7,3) + pkin(12) * (t196 ^ 2 + t250) * mrSges(3,3)) * qJD(1) - pkin(6) * (-mrSges(7,1) * t192 - mrSges(7,2) * t200) * qJD(6) + pkin(12) * (-mrSges(3,1) * t196 - mrSges(3,2) * t204) * qJD(2)) * qJD(1); mrSges(3,1) * t141 - mrSges(3,2) * t143 + Ifges(3,5) * t165 + Ifges(3,6) * t167 + Ifges(3,3) * qJDD(2) + t211 + t213 + (-t195 * (m(4) * t97 - mrSges(4,2) * t186 + mrSges(4,3) * t111 + t10 * t202 + t129 * t153 - t139 * t188 - t18 * t194) - t203 * (m(4) * t95 + mrSges(4,1) * t186 - mrSges(4,3) * t113 - t129 * t155 + t137 * t188 + t247) + t191 * (m(8) * t96 - mrSges(8,2) * t185 + mrSges(8,3) * t110 + t128 * t152 - t138 * t187 + t241) + t199 * (m(8) * t94 + mrSges(8,1) * t185 - mrSges(8,3) * t112 - t128 * t154 + t136 * t187 + t242)) * pkin(1) + (pkin(3) * (pkin(1) * (cos((-qJ(8) + t238)) - cos((qJ(8) + t238))) + (cos((-qJ(8) + t227)) - cos((qJ(8) + t227))) * pkin(2)) * t243 + (((cos(t221) - cos(t220)) * pkin(1) + (cos(t218) - cos(t219)) * pkin(2)) * pkin(3) + ((-cos((-qJ(8) + t221)) + cos((-qJ(8) + t220))) * pkin(1) + t252 * pkin(2)) * pkin(7)) * t240) / t252 * t210 + (t196 * t149 - t204 * t151) * qJD(1) + (pkin(1) * sin(t235) / pkin(5) * (mrSges(7,1) * t140 - mrSges(7,2) * t142 + Ifges(7,5) * t164 + Ifges(7,6) * t166 + Ifges(7,3) * qJDD(6) + (t192 * t148 - t200 * t150) * qJD(1)) + (-pkin(2) * t171 - pkin(1) * sin(t238)) * t210 * t213) / t171; -pkin(9) * t243 + t211 + (-(cos(0.2e1 * qJ(3) + t232) - cos((2 * qJ(7)) - 0.2e1 * pkin(15) + (2 * qJ(8)) + t232)) / (cos((2 * t230)) - cos(0.2e1 * pkin(15) - (2 * t237))) * t243 + sin(t232) / sin((t226 - t237)) * t240) * pkin(4); t215;];
tau = t6(:);
