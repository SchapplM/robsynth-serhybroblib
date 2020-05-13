% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% palh3m1IC
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m1IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1IC_invdynJ_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:32:06
% EndTime: 2020-04-20 17:32:09
% DurationCPUTime: 3.76s
% Computational Cost: add. (21704->436), mult. (47521->582), div. (16->6), fcn. (36580->26), ass. (0->181)
t204 = sin(qJ(7));
t209 = sin(qJ(2));
t212 = cos(qJ(7));
t217 = cos(qJ(2));
t157 = (t204 * t217 + t209 * t212) * qJD(1);
t234 = qJD(1) * qJD(2);
t172 = qJDD(1) * t209 + t217 * t234;
t232 = t209 * t234;
t174 = qJDD(1) * t217 - t232;
t110 = -qJD(7) * t157 - t172 * t204 + t174 * t212;
t208 = sin(qJ(3));
t216 = cos(qJ(3));
t158 = (-t208 * t217 - t209 * t216) * qJD(1);
t111 = -qJD(3) * t158 + t172 * t208 - t174 * t216;
t155 = (-t204 * t209 + t212 * t217) * qJD(1);
t112 = qJD(7) * t155 + t172 * t212 + t174 * t204;
t156 = (t208 * t209 - t216 * t217) * qJD(1);
t113 = qJD(3) * t156 - t172 * t216 - t174 * t208;
t210 = sin(qJ(1));
t218 = cos(qJ(1));
t231 = g(1) * t210 - t218 * g(2);
t165 = -qJDD(1) * pkin(12) - t231;
t132 = t165 + (-t174 + t232) * pkin(1);
t199 = qJD(2) + qJD(7);
t139 = -mrSges(8,2) * t199 + mrSges(8,3) * t155;
t200 = qJD(2) + qJD(3);
t140 = -mrSges(4,2) * t200 + mrSges(4,3) * t156;
t141 = mrSges(8,1) * t199 - mrSges(8,3) * t157;
t142 = mrSges(4,1) * t200 - mrSges(4,3) * t158;
t207 = sin(qJ(4));
t215 = cos(qJ(4));
t126 = t156 * t215 - t158 * t207;
t191 = qJD(4) + t200;
t114 = -mrSges(5,2) * t191 + mrSges(5,3) * t126;
t127 = t156 * t207 + t158 * t215;
t115 = mrSges(5,1) * t191 - mrSges(5,3) * t127;
t206 = sin(qJ(5));
t214 = cos(qJ(5));
t105 = t127 * t214 + t191 * t206;
t125 = qJD(5) - t126;
t67 = -qJD(4) * t127 + t111 * t215 - t113 * t207;
t68 = qJD(4) * t126 + t111 * t207 + t113 * t215;
t81 = t132 + (t158 * t200 - t111) * pkin(4);
t39 = (-t126 * t191 - t68) * pkin(10) + (t127 * t191 - t67) * pkin(8) + t81;
t187 = t191 ^ 2;
t198 = qJDD(2) + qJDD(3);
t189 = qJDD(4) + t198;
t219 = qJD(1) ^ 2;
t229 = -g(1) * t218 - g(2) * t210;
t169 = -pkin(12) * t219 + t229;
t144 = -g(3) * t217 - t169 * t209;
t135 = (t209 * t217 * t219 + qJDD(2)) * pkin(1) + t144;
t146 = -g(3) * t209 + t217 * t169;
t243 = t217 ^ 2;
t137 = (-qJD(2) ^ 2 - t219 * t243) * pkin(1) + t146;
t95 = -t135 * t216 + t208 * t137;
t78 = (t156 * t158 + t198) * pkin(4) + t95;
t97 = -t135 * t208 - t137 * t216;
t85 = (-t156 ^ 2 - t200 ^ 2) * pkin(4) + t97;
t50 = t207 * t78 + t215 * t85;
t91 = -pkin(8) * t126 - pkin(10) * t127;
t41 = -pkin(8) * t187 + pkin(10) * t189 + t126 * t91 + t50;
t30 = -t206 * t41 + t214 * t39;
t104 = -t127 * t206 + t191 * t214;
t47 = qJD(5) * t104 + t189 * t206 + t214 * t68;
t65 = qJDD(5) - t67;
t79 = -mrSges(6,1) * t104 + mrSges(6,2) * t105;
t82 = -mrSges(6,2) * t125 + mrSges(6,3) * t104;
t24 = m(6) * t30 + mrSges(6,1) * t65 - mrSges(6,3) * t47 - t105 * t79 + t125 * t82;
t31 = t206 * t39 + t214 * t41;
t46 = -qJD(5) * t105 + t189 * t214 - t206 * t68;
t83 = mrSges(6,1) * t125 - mrSges(6,3) * t105;
t25 = m(6) * t31 - mrSges(6,2) * t65 + mrSges(6,3) * t46 + t104 * t79 - t125 * t83;
t12 = t206 * t25 + t214 * t24;
t225 = -m(5) * t81 + t67 * mrSges(5,1) - t68 * mrSges(5,2) + t126 * t114 - t127 * t115 - t12;
t202 = sin(pkin(15));
t203 = cos(pkin(15));
t211 = cos(qJ(8));
t239 = sin(qJ(8));
t161 = t202 * t211 - t203 * t239;
t162 = -t202 * t239 - t203 * t211;
t100 = t155 * t162 - t157 * t161;
t101 = t155 * t161 + t157 * t162;
t56 = -qJD(8) * t101 + t110 * t162 - t112 * t161;
t57 = qJD(8) * t100 + t110 * t161 + t112 * t162;
t58 = (-t110 * t203 - t112 * t202 + (-t155 * t202 + t157 * t203) * t199) * pkin(3) + t132;
t190 = qJD(8) + t199;
t98 = -mrSges(9,2) * t190 + mrSges(9,3) * t100;
t99 = mrSges(9,1) * t190 - mrSges(9,3) * t101;
t227 = -m(9) * t58 + t56 * mrSges(9,1) - t57 * mrSges(9,2) + t100 * t98 - t101 * t99;
t244 = t111 * mrSges(4,1) + t110 * mrSges(8,1) - t113 * mrSges(4,2) - t112 * mrSges(8,2) + t155 * t139 + t156 * t140 - t157 * t141 - t158 * t142 + t225 + t227 + (-m(4) - m(8)) * t132;
t49 = -t207 * t85 + t215 * t78;
t40 = -pkin(8) * t189 - pkin(10) * t187 + t127 * t91 - t49;
t59 = Ifges(6,5) * t105 + Ifges(6,6) * t104 + Ifges(6,3) * t125;
t61 = Ifges(6,1) * t105 + Ifges(6,4) * t104 + Ifges(6,5) * t125;
t15 = -mrSges(6,1) * t40 + mrSges(6,3) * t31 + Ifges(6,4) * t47 + Ifges(6,2) * t46 + Ifges(6,6) * t65 - t105 * t59 + t125 * t61;
t60 = Ifges(6,4) * t105 + Ifges(6,2) * t104 + Ifges(6,6) * t125;
t16 = mrSges(6,2) * t40 - mrSges(6,3) * t30 + Ifges(6,1) * t47 + Ifges(6,4) * t46 + Ifges(6,5) * t65 + t104 * t59 - t125 * t60;
t226 = -m(6) * t40 + t46 * mrSges(6,1) - mrSges(6,2) * t47 + t104 * t82 - t105 * t83;
t230 = -t206 * t24 + t214 * t25;
t87 = Ifges(5,4) * t127 + Ifges(5,2) * t126 + Ifges(5,6) * t191;
t88 = Ifges(5,1) * t127 + Ifges(5,4) * t126 + Ifges(5,5) * t191;
t3 = pkin(8) * t226 + pkin(10) * t230 + mrSges(5,1) * t49 - mrSges(5,2) * t50 + Ifges(5,5) * t68 + Ifges(5,6) * t67 + Ifges(5,3) * t189 - t126 * t88 + t127 * t87 + t214 * t15 + t206 * t16;
t242 = pkin(7) * t3;
t197 = qJDD(2) + qJDD(7);
t188 = qJDD(8) + t197;
t124 = (-t155 * t203 - t157 * t202) * pkin(3);
t196 = t199 ^ 2;
t94 = t212 * t135 - t137 * t204;
t69 = -t124 * t157 + (t196 * t202 + t197 * t203) * pkin(3) + t94;
t96 = t204 * t135 + t212 * t137;
t70 = t124 * t155 + (-t196 * t203 + t197 * t202) * pkin(3) + t96;
t43 = -t161 * t70 + t162 * t69;
t44 = t161 * t69 + t162 * t70;
t72 = Ifges(9,4) * t101 + Ifges(9,2) * t100 + Ifges(9,6) * t190;
t73 = Ifges(9,1) * t101 + Ifges(9,4) * t100 + Ifges(9,5) * t190;
t26 = mrSges(9,1) * t43 - mrSges(9,2) * t44 + Ifges(9,5) * t57 + Ifges(9,6) * t56 + Ifges(9,3) * t188 - t100 * t73 + t101 * t72;
t241 = pkin(9) * t26;
t90 = -mrSges(5,1) * t126 + mrSges(5,2) * t127;
t10 = m(5) * t50 - mrSges(5,2) * t189 + mrSges(5,3) * t67 - t115 * t191 + t126 * t90 + t230;
t18 = m(5) * t49 + mrSges(5,1) * t189 - mrSges(5,3) * t68 + t114 * t191 - t127 * t90 + t226;
t240 = t207 * t10 + t215 * t18;
t238 = pkin(3) * t202;
t237 = pkin(3) * t203;
t76 = -mrSges(9,1) * t100 + mrSges(9,2) * t101;
t37 = m(9) * t43 + mrSges(9,1) * t188 - mrSges(9,3) * t57 - t101 * t76 + t190 * t98;
t38 = m(9) * t44 - mrSges(9,2) * t188 + mrSges(9,3) * t56 + t100 * t76 - t190 * t99;
t236 = t161 * t38 + t162 * t37;
t235 = -t161 * t37 + t162 * t38;
t201 = -qJ(7) + pkin(15);
t233 = qJD(1) * qJD(6);
t224 = mrSges(6,1) * t30 - mrSges(6,2) * t31 + Ifges(6,5) * t47 + Ifges(6,6) * t46 + Ifges(6,3) * t65 - t104 * t61 + t105 * t60;
t120 = Ifges(8,4) * t157 + Ifges(8,2) * t155 + Ifges(8,6) * t199;
t122 = Ifges(8,1) * t157 + Ifges(8,4) * t155 + Ifges(8,5) * t199;
t222 = mrSges(8,1) * t94 - mrSges(8,2) * t96 + Ifges(8,5) * t112 + Ifges(8,6) * t110 + Ifges(8,3) * t197 + t157 * t120 - t155 * t122 + t235 * t238 + t236 * t237 + t26;
t121 = Ifges(4,4) * t158 + Ifges(4,2) * t156 + Ifges(4,6) * t200;
t123 = Ifges(4,1) * t158 + Ifges(4,4) * t156 + Ifges(4,5) * t200;
t220 = pkin(4) * t240 + mrSges(4,1) * t95 - mrSges(4,2) * t97 + Ifges(4,5) * t113 + Ifges(4,6) * t111 + Ifges(4,3) * t198 + t158 * t121 - t156 * t123 + t3;
t213 = cos(qJ(6));
t205 = sin(qJ(6));
t194 = pkin(16) + qJ(7) + qJ(2);
t193 = -qJ(8) + t201;
t192 = qJ(3) + qJ(4) + pkin(14);
t184 = cos(t194);
t183 = cos(t193);
t182 = sin(t194);
t181 = sin(t193);
t179 = cos(t192);
t178 = sin(t192);
t173 = qJDD(1) * t213 - t205 * t233;
t171 = qJDD(1) * t205 + t213 * t233;
t170 = pkin(6) * t219 + t229;
t168 = pkin(1) * t217 + pkin(2) * t184;
t167 = -pkin(1) * t209 - pkin(2) * t182;
t166 = qJDD(1) * pkin(6) - t231;
t164 = pkin(4) * t216 + pkin(9) * t179;
t163 = -pkin(4) * t208 - pkin(9) * t178;
t160 = -t183 * pkin(7) + pkin(3) * cos(t201);
t159 = -t181 * pkin(7) + pkin(3) * sin(t201);
t154 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t209 + Ifges(3,4) * t217) * qJD(1);
t153 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t205 + Ifges(7,4) * t213) * qJD(1);
t152 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t209 + Ifges(3,2) * t217) * qJD(1);
t151 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t205 + Ifges(7,2) * t213) * qJD(1);
t145 = -g(3) * t205 + t170 * t213;
t143 = -g(3) * t213 - t170 * t205;
t136 = 0.1e1 / (t178 * t183 + t179 * t181) / pkin(9) / pkin(7);
t131 = -mrSges(4,1) * t156 + mrSges(4,2) * t158;
t130 = -mrSges(8,1) * t155 + mrSges(8,2) * t157;
t119 = Ifges(4,5) * t158 + Ifges(4,6) * t156 + Ifges(4,3) * t200;
t118 = Ifges(8,5) * t157 + Ifges(8,6) * t155 + Ifges(8,3) * t199;
t86 = Ifges(5,5) * t127 + Ifges(5,6) * t126 + Ifges(5,3) * t191;
t71 = Ifges(9,5) * t101 + Ifges(9,6) * t100 + Ifges(9,3) * t190;
t28 = mrSges(9,2) * t58 - mrSges(9,3) * t43 + Ifges(9,1) * t57 + Ifges(9,4) * t56 + Ifges(9,5) * t188 + t100 * t71 - t190 * t72;
t27 = -mrSges(9,1) * t58 + mrSges(9,3) * t44 + Ifges(9,4) * t57 + Ifges(9,2) * t56 + Ifges(9,6) * t188 - t101 * t71 + t190 * t73;
t8 = mrSges(8,2) * t132 - mrSges(8,3) * t94 + Ifges(8,1) * t112 + Ifges(8,4) * t110 + Ifges(8,5) * t197 + t118 * t155 - t120 * t199 - t161 * t27 + t162 * t28 + t227 * t238;
t7 = -mrSges(8,1) * t132 + mrSges(8,3) * t96 + Ifges(8,4) * t112 + Ifges(8,2) * t110 + Ifges(8,6) * t197 - t118 * t157 + t122 * t199 + t161 * t28 + t162 * t27 + t227 * t237;
t5 = -pkin(8) * t12 - mrSges(5,1) * t81 + mrSges(5,3) * t50 + Ifges(5,4) * t68 + Ifges(5,2) * t67 + Ifges(5,6) * t189 - t127 * t86 + t191 * t88 - t224;
t4 = -pkin(10) * t12 + mrSges(5,2) * t81 - mrSges(5,3) * t49 + Ifges(5,1) * t68 + Ifges(5,4) * t67 + Ifges(5,5) * t189 + t126 * t86 - t15 * t206 + t16 * t214 - t191 * t87;
t2 = mrSges(4,2) * t132 - mrSges(4,3) * t95 + Ifges(4,1) * t113 + Ifges(4,4) * t111 + Ifges(4,5) * t198 + t119 * t156 - t121 * t200 - t207 * t5 + t215 * t4;
t1 = pkin(4) * t225 - mrSges(4,1) * t132 + mrSges(4,3) * t97 + Ifges(4,4) * t113 + Ifges(4,2) * t111 + Ifges(4,6) * t198 - t158 * t119 + t200 * t123 + t207 * t4 + t215 * t5;
t6 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t231 - mrSges(2,2) * t229 + t209 * (mrSges(3,2) * t165 - mrSges(3,3) * t144 + Ifges(3,1) * t172 + Ifges(3,4) * t174 + Ifges(3,5) * qJDD(2) - qJD(2) * t152 + t1 * t208 - t2 * t216 - t204 * t7 + t212 * t8) + t217 * (t244 * pkin(1) - mrSges(3,1) * t165 + mrSges(3,3) * t146 + Ifges(3,4) * t172 + Ifges(3,2) * t174 + Ifges(3,6) * qJDD(2) + qJD(2) * t154 - t216 * t1 - t208 * t2 + t204 * t8 + t212 * t7) + pkin(12) * (-m(3) * t165 + t174 * mrSges(3,1) - t172 * mrSges(3,2) + t244) + t205 * (mrSges(7,2) * t166 - mrSges(7,3) * t143 + Ifges(7,1) * t171 + Ifges(7,4) * t173 + Ifges(7,5) * qJDD(6) - qJD(6) * t151) + t213 * (-mrSges(7,1) * t166 + mrSges(7,3) * t145 + Ifges(7,4) * t171 + Ifges(7,2) * t173 + Ifges(7,6) * qJDD(6) + qJD(6) * t153) - pkin(6) * (-m(7) * t166 + t173 * mrSges(7,1) - t171 * mrSges(7,2)) + ((-pkin(6) * (t205 ^ 2 + t213 ^ 2) * mrSges(7,3) + pkin(12) * (t209 ^ 2 + t243) * mrSges(3,3)) * qJD(1) - pkin(6) * (-mrSges(7,1) * t205 - mrSges(7,2) * t213) * qJD(6) + pkin(12) * (-mrSges(3,1) * t209 - mrSges(3,2) * t217) * qJD(2)) * qJD(1); (-t208 * (m(4) * t97 - mrSges(4,2) * t198 + mrSges(4,3) * t111 + t10 * t215 + t131 * t156 - t142 * t200 - t18 * t207) - t216 * (m(4) * t95 + mrSges(4,1) * t198 - mrSges(4,3) * t113 - t131 * t158 + t140 * t200 + t240) + t204 * (m(8) * t96 - mrSges(8,2) * t197 + mrSges(8,3) * t110 + t130 * t155 - t141 * t199 + t235) + t212 * (m(8) * t94 + mrSges(8,1) * t197 - mrSges(8,3) * t112 - t130 * t157 + t139 * t199 + t236)) * pkin(1) + (t209 * t152 - t217 * t154) * qJD(1) + t222 + t220 + ((-t167 * t184 - t168 * t182) * pkin(2) * (mrSges(7,1) * t143 - mrSges(7,2) * t145 + Ifges(7,5) * t171 + Ifges(7,6) * t173 + Ifges(7,3) * qJDD(6) + (t205 * t151 - t213 * t153) * qJD(1)) + (-t222 + ((-t159 * t183 + t160 * t181) * t242 + (-t159 * t179 - t160 * t178) * t241) * t136) * (t167 * t213 + t168 * t205) * pkin(5)) / (-t182 * t213 + t184 * t205) / pkin(5) / pkin(2) + Ifges(3,5) * t172 + Ifges(3,6) * t174 + mrSges(3,1) * t144 - mrSges(3,2) * t146 + Ifges(3,3) * qJDD(2); ((t163 * t179 + t164 * t178) * t241 + (t163 * t183 - t164 * t181) * t242) * t136 + t220; t224;];
tau = t6(:);
