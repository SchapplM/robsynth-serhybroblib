% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% picker2Dm2IC
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
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = picker2Dm2IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2IC_invdynJ_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:21:03
% EndTime: 2020-05-11 09:21:05
% DurationCPUTime: 1.27s
% Computational Cost: add. (5422->221), mult. (6552->238), div. (37->10), fcn. (3746->64), ass. (0->143)
t152 = 0.1e1 / pkin(5);
t123 = qJDD(1) + qJDD(8);
t135 = sin(qJ(8));
t144 = cos(qJ(8));
t142 = sin(qJ(1));
t151 = cos(qJ(1));
t170 = -g(1) * t142 + t151 * g(2);
t67 = qJDD(1) * pkin(1) + t170;
t175 = t151 * g(1) + t142 * g(2);
t68 = -qJD(1) ^ 2 * pkin(1) + t175;
t54 = t135 * t68 - t144 * t67;
t56 = -t135 * t67 - t144 * t68;
t165 = mrSges(9,1) * t54 - mrSges(9,2) * t56 + Ifges(9,3) * t123;
t187 = t152 * t165;
t139 = sin(qJ(4));
t148 = cos(qJ(4));
t126 = qJD(1) + qJD(2);
t109 = qJD(4) + t126;
t103 = t109 ^ 2;
t124 = qJDD(1) + qJDD(2);
t106 = qJDD(4) + t124;
t130 = sin(qJ(10));
t132 = cos(qJ(10));
t141 = sin(qJ(2));
t150 = cos(qJ(2));
t55 = -t141 * t68 + t150 * t67;
t47 = pkin(3) * t124 + t55;
t122 = t126 ^ 2;
t57 = t141 * t67 + t150 * t68;
t49 = -pkin(3) * t122 + t57;
t35 = -t139 * t49 + t148 * t47;
t27 = pkin(4) * t106 + t35;
t37 = t139 * t47 + t148 * t49;
t29 = -pkin(4) * t103 + t37;
t22 = t130 * t29 - t132 * t27;
t81 = qJDD(10) + t106;
t94 = qJD(10) + t109;
t85 = t94 ^ 2;
t16 = m(11) * t22 + mrSges(11,1) * t81 - mrSges(11,2) * t85;
t23 = -t130 * t27 - t132 * t29;
t17 = m(11) * t23 - mrSges(11,1) * t85 - mrSges(11,2) * t81;
t162 = -t130 * t17 - t132 * t16;
t8 = m(5) * t35 + mrSges(5,1) * t106 - mrSges(5,2) * t103 + t162;
t9 = m(5) * t37 - mrSges(5,1) * t103 - mrSges(5,2) * t106 + t130 * t16 - t132 * t17;
t186 = t139 * t9 + t148 * t8;
t185 = t151 * pkin(1);
t176 = qJ(3) - qJ(6);
t177 = qJ(2) + qJ(4);
t116 = qJ(1) + t177;
t99 = -qJ(9) + t116;
t77 = t99 - t176;
t98 = qJ(9) + t116;
t78 = t98 + t176;
t184 = sin(t78) - sin(t77);
t183 = cos(t78) - cos(t77);
t181 = sin(t99) - sin(t98);
t180 = cos(t99) - cos(t98);
t128 = qJ(1) + qJ(8);
t174 = pkin(8) + qJ(5);
t166 = t174 + t176;
t155 = t166 - t128;
t167 = t174 - t176;
t156 = t167 - t128;
t179 = t152 / (pkin(6) * (cos(t156) - cos(t155)) + (-cos(qJ(9) - t156) + cos(qJ(9) + t155)) * pkin(2));
t154 = 0.1e1 / pkin(3);
t172 = -qJ(9) - t176;
t53 = 0.1e1 / ((-cos(qJ(4) - t176) + cos(qJ(4) + t176)) * pkin(6) + (cos(qJ(4) + t172) - cos(qJ(4) - t172)) * pkin(2));
t178 = t154 * t53;
t127 = 0.1e1 / t139;
t153 = 0.1e1 / pkin(4);
t173 = pkin(1) * t127 * t153;
t171 = -qJ(1) + t174;
t140 = sin(qJ(3));
t149 = cos(qJ(3));
t48 = pkin(2) * t124 + t55;
t50 = -pkin(2) * t122 + t57;
t36 = t140 * t50 - t149 * t48;
t110 = qJD(3) + t126;
t107 = qJDD(3) + t124;
t105 = qJDD(6) + t124;
t137 = sin(qJ(6));
t146 = cos(qJ(6));
t40 = t137 * t57 - t146 * t55;
t41 = -t137 * t55 - t146 * t57;
t26 = mrSges(7,1) * t40 - mrSges(7,2) * t41 + Ifges(7,3) * t105;
t134 = sin(qJ(9));
t143 = cos(qJ(9));
t28 = pkin(6) * t107 + t36;
t104 = t110 ^ 2;
t38 = -t140 * t48 - t149 * t50;
t30 = -pkin(6) * t104 + t38;
t24 = t134 * t30 - t143 * t28;
t25 = -t134 * t28 - t143 * t30;
t92 = qJDD(9) + t107;
t15 = mrSges(10,1) * t24 - mrSges(10,2) * t25 + Ifges(10,3) * t92;
t169 = -qJ(4) + t171;
t168 = qJ(4) + t171;
t14 = mrSges(11,1) * t22 - mrSges(11,2) * t23 + Ifges(11,3) * t81;
t95 = qJD(9) + t110;
t93 = t95 ^ 2;
t18 = m(10) * t24 + mrSges(10,1) * t92 - mrSges(10,2) * t93;
t19 = m(10) * t25 - mrSges(10,1) * t93 - mrSges(10,2) * t92;
t161 = -t134 * t19 - t143 * t18;
t10 = m(4) * t36 + mrSges(4,1) * t107 - mrSges(4,2) * t104 + t161;
t11 = m(4) * t38 - mrSges(4,1) * t104 - mrSges(4,2) * t107 + t134 * t18 - t143 * t19;
t164 = -t10 * t149 - t11 * t140;
t129 = qJ(1) + qJ(2);
t112 = sin(t129);
t113 = cos(t129);
t96 = sin(t116);
t97 = cos(t116);
t163 = -t112 * t97 + t113 * t96;
t136 = sin(qJ(7));
t145 = cos(qJ(7));
t160 = t136 * t96 + t145 * t97;
t159 = -qJ(8) + t169;
t158 = -qJ(8) + t168;
t5 = pkin(6) * t161 + mrSges(4,1) * t36 - mrSges(4,2) * t38 + Ifges(4,3) * t107 + t15;
t4 = pkin(4) * t162 + mrSges(5,1) * t35 - mrSges(5,2) * t37 + Ifges(5,3) * t106 + t14;
t1 = pkin(2) * t164 + pkin(3) * t186 + mrSges(3,1) * t55 - mrSges(3,2) * t57 + Ifges(3,3) * t124 + t26 + t4 + t5;
t147 = cos(qJ(5));
t138 = sin(qJ(5));
t133 = cos(pkin(8));
t131 = sin(pkin(8));
t125 = qJD(1) + qJD(8);
t121 = t125 ^ 2;
t118 = t142 * pkin(1);
t115 = -qJ(9) + t174;
t114 = qJ(9) + t174;
t111 = sin(t177);
t108 = qJD(6) + t126;
t102 = t108 ^ 2;
t91 = qJ(9) + t166;
t90 = -qJ(9) + t167;
t70 = -pkin(5) * sin(t128) + t118;
t69 = t185 - pkin(5) * cos(t128);
t66 = -t131 * t138 + t133 * t147;
t65 = t131 * t147 + t133 * t138;
t61 = pkin(3) * t113 + pkin(4) * t97 + t185;
t60 = pkin(3) * t112 + pkin(4) * t96 + t118;
t32 = m(7) * t41 - mrSges(7,1) * t102 - mrSges(7,2) * t105;
t31 = m(7) * t40 + mrSges(7,1) * t105 - mrSges(7,2) * t102;
t2 = [-t141 * t14 * t173 - mrSges(2,2) * t175 - pkin(5) * t187 + t135 / sin(-qJ(8) + t171) * (Ifges(6,3) * qJDD(5) + mrSges(6,1) * (g(1) * t65 - g(2) * t66) - mrSges(6,2) * (-g(1) * t66 - g(2) * t65)) + t1 + mrSges(2,1) * t170 + Ifges(2,3) * qJDD(1) + t165 + (((-t183 * t61 - t184 * t60) * t178 + ((sin(t91) - sin(t90)) * t70 + (cos(t91) - cos(t90)) * t69) * t179) * t5 + ((t180 * t61 + t181 * t60) * t178 + (-(sin(t115) - sin(t114)) * t70 - (cos(t115) - cos(t114)) * t69) * t179) * t26) * pkin(2) + ((-cos(qJ(8) + 0.2e1 * qJ(1)) + cos(0.2e1 * pkin(8) + 0.2e1 * qJ(5) + qJ(8))) / (-cos(0.2e1 * t128) + cos(0.2e1 * t174)) * t187 - t135 * (m(9) * t56 - mrSges(9,1) * t121 - mrSges(9,2) * t123) - t144 * (m(9) * t54 + mrSges(9,1) * t123 - mrSges(9,2) * t121) + t141 * (m(3) * t57 - mrSges(3,1) * t122 - mrSges(3,2) * t124 + t10 * t140 - t11 * t149 + t137 * t31 - t139 * t8 - t146 * t32 + t148 * t9) + t150 * (m(3) * t55 + mrSges(3,1) * t124 - mrSges(3,2) * t122 - t137 * t32 - t146 * t31 + t164 + t186)) * pkin(1) + ((-pkin(1) * t111 - pkin(3) * t139) * t127 * t1 + (pkin(3) * t141 + pkin(4) * t111) * t4 * t173 - (pkin(3) * (cos(t169) - cos(t168)) + (cos(-qJ(2) + t159) - cos(qJ(2) + t158)) * pkin(5)) * pkin(1) * t152 / (cos(t159) - cos(t158)) * t15) * t154; Ifges(8,3) * qJDD(7) + mrSges(8,1) * (-g(1) * t145 - g(2) * t136) - mrSges(8,2) * (-g(1) * t136 + g(2) * t145) + ((-t181 * t26 + t184 * t5) * t145 + (t180 * t26 - t183 * t5) * t136) * t53 * pkin(2) + (-(t1 + t15) * t160 + ((t14 * t163 + t160 * t4) * pkin(4) + (-t14 + t4) * pkin(3) * (t112 * t136 + t113 * t145)) * t153) / t163;];
tau = t2(:);
