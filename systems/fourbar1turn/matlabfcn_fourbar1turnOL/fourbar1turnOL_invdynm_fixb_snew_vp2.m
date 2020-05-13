% Calculate vector of cutting torques with Newton-Euler for
% fourbar1turnOL
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = fourbar1turnOL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:50
% EndTime: 2020-04-12 19:40:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (2835->204), mult. (5914->277), div. (0->0), fcn. (3624->8), ass. (0->81)
t140 = sin(qJ(4));
t163 = qJD(1) * t140;
t142 = sin(qJ(2));
t162 = qJD(1) * t142;
t144 = cos(qJ(4));
t161 = qJD(1) * t144;
t146 = cos(qJ(2));
t160 = qJD(1) * t146;
t159 = qJD(1) * qJD(2);
t158 = qJD(1) * qJD(4);
t157 = t142 * t159;
t143 = sin(qJ(1));
t147 = cos(qJ(1));
t133 = t143 * g(1) - t147 * g(2);
t134 = -t147 * g(1) - t143 * g(2);
t108 = -t142 * g(3) + t146 * t134;
t107 = -t146 * g(3) - t142 * t134;
t141 = sin(qJ(3));
t145 = cos(qJ(3));
t115 = (t142 * t141 - t146 * t145) * qJD(1);
t138 = qJD(2) + qJD(3);
t103 = -t138 * mrSges(4,2) + t115 * mrSges(4,3);
t116 = (-t146 * t141 - t142 * t145) * qJD(1);
t137 = qJDD(2) + qJDD(3);
t148 = qJD(1) ^ 2;
t101 = (t142 * t146 * t148 + qJDD(2)) * pkin(2) + t107;
t102 = (-t146 ^ 2 * t148 - qJD(2) ^ 2) * pkin(2) + t108;
t89 = -t145 * t101 + t141 * t102;
t125 = t142 * qJDD(1) + t146 * t159;
t127 = t146 * qJDD(1) - t157;
t95 = t115 * qJD(3) - t145 * t125 - t141 * t127;
t99 = -t115 * mrSges(4,1) + t116 * mrSges(4,2);
t79 = m(4) * t89 + t137 * mrSges(4,1) - t95 * mrSges(4,3) + t138 * t103 - t116 * t99;
t104 = t138 * mrSges(4,1) - t116 * mrSges(4,3);
t90 = -t141 * t101 - t145 * t102;
t94 = -t116 * qJD(3) + t141 * t125 - t145 * t127;
t80 = m(4) * t90 - t137 * mrSges(4,2) + t94 * mrSges(4,3) - t138 * t104 + t115 * t99;
t75 = -t141 * t80 - t145 * t79;
t111 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t140 + Ifges(5,2) * t144) * qJD(1);
t113 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t140 + Ifges(5,4) * t144) * qJD(1);
t156 = t140 * t111 - t144 * t113;
t112 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t142 + Ifges(3,2) * t146) * qJD(1);
t114 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t142 + Ifges(3,4) * t146) * qJD(1);
t155 = t142 * t112 - t146 * t114;
t120 = -qJDD(1) * pkin(1) - t133;
t124 = t140 * qJDD(1) + t144 * t158;
t126 = t144 * qJDD(1) - t140 * t158;
t131 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t161;
t154 = -m(5) * t120 + t126 * mrSges(5,1) - t124 * mrSges(5,2) + t131 * t161;
t129 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t163;
t110 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t142 + Ifges(3,6) * t146) * qJD(1);
t100 = (-t127 + t157) * pkin(2) - t133;
t151 = m(4) * t100 - t94 * mrSges(4,1) + t95 * mrSges(4,2) - t115 * t103 + t116 * t104;
t96 = Ifges(4,5) * t116 + Ifges(4,6) * t115 + Ifges(4,3) * t138;
t98 = Ifges(4,1) * t116 + Ifges(4,4) * t115 + Ifges(4,5) * t138;
t77 = -mrSges(4,1) * t100 + mrSges(4,3) * t90 + Ifges(4,4) * t95 + Ifges(4,2) * t94 + Ifges(4,6) * t137 - t116 * t96 + t138 * t98;
t97 = Ifges(4,4) * t116 + Ifges(4,2) * t115 + Ifges(4,6) * t138;
t78 = mrSges(4,2) * t100 - mrSges(4,3) * t89 + Ifges(4,1) * t95 + Ifges(4,4) * t94 + Ifges(4,5) * t137 + t115 * t96 - t138 * t97;
t72 = mrSges(3,1) * t133 + mrSges(3,3) * t108 + Ifges(3,4) * t125 + Ifges(3,2) * t127 + Ifges(3,6) * qJDD(2) - pkin(2) * t151 + qJD(2) * t114 - t110 * t162 - t141 * t78 - t145 * t77;
t74 = -mrSges(3,2) * t133 - mrSges(3,3) * t107 + Ifges(3,1) * t125 + Ifges(3,4) * t127 + Ifges(3,5) * qJDD(2) - qJD(2) * t112 + t110 * t160 + t141 * t77 - t145 * t78;
t123 = -t148 * pkin(1) + t134;
t106 = -t140 * g(3) + t144 * t123;
t109 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t140 + Ifges(5,6) * t144) * qJD(1);
t83 = -mrSges(5,1) * t120 + mrSges(5,3) * t106 + Ifges(5,4) * t124 + Ifges(5,2) * t126 + Ifges(5,6) * qJDD(4) + qJD(4) * t113 - t109 * t163;
t105 = -t144 * g(3) - t140 * t123;
t84 = mrSges(5,2) * t120 - mrSges(5,3) * t105 + Ifges(5,1) * t124 + Ifges(5,4) * t126 + Ifges(5,5) * qJDD(4) - qJD(4) * t111 + t109 * t161;
t153 = -mrSges(2,2) * t134 + t140 * t84 + t142 * t74 + t144 * t83 + t146 * t72 + mrSges(2,1) * t133 + Ifges(2,3) * qJDD(1) + pkin(1) * (-t129 * t163 + t154);
t152 = -mrSges(4,1) * t89 + mrSges(4,2) * t90 - Ifges(4,5) * t95 - Ifges(4,6) * t94 - Ifges(4,3) * t137 + t115 * t98 - t116 * t97;
t150 = mrSges(5,1) * t105 - mrSges(5,2) * t106 + Ifges(5,5) * t124 + Ifges(5,6) * t126 + Ifges(5,3) * qJDD(4);
t149 = mrSges(3,1) * t107 - mrSges(3,2) * t108 + Ifges(3,5) * t125 + Ifges(3,6) * t127 + Ifges(3,3) * qJDD(2) + pkin(2) * t75 - t152;
t132 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t160;
t130 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t162;
t122 = (-mrSges(3,1) * t146 + mrSges(3,2) * t142) * qJD(1);
t121 = (-mrSges(5,1) * t144 + mrSges(5,2) * t140) * qJD(1);
t87 = m(5) * t106 - qJDD(4) * mrSges(5,2) + t126 * mrSges(5,3) - qJD(4) * t129 + t121 * t161;
t86 = m(5) * t105 + qJDD(4) * mrSges(5,1) - t124 * mrSges(5,3) + qJD(4) * t131 - t121 * t163;
t76 = qJDD(1) * mrSges(2,1) + t127 * mrSges(3,1) - t148 * mrSges(2,2) - t125 * mrSges(3,2) + (m(2) + m(3)) * t133 + (-t140 * t129 - t142 * t130 + t146 * t132) * qJD(1) - t151 + t154;
t70 = m(2) * t134 - qJDD(1) * mrSges(2,2) - t148 * mrSges(2,1) + t146 * (m(3) * t108 - qJDD(2) * mrSges(3,2) + t127 * mrSges(3,3) - qJD(2) * t130 + t122 * t160 + t141 * t79 - t145 * t80) - t142 * (m(3) * t107 + qJDD(2) * mrSges(3,1) - t125 * mrSges(3,3) + qJD(2) * t132 - t122 * t162 + t75) + t144 * t87 - t140 * t86;
t69 = -t149 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-t155 - t156) * qJD(1) + mrSges(2,3) * t134 - pkin(1) * (t140 * t87 + t144 * t86) + t148 * Ifges(2,5) - t150;
t68 = -mrSges(2,2) * g(3) - mrSges(2,3) * t133 + Ifges(2,5) * qJDD(1) - t148 * Ifges(2,6) - t140 * t83 - t142 * t72 + t144 * t84 + t146 * t74;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t147 * t68 - t143 * t69 - pkin(5) * (t143 * t70 + t147 * t76), t68, t74, t78, t84, 0, 0; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t143 * t68 + t147 * t69 + pkin(5) * (-t143 * t76 + t147 * t70), t69, t72, t77, t83, 0, 0; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t153, t153, t155 * qJD(1) + t149, -t152, t156 * qJD(1) + t150, 0, 0;];
m_new = t1;
