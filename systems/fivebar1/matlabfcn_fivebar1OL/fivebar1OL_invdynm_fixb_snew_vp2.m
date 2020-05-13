% Calculate vector of cutting torques with Newton-Euler for
% fivebar1OL
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
%   pkin=[AB,AE,BC,CD,ED]';
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
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = fivebar1OL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:11
% EndTime: 2020-04-27 06:13:12
% DurationCPUTime: 0.55s
% Computational Cost: add. (419->146), mult. (631->180), div. (0->0), fcn. (364->8), ass. (0->77)
t117 = sin(qJ(3));
t121 = cos(qJ(3));
t116 = sin(qJ(4));
t105 = t116 * mrSges(5,2);
t120 = cos(qJ(4));
t162 = t120 * mrSges(5,1);
t94 = -t105 + t162;
t80 = pkin(3) * m(5) + mrSges(4,1) + t94;
t93 = t116 * mrSges(5,1) + t120 * mrSges(5,2);
t88 = mrSges(4,2) + t93;
t149 = t117 * t80 + t88 * t121;
t123 = cos(qJ(1));
t118 = sin(qJ(2));
t122 = cos(qJ(2));
t96 = t118 * mrSges(3,1) + t122 * mrSges(3,2);
t171 = -mrSges(2,2) + t96;
t77 = t171 * t123;
t179 = -t149 + t77 - mrSges(1,2);
t119 = sin(qJ(1));
t128 = pkin(2) * m(3);
t107 = t118 * mrSges(3,2);
t161 = t122 * mrSges(3,1);
t175 = -t161 + t107;
t168 = -mrSges(2,1) - t128 - t175;
t178 = t119 * t168;
t177 = t119 * t171;
t176 = (2 * qJD(1) + qJD(2)) * qJD(2);
t166 = (-2 * qJD(3) - qJD(4)) * qJD(4);
t169 = -t88 * t117 + t80 * t121;
t72 = t168 * t123;
t174 = t169 + mrSges(1,1) - t72;
t115 = mrSges(3,3) + mrSges(2,3);
t114 = mrSges(4,3) + mrSges(5,3);
t103 = pkin(3) * t162;
t157 = Ifges(5,3) + t103;
t154 = pkin(3) * t105;
t147 = (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t107 - 0.2e1 * t161 + t128) * pkin(2)) * qJDD(1) + (t175 * pkin(2) + Ifges(3,3)) * qJDD(2) + t176 * pkin(2) * t96;
t95 = t122 * Ifges(3,5) - t118 * Ifges(3,6);
t92 = t120 * Ifges(5,5) - t116 * Ifges(5,6);
t142 = -t121 * g(1) - t117 * g(2);
t141 = -t117 * g(1) + t121 * g(2);
t140 = -t123 * g(1) - t119 * g(2);
t139 = -t119 * g(1) + t123 * g(2);
t97 = Ifges(3,5) * t118 + Ifges(3,6) * t122;
t91 = t116 * Ifges(5,5) + t120 * Ifges(5,6);
t78 = -pkin(3) * mrSges(5,3) + Ifges(4,5) + t92;
t87 = Ifges(4,6) + t91;
t138 = t117 * t78 + t87 * t121;
t136 = t117 * t92 + t91 * t121;
t68 = -t117 * t91 + t92 * t121;
t79 = -pkin(2) * mrSges(3,3) + Ifges(2,5) - t95;
t90 = -Ifges(2,6) + t97;
t135 = t119 * t79 - t90 * t123;
t69 = t119 * t97 - t95 * t123;
t134 = t119 * t95 + t97 * t123;
t133 = (pkin(3) ^ 2 * m(5) + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t103 - 0.2e1 * t154) * qJDD(3);
t132 = qJD(1) ^ 2;
t130 = qJD(3) ^ 2;
t126 = m(4) + m(5);
t113 = qJD(1) + qJD(2);
t112 = qJD(3) + qJD(4);
t111 = qJDD(1) + qJDD(2);
t110 = qJDD(3) + qJDD(4);
t109 = t113 ^ 2;
t108 = t112 ^ 2;
t100 = mrSges(1,3) + t114 + t115;
t86 = -t132 * pkin(2) + t140;
t85 = -t130 * pkin(3) + t142;
t84 = qJDD(1) * pkin(2) - t139;
t83 = qJDD(3) * pkin(3) - t141;
t67 = t119 * t90 + t79 * t123;
t66 = -t117 * t87 + t78 * t121;
t65 = -t118 * t84 - t122 * t86;
t64 = t116 * t83 + t120 * t85;
t63 = t118 * t86 - t122 * t84;
t62 = -t116 * t85 + t120 * t83;
t1 = [-t138 * t130 + t66 * qJDD(3) + t68 * qJDD(4) - (-t178 - t179) * g(3) + t69 * qJDD(2) - t135 * t132 + t67 * qJDD(1) + t100 * g(2) + t166 * t136 + t176 * t134, g(3) * t171 + t79 * qJDD(1) - t95 * qJDD(2) + t139 * t115 + t90 * t132 + t176 * t97, -mrSges(3,2) * g(3) - mrSges(3,3) * t63 + Ifges(3,5) * t111 - t109 * Ifges(3,6), -t88 * g(3) + t78 * qJDD(3) + t92 * qJDD(4) + t141 * t114 - t87 * t130 + t166 * t91, -mrSges(5,2) * g(3) - mrSges(5,3) * t62 + Ifges(5,5) * t110 - t108 * Ifges(5,6), 0, 0; t136 * qJDD(4) + t66 * t130 + t138 * qJDD(3) - (-pkin(1) * t126 - t174 - t177) * g(3) - t100 * g(1) - t134 * qJDD(2) + t67 * t132 + t135 * qJDD(1) + t176 * t69 - t166 * t68, -g(3) * t168 - t90 * qJDD(1) - t97 * qJDD(2) + t140 * t115 + t79 * t132 - t176 * t95, mrSges(3,1) * g(3) + mrSges(3,3) * t65 + t109 * Ifges(3,5) + Ifges(3,6) * t111, t80 * g(3) + t87 * qJDD(3) + t91 * qJDD(4) + t142 * t114 + t78 * t130 - t166 * t92, mrSges(5,1) * g(3) + mrSges(5,3) * t64 + t108 * Ifges(5,5) + Ifges(5,6) * t110, 0, 0; (-t154 + t157) * qJDD(4) + t133 - t179 * g(1) - t174 * g(2) - t166 * (-pkin(3) * t93 + (-t117 * t94 - t121 * t93) * pkin(1)) + (-g(1) * t168 - g(2) * t171) * t119 + ((-t117 * t93 + t121 * t94) * qJDD(4) + t169 * qJDD(3) - t149 * t130 - t126 * g(2)) * pkin(1) + t147, -(t77 + t178) * g(1) - (-t72 + t177) * g(2) + t147, mrSges(3,1) * t63 - mrSges(3,2) * t65 + Ifges(3,3) * t111, t149 * g(1) - t169 * g(2) + t133 + t157 * qJDD(4) + (-qJDD(4) * t105 + t166 * t93) * pkin(3), mrSges(5,1) * t62 - mrSges(5,2) * t64 + Ifges(5,3) * t110, 0, 0;];
m_new = t1;
