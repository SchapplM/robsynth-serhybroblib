% Calculate vector of inverse dynamics joint torques for
% palh2m2DE
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m2DE_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:32
% DurationCPUTime: 2.07s
% Computational Cost: add. (875->237), mult. (1862->355), div. (0->0), fcn. (909->8), ass. (0->120)
t74 = cos(qJ(3));
t146 = pkin(5) * t74;
t71 = sin(qJ(2));
t148 = pkin(4) * t71;
t75 = cos(qJ(2));
t147 = pkin(4) * t75;
t63 = qJDD(2) * t147;
t70 = sin(qJ(3));
t78 = qJD(3) ^ 2;
t79 = qJD(2) ^ 2;
t16 = -t79 * t148 + t63 + (qJDD(3) * t74 - t70 * t78) * pkin(5);
t165 = t16 * t146;
t164 = -t71 / 0.2e1;
t125 = pkin(4) * qJD(2);
t34 = t70 * t75 - t74 * t71;
t124 = qJD(1) * t70;
t44 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t124;
t121 = qJD(1) * t74;
t45 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t121;
t91 = t70 * t71 + t74 * t75;
t163 = (t34 * t44 - t45 * t91) * t125;
t53 = -pkin(1) - t147;
t105 = m(4) * t53 - mrSges(4,1);
t49 = -pkin(2) + t53;
t162 = m(5) * t49 + t105;
t160 = (-t34 * t70 - t74 * t91) * mrSges(5,3);
t159 = -mrSges(6,3) - mrSges(4,3);
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t158 = t73 * mrSges(7,1) - t69 * mrSges(7,2);
t156 = mrSges(7,1) * t69 + mrSges(7,2) * t73;
t155 = qJD(2) - qJD(3);
t141 = Ifges(3,4) * t71;
t154 = pkin(1) * (mrSges(3,1) * t71 + mrSges(3,2) * t75) + (Ifges(3,1) * t75 - t141) * t164;
t37 = t49 - t146;
t104 = m(6) * t37 - mrSges(6,1);
t68 = qJD(1) + qJD(4);
t153 = t156 * t68;
t80 = qJD(1) ^ 2;
t152 = t70 / 0.2e1;
t151 = t74 / 0.2e1;
t150 = t75 / 0.2e1;
t149 = m(6) + m(7);
t128 = t75 * t79;
t36 = (-qJDD(2) * t71 - t128) * pkin(4);
t15 = (-qJDD(3) * t70 - t74 * t78) * pkin(5) + t36;
t109 = t71 * t125;
t38 = -pkin(5) * qJD(3) * t70 - t109;
t127 = -t38 * qJD(1) + t37 * qJDD(1);
t8 = -qJDD(1) * pkin(3) + t127;
t116 = -pkin(3) + t37;
t19 = t116 * qJD(1);
t94 = t69 * t19 - t38 * t73;
t3 = t94 * qJD(4) - t15 * t69 - t73 * t8;
t67 = qJDD(1) + qJDD(4);
t144 = t3 * mrSges(7,1) + Ifges(7,3) * t67;
t140 = Ifges(5,4) * t70;
t139 = Ifges(3,5) * t75;
t138 = Ifges(5,5) * t74;
t137 = Ifges(3,6) * t71;
t136 = Ifges(5,6) * t70;
t133 = t68 * mrSges(7,1);
t131 = t71 * t80;
t129 = t74 * Ifges(5,2);
t54 = mrSges(7,1) * g(1) + mrSges(7,2) * g(2);
t55 = mrSges(7,1) * g(2) - mrSges(7,2) * g(1);
t126 = t54 * t69 - t55 * t73;
t123 = qJD(1) * t71;
t122 = qJD(1) * t73;
t120 = qJD(1) * t75;
t119 = qJD(2) * t75;
t118 = qJD(3) * t74;
t117 = qJD(4) * t68;
t115 = qJDD(1) * t71;
t114 = m(5) + t149;
t113 = qJD(1) * qJD(2);
t112 = qJD(1) * qJD(3);
t108 = pkin(4) * t119;
t106 = m(4) + t114;
t103 = t71 * t113;
t99 = t149 * pkin(5) + mrSges(5,1);
t72 = sin(qJ(1));
t76 = cos(qJ(1));
t98 = g(1) * t76 + g(2) * t72;
t97 = g(1) * t72 - g(2) * t76;
t6 = -t19 * t73 - t38 * t69;
t96 = t6 * t73 - t69 * t94;
t95 = t6 * t69 + t73 * t94;
t92 = -t54 * t73 - t55 * t69;
t50 = t106 * pkin(4) + mrSges(3,1);
t61 = Ifges(5,4) * t121;
t89 = -t70 * (Ifges(5,6) * qJD(3) + (t129 + t140) * qJD(1)) / 0.2e1 + (Ifges(5,1) * t124 + Ifges(5,5) * qJD(3) + t61) * t151;
t88 = t73 * t117 + t69 * t67;
t87 = t69 * t117 - t73 * t67;
t2 = t6 * qJD(4) + t15 * t73 - t69 * t8;
t83 = t96 * qJD(4) - t2 * t73 + t3 * t69;
t82 = t49 * (mrSges(5,1) * t70 + mrSges(5,2) * t74) + (Ifges(5,1) * t74 - t140) * t152;
t62 = Ifges(3,4) * t120;
t56 = mrSges(2,2) - mrSges(3,3) - mrSges(5,3) + t159;
t51 = pkin(4) * t103;
t43 = t75 * t113 + t115;
t42 = qJDD(1) * t75 - t103;
t41 = qJDD(1) * t70 + t74 * t112;
t40 = qJDD(1) * t74 - t70 * t112;
t39 = pkin(5) * t118 + t108;
t35 = (-mrSges(5,1) * t74 + mrSges(5,2) * t70) * qJD(1);
t29 = Ifges(3,1) * t123 + Ifges(3,5) * qJD(2) + t62;
t27 = Ifges(3,6) * qJD(2) + (t75 * Ifges(3,2) + t141) * qJD(1);
t24 = (-t73 * t119 - t69 * t123) * pkin(4);
t23 = (-t73 * t118 - t69 * t124) * pkin(5);
t22 = (t69 * t119 - t71 * t122) * pkin(4);
t21 = (t69 * t118 - t70 * t122) * pkin(5);
t20 = qJDD(1) * t49 + t51;
t17 = (m(3) + t106) * pkin(1) + pkin(3) * m(7) + mrSges(2,1) + mrSges(4,1) + mrSges(6,1) + t114 * pkin(2);
t14 = t16 * t147;
t10 = t155 * t91;
t9 = t155 * t34;
t5 = (-qJD(2) * t10 + qJDD(2) * t34) * pkin(4);
t4 = (qJD(2) * t9 + qJDD(2) * t91) * pkin(4);
t1 = [(g(1) * t17 + g(2) * t56 - t92) * t72 + (g(1) * t56 - g(2) * t17 + t126) * t76 - pkin(1) * (-mrSges(3,1) * t42 + mrSges(3,2) * t43) + t36 * mrSges(4,3) + t15 * mrSges(6,3) - t2 * mrSges(7,2) + (m(5) * t20 - mrSges(5,1) * t40 + mrSges(5,2) * t41) * t49 + t105 * (t53 * qJDD(1) + t51) + t104 * t127 + (m(7) * t96 - t104 * qJD(1) + t158 * t68) * t38 + (m(7) * (t95 * qJD(4) - t2 * t69 - t3 * t73) + t88 * mrSges(7,2) + t87 * mrSges(7,1)) * t116 + (Ifges(3,4) * t43 + Ifges(3,2) * t42 + Ifges(3,6) * qJDD(2) + t97 * t50) * t75 + (-t97 * mrSges(3,2) + Ifges(3,1) * t43 + Ifges(3,4) * t42 + Ifges(3,5) * qJDD(2)) * t71 + (-t53 * mrSges(4,1) - t37 * mrSges(6,1) + Ifges(4,2) + Ifges(6,2) + Ifges(2,3) + (m(3) * pkin(1) + mrSges(3,1) * t75 - mrSges(3,2) * t71) * pkin(1)) * qJDD(1) + (-t20 * mrSges(5,1) + t5 * mrSges(5,3) + Ifges(5,4) * t41 + Ifges(5,2) * t40 + Ifges(5,6) * qJDD(3) + t97 * t99) * t74 + (-t4 * mrSges(5,3) + Ifges(5,1) * t41 + Ifges(5,4) * t40 + Ifges(5,5) * qJDD(3) + (t20 - t97) * mrSges(5,2)) * t70 + ((t138 / 0.2e1 - t136 / 0.2e1) * qJD(3) + ((t74 * Ifges(5,4) - Ifges(5,2) * t70) * t151 + t82) * qJD(1) + t89) * qJD(3) + (t29 * t150 + t27 * t164 + (t139 / 0.2e1 - t137 / 0.2e1) * qJD(2) + ((Ifges(3,4) * t75 - Ifges(3,2) * t71) * t150 - t154) * qJD(1) + (qJD(3) * t160 + (t162 * qJD(1) + t35) * t71) * pkin(4)) * qJD(2) + t144; t27 * t123 / 0.2e1 - (-t137 + t139) * t113 / 0.2e1 - t22 * t133 + m(6) * t14 + (mrSges(3,2) * g(3) + t98 * t50) * t71 - (-t98 * mrSges(3,2) + g(3) * t50) * t75 + t24 * t68 * mrSges(7,2) + Ifges(3,6) * t42 + Ifges(3,5) * t43 + Ifges(3,3) * qJDD(2) + t154 * t80 - t163 - (-Ifges(3,2) * t123 + t29 + t62) * t120 / 0.2e1 + t108 * t153 + (t39 * t109 - t22 * t6 + t24 * t94 + t14) * m(7) + (t158 * t117 + t156 * t67) * t148 + (t159 * t115 + (mrSges(6,1) - t162) * t131 + m(5) * (t34 * t5 + t4 * t91) + m(7) * (t95 * t119 + (-qJD(2) * t39 + t83) * t71) + t91 * (qJDD(3) * mrSges(5,1) - mrSges(5,3) * t41) + t9 * t44 - t10 * t45 + t34 * (-qJDD(3) * mrSges(5,2) + mrSges(5,3) * t40) - t35 * t123 + m(4) * (-t36 * t71 + t63 * t75) + (-t37 * t131 - t15 * t71) * m(6) + ((-t10 * t34 + t9 * t91) * qJD(2) * m(5) - t128 * t71 * m(4)) * pkin(4)) * pkin(4); t4 * mrSges(5,1) + Ifges(5,5) * t41 + Ifges(5,6) * t40 + Ifges(5,3) * qJDD(3) + (-mrSges(7,1) * t21 + mrSges(7,2) * t23) * t68 + (-t74 * g(3) + t98 * t70) * t99 + (g(3) * t70 + t74 * t98 - t5) * mrSges(5,2) + t163 + m(6) * t165 + (-t74 * t61 / 0.2e1 - qJD(3) * (-t136 + t138) / 0.2e1 + (t129 * t152 - t82) * qJD(1) - t125 * t160 - t89) * qJD(1) + (t153 * t118 + (-m(6) * t15 + t88 * mrSges(7,1) - t87 * mrSges(7,2) - qJDD(1) * mrSges(6,3) - t104 * t80) * t70) * pkin(5) + (-t6 * t21 + t23 * t94 + t165 + (t95 * t118 + t83 * t70) * pkin(5)) * m(7); -t94 * t133 + t126 * t76 - t92 * t72 + (t6 * t68 - t2) * mrSges(7,2) + t144;];
tau = t1;
