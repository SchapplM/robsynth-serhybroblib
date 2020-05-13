% Calculate matrix of centrifugal and coriolis load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh2m1DE_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:07
% EndTime: 2020-05-02 23:52:09
% DurationCPUTime: 0.52s
% Computational Cost: add. (985->141), mult. (1264->203), div. (0->0), fcn. (473->22), ass. (0->101)
t90 = m(5) + m(6);
t65 = t90 * pkin(3) + mrSges(4,1);
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t25 = (t87 * mrSges(4,2) + t65 * t84) * pkin(2);
t79 = qJ(4) + qJ(2);
t76 = qJ(3) + t79;
t60 = sin(t76);
t138 = t60 / 0.2e1;
t80 = -qJ(4) + qJ(2);
t77 = qJ(3) + t80;
t61 = sin(t77);
t62 = cos(t76);
t63 = cos(t77);
t82 = qJ(2) + qJ(3);
t70 = sin(t82);
t74 = cos(t82);
t92 = pkin(3) ^ 2;
t12 = (t61 * t62 / 0.2e1 + t63 * t138 - t74 * t70) * m(6) * t92;
t72 = cos(t80);
t131 = pkin(2) * t72;
t71 = cos(t79);
t132 = pkin(2) * t71;
t54 = -t132 / 0.2e1;
t145 = t54 + t131 / 0.2e1;
t55 = mrSges(6,2) * t62;
t114 = pkin(3) * t55;
t129 = pkin(3) * t63;
t109 = t129 / 0.2e1;
t130 = pkin(3) * t61;
t40 = -t130 / 0.2e1;
t123 = mrSges(6,1) * t40 + mrSges(6,2) * t109;
t39 = pkin(3) * t138;
t106 = t114 / 0.2e1 + mrSges(6,1) * t39 + t123;
t105 = Ifges(4,6) * t70 + t106;
t135 = mrSges(5,3) * pkin(3);
t75 = -Ifges(4,5) + t135;
t125 = t75 * t74;
t144 = t105 + t125;
t142 = -Ifges(4,5) / 0.2e1 + t135 / 0.2e1;
t141 = -pkin(3) / 0.2e1;
t140 = pkin(3) / 0.2e1;
t139 = mrSges(6,1) / 0.4e1;
t137 = pkin(1) + pkin(4);
t67 = sin(t79);
t134 = pkin(2) * t67;
t68 = sin(t80);
t133 = pkin(2) * t68;
t56 = mrSges(6,1) * t60;
t86 = cos(qJ(4));
t128 = mrSges(6,1) * t86;
t127 = mrSges(6,2) * t63;
t83 = sin(qJ(4));
t126 = mrSges(6,2) * t83;
t88 = cos(qJ(2));
t124 = t84 * t88;
t122 = t40 + t39;
t41 = t62 * t141;
t121 = -t129 / 0.2e1 + t41;
t120 = t55 + t56;
t108 = m(6) * pkin(4) + mrSges(5,1);
t111 = -t133 / 0.2e1;
t91 = 0.2e1 * qJ(2);
t81 = t91 + qJ(3);
t69 = sin(t81);
t73 = cos(t81);
t78 = m(4) + t90;
t85 = sin(qJ(2));
t66 = 0.2e1 * t82;
t93 = mrSges(6,2) * t41 + t56 * t141 + Ifges(4,4) * cos(t66) - mrSges(4,2) * pkin(1) * t74 - (pkin(1) * t65 + pkin(3) * t108) * t70 - (t90 * t92 - Ifges(4,1) + Ifges(4,2)) * sin(t66) / 0.2e1 + t123;
t1 = t93 - ((pkin(1) * t78 + t108) * pkin(2) + mrSges(3,1) * pkin(1)) * t85 + Ifges(3,4) * cos(t91) - (t78 * pkin(2) ^ 2 - Ifges(3,1) + Ifges(3,2)) * sin(t91) / 0.2e1 - t65 * pkin(2) * t69 - mrSges(4,2) * pkin(2) * t73 - mrSges(3,2) * t88 * pkin(1) + t145 * mrSges(6,2) + (-t134 / 0.2e1 + t111) * mrSges(6,1);
t119 = t1 * qJD(1);
t2 = ((-t84 / 0.2e1 - t69 / 0.2e1) * t65 + (-t73 / 0.2e1 - t87 / 0.2e1) * mrSges(4,2)) * pkin(2) + t93;
t118 = t2 * qJD(1);
t64 = t87 * pkin(3) + pkin(2);
t23 = -pkin(3) * t124 - t64 * t85;
t3 = (-t23 * t83 / 0.2e1 + (t63 / 0.4e1 - t62 / 0.4e1) * pkin(3) + (t72 / 0.4e1 - t71 / 0.4e1) * pkin(2)) * mrSges(6,2) + (t23 * t86 / 0.2e1 + (-t61 / 0.4e1 - t60 / 0.4e1) * pkin(3) + (-t68 / 0.4e1 - t67 / 0.4e1) * pkin(2)) * mrSges(6,1);
t117 = t3 * qJD(1);
t103 = mrSges(6,1) * t83 + mrSges(6,2) * t86;
t94 = (t72 / 0.2e1 + t71 / 0.2e1) * mrSges(6,2) + (-t68 / 0.2e1 + t67 / 0.2e1) * mrSges(6,1);
t5 = t94 * pkin(2) + t103 * t137 + t106;
t116 = t5 * qJD(1);
t101 = -t85 * t87 - t124;
t100 = t101 * t128;
t107 = t114 / 0.4e1 + (-t127 / 0.4e1 + (t61 + t60) * t139) * pkin(3);
t99 = t101 * t126;
t7 = (-t100 / 0.2e1 + t99 / 0.2e1) * pkin(3) + t107;
t115 = t7 * qJD(1);
t104 = -mrSges(6,1) * t61 + t127;
t98 = pkin(3) * (t104 + t120);
t14 = (t75 / 0.2e1 - t142) * t74;
t95 = t122 * mrSges(6,1) / 0.2e1 - t121 * mrSges(6,2) / 0.2e1;
t9 = -t98 / 0.4e1 + t95;
t97 = t14 * qJD(1) - t25 * qJD(2) + t9 * qJD(4);
t53 = t134 / 0.2e1;
t11 = t12 * qJD(3);
t10 = t125 / 0.2e1 + t105 + t142 * t74;
t8 = t98 / 0.4e1 + t95;
t6 = t100 * t140 + t99 * t141 + t107;
t4 = (t128 / 0.2e1 - t126 / 0.2e1) * t23 + t107 + (t68 + t67) * pkin(2) * t139 + (-t131 / 0.4e1 + t132 / 0.4e1) * mrSges(6,2);
t13 = [t1 * qJD(2) + t2 * qJD(3) - t5 * qJD(4), t10 * qJD(3) + t4 * qJD(4) + t119 + (-Ifges(3,5) * t88 + Ifges(3,6) * t85 + ((mrSges(4,3) + mrSges(5,3)) * t88 + t94) * pkin(2) + t144) * qJD(2), t10 * qJD(2) + t144 * qJD(3) + t6 * qJD(4) + t118, -t116 + t4 * qJD(2) + t6 * qJD(3) - t103 * qJD(4) * (-t85 * t84 * pkin(3) + t64 * t88 + t137); t14 * qJD(3) - t3 * qJD(4) - t119, -t25 * qJD(3), t97 + (-t12 - t25) * qJD(3), -t117 + t9 * qJD(3) + (-(t109 + t41 + t145) * mrSges(6,2) + (t130 / 0.2e1 + t133 / 0.2e1 + t39 + t53) * mrSges(6,1)) * qJD(4); -t14 * qJD(2) + t7 * qJD(4) - t118, t11 - t97, t12 * qJD(2) + t11, t115 - t9 * qJD(2) + (-t104 + t120) * qJD(4) * t140; t3 * qJD(2) - t7 * qJD(3) + t116, t117 + (-(-t131 / 0.2e1 + t54 + t121) * mrSges(6,2) + (t111 + t53 + t122) * mrSges(6,1)) * qJD(2) + t8 * qJD(3), -t115 + t8 * qJD(2) + qJD(3) * t98 / 0.2e1, 0;];
Cq = t13;
