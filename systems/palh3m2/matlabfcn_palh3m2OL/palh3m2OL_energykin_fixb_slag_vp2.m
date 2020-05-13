% Calculate kinetic energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energykin_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_energykin_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energykin_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_energykin_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_energykin_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_energykin_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:27
% EndTime: 2020-05-07 04:33:28
% DurationCPUTime: 1.05s
% Computational Cost: add. (510->127), mult. (1146->220), div. (0->0), fcn. (872->16), ass. (0->54)
t150 = sin(qJ(8));
t126 = qJD(2) + qJD(7);
t149 = pkin(3) * t126;
t148 = pkin(1) * qJD(2);
t147 = qJD(2) ^ 2 * pkin(1) ^ 2;
t127 = qJD(2) + qJD(3);
t134 = sin(qJ(3));
t146 = t134 * t148;
t142 = cos(qJ(2));
t122 = (-pkin(1) * t142 - pkin(12)) * qJD(1);
t135 = sin(qJ(2));
t141 = cos(qJ(3));
t112 = (t134 * t135 - t141 * t142) * qJD(1);
t114 = (-t134 * t142 - t135 * t141) * qJD(1);
t133 = sin(qJ(4));
t140 = cos(qJ(4));
t104 = t112 * t140 - t114 * t133;
t120 = pkin(4) * t127 - t141 * t148;
t110 = t133 * t120 - t140 * t146;
t109 = t120 * t140 + t133 * t146;
t108 = -pkin(4) * t112 + t122;
t139 = cos(qJ(5));
t138 = cos(qJ(6));
t137 = cos(qJ(7));
t136 = cos(qJ(8));
t132 = sin(qJ(5));
t131 = sin(qJ(6));
t130 = sin(qJ(7));
t129 = cos(pkin(15));
t128 = sin(pkin(15));
t125 = qJD(4) + t127;
t124 = qJD(8) + t126;
t121 = t122 ^ 2;
t119 = -t128 * t150 - t129 * t136;
t118 = t128 * t136 - t129 * t150;
t117 = t129 * t149 + t137 * t148;
t116 = t128 * t149 + t130 * t148;
t113 = (t130 * t142 + t135 * t137) * qJD(1);
t111 = (-t130 * t135 + t137 * t142) * qJD(1);
t107 = pkin(10) * t125 + t110;
t106 = -pkin(8) * t125 - t109;
t105 = t112 * t133 + t114 * t140;
t103 = qJD(5) - t104;
t102 = t122 + (-t111 * t129 - t113 * t128) * pkin(3);
t101 = t105 * t139 + t125 * t132;
t100 = -t105 * t132 + t125 * t139;
t99 = t116 * t119 + t117 * t118;
t98 = -t116 * t118 + t117 * t119;
t97 = t111 * t118 + t113 * t119;
t96 = t111 * t119 - t113 * t118;
t95 = -pkin(8) * t104 - pkin(10) * t105 + t108;
t94 = t107 * t139 + t132 * t95;
t93 = -t107 * t132 + t139 * t95;
t1 = Ifges(7,3) * qJD(6) ^ 2 / 0.2e1 + m(6) * (t106 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(9) * (t102 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + Ifges(4,3) * t127 ^ 2 / 0.2e1 + m(8) * (t121 + (t130 ^ 2 + t137 ^ 2) * t147) / 0.2e1 + m(4) * (t121 + (t134 ^ 2 + t141 ^ 2) * t147) / 0.2e1 + Ifges(8,3) * t126 ^ 2 / 0.2e1 + (t102 * mrSges(9,2) - t98 * mrSges(9,3) + Ifges(9,1) * t97 / 0.2e1) * t97 + (t109 * mrSges(5,1) - t110 * mrSges(5,2) + Ifges(5,3) * t125 / 0.2e1) * t125 + (t93 * mrSges(6,1) - t94 * mrSges(6,2) + Ifges(6,3) * t103 / 0.2e1) * t103 + (-t102 * mrSges(9,1) + t99 * mrSges(9,3) + Ifges(9,4) * t97 + Ifges(9,2) * t96 / 0.2e1) * t96 + (t122 * mrSges(4,2) + Ifges(4,5) * t127 + Ifges(4,1) * t114 / 0.2e1) * t114 + (t122 * mrSges(8,2) + Ifges(8,5) * t126 + Ifges(8,1) * t113 / 0.2e1) * t113 + (t108 * mrSges(5,2) - t109 * mrSges(5,3) + Ifges(5,5) * t125 + Ifges(5,1) * t105 / 0.2e1) * t105 + (t106 * mrSges(6,2) - t93 * mrSges(6,3) + Ifges(6,5) * t103 + Ifges(6,1) * t101 / 0.2e1) * t101 + (t98 * mrSges(9,1) - t99 * mrSges(9,2) + Ifges(9,5) * t97 + Ifges(9,6) * t96 + Ifges(9,3) * t124 / 0.2e1) * t124 + (-t122 * mrSges(4,1) + Ifges(4,4) * t114 + Ifges(4,6) * t127 + Ifges(4,2) * t112 / 0.2e1) * t112 + (-t122 * mrSges(8,1) + Ifges(8,4) * t113 + Ifges(8,6) * t126 + Ifges(8,2) * t111 / 0.2e1) * t111 + (-t108 * mrSges(5,1) + t110 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,6) * t125 + Ifges(5,2) * t104 / 0.2e1) * t104 + (-t106 * mrSges(6,1) + t94 * mrSges(6,3) + Ifges(6,4) * t101 + Ifges(6,6) * t103 + Ifges(6,2) * t100 / 0.2e1) * t100 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t130 * (-mrSges(8,2) * t126 + mrSges(8,3) * t111) - t134 * (-mrSges(4,2) * t127 + mrSges(4,3) * t112) + t137 * (mrSges(8,1) * t126 - mrSges(8,3) * t113) - t141 * (mrSges(4,1) * t127 - t114 * mrSges(4,3))) * pkin(1)) * qJD(2) + (qJD(6) * (Ifges(7,5) * t131 + Ifges(7,6) * t138) + qJD(2) * (Ifges(3,5) * t135 + Ifges(3,6) * t142) + (m(7) * pkin(6) ^ 2 / 0.2e1 + m(3) * pkin(12) ^ 2 / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(12) * mrSges(3,1) + Ifges(3,2) * t142 / 0.2e1) * t142 + (-pkin(6) * mrSges(7,1) + Ifges(7,2) * t138 / 0.2e1) * t138 + (-pkin(12) * mrSges(3,2) + Ifges(3,4) * t142 + Ifges(3,1) * t135 / 0.2e1) * t135 + (pkin(6) * mrSges(7,2) + Ifges(7,4) * t138 + Ifges(7,1) * t131 / 0.2e1) * t131) * qJD(1)) * qJD(1);
T = t1;
