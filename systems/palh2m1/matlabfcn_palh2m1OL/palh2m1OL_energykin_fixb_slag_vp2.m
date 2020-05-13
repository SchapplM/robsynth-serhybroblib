% Calculate kinetic energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:59:02
% EndTime: 2020-05-02 23:59:03
% DurationCPUTime: 0.56s
% Computational Cost: add. (332->78), mult. (716->140), div. (0->0), fcn. (550->8), ass. (0->32)
t127 = sin(qJ(3));
t138 = t127 * pkin(2);
t137 = pkin(3) * qJD(3);
t128 = sin(qJ(2));
t136 = t128 * t127;
t131 = cos(qJ(3));
t120 = t131 * pkin(2) + pkin(3);
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t113 = (t120 * t126 + t130 * t138) * qJD(2) + t126 * t137;
t124 = qJD(2) + qJD(3);
t132 = cos(qJ(2));
t118 = (-t131 * t132 + t136) * qJD(1);
t119 = (-t127 * t132 - t128 * t131) * qJD(1);
t109 = t130 * t118 - t126 * t119;
t114 = (t120 * t130 - t126 * t138) * qJD(2) + t130 * t137;
t134 = qJD(1) ^ 2;
t129 = cos(qJ(5));
t125 = sin(qJ(5));
t123 = qJD(4) + t124;
t121 = t132 * pkin(2) + pkin(1);
t116 = (t131 * pkin(3) + pkin(2)) * t132 - pkin(3) * t136 + pkin(1);
t112 = -t123 * pkin(4) - t114;
t111 = t123 * pkin(6) + t113;
t110 = t126 * t118 + t130 * t119;
t108 = qJD(5) - t109;
t107 = t129 * t110 + t125 * t123;
t106 = -t125 * t110 + t129 * t123;
t105 = -t109 * pkin(4) - t110 * pkin(6) + t116 * qJD(1);
t104 = t125 * t105 + t129 * t111;
t103 = t129 * t105 - t125 * t111;
t1 = m(5) * (t116 ^ 2 * t134 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t121 ^ 2 * t134 + (t127 ^ 2 + t131 ^ 2) * pkin(2) ^ 2 * qJD(2) ^ 2) / 0.2e1 + Ifges(4,3) * t124 ^ 2 / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t112 ^ 2) / 0.2e1 + (-t114 * mrSges(5,3) + Ifges(5,1) * t110 / 0.2e1) * t110 + (Ifges(4,5) * t124 + Ifges(4,1) * t119 / 0.2e1) * t119 + (t103 * mrSges(6,1) - t104 * mrSges(6,2) + Ifges(6,3) * t108 / 0.2e1) * t108 + (Ifges(4,4) * t119 + Ifges(4,6) * t124 + Ifges(4,2) * t118 / 0.2e1) * t118 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t127 * (-t124 * mrSges(4,2) + t118 * mrSges(4,3)) + t131 * (t124 * mrSges(4,1) - t119 * mrSges(4,3))) * pkin(2)) * qJD(2) + (t114 * mrSges(5,1) - t113 * mrSges(5,2) + Ifges(5,5) * t110 + Ifges(5,3) * t123 / 0.2e1) * t123 + (t112 * mrSges(6,2) - t103 * mrSges(6,3) + Ifges(6,5) * t108 + Ifges(6,1) * t107 / 0.2e1) * t107 + (t113 * mrSges(5,3) + Ifges(5,4) * t110 + Ifges(5,6) * t123 + Ifges(5,2) * t109 / 0.2e1) * t109 + (-t112 * mrSges(6,1) + t104 * mrSges(6,3) + Ifges(6,4) * t107 + Ifges(6,6) * t108 + Ifges(6,2) * t106 / 0.2e1) * t106 + (t116 * (-t109 * mrSges(5,1) + t110 * mrSges(5,2)) + t121 * (-t118 * mrSges(4,1) + t119 * mrSges(4,2)) + ((pkin(1) * mrSges(3,1) + Ifges(3,2) * t132 / 0.2e1) * t132 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t132 + Ifges(3,1) * t128 / 0.2e1) * t128) * qJD(1) + qJD(2) * (-Ifges(3,5) * t128 - Ifges(3,6) * t132)) * qJD(1) + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * t134 / 0.2e1;
T = t1;
