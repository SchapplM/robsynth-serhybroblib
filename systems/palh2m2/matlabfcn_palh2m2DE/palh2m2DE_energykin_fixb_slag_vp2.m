% Calculate kinetic energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2DE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:23
% EndTime: 2020-05-03 01:06:23
% DurationCPUTime: 0.32s
% Computational Cost: add. (95->53), mult. (217->94), div. (0->0), fcn. (81->6), ass. (0->24)
t124 = m(4) / 0.2e1;
t123 = m(5) / 0.2e1;
t122 = m(6) / 0.2e1;
t121 = pkin(4) * qJD(2);
t120 = pkin(5) * qJD(3);
t116 = cos(qJ(2));
t109 = -t116 * pkin(4) - pkin(1);
t108 = -pkin(2) + t109;
t115 = cos(qJ(3));
t114 = cos(qJ(4));
t113 = sin(qJ(2));
t112 = sin(qJ(3));
t111 = sin(qJ(4));
t110 = qJD(1) + qJD(4);
t107 = t115 * t120 + t116 * t121;
t106 = -t112 * t120 - t113 * t121;
t105 = -t115 * pkin(5) + t108;
t104 = t107 ^ 2;
t103 = t112 * t113 + t115 * t116;
t102 = t112 * t116 - t115 * t113;
t101 = (-pkin(3) + t105) * qJD(1);
t100 = -t111 * t101 + t114 * t106;
t99 = -t114 * t101 - t111 * t106;
t1 = Ifges(5,3) * qJD(3) ^ 2 / 0.2e1 + (t106 ^ 2 + t104) * t122 + m(7) * (t100 ^ 2 + t99 ^ 2 + t104) / 0.2e1 + (t99 * mrSges(7,1) - t100 * mrSges(7,2) + Ifges(7,3) * t110 / 0.2e1) * t110 + ((mrSges(5,1) * t103 - mrSges(5,2) * t102) * qJD(3) * pkin(4) + (((t113 ^ 2 + t116 ^ 2) * t124 + (t102 ^ 2 + t103 ^ 2) * t123) * pkin(4) ^ 2 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + (t106 * mrSges(6,3) + qJD(3) * (Ifges(5,5) * t112 + Ifges(5,6) * t115) + (Ifges(3,5) * t113 + Ifges(3,6) * t116 + (-mrSges(4,3) * t113 + (t102 * t115 - t103 * t112) * mrSges(5,3)) * pkin(4)) * qJD(2) + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + t108 ^ 2 * t123 + m(3) * pkin(1) ^ 2 / 0.2e1 + (t109 * t124 - mrSges(4,1)) * t109 + (t105 * t122 - mrSges(6,1)) * t105 + (pkin(1) * mrSges(3,1) + Ifges(3,2) * t116 / 0.2e1) * t116 + (-t108 * mrSges(5,1) + Ifges(5,2) * t115 / 0.2e1) * t115 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t116 + Ifges(3,1) * t113 / 0.2e1) * t113 + (t108 * mrSges(5,2) + Ifges(5,4) * t115 + Ifges(5,1) * t112 / 0.2e1) * t112) * qJD(1)) * qJD(1);
T = t1;
