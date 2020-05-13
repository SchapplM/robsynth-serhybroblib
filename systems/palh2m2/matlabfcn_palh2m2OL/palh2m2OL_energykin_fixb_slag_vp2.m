% Calculate kinetic energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
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
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2OL_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:13:44
% EndTime: 2020-05-03 01:13:45
% DurationCPUTime: 0.90s
% Computational Cost: add. (638->110), mult. (1411->192), div. (0->0), fcn. (1202->10), ass. (0->46)
t169 = pkin(5) * qJD(4);
t153 = sin(qJ(4));
t154 = sin(qJ(3));
t168 = t153 * t154;
t155 = sin(qJ(2));
t167 = t155 * t154;
t166 = pkin(2) * t153;
t165 = pkin(2) * qJD(3);
t150 = qJD(2) + qJD(3);
t159 = cos(qJ(3));
t147 = t159 * pkin(4) + pkin(2);
t158 = cos(qJ(4));
t164 = -pkin(4) * t168 + t147 * t158;
t138 = pkin(5) + t164;
t139 = t158 * t154 * pkin(4) + t147 * t153;
t146 = t158 * pkin(2) + pkin(5);
t152 = sin(qJ(5));
t157 = cos(qJ(5));
t128 = (t138 * t157 - t139 * t152) * qJD(2) + (t146 * t157 - t152 * t166) * qJD(3) + t157 * t169;
t149 = qJD(4) + t150;
t160 = cos(qJ(2));
t140 = (t159 * t160 - t167) * qJD(1);
t141 = (t154 * t160 + t155 * t159) * qJD(1);
t132 = t158 * t140 - t153 * t141;
t133 = t153 * t140 + t158 * t141;
t124 = t157 * t132 - t152 * t133;
t162 = qJD(1) ^ 2;
t156 = cos(qJ(6));
t151 = sin(qJ(6));
t145 = t158 * pkin(5) + pkin(2);
t144 = qJD(5) + t149;
t142 = (-pkin(4) * t160 - pkin(1)) * qJD(1);
t137 = (-t159 * pkin(2) - pkin(4)) * t160 + pkin(2) * t167 - pkin(1);
t135 = t164 * qJD(2) + t158 * t165;
t134 = t139 * qJD(2) + t153 * t165;
t130 = (pkin(5) * t168 - t145 * t159 - pkin(4)) * t160 + (t159 * t153 * pkin(5) + t154 * t145) * t155 - pkin(1);
t127 = (t138 * t152 + t157 * t139) * qJD(2) + (t146 * t152 + t157 * t166) * qJD(3) + t152 * t169;
t126 = t144 * pkin(3) + t128;
t125 = t152 * t132 + t157 * t133;
t123 = qJD(6) + t124;
t122 = t156 * t125 - t151 * t144;
t121 = -t151 * t125 - t156 * t144;
t120 = -t124 * pkin(3) + t130 * qJD(1);
t119 = -t151 * t120 + t156 * t127;
t118 = -t156 * t120 - t151 * t127;
t1 = m(5) * (t137 ^ 2 * t162 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(6) * (t130 ^ 2 * t162 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + (t154 ^ 2 + t159 ^ 2) * pkin(4) ^ 2 * qJD(2) ^ 2) / 0.2e1 + Ifges(4,3) * t150 ^ 2 / 0.2e1 + m(7) * (t118 ^ 2 + t119 ^ 2 + t126 ^ 2) / 0.2e1 + (-t135 * mrSges(5,3) + Ifges(5,1) * t133 / 0.2e1) * t133 + (-t128 * mrSges(6,3) + Ifges(6,1) * t125 / 0.2e1) * t125 + (t118 * mrSges(7,1) - t119 * mrSges(7,2) + Ifges(7,3) * t123 / 0.2e1) * t123 + (t142 * mrSges(4,2) + Ifges(4,5) * t150 + Ifges(4,1) * t141 / 0.2e1) * t141 + (-t142 * mrSges(4,1) + Ifges(4,4) * t141 + Ifges(4,6) * t150 + Ifges(4,2) * t140 / 0.2e1) * t140 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t154 * (-t150 * mrSges(4,2) + t140 * mrSges(4,3)) + t159 * (t150 * mrSges(4,1) - t141 * mrSges(4,3))) * pkin(4)) * qJD(2) + (t135 * mrSges(5,1) - t134 * mrSges(5,2) + Ifges(5,5) * t133 + Ifges(5,3) * t149 / 0.2e1) * t149 + (t128 * mrSges(6,1) - t127 * mrSges(6,2) + Ifges(6,5) * t125 + Ifges(6,3) * t144 / 0.2e1) * t144 + (t126 * mrSges(7,2) - t118 * mrSges(7,3) + Ifges(7,5) * t123 + Ifges(7,1) * t122 / 0.2e1) * t122 + (t134 * mrSges(5,3) + Ifges(5,4) * t133 + Ifges(5,6) * t149 + Ifges(5,2) * t132 / 0.2e1) * t132 + (t127 * mrSges(6,3) + Ifges(6,4) * t125 + Ifges(6,6) * t144 + Ifges(6,2) * t124 / 0.2e1) * t124 + (-t126 * mrSges(7,1) + t119 * mrSges(7,3) + Ifges(7,4) * t122 + Ifges(7,6) * t123 + Ifges(7,2) * t121 / 0.2e1) * t121 + (t130 * (-t124 * mrSges(6,1) + t125 * mrSges(6,2)) + t137 * (-t132 * mrSges(5,1) + t133 * mrSges(5,2)) + ((pkin(1) * mrSges(3,1) + Ifges(3,2) * t160 / 0.2e1) * t160 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t160 + Ifges(3,1) * t155 / 0.2e1) * t155) * qJD(1) + qJD(2) * (Ifges(3,5) * t155 + Ifges(3,6) * t160)) * qJD(1) + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * t162 / 0.2e1;
T = t1;
