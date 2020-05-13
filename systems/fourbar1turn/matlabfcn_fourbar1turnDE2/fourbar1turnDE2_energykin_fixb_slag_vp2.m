% Calculate kinetic energy for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnDE2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:24
% EndTime: 2020-04-12 19:33:25
% DurationCPUTime: 0.74s
% Computational Cost: add. (5309->75), mult. (7589->167), div. (370->13), fcn. (2078->8), ass. (0->53)
t128 = cos(qJ(2));
t156 = pkin(2) * t128;
t137 = pkin(2) ^ 2;
t138 = pkin(1) ^ 2;
t150 = -0.2e1 * pkin(1) * t156 + t138;
t122 = t137 + t150;
t161 = pkin(3) ^ 2;
t162 = pkin(4) ^ 2;
t149 = t161 - t162;
t118 = t122 + t149;
t124 = pkin(1) * t128 - pkin(2);
t127 = sin(qJ(2));
t159 = -pkin(3) - pkin(4);
t116 = (pkin(2) - t159) * (pkin(2) + t159) + t150;
t158 = -pkin(3) + pkin(4);
t117 = (pkin(2) - t158) * (pkin(2) + t158) + t150;
t139 = sqrt(-t116 * t117);
t152 = t127 * t139;
t109 = -pkin(1) * t152 - t118 * t124;
t163 = t109 ^ 2;
t123 = pkin(1) - t156;
t119 = t122 - t149;
t153 = t127 * t119;
t111 = pkin(2) * t153 + t123 * t139;
t107 = t111 ^ 2;
t160 = -0.2e1 * pkin(2);
t126 = t127 ^ 2;
t157 = pkin(1) * t127;
t112 = t118 * t157 - t124 * t139;
t155 = t112 * t127;
t120 = 0.1e1 / t122;
t154 = t120 * t126;
t151 = t128 * t139;
t121 = 0.1e1 / t122 ^ 2;
t148 = t121 * t160;
t108 = t112 ^ 2;
t135 = 0.1e1 / pkin(3);
t147 = t120 * t135 * ((t108 + t163) / t161 * t121) ^ (-0.1e1 / 0.2e1);
t132 = 0.1e1 / pkin(4);
t110 = -pkin(2) * t152 + t119 * t123;
t140 = t110 ^ 2;
t146 = t132 * t120 * ((t107 + t140) / t162 * t121) ^ (-0.1e1 / 0.2e1);
t145 = (-t116 - t117) * qJD(2) * pkin(2) * t157 / t139 * t120;
t142 = qJD(1) * t147;
t141 = t145 / 0.2e1;
t130 = qJD(1) ^ 2;
t106 = 0.1e1 / t140;
t105 = 0.1e1 / t163;
t97 = (-t109 * t128 + t155) * t142;
t96 = (-t109 * t127 - t112 * t128) * t142;
t95 = 0.2e1 * (-(t123 * t141 + (t137 * pkin(1) * t154 + ((t119 * t128 + t152) * t120 / 0.2e1 - t111 * t121 * t157) * pkin(2)) * qJD(2)) / t110 - (t127 * t141 + (-(-t151 + t153) * t120 / 0.2e1 + (t110 * t121 - t120 * t123) * t157) * qJD(2)) * pkin(2) * t111 * t106) * pkin(4) / (t106 * t107 + 0.1e1) * t122 * t132;
t94 = qJD(2) + ((-t124 * t145 + (0.2e1 * t138 * pkin(2) * t154 + ((t118 * t128 + t152) * t120 + t148 * t155) * pkin(1)) * qJD(2)) / t109 - (-t127 * t145 + (-t120 * t151 + ((t124 * t160 + t118) * t120 + t109 * t148) * t127) * qJD(2)) * pkin(1) * t112 * t105) * pkin(3) / (t105 * t108 + 0.1e1) * t122 * t135;
t1 = Ifges(4,2) * t97 ^ 2 / 0.2e1 + Ifges(5,3) * t95 ^ 2 / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(5) * t138 / 0.2e1) * t130 + (Ifges(4,4) * t97 + Ifges(4,1) * t96 / 0.2e1) * t96 + m(4) * (t128 ^ 2 * t130 + qJD(2) ^ 2) * t137 / 0.2e1 + (Ifges(4,5) * t96 + Ifges(4,6) * t97 + Ifges(4,3) * t94 / 0.2e1) * t94 + (-(-mrSges(4,1) * t97 + t96 * mrSges(4,2)) * t156 + t95 * (Ifges(5,5) * t111 - Ifges(5,6) * t110) * t146 + (Ifges(3,1) * t126 / 0.2e1 + (Ifges(3,4) * t127 + Ifges(3,2) * t128 / 0.2e1) * t128 + (-pkin(1) * (mrSges(5,1) * t110 + mrSges(5,2) * t111) + (Ifges(5,1) * t107 / 0.2e1 + (-Ifges(5,4) * t111 + Ifges(5,2) * t110 / 0.2e1) * t110) * t146) * t146) * qJD(1)) * qJD(1) + (Ifges(3,3) * qJD(2) / 0.2e1 + (Ifges(3,5) * t127 + Ifges(3,6) * t128) * qJD(1) + (-t112 * (-mrSges(4,2) * t94 + mrSges(4,3) * t97) - t109 * (mrSges(4,1) * t94 - mrSges(4,3) * t96)) * pkin(2) * t147) * qJD(2);
T = t1;
