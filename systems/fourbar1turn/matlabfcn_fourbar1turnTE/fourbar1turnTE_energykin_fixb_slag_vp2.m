% Calculate kinetic energy for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnTE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:29
% EndTime: 2020-04-12 19:18:30
% DurationCPUTime: 0.68s
% Computational Cost: add. (3675->80), mult. (5197->186), div. (216->12), fcn. (1398->4), ass. (0->55)
t160 = pkin(3) ^ 2;
t159 = -2 * pkin(2);
t122 = qJD(1) ^ 2;
t158 = (-pkin(3) - pkin(4));
t157 = (pkin(4) - pkin(3));
t128 = pkin(2) ^ 2;
t129 = pkin(1) ^ 2;
t120 = cos(qJ(2));
t150 = pkin(2) * t120;
t143 = -0.2e1 * pkin(1) * t150 + t129;
t114 = t128 + t143;
t142 = -pkin(4) ^ 2 + t160;
t111 = t114 - t142;
t115 = pkin(1) - t150;
t119 = sin(qJ(2));
t108 = ((pkin(2) - t158) * (pkin(2) + t158)) + t143;
t109 = ((pkin(2) - t157) * (pkin(2) + t157)) + t143;
t130 = sqrt(-t108 * t109);
t145 = t119 * t130;
t102 = -pkin(2) * t145 + t111 * t115;
t156 = -t102 / 0.2e1;
t146 = t119 * t111;
t103 = pkin(2) * t146 + t115 * t130;
t155 = t103 / 0.2e1;
t112 = 0.1e1 / t114;
t154 = -t112 / 0.2e1;
t153 = -t120 / 0.2e1;
t151 = pkin(1) * t119;
t149 = pkin(2) * qJD(2);
t110 = t114 + t142;
t116 = pkin(1) * t120 - pkin(2);
t104 = t110 * t151 - t116 * t130;
t148 = t104 * t119;
t147 = t112 * t119 ^ 2;
t144 = t120 * t130;
t141 = qJD(1) * t112;
t140 = t120 * qJD(1);
t113 = 0.1e1 / t114 ^ 2;
t139 = t113 * t159;
t138 = (-t108 - t109) * t149 * t151 / t130 * t112;
t126 = 0.1e1 / pkin(3);
t137 = t126 * t141;
t124 = 0.1e1 / pkin(4);
t136 = t124 * t141;
t133 = t138 / 0.2e1;
t101 = -pkin(1) * t145 - t110 * t116;
t131 = t101 ^ 2;
t100 = t104 ^ 2;
t99 = 0.1e1 / t102 ^ 2;
t98 = 0.1e1 / t131;
t95 = (t148 / 0.2e1 + t101 * t153) * t137;
t94 = (-t101 * t119 / 0.2e1 + t104 * t153) * t137;
t93 = 0.2e1 * (-(t115 * t133 + (t128 * pkin(1) * t147 + ((t111 * t120 + t145) * t112 / 0.2e1 - t103 * t113 * t151) * pkin(2)) * qJD(2)) / t102 - (t119 * t133 + ((-t144 + t146) * t154 + (t102 * t113 - t112 * t115) * t151) * qJD(2)) * pkin(2) * t103 * t99) * pkin(4) * t114 * t124 / (t103 ^ 2 * t99 + 0.1e1);
t92 = qJD(2) + ((-t116 * t138 + (0.2e1 * t129 * pkin(2) * t147 + ((t110 * t120 + t145) * t112 + t139 * t148) * pkin(1)) * qJD(2)) / t101 - (-t119 * t138 + (-t112 * t144 + ((t116 * t159 + t110) * t112 + t101 * t139) * t119) * qJD(2)) * pkin(1) * t104 * t98) * pkin(3) * t114 * t126 / (t100 * t98 + 0.1e1);
t1 = t119 * qJD(1) * ((Ifges(3,5) * qJD(2)) + (Ifges(3,1) * t119 + Ifges(3,4) * t120) * qJD(1)) / 0.2e1 + ((Ifges(3,6) * qJD(2)) + (Ifges(3,4) * t119 + Ifges(3,2) * t120) * qJD(1)) * t140 / 0.2e1 + qJD(2) * ((Ifges(3,3) * qJD(2)) + (Ifges(3,5) * t119 + Ifges(3,6) * t120) * qJD(1)) / 0.2e1 + Ifges(4,1) * t94 ^ 2 / 0.2e1 + Ifges(4,2) * t95 ^ 2 / 0.2e1 + t103 * (Ifges(5,5) * t93 + (Ifges(5,1) * t155 + Ifges(5,4) * t156) * t136) * t136 / 0.4e1 - t102 * (Ifges(5,6) * t93 + (Ifges(5,4) * t155 + Ifges(5,2) * t156) * t136) * t136 / 0.4e1 + t93 * (Ifges(5,3) * t93 + (Ifges(5,5) * t155 + Ifges(5,6) * t156) * t136) / 0.2e1 - pkin(2) * (-mrSges(4,1) * t95 + mrSges(4,2) * t94) * t140 - t122 * pkin(1) * (mrSges(5,2) * t155 + t102 * mrSges(5,1) / 0.2e1) * t124 * t112 + m(4) * (t120 ^ 2 * t122 + (t100 / 0.4e1 + t131 / 0.4e1) / t160 * (qJD(2) ^ 2) * t113) * t128 / 0.2e1 + t94 * t95 * Ifges(4,4) + (m(5) * t129 + Ifges(2,3)) * t122 / 0.2e1 + (Ifges(4,5) * t94 + Ifges(4,6) * t95 + Ifges(4,3) * t92 / 0.2e1) * t92 + (t104 * (-mrSges(4,2) * t92 + mrSges(4,3) * t95) + t101 * (mrSges(4,1) * t92 - t94 * mrSges(4,3))) * t126 * t149 * t154;
T = t1;
