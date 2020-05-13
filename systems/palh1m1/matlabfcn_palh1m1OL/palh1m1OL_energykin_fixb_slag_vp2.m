% Calculate kinetic energy for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:23
% EndTime: 2020-04-15 19:28:26
% DurationCPUTime: 2.00s
% Computational Cost: add. (610->158), mult. (1383->278), div. (0->0), fcn. (1044->20), ass. (0->66)
t153 = qJD(2) + qJD(7);
t184 = pkin(4) * t153;
t183 = sin(qJ(10));
t182 = pkin(1) * qJD(2);
t181 = qJD(1) * pkin(15);
t180 = qJD(2) ^ 2 * pkin(1) ^ 2;
t154 = qJD(2) + qJD(3);
t152 = qJD(2) + qJD(8);
t172 = cos(qJ(3));
t179 = t172 * t182;
t165 = sin(qJ(2));
t145 = t165 * qJD(1) * pkin(1) - t181;
t164 = sin(qJ(3));
t173 = cos(qJ(2));
t134 = (t164 * t173 + t165 * t172) * qJD(1);
t137 = (-t164 * t165 + t172 * t173) * qJD(1);
t163 = sin(qJ(4));
t171 = cos(qJ(4));
t123 = -t134 * t163 + t137 * t171;
t143 = pkin(5) * t154 + t164 * t182;
t130 = t163 * t143 - t171 * t179;
t128 = -pkin(5) * t137 + t145;
t129 = t143 * t171 + t163 * t179;
t170 = cos(qJ(5));
t169 = cos(qJ(6));
t168 = cos(qJ(7));
t167 = cos(qJ(8));
t166 = cos(qJ(9));
t162 = sin(qJ(5));
t161 = sin(qJ(6));
t160 = sin(qJ(7));
t159 = sin(qJ(8));
t158 = sin(qJ(9));
t157 = cos(qJ(10));
t156 = cos(pkin(19));
t155 = sin(pkin(19));
t150 = qJD(4) + t154;
t149 = qJD(9) + t152;
t148 = qJD(10) + t153;
t144 = t145 ^ 2;
t142 = -t155 * t183 - t156 * t157;
t141 = t155 * t157 - t156 * t183;
t140 = t156 * t184 + t168 * t182;
t139 = t155 * t184 + t160 * t182;
t136 = (-t160 * t165 + t168 * t173) * qJD(1);
t135 = (-t159 * t165 + t167 * t173) * qJD(1);
t133 = (-t160 * t173 - t165 * t168) * qJD(1);
t132 = (-t159 * t173 - t165 * t167) * qJD(1);
t131 = -pkin(2) * t132 - t181;
t127 = pkin(11) * t150 + t130;
t126 = -pkin(9) * t150 - t129;
t125 = t134 * t171 + t137 * t163;
t124 = -t132 * t158 - t135 * t166;
t122 = -t132 * t166 + t135 * t158;
t121 = qJD(5) - t123;
t120 = (-t133 * t156 - t136 * t155) * pkin(4) + t145;
t119 = t125 * t170 + t150 * t162;
t118 = -t125 * t162 + t150 * t170;
t117 = t139 * t142 + t140 * t141;
t116 = -t139 * t141 + t140 * t142;
t115 = t133 * t141 + t136 * t142;
t114 = t133 * t142 - t136 * t141;
t113 = -pkin(9) * t123 - pkin(11) * t125 + t128;
t112 = t113 * t162 + t127 * t170;
t111 = t113 * t170 - t127 * t162;
t1 = Ifges(4,3) * t154 ^ 2 / 0.2e1 + (t111 * mrSges(6,1) - t112 * mrSges(6,2) + Ifges(6,3) * t121 / 0.2e1) * t121 + Ifges(10,3) * t149 ^ 2 / 0.2e1 + (t145 * mrSges(4,2) + Ifges(4,4) * t137 + Ifges(4,5) * t154 + Ifges(4,1) * t134 / 0.2e1) * t134 + (-t145 * mrSges(8,1) + Ifges(8,4) * t136 + Ifges(8,6) * t153 + Ifges(8,2) * t133 / 0.2e1) * t133 + (t145 * mrSges(8,2) + Ifges(8,5) * t153 + Ifges(8,1) * t136 / 0.2e1) * t136 + (t120 * mrSges(11,2) - t116 * mrSges(11,3) + Ifges(11,5) * t148 + Ifges(11,1) * t115 / 0.2e1) * t115 + Ifges(8,3) * t153 ^ 2 / 0.2e1 + (-t126 * mrSges(6,1) + t112 * mrSges(6,3) + Ifges(6,4) * t119 + Ifges(6,6) * t121 + Ifges(6,2) * t118 / 0.2e1) * t118 + (Ifges(9,4) * t135 + Ifges(9,2) * t132 / 0.2e1) * t132 + (-t120 * mrSges(11,1) + t117 * mrSges(11,3) + Ifges(11,4) * t115 + Ifges(11,6) * t148 + Ifges(11,2) * t114 / 0.2e1) * t114 + (-t131 * mrSges(10,1) + Ifges(10,4) * t124 + Ifges(10,6) * t149 + Ifges(10,2) * t122 / 0.2e1) * t122 + Ifges(7,3) * qJD(6) ^ 2 / 0.2e1 + (t116 * mrSges(11,1) - t117 * mrSges(11,2) + Ifges(11,3) * t148 / 0.2e1) * t148 + m(11) * (t116 ^ 2 + t117 ^ 2 + t120 ^ 2) / 0.2e1 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t160 * (-mrSges(8,2) * t153 + mrSges(8,3) * t133) + t164 * (mrSges(4,1) * t154 - mrSges(4,3) * t134) + t168 * (mrSges(8,1) * t153 - mrSges(8,3) * t136) - t172 * (-mrSges(4,2) * t154 + mrSges(4,3) * t137)) * pkin(1)) * qJD(2) + (t129 * mrSges(5,1) - t130 * mrSges(5,2) + Ifges(5,3) * t150 / 0.2e1) * t150 + m(10) * (t131 ^ 2 + (t158 ^ 2 + t166 ^ 2) * pkin(2) ^ 2 * t152 ^ 2) / 0.2e1 + (t131 * mrSges(10,2) + Ifges(10,5) * t149 + Ifges(10,1) * t124 / 0.2e1) * t124 + (t126 * mrSges(6,2) - t111 * mrSges(6,3) + Ifges(6,5) * t121 + Ifges(6,1) * t119 / 0.2e1) * t119 + (-t128 * mrSges(5,1) + t130 * mrSges(5,3) + Ifges(5,4) * t125 + Ifges(5,6) * t150 + Ifges(5,2) * t123 / 0.2e1) * t123 + (-pkin(15) * (-mrSges(9,1) * t132 + t135 * mrSges(9,2)) + qJD(6) * (Ifges(7,5) * t161 + Ifges(7,6) * t169) + qJD(2) * (Ifges(3,5) * t173 - Ifges(3,6) * t165) + ((-pkin(15) * mrSges(3,2) + Ifges(3,1) * t173 / 0.2e1) * t173 + (-pkin(14) * mrSges(7,1) + Ifges(7,2) * t169 / 0.2e1) * t169 + (-pkin(15) * mrSges(3,1) - Ifges(3,4) * t173 + Ifges(3,2) * t165 / 0.2e1) * t165 + (pkin(14) * mrSges(7,2) + Ifges(7,4) * t169 + Ifges(7,1) * t161 / 0.2e1) * t161 + Ifges(2,3) / 0.2e1 + m(7) * pkin(14) ^ 2 / 0.2e1 + (m(3) / 0.2e1 + m(9) / 0.2e1) * pkin(15) ^ 2) * qJD(1)) * qJD(1) + (Ifges(9,5) * t135 + Ifges(9,6) * t132 + Ifges(9,3) * t152 / 0.2e1 + (-t158 * (-mrSges(10,2) * t149 + mrSges(10,3) * t122) - t166 * (mrSges(10,1) * t149 - mrSges(10,3) * t124)) * pkin(2)) * t152 + Ifges(9,1) * t135 ^ 2 / 0.2e1 + m(4) * (t144 + (t164 ^ 2 + t172 ^ 2) * t180) / 0.2e1 + m(8) * (t144 + (t160 ^ 2 + t168 ^ 2) * t180) / 0.2e1 + (t128 * mrSges(5,2) - t129 * mrSges(5,3) + Ifges(5,5) * t150 + Ifges(5,1) * t125 / 0.2e1) * t125 + (-t145 * mrSges(4,1) + Ifges(4,6) * t154 + Ifges(4,2) * t137 / 0.2e1) * t137 + m(6) * (t111 ^ 2 + t112 ^ 2 + t126 ^ 2) / 0.2e1 + m(5) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1;
T = t1;
