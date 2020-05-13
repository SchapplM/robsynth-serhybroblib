% Calculate potential energy for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2DE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_energypot_fixb_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:27
% EndTime: 2020-05-07 02:14:28
% DurationCPUTime: 1.22s
% Computational Cost: add. (668->143), mult. (350->144), div. (0->0), fcn. (182->68), ass. (0->76)
t185 = 2 * g(1);
t188 = m(5) + m(6);
t157 = m(9) + m(4) + m(8) + t188;
t187 = m(3) + t157;
t118 = (pkin(1) * t157 + mrSges(3,1));
t186 = (m(7) * pkin(6) - pkin(12) * t187 - mrSges(2,1));
t184 = pkin(3) * m(9);
t181 = mrSges(6,1) * g(1);
t180 = mrSges(6,1) * g(2);
t179 = mrSges(7,1) * g(1);
t178 = mrSges(7,1) * g(2);
t177 = mrSges(8,1) * g(1);
t176 = mrSges(8,1) * g(2);
t175 = mrSges(9,1) * g(2);
t174 = mrSges(3,2) * g(1);
t173 = mrSges(3,2) * g(2);
t172 = mrSges(4,2) * g(1);
t171 = mrSges(4,2) * g(2);
t170 = mrSges(6,2) * g(1);
t168 = mrSges(7,2) * g(1);
t167 = mrSges(7,2) * g(2);
t166 = mrSges(8,2) * g(2);
t165 = mrSges(9,2) * g(2);
t160 = g(2) * t118;
t132 = (pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3));
t159 = g(2) * t132;
t137 = (pkin(8) * m(6) + mrSges(5,1));
t158 = g(2) * t137;
t156 = pkin(15) + qJ(2);
t140 = qJ(2) + qJ(3);
t133 = pkin(18) + t156;
t114 = qJ(2) + atan2(sin(t133), -cos(t133));
t134 = pkin(14) - t156;
t124 = pkin(17) + qJ(3) + t133;
t119 = pkin(16) + t124;
t107 = atan2(-sin(t119), cos(t119)) + t140;
t151 = pkin(17) - t114;
t104 = -qJ(4) + t107;
t103 = qJ(4) + t107;
t98 = -atan2(-sin(t124), -cos(t124)) + t151;
t148 = mrSges(9,1) * g(1);
t147 = mrSges(6,2) * g(2);
t146 = mrSges(8,2) * g(1);
t145 = mrSges(9,2) * g(1);
t144 = qJ(1) - qJ(2);
t143 = qJ(1) + qJ(2);
t142 = qJ(1) - qJ(4);
t141 = qJ(1) + qJ(4);
t136 = qJ(1) - t140;
t135 = qJ(1) + t140;
t131 = -qJ(1) + t134;
t130 = qJ(1) + t134;
t129 = -t170 + t180;
t128 = t170 + t180;
t127 = -t147 + t181;
t126 = t147 + t181;
t125 = (t188 * pkin(4) + mrSges(4,1));
t123 = t137 * t185;
t122 = (mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3));
t121 = g(2) * t125;
t120 = t132 * t185;
t117 = t125 * t185;
t116 = t118 * t185;
t113 = qJ(1) - t114;
t112 = qJ(1) + t114;
t111 = -qJ(1) + t151;
t110 = qJ(1) + t151;
t106 = qJ(1) - t107;
t105 = qJ(1) + t107;
t102 = qJ(1) - t103;
t101 = qJ(1) + t104;
t100 = qJ(1) - t104;
t99 = qJ(1) + t103;
t97 = -qJ(1) + t98;
t96 = qJ(1) + t98;
t1 = (sin(t100) + sin(t99)) * t129 / 0.4e1 + (sin(t102) + sin(t101)) * t128 / 0.4e1 + (cos(t102) + cos(t101)) * t127 / 0.4e1 + (cos(t100) + cos(t99)) * t126 / 0.4e1 + (-mrSges(7,2) * cos(t134) + mrSges(7,1) * sin(t134) + t137 * sin(t107) + t125 * sin(t140) + mrSges(4,2) * cos(t140) - mrSges(8,1) * sin(t114) - mrSges(8,2) * cos(t114) - t118 * sin(qJ(2)) - t132 * cos(t107) - mrSges(9,1) * sin(t98) + mrSges(9,2) * cos(t98) - pkin(13) * m(7) - mrSges(1,3) - mrSges(2,3) - (m(7) + m(2) + t187) * pkin(11) + sin(t151) * t184 - mrSges(3,2) * cos(qJ(2)) + (cos(t103) / 0.2e1 - cos(t104) / 0.2e1) * mrSges(6,2) + (sin(t104) + sin(t103)) * mrSges(6,1) / 0.2e1) * g(3) - (g(2) * sin(t110) + (cos(t111) + cos(t110)) * g(1)) * t184 / 0.2e1 + (-t116 + 2 * t173) * cos(t144) / 0.4e1 + (-t116 - 2 * t173) * cos(t143) / 0.4e1 + (t186 * g(1) - g(2) * t122) * cos(qJ(1)) + (-t160 - t174) * sin(t144) / 0.2e1 + (-t160 + t174) * sin(t143) / 0.2e1 + (t145 - t175) * sin(t97) / 0.2e1 + (t145 + t175) * sin(t96) / 0.2e1 + (t146 - t176) * sin(t112) / 0.2e1 - (t146 + t176) * sin(t113) / 0.2e1 + (-t166 - t177) * cos(t112) / 0.2e1 + (t166 - t177) * cos(t113) / 0.2e1 + (-t168 - t178) * sin(t130) / 0.2e1 + (-t168 + t178) * sin(t131) / 0.2e1 + (-t167 - t179) * cos(t131) / 0.2e1 + (t167 - t179) * cos(t130) / 0.2e1 + (t148 - t165) * cos(t96) / 0.2e1 + (t148 + t165) * cos(t97) / 0.2e1 + (t117 + 2 * t171) * cos(t135) / 0.4e1 + (t117 - 2 * t171) * cos(t136) / 0.4e1 + (t121 - t172) * sin(t135) / 0.2e1 + (t121 + t172) * sin(t136) / 0.2e1 + (t120 + 2 * t158) * sin(t105) / 0.4e1 - (t120 - 2 * t158) * sin(t106) / 0.4e1 + (t123 + 2 * t159) * cos(t106) / 0.4e1 + (t123 - 2 * t159) * cos(t105) / 0.4e1 + t126 * cos(t141) / 0.2e1 + t129 * sin(t141) / 0.2e1 - t127 * cos(t142) / 0.2e1 - t128 * sin(t142) / 0.2e1 + g(2) * sin(t111) * t184 / 0.2e1 + (t122 * g(1) + t186 * g(2)) * sin(qJ(1)) - (mrSges(1,1) * g(1)) - (mrSges(1,2) * g(2));
U = t1;
