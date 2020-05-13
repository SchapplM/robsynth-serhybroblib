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
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2DE2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energypot_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_energypot_fixb_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:03
% EndTime: 2020-05-07 02:14:04
% DurationCPUTime: 1.23s
% Computational Cost: add. (422->119), mult. (375->156), div. (0->0), fcn. (220->52), ass. (0->65)
t168 = (m(5) + m(6));
t176 = (m(4) + m(8) + m(9) + t168);
t174 = -(2 * rSges(2,2) * m(2)) + (2 * rSges(3,3) * m(3)) + (2 * m(4) * rSges(4,3)) + (2 * m(5) * rSges(5,3)) + (2 * rSges(7,3) * m(7)) + (2 * m(8) * rSges(8,3)) + (2 * m(9) * rSges(9,3)) + 0.2e1 * (rSges(6,1) * sin(qJ(4)) + rSges(6,2) * cos(qJ(4))) * m(6);
t169 = m(9) * g(3);
t167 = m(4) * rSges(4,2);
t166 = rSges(7,1) * g(2);
t165 = rSges(3,2) * m(3);
t164 = rSges(7,2) * g(2);
t134 = sin(qJ(1));
t139 = cos(qJ(1));
t108 = t139 * g(1) + t134 * g(2);
t163 = m(6) * t108;
t161 = m(9) * t108;
t133 = sin(qJ(2));
t160 = rSges(6,2) * t133;
t138 = cos(qJ(2));
t159 = rSges(6,2) * t138;
t126 = pkin(15) + pkin(18);
t120 = qJ(2) + t126;
t112 = pkin(17) + qJ(3) + t120;
t109 = pkin(16) + t112;
t98 = atan2(-sin(t109), cos(t109));
t97 = qJ(3) + t98;
t158 = pkin(14) - qJ(2);
t127 = qJ(2) + qJ(3);
t99 = pkin(17) - atan2(-sin(t112), -cos(t112));
t153 = t163 / 0.2e1;
t95 = -qJ(4) + t97;
t94 = qJ(4) + t97;
t151 = -m(3) - t176;
t103 = pkin(1) * t176 + rSges(3,1) * m(3);
t149 = (2 * m(7) * pkin(6)) - (2 * rSges(2,1) * m(2)) + (2 * t151 * pkin(12)) - 0.2e1 * t103 * t138 + 0.2e1 * t133 * t165;
t142 = rSges(7,1) * g(1);
t141 = rSges(7,2) * g(1);
t140 = cos(pkin(15));
t137 = cos(qJ(3));
t135 = sin(pkin(15));
t132 = sin(qJ(3));
t130 = cos(pkin(18));
t129 = sin(pkin(18));
t128 = -pkin(15) + pkin(14);
t125 = pkin(17) - qJ(2);
t124 = rSges(6,1) * t138;
t123 = rSges(6,1) * t133;
t122 = -qJ(1) + t158;
t121 = qJ(1) + t158;
t119 = pkin(8) * m(6) + (m(5) * rSges(5,1));
t118 = t141 - t166;
t117 = t141 + t166;
t116 = t142 - t164;
t115 = t142 + t164;
t114 = cos(t120);
t113 = sin(t120);
t111 = (-rSges(6,3) - pkin(10)) * m(6) + (m(5) * rSges(5,2));
t110 = t168 * pkin(4) + m(4) * rSges(4,1);
t107 = t138 * rSges(9,1) - t133 * rSges(9,2);
t106 = rSges(9,1) * t133 + t138 * rSges(9,2);
t105 = -t135 * t129 + t140 * t130;
t104 = t140 * t129 + t135 * t130;
t96 = t98 + t127;
t93 = qJ(2) + t95;
t92 = qJ(2) + t94;
t91 = -atan2(t113, -t114) + t99;
t90 = pkin(17) - atan2(t138 * t104 + t133 * t105, t133 * t104 - t138 * t105);
t1 = -m(1) * (rSges(1,1) * g(1) + rSges(1,2) * g(2)) + (-t174 * g(1) + g(2) * t149) * t134 / 0.2e1 + (g(1) * t149 + t174 * g(2)) * t139 / 0.2e1 + ((-t106 * t114 + t113 * t107) * cos(t99) + (t113 * t106 + t107 * t114) * sin(t99)) * t169 + (t106 * sin(t91) + t107 * cos(t91)) * t161 + (cos(t92) + cos(t93)) * rSges(6,1) * t153 + (-sin(t92) * t163 / 0.2e1 + sin(t93) * t153) * rSges(6,2) - ((t135 * t116 + t117 * t140) * sin(t121) + (t115 * t140 - t135 * t118) * cos(t122) + (t116 * t140 - t135 * t117) * cos(t121) + (t135 * t115 + t118 * t140) * sin(t122)) * m(7) / 0.2e1 + ((rSges(8,2) * sin(t126) + rSges(8,1) * cos(t126)) * m(8) - sin(t127) * t167 + t119 * cos(t96) + t110 * cos(t127) - t111 * sin(t96)) * t108 + ((-cos(t90) * t133 + sin(t90) * t138) * t169 + (t114 * cos(t125) - t113 * sin(t125)) * t161) * pkin(3) + (-rSges(1,3) * m(1) - (m(2) * rSges(2,3)) + ((-m(2) + t151) * pkin(11)) + ((t123 - t159) * cos(t95) + (t123 + t159) * cos(t94) + (t124 - t160) * sin(t94) + (t124 + t160) * sin(t95)) * m(6) / 0.2e1 + (-(t138 * rSges(8,1) - t133 * rSges(8,2)) * t113 + (rSges(8,1) * t133 + t138 * rSges(8,2)) * t114) * m(8) - (-t110 * t137 + t132 * t167 + t103) * t133 + (t110 * t132 + t137 * t167 - t165) * t138 + (t111 * t138 + t133 * t119) * cos(t97) + (-t133 * t111 + t119 * t138) * sin(t97) + (-pkin(13) - pkin(11) - (rSges(7,1) * t133 + rSges(7,2) * t138) * cos(t128) + (rSges(7,1) * t138 - rSges(7,2) * t133) * sin(t128)) * m(7)) * g(3);
U = t1;
