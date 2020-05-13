% Calculate potential energy for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE1_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:45
% EndTime: 2020-05-01 20:55:46
% DurationCPUTime: 0.59s
% Computational Cost: add. (289->115), mult. (499->165), div. (0->0), fcn. (312->24), ass. (0->57)
t106 = sin(qJ(4));
t112 = cos(qJ(4));
t133 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + (rSges(6,1) * t106 + rSges(6,2) * t112) * m(6) + rSges(7,3) * m(7);
t132 = m(5) + m(6);
t131 = m(4) * rSges(4,2);
t130 = pkin(13) * g(3);
t114 = cos(qJ(2));
t129 = t114 * pkin(1) + pkin(13);
t128 = m(3) + m(4) + m(5);
t108 = sin(qJ(2));
t127 = -pkin(1) * t108 + pkin(15);
t126 = (-rSges(6,3) - pkin(11)) * m(6) + m(5) * rSges(5,2);
t109 = sin(qJ(1));
t115 = cos(qJ(1));
t125 = g(1) * t115 + g(2) * t109;
t102 = cos(pkin(22));
t110 = sin(pkin(18));
t116 = cos(pkin(18));
t98 = sin(pkin(22));
t123 = t102 * t116 + t110 * t98;
t122 = t102 * t110 - t116 * t98;
t103 = cos(pkin(21));
t97 = qJ(3) + qJ(2);
t93 = sin(t97);
t94 = cos(t97);
t99 = sin(pkin(21));
t121 = rSges(11,1) * t94 - rSges(11,2) * t93 + (-t123 * t103 - t122 * t99) * pkin(4) + t127;
t120 = m(6) * (rSges(6,1) * t112 - rSges(6,2) * t106);
t107 = sin(qJ(3));
t113 = cos(qJ(3));
t111 = sin(pkin(17));
t117 = cos(pkin(17));
t89 = rSges(7,1) * t111 + rSges(7,2) * t117;
t90 = rSges(7,1) * t117 - rSges(7,2) * t111;
t91 = t132 * pkin(5) + m(4) * rSges(4,1);
t72 = t113 * t131 + rSges(3,1) * m(3) + t107 * t91 + (t110 * t89 + t116 * t90) * m(7) + (m(4) + t132) * pkin(1);
t73 = -t107 * t131 - rSges(3,2) * m(3) + t113 * t91 + (t110 * t90 - t116 * t89) * m(7);
t119 = -m(2) * rSges(2,1) + pkin(14) * m(7) + t108 * t72 - t114 * t73 + (-m(6) - t128) * pkin(15);
t105 = cos(pkin(19));
t104 = cos(pkin(20));
t101 = sin(pkin(19));
t100 = sin(pkin(20));
t96 = pkin(22) + pkin(21);
t92 = pkin(9) * m(6) + m(5) * rSges(5,1);
t88 = rSges(9,1) * t101 + rSges(9,2) * t105;
t87 = rSges(9,1) * t105 - rSges(9,2) * t101;
t82 = -rSges(10,2) + (-t101 * t107 + t105 * t113) * pkin(2);
t81 = rSges(10,1) + (t101 * t113 + t105 * t107) * pkin(2);
t80 = t100 * t92 - t126 * t104;
t79 = t100 * t126 + t92 * t104;
t78 = (-rSges(8,1) * t116 + rSges(8,2) * t110) * t102 + (-rSges(8,1) * t110 - rSges(8,2) * t116) * t98 + t127;
t77 = -t108 * t81 + t114 * t82 + pkin(15);
t76 = -g(3) * t87 + t125 * t88;
t75 = g(3) * t88 + t125 * t87;
t71 = t110 * t80 + t116 * t79 + (t100 * t110 + t104 * t116) * t120;
t70 = -t110 * t79 + t116 * t80 + (t100 * t116 - t104 * t110) * t120;
t1 = (-g(3) * t70 + t125 * t71) * cos(t96) + (-g(3) * t71 - t125 * t70) * sin(t96) + (t119 * g(1) + t133 * g(2)) * t115 + (-t133 * g(1) + t119 * g(2)) * t109 - t73 * g(3) * t108 - t72 * g(3) * t114 - m(6) * t130 - (m(2) * rSges(2,3) - pkin(16) * m(7) + (m(7) + m(2) + t128) * pkin(13)) * g(3) - m(8) * (g(1) * (rSges(8,3) * t109 + t115 * t78) + g(2) * (-rSges(8,3) * t115 + t109 * t78) + g(3) * (-(rSges(8,1) * t102 + rSges(8,2) * t98) * t110 + (rSges(8,1) * t98 - rSges(8,2) * t102) * t116 + t129)) - ((-t107 * t76 + t113 * t75) * t114 + (-rSges(9,3) * g(2) + pkin(15) * g(1)) * t115 + (rSges(9,3) * g(1) + g(2) * pkin(15)) * t109 + t130 + (-t107 * t75 - t113 * t76) * t108) * m(9) - m(10) * (g(1) * (rSges(10,3) * t109 + t115 * t77) + g(2) * (-rSges(10,3) * t115 + t109 * t77) + g(3) * (t108 * t82 + t114 * t81 + pkin(13))) - m(11) * ((-g(2) * rSges(11,3) + g(1) * t121) * t115 + (g(1) * rSges(11,3) + g(2) * t121) * t109 + (t93 * rSges(11,1) + t94 * rSges(11,2) + t129 + (-t122 * t103 + t123 * t99) * pkin(4)) * g(3)) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1);
U = t1;
