% Calculate potential energy for
% palh1m2DE2
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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:51
% EndTime: 2020-05-02 20:57:53
% DurationCPUTime: 0.70s
% Computational Cost: add. (359->108), mult. (276->105), div. (0->0), fcn. (219->45), ass. (0->57)
t100 = sin(qJ(4));
t106 = cos(qJ(4));
t135 = mrSges(6,1) * t100 + mrSges(6,2) * t106 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t102 = sin(qJ(2));
t108 = cos(qJ(2));
t132 = m(5) + m(6);
t133 = m(8) + m(4) + t132;
t123 = -m(9) - m(10) - m(3) - t133;
t101 = sin(qJ(3));
t105 = sin(pkin(17));
t107 = cos(qJ(3));
t111 = cos(pkin(17));
t104 = sin(pkin(18));
t110 = cos(pkin(18));
t114 = mrSges(7,1) * t104 - t110 * mrSges(7,2);
t115 = t110 * mrSges(7,1) + mrSges(7,2) * t104;
t88 = t132 * pkin(5) + mrSges(4,1);
t72 = t133 * pkin(1) + mrSges(4,2) * t107 + t88 * t101 + t114 * t105 + t115 * t111 + mrSges(3,1);
t73 = -t101 * mrSges(4,2) - t115 * t105 + t88 * t107 + t114 * t111 - mrSges(3,2);
t134 = pkin(14) * m(7) + t72 * t102 - t73 * t108 - mrSges(2,1) + (-m(11) + t123) * pkin(15);
t130 = m(11) * g(1);
t129 = m(11) * g(2);
t128 = g(3) * m(11);
t126 = mrSges(11,1) * g(2);
t125 = mrSges(11,2) * g(1);
t124 = mrSges(11,2) * g(2);
t97 = sin(pkin(19));
t98 = cos(pkin(19));
t83 = t101 * t98 + t107 * t97;
t84 = t101 * t97 - t107 * t98;
t77 = qJ(2) + atan2(t84, t83);
t94 = pkin(18) - pkin(22);
t122 = -t130 / 0.2e1;
t121 = -t129 / 0.2e1;
t120 = t129 / 0.2e1;
t92 = -qJ(2) + t94;
t119 = pkin(21) - atan2(cos(t92), -sin(t92)) - qJ(2);
t117 = -pkin(21) + t94;
t81 = -qJ(1) + t119;
t80 = qJ(1) + t119;
t103 = sin(qJ(1));
t109 = cos(qJ(1));
t116 = g(1) * t109 + g(2) * t103;
t99 = mrSges(11,1) * g(1);
t96 = qJ(1) - qJ(2);
t95 = qJ(1) + qJ(2);
t93 = pkin(2) * m(10) + mrSges(9,1);
t90 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t89 = -pkin(20) + t117;
t87 = -qJ(2) - qJ(3) + t117;
t85 = pkin(9) * m(6) + mrSges(6,1) * t106 - mrSges(6,2) * t100 + mrSges(5,1);
t79 = atan2(-sin(t87), cos(t87));
t76 = -t79 + t119;
t75 = -t79 + t81;
t74 = -t79 + t80;
t71 = atan2(t84, -t83) + t77;
t1 = (t99 - t124) * sin(t74) / 0.2e1 + (-t125 - t126) * cos(t74) / 0.2e1 + (t99 + t124) * sin(t75) / 0.2e1 + (-t125 + t126) * cos(t75) / 0.2e1 + (-g(3) * t90 + t116 * t85) * cos(t89) + (t116 * mrSges(9,2) - g(3) * t93) * cos(t77) + (g(3) * t85 + t116 * t90) * sin(t89) + (mrSges(9,2) * g(3) + t116 * t93) * sin(t77) + (t116 * mrSges(8,1) + g(3) * mrSges(8,2)) * cos(t94) + (-t116 * mrSges(10,1) - g(3) * mrSges(10,2)) * sin(t71) + (g(3) * mrSges(10,1) - t116 * mrSges(10,2)) * cos(t71) + (g(3) * mrSges(8,1) - t116 * mrSges(8,2)) * sin(t94) - pkin(13) * t128 + (-t108 * t128 + sin(t95) * t130 / 0.2e1 + sin(t96) * t122 + cos(t96) * t120 + cos(t95) * t121) * pkin(1) - (-pkin(16) * m(7) + mrSges(1,3) + mrSges(2,3) + (m(7) + m(2) - t123) * pkin(13)) * g(3) - g(3) * t73 * t102 - t72 * g(3) * t108 + mrSges(11,1) * cos(t76) * g(3) + mrSges(11,2) * sin(t76) * g(3) + (-t135 * g(1) + t134 * g(2)) * t103 + (t134 * g(1) + t135 * g(2)) * t109 - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (-cos(t119) * t128 + cos(t80) * t120 + cos(t81) * t121 + (sin(t81) + sin(t80)) * t122) * pkin(4);
U = t1;
