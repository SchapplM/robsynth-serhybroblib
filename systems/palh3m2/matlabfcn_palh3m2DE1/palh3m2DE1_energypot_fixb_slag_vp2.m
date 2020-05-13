% Calculate potential energy for
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2DE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_energypot_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_energypot_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE1_energypot_fixb_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:57:58
% EndTime: 2020-05-07 01:57:59
% DurationCPUTime: 0.31s
% Computational Cost: add. (237->83), mult. (319->93), div. (0->0), fcn. (239->22), ass. (0->45)
t89 = sin(qJ(4));
t95 = cos(qJ(4));
t114 = mrSges(6,1) * t89 + mrSges(6,2) * t95 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3);
t85 = sin(pkin(18));
t87 = cos(pkin(18));
t93 = sin(pkin(15));
t99 = cos(pkin(15));
t104 = -t85 * t93 + t87 * t99;
t105 = t85 * t99 + t87 * t93;
t112 = m(5) + m(6);
t110 = m(4) + t112;
t109 = -m(3) - t110;
t100 = cos(pkin(14));
t94 = sin(pkin(14));
t70 = mrSges(7,1) * t94 - mrSges(7,2) * t100;
t71 = mrSges(7,1) * t100 + mrSges(7,2) * t94;
t72 = pkin(4) * t112 + mrSges(4,1);
t90 = sin(qJ(3));
t96 = cos(qJ(3));
t61 = pkin(1) * t110 + t90 * mrSges(4,2) + t70 * t93 + t71 * t99 - t72 * t96 + mrSges(3,1);
t63 = t96 * mrSges(4,2) + t70 * t99 - t71 * t93 + t72 * t90 - mrSges(3,2);
t97 = cos(qJ(2));
t74 = pkin(1) * t97 + pkin(12);
t82 = qJ(3) + qJ(2);
t76 = sin(t82);
t77 = cos(t82);
t86 = sin(pkin(17));
t88 = cos(pkin(17));
t91 = sin(qJ(2));
t113 = -m(9) * ((-t104 * t88 + t105 * t86) * pkin(3) + t74) + mrSges(9,1) * t77 - mrSges(9,2) * t76 + m(7) * pkin(6) + pkin(12) * t109 - t61 * t97 - t63 * t91 - mrSges(2,1) - m(8) * t74 - (-mrSges(8,1) * t99 - mrSges(8,2) * t93) * t87 - t85 * (mrSges(8,1) * t93 - mrSges(8,2) * t99);
t111 = t91 * pkin(1) + pkin(11);
t92 = sin(qJ(1));
t98 = cos(qJ(1));
t107 = g(1) * t98 + g(2) * t92;
t106 = mrSges(6,1) * t95 - mrSges(6,2) * t89;
t84 = cos(pkin(16));
t83 = sin(pkin(16));
t81 = pkin(17) + pkin(18);
t78 = pkin(8) * m(6) + mrSges(5,1);
t75 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t66 = -t83 * t75 + t78 * t84;
t65 = t75 * t84 + t83 * t78;
t60 = -t93 * t65 + t66 * t99 + t106 * (-t83 * t93 + t84 * t99);
t59 = t65 * t99 + t93 * t66 + t106 * (t83 * t99 + t84 * t93);
t1 = (-g(3) * t59 + t107 * t60) * cos(t81) + (-g(3) * t60 - t107 * t59) * sin(t81) + g(3) * t63 * t97 - g(3) * t61 * t91 - ((m(2) + m(7) - t109) * pkin(11) + m(7) * pkin(13) + mrSges(1,3) + mrSges(2,3)) * g(3) - g(3) * (m(8) * t111 + (t85 * mrSges(8,1) - mrSges(8,2) * t87) * t99 + t93 * (mrSges(8,1) * t87 + mrSges(8,2) * t85)) - g(3) * (-t76 * mrSges(9,1) - t77 * mrSges(9,2) + (t111 + (t104 * t86 + t105 * t88) * pkin(3)) * m(9)) + (t113 * t92 + t114 * t98 - mrSges(1,2)) * g(2) + (t113 * t98 - t114 * t92 - mrSges(1,1)) * g(1);
U = t1;
