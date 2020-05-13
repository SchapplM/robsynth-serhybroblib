% Calculate potential energy for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:18
% EndTime: 2020-04-15 19:28:19
% DurationCPUTime: 0.58s
% Computational Cost: add. (251->81), mult. (241->61), div. (0->0), fcn. (184->22), ass. (0->39)
t120 = -m(5) - m(6);
t119 = pkin(4) * m(11);
t93 = sin(qJ(5));
t97 = cos(qJ(5));
t118 = -m(6) * pkin(9) - t97 * mrSges(6,1) + t93 * mrSges(6,2) - mrSges(5,1);
t117 = -m(10) - m(9) - m(3);
t116 = -m(11) - m(4) - m(8);
t115 = m(10) * pkin(2) + mrSges(9,1);
t114 = m(6) * pkin(11) - mrSges(5,2) + mrSges(6,3);
t113 = -t93 * mrSges(6,1) - t97 * mrSges(6,2) + mrSges(2,2) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3);
t90 = qJ(2) + qJ(7);
t79 = pkin(19) - t90;
t71 = -qJ(10) + t79;
t69 = sin(t71);
t70 = cos(t71);
t94 = sin(qJ(2));
t73 = -t94 * pkin(1) + pkin(15);
t89 = qJ(2) + qJ(8);
t86 = qJ(9) + t89;
t75 = sin(t86);
t91 = qJ(2) + qJ(3);
t87 = qJ(4) + t91;
t76 = sin(t87);
t77 = cos(t86);
t78 = cos(t87);
t80 = sin(t89);
t81 = sin(t90);
t82 = sin(t91);
t83 = cos(t89);
t84 = cos(t90);
t85 = cos(t91);
t92 = sin(qJ(6));
t96 = cos(qJ(6));
t98 = cos(qJ(2));
t112 = -mrSges(7,1) * t96 + mrSges(7,2) * t92 + mrSges(11,1) * t69 - mrSges(11,2) * t70 - mrSges(10,1) * t75 - mrSges(10,2) * t77 + mrSges(9,2) * t83 + mrSges(8,1) * t81 + mrSges(8,2) * t84 - mrSges(4,1) * t85 + mrSges(4,2) * t82 + mrSges(3,1) * t94 + mrSges(3,2) * t98 - sin(t79) * t119 + m(7) * pkin(14) - mrSges(2,1) + t115 * t80 + t118 * t78 + t120 * (pkin(5) * t85 + t73) - t114 * t76 + t116 * t73 + t117 * pkin(15);
t74 = t98 * pkin(1) + pkin(13);
t99 = cos(qJ(1));
t95 = sin(qJ(1));
t1 = (-cos(t79) * t119 + t94 * mrSges(3,2) - t96 * mrSges(7,2) - t98 * mrSges(3,1) + t80 * mrSges(9,2) + t81 * mrSges(8,2) - t82 * mrSges(4,1) - t84 * mrSges(8,1) - t85 * mrSges(4,2) - t92 * mrSges(7,1) + t69 * mrSges(11,2) + t70 * mrSges(11,1) - t75 * mrSges(10,2) + t77 * mrSges(10,1) + m(7) * pkin(16) - mrSges(1,3) - mrSges(2,3) - t115 * t83 + t120 * (pkin(5) * t82 + t74) + t114 * t78 + t118 * t76 + t116 * t74 + (-m(7) - m(2) + t117) * pkin(13)) * g(3) + (t112 * t95 - t113 * t99 - mrSges(1,2)) * g(2) + (t112 * t99 + t113 * t95 - mrSges(1,1)) * g(1);
U = t1;
