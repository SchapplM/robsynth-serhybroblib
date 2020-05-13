% Calculate potential energy for
% palh1m2OL
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:20
% EndTime: 2020-05-02 21:17:21
% DurationCPUTime: 0.38s
% Computational Cost: add. (263->80), mult. (240->72), div. (0->0), fcn. (134->24), ass. (0->36)
t88 = m(5) + m(6);
t92 = m(8) + m(4) + t88;
t60 = -qJ(7) + pkin(19) - qJ(2);
t58 = -qJ(10) + t60;
t56 = sin(t58);
t57 = cos(t58);
t66 = sin(qJ(6));
t70 = sin(qJ(2));
t74 = cos(qJ(6));
t85 = -m(9) - m(3) - t92;
t91 = pkin(14) * m(7) - t74 * mrSges(7,1) + mrSges(7,2) * t66 - mrSges(2,1) - m(11) * (-t70 * pkin(1) + pkin(4) * sin(t60)) + mrSges(11,1) * t56 - mrSges(11,2) * t57 + (t85 - m(11) - m(10)) * pkin(15);
t67 = sin(qJ(5));
t75 = cos(qJ(5));
t90 = mrSges(6,1) * t67 + mrSges(6,2) * t75 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t89 = m(10) * g(3);
t63 = qJ(2) + qJ(8);
t71 = sin(qJ(1));
t79 = cos(qJ(1));
t84 = g(1) * t79 + t71 * g(2);
t54 = m(6) * pkin(9) + mrSges(6,1) * t75 - mrSges(6,2) * t67 + mrSges(5,1);
t59 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t68 = sin(qJ(4));
t76 = cos(qJ(4));
t82 = t54 * t68 - t59 * t76 + mrSges(4,2);
t78 = cos(qJ(2));
t77 = cos(qJ(3));
t73 = cos(qJ(7));
t72 = cos(qJ(8));
t69 = sin(qJ(3));
t65 = sin(qJ(7));
t64 = sin(qJ(8));
t61 = qJ(9) + t63;
t52 = t88 * pkin(5) + t54 * t76 + t59 * t68 + mrSges(4,1);
t51 = -t65 * mrSges(8,1) - t64 * mrSges(9,1) - t73 * mrSges(8,2) - t72 * mrSges(9,2) + t52 * t77 - t82 * t69 - mrSges(3,2);
t50 = t92 * pkin(1) + t73 * mrSges(8,1) + t72 * mrSges(9,1) - mrSges(8,2) * t65 - t64 * mrSges(9,2) + t52 * t69 + t82 * t77 + mrSges(3,1);
t1 = (-t51 * g(3) + t84 * t50) * t70 + (-t50 * g(3) - t84 * t51) * t78 - mrSges(7,2) * t74 * g(3) - mrSges(7,1) * t66 * g(3) - (-m(7) * pkin(16) + mrSges(1,3) + mrSges(2,3) + (m(7) + m(2) - t85) * pkin(13)) * g(3) + (g(3) * mrSges(10,1) - t84 * mrSges(10,2)) * cos(t61) + (-t84 * mrSges(10,1) - g(3) * mrSges(10,2)) * sin(t61) - pkin(13) * t89 - g(3) * (m(11) * (pkin(4) * cos(t60) + t78 * pkin(1) + pkin(13)) - t57 * mrSges(11,1) - t56 * mrSges(11,2)) + (t91 * t71 + t90 * t79 - mrSges(1,2)) * g(2) + (-t90 * t71 + t91 * t79 - mrSges(1,1)) * g(1) + (m(10) * t84 * sin(t63) - cos(t63) * t89) * pkin(2);
U = t1;
