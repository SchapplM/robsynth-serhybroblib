% Calculate potential energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_energypot_fixb_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:23
% EndTime: 2020-05-07 04:33:24
% DurationCPUTime: 0.42s
% Computational Cost: add. (205->67), mult. (198->52), div. (0->0), fcn. (153->18), ass. (0->31)
t97 = pkin(3) * m(9);
t96 = -m(5) - m(6);
t74 = sin(qJ(5));
t78 = cos(qJ(5));
t95 = m(6) * pkin(8) + t78 * mrSges(6,1) - t74 * mrSges(6,2) + mrSges(5,1);
t94 = -m(9) - m(4) - m(8);
t93 = m(6) * pkin(10) - mrSges(5,2) + mrSges(6,3);
t71 = qJ(2) + qJ(7);
t67 = pkin(15) - t71;
t62 = -qJ(8) + t67;
t56 = sin(t62);
t57 = cos(t62);
t79 = cos(qJ(2));
t59 = t79 * pkin(1) + pkin(12);
t72 = qJ(2) + qJ(3);
t68 = qJ(4) + t72;
t60 = sin(t68);
t61 = cos(t68);
t63 = sin(t71);
t64 = sin(t72);
t65 = cos(t71);
t66 = cos(t72);
t73 = sin(qJ(6));
t75 = sin(qJ(2));
t77 = cos(qJ(6));
t92 = -cos(t67) * t97 + m(7) * pkin(6) - m(3) * pkin(12) - mrSges(3,1) * t79 + mrSges(4,1) * t66 - mrSges(7,1) * t77 - mrSges(8,1) * t65 + mrSges(9,1) * t57 + mrSges(3,2) * t75 - mrSges(4,2) * t64 + mrSges(7,2) * t73 + mrSges(8,2) * t63 + mrSges(9,2) * t56 - mrSges(2,1) + t95 * t61 + t96 * (-pkin(4) * t66 + t59) + t93 * t60 + t94 * t59;
t91 = -t74 * mrSges(6,1) - t78 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t58 = t75 * pkin(1) + pkin(11);
t80 = cos(qJ(1));
t76 = sin(qJ(1));
t1 = (sin(t67) * t97 - m(7) * pkin(13) - mrSges(3,1) * t75 + mrSges(4,1) * t64 - t73 * mrSges(7,1) - mrSges(8,1) * t63 - t56 * mrSges(9,1) - mrSges(3,2) * t79 + mrSges(4,2) * t66 - t77 * mrSges(7,2) - mrSges(8,2) * t65 + t57 * mrSges(9,2) - mrSges(1,3) - mrSges(2,3) + t96 * (-pkin(4) * t64 + t58) - t93 * t61 + t95 * t60 + t94 * t58 + (-m(7) - m(2) - m(3)) * pkin(11)) * g(3) + (t92 * t76 - t91 * t80 - mrSges(1,2)) * g(2) + (t91 * t76 + t92 * t80 - mrSges(1,1)) * g(1);
U = t1;
