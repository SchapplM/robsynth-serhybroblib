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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:44
% EndTime: 2020-05-01 20:55:45
% DurationCPUTime: 0.19s
% Computational Cost: add. (371->76), mult. (412->73), div. (0->0), fcn. (308->20), ass. (0->44)
t108 = m(5) + m(6);
t99 = m(11) + m(4) + m(8) + t108;
t107 = -m(3) - m(9) - m(10) - t99;
t84 = sin(qJ(4));
t89 = cos(qJ(4));
t106 = -mrSges(6,1) * t89 + mrSges(6,2) * t84;
t72 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t74 = pkin(9) * m(6) + mrSges(5,1);
t78 = sin(pkin(20));
t82 = cos(pkin(20));
t97 = t72 * t82 + t78 * t74;
t96 = mrSges(6,1) * t84 + mrSges(6,2) * t89 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t75 = pkin(2) * m(10) + mrSges(9,1);
t79 = sin(pkin(19));
t83 = cos(pkin(19));
t95 = mrSges(9,2) * t83 + t75 * t79 + mrSges(11,2) + mrSges(4,2);
t66 = m(11) * pkin(4) - t72 * t78 + t74 * t82;
t77 = sin(pkin(21));
t81 = cos(pkin(21));
t61 = t66 * t81 - t97 * t77 + mrSges(8,1);
t62 = t66 * t77 + t97 * t81 - mrSges(8,2);
t76 = sin(pkin(22));
t80 = cos(pkin(22));
t56 = t61 * t80 - t62 * t76;
t57 = t61 * t76 + t62 * t80;
t65 = t108 * pkin(5) - t79 * mrSges(9,2) + t75 * t83 + mrSges(11,1) + mrSges(4,1);
t88 = sin(pkin(17));
t93 = cos(pkin(17));
t70 = mrSges(7,1) * t88 + mrSges(7,2) * t93;
t71 = mrSges(7,1) * t93 - mrSges(7,2) * t88;
t85 = sin(qJ(3));
t87 = sin(pkin(18));
t90 = cos(qJ(3));
t92 = cos(pkin(18));
t58 = t99 * pkin(1) + t65 * t85 + t70 * t87 + t71 * t92 + t95 * t90 + mrSges(3,1) + mrSges(10,1);
t60 = t65 * t90 - t70 * t92 + t71 * t87 - t95 * t85 - mrSges(3,2) - mrSges(10,2);
t67 = t82 * t77 + t81 * t78;
t68 = -t78 * t77 + t82 * t81;
t63 = t67 * t80 + t76 * t68;
t64 = -t76 * t67 + t68 * t80;
t86 = sin(qJ(2));
t91 = cos(qJ(2));
t94 = pkin(14) * m(7) + t56 * t92 + t57 * t87 + t58 * t86 - t60 * t91 - mrSges(2,1) + t107 * pkin(15) - t106 * (t63 * t87 + t64 * t92);
t1 = (t94 * g(1) + t96 * g(2)) * cos(qJ(1)) + (-t96 * g(1) + t94 * g(2)) * sin(qJ(1)) - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (t106 * (t63 * t92 - t87 * t64) - t60 * t86 - t58 * t91 - t57 * t92 + t56 * t87 + m(7) * pkin(16) - mrSges(1,3) - mrSges(2,3) - (m(7) + m(2) - t107) * pkin(13)) * g(3);
U = t1;
