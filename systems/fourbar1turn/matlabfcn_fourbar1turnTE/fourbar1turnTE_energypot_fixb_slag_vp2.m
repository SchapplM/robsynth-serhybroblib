% Calculate potential energy for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnTE_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:28
% EndTime: 2020-04-12 19:18:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (424->49), mult. (608->60), div. (36->3), fcn. (195->6), ass. (0->30)
t60 = cos(qJ(2));
t81 = pkin(2) * t60;
t77 = (-0.2e1 * t81 + pkin(1)) * pkin(1);
t74 = pkin(2) ^ 2 + t77;
t76 = pkin(3) ^ 2 - pkin(4) ^ 2;
t52 = t74 + t76;
t56 = pkin(1) * t60 - pkin(2);
t83 = -pkin(3) + pkin(4);
t84 = -pkin(3) - pkin(4);
t51 = sqrt(-((pkin(2) - t84) * (pkin(2) + t84) + t77) * ((pkin(2) - t83) * (pkin(2) + t83) + t77));
t58 = sin(qJ(2));
t80 = t51 * t58;
t47 = -pkin(1) * t80 - t52 * t56;
t50 = pkin(1) * t52 * t58 - t51 * t56;
t54 = 0.1e1 / t74;
t78 = t54 / pkin(3);
t85 = t58 / 0.2e1;
t89 = (t60 * t50 / 0.2e1 + t47 * t85) * t78;
t88 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t53 = t74 - t76;
t55 = pkin(1) - t81;
t82 = pkin(2) * t58;
t49 = t51 * t55 + t53 * t82;
t71 = (t50 * t85 - t60 * t47 / 0.2e1) * t78;
t79 = t54 / pkin(4);
t86 = pkin(2) * t80 / 0.2e1 - t53 * t55 / 0.2e1;
t87 = -mrSges(2,1) - m(5) * pkin(1) - (mrSges(5,1) * t86 - t49 * mrSges(5,2) / 0.2e1) * t79 - mrSges(4,2) * t89 - mrSges(4,1) * t71 - mrSges(3,1) * t60 + mrSges(3,2) * t58 - m(4) * t81;
t61 = cos(qJ(1));
t59 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) - t58 * mrSges(3,1) - t60 * mrSges(3,2) - m(4) * t82 + mrSges(4,1) * t89 - mrSges(4,2) * t71 - (t49 * mrSges(5,1) / 0.2e1 + mrSges(5,2) * t86) * t79 + (-m(2) - m(3) - m(4) - m(5)) * pkin(5)) * g(3) + (t87 * t59 - t88 * t61 - mrSges(1,2)) * g(2) + (t88 * t59 + t87 * t61 - mrSges(1,1)) * g(1);
U = t1;
