% Calculate potential energy for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:25
% EndTime: 2020-04-12 19:25:25
% DurationCPUTime: 0.32s
% Computational Cost: add. (1198->51), mult. (1742->64), div. (108->6), fcn. (519->10), ass. (0->32)
t93 = pkin(4) ^ 2;
t92 = pkin(3) ^ 2;
t66 = cos(qJ(2));
t86 = pkin(2) * t66;
t84 = (-0.2e1 * t86 + pkin(1)) * pkin(1);
t60 = pkin(2) ^ 2 + t84;
t83 = t92 - t93;
t57 = t60 - t83;
t61 = pkin(1) - t86;
t88 = -pkin(3) + pkin(4);
t89 = -pkin(3) - pkin(4);
t55 = sqrt(-((pkin(2) - t89) * (pkin(2) + t89) + t84) * ((pkin(2) - t88) * (pkin(2) + t88) + t84));
t64 = sin(qJ(2));
t85 = t55 * t64;
t52 = -pkin(2) * t85 + t57 * t61;
t87 = pkin(2) * t64;
t53 = t55 * t61 + t57 * t87;
t56 = t60 + t83;
t62 = pkin(1) * t66 - pkin(2);
t51 = -pkin(1) * t85 - t56 * t62;
t54 = pkin(1) * t56 * t64 - t55 * t62;
t78 = t51 * t64 + t54 * t66;
t79 = -t51 * t66 + t54 * t64;
t58 = 0.1e1 / t60;
t59 = 0.1e1 / t60 ^ 2;
t81 = ((t51 ^ 2 + t54 ^ 2) / t92 * t59) ^ (-0.1e1 / 0.2e1) * t58 / pkin(3);
t82 = ((t52 ^ 2 + t53 ^ 2) / t93 * t59) ^ (-0.1e1 / 0.2e1) * t58 / pkin(4);
t91 = -mrSges(2,1) - m(4) * t86 - (t79 * mrSges(4,1) + t78 * mrSges(4,2)) * t81 - m(5) * pkin(1) - (-mrSges(5,1) * t52 - mrSges(5,2) * t53) * t82 - mrSges(3,1) * t66 + mrSges(3,2) * t64;
t90 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t67 = cos(qJ(1));
t65 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,1) * t64 - mrSges(3,2) * t66 - m(4) * t87 - (-t78 * mrSges(4,1) + t79 * mrSges(4,2)) * t81 - (mrSges(5,1) * t53 - mrSges(5,2) * t52) * t82 + (-m(2) - m(3) - m(4) - m(5)) * pkin(5)) * g(3) + (t91 * t65 - t90 * t67 - mrSges(1,2)) * g(2) + (t90 * t65 + t91 * t67 - mrSges(1,1)) * g(1);
U = t1;
