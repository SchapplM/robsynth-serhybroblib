% Calculate potential energy for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:23
% EndTime: 2020-04-12 19:33:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (682->48), mult. (958->53), div. (60->5), fcn. (299->11), ass. (0->29)
t82 = pkin(4) ^ 2;
t81 = -pkin(3) - pkin(4);
t80 = -pkin(3) + pkin(4);
t62 = cos(qJ(2));
t79 = pkin(2) * t62;
t76 = (-0.2e1 * t79 + pkin(1)) * pkin(1);
t51 = sqrt(-((pkin(2) - t81) * (pkin(2) + t81) + t76) * ((pkin(2) - t80) * (pkin(2) + t80) + t76));
t60 = sin(qJ(2));
t78 = t51 * t60;
t56 = pkin(2) ^ 2 + t76;
t54 = 0.1e1 / t56;
t77 = t54 / pkin(3);
t75 = pkin(3) ^ 2 - t82;
t53 = t56 - t75;
t57 = pkin(1) - t79;
t49 = -pkin(2) * t78 + t57 * t53;
t50 = pkin(2) * t60 * t53 + t57 * t51;
t74 = ((t49 ^ 2 + t50 ^ 2) / t82 / t56 ^ 2) ^ (-0.1e1 / 0.2e1) * t54 / pkin(4);
t73 = -m(4) * pkin(2) - mrSges(3,1);
t72 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t52 = t56 + t75;
t58 = pkin(1) * t62 - pkin(2);
t47 = qJ(2) + atan2((pkin(1) * t60 * t52 - t58 * t51) * t77, (-pkin(1) * t78 - t58 * t52) * t77);
t45 = sin(t47);
t46 = cos(t47);
t71 = -m(5) * pkin(1) + mrSges(4,1) * t46 + mrSges(3,2) * t60 - mrSges(4,2) * t45 - mrSges(2,1) + t73 * t62 + (mrSges(5,1) * t49 + mrSges(5,2) * t50) * t74;
t63 = cos(qJ(1));
t61 = sin(qJ(1));
t1 = (t45 * mrSges(4,1) - t62 * mrSges(3,2) + t46 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t73 * t60 + (-mrSges(5,1) * t50 + mrSges(5,2) * t49) * t74 + (-m(2) - m(3) - m(4) - m(5)) * pkin(5)) * g(3) + (t71 * t61 - t72 * t63 - mrSges(1,2)) * g(2) + (t72 * t61 + t71 * t63 - mrSges(1,1)) * g(1);
U = t1;
