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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:25
% EndTime: 2020-04-12 19:25:26
% DurationCPUTime: 0.35s
% Computational Cost: add. (1198->60), mult. (1739->93), div. (108->6), fcn. (519->10), ass. (0->30)
t89 = pkin(4) ^ 2;
t88 = pkin(3) ^ 2;
t62 = sin(qJ(1));
t64 = cos(qJ(1));
t87 = g(1) * t64 + g(2) * t62;
t86 = -pkin(3) - pkin(4);
t85 = -pkin(3) + pkin(4);
t63 = cos(qJ(2));
t82 = pkin(2) * t63;
t61 = sin(qJ(2));
t81 = t61 * pkin(2);
t79 = (-0.2e1 * t82 + pkin(1)) * pkin(1);
t52 = sqrt(-((pkin(2) - t86) * (pkin(2) + t86) + t79) * ((pkin(2) - t85) * (pkin(2) + t85) + t79));
t80 = t52 * t61;
t78 = t88 - t89;
t57 = pkin(2) ^ 2 + t79;
t77 = rSges(3,1) * t63 - rSges(3,2) * t61;
t53 = t57 + t78;
t59 = pkin(1) * t63 - pkin(2);
t48 = -pkin(1) * t80 - t53 * t59;
t51 = pkin(1) * t53 * t61 - t52 * t59;
t75 = -t48 * t63 + t51 * t61;
t74 = t48 * t61 + t51 * t63;
t58 = pkin(1) - t82;
t56 = 0.1e1 / t57 ^ 2;
t55 = 0.1e1 / t57;
t54 = t57 - t78;
t50 = t52 * t58 + t54 * t81;
t49 = -pkin(2) * t80 + t54 * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t64 - rSges(2,2) * t62) + g(2) * (rSges(2,1) * t62 + rSges(2,2) * t64) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t62 * rSges(3,3) + t64 * t77) + g(2) * (-t64 * rSges(3,3) + t62 * t77) + g(3) * (rSges(3,1) * t61 + rSges(3,2) * t63 + pkin(5))) - m(4) * (g(1) * (t62 * rSges(4,3) + t64 * t82) + g(2) * (-t64 * rSges(4,3) + t62 * t82) + g(3) * (pkin(5) + t81) + (g(3) * (-rSges(4,1) * t74 + rSges(4,2) * t75) + t87 * (rSges(4,1) * t75 + rSges(4,2) * t74)) / pkin(3) * t55 * ((t48 ^ 2 + t51 ^ 2) / t88 * t56) ^ (-0.1e1 / 0.2e1)) - m(5) * (g(1) * (t62 * rSges(5,3) + t64 * pkin(1)) + g(2) * (-t64 * rSges(5,3) + t62 * pkin(1)) + g(3) * pkin(5) + (g(3) * (rSges(5,1) * t50 - rSges(5,2) * t49) + t87 * (-rSges(5,1) * t49 - rSges(5,2) * t50)) / pkin(4) * t55 * ((t49 ^ 2 + t50 ^ 2) / t89 * t56) ^ (-0.1e1 / 0.2e1));
U = t1;
