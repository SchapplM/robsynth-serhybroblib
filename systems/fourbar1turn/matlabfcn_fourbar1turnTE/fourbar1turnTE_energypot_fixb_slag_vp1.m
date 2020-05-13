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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnTE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:29
% EndTime: 2020-04-12 19:18:29
% DurationCPUTime: 0.22s
% Computational Cost: add. (424->57), mult. (605->85), div. (36->3), fcn. (195->6), ass. (0->31)
t55 = cos(qJ(2));
t76 = pkin(2) * t55;
t71 = (-0.2e1 * t76 + pkin(1)) * pkin(1);
t69 = pkin(2) ^ 2 + t71;
t70 = pkin(3) ^ 2 - pkin(4) ^ 2;
t47 = t69 + t70;
t51 = pkin(1) * t55 - pkin(2);
t77 = -pkin(3) + pkin(4);
t78 = -pkin(3) - pkin(4);
t46 = sqrt(-((pkin(2) - t78) * (pkin(2) + t78) + t71) * ((pkin(2) - t77) * (pkin(2) + t77) + t71));
t53 = sin(qJ(2));
t74 = t46 * t53;
t42 = -pkin(1) * t74 - t51 * t47;
t45 = pkin(1) * t53 * t47 - t51 * t46;
t49 = 0.1e1 / t69;
t72 = t49 / pkin(3);
t79 = t53 / 0.2e1;
t82 = (t55 * t45 / 0.2e1 + t42 * t79) * t72;
t66 = (t45 * t79 - t55 * t42 / 0.2e1) * t72;
t81 = rSges(4,1) * t66 + rSges(4,2) * t82 + t76;
t48 = t69 - t70;
t50 = pkin(1) - t76;
t80 = pkin(2) * t74 / 0.2e1 - t50 * t48 / 0.2e1;
t75 = t53 * pkin(2);
t73 = t49 / pkin(4);
t68 = rSges(3,1) * t55 - rSges(3,2) * t53;
t44 = t50 * t46 + t48 * t75;
t63 = pkin(1) + (rSges(5,1) * t80 - t44 * rSges(5,2) / 0.2e1) * t73;
t56 = cos(qJ(1));
t54 = sin(qJ(1));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t56 * rSges(2,1) - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) + t56 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t54 * rSges(3,3) + t68 * t56) + g(2) * (-t56 * rSges(3,3) + t68 * t54) + g(3) * (t53 * rSges(3,1) + t55 * rSges(3,2) + pkin(5))) - m(4) * (g(3) * (-rSges(4,1) * t82 + rSges(4,2) * t66 + pkin(5) + t75) + (-g(2) * rSges(4,3) + g(1) * t81) * t56 + (g(1) * rSges(4,3) + g(2) * t81) * t54) - m(5) * (g(1) * (t54 * rSges(5,3) + t63 * t56) + g(2) * (-t56 * rSges(5,3) + t63 * t54) + g(3) * (pkin(5) + (t44 * rSges(5,1) / 0.2e1 + rSges(5,2) * t80) * t73));
U = t1;
