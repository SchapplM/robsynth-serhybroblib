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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:23
% EndTime: 2020-04-12 19:33:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (682->58), mult. (955->82), div. (60->5), fcn. (299->11), ass. (0->28)
t77 = pkin(4) ^ 2;
t76 = -pkin(3) - pkin(4);
t75 = -pkin(3) + pkin(4);
t57 = cos(qJ(2));
t74 = pkin(2) * t57;
t55 = sin(qJ(2));
t73 = t55 * pkin(2);
t70 = (-0.2e1 * t74 + pkin(1)) * pkin(1);
t46 = sqrt(-((pkin(2) - t76) * (pkin(2) + t76) + t70) * ((pkin(2) - t75) * (pkin(2) + t75) + t70));
t72 = t46 * t55;
t51 = pkin(2) ^ 2 + t70;
t49 = 0.1e1 / t51;
t71 = t49 / pkin(3);
t69 = pkin(3) ^ 2 - t77;
t68 = rSges(3,1) * t57 - rSges(3,2) * t55;
t47 = t51 + t69;
t53 = pkin(1) * t57 - pkin(2);
t42 = qJ(2) + atan2((pkin(1) * t55 * t47 - t53 * t46) * t71, (-pkin(1) * t72 - t53 * t47) * t71);
t40 = sin(t42);
t41 = cos(t42);
t66 = -rSges(4,1) * t41 + rSges(4,2) * t40 + t74;
t58 = cos(qJ(1));
t56 = sin(qJ(1));
t52 = pkin(1) - t74;
t48 = t51 - t69;
t45 = t52 * t46 + t48 * t73;
t44 = -pkin(2) * t72 + t52 * t48;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t58 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) + t58 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t56 * rSges(3,3) + t68 * t58) + g(2) * (-t58 * rSges(3,3) + t68 * t56) + g(3) * (t55 * rSges(3,1) + t57 * rSges(3,2) + pkin(5))) - m(4) * (g(3) * (-t40 * rSges(4,1) - t41 * rSges(4,2) + pkin(5) + t73) + (-g(2) * rSges(4,3) + g(1) * t66) * t58 + (g(1) * rSges(4,3) + g(2) * t66) * t56) - m(5) * (g(1) * (t56 * rSges(5,3) + t58 * pkin(1)) + g(2) * (-t58 * rSges(5,3) + t56 * pkin(1)) + g(3) * pkin(5) + (g(3) * (rSges(5,1) * t45 - rSges(5,2) * t44) + (g(1) * t58 + g(2) * t56) * (-rSges(5,1) * t44 - rSges(5,2) * t45)) / pkin(4) * t49 * ((t44 ^ 2 + t45 ^ 2) / t77 / t51 ^ 2) ^ (-0.1e1 / 0.2e1));
U = t1;
