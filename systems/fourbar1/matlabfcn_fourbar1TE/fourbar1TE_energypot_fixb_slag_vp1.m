% Calculate potential energy for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1TE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_energypot_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1TE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1TE_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:48:56
% EndTime: 2020-04-24 19:48:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (112->39), mult. (170->51), div. (8->3), fcn. (44->4), ass. (0->21)
t39 = cos(qJ(1));
t53 = pkin(2) * t39;
t51 = (-0.2e1 * t53 + pkin(1)) * pkin(1);
t49 = pkin(2) ^ 2 + t51;
t56 = 0.1e1 / t49 / 0.2e1;
t55 = -pkin(3) - pkin(4);
t54 = -pkin(3) + pkin(4);
t38 = sin(qJ(1));
t52 = t38 * pkin(2);
t50 = pkin(3) ^ 2 - pkin(4) ^ 2;
t48 = 0.1e1 / pkin(4) * t56;
t47 = 0.1e1 / pkin(3) * t56;
t46 = -pkin(1) + t53;
t35 = t49 - t50;
t34 = t49 + t50;
t33 = -t46 * rSges(3,1) + rSges(3,2) * t52;
t32 = rSges(3,1) * t52 + t46 * rSges(3,2);
t31 = rSges(4,1) * t52 + t46 * rSges(4,2);
t30 = -t46 * rSges(4,1) + rSges(4,2) * t52;
t29 = sqrt(-((pkin(2) - t55) * (pkin(2) + t55) + t51) * ((pkin(2) - t54) * (pkin(2) + t54) + t51));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - ((rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t39 + (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t38 + g(3) * rSges(2,3)) * m(2) - m(3) * (g(1) * (t53 + (t32 * t29 + t34 * t33) * t47) + g(2) * (t52 + (t33 * t29 - t34 * t32) * t47) + g(3) * rSges(3,3)) - m(4) * (g(1) * (pkin(1) + (t31 * t29 - t35 * t30) * t48) + g(2) * (t30 * t29 + t35 * t31) * t48 + g(3) * rSges(4,3));
U = t1;
