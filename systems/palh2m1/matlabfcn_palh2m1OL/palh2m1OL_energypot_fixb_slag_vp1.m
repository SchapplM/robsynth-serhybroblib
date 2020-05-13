% Calculate potential energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:00:28
% EndTime: 2020-05-03 00:00:28
% DurationCPUTime: 0.23s
% Computational Cost: add. (128->52), mult. (144->65), div. (0->0), fcn. (81->12), ass. (0->26)
t68 = -m(3) - m(4);
t51 = sin(qJ(5));
t55 = cos(qJ(5));
t64 = t55 * rSges(6,1) - t51 * rSges(6,2) + pkin(4);
t63 = pkin(6) + rSges(6,3);
t50 = qJ(2) + qJ(3);
t62 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t49 = qJ(4) + t50;
t45 = sin(t49);
t46 = cos(t49);
t61 = rSges(5,1) * t46 - rSges(5,2) * t45;
t53 = sin(qJ(2));
t60 = -t53 * pkin(2) - pkin(3) * sin(t50) + pkin(5);
t52 = sin(qJ(3));
t56 = cos(qJ(3));
t37 = rSges(3,1) * m(3) + (rSges(4,1) * t56 - rSges(4,2) * t52 + pkin(2)) * m(4);
t40 = rSges(3,2) * m(3) + (rSges(4,1) * t52 + rSges(4,2) * t56) * m(4);
t57 = cos(qJ(2));
t59 = -m(2) * rSges(2,1) + t68 * pkin(1) - t37 * t57 + t40 * t53;
t58 = cos(qJ(1));
t54 = sin(qJ(1));
t42 = t51 * rSges(6,1) + t55 * rSges(6,2);
t41 = pkin(3) * cos(t50) + t57 * pkin(2) + pkin(1);
t39 = t58 * t41;
t38 = t54 * t41;
t1 = (t59 * g(1) - t62 * g(2)) * t58 + (t62 * g(1) + t59 * g(2)) * t54 - m(5) * (g(1) * (-t54 * rSges(5,3) + t61 * t58 + t39) + g(2) * (t58 * rSges(5,3) + t61 * t54 + t38)) - m(6) * (g(1) * (-t54 * t42 + t39) + g(2) * (t58 * t42 + t38) + (t63 * t45 + t64 * t46) * (g(1) * t58 + g(2) * t54)) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (t40 * t57 + t37 * t53 - m(2) * rSges(2,3) - m(1) * rSges(1,3) - m(5) * (-t45 * rSges(5,1) - t46 * rSges(5,2) + t60) - m(6) * (-t64 * t45 + t63 * t46 + t60) + (-m(2) + t68) * pkin(5)) * g(3);
U = t1;
