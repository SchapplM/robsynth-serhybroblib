% Calculate potential energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:18:28
% EndTime: 2020-05-03 01:18:28
% DurationCPUTime: 0.21s
% Computational Cost: add. (167->58), mult. (165->75), div. (0->0), fcn. (101->14), ass. (0->29)
t66 = qJ(2) + qJ(3);
t73 = cos(qJ(2));
t49 = pkin(2) * cos(t66) + t73 * pkin(4) + pkin(1);
t63 = qJ(4) + t66;
t60 = cos(t63);
t84 = pkin(5) * t60 + t49;
t69 = sin(qJ(2));
t82 = pkin(2) * sin(t66) + t69 * pkin(4);
t59 = sin(t63);
t81 = pkin(5) * t59 + t82;
t80 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t79 = t67 * rSges(7,1) + t71 * rSges(7,2);
t78 = rSges(5,1) * t60 - rSges(5,2) * t59 + t49;
t61 = qJ(5) + t63;
t52 = sin(t61);
t53 = cos(t61);
t77 = rSges(6,1) * t53 - rSges(6,2) * t52 + t84;
t68 = sin(qJ(3));
t72 = cos(qJ(3));
t47 = m(3) * rSges(3,1) + (rSges(4,1) * t72 - rSges(4,2) * t68 + pkin(4)) * m(4);
t48 = rSges(3,2) * m(3) + (rSges(4,1) * t68 + rSges(4,2) * t72) * m(4);
t76 = -m(2) * rSges(2,1) - t47 * t73 + t48 * t69 + (-m(3) - m(4)) * pkin(1);
t50 = t71 * rSges(7,1) - rSges(7,2) * t67 + pkin(3);
t75 = -t52 * rSges(7,3) + t50 * t53 + t84;
t74 = cos(qJ(1));
t70 = sin(qJ(1));
t1 = (t76 * g(1) + t80 * g(2)) * t74 + (-t80 * g(1) + t76 * g(2)) * t70 - m(5) * ((-g(2) * rSges(5,3) + g(1) * t78) * t74 + (g(1) * rSges(5,3) + g(2) * t78) * t70) - m(6) * ((-g(2) * rSges(6,3) + g(1) * t77) * t74 + (g(1) * rSges(6,3) + g(2) * t77) * t70) - m(7) * ((g(1) * t75 + g(2) * t79) * t74 + (-g(1) * t79 + g(2) * t75) * t70) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (-t48 * t73 - t47 * t69 - m(1) * rSges(1,3) - m(2) * rSges(2,3) - m(5) * (t59 * rSges(5,1) + t60 * rSges(5,2) + t82) - m(6) * (t52 * rSges(6,1) + t53 * rSges(6,2) + t81) - m(7) * (t53 * rSges(7,3) + t50 * t52 + t81)) * g(3);
U = t1;
