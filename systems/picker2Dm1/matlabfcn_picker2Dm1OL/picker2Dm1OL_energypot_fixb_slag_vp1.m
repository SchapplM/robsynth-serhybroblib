% Calculate potential energy for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_energypot_fixb_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1OL_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:49
% EndTime: 2020-05-11 05:44:50
% DurationCPUTime: 0.26s
% Computational Cost: add. (157->85), mult. (112->96), div. (0->0), fcn. (68->22), ass. (0->35)
t60 = sin(qJ(1));
t68 = t60 * pkin(1);
t62 = cos(qJ(1));
t67 = t62 * pkin(1);
t58 = qJ(1) + qJ(2);
t55 = qJ(3) + t58;
t54 = qJ(4) + t58;
t50 = sin(t58);
t66 = -pkin(2) * t50 - t68;
t52 = cos(t58);
t65 = -pkin(2) * t52 - t67;
t64 = -pkin(3) * t50 - t68;
t63 = -pkin(3) * t52 - t67;
t61 = cos(qJ(7));
t59 = sin(qJ(7));
t57 = qJ(1) + qJ(8);
t56 = pkin(8) + qJ(5);
t53 = qJ(6) + t58;
t51 = cos(t57);
t49 = sin(t57);
t48 = cos(t56);
t47 = sin(t56);
t46 = qJ(9) + t55;
t45 = qJ(10) + t54;
t44 = cos(t55);
t43 = cos(t54);
t42 = cos(t53);
t41 = sin(t55);
t40 = sin(t54);
t39 = sin(t53);
t38 = cos(t46);
t37 = sin(t46);
t36 = cos(t45);
t35 = sin(t45);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (-t62 * rSges(2,1) + t60 * rSges(2,2)) + g(2) * (-t60 * rSges(2,1) - t62 * rSges(2,2)) + g(3) * rSges(2,3)) - m(3) * (g(1) * (-t52 * rSges(3,1) + t50 * rSges(3,2) - t67) + g(2) * (-t50 * rSges(3,1) - t52 * rSges(3,2) - t68) + g(3) * rSges(3,3)) - m(4) * (g(1) * (t44 * rSges(4,1) - t41 * rSges(4,2) + t65) + g(2) * (t41 * rSges(4,1) + t44 * rSges(4,2) + t66) + g(3) * rSges(4,3)) - m(5) * (g(1) * (-t43 * rSges(5,1) + t40 * rSges(5,2) + t63) + g(2) * (-t40 * rSges(5,1) - t43 * rSges(5,2) + t64) + g(3) * rSges(5,3)) - m(6) * (g(1) * (cos(pkin(8)) * pkin(5) + t48 * rSges(6,1) - t47 * rSges(6,2)) + g(2) * (sin(pkin(8)) * pkin(5) + t47 * rSges(6,1) + t48 * rSges(6,2)) + g(3) * rSges(6,3)) - m(7) * (g(1) * (t42 * rSges(7,1) - t39 * rSges(7,2) - t67) + g(2) * (t39 * rSges(7,1) + t42 * rSges(7,2) - t68) + g(3) * rSges(7,3)) - m(8) * (g(1) * (t59 * rSges(8,1) + t61 * rSges(8,2) + pkin(7)) + g(2) * (-t61 * rSges(8,1) + t59 * rSges(8,2)) + g(3) * rSges(8,3)) - m(9) * (g(1) * (t51 * rSges(9,1) - t49 * rSges(9,2) - t67) + g(2) * (t49 * rSges(9,1) + t51 * rSges(9,2) - t68) + g(3) * rSges(9,3)) - m(10) * (g(1) * (pkin(6) * t44 - t38 * rSges(10,1) + t37 * rSges(10,2) + t65) + g(2) * (pkin(6) * t41 - t37 * rSges(10,1) - t38 * rSges(10,2) + t66) + g(3) * rSges(10,3)) - m(11) * (g(1) * (-pkin(4) * t43 + t36 * rSges(11,1) - t35 * rSges(11,2) + t63) + g(2) * (-pkin(4) * t40 + t35 * rSges(11,1) + t36 * rSges(11,2) + t64) + g(3) * rSges(11,3));
U = t1;
