% Calculate inertial parameters regressor of potential energy for
% picker2Dm2OL
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
% 
% Output:
% U_reg [1x(12*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = picker2Dm2OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_energypot_fixb_reg2_slag_vp: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:19:12
% EndTime: 2020-05-09 23:19:12
% DurationCPUTime: 0.14s
% Computational Cost: add. (127->52), mult. (85->63), div. (0->0), fcn. (68->22), ass. (0->37)
t71 = sin(qJ(1));
t80 = t71 * pkin(1);
t73 = cos(qJ(1));
t79 = t73 * pkin(1);
t69 = qJ(1) + qJ(2);
t66 = qJ(3) + t69;
t65 = qJ(4) + t69;
t61 = sin(t69);
t78 = -pkin(2) * t61 - t80;
t63 = cos(t69);
t77 = -pkin(2) * t63 - t79;
t76 = -pkin(3) * t61 - t80;
t75 = -pkin(3) * t63 - t79;
t74 = g(1) * t73 + g(2) * t71;
t72 = cos(qJ(7));
t70 = sin(qJ(7));
t68 = qJ(1) + qJ(8);
t67 = pkin(8) + qJ(5);
t64 = qJ(6) + t69;
t62 = cos(t68);
t60 = sin(t68);
t59 = cos(t67);
t58 = sin(t67);
t57 = qJ(9) + t66;
t56 = qJ(10) + t65;
t55 = cos(t66);
t54 = cos(t65);
t53 = cos(t64);
t52 = sin(t66);
t51 = sin(t65);
t50 = sin(t64);
t49 = cos(t57);
t48 = sin(t57);
t47 = cos(t56);
t46 = sin(t56);
t45 = t74 * pkin(1);
t1 = [0, 0, 0, 0, 0, 0, t74, -g(1) * t71 + g(2) * t73, -g(3), 0, 0, 0, 0, 0, 0, 0, g(1) * t63 + g(2) * t61, -g(1) * t61 + g(2) * t63, -g(3), t45, 0, 0, 0, 0, 0, 0, -g(1) * t55 - g(2) * t52, g(1) * t52 - g(2) * t55, -g(3), -g(1) * t77 - g(2) * t78, 0, 0, 0, 0, 0, 0, g(1) * t54 + g(2) * t51, -g(1) * t51 + g(2) * t54, -g(3), -g(1) * t75 - g(2) * t76, 0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t58, g(1) * t58 - g(2) * t59, -g(3), (-sin(pkin(8)) * g(2) - cos(pkin(8)) * g(1)) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t50, g(1) * t50 - g(2) * t53, -g(3), t45, 0, 0, 0, 0, 0, 0, -g(1) * t70 + g(2) * t72, -g(1) * t72 - g(2) * t70, -g(3), -g(1) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t62 - g(2) * t60, g(1) * t60 - g(2) * t62, -g(3), t45, 0, 0, 0, 0, 0, 0, g(1) * t49 + g(2) * t48, -g(1) * t48 + g(2) * t49, -g(3), -g(2) * (pkin(6) * t52 + t78) - g(1) * (pkin(6) * t55 + t77), 0, 0, 0, 0, 0, 0, -g(1) * t47 - g(2) * t46, g(1) * t46 - g(2) * t47, -g(3), -g(2) * (-pkin(4) * t51 + t76) - g(1) * (-pkin(4) * t54 + t75);];
U_reg = t1;
