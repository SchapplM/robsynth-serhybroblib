% Calculate inertial parameters regressor of potential energy for
% fourbar1DE1
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
% 
% Output:
% U_reg [1x(1*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1DE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_energypot_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:17
% EndTime: 2020-04-24 19:57:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (104->29), mult. (156->31), div. (8->3), fcn. (46->4), ass. (0->22)
t71 = cos(qJ(1));
t84 = pkin(2) * t71;
t83 = (-0.2e1 * t84 + pkin(1)) * pkin(1);
t81 = pkin(2) ^ 2 + t83;
t88 = 0.1e1 / t81 / 0.2e1;
t87 = -pkin(3) - pkin(4);
t86 = -pkin(3) + pkin(4);
t70 = sin(qJ(1));
t85 = pkin(2) * t70;
t82 = pkin(3) ^ 2 - pkin(4) ^ 2;
t80 = 0.1e1 / pkin(4) * t88;
t79 = 0.1e1 / pkin(3) * t88;
t78 = -pkin(1) + t84;
t68 = g(1) * t71 + g(2) * t70;
t66 = t81 - t82;
t65 = t81 + t82;
t64 = t78 * g(1) + g(2) * t85;
t63 = g(1) * t85 - t78 * g(2);
t62 = sqrt(-((pkin(2) - t87) * (pkin(2) + t87) + t83) * ((pkin(2) - t86) * (pkin(2) + t86) + t83));
t61 = t64 * t62;
t60 = t63 * t62;
t1 = [0, 0, 0, 0, 0, 0, -t68, g(1) * t70 - g(2) * t71, -g(3), 0, 0, 0, 0, 0, 0, 0, (t65 * t64 - t60) * t79, (-t65 * t63 - t61) * t79, -g(3), -pkin(2) * t68, 0, 0, 0, 0, 0, 0, (-t66 * t64 - t60) * t80, (t66 * t63 - t61) * t80, -g(3), -g(1) * pkin(1);];
U_reg = t1;
