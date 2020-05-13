% Calculate inertial parameters regressor of potential energy for
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
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnDE2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:34:55
% EndTime: 2020-04-12 19:34:55
% DurationCPUTime: 0.24s
% Computational Cost: add. (667->39), mult. (943->50), div. (60->5), fcn. (299->11), ass. (0->30)
t125 = pkin(4) ^ 2;
t124 = g(3) * pkin(5);
t123 = -pkin(3) - pkin(4);
t122 = -pkin(3) + pkin(4);
t105 = cos(qJ(2));
t121 = pkin(2) * t105;
t103 = sin(qJ(2));
t118 = (-0.2e1 * t121 + pkin(1)) * pkin(1);
t92 = sqrt(-((pkin(2) - t123) * (pkin(2) + t123) + t118) * ((pkin(2) - t122) * (pkin(2) + t122) + t118));
t120 = t103 * t92;
t98 = pkin(2) ^ 2 + t118;
t95 = 0.1e1 / t98;
t119 = 0.1e1 / pkin(3) * t95;
t117 = pkin(3) ^ 2 - t125;
t94 = t98 - t117;
t99 = pkin(1) - t121;
t90 = -pkin(2) * t120 + t94 * t99;
t91 = pkin(2) * t103 * t94 + t92 * t99;
t116 = 0.1e1 / pkin(4) * ((t90 ^ 2 + t91 ^ 2) / t98 ^ 2 / t125) ^ (-0.1e1 / 0.2e1) * t95;
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t115 = g(1) * t106 + g(2) * t104;
t114 = -g(3) * t103 - t105 * t115;
t100 = pkin(1) * t105 - pkin(2);
t97 = g(1) * t104 - g(2) * t106;
t93 = t98 + t117;
t88 = qJ(2) + atan2((pkin(1) * t103 * t93 - t100 * t92) * t119, (-pkin(1) * t120 - t100 * t93) * t119);
t87 = cos(t88);
t86 = sin(t88);
t1 = [0, 0, 0, 0, 0, 0, -t115, t97, -g(3), -t124, 0, 0, 0, 0, 0, 0, t114, -g(3) * t105 + t103 * t115, -t97, -t124, 0, 0, 0, 0, 0, 0, g(3) * t86 + t115 * t87, g(3) * t87 - t115 * t86, -t97, pkin(2) * t114 - t124, 0, 0, 0, 0, 0, 0, (-g(3) * t91 + t115 * t90) * t116, (g(3) * t90 + t115 * t91) * t116, -t97, -pkin(1) * t115 - t124;];
U_reg = t1;
