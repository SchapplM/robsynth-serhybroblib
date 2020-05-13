% Calculate inertial parameters regressor of potential energy for
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
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnDE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:27:13
% EndTime: 2020-04-12 19:27:13
% DurationCPUTime: 0.21s
% Computational Cost: add. (1183->41), mult. (1727->60), div. (108->6), fcn. (519->10), ass. (0->33)
t131 = pkin(4) ^ 2;
t130 = pkin(3) ^ 2;
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t118 = g(1) * t106 + g(2) * t104;
t129 = g(3) * pkin(5);
t128 = -pkin(3) - pkin(4);
t127 = -pkin(3) + pkin(4);
t105 = cos(qJ(2));
t124 = pkin(2) * t105;
t103 = sin(qJ(2));
t122 = (-0.2e1 * t124 + pkin(1)) * pkin(1);
t92 = sqrt(-((pkin(2) - t128) * (pkin(2) + t128) + t122) * ((pkin(2) - t127) * (pkin(2) + t127) + t122));
t123 = t103 * t92;
t121 = t130 - t131;
t98 = pkin(2) ^ 2 + t122;
t94 = t98 - t121;
t99 = pkin(1) - t124;
t89 = -pkin(2) * t123 + t99 * t94;
t90 = pkin(2) * t103 * t94 + t99 * t92;
t95 = 0.1e1 / t98;
t96 = 0.1e1 / t98 ^ 2;
t120 = 0.1e1 / pkin(4) * ((t89 ^ 2 + t90 ^ 2) * t96 / t131) ^ (-0.1e1 / 0.2e1) * t95;
t100 = pkin(1) * t105 - pkin(2);
t93 = t98 + t121;
t88 = -pkin(1) * t123 - t100 * t93;
t91 = pkin(1) * t103 * t93 - t100 * t92;
t119 = 0.1e1 / pkin(3) * ((t88 ^ 2 + t91 ^ 2) * t96 / t130) ^ (-0.1e1 / 0.2e1) * t95;
t117 = t103 * t91 - t105 * t88;
t116 = t103 * t88 + t105 * t91;
t115 = -g(3) * t103 - t105 * t118;
t97 = g(1) * t104 - g(2) * t106;
t1 = [0, 0, 0, 0, 0, 0, -t118, t97, -g(3), -t129, 0, 0, 0, 0, 0, 0, t115, -g(3) * t105 + t103 * t118, -t97, -t129, 0, 0, 0, 0, 0, 0, (g(3) * t116 - t118 * t117) * t119, (-g(3) * t117 - t118 * t116) * t119, -t97, pkin(2) * t115 - t129, 0, 0, 0, 0, 0, 0, (-g(3) * t90 + t118 * t89) * t120, (g(3) * t89 + t118 * t90) * t120, -t97, -pkin(1) * t118 - t129;];
U_reg = t1;
