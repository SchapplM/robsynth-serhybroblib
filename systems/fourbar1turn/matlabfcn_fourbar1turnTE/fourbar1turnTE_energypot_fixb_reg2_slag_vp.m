% Calculate inertial parameters regressor of potential energy for
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
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnTE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:19:50
% EndTime: 2020-04-12 19:19:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (409->40), mult. (593->56), div. (36->3), fcn. (195->6), ass. (0->34)
t98 = cos(qJ(2));
t118 = pkin(2) * t98;
t114 = (-0.2e1 * t118 + pkin(1)) * pkin(1);
t112 = pkin(2) ^ 2 + t114;
t90 = 0.1e1 / t112;
t115 = 0.1e1 / pkin(3) * t90;
t96 = sin(qJ(2));
t123 = t96 / 0.2e1;
t121 = -pkin(3) + pkin(4);
t122 = -pkin(3) - pkin(4);
t87 = sqrt(-((pkin(2) - t122) * (pkin(2) + t122) + t114) * ((pkin(2) - t121) * (pkin(2) + t121) + t114));
t117 = t87 * t96;
t113 = pkin(3) ^ 2 - pkin(4) ^ 2;
t88 = t112 + t113;
t93 = pkin(1) * t98 - pkin(2);
t83 = -pkin(1) * t117 - t93 * t88;
t86 = pkin(1) * t96 * t88 - t93 * t87;
t125 = (t98 * t86 / 0.2e1 + t83 * t123) * t115;
t97 = sin(qJ(1));
t119 = g(2) * t97;
t99 = cos(qJ(1));
t120 = g(1) * t99;
t111 = t119 + t120;
t124 = g(3) * pkin(5);
t116 = 0.1e1 / pkin(4) * t90;
t110 = t120 / 0.2e1 + t119 / 0.2e1;
t107 = (t86 * t123 - t98 * t83 / 0.2e1) * t115;
t106 = -g(3) * t96 - t111 * t98;
t92 = pkin(1) - t118;
t91 = g(1) * t97 - g(2) * t99;
t89 = t112 - t113;
t85 = pkin(2) * t96 * t89 + t92 * t87;
t84 = -pkin(2) * t117 + t92 * t89;
t1 = [0, 0, 0, 0, 0, 0, -t111, t91, -g(3), -t124, 0, 0, 0, 0, 0, 0, t106, -g(3) * t98 + t111 * t96, -t91, -t124, 0, 0, 0, 0, 0, 0, g(3) * t125 - t111 * t107, -g(3) * t107 - t111 * t125, -t91, t106 * pkin(2) - t124, 0, 0, 0, 0, 0, 0, (-g(3) * t85 / 0.2e1 + t110 * t84) * t116, (g(3) * t84 / 0.2e1 + t110 * t85) * t116, -t91, -t111 * pkin(1) - t124;];
U_reg = t1;
