% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnTE_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energypot_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:30
% EndTime: 2020-06-27 16:22:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (398->32), mult. (572->53), div. (36->3), fcn. (182->6), ass. (0->31)
t86 = cos(qJ(2));
t105 = pkin(2) * t86;
t101 = (-0.2e1 * t105 + pkin(1)) * pkin(1);
t99 = pkin(2) ^ 2 + t101;
t80 = 0.1e1 / t99;
t102 = t80 / pkin(3);
t84 = sin(qJ(2));
t110 = t84 / 0.2e1;
t108 = -pkin(3) + pkin(4);
t109 = -pkin(3) - pkin(4);
t77 = sqrt(-((pkin(2) - t109) * (pkin(2) + t109) + t101) * ((pkin(2) - t108) * (pkin(2) + t108) + t101));
t104 = t77 * t84;
t100 = pkin(3) ^ 2 - pkin(4) ^ 2;
t78 = t99 + t100;
t82 = pkin(1) * t86 - pkin(2);
t73 = -pkin(1) * t104 - t82 * t78;
t76 = pkin(1) * t84 * t78 - t82 * t77;
t111 = (t86 * t76 / 0.2e1 + t73 * t110) * t102;
t85 = sin(qJ(1));
t106 = g(2) * t85;
t87 = cos(qJ(1));
t107 = g(1) * t87;
t98 = t106 + t107;
t103 = t80 / pkin(4);
t97 = t107 / 0.2e1 + t106 / 0.2e1;
t94 = (t76 * t110 - t86 * t73 / 0.2e1) * t102;
t81 = pkin(1) - t105;
t79 = t99 - t100;
t75 = pkin(2) * t84 * t79 + t81 * t77;
t74 = -pkin(2) * t104 + t81 * t79;
t1 = [0, -t98, g(1) * t85 - g(2) * t87, 0, 0, 0, 0, 0, -g(3) * t84 - t98 * t86, -g(3) * t86 + t98 * t84, 0, 0, 0, 0, 0, 0, g(3) * t111 - t98 * t94, -g(3) * t94 - t98 * t111, 0, 0, 0, 0, 0, (-g(3) * t75 / 0.2e1 + t97 * t74) * t103, (g(3) * t74 / 0.2e1 + t97 * t75) * t103;];
U_reg = t1;
