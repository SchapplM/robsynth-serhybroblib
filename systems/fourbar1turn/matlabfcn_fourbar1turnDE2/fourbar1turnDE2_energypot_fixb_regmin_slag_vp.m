% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnDE2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energypot_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:48:27
% EndTime: 2020-06-27 16:48:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (656->31), mult. (922->47), div. (60->5), fcn. (286->11), ass. (0->27)
t115 = pkin(4) ^ 2;
t114 = -pkin(3) - pkin(4);
t113 = -pkin(3) + pkin(4);
t97 = cos(qJ(2));
t112 = pkin(2) * t97;
t109 = (-0.2e1 * t112 + pkin(1)) * pkin(1);
t86 = sqrt(-((pkin(2) - t114) * (pkin(2) + t114) + t109) * ((pkin(2) - t113) * (pkin(2) + t113) + t109));
t95 = sin(qJ(2));
t111 = t86 * t95;
t91 = pkin(2) ^ 2 + t109;
t89 = 0.1e1 / t91;
t110 = 0.1e1 / pkin(3) * t89;
t108 = -pkin(3) ^ 2 + t115;
t88 = t91 + t108;
t92 = pkin(1) - t112;
t84 = -pkin(2) * t111 + t92 * t88;
t85 = pkin(2) * t95 * t88 + t92 * t86;
t107 = 0.1e1 / pkin(4) * ((t84 ^ 2 + t85 ^ 2) / t91 ^ 2 / t115) ^ (-0.1e1 / 0.2e1) * t89;
t96 = sin(qJ(1));
t98 = cos(qJ(1));
t106 = g(1) * t98 + g(2) * t96;
t93 = pkin(1) * t97 - pkin(2);
t87 = t91 - t108;
t82 = qJ(2) + atan2((pkin(1) * t95 * t87 - t93 * t86) * t110, (-pkin(1) * t111 - t93 * t87) * t110);
t81 = cos(t82);
t80 = sin(t82);
t1 = [0, -t106, g(1) * t96 - g(2) * t98, 0, 0, 0, 0, 0, -g(3) * t95 - t106 * t97, -g(3) * t97 + t106 * t95, 0, 0, 0, 0, 0, 0, g(3) * t80 + t106 * t81, g(3) * t81 - t106 * t80, 0, 0, 0, 0, 0, (-g(3) * t85 + t106 * t84) * t107, (g(3) * t84 + t106 * t85) * t107;];
U_reg = t1;
