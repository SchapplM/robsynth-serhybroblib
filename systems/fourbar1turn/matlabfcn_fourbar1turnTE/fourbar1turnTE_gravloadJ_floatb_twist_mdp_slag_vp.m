% Calculate Gravitation load on the joints for
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnTE_convert_par2_MPV_fixb.m
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnTE_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:37
% EndTime: 2020-06-27 16:22:44
% DurationCPUTime: 0.70s
% Computational Cost: add. (1836->80), mult. (2646->147), div. (120->5), fcn. (791->6), ass. (0->68)
t108 = cos(qJ(2));
t103 = pkin(1) * t108 - pkin(2);
t106 = sin(qJ(2));
t115 = pkin(1) ^ 2;
t145 = pkin(2) * t108;
t133 = -0.2e1 * pkin(1) * t145 + t115;
t153 = -pkin(3) - pkin(4);
t95 = (pkin(2) - t153) * (pkin(2) + t153) + t133;
t152 = -pkin(3) + pkin(4);
t96 = (pkin(2) - t152) * (pkin(2) + t152) + t133;
t116 = sqrt(-t95 * t96);
t135 = t106 * t116;
t114 = pkin(2) ^ 2;
t101 = t114 + t133;
t132 = pkin(3) ^ 2 - pkin(4) ^ 2;
t97 = t101 + t132;
t86 = -pkin(1) * t135 - t103 * t97;
t159 = t86 / 0.2e1;
t94 = pkin(1) * t106 * t97;
t89 = -t103 * t116 + t94;
t158 = -t89 / 0.2e1;
t157 = t106 * t158 + t108 * t159;
t107 = sin(qJ(1));
t148 = g(2) * t107;
t109 = cos(qJ(1));
t149 = g(1) * t109;
t123 = t148 + t149;
t134 = t108 * t116;
t146 = pkin(2) * t106;
t130 = pkin(1) * t146;
t144 = 0.1e1 / t116 * (-t95 - t96) * t130;
t78 = t94 + (-t134 + (-0.2e1 * t103 * pkin(2) - t144) * t106) * pkin(1);
t126 = t158 + t78 / 0.2e1;
t105 = t106 ^ 2;
t155 = 0.2e1 * t105;
t80 = -t103 * t144 + t115 * pkin(2) * t155 + (t108 * t97 + t135) * pkin(1);
t154 = t80 / 0.2e1;
t127 = t159 + t154;
t156 = t126 * t106 + t127 * t108;
t151 = t106 / 0.2e1;
t150 = g(1) * t107;
t147 = g(2) * t109;
t113 = 0.1e1 / pkin(3);
t99 = 0.1e1 / t101;
t137 = t113 * t99;
t129 = t107 * t137;
t143 = t157 * t129;
t128 = t109 * t137;
t138 = t108 * t89;
t141 = t106 * t86;
t142 = (t138 + t141) * t128 / 0.2e1;
t136 = t106 * t108;
t100 = 0.1e1 / t101 ^ 2;
t131 = pkin(1) * pkin(2) * t100;
t125 = t100 * t130;
t121 = t149 / 0.2e1 + t148 / 0.2e1;
t120 = t105 * t86 + t89 * t136;
t119 = -t105 * t89 + t86 * t136;
t117 = t127 * t106 - t126 * t108;
t111 = 0.1e1 / pkin(4);
t102 = pkin(1) - t145;
t98 = t101 - t132;
t93 = t98 * t146;
t88 = t102 * t116 + t93;
t87 = -pkin(2) * t135 + t102 * t98;
t81 = t102 * t144 + t114 * pkin(1) * t155 + (t108 * t98 + t135) * pkin(2);
t79 = t93 + (-t134 + (0.2e1 * t102 * pkin(1) - t144) * t106) * pkin(2);
t1 = [t123 * MDP(3) + (g(2) * t157 * t128 - g(1) * t143) * MDP(17) + (-g(2) * t142 - g(1) * (-t138 / 0.2e1 - t141 / 0.2e1) * t129) * MDP(18) + (t87 * MDP(24) + t88 * MDP(25)) * t111 * t99 * (-t150 / 0.2e1 + t147 / 0.2e1) + (-t106 * MDP(10) + t108 * MDP(9) + MDP(2)) * (-t147 + t150); (-g(3) * t108 + t123 * t106) * MDP(9) + (g(3) * t106 + t123 * t108) * MDP(10) + (-g(1) * t142 + ((-g(3) * t120 - t123 * t119) * t131 + (g(3) * t156 - (-t108 * t78 / 0.2e1 + t80 * t151) * t149 - t117 * t148) * t99) * t113) * MDP(17) + (-g(2) * t143 + ((-g(3) * t119 + t123 * t120) * t131 + (-g(3) * t117 - (t108 * t154 + t78 * t151) * t148 - t156 * t149) * t99) * t113) * MDP(18) + (((-g(3) * t81 / 0.2e1 + t121 * t79) * t99 + (g(3) * t88 - t123 * t87) * t125) * MDP(24) + ((g(3) * t79 / 0.2e1 + t121 * t81) * t99 + (-g(3) * t87 - t123 * t88) * t125) * MDP(25)) * t111;];
taug = t1;
