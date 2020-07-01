% Calculate Gravitation load on the joints for
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnDE2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnDE2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:48:34
% EndTime: 2020-06-27 16:48:40
% DurationCPUTime: 0.72s
% Computational Cost: add. (4834->71), mult. (6804->132), div. (304->10), fcn. (1917->11), ass. (0->59)
t162 = pkin(4) ^ 2;
t161 = 2 * pkin(1);
t126 = pkin(2) ^ 2;
t127 = pkin(1) ^ 2;
t119 = cos(qJ(2));
t148 = pkin(2) * t119;
t141 = -0.2e1 * pkin(1) * t148 + t127;
t112 = t126 + t141;
t109 = 0.1e1 / t112;
t140 = pkin(3) ^ 2 - t162;
t108 = t112 - t140;
t113 = pkin(1) - t148;
t117 = sin(qJ(2));
t155 = -pkin(3) - pkin(4);
t105 = (pkin(2) - t155) * (pkin(2) + t155) + t141;
t154 = -pkin(3) + pkin(4);
t106 = (pkin(2) - t154) * (pkin(2) + t154) + t141;
t129 = sqrt(-t105 * t106);
t143 = t117 * t129;
t97 = -pkin(2) * t143 + t113 * t108;
t149 = pkin(2) * t117;
t103 = t108 * t149;
t98 = t113 * t129 + t103;
t146 = t97 ^ 2 + t98 ^ 2;
t110 = 0.1e1 / t112 ^ 2;
t159 = t110 / t162;
t91 = t146 * t159;
t89 = t91 ^ (-0.1e1 / 0.2e1);
t160 = t109 * t89;
t118 = sin(qJ(1));
t150 = g(2) * t118;
t120 = cos(qJ(1));
t151 = g(1) * t120;
t158 = -t151 / 0.2e1 - t150 / 0.2e1;
t134 = t150 + t151;
t157 = -0.2e1 * pkin(2);
t156 = 0.2e1 * t117 ^ 2;
t153 = g(3) * t97;
t152 = g(3) * t98;
t147 = t117 * pkin(1);
t139 = pkin(2) * t147;
t145 = 0.1e1 / t129 * (-t105 - t106) * t139;
t125 = 0.1e1 / pkin(3);
t144 = t109 * t125;
t142 = t129 * t119;
t138 = t110 * t149;
t122 = 0.1e1 / pkin(4);
t114 = pkin(1) * t119 - pkin(2);
t107 = t112 + t140;
t104 = t107 * t147;
t99 = -t114 * t129 + t104;
t96 = -pkin(1) * t143 - t114 * t107;
t93 = 0.1e1 / t96 ^ 2;
t88 = qJ(2) + atan2(t99 * t144, t96 * t144);
t87 = cos(t88);
t86 = sin(t88);
t85 = t113 * t145 + t126 * pkin(1) * t156 + (t119 * t108 + t143) * pkin(2);
t84 = t103 + (-t142 + (t113 * t161 - t145) * t117) * pkin(2);
t1 = [t134 * MDP(3) + (MDP(2) - t87 * MDP(17) + t86 * MDP(18) + t119 * MDP(9) - t117 * MDP(10) + (-t97 * MDP(24) - t98 * MDP(25)) * t122 * t160) * (g(1) * t118 - g(2) * t120); (-g(3) * t119 + t134 * t117) * MDP(9) + (g(3) * t117 + t134 * t119) * MDP(10) + ((g(3) * t87 - t134 * t86) * MDP(17) + (-g(3) * t86 - t134 * t87) * MDP(18)) * (0.1e1 + (((t127 * pkin(2) * t156 - t114 * t145) * t109 + ((t119 * t107 + t143) * t109 - 0.2e1 * t99 * t138) * pkin(1)) / t96 - (t104 * t109 + (-t109 * t142 + ((t114 * t157 - t145) * t109 + t96 * t110 * t157) * t117) * pkin(1)) * t99 * t93) * pkin(3) * t112 * t125 / (t99 ^ 2 * t93 + 0.1e1)) + (((-t134 * t97 + t152) * MDP(24) + (-t134 * t98 - t153) * MDP(25)) * t89 * t138 * t161 + ((-g(3) * t85 + t134 * t84) * MDP(24) + (g(3) * t84 + t134 * t85) * MDP(25) + 0.2e1 * ((t152 / 0.2e1 + t158 * t97) * MDP(24) + (-t153 / 0.2e1 + t158 * t98) * MDP(25)) / t91 * (-0.2e1 * t109 * t139 * t146 + t84 * t97 + t85 * t98) * t159) * t160) * t122;];
taug = t1;
