% Calculate joint inertia matrix for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnDE1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_inertiaJ_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnDE1_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:27
% EndTime: 2020-06-27 16:35:29
% DurationCPUTime: 0.46s
% Computational Cost: add. (3844->66), mult. (5354->135), div. (258->14), fcn. (1498->8), ass. (0->54)
t141 = pkin(4) ^ 2;
t140 = pkin(3) ^ 2;
t127 = 2 * pkin(2);
t139 = -2 * pkin(2);
t138 = (-pkin(3) - pkin(4));
t137 = (-pkin(3) + pkin(4));
t108 = sin(qJ(2));
t136 = pkin(2) * t108;
t109 = cos(qJ(2));
t135 = pkin(2) * t109;
t134 = t108 * pkin(1);
t117 = pkin(1) ^ 2;
t129 = -0.2e1 * pkin(1) * t135 + t117;
t97 = ((pkin(2) - t138) * (pkin(2) + t138)) + t129;
t98 = ((pkin(2) - t137) * (pkin(2) + t137)) + t129;
t118 = sqrt(-t97 * t98);
t126 = pkin(2) * t134;
t133 = 0.1e1 / t118 * (-t97 - t98) * t126;
t116 = pkin(2) ^ 2;
t103 = t116 + t129;
t102 = 0.1e1 / t103 ^ 2;
t132 = t102 / t141;
t131 = t108 * t118;
t130 = t109 * t118;
t128 = t140 - t141;
t125 = t102 * t136;
t101 = 0.1e1 / t103;
t111 = 0.1e1 / pkin(4);
t100 = t103 - t128;
t104 = pkin(1) - t135;
t89 = -pkin(2) * t131 + t104 * t100;
t120 = t89 ^ 2;
t95 = t100 * t136;
t90 = t104 * t118 + t95;
t86 = t90 ^ 2;
t81 = (t120 + t86) * t132;
t124 = t101 * t111 * t81 ^ (-0.1e1 / 0.2e1);
t114 = 0.1e1 / pkin(3);
t105 = pkin(1) * t109 - pkin(2);
t99 = t103 + t128;
t88 = -pkin(1) * t131 - t105 * t99;
t119 = t88 ^ 2;
t96 = t99 * t134;
t91 = -t105 * t118 + t96;
t87 = t91 ^ 2;
t123 = t101 * t114 * ((t119 + t87) / t140 * t102) ^ (-0.1e1 / 0.2e1);
t107 = t108 ^ 2;
t85 = 0.1e1 / t120;
t84 = 0.1e1 / t119;
t77 = (t108 * t91 - t109 * t88) * t123;
t76 = (-t108 * t88 - t109 * t91) * t123;
t75 = 0.1e1 + (((t117 * t107 * t127 - t105 * t133) * t101 + ((t109 * t99 + t131) * t101 - 0.2e1 * t91 * t125) * pkin(1)) / t88 - (t96 * t101 + (-t101 * t130 + ((t105 * t139 - t133) * t101 + t88 * t102 * t139) * t108) * pkin(1)) * t91 * t84) * pkin(3) * t103 * t114 / (t87 * t84 + 0.1e1);
t74 = 0.2e1 * (-((t104 * t133 + (t109 * t100 + t131) * pkin(2)) * t101 / 0.2e1 + (t116 * t107 * t101 - t90 * t125) * pkin(1)) / t89 - (-(t95 + (-t108 * t133 - t130) * pkin(2)) * t101 / 0.2e1 + (-t101 * t104 + t102 * t89) * t126) * t90 * t85) * pkin(4) * t103 * t111 / (t86 * t85 + 0.1e1);
t1 = [t77 ^ 2 * MDP(14) + t107 * MDP(4) + MDP(1) + (MDP(11) * t76 + 0.2e1 * MDP(12) * t77) * t76 + (-0.2e1 * MDP(20) * t89 * t90 + MDP(19) * t86) / t81 * t132 + 0.2e1 * (-MDP(24) * t89 - MDP(25) * t90) * pkin(1) * t124 + (0.2e1 * MDP(5) * t108 + (MDP(17) * t77 - MDP(18) * t76) * t127) * t109; t108 * MDP(6) + t109 * MDP(7) + (MDP(13) * t76 + MDP(15) * t77) * t75 + (MDP(21) * t90 - MDP(22) * t89) * t74 * t124; t74 ^ 2 * MDP(23) + MDP(8) + (MDP(16) * t75 + (-MDP(17) * t88 + MDP(18) * t91) * t123 * t127) * t75;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t1(1), t1(2); t1(2), t1(3);];
Mq = res;
