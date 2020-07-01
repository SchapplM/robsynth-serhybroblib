% Calculate joint inertia matrix for
% fourbar1turnTE
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
%   see fourbar1turnTE_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnTE_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaJ_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnTE_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:37
% EndTime: 2020-06-27 16:22:41
% DurationCPUTime: 0.40s
% Computational Cost: add. (2812->64), mult. (3837->133), div. (160->12), fcn. (1070->4), ass. (0->50)
t124 = pkin(4) ^ 2;
t123 = -2 * pkin(2);
t102 = pkin(1) ^ 2;
t95 = cos(qJ(2));
t117 = pkin(2) * t95;
t108 = -0.2e1 * t117;
t109 = pkin(1) * t108 + t102;
t120 = (-pkin(3) - pkin(4));
t83 = ((pkin(2) - t120) * (pkin(2) + t120)) + t109;
t119 = (-pkin(3) + pkin(4));
t84 = ((pkin(2) - t119) * (pkin(2) + t119)) + t109;
t103 = sqrt(-t83 * t84);
t94 = sin(qJ(2));
t111 = t103 * t94;
t113 = pkin(3) ^ 2 - t124;
t101 = pkin(2) ^ 2;
t89 = t101 + t109;
t86 = t89 - t113;
t90 = pkin(1) - t117;
t75 = -pkin(2) * t111 + t86 * t90;
t122 = -t75 / 0.2e1;
t121 = -t95 / 0.2e1;
t118 = pkin(2) * t94;
t116 = t94 * pkin(1);
t107 = pkin(2) * t116;
t115 = 0.1e1 / t103 * (-t83 - t84) * t107;
t87 = 0.1e1 / t89;
t97 = 0.1e1 / pkin(4);
t114 = t87 * t97;
t100 = 0.1e1 / pkin(3);
t112 = t100 * t87;
t110 = t95 * t103;
t88 = 0.1e1 / t89 ^ 2;
t106 = t88 * t118;
t93 = t94 ^ 2;
t91 = t95 * pkin(1) - pkin(2);
t85 = t89 + t113;
t82 = t85 * t116;
t81 = t86 * t118;
t77 = -t103 * t91 + t82;
t76 = t103 * t90 + t81;
t74 = -pkin(1) * t111 - t85 * t91;
t73 = t76 ^ 2;
t72 = 0.1e1 / t75 ^ 2;
t71 = 0.1e1 / t74 ^ 2;
t68 = (t94 * t77 / 0.2e1 + t74 * t121) * t112;
t67 = (-t94 * t74 / 0.2e1 + t77 * t121) * t112;
t66 = 0.1e1 + (((0.2e1 * t102 * t93 * pkin(2) - t91 * t115) * t87 + ((t95 * t85 + t111) * t87 - 0.2e1 * t77 * t106) * pkin(1)) / t74 - (t82 * t87 + (-t87 * t110 + ((t91 * t123 - t115) * t87 + t74 * t88 * t123) * t94) * pkin(1)) * t77 * t71) * pkin(3) * t100 / (t71 * t77 ^ 2 + 0.1e1) * t89;
t65 = 0.2e1 * (-((t90 * t115 + (t95 * t86 + t111) * pkin(2)) * t87 / 0.2e1 + (t101 * t93 * t87 - t76 * t106) * pkin(1)) / t75 - (-(t81 + (-t94 * t115 - t110) * pkin(2)) * t87 / 0.2e1 + (t75 * t88 - t87 * t90) * t107) * t76 * t72) * pkin(4) / (t72 * t73 + 0.1e1) * t89 * t97;
t1 = [0.2e1 * t94 * t95 * MDP(5) + t93 * MDP(4) + MDP(1) + (t73 * MDP(19) / 0.4e1 + t76 * MDP(20) * t122) / t124 * t88 + (MDP(14) * t68 + 0.2e1 * MDP(17) * t117) * t68 + (-MDP(24) * t75 - MDP(25) * t76) * pkin(1) * t114 + (MDP(11) * t67 + 0.2e1 * t68 * MDP(12) + MDP(18) * t108) * t67; t94 * MDP(6) + t95 * MDP(7) + (t67 * MDP(13) + t68 * MDP(15)) * t66 + (t76 * MDP(21) / 0.2e1 + MDP(22) * t122) * t65 * t114; MDP(23) * t65 ^ 2 + MDP(8) + (MDP(16) * t66 + (-MDP(17) * t74 + MDP(18) * t77) * pkin(2) * t112) * t66;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t1(1), t1(2); t1(2), t1(3);];
Mq = res;
