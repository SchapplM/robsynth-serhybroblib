% Calculate joint inertia matrix for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1OL_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1OL_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'palh2m1OL_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:14
% EndTime: 2020-06-30 17:45:16
% DurationCPUTime: 0.32s
% Computational Cost: add. (389->92), mult. (747->120), div. (0->0), fcn. (838->8), ass. (0->51)
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t120 = t89 * MDP(27) + t93 * MDP(28);
t118 = MDP(30) * t93;
t102 = -MDP(31) * t89 + t118;
t127 = 0.2e1 * t102;
t126 = MDP(30) * t89 + MDP(31) * t93;
t91 = sin(qJ(3));
t92 = sin(qJ(2));
t95 = cos(qJ(3));
t96 = cos(qJ(2));
t105 = t91 * t92 - t95 * t96;
t81 = t96 * pkin(2) + pkin(1);
t66 = -t105 * pkin(3) + t81;
t125 = 0.2e1 * t66;
t124 = pkin(2) * t91;
t123 = pkin(4) * t89;
t94 = cos(qJ(4));
t122 = t94 * pkin(3);
t121 = t89 * t93;
t74 = -t91 * t96 - t95 * t92;
t90 = sin(qJ(4));
t63 = -t94 * t105 + t90 * t74;
t116 = t63 * MDP(29);
t80 = t95 * pkin(2) + pkin(3);
t75 = t94 * t80;
t69 = -t90 * t124 + t75;
t115 = t69 * MDP(23);
t70 = -t94 * t124 - t90 * t80;
t114 = t70 * MDP(24);
t113 = t90 * MDP(24);
t112 = t95 * MDP(16);
t111 = MDP(26) * t121;
t87 = t89 ^ 2;
t110 = t87 * MDP(25) + MDP(22) + 0.2e1 * t111;
t109 = MDP(15) + t110;
t64 = t90 * t105 + t94 * t74;
t88 = t93 ^ 2;
t104 = ((-t87 + t88) * MDP(26) + MDP(25) * t121 + MDP(20)) * t64 + (-MDP(21) + t120) * t63;
t103 = MDP(27) * t93 - MDP(28) * t89;
t100 = (t94 * MDP(23) - t113) * pkin(3);
t99 = t74 * MDP(13) + t105 * MDP(14) + t104;
t86 = pkin(4) * t93;
t79 = -pkin(4) - t122;
t78 = t90 * pkin(3) + pkin(6);
t76 = t79 * t89;
t68 = pkin(6) - t70;
t67 = -pkin(4) - t69;
t65 = t67 * t89;
t55 = t63 * pkin(4) - t64 * pkin(6) + t66;
t1 = [MDP(1) + t92 ^ 2 * MDP(4) + 0.2e1 * t92 * t96 * MDP(5) - 0.2e1 * t81 * t105 * MDP(16) + 0.2e1 * (-t92 * MDP(10) + t96 * MDP(9)) * pkin(1) + (MDP(11) * t74 + 0.2e1 * t105 * MDP(12) + 0.2e1 * t81 * MDP(17)) * t74 + (MDP(23) * t125 + t127 * t55 + t116) * t63 + (MDP(24) * t125 + 0.2e1 * (-MDP(19) + t103) * t63 + (t88 * MDP(25) + MDP(18) - 0.2e1 * t111) * t64) * t64; -t92 * MDP(6) - t96 * MDP(7) + t99 + t126 * (-t63 * t68 + t64 * t67); MDP(8) - t67 * t127 + 0.2e1 * (-t91 * MDP(17) + t112) * pkin(2) + 0.2e1 * t115 + 0.2e1 * t114 + t109; t99 + t126 * (-t63 * t78 + t64 * t79); (t75 + t122) * MDP(23) + (t76 + t65) * MDP(31) + (-t67 - t79) * t118 + (-pkin(3) - t80) * t113 + (t112 + (-MDP(23) * t90 - MDP(24) * t94 - MDP(17)) * t91) * pkin(2) + t109; -t127 * t79 + 0.2e1 * t100 + t109; t104 + t126 * (-pkin(4) * t64 - pkin(6) * t63); t115 + t114 + (-t67 * t93 + t86) * MDP(30) + (t65 - t123) * MDP(31) + t110; (-t79 * t93 + t86) * MDP(30) + (t76 - t123) * MDP(31) + t100 + t110; pkin(4) * t127 + t110; t102 * t55 + t103 * t64 + t116; -t126 * t68 + t120; -t126 * t78 + t120; -pkin(6) * t126 + t120; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
