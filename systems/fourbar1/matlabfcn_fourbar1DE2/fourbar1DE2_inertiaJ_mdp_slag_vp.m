% Calculate joint inertia matrix for
% fourbar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1DE2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:38
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1DE2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_inertiaJ_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1DE2_inertiaJ_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:38:32
% EndTime: 2020-06-26 17:38:33
% DurationCPUTime: 0.21s
% Computational Cost: add. (191->56), mult. (316->84), div. (10->4), fcn. (48->4), ass. (0->34)
t84 = pkin(3) - pkin(4);
t83 = pkin(3) + pkin(4);
t81 = sin(qJ(1));
t111 = pkin(1) * t81;
t82 = cos(qJ(1));
t110 = pkin(2) * t82;
t101 = pkin(1) * t110;
t92 = pkin(1) ^ 2;
t104 = -0.2e1 * t101 + t92;
t107 = (pkin(2) + t84) * (pkin(2) - t84);
t108 = (pkin(2) + t83) * (pkin(2) - t83);
t109 = (t104 + t108) * (t104 + t107);
t69 = pkin(1) * t82 - pkin(2);
t93 = sqrt(-t109);
t60 = t69 * t93;
t106 = t84 ^ 2 * t83 ^ 2;
t76 = t82 ^ 2;
t105 = t92 * t76;
t86 = pkin(4) ^ 2;
t90 = pkin(2) ^ 2;
t103 = -t86 / 0.2e1 + t90;
t87 = pkin(3) ^ 2;
t102 = t86 - t87;
t68 = t90 + t104;
t100 = t87 / 0.2e1 + t103;
t99 = t90 - t102;
t94 = pkin(1) * t92;
t91 = t92 ^ 2;
t89 = t90 ^ 2;
t85 = 0.2e1 * t90;
t71 = t86 + t87;
t59 = t60 + (t68 + t102) * t111;
t58 = t60 - (t99 + t104) * t111;
t1 = [MDP(1) + ((-MDP(4) * t59 ^ 2 - MDP(7) * t58 ^ 2) * t90 / t109 + (-t90 * ((0.4e1 * t100 * t105 - 0.4e1 * (t92 + t100) * t101 + t91 + (t85 + t102) * t92 + t90 * t99) * t93 + 0.2e1 * t69 * (t68 * t71 - t106) * t111) * MDP(5) - 0.2e1 * t90 * (0.6e1 * ((t90 - t87 / 0.6e1 - t86 / 0.6e1) * t92 + t89 + (-0.5e1 / 0.6e1 * t87 - 0.5e1 / 0.6e1 * t86) * t90 + t106 / 0.6e1) * t105 + 0.3e1 / 0.2e1 * t91 * t90 + (0.3e1 / 0.2e1 * t89 - t71 * t90 - t106 / 0.2e1) * t92 + t90 * t107 * t108 / 0.2e1 + (t81 * t84 * t83 * t60 - 0.3e1 * (t91 + (t85 - 0.2e1 / 0.3e1 * t87 - 0.2e1 / 0.3e1 * t86) * t92 + t89 + (-0.4e1 / 0.3e1 * t87 - 0.4e1 / 0.3e1 * t86) * t90 + t106 / 0.3e1) * t110) * pkin(1) + (-0.4e1 * (-t87 / 0.2e1 + t103) * t76 * t110 + t94 / 0.2e1) * t94) * MDP(6)) / pkin(3) / t93) / t68 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
