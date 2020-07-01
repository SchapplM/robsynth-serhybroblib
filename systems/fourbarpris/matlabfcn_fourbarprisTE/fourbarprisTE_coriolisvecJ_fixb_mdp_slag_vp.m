% Calculate minimal parameter regressor of Coriolis joint torque vector for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisTE_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbarprisTE_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:07:05
% EndTime: 2020-06-27 17:07:08
% DurationCPUTime: 0.47s
% Computational Cost: add. (656->41), mult. (405->81), div. (140->12), fcn. (10->2), ass. (0->38)
t95 = (qJ(1) + pkin(3));
t109 = -pkin(2) + t95;
t86 = pkin(1) + t109;
t110 = -pkin(2) - t95;
t87 = pkin(1) + t110;
t113 = t86 * t87;
t85 = pkin(1) - t109;
t108 = t85 * t113;
t84 = pkin(1) - t110;
t107 = t84 * t108;
t127 = (-t107) ^ (-0.1e1 / 0.2e1) * MDP(4);
t97 = (qJ(1) ^ 2);
t71 = pkin(1) ^ 2 - pkin(2) ^ 2 - t97 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t124 = t71 ^ 2;
t93 = 1 / (t95 ^ 2);
t126 = t124 * t93;
t91 = (t95 ^ 2);
t92 = 1 / t95;
t125 = t124 * t92 / t91;
t123 = 1 / t85;
t122 = 1 / t87;
t77 = 1 / (t85 ^ 2);
t81 = 1 / (t87 ^ 2);
t120 = 2 * MDP(7);
t119 = t71 * t92;
t118 = t71 * t93;
t117 = t123 * t122;
t116 = t123 * t81;
t115 = (t77 * t122);
t114 = t123 * t97;
t112 = MDP(5) * qJ(1);
t111 = -2 * t122 * t95;
t106 = MDP(6) * t122 * t114;
t105 = MDP(1) + 2 * t112;
t103 = (t106 * t126) / 0.2e1 + ((t91 * t120) + (MDP(1) / 0.2e1 + t112) * t126) * t117;
t78 = 0.1e1 / t86;
t74 = 0.1e1 / t84;
t1 = [(t103 * t78 / t84 ^ 2 + (t103 / t86 ^ 2 + (((-t91 * t115 + (-t81 * t91 + t111) * t123) * t120) + ((-2 * t119 + t125) * t106) + t105 * t117 * t125 + ((t105 * t123 * t71 * t111) + ((-(t81 * t114) / 0.2e1 - ((qJ(1) * t123) + (t97 * t77) / 0.2e1) * t122) * MDP(6) + (-t116 / 0.2e1 - t115 / 0.2e1) * MDP(1) + (-t117 + (-t115 - t116) * qJ(1)) * MDP(5)) * t124) * t93) * t78) * t74 + (0.1e1 / t107 * (-t108 + (t113 + (t86 - t87) * t85) * t84) * t119 / 0.2e1 - t118) * t127 + 0.2e1 * (-t92 * t127 + ((2 * t106) + (2 * MDP(1) + 4 * t112) * t117) * t78 * t74 * t118) * t95) * qJD(1) ^ 2;];
tauc = t1;
