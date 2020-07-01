% Calculate minimal parameter regressor of Coriolis joint torque vector for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2DE_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m2DE_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:47
% EndTime: 2020-06-30 17:56:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (168->63), mult. (412->109), div. (0->0), fcn. (198->6), ass. (0->42)
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t129 = pkin(1) * t90 * MDP(10) + (t87 ^ 2 - t90 ^ 2) * MDP(5);
t80 = qJD(1) + qJD(4);
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t126 = (t88 * MDP(21) - t85 * MDP(22)) * t80 ^ 2;
t124 = pkin(1) * MDP(9) - t90 * MDP(4);
t123 = qJD(2) - qJD(3);
t86 = sin(qJ(3));
t89 = cos(qJ(3));
t122 = -t86 * t89 * MDP(12) + (t86 ^ 2 - t89 ^ 2) * MDP(13);
t120 = pkin(4) * t87;
t119 = pkin(5) * t86;
t118 = t89 * pkin(5);
t79 = t90 * pkin(4);
t110 = pkin(4) * qJD(2);
t74 = -qJD(3) * t119 - t87 * t110;
t117 = t85 * t74;
t93 = qJD(1) ^ 2;
t116 = t86 * t93;
t115 = t87 * t93;
t114 = t88 * t74;
t113 = t89 * t93;
t108 = qJD(4) * t85;
t107 = qJD(4) * t88;
t106 = qJD(4) - t80;
t91 = qJD(3) ^ 2;
t92 = qJD(2) ^ 2;
t75 = t91 * t118 + t92 * t79;
t103 = qJD(1) * t114 - t74 * t107 + t85 * t75;
t76 = pkin(2) + t79 + pkin(1);
t71 = pkin(3) + t76 + t118;
t102 = t71 * t108;
t101 = t85 * t71 * qJD(1);
t98 = t90 * t86 - t89 * t87;
t97 = t87 * t86 + t90 * t89;
t94 = (-t71 * t107 - t117) * MDP(22);
t70 = t88 * t75;
t66 = t123 * t98;
t65 = t123 * t97;
t1 = [t103 * MDP(21) + (t74 * t108 + t70) * MDP(22) + (t90 * MDP(6) - t87 * MDP(7)) * t92 + (t89 * MDP(14) - t86 * MDP(15)) * t91 + ((-t102 + t114) * MDP(21) + t94) * t80 + (0.2e1 * t74 * MDP(19) - MDP(21) * t102 + t94 + 0.2e1 * ((-t86 * MDP(17) - t89 * MDP(18)) * t76 - t122) * qJD(3) + 0.2e1 * ((-pkin(4) * MDP(11) - t124) * t87 + (-t89 * MDP(17) + t86 * MDP(18)) * t120 - t129) * qJD(2)) * qJD(1); t129 * t93 + t126 * t120 + t124 * t115 + ((MDP(11) + MDP(19)) * t115 + (t87 * t113 + (-t98 * qJD(2) + t66) * qJD(3)) * MDP(17) + (-t86 * t115 + (-t97 * qJD(2) + t65) * qJD(3)) * MDP(18)) * pkin(4); (t76 * t116 + (t98 * qJD(3) + t66) * t110) * MDP(17) + (t76 * t113 + (t97 * qJD(3) + t65) * t110) * MDP(18) + pkin(5) * MDP(19) * t116 + t122 * t93 + t126 * t119; (-qJD(4) * t101 - t80 * (-t101 - t114) + t103) * MDP(21) + (t70 + t106 * t117 + (-t106 * t88 * t71 - t117) * qJD(1)) * MDP(22);];
tauc = t1;
