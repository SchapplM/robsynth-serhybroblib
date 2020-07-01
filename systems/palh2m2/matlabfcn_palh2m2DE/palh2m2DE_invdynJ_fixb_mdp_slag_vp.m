% Calculate vector of inverse dynamics joint torques for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2DE_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m2DE_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'palh2m2DE_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:47
% EndTime: 2020-06-30 17:56:49
% DurationCPUTime: 0.66s
% Computational Cost: add. (301->119), mult. (598->180), div. (0->0), fcn. (339->8), ass. (0->55)
t169 = 2 * qJD(1);
t121 = sin(qJ(1));
t125 = cos(qJ(1));
t102 = g(1) * t125 + g(2) * t121;
t119 = sin(qJ(3));
t120 = sin(qJ(2));
t123 = cos(qJ(3));
t124 = cos(qJ(2));
t126 = qJD(3) ^ 2;
t127 = qJD(2) ^ 2;
t168 = t102 + (qJDD(3) * t119 + t123 * t126) * pkin(5) + (qJDD(2) * t120 + t124 * t127) * pkin(4);
t113 = qJD(1) + qJD(4);
t167 = t113 ^ 2;
t115 = t120 ^ 2;
t166 = MDP(5) * (-t124 ^ 2 + t115);
t165 = g(1) * t121 - g(2) * t125;
t112 = qJDD(1) + qJDD(4);
t118 = sin(qJ(4));
t122 = cos(qJ(4));
t164 = MDP(21) * (t112 * t118 + t122 * t167) + MDP(22) * (t112 * t122 - t118 * t167);
t114 = t119 ^ 2;
t163 = MDP(13) * (-t123 ^ 2 + t114);
t162 = qJD(2) - qJD(3);
t147 = pkin(4) * t124 + pkin(1);
t103 = pkin(2) + t147;
t161 = -0.2e1 * t103;
t160 = 2 * qJDD(1);
t159 = pkin(4) * t120;
t158 = pkin(5) * t119;
t156 = MDP(4) * t124;
t128 = qJD(1) ^ 2;
t154 = t120 * t128;
t151 = MDP(12) * t123;
t144 = qJD(2) * t159;
t100 = -qJD(3) * t158 - t144;
t150 = t100 * qJD(1);
t139 = pkin(5) * t123 + t103;
t95 = pkin(3) + t139;
t149 = t95 * qJDD(1);
t148 = qJD(4) - t113;
t146 = qJD(1) * qJD(2);
t143 = -0.2e1 * t146;
t142 = t148 * t95;
t141 = t100 * t113 + t165;
t140 = qJD(4) * (-qJD(1) - t113);
t137 = t119 * t124 - t120 * t123;
t136 = t119 * t120 + t123 * t124;
t135 = t103 * t128 + t102;
t134 = pkin(1) * t128 + t102;
t133 = qJDD(1) * t161 - t165;
t132 = pkin(1) * t160 + t165;
t129 = t112 * MDP(20) + t168 * MDP(21) * t118 + (t168 * MDP(22) + (-qJD(4) * t100 + t149 + t150) * MDP(21)) * t122;
t88 = t162 * t137;
t87 = t162 * t136;
t1 = [(qJDD(1) * MDP(1)) + t165 * MDP(2) + t102 * MDP(3) + t115 * qJDD(1) * MDP(4) + (t147 * t160 + t165) * MDP(11) + t114 * qJDD(1) * MDP(12) + (t139 * t160 + 0.2e1 * t150 + t165) * MDP(19) + (pkin(1) * MDP(10) * t143 + t127 * MDP(6) + qJDD(2) * MDP(7) + MDP(9) * t132) * t124 + (qJDD(2) * MDP(6) - t127 * MDP(7) - t132 * MDP(10) + 0.2e1 * (-pkin(1) * MDP(9) - pkin(4) * MDP(11) + t156) * t146) * t120 + (t126 * MDP(14) + qJDD(3) * MDP(15) + (t143 * t159 - t133) * MDP(17) + qJD(1) * qJD(3) * MDP(18) * t161) * t123 + (qJDD(3) * MDP(14) - t126 * MDP(15) + t133 * MDP(18) + (MDP(18) * t144 + (-MDP(17) * t103 + t151) * qJD(3)) * t169) * t119 + (t119 * t123 * MDP(13) + t120 * t124 * MDP(5)) * t160 + (-qJD(2) * t166 - qJD(3) * t163) * t169 + ((t112 * t95 + t141) * MDP(21) + t95 * MDP(22) * t140) * t122 + ((-t165 + (-qJD(1) + t148) * t100) * MDP(22) + ((-qJDD(1) - t112) * MDP(22) + MDP(21) * t140) * t95) * t118 + t129; -t154 * t156 + t128 * t166 + qJDD(2) * MDP(8) + (-g(3) * t124 + t120 * t134) * MDP(9) + (g(3) * t120 + t124 * t134) * MDP(10) + t164 * t159 + (MDP(6) * t120 + MDP(7) * t124) * qJDD(1) + ((MDP(11) + MDP(19)) * t154 + (t123 * t154 + qJDD(3) * t136 + (-qJD(2) * t137 + t88) * qJD(3)) * MDP(17) + (-t119 * t154 - qJDD(3) * t137 + (-qJD(2) * t136 + t87) * qJD(3)) * MDP(18)) * pkin(4); qJDD(3) * MDP(16) + (-g(3) * t123 + t135 * t119 + (qJDD(2) * t136 + (qJD(3) * t137 + t88) * qJD(2)) * pkin(4)) * MDP(17) + (g(3) * t119 + t135 * t123 + (-qJDD(2) * t137 + (qJD(3) * t136 + t87) * qJD(2)) * pkin(4)) * MDP(18) + t164 * t158 + (MDP(14) * t119 + MDP(15) * t123) * qJDD(1) + (t163 + (pkin(5) * MDP(19) - t151) * t119) * t128; (-MDP(22) * qJD(1) * t142 + t141 * MDP(21)) * t122 + ((t100 * t148 - t149 - t165) * MDP(22) + (-MDP(21) * t142 - t100 * MDP(22)) * qJD(1)) * t118 + t129;];
tau = t1;
