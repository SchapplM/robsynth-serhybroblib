% Calculate minimal parameter regressor of Coriolis joint torque vector for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1DE_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m1DE_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:36
% EndTime: 2020-06-30 17:39:37
% DurationCPUTime: 0.48s
% Computational Cost: add. (339->91), mult. (807->160), div. (0->0), fcn. (552->6), ass. (0->70)
t115 = qJD(2) + qJD(3);
t122 = sin(qJ(3));
t126 = cos(qJ(2));
t160 = t126 * t122;
t123 = sin(qJ(2));
t125 = cos(qJ(3));
t161 = t125 * t123;
t133 = t160 + t161;
t98 = t115 * t133;
t116 = qJD(1) + qJD(4);
t124 = cos(qJ(4));
t173 = t116 * t124;
t107 = -0.2e1 * t160 * t161;
t119 = t125 ^ 2;
t120 = t126 ^ 2;
t172 = t107 + (-t122 ^ 2 + t119) * (t120 - 0.1e1 / 0.2e1);
t159 = t123 ^ 2 - t120;
t171 = t107 - t159 * (t119 - 0.1e1 / 0.2e1);
t170 = t123 * t126 * MDP(4) - t159 * MDP(5);
t163 = t123 * t122;
t106 = t125 * t126 - t163;
t167 = pkin(3) * qJD(3);
t103 = t106 * t167;
t112 = pkin(3) * t125 + pkin(2);
t135 = -pkin(3) * t163 + t112 * t126;
t130 = t106 * qJD(3);
t152 = qJD(2) * t126;
t153 = qJD(2) * t123;
t93 = t112 * t152 + (-t122 * t153 + t130) * pkin(3);
t168 = -t135 * qJD(2) - t103 + t93;
t121 = sin(qJ(4));
t102 = pkin(3) * t160 + t112 * t123;
t95 = qJD(2) * t102 + t133 * t167;
t166 = t121 * t95;
t101 = pkin(1) + pkin(4) + t135;
t165 = t101 * t124;
t113 = pkin(2) * t126 + pkin(1);
t128 = qJD(1) ^ 2;
t164 = t113 * t128;
t162 = t123 * t128;
t158 = MDP(11) * t133;
t157 = MDP(20) * t116;
t156 = MDP(21) * t116;
t155 = qJD(1) * t121;
t154 = qJD(1) * t124;
t151 = qJD(4) * t124;
t150 = t128 * MDP(18);
t149 = -qJD(1) - t116;
t148 = -qJD(4) + t116;
t146 = qJD(1) * qJD(2);
t131 = t106 * qJD(2);
t97 = t130 + t131;
t87 = qJD(2) * t93 + t97 * t167;
t92 = -t112 * t153 + (-t133 * qJD(3) - t122 * t152) * pkin(3);
t145 = t121 * t87 + t95 * t151 + t92 * t154;
t144 = -0.2e1 * MDP(12) * t128;
t143 = pkin(2) * t153;
t140 = t101 * t155;
t138 = -0.2e1 * pkin(1) * t146;
t137 = t149 * t121;
t136 = qJD(3) * (-qJD(2) - t115);
t134 = pkin(2) * qJD(2) * (-qJD(3) + t115);
t132 = -t128 * t106 * t158 + (t106 * t115 - t97) * qJD(1) * MDP(13);
t127 = qJD(2) ^ 2;
t104 = t133 * pkin(3);
t100 = t106 * t164;
t99 = t133 * t164;
t96 = pkin(3) * t131 + t103;
t86 = t124 * t87;
t1 = [(t101 * qJD(4) * t137 + t92 * t173 + t145) * MDP(20) + (t86 + t92 * t137 + (t149 * t165 - t166) * qJD(4)) * MDP(21) + 0.2e1 * t170 * t146 + (MDP(10) * t138 - t127 * MDP(6)) * t126 + (t127 * MDP(7) + MDP(9) * t138) * t123 + (-t97 * MDP(13) + t98 * MDP(14)) * t115 + (0.2e1 * t97 * t158 + 0.4e1 * (t171 * qJD(2) + t172 * qJD(3)) * MDP(12)) * qJD(1) + 0.2e1 * ((-t106 * t143 - t113 * t98) * MDP(16) + (-t113 * t97 + t133 * t143) * MDP(17) + t92 * MDP(18)) * qJD(1); pkin(1) * MDP(9) * t162 + t171 * t144 + (t99 + (t106 * t162 + t122 * t136) * pkin(2)) * MDP(16) + (t100 + (t125 * t136 - t133 * t162) * pkin(2)) * MDP(17) + t102 * t150 + (t102 * t173 + t168 * t121) * t157 + (-t116 * t121 * t102 + t168 * t124) * t156 + t132 + (pkin(1) * t126 * MDP(10) - t170) * t128; t172 * t144 + (t122 * t134 + t99) * MDP(16) + (t125 * t134 + t100) * MDP(17) + t104 * t150 + (t104 * t154 - t121 * t96 + (t121 * t97 + t133 * t151) * pkin(3)) * t157 + (-t104 * t155 - t124 * t96 + (-qJD(4) * t121 * t133 + t124 * t97) * pkin(3)) * t156 + t132; (-qJD(4) * t140 - (t124 * t95 - t140) * t116 + t145) * MDP(20) + (t86 + t148 * t166 + (-t121 * t92 + t148 * t165) * qJD(1)) * MDP(21);];
tauc = t1;
