% Calculate vector of inverse dynamics joint torques for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnOL_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnOL_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'fourbar1turnOL_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:31
% EndTime: 2020-06-27 16:56:33
% DurationCPUTime: 0.70s
% Computational Cost: add. (261->105), mult. (627->172), div. (0->0), fcn. (438->10), ass. (0->45)
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t105 = sin(qJ(1));
t109 = cos(qJ(1));
t118 = g(1) * t109 + g(2) * t105;
t121 = pkin(2) * qJD(2) * qJD(3);
t108 = cos(qJ(2));
t127 = qJD(1) * t108;
t122 = pkin(2) * t127;
t104 = sin(qJ(2));
t128 = qJD(1) * t104;
t120 = t103 * t128;
t80 = -t107 * t127 + t120;
t101 = qJ(2) + qJ(3);
t93 = sin(t101);
t94 = cos(t101);
t142 = t103 * qJDD(2) * pkin(2) - g(3) * t93 + t107 * t121 - t118 * t94 + t80 * t122;
t129 = t103 * t108;
t115 = t104 * t107 + t129;
t81 = t115 * qJD(1);
t141 = g(3) * t94 + t103 * t121 - t118 * t93 - t81 * t122;
t96 = qJD(2) + qJD(3);
t98 = t104 ^ 2;
t138 = MDP(5) * (t108 ^ 2 - t98);
t77 = t96 * t115;
t82 = t103 * t104 - t107 * t108;
t74 = qJD(1) * t77 + qJDD(1) * t82;
t106 = cos(qJ(4));
t102 = sin(qJ(4));
t97 = t102 ^ 2;
t137 = MDP(19) * (-t106 ^ 2 + t97);
t131 = qJD(2) * t96;
t130 = qJD(3) * t96;
t125 = qJD(1) * qJD(4);
t124 = qJDD(2) * t107;
t119 = t106 * t125;
t112 = qJD(1) ^ 2;
t114 = pkin(1) * t112 + t118;
t73 = -qJDD(1) * t129 + (-qJDD(1) * t104 - t127 * t96) * t107 + t96 * t120;
t95 = qJDD(2) + qJDD(3);
t113 = t81 * t80 * MDP(11) + (-t80 * t96 + t73) * MDP(13) + (-t81 * t96 + t74) * MDP(14) + (-t80 ^ 2 + t81 ^ 2) * MDP(12) + t95 * MDP(15);
t111 = qJD(2) ^ 2;
t110 = qJD(4) ^ 2;
t76 = t96 * t82;
t1 = [t118 * MDP(3) + (-t115 * t73 - t76 * t81) * MDP(11) + (-t115 * t74 + t73 * t82 + t76 * t80 - t77 * t81) * MDP(12) + (-t115 * t95 + t76 * t96) * MDP(13) + (t77 * t96 + t82 * t95) * MDP(14) + (t111 * MDP(6) + qJDD(2) * MDP(7)) * t108 + (0.2e1 * qJD(2) * MDP(4) * t127 + qJDD(2) * MDP(6) - t111 * MDP(7)) * t104 + (t97 * MDP(18) + t98 * MDP(4) + MDP(1)) * qJDD(1) + (t110 * MDP(20) + qJDD(4) * MDP(21)) * t106 + (0.2e1 * MDP(18) * t119 + qJDD(4) * MDP(20) - t110 * MDP(21)) * t102 + 0.2e1 * ((qJDD(1) * t106 - t102 * t125) * MDP(23) + (-qJDD(1) * t102 - t119) * MDP(24)) * pkin(1) + 0.2e1 * (t102 * t106 * MDP(19) + t104 * t108 * MDP(5)) * qJDD(1) + 0.2e1 * (qJD(2) * t138 - qJD(4) * t137) * qJD(1) + (((-qJD(1) * t82 - t80) * MDP(16) - 0.2e1 * t81 * MDP(17)) * t104 * qJD(2) + (0.2e1 * t74 * MDP(16) + (-qJD(1) * t76 + qJDD(1) * t115 - t73) * MDP(17)) * t108) * pkin(2) + (-MDP(10) * t104 - MDP(16) * t94 + MDP(17) * t93 + MDP(23) * t106 - MDP(24) * t102 + MDP(9) * t108 + MDP(2)) * (g(1) * t105 - g(2) * t109); qJDD(2) * MDP(8) + t141 * MDP(16) + t142 * MDP(17) - t112 * t138 + (MDP(10) * t118 + qJDD(1) * MDP(7) - g(3) * MDP(9)) * t108 + (-t112 * t108 * MDP(4) + g(3) * MDP(10) + qJDD(1) * MDP(6) + MDP(9) * t118) * t104 + ((t103 * t130 - t107 * t95 + t128 * t80 - t124) * MDP(16) + (t103 * t95 + t107 * t130 + t128 * t81) * MDP(17)) * pkin(2) + t113; ((-t103 * t131 - t124) * pkin(2) + t141) * MDP(16) + (-pkin(2) * t107 * t131 + t142) * MDP(17) + t113; qJDD(4) * MDP(22) + t112 * t137 + (qJDD(1) * MDP(21) - g(3) * MDP(23) + MDP(24) * t114) * t106 + (-t112 * t106 * MDP(18) + qJDD(1) * MDP(20) + MDP(23) * t114 + g(3) * MDP(24)) * t102; 0;];
tau = t1;
