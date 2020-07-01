% Calculate minimal parameter regressor of Coriolis joint torque vector for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnOL_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1turnOL_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:30
% EndTime: 2020-06-27 16:56:32
% DurationCPUTime: 0.29s
% Computational Cost: add. (148->49), mult. (456->97), div. (0->0), fcn. (294->6), ass. (0->30)
t80 = sin(qJ(2));
t100 = qJD(1) * t80;
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t83 = cos(qJ(2));
t99 = qJD(1) * t83;
t64 = t79 * t100 - t82 * t99;
t73 = qJD(2) + qJD(3);
t58 = t64 * t73;
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t109 = t78 * t81 * MDP(18) - (t78 ^ 2 - t81 ^ 2) * MDP(19) - (t78 * MDP(23) + t81 * MDP(24)) * pkin(1);
t107 = -t80 * t83 * MDP(4) + (t80 ^ 2 - t83 ^ 2) * MDP(5);
t101 = pkin(2) * qJD(2);
t90 = qJD(3) * t101;
t96 = pkin(2) * t99;
t105 = t64 * t96 + t82 * t90;
t89 = t79 * t83 + t80 * t82;
t65 = t89 * qJD(1);
t104 = -t65 * t96 + t79 * t90;
t98 = qJD(2) * t80;
t97 = qJD(3) * t73;
t95 = t73 * t101;
t66 = t79 * t80 - t82 * t83;
t61 = t73 * t89;
t59 = qJD(1) * t61;
t87 = t65 * t64 * MDP(11) + (-t65 * t73 + t59) * MDP(14) + (-t64 ^ 2 + t65 ^ 2) * MDP(12);
t86 = qJD(1) ^ 2;
t60 = t73 * t66;
t1 = [(-t58 * t89 - t60 * t65) * MDP(11) + (t58 * t66 - t59 * t89 + t60 * t64 - t61 * t65) * MDP(12) + (t83 * MDP(6) - t80 * MDP(7)) * qJD(2) ^ 2 + (t81 * MDP(20) - t78 * MDP(21)) * qJD(4) ^ 2 + (t60 * MDP(13) + t61 * MDP(14)) * t73 + ((t59 * t83 - t64 * t98) * MDP(16) + (-t58 * t83 - t65 * t98) * MDP(17)) * pkin(2) + (((t61 * t83 - t66 * t98) * MDP(16) + (-t60 * t83 - t89 * t98) * MDP(17)) * pkin(2) - 0.2e1 * t107 * qJD(2) + 0.2e1 * t109 * qJD(4)) * qJD(1); t104 * MDP(16) + t105 * MDP(17) + ((t64 * t100 + t79 * t97) * MDP(16) + (t65 * t100 + t82 * t97) * MDP(17)) * pkin(2) + t107 * t86 + t87; (-t79 * t95 + t104) * MDP(16) + (-t82 * t95 + t105) * MDP(17) + t87; -t109 * t86; 0;];
tauc = t1;
