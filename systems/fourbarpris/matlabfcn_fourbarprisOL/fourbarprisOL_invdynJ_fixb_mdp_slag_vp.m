% Calculate vector of inverse dynamics joint torques for
% fourbarprisOL
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisOL_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisOL_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisOL_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:58
% EndTime: 2020-06-27 17:29:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->18), mult. (42->28), div. (0->0), fcn. (16->4), ass. (0->10)
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t29 = -g(1) * t20 + g(2) * t22;
t28 = (2 * qJD(2) * qJD(1)) + t29;
t26 = -g(1) * t22 - g(2) * t20;
t25 = qJDD(2) - t26;
t23 = qJD(1) ^ 2;
t21 = cos(qJ(3));
t19 = sin(qJ(3));
t1 = [qJDD(1) * MDP(1) + t29 * MDP(2) + t26 * MDP(3) + t25 * MDP(4) + t28 * MDP(5) + (0.2e1 * qJDD(1) * MDP(5) + (qJDD(1) * qJ(2) + t28) * MDP(6)) * qJ(2); qJDD(1) * MDP(4) - t23 * MDP(5) + (-t23 * qJ(2) + t25) * MDP(6); qJDD(3) * MDP(7) + (-g(1) * t19 + g(2) * t21) * MDP(8) + (-g(1) * t21 - g(2) * t19) * MDP(9); 0;];
tau = t1;
