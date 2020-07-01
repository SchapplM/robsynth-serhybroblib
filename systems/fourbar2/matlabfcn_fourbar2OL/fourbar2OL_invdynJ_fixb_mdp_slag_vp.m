% Calculate vector of inverse dynamics joint torques for
% fourbar2OL
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar2OL_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 18:03
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar2OL_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar2OL_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 18:03:10
% EndTime: 2020-06-26 18:03:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (47->29), mult. (64->42), div. (0->0), fcn. (30->8), ass. (0->19)
t29 = qJD(1) + qJD(2);
t42 = qJD(1) * t29;
t41 = qJD(2) * t29;
t35 = cos(qJ(2));
t40 = qJDD(1) * t35;
t39 = pkin(2) * qJD(1) * qJD(2);
t30 = qJ(1) + qJ(2);
t26 = sin(t30);
t27 = cos(t30);
t32 = sin(qJ(2));
t38 = -g(1) * t26 + g(2) * t27 + t32 * t39;
t37 = t32 * qJDD(1) * pkin(2) - g(1) * t27 - g(2) * t26 + t35 * t39;
t36 = cos(qJ(1));
t34 = cos(qJ(3));
t33 = sin(qJ(1));
t31 = sin(qJ(3));
t28 = qJDD(1) + qJDD(2);
t25 = t28 * MDP(4);
t1 = [qJDD(1) * MDP(1) + (g(1) * t33 - g(2) * t36) * MDP(2) + (g(1) * t36 + g(2) * t33) * MDP(3) + t25 + t38 * MDP(5) + t37 * MDP(6) + ((-t28 * t35 + t32 * t41 - t40) * MDP(5) + (t28 * t32 + t35 * t41) * MDP(6)) * pkin(2); t25 + ((-t32 * t42 - t40) * pkin(2) + t38) * MDP(5) + (-t35 * pkin(2) * t42 + t37) * MDP(6); qJDD(3) * MDP(7) + (g(1) * t31 - g(2) * t34) * MDP(8) + (g(1) * t34 + g(2) * t31) * MDP(9); 0;];
tau = t1;
