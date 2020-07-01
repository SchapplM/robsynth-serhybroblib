% Calculate minimal parameter regressor of Coriolis joint torque vector for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1OL_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:43
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1OL_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1OL_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:43:19
% EndTime: 2020-06-26 17:43:19
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->7), mult. (28->14), div. (0->0), fcn. (8->2), ass. (0->10)
t18 = pkin(2) * qJD(2);
t12 = qJD(1) + qJD(2);
t17 = pkin(2) * qJD(1) * t12;
t16 = t12 * t18;
t15 = qJD(1) * t18;
t14 = cos(qJ(2));
t13 = sin(qJ(2));
t11 = t14 * t15;
t10 = t13 * t15;
t1 = [(t13 * t16 + t10) * MDP(5) + (t14 * t16 + t11) * MDP(6); (-t13 * t17 + t10) * MDP(5) + (-t14 * t17 + t11) * MDP(6); 0; 0;];
tauc = t1;
