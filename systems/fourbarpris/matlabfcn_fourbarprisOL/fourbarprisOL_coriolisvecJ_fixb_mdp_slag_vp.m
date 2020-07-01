% Calculate minimal parameter regressor of Coriolis joint torque vector for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisOL_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbarprisOL_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:58
% EndTime: 2020-06-27 17:29:58
% DurationCPUTime: 0.01s
% Computational Cost: add. (2->2), mult. (12->6), div. (0->0), fcn. (0->0), ass. (0->2)
t4 = (MDP(6) * qJ(2) + MDP(5));
t1 = [2 * t4 * qJD(2) * qJD(1); -t4 * qJD(1) ^ 2; 0; 0;];
tauc = t1;
