% Calculate vector of inverse dynamics joint torques for
% fourbar2DE2
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [3x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see fourbar2DE2_invdynJ_fixb_regmin2vec.m
% MDP [3x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar2DE2_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:59
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbar2DE2_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1), zeros(3,1)}
assert(isreal(MDP) && all(size(MDP) == [3 1]), ...
  'fourbar2DE2_invdynJ_fixb_mdp_slag_vr: MDP has to be [3x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:58:59
% EndTime: 2020-06-26 17:58:59
% DurationCPUTime: 0.01s
% Computational Cost: add. (2->2), mult. (3->3), div. (0->0), fcn. (3->3), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3);];
tauJ = t1;
