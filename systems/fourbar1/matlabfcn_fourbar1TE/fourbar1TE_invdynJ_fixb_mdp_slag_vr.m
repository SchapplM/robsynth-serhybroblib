% Calculate vector of inverse dynamics joint torques for
% fourbar1TE
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [9x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see fourbar1TE_invdynJ_fixb_regmin2vec.m
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1TE_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbar1TE_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(9,1), zeros(9,1)}
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1TE_invdynJ_fixb_mdp_slag_vr: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:21:07
% EndTime: 2020-06-26 17:21:07
% DurationCPUTime: 0.04s
% Computational Cost: add. (8->8), mult. (9->9), div. (0->0), fcn. (9->9), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(5) * MDP(5) + RV(6) * MDP(6) + RV(7) * MDP(7) + RV(8) * MDP(8) + RV(9) * MDP(9);];
tauJ = t1;
