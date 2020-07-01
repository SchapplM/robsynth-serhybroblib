% Calculate vector of inverse dynamics joint torques for
% fourbar1OL
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [12x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see fourbar1OL_invdynJ_fixb_regmin2vec.m
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1OL_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:43
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbar1OL_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1), zeros(9,1)}
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1OL_invdynJ_fixb_mdp_slag_vr: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:43:20
% EndTime: 2020-06-26 17:43:20
% DurationCPUTime: 0.01s
% Computational Cost: add. (9->9), mult. (12->12), div. (0->0), fcn. (12->12), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6); RV(10) * MDP(7) + RV(11) * MDP(8) + RV(12) * MDP(9); 0;];
tauJ = t1;
