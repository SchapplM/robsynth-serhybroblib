% Calculate vector of inverse dynamics joint torques for
% fourbarprisTE
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [9x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see fourbarprisTE_invdynJ_fixb_regmin2vec.m
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisTE_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbarprisTE_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(9,1), zeros(9,1)}
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisTE_invdynJ_fixb_mdp_slag_vr: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:07:10
% EndTime: 2020-06-27 17:07:10
% DurationCPUTime: 0.01s
% Computational Cost: add. (8->8), mult. (9->9), div. (0->0), fcn. (9->9), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(5) * MDP(5) + RV(6) * MDP(6) + RV(7) * MDP(7) + RV(8) * MDP(8) + RV(9) * MDP(9);];
tauJ = t1;
