% Calculate vector of inverse dynamics joint torques for
% fourbarprisDE1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [9x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see fourbarprisDE1_invdynJ_fixb_regmin2vec.m
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisDE1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbarprisDE1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(9,1), zeros(9,1)}
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisDE1_invdynJ_fixb_mdp_slag_vr: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:15:00
% EndTime: 2020-06-27 17:15:00
% DurationCPUTime: 0.01s
% Computational Cost: add. (8->8), mult. (9->9), div. (0->0), fcn. (9->9), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(5) * MDP(5) + RV(6) * MDP(6) + RV(7) * MDP(7) + RV(8) * MDP(8) + RV(9) * MDP(9);];
tauJ = t1;
