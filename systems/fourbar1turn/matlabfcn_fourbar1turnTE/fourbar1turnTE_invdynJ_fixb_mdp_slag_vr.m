% Calculate vector of inverse dynamics joint torques for
% fourbar1turnTE
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [44x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see fourbar1turnTE_invdynJ_fixb_regmin2vec.m
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnTE_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbar1turnTE_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(44,1), zeros(25,1)}
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vr: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:23:24
% EndTime: 2020-06-27 16:23:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (42->42), mult. (44->44), div. (0->0), fcn. (44->44), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(9) + RV(15) * MDP(10) + RV(17) * MDP(11) + RV(19) * MDP(12) + RV(21) * MDP(13) + RV(23) * MDP(14) + RV(25) * MDP(15) + RV(28) * MDP(17) + RV(30) * MDP(18) + RV(32) * MDP(19) + RV(34) * MDP(20) + RV(36) * MDP(21) + RV(38) * MDP(22) + RV(41) * MDP(24) + RV(43) * MDP(25); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(20) * MDP(12) + RV(22) * MDP(13) + RV(24) * MDP(14) + RV(26) * MDP(15) + RV(27) * MDP(16) + RV(29) * MDP(17) + RV(31) * MDP(18) + RV(33) * MDP(19) + RV(35) * MDP(20) + RV(37) * MDP(21) + RV(39) * MDP(22) + RV(40) * MDP(23) + RV(42) * MDP(24) + RV(44) * MDP(25);];
tauJ = t1;
