% Calculate vector of inverse dynamics joint torques for
% palh2m1DE
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [49x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see palh2m1DE_invdynJ_fixb_regmin2vec.m
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1DE_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = palh2m1DE_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(49,1), zeros(21,1)}
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vr: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:39
% EndTime: 2020-06-30 17:39:39
% DurationCPUTime: 0.02s
% Computational Cost: add. (45->45), mult. (49->49), div. (0->0), fcn. (49->49), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(9) + RV(15) * MDP(10) + RV(17) * MDP(11) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(26) * MDP(14) + RV(31) * MDP(16) + RV(34) * MDP(17) + RV(37) * MDP(18) + RV(40) * MDP(19) + RV(42) * MDP(20) + RV(46) * MDP(21); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(24) * MDP(13) + RV(27) * MDP(14) + RV(29) * MDP(15) + RV(32) * MDP(16) + RV(35) * MDP(17) + RV(38) * MDP(18) + RV(43) * MDP(20) + RV(47) * MDP(21); RV(19) * MDP(11) + RV(22) * MDP(12) + RV(25) * MDP(13) + RV(28) * MDP(14) + RV(30) * MDP(15) + RV(33) * MDP(16) + RV(36) * MDP(17) + RV(39) * MDP(18) + RV(44) * MDP(20) + RV(48) * MDP(21); RV(41) * MDP(19) + RV(45) * MDP(20) + RV(49) * MDP(21);];
tauJ = t1;
