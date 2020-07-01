% Convert inverse dynamics regressor matrix to a vector for
% fourbar2TE
% Use sparsity of the regressor matrix: 3/(1*3) elements are non-zero
%
% Input:
% RM [1x3]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [3x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:50
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = fourbar2TE_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3);];
RV = t1;
