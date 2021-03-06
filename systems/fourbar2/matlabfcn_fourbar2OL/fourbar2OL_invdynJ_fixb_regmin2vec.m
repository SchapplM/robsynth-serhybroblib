% Convert inverse dynamics regressor matrix to a vector for
% fourbar2OL
% Use sparsity of the regressor matrix: 12/(4*9) elements are non-zero
%
% Input:
% RM [4x9]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [12x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 18:03
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = fourbar2OL_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(3, 7); RM(3, 8); RM(3, 9);];
RV = t1;
