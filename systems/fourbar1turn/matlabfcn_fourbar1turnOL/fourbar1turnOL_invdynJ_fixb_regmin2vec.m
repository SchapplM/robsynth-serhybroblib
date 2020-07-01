% Convert inverse dynamics regressor matrix to a vector for
% fourbar1turnOL
% Use sparsity of the regressor matrix: 49/(5*24) elements are non-zero
%
% Input:
% RM [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [49x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = fourbar1turnOL_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(2, 8); RM(1, 9); RM(2, 9); RM(1, 10); RM(2, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(2, 15); RM(3, 15); RM(1, 16); RM(2, 16); RM(3, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(1, 18); RM(4, 18); RM(1, 19); RM(4, 19); RM(1, 20); RM(4, 20); RM(1, 21); RM(4, 21); RM(4, 22); RM(1, 23); RM(4, 23); RM(1, 24); RM(4, 24);];
RV = t1;
