% Convert inverse dynamics regressor matrix to a vector for
% fourbar1turnTE
% Use sparsity of the regressor matrix: 44/(2*25) elements are non-zero
%
% Input:
% RM [2x25]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [44x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = fourbar1turnTE_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(2, 8); RM(1, 9); RM(2, 9); RM(1, 10); RM(2, 10); RM(1, 11); RM(2, 11); RM(1, 12); RM(2, 12); RM(1, 13); RM(2, 13); RM(1, 14); RM(2, 14); RM(1, 15); RM(2, 15); RM(2, 16); RM(1, 17); RM(2, 17); RM(1, 18); RM(2, 18); RM(1, 19); RM(2, 19); RM(1, 20); RM(2, 20); RM(1, 21); RM(2, 21); RM(1, 22); RM(2, 22); RM(2, 23); RM(1, 24); RM(2, 24); RM(1, 25); RM(2, 25);];
RV = t1;
