% Convert inverse dynamics regressor matrix to a vector for
% palh2m2DE
% Use sparsity of the regressor matrix: 46/(4*22) elements are non-zero
%
% Input:
% RM [4x22]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [46x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = palh2m2DE_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(2, 8); RM(1, 9); RM(2, 9); RM(1, 10); RM(2, 10); RM(1, 11); RM(2, 11); RM(1, 12); RM(3, 12); RM(1, 13); RM(3, 13); RM(1, 14); RM(3, 14); RM(1, 15); RM(3, 15); RM(3, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(1, 20); RM(4, 20); RM(1, 21); RM(2, 21); RM(3, 21); RM(4, 21); RM(1, 22); RM(2, 22); RM(3, 22); RM(4, 22);];
RV = t1;
