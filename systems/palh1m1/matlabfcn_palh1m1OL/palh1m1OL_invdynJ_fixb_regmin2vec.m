% Convert inverse dynamics regressor matrix to a vector for
% palh1m1OL
% Use sparsity of the regressor matrix: 205/(13*66) elements are non-zero
%
% Input:
% RM [13x66]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [205x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = palh1m1OL_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(2, 8); RM(1, 9); RM(2, 9); RM(1, 10); RM(2, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(2, 15); RM(3, 15); RM(1, 16); RM(2, 16); RM(3, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(4, 19); RM(1, 20); RM(2, 20); RM(3, 20); RM(4, 20); RM(1, 21); RM(2, 21); RM(3, 21); RM(4, 21); RM(2, 22); RM(3, 22); RM(4, 22); RM(1, 23); RM(2, 23); RM(3, 23); RM(4, 23); RM(1, 24); RM(2, 24); RM(3, 24); RM(4, 24); RM(1, 25); RM(2, 25); RM(3, 25); RM(4, 25); RM(5, 25); RM(1, 26); RM(2, 26); RM(3, 26); RM(4, 26); RM(5, 26); RM(1, 27); RM(2, 27); RM(3, 27); RM(4, 27); RM(5, 27); RM(1, 28); RM(2, 28); RM(3, 28); RM(4, 28); RM(5, 28); RM(1, 29); RM(2, 29); RM(3, 29); RM(4, 29); RM(5, 29); RM(1, 30); RM(2, 30); RM(3, 30); RM(4, 30); RM(5, 30); RM(1, 31); RM(2, 31); RM(3, 31); RM(4, 31); RM(5, 31); RM(1, 32); RM(6, 32); RM(1, 33); RM(6, 33); RM(1, 34); RM(6, 34); RM(1, 35); RM(6, 35); RM(6, 36); RM(1, 37); RM(6, 37); RM(1, 38); RM(6, 38); RM(1, 39); RM(2, 39); RM(7, 39); RM(1, 40); RM(2, 40); RM(7, 40); RM(1, 41); RM(2, 41); RM(7, 41); RM(1, 42); RM(2, 42); RM(7, 42); RM(2, 43); RM(7, 43); RM(1, 44); RM(2, 44); RM(7, 44); RM(1, 45); RM(2, 45); RM(7, 45); RM(1, 46); RM(2, 46); RM(8, 46); RM(1, 47); RM(2, 47); RM(8, 47); RM(1, 48); RM(2, 48); RM(8, 48); RM(1, 49); RM(2, 49); RM(8, 49); RM(2, 50); RM(8, 50); RM(1, 51); RM(2, 51); RM(8, 51); RM(1, 52); RM(2, 52); RM(8, 52); RM(1, 53); RM(2, 53); RM(8, 53); RM(9, 53); RM(1, 54); RM(2, 54); RM(8, 54); RM(9, 54); RM(1, 55); RM(2, 55); RM(8, 55); RM(9, 55); RM(1, 56); RM(2, 56); RM(8, 56); RM(9, 56); RM(2, 57); RM(8, 57); RM(9, 57); RM(1, 58); RM(2, 58); RM(8, 58); RM(9, 58); RM(1, 59); RM(2, 59); RM(8, 59); RM(9, 59); RM(1, 60); RM(2, 60); RM(7, 60); RM(10, 60); RM(1, 61); RM(2, 61); RM(7, 61); RM(10, 61); RM(1, 62); RM(2, 62); RM(7, 62); RM(10, 62); RM(1, 63); RM(2, 63); RM(7, 63); RM(10, 63); RM(2, 64); RM(7, 64); RM(10, 64); RM(1, 65); RM(2, 65); RM(7, 65); RM(10, 65); RM(1, 66); RM(2, 66); RM(7, 66); RM(10, 66);];
RV = t1;
