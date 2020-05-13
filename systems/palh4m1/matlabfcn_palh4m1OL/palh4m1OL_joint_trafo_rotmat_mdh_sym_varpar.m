% Calculate homogenous joint transformation matrices for
% palh4m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% T_mdh [4x4x9]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh4m1OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:04:07
% EndTime: 2020-04-11 23:04:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (17->17), mult. (0->0), div. (0->0), fcn. (28->14), ass. (0->158)
unknown=NaN(36,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = cos(qJ(2));
t4 = sin(qJ(2));
t5 = sin(qJ(4));
t6 = cos(qJ(4));
t7 = cos(qJ(5));
t8 = sin(qJ(5));
t9 = cos(qJ(6));
t10 = sin(qJ(6));
t11 = sin(qJ(7));
t12 = cos(qJ(7));
t13 = cos(qJ(8));
t14 = sin(qJ(8));
unknown(1,1) = t1;
unknown(1,2) = -t2;
unknown(1,3) = 0.0e0;
unknown(1,4) = 0.0e0;
unknown(2,1) = t2;
unknown(2,2) = t1;
unknown(2,3) = 0.0e0;
unknown(2,4) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = 0.1e1;
unknown(3,4) = pkin(5);
unknown(4,1) = 0.0e0;
unknown(4,2) = 0.0e0;
unknown(4,3) = 0.0e0;
unknown(4,4) = 0.1e1;
unknown(5,1) = t3;
unknown(5,2) = -t4;
unknown(5,3) = 0.0e0;
unknown(5,4) = -pkin(7);
unknown(6,1) = 0.0e0;
unknown(6,2) = 0.0e0;
unknown(6,3) = -0.1e1;
unknown(6,4) = 0.0e0;
unknown(7,1) = t4;
unknown(7,2) = t3;
unknown(7,3) = 0.0e0;
unknown(7,4) = 0.0e0;
unknown(8,1) = 0.0e0;
unknown(8,2) = 0.0e0;
unknown(8,3) = 0.0e0;
unknown(8,4) = 0.1e1;
unknown(9,1) = 0.0e0;
unknown(9,2) = 0.0e0;
unknown(9,3) = 0.1e1;
unknown(9,4) = qJ(3);
unknown(10,1) = -0.1e1;
unknown(10,2) = 0.0e0;
unknown(10,3) = 0.0e0;
unknown(10,4) = 0.0e0;
unknown(11,1) = 0.0e0;
unknown(11,2) = -0.1e1;
unknown(11,3) = 0.0e0;
unknown(11,4) = 0.0e0;
unknown(12,1) = 0.0e0;
unknown(12,2) = 0.0e0;
unknown(12,3) = 0.0e0;
unknown(12,4) = 0.1e1;
unknown(13,1) = -t5;
unknown(13,2) = -t6;
unknown(13,3) = 0.0e0;
unknown(13,4) = 0.0e0;
unknown(14,1) = 0.0e0;
unknown(14,2) = 0.0e0;
unknown(14,3) = -0.1e1;
unknown(14,4) = 0.0e0;
unknown(15,1) = t6;
unknown(15,2) = -t5;
unknown(15,3) = 0.0e0;
unknown(15,4) = 0.0e0;
unknown(16,1) = 0.0e0;
unknown(16,2) = 0.0e0;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.1e1;
unknown(17,1) = t7;
unknown(17,2) = -t8;
unknown(17,3) = 0.0e0;
unknown(17,4) = pkin(3);
unknown(18,1) = t8;
unknown(18,2) = t7;
unknown(18,3) = 0.0e0;
unknown(18,4) = 0.0e0;
unknown(19,1) = 0.0e0;
unknown(19,2) = 0.0e0;
unknown(19,3) = 0.1e1;
unknown(19,4) = 0.0e0;
unknown(20,1) = 0.0e0;
unknown(20,2) = 0.0e0;
unknown(20,3) = 0.0e0;
unknown(20,4) = 0.1e1;
unknown(21,1) = 0.0e0;
unknown(21,2) = 0.0e0;
unknown(21,3) = 0.1e1;
unknown(21,4) = pkin(4);
unknown(22,1) = -t9;
unknown(22,2) = t10;
unknown(22,3) = 0.0e0;
unknown(22,4) = 0.0e0;
unknown(23,1) = -t10;
unknown(23,2) = -t9;
unknown(23,3) = 0.0e0;
unknown(23,4) = 0.0e0;
unknown(24,1) = 0.0e0;
unknown(24,2) = 0.0e0;
unknown(24,3) = 0.0e0;
unknown(24,4) = 0.1e1;
unknown(25,1) = -t11;
unknown(25,2) = -t12;
unknown(25,3) = 0.0e0;
unknown(25,4) = pkin(6);
unknown(26,1) = 0.0e0;
unknown(26,2) = 0.0e0;
unknown(26,3) = -0.1e1;
unknown(26,4) = 0.0e0;
unknown(27,1) = t12;
unknown(27,2) = -t11;
unknown(27,3) = 0.0e0;
unknown(27,4) = 0.0e0;
unknown(28,1) = 0.0e0;
unknown(28,2) = 0.0e0;
unknown(28,3) = 0.0e0;
unknown(28,4) = 0.1e1;
unknown(29,1) = -t13;
unknown(29,2) = t14;
unknown(29,3) = 0.0e0;
unknown(29,4) = -pkin(2);
unknown(30,1) = -t14;
unknown(30,2) = -t13;
unknown(30,3) = 0.0e0;
unknown(30,4) = 0.0e0;
unknown(31,1) = 0.0e0;
unknown(31,2) = 0.0e0;
unknown(31,3) = 0.1e1;
unknown(31,4) = 0.0e0;
unknown(32,1) = 0.0e0;
unknown(32,2) = 0.0e0;
unknown(32,3) = 0.0e0;
unknown(32,4) = 0.1e1;
unknown(33,1) = 0.1e1;
unknown(33,2) = 0.0e0;
unknown(33,3) = 0.0e0;
unknown(33,4) = pkin(1);
unknown(34,1) = 0.0e0;
unknown(34,2) = 0.1e1;
unknown(34,3) = 0.0e0;
unknown(34,4) = 0.0e0;
unknown(35,1) = 0.0e0;
unknown(35,2) = 0.0e0;
unknown(35,3) = 0.1e1;
unknown(35,4) = 0.0e0;
unknown(36,1) = 0.0e0;
unknown(36,2) = 0.0e0;
unknown(36,3) = 0.0e0;
unknown(36,4) = 0.1e1;
T_ges = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,9);             % numerisch
else,                         T_mdh = sym('xx', [4,4,9]); end % symbolisch

for i = 1:9
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
