% Calculate homogenous joint transformation matrices for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = palh2m1DE_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:38:53
% EndTime: 2020-06-30 17:38:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (12->9), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->12)
t31 = cos(qJ(1));
t30 = cos(qJ(2));
t29 = cos(qJ(3));
t28 = cos(qJ(4));
t27 = sin(qJ(1));
t26 = sin(qJ(2));
t25 = sin(qJ(3));
t24 = sin(qJ(4));
t23 = qJ(2) + qJ(3);
t22 = cos(t23);
t21 = sin(t23);
t1 = [t31, -t27, 0, 0; t27, t31, 0, 0; 0, 0, 1, pkin(5); t30, -t26, 0, pkin(1); 0, 0, 1, 0; -t26, -t30, 0, 0; t29, -t25, 0, pkin(2); t25, t29, 0, 0; 0, 0, 1, 0; t22, t21, 0, pkin(3); -t21, t22, 0, 0; 0, 0, 1, 0; t28, -t24, 0, pkin(4); 0, 0, -1, -pkin(6); t24, t28, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
