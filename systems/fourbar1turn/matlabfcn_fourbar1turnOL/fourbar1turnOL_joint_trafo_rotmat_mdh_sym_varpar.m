% Calculate homogenous joint transformation matrices for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = fourbar1turnOL_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:44
% EndTime: 2020-04-12 19:40:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (7->7), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t33 = cos(qJ(1));
t32 = cos(qJ(2));
t31 = cos(qJ(3));
t30 = cos(qJ(4));
t29 = cos(qJ(5));
t28 = sin(qJ(1));
t27 = sin(qJ(2));
t26 = sin(qJ(3));
t25 = sin(qJ(4));
t24 = sin(qJ(5));
t1 = [t33, -t28, 0, 0; t28, t33, 0, 0; 0, 0, 1, pkin(5); 0, 0, 0, 1; t32, -t27, 0, 0; 0, 0, -1, 0; t27, t32, 0, 0; 0, 0, 0, 1; -t31, t26, 0, pkin(2); -t26, -t31, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t25, 0, pkin(1); 0, 0, -1, 0; t25, t30, 0, 0; 0, 0, 0, 1; t29, -t24, 0, pkin(3); t24, t29, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
