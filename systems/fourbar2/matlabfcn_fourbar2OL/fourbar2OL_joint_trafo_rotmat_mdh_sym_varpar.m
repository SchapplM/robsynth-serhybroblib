% Calculate homogenous joint transformation matrices for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = fourbar2OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:12
% EndTime: 2020-04-24 20:32:12
% DurationCPUTime: 0.03s
% Computational Cost: add. (6->6), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t23 = cos(qJ(1));
t22 = cos(qJ(2));
t21 = cos(qJ(3));
t20 = cos(qJ(4));
t19 = sin(qJ(1));
t18 = sin(qJ(2));
t17 = sin(qJ(3));
t16 = sin(qJ(4));
t1 = [t23, -t19, 0, 0; t19, t23, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t22, t18, 0, pkin(2); -t18, -t22, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t17, 0, pkin(1); t17, t21, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t16, 0, pkin(1); t16, t20, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(2); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
