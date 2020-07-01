% Calculate homogenous joint transformation matrices for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = fourbarprisOL_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:51
% EndTime: 2020-06-27 17:29:51
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->9), mult. (0->0), div. (0->0), fcn. (12->6), ass. (0->7)
t16 = cos(qJ(1));
t15 = cos(qJ(3));
t14 = cos(qJ(4));
t13 = sin(qJ(1));
t12 = sin(qJ(3));
t11 = sin(qJ(4));
t1 = [-t16, t13, 0, pkin(1); -t13, -t16, 0, 0; 0, 0, 1, 0; 0, 0, 1, qJ(2); -1, 0, 0, 0; 0, -1, 0, 0; -t15, t12, 0, 0; -t12, -t15, 0, 0; 0, 0, 1, 0; -t14, t11, 0, pkin(2); -t11, -t14, 0, 0; 0, 0, 1, 0; 0, -1, 0, 0; 0, 0, -1, 0; 1, 0, 0, pkin(3);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
