% Calculate homogenous joint transformation matrices for
% fourbar2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
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
% Datum: 2020-04-24 20:17
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = fourbar2TE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2TE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2TE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:17:20
% EndTime: 2020-04-24 20:17:20
% DurationCPUTime: 0.02s
% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (16->2), ass. (0->3)
t8 = cos(qJ(1));
t7 = sin(qJ(1));
t1 = [t8, -t7, 0, 0; t7, t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t8, t7, 0, pkin(2); -t7, t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, pkin(1); t7, t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, pkin(1); t7, t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(2); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
