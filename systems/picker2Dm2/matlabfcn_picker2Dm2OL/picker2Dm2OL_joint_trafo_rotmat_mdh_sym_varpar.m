% Calculate homogenous joint transformation matrices for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% T_mdh [4x4x15]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = picker2Dm2OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:42
% EndTime: 2020-05-09 23:18:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (29->28), mult. (10->6), div. (0->0), fcn. (62->26), ass. (0->29)
t75 = cos(qJ(1));
t74 = cos(qJ(2));
t73 = cos(qJ(3));
t72 = cos(qJ(4));
t71 = cos(qJ(5));
t70 = cos(qJ(6));
t69 = cos(qJ(7));
t68 = cos(qJ(8));
t67 = cos(qJ(9));
t66 = sin(qJ(1));
t65 = sin(qJ(2));
t64 = sin(qJ(3));
t63 = sin(qJ(4));
t62 = sin(qJ(5));
t61 = sin(qJ(6));
t60 = sin(qJ(7));
t59 = sin(qJ(8));
t58 = sin(qJ(9));
t57 = cos(pkin(8));
t56 = cos(qJ(10));
t55 = cos(qJ(11));
t54 = cos(qJ(12));
t53 = sin(pkin(8));
t52 = sin(qJ(10));
t51 = sin(qJ(11));
t50 = sin(qJ(12));
t49 = -t53 * t62 + t57 * t71;
t48 = t53 * t71 + t57 * t62;
t1 = [-t75, t66, 0, 0; -t66, -t75, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t74, -t65, 0, pkin(1); t65, t74, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t73, t64, 0, pkin(2); -t64, -t73, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t72, -t63, 0, pkin(3); t63, t72, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t49, -t48, 0, t57 * pkin(5); t48, t49, 0, t53 * pkin(5); 0, 0, 1, 0; 0, 0, 0, 1; -t70, t61, 0, 0; -t61, -t70, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t60, t69, 0, pkin(7); -t69, t60, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t68, t59, 0, pkin(1); -t59, -t68, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t67, t58, 0, pkin(6); -t58, -t67, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t56, t52, 0, pkin(4); -t52, -t56, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t55, -t51, 0, pkin(5); t51, t55, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t54, t50, 0, pkin(6); -t50, -t54, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(3); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(1); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(2); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,15);             % numerisch
else,                         T_mdh = sym('xx', [4,4,15]); end % symbolisch

for i = 1:15
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
