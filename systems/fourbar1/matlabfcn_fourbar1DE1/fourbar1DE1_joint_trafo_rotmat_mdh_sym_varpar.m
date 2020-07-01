% Calculate homogenous joint transformation matrices for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:30
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = fourbar1DE1_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:26:31
% EndTime: 2020-06-26 17:26:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (203->24), mult. (274->42), div. (24->3), fcn. (78->4), ass. (0->35)
t60 = cos(qJ(1));
t79 = pkin(2) * t60;
t74 = pkin(1) * t79;
t58 = -0.2e1 * t74;
t66 = pkin(1) ^ 2;
t75 = pkin(2) ^ 2 + t66;
t73 = t58 + t75;
t55 = 0.1e1 / t73;
t83 = -t55 / 0.2e1;
t82 = t55 / 0.2e1;
t81 = -pkin(3) - pkin(4);
t80 = -pkin(3) + pkin(4);
t76 = t58 + t66;
t51 = sqrt(-((pkin(2) - t81) * (pkin(2) + t81) + t76) * ((pkin(2) - t80) * (pkin(2) + t80) + t76));
t59 = sin(qJ(1));
t78 = t51 * t59;
t62 = 0.1e1 / pkin(4);
t64 = 0.1e1 / pkin(3);
t77 = t62 * t64;
t63 = pkin(3) ^ 2;
t72 = -t63 + t75;
t71 = t62 * t82;
t70 = t64 * t82;
t69 = t77 / 0.2e1;
t61 = pkin(4) ^ 2;
t53 = -t61 + t63 + t73;
t57 = pkin(1) * t60 - pkin(2);
t68 = pkin(1) * t59 * t53 - t57 * t51;
t54 = t58 + t61 + t72;
t56 = -pkin(1) + t79;
t67 = pkin(2) * t59 * t54 - t56 * t51;
t52 = (t61 - t72 + 0.2e1 * t74) * t69;
t50 = (pkin(2) * t78 + t56 * t54) * t71;
t49 = (pkin(1) * t78 + t57 * t53) * t70;
t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, 0; t49, t68 * t70, 0, pkin(2); t68 * t64 * t83, t49, 0, 0; 0, 0, 1, 0; t50, t67 * t62 * t83, 0, pkin(1); t67 * t71, t50, 0, 0; 0, 0, 1, 0; t52, -t51 * t77 / 0.2e1, 0, pkin(3); t51 * t69, t52, 0, 0; 0, 0, 1, 0; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
