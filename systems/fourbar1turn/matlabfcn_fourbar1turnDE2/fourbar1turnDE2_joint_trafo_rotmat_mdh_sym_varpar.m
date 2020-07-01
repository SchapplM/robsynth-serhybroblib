% Calculate homogenous joint transformation matrices for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(6+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = fourbar1turnDE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:43:08
% EndTime: 2020-06-27 16:43:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (623->30), mult. (874->35), div. (72->3), fcn. (246->15), ass. (0->38)
t88 = -pkin(3) - pkin(4);
t87 = -pkin(3) + pkin(4);
t69 = cos(qJ(2));
t86 = pkin(2) * t69;
t79 = pkin(1) * t86;
t66 = -0.2e1 * t79;
t76 = pkin(1) ^ 2;
t81 = t66 + t76;
t60 = sqrt(-((pkin(2) - t88) * (pkin(2) + t88) + t81) * ((pkin(2) - t87) * (pkin(2) + t87) + t81));
t67 = sin(qJ(2));
t85 = t60 * t67;
t80 = pkin(2) ^ 2 + t76;
t78 = t66 + t80;
t63 = 0.1e1 / t78;
t72 = 0.1e1 / pkin(4);
t84 = t63 * t72;
t74 = 0.1e1 / pkin(3);
t83 = t63 * t74;
t82 = t72 * t74;
t73 = pkin(3) ^ 2;
t77 = -t73 + t80;
t71 = pkin(4) ^ 2;
t70 = cos(qJ(1));
t68 = sin(qJ(1));
t65 = pkin(1) * t69 - pkin(2);
t64 = pkin(1) - t86;
t62 = t66 + t71 + t77;
t61 = -t71 + t73 + t78;
t59 = atan2(t60 * t82, (t71 - t77 + 0.2e1 * t79) * t82);
t58 = cos(t59);
t57 = sin(t59);
t56 = atan2((pkin(1) * t67 * t61 - t65 * t60) * t83, (-pkin(1) * t85 - t65 * t61) * t83);
t55 = cos(t56);
t54 = sin(t56);
t53 = atan2((pkin(2) * t67 * t62 + t64 * t60) * t84 / 0.2e1, -(-pkin(2) * t85 + t64 * t62) * t84 / 0.2e1);
t52 = cos(t53);
t51 = sin(t53);
t1 = [t70, -t68, 0, 0; t68, t70, 0, 0; 0, 0, 1, pkin(5); t69, -t67, 0, 0; 0, 0, -1, 0; t67, t69, 0, 0; -t55, t54, 0, pkin(2); -t54, -t55, 0, 0; 0, 0, 1, 0; t52, -t51, 0, pkin(1); 0, 0, -1, 0; t51, t52, 0, 0; t58, -t57, 0, pkin(3); t57, t58, 0, 0; 0, 0, 1, 0; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
