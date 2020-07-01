% Calculate homogenous joint transformation matrices for
% fourbar1turnDE1
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
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = fourbar1turnDE1_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:30:48
% EndTime: 2020-06-27 16:30:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (623->30), mult. (874->35), div. (72->3), fcn. (246->15), ass. (0->38)
t100 = -pkin(3) - pkin(4);
t99 = -pkin(3) + pkin(4);
t81 = cos(qJ(2));
t98 = pkin(2) * t81;
t91 = pkin(1) * t98;
t78 = -0.2e1 * t91;
t88 = pkin(1) ^ 2;
t93 = t78 + t88;
t72 = sqrt(-((pkin(2) - t100) * (pkin(2) + t100) + t93) * ((pkin(2) - t99) * (pkin(2) + t99) + t93));
t79 = sin(qJ(2));
t97 = t72 * t79;
t92 = pkin(2) ^ 2 + t88;
t90 = t78 + t92;
t75 = 0.1e1 / t90;
t84 = 0.1e1 / pkin(4);
t96 = t75 * t84;
t86 = 0.1e1 / pkin(3);
t95 = t75 * t86;
t94 = t84 * t86;
t85 = pkin(3) ^ 2;
t89 = -t85 + t92;
t83 = pkin(4) ^ 2;
t82 = cos(qJ(1));
t80 = sin(qJ(1));
t77 = pkin(1) * t81 - pkin(2);
t76 = pkin(1) - t98;
t74 = t78 + t83 + t89;
t73 = -t83 + t85 + t90;
t71 = atan2(t72 * t94, (t83 - t89 + 0.2e1 * t91) * t94);
t70 = cos(t71);
t69 = sin(t71);
t68 = atan2((pkin(1) * t79 * t73 - t77 * t72) * t95, (-pkin(1) * t97 - t77 * t73) * t95);
t67 = cos(t68);
t66 = sin(t68);
t65 = atan2((pkin(2) * t79 * t74 + t76 * t72) * t96 / 0.2e1, -(-pkin(2) * t97 + t76 * t74) * t96 / 0.2e1);
t64 = cos(t65);
t63 = sin(t65);
t1 = [t82, -t80, 0, 0; t80, t82, 0, 0; 0, 0, 1, pkin(5); t81, -t79, 0, 0; 0, 0, -1, 0; t79, t81, 0, 0; -t67, t66, 0, pkin(2); -t66, -t67, 0, 0; 0, 0, 1, 0; t64, -t63, 0, pkin(1); 0, 0, -1, 0; t63, t64, 0, 0; t70, -t69, 0, pkin(3); t69, t70, 0, 0; 0, 0, 1, 0; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
