% Calculate homogenous joint transformation matrices for
% fourbar1turnTE
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
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = fourbar1turnTE_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:18:53
% EndTime: 2020-06-27 16:18:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (204->25), mult. (274->42), div. (24->3), fcn. (82->6), ass. (0->37)
t109 = -pkin(3) - pkin(4);
t108 = -pkin(3) + pkin(4);
t87 = cos(qJ(2));
t107 = pkin(2) * t87;
t100 = pkin(1) * t107;
t84 = -0.2e1 * t100;
t94 = pkin(1) ^ 2;
t102 = t84 + t94;
t77 = sqrt(-((pkin(2) - t109) * (pkin(2) + t109) + t102) * ((pkin(2) - t108) * (pkin(2) + t108) + t102));
t85 = sin(qJ(2));
t106 = t77 * t85;
t101 = pkin(2) ^ 2 + t94;
t99 = t84 + t101;
t81 = 0.1e1 / t99;
t90 = 0.1e1 / pkin(4);
t105 = t81 * t90;
t92 = 0.1e1 / pkin(3);
t104 = t81 * t92;
t103 = t90 * t92;
t91 = pkin(3) ^ 2;
t98 = -t91 + t101;
t97 = -t105 / 0.2e1;
t96 = -t104 / 0.2e1;
t95 = t103 / 0.2e1;
t89 = pkin(4) ^ 2;
t88 = cos(qJ(1));
t86 = sin(qJ(1));
t83 = pkin(1) * t87 - pkin(2);
t82 = pkin(1) - t107;
t80 = t84 + t89 + t98;
t79 = -t89 + t91 + t99;
t78 = (t89 - t98 + 0.2e1 * t100) * t95;
t76 = pkin(1) * t85 * t79 - t83 * t77;
t75 = pkin(2) * t85 * t80 + t82 * t77;
t74 = (-pkin(2) * t106 + t82 * t80) * t97;
t73 = (-pkin(1) * t106 - t83 * t79) * t96;
t1 = [t88, -t86, 0, 0; t86, t88, 0, 0; 0, 0, 1, pkin(5); t87, -t85, 0, 0; 0, 0, -1, 0; t85, t87, 0, 0; t73, t76 * t104 / 0.2e1, 0, pkin(2); t76 * t96, t73, 0, 0; 0, 0, 1, 0; t74, t75 * t97, 0, pkin(1); 0, 0, -1, 0; t75 * t105 / 0.2e1, t74, 0, 0; t78, -t77 * t103 / 0.2e1, 0, pkin(3); t77 * t95, t78, 0, 0; 0, 0, 1, 0; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
