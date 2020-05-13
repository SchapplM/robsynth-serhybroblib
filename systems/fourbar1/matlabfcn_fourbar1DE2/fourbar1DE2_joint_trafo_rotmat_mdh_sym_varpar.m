% Calculate homogenous joint transformation matrices for
% fourbar1DE2
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

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = fourbar1DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:02:25
% EndTime: 2020-04-24 20:02:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (203->24), mult. (274->42), div. (24->3), fcn. (78->4), ass. (0->35)
t69 = cos(qJ(1));
t88 = pkin(2) * t69;
t83 = pkin(1) * t88;
t67 = -0.2e1 * t83;
t75 = pkin(1) ^ 2;
t84 = pkin(2) ^ 2 + t75;
t82 = t67 + t84;
t64 = 0.1e1 / t82;
t92 = -t64 / 0.2e1;
t91 = t64 / 0.2e1;
t90 = -pkin(3) - pkin(4);
t89 = -pkin(3) + pkin(4);
t85 = t67 + t75;
t60 = sqrt(-((pkin(2) - t90) * (pkin(2) + t90) + t85) * ((pkin(2) - t89) * (pkin(2) + t89) + t85));
t68 = sin(qJ(1));
t87 = t60 * t68;
t71 = 0.1e1 / pkin(4);
t73 = 0.1e1 / pkin(3);
t86 = t71 * t73;
t72 = pkin(3) ^ 2;
t81 = -t72 + t84;
t80 = t71 * t91;
t79 = t73 * t91;
t78 = t86 / 0.2e1;
t70 = pkin(4) ^ 2;
t62 = -t70 + t72 + t82;
t66 = pkin(1) * t69 - pkin(2);
t77 = pkin(1) * t68 * t62 - t66 * t60;
t63 = t67 + t70 + t81;
t65 = -pkin(1) + t88;
t76 = pkin(2) * t68 * t63 - t65 * t60;
t61 = (t70 - t81 + 0.2e1 * t83) * t78;
t59 = (pkin(2) * t87 + t65 * t63) * t80;
t58 = (pkin(1) * t87 + t66 * t62) * t79;
t1 = [t69, -t68, 0, 0; t68, t69, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, t77 * t79, 0, pkin(2); t77 * t73 * t92, t58, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t59, t76 * t71 * t92, 0, pkin(1); t76 * t80, t59, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t61, -t60 * t86 / 0.2e1, 0, pkin(3); t60 * t78, t61, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
