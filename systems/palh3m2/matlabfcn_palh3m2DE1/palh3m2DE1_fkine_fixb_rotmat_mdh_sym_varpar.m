% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% T_c_mdh [4x4x(12+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   9:  mdh base (link 0) -> mdh frame (9-1), link (9-1)
%   ...
%   12+1:  mdh base (link 0) -> mdh frame (12)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m2DE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:57:23
% EndTime: 2020-05-07 01:57:24
% DurationCPUTime: 0.41s
% Computational Cost: add. (342->106), mult. (489->116), div. (0->0), fcn. (697->22), ass. (0->82)
t62 = sin(qJ(1));
t63 = sin(pkin(15));
t67 = cos(pkin(15));
t68 = cos(pkin(14));
t99 = sin(pkin(14));
t27 = t68 * t63 - t99 * t67;
t28 = t99 * t63 + t68 * t67;
t61 = sin(qJ(2));
t65 = cos(qJ(2));
t8 = t65 * t27 + t61 * t28;
t108 = t8 * t62;
t66 = cos(qJ(1));
t107 = t8 * t66;
t55 = sin(pkin(18));
t57 = cos(pkin(18));
t24 = -t55 * t63 + t57 * t67;
t77 = -t61 * t27 + t28 * t65;
t106 = t62 * t77;
t105 = t66 * t77;
t53 = sin(pkin(16));
t54 = cos(pkin(16));
t21 = t53 * t67 + t54 * t63;
t22 = -t53 * t63 + t54 * t67;
t51 = pkin(17) + pkin(18);
t44 = sin(t51);
t45 = cos(t51);
t78 = t21 * t45 + t22 * t44;
t104 = t78 * t62;
t103 = t78 * t66;
t102 = t44 * t21 - t45 * t22;
t100 = cos(qJ(3));
t60 = sin(qJ(3));
t95 = t60 * t61;
t94 = t60 * t65;
t25 = -t100 * t65 + t95;
t92 = t62 * t25;
t26 = t100 * t61 + t94;
t91 = t62 * t26;
t52 = qJ(3) + qJ(2);
t46 = sin(t52);
t90 = t62 * t46;
t47 = cos(t52);
t89 = t62 * t47;
t88 = t62 * t61;
t87 = t66 * t25;
t86 = t66 * t26;
t85 = t66 * t46;
t84 = t66 * t47;
t83 = t66 * t61;
t38 = pkin(1) * t65 + pkin(12);
t23 = t55 * t67 + t57 * t63;
t56 = sin(pkin(17));
t58 = cos(pkin(17));
t3 = (t23 * t56 - t24 * t58) * pkin(3) + t38;
t82 = t3 * t62 + 0;
t81 = t3 * t66 + 0;
t50 = pkin(11) + 0;
t43 = -pkin(4) * t100 + pkin(1);
t18 = pkin(4) * t95 + t43 * t65 + pkin(12);
t80 = t18 * t62 + 0;
t79 = t18 * t66 + 0;
t36 = t61 * pkin(1) + t50;
t76 = pkin(13) + t50;
t33 = t63 * pkin(8) + t67 * pkin(10);
t34 = t67 * pkin(8) - t63 * pkin(10);
t12 = t33 * t54 + t53 * t34;
t13 = -t53 * t33 + t34 * t54;
t75 = t12 * t44 - t13 * t45;
t73 = pkin(9) * (t44 * t63 - t45 * t67);
t59 = sin(qJ(4));
t72 = t59 * t102;
t64 = cos(qJ(4));
t71 = t102 * t64;
t70 = -pkin(4) * t94 + t43 * t61 + t50;
t69 = t36 + (t23 * t58 + t24 * t56) * pkin(3);
t42 = t66 * t65;
t39 = t62 * t65;
t31 = -pkin(2) * t67 + t38;
t30 = t38 * t66 + 0;
t29 = t38 * t62 + 0;
t7 = t77 * pkin(5) - pkin(6);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t66, -t62, 0, 0; t62, t66, 0, 0; 0, 0, 1, t50; 0, 0, 0, 1; t42, -t83, t62, t66 * pkin(12) + 0; t39, -t88, -t66, t62 * pkin(12) + 0; t61, t65, 0, t50; 0, 0, 0, 1; t87, t86, t62, t30; t92, t91, -t66, t29; -t26, t25, 0, t36; 0, 0, 0, 1; t102 * t66, -t103, t62, t79; t102 * t62, -t104, -t66, t80; t78, t102, 0, t70; 0, 0, 0, 1; t62 * t59 + t66 * t71, t62 * t64 - t66 * t72, t103, t75 * t66 + t79; -t66 * t59 + t62 * t71, -t62 * t72 - t66 * t64, t104, t75 * t62 + t80; t78 * t64, -t78 * t59, -t102, t12 * t45 + t13 * t44 + t70; 0, 0, 0, 1; t105, -t107, t62, -t66 * pkin(6) + 0; t106, -t108, -t66, -t62 * pkin(6) + 0; t8, t77, 0, t76; 0, 0, 0, 1; -t66 * t24, -t66 * t23, t62, t30; -t62 * t24, -t62 * t23, -t66, t29; t23, -t24, 0, t36; 0, 0, 0, 1; -t84, t85, t62, t81; -t89, t90, -t66, t82; -t46, -t47, 0, t69; 0, 0, 0, 1; t42, -t83, t62, t31 * t66 + 0; t39, -t88, -t66, t31 * t62 + 0; t61, t65, 0, pkin(2) * t63 + t36; 0, 0, 0, 1; -t87, -t86, t62, t66 * t73 + t79; -t92, -t91, -t66, t62 * t73 + t80; t26, -t25, 0, (t44 * t67 + t45 * t63) * pkin(9) + t70; 0, 0, 0, 1; t105, -t107, t62, t7 * t66 + 0; t106, -t108, -t66, t7 * t62 + 0; t8, t77, 0, t8 * pkin(5) + t76; 0, 0, 0, 1; t84, -t85, t62, -pkin(7) * t84 + t81; t89, -t90, -t66, -pkin(7) * t89 + t82; t46, t47, 0, -t46 * pkin(7) + t69; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,12+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,12+1]); end % symbolisch
for i = 1:12+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
