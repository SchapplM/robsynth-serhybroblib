% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m2TE
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m2TE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:30
% EndTime: 2020-05-07 01:41:30
% DurationCPUTime: 0.37s
% Computational Cost: add. (342->106), mult. (489->116), div. (0->0), fcn. (697->22), ass. (0->82)
t59 = sin(qJ(1));
t60 = sin(pkin(15));
t64 = cos(pkin(15));
t65 = cos(pkin(14));
t95 = sin(pkin(14));
t25 = t65 * t60 - t64 * t95;
t26 = t95 * t60 + t65 * t64;
t58 = sin(qJ(2));
t62 = cos(qJ(2));
t7 = t62 * t25 + t58 * t26;
t104 = t7 * t59;
t63 = cos(qJ(1));
t103 = t7 * t63;
t52 = sin(pkin(18));
t54 = cos(pkin(18));
t22 = -t52 * t60 + t54 * t64;
t74 = -t58 * t25 + t26 * t62;
t102 = t59 * t74;
t101 = t63 * t74;
t50 = sin(pkin(16));
t51 = cos(pkin(16));
t19 = t50 * t64 + t51 * t60;
t20 = -t50 * t60 + t51 * t64;
t48 = pkin(17) + pkin(18);
t41 = sin(t48);
t42 = cos(t48);
t100 = t41 * t19 - t42 * t20;
t4 = t19 * t42 + t20 * t41;
t98 = t59 * t4;
t97 = t63 * t4;
t96 = cos(qJ(3));
t57 = sin(qJ(3));
t91 = t57 * t58;
t90 = t57 * t62;
t23 = -t96 * t62 + t91;
t88 = t59 * t23;
t24 = t96 * t58 + t90;
t87 = t59 * t24;
t49 = qJ(3) + qJ(2);
t43 = sin(t49);
t86 = t59 * t43;
t44 = cos(t49);
t85 = t59 * t44;
t84 = t59 * t58;
t83 = t63 * t23;
t82 = t63 * t24;
t81 = t63 * t43;
t80 = t63 * t44;
t79 = t63 * t58;
t35 = pkin(1) * t62 + pkin(12);
t21 = t52 * t64 + t54 * t60;
t53 = sin(pkin(17));
t55 = cos(pkin(17));
t3 = (t21 * t53 - t22 * t55) * pkin(3) + t35;
t78 = t3 * t59 + 0;
t77 = t3 * t63 + 0;
t47 = pkin(11) + 0;
t40 = -pkin(4) * t96 + pkin(1);
t16 = pkin(4) * t91 + t40 * t62 + pkin(12);
t76 = t16 * t59 + 0;
t75 = t16 * t63 + 0;
t34 = t58 * pkin(1) + t47;
t73 = pkin(13) + t47;
t31 = t60 * pkin(8) + t64 * pkin(10);
t32 = t64 * pkin(8) - t60 * pkin(10);
t11 = t31 * t51 + t50 * t32;
t12 = -t50 * t31 + t32 * t51;
t72 = t11 * t41 - t12 * t42;
t70 = pkin(9) * (t41 * t60 - t42 * t64);
t56 = sin(qJ(4));
t69 = t56 * t100;
t61 = cos(qJ(4));
t68 = t100 * t61;
t67 = -pkin(4) * t90 + t40 * t58 + t47;
t66 = t34 + (t21 * t55 + t22 * t53) * pkin(3);
t38 = t63 * t62;
t36 = t59 * t62;
t29 = -pkin(2) * t64 + t35;
t28 = t35 * t63 + 0;
t27 = t35 * t59 + 0;
t6 = t74 * pkin(5) - pkin(6);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t63, -t59, 0, 0; t59, t63, 0, 0; 0, 0, 1, t47; 0, 0, 0, 1; t38, -t79, t59, t63 * pkin(12) + 0; t36, -t84, -t63, t59 * pkin(12) + 0; t58, t62, 0, t47; 0, 0, 0, 1; t83, t82, t59, t28; t88, t87, -t63, t27; -t24, t23, 0, t34; 0, 0, 0, 1; t100 * t63, -t97, t59, t75; t100 * t59, -t98, -t63, t76; t4, t100, 0, t67; 0, 0, 0, 1; t59 * t56 + t63 * t68, t59 * t61 - t63 * t69, t97, t72 * t63 + t75; -t63 * t56 + t59 * t68, -t59 * t69 - t63 * t61, t98, t72 * t59 + t76; t4 * t61, -t4 * t56, -t100, t11 * t42 + t12 * t41 + t67; 0, 0, 0, 1; t101, -t103, t59, -t63 * pkin(6) + 0; t102, -t104, -t63, -t59 * pkin(6) + 0; t7, t74, 0, t73; 0, 0, 0, 1; -t63 * t22, -t63 * t21, t59, t28; -t59 * t22, -t59 * t21, -t63, t27; t21, -t22, 0, t34; 0, 0, 0, 1; -t80, t81, t59, t77; -t85, t86, -t63, t78; -t43, -t44, 0, t66; 0, 0, 0, 1; t38, -t79, t59, t29 * t63 + 0; t36, -t84, -t63, t29 * t59 + 0; t58, t62, 0, pkin(2) * t60 + t34; 0, 0, 0, 1; -t83, -t82, t59, t63 * t70 + t75; -t88, -t87, -t63, t59 * t70 + t76; t24, -t23, 0, (t41 * t64 + t42 * t60) * pkin(9) + t67; 0, 0, 0, 1; t101, -t103, t59, t6 * t63 + 0; t102, -t104, -t63, t6 * t59 + 0; t7, t74, 0, t7 * pkin(5) + t73; 0, 0, 0, 1; t80, -t81, t59, -pkin(7) * t80 + t77; t85, -t86, -t63, -pkin(7) * t85 + t78; t43, t44, 0, -t43 * pkin(7) + t66; 0, 0, 0, 1;];
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
