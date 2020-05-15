% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m2DE2
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
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m2DE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:13:04
% EndTime: 2020-05-07 02:13:04
% DurationCPUTime: 0.61s
% Computational Cost: add. (1929->140), mult. (2750->145), div. (0->0), fcn. (4406->57), ass. (0->105)
t117 = sin(pkin(14));
t85 = sin(pkin(15));
t90 = cos(pkin(15));
t91 = cos(pkin(14));
t51 = -t90 * t117 + t85 * t91;
t54 = t85 * t117 + t90 * t91;
t83 = sin(qJ(2));
t88 = cos(qJ(2));
t101 = -t83 * t51 + t54 * t88;
t74 = qJ(2) + qJ(3);
t100 = pkin(16) + t74;
t73 = pkin(17) + pkin(18);
t59 = pkin(15) + t100 + t73;
t45 = atan2(-sin(t59), cos(t59));
t67 = qJ(1) + t74;
t42 = t45 + t67;
t123 = -sin(t42) / 0.2e1;
t68 = qJ(1) - t74;
t43 = t45 - t68;
t122 = -cos(t43) / 0.2e1;
t79 = sin(pkin(18));
t80 = cos(pkin(18));
t49 = t90 * t79 + t85 * t80;
t50 = -t85 * t79 + t90 * t80;
t27 = qJ(2) + atan2(t88 * t49 + t83 * t50, t83 * t49 - t88 * t50);
t25 = pkin(18) + t27;
t121 = pkin(2) * cos(t25);
t26 = pkin(17) - t27;
t120 = pkin(3) * cos(t26);
t77 = sin(pkin(16));
t78 = cos(pkin(16));
t47 = t77 * t90 + t78 * t85;
t48 = -t77 * t85 + t78 * t90;
t82 = sin(qJ(3));
t87 = cos(qJ(3));
t32 = t47 * t87 + t82 * t48;
t33 = -t82 * t47 + t48 * t87;
t19 = -t83 * t32 + t33 * t88;
t63 = sin(t73);
t64 = cos(t73);
t95 = t32 * t88 + t83 * t33;
t9 = atan2(-t19 * t63 - t95 * t64, t19 * t64 - t63 * t95);
t8 = t9 + t74;
t6 = sin(t8);
t84 = sin(qJ(1));
t119 = t84 * t6;
t89 = cos(qJ(1));
t118 = t89 * t6;
t52 = t82 * t90 + t87 * t85;
t53 = -t82 * t85 + t87 * t90;
t36 = t83 * t52 - t53 * t88;
t94 = t52 * t88 + t83 * t53;
t13 = atan2(t36 * t63 - t94 * t64, t36 * t64 + t63 * t94);
t12 = -t13 + t26;
t10 = sin(t12);
t115 = t84 * t10;
t11 = cos(t12);
t114 = t84 * t11;
t34 = t88 * t51 + t83 * t54;
t113 = t84 * t34;
t66 = cos(t74);
t112 = t84 * t66;
t81 = sin(qJ(4));
t111 = t84 * t81;
t86 = cos(qJ(4));
t110 = t84 * t86;
t109 = t89 * t10;
t108 = t89 * t11;
t107 = t89 * t34;
t106 = t89 * t66;
t105 = t89 * t81;
t104 = t89 * t86;
t72 = pkin(11) + 0;
t61 = t88 * pkin(1) + pkin(12);
t55 = t84 * t61 + 0;
t56 = t89 * t61 + 0;
t103 = t84 * pkin(12) + 0;
t102 = t89 * pkin(12) + 0;
t99 = t84 * t120 + t55;
t98 = t89 * t120 + t56;
t60 = t83 * pkin(1) + t72;
t97 = pkin(13) + t72;
t5 = t9 + t100;
t96 = -pkin(4) * t66 - pkin(9) * cos(t5);
t93 = -pkin(3) * sin(t26) + t60;
t65 = sin(t74);
t92 = -pkin(4) * t65 + t60;
t76 = qJ(1) - qJ(2);
t75 = qJ(1) + qJ(2);
t44 = t45 + t74;
t40 = cos(t42);
t39 = sin(t43);
t31 = t101 * pkin(5) - pkin(6);
t30 = t101 * t84;
t29 = t89 * t101;
t24 = cos(t27);
t23 = sin(t27);
t16 = atan2(t83 * t90 + t88 * t85, -t83 * t85 + t88 * t90) + t25;
t15 = cos(t16);
t14 = sin(t16);
t7 = cos(t8);
t3 = t13 + t5;
t2 = cos(t3);
t1 = sin(t3);
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t89, -t84, 0, 0; t84, t89, 0, 0; 0, 0, 1, t72; 0, 0, 0, 1; t89 * t88, -t89 * t83, t84, t102; t84 * t88, -t84 * t83, -t89, t103; t83, t88, 0, t72; 0, 0, 0, 1; -t106, t89 * t65, t84, t56; -t112, t84 * t65, -t89, t55; -t65, -t66, 0, t60; 0, 0, 0, 1; -t89 * t7, t118, t84, -pkin(4) * t106 + t56; -t84 * t7, t119, -t89, -pkin(4) * t112 + t55; -t6, -t7, 0, t92; 0, 0, 0, 1; -t7 * t104 + t111, t7 * t105 + t110, -t118, (t123 - t39 / 0.2e1) * pkin(10) + (t122 - t40 / 0.2e1) * pkin(8) + (-cos(t68) / 0.2e1 - cos(t67) / 0.2e1) * pkin(4) + (cos(t76) / 0.2e1 + cos(t75) / 0.2e1) * pkin(1) + t102; -t7 * t110 - t105, t7 * t111 - t104, -t119, (t122 + t40 / 0.2e1) * pkin(10) + (t123 + t39 / 0.2e1) * pkin(8) + (-sin(t67) / 0.2e1 - sin(t68) / 0.2e1) * pkin(4) + (sin(t75) / 0.2e1 + sin(t76) / 0.2e1) * pkin(1) + t103; -t6 * t86, t6 * t81, t7, -sin(t44) * pkin(8) + cos(t44) * pkin(10) + t92; 0, 0, 0, 1; t29, -t107, t84, -t89 * pkin(6) + 0; t30, -t113, -t89, -t84 * pkin(6) + 0; t34, t101, 0, t97; 0, 0, 0, 1; t89 * t24, -t89 * t23, t84, t56; t84 * t24, -t84 * t23, -t89, t55; t23, t24, 0, t60; 0, 0, 0, 1; -t108, -t109, t84, t98; -t114, -t115, -t89, t99; t10, -t11, 0, t93; 0, 0, 0, 1; -t89 * t15, t89 * t14, t84, t89 * t121 + t56; -t84 * t15, t84 * t14, -t89, t84 * t121 + t55; -t14, -t15, 0, pkin(2) * sin(t25) + t60; 0, 0, 0, 1; -t89 * t2, t89 * t1, t84, t96 * t89 + t56; -t84 * t2, t84 * t1, -t89, t96 * t84 + t55; -t1, -t2, 0, -pkin(9) * sin(t5) + t92; 0, 0, 0, 1; t29, -t107, t84, t31 * t89 + 0; t30, -t113, -t89, t31 * t84 + 0; t34, t101, 0, t34 * pkin(5) + t97; 0, 0, 0, 1; t108, t109, t84, -pkin(7) * t108 + t98; t114, t115, -t89, -pkin(7) * t114 + t99; -t10, t11, 0, pkin(7) * t10 + t93; 0, 0, 0, 1;];
T_ges = t4;
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
