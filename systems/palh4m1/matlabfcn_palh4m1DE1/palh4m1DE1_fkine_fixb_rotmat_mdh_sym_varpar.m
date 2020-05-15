% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh4m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% T_c_mdh [4x4x(9+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   8:  mdh base (link 0) -> mdh frame (8-1), link (8-1)
%   ...
%   9+1:  mdh base (link 0) -> mdh frame (9)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 22:26
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh4m1DE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 21:58:19
% EndTime: 2020-04-11 21:58:19
% DurationCPUTime: 0.78s
% Computational Cost: add. (18057->132), mult. (19980->156), div. (1605->9), fcn. (5446->16), ass. (0->258)
unknown=NaN(40,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = (pkin(7) + 0);
t4 = cos(qJ(5));
t5 = sin(qJ(5));
t6 = pkin(1) * t5;
t8 = 0.2e1 * pkin(2) * t6;
t9 = pkin(1) ^ 2;
t18 = (-t8 + t9 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3))) * (-t8 + t9 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3)));
t19 = sqrt(-t18);
t22 = t6 - pkin(2);
t23 = pkin(2) ^ 2;
t24 = qJ(2) + pkin(6);
t25 = t24 ^ 2;
t26 = pkin(3) ^ 2;
t27 = -t8 + t9 + t23 + t25 - t26;
t29 = -pkin(1) * t19 * t4 - t27 * t22;
t30 = t29 * t1;
t31 = 0.1e1 / t24;
t32 = -t8 + t9 + t23;
t33 = 0.1e1 / t32;
t36 = t4 * pkin(1);
t38 = -t19 * t22 + t27 * t36;
t39 = t38 ^ 2;
t40 = 0.1e1 / t25;
t42 = t32 ^ 2;
t43 = 0.1e1 / t42;
t45 = t29 ^ 2;
t49 = sqrt(t43 * t40 * t39 + t43 * t40 * t45);
t50 = 0.1e1 / t49;
t51 = t50 * t33 * t31;
t52 = t51 * t30;
t53 = t38 * t1;
t54 = t51 * t53;
t55 = pkin(9) * t1;
t57 = t29 * t2;
t58 = t51 * t57;
t59 = t38 * t2;
t60 = t51 * t59;
t61 = pkin(9) * t2;
t64 = t50 * t33;
t65 = t64 * t31 * t38;
t67 = t64 * t31 * t29;
t68 = t64 * t30;
t69 = t68 - t55 + 0;
t70 = t64 * t57;
t71 = t70 - t61 + 0;
t73 = t50 * t33 * t38;
t74 = t73 + pkin(7) + 0;
t75 = t33 * t40;
t76 = t75 * t53;
t77 = 0.1e1 / pkin(3);
t79 = 0.1e1 / t26;
t82 = t8 - t9 - t23 + t25 + t26;
t83 = t82 ^ 2;
t87 = sqrt(-t18 * t79 * t40 + t79 * t40 * t83);
t88 = 0.1e1 / t87;
t90 = t88 * t19 * t77 * t50;
t92 = t75 * t30;
t95 = t88 * t77 * t82 * t50;
t97 = -t90 * t76 + t95 * t92;
t100 = -t95 * t76 - t90 * t92;
t101 = t75 * t59;
t103 = t75 * t57;
t105 = -t90 * t101 + t95 * t103;
t108 = -t95 * t101 - t90 * t103;
t110 = t33 * t40 * t29;
t113 = t33 * t40 * t38;
t115 = t90 * t110 + t95 * t113;
t118 = t95 * t110 - t90 * t113;
t119 = cos(qJ(3));
t121 = sin(qJ(3));
t123 = t121 * t100 + t119 * t97;
t126 = t119 * t100 - t121 * t97;
t127 = pkin(4) * t97;
t131 = t119 * t105 + t121 * t108;
t134 = -t121 * t105 + t119 * t108;
t135 = pkin(4) * t105;
t139 = t119 * t115 + t121 * t118;
t142 = -t121 * t115 + t119 * t118;
t143 = pkin(4) * t115;
t145 = cos(qJ(4));
t147 = sin(qJ(4));
t167 = t5 * t1;
t168 = t4 * t1;
t169 = pkin(8) * t1;
t171 = t5 * t2;
t172 = t4 * t2;
t173 = pkin(8) * t2;
t175 = t4 * pkin(2);
t178 = -t5 * pkin(2) + pkin(1);
t179 = -t8 + t9 + t23 - t25 + t26;
t181 = -t19 * t175 + t179 * t178;
t186 = t179 * t175 + t19 * t178;
t187 = t186 ^ 2;
t190 = t181 ^ 2;
t194 = sqrt(t43 * t79 * t187 + t43 * t79 * t190);
t196 = 0.1e1 / t194 * t33 * t77;
unknown(1,1) = 1;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(2,1) = 0;
unknown(2,2) = 1;
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 1;
unknown(3,4) = 0;
unknown(4,1) = 0;
unknown(4,2) = 0;
unknown(4,3) = 0;
unknown(4,4) = 1;
unknown(5,1) = t1;
unknown(5,2) = -t2;
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(6,1) = t2;
unknown(6,2) = t1;
unknown(6,3) = 0;
unknown(6,4) = 0;
unknown(7,1) = 0;
unknown(7,2) = 0;
unknown(7,3) = 1;
unknown(7,4) = t3;
unknown(8,1) = 0;
unknown(8,2) = 0;
unknown(8,3) = 0;
unknown(8,4) = 1;
unknown(9,1) = t52;
unknown(9,2) = -t54;
unknown(9,3) = t2;
unknown(9,4) = (-t55 + 0);
unknown(10,1) = t58;
unknown(10,2) = -t60;
unknown(10,3) = -t1;
unknown(10,4) = (-t61 + 0);
unknown(11,1) = t65;
unknown(11,2) = t67;
unknown(11,3) = 0;
unknown(11,4) = t3;
unknown(12,1) = 0;
unknown(12,2) = 0;
unknown(12,3) = 0;
unknown(12,4) = 1;
unknown(13,1) = t54;
unknown(13,2) = -t2;
unknown(13,3) = t52;
unknown(13,4) = t69;
unknown(14,1) = t60;
unknown(14,2) = t1;
unknown(14,3) = t58;
unknown(14,4) = t71;
unknown(15,1) = -t67;
unknown(15,2) = 0;
unknown(15,3) = t65;
unknown(15,4) = t74;
unknown(16,1) = 0;
unknown(16,2) = 0;
unknown(16,3) = 0;
unknown(16,4) = 1;
unknown(17,1) = t97;
unknown(17,2) = t100;
unknown(17,3) = t2;
unknown(17,4) = t69;
unknown(18,1) = t105;
unknown(18,2) = t108;
unknown(18,3) = -t1;
unknown(18,4) = t71;
unknown(19,1) = t115;
unknown(19,2) = t118;
unknown(19,3) = 0;
unknown(19,4) = t74;
unknown(20,1) = 0;
unknown(20,2) = 0;
unknown(20,3) = 0;
unknown(20,4) = 1;
unknown(21,1) = t123;
unknown(21,2) = t126;
unknown(21,3) = t2;
unknown(21,4) = (t127 + t68 - t55 + 0);
unknown(22,1) = t131;
unknown(22,2) = t134;
unknown(22,3) = -t1;
unknown(22,4) = (t135 + t70 - t61 + 0);
unknown(23,1) = t139;
unknown(23,2) = t142;
unknown(23,3) = 0;
unknown(23,4) = (t143 + t73 + pkin(7) + 0);
unknown(24,1) = 0;
unknown(24,2) = 0;
unknown(24,3) = 0;
unknown(24,4) = 1;
unknown(25,1) = (-t145 * t126 - t147 * t2);
unknown(25,2) = (t147 * t126 - t145 * t2);
unknown(25,3) = t123;
unknown(25,4) = (pkin(5) * t123 + t127 - t55 + t68 + 0);
unknown(26,1) = (t147 * t1 - t145 * t134);
unknown(26,2) = (t145 * t1 + t147 * t134);
unknown(26,3) = t131;
unknown(26,4) = (pkin(5) * t131 + t135 - t61 + t70 + 0);
unknown(27,1) = -(t145 * t142);
unknown(27,2) = (t147 * t142);
unknown(27,3) = t139;
unknown(27,4) = (pkin(5) * t139 + pkin(7) + t143 + t73 + 0);
unknown(28,1) = 0;
unknown(28,2) = 0;
unknown(28,3) = 0;
unknown(28,4) = 1;
unknown(29,1) = -t167;
unknown(29,2) = -t168;
unknown(29,3) = t2;
unknown(29,4) = (t169 + 0);
unknown(30,1) = -t171;
unknown(30,2) = -t172;
unknown(30,3) = -t1;
unknown(30,4) = (t173 + 0);
unknown(31,1) = t4;
unknown(31,2) = -t5;
unknown(31,3) = 0;
unknown(31,4) = t3;
unknown(32,1) = 0;
unknown(32,2) = 0;
unknown(32,3) = 0;
unknown(32,4) = 1;
unknown(33,1) = (-t196 * t186 * t100 + t196 * t181 * t97);
unknown(33,2) = (t196 * t181 * t100 + t196 * t186 * t97);
unknown(33,3) = t2;
unknown(33,4) = (-pkin(3) * t97 - t55 + t68 + 0);
unknown(34,1) = (t196 * t181 * t105 - t196 * t186 * t108);
unknown(34,2) = (t196 * t186 * t105 + t196 * t181 * t108);
unknown(34,3) = -t1;
unknown(34,4) = (-pkin(3) * t105 - t61 + t70 + 0);
unknown(35,1) = (t196 * t181 * t115 - t196 * t186 * t118);
unknown(35,2) = (t196 * t186 * t115 + t196 * t181 * t118);
unknown(35,3) = 0;
unknown(35,4) = (-pkin(3) * t115 + pkin(7) + t73 + 0);
unknown(36,1) = 0;
unknown(36,2) = 0;
unknown(36,3) = 0;
unknown(36,4) = 1;
unknown(37,1) = -t167;
unknown(37,2) = -t168;
unknown(37,3) = t2;
unknown(37,4) = (-pkin(1) * t167 + t169 + 0);
unknown(38,1) = -t171;
unknown(38,2) = -t172;
unknown(38,3) = -t1;
unknown(38,4) = (-pkin(1) * t171 + t173 + 0);
unknown(39,1) = t4;
unknown(39,2) = -t5;
unknown(39,3) = 0;
unknown(39,4) = (t36 + pkin(7) + 0);
unknown(40,1) = 0;
unknown(40,2) = 0;
unknown(40,3) = 0;
unknown(40,4) = 1;
T_ges = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,9+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,9+1]); end % symbolisch
for i = 1:9+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
