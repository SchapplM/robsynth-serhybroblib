% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh4m1TE
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
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh4m1TE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 21:35:23
% EndTime: 2020-04-11 21:35:24
% DurationCPUTime: 0.24s
% Computational Cost: add. (5979->129), mult. (6435->171), div. (441->4), fcn. (1846->10), ass. (0->240)
unknown=NaN(40,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = (pkin(7) + 0);
t4 = cos(qJ(5));
t5 = sin(qJ(5));
t6 = t5 * pkin(1);
t8 = 0.2e1 * t6 * pkin(2);
t9 = pkin(1) ^ 2;
t19 = sqrt(-(-t8 + t9 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3))) * (-t8 + t9 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3))));
t22 = t6 - pkin(2);
t23 = pkin(2) ^ 2;
t24 = qJ(2) + pkin(6);
t25 = t24 ^ 2;
t26 = pkin(3) ^ 2;
t27 = -t8 + t9 + t23 + t25 - t26;
t29 = -t4 * t19 * pkin(1) - t22 * t27;
t30 = t1 * t29;
t31 = 0.1e1 / t24;
t33 = 0.1e1 / (-t8 + t9 + t23);
t34 = t31 * t33;
t36 = t30 * t34 / 0.2e1;
t38 = pkin(1) * t4;
t40 = -t22 * t19 + t38 * t27;
t41 = t1 * t40;
t43 = t41 * t34 / 0.2e1;
t44 = t1 * pkin(9);
t46 = t2 * t29;
t48 = t46 * t34 / 0.2e1;
t49 = t2 * t40;
t51 = t49 * t34 / 0.2e1;
t52 = t2 * pkin(9);
t56 = t40 * t31 * t33 / 0.2e1;
t59 = t29 * t31 * t33 / 0.2e1;
t61 = t30 * t33 / 0.2e1;
t62 = t61 - t44 + 0;
t64 = t46 * t33 / 0.2e1;
t65 = t64 - t52 + 0;
t67 = t40 * t33 / 0.2e1;
t68 = t67 + pkin(7) + 0;
t69 = 0.1e1 / t25;
t70 = t41 * t69;
t71 = 0.1e1 / pkin(3);
t72 = t33 * t71;
t73 = t72 * t19;
t75 = t30 * t69;
t78 = t33 * (t8 - t9 - t23 + t25 + t26) * t71;
t80 = -t70 * t73 + t75 * t78;
t83 = -t70 * t78 - t75 * t73;
t84 = t49 * t69;
t86 = t46 * t69;
t88 = -t84 * t73 + t86 * t78;
t91 = -t86 * t73 - t84 * t78;
t92 = t29 * t69;
t94 = t40 * t69;
t96 = t92 * t73 + t94 * t78;
t99 = -t94 * t73 + t92 * t78;
t100 = cos(qJ(3));
t102 = sin(qJ(3));
t104 = t80 * t100 / 0.4e1 + t83 * t102 / 0.4e1;
t107 = -t80 * t102 / 0.4e1 + t83 * t100 / 0.4e1;
t108 = t80 * pkin(4) / 0.4e1;
t112 = t88 * t100 / 0.4e1 + t91 * t102 / 0.4e1;
t115 = -t88 * t102 / 0.4e1 + t91 * t100 / 0.4e1;
t116 = t88 * pkin(4) / 0.4e1;
t120 = t96 * t100 / 0.4e1 + t99 * t102 / 0.4e1;
t123 = -t96 * t102 / 0.4e1 + t99 * t100 / 0.4e1;
t124 = t96 * pkin(4) / 0.4e1;
t126 = cos(qJ(4));
t128 = sin(qJ(4));
t148 = t1 * t5;
t149 = t1 * t4;
t150 = t1 * pkin(8);
t152 = t2 * t5;
t153 = t2 * t4;
t154 = t2 * pkin(8);
t156 = pkin(2) * t4;
t159 = -pkin(2) * t5 + pkin(1);
t160 = -t8 + t9 + t23 - t25 + t26;
t162 = -t156 * t19 + t159 * t160;
t167 = t156 * t160 + t159 * t19;
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
unknown(9,1) = t36;
unknown(9,2) = -t43;
unknown(9,3) = t2;
unknown(9,4) = (-t44 + 0);
unknown(10,1) = t48;
unknown(10,2) = -t51;
unknown(10,3) = -t1;
unknown(10,4) = (-t52 + 0);
unknown(11,1) = t56;
unknown(11,2) = t59;
unknown(11,3) = 0;
unknown(11,4) = t3;
unknown(12,1) = 0;
unknown(12,2) = 0;
unknown(12,3) = 0;
unknown(12,4) = 1;
unknown(13,1) = t43;
unknown(13,2) = -t2;
unknown(13,3) = t36;
unknown(13,4) = t62;
unknown(14,1) = t51;
unknown(14,2) = t1;
unknown(14,3) = t48;
unknown(14,4) = t65;
unknown(15,1) = -t59;
unknown(15,2) = 0;
unknown(15,3) = t56;
unknown(15,4) = t68;
unknown(16,1) = 0;
unknown(16,2) = 0;
unknown(16,3) = 0;
unknown(16,4) = 1;
unknown(17,1) = (t80 / 0.4e1);
unknown(17,2) = (t83 / 0.4e1);
unknown(17,3) = t2;
unknown(17,4) = t62;
unknown(18,1) = (t88 / 0.4e1);
unknown(18,2) = (t91 / 0.4e1);
unknown(18,3) = -t1;
unknown(18,4) = t65;
unknown(19,1) = (t96 / 0.4e1);
unknown(19,2) = (t99 / 0.4e1);
unknown(19,3) = 0;
unknown(19,4) = t68;
unknown(20,1) = 0;
unknown(20,2) = 0;
unknown(20,3) = 0;
unknown(20,4) = 1;
unknown(21,1) = t104;
unknown(21,2) = t107;
unknown(21,3) = t2;
unknown(21,4) = (t108 + t61 - t44 + 0);
unknown(22,1) = t112;
unknown(22,2) = t115;
unknown(22,3) = -t1;
unknown(22,4) = (t116 + t64 - t52 + 0);
unknown(23,1) = t120;
unknown(23,2) = t123;
unknown(23,3) = 0;
unknown(23,4) = (t124 + t67 + pkin(7) + 0);
unknown(24,1) = 0;
unknown(24,2) = 0;
unknown(24,3) = 0;
unknown(24,4) = 1;
unknown(25,1) = (-t107 * t126 - t2 * t128);
unknown(25,2) = (t107 * t128 - t2 * t126);
unknown(25,3) = t104;
unknown(25,4) = (t104 * pkin(5) + t108 - t44 + t61 + 0);
unknown(26,1) = (t1 * t128 - t115 * t126);
unknown(26,2) = (t1 * t126 + t115 * t128);
unknown(26,3) = t112;
unknown(26,4) = (t112 * pkin(5) + t116 - t52 + t64 + 0);
unknown(27,1) = -(t123 * t126);
unknown(27,2) = (t123 * t128);
unknown(27,3) = t120;
unknown(27,4) = (t120 * pkin(5) + pkin(7) + t124 + t67 + 0);
unknown(28,1) = 0;
unknown(28,2) = 0;
unknown(28,3) = 0;
unknown(28,4) = 1;
unknown(29,1) = -t148;
unknown(29,2) = -t149;
unknown(29,3) = t2;
unknown(29,4) = (t150 + 0);
unknown(30,1) = -t152;
unknown(30,2) = -t153;
unknown(30,3) = -t1;
unknown(30,4) = (t154 + 0);
unknown(31,1) = t4;
unknown(31,2) = -t5;
unknown(31,3) = 0;
unknown(31,4) = t3;
unknown(32,1) = 0;
unknown(32,2) = 0;
unknown(32,3) = 0;
unknown(32,4) = 1;
unknown(33,1) = (t80 * t162 * t72 / 0.8e1 - t83 * t167 * t72 / 0.8e1);
unknown(33,2) = (t80 * t167 * t72 / 0.8e1 + t83 * t162 * t72 / 0.8e1);
unknown(33,3) = t2;
unknown(33,4) = (-t80 * pkin(3) / 0.4e1 + t61 - t44 + 0);
unknown(34,1) = (t88 * t162 * t72 / 0.8e1 - t91 * t167 * t72 / 0.8e1);
unknown(34,2) = (t88 * t167 * t72 / 0.8e1 + t91 * t162 * t72 / 0.8e1);
unknown(34,3) = -t1;
unknown(34,4) = (-t88 * pkin(3) / 0.4e1 + t64 - t52 + 0);
unknown(35,1) = (t96 * t162 * t72 / 0.8e1 - t99 * t167 * t72 / 0.8e1);
unknown(35,2) = (t96 * t167 * t72 / 0.8e1 + t99 * t162 * t72 / 0.8e1);
unknown(35,3) = 0;
unknown(35,4) = (-t96 * pkin(3) / 0.4e1 + t67 + pkin(7) + 0);
unknown(36,1) = 0;
unknown(36,2) = 0;
unknown(36,3) = 0;
unknown(36,4) = 1;
unknown(37,1) = -t148;
unknown(37,2) = -t149;
unknown(37,3) = t2;
unknown(37,4) = (-t148 * pkin(1) + t150 + 0);
unknown(38,1) = -t152;
unknown(38,2) = -t153;
unknown(38,3) = -t1;
unknown(38,4) = (-t152 * pkin(1) + t154 + 0);
unknown(39,1) = t4;
unknown(39,2) = -t5;
unknown(39,3) = 0;
unknown(39,4) = (t38 + pkin(7) + 0);
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
