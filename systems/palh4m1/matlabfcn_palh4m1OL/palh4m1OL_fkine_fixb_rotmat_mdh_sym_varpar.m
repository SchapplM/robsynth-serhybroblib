% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh4m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
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
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh4m1OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:04:07
% EndTime: 2020-04-11 23:04:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (168->101), mult. (221->84), div. (0->0), fcn. (324->14), ass. (0->212)
unknown=NaN(40,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = (pkin(5) + 0);
t4 = cos(qJ(2));
t5 = t1 * t4;
t6 = sin(qJ(2));
t7 = t1 * t6;
t8 = t1 * pkin(7);
t10 = t2 * t4;
t11 = t2 * t6;
t12 = t2 * pkin(7);
t14 = t5 * qJ(3);
t15 = t14 - t8 + 0;
t16 = t10 * qJ(3);
t17 = t16 - t12 + 0;
t18 = t6 * qJ(3);
t19 = t18 + pkin(5) + 0;
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t34 = t4 * t20;
t35 = t6 * t22;
t40 = qJ(4) + qJ(5);
t41 = sin(t40);
t43 = cos(t40);
t45 = -t7 * t41 + t5 * t43;
t48 = -t5 * t41 - t7 * t43;
t49 = t20 * pkin(3);
t50 = t7 * t49;
t51 = t22 * pkin(3);
t52 = t5 * t51;
t56 = t10 * t43 - t11 * t41;
t59 = -t10 * t41 - t11 * t43;
t60 = t11 * t49;
t61 = t10 * t51;
t65 = t4 * t41 + t6 * t43;
t68 = t4 * t43 - t6 * t41;
t69 = t34 * pkin(3);
t70 = t35 * pkin(3);
t72 = cos(qJ(6));
t74 = sin(qJ(6));
t94 = sin(qJ(7));
t95 = t1 * t94;
t96 = cos(qJ(7));
t97 = t1 * t96;
t100 = t2 * t94;
t101 = t2 * t96;
t104 = qJ(4) + qJ(8);
t105 = sin(t104);
t107 = cos(t104);
t113 = t20 * pkin(2);
t115 = t22 * pkin(2);
t137 = -t94 * pkin(1) + pkin(6);
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
unknown(9,1) = t5;
unknown(9,2) = -t7;
unknown(9,3) = t2;
unknown(9,4) = (-t8 + 0);
unknown(10,1) = t10;
unknown(10,2) = -t11;
unknown(10,3) = -t1;
unknown(10,4) = (-t12 + 0);
unknown(11,1) = t6;
unknown(11,2) = t4;
unknown(11,3) = 0;
unknown(11,4) = t3;
unknown(12,1) = 0;
unknown(12,2) = 0;
unknown(12,3) = 0;
unknown(12,4) = 1;
unknown(13,1) = t7;
unknown(13,2) = -t2;
unknown(13,3) = t5;
unknown(13,4) = t15;
unknown(14,1) = t11;
unknown(14,2) = t1;
unknown(14,3) = t10;
unknown(14,4) = t17;
unknown(15,1) = -t4;
unknown(15,2) = 0;
unknown(15,3) = t6;
unknown(15,4) = t19;
unknown(16,1) = 0;
unknown(16,2) = 0;
unknown(16,3) = 0;
unknown(16,4) = 1;
unknown(17,1) = (-t7 * t20 + t5 * t22);
unknown(17,2) = (-t5 * t20 - t7 * t22);
unknown(17,3) = t2;
unknown(17,4) = t15;
unknown(18,1) = (t10 * t22 - t11 * t20);
unknown(18,2) = (-t10 * t20 - t11 * t22);
unknown(18,3) = -t1;
unknown(18,4) = t17;
unknown(19,1) = (t34 + t35);
unknown(19,2) = (-t6 * t20 + t4 * t22);
unknown(19,3) = 0;
unknown(19,4) = t19;
unknown(20,1) = 0;
unknown(20,2) = 0;
unknown(20,3) = 0;
unknown(20,4) = 1;
unknown(21,1) = t45;
unknown(21,2) = t48;
unknown(21,3) = t2;
unknown(21,4) = (-t50 + t52 + t14 - t8 + 0);
unknown(22,1) = t56;
unknown(22,2) = t59;
unknown(22,3) = -t1;
unknown(22,4) = (-t60 + t61 + t16 - t12 + 0);
unknown(23,1) = t65;
unknown(23,2) = t68;
unknown(23,3) = 0;
unknown(23,4) = (t69 + t70 + t18 + pkin(5) + 0);
unknown(24,1) = 0;
unknown(24,2) = 0;
unknown(24,3) = 0;
unknown(24,4) = 1;
unknown(25,1) = (-t2 * t74 - t48 * t72);
unknown(25,2) = (-t2 * t72 + t48 * t74);
unknown(25,3) = t45;
unknown(25,4) = (t45 * pkin(4) + t14 - t50 + t52 - t8 + 0);
unknown(26,1) = (t1 * t74 - t59 * t72);
unknown(26,2) = (t1 * t72 + t59 * t74);
unknown(26,3) = t56;
unknown(26,4) = (t56 * pkin(4) - t12 + t16 - t60 + t61 + 0);
unknown(27,1) = -(t68 * t72);
unknown(27,2) = (t68 * t74);
unknown(27,3) = t65;
unknown(27,4) = (t65 * pkin(4) + pkin(5) + t18 + t69 + t70 + 0);
unknown(28,1) = 0;
unknown(28,2) = 0;
unknown(28,3) = 0;
unknown(28,4) = 1;
unknown(29,1) = -t95;
unknown(29,2) = -t97;
unknown(29,3) = t2;
unknown(29,4) = (t1 * pkin(6) + 0);
unknown(30,1) = -t100;
unknown(30,2) = -t101;
unknown(30,3) = -t1;
unknown(30,4) = (t2 * pkin(6) + 0);
unknown(31,1) = t96;
unknown(31,2) = -t94;
unknown(31,3) = 0;
unknown(31,4) = t3;
unknown(32,1) = 0;
unknown(32,2) = 0;
unknown(32,3) = 0;
unknown(32,4) = 1;
unknown(33,1) = (t7 * t105 - t5 * t107);
unknown(33,2) = (t5 * t105 + t7 * t107);
unknown(33,3) = t2;
unknown(33,4) = (t7 * t113 - t5 * t115 + t14 - t8 + 0);
unknown(34,1) = (-t10 * t107 + t11 * t105);
unknown(34,2) = (t10 * t105 + t11 * t107);
unknown(34,3) = -t1;
unknown(34,4) = (-t10 * t115 + t11 * t113 - t12 + t16 + 0);
unknown(35,1) = (-t4 * t105 - t6 * t107);
unknown(35,2) = (t6 * t105 - t4 * t107);
unknown(35,3) = 0;
unknown(35,4) = (-t34 * pkin(2) - t35 * pkin(2) + pkin(5) + t18 + 0);
unknown(36,1) = 0;
unknown(36,2) = 0;
unknown(36,3) = 0;
unknown(36,4) = 1;
unknown(37,1) = -t95;
unknown(37,2) = -t97;
unknown(37,3) = t2;
unknown(37,4) = (t1 * t137 + 0);
unknown(38,1) = -t100;
unknown(38,2) = -t101;
unknown(38,3) = -t1;
unknown(38,4) = (t2 * t137 + 0);
unknown(39,1) = t96;
unknown(39,2) = -t94;
unknown(39,3) = 0;
unknown(39,4) = (t96 * pkin(1) + pkin(5) + 0);
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
