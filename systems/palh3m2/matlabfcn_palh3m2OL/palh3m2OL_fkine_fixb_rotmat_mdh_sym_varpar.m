% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m2OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:16
% EndTime: 2020-05-07 04:33:16
% DurationCPUTime: 0.26s
% Computational Cost: add. (318->96), mult. (139->80), div. (0->0), fcn. (236->26), ass. (0->61)
t39 = qJ(2) + qJ(7);
t34 = pkin(15) - t39;
t28 = -qJ(8) + t34;
t16 = sin(t28);
t44 = sin(qJ(1));
t66 = t44 * t16;
t18 = cos(t28);
t65 = t44 * t18;
t40 = qJ(2) + qJ(3);
t35 = qJ(4) + t40;
t25 = sin(t35);
t64 = t44 * t25;
t41 = sin(qJ(6));
t63 = t44 * t41;
t42 = sin(qJ(5));
t62 = t44 * t42;
t46 = cos(qJ(5));
t61 = t44 * t46;
t48 = cos(qJ(1));
t60 = t48 * t16;
t59 = t48 * t18;
t58 = t48 * t25;
t57 = t48 * t41;
t56 = t48 * t42;
t55 = t48 * t46;
t47 = cos(qJ(2));
t22 = t47 * pkin(1) + pkin(12);
t32 = cos(t40);
t7 = -pkin(4) * t32 + t22;
t54 = t44 * t7 + 0;
t53 = t48 * t7 + 0;
t38 = pkin(11) + 0;
t5 = pkin(3) * cos(t34) + t22;
t33 = pkin(16) + t39;
t43 = sin(qJ(2));
t12 = t43 * pkin(1) + t38;
t52 = pkin(13) + t38;
t20 = pkin(14) + t35;
t26 = cos(t35);
t51 = -pkin(8) * t26 - pkin(10) * t25;
t50 = -pkin(3) * sin(t34) + t12;
t30 = sin(t40);
t49 = -pkin(4) * t30 + t12;
t45 = cos(qJ(6));
t31 = cos(t39);
t29 = sin(t39);
t27 = qJ(9) + t33;
t24 = t48 * t45;
t23 = t44 * t45;
t21 = t45 * pkin(5) - pkin(6);
t17 = cos(t27);
t15 = sin(t27);
t14 = qJ(10) + t20;
t11 = cos(t14);
t10 = sin(t14);
t9 = t48 * t22 + 0;
t8 = t44 * t22 + 0;
t6 = pkin(2) * cos(t33) + t22;
t2 = -pkin(9) * cos(t20) + t7;
t1 = -pkin(7) * t18 + t5;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t48, -t44, 0, 0; t44, t48, 0, 0; 0, 0, 1, t38; 0, 0, 0, 1; t48 * t47, -t48 * t43, t44, t48 * pkin(12) + 0; t44 * t47, -t44 * t43, -t48, t44 * pkin(12) + 0; t43, t47, 0, t38; 0, 0, 0, 1; -t48 * t32, t48 * t30, t44, t9; -t44 * t32, t44 * t30, -t48, t8; -t30, -t32, 0, t12; 0, 0, 0, 1; -t48 * t26, t58, t44, t53; -t44 * t26, t64, -t48, t54; -t25, -t26, 0, t49; 0, 0, 0, 1; -t26 * t55 + t62, t26 * t56 + t61, -t58, t48 * t51 + t53; -t26 * t61 - t56, t26 * t62 - t55, -t64, t44 * t51 + t54; -t25 * t46, t25 * t42, t26, -t25 * pkin(8) + t26 * pkin(10) + t49; 0, 0, 0, 1; t24, -t57, t44, -t48 * pkin(6) + 0; t23, -t63, -t48, -t44 * pkin(6) + 0; t41, t45, 0, t52; 0, 0, 0, 1; t48 * t31, -t48 * t29, t44, t9; t44 * t31, -t44 * t29, -t48, t8; t29, t31, 0, t12; 0, 0, 0, 1; -t59, -t60, t44, t48 * t5 + 0; -t65, -t66, -t48, t44 * t5 + 0; t16, -t18, 0, t50; 0, 0, 0, 1; -t48 * t17, t48 * t15, t44, t48 * t6 + 0; -t44 * t17, t44 * t15, -t48, t44 * t6 + 0; -t15, -t17, 0, pkin(2) * sin(t33) + t12; 0, 0, 0, 1; -t48 * t11, t48 * t10, t44, t48 * t2 + 0; -t44 * t11, t44 * t10, -t48, t44 * t2 + 0; -t10, -t11, 0, -pkin(9) * sin(t20) + t49; 0, 0, 0, 1; t24, -t57, t44, t48 * t21 + 0; t23, -t63, -t48, t44 * t21 + 0; t41, t45, 0, t41 * pkin(5) + t52; 0, 0, 0, 1; t59, t60, t44, t48 * t1 + 0; t65, t66, -t48, t44 * t1 + 0; -t16, t18, 0, pkin(7) * t16 + t50; 0, 0, 0, 1;];
T_ges = t3;
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
