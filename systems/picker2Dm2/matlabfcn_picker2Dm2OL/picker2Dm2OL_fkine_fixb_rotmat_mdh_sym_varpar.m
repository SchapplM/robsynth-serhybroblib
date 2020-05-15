% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% T_c_mdh [4x4x(15+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   11:  mdh base (link 0) -> mdh frame (11-1), link (11-1)
%   ...
%   15+1:  mdh base (link 0) -> mdh frame (15)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = picker2Dm2OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:42
% EndTime: 2020-05-09 23:18:42
% DurationCPUTime: 0.24s
% Computational Cost: add. (211->64), mult. (50->22), div. (0->0), fcn. (110->26), ass. (0->46)
t36 = qJ(1) + qJ(2);
t35 = qJ(1) + qJ(8);
t47 = pkin(7) + 0;
t46 = sin(pkin(8)) * pkin(5) + 0;
t45 = cos(pkin(8)) * pkin(5) + 0;
t31 = qJ(3) + t36;
t30 = qJ(4) + t36;
t29 = qJ(6) + t36;
t38 = sin(qJ(1));
t9 = -t38 * pkin(1) + 0;
t40 = cos(qJ(1));
t10 = -t40 * pkin(1) + 0;
t25 = sin(t36);
t44 = -pkin(2) * t25 + t9;
t27 = cos(t36);
t43 = -pkin(2) * t27 + t10;
t42 = -pkin(3) * t25 + t9;
t41 = -pkin(3) * t27 + t10;
t39 = cos(qJ(7));
t37 = sin(qJ(7));
t34 = pkin(8) + qJ(5);
t28 = qJ(11) + t35;
t26 = cos(t35);
t24 = sin(t35);
t23 = cos(t34);
t22 = sin(t34);
t21 = qJ(9) + t31;
t20 = qJ(10) + t30;
t19 = qJ(12) + t29;
t18 = cos(t31);
t17 = cos(t30);
t16 = cos(t29);
t15 = sin(t31);
t14 = sin(t30);
t13 = sin(t29);
t12 = cos(t28);
t11 = sin(t28);
t8 = cos(t21);
t7 = sin(t21);
t6 = cos(t20);
t5 = cos(t19);
t4 = sin(t20);
t3 = sin(t19);
t2 = pkin(6) * t18;
t1 = pkin(6) * t15;
t32 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t40, t38, 0, 0; -t38, -t40, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t27, t25, 0, t10; -t25, -t27, 0, t9; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t15, 0, t43; t15, t18, 0, t44; 0, 0, 1, 0; 0, 0, 0, 1; -t17, t14, 0, t41; -t14, -t17, 0, t42; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t22, 0, t45; t22, t23, 0, t46; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t13, 0, t10; t13, t16, 0, t9; 0, 0, 1, 0; 0, 0, 0, 1; t37, t39, 0, t47; -t39, t37, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t24, 0, t10; t24, t26, 0, t9; 0, 0, 1, 0; 0, 0, 0, 1; -t8, t7, 0, t2 + t43; -t7, -t8, 0, t1 + t44; 0, 0, 1, 0; 0, 0, 0, 1; t6, -t4, 0, -pkin(4) * t17 + t41; t4, t6, 0, -pkin(4) * t14 + t42; 0, 0, 1, 0; 0, 0, 0, 1; t12, -t11, 0, pkin(5) * t26 + t10; t11, t12, 0, pkin(5) * t24 + t9; 0, 0, 1, 0; 0, 0, 0, 1; -t5, t3, 0, pkin(6) * t16 + t10; -t3, -t5, 0, pkin(6) * t13 + t9; 0, 0, 1, 0; 0, 0, 0, 1; t37, t39, 0, t37 * pkin(3) + t47; -t39, t37, 0, -t39 * pkin(3) + 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t22, 0, t23 * pkin(1) + t45; t22, t23, 0, t22 * pkin(1) + t46; 0, 0, 1, 0; 0, 0, 0, 1; -t8, t7, 0, t2 + (-t27 - t8) * pkin(2) + t10; -t7, -t8, 0, t1 + (-t25 - t7) * pkin(2) + t9; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t32;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,15+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,15+1]); end % symbolisch
for i = 1:15+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
