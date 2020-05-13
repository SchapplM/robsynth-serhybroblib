% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar1turnDE2
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
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = fourbar1turnDE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:18
% EndTime: 2020-04-12 19:33:19
% DurationCPUTime: 0.55s
% Computational Cost: add. (1761->59), mult. (2371->63), div. (171->3), fcn. (750->15), ass. (0->49)
t50 = -pkin(3) - pkin(4);
t49 = -pkin(3) + pkin(4);
t26 = cos(qJ(2));
t48 = pkin(2) * t26;
t24 = sin(qJ(2));
t22 = t24 * pkin(2);
t25 = sin(qJ(1));
t37 = pkin(1) * t48;
t21 = -0.2e1 * t37;
t33 = pkin(1) ^ 2;
t39 = t21 + t33;
t15 = sqrt(-((pkin(2) - t50) * (pkin(2) + t50) + t39) * ((pkin(2) - t49) * (pkin(2) + t49) + t39));
t28 = pkin(4) ^ 2;
t30 = pkin(3) ^ 2;
t38 = pkin(2) ^ 2 + t33;
t35 = -t30 + t38;
t17 = t21 + t28 + t35;
t19 = pkin(1) - t48;
t36 = t21 + t38;
t18 = 0.1e1 / t36;
t29 = 0.1e1 / pkin(4);
t44 = t18 * t29;
t45 = t15 * t24;
t12 = atan2((t19 * t15 + t17 * t22) * t44 / 0.2e1, -(-pkin(2) * t45 + t19 * t17) * t44 / 0.2e1);
t8 = sin(t12);
t47 = t25 * t8;
t27 = cos(qJ(1));
t46 = t27 * t8;
t31 = 0.1e1 / pkin(3);
t43 = t18 * t31;
t42 = t25 * t26;
t41 = t27 * t26;
t40 = t29 * t31;
t16 = -t28 + t30 + t36;
t20 = pkin(1) * t26 - pkin(2);
t13 = qJ(2) + atan2((pkin(1) * t24 * t16 - t20 * t15) * t43, (-pkin(1) * t45 - t20 * t16) * t43);
t23 = pkin(5) + 0;
t34 = t22 + t23;
t11 = cos(t13);
t10 = sin(t13);
t9 = cos(t12);
t7 = t9 * pkin(4) + pkin(1);
t6 = t27 * t9;
t5 = t25 * t9;
t4 = -pkin(3) * t11 + t48;
t3 = atan2(t15 * t40, (t28 - t35 + 0.2e1 * t37) * t40) + t13;
t2 = cos(t3);
t1 = sin(t3);
t14 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t25, 0, 0; t25, t27, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t41, -t27 * t24, t25, 0; t42, -t25 * t24, -t27, 0; t24, t26, 0, t23; 0, 0, 0, 1; -t27 * t11, t27 * t10, t25, pkin(2) * t41 + 0; -t25 * t11, t25 * t10, -t27, pkin(2) * t42 + 0; -t10, -t11, 0, t34; 0, 0, 0, 1; t6, -t46, t25, t27 * pkin(1) + 0; t5, -t47, -t27, t25 * pkin(1) + 0; t8, t9, 0, t23; 0, 0, 0, 1; -t27 * t2, t27 * t1, t25, t27 * t4 + 0; -t25 * t2, t25 * t1, -t27, t25 * t4 + 0; -t1, -t2, 0, -pkin(3) * t10 + t34; 0, 0, 0, 1; t6, -t46, t25, t27 * t7 + 0; t5, -t47, -t27, t25 * t7 + 0; t8, t9, 0, t8 * pkin(4) + t23; 0, 0, 0, 1;];
T_ges = t14;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
