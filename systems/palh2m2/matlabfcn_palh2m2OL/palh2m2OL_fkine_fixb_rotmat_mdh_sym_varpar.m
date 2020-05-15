% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh2m2OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:12:11
% EndTime: 2020-05-03 01:12:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (153->40), mult. (74->44), div. (0->0), fcn. (123->12), ass. (0->31)
t22 = sin(qJ(1));
t19 = qJ(2) + qJ(3);
t16 = qJ(4) + t19;
t13 = qJ(5) + t16;
t6 = sin(t13);
t37 = t22 * t6;
t25 = cos(qJ(1));
t36 = t25 * t6;
t20 = sin(qJ(6));
t35 = t22 * t20;
t23 = cos(qJ(6));
t34 = t22 * t23;
t33 = t25 * t20;
t32 = t25 * t23;
t7 = cos(t13);
t31 = pkin(3) * t7;
t24 = cos(qJ(2));
t10 = t24 * pkin(4) + pkin(1);
t12 = cos(t16);
t15 = cos(t19);
t4 = pkin(2) * t15 + t10;
t3 = pkin(5) * t12 + t4;
t30 = t22 * t3 + 0;
t29 = t25 * t3 + 0;
t21 = sin(qJ(2));
t28 = t21 * pkin(4) + 0;
t14 = sin(t19);
t27 = pkin(2) * t14 + t28;
t11 = sin(t16);
t26 = pkin(5) * t11 + t27;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t22, 0, 0; t22, t25, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25 * t24, -t25 * t21, t22, t25 * pkin(1) + 0; t22 * t24, -t22 * t21, -t25, t22 * pkin(1) + 0; t21, t24, 0, 0; 0, 0, 0, 1; t25 * t15, -t25 * t14, t22, t25 * t10 + 0; t22 * t15, -t22 * t14, -t25, t22 * t10 + 0; t14, t15, 0, t28; 0, 0, 0, 1; t25 * t12, -t25 * t11, t22, t25 * t4 + 0; t22 * t12, -t22 * t11, -t25, t22 * t4 + 0; t11, t12, 0, t27; 0, 0, 0, 1; t25 * t7, -t36, t22, t29; t22 * t7, -t37, -t25, t30; t6, t7, 0, t26; 0, 0, 0, 1; t32 * t7 - t35, -t33 * t7 - t34, -t36, t25 * t31 + t29; t34 * t7 + t33, -t35 * t7 + t32, -t37, t22 * t31 + t30; t6 * t23, -t6 * t20, t7, t6 * pkin(3) + t26; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
