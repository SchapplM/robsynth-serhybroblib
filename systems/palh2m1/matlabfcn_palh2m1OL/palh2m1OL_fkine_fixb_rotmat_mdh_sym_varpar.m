% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = palh2m1OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:58:22
% EndTime: 2020-05-02 23:58:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (109->39), mult. (61->38), div. (0->0), fcn. (102->10), ass. (0->27)
t15 = sin(qJ(1));
t12 = qJ(2) + qJ(3);
t9 = qJ(4) + t12;
t4 = sin(t9);
t29 = t15 * t4;
t18 = cos(qJ(1));
t28 = t18 * t4;
t17 = cos(qJ(2));
t6 = t17 * pkin(2) + pkin(1);
t13 = sin(qJ(5));
t27 = t15 * t13;
t16 = cos(qJ(5));
t26 = t15 * t16;
t25 = t18 * t13;
t24 = t18 * t16;
t11 = pkin(5) + 0;
t8 = cos(t12);
t3 = pkin(3) * t8 + t6;
t23 = t15 * t3 + 0;
t22 = t18 * t3 + 0;
t5 = cos(t9);
t21 = pkin(4) * t5 + pkin(6) * t4;
t14 = sin(qJ(2));
t20 = -pkin(2) * t14 + t11;
t7 = sin(t12);
t19 = -pkin(3) * t7 + t20;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t15, 0, 0; t15, t18, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t18 * t17, -t18 * t14, -t15, pkin(1) * t18 + 0; t15 * t17, -t15 * t14, t18, pkin(1) * t15 + 0; -t14, -t17, 0, t11; 0, 0, 0, 1; t18 * t8, -t18 * t7, -t15, t18 * t6 + 0; t15 * t8, -t15 * t7, t18, t15 * t6 + 0; -t7, -t8, 0, t20; 0, 0, 0, 1; t18 * t5, -t28, -t15, t22; t15 * t5, -t29, t18, t23; -t4, -t5, 0, t19; 0, 0, 0, 1; t24 * t5 - t27, -t25 * t5 - t26, t28, t18 * t21 + t22; t26 * t5 + t25, -t27 * t5 + t24, t29, t15 * t21 + t23; -t4 * t16, t4 * t13, t5, -pkin(4) * t4 + pkin(6) * t5 + t19; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
