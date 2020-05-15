% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh2m1DE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:35
% EndTime: 2020-05-02 23:51:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (59->31), mult. (41->24), div. (0->0), fcn. (73->8), ass. (0->18)
t15 = cos(qJ(2));
t5 = t15 * pkin(2) + pkin(1);
t10 = qJ(2) + qJ(3);
t7 = cos(t10);
t3 = pkin(3) * t7 + t5;
t9 = pkin(5) + 0;
t19 = pkin(4) + t3;
t12 = sin(qJ(2));
t18 = -pkin(2) * t12 + t9;
t6 = sin(t10);
t17 = -pkin(3) * t6 + t18;
t16 = cos(qJ(1));
t14 = cos(qJ(4));
t13 = sin(qJ(1));
t11 = sin(qJ(4));
t2 = -t11 * t13 + t14 * t16;
t1 = t11 * t16 + t13 * t14;
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t13, 0, 0; t13, t16, 0, 0; 0, 0, 1, t9; 0, 0, 0, 1; t16 * t15, -t16 * t12, -t13, pkin(1) * t16 + 0; t13 * t15, -t13 * t12, t16, pkin(1) * t13 + 0; -t12, -t15, 0, t9; 0, 0, 0, 1; t16 * t7, -t16 * t6, -t13, t16 * t5 + 0; t13 * t7, -t13 * t6, t16, t13 * t5 + 0; -t6, -t7, 0, t18; 0, 0, 0, 1; t16, 0, -t13, t16 * t3 + 0; t13, 0, t16, t13 * t3 + 0; 0, -1, 0, t17; 0, 0, 0, 1; t2, -t1, 0, t16 * t19 + 0; t1, t2, 0, t13 * t19 + 0; 0, 0, 1, pkin(6) + t17; 0, 0, 0, 1;];
T_ges = t4;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
