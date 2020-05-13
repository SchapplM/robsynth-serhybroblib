% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = palh2m2DE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:09
% EndTime: 2020-05-03 01:06:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (49->28), mult. (44->26), div. (0->0), fcn. (80->8), ass. (0->17)
t18 = cos(qJ(2));
t8 = t18 * pkin(4) + pkin(1);
t14 = sin(qJ(2));
t7 = t14 * pkin(4) + 0;
t6 = pkin(2) + t8;
t17 = cos(qJ(3));
t4 = t17 * pkin(5) + t6;
t19 = cos(qJ(1));
t16 = cos(qJ(4));
t15 = sin(qJ(1));
t13 = sin(qJ(3));
t12 = sin(qJ(4));
t5 = pkin(5) * t13 + t7;
t3 = -t12 * t15 + t16 * t19;
t2 = t12 * t19 + t15 * t16;
t1 = pkin(3) + t4;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t15, 0, 0; t15, t19, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19 * t18, -t19 * t14, t15, pkin(1) * t19 + 0; t15 * t18, -t15 * t14, -t19, pkin(1) * t15 + 0; t14, t18, 0, 0; 0, 0, 0, 1; t19, 0, t15, t19 * t8 + 0; t15, 0, -t19, t15 * t8 + 0; 0, 1, 0, t7; 0, 0, 0, 1; t19 * t17, -t19 * t13, t15, t19 * t6 + 0; t15 * t17, -t15 * t13, -t19, t15 * t6 + 0; t13, t17, 0, t7; 0, 0, 0, 1; t19, 0, t15, t19 * t4 + 0; t15, 0, -t19, t15 * t4 + 0; 0, 1, 0, t5; 0, 0, 0, 1; t3, -t2, 0, t1 * t19 + 0; t2, t3, 0, t1 * t15 + 0; 0, 0, 1, t5; 0, 0, 0, 1;];
T_ges = t9;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
