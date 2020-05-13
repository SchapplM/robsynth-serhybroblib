% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
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
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = fivebar1OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:49
% EndTime: 2020-04-27 06:12:49
% DurationCPUTime: 0.11s
% Computational Cost: add. (51->24), mult. (12->8), div. (0->0), fcn. (36->10), ass. (0->19)
t12 = qJ(3) + qJ(4);
t22 = pkin(1) + 0;
t14 = sin(qJ(3));
t21 = t14 * pkin(3) + 0;
t15 = sin(qJ(1));
t20 = t15 * pkin(2) + 0;
t17 = cos(qJ(1));
t19 = t17 * pkin(2) + 0;
t16 = cos(qJ(3));
t18 = t16 * pkin(3) + t22;
t13 = qJ(1) + qJ(2);
t7 = qJ(5) + t12;
t6 = cos(t13);
t5 = cos(t12);
t4 = sin(t13);
t3 = sin(t12);
t2 = cos(t7);
t1 = sin(t7);
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t15, 0, 0; t15, t17, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t6, t4, 0, t19; -t4, -t6, 0, t20; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, t22; t14, t16, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t5, -t3, 0, t18; t3, t5, 0, t21; 0, 0, 1, 0; 0, 0, 0, 1; -t2, t1, 0, pkin(4) * t5 + t18; -t1, -t2, 0, pkin(4) * t3 + t21; 0, 0, 1, 0; 0, 0, 0, 1; -t6, t4, 0, -pkin(5) * t6 + t19; -t4, -t6, 0, -pkin(5) * t4 + t20; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t8;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
