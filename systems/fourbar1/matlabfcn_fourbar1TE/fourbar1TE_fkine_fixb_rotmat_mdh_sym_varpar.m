% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   4:  mdh base (link 0) -> mdh frame (4-1), link (4-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = fourbar1TE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:46:45
% EndTime: 2020-04-24 19:46:45
% DurationCPUTime: 0.17s
% Computational Cost: add. (436->39), mult. (605->34), div. (36->4), fcn. (167->4), ass. (0->32)
t19 = cos(qJ(1));
t33 = -0.2e1 * pkin(1) * t19;
t17 = pkin(2) * t33;
t25 = pkin(2) ^ 2;
t26 = pkin(1) ^ 2;
t34 = t25 + t26;
t46 = 0.1e1 / (t17 + t34) / 0.2e1;
t29 = 0.1e1 / pkin(4) * t46;
t35 = pkin(3) ^ 2 - pkin(4) ^ 2;
t30 = t25 - t35;
t36 = t17 + t26;
t13 = t30 + t36;
t41 = -pkin(3) + pkin(4);
t42 = -pkin(3) - pkin(4);
t10 = sqrt(-((pkin(2) - t42) * (pkin(2) + t42) + t36) * ((pkin(2) - t41) * (pkin(2) + t41) + t36));
t40 = pkin(2) * t19;
t16 = -pkin(1) + t40;
t37 = t16 * t10;
t18 = sin(qJ(1));
t39 = t18 * pkin(2);
t45 = -t13 * t39 + t37;
t48 = t45 * t29;
t27 = t34 + t35;
t12 = t17 + t27;
t28 = 0.1e1 / pkin(3) * t46;
t47 = (t12 * t39 + t37) * t28;
t32 = pkin(1) + 0;
t9 = t10 * t39;
t5 = t16 * t13 + t9;
t2 = t5 * t29;
t1 = (-t16 * t12 + t9) * t28;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t18, 0, 0; t18, t19, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t1, t47, 0, 0 + t40; -t47, t1, 0, 0 + t39; 0, 0, 1, 0; 0, 0, 0, 1; t2, t48, 0, t32; -t48, t2, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t2, t48, 0, (t9 + 0.2e1 * t34 * 0 + pkin(1) * t27 + (-pkin(2) * (0.4e1 * pkin(1) * 0 + t26 - t30) + t25 * t33) * t19) * t46; -t48, t2, 0, (-t45 + 0 / t46) * t46; 0, 0, 1, 0; 0, 0, 0, 1; t2, t48, 0, t46 * t5 + t32; -t48, t2, 0, -t45 * t46 + 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
