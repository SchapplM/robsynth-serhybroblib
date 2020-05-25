% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbarprisDE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:08:37
% EndTime: 2020-05-07 09:08:38
% DurationCPUTime: 0.13s
% Computational Cost: add. (247->18), mult. (166->28), div. (48->3), fcn. (26->2), ass. (0->24)
t29 = 2 * pkin(1);
t16 = 1 / pkin(1);
t28 = t16 / 0.2e1;
t25 = (qJ(1) + pkin(3));
t27 = 1 / t25 * t16;
t26 = 1 / pkin(2) * t16;
t24 = pkin(1) + 0;
t23 = -pkin(2) - t25;
t22 = -pkin(2) + t25;
t21 = -t27 / 0.2e1;
t20 = t27 / 0.2e1;
t19 = t26 / 0.2e1;
t18 = pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t4 = sqrt(-((pkin(1) + t23) * (pkin(1) + t22) * (pkin(1) - t22) * (pkin(1) - t23)));
t3 = t4 * t21;
t15 = pkin(1) ^ 2;
t8 = -t15 + t18;
t6 = t8 * t20;
t17 = t15 + t18;
t7 = t17 * t19;
t5 = (0 * t29 + t17) * t28;
t2 = t4 * t20;
t1 = ((0 * t29) - t4) * t28;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t6, t2, 0, t24; t3, t6, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t3, 0, t6, qJ(1) * t6 + t24; t8 * t21, 0, t3, qJ(1) * t3 + 0; 0, -1, 0, 0; 0, 0, 0, 1; t7, t4 * t19, 0, 0; -t4 * t26 / 0.2e1, t7, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t6, t2, 0, t5; t3, t6, 0, t1; 0, 0, 1, 0; 0, 0, 0, 1; t6, t2, 0, t5; t3, t6, 0, t1; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t9;
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
