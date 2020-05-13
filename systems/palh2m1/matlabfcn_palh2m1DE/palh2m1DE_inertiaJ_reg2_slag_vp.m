% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m1DE_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t53 = 2 * pkin(1);
t35 = cos(qJ(2));
t50 = pkin(2) * t35;
t21 = pkin(1) + t50;
t52 = 0.2e1 * t21;
t32 = sin(qJ(2));
t34 = cos(qJ(3));
t31 = sin(qJ(3));
t43 = t35 * t31;
t7 = -t34 * t32 - t43;
t51 = pkin(3) * t7;
t38 = pkin(3) ^ 2;
t28 = t38 / 0.2e1;
t39 = pkin(2) ^ 2;
t29 = t39 / 0.2e1;
t27 = qJ(2) + qJ(3);
t49 = pkin(3) * cos(t27);
t20 = pkin(3) * t34 + pkin(2);
t45 = t32 * t31;
t42 = -pkin(3) * t45 + t20 * t35 + pkin(1);
t3 = pkin(4) + t42;
t33 = cos(qJ(4));
t48 = t3 * t33;
t30 = sin(qJ(4));
t47 = t30 * t3;
t46 = t31 * pkin(2);
t25 = t34 * pkin(2);
t44 = t32 * t35;
t22 = pkin(3) * t25;
t37 = 0.2e1 * qJ(2);
t41 = cos(0.2e1 * t27) * t28 + pkin(2) * pkin(3) * cos(qJ(3) + t37) + cos(t37) * t29 + t22 + t28 + t29;
t40 = pkin(1) ^ 2;
t36 = pkin(1) + pkin(4);
t26 = t35 ^ 2;
t24 = t32 * pkin(2);
t16 = pkin(3) * sin(t27);
t13 = pkin(3) * (t25 + pkin(3));
t12 = 0.2e1 * t22 + t38 + t39;
t9 = t34 * t35 - t45;
t6 = pkin(3) * t43 + t20 * t32;
t5 = t33 * t51;
t4 = t30 * t51;
t2 = t6 * t33;
t1 = t30 * t6;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32 ^ 2, 0.2e1 * t44, 0, t26, 0, 0, t35 * t53, -0.2e1 * pkin(1) * t32, 0, t40, t7 ^ 2, 0.4e1 * (t34 ^ 2 - 0.1e1 / 0.2e1) * t44 + (0.4e1 * t26 - 0.2e1) * t34 * t31, 0, t9 ^ 2, 0, 0, t9 * t52, t7 * t52, 0, t21 ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t42, 0, 0, t40 + (t49 + t50) * t53 + t41, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t47, 0, (t36 + 0.2e1 * t49 + 0.2e1 * t50) * t36 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t35, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t9, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 + t24, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t46, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t48, -t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
