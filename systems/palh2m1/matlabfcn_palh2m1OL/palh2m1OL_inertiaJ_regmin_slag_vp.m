% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m1OL_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:12
% EndTime: 2020-06-30 17:45:14
% DurationCPUTime: 0.40s
% Computational Cost: add. (300->44), mult. (643->98), div. (0->0), fcn. (838->8), ass. (0->51)
t37 = cos(qJ(2));
t23 = t37 * pkin(2) + pkin(1);
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t36 = cos(qJ(3));
t38 = t32 * t33 - t36 * t37;
t10 = -t38 * pkin(3) + t23;
t54 = 0.2e1 * t10;
t16 = -t32 * t37 - t36 * t33;
t53 = 0.2e1 * t16;
t52 = 0.2e1 * t37;
t30 = sin(qJ(5));
t51 = pkin(4) * t30;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t7 = t31 * t16 - t35 * t38;
t4 = t30 * t7;
t50 = t31 * pkin(3);
t49 = t32 * pkin(2);
t34 = cos(qJ(5));
t5 = t34 * t7;
t25 = t36 * pkin(2);
t22 = t25 + pkin(3);
t42 = -t35 * t22 + t31 * t49;
t11 = -pkin(4) + t42;
t48 = t11 * t34;
t24 = t35 * pkin(3);
t21 = -t24 - pkin(4);
t47 = t21 * t34;
t46 = t30 * t34;
t45 = -0.2e1 * t4;
t44 = 0.2e1 * t5;
t43 = t35 * t49;
t8 = t35 * t16 + t31 * t38;
t41 = -pkin(4) * t8 - pkin(6) * t7;
t14 = -t31 * t22 - t43;
t12 = pkin(6) - t14;
t40 = t11 * t8 - t12 * t7;
t20 = pkin(6) + t50;
t39 = -t20 * t7 + t21 * t8;
t29 = t34 ^ 2;
t28 = t30 ^ 2;
t27 = pkin(4) * t34;
t19 = 0.2e1 * t46;
t18 = t21 * t30;
t9 = t11 * t30;
t6 = t8 ^ 2;
t3 = t8 * t46;
t2 = (-t28 + t29) * t8;
t1 = t7 * pkin(4) - t8 * pkin(6) + t10;
t13 = [1, 0, 0, t33 ^ 2, t33 * t52, 0, 0, 0, pkin(1) * t52, -0.2e1 * pkin(1) * t33, t16 ^ 2, t38 * t53, 0, 0, 0, -0.2e1 * t23 * t38, t23 * t53, t6, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t54, t8 * t54, t29 * t6, -0.2e1 * t6 * t46, t8 * t44, t8 * t45, t7 ^ 2, t1 * t44, t1 * t45; 0, 0, 0, 0, 0, -t33, -t37, 0, 0, 0, 0, 0, t16, t38, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, t3, t2, t4, t5, 0, t40 * t30, t40 * t34; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t49, 0, 0, 0, 0, 1, -0.2e1 * t42, 0.2e1 * t14, t28, t19, 0, 0, 0, -0.2e1 * t48, 0.2e1 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t38, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, t3, t2, t4, t5, 0, t39 * t30, t39 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t49, 0, 0, 0, 0, 1, t24 - t42, -t43 + (-pkin(3) - t22) * t31, t28, t19, 0, 0, 0, (-t11 - t21) * t34, t18 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, -0.2e1 * t50, t28, t19, 0, 0, 0, -0.2e1 * t47, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, t3, t2, t4, t5, 0, t41 * t30, t41 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t42, t14, t28, t19, 0, 0, 0, t27 - t48, t9 - t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, -t50, t28, t19, 0, 0, 0, t27 - t47, t18 - t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, t19, 0, 0, 0, 0.2e1 * t27, -0.2e1 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t8, -t30 * t8, t7, t34 * t1, -t30 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34, 0, -t30 * t12, -t34 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34, 0, -t30 * t20, -t34 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34, 0, -t30 * pkin(6), -t34 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t13;
