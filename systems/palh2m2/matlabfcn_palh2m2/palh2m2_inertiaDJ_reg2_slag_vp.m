% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% palh2m2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 14:45:19
% EndTime: 2019-06-06 14:45:21
% DurationCPUTime: 0.41s
% Computational Cost: add. (413->42), mult. (1366->126), div. (0->0), fcn. (957->6), ass. (0->58)
t24 = sin(qJ(3));
t25 = sin(qJ(2));
t27 = cos(qJ(3));
t28 = cos(qJ(2));
t13 = (t24 * t25 + t27 * t28) * pkin(4);
t14 = (t24 * t28 - t25 * t27) * pkin(4);
t62 = qJD(2) - qJD(3);
t19 = t25 ^ 2;
t22 = t28 ^ 2;
t16 = t19 + t22;
t18 = t24 ^ 2;
t21 = t27 ^ 2;
t10 = (t18 + t21) * t16;
t61 = 0.2e1 * t10;
t52 = qJD(3) * t24;
t44 = t16 * t52;
t48 = t25 * qJD(2);
t46 = pkin(4) * t48;
t11 = pkin(5) * t44 + t46;
t60 = 0.2e1 * t11;
t59 = 0.2e1 * t16;
t6 = t27 * t13 + t24 * t14;
t8 = t62 * t13;
t9 = t62 * t14;
t3 = -qJD(3) * t6 - t24 * t9 - t27 * t8;
t7 = -t24 * t13 + t27 * t14;
t58 = t7 * t3;
t57 = pkin(5) * t27;
t56 = t24 * t3;
t55 = t24 * t6;
t54 = t27 * t7;
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t53 = t23 ^ 2 + t26 ^ 2;
t51 = qJD(3) * t27;
t50 = qJD(4) * t23;
t49 = qJD(4) * t26;
t47 = t28 * qJD(2);
t45 = t24 * t51;
t43 = t25 * t47;
t42 = -0.2e1 * pkin(1) * qJD(2);
t41 = 0.2e1 * t46;
t40 = -t28 * pkin(4) - pkin(1);
t15 = t16 ^ 2;
t39 = t15 * t45;
t36 = -t23 * t3 - t7 * t49;
t35 = -t26 * t3 + t7 * t50;
t30 = (-pkin(2) - t57) * t16 + t40;
t5 = -t10 * pkin(3) + t30;
t34 = t23 * t11 + t5 * t49;
t33 = -t26 * t11 + t5 * t50;
t32 = (-t24 * t50 + t26 * t51) * pkin(5);
t31 = (t23 * t51 + t24 * t49) * pkin(5);
t12 = -t16 * pkin(2) + t40;
t4 = t7 * qJD(3) - t24 * t8 + t27 * t9;
t2 = t4 * t57;
t1 = t6 * t4;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, 0.2e1 * (-t19 + t22) * qJD(2), 0, -0.2e1 * t43, 0, 0, t25 * t42, t28 * t42, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t16 * t46, 0, 0, t40 * t41, 0.2e1 * t39, 0.2e1 * (-t18 + t21) * t15 * qJD(3), 0, -0.2e1 * t39, 0, 0, (t12 * t52 - t27 * t46) * t59, (t12 * t51 + t24 * t46) * t59, 0, t12 * t41, 0, 0, 0, 0, 0, 0, -0.2e1 * t11 * t10, 0, 0, t30 * t60, 0, 0, 0, 0, 0, 0, t33 * t61, t34 * t61, 0, t53 * t5 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(4) * t16 * t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t10, 0, 0, 0, 0, 0, 0, 0, t36 * t10, t35 * t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t13 * t9 - 0.2e1 * t14 * t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t1 + 0.2e1 * t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t53 * t58 + 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t51, 0, -t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(5) * t10 * t51, 0, 0, 0, 0, 0, 0, 0, t10 * t31, t10 * t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 + (-t56 + (-t54 - t55) * qJD(3)) * pkin(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 + (-t53 * t56 + (-t53 * t54 - t55) * qJD(3)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t53) * pkin(5) ^ 2 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t17;
