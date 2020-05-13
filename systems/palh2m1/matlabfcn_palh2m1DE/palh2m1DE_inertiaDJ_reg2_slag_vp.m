% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m1DE_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:32
% EndTime: 2020-05-02 23:52:34
% DurationCPUTime: 0.34s
% Computational Cost: add. (165->47), mult. (471->107), div. (0->0), fcn. (314->11), ass. (0->56)
t27 = qJD(2) + qJD(3);
t66 = 2 * qJD(2);
t40 = pkin(1) + pkin(4);
t65 = pkin(3) * t27;
t35 = sin(qJ(3));
t36 = sin(qJ(2));
t64 = t36 * t35;
t39 = cos(qJ(2));
t63 = t39 * t35;
t38 = cos(qJ(3));
t62 = t39 * t38;
t31 = t38 ^ 2;
t61 = -t35 ^ 2 + t31;
t32 = t39 ^ 2;
t60 = -t36 ^ 2 + t32;
t59 = qJD(2) * t36;
t58 = qJD(2) * t39;
t57 = qJD(3) * t35;
t56 = qJD(3) * t38;
t34 = sin(qJ(4));
t55 = qJD(4) * t34;
t37 = cos(qJ(4));
t54 = qJD(4) * t37;
t53 = -0.2e1 * t59;
t52 = pkin(2) * t59;
t51 = pkin(2) * t57;
t50 = pkin(2) * t56;
t49 = t36 * t58;
t33 = qJ(2) + qJ(3);
t48 = -0.2e1 * sin(t33) * t65;
t47 = pkin(1) * t53;
t46 = pkin(3) * t51;
t45 = -0.4e1 * t62 * t64;
t41 = 0.2e1 * qJ(2);
t44 = -qJD(2) * sin(t41) * pkin(2) ^ 2 - pkin(3) ^ 2 * t27 * sin(0.2e1 * t33);
t14 = t38 * t36 + t63;
t15 = t62 - t64;
t43 = (-sin(qJ(3) + t41) * (qJD(3) + t66) - t57) * pkin(3);
t22 = pkin(2) * t58;
t21 = t39 * pkin(2) + pkin(1);
t20 = t38 * pkin(3) + pkin(2);
t19 = -0.2e1 * t46;
t16 = cos(t33) * t65;
t12 = pkin(3) * t63 + t36 * t20;
t11 = -pkin(3) * t64 + t20 * t39 + t40;
t10 = t27 * t14;
t9 = t27 * t64 - t38 * t58 - t39 * t56;
t8 = t20 * t58 + (t15 * qJD(3) - t35 * t59) * pkin(3);
t7 = -t20 * t59 + (-t14 * qJD(3) - t35 * t58) * pkin(3);
t6 = (-t14 * t55 - t37 * t9) * pkin(3);
t5 = (t14 * t54 - t34 * t9) * pkin(3);
t4 = -t12 * t55 + t8 * t37;
t3 = t12 * t54 + t34 * t8;
t2 = -t11 * t55 + t7 * t37;
t1 = -t11 * t54 - t34 * t7;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t49, t60 * t66, 0, -0.2e1 * t49, 0, 0, t47, -0.2e1 * pkin(1) * t58, 0, 0, -0.2e1 * t14 * t9, (t45 + 0.2e1 * t60 * (t31 - 0.1e1 / 0.2e1)) * t66 + 0.2e1 * (0.2e1 * t32 * t61 + t45 - t61) * qJD(3), 0, -0.2e1 * t15 * t10, 0, 0, -0.2e1 * t21 * t10 - 0.2e1 * t15 * t52, 0.2e1 * t14 * t52 + 0.2e1 * t21 * t9, 0, -0.2e1 * t21 * t52, 0, 0, 0, 0, 0, 0, 0.2e1 * t7, 0, 0, pkin(1) * t48 + (t43 + t47) * pkin(2) + t44, 0, 0, 0, 0, 0, 0, 0.2e1 * t2, 0.2e1 * t1, 0, t40 * t48 + (t40 * t53 + t43) * pkin(2) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, t59, 0, 0, 0, 0, 0, 0, 0, t9, 0, t10, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 + t22, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t51, -0.2e1 * t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
