% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% palh2m2DE
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
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m2DE_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:46
% EndTime: 2020-05-03 01:06:48
% DurationCPUTime: 0.35s
% Computational Cost: add. (71->24), mult. (219->69), div. (0->0), fcn. (133->6), ass. (0->34)
t38 = pkin(4) * (qJD(2) - qJD(3));
t16 = sin(qJ(2));
t35 = qJD(2) * t16;
t28 = pkin(4) * t35;
t15 = sin(qJ(3));
t33 = qJD(3) * t15;
t10 = -pkin(5) * t33 - t28;
t37 = 0.2e1 * t10;
t18 = cos(qJ(3));
t19 = cos(qJ(2));
t36 = (t15 * t19 - t16 * t18) * t38;
t34 = qJD(2) * t19;
t32 = qJD(3) * t18;
t14 = sin(qJ(4));
t31 = qJD(4) * t14;
t17 = cos(qJ(4));
t30 = qJD(4) * t17;
t29 = t19 * pkin(4) + pkin(1);
t27 = pkin(4) * t34;
t26 = t15 * t32;
t25 = t16 * t34;
t24 = -0.2e1 * pkin(1) * qJD(2);
t11 = pkin(2) + t29;
t23 = -0.2e1 * t28;
t22 = t18 * pkin(5) + t11;
t9 = pkin(3) + t22;
t8 = (t14 * t34 + t16 * t30) * pkin(4);
t7 = (t14 * t32 + t15 * t30) * pkin(5);
t6 = (-t16 * t31 + t17 * t34) * pkin(4);
t5 = (-t15 * t31 + t17 * t32) * pkin(5);
t3 = pkin(5) * t36;
t2 = -t14 * t10 - t30 * t9;
t1 = t17 * t10 - t31 * t9;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t25, 0.2e1 * (-t16 ^ 2 + t19 ^ 2) * qJD(2), 0, -0.2e1 * t25, 0, 0, t16 * t24, t19 * t24, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, t29 * t23, 0.2e1 * t26, 0.2e1 * (-t15 ^ 2 + t18 ^ 2) * qJD(3), 0, -0.2e1 * t26, 0, 0, -0.2e1 * t11 * t33 - 0.2e1 * t18 * t28, -0.2e1 * t11 * t32 + 0.2e1 * t15 * t28, 0, t11 * t23, 0, 0, 0, 0, 0, 0, t37, 0, 0, t22 * t37, 0, 0, 0, 0, 0, 0, 0.2e1 * t1, 0.2e1 * t2, 0, t9 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, 0, 0, 0, 0, t8, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(5) * t32, 0, 0, 0, 0, 0, 0, 0, t7, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, (t15 * t16 + t18 * t19) * t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
