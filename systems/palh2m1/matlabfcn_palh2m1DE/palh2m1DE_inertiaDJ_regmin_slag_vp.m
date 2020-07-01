% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m1DE_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:35
% EndTime: 2020-06-30 17:39:35
% DurationCPUTime: 0.20s
% Computational Cost: add. (129->33), mult. (373->78), div. (0->0), fcn. (278->6), ass. (0->42)
t46 = qJD(2) + qJD(3);
t45 = 2 * qJD(2);
t24 = sin(qJ(3));
t25 = sin(qJ(2));
t44 = t25 * t24;
t28 = cos(qJ(2));
t43 = t28 * t24;
t27 = cos(qJ(3));
t42 = t28 * t27;
t21 = t27 ^ 2;
t41 = -t24 ^ 2 + t21;
t22 = t28 ^ 2;
t40 = -t25 ^ 2 + t22;
t39 = qJD(2) * t25;
t38 = qJD(2) * t28;
t37 = qJD(3) * t27;
t23 = sin(qJ(4));
t36 = qJD(4) * t23;
t26 = cos(qJ(4));
t35 = qJD(4) * t26;
t34 = -2 * pkin(1) * qJD(2);
t33 = pkin(2) * t39;
t32 = qJD(3) * t24 * pkin(2);
t31 = pkin(2) * t37;
t30 = -0.4e1 * t42 * t44;
t13 = t27 * t25 + t43;
t29 = -t42 + t44;
t17 = t28 * pkin(2) + pkin(1);
t16 = t27 * pkin(3) + pkin(2);
t12 = pkin(3) * t43 + t25 * t16;
t11 = -pkin(3) * t44 + t16 * t28 + pkin(1) + pkin(4);
t10 = t46 * t13;
t9 = -t27 * t38 - t28 * t37 + t46 * t44;
t8 = t16 * t38 + (-qJD(3) * t29 - t24 * t39) * pkin(3);
t7 = -t16 * t39 + (-qJD(3) * t13 - t24 * t38) * pkin(3);
t6 = (t13 * t35 - t23 * t9) * pkin(3);
t5 = (-t13 * t36 - t26 * t9) * pkin(3);
t4 = t12 * t35 + t23 * t8;
t3 = -t12 * t36 + t26 * t8;
t2 = -t11 * t35 - t23 * t7;
t1 = -t11 * t36 + t26 * t7;
t14 = [0, 0, 0, 0.2e1 * t25 * t38, t40 * t45, 0, 0, 0, t25 * t34, t28 * t34, -0.2e1 * t13 * t9, (t30 + 0.2e1 * t40 * (t21 - 0.1e1 / 0.2e1)) * t45 + 0.2e1 * (0.2e1 * t22 * t41 + t30 - t41) * qJD(3), 0, 0, 0, -0.2e1 * t17 * t10 + 0.2e1 * t29 * t33, 0.2e1 * t13 * t33 + 0.2e1 * t17 * t9, 0.2e1 * t7, 0, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 0, 0, -t38, t39, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t32, -0.2e1 * t31, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
