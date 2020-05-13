% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% T_reg [1x(6*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m2OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 02:56:58
% EndTime: 2020-05-03 02:56:59
% DurationCPUTime: 0.27s
% Computational Cost: add. (620->68), mult. (1429->176), div. (0->0), fcn. (1228->10), ass. (0->54)
t49 = qJD(1) ^ 2;
t65 = t49 / 0.2e1;
t64 = cos(qJ(5));
t63 = cos(qJ(6));
t42 = sin(qJ(4));
t43 = sin(qJ(3));
t62 = t42 * t43;
t44 = sin(qJ(2));
t61 = t44 * t43;
t60 = pkin(4) * qJD(2);
t59 = pkin(5) * qJD(4);
t58 = pkin(1) * t49;
t57 = pkin(2) * t42;
t46 = cos(qJ(3));
t47 = cos(qJ(2));
t23 = (-pkin(2) * t46 - pkin(4)) * t47 + pkin(2) * t61 - pkin(1);
t56 = qJD(1) * t23;
t45 = cos(qJ(4));
t34 = pkin(5) * t45 + pkin(2);
t16 = (pkin(5) * t62 - t34 * t46 - pkin(4)) * t47 + (pkin(5) * t42 * t46 + t34 * t43) * t44 - pkin(1);
t55 = t16 * qJD(1);
t54 = pkin(2) * qJD(3);
t53 = qJD(1) * qJD(2);
t39 = qJD(2) + qJD(3);
t36 = pkin(4) * t46 + pkin(2);
t51 = -pkin(4) * t62 + t36 * t45;
t24 = pkin(5) + t51;
t25 = pkin(4) * t43 * t45 + t36 * t42;
t35 = pkin(2) * t45 + pkin(5);
t41 = sin(qJ(5));
t14 = (t24 * t64 - t25 * t41) * qJD(2) + (t35 * t64 - t41 * t57) * qJD(3) + t64 * t59;
t52 = t39 * t60;
t26 = (-t46 * t47 + t61) * qJD(1);
t28 = (t43 * t47 + t44 * t46) * qJD(1);
t18 = -t26 * t45 - t28 * t42;
t19 = -t26 * t42 + t28 * t45;
t9 = -t64 * t18 + t19 * t41;
t38 = qJD(4) + t39;
t48 = qJD(2) ^ 2;
t40 = sin(qJ(6));
t33 = qJD(5) + t38;
t29 = (-pkin(4) * t47 - pkin(1)) * qJD(1);
t21 = qJD(2) * t51 + t45 * t54;
t20 = qJD(2) * t25 + t42 * t54;
t13 = (t24 * t41 + t25 * t64) * qJD(2) + (t35 * t41 + t57 * t64) * qJD(3) + t41 * t59;
t12 = pkin(3) * t33 + t14;
t11 = t41 * t18 + t19 * t64;
t7 = -qJD(6) + t9;
t6 = t11 * t63 - t40 * t33;
t4 = t40 * t11 + t33 * t63;
t3 = pkin(3) * t9 + t55;
t2 = t13 * t63 - t40 * t3;
t1 = -t40 * t13 - t3 * t63;
t5 = [0, 0, 0, 0, 0, t65, 0, 0, 0, 0, t44 ^ 2 * t65, t44 * t49 * t47, t44 * t53, t47 ^ 2 * t65, t47 * t53, t48 / 0.2e1, t47 * t58, -t44 * t58, 0, pkin(1) ^ 2 * t65, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t39, t26 ^ 2 / 0.2e1, -t26 * t39, t39 ^ 2 / 0.2e1, t26 * t29 + t46 * t52, t28 * t29 - t43 * t52, (-t26 * t43 - t28 * t46) * t60, t29 ^ 2 / 0.2e1 + (t43 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * pkin(4) ^ 2 * t48, t19 ^ 2 / 0.2e1, t19 * t18, t19 * t38, t18 ^ 2 / 0.2e1, t18 * t38, t38 ^ 2 / 0.2e1, -t18 * t56 + t21 * t38, t19 * t56 - t20 * t38, t18 * t20 - t19 * t21, t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t23 ^ 2 * t65, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t33, t9 ^ 2 / 0.2e1, -t9 * t33, t33 ^ 2 / 0.2e1, t14 * t33 + t55 * t9, t11 * t55 - t13 * t33, -t11 * t14 - t13 * t9, t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t16 ^ 2 * t65, t6 ^ 2 / 0.2e1, -t6 * t4, -t6 * t7, t4 ^ 2 / 0.2e1, t4 * t7, t7 ^ 2 / 0.2e1, -t1 * t7 + t12 * t4, t12 * t6 + t2 * t7, -t1 * t6 - t2 * t4, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t5;
