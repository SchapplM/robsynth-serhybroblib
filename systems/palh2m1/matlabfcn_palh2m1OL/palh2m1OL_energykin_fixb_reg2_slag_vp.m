% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m1OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:27:43
% EndTime: 2020-05-03 00:27:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (309->46), mult. (729->126), div. (0->0), fcn. (570->8), ass. (0->40)
t36 = qJD(1) ^ 2;
t48 = t36 / 0.2e1;
t30 = sin(qJ(3));
t47 = pkin(2) * t30;
t46 = cos(qJ(5));
t31 = sin(qJ(2));
t45 = t31 * t30;
t34 = cos(qJ(2));
t44 = t34 * t36;
t33 = cos(qJ(3));
t23 = t33 * pkin(2) + pkin(3);
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t42 = pkin(3) * qJD(3);
t13 = (t23 * t29 + t32 * t47) * qJD(2) + t29 * t42;
t43 = pkin(2) * qJD(2);
t24 = t34 * pkin(2) + pkin(1);
t41 = qJD(1) * t24;
t17 = (t33 * pkin(3) + pkin(2)) * t34 - pkin(3) * t45 + pkin(1);
t40 = t17 * qJD(1);
t39 = qJD(1) * qJD(2);
t27 = qJD(2) + qJD(3);
t38 = t27 * t43;
t19 = (-t33 * t34 + t45) * qJD(1);
t20 = (-t30 * t34 - t31 * t33) * qJD(1);
t8 = -t32 * t19 + t29 * t20;
t14 = (t23 * t32 - t29 * t47) * qJD(2) + t32 * t42;
t35 = qJD(2) ^ 2;
t28 = sin(qJ(5));
t26 = qJD(4) + t27;
t12 = -t26 * pkin(4) - t14;
t11 = t26 * pkin(6) + t13;
t10 = t29 * t19 + t32 * t20;
t7 = qJD(5) + t8;
t6 = t46 * t10 + t28 * t26;
t4 = t28 * t10 - t46 * t26;
t3 = t8 * pkin(4) - t10 * pkin(6) + t40;
t2 = t46 * t11 + t28 * t3;
t1 = -t28 * t11 + t46 * t3;
t5 = [0, 0, 0, 0, 0, t48, 0, 0, 0, 0, t31 ^ 2 * t48, t31 * t44, -t31 * t39, t34 ^ 2 * t48, -t34 * t39, t35 / 0.2e1, pkin(1) * t44, -t36 * pkin(1) * t31, 0, pkin(1) ^ 2 * t48, t20 ^ 2 / 0.2e1, t20 * t19, t20 * t27, t19 ^ 2 / 0.2e1, t19 * t27, t27 ^ 2 / 0.2e1, -t19 * t41 + t33 * t38, t20 * t41 - t30 * t38, (t19 * t30 - t20 * t33) * t43, t24 ^ 2 * t48 + (t30 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t35, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t26, t8 ^ 2 / 0.2e1, -t8 * t26, t26 ^ 2 / 0.2e1, t14 * t26 + t8 * t40, t10 * t40 - t13 * t26, -t14 * t10 - t13 * t8, t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t17 ^ 2 * t48, t6 ^ 2 / 0.2e1, -t6 * t4, t6 * t7, t4 ^ 2 / 0.2e1, -t4 * t7, t7 ^ 2 / 0.2e1, t1 * t7 + t12 * t4, t12 * t6 - t2 * t7, -t1 * t6 - t2 * t4, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t5;
