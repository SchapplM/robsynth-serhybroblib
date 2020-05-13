% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m1DE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:29
% EndTime: 2020-05-02 23:52:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (90->44), mult. (249->99), div. (0->0), fcn. (96->11), ass. (0->43)
t53 = qJD(2) ^ 2 / 0.2e1;
t29 = qJD(1) ^ 2;
t18 = t29 / 0.2e1;
t52 = t29 / 0.4e1;
t51 = pkin(1) * t29;
t25 = cos(qJ(3));
t50 = t25 * pkin(3);
t26 = cos(qJ(2));
t49 = t26 * pkin(2);
t12 = pkin(1) + t49;
t48 = t12 * t29;
t22 = sin(qJ(3));
t23 = sin(qJ(2));
t47 = t23 * t22;
t46 = t26 * t22;
t45 = t26 * t29;
t44 = pkin(2) * qJD(2);
t43 = pkin(3) * qJD(3);
t11 = pkin(2) + t50;
t35 = -pkin(3) * t47 + t11 * t26 + pkin(1);
t42 = qJD(1) * (pkin(4) + t35);
t16 = qJD(2) + qJD(3);
t41 = qJD(1) * t16;
t40 = qJD(1) * qJD(2);
t31 = pkin(3) ^ 2;
t32 = pkin(2) ^ 2;
t39 = 0.2e1 * pkin(2) * t50 + t31 + t32;
t38 = t16 * t44;
t37 = t23 * t40;
t20 = qJ(2) + qJ(3);
t36 = pkin(3) * cos(t20) + t49;
t30 = 0.2e1 * qJ(2);
t34 = t39 * t53 + cos(t30) * t32 * t52 + (t25 * pkin(2) + pkin(3)) * qJD(2) * t43 + pkin(2) * pkin(3) * cos(qJ(3) + t30) * t18 + (qJD(3) ^ 2 / 0.2e1 + cos(0.2e1 * t20) * t52) * t31;
t33 = pkin(1) ^ 2;
t27 = pkin(1) + pkin(4);
t24 = cos(qJ(4));
t21 = sin(qJ(4));
t19 = t26 ^ 2;
t17 = qJD(1) + qJD(4);
t5 = t26 * t25 - t47;
t4 = t25 * t23 + t46;
t1 = (pkin(3) * t46 + t23 * t11) * qJD(2) + t4 * t43;
t2 = [0, 0, 0, 0, 0, t18, 0, 0, 0, 0, t23 ^ 2 * t18, t23 * t45, -t37, t19 * t18, -t26 * t40, t53, pkin(1) * t45, -t23 * t51, 0, t33 * t18, t4 ^ 2 * t18, 0.2e1 * (t23 * (t25 ^ 2 - 0.1e1 / 0.2e1) * t26 + (t19 - 0.1e1 / 0.2e1) * t25 * t22) * t29, -t4 * t41, t5 ^ 2 * t18, -t5 * t41, t16 ^ 2 / 0.2e1, t25 * t38 + t5 * t48, -t22 * t38 - t4 * t48, pkin(2) * t37, t12 ^ 2 * t18 + t32 * t53, 0, 0, 0, t18, 0, 0, t35 * t29, 0, (pkin(3) * t16 * sin(t20) + t23 * t44) * qJD(1), (0.2e1 * t33 + t39) * t52 + t36 * t51 + t34, 0, 0, 0, 0, 0, t17 ^ 2 / 0.2e1, (t21 * t1 + t24 * t42) * t17, -(-t1 * t24 + t21 * t42) * t17, 0, t39 * t52 + (t27 / 0.2e1 + t36) * t29 * t27 + t34;];
T_reg = t2;
