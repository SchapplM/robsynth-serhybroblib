% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m2DE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:44
% EndTime: 2020-05-03 01:06:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (52->23), mult. (166->67), div. (0->0), fcn. (59->6), ass. (0->29)
t35 = qJD(3) ^ 2 / 0.2e1;
t34 = qJD(2) ^ 2 / 0.2e1;
t22 = qJD(1) ^ 2;
t13 = t22 / 0.2e1;
t18 = cos(qJ(3));
t33 = t18 * t22;
t32 = pkin(4) * qJD(2);
t19 = cos(qJ(2));
t7 = t19 * pkin(4) + pkin(1);
t6 = pkin(2) + t7;
t4 = t18 * pkin(5) + t6;
t2 = pkin(3) + t4;
t31 = qJD(1) * t2;
t30 = pkin(1) * t22;
t15 = sin(qJ(3));
t16 = sin(qJ(2));
t26 = qJD(3) * t32;
t24 = (t15 * t16 + t18 * t19) * t26;
t9 = pkin(4) ^ 2 * t34;
t29 = t9 + (t35 * pkin(5) + t24) * pkin(5);
t28 = qJD(1) * qJD(2);
t27 = qJD(1) * qJD(3);
t25 = t16 * t28;
t23 = pkin(4) * t25;
t17 = cos(qJ(4));
t14 = sin(qJ(4));
t12 = qJD(1) + qJD(4);
t5 = pkin(5) * qJD(3) * t15 + t16 * t32;
t1 = [0, 0, 0, 0, 0, t13, 0, 0, 0, 0, t16 ^ 2 * t13, t16 * t22 * t19, t25, t19 ^ 2 * t13, t19 * t28, t34, t19 * t30, -t16 * t30, 0, pkin(1) ^ 2 * t13, 0, 0, 0, t13, 0, 0, t7 * t22, 0, -t23, t7 ^ 2 * t13 + t9, t15 ^ 2 * t13, t15 * t33, t15 * t27, t18 ^ 2 * t13, t18 * t27, t35, t6 * t33 + t24, -t15 * t6 * t22 + (-t15 * t19 + t18 * t16) * t26, -t23, t6 ^ 2 * t13 + t9, 0, 0, 0, t13, 0, 0, t4 * t22, 0, -t5 * qJD(1), t4 ^ 2 * t13 + t29, 0, 0, 0, 0, 0, t12 ^ 2 / 0.2e1, t12 * (t14 * t5 + t17 * t31), -t12 * (t14 * t31 - t17 * t5), 0, t2 ^ 2 * t13 + t29;];
T_reg = t1;
