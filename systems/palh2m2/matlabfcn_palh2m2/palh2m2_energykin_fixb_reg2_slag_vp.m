% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 14:45:19
% EndTime: 2019-06-06 14:45:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (148->31), mult. (552->95), div. (0->0), fcn. (264->6), ass. (0->38)
t27 = sin(qJ(2));
t22 = t27 ^ 2;
t30 = cos(qJ(2));
t24 = t30 ^ 2;
t17 = (t22 + t24) * qJD(1);
t16 = t17 ^ 2;
t42 = t16 / 0.2e1;
t32 = qJD(1) ^ 2;
t41 = t32 / 0.2e1;
t18 = (-pkin(4) * t30 - pkin(1)) * qJD(1);
t10 = -t17 * pkin(2) + t18;
t40 = t10 * t17;
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t38 = pkin(4) * qJD(2);
t13 = (t26 * t30 - t27 * t29) * t38;
t39 = t29 * t13;
t21 = t26 ^ 2;
t23 = t29 ^ 2;
t9 = (t21 + t23) * t17;
t34 = t27 * t38;
t12 = t29 * t30 * t38 + t26 * t34;
t37 = pkin(1) * t32;
t36 = qJD(3) * t17;
t35 = qJD(1) * qJD(2);
t7 = (-pkin(5) * t29 - pkin(2)) * t17 + t18;
t31 = qJD(2) ^ 2;
t28 = cos(qJ(4));
t25 = sin(qJ(4));
t11 = qJD(3) * pkin(5) + t12;
t8 = qJD(4) + t9;
t6 = -t26 * t11 + t39;
t5 = t29 * t11 + t26 * t13;
t4 = t5 ^ 2 / 0.2e1;
t3 = -t9 * pkin(3) + t7;
t2 = -t25 * t3 + t28 * t6;
t1 = -t25 * t6 - t28 * t3;
t14 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, t22 * t41, t27 * t32 * t30, t27 * t35, t24 * t41, t30 * t35, t31 / 0.2e1, t30 * t37, -t27 * t37, 0, pkin(1) ^ 2 * t41, 0, 0, 0, t42, 0, 0, -t18 * t17, 0, -t17 * t34, t18 ^ 2 / 0.2e1 + (t22 / 0.2e1 + t24 / 0.2e1) * pkin(4) ^ 2 * t31, t21 * t42, t26 * t16 * t29, t26 * t36, t23 * t42, t29 * t36, qJD(3) ^ 2 / 0.2e1, t12 * qJD(3) - t29 * t40, -t13 * qJD(3) + t26 * t40, (-t12 * t26 + t39) * t17, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, 0, 0, 0, t9 ^ 2 / 0.2e1, 0, 0, -t7 * t9, 0, t6 * t9, t6 ^ 2 / 0.2e1 + t4 + t7 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t8 ^ 2 / 0.2e1, t1 * t8, -t2 * t8, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4;];
T_reg  = t14;
