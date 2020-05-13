% Calculate inertial parameters regressor of fixed base kinetic energy for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% T_reg [1x(12*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = picker2Dm1OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_energykin_fixb_reg2_slag_vp: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_energykin_fixb_reg2_slag_vp: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:45:11
% EndTime: 2020-05-11 05:45:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (164->35), mult. (302->95), div. (0->0), fcn. (136->14), ass. (0->42)
t39 = qJD(1) ^ 2;
t44 = t39 * pkin(1) ^ 2;
t43 = pkin(1) * qJD(1);
t24 = qJD(1) + qJD(2);
t23 = qJD(1) + qJD(8);
t42 = t23 * t43;
t32 = sin(qJ(2));
t41 = t32 * t43;
t38 = cos(qJ(2));
t19 = t38 * t43;
t15 = t24 * pkin(2) + t19;
t31 = sin(qJ(3));
t37 = cos(qJ(3));
t8 = -t37 * t15 + t31 * t41;
t22 = qJD(3) + t24;
t21 = qJD(4) + t24;
t14 = t24 * pkin(3) + t19;
t30 = sin(qJ(4));
t36 = cos(qJ(4));
t7 = t36 * t14 - t30 * t41;
t35 = cos(qJ(6));
t34 = cos(qJ(8));
t33 = cos(qJ(9));
t29 = sin(qJ(6));
t28 = sin(qJ(8));
t27 = sin(qJ(9));
t26 = cos(qJ(10));
t25 = sin(qJ(10));
t20 = qJD(6) + t24;
t18 = qJD(9) + t22;
t17 = qJD(10) + t21;
t12 = (-t29 * t38 - t32 * t35) * t43;
t11 = (t29 * t32 - t35 * t38) * t43;
t10 = -t31 * t15 - t37 * t41;
t9 = t30 * t14 + t36 * t41;
t6 = t22 * pkin(6) + t8;
t5 = t21 * pkin(4) + t7;
t4 = -t33 * t10 - t27 * t6;
t3 = t27 * t10 - t33 * t6;
t2 = -t25 * t5 - t26 * t9;
t1 = t25 * t9 - t26 * t5;
t13 = [0, 0, 0, 0, 0, t39 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 ^ 2 / 0.2e1, t24 * t19, -t24 * t41, 0, (t32 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * t44, 0, 0, 0, 0, 0, t22 ^ 2 / 0.2e1, t8 * t22, -t10 * t22, 0, t10 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t21 ^ 2 / 0.2e1, t7 * t21, -t9 * t21, 0, t9 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, qJD(5) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 ^ 2 / 0.2e1, t11 * t20, -t12 * t20, 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, qJD(7) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 ^ 2 / 0.2e1, -t34 * t42, t28 * t42, 0, (t28 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * t44, 0, 0, 0, 0, 0, t18 ^ 2 / 0.2e1, t3 * t18, -t4 * t18, 0, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t17 ^ 2 / 0.2e1, t1 * t17, -t2 * t17, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t13;
