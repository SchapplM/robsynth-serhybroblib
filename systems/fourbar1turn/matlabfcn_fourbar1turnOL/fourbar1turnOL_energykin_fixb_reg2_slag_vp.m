% Calculate inertial parameters regressor of fixed base kinetic energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1turnOL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:04
% EndTime: 2020-04-12 19:41:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (21->9), mult. (131->62), div. (0->0), fcn. (63->6), ass. (0->19)
t13 = qJD(1) ^ 2;
t21 = t13 / 0.2e1;
t9 = cos(qJ(4));
t20 = t13 * t9;
t4 = qJD(2) + qJD(3);
t19 = qJD(2) * t4;
t11 = cos(qJ(2));
t18 = qJD(1) * t11;
t17 = qJD(1) * qJD(2);
t16 = qJD(1) * qJD(4);
t15 = t11 ^ 2 * t21;
t12 = qJD(2) ^ 2;
t10 = cos(qJ(3));
t8 = sin(qJ(2));
t7 = sin(qJ(3));
t6 = sin(qJ(4));
t2 = (-t10 * t8 - t11 * t7) * qJD(1);
t1 = (-t10 * t11 + t7 * t8) * qJD(1);
t3 = [0, 0, 0, 0, 0, t21, 0, 0, 0, 0, t8 ^ 2 * t21, t8 * t13 * t11, t8 * t17, t15, t11 * t17, t12 / 0.2e1, 0, 0, 0, 0, t2 ^ 2 / 0.2e1, t2 * t1, t2 * t4, t1 ^ 2 / 0.2e1, t1 * t4, t4 ^ 2 / 0.2e1, (t1 * t18 - t10 * t19) * pkin(2), (-t2 * t18 + t7 * t19) * pkin(2), (-t1 * t7 + t10 * t2) * qJD(2) * pkin(2), (t15 + (t7 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1) * t12) * pkin(2) ^ 2, t6 ^ 2 * t21, t6 * t20, t6 * t16, t9 ^ 2 * t21, t9 * t16, qJD(4) ^ 2 / 0.2e1, pkin(1) * t20, -t13 * pkin(1) * t6, 0, pkin(1) ^ 2 * t21;];
T_reg = t3;
