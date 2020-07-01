% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m1DE_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:34
% EndTime: 2020-06-30 17:39:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (47->24), mult. (108->58), div. (0->0), fcn. (71->6), ass. (0->24)
t25 = qJD(1) ^ 2;
t35 = t25 / 0.2e1;
t24 = cos(qJ(2));
t34 = (t24 * pkin(2) + pkin(1)) * t25;
t20 = sin(qJ(3));
t21 = sin(qJ(2));
t33 = t21 * t20;
t32 = t24 * t20;
t31 = t24 * t25;
t23 = cos(qJ(3));
t15 = t23 * pkin(3) + pkin(2);
t26 = -pkin(3) * t33 + t15 * t24 + pkin(1);
t30 = qJD(1) * (pkin(4) + t26);
t17 = qJD(2) + qJD(3);
t29 = qJD(1) * t17;
t28 = qJD(1) * qJD(2);
t27 = pkin(2) * qJD(2) * t17;
t22 = cos(qJ(4));
t19 = sin(qJ(4));
t18 = qJD(1) + qJD(4);
t13 = t24 * t23 - t33;
t12 = t23 * t21 + t32;
t10 = (pkin(3) * t32 + t21 * t15) * qJD(2) + qJD(3) * pkin(3) * t12;
t1 = [t35, 0, 0, t21 ^ 2 * t35, t21 * t31, -t21 * t28, -t24 * t28, qJD(2) ^ 2 / 0.2e1, pkin(1) * t31, -t25 * pkin(1) * t21, t12 ^ 2 * t35, 0.2e1 * (t21 * (t23 ^ 2 - 0.1e1 / 0.2e1) * t24 + (t24 ^ 2 - 0.1e1 / 0.2e1) * t23 * t20) * t25, -t12 * t29, -t13 * t29, t17 ^ 2 / 0.2e1, t13 * t34 + t23 * t27, -t12 * t34 - t20 * t27, t26 * t25, t18 ^ 2 / 0.2e1, (t19 * t10 + t22 * t30) * t18, -(-t22 * t10 + t19 * t30) * t18;];
T_reg = t1;
