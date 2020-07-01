% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m2DE_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:45
% EndTime: 2020-06-30 17:56:45
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->14), mult. (81->50), div. (0->0), fcn. (39->6), ass. (0->21)
t24 = qJD(1) ^ 2;
t34 = t24 / 0.2e1;
t22 = cos(qJ(3));
t33 = t22 * t24;
t32 = pkin(4) * qJD(2);
t31 = pkin(1) * t24;
t23 = cos(qJ(2));
t29 = t23 * pkin(4) + pkin(1);
t14 = pkin(2) + t29;
t25 = t22 * pkin(5) + t14;
t30 = qJD(1) * (pkin(3) + t25);
t28 = qJD(1) * qJD(2);
t27 = qJD(1) * qJD(3);
t26 = qJD(3) * t32;
t21 = cos(qJ(4));
t20 = sin(qJ(2));
t19 = sin(qJ(3));
t18 = sin(qJ(4));
t17 = qJD(1) + qJD(4);
t13 = pkin(5) * qJD(3) * t19 + t20 * t32;
t1 = [t34, 0, 0, t20 ^ 2 * t34, t20 * t24 * t23, t20 * t28, t23 * t28, qJD(2) ^ 2 / 0.2e1, t23 * t31, -t20 * t31, t29 * t24, t19 ^ 2 * t34, t19 * t33, t19 * t27, t22 * t27, qJD(3) ^ 2 / 0.2e1, t14 * t33 + (t20 * t19 + t23 * t22) * t26, -t19 * t14 * t24 + (-t23 * t19 + t22 * t20) * t26, t25 * t24, t17 ^ 2 / 0.2e1, t17 * (t18 * t13 + t21 * t30), -(-t21 * t13 + t18 * t30) * t17;];
T_reg = t1;
