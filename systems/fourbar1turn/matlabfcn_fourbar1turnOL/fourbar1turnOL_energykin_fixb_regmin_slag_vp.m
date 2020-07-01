% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1turnOL_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:28
% EndTime: 2020-06-27 16:56:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->6), mult. (81->44), div. (0->0), fcn. (44->6), ass. (0->17)
t13 = qJD(1) ^ 2;
t19 = t13 / 0.2e1;
t10 = cos(qJ(4));
t18 = t10 * t13;
t6 = qJD(2) + qJD(3);
t17 = qJD(2) * t6;
t12 = cos(qJ(2));
t16 = qJD(1) * t12;
t15 = qJD(1) * qJD(2);
t14 = qJD(1) * qJD(4);
t11 = cos(qJ(3));
t9 = sin(qJ(2));
t8 = sin(qJ(3));
t7 = sin(qJ(4));
t5 = (-t11 * t9 - t12 * t8) * qJD(1);
t4 = (-t11 * t12 + t8 * t9) * qJD(1);
t1 = [t19, 0, 0, t9 ^ 2 * t19, t9 * t13 * t12, t9 * t15, t12 * t15, qJD(2) ^ 2 / 0.2e1, 0, 0, t5 ^ 2 / 0.2e1, t5 * t4, t5 * t6, t4 * t6, t6 ^ 2 / 0.2e1, (-t11 * t17 + t4 * t16) * pkin(2), (-t5 * t16 + t8 * t17) * pkin(2), t7 ^ 2 * t19, t7 * t18, t7 * t14, t10 * t14, qJD(4) ^ 2 / 0.2e1, pkin(1) * t18, -t13 * pkin(1) * t7;];
T_reg = t1;
