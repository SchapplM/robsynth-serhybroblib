% Calculate minimal parameter regressor of fixed base kinetic energy for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_reg [1x9]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbarprisTE_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energykin_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_energykin_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energykin_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:07:01
% EndTime: 2020-06-27 17:07:02
% DurationCPUTime: 0.04s
% Computational Cost: add. (84->14), mult. (54->24), div. (20->3), fcn. (2->2), ass. (0->13)
t13 = (qJ(1) + pkin(3));
t19 = -pkin(2) + t13;
t20 = -pkin(2) - t13;
t26 = ((pkin(1) + t20) * (pkin(1) - t19) * (pkin(1) - t20) * (pkin(1) + t19));
t25 = 1 / t26;
t16 = (t13 ^ 2);
t15 = (qJ(1) ^ 2);
t2 = pkin(1) ^ 2 - pkin(2) ^ 2 - t15 - (2 * qJ(1) + pkin(3)) * pkin(3);
t23 = t2 ^ 2 / t16;
t14 = (qJD(1) ^ 2);
t18 = (t14 * t25);
t17 = t18 * t23;
t1 = [-t17 / 0.2e1, 0, 0, t14 * t2 / t13 * (-t26) ^ (-0.1e1 / 0.2e1), -qJ(1) * t17, (-(t15 * t23 * t25) / 0.2e1 + 0.1e1 / 0.2e1) * t14, -2 * t16 * t18, 0, 0;];
T_reg = t1;
