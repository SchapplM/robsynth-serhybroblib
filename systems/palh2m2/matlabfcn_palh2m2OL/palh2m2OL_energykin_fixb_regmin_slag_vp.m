% Calculate minimal parameter regressor of fixed base kinetic energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m2OL_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:35
% EndTime: 2020-06-30 18:06:35
% DurationCPUTime: 0.14s
% Computational Cost: add. (361->42), mult. (865->99), div. (0->0), fcn. (727->10), ass. (0->41)
t146 = qJD(1) ^ 2;
t156 = t146 / 0.2e1;
t155 = cos(qJ(3));
t154 = cos(qJ(5));
t153 = pkin(4) * qJD(2);
t152 = pkin(1) * t146;
t151 = qJD(1) * qJD(2);
t137 = qJD(2) + qJD(3);
t141 = sin(qJ(3));
t150 = t141 * t153;
t136 = qJD(4) + t137;
t142 = sin(qJ(2));
t145 = cos(qJ(2));
t129 = (t141 * t142 - t145 * t155) * qJD(1);
t130 = (t141 * t145 + t155 * t142) * qJD(1);
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t123 = t144 * t129 + t140 * t130;
t124 = -t140 * t129 + t144 * t130;
t139 = sin(qJ(5));
t117 = t154 * t123 + t139 * t124;
t148 = t155 * t153;
t132 = t137 * pkin(2) + t148;
t147 = t144 * t132 - t140 * t150;
t126 = t136 * pkin(5) + t147;
t128 = t140 * t132 + t144 * t150;
t149 = t154 * t126 - t139 * t128;
t133 = (-pkin(4) * t145 - pkin(1)) * qJD(1);
t127 = t129 * pkin(2) + t133;
t120 = t123 * pkin(5) + t127;
t143 = cos(qJ(6));
t138 = sin(qJ(6));
t135 = qJD(5) + t136;
t121 = t139 * t126 + t154 * t128;
t119 = t135 * pkin(3) + t149;
t118 = -t139 * t123 + t154 * t124;
t115 = -qJD(6) + t117;
t114 = t143 * t118 - t138 * t135;
t113 = t138 * t118 + t143 * t135;
t112 = t117 * pkin(3) + t120;
t1 = [t156, 0, 0, t142 ^ 2 * t156, t142 * t146 * t145, t142 * t151, t145 * t151, qJD(2) ^ 2 / 0.2e1, t145 * t152, -t142 * t152, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t137, -t129 * t137, t137 ^ 2 / 0.2e1, t133 * t129 + t137 * t148, t133 * t130 - t137 * t150, t124 ^ 2 / 0.2e1, -t124 * t123, t124 * t136, -t123 * t136, t136 ^ 2 / 0.2e1, t127 * t123 + t136 * t147, t127 * t124 - t128 * t136, t118 ^ 2 / 0.2e1, -t118 * t117, t118 * t135, -t117 * t135, t135 ^ 2 / 0.2e1, t120 * t117 + t135 * t149, t120 * t118 - t121 * t135, t114 ^ 2 / 0.2e1, -t114 * t113, -t114 * t115, t113 * t115, t115 ^ 2 / 0.2e1, -(-t143 * t112 - t138 * t121) * t115 + t119 * t113, (-t138 * t112 + t143 * t121) * t115 + t119 * t114;];
T_reg = t1;
