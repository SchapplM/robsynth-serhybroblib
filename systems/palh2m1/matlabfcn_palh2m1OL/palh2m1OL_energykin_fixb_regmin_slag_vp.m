% Calculate minimal parameter regressor of fixed base kinetic energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh2m1OL_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:10
% EndTime: 2020-06-30 17:45:10
% DurationCPUTime: 0.12s
% Computational Cost: add. (176->33), mult. (422->80), div. (0->0), fcn. (331->8), ass. (0->33)
t121 = qJD(1) ^ 2;
t129 = t121 / 0.2e1;
t128 = pkin(2) * qJD(2);
t120 = cos(qJ(2));
t127 = t120 * t121;
t111 = qJD(2) + qJD(3);
t119 = cos(qJ(3));
t123 = t119 * t128;
t105 = t111 * pkin(3) + t123;
t114 = sin(qJ(4));
t118 = cos(qJ(4));
t115 = sin(qJ(3));
t124 = t115 * t128;
t126 = t114 * t105 + t118 * t124;
t106 = (pkin(2) * t120 + pkin(1)) * qJD(1);
t125 = qJD(1) * qJD(2);
t116 = sin(qJ(2));
t102 = (t115 * t116 - t119 * t120) * qJD(1);
t103 = (-t115 * t120 - t116 * t119) * qJD(1);
t96 = -t118 * t102 + t103 * t114;
t100 = -pkin(3) * t102 + t106;
t122 = t118 * t105 - t114 * t124;
t117 = cos(qJ(5));
t113 = sin(qJ(5));
t110 = qJD(4) + t111;
t99 = pkin(6) * t110 + t126;
t98 = -t110 * pkin(4) - t122;
t97 = t114 * t102 + t118 * t103;
t95 = qJD(5) + t96;
t94 = t113 * t110 + t117 * t97;
t93 = -t117 * t110 + t113 * t97;
t92 = pkin(4) * t96 - pkin(6) * t97 + t100;
t1 = [t129, 0, 0, t116 ^ 2 * t129, t116 * t127, -t116 * t125, -t120 * t125, qJD(2) ^ 2 / 0.2e1, pkin(1) * t127, -t121 * pkin(1) * t116, t103 ^ 2 / 0.2e1, t103 * t102, t103 * t111, t102 * t111, t111 ^ 2 / 0.2e1, -t106 * t102 + t111 * t123, t106 * t103 - t111 * t124, t97 ^ 2 / 0.2e1, -t97 * t96, t97 * t110, -t96 * t110, t110 ^ 2 / 0.2e1, t100 * t96 + t122 * t110, t100 * t97 - t126 * t110, t94 ^ 2 / 0.2e1, -t94 * t93, t94 * t95, -t93 * t95, t95 ^ 2 / 0.2e1, (-t113 * t99 + t117 * t92) * t95 + t98 * t93, -(t113 * t92 + t117 * t99) * t95 + t98 * t94;];
T_reg = t1;
